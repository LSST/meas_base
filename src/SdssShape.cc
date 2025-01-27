// -*- lsst-c++ -*-
/*
 * LSST Data Management System
 * Copyright 2008-2016 AURA/LSST.
 *
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the LSST License Statement and
 * the GNU General Public License along with this program.  If not,
 * see <https://www.lsstcorp.org/LegalNotices/>.
 */

#include <cmath>
#include <tuple>

#include "boost/format.hpp"
#include "Eigen/LU"
#include "lsst/geom/Angle.h"
#include "lsst/geom/Box.h"
#include "lsst/geom/LinearTransform.h"
#include "lsst/afw/image.h"
#include "lsst/afw/detection/Psf.h"
#include "lsst/afw/geom/ellipses.h"
#include "lsst/afw/table/Source.h"
#include "lsst/meas/base/exceptions.h"
#include "lsst/meas/base/SdssShape.h"

namespace lsst {
namespace meas {
namespace base {
namespace {
FlagDefinitionList flagDefinitions;
}  // namespace

FlagDefinition const SdssShapeAlgorithm::FAILURE = flagDefinitions.addFailureFlag();
FlagDefinition const SdssShapeAlgorithm::UNWEIGHTED_BAD =
        flagDefinitions.add("flag_unweightedBad", "Both weighted and unweighted moments were invalid");
FlagDefinition const SdssShapeAlgorithm::UNWEIGHTED = flagDefinitions.add(
        "flag_unweighted", "Weighted moments converged to an invalid value; using unweighted moments");
FlagDefinition const SdssShapeAlgorithm::SHIFT =
        flagDefinitions.add("flag_shift", "centroid shifted by more than the maximum allowed amount");
FlagDefinition const SdssShapeAlgorithm::MAXITER =
        flagDefinitions.add("flag_maxIter", "Too many iterations in adaptive moments");
FlagDefinition const SdssShapeAlgorithm::PSF_SHAPE_BAD =
        flagDefinitions.add("flag_psf", "Failure in measuring PSF model shape");

FlagDefinitionList const &SdssShapeAlgorithm::getFlagDefinitions() { return flagDefinitions; }

namespace {  // anonymous

typedef Eigen::Matrix<double, 4, 4, Eigen::DontAlign> Matrix4d;

/*****************************************************************************/
/*
 * Error analysis, courtesy of David Johnston, University of Chicago
 */
/*
 * This function takes the 4 Gaussian parameters A, sigmaXXW and the
 * sky variance and fills in the Fisher matrix from the least squares fit.
 *
 * Following "Numerical Recipes in C" section 15.5, it ignores the 2nd
 * derivative parts and so the fisher matrix is just a function of these
 * best fit model parameters. The components are calculated analytically.
 */
Matrix4d calc_fisher(SdssShapeResult const &shape,  // the Shape that we want the the Fisher matrix for
                     float bkgd_var                 // background variance level for object
                     ) {
    float const A = shape.instFlux;  // amplitude; will be converted to instFlux later
    float const sigma11W = shape.xx;
    float const sigma12W = shape.xy;
    float const sigma22W = shape.yy;

    double const D = sigma11W * sigma22W - sigma12W * sigma12W;

    if (D <= std::numeric_limits<double>::epsilon()) {
        throw LSST_EXCEPT(pex::exceptions::DomainError, "Determinant is too small calculating Fisher matrix");
    }
    /*
     * a normalization factor
     */
    if (bkgd_var <= 0.0) {
        throw LSST_EXCEPT(pex::exceptions::DomainError,
                          (boost::format("Background variance must be positive (saw %g)") % bkgd_var).str());
    }
    double const F = geom::PI * sqrt(D) / bkgd_var;
    /*
     * Calculate the 10 independent elements of the 4x4 Fisher matrix
     */
    Matrix4d fisher;

    double fac = F * A / (4.0 * D);
    fisher(0, 0) = F;
    fisher(0, 1) = fac * sigma22W;
    fisher(1, 0) = fisher(0, 1);
    fisher(0, 2) = fac * sigma11W;
    fisher(2, 0) = fisher(0, 2);
    fisher(0, 3) = -fac * 2 * sigma12W;
    fisher(3, 0) = fisher(0, 3);

    fac = 3.0 * F * A * A / (16.0 * D * D);
    fisher(1, 1) = fac * sigma22W * sigma22W;
    fisher(2, 2) = fac * sigma11W * sigma11W;
    fisher(3, 3) = fac * 4.0 * (sigma12W * sigma12W + D / 3.0);

    fisher(1, 2) = fisher(3, 3) / 4.0;
    fisher(2, 1) = fisher(1, 2);
    fisher(1, 3) = fac * (-2 * sigma22W * sigma12W);
    fisher(3, 1) = fisher(1, 3);
    fisher(2, 3) = fac * (-2 * sigma11W * sigma12W);
    fisher(3, 2) = fisher(2, 3);

    return fisher;
}
//
// Here's a class to allow us to get the Image and variance from an Image or MaskedImage
//
template <typename ImageT>  // general case
struct ImageAdaptor {
    typedef ImageT Image;

    static bool const hasVariance = false;

    Image const &getImage(ImageT const &image) const { return image; }

    double getVariance(ImageT const &, int, int) { return std::numeric_limits<double>::quiet_NaN(); }
};

template <typename T>  // specialise to a MaskedImage
struct ImageAdaptor<afw::image::MaskedImage<T> > {
    typedef typename afw::image::MaskedImage<T>::Image Image;

    static bool const hasVariance = true;

    Image const &getImage(afw::image::MaskedImage<T> const &mimage) const { return *mimage.getImage(); }

    double getVariance(afw::image::MaskedImage<T> const &mimage, int ix, int iy) {
        return mimage.at(ix, iy).variance();
    }
};

/// Calculate weights from moments
std::tuple<std::pair<bool, double>, double, double, double> getWeights(double sigma11, double sigma12,
                                                                       double sigma22) {
    double const NaN = std::numeric_limits<double>::quiet_NaN();
    if (std::isnan(sigma11) || std::isnan(sigma12) || std::isnan(sigma22)) {
        return std::make_tuple(std::make_pair(false, NaN), NaN, NaN, NaN);
    }
    double const det = sigma11 * sigma22 - sigma12 * sigma12;              // determinant of sigmaXX matrix
    if (std::isnan(det) || det < std::numeric_limits<float>::epsilon()) {  // a suitably small number
        return std::make_tuple(std::make_pair(false, det), NaN, NaN, NaN);
    }
    return std::make_tuple(std::make_pair(true, det), sigma22 / det, -sigma12 / det, sigma11 / det);
}

/// Should we be interpolating?
bool shouldInterp(double sigma11, double sigma22, double det) {
    float const xinterp = 0.25;  // I.e. 0.5*0.5
    return (sigma11 < xinterp || sigma22 < xinterp || det < xinterp * xinterp);
}

// Decide on the bounding box for the region to examine while calculating the adaptive moments
// This routine will work in either LOCAL or PARENT coordinates (but of course which you pass
// determine which you will get back).
geom::Box2I computeAdaptiveMomentsBBox(geom::Box2I const &bbox,      // full image bbox
                                       geom::Point2D const &center,  // centre of object
                                       double sigma11_w,             // quadratic moments of the
                                       double,                       //         weighting function
                                       double sigma22_w,             //                    xx, xy, and yy
                                       double maxRadius = 1000       // Maximum radius of area to use
                                       ) {
    double radius = std::min(4 * std::sqrt(std::max(sigma11_w, sigma22_w)), maxRadius);
    geom::Extent2D offset(radius, radius);
    geom::Box2I result(geom::Box2D(center - offset, center + offset));
    result.clip(bbox);
    return result;
}

/*****************************************************************************/
/*
 * Calculate weighted moments of an object up to 2nd order
 */
template <typename ImageT>
static void calcmom(ImageT const &image,                             // the image data
                   float xcen, float ycen,                          // centre of object
                   const geom::BoxI &bbox,                          // bounding box to consider
                   float bkgd,                                      // data's background level
                   bool interpflag,                                 // interpolate within pixels?
                   double w11, double w12, double w22,              // weights
                   double &psum) {                                    // sum w*I

    float tmod, ymod;
    float X, Y;  // sub-pixel interpolated [xy]
    float weight;
    float tmp;
    double sum, sumx, sumy, sumxx, sumyy, sumxy, sums4;

    if (w11 < 0 ||  w11 > 1e6 || fabs(w12) > 1E6 || w22 < 0 || w22 > 1e6) {
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterError, "Invalid weight parameter(s)");
    }

    sum = sumx = sumy = sumxx = sumxy = sumyy = sums4 = 0;

    int const ix0 = bbox.getMinX();  // corners of the box being analyzed
    int const ix1 = bbox.getMaxX();
    int const iy0 = bbox.getMinY();  // corners of the box being analyzed
    int const iy1 = bbox.getMaxY();

    if (ix0 < 0 || ix1 >= image.getWidth() || iy0 < 0 || iy1 >= image.getHeight()) {
        throw LSST_EXCEPT(pex::exceptions::LengthError, "Invalid image dimensions");
    }

    for (int i = iy0; i <= iy1; ++i) {
        typename ImageT::x_iterator ptr = image.x_at(ix0, i);
        float const y = i - ycen;
        float const y2 = y * y;
        float const yl = y - 0.375;
        float const yh = y + 0.375;
        for (int j = ix0; j <= ix1; ++j, ++ptr) {
            float x = j - xcen;
            if (interpflag) {
                float const xl = x - 0.375;
                float const xh = x + 0.375;

                float expon = xl * xl * w11 + yl * yl * w22 + 2.0 * xl * yl * w12;
                tmp = xh * xh * w11 + yh * yh * w22 + 2.0 * xh * yh * w12;
                expon = (expon > tmp) ? expon : tmp;
                tmp = xl * xl * w11 + yh * yh * w22 + 2.0 * xl * yh * w12;
                expon = (expon > tmp) ? expon : tmp;
                tmp = xh * xh * w11 + yl * yl * w22 + 2.0 * xh * yl * w12;
                expon = (expon > tmp) ? expon : tmp;

                if (expon <= 9.0) {
                    tmod = *ptr - bkgd;
                    for (Y = yl; Y <= yh; Y += 0.25) {
                        double const interpY2 = Y * Y;
                        for (X = xl; X <= xh; X += 0.25) {
                            double const interpX2 = X * X;
                            double const interpXy = X * Y;
                            expon = interpX2 * w11 + 2 * interpXy * w12 + interpY2 * w22;
                            weight = std::exp(-0.5 * expon);

                            ymod = tmod * weight;
                            sum += ymod;
                        }
                    }
                }
            } else {
                float x2 = x * x;
                float xy = x * y;
                float expon = x2 * w11 + 2 * xy * w12 + y2 * w22;

                if (expon <= 14.0) {
                    weight = std::exp(-0.5 * expon);
                    tmod = *ptr - bkgd;
                    ymod = tmod * weight;
                    sum += ymod;
                }
            }
        }
    }

    psum = sum;

}

template <typename ImageT>
static int calcmom(ImageT const &image,                             // the image data
                   float xcen, float ycen,                          // centre of object
                   const geom::BoxI& bbox,                                 // bounding box to consider
                   float bkgd,                                      // data's background level
                   bool interpflag,                                 // interpolate within pixels?
                   double w11, double w12, double w22,              // weights
                   double &pI0,                                     // amplitude of fit (if !NULL)
                   double &psum,                                    // sum w*I (if !NULL)
                   double &psumx, double &psumy,                    // sum [xy]*w*I (if !instFluxOnly)
                   double &psumxx, double &psumxy, double &psumyy,  // sum [xy]^2*w*I (if !instFluxOnly)
                   double &psums4,  // sum w*I*weight^2 (if !instFluxOnly && !NULL)
                   bool negative = false) {
    float tmod, ymod;
    float X, Y;  // sub-pixel interpolated [xy]
    float weight;
    float tmp;
    double sum, sumx, sumy, sumxx, sumyy, sumxy, sums4;

    if (w11 < 0 ||  w11 > 1e6 || fabs(w12) > 1E6 || w22 < 0 || w22 > 1e6) {
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterError, "Invalid weight parameter(s)");
    }

    sum = sumx = sumy = sumxx = sumxy = sumyy = sums4 = 0;

    int const ix0 = bbox.getMinX();  // corners of the box being analyzed
    int const ix1 = bbox.getMaxX();
    int const iy0 = bbox.getMinY();  // corners of the box being analyzed
    int const iy1 = bbox.getMaxY();

    if (ix0 < 0 || ix1 >= image.getWidth() || iy0 < 0 || iy1 >= image.getHeight()) {
        throw LSST_EXCEPT(pex::exceptions::LengthError, "Invalid image dimensions");
    }

    for (int i = iy0; i <= iy1; ++i) {
        typename ImageT::x_iterator ptr = image.x_at(ix0, i);
        float const y = i - ycen;
        float const y2 = y * y;
        float const yl = y - 0.375;
        float const yh = y + 0.375;
        for (int j = ix0; j <= ix1; ++j, ++ptr) {
            float x = j - xcen;
            if (interpflag) {
                float const xl = x - 0.375;
                float const xh = x + 0.375;

                float expon = xl * xl * w11 + yl * yl * w22 + 2.0 * xl * yl * w12;
                tmp = xh * xh * w11 + yh * yh * w22 + 2.0 * xh * yh * w12;
                expon = (expon > tmp) ? expon : tmp;
                tmp = xl * xl * w11 + yh * yh * w22 + 2.0 * xl * yh * w12;
                expon = (expon > tmp) ? expon : tmp;
                tmp = xh * xh * w11 + yl * yl * w22 + 2.0 * xh * yl * w12;
                expon = (expon > tmp) ? expon : tmp;

                if (expon <= 9.0) {
                    tmod = *ptr - bkgd;
                    for (Y = yl; Y <= yh; Y += 0.25) {
                        double const interpY2 = Y * Y;
                        for (X = xl; X <= xh; X += 0.25) {
                            double const interpX2 = X * X;
                            double const interpXy = X * Y;
                            expon = interpX2 * w11 + 2 * interpXy * w12 + interpY2 * w22;
                            weight = std::exp(-0.5 * expon);

                            ymod = tmod * weight;
                            sum += ymod;
                            sumx += ymod * (X + xcen);
                            sumy += ymod * (Y + ycen);
                            sumxx += interpX2 * ymod;
                            sumxy += interpXy * ymod;
                            sumyy += interpY2 * ymod;
                            sums4 += expon * expon * ymod;
                        }
                    }
                }
            } else {
                float x2 = x * x;
                float xy = x * y;
                float expon = x2 * w11 + 2 * xy * w12 + y2 * w22;

                if (expon <= 14.0) {
                    weight = std::exp(-0.5 * expon);
                    tmod = *ptr - bkgd;
                    ymod = tmod * weight;
                    sum += ymod;
                    sumx += ymod * j;
                    sumy += ymod * i;
                    sumxx += x2 * ymod;
                    sumxy += xy * ymod;
                    sumyy += y2 * ymod;
                    sums4 += expon * expon * ymod;
                }
            }
        }
    }


    std::tuple<std::pair<bool, double>, double, double, double> const weights = getWeights(w11, w12, w22);
    double const detW = std::get<1>(weights) * std::get<3>(weights) - std::pow(std::get<2>(weights), 2);
    pI0 = sum / (geom::PI * sqrt(detW));


    psum = sum;

    psumx = sumx;
    psumy = sumy;
    psumxx = sumxx;
    psumxy = sumxy;
    psumyy = sumyy;

    psums4 = sums4;

    if (negative) {
        return (sum < 0 && sumxx < 0 && sumyy < 0) ? 0 : -1;
    } else {
        return (sum > 0 && sumxx > 0 && sumyy > 0) ? 0 : -1;
    }
}

/*
 * Workhorse for adaptive moments
 *
 * All inputs are expected to be in LOCAL image coordinates
 */
template <typename ImageT>
bool getAdaptiveMoments(ImageT const &mimage, double bkgd, double xcen, double ycen, double shiftmax,
                        SdssShapeResult *shape, int maxIter, float tol1, float tol2, bool negative) {
    double I0 = 0;               // amplitude of best-fit Gaussian
    double sum;                  // sum of intensity*weight
    double sumx, sumy;           // sum ((int)[xy])*intensity*weight
    double sumxx, sumxy, sumyy;  // sum {x^2,xy,y^2}*intensity*weight
    double sums4;                // sum intensity*weight*exponent^2
    float const xcen0 = xcen;    // initial centre
    float const ycen0 = ycen;    //                of object

    double sigma11W = 1.5;  // quadratic moments of the
    double sigma12W = 0.0;  //     weighting fcn;
    double sigma22W = 1.5;  //               xx, xy, and yy

    double w11 = -1, w12 = -1, w22 = -1;  // current weights for moments; always set when iter == 0
    float e1_old = 1e6, e2_old = 1e6;     // old values of shape parameters e1 and e2
    float sigma11_ow_old = 1e6;           // previous version of sigma11_ow

    typename ImageAdaptor<ImageT>::Image const &image = ImageAdaptor<ImageT>().getImage(mimage);

    if (std::isnan(xcen) || std::isnan(ycen)) {
        // Can't do anything
        shape->flags[SdssShapeAlgorithm::UNWEIGHTED_BAD.number] = true;
        return false;
    }

    bool interpflag = false;  // interpolate finer than a pixel?
    geom::BoxI bbox;
    int iter = 0;  // iteration number
    for (; iter < maxIter; iter++) {
        bbox = computeAdaptiveMomentsBBox(image.getBBox(afw::image::LOCAL), geom::Point2D(xcen, ycen),
                                          sigma11W, sigma12W, sigma22W);
        std::tuple<std::pair<bool, double>, double, double, double> weights =
                getWeights(sigma11W, sigma12W, sigma22W);
        if (!std::get<0>(weights).first) {
            shape->flags[SdssShapeAlgorithm::UNWEIGHTED.number] = true;
            break;
        }

        double const detW = std::get<0>(weights).second;

        if( sigma11W * sigma22W < sigma12W * sigma12W - std::numeric_limits<float>::epsilon()) return false;

        {
            const double ow11 = w11;  // old
            const double ow12 = w12;  //     values
            const double ow22 = w22;  //            of w11, w12, w22

            w11 = std::get<1>(weights);
            w12 = std::get<2>(weights);
            w22 = std::get<3>(weights);

            if (shouldInterp(sigma11W, sigma22W, detW)) {
                if (!interpflag) {
                    interpflag = true;  // N.b.: stays set for this object
                    if (iter > 0) {
                        sigma11_ow_old = 1.e6;  // force at least one more iteration
                        w11 = ow11;
                        w12 = ow12;
                        w22 = ow22;
                        iter--;  // we didn't update wXX
                    }
                }
            }
        }
         if (calcmom(image, xcen, ycen, bbox, bkgd, interpflag, w11, w12, w22, I0, sum, sumx, sumy,
                           sumxx, sumxy, sumyy, sums4, negative) < 0) {
            shape->flags[SdssShapeAlgorithm::UNWEIGHTED.number] = true;
            break;
        }

        shape->x = sumx / sum;  // update centroid.  N.b. we're not setting errors here
        shape->y = sumy / sum;

        if (fabs(shape->x - xcen0) > shiftmax || fabs(shape->y - ycen0) > shiftmax) {
            shape->flags[SdssShapeAlgorithm::SHIFT.number] = true;
        }
        /*
         * OK, we have the centre. Proceed to find the second moments.
         */
        float const sigma11_ow = sumxx / sum;  // quadratic moments of
        float const sigma22_ow = sumyy / sum;  //          weight*object
        float const sigma12_ow = sumxy / sum;  //                 xx, xy, and yy

        if (sigma11_ow <= 0 || sigma22_ow <= 0) {
            shape->flags[SdssShapeAlgorithm::UNWEIGHTED.number] = true;
            break;
        }

        float const d = sigma11_ow + sigma22_ow;  // current values of shape parameters
        float const e1 = (sigma11_ow - sigma22_ow) / d;
        float const e2 = 2.0 * sigma12_ow / d;
        /*
         * Did we converge?
         */
        if (iter > 0 && fabs(e1 - e1_old) < tol1 && fabs(e2 - e2_old) < tol1 &&
            fabs(sigma11_ow / sigma11_ow_old - 1.0) < tol2) {
            break;  // yes; we converged
        }

        e1_old = e1;
        e2_old = e2;
        sigma11_ow_old = sigma11_ow;
        /*
         * Didn't converge, calculate new values for weighting function
         *
         * The product of two Gaussians is a Gaussian:
         * <x^2 exp(-a x^2 - 2bxy - cy^2) exp(-Ax^2 - 2Bxy - Cy^2)> =
         *                            <x^2 exp(-(a + A) x^2 - 2(b + B)xy - (c + C)y^2)>
         * i.e. the inverses of the covariances matrices add.
         *
         * We know sigmaXX_ow and sigmaXXW, the covariances of the weighted object
         * and of the weights themselves.  We can estimate the object's covariance as
         *   sigmaXX_ow^-1 - sigmaXXW^-1
         * and, as we want to find a set of weights with the _same_ covariance as the
         * object we take this to be the an estimate of our correct weights.
         *
         * N.b. This assumes that the object is roughly Gaussian.
         * Consider the object:
         *   O == delta(x + p) + delta(x - p)
         * the covariance of the weighted object is equal to that of the unweighted
         * object, and this prescription fails badly.  If we detect this, we set
         * the UNWEIGHTED flag, and calculate the UNweighted moments
         * instead.
         */
        {
            float n11, n12, n22;     // elements of inverse of next guess at weighting function
            float ow11, ow12, ow22;  // elements of inverse of sigmaXX_ow

            std::tuple<std::pair<bool, double>, double, double, double> weights =
                    getWeights(sigma11_ow, sigma12_ow, sigma22_ow);
            if (!std::get<0>(weights).first) {
                shape->flags[SdssShapeAlgorithm::UNWEIGHTED.number] = true;
                break;
            }

            ow11 = std::get<1>(weights);
            ow12 = std::get<2>(weights);
            ow22 = std::get<3>(weights);

            n11 = ow11 - w11;
            n12 = ow12 - w12;
            n22 = ow22 - w22;

            weights = getWeights(n11, n12, n22);
            if (!std::get<0>(weights).first) {
                // product-of-Gaussians assumption failed
                shape->flags[SdssShapeAlgorithm::UNWEIGHTED.number] = true;
                break;
            }

            sigma11W = std::get<1>(weights);
            sigma12W = std::get<2>(weights);
            sigma22W = std::get<3>(weights);
        }

        if (sigma11W <= 0 || sigma22W <= 0) {
            shape->flags[SdssShapeAlgorithm::UNWEIGHTED.number] = true;
            break;
        }
    }

    if (iter == maxIter) {
        shape->flags[SdssShapeAlgorithm::UNWEIGHTED.number] = true;
        shape->flags[SdssShapeAlgorithm::MAXITER.number] = true;
    }

    if (sumxx + sumyy == 0.0) {
        shape->flags[SdssShapeAlgorithm::UNWEIGHTED.number] = true;
    }
    /*
     * Problems; try calculating the un-weighted moments
     */
    if (shape->flags[SdssShapeAlgorithm::UNWEIGHTED.number]) {
        w11 = w22 = w12 = 0;
        double ignored;
        if (calcmom(image, xcen, ycen, bbox, bkgd, interpflag, w11, w12, w22, I0, sum, sumx, sumy,
                           sumxx, sumxy, sumyy, ignored, negative) < 0 ||
            (!negative && sum <= 0) || (negative && sum >= 0)) {
            shape->flags[SdssShapeAlgorithm::UNWEIGHTED.number] = false;
            shape->flags[SdssShapeAlgorithm::UNWEIGHTED_BAD.number] = true;

            if (sum > 0) {
                shape->xx = 1 / 12.0;  // a single pixel
                shape->xy = 0.0;
                shape->yy = 1 / 12.0;
            }

            return false;
        }

        sigma11W = sumxx / sum;  // estimate of object moments
        sigma12W = sumxy / sum;  //   usually, object == weight
        sigma22W = sumyy / sum;  //      at this point
    }

    shape->instFlux = I0;
    shape->xx = sigma11W;
    shape->xy = sigma12W;
    shape->yy = sigma22W;

    if (shape->xx + shape->yy != 0.0) {
        int const ix = afw::image::positionToIndex(xcen);
        int const iy = afw::image::positionToIndex(ycen);

        if (ix >= 0 && ix < mimage.getWidth() && iy >= 0 && iy < mimage.getHeight()) {
            float const bkgd_var = ImageAdaptor<ImageT>().getVariance(
                    mimage, ix, iy);  // XXX Overestimate as it includes object

            if (bkgd_var > 0.0) {  // NaN is not > 0.0
                if (!(shape->flags[SdssShapeAlgorithm::UNWEIGHTED.number])) {
                    Matrix4d fisher = calc_fisher(*shape, bkgd_var);  // Fisher matrix
                    Matrix4d cov = fisher.inverse();
                    // convention is to order moments (xx, yy, xy)
                    shape->instFluxErr = std::sqrt(cov(0, 0));
                    shape->xxErr = std::sqrt(cov(1, 1));
                    shape->yyErr = std::sqrt(cov(2, 2));
                    shape->xyErr = std::sqrt(cov(3, 3));
                    shape->instFlux_xx_Cov = cov(0, 1);
                    shape->instFlux_yy_Cov = cov(0, 2);
                    shape->instFlux_xy_Cov = cov(0, 3);
                    shape->xx_yy_Cov = cov(1, 2);
                    shape->xx_xy_Cov = cov(1, 3);
                    shape->yy_xy_Cov = cov(2, 3);
                }
            }
        }
    }

    return true;
}

}  // namespace

SdssShapeResult::SdssShapeResult()
        : instFlux_xx_Cov(std::numeric_limits<ErrElement>::quiet_NaN()),
          instFlux_yy_Cov(std::numeric_limits<ErrElement>::quiet_NaN()),
          instFlux_xy_Cov(std::numeric_limits<ErrElement>::quiet_NaN()) {}

SdssShapeResultKey SdssShapeResultKey::addFields(afw::table::Schema &schema, std::string const &name,
                                                 bool doMeasurePsf) {
    SdssShapeResultKey r;
    r._shapeResult =
            ShapeResultKey::addFields(schema, name, "elliptical Gaussian adaptive moments", SIGMA_ONLY);
    r._centroidResult = CentroidResultKey::addFields(schema, name, "elliptical Gaussian adaptive moments",
                                                     NO_UNCERTAINTY);
    r._instFluxResult = FluxResultKey::addFields(schema, name, "elliptical Gaussian adaptive moments");

    // Only include PSF shape fields if doMeasurePsf = True
    if (doMeasurePsf) {
        r._includePsf = true;
        r._psfShapeResult = afw::table::QuadrupoleKey::addFields(
                schema, schema.join(name, "psf"), "adaptive moments of the PSF model at the object position");
    } else {
        r._includePsf = false;
    }

    r._instFlux_xx_Cov =
            schema.addField<ErrElement>(schema.join(name, "instFlux", "xx", "Cov"),
                                        (boost::format("uncertainty covariance between %s and %s") %
                                         schema.join(name, "instFlux") % schema.join(name, "xx"))
                                                .str(),
                                        "count*pixel^2");
    r._instFlux_yy_Cov =
            schema.addField<ErrElement>(schema.join(name, "instFlux", "yy", "Cov"),
                                        (boost::format("uncertainty covariance between %s and %s") %
                                         schema.join(name, "instFlux") % schema.join(name, "yy"))
                                                .str(),
                                        "count*pixel^2");
    r._instFlux_xy_Cov =
            schema.addField<ErrElement>(schema.join(name, "instFlux", "xy", "Cov"),
                                        (boost::format("uncertainty covariance between %s and %s") %
                                         schema.join(name, "instFlux") % schema.join(name, "xy"))
                                                .str(),
                                        "count*pixel^2");

    // Skip the psf flag if not recording the PSF shape.
    if (r._includePsf) {
        r._flagHandler = FlagHandler::addFields(schema, name, SdssShapeAlgorithm::getFlagDefinitions());
    } else {
        r._flagHandler = FlagHandler::addFields(schema, name, SdssShapeAlgorithm::getFlagDefinitions(),
                                                {SdssShapeAlgorithm::PSF_SHAPE_BAD});
    }
    return r;
}

SdssShapeResultKey::SdssShapeResultKey(afw::table::SubSchema const &s)
        : _shapeResult(s),
          _centroidResult(s),
          _instFluxResult(s),
          _instFlux_xx_Cov(s["instFlux"]["xx"]["Cov"]),
          _instFlux_yy_Cov(s["instFlux"]["yy"]["Cov"]),
          _instFlux_xy_Cov(s["instFlux"]["xy"]["Cov"]) {
    // The input SubSchema may optionally provide for a PSF.
    try {
        _psfShapeResult = afw::table::QuadrupoleKey(s["psf"]);
        _flagHandler = FlagHandler(s, SdssShapeAlgorithm::getFlagDefinitions());
        _includePsf = true;
    } catch (pex::exceptions::NotFoundError &e) {
        _flagHandler =
                FlagHandler(s, SdssShapeAlgorithm::getFlagDefinitions(), {SdssShapeAlgorithm::PSF_SHAPE_BAD});
        _includePsf = false;
    }
}

SdssShapeResult SdssShapeResultKey::get(afw::table::BaseRecord const &record) const {
    SdssShapeResult result;
    static_cast<ShapeResult &>(result) = record.get(_shapeResult);
    static_cast<CentroidResult &>(result) = record.get(_centroidResult);
    static_cast<FluxResult &>(result) = record.get(_instFluxResult);
    result.instFlux_xx_Cov = record.get(_instFlux_xx_Cov);
    result.instFlux_yy_Cov = record.get(_instFlux_yy_Cov);
    result.instFlux_xy_Cov = record.get(_instFlux_xy_Cov);
    for (size_t n = 0; n < SdssShapeAlgorithm::N_FLAGS; ++n) {
        if (n == SdssShapeAlgorithm::PSF_SHAPE_BAD.number && !_includePsf) continue;
        result.flags[n] = _flagHandler.getValue(record, n);
    }
    return result;
}

afw::geom::ellipses::Quadrupole SdssShapeResultKey::getPsfShape(afw::table::BaseRecord const &record) const {
    return record.get(_psfShapeResult);
}

void SdssShapeResultKey::set(afw::table::BaseRecord &record, SdssShapeResult const &value) const {
    record.set(_shapeResult, value);
    record.set(_centroidResult, value);
    record.set(_instFluxResult, value);
    record.set(_instFlux_xx_Cov, value.instFlux_xx_Cov);
    record.set(_instFlux_yy_Cov, value.instFlux_yy_Cov);
    record.set(_instFlux_xy_Cov, value.instFlux_xy_Cov);
    for (size_t n = 0; n < SdssShapeAlgorithm::N_FLAGS; ++n) {
        if (n == SdssShapeAlgorithm::PSF_SHAPE_BAD.number && !_includePsf) continue;
        _flagHandler.setValue(record, n, value.flags[n]);
    }
}

void SdssShapeResultKey::setPsfShape(afw::table::BaseRecord &record,
                                     afw::geom::ellipses::Quadrupole const &value) const {
    record.set(_psfShapeResult, value);
}

bool SdssShapeResultKey::operator==(SdssShapeResultKey const &other) const {
    return _shapeResult == other._shapeResult && _centroidResult == other._centroidResult &&
           _instFluxResult == other._instFluxResult && _psfShapeResult == other._psfShapeResult &&
           _instFlux_xx_Cov == other._instFlux_xx_Cov && _instFlux_yy_Cov == other._instFlux_yy_Cov &&
           _instFlux_xy_Cov == other._instFlux_xy_Cov;
    // don't bother with flags - if we've gotten this far, it's basically impossible the flags don't match
}

bool SdssShapeResultKey::isValid() const {
    return _shapeResult.isValid() && _centroidResult.isValid() && _instFluxResult.isValid() &&
           _psfShapeResult.isValid() && _instFlux_xx_Cov.isValid() && _instFlux_yy_Cov.isValid() &&
           _instFlux_xy_Cov.isValid();
    // don't bother with flags - if we've gotten this far, it's basically impossible the flags are invalid
}

SdssShapeAlgorithm::SdssShapeAlgorithm(Control const &ctrl, std::string const &name,
                                       afw::table::Schema &schema)
        : _ctrl(ctrl),
          _resultKey(ResultKey::addFields(schema, name, ctrl.doMeasurePsf)),
          _centroidExtractor(schema, name) {}

template <typename ImageT>
SdssShapeResult SdssShapeAlgorithm::computeAdaptiveMoments(ImageT const &image, geom::Point2D const &center,
                                                           bool negative, Control const &control) {
    double xcen = center.getX();  // object's column position
    double ycen = center.getY();  // object's row position

    xcen -= image.getX0();  // work in image Pixel coordinates
    ycen -= image.getY0();

    float shiftmax = control.maxShift;  // Max allowed centroid shift
    if (shiftmax < 2) {
        shiftmax = 2;
    } else if (shiftmax > 10) {
        shiftmax = 10;
    }

    SdssShapeResult result;
    try {
        result.flags[FAILURE.number] =
                !getAdaptiveMoments(image, control.background, xcen, ycen, shiftmax, &result, control.maxIter,
                                    control.tol1, control.tol2, negative);
    } catch (pex::exceptions::Exception &err) {
        result.flags[FAILURE.number] = true;
    }
    if (result.flags[UNWEIGHTED.number] || result.flags[SHIFT.number]) {
        // These are also considered fatal errors in terms of the quality of the results,
        // even though they do produce some results.
        result.flags[FAILURE.number] = true;
    }
    double IxxIyy = result.getQuadrupole().getIxx() * result.getQuadrupole().getIyy();
    double Ixy_sq = result.getQuadrupole().getIxy();
    Ixy_sq *= Ixy_sq;
    double epsilon = 1.0e-6;
    if (IxxIyy < (1.0 + epsilon) * Ixy_sq)
    // We are checking that Ixx*Iyy > (1 + epsilon)*Ixy*Ixy where epsilon is suitably small. The
    // value of epsilon used here is a magic number subject to future review (DM-5801 was marked won't fix).
    {
        if (!result.flags[FAILURE.number]) {
            throw LSST_EXCEPT(pex::exceptions::LogicError,
                              (boost::format("computeAdaptiveMoments IxxIxy %d < (1 + eps=%d)*(Ixy^2=%d);"
                                             " implying singular moments without any flag set")
                              % IxxIyy % epsilon % Ixy_sq).str());
        }
    }

    // getAdaptiveMoments() just computes the zeroth moment in result.instFlux (and its error in
    // result.instFluxErr, result.instFlux_xx_Cov, etc.)  That's related to the instFlux by some geometric
    // factors, which we apply here.
    // This scale factor is just the inverse of the normalization constant in a bivariate normal distribution:
    // 2*pi*sigma_x*sigma_y*sqrt(1-rho^2), where rho is the correlation.
    // This happens to be twice the ellipse area pi*sigma_maj*sigma_min, which is identically equal to:
    // pi*sqrt(I_xx*I_yy - I_xy^2); i.e. pi*sqrt(determinant(I))
    double instFluxScale = geom::TWOPI * std::sqrt(IxxIyy - Ixy_sq);

    result.instFlux *= instFluxScale;
    result.instFluxErr *= instFluxScale;
    result.x += image.getX0();
    result.y += image.getY0();

    if (ImageAdaptor<ImageT>::hasVariance) {
        result.instFlux_xx_Cov *= instFluxScale;
        result.instFlux_yy_Cov *= instFluxScale;
        result.instFlux_xy_Cov *= instFluxScale;
    }

    return result;
}

template <typename ImageT>
FluxResult SdssShapeAlgorithm::computeFixedMomentsFlux(ImageT const &image,
                                                       afw::geom::ellipses::Quadrupole const &shape,
                                                       geom::Point2D const &center) {
    // while arguments to computeFixedMomentsFlux are in PARENT coordinates, the implementation is LOCAL.
    geom::Point2D localCenter = center - geom::Extent2D(image.getXY0());

    geom::BoxI const bbox = computeAdaptiveMomentsBBox(image.getBBox(afw::image::LOCAL), localCenter,
                                                       shape.getIxx(), shape.getIxy(), shape.getIyy());

    std::tuple<std::pair<bool, double>, double, double, double> weights =
            getWeights(shape.getIxx(), shape.getIxy(), shape.getIyy());

    FluxResult result;

    if (!std::get<0>(weights).first) {
        throw pex::exceptions::InvalidParameterError("Input shape is singular");
    }

    double const w11 = std::get<1>(weights);
    double const w12 = std::get<2>(weights);
    double const w22 = std::get<3>(weights);
    bool const interp = shouldInterp(shape.getIxx(), shape.getIyy(), std::get<0>(weights).second);

    double sum0 = 0;  //  sum of pixel values weighted by a Gaussian
    calcmom(ImageAdaptor<ImageT>().getImage(image), localCenter.getX(), localCenter.getY(), bbox,
                      0.0, interp, w11, w12, w22, sum0);

    result.instFlux = sum0 * 2.0;

    if (ImageAdaptor<ImageT>::hasVariance) {
        int ix = static_cast<int>(center.getX() - image.getX0());
        int iy = static_cast<int>(center.getY() - image.getY0());
        if (!image.getBBox(afw::image::LOCAL).contains(geom::Point2I(ix, iy))) {
            throw LSST_EXCEPT(pex::exceptions::RuntimeError,
                              (boost::format("Center (%d,%d) not in image (%dx%d)") % ix % iy %
                               image.getWidth() % image.getHeight())
                                      .str());
        }
        double var = ImageAdaptor<ImageT>().getVariance(image, ix, iy);
        // 0th moment (i0) error = sqrt(var / wArea); instFlux (error) = 2 * wArea * i0 (error)
        double const wArea = geom::PI * std::sqrt(shape.getDeterminant());
        result.instFluxErr = 2 * std::sqrt(var * wArea);
    }

    return result;
}

void SdssShapeAlgorithm::measure(afw::table::SourceRecord &measRecord,
                                 afw::image::Exposure<float> const &exposure) const {
    bool negative = false;

    try {
        negative = measRecord.get(measRecord.getSchema().find<afw::table::Flag>("is_negative").key);
    } catch (pexExcept::Exception &e) {
    }
    SdssShapeResult result = computeAdaptiveMoments(
            exposure.getMaskedImage(), _centroidExtractor(measRecord, _resultKey.getFlagHandler()), negative,
            _ctrl);

    if (_ctrl.doMeasurePsf) {
        // Compute moments of Psf model.  In the interest of implementing this quickly, we're just
        // calling Psf::computeShape(), which delegates to SdssShapeResult::computeAdaptiveMoments
        // for all nontrivial Psf classes.  But this could in theory save the results of a shape
        // computed some other way as part of base_SdssShape, which might be confusing.  We should
        // fix this eventually either by making Psf shape measurement not part of base_SdssShape, or
        // by making the measurements stored with shape.sdss always computed via the
        // SdssShapeAlgorithm instead of delegating to the Psf class.
        try {
            std::shared_ptr<afw::detection::Psf const> psf = exposure.getPsf();
            if (!psf) {
                result.flags[PSF_SHAPE_BAD.number] = true;
            } else {
                _resultKey.setPsfShape(measRecord, psf->computeShape(geom::Point2D(result.x, result.y)));
            }
        } catch (pex::exceptions::Exception &err) {
            result.flags[PSF_SHAPE_BAD.number] = true;
        }
    }

    measRecord.set(_resultKey, result);
}

void SdssShapeAlgorithm::fail(afw::table::SourceRecord &measRecord, MeasurementError *error) const {
    _resultKey.getFlagHandler().handleFailure(measRecord, error);
}

#define INSTANTIATE_IMAGE(IMAGE)                                          \
    template SdssShapeResult SdssShapeAlgorithm::computeAdaptiveMoments(  \
            IMAGE const &, geom::Point2D const &, bool, Control const &); \
    template FluxResult SdssShapeAlgorithm::computeFixedMomentsFlux(      \
            IMAGE const &, afw::geom::ellipses::Quadrupole const &, geom::Point2D const &)

#define INSTANTIATE_PIXEL(PIXEL)                 \
    INSTANTIATE_IMAGE(afw::image::Image<PIXEL>); \
    INSTANTIATE_IMAGE(afw::image::MaskedImage<PIXEL>);

INSTANTIATE_PIXEL(int);
INSTANTIATE_PIXEL(float);
INSTANTIATE_PIXEL(double);

SdssShapeTransform::SdssShapeTransform(Control const &ctrl, std::string const &name,
                                       afw::table::SchemaMapper &mapper)
        : BaseTransform{name}, _instFluxTransform{name, mapper}, _centroidTransform{name, mapper} {
    // If the input schema has a PSF flag -- it's optional --  assume we are also transforming the PSF.
    _transformPsf = mapper.getInputSchema().getNames().count("sdssShape_flag_psf") ? true : false;

    // Skip the last flag if not transforming the PSF shape.
    for (std::size_t i = 0; i < SdssShapeAlgorithm::getFlagDefinitions().size(); i++) {
        FlagDefinition const &flag = SdssShapeAlgorithm::getFlagDefinitions()[i];
        if (flag == SdssShapeAlgorithm::FAILURE) continue;
        if (mapper.getInputSchema().getNames().count(mapper.getInputSchema().join(name, flag.name)) == 0)
            continue;
        afw::table::Key<afw::table::Flag> key =
                mapper.getInputSchema().find<afw::table::Flag>(name + "_" + flag.name).key;
        mapper.addMapping(key);
    }

    _outShapeKey = ShapeResultKey::addFields(mapper.editOutputSchema(), name, "Shape in celestial moments",
                                             FULL_COVARIANCE, afw::table::CoordinateType::CELESTIAL);
    if (_transformPsf) {
        _outPsfShapeKey = afw::table::QuadrupoleKey::addFields(mapper.editOutputSchema(), name + "_psf",
                                                               "PSF shape in celestial moments",
                                                               afw::table::CoordinateType::CELESTIAL);
    }
}

void SdssShapeTransform::operator()(afw::table::SourceCatalog const &inputCatalog,
                                    afw::table::BaseCatalog &outputCatalog, afw::geom::SkyWcs const &wcs,
                                    afw::image::PhotoCalib const &photoCalib) const {
    // The instFlux and cetroid transforms will check that the catalog lengths
    // match and throw if not, so we don't repeat the test here.
    _instFluxTransform(inputCatalog, outputCatalog, wcs, photoCalib);
    _centroidTransform(inputCatalog, outputCatalog, wcs, photoCalib);

    CentroidResultKey centroidKey(inputCatalog.getSchema()[_name]);
    ShapeResultKey inShapeKey(inputCatalog.getSchema()[_name]);
    afw::table::QuadrupoleKey inPsfShapeKey;
    if (_transformPsf) {
        inPsfShapeKey = afw::table::QuadrupoleKey(
                inputCatalog.getSchema()[inputCatalog.getSchema().join(_name, "psf")]);
    }

    afw::table::SourceCatalog::const_iterator inSrc = inputCatalog.begin();
    afw::table::BaseCatalog::iterator outSrc = outputCatalog.begin();
    for (; inSrc != inputCatalog.end(); ++inSrc, ++outSrc) {
        ShapeResult inShape = inShapeKey.get(*inSrc);
        ShapeResult outShape;

        // The transformation from the (x, y) to the (Ra, Dec) basis.
        geom::AffineTransform crdTr =
                wcs.linearizePixelToSky(centroidKey.get(*inSrc).getCentroid(), geom::radians);
        outShape.setShape(inShape.getShape().transform(crdTr.getLinear()));

        // Transformation matrix from pixel to celestial basis.
        ShapeTrMatrix m = makeShapeTransformMatrix(crdTr.getLinear());
        outShape.setShapeErr((m * inShape.getShapeErr().cast<double>() * m.transpose()).cast<ErrElement>());

        _outShapeKey.set(*outSrc, outShape);

        if (_transformPsf) {
            _outPsfShapeKey.set(*outSrc, inPsfShapeKey.get(*inSrc).transform(crdTr.getLinear()));
        }
    }
}

}  // namespace base
}  // namespace meas
}  // namespace lsst
