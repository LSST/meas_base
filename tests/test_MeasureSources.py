# This file is part of meas_base.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import math
import unittest

import numpy as np

import lsst.pex.exceptions
import lsst.daf.base as dafBase
import lsst.geom
import lsst.afw.detection as afwDetection
import lsst.afw.math as afwMath
import lsst.afw.geom as afwGeom
import lsst.afw.table as afwTable
import lsst.afw.image as afwImage
import lsst.meas.base as measBase
import lsst.utils.tests

FwhmPerSigma = 2*math.sqrt(2*math.log(2))  # FWHM for an N(0, 1) Gaussian


def makePluginAndCat(alg, name, control, metadata=False, centroid=None):
    schema = afwTable.SourceTable.makeMinimalSchema()
    if centroid:
        schema.addField(centroid + "_x", type=np.float64)
        schema.addField(centroid + "_y", type=np.float64)
        schema.addField(centroid + "_flag", type='Flag')
        schema.getAliasMap().set("slot_Centroid", centroid)
    if metadata:
        plugin = alg(control, name, schema, dafBase.PropertySet())
    else:
        plugin = alg(control, name, schema)
    cat = afwTable.SourceCatalog(schema)
    return plugin, cat


class MeasureSourcesTestCase(lsst.utils.tests.TestCase):

    def testCircularApertureMeasure(self):
        mi = afwImage.MaskedImageF(lsst.geom.ExtentI(100, 200))
        mi.set(10)
        #
        # Create our measuring engine
        #

        radii = (1.0, 5.0, 10.0)  # radii to use

        control = measBase.ApertureFluxControl()
        control.radii = radii

        exp = afwImage.makeExposure(mi)
        x0, y0 = 1234, 5678
        exp.setXY0(lsst.geom.Point2I(x0, y0))

        plugin, cat = makePluginAndCat(measBase.CircularApertureFluxAlgorithm,
                                       "test", control, True, centroid="centroid")
        source = cat.makeRecord()
        source.set("centroid_x", 30+x0)
        source.set("centroid_y", 50+y0)
        plugin.measure(source, exp)

        for r in radii:
            currentFlux = source.get("%s_instFlux" %
                                     measBase.CircularApertureFluxAlgorithm.makeFieldPrefix("test", r))
            self.assertAlmostEqual(10.0*math.pi*r*r/currentFlux, 1.0, places=4)

    def testPeakLikelihoodFlux(self):
        """Test measurement with PeakLikelihoodFlux.

        Notes
        -----
        This test makes and measures a series of exposures containing just one
        star, approximately centered.
        """

        bbox = lsst.geom.Box2I(lsst.geom.Point2I(0, 0), lsst.geom.Extent2I(100, 101))
        kernelWidth = 35
        var = 100
        fwhm = 3.0
        sigma = fwhm/FwhmPerSigma
        convolutionControl = afwMath.ConvolutionControl()
        psf = afwDetection.GaussianPsf(kernelWidth, kernelWidth, sigma)
        psfKernel = psf.getLocalKernel(psf.getAveragePosition())
        psfImage = psf.computeKernelImage(psf.getAveragePosition())
        sumPsfSq = np.sum(psfImage.array**2)
        psfSqArr = psfImage.array**2

        for instFlux in (1000, 10000):
            ctrInd = lsst.geom.Point2I(50, 51)
            ctrPos = lsst.geom.Point2D(ctrInd)

            kernelBBox = psfImage.getBBox()
            kernelBBox.shift(lsst.geom.Extent2I(ctrInd))

            # compute predicted instFlux error
            unshMImage = makeFakeImage(bbox, [ctrPos], [instFlux], fwhm, var)

            # filter image by PSF
            unshFiltMImage = afwImage.MaskedImageF(unshMImage.getBBox())
            afwMath.convolve(unshFiltMImage, unshMImage, psfKernel, convolutionControl)

            # compute predicted instFlux = value of image at peak / sum(PSF^2)
            # this is a sanity check of the algorithm, as much as anything
            predFlux = unshFiltMImage.image[ctrInd, afwImage.LOCAL] / sumPsfSq
            self.assertLess(abs(instFlux - predFlux), instFlux * 0.01)

            # compute predicted instFlux error based on filtered pixels
            # = sqrt(value of filtered variance at peak / sum(PSF^2)^2)
            predFluxErr = math.sqrt(unshFiltMImage.variance[ctrInd, afwImage.LOCAL]) / sumPsfSq

            # compute predicted instFlux error based on unfiltered pixels
            # = sqrt(sum(unfiltered variance * PSF^2)) / sum(PSF^2)
            # and compare to that derived from filtered pixels;
            # again, this is a test of the algorithm
            varView = afwImage.ImageF(unshMImage.variance, kernelBBox)
            varArr = varView.array
            unfiltPredFluxErr = math.sqrt(np.sum(varArr*psfSqArr)) / sumPsfSq
            self.assertLess(abs(unfiltPredFluxErr - predFluxErr), predFluxErr * 0.01)

            for fracOffset in (lsst.geom.Extent2D(0, 0), lsst.geom.Extent2D(0.2, -0.3)):
                adjCenter = ctrPos + fracOffset
                if fracOffset == lsst.geom.Extent2D(0, 0):
                    maskedImage = unshMImage
                    filteredImage = unshFiltMImage
                else:
                    maskedImage = makeFakeImage(bbox, [adjCenter], [instFlux], fwhm, var)
                    # filter image by PSF
                    filteredImage = afwImage.MaskedImageF(maskedImage.getBBox())
                    afwMath.convolve(filteredImage, maskedImage, psfKernel, convolutionControl)

                exp = afwImage.makeExposure(filteredImage)
                exp.setPsf(psf)
                control = measBase.PeakLikelihoodFluxControl()
                plugin, cat = makePluginAndCat(measBase.PeakLikelihoodFluxAlgorithm, "test",
                                               control, centroid="centroid")
                source = cat.makeRecord()
                source.set("centroid_x", adjCenter.getX())
                source.set("centroid_y", adjCenter.getY())
                plugin.measure(source, exp)
                measFlux = source.get("test_instFlux")
                measFluxErr = source.get("test_instFluxErr")
                self.assertLess(abs(measFlux - instFlux), instFlux * 0.003)

                self.assertLess(abs(measFluxErr - predFluxErr), predFluxErr * 0.2)

                # try nearby points and verify that the instFlux is smaller;
                # this checks that the sub-pixel shift is performed in the
                # correct direction
                for dx in (-0.2, 0, 0.2):
                    for dy in (-0.2, 0, 0.2):
                        if dx == dy == 0:
                            continue
                        offsetCtr = lsst.geom.Point2D(adjCenter[0] + dx, adjCenter[1] + dy)
                        source = cat.makeRecord()
                        source.set("centroid_x", offsetCtr.getX())
                        source.set("centroid_y", offsetCtr.getY())
                        plugin.measure(source, exp)
                        self.assertLess(source.get("test_instFlux"), measFlux)

        # source so near edge of image that PSF does not overlap exposure
        # should result in failure
        for edgePos in (
            (1, 50),
            (50, 1),
            (50, bbox.getHeight() - 1),
            (bbox.getWidth() - 1, 50),
        ):
            source = cat.makeRecord()
            source.set("centroid_x", edgePos[0])
            source.set("centroid_y", edgePos[1])
            with self.assertRaises(lsst.pex.exceptions.RangeError):
                plugin.measure(source, exp)

        # no PSF should result in failure: flags set
        noPsfExposure = afwImage.ExposureF(filteredImage)
        source = cat.makeRecord()
        source.set("centroid_x", edgePos[0])
        source.set("centroid_y", edgePos[1])
        with self.assertRaises(lsst.pex.exceptions.InvalidParameterError):
            plugin.measure(source, noPsfExposure)

    def testPixelFlags(self):
        width, height = 100, 100
        mi = afwImage.MaskedImageF(width, height)
        exp = afwImage.makeExposure(mi)
        mi.image.set(0)
        mask = mi.mask
        sat = mask.getPlaneBitMask('SAT')
        interp = mask.getPlaneBitMask('INTRP')
        edge = mask.getPlaneBitMask('EDGE')
        bad = mask.getPlaneBitMask('BAD')
        nodata = mask.getPlaneBitMask('NO_DATA')
        mask.addMaskPlane('CLIPPED')
        clipped = mask.getPlaneBitMask('CLIPPED')
        mask.set(0)
        mask[20, 20, afwImage.LOCAL] = sat
        mask[60, 60, afwImage.LOCAL] = interp
        mask[40, 20, afwImage.LOCAL] = bad
        mask[20, 80, afwImage.LOCAL] = nodata
        mask[30, 30, afwImage.LOCAL] = clipped
        mask.Factory(mask, lsst.geom.Box2I(lsst.geom.Point2I(0, 0), lsst.geom.Extent2I(3, height))).set(edge)
        x0, y0 = 1234, 5678
        exp.setXY0(lsst.geom.Point2I(x0, y0))
        control = measBase.PixelFlagsControl()
        # Change the configuration of control to test for clipped mask
        control.masksFpAnywhere = ['CLIPPED']
        plugin, cat = makePluginAndCat(measBase.PixelFlagsAlgorithm, "test", control, centroid="centroid")
        allFlags = [
            "",
            "edge",
            "interpolated",
            "interpolatedCenter",
            "saturated",
            "saturatedCenter",
            "cr",
            "crCenter",
            "bad",
            "clipped",
        ]
        for x, y, setFlags in [(1, 50, ['edge']),
                               (40, 20, ['bad']),
                               (20, 20, ['saturatedCenter',
                                         'saturated']),
                               (20, 22, ['saturated']),
                               (60, 60, ['interpolatedCenter',
                                         'interpolated']),
                               (60, 62, ['interpolated']),
                               (20, 80, ['nodata']),
                               (30, 30, ['clipped']),
                               ]:
            spans = afwGeom.SpanSet.fromShape(5).shiftedBy(x + x0,
                                                           y + y0)
            foot = afwDetection.Footprint(spans)
            source = cat.makeRecord()
            source.setFootprint(foot)
            source.set("centroid_x", x+x0)
            source.set("centroid_y", y+y0)
            plugin.measure(source, exp)
            for flag in allFlags[1:]:
                value = source.get("test_flag_" + flag)
                if flag in setFlags:
                    self.assertTrue(value, "Flag %s should be set for %f,%f" % (flag, x, y))
                else:
                    self.assertFalse(value, "Flag %s should not be set for %f,%f" % (flag, x, y))

        # the new code which grabs the center of a record throws when a NaN is
        # set in the centroid slot and the algorithm attempts to get the
        # default center position
        source = cat.makeRecord()
        source.set("centroid_x", float("NAN"))
        source.set("centroid_y", 40)
        source.set("centroid_flag", True)
        tmpSpanSet = afwGeom.SpanSet.fromShape(5).shiftedBy(x + x0,
                                                            y + y0)
        source.setFootprint(afwDetection.Footprint(tmpSpanSet))
        with self.assertRaises(lsst.pex.exceptions.RuntimeError):
            plugin.measure(source, exp)

        # Test that if there is no center and centroider that the object
        # should look at the footprint
        plugin, cat = makePluginAndCat(measBase.PixelFlagsAlgorithm, "test", control)
        # The first test should raise exception because there is no footprint
        source = cat.makeRecord()
        with self.assertRaises(lsst.pex.exceptions.RuntimeError):
            plugin.measure(source, exp)
        # The second test will raise an error because no peaks are present
        tmpSpanSet2 = afwGeom.SpanSet.fromShape(5).shiftedBy(x + x0,
                                                             y + y0)
        source.setFootprint(afwDetection.Footprint(tmpSpanSet2))
        with self.assertRaises(lsst.pex.exceptions.RuntimeError):
            plugin.measure(source, exp)
        # The final test should pass because it detects a peak, we are reusing
        # the location of the clipped bit in the mask plane, so we will check
        # first that it is False, then True
        source.getFootprint().addPeak(x+x0, y+y0, 100)
        self.assertFalse(source.get("test_flag_clipped"), "The clipped flag should be set False")
        plugin.measure(source, exp)
        self.assertTrue(source.get("test_flag_clipped"), "The clipped flag should be set True")


def addStar(image, center, instFlux, fwhm):
    """Add a perfect single Gaussian star to an image.

    Parameters
    ----------
    image : `lsst.afw.image.ImageF`
        Image to which the star will be added.
    center : `list` or `tuple` of `float`, length 2
        Position of the center of the star on the image.
    instFlux : `float`
        instFlux of the Gaussian star, in counts.
    fwhm : `float`
        FWHM of the Gaussian star, in pixels.

    Notes
    -----
    Uses Python to iterate over all pixels (because there is no C++
    function that computes a Gaussian offset by a non-integral amount).
    """
    sigma = fwhm/FwhmPerSigma
    func = afwMath.GaussianFunction2D(sigma, sigma, 0)
    starImage = afwImage.ImageF(image.getBBox())
    # The instFlux in the region of the image will not be exactly the desired
    # instFlux because the Gaussian does not extend to infinity, so keep track
    # of the actual instFlux and correct for it
    actFlux = 0
    # No function exists that has a fractional x and y offset, so set the
    # image the slow way
    for i in range(image.getWidth()):
        x = center[0] - i
        for j in range(image.getHeight()):
            y = center[1] - j
            pixVal = instFlux * func(x, y)
            actFlux += pixVal
            starImage[i, j, afwImage.LOCAL] += pixVal
    starImage *= instFlux / actFlux

    image += starImage


def makeFakeImage(bbox, centerList, instFluxList, fwhm, var):
    """Make a fake image containing a set of stars with variance = image + var.

    Paramters
    ---------
    bbox : `lsst.afw.image.Box2I`
        Bounding box for image.
    centerList : iterable of pairs of `float`
        list of positions of center of star on image.
    instFluxList : `list` of `float`
        instFlux of each star, in counts.
    fwhm : `float`
        FWHM of Gaussian star, in pixels.
    var : `float`
        Value of variance plane, in counts.

    Returns
    -------
    maskedImage : `lsst.afw.image.MaskedImageF`
        Resulting fake image.

    Notes
    -----
    It is trivial to add Poisson noise, which would be more accurate, but
    hard to make a unit test that can reliably determine whether such an
    image passes a test.
    """
    if len(centerList) != len(instFluxList):
        raise RuntimeError("len(centerList) != len(instFluxList)")
    maskedImage = afwImage.MaskedImageF(bbox)
    image = maskedImage.image
    for center, instFlux in zip(centerList, instFluxList):
        addStar(image, center=center, instFlux=instFlux, fwhm=fwhm)
    variance = maskedImage.variance
    variance[:] = image
    variance += var
    return maskedImage


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
