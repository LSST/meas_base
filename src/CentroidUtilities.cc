/*
 * LSST Data Management System
 * Copyright 2008-2014 LSST Corporation.
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
 * see <http://www.lsstcorp.org/LegalNotices/>.
 */

#include "lsst/meas/base/CentroidUtilities.h"
#include "lsst/afw/table/BaseRecord.h"

namespace lsst { namespace meas { namespace base {

CentroidResult::CentroidResult() :
    x(std::numeric_limits<CentroidElement>::quiet_NaN()),
    y(std::numeric_limits<CentroidElement>::quiet_NaN()),
    xSigma(std::numeric_limits<ErrElement>::quiet_NaN()),
    ySigma(std::numeric_limits<ErrElement>::quiet_NaN()),
    x_y_Cov(std::numeric_limits<ErrElement>::quiet_NaN())
{}

Centroid const CentroidResult::getCentroid() const { return Centroid(x, y); }

void CentroidResult::setCentroid(Centroid const & centroid) {
    x = centroid.getX();
    y = centroid.getY();
}

CentroidCov const CentroidResult::getCentroidErr() const {
    CentroidCov m;
    m <<
        xSigma*xSigma, x_y_Cov,
        x_y_Cov, ySigma*ySigma;
    return m;
}

void CentroidResult::setCentroidErr(CentroidCov const & matrix) {
    xSigma = std::sqrt(matrix(0, 0));
    ySigma = std::sqrt(matrix(1, 1));
    x_y_Cov = matrix(0, 1);
}

CentroidResultKey CentroidResultKey::addFields(
    afw::table::Schema & schema,
    std::string const & name,
    std::string const & doc,
    UncertaintyEnum uncertainty
) {
    CentroidResultKey r;
    r._centroid = afw::table::PointKey<CentroidElement>::addFields(
        schema,
        name,
        doc,
        "pixels"
    );
    if (uncertainty != NO_UNCERTAINTY) {
        std::vector< afw::table::Key<ErrElement> > sigma(2);
        std::vector< afw::table::Key<ErrElement> > cov;
        sigma[0] = schema.addField<ErrElement>(
            schema.join(name, "xSigma"), "1-sigma uncertainty on x position", "pixels"
        );
        sigma[1] = schema.addField<ErrElement>(
            schema.join(name, "ySigma"), "1-sigma uncertainty on y position", "pixels"
        );
        if (uncertainty == FULL_COVARIANCE) {
            cov.push_back(
                schema.addField<ErrElement>(
                    schema.join(name, "x_y_Cov"), "uncertainty covariance in x and y", "pixels^2"
                )
            );
        }
        r._centroidErr = afw::table::CovarianceMatrixKey<ErrElement,2>(sigma, cov);
    }
    return r;
}

namespace {

std::vector<std::string> getNameVector() {
    std::vector<std::string> v;
    v.push_back("x");
    v.push_back("y");
    return v;
}

} // anonymous

CentroidResultKey::CentroidResultKey(afw::table::SubSchema const & s) :
    _centroid(s)
{
    static std::vector<std::string> names = getNameVector(); // C++11 TODO: just use initializer list
    try {
        _centroidErr = afw::table::CovarianceMatrixKey<ErrElement,2>(s, names);
    } catch (pex::exceptions::NotFoundError &) {}
}

CentroidResult CentroidResultKey::get(afw::table::BaseRecord const & record) const {
    CentroidResult r;
    r.setCentroid(record.get(_centroid));
    if (_centroidErr.isValid()) {
        r.setCentroidErr(record.get(_centroidErr));
    }
    return r;
}

void CentroidResultKey::set(afw::table::BaseRecord & record, CentroidResult const & value) const {
    record.set(_centroid, value.getCentroid());
    if (_centroidErr.isValid()) {
        record.set(_centroidErr, value.getCentroidErr());
    }
}

}}} // lsst::meas::base
