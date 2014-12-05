// -*- lsst-c++ -*-
/*
 * LSST Data Management System
 * Copyright 2008-2013 LSST Corporation.
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

#include "ndarray/eigen.h"

#include "lsst/afw/detection/Psf.h"
#include "lsst/afw/geom/Box.h"
#include "lsst/afw/geom/ellipses/Ellipse.h"
#include "lsst/afw/detection/FootprintArray.h"
#include "lsst/afw/detection/FootprintArray.cc"
#include "lsst/afw/table/Source.h"
#include "lsst/meas/base/GaussianFlux.h"
#include "lsst/meas/base/detail/SdssShapeImpl.h"

namespace lsst { namespace meas { namespace base {

GaussianFluxAlgorithm::GaussianFluxAlgorithm(
    Control const & ctrl,
    std::string const & name,
    afw::table::Schema & schema
) : _ctrl(ctrl),
    _fluxResultKey(
        FluxResultKey::addFields(schema, name, "flux from Gaussian Flux algorithm")
    ),
    _centroidExtractor(schema, name),
    _shapeExtractor(schema, name)
{
    static boost::array<FlagDefinition,N_FLAGS> const flagDefs = {{
        {"flag", "general failure flag, set if anything went wrong"}
    }};
    _flagHandler = FlagHandler::addFields(schema, name, flagDefs.begin(), flagDefs.end());
}

void GaussianFluxAlgorithm::measure(
    afw::table::SourceRecord & measRecord,
    afw::image::Exposure<float> const & exposure
) const {
    // get the value from the centroid slot only
    afw::geom::Point2D centroid = _centroidExtractor(measRecord, _flagHandler);
    afw::geom::ellipses::Quadrupole shape = _shapeExtractor(measRecord, _flagHandler);
    FluxResult result;

    //  This code came straight out of the GaussianFlux.apply() in meas_algorithms with few changes
    afw::image::Exposure<float>::MaskedImageT const& mimage = exposure.getMaskedImage();

    detail::SdssShapeImpl sdss(centroid, shape);
    std::pair<double, double> fluxResult
        = detail::getFixedMomentsFlux(mimage, _ctrl.background, centroid.getX()-mimage.getX0(), centroid.getY()-mimage.getY0(), sdss);
    result.flux =  fluxResult.first;
    result.fluxSigma = fluxResult.second;
    measRecord.set(_fluxResultKey, result);
    _flagHandler.setValue(measRecord, FAILURE, false);
}


void GaussianFluxAlgorithm::fail(afw::table::SourceRecord & measRecord, MeasurementError * error) const {
    _flagHandler.handleFailure(measRecord, error);
}

}}} // namespace lsst::meas::base
