#if 0
// -*- lsst-c++ -*-
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

#include "lsst/meas/base/CircularApertureFlux.h"
#include "lsst/meas/base/SincCoeffs.h"

namespace lsst { namespace meas { namespace base {

CircularApertureFluxAlgorithm::CircularApertureFluxAlgorithm(
    Control const & ctrl,
    std::string const & name,
    afw::table::Schema & schema
) : ApertureFluxAlgorithm(ctrl, name, schema)
{
    for (std::size_t i = 0; i < ctrl.radii.size(); ++i) {
        if (ctrl.radii[i] > ctrl.maxSincRadius) break;
        SincCoeffs<float>::cache(0.0, ctrl.radii[i]);
    }
}

void CircularApertureFluxAlgorithm::measure(
    afw::table::SourceRecord & record,
    afw::image::Exposure<float> const & exposure
) const {
    afw::geom::ellipses::Ellipse ellipse(afw::geom::ellipses::Axes(1.0, 1.0, 0.0), record.getCentroid());
    PTR(afw::geom::ellipses::Axes) axes
        = boost::static_pointer_cast<afw::geom::ellipses::Axes>(ellipse.getCorePtr());
    for (std::size_t i = 0; i < _ctrl.radii.size(); ++i) {
        axes->setA(_ctrl.radii[i]);
        axes->setB(_ctrl.radii[i]);
        Result result = computeFlux(exposure.getMaskedImage(), ellipse, _ctrl);
        copyResultToRecord(result, record, i);
    }
}

}}} // namespace lsst::meas::base
#endif
