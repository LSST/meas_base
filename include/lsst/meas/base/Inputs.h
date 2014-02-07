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

#ifndef LSST_MEAS_BASE_Inputs_h_INCLUDED
#define LSST_MEAS_BASE_Inputs_h_INCLUDED

#include <vector>

#include "lsst/afw/geom/Point.h"
#include "lsst/afw/geom/ellipses/Quadrupole.h"
#include "lsst/afw/detection/Footprint.h"
#include "lsst/afw/table/Source.h"

namespace lsst { namespace meas { namespace base {

/**
 *  @brief Empty control object used by algorithm classes that don't have any configuration parameters.
 *
 *  It'd be a bit cleaner in C++ just to not have one, but having a null one makes the Python side
 *  much cleaner.
 */
struct NullControl {};

//@{
/**
 *  These reusable structs represent the inputs for most algorithms.
 *
 *  Each of these can be constructed directly from its constituents, or constructed from a SourceRecord.
 *  The latter lets us use these interchangeably in Python.  Similarly, they all provide a static method
 *  to create a std::vector of Inputs from a SourceCatalog.
 */

/// An Input struct for algorithms that require only a Footprint
struct AlgorithmInput1 {
    typedef std::vector<AlgorithmInput1> Vector;

    PTR(afw::detection::Footprint) footprint;

    explicit AlgorithmInput1(PTR(afw::detection::Footprint) footprint_) : footprint(footprint_) {}

    explicit AlgorithmInput1(afw::table::SourceRecord const & record) : footprint(record.getFootprint()) {}

    static Vector makeVector(afw::table::SourceCatalog const & catalog);

};

/// An Input struct for algorithms that require a position as well as a Footprint
struct AlgorithmInput2 : public AlgorithmInput1 {
    typedef std::vector<AlgorithmInput2> Vector;

    afw::geom::Point2D position;

    explicit AlgorithmInput2(
        PTR(afw::detection::Footprint) footprint_, afw::geom::Point2D const & position_
    ) : AlgorithmInput1(footprint_), position(position_) {}

    AlgorithmInput2(afw::table::SourceRecord const & record) :
        AlgorithmInput1(record), position(record.getCentroid())
    {}

    static Vector makeVector(afw::table::SourceCatalog const & catalog);

};

/// An Input struct for algorithms that require a position and shape as well as a Footprint
struct AlgorithmInput3 : public AlgorithmInput2 {
    typedef std::vector<AlgorithmInput3> Vector;

    afw::geom::ellipses::Quadrupole shape;

    AlgorithmInput3(
        PTR(afw::detection::Footprint) footprint_,
        afw::geom::Point2D const & position_,
        afw::geom::ellipses::Quadrupole const & shape_
    ) : AlgorithmInput2(footprint_, position_), shape(shape_) {}

    explicit AlgorithmInput3(afw::table::SourceRecord const & record) :
        AlgorithmInput2(record), shape(record.getShape())
    {}

    static Vector makeVector(afw::table::SourceCatalog const & catalog);

};

//@}

}}} // lsst::meas::base

#endif // !LSST_MEAS_BASE_Inputs_h_INCLUDED
