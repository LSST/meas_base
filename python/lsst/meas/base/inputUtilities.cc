/*
 * LSST Data Management System
 * Copyright 2008-2017  AURA/LSST.
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

#include "pybind11/pybind11.h"

#include "lsst/afw/table/Source.h"
#include "lsst/meas/base/InputUtilities.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace meas {
namespace base {

PYBIND11_PLUGIN(inputUtilities) {
    py::module mod("inputUtilities");

    /* Module level */
    py::class_<SafeCentroidExtractor> clsSafeCentroidExtractor(mod, "SafeCentroidExtractor");

    /* Constructors */
    clsSafeCentroidExtractor.def(py::init<afw::table::Schema &,
                                          std::string const &,
                                          bool>(),
                                 "schema"_a, "name"_a, "isCentroider"_a=false);

    /* Operators */
    clsSafeCentroidExtractor.def("__call__", [](SafeCentroidExtractor const & self,
                           afw::table::SourceRecord & record,
                           FlagHandler const & flags) {
            return self(record, flags);
        }, "record"_a, "flags"_a);

    return mod.ptr();
}

}}}     // lsst::meas::base
