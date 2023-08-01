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
#include "lsst/cpputils/python.h"

#include <memory>

#include "lsst/pex/config/python.h"
#include "lsst/meas/base/python.h"

#include "lsst/meas/base/SdssCentroid.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace meas {
namespace base {

namespace {

using PyCentroidAlgorithm =
        py::class_<SdssCentroidAlgorithm, SimpleAlgorithm>;
using PyCentroidControl = py::class_<SdssCentroidControl>;
using PyCentroidTransform =
        py::class_<SdssCentroidTransform, BaseTransform>;

PyCentroidControl declareCentroidControl(lsst::cpputils::python::WrapperCollection &wrappers) {
    return wrappers.wrapType(PyCentroidControl(wrappers.module, "SdssCentroidControl"), [](auto &mod, auto &cls) {
        LSST_DECLARE_CONTROL_FIELD(cls, SdssCentroidControl, binmax);
        LSST_DECLARE_CONTROL_FIELD(cls, SdssCentroidControl, peakMin);
        LSST_DECLARE_CONTROL_FIELD(cls, SdssCentroidControl, wfac);
        LSST_DECLARE_CONTROL_FIELD(cls, SdssCentroidControl, doFootprintCheck);
        LSST_DECLARE_CONTROL_FIELD(cls, SdssCentroidControl, maxDistToPeak);

        cls.def(py::init<>());
    });
}

PyCentroidAlgorithm declareCentroidAlgorithm(lsst::cpputils::python::WrapperCollection &wrappers) {
    return wrappers.wrapType(PyCentroidAlgorithm(wrappers.module, "SdssCentroidAlgorithm"), [](auto &mod, auto &cls) {
        cls.attr("FAILURE") = py::cast(SdssCentroidAlgorithm::FAILURE);
        cls.attr("EDGE") = py::cast(SdssCentroidAlgorithm::EDGE);
        cls.attr("NO_SECOND_DERIVATIVE") = py::cast(SdssCentroidAlgorithm::NO_SECOND_DERIVATIVE);
        cls.attr("ALMOST_NO_SECOND_DERIVATIVE") = py::cast(SdssCentroidAlgorithm::ALMOST_NO_SECOND_DERIVATIVE);
        cls.attr("NOT_AT_MAXIMUM") = py::cast(SdssCentroidAlgorithm::NOT_AT_MAXIMUM);

        cls.def(py::init<SdssCentroidAlgorithm::Control const &, std::string const &, afw::table::Schema &>(),
                "ctrl"_a, "name"_a, "schema"_a);

        cls.def("measure", &SdssCentroidAlgorithm::measure, "measRecord"_a, "exposure"_a);
        cls.def("fail", &SdssCentroidAlgorithm::fail, "measRecord"_a, "error"_a = nullptr);
    });
}

PyCentroidTransform declareCentroidTransform(lsst::cpputils::python::WrapperCollection &wrappers) {
    return wrappers.wrapType(PyCentroidTransform(wrappers.module, "SdssCentroidTransform"), [](auto &mod, auto &cls) {
        cls.def(py::init<SdssCentroidTransform::Control const &, std::string const &,
                        afw::table::SchemaMapper &>(),
                "ctrl"_a, "name"_a, "mapper"_a);
    });
}

}  // namespace

void wrapSddsCentroid(lsst::cpputils::python::WrapperCollection &wrappers) {
    auto clsCentroidControl = declareCentroidControl(wrappers);
    auto clsCentroidAlgorithm = declareCentroidAlgorithm(wrappers);
    auto clsCentroidTransform = declareCentroidTransform(wrappers);

    clsCentroidAlgorithm.attr("Control") = clsCentroidControl;
    clsCentroidTransform.attr("Control") = clsCentroidControl;

    python::declareAlgorithm<SdssCentroidAlgorithm, SdssCentroidControl, SdssCentroidTransform>(
            clsCentroidAlgorithm, clsCentroidControl, clsCentroidTransform);
}

}  // namespace base
}  // namespace meas
}  // namespace lsst
