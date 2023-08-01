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

#include "lsst/meas/base/ScaledApertureFlux.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace meas {
namespace base {

namespace {

using PyFluxAlgorithm = py::class_<ScaledApertureFluxAlgorithm,
                                   SimpleAlgorithm>;
using PyFluxControl = py::class_<ScaledApertureFluxControl>;
using PyFluxTransform =
        py::class_<ScaledApertureFluxTransform, BaseTransform>;

PyFluxControl declareFluxControl(lsst::cpputils::python::WrapperCollection &wrappers) {
    return wrappers.wrapType(PyFluxControl(wrappers.module, "ScaledApertureFluxControl"), [](auto &mod, auto &cls) {
        LSST_DECLARE_CONTROL_FIELD(cls, ScaledApertureFluxControl, scale);
        LSST_DECLARE_CONTROL_FIELD(cls, ScaledApertureFluxControl, shiftKernel);

        cls.def(py::init<>());
    });
}

PyFluxAlgorithm declareFluxAlgorithm(lsst::cpputils::python::WrapperCollection &wrappers) {
    return wrappers.wrapType(PyFluxAlgorithm(wrappers.module, "ScaledApertureFluxAlgorithm"), [](auto &mod, auto &cls) {
        cls.def(py::init<ScaledApertureFluxAlgorithm::Control const &, std::string const &,
                        afw::table::Schema &>(),
                "ctrl"_a, "name"_a, "schema"_a);

        cls.def("measure", &ScaledApertureFluxAlgorithm::measure, "measRecord"_a, "exposure"_a);
        cls.def("fail", &ScaledApertureFluxAlgorithm::fail, "measRecord"_a, "error"_a = nullptr);
    });
}

PyFluxTransform declareFluxTransform(lsst::cpputils::python::WrapperCollection &wrappers) {
    return wrappers.wrapType(PyFluxTransform(wrappers.module, "ScaledApertureFluxTransform"), [](auto &mod, auto &cls) {


        cls.def(py::init<ScaledApertureFluxTransform::Control const &, std::string const &,
                        afw::table::SchemaMapper &>(),
                "ctrl"_a, "name"_a, "mapper"_a);
    });
}

}  // namespace

void wrapScaledApertureFlux(lsst::cpputils::python::WrapperCollection &wrappers) {
    auto clsFluxControl = declareFluxControl(wrappers);
    auto clsFluxAlgorithm = declareFluxAlgorithm(wrappers);
    auto clsFluxTransform = declareFluxTransform(wrappers);

    clsFluxAlgorithm.attr("Control") = clsFluxControl;
    clsFluxTransform.attr("Control") = clsFluxControl;

    python::declareAlgorithm<ScaledApertureFluxAlgorithm, ScaledApertureFluxControl,
            ScaledApertureFluxTransform>(clsFluxAlgorithm, clsFluxControl, clsFluxTransform);
}

}  // namespace base
}  // namespace meas
}  // namespace lsst
