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

#include "lsst/meas/base/GaussianFlux.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace meas {
namespace base {

namespace {

using PyFluxAlgorithm =
        py::class_<GaussianFluxAlgorithm, SimpleAlgorithm>;
using PyFluxControl = py::class_<GaussianFluxControl>;
using PyFluxTransform =
        py::class_<GaussianFluxTransform, BaseTransform>;

PyFluxControl declareFluxControl(lsst::cpputils::python::WrapperCollection &wrappers) {
    return wrappers.wrapType(PyFluxControl(wrappers.module, "GaussianFluxControl"), [](auto &mod, auto &cls) {
        LSST_DECLARE_CONTROL_FIELD(cls, GaussianFluxControl, background);
    });
}

PyFluxAlgorithm declareFluxAlgorithm(lsst::cpputils::python::WrapperCollection &wrappers) {
    return wrappers.wrapType(PyFluxAlgorithm(wrappers.module, "GaussianFluxAlgorithm"), [](auto &mod, auto &cls) {
        cls.def(py::init<GaussianFluxAlgorithm::Control const &, std::string const &, afw::table::Schema &>(),
                "ctrl"_a, "name"_a, "schema"_a);

        cls.attr("FAILURE") = py::cast(GaussianFluxAlgorithm::FAILURE);

        cls.def("measure", &GaussianFluxAlgorithm::measure, "measRecord"_a, "exposure"_a);
        cls.def("fail", &GaussianFluxAlgorithm::fail, "measRecord"_a, "error"_a = nullptr);
    });
}

PyFluxTransform declareFluxTransform(lsst::cpputils::python::WrapperCollection &wrappers) {
    return wrappers.wrapType(PyFluxTransform(wrappers.module, "GaussianFluxTransform"), [](auto &mod, auto &cls) {
        cls.def(py::init<GaussianFluxTransform::Control const &, std::string const &,
                        afw::table::SchemaMapper &>(),
                "ctrl"_a, "name"_a, "mapper"_a);
    });
}

}  // namespace

void wrapGaussianFlux(lsst::cpputils::python::WrapperCollection &wrappers) {
    auto clsFluxControl = declareFluxControl(wrappers);
    auto clsFluxAlgorithm = declareFluxAlgorithm(wrappers);
    auto clsFluxTransform = declareFluxTransform(wrappers);

    clsFluxAlgorithm.attr("Control") = clsFluxControl;
    clsFluxTransform.attr("Control") = clsFluxControl;

    python::declareAlgorithm<GaussianFluxAlgorithm, GaussianFluxControl, GaussianFluxTransform>(
            clsFluxAlgorithm, clsFluxControl, clsFluxTransform);
}

}  // namespace base
}  // namespace meas
}  // namespace lsst
