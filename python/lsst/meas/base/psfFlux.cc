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
#include "pybind11/stl.h"
#include "lsst/cpputils/python.h"

#include <memory>

#include "lsst/pex/config/python.h"
#include "lsst/meas/base/python.h"

#include "lsst/meas/base/PsfFlux.h"
#include "lsst/meas/base/FluxUtilities.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace meas {
namespace base {

namespace {

using PyFluxAlgorithm = py::class_<PsfFluxAlgorithm, SimpleAlgorithm>;
using PyFluxControl = py::class_<PsfFluxControl>;
using PyFluxTransform = py::class_<PsfFluxTransform, BaseTransform>;

PyFluxControl declareFluxControl(lsst::cpputils::python::WrapperCollection &wrappers) {
     return wrappers.wrapType(PyFluxControl(wrappers.module, "PsfFluxControl"), [](auto &mod, auto &cls) {

         LSST_DECLARE_CONTROL_FIELD(cls, PsfFluxControl, badMaskPlanes);

         cls.def(py::init<>());
     });
}

PyFluxAlgorithm declareFluxAlgorithm(lsst::cpputils::python::WrapperCollection &wrappers) {
    return wrappers.wrapType(PyFluxAlgorithm(wrappers.module, "PsfFluxAlgorithm"), [](auto &mod, auto &cls) {
        cls.attr("FAILURE") = py::cast(PsfFluxAlgorithm::FAILURE);
        cls.attr("NO_GOOD_PIXELS") = py::cast(PsfFluxAlgorithm::NO_GOOD_PIXELS);
        cls.attr("EDGE") = py::cast(PsfFluxAlgorithm::EDGE);

        cls.def(py::init<PsfFluxAlgorithm::Control const &, std::string const &, afw::table::Schema &>(),
                "ctrl"_a, "name"_a, "schema"_a);

        cls.def(py::init<PsfFluxAlgorithm::Control const &, std::string const &, afw::table::Schema &,
                        std::string const &>(),
                "ctrl"_a, "name"_a, "schema"_a, "logName"_a);
    });
}

PyFluxTransform declareFluxTransform(lsst::cpputils::python::WrapperCollection &wrappers) {
    return wrappers.wrapType(PyFluxTransform(wrappers.module, "PsfFluxTransform"), [](auto &mod, auto &cls) {
        cls.def(py::init<PsfFluxTransform::Control const &, std::string const &, afw::table::SchemaMapper &>(),
                "ctrl"_a, "name"_a, "mapper"_a);
    });
}

}  // namespace

void wrapPsfFlux(lsst::cpputils::python::WrapperCollection &wrappers) {
    auto clsFluxControl = declareFluxControl(wrappers);
    auto clsFluxAlgorithm = declareFluxAlgorithm(wrappers);
    auto clsFluxTransform = declareFluxTransform(wrappers);
    clsFluxAlgorithm.attr("Control") = clsFluxControl;
    clsFluxTransform.attr("Control") = clsFluxControl;

    python::declareAlgorithm<PsfFluxAlgorithm, PsfFluxControl, PsfFluxTransform>(
            clsFluxAlgorithm, clsFluxControl, clsFluxTransform);
}

}  // namespace base
}  // namespace meas
}  // namespace lsst
