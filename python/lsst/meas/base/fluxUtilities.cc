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

#include "lsst/meas/base/ApertureFlux.h"

namespace py = pybind11;
using namespace py::literals;

namespace lsst {
namespace meas {
namespace base {

namespace {

using PyFluxResult = py::class_<FluxResult>;
using PyFluxResultKey = py::class_<FluxResultKey>;
using PyMagResult = py::class_<MagResult>;
using PyMagResultKey = py::class_<MagResultKey>;

void declareFluxResult(lsst::cpputils::python::WrapperCollection &wrappers) {
    wrappers.wrapType(PyFluxResult(wrappers.module, "FluxResult"), [](auto &mod, auto &cls) {
        cls.def_readwrite("instFlux", &FluxResult::instFlux);
        cls.def_readwrite("instFluxErr", &FluxResult::instFluxErr);
    });
}

void declareFluxResultKey(lsst::cpputils::python::WrapperCollection &wrappers) {
    wrappers.wrapType(PyFluxResultKey(wrappers.module, "FluxResultKey"), [](auto &mod, auto &cls) {
        cls.def(py::init<>());
        cls.def(py::init<afw::table::Key<meas::base::Flux> const &,
                        afw::table::Key<meas::base::FluxErrElement> const &>(),
                "instFlux"_a, "instFluxErr"_a);
        cls.def(py::init<afw::table::SubSchema const &>());

        cls.def("__eq__", &FluxResultKey::operator==, py::is_operator());
        cls.def("__ne__", &FluxResultKey::operator!=, py::is_operator());

        cls.def("get", &FluxResultKey::get);
        cls.def("set", &FluxResultKey::set);
        cls.def_static("addFields", &FluxResultKey::addFields, "schema"_a, "name"_a, "doc"_a);
        cls.def("isValid", &FluxResultKey::isValid);
        cls.def("getInstFlux", &FluxResultKey::getInstFlux);
        cls.def("getInstFluxErr", &FluxResultKey::getInstFluxErr);
    });
}

void declareMagResult(lsst::cpputils::python::WrapperCollection &wrappers) {
    wrappers.wrapType(PyMagResult(wrappers.module, "MagResult"), [](auto &mod, auto &cls) {
        cls.def_readwrite("mag", &MagResult::mag);
        cls.def_readwrite("magErr", &MagResult::magErr);
    });
}

void declareMagResultKey(lsst::cpputils::python::WrapperCollection &wrappers) {
    wrappers.wrapType(PyMagResultKey(wrappers.module, "MagResultKey"), [](auto &mod, auto &cls) {
        cls.def(py::init<>());
        cls.def(py::init<afw::table::SubSchema const &>());

        cls.def("get", &MagResultKey::get);
        cls.def("set",
                (void (MagResultKey::*)(afw::table::BaseRecord &, MagResult const &) const) &MagResultKey::set);
        cls.def("set", (void (MagResultKey::*)(afw::table::BaseRecord &, afw::image::Measurement const &) const) &
                MagResultKey::set);
        cls.def_static("addFields", &MagResultKey::addFields, "schema"_a, "name"_a);
    });
}
}  // namespace

void wrapFluxUtilities(lsst::cpputils::python::WrapperCollection &wrappers) {
    declareFluxResult(wrappers);
    declareFluxResultKey(wrappers);
    declareMagResult(wrappers);
    declareMagResultKey(wrappers);
}

}  // namespace base
}  // namespace meas
}  // namespace lsst
