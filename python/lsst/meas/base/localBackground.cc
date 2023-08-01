/*
 * LSST Data Management System
 * Copyright 2008-2018  AURA/LSST.
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

#include "lsst/meas/base/LocalBackground.h"
#include "lsst/meas/base/FluxUtilities.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace meas {
namespace base {

namespace {

using PyAlgorithm =
        py::class_<LocalBackgroundAlgorithm, SimpleAlgorithm>;
using PyControl = py::class_<LocalBackgroundControl>;
using PyTransform =
        py::class_<LocalBackgroundTransform, BaseTransform>;

PyControl declareControl(lsst::cpputils::python::WrapperCollection &wrappers) {
    return wrappers.wrapType(PyControl(wrappers.module, "LocalBackgroundControl"), [](auto &mod, auto &cls) {
        LSST_DECLARE_CONTROL_FIELD(cls, LocalBackgroundControl, badMaskPlanes);
        LSST_DECLARE_CONTROL_FIELD(cls, LocalBackgroundControl, annulusInner);
        LSST_DECLARE_CONTROL_FIELD(cls, LocalBackgroundControl, annulusOuter);
        LSST_DECLARE_CONTROL_FIELD(cls, LocalBackgroundControl, bgRej);
        LSST_DECLARE_CONTROL_FIELD(cls, LocalBackgroundControl, bgIter);
        cls.def(py::init<>());
    });
}

PyAlgorithm declareAlgorithm(lsst::cpputils::python::WrapperCollection &wrappers) {
    return wrappers.wrapType(PyAlgorithm(wrappers.module, "LocalBackgroundAlgorithm"), [](auto &mod, auto &cls) {
        cls.attr("FAILURE") = py::cast(LocalBackgroundAlgorithm::FAILURE);
        cls.attr("NO_GOOD_PIXELS") = py::cast(LocalBackgroundAlgorithm::NO_GOOD_PIXELS);

        cls.def(py::init<LocalBackgroundAlgorithm::Control const &, std::string const &, afw::table::Schema &>(),
                "ctrl"_a, "name"_a, "schema"_a);

        cls.def(py::init<LocalBackgroundAlgorithm::Control const &, std::string const &, afw::table::Schema &,
                        std::string const &>(),
                "ctrl"_a, "name"_a, "schema"_a, "logName"_a);
    });
}

PyTransform declareTransform(lsst::cpputils::python::WrapperCollection &wrappers) {
    return wrappers.wrapType(PyTransform(wrappers.module, "LocalBackgroundTransform"), [](auto &mod, auto &cls) {
        cls.def(py::init<LocalBackgroundTransform::Control const &, std::string const &,
                        afw::table::SchemaMapper &>(),
                "ctrl"_a, "name"_a, "mapper"_a);
    });
}

}  // namespace

void wrapLocalBackground(lsst::cpputils::python::WrapperCollection &wrappers) {
    auto clsControl = declareControl(wrappers);
    auto clsAlgorithm = declareAlgorithm(wrappers);
    auto clsTransform = declareTransform(wrappers);

    clsAlgorithm.attr("Control") = clsControl;
    clsTransform.attr("Control") = clsControl;

    python::declareAlgorithm<LocalBackgroundAlgorithm, LocalBackgroundControl, LocalBackgroundTransform>(
            clsAlgorithm, clsControl, clsTransform);
}

}  // namespace base
}  // namespace meas
}  // namespace lsst
