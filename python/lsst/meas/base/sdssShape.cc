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

#include "lsst/afw/table/FunctorKey.h"
#include "lsst/meas/base/ShapeUtilities.h"
#include "lsst/meas/base/SdssShape.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace meas {
namespace base {

namespace {

using PyShapeControl = py::class_<SdssShapeControl>;
// TODO decide if we need to mention afw::table::FunctorKey<SdssShapeResult>
// and if so, wrap it; if not, document that here
using PyShapeResultKey = py::class_<SdssShapeResultKey>;
using PyShapeResult = py::class_<SdssShapeResult, ShapeResult,
                                 CentroidResult, FluxResult>;
using PyShapeAlgorithm = py::class_<SdssShapeAlgorithm, SimpleAlgorithm>;
using PyShapeTransform = py::class_<SdssShapeTransform, BaseTransform>;

PyShapeControl declareShapeControl(lsst::cpputils::python::WrapperCollection &wrappers) {
    return wrappers.wrapType(PyShapeControl(wrappers.module, "SdssShapeControl"), [](auto &mod, auto &cls) {
        
        LSST_DECLARE_CONTROL_FIELD(cls, SdssShapeControl, background);
        LSST_DECLARE_CONTROL_FIELD(cls, SdssShapeControl, maxIter);
        LSST_DECLARE_CONTROL_FIELD(cls, SdssShapeControl, maxShift);
        LSST_DECLARE_CONTROL_FIELD(cls, SdssShapeControl, tol1);
        LSST_DECLARE_CONTROL_FIELD(cls, SdssShapeControl, tol2);
        LSST_DECLARE_CONTROL_FIELD(cls, SdssShapeControl, doMeasurePsf);
        
cls.def(py::init<>());
    });
}

void declareShapeResultKey(lsst::cpputils::python::WrapperCollection &wrappers) {
    wrappers.wrapType(PyShapeResultKey(wrappers.module, "SdssShapeResultKey"), [](auto &mod, auto &cls) {
        // TODO decide whether to wrap default constructor and do it or document why not
        cls.def(py::init<afw::table::SubSchema const &>(), "subSchema"_a);

        cls.def_static("addFields", &FluxResultKey::addFields, "schema"_a, "name"_a, "doMeasurePsf"_a);

        cls.def("__eq__", &SdssShapeResultKey::operator==, py::is_operator());
        cls.def("__ne__", &SdssShapeResultKey::operator!=, py::is_operator());

        cls.def("get", &SdssShapeResultKey::get, "record"_a);
        cls.def("set", &SdssShapeResultKey::set, "record"_a, "value"_a);
        cls.def("getPsfShape", &SdssShapeResultKey::getPsfShape, "record"_a);
        cls.def("setPsfShape", &SdssShapeResultKey::setPsfShape, "record"_a, "value"_a);
        cls.def("isValid", &SdssShapeResultKey::isValid);
        cls.def("getFlagHandler", &SdssShapeResultKey::getFlagHandler);
    });
}

template <typename ImageT>
void declareComputeMethods(PyShapeAlgorithm &cls) {
    cls.def_static("computeAdaptiveMoments",(SdssShapeResult(*)(ImageT const &, geom::Point2D const &, bool, SdssShapeControl const &)) &
                    SdssShapeAlgorithm::computeAdaptiveMoments,
            "image"_a, "position"_a, "negative"_a = false, "ctrl"_a = SdssShapeControl());
    cls.def_static("computeFixedMomentsFlux",(FluxResult(*)(ImageT const &, afw::geom::ellipses::Quadrupole const &, geom::Point2D const &)) &
                    SdssShapeAlgorithm::computeFixedMomentsFlux,
            "image"_a, "shape"_a, "position"_a);
}

PyShapeAlgorithm declareShapeAlgorithm(lsst::cpputils::python::WrapperCollection &wrappers) {
    return wrappers.wrapType(PyShapeAlgorithm(wrappers.module, "SdssShapeAlgorithm"), [](auto &mod, auto &cls) {
        cls.attr("FAILURE") = py::cast(SdssShapeAlgorithm::FAILURE);
        cls.attr("UNWEIGHTED_BAD") = py::cast(SdssShapeAlgorithm::UNWEIGHTED_BAD);
        cls.attr("UNWEIGHTED") = py::cast(SdssShapeAlgorithm::UNWEIGHTED);
        cls.attr("SHIFT") = py::cast(SdssShapeAlgorithm::SHIFT);
        cls.attr("MAXITER") = py::cast(SdssShapeAlgorithm::MAXITER);
        cls.attr("PSF_SHAPE_BAD") = py::cast(SdssShapeAlgorithm::PSF_SHAPE_BAD);

        cls.def(py::init<SdssShapeAlgorithm::Control const &, std::string const &, afw::table::Schema &>(),
                "ctrl"_a, "name"_a, "schema"_a);

        declareComputeMethods<afw::image::Image<int>>(cls);
        declareComputeMethods<afw::image::Image<float>>(cls);
        declareComputeMethods<afw::image::Image<double>>(cls);
        declareComputeMethods<afw::image::MaskedImage<int>>(cls);
        declareComputeMethods<afw::image::MaskedImage<float>>(cls);
        declareComputeMethods<afw::image::MaskedImage<double>>(cls);

        cls.def("measure", &SdssShapeAlgorithm::measure, "measRecord"_a, "exposure"_a);
        cls.def("fail", &SdssShapeAlgorithm::fail, "measRecord"_a, "error"_a = nullptr);
    });
}

void declareShapeResult(lsst::cpputils::python::WrapperCollection &wrappers) {
    wrappers.wrapType(PyShapeResult(wrappers.module, "SdssShapeResult"), [](auto &mod, auto &cls) {
        cls.def(py::init<>());
        cls.def_readwrite("instFlux_xx_Cov", &SdssShapeResult::instFlux_xx_Cov);
        cls.def_readwrite("instFlux_yy_Cov", &SdssShapeResult::instFlux_yy_Cov);
        cls.def_readwrite("instFlux_xy_Cov", &SdssShapeResult::instFlux_xy_Cov);
        cls.def_readwrite("flags", &SdssShapeResult::flags);

        // TODO this method says it's a workaround for Swig which doesn't understand std::bitset
        cls.def("getFlag", (bool (SdssShapeResult::*)(unsigned int) const) &SdssShapeResult::getFlag,
                "index"_a);
        cls.def("getFlag", (bool (SdssShapeResult::*)(std::string const &name) const) &SdssShapeResult::getFlag,
                "name"_a);
    });
}

PyShapeTransform declareShapeTransform(lsst::cpputils::python::WrapperCollection &wrappers) {
    return wrappers.wrapType(PyShapeTransform(wrappers.module, "SdssShapeTransform"), [](auto &mod, auto &cls) {
        cls.def(py::init<SdssShapeTransform::Control const &, std::string const &, afw::table::SchemaMapper &>(),
                "ctrl"_a, "name"_a, "mapper"_a);

        cls.def("__call__", &SdssShapeTransform::operator(), "inputCatalog"_a, "outputCatalog"_a, "wcs"_a,
                "photoCalib"_a);
    });
}

}  // namespace

void wrapSsdsShape(lsst::cpputils::python::WrapperCollection &wrappers) {
    auto clsShapeControl = declareShapeControl(wrappers);
    declareShapeResultKey(wrappers);
    auto clsShapeAlgorithm = declareShapeAlgorithm(wrappers);
    declareShapeResult(wrappers);
    auto clsShapeTransform = declareShapeTransform(wrappers);

    clsShapeAlgorithm.attr("Control") = clsShapeControl;
    clsShapeTransform.attr("Control") = clsShapeControl;

    python::declareAlgorithm<SdssShapeAlgorithm, SdssShapeControl, SdssShapeTransform>(
            clsShapeAlgorithm, clsShapeControl, clsShapeTransform);
}

}  // namespace base
}  // namespace meas
}  // namespace lsst
