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
#include "pybind11/eigen.h"
#include "lsst/cpputils/python.h"

#include <memory>

#include "lsst/afw/table/BaseRecord.h"
#include "lsst/meas/base/ShapeUtilities.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace meas {
namespace base {

namespace {

using PyShapeResult = py::class_<ShapeResult>;
using PyShapeResultKey = py::class_<ShapeResultKey>;

void declareShapeResult(lsst::cpputils::python::WrapperCollection &wrappers) {
    wrappers.wrapType(PyShapeResult(wrappers.module, "ShapeResult"), [](auto &mod, auto &cls) {
        cls.def(py::init<>());
        cls.def(py::init<ShapeElement, ShapeElement, ShapeElement, ShapeCov const &>(), "xx"_a, "yy"_a, "xy"_a,
                "matrix"_a);
        cls.def(py::init<ShapeElement, ShapeElement, ShapeElement, ErrElement, ErrElement, ErrElement>(), "xx"_a,
                "yy"_a, "xy"_a, "xxErr"_a, "yyErr"_a, "xyErr"_a);
        cls.def("getShape", &ShapeResult::getShape);
        cls.def("getQuadrupole", &ShapeResult::getQuadrupole);
        cls.def("setShape", &ShapeResult::setShape, "shape"_a);
        cls.def("getShapeErr", &ShapeResult::getShapeErr);
        cls.def("setShapeErr", (void (ShapeResult::*)(ShapeCov const &)) &ShapeResult::setShapeErr, "matrix"_a);
        cls.def("setShapeErr",
                (void (ShapeResult::*)(ErrElement, ErrElement, ErrElement)) &ShapeResult::setShapeErr,
                "xxErr"_a, "yyErr"_a, "xyErr"_a);
        cls.def_readwrite("xx", &ShapeResult::xx);
        cls.def_readwrite("yy", &ShapeResult::yy);
        cls.def_readwrite("xy", &ShapeResult::xy);
        cls.def_readwrite("xxErr", &ShapeResult::xxErr);
        cls.def_readwrite("yyErr", &ShapeResult::yyErr);
        cls.def_readwrite("xyErr", &ShapeResult::xyErr);
        cls.def_readwrite("xx_yy_Cov", &ShapeResult::xx_yy_Cov);
        cls.def_readwrite("xx_xy_Cov", &ShapeResult::xx_xy_Cov);
        cls.def_readwrite("yy_xy_Cov", &ShapeResult::yy_xy_Cov);
    });
}

void declareShapeResultKey(lsst::cpputils::python::WrapperCollection &wrappers) {
    wrappers.wrapType(PyShapeResultKey(wrappers.module, "ShapeResultKey"), [](auto &mod, auto &cls) {
        cls.def_static("addFields", &ShapeResultKey::addFields, "schema"_a, "name"_a, "doc"_a, "uncertainty"_a,
                       "coordType"_a = afw::table::CoordinateType::PIXEL);

        cls.def(py::init<>());
        cls.def(py::init<afw::table::QuadrupoleKey const &,
                        afw::table::CovarianceMatrixKey<ErrElement, 3> const &>(),
                "shape"_a, "shapeErr"_a);
        cls.def(py::init<afw::table::SubSchema const &>(), "subSchema"_a);

        cls.def("__eq__", &ShapeResultKey::operator==, py::is_operator());
        cls.def("__ne__", &ShapeResultKey::operator!=, py::is_operator());

        cls.def("get", &ShapeResultKey::get, "record"_a);
        cls.def("set", &ShapeResultKey::set, "record"_a, "value"_a);
        cls.def("isValid", &ShapeResultKey::isValid);
        cls.def("getShape", &ShapeResultKey::getShape);
        cls.def("getShapeErr", &ShapeResultKey::getShapeErr);
        cls.def("getIxx", &ShapeResultKey::getIxx);
        cls.def("getIyy", &ShapeResultKey::getIyy);
        cls.def("getIxy", &ShapeResultKey::getIxy);
    });
}

}  // namespace

void wrapShapeUtilities(lsst::cpputils::python::WrapperCollection &wrappers) {
    declareShapeResult(wrappers);
    declareShapeResultKey(wrappers);
    wrappers.wrap([](auto &mod) {
        mod.def("makeShapeTransformMatrix", &makeShapeTransformMatrix, "xform"_a);
    });
}

}  // namespace base
}  // namespace meas
}  // namespace lsst
