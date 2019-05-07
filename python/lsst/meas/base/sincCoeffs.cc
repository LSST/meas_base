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

#include "lsst/meas/base/SincCoeffs.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace meas {
namespace base {

namespace {

template <typename T>
void declareSincCoeffs(py::module& mod, std::string const& suffix) {
    py::class_<SincCoeffs<T>> cls(mod, ("SincCoeffs" + suffix).c_str());

    cls.def_static("cache", &SincCoeffs<T>::cache, "rInner"_a, "rOuter"_a);
    cls.def_static("get", &SincCoeffs<T>::get, "outerEllipse"_a, "innerRadiusFactor"_a);

    cls.attr("DISABLED_AT_COMPILE_TIME") =
#ifdef LSST_DISABLE_SINC_PHOTOMETRY
        true;
#else
        false;
#endif
}

}  // namespace

PYBIND11_MODULE(sincCoeffs, mod) {
    py::module::import("lsst.afw.geom");
    py::module::import("lsst.afw.image");

    declareSincCoeffs<float>(mod, "F");
    declareSincCoeffs<double>(mod, "D");
}

}  // namespace base
}  // namespace meas
}  // namespace lsst
