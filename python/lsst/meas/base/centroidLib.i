// -*- lsst-c++ -*-
/*
 * LSST Data Management System
 * Copyright 2008-2013 LSST Corporation.
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
 * see <http://www.lsstcorp.org/LegalNotices/>.
 */

%{
#include "lsst/meas/base/SdssCentroid.h"
%}

%include "lsst/meas/base/SdssCentroid.h"

%template(apply) lsst::meas::base::SdssCentroidAlgorithm::apply<float>;
%template(apply) lsst::meas::base::SdssCentroidAlgorithm::apply<double>;
%wrapMeasurementAlgorithm1(lsst::meas::base, SdssCentroidAlgorithm,
                           SdssCentroidControl, FootprintCentroidInput, CentroidComponent)

%include "lsst/meas/base/GaussianCentroid.h"
%template(apply) lsst::meas::base::GaussianCentroidAlgorithm::apply<float>;
%template(apply) lsst::meas::base::GaussianCentroidAlgorithm::apply<double>;
%wrapMeasurementAlgorithm1(lsst::meas::base, GaussianCentroidAlgorithm, GaussianCentroidControl,
                           FootprintCentroidInput, CentroidComponent)

%include "lsst/meas/base/NaiveCentroid.h"
%template(apply) lsst::meas::base::NaiveCentroidAlgorithm::apply<float>;
%template(apply) lsst::meas::base::NaiveCentroidAlgorithm::apply<double>;
%wrapMeasurementAlgorithm1(lsst::meas::base, NaiveCentroidAlgorithm,
                           NaiveCentroidControl, FootprintCentroidInput, CentroidComponent)

%include "lsst/meas/base/SdssCentroid.h"
%template(apply) lsst::meas::base::SdssCentroidAlgorithm::apply<float>;
%template(apply) lsst::meas::base::SdssCentroidAlgorithm::apply<double>;
%wrapMeasurementAlgorithm1(lsst::meas::base, SdssCentroidAlgorithm,
                           SdssCentroidControl, FootprintCentroidInput, CentroidComponent)
