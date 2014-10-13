#!/usr/bin/env python
#
# LSST Data Management System
# Copyright 2008-2013 LSST Corporation.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#

import math
import os
from lsst.afw.table import Schema,SchemaMapper,SourceCatalog,SourceTable
from lsst.meas.base.sfm import SingleFramePluginConfig, SingleFramePlugin, SingleFrameMeasurementTask
from lsst.meas.base.base import *
from lsst.meas.base.tests import *
import unittest
import lsst.utils.tests
import numpy

numpy.random.seed(1234)


DATA_DIR = os.path.join(os.environ["MEAS_BASE_DIR"], "tests")

class SFMTestCase(lsst.utils.tests.TestCase):

    def testAlgorithm(self):

        srccat, bbox = MakeTestData.makeCatalog()
        exposure = MakeTestData.makeEmptyExposure(bbox)
        MakeTestData.fillImages(srccat, exposure)
       
        sfm_config = lsst.meas.base.sfm.SingleFrameMeasurementConfig()
        mapper = SchemaMapper(srccat.getSchema())
        mapper.addMinimalSchema(srccat.getSchema())
        outschema = mapper.getOutputSchema()

        #  Basic test of SincFlux algorithm, no C++ slots
        sfm_config.plugins = ["base_PeakCentroid", "base_SincFlux"]
        sfm_config.slots.centroid = "base_PeakCentroid"
        sfm_config.slots.shape = None
        sfm_config.slots.psfFlux = None
        sfm_config.slots.modelFlux = None
        sfm_config.slots.apFlux = None
        sfm_config.slots.instFlux = "base_SincFlux"
        sfm_config.plugins["base_SincFlux"].radius1 = 0.0
        sfm_config.plugins["base_SincFlux"].radius2 = 16.0
        task = SingleFrameMeasurementTask(outschema, config=sfm_config)
        measCat = SourceCatalog(outschema)
        measCat.extend(srccat, mapper=mapper)
        # now run the SFM task with the test plugin
        task.run(measCat, exposure)
        for i in range(len(measCat)):
            record = measCat[i]
            srcRec = srccat[i]
            # check all the flags
            self.assertFalse(record.get("base_SincFlux_flag"))
            # check the slots
            centroid = record.getCentroid() 
            flux = record.get("base_SincFlux_flux")
            fluxerr = record.get("base_SincFlux_fluxSigma")
            truthFlux = srcRec.get("truth_flux")
            # if a star, see if the flux measured is decent
            if srcRec.get("truth_isStar"):
                self.assertClose(truthFlux, flux, atol=None, rtol=.1)
            self.assertEqual(record.getInstFlux(), flux) 
            self.assertEqual(record.getInstFluxErr(), fluxerr) 


def suite():
    """Returns a suite containing all the test cases in this module."""

    lsst.utils.tests.init()

    suites = []
    suites += unittest.makeSuite(SFMTestCase)
    suites += unittest.makeSuite(lsst.utils.tests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(shouldExit=False):
    """Run the tests"""
    lsst.utils.tests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)
