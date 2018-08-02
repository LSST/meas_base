#
# LSST Data Management System
# Copyright 2008-2017 LSST Corporation.
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

import unittest

import numpy as np

import lsst.utils.tests
import lsst.geom
import lsst.meas.base.tests
import lsst.afw.image as afwImage
import lsst.afw.table as afwTable
import lsst.meas.base.applyApCorr as applyApCorr
from lsst.afw.math import ChebyshevBoundedField
from lsst.meas.base.apCorrRegistry import addApCorrName


def initializeSourceCatalog(schema=None, name=None, flux=None, sigma=None, centroid=None):
    fluxName = name + "_flux"
    fluxSigmaName = name + "_fluxSigma"
    fluxKey = schema.find(fluxName).key
    centroidKey = afwTable.Point2DKey(schema["slot_Centroid"])
    sourceCat = afwTable.SourceCatalog(schema)
    source = sourceCat.addNew()
    source.set(fluxKey, flux)
    source.set(fluxSigmaName, sigma)
    source.set(centroidKey, centroid)
    return(sourceCat)


class ApplyApCorrTestCase(lsst.meas.base.tests.AlgorithmTestCase, lsst.utils.tests.TestCase):

    def setUp(self):
        schema = afwTable.SourceTable.makeMinimalSchema()
        name = "test"
        addApCorrName(name)
        schema.addField(name + "_flux", type=np.float64)
        schema.addField(name + "_fluxSigma", type=np.float64)
        schema.addField(name + "_flag", type=np.float64)
        schema.addField(name + "_Centroid_x", type=np.float64)
        schema.addField(name + "_Centroid_y", type=np.float64)
        schema.getAliasMap().set('slot_Centroid', name + '_Centroid')
        self.ap_corr_task = applyApCorr.ApplyApCorrTask(schema=schema)
        self.name = name
        self.schema = schema

    def tearDown(self):
        del self.schema
        del self.ap_corr_task

    def testAddFields(self):
        # Check that the required fields have been added to the schema
        self.assertIn(self.name + "_apCorr", self.schema.getNames())
        self.assertIn(self.name + "_apCorrSigma", self.schema.getNames())
        self.assertIn(self.name + "_flag_apCorr", self.schema.getNames())

    def testSuccessUnflagged(self):
        # Check that the aperture correction flag is set to False if aperture correction was successfully run
        flagName = self.name + "_flag_apCorr"
        flagKey = self.schema.find(flagName).key
        source_test_flux = 5.1
        source_test_centroid = lsst.geom.Point2D(5, 7.1)
        sourceCat = initializeSourceCatalog(schema=self.schema, name=self.name, flux=source_test_flux,
                                            sigma=0, centroid=source_test_centroid)
        fluxName = self.name + "_flux"
        fluxSigmaName = self.name + "_fluxSigma"

        apCorrMap = afwImage.ApCorrMap()
        bbox = lsst.geom.Box2I(lsst.geom.Point2I(0, 0), lsst.geom.ExtentI(10, 10), invert=False)
        coefficients = np.ones((1, 1), dtype=np.float64)
        coefficients_sigma = np.zeros((1, 1), dtype=np.float64)
        apCorrMap[fluxName] = ChebyshevBoundedField(bbox, coefficients)
        apCorrMap[fluxSigmaName] = ChebyshevBoundedField(bbox, coefficients_sigma)
        self.ap_corr_task.run(sourceCat, apCorrMap)
        self.assertFalse(sourceCat[flagKey])

    def testFailureFlagged(self):
        # Check that aperture correction flag is set to True if aperture correction is invalid (negative)
        flagName = self.name + "_flag_apCorr"
        flagKey = self.schema.find(flagName).key
        source_test_flux = 5.2
        source_test_centroid = lsst.geom.Point2D(5, 7.1)
        sourceCat = initializeSourceCatalog(schema=self.schema, name=self.name, flux=source_test_flux,
                                            sigma=0, centroid=source_test_centroid)
        fluxName = self.name + "_flux"
        fluxSigmaName = self.name + "_fluxSigma"

        apCorrMap = afwImage.ApCorrMap()
        bbox = lsst.geom.Box2I(lsst.geom.Point2I(0, 0), lsst.geom.ExtentI(10, 10), invert=False)
        coefficients = -(np.ones((1, 1), dtype=np.float64))
        coefficients_sigma = np.zeros((1, 1), dtype=np.float64)
        apCorrMap[fluxName] = ChebyshevBoundedField(bbox, coefficients)
        apCorrMap[fluxSigmaName] = ChebyshevBoundedField(bbox, coefficients_sigma)
        self.ap_corr_task.run(sourceCat, apCorrMap)
        self.assertTrue(sourceCat[flagKey])

    def testCatFluxUnchanged(self):
        # Pick arbitrary but unique values for the test case
        source_test_flux = 5.3
        source_test_centroid = lsst.geom.Point2D(5, 7.1)
        sourceCat = initializeSourceCatalog(schema=self.schema, name=self.name, flux=source_test_flux,
                                            sigma=0, centroid=source_test_centroid)
        fluxName = self.name + "_flux"
        fluxSigmaName = self.name + "_fluxSigma"
        fluxKey = self.schema.find(fluxName).key

        apCorrMap = afwImage.ApCorrMap()
        bbox = lsst.geom.Box2I(lsst.geom.Point2I(0, 0), lsst.geom.ExtentI(10, 10), invert=False)
        coefficients = np.ones((1, 1), dtype=np.float64)
        coefficients_sigma = np.zeros((1, 1), dtype=np.float64)
        apCorrMap[fluxName] = ChebyshevBoundedField(bbox, coefficients)
        apCorrMap[fluxSigmaName] = ChebyshevBoundedField(bbox, coefficients_sigma)
        self.ap_corr_task.run(sourceCat, apCorrMap)

        self.assertEqual(sourceCat[fluxKey], source_test_flux)

    def testCatFluxHalf(self):
        # Pick arbitrary but unique values for the test case
        source_test_flux = 5.4
        source_test_centroid = lsst.geom.Point2D(5, 7.1)
        sourceCat = initializeSourceCatalog(schema=self.schema, name=self.name, flux=source_test_flux,
                                            sigma=0, centroid=source_test_centroid)
        fluxName = self.name + "_flux"
        fluxSigmaName = self.name + "_fluxSigma"
        fluxKey = self.schema.find(fluxName).key

        apCorrMap = afwImage.ApCorrMap()
        bbox = lsst.geom.Box2I(lsst.geom.Point2I(0, 0), lsst.geom.ExtentI(10, 10), invert=False)
        coefficients = np.ones((1, 1), dtype=np.float64)
        coefficients /= 2.
        coefficients_sigma = np.zeros((1, 1), dtype=np.float64)
        apCorrMap[fluxName] = ChebyshevBoundedField(bbox, coefficients)
        apCorrMap[fluxSigmaName] = ChebyshevBoundedField(bbox, coefficients_sigma)
        self.ap_corr_task.run(sourceCat, apCorrMap)

        self.assertAlmostEqual(sourceCat[fluxKey], source_test_flux / 2)

    def testCatFluxSigma(self):
        """
        Important note! This test will break if UseNaiveFluxSigma = False
        The alternate method significantly overestimates noise, causing this test to fail.
        It is likely that this test will need to be modified if the noise calculation is updated.
        """
        # Pick arbitrary but unique values for the test case
        source_test_flux = 5.5
        source_test_sigma = 0.23
        source_test_centroid = lsst.geom.Point2D(5, 7.3)
        sourceCat = initializeSourceCatalog(schema=self.schema, name=self.name, flux=source_test_flux,
                                            sigma=source_test_sigma, centroid=source_test_centroid)

        fluxName = self.name + "_flux"
        fluxSigmaName = self.name + "_fluxSigma"
        fluxSigmaKey = self.schema.find(fluxSigmaName).key

        apCorrMap = afwImage.ApCorrMap()
        bbox = lsst.geom.Box2I(lsst.geom.Point2I(0, 0), lsst.geom.ExtentI(10, 10), invert=False)
        coefficients = np.ones((1, 1), dtype=np.float64)
        coefficients_sigma = np.ones((1, 1), dtype=np.float64)
        apCorrMap[fluxName] = ChebyshevBoundedField(bbox, coefficients)
        apCorrMap[fluxSigmaName] = ChebyshevBoundedField(bbox, coefficients_sigma)
        self.ap_corr_task.run(sourceCat, apCorrMap)

        self.assertAlmostEqual(sourceCat[fluxSigmaKey], source_test_sigma)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
