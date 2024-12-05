# This file is part of meas_base.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
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
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import unittest

import numpy as np

import lsst.geom
from lsst.meas.base.tests import (AlgorithmTestCase, CentroidTransformTestCase,
                                  SingleFramePluginTransformSetupHelper)
import lsst.afw.geom
import lsst.utils.tests

# N.B. Some tests here depend on the noise realization in the test data
# or from the numpy random number generator.
# For the current test data and seed value, they pass, but they may not
# if the test data is regenerated or the seed value changes.  I've marked
# these with an "rng dependent" comment.  In most cases, they test that
# the measured instFlux lies within 2 sigma of the correct value, which we
# should expect to fail sometimes.


class SdssCentroidTestCase(AlgorithmTestCase, lsst.utils.tests.TestCase):

    def setUp(self):
        self.center = lsst.geom.Point2D(50.1, 49.8)
        self.bbox = lsst.geom.Box2I(lsst.geom.Point2I(-20, -30),
                                    lsst.geom.Extent2I(140, 160))
        self.dataset = lsst.meas.base.tests.TestDataset(self.bbox)
        self.dataset.addSource(100000.0, self.center)

    def makeAlgorithm(self, ctrl=None):
        """Construct an algorithm and return both it and its schema.
        """
        if ctrl is None:
            ctrl = lsst.meas.base.SdssCentroidControl()
        schema = lsst.meas.base.tests.TestDataset.makeMinimalSchema()
        algorithm = lsst.meas.base.SdssCentroidAlgorithm(ctrl, "base_SdssCentroid", schema)
        return algorithm, schema

    def testSingleFramePlugin(self):
        """Test calling the algorithm through the plugin interface.
        """
        task = self.makeSingleFrameMeasurementTask("base_SdssCentroid")
        exposure, catalog = self.dataset.realize(10.0, task.schema, randomSeed=0)
        task.run(catalog, exposure)
        record = catalog[0]
        self.assertFalse(record.get("base_SdssCentroid_flag"))
        self.assertFalse(record.get("base_SdssCentroid_flag_edge"))
        self.assertFalse(record.get("base_SdssCentroid_flag_notAtMaximum"))
        self.assertFalse(record.get("base_SdssCentroid_flag_near_edge"))
        self.assertFalse(record.get("base_SdssCentroid_flag_resetToPeak"))
        self.assertFalse(record.get("base_SdssCentroid_flag_badError"))
        self.assertFloatsAlmostEqual(record.get("base_SdssCentroid_x"), record.get("truth_x"), rtol=0.005)
        self.assertFloatsAlmostEqual(record.get("base_SdssCentroid_y"), record.get("truth_y"), rtol=0.005)

    def testMonteCarlo(self):
        """Test an ideal simulation, with no noise.

        Demonstrate that:

        - We get exactly the right answer, and
        - The reported uncertainty agrees with a Monte Carlo test of the noise.
        """
        algorithm, schema = self.makeAlgorithm()
        exposure, catalog = self.dataset.realize(0.0, schema, randomSeed=1)
        record = catalog[0]
        x = record.get("truth_x")
        y = record.get("truth_y")
        instFlux = record.get("truth_instFlux")
        algorithm.measure(record, exposure)
        self.assertFloatsAlmostEqual(record.get("base_SdssCentroid_x"), x, rtol=1E-4)
        self.assertFloatsAlmostEqual(record.get("base_SdssCentroid_y"), y, rtol=1E-4)
        for noise in (0.001, 0.01):
            xList = []
            yList = []
            xErrList = []
            yErrList = []
            nSamples = 1000
            for repeat in range(nSamples):
                # By using ``repeat`` to seed the RNG, we get results which
                # fall within the tolerances defined below. If we allow this
                # test to be truly random, passing becomes RNG-dependent.
                exposure, catalog = self.dataset.realize(noise*instFlux, schema, randomSeed=repeat)
                record = catalog[0]
                algorithm.measure(record, exposure)
                xList.append(record.get("base_SdssCentroid_x"))
                yList.append(record.get("base_SdssCentroid_y"))
                xErrList.append(record.get("base_SdssCentroid_xErr"))
                yErrList.append(record.get("base_SdssCentroid_yErr"))
            xMean = np.mean(xList)
            yMean = np.mean(yList)
            xErrMean = np.mean(xErrList)
            yErrMean = np.mean(yErrList)
            xStandardDeviation = np.std(xList)
            yStandardDeviation = np.std(yList)
            self.assertFloatsAlmostEqual(xErrMean, xStandardDeviation, rtol=0.2)   # rng dependent
            self.assertFloatsAlmostEqual(yErrMean, yStandardDeviation, rtol=0.2)   # rng dependent
            self.assertLess(abs(xMean - x), 3.0*xErrMean / nSamples**0.5)   # rng dependent
            self.assertLess(abs(yMean - y), 3.0*yErrMean / nSamples**0.5)   # rng dependent

    def testEdge(self):
        task = self.makeSingleFrameMeasurementTask("base_SdssCentroid")
        exposure, catalog = self.dataset.realize(10.0, task.schema, randomSeed=2)
        psfImage = exposure.getPsf().computeImage(self.center)
        # construct a box that won't fit the full PSF model
        bbox = psfImage.getBBox()
        bbox.grow(-5)
        subImage = lsst.afw.image.ExposureF(exposure, bbox)
        # we also need to install a smaller footprint, or NoiseReplacer
        # complains before we even get to measuring the centroid
        record = catalog[0]
        spanSet = lsst.afw.geom.SpanSet(bbox)
        newFootprint = lsst.afw.detection.Footprint(spanSet)
        peak = record.getFootprint().getPeaks()[0]
        newFootprint.addPeak(peak.getFx(), peak.getFy(), peak.getPeakValue())
        record.setFootprint(newFootprint)
        # just measure the one object we've prepared for
        task.measure(catalog, subImage)
        self.assertTrue(record.get("base_SdssCentroid_flag"))
        self.assertTrue(record.get("base_SdssCentroid_flag_edge"))

    def testNearEdge(self):
        task = self.makeSingleFrameMeasurementTask("base_SdssCentroid")
        exposure, catalog = self.dataset.realize(10.0, task.schema, randomSeed=2)
        psfImage = exposure.getPsf().computeImage(self.center)
        # construct a box that won't fit the full PSF model
        bbox = psfImage.getBBox()
        bbox.grow(-3)
        subImage = lsst.afw.image.ExposureF(exposure, bbox)
        # we also need to install a smaller footprint, or NoiseReplacer
        # complains before we even get to measuring the centroid
        record = catalog[0]
        spanSet = lsst.afw.geom.SpanSet(bbox)
        newFootprint = lsst.afw.detection.Footprint(spanSet)
        peak = record.getFootprint().getPeaks()[0]
        newFootprint.addPeak(peak.getFx(), peak.getFy(), peak.getPeakValue())
        record.setFootprint(newFootprint)
        # just measure the one object we've prepared for
        task.measure(catalog, subImage)
        self.assertTrue(record.get("base_SdssCentroid_flag_near_edge"))
        self.assertTrue(record.get("base_SdssCentroid_flag"))

    def testNo2ndDerivative(self):
        task = self.makeSingleFrameMeasurementTask("base_SdssCentroid")
        exposure, catalog = self.dataset.realize(10.0, task.schema, randomSeed=3)
        # cutout a subimage around object in the test image
        bbox = lsst.geom.Box2I(lsst.geom.Point2I(self.center), lsst.geom.Extent2I(1, 1))
        bbox.grow(20)
        subImage = lsst.afw.image.ExposureF(exposure, bbox)
        # A completely flat image will trigger the no 2nd derivative error
        subImage.image.array[:] = 0
        task.measure(catalog, subImage)
        self.assertTrue(catalog[0].get("base_SdssCentroid_flag"))
        self.assertTrue(catalog[0].get("base_SdssCentroid_flag_noSecondDerivative"))

    def testNotAtMaximum(self):
        task = self.makeSingleFrameMeasurementTask("base_SdssCentroid")
        exposure, catalog = self.dataset.realize(10.0, task.schema, randomSeed=4)
        # cutout a subimage around the object in the test image
        bbox = lsst.geom.Box2I(lsst.geom.Point2I(self.center), lsst.geom.Extent2I(1, 1))
        bbox.grow(20)
        subImage = lsst.afw.image.ExposureF(exposure, bbox)
        # zero out the central region, which will destroy the maximum
        subImage.image.array[18:22, 18:22] = 0
        task.measure(catalog, subImage)
        self.assertTrue(catalog[0].get("base_SdssCentroid_flag"))
        self.assertTrue(catalog[0].get("base_SdssCentroid_flag_notAtMaximum"))

    def testNegative(self):
        """Test that negative sources are well measured, without error flags.
        """
        dataset = lsst.meas.base.tests.TestDataset(self.bbox)
        dataset.addSource(-10000.0, self.center, negative=True)
        task = self.makeSingleFrameMeasurementTask("base_SdssCentroid")
        exposure, catalog = dataset.realize(10.0, task.schema, randomSeed=4)

        print("!!!!!!!")
        print("!!!!!!!")
        task.run(catalog, exposure)
        record = catalog[0]
        print(record)
        # import os; print(os.getpid()); import ipdb; ipdb.set_trace();
        self.assertFalse(record.get("base_SdssCentroid_flag"))
        self.assertFalse(record.get("base_SdssCentroid_flag_edge"))
        self.assertFalse(record.get("base_SdssCentroid_flag_notAtMaximum"))
        self.assertFalse(record.get("base_SdssCentroid_flag_near_edge"))
        self.assertFalse(record.get("base_SdssCentroid_flag_resetToPeak"))
        self.assertFalse(record.get("base_SdssCentroid_flag_badError"))
        self.assertFloatsAlmostEqual(record.get("base_SdssCentroid_x"), record.get("truth_x"), rtol=0.005)
        self.assertFloatsAlmostEqual(record.get("base_SdssCentroid_y"), record.get("truth_y"), rtol=0.005)


class SdssCentroidTransformTestCase(CentroidTransformTestCase,
                                    SingleFramePluginTransformSetupHelper,
                                    lsst.utils.tests.TestCase):
    controlClass = lsst.meas.base.SdssCentroidControl
    algorithmClass = lsst.meas.base.SdssCentroidAlgorithm
    transformClass = lsst.meas.base.SdssCentroidTransform
    flagNames = ('flag', 'flag_edge', 'flag_badData')
    singleFramePlugins = ('base_SdssCentroid',)
    forcedPlugins = ('base_SdssCentroid',)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
