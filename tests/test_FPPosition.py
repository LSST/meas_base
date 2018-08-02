#
# LSST Data Management System
# Copyright 2008-2015 AURA/LSST.
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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#

import unittest

import lsst.geom
import lsst.meas.base.tests
import lsst.utils.tests

from lsst.afw.cameraGeom.testUtils import DetectorWrapper
from lsst.meas.base.tests import AlgorithmTestCase


class FPPositionTestCase(AlgorithmTestCase, lsst.utils.tests.TestCase):

    def setUp(self):
        # Define point2D object which are distributed about a detector
        self.positions = [lsst.geom.Point2D(*x) for x in ((50.1, 49.8), (12, 15.6), (13.4, 100.0))]
        # Define a box which will be used to as boundaries to construct an detector object
        self.bbox = lsst.geom.Box2I(lsst.geom.Point2I(-20, -30), lsst.geom.Extent2I(140, 160), invert=False)
        self.dataset = lsst.meas.base.tests.TestDataset(self.bbox)
        # Add in sources to synthetic dataset at defined positions with an arbitrary value
        for pos in self.positions:
            self.dataset.addSource(100000.0, pos)
        self.dw = DetectorWrapper()

    def tearDown(self):
        del self.positions
        del self.bbox
        del self.dataset
        del self.dw

    def testFPPosition(self):
        task = self.makeSingleFrameMeasurementTask("base_FPPosition")
        exposure, catalog = self.dataset.realize(10.0, task.schema, randomSeed=0)
        exposure.setDetector(self.dw.detector)
        task.run(catalog, exposure)
        pointKey = lsst.afw.table.Point2DKey(catalog.schema["base_FPPosition"])
        # Compare the derived focal plane position to the true position.
        # True position is calculated as pixel coordinate plus half pixel offset
        # to the center of a pixel, times the pixel scale.
        for record, pos in zip(catalog, self.positions):
            self.assertFalse(record.get("base_FPPosition_flag"))
            point = record.get(pointKey)
            self.assertAlmostEqual(point.getX(), self.dw.pixelSize[0] * (0.5 + pos[0]), 3)
            self.assertAlmostEqual(point.getY(), self.dw.pixelSize[1] * (0.5 + pos[1]), 3)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
