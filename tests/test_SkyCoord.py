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

import lsst.geom
import lsst.meas.base.tests
import lsst.utils.tests


class SkyCoordTestCase(lsst.meas.base.tests.AlgorithmTestCase, lsst.utils.tests.TestCase):

    def setUp(self):
        self.center = lsst.geom.Point2D(50.1, 49.8)
        self.bbox = lsst.geom.Box2I(lsst.geom.Point2I(-20, -30), lsst.geom.Extent2I(140, 160), invert=False)
        self.dataset = lsst.meas.base.tests.TestDataset(self.bbox)
        self.dataset.addSource(100000.0, self.center)

    def tearDown(self):
        del self.center
        del self.bbox
        del self.dataset

    def testSingleFramePlugin(self):
        task = self.makeSingleFrameMeasurementTask("base_SkyCoord")
        exposure, catalog = self.dataset.realize(10.0, task.schema, randomSeed=0)
        task.run(catalog, exposure)
        record = catalog[0]
        position = exposure.getWcs().skyToPixel(record.getCoord())
        self.assertFloatsAlmostEqual(position.getX(), record.get("truth_x"), rtol=1E-8)
        self.assertFloatsAlmostEqual(position.getY(), record.get("truth_y"), rtol=1E-8)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
