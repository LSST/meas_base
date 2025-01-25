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
import lsst.utils.tests
import lsst.meas.base.tests


class PixelFlagsTestCase(lsst.meas.base.tests.AlgorithmTestCase, lsst.utils.tests.TestCase):

    def setUp(self):
        self.center = lsst.geom.Point2D(50.1, 49.8)
        self.bbox = lsst.geom.Box2I(lsst.geom.Point2I(-20, -30),
                                    lsst.geom.Extent2I(140, 160))
        self.dataset = lsst.meas.base.tests.TestDataset(self.bbox)

    def testNoFlags(self):
        self.dataset.addSource(100000.0, self.center)
        task = self.makeSingleFrameMeasurementTask("base_PixelFlags")
        exposure, catalog = self.dataset.realize(10.0, task.schema, randomSeed=0)
        task.run(catalog, exposure)
        record = catalog[0]
        self.assertFalse(record.get("base_PixelFlags_flag"))
        self.assertFalse(record.get("base_PixelFlags_flag_edge"))
        self.assertFalse(record.get("base_PixelFlags_flag_edgeCenter"))
        self.assertFalse(record.get("base_PixelFlags_flag_edgeCenterAll"))
        self.assertFalse(record.get("base_PixelFlags_flag_nodata"))
        self.assertFalse(record.get("base_PixelFlags_flag_nodataCenter"))
        self.assertFalse(record.get("base_PixelFlags_flag_nodataCenterAll"))
        self.assertFalse(record.get("base_PixelFlags_flag_interpolated"))
        self.assertFalse(record.get("base_PixelFlags_flag_interpolatedCenter"))
        self.assertFalse(record.get("base_PixelFlags_flag_interpolatedCenterAll"))
        self.assertFalse(record.get("base_PixelFlags_flag_saturated"))
        self.assertFalse(record.get("base_PixelFlags_flag_saturatedCenter"))
        self.assertFalse(record.get("base_PixelFlags_flag_saturatedCenterAll"))
        self.assertFalse(record.get("base_PixelFlags_flag_cr"))
        self.assertFalse(record.get("base_PixelFlags_flag_crCenter"))
        self.assertFalse(record.get("base_PixelFlags_flag_crCenterAll"))
        self.assertFalse(record.get("base_PixelFlags_flag_bad"))
        self.assertFalse(record.get("base_PixelFlags_flag_badCenter"))
        self.assertFalse(record.get("base_PixelFlags_flag_badCenterAll"))

    def testSomeFlags(self):
        self.dataset.addSource(100000.0, self.center)
        task = self.makeSingleFrameMeasurementTask("base_PixelFlags")
        exposure, catalog = self.dataset.realize(10.0, task.schema, randomSeed=0)
        # one cr pixel outside the center
        cosmicray = exposure.mask.getPlaneBitMask("CR")
        x = round(self.center.x)
        y = round(self.center.y)
        exposure.mask[x+3, y+4] |= cosmicray
        # one interpolated pixel near the center
        interpolated = exposure.mask.getPlaneBitMask("INTRP")
        exposure.mask[self.center] |= interpolated
        # all pixels in the center are bad
        bad = exposure.mask.getPlaneBitMask("BAD")
        exposure.mask[x-1:x+2, y-1:y+2] |= bad
        task.run(catalog, exposure)
        record = catalog[0]

        self.assertFalse(record.get("base_PixelFlags_flag_edge"))
        self.assertFalse(record.get("base_PixelFlags_flag_edgeCenter"))
        self.assertFalse(record.get("base_PixelFlags_flag_edgeCenterAll"))
        self.assertFalse(record.get("base_PixelFlags_flag_nodata"))
        self.assertFalse(record.get("base_PixelFlags_flag_nodataCenter"))
        self.assertFalse(record.get("base_PixelFlags_flag_nodataCenterAll"))
        self.assertTrue(record.get("base_PixelFlags_flag_cr"))
        self.assertFalse(record.get("base_PixelFlags_flag_crCenter"))
        self.assertFalse(record.get("base_PixelFlags_flag_crCenterAll"))
        self.assertTrue(record.get("base_PixelFlags_flag_interpolated"))
        self.assertTrue(record.get("base_PixelFlags_flag_interpolatedCenter"))
        self.assertFalse(record.get("base_PixelFlags_flag_interpolatedCenterAll"))
        self.assertTrue(record.get("base_PixelFlags_flag_bad"))
        self.assertTrue(record.get("base_PixelFlags_flag_badCenter"))
        self.assertTrue(record.get("base_PixelFlags_flag_badCenterAll"))

    def testOffimageNonfinite(self):
        """Test that flag_offimage and flag_edge are set correctly for
        nonfinite sources.
        """
        # These three will be explicitly set to have a non-finite centroid; we
        # cannot add sources off the image in TestDataset.
        self.dataset.addSource(100000, lsst.geom.Point2D(20, 20))
        self.dataset.addSource(100000, lsst.geom.Point2D(20, 100))
        self.dataset.addSource(100000, lsst.geom.Point2D(100, 30))
        task = self.makeSingleFrameMeasurementTask("base_PixelFlags")
        exposure, catalog = self.dataset.realize(10.0, task.schema, randomSeed=0)
        catalog[0]["slot_Centroid_x"] = np.inf
        catalog[1]["slot_Centroid_x"] = -np.inf
        catalog[2]["slot_Centroid_y"] = np.nan
        task.run(catalog, exposure)
        np.testing.assert_array_equal(catalog["base_PixelFlags_flag_offimage"],
                                      np.array([True, True, True]))
        np.testing.assert_array_equal(catalog["base_PixelFlags_flag_edge"],
                                      np.array([True, True, True]))
        np.testing.assert_array_equal(catalog["base_PixelFlags_flag_edgeCenter"],
                                      np.array([True, True, True]))
        np.testing.assert_array_equal(catalog["base_PixelFlags_flag_edgeCenterAll"],
                                      np.array([True, True, True]))

    def testOffimage(self):
        """Test that sources at the boundary of the image get flag_offimage
        and flag_edge[Center[All]] set appropriately.
        """
        # These four will be explicitly set to values at the boundary; we
        # cannot add sources off the image in TestDataset.
        self.dataset.addSource(100000, lsst.geom.Point2D(20, 20))
        self.dataset.addSource(100000, lsst.geom.Point2D(20, 100))
        self.dataset.addSource(100000, lsst.geom.Point2D(40, 100))
        self.dataset.addSource(100000, lsst.geom.Point2D(80, 30))

        task = self.makeSingleFrameMeasurementTask("base_PixelFlags")
        exposure, catalog = self.dataset.realize(10.0, task.schema, randomSeed=0)

        # Mask 4 pixels at the y-edges (note numpy vs. afw flipped indices!)
        edgeBit = exposure.mask.getPlaneBitMask("EDGE")
        exposure.mask.array[:4, :] |= edgeBit
        exposure.mask.array[-4:, ] |= edgeBit
        # Mask 4 pixels at the x-edges
        nodataBit = exposure.mask.getPlaneBitMask("NO_DATA")
        exposure.mask.array[:, :4] |= nodataBit
        exposure.mask.array[:, -4:] |= nodataBit

        # All of these should get Center and CenterAll.
        catalog[0]["slot_Centroid_x"] = -20.5  # on image
        catalog[1]["slot_Centroid_x"] = -20.51  # off image
        catalog[2]["slot_Centroid_y"] = 129.4  # on image
        catalog[3]["slot_Centroid_y"] = 129.50  # off image

        task.run(catalog, exposure)

        np.testing.assert_array_equal(catalog["base_PixelFlags_flag_offimage"],
                                      np.array([False, True, False, True]))

        # Centroid and footprints don't line up, so edge is only set for the
        # sources that have off-image centroids.
        np.testing.assert_array_equal(catalog["base_PixelFlags_flag_edge"],
                                      np.array([False, True, False, True]))
        # Center[All] is set for off-image and on-EDGE sources, because it's
        # based on the centroid location.
        np.testing.assert_array_equal(catalog["base_PixelFlags_flag_edgeCenter"],
                                      np.array([False, True, True, True]))
        np.testing.assert_array_equal(catalog["base_PixelFlags_flag_edgeCenterAll"],
                                      np.array([False, True, True, True]))

        # Centroid and footprints don't line up, so nodata is not set.
        np.testing.assert_array_equal(catalog["base_PixelFlags_flag_nodata"],
                                      np.array([False, False, False, False]))
        # Center[All] is set for on-NO_DATA sources, because it's based on the
        # centroid location.
        np.testing.assert_array_equal(catalog["base_PixelFlags_flag_nodataCenter"],
                                      np.array([True, True, False, False]))
        np.testing.assert_array_equal(catalog["base_PixelFlags_flag_nodataCenterAll"],
                                      np.array([True, True, False, False]))

    def testEdgeFlags(self):
        """Test that edge/edgeCenter/edgeCenterAll and the equivalent nodata
        flags are set appropriately.
        """
        # Half of the Central 3x3 pixels are on the masked NO_DATA; no offimage or
        # nodataCenterAll, but nodata and nodataCenter.
        self.dataset.addSource(100000, lsst.geom.Point2D(115, 30))
        # Central 3x3 pixels all masked EDGE; no offimage but all edge flags.
        self.dataset.addSource(100000, lsst.geom.Point2D(10, 127))

        task = self.makeSingleFrameMeasurementTask("base_PixelFlags")
        exposure, catalog = self.dataset.realize(10.0, task.schema, randomSeed=0)

        # Mask 4 pixels at the y-edges (note numpy vs. afw flipped indices!)
        edgeBit = exposure.mask.getPlaneBitMask("EDGE")
        exposure.mask.array[:4, :] |= edgeBit
        exposure.mask.array[-4:, ] |= edgeBit
        # Mask 4 pixels at the x-edges
        nodataBit = exposure.mask.getPlaneBitMask("NO_DATA")
        exposure.mask.array[:, :4] |= nodataBit
        exposure.mask.array[:, -4:] |= nodataBit

        task.run(catalog, exposure)

        np.testing.assert_array_equal(catalog["base_PixelFlags_flag_offimage"],
                                      np.array([False, False]))

        np.testing.assert_array_equal(catalog["base_PixelFlags_flag_nodata"],
                                      np.array([True, False]))
        np.testing.assert_array_equal(catalog["base_PixelFlags_flag_nodataCenter"],
                                      np.array([True, False]))
        np.testing.assert_array_equal(catalog["base_PixelFlags_flag_nodataCenterAll"],
                                      np.array([False, False]))

        np.testing.assert_array_equal(catalog["base_PixelFlags_flag_edge"],
                                      np.array([False, True]))
        np.testing.assert_array_equal(catalog["base_PixelFlags_flag_edgeCenter"],
                                      np.array([False, True]))
        np.testing.assert_array_equal(catalog["base_PixelFlags_flag_edgeCenterAll"],
                                      np.array([False, True]))


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
