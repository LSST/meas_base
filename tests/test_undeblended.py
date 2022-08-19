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

"""Tests for measuring sources on undeblended images.
"""

import sys
import unittest
import warnings

import numpy as np

import lsst.geom
import lsst.afw.image as afwImage
import lsst.afw.table as afwTable
import lsst.afw.geom as afwGeom
import lsst.afw.detection as afwDetection
import lsst.afw.math as afwMath
import lsst.meas.base as measBase
import lsst.utils.tests


class UndeblendedTestCase(lsst.utils.tests.TestCase):
    def testUndeblendedMeasurement(self):
        """Check undeblended measurement and aperture correction.
        """
        width, height = 100, 100  # Dimensions of image
        x0, y0 = 1234, 5678  # Offset of image
        radius = 3.0  # Aperture radius

        # Position of first source; integer values, for convenience
        xCenter, yCenter = width//2, height//2
        xOffset, yOffset = 1, 1  # Offset from first source to second source
        instFlux1, instFlux2 = 1000, 1  # Flux of sources
        apCorrValue = 3.21  # Aperture correction value to apply

        image = afwImage.MaskedImageF(lsst.geom.ExtentI(width, height))
        image.setXY0(x0, y0)
        image.getVariance().set(1.0)

        schema = afwTable.SourceTable.makeMinimalSchema()
        schema.addField("centroid_x", type=np.float64)
        schema.addField("centroid_y", type=np.float64)
        schema.addField("centroid_flag", type='Flag')
        schema.getAliasMap().set("slot_Centroid", "centroid")

        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message="ignoreSlotPluginChecks", category=FutureWarning)
            sfmConfig = measBase.SingleFrameMeasurementConfig(ignoreSlotPluginChecks=True)
        algName = "base_CircularApertureFlux"

        for subConfig in (sfmConfig.plugins, sfmConfig.undeblended):
            subConfig.names = [algName]
            subConfig[algName].radii = [radius]
            # Disable sinc photometry because we're undersampled
            subConfig[algName].maxSincRadius = 0
        slots = sfmConfig.slots
        slots.centroid = "centroid"
        slots.shape = None
        slots.psfShape = None
        slots.apFlux = None
        slots.modelFlux = None
        slots.psfFlux = None
        slots.gaussianFlux = None
        slots.calibFlux = None

        fieldName = lsst.meas.base.CircularApertureFluxAlgorithm.makeFieldPrefix(algName, radius)
        measBase.addApCorrName(fieldName)

        apCorrConfig = measBase.ApplyApCorrConfig()
        apCorrConfig.proxies = {"undeblended_" + fieldName: fieldName}

        sfm = measBase.SingleFrameMeasurementTask(config=sfmConfig, schema=schema)
        apCorr = measBase.ApplyApCorrTask(config=apCorrConfig, schema=schema)

        cat = afwTable.SourceCatalog(schema)
        parent = cat.addNew()
        parent.set("centroid_x", x0 + xCenter)
        parent.set("centroid_y", y0 + yCenter)
        spanSetParent = afwGeom.SpanSet.fromShape(int(radius))
        spanSetParent = spanSetParent.shiftedBy(x0 + xCenter, y0 + yCenter)
        parent.setFootprint(afwDetection.Footprint(spanSetParent))

        # First child is bright, dominating the blend
        child1 = cat.addNew()
        child1.set("centroid_x", parent.get("centroid_x"))
        child1.set("centroid_y", parent.get("centroid_y"))
        child1.setParent(parent.getId())
        image[xCenter, yCenter, afwImage.LOCAL] = (instFlux1, 0, 0)
        spanSetChild1 = afwGeom.SpanSet.fromShape(1)
        spanSetChild1 = spanSetChild1.shiftedBy(x0 + xCenter, y0 + yCenter)
        foot1 = afwDetection.Footprint(spanSetChild1)
        child1.setFootprint(afwDetection.HeavyFootprintF(foot1, image))

        # Second child is fainter, but we want to be able to measure it!
        child2 = cat.addNew()
        child2.set("centroid_x", parent.get("centroid_x") + xOffset)
        child2.set("centroid_y", parent.get("centroid_y") + yOffset)
        child2.setParent(parent.getId())
        image[xCenter + xOffset, yCenter + yOffset, afwImage.LOCAL] = (instFlux2, 0, 0)
        spanSetChild2 = afwGeom.SpanSet.fromShape(1)
        tmpPoint = (x0 + xCenter + xOffset, y0 + yCenter + yOffset)
        spanSetChild2 = spanSetChild2.shiftedBy(*tmpPoint)
        foot2 = afwDetection.Footprint(spanSetChild2)
        child2.setFootprint(afwDetection.HeavyFootprintF(foot2, image))

        spans = foot1.spans.union(foot2.spans)
        bbox = lsst.geom.Box2I()
        bbox.include(foot1.getBBox())
        bbox.include(foot2.getBBox())
        parent.setFootprint(afwDetection.Footprint(spans, bbox))

        exposure = afwImage.makeExposure(image)

        sfm.run(cat, exposure)

        def checkSource(source, baseName, expectedFlux):
            """Check that we get the expected results.
            """
            self.assertEqual(source.get(baseName + "_instFlux"), expectedFlux)
            self.assertGreater(source.get(baseName + "_instFluxErr"), 0)
            self.assertFalse(source.get(baseName + "_flag"))

        # Deblended
        checkSource(child1, fieldName, instFlux1)
        checkSource(child2, fieldName, instFlux2)

        # Undeblended
        checkSource(child1, "undeblended_" + fieldName, instFlux1 + instFlux2)
        checkSource(child2, "undeblended_" + fieldName, instFlux1 + instFlux2)

        # Apply aperture correction
        apCorrMap = afwImage.ApCorrMap()
        apCorrMap[fieldName + "_instFlux"] = afwMath.ChebyshevBoundedField(
            image.getBBox(),
            apCorrValue*np.ones((1, 1), dtype=np.float64)
        )
        apCorrMap[fieldName + "_instFluxErr"] = afwMath.ChebyshevBoundedField(
            image.getBBox(),
            apCorrValue*np.zeros((1, 1), dtype=np.float64)
        )

        apCorr.run(cat, apCorrMap)

        # Deblended
        checkSource(child1, fieldName, instFlux1*apCorrValue)
        checkSource(child2, fieldName, instFlux2*apCorrValue)

        # Undeblended
        checkSource(child1, "undeblended_" + fieldName, (instFlux1 + instFlux2)*apCorrValue)
        checkSource(child2, "undeblended_" + fieldName, (instFlux1 + instFlux2)*apCorrValue)

        self.assertIn(fieldName + "_apCorr", schema)
        self.assertIn(fieldName + "_apCorrErr", schema)
        self.assertIn("undeblended_" + fieldName + "_apCorr", schema)
        self.assertIn("undeblended_" + fieldName + "_apCorrErr", schema)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    setup_module(sys.modules[__name__])
    unittest.main()
