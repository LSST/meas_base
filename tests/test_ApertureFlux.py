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
import lsst.afw.geom
import lsst.afw.image
import lsst.utils.tests
from lsst.meas.base import ApertureFluxAlgorithm
from lsst.meas.base.tests import (AlgorithmTestCase, FluxTransformTestCase,
                                  SingleFramePluginTransformSetupHelper)


class ApertureFluxTestCase(lsst.utils.tests.TestCase):
    """Test case for the ApertureFlux algorithm base class.
    """

    def setUp(self):
        self.bbox = lsst.geom.Box2I(lsst.geom.Point2I(20, -100), lsst.geom.Point2I(100, -20))
        self.exposure = lsst.afw.image.ExposureF(self.bbox)
        self.exposure.getMaskedImage().getImage().set(1.0)
        self.exposure.getMaskedImage().getVariance().set(0.25)
        self.ctrl = ApertureFluxAlgorithm.Control()

    def tearDown(self):
        del self.bbox
        del self.exposure

    def computeNaiveArea(self, position, radius):
        """Computes the area of a circular aperture.

        Calculates the area of the aperture by the "naive" approach of testing
        each pixel to see whether its center lies within the aperture.
        """
        x, y = np.meshgrid(np.arange(self.bbox.getBeginX(), self.bbox.getEndX()),
                           np.arange(self.bbox.getBeginY(), self.bbox.getEndY()))
        return ((x - position.getX())**2 + (y - position.getY())**2 <= radius**2).sum()

    def testNaive(self):
        positions = [lsst.geom.Point2D(60.0, -60.0),
                     lsst.geom.Point2D(60.5, -60.0),
                     lsst.geom.Point2D(60.0, -60.5),
                     lsst.geom.Point2D(60.5, -60.5)]
        radii = [12.0, 17.0]
        for position in positions:
            for radius in radii:
                ellipse = lsst.afw.geom.Ellipse(lsst.afw.geom.ellipses.Axes(radius, radius, 0.0), position)
                area = self.computeNaiveArea(position, radius)
                # test that this isn't the same as the sinc instFlux
                self.assertFloatsNotEqual(
                    ApertureFluxAlgorithm.computeSincFlux(self.exposure.getMaskedImage().getImage(),
                                                          ellipse, self.ctrl).instFlux, area)

                def check(method, image):
                    """Test that all instFlux measurement invocations work.

                    That is, that they return the expected value.
                    """
                    result = method(image, ellipse, self.ctrl)
                    self.assertFloatsAlmostEqual(result.instFlux, area)
                    self.assertFalse(result.getFlag(ApertureFluxAlgorithm.APERTURE_TRUNCATED.number))
                    self.assertFalse(result.getFlag(ApertureFluxAlgorithm.SINC_COEFFS_TRUNCATED.number))
                    if hasattr(image, "getVariance"):
                        self.assertFloatsAlmostEqual(result.instFluxErr, (area*0.25)**0.5)
                    else:
                        self.assertTrue(np.isnan(result.instFluxErr))
                check(ApertureFluxAlgorithm.computeNaiveFlux, self.exposure.getMaskedImage())
                check(ApertureFluxAlgorithm.computeNaiveFlux, self.exposure.getMaskedImage().getImage())
                check(ApertureFluxAlgorithm.computeFlux, self.exposure.getMaskedImage())
                check(ApertureFluxAlgorithm.computeFlux, self.exposure.getMaskedImage().getImage())
        # test failure conditions when the aperture itself is truncated
        invalid = ApertureFluxAlgorithm.computeNaiveFlux(
            self.exposure.getMaskedImage().getImage(),
            lsst.afw.geom.Ellipse(lsst.afw.geom.ellipses.Axes(12.0, 12.0),
                                  lsst.geom.Point2D(25.0, -60.0)),
            self.ctrl)
        self.assertTrue(invalid.getFlag(ApertureFluxAlgorithm.APERTURE_TRUNCATED.number))
        self.assertFalse(invalid.getFlag(ApertureFluxAlgorithm.SINC_COEFFS_TRUNCATED.number))
        self.assertTrue(np.isnan(invalid.instFlux))

    def testSinc(self):
        positions = [lsst.geom.Point2D(60.0, -60.0),
                     lsst.geom.Point2D(60.5, -60.0),
                     lsst.geom.Point2D(60.0, -60.5),
                     lsst.geom.Point2D(60.5, -60.5)]
        radii = [7.0, 9.0]
        for position in positions:
            for radius in radii:
                ellipse = lsst.afw.geom.Ellipse(lsst.afw.geom.ellipses.Axes(radius, radius, 0.0), position)
                area = ellipse.getCore().getArea()
                # test that this isn't the same as the naive instFlux
                self.assertFloatsNotEqual(
                    ApertureFluxAlgorithm.computeNaiveFlux(self.exposure.getMaskedImage().getImage(),
                                                           ellipse, self.ctrl).instFlux, area)

                def check(method, image):
                    # test that all the ways we could invoke sinc flux
                    # measurement produce the expected result
                    result = method(image, ellipse, self.ctrl)
                    self.assertFloatsAlmostEqual(result.instFlux, area, rtol=1E-3)
                    self.assertFalse(result.getFlag(ApertureFluxAlgorithm.APERTURE_TRUNCATED.number))
                    self.assertFalse(result.getFlag(ApertureFluxAlgorithm.SINC_COEFFS_TRUNCATED.number))
                    if hasattr(image, "getVariance"):
                        self.assertFalse(np.isnan(result.instFluxErr))
                    else:
                        self.assertTrue(np.isnan(result.instFluxErr))
                check(ApertureFluxAlgorithm.computeSincFlux, self.exposure.getMaskedImage())
                check(ApertureFluxAlgorithm.computeSincFlux, self.exposure.getMaskedImage().getImage())
                check(ApertureFluxAlgorithm.computeFlux, self.exposure.getMaskedImage())
                check(ApertureFluxAlgorithm.computeFlux, self.exposure.getMaskedImage().getImage())
        # test failure conditions when the aperture itself is truncated
        invalid1 = ApertureFluxAlgorithm.computeSincFlux(
            self.exposure.getMaskedImage().getImage(),
            lsst.afw.geom.Ellipse(lsst.afw.geom.ellipses.Axes(9.0, 9.0), lsst.geom.Point2D(25.0, -60.0)),
            self.ctrl)
        self.assertTrue(invalid1.getFlag(ApertureFluxAlgorithm.APERTURE_TRUNCATED.number))
        self.assertTrue(invalid1.getFlag(ApertureFluxAlgorithm.SINC_COEFFS_TRUNCATED.number))
        self.assertTrue(np.isnan(invalid1.instFlux))
        # test failure conditions when the aperture is not truncated, but the
        # sinc coeffs are
        invalid2 = ApertureFluxAlgorithm.computeSincFlux(
            self.exposure.getMaskedImage().getImage(),
            lsst.afw.geom.Ellipse(lsst.afw.geom.ellipses.Axes(9.0, 9.0), lsst.geom.Point2D(30.0, -60.0)),
            self.ctrl)
        self.assertFalse(invalid2.getFlag(ApertureFluxAlgorithm.APERTURE_TRUNCATED.number))
        self.assertTrue(invalid2.getFlag(ApertureFluxAlgorithm.SINC_COEFFS_TRUNCATED.number))
        self.assertFalse(np.isnan(invalid2.instFlux))


class CircularApertureFluxTestCase(AlgorithmTestCase, lsst.utils.tests.TestCase):
    """Test case for the CircularApertureFlux algorithm/plugin.
    """

    def setUp(self):
        self.bbox = lsst.geom.Box2I(lsst.geom.Point2I(0, 0),
                                    lsst.geom.Extent2I(100, 100))
        self.dataset = lsst.meas.base.tests.TestDataset(self.bbox)
        # first source is a point
        self.dataset.addSource(100000.0, lsst.geom.Point2D(49.5, 49.5))

    def tearDown(self):
        del self.bbox
        del self.dataset

    def testSingleFramePlugin(self):
        baseName = "base_CircularApertureFlux"
        config = self.makeSingleFrameMeasurementConfig(baseName)
        config.plugins[baseName].maxSincRadius = 20
        ctrl = config.plugins[baseName].makeControl()
        algMetadata = lsst.daf.base.PropertyList()
        task = self.makeSingleFrameMeasurementTask(config=config, algMetadata=algMetadata)
        exposure, catalog = self.dataset.realize(10.0, task.schema, randomSeed=0)
        task.run(catalog, exposure)
        radii = algMetadata.getArray("%s_RADII" % (baseName.upper(),))
        self.assertEqual(list(radii), list(ctrl.radii))
        for record in catalog:
            lastFlux = 0.0
            lastFluxErr = 0.0
            for n, radius in enumerate(radii):
                # Test that the flags are what we expect
                prefix = ApertureFluxAlgorithm.makeFieldPrefix(baseName, radius)
                if radius <= ctrl.maxSincRadius:
                    self.assertFalse(record.get(record.schema.join(prefix, "flag")))
                    self.assertFalse(record.get(record.schema.join(prefix, "flag_apertureTruncated")))
                    self.assertEqual(
                        record.get(record.schema.join(prefix, "flag_sincCoeffsTruncated")),
                        radius > 12
                    )
                else:
                    self.assertTrue(record.schema.join(prefix, "flag_sincCoeffsTruncated")
                                    not in record.getSchema())
                    self.assertEqual(record.get(record.schema.join(prefix, "flag")), radius >= 50)
                    self.assertEqual(record.get(record.schema.join(prefix, "flag_apertureTruncated")),
                                     radius >= 50)
                # Test that the instFluxes and uncertainties increase as we
                # increase the apertures, or that they match the true instFlux
                # within 3 sigma.  This is just a test as to whether the
                # values are reasonable.  As to whether the values are exactly
                # correct, we rely on the tests on ApertureFluxAlgorithm's
                # static methods, as the way the plugins code calls that is
                # extremely simple, so if the results we get are reasonable,
                # it's hard to imagine how they could be incorrect if
                # ApertureFluxAlgorithm's tests are valid.
                currentFlux = record.get(record.schema.join(prefix, "instFlux"))
                currentFluxErr = record.get(record.schema.join(prefix, "instFluxErr"))
                if not record.get(record.schema.join(prefix, "flag")):
                    self.assertTrue(currentFlux > lastFlux
                                    or (record.get("truth_instFlux") - currentFlux) < 3*currentFluxErr)
                    self.assertGreater(currentFluxErr, lastFluxErr)
                    lastFlux = currentFlux
                    lastFluxErr = currentFluxErr
                else:
                    self.assertTrue(np.isnan(currentFlux))
                    self.assertTrue(np.isnan(currentFluxErr))
            # When measuring an isolated point source with a sufficiently
            # large aperture, we should recover the known input instFlux.
            if record.get("truth_isStar") and record.get("parent") == 0:
                self.assertFloatsAlmostEqual(record.get("base_CircularApertureFlux_25_0_instFlux"),
                                             record.get("truth_instFlux"), rtol=0.02)

    def testForcedPlugin(self):
        baseName = "base_CircularApertureFlux"
        algMetadata = lsst.daf.base.PropertyList()
        task = self.makeForcedMeasurementTask(baseName, algMetadata=algMetadata)
        radii = algMetadata.getArray("%s_RADII" % (baseName.upper(),))
        measWcs = self.dataset.makePerturbedWcs(self.dataset.exposure.getWcs(), randomSeed=1)
        measDataset = self.dataset.transform(measWcs)
        exposure, truthCatalog = measDataset.realize(10.0, measDataset.makeMinimalSchema(), randomSeed=1)
        refCat = self.dataset.catalog
        refWcs = self.dataset.exposure.getWcs()
        measCat = task.generateMeasCat(exposure, refCat, refWcs)
        task.attachTransformedFootprints(measCat, refCat, exposure, refWcs)
        task.run(measCat, exposure, refCat, refWcs)
        for measRecord, truthRecord in zip(measCat, truthCatalog):
            # Centroid tolerances set to ~ single precision epsilon
            self.assertFloatsAlmostEqual(measRecord.get("slot_Centroid_x"),
                                         truthRecord.get("truth_x"), rtol=1E-7)
            self.assertFloatsAlmostEqual(measRecord.get("slot_Centroid_y"),
                                         truthRecord.get("truth_y"), rtol=1E-7)
            for n, radius in enumerate(radii):
                prefix = ApertureFluxAlgorithm.makeFieldPrefix(baseName, radius)
                self.assertFalse(measRecord.get(measRecord.schema.join(prefix, "flag")))
                # CircularApertureFlux isn't designed to do a good job in
                # forced mode, because it doesn't account for changes in the
                # PSF or changes in the WCS.  Hence, this is really just a
                # test to make sure the values are reasonable and that it runs
                # with no unexpected errors.
                self.assertFloatsAlmostEqual(measRecord.get(measRecord.schema.join(prefix, "instFlux")),
                                             truthCatalog.get("truth_instFlux"), rtol=1.0)
                self.assertLess(measRecord.get(measRecord.schema.join(prefix, "instFluxErr")), (n+1)*150.0)


class ApertureFluxTransformTestCase(FluxTransformTestCase, SingleFramePluginTransformSetupHelper,
                                    lsst.utils.tests.TestCase):

    class CircApFluxAlgorithmFactory:
        """Supply an empty ``PropertyList`` to `CircularApertureFluxAlgorithm`.

        This is a helper class to make testing more convenient.
        """

        def __call__(self, control, name, inputSchema):
            return lsst.meas.base.CircularApertureFluxAlgorithm(control, name, inputSchema,
                                                                lsst.daf.base.PropertyList())

    controlClass = lsst.meas.base.ApertureFluxAlgorithm.Control
    algorithmClass = CircApFluxAlgorithmFactory()
    transformClass = lsst.meas.base.ApertureFluxTransform
    flagNames = ('flag', 'flag_apertureTruncated', 'flag_sincCoeffsTruncated')
    singleFramePlugins = ('base_CircularApertureFlux',)
    forcedPlugins = ('base_CircularApertureFlux',)

    def testTransform(self):
        """Test `ApertureFluxTransform` with a synthetic catalog.
        """
        FluxTransformTestCase.testTransform(self, [ApertureFluxAlgorithm.makeFieldPrefix(self.name, r)
                                                   for r in self.control.radii])


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
