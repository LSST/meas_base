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
import lsst.meas.base
import lsst.utils.tests
from lsst.meas.base.tests import (AlgorithmTestCase, FluxTransformTestCase,
                                  SingleFramePluginTransformSetupHelper)


class GaussianFluxTestCase(AlgorithmTestCase, lsst.utils.tests.TestCase):

    def setUp(self):
        self.bbox = lsst.geom.Box2I(lsst.geom.Point2I(-20, -30),
                                    lsst.geom.Extent2I(240, 1600))
        self.dataset = lsst.meas.base.tests.TestDataset(self.bbox)
        # first source is a point
        self.dataset.addSource(100000.0, lsst.geom.Point2D(50.1, 49.8))
        # second source is extended
        self.dataset.addSource(100000.0, lsst.geom.Point2D(149.9, 50.3),
                               lsst.afw.geom.Quadrupole(8, 9, 3))

    def tearDown(self):
        del self.bbox
        del self.dataset

    def makeAlgorithm(self, ctrl=None):
        """Construct an algorithm and return both it and its schema.
        """
        if ctrl is None:
            ctrl = lsst.meas.base.GaussianFluxControl()
        schema = lsst.meas.base.tests.TestDataset.makeMinimalSchema()
        algorithm = lsst.meas.base.GaussianFluxAlgorithm(ctrl, "base_GaussianFlux", schema)
        return algorithm, schema

    def testGaussians(self):
        """Test for correct instFlux given known position and shape.
        """
        task = self.makeSingleFrameMeasurementTask("base_GaussianFlux")
        # Results are RNG dependent; we choose a seed that is known to pass.
        exposure, catalog = self.dataset.realize(10.0, task.schema, randomSeed=0)
        task.run(catalog, exposure)
        for measRecord in catalog:
            self.assertFloatsAlmostEqual(measRecord.get("base_GaussianFlux_instFlux"),
                                         measRecord.get("truth_instFlux"), rtol=3E-3)

    def testMonteCarlo(self):
        """Test an ideal simulation, with no noise.

        Demonstrate that:

        - We get exactly the right answer, and
        - The reported uncertainty agrees with a Monte Carlo test of the noise.
        """
        algorithm, schema = self.makeAlgorithm()
        # Results are RNG dependent; we choose a seed that is known to pass.
        exposure, catalog = self.dataset.realize(1E-8, schema, randomSeed=1)
        record = catalog[0]
        instFlux = record.get("truth_instFlux")
        algorithm.measure(record, exposure)
        self.assertFloatsAlmostEqual(record.get("base_GaussianFlux_instFlux"), instFlux, rtol=1E-3)
        self.assertLess(record.get("base_GaussianFlux_instFluxErr"), 1E-3)
        for noise in (0.001, 0.01, 0.1):
            instFluxes = []
            instFluxErrs = []
            nSamples = 1000
            for repeat in range(nSamples):
                # By using ``repeat`` to seed the RNG, we get results which
                # fall within the tolerances defined below. If we allow this
                # test to be truly random, passing becomes RNG-dependent.
                exposure, catalog = self.dataset.realize(noise*instFlux, schema, randomSeed=repeat)
                record = catalog[1]
                algorithm.measure(record, exposure)
                instFluxes.append(record.get("base_GaussianFlux_instFlux"))
                instFluxErrs.append(record.get("base_GaussianFlux_instFluxErr"))
            instFluxMean = np.mean(instFluxes)
            instFluxErrMean = np.mean(instFluxErrs)
            instFluxStandardDeviation = np.std(instFluxes)
            self.assertFloatsAlmostEqual(instFluxErrMean, instFluxStandardDeviation, rtol=0.10)
            self.assertLess(abs(instFluxMean - instFlux), 2.0*instFluxErrMean / nSamples**0.5)

    def testForcedPlugin(self):
        task = self.makeForcedMeasurementTask("base_GaussianFlux")
        # Results of this test are RNG dependent: we choose seeds that are
        # known to pass.
        measWcs = self.dataset.makePerturbedWcs(self.dataset.exposure.getWcs(), randomSeed=2)
        measDataset = self.dataset.transform(measWcs)
        exposure, truthCatalog = measDataset.realize(10.0, measDataset.makeMinimalSchema(), randomSeed=2)
        refWcs = self.dataset.exposure.getWcs()
        refCat = self.dataset.catalog
        measCat = task.generateMeasCat(exposure, refCat, refWcs)
        task.attachTransformedFootprints(measCat, refCat, exposure, refWcs)
        task.run(measCat, exposure, refCat, refWcs)
        for measRecord, truthRecord in zip(measCat, truthCatalog):
            # Centroid tolerances set to ~ single precision epsilon
            self.assertFloatsAlmostEqual(measRecord.get("slot_Centroid_x"),
                                         truthRecord.get("truth_x"), rtol=1E-7)
            self.assertFloatsAlmostEqual(measRecord.get("slot_Centroid_y"),
                                         truthRecord.get("truth_y"), rtol=1E-7)
            self.assertFalse(measRecord.get("base_GaussianFlux_flag"))
            # GaussianFlux isn't designed to do a good job in forced mode,
            # because it doesn't account for changes in the PSF (and in fact
            # confuses them with changes in the WCS).  Hence, this is really
            # just a regression test, with the initial threshold set to just a
            # bit more than what it was found to be at one point.
            self.assertFloatsAlmostEqual(measRecord.get("base_GaussianFlux_instFlux"),
                                         truthCatalog.get("truth_instFlux"), rtol=0.3)
            self.assertLess(measRecord.get("base_GaussianFlux_instFluxErr"), 500.0)


class GaussianFluxTransformTestCase(FluxTransformTestCase, SingleFramePluginTransformSetupHelper,
                                    lsst.utils.tests.TestCase):
    controlClass = lsst.meas.base.GaussianFluxControl
    algorithmClass = lsst.meas.base.GaussianFluxAlgorithm
    transformClass = lsst.meas.base.GaussianFluxTransform
    singleFramePlugins = ('base_GaussianFlux',)
    forcedPlugins = ('base_GaussianFlux',)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
