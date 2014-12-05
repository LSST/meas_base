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

import numpy
import unittest

import lsst.afw.geom
import lsst.afw.image
import lsst.afw.detection
import lsst.afw.table
import lsst.meas.algorithms
import lsst.meas.base
import lsst.utils.tests

def makePluginAndCat(alg, name, control, metadata=False, centroid=None):
    schema = afwTable.SourceTable.makeMinimalSchema()
    if centroid:
        schema.addField(centroid + "_x", type=float)
        schema.addField(centroid + "_y", type=float)
        schema.addField(centroid + "_flag", type='Flag')
    if metadata:
        plugin = alg(control, name, schema, dafBase.PropertySet())
    else:
        plugin = alg(control, name, schema)
    cat = afwTable.SourceCatalog(schema)
    if centroid:
        cat.defineCentroid(centroid)
    return plugin, cat

class FluxTestCase(unittest.TestCase):

    def assertClose(self, a, b, rtol=1E-5, atol=1E-8):
        self.assert_(numpy.allclose(a, b, rtol=rtol, atol=atol), "\n%s\n!=\n%s" % (a, b))

    def assertNotClose(self, a, b, rtol=1E-5, atol=1E-8):
        self.assertFalse(numpy.allclose(a, b, rtol=rtol, atol=atol), "\n%s\n==\n%s" % (a, b))

    def setUp(self):
        self.exposure = lsst.afw.image.ExposureF(201, 201)
        crval = lsst.afw.coord.IcrsCoord(45.0*lsst.afw.geom.degrees, 45.0*lsst.afw.geom.degrees)
        crpix = lsst.afw.geom.Point2D(0,0)
        cdelt = (0.2 * lsst.afw.geom.arcseconds).asDegrees()
        wcs = lsst.afw.image.makeWcs(crval, crpix, cdelt, 0.0, 0.0, cdelt)
        self.exposure.setWcs(wcs)
        # for convenience, we'll put the source at the origin
        self.exposure.setXY0(lsst.afw.geom.Point2I(-100,-100))
        self.exposure.getMaskedImage().getVariance()[:] = 1.0
        self.psf = lsst.meas.algorithms.DoubleGaussianPsf(71, 71, 8.0, 15.0, 1.0)
        self.exposure.setPsf(self.psf)
        self.flux = 50.0
        psfImage = self.psf.computeImage()
        box = psfImage.getBBox()
        image = self.exposure.getMaskedImage().getImage()
        subImage = image.Factory(image, box, lsst.afw.image.PARENT, False)
        subImage.scaledPlus(self.flux, psfImage.convertF())
        self.footprint = lsst.afw.detection.Footprint(box)
        self.footprint.getPeaks().append(lsst.afw.detection.Peak(0,0))
        self.config = lsst.meas.base.SingleFrameMeasurementConfig()
        self.config.doReplaceWithNoise = False
    def tearDown(self):
        del self.psf
        del self.exposure
        del self.footprint

    def measure(self, radius=None):
        if radius is not None:
            self.config.algorithms["base_SincFlux"].radius = radius
         #   self.config.algorithms["correctfluxes"].apCorrRadius = radius
        schema = lsst.afw.table.SourceTable.makeMinimalSchema()
        task = lsst.meas.base.SingleFrameMeasurementTask(config=self.config, schema=schema)
        catalog = lsst.afw.table.SourceCatalog(schema)
        source = catalog.addNew()
        source.setFootprint(self.footprint)
        task.run(self.exposure, catalog)
        # flux.psf.psffactor should be 1.0 because it's just the dot product of the PSF with itself
        #self.assertClose(source.get("base_PsfFlux_psfFactor"), 1.0)
        # flux.gaussian.psffactor should be the dot product of a Gaussian with a double-Gaussian PSF.
        #self.assertNotClose(source.get("base_GaussianFlux_psfFactor"), 1.0, rtol=1E-2, atol=1E-2)
        return source

    def testGaussian(self):
        """Test that we can measure a Gaussian flux"""

        #self.config.algorithms["correctfluxes"].doApCorr = False
        source = self.measure()

        self.assertClose(self.flux, source.get("base_GaussianFlux_flux"), rtol=.05)
        self.assertTrue(numpy.isfinite(source.get("base_GaussianFlux_fluxSigma")))

def suite():
    """Returns a suite containing all the test cases in this module."""
    lsst.utils.tests.init()

    suites = []
    suites += unittest.makeSuite(FluxTestCase)
    suites += unittest.makeSuite(lsst.utils.tests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(exit=False):
    lsst.utils.tests.run(suite(), exit)

if __name__ == "__main__":
    run(True)
