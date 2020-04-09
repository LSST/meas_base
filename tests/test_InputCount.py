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

"""Tests for InputCounts measurement algorithm.
"""
import numpy as np
import itertools
from collections import namedtuple

import unittest
import lsst.utils.tests

import lsst.geom
import lsst.afw.detection as afwDetection
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.table as afwTable
import lsst.meas.base as measBase


try:
    display
    import matplotlib.pyplot as plt
    import matplotlib.patches as patches
except NameError:
    display = False


def ccdVennDiagram(exp, showImage=True, legendLocation='best'):
    """Display and exposure & bounding boxes for images which go into a coadd.

    Parameters
    ----------
    exp : `lsst.afw.image.Exposure`
        The exposure object to plot, must be the product of a coadd
    showImage : `bool`, optional
        Plot image data in addition to its bounding box
    legendLocation : `str`, optional
        Matplotlib legend location code. Can be: ``'best'``, ``'upper
        right'``, ``'upper left'``, ``'lower left'``, ``'lower right'``,
        ``'right'``, ``'center left'``, ``'center right'``, ``'lower
        center'``, ``'upper center'``, ``'center'``

    """
    # Create the figure object
    fig = plt.figure()
    # Use all the built in matplotib line style attributes to create a list of
    # the possible styles
    linestyles = ['solid', 'dashed', 'dashdot', 'dotted']
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
    # Calculate the cartisian product of the styles, and randomize the order,
    # to help each CCD get it's own color
    pcomb = np.random.permutation(list(itertools.product(colors, linestyles)))
    # Filter out a black solid box, as that will be the style of the given exp
    # object
    pcomb = pcomb[((pcomb[:, 0] == 'k') * (pcomb[:, 1] == 'solid')) is False]
    # Get the image properties
    origin = lsst.geom.PointD(exp.getXY0())
    mainBox = exp.getBBox().getCorners()
    # Plot the exposure
    plt.gca().add_patch(patches.Rectangle((0, 0), *list(mainBox[2]-mainBox[0]), fill=False, label="exposure"))
    # Grab all of the CCDs that went into creating the exposure
    ccds = exp.getInfo().getCoaddInputs().ccds
    # Loop over and plot the extents of each ccd
    for i, ccd in enumerate(ccds):
        ccdBox = lsst.geom.Box2D(ccd.getBBox())
        ccdCorners = ccdBox.getCorners()
        coaddCorners = [exp.getWcs().skyToPixel(ccd.getWcs().pixelToSky(point))
                        + (lsst.geom.PointD() - origin) for point in ccdCorners]
        plt.gca().add_patch(patches.Rectangle(coaddCorners[0], *list(coaddCorners[2]-coaddCorners[0]),
                                              fill=False, color=pcomb[i][0], ls=pcomb[i][1],
                                              label="CCD{}".format(i)))
    # If showImage is true, plot the data contained in exp as well as the
    # boundaries
    if showImage:
        plt.imshow(exp.getMaskedImage().getArrays()[0], cmap='Greys', origin='lower')
        plt.colorbar()
    # Adjust plot parameters and plot
    plt.gca().relim()
    plt.gca().autoscale_view()
    ylim = plt.gca().get_ylim()
    xlim = plt.gca().get_xlim()
    plt.gca().set_ylim(1.5*ylim[0], 1.5*ylim[1])
    plt.gca().set_xlim(1.5*xlim[0], 1.5*xlim[1])
    plt.legend(loc=legendLocation)
    fig.canvas.draw()
    plt.show()


class InputCountTest(lsst.utils.tests.TestCase):

    def testInputCounts(self, showPlot=False):
        # Generate a simulated coadd of four overlapping-but-offset CCDs.
        # Populate it with three sources.
        # Demonstrate that we can correctly recover the number of images which
        # contribute to each source.

        size = 20  # Size of images (pixels)
        value = 100.0  # Source flux

        ccdPositions = [
            lsst.geom.Point2D(8, 0),
            lsst.geom.Point2D(10, 10),
            lsst.geom.Point2D(-8, -8),
            lsst.geom.Point2D(-8, 8)
        ]

        # Represent sources by a tuple of position and expected number of
        # contributing CCDs (based on the size/positions given above).
        Source = namedtuple("Source", ["pos", "count"])
        sources = [
            Source(pos=lsst.geom.Point2D(6, 6), count=2),
            Source(pos=lsst.geom.Point2D(10, 10), count=3),
            Source(pos=lsst.geom.Point2D(14, 14), count=1)
        ]

        # These lines are used in the creation of WCS information
        scale = 1.0e-5 * lsst.geom.degrees
        cdMatrix = afwGeom.makeCdMatrix(scale=scale)
        crval = lsst.geom.SpherePoint(0.0, 0.0, lsst.geom.degrees)

        # Construct the info needed to set the exposure object
        imageBox = lsst.geom.Box2I(lsst.geom.Point2I(0, 0), lsst.geom.Extent2I(size, size))
        wcsRef = afwGeom.makeSkyWcs(crpix=lsst.geom.Point2D(0, 0), crval=crval, cdMatrix=cdMatrix)

        # Create the exposure object, and set it up to be the output of a coadd
        exp = afwImage.ExposureF(size, size)
        exp.setWcs(wcsRef)
        exp.getInfo().setCoaddInputs(afwImage.CoaddInputs(afwTable.ExposureTable.makeMinimalSchema(),
                                                          afwTable.ExposureTable.makeMinimalSchema()))

        # Set the fake CCDs that "went into" making this coadd, using the
        # differing wcs objects created above.
        ccds = exp.getInfo().getCoaddInputs().ccds
        for pos in ccdPositions:
            record = ccds.addNew()
            record.setWcs(afwGeom.makeSkyWcs(crpix=pos, crval=crval, cdMatrix=cdMatrix))
            record.setBBox(imageBox)
            record.setValidPolygon(afwGeom.Polygon(lsst.geom.Box2D(imageBox)))

        # Configure a SingleFrameMeasurementTask to run InputCounts.
        measureSourcesConfig = measBase.SingleFrameMeasurementConfig()
        measureSourcesConfig.plugins.names = ["base_PeakCentroid", "base_InputCount"]
        measureSourcesConfig.slots.centroid = "base_PeakCentroid"
        measureSourcesConfig.slots.psfFlux = None
        measureSourcesConfig.slots.apFlux = None
        measureSourcesConfig.slots.modelFlux = None
        measureSourcesConfig.slots.gaussianFlux = None
        measureSourcesConfig.slots.calibFlux = None
        measureSourcesConfig.slots.shape = None
        measureSourcesConfig.validate()
        schema = afwTable.SourceTable.makeMinimalSchema()
        task = measBase.SingleFrameMeasurementTask(schema, config=measureSourcesConfig)
        catalog = afwTable.SourceCatalog(schema)

        # Add simulated sources to the measurement catalog.
        for src in sources:
            spans = afwGeom.SpanSet.fromShape(1)
            spans = spans.shiftedBy(int(src.pos.getX()), int(src.pos.getY()))
            foot = afwDetection.Footprint(spans)
            peak = foot.getPeaks().addNew()
            peak.setFx(src.pos[0])
            peak.setFy(src.pos[1])
            peak.setPeakValue(value)
            catalog.addNew().setFootprint(foot)

        task.run(catalog, exp)

        for src, rec in zip(sources, catalog):
            self.assertEqual(rec.get("base_InputCount_value"), src.count)

        if display:
            ccdVennDiagram(exp)

    def _preparePlugin(self, addCoaddInputs):
        """Prepare a `SingleFrameInputCountPlugin` for running.

        Sets up the plugin to run on an empty catalog together with a
        synthetic, content-free `~lsst.afw.image.ExposureF`.

        Parameters
        ----------
        addCoaddInputs : `bool`
            Should we add the coadd inputs?

        Returns
        -------
        inputCount : `SingleFrameInputCountPlugin`
            Initialized measurement plugin.
        catalog : `lsst.afw.table.SourceCatalog`
            Empty Catalog.
        exp : `lsst.afw.image.ExposureF`
            Synthetic exposure.
        """
        exp = afwImage.ExposureF(20, 20)
        scale = 1.0e-5*lsst.geom.degrees
        wcs = afwGeom.makeSkyWcs(crpix=lsst.geom.Point2D(0, 0),
                                 crval=lsst.geom.SpherePoint(0.0, 0.0, lsst.geom.degrees),
                                 cdMatrix=afwGeom.makeCdMatrix(scale=scale))
        exp.setWcs(wcs)
        if addCoaddInputs:
            exp.getInfo().setCoaddInputs(afwImage.CoaddInputs(afwTable.ExposureTable.makeMinimalSchema(),
                                                              afwTable.ExposureTable.makeMinimalSchema()))
            ccds = exp.getInfo().getCoaddInputs().ccds
            record = ccds.addNew()
            record.setWcs(wcs)
            record.setBBox(exp.getBBox())
            record.setValidPolygon(afwGeom.Polygon(lsst.geom.Box2D(exp.getBBox())))

        schema = afwTable.SourceTable.makeMinimalSchema()
        measBase.SingleFramePeakCentroidPlugin(measBase.SingleFramePeakCentroidConfig(),
                                               "centroid", schema, None)
        schema.getAliasMap().set("slot_Centroid", "centroid")
        inputCount = measBase.SingleFrameInputCountPlugin(measBase.InputCountConfig(),
                                                          "inputCount", schema, None)
        catalog = afwTable.SourceCatalog(schema)
        return inputCount, catalog, exp

    def testBadCentroid(self):
        """Test that a NaN centroid raises and results in a correct flag.

        The flag from the centroid slot should propagate to the badCentroid
        flag on InputCount and the algorithm should throw a MeasurementError
        when it encounters a NaN position.
        """
        inputCount, catalog, exp = self._preparePlugin(True)
        record = catalog.addNew()

        # The inputCount's badCentroid flag is an alias to the centroid's
        # global flag, so it should be set immediately.
        record.set("centroid_flag", True)
        self.assertTrue(record.get("inputCount_flag_badCentroid"))

        # Even though the source is flagged as bad, if the position is good we
        # should still get a measurement.
        record.set("slot_Centroid_x", 10)
        record.set("slot_Centroid_y", 10)
        inputCount.measure(record, exp)
        self.assertTrue(record.get("inputCount_flag_badCentroid"))
        self.assertEqual(record.get("inputCount_value"), 1)

        # The centroid is bad (with a NaN) even though the centroid isn't
        # flagged, so we should get a MeasurementError indicating an expected
        # failure.
        record = catalog.addNew()
        record.set("slot_Centroid_x", float("nan"))
        record.set("slot_Centroid_y", 12.345)
        record.set("centroid_flag", False)
        with self.assertRaises(measBase.MeasurementError) as measErr:
            inputCount.measure(record, exp)

        # Calling the fail() method should set the global flag.
        record = catalog.addNew()
        record.set("inputCount_flag", False)
        inputCount.fail(record, measErr.exception)
        self.assertTrue(record.get("inputCount_flag"))

    def testBadCoaddInputs(self):
        """Test that no coadd inputs raises and results in a correct flag.

        When there are no coadd inputs on the input exposure we should throw a
        MeasurementError and set both the global flag and flag_noInputs.
        """
        inputCount, catalog, exp = self._preparePlugin(False)
        record = catalog.addNew()

        # Initially, the record is not flagged.
        self.assertFalse(record.get("inputCount_flag"))
        self.assertFalse(record.get("inputCount_flag_noInputs"))

        # There are no coadd inputs, so we should get a MeasurementError
        # indicating an expected failure.
        with self.assertRaises(measBase.MeasurementError) as measErr:
            inputCount.measure(record, exp)

        # Calling the fail() method should set the noInputs and global flags.
        inputCount.fail(record, measErr.exception)
        self.assertTrue(record.get("inputCount_flag"))
        self.assertTrue(record.get("inputCount_flag_noInputs"))


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
