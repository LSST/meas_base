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

import warnings

import numpy as np

import lsst.geom
import lsst.afw.table
import lsst.afw.image
import lsst.afw.detection
import lsst.afw.cameraGeom.testUtils
import lsst.afw.geom
import lsst.afw.coord
import lsst.daf.base
import lsst.pex.exceptions

from .sfm import SingleFrameMeasurementTask
from .forcedMeasurement import ForcedMeasurementTask
from ._measBaseLib import CentroidResultKey

__all__ = ("BlendContext", "TestDataset", "AlgorithmTestCase", "TransformTestCase",
           "SingleFramePluginTransformSetupHelper", "ForcedPluginTransformSetupHelper",
           "FluxTransformTestCase", "CentroidTransformTestCase")


class BlendContext:
    """Context manager which adds multiple overlapping sources and a parent.

    Notes
    -----
    This is used as the return value for `TestDataset.addBlend`, and this is
    the only way it should be used.
    """

    def __init__(self, owner):
        self.owner = owner
        self.parentRecord = self.owner.catalog.addNew()
        self.parentImage = lsst.afw.image.ImageF(self.owner.exposure.getBBox())
        self.children = []

    def __enter__(self):
        # BlendContext is its own context manager, so we just return self.
        return self

    def addChild(self, instFlux, centroid, shape=None):
        """Add a child to the blend; return corresponding truth catalog record.

        instFlux : `float`
            Total instFlux of the source to be added.
        centroid : `lsst.geom.Point2D`
            Position of the source to be added.
        shape : `lsst.afw.geom.Quadrupole`
            Second moments of the source before PSF convolution.  Note that
            the truth catalog records post-convolution moments)
        """
        record, image = self.owner.addSource(instFlux, centroid, shape)
        record.set(self.owner.keys["parent"], self.parentRecord.getId())
        self.parentImage += image
        self.children.append((record, image))
        return record

    def __exit__(self, type_, value, tb):
        # We're not using the context manager for any kind of exception safety
        # or guarantees; we just want the nice "with" statement syntax.

        if type_ is not None:
            # exception was raised; just skip all this and let it propagate
            return

        # On exit, compute and set the truth values for the parent object.
        self.parentRecord.set(self.owner.keys["nChild"], len(self.children))
        # Compute instFlux from sum of component fluxes
        instFlux = 0.0
        for record, image in self.children:
            instFlux += record.get(self.owner.keys["instFlux"])
        self.parentRecord.set(self.owner.keys["instFlux"], instFlux)
        # Compute centroid from instFlux-weighted mean of component centroids
        x = 0.0
        y = 0.0
        for record, image in self.children:
            w = record.get(self.owner.keys["instFlux"])/instFlux
            x += record.get(self.owner.keys["centroid"].getX())*w
            y += record.get(self.owner.keys["centroid"].getY())*w
        self.parentRecord.set(self.owner.keys["centroid"], lsst.geom.Point2D(x, y))
        # Compute shape from instFlux-weighted mean of offset component shapes
        xx = 0.0
        yy = 0.0
        xy = 0.0
        for record, image in self.children:
            w = record.get(self.owner.keys["instFlux"])/instFlux
            dx = record.get(self.owner.keys["centroid"].getX()) - x
            dy = record.get(self.owner.keys["centroid"].getY()) - y
            xx += (record.get(self.owner.keys["shape"].getIxx()) + dx**2)*w
            yy += (record.get(self.owner.keys["shape"].getIyy()) + dy**2)*w
            xy += (record.get(self.owner.keys["shape"].getIxy()) + dx*dy)*w
        self.parentRecord.set(self.owner.keys["shape"], lsst.afw.geom.Quadrupole(xx, yy, xy))
        # Run detection on the parent image to get the parent Footprint.
        self.owner._installFootprint(self.parentRecord, self.parentImage)
        # Create perfect HeavyFootprints for all children; these will need to
        # be modified later to account for the noise we'll add to the image.
        deblend = lsst.afw.image.MaskedImageF(self.owner.exposure.maskedImage, True)
        for record, image in self.children:
            deblend.image.array[:, :] = image.array
            heavyFootprint = lsst.afw.detection.HeavyFootprintF(self.parentRecord.getFootprint(), deblend)
            record.setFootprint(heavyFootprint)


class TestDataset:
    """A simulated dataset consisuting of test image and truth catalog.

    TestDataset creates an idealized image made of pure Gaussians (including a
    Gaussian PSF), with simple noise and idealized Footprints/HeavyFootprints
    that simulated the outputs of detection and deblending.  Multiple noise
    realizations can be created from the same underlying sources, allowing
    uncertainty estimates to be verified via Monte Carlo.

    Parameters
    ----------
    bbox : `lsst.geom.Box2I` or `lsst.geom.Box2D`
        Bounding box of the test image.
    threshold : `float`
        Threshold absolute value used to determine footprints for
        simulated sources.  This thresholding will be applied before noise is
        actually added to images (or before the noise level is even known), so
        this will necessarily produce somewhat artificial footprints.
    exposure : `lsst.afw.image.ExposureF`
        The image to which test sources should be added. Ownership should
        be considered transferred from the caller to the TestDataset.
        Must have a Gaussian PSF for truth catalog shapes to be exact.
    **kwds
        Keyword arguments forwarded to makeEmptyExposure if exposure is `None`.

    Notes
    -----
    Typical usage:

    .. code-block: py

        bbox = lsst.geom.Box2I(lsst.geom.Point2I(0,0), lsst.geom.Point2I(100,
                                                                         100))
        dataset = TestDataset(bbox)
        dataset.addSource(instFlux=1E5, centroid=lsst.geom.Point2D(25, 26))
        dataset.addSource(instFlux=2E5, centroid=lsst.geom.Point2D(75, 24),
                        shape=lsst.afw.geom.Quadrupole(8, 7, 2))
        with dataset.addBlend() as family:
            family.addChild(instFlux=2E5, centroid=lsst.geom.Point2D(50, 72))
            family.addChild(instFlux=1.5E5, centroid=lsst.geom.Point2D(51, 74))
        exposure, catalog = dataset.realize(noise=100.0,
                                            schema=TestDataset.makeMinimalSchema())
    """

    def __init__(self, bbox, threshold=10.0, exposure=None, **kwds):
        if exposure is None:
            exposure = self.makeEmptyExposure(bbox, **kwds)
        self.threshold = lsst.afw.detection.Threshold(threshold, lsst.afw.detection.Threshold.VALUE)
        self.exposure = exposure
        self.psfShape = self.exposure.getPsf().computeShape(bbox.getCenter())
        self.schema = self.makeMinimalSchema()
        self.catalog = lsst.afw.table.SourceCatalog(self.schema)

    @classmethod
    def makeMinimalSchema(cls):
        """Return the minimal schema needed to hold truth catalog fields.

        Notes
        -----
        When `TestDataset.realize` is called, the schema must include at least
        these fields.  Usually it will include additional fields for
        measurement algorithm outputs, allowing the same catalog to be used
        for both truth values (the fields from the minimal schema) and the
        measurements.
        """
        if not hasattr(cls, "_schema"):
            schema = lsst.afw.table.SourceTable.makeMinimalSchema()
            cls.keys = {}
            cls.keys["coordErr"] = lsst.afw.table.CoordKey.addErrorFields(schema)
            cls.keys["parent"] = schema.find("parent").key
            cls.keys["nChild"] = schema.addField("deblend_nChild", type=np.int32)
            cls.keys["instFlux"] = schema.addField("truth_instFlux", type=np.float64,
                                                   doc="true instFlux", units="count")
            cls.keys["instFluxErr"] = schema.addField("truth_instFluxErr", type=np.float64,
                                                      doc="true instFluxErr", units="count")
            cls.keys["centroid"] = lsst.afw.table.Point2DKey.addFields(
                schema, "truth", "true simulated centroid", "pixel"
            )
            cls.keys["centroid_sigma"] = lsst.afw.table.CovarianceMatrix2fKey.addFields(
                schema, "truth", ['x', 'y'], "pixel"
            )
            cls.keys["centroid_flag"] = schema.addField("truth_flag", type="Flag",
                                                        doc="set if the object is a star")
            cls.keys["shape"] = lsst.afw.table.QuadrupoleKey.addFields(
                schema, "truth", "true shape after PSF convolution", lsst.afw.table.CoordinateType.PIXEL
            )
            cls.keys["isStar"] = schema.addField("truth_isStar", type="Flag",
                                                 doc="set if the object is a star")
            cls.keys["negative"] = schema.addField("flags_negative", type="Flag",
                                                   doc="set if source was detected as significantly negative"
                                                   )
            schema.getAliasMap().set("slot_Shape", "truth")
            schema.getAliasMap().set("slot_Centroid", "truth")
            schema.getAliasMap().set("slot_ModelFlux", "truth")
            cls._schema = schema
        schema = lsst.afw.table.Schema(cls._schema)
        schema.disconnectAliases()
        return schema

    @staticmethod
    def makePerturbedWcs(oldWcs, minScaleFactor=1.2, maxScaleFactor=1.5,
                         minRotation=None, maxRotation=None,
                         minRefShift=None, maxRefShift=None,
                         minPixShift=2.0, maxPixShift=4.0, randomSeed=1):
        """Return a perturbed version of the input WCS.

        Create a new undistorted TAN WCS that is similar but not identical to
        another, with random scaling, rotation, and offset (in both pixel
        position and reference position).

        Parameters
        ----------
        oldWcs : `lsst.afw.geom.SkyWcs`
            The input WCS.
        minScaleFactor : `float`
            Minimum scale factor to apply to the input WCS.
        maxScaleFactor : `float`
            Maximum scale factor to apply to the input WCS.
        minRotation : `lsst.geom.Angle` or `None`
            Minimum rotation to apply to the input WCS. If `None`, defaults to
            30 degrees.
        maxRotation : `lsst.geom.Angle` or `None`
            Minimum rotation to apply to the input WCS. If `None`, defaults to
            60 degrees.
        minRefShift : `lsst.geom.Angle` or `None`
            Miniumum shift to apply to the input WCS reference value. If
            `None`, defaults to 0.5 arcsec.
        maxRefShift : `lsst.geom.Angle` or `None`
            Miniumum shift to apply to the input WCS reference value. If
            `None`, defaults to 1.0 arcsec.
        minPixShift : `float`
            Minimum shift to apply to the input WCS reference pixel.
        maxPixShift : `float`
            Maximum shift to apply to the input WCS reference pixel.
        randomSeed : `int`
            Random seed.

        Returns
        -------
        newWcs : `lsst.afw.geom.SkyWcs`
            A perturbed version of the input WCS.

        Notes
        -----
        The maximum and minimum arguments are interpreted as absolute values
        for a split range that covers both positive and negative values (as
        this method is used in testing, it is typically most important to
        avoid perturbations near zero). Scale factors are treated somewhat
        differently: the actual scale factor is chosen between
        ``minScaleFactor`` and ``maxScaleFactor`` OR (``1/maxScaleFactor``)
        and (``1/minScaleFactor``).

        The default range for rotation is 30-60 degrees, and the default range
        for reference shift is 0.5-1.0 arcseconds (these cannot be safely
        included directly as default values because Angle objects are
        mutable).

        The random number generator is primed with the seed given. If
        `None`, a seed is automatically chosen.
        """
        random_state = np.random.RandomState(randomSeed)
        if minRotation is None:
            minRotation = 30.0*lsst.geom.degrees
        if maxRotation is None:
            maxRotation = 60.0*lsst.geom.degrees
        if minRefShift is None:
            minRefShift = 0.5*lsst.geom.arcseconds
        if maxRefShift is None:
            maxRefShift = 1.0*lsst.geom.arcseconds

        def splitRandom(min1, max1, min2=None, max2=None):
            if min2 is None:
                min2 = -max1
            if max2 is None:
                max2 = -min1
            if random_state.uniform() > 0.5:
                return float(random_state.uniform(min1, max1))
            else:
                return float(random_state.uniform(min2, max2))
        # Generate random perturbations
        scaleFactor = splitRandom(minScaleFactor, maxScaleFactor, 1.0/maxScaleFactor, 1.0/minScaleFactor)
        rotation = splitRandom(minRotation.asRadians(), maxRotation.asRadians())*lsst.geom.radians
        refShiftRa = splitRandom(minRefShift.asRadians(), maxRefShift.asRadians())*lsst.geom.radians
        refShiftDec = splitRandom(minRefShift.asRadians(), maxRefShift.asRadians())*lsst.geom.radians
        pixShiftX = splitRandom(minPixShift, maxPixShift)
        pixShiftY = splitRandom(minPixShift, maxPixShift)
        # Compute new CD matrix
        oldTransform = lsst.geom.LinearTransform(oldWcs.getCdMatrix())
        rTransform = lsst.geom.LinearTransform.makeRotation(rotation)
        sTransform = lsst.geom.LinearTransform.makeScaling(scaleFactor)
        newTransform = oldTransform*rTransform*sTransform
        matrix = newTransform.getMatrix()
        # Compute new coordinate reference pixel (CRVAL)
        oldSkyOrigin = oldWcs.getSkyOrigin()
        newSkyOrigin = lsst.geom.SpherePoint(oldSkyOrigin.getRa() + refShiftRa,
                                             oldSkyOrigin.getDec() + refShiftDec)
        # Compute new pixel reference pixel (CRPIX)
        oldPixOrigin = oldWcs.getPixelOrigin()
        newPixOrigin = lsst.geom.Point2D(oldPixOrigin.getX() + pixShiftX,
                                         oldPixOrigin.getY() + pixShiftY)
        return lsst.afw.geom.makeSkyWcs(crpix=newPixOrigin, crval=newSkyOrigin, cdMatrix=matrix)

    @staticmethod
    def makeEmptyExposure(bbox, wcs=None, crval=None, cdelt=None, psfSigma=2.0, psfDim=17, calibration=4,
                          visitId=1234, mjd=60000.0, detector=1):
        """Create an Exposure, with a PhotoCalib, Wcs, and Psf, but no pixel values.

        Parameters
        ----------
        bbox : `lsst.geom.Box2I` or `lsst.geom.Box2D`
            Bounding box of the image in image coordinates.
        wcs : `lsst.afw.geom.SkyWcs`, optional
            New WCS for the exposure (created from CRVAL and CDELT if `None`).
        crval : `lsst.afw.geom.SpherePoint`, optional
            ICRS center of the TAN WCS attached to the image. If `None`, (45
            degrees, 45 degrees) is assumed.
        cdelt : `lsst.geom.Angle`, optional
            Pixel scale of the image. If `None`, 0.2 arcsec is assumed.
        psfSigma : `float`, optional
            Radius (sigma) of the Gaussian PSF attached to the image
        psfDim : `int`, optional
            Width and height of the image's Gaussian PSF attached to the image
        calibration : `float`, optional
            The spatially-constant calibration (in nJy/count) to set the
            PhotoCalib of the exposure.
        visitId : `int`, optional
            Visit id to store in VisitInfo.
        mjd : `float`, optional
            Modified Julian Date of this exposure to store in VisitInfo.
        detector: `int`, optional
            Detector id to assign to the attached Detector object.

        Returns
        -------
        exposure : `lsst.age.image.ExposureF`
            An empty image.
        """
        if wcs is None:
            if crval is None:
                crval = lsst.geom.SpherePoint(45.0, 45.0, lsst.geom.degrees)
            if cdelt is None:
                cdelt = 0.2*lsst.geom.arcseconds
            crpix = lsst.geom.Box2D(bbox).getCenter()
            wcs = lsst.afw.geom.makeSkyWcs(crpix=crpix, crval=crval,
                                           cdMatrix=lsst.afw.geom.makeCdMatrix(scale=cdelt))
        exposure = lsst.afw.image.ExposureF(bbox)
        psf = lsst.afw.detection.GaussianPsf(psfDim, psfDim, psfSigma)
        photoCalib = lsst.afw.image.PhotoCalib(calibration)
        visitInfo = lsst.afw.image.VisitInfo(id=visitId,
                                             exposureTime=30.0,
                                             date=lsst.daf.base.DateTime(mjd),
                                             observatory=lsst.afw.coord.Observatory(11.1*lsst.geom.degrees,
                                                                                    22.2*lsst.geom.degrees,
                                                                                    0.333),
                                             hasSimulatedContent=True)
        exposure.setWcs(wcs)
        exposure.setPsf(psf)
        exposure.setPhotoCalib(photoCalib)
        exposure.info.setVisitInfo(visitInfo)
        exposure.info.setDetector(lsst.afw.cameraGeom.testUtils.DetectorWrapper(id=detector).detector)
        return exposure

    @staticmethod
    def drawGaussian(bbox, instFlux, ellipse):
        """Create an image of an elliptical Gaussian.

        Parameters
        ----------
        bbox : `lsst.geom.Box2I` or `lsst.geom.Box2D`
            Bounding box of image to create.
        instFlux : `float`
            Total instrumental flux of the Gaussian (normalized analytically,
            not using pixel values).
        ellipse : `lsst.afw.geom.Ellipse`
            Defines the centroid and shape.

        Returns
        -------
        image : `lsst.afw.image.ImageF`
            An image of the Gaussian.
        """
        x, y = np.meshgrid(np.arange(bbox.getBeginX(), bbox.getEndX()),
                           np.arange(bbox.getBeginY(), bbox.getEndY()))
        t = ellipse.getGridTransform()
        xt = t[t.XX] * x + t[t.XY] * y + t[t.X]
        yt = t[t.YX] * x + t[t.YY] * y + t[t.Y]
        image = lsst.afw.image.ImageF(bbox)
        image.array[:, :] = np.exp(-0.5*(xt**2 + yt**2))*instFlux/(2.0*ellipse.getCore().getArea())
        return image

    def _installFootprint(self, record, image, setPeakSignificance=True, negative=False):
        """Create simulated Footprint and add it to a truth catalog record.
        """
        schema = lsst.afw.detection.PeakTable.makeMinimalSchema()
        if setPeakSignificance:
            schema.addField("significance", type=float,
                            doc="Ratio of peak value to configured standard deviation.")
        if negative:
            self.threshold.setPolarity(False)
        # Run detection on the single-source image
        fpSet = lsst.afw.detection.FootprintSet(image, self.threshold, peakSchema=schema)
        # the call below to the FootprintSet ctor is actually a grow operation
        fpSet = lsst.afw.detection.FootprintSet(fpSet, int(self.psfShape.getDeterminantRadius() + 1.0), True)
        if setPeakSignificance:
            # This isn't a traditional significance, since we're using the VALUE
            # threshold type, but it's the best we can do in that case.
            for footprint in fpSet.getFootprints():
                footprint.updatePeakSignificance(self.threshold.getValue())
        if negative:
            self.threshold.setPolarity(True)
        # Update the full exposure's mask plane to indicate the detection
        fpSet.setMask(self.exposure.mask, "DETECTED")
        # Attach the new footprint to the exposure
        if len(fpSet.getFootprints()) > 1:
            raise RuntimeError("Threshold value results in multiple Footprints for a single object")
        if len(fpSet.getFootprints()) == 0:
            raise RuntimeError("Threshold value results in zero Footprints for object")
        record.setFootprint(fpSet.getFootprints()[0])

    def addSource(self, instFlux, centroid, shape=None, setPeakSignificance=True,
                  negative=False):
        """Add a source to the simulation.

        To insert a point source with a given signal-to-noise (sn), the total
        ``instFlux`` should be: ``sn*noise*psf_scale``, where ``noise`` is the
        noise you will pass to ``realize()``, and
        ``psf_scale=sqrt(4*pi*r^2)``, where ``r`` is the width of the PSF.

        Parameters
        ----------
        instFlux : `float`
            Total instFlux of the source to be added.
        centroid : `lsst.geom.Point2D`
            Position of the source to be added.
        shape : `lsst.afw.geom.Quadrupole`
            Second moments of the source before PSF convolution. Note that the
            truth catalog records post-convolution moments. If `None`, a point
            source will be added.
        setPeakSignificance : `bool`
            Set the ``significance`` field for peaks in the footprints?
            See ``lsst.meas.algorithms.SourceDetectionTask.setPeakSignificance``
            for how this field is computed for real datasets.
        negative : `bool`
            Treat this as a negative source, using a negative polarity
            Threshold to create the footprint.

        Returns
        -------
        record : `lsst.afw.table.SourceRecord`
            A truth catalog record.
        image : `lsst.afw.image.ImageF`
            Single-source image corresponding to the new source.
        """
        # Create and set the truth catalog fields
        record = self.catalog.addNew()
        record.set(self.keys["instFlux"], instFlux)
        record.set(self.keys["instFluxErr"], 0)
        record.set(self.keys["centroid"], centroid)
        record.set(self.keys["negative"], negative)
        covariance = np.random.normal(0, 0.1, 4).reshape(2, 2)
        covariance[0, 1] = covariance[1, 0]  # CovarianceMatrixKey assumes symmetric x_y_Cov
        record.set(self.keys["centroid_sigma"], covariance.astype(np.float32))
        if shape is None:
            record.set(self.keys["isStar"], True)
            fullShape = self.psfShape
        else:
            record.set(self.keys["isStar"], False)
            fullShape = shape.convolve(self.psfShape)
        record.set(self.keys["shape"], fullShape)
        # Create an image containing just this source
        image = self.drawGaussian(self.exposure.getBBox(), instFlux,
                                  lsst.afw.geom.Ellipse(fullShape, centroid))
        # Generate a footprint for this source
        self._installFootprint(record, image, setPeakSignificance, negative=negative)
        # Actually add the source to the full exposure
        self.exposure.image.array[:, :] += image.array
        return record, image

    def addBlend(self):
        """Return a context manager which can add a blend of multiple sources.

        Notes
        -----
        Note that nothing stops you from creating overlapping sources just using the addSource() method,
        but addBlend() is necesssary to create a parent object and deblended HeavyFootprints of the type
        produced by the detection and deblending pipelines.

        Examples
        --------
        .. code-block: py
            d = TestDataset(...)
            with d.addBlend() as b:
                b.addChild(flux1, centroid1)
                b.addChild(flux2, centroid2, shape2)
        """
        return BlendContext(self)

    def transform(self, wcs, **kwds):
        """Copy this dataset transformed to a new WCS, with new Psf and PhotoCalib.

        Parameters
        ----------
        wcs : `lsst.afw.geom.SkyWcs`
            WCS for the new dataset.
        **kwds
            Additional keyword arguments passed on to
            `TestDataset.makeEmptyExposure`.  If not specified, these revert
            to the defaults for `~TestDataset.makeEmptyExposure`, not the
            values in the current dataset.

        Returns
        -------
        newDataset : `TestDataset`
            Transformed copy of this dataset.
        """
        bboxD = lsst.geom.Box2D()
        xyt = lsst.afw.geom.makeWcsPairTransform(self.exposure.getWcs(), wcs)
        for corner in lsst.geom.Box2D(self.exposure.getBBox()).getCorners():
            bboxD.include(xyt.applyForward(lsst.geom.Point2D(corner)))
        bboxI = lsst.geom.Box2I(bboxD)
        result = TestDataset(bbox=bboxI, wcs=wcs, **kwds)
        oldPhotoCalib = self.exposure.getPhotoCalib()
        newPhotoCalib = result.exposure.getPhotoCalib()
        oldPsfShape = self.exposure.getPsf().computeShape(bboxD.getCenter())
        for record in self.catalog:
            if record.get(self.keys["nChild"]):
                raise NotImplementedError("Transforming blended sources in TestDatasets is not supported")
            magnitude = oldPhotoCalib.instFluxToMagnitude(record.get(self.keys["instFlux"]))
            newFlux = newPhotoCalib.magnitudeToInstFlux(magnitude)
            oldCentroid = record.get(self.keys["centroid"])
            newCentroid = xyt.applyForward(oldCentroid)
            if record.get(self.keys["isStar"]):
                newDeconvolvedShape = None
            else:
                affine = lsst.afw.geom.linearizeTransform(xyt, oldCentroid)
                oldFullShape = record.get(self.keys["shape"])
                oldDeconvolvedShape = lsst.afw.geom.Quadrupole(
                    oldFullShape.getIxx() - oldPsfShape.getIxx(),
                    oldFullShape.getIyy() - oldPsfShape.getIyy(),
                    oldFullShape.getIxy() - oldPsfShape.getIxy(),
                    False
                )
                newDeconvolvedShape = oldDeconvolvedShape.transform(affine.getLinear())
            result.addSource(newFlux, newCentroid, newDeconvolvedShape)
        return result

    def realize(self, noise, schema, randomSeed=1):
        r"""Simulate an exposure and detection catalog for this dataset.

        The simulation includes noise, and the detection catalog includes
        `~lsst.afw.detection.heavyFootprint.HeavyFootprint`\ s.

        Parameters
        ----------
        noise : `float`
            Standard deviation of noise to be added to the exposure.  The
            noise will be Gaussian and constant, appropriate for the
            sky-limited regime.
        schema : `lsst.afw.table.Schema`
            Schema of the new catalog to be created.  Must start with
            ``self.schema`` (i.e. ``schema.contains(self.schema)`` must be
            `True`), but typically contains fields for already-configured
            measurement algorithms as well.
        randomSeed : `int`, optional
            Seed for the random number generator.
            If `None`, a seed is chosen automatically.

        Returns
        -------
        `exposure` : `lsst.afw.image.ExposureF`
            Simulated image.
        `catalog` : `lsst.afw.table.SourceCatalog`
            Simulated detection catalog.
        """
        random_state = np.random.RandomState(randomSeed)
        assert schema.contains(self.schema)
        mapper = lsst.afw.table.SchemaMapper(self.schema)
        mapper.addMinimalSchema(self.schema, True)
        exposure = self.exposure.clone()
        exposure.variance.array[:, :] = noise**2
        exposure.image.array[:, :] += random_state.randn(exposure.height, exposure.width)*noise
        catalog = lsst.afw.table.SourceCatalog(schema)
        catalog.extend(self.catalog, mapper=mapper)
        # Loop over sources and generate new HeavyFootprints that divide up
        # the noisy pixels, not the ideal no-noise pixels.
        for record in catalog:
            # parent objects have non-Heavy Footprints, which don't need to be
            # updated after adding noise.
            if record.getParent() == 0:
                continue
            # get flattened arrays that correspond to the no-noise and noisy
            # parent images
            parent = catalog.find(record.getParent())
            footprint = parent.getFootprint()
            parentFluxArrayNoNoise = np.zeros(footprint.getArea(), dtype=np.float32)
            footprint.spans.flatten(parentFluxArrayNoNoise, self.exposure.image.array, self.exposure.getXY0())
            parentFluxArrayNoisy = np.zeros(footprint.getArea(), dtype=np.float32)
            footprint.spans.flatten(parentFluxArrayNoisy, exposure.image.array, exposure.getXY0())
            oldHeavy = record.getFootprint()
            fraction = (oldHeavy.getImageArray() / parentFluxArrayNoNoise)
            # N.B. this isn't a copy ctor - it's a copy from a vanilla
            # Footprint, so it doesn't copy the arrays we don't want to
            # change, and hence we have to do that ourselves below.
            newHeavy = lsst.afw.detection.HeavyFootprintF(oldHeavy)
            newHeavy.getImageArray()[:] = parentFluxArrayNoisy*fraction
            newHeavy.getMaskArray()[:] = oldHeavy.getMaskArray()
            newHeavy.getVarianceArray()[:] = oldHeavy.getVarianceArray()
            record.setFootprint(newHeavy)
        lsst.afw.table.updateSourceCoords(exposure.wcs, catalog)
        return exposure, catalog


class AlgorithmTestCase:
    """Base class for tests of measurement tasks.
    """
    def makeSingleFrameMeasurementConfig(self, plugin=None, dependencies=()):
        """Create an instance of `SingleFrameMeasurementTask.ConfigClass`.

        Only the specified plugin and its dependencies will be run; the
        Centroid, Shape, and ModelFlux slots will be set to the truth fields
        generated by the `TestDataset` class.

        Parameters
        ----------
        plugin : `str`
            Name of measurement plugin to enable.
        dependencies : iterable of `str`, optional
            Names of dependencies of the measurement plugin.

        Returns
        -------
        config : `SingleFrameMeasurementTask.ConfigClass`
            The resulting task configuration.
        """
        config = SingleFrameMeasurementTask.ConfigClass()
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message="ignoreSlotPluginChecks", category=FutureWarning)
            config = SingleFrameMeasurementTask.ConfigClass(ignoreSlotPluginChecks=True)
        config.slots.centroid = "truth"
        config.slots.shape = "truth"
        config.slots.modelFlux = None
        config.slots.apFlux = None
        config.slots.psfFlux = None
        config.slots.gaussianFlux = None
        config.slots.calibFlux = None
        config.plugins.names = (plugin,) + tuple(dependencies)
        return config

    def makeSingleFrameMeasurementTask(self, plugin=None, dependencies=(), config=None, schema=None,
                                       algMetadata=None):
        """Create a configured instance of `SingleFrameMeasurementTask`.

        Parameters
        ----------
        plugin : `str`, optional
            Name of measurement plugin to enable. If `None`, a configuration
            must be supplied as the ``config`` parameter. If both are
            specified, ``config`` takes precedence.
        dependencies : iterable of `str`, optional
            Names of dependencies of the specified measurement plugin.
        config : `SingleFrameMeasurementTask.ConfigClass`, optional
            Configuration for the task. If `None`, a measurement plugin must
            be supplied as the ``plugin`` paramter. If both are specified,
            ``config`` takes precedence.
        schema : `lsst.afw.table.Schema`, optional
            Measurement table schema. If `None`, a default schema is
            generated.
        algMetadata : `lsst.daf.base.PropertyList`, optional
            Measurement algorithm metadata. If `None`, a default container
            will be generated.

        Returns
        -------
        task : `SingleFrameMeasurementTask`
            A configured instance of the measurement task.
        """
        if config is None:
            if plugin is None:
                raise ValueError("Either plugin or config argument must not be None")
            config = self.makeSingleFrameMeasurementConfig(plugin=plugin, dependencies=dependencies)
        if schema is None:
            schema = TestDataset.makeMinimalSchema()
            # Clear all aliases so only those defined by config are set.
            schema.setAliasMap(None)
        if algMetadata is None:
            algMetadata = lsst.daf.base.PropertyList()
        return SingleFrameMeasurementTask(schema=schema, algMetadata=algMetadata, config=config)

    def makeForcedMeasurementConfig(self, plugin=None, dependencies=()):
        """Create an instance of `ForcedMeasurementTask.ConfigClass`.

        In addition to the plugins specified in the plugin and dependencies
        arguments, the `TransformedCentroid` and `TransformedShape` plugins
        will be run and used as the centroid and shape slots; these simply
        transform the reference catalog centroid and shape to the measurement
        coordinate system.

        Parameters
        ----------
        plugin : `str`
            Name of measurement plugin to enable.
        dependencies : iterable of `str`, optional
            Names of dependencies of the measurement plugin.

        Returns
        -------
        config : `ForcedMeasurementTask.ConfigClass`
            The resulting task configuration.
        """

        config = ForcedMeasurementTask.ConfigClass()
        config.slots.centroid = "base_TransformedCentroid"
        config.slots.shape = "base_TransformedShape"
        config.slots.modelFlux = None
        config.slots.apFlux = None
        config.slots.psfFlux = None
        config.slots.gaussianFlux = None
        config.plugins.names = (plugin,) + tuple(dependencies) + ("base_TransformedCentroid",
                                                                  "base_TransformedShape")
        return config

    def makeForcedMeasurementTask(self, plugin=None, dependencies=(), config=None, refSchema=None,
                                  algMetadata=None):
        """Create a configured instance of `ForcedMeasurementTask`.

        Parameters
        ----------
        plugin : `str`, optional
            Name of measurement plugin to enable. If `None`, a configuration
            must be supplied as the ``config`` parameter. If both are
            specified, ``config`` takes precedence.
        dependencies : iterable of `str`, optional
            Names of dependencies of the specified measurement plugin.
        config : `SingleFrameMeasurementTask.ConfigClass`, optional
            Configuration for the task. If `None`, a measurement plugin must
            be supplied as the ``plugin`` paramter. If both are specified,
            ``config`` takes precedence.
        refSchema : `lsst.afw.table.Schema`, optional
            Reference table schema. If `None`, a default schema is
            generated.
        algMetadata : `lsst.daf.base.PropertyList`, optional
            Measurement algorithm metadata. If `None`, a default container
            will be generated.

        Returns
        -------
        task : `ForcedMeasurementTask`
            A configured instance of the measurement task.
        """
        if config is None:
            if plugin is None:
                raise ValueError("Either plugin or config argument must not be None")
            config = self.makeForcedMeasurementConfig(plugin=plugin, dependencies=dependencies)
        if refSchema is None:
            refSchema = TestDataset.makeMinimalSchema()
        if algMetadata is None:
            algMetadata = lsst.daf.base.PropertyList()
        return ForcedMeasurementTask(refSchema=refSchema, algMetadata=algMetadata, config=config)


class TransformTestCase:
    """Base class for testing measurement transformations.

    Notes
    -----
    We test both that the transform itself operates successfully (fluxes are
    converted to magnitudes, flags are propagated properly) and that the
    transform is registered as the default for the appropriate measurement
    algorithms.

    In the simple case of one-measurement-per-transformation, the developer
    need not directly write any tests themselves: simply customizing the class
    variables is all that is required. More complex measurements (e.g.
    multiple aperture fluxes) require extra effort.
    """
    name = "MeasurementTransformTest"
    """The name used for the measurement algorithm (str).

    Notes
    -----
    This determines the names of the fields in the resulting catalog. This
    default should generally be fine, but subclasses can override if
    required.
    """

    # These should be customized by subclassing.
    controlClass = None
    algorithmClass = None
    transformClass = None

    flagNames = ("flag",)
    """Flags which may be set by the algorithm being tested (iterable of `str`).
    """

    # The plugin being tested should be registered under these names for
    # single frame and forced measurement. Should be customized by
    # subclassing.
    singleFramePlugins = ()
    forcedPlugins = ()

    def setUp(self):
        bbox = lsst.geom.Box2I(lsst.geom.Point2I(0, 0), lsst.geom.Point2I(200, 200))
        self.calexp = TestDataset.makeEmptyExposure(bbox)
        self._setupTransform()

    def tearDown(self):
        del self.calexp
        del self.inputCat
        del self.mapper
        del self.transform
        del self.outputCat

    def _populateCatalog(self, baseNames):
        records = []
        for flagValue in (True, False):
            records.append(self.inputCat.addNew())
            for baseName in baseNames:
                for flagName in self.flagNames:
                    if records[-1].schema.join(baseName, flagName) in records[-1].schema:
                        records[-1].set(records[-1].schema.join(baseName, flagName), flagValue)
        self._setFieldsInRecords(records, baseName)

    def _checkOutput(self, baseNames):
        for inSrc, outSrc in zip(self.inputCat, self.outputCat):
            for baseName in baseNames:
                self._compareFieldsInRecords(inSrc, outSrc, baseName)
                for flagName in self.flagNames:
                    keyName = outSrc.schema.join(baseName, flagName)
                    if keyName in inSrc.schema:
                        self.assertEqual(outSrc.get(keyName), inSrc.get(keyName))
                    else:
                        self.assertFalse(keyName in outSrc.schema)

    def _runTransform(self, doExtend=True):
        if doExtend:
            self.outputCat.extend(self.inputCat, mapper=self.mapper)
        self.transform(self.inputCat, self.outputCat, self.calexp.getWcs(), self.calexp.getPhotoCalib())

    def testTransform(self, baseNames=None):
        """Test the transformation on a catalog containing random data.

        Parameters
        ----------
        baseNames : iterable of `str`
            Iterable of the initial parts of measurement field names.

        Notes
        -----
        We check that:

        - An appropriate exception is raised on an attempt to transform
          between catalogs with different numbers of rows;
        - Otherwise, all appropriate conversions are properly appled and that
          flags have been propagated.

        The ``baseNames`` argument requires some explanation. This should be
        an iterable of the leading parts of the field names for each
        measurement; that is, everything that appears before ``_instFlux``,
        ``_flag``, etc. In the simple case of a single measurement per plugin,
        this is simply equal to ``self.name`` (thus measurements are stored as
        ``self.name + "_instFlux"``, etc). More generally, the developer may
        specify whatever iterable they require. For example, to handle
        multiple apertures, we could have ``(self.name + "_0", self.name +
        "_1", ...)``.
        """
        baseNames = baseNames or [self.name]
        self._populateCatalog(baseNames)
        self.assertRaises(lsst.pex.exceptions.LengthError, self._runTransform, False)
        self._runTransform()
        self._checkOutput(baseNames)

    def _checkRegisteredTransform(self, registry, name):
        # If this is a Python-based transform, we can compare directly; if
        # it's wrapped C++, we need to compare the wrapped class.
        self.assertEqual(registry[name].PluginClass.getTransformClass(), self.transformClass)

    def testRegistration(self):
        """Test that the transformation is appropriately registered.
        """
        for pluginName in self.singleFramePlugins:
            self._checkRegisteredTransform(lsst.meas.base.SingleFramePlugin.registry, pluginName)
        for pluginName in self.forcedPlugins:
            self._checkRegisteredTransform(lsst.meas.base.ForcedPlugin.registry, pluginName)


class SingleFramePluginTransformSetupHelper:

    def _setupTransform(self):
        self.control = self.controlClass()
        inputSchema = lsst.afw.table.SourceTable.makeMinimalSchema()
        # Trick algorithms that depend on the slot centroid or alias into thinking they've been defined;
        # it doesn't matter for this test since we won't actually use the plugins for anything besides
        # defining the schema.
        inputSchema.getAliasMap().set("slot_Centroid", "dummy")
        inputSchema.getAliasMap().set("slot_Shape", "dummy")
        self.algorithmClass(self.control, self.name, inputSchema)
        inputSchema.getAliasMap().erase("slot_Centroid")
        inputSchema.getAliasMap().erase("slot_Shape")
        self.inputCat = lsst.afw.table.SourceCatalog(inputSchema)
        self.mapper = lsst.afw.table.SchemaMapper(inputSchema)
        self.transform = self.transformClass(self.control, self.name, self.mapper)
        self.outputCat = lsst.afw.table.BaseCatalog(self.mapper.getOutputSchema())


class ForcedPluginTransformSetupHelper:

    def _setupTransform(self):
        self.control = self.controlClass()
        inputMapper = lsst.afw.table.SchemaMapper(lsst.afw.table.SourceTable.makeMinimalSchema(),
                                                  lsst.afw.table.SourceTable.makeMinimalSchema())
        # Trick algorithms that depend on the slot centroid or alias into thinking they've been defined;
        # it doesn't matter for this test since we won't actually use the plugins for anything besides
        # defining the schema.
        inputMapper.editOutputSchema().getAliasMap().set("slot_Centroid", "dummy")
        inputMapper.editOutputSchema().getAliasMap().set("slot_Shape", "dummy")
        self.algorithmClass(self.control, self.name, inputMapper, lsst.daf.base.PropertyList())
        inputMapper.editOutputSchema().getAliasMap().erase("slot_Centroid")
        inputMapper.editOutputSchema().getAliasMap().erase("slot_Shape")
        self.inputCat = lsst.afw.table.SourceCatalog(inputMapper.getOutputSchema())
        self.mapper = lsst.afw.table.SchemaMapper(inputMapper.getOutputSchema())
        self.transform = self.transformClass(self.control, self.name, self.mapper)
        self.outputCat = lsst.afw.table.BaseCatalog(self.mapper.getOutputSchema())


class FluxTransformTestCase(TransformTestCase):

    def _setFieldsInRecords(self, records, name):
        for record in records:
            record[record.schema.join(name, 'instFlux')] = np.random.random()
            record[record.schema.join(name, 'instFluxErr')] = np.random.random()

        # Negative instFluxes should be converted to NaNs.
        assert len(records) > 1
        records[0][record.schema.join(name, 'instFlux')] = -1

    def _compareFieldsInRecords(self, inSrc, outSrc, name):
        instFluxName = inSrc.schema.join(name, 'instFlux')
        instFluxErrName = inSrc.schema.join(name, 'instFluxErr')
        if inSrc[instFluxName] > 0:
            mag = self.calexp.getPhotoCalib().instFluxToMagnitude(inSrc[instFluxName],
                                                                  inSrc[instFluxErrName])
            self.assertEqual(outSrc[outSrc.schema.join(name, 'mag')], mag.value)
            self.assertEqual(outSrc[outSrc.schema.join(name, 'magErr')], mag.error)
        else:
            # negative instFlux results in NaN magnitude, but can still have finite error
            self.assertTrue(np.isnan(outSrc[outSrc.schema.join(name, 'mag')]))
            if np.isnan(inSrc[instFluxErrName]):
                self.assertTrue(np.isnan(outSrc[outSrc.schema.join(name, 'magErr')]))
            else:
                mag = self.calexp.getPhotoCalib().instFluxToMagnitude(inSrc[instFluxName],
                                                                      inSrc[instFluxErrName])
                self.assertEqual(outSrc[outSrc.schema.join(name, 'magErr')], mag.error)


class CentroidTransformTestCase(TransformTestCase):

    def _setFieldsInRecords(self, records, name):
        for record in records:
            record[record.schema.join(name, 'x')] = np.random.random()
            record[record.schema.join(name, 'y')] = np.random.random()
            # Some algorithms set no errors; some set only sigma on x & y; some provide
            # a full covariance matrix. Set only those which exist in the schema.
            for fieldSuffix in ('xErr', 'yErr', 'x_y_Cov'):
                fieldName = record.schema.join(name, fieldSuffix)
                if fieldName in record.schema:
                    record[fieldName] = np.random.random()

    def _compareFieldsInRecords(self, inSrc, outSrc, name):
        centroidResultKey = CentroidResultKey(inSrc.schema[self.name])
        centroidResult = centroidResultKey.get(inSrc)

        coord = lsst.afw.table.CoordKey(outSrc.schema[self.name]).get(outSrc)
        coordTruth = self.calexp.getWcs().pixelToSky(centroidResult.getCentroid())
        self.assertEqual(coordTruth, coord)

        # If the centroid has an associated uncertainty matrix, the coordinate
        # must have one too, and vice versa.
        try:
            coordErr = lsst.afw.table.CovarianceMatrix2fKey(outSrc.schema[self.name],
                                                            ["ra", "dec"]).get(outSrc)
        except lsst.pex.exceptions.NotFoundError:
            self.assertFalse(centroidResultKey.getCentroidErr().isValid())
        else:
            transform = self.calexp.getWcs().linearizePixelToSky(coordTruth, lsst.geom.radians)
            coordErrTruth = np.dot(np.dot(transform.getLinear().getMatrix(),
                                          centroidResult.getCentroidErr()),
                                   transform.getLinear().getMatrix().transpose())
            np.testing.assert_array_almost_equal(np.array(coordErrTruth), coordErr)
