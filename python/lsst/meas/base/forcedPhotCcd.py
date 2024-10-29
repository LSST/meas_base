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

import pandas as pd
import numpy as np

import lsst.pex.config
import lsst.pex.exceptions
import lsst.pipe.base
import lsst.geom
import lsst.afw.detection
import lsst.afw.geom
import lsst.afw.image
import lsst.afw.table
import lsst.sphgeom

from lsst.pipe.base import PipelineTaskConnections
import lsst.pipe.base.connectionTypes as cT

import lsst.pipe.base as pipeBase
from lsst.skymap import BaseSkyMap

from .forcedMeasurement import ForcedMeasurementTask
from .applyApCorr import ApplyApCorrTask
from .catalogCalculation import CatalogCalculationTask
from ._id_generator import DetectorVisitIdGeneratorConfig

__all__ = (
    "ForcedPhotCcdConfig",
    "ForcedPhotCcdTask",
    "ForcedPhotCcdFromDataFrameTask",
    "ForcedPhotCcdFromDataFrameConfig",
)


class ForcedPhotCcdConnections(
    PipelineTaskConnections,
    dimensions=("instrument", "visit", "detector", "skymap", "tract"),
    defaultTemplates={"inputCoaddName": "deep", "inputName": "pvi"},
):
    inputSchema = cT.InitInput(
        doc="Schema for the input measurement catalogs.",
        name="{inputCoaddName}Coadd_ref_schema",
        storageClass="SourceCatalog",
    )
    outputSchema = cT.InitOutput(
        doc="Schema for the output forced measurement catalogs.",
        name="forced_src_schema",
        storageClass="SourceCatalog",
    )
    exposure = cT.Input(
        doc="Input exposure to perform photometry on.",
        name="{inputName}",
        storageClass="ExposureF",
        dimensions=["instrument", "visit", "detector"],
    )
    refCat = cT.Input(
        doc="Catalog of shapes and positions at which to force photometry.",
        name="{inputCoaddName}Coadd_ref",
        storageClass="SourceCatalog",
        dimensions=["skymap", "tract", "patch"],
        multiple=True,
        deferLoad=True,
    )
    skyMap = cT.Input(
        doc="SkyMap dataset that defines the coordinate system of the reference catalog.",
        name=BaseSkyMap.SKYMAP_DATASET_TYPE_NAME,
        storageClass="SkyMap",
        dimensions=["skymap"],
    )
    skyCorr = cT.Input(
        doc="Input Sky Correction to be subtracted from the exposure if doApplySkyCorr=True",
        name="skyCorr",
        storageClass="Background",
        dimensions=("instrument", "visit", "detector"),
    )
    visitSummary = cT.Input(
        doc="Input visit-summary catalog with updated calibration objects.",
        name="finalVisitSummary",
        storageClass="ExposureCatalog",
        dimensions=("instrument", "visit"),
    )
    measCat = cT.Output(
        doc="Output forced photometry catalog.",
        name="forced_src",
        storageClass="SourceCatalog",
        dimensions=["instrument", "visit", "detector", "skymap", "tract"],
    )

    def __init__(self, *, config=None):
        super().__init__(config=config)
        if not config.doApplySkyCorr:
            del self.skyCorr
        if not config.useVisitSummary:
            del self.visitSummary


class ForcedPhotCcdConfig(pipeBase.PipelineTaskConfig, pipelineConnections=ForcedPhotCcdConnections):
    """Config class for forced measurement driver task."""

    measurement = lsst.pex.config.ConfigurableField(
        target=ForcedMeasurementTask, doc="subtask to do forced measurement"
    )
    coaddName = lsst.pex.config.Field(
        doc="coadd name: typically one of deep or goodSeeing",
        dtype=str,
        default="deep",
    )
    doApCorr = lsst.pex.config.Field(
        dtype=bool, default=True, doc="Run subtask to apply aperture corrections"
    )
    applyApCorr = lsst.pex.config.ConfigurableField(
        target=ApplyApCorrTask, doc="Subtask to apply aperture corrections"
    )
    catalogCalculation = lsst.pex.config.ConfigurableField(
        target=CatalogCalculationTask, doc="Subtask to run catalogCalculation plugins on catalog"
    )
    doApplySkyCorr = lsst.pex.config.Field(
        dtype=bool,
        default=False,
        doc="Apply sky correction?",
    )
    useVisitSummary = lsst.pex.config.Field(
        dtype=bool,
        default=False,
        doc=(
            "Use updated WCS, PhotoCalib, ApCorr, and PSF from visit summary? "
            "This should be False if and only if the input image already has the best-available calibration "
            "objects attached."
        ),
    )
    includePhotoCalibVar = lsst.pex.config.Field(
        dtype=bool,
        default=False,
        doc="Add photometric calibration variance to warp variance plane?",
    )
    footprintSource = lsst.pex.config.ChoiceField(
        dtype=str,
        doc="Where to obtain footprints to install in the measurement catalog, prior to measurement.",
        allowed={
            "transformed": "Transform footprints from the reference catalog (downgrades HeavyFootprints).",
            "psf": (
                "Use the scaled shape of the PSF at the position of each source (does not generate "
                "HeavyFootprints)."
            ),
        },
        optional=True,
        default="transformed",
    )
    psfFootprintScaling = lsst.pex.config.Field(
        dtype=float,
        doc="Scaling factor to apply to the PSF shape when footprintSource='psf' (ignored otherwise).",
        default=3.0,
    )
    idGenerator = DetectorVisitIdGeneratorConfig.make_field()

    def setDefaults(self):
        # Docstring inherited.
        super().setDefaults()
        # Footprints here will not be entirely correct, so don't try to make
        # a biased correction for blended neighbors.
        self.measurement.doReplaceWithNoise = False
        # Only run a minimal set of plugins, as these measurements are only
        # needed for PSF-like sources.
        self.measurement.plugins.names = [
            "base_PixelFlags",
            "base_TransformedCentroid",
            "base_PsfFlux",
            "base_LocalBackground",
            "base_LocalPhotoCalib",
            "base_LocalWcs",
        ]
        self.measurement.slots.shape = None
        # Make catalogCalculation a no-op by default as no modelFlux is setup
        # by default in ForcedMeasurementTask.
        self.catalogCalculation.plugins.names = []


class ForcedPhotCcdTask(pipeBase.PipelineTask):
    """A pipeline task for performing forced measurement on CCD images.

    Parameters
    ----------
    refSchema : `lsst.afw.table.Schema`, optional
        The schema of the reference catalog, passed to the constructor of the
        references subtask. Optional, but must be specified if ``initInputs``
        is not; if both are specified, ``initInputs`` takes precedence.
    initInputs : `dict`
        Dictionary that can contain a key ``inputSchema`` containing the
        schema. If present will override the value of ``refSchema``.
    **kwargs
        Keyword arguments are passed to the supertask constructor.
    """

    ConfigClass = ForcedPhotCcdConfig
    _DefaultName = "forcedPhotCcd"
    dataPrefix = ""

    def __init__(self, refSchema=None, initInputs=None, **kwargs):
        super().__init__(**kwargs)

        if initInputs is not None:
            refSchema = initInputs["inputSchema"].schema

        if refSchema is None:
            raise ValueError("No reference schema provided.")

        self.makeSubtask("measurement", refSchema=refSchema)
        # It is necessary to get the schema internal to the forced measurement
        # task until such a time that the schema is not owned by the
        # measurement task, but is passed in by an external caller.
        if self.config.doApCorr:
            self.makeSubtask("applyApCorr", schema=self.measurement.schema)
        self.makeSubtask("catalogCalculation", schema=self.measurement.schema)
        self.outputSchema = lsst.afw.table.SourceCatalog(self.measurement.schema)

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        inputs = butlerQC.get(inputRefs)

        tract = butlerQC.quantum.dataId["tract"]
        skyMap = inputs.pop("skyMap")
        inputs["refWcs"] = skyMap[tract].getWcs()

        # Connections only exist if they are configured to be used.
        skyCorr = inputs.pop("skyCorr", None)

        inputs["exposure"] = self.prepareCalibratedExposure(
            inputs["exposure"],
            skyCorr=skyCorr,
            visitSummary=inputs.pop("visitSummary", None),
        )

        inputs["refCat"] = self.mergeAndFilterReferences(
            inputs["exposure"], inputs["refCat"], inputs["refWcs"]
        )

        if inputs["refCat"] is None:
            self.log.info(
                "No WCS for exposure %s.  No %s catalog will be written.",
                butlerQC.quantum.dataId,
                outputRefs.measCat.datasetType.name,
            )
        else:
            inputs["measCat"], inputs["exposureId"] = self.generateMeasCat(
                inputRefs.exposure.dataId, inputs["exposure"], inputs["refCat"], inputs["refWcs"]
            )
            self.attachFootprints(inputs["measCat"], inputs["refCat"], inputs["exposure"], inputs["refWcs"])
            outputs = self.run(**inputs)
            butlerQC.put(outputs, outputRefs)

    def prepareCalibratedExposure(self, exposure, skyCorr=None, visitSummary=None):
        """Prepare a calibrated exposure and apply external calibrations
        and sky corrections if so configured.

        Parameters
        ----------
        exposure : `lsst.afw.image.exposure.Exposure`
            Input exposure to adjust calibrations.
        skyCorr : `lsst.afw.math.backgroundList`, optional
            Sky correction frame to apply if doApplySkyCorr=True.
        visitSummary : `lsst.afw.table.ExposureCatalog`, optional
            Exposure catalog with update calibrations; any not-None calibration
            objects attached will be used.  These are applied first and may be
            overridden by other arguments.

        Returns
        -------
        exposure : `lsst.afw.image.exposure.Exposure`
            Exposure with adjusted calibrations.
        """
        detectorId = exposure.getInfo().getDetector().getId()

        if visitSummary is not None:
            row = visitSummary.find(detectorId)
            if row is None:
                raise RuntimeError(f"Detector id {detectorId} not found in visitSummary.")
            if (photoCalib := row.getPhotoCalib()) is not None:
                exposure.setPhotoCalib(photoCalib)
            if (skyWcs := row.getWcs()) is not None:
                exposure.setWcs(skyWcs)
            if (psf := row.getPsf()) is not None:
                exposure.setPsf(psf)
            if (apCorrMap := row.getApCorrMap()) is not None:
                exposure.info.setApCorrMap(apCorrMap)

        if skyCorr is not None:
            exposure.maskedImage -= skyCorr.getImage()

        return exposure

    def mergeAndFilterReferences(self, exposure, refCats, refWcs):
        """Filter reference catalog so that all sources are within the
        boundaries of the exposure.

        Parameters
        ----------
        exposure : `lsst.afw.image.exposure.Exposure`
            Exposure to generate the catalog for.
        refCats : sequence of `lsst.daf.butler.DeferredDatasetHandle`
            Handles for catalogs of shapes and positions at which to force
            photometry.
        refWcs : `lsst.afw.image.SkyWcs`
            Reference world coordinate system.

        Returns
        -------
        refSources : `lsst.afw.table.SourceCatalog`
            Filtered catalog of forced sources to measure.

        Notes
        -----
        The majority of this code is based on the methods of
        lsst.meas.algorithms.loadReferenceObjects.ReferenceObjectLoader

        """
        mergedRefCat = None

        # Step 1: Determine bounds of the exposure photometry will
        # be performed on.
        expWcs = exposure.getWcs()
        if expWcs is None:
            self.log.info("Exposure has no WCS.  Returning None for mergedRefCat.")
        else:
            expRegion = exposure.getBBox(lsst.afw.image.PARENT)
            expBBox = lsst.geom.Box2D(expRegion)
            expBoxCorners = expBBox.getCorners()
            expSkyCorners = [expWcs.pixelToSky(corner).getVector() for corner in expBoxCorners]
            expPolygon = lsst.sphgeom.ConvexPolygon(expSkyCorners)

            # Step 2: Filter out reference catalog sources that are
            # not contained within the exposure boundaries, or whose
            # parents are not within the exposure boundaries.  Note
            # that within a single input refCat, the parents always
            # appear before the children.
            for refCat in refCats:
                refCat = refCat.get()
                if mergedRefCat is None:
                    mergedRefCat = lsst.afw.table.SourceCatalog(refCat.table)
                    containedIds = {0}  # zero as a parent ID means "this is a parent"
                for record in refCat:
                    if (
                        expPolygon.contains(record.getCoord().getVector())
                        and record.getParent() in containedIds
                    ):
                        record.setFootprint(record.getFootprint())
                        mergedRefCat.append(record)
                        containedIds.add(record.getId())
            if mergedRefCat is None:
                raise RuntimeError("No reference objects for forced photometry.")
            mergedRefCat.sort(lsst.afw.table.SourceTable.getParentKey())
        return mergedRefCat

    def generateMeasCat(self, dataId, exposure, refCat, refWcs):
        """Generate a measurement catalog.

        Parameters
        ----------
        dataId : `lsst.daf.butler.DataCoordinate`
            Butler data ID for this image, with ``{visit, detector}`` keys.
        exposure : `lsst.afw.image.exposure.Exposure`
            Exposure to generate the catalog for.
        refCat : `lsst.afw.table.SourceCatalog`
            Catalog of shapes and positions at which to force photometry.
        refWcs : `lsst.afw.image.SkyWcs`
            Reference world coordinate system.
            This parameter is not currently used.

        Returns
        -------
        measCat : `lsst.afw.table.SourceCatalog`
            Catalog of forced sources to measure.
        expId : `int`
            Unique binary id associated with the input exposure
        """
        id_generator = self.config.idGenerator.apply(dataId)
        measCat = self.measurement.generateMeasCat(
            exposure, refCat, refWcs, idFactory=id_generator.make_table_id_factory()
        )
        return measCat, id_generator.catalog_id

    def run(self, measCat, exposure, refCat, refWcs, exposureId=None):
        """Perform forced measurement on a single exposure.

        Parameters
        ----------
        measCat : `lsst.afw.table.SourceCatalog`
            The measurement catalog, based on the sources listed in the
            reference catalog.
        exposure : `lsst.afw.image.Exposure`
            The measurement image upon which to perform forced detection.
        refCat : `lsst.afw.table.SourceCatalog`
            The reference catalog of sources to measure.
        refWcs : `lsst.afw.image.SkyWcs`
            The WCS for the references.
        exposureId : `int`
            Optional unique exposureId used for random seed in measurement
            task.

        Returns
        -------
        result : `lsst.pipe.base.Struct`
            Structure with fields:

            ``measCat``
                Catalog of forced measurement results
                (`lsst.afw.table.SourceCatalog`).
        """
        self.measurement.run(measCat, exposure, refCat, refWcs, exposureId=exposureId)
        if self.config.doApCorr:
            apCorrMap = exposure.getInfo().getApCorrMap()
            if apCorrMap is None:
                self.log.warning("Forced exposure image does not have valid aperture correction; skipping.")
            else:
                self.applyApCorr.run(
                    catalog=measCat,
                    apCorrMap=apCorrMap,
                )
        self.catalogCalculation.run(measCat)

        return pipeBase.Struct(measCat=measCat)

    def attachFootprints(self, sources, refCat, exposure, refWcs):
        """Attach footprints to blank sources prior to measurements.

        Notes
        -----
        `~lsst.afw.detection.Footprint` objects for forced photometry must
        be in the pixel coordinate system of the image being measured, while
        the actual detections may start out in a different coordinate system.

        Subclasses of this class may implement this method to define how
        those `~lsst.afw.detection.Footprint` objects should be generated.

        This default implementation transforms depends on the
        ``footprintSource`` configuration parameter.
        """
        if self.config.footprintSource == "transformed":
            return self.measurement.attachTransformedFootprints(sources, refCat, exposure, refWcs)
        elif self.config.footprintSource == "psf":
            return self.measurement.attachPsfShapeFootprints(
                sources, exposure, scaling=self.config.psfFootprintScaling
            )


class ForcedPhotCcdFromDataFrameConnections(
    PipelineTaskConnections,
    dimensions=("instrument", "visit", "detector", "skymap", "tract"),
    defaultTemplates={
        "inputCoaddName": "goodSeeing",
        "inputName": "pvi",
    },
):
    refCat = cT.Input(
        doc="Catalog of positions at which to force photometry.",
        name="{inputCoaddName}Diff_fullDiaObjTable",
        storageClass="DataFrame",
        dimensions=["skymap", "tract", "patch"],
        multiple=True,
        deferLoad=True,
    )
    exposure = cT.Input(
        doc="Input exposure to perform photometry on.",
        name="{inputName}",
        storageClass="ExposureF",
        dimensions=["instrument", "visit", "detector"],
    )
    skyCorr = cT.Input(
        doc="Input Sky Correction to be subtracted from the exposure if doApplySkyCorr=True",
        name="skyCorr",
        storageClass="Background",
        dimensions=("instrument", "visit", "detector"),
    )
    visitSummary = cT.Input(
        doc="Input visit-summary catalog with updated calibration objects.",
        name="finalVisitSummary",
        storageClass="ExposureCatalog",
        dimensions=("instrument", "visit"),
    )
    skyMap = cT.Input(
        doc="SkyMap dataset that defines the coordinate system of the reference catalog.",
        name=BaseSkyMap.SKYMAP_DATASET_TYPE_NAME,
        storageClass="SkyMap",
        dimensions=["skymap"],
    )
    measCat = cT.Output(
        doc="Output forced photometry catalog.",
        name="forced_src_diaObject",
        storageClass="SourceCatalog",
        dimensions=["instrument", "visit", "detector", "skymap", "tract"],
    )
    outputSchema = cT.InitOutput(
        doc="Schema for the output forced measurement catalogs.",
        name="forced_src_diaObject_schema",
        storageClass="SourceCatalog",
    )

    def __init__(self, *, config=None):
        super().__init__(config=config)
        if not config.doApplySkyCorr:
            del self.skyCorr
        if not config.useVisitSummary:
            del self.visitSummary


class ForcedPhotCcdFromDataFrameConfig(
    ForcedPhotCcdConfig, pipelineConnections=ForcedPhotCcdFromDataFrameConnections
):
    def setDefaults(self):
        super().setDefaults()
        self.footprintSource = "psf"
        self.measurement.doReplaceWithNoise = False
        # Only run a minimal set of plugins, as these measurements are only
        # needed for PSF-like sources.
        self.measurement.plugins.names = [
            "base_PixelFlags",
            "base_TransformedCentroidFromCoord",
            "base_PsfFlux",
            "base_LocalBackground",
            "base_LocalPhotoCalib",
            "base_LocalWcs",
        ]
        self.measurement.slots.shape = None
        # Make catalogCalculation a no-op by default as no modelFlux is setup
        # by default in ForcedMeasurementTask.
        self.catalogCalculation.plugins.names = []

        self.measurement.copyColumns = {"id": "diaObjectId", "coord_ra": "coord_ra", "coord_dec": "coord_dec"}
        self.measurement.slots.centroid = "base_TransformedCentroidFromCoord"
        self.measurement.slots.psfFlux = "base_PsfFlux"

    def validate(self):
        super().validate()
        if self.footprintSource == "transformed":
            raise ValueError(
                "Cannot transform footprints from reference catalog, "
                "because DataFrames can't hold footprints."
            )


class ForcedPhotCcdFromDataFrameTask(ForcedPhotCcdTask):
    """Force Photometry on a per-detector exposure with coords from a DataFrame

    Uses input from a DataFrame instead of SourceCatalog
    like the base class ForcedPhotCcd does.
    Writes out a SourceCatalog so that the downstream
    WriteForcedSourceTableTask can be reused with output from this Task.
    """

    _DefaultName = "forcedPhotCcdFromDataFrame"
    ConfigClass = ForcedPhotCcdFromDataFrameConfig

    def __init__(self, refSchema=None, initInputs=None, **kwargs):
        # Parent's init assumes that we have a reference schema; Cannot reuse
        pipeBase.PipelineTask.__init__(self, **kwargs)

        self.makeSubtask("measurement", refSchema=lsst.afw.table.SourceTable.makeMinimalSchema())

        if self.config.doApCorr:
            self.makeSubtask("applyApCorr", schema=self.measurement.schema)
        self.makeSubtask("catalogCalculation", schema=self.measurement.schema)
        self.outputSchema = lsst.afw.table.SourceCatalog(self.measurement.schema)

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        inputs = butlerQC.get(inputRefs)

        tract = butlerQC.quantum.dataId["tract"]
        skyMap = inputs.pop("skyMap")
        inputs["refWcs"] = skyMap[tract].getWcs()

        # Connections only exist if they are configured to be used.
        skyCorr = inputs.pop("skyCorr", None)

        inputs["exposure"] = self.prepareCalibratedExposure(
            inputs["exposure"],
            skyCorr=skyCorr,
            visitSummary=inputs.pop("visitSummary", None),
        )

        self.log.info("Filtering ref cats: %s", ",".join([str(i.dataId) for i in inputs["refCat"]]))
        if inputs["exposure"].getWcs() is not None:
            refCat = self.df2RefCat(
                [i.get(parameters={"columns": ["diaObjectId", "ra", "dec"]}) for i in inputs["refCat"]],
                inputs["exposure"].getBBox(),
                inputs["exposure"].getWcs(),
            )
            inputs["refCat"] = refCat
            # generateMeasCat does not use the refWcs.
            inputs["measCat"], inputs["exposureId"] = self.generateMeasCat(
                inputRefs.exposure.dataId, inputs["exposure"], inputs["refCat"], inputs["refWcs"]
            )
            # attachFootprints only uses refWcs in ``transformed`` mode, which is not
            # supported in the DataFrame-backed task.
            self.attachFootprints(inputs["measCat"], inputs["refCat"], inputs["exposure"], inputs["refWcs"])
            outputs = self.run(**inputs)

            butlerQC.put(outputs, outputRefs)
        else:
            self.log.info(
                "No WCS for %s.  Skipping and no %s catalog will be written.",
                butlerQC.quantum.dataId,
                outputRefs.measCat.datasetType.name,
            )

    def df2RefCat(self, dfList, exposureBBox, exposureWcs):
        """Convert list of DataFrames to reference catalog

        Concatenate list of DataFrames presumably from multiple patches and
        downselect rows that overlap the exposureBBox using the exposureWcs.

        Parameters
        ----------
        dfList : `list` of `pandas.DataFrame`
            Each element containst diaObjects with ra/dec position in degrees
            Columns 'diaObjectId', 'ra', 'dec' are expected
        exposureBBox :   `lsst.geom.Box2I`
            Bounding box on which to select rows that overlap
        exposureWcs : `lsst.afw.geom.SkyWcs`
            World coordinate system to convert sky coords in ref cat to
            pixel coords with which to compare with exposureBBox

        Returns
        -------
        refCat : `lsst.afw.table.SourceTable`
            Source Catalog with minimal schema that overlaps exposureBBox
        """
        df = pd.concat(dfList)
        # translate ra/dec coords in dataframe to detector pixel coords
        # to down select rows that overlap the detector bbox
        mapping = exposureWcs.getTransform().getMapping()
        x, y = mapping.applyInverse(np.array(df[["ra", "dec"]].values * 2 * np.pi / 360).T)
        inBBox = lsst.geom.Box2D(exposureBBox).contains(x, y)
        refCat = self.df2SourceCat(df[inBBox])
        return refCat

    def df2SourceCat(self, df):
        """Create minimal schema SourceCatalog from a pandas DataFrame.

        The forced measurement subtask expects this as input.

        Parameters
        ----------
        df : `pandas.DataFrame`
            DiaObjects with locations and ids.

        Returns
        -------
        outputCatalog : `lsst.afw.table.SourceTable`
            Output catalog with minimal schema.
        """
        schema = lsst.afw.table.SourceTable.makeMinimalSchema()
        outputCatalog = lsst.afw.table.SourceCatalog(schema)
        outputCatalog.reserve(len(df))

        for diaObjectId, ra, dec in df[["ra", "dec"]].itertuples():
            outputRecord = outputCatalog.addNew()
            outputRecord.setId(diaObjectId)
            outputRecord.setCoord(lsst.geom.SpherePoint(ra, dec, lsst.geom.degrees))
        return outputCatalog
