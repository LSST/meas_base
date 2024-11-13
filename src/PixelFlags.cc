// -*- lsst-c++ -*-
/*
 * LSST Data Management System
 * Copyright 2008-2013 LSST Corporation.
 *
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the LSST License Statement and
 * the GNU General Public License along with this program.  If not,
 * see <http://www.lsstcorp.org/LegalNotices/>.
 */

#include <cctype>     // ::tolower
#include <algorithm>  // std::transform
#include <cmath>

#include "ndarray/eigen.h"

#include "lsst/geom/Box.h"
#include "lsst/geom/Point.h"
#include "lsst/afw/detection/Psf.h"
#include "lsst/afw/table/Source.h"
#include "lsst/afw/geom/SpanSet.h"
#include "lsst/meas/base/PixelFlags.h"

namespace lsst {
namespace meas {
namespace base {
namespace {
template <typename MaskedImageT>
class FootprintBits {
public:
    explicit FootprintBits() : _anyBits(0), _allBits(~static_cast<typename MaskedImageT::Mask::Pixel>(0x0)) {}

    /// \brief Reset everything for a new Footprint
    void reset() {
        _anyBits = 0x0;
        _allBits = ~static_cast<typename MaskedImageT::Mask::Pixel>(0x0);
    }

    void operator()(geom::Point2I const& point, typename MaskedImageT::Mask::Pixel const& value) {
        _anyBits |= value;
        _allBits &= value;
    }

    /// Return the union of the bits set anywhere in the Footprint.
    typename MaskedImageT::Mask::Pixel getAnyBits() const { return _anyBits; }
    /// Return the union of the bits set everywhere in the Footprint.
    typename MaskedImageT::Mask::Pixel getAllBits() const { return _allBits; }

private:
    typename MaskedImageT::Mask::Pixel _anyBits;
    typename MaskedImageT::Mask::Pixel _allBits;
};

typedef afw::image::MaskedImage<float> MaskedImageF;

// Set flags when any pixel in func has the mask bit set.
void updateFlags(PixelFlagsAlgorithm::KeyMap const& maskFlagToPixelFlag,
                 const FootprintBits<MaskedImageF>& func, afw::table::SourceRecord& measRecord) {
    for (auto const& i : maskFlagToPixelFlag) {
        try {
            if (func.getAnyBits() & MaskedImageF::Mask::getPlaneBitMask(i.first)) {
                measRecord.set(i.second, true);
            }
        } catch (pex::exceptions::InvalidParameterError& err) {
            throw LSST_EXCEPT(FatalAlgorithmError, err.what());
        }
    }
}

// Set flags when all pixels in func have the mask bit set.
void updateFlagsAll(PixelFlagsAlgorithm::KeyMap const& maskFlagToPixelFlag,
                    const FootprintBits<MaskedImageF>& func, afw::table::SourceRecord& measRecord) {
    for (auto const& i : maskFlagToPixelFlag) {
        try {
            if (func.getAllBits() & MaskedImageF::Mask::getPlaneBitMask(i.first)) {
                measRecord.set(i.second, true);
            }
        } catch (pex::exceptions::InvalidParameterError& err) {
            throw LSST_EXCEPT(FatalAlgorithmError, err.what());
        }
    }
}

}  // end anonymous namespace

PixelFlagsAlgorithm::PixelFlagsAlgorithm(Control const& ctrl, std::string const& name,
                                         afw::table::Schema& schema)
        : _ctrl(ctrl) {
    // Add generic keys first, which don't correspond to specific mask planes
    _generalFailureKey = schema.addField<afw::table::Flag>(
            name + "_flag", "General failure flag, set if anything went wrong setting flags.");
    _offImageKey =
            schema.addField<afw::table::Flag>(name + "_flag" + "_offimage", "Source center is off image");
    // Set all the flags that correspond to mask planes anywhere in the footprint
    _anyKeys["EDGE"] = schema.addField<afw::table::Flag>(
            name + "_flag_edge",
            "Pixel in source outside usable exposure region (masked EDGE or centroid off image).");
    _anyKeys["NO_DATA"] = schema.addField<afw::table::Flag>(name + "_flag_nodata",
                                                            "NO_DATA pixel in the source footprint.");
    _anyKeys["INTRP"] = schema.addField<afw::table::Flag>(name + "_flag_interpolated",
                                                          "Interpolated pixel in the Source footprint");
    _anyKeys["SAT"] = schema.addField<afw::table::Flag>(name + "_flag_saturated",
                                                        "Saturated pixel in the Source footprint");
    _anyKeys["CR"] =
            schema.addField<afw::table::Flag>(name + "_flag_cr", "Cosmic ray in the Source footprint");
    _anyKeys["BAD"] =
            schema.addField<afw::table::Flag>(name + "_flag_bad", "Bad pixel in the Source footprint");
    _anyKeys["SUSPECT"] = schema.addField<afw::table::Flag>(name + "_flag_suspect",
                                                            "Source's footprint includes suspect pixels");
    // Flags that correspond to mask bits which are set anywhere in the 3x3 central region of the object.
    _centerKeys["EDGE"] = schema.addField<afw::table::Flag>(
            name + "_flag_edgeCenter", "EDGE pixel in the 3x3 region around the centroid.");
    _centerKeys["NO_DATA"] = schema.addField<afw::table::Flag>(
            name + "_flag_nodataCenter", "NO_DATA pixel in the 3x3 region around the centroid.");
    _centerKeys["INTRP"] = schema.addField<afw::table::Flag>(
            name + "_flag_interpolatedCenter", "Interpolated pixel in the 3x3 region around the centroid.");
    _centerKeys["SAT"] = schema.addField<afw::table::Flag>(
            name + "_flag_saturatedCenter", "Saturated pixel in the 3x3 region around the centroid.");
    _centerKeys["CR"] = schema.addField<afw::table::Flag>(
            name + "_flag_crCenter", "Cosmic ray in the 3x3 region around the centroid.");
    _centerKeys["BAD"] = schema.addField<afw::table::Flag>(name + "_flag_badCenter",
                                                           "Bad pixel in the 3x3 region around the centroid");
    _centerKeys["SUSPECT"] = schema.addField<afw::table::Flag>(
            name + "_flag_suspectCenter", "Suspect pixel in the 3x3 region around the centroid.");

    // Flags that correspond to mask bits which are set on all of the 3x3 central pixels of the object.
    _centerAllKeys["EDGE"] = schema.addField<afw::table::Flag>(
            name + "_flag_edgeCenterAll",
            "All pixels in the 3x3 region around the centroid are marked EDGE.");
    _centerAllKeys["NO_DATA"] = schema.addField<afw::table::Flag>(
            name + "_flag_nodataCenterAll",
            "All pixels in the 3x3 region around the centroid are marked NO_DATA");
    _centerAllKeys["INTRP"] = schema.addField<afw::table::Flag>(
            name + "_flag_interpolatedCenterAll",
            "All pixels in the 3x3 region around the centroid are interpolated.");
    _centerAllKeys["SAT"] = schema.addField<afw::table::Flag>(
            name + "_flag_saturatedCenterAll",
            "All pixels in the 3x3 region around the centroid are saturated.");
    _centerAllKeys["CR"] = schema.addField<afw::table::Flag>(
            name + "_flag_crCenterAll",
            "All pixels in the 3x3 region around the centroid have the cosmic ray mask bit.");
    _centerAllKeys["BAD"] = schema.addField<afw::table::Flag>(
            name + "_flag_badCenterAll", "All pixels in the 3x3 region around the centroid are bad.");
    _centerAllKeys["SUSPECT"] = schema.addField<afw::table::Flag>(
            name + "_flag_suspectCenterAll", "All pixels in the 3x3 region around the centroid are suspect.");

    // Read in the flags passed from the configuration, and add them to the schema
    for (auto const& i : _ctrl.masksFpCenter) {
        std::string maskName(i);
        std::transform(maskName.begin(), maskName.end(), maskName.begin(), ::tolower);
        _centerKeys[i] = schema.addField<afw::table::Flag>(
                name + "_flag_" + maskName + "Center", "3x3 region around the centroid has " + i + " pixels");
    }
    for (auto const& i : _ctrl.masksFpCenter) {
        std::string maskName(i);
        std::transform(maskName.begin(), maskName.end(), maskName.begin(), ::tolower);
        _centerAllKeys[i] = schema.addField<afw::table::Flag>(
                name + "_flag_" + maskName + "CenterAll",
                "All pixels in the 3x3 region around the source centroid are " + i + " pixels");
    }

    for (auto const& i : _ctrl.masksFpAnywhere) {
        std::string maskName(i);
        std::transform(maskName.begin(), maskName.end(), maskName.begin(), ::tolower);
        _anyKeys[i] = schema.addField<afw::table::Flag>(name + "_flag_" + maskName,
                                                        "Source footprint includes " + i + " pixels");
    }
}

void PixelFlagsAlgorithm::measure(afw::table::SourceRecord& measRecord,
                                  afw::image::Exposure<float> const& exposure) const {
    MaskedImageF mimage = exposure.getMaskedImage();
    FootprintBits<MaskedImageF> func;

    // Check if the measRecord has a valid centroid key, i.e. it was centroided
    geom::Point2D center;
    if (measRecord.getTable()->getCentroidSlot().getMeasKey().isValid()) {
        center = measRecord.getCentroid();
    } else {
        // Set the general failure flag because using the Peak might affect
        // the current measurement
        measRecord.set(_generalFailureKey, true);
        // Attempting to set pixel flags with no centroider, the center will
        // be determined though the peak pixel. The first peak in the
        // footprint (supplied by the measurement framework) should be the
        // highest peak so that one will be used as a proxy for the central
        // tendency of the distribution of flux for the record.
        std::shared_ptr<afw::detection::Footprint> footprint = measRecord.getFootprint();
        // If there is no footprint or the footprint contains no peaks, throw
        // a runtime error.
        if (!footprint || footprint->getPeaks().empty()) {
            throw LSST_EXCEPT(pex::exceptions::RuntimeError, "No footprint, or no footprint peaks detected");
        } else {
            center.setX(footprint->getPeaks().front().getFx());
            center.setY(footprint->getPeaks().front().getFy());
        }
    }

    // Flag centroids off the image.
    // Float bbox to properly handle non-finite centroids (casting Point2D->Point2I breaks NaN/inf).
    geom::Box2D bbox(mimage.getBBox());
    if (!bbox.contains(center)) {
        measRecord.set(_offImageKey, true);
        measRecord.set(_anyKeys.at("EDGE"), true);
        measRecord.set(_centerKeys.at("EDGE"), true);
        measRecord.set(_centerAllKeys.at("EDGE"), true);
    }

    // Check for bits set in the source's Footprint
    auto footprint = measRecord.getFootprint();
    if (!footprint) {
        throw LSST_EXCEPT(pex::exceptions::RuntimeError,
                          (boost::format("Source id %d has no footprint.") % measRecord.getId()).str());
    }
    auto fullSpans = footprint->getSpans();
    if (!fullSpans) {
        throw LSST_EXCEPT(
                pex::exceptions::RuntimeError,
                (boost::format("Source id %d has no spans in footprint.") % measRecord.getId()).str());
    }
    fullSpans->clippedTo(mimage.getBBox())->applyFunctor(func, *(mimage.getMask()));

    // update the source record for the any keys
    updateFlags(_anyKeys, func, measRecord);

    if (!(std::isfinite(center.getX()) && std::isfinite(center.getY()))) {
        auto msg =
                (boost::format("Centroid of source id %d passed to PixelFlags is non-finite; "
                               "footprint-based flags have been set, but centroid-based flags will not be.") %
                 measRecord.getId())
                        .str();
        throw LSST_EXCEPT(pex::exceptions::RuntimeError, msg);
    } else {
        // Check for bits set in the 3x3 box around the center
        geom::Point2I llc(afw::image::positionToIndex(center.getX()) - 1,
                          afw::image::positionToIndex(center.getY()) - 1);

        func.reset();
        auto spans = std::make_shared<afw::geom::SpanSet>(geom::Box2I(llc, geom::ExtentI(3)));
        afw::detection::Footprint const middle(spans);  // central 3x3
        middle.getSpans()->clippedTo(mimage.getBBox())->applyFunctor(func, *(mimage.getMask()));

        // Update the flags which have to do with the center of the footprint
        updateFlags(_centerKeys, func, measRecord);
        updateFlagsAll(_centerAllKeys, func, measRecord);
    }
}

void PixelFlagsAlgorithm::fail(afw::table::SourceRecord& measRecord, MeasurementError* error) const {
    measRecord.set(_generalFailureKey, true);
}

}  // namespace base
}  // namespace meas
}  // namespace lsst
