// -*- LSST-C++ -*-
/*
 * This file is part of meas_base.
 *
 * Developed for the LSST Data Management System.
 * This product includes software developed by the LSST Project
 * (https://www.lsst.org).
 * See the COPYRIGHT file at the top-level directory of this distribution
 * for details of code ownership.
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
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef LSST_MEAS_BASE_SillyCentroid_h_INCLUDED
#define LSST_MEAS_BASE_SillyCentroid_h_INCLUDED

#include "lsst/geom/Point.h"
#include "lsst/afw/table.h"
#include "lsst/meas/base.h"
#include "lsst/pex/config.h"

namespace test { namespace foo { namespace bar {
/**
 *  @brief A C++ control class to handle SillyCentroidAlgorithm's configuration
 */
class SillyCentroidControl {
public:
    LSST_CONTROL_FIELD(param, double, "Value to offset the centroid position in both x and y");

    /**
     *  @brief Default constructor
     *
     *  All control classes should define a default constructor that sets all fields to their default values.
     */
    SillyCentroidControl() : param(0.0) {}
};

/**
 *  @brief A class that calculates a centroid as a simple unweighted first moment
 *         of the 3x3 region around a pixel; this should only be used in tests.
 *
 *   A fixed background (set via config) may optionally be subtracted.
 *   This algorithm does not report an error.
 */

class SillyCentroidAlgorithm : public lsst::meas::base::SimpleAlgorithm {
public:

    /**
     *  @brief Flag bits to be used with the 'flags' data member of the Result object.
     *
     *  Inspect getFlagDefinitions() for more detailed explanations of each flag.
     */
    typedef SillyCentroidControl Control;

    SillyCentroidAlgorithm(
        Control const & ctrl,
        std::string const & name,
        lsst::afw::table::Schema & schema
    ) : _ctrl(ctrl),
        _centroidKey(
            lsst::meas::base::CentroidResultKey::addFields(schema, name,
            "centroid from Silly Centroid algorithm", lsst::meas::base::SIGMA_ONLY)
        ),
        _centroidExtractor(schema, name, true)
    {
        static lsst::meas::base::FlagDefinitionList flagDefs = {{
            {"flag", "general failure flag, set if anything went wrong"},
            {"flag_noCounts", "Object to be centroided has no counts"},
            {"flag_edge", "Object too close to edge"}
        }};
        _flagHandler = lsst::meas::base::FlagHandler::addFields(schema, name, flagDefs);
    }
    
    void measure(
        lsst::afw::table::SourceRecord & measRecord,
        lsst::afw::image::Exposure<float> const & exposure
    ) const {
        
        lsst::geom::Point2D center = _centroidExtractor(measRecord, _flagHandler);
        lsst::meas::base::CentroidResult result;
        result.x = center.getX() + _ctrl.param;
        result.y = center.getY() + _ctrl.param; 
        measRecord.set(_centroidKey, result); // better than NaN
    
        measRecord.set(_centroidKey, result);
        _flagHandler.setValue(measRecord, "flag", false);  // if we had a suspect flag, we'd set that instead
    }

    void fail(lsst::afw::table::SourceRecord & measRecord, lsst::meas::base::MeasurementError * error) const {
        _flagHandler.handleFailure(measRecord, error);
    }

private:

    Control _ctrl;
    
    lsst::meas::base::CentroidResultKey _centroidKey;
    lsst::meas::base::FlagHandler _flagHandler;
    lsst::meas::base::SafeCentroidExtractor _centroidExtractor;
};

/**
 *  @brief A trivial example of a transformation.
 */
class SillyTransform :
    public lsst::meas::base::BaseTransform
{
public:
    typedef SillyCentroidControl Control;

    SillyTransform(Control const & ctrl, std::string const & name, lsst::afw::table::SchemaMapper & mapper)
                   : BaseTransform(name), _ctrl(ctrl) {
        // Map these fields from the input to the output
        lsst::afw::table::Key<double> xkey = mapper.getInputSchema()[name + "_x"];
        lsst::afw::table::Key<double> ykey = mapper.getInputSchema()[name + "_y"];
        mapper.addMapping(xkey);
        mapper.addMapping(ykey);

        // Add these fields to the output schema to fill in later; the keys
        // are stored as members of this class
        _key_revX = mapper.editOutputSchema().addField<double>(name + "_reverse_x", "reversed X");
        _key_revY = mapper.editOutputSchema().addField<double>(name + "_reverse_y", "reversed Y");
    }

    virtual void operator()(lsst::afw::table::SourceCatalog const& inputCatalog,
                            lsst::afw::table::BaseCatalog& outputCatalog, lsst::afw::geom::SkyWcs const& wcs,
                            lsst::afw::image::PhotoCalib const& photoCalib) const {
        checkCatalogSize(inputCatalog, outputCatalog);
        lsst::afw::table::Key<double> xkey = inputCatalog.getSchema()[_name + "_x"];
        lsst::afw::table::Key<double> ykey = inputCatalog.getSchema()[_name + "_y"];
        lsst::afw::table::SourceCatalog::const_iterator inSrc = inputCatalog.begin();
        lsst::afw::table::BaseCatalog::iterator outSrc = outputCatalog.begin();
        for (; inSrc < inputCatalog.end() && outSrc < outputCatalog.end(); ++inSrc, ++outSrc) {
            // Store the "reversed" versions of the x and y positions in the
            // output catalog.
            outSrc->set(_key_revX, -1.0 * inSrc->get(xkey));
            outSrc->set(_key_revY, -1.0 * inSrc->get(ykey));
        }
    }  // operator()

private:
    Control _ctrl;
    lsst::afw::table::Key<double> _key_revX;
    lsst::afw::table::Key<double> _key_revY;
};

}  // namespace bar
}  // namespace foo
}  // namespace test
#endif  // !LSST_MEAS_BASE_SillyCentroid_h_INCLUDED
