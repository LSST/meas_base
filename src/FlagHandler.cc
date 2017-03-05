// -*- lsst-c++ -*-
/*
 * LSST Data Management System
 * Copyright 2008-2014 LSST Corporation.
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

#include "lsst/meas/base/FlagHandler.h"

namespace lsst { namespace meas { namespace base {

// so pybind11 can get the address of number_undefined
constexpr std::size_t FlagDefinition::number_undefined;

FlagDefinition FlagDefinitionList::addFailureFlag(
    std::string const & doc
) {
        return add(FlagHandler::getFailureFlagName(), doc);
}

FlagHandler FlagHandler::addFields(
    afw::table::Schema & schema,
    std::string const & prefix,
    FlagDefinitionList const & flagDefs,
    FlagDefinitionList const & exclDefs
) {
    FlagHandler r;
    r._vector.reserve(flagDefs.size());
    for (std::size_t i = 0; i < flagDefs.size(); i++) {
        FlagDefinition const & flagDef = flagDefs[i];
        if (exclDefs.hasDefinition(flagDef.name)) {
            afw::table::Key<afw::table::Flag> key;
            r._vector.push_back( std::make_pair( flagDef.name, key));
        }
        else {
            afw::table::Key<afw::table::Flag> key(schema.addField<afw::table::Flag>(schema.join(prefix, flagDef.name), flagDef.doc));
            r._vector.push_back( std::make_pair( flagDef.name, key));
            if (flagDef.name == FlagHandler::getFailureFlagName()) {
                r.failureFlagNumber = i;
     }
        }
    }
    return r;
}

FlagHandler::FlagHandler(
    afw::table::SubSchema const & s,
    FlagDefinitionList const & flagDefs,
    FlagDefinitionList const & exclDefs
) : failureFlagNumber(FlagDefinition::number_undefined) {
    _vector.reserve(flagDefs.size());
    for (std::size_t i = 0; i < flagDefs.size(); i++ ) {
        FlagDefinition const & flagDef = flagDefs[i];
        if (exclDefs.hasDefinition(flagDef.name)) {
            afw::table::Key<afw::table::Flag> key;
            _vector.push_back(
                std::make_pair(
                    flagDef.name,
                    key
                )
            );
        }
        else {
            _vector.push_back(
                std::make_pair(
                    flagDef.name,
                    s[flagDef.name]
                )
            );
            if (flagDef.name == FlagHandler::getFailureFlagName()) {
                failureFlagNumber = i;
            }
        }
    }
}

void FlagHandler::handleFailure(afw::table::BaseRecord & record, MeasurementError const * error) const {
    std::size_t const numFlags = _vector.size();
    if (failureFlagNumber != FlagDefinition::number_undefined) {
        record.set(_vector[failureFlagNumber].second, true);
    }
    if (error && error->getFlagBit() != FlagDefinition::number_undefined) {
        assert(numFlags > error->getFlagBit());  // We need the particular flag
        record.set(_vector[error->getFlagBit()].second, true);
    }
}

}}} // lsst::meas::base
