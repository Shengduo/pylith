// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2017 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "SlipWeakeningFH.hh" // implementation of object methods

#include "pylith/materials/Metadata.hh" // USES Metadata

#include "pylith/utils/array.hh" // USES scalar_array
#include "pylith/utils/constdefs.h" // USES MAXSCALAR

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "petsc.h" // USES PetscLogFlops

#include <cassert> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error
// ----------------------------------------------------------------------
namespace pylith {
  namespace friction {
    namespace _SlipWeakeningFH {

      // Number of physical properties.
      const int numProperties = 6;

      // Physical properties.
      const pylith::materials::Metadata::ParamDescription properties[] = {
	{ "static_coefficient", 1, pylith::topology::FieldBase::SCALAR },
	{ "dynamic_coefficient", 1, pylith::topology::FieldBase::SCALAR },
        { "slip_weakening_parameter", 1, pylith::topology::FieldBase::SCALAR },
	{ "cohesion", 1, pylith::topology::FieldBase::SCALAR },
  { "fw", 1, pylith::topology::FieldBase::SCALAR },
  { "Vw", 1, pylith::topology::FieldBase::SCALAR },
      };

      // Number of State Variables.
      const int numStateVars = 2;

      // State Variables.
      const pylith::materials::Metadata::ParamDescription stateVars[] = {
	{ "cumulative_slip", 1, pylith::topology::FieldBase::SCALAR },
	{ "previous_slip", 1, pylith::topology::FieldBase::SCALAR },
      };

      // Values expected in spatial database
      const int numDBProperties = 6;
      const char* dbProperties[6] = { "static-coefficient",
				      "dynamic-coefficient",
				      "slip-weakening-parameter",
				      "cohesion",
              "fw", 
              "Vw", 
      };

      const int numDBStateVars = 2;
      const char* dbStateVars[2] = { "cumulative-slip",
				     "previous-slip",
      };      
      
    } // _SlipWeakeningFH
  } // friction
} // pylith

// Indices of physical properties
const int pylith::friction::SlipWeakeningFH::p_coefS = 0;
const int pylith::friction::SlipWeakeningFH::p_coefD = 
  pylith::friction::SlipWeakeningFH::p_coefS + 1;
const int pylith::friction::SlipWeakeningFH::p_d0 = 
  pylith::friction::SlipWeakeningFH::p_coefD + 1;
const int pylith::friction::SlipWeakeningFH::p_cohesion =
  pylith::friction::SlipWeakeningFH::p_d0 + 1;
const int pylith::friction::SlipWeakeningFH::p_fw =
  pylith::friction::SlipWeakeningFH::p_cohesion + 1;
const int pylith::friction::SlipWeakeningFH::p_Vw =
  pylith::friction::SlipWeakeningFH::p_fw + 1;

// Indices of database values (order must match dbProperties)
const int pylith::friction::SlipWeakeningFH::db_coefS = 0;
const int pylith::friction::SlipWeakeningFH::db_coefD = 
  pylith::friction::SlipWeakeningFH::db_coefS + 1;
const int pylith::friction::SlipWeakeningFH::db_d0 = 
  pylith::friction::SlipWeakeningFH::db_coefD + 1;
const int pylith::friction::SlipWeakeningFH::db_cohesion =
  pylith::friction::SlipWeakeningFH::db_d0 + 1;
const int pylith::friction::SlipWeakeningFH::db_fw =
  pylith::friction::SlipWeakeningFH::db_cohesion + 1;
const int pylith::friction::SlipWeakeningFH::db_Vw =
  pylith::friction::SlipWeakeningFH::db_fw + 1;

// Indices of state variables.
const int pylith::friction::SlipWeakeningFH::s_slipCum = 0;
const int pylith::friction::SlipWeakeningFH::s_slipPrev = 
  pylith::friction::SlipWeakeningFH::s_slipCum + 1;

// Indices of database values (order must match dbProperties)
const int pylith::friction::SlipWeakeningFH::db_slipCum = 0;
const int pylith::friction::SlipWeakeningFH::db_slipPrev = 
  pylith::friction::SlipWeakeningFH::db_slipCum + 1;

// ----------------------------------------------------------------------
// Default constructor.
pylith::friction::SlipWeakeningFH::SlipWeakeningFH(void) :
  FrictionModel(materials::Metadata(_SlipWeakeningFH::properties,
				    _SlipWeakeningFH::numProperties,
				    _SlipWeakeningFH::dbProperties,
				    _SlipWeakeningFH::numDBProperties,
				    _SlipWeakeningFH::stateVars,
				    _SlipWeakeningFH::numStateVars,
				    _SlipWeakeningFH::dbStateVars,
				    _SlipWeakeningFH::numDBStateVars)),
  _forceHealing(false)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::friction::SlipWeakeningFH::~SlipWeakeningFH(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Compute properties from values in spatial database.
void
pylith::friction::SlipWeakeningFH::forceHealing(const bool flag)
{ // forceHealing
  _forceHealing = flag;
} // forceHealing

// ----------------------------------------------------------------------
// Compute properties from values in spatial database.
void
pylith::friction::SlipWeakeningFH::_dbToProperties(
					   PylithScalar* const propValues,
					   const scalar_array& dbValues) const
{ // _dbToProperties
  assert(propValues);
  const int numDBValues = dbValues.size();
  assert(_SlipWeakeningFH::numDBProperties == numDBValues);

  const PylithScalar db_static = dbValues[db_coefS];
  const PylithScalar db_dynamic = dbValues[db_coefD];
  const PylithScalar db_do = dbValues[db_d0];
  const PylithScalar db_c = dbValues[db_cohesion];
  const PylithScalar db_FW = dbValues[db_fw];
  const PylithScalar db_VW = dbValues[db_Vw];

  if (db_static < 0.0) {
    std::ostringstream msg;
    msg << "Spatial database returned negative value for static coefficient "
	<< "of friction.\n"
	<< "static coefficient of friction: " << db_static << "\n";
    throw std::runtime_error(msg.str());
  } // if

  if (db_dynamic < 0.0) {
    std::ostringstream msg;
    msg << "Spatial database returned negative value for dynamic coefficient "
	<< "of friction.\n"
	<< "dynamic coefficient of friction: " << db_dynamic << "\n";
    throw std::runtime_error(msg.str());
  } // if

  if (db_d0 <= 0.0) {
    std::ostringstream msg;
    msg << "Spatial database returned nonpositive value for slip weakening parameter "
	<< "of friction.\n"
	<< "slip weakening parameter of friction: " << db_do << "\n";
    throw std::runtime_error(msg.str());
  } // if

  propValues[p_coefS] = db_static;
  propValues[p_coefD] = db_dynamic;
  propValues[p_d0] = db_do;
  propValues[p_cohesion] = db_c;
  propValues[p_fw] = db_FW;
  propValues[p_Vw] = db_VW;
} // _dbToProperties

// ----------------------------------------------------------------------
// Nondimensionalize properties.
void
pylith::friction::SlipWeakeningFH::_nondimProperties(PylithScalar* const values,
						    const int nvalues) const
{ // _nondimProperties
  assert(_normalizer);
  assert(values);
  assert(nvalues == _SlipWeakeningFH::numProperties);

  const PylithScalar lengthScale = _normalizer->lengthScale();
  const PylithScalar pressureScale = _normalizer->pressureScale();
  const PylithScalar timeScale = _normalizer->timeScale();
  values[p_d0] /= lengthScale;
  values[p_cohesion] /= pressureScale;
  values[p_Vw] /= (lengthScale / )
} // _nondimProperties

// ----------------------------------------------------------------------
// Dimensionalize properties.
void
pylith::friction::SlipWeakeningFH::_dimProperties(PylithScalar* const values,
						      const int nvalues) const
{ // _dimProperties
  assert(_normalizer);
  assert(values);
  assert(nvalues == _SlipWeakeningFH::numProperties);

  const PylithScalar lengthScale = _normalizer->lengthScale();
  const PylithScalar pressureScale = _normalizer->pressureScale();

  values[p_d0] *= lengthScale;
  values[p_cohesion] *= pressureScale;
} // _dimProperties

// ----------------------------------------------------------------------
// Compute state variables from values in spatial database.
void
pylith::friction::SlipWeakeningFH::_dbToStateVars(
					   PylithScalar* const stateValues,
					   const scalar_array& dbValues) const
{ // _dbToStateVars
  assert(stateValues);
  const int numDBValues = dbValues.size();
  assert(_SlipWeakeningFH::numDBStateVars == numDBValues);

  const PylithScalar cumulativeSlip = dbValues[db_slipCum];
  const PylithScalar previousSlip = dbValues[db_slipPrev];
 
  stateValues[s_slipCum] = cumulativeSlip;
  stateValues[s_slipPrev] = previousSlip;
} // _dbToStateVars

// ----------------------------------------------------------------------
// Nondimensionalize state variables.
void
pylith::friction::SlipWeakeningFH::_nondimStateVars(PylithScalar* const values,
						  const int nvalues) const
{ // _nondimStateVars
  assert(_normalizer);
  assert(values);
  assert(nvalues == _SlipWeakeningFH::numStateVars);

  const PylithScalar lengthScale = _normalizer->lengthScale();

  values[s_slipCum] /= lengthScale;
  values[s_slipPrev] /= lengthScale;
} // _nondimStateVars

// ----------------------------------------------------------------------
// Dimensionalize state variables.
void
pylith::friction::SlipWeakeningFH::_dimStateVars(PylithScalar* const values,
					       const int nvalues) const
{ // _dimStateVars
  assert(_normalizer);
  assert(values);
  assert(nvalues == _SlipWeakeningFH::numStateVars);

  const PylithScalar lengthScale = _normalizer->lengthScale();

  values[s_slipCum] *= lengthScale;
  values[s_slipPrev] *= lengthScale;
} // _dimStateVars

// ----------------------------------------------------------------------
// Compute friction from properties and state variables.
PylithScalar
pylith::friction::SlipWeakeningFH::_calcFriction(const PylithScalar t,
					       const PylithScalar slip,
					       const PylithScalar slipRate,
					       const PylithScalar normalTraction,
					       const PylithScalar* properties,
					       const int numProperties,
					       const PylithScalar* stateVars,
					       const int numStateVars)
{ // _calcFriction
  assert(properties);
  assert(_SlipWeakeningFH::numProperties == numProperties);
  assert(stateVars);
  assert(_SlipWeakeningFH::numStateVars == numStateVars);

  PylithScalar friction = 0.0;
  PylithScalar mu_f = 0.0;
  if (normalTraction <= 0.0) {
    // if fault is in compression
    const PylithScalar slipPrev = stateVars[s_slipPrev];
    const PylithScalar slipCum = stateVars[s_slipCum] + fabs(slip - slipPrev);

    if (slipCum < properties[p_d0]) {
	// if/else linear slip-weakening form of mu_f 
	mu_f = properties[p_coefS] -
	  (properties[p_coefS] - properties[p_coefD]) * 
	  slipCum / properties[p_d0];
      } else {
	mu_f = properties[p_coefD];
      } // if/else
    friction = -mu_f * normalTraction + properties[p_cohesion];
  } else { // else
    friction = properties[p_cohesion];
  } // if/else

  PetscLogFlops(10);

  return friction;
} // _calcFriction

// ----------------------------------------------------------------------
// Compute derivative of friction with slip from properties and
// state variables.
PylithScalar
pylith::friction::SlipWeakeningFH::_calcFrictionDeriv(const PylithScalar t,
						    const PylithScalar slip,
						    const PylithScalar slipRate,
						    const PylithScalar normalTraction,
						    const PylithScalar* properties,
						    const int numProperties,
						    const PylithScalar* stateVars,
						    const int numStateVars)
{ // _calcFrictionDeriv
  assert(properties);
  assert(_SlipWeakeningFH::numProperties == numProperties);
  assert(stateVars);
  assert(_SlipWeakeningFH::numStateVars == numStateVars);

  PylithScalar frictionDeriv = 0.0;
  if (normalTraction <= 0.0) {
    // if fault is in compression
    const PylithScalar slipPrev = stateVars[s_slipPrev];
    const PylithScalar slipCum = stateVars[s_slipCum] + fabs(slip - slipPrev);

    if (slipCum < properties[p_d0]) {
      frictionDeriv = normalTraction * (properties[p_coefS] - properties[p_coefD]) / properties[p_d0];
    } // if
  } // if

  PetscLogFlops(6);

  return frictionDeriv;
} // _calcFrictionDeriv


// ----------------------------------------------------------------------
// Update state variables (for next time step).
void
pylith::friction::SlipWeakeningFH::_updateStateVars(const PylithScalar t,
						  const PylithScalar slip,
						  const PylithScalar slipRate,
						  const PylithScalar normalTraction,
						  PylithScalar* const stateVars,
						  const int numStateVars,
						  const PylithScalar* properties,
						  const int numProperties)
{ // _updateStateVars
  assert(properties);
  assert(_SlipWeakeningFH::numProperties == numProperties);
  assert(stateVars);
  assert(_SlipWeakeningFH::numStateVars == numStateVars);

  const PylithScalar tolerance = 1.0e-12;
  if (slipRate > tolerance && !_forceHealing) {
    const PylithScalar slipPrev = stateVars[s_slipPrev];

    stateVars[s_slipPrev] = slip;
    stateVars[s_slipCum] += fabs(slip - slipPrev);
  } else {
    // Sliding has stopped, so reset state variables.
    stateVars[s_slipPrev] = slip;
    stateVars[s_slipCum] = 0.0;
  } // else

  PetscLogFlops(3);
} // _updateStateVars


// End of file 
