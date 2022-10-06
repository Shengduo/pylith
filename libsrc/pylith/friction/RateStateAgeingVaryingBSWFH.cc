// -*- C++ -*-
//
// ----------------------------------------------------------------------
// Shengduo Liu, Caltech
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

#include "RateStateAgeingVaryingBSWFH.hh" // implementation of object methods

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
    namespace _RateStateAgeingVaryingBSWFH {

      // Number of physical properties.
      const int numProperties = 13;

      // Physical properties.
      const pylith::materials::Metadata::ParamDescription properties[] = {
        { "reference_friction_coefficient", 1, pylith::topology::FieldBase::SCALAR },
        { "reference_slip_rate", 1, pylith::topology::FieldBase::SCALAR },
        { "characteristic_slip_distance", 1, pylith::topology::FieldBase::SCALAR },
        { "constitutive_parameter_a", 1, pylith::topology::FieldBase::SCALAR },
        { "constitutive_parameter_b", 1, pylith::topology::FieldBase::SCALAR },
	      { "cohesion", 1, pylith::topology::FieldBase::SCALAR },
        { "flash_heating_coefficient", 1, pylith::topology::FieldBase::SCALAR }, 
        { "flash_heating_slip_rate", 1, pylith::topology::FieldBase::SCALAR }, 
        { "constitutive_parameter_l_low", 1, pylith::topology::FieldBase::SCALAR }, 
        { "constitutive_parameter_l_high", 1, pylith::topology::FieldBase::SCALAR }, 
        { "constitutive_parameter_si_low", 1, pylith::topology::FieldBase::SCALAR }, 
        { "constitutive_parameter_si_high", 1, pylith::topology::FieldBase::SCALAR }, 
        { "constitutive_parameter_fwLexp", 1, pylith::topology::FieldBase::SCALAR }
      };

      // Number of State Variables.
      const int numStateVars = 1;

      // State Variables.
      const pylith::materials::Metadata::ParamDescription stateVars[] = {
        { "state_variable", 1, pylith::topology::FieldBase::SCALAR },
      };

      // Values expected in spatial database
      const int numDBProperties = 13;
      const char* dbProperties[13] = {
	"reference-friction-coefficient",
	"reference-slip-rate",
	"characteristic-slip-distance",
	"constitutive-parameter-a",
	"constitutive-parameter-b",
	"cohesion",
  "flash-heating-coefficient", 
  "flash-heating-slip-rate", 
  "constitutive-parameter-l-low", 
  "constitutive-parameter-l-high", 
  "constitutive-parameter-si-low", 
  "constitutive-parameter-si-high", 
  "constitutive-parameter-fwLexp"
      };

      const int numDBStateVars = 1;
      const char* dbStateVars[1] = {
	"state-variable",
      };
      
    } // _RateStateAgeingVaryingBSWFH
  } // friction
} // pylith

// Indices of physical properties
const int pylith::friction::RateStateAgeingVaryingBSWFH::p_coef = 0;
const int pylith::friction::RateStateAgeingVaryingBSWFH::p_slipRate0 = 
  pylith::friction::RateStateAgeingVaryingBSWFH::p_coef + 1;
const int pylith::friction::RateStateAgeingVaryingBSWFH::p_L = 
  pylith::friction::RateStateAgeingVaryingBSWFH::p_slipRate0 + 1;
const int pylith::friction::RateStateAgeingVaryingBSWFH::p_a = 
  pylith::friction::RateStateAgeingVaryingBSWFH::p_L + 1;
const int pylith::friction::RateStateAgeingVaryingBSWFH::p_b = 
  pylith::friction::RateStateAgeingVaryingBSWFH::p_a + 1;
const int pylith::friction::RateStateAgeingVaryingBSWFH::p_cohesion =
  pylith::friction::RateStateAgeingVaryingBSWFH::p_b + 1;
const int pylith::friction::RateStateAgeingVaryingBSWFH::p_fwcoef =
  pylith::friction::RateStateAgeingVaryingBSWFH::p_cohesion + 1;
const int pylith::friction::RateStateAgeingVaryingBSWFH::p_fwSlipRate =
  pylith::friction::RateStateAgeingVaryingBSWFH::p_fwcoef + 1;
const int pylith::friction::RateStateAgeingVaryingBSWFH::p_fwLLow =
  pylith::friction::RateStateAgeingVaryingBSWFH::p_fwSlipRate + 1;
const int pylith::friction::RateStateAgeingVaryingBSWFH::p_fwLHigh =
  pylith::friction::RateStateAgeingVaryingBSWFH::p_fwLLow + 1;
const int pylith::friction::RateStateAgeingVaryingBSWFH::p_fwSiLow =
  pylith::friction::RateStateAgeingVaryingBSWFH::p_fwLHigh + 1;
const int pylith::friction::RateStateAgeingVaryingBSWFH::p_fwSiHigh =
  pylith::friction::RateStateAgeingVaryingBSWFH::p_fwSiLow + 1;
const int pylith::friction::RateStateAgeingVaryingBSWFH::p_fwLexp =
  pylith::friction::RateStateAgeingVaryingBSWFH::p_fwSiHigh + 1;

// Indices of database values (order must match dbProperties)
const int pylith::friction::RateStateAgeingVaryingBSWFH::db_coef = 0;
const int pylith::friction::RateStateAgeingVaryingBSWFH::db_slipRate0 = 
  pylith::friction::RateStateAgeingVaryingBSWFH::db_coef + 1;
const int pylith::friction::RateStateAgeingVaryingBSWFH::db_L = 
  pylith::friction::RateStateAgeingVaryingBSWFH::db_slipRate0 + 1;
const int pylith::friction::RateStateAgeingVaryingBSWFH::db_a = 
  pylith::friction::RateStateAgeingVaryingBSWFH::db_L + 1;
const int pylith::friction::RateStateAgeingVaryingBSWFH::db_b = 
  pylith::friction::RateStateAgeingVaryingBSWFH::db_a + 1;
const int pylith::friction::RateStateAgeingVaryingBSWFH::db_cohesion =
  pylith::friction::RateStateAgeingVaryingBSWFH::db_b + 1;
const int pylith::friction::RateStateAgeingVaryingBSWFH::db_fwcoef =
  pylith::friction::RateStateAgeingVaryingBSWFH::db_cohesion + 1;
const int pylith::friction::RateStateAgeingVaryingBSWFH::db_fwSlipRate =
  pylith::friction::RateStateAgeingVaryingBSWFH::db_fwcoef + 1;
const int pylith::friction::RateStateAgeingVaryingBSWFH::db_fwLLow =
  pylith::friction::RateStateAgeingVaryingBSWFH::db_fwSlipRate + 1;
const int pylith::friction::RateStateAgeingVaryingBSWFH::db_fwLHigh =
  pylith::friction::RateStateAgeingVaryingBSWFH::db_fwLLow + 1;
const int pylith::friction::RateStateAgeingVaryingBSWFH::db_fwSiLow =
  pylith::friction::RateStateAgeingVaryingBSWFH::db_fwLHigh + 1;
const int pylith::friction::RateStateAgeingVaryingBSWFH::db_fwSiHigh =
  pylith::friction::RateStateAgeingVaryingBSWFH::db_fwSiLow + 1;
const int pylith::friction::RateStateAgeingVaryingBSWFH::db_fwLexp =
  pylith::friction::RateStateAgeingVaryingBSWFH::db_fwSiHigh + 1;

// Indices of state variables.
const int pylith::friction::RateStateAgeingVaryingBSWFH::s_state = 0;

// Indices of database values (order must match dbProperties)
const int pylith::friction::RateStateAgeingVaryingBSWFH::db_state = 0;

// ----------------------------------------------------------------------
// Default constructor.
pylith::friction::RateStateAgeingVaryingBSWFH::RateStateAgeingVaryingBSWFH(void) :
  FrictionModel(materials::Metadata(_RateStateAgeingVaryingBSWFH::properties,
				    _RateStateAgeingVaryingBSWFH::numProperties,
				    _RateStateAgeingVaryingBSWFH::dbProperties,
				    _RateStateAgeingVaryingBSWFH::numDBProperties,
				    _RateStateAgeingVaryingBSWFH::stateVars,
				    _RateStateAgeingVaryingBSWFH::numStateVars,
				    _RateStateAgeingVaryingBSWFH::dbStateVars,
				    _RateStateAgeingVaryingBSWFH::numDBStateVars)),
  _linearSlipRate(1.0e-12)
{ // constructor
  _normalTraction = -10.0e6;
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::friction::RateStateAgeingVaryingBSWFH::~RateStateAgeingVaryingBSWFH(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Set nondimensional slip rate below which friction varies
//  linearly with slip rate.
void
pylith::friction::RateStateAgeingVaryingBSWFH::linearSlipRate(const PylithScalar value)
{ // linearSlipRate
  if (value < 0.0) {
    std::ostringstream msg;
    msg << "Nondimensional linear slip rate threshold (" << value << ") for rate state friction model with flash heating and varying ratestate-b"
	<< label() << " must be nonnegative.";
    throw std::runtime_error(msg.str());
  } // if
  
  _linearSlipRate = value;
} // linearSlipRate

// ----------------------------------------------------------------------
// Compute properties from values in spatial database.
void
pylith::friction::RateStateAgeingVaryingBSWFH::_dbToProperties(PylithScalar* const propValues,
						   const scalar_array& dbValues) const
{ // _dbToProperties
  assert(propValues);
  const int numDBValues = dbValues.size();
  assert(_RateStateAgeingVaryingBSWFH::numDBProperties == numDBValues);

  const PylithScalar frictionCoef = dbValues[db_coef];
  const PylithScalar slipRate0 = dbValues[db_slipRate0];
  const PylithScalar dc = dbValues[db_L];
  const PylithScalar a = dbValues[db_a];
  const PylithScalar b = dbValues[db_b];
  const PylithScalar cohesion = dbValues[db_cohesion];
  const PylithScalar FHfrictionCoef = dbValues[db_fwcoef];
  const PylithScalar FHSlipRate = dbValues[db_fwSlipRate];
  const PylithScalar FHLLow = dbValues[db_fwLLow];
  const PylithScalar FHLHigh = dbValues[db_fwLHigh];
  const PylithScalar FHSiLow = dbValues[db_fwSiLow];
  const PylithScalar FHSiHigh = dbValues[db_fwSiHigh];
  const PylithScalar FHLexp = dbValues[db_fwLexp];

  if (frictionCoef < 0.0) {
    std::ostringstream msg;
    msg << "Spatial database returned negative value for reference coefficient "
	<< "of Rate and State friction Ageing Law.\n"
	<< "reference coefficient of friction: " << frictionCoef << "\n";
    throw std::runtime_error(msg.str());
  } // if

  if (dc <= 0.0) {
    std::ostringstream msg;
    msg << "Spatial database returned nonpositive value for characteristic"
	<< "slip distance of Rate and State friction Ageing Law.\n"
	<< "characteristic slip distance: " << dc << "\n";
    throw std::runtime_error(msg.str());
  } // if

  if (a <= 0.0) {
    std::ostringstream msg;
    msg << "Spatial database returned nonpositive value for constitutive "
	<< "parameter 'a' of Rate and State friction Ageing Law.\n"
	<< "Rate and State parameter 'a' of Ageing Law of friction: " << a << "\n";
    throw std::runtime_error(msg.str());
  } // if

  if (b <= 0.0) {
    std::ostringstream msg;
    msg << "Spatial database returned nonpositive value for constitutive "
	<< "parameter 'b' of Rate and State friction Ageing Law.\n"
	<< "Rate and State parameter 'b' of Ageing Law of friction: " << b << "\n";
    throw std::runtime_error(msg.str());
  } // if

  if (FHfrictionCoef <= 0.0) {
    std::ostringstream msg;
    msg << "Spatial database returned nonpositive value for constitutive "
  << "parameter 'flash heating friction coefficient' of Rate and State friction Ageing Law.\n"
  << "Flash heating friction coefficient 'fw' of Ageing Law of friction: " << FHfrictionCoef << "\n";
    throw std::runtime_error(msg.str());
  } // if

  if (FHLLow <= 0.0) {
    std::ostringstream msg;
    msg << "Spatial database returned nonpositive value for constitutive "
  << "parameter 'flash heating slip weakening Llow' of Rate and State friction Ageing Law.\n"
  << "Flash heating slip weakening 'Llow' of Ageing Law of friction: " << FHLLow << "\n";
    throw std::runtime_error(msg.str());
  } // if

  if (FHLHigh <= 0.0) {
    std::ostringstream msg;
    msg << "Spatial database returned nonpositive value for constitutive "
  << "parameter 'flash heating slip weakening Lhigh of Rate and State friction Ageing Law.\n"
  << "Flash heating slip weakening 'Lhigh' of Ageing Law of friction: " << FHLHigh << "\n";
    throw std::runtime_error(msg.str());
  } // if

  if (FHSiLow <= 0.0) {
    std::ostringstream msg;
    msg << "Spatial database returned nonpositive value for constitutive "
  << "parameter 'flash heating slip weakening silow' of Rate and State friction Ageing Law.\n"
  << "Flash heating slip weakening 'silow' of Ageing Law of friction: " << FHSiLow << "\n";
    throw std::runtime_error(msg.str());
  } // if

  if (FHSiHigh <= 0.0) {
    std::ostringstream msg;
    msg << "Spatial database returned nonpositive value for constitutive "
  << "parameter 'flash heating slip weakening sihigh' of Rate and State friction Ageing Law.\n"
  << "Flash heating slip weakening 'sihigh' of Ageing Law of friction: " << FHSiHigh << "\n";
    throw std::runtime_error(msg.str());
  } // if

  if (FHSlipRate <= 0.0) {
    std::ostringstream msg;
    msg << "Spatial database returned nonpositive value for constitutive "
  << "parameter 'flash heating slip rate of Rate and State friction Ageing Law.\n"
  << "Flash heating friction slip rate 'Vw' of Ageing Law of friction: " << FHSlipRate << "\n";
    throw std::runtime_error(msg.str());
  } // if

  if (FHLexp <= 0.0) {
    std::ostringstream msg;
    msg << "Spatial database returned nonpositive value for constitutive "
  << "parameter 'flash heating exponential decay distance fwLexp.\n"
  << "Flash heating decay distance 'fwLexp' of Ageing Law of friction: " << FHLexp << "\n";
    throw std::runtime_error(msg.str());
  } // if

  propValues[p_coef] = frictionCoef;
  propValues[p_slipRate0] = slipRate0;
  propValues[p_L] = dc;
  propValues[p_a] = a;
  propValues[p_b] = b;
  propValues[p_cohesion] = cohesion;
  propValues[p_fwcoef] = FHfrictionCoef;
  propValues[p_fwSlipRate] = FHSlipRate;
  propValues[p_fwLLow] = FHLLow;
  propValues[p_fwLHigh] = FHLHigh;
  propValues[p_fwSiLow] = FHSiLow;
  propValues[p_fwSiHigh] = FHSiHigh;
  propValues[p_fwLexp] = FHLexp;

} // _dbToProperties

// ----------------------------------------------------------------------
// Nondimensionalize properties.
void
pylith::friction::RateStateAgeingVaryingBSWFH::_nondimProperties(PylithScalar* const values,
						    const int nvalues) const
{ // _nondimProperties
  assert(_normalizer);
  assert(values);
  assert(nvalues == _RateStateAgeingVaryingBSWFH::numProperties);

  const PylithScalar lengthScale = _normalizer->lengthScale();
  const PylithScalar timeScale = _normalizer->timeScale();
  const PylithScalar pressureScale = _normalizer->pressureScale();

  values[p_slipRate0] /= lengthScale / timeScale;
  values[p_L] /= lengthScale;
  values[p_cohesion] /= pressureScale;
  values[p_fwSlipRate] /= lengthScale / timeScale;
  values[p_fwLLow] /= lengthScale;
  values[p_fwLHigh] /= lengthScale;
  values[p_fwSiLow] /= pressureScale;
  values[p_fwSiHigh] /= pressureScale;
  values[p_fwLexp] /= lengthScale;

} // _nondimProperties

// ----------------------------------------------------------------------
// Dimensionalize properties.
void
pylith::friction::RateStateAgeingVaryingBSWFH::_dimProperties(PylithScalar* const values,
						  const int nvalues) const
{ // _dimProperties
  assert(_normalizer);
  assert(values);
  assert(nvalues == _RateStateAgeingVaryingBSWFH::numProperties);

  const PylithScalar lengthScale = _normalizer->lengthScale();
  const PylithScalar timeScale = _normalizer->timeScale();
  const PylithScalar pressureScale = _normalizer->pressureScale();

  values[p_slipRate0] *= lengthScale / timeScale;
  values[p_L] *= lengthScale;
  values[p_cohesion] *= pressureScale;
  values[p_fwSlipRate] *= lengthScale / timeScale;
  values[p_fwLLow] *= lengthScale;
  values[p_fwLHigh] *= lengthScale;
  values[p_fwSiLow] *= pressureScale;
  values[p_fwSiHigh] *= pressureScale;
  values[p_fwLexp] *= lengthScale;
} // _dimProperties

// ----------------------------------------------------------------------
// Compute state variables from values in spatial database.
void
pylith::friction::RateStateAgeingVaryingBSWFH::_dbToStateVars(PylithScalar* const stateValues,
						  const scalar_array& dbValues) const
{ // _dbToStateVars
  assert(stateValues);
  const int numDBValues = dbValues.size();
  assert(_RateStateAgeingVaryingBSWFH::numDBStateVars == numDBValues);

  stateValues[s_state] = dbValues[db_state];
} // _dbToStateVars

// ----------------------------------------------------------------------
// Nondimensionalize state variables.
void
pylith::friction::RateStateAgeingVaryingBSWFH::_nondimStateVars(PylithScalar* const values,
						    const int nvalues) const
{ // _nondimStateVars
  assert(_normalizer);
  assert(values);
  assert(nvalues == _RateStateAgeingVaryingBSWFH::numStateVars);

  const PylithScalar timeScale = _normalizer->timeScale();

  values[s_state] /= timeScale;
} // _nondimStateVars

// ----------------------------------------------------------------------
// Dimensionalize state variables.
void
pylith::friction::RateStateAgeingVaryingBSWFH::_dimStateVars(PylithScalar* const values,
						 const int nvalues) const
{ // _dimStateVars
  assert(_normalizer);
  assert(values);
  assert(nvalues == _RateStateAgeingVaryingBSWFH::numStateVars);

  const PylithScalar timeScale = _normalizer->timeScale();

  values[s_state] *= timeScale;
} // _dimStateVars

// ----------------------------------------------------------------------
// Compute friction from properties and state variables.
PylithScalar
pylith::friction::RateStateAgeingVaryingBSWFH::_calcFriction(const PylithScalar t,
						 const PylithScalar slip,
						 const PylithScalar slipRate,
						 const PylithScalar normalTraction,
						 const PylithScalar* properties,
						 const int numProperties,
						 const PylithScalar* stateVars,
						 const int numStateVars)
{ // _calcFriction
  assert(properties);
  assert(_RateStateAgeingVaryingBSWFH::numProperties == numProperties);
  assert(numStateVars);
  assert(_RateStateAgeingVaryingBSWFH::numStateVars == numStateVars);

  PylithScalar friction = 0.0;
  PylithScalar mu_f = 0.0;
  if (normalTraction <= 0.0) {
    // if fault is in compression
    if(t >= 0. && _normalTraction < 0.) _normalTraction = -normalTraction;

    const PylithScalar slipRateLinear = _linearSlipRate;

    const PylithScalar f0= properties[p_coef];
    const PylithScalar a = properties[p_a];
    const PylithScalar b = properties[p_b];
    const PylithScalar L = properties[p_L];
    const PylithScalar slipRate0 = properties[p_slipRate0];
    const PylithScalar fw = properties[p_fwcoef];
    const PylithScalar fwSlipRate = properties[p_fwSlipRate];
    const PylithScalar fwLLow = properties[p_fwLLow];
    const PylithScalar fwLHigh = properties[p_fwLHigh];
    const PylithScalar fwSiLow = properties[p_fwSiLow];
    const PylithScalar fwSiHigh = properties[p_fwSiHigh];
    const PylithScalar fwLexp = properties[p_fwLexp];

    // Prevent zero value for theta, reasonable value is L / slipRate0
    const PylithScalar si0 = 10.0e6;
    const PylithScalar theta = (stateVars[s_state] > 0.0) ? stateVars[s_state] : L / slipRate0;
    const PylithScalar beta = log(fwLLow / fwLHigh) / log(fwSiLow / fwSiHigh);
    const PylithScalar alp = fwLLow / pow(fwSiLow / si0, beta);
    const PylithScalar Lw = t >= 0. ? alp * pow(_normalTraction / si0, beta) : alp;

    // const PylithScalar beta = (fwLLow - fwLHigh) / (fwSiLow - fwSiHigh);
    // const PylithScalar alp = fwLLow - beta * fwSiLow;
    // PylithScalar Lw = beta * (-normalTraction) + alp;
    // if (-normalTraction < fwSiLow) Lw = fwLLow;
    // else if (-normalTraction > fwSiHigh) Lw = fwLHigh;

    if (slipRate >= slipRateLinear) {
      mu_f = f0 + a*log(slipRate / slipRate0) + b*log(slipRate0*theta/L);
      mu_f = fw + (mu_f - fw) / (1. + (L / theta / fwSlipRate + slip / Lw) * exp(-slip / fwLexp));
    } else {
      mu_f = f0 + a*log(slipRateLinear / slipRate0) + b*log(slipRate0*theta/L) -
	a*(1.0 - slipRate/slipRateLinear);
      mu_f = fw + (mu_f - fw) / (1. + (L / theta / fwSlipRate + slip / Lw) * exp(-slip / fwLexp));
    } // else
    friction = -mu_f * normalTraction + properties[p_cohesion];
    
  } else {
    friction = properties[p_cohesion];
  } // if/else

  PetscLogFlops(12);

  return friction;
} // _calcFriction


// ----------------------------------------------------------------------
// Compute derivative of friction with slip from properties and
// state variables.
PylithScalar
pylith::friction::RateStateAgeingVaryingBSWFH::_calcFrictionDeriv(const PylithScalar t,
						      const PylithScalar slip,
						      const PylithScalar slipRate,
						      const PylithScalar normalTraction,
						      const PylithScalar* properties,
						      const int numProperties,
						      const PylithScalar* stateVars,
						      const int numStateVars)
{ // _calcFrictionDeriv
  assert(properties);
  assert(_RateStateAgeingVaryingBSWFH::numProperties == numProperties);
  assert(numStateVars);
  assert(_RateStateAgeingVaryingBSWFH::numStateVars == numStateVars);

  PylithScalar frictionDeriv = 0.0;
  if (normalTraction <= 0.0) {
    // if fault is in compression

    // initialize _normaltraction
    if(t >= 0. && _normalTraction < 0.) _normalTraction = -normalTraction;

    const PylithScalar slipRateLinear = _linearSlipRate;
    const PylithScalar f0= properties[p_coef];
    const PylithScalar a = properties[p_a];
    const PylithScalar b = properties[p_b];
    const PylithScalar L = properties[p_L];
    const PylithScalar fw = properties[p_fwcoef];
    const PylithScalar Vw = properties[p_fwSlipRate];
    const PylithScalar slipRate0 = properties[p_slipRate0];
    const PylithScalar fwLLow = properties[p_fwLLow];
    const PylithScalar fwLHigh = properties[p_fwLHigh];
    const PylithScalar fwSiLow = properties[p_fwSiLow];
    const PylithScalar fwSiHigh = properties[p_fwSiHigh];
    const PylithScalar fwLexp = properties[p_fwLexp];

    // Prevent zero value for theta, reasonable value is L / slipRate0
    // lw = alpha * (si / si0)^beta
    const PylithScalar si0 = 10.0e6;
    const PylithScalar theta = (stateVars[s_state] > 0.0) ? stateVars[s_state] : L / slipRate0;
    const PylithScalar beta = log(fwLLow / fwLHigh) / log(fwSiLow / fwSiHigh);
    const PylithScalar alp = fwLLow / pow(fwSiLow / si0, beta);
    const PylithScalar Lw = t >= 0.? alp * pow(_normalTraction / si0, beta) : alp;

    // const PylithScalar beta = (fwLLow - fwLHigh) / (fwSiLow - fwSiHigh);
    // const PylithScalar alp = fwLLow - beta * fwSiLow;
    // PylithScalar Lw = beta * (-normalTraction) + alp;
    // if (-normalTraction < fwSiLow) Lw = fwLLow;
    // else if (-normalTraction > fwSiHigh) Lw = fwLHigh;
    
    PylithScalar fRS;

    if (slipRate >= slipRateLinear) {
      // frictionDeriv = -normalTraction * a / (slipRate * _dt);
      fRS = f0 + a * log(slipRate / slipRate0) + b * log(slipRate0 * theta / L);
      frictionDeriv = -normalTraction * ((a / (slipRate * _dt)) / (1. + (L / theta / Vw + slip / Lw) * exp(-slip / fwLexp))
                                          - (fRS - fw) / (1. + (L / theta / Vw + slip / Lw) * exp(-slip / fwLexp)) / (1. + (L / theta / Vw + slip / Lw) * exp(-slip / fwLexp))
                                            * ((1 - slip / fwLexp) / Lw - L / theta / Vw / fwLexp) * exp(-slip / fwLexp));
    } else {
      fRS = f0 + a*log(slipRateLinear / slipRate0) + b*log(slipRate0*theta/L) - a*(1.0 - slipRate/slipRateLinear);
      frictionDeriv = -normalTraction * ((a / (slipRateLinear * _dt)) / (1. + L / theta / Vw + slip * exp(-slip / fwLexp) / Lw)
                                          - (fRS - fw) / (1. + (L / theta / Vw + slip / Lw) * exp(-slip / fwLexp)) / (1. + (L / theta / Vw + slip / Lw) * exp(-slip / fwLexp))
                                            * ((1 - slip / fwLexp) / Lw - L / theta / Vw / fwLexp) * exp(-slip / fwLexp));
    }
  } // if    

  PetscLogFlops(12);

  return frictionDeriv;
} // _calcFrictionDeriv


// ----------------------------------------------------------------------
// Update state variables (for next time step).
void
pylith::friction::RateStateAgeingVaryingBSWFH::_updateStateVars(const PylithScalar t,
						    const PylithScalar slip,
						    const PylithScalar slipRate,
						    const PylithScalar normalTraction,
						    PylithScalar* const stateVars,
						    const int numStateVars,
						    const PylithScalar* properties,
						    const int numProperties)
{ // _updateStateVars
  assert(properties);
  assert(_RateStateAgeingVaryingBSWFH::numProperties == numProperties);
  assert(numStateVars);
  assert(_RateStateAgeingVaryingBSWFH::numStateVars == numStateVars);

  // d(theta)/dt = (1 - slipRate * theta / L)
  //
  // Use separation of variables to integrate above ODE from t->t+dt,
  // keeping slip rate constant.
  //
  // thetaTpdt = thetaT * exp(-slipRate/L * dt)
  //             + L/slipRate * (1 -  exp(-slipRate/L * dt))
  //
  // As slipRate --> 0, L/sliprate --> infinity and
  // exp(-sliprate/L*dt) --> 1.  To determine, d(theta)/dt near
  // sliprate == 0, we expand the exponential term in a Taylor series:
  //
  // exp(-x) = 1 - x +1/2*x**2 + 1/6*x**3
  //
  // This leads to (in the vicinity of slipRate == 0):
  //
  // thetaTpdt = thetaT * exp(-slipRate/L * dt)
  //             + dt - 0.5*(sliprate/L)*dt**2 + 1.0/6.0*(slipRate/L)*dt**3;

  const PylithScalar dt = _dt;
  const PylithScalar thetaTVertex = stateVars[s_state];
  const PylithScalar L = properties[p_L];
  const PylithScalar vDtL = slipRate * dt / L;
  const PylithScalar expTerm = exp(-vDtL);

  PylithScalar thetaTpdtVertex = 0.0;
  if (vDtL > 1.0e-20) {
    thetaTpdtVertex = thetaTVertex * expTerm + L / slipRate * (1 - expTerm);
    PetscLogFlops(7);
  } else {
    thetaTpdtVertex = thetaTVertex * expTerm + dt - 0.5 * slipRate/L * dt*dt;
    PetscLogFlops(9);
  } // if/else
  
  stateVars[s_state] = thetaTpdtVertex;

} // _updateStateVars


// End of file 
