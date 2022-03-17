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

#include "RateStateAgeingVaryingBFH.hh" // implementation of object methods

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
    namespace _RateStateAgeingVaryingBFH {

      // Number of physical properties.
      const int numProperties = 10;

      // Physical properties.
      const pylith::materials::Metadata::ParamDescription properties[] = {
        { "reference_friction_coefficient", 1, pylith::topology::FieldBase::SCALAR },
        { "reference_slip_rate", 1, pylith::topology::FieldBase::SCALAR },
        { "characteristic_slip_distance", 1, pylith::topology::FieldBase::SCALAR },
        { "constitutive_parameter_a", 1, pylith::topology::FieldBase::SCALAR },
        { "constitutive_parameter_b_initial", 1, pylith::topology::FieldBase::SCALAR },
        { "constitutive_parameter_b_final", 1, pylith::topology::FieldBase::SCALAR },
        { "constitutive_parameter_b_distance", 1, pylith::topology::FieldBase::SCALAR },
	      { "cohesion", 1, pylith::topology::FieldBase::SCALAR },
        { "flash_heating_coefficient", 1, pylith::topology::FieldBase::SCALAR }, 
        { "flash_heating_slip_rate", 1, pylith::topology::FieldBase::SCALAR }
      };

      // Number of State Variables.
      const int numStateVars = 1;

      // State Variables.
      const pylith::materials::Metadata::ParamDescription stateVars[] = {
        { "state_variable", 1, pylith::topology::FieldBase::SCALAR },
      };

      // Values expected in spatial database
      const int numDBProperties = 10;
      const char* dbProperties[10] = {
	"reference-friction-coefficient",
	"reference-slip-rate",
	"characteristic-slip-distance",
	"constitutive-parameter-a",
	"constitutive-parameter-b-initial",
  "constitutive-parameter-b-final", 
  "constitutive-parameter-b-distance",
	"cohesion",
  "flash_heating_coefficient", 
  "flash_heating_slip_rate"
      };

      const int numDBStateVars = 1;
      const char* dbStateVars[1] = {
	"state-variable",
      };
      
    } // _RateStateAgeingVaryingBFH
  } // friction
} // pylith

// Indices of physical properties
const int pylith::friction::RateStateAgeingVaryingBFH::p_coef = 0;
const int pylith::friction::RateStateAgeingVaryingBFH::p_slipRate0 = 
  pylith::friction::RateStateAgeingVaryingBFH::p_coef + 1;
const int pylith::friction::RateStateAgeingVaryingBFH::p_L = 
  pylith::friction::RateStateAgeingVaryingBFH::p_slipRate0 + 1;
const int pylith::friction::RateStateAgeingVaryingBFH::p_a = 
  pylith::friction::RateStateAgeingVaryingBFH::p_L + 1;
const int pylith::friction::RateStateAgeingVaryingBFH::p_b_initial = 
  pylith::friction::RateStateAgeingVaryingBFH::p_a + 1;
const int pylith::friction::RateStateAgeingVaryingBFH::p_b_final = 
  pylith::friction::RateStateAgeingVaryingBFH::p_b_initial + 1;
const int pylith::friction::RateStateAgeingVaryingBFH::p_b_distance = 
  pylith::friction::RateStateAgeingVaryingBFH::p_b_final + 1;
const int pylith::friction::RateStateAgeingVaryingBFH::p_cohesion =
  pylith::friction::RateStateAgeingVaryingBFH::p_b_distance + 1;
const int pylith::friction::RateStateAgeingVaryingBFH::p_fwcoef =
  pylith::friction::RateStateAgeingVaryingBFH::p_cohesion + 1;
const int pylith::friction::RateStateAgeingVaryingBFH::p_fwSlipRate =
  pylith::friction::RateStateAgeingVaryingBFH::p_fwcoef + 1;

// Indices of database values (order must match dbProperties)
const int pylith::friction::RateStateAgeingVaryingBFH::db_coef = 0;
const int pylith::friction::RateStateAgeingVaryingBFH::db_slipRate0 = 
  pylith::friction::RateStateAgeingVaryingBFH::db_coef + 1;
const int pylith::friction::RateStateAgeingVaryingBFH::db_L = 
  pylith::friction::RateStateAgeingVaryingBFH::db_slipRate0 + 1;
const int pylith::friction::RateStateAgeingVaryingBFH::db_a = 
  pylith::friction::RateStateAgeingVaryingBFH::db_L + 1;
const int pylith::friction::RateStateAgeingVaryingBFH::db_b_initial = 
  pylith::friction::RateStateAgeingVaryingBFH::db_a + 1;
const int pylith::friction::RateStateAgeingVaryingBFH::db_b_final = 
  pylith::friction::RateStateAgeingVaryingBFH::db_b_initial + 1;
const int pylith::friction::RateStateAgeingVaryingBFH::db_b_distance = 
  pylith::friction::RateStateAgeingVaryingBFH::db_b_final + 1;
const int pylith::friction::RateStateAgeingVaryingBFH::db_cohesion =
  pylith::friction::RateStateAgeingVaryingBFH::db_b_distance + 1;
const int pylith::friction::RateStateAgeingVaryingBFH::db_fwcoef =
  pylith::friction::RateStateAgeingVaryingBFH::db_cohesion + 1;
const int pylith::friction::RateStateAgeingVaryingBFH::db_fwSlipRate =
  pylith::friction::RateStateAgeingVaryingBFH::db_fwcoef + 1;

// Indices of state variables.
const int pylith::friction::RateStateAgeingVaryingBFH::s_state = 0;

// Indices of database values (order must match dbProperties)
const int pylith::friction::RateStateAgeingVaryingBFH::db_state = 0;

// ----------------------------------------------------------------------
// Default constructor.
pylith::friction::RateStateAgeingVaryingBFH::RateStateAgeingVaryingBFH(void) :
  FrictionModel(materials::Metadata(_RateStateAgeingVaryingBFH::properties,
				    _RateStateAgeingVaryingBFH::numProperties,
				    _RateStateAgeingVaryingBFH::dbProperties,
				    _RateStateAgeingVaryingBFH::numDBProperties,
				    _RateStateAgeingVaryingBFH::stateVars,
				    _RateStateAgeingVaryingBFH::numStateVars,
				    _RateStateAgeingVaryingBFH::dbStateVars,
				    _RateStateAgeingVaryingBFH::numDBStateVars)),
  _linearSlipRate(1.0e-12)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::friction::RateStateAgeingVaryingBFH::~RateStateAgeingVaryingBFH(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Set nondimensional slip rate below which friction varies
//  linearly with slip rate.
void
pylith::friction::RateStateAgeingVaryingBFH::linearSlipRate(const PylithScalar value)
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
pylith::friction::RateStateAgeingVaryingBFH::_dbToProperties(PylithScalar* const propValues,
						   const scalar_array& dbValues) const
{ // _dbToProperties
  assert(propValues);
  const int numDBValues = dbValues.size();
  assert(_RateStateAgeingVaryingBFH::numDBProperties == numDBValues);

  const PylithScalar frictionCoef = dbValues[db_coef];
  const PylithScalar slipRate0 = dbValues[db_slipRate0];
  const PylithScalar dc = dbValues[db_L];
  const PylithScalar a = dbValues[db_a];
  const PylithScalar b_initial = dbValues[db_b_initial];
  const PylithScalar b_final = dbValues[db_b_final];
  const PylithScalar b_distance = dbValues[db_b_distance];
  const PylithScalar cohesion = dbValues[db_cohesion];
  const PylithScalar FHfrictionCoef = dbValues[db_fwcoef];
  const PylithScalar FHSlipRate = dbValues[db_fwSlipRate];
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

  if (b_initial <= 0.0) {
    std::ostringstream msg;
    msg << "Spatial database returned nonpositive value for constitutive "
	<< "parameter 'b_initial' of Rate and State friction Ageing Law.\n"
	<< "Rate and State parameter 'b_initial' of Ageing Law of friction: " << b_initial << "\n";
    throw std::runtime_error(msg.str());
  } // if

  if (b_final <= 0.0) {
    std::ostringstream msg;
    msg << "Spatial database returned nonpositive value for constitutive "
	<< "parameter 'b_final' of Rate and State friction Ageing Law.\n"
	<< "Rate and State parameter 'b_final' of Ageing Law of friction: " << b_final << "\n";
    throw std::runtime_error(msg.str());
  } // if

  if (b_initial <= 0.0) {
    std::ostringstream msg;
    msg << "Spatial database returned nonpositive value for constitutive "
	<< "parameter 'b_distance' of Rate and State friction Ageing Law.\n"
	<< "Rate and State parameter 'b_distance' of Ageing Law of friction: " << b_distance << "\n";
    throw std::runtime_error(msg.str());
  } // if

  if (FHfrictionCoef <= 0.0) {
    std::ostringstream msg;
    msg << "Spatial database returned nonpositive value for constitutive "
  << "parameter 'flash heating friction coefficient' of Rate and State friction Ageing Law.\n"
  << "Flash heating friction coefficient 'fw' of Ageing Law of friction: " << FHfrictionCoef << "\n";
    throw std::runtime_error(msg.str());
  } // if

  propValues[p_coef] = frictionCoef;
  propValues[p_slipRate0] = slipRate0;
  propValues[p_L] = dc;
  propValues[p_a] = a;
  propValues[p_b_initial] = b_initial;
  propValues[p_b_final] = b_final;
  propValues[p_b_distance] = b_distance;
  propValues[p_cohesion] = cohesion;
  propValues[p_fwcoef] = FHfrictionCoef;
  propValues[p_fwSlipRate] = FHSlipRate;

} // _dbToProperties

// ----------------------------------------------------------------------
// Nondimensionalize properties.
void
pylith::friction::RateStateAgeingVaryingBFH::_nondimProperties(PylithScalar* const values,
						    const int nvalues) const
{ // _nondimProperties
  assert(_normalizer);
  assert(values);
  assert(nvalues == _RateStateAgeingVaryingBFH::numProperties);

  const PylithScalar lengthScale = _normalizer->lengthScale();
  const PylithScalar timeScale = _normalizer->timeScale();
  const PylithScalar pressureScale = _normalizer->pressureScale();

  values[p_slipRate0] /= lengthScale / timeScale;
  values[p_L] /= lengthScale;
  values[p_cohesion] /= pressureScale;
  values[p_fwSlipRate] /= lengthScale / timeScale;
  values[p_b_distance] /= lengthScale;
} // _nondimProperties

// ----------------------------------------------------------------------
// Dimensionalize properties.
void
pylith::friction::RateStateAgeingVaryingBFH::_dimProperties(PylithScalar* const values,
						  const int nvalues) const
{ // _dimProperties
  assert(_normalizer);
  assert(values);
  assert(nvalues == _RateStateAgeingVaryingBFH::numProperties);

  const PylithScalar lengthScale = _normalizer->lengthScale();
  const PylithScalar timeScale = _normalizer->timeScale();
  const PylithScalar pressureScale = _normalizer->pressureScale();

  values[p_slipRate0] *= lengthScale / timeScale;
  values[p_L] *= lengthScale;
  values[p_cohesion] *= pressureScale;
  values[p_fwSlipRate] *= lengthScale / timeScale;
  values[p_b_distance] *= lengthScale;
} // _dimProperties

// ----------------------------------------------------------------------
// Compute state variables from values in spatial database.
void
pylith::friction::RateStateAgeingVaryingBFH::_dbToStateVars(PylithScalar* const stateValues,
						  const scalar_array& dbValues) const
{ // _dbToStateVars
  assert(stateValues);
  const int numDBValues = dbValues.size();
  assert(_RateStateAgeingVaryingBFH::numDBStateVars == numDBValues);

  stateValues[s_state] = dbValues[db_state];
} // _dbToStateVars

// ----------------------------------------------------------------------
// Nondimensionalize state variables.
void
pylith::friction::RateStateAgeingVaryingBFH::_nondimStateVars(PylithScalar* const values,
						    const int nvalues) const
{ // _nondimStateVars
  assert(_normalizer);
  assert(values);
  assert(nvalues == _RateStateAgeingVaryingBFH::numStateVars);

  const PylithScalar timeScale = _normalizer->timeScale();

  values[s_state] /= timeScale;
} // _nondimStateVars

// ----------------------------------------------------------------------
// Dimensionalize state variables.
void
pylith::friction::RateStateAgeingVaryingBFH::_dimStateVars(PylithScalar* const values,
						 const int nvalues) const
{ // _dimStateVars
  assert(_normalizer);
  assert(values);
  assert(nvalues == _RateStateAgeingVaryingBFH::numStateVars);

  const PylithScalar timeScale = _normalizer->timeScale();

  values[s_state] *= timeScale;
} // _dimStateVars

// ----------------------------------------------------------------------
// Compute friction from properties and state variables.
PylithScalar
pylith::friction::RateStateAgeingVaryingBFH::_calcFriction(const PylithScalar t,
						 const PylithScalar slip,
						 const PylithScalar slipRate,
						 const PylithScalar normalTraction,
						 const PylithScalar* properties,
						 const int numProperties,
						 const PylithScalar* stateVars,
						 const int numStateVars)
{ // _calcFriction
  assert(properties);
  assert(_RateStateAgeingVaryingBFH::numProperties == numProperties);
  assert(numStateVars);
  assert(_RateStateAgeingVaryingBFH::numStateVars == numStateVars);

  PylithScalar friction = 0.0;
  PylithScalar mu_f = 0.0;
  if (normalTraction <= 0.0) {
    // if fault is in compression

    const PylithScalar slipRateLinear = _linearSlipRate;

    const PylithScalar f0 = properties[p_coef];
    const PylithScalar a = properties[p_a];
    const PylithScalar b_initial = properties[p_b_initial];
    const PylithScalar b_final = properties[p_b_final];
    const PylithScalar b_distance = properties[p_b_distance];
    const PylithScalar L = properties[p_L];
    const PylithScalar slipRate0 = properties[p_slipRate0];
    const PylithScalar fw = properties[p_fwcoef];
    const PylithScalar fwSlipRate = properties[p_fwSlipRate];

    // Prevent zero value for theta, reasonable value is L / slipRate0
    const PylithScalar theta = (stateVars[s_state] > 0.0) ? stateVars[s_state] : L / slipRate0;
    PylithScalar b;
    if (slip > b_distance) {
      b = b_final;
    } else {
      b = b_initial + (b_final - b_initial) / b_distance * slip;
    }
    
    if (slipRate >= slipRateLinear) {
      mu_f = f0 + a*log(slipRate / slipRate0) + b*log(slipRate0*theta/L);
      mu_f = fw + (mu_f - fw) / (1. + L / theta / fwSlipRate);
    } else {
      mu_f = f0 + a*log(slipRateLinear / slipRate0) + b*log(slipRate0*theta/L) -
	a*(1.0 - slipRate/slipRateLinear);
      mu_f = fw + (mu_f - fw) / (1. + L / theta / fwSlipRate);
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
pylith::friction::RateStateAgeingVaryingBFH::_calcFrictionDeriv(const PylithScalar t,
						      const PylithScalar slip,
						      const PylithScalar slipRate,
						      const PylithScalar normalTraction,
						      const PylithScalar* properties,
						      const int numProperties,
						      const PylithScalar* stateVars,
						      const int numStateVars)
{ // _calcFrictionDeriv
  assert(properties);
  assert(_RateStateAgeingVaryingBFH::numProperties == numProperties);
  assert(numStateVars);
  assert(_RateStateAgeingVaryingBFH::numStateVars == numStateVars);

  PylithScalar frictionDeriv = 0.0;
  if (normalTraction <= 0.0) {
    // if fault is in compression

    const PylithScalar slipRateLinear = _linearSlipRate;

    const PylithScalar a = properties[p_a];
    const PylithScalar b_initial = properties[p_b_initial];
    const PylithScalar b_final = properties[p_b_final];
    const PylithScalar b_distance = properties[p_b_distance];
    const PylithScalar L = properties[p_L];
    const PylithScalar fw = properties[p_fwcoef];
    const PylithScalar Vw = properties[p_fwSlipRate];
    const PylithScalar slipRate0 = properties[p_slipRate0];

    // Prevent zero value for theta, reasonable value is L / slipRate0
    const PylithScalar theta = (stateVars[s_state] > 0.0) ? stateVars[s_state] : L / slipRate0;

    if (slipRate >= slipRateLinear) {
      // frictionDeriv = -normalTraction * a / (slipRate * _dt);
      if (slip > b_distance) {
        frictionDeriv = -normalTraction * (a / (slipRate * _dt))/ (1. + L / theta / Vw);
      } else {
        frictionDeriv = -normalTraction * (a / (slipRate * _dt) + (b_final - b_initial) / b_distance * log(slipRate0*theta/L))/ (1. + L / theta / Vw);
      }
    } else {
      if (slip > b_distance) {
        frictionDeriv = -normalTraction * (a / (slipRateLinear * _dt)) / (1. + L / theta / Vw);
      } else {
        frictionDeriv = -normalTraction * (a / (slipRateLinear * _dt) + (b_final - b_initial) / b_distance * log(slipRate0*theta/L)) / (1. + L / theta / Vw);
      }
      
    } // else
  } // if    

  PetscLogFlops(12);

  return frictionDeriv;
} // _calcFrictionDeriv


// ----------------------------------------------------------------------
// Update state variables (for next time step).
void
pylith::friction::RateStateAgeingVaryingBFH::_updateStateVars(const PylithScalar t,
						    const PylithScalar slip,
						    const PylithScalar slipRate,
						    const PylithScalar normalTraction,
						    PylithScalar* const stateVars,
						    const int numStateVars,
						    const PylithScalar* properties,
						    const int numProperties)
{ // _updateStateVars
  assert(properties);
  assert(_RateStateAgeingVaryingBFH::numProperties == numProperties);
  assert(numStateVars);
  assert(_RateStateAgeingVaryingBFH::numStateVars == numStateVars);

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
