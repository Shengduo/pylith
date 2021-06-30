// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University at Buffalo
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2021 University of California, Davis
//
// See LICENSE.md.md for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "AuxiliaryFactoryKinematic.hh" // implementation of object methods

#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldQuery.hh" // HOLDSA FieldQuery

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "pylith/utils/error.hh" // USES PYLITH_METHOD*
#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL*

#include <cmath> // USES pow()
#include <cassert>

// ---------------------------------------------------------------------------------------------------------------------
// Default constructor.
pylith::faults::AuxiliaryFactoryKinematic::AuxiliaryFactoryKinematic(void) {
    GenericComponent::setName("auxiliaryfactorykinematic");
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor.
pylith::faults::AuxiliaryFactoryKinematic::~AuxiliaryFactoryKinematic(void) {}


// ---------------------------------------------------------------------------------------------------------------------
// Add fault strike direction subfield to auxiliary fields.
void
pylith::faults::AuxiliaryFactoryKinematic::addStrikeDir(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addStrikeDir(void)");

    const char* fieldName = "strike_dir";
    const char* componentNames[3] = { "strike_dir_x", "strike_dir_y", "strike_dir_z" };

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.alias = fieldName;
    description.vectorFieldType = pylith::topology::Field::VECTOR;
    description.numComponents = _spaceDim;
    description.componentNames.resize(_spaceDim);
    for (int i = 0; i < _spaceDim; ++i) {
        description.componentNames[i] = componentNames[i];
    } // for
    description.scale = 1.0;
    description.validator = NULL;

    _field->subfieldAdd(description, getSubfieldDiscretization(fieldName));
    // No query; computed during initialization.

    PYLITH_METHOD_END;
} // addStrikeDir


// ---------------------------------------------------------------------------------------------------------------------
// Add fault up-dip direction subfield to auxiliary fields.
void
pylith::faults::AuxiliaryFactoryKinematic::addUpDipDir(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addUpDipDir(void)");

    const char* fieldName = "up_dip_dir";
    const char* componentNames[3] = { "up_dip_dir_x", "up_dip_dir_y", "up_dip_dir_z" };

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.alias = fieldName;
    description.vectorFieldType = pylith::topology::Field::VECTOR;
    description.numComponents = _spaceDim;
    description.componentNames.resize(_spaceDim);
    for (int i = 0; i < _spaceDim; ++i) {
        description.componentNames[i] = componentNames[i];
    } // for
    description.scale = 1.0;
    description.validator = NULL;

    _field->subfieldAdd(description, getSubfieldDiscretization(fieldName));
    // No query; computed during initialization.

    PYLITH_METHOD_END;
} // addUpDipDir


// ---------------------------------------------------------------------------------------------------------------------
// Add fault normal direction subfield to auxiliary fields.
void
pylith::faults::AuxiliaryFactoryKinematic::addNormalDir(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addNormalDir(void)");

    const char* fieldName = "normal_dir";
    const char* componentNames[3] = { "normal_dir_x", "normal_dir_y", "normal_dir_z" };

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.alias = fieldName;
    description.vectorFieldType = pylith::topology::Field::VECTOR;
    description.numComponents = _spaceDim;
    description.componentNames.resize(_spaceDim);
    for (int i = 0; i < _spaceDim; ++i) {
        description.componentNames[i] = componentNames[i];
    } // for
    description.scale = 1.0;
    description.validator = NULL;

    _field->subfieldAdd(description, getSubfieldDiscretization(fieldName));
    // No query; computed during initialization.

    PYLITH_METHOD_END;
} // addNormalDir


// ---------------------------------------------------------------------------------------------------------------------
// Add fault slip subfield to auxiliary fields.
void
pylith::faults::AuxiliaryFactoryKinematic::addSlip(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addSlip(void)");

    const char* fieldName = "slip";
    const char* componentNames[3] = { "slip_opening", "slip_left_lateral", "slip_reverse" };

    const PylithReal lengthScale = _normalizer->getLengthScale();

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.alias = fieldName;
    description.vectorFieldType = pylith::topology::Field::VECTOR;
    description.numComponents = _spaceDim;
    description.componentNames.resize(_spaceDim);
    for (int i = 0; i < _spaceDim; ++i) {
        description.componentNames[i] = componentNames[i];
    } // for
    description.scale = lengthScale;
    description.validator = NULL;

    _field->subfieldAdd(description, getSubfieldDiscretization(fieldName));
    // No query; populated by kinematic source at beginning of time step.

    PYLITH_METHOD_END;
} // addSlip


// ---------------------------------------------------------------------------------------------------------------------
// Add fault slip rate subfield to auxiliary fields.
void
pylith::faults::AuxiliaryFactoryKinematic::addSlipRate(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addSlipRate(void)");

    const char* fieldName = "slip_rate";
    const char* componentNames[3] = { "slip_rate_opening", "slip_rate_left_lateral", "slip_rate_reverse" };

    const PylithReal velocityScale = _normalizer->getLengthScale() / _normalizer->getTimeScale();

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.alias = fieldName;
    description.vectorFieldType = pylith::topology::Field::VECTOR;
    description.numComponents = _spaceDim;
    description.componentNames.resize(_spaceDim);
    for (int i = 0; i < _spaceDim; ++i) {
        description.componentNames[i] = componentNames[i];
    } // for
    description.scale = velocityScale;
    description.validator = NULL;

    _field->subfieldAdd(description, getSubfieldDiscretization(fieldName));
    // No query; populated by kinematic source at beginning of time step.

    PYLITH_METHOD_END;
} // addSlipRate


// ---------------------------------------------------------------------------------------------------------------------
// Add fault slip acceleration subfield to auxiliary fields.
void
pylith::faults::AuxiliaryFactoryKinematic::addSlipAcceleration(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addSlipAcceleration(void)");

    const char* fieldName = "slip_acceleration";
    const char* componentNames[3] = {
        "slip_acceleration_opening",
        "slip_acceleration_left_lateral",
        "slip_acceleration_reverse",
    };

    const PylithReal accelerationScale = _normalizer->getLengthScale() / pow(_normalizer->getTimeScale(), 2);

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.alias = fieldName;
    description.vectorFieldType = pylith::topology::Field::VECTOR;
    description.numComponents = _spaceDim;
    description.componentNames.resize(_spaceDim);
    for (int i = 0; i < _spaceDim; ++i) {
        description.componentNames[i] = componentNames[i];
    } // for
    description.scale = accelerationScale;
    description.validator = NULL;

    _field->subfieldAdd(description, getSubfieldDiscretization(fieldName));
    // No query; populated by kinematic source at beginning of time step.

    PYLITH_METHOD_END;
} // addSlipAcc
// --------------------------------------------------------------------
// Add undrained bulk modulus subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryKinematic::addUndrainedBulkModulus(void) { // UndrainedBulkModulus
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addUndrainedBulkModulus(void)");

    const char* subfieldName = "undrained_bulk_modulus";
    const PylithReal pressureScale = _normalizer->getPressureScale();

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = pressureScale;
    description.validator = pylith::topology::FieldQuery::validatorPositive;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addUndrainedBulkModulus


// ---------------------------------------------------------------------------------------------------------------------
// Add shear modulus subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryKinematic::addShearModulus(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addShearModulus(void)");

    const char* subfieldName = "shear_modulus";
    const PylithReal pressureScale = _normalizer->getPressureScale();

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = pressureScale;
    description.validator = pylith::topology::FieldQuery::validatorNonnegative;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addShearModulus


// ---------------------------------------------------------------------
// Add Skempton coefficient subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryKinematic::addSkemptonCoefficient(void) { // skemptonCoefficient
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addSkemptonCoefficient(void)");

    const char* subfieldName = "skempton_coefficient";
    const PylithReal noScale = 1;

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = noScale;
    description.validator = pylith::topology::FieldQuery::validatorPositive;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    pylith::faults::Query::skemptonCoefficientFromInput(subfieldName, this);

    PYLITH_METHOD_END;
} // addSkemptonCoefficient
// Add auxilliary fields for FaultPoroDiffusionCohesiveKin
// Add auxilliary fields for FaultPoroDiffusionCohesiveKin

// Add layer thickness to auxiliary fields.
void
pylith::faults::AuxiliaryFactoryKinematic::addThickness(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addThickness(void)");

    const char* subfieldName = "thickness";

    const PylithReal lengthScale = _normalizer->getLengthScale();

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName
    description.scale = lengthScale;
    description.validator = pylith::topology::FieldQuery::validatorPositive;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));

    // ** TO DO **
    // Is this step correct?
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addThickness

// Add layer porosity to auxiliary fields.
void
pylith::faults::AuxiliaryFactoryKinematic::addPorosity(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addPorosity(void)");

    const char* subfieldName = "porosity";

    const PylithReal noScale = 1;

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName
    description.scale = noScale;
    description.validator = pylith::topology::FieldQuery::validatorPositive;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));

    // ** TO DO **
    // Is this step correct?
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addPorosity

// Add beta p to auxiliary fields.
void
pylith::faults::AuxiliaryFactoryKinematic::addBetaP(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addBetaP(void)");

    const char* subfieldName = "beta_p";

    // ** TO DO **
    // This works? The scale should be 1/Pa
    const PylithReal betaScale = 1. / _normalizer->getPressureScale();

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName
    description.scale = betaScale;
    description.validator = pylith::topology::FieldQuery::validatorPositive;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));

    // ** TO DO **
    // Is this step correct?
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addBetaP

// Add beta_sigma to auxiliary fields.
void
pylith::faults::AuxiliaryFactoryKinematic::addBetaSigma(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addBetaSigma(void)");

    const char* subfieldName = "beta_sigma";

    const PylithReal betaScale = 1. / _normalizer->getPressureScale();

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName
    description.scale = betaScale;
    description.validator = pylith::topology::FieldQuery::validatorPositive;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));

    // ** TO DO **
    // Is this step correct?
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addBetaSigma

// Add tangential permeability to auxiliary fields.
void
pylith::faults::AuxiliaryFactoryKinematic::addPermeabilityTangential(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addPermeabilityTangential(void)");

    const char* subfieldName = "permeability_tangential";
    
    // ** TO DO **
    // Please verify this following line
    // Permeability scale should be m^2
    const PylithReal permeabilityScale = (_normalizer->getLengthScale())^2;

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName
    description.scale = permeabilityScale;
    description.validator = pylith::topology::FieldQuery::validatorPositive;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));

    // ** TO DO **
    // Is this step correct?
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addPermeabilityTangential

// Add normal permeability to auxiliary fields.
void
pylith::faults::AuxiliaryFactoryKinematic::addPermeabilityNormal(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addPermeabilityNormal(void)");

    const char* subfieldName = "permeability_normal";

    // ** TO DO **
    // Please verify this following line
    // Permeability scale should be m^2
    const PylithReal permeabilityScale = (_normalizer->getLengthScale())^2;

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName
    description.scale = permeabilityScale;
    description.validator = pylith::topology::FieldQuery::validatorPositive;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));

    // ** TO DO **
    // Is this step correct?
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addPermeabilityNormal

// Add fluid viscosity to auxiliary fields.
void
pylith::faults::AuxiliaryFactoryKinematic::addFluidViscosity(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addFluidViscosity(void)");

    const char* subfieldName = "fluid_viscosity";

    // ** TO DO **
    // Please verify this following line
    // viscosity scale should be pa s
    const PylithReal fluidViscosityScale = _normalizer->getPressureScale()
                                           * _normalizer->getTimeScale();

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName
    description.scale = fluiViscosityScale;
    description.validator = pylith::topology::FieldQuery::validatorPositive;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));

    // ** TO DO **
    // Is this step correct?
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addFluidViscosity

// Add negative side bulk modulus to auxiliary fields.
void
pylith::faults::AuxiliaryFactoryKinematic::addBulkModulusNegative(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addBulkModulusNegative(void)");

    const char* subfieldName = "bulk_modulus_negative";

    const PylithReal bulkModulusScale = _normalizer->getPressureScale();

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName
    description.scale = bulkModulusScale;
    description.validator = pylith::topology::FieldQuery::validatorPositive;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));

    // ** TO DO **
    // Is this step correct?
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addBulkModulusNegative

// Add positive side bulk modulus to auxiliary fields.
void
pylith::faults::AuxiliaryFactoryKinematic::addBulkModulusPositive(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addBulkModulusPositive(void)");

    const char* subfieldName = "bulk_modulus_positive";

    const PylithReal bulkModulusScale = _normalizer->getPressureScale();

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName
    description.scale = bulkModulusScale;
    description.validator = pylith::topology::FieldQuery::validatorPositive;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));

    // ** TO DO **
    // Is this step correct?
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addBulkModulusPositive

// Add negative side shear modulus to auxiliary fields.
void
pylith::faults::AuxiliaryFactoryKinematic::addShearModulusNegative(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addShearModulusNegative(void)");

    const char* subfieldName = "shear_modulus_negative";

    const PylithReal shearModulusScale = _normalizer->getPressureScale();

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName
    description.scale = shearModulusScale;
    description.validator = pylith::topology::FieldQuery::validatorPositive;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));

    // ** TO DO **
    // Is this step correct?
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addShearModulusNegative

// Add positive side shear modulus to auxiliary fields.
void
pylith::faults::AuxiliaryFactoryKinematic::addShearModulusPositive(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addShearModulusPositive(void)");

    const char* subfieldName = "shear_modulus_positive";

    const PylithReal shearModulusScale = _normalizer->getPressureScale();

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName
    description.scale = shearModulusScale;
    description.validator = pylith::topology::FieldQuery::validatorPositive;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));

    // ** TO DO **
    // Is this step correct?
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addShearModulusPositive

// Add body force subfield to auxiliary fields.
void
pylith::faults::AuxiliaryFactoryKinematic::addBodyForce(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addBodyForce(void)");

    const char* fieldName = "body_force";
    const char* componentNames[3] = { "body_force_x", "body_force_y", "body_force_z" };

    // ** TO DO **
    // Verify this line
    // The scale should be pa / m
    const PylithReal bodyForceScale = normalizer->getPressureScale() / _normalizer->getLengthScale();

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.alias = fieldName;
    description.vectorFieldType = pylith::topology::Field::VECTOR;
    description.numComponents = _spaceDim;
    description.componentNames.resize(_spaceDim);
    for (int i = 0; i < _spaceDim; ++i) {
        description.componentNames[i] = componentNames[i];
    } // for
    description.scale = bodyForceScale;
    description.validator = NULL;

    _field->subfieldAdd(description, getSubfieldDiscretization(fieldName));
    // No query; populated by kinematic source at beginning of time step.

    // ** TO DO **
    // Is this step correct?
    this->setSubfieldQuery(fieldName);

    PYLITH_METHOD_END;
} // addBodyForce

// Add source to auxiliary fields.
void
pylith::faults::AuxiliaryFactoryKinematic::addSource(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addSource(void)");

    const char* subfieldName = "source";
    // ** TO DO **
    // Verify this line
    // The scale for source is 1 / s
    const PylithReal sourceScale = 1. / _normalizer->getTimeScale();

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName
    description.scale = sourceScale;
    description.validator = pylith::topology::FieldQuery::validatorNonnegative;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));

    // ** TO DO **
    // Is this step correct?
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addSource

// End of file
