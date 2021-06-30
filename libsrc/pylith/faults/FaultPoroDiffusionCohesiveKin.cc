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

#include "FaultPoroDiffusionCohesiveKin.hh" // implementation of object methods

#include "pylith/faults/KinSrc.hh" // USES KinSrc
#include "pylith/faults/AuxiliaryFactoryKinematic.hh" // USES AuxiliaryFactoryKinematic
#include "pylith/feassemble/IntegratorInterface.hh" // USES IntegratorInterface
#include "pylith/feassemble/ConstraintSimple.hh" // USES ConstraintSimple

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldOps.hh" // USES FieldOps::checkDiscretization()
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh

#include "pylith/fekernels/FaultPoroDiffusionCohesiveKin.hh" // USES FaultPoroDiffusionCohesiveKin

#include "pylith/utils/EventLogger.hh" // USES EventLogger
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensionalizer
#include "spatialdata/spatialdb/SpatialDB.hh" // USES SpatialDB

#include <cmath> // USES pow(), sqrt()
#include <strings.h> // USES strcasecmp()
#include <cstring> // USES strlen()
#include <cstdlib> // USES atoi()
#include <cassert> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error
#include <typeinfo> // USES typeid()

// ---------------------------------------------------------------------------------------------------------------------
typedef pylith::feassemble::IntegratorInterface::ResidualKernels ResidualKernels;
typedef pylith::feassemble::IntegratorInterface::JacobianKernels JacobianKernels;

// ---------------------------------------------------------------------------------------------------------------------
namespace pylith {
    namespace faults {
        class _FaultPoroDiffusionCohesiveKin {
            // PUBLIC MEMBERS //////////////////////////////////////////////////////////////////////////////////////////
public:

            /** Set kernels for LHS residual.
             *
             * @param[out] integrator Integrator for interface.
             * @param[in] fault Fault object for kinematic ruptures.
             * @param[in] solution Solution field.
             * @param[in] formulation Formulation for equations.
             */
            static void setKernelsLHSResidual(pylith::feassemble::IntegratorInterface *integrator,
                                              const pylith::faults::FaultPoroDiffusionCohesiveKin &fault,
                                              const pylith::topology::Field &solution,
                                              const pylith::problems::Physics::FormulationEnum formulation);

            /** Set kernels for LHS Jacobian.
             *
             * @param[out] integrator Integrator for interface.
             * @param[in] fault Fault object for kinematic ruptures.
             * @param[in] solution Solution field.
             * @param[in] formulation Formulation for equations.
             */
            static void setKernelsLHSJacobian(pylith::feassemble::IntegratorInterface *integrator,
                                              const pylith::faults::FaultPoroDiffusionCohesiveKin &fault,
                                              const pylith::topology::Field &solution,
                                              const pylith::problems::Physics::FormulationEnum formulation);

            static const char *pyreComponent;
        };
        const char *_FaultPoroDiffusionCohesiveKin::pyreComponent = "faultporodiffusioncohesivekin";

        // _FaultPoroDiffusionCohesiveKin

    } // faults
} // pylith

// ---------------------------------------------------------------------------------------------------------------------
// Default constructor.
pylith::faults::FaultPoroDiffusionCohesiveKin::FaultPoroDiffusionCohesiveKin(void) : _auxiliaryFactory(new pylith::faults::AuxiliaryFactoryKinematic),
    _slipVecRupture(NULL),
    _slipVecTotal(NULL) {
    pylith::utils::PyreComponent::setName(_FaultPoroDiffusionCohesiveKin::pyreComponent);
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor.
pylith::faults::FaultPoroDiffusionCohesiveKin::~FaultPoroDiffusionCohesiveKin(void) {
    deallocate();
} // destructor


// ---------------------------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::faults::FaultPorodiffusionCohesiveKin::deallocate(void) {
    FaultCohesive::deallocate();

    PetscErrorCode err = VecDestroy(&_slipVecRupture);
    PYLITH_CHECK_ERROR(err);
    err = VecDestroy(&_slipVecTotal);
    PYLITH_CHECK_ERROR(err);
    delete _auxiliaryFactory;
    _auxiliaryFactory = NULL;
    _ruptures.clear(); // :TODO: Use shared pointers for earthquake ruptures
} // deallocate


// ---------------------------------------------------------------------------------------------------------------------
// Set kinematic earthquake ruptures.
void
pylith::faults::FaultPoroDiffusionCohesiveKin::setEqRuptures(const char *const *names,
                                                             const int numNames,
                                                             KinSrc **ruptures,
                                                             const int numRuptures) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("setEqRuptures(names=" << names << ", numNames=" << numNames << ", ruptures=" << ruptures << ", numRuptures=" << numRuptures << ")");

    assert(numNames == numRuptures);

    // :TODO: Use shared pointers for earthquake ruptures
    _ruptures.clear();
    for (int i = 0; i < numRuptures; ++i) {
        if (!ruptures[i]) {
            std::ostringstream msg;
            msg << "Null earthquake rupture object for earthquake rupture '" << names[i] << "'.";
            throw std::runtime_error(msg.str());
        } // if
        _ruptures[std::string(names[i])] = ruptures[i];
    } // for

    PYLITH_METHOD_END;
} // setEqRuptures


// ---------------------------------------------------------------------------------------------------------------------
// Include body force?
void
pylith::faults::FaultPoroDiffusionCohesiveKin::useBodyForce(const bool value) {
    PYLITH_COMPONENT_DEBUG("useBodyForce(value=" << value << ")");

    _useBodyForce = value;
} // useBodyForce


// ---------------------------------------------------------------------------------------------------------------------
// Include body force?
bool
pylith::faults::FaultPoroDiffusionCohesiveKin::useBodyForce(void) const {
    return _useBodyForce;
} // useBodyForce


// ----------------------------------------------------------------------
// Include source?
void
pylith::faults::FaultPoroDiffusionCohesiveKin::useSource(const bool value) {
    PYLITH_COMPONENT_DEBUG("useSource(value=" << value << ")");

    _useSource = value;
} // useSource


// ----------------------------------------------------------------------
// Include source density?
bool
pylith::faults::FaultPoroDiffusionCohesiveKin::useSource(void) const {
    return _useSource;
} // useSource


// ----------------------------------------------------------------------
// Include constant pressure source?
void
pylith::faults::FaultPoroDiffusionCohesiveKin::useConstantPressureSource(const bool value) {
    PYLITH_COMPONENT_DEBUG("useConstantPressureSource(value=" << value << ")");

    _useConstantPressureSource = value;
} // useConstantPressureSource


// ----------------------------------------------------------------------
// Include constant pressure source?
bool
pylith::faults::FaultPoroDiffusionCohesiveKin::useConstantPressureSource(void) const {
    return _useConstantPressureSource;
} // useConstantPressureSource


// ---------------------------------------------------------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::faults::FaultPoroDiffusionCohesiveKin::verifyConfiguration(const pylith::topology::Field &solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("verifyConfiguration(solution=" << solution.getLabel() << ")");

    if (!solution.hasSubfield("lagrange_multiplier_fault")) {
        std::ostringstream msg;
        msg << "Cannot find 'lagrange_multiplier_fault' subfield in solution field for poroelastic diffusive fault implementation in component '"
            << PyreComponent::getIdentifier() << "'.";
        throw std::runtime_error(msg.str());
    } // if

    /** TO DO **
     * Implement or verify field::hasSubfield("pressure", "trace_strain", "fault_pressure") is implemented
     * Seems that this is fine, at least no need to change the function field::hasSubfield
     */
    if (!solution.hasSubfield("pressure")) {
        std::ostringstream msg;
        msg << "Cannot find 'pressure' subfield in solution field for poroelastic diffusive fault implementation in component '"
            << PyreComponent::getIdentifier() << "'.";
        throw std::runtime_error(msg.str());
    } // if

    if (!solution.hasSubfield("fault_pressure")) {
        std::ostringstream msg;
        msg << "Cannot find 'fault_pressure' subfield in solution field for poroelastic diffusive fault implementation in component '"
            << PyreComponent::getIdentifier() << "'.";
        throw std::runtime_error(msg.str());
    } // if

    switch (_formulation) {
    case QUASISTATIC:
        if (!solution.hasSubfield("displacement")) {
            std::ostringstream msg;
            msg << "Cannot find 'displacement' subfield in solution field for fault implementation in component '"
                << PyreComponent::getIdentifier() << "'.";
            throw std::runtime_error(msg.str());
        } // if
        break;
    case DYNAMIC_IMEX:
        if (!solution.hasSubfield("velocity")) {
            std::ostringstream msg;
            msg << "Cannot find 'velocity' subfield in solution field for fault implementation in component '"
                << PyreComponent::getIdentifier() << "'.";
            throw std::runtime_error(msg.str());
        } // if
        break;
    case DYNAMIC:
        PYLITH_COMPONENT_LOGICERROR("Fault implementation is incompatible with 'dynamic' formulation. Use 'dynamic_imex'.");
        break;
    default:
        PYLITH_COMPONENT_LOGICERROR("Unknown formulation for equations (" << _formulation << ").");
    } // switch

    PYLITH_METHOD_END;
} // verifyConfiguration


// ---------------------------------------------------------------------------------------------------------------------
// Create integrator and set kernels.
pylith::feassemble::Integrator *
pylith::faults::FaultPoroDiffusionCohesiveKin::createIntegrator(const pylith::topology::Field &solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("createIntegrator(solution=" << solution.getLabel() << ")");

    pylith::feassemble::IntegratorInterface *integrator = new pylith::feassemble::IntegratorInterface(this);
    assert(integrator);
    integrator->setLabelValue(getInterfaceId());
    integrator->setSurfaceMarkerLabel(getSurfaceMarkerLabel());

    _FaultPoroDiffusionCohesiveKin::setKernelsLHSResidual(integrator, *this, solution, _formulation);
    _FaultPoroDiffusionCohesiveKin::setKernelsLHSJacobian(integrator, *this, solution, _formulation);
    // No state variables.
    // _FaultPoroDiffusionCohesiveKin::setKernelsDerivedFields(integrator, *this, solution);

    PYLITH_METHOD_RETURN(integrator);
} // createIntegrator


// ---------------------------------------------------------------------------------------------------------------------
// Create constraint for buried fault edges and faces.
// ** TO DO **
// Check if new constraint needed for poro-fault
pylith::feassemble::Constraint *
pylith::faults::FaultPoroDiffusionCohesiveKin::createConstraint(const pylith::topology::Field &solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("createConstraint(solution=" << solution.getLabel() << ")");

    if (0 == strlen(getBuriedEdgesMarkerLabel())) {
        PYLITH_METHOD_RETURN(NULL);
    } // if

    const char *lagrangeName = "lagrange_multiplier_fault";
    // const PylithInt numComponents = solution.subfieldInfo(lagrangeName).fe.numComponents;
    const PylithInt numComponents = solution.getSpaceDim();

    pylith::int_array constrainedDOF;
    constrainedDOF.resize(numComponents);
    for (int c = 0; c < numComponents; ++c) {
        constrainedDOF[c] = c;
    }
    // Make new label for cohesive edges and faces
    PetscDM dm = solution.dmMesh();
    PetscDMLabel buriedLabel = NULL;
    PetscDMLabel buriedCohesiveLabel = NULL;
    PetscIS pointIS = NULL;
    const PetscInt *points = NULL;
    PetscInt n;
    std::ostringstream labelstream;
    labelstream << getBuriedEdgesMarkerLabel() << "_cohesive";
    std::string labelname = labelstream.str();
    PetscErrorCode err;

    err = DMCreateLabel(dm, labelname.c_str());
    PYLITH_CHECK_ERROR(err);
    err = DMGetLabel(dm, getBuriedEdgesMarkerLabel(), &buriedLabel);
    PYLITH_CHECK_ERROR(err);
    err = DMGetLabel(dm, labelname.c_str(), &buriedCohesiveLabel);
    PYLITH_CHECK_ERROR(err);
    err = DMLabelGetStratumIS(buriedLabel, 1, &pointIS);
    PYLITH_CHECK_ERROR(err);
    err = ISGetLocalSize(pointIS, &n);
    PYLITH_CHECK_ERROR(err);
    err = ISGetIndices(pointIS, &points);
    PYLITH_CHECK_ERROR(err);
    for (int p = 0; p < n; ++p) {
        const PetscInt *support = NULL;
        PetscInt supportSize;

        err = DMPlexGetSupportSize(dm, points[p], &supportSize);
        PYLITH_CHECK_ERROR(err);
        err = DMPlexGetSupport(dm, points[p], &support);
        PYLITH_CHECK_ERROR(err);
        for (int s = 0; s < supportSize; ++s) {
            DMPolytopeType ct;
            const PetscInt spoint = support[s];

            err = DMPlexGetCellType(dm, spoint, &ct);
            PYLITH_CHECK_ERROR(err);
            if ((ct == DM_POLYTOPE_SEG_PRISM_TENSOR) || (ct == DM_POLYTOPE_POINT_PRISM_TENSOR)) {
                const PetscInt *cone = NULL;
                PetscInt coneSize;

                err = DMPlexGetConeSize(dm, spoint, &coneSize);
                PYLITH_CHECK_ERROR(err);
                err = DMPlexGetCone(dm, spoint, &cone);
                PYLITH_CHECK_ERROR(err);
                for (int c = 0; c < coneSize; ++c) {
                    PetscInt val;
                    err = DMLabelGetValue(buriedLabel, cone[c], &val);
                    PYLITH_CHECK_ERROR(err);
                    if (val >= 0) {
                        err = DMLabelSetValue(buriedCohesiveLabel, spoint, 1);
                        PYLITH_CHECK_ERROR(err);
                        break;
                    }
                } // for
            } // if
        } // for
    } // for

    pylith::feassemble::ConstraintSimple *constraint = new pylith::feassemble::ConstraintSimple(this);
    assert(constraint);
    constraint->setMarkerLabel(labelname.c_str());
    err = PetscObjectViewFromOptions((PetscObject)buriedLabel, NULL, "-buried_edge_label_view");
    err = PetscObjectViewFromOptions((PetscObject)buriedCohesiveLabel, NULL, "-buried_cohesive_edge_label_view");
    constraint->setConstrainedDOF(&constrainedDOF[0], constrainedDOF.size());
    constraint->setSubfieldName(lagrangeName);
    constraint->setUserFn(_zero);

    PYLITH_METHOD_RETURN(constraint);
} // createConstraint


// ---------------------------------------------------------------------------------------------------------------------
// Create auxiliary field.
pylith::topology::Field *
pylith::faults::FaultPoroDiffusionCohesiveKin::createAuxiliaryField(const pylith::topology::Field &solution,
                                                                    const pylith::topology::Mesh &domainMesh) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("createAuxiliaryField(solution=" << solution.getLabel() << ", domainMesh=)" << typeid(domainMesh).name() << ")");

    assert(_normalizer);

    pylith::topology::Field *auxiliaryField = new pylith::topology::Field(domainMesh);
    assert(auxiliaryField);
    auxiliaryField->setLabel("FaultPoroDiffusionCohesiveKin auxiliary field");

    // Set default discretization of auxiliary subfields to match lagrange_multiplier_fault subfield in solution.
    assert(_auxiliaryFactory);
    const pylith::topology::FieldBase::Discretization &discretization = solution.subfieldInfo("lagrange_multiplier_fault").fe;
    _auxiliaryFactory->setSubfieldDiscretization("default", discretization.basisOrder, discretization.quadOrder, -1,
                                                 discretization.cellBasis, discretization.isBasisContinuous, discretization.feSpace);

    // ** TO DO **
    // Set discretization of aux fields to match also
    // pressure (\Gamma^+ and - sides), trace_strain (\Gamma^+ and - sides), fault_pressure (only \Gamma^f)

    assert(_auxiliaryFactory);
    assert(_normalizer);
    _auxiliaryFactory->initialize(auxiliaryField, *_normalizer, solution.getSpaceDim());

    // ** TO DO **
    // implement the following add* functions
    /** -- numA : number of auxiliary fields
     ***** Required fields
     * - 0: thickness(1)
     * - 1: porosity(1)
     * - 2: beta_p(1)
     * - 3: beta_sigma(1)
     * - 4: permeability_tangential(1)
     * - 5: permeability_normal(1)
     * - 6: fluid_viscosity(1)
     * - 7: bulk_modulus_negative(1)
     * - 8: shear_modulus_negative(1)
     * - 9: bulk_modulus_positive(1)
     * - 10: shear_modulus_positive(1)
     * - numA - 3: body_force(dim)
     * - numA - 2: source(1)
     * - numA - 1: slip(dim)
     */
    _auxiliaryFactory->addThickness(); // 0
    _auxiliaryFactory->addPorosity(); // 1
    _auxiliaryFactory->addBetaP(); // 2
    _auxiliaryFactory->addBetaSigma(); // 3
    _auxiliaryFactory->addPermeabilityTangential(); // 4
    _auxiliaryFactory->addPermeabilityNormal(); // 5
    _auxiliaryFactory->addFluidViscosity(); // 6
    _auxiliaryFactory->addBulkModulusNegative(); // 7
    _auxiliaryFactory->addShearModulusNegative(); // 8
    _auxiliaryFactory->addBulkModulusPositive(); // 9
    _auxiliaryFactory->addShearModulusPositive(); // 10
    if (_useBodyForce) {
        _auxiliaryFactory->addBodyForce(); // +1
    } // if
    if (_useSource) {
        _auxiliaryFactory->addSource(); // +1
    } // if
    if (_useConstantPressureSource) {
        auxiliaryFactory->addConstantPressureSource(); // +1
    }
    // :ATTENTION: The order for adding subfields must match the order of the auxiliary fields in the FE kernels.

    switch (_formulation) {
    case QUASISTATIC:
        _auxiliaryFactory->addSlip(); // numA - 1
        break;
    case DYNAMIC_IMEX:
        _auxiliaryFactory->addSlipRate(); // numA - 1
        break;
    case DYNAMIC:
        PYLITH_COMPONENT_LOGICERROR("Fault implementation is incompatible with 'dynamic' time-stepping formulation. Use 'dynamic_imex'.");
        break;
    default:
        PYLITH_COMPONENT_LOGICERROR("Unknown formulation for equations (" << _formulation << ").");
    } // switch

    auxiliaryField->subfieldsSetup();
    auxiliaryField->createDiscretization();
    pylith::topology::FieldOps::checkDiscretization(solution, *auxiliaryField);
    auxiliaryField->allocate();
    auxiliaryField->createOutputVector();

    // We don't populate the auxiliary field via a spatial database, because they will be set from the earthquake
    // rupture.

    // Initialize auxiliary fields for kinematic ruptures.
    assert(auxiliaryField);
    const srcs_type::const_iterator rupturesEnd = _ruptures.end();
    for (srcs_type::iterator r_iter = _ruptures.begin(); r_iter != rupturesEnd; ++r_iter) {
        KinSrc *src = r_iter->second;
        assert(src);
        src->initialize(*auxiliaryField, *_normalizer, solution.mesh().getCoordSys());
    } // for

    // Create local PETSc vector to hold current slip.
    PetscErrorCode err = 0;
    err = DMCreateLocalVector(auxiliaryField->dmMesh(), &_slipVecRupture);
    PYLITH_CHECK_ERROR(err);
    err = DMCreateLocalVector(auxiliaryField->dmMesh(), &_slipVecTotal);
    PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_RETURN(auxiliaryField);
} // createAuxiliaryField


// ---------------------------------------------------------------------------------------------------------------------
// Create derived field.
pylith::topology::Field *
pylith::faults::FaultPoroDiffusionCohesiveKin::createDerivedField(const pylith::topology::Field &solution,
                                                                  const pylith::topology::Mesh &domainMesh) {
    return NULL;
} // createDerivedField


// ---------------------------------------------------------------------------------------------------------------------
// Update auxiliary fields at beginning of time step.
void
pylith::faults::FaultPoroDiffusionCohesiveKin::updateAuxiliaryField(pylith::topology::Field *auxiliaryField,
                                                                    const double t) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("updateAuxiliaryField(auxiliaryField=" << auxiliaryField << ", t=" << t << ")");

    assert(auxiliaryField);
    assert(_normalizer);

    switch (_formulation) {
    case QUASISTATIC:
        this->_updateSlip(auxiliaryField, t);
        break;
    case DYNAMIC_IMEX:
        this->_updateSlipRate(auxiliaryField, t);
        break;
    case DYNAMIC:
        PYLITH_COMPONENT_LOGICERROR("Fault implementation is incompatible with 'dynamic' formulation. Use 'dynamic_imex'.");
        break;
    default:
        PYLITH_COMPONENT_LOGICERROR("Unknown formulation for equations (" << _formulation << ").");
    } // switch

    PYLITH_METHOD_END;
} // updateAuxiliaryField


// ---------------------------------------------------------------------------------------------------------------------
// Get auxiliary factory associated with physics.
pylith::feassemble::AuxiliaryFactory *
pylith::faults::FaultPoroDiffusionCohesiveKin::_getAuxiliaryFactory(void) {
    return _auxiliaryFactory;
} // _getAuxiliaryFactory


// ---------------------------------------------------------------------------------------------------------------------
// Update kernel constants.
void
pylith::faults::FaultPoroDiffusionCohesiveKin::_updateKernelConstants(const PylithReal dt) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_setKernelConstants(dt=" << dt << ")");

    if (6 != _kernelConstants.size()) {
        _kernelConstants.resize(6);
    }
    _kernelConstants[0] = _refDir1[0];
    _kernelConstants[1] = _refDir1[1];
    _kernelConstants[2] = _refDir1[2];
    _kernelConstants[3] = _refDir2[0];
    _kernelConstants[4] = _refDir2[1];
    _kernelConstants[5] = _refDir2[2];

    PYLITH_METHOD_END;
} // _updateKernelConstants


// ---------------------------------------------------------------------------------------------------------------------
// Update slip subfield in auxiliary field at beginning of time step.
void
pylith::faults::FaultPoroDiffusionCohesiveKin::_updateSlip(pylith::topology::Field *auxiliaryField,
                                                           const double t) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("updateSlip(auxiliaryField=" << auxiliaryField << ", t=" << t << ")");

    assert(auxiliaryField);
    assert(_normalizer);

    // Update slip subfield at current time step
    PetscErrorCode err = VecSet(_slipVecTotal, 0.0);
    PYLITH_CHECK_ERROR(err);
    const srcs_type::const_iterator rupturesEnd = _ruptures.end();
    for (srcs_type::iterator r_iter = _ruptures.begin(); r_iter != rupturesEnd; ++r_iter) {
        err = VecSet(_slipVecRupture, 0.0);
        PYLITH_CHECK_ERROR(err);

        KinSrc *src = r_iter->second;
        assert(src);
        src->updateSlip(_slipVecRupture, auxiliaryField, t, _normalizer->getTimeScale());
        err = VecAYPX(_slipVecTotal, 1.0, _slipVecRupture);
    } // for

    // Transfer slip values from local PETSc slip vector to fault auxiliary field.
    PetscInt pStart = 0, pEnd = 0;
    err = PetscSectionGetChart(auxiliaryField->localSection(), &pStart, &pEnd);
    PYLITH_CHECK_ERROR(err);

    pylith::topology::VecVisitorMesh auxiliaryVisitor(*auxiliaryField, "slip");
    PylithScalar *auxiliaryArray = auxiliaryVisitor.localArray();

    const PylithScalar *slipArray = NULL;
    err = VecGetArrayRead(_slipVecTotal, &slipArray);
    PYLITH_CHECK_ERROR(err);

    for (PetscInt p = pStart, iSlip = 0; p < pEnd; ++p) {
        const PetscInt slipDof = auxiliaryVisitor.sectionDof(p);
        const PetscInt slipOff = auxiliaryVisitor.sectionOffset(p);
        for (PetscInt iDof = 0; iDof < slipDof; ++iDof, ++iSlip) {
            auxiliaryArray[slipOff + iDof] = slipArray[iSlip];
        } // for
    } // for
    err = VecRestoreArrayRead(_slipVecTotal, &slipArray);
    PYLITH_CHECK_ERROR(err);

    pythia::journal::debug_t debug(pylith::utils::PyreComponent::getName());
    if (debug.state()) {
        auxiliaryField->view("Fault auxiliary field after setting slip.");
    } // if

    PYLITH_METHOD_END;
} // _updateSlip


// ---------------------------------------------------------------------------------------------------------------------
// Update slip rate subfield in auxiliary field at beginning of time step.
void
pylith::faults::FaultPoroDiffusionCohesiveKin::_updateSlipRate(pylith::topology::Field *auxiliaryField,
                                                               const double t) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_updateSlipRate(auxiliaryField=" << auxiliaryField << ", t=" << t << ")");

    assert(auxiliaryField);
    assert(_normalizer);

    // Update slip rate subfield at current time step
    PetscErrorCode err = VecSet(_slipVecTotal, 0.0);
    PYLITH_CHECK_ERROR(err);
    const srcs_type::const_iterator rupturesEnd = _ruptures.end();
    for (srcs_type::iterator r_iter = _ruptures.begin(); r_iter != rupturesEnd; ++r_iter) {
        err = VecSet(_slipVecRupture, 0.0);
        PYLITH_CHECK_ERROR(err);

        KinSrc *src = r_iter->second;
        assert(src);
        src->updateSlipRate(_slipVecRupture, auxiliaryField, t, _normalizer->getTimeScale());
        err = VecAYPX(_slipVecTotal, 1.0, _slipVecRupture);
    } // for

    // Transfer slip values from local PETSc slip vector to fault auxiliary field.
    PetscInt pStart = 0, pEnd = 0;
    err = PetscSectionGetChart(auxiliaryField->localSection(), &pStart, &pEnd);
    PYLITH_CHECK_ERROR(err);

    pylith::topology::VecVisitorMesh auxiliaryVisitor(*auxiliaryField, "slip_rate");
    PylithScalar *auxiliaryArray = auxiliaryVisitor.localArray();

    const PylithScalar *slipRateArray = NULL;
    err = VecGetArrayRead(_slipVecTotal, &slipRateArray);
    PYLITH_CHECK_ERROR(err);

    for (PetscInt p = pStart, iSlipRate = 0; p < pEnd; ++p) {
        const PetscInt slipRateDof = auxiliaryVisitor.sectionDof(p);
        const PetscInt slipRateOff = auxiliaryVisitor.sectionOffset(p);
        for (PetscInt iDof = 0; iDof < slipRateDof; ++iDof, ++iSlipRate) {
            auxiliaryArray[slipRateOff + iDof] = slipRateArray[iSlipRate];
        } // for
    } // for
    err = VecRestoreArrayRead(_slipVecTotal, &slipRateArray);
    PYLITH_CHECK_ERROR(err);

    pythia::journal::debug_t debug(pylith::utils::PyreComponent::getName());
    if (debug.state()) {
        auxiliaryField->view("Fault auxiliary field after setting slip rate.");
    } // if

    PYLITH_METHOD_END;
} // _updateSlipRate


// ---------------------------------------------------------------------------------------------------------------------
// Set kernels for LHS residual.
void
pylith::faults::_FaultPoroDiffusionCohesiveKin::setKernelsLHSResidual(pylith::feassemble::IntegratorInterface *integrator,
                                                                      const pylith::faults::FaultPoroDiffusionCohesiveKin &fault,
                                                                      const pylith::topology::Field &solution,
                                                                      const pylith::problems::Physics::FormulationEnum formulation) {
    PYLITH_METHOD_BEGIN;
    pythia::journal::debug_t debug(_FaultPoroDiffusionCohesiveKin::pyreComponent);
    debug << pythia::journal::at(__HERE__)
          << "setKernelsLHSResidual(integrator=" << integrator << ", fault=" << typeid(fault).name()
          << ", solution=" << solution.getLabel() << ")"
          << pythia::journal::endl;

    std::vector<ResidualKernels> kernels(5);

    // Current solution field is
    // [displacement, pressure, trace_strain, lagrange, fault_pressure]

    // Elasticity equation (displacement/velocity).
    const PetscBdPointFunc f0u = pylith::fekernels::FaultPoroDiffusionCohesiveKin::f0u;
    const PetscBdPointFunc f1u = NULL;

    /* ** TO DO **
     * The pressure constraint and trace_strain constraint maybe unnecessary
     * since both f0 and f1 are NULLs, but I'm not sure so I kept it like this
     * ----Shengduo
     */

    // Pressure constraint
    const PetscBdPointFunc f0p = pylith::fekernels::FaultPoroDiffusionCohesiveKin::f0p;
    const PetscBdPointFunc f1p = NULL;

    // Trace_strain constraint
    const PetscBdPointFunc f0e = NULL;
    const PetscBdPointFunc f1e = NULL;

    // Fault_pressure constraint
    const PetscBdPointFunc f0p_fault = pylith::fekernels::FaultPoroDiffusionCohesiveKin::f0p_fault;
    const PetscBdPointFunc f1p_fault = pylith::fekernels::FaultPoroDiffusionCohesiveKin::f1p_fault;

    switch (formulation) {
    case pylith::problems::Physics::QUASISTATIC:
    {
        // Fault slip constraint equation.
        const PetscBdPointFunc f0l = pylith::fekernels::FaultPoroDiffusionCohesiveKin::f0l_u;
        const PetscBdPointFunc f1l = NULL;

        kernels[0] = ResidualKernels("displacement", f0u, f1u);
        kernels[1] = ResidualKernels("pressure", f0p, f1p);
        kernels[2] = ResidualKernels("trace_strain", f0e, f1e);
        kernels[3] = ResidualKernels("lagrange_multiplier_fault", f0l, f1l);
        kernels[4] = ResidualKernels("fault_pressure", f0p_fault, f1p_fault);

        break;
    } // QUASISTATIC
    case pylith::problems::Physics::DYNAMIC_IMEX:
    {
        // Fault slip rate constraint equation.
        const PetscBdPointFunc f0l = pylith::fekernels::FaultPoroDiffusionCohesiveKin::f0l_v;
        const PetscBdPointFunc f1l = NULL;

        kernels[0] = ResidualKernels("velocity", f0u, f1u);
        kernels[1] = ResidualKernels("pressure", f0p, f1p);
        kernels[2] = ResidualKernels("trace_strain", f0e, f1e);
        kernels[3] = ResidualKernels("lagrange_multiplier_fault", f0l, f1l);
        kernels[4] = ResidualKernels("fault_pressure", f0p_fault, f1p_fault);
        break;

    } // DYNAMIC_IMEX
    case pylith::problems::Physics::DYNAMIC:
    {
        pythia::journal::firewall_t firewall(_FaultPoroDiffusionCohesiveKin::pyreComponent);
        firewall << pythia::journal::at(__HERE__)
                 << "Poroelastic Fault implementation is incompatible with 'dynamic' formulation. Use 'dynamic_imex'."
                 << pythia::journal::endl;

    } // DYNAMIC
    default:
    {
        pythia::journal::firewall_t firewall(_FaultPoroDiffusionCohesiveKin::pyreComponent);
        firewall << pythia::journal::at(__HERE__)
                 << "Unknown formulation for equations (" << formulation << ")." << pythia::journal::endl;
    } // default
    } // switch

    assert(integrator);
    integrator->setKernelsLHSResidual(kernels);

    PYLITH_METHOD_END;
} // setKernelsLHSResidual


// ---------------------------------------------------------------------------------------------------------------------
// Set kernels for LHS Jacobian.
void
pylith::faults::_FaultPoroDiffusionCohesiveKin::setKernelsLHSJacobian(pylith::feassemble::IntegratorInterface *integrator,
                                                                      const pylith::faults::FaultPoroDiffusionCohesiveKin &fault,
                                                                      const pylith::topology::Field &solution,
                                                                      const pylith::problems::Physics::FormulationEnum formulation) {
    PYLITH_METHOD_BEGIN;
    pythia::journal::debug_t debug(_FaultPoroDiffusionCohesiveKin::pyreComponent);
    debug << pythia::journal::at(__HERE__)
          << "setKernelsLHSJacobian(integrator=" << integrator << ", fault=" << typeid(fault).name()
          << ", solution=" << solution.getLabel() << ")" << pythia::journal::endl;

    // Eight non-trivial kernels,
    // Jful; Jfp_fp_f, Jfp_fe, Jfp_fu, Jfp_fp; Jfpp, Jfpp_f; Jflu
    std::vector<JacobianKernels> kernels(8);

    const PetscBdPointJac Jf0ul = pylith::fekernels::FaultPoroDiffusionCohesiveKin::Jf0ul;
    const PetscBdPointJac Jf1ul = NULL;
    const PetscBdPointJac Jf2ul = NULL;
    const PetscBdPointJac Jf3ul = NULL;

    const PetscBdPointJac Jf0p_fp_f = pylith::fekernels::FaultPoroDiffusionCohesiveKin::Jf0p_fp_f;
    const PetscBdPointJac Jf1p_fp_f = NULL;
    const PetscBdPointJac Jf2p_fp_f = NULL;
    const PetscBdPointJac Jf3p_fp_f = pylith::fekernels::FaultPoroDiffusionCohesiveKin::Jf3p_fp_f;

    const PetscBdPointJac Jf0p_fe = pylith::fekernels::FaultPoroDiffusionCohesiveKin::Jf0p_fe;
    const PetscBdPointJac Jf1p_fe = NULL;
    const PetscBdPointJac Jf2p_fe = NULL;
    const PetscBdPointJac Jf3p_fe = NULL;

    const PetscBdPointJac Jf0p_fu = NULL;
    const PetscBdPointJac Jf1p_fu = pylith::fekernels::FaultPoroDiffusionCohesiveKin::Jf1p_fu;
    const PetscBdPointJac Jf2p_fu = NULL;
    const PetscBdPointJac Jf3p_fu = NULL;

    const PetscBdPointJac Jf0p_fp = pylith::fekernels::FaultPoroDiffusionCohesiveKin::Jf0p_fp;
    const PetscBdPointJac Jf1p_fp = NULL;
    const PetscBdPointJac Jf2p_fp = NULL;
    const PetscBdPointJac Jf3p_fp = pylith::fekernels::FaultPoroDiffusionCohesiveKin::Jf3p_fp;

    const PetscBdPointJac Jf0pp = pylith::fekernels::FaultPoroDiffusionCohesiveKin::Jf0pp;
    const PetscBdPointJac Jf1pp = NULL;
    const PetscBdPointJac Jf2pp = NULL;
    const PetscBdPointJac Jf3pp = NULL;

    const PetscBdPointJac Jf0pp_f = pylith::fekernels::FaultPoroDiffusionCohesiveKin::Jf0pp_f;
    const PetscBdPointJac Jf1pp_f = NULL;
    const PetscBdPointJac Jf2pp_f = NULL;
    const PetscBdPointJac Jf3pp_f = NULL;

    const PetscBdPointJac Jf0lu = pylith::fekernels::FaultPoroDiffusionCohesiveKin::Jf0lu;
    const PetscBdPointJac Jf1lu = NULL;
    const PetscBdPointJac Jf2lu = NULL;
    const PetscBdPointJac Jf3lu = NULL;

    const char *nameDispVel = NULL;
    const char *nameLagrangeMultiplier = "lagrange_multiplier_fault";
    const char *namePressure = "pressure";
    const char *nameTraceStrain = "trace_strain";
    const char *nameFaultPressure = "fault_pressure";

    switch (formulation) {
    case pylith::problems::Physics::QUASISTATIC:
        nameDispVel = "displacement";
        break;
    case pylith::problems::Physics::DYNAMIC_IMEX:
        nameDispVel = "velocity";
        break;
    case pylith::problems::Physics::DYNAMIC:
    {
        pythia::journal::firewall_t firewall(_FaultPoroDiffusionCohesiveKin::pyreComponent);
        firewall << pythia::journal::at(__HERE__)
                 << "Fault implementation is incompatible with 'dynamic' formulation. Use 'dynamic_imex'."
                 << pythia::journal::endl;
    } // DYNAMIC
    default:
    {
        pythia::journal::firewall_t firewall(_FaultPoroDiffusionCohesiveKin::pyreComponent);
        firewall << pythia::journal::at(__HERE__)
                 << "Unknown formulation for equations (" << formulation << ")."
                 << pythia::journal::endl;
    } // default
    } // switch

    kernels[0] = JacobianKernels(nameDispVel, nameLagrangeMultiplier, Jf0ul, Jf1ul, Jf2ul, Jf3ul);
    kernels[1] = JacobianKernels(nameFaultPressure, nameFaultPressure, Jf0p_fp_f, Jf1p_fp_f, Jf2p_fp_f, Jf3p_fp_f);
    kernels[2] = JacobianKernels(nameFaultPressure, nameTraceStrain, Jf0p_fe, Jf1p_fe, Jf2p_fe, Jf3p_fe);
    kernels[3] = JacobianKernels(nameFaultPressure, nameDispVel, Jf0p_fu, Jf1p_fu, Jf2p_fu, Jf3p_fu);
    kernels[4] = JacobianKernels(nameLagrangeMultiplier, nameDispVel, Jf0lu, Jf1lu, Jf2lu, Jf3lu);
    kernels[5] = JacobianKernels(namePressure, namePressure, Jf0pp, Jf1pp, Jf2pp, Jf3pp);
    kernels[6] = JacobianKernels(nameFaultPressure, namePressure, Jf0p_fp, Jf1p_fp, Jf2p_fp, Jf3p_fp);
    kernels[7] = JacobianKernels(namePressure, nameFaultPressure, Jf0pp_f, Jf1pp_f, Jf2pp_f, Jf3pp_f);

    assert(integrator);
    integrator->setKernelsLHSJacobian(kernels);

    PYLITH_METHOD_END;
} // setKernelsLHSJacobian


// End of file
