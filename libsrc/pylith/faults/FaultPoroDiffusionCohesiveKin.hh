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

/** @file libsrc/faults/FaultPoroDiffusionCohesiveKin.hh
 *
 * @brief C++ implementation for a undrained poroelastic fault surface with
 * fluid diffusion and kinematic
 * (prescribed) slip implemented with cohesive elements.
 */

#if !defined(pylith_faults_faultporodiffusioncohesivekin_hh)
#define pylith_faults_faultporodiffusioncohesivekin_hh

#include "FaultCohesive.hh" // ISA FaultCohesive

#include <string> // HASA std::string
#include <map> // HASA std::map

class pylith::faults::FaultPoroDiffusionCohesiveKin : public pylith::faults::FaultCohesive {
    // ** TO DO **
    // Implement this friend class
    friend class TestFaultPoroDiffusionCohesiveKin; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    FaultPoroDiffusionCohesiveKin(void);

    /// Destructor.
    ~FaultPoroDiffusionCohesiveKin(void);

    /// Deallocate PETSc and local data structures.
    void deallocate(void);

    /** Set kinematic earthquake ruptures.
     *
     * @param names Array of kinematic earthquake rupture names.
     * @param numNames Number of earthquake rupture names.
     * @param ruptures Array of kinematic earthquake ruptures.
     * @param numRuptures Number of earthquake ruptures.
     */
    void setEqRuptures(const char* const* names,
                       const int numNames,
                       KinSrc** ruptures,
                       const int numRuptures);

    /** Include body force?
     *
     * @param[in] value Flag indicating to include body force term.
     */
    void useBodyForce(const bool value);

    /** Include body force?
     *
     * @returns True if including body force term, false otherwise.
     */
    bool useBodyForce(void) const;

    /** Include source?
     *
     * @param[in] value Flag indicating to include source term.
     */
    void useSource(const bool value);

    /** Include source?
     *
     * @returns True if including source term, false otherwise.
     */
    bool useSource(void) const;

    /** Include constant pressure source?
     *
     * @param[in] value Flag indicating to include constant pressure source term.
     */
    void useConstantPressureSource(const bool value);

    /** Include constant pressure source?
     *
     * @returns True if including constant pressure source term, false otherwise.
     */
    bool useConstantPressureSource(void) const;

    /** Verify configuration is acceptable.
     *
     * @param[in] solution Solution field.
     */
    void verifyConfiguration(const pylith::topology::Field& solution) const;

    /** Create integrator and set kernels.
     *
     * @param[in] solution Solution field.
     * @returns Integrator if applicable, otherwise NULL.
     */
    pylith::feassemble::Integrator* createIntegrator(const pylith::topology::Field& solution);

    /** Create constraint and set kernels.
     *
     * @param[in] solution Solution field.
     * @returns Constraint if applicable, otherwise NULL.
     */
    pylith::feassemble::Constraint* createConstraint(const pylith::topology::Field& solution);

    /** Create auxiliary field.
     *
     * @param[in] solution Solution field.
     * @param[in\ domainMesh Finite-element mesh associated with integration domain.
     *
     * @returns Auxiliary field if applicable, otherwise NULL.
     */
    pylith::topology::Field* createAuxiliaryField(const pylith::topology::Field& solution,
                                                  const pylith::topology::Mesh& domainMesh);

    /** Create derived field.
     *
     * @param[in] solution Solution field.
     * @param[in\ domainMesh Finite-element mesh associated with integration domain.
     *
     * @returns Derived field if applicable, otherwise NULL.
     */
    pylith::topology::Field* createDerivedField(const pylith::topology::Field& solution,
                                                const pylith::topology::Mesh& domainMesh);

    /** Update auxiliary subfields at beginning of time step.
     *
     * @param[out] auxiliaryField Auxiliary field.
     * @param[in] t Current time.
     */
    void updateAuxiliaryField(pylith::topology::Field* auxiliaryField,
                              const double t);

    // PROTECTED METHODS ///////////////////////////////////////////////////////////////////////////////////////////////
protected:

    /** Get auxiliary factory associated with physics.
     *
     * @return Auxiliary factory for physics object.
     */
    pylith::feassemble::AuxiliaryFactory* _getAuxiliaryFactory(void);

    /** Update kernel constants.
     *
     * @param[in] dt Current time step.
     */
    void _updateKernelConstants(const PylithReal dt);

    /** Update slip subfield in auxiliary field at beginning of time step.
     *
     * @param[out] auxiliaryField Auxiliary field.
     * @param[in] t Current time.
     */
    void _updateSlip(pylith::topology::Field* auxiliaryField,
                     const double t);

    /** Update slip rate subfield in auxiliary field at beginning of time step.
     *
     * @param[out] auxiliaryField Auxiliary field.
     * @param[in] t Current time.
     */
    void _updateSlipRate(pylith::topology::Field* auxiliaryField,
                         const double t);

    // PRIVATE TYPEDEFS ////////////////////////////////////////////////////////////////////////////////////////////////
private:

    typedef std::map<std::string, KinSrc*> srcs_type;

    // PRIVATE MEMBERS /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    pylith::faults::AuxiliaryFactoryKinematic* _auxiliaryFactory; ///< Factory for auxiliary subfields.
    srcs_type _ruptures; ///< Array of kinematic earthquake ruptures.
    PetscVec _slipVecRupture; ///< PETSc local Vec to hold slip for one kinematic rupture.
    PetscVec _slipVecTotal; ///< PETSc local Vec to hold slip for all kinematic ruptures.

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    FaultPoroDiffusionCohesiveKin(const FaultPoroDiffusionCohesiveKin&); ///< Not implemented
    const FaultPoroDiffusionCohesiveKin& operator=(const FaultPoroDiffusionCohesiveKin&); ///< Not implemented.

    static PetscErrorCode _zero(PetscInt dim,
                                PetscReal t,
                                const PetscReal x[],
                                PetscInt Nc,
                                PetscScalar *u,
                                void *ctx) {
        for (int c = 0; c < Nc; ++c) {
            u[c] = 0.0;
        }
        return 0;
    }

}; // class FaultPoroDiffusionCohesiveKin

#endif // pylith_faults_faultporodiffusioncohesivekin_hh

// End of file
