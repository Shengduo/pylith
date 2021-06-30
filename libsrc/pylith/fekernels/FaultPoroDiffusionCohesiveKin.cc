/* -*- C++ -*-
 *
 * ----------------------------------------------------------------------
 *
 * Brad T. Aagaard, U.S. Geological Survey
 * Charles A. Williams, GNS Science
 * Matthew G. Knepley, University of Chicago
 *
 * This code was developed as part of the Computational Infrastructure
 * for Geodynamics (http:*geodynamics.org).
 *
 * Copyright (c) 2010-2015 University of California, Davis
 *
 * See COPYING for license information.
 *
 * ----------------------------------------------------------------------
 */

#include <portinfo>

#include "pylith/fekernels/FaultPoroDiffusionCohesiveKin.hh"

#include <cassert> // USES assert()

/* ======================================================================
 * Kernels for prescribed fault slip.
 *
 * Solution fields: [disp(dim), vel(dim, if dynamic), pressure(1), trace_strain(1), lagrange(dim), fault_pressure(1)]
 *
 * Auxiliary fields:
 * -- numA : number of auxiliary fields
 *
 ***** Required fields
 * thickness(1),                       0
 * porosity(1),                        1
 * beta_p(1),                          2
 * beta_sigma(1),                      3
 * permeability_tangential(1),         4
 * permeability_normal(1),             5
 * fluid_viscosity(1),                 6
 * bulk_modulus_negative(1),           7
 * shear_modulus_negative(1),          8
 * bulk_modulus_positive(1),           9
 * shear_modulus_positive(1),          10
 * body_force(dim),                    numA - 3
 * source (1),                         numA - 2
 * slip(dim)                           numA - 1 
 * ** TO DO **
 * (**NEED IMPLEMENTATION**) *
 * ======================================================================
 */

namespace pylith
{
    namespace fekernels
    {
        class _FaultPoroDiffusionCohesiveKin
        {
        public:
            /** Get offset in s where velocity subfield starts.
             *
             * Normally this would be sOff, but sOff doesn't account for having DOF for the two sides of the fault
             * passed to the hybrid kernels. This functions computes the correct offset into s for the velocity
             * subfield.
             *
             * @param[in] sOff Offset of registered subfields in solution field [numS].
             * @param[in] numS Number of registered subfields in solution field.
             *
             * @returns Offset of velocity subfield in s.
             */
            static PylithInt velocity_sOff(const PylithInt sOff[],
                                           const PylithInt numS);

            /** Get offset in s where pressure subfield starts.
             *
             * Normally this would be sOff, but sOff doesn't account for having DOF for the two sides of the fault
             * passed to the hybrid kernels. This functions computes the correct offset into s for the velocity
             * subfield.
             *
             * @param[in] sOff Offset of registered subfields in solution field [numS].
             * @param[in] numS Number of registered subfields in solution field.
             *
             * @returns Offset of pressure subfield in s.
             */
            static PylithInt pressure_sOff(const PylithInt sOff[],
                                           const PylithInt numS);

            /** Get offset in s where trace_strain subfield starts.
             *
             * Normally this would be sOff, but sOff doesn't account for having DOF for the two sides of the fault
             * passed to the hybrid kernels. This functions computes the correct offset into s for the velocity
             * subfield.
             *
             * @param[in] sOff Offset of registered subfields in solution field [numS].
             * @param[in] numS Number of registered subfields in solution field.
             *
             * @returns Offset of trace_strain subfield in s.
             */
            static PylithInt trace_strain_sOff(const PylithInt sOff[],
                                               const PylithInt numS);

            /** Get offset in s where Lagrange multiplier subfield starts.
             *
             * Normally this would be sOff, but sOff doesn't account for having DOF for the two sides of the fault
             * passed to the hybrid kernels. This functions computes the correct offset into s for the Lagrange
             * multiplier subfield.
             *
             * @param[in] sOff Offset of registered subfields in solution field [numS].
             * @param[in] numS Number of registered subfields in solution field.
             *
             * @returns Offset of Lagrange multiplier subfield in s.
             */
            static PylithInt lagrange_sOff(const PylithInt sOff[],
                                           const PylithInt numS);

            /** Get offset in s where fault_pressure subfield starts.
             *
             * Normally this would be sOff, but sOff doesn't account for having DOF for the two sides of the fault
             * passed to the hybrid kernels. This functions computes the correct offset into s for the Lagrange
             * multiplier subfield.
             *
             * @param[in] sOff Offset of registered subfields in solution field [numS].
             * @param[in] numS Number of registered subfields in solution field.
             *
             * @returns Offset of fault_pressure subfield in s.
             */
            static PylithInt fault_pressure_sOff(const PylithInt sOff[],
                                                 const PylithInt numS);

            /* Compute tangential directions for 3-D fault.
             *
             * @param[in] dim Spatial dimension.
             * @param[in] refDir1 First choice for reference direction.
             * @param[in] refDir2 Second choice for reference direction if first fails.
             * @param[in] normDir Normal direction.
             * @param[out] tanDir1 First tangential direction.
             * @param[out] tanDIr2 Second tangential direction.
             */
            static void tangential_directions(const PylithInt dim,
                                              const PylithScalar refDir1[],
                                              const PylithScalar refDir2[],
                                              const PylithScalar normDir[],
                                              PylithScalar tanDir1[],
                                              PylithScalar tanDir2[]);

        }; // _FaultPoroDiffusionCohesiveKin
    }      // fekernels
} // pylith

// ----------------------------------------------------------------------
// Get offset in s where velocity subfield starts.
PylithInt
pylith::fekernels::_FaultPoroDiffusionCohesiveKin::velocity_sOff(const PylithInt sOff[],
                                                                 const PylithInt numS)
{
    PylithInt off = 0;
    const PylithInt numCount = 1; // [displacement, velocity, ...]
    for (PylithInt i = 0; i < numCount; ++i)
    {
        off += 2 * (sOff[i + 1] - sOff[i]);
    } // for
    return off;
} // velocity_sOff

// Get offset in s where pressure subfield starts.
PylithInt
pylith::fekernels::_FaultPoroDiffusionCohesiveKin::pressure_sOff(const PylithInt sOff[],
                                                                 const PylithInt numS)
{
    PylithInt off = 0;
    const PylithInt numCount = 2; // [displacement, velocity, pressure...]
    for (PylithInt i = 0; i < numCount; ++i)
    {
        off += 2 * (sOff[i + 1] - sOff[i]);
    } // for
    return off;
} // pressure_sOff

// Get offset in s where trace_strain subfield starts.
PylithInt
pylith::fekernels::_FaultPoroDiffusionCohesiveKin::trace_strain_sOff(const PylithInt sOff[],
                                                                     const PylithInt numS)
{
    PylithInt off = 0;
    const PylithInt numCount = 3; // [displacement, velocity, pressure, trace_strain...]
    for (PylithInt i = 0; i < numCount; ++i)
    {
        off += 2 * (sOff[i + 1] - sOff[i]);
    } // for
    return off;
} // trace_strain_sOff

// ----------------------------------------------------------------------
// Get offset in s where Lagrange multiplier field starts.
PylithInt
pylith::fekernels::_FaultPoroDiffusionCohesiveKin::lagrange_sOff(const PylithInt sOff[],
                                                                 const PylithInt numS)
{
    PylithInt off = 0;
    const PylithInt numCount = numS - 2; // Don't include last 2 field (Lagrange multiplier, fault_pressure)
    for (PylithInt i = 0; i < numCount; ++i)
    {
        off += 2 * (sOff[i + 1] - sOff[i]);
    } // for
    return off;
} // lagrange_sOff

// ----------------------------------------------------------------------
// Get offset in s where pressure_fault subfield starts.
// Seems that the offset contributed by Lagrange multiplier does not need *2
// ** TO DO **
// ** NEEDS VERFICATION **
PylithInt
pylith::fekernels::_FaultPoroDiffusionCohesiveKin::fault_pressure_sOff(const PylithInt sOff[],
                                                                       const PylithInt numS)
{
    PylithInt off = 0;
    const PylithInt numCount = numS - 2; // [..., Lagrange multiplier, fault_pressure]
    for (PylithInt i = 0; i < numCount; ++i)
    {
        off += 2 * (sOff[i + 1] - sOff[i]);
    } // for

    off += (sOff[numS - 1] - sOff[numS - 2]; 
    return off;
} // fault_pressure_sOff

// ----------------------------------------------------------------------
// Compute tangential directions from reference direction (first and second choice) and normal direction in 3-D.
void pylith::fekernels::_FaultPoroDiffuionCohesiveKin::tangential_directions(const PylithInt dim,
                                                                             const PylithScalar refDir1[],
                                                                             const PylithScalar refDir2[],
                                                                             const PylithScalar normDir[],
                                                                             PylithScalar tanDir1[],
                                                                             PylithScalar tanDir2[])
{
    assert(3 == dim);
    assert(refDir1);
    assert(refDir2);
    assert(normDir);
    assert(tanDir1);
    assert(tanDir2);

    const PylithInt _dim = 3;
    PylithScalar refDir[3] = {refDir1[0], refDir1[1], refDir1[2]};
    if (fabs(refDir[0] * normDir[0] + refDir[1] * normDir[1] + refDir[2] * normDir[2]) > 0.98)
    {
        for (PylithInt i = 0; i < _dim; ++i)
        {
            refDir[i] = refDir2[i];
        } // for
    }     // if

    // refDir x normDir, normalization required
    tanDir1[0] = +refDir[1] * normDir[2] - refDir[2] * normDir[1];
    tanDir1[1] = +refDir[2] * normDir[0] - refDir[0] * normDir[2];
    tanDir1[2] = +refDir[0] * normDir[1] - refDir[1] * normDir[0];

    PylithScalar _norm = sqrt(tanDir1[0] * tanDir1[0] + tanDir1[1] * tanDir1[1] + tanDir1[2] * tanDir1[2]);
    tanDir1[0] /= _norm;
    tanDir1[1] /= _norm;
    tanDir1[2] /= _norm;

    // normDir x tanDir1, normalization not required
    tanDir2[0] = +normDir[1] * tanDir1[2] - normDir[2] * tanDir1[1];
    tanDir2[1] = +normDir[2] * tanDir1[0] - normDir[0] * tanDir1[2];
    tanDir2[2] = +normDir[0] * tanDir1[1] - normDir[1] * tanDir1[0];
} // _tangential_directions

// ----------------------------------------------------------------------
// f0 function for elasticity equation: f0u = -\lambda (pos side), +\lambda (neg side).
void pylith::fekernels::FaultPoroDiffusionCohesiveKin::f0u(const PylithInt dim,
                                                           const PylithInt numS,
                                                           const PylithInt numA,
                                                           const PylithInt sOff[],
                                                           const PylithInt sOff_x[],
                                                           const PylithScalar s[],
                                                           const PylithScalar s_t[],
                                                           const PylithScalar s_x[],
                                                           const PylithInt aOff[],
                                                           const PylithInt aOff_x[],
                                                           const PylithScalar a[],
                                                           const PylithScalar a_t[],
                                                           const PylithScalar a_x[],
                                                           const PylithReal t,
                                                           const PylithScalar x[],
                                                           const PylithReal n[],
                                                           const PylithInt numConstants,
                                                           const PylithScalar constants[],
                                                           PylithScalar f0[])
{
    assert(sOff);
    assert(s);
    assert(f0);

    assert(numS >= 5);

    const PylithInt spaceDim = dim + 1; // :KLUDGE: dim passed in is spaceDim-1

    const PylithInt fOffN = 0;
    const PylithInt fOffP = fOffN + spaceDim;
    const PylithInt sOffLagrange = pylith::fekernels::_FaultPoroDiffusionCohesiveKin::lagrange_sOff(sOff, numS);
    const PylithScalar *lagrange = &s[sOffLagrange];

    for (PylithInt i = 0; i < spaceDim; ++i)
    {
        f0[fOffN + i] += -lagrange[i];
        f0[fOffP + i] += +lagrange[i];
    } // for
} // f0u

// ----------------------------------------------------------------------
// f0 function for bulk pressure: f0p = [\kappa_{cz} / \mu * ((p^+ - p^f)/h - n \cdot f_f),
//                                       \kappa_{cz} / \mu * ((p^- - p^f)/h + n \cdot f_f)]
void pylith::fekernels::FaultPoroDiffusionCohesiveKin::f0p(const PylithInt dim,
                                                           const PylithInt numS,
                                                           const PylithInt numA,
                                                           const PylithInt sOff[],
                                                           const PylithInt sOff_x[],
                                                           const PylithScalar s[],
                                                           const PylithScalar s_t[],
                                                           const PylithScalar s_x[],
                                                           const PylithInt aOff[],
                                                           const PylithInt aOff_x[],
                                                           const PylithScalar a[],
                                                           const PylithScalar a_t[],
                                                           const PylithScalar a_x[],
                                                           const PylithReal t,
                                                           const PylithScalar x[],
                                                           const PylithReal n[],
                                                           const PylithInt numConstants,
                                                           const PylithScalar constants[],
                                                           PylithScalar f0[])
{
    assert(sOff);
    assert(s);
    assert(f0);

    assert(numS >= 5);
    assert(a);
    assert(numA >= 10);

    const PylithInt spaceDim = dim + 1; // :KLUDGE: dim passed in is spaceDim-1
    const PylithInt i_thickness = 0;
    const PylithInt i_permeability_normal = 5;
    const PylithInt i_fluid_viscosity = 6;
    const PylithInt sOffPressureN = _FaultPoroDiffusionCohesiveKin::pressure_sOff(sOff, numS);
    const PylithInt sOffPressureP = sOffPressureN + 1;
    const PylithInt sOffPressureFault = _FaultPoroDiffusionCohesiveKin::fault_pressure_sOff(sOff, numS);
    const PylithInt fOffN = 0;
    const PylithInt fOffP = fOffN + 1;

    const PylithScalar thickness = a[aOff[i_thickness]];
    const PylithScalar permeabilityNormal = a[aOff[i_permeability_normal]];
    const PylithScalar fluidViscosity = a[aOff[i_fluid_viscosity]];
    const PylithScalar pressureN = s[sOffPressureN];
    const PylithScalar pressureP = s[sOffPressureP];
    const PylithScalar pressureFault = s[sOffPressureFault];

    f0[fOffN] += permeabilityNormal / fluidViscosity *
                 ((pressureN - pressureFault) / thickness);
    f0[fOffP] += permeabilityNormal / fluidViscosity *
                 ((pressureP - pressureFault) / thickness);
} // f0p

// ----------------------------------------------------------------------
// f0 function for bulk pressure: f0p = [\kappa_{cz} / \mu * ((p^+ - p^f)/h - n \cdot f_f),
//                                       \kappa_{cz} / \mu * ((p^- - p^f)/h + n \cdot f_f)]
void pylith::fekernels::FaultPoroDiffusionCohesiveKin::f0p_body(const PylithInt dim,
                                                           const PylithInt numS,
                                                           const PylithInt numA,
                                                           const PylithInt sOff[],
                                                           const PylithInt sOff_x[],
                                                           const PylithScalar s[],
                                                           const PylithScalar s_t[],
                                                           const PylithScalar s_x[],
                                                           const PylithInt aOff[],
                                                           const PylithInt aOff_x[],
                                                           const PylithScalar a[],
                                                           const PylithScalar a_t[],
                                                           const PylithScalar a_x[],
                                                           const PylithReal t,
                                                           const PylithScalar x[],
                                                           const PylithReal n[],
                                                           const PylithInt numConstants,
                                                           const PylithScalar constants[],
                                                           PylithScalar f0[])
{
    assert(sOff);
    assert(s);
    assert(f0);

    assert(numS >= 5);
    assert(a);
    assert(numA >= 10);

    const PylithInt spaceDim = dim + 1; // :KLUDGE: dim passed in is spaceDim-1
    const PylithInt i_thickness = 0;
    const PylithInt i_permeability_normal = 5;
    const PylithInt i_fluid_viscosity = 6;
    const PylithInt i_body_force = numA - 3;
    const PylithInt sOffPressureN = _FaultPoroDiffusionCohesiveKin::pressure_sOff(sOff, numS);
    const PylithInt sOffPressureP = sOffPressureN + 1;
    const PylithInt sOffPressureFault = _FaultPoroDiffusionCohesiveKin::fault_pressure_sOff(sOff, numS);
    const PylithInt fOffN = 0;
    const PylithInt fOffP = fOffN + 1;

    const PylithScalar thickness = a[aOff[i_thickness]];
    const PylithScalar permeabilityNormal = a[aOff[i_permeability_normal]];
    const PylithScalar fluidViscosity = a[aOff[i_fluid_viscosity]];
    const PylithScalar *bodyForce = &a[aOff[i_body_force]];
    const PylithScalar pressureN = s[sOffPressureN];
    const PylithScalar pressureP = s[sOffPressureP];
    const PylithScalar pressureFault = s[sOffPressureFault];

    PylithScalar nDotBodyForce = 0.;
    for (PylithInt i = 0; i < spaceDim; ++i)
    {
        nDotBodyForce += n[i] * bodyForce[i];
    }

    f0[fOffN] += permeabilityNormal / fluidViscosity *
                 ((pressureN - pressureFault) / thickness + nDotBodyForce);
    f0[fOffP] += permeabilityNormal / fluidViscosity *
                 ((pressureP - pressureFault) / thickness - nDotBodyForce);
} // f0p_body

// ----------------------------------------------------------------------
// f0 function for slip constraint equation: f0\lambda = (u^+ - u^-) - d
void pylith::fekernels::FaultPoroDiffusionCohesiveKin::f0l_u(const PylithInt dim,
                                                             const PylithInt numS,
                                                             const PylithInt numA,
                                                             const PylithInt sOff[],
                                                             const PylithInt sOff_x[],
                                                             const PylithScalar s[],
                                                             const PylithScalar s_t[],
                                                             const PylithScalar s_x[],
                                                             const PylithInt aOff[],
                                                             const PylithInt aOff_x[],
                                                             const PylithScalar a[],
                                                             const PylithScalar a_t[],
                                                             const PylithScalar a_x[],
                                                             const PylithReal t,
                                                             const PylithScalar x[],
                                                             const PylithReal n[],
                                                             const PylithInt numConstants,
                                                             const PylithScalar constants[],
                                                             PylithScalar f0[])
{
    assert(sOff);
    assert(aOff);
    assert(s);
    assert(a);
    assert(f0);

    assert(numS >= 5);
    assert(numA >= 6);

    const PylithInt spaceDim = dim + 1; // :KLUDGE: dim passed in is spaceDim-1
    const PylithInt i_slip = numA - 1;
    const PylithInt i_disp = 0;

    const PylithScalar *slip = &a[aOff[i_slip]];

    const PylithInt sOffDispN = sOff[i_disp];
    const PylithInt sOffDispP = sOffDispN + spaceDim;
    const PylithInt fOffLagrange = 0;

    const PylithScalar *dispN = &s[sOffDispN];
    const PylithScalar *dispP = &s[sOffDispP];

    switch (spaceDim)
    {
    case 2:
    {
        const PylithInt _spaceDim = 2;
        const PylithScalar tanDir[2] = {-n[1], n[0]};
        for (PylithInt i = 0; i < _spaceDim; ++i)
        {
            const PylithScalar slipXY = n[i] * slip[0] + tanDir[i] * slip[1];
            f0[fOffLagrange + i] += dispP[i] - dispN[i] - slipXY;
        } // for
        break;
    } // case 2
    case 3:
    {
        const PylithInt _spaceDim = 3;
        const PylithScalar *refDir1 = &constants[0];
        const PylithScalar *refDir2 = &constants[3];
        PylithScalar tanDir1[3], tanDir2[3];
        pylith::fekernels::_FaultPoroDiffusionCohesiveKin::tangential_directions(_spaceDim, refDir1, refDir2, n, tanDir1, tanDir2);

        for (PylithInt i = 0; i < _spaceDim; ++i)
        {
            const PylithScalar slipXYZ = n[i] * slip[0] + tanDir1[i] * slip[1] + tanDir2[i] * slip[2];
            f0[fOffLagrange + i] += dispP[i] - dispN[i] - slipXYZ;
        } // for
        break;
    } // case 3
    default:
        assert(0);
    } // switch
} // f0l_u

// ----------------------------------------------------------------------
// f0 function for slip rate constraint equation: f0\lambda = (v^+ - v^-) - \dot{d}
void pylith::fekernels::FaultPoroDiffusionCohesiveKin::f0l_v(const PylithInt dim,
                                                             const PylithInt numS,
                                                             const PylithInt numA,
                                                             const PylithInt sOff[],
                                                             const PylithInt sOff_x[],
                                                             const PylithScalar s[],
                                                             const PylithScalar s_t[],
                                                             const PylithScalar s_x[],
                                                             const PylithInt aOff[],
                                                             const PylithInt aOff_x[],
                                                             const PylithScalar a[],
                                                             const PylithScalar a_t[],
                                                             const PylithScalar a_x[],
                                                             const PylithReal t,
                                                             const PylithScalar x[],
                                                             const PylithReal n[],
                                                             const PylithInt numConstants,
                                                             const PylithScalar constants[],
                                                             PylithScalar f0[])
{
    assert(sOff);
    assert(aOff);
    assert(s);
    assert(a);
    assert(f0);

    assert(numS >= 5);
    assert(numA >= 6);

    const PylithInt spaceDim = dim + 1; // :KLUDGE: dim passed in is spaceDim-1
    const PylithInt i_slipRate = numA - 1;

    const PylithScalar *slipRate = &a[aOff[i_slipRate]];

    const PylithInt sOffVelN = _FaultPoroDiffusionCohesiveKin::velocity_sOff(sOff, numS);
    const PylithInt sOffVelP = sOffVelN + spaceDim;
    const PylithInt fOffLagrange = 0;

    const PylithScalar *velN = &s[sOffVelN];
    const PylithScalar *velP = &s[sOffVelP];

    switch (spaceDim)
    {
    case 2:
    {
        const PylithInt _spaceDim = 2;
        const PylithScalar tanDir[2] = {-n[1], n[0]};
        for (PylithInt i = 0; i < _spaceDim; ++i)
        {
            const PylithScalar slipRateXY = n[i] * slipRate[0] + tanDir[i] * slipRate[1];
            f0[fOffLagrange + i] += velP[i] - velN[i] - slipRateXY;
        } // for
        break;
    } // case 2
    case 3:
    {
        const PylithInt _spaceDim = 3;
        const PylithScalar *refDir1 = &constants[0];
        const PylithScalar *refDir2 = &constants[3];
        PylithScalar tanDir1[3], tanDir2[3];
        pylith::fekernels::_FaultPoroDiffusionCohesiveKin::tangential_directions(_spaceDim, refDir1, refDir2, n, tanDir1, tanDir2);

        for (PylithInt i = 0; i < _spaceDim; ++i)
        {
            const PylithScalar slipRateXYZ = n[i] * slipRate[0] + tanDir1[i] * slipRate[1] + tanDir2[i] * slipRate[2];
            f0[fOffLagrange + i] += velP[i] - velN[i] - slipRateXYZ;
        } // for
        break;
    } // case 3
    default:
        assert(0);
    } // switch
} // f0l_v

// ----------------------------------------------------------------------
// f0 function for p_fault constraint equation:
// f0p_fault = porosity * (\beta^p * (\dot{p}^+ + 2 \dot{p}^f + \dot{p}^-)/4
//                         + \beta^\sigma * (\dot{\sigma}^+_{nn} + \dot{\sigma}^-_{nn}))
//             + \kappa_{fx} / \mu * \vnabla(2D) \cdot body_force
//             - \kappa_{fz} / \mu * (p^+ - 2p_f + p^-) / h^2 - source
void pylith::fekernels::FaultPoroDiffusionCohesiveKin::f0p_fault(const PylithInt dim,
                                                                 const PylithInt numS,
                                                                 const PylithInt numA,
                                                                 const PylithInt sOff[],
                                                                 const PylithInt sOff_x[],
                                                                 const PylithScalar s[],
                                                                 const PylithScalar s_t[],
                                                                 const PylithScalar s_x[],
                                                                 const PylithInt aOff[],
                                                                 const PylithInt aOff_x[],
                                                                 const PylithScalar a[],
                                                                 const PylithScalar a_t[],
                                                                 const PylithScalar a_x[],
                                                                 const PylithReal t,
                                                                 const PylithScalar x[],
                                                                 const PylithReal n[],
                                                                 const PylithInt numConstants,
                                                                 const PylithScalar constants[],
                                                                 PylithScalar f0[])
{
    assert(sOff);
    assert(aOff);
    assert(sOff_x);
    assert(a_x);
    assert(aOff_x);
    assert(s);
    assert(a);
    assert(s_x);
    assert(s_t);
    assert(f0);

    assert(numS >= 5);
    assert(numA >= 10);

    const PylithInt spaceDim = dim + 1; // :KLUDGE: dim passed in is spaceDim-1
    const PylithInt i_disp_x = 0;

    // Index for auxiliary fields
    const PylithInt i_thickness = 0;
    const PylithInt i_porosity = 1;
    const PylithInt i_beta_p = 2;
    const PylithInt i_beta_sigma = 3;
    const PylithInt i_permeabilility_tangential = 4;
    const PylithInt i_permeabilility_normal = 5;
    const PylithInt i_fluid_viscosity = 6;
    const PylithInt i_bulk_modulus_negative = 7;
    const PylithInt i_shear_modulus_negative = 8;
    const PylithInt i_bulk_modulus_positive = 9;
    const PylithInt i_shear_modulus_positive = 10;

    const PylithScalar thickness = a[aOff[i_thickness]];
    const PylithScalar porosity = a[aOff[i_porosity]];
    const PylithScalar betaP = a[aOff[i_beta_p]];
    const PylithScalar betaSigma = a[aOff[i_beta_sigma]];
    const PylithScalar permeabilityTangential = a[aOff[i_permeabilility_tangential]];
    const PylithScalar permeabilityNormal = a[aOff[i_permeabilility_normal]];
    const PylithScalar fluidViscosity = a[aOff[i_fluid_viscosity]];

    // Pressure and pressure_t
    const PylithInt sOffpressureN = _FaultPoroDiffusionCohesiveKin::pressure_sOff(sOff, numS);
    const PylithInt sOffpressureP = sOffpressureN + 1;
    const PylithInt sOffpressureFault = _FaultPoroDiffusionCohesiveKin::fault_pressure_sOff(sOff, numS);

    const PylithScalar pressureN = s[sOffpressureN];
    const PylithScalar pressureP = s[sOffpressureP];
    const PylithScalar pressureFault = s[sOffpressureFault];
    const PylithScalar pressureN_t = s_t[sOffpressureN];
    const PylithScalar pressureP_t = s_t[sOffpressureP];
    const PylithScalar pressureFault_t = s_t[sOffpressureFault];

    // ** TO DO **
    // Pull out drained bulk modulus and shear modulus from the surrounding bulks
    const PylithScalar bulkModulusN = a[aOff[i_bulk_modulus_negative]];
    const PylithScalar bulkModulusP = a[aOff[i_bulk_modulus_positive]];
    const PylithScalar shearModulusN = a[aOff[i_shear_modulus_negative]];
    const PylithScalar shearModulusP = a[aOff[i_shear_modulus_positive]];

    // Strain components
    const PylithInt sOffDispN_x = sOff_x[i_disp_x];
    const PylithInt sOffDispP_x = sOffDispN_x + spaceDim ^ 2;
    const PylithInt fOffp_fault = 0;

    const PylithScalar *dispN_x = &s_x[sOffDispN_x];
    const PylithScalar *dispP_x = &s_x[sOffDispP_x];

    // Trace_strain, no transformation required

    const PylithInt sOffTraceStrainN = _FaultPoroDiffusionCohesiveKin::trace_strain_sOff(sOff, numS);
    const PylithInt sOffTraceStrainP = sOffTraceStrainN + 1;
    const PylithScalar traceStrainN = s[sOffTraceStrainN];
    const PylithScalar traceStrainP = s[sOffTraceStrainP];

    // \sigma_nn, requires transformation
    // \sigma_nn = (K_u - 2G/3) \epsilon_v + 2G strain_nn
    // strain_nn = n_i u_{i,j} n_j
    PylithScalar strain_nnN = 0.;
    PylithScalar strain_nnP = 0.;
    for (PylithInt i = 0; i < spaceDim; ++i)
    {
        for (PylithInt j = 0; j < spaceDim; ++j)
        {
            strain_nnN += n[i] * dispN_x[i * spaceDim + j] * n[j];
            strain_nnP += n[i] * dispP_x[i * spaceDim + j] * n[j];
        }
    }

    const PylithScalar stress_nnN = (bulkModulusN - 2. * shearModulusN / 3.) * traceStrainN + 2. * shearModulusN * strain_nnN;
    const PylithScalar stress_nnP = (bulkModulusP - 2. * shearModulusP / 3.) * traceStrainP + 2. * shearModulusP * strain_nnP;

    f0[fOffp_fault] += porosity * (betaP * (pressureN_t + 2. * pressureFault_t + pressureP_t) / 4. +
                                   betaSigma * (stress_nnN + stress_nnP) / 2.) +
                           permeabilityTangential / fluidViscosity - permeabilityNormal / fluidViscosity * (pressureP - 2. * pressureFault + pressureN) / (thickness*thickness);

} // f0p_fault

// ----------------------------------------------------------------------
// f0 function for p_fault constraint equation:
// f0p_fault = porosity * (\beta^p * (\dot{p}^+ + 2 \dot{p}^f + \dot{p}^-)/4
//                         + \beta^\sigma * (\dot{\sigma}^+_{nn} + \dot{\sigma}^-_{nn}))
//             + \kappa_{fx} / \mu * \vnabla(2D) \cdot body_force
//             - \kappa_{fz} / \mu * (p^+ - 2p_f + p^-) / h^2 - source
void pylith::fekernels::FaultPoroDiffusionCohesiveKin::f0p_fault_body(const PylithInt dim,
                                                                 const PylithInt numS,
                                                                 const PylithInt numA,
                                                                 const PylithInt sOff[],
                                                                 const PylithInt sOff_x[],
                                                                 const PylithScalar s[],
                                                                 const PylithScalar s_t[],
                                                                 const PylithScalar s_x[],
                                                                 const PylithInt aOff[],
                                                                 const PylithInt aOff_x[],
                                                                 const PylithScalar a[],
                                                                 const PylithScalar a_t[],
                                                                 const PylithScalar a_x[],
                                                                 const PylithReal t,
                                                                 const PylithScalar x[],
                                                                 const PylithReal n[],
                                                                 const PylithInt numConstants,
                                                                 const PylithScalar constants[],
                                                                 PylithScalar f0[])
{
    assert(sOff);
    assert(aOff);
    assert(sOff_x);
    assert(a_x);
    assert(aOff_x);
    assert(s);
    assert(a);
    assert(s_x);
    assert(s_t);
    assert(f0);

    assert(numS >= 5);
    assert(numA >= 10);

    const PylithInt spaceDim = dim + 1; // :KLUDGE: dim passed in is spaceDim-1
    const PylithInt i_disp_x = 0;

    // Index for auxiliary fields
    const PylithInt i_thickness = 0;
    const PylithInt i_porosity = 1;
    const PylithInt i_beta_p = 2;
    const PylithInt i_beta_sigma = 3;
    const PylithInt i_permeabilility_tangential = 4;
    const PylithInt i_permeabilility_normal = 5;
    const PylithInt i_fluid_viscosity = 6;
    const PylithInt i_bulk_modulus_negative = 7;
    const PylithInt i_shear_modulus_negative = 8;
    const PylithInt i_bulk_modulus_positive = 9;
    const PylithInt i_shear_modulus_positive = 10;
    const PylithInt i_body_force = numA - 3;

    const PylithScalar thickness = a[aOff[i_thickness]];
    const PylithScalar porosity = a[aOff[i_porosity]];
    const PylithScalar betaP = a[aOff[i_beta_p]];
    const PylithScalar betaSigma = a[aOff[i_beta_sigma]];
    const PylithScalar permeabilityTangential = a[aOff[i_permeabilility_tangential]];
    const PylithScalar permeabilityNormal = a[aOff[i_permeabilility_normal]];
    const PylithScalar fluidViscosity = a[aOff[i_fluid_viscosity]];
    const PylithScalar *bodyForce_x = &a_x[aOff_x[i_body_force]];

    // Divergence of body force in the first spaceDim - 1 directions
    PylithScalar bodyForce_div = 0.;

    for (PylithInt i = 0; i < spaceDim - 1; ++i)
    {
        bodyForce_div += bodyForce_x[i * spaceDim + i];
    }

    // Pressure and pressure_t
    const PylithInt sOffpressureN = _FaultPoroDiffusionCohesiveKin::pressure_sOff(sOff, numS);
    const PylithInt sOffpressureP = sOffpressureN + 1;
    const PylithInt sOffpressureFault = _FaultPoroDiffusionCohesiveKin::fault_pressure_sOff(sOff, numS);

    const PylithScalar pressureN = s[sOffpressureN];
    const PylithScalar pressureP = s[sOffpressureP];
    const PylithScalar pressureFault = s[sOffpressureFault];
    const PylithScalar pressureN_t = s_t[sOffpressureN];
    const PylithScalar pressureP_t = s_t[sOffpressureP];
    const PylithScalar pressureFault_t = s_t[sOffpressureFault];

    // ** TO DO **
    // Pull out drained bulk modulus and shear modulus from the surrounding bulks
    const PylithScalar bulkModulusN = a[aOff[i_bulk_modulus_negative]];
    const PylithScalar bulkModulusP = a[aOff[i_bulk_modulus_positive]];
    const PylithScalar shearModulusN = a[aOff[i_shear_modulus_negative]];
    const PylithScalar shearModulusP = a[aOff[i_shear_modulus_positive]];

    // Strain components
    const PylithInt sOffDispN_x = sOff_x[i_disp_x];
    const PylithInt sOffDispP_x = sOffDispN_x + spaceDim ^ 2;
    const PylithInt fOffp_fault = 0;

    const PylithScalar *dispN_x = &s_x[sOffDispN_x];
    const PylithScalar *dispP_x = &s_x[sOffDispP_x];

    // Trace_strain, no transformation required

    const PylithInt sOffTraceStrainN = _FaultPoroDiffusionCohesiveKin::trace_strain_sOff(sOff, numS);
    const PylithInt sOffTraceStrainP = sOffTraceStrainN + 1;
    const PylithScalar traceStrainN = s[sOffTraceStrainN];
    const PylithScalar traceStrainP = s[sOffTraceStrainP];

    // \sigma_nn, requires transformation
    // \sigma_nn = (K_u - 2G/3) \epsilon_v + 2G strain_nn
    // strain_nn = n_i u_{i,j} n_j
    PylithScalar strain_nnN = 0.;
    PylithScalar strain_nnP = 0.;
    for (PylithInt i = 0; i < spaceDim; ++i)
    {
        for (PylithInt j = 0; j < spaceDim; ++j)
        {
            strain_nnN += n[i] * dispN_x[i * spaceDim + j] * n[j];
            strain_nnP += n[i] * dispP_x[i * spaceDim + j] * n[j];
        }
    }

    const PylithScalar stress_nnN = (bulkModulusN - 2. * shearModulusN / 3.) * traceStrainN + 2. * shearModulusN * strain_nnN;
    const PylithScalar stress_nnP = (bulkModulusP - 2. * shearModulusP / 3.) * traceStrainP + 2. * shearModulusP * strain_nnP;

    f0[fOffp_fault] += porosity * (betaP * (pressureN_t + 2. * pressureFault_t + pressureP_t) / 4. +
                                   betaSigma * (stress_nnN + stress_nnP) / 2.) +
                           permeabilityTangential / fluidViscosity * bodyForce_div - permeabilityNormal / fluidViscosity * (pressureP - 2. * pressureFault + pressureN) / (thickness*thickness);

} // f0p_fault_body_source


// ----------------------------------------------------------------------
// f0 function for p_fault constraint equation:
// f0p_fault = porosity * (\beta^p * (\dot{p}^+ + 2 \dot{p}^f + \dot{p}^-)/4
//                         + \beta^\sigma * (\dot{\sigma}^+_{nn} + \dot{\sigma}^-_{nn}))
//             + \kappa_{fx} / \mu * \vnabla(2D) \cdot body_force
//             - \kappa_{fz} / \mu * (p^+ - 2p_f + p^-) / h^2 - source
void pylith::fekernels::FaultPoroDiffusionCohesiveKin::f0p_fault_source(const PylithInt dim,
                                                                 const PylithInt numS,
                                                                 const PylithInt numA,
                                                                 const PylithInt sOff[],
                                                                 const PylithInt sOff_x[],
                                                                 const PylithScalar s[],
                                                                 const PylithScalar s_t[],
                                                                 const PylithScalar s_x[],
                                                                 const PylithInt aOff[],
                                                                 const PylithInt aOff_x[],
                                                                 const PylithScalar a[],
                                                                 const PylithScalar a_t[],
                                                                 const PylithScalar a_x[],
                                                                 const PylithReal t,
                                                                 const PylithScalar x[],
                                                                 const PylithReal n[],
                                                                 const PylithInt numConstants,
                                                                 const PylithScalar constants[],
                                                                 PylithScalar f0[])
{
    assert(sOff);
    assert(aOff);
    assert(sOff_x);
    assert(a_x);
    assert(aOff_x);
    assert(s);
    assert(a);
    assert(s_x);
    assert(s_t);
    assert(f0);

    assert(numS >= 5);
    assert(numA >= 10);

    const PylithInt spaceDim = dim + 1; // :KLUDGE: dim passed in is spaceDim-1
    const PylithInt i_disp_x = 0;

    // Index for auxiliary fields
    const PylithInt i_thickness = 0;
    const PylithInt i_porosity = 1;
    const PylithInt i_beta_p = 2;
    const PylithInt i_beta_sigma = 3;
    const PylithInt i_permeabilility_tangential = 4;
    const PylithInt i_permeabilility_normal = 5;
    const PylithInt i_fluid_viscosity = 6;
    const PylithInt i_bulk_modulus_negative = 7;
    const PylithInt i_shear_modulus_negative = 8;
    const PylithInt i_bulk_modulus_positive = 9;
    const PylithInt i_shear_modulus_positive = 10;
    const PylithInt i_source = numA - 2;

    const PylithScalar thickness = a[aOff[i_thickness]];
    const PylithScalar porosity = a[aOff[i_porosity]];
    const PylithScalar betaP = a[aOff[i_beta_p]];
    const PylithScalar betaSigma = a[aOff[i_beta_sigma]];
    const PylithScalar permeabilityTangential = a[aOff[i_permeabilility_tangential]];
    const PylithScalar permeabilityNormal = a[aOff[i_permeabilility_normal]];
    const PylithScalar fluidViscosity = a[aOff[i_fluid_viscosity]];
    const PylithScalar *bodyForce_x = &a_x[aOff_x[i_body_force]];
    const PylithScalar source = a[aOff[i_source]];


    // Pressure and pressure_t
    const PylithInt sOffpressureN = _FaultPoroDiffusionCohesiveKin::pressure_sOff(sOff, numS);
    const PylithInt sOffpressureP = sOffpressureN + 1;
    const PylithInt sOffpressureFault = _FaultPoroDiffusionCohesiveKin::fault_pressure_sOff(sOff, numS);

    const PylithScalar pressureN = s[sOffpressureN];
    const PylithScalar pressureP = s[sOffpressureP];
    const PylithScalar pressureFault = s[sOffpressureFault];
    const PylithScalar pressureN_t = s_t[sOffpressureN];
    const PylithScalar pressureP_t = s_t[sOffpressureP];
    const PylithScalar pressureFault_t = s_t[sOffpressureFault];

    // ** TO DO **
    // Pull out drained bulk modulus and shear modulus from the surrounding bulks
    const PylithScalar bulkModulusN = a[aOff[i_bulk_modulus_negative]];
    const PylithScalar bulkModulusP = a[aOff[i_bulk_modulus_positive]];
    const PylithScalar shearModulusN = a[aOff[i_shear_modulus_negative]];
    const PylithScalar shearModulusP = a[aOff[i_shear_modulus_positive]];

    // Strain components
    const PylithInt sOffDispN_x = sOff_x[i_disp_x];
    const PylithInt sOffDispP_x = sOffDispN_x + spaceDim ^ 2;
    const PylithInt fOffp_fault = 0;

    const PylithScalar *dispN_x = &s_x[sOffDispN_x];
    const PylithScalar *dispP_x = &s_x[sOffDispP_x];

    // Trace_strain, no transformation required

    const PylithInt sOffTraceStrainN = _FaultPoroDiffusionCohesiveKin::trace_strain_sOff(sOff, numS);
    const PylithInt sOffTraceStrainP = sOffTraceStrainN + 1;
    const PylithScalar traceStrainN = s[sOffTraceStrainN];
    const PylithScalar traceStrainP = s[sOffTraceStrainP];

    // \sigma_nn, requires transformation
    // \sigma_nn = (K_u - 2G/3) \epsilon_v + 2G strain_nn
    // strain_nn = n_i u_{i,j} n_j
    PylithScalar strain_nnN = 0.;
    PylithScalar strain_nnP = 0.;
    for (PylithInt i = 0; i < spaceDim; ++i)
    {
        for (PylithInt j = 0; j < spaceDim; ++j)
        {
            strain_nnN += n[i] * dispN_x[i * spaceDim + j] * n[j];
            strain_nnP += n[i] * dispP_x[i * spaceDim + j] * n[j];
        }
    }

    const PylithScalar stress_nnN = (bulkModulusN - 2. * shearModulusN / 3.) * traceStrainN + 2. * shearModulusN * strain_nnN;
    const PylithScalar stress_nnP = (bulkModulusP - 2. * shearModulusP / 3.) * traceStrainP + 2. * shearModulusP * strain_nnP;

    f0[fOffp_fault] += porosity * (betaP * (pressureN_t + 2. * pressureFault_t + pressureP_t) / 4. +
                                   betaSigma * (stress_nnN + stress_nnP) / 2.) +
                           permeabilityTangential / fluidViscosity  - permeabilityNormal / fluidViscosity * (pressureP - 2. * pressureFault + pressureN) / thickness ^
                       2 - source;

} // f0p_fault_source

// ----------------------------------------------------------------------
// f0 function for p_fault constraint equation:
// f0p_fault = porosity * (\beta^p * (\dot{p}^+ + 2 \dot{p}^f + \dot{p}^-)/4
//                         + \beta^\sigma * (\dot{\sigma}^+_{nn} + \dot{\sigma}^-_{nn}))
//             + \kappa_{fx} / \mu * \vnabla(2D) \cdot body_force
//             - \kappa_{fz} / \mu * (p^+ - 2p_f + p^-) / h^2 - source
void pylith::fekernels::FaultPoroDiffusionCohesiveKin::f0p_fault_body_source(const PylithInt dim,
                                                                 const PylithInt numS,
                                                                 const PylithInt numA,
                                                                 const PylithInt sOff[],
                                                                 const PylithInt sOff_x[],
                                                                 const PylithScalar s[],
                                                                 const PylithScalar s_t[],
                                                                 const PylithScalar s_x[],
                                                                 const PylithInt aOff[],
                                                                 const PylithInt aOff_x[],
                                                                 const PylithScalar a[],
                                                                 const PylithScalar a_t[],
                                                                 const PylithScalar a_x[],
                                                                 const PylithReal t,
                                                                 const PylithScalar x[],
                                                                 const PylithReal n[],
                                                                 const PylithInt numConstants,
                                                                 const PylithScalar constants[],
                                                                 PylithScalar f0[])
{
    assert(sOff);
    assert(aOff);
    assert(sOff_x);
    assert(a_x);
    assert(aOff_x);
    assert(s);
    assert(a);
    assert(s_x);
    assert(s_t);
    assert(f0);

    assert(numS >= 5);
    assert(numA >= 10);

    const PylithInt spaceDim = dim + 1; // :KLUDGE: dim passed in is spaceDim-1
    const PylithInt i_disp_x = 0;

    // Index for auxiliary fields
    const PylithInt i_thickness = 0;
    const PylithInt i_porosity = 1;
    const PylithInt i_beta_p = 2;
    const PylithInt i_beta_sigma = 3;
    const PylithInt i_permeabilility_tangential = 4;
    const PylithInt i_permeabilility_normal = 5;
    const PylithInt i_fluid_viscosity = 6;
    const PylithInt i_bulk_modulus_negative = 7;
    const PylithInt i_shear_modulus_negative = 8;
    const PylithInt i_bulk_modulus_positive = 9;
    const PylithInt i_shear_modulus_positive = 10;
    const PylithInt i_body_force = numA - 3;
    const PylithInt i_source = numA - 2;

    const PylithScalar thickness = a[aOff[i_thickness]];
    const PylithScalar porosity = a[aOff[i_porosity]];
    const PylithScalar betaP = a[aOff[i_beta_p]];
    const PylithScalar betaSigma = a[aOff[i_beta_sigma]];
    const PylithScalar permeabilityTangential = a[aOff[i_permeabilility_tangential]];
    const PylithScalar permeabilityNormal = a[aOff[i_permeabilility_normal]];
    const PylithScalar fluidViscosity = a[aOff[i_fluid_viscosity]];
    const PylithScalar *bodyForce_x = &a_x[aOff_x[i_body_force]];
    const PylithScalar source = a[aOff[i_source]];

    // Divergence of body force in the first spaceDim - 1 directions
    PylithScalar bodyForce_div = 0.;

    for (PylithInt i = 0; i < spaceDim - 1; ++i)
    {
        bodyForce_div += bodyForce_x[i * spaceDim + i];
    }

    // Pressure and pressure_t
    const PylithInt sOffpressureN = _FaultPoroDiffusionCohesiveKin::pressure_sOff(sOff, numS);
    const PylithInt sOffpressureP = sOffpressureN + 1;
    const PylithInt sOffpressureFault = _FaultPoroDiffusionCohesiveKin::fault_pressure_sOff(sOff, numS);

    const PylithScalar pressureN = s[sOffpressureN];
    const PylithScalar pressureP = s[sOffpressureP];
    const PylithScalar pressureFault = s[sOffpressureFault];
    const PylithScalar pressureN_t = s_t[sOffpressureN];
    const PylithScalar pressureP_t = s_t[sOffpressureP];
    const PylithScalar pressureFault_t = s_t[sOffpressureFault];

    // ** TO DO **
    // Pull out drained bulk modulus and shear modulus from the surrounding bulks
    const PylithScalar bulkModulusN = a[aOff[i_bulk_modulus_negative]];
    const PylithScalar bulkModulusP = a[aOff[i_bulk_modulus_positive]];
    const PylithScalar shearModulusN = a[aOff[i_shear_modulus_negative]];
    const PylithScalar shearModulusP = a[aOff[i_shear_modulus_positive]];

    // Strain components
    const PylithInt sOffDispN_x = sOff_x[i_disp_x];
    const PylithInt sOffDispP_x = sOffDispN_x + spaceDim ^ 2;
    const PylithInt fOffp_fault = 0;

    const PylithScalar *dispN_x = &s_x[sOffDispN_x];
    const PylithScalar *dispP_x = &s_x[sOffDispP_x];

    // Trace_strain, no transformation required

    const PylithInt sOffTraceStrainN = _FaultPoroDiffusionCohesiveKin::trace_strain_sOff(sOff, numS);
    const PylithInt sOffTraceStrainP = sOffTraceStrainN + 1;
    const PylithScalar traceStrainN = s[sOffTraceStrainN];
    const PylithScalar traceStrainP = s[sOffTraceStrainP];

    // \sigma_nn, requires transformation
    // \sigma_nn = (K_u - 2G/3) \epsilon_v + 2G strain_nn
    // strain_nn = n_i u_{i,j} n_j
    PylithScalar strain_nnN = 0.;
    PylithScalar strain_nnP = 0.;
    for (PylithInt i = 0; i < spaceDim; ++i)
    {
        for (PylithInt j = 0; j < spaceDim; ++j)
        {
            strain_nnN += n[i] * dispN_x[i * spaceDim + j] * n[j];
            strain_nnP += n[i] * dispP_x[i * spaceDim + j] * n[j];
        }
    }

    const PylithScalar stress_nnN = (bulkModulusN - 2. * shearModulusN / 3.) * traceStrainN + 2. * shearModulusN * strain_nnN;
    const PylithScalar stress_nnP = (bulkModulusP - 2. * shearModulusP / 3.) * traceStrainP + 2. * shearModulusP * strain_nnP;

    f0[fOffp_fault] += porosity * (betaP * (pressureN_t + 2. * pressureFault_t + pressureP_t) / 4. +
                                   betaSigma * (stress_nnN + stress_nnP) / 2.) +
                           permeabilityTangential / fluidViscosity * bodyForce_div - permeabilityNormal / fluidViscosity * (pressureP - 2. * pressureFault + pressureN) / thickness ^
                       2 - source;

} // f0p_fault_body_source

// ----------------------------------------------------------------------
// f1 function for p_fault constraint equation:
// f1p_fault = \kappa_{fx} / (4\mu) \vnabla (p^+ + 2 p^f + p^-)
void pylith::fekernels::FaultPoroDiffusionCohesiveKin::f1p_fault(const PylithInt dim,
                                                                 const PylithInt numS,
                                                                 const PylithInt numA,
                                                                 const PylithInt sOff[],
                                                                 const PylithInt sOff_x[],
                                                                 const PylithScalar s[],
                                                                 const PylithScalar s_t[],
                                                                 const PylithScalar s_x[],
                                                                 const PylithInt aOff[],
                                                                 const PylithInt aOff_x[],
                                                                 const PylithScalar a[],
                                                                 const PylithScalar a_t[],
                                                                 const PylithScalar a_x[],
                                                                 const PylithReal t,
                                                                 const PylithScalar x[],
                                                                 const PylithReal n[],
                                                                 const PylithInt numConstants,
                                                                 const PylithScalar constants[],
                                                                 PylithScalar f1[])
{
    assert(aOff);
    assert(sOff_x);
    assert(s);
    assert(a);
    assert(s_x);
    assert(f1);

    assert(numS >= 5);
    assert(numA >= 10);

    const PylithInt spaceDim = dim + 1; // :KLUDGE: dim passed in is spaceDim-1

    // Index for auxiliary fields
    const PylithInt i_permeabilility_tangential = 4;
    const PylithInt i_fluid_viscosity = 6;

    const PylithScalar permeabilityTangential = a[aOff[i_permeabilility_tangential]];
    const PylithScalar fluidViscosity = a[aOff[i_fluid_viscosity]];

    // Pressure_x
    const PylithInt sOffPressureN = _FaultPoroDiffusionCohesiveKin::pressure_sOff(sOff, numS);
    const PylithInt sOffPressureP = sOffPressureN + 1;
    const PylithInt sOffPressureFault = _FaultPoroDiffusionCohesiveKin::fault_pressure_sOff(sOff, numS);

    const PylithScalar *pressureN_x = &s_x[sOff_x[sOffPressureN]];
    const PylithScalar *pressureP_x = &s_x[sOff_x[sOffPressureP]];
    const PylithScalar *pressureFault_x = &s_x[sOff_x[sOffPressureFault]];
    const PylithInt fOffp_fault = 0;

    for (PylithInt i = 0; i < spaceDim; ++i)
    {
        f1[fOffp_fault + i] += permeabilityTangential / fluidViscosity / 4. *
                               (pressureN_x[i] + 2. * pressureFault_x[i] + pressureP_x[i]);
    }
} // f1p_fault

// ----------------------------------------------------------------------
/* Jf0 function for integration of the displacement equation.
 *
 * Solution fields = [disp(dim), ..., lagrange(dim)]
 * Auxiliary fields = None
 */
void pylith::fekernels::FaultPoroDiffusionCohesiveKin::Jf0ul(const PylithInt dim,
                                                             const PylithInt numS,
                                                             const PylithInt numA,
                                                             const PylithInt sOff[],
                                                             const PylithInt sOff_x[],
                                                             const PylithScalar s[],
                                                             const PylithScalar s_t[],
                                                             const PylithScalar s_x[],
                                                             const PylithInt aOff[],
                                                             const PylithInt aOff_x[],
                                                             const PylithScalar a[],
                                                             const PylithScalar a_t[],
                                                             const PylithScalar a_x[],
                                                             const PylithReal t,
                                                             const PylithReal s_tshift,
                                                             const PylithScalar x[],
                                                             const PylithReal n[],
                                                             const PylithInt numConstants,
                                                             const PylithScalar constants[],
                                                             PylithScalar Jf0[])
{
    assert(numS >= 5);
    assert(numA >= 10);
    assert(Jf0);
    assert(sOff);
    assert(aOff);
    assert(n);

    const PylithInt spaceDim = dim + 1; // :KLUDGE: dim passed in is spaceDim-1

    const PylithInt gOffN = 0;
    const PylithInt gOffP = gOffN + spaceDim;
    const PylithInt ncols = spaceDim;

    for (PylithInt i = 0; i < spaceDim; ++i)
    {
        Jf0[(gOffN + i) * ncols + i] += -1.0;
        Jf0[(gOffP + i) * ncols + i] += +1.0;
    } // for
} // Jg0ul

// ----------------------------------------------------------------------
/* Jf0p_fp_f function for integration of the displacement equation:
 * 2 \kappa_{fz} / (\mu h^2) + 2 \phi_f \beta^p t_shift
 * Solution fields = [disp(dim), ..., lagrange(dim), fault_pressure(1)]
 * Auxiliary fields
 */
void pylith::fekernels::FaultPoroDiffusionCohesiveKin::Jf0p_fp_f(const PylithInt dim,
                                                                 const PylithInt numS,
                                                                 const PylithInt numA,
                                                                 const PylithInt sOff[],
                                                                 const PylithInt sOff_x[],
                                                                 const PylithScalar s[],
                                                                 const PylithScalar s_t[],
                                                                 const PylithScalar s_x[],
                                                                 const PylithInt aOff[],
                                                                 const PylithInt aOff_x[],
                                                                 const PylithScalar a[],
                                                                 const PylithScalar a_t[],
                                                                 const PylithScalar a_x[],
                                                                 const PylithReal t,
                                                                 const PylithReal s_tshift,
                                                                 const PylithScalar x[],
                                                                 const PylithReal n[],
                                                                 const PylithInt numConstants,
                                                                 const PylithScalar constants[],
                                                                 PylithScalar Jf0[])
{
    // Check data fields
    assert(numS >= 5);
    assert(numA >= 10);
    assert(Jf0);
    assert(sOff);
    assert(aOff);
    assert(n);
    const PylithInt gOff = 0;

    // Index for auxiliary fields
    const PylithInt i_thickness = 0;
    const PylithInt i_porosity = 1;
    const PylithInt i_beta_p = 2;
    const PylithInt i_permeabilility_normal = 5;
    const PylithInt i_fluid_viscosity = 6;

    const PylithScalar thickness = a[aOff[i_thickness]];
    const PylithScalar porosity = a[aOff[i_porosity]];
    const PylithScalar betaP = a[aOff[i_beta_p]];
    const PylithScalar permeabilityNormal = a[aOff[i_permeabilility_normal]];
    const PylithScalar fluidViscosity = a[aOff[i_fluid_viscosity]];

    Jf0[gOff] += 2 permeabilityNormal / (fluidViscosity * thickness ^ 2) + 2. * porosity * betaP * s_t;

} // Jf0p_fp_f

// ----------------------------------------------------------------------
/* Jf3p_fp_f function for integration of the displacement equation:
 * \kappa_{fx} / (2 \mu) \te{I}
 * Solution fields = [disp(dim), ..., lagrange(dim), fault_pressure(1)]
 * Auxiliary fields
 */
void pylith::fekernels::FaultPoroDiffusionCohesiveKin::Jf3p_fp_f(const PylithInt dim,
                                                                 const PylithInt numS,
                                                                 const PylithInt numA,
                                                                 const PylithInt sOff[],
                                                                 const PylithInt sOff_x[],
                                                                 const PylithScalar s[],
                                                                 const PylithScalar s_t[],
                                                                 const PylithScalar s_x[],
                                                                 const PylithInt aOff[],
                                                                 const PylithInt aOff_x[],
                                                                 const PylithScalar a[],
                                                                 const PylithScalar a_t[],
                                                                 const PylithScalar a_x[],
                                                                 const PylithReal t,
                                                                 const PylithReal s_tshift,
                                                                 const PylithScalar x[],
                                                                 const PylithReal n[],
                                                                 const PylithInt numConstants,
                                                                 const PylithScalar constants[],
                                                                 PylithScalar Jf3[])
{
    // Check data fields
    assert(numS >= 5);
    assert(numA >= 10);
    assert(Jf3);
    assert(sOff);
    assert(aOff);
    assert(n);
    const PylithInt gOff = 0;

    const PylithInt spaceDim = dim + 1;

    // Index for auxiliary fields
    const PylithInt i_permeabilility_tangential = 4;
    const PylithInt i_fluid_viscosity = 6;

    const PylithScalar permeabilityTangential = a[aOff[i_permeabilility_tangential]];
    const PylithScalar fluidViscosity = a[aOff[i_fluid_viscosity]];

    const PylithInt ncols = spaceDim;
    for (PylithInt i = 0; i < spaceDim; ++i)
    {
        Jf3[i * ncols + i] += permeabilityTangential / 2. / fluidViscosity;
    }

} // Jf3p_fp_f

// ----------------------------------------------------------------------
/* Jf1p_fu function for integration of the displacement equation.
 * [t_shift \phi_f \beta^\sigma G \ve{n} \ve{n}, t_shift \phi_f \beta^\sigma G \ve{n} \ve{n}]
 * Solution fields = [disp(dim), ..., lagrange(dim), fault_pressure(1)]
 * Auxiliary fields
 */
void pylith::fekernels::FaultPoroDiffusionCohesiveKin::Jf1p_fu(const PylithInt dim,
                                                               const PylithInt numS,
                                                               const PylithInt numA,
                                                               const PylithInt sOff[],
                                                               const PylithInt sOff_x[],
                                                               const PylithScalar s[],
                                                               const PylithScalar s_t[],
                                                               const PylithScalar s_x[],
                                                               const PylithInt aOff[],
                                                               const PylithInt aOff_x[],
                                                               const PylithScalar a[],
                                                               const PylithScalar a_t[],
                                                               const PylithScalar a_x[],
                                                               const PylithReal t,
                                                               const PylithReal s_tshift,
                                                               const PylithScalar x[],
                                                               const PylithReal n[],
                                                               const PylithInt numConstants,
                                                               const PylithScalar constants[],
                                                               PylithScalar Jf1[])
{
    // Check data fields
    assert(numS >= 5);
    assert(numA >= 10);
    assert(Jf1);
    assert(sOff);
    assert(aOff);
    assert(n);

    const PylithInt spaceDim = dim + 1;
    const PylithInt gOffN = 0;
    const PylithInt gOffP = gOffN + spaceDim;

    // Index for auxiliary fields
    const PylithInt i_porosity = 1;
    const PylithInt i_beta_sigma = 3;
    const PylithInt i_shear_modulus_negative = 8;
    const PylithInt i_shear_modulus_positive = 10;

    // ** TO DO **
    // Pull out shear modulus from the bulks
    const PylithScalar shearModulusN = a[aOff[i_shear_modulus_negative]];
    const PylithScalar shearModulusP = a[aOff[i_shear_modulus_positive]];

    const PylithScalar porosity = a[aOff[i_porosity]];
    const PylithScalar betaSigma = a[aOff[i_beta_sigma]];

    const PylithInt ncols = 2 * spaceDim;
    for (PylithInt i = 0; i < spaceDim; ++i)
    {
        for (PylithInt j = 0; j < spaceDim; ++j)
        {
            Jf1[i * ncols + gOffN + j] += s_tshift * porosity * betaSigma * shearModulusN * n[i] * n[j];
            Jf1[i * ncols + gOffP + j] += s_tshift * porosity * betaSigma * shearModulusP * n[i] * n[j];
        }
    }
} // Jf1p_fu

// ----------------------------------------------------------------------
/* Jf0p_fe function for integration of the displacement equation.
 *
 * Solution fields = [disp(dim), ..., lagrange(dim), fault_pressure(1)]
 * Auxiliary fields
 */
void pylith::fekernels::FaultPoroDiffusionCohesiveKin::Jf0p_fe(const PylithInt dim,
                                                               const PylithInt numS,
                                                               const PylithInt numA,
                                                               const PylithInt sOff[],
                                                               const PylithInt sOff_x[],
                                                               const PylithScalar s[],
                                                               const PylithScalar s_t[],
                                                               const PylithScalar s_x[],
                                                               const PylithInt aOff[],
                                                               const PylithInt aOff_x[],
                                                               const PylithScalar a[],
                                                               const PylithScalar a_t[],
                                                               const PylithScalar a_x[],
                                                               const PylithReal t,
                                                               const PylithReal s_tshift,
                                                               const PylithScalar x[],
                                                               const PylithReal n[],
                                                               const PylithInt numConstants,
                                                               const PylithScalar constants[],
                                                               PylithScalar Jf0[])
{
    assert(numS >= 5);
    assert(numA >= 10);
    assert(Jf0);
    assert(sOff);
    assert(aOff);
    assert(n);

    const PylithInt spaceDim = dim + 1;
    const PylithInt gOffN = 0;
    const PylithInt gOffP = gOffN + 1;

    // Index for auxiliary fields
    const PylithInt i_porosity = 1;
    const PylithInt i_beta_sigma = 3;

    // ** TO DO **
    // Pull out shear and bulk modulus from the bulks
    const PylithInt i_bulk_modulus_negative = 7;
    const PylithInt i_shear_modulus_negative = 8;
    const PylithInt i_bulk_modulus_positive = 9;
    const PylithInt i_shear_modulus_positive = 10;
    
    const PylithScalar bulkModulusN = a[aOff[i_bulk_modulus_negative]];
    const PylithScalar bulkModulusP = a[aOff[i_bulk_modulus_positive]];
    const PylithScalar shearModulusN = a[aOff[i_shear_modulus_negative]];
    const PylithScalar shearModulusP = a[aOff[i_shear_modulus_positive]];

    const PylithScalar porosity = a[aOff[i_porosity]];
    const PylithScalar betaSigma = a[aOff[i_beta_sigma]];

    Jf0[gOffN] += s_tshift * porosity * betaSigma * (bulkModulusN - 2. * shearModulusN / 3.) / 2.;
    Jf0[gOffP] += s_tshift * porosity * betaSigma * (bulkModulusP - 2. * shearModulusP / 3.) / 2.;

} // Jf0p_fe

// ----------------------------------------------------------------------
/* Jf0p_fp function for integration of the displacement equation.
 * [\phi_f beta^p t_shift /4 - \kappa_{fz}/\mu h^2,\phi_f beta^p t_shift /4 - \kappa_{fz}/\mu h^2]
 * Solution fields = [disp(dim), ..., lagrange(dim), fault_pressure(1)]
 * Auxiliary fields
 */
void pylith::fekernels::FaultPoroDiffusionCohesiveKin::Jf0p_fp(const PylithInt dim,
                                                               const PylithInt numS,
                                                               const PylithInt numA,
                                                               const PylithInt sOff[],
                                                               const PylithInt sOff_x[],
                                                               const PylithScalar s[],
                                                               const PylithScalar s_t[],
                                                               const PylithScalar s_x[],
                                                               const PylithInt aOff[],
                                                               const PylithInt aOff_x[],
                                                               const PylithScalar a[],
                                                               const PylithScalar a_t[],
                                                               const PylithScalar a_x[],
                                                               const PylithReal t,
                                                               const PylithReal s_tshift,
                                                               const PylithScalar x[],
                                                               const PylithReal n[],
                                                               const PylithInt numConstants,
                                                               const PylithScalar constants[],
                                                               PylithScalar Jf0[])
{
    assert(numS >= 5);
    assert(numA >= 10);
    assert(Jf0);
    assert(sOff);
    assert(aOff);
    assert(n);

    // Index for auxiliary fields
    const PylithInt i_thickness = 0;
    const PylithInt i_porosity = 1;
    const PylithInt i_beta_p = 2;
    const PylithInt i_permeabilility_normal = 5;
    const PylithInt i_fluid_viscosity = 6;

    const PylithScalar thickness = a[aOff[i_thickness]];
    const PylithScalar porosity = a[aOff[i_porosity]];
    const PylithScalar betaP = a[aOff[i_beta_p]];
    const PylithScalar permeabilityNormal = a[aOff[i_permeabilility_normal]];
    const PylithScalar fluidViscosity = a[aOff[i_fluid_viscosity]];

    const PylithInt spaceDim = dim + 1; // :KLUDGE: dim passed in is spaceDim-1

    const PylithInt gOffN = 0;
    const PylithInt gOffP = gOffN + 1;

    Jf0[gOffN] += s_tshift * porosity * betaP / 4. - permeabilityNormal / fluidViscosity / thickness ^ 2;
    Jf0[gOffP] += s_tshift * porosity * betaP / 4. - permeabilityNormal / fluidViscosity / thickness ^ 2;

} // Jf0p_fp

// ----------------------------------------------------------------------
/* Jf3p_fp function for integration of the displacement equation.
 * [\kappa_{fx} / 4 / \mu \te{I}, \kappa_{fx} / 4 / \mu \te{I}]
 * Solution fields = [disp(dim), ..., lagrange(dim), fault_pressure(1)]
 * Auxiliary fields
 */
void pylith::fekernels::FaultPoroDiffusionCohesiveKin::Jf3p_fp(const PylithInt dim,
                                                               const PylithInt numS,
                                                               const PylithInt numA,
                                                               const PylithInt sOff[],
                                                               const PylithInt sOff_x[],
                                                               const PylithScalar s[],
                                                               const PylithScalar s_t[],
                                                               const PylithScalar s_x[],
                                                               const PylithInt aOff[],
                                                               const PylithInt aOff_x[],
                                                               const PylithScalar a[],
                                                               const PylithScalar a_t[],
                                                               const PylithScalar a_x[],
                                                               const PylithReal t,
                                                               const PylithReal s_tshift,
                                                               const PylithScalar x[],
                                                               const PylithReal n[],
                                                               const PylithInt numConstants,
                                                               const PylithScalar constants[],
                                                               PylithScalar Jf3[])
{
    assert(numS >= 5);
    assert(numA >= 10);
    assert(Jf3);
    assert(sOff);
    assert(aOff);
    assert(n);

    // Index for auxiliary fields
    const PylithInt i_permeabilility_tangential = 4;
    const PylithInt i_fluid_viscosity = 6;

    const PylithScalar permeabilityTangential = a[aOff[i_permeabilility_tangential]];
    const PylithScalar fluidViscosity = a[aOff[i_fluid_viscosity]];

    const PylithInt spaceDim = dim + 1; // :KLUDGE: dim passed in is spaceDim-1

    const PylithInt gOffN = 0;
    const PylithInt gOffP = gOffN + spaceDim;

    const PylithInt ncols = 2 * spaceDim;

    for (PylithInt i = 0; i < spaceDim; ++i)
    {
        Jf3[i * ncols + gOffN + i] += permeabilityTangential / 4. / fluidViscosity;
        Jf3[i * ncols + gOffP + i] += permeabilityTangential / 4. / fluidViscosity;
    }
} // Jf3p_fp

// ----------------------------------------------------------------------
/* Jf0pp function for integration of the slip constraint equation.
 *
 * Solution fields = [disp(dim), ..., lagrange(dim)]
 * Auxiliary fields = None
 */
void pylith::fekernels::FaultPoroDiffusionCohesiveKin::Jf0pp(const PylithInt dim,
                                                             const PylithInt numS,
                                                             const PylithInt numA,
                                                             const PylithInt sOff[],
                                                             const PylithInt sOff_x[],
                                                             const PylithScalar s[],
                                                             const PylithScalar s_t[],
                                                             const PylithScalar s_x[],
                                                             const PylithInt aOff[],
                                                             const PylithInt aOff_x[],
                                                             const PylithScalar a[],
                                                             const PylithScalar a_t[],
                                                             const PylithScalar a_x[],
                                                             const PylithReal t,
                                                             const PylithReal s_tshift,
                                                             const PylithScalar x[],
                                                             const PylithReal n[],
                                                             const PylithInt numConstants,
                                                             const PylithScalar constants[],
                                                             PylithScalar Jf0[])
{
    assert(numS >= 5);
    assert(numA >= 10);
    assert(Jf0);
    assert(sOff);
    assert(aOff);
    assert(n);

    // Index for auxiliary fields
    const PylithInt i_thickness = 0;
    const PylithInt i_permeabilility_normal = 5;
    const PylithInt i_fluid_viscosity = 6;

    const PylithScalar thickness = a[aOff[i_thickness]];
    const PylithScalar permeabilityNormal = a[aOff[i_permeabilility_normal]];
    const PylithScalar fluidViscosity = a[aOff[i_fluid_viscosity]];

    const PylithInt gOffN = 0;
    const PylithInt gOffP = gOffN + 3;

    Jf0[gOffN] += permeabilityNormal / thickness / fluidViscosity;
    Jf0[gOffP] += permeabilityNormal / thickness / fluidViscosity;

} // Jf0pp

// ----------------------------------------------------------------------
/* Jf0pp_f function for integration of the slip constraint equation.
 * [-\kappa_cz / \mu / h, -\kappa_cz / \mu / h]
 * Solution fields = [disp(dim), ..., lagrange(dim)]
 * Auxiliary fields = None
 */
void pylith::fekernels::FaultPoroDiffusionCohesiveKin::Jf0pp_f(const PylithInt dim,
                                                               const PylithInt numS,
                                                               const PylithInt numA,
                                                               const PylithInt sOff[],
                                                               const PylithInt sOff_x[],
                                                               const PylithScalar s[],
                                                               const PylithScalar s_t[],
                                                               const PylithScalar s_x[],
                                                               const PylithInt aOff[],
                                                               const PylithInt aOff_x[],
                                                               const PylithScalar a[],
                                                               const PylithScalar a_t[],
                                                               const PylithScalar a_x[],
                                                               const PylithReal t,
                                                               const PylithReal s_tshift,
                                                               const PylithScalar x[],
                                                               const PylithReal n[],
                                                               const PylithInt numConstants,
                                                               const PylithScalar constants[],
                                                               PylithScalar Jf0[])
{
    assert(numS >= 5);
    assert(numA >= 10);
    assert(Jf0);
    assert(sOff);
    assert(aOff);
    assert(n);

    // Index for auxiliary fields
    const PylithInt i_thickness = 0;
    const PylithInt i_permeabilility_normal = 5;
    const PylithInt i_fluid_viscosity = 6;

    const PylithScalar thickness = a[aOff[i_thickness]];
    const PylithScalar permeabilityNormal = a[aOff[i_permeabilility_normal]];
    const PylithScalar fluidViscosity = a[aOff[i_fluid_viscosity]];

    const PylithInt gOffN = 0;
    const PylithInt gOffP = gOffN + 1;

    Jf0[gOffN] += -permeabilityNormal / thickness / fluidViscosity;
    Jf0[gOffP] += -permeabilityNormal / thickness / fluidViscosity;

} // Jf0pp_f

// ----------------------------------------------------------------------
/* Jg0 function for integration of the slip constraint equation.
 *
 * Solution fields = [disp(dim), ..., lagrange(dim)]
 * Auxiliary fields = None
 */
void pylith::fekernels::FaultPoroDiffusionCohesiveKin::Jf0lu(const PylithInt dim,
                                                             const PylithInt numS,
                                                             const PylithInt numA,
                                                             const PylithInt sOff[],
                                                             const PylithInt sOff_x[],
                                                             const PylithScalar s[],
                                                             const PylithScalar s_t[],
                                                             const PylithScalar s_x[],
                                                             const PylithInt aOff[],
                                                             const PylithInt aOff_x[],
                                                             const PylithScalar a[],
                                                             const PylithScalar a_t[],
                                                             const PylithScalar a_x[],
                                                             const PylithReal t,
                                                             const PylithReal s_tshift,
                                                             const PylithScalar x[],
                                                             const PylithReal n[],
                                                             const PylithInt numConstants,
                                                             const PylithScalar constants[],
                                                             PylithScalar Jf0[])
{
    assert(numS >= 5);
    assert(numA >= 10);
    assert(Jf0);
    assert(sOff);
    assert(aOff);
    assert(n);

    const PylithInt spaceDim = dim + 1; // :KLUDGE: dim passed in is spaceDim-1

    const PylithInt gOffN = 0;
    const PylithInt gOffP = gOffN + spaceDim;
    const PylithInt ncols = 2 * spaceDim;

    for (PylithInt i = 0; i < spaceDim; ++i)
    {
        Jf0[i * ncols + gOffN + i] += -1.0;
        Jf0[i * ncols + gOffP + i] += +1.0;
    } // for
} // Jg0lu

// End of file
