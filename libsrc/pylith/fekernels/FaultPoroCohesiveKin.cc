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

#include "pylith/fekernels/FaultPoroCohesiveKin.hh"

#include <cassert> // USES assert()

/* ======================================================================
 * Kernels for prescribed fault slip.
 *
 * Solution fields: [disp(dim), vel(dim, if dynamic), pressure(1), trace_strain(1), lagrange(dim), fault_pressure(1)]
 * 
 * Auxiliary fields = [fault_undrained_bulk_modulus, 
 *                     fault_shear_modulus, 
 *                     fault_skempton_coefficient (B), 
 *                     solid_undrained_bulk_modulus,
 *                     solid_shear_modulus, 
 *                     slip (dim)]                        (**NEED IMPLEMENTATION)
 * 
 * -- numA : number of auxiliary fields
 ***** Required fields
 * - 0: fault_undrained_bulk_modulus(1)
 * - 1: fault_shear_modulus(1)
 * - 2: fault_skempton_coefficient[B](1) 
 * - 3: solid_undrained_bulk_modulus(1)
 * - 4: solid_shear_modulus(1)
 * - numA - 1: slip(dim)
 * 
 * ======================================================================
 */

namespace pylith
{
    namespace fekernels
    {
        class _FaultPoroCohesiveKin
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

        }; // _FaultCohesiveKin
    }      // fekernels
} // pylith

// ----------------------------------------------------------------------
// Get offset in s where velocity subfield starts.
PylithInt
pylith::fekernels::_FaultPoroCohesiveKin::velocity_sOff(const PylithInt sOff[],
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
pylith::fekernels::_FaultPoroCohesiveKin::pressure_sOff(const PylithInt sOff[],
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
pylith::fekernels::_FaultPoroCohesiveKin::trace_strain_sOff(const PylithInt sOff[],
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
pylith::fekernels::_FaultPoroCohesiveKin::lagrange_sOff(const PylithInt sOff[],
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

// Get offset in s where pressure_fault subfield starts.
PylithInt
pylith::fekernels::_FaultPoroCohesiveKin::fault_pressure_sOff(const PylithInt sOff[],
                                                              const PylithInt numS)
{
    PylithInt off = 0;
    const PylithInt numCount = numS - 1; // [..., fault_pressure]
    for (PylithInt i = 0; i < numCount; ++i)
    {
        off += 2 * (sOff[i + 1] - sOff[i]);
    } // for
    return off;
} // fault_pressure_sOff
// ----------------------------------------------------------------------
// Compute tangential directions from reference direction (first and second choice) and normal direction in 3-D.
void pylith::fekernels::_FaultPoroCohesiveKin::tangential_directions(const PylithInt dim,
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

    PylithScalar _norm = sqrt(tanDir1[0] ^ 2 + tanDir1[1] ^ 2 + tanDir1[2] ^ 2);
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
void pylith::fekernels::FaultPoroCohesiveKin::f0u(const PylithInt dim,
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
    const PylithInt sOffLagrange = pylith::fekernels::_FaultPoroCohesiveKin::lagrange_sOff(sOff, numS);
    const PylithScalar *lagrange = &s[sOffLagrange];

    for (PylithInt i = 0; i < spaceDim; ++i)
    {
        f0[fOffN + i] += -lagrange[i];
        f0[fOffP + i] += +lagrange[i];
    } // for
} // f0u

// ----------------------------------------------------------------------
// f0 function for slip constraint equation: f0\lambda = (u^+ - u^-) - d
void pylith::fekernels::FaultPoroCohesiveKin::f0l_u(const PylithInt dim,
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
        pylith::fekernels::_FaultPoroCohesiveKin::tangential_directions(_spaceDim, refDir1, refDir2, n, tanDir1, tanDir2);

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
// f0 function for p' constraint equation: f0p_fault = p' - B'K'_u / M_u'
// * [G' M_u / (G K_u) * 3K_u(trace_strain^- + trace_strain^-) / 2 +
// (G - G') / G * ((K_u - 2G/3) * (trace_strain^- + trace_strain^-) + 2G (u3,3^- + u3,3^+))/2]
void pylith::fekernels::FaultPoroCohesiveKin::f0p_fault(const PylithInt dim,
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
    assert(s);
    assert(a);
    assert(s_x);
    assert(f0);

    assert(numS >= 5);
    assert(numA >= 6);

    const PylithInt spaceDim = dim + 1; // :KLUDGE: dim passed in is spaceDim-1
    const PylithInt i_disp_x = 0;

    // Index for auxiliary fields
    const PylithInt i_fault_undrained_bulk_modulus = 0;
    const PylithInt i_fault_shear_modulus = 1;
    const PylithInt i_fault_skempton_coefficient = 2;
    const PylithInt i_solid_undrained_bulk_modulus = 3;
    const PylithInt i_solid_shear_modulus = 4;

    const PylithScalar faultUndrainedBulkModulus = a[aOff[i_fault_undrained_bulk_modulus]];
    const PylithScalar faultShearModulus = a[aOff[i_fault_shear_modulus]];
    const PylithScalar faultSkemptonCoef = a[aOff[i_fault_skempton_coefficient]];
    const PylithScalar undrainedBulkModulus = a[aOff[i_solid_undrained_bulk_modulus]];
    const PylithScalar shearModulus = a[aOff[i_solid_shear_modulus]];

    // Strain components
    const PylithInt sOffDispN_x = sOff_x[i_disp_x];
    const PylithInt sOffDispP_x = sOffDispN_x + spaceDim ^ 2;
    const PylithInt fOffp_fault = 0;

    const PylithScalar *dispN_x = &s_x[sOffDispN_x];
    const PylithScalar *dispP_x = &s_x[sOffDispP_x];

    // Trace_strain, no transformation required
    const PylithInt sOffTraceStrainN = _FaultPoroCohesiveKin::trace_strain_sOff(sOff, numS);
    const PylithInt sOffTraceStrainP = sOffTraceStrainN + 1;
    const PylithScalar traceStrainN = s[sOffTraceStrainN];
    const PylithScalar traceStrainP = s[sOffTraceStrainP];

    // fault pressure
    const PylithInt sOffFaultPressure = _FaultPoroCohesiveKin::fault_pressure_sOff(sOff, numS);
    const PylithScalar faultPressure = s[sOffFaultPressure];

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

    const PylithScalar stress_nnN = (undrainedBulkModulus - 2. * shearModulus / 3.) * traceStrainN + 2. * shearModulus * strain_nnN;
    const PylithScalar stress_nnP = (undrainedBulkModulus - 2. * shearModulus / 3.) * traceStrainP + 2. * shearModulus * strain_nnP;

    // Constants M_u and M_u' for f0o_fault
    const PylithScalar M_u = undrainedBulkModulus + 4. * shearModulus / 3.;
    const PylithScalar M_u_prime = faultUndrainedBulkModulus + 4. * faultShearModulus / 3.;

    f0[0] = faultPressure - faultSkemptonCoef * faultUndrainedBulkModulus / M_u_prime *
                                (faultShearModulus * M_u / shearModulus * (traceStrainN + traceStrainP) * 3. / 2. + (shearModulus - faultShearModulus) / shearModulus * (stress_nnN + stress_nnP) / 2.);

} // f0p_fault

// ----------------------------------------------------------------------
// f0 function for slip rate constraint equation: f0\lambda = (v^+ - v^-) - \dot{d}
void pylith::fekernels::FaultPoroCohesiveKin::f0l_v(const PylithInt dim,
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

    const PylithInt sOffVelN = _FaultPoroCohesiveKin::velocity_sOff(sOff, numS);
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
        pylith::fekernels::_FaultCohesiveKin::tangential_directions(_spaceDim, refDir1, refDir2, n, tanDir1, tanDir2);

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
/* Jf0 function for integration of the displacement equation.
 *
 * Solution fields = [disp(dim), ..., lagrange(dim)]
 * Auxiliary fields = None
 */
void pylith::fekernels::FaultPoroCohesiveKin::Jf0ul(const PylithInt dim,
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
    assert(numA >= 6);
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
/* Jf0p_fp_f function for integration of the displacement equation.
 *
 * Solution fields = [disp(dim), ..., lagrange(dim), fault_pressure(1)]
 * Auxiliary fields
 */
void pylith::fekernels::FaultPoroCohesiveKin::Jf0p_fp_f(const PylithInt dim,
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
    assert(numA >= 6);
    assert(Jf0);
    assert(sOff);
    assert(aOff);
    assert(n);

    Jf0[0] += 1.;

} // Jf0p_fp_f

// ----------------------------------------------------------------------
/* Jf0p_fe function for integration of the displacement equation.
 *
 * Solution fields = [disp(dim), ..., lagrange(dim), fault_pressure(1)]
 * Auxiliary fields
 */
void pylith::fekernels::FaultPoroCohesiveKin::Jf0p_fe(const PylithInt dim,
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
    assert(numA >= 6);
    assert(Jf0);
    assert(sOff);
    assert(aOff);
    assert(n);

    // Index for auxiliary fields
    const PylithInt i_fault_undrained_bulk_modulus = 0;
    const PylithInt i_fault_shear_modulus = 1;
    const PylithInt i_fault_skempton_coefficient = 2;
    const PylithInt i_solid_undrained_bulk_modulus = 3;
    const PylithInt i_solid_shear_modulus = 4;

    const PylithScalar faultUndrainedBulkModulus = a[aOff[i_fault_undrained_bulk_modulus]];
    const PylithScalar faultShearModulus = a[aOff[i_fault_shear_modulus]];
    const PylithScalar faultSkemptonCoef = a[aOff[i_fault_skempton_coefficient]];
    const PylithScalar undrainedBulkModulus = a[aOff[i_solid_undrained_bulk_modulus]];
    const PylithScalar shearModulus = a[aOff[i_solid_shear_modulus]];

    // Constants M_u and M_u' for f0o_fault
    const PylithScalar M_u = undrainedBulkModulus + 4. * shearModulus / 3.;
    const PylithScalar M_u_prime = faultUndrainedBulkModulus + 4. * faultShearModulus / 3.;

    Jf0[0] += -faultSkemptonCoef * faultUndrainedBulkModulus / M_u_prime *
              (faultShearModulus * M_u * 3. / 2. / shearModulus + (shearModulus - faultShearModulus) / shearModulus *
                                                                      (undrainedBulkModulus - 2. * shearModulus / 3.) / 2.);

    Jf0[1] += -faultSkemptonCoef * faultUndrainedBulkModulus / M_u_prime *
              (faultShearModulus * M_u * 3. / 2. / shearModulus + (shearModulus - faultShearModulus) / shearModulus *
                                                                      (undrainedBulkModulus - 2. * shearModulus / 3.) / 2.);

} // Jf0p_fe

// ----------------------------------------------------------------------
/* Jf1p_fu function for integration of the displacement equation.
 *
 * Solution fields = [disp(dim), ..., lagrange(dim), fault_pressure(1)]
 * Auxiliary fields
 */
void pylith::fekernels::FaultPoroCohesiveKin::Jf1p_fu(const PylithInt dim,
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
    assert(numS >= 5);
    assert(numA >= 6);
    assert(Jf1);
    assert(sOff);
    assert(aOff);
    assert(n);

    // Index for auxiliary fields
    const PylithInt i_fault_undrained_bulk_modulus = 0;
    const PylithInt i_fault_shear_modulus = 1;
    const PylithInt i_fault_skempton_coefficient = 2;
    // const PylithInt i_solid_undrained_bulk_modulus = 3;
    const PylithInt i_solid_shear_modulus = 4;

    const PylithScalar faultUndrainedBulkModulus = a[aOff[i_fault_undrained_bulk_modulus]];
    const PylithScalar faultShearModulus = a[aOff[i_fault_shear_modulus]];
    const PylithScalar faultSkemptonCoef = a[aOff[i_fault_skempton_coefficient]];
    // const PylithScalar undrainedBulkModulus = a[aOff[i_solid_undrained_bulk_modulus]];
    const PylithScalar shearModulus = a[aOff[i_solid_shear_modulus]];

    // Constants M_u and M_u' for f0o_fault
    // const PylithScalar M_u = undrainedBulkModulus + 4. * shearModulus / 3.;
    const PylithScalar M_u_prime = faultUndrainedBulkModulus + 4. * faultShearModulus / 3.;

    // Start computing Jf0
    const PylithInt spaceDim = dim + 1; // :KLUDGE: dim passed in is spaceDim-1

    const PylithInt gOffN = 0;
    const PylithInt gOffP = gOffN + spaceDim;
    const PylithInt ncols = spaceDim;

    Jf1[(gOffN + spaceDim - 1) * ncols + spaceDim - 1] += -faultSkemptonCoef *
                                                          faultUndrainedBulkModulus *
                                                          (shearModulus - faultShearModulus) /
                                                          M_u_prime;

    Jf1[(gOffP + spaceDim - 1) * ncols + spaceDim - 1] += -faultSkemptonCoef *
                                                          faultUndrainedBulkModulus *
                                                          (shearModulus - faultShearModulus) /
                                                          M_u_prime;

} // Jf1p_fu

// ----------------------------------------------------------------------
/* Jg0 function for integration of the slip constraint equation.
 *
 * Solution fields = [disp(dim), ..., lagrange(dim)]
 * Auxiliary fields = None
 */
void pylith::fekernels::FaultPoroCohesiveKin::Jf0lu(const PylithInt dim,
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
    assert(numA >= 6);
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
