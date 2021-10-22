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

/** @file libsrc/fekernels/FaultPoroCohesiveKin.hh
 *
 * Kernels for poroelastic faults with prescribed slip.
 *
 * Solution fields: [disp(dim), vel(dim, optional), pressure(1), trace_strain(1), lagrange(dim), fault_pressure(1)]
 *
 * Auxiliary fields: [fault_undrained_bulk_modulus,
 *                    fault_shear_modulus,
 *                    fault_skempton_coefficient (B),
 *                    solid_undrained_bulk_modulus,
 *                    solid_shear_modulus,
 *                    slip]                   (**NEED IMPLEMENTATION**)
 *
 * LHS Residual: no contribution
 *
 * RHS Residual
 *
 *  - g0u^+ = -\lambda
 *  - g0u^- = +\lambda
 *  - g0\lambda = d - u^+ + u^-
 *
 * LHS Jacobian: no contribution
 *
 * RHS Jacobian
 *
 *  - Jg0^{u \lambda} = -1 for u^+, +1 for u^-
 *  - Jg0^{\lambda u} = -1 for u^+, +1 for u^-
 *
 * ======================================================================
 */

#if !defined(pylith_fekernels_faultporocohesivekin_hh)
#define pylith_fekernels_faultporocohesivekin_hh

// Include directives ---------------------------------------------------
#include "fekernelsfwd.hh" // forward declarations

#include "pylith/utils/types.hh"

class pylith::fekernels::FaultPoroCohesiveKin {
    // PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

    /** Kernel interface.
     *
     * @param[in] dim Spatial dimension.
     * @param[in] numS Number of registered subfields in solution field.
     * @param[in] numA Number of registered subfields in auxiliary field.
     * @param[in] sOff Offset of registered subfields in solution field [numS].
     * @param[in] sOff_x Offset of registered subfields in gradient of the solution field [numS].
     * @param[in] s Solution field with all subfields.
     * @param[in] s_t Time derivative of solution field.
     * @param[in] s_x Gradient of solution field.
     * @param[in] aOff Offset of registered subfields in auxiliary field [numA]
     * @param[in] aOff_x Offset of registered subfields in gradient of auxiliary field [numA]
     * @param[in] a Auxiliary field with all subfields.
     * @param[in] a_t Time derivative of auxiliary field.
     * @param[in] a_x Gradient of auxiliary field.
     * @param[in] t Time for residual evaluation.
     * @param[in] x Coordinates of point evaluation.
     * @param[in] numConstants Number of registered constants.
     * @param[in] constants Array of registered constants.
     * @param[out] f0 [dim].
     */

    /** f0 function for poroelasticity equation: f0u = -\lambda (pos side), +\lambda (neg side).
     *
     * Solution fields: [disp(dim), ..., lagrange(dim), fault_pressure(1)]
     */
    static void f0u(const PylithInt dim,
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
                    PylithScalar f0[]);

    /** f0 function for slip constraint equation: f0\lambda = (u^+ - u^-) - d
     *
     * Solution fields: [disp(dim), ..., lagrange(dim), fault_pressure(1)]
     */
    static void f0l_u(const PylithInt dim,
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
                      PylithScalar f0[]);

    /** f0 function for slip rate constraint equation: f0\lambda = (v^+ - v^-) - \dot{d}
     *
     * Solution fields: [disp(dim), vel(dim), ..., lagrange(dim), fault_pressure(1)]
     */
    static void f0l_v(const PylithInt dim,
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
                      PylithScalar f0[]);

    /** f0 function for fault pressure constraint equation
     *  f0p_fault = p' - B'K'_u / M_u'
     *  * [G' M_u / (G K_u) * 3K_u(trace_strain^- + trace_strain^-) / 2 +
     *  (G - G') / G * ((K_u - 2G/3) * (trace_strain^- + trace_strain^-)
     *  + 2G (u3,3^- + u3,3^+))/2]
     *
     *  Solution fields: [disp(dim), vel(dim), ..., lagrange(dim), fault_pressure(1)]
     */
    static void f0p_fault(const PylithInt dim,
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
                          PylithScalar f0[]);

    /** Jf0 function for displacement equation: +\lambda (pos side), -\lambda (neg side).
     */
    static void Jf0ul(const PylithInt dim,
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
                      PylithScalar Jf0[]);

    /** Jf0 function for p' p': 1.
     *
     * Solution fields: [disp(dim), ..., lagrange(dim), fault_pressure]
     */
    static void Jf0p_fp_f(const PylithInt dim,
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
                          PylithScalar Jf0[]);

    /** Jf0 function for p' trace_strain.
     *
     * Solution fields: [disp(dim), ..., lagrange(dim), fault_pressure(1)]
     */
    static void Jf0p_fe(const PylithInt dim,
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
                        PylithScalar Jf0[]);

    /** Jf1 function for fault pressure - u.
     *
     * Solution fields: [disp(dim), ..., lagrange(dim), fault_pressure(1)]
     */
    static void Jf1p_fu(const PylithInt dim,
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
                        PylithScalar Jf1[]);

    /** Jf0 function for slip constraint equation: +\lambda (pos side), -\lambda (neg side).
     *
     * Solution fields: [disp(dim), ..., lagrange(dim), fault_pressure(1)]
     */
    static void Jf0lu(const PylithInt dim,
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
                      PylithScalar Jf0[]);

}; // FaultPoroCohesiveKin

#endif // pylith_fekernels_faultporocohesivekin_hh

/* End of file */
