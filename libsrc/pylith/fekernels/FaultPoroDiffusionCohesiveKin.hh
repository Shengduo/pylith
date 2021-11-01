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

/** @file libsrc/fekernels/FaultPoroDiffusionCohesiveKin.hh
 *
 * Kernels for poroelastic faults 3D fluid diffusion and prescribed slip.
 *
 * Solution fields: [disp(dim), vel(dim, optional), pressure(1), trace_strain(1), lagrange(dim), fault_pressure(1)]
 *
 * Auxiliary fields: [
 *                    thickness(1)                        0
 *                    porosity(1),                        1
 *                    beta_p(1),                          2
 *                    beta_sigma(1),                      3
 *                    permeability_tangential(1),         4
 *                    permeability_normal(1),             5
 *                    fluid_viscosity(1),                 6
 *                    body_force(dim),                    numA - 4
 *                    source (1),                         numA - 3
 *                    constant_pressure_source(1)         numA - 2
 *                    slip(dim)                           numA - 1
 *                    ]                                   (**NEED IMPLEMENTATION**)
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

#if !defined(pylith_fekernels_faultporodiffusioncohesivekin_hh)
#define pylith_fekernels_faultporodiffusioncohesivekin_hh

// Include directives ---------------------------------------------------
#include "fekernelsfwd.hh" // forward declarations

#include "pylith/utils/types.hh"

class pylith::fekernels::FaultPoroDiffusionCohesiveKin {
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

    /** f0 function for poroelasticity equation: f0u = -\lambda (neg side), +\lambda (pos side).
     *
     * Solution fields: [disp(dim), ..., lagrange(dim), fault_pressure(1)]
     */
    static void f0u_neg(const PylithInt dim,
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

    /** f0 function for poroelasticity equation: f0u = -\lambda (neg side), +\lambda (pos side).
     *
     * Solution fields: [disp(dim), ..., lagrange(dim), fault_pressure(1)]
     */
    static void f0u_pos(const PylithInt dim,
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
    
    /** f0 function for bulk pressure equation: f0p = [\kappa_{cz} / \mu * ((p^+ - p^f)/h - n \cdot f_f),
     *                                                 \kappa_{cz} / \mu * ((p^- - p^f)/h + n \cdot f_f)]
     *
     * Solution fields: [disp(dim), vel(dim), ..., lagrange(dim), fault_pressure(1)]
     */
    static void f0p_neg(const PylithInt dim,
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

    /** f0 function for bulk pressure equation: f0p = [\kappa_{cz} / \mu * ((p^+ - p^f)/h - n \cdot f_f),
     *                                                 \kappa_{cz} / \mu * ((p^- - p^f)/h + n \cdot f_f)]
     *
     * Solution fields: [disp(dim), vel(dim), ..., lagrange(dim), fault_pressure(1)]
     */
    static void f0p_pos(const PylithInt dim,
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

    /** f0 function for bulk pressure equation: f0p = [\kappa_{cz} / \mu * ((p^+ - p^f)/h - n \cdot f_f),
     *                                                 \kappa_{cz} / \mu * ((p^- - p^f)/h + n \cdot f_f)]
     *
     * Solution fields: [disp(dim), vel(dim), ..., lagrange(dim), fault_pressure(1)]
     */
    static void f0p_body_neg(const PylithInt dim,
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
    
    /** f0 function for bulk pressure equation: f0p = [\kappa_{cz} / \mu * ((p^+ - p^f)/h - n \cdot f_f),
     *                                                 \kappa_{cz} / \mu * ((p^- - p^f)/h + n \cdot f_f)]
     *
     * Solution fields: [disp(dim), vel(dim), ..., lagrange(dim), fault_pressure(1)]
     */
    static void f0p_body_pos(const PylithInt dim,
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
     *  f0p_fault = porosity * (\beta^p * (\dot{p}^+ + 2 \dot{p}^f + \dot{p}^-)/4
     *                        + \beta^\sigma * (\dot{\sigma}^+_{nn} + \dot{\sigma}^-_{nn}))
     *              + \kappa_{fx} / \mu * \vnabla(2D) \cdot body_force
     *              - \kappa_{fz} / \mu * (p^+ - 2p_f + p^-) / h^2 - source
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

    /** f0 function for fault pressure constraint equation
     *  f0p_fault = porosity * (\beta^p * (\dot{p}^+ + 2 \dot{p}^f + \dot{p}^-)/4
     *                        + \beta^\sigma * (\dot{\sigma}^+_{nn} + \dot{\sigma}^-_{nn}))
     *              + \kappa_{fx} / \mu * \vnabla(2D) \cdot body_force
     *              - \kappa_{fz} / \mu * (p^+ - 2p_f + p^-) / h^2 - source
     *
     *  Solution fields: [disp(dim), vel(dim), ..., lagrange(dim), fault_pressure(1)]
     */
    static void f0p_fault_body(const PylithInt dim,
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
     *  f0p_fault = porosity * (\beta^p * (\dot{p}^+ + 2 \dot{p}^f + \dot{p}^-)/4
     *                        + \beta^\sigma * (\dot{\sigma}^+_{nn} + \dot{\sigma}^-_{nn}))
     *              + \kappa_{fx} / \mu * \vnabla(2D) \cdot body_force
     *              - \kappa_{fz} / \mu * (p^+ - 2p_f + p^-) / h^2 - source
     *
     *  Solution fields: [disp(dim), vel(dim), ..., lagrange(dim), fault_pressure(1)]
     */
    static void f0p_fault_source(const PylithInt dim,
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
     *  f0p_fault = porosity * (\beta^p * (\dot{p}^+ + 2 \dot{p}^f + \dot{p}^-)/4
     *                        + \beta^\sigma * (\dot{\sigma}^+_{nn} + \dot{\sigma}^-_{nn}))
     *              + \kappa_{fx} / \mu * \vnabla(2D) \cdot body_force
     *              - \kappa_{fz} / \mu * (p^+ - 2p_f + p^-) / h^2 - source
     *
     *  Solution fields: [disp(dim), vel(dim), ..., lagrange(dim), fault_pressure(1)]
     */
    static void f0p_fault_body_source(const PylithInt dim,
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

    /** f1 function for fault pressure constraint equation
     *  f1p_fault = \kappa_{fx} / (4\mu) \vnabla (p^+ + 2 p^f + p^-)
     *
     *  Solution fields: [disp(dim), vel(dim), ..., lagrange(dim), fault_pressure(1)]
     */
    static void f1p_fault(const PylithInt dim,
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
                          PylithScalar f1[]);

    /** Jf0 function for displacement equation: +\lambda (pos side), -\lambda (neg side).
     */
    static void Jf0ul_neg(const PylithInt dim,
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

    /** Jf0 function for displacement equation: +\lambda (pos side), -\lambda (neg side).
     */
    static void Jf0ul_pos(const PylithInt dim,
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

    /** Jf0 function for p_f p_f :
     * 2 \kappa_{fz} / (\mu h^2) + 2 \phi_f \beta^p t_shift
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

    /** Jf3 function for p_f p_f: \kappa_{fx}/(2\mu) \te{I}.
     *
     * Solution fields: [disp(dim), ..., lagrange(dim), fault_pressure]
     */
    static void Jf3p_fp_f(const PylithInt dim,
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
                          PylithScalar Jf3[]);

    /** Jf0 function for fault pressure - lambda.
     * s_tshift \phi_f \beta^\sigma \ve{n} 
     * Solution fields: [disp(dim), ..., lagrange(dim), fault_pressure(1)]
     */
    static void Jf0p_fl(const PylithInt dim,
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


    /** Jf0 function for fault pressure - pressure :
     * [\phi_f beta^p t_shift /4 - \kappa_{fz}/\mu h^2,\phi_f beta^p t_shift /4 - \kappa_{fz}/\mu h^2]
     * Solution fields: [disp(dim), ..., lagrange(dim), fault_pressure(1)]
     */
    static void Jf0p_fp(const PylithInt dim,
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

    /** Jf3 function for fault pressure - p :
     *  [\kappa_{fx}/4\mu \te{I}, \kappa_{fx}/4\mu \te{I}]
     * Solution fields: [disp(dim), ..., lagrange(dim), fault_pressure(1)]
     */
    static void Jf3p_fp(const PylithInt dim,
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
                        PylithScalar Jf3[]);

    /** Jf0 function for pressure pressure:
     * [\kappa_{cz} / \mu / h, 0;
     *  0, \kappa_{cz} / \mu / h]
     *
     * Solution fields: [disp(dim), ..., lagrange(dim), fault_pressure(1)]
     */
    static void Jf0pp_neg(const PylithInt dim,
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
    /** Jf0 function for pressure pressure:
     * [\kappa_{cz} / \mu / h, 0;
     *  0, \kappa_{cz} / \mu / h]
     *
     * Solution fields: [disp(dim), ..., lagrange(dim), fault_pressure(1)]
     */
    static void Jf0pp_pos(const PylithInt dim,
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


    /** Jf0 function for pressure - fault_pressure.
     * [-\kappa_cz / \mu / h, -\kappa_cz / \mu / h]
     * Solution fields: [disp(dim), ..., lagrange(dim), fault_pressure(1)]
     */
    static void Jf0pp_f_neg(const PylithInt dim,
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

    /** Jf0 function for pressure - fault_pressure.
     * [-\kappa_cz / \mu / h, -\kappa_cz / \mu / h]
     * Solution fields: [disp(dim), ..., lagrange(dim), fault_pressure(1)]
     */
    static void Jf0pp_f_pos(const PylithInt dim,
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

}; // FaultPoroDiffusionCohesiveKin

#endif // pylith_fekernels_faultporodiffusioncohesivekin_hh

/* End of file */
