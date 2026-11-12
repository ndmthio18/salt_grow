// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright (c) DuMux Project contributors
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Chemistry
 * \brief Simplified Pitzer model for aqueous NaCl: single shared
 *        source of truth for the ion activity coefficient
 *        \f$\gamma_{\pm}\f$, the osmotic coefficient \f$\phi\f$, and
 *        the water activity \f$a_w\f$. The \f$\gamma_{\pm}\f$
 *        expression matches the inline formula previously used in the
 *        2pncmin precipitation source term verbatim; the \f$\phi\f$
 *        expression is its Gibbs--Duhem-consistent counterpart, so
 *        the two are thermodynamically coupled by construction.
 *
 *        Parameters are hard-coded for NaCl at T = 298.15 K. The
 *        TODO at the end of the parameter block flags the missing
 *        temperature dependence.
 */
#ifndef DUMUX_PITZER_NACL_HH
#define DUMUX_PITZER_NACL_HH

#include <cmath>
#include <dumux/material/components/h2o.hh>

namespace Dumux {

/*!
 * \ingroup Chemistry
 * \brief Simplified Pitzer model for aqueous NaCl.
 *
 *  The \f$\gamma_{\pm}\f$ formula is bit-for-bit equivalent to the
 *  inline expression that was previously hard-coded in
 *  appl/2pncmin/problem.hh; refactoring it into this header guarantees
 *  that the corresponding water activity used in the vapor pressure
 *  uses the same numerical parameters.
 *
 *  The \f$\phi\f$ formula is derived from the same \f$\gamma_{\pm}\f$
 *  expression by integrating the electrolyte Gibbs--Duhem relation
 *  \f$\,m\,d(\phi-1) = m\,d\ln\gamma_{\pm} - (\phi-1)\,dm\f$, term by
 *  term. For NaCl (1:1 electrolyte, \f$\nu = 2\f$) this gives
 *  \f$a_w = \exp(-\phi(m)\,M_w\,2m)\f$.
 */
template<class Scalar>
struct PitzerNaCl
{
    // Pitzer parameters for NaCl at T = 298.15 K. The values match
    // appl/2pncmin/problem.hh:408-414 verbatim. TODO: replace by a
    // temperature-dependent expression for non-isothermal use.
    static constexpr Scalar A      = 0.39159;   // Debye-Hueckel constant [(kg/mol)^0.5]
    static constexpr Scalar b      = 2.22718;   // [(kg/mol)^0.5]
    static constexpr Scalar betaMX = 0.05383;   // single-B Pitzer parameter
    static constexpr Scalar CMX    = 0.00134;   // third-order Pitzer parameter
    // Charge / stoichiometry constants for NaCl (1:1 electrolyte).
    static constexpr Scalar zM     = 1.0;
    static constexpr Scalar zX     = -1.0;
    static constexpr Scalar v      = 2.0;
    static constexpr Scalar vM     = 1.0;
    static constexpr Scalar vX     = 1.0;

    /*!
     * \brief Mean ionic activity coefficient \f$\gamma_{\pm}\f$ of
     *        aqueous NaCl as a function of molality.
     *
     * \param molality NaCl molality [mol/kg H2O]; for NaCl the ionic
     *                 strength equals the molality.
     */
    static Scalar activityCoefficient(Scalar molality)
    {
        using std::sqrt;
        using std::log;
        using std::exp;
        using std::pow;
        using std::abs;

        const Scalar I  = molality;            // ionic strength = molality for NaCl
        const Scalar sI = sqrt(I);

        // Bit-for-bit equivalent to the original inline expression at
        // appl/2pncmin/problem.hh:415.
        return exp(
            -abs(zM*zX) * A * (sI/(1 + b*sI) + 2.0/b * log(1 + b*sI))
            + (4.0*vM*vX/v) * molality * betaMX
            + (6.0*pow(vM*vX, 1.5)/v) * sqrt(abs(zM*zX)) * molality*molality * CMX
        );
    }

    /*!
     * \brief Osmotic coefficient \f$\phi\f$ of aqueous NaCl as a
     *        function of molality, derived from the
     *        \f$\gamma_{\pm}\f$ above by Gibbs--Duhem integration.
     *
     *  For NaCl (1:1, \f$\nu = 2\f$) the derivation reduces to
     *  halving the linear-in-m coefficient and applying a 2/3
     *  factor to the m^2 coefficient relative to \f$\ln\gamma_{\pm}\f$,
     *  while replacing the Debye-Hueckel block
     *  \f$\sqrt{I}/(1+b\sqrt{I}) + (2/b)\ln(1+b\sqrt{I})\f$ by the
     *  bare \f$\sqrt{I}/(1+b\sqrt{I})\f$.
     *
     *  Numerical check: at saturation (m = 6.14 mol/kg) this returns
     *  \f$\phi \approx 1.283\f$, giving \f$a_w \approx 0.753\f$ as
     *  reported in the manuscript.
     *
     * \param molality NaCl molality [mol/kg H2O].
     */
    static Scalar osmoticCoefficient(Scalar molality)
    {
        using std::sqrt;
        using std::abs;
        using std::max;
        using std::min;

        // Defensive clamp: Newton may probe non-physical states with
        // molality < 0 (would NaN sqrt) or molality far above NaCl
        // saturation (~6.14 mol/kg) where the parameter set isn't
        // calibrated. Out-of-range trial states must NOT return NaN, or
        // the Jacobian collapses and the linear solve aborts. Inside
        // the calibrated range the result is bit-identical.
        const Scalar m = max(Scalar(0), min(molality, Scalar(6.5)));

        const Scalar sI = sqrt(m);

        return 1.0
             - abs(zM*zX) * A * sI/(1 + b*sI)
             + m * betaMX
             + 2.0 * m * m * CMX;
    }

    /*!
     * \brief Water activity \f$a_w\f$ of aqueous NaCl as a function of
     *        molality, computed as \f$a_w = \exp(-\phi\,M_w\,\nu\,m)\f$.
     *
     *  This is the Gibbs--Duhem-consistent counterpart to the ideal
     *  Raoult approximation \f$a_w = x_l^w\f$. At NaCl saturation
     *  (m = 6.14 mol/kg) it returns approximately 0.753, vs the
     *  ideal value 0.900.
     *
     * \param molality NaCl molality [mol/kg H2O].
     */
    static Scalar waterActivity(Scalar molality)
    {
        using std::exp;
        using std::max;
        using std::min;
        const Scalar Mw = Components::H2O<Scalar>::molarMass(); // kg/mol
        // Mirror the clamp inside osmoticCoefficient so the exponent
        // m·φ·v·Mw stays finite for any Newton trial state.
        const Scalar m = max(Scalar(0), min(molality, Scalar(6.5)));
        const Scalar phi = osmoticCoefficient(m);
        // sum of solute molalities for fully dissociated NaCl is 2*m.
        return exp(-phi * Mw * v * m);
    }
};

} // end namespace Dumux

#endif
