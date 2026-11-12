// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup SolidStates
 * \brief Represents all relevant thermodynamic quantities of a noninert solid system.
 */
#ifndef DUMUX_NONINERT_SOLID_STATE_HH
#define DUMUX_NONINERT_SOLID_STATE_HH

#include <dumux/material/solidstates/inertsolidstate.hh>

namespace Dumux {

/*!
 * \ingroup SolidStates
 * \brief Represents all relevant thermodynamic quantities of a inert solid system.
 */
template <class Scalar, class SolidSystemType>
class NonInertSolidState: public InertSolidState<Scalar, SolidSystemType>
{
public:
    using SolidSystem = SolidSystemType;
    using ParentType = InertSolidState<Scalar, SolidSystemType>;

    enum
    {
        numComponents = SolidSystem::numComponents,
        numInertComponents = SolidSystem::numInertComponents,
    };

    /*!
     * \brief Allows compile-time evaluation of if the solid system
     *        is inert or takes part in any kind of reactions.
     */
    static constexpr bool isInert()
    {
        static_assert(!SolidSystem::isInert(), "Only noninert solid systems are allowed with the NonInertSolidState");
        return false;
    }


    /*!
     * \brief The porosity of the porous medium
     */
    Scalar porosity() const
    {
        bool poreVolumeChange = getParam<bool>("Problem.PoreVolumeChange", true);

        Scalar porosity = 1.0;
        if(poreVolumeChange){
            Scalar sumVolumeFraction = 0.0;
            for (int compIdx =0; compIdx < numComponents; ++compIdx)
                sumVolumeFraction += this->volumeFraction(compIdx);
            porosity = 1-sumVolumeFraction;
        }
        return porosity;
    }

};

} // end namespace Dumux

#endif
