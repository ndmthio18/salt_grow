// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PNMTwoPNCModel
 * \copydoc Dumux::PoreNetwork::TwoPNCIOFields
 */
#ifndef DUMUX_PNM_2P_NC_CUSTOM_IO_FIELDS_HH
#define DUMUX_PNM_2P_NC_CUSTOM_IO_FIELDS_HH

#include <dumux/porousmediumflow/2pnc/iofields.hh>
#include <dumux/porenetwork/common/iofields.hh>

namespace Dumux::PoreNetwork {

/*!
 * \ingroup PNMTwoPNCModel
 * \brief Adds output fields specific to the PNM 2pnc model
 */
class TwoPNCCustomIOFields
{
public:
    template <class OutputModule>
    static void initOutputModule(OutputModule& out)
    {
         // use default fields from the 2pnc model
        Dumux::TwoPNCIOFields::initOutputModule(out);

        using VolumeVariables = typename OutputModule::VolumeVariables;
        using FS = typename VolumeVariables::FluidSystem;

        CommonIOFields::initOutputModule(out);

        out.addFluxVariable([](const auto& fluxVars, const auto& fluxVarsCache)
                              { return fluxVarsCache.pcEntry(); }, "pcEntry");
        out.addFluxVariable([](const auto& fluxVars, const auto& fluxVarsCache)
                              { return fluxVarsCache.pcSnapoff(); }, "pcSnapoff");

        out.addFluxVariable([](const auto& fluxVars, const auto& fluxVarsCache)
                              { return fluxVarsCache.transmissibility(FS::phase0Idx); }, "transmissibilityW");

        out.addFluxVariable([](const auto& fluxVars, const auto& fluxVarsCache)
                              { return fluxVarsCache.transmissibility(FS::phase1Idx); }, "transmissibilityN");



        out.addVolumeVariable([](const auto& v){ return v.pressure(1)-1e5; }, "delPN");
        out.addVolumeVariable([](const auto& v){ return v.pressure(0)-1e5; }, "delPW");
        out.addVolumeVariable([](const auto& v){ return v.surfaceTension(); }, "surfaceTension");

        out.addFluxVariable([](const auto& fluxVars, const auto& fluxVarsCache) { return fluxVars.molecularDiffusionFlux(0)[1]; }, "diffusiveFlux_0^1");
        out.addFluxVariable([](const auto& fluxVars, const auto& fluxVarsCache) { return fluxVars.molecularDiffusionFlux(1)[0]; }, "diffusiveFlux_1^0");
        out.addFluxVariable([](const auto& fluxVars, const auto& fluxVarsCache) { return fluxVars.heatConductionFlux(); }, "heatConductionFlux");
        out.addFluxVariable([](const auto& fluxVars, const auto& fluxVarsCache) { return fluxVarsCache.throatCrossSectionalArea(0); }, "throatCrossSectionalArea0");
        out.addFluxVariable([](const auto& fluxVars, const auto& fluxVarsCache) { return fluxVarsCache.throatCrossSectionalArea(1); }, "throatCrossSectionalArea1");
        out.addFluxVariable([](const auto& fluxVars, const auto& fluxVarsCache) { return fluxVarsCache.throatCrossSectionalArea(0)/(fluxVarsCache.throatCrossSectionalArea(0)+fluxVarsCache.throatCrossSectionalArea(1)); }, "throatCrossSectionalSaturation0");
        out.addFluxVariable([](const auto& fluxVars, const auto& fluxVarsCache) { return fluxVarsCache.throatCrossSectionalArea(1)/(fluxVarsCache.throatCrossSectionalArea(0)+fluxVarsCache.throatCrossSectionalArea(1)); }, "throatCrossSectionalSaturation1");

        out.addVolumeVariable([](const auto& volVars){ return volVars.poreVolume(); }, "poreVolume");




        auto volumeFluxW = [](const auto& fluxVars, const auto& fluxVarsCache)
        {
            auto upwindTerm = [](const auto& volVars) { return volVars.mobility(FS::phase0Idx); };
            using std::abs;
            return abs(fluxVars.advectiveFlux(FS::phase0Idx, upwindTerm));
        };
        out.addFluxVariable(volumeFluxW, "volumeFluxW");

        auto volumeFluxN = [](const auto& fluxVars, const auto& fluxVarsCache)
        {
            auto upwindTerm = [](const auto& volVars) { return volVars.mobility(FS::phase1Idx); };
            using std::abs;
            return abs(fluxVars.advectiveFlux(FS::phase1Idx, upwindTerm));
        };
        out.addFluxVariable(volumeFluxN, "volumeFluxN");
    }
};

} // end namespace Dumux::PoreNetwork

#endif
