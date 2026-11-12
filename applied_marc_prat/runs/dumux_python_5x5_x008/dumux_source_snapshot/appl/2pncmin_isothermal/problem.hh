// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 *
 * \brief A test problem for the two-phase n-components pore network model.
 */
#ifndef DUMUX_PNM_2P_NC_PROBLEM_HH
#define DUMUX_PNM_2P_NC_PROBLEM_HH

#include <dumux/common/boundarytypes.hh>
#include <dumux/common/parameters.hh>
#include <dumux/2pncmin/model.hh>
#include <dumux/porousmediumflow/problem.hh>
#include <dumux/material/components/air.hh>
#include <test/multidomain/porenetwork/2p/regularization.hh>

#include "../../dumux/material/pitzernacl.hh"

namespace Dumux {

template <class TypeTag>
class DrainageProblem;

template <class TypeTag>
class DrainageProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using FluidState = GetPropType<TypeTag, Properties::FluidState>;

    // copy some indices for convenience
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using Labels = GetPropType<TypeTag, Properties::Labels>;
    using Element = typename GridView::template Codim<0>::Entity;
    using Vertex = typename GridView::template Codim<GridView::dimension>::Entity;

    using GridFluxVariablesCache = GetPropType<TypeTag, Properties::GridFluxVariablesCache>;
    using InvasionState = std::decay_t<decltype(std::declval<GridFluxVariablesCache>().invasionState())>;
    static constexpr bool useThetaRegularization = InvasionState::stateMethod == Dumux::PoreNetwork::StateSwitchMethod::theta;

    enum
    {
        gasPhaseIdx = FluidSystem::gasPhaseIdx, // =1
        liquidPhaseIdx = FluidSystem::liquidPhaseIdx, //=0

        // component indices in fluid system
        H2OIdx = FluidSystem::H2OIdx, //=0
        AirIdx = FluidSystem::AirIdx, //=1
#if USEBRINE
        NaClIdx = FluidSystem::NaClIdx, //=2

        // component indices in the solid sysem
        NaClSolidSystemIdx = 0/*SolidSystem::comp0Idx*/,
#endif

        // Index of the balance equations
        conti0EqIdx = Indices::conti0EqIdx, //=0 water component mole balance equation
        conti1EqIdx = Indices::conti0EqIdx + AirIdx, //=1 air component mole balance equation
#if USEBRINE
        contiNaClEqIdx = Indices::conti0EqIdx + NaClIdx, //=2 NaCl component mole balance equation
        contiSolidEqIdx = Indices::conti0EqIdx + FluidSystem::numComponents, //=3 solid NaCl component mole balance equation
#endif
#if !ISOTHERMAL
        energyEqIdx = Indices::energyEqIdx, // =4 energy balance equation
#endif

        //FluidSystem
        numPhases = FluidSystem::numPhases,
        numComponents = FluidSystem::numComponents,

        // primary variable indices
        pressureIdx = Indices::pressureIdx, //=0
        switchIdx = Indices::switchIdx, //=1
#if USEBRINE
        solidNaClIdx = FluidSystem::numComponents, //=3
#endif
#if !ISOTHERMAL
        temperatureIdx = Indices::temperatureIdx, //=4
#endif

        // Phase State
        bothPhases = Indices::bothPhases,

        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,
    };


public:
    template<class SpatialParams>
    DrainageProblem(std::shared_ptr<const GridGeometry> gridGeometry, std::shared_ptr<SpatialParams> spatialParams)
    : ParentType(gridGeometry, spatialParams),
    reg_(getParam<std::string>("Problem.RegularizationFunction", "Sine"),
           getParam<Scalar>("Problem.RegularizationDelta", 1e-8))
    {
        episodeLength_ = getParam<Scalar>("TimeLoop.EpisodeLength", 0.0);
        useSmallGrid_ = getParamFromGroup<bool>(this->paramGroup(), "Grid.UseSmallGrid", true);
        vtpOutputFrequency_ = getParam<int>("Problem.VtpOutputFrequency");
        useKineticPrecipitationReaction_ = getParam<bool>("Problem.UseKineticPrecipitationReaction", false);
        activityCoefficientPitzerApproach_ = getParam<bool>("Problem.UseActivityCoefficientPitzerApproach", true);
#if USEBRINE
        // Toggle Pitzer-derived water activity in the brine vapor pressure
        // (Gibbs--Duhem-consistent with the salt ion activity coefficient).
        // Default true to match the manuscript claim; set false in a case's
        // .input file to recover the legacy a_w = x_l^w behaviour.
        FluidSystem::Brine::setUsePitzerWaterActivity(
            getParam<bool>("Brine.UsePitzerWaterActivity", true));
#endif
        useDirichletAtBottom_ = getParam<bool>("Problem.UseDirichletAtBottom", false);

        //initial conditions
        initialSaturation_ = getParamFromGroup<Scalar>(this->paramGroup(), "Problem.InitialSaturation");
        initialSaturationBottom_ = getParamFromGroup<Scalar>(this->paramGroup(), "Problem.InitialSaturationBottom");
        initialSaturationTop_ = getParamFromGroup<Scalar>(this->paramGroup(), "Problem.InitialSaturationTop");
        initialPressure_ = getParamFromGroup<Scalar>(this->paramGroup(), "Problem.InitialPressure", 1e5);
#if USEBRINE
        initialSaltMoleFraction_ = getParamFromGroup<Scalar>(this->paramGroup(), "Problem.InitialSaltMoleFraction", 0.0);
        initialSolidVolFrac_ = getParamFromGroup<Scalar>(this->paramGroup(), "Problem.InitialSolidVolumeFraction", 0.0);
        throatApproach_ = getParam<Scalar>("Problem.ThroatApproach", 2);
#endif

        //evaporation boundary condition at the top
        boundaryLayerThickness_ = getParam<Scalar>("FreeFlow.BoundaryLayerThickness", 1.6e-4);
        refT_ = getParam<Scalar>("FreeFlow.ReferenceTemperature", 293.15);
        refWaterVaporMoleFrac_ = getParam<Scalar>("FreeFlow.ReferenceWaterVaporMoleFraction", 0.0);
        refWaterVaporMassFrac_ = moleToMassFrac_(refWaterVaporMoleFrac_, gasPhaseIdx, H2OIdx);
        refGasPressure_ = getParam<Scalar>("FreeFlow.ReferenceGasPressure", 1e5);
        transmissibilityFactor_ = getParam<Scalar>("FreeFlow.TransmissibilityFactor", 1.0);


#if !ISOTHERMAL
        initialTemperature_ = getParamFromGroup<Scalar>(this->paramGroup(), "Problem.InitialTemperature", 273.15 + 20.0);
#endif
        //set initial conditions
        getInitial();
    }

    /*!
     * \name Simulation steering
     */
    // \{

    /*!
     * \name Problem parameters
     */
    // \{

    bool shouldWriteOutput(const int timeStepIndex, const GridVariables& gridVariables) const
    {
        if (vtpOutputFrequency_ < 0)
            return true;

        if (vtpOutputFrequency_ == 0)
            return (timeStepIndex == 0 || gridVariables.gridFluxVarsCache().invasionState().hasChanged());
        else
            return (timeStepIndex % vtpOutputFrequency_ == 0 || gridVariables.gridFluxVarsCache().invasionState().hasChanged());
    }



        //! Called after every time step
        //! Output the total global exchange term
        template<class GridVariables, class SolutionVector>
        void postTimeStep(const SolutionVector& curSol,
                          const GridVariables& gridVariables,
                          std::ofstream& evaporationFile)
        {
            Scalar evaporation_mol = 0.0;
            Scalar evaporation_mass = 0.0;
            Scalar evaporation_volume = 0.0;

            Scalar advectiveAirFlux_mol = 0.0;
            Scalar advectiveWaterFlux_mol = 0.0;
            Scalar diffusiveAirFlux_mol = 0.0;
            Scalar diffusiveWaterFlux_mol = 0.0;

            Scalar advectiveAirFlux_volume = 0.0;
            Scalar advectiveWaterFlux_volume = 0.0;
            Scalar diffusiveAirFlux_volume = 0.0;
            Scalar diffusiveWaterFlux_volume = 0.0;

            Scalar airFlux_mol = 0.0;
            Scalar airFlux_volume = 0.0;
#if !ISOTHERMAL
            Scalar energyFlux = 0.0;
            Scalar energyConductFlux = 0.0;
            Scalar energyConvectFlux = 0.0;
#endif

            //Free flow fluidstate
            FluidState fsRef;
            fsRef.setTemperature(refT_);
            fsRef.setSaturation(gasPhaseIdx, 1.0);
            fsRef.setPressure(gasPhaseIdx, refGasPressure_);
            fsRef.setMoleFraction(gasPhaseIdx, H2OIdx, refWaterVaporMoleFrac_);
            fsRef.setMoleFraction(gasPhaseIdx, AirIdx, 1-refWaterVaporMoleFrac_);

            auto hmax = this->gridGeometry().bBoxMax()[1];

                for (const auto& element : elements(this->gridGeometry().gridView())){
                    auto fvGeometry = localView(this->gridGeometry());
                    auto elemVolVars = localView(gridVariables.curGridVolVars());
                    auto elemFluxVarsCache = localView(gridVariables.gridFluxVarsCache());
                    fvGeometry.bind(element);
                    elemVolVars.bind(element, fvGeometry, curSol);
                    elemFluxVarsCache.bind(element, fvGeometry, elemVolVars);

                    for (auto&& scv : scvs(fvGeometry)){
                        if (scv.dofPosition()[1] > hmax - eps_){
                            auto fluxes = source(element,  fvGeometry, elemVolVars, scv);
                            evaporation_mol += fluxes[conti0EqIdx] * scv.volume();
                            evaporation_mass += fluxes[conti0EqIdx] * scv.volume() * Dumux::Components::H2O<Scalar>::molarMass();
                            evaporation_volume += fluxes[conti0EqIdx] * scv.volume() / Dumux::Components::H2O<Scalar>::gasMolarDensity(refT_, fsRef.partialPressure(gasPhaseIdx, H2OIdx));
                            airFlux_mol += fluxes[conti1EqIdx] * scv.volume();
                            airFlux_volume += fluxes[conti1EqIdx] * scv.volume() / Dumux::Components::Air<Scalar>::gasMolarDensity(refT_, fsRef.partialPressure(gasPhaseIdx, AirIdx));
#if !ISOTHERMAL
                            energyFlux += fluxes[Indices::energyEqIdx] * scv.volume();
                            energyConductFlux += energyConductionFlux(element,  fvGeometry, elemVolVars, scv) * scv.volume();
                            energyConvectFlux += energyConvectionFlux(element,  fvGeometry, elemVolVars, scv) * scv.volume();
#endif
                            advectiveAirFlux_mol += advAirFlux(element,  fvGeometry, elemVolVars, scv) * scv.volume();
                            advectiveWaterFlux_mol += advWaterFlux(element,  fvGeometry, elemVolVars, scv) * scv.volume();
                            diffusiveAirFlux_mol += diffAirFlux(element,  fvGeometry, elemVolVars, scv) * scv.volume();
                            diffusiveWaterFlux_mol += diffWaterFlux(element,  fvGeometry, elemVolVars, scv) * scv.volume();
                            advectiveAirFlux_volume += advAirFlux(element,  fvGeometry, elemVolVars, scv)  * scv.volume() / Dumux::Components::Air<Scalar>::gasMolarDensity(refT_, fsRef.partialPressure(gasPhaseIdx, AirIdx));
                            advectiveWaterFlux_volume += advWaterFlux(element,  fvGeometry, elemVolVars, scv)  * scv.volume() / Dumux::Components::H2O<Scalar>::gasMolarDensity(refT_, fsRef.partialPressure(gasPhaseIdx, H2OIdx));
                            diffusiveAirFlux_volume += diffAirFlux(element,  fvGeometry, elemVolVars, scv)  * scv.volume() / Dumux::Components::Air<Scalar>::gasMolarDensity(refT_, fsRef.partialPressure(gasPhaseIdx, AirIdx));
                            diffusiveWaterFlux_volume += diffWaterFlux(element,  fvGeometry, elemVolVars, scv)  * scv.volume() / Dumux::Components::H2O<Scalar>::gasMolarDensity(refT_, fsRef.partialPressure(gasPhaseIdx, H2OIdx));
                        }
                    }
                }

#if !ISOTHERMAL
            evaporationFile << time_ << ", " << evaporation_mol << ", " << evaporation_mass << ", " << evaporation_volume << ", " << airFlux_mol << ", " << airFlux_volume << ", " << energyFlux << ", " << energyConductFlux << ", " << energyConvectFlux << ", " << advectiveWaterFlux_mol << ", " << diffusiveWaterFlux_mol << ", " << advectiveAirFlux_mol << ", " << diffusiveAirFlux_mol << ", " << advectiveWaterFlux_volume << ", " << diffusiveWaterFlux_volume << ", " << advectiveAirFlux_volume << ", " << diffusiveAirFlux_volume;
#else
            evaporationFile << time_ << ", " << evaporation_mol << ", " << evaporation_mass << ", " << evaporation_volume << ", " << airFlux_mol << ", " << airFlux_volume << ", " << advectiveWaterFlux_mol << ", " << diffusiveWaterFlux_mol << ", " << advectiveAirFlux_mol << ", " << diffusiveAirFlux_mol << ", " << advectiveWaterFlux_volume << ", " << diffusiveWaterFlux_volume << ", " << advectiveAirFlux_volume << ", " << diffusiveAirFlux_volume;
#endif


            for (const auto& element : elements(this->gridGeometry().gridView())){
                auto fvGeometry = localView(this->gridGeometry());
                auto elemVolVars = localView(gridVariables.curGridVolVars());
                auto elemFluxVarsCache = localView(gridVariables.gridFluxVarsCache());
                fvGeometry.bind(element);
                elemVolVars.bind(element, fvGeometry, curSol);
                elemFluxVarsCache.bind(element, fvGeometry, elemVolVars);
                for (auto&& scv : scvs(fvGeometry)){
                    if (scv.dofPosition()[1] > hmax - eps_){
                        auto fluxes = source(element,  fvGeometry, elemVolVars, scv);
                        evaporation_mol = fluxes[Indices::conti0EqIdx] * scv.volume();
                        evaporationFile << ", " << scv.dofPosition()[0] << ", " << evaporation_mol;
                    }
                }
            }
            evaporationFile << std::endl;
        }


    /*!
     * \brief Grow of salt in the pore throats based on the sourceterm
     * in the adjacent pore bodies (approach 2).
     * Only used for approach 2.
     */
#if USEBRINE
    template<class SolutionVector, class GridVariables>
    void calculateThroatPrecipitatedAreaAverageSource(const SolutionVector& sol, const GridVariables& gridVariables)
    {
         if(throatApproach_ == 2){
            auto fvGeometry = localView(this->gridGeometry());
            auto elemVolVars = localView(gridVariables.curGridVolVars());

            const auto numElements = this->gridGeometry().gridView().size(0);
            std::vector<Scalar> throatPrecipitatedAreaInTimeStep(numElements, 0.0);


            for(const auto& element:elements(this->gridGeometry().gridView())){
                fvGeometry.bind(element);
                elemVolVars.bind(element, fvGeometry, sol);
                const auto eIdx = this->gridGeometry().elementMapper().index(element);

                if (this->gridGeometry().useSameShapeForAllThroats()){

                        auto sourceterm = 0.0;
                        for (auto&& scv : scvs(fvGeometry)){
                            const Scalar dofIdx = scv.dofIndex();
                            const auto& volVars = elemVolVars[scv];
                            BoundaryTypes bcTypes = boundaryTypes(element, scv);
                            if(!bcTypes.isDirichlet(contiSolidEqIdx)){
                                if(useKineticPrecipitationReaction_){
                                    sourceterm += this->source(element, fvGeometry, elemVolVars, scv)[contiSolidEqIdx] * this->spatialParams().initialPoreVolume(dofIdx)/this->spatialParams().poreReactiveSurfaceArea(scv, volVars);
                                }else{
                                    Scalar precipVolFrac = sol[scv.localDofIndex()][solidNaClIdx];
                                    sourceterm += this->source(element, fvGeometry, elemVolVars, scv)[contiSolidEqIdx] / this->spatialParams().porosity(element, scv, precipVolFrac);
                                }
                            }
                        }

                        sourceterm  /= 2.0;

                        Scalar A = 0.0;
                        if(useKineticPrecipitationReaction_){
                            A = this->spatialParams().throatReactiveSurfaceArea(element, elemVolVars) / this->spatialParams().throatLength(element, elemVolVars); //m^2/m reative area per throat length
                        }else{
                            A = this->spatialParams().throatCrossSectionalArea(element, elemVolVars);
                        }

                        auto rho_solid = Dumux::Components::NaCl<Scalar>::solidMolarDensity(/*temperature=*/293); //TODO
                        throatPrecipitatedAreaInTimeStep[eIdx] = sourceterm * A * timeStepSize_ / rho_solid; //[m^3/ m length porethroat]
                }
                else {
                        DUNE_THROW(Dune::InvalidStateException, "Varying pore-throat geometries are not implemented for the mineralization model.");
                }
            }
            this->spatialParams().setThroatPrecipitatedAreaAverageSource(throatPrecipitatedAreaInTimeStep);
         }
         else{
             DUNE_THROW(Dune::InvalidStateException, "Function calculateSourceApprochPrecipitation just used for the averaged source approach (approach 2).");
         }
    }
#endif

     /*!
     * \name Boundary conditions
     */
    // \{
    //! Specifies which kind of boundary condition should be used for
    //! which equation for a finite volume on the boundary.
    BoundaryTypes boundaryTypes(const Element& element, const SubControlVolume& scv) const
    {
        BoundaryTypes bcTypes;
        bcTypes.setAllNeumann();
        if(useDirichletAtBottom_ && scv.dofPosition()[1] < eps_){
            bcTypes.setAllDirichlet();
        }
        return bcTypes;
    }


    //! Evaluate the boundary conditions for a Dirichlet control volume.
    PrimaryVariables dirichlet(const Element& element,
                               const SubControlVolume& scv) const
    {
        PrimaryVariables values(0.0);
        values = initialConditions_[scv.dofIndex()];
        return values;
    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{

    //! Evaluate the source term for all phases within a given sub-control-volume.
    PrimaryVariables source(const Element& element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const SubControlVolume& scv) const
    {
        PrimaryVariables values(0.0);

        ////////////////////////////////////////////
        // EVAPORATION SOURCE TERM
        ////////////////////////////////////////////

        const Scalar hmax = this->gridGeometry().bBoxMax()[1];

        if(scv.dofPosition()[1] > hmax - eps_){

            Scalar advectiveAirFlux = advAirFlux(element, fvGeometry, elemVolVars, scv);
            Scalar advectiveWaterFlux = 0.0;
            Scalar diffusiveAirFlux = 0.0;
            Scalar diffusiveWaterFlux = 0.0;
#if !ISOTHERMAL
            Scalar energyConduction = 0.0;
            Scalar energyConvection = 0.0;
#endif

            if(time_ > episodeLength_){
                advectiveWaterFlux = advWaterFlux(element, fvGeometry, elemVolVars, scv);
                diffusiveAirFlux = diffAirFlux(element, fvGeometry, elemVolVars, scv);
                diffusiveWaterFlux = diffWaterFlux(element, fvGeometry, elemVolVars, scv);
#if !ISOTHERMAL
                energyConduction = energyConductionFlux(element, fvGeometry, elemVolVars, scv);
                energyConvection = energyConvectionFlux(element, fvGeometry, elemVolVars, scv);
#endif
            }

            values[conti0EqIdx] = advectiveWaterFlux + diffusiveWaterFlux;
            values[conti1EqIdx] = advectiveAirFlux + diffusiveAirFlux;
#if !ISOTHERMAL
            values[energyEqIdx] = energyConduction + energyConvection;
#endif
        }

#if USEBRINE
        ////////////////////////////////////////////
        // PRECIPITATION SOURCE TERM
        ////////////////////////////////////////////
        Scalar saltSource = 0.0;

        const auto& volVars = elemVolVars[scv];
        const Scalar dofIdx = scv.dofIndex();

        // some parameters for convenience
        Scalar molFraction = volVars.moleFraction(liquidPhaseIdx, NaClIdx); //actual mole fraction of salt in the liquid phase
        Scalar liquidSaturation = volVars.saturation(liquidPhaseIdx);
        Scalar precipVolFrac = volVars.solidVolumeFraction(NaClSolidSystemIdx); //solidVolumeFraction(0) is the solid volume of NaCl per volume of the respective scv, which is the initial pore-body volume in case of a pnm
        Scalar porosity = this->spatialParams().porosity(element, scv, precipVolFrac);
        Scalar solubilityMolFraction = 0.0996; //corresponds to molality = 6.14 mol NaCl/kgH2O

        if (useKineticPrecipitationReaction_){
            auto molality = molFraction / ((1-molFraction) * FluidSystem::molarMass(H2OIdx)); //mol NaCl/kg H20
            auto ionicStrength = molality; //assumption that molarität Na = molarität Cl
            auto saturationIndex = 0.0;
            if (activityCoefficientPitzerApproach_){
                auto equilibriumSolubilityProduct = 38.0464; //mol²/kg² 38.04643300531410 corrsponds to mole fraction 0.0996 and a molality of 6.14 mol NaCl/ kg H2O (38.0410997957714 corresponds exactly to molality 6.14 mol NaCl/kg H2O)
                auto activityCoefficientNaCl = PitzerNaCl<Scalar>::activityCoefficient(molality);
                saturationIndex = activityCoefficientNaCl * molality * activityCoefficientNaCl * molality/ equilibriumSolubilityProduct;
            }else{
                auto equilibriumSolubilityProduct = 32.5695; //mol²/kg² 32.5694616843 corrsponds to mole fraction 0.0996
                auto activityCoefficientNa = std::exp(-0.5085 * std::sqrt(ionicStrength) / (1 + 3.285e9 * 4e-10 * std::sqrt(ionicStrength)) + 0.075 * ionicStrength);
                auto activityCoefficientCl = std::exp(-0.5085 * std::sqrt(ionicStrength) / (1 + 3.285e9 * 3.5e-10 * std::sqrt(ionicStrength)) + 0.015 * ionicStrength);
                saturationIndex = activityCoefficientNa * molality * activityCoefficientCl * molality/ equilibriumSolubilityProduct;
            }
//             auto reactiveSurfaceArea = 4 * M_PI * this->spatialParams().poreInscribedRadius(dofIdx, volVars) * this->spatialParams().poreInscribedRadius(dofIdx, volVars);
            //TODO: is just for cubic bodies
            auto reactiveSurfaceArea = this->spatialParams().poreReactiveSurfaceArea(scv, volVars);
            auto reactionRate = 1e-3 * reactiveSurfaceArea * (saturationIndex - 1.0); //mol/s
            auto precipitingSaltRate = 0.0;

            Scalar initPoreVolume = this ->spatialParams().initialPoreVolume(dofIdx);

            //check that not more than the available amount of mol in the solid phase dissolves
            if (saturationIndex < 1.0){
                Scalar mol = - reactionRate * timeStepSize_;//[molNaCl] total number of molecules which would dissolve in the pore body so that the solubility limit is reached
                Scalar availableMol = precipVolFrac * initPoreVolume * Dumux::Components::NaCl<Scalar>::solidMolarDensity(/*temperature=*/293);//[molNaCl] available number of molecules for dissolution in the solid phase
                precipitingSaltRate = - std::min(mol, availableMol) / (timeStepSize_ * initPoreVolume);//[molNaCl/m^3]
            }
            else{ //check that not more salt is precipitating than there is available in the solution
                auto availableMolInSolution = (molFraction - solubilityMolFraction) * liquidSaturation * this->spatialParams().poreVolume(scv, volVars) * volVars.molarDensity(liquidPhaseIdx);
                auto potentialMolPrecipitate = reactionRate * timeStepSize_;
                precipitingSaltRate = std::min(availableMolInSolution, potentialMolPrecipitate) / (timeStepSize_ * initPoreVolume);
            }


            values[contiNaClEqIdx] += -precipitingSaltRate;
            values[contiSolidEqIdx] += precipitingSaltRate;

            // ---------------------------------------------------------------
            // Low-saturation regularization of the dissolved-salt equation.
            // As the brine vanishes (sw -> 0), both the dissolved-NaCl storage
            // (~ sw) and the precipitation-source coupling (~ sw) go to zero,
            // so the Jacobian row/column for the dissolved mole fraction x_NaCl
            // degenerates and the linear system becomes singular at the
            // bothPhases -> gas-only transition. We restore a non-vanishing
            // diagonal by relaxing x_NaCl towards the solubility limit (the
            // physical state of an evaporating brine film) with a strength that
            // is only active below swSaltReg and matches the storage scale.
            static const Scalar swSaltReg = getParam<Scalar>("Problem.SaltRegularizationSaturation", 1e-3);
            static const Scalar regStrength = getParam<Scalar>("Problem.SaltRegularizationStrength", 1.0);
            if (liquidSaturation < swSaltReg)
            {
                const Scalar w = (swSaltReg - liquidSaturation) / swSaltReg; // 1 at sw=0, 0 at sw=swSaltReg
                const Scalar regRate = regStrength * volVars.molarDensity(liquidPhaseIdx) * swSaltReg / timeStepSize_; //[molNaCl/(m^3 s)]
                const Scalar regTerm = regRate * w * (solubilityMolFraction - molFraction); // <0 if supersaturated -> precipitate
                values[contiNaClEqIdx] += regTerm;        // pin dissolved salt to solubility as brine dries
                values[contiSolidEqIdx] -= regTerm;       // conserve total salt (dissolved <-> solid)
            }

        }else{ //equilibrium reaction

            //precipitates/ dissolves so that equilibrium state is reached
            saltSource = porosity * liquidSaturation * volVars.molarDensity(liquidPhaseIdx) * (molFraction - solubilityMolFraction);// [molNaCl/m^3] positiv for precipitation, negativ for dissolution

            //check that not more than the available amount of mol in the solid phase dissolves
            if (molFraction < solubilityMolFraction){

                Scalar initPoreVolume = this ->spatialParams().initialPoreVolume(dofIdx);

                Scalar mol = - saltSource * initPoreVolume;//[molNaCl] total number of molecules which would dissolve in the pore body so that the solubility limit is reached
                Scalar availableMol = precipVolFrac * initPoreVolume * Dumux::Components::NaCl<Scalar>::solidMolarDensity(/*temperature=*/293);//[molNaCl] available number of molecules for dissolution in the solid phase (solid molar density is constant and not temperature dependent)
                saltSource = - std::min(mol, availableMol) / initPoreVolume;//[molNaCl/m^3]
            }

            //TODO handle case of fully blocked pore-bodies

            const Scalar time = timeStepSize_;
            saltSource /= time;//[molNaCl/(m^3*s)]

            values[contiNaClEqIdx] += -saltSource;
            values[contiSolidEqIdx] += saltSource;
        }
#endif


        return values;
    }
    // \}


    ///////////////////////////////////////////////////////////////////////
    // Phase fluxes
    //////////////////////////////////////////////////////////////////////

    // advective flux of the air component in the gas phase
    Scalar advAirFlux(const Element& element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const SubControlVolume& scv) const
    {
        Scalar advectiveAirFlux = 0.0;
        const auto& volVars = elemVolVars[scv];

        //Free flow fluidstate
        FluidState fsRef;
        fsRef.setTemperature(refT_);
        fsRef.setSaturation(gasPhaseIdx, 1.0);
        fsRef.setPressure(gasPhaseIdx, refGasPressure_);
        fsRef.setMoleFraction(gasPhaseIdx, H2OIdx, refWaterVaporMoleFrac_);
        fsRef.setMoleFraction(gasPhaseIdx, AirIdx, 1-refWaterVaporMoleFrac_);

        //Freeflow properties
        const Scalar molarGasDensityRef = FluidSystem::molarDensity(fsRef, gasPhaseIdx);
        const Scalar molarGasViscosityRef = FluidSystem::viscosity(fsRef, gasPhaseIdx);

        // gas flows out
        if (volVars.pressure(gasPhaseIdx) - refGasPressure_ > 0) {
            advectiveAirFlux = -1.0 * transmissibilityFactor_ * (volVars.pressure(gasPhaseIdx) - refGasPressure_)
                                    /(boundaryLayerThickness_ * volVars.viscosity(gasPhaseIdx))
                                    *volVars.molarDensity(gasPhaseIdx)
                                    *volVars.moleFraction(gasPhaseIdx, AirIdx);
        }
        //gas flows in
        else {
            advectiveAirFlux = -1.0 * transmissibilityFactor_ * (volVars.pressure(gasPhaseIdx) - refGasPressure_)
                                    /(boundaryLayerThickness_ * molarGasViscosityRef)
                                    *molarGasDensityRef * (1-refWaterVaporMoleFrac_);
        }
        advectiveAirFlux /= scv.volume();
        return advectiveAirFlux;
    }


    //Advective transport of water component in the gas phase
    Scalar advWaterFlux(const Element& element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const SubControlVolume& scv) const
    {
        Scalar advectiveWaterFlux = 0.0;
        const auto& volVars = elemVolVars[scv];

        //Free flow fluidstate
        FluidState fsRef;
        fsRef.setTemperature(refT_);
        fsRef.setSaturation(gasPhaseIdx, 1.0);
        fsRef.setPressure(gasPhaseIdx, refGasPressure_);
        fsRef.setMoleFraction(gasPhaseIdx, H2OIdx, refWaterVaporMoleFrac_);
        fsRef.setMoleFraction(gasPhaseIdx, AirIdx, 1-refWaterVaporMoleFrac_);

        //Freeflow properties
        const Scalar molarGasDensityRef = FluidSystem::molarDensity(fsRef, gasPhaseIdx);
        const Scalar molarGasViscosityRef = FluidSystem::viscosity(fsRef, gasPhaseIdx);

        if (volVars.pressure(gasPhaseIdx) - refGasPressure_ > 0) {
            advectiveWaterFlux = -1.0 * transmissibilityFactor_ * (volVars.pressure(gasPhaseIdx) - refGasPressure_)
                                    /(boundaryLayerThickness_ * volVars.viscosity(gasPhaseIdx))
                                    *volVars.molarDensity(gasPhaseIdx)
                                    *volVars.moleFraction(gasPhaseIdx, H2OIdx);
        }
        //gas flows in
        else {
            advectiveWaterFlux = -1.0 * transmissibilityFactor_ * (volVars.pressure(gasPhaseIdx) - refGasPressure_)
                                    /(boundaryLayerThickness_ * molarGasViscosityRef)
                                    *molarGasDensityRef * refWaterVaporMoleFrac_;
        }
        advectiveWaterFlux /= scv.volume();
        return advectiveWaterFlux;
    }

    //Diffusion of air component in gas phase
    Scalar diffAirFlux(const Element& element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const SubControlVolume& scv) const
    {
        Scalar diffusiveAirFlux = 0.0;
        const auto& volVars = elemVolVars[scv];

        //Free flow fluidstate
        FluidState fsRef;
        fsRef.setTemperature(refT_);
        fsRef.setSaturation(gasPhaseIdx, 1.0);
        fsRef.setPressure(gasPhaseIdx, refGasPressure_);
        fsRef.setMoleFraction(gasPhaseIdx, H2OIdx, refWaterVaporMoleFrac_);
        fsRef.setMoleFraction(gasPhaseIdx, AirIdx, 1-refWaterVaporMoleFrac_);

        //Freeflow properties
        const Scalar massGasDensityRef = FluidSystem::density(fsRef, gasPhaseIdx);
        const Scalar DRef = FluidSystem::binaryDiffusionCoefficient(fsRef, gasPhaseIdx, H2OIdx, AirIdx);

        const Scalar poreCrossSectionalArea = this ->spatialParams().poreCrossSectionalArea(scv, volVars);//M_PI * 0.0002 *0.0002;//

        //Diffusive transport
        Scalar DAvg = 0.5*(volVars.diffusionCoefficient(gasPhaseIdx, H2OIdx, AirIdx) + DRef);
        const Scalar gasDensityAvg = 0.5*(volVars.density(gasPhaseIdx) + massGasDensityRef);
        Scalar airMolarMass = Dumux::Components::Air<Scalar>::molarMass();
        Scalar massFracAir = volVars.massFraction(gasPhaseIdx, AirIdx);

        diffusiveAirFlux = - gasDensityAvg / airMolarMass * (massFracAir - (1-refWaterVaporMassFrac_)) / boundaryLayerThickness_ * DAvg * poreCrossSectionalArea;

        diffusiveAirFlux /= scv.volume();
        return diffusiveAirFlux;
    }

    //Diffusion of water component in gas phase
    Scalar diffWaterFlux(const Element& element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const SubControlVolume& scv) const
    {
        Scalar diffusiveWaterFlux = 0.0;
        const auto& volVars = elemVolVars[scv];

        //Free flow fluidstate
        FluidState fsRef;
        fsRef.setTemperature(refT_);
        fsRef.setSaturation(gasPhaseIdx, 1.0);
        fsRef.setPressure(gasPhaseIdx, refGasPressure_);
        fsRef.setMoleFraction(gasPhaseIdx, H2OIdx, refWaterVaporMoleFrac_);
        fsRef.setMoleFraction(gasPhaseIdx, AirIdx, 1-refWaterVaporMoleFrac_);

        //Freeflow properties
        const Scalar massGasDensityRef = FluidSystem::density(fsRef, gasPhaseIdx);
        const Scalar DRef = FluidSystem::binaryDiffusionCoefficient(fsRef, gasPhaseIdx, H2OIdx, AirIdx);

        const Scalar poreCrossSectionalArea = this ->spatialParams().poreCrossSectionalArea(scv, volVars);//M_PI * 0.0002 *0.0002;//

        //Diffusive transport
        Scalar DAvg = 0.5*(volVars.diffusionCoefficient(gasPhaseIdx, H2OIdx, AirIdx) + DRef);
        const Scalar gasDensityAvg = 0.5*(volVars.density(gasPhaseIdx) + massGasDensityRef);
        Scalar massFracWaterVapor = volVars.massFraction(gasPhaseIdx, H2OIdx);
        Scalar waterMolarMass = Dumux::Components::H2O<Scalar>::molarMass();

        diffusiveWaterFlux = - gasDensityAvg / waterMolarMass * (massFracWaterVapor - refWaterVaporMassFrac_) / boundaryLayerThickness_ * DAvg * poreCrossSectionalArea;

        diffusiveWaterFlux /= scv.volume();
        return diffusiveWaterFlux;
    }


    ///////////////////////////////////////////////////////////////////////
    // Energy fluxes
    //////////////////////////////////////////////////////////////////////
#if !ISOTHERMAL
    Scalar energyConductionFlux(const Element& element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const SubControlVolume& scv) const
    {
        Scalar conductionFlux = 0.0;
        const auto& volVars = elemVolVars[scv];

        //Free flow fluidstate
        FluidState fsRef;
        fsRef.setTemperature(refT_);
        fsRef.setSaturation(gasPhaseIdx, 1.0);
        fsRef.setPressure(gasPhaseIdx, refGasPressure_);
        fsRef.setMoleFraction(gasPhaseIdx, H2OIdx, refWaterVaporMoleFrac_);
        fsRef.setMoleFraction(gasPhaseIdx, AirIdx, 1-refWaterVaporMoleFrac_);

        //Energy flux
        const Scalar poreCrossSectionalArea = this ->spatialParams().poreCrossSectionalArea(scv, volVars);
        Scalar lambdaGas = FluidSystem::thermalConductivity(volVars.fluidState(), gasPhaseIdx);
        Scalar T = volVars.temperature();
        Scalar lambdaGasRef = FluidSystem::thermalConductivity(fsRef, gasPhaseIdx);

        //heat conduction
        conductionFlux = 0.5 * (lambdaGas + lambdaGasRef) * (refT_ - T)/boundaryLayerThickness_ * poreCrossSectionalArea/scv.volume();
        return conductionFlux;
    }


    Scalar energyConvectionFlux(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const SubControlVolume& scv) const
    {
        Scalar convectiveFlux = 0.0;
        const auto& volVars = elemVolVars[scv];

        //Free flow fluidstate
        FluidState fsRef;
        fsRef.setTemperature(refT_);
        fsRef.setSaturation(gasPhaseIdx, 1.0);
        fsRef.setPressure(gasPhaseIdx, refGasPressure_);
        fsRef.setMoleFraction(gasPhaseIdx, H2OIdx, refWaterVaporMoleFrac_);
        fsRef.setMoleFraction(gasPhaseIdx, AirIdx, 1-refWaterVaporMoleFrac_);

        Scalar waterMolarMass = Dumux::Components::H2O<Scalar>::molarMass();
        Scalar airMolarMass = Dumux::Components::Air<Scalar>::molarMass();

        Scalar hgw = FluidSystem::componentEnthalpy(volVars.fluidState(), gasPhaseIdx, H2OIdx);
        Scalar hga = FluidSystem::componentEnthalpy(volVars.fluidState(), gasPhaseIdx, AirIdx);
//         Scalar T = volVars.temperature();

        Scalar hgwRef = FluidSystem::componentEnthalpy(fsRef, gasPhaseIdx, H2OIdx);
        Scalar hgaRef = FluidSystem::componentEnthalpy(fsRef, gasPhaseIdx, AirIdx);

        Scalar diffusiveWaterFlux = diffWaterFlux(element, fvGeometry, elemVolVars, scv);
        Scalar diffusiveAirFlux = diffAirFlux(element, fvGeometry, elemVolVars, scv);
        Scalar advectiveWaterFlux = advWaterFlux(element, fvGeometry, elemVolVars, scv);
        Scalar advectiveAirFlux = advAirFlux(element, fvGeometry, elemVolVars, scv);

        //diffusive energy transport
        if(diffusiveWaterFlux < 0){ //outflow
            convectiveFlux += diffusiveWaterFlux * hgw * waterMolarMass;
        }
        else{ //inflow
            convectiveFlux += diffusiveWaterFlux * hgwRef * waterMolarMass;
        }

        if(diffusiveAirFlux < 0){ //outflow
            convectiveFlux += diffusiveAirFlux * hga * airMolarMass;
        }
        else{ //inflow
            convectiveFlux += diffusiveAirFlux * hgaRef * airMolarMass;
        }


        //advective energy transport
        if (volVars.pressure(gasPhaseIdx) - refGasPressure_ > 0) { //outflow
            convectiveFlux += advectiveWaterFlux * hgw * waterMolarMass + advectiveAirFlux * hga * airMolarMass;
        }
        else{ //inflow
            convectiveFlux += advectiveWaterFlux * hgwRef * waterMolarMass + advectiveAirFlux * hgaRef * airMolarMass;
        }

        return convectiveFlux;
    }
#endif // !ISOTHERMAL



    //! Evaluate the initial value for a control volume.
    PrimaryVariables initial(const Vertex& vertex) const
    {
        return initialConditions_[this->gridGeometry().vertexMapper().index(vertex)];
    }

    void getInitial()
    {

        initialConditions_.resize(this->gridGeometry().numDofs());

        PrimaryVariables priVars(0.0);
        priVars.setState(Indices::bothPhases);

        auto fvGeometry = localView(this->gridGeometry());
        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            fvGeometry.bind(element);
            for (const auto& scv : scvs(fvGeometry)){
                priVars[pressureIdx] = initialPressure_;
                // set water saturation
                if (scv.dofPosition()[dimWorld-1] > this->gridGeometry().bBoxMax()[dimWorld-1] - eps_)
                {
                    priVars[switchIdx] = initialSaturationTop_;
                }
                else if (scv.dofPosition()[dimWorld-1] < this->gridGeometry().bBoxMin()[dimWorld-1] + eps_)
                {
                    priVars[switchIdx] =  initialSaturationBottom_;
                }
                else{
                    priVars[switchIdx] = initialSaturation_;//Sw for bothPhases
                }

#if USEBRINE
                priVars[NaClIdx] = initialSaltMoleFraction_;
                priVars[solidNaClIdx] = initialSolidVolFrac_;
#endif

#if !ISOTHERMAL
                priVars[temperatureIdx] = initialTemperature_;
#endif

                initialConditions_[scv.dofIndex()] = priVars;
            }
        }
    }

    //!  Evaluate the initial invasion state of a pore throat
    bool initialInvasionState(const Element& element) const
    {
        const auto eIdx = this->gridGeometry().elementMapper().index(element);
        if (this->gridGeometry().throatLabel(eIdx) == 3){
            return true;
        }else{
            return false;
        }
    }

    // \}

    // throat parameter indicating invasion state
    template<class FluxVariablesCache, class SubControlVolumeFace>
    Scalar theta(const Element& element,
                 const FVElementGeometry& fvGeometry,
                 const ElementVolumeVariables& elemVolVars,
                 const FluxVariablesCache& fluxVarsCache,
                 const SubControlVolumeFace& scvf) const
    {
        const auto invaded = fluxVarsCache.invaded();
        if constexpr (useThetaRegularization)
        {
            if(!invaded)
            {
                using std::max;
                auto pcEntry = this->spatialParams().pcEntry(element, elemVolVars);
                auto dp = max(elemVolVars[0].capillaryPressure(),
                                elemVolVars[1].capillaryPressure()) / pcEntry - 1.0;
                // Use a regularization function for theta
                return reg_.eval(dp,invaded);
            }
            else
            {
                using std::min; using std::abs;
                auto pcSnapoff = this->spatialParams().pcSnapoff(element, elemVolVars);
                auto dp = min(elemVolVars[0].capillaryPressure(),
                                elemVolVars[1].capillaryPressure()) / abs(pcSnapoff) - sign(pcSnapoff);
                // Use a regularization function for theta
                return reg_.eval(dp,invaded);
            }
        }
        else
        {
            return invaded ? 1.0 : 0.0;
        }
    }

    void setTime(const Scalar t)
    { time_ = t; }

    void setTimeStepSize(const Scalar t)
    { timeStepSize_ = t; }

private:

    bool isInletPore_(const SubControlVolume& scv) const
    {
        return isInletPore_(scv.dofIndex());
    }

    bool isInletPore_(const std::size_t dofIdxGlobal) const
    {
        return this->gridGeometry().poreLabel(dofIdxGlobal) == Labels::inlet;
    }

    bool isOutletPore_(const SubControlVolume& scv) const
    {
        return this->gridGeometry().poreLabel(scv.dofIndex()) == Labels::outlet;
    }

    //Conversion functions for mass and mole fractions
    static const Scalar massToMoleFrac_(Scalar X, Scalar phaseIdx, Scalar compIdx, Scalar Xaliq = 0.0)
    {
        if (phaseIdx == gasPhaseIdx && compIdx == H2OIdx){
            const Scalar Mw = FluidSystem::molarMass(H2OIdx); /* molecular weight of water [kg/mol] */
            const Scalar Ma = FluidSystem::molarMass(AirIdx); /* molecular weight of air  [kg/mol] */

            /* XwNaCl: conversion from mass fraction to mol fraction */
            auto xw = Ma * X / ((Ma - Mw) * X + Mw);
            return xw;
        }

#if USEBRINE
        else if (phaseIdx == liquidPhaseIdx && compIdx == NaClIdx){
                const Scalar Mw = FluidSystem::molarMass(H2OIdx); /* molecular weight of water [kg/mol] */
                const Scalar MNaCl = FluidSystem::molarMass(NaClIdx); /* molecular weight of NaCl  [kg/mol] */
                const Scalar Ma = FluidSystem::molarMass(AirIdx); /* molecular weight of air  [kg/mol] */

                /* XwNaCl: conversion from mass fraction to mol fraction */
                auto xNaCl = (X / MNaCl)  / (X / MNaCl + (1 - X - Xaliq) / Mw + Xaliq / Ma);
                return xNaCl;
        }
#endif
        else {
                DUNE_THROW(Dune::InvalidStateException, "Mass to mole fraction conversion not implemented for this phaseIdx compIdx combination.");
        }
    }

    static const Scalar moleToMassFrac_(Scalar x, Scalar phaseIdx, Scalar compIdx, Scalar xaliq = 0.0)
    {
        if (phaseIdx == gasPhaseIdx && compIdx == H2OIdx){
            const Scalar Mw = FluidSystem::molarMass(H2OIdx); /* molecular weight of water [kg/mol] */
            const Scalar Ma = FluidSystem::molarMass(AirIdx); /* molecular weight of air  [kg/mol] */

            /* XwNaCl: conversion from mass fraction to mol fraction */
            auto Xw = Mw * x / ((Mw - Ma) * x + Ma);
            return Xw;
        }
#if USEBRINE
        else if (phaseIdx == liquidPhaseIdx && compIdx == NaClIdx){
            const Scalar Mw = FluidSystem::molarMass(H2OIdx); /* molecular weight of water [kg/mol] */
            const Scalar MNaCl = FluidSystem::molarMass(NaClIdx); /* molecular weight of NaCl  [kg/mol] */
            const Scalar Ma = FluidSystem::molarMass(AirIdx); /* molecular weight of air  [kg/mol] */

            /* XwNaCl: conversion from mass fraction to mol fraction */
            auto XNaCl = MNaCl * x / (x * MNaCl + (1 - x - xaliq) * Mw + xaliq * Ma);
            return XNaCl;
        }
#endif
        else {
            DUNE_THROW(Dune::InvalidStateException, "Mole to mass fraction conversion not implemented for this phaseIdx compIdx combination.");
        }
    }



    Scalar episodeLength_;
    Scalar time_ = 0.0;
    Scalar timeStepSize_ = 0.0;
    Scalar eps_ = 1e-6;
    bool useKineticPrecipitationReaction_;
    bool activityCoefficientPitzerApproach_;

    //io
    bool useSmallGrid_;
    int vtpOutputFrequency_;

    std::vector<PrimaryVariables> initialConditions_;
    Scalar initialSaturation_;
    Scalar initialSaturationBottom_;
    Scalar initialSaturationTop_;
    Scalar initialPressure_;
#if USEBRINE
    Scalar initialSaltMoleFraction_;
    Scalar initialSolidVolFrac_;
    Scalar throatApproach_;
#endif

    Scalar boundaryLayerThickness_;
    Scalar refT_;
    Scalar refWaterVaporMoleFrac_;
    Scalar refWaterVaporMassFrac_;
    Scalar refGasPressure_;
    Scalar transmissibilityFactor_;
    bool useDirichletAtBottom_;
    Dumux::PoreNetwork::Throat::Regularization<Scalar> reg_;



#if !ISOTHERMAL
    Scalar initialTemperature_;
#endif
};
} //end namespace Dumux

#endif
