// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 *
 * \brief Test for the two-phase n-components pore-network model
 */
#include <config.h>

#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>
#include <dune/grid/io/file/dgfparser/dgfexception.hh>
#include <dune/grid/io/file/vtk.hh>

#include <dumux/assembly/fvassembler.hh>
#include <dumux/common/initialize.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/io/grid/porenetwork/gridmanager.hh>
#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/porenetwork/common/pnmvtkoutputmodule.hh>
#include <dumux/nonlinear/newtonsolver.hh>
#include <dumux/porenetwork/2p/newtonsolver.hh>

#include "properties.hh"

int main(int argc, char** argv)
{
    using namespace Dumux;

    using TypeTag = Properties::TTag::DrainageProblem;

    // maybe initialize MPI and/or multithreading backend
    Dumux::initialize(argc, argv);
    const auto& mpiHelper = Dune::MPIHelper::instance();

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    ////////////////////////////////////////////////////////////
    // parse the command line arguments and input file
    ////////////////////////////////////////////////////////////

    // parse command line arguments
    Parameters::init(argc, argv);

    //////////////////////////////////////////////////////////////////////
    // try to create a grid (from the given grid file or the input file)
    /////////////////////////////////////////////////////////////////////

    using GridManager = Dumux::PoreNetwork::GridManager<3>;
    GridManager gridManager;
    gridManager.init();

    // we compute on the leaf grid view
    const auto& leafGridView = gridManager.grid().leafGridView();
    auto gridData = gridManager.getGridData();

    // create the finite volume grid geometry
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(leafGridView, *gridData);

    // the spatial parameters
    using SpatialParams = GetPropType<TypeTag, Properties::SpatialParams>;
    auto spatialParams = std::make_shared<SpatialParams>(gridGeometry);

    // the problem (boundary conditions)
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry, spatialParams);

    // the solution vector
    using GridView = typename GridGeometry::GridView;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector x(leafGridView.size(GridView::dimension));
    problem->applyInitialSolution(x);
    auto xOld = x;

    // the grid variables
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);

    // Validation-only diagnostics.  The standard VTP writer stores these small pressure drops
    // as Float32 and writes phase flux magnitudes.  Keep the governing equations untouched and
    // additionally write full-precision state and signed throat/component fluxes.
    using FluxVariables = GetPropType<TypeTag, Properties::FluxVariables>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    static constexpr int liquidPhaseIdx = FluidSystem::liquidPhaseIdx;
    static constexpr int waterCompIdx = FluidSystem::H2OIdx;
    static constexpr int saltCompIdx = FluidSystem::NaClIdx;
    std::ofstream validationPores("validation_pores.csv");
    std::ofstream validationThroats("validation_throats.csv");
    validationPores << std::scientific << std::setprecision(17);
    validationThroats << std::scientific << std::setprecision(17);
    validationPores << "time_s,pore_index,x_m,y_m,z_m,p_liq_pa,s_liq,x_nacl_liq,"
                     << "solid_volume_fraction,pore_volume_m3,water_inventory_mol,"
                     << "dissolved_nacl_inventory_mol,solid_nacl_inventory_mol\n";
    validationThroats << "time_s,throat_index,pore_a,pore_b,q_liq_a_to_b_m3_s,"
                       << "transmissibility_liq_m3_pa_s,diff_water_mol_s,diff_nacl_mol_s,"
                       << "adv_water_mol_s,adv_nacl_mol_s,total_water_mol_s,"
                       << "total_nacl_mol_s,total_water_all_phases_mol_s,"
                       << "total_nacl_all_phases_mol_s,throat_radius_m\n";

    auto writeValidationDiagnostics = [&](const Scalar time)
    {
        auto fvGeometry = localView(*gridGeometry);
        auto elemVolVars = localView(gridVariables->curGridVolVars());
        auto elemFluxVarsCache = localView(gridVariables->gridFluxVarsCache());
        std::vector<bool> poreWritten(leafGridView.size(GridView::dimension), false);

        for (const auto& element : elements(leafGridView, Dune::Partitions::interior))
        {
            fvGeometry.bindElement(element);
            elemVolVars.bind(element, fvGeometry, x);
            elemFluxVarsCache.bind(element, fvGeometry, elemVolVars);

            for (const auto& scv : scvs(fvGeometry))
            {
                const auto dofIdx = scv.dofIndex();
                if (poreWritten[dofIdx])
                    continue;
                poreWritten[dofIdx] = true;
                const auto& v = elemVolVars[scv];
                const auto& pos = scv.dofPosition();
                Scalar waterInventory = 0.0;
                Scalar dissolvedSaltInventory = 0.0;
                for (int phaseIdx = 0; phaseIdx < FluidSystem::numPhases; ++phaseIdx)
                {
                    const auto phaseMoles = scv.volume() * v.porosity()
                                            * v.saturation(phaseIdx)
                                            * v.molarDensity(phaseIdx);
                    waterInventory += phaseMoles * v.moleFraction(phaseIdx, waterCompIdx);
                    dissolvedSaltInventory += phaseMoles * v.moleFraction(phaseIdx, saltCompIdx);
                }
                const auto solidSaltInventory = scv.volume() * v.solidVolumeFraction(0)
                                                * v.solidComponentMolarDensity(0);
                validationPores << time << ',' << dofIdx << ','
                                << pos[0] << ',' << pos[1] << ',' << pos[2] << ','
                                << v.pressure(liquidPhaseIdx) << ','
                                << v.saturation(liquidPhaseIdx) << ','
                                << v.moleFraction(liquidPhaseIdx, saltCompIdx) << ','
                                << v.solidVolumeFraction(0) << ',' << v.poreVolume() << ','
                                << waterInventory << ',' << dissolvedSaltInventory << ','
                                << solidSaltInventory << '\n';
            }

            for (const auto& scvf : scvfs(fvGeometry))
            {
                if (scvf.boundary())
                    continue;
                FluxVariables fluxVars;
                fluxVars.init(*problem, element, fvGeometry, elemVolVars, scvf,
                              elemFluxVarsCache);
                const auto& cache = elemFluxVarsCache[scvf];
                const auto& inside = fvGeometry.scv(scvf.insideScvIdx());
                const auto& outside = fvGeometry.scv(scvf.outsideScvIdx());
                auto volumeUpwind = [](const auto& v)
                { return v.mobility(liquidPhaseIdx); };
                auto waterUpwind = [](const auto& v)
                {
                    return v.molarDensity(liquidPhaseIdx)
                           * v.moleFraction(liquidPhaseIdx, waterCompIdx)
                           * v.mobility(liquidPhaseIdx);
                };
                auto saltUpwind = [](const auto& v)
                {
                    return v.molarDensity(liquidPhaseIdx)
                           * v.moleFraction(liquidPhaseIdx, saltCompIdx)
                           * v.mobility(liquidPhaseIdx);
                };
                const auto q = fluxVars.advectiveFlux(liquidPhaseIdx, volumeUpwind);
                const auto advWater = fluxVars.advectiveFlux(liquidPhaseIdx, waterUpwind);
                const auto advSalt = fluxVars.advectiveFlux(liquidPhaseIdx, saltUpwind);
                const auto diff = fluxVars.molecularDiffusionFlux(liquidPhaseIdx);
                static constexpr auto referenceSystem =
                    FluxVariables::MolecularDiffusionType::referenceSystemFormulation();
                auto molarDiffusion = [](const auto rawFlux, const int compIdx)
                {
                    if constexpr (referenceSystem == ReferenceSystemFormulation::massAveraged)
                        return rawFlux / FluidSystem::molarMass(compIdx);
                    else
                        return rawFlux;
                };
                const auto diffWater = molarDiffusion(diff[waterCompIdx], waterCompIdx);
                const auto diffSalt = molarDiffusion(diff[saltCompIdx], saltCompIdx);
                Scalar totalWaterAllPhases = 0.0;
                Scalar totalSaltAllPhases = 0.0;
                for (int phaseIdx = 0; phaseIdx < FluidSystem::numPhases; ++phaseIdx)
                {
                    auto waterPhaseUpwind = [phaseIdx](const auto& v)
                    {
                        return v.molarDensity(phaseIdx)
                               * v.moleFraction(phaseIdx, waterCompIdx)
                               * v.mobility(phaseIdx);
                    };
                    auto saltPhaseUpwind = [phaseIdx](const auto& v)
                    {
                        return v.molarDensity(phaseIdx)
                               * v.moleFraction(phaseIdx, saltCompIdx)
                               * v.mobility(phaseIdx);
                    };
                    const auto phaseDiff = fluxVars.molecularDiffusionFlux(phaseIdx);
                    totalWaterAllPhases += fluxVars.advectiveFlux(phaseIdx, waterPhaseUpwind)
                                           + molarDiffusion(phaseDiff[waterCompIdx], waterCompIdx);
                    totalSaltAllPhases += fluxVars.advectiveFlux(phaseIdx, saltPhaseUpwind)
                                          + molarDiffusion(phaseDiff[saltCompIdx], saltCompIdx);
                }
                const auto eIdx = gridGeometry->elementMapper().index(element);
                validationThroats << time << ',' << eIdx << ',' << inside.dofIndex() << ','
                                  << outside.dofIndex() << ',' << q << ','
                                  << cache.transmissibility(liquidPhaseIdx) << ','
                                  << diffWater << ',' << diffSalt << ','
                                  << advWater << ',' << advSalt << ','
                                  << advWater + diffWater << ',' << advSalt + diffSalt << ','
                                  << totalWaterAllPhases << ',' << totalSaltAllPhases << ','
                                  << problem->spatialParams().throatInscribedRadius(
                                         element, elemVolVars) << '\n';
            }
        }
        validationPores.flush();
        validationThroats.flush();
    };

    // get some time loop parameters
    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    const auto maxDt = getParam<Scalar>("TimeLoop.MaxTimeStepSize");
    auto dt = getParam<Scalar>("TimeLoop.DtInitial");
    auto checkPointInterval = getParam<Scalar>("TimeLoop.CheckPointInterval");

    // check if we are about to restart a previously interrupted simulation
    Scalar restartTime = 0;
    if (Parameters::getTree().hasKey("Restart") || Parameters::getTree().hasKey("TimeLoop.Restart"))
        restartTime = getParam<Scalar>("TimeLoop.Restart");

    // intialize the vtk output module
    using IOFields = GetPropType<TypeTag, Properties::IOFields>;
    PoreNetwork::VtkOutputModule<GridVariables, GetPropType<TypeTag, Properties::FluxVariables>, SolutionVector> vtkWriter(*gridVariables, x, problem->name());
    IOFields::initOutputModule(vtkWriter); //! Add model specific output fields

    vtkWriter.write(0.0);
    writeValidationDiagnostics(0.0);

    // instantiate time loop
    auto timeLoop =std::make_shared<CheckPointTimeLoop<Scalar>>(restartTime, dt, tEnd);
    timeLoop->setMaxTimeStepSize(maxDt);
//     timeLoop->setCheckPoint(checkPoints.begin(), checkPoints.end());
    timeLoop->setPeriodicCheckPoint(checkPointInterval);

    // the assembler with time loop for instationary problem
    using Assembler = FVAssembler<TypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables, timeLoop, xOld);

    // the linear solver
    using LinearSolver = Dumux::UMFPackIstlSolver<SeqLinearSolverTraits, LinearAlgebraTraitsFromAssembler<Assembler>>;
    auto linearSolver = std::make_shared<LinearSolver>();

    // the non-linear solver
#if !USETHETAREGULARIZATION
    using NewtonSolver = PoreNetwork::TwoPNewtonSolver<Assembler, LinearSolver>;
#else
    using NewtonSolver = Dumux::NewtonSolver<Assembler, LinearSolver>;
#endif
    NewtonSolver nonLinearSolver(assembler, linearSolver);

#if USEBRINE
    int throatApproach = getParam<Scalar>("Problem.ThroatApproach");
#endif

    std::ofstream outFile("evaporation.csv");
#if !ISOTHERMAL
    outFile << "#time [s], evaporation rate [mol/s], evaporation rate [kg/s], evaporation rate [m^3/s], air flux [mol/s], air flux [m^3/s], energy flux [J/s], energy conduction flux [J/s], energy convection flux [J/s], advective water flux [mol/s], diffusive water flux [mol/s], advective air flux [mol/s], diffusive air flux [mol/s], advective water flux [m^3/s], diffusive water flux [m^3/s], advective air flux [m^3/s], diffusive air flux [m^3/s], x-coordinate pb1, evaporation rate pb1 [mol/s], x-coordinate pb2, evaporation rate pb2 [mol/s], x-coordinate pb3, evaporation rate pb3 [mol/s], x-coordinate pb4, evaporation rate pb4 [mol/s], x-coordinate pb5, evaporation rate pb5 [mol/s], x-coordinate pb6, evaporation rate pb6 [mol/s] "<< std::endl;
#else
    outFile << "#time [s], evaporation rate [mol/s], evaporation rate [kg/s], evaporation rate [m^3/s], air flux [mol/s], air flux [m^3/s], advective water flux [mol/s], diffusive water flux [mol/s], advective air flux [mol/s], diffusive air flux [mol/s], advective water flux [m^3/s], diffusive water flux [m^3/s], advective air flux [m^3/s], diffusive air flux [m^3/s], x-coordinate pb1, evaporation rate pb1 [mol/s], x-coordinate pb2, evaporation rate pb2 [mol/s], x-coordinate pb3, evaporation rate pb3 [mol/s], x-coordinate pb4, evaporation rate pb4 [mol/s], x-coordinate pb5, evaporation rate pb5 [mol/s], x-coordinate pb6, evaporation rate pb6 [mol/s] "<< std::endl;
#endif

    // time loop
    timeLoop->start(); do
    {
        problem->setTime(timeLoop->time()+timeLoop->timeStepSize());
        problem->setTimeStepSize(timeLoop->timeStepSize());
        // try solving the non-linear system
        nonLinearSolver.solve(x, *timeLoop);

        // average source approach
#if USEBRINE
        if (throatApproach == 2){
            problem->calculateThroatPrecipitatedAreaAverageSource(x, *gridVariables);
        }
#endif

        // make the new solution the old solution
        xOld = x;
        gridVariables->advanceTimeStep();
        if (gridVariables->gridFluxVarsCache().invasionState().updateAfterTimeStep()){
            gridVariables->gridFluxVarsCache().invasionState().update(x, gridVariables->curGridVolVars(), gridVariables->gridFluxVarsCache());
        }

        // advance to the time loop to the next step
        timeLoop->advanceTimeStep();

        // write vtk output
//         if(problem->shouldWriteOutput(timeLoop->timeStepIndex(), *gridVariables))
//         {

//         if (timeLoop->time() > 3400 && timeLoop->time() < 3500){
//             vtkWriter.write(timeLoop->time());
//         }
//         else if(timeLoop->isCheckPoint()){
//             vtkWriter.write(timeLoop->time());
//         }

        if(timeLoop->isCheckPoint()){
            vtkWriter.write(timeLoop->time());
            writeValidationDiagnostics(timeLoop->time());
        }
        problem->postTimeStep(x, *gridVariables, outFile);
//         }

        // report statistics of this time step
        timeLoop->reportTimeStep();

        // set new dt as suggested by newton solver
        timeLoop->setTimeStepSize(nonLinearSolver.suggestTimeStepSize(timeLoop->timeStepSize()));

    } while (!timeLoop->finished());

    ////////////////////////////////////////////////////////////
    // finalize, print dumux message to say goodbye
    ////////////////////////////////////////////////////////////
    timeLoop->finalize(leafGridView.comm());

    // print dumux end message
    if (mpiHelper.rank() == 0)
    {
        Parameters::print();
        DumuxMessage::print(/*firstCall=*/false);
    }
}
