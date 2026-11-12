// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PoreNetworkModels
 * \ingroup SpatialParameters
 * \brief Spatial parameters for a 2pnc pore-network model
 */

#ifndef DUMUX_PNM_2P_NC_MIN_SPATIAL_PARAMS_HH
#define DUMUX_PNM_2P_NC_MIN_SPATIAL_PARAMS_HH

#include <dumux/porenetwork/2p/spatialparams.hh>

namespace Dumux::PoreNetwork {

template<class GridGeometry, class Scalar, class MaterialLawT>
class TwoPNCMinSpatialParams
: public TwoPSpatialParams<GridGeometry, Scalar, MaterialLawT, TwoPNCMinSpatialParams<GridGeometry, Scalar, MaterialLawT>>
{
    using ParentType = TwoPSpatialParams<GridGeometry, Scalar, MaterialLawT, TwoPNCMinSpatialParams<GridGeometry, Scalar, MaterialLawT>>;

    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    //some Indexes for convinience
    static const int sNaClIdx  = 0;
    static const int sPhaseIdx = 3;

    enum
    {
         // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,
    };
public:
    using PermeabilityType = Scalar;
    using ParentType::ParentType;

    TwoPNCMinSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        temperature_ = getParam<Scalar>("SpatialParams.Temperature", 283.15);
        gamma_ = getParam<Scalar>("SpatialParams.SurfaceTension", 0.0725); // default to surface tension of water/air
        theta_ = getParam<Scalar>("SpatialParams.ContactAngle", 0.0);

#if USEBRINE
        poreVolumeChange_ = getParam<bool>("Problem.PoreVolumeChange", true);
        reducedEvaporativeArea_ = getParam<bool>("Problem.ReducedEvaporativeArea", false);
        throatRadiusChange_ = getParam<bool>("Problem.ThroatRadiusChange", true);
        throatApproach_ = getParam<Scalar>("Problem.ThroatApproach", 2);
#else
        poreVolumeChange_ = getParam<bool>("Problem.PoreVolumeChange", false);
        reducedEvaporativeArea_ = getParam<bool>("Problem.ReducedEvaporativeArea", false);
        throatRadiusChange_ = getParam<bool>("Problem.ThroatRadiusChange", false);
        throatApproach_ = getParam<Scalar>("Problem.ThroatApproach", 2);
#endif

        if(this->gridGeometry().throatCrossSectionShape()[0] == Throat::Shape::rectangle)
        {
            rectangleApproach_ = getParam<Scalar>("Grid.RectangleApproach");
            initialRectangleWidth_ = 2 * initialThroatRadius(0); //TODO just valid if all cross section are the same initially
            initialRectangleHeight_ = getParam<Scalar>("Grid.ThroatHeight");
            rectangleRatio_ = initialRectangleHeight_ / initialRectangleWidth_;
        }

        if (throatApproach_ == 2){
            int numElements = gridGeometry->gridView().size(0);

            throatPrecipitatedAreaInTimeStep_.resize(numElements);
            for (int i = 0; i < throatPrecipitatedAreaInTimeStep_.size(); ++i)
            throatPrecipitatedAreaInTimeStep_[i] = 0.0;

            throatPrecipitatedAreaTotal_.resize(numElements);
            for (int i = 0; i < throatPrecipitatedAreaTotal_.size(); ++i)
            throatPrecipitatedAreaTotal_[i] = 0.0;
        }
    }

    //! \brief Function for defining the temperature
    Scalar temperatureAtPos(const GlobalPosition& globalPos) const
    { return temperature_; }

    //! \brief Function for defining the wetting phase
    template<class FluidSystem>
    int wettingPhaseAtPos(const GlobalPosition& globalPos) const
    { return FluidSystem::phase0Idx; }


    //! \brief Function for defining the contact angle
    int contactAngleAtPos(const GlobalPosition& globalPos) const
    { return theta_; }

    //! \brief Function for defining the surface tension
    Scalar surfaceTensionAtPos(const GlobalPosition& globalPos) const
    { return gamma_; }

    template<class SolidSystem, class ElementSolution>
    Scalar inertVolumeFraction(const Element& element,
                               const SubControlVolume& scv,
                               const ElementSolution& elemSol,
                               int compIdx) const
    {
        // there are no inert components at all
        return 0.0;
    }

    template<class ElementSolution>
    Scalar porosity(const Element& element,
                    const SubControlVolume& scv,
                    const ElementSolution& elemSol) const
    {
            const Scalar precipVolFrac = elemSol[scv.localDofIndex()][sPhaseIdx];
            return porosity(element, scv, precipVolFrac);
    }


    Scalar porosity(const Element& element,
                    const SubControlVolume& scv,
                    const Scalar precipVolFrac) const
    {
        if(!poreVolumeChange_){
            return 1.0;
        }
        else {
            const Scalar initPoreVolume = initialPoreVolume(scv.dofIndex());
            const Scalar poreVol = poreVolume(initPoreVolume, precipVolFrac);
            return poreVol/initPoreVolume;
        }
    }

    //! Required for compatibility reasons with porous medium-flow models.
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    {   if(!poreVolumeChange_){
            return 1.0;
        }
        else {DUNE_THROW(Dune::InvalidStateException, "Porosity change not implemented for porosityAtPos function.");}
    }

///////////////////////////////////////////////////////////////////////////
// Pore body
///////////////////////////////////////////////////////////////////////////
    /*!
     * \brief Returns the actual pore volume.
     */
    template<class VolumeVariables>
    Scalar poreVolume(const SubControlVolume& scv,
                      const VolumeVariables& volVars) const
    {
        return poreVolume(scv.dofIndex(), volVars);
    }

    /*!
     * \brief Returns the actual pore volume.
     */
    template<class VolumeVariables>
    Scalar poreVolume(const int dofIdx,
                      const VolumeVariables& volVars) const
    {
        const Scalar precipVolFrac = volVars.solidVolumeFraction(sNaClIdx);
        const Scalar initPoreVolume = initialPoreVolume(dofIdx);
        return poreVolume(initPoreVolume, precipVolFrac);
    }

    /*!
     * \brief Returns the actual pore volume.
     *
     * \param initialPoreVolume
     * \param precipitateVolumeFraction
     */
    Scalar poreVolume(const Scalar initialPoreVolume, const Scalar precipitateVolumeFraction) const
    {
        if(!poreVolumeChange_){
            return initialPoreVolume;
        }
        else {
            return (1.0 - precipitateVolumeFraction) * initialPoreVolume;
        }
    }

    /*!
     * \brief Returns the initial pore volume without any volume fraction of the solid phase.
     */
    Scalar initialPoreVolume(const int dofIdx) const
    { return this->gridGeometry().poreVolume(dofIdx); }

    /*!
     * \brief Returns the actual pore radius.
     */
    template<class ElementSolutionVector>
    Scalar poreInscribedRadius(const Element& element,
                               const SubControlVolume& scv,
                               const ElementSolutionVector& elemSol) const
    {
        const Scalar initPoreVolume = initialPoreVolume(scv.dofIndex());
        const Scalar precipVolFrac = elemSol[scv.localDofIndex()][sPhaseIdx];
        const Scalar factor = volumeFactor_(scv);
        return poreInscribedRadius(initPoreVolume, precipVolFrac, factor);
    }


    /*!
     * \brief Returns the actual pore radius.
     */
    template<class VolumeVariables>
    Scalar poreInscribedRadius(const SubControlVolume& scv,
                               const VolumeVariables volVars) const
    {
        const Scalar precipVolFrac = volVars.solidVolumeFraction(sNaClIdx);
        const Scalar initPoreVolume = initialPoreVolume(scv.dofIndex());
        const Scalar factor = volumeFactor_(scv);
        return poreInscribedRadius(initPoreVolume, precipVolFrac, factor);
    }


    Scalar poreInscribedRadius(const Scalar initialPoreVolume, const Scalar precipitateVolumeFraction, const Scalar factor) const
    {
        const Scalar poreVol = poreVolume(initialPoreVolume, precipitateVolumeFraction);
        Scalar volume = poreVol * factor;

        return poreInscribedRadius_(volume);
    }


    /*!
     * \brief Returns the initial pore radius without any volume fraction of the solid phase.
     */
    Scalar initialPoreRadius(const int dofIndex) const
    { return this->gridGeometry().poreInscribedRadius(dofIndex); }


    template<class VolumeVariables>
    Scalar poreCrossSectionalArea(SubControlVolume scv, const VolumeVariables volVars) const
	{
        Scalar r = 0.0;
        if(reducedEvaporativeArea_){
            r = poreInscribedRadius(scv, volVars);
//             std::cout << "reducedEvaporativeArea_: inscribed radius" << r << std::endl;
        }else{
            r = initialPoreRadius(scv.dofIndex());
//             std::cout << "reducedEvaporativeArea_: initial radius" << r << std::endl;
        }

        // TODO find formulation for cylindrical pore-body and aspectRatio function
        const auto& poreGeometry = this->gridGeometry().poreGeometry();

        if (this->gridGeometry().useSameGeometryForAllPores())
        {
            switch (poreGeometry[0])
            {
                case PoreNetwork::Pore::Shape::square: return 4 * r * r; break; //could also be a 1-dimensional cross-section with A = 2*r
                case PoreNetwork::Pore::Shape::circle: return M_PI * r * r; break; //could also be a 1-dimensional cross-section with A = 2*r
                case PoreNetwork::Pore::Shape::cube: return 4 * r * r; break;
                case PoreNetwork::Pore::Shape::sphere: return M_PI * r * r; break;
                case PoreNetwork::Pore::Shape::cylinder: DUNE_THROW(Dune::InvalidStateException, "Cylindrical pore-body geometry not implemented for mineralization in pore-network model.");
                default: DUNE_THROW(Dune::InvalidStateException, "Unknown shape or shape not implemented for mineralization pore-network model!");
            }
        }
        else
            DUNE_THROW(Dune::InvalidStateException, "Varying pore-body geometries are not implemented for the mineralization model.");
    }

    template<class VolumeVariables>
    Scalar poreReactiveSurfaceArea(SubControlVolume scv, const VolumeVariables volVars) const
	{
        const Scalar r = poreInscribedRadius(scv, volVars);

        // TODO find formulation for cylindrical pore-body and aspectRatio function
        const auto& poreGeometry = this->gridGeometry().poreGeometry();

        if (this->gridGeometry().useSameGeometryForAllPores())
        {
            switch (poreGeometry[0])
            {
                case PoreNetwork::Pore::Shape::cube: return 6 * 2 * r * 2 * r; break;
                default: DUNE_THROW(Dune::InvalidStateException, "Unknown shape or shape not implemented for mineralization pore-network model!");
            }
        }
        else
            DUNE_THROW(Dune::InvalidStateException, "Varying pore-body geometries are not implemented for the mineralization model.");
    }

///////////////////////////////////////////////////////////////////////////
// Pore throat
///////////////////////////////////////////////////////////////////////////

    /*!
     * \brief Returns the actual pore-throat radius.
     */
    template<class ElementVolumeVariables>
    Scalar throatInscribedRadius(const Element& element, const ElementVolumeVariables& elemVolVars) const
    {
        const auto eIdx = this->gridGeometry().elementMapper().index(element);
        const auto& throatCrossSectionShape = this->gridGeometry().throatCrossSectionShape(eIdx);

        if(!throatRadiusChange_){
            return initialThroatRadius(eIdx);
        }else{
            //average source terms approach
            if (throatApproach_ == 2){
                    const auto Aini = initialThroatCrossSectionalArea(element, elemVolVars);
                    const auto A = Aini-throatPrecipitatedAreaTotal_[eIdx];

                    switch (throatCrossSectionShape){
                        case Throat::Shape::square: return 0.5 * pow(A, 0.5); break;
                        case Throat::Shape::circle: return pow(A/M_PI, 0.5); break;
                        case Throat::Shape::rectangle:
                            if(rectangleApproach_ == 0){
                                return pow(A / (4 * rectangleRatio_), 0.5);
                            }else if (rectangleApproach_ == 1){
                                Scalar a = 4.0;
                                Scalar b = 2.0 * (initialRectangleWidth_-initialRectangleHeight_);
                                Scalar c = A;
                                Scalar radius = (b+sqrt(b*b + 4*a*c))/(2*a);
                                return radius;
                            }else{
                                DUNE_THROW(Dune::InvalidStateException, "Rectangle approach not available.");
                            }
                            break;
                        default: DUNE_THROW(Dune::InvalidStateException, "Pore cross section shape is not implemented for average source terms approach (approach 2).");
                    }
            }
            else if(throatApproach_ == 1){
                auto fvGeometry = localView(this->gridGeometry());
                fvGeometry.bindElement(element);
                if (throatCrossSectionShape==Throat::Shape::circle){
                    std::vector<Scalar> r = throatEndRadiiEvenGrowth(element, elemVolVars);
                    Scalar averageThroatRadius = pow((3*pow(r[0], 3.0)*pow(r[1], 3.0))/(pow(r[0], 2.0) + r[0]*r[1] + pow(r[1], 2.0)), 0.25);
                    return averageThroatRadius;
                 }else if(throatCrossSectionShape==Throat::Shape::square){
                    Scalar wi = 2.0*throatEndRadiiEvenGrowth(element, elemVolVars)[0];
                    Scalar wj = 2.0*throatEndRadiiEvenGrowth(element, elemVolVars)[1];
                    Scalar averageThroatWidth = pow((3*pow(wi, 3.0)*pow(wj, 3.0))/(pow(wi, 2.0) + wi*wj + pow(wj, 2.0)), 0.25);
                    return averageThroatWidth/2.0;
                 }

//                     else if(throatCrossSectionShape[0]==Throat::Shape::rectangle){
//                     std::vector<Scalar> r = throatEndRadiiEvenGrowth(element, elemVolVars);
//                     if(rectangleApproach_ == 0 || rectangleApproach_ == 2){
//                         std::vector<Scalar> heights(fvGeometry.numScv(), 0);
//                         std::vector<Scalar> widths(fvGeometry.numScv(), 0);
//                         for (int i = 0; i < r.size(); ++i){
//                             widths[i]  = 2.0 * r[i];
//                             heights[i]  = rectangleRatio_ * widths[i];
//                         }
//                         Scalar averageHeight = pow((3*pow(heights[0], 3.0)*pow(heights[1], 3.0))/(pow(heights[0], 2.0) + heights[0]*heights[1] + pow(heights[1], 2.0)), 0.25);
//                         Scalar averageWidth = averageHeight / rectangleRatio_;
//                         return 0.5 * averageWidth;
//                     }else if (rectangleApproach_ == 1){
//                         std::vector<Scalar> heights(fvGeometry.numScv(), 0);
//                         std::vector<Scalar> widths(fvGeometry.numScv(), 0);
//                         for (int i = 0; i < r.size(); ++i){
//                             widths[i]  = 2.0 * r[i];
//                             heights[i]  = initialRectangleHeight_ - initialRectangleWidth_ + widths[i];
//                         }
//                         Scalar averageHeight = pow((3*pow(heights[0], 3.0)*pow(heights[1], 3.0))/(pow(heights[0], 2.0) + heights[0]*heights[1] + pow(heights[1], 2.0)), 0.25);
//
//                         Scalar area = averageHeight * (initialRectangleWidth_-initialRectangleHeight_+averageHeight);
//                         return 0.5 * area / averageHeight;
//                     }else{
//                         DUNE_THROW(Dune::InvalidStateException, "Rectangle approach not available.");
//                     }
//                 }
                else
                    DUNE_THROW(Dune::InvalidStateException, "Pore cross section shape is not implemented for equal growth approach (approach 1).");
            }
            else{
                DUNE_THROW(Dune::InvalidStateException, "Throat approch not implemented.");
            }
        }
    }

    /*!
     * \brief Returns the initial pore-throat radius without any volume fraction of the solid phase.
     */
    Scalar initialThroatRadius(const int eIdx) const
	{ return this->gridGeometry().throatInscribedRadius(eIdx); }

	/*!
     * \brief Returns the actual pore-throat radii at each throat end for the even growth approach (approach 1).
     * Only used for approach 1.
     */
    template<class ElementVolumeVariables>
    std::vector<Scalar> throatEndRadiiEvenGrowth(const Element& element, const ElementVolumeVariables& elemVolVars) const
    {
        if (throatApproach_ == 1){
            const auto eIdx = ParentType::gridView().indexSet().index(element);
            Scalar throatInitialRadius = initialThroatRadius(eIdx);
            const auto& throatCrossSectionShape = this->gridGeometry().throatCrossSectionShape(eIdx);

            auto fvGeometry = localView(this->gridGeometry());
            fvGeometry.bindElement(element);
            std::vector<Scalar> throatRadii(fvGeometry.numScv(), 0);

//             if (this->gridGeometry().throatCrossSectionShape()[0] == Throat::Shape::rectangle && rectangleApproach_ == 2){
//                 std::vector<Scalar> throatWidths(fvGeometry.numScv(), 0);
//                 for (const auto& scv : scvs(fvGeometry)){
//                     Scalar poreRadius = poreInscribedRadius(scv, elemVolVars[scv.indexInElement()]);
//                     Scalar poreInitialRadius = initialPoreRadius(scv.dofIndex());
//                     Scalar porePrecipitateThickness = poreInitialRadius - poreRadius;
//                     throatWidths[scv.indexInElement()] = pow(2*throatInitialRadius - 2*porePrecipitateThickness / rectangleRatio_, 0.5) * pow(2*throatInitialRadius - 2*porePrecipitateThickness, 0.5);
//                     throatRadii[scv.indexInElement()] = 0.5 * throatWidths[scv.indexInElement()];
//                 }
//             }
            if(throatCrossSectionShape==Throat::Shape::square || throatCrossSectionShape==Throat::Shape::circle || throatCrossSectionShape==Throat::Shape::equilateralTriangle){
                for (const auto& scv : scvs(fvGeometry)){
                    Scalar poreRadius = poreInscribedRadius(scv, elemVolVars[scv.indexInElement()]);
                    Scalar poreInitialRadius = initialPoreRadius(scv.dofIndex());
                    Scalar porePrecipitateThickness = poreInitialRadius - poreRadius;
                    throatRadii[scv.indexInElement()] = throatInitialRadius - porePrecipitateThickness;
                }
            }
            return throatRadii;
        }
        else{
            DUNE_THROW(Dune::InvalidStateException, "This function is just necessary for the even growth approach (approach 1) and should not be used by the other approaches.");
        }
    }


    /*!
     * \brief Returns the actual pore-throat cross-sectional area.
     */
   template<class ElementVolumeVariables>
   Scalar throatCrossSectionalArea(const Element& element,
                                   const ElementVolumeVariables& elemVolVars) const
   {
        const Scalar r = throatInscribedRadius(element, elemVolVars);
        const auto& throatCrossSectionShape = this->gridGeometry().throatCrossSectionShape();

        if (this->gridGeometry().useSameShapeForAllThroats())
        {
            switch (throatCrossSectionShape[0])
            {
                case PoreNetwork::Throat::Shape::equilateralTriangle: return pow(3, 0.5) * r * r; break;
                case PoreNetwork::Throat::Shape::square: return 4 * r * r; break;
                case PoreNetwork::Throat::Shape::circle: return M_PI * r * r ; break;
                case PoreNetwork::Throat::Shape::twoPlates: return 2 * r; break;
                case PoreNetwork::Throat::Shape::rectangle: DUNE_THROW(Dune::InvalidStateException, "Rectangle throat cross-section geometry not implemented for mineralization in pore-network model.");
                case PoreNetwork::Throat::Shape::scaleneTriangle: DUNE_THROW(Dune::InvalidStateException, "Scalene triangle throat cross-section geometry not implemented for mineralization in pore-network model.");
                case PoreNetwork::Throat::Shape::polygon: DUNE_THROW(Dune::InvalidStateException, "Polygon throat cross-section geometry not implemented for mineralization in pore-network model.");
                default: DUNE_THROW(Dune::InvalidStateException, "Unknown shape!");
            }
        }
        else
            DUNE_THROW(Dune::InvalidStateException, "Varying throat geometries are not implemented for the mineralization model.");
   }

    /*!
     * \brief Returns the initial pore-throat cross-sectional area without any volume fraction of the solid phase.
     */
   template<class ElementVolumeVariables>
   Scalar initialThroatCrossSectionalArea(const Element& element,
                                    const ElementVolumeVariables& elemVolVars) const
   {
       const auto eIdx = this->gridGeometry().elementMapper().index(element);
       return this->gridGeometry().throatCrossSectionalArea(eIdx);
   }


    template<class ElementVolumeVariables>
    Scalar throatReactiveSurfaceArea(const Element& element,
                                    const ElementVolumeVariables& elemVolVars) const
	{
        const Scalar r = throatInscribedRadius(element, elemVolVars);
        const auto& throatCrossSectionShape = this->gridGeometry().throatCrossSectionShape();
        const Scalar throatL = throatLength(element, elemVolVars);

        if (this->gridGeometry().useSameShapeForAllThroats())
        {
            switch (throatCrossSectionShape[0])
            {
                case PoreNetwork::Throat::Shape::square: return 8 * r * throatL; break;
                case PoreNetwork::Throat::Shape::circle: return 2 * M_PI * r * throatL; break;
                case PoreNetwork::Throat::Shape::equilateralTriangle: DUNE_THROW(Dune::InvalidStateException, "equilateral triangle throat cross-section geometry not implemented for mineralization in pore-network model.");
                case PoreNetwork::Throat::Shape::twoPlates: DUNE_THROW(Dune::InvalidStateException, "Two plates throat cross-section geometry not implemented for mineralization in pore-network model.");
                case PoreNetwork::Throat::Shape::rectangle: DUNE_THROW(Dune::InvalidStateException, "Rectangle throat cross-section geometry not implemented for mineralization in pore-network model.");
                case PoreNetwork::Throat::Shape::scaleneTriangle: DUNE_THROW(Dune::InvalidStateException, "Scalene triangle throat cross-section geometry not implemented for mineralization in pore-network model.");
                case PoreNetwork::Throat::Shape::polygon: DUNE_THROW(Dune::InvalidStateException, "Polygon throat cross-section geometry not implemented for mineralization in pore-network model.");
                default: DUNE_THROW(Dune::InvalidStateException, "Unknown shape!");
            }
        }
        else
            DUNE_THROW(Dune::InvalidStateException, "Varying throat geometries are not implemented for the mineralization model.");
    }

    /*!
     * \brief Returns the actual pore-throat length.
     */
    template<class ElementVolumeVariables>
    Scalar throatLength(const Element& element,
                        const ElementVolumeVariables& elemVolVars) const
    {
//         auto fvGeometry = localView(this->gridGeometry());
//         fvGeometry.bindElement(element);
//
//         //add reduced pore-body radii to pore throat length
//         Scalar additionalPoreLength = 0;
//         for(int i = 0; i < fvGeometry.numScv() ; ++i)
//         {
//             const auto dofIdxGlobal = ParentType::gridView().indexSet().subIndex(element, i, dim);
//             additionalPoreLength += initialPoreRadius(dofIdxGlobal) - poreInscribedRadius(dofIdxGlobal, elemVolVars[i]);
//         }
        return initialThroatLength(element);/* + additionalPoreLength;*/

    }

    /*!
     * \brief Returns the initial pore-throat length without any volume fraction of the solid phase.
     */
    Scalar initialThroatLength(const Element& element) const
    {
        const auto eIdx = this->gridGeometry().elementMapper().index(element);
        return this->gridGeometry().throatLength(eIdx);
    }

    /*!
     * \brief Used to set the precipitated area for the averaged source approach calculated in the problem where the source term is available (approach 2).
     * Only used for approach 2.
     */
    void setThroatPrecipitatedAreaAverageSource(const std::vector<Scalar>& timeStepPrecipitatedArea){
        if(throatApproach_ == 2){
            for (const auto& element : elements(this->gridGeometry().gridView())){
                const int eIdx = ParentType::gridView().indexSet().index(element);
                throatPrecipitatedAreaInTimeStep_[eIdx] = timeStepPrecipitatedArea[eIdx];
                throatPrecipitatedAreaTotal_[eIdx] += timeStepPrecipitatedArea[eIdx];
            }
        }
        else{
            DUNE_THROW(Dune::InvalidStateException, "Function just used for the averaged source approach (approach 2).");
        }
    }


private:
    bool onBoundary_(const GlobalPosition& pos, std::size_t boundaryIdx) const
    {
        constexpr auto eps = 1e-8; //TODO
        //! the bounding box of the whole domain
        GlobalPosition bBoxMin_ = this->gridGeometry().bBoxMin();
        GlobalPosition bBoxMax_ = this->gridGeometry().bBoxMax();

        static constexpr std::array boundaryIdxToCoordinateIdx{0, 0, 1, 1, 2, 2};
        const auto coordinateIdx = boundaryIdxToCoordinateIdx[boundaryIdx];
        auto isMaxBoundary = [] (int n) { return (n % 2 != 0); };

        if (isMaxBoundary(boundaryIdx))
            return pos[coordinateIdx] > bBoxMax_[coordinateIdx] - eps;
        else
            return pos[coordinateIdx] < bBoxMin_[coordinateIdx] + eps;
    }


    /*!
     * \brief Returns the actual pore radius.
     */
    Scalar poreInscribedRadius_(const Scalar volume) const
	{
        // TODO find formulation for cylindrical pore-body and aspectRatio function
        const auto& poreGeometry = this->gridGeometry().poreGeometry();

        if (this->gridGeometry().useSameGeometryForAllPores())
        {
            switch (poreGeometry[0])
            {
                case PoreNetwork::Pore::Shape::square: return 0.5 * pow(volume, 0.5); break;
                case PoreNetwork::Pore::Shape::circle: return pow(volume/M_PI, 0.5); break;
                case PoreNetwork::Pore::Shape::cube: return 0.5 * pow(volume, (1.0/3.0)); break;
                case PoreNetwork::Pore::Shape::sphere: return pow((3*volume)/(4*M_PI), (1.0/3.0)); break;
                case PoreNetwork::Pore::Shape::cylinder: DUNE_THROW(Dune::InvalidStateException, "Cylindrical pore-body geometry not implemented for mineralization in pore-network model.");
                default: DUNE_THROW(Dune::InvalidStateException, "Unknown shape or shape not implemented for mineralization pore-network model!");
            }
        }
        else
            DUNE_THROW(Dune::InvalidStateException, "Varying pore-body geometries are not implemented for the mineralization model.");
    }

    Scalar volumeFactor_(const SubControlVolume& scv) const
    {
        const auto capPoresOnBoundaries = getParam<std::vector<int>>("Grid.CapPoresOnBoundaries", std::vector<int>{});

        Scalar factor = 1.0;
        if (capPoresOnBoundaries.empty())
                factor = 1.0;
        else
        {
            std::size_t numCaps = 0;
//             const auto& pos = scv.center();
            const auto& pos = scv.dofPosition();
            for (auto boundaryIdx : capPoresOnBoundaries)
            {
                if (onBoundary_(pos, boundaryIdx))
                {
                    factor = 2.0;
                    ++numCaps;
                }
            }
            if (numCaps > dimWorld)
                DUNE_THROW(Dune::InvalidStateException, "Pore at " << pos << " capped " << numCaps << " times. Capping should not happen more than " << dimWorld << " times");
        }
        return factor;
    }


    Scalar temperature_;
    Scalar gamma_;
    Scalar theta_;
    bool poreVolumeChange_;
    bool throatRadiusChange_;
    Scalar throatApproach_;
    bool kineticPrecipitationReaction_;
    bool reducedEvaporativeArea_;

    std::vector<Scalar> throatPrecipitatedAreaInTimeStep_;
    std::vector<Scalar> throatPrecipitatedAreaTotal_;

//     if(this->gridGeometry().throatCrossSectionShape()[0] == Throat::Shape::rectangle)
//     {
        Scalar rectangleApproach_;
        Scalar initialRectangleWidth_; //TODO just valid if all cross section are the same initially
        Scalar initialRectangleHeight_;
        Scalar rectangleRatio_;
//     }

};
} // end namespace Dumux::PoreNetwork

#endif
