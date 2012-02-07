// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     FeSpaceL2Nodal.hh
 *       \brief    <Description>
 *
 *       \date     Feb 2, 2012
 *       \author   ahueppe
 */
//================================================================================================

#ifndef FESPACEL2NODAL_HH_
#define FESPACEL2NODAL_HH_

#include "FeSpaceL2.hh"
#include "FeBasis/H1/H1ElemsLagExpl.hh"
#include "FeBasis/H1/H1ElemsLagVar.hh"

namespace CoupledField {

class FeSpaceL2Nodal : public FeSpaceL2 {
public:
     //! Constructor
     FeSpaceL2Nodal(PtrParamNode aNode, PtrParamNode infoNode );

     //! Destructor
     virtual ~FeSpaceL2Nodal();

     //! Initialize class (read order etc.)
     void Init();
     //! \copydoc FeSpace::GetFe(EntityIterator,shared_ptr<IntScheme>&)
    virtual BaseFE* GetFe( const EntityIterator ent ,
                           shared_ptr<IntScheme>& intScheme );

    //! \copydoc FeSpace::GetFe(EntityIterator)
    virtual BaseFE* GetFe( const EntityIterator ent );

    //! \copydoc FeSpace::GetFe(UInt)
    virtual BaseFE* GetFe( UInt elemNum );

    ////! Return equation numbers
    //virtual void GetEqns( StdVector<Integer>& eqns, const EntityIterator ent );

    ////! Return equation numbers for a specific dof
    //virtual void GetEqns( StdVector<Integer>& eqns, const EntityIterator ent
    //                      , UInt dof );

    //! \copydoc FeSpace::GetEntityOrder
    UInt GetEntityOrder( UInt elemNum, BaseFE::EntityType type,
                           UInt entityNum, UInt comp = 1 );

    //! \copydoc FeSpace:: GetMaxEntityOrder
    UInt GetMaxEntityOrder( UInt elemNum, BaseFE::EntityType type,
                            UInt entityNum );

    //! Precalculate integration points
    virtual void PreCalcShapeFncs();

    //! Map equations i.e. initialize object
    virtual void Finalize();

  protected:
    // ====================================================================
    // INTERNAL INITIALIZATION
    // ====================================================================

    //! Set the order and mapping type of a specific region
    virtual void SetRegionElements( RegionIdType region,
                                    MappingType mType,
                                    const Matrix<Integer>& order,
                                    PtrParamNode infoNode );

    //! Here the spaces have the possibility to check if user definitions makes sense
    //! e.g. if the chosen integration is correct or the element order is nice
    //! here one could e.g. adjust the integration oder according to the element order
    virtual void CheckConsistency();

    //! sets the default integration scheme and order
    virtual void SetDefaultIntegration( PtrParamNode infoNode );

    //! Create default finite elements to be used if nothing else is requested
    virtual void SetDefaultElements( PtrParamNode infoNode );

    // ====================================================================
    // PROCESS USER INPUT
    // ====================================================================
    //! reads in special options for the space under consideration
    //! here it cosideres the spectral flag
    virtual void ReadCustomAttributes(PtrParamNode pNode,RegionIdType region);

    //! This array stores all nodes, including the virtual ones, which are available in the feSpace
    //StdVector<UInt> nodes;

    //! Set with all regions being treated as spectral
    std::set<RegionIdType> spectralRegions_;

    //! Map for reference elements by region
    std::map< RegionIdType, std::map<Elem::FEType, FeH1* > > refElems_;
};

}

#endif /* FESPACEL2NODAL_HH_ */
