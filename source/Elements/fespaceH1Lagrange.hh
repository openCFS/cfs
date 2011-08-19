// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
#ifndef FILE_CFS_FESPACE_H1_LAGRANGE_HH
#define FILE_CFS_FESPACE_H1_LAGRANGE_HH

// =====================================================================================
// 
//       Filename:  fespaceH1Lagrange.hh
// 
//    Description:  This implements the finite element space for Lagrange polinomial based
//                  approximations.
//                  This space is capabale to treat the mapping of equation numbers for 
//                  elements of first and second order as well as arbitrary order tesnsor 
//                  product elements.
//                  Virtual Node Array: Simplify the implementation in the referenece 
//                                      by taking element edge and face orientations
//                                      into account
//                  
// 
//        Version:  0.1
//        Created:  01/21/2010 10:31:09 AM
//       Revision:  none
// 
//         Author:  Andreas Hueppe (AHU), andreas.hueppe@uni-klu.ac.at
//        Company:  Universitaet Klagenfurt
// 
// =====================================================================================

#include "Elements/fespaceH1.hh"
#include "Elements/H1ElemsLagExpl.hh"
#include "Elements/H1ElemsLagVar.hh"


namespace CoupledField {

class FeSpaceH1Lagrange : public FeSpaceH1 {

  public:

    //! Constructor
    FeSpaceH1Lagrange(PtrParamNode aNode );

    //! Destructor
    ~FeSpaceH1Lagrange();
    
    //! Initialize class (read order etc.)
    void Init();

    //! Return pointer to reference element
    virtual BaseFE* GetFe( const EntityIterator ent );
    
    //! Return pointer to reference element (by element number)
    virtual BaseFE* GetFe( UInt elemNum );
    
    ////! Return equation numbers
    //virtual void GetEqns( StdVector<Integer>& eqns, const EntityIterator ent ); 

    ////! Return equation numbers for a specific dof
    //virtual void GetEqns( StdVector<Integer>& eqns, const EntityIterator ent
    //                      , UInt dof );
    
    //! @see FeSpace::GetEntityOrder
    UInt GetEntityOrder( UInt elemNum, BaseFE::EntityType type, 
                           UInt entityNum, UInt comp = 1 );
    
    //! @see FeSpace:: GetMaxEntityOrder
    UInt GetMaxEntityOrder( UInt elemNum, BaseFE::EntityType type, 
                            UInt entityNum );

    //! Map equations i.e. initialize object
    virtual void Finalize();

  protected:
    // ====================================================================
    // INTERNAL INITIALIZATION
    // ====================================================================
    //! read in integration data and set defaults
    virtual void SetRegionIntegration(RegionIdType region, IntScheme::IntegMethod method, Matrix<Integer> order);

    //! Set the order and mapping type of a specific region
    virtual void SetRegionElements(RegionIdType region, MappingType mType,Matrix<Integer> order);

    //! Here the spaces have the possibility to check if user definitions makes sense
    //! e.g. if the chosen integration is correct or the element order is nice
    //! here one could e.g. adjust the integration oder according to the element order
    virtual void CheckConsistency();

    //! sets the default integration scheme and order
    virtual void SetDefaultIntegration();

    //! Create default finite elements to be used if nothing else is requested
    virtual void CreateDefaultElements();

    // ====================================================================
    // PROCESS USER INPUT
    // ====================================================================
    //! Here we pass a fePolynomial parameter node such that the feSpace can extract the information
    //! which is important for the specific space
    virtual void ProcessPolyRegionNode(PtrParamNode node, RegionIdType region);

    //! This array stores all nodes, including the virtual ones, which are available in the feSpace
    //StdVector<UInt> nodes;


  private:
    std::map<RegionIdType,bool> spectralRegions_;
};
}
#endif //
