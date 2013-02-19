// =====================================================================================
// 
//       Filename:  FeSpaceH1Nodal.hh
// 
//    Description:  This implements the finite element space for nodal finite elements
//                  In our understanding this includes every element which features local
//                  or "nodal" support. i.e. the basis functions become one only at a single dof!
//
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

#ifndef FILE_CFS_FESPACE_H1_NODAL_HH
#define FILE_CFS_FESPACE_H1_NODAL_HH

#include "FeBasis/FeSpace.hh"
#include "H1ElemsLagExpl.hh"
#include "H1ElemsLagVar.hh"


namespace CoupledField {

class FeSpaceH1Nodal : public FeSpace{

  public:

    //! Constructor
    FeSpaceH1Nodal(PtrParamNode aNode, PtrParamNode infoNode, Grid* ptGrid );

    //! Destructor
    virtual ~FeSpaceH1Nodal();
    
    //! \copydoc FeSpace::Init()
    void Init( shared_ptr<SolStrategy> solStrat );

    //! \copydoc FeSpace::GetFe(EntityIterator,shared_ptr<IntScheme>&)
    virtual BaseFE* GetFe( const EntityIterator ent ,
                           IntScheme::IntegMethod& method,
                           IntegOrder & order );
    
    //! \copydoc FeSpace::GetFe(EntityIterator)
    virtual BaseFE* GetFe( const EntityIterator ent );
    
    //! \copydoc FeSpace::GetFe(UInt)
    virtual BaseFE* GetFe( UInt elemNum );
    
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
                                    const ApproxOrder& order,
                                    PtrParamNode infoNode );

    //! Here the spaces have the possibility to check if user definitions makes sense
    //! e.g. if the chosen integration is correct or the element order is nice
    //! here one could e.g. adjust the integration oder according to the element order
    virtual void CheckConsistency();

    //! sets the default integration scheme and order
    virtual void SetDefaultIntegration( PtrParamNode infoNode );

    //! Create default finite elements to be used if nothing else is requested
    virtual void SetDefaultElements( PtrParamNode infoNode );

    //! \copydoc FeSpace::MapNodalBCs
    virtual void MapNodalBCs();
    
    // ====================================================================
    // PROCESS USER INPUT
    // ====================================================================
    //! reads in special options for the space under consideration
    //! here it cosideres the spectral flag
    virtual void ReadCustomAttributes(PtrParamNode pNode,RegionIdType region);

    //! Set with all regions being treated as spectral
    std::set<RegionIdType> spectralRegions_;
    
    //! Map for reference elements by region
    std::map< RegionIdType, std::map<Elem::FEType, FeH1* > > refElems_;
    
    //! Order of the space
 
    //! As the Lagrange shape functions have to be continuous, we only 
    //! allow one specific order of the shape functions. 
    UInt order_;
};
}
#endif //
