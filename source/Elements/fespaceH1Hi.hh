// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
#ifndef FILE_CFS_FESPACE_H1_HI_HH
#define FILE_CFS_FESPACE_H1_HI_HH

#include "Elements/fespaceH1.hh"

namespace CoupledField {

// forward class declaration
class FeH1Hi;

//! Finite element space for hierarchical H1 elements
class FeSpaceH1Hi : public FeSpaceH1 {

  public:

    //! Constructor
    FeSpaceH1Hi(PtrParamNode aNode, PtrParamNode infoNode);

    //! Destructor
    ~FeSpaceH1Hi();

    //! Initialize class (read order etc.)
    void Init();
        
    //! \copydoc FeSpace::GetFe(EntityIterator,shared_ptr<IntScheme>&)
    virtual BaseFE* GetFe( const EntityIterator ent ,
                           shared_ptr<IntScheme>& intScheme );

    //! \copydoc FeSpace::GetFe(EntityIterator)
    virtual BaseFE* GetFe( const EntityIterator ent );

    //! \copydoc FeSpace::GetFe(UInt)
    virtual BaseFE* GetFe( UInt elemNum );

    //! @copydoc FeSpace::GetNumEntityOrder
    UInt GetEntityOrder( UInt elemNum, BaseFE::EntityType type, 
                         UInt entityNum, UInt comp = 1 );
    
    //! @copydoc FeSpace::GetMaxEntityOrder
    UInt GetMaxEntityOrder( UInt elemNum, BaseFE::EntityType type, 
                            UInt entityNum );
    
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

    //! Adjust orders of edges / faces according to min/max rule
    void AdjustEntityOrder();
    
    //! Map for reference elements by region
    std::map< RegionIdType, std::map<Elem::FEType, FeH1Hi* > > refElems_;
    
    //! Map containing the polynomial order for every edge
    std::map<UInt, StdVector<UInt> > edgeOrder_;
    
    // For the polynomial degrees of faces / interior we have to think
    //! Map containing the polynomial order for every face
//    std::map<UInt, UInt> faceOrder_;
    
    //! Map containing the polynomial order for every element interior
//    std::map<UInt, UInt> faceOrder_;
    
    //! Test continuity of shape functions
    void TestContinuity();


  private:
};
} // end of namespace
#endif // header guard
