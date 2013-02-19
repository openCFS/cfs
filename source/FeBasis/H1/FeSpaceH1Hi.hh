#ifndef FILE_CFS_FESPACE_H1_HI_HH
#define FILE_CFS_FESPACE_H1_HI_HH

#include "FeBasis/FeSpaceHi.hh"


namespace CoupledField {

// forward class declaration
class FeH1Hi;

//! Finite element space for hierarchical H1 elements
class FeSpaceH1Hi : public FeSpaceHi {

  public:

    //! Constructor
    FeSpaceH1Hi(PtrParamNode aNode, PtrParamNode infoNode, Grid* ptGrid);

    //! Destructor
    ~FeSpaceH1Hi();

    //! \copydoc FeSpace::Init()
    void Init( shared_ptr<SolStrategy> solStrat );
        
    //! \copydoc FeSpace::GetFe(EntityIterator,shared_ptr<IntScheme>&)
    virtual BaseFE* GetFe( const EntityIterator ent ,
                           IntScheme::IntegMethod& method,
                           IntegOrder & order  );

    //! \copydoc FeSpace::GetFe(EntityIterator)
    virtual BaseFE* GetFe( const EntityIterator ent );

    //! \copydoc FeSpace::GetFe(UInt)
    virtual BaseFE* GetFe( UInt elemNum );

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

    //! \copydoc FeSpaceHi::GetFeHi
    virtual FeHi* GetFeHi( RegionIdType region, Elem::FEType type );
    
    //! Map for reference elements by region
    std::map< RegionIdType, std::map<Elem::FEType, FeH1Hi* > > refElems_;
    
  private:
};
} // end of namespace
#endif // header guard
