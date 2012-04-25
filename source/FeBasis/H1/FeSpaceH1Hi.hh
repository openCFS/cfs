#ifndef FILE_CFS_FESPACE_H1_HI_HH
#define FILE_CFS_FESPACE_H1_HI_HH

#include<boost/unordered_set.hpp>
#include<boost/unordered_map.hpp>


#include "FeSpaceH1.hh"



namespace CoupledField {

// forward class declaration
class FeH1Hi;

//! Finite element space for hierarchical H1 elements
class FeSpaceH1Hi : public FeSpaceH1 {

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

    //! Map for reference elements by region
    std::map< RegionIdType, std::map<Elem::FEType, FeH1Hi* > > refElems_;
    
    //! Map for each region the element order
    std::map< RegionIdType, ApproxOrder > regionOrder_;
    
    // ====================================================================
    // VARIABLE ENTITY ORDER
    // ====================================================================
    
    //! Initialize a finite element with the correct order
    void SetElemOrder( const Elem* ptEl, FeH1Hi* ptFe, 
                       const ApproxOrder& order ,
                       bool applyMinMaxRule );
    
    //! Adjust order of edges / faces which have mixed order neighbours
    void AdjustEntityOrder();
    
    //! Set polynomial order of vectorial unknowns for for anisotropic order
    
    //! In case of vectorial unknowns and anisotropic polynomial order
    //! (= different order for every component of the vector), the
    //! higher order nodes might have to be fixed by a homogeneous Dirichlet
    //! BC.
    void FixHigherOrderAnisoDofs();

    //! Set containing all edges, where the min/max rule gets applied
    boost::unordered_set<UInt> adjustedEdges_;
    
    //! Set containing all faces, where the min/max rule gets applied
    boost::unordered_set<UInt> adjustedFaces_;
    
    //! Map containing the order of adjusted edges (key: edge number)
    boost::unordered_map<UInt, UInt> orderEdges_;
    
    //! Map containing the order of adjusted faces (key: face number)
    boost::unordered_map<UInt, boost::array<UInt,2> > orderFaces_;

  private:
};
} // end of namespace
#endif // header guard
