#ifndef FILE_CFS_FESPACE_HCURL_HI_HH
#define FILE_CFS_FESPACE_HCURL_HI_HH

#include "FeBasis/H1/FeSpaceH1.hh"
#include <boost/array.hpp>

namespace CoupledField {

// forward class declaration
class FeHCurlHi;

//! Finite element space for hierarchical H1 elements
class FeSpaceHCurlHi : public FeSpaceH1 {

public:

  //! Constructor
  FeSpaceHCurlHi(PtrParamNode aNode, PtrParamNode infoNode, Grid* ptGrid);

  //! Destructor
  ~FeSpaceHCurlHi();
  
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
  
  //! \copydoc FeSpace::UpdateToSolStrategy()
  virtual void UpdateToSolStrategy();
  
  //! Return SBM-block and Matrix-SubBlockdefinition according to strategy
  virtual void GetOlasMappings( StdVector<AlgebraicSys::SBMBlockDef>& sbmBlocks,
                                std::map<UInt,StdVector<std::set<Integer> > >&
                                minorBlocks );


protected:

  // ====================================================================
  // INTERNAL INITIALIZATION
  // ====================================================================

  //! Set the order and mapping type of a specific region
  virtual void SetRegionElements( RegionIdType region, 
                                  MappingType mType,
                                  const Matrix<Integer>& order,
                                  PtrParamNode infoNode );

  //! \copydoc FeSpace::GetNodesOfElement
  virtual void GetNodesOfElement( StdVector<UInt>& nodes,
                                  const Elem* ptElem,
                                  BaseFE::EntityType entType = BaseFE::ALL);
  
  //! Here the spaces have the possibility to check if user definitions makes sense
  //! e.g. if the chosen integration is correct or the element order is nice
  //! here one could e.g. adjust the integration oder according to the element order
  virtual void CheckConsistency();

  //! sets the default integration scheme and order
  virtual void SetDefaultIntegration( PtrParamNode infoNode );

  //! Create default finite elements to be used if nothing else is requested
  virtual void SetDefaultElements( PtrParamNode infoNode );

  //! Specialized version of the method of the base class
  //! We number the lowest order dofs consecutively 
  virtual void CreateVirtualNodes();
  
  //! Adjust orders of edges / faces according to min/max rule
  void AdjustEntityOrder();

  //! Virtual node numbers for lowest order edge functions
  boost::array<UInt,2> nedelecNodeRange_;
  
  //! Virtual node numbers for higher order edge functions
  boost::array<UInt,2> edgeHiNodeRange_;
  
  //! Virtual node numbers for higher order face functions
  boost::array<UInt,2> faceNodeRange_;
  
  //! Virtual node numbers for higher order inner functions
  boost::array<UInt,2> innerNodeRange_;
    
  //! Map containing the polynomial order for every edge
  std::map<UInt, StdVector<UInt> > edgeOrder_;
  
  //! Map for reference elements by region
  std::map< RegionIdType, std::map<Elem::FEType, FeHCurlHi* > > refElems_;
  
  //! Map for reference elements by region (1st order)
  std::map< RegionIdType, std::map<Elem::FEType, FeHCurlHi* > > refElems1St_;
  
  //! Flag, if only lowest order should be used
  bool onlyLowestOrder_;
  
  // ====================================================================
  // PROCESS USER INPUT
  // ====================================================================


private:
};
} // end of namespace
#endif // header guard
