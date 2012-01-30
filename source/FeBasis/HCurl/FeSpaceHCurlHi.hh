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
  FeSpaceHCurlHi(PtrParamNode aNode, PtrParamNode infoNode);

  //! Destructor
  ~FeSpaceHCurlHi();
  
  //! Initialize class (read order etc.)
  void Init();

  //! \copydoc FeSpace::GetFe(EntityIterator,shared_ptr<IntScheme>&)
  virtual BaseFE* GetFe( const EntityIterator ent ,
                         shared_ptr<IntScheme>& intScheme );

  //! \copydoc FeSpace::GetFe(EntityIterator)
  virtual BaseFE* GetFe( const EntityIterator ent );

  //! \copydoc FeSpace::GetFe(UInt)
  virtual BaseFE* GetFe( UInt elemNum );


  //! \copydoc FeSpace::GetNumEntityOrder
  UInt GetEntityOrder( UInt elemNum, BaseFE::EntityType type, 
                       UInt entityNum, UInt comp = 1 );

  //! \copydoc FeSpace::GetMaxEntityOrder
  UInt GetMaxEntityOrder( UInt elemNum, BaseFE::EntityType type, 
                          UInt entityNum );

  //! Map equations i.e. initialize object
  virtual void Finalize();
  
  //! Return SBM-block and Matrix-SubBlockdefinition according to strategy
  virtual void GetOlasMappings( shared_ptr<SolStrategy> solStrat, 
                                StdVector<AlgebraicSys::SBMBlockDef>& sbmBlocks,
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

  // ====================================================================
  // PROCESS USER INPUT
  // ====================================================================


private:
};
} // end of namespace
#endif // header guard
