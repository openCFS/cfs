#ifndef FILE_CFS_FESPACE_HCURL_HI_HH
#define FILE_CFS_FESPACE_HCURL_HI_HH

#include "FeBasis/FeSpaceHi.hh"
#include <array>
#include <unordered_map>
#include <unordered_set>

#include "Utils/ThreadLocalStorage.hh"

namespace CoupledField {

// forward class declaration
class FeHCurlHi;

//! Finite element space for hierarchical H1 elements
class FeSpaceHCurlHi : public FeSpaceHi {

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
  
  //! Ensure usage of gradient (default: omit)
  virtual void SetUseGradients(RegionIdType region);
  
  //! Treat thin regions specially
  virtual void TreatThinElements(Double maxAspectRatio );
  
  //! \copydoc FeSpace::IsSameEntityApproximation
  virtual bool IsSameEntityApproximation( shared_ptr<EntityList> list,
                                          shared_ptr<FeSpace> space );

protected:

  // ====================================================================
  // INTERNAL INITIALIZATION
  // ====================================================================

  //! Set the order and mapping type of a specific region
  virtual void SetRegionElements( RegionIdType region, 
                                  MappingType mType,
                                  const ApproxOrder& order,
                                  PtrParamNode infoNode );

  //! \copydoc FeSpace::GetNodesOfElement
  virtual void GetNodesOfElement( StdVector<UInt>& nodes,
                                  const Elem* ptElem,
                                  BaseFE::EntityType entType = BaseFE::ALL);
  
  //! Here the spaces have the possibility to check if user definitions makes sense
  //! e.g. if the chosen integration is correct or the element order is nice
  //! here one could e.g. adjust the integration oder according to the element order
  virtual void CheckConsistency();

  //! Sets the default integration scheme and order
  virtual void SetDefaultIntegration( PtrParamNode infoNode );

  //! Create default finite elements to be used if nothing else is requested
  virtual void SetDefaultElements( PtrParamNode infoNode );
  
  //! \copydoc FeSpaceHi::GetFeHi
  virtual FeHi* GetFeHi( RegionIdType region, Elem::FEType type );
  
  //! \copydoc FeSpace::GetNumDofs()
  virtual UInt GetNumDofs() const;

  //! Return groups with faces according for block mapping
  
  //! This method returns the groups of virtual face nodes, which belong
  //! together due to "thin" elements (high aspect ratio).
  void GetThinFaceGroups( StdVector<StdVector<UInt> >& faceGroups );
  
  //! Map for reference elements by region
  std::map< RegionIdType, std::map<Elem::FEType, FeHCurlHi* > > refElems_;
  
  //! Map for reference elements by region (1st order)
  std::map< RegionIdType, std::map<Elem::FEType, FeHCurlHi* > > refElems1St_;
  
  // Cached reference elements used during element matrix computation
  //! Map for reference elements by region
  std::map< RegionIdType, TLMap<Elem::FEType, FeHCurlHi* > > TL_refElems_;

  //! Map for reference elements by region (1st order)
  std::map< RegionIdType, TLMap<Elem::FEType, FeHCurlHi* > > TL_refElems1St_;

  //! Flag, if only lowest order should be used
  bool onlyLowestOrder_;
  
  //! Flag, if anisotropic elements get smoothed in one block
  bool groupAnisoEdges_;
  
  //! Maximum aspect ratio of elements in case of anisotropic element smoothing
  Double maxAspectRatio_;
  
  // ====================================================================
  // GRADIENT HANDLING (only for HCurl)
  // ====================================================================
  
  
  //! Set handling of gradients  for edges / faced interior
  void AdjustGradients( );
  
  //! Set the usage of gradients for one element
  void SetElemGrad( const Elem* ptEl, FeHCurlHi* ptFe, 
                    RegionIdType regionId, bool applyMaxRule );
  
  //! Map for each region, if gradients should be used
  std::map< RegionIdType, bool> useGradients_;
  
  //! Set containing all edges, where use of gradient was adjusted
  std::unordered_set<UInt> adjustedGradEdges_;

  //! Set containing all faces, where use of gradient was adjusted
  std::unordered_set<UInt> adjustedGradFaces_;

  //! Map usage of gradients of adjusted edges (key: edge number)
  std::unordered_map<UInt, bool> gradEdges_;

  //! Map usage of gradients of adjusted faces (key: face number)
  std::unordered_map<UInt, bool> gradFaces_;

  
private:
};
} // end of namespace
#endif // header guard
