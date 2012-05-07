#ifndef FILE_CFS_FESPACE_HI_HH
#define FILE_CFS_FESPACE_HI_HH

#include<boost/unordered_set.hpp>
#include<boost/unordered_map.hpp>

#include "FeSpace.hh"

namespace CoupledField {

// forward class declarations
class FeHi;

//! Class for hierarchical finite elements
class FeSpaceHi : public FeSpace {
public:
  //! Constructor
  FeSpaceHi (PtrParamNode paramNode, PtrParamNode infoNode, Grid* ptGrid );

  //! Destructor
  virtual ~FeSpaceHi();
  
protected:
  
  //! Get unmodified higher order finite element
  virtual FeHi* GetFeHi( RegionIdType region, Elem::FEType type ) = 0;
  
  //! \copydoc FeSpace::MapNodalBCs()
  virtual void MapNodalBCs();

  //! Map for each region the element order
  std::map< RegionIdType, ApproxOrder > regionOrder_;
    
  // ====================================================================
  // VARIABLE ENTITY ORDER
  // ====================================================================
  
  //! Set region order in the beginning
  void SetRegionOrder( RegionIdType region, 
                       const ApproxOrder& order );
  
  //! Initialize a finite element with the correct order
  void SetElemOrder( const Elem* ptEl, FeHi* ptFe, 
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
};

} // end of namespace
#endif
