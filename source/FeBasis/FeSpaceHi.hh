#ifndef FILE_CFS_FESPACE_HI_HH
#define FILE_CFS_FESPACE_HI_HH

#include<boost/unordered_set.hpp>
#include<boost/unordered_map.hpp>

#include "FeSpace.hh"

namespace CoupledField {

// forward class declarations
class FeHi;
class Assemble;

#include <unordered_set>
#include <unordered_map>

//! Class for hierarchical finite elements
class FeSpaceHi: public FeSpace {

public:

  // ========================================================================
  //  Helper struct for caching the result of a mapping of coefs to space
  // ========================================================================
  //! Helper struct for caching the result of a function mapping

  //! This struct is used to store the auxiliary data related to the
  //! the mapping of a general coefficient function to this feFunction /
  //! its associated function space.
  struct MapContext {

    //! Constructor
    MapContext( );

    //! Destructor
    virtual ~MapContext( );

    //! Parameter node algebraic system
    PtrParamNode olasNode;

    //! Infonode for algebraic system
    PtrParamNode infoNode;

    //! Pointer to algebraic system
    AlgebraicSys* algSys;

    //! Pointer to assemble class
    Assemble* assemble;

    //! Equation numbers of entity to be mapped
    StdVector<Integer> entityEqns;

    //! Vector containing the solution
    SBM_Vector * sol;
  };

  //! Constructor
  FeSpaceHi( PtrParamNode paramNode, PtrParamNode infoNode, Grid* ptGrid );

  //! Destructor
  virtual ~FeSpaceHi( );

  //! \copydoc FeSpace::MapCoefFctToSpace
  void MapCoefFctToSpace( StdVector<shared_ptr<EntityList> > support,
                          shared_ptr<CoefFunction> coefFct,
                          shared_ptr<BaseFeFunction> feFct,
                          std::map<Integer, Double>& vals, bool cache,
                          const std::set<UInt>& comp = std::set<UInt>(),
                          bool updatedGeo = true );

  //! \copydoc FeSpace::MapCoefFctToSpace
  void MapCoefFctToSpace( StdVector<shared_ptr<EntityList> > support,
                          shared_ptr<CoefFunction> coefFct,
                          shared_ptr<BaseFeFunction> feFct,
                          std::map<Integer, Complex>& vals, bool cache,
                          const std::set<UInt>& comp = std::set<UInt>(),
                          bool updatedGeo = true);

protected:

  template<typename T>
  void MapCoefFctToSpacePriv( StdVector<shared_ptr<EntityList> > support,
                              shared_ptr<CoefFunction> coefFct,
                              shared_ptr<BaseFeFunction> feFct,
                              std::map<Integer, T>& vals, bool cache,
                              const std::set<UInt>& comp = std::set<UInt>() );

  //! Get unmodified higher order finite element
  virtual FeHi* GetFeHi( RegionIdType region, Elem::FEType type ) = 0;

  //! \copydoc FeSpace::MapNodalBCs()
  virtual void MapNodalBCs( );

  //! Map for each region the element order
  std::map<RegionIdType, ApproxOrder> regionOrder_;

  // ====================================================================
  // VARIABLE ENTITY ORDER
  // ====================================================================

  //! Set region order in the beginning
  void SetRegionOrder( RegionIdType region, const ApproxOrder& order );

  //! Initialize a finite element with the correct order
  void SetElemOrder( const Elem* ptEl, FeHi* ptFe, const ApproxOrder& order,
                     bool applyMinMaxRule );

  //! Adjust order of edges / faces which have mixed order neighbours
  void AdjustEntityOrder( );

  //! Set polynomial order of vectorial unknowns for for anisotropic order

  //! In case of vectorial unknowns and anisotropic polynomial order
  //! (= different order for every component of the vector), the
  //! higher order nodes might have to be fixed by a homogeneous Dirichlet
  //! BC.
  void FixHigherOrderAnisoDofs( );

  //! Set containing all edges, where the min/max rule gets applied
  std::unordered_set<UInt> adjustedEdges_;

  //! Set containing all faces, where the min/max rule gets applied
  std::unordered_set<UInt> adjustedFaces_;

  //! Map containing the order of adjusted edges (key: edge number)
  std::unordered_map<UInt, UInt> orderEdges_;

  //! Map containing the order of adjusted faces (key: face number)
  std::unordered_map<UInt, boost::array<UInt, 2> > orderFaces_;

  //! Associate entity name with mapping context
  std::map<std::string, MapContext*> entityCtx_;
};

} // end of namespace
#endif
