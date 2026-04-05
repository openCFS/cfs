#ifndef FILE_FE_HI_HH
#define FILE_FE_HI_HH

#include <array>

#include "BaseFE.hh"
#include "Domain/ElemMapping/Elem.hh"

namespace CoupledField {

//! Base class for all hierarchical finite elements 

//! This class encapsulates the functionality of hierarchic, arbitrary order
//! finite elements. It implements the assignment and query of polynomial
//! order (isotropic or anisotropic) to the different entities (edges /
//! faces / interior). 
class FeHi{
public:

  //! Constructor
  FeHi( Elem::FEType feType  );
  
  //! Copy constructor
  FeHi(const FeHi& other);

  //! Destructor
  virtual ~FeHi();
  
  //! Deep Copy
  virtual FeHi* Clone() = 0;

  // ----------------------------------------------------------------------
  //  Polynomial Order
  // ----------------------------------------------------------------------
  //! Set isotropic polynomial order
  void SetIsoOrder( UInt order );

  //! Set anisotropic / directional dependent order

  //! Set the polynomial order direction dependent, i.e. different for each 
  //! direction. This can be used to prefer certain directions, i.e. for 
  //! thin structures or anisotropic field distribution.
  //! \param order Vector containing the polynomial order in each local 
  //!              direction (xi,eta,zeta).
  //! \param ptElem Pointer to geometric element
  //!
  //! \note Only simplex-based elements (quad, wedge, hex) support the 
  //!       prescription of anisotropic polynomial order.

  void SetAnisoOrder( const StdVector<UInt>& order,
                      const Elem* ptElem );

  //! Return set of basis functions, which exceeds a given polynomial order

  //! This method returns a set of indices of basis functions, which have
  //! a higher polynomial degree than the (spatially dependent) order
  //! provided as argument. The indices correspond to the positions of the
  //! basis functions / derivatives obtainable e.g. by the method
  //! FeH1::GetShFnc.
  //! \param nodes Contains all indices of basis functions, which have higher
  //!              polynomial degree than the one provided in the order
  //!              vector
  //! \param order Vector containing the polynomial order in each local 
  //!              direction (xi,eta,zeta) (see FeH1Hi::SetAnisoOrder)
  //! \param ptElem Pointer to geometric element
  void GetNodesExceedingOrder( std::set<UInt>& nodes, 
                               const StdVector<UInt>& order,
                               const Elem* ptElem );

  //! Set polynomial order for an edge

  //! Sets the order of a given edge number (0-based) to a given order. 
  void SetEdgeOrder( UInt edgeNum, UInt order  );

  //! Set polynomial order for a face

  //! Sets the order of a given face number (0-based) to a given order.
  //! The order is specified in terms of face-local directions.
  void SetFaceOrder( UInt faceNum, const std::array<UInt,2>& order );

  //! Set polynomial order for element interior
  void SetInteriorOrder( const std::array<UInt,3>& order ); 

  //! Return edge order
  const StdVector<UInt>& GetEdgeOrder( ) const {
    return orderEdge_; 
  }

  //! Return face order 
  const StdVector<std::array<UInt,2> >& GetFaceOrder( ) const {
    return orderFace_;
  }

  //! Return interior order
  const std::array<UInt,3>& GetInnerOrder( ) const {
    return orderInner_;
  }
  
protected:
  
  //! Calculate number of unknowns

  //! This method calculates the number of unknowns functions
  //! for real tensor-product elements, i.e. line, quadrilateral
  //! and hexahedral elements.
  virtual void CalcNumUnknowns() = 0;

  //! Return matrix with polynomial degrees for each shape function

  //! This method returns a matrix, which contains for every basis function
  //! and every local direction the polynomial degree. This can be used e.g.
  //! to determine the spatial order of a given node.
  virtual void GetPolyOrderOfNodes( Matrix<UInt>& polyOrder,
                                    const Elem* ptElem ) {
    EXCEPTION( "Not implemented" );
  }

  //! Element shape of reference element
  ElemShape elemShape_;
  
  //! Geometric type of finite element (line, quad, hex ...)
  Elem::FEType myFeType_ = Elem::ET_UNDEF;
  
  //! Flag if re-calculation of number of unknowns is needed

  //! After changing the order of the element, a re-calculation of 
  //! the number of unknowns (actNumFncs_) is necessary.
  bool updateUnknowns_ = false;

  //! Flag for isotropic polynomial order
  bool isIsotropic_ = false;
  
  //! Number of shape functions per entity
  std::map<BaseFE::EntityType,StdVector<UInt> > entityFncs_;

  //! Isotropic order. 0 if anisotropic
  UInt isoOrder_ = 0;
  
  //! Directional-dependent anisotropic order (w.r.t. to local directions)
  StdVector<UInt> anisoOrder_;

  //! Maximum polynomial degree of element
  UInt maxOrder_ = 0;

  // ========================================================================
  // DEFINITION OF (ANISOTROPIC) ORDER
  // ========================================================================
  //@{ \name Definition of (anisotropic) polynomial order

  //! Polynomial order of edges (#edges x 1 local direction)
  StdVector<UInt> orderEdge_;

  //! Polynomial order of faces (#faces x 2 local directions)
  StdVector<std::array<UInt,2> > orderFace_;

  //! Polynomial order of inner (1 x 3 local directions)
  std::array<UInt,3> orderInner_;
  //@}
}; // FeHi class

} // end of namespace


#endif // header guard
