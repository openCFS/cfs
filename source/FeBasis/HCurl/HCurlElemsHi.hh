#ifndef FILE_CFS_HCURL_ELEMENTS_HI_HH
#define FILE_CFS_HCURL_ELEMENTS_HI_HH

#include "HCurlElems.hh"
#include <boost/array.hpp>

namespace CoupledField {


//! This is the base class for all HCurl elements of arbitrary order using 
//! hierarchical polynomials. Currently we use the integrated Legendre polynomials,
//! but in general we can plug in any arbitrary orthogonal polynomial (e.h. Jacobi,
//! Gegenbauer, etc.) 
class FeHCurlHi : public FeHCurl {
public:
  //! Constructor
  FeHCurlHi();

  //! Destructor
  virtual ~FeHCurlHi();

  //! Set isotropic polynomial order
  void SetIsoOrder(UInt order);
  
  //! Set usage of gradients
  void UseGradient(EntityType entity, bool usage);
  
  //! Set Usage of only lowest order functions
  void SetOnlyLowestOrder( bool flag) { onlyLowestOrder_ = flag; }
  
  
  //! Get number of shape functions for a given type (NODE/EDGE/FACE/ELEM)
  void GetNumFncs( StdVector<UInt>& numFcns,
                   EntityType fctEntityType,
                   UInt dof = 1 );
  
  //! This holds only for line,quad and hex,
  //! if other types are available the has to be reimplemented!
  //! Get the permutation Vector for a given Face or Edge
  //! e.g. If asked for a face, the element will check the flags
  //! of this face and return a vector of size NumberOfFncs on the Face
  //! holding the correct ordering 
  /*!
  \param fncPermutation (output) The Permuation Vector 
  \param ptElem (input) pointer to Grid Element to get grip of flags 
  \param fctEntityType (input) The Entity type, Node/Edge/Face where the 
  nodes are located at
  \param entNumber (input) The local entity number 
  */
  virtual void GetNodalPermutation( StdVector<UInt>& fncPermutation,
                                    const Elem* ptElem,
                                    EntityType fctEntityType,
                                    UInt entNumber);
  //! \copydoc BaseFE::GetIsoOrder
  virtual UInt GetIsoOrder() const;

  //! \copydoc BaseFE::GetMaxOrder
  virtual UInt GetMaxOrder() const;

  //! \copydoc BaseFE::GetMaxOrderLocDir
  virtual void GetMaxOrderLocDir(StdVector<UInt>& order );
       
protected:
  
  //! Calculate number of unknowns
  virtual void CalcNumUnknowns() = 0;

  //! Isotropic order. 0 if anisotropic
  UInt isoOrder_;
    
  //! Number of shape functions per entity
  std::map<EntityType,StdVector<UInt> > entityFncs_;

  // ========================================================================
  // DEFINITION OF (ANISOTROPIC) ORDER
  // ========================================================================
  //@{ \name Definition of (anisotropic) polynomial order

  //! Polynomial order of edges (#edges x 1 local direction)
  StdVector<UInt> orderEdge_;

  //! Polynomial order of faces (#faces x 2 local directions)
  StdVector<boost::array<UInt,2> > orderFace_;

  //! Polynomial order of inner (1 x 3 local directions)
  boost::array<UInt,3> orderInner_;
  //@}
  
  // ========================================================================
  // Usage Of Gradient function
  // ========================================================================
  
  //! Usage of gradient for given entitytype
  std::map<EntityType, bool> useGrad_;
  
  //! Flag for using only lowst order shape functions
  bool onlyLowestOrder_; 
};



//! HCurl conforming hierarchical higher order quadrilateral element
class  FeHCurlHiQuad : public FeHCurlHi {

public:

  //! Constructor
  FeHCurlHiQuad();

  //! Destructor
  virtual ~FeHCurlHiQuad();

  //! @copydoc FeHCurl::CalcLocShFnc
  void CalcLocShFnc( Matrix<Double>& shape, const LocPointMapped& lp,
                     const Elem* elem, UInt comp = 1 );

  //! @copydoc FeHCurl::CalcLocCurlShFnc
  void CalcLocCurlShFnc( Matrix<Double>& curl, const LocPointMapped& lp,
                         const Elem* elem, UInt comp = 1 );
protected:
  
  //! Calculate number of unknowns
  void CalcNumUnknowns();
  
};



//! HCurl conforming hierarchical higher order hexahedral element
class  FeHCurlHiHex : public FeHCurlHi {

public:

  //! Constructor
  FeHCurlHiHex();

  //! Destructor
  virtual ~FeHCurlHiHex();

  //! @copydoc FeHCurl::CalcLocShFnc
  void CalcLocShFnc( Matrix<Double>& shape, const LocPointMapped& lp,
                     const Elem* elem, UInt comp = 1 );

  //! @copydoc FeHCurl::CalcLocCurlShFnc
  void CalcLocCurlShFnc( Matrix<Double>& curl, const LocPointMapped& lp,
                         const Elem* elem, UInt comp = 1 );

protected:

  //! Calculate number of unknowns
  void CalcNumUnknowns();
};

} // namespace

#endif

