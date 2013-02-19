#ifndef FILE_CFS_HCURL_ELEMENTS_HI_HH
#define FILE_CFS_HCURL_ELEMENTS_HI_HH

#include "HCurlElems.hh"
#include "FeBasis/FeHi.hh"

namespace CoupledField {


//! This is the base class for all HCurl elements of arbitrary order using 
//! hierarchical polynomials. Currently we use the integrated Legendre polynomials,
//! but in general we can plug in any arbitrary orthogonal polynomial (e.h. Jacobi,
//! Gegenbauer, etc.) 
class FeHCurlHi : public FeHCurl, public FeHi  {
public:

  //! Denotes the derivative type of the shape function to be calculated
  typedef enum { ID = 0, CURL =1} DiffType;
  
  //! Constructor
  FeHCurlHi( Elem::FEType feType  );

  //! Destructor
  virtual ~FeHCurlHi();

  //! Set usage of gradients
  void UseGradient(EntityType entity, bool usage);
  
  //! Set Usage of only lowest order functions
  void SetOnlyLowestOrder( bool flag) { onlyLowestOrder_ = flag; }
  
  //! Get total number of functions
  UInt GetNumFncs();
  
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

protected:
  
  //! Calculate number of unknowns
  virtual void CalcNumUnknowns() = 0;

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


//! HCurl conforming hierarchical higher order wedge element
class  FeHCurlHiWedge : public FeHCurlHi {

public:

  //! Constructor
  FeHCurlHiWedge();

  //! Destructor
  virtual ~FeHCurlHiWedge();

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

  //! Return HCurl shape functions 
  virtual void GetShFnc( Matrix<Double>& shape, 
                         const LocPointMapped& lp,
                         const Elem* elem, UInt comp = 1 );

  //! Return global curl of shape functions
  virtual void GetCurlShFnc( Matrix<Double>& curl, 
                             const LocPointMapped& lp,
                             const Elem* elem, UInt comp = 1 );
  
  //! Second version using the generalized curl-mechanism
  template<DiffType DIFF_TYPE>
  void CalcLocShFnc2( Matrix<Double>& curl, 
                      const LocPointMapped& lp,
                      const Elem* elem, UInt comp = 1 );
  
protected:

  //! Calculate number of unknowns
  void CalcNumUnknowns();
};

} // namespace

#endif

