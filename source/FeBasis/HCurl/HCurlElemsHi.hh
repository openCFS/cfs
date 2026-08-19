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

  //! Copy Constructor
  FeHCurlHi(const FeHCurlHi & other);

  //! Destructor
  virtual ~FeHCurlHi();

  //! Deep Copy
  virtual FeHCurlHi* Clone() = 0;

  //! Set Usage of only lowest order functions
  void SetOnlyLowestOrder( bool flag);
  
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
  //! \copydoc BaseFE::IsIsotropic 
  virtual bool IsIsotropic() {
    return isIsotropic_;
  }
  
  //! \copydoc BaseFE::GetIsoOrder
  virtual UInt GetIsoOrder() const;

  //! \copydoc BaseFE::GetMaxOrder
  virtual UInt GetMaxOrder() const;
  
  //! \copydoc BaseFE::GetAnsiOrder
  virtual void GetAnisoOrder(StdVector<UInt>& order ) const;
  
  //! Compare two element for equality (= same shape and approximation);
  bool operator==( const FeHCurlHi& comp) const;
  
  //! Set general usage of gradient shape functions
  void SetUseGradients(bool useGrad);
  
  //! Set usage of gradient shape functions on specified edge
  void SetEdgeGradient(UInt edgeNum, bool useGrad);
  
  //! Set usage of gradient shape functions on specified face
  void SetFaceGradient(UInt faceNum, bool useGrad);
  
  //! Get gradient usage on edges
  const StdVector<bool>& GetEdgeGradient() const {
    return useEdgeGrad_;
  }
  
  //! Get gradient usage on faces
  const StdVector<bool>& GetFaceGradient() const {
    return useFaceGrad_;
  }
  
protected:
  
  //! Calculate number of unknowns
  virtual void CalcNumUnknowns() = 0;

  //! Flag for using only lowest order shape functions
  bool onlyLowestOrder_; 
  
  // ========================================================================
  // Usage Of Gradient function
  // ========================================================================
  //@{ \name Usage of gradient functions

  //! Flag vector (#edges) if gradient shape functions are used on edges
  StdVector<bool> useEdgeGrad_;
  
  //! Flag vector (#faces) if gradient shape functions are used on face
  StdVector<bool> useFaceGrad_;
  
  //! Flag if gradient shape functions are used in interior
  bool useInteriorGrad_;
  //@}
};


class  FeHCurlHiLine : public FeHCurlHi {

public:

  //! Constructor
  FeHCurlHiLine();

  //! Copy Constructor
  FeHCurlHiLine(const FeHCurlHiLine & other)
    : FeHCurlHi(other){
  }

  //! Create deep copy
  virtual FeHCurlHiLine* Clone(){
    return new FeHCurlHiLine(*this);
  }

  //! Destructor
  virtual ~FeHCurlHiLine();



  //! Return HCurl shape functions
  virtual void GetShFnc( Matrix<Double>& shape,
                         const LocPointMapped& lp,
                         const Elem* elem, UInt comp = 1 )
                         {};

  //! Return global curl of shape functions
  virtual void GetCurlShFnc( Matrix<Double>& curl,
                             const LocPointMapped& lp,
                             const Elem* elem, UInt comp = 1 )
                            {};

  //! Internal method for calculating generalized curl functions
  template<DiffType DIFF_TYPE>
  void CalcLocShFnc2( Matrix<Double>& curl,
                      const LocPointMapped& lp,
                      const Elem* elem, UInt comp = 1 )
                      {};

protected:
  
  //! Calculate number of unknowns
  void CalcNumUnknowns(){};
  
};


//! HCurl conforming hierarchical higher order triangular element
class  FeHCurlHiTria : public FeHCurlHi {

public:

  //! Constructor
  FeHCurlHiTria();

  //! Copy Constructor
  FeHCurlHiTria(const FeHCurlHiTria & other)
    : FeHCurlHi(other){
  }

  //! Create deep copy
  virtual FeHCurlHiTria* Clone(){
    return new FeHCurlHiTria(*this);
  }

  //! Destructor
  virtual ~FeHCurlHiTria();


  //! Return HCurl shape functions
  virtual void GetShFnc( Matrix<Double>& shape,
                         const LocPointMapped& lp,
                         const Elem* elem, UInt comp = 1 );

  //! Return global curl of shape functions
  virtual void GetCurlShFnc( Matrix<Double>& curl,
                             const LocPointMapped& lp,
                             const Elem* elem, UInt comp = 1 );

  //! Internal method for calculating generalized curl functions
  template<DiffType DIFF_TYPE>
  void CalcLocShFnc2( Matrix<Double>& curl,
                      const LocPointMapped& lp,
                      const Elem* elem, UInt comp = 1 );

protected:
  
  //! Calculate number of unknowns
  void CalcNumUnknowns();
  
};


//! HCurl conforming hierarchical higher order quadrilateral element
class  FeHCurlHiQuad : public FeHCurlHi {

public:

  //! Constructor
  FeHCurlHiQuad();

  //! Copy Constructor
  FeHCurlHiQuad(const FeHCurlHiQuad & other)
    : FeHCurlHi(other){
  }

  //! Destructor
  virtual ~FeHCurlHiQuad();

  //! Create deep copy
  virtual FeHCurlHiQuad* Clone(){
    return new FeHCurlHiQuad(*this);
  }

  //! Return HCurl shape functions 
  virtual void GetShFnc( Matrix<Double>& shape, 
                         const LocPointMapped& lp,
                         const Elem* elem, UInt comp = 1 );

  //! Return global curl of shape functions
  virtual void GetCurlShFnc( Matrix<Double>& curl, 
                             const LocPointMapped& lp,
                             const Elem* elem, UInt comp = 1 );
  
  //! Internal method for calculating generalized curl functions
  template<DiffType DIFF_TYPE>
  void CalcLocShFnc2( Matrix<Double>& curl, 
                      const LocPointMapped& lp,
                      const Elem* elem, UInt comp = 1 );
protected:
  
  //! Calculate number of unknowns
  void CalcNumUnknowns();
  
};



// ============
//  HEXAHEDRAL 
// ============

//! HCurl conforming hierarchical higher order hexahedral element
class  FeHCurlHiHex : public FeHCurlHi {

public:

  //! Constructor
  FeHCurlHiHex();

  //! Copy Constructor
  FeHCurlHiHex(const FeHCurlHiHex & other)
    : FeHCurlHi(other){
  }

  //! Create deep copy
  virtual FeHCurlHiHex* Clone(){
    return new FeHCurlHiHex(*this);
  }

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

// =======
//  WEDGE
// =======
//! HCurl conforming hierarchical higher order wedge / prismatic element
class  FeHCurlHiWedge : public FeHCurlHi {

public:

  //! Constructor
  FeHCurlHiWedge();

  //! Copy Constructor
  FeHCurlHiWedge(const FeHCurlHiWedge & other)
    : FeHCurlHi(other){
  }

  //! Create deep copy
  virtual FeHCurlHiWedge* Clone(){
    return new FeHCurlHiWedge(*this);
  }


  //! Destructor
  virtual ~FeHCurlHiWedge();

  //! Return HCurl shape functions 
  virtual void GetShFnc( Matrix<Double>& shape, 
                         const LocPointMapped& lp,
                         const Elem* elem, UInt comp = 1 );

  //! Return global curl of shape functions
  virtual void GetCurlShFnc( Matrix<Double>& curl, 
                             const LocPointMapped& lp,
                             const Elem* elem, UInt comp = 1 );
  
  //! Internal method for calculating generalized curl functions
  template<DiffType DIFF_TYPE>
  void CalcLocShFnc2( Matrix<Double>& curl, 
                      const LocPointMapped& lp,
                      const Elem* elem, UInt comp = 1 );
protected:

  //! Calculate number of unknowns
  void CalcNumUnknowns();
};

// =======
//  TET  
// =======
//! HCurl conforming hierarchical higher order tetrahedral element
class  FeHCurlHiTet : public FeHCurlHi {

public:

  //! Constructor
  FeHCurlHiTet();

  //! Copy Constructor
  FeHCurlHiTet(const FeHCurlHiTet & other)
    : FeHCurlHi(other){
  }

  //! Create deep copy
  virtual FeHCurlHiTet* Clone(){
    return new FeHCurlHiTet(*this);
  }

  //! Destructor
  virtual ~FeHCurlHiTet();

  //! Return HCurl shape functions 
  virtual void GetShFnc( Matrix<Double>& shape, 
                         const LocPointMapped& lp,
                         const Elem* elem, UInt comp = 1 );

  //! Return global curl of shape functions
  virtual void GetCurlShFnc( Matrix<Double>& curl, 
                             const LocPointMapped& lp,
                             const Elem* elem, UInt comp = 1 );
  
  //! Internal method for calculating generalized curl functions
  template<DiffType DIFF_TYPE>
  void CalcLocShFnc2( Matrix<Double>& curl, 
                      const LocPointMapped& lp,
                      const Elem* elem, UInt comp = 1 );
protected:

  //! Calculate number of unknowns
  void CalcNumUnknowns();
};


// =======
//  PYRA
// =======
//! HCurl conforming hierarchical higher order pyramidal element
class  FeHCurlHiPyra : public FeHCurlHi {

public:

  //! Constructor
  FeHCurlHiPyra();

  //! Copy Constructor
  FeHCurlHiPyra(const FeHCurlHiPyra & other)
    : FeHCurlHi(other){
  }

  //! Create deep copy
  virtual FeHCurlHiPyra* Clone(){
    return new FeHCurlHiPyra(*this);
  }

  //! Destructor
  virtual ~FeHCurlHiPyra();

  //! Return HCurl shape functions 
  virtual void GetShFnc( Matrix<Double>& shape, 
                         const LocPointMapped& lp,
                         const Elem* elem, UInt comp = 1 );

  //! Return global curl of shape functions
  virtual void GetCurlShFnc( Matrix<Double>& curl, 
                             const LocPointMapped& lp,
                             const Elem* elem, UInt comp = 1 );
  
  //! Internal method for calculating generalized curl functions
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

