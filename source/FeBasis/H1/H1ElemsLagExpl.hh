#ifndef FILE_CFS_H1_ELEMENTS_LAGRANGE_EXPLICIT_HH
#define FILE_CFS_H1_ELEMENTS_LAGRANGE_EXPLICIT_HH

#include <array>
#include "H1Elems.hh"
#include "FeBasis/FeNodal.hh"

namespace CoupledField {

// ========================================================================
//  H1 Fe Lagrangian Elements of lowest order (1st / 2nd)
// ========================================================================

//! Lagrangian based elements of lowest order (explicitly coded)

//! This serves as a base class for all Lagrangian bases Elements of 
//! lowest order, i.e. the  shape functions can be pre-calculated without 
//! knowledge of the global element.
//! Currently these elements are all 1st / 2nd order Lagrangian elements
class FeH1LagrangeExpl : public FeH1, public FeNodal {

public:
  //! Constructor
  FeH1LagrangeExpl();

  //! Copy Constructor
  FeH1LagrangeExpl(const FeH1LagrangeExpl & other);

  //! Destructor
  virtual ~FeH1LagrangeExpl();

  //! Get number of shape functions for a given type (NODE/EDGE/FACE/ELEM)
  void GetNumFncs( StdVector<UInt>& numFcns,
                   EntityType fctEntityType,
                   UInt dof = 1 );

  //! Get the permutation Vector for a given Face or Edge
  //! e.g. If asked for a face, the element will check the flags
  //! of this face and return a vector of size NumberOfFncs on the Face
  //! holding the correct ordering. BUT: this only works for tensor product elements!
  //! For Explicit elements this needs to be overloaded!
  /*!
  \param fncPermutation (output) The Permuation Vector 
  \param ptElem (input) pointer to Grid Element to get grip of flags 
  \param fctEntityType (input) The Entity type, Node/Edge/Face
  \param entNumber (input) The local entity number 
  */
  virtual void GetNodalPermutation( StdVector<UInt>& fncPermutation,
                                  const Elem* ptElem,
                                  EntityType fctEntityType,
                                  UInt entNumber);

  //! \copydoc BaseFE::GetIsoOrder
  virtual UInt GetIsoOrder() const {
    return order_;
  }
  //! \copydoc BaseFE::GetMaxOrder
  virtual UInt GetMaxOrder() const {
    return order_;
  }
  
  //! Compare two element for equality (= same shape and approximation);
  bool operator==( const FeH1LagrangeExpl& comp) const;
  
  //! =======================================================================
  //! G E O M E T R I C   I N F O R M A T I O N
  //! =======================================================================
  //! Returns whether a given local coordinate is within this element
  //! @param point input Local point
  //! @param tolerance input Additioanl (relative) tolerance
  //! \return flag if point is inside the element
  virtual bool CoordIsInsideElem( const Vector<Double>& point,
                                  Double tolerance ) = 0;
  
  //! Calculates corresponding volume point of neighboring surfaces

  //! For a given surface element and a neighboring volume element this
  //! method calculates the local volume-coordinates out of the given
  //! local surface-coordinates, which have one less dimension.
  //! This can be used to get the corresponding volume coordinates of
  //! the integration points of a surface. Therefore it calculates
  //! on which side of the volume element the surface element lies
  //! and creates the according volume point.
  //! 
  //! \note The resulting normal vector will always point OUT of the
  //!       first volume neighbor of the surface element 
  //!       (SurfElem.ptVolElems[0]).
  /*!
       \param surfConnect (input) Node numbers of surface element
       \param volConnect (input) Node numbers of volume element
       \param surfIntPoint (input) Surface integration point, which gets mapped
       onto the volume element
       \param volIntPoint (output) Corresponding volume integration point
       \param locNormal (output) Normal direction of surface element in local 
       coordinates of the volume element, pointint OUT of the first volume element
   */
  virtual void GetLocalIntPoints4Surface(const StdVector<UInt> & surfConnect,
                                         const StdVector<UInt> & volConnect,
                                         const LocPoint & surfIntPoint,
                                         LocPoint & volIntPoint,
                                         Vector<Double>& locNormal ) = 0;

  //! @copydoc BaseFE::ComputeMonomialCoefficients
  //! Overloaded method for lagrange Elements
  //virtual void ComputeMonomialCoefficients(Matrix<Integer>& P, Matrix<Double>& C);

  //! @copydoc BaseFE::GetLocalDOFCoordinates
  virtual void GetLocalDOFCoordinates(Matrix<Double> & coordMat);
  
  //! return index of triangulation
  virtual void Triangulate(StdVector< StdVector<UInt> > & triConnect){
    EXCEPTION("Triangulated not available for this element.");
  }

protected:

  //! @copydoc FeNodal::SetFunctionsAtIp
  void SetFunctionsAtIp( const StdVector<LocPoint>& iPoints );

  //! @copydoc FeNodal::SetFunctionsAtIp
  void SetFunctionsAtIp( const std::map<Integer, LocPoint >& 
                         iPoints);

  //! @copydoc FeH1::CalcShFnc
  virtual void CalcShFnc( Vector<Double>& shape,
                          const Vector<Double>& point,
                          const Elem* ptElem,
                          UInt comp = 1 ) = 0;

  //! @copydoc FeH1::CalcLocDerivShFnc
  virtual void CalcLocDerivShFnc( Matrix<Double> & deriv, 
                                  const Vector<Double>& point,
                                  const Elem* ptElem,
                                  UInt comp = 1 ) = 0;

  /*! mapping function that checks the orientation of quad elements.
    The function is used in GetLocalIntPoints4Surface() for the Hex, Wedge and Pyra elements.
    It is put here into the base class to avoid code duplication.
    \param commonIndexMap (input) maps the vol elem indices of the surface to the unitary [1,2,3,4] indices
    \param commonIndex (input) common indices
    \param surfIntPoint (input) surface integration point
    \param volIntPoint (output) volume integration point
  */
  void MapQuadSurfOrientation(const std::map<UInt, UInt>& commonIndexMap, 
                              const StdVector<UInt>& commonIndex, 
                              const LocPoint& surfIntPoint, 
                              std::array<double, 2>& volIntPoint);

  /*! mapping function that checks the orientation of triangular elements.
    The function is used in GetLocalIntPoints4Surface() for the Tet, Wedge and Pyra elements.
    It is put here into the base class to avoid code duplication.
    \param commonIndexMap (input) maps the vol elem indices of the surface to the unitary [1,2,3] indices
    \param commonIndex (input) common indices
    \param surfIntPoint (input) surface integration point
    \param volIntPoint (output) volume integration point
    \param isEquilateral (input) flag to indicate if the triangle is equilateral, as occurs in one of the tetrhedron faces
  */
  void MapTriaSurfOrientation(const std::map<UInt, UInt>& commonIndexMap, 
                              const StdVector<UInt>& commonIndex, 
                              const LocPoint& surfIntPoint, 
                              std::array<double, 2>& volIntPoint,
                              bool isEquilateral = false);

  //! Polynomial order of the finite element
  UInt order_;
};


// ========================================================================
// Common base classes for elements of 1st and 2nd order
// ========================================================================

//! Lagrangian line element
class FeH1LagrangeLine : public FeH1LagrangeExpl {

public:

  //! Constructor 
  FeH1LagrangeLine() : FeH1LagrangeExpl() {};

  //! Copy Constructor
  FeH1LagrangeLine(const FeH1LagrangeLine& other)
   : FeH1LagrangeExpl(other){
  }

  //! Destructor
  virtual ~FeH1LagrangeLine() {};

protected:

  //! @copydoc FeH1::CalcShFnc
  void CalcShFnc( Vector<Double>& shape,
                  const Vector<Double>& point,
                  const Elem* ptElem,
                  UInt comp = 1 ) = 0;

  //! @copydoc FeH1::CalcLocDerivShFnc
  void CalcLocDerivShFnc( Matrix<Double> & deriv, 
                          const Vector<Double>& point,
                          const Elem* ptElem,
                          UInt comp = 1 ) = 0;
  
  //! @copydoc FeH1LagrangeExpl::CoordIsInsideElem
  bool CoordIsInsideElem( const Vector<Double>& point,
                          Double tolerance );
  
  //! @copydoc FeH1LagrangeExpl::GetLocalIntPoints4Surface
  void GetLocalIntPoints4Surface(const StdVector<UInt> & surfConnect,
                                 const StdVector<UInt> & volConnect,
                                 const LocPoint & surfIntPoint,
                                 LocPoint & volIntPoint,
                                 Vector<Double>& locNormal );
};

//! Lagrangian triangular element
class FeH1LagrangeTria : public FeH1LagrangeExpl {

public:

  //! Constructor
  FeH1LagrangeTria() : FeH1LagrangeExpl() {};

  //! Copy Constructor
  FeH1LagrangeTria(const FeH1LagrangeTria& other)
   : FeH1LagrangeExpl(other){
  }

  //! Destructor
  virtual ~FeH1LagrangeTria() {};

protected:

  //! @copydoc FeH1::CalcShFnc
  void CalcShFnc( Vector<Double>& shape,
                  const Vector<Double>& point,
                  const Elem* ptElem,
                  UInt comp = 1 ) = 0;

  //! @copydoc FeH1::CalcLocDerivShFnc
  void CalcLocDerivShFnc( Matrix<Double> & deriv, 
                          const Vector<Double>& point,
                          const Elem* ptElem,
                          UInt comp = 1 ) = 0;
  
  //! @copydoc FeH1LagrangeExpl::CoordIsInsideElem
  bool CoordIsInsideElem( const Vector<Double>& point,
                          Double tolerance );
  
  //! @copydoc FeH1LagrangeExpl::GetLocalIntPoints4Surface
  void GetLocalIntPoints4Surface(const StdVector<UInt> & surfConnect,
                                 const StdVector<UInt> & volConnect,
                                 const LocPoint & surfIntPoint,
                                 LocPoint & volIntPoint,
                                 Vector<Double>& locNormal );

  //! @copydoc BaseFE::ComputeMonomialCoefficients
  //! Overloaded method for lagrange Elements
  virtual void ComputeMonomialCoefficients(Matrix<Integer>& P, Matrix<Double>& C);

};

//! Lagrangian quadrilateral element
class FeH1LagrangeQuad : public FeH1LagrangeExpl {

public:

  //! Constructor
  FeH1LagrangeQuad() : FeH1LagrangeExpl() {};

  //! Copy Constructor
  FeH1LagrangeQuad(const FeH1LagrangeQuad& other)
   : FeH1LagrangeExpl(other){
  }

  //! Destructor
  virtual ~FeH1LagrangeQuad() {};

protected:

  //! @copydoc FeH1::CalcShFnc
  void CalcShFnc( Vector<Double>& shape,
                  const Vector<Double>& point,
                  const Elem* ptElem,
                  UInt comp = 1 ) = 0;

  //! @copydoc FeH1::CalcLocDerivShFnc
  void CalcLocDerivShFnc( Matrix<Double> & deriv, 
                          const Vector<Double>& point,
                          const Elem* ptElem,
                          UInt comp = 1 ) = 0;
  
  //! @copydoc FeH1LagrangeExpl::CoordIsInsideElem
  bool CoordIsInsideElem( const Vector<Double>& point,
                          Double tolerance );
  
  //! @copydoc FeH1LagrangeExpl::GetLocalIntPoints4Surface
  void GetLocalIntPoints4Surface(const StdVector<UInt> & surfConnect,
                                 const StdVector<UInt> & volConnect,
                                 const LocPoint & surfIntPoint,
                                 LocPoint & volIntPoint,
                                 Vector<Double>& locNormal );

};

//! Lagrangian hexahedral element
class FeH1LagrangeHex : public FeH1LagrangeExpl {

public:

  //! Constructor  
  FeH1LagrangeHex() : FeH1LagrangeExpl() {};

  //! Copy Constructor
  FeH1LagrangeHex(const FeH1LagrangeHex& other)
   : FeH1LagrangeExpl(other){
  }

  //! Destructor
  virtual  ~FeH1LagrangeHex() {};

protected:

  //! @copydoc FeH1::CalcShFnc
  void CalcShFnc( Vector<Double>& shape,
                  const Vector<Double>& point,
                  const Elem* ptElem,
                  UInt comp = 1 ) = 0;

  //! @copydoc FeH1::CalcLocDerivShFnc
  void CalcLocDerivShFnc( Matrix<Double> & deriv, 
                          const Vector<Double>& point,
                          const Elem* ptElem,
                          UInt comp = 1 ) = 0;
  
  //! @copydoc FeH1LagrangeExpl::CoordIsInsideElem
  bool CoordIsInsideElem( const Vector<Double>& point,
                          Double tolerance );
  
  //! @copydoc FeH1LagrangeExpl::GetLocalIntPoints4Surface
  void GetLocalIntPoints4Surface(const StdVector<UInt> & surfConnect,
                                 const StdVector<UInt> & volConnect,
                                 const LocPoint & surfIntPoint,
                                 LocPoint & volIntPoint,
                                 Vector<Double>& locNormal );

};


//! Lagrangian wedge element
class FeH1LagrangeWedge : public FeH1LagrangeExpl {

public:

  //! Constructor  
  FeH1LagrangeWedge() : FeH1LagrangeExpl() {};

  //! Copy Constructor
  FeH1LagrangeWedge(const FeH1LagrangeWedge& other)
   : FeH1LagrangeExpl(other){
  }

  //! Destructor
  virtual  ~FeH1LagrangeWedge() {};

protected:

  //! @copydoc FeH1::CalcShFnc
  void CalcShFnc( Vector<Double>& shape,
                  const Vector<Double>& point,
                  const Elem* ptElem,
                  UInt comp = 1 ) = 0;

  //! @copydoc FeH1::CalcLocDerivShFnc
  void CalcLocDerivShFnc( Matrix<Double> & deriv, 
                          const Vector<Double>& point,
                          const Elem* ptElem,
                          UInt comp = 1 ) = 0;

  //! @copydoc FeH1LagrangeExpl::CoordIsInsideElem
  bool CoordIsInsideElem( const Vector<Double>& point,
                          Double tolerance );

  //! @copydoc FeH1LagrangeExpl::GetLocalIntPoints4Surface
  void GetLocalIntPoints4Surface(const StdVector<UInt> & surfConnect,
                                 const StdVector<UInt> & volConnect,
                                 const LocPoint & surfIntPoint,
                                 LocPoint & volIntPoint,
                                 Vector<Double>& locNormal );

};

//! Lagrangian tetrahedron element
class FeH1LagrangeTet : public FeH1LagrangeExpl {

public:

  //! Constructor
  FeH1LagrangeTet() : FeH1LagrangeExpl() {};

  //! Copy Constructor
  FeH1LagrangeTet(const FeH1LagrangeTet& other)
   : FeH1LagrangeExpl(other){
  }

  //! Destructor
  virtual  ~FeH1LagrangeTet() {};

protected:

  //! @copydoc FeH1::CalcShFnc
  void CalcShFnc( Vector<Double>& shape,
                  const Vector<Double>& point,
                  const Elem* ptElem,
                  UInt comp = 1 ) = 0;

  //! @copydoc FeH1::CalcLocDerivShFnc
  void CalcLocDerivShFnc( Matrix<Double> & deriv,
                          const Vector<Double>& point,
                          const Elem* ptElem,
                          UInt comp = 1 ) = 0;

  //! @copydoc FeH1LagrangeExpl::CoordIsInsideElem
  bool CoordIsInsideElem( const Vector<Double>& point,
                          Double tolerance );

  //! @copydoc FeH1LagrangeExpl::GetLocalIntPoints4Surface
  void GetLocalIntPoints4Surface(const StdVector<UInt> & surfConnect,
                                 const StdVector<UInt> & volConnect,
                                 const LocPoint & surfIntPoint,
                                 LocPoint & volIntPoint,
                                 Vector<Double>& locNormal );

};

//! Lagrangian pyramid element
class FeH1LagrangePyra : public FeH1LagrangeExpl {

public:

  //! Constructor
  FeH1LagrangePyra() : FeH1LagrangeExpl()  {};

  //! Copy Constructor
  FeH1LagrangePyra(const FeH1LagrangePyra& other)
   : FeH1LagrangeExpl(other){
  }

  //! Destructor
  virtual  ~FeH1LagrangePyra() {};

protected:

  //! @copydoc FeH1::CalcShFnc
  void CalcShFnc( Vector<Double>& shape,
                  const Vector<Double>& point,
                  const Elem* ptElem,
                  UInt comp = 1 ) = 0;

  //! @copydoc FeH1::CalcLocDerivShFnc
  void CalcLocDerivShFnc( Matrix<Double> & deriv,
                          const Vector<Double>& point,
                          const Elem* ptElem,
                          UInt comp = 1 ) = 0;

  //! @copydoc FeH1LagrangeExpl::CoordIsInsideElem
  bool CoordIsInsideElem( const Vector<Double>& point,
                          Double tolerance );

  //! @copydoc FeH1LagrangeExpl::GetLocalIntPoints4Surface
  void GetLocalIntPoints4Surface(const StdVector<UInt> & surfConnect,
                                 const StdVector<UInt> & volConnect,
                                 const LocPoint & surfIntPoint,
                                 LocPoint & volIntPoint,
                                 Vector<Double>& locNormal );

};

// ========================================================================
// Lagrange H1 elements of 1st order
// ========================================================================

//! Lagrangian line element of 1st order (ET_LINE2)
class FeH1LagrangeLine1 : public FeH1LagrangeLine {

public:

  //! Constructor 
  FeH1LagrangeLine1();

  //! Copy Constructor
  FeH1LagrangeLine1(const FeH1LagrangeLine1& other)
   : FeH1LagrangeLine(other){
  }

  //! Destructor
  virtual ~FeH1LagrangeLine1();

  //! Create deep copy
  virtual FeH1LagrangeLine1* Clone(){
    return new FeH1LagrangeLine1(*this);
  }
protected:

  //! @copydoc FeH1::CalcShFnc
  void CalcShFnc( Vector<Double>& shape,
                  const Vector<Double>& point,
                  const Elem* ptElem,
                  UInt comp = 1 );

  //! @copydoc FeH1::CalcLocDerivShFnc
  void CalcLocDerivShFnc( Matrix<Double> & deriv, 
                          const Vector<Double>& point,
                          const Elem* ptElem,
                          UInt comp = 1 );
};

//! Lagrangian triangular element of 1st order (ET_TRIA3)
class FeH1LagrangeTria1 : public FeH1LagrangeTria {

public:

  //! Constructor
  FeH1LagrangeTria1();

  //! Copy Constructor
  FeH1LagrangeTria1(const FeH1LagrangeTria1& other)
   : FeH1LagrangeTria(other){
  }

  //! Destructor
  virtual ~FeH1LagrangeTria1();

  //! Create deep copy
  virtual FeH1LagrangeTria1* Clone(){
    return new FeH1LagrangeTria1(*this);
  }

protected:

  //! @copydoc FeH1::CalcShFnc
  void CalcShFnc( Vector<Double>& shape,
                  const Vector<Double>& point,
                  const Elem* ptElem,
                  UInt comp = 1 );

  //! @copydoc FeH1::CalcLocDerivShFnc
  void CalcLocDerivShFnc( Matrix<Double> & deriv, 
                          const Vector<Double>& point,
                          const Elem* ptElem,
                          UInt comp = 1 );

  //! @copydoc BaseFE::Triangulate
  virtual void Triangulate(StdVector< StdVector<UInt> > & triConnect);
};

//! Lagrangian quadrilateral element of 1st order (ET_QUAD4)
class FeH1LagrangeQuad1 : public FeH1LagrangeQuad {

public:

  //! Constructor
  FeH1LagrangeQuad1();

  //! Copy Constructor
  FeH1LagrangeQuad1(const FeH1LagrangeQuad1& other)
   : FeH1LagrangeQuad(other){
  }

  //! Destructor
  virtual ~FeH1LagrangeQuad1();

  //! Create deep copy
  virtual FeH1LagrangeQuad1* Clone(){
    return new FeH1LagrangeQuad1(*this);
  }


  virtual void ComputeMonomialCoefficients(Matrix<Integer>& P, Matrix<Double>& C);
protected:

  //! @copydoc FeH1::CalcShFnc
  void CalcShFnc( Vector<Double>& shape,
                  const Vector<Double>& point,
                  const Elem* ptElem,
                  UInt comp = 1 );

  //! @copydoc FeH1::CalcLocDerivShFnc
  void CalcLocDerivShFnc( Matrix<Double> & deriv, 
                          const Vector<Double>& point,
                          const Elem* ptElem,
                          UInt comp = 1 );
  

  //! @copydoc FeH1::GetNumICModes
  virtual UInt GetNumICModes() { return 2;}
    
  //! @copydoc FeH1::CalcShFncICModes
  void CalcShFncICModes( Vector<Double>& shape,
                         const Vector<Double>& point,
                         const Elem* ptElem,
                         UInt comp = 1 );

  //! @copydoc FeH1::CalcLocDerivShFncICModes
  void CalcLocDerivShFncICModes( Matrix<Double> & deriv, 
                                 const Vector<Double>& point,
                                 const Elem* ptElem,
                                 UInt comp = 1 );


  //! @copydoc BaseFE::Triangulate
  virtual void Triangulate(StdVector< StdVector<UInt> > & triConnect);
};

//! Lagrangian hexahedral element of 1st order (ET_HEXA8)
class FeH1LagrangeHex1 : public FeH1LagrangeHex {

public:

  //! Constructor  
  FeH1LagrangeHex1();

  //! Copy Constructor
  FeH1LagrangeHex1(const FeH1LagrangeHex1& other)
   : FeH1LagrangeHex(other){
  }

  //! Destructor
  virtual  ~FeH1LagrangeHex1();

  //! Create deep copy
  virtual FeH1LagrangeHex1* Clone(){
    return new FeH1LagrangeHex1(*this);
  }

  void ComputeMonomialCoefficients(Matrix<Integer>& P, Matrix<Double>& C);

protected:

  //! @copydoc FeH1::CalcShFnc
  void CalcShFnc( Vector<Double>& shape,
                  const Vector<Double>& point,
                  const Elem* ptElem,
                  UInt comp = 1 );

  //! @copydoc FeH1::CalcLocDerivShFnc
  void CalcLocDerivShFnc( Matrix<Double> & deriv, 
                          const Vector<Double>& point,
                          const Elem* ptElem,
                          UInt comp = 1 );
  
  //! @copydoc FeH1::GetNumICModes
  virtual UInt GetNumICModes() { return 3;}
  
  //! @copydoc FeH1::CalcShFncICModes
  void CalcShFncICModes( Vector<Double>& shape,
                         const Vector<Double>& point,
                         const Elem* ptElem,
                         UInt comp = 1 );

  //! @copydoc FeH1::CalcLocDerivShFncICModes
  void CalcLocDerivShFncICModes( Matrix<Double> & deriv, 
                                 const Vector<Double>& point,
                                 const Elem* ptElem,
                                 UInt comp = 1 );

  //! @copydoc BaseFE::Triangulate
  virtual void Triangulate(StdVector< StdVector<UInt> > & triConnect);
};


//! Lagrangian wedge element of 1st order (ET_WEDGE6)
class FeH1LagrangeWedge1 : public FeH1LagrangeWedge {

public:

  //! Constructor  
  FeH1LagrangeWedge1();

  //! Copy Constructor
  FeH1LagrangeWedge1(const FeH1LagrangeWedge1& other)
   : FeH1LagrangeWedge(other){
  }

  //! Destructor
  virtual  ~FeH1LagrangeWedge1();

  //! Create deep copy
  virtual FeH1LagrangeWedge1* Clone(){
    return new FeH1LagrangeWedge1(*this);
  }
protected:

  //! @copydoc FeH1::CalcShFnc
  void CalcShFnc( Vector<Double>& shape,
                  const Vector<Double>& point,
                  const Elem* ptElem,
                  UInt comp = 1 );

  //! @copydoc FeH1::CalcLocDerivShFnc
  void CalcLocDerivShFnc( Matrix<Double> & deriv, 
                          const Vector<Double>& point,
                          const Elem* ptElem,
                          UInt comp = 1 );
  
  //! @copydoc BaseFE::Triangulate
  virtual void Triangulate(StdVector< StdVector<UInt> > & triConnect);

};

//! Lagrangian tetrahedron element of 1st order (ET_TET4)
class FeH1LagrangeTet1 : public FeH1LagrangeTet {

public:

  //! Constructor
  FeH1LagrangeTet1();

  //! Copy Constructor
  FeH1LagrangeTet1(const FeH1LagrangeTet1& other)
   : FeH1LagrangeTet(other){
  }

  //! Destructor
  virtual  ~FeH1LagrangeTet1();

  //! Create deep copy
  virtual FeH1LagrangeTet1* Clone(){
    return new FeH1LagrangeTet1(*this);
  }

  void ComputeMonomialCoefficients(Matrix<Integer>& P, Matrix<Double>& C);

protected:

  //! @copydoc FeH1::CalcShFnc
  void CalcShFnc( Vector<Double>& shape,
                  const Vector<Double>& point,
                  const Elem* ptElem,
                  UInt comp = 1 );

  //! @copydoc FeH1::CalcLocDerivShFnc
  void CalcLocDerivShFnc( Matrix<Double> & deriv,
                          const Vector<Double>& point,
                          const Elem* ptElem,
                          UInt comp = 1 );

  //! @copydoc BaseFE::Triangulate
  virtual void Triangulate(StdVector< StdVector<UInt> > & triConnect);
};

//! Lagrangian pyramid element of 1st order (ET_PYRA5)
class FeH1LagrangePyra1 : public FeH1LagrangePyra {

public:

  //! Constructor
  FeH1LagrangePyra1();

  //! Copy Constructor
  FeH1LagrangePyra1(const FeH1LagrangePyra1& other)
   : FeH1LagrangePyra(other){
  }

  //! Destructor
  virtual  ~FeH1LagrangePyra1();

  //! Create deep copy
  virtual FeH1LagrangePyra1* Clone(){
    return new FeH1LagrangePyra1(*this);
  }

protected:

  //! @copydoc FeH1::CalcShFnc
  void CalcShFnc( Vector<Double>& shape,
                  const Vector<Double>& point,
                  const Elem* ptElem,
                  UInt comp = 1 );

  //! @copydoc FeH1::CalcLocDerivShFnc
  void CalcLocDerivShFnc( Matrix<Double> & deriv,
                          const Vector<Double>& point,
                          const Elem* ptElem,
                          UInt comp = 1 );

  //! @copydoc BaseFE::Triangulate
  virtual void Triangulate(StdVector< StdVector<UInt> > & triConnect);
};

// ========================================================================
// Lagrange H1 elements of 2nd order
// ========================================================================
//! Lagrangian line element of 2nd order (ET_LINE3)
class FeH1LagrangeLine2 : public FeH1LagrangeLine {

public:

  //! Constructor 
  FeH1LagrangeLine2();

  //! Copy Constructor
  FeH1LagrangeLine2(const FeH1LagrangeLine2& other)
   : FeH1LagrangeLine(other){
  }

  //! Destructor
  virtual ~FeH1LagrangeLine2();

  //! Create deep copy
  virtual FeH1LagrangeLine2* Clone(){
    return new FeH1LagrangeLine2(*this);
  }

protected:

  //! @copydoc FeH1::CalcShFnc
  void CalcShFnc( Vector<Double>& shape,
                  const Vector<Double>& point,
                  const Elem* ptElem,
                  UInt comp = 1 );

  //! @copydoc FeH1::CalcLocDerivShFnc
  void CalcLocDerivShFnc( Matrix<Double> & deriv, 
                          const Vector<Double>& point,
                          const Elem* ptElem,
                          UInt comp = 1 );
};

//! Lagrangian triangle element of 2nd order (ET_TRIA6)
class FeH1LagrangeTria2 : public FeH1LagrangeTria {

public:

  //! Constructor
  FeH1LagrangeTria2();

  //! Copy Constructor
  FeH1LagrangeTria2(const FeH1LagrangeTria2& other)
   : FeH1LagrangeTria(other){
  }

  //! Destructor
  virtual ~FeH1LagrangeTria2();

  //! Create deep copy
  virtual FeH1LagrangeTria2* Clone(){
    return new FeH1LagrangeTria2(*this);
  }

protected:

  //! @copydoc FeH1::CalcShFnc
  void CalcShFnc( Vector<Double>& shape,
                  const Vector<Double>& point,
                  const Elem* ptElem,
                  UInt comp = 1 );

  //! @copydoc FeH1::CalcLocDerivShFnc
  void CalcLocDerivShFnc( Matrix<Double> & deriv,
                          const Vector<Double>& point,
                          const Elem* ptElem,
                          UInt comp = 1 );
};

//! Lagrangian quadrilateral element of 1st order (ET_QUAD8)
class FeH1LagrangeQuad2 : public FeH1LagrangeQuad {

public:

  //! Constructor
  FeH1LagrangeQuad2();

  //! Copy Constructor
  FeH1LagrangeQuad2(const FeH1LagrangeQuad2& other)
   : FeH1LagrangeQuad(other){
  }

  //! Destructor
  virtual ~FeH1LagrangeQuad2();

  //! Create deep copy
  virtual FeH1LagrangeQuad2* Clone(){
    return new FeH1LagrangeQuad2(*this);
  }

protected:

  //! @copydoc FeH1::CalcShFnc
  void CalcShFnc( Vector<Double>& shape,
                  const Vector<Double>& point,
                  const Elem* ptElem,
                  UInt comp = 1 );

  //! @copydoc FeH1::CalcLocDerivShFnc
  void CalcLocDerivShFnc( Matrix<Double> & deriv, 
                          const Vector<Double>& point,
                          const Elem* ptElem,
                          UInt comp = 1 );
};

//! Lagrangian quadrilateral tensor product element of 2nd order (ET_QUAD9)
class FeH1LagrangeQuad9 : public FeH1LagrangeQuad {

public:

  //! Constructor
  FeH1LagrangeQuad9();

  //! Copy Constructor
  FeH1LagrangeQuad9(const FeH1LagrangeQuad9& other)
   : FeH1LagrangeQuad(other){
  }

  //! Destructor
  virtual ~FeH1LagrangeQuad9();

  //! Create deep copy
  virtual FeH1LagrangeQuad9* Clone(){
    return new FeH1LagrangeQuad9(*this);
  }

protected:

  //! @copydoc FeH1::CalcShFnc
  void CalcShFnc( Vector<Double>& shape,
                  const Vector<Double>& point,
                  const Elem* ptElem,
                  UInt comp = 1 );

  //! @copydoc FeH1::CalcLocDerivShFnc
  void CalcLocDerivShFnc( Matrix<Double> & deriv,
                          const Vector<Double>& point,
                          const Elem* ptElem,
                          UInt comp = 1 );
};

//! Lagrangian hexahedral element of 2nd order (ET_HEX20)
class FeH1LagrangeHex2 : public FeH1LagrangeHex {

public:

  //! Constructor  
  FeH1LagrangeHex2();

  //! Copy Constructor
  FeH1LagrangeHex2(const FeH1LagrangeHex2& other)
   : FeH1LagrangeHex(other){
  }

  //! Destructor
  virtual  ~FeH1LagrangeHex2();

  //! Create deep copy
  virtual FeH1LagrangeHex2* Clone(){
    return new FeH1LagrangeHex2(*this);
  }

protected:

  //! @copydoc FeH1::CalcShFnc
  void CalcShFnc( Vector<Double>& shape,
                  const Vector<Double>& point,
                  const Elem* ptElem,
                  UInt comp = 1 );

  //! @copydoc FeH1::CalcLocDerivShFnc
  void CalcLocDerivShFnc( Matrix<Double> & deriv, 
                          const Vector<Double>& point,
                          const Elem* ptElem,
                          UInt comp = 1 );
};

//! Lagrangian hexahedral element of 2nd order (ET_HEX27)
class FeH1LagrangeHex27 : public FeH1LagrangeHex {

public:

  //! Constructor
  FeH1LagrangeHex27();

  //! Copy Constructor
  FeH1LagrangeHex27(const FeH1LagrangeHex27& other)
   : FeH1LagrangeHex(other){
  }

  //! Destructor
  virtual  ~FeH1LagrangeHex27();

  //! Create deep copy
  virtual FeH1LagrangeHex27* Clone(){
    return new FeH1LagrangeHex27(*this);
  }

protected:

  //! @copydoc FeH1::CalcShFnc
  void CalcShFnc( Vector<Double>& shape,
                  const Vector<Double>& point,
                  const Elem* ptElem,
                  UInt comp = 1 );

  //! @copydoc FeH1::CalcLocDerivShFnc
  void CalcLocDerivShFnc( Matrix<Double> & deriv,
                          const Vector<Double>& point,
                          const Elem* ptElem,
                          UInt comp = 1 );
};

//! Lagrangian wedge element of 2nd order (ET_WEDGE15)
class FeH1LagrangeWedge2 : public FeH1LagrangeWedge {

public:

  //! Constructor  
  FeH1LagrangeWedge2();

  //! Copy Constructor
  FeH1LagrangeWedge2(const FeH1LagrangeWedge2& other)
   : FeH1LagrangeWedge(other){
  }

  //! Destructor
  virtual  ~FeH1LagrangeWedge2();

  //! Create deep copy
  virtual FeH1LagrangeWedge2* Clone(){
    return new FeH1LagrangeWedge2(*this);
  }

protected:

  //! @copydoc FeH1::CalcShFnc
  void CalcShFnc( Vector<Double>& shape,
                  const Vector<Double>& point,
                  const Elem* ptElem,
                  UInt comp = 1 );

  //! @copydoc FeH1::CalcLocDerivShFnc
  void CalcLocDerivShFnc( Matrix<Double> & deriv, 
                          const Vector<Double>& point,
                          const Elem* ptElem,
                          UInt comp = 1 );
};

//! Lagrangian wedge element of 2nd order (ET_WEDGE18)
class FeH1LagrangeWedge18 : public FeH1LagrangeWedge {

public:

  //! Constructor
  FeH1LagrangeWedge18();

  //! Copy Constructor
  FeH1LagrangeWedge18(const FeH1LagrangeWedge18& other)
   : FeH1LagrangeWedge(other){
  }

  //! Destructor
  virtual  ~FeH1LagrangeWedge18();

  //! Create deep copy
  virtual FeH1LagrangeWedge18* Clone(){
    return new FeH1LagrangeWedge18(*this);
  }

protected:

  //! @copydoc FeH1::CalcShFnc
  void CalcShFnc( Vector<Double>& shape,
                  const Vector<Double>& point,
                  const Elem* ptElem,
                  UInt comp = 1 );

  //! @copydoc FeH1::CalcLocDerivShFnc
  void CalcLocDerivShFnc( Matrix<Double> & deriv,
                          const Vector<Double>& point,
                          const Elem* ptElem,
                          UInt comp = 1 );

};
//! Lagrangian tetrahedron element of 2nd order (ET_TET10)
class FeH1LagrangeTet2 : public FeH1LagrangeTet {

public:

  //! Constructor
  FeH1LagrangeTet2();

  //! Copy Constructor
  FeH1LagrangeTet2(const FeH1LagrangeTet2& other)
   : FeH1LagrangeTet(other){
  }

  //! Destructor
  virtual  ~FeH1LagrangeTet2();

  //! Create deep copy
  virtual FeH1LagrangeTet2* Clone(){
    return new FeH1LagrangeTet2(*this);
  }

protected:

  //! @copydoc FeH1::CalcShFnc
  void CalcShFnc( Vector<Double>& shape,
                  const Vector<Double>& point,
                  const Elem* ptElem,
                  UInt comp = 1 );

  //! @copydoc FeH1::CalcLocDerivShFnc
  void CalcLocDerivShFnc( Matrix<Double> & deriv,
                          const Vector<Double>& point,
                          const Elem* ptElem,
                          UInt comp = 1 );
};

//! Lagrangian pyramid element of 2nd order (ET_PYRA13)
class FeH1LagrangePyra2 : public FeH1LagrangePyra {

public:

  //! Constructor
  FeH1LagrangePyra2();

  //! Copy Constructor
  FeH1LagrangePyra2(const FeH1LagrangePyra2& other)
   : FeH1LagrangePyra(other){
  }

  //! Destructor
  virtual  ~FeH1LagrangePyra2();

  //! Create deep copy
  virtual FeH1LagrangePyra2* Clone(){
    return new FeH1LagrangePyra2(*this);
  }

protected:

  //! @copydoc FeH1::CalcShFnc
  void CalcShFnc( Vector<Double>& shape,
                  const Vector<Double>& point,
                  const Elem* ptElem,
                  UInt comp = 1 );

  //! @copydoc FeH1::CalcLocDerivShFnc
  void CalcLocDerivShFnc( Matrix<Double> & deriv,
                          const Vector<Double>& point,
                          const Elem* ptElem,
                          UInt comp = 1 );
};

//! Lagrangian pyramid element of 2nd order (ET_PYRA14)
class FeH1LagrangePyra14 : public FeH1LagrangePyra {

public:

  //! Constructor
  FeH1LagrangePyra14();

  //! Copy Constructor
  FeH1LagrangePyra14(const FeH1LagrangePyra14& other)
   : FeH1LagrangePyra(other){
  }

  //! Destructor
  virtual  ~FeH1LagrangePyra14();

  //! Create deep copy
  virtual FeH1LagrangePyra14* Clone(){
    return new FeH1LagrangePyra14(*this);
  }

protected:

  //! @copydoc FeH1::CalcShFnc
  void CalcShFnc( Vector<Double>& shape,
                  const Vector<Double>& point,
                  const Elem* ptElem,
                  UInt comp = 1 );

  //! @copydoc FeH1::CalcLocDerivShFnc
  void CalcLocDerivShFnc( Matrix<Double> & deriv,
                          const Vector<Double>& point,
                          const Elem* ptElem,
                          UInt comp = 1 );
};

} // namespace CoupledField

#endif
