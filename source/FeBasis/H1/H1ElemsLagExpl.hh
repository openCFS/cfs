#ifndef FILE_CFS_H1_ELEMENTS_LAGRANGE_EXPLICIT_HH
#define FILE_CFS_H1_ELEMENTS_LAGRANGE_EXPLICIT_HH

#include "H1Elems.hh"

namespace CoupledField {

// ========================================================================
//  H1 Fe Lagrangian Elements of lowest order (1st / 2nd)
// ========================================================================

//! Lagrangian based elements of lowest order (explicitly coded)

//! This serves as a base class for all Lagrangian bases Elements of 
//! lowest order, i.e. the  shape functions can be pre-calculated without 
//! knowledge of the global element.
//! Currently these elements are all 1st / 2nd order Lagrangian elements
class FeH1LagrangeExpl : public FeH1 {

public:
  //! Constructor
  FeH1LagrangeExpl();

  //! Destructor
  virtual ~FeH1LagrangeExpl();

  //! Pre-calculate values at integration points
  void SetIntPoints( StdVector<LocPoint>& intPoints );

  //! Get number of shape functions for a given type (NODE/EDGE/FACE/ELEM)
  void GetNumFncs( StdVector<UInt>& numFcns,
                   EntityType fctEntityType,
                   UInt dof = 1 );

  //! Get the permutation Vector for a given Face or Edge
  //! e.g. If asked for a face, the element will check the flags
  //! of this face and return a vector of size NumberOfFncs on the Face
  //! holding the correct ordering 
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

  //! returns the number of functions for a single edge or face
  UInt GetNumFncsPerEntType( EntityType fctEntityType, UInt dof = 1);

  //! \copydoc BaseFE::GetIsoOrder
  virtual UInt GetIsoOrder() const {
    return order_;
  }
  //! \copydoc BaseFE::GetMaxOrder
  virtual UInt GetMaxOrder() const {
    return order_;
  }
  
  //! \copydoc BaseFE::GetMaxOrderLocDir
  virtual void GetMaxOrderLocDir(StdVector<UInt>& order );
  
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
  /*!
       \param surfConnect (input) Node numbers of surface element
       \param volConnect (input) Node numbers of volume element
       \param surfIntPoint (input) Surface integration point, which gets mapped
       onto the volume element
       \param volIntPoint (output) Corresponding volume integration point
       \param locNormal (output) Normal direction of surface element in local 
       coordinates of the volume element
   */
  virtual void GetLocalIntPoints4Surface(const StdVector<UInt> & surfConnect,
                                         const StdVector<UInt> & volConnect,
                                         const LocPoint & surfIntPoint,
                                         LocPoint & volIntPoint,
                                         Vector<Double>& locNormal ) = 0;

  //! @copydoc BaseFE::ComputeMonomialCoefficients
  //! Overloaded method for lagrange Elements
  //virtual void ComputeMonomialCoefficients(Matrix<Integer>& P, Matrix<Double>& C);

protected:

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
  
  //! Polynomial order of the finite element
  UInt order_;
};


// ========================================================================
// Lagrange H1 elements of 1st order
// ========================================================================

//! Lagrangian line element of 1st order (ET_LINE2)
class FeH1LagrangeLine1 : public FeH1LagrangeExpl {

public:

  //! Constructor 
  FeH1LagrangeLine1();

  //! Destructor
  virtual ~FeH1LagrangeLine1();

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

//! Lagrangian triangular element of 1st order (ET_TRIA3)
class FeH1LagrangeTria1 : public FeH1LagrangeExpl {

public:

  //! Constructor
  FeH1LagrangeTria1();

  //! Destructor
  virtual ~FeH1LagrangeTria1();

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

//! Lagrangian quadrilateral element of 1st order (ET_QUAD4)
class FeH1LagrangeQuad1 : public FeH1LagrangeExpl {

public:

  //! Constructor
  FeH1LagrangeQuad1();

  //! Destructor
  virtual ~FeH1LagrangeQuad1();

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

//! Lagrangian hexahedral element of 1st order (ET_HEX8)
class FeH1LagrangeHex1 : public FeH1LagrangeExpl {

public:

  //! Constructor  
  FeH1LagrangeHex1();

  //! Destructor
  virtual  ~FeH1LagrangeHex1();

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


//! Lagrangian wedge element of 1st order (ET_WEDGE6)
class FeH1LagrangeWedge1 : public FeH1LagrangeExpl {

public:

  //! Constructor  
  FeH1LagrangeWedge1();

  //! Destructor
  virtual  ~FeH1LagrangeWedge1();

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
// Lagrange H1 elements of 2nd order
// ========================================================================
//! Lagrangian line element of 2nd order (ET_LINE3)
class FeH1LagrangeLine2 : public FeH1LagrangeExpl {

public:

  //! Constructor 
  FeH1LagrangeLine2();

  //! Destructor
  virtual ~FeH1LagrangeLine2();

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

//! Lagrangian quadrilateral element of 1st order (ET_QUAD8)
class FeH1LagrangeQuad2 : public FeH1LagrangeExpl {

public:

  //! Constructor
  FeH1LagrangeQuad2();

  //! Destructor
  virtual ~FeH1LagrangeQuad2();

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

//! Lagrangian hexahedral element of 2nd order (ET_HEX20)
class FeH1LagrangeHex2 : public FeH1LagrangeExpl {

public:

  //! Constructor  
  FeH1LagrangeHex2();

  //! Destructor
  virtual  ~FeH1LagrangeHex2();

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

//! Lagrangian wedge element of 1st order (ET_WEDGE15)
class FeH1LagrangeWedge2 : public FeH1LagrangeExpl {

public:

  //! Constructor  
  FeH1LagrangeWedge2();

  //! Destructor
  virtual  ~FeH1LagrangeWedge2();

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



} // namespace CoupledField

#endif
