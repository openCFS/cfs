// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
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
  
  //! =======================================================================
  //! G E O M E T R I C   I N F O R M A T I O N
  //! =======================================================================
  //! Returns whether a given local coordinate is within this element
  //! @param point input Local point
  //! @param tolerance input Additioanl (relative) tolerance
  //! \return flag if point is inside the element
  virtual bool CoordIsInsideElem( const Vector<Double>& point,
                                  Double tolerance ) = 0;

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

};


} // namespace CoupledField

#endif
