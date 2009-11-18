// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
#ifndef FILE_CFS_H1_ELEMENTS_HH
#define FILE_CFS_H1_ELEMENTS_HH

#include "basefe.hh"

namespace CoupledField {

  //! Base class for H1-conforming reference elements
  
  //! The H1 element represents all elements with H1 conforming shape function
  class FeH1 : public BaseFE {
  public:

    //! Constructor
    FeH1();

    //! Destructor
    virtual ~FeH1();

    //! Evaluate Lagrange polynomial (to be implemented by
    //! Andi Hueppe) 
    void EvaluatePolynomial( Vector<Double> & shape, Double coord );

    //! Evaluate derivative of Lagrange polynomials (to be implemented by
    //! Andi Hueppe)
    void EvaluateDerivPolynomial( Vector<Double> & deriv, Double coord );

    //! Return global derivative of shape functions
    void GetGlobDerivShFnc( Matrix<Double>& deriv, LocPointMapped& lp,
                            const Elem* elem, UInt comp = 1 );

  protected:

    //! Node numbers (to determine orientation)
    StdVector<UInt> nodeNums;
  };

  
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
                     const shared_ptr<AnsatzFct>& fcnType,
                     AnsatzFct::FctEntityType fctEntityType,
                     UInt dof = 1 );

    //! Return shape function
    void GetShFnc( Vector<Double> & S, const LocPoint& lp,
                   const Elem* ptElem,  UInt comp = 1 );

    //! Return local derivative of shape function
    void GetDerivShFnc( Matrix<Double> & deriv, 
                        const LocPoint& lp,
                        const Elem * elem, 
                        UInt comp = 1 );
  protected:

    //! Compute shape function at given position
    virtual void CalcShFnc( Vector<Double>& shape,
                            const Vector<Double>& point ) = 0;

    //! Compute local derivative at of shape function at given position
    virtual void CalcDerivShFnc( Matrix<Double> & deriv, 
                                 const Vector<Double>& point ) = 0;

    //! Valued of shape functions at integration points
    StdVector<Vector<Double> > shapeAtIp_;

    //! Value of local derivatives of shape functions at integration points
    StdVector<Matrix<Double> > shapeDerivAtIp_;
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

    //! @see FeH1LagrangeExpl::CalcShapeFnc
    void CalcShFnc( Vector<Double>& shape,
                    const Vector<Double>& point );

    //! @see FeH1LagrangeExpl::CalcDerivShFnc
    void CalcDerivShFnc( Matrix<Double> & deriv, 
                         const Vector<Double>& point );
  };

  //! Lagrangian quadrilateral element of 1st order (ET_QUAD4)
  class FeH1LagrangeQuad1 : public FeH1LagrangeExpl {

  public:

    //! Constructor
    FeH1LagrangeQuad1();

    //! Destructor
    virtual ~FeH1LagrangeQuad1();

  protected:

    //! @see FeH1LagrangeExpl::CalcShapeFnc 
    void CalcShFnc( Vector<Double>& shape,
                    const Vector<Double>& point );

    //! @see FeH1LagrangeExpl::CalcDerivShFnc 
    void CalcDerivShFnc( Matrix<Double> & deriv, 
                         const Vector<Double>& point );
  };

  //! Lagrangian hexahedral element of 1st order (ET_HEX8)
  class FeH1LagrangeHex1 : public FeH1LagrangeExpl {

  public:

    //! Constructor  
    FeH1LagrangeHex1();

    //! Destructor
    virtual  ~FeH1LagrangeHex1();

  protected:

    //! @see FeH1LagrangeExpl::CalcShapeFnc 
    void CalcShFnc( Vector<Double>& shape,
                    const Vector<Double>& point );

    //! @see FeH1LagrangeExpl::CalcDerivShFnc 
    void CalcDerivShFnc( Matrix<Double> & deriv, 
                         const Vector<Double>& point );
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

     //! @see FeH1LagrangeExpl::CalcShapeFnc
     void CalcShFnc( Vector<Double>& shape,
                     const Vector<Double>& point );

     //! @see FeH1LagrangeExpl::CalcDerivShFnc
     void CalcDerivShFnc( Matrix<Double> & deriv, 
                          const Vector<Double>& point );
   };

   //! Lagrangian quadrilateral element of 1st order (ET_QUAD8)
   class FeH1LagrangeQuad2 : public FeH1LagrangeExpl {

   public:

     //! Constructor
     FeH1LagrangeQuad2();

     //! Destructor
     virtual ~FeH1LagrangeQuad2();

   protected:

     //! @see FeH1LagrangeExpl::CalcShapeFnc 
     void CalcShFnc( Vector<Double>& shape,
                     const Vector<Double>& point );

     //! @see FeH1LagrangeExpl::CalcDerivShFnc 
     void CalcDerivShFnc( Matrix<Double> & deriv, 
                          const Vector<Double>& point );
   };

   //! Lagrangian hexahedral element of 2nd order (ET_HEX20)
   class FeH1LagrangeHex2 : public FeH1LagrangeExpl {

   public:

     //! Constructor  
     FeH1LagrangeHex2();

     //! Destructor
     virtual  ~FeH1LagrangeHex2();

   protected:

     //! @see FeH1LagrangeExpl::CalcShapeFnc 
     void CalcShFnc( Vector<Double>& shape,
                     const Vector<Double>& point );

     //! @see FeH1LagrangeExpl::CalcDerivShFnc 
     void CalcDerivShFnc( Matrix<Double> & deriv, 
                          const Vector<Double>& point );
   };


  // ========================================================================
  // Lagrange H1 elements of arbitrary order (evaluation using
  // general Lagrange polynomials)
  // ========================================================================

   //! Lagrangian line element of variable order
   class FeH1LoLineVar : public FeH1 {

   };
   //! Lagrangian quadrilateral element of variable order 
   class FeH1LoQuadVar : public FeH1 {

   };
   //! Lagrangian hexahedral element of varaiable order
   class FeH1LoHexVar : public FeH1 {

   };  
  
} // namespace CoupledField

#endif
