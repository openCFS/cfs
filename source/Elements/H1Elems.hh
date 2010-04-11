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

    //! Evaluate Lagrange polynomial 
    void EvaluateLagrangePolynomial( Vector<Double> & shape, Double coord );

    //! Evaluate derivative of Lagrange polynomials 
    void EvaluateDerivLagrangePolynomial( Vector<Double> & deriv, Double coord );

    //! Return global derivative of shape functions
    void GetGlobDerivShFnc( Matrix<Double>& deriv, LocPointMapped& lp,
                            const Elem* elem, UInt comp = 1 );

    //! Set the isotropic order of the Element. This methods gets overwritten 
    //! by the child classes to calculate the number of functions according to
    //! the given order
    //! \param order (input) The desired order of the element
    virtual void SetIsoOrder(UInt order){
      // after debugging phase, this EXCEPTION can become a warning
      EXCEPTION("Trying to set the ISOTROPIC order of an element which does not support this opeeration.\
                The element will remain unchanged!");
    }

    //! Set the Anisotropic order of the Element. This methods gets overwritten 
    //! by the child classes to calculate the number of functions according to
    //! the given order
    //! \param order (input) vector of element orders for each space direction 
    virtual void SetAnisoOrder(StdVector<UInt> order){
      // after debugging phase, this EXCEPTION can become a warning
      EXCEPTION("Trying to set the ANISOTROPIC order of an element which does not support this operation.\
                The element will remain unchanged!");
    }

  protected:

    //! Calculate the location of unknowns for a line up to given order
    //! Righ now the calculation is done here but has to move into the 
    //! Integration Scheme class!
    //! This method gets reimplemented in the Spectral classes to ask the 
    //! Integration scheme to provide the points
    void CalcAllSupportingPoints(UInt maxOrder);

    //! Calculate the location of unknowns for a given order
    void CalcSupportingPoints(UInt order);

    //! Calculates the Gauss Lobatto Points for a given order
    //! Is to be moved into Integration Scheme Class
    StdVector<Double> CalcGaussLobattoPoints(UInt order);

    //! Stores the Locations of the Element DOFs for a line for every order
    std::map<UInt,StdVector<Double> > supportingPoints_;

  private:
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
    virtual void GetFncPermutation( StdVector<UInt>& fncPermutation,
                                    const Elem* ptElem,
                                    EntityType fctEntityType,
                                    UInt entNumber);

    ////! Return shape function
    //void GetShFnc( Vector<Double> & S, const LocPoint& lp,
    //               const Elem* ptElem,  UInt comp = 1 );

    ////! Return local derivative of shape function
    //void GetDerivShFnc( Matrix<Double> & deriv, 
    //                    const LocPoint& lp,
    //                    const Elem * elem, 
    //                    UInt comp = 1 );

    //! returns the number of functions for a single edge or face
    UInt GetNumFncsPerEntType( EntityType fctEntityType, UInt dof = 1);

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

  
} // namespace CoupledField

#endif
