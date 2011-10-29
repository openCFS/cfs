// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode c
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
    void EvaluateLagrangePolynomial( Vector<Double> & shape, Double coord,
                                     UInt order);

    //! Evaluate derivative of Lagrange polynomials 
    void EvaluateDerivLagrangePolynomial( Vector<Double> & deriv, Double coord,
                                          UInt order);
    
    //! Get value of all shape fnc at local poiont lp
    /*!
    \param S (output) Vector of shape fnc values \f$ (N_{1},\cdots\,N_{NumNodes})^T \f$
    \param ip (input) Integration point
    */
    virtual void GetShFnc( Vector<Double>& S, const LocPoint& lp,
                           const Elem* ptElem,  UInt comp = 1 );
    
    //! Return global derivative of shape functions
    void GetGlobDerivShFnc( Matrix<Double>& deriv, const LocPointMapped& lp,
                            const Elem* elem, UInt comp = 1 );

    //! Return global derivative of shape functions
    void GetLocDerivShFnc( Matrix<Double>& deriv, const LocPoint& lp,
                           const Elem* elem, UInt comp = 1 );

    //! @copydoc BaseFE::SetFunctionsAtIp
    void SetFunctionsAtIp(const StdVector<LocPoint>& iPoints);
        


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
    
    //! Compute shape function at given position
    virtual void CalcShFnc( Vector<Double>& shape,
                            const Vector<Double>& point,
                            const Elem* ptElem,
                            UInt comp = 1 ) = 0;

    //! Get local derivatives of all shape fnc at arbitrary local point
    //! Local means here on the reference element
    /*! 
    \param S (output) Matrix with global derivatives of all shape functions
    \f [ \left( \begin{array}{ccc} N_{1,dx} & N_{1,dy} & \cdots \\
    N_{2,dx} & N_{2,dy} & \cdots \\
    \cdots     & \cdots      & \cdots \end{array}\right) \f ]
    \param LCoord (input) Local Coordinates of evalutaion point
    \param CornerCoords (input) Coordinates of element corners
    \f [ \left( \begin{array}{ccc} x_{1} & x_{2} & \cdots \\ y_{1} & y_{2} & \cdots \\
    \cdots & \cdots & \cdots \end{array} \right) \f ]       
    */
    virtual void CalcLocDerivShFnc( Matrix<Double> & deriv, 
                                    const Vector<Double>& point,
                                    const Elem* ptElem,
                                    UInt comp = 1 ) = 0;
    
    // =======================================================================
    //  PRE CALCULATION OF SHAPE FUNCTIONS AT INTEGRATION POINTS
    // =======================================================================

    //! Stores Shape Functions for each integration point definied
    StdVector< Vector<Double> > shapeFncsAtIp_;

    //! Stores shape function derivatives for each integration point
    StdVector< Matrix<Double> > shapeFncDerivsAtIp_;
    
    //! Stores the Locations of the Element DOFs for a line for every order
    std::map<UInt,StdVector<Double> > supportingPoints_;


  private:
  };

} // namespace CoupledField

#endif
