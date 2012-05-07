#ifndef FILE_CFS_H1_ELEMENTS_HH
#define FILE_CFS_H1_ELEMENTS_HH

#include "FeBasis/BaseFE.hh"


namespace CoupledField {

  //! Base class for H1-conforming reference elements
  
  //! The H1 element represents all elements with H1 conforming shape function
  class FeH1 : public BaseFE {
  public:

    //! Constructor
    FeH1();

    //! Destructor
    virtual ~FeH1();

    //! Get value of all shape fnc at local poiont lp
    /*!
    \param S (output) Vector of shape fnc values \f$ (N_{1},\cdots\,N_{NumNodes})^T \f$
    \param ip (input) Integration point
    */
    virtual void GetShFnc( Vector<Double>& S, const LocPoint& lp,
                           const Elem* ptElem,  UInt comp = 1 );
    
    //! Return global derivative of shape functions
    void GetGlobDerivShFnc( Matrix<Double>& deriv, 
                            const LocPointMapped& lp,
                            const Elem* elem, UInt comp = 1 );

    //! Return global derivative of shape functions
    void GetLocDerivShFnc( Matrix<Double>& deriv, 
                           const LocPoint& lp,
                           const Elem* elem, UInt comp = 1 );

  protected:

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
    //  Changed here to map data structure as it is more versatile and
    //  seeing the fact, that we do not have millions of integration points,
    //  it should be ok
    // =======================================================================

    //! Stores Shape Functions for each integration point definied
    std::map<Integer, Vector<Double> > shapeFncsAtIp_;

    //! Stores shape function derivatives for each integration point
    std::map<Integer, Matrix<Double> > shapeFncDerivsAtIp_;

  private:
  };

} // namespace CoupledField

#endif
