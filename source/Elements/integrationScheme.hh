/*
 * inttable.hh
 *
 *  Created on: 23.07.2009
 *      Author: ahauck
 */

#ifndef FILE_INTEGRATION_SCHEME_HH_
#define FILE_INTEGRATION_SCHEME_HH_

#include "Domain/elem.hh"
#include "Elements/elemShapeMap.hh"

namespace CoupledField {

  //! Class defininig Numerical Integration

  //! The integration scheme returns the integration points / elements for a 
  //! given element shape. As additional parameters one cann pass the 
  //! order (maybe anisotropic) and the type of integration used
  class IntScheme {
    
  public:

    /** enumeration with integration types - to be renamed into IntegrationType
    * UNDEFINED: only interal use / not define
    * GAUSS: Classic Gauss-Legendre integration points (quad, hex: tensor product of 1D points)
    * GAUSS_CART: As GAUSS, only with anisotropic weightning in 3 directions. 
    *             The order is encoded in xxyyzz form x1 =0-9, x2 = 1x-9x, x3 = 1xx-9xx
    * GAUSS_ECO: are the "efficient" gaussian quadrature weights
    *            (->Solin, Segeth, Dolezel, Higher-Order Finite Element Methods)
    * LOBATTO: Gauss-Lobatto integration points (used for spectral method)
    * CHEBYSHEV: special line intgration methods (->Solin, Segeth, Dolezel, 
    *            Higher-Order Finite Element Methods)
    * SPECIAL: Special points (mainly used for HEX27 element)*/
    typedef enum { UNDEFINED, GAUSS, GAUSS_CART, GAUSS_ECO, 
                   LOBATTO, CHEBYSHEV, SPECIAL } IntegMethod;
    static Enum<IntegMethod> IntegMethodEnum;
    
    typedef std::map<Elem::ShapeType,StdVector<LocPoint> > IntegrationPoints;
    typedef std::map<Elem::ShapeType,StdVector<Double> > IntegrationWeights;
    
    //! Constructor
    IntScheme();

    //! Destructor
    ~IntScheme();
  
    //! Set method and order
    void SetOrder( std::string method, UInt order);

    //! Set method and order
    void SetOrder( IntegMethod method, UInt order);
    
    //! Obtain method
    IntegMethod GetMethod() { return integMethod_; }
    
    //! obtain order
    UInt GetOrder() { return order_;}

    //! Get integration points / weights
    void GetIntPoints( Elem::ShapeType elemType,
                       StdVector<LocPoint>& intPts, 
                       StdVector<Double>& weights );  

    //! Returns all definied integration points in a single vector
    void GetAllIntegrationPoints(StdVector< LocPoint >& points, Elem::ShapeType type); 

private:
    //! Adds the Gauss Lobatto/Legendre  Points up to the given order to the Integration maps
    void FillIntegPoints(UInt order);

    //! Calculate the Gauss-Lobatto points an weights for the given order
    void CalcGaussLobattoPointsWeights(UInt order,StdVector<Double>& points, StdVector<Double>& weights);
    
    //! Calculate the Gauss-Lobatto points an weights for the given order
    void CalcGaussLegendrePointsWeights(UInt order, StdVector<Double>& points, 
                                        StdVector<Double>& weights );
  
    //! Map with integration points for each element type According to the Template Paramter Integration Scheme
    std::map<IntegMethod, std::map< UInt, IntegrationPoints > > intPoints_;
    
    //! Map with integration weights for each element type
    std::map<IntegMethod, std::map< UInt, IntegrationWeights > > intWeights_;

    //! the current method of integration
    IntegMethod integMethod_;

    UInt order_;

    //! stores the overall number of integration points definied in this class
    //! for each element type, zero based
    std::map<Elem::ShapeType,UInt> numIntPts_;
      
};



} // namespace CoupledField

#endif /* INTTABLE_HH_ */
