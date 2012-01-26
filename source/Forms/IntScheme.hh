/*
 * inttable.hh
 *
 *  Created on: 23.07.2009
 *      Author: ahauck
 */

#ifndef FILE_INTEGRATION_SCHEME_HH_
#define FILE_INTEGRATION_SCHEME_HH_

#include "Domain/ElemMapping/Elem.hh"
#include "Domain/ElemMapping/ElemShapeMap.hh"

namespace CoupledField {

  //! Class defining Numerical Integration

  //! The integration scheme returns the integration points / elements for a 
  //! given element shape. As additional parameters one can pass the 
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
    //! \param order Integration order (-1 defaults to 2) 
    void SetOrder( std::string method, Integer order);

    //! Set method and order
    void SetOrder( IntegMethod method, Integer order);
    
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

    //! Returns definied integration points in a single vector for a given order
    void GetIntegrationPoints(std::map<Integer, LocPoint >& points, Elem::ShapeType type,UInt order,IntegMethod method);

private:
    //! Adds the Gauss Lobatto/Legendre  Points up to the given order to the Integration maps
    void FillIntegPoints(UInt order);

    //! Calculate the Gauss-Lobatto points an weights for the given order
    void CalcGaussLobattoPointsWeights( UInt order,StdVector<Double>& points, 
                                        StdVector<Double>& weights );
    
    //! Calculate the Gauss-Lobatto points an weights for the given order
    void CalcGaussLegendrePointsWeights( UInt order, StdVector<Double>& points, 
                                         StdVector<Double>& weights );

    //! Add single integration point set
    
    //! This method allows to define by hand integrations points for a 
    //! given element.
    //! \param method Integration method  
    //! \param shape Shape of element
    //! \param order Integration order
    //! \param numPoints Number of points
    //! \param data Actually [][2] for 1D (coord + weight) and [][3] 
    //!             for 2D and  [][4] for 3D
    void AddIntegrationSet( IntegMethod method, Elem::ShapeType shape, 
                            UInt order, UInt numPoints, Double* data );
    
    //! Define integrations points / weights for triangular elements
    void DefineTriaPoints();
    
    //! Calculate integration points for triangles with Duffy transformation
    void CalcIntTria( IntegMethod, UInt order,StdVector<LocPoint>& points,
                      StdVector<Double>& weights, bool numberPoints );
    
   //! Define integration points for wedge
   void DefineWedgePoints();

   //! Define integration points for tetrahedron
   void DefineTetPoints();

   //! Define integration points for pyramid
   void DefinePyraPoints();

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
