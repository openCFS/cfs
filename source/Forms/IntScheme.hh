/*
 * inttable.hh
 *
 *  Created on: 23.07.2009
 *      Author: ahauck
 */

#ifndef FILE_INTEGRATION_SCHEME_HH_
#define FILE_INTEGRATION_SCHEME_HH_

#include <boost/functional/hash.hpp>
#include <boost/array.hpp>
#include <boost/unordered_map.hpp>
// needed for defining a hash function for integration orders


#include "Domain/ElemMapping/Elem.hh"
#include "Domain/ElemMapping/ElemShapeMap.hh"



namespace CoupledField {


  //! Lightweight struct for defining the integration order
  
  //! This struct encapsulates the integration order of an element.
  //! The integration order can be defined either isotropic by just
  //! an integer value or anisotropic, allowing a different order
  //! in each element-local direction (xi, eta, zeta).
  class IntegOrder {
    
  public:
    
    
    // ----------------------------------------------------------------------
    //  Initialization 
    // ----------------------------------------------------------------------
    //! Default constructor
    IntegOrder();
    
    //! Constructor for isotropic order
    IntegOrder( UInt order);
    
    //! Constructor for anisotropic order
    IntegOrder( const StdVector<UInt>& order );
    
    // ----------------------------------------------------------------------
    //  Set Methods
    // ----------------------------------------------------------------------

    //! Set isotropic order
    void SetIsoOrder( UInt order );

    //! Set anisotropic order
    void SetAnisoOrder( const StdVector<UInt>& order );

    //! Adds an isotropic offset to the integration order
    IntegOrder& operator+=( const UInt add );
    
    //! Adds a second IntegOrder struct to itself
    IntegOrder& operator+=( const IntegOrder& other );
    
    // ----------------------------------------------------------------------
    //  Query Methods
    // ----------------------------------------------------------------------
    
    //! Query if integration is set at all
    bool IsSet() const;
    
    //! Query if integration order is isotropic
    bool IsIsotropic() const;
    
    //! Return isotropic order. Returns 0 if not set
    UInt GetIsoOrder() const;
    
    //! Return anisotropic order. Returns empty vector if not set.
    void GetAnisoOrder( StdVector<UInt>& order ) const;
    
    //! Convert to string representation
    std::string ToString() const;
    
    //! Returns the maximum of two integration order structs
    static IntegOrder GetMax( const IntegOrder& order1, 
                              const IntegOrder& order2 );
    
    // ----------------------------------------------------------------------
    //  Hash method
    // ----------------------------------------------------------------------
    
    //! Test for equality
    bool operator== (const IntegOrder &) const;
    
    //! Helper struct for hash function
    struct Hash
    : std::unary_function<IntegOrder, std::size_t> {
      std::size_t operator()(IntegOrder const& p) const {
        return boost::hash_range(&(p.order_[0]), &(p.order_[2]));
      }
    };

    
  private:

    //! Vector containing the integration order
    boost::array<UInt,3> order_;
    
    //! Flag if integration order is isotropic
    bool isIsotropic_;
    
    //! Flag if order is set
    bool isSet_;
    
  };


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

    //! Get integration points / weights
    void GetIntPoints( Elem::ShapeType elemType,
                       IntegMethod method,
                       const IntegOrder& order,
                       StdVector<LocPoint>& intPts, 
                       StdVector<Double>& weights );
    
    //! Get integration points / weights
    void GetIntPoints( Elem::ShapeType elemType,
                       IntegMethod method1,
                       const IntegOrder& order1,
                       IntegMethod method2,
                       const IntegOrder& order2,
                       StdVector<LocPoint>& intPts, 
                       StdVector<Double>& weights );  

    //! Returns all definied integration points in a single vector
    void GetAllIntegrationPoints(StdVector< LocPoint >& points, Elem::ShapeType type); 

    //! Returns definied integration points in a single vector for a given order
    void GetIntegrationPoints(std::map<Integer, LocPoint >& points, Elem::ShapeType type,
                              IntegMethod method, const IntegOrder& order );

private:
    //! Adds the Gauss Lobatto/Legendre  Points up to the given order to the Integration maps
    void FillIntegPoints(UInt order);

    //! Calculate the Gauss-Lobatto points an weights for the given order
    void CalcGaussLobattoPointsWeights( UInt order,StdVector<Double>& points, 
                                        StdVector<Double>& weights );
    
    //! Calculate the Gauss-Lobatto points an weights for the given order
    void CalcGaussLegendrePointsWeights( UInt order, StdVector<Double>& points, 
                                         StdVector<Double>& weights );

    //! Generate set of integration points / weight for gien element, method and order 
    
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

    //! typdef for map of integration points
    typedef boost::unordered_map< IntegOrder, IntegrationPoints, IntegOrder::Hash > IntPointMap;
    
   //! typedef for map of integration weights
    typedef boost::unordered_map< IntegOrder, IntegrationWeights, IntegOrder::Hash > IntWeightMap;
    
    //! Map with integration points for each element type According to the Template Parameter Integration Scheme
    std::map<IntegMethod, IntPointMap > intPoints_;
    
    //! Map with integration weights for each element type
    std::map<IntegMethod, IntWeightMap > intWeights_;

    //! stores the overall number of integration points definied in this class
    //! for each element type, zero based
    std::map<Elem::ShapeType,UInt> numIntPts_;
      
};



} // namespace CoupledField

#endif /* INTTABLE_HH_ */
