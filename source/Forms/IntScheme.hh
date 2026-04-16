/*
 * inttable.hh
 *
 *  Created on: 23.07.2009
 *      Author: ahauck
 */

#ifndef FILE_INTEGRATION_SCHEME_HH_
#define FILE_INTEGRATION_SCHEME_HH_

#include <array>
#include <unordered_map>
#include <functional> // for std::hash

#include "Domain/ElemMapping/Elem.hh"
#include "Domain/ElemMapping/ElemShapeMap.hh"



namespace CoupledField {


  //! Lightweight struct for defining the integration order
  
  //! This struct encapsulates the integration order of an element.
  //! The integration order can be defined either isotropic by just
  //! an integer value or anisotropic, allowing a different order
  //! in each element-local direction (xi, eta, zeta).
  class IntegOrder {
    
    friend bool operator< (const IntegOrder&, const IntegOrder&);
    friend std::size_t hash_value(const IntegOrder& p);
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
    
    //! Multiplication of integration order by scalar factor
    IntegOrder operator*(const UInt mult) const;
    
    //! Add scalar, isotropic offset to integration order 
    IntegOrder operator+(const UInt add) const;
    
    //! Subtract scalar, isotropic offset to integration order 
    IntegOrder operator-(const UInt add) const;

    // ----------------------------------------------------------------------
    //  Query Methods
    // ----------------------------------------------------------------------
    
    //! Query if integration is set at all
    bool IsSet() const;
    
    //! Query if integration order is isotropic
    bool IsIsotropic() const;
    
    //! Return isotropic order. Returns 0 if not set
    UInt GetIsoOrder() const;
    
    //! Return anisotropic order.
    void GetAnisoOrder( StdVector<UInt>& order ) const;
    
    //! Return maximum integration order
    UInt GetMaxOrder() const;
    
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
    
  private:

    //! Vector containing the integration order
    std::array<UInt,3> order_;
    
    //! Flag if integration order is isotropic
    bool isIsotropic_;
    
    //! Flag if order is set
    bool isSet_;
    
  };
  
  
  //! external operator for comparing two IntegOrders
  bool operator<( const IntegOrder& a, const IntegOrder& b );

  /** necessary hash function for unordered_maps of key type involving IntegOrder */
  std::size_t hash_value(const IntegOrder& p);
} // namespace CoupledField

// Explicit specialization of std::hash must be injected into namespace std.
// This allows IntegOrder to be used as a key in std::unordered_map.
namespace std {
  template<> struct hash<CoupledField::IntegOrder> {
    std::size_t operator()(const CoupledField::IntegOrder& p) const {
      return CoupledField::hash_value(p);
    }
  };
}

namespace CoupledField {

  //! Class defining Numerical Integration

  //! The integration scheme returns the integration points / elements for a 
  //! given element shape. As additional parameters one can pass the 
  //! order (maybe anisotropic) and the type of integration used
  class IntScheme {
    
  public:

    /** enumeration with integration types - to be renamed into IntegrationType
    * UNDEFINED: only internal use / not define
    * GAUSS: Classic Gauss-Legendre integration points (quad, hex: tensor product of 1D points),
    *        which also supports anisotropic integration order.
    * GAUSS_ECO: are the "efficient" Gaussian quadrature weights
    *            (->Solin, Segeth, Dolezel, Higher-Order Finite Element Methods)
    * LOBATTO: Gauss-Lobatto integration points (used for spectral method)
    *            Higher-Order Finite Element Methods */
    typedef enum { UNDEFINED, GAUSS,  GAUSS_ECO, 
                   LOBATTO, CHEBYSHEV  } IntegMethod;
    static Enum<IntegMethod> IntegMethodEnum;
    
    typedef std::map<Elem::ShapeType,StdVector<LocPoint> > IntegrationPoints;
    typedef std::map<Elem::ShapeType,StdVector<Double> > IntegrationWeights;
    
    //! Constructor
    IntScheme();

    //! Destructor
    ~IntScheme();
  
    //! Get integration points / weights
    void GetIntPoints( Elem::ShapeType elemType,
                       IntegMethod method,
                       const IntegOrder& order,
                       StdVector<LocPoint>& intPts, 
                       StdVector<Double>& weights );
    
    //! Get integration points / weights for combinations of two orders
    void GetIntPoints( Elem::ShapeType elemType,
                       IntegMethod method1,
                       const IntegOrder& order1,
                       IntegMethod method2,
                       const IntegOrder& order2,
                       StdVector<LocPoint>& intPts, 
                       StdVector<Double>& weights );  

    //! Returns all defined integration points in a single vector
    void GetAllIntegrationPoints(StdVector< LocPoint >& points, Elem::ShapeType type); 

    //! Returns defined integration points in a single vector for a given order
    void GetIntegrationPoints(std::map<Integer, LocPoint >& points, Elem::ShapeType type,
                              IntegMethod method, const IntegOrder& order );

    //! Print list of available integration points
    std::string PrintList( bool detailed ) const;
    
private:

    //! Define specific integration set and save it to internal map (optional)
    void DefineIntPoints( Elem::ShapeType shapeType,
                          IntegMethod method, const IntegOrder& order,
                          StdVector<LocPoint>& points, 
                          StdVector<Double>& weights,
                          bool saveInternal = true);

    // ======================================================================
    //  Element Shape Specific Integration Points
    // ======================================================================
    //@{ \name Element Shape Specific Integration Rule Methods
    
    //! Define integration points / weights for line elements
    void DefineLinePoints( IntegMethod method, const IntegOrder& order,
                           StdVector<LocPoint>& points, 
                           StdVector<Double>& weights );

    //! Define integration points / weights for triangular elements
    void DefineTriaPoints( IntegMethod method, const IntegOrder& order,
                           StdVector<LocPoint>& points, 
                           StdVector<Double>& weights );

    //! Define integration points / weights for quadrilateral elements
    void DefineQuadPoints( IntegMethod method, const IntegOrder& order,
                           StdVector<LocPoint>& points, 
                           StdVector<Double>& weights );
    
    //! Define integration points / weights for tetrahedral elements
    void DefineTetPoints( IntegMethod method, const IntegOrder& order,
                           StdVector<LocPoint>& points, 
                           StdVector<Double>& weights );
    
    //! Define integration points / weights for hexahderal elements
    void DefineHexPoints( IntegMethod method, const IntegOrder& order,
                           StdVector<LocPoint>& points, 
                           StdVector<Double>& weights );
    
    //! Define integration points / weights for pyramidal elements
    void DefinePyraPoints( IntegMethod method, const IntegOrder& order,
                            StdVector<LocPoint>& points, 
                            StdVector<Double>& weights );
    
    //! Define integration points / weights for wedge / prism elements
    void DefineWedgePoints( IntegMethod method, const IntegOrder& order,
                           StdVector<LocPoint>& points, 
                           StdVector<Double>& weights );
    

    //@}

    // ======================================================================
    //  Auxilliary Methods
    // ======================================================================  
    //! Calculate integration points for triangles with Duffy transformation
    void CalcIntTria( IntegMethod, UInt order,StdVector<LocPoint>& points,
                      StdVector<Double>& weights );

    //! Calculate gauss integration points for tetrahedras by applying 3D Duffy transformation to hexahedrals
    void CalcIntTet( IntegMethod, UInt order,StdVector<LocPoint>& points,
                      StdVector<Double>& weights );

    //! Fill initial set of integration points up to a given order.
    void FillInitialIntegPoints(UInt order);

    //! Convert array representation of LocPoints/Weight-array to vectors

    //! This method converts the raw double-array representation of 
    //! integration weights and numbers to a vector representation of 
    //! LocPoints and weights.
    //! \param dim Spatial dimension of the element
    //! \param shape Shape of element
    //! \param data Actually [][2] for 1D (coord + weight) and [][3] 
    //!             for 2D and  [][4] for 3D
    //! \param points Integration points
    //! \param weights Integration weights

    void Convert( Elem::ShapeType shape, UInt nPoints, Double *data,
                  StdVector<LocPoint>& points, StdVector<Double>& weights);

    //! Calculate the Gauss-Lobatto points an weights for the given order
    void CalcGaussLobattoPointsWeights( UInt order,StdVector<Double>& points, 
                                        StdVector<Double>& weights );

    //! Calculate the Gauss-Legendre points an weights for the given order
    void CalcGaussLegendrePointsWeights( UInt order, StdVector<Double>& points, 
                                         StdVector<Double>& weights );

    /** typedef for map of integration points. We need a hash_value() function for IntegOrder */
    typedef std::unordered_map< IntegOrder, IntegrationPoints > IntPointMap;

    //! typedef for map of integration weights
    typedef std::unordered_map< IntegOrder, IntegrationWeights > IntWeightMap;

    //! Map with integration points for each element type According to the Template Parameter Integration Scheme
    std::map<IntegMethod, IntPointMap > intPoints_;

    //! Map with integration weights for each element type
    std::map<IntegMethod, IntWeightMap > intWeights_;

    //! stores the overall number of integration points defined in this class
    //! for each element type, zero based
    std::map<Elem::ShapeType,UInt> numIntPts_;

  };

} // namespace CoupledField

#endif /* INTTABLE_HH_ */
