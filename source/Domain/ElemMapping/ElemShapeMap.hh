#ifndef FILE_ELEMSHAPEMAP_HH_
#define FILE_ELEMSHAPEMAP_HH_

#include "boost/scoped_ptr.hpp"
#include "General/defs.hh"
#include "MatVec/Vector.hh"
#include "MatVec/Matrix.hh"
#include "Elem.hh"
#include "SurfElem.hh"
#include <climits>


namespace CoupledField {


  //! forward class declarations
  class Grid;
  class ElemShapeMap;
  class FeH1LagrangeExpl;
  class FeH1;
  class CoefFunction;
  
  // ===========================================================================
  //  C L A S S   LocPoint
  // ===========================================================================

  //! Simple struct representing an element local point

  //! This struct represents an element local point (xi, eta, zeta)
  //! which is either defined by its coordinates or by an integration 
  //! point
  struct LocPoint {

  public:

    //! Constructor
    LocPoint();
    
    //! Conversion from standard vector
    LocPoint( const Vector<Double>& vec);
    
    //! Return coordinate component of point
    inline Double& operator[]( UInt i) 
    { return coord[i];}

    //! Return coordinate component of point
    inline const Double& operator[]( UInt i) const 
    { return coord[i];}

    //! Alias for point, which is not explicitly represented by 
    //! a integration point
    enum {NOT_SET = INT_MAX};

    //! Number of corresponding integration point 
    Integer number;

    //! Coordinate of local point
    Vector<Double> coord;
  };
  
  // free outstream operator
  //! Overloading << for class LocPoint
  std::ostream& operator << ( std::ostream & , const LocPoint &);

  // ===========================================================================
  //  C L A S S   LocPointMapped
  // ===========================================================================

  //! Struct representing mapped element local information at a given point
  
  //! This class represents all geometry related information of a finite element
  //! at a given local point. The idea is to bundel up any geometry related
  //! information when it is computed first (e.g. within the integration loop
  //! of an integrator) and to pass it to any succeeding functions (e.g. within
  //! the bilinearforms, the mapped point is passed to the differential operators,
  //! which need the Jacobian matrix).
  //! 
  //! Thus, this struct basically contains the following data:
  //! 
  //! - LocPoint (point in element-local xi/eta/zeta coordinates)
  //! - ElementTransformation (mapping for reference -> physical domain)
  //! - Pointer to geometrical element
  //! - Jacobian Matrix, its inverse and determinant for given local point
  //! 
  //! If the mapped element is a surface element, it additionally contains
  //! the following information:
  //! 
  //! - Normal direction of the surface element in this point
  //! - Mapped local point of the surface element coordinates in the 
  //!   neighboring volume element
  //!
  //! Upon initialization, the following things are performed:
  //! - Calculation of Jacobian matrix and its inverse (volume elements only)
  //! - Calculation of its determinant, incorporating information about 
  //!   axi-symmetry or global model depth (only 2D)
  //!
  //! If the object gets initialized with a surface element, additionally the 
  //! following steps are performed:
  //! 
  //! - The correct volume neighbor element is selected, depending on the list
  //!   of regionIds passed to the Set()-method
  //! - Projection of surface element coordinates to the neighboring volume 
  //!   element (which is a LocPointMapped struct itself)
  //! - Calculation of normal direction, pointing out of the neighboring
  //!   volume element.  
  struct LocPointMapped {

  public:

    //! Constructor
    LocPointMapped();

    //! Initialize shape map with with given local point and shape map
    
    //! This constructor initializes the struct for a local point in a 
    //! general (volume) element.
    //! \param lp Local point to bet set
    //! \param esm ElemShapeMap, representing the mapping from reference to
    //!            physical domain    
    void Set( const LocPoint& lp, shared_ptr<ElemShapeMap> esm);
    
    //! Set method for a local point of a surface element
    
    //! This constructor initializes the struct for a local point in a
    //! surface element. To determine the "correct" volume neighbor, a set of
    //! regionIds of neighboring volume regions is passed. The normal
    //! direction will then be calculated to point out of those volume
    //! elements.
    //! \param lp Local point to bet set
    //! \param esm ElemShapeMap, representing the mapping from reference to
    //!            physical domain    
    //! \param volRegions Set containing the regionIds of the neighboring
    //!                   volume regions.
    void Set( const LocPoint& lp, shared_ptr<ElemShapeMap> esm,
              const std::set<RegionIdType>& volRegions ); 

    //! Shape map for this element
    shared_ptr<ElemShapeMap> shapeMap;

    //! Pointer to element
    const Elem * ptEl;

    //! Element local point
    LocPoint lp;

    //! Jacobian matrix in local point @ref lp
    Matrix<Double> jac;

    //! Inverse of Jacobian matrix in local point @ref lp
    Matrix<Double> jacInv;

    //! Jacobian determinant in in local point @ref lp
    //! @note: If the underlying geometry is axi-symmetric, the Jacobian 
    //! determinant is scaled by \f$ 2\pi r \f$ with r being the 
    //! global x coordinate of the local integration point @ref lp.
    Double jacDet;
    
    // =======================================
    //  Surface element specific section
    // =======================================
    //@{ \name section for surface elements
    
    //! Flag indicating, if surface element information is set
    bool isSurface;
    
    //! Normal direction of surface element (only set in case of a surface element)
    Vector<Double> normal;
    
    //! Mapped local point of the correct neighboring volume element
    
    //! If the struct was initialized with a surface element and the set of
    //! neighboring volume regions, this object represents the projected 
    //! local point of the surface element in the neighboring volume
    //! element. This is necessary, e.g. to calculate normal directions or
    //! derivatives in normal direction.
    shared_ptr<LocPointMapped> lpmVol;
    
    //@}
  };


  
  // ===========================================================================
  //  C L A S S   ElemshapeMap
  // ===========================================================================
  
  //! Base class for representing geometric shapes of elements 
  //!
  //! This class represents the geometric representation of 
  //! mesh elements, i.e. it can map local/global coordinates,
  //! calculate the Jacobian matrix / determinant and the normal direction
  //! for surface elements
  //! In most cases the concrete implementation will be a 1st/2nd order
  //! Lagrangian element, as provided by most mesh generators.
  //!
  //! Collaboration: 
  //! The class gets instantiated in the Grid-class, so the grid class needs
  //! a method GetShapeMap(const Elem* ). The Elem class therefore just receives
  //! an additional flag of type ShapeType. This may not be necessary at
  //! all, as the grid could store this information ....
  //! 
  
  class ElemShapeMap {
  
  public:

    //! Define enumeration data type
    typedef 
    enum { NO_TYPE, LAGRANGE, LAGRANGE_BLENDED, 
      ANALYTICAL } ShapeMapType;  

    //! Enum for shape map type
    static Enum<ShapeMapType> shapeMapType;

    //! Constructor
    //! \param ptGrid input Pointer to Grid
    //! \param isUpdated input Use of updated Lagrangian formulation
    ElemShapeMap( Grid* ptGrid );

    //! Destructor
    virtual ~ElemShapeMap();
    
    //! Return pointer to grid
    Grid* GetGrid() {
      return ptGrid_;
    }

    //! Set current element
    //! \param ptElem output Current element
    virtual void SetElem( const Elem* ptElem,
                          bool isUpdated = false );
    
    //! Get current element
    //! \param ptElem output Current element
    virtual const Elem* GetElem( ) {
      return ptElem_;
    }
    
    //! Get current surface element
    
    //! This method returns a pointer to the internal element in form of a
    //! surface element pointer. This pointer only is valid, if the set element
    //! was in fact a surface element
    //! \param ptElem output Current element
    virtual const SurfElem* GetSurfElem( ) {
      return ptSurfElem_;
    }

    //! Query type of shape-type
    ShapeMapType GetShapeMapType() const {return type_;}

    // ---------------------------------------------------
    //   Coordinate Mapping
    // ---------------------------------------------------

    //! Return axi-symmetric flag
    bool IsAxi() const {return isAxi_;}
    
    //! Return if updated geometry is used
    bool IsUpadet() const {return isUpdated_;}

    //! Calculate Local -> Global Mapping
    //! \param globPoint output Global point in cartesian coordinates
    //! \param locPoint input Element local point
    virtual void Local2Global( Vector<Double>& globPoint, 
                               const LocPoint& locPoint ) = 0;

    //! Calculate Global -> Local Mapping

    //! Calculates the mapping global->local using a Newton Raphson method
    //! \param locPoint output Element local point
    //! \param globPoint input Global point in cartesian coordinates
    virtual void Global2Local( Vector<Double>& locPoint, 
                               const Vector<Double>& globPoint ) = 0;

    //! Calculate global midpoint of element
    //! \param midPoint output Global element midpoint
    virtual void GetGlobMidPoint( Vector<Double>& midPoint ) = 0;

    //! Calculate volume of element
    virtual Double CalcVolume( ) = 0;

    //! Calculate normal of element
    //! \param normal output Normal vector in global coordinates
    //! \param lp input Element local point
    //! \note The direction of the normal has undefined sign!
    virtual void CalcNormal( Vector<Double>& normal, 
                             const LocPoint& lp ) = 0;
    

    //! Returns surface element normal with defined orientation

    //! Calculates the surface normal pointing OUT OF the neighboring
    //! volume element
    //! \param n (out) normal vector
    //! \param lp (input) Element local point
    //! \param volElem (in) volume element
    virtual void CalcNormalOutOfVol( Vector<Double> & normal,
                                     const LocPoint& lp,
                                     const Elem & volElem ) = 0;

    //! Calculates corresponding volume point of neighboring surfaces

    //! For a given surface element  this
    //! method calculates the local volume-coordinates out of the given
    //! local surface-coordinates, which have one less dimension.
    //! This can be used to get the corresponding volume coordinates of
    //! the integration points of a surface. Therefore it calculates
    //! on which side of the volume element the surface element lies
    //! and creates the according volume point.
    //! In addition it returns the local normal direction of the side
    //! of the volume element, which this element is neighhoring.
    /*!
       \param surfConnect (input) Node numbers of surface element
       \param volConnect (input) Node numbers of volume element
       \param surfIntPoint (input) Surface integration point, which gets mapped
                                   onto the volume element
       \param volIntPoint (output) Corresponding volume integration point
       \param locNormal (output) Normal direction of surface element in local 
                                 coordinates of the volume element
     */
    virtual void GetLocalIntPoints4Surface( const StdVector<UInt> & surfConnect,
                                            const LocPoint & surfIntPoint,
                                            LocPoint & volIntPoint,
                                            Vector<Double>& locNormal ) = 0;
    
    //! Returns whether a given local coordinate is within this element
    //! \param point input Local point
    //! \param tolerance input Additioanl (relative) tolerance
    //! \return flag if point is inside the element
    virtual bool CoordIsInsideElem( const Vector<Double>& point,Double tolerance = 0.0 ) = 0;

    //! Calculate the diameter vector of the element.
    //! Handles the element by itself and no axis-symmetric case.
    //! \param diameter output vector with only positive entries. */
    virtual void CalcDiameter( Vector<Double>& diameter ) = 0;

    //! Calculates the barycenter of the element. For 3D and 2D (z=0 then)
    //! This is the average of every coordinate.
    //! \param barycenter the output. In 2D z=0. If you want 2D then overload */
    virtual void CalcBarycenter( Vector<Double>& baryCenter )  = 0;

    //! Compute length of edge with maximal/minimal size
    virtual void GetMaxMinEdgeLength( Double& max, Double& min)  = 0;

    //! Compute the average side length for all dimensions. Clearly works only 
    //! for quadrilaterals and hexahedrons.
    //! \param edges_out array with index 0 for x, ... */
    virtual void GetEdgeLength( StdVector<Double>& edges_out) = 0;

    // ---------------------------------------------------
    //   Jacobian 
    // ---------------------------------------------------

    //! Calculate Jacobian Matrix
    /*! \param jac output Jacobian Matrix
    \f$ J = \left( \begin{array}{ccc} x_{\xi} & x_{\eta} & x_{\zeta} \\
    y_{\xi} & y_{\eta} & y_{\zeta} \\
    z_{\xi} & z_{\eta} & z_{\zeta} 
    \end{array}\right) \f$ 
    */
    virtual void CalcJ( Matrix<Double>& jac, 
                        const LocPoint& ip ) = 0;

    //! Calculation of Jacobian determinant
    /*!
         \param LCoord (input) Local Coordinates of evaluation point
         \param CornerCoords (input) Coordinates of element corners
          \f$ [ \left( \begin{array}{ccc} x_{1} & x_{2} & \cdots \\ y_{1} & y_{2} & \cdots \\
         \cdots & \cdots & \cdots \end{array} \right)  \f$ ]
     */
    virtual Double CalcJDet( Matrix<Double>& jac, 
                             const LocPoint& ip ) = 0;
    
  protected:

    //! Type of shape mapping
    ShapeMapType type_;

    //! Pointer to grid
    Grid *ptGrid_;

    //! Flag for axisymmtry
    bool isAxi_;

    //! Flag if updated coordinates are uses (updated Lagrange)
    bool isUpdated_;

    //! Depth of elements (only for 2D)
    Double depth_;

    //! Pointer to current element
    const Elem* ptElem_;
    
    //! Pointer to current surface element (if valid)
    const SurfElem* ptSurfElem_;

  };
  
  // ========================================================================
  //  Lagrangian Element Shape Map
  // ========================================================================
  
  
  //! Lagrangian representation for finite elements
  
  //! This class represents the classical 1st/2nd order Lagrangian element
  //! shape map, i.e. the shape of the element is transformed using the 
  //! bi/trilinear shape functions of the Lagrangian elements.
  class LagrangeElemShapeMap : public ElemShapeMap {

  public:
    //! Constructor
    LagrangeElemShapeMap( Grid* ptGrid );

    //! Destructor
    ~LagrangeElemShapeMap();

    //! @copydoc ElemShapeMap::SetElem
    virtual void SetElem( const Elem* ptElem, bool isUpdated = false );

    // ---------------------------------------------------
    //   Coordinate Mapping
    // ---------------------------------------------------

    //! @copydoc ElemShapeMap::Local2Global
    void Local2Global( Vector<Double>& globPoint, 
                       const LocPoint& lp );

    //! @copydoc ElemShapeMap::Global2Local
    void Global2Local( Vector<Double>& locPoint, 
                       const Vector<Double>& glob );
    
    //! Specialized version for quad4
    void Global2LocalQuad4( Vector<Double>& locPoint,
                               const Vector<Double>& glob );

    //! General version
    void Global2LocalGeneral( Vector<Double>& locPoint,
                               const Vector<Double>& glob );

    //! @copydoc ElemShapeMap::GetGlobMidPoint
    void GetGlobMidPoint( Vector<Double>& midPoint );

    //! @copydoc ElemShapeMap::CalcVolume
    Double CalcVolume( );

    //! @copydoc ElemShapeMap::CalcNormal
    void CalcNormal( Vector<Double>& normal, 
                     const LocPoint& lp );
    
    //! @copydoc ElemShapeMap::CalcNormalOutOfVol
    void CalcNormalOutOfVol( Vector<Double> & normal,
                             const LocPoint& lp,
                             const Elem & volElem );
    
    //! @copydoc ElemShapeMap::GetLocalIntPoints4Surface
    void GetLocalIntPoints4Surface( const StdVector<UInt> & surfConnect,
                                    const LocPoint & surfIntPoint,
                                    LocPoint & volIntPoint,
                                    Vector<Double>& locNormal );
    
    //! @copydoc ElemShapeMap::CoordIsInsideElem
    bool CoordIsInsideElem( const Vector<Double>& point,Double tolerance = 0.0 );

    //! @copydoc ElemShapeMap::CalcDiameter
    void CalcDiameter( Vector<Double>& diameter );

    //! @copydoc ElemShapeMap::CalcBarycenter
    void CalcBarycenter( Vector<Double>& baryCenter );

    //! @copydoc ElemShapeMap::GetMaxMinEdgeLength
    void GetMaxMinEdgeLength( Double& max, Double& min);

    //! @copydoc ElemShapeMap::GetEdgeLength
    void GetEdgeLength( StdVector<Double>& edges_out);

    // ---------------------------------------------------
    //   Jacobian 
    // ---------------------------------------------------

    //! @copydoc ElemShapeMap::CalcJ
    void CalcJ( Matrix<Double>& jac, 
                const LocPoint& ip );
    
    //! @copydoc ElemShapeMap::CalcJ
    Double CalcJDet( Matrix<Double>& jac, 
                     const LocPoint& ip );
    
    //! Initialize static members
    static void InitStaticMembers();
    
  protected:

    //! Pointer to H1 element of lower order
    FeH1LagrangeExpl* ptFe_;
    
    //! Map with pointers to reference elements
    static std::map<Elem::FEType, FeH1LagrangeExpl* > feMap_;

    //! Nodal coordinates
    Matrix<Double> coords_;
    
    //! Shape of reference element
    ElemShape shape_;
    
  };

}

#endif /* ELEMSHAPEMAP_HH_ */
