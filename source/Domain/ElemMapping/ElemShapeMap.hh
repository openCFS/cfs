#ifndef FILE_ELEMSHAPEMAP_HH_
#define FILE_ELEMSHAPEMAP_HH_

#include "boost/scoped_ptr.hpp"
#include "General/defs.hh"
#include "MatVec/Vector.hh"
#include "MatVec/Matrix.hh"
#include "Elem.hh"
#include "SurfElem.hh"

#include "Utils/ThreadLocalStorage.hh"

#include <climits>


namespace CoupledField {


  //! forward class declarations
  class Grid;
  class ElemShapeMap;
  class FeH1LagrangeExpl;
  class FeH1;
  class BaseFE;
  class CoefFunction;
  class IntScheme;

  
  // ==========================================================================
  //  C L A S S   LocPoint
  // ==========================================================================

  //! Multiple usage struct:
  //! - standard case: representing an element local point
  //! - rare case: transport coordinate for nodal points in global system
  //! - rare case: transport number of element or node

  //! This struct represents an element local point (xi, eta, zeta)
  //! which is either defined by its coordinates or by an integration 
  //! point created by the integration scheme class (\see IntScheme).
  //! In this case the point is uniquely defined by its number 
  //! (IntScheme::number).
  //! If an arbitrary point is set, which is not defined by a numbered
  //! integration point, the number is set to NOT_SET. 
  struct LocPoint {

  public:

    //! Constructor
    LocPoint() {};
    
    //! Conversion from standard vector
    LocPoint( const Vector<Double>& vec);
    
    //! Return coordinate component of point
    inline Double& operator[]( UInt i) { return coord[i];}

    //! Return coordinate component of point
    inline const Double& operator[]( UInt i) const { return coord[i];}

    //! Alias for point, which is not explicitly represented by 
    //! a integration point
    enum {NOT_SET = -1};

    //! Number of corresponding integration point 
    Integer number = NOT_SET;

    //! Coordinate of local point
    Vector<Double> coord;
  };
  
  //! Overloading << for class LocPoint
  std::ostream& operator << ( std::ostream & , const LocPoint &);

  // ===========================================================================
  //  C L A S S   LocPointMapped
  // ===========================================================================

  //! Struct representing mapped element local information at a given point
  
  //! The struct servers the following usage:
  //! - provide FE relatated information for an element, including local/global transformation
  //! - (rarely) transport information about a nodal point of interest in global coordinate space
  //! - (rarely) transport about an element or nodal number of interest.
  //!
  //! For the rarely usage, lp is misused, lpm is overly complicated and we do not
  //! identify what the intended usage is, just having shapeMap not set gives an indication

  //! This class represents all geometry related information of a finite element
  //! at a given local point. The idea is to bundle up any geometry related
  //! information when it is computed first (e.g. within the integration loop
  //! of an integrator) and to pass it to any succeeding functions (e.g. within
  //! the bilinearforms, the mapped point is passed to the differential operators,
  //! which need the Jacobian matrix).
  //!
  //! Thus, this struct basically contains the following data:
  //! 
  //! - LocPoint (point in element-local xi/eta/zeta coordinates)
  //!   Alternatively number and coord are used for the rarely cases mentioned above
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
    //! \param weight Integration weight (should be set to 0.0 if not used 
    //!               within an integration loop)
    void Set(const LocPoint& lp, shared_ptr<ElemShapeMap> esm, Double weight = 0.0);
    
    //! Initialize shape map with with given local point and shape map

    //! This constructor initializes the struct for a local point in a
    //! general (volume) element.
    //! \param lp Local point to bet set
    //! \param esm ElemShapeMap, representing the mapping from reference to
    //!            physical domain
    //! \param weight Integration weight (should be set to 0.0 if not used
    //!               within an integration loop)
    //! \param cornerCoord explicit definition of corner coordiantes of element
    void Set(const LocPoint& lp, shared_ptr<ElemShapeMap> esm,
                 Double weight, Matrix<Double>& cornerCoord);

    //! This method allows to set information specific to a surface
    //! element, in case it was initialized with only the volume information;
    //! This constructor initializes the struct for a local point in a
    //! surface element. To determine the "correct" volume neighbor, a set of
    //! regionIds of neighboring volume regions is passed. The normal
    //! direction will then be calculated to point out of those volume
    //! elements.
    //! \note The search for the correct volume neighbor is just performed,
    //!       if a new element is set. If only the local point within one
    //!       element is set, the neighbor information is re-used. 
    //! \param lp Local point to bet set
    //! \param esm ElemShapeMap, representing the mapping from reference to
    //!            physical domain    
    //! \param volRegions Set containing the regionIds of the neighboring
    //!                   volume regions.
    //! \param weight Integration weight (should be set to 0.0 if not used 
    //!               within an integration loop)
    void SetWithSurface( const LocPoint& lp, shared_ptr<ElemShapeMap> esm, const std::set<RegionIdType>& volRegions, Double weight );

    //! This method allows to set information specific to a surface
    //! element, in case it was initialized with only the volume information;
    //! \param volRegions (in) Set containing the regionIds of the neighboring volume regions.
    //! \param volRegid (in) Allows to explicitly specify a volume region in case of the surface being an interface
    void SetSurfInfo( const std::set<RegionIdType>& volRegions, const RegionIdType volRegid = NO_REGION_ID );

    //! set, if jacobi determinat should be checked; standard is YES
    void SetCheckJacobi(bool check) {
    		checkJacobi_ = check;
    };

    //! Specialized version for NMG points
    void SetMortar( const LocPoint& lp, shared_ptr<ElemShapeMap> esm,
                                          Double weight, bool useMaster);

    //! Return shape map
    const shared_ptr<ElemShapeMap> GetShapeMap() const {return shapeMap;}

    /** service function which obtains a global point based on a LocPoint
     * @parm loc if not given the own lp is used
     * @param fallback if shapeMap is not set, assue we want the rare nodal case and use lp->number
     *        as nodal number and get the coordinate from the mesh
     * @param update for the fallback case, obtain the updated coordinate from the mesh
     * @return the given coord as convenience */
    Vector<double>& GetGlobal(Vector<double>& coord, const LocPoint* loc = NULL, bool fallback = true, bool update = false) const;

    /** same as the other GetGlobal() helper, but we (abuse) lp::coord, either given or internal lp */
    Vector<double>& GetGlobal(LocPoint* loc = NULL, bool fallback = true, bool update = false) {
      return GetGlobal(lp.coord, loc, fallback, update);
    }

    //! Shape map for this element
    shared_ptr<ElemShapeMap> shapeMap;

    //! Pointer to element
    const Elem* ptEl;

    //! Element local point, which is of multiple usage
    LocPoint lp;
    
    //! Integration weight
    Double weight;

    //! Jacobian matrix in local point @ref lp
    Matrix<Double> jac;

    //! Inverse of Jacobian matrix in local point @ref lp
    Matrix<Double> jacInv;

    //! Jacobian determinant in local point @ref lp
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
    
    //! check Jacobi detrminat
    bool checkJacobi_;

    //! Normal direction of surface element (only set in case of a surface element)
    
    //! In case the struct gets initialized with a surface element, this member
    //! contains the surface normal, which points out of the "correct" neighbouring
    //! volume element (\see LocPointMapped::lpmVol).
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
  
  class ElemShapeMap : public boost::enable_shared_from_this<ElemShapeMap>{
  public:

    //! Define enumeration data type
    typedef 
    enum { NO_TYPE, LAGRANGE, LAGRANGE_BLENDED, ANALYTICAL } ShapeMapType;

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

    //! Set current element.
    //! Be very careful with it,
    //! @see Grid::GetElemShapeMap()
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
    bool IsUpdated() const {return isUpdated_;}

    //! Calculate Local -> Global Mapping
    //! \param globPoint output Global point in Cartesian coordinates
    //! \param locPoint input Element local point
    virtual void Local2Global( Vector<Double>& globPoint, 
                               const LocPoint& locPoint ) = 0;

    //! Calculate Global -> Local Mapping

    //! Calculates the mapping global->local using a Newton Raphson method
    //! \param locPoint output Element local point
    //! \param globPoint input Global point in Cartesian coordinates
    virtual void Global2Local( Vector<Double>& locPoint, 
                               const Vector<Double>& globPoint ) = 0;

    //! Calculate global midpoint of element
    //! \param midPoint output Global element midpoint
    virtual void GetGlobMidPoint( Vector<Double>& midPoint ) = 0;

    //! Calculate volume of element
    virtual Double CalcVolume(bool useDepth = false) = 0;

    //! Calculate normal of element
    
    //! This method calculate the normal vector of a surface element in the
    //! given local point. The resulting vector will point OUT of the 
    //! first volume neighbor of the surface element!
    //! \param normal output Normal vector in global coordinates
    //! \param lp input Element local point
    virtual void CalcNormal( Vector<Double>& normal, 
                             const LocPoint& lp ) = 0;
    
    //!This method checks if a given local point is on the surface
    //! of the element. If so, it returns true and calculates the
    //! surface normal according to the point. If not, it returns
    //! false
    virtual bool CalcNormalOutOfVolume(Vector<Double> & normal,
                                           const LocPoint & lp,
                                           const Elem* volElem,
                                           const SurfElem* edgeFaceElem)=0;

    //! This method calculates a parallel projection of the given point onto the surface element
    //! with respect to the given direction
    //! \param direction vector detrmining the direction of the projection
    //! \param point point to be projected
    virtual void TranslatePointOntoSurface(const Vector<Double>& direction, Vector<Double>& point) = 0;


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
    
    //! Calculate local directions of surface / edge w.r.t. to volume element
    
    //! This method maps a surface element to an edge / surface of the 
    //! the neighboring volume element (= this instance) and returns
    //! a vector containing the local direction(s) of this surface element
    //! w.r.t. to the volume.
    virtual void MapSurfLocDirs( const Elem* ptSurfElem, 
                                 StdVector<UInt>& surfLocDirs ) = 0;
    
    //! Returns whether a given local coordinate is within this element
    //! \param point input Local point
    //! \param tolerance input Additional (relative) tolerance
    //! \return flag if point is inside the element
    virtual bool CoordIsInsideElem( const Vector<Double>& point,Double tolerance = 0.0 ) = 0;

    //! Calculate the diameter vector of the element.
    //! Handles the element by itself and no axi-symmetric case.
    //! \param diameter output vector with only positive entries. */
    virtual void CalcDiameter( Vector<Double>& diameter ) = 0;

    //! Calculates the barycenter of the element. For 3D and 2D (z=0 then)
    //! This is the average of every coordinate.
    //! \param barycenter the output. In 2D z=0. If you want 2D then overload */
    virtual void CalcBarycenter(Point& baryCenter )  = 0;

    //! Compute length of edge with maximal/minimal size
    virtual void GetMaxMinEdgeLength( Double& max, Double& min)  = 0;

    //! Compute the average side length for all dimensions. Clearly works only 
    //! for quadrilaterals and hexahedrons.
    //! \param edges_out array with index 0 for x, ... */
    virtual void GetEdgeLength( StdVector<Double>& edges_out) = 0;
    
    //! Compute extension of the element in local directions
    
    //! This method computes the extension of the element (w.r.t. to
    //! global coordinates) and returns a them in a vector, mapped to 
    //! local coordinate directions. Internally, the extension is determined
    //! using the edge sizes and for each edge the according local coordinate
    //! direction is taken. Thus the diagonal extension of the element is not
    //! considered.
    //! \param extension Contains the extension (in global coordinates) of the
    //!                  element w.r.t. local coordinate directions.
    virtual void GetExtensionLocalDir( Vector<Double>& extension ) = 0;

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

    //! Calculation of Jacobian with given coordinates
    virtual void CalcJ( Matrix<Double>& jac,
       		            const LocPoint& ip,
   				        Matrix<Double>& cornerCoords) = 0;

    //! Calculation of Jacobian determinant
    /*!
         \param LCoord (input) Local Coordinates of evaluation point
         \param CornerCoords (input) Coordinates of element corners
          \f$ [ \left( \begin{array}{ccc} x_{1} & x_{2} & \cdots \\ y_{1} & y_{2} & \cdots \\
         \cdots & \cdots & \cdots \end{array} \right)  \f$ ]
     */
    virtual Double CalcJDet( Matrix<Double>& jac, 
                             const LocPoint& ip, bool useDepth = false ) = 0;

    virtual void SetModelDepth(Double depth) = 0;
    
    virtual Double GetModelDepth() = 0;
    
    //! obtain pointer to geometric reference element
    virtual BaseFE* GetBaseFE()  = 0;

    /** for debugging purpos */
    virtual std::string ToString() const;

  protected:

    //! Type of shape mapping
    ShapeMapType type_ = NO_TYPE;

    //! Pointer to grid
    Grid *ptGrid_ = NULL;

    //! Flag for axisymmetry
    bool isAxi_ = false;

    //! Flag if updated coordinates are uses (updated Lagrange)
    bool isUpdated_ = false;

    //! Depth of elements (only for 2D)
    Double depth_ = 1.0;

    //! Pointer to current element
    const Elem* ptElem_ = NULL;
    
    //! Pointer to current surface element (if valid)
    const SurfElem* ptSurfElem_ = NULL;

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

    /** @copydoc ElemShapeMap::SetElem
     * this sets coords_ */
    virtual void SetElem( const Elem* ptElem, bool isUpdated = false );

    /** take some info from the element but use the given coordinates, not the one from elem */
    void SetElem(const Elem* elem, const Matrix<double>& coords);

    // ---------------------------------------------------
    //   Coordinate Mapping
    // ---------------------------------------------------

    //! @copydoc ElemShapeMap::Local2Global
    void Local2Global( Vector<Double>& globPoint, 
                       const LocPoint& lp );

    //! @copydoc ElemShapeMap::Global2Local
    void Global2Local( Vector<Double>& locPoint, 
                       const Vector<Double>& glob );
    
    //! @copydoc ElemShapeMap::GetGlobMidPoint
    void GetGlobMidPoint( Vector<Double>& midPoint );

    //! @copydoc ElemShapeMap::CalcVolume
    Double CalcVolume(bool useDepth = false);

    //! @copydoc ElemShapeMap::CalcNormal
    void CalcNormal( Vector<Double>& normal, 
                     const LocPoint& lp );
    
    //! @copydoc ElemShapeMap::CalcNormalOutOfVolume
    bool CalcNormalOutOfVolume(Vector<Double> & normal,
                               const LocPoint & lp,
                               const Elem* volElem,
                               const SurfElem* edgeFaceElem);

    //! @copydoc ElemShapeMap::TranslatePointOntoSurface
    void TranslatePointOntoSurface(const Vector<Double>& direction, Vector<Double>& point);

    //! @copydoc ElemShapeMap::GetLocalIntPoints4Surface
    void GetLocalIntPoints4Surface( const StdVector<UInt> & surfConnect,
                                    const LocPoint & surfIntPoint,
                                    LocPoint & volIntPoint,
                                    Vector<Double>& locNormal );
    
    //! @copydoc ElemShapeMap::MapSurfLocDirs
    void MapSurfLocDirs( const Elem* ptSurfElem, 
                         StdVector<UInt>& surfLocDirs );
       
    
    //! @copydoc ElemShapeMap::CoordIsInsideElem
    bool CoordIsInsideElem( const Vector<Double>& point,Double tolerance = 0.0 );

    //! @copydoc ElemShapeMap::CalcDiameter
    void CalcDiameter( Vector<Double>& diameter );

    //! @copydoc ElemShapeMap::CalcBarycenter
    void CalcBarycenter(Point& baryCenter );

    //! @copydoc ElemShapeMap::GetMaxMinEdgeLength
    void GetMaxMinEdgeLength( Double& max, Double& min);

    //! @copydoc ElemShapeMap::GetEdgeLength
    void GetEdgeLength( StdVector<Double>& edges_out);
    
    //! @copydoc ElemShapeMap::GetMinExtensionLocalDir
    void GetExtensionLocalDir( Vector<Double>& extension );

    // ---------------------------------------------------
    //   Jacobian 
    // ---------------------------------------------------

    //! @copydoc ElemShapeMap::CalcJ
    void CalcJ( Matrix<Double>& jac, 
                const LocPoint& ip );
    
    //! @copydoc ElemShapeMap::CalcJ with given coordinates
    void CalcJ( Matrix<Double>& jac,
    		    const LocPoint& ip,
				Matrix<Double>& cornerCoords);

    //! @copydoc ElemShapeMap::CalcJ
    Double CalcJDet( Matrix<Double>& jac, 
                     const LocPoint& ip, bool useDepth = false );
    
    void SetModelDepth(Double depth);
    
    Double GetModelDepth();
    
    //! @copydoc ElemShapeMap::GetBaseFE
    virtual BaseFE* GetBaseFE();

    virtual std::string ToString() const;

  protected:

    
    // ---------------------------------------------------
    //   Specialized Global -> Local Algorithms 
    // ---------------------------------------------------
    //@{ \name Specialized Global -> Local Mapping Algorithms 
    
    //! Return local search directions and Jacobian determinant
    
    //! This method calculates for the given, full Jacobian matrix
    //! the determinant and returns the local search direction. 
    //! It takes into account the possible combinations (globDir,elemDir).
    //! \param delta_xi (out) Local search direction
    //! \param delta_x (in) Global search direction
    //! \param jac (in) Jacobian matrix
    Double GetLocDirJac( Vector<Double>& delta_xi, 
                         const Vector<Double>& delta_x,
                         const Matrix<Double>& jac );
    
    //! check if the point coincides with corner nodes
    //! this function is introduced as the usual non-linear iterations
    //! sometimes mess up if the requested point coincides with a element node
    //! tolerance is 1e-10!
    bool Global2LocalOnNode(Vector<Double>& locPoint,
                            const Vector<Double>& glob);

    //! Specialized version for quad4
    void Global2LocalQuad4( Vector<Double>& locPoint,
                            const Vector<Double>& glob );

    //! Orignial version from duester skript
    void Global2LocalDuester(Vector<Double>& locPoint,
                             const Vector<Double>& globalPoint);

    //! For elements using barycentric coordinates right now limited to straight edges
    void Global2LocalBarycentric( Vector<Double>& locPoint,
                                  const Vector<Double>& globalPoint );

    //! Version for linear line elements
    void Global2LocalLine2(Vector<Double>& locPoint,
                           const Vector<Double>& globalPoint);

    //! General version
    void Global2LocalGeneral( Vector<Double>& locPoint,
                              const Vector<Double>& glob );

    //@}
    
    //! Pointer to H1 element of lower order
    FeH1LagrangeExpl* ptFe_;

    //! Helper struct containing the reference element
    class LagrangeMapSingleton
    { 
    private: 
      LagrangeMapSingleton();   
      LagrangeMapSingleton(const LagrangeMapSingleton&);            
      LagrangeMapSingleton& operator=(const LagrangeMapSingleton&);
      ~LagrangeMapSingleton();
    public: 
      static LagrangeMapSingleton& getInstance();

      //! Map containing the reference elements
      TLMap<Elem::FEType, FeH1LagrangeExpl* > feMap_;
    }; 

    //! Nodal coordinates
    Matrix<Double> coords_;
    
    //! Shape of reference element
    const ElemShape* shape_;
    
    //! Pointer to integration scheme
    shared_ptr<IntScheme> intScheme_;
    
    //! Pointer to class containing the reference elements
    LagrangeMapSingleton & elems_;

    //! Cached entry for local derivative (used for Jacobian)
    Matrix<Double> deriv_;

    //! Cached entry for shape functions (used for local -> global)
    Vector<Double> shFnc_;
  };

}

#endif /* ELEMSHAPEMAP_HH_ */
