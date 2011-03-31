// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
#ifndef FILE_ELEMSHAPEMAP_HH_
#define FILE_ELEMSHAPEMAP_HH_

#include "General/defs.hh"
#include "MatVec/vector.hh"
#include "MatVec/matrix.hh"
#include "Domain/elem.hh"


namespace CoupledField {


  //! forward class declarations
  class Grid;
  class ElemShapeMap;
  class FeH1LagrangeExpl;
  class FeH1;
  
  // ===========================================================================
  //  C L A S S   LocPoint
  // ===========================================================================

  //! Simple struct representing an element local point

  //! This struct represents an element local point (xi, eta, zeta)
  //! which is either defined by its coordinates or by an ingegration 
  //! point
  struct LocPoint {

  public:

    //! Constructor
    LocPoint();
    
    //! Conversion from standard vector
    LocPoint( const Vector<Double>& vec);

    //! Alias for point, which is not explicitly represented by 
    //! a integration point
    enum {NOT_SET = -1};

    //! Number of corresponding integration point 
    Integer number;

    //! Coordinate of local point
    Vector<Double> coord;
  };


  // ===========================================================================
  //  C L A S S   LocPointMapped
  // ===========================================================================

  //! Represents a local point in a mapped finite element

  //! This class contains the following information:
  //! - LocPoint
  //! - ElementTransformation
  //! - Jacobian Matrix for local point 
  //! 
  //! Upon creation it calculates:
  //! - Jacobian matrix
  //! - Its inverse
  //! - Its determinant
  //! 
  //! The idea is as follows:
  //! - This struct gets created within CalcElementMatrix() using the Element
  //!   transformation, the real element and the current integration point
  //! - After that it can be passed to the B-Operator  
  struct LocPointMapped {

  public:

    //! Constructor
    LocPointMapped();

    //! Initialize shape map with with given local point and shape map    
    void Set( const LocPoint& lp, shared_ptr<ElemShapeMap> esm);

    //! Shape map for this element
    shared_ptr<ElemShapeMap> shapeMap;

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
    //! @param ptGrid input Pointer to Grid
    //! @param isUpdated input Use of updated Lagrangian formulation
    ElemShapeMap( Grid* ptGrid );

    //! Destructor
    virtual ~ElemShapeMap();

    //! Set current element
    //! @param ptElem output Current element
    virtual void SetElem( const Elem* ptElem,bool isUpdated = false  );
    
    //! Set current element
    //! @param ptElem output Current element
    virtual const Elem* GetElem( ) {
      return ptElem_;
    }

    //! Query type of shape-type
    ShapeMapType GetShapeMapType() const {return type_;}

    // ---------------------------------------------------
    //   Coordinate Mapping
    // ---------------------------------------------------

    //! Return axi-symmetric flag
    bool IsAxi() const {return isAxi_;}

    //! Calculate Local -> Global Mapping
    //! @param globPoint output Global point in cartesian coordinates
    //! @param locPoint input Element local point
    virtual void Local2Global( Vector<Double>& globPoint, 
                               const LocPoint& locPoint ) = 0;

    //! Calculate Global -> Local Mapping

    //! Calculates the mapping global->local using a Newton Raphson method
    //! @param locPoint output Element local point
    //! @param globPoint input Global point in cartesian coordinates
    virtual void Global2Local( Vector<Double>& locPoint, 
                               const Vector<Double>& globPoint ) = 0;

    //! Calculate global midpoint of element
    //! @param midPoint output Global element midpoint
    virtual void GetGlodMidPoint( Vector<Double>& midPoint ) = 0;

    //! Calculate volume of element
    virtual Double CalcVolume( ) = 0;

    //! Calculate normal of element
    //! @param normal output Normal vector in global coordinates
    //! @param lp input Element local point
    virtual void CalcNormal( Vector<Double>& normal, 
                             const LocPoint& lp ) = 0;

    //! Returns whether a given local coordinate is within this element
    //! @param lps input Vector containing element local points
    //! @param tolerance input Additioanl (relative) tolerance
    //! @param coordsInside output Vector with flags
    virtual bool CoordIsInsideElem( const StdVector<LocPoint>& lps,
                                    const Double tolerance,
                                    StdVector<bool>& coordsInside ) = 0;

    //! Calculate the diameter vector of the element.
    //! Handles the element by itself and no axis-symmetric case.
    //! @param diameter output vector with only positive entries. */
    virtual void CalcDiameter( Vector<Double>& diameter ) = 0;

    //! Calculates the barycenter of the element. For 3D and 2D (z=0 then)
    //! This is the average of every coordinate.
    //! @param barycenter the output. In 2D z=0. If you want 2D then overload */
    virtual void CalcBarycenter( Vector<Double>& baryCenter )  = 0;

    //! Compute length of edge with maximal/minimal size
    virtual void GetMaxMinEdgeLength( Double& max, Double& min)  = 0;

    //! Compute the average side length for all dimenstions. Clearly works only 
    //! for quadrilaterals ans hexahedrons.
    //! @param edges_out array with index 0 for x, ... */
    virtual void GetEdgeLength( StdVector<Double>& edges_out) = 0;

    // ---------------------------------------------------
    //   Jacobian 
    // ---------------------------------------------------

    //! Calculate Jacobian Matrix
    /*! @param jac output Jacobian Matrix
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

    //! @see ElemShapeMap::SetElem
    virtual void SetElem( const Elem* ptElem, bool isUpdated = false );

    // ---------------------------------------------------
    //   Coordinate Mapping
    // ---------------------------------------------------

    //! @see ElemShapeMap::Local2Global
    void Local2Global( Vector<Double>& globPoint, 
                       const LocPoint& lp );

    //! @see ElemShapeMap::Global2Local
    void Global2Local( Vector<Double>& locPoint, 
                       const Vector<Double>& glob );
    
    //! @see ElemShapeMap::GetGlobMidPoint
    void GetGlodMidPoint( Vector<Double>& midPoint );

    //! @see ElemShapeMap::CalcVolume
    Double CalcVolume( );

    //! @see ElemShapeMap::CalcNormal
    void CalcNormal( Vector<Double>& normal, 
                     const LocPoint& lp );

    //! @see ElemShapeMap::CoordIsInsideElem
    bool CoordIsInsideElem( const StdVector<LocPoint>& lps,
                            const Double tolerance,
                            StdVector<bool>& coordsInside );

    //! @see ElemShapeMap::CalcDiameter
    void CalcDiameter( Vector<Double>& diameter );

    //! @see ElemShapeMap::CalcBarycenter
    void CalcBarycenter( Vector<Double>& baryCenter );

    //! @see ElemShapeMap::GetMaxMinEdgeLength
    void GetMaxMinEdgeLength( Double& max, Double& min);

    //! @see ElemShapeMap::GetEdgeLength
    void GetEdgeLength( StdVector<Double>& edges_out);

    // ---------------------------------------------------
    //   Jacobian 
    // ---------------------------------------------------

    //! @see ElemShapeMap::CalcJ
    void CalcJ( Matrix<Double>& jac, 
                const LocPoint& ip );
    
    //! @see ElemShapeMap::CalcJ
    Double CalcJDet( Matrix<Double>& jac, 
                     const LocPoint& ip );

  protected:

    //! Pointer to H1 element of lower order
    FeH1* ptFe_;
    
    //! Map with pointers to reference elements
    std::map<Elem::FEType, FeH1 *> feMap_;

    //! Nodal coordinates
    Matrix<Double> coords_;
    
    //! Shape of reference element
    ElemShape shape_;
    
  };

}

#endif /* ELEMSHAPEMAP_HH_ */
