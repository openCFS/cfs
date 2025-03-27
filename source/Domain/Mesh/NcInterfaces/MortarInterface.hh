// ================================================================================================
/*!
 *       \file     MortarInterface.cc
 *       \brief    This class implements all grid operations for a non-conforming interface for 
 *                 the coupling of two regions. The interface can be coplanar or non-coplanar, 
 *                 as well as moving or stationary. In any case, intersection elements are computed.
 *                 For moving interfaces, the intersections are re-computed at every timestep. 
 *                 In the non-coplanar case, intersection elements are computed after orthogonal
 *                 projection of the primary surface onto the secondary surface. 
 *                 The class was created by "jens" in 2013. Edit by "pheidegger" in 2024. 
 *       \date     Mar 2024
 *       \author   pheidegger
 */
//================================================================================================

#ifndef _MORTARINTERFACE_HH_
#define _MORTARINTERFACE_HH_

#include "BaseNcInterface.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#ifdef USE_OPENMP
#include <omp.h>
#endif

namespace CoupledField {

// forward declarations
class CoordSystem;
class MathParser;
class NcSurfElemList;
struct MortarNcSurfElem;
struct SurfElem;
template<class TYPE> class StdVector;
template<class TYPE> class Vector;

class MortarInterface : public BaseNcInterface {
    
  public:
    
    MortarInterface(Grid* grid, PtrParamNode nciNode);

    virtual ~MortarInterface();
    
    // getters for surface regions
    RegionIdType GetPrimarySurfRegion() const { return primarySurfRegion_; }
    RegionIdType GetSecondarySurfRegion() const { return secondarySurfRegion_; }
    // getters for volume regions
    RegionIdType GetPrimaryVolRegion() const { return primaryVolRegion_; }
    RegionIdType GetSecondaryVolRegion() const { return secondaryVolRegion_; }
    RegionIdType GetMovingVolRegion() const { return movingVolRegion_; }
    StdVector<RegionIdType> GetConnectedVolRegions() const { return connectedVolRegions_; }

    //! return if the interface is a coplanar interface
    bool IsCoplanar() const { return isCoplanar_; }
    //! return if the interface should assign counter-rotating grid velocities
    bool IsEulerian() const { return isEulerian_; }
    // return if the interface is a moving interface
    bool IsMoving() const override { return isMoving_; }

    // return the name (id) of the local coordinate system
    const std::string& GetCoordSysName() const { return coordSysName_; }
    
    //! return the coefFunction describing the cylindrical grid velocity.
    //! used in the singlePDE when useEulerian is activated
    PtrCoefFct GetGridVelocity() const { return gridVelo_; };

    //! Resets the interface.
    //! 1) Clears the interface surface region of
    //! the grid. This also deletes the intersection elements; 
    //! 2) Deletes the corresponding named nodes of the grid, 
    //! i.e., the nodes of the intersection elements.
    //! Checks if the interface needs to be reset.
    //! ATTENTION: GridCFS::DeleteNamedNodes() will lead to problems if
    //! other nodes / elements were inserted after the interface that is 
    //! being reset, as this would change the element connectivity.
    //! Thus, moving interfaces must always be created after
    //! static ones. This is ensured in Grid::InitNcInterfacesFromXML().
    //! Also, multiple moving interfaces must be moved all together (in every timestep).
    void ResetInterface() override;

    //! Initializes the nci. Is repeatedly called by the TransientDriver 
    //! via Grid::MoveNcInterfaces() case of motion
    void UpdateInterface() override;
    
  private:

    //! Defines the volume and surf regions belonging to the interface:
    //! primarySurfRegion_, secondarySurfRegion_
    //! primaryVolRegion_, secondaryVolRegion_
    //! \param nciNode (in) the param node defining the non conforming interface
    void SetRegions(const PtrParamNode nciNode);

    //! Defines the intersection algorithm: intersectAlgo_;
    //! \param nciNode (in) the param node defining the non conforming interface
    void SetIntersectionType(const PtrParamNode nciNode);

    //! Set tolerances for line intersection algorithms and the minimum element size below which elements/nodes are discarded
    void SetTolerances(const PtrParamNode nciNode);

    //! Set the interface coordiates to forced values. 
    //! This should be used carefully as the coordinates of the actual grid are 
    //! changed here. However, this can speed up computation.
    //! \param nciNode (in) the param node defining the non conforming interface
    void SetForcedCoords(const PtrParamNode nciNode);

    //! kirill: translational p.b.c.
    //! A common interface for the primary and the secondary surfaces is for now created throug a parallel projection
    //! of the primary onto the secondary. Namely, we have to find a vector along which the primary is translated from
    //! the secondary. The subtraction of this vector from the primary grid nodes will give us the desired coordinates
    //! near the secondary side. These translated grid is to be used instead of the original primary grid during
    //! the grid intersection procedures
    //! \param nciNode (in) the param node defining the non conforming interface
    void SetMutualProjection(const PtrParamNode nciNode);

    //! Define motion settings. Checks if translational, rotational motion, or mesh smoothing is specified.
    //! If yes calls the needed functions to set the motion to the moving region nodes.
    //! \param nciNode (in) the param node defining the non conforming interface
    void SetMotion(const PtrParamNode nciNode);

    //! Defines the math parser offset expressions to define rotational motion.
    //! \param nciNode (in) the param node defining the non conforming interface
    void SetRotation(const PtrParamNode motionNode);

    //! Defines the math parser offset expressions to define translational motion.
    //! \param nciNode (in) the param node defining the non conforming interface
    void SetTranslation(const PtrParamNode motionNode);

    //! Assigns the motion math parser expression to the math parser
    //! \param offsetExpression (in) math parser expression that describes motion of the interface in every direction.
    void SetMotionMathParser(const StdVector<std::string>& offsetExpression);

    //! Check if connected regions are specified and set them in case
    //! \param motionNode (in) the param node describing the motion
    void SetConnectedRegions(const PtrParamNode motionNode);

    //! Check if the surface is coplanar and set isCoplanar_ accordingly
    //! \param nciNode (in) the param node defining the non conforming interface
    void SetCoplanar(const PtrParamNode nciNode);

    //! Checks the type of every primary / secondary element and determines if the meshes are valid for the given
    //! configurations. Is only called when geometry warnings are activated via geometryWarnings="yes" (=default)
    void CheckMeshValidity();

    //! Export the initial intersection grid. If isMoving_ is specified, this function exports
    //! the intersection grid of every timestep. 
    //! This is done by adding additional copies of the intersection elements to the grid.
    void ExportToGrid();

    //! Moves the interface by computing the necessary node offsets and assigning them
    //! to the grid. For non-default coordinate systems, offsets are computed by 
    //! retrieving the local offsets from the specified coordinate system and 
    //! assigning them to the local coordinate of the node. Eventually, the 
    //! shifted local coordinate is transformed to the global coordinate.
    //! The offsets are assigned to the grid via GridCFS::SetNodeOffset().
    //! Checks if the mesh needs to be moved before computation.
    void MoveInterface();

    // =======================================================================
    // Non-matching grid interface calculation
    // =======================================================================

    enum NcIntersectAlgo {
      NCI_INTERSECT_NONE,
      NCI_INTERSECT_LINE,
      NCI_INTERSECT_RECT,
      NCI_INTERSECT_POLYGON
    };


    //! return codes of function CutLines
    enum LineIntersectType {
      INTERSECT_NONE,
      INTERSECT_OUTSIDE,
      INTERSECT_ON_LINE2,
      INTERSECT_CROSS,
      INTERSECT_IN_A,
      INTERSECT_IN_B,
      INTERSECT_IN_C,
      INTERSECT_IN_D,
      INTERSECT_A_EQ_C,
      INTERSECT_A_AND_C
    };


    //! intersect two line elements
    bool IntersectLines( SurfElem *ifaceElem1, SurfElem *ifaceElem2,
                         StdVector<UInt> &newNodes );


    //! intersect two axiparallel quads
    bool IntersectRects( SurfElem *ifaceElem1, SurfElem *ifaceElem2,
                         StdVector<UInt> &newNodes );


    //! get interface element coordinates for polygon intersection
    void GetInterfaceElemCoordinates(SurfElem *ifElem, StdVector< Vector<Double> >& coordinates);


    //! Calculates the intersection of two lines [a,b] and [c,d] and stores
    //! the intersection point (if any) in e.
    //! return codes:
    //! INTERSECT_NONE: no intersection
    //! INTERSECT_OUTSIDE: lines [a,inf) and [c,inf) intersect outside [c,d]
    //! INTERSECT_ON_LINE2: lines [a,inf) and [c,d] intersect on [c,d]
    //! INTERSECT_CROSS: lines (a,b) and (c,d) intersect in e
    //! INTERSECT_IN_A: intersection in a
    //! INTERSECT_IN_B: intersection in b
    //! INTERSECT_IN_C: intersection in c
    //! INTERSECT_IN_D: intersection in d
    //! INTERSECT_A_EQ_C: a equals c
    //! \param a (in) starting point of line 1
    //! \param b (in) end point of line 1
    //! \param c (in) starting point of line 2
    //! \param d (in) end point of line 2
    //! \param e (out) intersection point
    LineIntersectType CutLines(const Vector<Double> &a,
        const Vector<Double> &b, const Vector<Double> &c,
        const Vector<Double> &d, Vector<Double> &e) const;


    //! Returns if a point lies inside a convex polygon. Optionally the
    //! centroid of the polygon may be given in order to save computational
    //! time.
    //! \param p (in) coordinates of the point of interest
    //! \param poly (in) polygon to be tested
    //! \param c (in) pointer to the centroid of the polygon (may be NULL)
    bool PointInsidePoly(const Vector<Double> &p,
      const StdVector< Vector<Double> > &poly,
      const Vector<Double> *const c);


    //! Computes the centroid of a polygon.
    //! \param p (in) polygon for which the centroid is to be computed.
    //! \param c (out) coordinates of the centroid
    void PolyCentroid(const StdVector< Vector<Double> > &p,
                        Vector<Double> &c);


    //! Returns the radius of the surrounding circle.
    //! \param p (in) polygon for which the centroid is to be computed.
    //! \param c (in) coordinates of the centroid
    Double PolyCircumcircle(const StdVector< Vector<Double> > &p,
                            const Vector<Double> &c);


    //! Computes a triangulation for a polygon.
    //! 
    //! \param inputPolygon (in) polygon to be triangulized
    //! \param outputTriangles (out) list of triangles
    //! \param newNodeIds (out) node Ids of the generated nodes
    void TriangulatePoly(const StdVector<Vector<Double>>& inputPolygon,
                         StdVector<MortarNcSurfElem*>& outputTriangles,
                         StdVector<UInt>& newNodeIds);


    //! Adds a node to the grid and adds the node number to the new nodes vector
    //! Inherently calls the Grid::AddNode() but ensures openmp locking.
    //! \param nodeCoordinate (in) the coordinates of the new node
    //! \param newNodeId (out) is the node number the new node obtained from the grid
    //! \param newNodeIds (in,out) vector containing multiple new nodes.
    void AddNodeToGrid(const Vector<Double>& nodeCoordinate, UInt& newNodeId, StdVector<UInt>& newNodeIds);
    

    //! Computes the possible candidates for element intersections. In default build this just every possible permutation 
    //! of surface elements.
    //! When using CGAL, use bounding boxes to check for intersection candidates to speed up computation.
    //! Allows to get rid of impossible combinations beforehand. 
    //! \param primarySurfElems (in) surface elements of the primary side
    //! \param secondarySurfElems (in) surface elements of the secondary side
    void ComputeIntersectionCandidates();


    //! Parses the intersectAlgo_ and isCoplanar_ and calles the correct type of intersection.
    //! Returns the nodes of the new intersection-element nodes
    void ComputeIntersection(StdVector<UInt>& newNodeIds);

    //! Projects the a polygon onto the plane of a donor polygon. Distinguishes between triangle and quadrilateral elements
    //! \param projectionPolygonCoords (in,out) takes the coordinates of the element to be projected and outputs the projected coordinates
    //! \param donorPolygonCoords (in) the coordinates of the polygon that donates the orientation
    //! \param projectionNormVec (out) the normal vector of the projection plane
    //! \param tempVec1, tempVec2 () helper vectors for mathematical operations
    void ProjectPolygon(StdVector<Vector<Double>>& projectionPolygonCoords, StdVector<Vector<Double>>& donorPolygonCoords,
                        Vector<Double>& projectionNormVec, Vector<Double>& tempVec1, Vector<Double>& tempVec2);


    //! Computes the rotation matrix to rotate elements into the (x,y,0) plane to use CGAL's 2D intersection methods
    //! Only computes the matrix for the first element if the interface isCoplanar_
    //! \param rotationMatrix (out) the computed rotation matrix
    //! \param rotationMatrixTrans (out) the transposed rotation matrix for the inverse rotation
    //! \param normVec (in) the normal vector of the element with ||normVec||=1
    //! \param rotationAxis (in) the rotation axis vector computed as cross(normVec, ez)/||cross(normVec, ez)||
    //! \param dirVec (in) direction vector into which the rotated normal should point
    bool GetRotationMatrix(Matrix<Double>& rotationMatrix, Matrix<Double>& rotationMatrixTrans,
                           Vector<Double>& normVec, Vector<Double>& rotationAxis, Vector<Double>& dirVec);

    //! =======================================================================
    //! Class variables
    //! =======================================================================
    //! surface region ids of both sides
    RegionIdType primarySurfRegion_;
    RegionIdType secondarySurfRegion_;
    //! corresponding volume region ids of both sides
    RegionIdType primaryVolRegion_;
    RegionIdType secondaryVolRegion_;
    //! specify the moving volume region. Is set to -1 by default (if no moving NCI is specified)
    RegionIdType movingVolRegion_;
    //! the surface region assigned to the non conforming interface
    RegionIdType nciRegion_;
    //! Store the number of calls of UpdateInterface().
    UInt callCounter_;
    //! additional volume regions that are statically connected to the moving region
    StdVector<RegionIdType> connectedVolRegions_;
    //! vector holding pairs of node ids for possible element intersections
    std::vector<std::pair<UInt,UInt>> intersectionCandiatesIdx_;
    //! caching of node ids and coordinates for moving interfaces
    StdVector<UInt> movingNodeIds_;                     //! ids of moving nodes on the (original) grid
    StdVector<Vector<Double>> movingNodeLocalCoords_;   //! original local coordinates of the moving nodes
    StdVector<Vector<Double>> movingNodeCoords_;        //! original global coordinates of the moving nodes
    //! caching of elements for rotating interface
    StdVector<SurfElem*> primarySurfElems_;    //! elements of the primary side
    StdVector<SurfElem*> secondarySurfElems_;  //! elements of the secondary side

    //! specified id of the coordinate system
    std::string coordSysName_;
    //! pointer to the local coordinate system corresponding to coordSysName_
    CoordSystem* coordSys_;

    //! pointer to the domain's mathParser
    MathParser* mParser_;
    //! pointer to the coefFunction describing the moving grid velocity
    PtrCoefFct gridVelo_;
    //! mathParser handles for the offset (motion) expressions
    UInt offsetMpHandles_[3];
    //! Type of the used intersection algorithm
    NcIntersectAlgo intersectAlgo_;
    //! absolute tolerance for the intersection algorithms.
    //! used to check for duplicate cuts in CutPolygon()
    //! used to check PointInsidePoly()
    //! used to check for line intersections in CutLines()
    Double absoluteTolerance_;
    //! relative tolerance used as an additional tolerance in CutLines()
    Double relativeTolerance_;
    //! minimum distance below which intersection nodes are discarded
    Double coincidenceRadius_;
    //! minimum relative side length of the intersection polygon, relative to the size of the smaller interface element.
    //! Below the threshold, which intersection elements are discarded.
    Double minRelativeSideLength_;
    //! used to create a bounding box of elements in ComputeIntersectionCandidates()
    Double boundingBoxTolerance_;

    //! Translation vector is only set for the mutualProjection method.
    //! Is 1/2 of the difference vector between centers of bounding boxes of both sides,
    //! i.e., the centers of the two coplanar surfaces
    Vector<Double> translationVector_;
    //! caching for Line intersection
    Vector<Double> c0_Line_, c1_Line_, d0_Line_, d1_Line_, tmp_Line_;
    Vector<Double> diff0_Line_, diff1_Line_;
    Vector<Double> s_Line_, t_Line_;
    Vector<Double> normal_Line_;

    //! specify if the surface is coplanar. In this case we can used a simplified 
    //! computation
    bool isCoplanar_;
    //! isEulerian_ can be optionally set for cylindrically moving interfaces.
    //! If specified, the grid velocity is computed by the gridVelo_ function and
    //! set against the grid motion. Hence, the physical field is not rotated.
    bool isEulerian_;
    //! specify if the NCI is moving
    bool isMoving_;
    //! export (initial) intersection grid to output grid
    bool exportToGrid_;
    //! print a warning if an intersection element is too small
    bool geoWarn_;
    //! Check if the surface has been reset after movement
    bool isReset_;
    //! flag if mutualProjection is used. The XSD scheme describes the mutualProjection method as:
    //! "Used for creating an intersection grid handling periodic boundary conditions.
    //! In this case, the primary and the secondary surfaces don't lie on the same interface. Therefore,
    //! there will be no intersection grid found unless they are projected on a common interface."
    bool mutualProjection_;
    //! If true, nc interface will be updated, but no mesh update will be done
    //! the mesh update will be done by iterative geometry updates (and there it will also be triggered)
    //! cannot be used in combination with a moving or rotating interface
    bool useMeshSmoothing_;


  #ifdef USE_OPENMP
    // lock for the grid for parallel evaluation of intersection faces
    omp_lock_t gridLock_;
    // lock for the new nodes vector for parallel evaluation of intersection faces
    omp_lock_t newNodesLock_;
    // lock for the new element vector for parallel evaluation of intersection faces
    omp_lock_t elemListLock_;
  #endif

  #ifdef USE_CGAL
    //! Compute the intersection grid for non-coplanar interfaces using CGAL
    //! \param newNodeIds (out) vector containing the node ids for the intersection grid
    void IntersectCoplanarPolygons(StdVector<UInt>& newNodeIds);


    //! Compute the intersection grid for non-coplanar interfaces using CGAL
    //! \param newNodeIds (out) vector containing the node ids for the intersection grid
    void IntersectGeneralPolygons(StdVector<UInt>& newNodeIds);


    //! Compute the intersection of two polygons using CGAL. The polygons are rotated first so that they are oriented along the (x,y) plane.
    //! The temporary vectors are initialized externally so they are not required to be re-initialized in every loop.
    //! \param projectedPrimaryNodeCoords     (in) coordinates of the primary polygon
    //! \param secondaryNodeCoords            (in) coordinates of the secondary polygon
    //! \param intersectionPolygonCoords      (out) coordinates of the intersection polygon
    //! \param tempVec1                       (in) helper vector for mathematical operations
    //! \param tempVec2                       (in) helper vector for mathematical operations
    //! \param rotationAxis                   (in) the rotation axis vector computed as cross(normVec, ez)/||cross(normVec, ez)||
    //! \param rotationMatrix                 (in) the rotation matrix to rotate the polygons into the (x,y) plane
    //! \param rotationMatrixTrans            (in) the transposed rotation matrix for the inverse rotation
    //! \param rotatedPrimaryNodeCoords       (in) coordinates of the rotated primary polygon
    //! \param rotatedSecondaryNodeCoords     (in) coordinates of the rotated secondary polygon
    //! \param rotatePolys                    (in) boolean to state if the polygons should be rotated
    bool CGAL_Cut2DPolygons(StdVector<Vector<Double>>& projectedPrimaryNodeCoords,
                            StdVector<Vector<Double>>& secondaryNodeCoords,
                            StdVector<Vector<Double>>& intersectionPolygonCoords,
                            Vector<Double>& tempVec1,
                            Vector<Double>& tempVec2,
                            Vector<Double>& rotationAxis,
                            Matrix<Double>& rotationMatrix,
                            Matrix<Double>& rotationMatrixTrans,
                            StdVector<Vector<Double>>& rotatedPrimaryNodeCoords,
                            StdVector<Vector<Double>>& rotatedSecondaryNodeCoords,
                            const bool& rotatePolys);


    //! List containing the indexes of the surface-element bounding boxes. Needed to precompute intersection candidates
    StdVector<UInt> primaryBoxIndexes_;
    StdVector<UInt> secondaryBoxIndexes_;


    //! Vectors storing the handles to the bounding boxes of the surface elements of the NCI surfaces
    std::vector<Grid::HandleBox> primarySurfElemBoxes_;
    std::vector<Grid::HandleBox> secondarySurfElemBoxes_;
  #else // USE_CGAL
    //! Computes intersection elements between one primary (primary) and one secondary (secondary) element
    //! of an interface and adds them to the grid.
    //! Internally calls functions CutPolygon() or CutPolygon() for computing the actual element intersections
    //! The temporary vectors are initialized externally so they are not required to be re-initialized in every loop.
    //! \param primaryElement             (in) pointer to an element on the primary side of the non-conforming interface
    //! \param secondaryElement           (in) pointer to an element on the secondary side of the non-conforming interface
    //! \param newNodeIds                 (in, out) reference to a StdVector holding all Ids of newly generated nodes
    //! \param primaryNodeCoords          (in) reference to all node coordinates of the primary element
    //! \param secondaryNodeCoords        (in) reference to all node coordinates of the secondary element
    //! \param intersectionPolygonCoords  (in, out) reference to all computed intersection-polygon node coordinates
    //! \param tempVec1                   (in) reference to memory used for multithreading purpose. Gets passed to CutPolygon()
    //! \param tempVec2                   (in) reference to memory used for multithreading purpose. Gets passed to CutPolygon()
    //! \param tempVec3                   (in) reference to memory used for multithreading purpose. Gets passed to CutPolygon()
    //! \param tempVec4                   (in) reference to memory used for multithreading purpose. Gets passed to CutPolygon()
    //! \param tempVec5                   (in) reference to memory used for multithreading purpose. Gets passed to CutPolygon()
    //! \param tempRotMat                 (in) reference to memory used for multithreading purpose. Gets passed to CutPolygon()
    //!                                        Stores the roatation matrix.
    //! \param tempRotMatTrans            (in) reference to memory used for multithreading purpose. Gets passed to CutPolygon()
    //!                                        Stores the transposed roatation matrix.
    //! \param tempRotPolyPrim            (in) reference to memory used for multithreading purpose. Gets passed to CutPolygon()
    //!                                        Stores the coordinates of the rotated polygon of the primary side.
    //! \param tempRotPolySec             (in) reference to memory used for multithreading purpose. Gets passed to CutPolygon()
    //!                                        Stores the coordinates of the rotated polygon of the secondary side.
    bool IntersectPolygons(SurfElem *primaryElement,
                           SurfElem *secondaryElement,
                           StdVector<UInt> &newNodeIds,
                           StdVector<Vector<Double>>& primaryNodeCoords,
                           StdVector<Vector<Double>>& secondaryNodeCoords,
                           StdVector<Vector<Double>>& intersectionPolygonCoords,
                           Vector<Double>& tempVec1,
                           Vector<Double>& tempVec2,
                           Vector<Double>& tempVec3,
                           Vector<Double>& tempVec4,
                           Vector<Double>& tempVec5,
                           Matrix<Double>& rotMat,
                           Matrix<Double>& rotMatTrans,
                           StdVector<Vector<Double>>& rotPolyPrim,
                           StdVector<Vector<Double>>& rotPolySec);


    //! Calculates the intersection between two polygons.
    //! The temporary vectors are initialized externally so they are not required to be re-initialized in every loop.
    //! \param primaryNodeCoords          (in) reference to all node coordinates of the primary element
    //! \param secondaryNodeCoords        (in) reference to all node coordinates of the secondary element
    //! \param isCoplanar                 (in) boolean to state if the interface is coplanar
    //! \param intersectionPolygonCoords  (in, out) reference to all computed intersection-polygon node coordinates
    //! \param primaryCentroid            (in) temporary vector storing the centroid of the primary polygon.
    //! \param secondaryCentroid          (in) temporary vector storing the centroid of the secondary polygon.
    //! \param tempVec1                   (in) temporary vector for mixed purposes
    //! \param tempVec2                   (in) temporary vector for mixed purposes
    //! \param intersectionPointCoords    (in) temporary vector storing intersection coordinates of lines.
    bool CutPolygon(StdVector<Vector<Double>>& primaryNodeCoords,
                    StdVector<Vector<Double>>& secondaryNodeCoords,
                    StdVector<Vector<Double>>& intersectionPolygonCoords,
                    Vector<Double>& primaryCentroid,
                    Vector<Double>& secondaryCentroid,
                    Vector<Double>& tempVec1,
                    Vector<Double>& tempVec2,
                    Vector<Double>& intersectionPointCoords);
  #endif
};

} /* namespace CoupledField */

#endif /* _MORTARINTERFACE_HH_ */

