/*
 * MortarInterface.hh
 *
 *  Created on: 23.02.2013
 *      Author: jens
 */

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
    
    RegionIdType GetMasterSurfRegion() const { return masterSurfRegion_; }
    
    RegionIdType GetSlaveSurfRegion() const { return slaveSurfRegion_; }
    
    RegionIdType GetMasterVolRegion() const { return masterVolRegion_; }
    
    RegionIdType GetSlaveVolRegion() const { return slaveVolRegion_; }
    
    bool IsPlanar() const { return isCoplanar_; }
    
    bool IsEulerian() const { return isEulerian_; }

    bool NeedsUpdate() const { return isMoving_; }
    
    const std::string& GetCoordSys() const { return coordSysId_; }
    
    PtrCoefFct GetGridVelocity() const { return gridVelo_; };

    void ResetInterface();

    void UpdateInterface();
    
  protected:
    
    void SetRotation(const std::string &coordSysId,
                     Double rpm,
                     const std::string &connectedRegions = "");
    
    void SetMotion( const StdVector<std::string> &offsetExpr,
                    const std::string &coordSysId = "default" );
    
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

    //! intersect two elements of arbitrary type
    bool IntersectPolygons( SurfElem *ifElem1, SurfElem *ifElem2,
                            StdVector<UInt> &newNodes, 
                            StdVector< Vector<Double> >& p1, StdVector< Vector<Double> >& p2, StdVector< Vector<Double> >& r,    
                            Vector<Double>& temp1, Vector<Double>& temp2, 
                            Vector<Double>& ez, Vector<Double>& v, Vector<Double>& e,
                            Matrix<Double>& rMat, Matrix<Double>& rMatTrans,
                            StdVector< Vector<Double> >& p1Rot, StdVector< Vector<Double> >& p2Rot);

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

    //! Calculates the intersection between two polygons
    //! \param p1 (in) first polygon to be intersected
    //! \param p2 (in) second polygon to be intersected
    //! \param r (out) resulting polygon
    //! \param c1 temporary vector for each thread
    //! \param c2 temporary vector for each thread
    //! \param e temporary vector for each thread
    //! \param temp1 temporary vector for each thread
    //! \param temp2 temporary vector for each thread
    bool CutPolys(StdVector< Vector<Double> > &p1,
      StdVector< Vector<Double> > &p2, const bool coPlanarIface,
      StdVector< Vector<Double> > &r,
      Vector<Double> & c1, Vector<Double> & c2, Vector<Double> & e,
      Vector<Double> & temp1, Vector<Double> & temp2);

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
    //! \param p (in) polygon to be triangulized
    //! \param tri (out) list of triangles
    UInt TriangulatePoly( const StdVector< Vector<Double> > &p,
                          StdVector<MortarNcSurfElem*> &tri,
                          StdVector<UInt> &newNodes );

    //! Adds a node to the grid and adds the node number to the new nodes vector
    //! \param coordinate is the position of the new node
    //! \param nodeNo is the node number the new node obtains from the grid
    //! \param newNodes is the new nodes vector
    void AddNodeToGrid(const Vector<Double>& coordinate, UInt& nodeNo, StdVector<UInt>& newNodes);


    // =======================================================================
    // Class variables
    // =======================================================================

    RegionIdType masterSurfRegion_;
    RegionIdType slaveSurfRegion_;
    RegionIdType masterVolRegion_;
    StdVector<RegionIdType> additionalVolRegions_;
    RegionIdType slaveVolRegion_;
    bool isCoplanar_;
    bool isEulerian_;
    bool isMoving_;
    bool moveMaster_;
    bool exportToGrid_;
    bool geoWarn_;
    std::string coordSysId_;
    CoordSystem* coordSys_;
    MathParser* mParser_;
    PtrCoefFct gridVelo_;
    StdVector<std::string> offsetExpr_;
    UInt mphOffset_[3];
    NcIntersectAlgo intersectAlgo_;
    Double tolAbs_;
    Double tolRel_;
    RegionIdType region_;
    bool mutualProjection_; // 2D or 3D projection of parallel surfaces
    bool cakePieceProjection_; // 2D projection of a wedge (cake piece)

    //caching for Line intersection
    Vector<Double> c0_Line_, c1_Line_, d0_Line_, d1_Line_, tmp_Line_;
    Vector<Double> diff0_Line_, diff1_Line_;
    Vector<Double> s_Line_, t_Line_;
    Vector<Double> normal_Line_;

    bool isReset_;
    
    //if true, nc interface will be updated, but no mesh update will be done
    //the mesh update will be done by iterative geometry updates (and there it will also be triggered)
    bool useMeshSmoothing_;

    bool precomputeIntersectionCandiates_;

    std::vector<Grid::ElemElemMatch> intersectionCandiatesIdx_;

    Vector<Double> translationVector_;
    // rotation angle for bloch periodic BC of a circle in radians
    Double cakePieceRotationAngle_;
    // center of rotation
    Vector<Double> cakePieceRotationCenter_;

    // caching of node numbers for rotating interfaces
    bool hasNodeNums_;
    StdVector<UInt> nodeNums_;

    // caching of moved interface coordinated for rotating interfaces
    bool hasMoveInterfaceCoords_;
    StdVector< Vector<Double> > moveInterfaceLocalCoords_;
    StdVector< Vector<Double> > moveInterfaceOrigCoords_;

    // caching of elements for rotating interface
    bool hasSurfElems_;
    StdVector<SurfElem*> masterElems_;
    StdVector<SurfElem*> slaveElems_;

#ifdef USE_OPENMP
    // lock for the grid for parallel evaluation of intersection faces
    omp_lock_t gridLock_;
    // lock for the new nodes vector for parallel evaluation of intersection faces
    omp_lock_t newNodesLock_;
    // lock for the new element vector for parallel evaluation of intersection faces
    omp_lock_t elemListLock_;
#endif


#ifdef USE_CGAL
    //! Calculates the intersection between two polygons using CGAL
    //! \param p1 (in) first polygon to be intersected
    //! \param p2 (in) second polygon to be intersected
    //! \param r (out) resulting polygon
    //! \param temp1 temporary vector for each thread
    //! \param temp2 temporary vector for each thread
    //! \param ez temporary vector for each thread
    //! \param v temporary vector for each thread
    //! \param e temporary vector for each thread
    //! \param rMat temporary matrix for each thread
    //! \param rMatTrans temporary matrix for each thread
    //! \param p1rot temporary vector for each thread
    //! \param p2rot temporary vector for each thread
    bool CutPolysCGAL(StdVector<Vector<Double> > &p1, StdVector<Vector<Double> > &p2,
                      const bool coPlanarIface, StdVector<Vector<Double> > &r, 
                      Vector<Double>& temp1, Vector<Double>& temp2,
                      Vector<Double>& ez, Vector<Double>& v, Vector<Double>& e,
                      Matrix<Double>& rMat, Matrix<Double>& rMatTrans,
                      StdVector< Vector<Double> >& p1Rot, StdVector< Vector<Double> >& p2Rot);


    void PreComputeIntersectionCandidatesCGAL(const StdVector<SurfElem*>& masterElems,const StdVector<SurfElem*>& slaveElems);

    //! Map for each dimension (key) a list containing the "boxes" of elements (value)
    StdVector<UInt> uniqueIdxMaster_;
    StdVector<UInt> uniqueIdxSlave_;
    std::vector<Grid::HandleBox> slaveBoxes_;
    std::vector<Grid::HandleBox> masterBoxes_;
#endif
};

} /* namespace CoupledField */

#endif /* _MORTARINTERFACE_HH_ */

