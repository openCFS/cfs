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
                     Double rpm);
    
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

    //! intersect two elements of arbitrary type
    bool IntersectPolygons( SurfElem *ifElem1, SurfElem *ifElem2,
                            StdVector<UInt> &newNodes );


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
    bool CutPolys(StdVector< Vector<Double> > &p1,
      StdVector< Vector<Double> > &p2, const bool coPlanarIface,
      StdVector< Vector<Double> > &r) const;

    //! Returns if a point lies inside a convex polygon. Optionally the
    //! centroid of the polygon may be given in order to save computational
    //! time.
    //! \param p (in) coordinates of the point of interest
    //! \param poly (in) polygon to be tested
    //! \param c (in) pointer to the centroid of the polygon (may be NULL)
    bool PointInsidePoly(const Vector<Double> &p,
      const StdVector< Vector<Double> > &poly,
      const Vector<Double> *const c) const;

    //! Computes the centroid of a polygon and returns the radius of the
    //! surrounding circle.
    //! \param p (in) polygon for which the centroid is to be computed.
    //! \param c (out) coordinates of the centroid
    Double PolyCentroid(const StdVector< Vector<Double> > &p,
      Vector<Double> &c) const;

    //! Computes a triangulation for a polygon.
    //! \param p (in) polygon to be triangulized
    //! \param tri (out) list of triangles
    UInt TriangulatePoly( const StdVector< Vector<Double> > &p,
                          StdVector<MortarNcSurfElem*> &tri,
                          StdVector<UInt> &newNodes );

    // =======================================================================
    // Class variables
    // =======================================================================

    RegionIdType masterSurfRegion_;
    RegionIdType slaveSurfRegion_;
    RegionIdType masterVolRegion_;
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

    bool isReset_;

};

} /* namespace CoupledField */

#endif /* _MORTARINTERFACE_HH_ */
