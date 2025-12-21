// ================================================================================================
/*!
 *       \file     MortarInterface.cc
 *       \brief    This class implements all grid operations for a non-conforming interface for 
 *                 the coupling of two regions. The interface can be coplanar or non-coplanar, 
 *                 as well as moving or stationary. In any case, intersection elements are computed.
 *                 For moving interfaces, the intersections are re-computed at every timestep. 
 *                 In the non-coplanar case, intersection elements are computed afer orthogonal
 *                 projection of the primary surface onto the secondary surface. 
 *                 The class was created by "jens" in 2013. Edit by "pheidegger" in 2024. 
 *       \date     Mar 2024
 *       \author   pheidegger, jens
 */
//================================================================================================

#include <def_use_cgal.hh>

#include "MortarInterface.hh"
#include "PolygonIterators.hh"
#include "Domain/Domain.hh"
#include "Domain/ElemMapping/SurfElem.hh"
#include "Domain/ElemMapping/EntityLists.hh"
#include "Domain/CoefFunction/CoefFunction.hh"
#include "Domain/CoordinateSystems/CoordSystem.hh"
#include "Driver/TransientDriver.hh"
#include "Utils/StdVector.hh"
#include "MatVec/Vector.hh"
#include "Utils/mathParser/mathParser.hh"

#include "Utils/Timer.hh"

#include <sstream>
#include <boost/shared_ptr.hpp>

#ifdef USE_CGAL

#pragma GCC diagnostic push
// for gcc7 and CGAL 4.9.1
#pragma GCC diagnostic ignored "-Wuninitialized"

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Boolean_set_operations_2.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef Kernel::Point_2 CGALPoint2;
typedef CGAL::Polygon_2<Kernel> CGALPolygon2;
typedef CGAL::Polygon_with_holes_2<Kernel> CGALPolygonWithHoles2;

#pragma GCC diagnostic pop

#endif

#include "DataInOut/Logging/LogConfigurator.hh"
DEFINE_LOG(mortarInterface, "mortarInterface")
namespace CoupledField {

  //=================================================
  // CGAL Presort of intersection candidates
  //=================================================

  #ifdef USE_CGAL
  // Iterator reporter class, returning the two ids of the CGAL-Boxes
  template <class OutputIterator>
  struct CGAL_ElemElemIdReporter {
    OutputIterator it;
    CGAL_ElemElemIdReporter(OutputIterator i  )
    : it(i) {} // store iterator in object

    // We write the id-number of box a to the output iterator assuming
    // that box b (the query box) is not interesting in the result.
    void operator()( const Grid::HandleBox& a, const Grid::HandleBox& b) {
      std::pair<UInt, UInt > pair;
      //ids seems to be one based
      pair.first = *a.handle();
      pair.second = *b.handle();
      *it++ = pair;
    }
  };
  // helper function to create the function object
  template <class Iter>
  CGAL_ElemElemIdReporter<Iter> elemElemIdReporter(Iter it)
  { return CGAL_ElemElemIdReporter<Iter>(it); }
  #endif
  //=================================================
  // END: CGAL Presort of intersection candidates
  //=================================================


  MortarInterface::MortarInterface(Grid* grid, PtrParamNode nciNode) :
    BaseNcInterface(grid),
    primarySurfRegion_(NO_REGION_ID),
    secondarySurfRegion_(NO_REGION_ID),
    primaryVolRegion_(NO_REGION_ID),
    secondaryVolRegion_(NO_REGION_ID),
    movingVolRegion_(NO_REGION_ID),
    nciRegion_(NO_REGION_ID),
    callCounter_(0),
    coordSysName_(""),
    coordSys_(nullptr),
    mParser_(nullptr),
    intersectAlgo_(NCI_INTERSECT_NONE),
    absoluteTolerance_(-1),
    relativeTolerance_(-1),
    isCoplanar_(false),
    isEulerian_(false),
    isMoving_(false),
    exportToGrid_(true),
    geoWarn_(true),
    isReset_(false),
    mutualProjection_(false),
    useMeshSmoothing_(false) {
  #ifdef USE_OPENMP
    // initialize locks to synchronize parallel threads
    omp_init_lock(&gridLock_);
    omp_init_lock(&newNodesLock_);
    omp_init_lock(&elemListLock_);
  #endif
    // get info from the xml...
    // name of the nci
    this->SetName(nciNode->Get("name")->As<std::string>());
    // set the name of the element list
    elemList_->SetName(name_);
    // state that nci is being set
    std::cout << std::endl << "Assigning non-conforming interface '" << name_ <<"'." << std::endl;
    // assign regions
    SetRegions(nciNode);
    // check for forced coords
    SetForcedCoords(nciNode);
    // check if intersection grid should be exported to the result
    nciNode->GetValue("storeIntegrationGrid", exportToGrid_, ParamNode::PASS);
    // check for warning outputs. Print warning if small elements are rejected
    nciNode->GetValue("geometryWarnings", geoWarn_, ParamNode::PASS);
    // set tolerances dependent on intersection method
    SetTolerances(nciNode);
    // check for mutual projection method
    SetMutualProjection(nciNode);
    // check for general motion
    SetMotion(nciNode);
    // check if interface is coplanar
    SetCoplanar(nciNode);
    // assign intersection type
    SetIntersectionType(nciNode);
    // check if we have valid meshes on both interface sides
    if (geoWarn_==true)
      CheckMeshValidity();
    // Finally, initialize the interface by calling UpdateInterface(). Needed for Static and moving meshes.
    // If there is motion, UpdateInterface will be called repeatedly by TransientDriver.
    UpdateInterface();
    // Output debug message that NCI was assigned
    LOG_DBG(mortarInterface) << "Assigned non-conforming interface with the following properties:" << std::endl
                             << "Name: " << name_ << std::endl
                             << "Primary surfregion: " << ptGrid_->GetRegionName(primarySurfRegion_) << std::endl
                             << "Secondary surf region: " << ptGrid_->GetRegionName(secondarySurfRegion_) << std::endl
                             << "Primary volume region: " << ptGrid_->GetRegionName(primaryVolRegion_) << std::endl
                             << "Secondary volume region: " << ptGrid_->GetRegionName(secondaryVolRegion_) << std::endl
                             << "Absolute intersection tolerance: " << absoluteTolerance_ << std::endl
                             << "Relative intersection tolerance: " << relativeTolerance_ << std::endl
                             << "Is moving: " << isMoving_ << std::endl
                             << "Is coplanar: " << isCoplanar_ << std::endl
                             << "Is Eulerian: " << isEulerian_ << std::endl
                             << "Do grid export: " << exportToGrid_ << std::endl
                             << "Print geometry warnings: " << geoWarn_ << std::endl
                             << "Use mutual projection method: " << mutualProjection_ << std::endl;
    if (isMoving_) {
      LOG_DBG(mortarInterface) << "Non-conforming interface is moving and has the following attributes:" << std::endl
                               << "Name of the coordinate system: " << coordSysName_ << std::endl
                               << "Moving volume region: " << ptGrid_->GetRegionName(movingVolRegion_) << std::endl
                               << "Has connected regions: " << (connectedVolRegions_.GetSize()>0) << std::endl
                               << "Updated geometry movement (mesh smoothing): " << useMeshSmoothing_ << std::endl;
    }
  }

  MortarInterface::~MortarInterface() {
  #ifdef USE_OPENMP
    // destroy locks
    omp_destroy_lock(&gridLock_);
    omp_destroy_lock(&newNodesLock_);
    omp_destroy_lock(&elemListLock_);
  #endif
    // reset pointers
    ptGrid_ = nullptr;
    coordSys_ = nullptr;
    elemList_ = nullptr;
    // reset math parser and handles
    if (mParser_) {
      for (UInt iDim = 0; iDim < 3; ++iDim ) {
        if (offsetMpHandles_[iDim] != MathParser::GLOB_HANDLER)
          mParser_->ReleaseHandle(offsetMpHandles_[iDim]);
      }
      mParser_ = nullptr;
    }
  }


  void MortarInterface::SetRegions(const PtrParamNode nciNode) {
    // get primary / secondary region ids
    primarySurfRegion_ = ptGrid_->GetRegion().Parse(nciNode->Get("primarySide")->As<std::string>(), NO_REGION_ID);
    secondarySurfRegion_ = ptGrid_->GetRegion().Parse(nciNode->Get("secondarySide")->As<std::string>(), NO_REGION_ID);
    if (primarySurfRegion_ == NO_REGION_ID)
      EXCEPTION("Cannot find primary region of ncInterface: '" << name_ << "'.")
    if (secondarySurfRegion_ == NO_REGION_ID)
      EXCEPTION("Cannot find secondary region of ncInterface '" << name_ << "'.");
    // get the interface surface elements from the grid
    ptGrid_->GetSurfElems(primarySurfElems_, primarySurfRegion_);
    ptGrid_->GetSurfElems(secondarySurfElems_, secondarySurfRegion_);
    if (primarySurfElems_.GetSize() == 0)
      EXCEPTION("There are no surface elements in the primary region of " << "ncInterface '" << name_ << "'.");
    if (secondarySurfElems_.GetSize() == 0)
      EXCEPTION("There are no surface elements in the secondary region of " << "ncInterface '" << name_ << "'.");
    // get corresponding volume ids
    primaryVolRegion_ = primarySurfElems_[0]->ptVolElems[0]->regionId;
    secondaryVolRegion_ = secondarySurfElems_[0]->ptVolElems[0]->regionId;
  }


  void MortarInterface::SetIntersectionType(const PtrParamNode nciNode) {
    // define the intersection algorithm
    std::string intersectionMethod;
    nciNode->GetValue("intersectionMethod", intersectionMethod, ParamNode::PASS);
    if (ptGrid_->GetDim() == 2) {
      intersectAlgo_ = NCI_INTERSECT_LINE;
      if ( !intersectionMethod.empty() ) {
        WARN("ncInterface '" << name_ << "': intersection algorithm '" << intersectionMethod << "' is ignored in a 2D simulation and always set to line intersectoins.");
      }
    } else { // 3D
      if (intersectionMethod == "coaxi") {
        WARN("You provided the intersection method 'coaxi' for the non-conforming interface '" << name_ << "'." << std::endl
               << "This is currently untested!" << std::endl
               << "Please consider adding a testcase to the testsuite." << std::endl);
        if (!isCoplanar_)
          EXCEPTION("Only coplanar interfaces are supported with coaxial rectangle algorithm.");
        intersectAlgo_ = NCI_INTERSECT_RECT;
      } else {
        intersectAlgo_ = NCI_INTERSECT_POLYGON;
      }
    }
  }


  void MortarInterface::SetTolerances(const PtrParamNode nciNode) {
    // get tolerances for intersection computation if specified
    nciNode->GetValue("tolAbs", absoluteTolerance_);
    nciNode->GetValue("tolRel", relativeTolerance_);
    nciNode->GetValue("minRelSideLength", minRelativeSideLength_);
    nciNode->GetValue("coincidenceRadius", coincidenceRadius_);
    nciNode->GetValue("tolBbox", boundingBoxTolerance_);
  }


  void MortarInterface::SetForcedCoords(const PtrParamNode nciNode) {
    // check if forced x/y/z coordinates are specified
    bool doForceX = (nciNode->Get("forceXValue")->As<std::string>()!="");
    bool doForceY = (nciNode->Get("forceYValue")->As<std::string>()!="");
    bool doForceZ = (nciNode->Get("forceZValue")->As<std::string>()!="");
    if(doForceX||doForceY||doForceZ) {
      Double fX=0,fY=0,fZ=0;
      StdVector<UInt> ifNodeList;
      Vector<Double> curCoord(ptGrid_->GetDim());
      Vector<Double> newCoord(ptGrid_->GetDim());
      // get forced coords
      if (doForceX)
        nciNode->GetValue("forceXValue", fX, ParamNode::PASS);
      if (doForceY)
        nciNode->GetValue("forceYValue", fY, ParamNode::PASS);
      if (doForceZ)
        nciNode->GetValue("forceZValue", fZ, ParamNode::PASS);
      // set the coords of primary side
      ptGrid_->GetNodesByRegion(ifNodeList,primarySurfRegion_);
      for(UInt iNodes = 0; iNodes < ifNodeList.GetSize(); ++iNodes) {
        ptGrid_->GetNodeCoordinate(curCoord,ifNodeList[iNodes],false);
        newCoord[0] = (doForceX)? fX : curCoord[0];
        newCoord[1] = (doForceY)? fY : curCoord[1];
        if(ptGrid_->GetDim() == 3)
          newCoord[2] = (doForceZ)? fZ : curCoord[2];
        ptGrid_->SetNodeCoordinate(ifNodeList[iNodes],newCoord);
      }
      // clear
      ifNodeList.Clear(true);
      // set the coords of secondary side
      ptGrid_->GetNodesByRegion(ifNodeList,secondarySurfRegion_);
      for(UInt iNodes = 0; iNodes < ifNodeList.GetSize(); ++iNodes) {
        ptGrid_->GetNodeCoordinate(curCoord,ifNodeList[iNodes],false);
        newCoord[0] = (doForceX)? fX : curCoord[0];
        newCoord[1] = (doForceY)? fY : curCoord[1];
        if (ptGrid_->GetDim() == 3)
          newCoord[2] = (doForceZ)? fZ : curCoord[2];
        ptGrid_->SetNodeCoordinate(ifNodeList[iNodes],newCoord);
      }
    }
  }

  void MortarInterface::SetMutualProjection(const PtrParamNode nciNode)
  {
    nciNode->GetValue("mutualProjection", mutualProjection_, ParamNode::PASS);
    nciNode->GetValue("cakePieceProjection", cakePieceProjection_, ParamNode::PASS);
    if (mutualProjection_ && cakePieceProjection_) {
      EXCEPTION("Both mutual and cake-piece projection methods are specified. Only one of them can be used at a time.");
    }

    if (mutualProjection_) {
      translationVector_.Resize(ptGrid_->GetDim());
      Matrix<Double> bboxPrim, bboxSec;
      ptGrid_->CalcBoundingBoxOfRegion(primarySurfRegion_, bboxPrim);
      ptGrid_->CalcBoundingBoxOfRegion(secondarySurfRegion_, bboxSec);
      // calculate the vector along which the primary side was translated from the secondary side
      // as the difference between the centres of primary's and secondary's bounding boxes
      for (UInt iDim = 0; iDim < ptGrid_->GetDim(); ++iDim)
        translationVector_[iDim] = bboxPrim[iDim][0] + bboxPrim[iDim][1] - bboxSec[iDim][0] - bboxSec[iDim][1];
      translationVector_ *= 0.5;
    }
    else if (cakePieceProjection_) {
      // consider 2D-rotated mutual projection (cake-piece projection)
      if (cakePieceProjection_) {
        if (ptGrid_->GetDim() != 2) {
          EXCEPTION("Rotated mutual projection is only supported for 2D grids");
        }
        else {
          Matrix<Double> bboxPrim, bboxSec;
          // get bounding boxes of lines [min(x),max(x);min(y),max(y)]
          ptGrid_->CalcBoundingBoxOfRegion(primarySurfRegion_, bboxPrim);
          ptGrid_->CalcBoundingBoxOfRegion(secondarySurfRegion_, bboxSec);
          // normalize tolerances
          Double tol = fmax(fmax(fabs(bboxPrim[0][0] - bboxPrim[0][1]),
                                 fabs(bboxPrim[1][0] - bboxPrim[1][1])),
                            fmax(fabs(bboxSec[0][0] - bboxSec[0][1]),
                                 fabs(bboxSec[1][0] - bboxSec[1][1]))) *
                       relativeTolerance_;

          // find intersection point of the bounding boxes
          if (bboxPrim[0][0] > bboxSec[0][1] + tol || bboxPrim[1][1] < bboxSec[1][0] - tol)
            EXCEPTION("No intersecting bounding boxes found for cake-piece projection!");
          if (fabs(bboxPrim[0][0] - bboxSec[0][0]) < tol &&
              fabs(bboxPrim[0][1] - bboxSec[0][1]) < tol &&
              fabs(bboxPrim[1][0] - bboxSec[1][0]) < tol &&
              fabs(bboxPrim[1][1] - bboxSec[1][1]) < tol)
            EXCEPTION("Bounding boxes for cake-piece projection entirely overlap. Intersection is not unique!");
          // check if both vectors have equal lengths
          Double lenPrim, lenSec;
          lenPrim = sqrt(pow(bboxPrim[0][1] - bboxPrim[0][0], 2) + pow(bboxPrim[1][1] - bboxPrim[1][0], 2));
          lenSec = sqrt(pow(bboxSec[0][1] - bboxSec[0][0], 2) + pow(bboxSec[1][1] - bboxSec[1][0], 2));
          if (fabs(lenPrim - lenSec) > tol)
            EXCEPTION("Rotated mutual projection is only supported for lines of equal lengths");

          // get two nodes of each line..
          StdVector<UInt> nodesPrim, nodesSec;
          StdVector<Vector<Double>> coordsPrim(2), coordsSec(2);
          ptGrid_->GetNodesByRegion(nodesPrim, primarySurfRegion_);
          ptGrid_->GetNodesByRegion(nodesSec, secondarySurfRegion_);
          for (UInt iNode = 0; iNode < 2; ++iNode) {
            ptGrid_->GetNodeCoordinate(coordsPrim[iNode], nodesPrim[iNode], false);
            ptGrid_->GetNodeCoordinate(coordsSec[iNode], nodesSec[iNode], false);
          }

          // get line parameters and find intersections
          Vector<Double> intersectionPoint(2);
          if (fabs(coordsPrim[1][0] - coordsPrim[0][0]) > tol &&
              fabs(coordsSec[1][0] - coordsSec[0][0])) {
            // if both lines can be parametrized as y = k*x + d, find the params
            Double kPrim, kSec, dPrim, dSec;
            kPrim = (coordsPrim[1][1] - coordsPrim[0][1]) / (coordsPrim[1][0] - coordsPrim[0][0]);
            kSec = (coordsSec[1][1] - coordsSec[0][1]) / (coordsSec[1][0] - coordsSec[0][0]);
            dPrim = coordsPrim[0][1] - kPrim * coordsPrim[0][0];
            dSec = coordsSec[0][1] - kSec * coordsSec[0][0];
            intersectionPoint[0] = (dSec - dPrim) / (kPrim - kSec);
            intersectionPoint[1] = kPrim * intersectionPoint[0] + dPrim;
          }
          else if (fabs(coordsPrim[1][0] - coordsPrim[0][0]) <= tol) {
            // if the primary line is vertical (would lead to infinite kPrim)
            // find the params of secondary side
            Double kSec, dSec;
            kSec = (coordsSec[1][1] - coordsSec[0][1]) / (coordsSec[1][0] - coordsSec[0][0]);
            dSec = coordsSec[0][1] - kSec * coordsSec[0][0];
            intersectionPoint[0] = (coordsPrim[1][0] + coordsPrim[0][0]) / 2;
            intersectionPoint[1] = kSec * intersectionPoint[0] + dSec;
          }
          else if (fabs(coordsSec[1][0] - coordsSec[0][0]) <= tol) {
            // if the secondary line is vertical (would lead to infinite kSec)
            // find the params of primary side
            Double kPrim, dPrim;
            kPrim = (coordsPrim[1][1] - coordsPrim[0][1]) / (coordsPrim[1][0] - coordsPrim[0][0]);
            dPrim = coordsPrim[0][1] - kPrim * coordsPrim[0][0];
            intersectionPoint[0] = (coordsSec[1][0] + coordsSec[0][0]) / 2;
            intersectionPoint[1] = kPrim * intersectionPoint[0] + dPrim;
          }
          else {
            EXCEPTION("Unexpected case.");
          }

          // check if intersection point is in bounding box. We need to extend the bounding box to include the intersection point
          if (intersectionPoint[0] < bboxPrim[0][0] - tol || intersectionPoint[0] > bboxPrim[0][1] + tol ||
              intersectionPoint[1] < bboxPrim[1][0] - tol || intersectionPoint[1] > bboxPrim[1][1] + tol) {
            WARN("Intersection point :" << intersectionPoint << "\nfor cake-piece projection is not in bounding box of primary surface:\n"
                                        << bboxPrim << "\nThe bounding box will be extended to include the intersection point.");
            bboxPrim[0][0] = fmin(bboxPrim[0][0], intersectionPoint[0]);
            bboxPrim[0][1] = fmax(bboxPrim[0][1], intersectionPoint[0]);
            bboxPrim[1][0] = fmin(bboxPrim[1][0], intersectionPoint[1]);
            bboxPrim[1][1] = fmax(bboxPrim[1][1], intersectionPoint[1]);
          }
          if (intersectionPoint[0] < bboxSec[0][0] - tol || intersectionPoint[0] > bboxSec[0][1] + tol ||
              intersectionPoint[1] < bboxSec[1][0] - tol || intersectionPoint[1] > bboxSec[1][1] + tol) {
            WARN("Intersection point :" << intersectionPoint << "\nfor cake-piece projection is not in bounding box of secondary surface:\n"
                                        << bboxSec << "\nThe bounding box will be extended to include the intersection point.");
            bboxSec[0][0] = fmin(bboxSec[0][0], intersectionPoint[0]);
            bboxSec[0][1] = fmax(bboxSec[0][1], intersectionPoint[0]);
            bboxSec[1][0] = fmin(bboxSec[1][0], intersectionPoint[1]);
            bboxSec[1][1] = fmax(bboxSec[1][1], intersectionPoint[1]);
          }

          // extract line vectors from bounding boxes, knowing the intersection point.
          Vector<Double> linePrim(2), lineSec(2);
          linePrim[0] = (fabs(bboxPrim[0][0] - intersectionPoint[0]) < tol) ? bboxPrim[0][1] - intersectionPoint[0] : bboxPrim[0][0] - intersectionPoint[0];
          linePrim[1] = (fabs(bboxPrim[1][0] - intersectionPoint[1]) < tol) ? bboxPrim[1][1] - intersectionPoint[1] : bboxPrim[1][0] - intersectionPoint[1];
          lineSec[0] = (fabs(bboxSec[0][0] - intersectionPoint[0]) < tol) ? bboxSec[0][1] - intersectionPoint[0] : bboxSec[0][0] - intersectionPoint[0];
          lineSec[1] = (fabs(bboxSec[1][0] - intersectionPoint[1]) < tol) ? bboxSec[1][1] - intersectionPoint[1] : bboxSec[1][0] - intersectionPoint[1];

          // calculate the angle between the lines
          cakePieceRotationAngle_ = acos((linePrim * lineSec) / (linePrim.NormL2() * lineSec.NormL2()));
          // get the direction of rotation via cross product
          Double crossProd;
          crossProd = linePrim[0] * lineSec[1] - linePrim[1] * lineSec[0];
          if (crossProd < 0)
            cakePieceRotationAngle_ = -cakePieceRotationAngle_;

          // store the center of the secondary bounding box as the center of rotation
          cakePieceRotationCenter_.Resize(ptGrid_->GetDim());
          for (UInt iDim = 0; iDim < ptGrid_->GetDim(); ++iDim)
            cakePieceRotationCenter_[iDim] = (bboxPrim[iDim][0] + bboxPrim[iDim][1]) / 2;

          // calculate the vector along which the master side was translated from the slave side
          // as the difference between the centres of master's and slave's bounding boxes
          translationVector_.Resize(ptGrid_->GetDim());
          for (UInt iDim = 0; iDim < ptGrid_->GetDim(); ++iDim)
            translationVector_[iDim] = bboxPrim[iDim][0] + bboxPrim[iDim][1] - bboxSec[iDim][0] - bboxSec[iDim][1];
          translationVector_ *= 0.5;
          // debug messages (convert to debug logs before merge)
          // std::cout << "bboxPrim: \n" << bboxPrim << std::endl;
          // std::cout << "bboxSec: \n" << bboxSec << std::endl;
          // std::cout << "lenPrim: \n" << lenPrim << std::endl;
          // std::cout << "lenSec: \n" << lenSec << std::endl;
          // std::cout << "coordsPrim: \n" << coordsPrim << std::endl;
          // std::cout << "coordsSec: \n" << coordsSec << std::endl;
          // std::cout << "Intersection point: \n" << intersectionPoint << std::endl;
          // std::cout << "linePrim: \n" << linePrim << std::endl;
          // std::cout << "lineSec: \n" << lineSec << std::endl;
          // std::cout << "crossProd" << crossProd << std::endl;
          // Double rotDeg = cakePieceRotationAngle_ * 180.0 / M_PI;
          // std::cout << "cakePieceRotationAngle: " << rotDeg << std::endl;
          // std::cout << "cakePieceRotationCenter_: \n" << cakePieceRotationCenter_ << std::endl;
          // std::cout << "translationVector: \n" << translationVector_ << std::endl;
        }
      }
    }
    else {
      translationVector_.Resize(0);
    }
  }

  void MortarInterface::SetMotion(const PtrParamNode nciNode) {
    bool isRotating = false;
    // check if rotational motion is specified
    PtrParamNode motionNode = nciNode->Get("rotation", ParamNode::PASS);
    if (motionNode)
      isRotating = true;
    else
      // check if translational motion is specified
      motionNode = nciNode->Get("generalMotion", ParamNode::PASS);

    if (motionNode) {
      // set the motion coordinate system
      // gather the specified coord system
      coordSysName_ = motionNode->Get("coordSysId")->As<std::string>();
      coordSys_ = domain->GetCoordSystem(coordSysName_);
      // specify the moving region
      if (motionNode->Get("movingSide", ParamNode::INSERT)->As<std::string>() == "primary")
        movingVolRegion_ = primaryVolRegion_;
      else if (motionNode->Get("movingSide", ParamNode::INSERT)->As<std::string>() == "secondary")
        movingVolRegion_ = secondaryVolRegion_;
      // get info if Eulerian system is used
      isEulerian_ = (motionNode->Get("ALESystem")->As<std::string>()=="yes");
      if (isEulerian_==true && isRotating==false) {
        WARN("You activated the ALESystem in a non-rotational domain. This is not supported yet." <<
             "ALESystem will be deactivated!")
        isEulerian_ = false;
      }
    } else {
      isEulerian_ = false;
    }

    // check if mesh smoothing is specified
    useMeshSmoothing_ = (nciNode->Get("passiveGeomUpdate")->As<std::string>()=="yes");
    if(useMeshSmoothing_ == true) {
      if (motionNode) {
        WARN("You activated passiveGeomUpdate. Mortar interface will not move actively" <<
             " but will be updated at every timestep. The specified motion will be ignored!");
      } else if (!motionNode) {
        EXCEPTION("You must specify which side is moved when passiveGeomUpdate is activated!");
      }
      isMoving_ = true;
    }

    // if no motion is defined, return
    if (motionNode == nullptr)
      return;

    // set the motion to the mathparser
    if (!isRotating) {
      SetTranslation(motionNode);
    } else {
      SetRotation(motionNode);
    }

    // initialize the moving nodes by assigning a 0.0 node offset to every node
    if (isMoving_) {
      Vector<Double> nullOffsets;
      // cache all nodes that are to be moved
      ptGrid_->GetNodesByRegion(movingNodeIds_, movingVolRegion_);
      // check for connected regions
      SetConnectedRegions(motionNode);
      // set 0 offset  vector
      nullOffsets.Resize(movingNodeIds_.GetSize() * coordSys_->GetDim(), 0.0);
      // set the offset
      ptGrid_->SetNodeOffset(movingNodeIds_, nullOffsets);
      // cache the global and local original node coordinates
      movingNodeCoords_.Resize(movingNodeIds_.GetSize());
      movingNodeLocalCoords_.Resize(movingNodeIds_.GetSize());
      #pragma omp parallel for num_threads(CFS_NUM_THREADS)
      for (Integer iNodes = 0; iNodes < (Integer)movingNodeIds_.GetSize(); ++iNodes) { // omp only allows us to use int here
        // get the node coordinates
        ptGrid_->GetNodeCoordinate(movingNodeCoords_[iNodes], movingNodeIds_[iNodes], false);
        // get the corresponding local coordinates
        coordSys_->Global2LocalCoord(movingNodeLocalCoords_[iNodes], movingNodeCoords_[iNodes]);
      }
    } else {
      WARN("You supplied constant expressions for displacements for moving ncInterface '" << name_ <<
           "'. The interface is assumed stationary.");
    }
  }


  void MortarInterface::SetTranslation(const PtrParamNode motionNode) {
    // set expression for translational motion...
    // extract the math parser expressions for the displacement functions
    StdVector<std::string> offsetExpression;
    if (motionNode->Has("displace3")) {
      offsetExpression.Resize(3);
      offsetExpression[2] = motionNode->Get("displace3")->As<std::string>();
    } else {
      offsetExpression.Resize(2);
    }
    offsetExpression[0] = motionNode->Get("displace1")->As<std::string>();
    offsetExpression[1] = motionNode->Get("displace2")->As<std::string>();
    if ( coordSys_->GetDim() != offsetExpression.GetSize() ) {
      EXCEPTION("You specified "<< offsetExpression.GetSize()
             << " displace expressions for a moving interface, for a " << coordSys_->GetDim() << " dimensional problem." << std::endl
             << "You must provide exactly as many math parser expressions as there are dimensions in the coordinate system.")
    }
    // assign the math parser
    SetMotionMathParser(offsetExpression);
    // check if the correct number of expressions was specified
  }


  void MortarInterface::SetRotation(const PtrParamNode motionNode) {
    // set expression for rotational motion...
    // check if the correct type of coordinate system was passed
    StdVector<std::string> offsetExpression;
    if ( coordSys_->GetDofName(2) != "phi" ) {
      EXCEPTION("For a rotating ncInterface the coordinate system must be either polar (in 2D) of cylindrical (3D).");
    }
    // get the rpm
    Double rpm = motionNode->Get("rpm")->As<Double>();
    // set the expression...
    std::ostringstream exprSstr("");
    // generate math parser expression for roatation, where 1 rpm = 360 / 60 s = 6 /s and we use a 10 digit precision
    exprSstr << std::setprecision(std::numeric_limits<Double>::digits10) << std::scientific << (6.0*rpm) << "*t";
    // assign to the offset expression
    offsetExpression.Resize(coordSys_->GetDim());
    offsetExpression[1] = exprSstr.str(); // second coordinate is always "phi" in cylindrical coord system
    // set math parser
    SetMotionMathParser(offsetExpression);
    // set the grid velocity by mathParser expression for counter rotation. This is only needed if isEulerian_ true
    if (isEulerian_){
      StdVector<std::string> veloExpr(coordSys_->GetDim());
      veloExpr[0] = "0.0";
      if (veloExpr.GetSize() == 3)
        veloExpr[2] = "0.0";
      exprSstr.str("");
      exprSstr.clear();
      // generate math parser expression for roatation, where 1 rpm = 360 / 60 s = 6 /s and we use a 10 digit precision
      // to make things generic, we pick up the base vectors of the coordinate system and use them to generate the expression
      // this is a bit of a hack, as the math parser does not handle 'r', so we express it in terms of the rotated Cartesian coordinates
      // and the origin 
      Vector<Double> origin = *coordSys_->GetOrigin();
      Matrix<Double> rotMat = *coordSys_->GetRotationMatrix();
      // -2*pi*rpm/60*sqrt(((x-x0)*a1 + (y-y0)*a2 + (z-z0)*a3)^2 + ((x-x0)*b1 + (y-y0)*b2 + (z-z0)*b3)^2)
      exprSstr << std::setprecision(std::numeric_limits<Double>::digits10) << std::scientific << "-2*pi*" << (rpm/60.0) << "*sqrt((" 
               << "(x-" << origin[0] << ")*" << rotMat[0][0] << "+(y-" << origin[1] << ")*" << rotMat[0][1] << "+(z-" << origin[2] << ")*" << rotMat[0][2] << ")^2 + (" 
               << "(x-" << origin[0] << ")*" << rotMat[1][0] << "+(y-" << origin[1] << ")*" << rotMat[1][1] << "+(z-" << origin[2] << ")*" << rotMat[1][2] << ")^2 )";
      veloExpr[1] = exprSstr.str();
      // generate the coefFunction
      gridVelo_ = CoefFunction::Generate(mParser_, Global::REAL, veloExpr);
      // assign the specified coordinate system
      gridVelo_->SetCoordinateSystem(coordSys_);
      LOG_DBG2(mortarInterface) << "Assigned grid velocity expression '" << veloExpr[1] << std::endl;
    }
  }


  void MortarInterface::SetConnectedRegions(const PtrParamNode motionNode) {
    // check for connected regions
    if(motionNode->Has("connectedRegions")) {
      std::string connectedRegions = motionNode->Get("connectedRegions")->As<std::string>();
      // Read "connectedRegions" tag from xml file to specify statically connected volume regions that
      // must be moved with the rotating region
      std::stringstream regionStream(connectedRegions);
      std::string regionName;
      // Read the passed regions line by line, tokenizing w.r.t. space ' '
      while(getline(regionStream, regionName, ' ')) {
        connectedVolRegions_.Push_back(ptGrid_->GetRegionId(regionName));
        LOG_DBG2(mortarInterface) << "Added region '" << regionName << "' to move with '" << ptGrid_->GetRegionName(movingVolRegion_) << "'." << std::endl;
      }
    }
    // if connected regions are specified, include these nodes as well
    if(connectedVolRegions_.GetSize() != 0) {
      for(auto itVolRegion : connectedVolRegions_) {
        StdVector<UInt> nodeIds;
        ptGrid_->GetNodesByRegion(nodeIds, itVolRegion);
        movingNodeIds_.Append(nodeIds);
      }
    }
  }


  void MortarInterface::SetMotionMathParser(const StdVector<std::string>& offsetExpression) {
    // set math parser expressions...
    mParser_ = domain->GetMathParser();
    for (UInt iDim = 0; iDim < 3; ++iDim) {
      if (iDim < coordSys_->GetDim() && !offsetExpression[iDim].empty()) {
        // get the handles for the mathparser for all specified offset expressions
        offsetMpHandles_[iDim] = mParser_->GetNewHandle(true);
        // assign the expression to the handle
        mParser_->SetExpr(offsetMpHandles_[iDim], offsetExpression[iDim]);
        LOG_DBG2(mortarInterface) << "Added offset expression '" << offsetExpression[iDim] << "' to the mathParser handle of dimension " << iDim << std::endl;
        // check if the expression is really time dependent.
        // Otherwise there will be no motion
        if ( mParser_->IsExprVariable( offsetMpHandles_[iDim], "t") )
          isMoving_ = true;
        else
          WARN("Offset expression '" << offsetExpression[iDim] << "' is not time dependent. Motion remains deactivated.");
      } else {
        // assign the global handle to the 3rd entry for 2D coord systems and for non-specified dimensions
        offsetMpHandles_[iDim] = MathParser::GLOB_HANDLER;
      }
    }
  }


  void MortarInterface::SetCoplanar(const PtrParamNode nciNode) {
    // check if interface is coplanar
    // if both surfaces are planar we can use 2 points from each surface to check if they are coplanar
    bool isPrimaryPlanar = ptGrid_->IsSurfacePlanar(primarySurfElems_);
    bool isSecondaryPlanar = ptGrid_->IsSurfacePlanar(secondarySurfElems_);
    if (isPrimaryPlanar && isSecondaryPlanar) {
      StdVector<SurfElem*> mixElems = primarySurfElems_;
      mixElems.Append(secondarySurfElems_);
      if (ptGrid_->IsSurfacePlanar(mixElems))
        isCoplanar_ = true;
    }
    // if only one interface is planar we have a problem
    if (!isCoplanar_ && (isPrimaryPlanar && isSecondaryPlanar))
      WARN("Non-conforming interface '" << name_ << "' is assumed to be curved even as the two planar surfaces that are not coplanar!" << std::endl <<
           "This might lead to errors!");
    if ((!isPrimaryPlanar && isSecondaryPlanar) || (isPrimaryPlanar && !isSecondaryPlanar))
      WARN("Non-conforming interface '" << name_ << "' is assumed to be curved even though " << std::endl <<
           "isPrimaryPlanar=" << isPrimaryPlanar << " and isSecondaryPlanar=" << isSecondaryPlanar << "!" << std::endl <<
           "This might lead to errors!");
  }


  void MortarInterface::CheckMeshValidity() {
    if (ptGrid_->GetDim() == 3) {
      // lambda function to call primary and secondary side individually
      auto CheckElemTypes = [&] (StdVector<SurfElem*>& elemVec) {
        UInt numElems = elemVec.GetSize();
        UInt elemTypeCounterQuad4 = 0;
        UInt elemTypeCounterQuad8 = 0;
        UInt elemTypeCounterQuad9 = 0;
        UInt elemTypeCounterTria3 = 0;
        UInt elemTypeCounterTria6 = 0;
        // loop over primary interface elements and check for type
        for (UInt iElem = 0; iElem < numElems; ++iElem) {
          switch(elemVec[iElem]->type) {
            case Elem::ET_QUAD4: {
              ++elemTypeCounterQuad4;
              break;
            }
            case Elem::ET_QUAD8: {
              ++elemTypeCounterQuad8;
              break;
            }
            case Elem::ET_QUAD9: {
              ++elemTypeCounterQuad9;
              break;
            }
            case Elem::ET_TRIA3: {
              ++elemTypeCounterTria3;
              break;
            }
            case Elem::ET_TRIA6: {
              ++elemTypeCounterTria6;
              break;
            }
            default: {
              WARN("Found unexpected element type " << elemVec[iElem]->type << " in non-conforming interface '" << name_ << "'.");
            }
          }
        }
        // check for quadratic elements
        if (elemTypeCounterQuad8 > 0 || elemTypeCounterQuad9 > 0 || elemTypeCounterTria6 > 0)
          WARN("You provided a surface with quadratic element type for the non-conforming interface '" << name_ << "'." << std::endl
            << "OpenCFS will treat them as linear for computing intersection elements!");
        // log
        LOG_DBG2(mortarInterface) << "The interface side consists of elements:" << std::endl
                                  << elemTypeCounterQuad4 << " QUAD4" << std::endl
                                  << elemTypeCounterQuad8 << " QUAD8" << std::endl
                                  << elemTypeCounterQuad9 << " QUAD9" << std::endl
                                  << elemTypeCounterTria3 << " TRIA3" << std::endl
                                  << elemTypeCounterTria6 << " TRIA6" << std::endl;
        // the NCI_INTERSECT_RECT works only with linear quadrilaterals
        if(intersectAlgo_ == NCI_INTERSECT_RECT && elemTypeCounterQuad4 < numElems)
          WARN("Only linear quadrilaterals can be intersected with coaxial rectangle algorithm!");
      };
      // check primary side:
      LOG_DBG2(mortarInterface) << "Checking element types of primary side." << std::endl;
      CheckElemTypes(primarySurfElems_);
      // check secondary side:
      LOG_DBG2(mortarInterface) << "Checking element types of secondary side." << std::endl;
      CheckElemTypes(secondarySurfElems_);
    }
    else if (ptGrid_->GetDim() == 2) {
      // lambda function to call primary and secondary side individually
      auto CheckElemTypes = [&] (StdVector<SurfElem*>& elemVec) {
        UInt numElems = elemVec.GetSize();
        UInt elemTypeCounterLine2 = 0;
        UInt elemTypeCounterLine3 = 0;
        // loop over primary interface elements and check for type
        for (UInt iElem = 0; iElem < numElems; ++iElem) {
          switch(elemVec[iElem]->type) {
            case Elem::ET_LINE2: {
              ++elemTypeCounterLine2;
              break;
            }
            case Elem::ET_LINE3: {
              ++elemTypeCounterLine3;
              break;
            }
            default: {
              EXCEPTION("Found unexpected element type " << elemVec[iElem]->type << " in non-conforming interface '" << name_ << "'.");
            }
          }
        }
        // check for quadratic elements
        if (elemTypeCounterLine3 > 0)
          WARN("You provided a surface with quadratic element type for the non-conforming interface '" << name_ << "'." << std::endl
                << "OpenCFS will treat them as linear for computing intersection elements!");
        // log
        LOG_DBG2(mortarInterface) << "The interface side consists of elements:" << std::endl
                                  << elemTypeCounterLine2 << " ET_LINE2" << std::endl
                                  << elemTypeCounterLine3 << " ET_LINE3" << std::endl;
      };
      // check primary side:
      LOG_DBG2(mortarInterface) << "Checking element types of primary side." << std::endl;
      CheckElemTypes(primarySurfElems_);
      // check secondary side:
      LOG_DBG2(mortarInterface) << "Checking element types of secondary side." << std::endl;
      CheckElemTypes(secondarySurfElems_);
    } else {
      EXCEPTION("Invalid Grid Dimension!!");
    }
  }


  void MortarInterface::MoveInterface() {
    // if no motion is specified or we use iterative coupling or at initialization, return
    if (!isMoving_ || useMeshSmoothing_ || !isReset_)
      return;
    // coordinates and offsets
    Vector<Double> movedLocalCoord;
    Vector<Double> movedGlobalCoord;
    Vector<Double> globalNodeOffset;
    Vector<Double> localNodeOffset;

    // compute node offsets...
    globalNodeOffset.Resize(movingNodeIds_.GetSize() * coordSys_->GetDim());
    Integer numNodes = movingNodeIds_.GetSize(); // sanity check for type conversion
    if (numNodes < 0)
      EXCEPTION("MortarInterface: the number of nodes in the grid exceeds the integer limits!");
    if (coordSysName_ == "default") {
      // for the default system we can assign directly
      #pragma omp parallel for num_threads(CFS_NUM_THREADS)
      for (Integer iNodes = 0; iNodes < numNodes; ++iNodes) { // omp only allows us to use int here
        for (UInt iDim = 0; iDim < coordSys_->GetDim(); ++iDim) {
          globalNodeOffset[iNodes * coordSys_->GetDim() + iDim] = mParser_->Eval(offsetMpHandles_[iDim]);
        }
      }
    } else {
      // for other coordinate systems we need to take further steps...
      // evaluate the node offsets via math parser
      localNodeOffset.Resize(coordSys_->GetDim());
      for (UInt iDim = 0; iDim < coordSys_->GetDim(); ++iDim) {
        localNodeOffset[iDim] = mParser_->Eval(offsetMpHandles_[iDim]);
      }

      // compute the coordinate offset
      #pragma omp parallel for private (movedGlobalCoord, movedLocalCoord) num_threads(CFS_NUM_THREADS)
      for (Integer iNodes = 0; iNodes < numNodes; ++iNodes) { // omp only allows us to use int here
        movedGlobalCoord.Resize(coordSys_->GetDim());
        movedLocalCoord.Resize(coordSys_->GetDim());
        // get node offsets via the local coord system
        for (UInt iDim = 0; iDim < coordSys_->GetDim(); ++iDim) {
          movedLocalCoord[iDim] = movingNodeLocalCoords_[iNodes][iDim] + localNodeOffset[iDim];
        }
        // transform to moved global coordinate
        coordSys_->Local2GlobalCoord(movedGlobalCoord, movedLocalCoord);
        // compute offsets
        globalNodeOffset[iNodes * coordSys_->GetDim()    ] = movedGlobalCoord[0] - movingNodeCoords_[iNodes][0];
        globalNodeOffset[iNodes * coordSys_->GetDim() + 1] = movedGlobalCoord[1] - movingNodeCoords_[iNodes][1];
        if (coordSys_->GetDim() == 3)
          globalNodeOffset[iNodes * coordSys_->GetDim() + 2] = movedGlobalCoord[2] - movingNodeCoords_[iNodes][2];
      }
    }
    // assign offsets to the grid.
    ptGrid_->SetNodeOffset(movingNodeIds_, globalNodeOffset);
  }


  void MortarInterface::ResetInterface(){
    // we only reset in the initial run or if the interface is moving
    if (!isMoving_ && elemList_->GetSize() > 0)
      return;
    // set variables
    StdVector<std::string> namedNodesLists;
    std::string nciNodesName = name_ + "_nodes";
    // clear all intersection elements
    elemList_->Clear(true);
    // delete all (surface) elements of the current region
    // these are only set if exportToGrid_ is true
    if (nciRegion_ != NO_REGION_ID && exportToGrid_)
      ptGrid_->ClearRegion(nciRegion_);
    // get the names of stored node lists containing the nci nodes
    ptGrid_->GetListNodeNames(namedNodesLists);
    // delete the nci nodes
    if (namedNodesLists.Find(nciNodesName) != -1)
      ptGrid_->DeleteNamedNodes(nciNodesName);
    // set flag
    isReset_ = true;
  }


  void MortarInterface::UpdateInterface() {
    // either we set motion, initialize, or update the nci after reset. Otherwise return
    if (!isMoving_ && (elemList_->GetSize() > 0 || isReset_))
      return;
    // name to add the new nodeset as region to the grid
    std::string newNodesName = name_ + "_nodes";

    // try moving the interface one step
    MoveInterface();
    
    // reset the isReset_ flag to indicate that the interface was updated
    isReset_ = false;

    // create lists of candidates for element intersection...
    // we can pre compute the intersection candidates to speed up computation if we use the CGAL build
    // without CGAL we still assign every possible combination
    ComputeIntersectionCandidates();

    // call the specified intersection algorithm...
    // node ids of the intersection elements
    StdVector<UInt> newNodeIds;
    // ATTENTION: the added node ids are not consistent when using multi threading! check this..
    ComputeIntersection(newNodeIds);

    // add named node list of the nodes of the computed intersection elements as named nodes to the grid
    if (newNodeIds.GetSize() > 0)
      ptGrid_->AddNamedNodes(newNodesName, newNodeIds);
    LOG_DBG2(mortarInterface) << "Added " << newNodeIds.GetSize() << " nodes to named nodes '" << newNodesName << "'." << std::endl;
    LOG_DBG3(mortarInterface) << "The new node Ids are:" << std::endl << newNodeIds << std::endl;

    // if elements were found, update the surface integrators and check for grid export
    if(elemList_->GetSize() > 0) {
      // add the intersection elements to the BiLinFormContext
      UpdateIntegrators();
      // export the interface to the grid if specified
      ExportToGrid();
      // increase the call counter
      callCounter_++;
    } else {
      EXCEPTION("No intersection elements were computed for non-conforming interface '" << name_
             << "'. Please check your mesh file. " << "Different precision of interface neighbors?");
    }
  }


  void MortarInterface::ExportToGrid() {
    // So far the ncis for every timestep are added to the grid. 
    // However, they are not visible in the output hdf5 as the grid gets only printed at initialization. See SimOutputHDF5::WriteGrid().
    // if not specified, return
    if (!exportToGrid_)
      return;
    std::string regionName;
    // only store the initial grid if we do not move. Store all if we do.
    if (!isMoving_ && callCounter_ ==0) {
      regionName = name_;
    }
    else if (isMoving_) {
      regionName = name_;
      regionName.append("_t");
      regionName.append(std::to_string(callCounter_));
    } else {
      return;
    }
    UInt numElems = elemList_->GetSize();
    StdVector<UInt> ncElemIds;
    // helper to store explicit copies needed to output the intersection grid
    StdVector<SurfElem*> ncElemsHelper;
    // add the new region to the grid
    nciRegion_ = ptGrid_->AddSurfaceRegion(regionName);
    // create helper elements.
    // We need to make explicit copies of the NcSurfElems, because the
    // Grid deletes all its elements when it gets destroyed.
    ncElemsHelper.Resize(numElems);
    for (UInt iElems = 0; iElems < numElems; ++iElems) {
      ncElemsHelper[iElems] = new SurfElem(*(elemList_->GetSurfElem(iElems)));
    }
    ptGrid_->AddSurfaceElems(nciRegion_, ncElemsHelper, ncElemIds);

    // assign the dimension of the output intersection grid to the grid and check if all elements have the same dimension
    for (UInt iElems = 0; iElems < numElems; ++iElems) {
      // search for entry in the entityDim
      std::map<std::string, UInt>::iterator it = ptGrid_->entityDim_.find(regionName);
      if( it != ptGrid_->entityDim_.end() ) {
        if( it->second != Elem::shapes[elemList_->GetSurfElem(iElems)->type].dim ) {
          // if the dimension of the current element is different raise an exception
          EXCEPTION( "Region '" << regionName << "' contains elements of different dimensions!");
        }
        // if nothing has been assigned yet assign
      } else {
        ptGrid_->entityDim_[regionName] = Elem::shapes[elemList_->GetSurfElem(iElems)->type].dim;
      }
    }
    LOG_DBG2(mortarInterface) << "Added surface region " << nciRegion_ << " with name '" << regionName << "' to the grid." << std::endl;
  }


  void MortarInterface::ComputeIntersection(StdVector<UInt>& newNodeIds) {
    // parse the intersectAlgo_ and call the correct intersection operation to compute the new intersection-element nodes
    LOG_DBG2(mortarInterface) << "Computing element intersections using algorithm '" << intersectAlgo_  << "'" << std::endl;

    UInt primaryIndex, secondaryIndex;
    bool doIntersect;
    switch (intersectAlgo_) {
      case NCI_INTERSECT_LINE:
        for (UInt iCandidate = 0; iCandidate < intersectionCandiatesIdx_.size(); ++iCandidate) {
          primaryIndex = intersectionCandiatesIdx_[iCandidate].first;
          secondaryIndex = intersectionCandiatesIdx_[iCandidate].second;
          // compute intersection
          doIntersect = IntersectLines(primarySurfElems_[primaryIndex], secondarySurfElems_[secondaryIndex], newNodeIds);
          LOG_DBG3(mortarInterface) << "Elements '" << primarySurfElems_[primaryIndex]->elemNum
                                      << "' and '" << secondarySurfElems_[secondaryIndex]->elemNum
                                      << "' intersect: "<< doIntersect;
        }
        break;
      case NCI_INTERSECT_RECT:
        for (UInt iCandidate = 0; iCandidate < intersectionCandiatesIdx_.size(); ++iCandidate) {
          primaryIndex = intersectionCandiatesIdx_[iCandidate].first;
          secondaryIndex = intersectionCandiatesIdx_[iCandidate].second;
          // compute intersection
          doIntersect = IntersectRects( primarySurfElems_[primaryIndex], secondarySurfElems_[secondaryIndex], newNodeIds);
          LOG_DBG3(mortarInterface) << "Elements '" << primarySurfElems_[primaryIndex]->elemNum
                                      << "' and '" << secondarySurfElems_[secondaryIndex]->elemNum
                                      << "' intersect: "<< doIntersect;
        }
        break;
      case NCI_INTERSECT_POLYGON: {
        Integer numCandidates = intersectionCandiatesIdx_.size(); // sanity check for type conversion
        if (numCandidates < 0)
          EXCEPTION("MortarInterface: the number of intersection candidates in the grid exceeds the integer limits!");
      #ifdef USE_CGAL
        if (isCoplanar_)
          IntersectCoplanarPolygons(newNodeIds);
        else
          IntersectGeneralPolygons(newNodeIds);
      #else
        #pragma omp parallel num_threads(CFS_NUM_THREADS) private(primaryIndex, secondaryIndex, doIntersect)
        {
          // caching of data for each thread separatly
          // create temporary vectors and matrices here (to speed up computation?)
          // these are private variables as they are defined inside the parallel region
          StdVector<UInt> tempNewNodeIds;
          StdVector<Vector<Double>> primaryNodeCoords;
          StdVector<Vector<Double>> secondaryNodeCoords;
          StdVector<Vector<Double>> intersectionPolygonCoords;
          Vector<Double> tempVec1;
          Vector<Double> tempVec2;
          Vector<Double> tempVec3;
          Vector<Double> tempVec4;
          Vector<Double> tempVec5;
          Matrix<Double> rotMat;
          Matrix<Double> rotMatTrans;
          StdVector<Vector<Double>> rotPolyPrim;
          StdVector<Vector<Double>> rotPolySec;
          tempVec1.Resize(ptGrid_->GetDim());
          tempVec2.Resize(ptGrid_->GetDim());
          tempVec3.Resize(ptGrid_->GetDim());
          tempVec4.Resize(ptGrid_->GetDim());
          tempVec5.Resize(ptGrid_->GetDim());
          rotMat.Resize(ptGrid_->GetDim(),ptGrid_->GetDim());
          rotMatTrans.Resize(ptGrid_->GetDim(),ptGrid_->GetDim());
          rotPolyPrim.Resize(ptGrid_->GetDim());
          rotPolySec.Resize(ptGrid_->GetDim());
          #pragma omp for
          for (Integer iCandidate = 0; iCandidate < (Integer)intersectionCandiatesIdx_.size(); ++iCandidate) {
            primaryIndex = intersectionCandiatesIdx_[iCandidate].first;
            secondaryIndex = intersectionCandiatesIdx_[iCandidate].second;
            // compute intersection
            doIntersect = IntersectPolygons(primarySurfElems_[primaryIndex],
                                            secondarySurfElems_[secondaryIndex],
                                            tempNewNodeIds,
                                            primaryNodeCoords,
                                            secondaryNodeCoords,
                                            intersectionPolygonCoords,
                                            tempVec1, tempVec2, tempVec3, tempVec4, tempVec5,
                                            rotMat, rotMatTrans, rotPolyPrim, rotPolySec);
            LOG_DBG3(mortarInterface) << "Elements '" << primarySurfElems_[primaryIndex]->elemNum
                                      << "' and '" << secondarySurfElems_[secondaryIndex]->elemNum
                                      << "' intersect: "<< doIntersect;
            }
          #ifdef USE_OPENMP
            // add computed nodes to the general node vector
            omp_set_lock(&newNodesLock_);
          #endif
            for (UInt iNodes = 0; iNodes < tempNewNodeIds.GetSize(); ++iNodes) {
              newNodeIds.Push_back(tempNewNodeIds[iNodes]);
            }
          #ifdef USE_OPENMP
            omp_unset_lock(&newNodesLock_);
          #endif
        }
      #endif
        break;
      }
      default:
        EXCEPTION("Unknown intersection algorithm '" << intersectAlgo_ <<"'.");
        break;
    }
  }

#ifdef USE_CGAL
  void MortarInterface::ComputeIntersectionCandidates()
  {
    if (cakePieceProjection_) {
      // in the case caePieceProjection The intersection simply check all possible permutations
      intersectionCandiatesIdx_.resize(primarySurfElems_.GetSize() * secondarySurfElems_.GetSize());
      UInt position = 0;
      for (UInt iPrimElems = 0; iPrimElems < primarySurfElems_.GetSize(); ++iPrimElems) {
        for (UInt iSecElems = 0; iSecElems < secondarySurfElems_.GetSize(); ++iSecElems) {
          intersectionCandiatesIdx_[position].first = iPrimElems;
          intersectionCandiatesIdx_[position].second = iSecElems;
          position++;
        }
      }
    }
    else {
      // clear intersection candidates
      intersectionCandiatesIdx_.clear();
      // (Re)compute bounding boxes if not done yet or if the region has moved...
      // primary region
      if (primarySurfElemBoxes_.size() != primarySurfElems_.GetSize() || movingVolRegion_ == primaryVolRegion_) {
        primarySurfElemBoxes_.resize(primarySurfElems_.GetSize());
        primaryBoxIndexes_.Resize(primarySurfElems_.GetSize());
#pragma omp parallel for num_threads(CFS_NUM_THREADS)
        for (Integer iBox = 0; iBox < (Integer)primarySurfElems_.GetSize(); ++iBox) {
          // create the box
          boost::array<Double, 6> bbox;
          ptGrid_->CreateBBoxFromElement(primarySurfElems_[iBox], boundingBoxTolerance_, &bbox[0], isMoving_);
          // assign the box index
          primaryBoxIndexes_[iBox] = iBox;
          // kirill: translational p.b.c.
          // a transformation must be applied to the bounding box
          if (mutualProjection_) {
            // translation
            for (UInt iDim = 0; iDim < ptGrid_->GetDim(); ++iDim) {
              bbox[iDim] -= translationVector_[iDim];
              bbox[iDim + 3] -= translationVector_[iDim];
            }
          } 
          // get a handle to the bounding box
          Grid::HandleBox boxHandle(Grid::BBox3D(bbox[0], bbox[1], bbox[2],
                                                 bbox[3], bbox[4], bbox[5]),
                                    &primaryBoxIndexes_[iBox]);
          // assign the handle
          primarySurfElemBoxes_[iBox] = boxHandle;
        }
      }
      // secondary region
      if (secondarySurfElemBoxes_.size() != secondarySurfElems_.GetSize() || movingVolRegion_ == secondaryVolRegion_) {
        secondarySurfElemBoxes_.resize(secondarySurfElems_.GetSize());
        secondaryBoxIndexes_.Resize(secondarySurfElems_.GetSize());
#pragma omp parallel for num_threads(CFS_NUM_THREADS)
        for (Integer iBox = 0; iBox < (Integer)secondarySurfElems_.GetSize(); ++iBox) {
          // create the box
          boost::array<Double, 6> bbox;
          ptGrid_->CreateBBoxFromElement(secondarySurfElems_[iBox], boundingBoxTolerance_, &bbox[0], isMoving_);
          // assign the box index
          secondaryBoxIndexes_[iBox] = iBox;
          // get a handle to the bounding box
          Grid::HandleBox boxHandle(Grid::BBox3D(bbox[0], bbox[1], bbox[2],
                                                 bbox[3], bbox[4], bbox[5]),
                                    &secondaryBoxIndexes_[iBox]);
          // assign the handle
          secondarySurfElemBoxes_[iBox] = boxHandle;
        }
      }
      CGAL::box_intersection_d(primarySurfElemBoxes_.begin(), primarySurfElemBoxes_.end(),
                               secondarySurfElemBoxes_.begin(), secondarySurfElemBoxes_.end(),
                               elemElemIdReporter(std::back_inserter(intersectionCandiatesIdx_)));
    }
  }
#else
  void MortarInterface::ComputeIntersectionCandidates() {
    // in the default build we still have no proper pre computation. The intersection candidates consist of all possible permuations
    intersectionCandiatesIdx_.resize(primarySurfElems_.GetSize() * secondarySurfElems_.GetSize());
    UInt position = 0;
    for (UInt iPrimElems = 0; iPrimElems < primarySurfElems_.GetSize(); ++iPrimElems) {
      for (UInt iSecElems = 0; iSecElems < secondarySurfElems_.GetSize(); ++iSecElems) {
        intersectionCandiatesIdx_[position].first = iPrimElems;
        intersectionCandiatesIdx_[position].second = iSecElems;
        position++;
      }
    }
  }
#endif


  /****************************************************************************
   **
  ** IntersectLines
  **
  **   computes the local coordinates of the overlap of the primary and secondary
  **   element with respect to the primary side in the order of the
  **   orientation of the secondary side element. It pushes back the intersection
  **   element to elemList.
  **
  ** Input Parameters:
  **   ifaceElem1:  Master Side
  **   ifaceElem2:  Slave Side
  **
  ** Output Parameters:
  **   elemList: the found intersection NCElems will be pushed
  **                     back to this vector
  **
  */


  bool MortarInterface::IntersectLines( SurfElem *ifaceElem1,
                                        SurfElem *ifaceElem2,
                                        StdVector<UInt> &newNodes )
  {
    // c0, c1, d0 and d1 are the endpoints of the two line elements
    //
    //           d0 x-----------+--------------x d1
    // c0 x---------+-----------x c1

    // primary elem node coordinates
    Vector<Double> c0(this->ptGrid_->GetDim());
    Vector<Double> c1(this->ptGrid_->GetDim());
    // secondary elem node coordinates
    Vector<Double> d0(this->ptGrid_->GetDim());
    Vector<Double> d1(this->ptGrid_->GetDim());
    // normal vector for projection
    Vector<Double> normal(this->ptGrid_->GetDim());
    // temporary vector
    Vector<Double> tmp(this->ptGrid_->GetDim());
    Vector<Double> tmp2(2);
    // local points in line1 representing the nodes of line2
    Vector<Double> intersectionCoords(2);
    // connectivity of intersection nodes
    StdVector<UInt> connect2(2);
    // axis representing the orientation of line1
    Vector<Double> xUnitVec(this->ptGrid_->GetDim());

    bool wasProjected0 = false;
    bool wasProjected1 = false;

    Double dist;
    UInt nodenum_c0, nodenum_c1, nodenum_d0, nodenum_d1;
    Double relativeElemVol = 0.0;

    // Get coordinates of the endpoints
    nodenum_c0 = ifaceElem1->connect[0];
    nodenum_c1 = ifaceElem1->connect[1];
    nodenum_d0 = ifaceElem2->connect[0];
    nodenum_d1 = ifaceElem2->connect[1];
    ptGrid_->GetNodeCoordinate(c0, nodenum_c0, isMoving_);
    ptGrid_->GetNodeCoordinate(c1, nodenum_c1, isMoving_);
    ptGrid_->GetNodeCoordinate(d0, nodenum_d0, isMoving_);
    ptGrid_->GetNodeCoordinate(d1, nodenum_d1, isMoving_);

    // for translational p.b.c., project the primary grid nodes onto the secondary interface
    if (mutualProjection_) {
      c0 -= translationVector_;
      c1 -= translationVector_;
      wasProjected0 = true;
      wasProjected1 = true;
    }

    // for cake-piece projection, perform translation and rotation of coordinates
    if (cakePieceProjection_) {
      // debug messages. Remove or make as debug logs before merge
      // std::cout << "c0: " << c0 << std::endl;
      // std::cout << "c1: " << c1 << std::endl;
      // subtract rotation center
      c0 -= cakePieceRotationCenter_;
      c1 -= cakePieceRotationCenter_;
      // rotate around z axis
      Vector<Double> tmp(2);
      tmp[0] = c0[0]*cos(cakePieceRotationAngle_) - c0[1]*sin(cakePieceRotationAngle_);
      tmp[1] = c0[0]*sin(cakePieceRotationAngle_) + c0[1]*cos(cakePieceRotationAngle_);
      c0 = tmp;
      tmp[0] = c1[0]*cos(cakePieceRotationAngle_) - c1[1]*sin(cakePieceRotationAngle_);
      tmp[1] = c1[0]*sin(cakePieceRotationAngle_) + c1[1]*cos(cakePieceRotationAngle_);
      c1 = tmp;
      // translate back to original position
      c0 += cakePieceRotationCenter_;
      c1 += cakePieceRotationCenter_;
      // translate to primary interface position
      c0 -= translationVector_;
      c1 -= translationVector_;
      // set projection bools to true since we already projected a line and we do not need
      // to wait for the check below
      wasProjected0 = true;
      wasProjected1 = true;
      // debug messages. Remove or make as debug logs before merge
      // std::cout << "c0: " << c0 << std::endl;
      // std::cout << "c1: " << c1 << std::endl;
      // std::cout << "c0-center: " << c0 << std::endl;
      // std::cout << "c1-center: " << c1 << std::endl;
      // std::cout << "c0-center+rot: " << c0 << std::endl;
      // std::cout << "c1-center+rot: " << c1 << std::endl;
      // std::cout << "c0-center+rot+center: " << c0 << std::endl;
      // std::cout << "c1-center+rot+center: " << c1 << std::endl;
      // std::cout << "c0-final: " << c0 << std::endl;
      // std::cout << "c1-final: " << c1 << std::endl;
    }

    // Project primary nodes onto secondary element, if interface is not coplanar
    if (!isCoplanar_) {
      shared_ptr<ElemShapeMap> esmSecondary = ptGrid_->GetElemShapeMap(ifaceElem2, isMoving_);
      LocPoint midPointSecondary = Elem::shapes[ifaceElem2->type].midPointCoord;

      // compute maximal allowed distance as sum of lengths of both lines
      tmp = c1 - c0;
      Double maxDist = tmp.NormL2();
      tmp = d1 - d0;
      maxDist += tmp.NormL2();
      // compute normal vector of secondary element
      esmSecondary->CalcNormal(normal, midPointSecondary);
      // compute distance of c0 to plane of secondary element
      tmp = c0 - d0;
      dist = normal.Inner(tmp);
      // make sure that distance does not exceed maximum distance
      if (fabs(dist) > maxDist)
        return false;
      // do the projection if necessary
      if (fabs(dist) > absoluteTolerance_) {
        c0 -= normal * dist;
        wasProjected0 = true;
      }
      // do the same for c1
      tmp = c1 - d0;
      dist = normal.Inner(tmp);
      if (fabs(dist) > absoluteTolerance_) {
        c1 -= normal * dist;
        wasProjected1 = true;
      }
    }

    // Compute and normalize vector from c0 to c1.
    // This becomes the new x-unit vector.
    xUnitVec = c1 - c0;
    dist = xUnitVec.NormL2();
    // check if both elements were perpendicular
    if (dist < absoluteTolerance_)
      return false;
    xUnitVec /= dist;

    // Compute x1 coordinate of line2 in respect to line1.
    tmp = d0 - c0;
    tmp.Inner(xUnitVec, tmp2[0]);
    tmp2[0] /= dist;

    // Compute x2 coordinate of line2 in respect to line1.
    tmp = d1 - c0;
    tmp.Inner(xUnitVec, tmp2[1]);
    tmp2[1] /= dist;

    // Bring line2's endpoints into ascending order.
    if(tmp2[1] < tmp2[0]) {
      intersectionCoords[0] = tmp2[1];
      intersectionCoords[1] = tmp2[0];
      connect2[0] = nodenum_d1;
      connect2[1] = nodenum_d0;
    } else {
      intersectionCoords[0] = tmp2[0];
      intersectionCoords[1] = tmp2[1];
      connect2[0] = nodenum_d0;
      connect2[1] = nodenum_d1;
    }

    // Check if an intersection between line1 and line2 exists.
    if((intersectionCoords[0] >= 1.0) || intersectionCoords[1] <= 0.0)
      return false;

    // create new mortar element
    shared_ptr<MortarNcSurfElem> ncElem(new MortarNcSurfElem());
    ncElem->connect.Resize(2);

    // Here we already know that an intersection exists, but we must distinguish 4 different cases.
    tmp.Resize(3); // use tmp for assigning the new node
    if(intersectionCoords[0] <= 0) {
      // connect2[0] x--------|------------...
      //                   c0 x--------------------x c1
      if (wasProjected0) {
        // create a new node for the projection of c0 onto d
        Vector<Double> new_node;
        new_node.Resize(3);
        new_node[0] = c0[0];
        new_node[1] = c0[1];
        if (this->ptGrid_->GetDim() == 2)
          new_node[2] = 0.0;
        else
          new_node[2] = c0[2];
        ptGrid_->AddNode(new_node, nodenum_c0);
        newNodes.Push_back(nodenum_c0);
        ncElem->connect[0] = nodenum_c0;
      } else {
        ncElem->connect[0] = nodenum_c0;
      }

      if(intersectionCoords[1] >= 1) {
        // connect2[0] x--------|--------------------|-----x connect2[1]
        //                   c0 x--------------------x c1
        if (wasProjected1) {
          // create a new node for the projection of c1 onto d
          Vector<Double> new_node;
          new_node.Resize(3);
          new_node[0] = c1[0];
          new_node[1] = c1[1];
          if (this->ptGrid_->GetDim() == 2)
            new_node[2] = 0.0;
          else
            new_node[2] = c1[2];
          ptGrid_->AddNode(new_node, nodenum_c1);
          newNodes.Push_back(nodenum_c1);
          ncElem->connect[1] = nodenum_c1;
        } else {
          ncElem->connect[1] = nodenum_c1;
        }
        relativeElemVol = 1;
      } else {
        // connect2[0] x-----|-------x connect2[1]
        //                c0 x-------|------x c1
        relativeElemVol = intersectionCoords[1];
        ncElem->connect[1] = connect2[1];
      }
    } else { // if(intersectionCoords[0] <= 0)
      //   connect2[0] x---------...
      //      c0 x-----|----------------x c1
      ncElem->connect[0] = connect2[0];

      if(intersectionCoords[1] >= 1) {
        // connect2[0] x----------------|------x connect2[1]
        //    c0 x-----|----------------x c1
        if (wasProjected1) {
          // create a new node for the projection of c1 onto d
          Vector<Double> new_node;
          new_node.Resize(3);
          new_node[0] = c1[0];
          new_node[1] = c1[1];
          if (this->ptGrid_->GetDim() == 2)
            new_node[2] = 0.0;
          else
            new_node[2] = c1[2];
          ptGrid_->AddNode(new_node, nodenum_c1);
          newNodes.Push_back(nodenum_c1);
          ncElem->connect[1] = nodenum_c1;
        } else {
          ncElem->connect[1] = nodenum_c1;
        }
        relativeElemVol = 1.0 - intersectionCoords[0];
      } else {
        // connect2[0] x------------x connect2[1]
        //      c0 x---|------------|---x c1
        ncElem->connect[1] = connect2[1];
        relativeElemVol = intersectionCoords[1] - intersectionCoords[0];
      }
    }
    // reject very small line elements.
    if (relativeElemVol < minRelativeSideLength_) {
      LOG_DBG3(mortarInterface) << "Rejecting line intersection element due to a relative size of " << relativeElemVol << std::endl
            << "  for intersection of elements " << ifaceElem1->elemNum << " and " << ifaceElem2->elemNum;
      return false;
    }

    // In case of a curved interface store the projected primary element.
    // This is needed for coordinate transform of integration points.
    shared_ptr<SurfElem> projElem = nullptr;
    if (!isCoplanar_) {
      Vector<Double> new_node(3);
      projElem.reset(new SurfElem());
      projElem->type = Elem::ET_LINE2;
      projElem->connect.Resize(2);
      if (wasProjected0) {
        // create a new node for the projection of c0 onto d
        new_node[0] = c0[0];
        new_node[1] = c0[1];
        if (this->ptGrid_->GetDim() == 2)
          new_node[2] = 0.0;
        else
          new_node[2] = c0[2];
        ptGrid_->AddNode(new_node, projElem->connect[0]);
        newNodes.Push_back(projElem->connect[0]);
      } else {
        projElem->connect[0] = nodenum_c0;
      }
      if (wasProjected1) {
        // create a new node for the projection of c1 onto d
        new_node[0] = c1[0];
        new_node[1] = c1[1];
        if (this->ptGrid_->GetDim() == 2)
          new_node[2] = 0.0;
        else
          new_node[2] = c1[2];
        ptGrid_->AddNode(new_node, projElem->connect[1]);
        newNodes.Push_back(projElem->connect[1]);
      } else {
        projElem->connect[1] = nodenum_c1;
      }
    }

    // Finally assign intersection element
    ncElem->type = Elem::ET_LINE2;
    ncElem->ptPrimary = ifaceElem1;
    ncElem->ptSecondary = ifaceElem2;
    ncElem->transVect = (mutualProjection_ || cakePieceProjection_) ? &GetTranslationVector() : nullptr;
    ncElem->rotationAngle = cakePieceRotationAngle_;
    ncElem->rotationCenter = cakePieceRotationCenter_;
    ncElem->projectedPrimary = projElem;
    elemList_->AddElement(ncElem);
    return true;
  }


  bool MortarInterface::IntersectRects( SurfElem *ifaceElem1, SurfElem *ifaceElem2, StdVector<UInt> &newNodes ) {
    Vector<Double> c0, c1, c2, d0, d1, d2;
    Vector<Double> diffS, diffX, diffY, diffX2;
    Vector<Double> s, t;
    StdVector<UInt> connect2;
    Double distX, distY, distX2, facX, facY, r;
    UInt nodeNr;
    // Introduce a tolerance to account for roundoff errors during the calculation of
    // normed new x base vector.
    Double tol_r;

    s.Resize(4);
    t.Resize(4);
    connect2.Resize(4);

    // The meaning of the points c0, c1, c2, d0 and d2
    // is as follows:
    //                x------------------x d2
    //                |                  |
    //                |                  |
    // c2 x-----------+----------x       |
    //    |           |          |       |
    //    |           |          |       |
    //    |        d0 x----------+-------x d1
    //    |                      |
    //    |                      |
    // c0 x----------------------x c1


    // Get coordinates of the endpoints
    ptGrid_->GetNodeCoordinate(c0, ifaceElem1->connect[0], isMoving_);
    ptGrid_->GetNodeCoordinate(c1, ifaceElem1->connect[1], isMoving_);
    ptGrid_->GetNodeCoordinate(c2, ifaceElem1->connect[3], isMoving_);
    ptGrid_->GetNodeCoordinate(d0, ifaceElem2->connect[0], isMoving_);
    ptGrid_->GetNodeCoordinate(d1, ifaceElem2->connect[1], isMoving_);
    ptGrid_->GetNodeCoordinate(d2, ifaceElem2->connect[2], isMoving_);

    // Compute and normalize vector from c0 to c1.
    // This becomes the new x-unit vector.
    diffX = c1 - c0;
    distX = diffX.NormL2();
    facX = 1.0 / distX;
    diffX *= facX;

    // Compute and normalize vector from c0 to c2
    // This becomes the new y-unit vector.
    diffY = c2 - c0;
    distY = diffY.NormL2();
    facY = 1.0 / distY;
    diffY *= facY;

    // Now compute vector from c0 to d0 and project
    // the result onto the new x- and y-axis.
    diffS = d0 - c0;
    diffS.Inner(diffX, s[0]);
    diffS.Inner(diffY, s[1]);
    s[0] *= facX;
    s[1] *= facY;

    // Now compute vector from c0 to d2 and project
    // the result onto the new x- and y-axis.
    diffS = d2 - c0;
    diffS.Inner(diffX, s[2]);
    diffS.Inner(diffY, s[3]);
    s[2] *= facX;
    s[3] *= facY;

    // Determine the orientation of the second rectangle
    // to make sure that the edges which connect c0 and c1
    // are parallel to the edges which connect d0 and d1.
    diffX2 = d1 - d0;
    distX2 = diffX2.NormL2();
    diffX.Inner(diffX2, r);

    // Set the tolerance for determining if the edges
    // mentioned in the last comment are parallel.
    tol_r = distX2 < distX ? distX2 / 10 : distX / 10;

    // Bring the x- and y-coordinates of the intersection
    // into an order, where the smaller coordinates come
    // first.
    if(s[2] < s[0])
    {
      t[0] = s[2];
      t[2] = s[0];

      if(s[3] < s[1])
      {
        t[1] = s[3];
        t[3] = s[1];
        connect2[0] = ifaceElem2->connect[2];
        connect2[2] = ifaceElem2->connect[0];
        if (fabs(r) < tol_r) {
          connect2[1] = ifaceElem2->connect[1];
          connect2[3] = ifaceElem2->connect[3];
        }
        else {
          connect2[1] = ifaceElem2->connect[3];
          connect2[3] = ifaceElem2->connect[1];
        }
      }
      else
      {
        t[1] = s[1];
        t[3] = s[3];
        connect2[1] = ifaceElem2->connect[0];
        connect2[3] = ifaceElem2->connect[2];
        if (fabs(r) < tol_r) {
          connect2[0] = ifaceElem2->connect[3];
          connect2[2] = ifaceElem2->connect[1];
        }
        else {
          connect2[0] = ifaceElem2->connect[1];
          connect2[2] = ifaceElem2->connect[3];
        }
      }

    }
    else
    {
      t[0] = s[0];
      t[2] = s[2];

      if(s[3] < s[1])
      {
        t[1] = s[3];
        t[3] = s[1];
        connect2[1] = ifaceElem2->connect[2];
        connect2[3] = ifaceElem2->connect[0];
        if (fabs(r) < tol_r) {
          connect2[0] = ifaceElem2->connect[1];
          connect2[2] = ifaceElem2->connect[3];
        }
        else {
          connect2[0] = ifaceElem2->connect[3];
          connect2[2] = ifaceElem2->connect[1];
        }
      }
      else
      {
        t[1] = s[1];
        t[3] = s[3];
        connect2[0] = ifaceElem2->connect[0];
        connect2[2] = ifaceElem2->connect[2];
        if (fabs(r) < tol_r) {
          connect2[1] = ifaceElem2->connect[3];
          connect2[3] = ifaceElem2->connect[1];
        }
        else {
          connect2[1] = ifaceElem2->connect[1];
          connect2[3] = ifaceElem2->connect[3];
        }
      }
    }

    // Check if an intersection between rectangle1
    // and rectangle2 exists.
    if(t[0] >= 1.0)
      return false;

    if(t[2] <= 0.0)
      return false;

    if(t[1] >= 1.0)
      return false;

    if(t[3] <= 0.0)
      return false;

    shared_ptr<MortarNcSurfElem> ncElem(new MortarNcSurfElem());
    ncElem->connect.Resize(4);

    diffX *= distX;
    diffY *= distY;

    // If an intersection actually exist, we eventually
    // have to compute the intersection points.
    // There exist 16 different cases how two axiparallel
    // rectangles can intersect each other.

    Vector<Double> tmp;

    if(t[0] <= 0)
    {
      if(t[2] >= 1)
      {
        if(t[1] <= 0)
        {
          ncElem->connect[0] = ifaceElem1->connect[0];
          ncElem->connect[1] = ifaceElem1->connect[1];

          if(t[3] >= 1)
          {
            ncElem->connect[2] = ifaceElem1->connect[2];
            ncElem->connect[3] = ifaceElem1->connect[3];
          }
          else
          {
            tmp = c0 + diffX     + diffY*t[3];
            ptGrid_->AddNode(tmp, nodeNr);
            newNodes.Push_back(nodeNr);
            ncElem->connect[2] = nodeNr;
            tmp = c0 + diffX*0.0 + diffY*t[3];
            ptGrid_->AddNode(tmp, nodeNr);
            newNodes.Push_back(nodeNr);
            ncElem->connect[3] = nodeNr;
          }

        }
        else
        {
          if(t[3] >= 1)
          {
            tmp = c0 + diffX*0.0 + diffY*t[1];
            ptGrid_->AddNode(tmp, nodeNr);
            newNodes.Push_back(nodeNr);
            ncElem->connect[0] = nodeNr;
            tmp = c0 + diffX     + diffY*t[1];
            ptGrid_->AddNode(tmp, nodeNr);
            newNodes.Push_back(nodeNr);
            ncElem->connect[1] = nodeNr;
            ncElem->connect[2] = ifaceElem1->connect[2];
            ncElem->connect[3] = ifaceElem1->connect[3];
          }
          else
          {
            tmp = c0 + diffX*0.0 + diffY*t[1];
            ptGrid_->AddNode(tmp, nodeNr);
            newNodes.Push_back(nodeNr);
            ncElem->connect[0] = nodeNr;
            tmp = c0 + diffX     + diffY*t[1];
            ptGrid_->AddNode(tmp, nodeNr);
            newNodes.Push_back(nodeNr);
            ncElem->connect[1] = nodeNr;
            tmp = c0 + diffX     + diffY*t[3];
            ptGrid_->AddNode(tmp, nodeNr);
            newNodes.Push_back(nodeNr);
            ncElem->connect[2] = nodeNr;
            tmp = c0 + diffX*0.0 + diffY*t[3];
            ptGrid_->AddNode(tmp, nodeNr);
            newNodes.Push_back(nodeNr);
            ncElem->connect[3] = nodeNr;
          }
        }
      }
      else
      {
        if(t[1] <= 0)
        {
          if(t[3] >= 1)
          {
            ncElem->connect[0] = ifaceElem1->connect[0];
            tmp = c0 + diffX*t[2] + diffY*0.0;
            ptGrid_->AddNode(tmp, nodeNr);
            newNodes.Push_back(nodeNr);
            ncElem->connect[1] = nodeNr;
            tmp = c0 + diffX*t[2] + diffY;
            ptGrid_->AddNode(tmp, nodeNr);
            newNodes.Push_back(nodeNr);
            ncElem->connect[2] = nodeNr;
            ncElem->connect[3] = ifaceElem1->connect[3];
          }
          else
          {
            ncElem->connect[0] = ifaceElem1->connect[0];
            tmp = c0 + diffX*t[2] + diffY*0.0;
            ptGrid_->AddNode(tmp, nodeNr);
            newNodes.Push_back(nodeNr);
            ncElem->connect[1] = nodeNr;
            ncElem->connect[2] = connect2[2];
            tmp = c0 + diffX*0.0  + diffY*t[3];
            ptGrid_->AddNode(tmp, nodeNr);
            newNodes.Push_back(nodeNr);
            ncElem->connect[3] = nodeNr;
          }

        }
        else
        {
          if(t[3] >= 1)
          {
            tmp = c0 + diffX*0.0  + diffY*t[1];
            ptGrid_->AddNode(tmp, nodeNr);
            newNodes.Push_back(nodeNr);
            ncElem->connect[0] = nodeNr;
            ncElem->connect[1] = connect2[1];
            tmp = c0 + diffX*t[2] + diffY;
            ptGrid_->AddNode(tmp, nodeNr);
            newNodes.Push_back(nodeNr);
            ncElem->connect[2] = nodeNr;
            ncElem->connect[3] = ifaceElem1->connect[3];
          }
          else
          {
            tmp = c0 + diffX*0.0  + diffY*t[1];
            ptGrid_->AddNode(tmp, nodeNr);
            newNodes.Push_back(nodeNr);
            ncElem->connect[0] = nodeNr;
            ncElem->connect[1] = connect2[1];
            ncElem->connect[2] = connect2[2];
            tmp = c0 + diffX*0.0  + diffY*t[3];
            ptGrid_->AddNode(tmp, nodeNr);
            newNodes.Push_back(nodeNr);
            ncElem->connect[3] = nodeNr;
          }
        }
      }
    }
    else
    {
      if(t[2] >= 1)
      {
        if(t[1] <= 0)
        {
          if(t[3] >= 1)
          {
            tmp = c0 + diffX*t[0] + diffY*0.0;
            ptGrid_->AddNode(tmp, nodeNr);
            newNodes.Push_back(nodeNr);
            ncElem->connect[0] = nodeNr;
            ncElem->connect[1] = ifaceElem1->connect[1];
            ncElem->connect[2] = ifaceElem1->connect[2];
            tmp = c0 + diffX*t[0] + diffY;
            ptGrid_->AddNode(tmp, nodeNr);
            newNodes.Push_back(nodeNr);
            ncElem->connect[3] = nodeNr;
          }
          else
          {
            tmp = c0 + diffX*t[0] + diffY*0.0;
            ptGrid_->AddNode(tmp, nodeNr);
            newNodes.Push_back(nodeNr);
            ncElem->connect[0] = nodeNr;
            ncElem->connect[1] = ifaceElem1->connect[1];
            tmp = c0 + diffX      + diffY*t[3];
            ptGrid_->AddNode(tmp, nodeNr);
            newNodes.Push_back(nodeNr);
            ncElem->connect[2] = nodeNr;
            ncElem->connect[3] = connect2[3];
          }

        }
        else
        {
          if(t[3] >= 1)
          {
            ncElem->connect[0] = connect2[0];
            tmp = c0 + diffX      + diffY*t[1];
            ptGrid_->AddNode(tmp, nodeNr);
            newNodes.Push_back(nodeNr);
            ncElem->connect[1] = nodeNr;
            ncElem->connect[2] = ifaceElem1->connect[2];
            tmp = c0 + diffX*t[0] + diffY;
            ptGrid_->AddNode(tmp, nodeNr);
            newNodes.Push_back(nodeNr);
            ncElem->connect[3] = nodeNr;
          }
          else
          {
            ncElem->connect[0] = connect2[0];
            tmp = c0 + diffX      + diffY*t[1];
            ptGrid_->AddNode(tmp, nodeNr);
            newNodes.Push_back(nodeNr);
            ncElem->connect[1] = nodeNr;
            tmp = c0 + diffX      + diffY*t[3];
            ptGrid_->AddNode(tmp, nodeNr);
            newNodes.Push_back(nodeNr);
            ncElem->connect[2] = nodeNr;
            ncElem->connect[3] = connect2[3];
          }
        }
      }
      else
      {
        if(t[1] <= 0)
        {
          if(t[3] >= 1)
          {
            tmp = c0 + diffX*t[0] + diffY*0.0;
            ptGrid_->AddNode(tmp, nodeNr);
            newNodes.Push_back(nodeNr);
            ncElem->connect[0] = nodeNr;
            tmp = c0 + diffX*t[2] + diffY*0.0;
            ptGrid_->AddNode(tmp, nodeNr);
            newNodes.Push_back(nodeNr);
            ncElem->connect[1] = nodeNr;
            tmp = c0 + diffX*t[2] + diffY;
            ptGrid_->AddNode(tmp, nodeNr);
            newNodes.Push_back(nodeNr);
            ncElem->connect[2] = nodeNr;
            tmp = c0 + diffX*t[0] + diffY;
            ptGrid_->AddNode(tmp, nodeNr);
            newNodes.Push_back(nodeNr);
            ncElem->connect[3] = nodeNr;
          }
          else
          {
            tmp = c0 + diffX*t[0] + diffY*0.0;
            ptGrid_->AddNode(tmp, nodeNr);
            newNodes.Push_back(nodeNr);
            ncElem->connect[0] = nodeNr;
            tmp = c0 + diffX*t[2] + diffY*0.0;
            ptGrid_->AddNode(tmp, nodeNr);
            newNodes.Push_back(nodeNr);
            ncElem->connect[1] = nodeNr;
            ncElem->connect[2] = connect2[2];
            ncElem->connect[3] = connect2[3];
          }

        }
        else
        {
          if(t[3] >= 1)
          {
            ncElem->connect[0] = connect2[0];
            ncElem->connect[1] = connect2[1];
            tmp = c0 + diffX*t[2] + diffY;
            ptGrid_->AddNode(tmp, nodeNr);
            newNodes.Push_back(nodeNr);
            ncElem->connect[2] = nodeNr;
            tmp = c0 + diffX*t[0] + diffY;
            ptGrid_->AddNode(tmp, nodeNr);
            newNodes.Push_back(nodeNr);
            ncElem->connect[3] = nodeNr;
          }
          else
          {
            ncElem->connect[0] = connect2[0];
            ncElem->connect[1] = connect2[1];
            ncElem->connect[2] = connect2[2];
            ncElem->connect[3] = connect2[3];
          }
        }
      }
    }

    ncElem->type = Elem::ET_QUAD4;
    ncElem->ptPrimary = ifaceElem1;
    ncElem->ptSecondary = ifaceElem2;

    elemList_->AddElement(ncElem);

    return true;
  }

  void MortarInterface::GetInterfaceElemCoordinates(SurfElem* ifElem, StdVector<Vector<Double>>& coordinates) {
    // we still need this querry here as we need to treat quadratic elements like linear ones!
    UInt pSize = 0;
    switch(ifElem->type) {
      case Elem::ET_TRIA3:
      case Elem::ET_TRIA6:
        pSize = 3;
        break;
      case Elem::ET_QUAD4:
      case Elem::ET_QUAD8:
      case Elem::ET_QUAD9:
        pSize = 4;
        break;
      default:
        EXCEPTION("First argument to PolygonOnPolygon may not be of type '"
                << Elem::feType.ToString(ifElem->type) << "!");
    }
    coordinates.Resize(pSize);
  #ifdef USE_OPENMP
    omp_set_lock(&gridLock_);
  #endif
    for (UInt iNodes = 0; iNodes < pSize; ++iNodes) {
      ptGrid_->GetNodeCoordinate(coordinates[iNodes], ifElem->connect[iNodes], isMoving_);
    }
  #ifdef USE_OPENMP
    omp_unset_lock(&gridLock_);
  #endif
  }


  MortarInterface::LineIntersectType MortarInterface::CutLines(const Vector<Double> &a,
                                    const Vector<Double> &b, const Vector<Double> &c,
                                    const Vector<Double> &d, Vector<Double> &e) const
  {
    Double l1, l2;
    Vector<Double> v1, v2, temp;

  #ifdef CHECK_INDEX
    if ((a.GetSize() != 3) || (b.GetSize() != 3) || (c.GetSize() != 3) ||
        (d.GetSize() != 3)) {
      EXCEPTION("Points must be given as 3D coordinates");
      return INTERSECT_NONE;
    }
  #endif

    v1.Resize(3);
    v2.Resize(3);
    e.Resize(3);

    // calculate vectors of both lines
    v1 = b - a;
    v2 = d - c;
    // calculate lengths of both lines
    l1 = v1.NormL2();
    l2 = v2.NormL2();

    if ( v1.Collinear(v2) ) { // lines are parallel
      // if line from a to d is also parallel then lines may intersect
      e = d - a;
      if ( v1.Collinear(e) ) {
        Double l_ac, l_ad, l_bc, l_bd;
        // calculate distances between points
        temp = (c - a);
        l_ac = temp.NormL2();

        temp = (d - a);
        l_ad = temp.NormL2();

        temp = (c - b);
        l_bc = temp.NormL2();

        temp = (d - b);
        l_bd = temp.NormL2();

        // does a lie on [c,d]?
        if (fabs(l_ac + l_ad - l2) < absoluteTolerance_) {
          e = a;
          if (l_ac < absoluteTolerance_) // is a=c?
            return INTERSECT_A_EQ_C;
          if (l_ad < absoluteTolerance_) { // is a=d?
            // usually a cut at d wins over a cut at a,
            // but c might be an endpoint here, too
            if (fabs(l_ac + l_bc - l1) < absoluteTolerance_) {
              e = c;
              return INTERSECT_IN_C;
            }
            return INTERSECT_IN_D;
          }
          // does c lie on [a,b]?
          if (fabs(l_ac + l_bc - l1) < absoluteTolerance_)
            return INTERSECT_A_AND_C; // intersection is [a,c]

          return INTERSECT_IN_A;
        }
        if (fabs(l_bc + l_bd - l2) < absoluteTolerance_) {
          e = b;
          return INTERSECT_IN_B;
        }
        if (fabs(l_ac + l_bc - l1) < absoluteTolerance_) {
          e = c;
          return INTERSECT_IN_C;
        }
        if (fabs(l_ad + l_bd - l1) < absoluteTolerance_) {
          e = d;
          return INTERSECT_IN_D;
        }

        // both lines lie on one infinite virtual line
        if (l_bc < l_ac) { // we only consider [a,inf) for type 0
          // return point closest to a
          if (l_ad < l_ac)
            e = d;
          else
            e = c;
          return INTERSECT_OUTSIDE;
        }
        return INTERSECT_NONE; // line a->b does not point to line [c,d]
      }
      return INTERSECT_NONE; // lines are parallel but do not intersect
    }

    /* At this point we know that the lines are not parallel,
    * so compute intersection.
    *
    * a + h * v1 = c + k * v2
    *
    * This is a system with 2 unknowns (h,k) and 3 equations. Compute k1
    * from equations 1 and 2, k2 from equations 1 and 3, and k3 from
    * equations 2 and 3. Depending on the orientation of the lines in 3D
    * space, up to two values out of (k1,k2,k3) may be undefined, because
    * the denominator is zero. Therefore we need to select the right k.
    */

    Double h, k, k1 = 0.0, k2 = 0.0, k3 = 0.0, denom1, denom2, denom3;

    denom1 = v1[1] * v2[0] - v1[0] * v2[1];
    if (fabs(denom1) > absoluteTolerance_)
      k1 = (v1[0] * (c[1] - a[1]) + v1[1] * (a[0] - c[0])) / denom1;
    denom2 = v1[2] * v2[0] - v1[0] * v2[2];
    if (fabs(denom2) > absoluteTolerance_)
      k2 = (v1[0] * (c[2] - a[2]) + v1[2] * (a[0] - c[0])) / denom2;
    denom3 = v1[2] * v2[1] - v1[1] * v2[2];
    if (fabs(denom3) > absoluteTolerance_)
      k3 = (v1[1] * (c[2] - a[2]) + v1[2] * (a[1] - c[1])) /denom3;

    // If this system has no solution, lines do not intersect.
    if ((fabs(denom1) <= absoluteTolerance_)
        && (fabs(denom2) <= absoluteTolerance_)
        && (fabs(denom3) <= absoluteTolerance_))
      return INTERSECT_NONE;

    /* TODO: jens
    * This check makes no sense for 3 k's. Maybe add a check based on the
    * standard deviation of k.
    */
    /*if ((fabs(denom1) > absoluteTolerance_)
        && (fabs(denom2) > absoluteTolerance_)) {
      if (fabs(k1 - k2) > relativeTolerance_)
        return INTERSECT_NONE;
    }*/

    // select the right one out of (k1,k2,k3)
    if (fabs(denom1) > absoluteTolerance_)
      k = k1;
    else if (fabs(denom2) > absoluteTolerance_)
      k = k2;
    else
      k = k3;

    // compute second unknown
    if (fabs(v1[0]) > absoluteTolerance_)
      h = (c[0] - a[0] + v2[0] * k) / v1[0];
    else if (fabs(v1[1]) > absoluteTolerance_)
      h = (c[1] - a[1] + v2[1] * k) / v1[1];
    else
      h = (c[2] - a[2] + v2[2] * k) / v1[2];

    // compute point of intersection
    e = c + v2 * k; // do not use h, because it was computed from k

    if (h > -relativeTolerance_) { // we consider only [a,inf)
      if ((k > -relativeTolerance_) && (k < 1.0 + relativeTolerance_)) { // intersection on [c,d]?
        if (h < 1.0 + relativeTolerance_) { // intersection on [a,b]?
          // treat special cases
          if (fabs(h - 1.0) < relativeTolerance_) // h=1 means intersection in b
            return INTERSECT_IN_B;
          if (fabs(k - 1.0) < relativeTolerance_) // k=1 means intersection in d
            return INTERSECT_IN_D;
          if (fabs(k) < relativeTolerance_) { // k=0 means intersection in c
            if (fabs(h) < relativeTolerance_) // h=0 means intersection in a
              return INTERSECT_A_EQ_C;
            return INTERSECT_IN_C;
          }
          if (fabs(h) < relativeTolerance_) // h=0 means intersection in a
            return INTERSECT_IN_A;
          return INTERSECT_CROSS; // X intersection
        }
        return INTERSECT_ON_LINE2; // [a,inf) with [c,d]
      }
      return INTERSECT_OUTSIDE; // intersection not on any line
    }

    return INTERSECT_NONE; // no intersection (with [a,inf))
  }


  void MortarInterface::ProjectPolygon(StdVector<Vector<Double>>& projectionPolygonCoords,
                                       StdVector<Vector<Double>>& donorPolygonCoords,
                                       Vector<Double>& projectionNormVec,
                                       Vector<Double>& tempVec1,
                                       Vector<Double>& tempVec2) {
    Double projectionLength;
    // compute surface normal of the secondary element
    switch (donorPolygonCoords.GetSize()) {
      case 3:  //triangle
      case 6:  //triangle
      case 4:  //quadrilateral
      case 8:  //quadrilateral
      case 9:{ //quadrilateral
        tempVec1 = donorPolygonCoords[1] - donorPolygonCoords[0];
        tempVec2 = donorPolygonCoords[2] - donorPolygonCoords[0];
        tempVec1.CrossProduct(tempVec2, projectionNormVec);
        projectionNormVec.Normalize();
        break;
      }
      default:
        EXCEPTION("Unsupported element type of donor element!");
    }
    // print a warning if the secondary element is a non-planar quadrilateral.
    // In this case, the projection is not well-defined.
    // We only need to check for the initial run of the function.
    if (geoWarn_ && callCounter_ == 0) {
      switch (donorPolygonCoords.GetSize()) {
        case 4:
        case 8:
        case 9: {
          Vector<Double> p1, p2;
          Vector<Double> uniquePlaneNormVec;
          // in the case of quadrilaterals, we can find a unique plane by taking the cross product of the lines connecting the centers 
          // two opposing lines of the element
          p1 = donorPolygonCoords[1] - donorPolygonCoords[0];
          p1.ScalarDiv(2.0);
          p1 += donorPolygonCoords[0];
          p2 = donorPolygonCoords[3] - donorPolygonCoords[2];
          p2.ScalarDiv(2.0);
          p2 += donorPolygonCoords[2];
          tempVec1 = p2 - p1;
          p1 = donorPolygonCoords[2] - donorPolygonCoords[1];
          p1.ScalarDiv(2.0);
          p1 += donorPolygonCoords[1];
          p2 = donorPolygonCoords[0] - donorPolygonCoords[3];
          p2.ScalarDiv(2.0);
          p2 += donorPolygonCoords[3];
          tempVec2 = p2 - p1;
          tempVec1.CrossProduct(tempVec2, uniquePlaneNormVec);
          // tempVec1 is the normal vector of the unique plane
          uniquePlaneNormVec.Normalize();
          // check if the normals coincide
          uniquePlaneNormVec.CrossProduct(projectionNormVec, tempVec1);
          if (tempVec1.NormL2() > NORM_EPS)
            WARN("The donor element in point "
                 << donorPolygonCoords << " on region " << secondarySurfRegion_
                 << " is a non-planar quadrilateral. The projection is not "
                    "well-defined. Switching primary/secondary side could "
                    "help.");
          break;
        }
        default:
          break;
      }
    }
    // project each node of the primary element
    for (UInt iNodes = 0; iNodes < projectionPolygonCoords.GetSize(); ++iNodes) {
      // compute the projection distance from a node to the projection surface
      tempVec1 = projectionPolygonCoords[iNodes] - donorPolygonCoords[0];
      projectionLength = projectionNormVec.Inner(tempVec1);
      // modify the coordinates
      projectionPolygonCoords[iNodes] -= projectionNormVec * projectionLength;
    }
  }


  bool MortarInterface::GetRotationMatrix(Matrix<Double>& rotationMatrix,
                                          Matrix<Double>& rotationMatrixTrans,
                                          Vector<Double>& normVec,
                                          Vector<Double>& rotationAxis,
                                          Vector<Double>& dirVec) {
    // get the rotation axis 
    normVec.CrossProduct(dirVec, rotationAxis);
    // if the vectors coincide we do not need to rotate
    if (rotationAxis.NormL2() == 0)
      return false;
    // normalize the rotation axis vector
    rotationAxis.Normalize();
    // ca = cos(theta) = cos(n ^ ez) = n_z/|n| = n_z
    Double ca = normVec[2];
    // sa = sin(theta) = sin(n ^ ez) = sqrt(n_x^2 + n_y^2)/|n| = sqrt(n_x^2 + n_y^2)
    Double sa = sqrt(normVec[0] * normVec[0] + normVec[1] * normVec[1]);
    // ci = 1-cos(theta) = 1 - cos(n ^ ez)
    Double ci = 1 - ca;
    // rotation matrix
    rotationMatrix[0][0] = ca + rotationAxis[0] * rotationAxis[0] * ci;
    rotationMatrix[1][0] = rotationAxis[0] * rotationAxis[1] * ci + rotationAxis[2] * sa;
    rotationMatrix[2][0] = rotationAxis[0] * rotationAxis[2] * ci - rotationAxis[1] * sa;
    rotationMatrix[0][1] = rotationAxis[0] * rotationAxis[1] * ci - rotationAxis[2] * sa;
    rotationMatrix[1][1] = ca + rotationAxis[1] * rotationAxis[1] * ci;
    rotationMatrix[2][1] = rotationAxis[1] * rotationAxis[2] * ci + rotationAxis[0] * sa;
    rotationMatrix[0][2] = rotationAxis[0] * rotationAxis[2] * ci + rotationAxis[1] * sa;
    rotationMatrix[1][2] = rotationAxis[1] * rotationAxis[2] * ci - rotationAxis[0] * sa;
    rotationMatrix[2][2] = ca + rotationAxis[2] * rotationAxis[2] * ci;
    // store the transposed rotation matrix for the inverse rotation afterwards.
    rotationMatrix.Transpose(rotationMatrixTrans);
    return true;
  }


#ifdef USE_CGAL
  void MortarInterface::IntersectCoplanarPolygons(StdVector<UInt>& newNodeIds) {
    // these variables will be shared
    bool rotatePolys;
    Vector<Double> normVec;
    Vector<Double> rotationAxis;
    Matrix<Double> rotationMatrix;
    Matrix<Double> rotationMatrixTrans;
    // these variables will be private
    UInt primaryIndex, secondaryIndex;
    StdVector<Vector<Double>> primaryNodeCoords;
    StdVector<Vector<Double>> secondaryNodeCoords;
    Vector<Double> tempVec1;
    Vector<Double> tempVec2;
    // resize
    normVec.Resize(ptGrid_->GetDim());
    rotationAxis.Resize(ptGrid_->GetDim());
    tempVec1.Resize(ptGrid_->GetDim());
    tempVec2.Resize(ptGrid_->GetDim());
    rotationMatrix.Resize(ptGrid_->GetDim(),ptGrid_->GetDim());
    rotationMatrixTrans.Resize(ptGrid_->GetDim(),ptGrid_->GetDim());
    // get the coordinates and the norm vector of the first secondary element
    secondaryIndex = intersectionCandiatesIdx_[0].second;
    GetInterfaceElemCoordinates(secondarySurfElems_[secondaryIndex], secondaryNodeCoords);
    // compute surface normal of the secondary element
    tempVec1 = secondaryNodeCoords[1] - secondaryNodeCoords[0];
    tempVec2 = secondaryNodeCoords[2] - secondaryNodeCoords[0];
    tempVec1.CrossProduct(tempVec2, normVec);
    normVec.Normalize();
    // get the z-unitary vector
    tempVec1[0] = 0.0;
    tempVec1[1] = 0.0;
    tempVec1[2] = 1.0;
    // compute the rotation matrix only once
    rotatePolys = GetRotationMatrix(rotationMatrix, rotationMatrixTrans, normVec, rotationAxis, tempVec1);
    // parallelize
    #pragma omp parallel num_threads(CFS_NUM_THREADS) private(primaryIndex, secondaryIndex, primaryNodeCoords, secondaryNodeCoords, tempVec1, tempVec2)
    {
      // additional private variables
      bool doIntersect;
      UInt numElemsStart, numElemsUpdated;
      StdVector<UInt> tempNewNodeIds;
      StdVector<Vector<Double>> intersectionPolygonCoords;
      StdVector<Vector<Double>> rotatedPrimaryNodeCoords;
      StdVector<Vector<Double>> rotatedSecondaryNodeCoords;
      // additional resize of private copies
      tempVec1.Resize(ptGrid_->GetDim());
      tempVec2.Resize(ptGrid_->GetDim());
      #pragma omp for
      for (Integer iCandidate = 0; iCandidate < (Integer)intersectionCandiatesIdx_.size(); ++iCandidate) {
        primaryIndex = intersectionCandiatesIdx_[iCandidate].first;
        secondaryIndex = intersectionCandiatesIdx_[iCandidate].second;
        // get node coordinates of both elements
        GetInterfaceElemCoordinates(primarySurfElems_[primaryIndex], primaryNodeCoords);
        GetInterfaceElemCoordinates(secondarySurfElems_[secondaryIndex], secondaryNodeCoords);
        // compute the intersection element using the method depending on the software build
        doIntersect = CGAL_Cut2DPolygons(primaryNodeCoords, secondaryNodeCoords, intersectionPolygonCoords, tempVec1, tempVec2, 
                                         rotationAxis, rotationMatrix, rotationMatrixTrans, rotatedPrimaryNodeCoords, rotatedSecondaryNodeCoords, rotatePolys);
        // triangulate the intersection polygon and add surface elements to the grid...
        if (!doIntersect) {
          // without intersection we move to the next combination
          continue;
        } else {
          // add the new polygon to the grid...
          StdVector<MortarNcSurfElem*> triangulatedElements;
          // check how many nodes and elements are already added to the current interface
          numElemsStart = triangulatedElements.GetSize();
          // triangulate the intersection polygon
          TriangulatePoly(intersectionPolygonCoords, triangulatedElements, tempNewNodeIds);
          // TriangulatePoly() might have added nodes and elements to the grid, so
          // get the updated numbers after triangulation
          numElemsUpdated = triangulatedElements.GetSize();
          // now loop over all triangulated elements and add to the NC elemList_
          for (UInt iElems = numElemsStart; iElems < numElemsUpdated; ++iElems) {
            shared_ptr<MortarNcSurfElem> elemToAdd(triangulatedElements[iElems]);
            elemToAdd->ptPrimary = primarySurfElems_[primaryIndex];
            elemToAdd->ptSecondary = secondarySurfElems_[secondaryIndex];
            elemToAdd->projectedPrimary = nullptr;
            elemToAdd->transVect = (mutualProjection_ || cakePieceProjection_) ? &GetTranslationVector() : nullptr;
          #ifdef USE_OPENMP
            // set lock to add elements to the grid
            omp_set_lock(&elemListLock_);
          #endif
            // add elements to the elem list
            elemList_->AddElement(elemToAdd);
          #ifdef USE_OPENMP
            omp_unset_lock(&elemListLock_);
          #endif
          }
        } // (!doIntersect)
      }
    #ifdef USE_OPENMP
      // add computed nodes to the general node vector
      omp_set_lock(&newNodesLock_);
    #endif
      for (UInt iNodes = 0; iNodes < tempNewNodeIds.GetSize(); ++iNodes) {
        newNodeIds.Push_back(tempNewNodeIds[iNodes]);
      }
    #ifdef USE_OPENMP
      omp_unset_lock(&newNodesLock_);
    #endif
    }
  }


  void MortarInterface::IntersectGeneralPolygons(StdVector<UInt>& newNodeIds) {
    // parallelize
    #pragma omp parallel num_threads(CFS_NUM_THREADS)
    {
      // declare private variables
      bool doIntersect, rotatePolys;
      StdVector<UInt> tempNewNodeIds;
      StdVector<Vector<Double>> intersectionPolygonCoords;
      StdVector<Vector<Double>> primaryNodeCoords;
      StdVector<Vector<Double>> secondaryNodeCoords;
      Vector<Double> projectionNormVec;
      Vector<Double> tempVec1;
      Vector<Double> tempVec2;
      Vector<Double> rotationAxis;
      Matrix<Double> rotationMatrix;
      Matrix<Double> rotationMatrixTrans;
      StdVector<Vector<Double>> rotatedPrimaryNodeCoords;
      StdVector<Vector<Double>> rotatedSecondaryNodeCoords;
      UInt primaryIndex, secondaryIndex;
      // resize
      projectionNormVec.Resize(ptGrid_->GetDim());
      rotationAxis.Resize(ptGrid_->GetDim());
      tempVec1.Resize(ptGrid_->GetDim());
      tempVec2.Resize(ptGrid_->GetDim());
      rotationMatrix.Resize(ptGrid_->GetDim(),ptGrid_->GetDim());
      rotationMatrixTrans.Resize(ptGrid_->GetDim(),ptGrid_->GetDim());
      #pragma omp for
      for (Integer iCandidate = 0; iCandidate < (Integer)intersectionCandiatesIdx_.size(); ++iCandidate) {
        primaryIndex = intersectionCandiatesIdx_[iCandidate].first;
        secondaryIndex = intersectionCandiatesIdx_[iCandidate].second;
        // get node coordinates of both elements
        GetInterfaceElemCoordinates(primarySurfElems_[primaryIndex], primaryNodeCoords);
        GetInterfaceElemCoordinates(secondarySurfElems_[secondaryIndex], secondaryNodeCoords);
        // project the primary element onto the secondary-element plane
        // note: afer calling this function, the primaryNodeCoords will not be the original, 
        // but the projected coordinates
        ProjectPolygon(primaryNodeCoords, secondaryNodeCoords, projectionNormVec, tempVec1, tempVec2);
        // we need to rotate the elements of both sides into the (x,y) plane 
        // so that we are able to use cgal's 2D intersection methods...
        // get the z-unitary vector
        tempVec1[0] = 0.0;
        tempVec1[1] = 0.0;
        tempVec1[2] = 1.0;
        // compute the roation matrix to use cgal 2d algorithms
        rotatePolys = GetRotationMatrix(rotationMatrix, rotationMatrixTrans, projectionNormVec, rotationAxis, tempVec1);
        // compute the intersection element using the method depending on the software build
        doIntersect = CGAL_Cut2DPolygons(primaryNodeCoords, secondaryNodeCoords, intersectionPolygonCoords, tempVec1, tempVec2, 
                                         rotationAxis, rotationMatrix, rotationMatrixTrans, rotatedPrimaryNodeCoords, rotatedSecondaryNodeCoords, rotatePolys);
        // triangulate the intersection polygon and add surface elements to the grid...
        if (!doIntersect) {
          // without intersection we move to the next combination
          LOG_DBG3(mortarInterface) << "No intersection element stored for:" << std::endl << "Primary element: " << primaryIndex << std::endl << "Secondary element: " << secondaryIndex << std::endl;
          continue;
        } else {
          // add the new polygon to the grid...
          StdVector<MortarNcSurfElem*> triangulatedElements;
          // check how many elements are already added to the current interface
          UInt numElemsStart = triangulatedElements.GetSize();
          // triangulate the intersection polygon
          TriangulatePoly(intersectionPolygonCoords, triangulatedElements, tempNewNodeIds);
          // TriangulatePoly() might have added elements to the grid, so
          // get the updated numbers after triangulation
          UInt numElemsUpdated = triangulatedElements.GetSize();
          // pointer to projected elements
          shared_ptr<SurfElem> projectedPrimaryElement = nullptr;
          // if the interface is not coplanar we need to do additional work
          // variable to loop over nodes of the intersection polygon
          Vector<Double> currNodeCoord;
          // create the projected elements on the heap to avoid another global-local transformation in the ABInt
          projectedPrimaryElement.reset(new SurfElem());
          projectedPrimaryElement->type = primarySurfElems_[primaryIndex]->type;
          projectedPrimaryElement->connect.Resize(primaryNodeCoords.GetSize());
          // Assign connectivity list to the projected elements.
          for (UInt iPrimaryNodes = 0; iPrimaryNodes < primaryNodeCoords.GetSize(); ++iPrimaryNodes) {
            // in the non-coplanar case, ProjectPolygon() has already projected the primaryNodeCoords.
            // so assign the new projected nodes to the grid and to the projected primary element.
            // AddNodeToGrid directly assigns the new node id to the projected element in this case.
            AddNodeToGrid(primaryNodeCoords[iPrimaryNodes], projectedPrimaryElement->connect[iPrimaryNodes], tempNewNodeIds);
          }
          // now loop over all triangulated elements and add to the NC elemList_
          for (UInt iElems = numElemsStart; iElems < numElemsUpdated; ++iElems) {
             // we only need the projected primary element in the non-coplanar case
            shared_ptr<MortarNcSurfElem> elemToAdd(triangulatedElements[iElems]);
            elemToAdd->ptPrimary = primarySurfElems_[primaryIndex];
            elemToAdd->ptSecondary = secondarySurfElems_[secondaryIndex];
            elemToAdd->projectedPrimary = projectedPrimaryElement;
            elemToAdd->transVect = (mutualProjection_ || cakePieceProjection_) ? &GetTranslationVector() : nullptr;
          #ifdef USE_OPENMP
            // set lock to add elements to the grid
            omp_set_lock(&elemListLock_);
          #endif
            // add elements to the elem list
            elemList_->AddElement(elemToAdd);
          #ifdef USE_OPENMP
            omp_unset_lock(&elemListLock_);
          #endif
          }
        } // (!doIntersect)
      }
    #ifdef USE_OPENMP
      // add computed nodes to the general node vector
      omp_set_lock(&newNodesLock_);
    #endif
      for (UInt iNodes = 0; iNodes < tempNewNodeIds.GetSize(); ++iNodes) {
        newNodeIds.Push_back(tempNewNodeIds[iNodes]);
      }
    #ifdef USE_OPENMP
      omp_unset_lock(&newNodesLock_);
    #endif
    }
  }


  bool MortarInterface::CGAL_Cut2DPolygons(StdVector<Vector<Double>>& primaryNodeCoords,
                                           StdVector<Vector<Double>>& secondaryNodeCoords,
                                           StdVector<Vector<Double>>& intersectionPolygonCoords,
                                           Vector<Double>& tempVec1,
                                           Vector<Double>& tempVec2,
                                           Vector<Double>& rotationAxis,
                                           Matrix<Double>& rotationMatrix,
                                           Matrix<Double>& rotationMatrixTrans,
                                           StdVector<Vector<Double>>& rotatedPrimaryNodeCoords,
                                           StdVector<Vector<Double>>& rotatedSecondaryNodeCoords,
                                           const bool& rotatePolys) {
    // iterate over nodes
    UInt iNodes;
    // reset the intersection polygon
    intersectionPolygonCoords.Clear(true);
    // the cgal polygons to use CGAL 2d intersection algorithms
    CGALPolygon2 cgalPolyPrim, cgalPolySec;
    // we need to rotate the elements of both sides into the (x,y) plane 
    // so that we are able to use cgal's 2D intersection methods...
    // if the secondary element's normal already coincides with the z-unitary vector, we can directly compute
    // the intersection using cgal's 2d algorithms. If not, we need to rotate the elements first into the (x,y) plane.
    if (!rotatePolys) {
      for (iNodes = 0; iNodes < primaryNodeCoords.GetSize(); ++iNodes) {
        cgalPolyPrim.push_back(CGALPoint2(primaryNodeCoords[iNodes][0], primaryNodeCoords[iNodes][1]));
      }
      for (iNodes = 0; iNodes < secondaryNodeCoords.GetSize(); ++iNodes) {
        cgalPolySec.push_back(CGALPoint2(secondaryNodeCoords[iNodes][0], secondaryNodeCoords[iNodes][1]));
      }
    } else {
      // Perform rotations and create CGAL 2D polygons.
      // We do not change the initial polygons, because, for
      // example, p1 is used further to create a projected primary element.
      rotatedPrimaryNodeCoords.Resize(primaryNodeCoords.GetSize());
      rotatedSecondaryNodeCoords.Resize(secondaryNodeCoords.GetSize());
      for (iNodes = 0; iNodes < primaryNodeCoords.GetSize(); ++iNodes) {
        rotatedPrimaryNodeCoords[iNodes] = rotationMatrix * primaryNodeCoords[iNodes];
        cgalPolyPrim.push_back(CGALPoint2(rotatedPrimaryNodeCoords[iNodes][0], rotatedPrimaryNodeCoords[iNodes][1]));
      }
      for (iNodes = 0; iNodes < secondaryNodeCoords.GetSize(); ++iNodes) {
        rotatedSecondaryNodeCoords[iNodes] = rotationMatrix * secondaryNodeCoords[iNodes];
        cgalPolySec.push_back(CGALPoint2(rotatedSecondaryNodeCoords[iNodes][0], rotatedSecondaryNodeCoords[iNodes][1]));
      }
    } // (!rotatePolys)
    // check orientation of polygons
    if (cgalPolyPrim.is_clockwise_oriented())
      cgalPolyPrim.reverse_orientation();
    if (cgalPolySec.is_clockwise_oriented())
      cgalPolySec.reverse_orientation();

    // compute intersections
    std::vector<CGALPolygonWithHoles2> intersectionList;
    CGAL::intersection(cgalPolyPrim, cgalPolySec, std::back_inserter(intersectionList));

    // do some sanity checks...
    // If the intersection list has the size = 0, the polygons don't
    // intersect with each other. If the size is > 1, then something
    // must be wrong, because the intersection two convex sets is a convex set.
    if (intersectionList.size() == 0) {
      LOG_DBG3(mortarInterface) << "No intersection found on element.";
      return false;
    }
    else if (intersectionList.size() > 1) {
      // here should be an exception as this is definitely not leading to valid simulation
      EXCEPTION("More than one intersection polygon found for elements at (" 
                << primaryNodeCoords << ") and (" << secondaryNodeCoords  << std::endl 
                <<  ") in non-conforming interface '" << name_ << "'." << std::endl 
                << "Something must be wrong with your mesh!");
      return false;
    }

    // We need to transform the vertices of the CGAL-polygon, back to a CFS-Vector
    for (CGALPolygon2::Edge_const_iterator edgeIt = intersectionList[0].outer_boundary().edges_begin();
        edgeIt != intersectionList[0].outer_boundary().edges_end(); ++edgeIt) {
      // In some cases the polygon can contain a negligibly short edges which
      // must be omitted. Therefore, we iterate over the edges, check their
      // lengths, and take the starting points of those having considerable lengths.
      if (edgeIt->squared_length() < pow(coincidenceRadius_,2)) {
        LOG_DBG3(mortarInterface) << "Omitting coincident intersection nodes.";
        continue;
      }
      // We omit the intersection polygon if its relative area is too little
      if (edgeIt->squared_length() < pow(minRelativeSideLength_,2) * min(cgalPolyPrim.area(),cgalPolySec.area())) {
        LOG_DBG3(mortarInterface) << "Rejecting intersection polygon due to small length detection.";
        return false;
      }
      // tempVec1 - a vertex of the intersection polygon; tempVec2 - a vertex rotated back
      tempVec1[0] = CGAL::to_double(edgeIt->vertex(0).x());
      tempVec1[1] = CGAL::to_double(edgeIt->vertex(0).y());
      if (rotatePolys) {
        // Restore the z-coordinate of the rotated polygons.
        // The cgal intersection polygon is parallel to the (x,y) plane, so we can choose any 
        // z-coordinate. The only condition - it must be the same for all vertices.
        tempVec1[2] = rotatedPrimaryNodeCoords[0][2];
        // Perform the back-rotation, into the original plane...
        // It is required in order to get correct results transforming Local-to-Global and back
        // within SurfaceMortarABInt::CalcElementMatrix method.
        tempVec2 = rotationMatrixTrans * tempVec1;
        intersectionPolygonCoords.Push_back(tempVec2);
      } else {
        // if we did not rotate we can directly assign the z coordinate of the projected primary element
        tempVec1[2] = primaryNodeCoords[0][2];
        // and we can directly assign the (x,y) coords of the intersection vertices
        intersectionPolygonCoords.Push_back(tempVec1);
      }
    }
    // feed back that we succeeded, if we have at least 3 intersection points. Otherwise omit the element
    if (intersectionPolygonCoords.GetSize() > 2) {
      return true;
    } else {
      LOG_DBG3(mortarInterface) << "Rejecting intersection polygon as less than 3 nodes are left.";
      return false;
    }
  }
#else
  bool MortarInterface::IntersectPolygons(SurfElem* primaryElement,
                                          SurfElem* secondaryElement,
                                          StdVector<UInt>& newNodeIds,
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
                                          StdVector<Vector<Double>>& rotPolySec) {
    // get node coordinates of both elements
    GetInterfaceElemCoordinates(primaryElement, primaryNodeCoords);
    GetInterfaceElemCoordinates(secondaryElement, secondaryNodeCoords);
    UInt numPrimaryNodes = primaryNodeCoords.GetSize();
    // boolean to check if an intersection element has been computed
    bool intersectionComputed;
    // compute the intersection element using the method depending on the software build
    intersectionComputed = CutPolygon(primaryNodeCoords, secondaryNodeCoords, intersectionPolygonCoords, tempVec1, tempVec2, tempVec3, tempVec4, tempVec5);
    // triangulate the intersection polygon and add surface elements to the grid...
    if (!intersectionComputed) {
      // without intersection, return
      return false;
    } else {
      // add the new polygon to the grid...
      StdVector<MortarNcSurfElem*> triangulatedElements;
      // check how many nodes and elements are already added to the current interface
      UInt numElemsStart = triangulatedElements.GetSize();
      // triangulate the intersection polygon
      TriangulatePoly(intersectionPolygonCoords, triangulatedElements, newNodeIds);
      // TriangulatePoly() might have added nodes and elements to the grid, so
      // get the updated numbers after triangulation
      UInt numElemsUpdated = triangulatedElements.GetSize();

      // pointer to projected primary element
      shared_ptr<SurfElem> projectedPrimaryElement = nullptr;
      // if the interface is not coplanar we need to do additional work
      if (!isCoplanar_) {
        // variable to loop over nodes of the intersection polygon
        Vector<Double> currNodeCoord;
        // create the projected primary element on the heap
        projectedPrimaryElement.reset(new SurfElem());
        projectedPrimaryElement->type = primaryElement->type;
        projectedPrimaryElement->connect.Resize(numPrimaryNodes);
        // Assign connectivity list to the projected primary element.
        for (UInt iPrimaryNodes = 0; iPrimaryNodes < numPrimaryNodes; ++iPrimaryNodes) {
          // in the non-coplanar case, CutPolygon() has already projected the primaryNodeCoords.
          // so assign the new projected nodes to the grid and to the projected primary element.
          // AddNodeToGrid directly assigns the new node id to the projected element in this case.
          AddNodeToGrid(primaryNodeCoords[iPrimaryNodes], projectedPrimaryElement->connect[iPrimaryNodes], newNodeIds);
        }
      } // (!isCoplanar_)
      // now loop over all triangulated elements and add to the NC elemList_
      for (UInt iElems = numElemsStart; iElems < numElemsUpdated; ++iElems) {
        shared_ptr<MortarNcSurfElem> elemToAdd(triangulatedElements[iElems]);
        elemToAdd->ptPrimary = primaryElement;
        elemToAdd->ptSecondary = secondaryElement;
        elemToAdd->projectedPrimary = projectedPrimaryElement;
        elemToAdd->transVect = (mutualProjection_ || cakePieceProjection_) ? &GetTranslationVector() : nullptr;
      #ifdef USE_OPENMP
        // set lock to add elements to the grid
        omp_set_lock(&elemListLock_);
      #endif
        // add elements to the elem list
        elemList_->AddElement(elemToAdd);
      #ifdef USE_OPENMP
        omp_unset_lock(&elemListLock_);
      #endif
      }
      // feed back that we succeeded
      return true;
    }
  }


  bool MortarInterface::CutPolygon(StdVector<Vector<Double>>& primaryNodeCoords,
                                   StdVector<Vector<Double>>& secondaryNodeCoords,
                                   StdVector<Vector<Double>>& intersectionPolygonCoords,
                                   Vector<Double>& primaryCentroid,
                                   Vector<Double>& secondaryCentroid,
                                   Vector<Double>& tempVec1,
                                   Vector<Double>& tempVec2,
                                   Vector<Double>& intersectionPointCoords) {
    // declare variables...
    // circumcircles of elements
    Double primaryElementCircle;
    Double secondaryElementCircle;
    // number of nodes in the elements
    UInt numPrimaryNodes = primaryNodeCoords.GetSize();
    UInt numSecondaryNodes = secondaryNodeCoords.GetSize();
    // start position of the polygon iterator
    UInt startPolyIter = numPrimaryNodes;
    // number of primary points lying in the secondary element
    UInt numContainedPoints = 0;
    // normal vectors of the (triangle) element
    Vector<Double> primaryNormVec, secondaryNormVec;
    // reset the intersection polygon
    intersectionPolygonCoords.Clear(true);
    // count the number of line cuts between the polygons
    UInt nCuts = 0;
    // counter to loop over nodes
    UInt iNodes = 0;

    // data type to store intersection info
    struct Intersection {
      UInt index;         // position of the polygon iterator
      UInt type;          // intersection type LineIntersectType
      bool swap;          // indicate if elements must be swapped for the intersection
      Vector<Double> loc; // coordinates of the intersection point;
    } cuts[2];

    // compute surrounding circles of both polygons
    PolyCentroid(primaryNodeCoords, primaryCentroid);
    primaryElementCircle = PolyCircumcircle(primaryNodeCoords,primaryCentroid);
    PolyCentroid(secondaryNodeCoords, secondaryCentroid);
    secondaryElementCircle = PolyCircumcircle(secondaryNodeCoords, secondaryCentroid);

    // quit, if surrounding circles do not intersect
    tempVec1 = primaryCentroid - secondaryCentroid;
    if ((primaryElementCircle + secondaryElementCircle) < sqrt(tempVec1*tempVec1))
      return false;

    // compute surface normal of the secondary element
    tempVec1 = secondaryNodeCoords[1] - secondaryNodeCoords[0];
    tempVec2 = secondaryNodeCoords[2] - secondaryNodeCoords[0];
    tempVec1.CrossProduct(tempVec2, secondaryNormVec);
    secondaryNormVec.Normalize();

    // if interface is not coplanar then project the primary element onto the secondary element
    if (!isCoplanar_) {
      // distance between a node and the projection plane
      Double projectionLength;
      // project each node of the primary element
      for (iNodes = 0; iNodes < numPrimaryNodes; ++iNodes) {
        // compute the projection distance from a node to the projection surface
        tempVec1 = primaryNodeCoords[iNodes] - secondaryNodeCoords[0];
        projectionLength = secondaryNormVec.Inner(tempVec1);
        // modify the coordinates
        primaryNodeCoords[iNodes] -= secondaryNormVec * projectionLength;
      }
    }

    // Count those points of p1 that are contained in p2. Choose a point
    // that lies outside of p2 as starting point.
    for (iNodes = 0; iNodes < numPrimaryNodes; ++iNodes) {
      if (PointInsidePoly(primaryNodeCoords[iNodes], secondaryNodeCoords, &secondaryCentroid))
        ++numContainedPoints;
      // set start position for the polygon iterator
      else if ((numContainedPoints == 0) || (startPolyIter == numPrimaryNodes))
        startPolyIter = iNodes;

    }
    // if the primary element is entirely contained in the secondary element,
    // we can directly set the intersection polygon as the primary element and return
    // (for convex polygons)
    if (numContainedPoints == numPrimaryNodes) {
      intersectionPolygonCoords = primaryNodeCoords;
      return true;
    }
    // no points are contained, check the opposite way
    if (numContainedPoints == 0) {
      for (iNodes = 0; iNodes < numSecondaryNodes; ++iNodes) {
        if (PointInsidePoly(secondaryNodeCoords[iNodes], primaryNodeCoords, &primaryCentroid))
          ++numContainedPoints;
      }
      // if the secondary element is entirely contained in the primary element,
      // directly assign it to the output
      if (numContainedPoints == numSecondaryNodes) {
        intersectionPolygonCoords = secondaryNodeCoords;
        return true;
      }
    }

    // WARNING: One can not conclude that two polygons do not intersect from
    // the fact that no point lies inside the other polygon.

    // create polygon iterators and make to sure that both polygons have the same orientation
    PolygonIterator primaryIterator(primaryNodeCoords, startPolyIter);
    PolygonIterator secondaryIterator(secondaryNodeCoords);
    tempVec1 = primaryNodeCoords[1] - primaryNodeCoords[0];
    tempVec2 = primaryNodeCoords[2] - primaryNodeCoords[0];
    tempVec1.CrossProduct(tempVec2, primaryNormVec);

    if (primaryNormVec * secondaryNormVec < 0.0)
      secondaryIterator.Reverse();

    // find the first cut of two edges of the polygons...
    // iterate over primary element
    do {
      // iterate over secondaty element
      do {
        // check for the intersection type
        LineIntersectType cutType = CutLines(*primaryIterator, primaryIterator.Next(),
                                            *secondaryIterator, secondaryIterator.Next(),
                                            intersectionPointCoords);
        Intersection cut = {secondaryIterator.GetPos(), cutType, false, intersectionPointCoords};
        // See what kind of cut we have found.
        // This section is different from the main loop, because we do
        // not know if we cut from outside into p2 or vice versa. This
        // can happen despite the starting point lying outside, because
        // an edge of p1 might cut p2 into halves.
        switch (cutType) {
          case INTERSECT_CROSS: {
            // always store the cut for cross intersection
            break;
          }
          case INTERSECT_IN_A: {
            // see if [a,b] lies inside of the secondary element or not.
            if ((CutLines(*primaryIterator, primaryIterator.Next(),
                          *secondaryIterator, secondaryIterator.Next(2),
                          intersectionPointCoords) >= INTERSECT_ON_LINE2) ||
                (CutLines(*primaryIterator, primaryIterator.Next(),
                          secondaryIterator.Prev(), secondaryIterator.Next(),
                          intersectionPointCoords) >= INTERSECT_ON_LINE2)) {
              break;
            }
            continue; // break without storing as [a,b] lies outside of p2
          }
          // do the same for INTERSECT_IN_C and INTERSECT_A_EQ_C
          case INTERSECT_IN_C:
          case INTERSECT_A_EQ_C: {
            // if [a,b] cuts into the secondary element, store the cut
            if (CutLines(*primaryIterator, primaryIterator.Next(),
                          secondaryIterator.Prev(), secondaryIterator.Next(),
                          intersectionPointCoords) >= INTERSECT_ON_LINE2) {
              break;
            }
            // if [c,d] lies inside of the primary element, swap the element and store cut
            if (CutLines(*secondaryIterator, secondaryIterator.Next(),
                          primaryIterator.Prev(), primaryIterator.Next(),
                          intersectionPointCoords) >= INTERSECT_ON_LINE2) {
              cut.swap = true;
              break;
            }
            // if line a is identical to line c, continue as polygons touch in c only
            if (cutType != INTERSECT_A_EQ_C) {
              if (CutLines(*secondaryIterator, secondaryIterator.Next(),
                          *primaryIterator, primaryIterator.Next(2),
                          intersectionPointCoords) >= INTERSECT_ON_LINE2) {
                cut.swap = true;
                break;
              }
            }
            continue; // break without storing as polygons touch in c only
          }
          case INTERSECT_A_AND_C: {
            // if we have an intersection through a and c, assign info to the intersection struct
            nCuts = 2;
            cuts[0] = cut;
            cuts[0].type = INTERSECT_IN_A;
            cuts[1].index = cuts[0].index;
            cuts[1].type = INTERSECT_IN_C;
            cuts[1].swap = true;
            cuts[1].loc = *secondaryIterator;
            continue;
          }
          default: {
            // cases for cuts in b and d are not stored, because
            // they would give duplicate cuts (polygons are closed!)
            continue; // break without storing
          }
        }
        // here we know we have an additional cut.
        // we store max. 2 cuts as a linear edge cannot cut through more than two edges
        if (nCuts < 2)
          cuts[nCuts] = cut; // store the cut
        ++nCuts;
        // next line of passive polygon
      } while (!(++secondaryIterator).AtBegin());

      // exit loop if first active line with cut is found
      if (nCuts > 0)
        break;
      // next line of active polygon
    } while (!(++primaryIterator).AtBegin());

    // do not proceed if there is no cut to start with and
    // make sure there are not more cuts than possible
    if (nCuts == 0) {
      return false;
    }
    else if (nCuts > 2) {
      WARN("Detected more than 2 line cuts!" << std::endl <<
           "A line cannot cut more than two edges of a convex polygon." << std::endl <<
           "This can occur, e.g. if two elements touch on a node or a line(2D)." << std::endl <<
           "Ignoring this pair of elements. Please check the intersection grid.");
      return false;
    }

    // save the position of the first cut in the active polygon
    primaryIterator.SetBegin(primaryIterator.GetPos());

    if (nCuts == 2) {
      // make sure we do not treat a duplicate cut
      tempVec1 = (cuts[1].loc - cuts[0].loc);
      if (tempVec1.NormL2() < absoluteTolerance_) {
        nCuts = 1;
      } else {
        // Here we can assume that we have found two "real" cuts. In
        // this case [a,b] runs completely through p2.
        // => sort cuts by distance to a
        tempVec1 = (cuts[0].loc - *primaryIterator);
        tempVec2 = (cuts[1].loc - *primaryIterator);
        if (tempVec1.NormL2() < tempVec2.NormL2())
        {
          intersectionPolygonCoords.Push_back(cuts[0].loc);
          intersectionPolygonCoords.Push_back(cuts[1].loc);
          secondaryIterator.Seek(cuts[0].index);
          if ((cuts[0].type != INTERSECT_IN_C) && (cuts[0].type != INTERSECT_A_EQ_C)) {
            ++secondaryIterator;
          }
          secondaryIterator.SetBegin();
          secondaryIterator.Seek(cuts[1].index);
        } else {
          intersectionPolygonCoords.Push_back(cuts[1].loc);
          intersectionPolygonCoords.Push_back(cuts[0].loc);
          secondaryIterator.Seek(cuts[1].index);
          if ((cuts[1].type != INTERSECT_IN_C) && (cuts[1].type != INTERSECT_A_EQ_C)) {
            ++secondaryIterator;
          }
          secondaryIterator.SetBegin();
          secondaryIterator.Seek(cuts[0].index);
        }
        ++primaryIterator; // avoid finding the same cut twice
        primaryIterator.Swap(secondaryIterator); // continue with p2
      }
    } // NO else clause here in order to catch duplicate cuts

    if (nCuts == 1) {
      // save the position of the first cut with the passive polygon
      secondaryIterator.Seek(cuts[0].index);
      if ((cuts[0].type != INTERSECT_IN_C) && (cuts[0].type != INTERSECT_A_EQ_C))
        ++secondaryIterator;
      secondaryIterator.SetBegin();

      // store first point of intersection polygon
      intersectionPolygonCoords.Push_back(cuts[0].loc);
      // avoid finding the same cut twice
      ++primaryIterator;
      // continue with p2, if indicated
      if (cuts[0].swap) {
        primaryIterator.Swap(secondaryIterator);
      } else { // [a,b] cuts into p2, so add b
        intersectionPolygonCoords.Push_back(*primaryIterator);
      }
    }
    // indicate if we need to apply swap in the primary iterator
    bool doSwap = false;
    // indicate if the primary iterator is swapped
    bool isSwapped = false;
    // get current positions of primary and secondary polygon iterators
    UInt primaryStartPos = primaryIterator.GetPos();
    UInt secondaryStartPos = secondaryIterator.GetPos();
    // main loop
    do {
      doSwap = false;
      if (!secondaryIterator.AtBegin() || !isSwapped ) {
        do {
          switch (CutLines(*primaryIterator, primaryIterator.Next(),
                          *secondaryIterator, secondaryIterator.Next(),
                          intersectionPointCoords)) {
            case INTERSECT_CROSS:
            case INTERSECT_IN_C: {
              intersectionPolygonCoords.Push_back(intersectionPointCoords);
              doSwap = true;
              break; // break the switch
            }
            case INTERSECT_IN_A:
            case INTERSECT_A_AND_C: {
              if (CutLines(*primaryIterator, primaryIterator.Next(),
                          *secondaryIterator, secondaryIterator.Next(2),
                          intersectionPointCoords) >= INTERSECT_ON_LINE2)
                continue; // continue to next loop
            }
            case INTERSECT_A_EQ_C: {
              if (CutLines(*primaryIterator, primaryIterator.Next(),
                          secondaryIterator.Prev(), secondaryIterator.Next(),
                          intersectionPointCoords) >= INTERSECT_ON_LINE2)
                continue; // continue to next loop
              doSwap = true;
              break; // break the switch
            }
            default:
              break; // break the switch
          }
          if (doSwap)
            break;
        } while ( ! (++secondaryIterator).AtBegin() );
      }
      ++primaryIterator;
      if (doSwap) {
        primaryIterator.Swap(secondaryIterator);
        primaryStartPos = primaryIterator.GetPos();
        primaryIterator.Seek(primaryStartPos);
        secondaryStartPos = secondaryIterator.GetPos();
        isSwapped = true;
      } else {
        intersectionPolygonCoords.Push_back(*primaryIterator);
        // Return to the point directly after the last cut (we can do
        // this due to the polygons being convex and having the same
        // orientation).
        secondaryIterator.Seek(secondaryStartPos);
      }
    } while ( !primaryIterator.AtEnd() );

    intersectionPolygonCoords.Erase(intersectionPolygonCoords.GetSize() - 1);
    // we succeeded if more than two intersection points are found
    return (intersectionPolygonCoords.GetSize() > 2);
  }
#endif


  bool MortarInterface::PointInsidePoly(const Vector<Double> &p,
                                        const StdVector<Vector<Double>> &poly,
                                        const Vector<Double> *const c) {
    bool result = false;
    LineIntersectType s;
    Vector<Double> center, e, temp;
    ConstPolygonIterator pi(poly);

    // compute centroid of polygon, if not given
    if (c == nullptr)
      PolyCentroid(poly, center);
    else
      center = *c;

    // Test if p is the centroid of the polygon (should always lie inside of a
    // convex polygon). In this case the algorithm below will not work.
    temp = (p - center);
    if ( temp.NormL2() < absoluteTolerance_)
      return true;

    // try intersecting [c,p] with each edge of the polygon
    do {
      s = CutLines(center, p, *pi, pi.Next(), e);
      if (s <= INTERSECT_OUTSIDE)
        continue;
      if ((s == INTERSECT_ON_LINE2) || (s == INTERSECT_IN_B)) {
        result = true;
        break;
      }
      if ((s == INTERSECT_CROSS) || (s >= INTERSECT_IN_C)) {
        result = false;
        break;
      }
    } while ( ! (++pi).AtBegin() );

    return result;
  }

  Double MortarInterface::PolyCircumcircle(const StdVector< Vector<Double> > &p,
                          const Vector<Double> &c){
    UInt i,j, d=c.GetSize(), n = p.GetSize();
    Double r = 0.0, r_max = 0.0, tmp=0.0;

    // find point with maximum distance from centroid
    for (i = 0; i < n; ++i) {
      tmp = 0;
      for(j=0;j<d;j++){
        tmp += (p[i][j] - c[j])*(p[i][j] - c[j]);
      }
      r = sqrt(tmp);
      if (r > r_max){
        r_max = r;
      }
    }

    return r_max;
  }

  void MortarInterface::PolyCentroid(const StdVector< Vector<Double> > &p,
                            Vector<Double> &c)
  {
    UInt i, n = p.GetSize();
    // set c to 0
    c.Resize(3);
    c.Init(0.0);

    // compute center of gravity
    for (i = 0; i < n; ++i){
      c[0] += p[i][0];
      c[1] += p[i][1];
      c[2] += p[i][2];
    }
    for(i=0;i<3;i++)
      c[i] /= (Double) n;
  }

  void MortarInterface::AddNodeToGrid(const Vector<Double>& nodeCoordinate, UInt& newNodeId, StdVector<UInt>& newNodeIds) {
  #ifdef USE_OPENMP
    omp_set_lock(&gridLock_);
  #endif
    ptGrid_->AddNode(nodeCoordinate, newNodeId);
    newNodeIds.Push_back(newNodeId);
  #ifdef USE_OPENMP
    omp_unset_lock(&gridLock_);
  #endif
  }


  void MortarInterface::TriangulatePoly(const StdVector<Vector<Double>>& inputPolygon,
                                        StdVector<MortarNcSurfElem*>& outputTriangles,
                                        StdVector<UInt>& newNodeIds) {
    // set variables...
    // the id of the new node on the grid
    UInt newNodeId;
    // number of node coords passed
    UInt numNodes = inputPolygon.GetSize();
    // pointer to the new elements
    MortarNcSurfElem* ncElem;
    if (numNodes > 4) {
      UInt centerNodeId, firstNodeId;
      // get the centroid of the intersection polygon and add it to grid
      Vector<Double> centerNodeCoords;
      PolyCentroid(inputPolygon, centerNodeCoords);
      AddNodeToGrid(centerNodeCoords, centerNodeId, newNodeIds);
      // add the first node of the intersection polygon to the grid
      // AddNodeToGrid() returns the newNodeId as id of the newly assigned node.
      // so newNodeId gets incremented by each call
      AddNodeToGrid(inputPolygon[0], newNodeId, newNodeIds);
      firstNodeId = newNodeId;
      // triangulate the polygon using the additional center node
      for (UInt iNodes = 1; iNodes < numNodes; ++iNodes) {
        ncElem = new MortarNcSurfElem;
        ncElem->type = Elem::ET_TRIA3;
        ncElem->connect.Resize(3);
        ncElem->connect[0] = newNodeId;
        AddNodeToGrid(inputPolygon[iNodes], newNodeId, newNodeIds);
        ncElem->connect[1] = newNodeId;
        ncElem->connect[2] = centerNodeId;
        outputTriangles.Push_back(ncElem);
      }
      // assign the last trianlge connecting the last and first nodes
      ncElem = new MortarNcSurfElem;
      ncElem->type = Elem::ET_TRIA3;
      ncElem->connect.Resize(3);
      ncElem->connect[0] = newNodeId;
      ncElem->connect[1] = firstNodeId;
      ncElem->connect[2] = centerNodeId;
      outputTriangles.Push_back(ncElem);
    }
    // if we have only 3 or 4 intersection nodes, we can directly assign a linear tri or quad element
    else if (numNodes == 4) {
      ncElem = new MortarNcSurfElem;
      ncElem->type = Elem::ET_QUAD4;
      ncElem->connect.Resize(numNodes);
      // add the nodes to the grid
      for (UInt iNodes = 0; iNodes < numNodes; ++iNodes) {
        AddNodeToGrid(inputPolygon[iNodes], newNodeId, newNodeIds);
        ncElem->connect[iNodes] = newNodeId;
      }
      outputTriangles.Push_back(ncElem);
    }
    else if (numNodes == 3) {
      // if we have only 3 or 4 intersection nodes, we can directly assign a linear tri or quad element
      ncElem = new MortarNcSurfElem;
      ncElem->type = Elem::ET_TRIA3;
      ncElem->connect.Resize(numNodes);
      // add the nodes to the grid
      for (UInt iNodes = 0; iNodes < numNodes; ++iNodes) {
        AddNodeToGrid(inputPolygon[iNodes], newNodeId, newNodeIds);
        ncElem->connect[iNodes] = newNodeId;
      }
      outputTriangles.Push_back(ncElem);
    } else {
      EXCEPTION("Cannot triangulate less then 3 nodes.")
    }
  }

} /* namespace CoupledField */

