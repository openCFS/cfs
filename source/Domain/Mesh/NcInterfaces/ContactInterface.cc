// ================================================================================================
/*!
 *       \file     ContactInterface.cc
 *       \brief    Geometry layer for genuine two-body mechanical contact.
 *                 See ContactInterface.hh for the frozen sign/naming conventions.
 *       \date     2026
 */
//================================================================================================

#include "ContactInterface.hh"

#include "Domain/ElemMapping/Elem.hh"
#include "Domain/ElemMapping/SurfElem.hh"
#include "Domain/ElemMapping/ElemShapeMap.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "MatVec/Matrix.hh"
#include "FeBasis/FeFunctions.hh"

DEFINE_LOG(contactinterface, "contactinterface")

namespace CoupledField {

namespace {
  //! Slack applied when clamping the projection iterate. Generous on purpose: it only
  //! prevents divergence, it must not decide whether the point belongs to the element.
  const Double PROJECTION_CLAMP_SLACK = 0.5;
  //! Absolute slack on the reference domain when deciding whether a converged projection
  //! lies on the element. Sized so a point on a shared edge is claimed by a neighbour.
  const Double PROJECTION_INSIDE_TOL = 1e-8;

  //! Slack used when projecting a PRIMARY corner node onto a secondary element during
  //! segmentation. Much larger than PROJECTION_CLAMP_SLACK because here the answer is wanted
  //! precisely when it lies outside the element: a primary element several times larger than
  //! the secondary one has corners that project far outside, and truncating them would move
  //! the cut. A corner that still hits this bound is reported as `clamped` and the caller
  //! refuses to segment that element rather than cut it in the wrong place.
  const Double SEGMENTATION_CLAMP_SLACK = 10.0;

  const Double SEGMENTATION_MERGE_TOL = 1e-8;

  const Double SEGMENTATION_AREA_TOL = 1e-8;
}

// ================================================================================================
//  Construction
// ================================================================================================

ContactInterface::ContactInterface(Grid* grid, PtrParamNode ctNode)
  : BaseNcInterface(grid),
    primarySurfRegion_(-1), secondarySurfRegion_(-1),
    primaryVolRegion_(-1), secondaryVolRegion_(-1),
    integrationOrder_(2),
    slidingType_(SLIDING_SMALL),
    normalSmoothing_(NORMAL_FACET),
    formulation_(FORM_PENALTY),
    augmentedMaxIter_(20),
    augmentedTol_(1e-8),
    mortarSegmentation_(true),
    numUnsegmented_(0),
    normalPenalty_(100.0),
    frictionType_(FRICTION_NONE),
    frictionMu_(0.0),
    frictionTangent_(FRICTION_TANGENT_REDUCED),
    tangentialPenalty_(100.0),
    stabilization_(0.0),
    refModulus_(0.0),
    applyForce_(true),
    zeroInitialGap_(false),
    activationTolerance_(1e-12),
    projectionTolerance_(1e-10),
    projectionMaxIter_(20),
    boundingBoxTolerance_(0.1),
    searchDistanceFactor_(2.0),
    segmentsBuilt_(false),
    projected_(false)
{
  SetName(ctNode->Get("name")->As<std::string>());
  SetRegions(ctNode);
  SetOptions(ctNode);
  CollectSurfaceElements();
}

ContactInterface::~ContactInterface() {
  segments_.Clear();
  primarySurfElems_.Clear();
  secondarySurfElems_.Clear();
}

// ================================================================================================
//  XML input
// ================================================================================================

void ContactInterface::SetRegions(const PtrParamNode ctNode) {

  const std::string primaryName   = ctNode->Get("primarySurface")->As<std::string>();
  const std::string secondaryName = ctNode->Get("secondarySurface")->As<std::string>();

  if (primaryName == secondaryName) {
    EXCEPTION("Contact pair '" << GetName() << "': primary and secondary surface are the "
              << "same region '" << primaryName << "'. Self-contact is not supported.");
  }

  primarySurfRegion_   = ptGrid_->GetRegionId(primaryName);
  secondarySurfRegion_ = ptGrid_->GetRegionId(secondaryName);

  // both must be surfaces, i.e. of dimension dim-1
  const UInt surfDim = ptGrid_->GetDim() - 1;
  if (ptGrid_->GetEntityDim(primaryName) != surfDim) {
    EXCEPTION("Contact pair '" << GetName() << "': primary surface '" << primaryName
              << "' has dimension " << ptGrid_->GetEntityDim(primaryName)
              << " but " << surfDim << " was expected.");
  }
  if (ptGrid_->GetEntityDim(secondaryName) != surfDim) {
    EXCEPTION("Contact pair '" << GetName() << "': secondary surface '" << secondaryName
              << "' has dimension " << ptGrid_->GetEntityDim(secondaryName)
              << " but " << surfDim << " was expected.");
  }
}

void ContactInterface::SetOptions(const PtrParamNode ctNode) {

  std::string slidingStr = "small";
  ctNode->GetValue("slidingType", slidingStr, ParamNode::PASS);
  if (slidingStr == "small") {
    slidingType_ = SLIDING_SMALL;
  } else if (slidingStr == "large") {
    slidingType_ = SLIDING_LARGE;
  } else {
    EXCEPTION("Contact pair '" << GetName() << "': unknown slidingType '" << slidingStr
              << "'. Use 'small' or 'large'.");
  }

  // Faceted normals jump at element boundaries, which makes g_N non-differentiable exactly
  // where a point crosses from one primary element to the next. "nodal" removes the jump; it
  // is not the default because on a flat interface it is provably identical and on a faceted
  // one it changes the discrete normal, hence every existing answer.
  std::string smoothStr = "none";
  ctNode->GetValue("normalSmoothing", smoothStr, ParamNode::PASS);
  if (smoothStr == "none") {
    normalSmoothing_ = NORMAL_FACET;
  } else if (smoothStr == "nodal") {
    normalSmoothing_ = NORMAL_NODAL;
  } else {
    EXCEPTION("Contact pair '" << GetName() << "': unknown normalSmoothing '" << smoothStr
              << "'. Use 'none' or 'nodal'.");
  }

  ctNode->GetValue("normalPenalty", normalPenalty_, ParamNode::PASS);
  if (normalPenalty_ <= 0.0) {
    EXCEPTION("Contact pair '" << GetName() << "': normalPenalty must be positive, got "
              << normalPenalty_);
  }

  ctNode->GetValue("stabilization", stabilization_, ParamNode::PASS);
  if (stabilization_ < 0.0) {
    EXCEPTION("Contact pair '" << GetName() << "': stabilization must not be negative, got "
              << stabilization_);
  }

  std::string applyStr = "yes";
  ctNode->GetValue("applyContactForce", applyStr, ParamNode::PASS);
  applyForce_ = (applyStr != "no");

  std::string gapStr = "fromMesh";
  ctNode->GetValue("initialGapHandling", gapStr, ParamNode::PASS);
  zeroInitialGap_ = (gapStr == "zero");

  std::string formStr = "penalty";
  ctNode->GetValue("formulation", formStr, ParamNode::PASS);
  if (formStr == "penalty") {
    formulation_ = FORM_PENALTY;
  } else if (formStr == "augmentedLagrange") {
    formulation_ = FORM_AUGMENTED_LAGRANGE;
  } else {
    EXCEPTION("Contact pair '" << GetName() << "': unknown formulation '" << formStr
              << "'. Use 'penalty' or 'augmentedLagrange'.");
  }

  PtrParamNode fricNode = ctNode->Get("friction", ParamNode::PASS);
  if (fricNode) {

    std::string fricStr = "coulomb";
    fricNode->GetValue("type", fricStr, ParamNode::PASS);
    if (fricStr == "coulomb") {
      frictionType_ = FRICTION_COULOMB;
    } else if (fricStr == "none") {
      frictionType_ = FRICTION_NONE;
    } else {
      EXCEPTION("Contact pair '" << GetName() << "': unknown friction type '" << fricStr
                << "'. Use 'coulomb' or 'none'.");
    }

    fricNode->GetValue("mu", frictionMu_, ParamNode::PASS);
    if (frictionMu_ < 0.0) {
      EXCEPTION("Contact pair '" << GetName() << "': friction mu must not be negative, got "
                << frictionMu_);
    }

    std::string tangStr = "reduced";
    fricNode->GetValue("tangent", tangStr, ParamNode::PASS);
    if (tangStr == "reduced") {
      frictionTangent_ = FRICTION_TANGENT_REDUCED;
    } else if (tangStr == "consistent") {
      frictionTangent_ = FRICTION_TANGENT_CONSISTENT;
    } else {
      EXCEPTION("Contact pair '" << GetName() << "': unknown friction tangent '" << tangStr
                << "'. Use 'reduced' or 'consistent'.");
    }

    fricNode->GetValue("tangentialPenalty", tangentialPenalty_, ParamNode::PASS);
    if (tangentialPenalty_ <= 0.0) {
      EXCEPTION("Contact pair '" << GetName() << "': friction tangentialPenalty must be "
                << "positive, got " << tangentialPenalty_);
    }

  }

  PtrParamNode augNode = ctNode->Get("augmentation", ParamNode::PASS);
  if (augNode) {
    augNode->GetValue("maxIter", augmentedMaxIter_, ParamNode::PASS);
    augNode->GetValue("tol", augmentedTol_, ParamNode::PASS);
    if (augmentedMaxIter_ < 1) {
      EXCEPTION("Contact pair '" << GetName() << "': augmentation maxIter must be at least 1.");
    }
    if (augmentedTol_ <= 0.0) {
      EXCEPTION("Contact pair '" << GetName() << "': augmentation tol must be positive.");
    }
  }

  std::string segStr = "mortar";
  ctNode->GetValue("segmentation", segStr, ParamNode::PASS);
  if (segStr == "mortar") {
    mortarSegmentation_ = true;
  } else if (segStr == "gpts") {
    mortarSegmentation_ = false;
  } else {
    EXCEPTION("Contact pair '" << GetName() << "': unknown segmentation '" << segStr
              << "'. Use 'mortar' or 'gpts'.");
  }

  // --------------------------------------------------------------------------------------
  //  Options that carry PER-POINT STATE hold the partition fixed inside a solve
  // --------------------------------------------------------------------------------------
  // Under SLIDING_LARGE with segmentation="mortar" the partition is re-cut on every
  // re-projection (RebuildSegmentPoints()), because its sub-cell boundaries are the projected
  // primary element outlines and those move with the body. That replaces the integration
  // points: a point is defined by where it sits, so a new one has no past.
  //
  partitionFrozenInSolve_ = (slidingType_ == SLIDING_LARGE) && mortarSegmentation_
                            && NeedsStablePointSet();

  ctNode->GetValue("integrationOrder", integrationOrder_, ParamNode::PASS);
  ctNode->GetValue("projectionTolerance", projectionTolerance_, ParamNode::PASS);
  ctNode->GetValue("boundingBoxTolerance", boundingBoxTolerance_, ParamNode::PASS);
  ctNode->GetValue("searchDistanceFactor", searchDistanceFactor_, ParamNode::PASS);
}

// ================================================================================================
//  Surface element collection
// ================================================================================================

void ContactInterface::CollectSurfaceElements() {

  ptGrid_->GetSurfElems(primarySurfElems_, primarySurfRegion_);
  ptGrid_->GetSurfElems(secondarySurfElems_, secondarySurfRegion_);

  if (primarySurfElems_.GetSize() == 0) {
    EXCEPTION("Contact pair '" << GetName() << "': primary surface has no elements.");
  }
  if (secondarySurfElems_.GetSize() == 0) {
    EXCEPTION("Contact pair '" << GetName() << "': secondary surface has no elements.");
  }

  // Determine the adjacent volume regions and, at the same time, verify the normal
  // orientation assumption: ElemShapeMap::CalcNormal returns the normal
  // pointing out of ptVolElems[0], so every primary surface element must have the primary
  // body as its FIRST volume neighbour. If a surface sits between two meshed bodies this
  // may not hold, and the resulting normals would be inconsistent across the surface.
  primaryVolRegion_ = -1;
  for (UInt i = 0; i < primarySurfElems_.GetSize(); ++i) {
    Elem* vol = primarySurfElems_[i]->ptVolElems[0];
    if (vol == nullptr) {
      EXCEPTION("Contact pair '" << GetName() << "': primary surface element "
                << primarySurfElems_[i]->elemNum << " has no volume neighbour.");
    }
    if (primaryVolRegion_ == (RegionIdType) -1) {
      primaryVolRegion_ = vol->regionId;
    } else if (primaryVolRegion_ != vol->regionId) {
      EXCEPTION("Contact pair '" << GetName() << "': the primary surface borders more than "
                << "one volume region. The surface normal orientation would be ambiguous. "
                << "Split the contact pair into one pair per volume region.");
    }
  }

  secondaryVolRegion_ = -1;
  for (UInt i = 0; i < secondarySurfElems_.GetSize(); ++i) {
    Elem* vol = secondarySurfElems_[i]->ptVolElems[0];
    if (vol == nullptr) {
      EXCEPTION("Contact pair '" << GetName() << "': secondary surface element "
                << secondarySurfElems_[i]->elemNum << " has no volume neighbour.");
    }
    if (secondaryVolRegion_ == (RegionIdType) -1) {
      secondaryVolRegion_ = vol->regionId;
    } else if (secondaryVolRegion_ != vol->regionId) {
      EXCEPTION("Contact pair '" << GetName() << "': the secondary surface borders more than "
                << "one volume region.");
    }
  }

  if (primaryVolRegion_ == secondaryVolRegion_) {
    EXCEPTION("Contact pair '" << GetName() << "': both surfaces belong to the same volume "
              << "region. Two-body contact needs two distinct bodies.");
  }

  primaryVolRegionSet_.clear();
  primaryVolRegionSet_.insert(primaryVolRegion_);
  secondaryVolRegionSet_.clear();
  secondaryVolRegionSet_.insert(secondaryVolRegion_);

  LOG_DBG(contactinterface) << "contact pair '" << GetName() << "': "
      << primarySurfElems_.GetSize() << " primary / "
      << secondarySurfElems_.GetSize() << " secondary surface elements";
}

// ================================================================================================
//  Lifecycle
// ================================================================================================

void ContactInterface::ResetInterface() {
  ClearGeometryCaches();
  segments_.Clear();
  graphPrimaries_.Clear();
  elemPairs_.clear();
  elemPairsFixed_ = false;
  segmentsBuilt_ = false;
  projected_ = false;
}

void ContactInterface::UpdateInterface() {
  if (!segmentsBuilt_) {
    BuildSegments();
  }
  UpdateProjection();
}

// ================================================================================================
//  Mortar segmentation
// ================================================================================================

void ContactInterface::GetCornerCoords(SurfElem* elem, RegionIdType volRegion,
                                       StdVector<Vector<Double> >& coords) const {

  const UInt nv = Elem::shapes[elem->type].numVertices;

  if (slidingType_ == SLIDING_LARGE && dispFct_) {
    Matrix<Double> ec;
    GetDeformedElemCoords(elem, volRegion, ec);

    const UInt sdim = ptGrid_->GetDim();
    coords.Resize(nv);
    for (UInt i = 0; i < nv; ++i) {
      coords[i].Resize(sdim);
      for (UInt d = 0; d < sdim; ++d) {
        coords[i][d] = ec[d][i];
      }
    }
    return;
  }

  StdVector<UInt> nodes;
  ptGrid_->GetElemNodes(nodes, elem->elemNum);

  StdVector<UInt> corners(nv);
  for (UInt i = 0; i < nv; ++i) {
    corners[i] = nodes[i];
  }
  ptGrid_->GetNodeCoordinates(coords, corners, false);
}

void ContactInterface::GetReferencePolygon(Elem::ShapeType shape,
                                           StdVector<Vector<Double> >& poly) {

  poly.Clear();
  auto add = [&poly](Double a, Double b) {
    Vector<Double> v(2);
    v[0] = a; v[1] = b;
    poly.Push_back(v);
  };

  switch (shape) {
    case Elem::ST_QUAD:                       // [-1,1] x [-1,1], counter-clockwise
      add(-1.0, -1.0); add(1.0, -1.0); add(1.0, 1.0); add(-1.0, 1.0);
      break;
    case Elem::ST_TRIA:                       // { xi >= 0, eta >= 0, xi + eta <= 1 }
      add(0.0, 0.0); add(1.0, 0.0); add(0.0, 1.0);
      break;
    default:
      EXCEPTION("ContactInterface: no reference polygon for shape '"
                << Elem::shapeType.ToString(shape) << "'.");
  }
}

Double ContactInterface::PolygonArea(const StdVector<Vector<Double> >& poly) {

  const UInt n = poly.GetSize();
  if (n < 3) {
    return 0.0;
  }
  Double a = 0.0;
  for (UInt i = 0; i < n; ++i) {
    const Vector<Double>& p = poly[i];
    const Vector<Double>& q = poly[(i + 1) % n];
    a += p[0] * q[1] - q[0] * p[1];
  }
  return 0.5 * a;
}

void ContactInterface::ClipPolygon(const StdVector<Vector<Double> >& subject,
                                   const StdVector<Vector<Double> >& clip,
                                   StdVector<Vector<Double> >& result) {

  // Sutherland-Hodgman: clip the subject successively against the half-plane left of each
  // clip edge. Correct for a CONVEX clip polygon, which the caller guarantees; the subject
  // (a reference square or triangle) is convex too, so the result is a convex polygon.
  result = subject;

  const UInt nc = clip.GetSize();
  StdVector<Vector<Double> > input;

  for (UInt e = 0; e < nc && result.GetSize() > 0; ++e) {

    input = result;
    result.Clear();

    const Vector<Double>& a = clip[e];
    const Vector<Double>& b = clip[(e + 1) % nc];
    const Double ex = b[0] - a[0];
    const Double ey = b[1] - a[1];

    // signed distance to the edge line, positive on the inside (left) for a CCW clip polygon
    auto side = [&](const Vector<Double>& p) {
      return ex * (p[1] - a[1]) - ey * (p[0] - a[0]);
    };

    const UInt ni = input.GetSize();
    for (UInt i = 0; i < ni; ++i) {

      const Vector<Double>& cur  = input[i];
      const Vector<Double>& prev = input[(i + ni - 1) % ni];
      const Double dCur  = side(cur);
      const Double dPrev = side(prev);

      if (dCur >= 0.0) {
        if (dPrev < 0.0) {
          Vector<Double> x(2);
          const Double t = dPrev / (dPrev - dCur);
          x[0] = prev[0] + t * (cur[0] - prev[0]);
          x[1] = prev[1] + t * (cur[1] - prev[1]);
          result.Push_back(x);
        }
        result.Push_back(cur);
      } else if (dPrev >= 0.0) {
        Vector<Double> x(2);
        const Double t = dPrev / (dPrev - dCur);
        x[0] = prev[0] + t * (cur[0] - prev[0]);
        x[1] = prev[1] + t * (cur[1] - prev[1]);
        result.Push_back(x);
      }
    }
  }
}

void ContactInterface::FindElementCandidates(SurfElem* secElem,
                                             StdVector<UInt>& candidates) const {

  const UInt dim = ptGrid_->GetDim();
  candidates.Clear();

  // Same configuration as the primary boxes this is tested against, and as the projection
  // that follows: reference under small sliding, deformed under large. Boxes drawn in two
  // different configurations would pre-select candidates by a distance that is neither.
  Matrix<Double> ec;
  if (slidingType_ == SLIDING_LARGE && dispFct_) {
    GetDeformedElemCoords(secElem, secondaryVolRegion_, ec);
  } else {
    ptGrid_->GetElemNodesCoord(ec, secElem->connect, false);
  }
  const UInt numNodes = ec.GetNumCols();

  Vector<Double> lo(dim), hi(dim);
  for (UInt d = 0; d < dim; ++d) {
    lo[d] = ec[d][0];
    hi[d] = ec[d][0];
  }
  for (UInt n = 1; n < numNodes; ++n) {
    for (UInt d = 0; d < dim; ++d) {
      lo[d] = std::min(lo[d], ec[d][n]);
      hi[d] = std::max(hi[d], ec[d][n]);
    }
  }
  for (UInt i = 0; i < primaryBoxMin_.GetSize(); ++i) {
    bool overlap = true;
    for (UInt d = 0; d < dim; ++d) {
      if (hi[d] < primaryBoxMin_[i][d] || lo[d] > primaryBoxMax_[i][d]) {
        overlap = false;
        break;
      }
    }
    if (overlap) {
      candidates.Push_back(i);
    }
  }
}

void ContactInterface::BuildSegmentCells(SurfElem* secElem, StdVector<ContactCell>& cells) {

  cells.Clear();

  const UInt dim = ptGrid_->GetDim();
  const Elem::ShapeType secShape = Elem::GetShapeType(secElem->type);

  StdVector<UInt> candidates;
  FindElementCandidates(secElem, candidates);
  if (candidates.GetSize() == 0) {
    return;                                   // no partner at all -- integrate as a whole
  }

  StdVector<Vector<Double> > corners;
  LocPoint lp;
  bool clamped = false;

  if (dim == 2) {

    // ------------------------------------------------------------------------------------
    // 2D: the cut is a set of break points on xi in [-1,1], and tiling the interval by them
    // is exact and complete by construction -- there is no coverage or overlap question, so
    // this branch can never fall back.
    // ------------------------------------------------------------------------------------
    StdVector<Double> breaks;
    breaks.Push_back(-1.0);
    breaks.Push_back(1.0);

    for (UInt iC = 0; iC < candidates.GetSize(); ++iC) {

      GetCornerCoords(primarySurfElems_[candidates[iC]], primaryVolRegion_, corners);

      for (UInt k = 0; k < corners.GetSize(); ++k) {
        if (!ClosestPointLocal(corners[k], secElem, secondaryVolRegion_,
                               SEGMENTATION_CLAMP_SLACK, lp, clamped)) {
          continue;
        }
        // A clamped coordinate is a guard value, and one that lies outside the element does
        // not cut it. Both are simply not break points.
        if (clamped) {
          continue;
        }
        const Double xi = lp.coord[0];
        if (xi > -1.0 + SEGMENTATION_MERGE_TOL && xi < 1.0 - SEGMENTATION_MERGE_TOL) {
          breaks.Push_back(xi);
        }
      }
    }

    if (breaks.GetSize() == 2) {
      return;
    }

    std::sort(breaks.Begin(), breaks.End());

    for (UInt i = 0; i + 1 < breaks.GetSize(); ++i) {
      const Double a = breaks[i];
      const Double b = breaks[i + 1];
      if (b - a < SEGMENTATION_MERGE_TOL) {
        continue;                             // duplicate cut
      }
      ContactCell cell;
      cell.vtx.Resize(2);
      cell.vtx[0].Resize(1);
      cell.vtx[1].Resize(1);
      cell.vtx[0][0] = a;
      cell.vtx[1][0] = b;
      cells.Push_back(cell);
    }
    return;
  }

  // --------------------------------------------------------------------------------------
  // 3D: clip the secondary reference polygon against each projected primary outline and
  // fan-triangulate what survives.
  // --------------------------------------------------------------------------------------
  StdVector<Vector<Double> > refPoly;
  GetReferencePolygon(secShape, refPoly);
  const Double refArea = PolygonArea(refPoly);

  StdVector<Vector<Double> > primPoly, clipped;
  StdVector<StdVector<Vector<Double> > > pieces;
  Double coveredArea = 0.0;

  for (UInt iC = 0; iC < candidates.GetSize(); ++iC) {

    GetCornerCoords(primarySurfElems_[candidates[iC]], primaryVolRegion_, corners);

    primPoly.Clear();
    bool usable = true;

    for (UInt k = 0; k < corners.GetSize(); ++k) {
      if (!ClosestPointLocal(corners[k], secElem, secondaryVolRegion_,
                               SEGMENTATION_CLAMP_SLACK, lp, clamped)
          || clamped) {
        usable = false;
        break;
      }
      Vector<Double> v(2);
      v[0] = lp.coord[0];
      v[1] = lp.coord[1];
      primPoly.Push_back(v);
    }

    if (!usable || primPoly.GetSize() < 3) {
      continue;
    }

    // Sutherland-Hodgman needs a counter-clockwise, convex clip polygon. The two surfaces
    // face each other, so a primary outline generally appears mirrored in the secondary
    // reference domain and has to be reversed.
    if (PolygonArea(primPoly) < 0.0) {
      StdVector<Vector<Double> > rev;
      for (UInt k = primPoly.GetSize(); k > 0; --k) {
        rev.Push_back(primPoly[k - 1]);
      }
      primPoly = rev;
    }

    bool convex = true;
    const UInt np = primPoly.GetSize();
    for (UInt k = 0; k < np; ++k) {
      const Vector<Double>& p0 = primPoly[k];
      const Vector<Double>& p1 = primPoly[(k + 1) % np];
      const Vector<Double>& p2 = primPoly[(k + 2) % np];
      const Double cross = (p1[0] - p0[0]) * (p2[1] - p1[1])
                         - (p1[1] - p0[1]) * (p2[0] - p1[0]);
      if (cross < -SEGMENTATION_MERGE_TOL) {
        convex = false;
        break;
      }
    }
    if (!convex) {
      continue;
    }

    ClipPolygon(refPoly, primPoly, clipped);
    const Double area = PolygonArea(clipped);
    if (clipped.GetSize() < 3 || area < SEGMENTATION_MERGE_TOL * std::abs(refArea)) {
      continue;
    }
    coveredArea += area;
    pieces.Push_back(clipped);
  }

  // Sub-cells covering MORE than the element would integrate part of the interface twice and
  // transmit too much load. That can only come from projected primary outlines overlapping,
  // i.e. from a projection this code cannot trust, so refuse the segmentation outright rather
  // than assemble a silently wrong force -- the failure mode this whole phase exists to avoid.
  if (coveredArea > std::abs(refArea) * (1.0 + SEGMENTATION_AREA_TOL)) {
    ++numUnsegmented_;
    return;
  }

  if (pieces.GetSize() == 1
      && std::abs(PolygonArea(pieces[0]) - refArea)
           <= SEGMENTATION_AREA_TOL * std::abs(refArea)) {
    return;
  }

  // fan-triangulate each convex clip polygon about its first vertex
  for (UInt p = 0; p < pieces.GetSize(); ++p) {
    const StdVector<Vector<Double> >& poly = pieces[p];
    for (UInt k = 1; k + 1 < poly.GetSize(); ++k) {
      ContactCell cell;
      cell.vtx.Resize(3);
      cell.vtx[0] = poly[0];
      cell.vtx[1] = poly[k];
      cell.vtx[2] = poly[k + 1];
      cells.Push_back(cell);
    }
  }
}

// ================================================================================================
//  Segment generation
// ================================================================================================

void ContactInterface::GenerateSegmentPoints(ContactSegment& seg) {

  const UInt dim = ptGrid_->GetDim();
  const IntegOrder order(integrationOrder_);

  StdVector<LocPoint> intPoints;
  StdVector<Double> intWeights;
  Matrix<Double> jac;
  StdVector<ContactCell> cells;

  const Elem::ShapeType cellShape = (dim == 3) ? Elem::ST_TRIA : Elem::ST_LINE;

  SurfElem* elem = seg.secElem;
  const Elem::ShapeType shape = Elem::GetShapeType(elem->type);

  // Weights stay on the REFERENCE configuration, consistent with the total-Lagrangian
  // formulation MechPDE uses. This is deliberately NOT the configuration the
  // partition is cut in: where the cut falls is a pairing question and follows the
  // deformation, while what a weight measures is an integration question and does not.
  shared_ptr<ElemShapeMap> esm = ptGrid_->GetElemShapeMap(elem, false);

  intScheme_.GetIntPoints(shape, IntScheme::GAUSS, order, intPoints, intWeights);

  // ------------------------------------------------------------------------------------
  // Integration points: one Gauss set per mortar sub-cell, or the element's own rule when
  // the element is not cut.
  // ------------------------------------------------------------------------------------
  cells.Clear();
  if (mortarSegmentation_) {
    BuildSegmentCells(elem, cells);
  }

  seg.points.Clear();

  if (cells.GetSize() == 0) {

      seg.numCells = 1;

      for (UInt iP = 0; iP < intPoints.GetSize(); ++iP) {
        ContactPoint cp;
        cp.secElem  = elem;
        cp.secLocal = intPoints[iP];
        esm->CalcJ(jac, cp.secLocal);
        cp.weight = intWeights[iP] * esm->CalcJDet(jac, cp.secLocal);
        seg.points.Push_back(cp);
      }

  } else {

      seg.numCells = cells.GetSize();

      // --------------------------------------------------------------------------------
      // The order of the SUB-CELL rule is not the user's integrationOrder.
      //
      // A sub-cell rule must integrate the secondary element's surface shape functions
      // exactly. That is what reproduces a constant contact pressure -- the contact patch
      // test -- and it is a property of the shape functions, not of how accurately the user
      // wants the traction resolved.
      //
      // It matters because the 3D sub-cells are TRIANGLES while the element's own rule is a
      // tensor-product one, and the two are not comparable at equal nominal order: on a
      // triangle IntegOrder(n) is exact for total degree n, so order 2 gives 3 points and
      // total degree 2, while the 2x2 tensor rule it replaces is exact for degree 3 in each
      // direction. A QUAD9 face carries biquadratic shape functions, total degree
      // (dim-1)*2 = 4, so a 3-point rule integrates them wrongly
      //
      const UInt shapeFncDegree = (dim - 1) * Elem::shapes[elem->type].order;
      const IntegOrder cellOrder(std::max(integrationOrder_, shapeFncDegree));

      StdVector<LocPoint> cellPoints;
      StdVector<Double> cellWeights;
      intScheme_.GetIntPoints(cellShape, IntScheme::GAUSS, cellOrder, cellPoints, cellWeights);

      Double wsum = 0.0;
      for (UInt iP = 0; iP < cellWeights.GetSize(); ++iP) {
        wsum += cellWeights[iP];
      }
      if (wsum <= 0.0) {
        EXCEPTION("Contact pair '" << GetName() << "': the integration rule for the mortar "
                  << "sub-cells has non-positive total weight.");
      }

      for (UInt iC = 0; iC < cells.GetSize(); ++iC) {

        const ContactCell& cell = cells[iC];

        // Measure of the cell in the secondary reference domain: the interval length in 2D,
        // the triangle area in 3D. The PHYSICAL weight is this times the pointwise surface
        // Jacobian, so a curved element is still integrated exactly -- only the partition
        // boundary is linearized.
        Double measure = 0.0;
        if (dim == 2) {
          measure = cell.vtx[1][0] - cell.vtx[0][0];
        } else {
          measure = 0.5 * std::abs((cell.vtx[1][0] - cell.vtx[0][0])
                                     * (cell.vtx[2][1] - cell.vtx[0][1])
                                 - (cell.vtx[2][0] - cell.vtx[0][0])
                                     * (cell.vtx[1][1] - cell.vtx[0][1]));
        }
        if (measure <= 0.0) {
          continue;
        }

        for (UInt iP = 0; iP < cellPoints.GetSize(); ++iP) {

          ContactPoint cp;
          cp.secElem = elem;
          cp.secLocal.number = LocPoint::NOT_SET;
          cp.secLocal.coord.Resize(dim - 1);

          if (dim == 2) {
            // map [-1,1] onto [a,b]
            const Double t = 0.5 * (cellPoints[iP].coord[0] + 1.0);
            cp.secLocal.coord[0] = cell.vtx[0][0] + t * measure;
          } else {
            // affine map of the reference triangle onto the cell
            const Double r = cellPoints[iP].coord[0];
            const Double s = cellPoints[iP].coord[1];
            for (UInt d = 0; d < 2; ++d) {
              cp.secLocal.coord[d] = cell.vtx[0][d]
                                   + r * (cell.vtx[1][d] - cell.vtx[0][d])
                                   + s * (cell.vtx[2][d] - cell.vtx[0][d]);
            }
          }

          esm->CalcJ(jac, cp.secLocal);
          cp.weight = (cellWeights[iP] / wsum) * measure
                    * esm->CalcJDet(jac, cp.secLocal);

          seg.points.Push_back(cp);
        }
      }
  }

  for (UInt iP = 0; iP < seg.points.GetSize(); ++iP) {
    ContactPoint& cp = seg.points[iP];
    cp.secGlobal.Resize(dim);
    esm->Local2Global(cp.secGlobal, cp.secLocal);
    cp.normal.Resize(dim);
    cp.normal.Init(0.0);
    cp.primGlobal.Resize(dim);
    cp.primGlobal.Init(0.0);
    cp.slipPlastic.Resize(dim);
    cp.slipPlastic.Init(0.0);
    cp.slipOffset.Resize(dim);
    cp.slipOffset.Init(0.0);
    cp.tracTang.Resize(dim);
    cp.tracTang.Init(0.0);
  }
}


void ContactInterface::BuildSegments() {

  if (segmentsBuilt_) {
    return;
  }

  const UInt dim = ptGrid_->GetDim();
  const IntegOrder order(integrationOrder_);

  segments_.Resize(secondarySurfElems_.GetSize());

  StdVector<LocPoint> intPoints;
  StdVector<Double> intWeights;
  Matrix<Double> jac;

  numUnsegmented_ = 0;
  if (mortarSegmentation_) {
    UpdatePrimaryBoundingBoxes();
  }

  UInt numPoints = 0;

  for (UInt iEl = 0; iEl < secondarySurfElems_.GetSize(); ++iEl) {

    SurfElem* elem = secondarySurfElems_[iEl];
    ContactSegment& seg = segments_[iEl];
    seg.secElem = elem;

    const Elem::ShapeType shape = Elem::GetShapeType(elem->type);
    shared_ptr<ElemShapeMap> esm = ptGrid_->GetElemShapeMap(elem, false);

    intScheme_.GetIntPoints(shape, IntScheme::GAUSS, order, intPoints, intWeights);

    seg.elemArea = 0.0;
    for (UInt iP = 0; iP < intPoints.GetSize(); ++iP) {
      esm->CalcJ(jac, intPoints[iP]);
      seg.elemArea += intWeights[iP] * esm->CalcJDet(jac, intPoints[iP]);
    }
    seg.elemSize = (dim == 3) ? std::sqrt(seg.elemArea) : seg.elemArea;

    GenerateSegmentPoints(seg);
    numPoints += seg.points.GetSize();
  }

  ReserveGraphPrimaries();

  segmentsBuilt_ = true;

  LOG_DBG(contactinterface) << "contact pair '" << GetName() << "': built "
      << segments_.GetSize() << " segments (" << GetNumCells() << " mortar cells, "
      << numUnsegmented_ << " not segmented) with " << numPoints << " integration points, "
      << "integrated area = " << GetIntegratedArea()
      << ", secondary area = " << GetSecondaryArea();
}


void ContactInterface::ReserveGraphPrimaries() {

  graphPrimaries_.Clear();

  if (slidingType_ != SLIDING_LARGE || !applyForce_) {
    return;
  }

  UpdatePrimaryBoundingBoxes();

  graphPrimaries_.Resize(segments_.GetSize());

  StdVector<UInt> candidates;
  UInt total = 0;

  for (UInt iSeg = 0; iSeg < segments_.GetSize(); ++iSeg) {

    FindElementCandidates(segments_[iSeg].secElem, candidates);

    graphPrimaries_[iSeg].Clear();
    for (UInt i = 0; i < candidates.GetSize(); ++i) {
      graphPrimaries_[iSeg].Push_back(primarySurfElems_[candidates[i]]);
    }
    total += graphPrimaries_[iSeg].GetSize();
  }

  LOG_DBG(contactinterface) << "contact pair '" << GetName() << "': reserved " << total
      << " (secondary, primary) blocks in the matrix graph for " << segments_.GetSize()
      << " segments, so the pairing may move without leaving the sparsity pattern";
}


bool ContactInterface::NeedsStablePointSet() const {
  return formulation_ == FORM_AUGMENTED_LAGRANGE
      || (frictionType_ == FRICTION_COULOMB && frictionMu_ > 0.0)
      || zeroInitialGap_;
}

Double ContactInterface::EvalFit(const Vector<Double>& coeff, const LocPoint& lp) const {

  Double v = coeff[0];
  for (UInt i = 1; i < coeff.GetSize(); ++i) {
    v += coeff[i] * lp.coord[i - 1];
  }
  return v;
}

void ContactInterface::FitWeightedLinear(const StdVector<ContactPoint>& pts,
                                         const StdVector<Double>& value,
                                         Vector<Double>& coeff) const {

  const UInt sdim = ptGrid_->GetDim() - 1;      // surface parameters: 1 in 2D, 2 in 3D
  const UInt nb   = sdim + 1;                   // basis (1, xi, eta)

  coeff.Resize(nb);
  coeff.Init(0.0);

  // Normal equations of the weighted least-squares problem, A c = r, with
  //   A_mn = SUM_k w_k b_m(xi_k) b_n(xi_k),   r_m = SUM_k w_k b_m(xi_k) f_k.
  Matrix<Double> A(nb, nb);
  Vector<Double> r(nb), b(nb);
  A.Init();
  r.Init(0.0);

  Double wSum = 0.0, wfSum = 0.0;

  for (UInt k = 0; k < pts.GetSize(); ++k) {

    const Double w = pts[k].weight;
    if (w <= 0.0) {
      continue;
    }

    b[0] = 1.0;
    for (UInt i = 0; i < sdim; ++i) {
      b[i + 1] = pts[k].secLocal.coord[i];
    }

    for (UInt m = 0; m < nb; ++m) {
      for (UInt n = 0; n < nb; ++n) {
        A[m][n] += w * b[m] * b[n];
      }
      r[m] += w * b[m] * value[k];
    }

    wSum  += w;
    wfSum += w * value[k];
  }

  if (wSum <= 0.0) {
    return;
  }

  // Cholesky with a rank guard. A pivot collapsing means the points do not determine a slope
  // in that direction -- an uncut segment in 2D has as few as two points, and a degenerate
  // sub-cell fewer still. Degrade to the weighted MEAN rather than invent a gradient: it is
  // the fit of the same family in the space the data does determine, and it keeps the
  // conservation property, which is the one thing this transfer must not lose.
  Matrix<Double> L(nb, nb);
  L.Init();
  bool ok = true;

  for (UInt m = 0; m < nb && ok; ++m) {
    Double d = A[m][m];
    for (UInt p = 0; p < m; ++p) {
      d -= L[m][p] * L[m][p];
    }
    // Relative to the leading entry, which is SUM w: this is a scale-free rank test.
    if (d <= 1e-12 * A[0][0]) {
      ok = false;
      break;
    }
    L[m][m] = std::sqrt(d);
    for (UInt n = m + 1; n < nb; ++n) {
      Double s = A[n][m];
      for (UInt p = 0; p < m; ++p) {
        s -= L[n][p] * L[m][p];
      }
      L[n][m] = s / L[m][m];
    }
  }

  if (!ok) {
    coeff[0] = wfSum / wSum;                    // weighted mean
    return;
  }

  Vector<Double> y(nb);
  for (UInt m = 0; m < nb; ++m) {
    Double s = r[m];
    for (UInt p = 0; p < m; ++p) {
      s -= L[m][p] * y[p];
    }
    y[m] = s / L[m][m];
  }
  for (UInt m = nb; m-- > 0; ) {
    Double s = y[m];
    for (UInt p = m + 1; p < nb; ++p) {
      s -= L[p][m] * coeff[p];
    }
    coeff[m] = s / L[m][m];
  }
}

void ContactInterface::SnapshotPointHistory(StdVector<PointHistoryFit>& fits) const {

  const UInt dim = ptGrid_->GetDim();

  fits.Resize(segments_.GetSize());

  StdVector<ContactPoint> pts;
  StdVector<Double> val;

  for (UInt iSeg = 0; iSeg < segments_.GetSize(); ++iSeg) {

    const ContactSegment& seg = segments_[iSeg];
    PointHistoryFit& fit = fits[iSeg];
    fit.valid = false;

    // Only points with a projection have a history worth carrying: an unprojected point never
    // had a multiplier, never slipped and never had its as-meshed gap frozen. Including them
    // would drag the fit towards zero over the part of the element that is not in contact.
    pts.Clear();
    for (UInt iP = 0; iP < seg.points.GetSize(); ++iP) {
      if (seg.points[iP].isProjected) {
        pts.Push_back(seg.points[iP]);
      }
    }
    if (pts.GetSize() == 0) {
      continue;
    }

    val.Resize(pts.GetSize());

    for (UInt k = 0; k < pts.GetSize(); ++k) {
      val[k] = pts[k].lambda;
    }
    FitWeightedLinear(pts, val, fit.lambda);

    for (UInt k = 0; k < pts.GetSize(); ++k) {
      val[k] = pts[k].gapOffset;
    }
    FitWeightedLinear(pts, val, fit.gapOffset);
    fit.gapOffsetSet = pts[0].gapOffsetSet;

    fit.slipPlastic.Resize(dim);
    fit.slipOffset.Resize(dim);
    for (UInt d = 0; d < dim; ++d) {

      for (UInt k = 0; k < pts.GetSize(); ++k) {
        val[k] = (pts[k].slipPlastic.GetSize() > d) ? pts[k].slipPlastic[d] : 0.0;
      }
      FitWeightedLinear(pts, val, fit.slipPlastic[d]);

      for (UInt k = 0; k < pts.GetSize(); ++k) {
        val[k] = (pts[k].slipOffset.GetSize() > d) ? pts[k].slipOffset[d] : 0.0;
      }
      FitWeightedLinear(pts, val, fit.slipOffset[d]);
    }

    fit.valid = true;
  }
}

void ContactInterface::RestorePointHistory(const StdVector<PointHistoryFit>& fits) {

  const UInt dim = ptGrid_->GetDim();

  if (fits.GetSize() != segments_.GetSize()) {
    EXCEPTION("ContactInterface '" << GetName() << "': the history snapshot has "
              << fits.GetSize() << " segments but the interface has " << segments_.GetSize()
              << ". A re-cut may change the POINTS of a segment, never the segments "
              << "themselves -- those are the secondary surface elements.");
  }

  for (UInt iSeg = 0; iSeg < segments_.GetSize(); ++iSeg) {

    const PointHistoryFit& fit = fits[iSeg];
    ContactSegment& seg = segments_[iSeg];

    if (!fit.valid) {
      continue;                                 // new points keep their zero-initialised state
    }

    for (UInt iP = 0; iP < seg.points.GetSize(); ++iP) {

      ContactPoint& cp = seg.points[iP];

      // lambda <= 0 by convention (it is a traction, not a pressure). A linear fit evaluated
      // outside the hull of the points it was fitted to can overshoot; clamping keeps the
      // Uzawa update on the admissible side, where the whole formulation lives.
      cp.lambda = std::min(0.0, EvalFit(fit.lambda, cp.secLocal));

      cp.gapOffset    = EvalFit(fit.gapOffset, cp.secLocal);
      cp.gapOffsetSet = fit.gapOffsetSet;

      cp.slipPlastic.Resize(dim);
      cp.slipOffset.Resize(dim);
      for (UInt d = 0; d < dim; ++d) {
        cp.slipPlastic[d] = EvalFit(fit.slipPlastic[d], cp.secLocal);
        cp.slipOffset[d]  = EvalFit(fit.slipOffset[d], cp.secLocal);
      }
    }
  }
}

void ContactInterface::ReSegment() {

  // Only meaningful where the partition is otherwise held fixed: without per-point history it
  // is re-cut every Newton iteration already, and there would be nothing for this to do but
  // repeat the last one.
  if (!partitionFrozenInSolve_ || !segmentsBuilt_ || !projected_ || !dispFct_) {
    return;
  }

  resegmentNow_ = true;
  UpdateProjection();
  resegmentNow_ = false;
}

void ContactInterface::RebuildSegmentPoints() {

  // Caller's contract (UpdateProjection): SLIDING_LARGE, mortar segmentation, a displacement
  // field, and the primary bounding boxes already refreshed for the current configuration --
  // BuildSegmentCells() searches candidates against them.
  numUnsegmented_ = 0;

  UInt numPoints = 0;
  for (UInt iEl = 0; iEl < segments_.GetSize(); ++iEl) {
    GenerateSegmentPoints(segments_[iEl]);
    numPoints += segments_[iEl].points.GetSize();
  }

  LOG_DBG(contactinterface) << "contact pair '" << GetName() << "': re-segmented into "
      << GetNumCells() << " mortar cells (" << numUnsegmented_ << " not segmented), "
      << numPoints << " integration points, integrated area = " << GetIntegratedArea();
}

// ================================================================================================
//  Bounding boxes and candidate search
// ================================================================================================

void ContactInterface::UpdatePrimaryBoundingBoxes() {

  const UInt dim = ptGrid_->GetDim();
  const UInt numEl = primarySurfElems_.GetSize();

  primaryBoxMin_.Resize(numEl);
  primaryBoxMax_.Resize(numEl);

  StdVector<Vector<Double> > nodeCoords;
  StdVector<UInt> nodes;

  for (UInt i = 0; i < numEl; ++i) {

    // Same configuration the projection uses, for the same reason: a box drawn around the
    // reference element would pre-select the wrong candidates for a query point that has
    // moved with the deformation. Under small sliding this is the reference configuration.
    Matrix<Double> ec;
    if (slidingType_ == SLIDING_LARGE && dispFct_) {
      GetDeformedElemCoords(primarySurfElems_[i], primaryVolRegion_, ec);
    } else {
      ptGrid_->GetElemNodesCoord(ec, primarySurfElems_[i]->connect, false);
    }
    const UInt numNodes = ec.GetNumCols();

    Vector<Double>& lo = primaryBoxMin_[i];
    Vector<Double>& hi = primaryBoxMax_[i];
    lo.Resize(dim);
    hi.Resize(dim);

    for (UInt d = 0; d < dim; ++d) {
      lo[d] = ec[d][0];
      hi[d] = ec[d][0];
    }
    for (UInt n = 1; n < numNodes; ++n) {
      for (UInt d = 0; d < dim; ++d) {
        lo[d] = std::min(lo[d], ec[d][n]);
        hi[d] = std::max(hi[d], ec[d][n]);
      }
    }

    // inflate, so that a point sitting exactly on the surface is not missed and so that
    // a separated but nearby point still finds its partner
    Double diag = 0.0;
    for (UInt d = 0; d < dim; ++d) {
      diag += (hi[d] - lo[d]) * (hi[d] - lo[d]);
    }
    diag = std::sqrt(diag);
    const Double pad = boundingBoxTolerance_ * diag + searchDistanceFactor_ * diag;
    for (UInt d = 0; d < dim; ++d) {
      lo[d] -= pad;
      hi[d] += pad;
    }
  }
}

void ContactInterface::FindCandidates(const Vector<Double>& point,
                                      StdVector<UInt>& candidates) const {

  const UInt dim = ptGrid_->GetDim();
  candidates.Clear();

  for (UInt i = 0; i < primaryBoxMin_.GetSize(); ++i) {
    bool inside = true;
    for (UInt d = 0; d < dim; ++d) {
      if (point[d] < primaryBoxMin_[i][d] || point[d] > primaryBoxMax_[i][d]) {
        inside = false;
        break;
      }
    }
    if (inside) {
      candidates.Push_back(i);
    }
  }

  if (candidates.GetSize() > 0) {
    return;
  }

  // Fallback: nearest box centre. A point that has no box hit is far from the primary
  // surface and will almost certainly be rejected by ProjectOntoElement, but returning
  // one candidate keeps the "no silent empty search" invariant and makes debugging easier.
  Double bestDist = std::numeric_limits<Double>::max();
  UInt bestIdx = 0;
  for (UInt i = 0; i < primaryBoxMin_.GetSize(); ++i) {
    Double d2 = 0.0;
    for (UInt d = 0; d < dim; ++d) {
      const Double c = 0.5 * (primaryBoxMin_[i][d] + primaryBoxMax_[i][d]);
      d2 += (point[d] - c) * (point[d] - c);
    }
    if (d2 < bestDist) {
      bestDist = d2;
      bestIdx = i;
    }
  }
  if (primaryBoxMin_.GetSize() > 0) {
    candidates.Push_back(bestIdx);
  }
}

// ================================================================================================
//  Geometry helpers
// ================================================================================================

void ContactInterface::EvalReferencePosition(SurfElem* elem, const LocPoint& lp,
                                             Vector<Double>& global) const {
  shared_ptr<ElemShapeMap> esm = ptGrid_->GetElemShapeMap(elem, false);
  global.Resize(ptGrid_->GetDim());
  esm->Local2Global(global, lp);
}

void ContactInterface::EvalCurrentPosition(SurfElem* elem, const LocPoint& lp,
                                           RegionIdType volRegion,
                                           Vector<Double>& global) const {

  EvalReferencePosition(elem, lp, global);

  // No displacement field yet -> we are still in the reference configuration. This is the
  // correct state for the very first pairing, which happens before any solve.
  if (!dispFct_) {
    return;
  }

  Vector<Double> u;
  EvalDisplacement(elem, lp, volRegion, u);

  const UInt dim = ptGrid_->GetDim();
  for (UInt d = 0; d < dim; ++d) {
    global[d] += u[d];
  }
}

void ContactInterface::EvalDisplacement(SurfElem* elem, const LocPoint& lp,
                                        RegionIdType volRegion,
                                        Vector<Double>& u) const {

  const UInt dim = ptGrid_->GetDim();
  u.Resize(dim);
  u.Init(0.0);

  if (!dispFct_) {
    return;
  }

  const std::set<RegionIdType>& volRegions =
      (volRegion == primaryVolRegion_) ? primaryVolRegionSet_ : secondaryVolRegionSet_;

  // Sampling the displacement field needs the surface point resolved onto its volume
  // neighbour, which is what SetWithSurface() does.
  shared_ptr<ElemShapeMap> esm = ptGrid_->GetElemShapeMap(elem, false);
  LocPointMapped lpm;
  lpm.SetWithSurface(lp, esm, volRegions, 0.0);

  dispFct_->GetVector(u, lpm);
}

void ContactInterface::GetDeformedElemCoords(SurfElem* elem, RegionIdType volRegion,
                                             Matrix<Double>& coords) const {

  // (sdim x numNodes), reference configuration to start from
  ptGrid_->GetElemNodesCoord(coords, elem->connect, false);

  if (!dispFct_) {
    return;                             // no solve yet: reference IS the current config
  }

  const UInt sdim = ptGrid_->GetDim();
  const ElemShape& shape = Elem::shapes[elem->type];
  LocPoint lp;
  lp.number = LocPoint::NOT_SET;
  Vector<Double> u;

  for (UInt n = 0; n < shape.numNodes; ++n) {
    lp.coord = shape.nodeCoords[n];
    EvalDisplacement(elem, lp, volRegion, u);
    for (UInt d = 0; d < sdim; ++d) {
      coords[d][n] += u[d];
    }
  }
}

shared_ptr<ElemShapeMap> ContactInterface::GetProjectionShapeMap(
    SurfElem* elem, RegionIdType volRegion) const {

  // Under small sliding the pairing is frozen at u = 0, so the reference configuration is
  // the configuration it was frozen in and there is nothing to deform.
  if (slidingType_ == SLIDING_SMALL || !dispFct_) {
    return ptGrid_->GetElemShapeMap(elem, false);
  }

  std::map<UInt, shared_ptr<ElemShapeMap> >::const_iterator it =
      projEsmCache_.find(elem->elemNum);
  if (it != projEsmCache_.end()) {
    return it->second;
  }

  Matrix<Double> coords;
  GetDeformedElemCoords(elem, volRegion, coords);

  shared_ptr<LagrangeElemShapeMap> esm(new LagrangeElemShapeMap(ptGrid_));
  esm->SetElem(elem, false);
  esm->SetElem(elem, coords);

  shared_ptr<ElemShapeMap> ret = esm;
  projEsmCache_[elem->elemNum] = ret;
  return ret;
}

Double ContactInterface::GetFacetOrientation(SurfElem* elem) const {

  std::map<UInt, Double>::const_iterator it = facetOrientation_.find(elem->elemNum);
  if (it != facetOrientation_.end()) {
    return it->second;
  }

  const UInt dim = ptGrid_->GetDim();
  shared_ptr<ElemShapeMap> esm = ptGrid_->GetElemShapeMap(elem, false);

  LocPoint mid;
  mid.coord  = Elem::shapes[elem->type].midPointCoord;
  mid.number = LocPoint::NOT_SET;

  Vector<Double> nRef(dim);
  esm->CalcNormal(nRef, mid);        // unit, out of ptVolElems[0] -> convention 2

  Matrix<Double> jac;
  esm->CalcJ(jac, mid);

  Vector<Double> nPar(dim);
  if (dim == 2) {
    nPar[0] =  jac[1][0];
    nPar[1] = -jac[0][0];
  } else {
    nPar[0] = jac[1][0] * jac[2][1] - jac[2][0] * jac[1][1];
    nPar[1] = jac[2][0] * jac[0][1] - jac[0][0] * jac[2][1];
    nPar[2] = jac[0][0] * jac[1][1] - jac[1][0] * jac[0][1];
  }

  Double dot = 0.0;
  for (UInt d = 0; d < dim; ++d) {
    dot += nPar[d] * nRef[d];
  }

  if (std::abs(dot) < 1e-30) {
    EXCEPTION("Contact pair '" << GetName() << "': primary surface element " << elem->elemNum
              << " has a degenerate parametrization -- its parametric normal is orthogonal to "
              << "its outward normal, so the contact normal orientation cannot be determined.");
  }

  const Double s = (dot > 0.0) ? 1.0 : -1.0;
  facetOrientation_[elem->elemNum] = s;
  return s;
}

void ContactInterface::EvalFacetNormalRaw(SurfElem* elem, const LocPoint& lp,
                                          Vector<Double>& normal) const {

  const UInt dim = ptGrid_->GetDim();

  // Deformed under large sliding: the contact normal is a property of the CURRENT primary
  // surface, and using a reference normal there is exactly the approximation Phase 8 removes.
  shared_ptr<ElemShapeMap> esm = GetProjectionShapeMap(elem, primaryVolRegion_);

  Matrix<Double> jac;
  esm->CalcJ(jac, lp);               // dim x (dim-1): the surface tangents a_i, deformed

  normal.Resize(dim);
  if (dim == 2) {
    normal[0] =  jac[1][0];
    normal[1] = -jac[0][0];
  } else {
    normal[0] = jac[1][0] * jac[2][1] - jac[2][0] * jac[1][1];
    normal[1] = jac[2][0] * jac[0][1] - jac[0][0] * jac[2][1];
    normal[2] = jac[0][0] * jac[1][1] - jac[1][0] * jac[0][1];
  }

  normal *= GetFacetOrientation(elem);
}

void ContactInterface::BuildNodalNormals() const {

  nodalNormals_.clear();

  const UInt dim = ptGrid_->GetDim();

  LocPoint lp;
  lp.number = LocPoint::NOT_SET;
  Vector<Double> n;

  for (UInt i = 0; i < primarySurfElems_.GetSize(); ++i) {

    SurfElem* elem = primarySurfElems_[i];
    const ElemShape& shape = Elem::shapes[elem->type];

    for (UInt a = 0; a < shape.numNodes; ++a) {

      lp.coord = shape.nodeCoords[a];
      EvalFacetNormalRaw(elem, lp, n);

      const UInt node = elem->connect[a];
      std::map<UInt, Vector<Double> >::iterator it = nodalNormals_.find(node);
      if (it == nodalNormals_.end()) {
        nodalNormals_[node] = n;
      } else {
        for (UInt d = 0; d < dim; ++d) {
          it->second[d] += n[d];
        }
      }
    }
  }

  for (std::map<UInt, Vector<Double> >::iterator it = nodalNormals_.begin();
       it != nodalNormals_.end(); ++it) {
    const Double len = it->second.NormL2();
    if (len > 1e-30) {
      it->second /= len;
    }
  }

  nodalNormalsBuilt_ = true;
}

shared_ptr<ElemShapeMap> ContactInterface::GetNormalShapeMap(SurfElem* elem) const {

  std::map<UInt, shared_ptr<ElemShapeMap> >::const_iterator it =
      normalEsmCache_.find(elem->elemNum);
  if (it != normalEsmCache_.end()) {
    return it->second;
  }

  if (!nodalNormalsBuilt_) {
    BuildNodalNormals();
  }

  const UInt dim = ptGrid_->GetDim();
  const ElemShape& shape = Elem::shapes[elem->type];

  // The "coordinates" of this map are normals, one column per node. Local2Global() then
  // computes coords_ * N(xi) = SUM_a N_a(xi) n_a, which is the interpolation wanted -- with
  // the element's own shape functions, so nothing can drift apart from the geometry.
  Matrix<Double> nodalN(dim, shape.numNodes);
  for (UInt a = 0; a < shape.numNodes; ++a) {
    std::map<UInt, Vector<Double> >::const_iterator nIt =
        nodalNormals_.find(elem->connect[a]);
    if (nIt == nodalNormals_.end()) {
      EXCEPTION("Contact pair '" << GetName() << "': node " << elem->connect[a]
                << " of primary surface element " << elem->elemNum
                << " has no averaged normal. BuildNodalNormals() and the primary surface "
                << "element list disagree.");
    }
    for (UInt d = 0; d < dim; ++d) {
      nodalN[d][a] = nIt->second[d];
    }
  }

  shared_ptr<LagrangeElemShapeMap> esm(new LagrangeElemShapeMap(ptGrid_));
  // Same two-step as GetProjectionShapeMap(): the coordinate overload reads ptElem_ from the
  // object and never assigns it, so the plain SetElem() has to run first.
  esm->SetElem(elem, false);
  esm->SetElem(elem, nodalN);

  shared_ptr<ElemShapeMap> ret = esm;
  normalEsmCache_[elem->elemNum] = ret;
  return ret;
}

void ContactInterface::ClearGeometryCaches() const {
  projEsmCache_.clear();
  normalEsmCache_.clear();
  nodalNormals_.clear();
  nodalNormalsBuilt_ = false;
}

void ContactInterface::EvalOutwardNormal(SurfElem* elem, const LocPoint& lp,
                                         Vector<Double>& normal) const {

  if (normalSmoothing_ == NORMAL_NODAL) {

    shared_ptr<ElemShapeMap> esm = GetNormalShapeMap(elem);
    esm->Local2Global(normal, lp);

    const Double len = normal.NormL2();
    if (len > 1e-30) {
      normal /= len;
      return;
    }
    // The interpolated normal can only cancel where the averaged nodal normals of one element
    // point in opposing directions, i.e. across a fold sharper than 90 degrees. Smoothing has
    // no meaning there; fall through to the facet normal, which at least still describes the
    // element the point was projected onto.
  }

  EvalFacetNormalRaw(elem, lp, normal);
  const Double len = normal.NormL2();
  if (len > 1e-30) {
    normal /= len;
  }
}

bool ContactInterface::IsInsideElement(Elem::ShapeType shape, const Vector<Double>& local,
                                       Double tol) const {
  switch (shape) {

    case Elem::ST_LINE:
      // reference domain [-1, 1]
      return (local[0] >= -1.0 - tol) && (local[0] <= 1.0 + tol);

    case Elem::ST_QUAD:
      // reference domain [-1, 1] x [-1, 1]
      return (local[0] >= -1.0 - tol) && (local[0] <= 1.0 + tol)
          && (local[1] >= -1.0 - tol) && (local[1] <= 1.0 + tol);

    case Elem::ST_TRIA:
      // reference domain { xi >= 0, eta >= 0, xi + eta <= 1 }
      return (local[0] >= -tol) && (local[1] >= -tol)
          && (local[0] + local[1] <= 1.0 + tol);

    default:
      EXCEPTION("ContactInterface: unsupported surface element shape '"
                << Elem::shapeType.ToString(shape) << "'. Supported: line, tria, quad.");
  }
  return false;
}

void ContactInterface::ClampToElement(Elem::ShapeType shape, Vector<Double>& local,
                                      Double slack) const {
  assert(slack > 0.0);

  switch (shape) {

    case Elem::ST_LINE:
      local[0] = std::max(-1.0 - slack, std::min(1.0 + slack, local[0]));
      break;

    case Elem::ST_QUAD:
      for (UInt d = 0; d < 2; ++d) {
        local[d] = std::max(-1.0 - slack, std::min(1.0 + slack, local[d]));
      }
      break;

    case Elem::ST_TRIA: {
      local[0] = std::max(-slack, local[0]);
      local[1] = std::max(-slack, local[1]);
      const Double s = local[0] + local[1];
      const Double smax = 1.0 + slack;
      if (s > smax) {
        local[0] *= smax / s;
        local[1] *= smax / s;
      }
      break;
    }

    default:
      EXCEPTION("ContactInterface: unsupported surface element shape '"
                << Elem::shapeType.ToString(shape) << "'. Supported: line, tria, quad.");
  }
}

// ================================================================================================
//  Closest-point projection
// ================================================================================================

bool ContactInterface::ClosestPointLocal(const Vector<Double>& point,
                                         SurfElem* elem,
                                         RegionIdType volRegion,
                                         Double clampSlack,
                                         LocPoint& localOut,
                                         bool& clamped) const {

  const UInt dim = ptGrid_->GetDim();
  const UInt sdim = dim - 1;                 // number of surface parameters
  const Elem::ShapeType shape = Elem::GetShapeType(elem->type);

  // Which configuration the projection runs in is decided by GetProjectionShapeMap():
  // reference under SLIDING_SMALL, where the pairing is frozen at u = 0 and reference and
  // current therefore coincide by construction; the DEFORMED surface under SLIDING_LARGE.
  //
  // This matters for the iteration below and not just for the answer's accuracy. The
  // stationarity condition it solves,
  //
  //     ( x_s - x_m(xi) ) . dx_m/dxi_i = 0
  //
  // involves the tangent of the surface being projected onto. Feeding it reference tangents
  // while measuring distances to a deformed surface would not converge slowly -- it would
  // converge to a different point. The deformed shape map supplies both x_m and dx_m/dxi
  // consistently, which is the whole reason it is built from deformed NODAL coordinates
  // rather than by correcting a reference map after the fact.
  shared_ptr<ElemShapeMap> esm = GetProjectionShapeMap(elem, volRegion);

  // start from the element midpoint
  localOut.coord = Elem::shapes[elem->type].midPointCoord;
  localOut.number = LocPoint::NOT_SET;
  clamped = false;

  Vector<Double> xm(dim), diff(dim), before(sdim);
  Matrix<Double> jac;

  // Gauss-Newton on  f(xi) = 0.5 * || point - x_m(xi) ||^2
  //   grad_i    = -( point - x_m ) . a_i           with a_i = dx_m/dxi_i
  //   Hess_ij  ~=  a_i . a_j                        (curvature term dropped)
  bool converged = false;

  for (UInt it = 0; it < projectionMaxIter_; ++it) {

    esm->Local2Global(xm, localOut);
    for (UInt d = 0; d < dim; ++d) {
      diff[d] = point[d] - xm[d];
    }

    esm->CalcJ(jac, localOut);

    // tangent vectors a_i are the columns of the Jacobian
    Vector<Double> grad(sdim);
    Matrix<Double> hess(sdim, sdim);

    for (UInt i = 0; i < sdim; ++i) {
      Double g = 0.0;
      for (UInt d = 0; d < dim; ++d) {
        g += diff[d] * jac[d][i];
      }
      grad[i] = g;
      for (UInt j = 0; j < sdim; ++j) {
        Double h = 0.0;
        for (UInt d = 0; d < dim; ++d) {
          h += jac[d][i] * jac[d][j];
        }
        hess[i][j] = h;
      }
    }

    // solve hess * delta = grad  (1x1 or 2x2)
    Vector<Double> delta(sdim);
    if (sdim == 1) {
      if (std::abs(hess[0][0]) < 1e-30) {
        return false;
      }
      delta[0] = grad[0] / hess[0][0];
    } else {
      const Double det = hess[0][0] * hess[1][1] - hess[0][1] * hess[1][0];
      if (std::abs(det) < 1e-30) {
        return false;
      }
      delta[0] = ( hess[1][1] * grad[0] - hess[0][1] * grad[1]) / det;
      delta[1] = (-hess[1][0] * grad[0] + hess[0][0] * grad[1]) / det;
    }

    for (UInt i = 0; i < sdim; ++i) {
      localOut.coord[i] += delta[i];
      before[i] = localOut.coord[i];
    }
    // Only a runaway guard -- the slack must stay well clear of the true element
    // boundary so that a point belonging to a neighbouring element is still allowed to
    // converge outside and be rejected by the caller.
    ClampToElement(shape, localOut.coord, clampSlack);
    // Only the last iteration decides: an early overshoot that the iteration then walks back
    // into the element is not a guard value.
    clamped = false;
    for (UInt i = 0; i < sdim; ++i) {
      if (std::abs(localOut.coord[i] - before[i]) > SEGMENTATION_MERGE_TOL) {
        clamped = true;
      }
    }

    if (delta.NormL2() < projectionTolerance_) {
      converged = true;
      break;
    }
  }

  return converged;
}

bool ContactInterface::ProjectOntoElement(const Vector<Double>& point,
                                          SurfElem* elem,
                                          LocPoint& localOut,
                                          Vector<Double>& globalOut,
                                          Vector<Double>& normalOut,
                                          Double& distOut) const {

  const UInt dim = ptGrid_->GetDim();
  const Elem::ShapeType shape = Elem::GetShapeType(elem->type);

  bool clamped = false;
  if (!ClosestPointLocal(point, elem, primaryVolRegion_, PROJECTION_CLAMP_SLACK,
                            localOut, clamped)) {
    return false;
  }

  shared_ptr<ElemShapeMap> esm = GetProjectionShapeMap(elem, primaryVolRegion_);
  Vector<Double> xm(dim);

  if (!IsInsideElement(shape, localOut.coord, PROJECTION_INSIDE_TOL)) {
    return false;
  }
  ClampToElement(shape, localOut.coord, 0.0 + 1e-14);

  esm->Local2Global(xm, localOut);
  globalOut = xm;
  EvalOutwardNormal(elem, localOut, normalOut);

  // convention 3:  g_N = (x_s - x_m) . n
  Double g = 0.0;
  for (UInt d = 0; d < dim; ++d) {
    g += (point[d] - xm[d]) * normalOut[d];
  }
  distOut = g;

  // Reject projections that are absurdly far away. Distance is measured against the
  // element diameter so the criterion is mesh-size independent.
  Vector<Double> diam;
  esm->CalcDiameter(diam);
  const Double elemSize = diam.NormL2();
  if (elemSize > 1e-30 && std::abs(g) > searchDistanceFactor_ * elemSize) {
    return false;
  }

  return true;
}

// ================================================================================================
//  Projection driver
// ================================================================================================

void ContactInterface::UpdateProjection() {

  if (!segmentsBuilt_) {
    EXCEPTION("ContactInterface '" << GetName()
              << "': UpdateProjection() called before BuildSegments().");
  }

  // Under small-sliding kinematics the pairing and the normals are frozen after the first
  // projection -- only the gap follows the deformation. This is the whole point of the
  // small-sliding assumption and is what keeps a Newton iteration cheap.
  if (projected_ && slidingType_ == SLIDING_SMALL) {
    UpdateGaps();
    return;
  }

  if (HasFriction()) {
    AnchorTangentialSlip();
  }

  ClearGeometryCaches();

  UpdatePrimaryBoundingBoxes();

  // Under large sliding the mortar partition goes stale exactly as the pairing does: its
  // sub-cell boundaries are the projected outlines of the primary elements, and those move
  // with the body. Rebuilding it here, right after the bounding boxes and before the
  // projection, means the points that get projected below are already the ones belonging to
  // the current partition -- so the whole update is one pass and no point is ever projected
  // twice.
  if (slidingType_ == SLIDING_LARGE && mortarSegmentation_ && dispFct_ && projected_
      && (!partitionFrozenInSolve_ || resegmentNow_)) {

    if (NeedsStablePointSet()) {
      // Carry the history across the point set that is about to be replaced. Cheap: one
      // small least-squares fit per secondary element per field, and only at step boundaries.
      StdVector<PointHistoryFit> fits;
      SnapshotPointHistory(fits);
      RebuildSegmentPoints();
      RestorePointHistory(fits);
    } else {
      RebuildSegmentPoints();
    }
  }

  const UInt dim = ptGrid_->GetDim();

  StdVector<UInt> candidates;
  LocPoint bestLocal, trialLocal;
  Vector<Double> bestGlobal(dim), trialGlobal(dim);
  Vector<Double> bestNormal(dim), trialNormal(dim);

  UInt pointIdx = 0;
  UInt numProjected = 0;
  Vector<Double> refPoint(dim);

  for (UInt iSeg = 0; iSeg < segments_.GetSize(); ++iSeg) {

    ContactSegment& seg = segments_[iSeg];

    for (UInt iP = 0; iP < seg.points.GetSize(); ++iP, ++pointIdx) {

      ContactPoint& cp = seg.points[iP];

      // The query point must live in the same configuration as the surface it is projected
      // onto: reference under small sliding, CURRENT under large sliding. Mixing the two
      // would measure the gap between a reference point and a deformed surface, which is
      // neither quantity.
      if (slidingType_ == SLIDING_LARGE) {
        EvalCurrentPosition(seg.secElem, cp.secLocal, secondaryVolRegion_, refPoint);
      } else {
        EvalReferencePosition(seg.secElem, cp.secLocal, refPoint);
      }

      FindCandidates(refPoint, candidates);

      // Among all candidates keep the one whose projection is closest, i.e. smallest |gap|.
      Double bestAbsDist = std::numeric_limits<Double>::max();
      Double bestDist = 0.0;
      SurfElem* bestElem = nullptr;

      for (UInt iC = 0; iC < candidates.GetSize(); ++iC) {

        SurfElem* cand = primarySurfElems_[candidates[iC]];
        Double trialDist = 0.0;

        if (!ProjectOntoElement(refPoint, cand, trialLocal, trialGlobal,
                                trialNormal, trialDist)) {
          continue;
        }

        if (std::abs(trialDist) < bestAbsDist) {
          bestAbsDist = std::abs(trialDist);
          bestDist    = trialDist;
          bestElem    = cand;
          bestLocal   = trialLocal;
          bestGlobal  = trialGlobal;
          bestNormal  = trialNormal;
        }
      }

      if (bestElem != nullptr) {
        cp.primElem    = bestElem;
        cp.primLocal   = bestLocal;
        cp.primGlobal  = bestGlobal;
        cp.normal      = bestNormal;
        cp.gap         = bestDist;      // reference-config value, overwritten below
        cp.isProjected = true;
        ++numProjected;
      } else {
        cp.primElem    = nullptr;
        cp.isProjected = false;
        cp.gap         = 0.0;
        cp.isActive    = false;
      }
    }
  }

  projected_ = true;

  if (HasFriction()) {
    RebaseTangentialSlip();
  }

  UpdateGaps();

  // Freeze the as-meshed gap so that an imperfectly meshed assembly starts from a state
  // that is exactly touching rather than pre-stressed. Done after the first UpdateGaps()
  // so the offset is captured from the same quantity it will later be subtracted from.
  //
  // The offset lives on the POINT and is guarded by a per-point flag rather than by "is this
  // the first projection", so that it cannot be mismatched to the wrong point. That mattered:
  // it used to be a parallel array indexed by a running point counter, which is only safe as
  // long as the point set never changes -- and under SLIDING_LARGE the mortar partition now
  // rebuilds it (RebuildSegmentPoints()).
  //
  // The flag is normally only false on the first projection, where the displacement field
  // does not exist yet and cp.gap IS the as-meshed gap. A point created by a step-boundary
  // re-cut is not such a case and must not be treated as one -- re-deriving its offset from
  // the current (deformed) gap would freeze the deformation into it. RestorePointHistory()
  // therefore hands the new point both the interpolated offset and the flag, so it arrives
  // here already set and this loop passes it over. That is the whole reason gapOffsetSet is
  // transferred as a boolean rather than recomputed.
  if (zeroInitialGap_) {

    bool anyNew = false;

    for (UInt iSeg = 0; iSeg < segments_.GetSize(); ++iSeg) {
      ContactSegment& seg = segments_[iSeg];
      for (UInt iP = 0; iP < seg.points.GetSize(); ++iP) {

        ContactPoint& cp = seg.points[iP];
        if (cp.gapOffsetSet) {
          continue;
        }
        cp.gapOffsetSet = true;
        anyNew = true;
        cp.gapOffset = cp.isProjected ? cp.gap : 0.0;
      }
    }

    if (anyNew) {
      UpdateGaps();
    }
  }

  // Freeze the reference gap: g_N = gap0 + (u_s - u_m).n from here on. At this moment the
  // displacement is whatever it was during this projection, so gap0 is captured as the gap
  // MINUS the displacement part, keeping the decomposition exact even for a re-projection
  // that happens at u != 0.
  {
    Vector<Double> us, um;
    for (UInt iSeg = 0; iSeg < segments_.GetSize(); ++iSeg) {
      ContactSegment& seg = segments_[iSeg];
      for (UInt iP = 0; iP < seg.points.GetSize(); ++iP) {
        ContactPoint& cp = seg.points[iP];
        if (!cp.isProjected) {
          cp.gap0 = 0.0;
          continue;
        }
        EvalDisplacement(cp.secElem, cp.secLocal, secondaryVolRegion_, us);
        EvalDisplacement(cp.primElem, cp.primLocal, primaryVolRegion_, um);
        Double du = 0.0;
        for (UInt d = 0; d < dim; ++d) {
          du += (us[d] - um[d]) * cp.normal[d];
        }
        cp.gap0 = cp.gap - du;
      }
    }
  }

  BuildElementList();

  LOG_DBG(contactinterface) << "contact pair '" << GetName() << "': projected "
      << numProjected << " of " << pointIdx << " points, "
      << GetNumActivePoints() << " active, min gap = " << GetMinGap();
}


// ================================================================================================
//  Assembler-facing element list
// ================================================================================================

void ContactInterface::BuildElementList() {

  elemList_->Clear();

  for (UInt iSeg = 0; iSeg < segments_.GetSize(); ++iSeg) {

    ContactSegment& seg = segments_[iSeg];
    StdVector<SurfElem*> prims;
    StdVector<StdVector<UInt> > groups;

    for (UInt iP = 0; iP < seg.points.GetSize(); ++iP) {

      const ContactPoint& cp = seg.points[iP];
      if (!cp.isProjected || cp.primElem == nullptr) {
        continue;
      }

      UInt g = prims.GetSize();
      for (UInt k = 0; k < prims.GetSize(); ++k) {
        if (prims[k] == cp.primElem) {
          g = k;
          break;
        }
      }
      if (g == prims.GetSize()) {
        prims.Push_back(cp.primElem);
        groups.Push_back(StdVector<UInt>());
      }
      groups[g].Push_back(iP);
    }

    // --------------------------------------------------------------------------------
    //  Graph reservation (large sliding, force-carrying pairs only)
    // --------------------------------------------------------------------------------
    // ReserveGraphPrimaries() recorded, once and on the reference configuration, every
    // primary element this secondary element could plausibly reach. Emitting an assembler
    // element for each of them -- with an EMPTY point group where nothing projects there --
    // makes the list independent of the current pairing, so the sparsity pattern built from
    // it at setup still covers the pairing after the body has slid.
    if (iSeg < graphPrimaries_.GetSize()) {
      for (UInt c = 0; c < graphPrimaries_[iSeg].GetSize(); ++c) {
        SurfElem* cand = graphPrimaries_[iSeg][c];
        bool have = false;
        for (UInt k = 0; k < prims.GetSize(); ++k) {
          if (prims[k] == cand) {
            have = true;
            break;
          }
        }
        if (!have) {
          prims.Push_back(cand);
          groups.Push_back(StdVector<UInt>());
        }
      }
    }

    for (UInt g = 0; g < prims.GetSize(); ++g) {

      shared_ptr<ContactNcSurfElem> ce(new ContactNcSurfElem());
      ce->segment     = &seg;
      ce->ptSecondary = seg.secElem;
      ce->ptPrimary   = prims[g];
      ce->pointIdx    = groups[g];

      ce->elemNum  = seg.secElem->elemNum;
      ce->type     = seg.secElem->type;
      ce->regionId = seg.secElem->regionId;
      ce->connect  = seg.secElem->connect;
      ce->ptVolElems[0] = seg.secElem->ptVolElems[0];
      ce->ptVolElems[1] = nullptr;

      elemList_->AddElement(ce);

      // ------------------------------------------------------------------------------
      //  The matrix graph is built ONCE, from this list.
      // ------------------------------------------------------------------------------
      // A (secondary, primary) pair that appears only later has no reserved entries in the
      // sparsity pattern, and the assembler aborts inside CRS_Matrix with
      //     AddToMatrixEntry: Index pair = (i , j) not found
      // which says nothing about contact.
      //
      // ReserveGraphPrimaries() is what normally prevents this, by putting the whole
      // reference-configuration candidate set into the list up front. This check is the
      // backstop for the case it cannot cover: a body that slides beyond the inflated
      // candidate boxes entirely, which is a genuine "the search radius was too small" and
      // has a genuine remedy. Diagnose it here, where the cause is known.
      const std::pair<UInt, UInt> key(seg.secElem->elemNum, prims[g]->elemNum);

      if (!elemPairsFixed_) {
        elemPairs_.insert(key);
      } else if (applyForce_ && elemPairs_.find(key) == elemPairs_.end()) {
        EXCEPTION("Contact pair '" << GetName() << "': secondary element "
                  << seg.secElem->elemNum << " has slid onto primary element "
                  << prims[g]->elemNum << ", which lies outside the candidate set reserved "
                  << "in the matrix sparsity pattern at setup, so there is nowhere to "
                  << "assemble its contact stiffness. The reservation covers every primary "
                  << "element within searchDistanceFactor (currently "
                  << searchDistanceFactor_ << ") element diameters of the secondary element "
                  << "in the REFERENCE configuration; this body has slid further. Raise "
                  << "searchDistanceFactor so the reserved set covers the primary elements "
                  << "the body will reach, or reduce the sliding per step with a transient "
                  << "or ramped run.");
      }
    }
  }

  // Everything built before the first solve belongs to the pattern the graph is built from.
  elemPairsFixed_ = true;

  LOG_DBG(contactinterface) << "contact pair '" << GetName() << "': element list has "
      << elemList_->GetSize() << " contact segments";
}


// ================================================================================================
//  Per-point gap against the current displacement
// ================================================================================================

void ContactInterface::EvalRelativeDisplacement(const ContactPoint& cp,
                                                Vector<Double>& du) const {

  const UInt dim = ptGrid_->GetDim();
  du.Resize(dim);
  du.Init(0.0);

  if (!cp.isProjected) {
    return;
  }

  Vector<Double> us(dim), um(dim);
  EvalDisplacement(cp.secElem, cp.secLocal, secondaryVolRegion_, us);
  EvalDisplacement(cp.primElem, cp.primLocal, primaryVolRegion_, um);

  for (UInt d = 0; d < dim; ++d) {
    du[d] = us[d] - um[d];
  }
}

Double ContactInterface::EvalCurrentGap(const ContactPoint& cp) const {

  if (!cp.isProjected) {
    return 0.0;
  }

  // g_N = gap0 + (u_s - u_m).n   -- exactly the decomposition the contact stiffness relies
  // on, see ContactABInt.hh. Evaluating the displacement part directly (rather than
  // differencing two current positions) keeps this consistent with the integrators by
  // construction.
  const UInt dim = ptGrid_->GetDim();
  Vector<Double> du;
  EvalRelativeDisplacement(cp, du);

  Double g = cp.gap0;
  for (UInt d = 0; d < dim; ++d) {
    g += du[d] * cp.normal[d];
  }
  return g;
}

void ContactInterface::EvalCurrentSlip(const ContactPoint& cp, Vector<Double>& gT) const {

  const UInt dim = ptGrid_->GetDim();
  Vector<Double> du;
  EvalRelativeDisplacement(cp, du);

  // g_T = slipOffset + P_T (u_s - u_m) = slipOffset + du - (du.n) n,
  // the tangential twin of  g_N = gap0 + (u_s - u_m).n.
  //
  // The closest-point projection makes the relative position purely normal at the moment the
  // pairing is made, so P_T (u_s - u_m) is the whole slip FROM THAT MOMENT. slipOffset carries
  // what came before it, and is identically zero whenever the pairing has never moved -- so
  // this reduces to the original expression under SLIDING_SMALL exactly, not approximately.
  // See ContactPoint::slipOffset and RebaseTangentialSlip().
  Double dn = 0.0;
  for (UInt d = 0; d < dim; ++d) {
    dn += du[d] * cp.normal[d];
  }

  gT.Resize(dim);
  const bool hasOffset = (cp.slipOffset.GetSize() == dim);
  for (UInt d = 0; d < dim; ++d) {
    gT[d] = du[d] - dn * cp.normal[d] + (hasOffset ? cp.slipOffset[d] : 0.0);
  }
}

void ContactInterface::AnchorTangentialSlip() {

  const UInt dim = ptGrid_->GetDim();
  Vector<Double> gT;

  for (UInt iSeg = 0; iSeg < segments_.GetSize(); ++iSeg) {
    ContactSegment& seg = segments_[iSeg];
    for (UInt iP = 0; iP < seg.points.GetSize(); ++iP) {

      ContactPoint& cp = seg.points[iP];
      if (!cp.isProjected) {
        cp.slipAnchored = false;
        continue;
      }

      EvalCurrentSlip(cp, gT);              // uses the pairing that is about to be discarded
      cp.slipOffset.Resize(dim);
      for (UInt d = 0; d < dim; ++d) {
        cp.slipOffset[d] = gT[d];
      }
      cp.slipAnchored = true;
    }
  }
}

void ContactInterface::RebaseTangentialSlip() {

  const UInt dim = ptGrid_->GetDim();
  Vector<Double> du;

  for (UInt iSeg = 0; iSeg < segments_.GetSize(); ++iSeg) {
    ContactSegment& seg = segments_[iSeg];
    for (UInt iP = 0; iP < seg.points.GetSize(); ++iP) {

      ContactPoint& cp = seg.points[iP];

      cp.slipOffset.Resize(dim);
      if (cp.slipPlastic.GetSize() != dim) {
        cp.slipPlastic.Resize(dim);
      }

      if (!cp.isProjected) {
        cp.slipAnchored = false;
        continue;
      }

      if (!cp.slipAnchored) {
        Double oN = 0.0, pN = 0.0;
        for (UInt d = 0; d < dim; ++d) {
          oN += cp.slipOffset[d] * cp.normal[d];
          pN += cp.slipPlastic[d] * cp.normal[d];
        }
        for (UInt d = 0; d < dim; ++d) {
          cp.slipOffset[d]  -= oN * cp.normal[d];
          cp.slipPlastic[d] -= pN * cp.normal[d];
        }
        continue;
      }
      cp.slipAnchored = false;

      EvalRelativeDisplacement(cp, du);     // against the NEW pairing

      // slipOffset <- P_T (anchor - du),  slipPlastic <- P_T slipPlastic, with the NEW normal.
      Double aN = 0.0, pN = 0.0;
      for (UInt d = 0; d < dim; ++d) {
        aN += (cp.slipOffset[d] - du[d]) * cp.normal[d];
        pN += cp.slipPlastic[d] * cp.normal[d];
      }
      for (UInt d = 0; d < dim; ++d) {
        cp.slipOffset[d]  = (cp.slipOffset[d] - du[d]) - aN * cp.normal[d];
        cp.slipPlastic[d] = cp.slipPlastic[d] - pN * cp.normal[d];
      }
    }
  }
}

bool ContactInterface::EvalFrictionTraction(const ContactSegment& seg,
                                            const ContactPoint& cp,
                                            Double gap,
                                            const Vector<Double>& gT,
                                            Vector<Double>& tracT) const {

  const UInt dim = ptGrid_->GetDim();
  tracT.Resize(dim);
  tracT.Init(0.0);

  if (!HasFriction() || !cp.isProjected) {
    return false;
  }

  // An open point is not in contact, so there is nothing for friction to act through.
  if (!IsGapActive(seg, cp, gap)) {
    return false;
  }

  const Double epsT = GetTangentialPenalty(seg);
  if (epsT <= 0.0) {
    return false;
  }

  // Elastic (trial) slip and trial traction. slipPlastic is HISTORY -- constant within a
  // Newton solve -- which is what makes this a well-defined function of u.
  Double normTrial = 0.0;
  for (UInt d = 0; d < dim; ++d) {
    tracT[d] = epsT * (gT[d] - cp.slipPlastic[d]);
    normTrial += tracT[d] * tracT[d];
  }
  normTrial = std::sqrt(normTrial);

  // Coulomb limit. The pressure is the same augmented traction the normal part uses, so a
  // pair under augmented Lagrange bounds friction by the augmented pressure automatically.
  const Double p = -EvalTraction(seg, cp, gap);
  const Double limit = frictionMu_ * p;

  if (normTrial <= limit) {
    return true;                        // stick: the trial traction is admissible as it is
  }

  // Slip: project the trial traction radially back onto the Coulomb cone. The DIRECTION is
  // that of the trial traction, which is why normTrial > limit >= 0 guarantees it is well
  // defined here -- a zero trial traction can never reach this branch.
  const Double scale = limit / normTrial;
  for (UInt d = 0; d < dim; ++d) {
    tracT[d] *= scale;
  }
  return false;
}

void ContactInterface::UpdateSlipHistory() {

  if (!HasFriction()) {
    return;
  }

  const UInt dim = ptGrid_->GetDim();
  Vector<Double> gT, tracT;


  for (UInt iSeg = 0; iSeg < segments_.GetSize(); ++iSeg) {

    ContactSegment& seg = segments_[iSeg];
    const Double epsT = GetTangentialPenalty(seg);
    if (epsT <= 0.0) {
      continue;
    }

    for (UInt iP = 0; iP < seg.points.GetSize(); ++iP) {

      ContactPoint& cp = seg.points[iP];
      if (!cp.isProjected) {
        continue;
      }

      const Double gap = EvalCurrentGap(cp);
      EvalCurrentSlip(cp, gT);

      if (!IsGapActive(seg, cp, gap)) {
        // Not in contact: re-engage from zero elastic slip when it closes again, rather than
        // from whatever it happened to carry when it separated.
        for (UInt d = 0; d < dim; ++d) {
          cp.slipPlastic[d] = gT[d];
        }
        continue;
      }

      EvalFrictionTraction(seg, cp, gap, gT, tracT);

      // g_T^p <- g_T - t_T/eps_T. A no-op while the point sticks, because there
      // t_T = eps_T (g_T - g_T^p) identically; it absorbs the slipped part once it does not.
      for (UInt d = 0; d < dim; ++d) {
        cp.slipPlastic[d] = gT[d] - tracT[d] / epsT;
      }
    }
  }

  LOG_DBG(contactinterface) << "contact pair '" << GetName() << "': slip history advanced, "
      << GetNumStickPoints() << " stick / " << GetNumSlipPoints() << " slip";
}

void ContactInterface::ResetSlipHistory() {
  for (UInt iSeg = 0; iSeg < segments_.GetSize(); ++iSeg) {
    ContactSegment& seg = segments_[iSeg];
    for (UInt iP = 0; iP < seg.points.GetSize(); ++iP) {
      seg.points[iP].slipPlastic.Init(0.0);
    }
  }
}

void ContactInterface::UpdateGaps() {

  const UInt dim = ptGrid_->GetDim();

  for (UInt iSeg = 0; iSeg < segments_.GetSize(); ++iSeg) {

    ContactSegment& seg = segments_[iSeg];

    for (UInt iP = 0; iP < seg.points.GetSize(); ++iP) {

      ContactPoint& cp = seg.points[iP];
      if (!cp.isProjected) {
        cp.isActive = false;
        continue;
      }

      // both bodies deform, so both positions have to be re-evaluated
      EvalCurrentPosition(seg.secElem, cp.secLocal, secondaryVolRegion_, cp.secGlobal);
      EvalCurrentPosition(cp.primElem, cp.primLocal, primaryVolRegion_, cp.primGlobal);

      Double g = 0.0;
      for (UInt d = 0; d < dim; ++d) {
        g += (cp.secGlobal[d] - cp.primGlobal[d]) * cp.normal[d];
      }
      cp.gap = g - cp.gapOffset;
      // KKT-consistent activation, with a tolerance so the exactly-touching case is not
      // decided by rounding -- see IsGapActive(). The penalty has to exist before this,
      // because under augmented Lagrange the test weighs lambda against eps_N * g.
      seg.penalty = (seg.elemSize > 0.0)
          ? normalPenalty_ * refModulus_ / seg.elemSize
          : 0.0;
      cp.isActive = IsGapActive(seg, cp, cp.gap);
    }
  }

  UpdateTractions();
}

// ================================================================================================
//  Contact traction
// ================================================================================================

void ContactInterface::UpdateTractions() {

  for (UInt iSeg = 0; iSeg < segments_.GetSize(); ++iSeg) {

    ContactSegment& seg = segments_[iSeg];

    // seg.elemSize is set once in BuildSegments(); it has to exist before the first gap is
    // evaluated because the activation tolerance is relative to it.
    seg.penalty = (seg.elemSize > 0.0)
        ? normalPenalty_ * refModulus_ / seg.elemSize
        : 0.0;

    for (UInt iP = 0; iP < seg.points.GetSize(); ++iP) {
      ContactPoint& cp = seg.points[iP];
      // convention: t_N = <lambda + eps_N * g_N>_- <= 0, and the pressure is p = -t_N >= 0.
      // Under FORM_PENALTY lambda is zero and this is eps_N * <g_N>_- exactly as before.
      cp.pressure = cp.isActive ? (-EvalTraction(seg, cp, cp.gap)) : 0.0;
    }
  }

  UpdateFrictionState();
}

void ContactInterface::UpdateFrictionState() {

  if (!HasFriction()) {
    return;
  }

  Vector<Double> gT;

  for (UInt iSeg = 0; iSeg < segments_.GetSize(); ++iSeg) {

    ContactSegment& seg = segments_[iSeg];

    for (UInt iP = 0; iP < seg.points.GetSize(); ++iP) {

      ContactPoint& cp = seg.points[iP];
      if (!cp.isProjected || !cp.isActive) {
        cp.tracTang.Init(0.0);
        cp.isSticking = false;
        continue;
      }

      EvalCurrentSlip(cp, gT);
      cp.isSticking = EvalFrictionTraction(seg, cp, cp.gap, gT, cp.tracTang);
    }
  }
}

void ContactInterface::UpdateMultipliers(Double& maxDelta, Double& maxLambda) {

  maxDelta = 0.0;
  maxLambda = 0.0;

  for (UInt iSeg = 0; iSeg < segments_.GetSize(); ++iSeg) {

    ContactSegment& seg = segments_[iSeg];

    for (UInt iP = 0; iP < seg.points.GetSize(); ++iP) {

      ContactPoint& cp = seg.points[iP];
      if (!cp.isProjected) {
        continue;
      }

      // Uzawa:  lambda <- < lambda + eps_N * g_N >_-
      //
      // The Macaulay bracket is what makes this an inequality method: a point whose augmented
      // traction has gone positive is separating, and its multiplier is reset to zero rather
      // than allowed to become a tensile (sticking) traction.
      const Double before = cp.lambda;
      cp.lambda = std::min(cp.lambda + seg.penalty * cp.gap, 0.0);

      maxDelta = std::max(maxDelta, std::abs(cp.lambda - before));
      maxLambda = std::max(maxLambda, std::abs(cp.lambda));
    }
  }

  UpdateTractions();
  for (UInt iSeg = 0; iSeg < segments_.GetSize(); ++iSeg) {
    ContactSegment& seg = segments_[iSeg];
    for (UInt iP = 0; iP < seg.points.GetSize(); ++iP) {
      ContactPoint& cp = seg.points[iP];
      cp.isActive = cp.isProjected && IsGapActive(seg, cp, cp.gap);
    }
  }

  LOG_DBG(contactinterface) << "contact pair '" << GetName() << "': Uzawa update, "
      << "maxDelta = " << maxDelta << ", maxLambda = " << maxLambda
      << ", minGap = " << GetMinGap();
}

void ContactInterface::ResetMultipliers() {
  for (UInt iSeg = 0; iSeg < segments_.GetSize(); ++iSeg) {
    ContactSegment& seg = segments_[iSeg];
    for (UInt iP = 0; iP < seg.points.GetSize(); ++iP) {
      seg.points[iP].lambda = 0.0;
    }
  }
}

// ================================================================================================
//  Diagnostics
// ================================================================================================

UInt ContactInterface::GetNumProjectedPoints() const {
  UInt n = 0;
  for (UInt i = 0; i < segments_.GetSize(); ++i) {
    for (UInt j = 0; j < segments_[i].points.GetSize(); ++j) {
      if (segments_[i].points[j].isProjected) ++n;
    }
  }
  return n;
}

UInt ContactInterface::GetNumActivePoints() const {
  UInt n = 0;
  for (UInt i = 0; i < segments_.GetSize(); ++i) {
    for (UInt j = 0; j < segments_[i].points.GetSize(); ++j) {
      if (segments_[i].points[j].isActive) ++n;
    }
  }
  return n;
}

Double ContactInterface::GetMinGap() const {
  Double g = 0.0;
  bool any = false;
  for (UInt i = 0; i < segments_.GetSize(); ++i) {
    for (UInt j = 0; j < segments_[i].points.GetSize(); ++j) {
      const ContactPoint& cp = segments_[i].points[j];
      if (!cp.isProjected) continue;
      if (!any || cp.gap < g) {
        g = cp.gap;
        any = true;
      }
    }
  }
  return g;
}

Double ContactInterface::GetTotalContactForce() const {
  Double f = 0.0;
  for (UInt i = 0; i < segments_.GetSize(); ++i) {
    for (UInt j = 0; j < segments_[i].points.GetSize(); ++j) {
      f += segments_[i].points[j].weight * segments_[i].points[j].pressure;
    }
  }
  return f;
}

void ContactInterface::GetContactForceResultant(Vector<Double>& f) const {

  const UInt dim = ptGrid_->GetDim();
  f.Resize(dim);
  f.Init(0.0);

  for (UInt i = 0; i < segments_.GetSize(); ++i) {
    for (UInt j = 0; j < segments_[i].points.GetSize(); ++j) {

      const ContactPoint& cp = segments_[i].points[j];

      // An inactive point carries pressure 0 and traction 0, so it contributes nothing and
      // needs no separate test -- but its normal may be stale or unset, and 0 * garbage is
      // not reliably 0. Skip it explicitly.
      if (!cp.isProjected || cp.normal.GetSize() < dim) continue;

      for (UInt d = 0; d < dim; ++d) {
        f[d] += cp.weight * cp.pressure * cp.normal[d];
      }
      if (cp.tracTang.GetSize() >= dim) {
        // MINUS, and the two terms above and below genuinely differ in sign convention.
        // cp.tracTang is the elastic restoring traction eps_T*(g_T - g_T^p), so it carries the
        // sign of the SLIP; the force it exerts on the secondary body opposes that slip, which
        // is why ContactLinInt assembles it as -rowSign*w*(N^T t_T). cp.pressure*normal is
        // already the force on the secondary (+n points from the primary towards it), so only
        // this term needs negating to put both in the documented convention.
        //
        // Getting this wrong is invisible on a FLAT interface -- t_T is then purely tangential
        // and leaks nothing into the component equilibrium constrains
        for (UInt d = 0; d < dim; ++d) {
          f[d] -= cp.weight * cp.tracTang[d];
        }
      }
    }
  }
}

Double ContactInterface::GetMaxPressure() const {
  Double p = 0.0;
  for (UInt i = 0; i < segments_.GetSize(); ++i) {
    for (UInt j = 0; j < segments_[i].points.GetSize(); ++j) {
      p = std::max(p, segments_[i].points[j].pressure);
    }
  }
  return p;
}

Double ContactInterface::GetReferenceGapForce() const {
  Double f = 0.0;
  for (UInt i = 0; i < segments_.GetSize(); ++i) {
    const ContactSegment& seg = segments_[i];
    for (UInt j = 0; j < seg.points.GetSize(); ++j) {
      const ContactPoint& cp = seg.points[j];
      if (cp.isProjected && IsGapActive(seg, cp, cp.gap0)
          && std::abs(cp.gap0) > activationTolerance_ * seg.elemSize) {
        f += seg.penalty * std::abs(cp.gap0) * cp.weight;
      }
    }
  }
  return f;
}

UInt ContactInterface::GetNumStickPoints() const {
  UInt n = 0;
  for (UInt i = 0; i < segments_.GetSize(); ++i) {
    for (UInt j = 0; j < segments_[i].points.GetSize(); ++j) {
      const ContactPoint& cp = segments_[i].points[j];
      if (cp.isActive && cp.isSticking) ++n;
    }
  }
  return n;
}

UInt ContactInterface::GetNumSlipPoints() const {
  UInt n = 0;
  for (UInt i = 0; i < segments_.GetSize(); ++i) {
    for (UInt j = 0; j < segments_[i].points.GetSize(); ++j) {
      const ContactPoint& cp = segments_[i].points[j];
      if (cp.isActive && !cp.isSticking) ++n;
    }
  }
  return n;
}

Double ContactInterface::GetMaxTangentialTraction() const {
  const UInt dim = ptGrid_->GetDim();
  Double t = 0.0;
  for (UInt i = 0; i < segments_.GetSize(); ++i) {
    for (UInt j = 0; j < segments_[i].points.GetSize(); ++j) {
      const ContactPoint& cp = segments_[i].points[j];
      if (cp.tracTang.GetSize() < dim) continue;
      Double n2 = 0.0;
      for (UInt d = 0; d < dim; ++d) {
        n2 += cp.tracTang[d] * cp.tracTang[d];
      }
      t = std::max(t, std::sqrt(n2));
    }
  }
  return t;
}

Double ContactInterface::GetTotalTangentialForce() const {
  const UInt dim = ptGrid_->GetDim();
  Vector<Double> f(dim);
  f.Init(0.0);
  for (UInt i = 0; i < segments_.GetSize(); ++i) {
    for (UInt j = 0; j < segments_[i].points.GetSize(); ++j) {
      const ContactPoint& cp = segments_[i].points[j];
      if (cp.tracTang.GetSize() < dim) continue;
      for (UInt d = 0; d < dim; ++d) {
        f[d] += cp.weight * cp.tracTang[d];
      }
    }
  }
  // The magnitude of the RESULTANT, not the integral of the magnitude -- see the header.
  Double n2 = 0.0;
  for (UInt d = 0; d < dim; ++d) {
    n2 += f[d] * f[d];
  }
  return std::sqrt(n2);
}

Double ContactInterface::GetIntegratedArea() const {
  Double a = 0.0;
  for (UInt i = 0; i < segments_.GetSize(); ++i) {
    for (UInt j = 0; j < segments_[i].points.GetSize(); ++j) {
      a += segments_[i].points[j].weight;
    }
  }
  return a;
}

Double ContactInterface::GetSecondaryArea() const {
  Double a = 0.0;
  for (UInt i = 0; i < segments_.GetSize(); ++i) {
    a += segments_[i].elemArea;
  }
  return a;
}

UInt ContactInterface::GetNumContactElements() const {

  // The elements that actually carry points, NOT elemList_->GetSize().
  //
  UInt n = 0;
  for (UInt i = 0; i < elemList_->GetSize(); ++i) {
    const ContactNcSurfElem* ce =
        dynamic_cast<const ContactNcSurfElem*>(elemList_->GetNcSurfElem(i));
    if (ce != nullptr && ce->pointIdx.GetSize() > 0) {
      ++n;
    }
  }
  return n;
}

UInt ContactInterface::GetNumCells() const {
  UInt n = 0;
  for (UInt i = 0; i < segments_.GetSize(); ++i) {
    n += segments_[i].numCells;
  }
  return n;
}

} /* namespace CoupledField */
