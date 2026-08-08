// ================================================================================================
/*!
 *       \file     ContactInterface.hh
 *       \brief    Geometry layer for genuine two-body mechanical contact.
 *
 *                 This class owns the *geometric* side of a contact pair: it pairs a
 *                 secondary ("slave") surface with a primary ("master") surface, generates
 *                 contact integration points on the secondary surface, projects each of them
 *                 onto the primary surface by closest-point projection, and evaluates the
 *                 normal gap.
 *
 *                 It deliberately does NOT touch the grid. Unlike MortarInterface, which
 *                 creates intersection nodes/elements via Grid::AddNode(), a contact interface
 *                 must be re-projected inside the Newton loop, where mutating the grid would
 *                 churn DOF numbering. All contact data therefore lives in this object.
 *
 *  ------------------------------------------------------------------------------------------
 *  FROZEN CONVENTIONS
 *  ------------------------------------------------------------------------------------------
 *
 *  1. NAMING       We reuse the MortarInterface vocabulary:
 *                    primary   == master  == the surface that is projected ONTO
 *                    secondary == slave   == the surface that carries the integration points
 *                  Contact tractions are integrated over the SECONDARY surface.
 *
 *  2. NORMAL       n is the outward unit normal of the PRIMARY (master) surface at the
 *                  projection point, i.e. it points out of the primary body into the space
 *                  where the secondary body sits.
 *
 *  3. GAP          g_N = (x_s - x_m) . n     with x = X + u  (current configuration)
 *                    g_N > 0  ->  separated  (no contact)
 *                    g_N = 0  ->  touching
 *                    g_N < 0  ->  penetration  (contact active)
 *                  This is the classical Signorini sign convention. Contact is active iff
 *                  g_N < 0 in the penalty setting.
 *
 *  4. TRACTION     t_N = eps_N * <g_N>_-  with <x>_- = min(x,0), so t_N <= 0 always.
 *                  The force on the secondary body is  -t_N * n  (i.e. pushing it away from
 *                  the primary body, along +n). Sign of the residual contribution must be
 *                  cross-checked against MechPDE::CoefFunction2ndPiolaTensor before use.
 *
 *  4b. TANGENTIAL  g_T = P_T (u_s - u_m)   with the tangential projector  P_T = I - n n^T,
 *                  i.e. the relative TANGENTIAL DISPLACEMENT of the frozen pairing.
 *
 *                  Unlike the normal gap there is NO reference offset: the closest-point
 *                  projection makes (x_s - x_m) purely normal at u = 0, so g_T is exactly
 *                  zero there and is exactly linear in u. That is why friction needs no
 *                  counterpart of ContactLinInt's g_0 term.
 *
 *                  t_T is the traction conjugate to g_T (virtual work INT t_T . dg_T), so
 *                  the friction force on the SECONDARY body is -t_T and on the primary +t_T
 *                  -- the same rowSign pattern the normal part uses. Coulomb bounds it by
 *                  ||t_T|| <= mu * p with p = -t_N >= 0 the contact pressure.
 *
 *  5. CONFIGURATION All projections and gaps are evaluated in the CURRENT (deformed)
 *                  configuration, x = X + u.
 *
 *                  X comes from the grid. u comes from the MECH_DISPLACEMENT FeFunction,
 *                  which must be handed in via SetDisplacementFunction() before the first
 *                  gap evaluation.
 *
                   Integration WEIGHTS are deliberately kept on the reference configuration,
 *                  consistent with the total-Lagrangian formulation MechPDE uses for
 *                  geometric nonlinearity.
 *
 */
//================================================================================================

#ifndef _CONTACTINTERFACE_HH_
#define _CONTACTINTERFACE_HH_

#include <map>
#include <set>
#include <utility>

#include "BaseNcInterface.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "Forms/IntScheme.hh"

namespace CoupledField {

// forward declarations
struct SurfElem;
class ElemShapeMap;
class BaseFeFunction;
template<class TYPE> class StdVector;
template<class TYPE> class Vector;
struct ContactPoint {

  // -------- secondary (slave) side --------
  //! Surface element on the secondary surface carrying this point
  SurfElem* secElem = nullptr;
  //! Local coordinate of the point on secElem
  LocPoint secLocal;
  //! Global coordinate of the point in the CURRENT configuration
  Vector<Double> secGlobal;
  //! Integration weight, already multiplied by the surface Jacobian determinant
  Double weight = 0.0;

  // -------- primary (master) side --------
  //! Surface element on the primary surface this point projects onto (nullptr if none)
  SurfElem* primElem = nullptr;
  //! Local coordinate of the projection point on primElem
  LocPoint primLocal;
  //! Global coordinate of the projection point in the CURRENT configuration
  Vector<Double> primGlobal;

  Vector<Double> normal;
  Double gap = 0.0;

  Double gap0 = 0.0;

  Double pressure = 0.0;

  //! Lagrange multiplier estimate for the normal traction, lambda <= 0 like t_N itself.
  //!
  //! Zero under FORM_PENALTY, and then every formula below degenerates to the pure penalty
  //! one -- which is exactly why augmented Lagrange reuses the penalty kernel unchanged.
  //!
  //! Under FORM_AUGMENTED_LAGRANGE the augmented traction is
  //!     t_N = < lambda + eps_N * g_N >_-
  //! and lambda is advanced ONLY between complete Newton solves, by UpdateMultipliers().
  //! Holding it fixed within a solve is what keeps the tangent the same positive definite
  //! penalty matrix and the Newton convergence quadratic.
  Double lambda = 0.0;

  //! Reference tangential slip, the tangential twin of gap0. The current
  //! slip decomposes as
  //!     g_T = slipOffset + P_T (u_s - u_m),      P_T = I - n n^T
  //! exactly as g_N = gap0 + (u_s - u_m).n does for the normal part, and for the same reason.
  //!
  //! NEVER CLEAR THIS WITHOUT CLEARING slipPlastic IN THE SAME PLACE. The two are measured in
  //! the same reference and only their difference, the elastic slip g_T - g_T^p, drives the
  //! traction. Zeroing this half alone, which RebaseTangentialSlip() used to do for a point
  //! that lost its projection, restarts g_T from zero on re-engagement while g_T^p keeps
  //! what it accumulated, so the elastic slip comes back as -g_T^p and the friction traction
  //! points ALONG the sliding direction.
  Vector<Double> slipOffset;

  Vector<Double> slipPlastic;

  Vector<Double> tracTang;

  //! True if a valid closest-point projection onto the primary surface was found.
  //! Points without a projection are inactive and contribute nothing.
  bool isProjected = false;

  bool slipAnchored = false;

  //! True if the point is currently in contact (gap < 0). Recomputed by UpdateGaps().
  bool isActive = false;

  //! True if the point is active AND inside the Coulomb cone, i.e. sticking. An active point
  //! is either sticking or slipping; an inactive one is neither.
  bool isSticking = false;

  //! As-meshed gap subtracted from every later gap evaluation under
  //! initialGapHandling="zero", so that an imperfectly meshed assembly starts out exactly
  //! touching instead of pre-stressed. Zero, and never set, under "fromMesh".
  //!
  //! It lives on the POINT rather than in a parallel per-point array on the interface
  //! because the point set is not fixed: under SLIDING_LARGE the mortar partition is rebuilt
  //! as the body slides (RebuildSegmentPoints()), which changes both the number of points
  //! and their order. An index-based array would silently pair each new point with some
  //! other point's offset.
  Double gapOffset = 0.0;

  //! True once gapOffset has been frozen for this point. A point created by a
  //! re-segmentation starts without one and gets it from its first projection, exactly as
  //! the original points got theirs from the first projection of all.
  bool gapOffsetSet = false;
};

struct ContactSegment {
  SurfElem* secElem = nullptr;
  StdVector<ContactPoint> points;
  Double elemArea = 0.0;
  Double elemSize = 0.0;
  Double penalty = 0.0;
  UInt numCells = 1;
};

struct ContactCell {
  //! 2 vertices (2D) or 3 vertices (3D), each of dimension dim-1
  StdVector<Vector<Double> > vtx;
};

struct ContactNcSurfElem : public MortarNcSurfElem {
  ContactSegment* segment = nullptr;
  StdVector<UInt> pointIdx;
};

// ================================================================================================
//  C L A S S   ContactInterface
// ================================================================================================
class ContactInterface : public BaseNcInterface {

  public:

    //! How the non-penetration constraint is enforced.
    enum Formulation {
      FORM_PENALTY,
      FORM_AUGMENTED_LAGRANGE
    };

    enum FrictionType {
      FRICTION_NONE,
      FRICTION_COULOMB
    };

    enum FrictionTangent {
      FRICTION_TANGENT_REDUCED,
      FRICTION_TANGENT_CONSISTENT
    };

    enum NormalSmoothing {
      NORMAL_FACET,
      NORMAL_NODAL
    };

    enum SlidingType {
      //! Pairing and normals frozen after the first projection; only the gap is updated.
      //! Valid while relative tangential motion stays well inside one primary element.
      SLIDING_SMALL,
      //! Full re-projection (pairing, normals, gap) every update.
      SLIDING_LARGE
    };

    ContactInterface(Grid* grid, PtrParamNode ctNode);

    virtual ~ContactInterface();

    // ======================================================================
    //  BaseNcInterface interface
    // ======================================================================

    bool IsMoving() const override { return false; }

    void ResetInterface() override;

    void UpdateInterface() override;

    // ======================================================================
    //  Contact-specific interface
    // ======================================================================
    void SetDisplacementFunction(shared_ptr<BaseFeFunction> dispFct) { dispFct_ = dispFct; }

    //! Generates the contact integration points on the secondary surface.
    //! Must be called once before UpdateProjection(). Idempotent.
    void BuildSegments();

    //! Recomputes, for every contact point, the closest-point projection onto the primary
    //! surface, the primary outward normal there, and the normal gap.
    //! For SLIDING_SMALL the pairing and normal are only computed on the first call; later
    //! calls just refresh the gap from the current geometry.
    void UpdateProjection();

    //! Re-cuts the mortar partition at a STEP boundary and carries the per-point history
    //! across the new point set
    //!
    void ReSegment();

    bool NeedsStablePointSet() const;

    void UpdateGaps();

    void SetReferenceModulus(Double e) { refModulus_ = e; }

    Double GetReferenceModulus() const { return refModulus_; }

    void UpdateTractions();

    Double EvalCurrentGap(const ContactPoint& cp) const;

    bool IsGapActive(const ContactSegment& seg, const ContactPoint& cp, Double gap) const {
      if (seg.penalty <= 0.0) {
        return gap <= activationTolerance_ * seg.elemSize;
      }

      const Double tol = activationTolerance_ * seg.elemSize * seg.penalty;
      return cp.lambda + seg.penalty * gap <= tol;
    }

    //! Artificial normal stiffness applied to projected-but-open points, to remove the
    //! singular direction of a body that is supported only through a contact interface and
    //! starts separated. Scaled exactly like the penalty, eps_stab = stabilization * E / h,
    //! so the dimensionless factor means the same thing across meshes and materials.
    //!
    //! Zero by default, i.e. off: it must be asked for, because it deliberately makes the
    //! tangent inconsistent (see ContactABInt.hh, MODE_STABILIZATION) and that costs the
    //! quadratic Newton convergence while any point is open.
    Double GetStabilizationStiffness(const ContactSegment& seg) const {
      return (seg.elemSize > 0.0) ? stabilization_ * refModulus_ / seg.elemSize : 0.0;
    }

    //! The dimensionless stabilisation factor as given by the user; 0 means off.
    Double GetStabilization() const { return stabilization_; }

    //! The augmented normal traction at a point, t_N = <lambda + eps_N * g_N>_- <= 0.
    //! Reduces to eps_N * <g_N>_- under FORM_PENALTY, where lambda == 0.
    Double EvalTraction(const ContactSegment& seg, const ContactPoint& cp, Double gap) const {
      return std::min(cp.lambda + seg.penalty * gap, 0.0);
    }

    // ======================================================================
    //  Friction
    // ======================================================================
    bool HasFriction() const {
      return frictionType_ == FRICTION_COULOMB && frictionMu_ > 0.0;
    }

    Double GetFrictionCoefficient() const { return frictionMu_; }

    FrictionTangent GetFrictionTangent() const { return frictionTangent_; }

    Double GetTangentialPenaltyFactor() const { return tangentialPenalty_; }

    Double GetTangentialPenalty(const ContactSegment& seg) const {
      return (seg.elemSize > 0.0) ? tangentialPenalty_ * refModulus_ / seg.elemSize : 0.0;
    }
  
    void EvalRelativeDisplacement(const ContactPoint& cp, Vector<Double>& du) const;

    void EvalCurrentSlip(const ContactPoint& cp, Vector<Double>& gT) const;

    //! Coulomb return map at one contact point.
    //!
    //!     s      = g_T - g_T^p                    elastic (trial) slip
    //!     t_tr   = eps_T s                        trial traction
    //!     Phi    = ||t_tr|| - mu p                Coulomb yield function
    //!
    //!     Phi <= 0   STICK   t_T = t_tr
    //!     Phi >  0   SLIP    t_T = mu p * t_tr/||t_tr||
    //!
    //! This is the exact tangential analogue of the normal Macaulay bracket, and of a radial
    //! return in plasticity: the stick branch is the elastic predictor, the slip branch the
    //! projection back onto the (pressure-dependent) Coulomb cone.
    bool EvalFrictionTraction(const ContactSegment& seg, const ContactPoint& cp,
                              Double gap, const Vector<Double>& gT,
                              Vector<Double>& tracT) const;

    void UpdateFrictionState();

    void UpdateSlipHistory();

    void ResetSlipHistory();

    UInt GetNumStickPoints() const;

    UInt GetNumSlipPoints() const;

    //! Largest tangential traction magnitude ||t_T|| over all points.
    Double GetMaxTangentialTraction() const;

    Double GetTotalTangentialForce() const;

    //! One Uzawa update of every multiplier: lambda <- <lambda + eps_N * g_N>_-
    //!
    //! Must be called only between complete Newton solves. Uses the currently cached gaps, so
    //! the caller has to have refreshed them (UpdateGaps()) against the converged solution.
    void UpdateMultipliers(Double& maxDelta, Double& maxLambda);

    void ResetMultipliers();

    Formulation GetFormulation() const { return formulation_; }
    UInt GetAugmentedMaxIter() const { return augmentedMaxIter_; }
    Double GetAugmentedTol() const { return augmentedTol_; }

    //! Grid this interface lives on; the integrators need it for shape maps.
    Grid* GetGrid() const { return ptGrid_; }

    //! Volume region sets, needed by the integrators to resolve a surface point onto the
    //! correct volume neighbour.
    const std::set<RegionIdType>& GetPrimaryVolRegions() const { return primaryVolRegionSet_; }
    const std::set<RegionIdType>& GetSecondaryVolRegions() const { return secondaryVolRegionSet_; }

    // ---------------- accessors ----------------

    RegionIdType GetPrimarySurfRegion() const { return primarySurfRegion_; }
    RegionIdType GetSecondarySurfRegion() const { return secondarySurfRegion_; }
    RegionIdType GetPrimaryVolRegion() const { return primaryVolRegion_; }
    RegionIdType GetSecondaryVolRegion() const { return secondaryVolRegion_; }

    const StdVector<ContactSegment>& GetSegments() const { return segments_; }
    StdVector<ContactSegment>& GetSegments() { return segments_; }

    SlidingType GetSlidingType() const { return slidingType_; }

    NormalSmoothing GetNormalSmoothing() const { return normalSmoothing_; }

    //! Dimensionless normal penalty factor as given by the user. The integrator is
    //! responsible for scaling it by E*h
    Double GetNormalPenalty() const { return normalPenalty_; }

    //! Whether the contact force is assembled into the equation system (default yes).
    //! "no" gives a geometry-only pair that reports gap and pressure but changes nothing.
    bool GetApplyForce() const { return applyForce_; }

    //! Number of contact points with a valid projection.
    UInt GetNumProjectedPoints() const;

    //! Number of contact points currently in contact (gap < 0).
    UInt GetNumActivePoints() const;

    //! Smallest (most negative) gap over all projected points. Returns 0 if none.
    Double GetMinGap() const;

    Double GetIntegratedArea() const;

    Double GetSecondaryArea() const;

    UInt GetNumContactElements() const;

    UInt GetNumCells() const;

    UInt GetNumUnsegmented() const { return numUnsegmented_; }

    //! Whether mortar segmentation is switched on.
    bool GetMortarSegmentation() const { return mortarSegmentation_; }

    //! Resultant normal contact force, integral of the pressure over the secondary surface.
    //!
    //! This is INT p dG, the integral of a magnitude -- it does not resolve direction and it
    //! is only the force the pair transmits when the normal is the same everywhere on the
    //! patch. On a curved interface use GetContactForceResultant() instead.
    Double GetTotalContactForce() const;

    //! Vector resultant of the whole contact traction on the SECONDARY body,
    //!   F = INT (p*n + t_T) dG,
    //! with n the primary outward normal at the projection point, i.e. the traction the pair
    //! actually applies. The friction part is included when the pair has friction, so this
    //! is the complete interface load, not just its normal share.
    //!
    //! This is the quantity global equilibrium constrains, and unlike GetTotalContactForce()
    //! it stays exact when the normal varies along the patch. For a secondary body loaded by
    //! an applied pressure p_a over a flat top of length L and held tangentially by a
    //! displacement boundary condition (whose reaction is purely in the constrained
    //! component), the transverse component of F must equal p_a*L to solver tolerance,
    //! whatever the curvature does to the pressure distribution or to the segmentation. That
    //! makes it the one contact number a curved-interface testcase can predict analytically.
    //!
    //! \param f (out) resized to the grid dimension
    void GetContactForceResultant(Vector<Double>& f) const;

    //! Largest contact pressure over all points.
    Double GetMaxPressure() const;

    //! Resultant of the reference-gap force, sum over active points of eps_N * |g_0| * w.
    //!
    Double GetReferenceGapForce() const;

  private:

    //! Resolves primary/secondary surface regions and their adjacent volume regions.
    void SetRegions(const PtrParamNode ctNode);

    //! Reads tolerances, penalty, sliding type, initial-gap handling.
    void SetOptions(const PtrParamNode ctNode);

    //! Builds the assembler-facing list of ContactNcSurfElem, one per segment that has at
    //! least one projected point. These elements are not added to the grid.
    void BuildElementList();

    //! Collects the SurfElem pointers of both surfaces into primarySurfElems_ /
    //! secondarySurfElems_ and precomputes primary bounding boxes for the candidate search.
    void CollectSurfaceElements();

    //! Recomputes the axis-aligned bounding box of every primary surface element in the
    //! CURRENT configuration. Cheap; called at the start of a full re-projection.
    void UpdatePrimaryBoundingBoxes();

    //! Cuts one secondary surface element into mortar sub-cells along the projected outlines
    //! of the primary elements that overlap it.
    void BuildSegmentCells(SurfElem* secElem, StdVector<ContactCell>& cells);

    void GenerateSegmentPoints(ContactSegment& seg);

    void RebuildSegmentPoints();

    struct PointHistoryFit {
      Vector<Double> lambda;
      Vector<Double> gapOffset;
      //! one coefficient vector per component of slipPlastic and of slipOffset
      StdVector<Vector<Double> > slipPlastic;
      StdVector<Vector<Double> > slipOffset;
      bool valid = false;
      bool gapOffsetSet = false;
    };

    //! Fits the history fields of the CURRENT points of every segment, before a re-cut.
    void SnapshotPointHistory(StdVector<PointHistoryFit>& fits) const;

    //! Evaluates the fits at the NEW points of every segment, after a re-cut.
    void RestorePointHistory(const StdVector<PointHistoryFit>& fits);

    //! Weighted least-squares fit of a linear polynomial in the secondary reference
    //! coordinates: minimise SUM_k w_k (c.b(xi_k) - f_k)^2 over c, with b = (1, xi, eta).
    //!
    //! Why weighted least squares on a linear basis and not an L2 projection onto the
    //! element's own shape functions, which is what the plan proposed. Because the data
    //! cannot support it: an uncut LINE3 segment carries TWO integration points against
    //! THREE nodal shape functions, so the projection's mass matrix is singular exactly in
    //! the most common case. A linear fit is the largest space the point set is guaranteed
    //! to determine (sdim+1 coefficients against at least sdim+1 points per sub-cell).
    //!
    void FitWeightedLinear(const StdVector<ContactPoint>& pts,
                           const StdVector<Double>& value,
                           Vector<Double>& coeff) const;

    Double EvalFit(const Vector<Double>& coeff, const LocPoint& lp) const;

    //! Step 1 of keeping the tangential slip continuous across a re-pairing: replaces
    //! slipOffset by the CURRENT total slip g_T, evaluated against the pairing that is about
    //! to be discarded. Called at the top of UpdateProjection(), before anything moves.
    //!
    //! Storing the anchor in slipOffset itself, rather than in a parallel array, is what lets
    //! a step-boundary re-cut carry it: it is then just another per-point field for
    //! SnapshotPointHistory() to fit and RestorePointHistory() to evaluate, on exactly the
    //! same footing as slipPlastic. A parallel array would be indexed by points that the
    //! re-cut has already replaced.
    void AnchorTangentialSlip();

    //! Step 2: turns the anchor back into an offset against the NEW pairing, and rotates both
    //! it and slipPlastic into the new tangent plane,
    //!     slipOffset  <- P_T^new ( anchor - (u_s - u_m)^new )
    //!     slipPlastic <- P_T^new slipPlastic
    //! Called from UpdateProjection() as soon as the new pairing and normal exist, and BEFORE
    //! UpdateGaps() -- which ends in UpdateTractions() and so reads the slip. Between
    //! AnchorTangentialSlip() and this call slipOffset holds a total, not an offset.
    //!
    //! Both are exact no-ops when the pairing did not move: the anchor is then P_T du, so the
    //! first line gives P_T du - P_T du = 0, and slipPlastic already lies in the plane it is
    //! being projected onto
    void RebaseTangentialSlip();

    //! Records, once and on the REFERENCE configuration, every primary element each secondary
    //! element could plausibly slide onto, so that BuildElementList() can emit an assembler
    //! element for all of them and the matrix sparsity pattern covers the pairing wherever it
    //! ends up.
    void ReserveGraphPrimaries();

    //! Primary elements whose inflated bounding box overlaps that of a secondary element.
    //! Element-level counterpart of FindCandidates(), used by the segmentation.
    void FindElementCandidates(SurfElem* secElem, StdVector<UInt>& candidates) const;

    //! Clips the convex polygon `subject` against the convex polygon `clip` (both counter-
    //! clockwise, in the secondary reference domain) by the Sutherland-Hodgman algorithm.
    //! Both polygons are 2D; used for the 3D case only.
    static void ClipPolygon(const StdVector<Vector<Double> >& subject,
                            const StdVector<Vector<Double> >& clip,
                            StdVector<Vector<Double> >& result);

    //! Signed area of a 2D polygon, positive for counter-clockwise ordering.
    static Double PolygonArea(const StdVector<Vector<Double> >& poly);

    void GetCornerCoords(SurfElem* elem, RegionIdType volRegion,
                         StdVector<Vector<Double> >& coords) const;

    static void GetReferencePolygon(Elem::ShapeType shape,
                                    StdVector<Vector<Double> >& poly);

    void FindCandidates(const Vector<Double>& point, StdVector<UInt>& candidates) const;

    bool ProjectOntoElement(const Vector<Double>& point,
                            SurfElem* elem,
                            LocPoint& localOut,
                            Vector<Double>& globalOut,
                            Vector<Double>& normalOut,
                            Double& distOut) const;

    bool ClosestPointLocal(const Vector<Double>& point,
                           SurfElem* elem,
                           RegionIdType volRegion,
                           Double clampSlack,
                           LocPoint& localOut,
                           bool& clamped) const;

    //! Shape map of a surface element on the configuration the PROJECTION runs in.
    //!
    //! SLIDING_SMALL -> the reference configuration. Exact there, and not an approximation:
    //!   the pairing is frozen at u = 0 by definition, so the reference geometry IS the
    //!   configuration it was frozen in.
    //!
    //! SLIDING_LARGE -> the CURRENT configuration, built from the deformed nodal coordinates
    //!   X_I + u_I. Because the geometry and the displacement of a Lagrange element share
    //!   their shape functions, a shape map on those coordinates reproduces the deformed
    //!   surface exactly -- its Local2Global is x_m(xi), its CalcJ the deformed tangents, and
    //!   its CalcNormal the deformed normal. That last point is what the closest-point
    //!   projection actually needs: the stationarity condition (x_s - x_m).dx_m/dxi = 0
    //!   involves the DEFORMED tangent, so projecting with reference tangents would converge
    //!   to the wrong point rather than merely converge slowly.
    //!
    //! Deliberately does NOT go through Grid::SetNodeOffset()/updated geometry, which is the
    //! other way to get a deformed shape map
    shared_ptr<ElemShapeMap> GetProjectionShapeMap(SurfElem* elem,
                                                   RegionIdType volRegion) const;

    //! Deformed nodal coordinates (sdim x numNodes) of a surface element, X_I + u_I.
    void GetDeformedElemCoords(SurfElem* elem, RegionIdType volRegion,
                               Matrix<Double>& coords) const;

    void EvalOutwardNormal(SurfElem* elem, const LocPoint& lp,
                           Vector<Double>& normal) const;

    void EvalFacetNormalRaw(SurfElem* elem, const LocPoint& lp,
                            Vector<Double>& normal) const;

    Double GetFacetOrientation(SurfElem* elem) const;

    void BuildNodalNormals() const;

    //! Shape map of a primary element whose "coordinates" are its nodal normals, so that
    //! Local2Global() evaluates SUM_a N_a(xi) n_a -- the interpolation of the averaged
    //! normals, using exactly the element's own shape functions. Cached like
    //! projEsmCache_ and invalidated with it.
    shared_ptr<ElemShapeMap> GetNormalShapeMap(SurfElem* elem) const;

    void ClearGeometryCaches() const;

    //! Maps a local point to global coordinates in the CURRENT configuration, i.e.
    //! x = X(lp) + u(lp)
    void EvalCurrentPosition(SurfElem* elem, const LocPoint& lp,
                             RegionIdType volRegion, Vector<Double>& global) const;

    //! Maps a local point to global coordinates in the REFERENCE configuration.
    void EvalReferencePosition(SurfElem* elem, const LocPoint& lp,
                               Vector<Double>& global) const;

    //! Samples MECH_DISPLACEMENT at a point on a surface element. Returns zero when no
    //! displacement function has been supplied.
    void EvalDisplacement(SurfElem* elem, const LocPoint& lp, RegionIdType volRegion,
                          Vector<Double>& u) const;

    bool IsInsideElement(Elem::ShapeType shape, const Vector<Double>& local,
                         Double tol) const;

    void ClampToElement(Elem::ShapeType shape, Vector<Double>& local, Double slack) const;

    // =======================================================================
    //  Class variables
    // =======================================================================

    //! surface region ids of both sides
    RegionIdType primarySurfRegion_;
    RegionIdType secondarySurfRegion_;
    //! adjacent volume region ids of both sides
    RegionIdType primaryVolRegion_;
    RegionIdType secondaryVolRegion_;

    //! surface elements of both sides
    StdVector<SurfElem*> primarySurfElems_;
    StdVector<SurfElem*> secondarySurfElems_;

    //! axis-aligned bounding boxes of the primary elements, current configuration.
    //! Layout: [i] = (min_0..min_{d-1}, max_0..max_{d-1})
    StdVector<Vector<Double> > primaryBoxMin_;
    StdVector<Vector<Double> > primaryBoxMax_;

    //! the contact points, grouped per secondary surface element
    StdVector<ContactSegment> segments_;

    //! MECH_DISPLACEMENT field, used to reach the current configuration. May be null.
    shared_ptr<BaseFeFunction> dispFct_;

    //! single-element region sets, cached for LocPointMapped::SetWithSurface()
    std::set<RegionIdType> primaryVolRegionSet_;
    std::set<RegionIdType> secondaryVolRegionSet_;

    //! integration scheme used to place the contact points
    IntScheme intScheme_;
    //! integration order on the secondary surface
    UInt integrationOrder_;

    //! sliding kinematics
    SlidingType slidingType_;

    //! where the primary contact normal comes from; facet normal by default
    NormalSmoothing normalSmoothing_;

    //! constraint enforcement, penalty or augmented Lagrange
    Formulation formulation_;
    //! maximum Uzawa outer iterations; the solve is accepted with a warning if it runs out
    UInt augmentedMaxIter_;
    //! relative convergence tolerance of the Uzawa loop, on |dLambda| against |lambda|
    Double augmentedTol_;

    //! true  -> mortar segmentation: the secondary element is cut along the projected primary
    //!          element outlines and each sub-cell is integrated separately
    //! false -> Gauss-point-to-segment: the secondary element's own quadrature is used
    //! Kept switchable because it is the one knob that changes the discretization, so a
    //! result can always be compared against the simpler scheme.
    bool mortarSegmentation_;

    //! Number of secondary elements for which the segmentation could not be trusted and the
    //! plain element quadrature was used instead. Reported, because it says how much of the
    //! interface is still integrated the GPTS way.
    UInt numUnsegmented_;

    //! dimensionless normal penalty factor, scaled by E/h in UpdateTractions()
    Double normalPenalty_;

    //! tangential contact law; FRICTION_NONE unless a <friction> element is present
    FrictionType frictionType_;
    //! Coulomb friction coefficient
    Double frictionMu_;
    //! which linearisation of the slip branch enters the tangent; reduced by default
    FrictionTangent frictionTangent_;
    //! dimensionless tangential (stick) penalty, scaled by E/h in GetTangentialPenalty()
    Double tangentialPenalty_;

    //! dimensionless stabilisation factor for a separated start, scaled by E/h.
    //! 0 = off, which is the default
    Double stabilization_;

    //! assemble the contact force into the system? default true
    bool applyForce_;

    //! representative stiffness used to give the penalty its dimension, see
    //! SetReferenceModulus(). Zero until the PDE supplies it.
    Double refModulus_;

    //! true  -> the gap measured at build time is subtracted, so the as-meshed state is
    //!          treated as exactly touching (removes overclosure from imperfect meshing)
    //! false -> the measured gap is used as-is
    bool zeroInitialGap_;

    //! relative tolerance of the active-set test, scaled by the segment size. See
    //! IsGapActive(); it exists only to keep the exactly-touching case off the rounding
    //! lottery, not to widen contact.
    Double activationTolerance_;

    //! Newton tolerance for the closest-point projection (on the local coordinate)
    Double projectionTolerance_;
    //! maximum Newton iterations for the closest-point projection
    UInt projectionMaxIter_;
    //! relative inflation of the primary bounding boxes for the candidate search
    Double boundingBoxTolerance_;
    //! projections farther away than this (relative to the element size) are rejected
    Double searchDistanceFactor_;

    //! Deformed shape maps of primary surface elements, keyed by element number. Valid only
    //! for the displacement state of the UpdateProjection() pass that filled it, which is
    //! why that method clears it on entry. Empty under SLIDING_SMALL.
    mutable std::map<UInt, shared_ptr<ElemShapeMap> > projEsmCache_;

    //! Averaged unit normal per PRIMARY surface node, keyed by global node number. Built by
    //! BuildNodalNormals() in the same configuration as projEsmCache_ and invalidated with
    //! it. Empty unless normalSmoothing_ == NORMAL_NODAL.
    mutable std::map<UInt, Vector<Double> > nodalNormals_;
    //! true while nodalNormals_ describes the current configuration
    mutable bool nodalNormalsBuilt_ = false;

    //! Shape maps carrying the nodal normals as coordinates, see GetNormalShapeMap().
    //! Same lifetime as nodalNormals_.
    mutable std::map<UInt, shared_ptr<ElemShapeMap> > normalEsmCache_;

    mutable std::map<UInt, Double> facetOrientation_;

    StdVector<StdVector<SurfElem*> > graphPrimaries_;

    std::set<std::pair<UInt, UInt> > elemPairs_;
    //! true once elemPairs_ describes the graph rather than being collected into
    bool elemPairsFixed_ = false;

    bool partitionFrozenInSolve_ = false;

    //! Set only for the duration of one ReSegment() call, to let that one UpdateProjection()
    //! re-cut a partition that is otherwise frozen.
    bool resegmentNow_ = false;

    //! true once BuildSegments() has run
    bool segmentsBuilt_;
    //! true once a full projection has been performed at least once
    bool projected_;
};

} /* namespace CoupledField */

#endif /* _CONTACTINTERFACE_HH_ */
