#ifndef OPTIMIZATION_DESIGN_SPLINEBOXDESIGN_HH_
#define OPTIMIZATION_DESIGN_SPLINEBOXDESIGN_HH_

#include "Optimization/Design/FeaturedDesign.hh"
#include "Utils/BiCubicInterpolate.hh"
#include "Utils/TriCubicInterpolate.hh"
#include "Utils/BSpline.hh"

namespace CoupledField {

class Optimization;

class SplineBoxDesign : public FeaturedDesign
{
public:

  SplineBoxDesign(StdVector<RegionIdType>& regionIds, PtrParamNode pn, ErsatzMaterial::Method method = ErsatzMaterial::NO_METHOD);

  virtual ~SplineBoxDesign() { } ;

  void PostInit(int objectives, int constraints) override;

  /** conditionally calls UpdateCoordinates()
   *  @see AuxDesign::ReadDesignFromExtern() */
  int ReadDesignFromExtern(const double* space_in, bool setAndWriteCurrent = true) override;

  void ToInfo(ErsatzMaterial* em) override;

  /** Set or move one control point */
  void SetControlPoint(int idx, Point coords, bool add = false);

  /** Set or move all control points*/
  void SetControlPoints(StdVector<Point> coords, bool add = false);

  /** Returns the current control points of the BSpline curve */
  StdVector<Point> GetControlPoints() { return this->control_points_; }

  /** Returns the initial control points of the BSpline curve */
  StdVector<Point> GetInitialControlPoints() { return this->initial_control_points_; }

  /** Returns one initial control point of the BSpline curve */
  Point GetInitialControlPoint(unsigned int index) { return this->initial_control_points_[index]; }

  /** Evaluates basis
   *  @param Matrix<double> transformed nodes, size: num_nodes_inside x 3
   */
  void EvalAll(Matrix<double>& out);

  /** Returns the transformed coordinates of point */
  Vector<double> Eval(Vector<double> point);

  /**
   * @return StdVector<Vector<double>> derivative at point w.r.t all distortion values
   *                                   size = dim x distortion.GetSize()
   */
  StdVector<Vector<double>> EvalDerivative(Vector<double> point);

  /** Write splinebox to VTK file */
  void WriteBoxToVTK();

  /** Called from DensityFile::ReadErsatzMaterial() with load ersatz material (-x)
   * @param set the set from the density.xml
   * @param lower_violation the maximal violation */
  void ReadDensityXml(PtrParamNode set, double& lower_violation, double& upper_violation) override;

  /** This is the variant of Function::Local::SetupVirtualElementMap() for slope constraints on ShapeParamElements.
   * This function is called within Function::Local() constructor, therefore Function::GetLocal() cannot work yet!
   * @param locality the local type */
  void SetupVirtualShapeElementMap(Function* f, StdVector<Function::Local::Identifier>& virtual_element_map, Function::Local::Locality locality) override;

  /** For SHAPE_MAP design. Combines NODE and PROFILE. Simple implementation, does not handle symmetry */
  void SetupVirtualMultiShapeElementMap(Function* f, StdVector<Function::Local::Identifier>& virtual_element_map, Function::Local::Locality locality) override;

  /** Returns the matrix for linear injectivity constraints
   *  such that A*cp < 0 assures injectivity. The matrix
   *  applies to the control points and not the deformation. */
  Matrix<double> GetInjectivityMatrix();

  typedef enum { NONE, CUBIC } Interpolation;

  static Enum<Interpolation> interpolation;

private:

  /** Setup distortion variables from xml file */
  void SetupDesign(PtrParamNode pn) override;

  /** Extract optimization variables from distortion */
  void SetupOptParam();

  /** some sanity checks, e.g. volume shall not be linear */
  void CheckPlausibility() override;

  /** bounding_box is assumed to be rectangular */
  StdVector<Vector<Double>> GetPointsForBasis();

  void GenerateBasis(StdVector<Vector<Double>> points);

  /** Returns wether a point is inside the bounding box of the spline box */
  bool IsInside(Vector<double> point) const;

  /** Evaluate the splinebox with current control points and set nodes in FE mesh accordingly */
  void UpdateFEMesh(Matrix<double> new_coords);

  /** Generate a matrix from the bspline basis, s.t. deformed_mesh = matrix * control_points
   *  Basis has to be set by calling GenerateBasis before.
   *  @return Matrix<double> of size nodes inside splinebox x total_num_cp */
  Matrix<double> GetValuesOfBasisFunctions();

  /** Generate a matrix from the bspline basis, s.t. deformed_point = matrix * control_points
   *  Basis has to be set by calling GenerateBasis before.
   *  @return Matrix<double> of size 1 x total_num_cp */
  Matrix<double> GetValuesOfBasisFunctions(Vector<double> point);

  /** Map (distorted) structure to rho (DesignSpace::data). Sets DesignSpace::data.
   *  Shall be called by ReadDesignFromExtern(). */
  void MapFeatureToDensity() override;

  /** Takes the density gradients and sums it up on the shape variables using map.
   *  To be called within WriteGradientToExtern().
   *  @param f the function we add the stuff to the gradient. */
  void MapFeatureGradient(Function* f) override;

  /** assumes rectangular grid */
  void ReadFeature(string file_in, string key = "last");

  void InterpolateFeature();

  void EvalAllCornerValues();

  double GetDensityAtCoord(Vector<double> point) const;

  Vector<double> GetDensityDerivativeAtCoord(Vector<double> point) const;

  /** degree of splines */
  unsigned int degree_;

  /** number of control points in each dimension */
  Vector<unsigned int> num_cp_;

  /** number of total control points = prod(num_cp) */
  unsigned int total_num_cp_;

  /** vector of current control points */
  StdVector<Point> control_points_;

  /** vector of initial control points */
  StdVector<Point> initial_control_points_;

  /** bounding box of spline box */
  Matrix<double> bounding_box_;

  /** if outermost layer of splinebox control points is fixed */
  bool fixed_boundary_;

  /** if forward (FE mesh optimization) or backward (feature mapping) thinking */
  bool forward_;

  /** If design variable is subject to optimization. same size as param_ */
  StdVector<bool> is_opt_;

  StdVector<BSpline*> bsplines_;

  StdVector<Matrix<double>> basis_;

  StdVector<int> local_index_;

  StdVector<StdVector<double>*> unique_coords_;

  StdVector<StdVector<unsigned int>> unique_inverse_;


  /** Feature density field */

  /** resolution of density field */
  Vector<unsigned int> density_resolution_;

  /** values of density field */
  StdVector<double> density_;

  /** derivatives/finite differences of density field */
  StdVector<StdVector<double>> density_derivative_;

  /** if density field should be assumed to be periodic */
  bool periodic_;

  /** interpolation method for density field */
  Interpolation interpolation_;

  /** 2D interpolator for density field */
  BiCubicInterpolate* bicubicInterpolator_;

  /** 3D interpolator for density field */
  TriCubicInterpolate* tricubicInterpolator_;

  /** bounding box of spline box w.r.t. density field */
  Matrix<double> cover_box_;

  typedef enum { FILE, SUM_OF_SINE, MAX_SINE, SINE_X, SINE_Y, SINE_Z } AnalyticFunc;

  static Enum<AnalyticFunc> analyticFunc;

  AnalyticFunc analyticFunc_;

  /** factor for smooth max or min function */
  double beta_;

  /** factor to scale the feature */
  double feature_scale_;

  Matrix<double> injMatrix_;

  NumInt numInt_;
};

} //end of namespace

#endif /* OPTIMIZATION_DESIGN_SPLINEBOXDESIGN_HH_ */
