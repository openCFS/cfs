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

  virtual void PostInit(int objectives, int constraints);

  /** conditionally calls UpdateCoordinates()
   *  @see AuxDesign::ReadDesignFromExtern() */
  virtual int ReadDesignFromExtern(const double* space_in);

  /** overwrites DesignSpace::CompareDesign() */
  virtual bool CompareDesign(const double* space_in);

  /** writes design to the vector, beginning with shape variables (shape_param_) and then aux_design_ */
  virtual int WriteDesignToExtern(double* space_out, bool scaling = true) const;

  /** write gradient out to the vector, appending with shape gradient
   * Sparse and dense! */
  virtual void WriteGradientToExtern(StdVector<double>& out, DesignElement::ValueSpecifier vs, DesignElement::Access access, Function* f, bool scaling = true);

  /** same as in DesignSpace, setting elements to zero, but also aux elements */
  virtual void Reset(DesignElement::ValueSpecifier vs, DesignElement::Type design = DesignElement::DEFAULT);

  virtual void WriteBoundsToExtern(double* x_l, double* x_u) const;

  virtual unsigned int GetNumberOfVariables() const;

  /** In case DesignSpace::FindDesign() searches for SPLINE_BOX and CP.
   * @return either DesignSpace::FindDesign() or the index within param_ */
  virtual int FindDesign(DesignElement::Type dt, bool throw_exception = true) const;

  /** goes on opt_shape_param_ only! */
  virtual BaseDesignElement* GetDesignElement(unsigned int idx);

  virtual void ToInfo(ErsatzMaterial* em);

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
   *  @param Matrix<double> output of size num_nodes_inside x 3
   */
  void EvalAll(Matrix<double>& out);

  Vector<double> Eval(Vector<double> point);

  /**
   * @return StdVector<Vector<double>> derivative at point w.r.t all distortion values
   *                                   size = dim x distortion.GetSize()
   */
  StdVector<Vector<double>> EvalDerivative(Vector<double> point);

  /** Write splinebox to VTK file */
  void Write();

  /** Called from DensityFile::ReadErsatzMaterial() with load ersatz material (-x)
   * @param set the set from the density.xml
   * @param lower_violation the maximal violation */
  void ReadDensityXml(PtrParamNode set, double& lower_violation, double& upper_violation);

  /** This is the variant of Function::Local::SetupVirtualElementMap() for slope constraints on ShapeParamElements.
   * This function is called within Function::Local() constructor, therefore Function::GetLocal() cannot work yet!
   * @param locality the local type */
  void SetupVirtualShapeElementMap(Function* f, StdVector<Function::Local::Identifier>& virtual_element_map, Function::Local::Locality locality);

  /** For SHAPE_MAP design. Combines NODE and PROFILE. Simple implementation, does not handle symmetry */
  void SetupVirtualMultiShapeElementMap(Function* f, StdVector<Function::Local::Identifier>& virtual_element_map, Function::Local::Locality locality);

  /** Returns the matrix for linear injectivity constraints
   *  such that A*cp < 0 assures injectivity. The matrix
   *  applies to the control points and not the deformation. */
  Matrix<double> GetInjectivityMatrix();

  typedef enum { NONE, CUBIC } Interpolation;

  static Enum<Interpolation> interpolation;

private:

  /** Setup distortion variables from xml file */
  void SetupDesign(PtrParamNode pn);

  /** Extract optimization variables from distortion */
  void SetupOptParam();

  /** some sanity checks, e.g. volume shall not be linear */
  void CheckPlausibility();

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
  Matrix<double> GetBasisMatrix();

  /** Generate a matrix from the bspline basis, s.t. deformed_point = matrix * control_points
   *  Basis has to be set by calling GenerateBasis before.
   *  @return Matrix<double> of size 1 x total_num_cp */
  Matrix<double> GetBasisMatrix(Vector<double> point);

  /** Map (distorted) structure to rho (DesignSpace::data). Sets DesignSpace::data.
   *  Shall be called by ReadDesignFromExtern(). */
  void MapFeatureToDensity();

  /** Takes the density gradients and sums it up on the shape variables using map_.
   *  To be called within WriteGradientToExtern().
   *  @param f the function we add the stuff to the gradient. */
  void MapFeatureGradient(const Function* f);

  /** assumes rectangular grid */
  void ReadFeature(string file_in, string key = "last");

  void InterpolateFeature();

  void EvalAllCornerValues();

  double EvalAtCoord(Vector<double> point) const;

  Vector<double> EvalDerivativeAtCoord(Vector<double> point) const;

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


  /** reference to optimization as we need it in MapFeatureGradient() to get the functions */
  Optimization* opt_ = NULL; // set in PostInit() if we have optimization and not only external design for sim

  /** this is the design_id for the last MapFeatureToDensity() run */
  int mapped_design_ = -1;

  Matrix<double> injMatrix_;

  NumInt numInt_;
};

} //end of namespace

#endif /* OPTIMIZATION_DESIGN_SPLINEBOXDESIGN_HH_ */
