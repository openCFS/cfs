#ifndef OPTIMIZATION_DESIGN_SPLINEBOXDESIGN_HH_
#define OPTIMIZATION_DESIGN_SPLINEBOXDESIGN_HH_

#include "Optimization/Design/AuxDesign.hh"
#include "Optimization/Design/SplineBoxDesign.hh"
#include "Utils/BiCubicInterpolate.hh"
#include "Utils/TriCubicInterpolate.hh"
#include "Utils/BSpline.hh"

namespace CoupledField {

class Optimization;

class SplineBoxDesign : public AuxDesign
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

  virtual int GetNumberOfFeatureMappingVariables() const { return distortion_.GetSize(); }

  /** goes on opt_distortion only! */
  virtual BaseDesignElement* GetDesignElement(unsigned int idx);

  /** does not go on opt_distortion but on distortion */
  ShapeParamElement* GetSplineBoxDesignElement(unsigned int idx) { return &(distortion_[idx]);}

  virtual void ToInfo(ErsatzMaterial* em);

  /** Set or move one control point */
  void SetControlPoint(int idx, Point coords, bool add = false);

  /** Set or move all control points*/
  void SetControlPoints(StdVector<Point> coords, bool add = false);

  /** Returns the control points of the BSpline curve */
  StdVector<Point> GetControlPoints() { return this->control_points_; }

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

  /** The integration strategy applied forItem::GetOrder() */
  typedef enum { CONSTANT_FULL, FULL_OR_NOTHING } IntStrategy;

  static Enum<IntStrategy> intStrategy;

  /** combines our settings for numerical integration */
  struct NumInt
  {
    /** the constructor with the ShapeMap pn. Calls SetTailored conditionally */
    void Init(SplineBoxDesign* sb, PtrParamNode pn, PtrParamNode info);

    /** standard output, does not print warning here */
    void ToInfo(PtrParamNode info) const;

    /** for Item::GetOrder(), however it makes only sense in debug to reduce it. */
    int max_order;

    IntStrategy strategy;

    /** this is the decision value for CloseEnough() */
    double sensitivity;

    /** for statistics: cells which need integration */
    unsigned int int_cells_cnt = 0;

    /** for statistics: sum of 1D integration on cells we integrate (-> avg. order)  */
    unsigned int int_cells_order_sum = 0;

  private:
    double beta = -1.0;

    /** number of elements in FE mesh */
    unsigned int cells_;

    /** what GetTailoredOrder() returns in the LINEAR case with 1st order integration */
    int linear_int_order_ = -1;
  };

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
  void ReadFeature(string file_in);

  void InterpolateFeature();

  void EvalAllCornerValues();

  double EvalAtCoord(Vector<double> point) const;

  Vector<double> EvalDerivativeAtCoord(Vector<double> point) const;

  /** Multiple subscripts from linear index
   *  This function corresponds to ShapeMapDesign::DensityIdx(int x, int y, int z)*/
  void Sub2Ind(Vector<unsigned int> size, StdVector<int> sub, unsigned int &ind) const;

  /** Linear index from multiple subscripts
   *  This function corresponds to ShapeMapDesign::DensityIdx(unsigned int i, Vector<unsigned int>& idx)*/
  void Ind2Sub(Vector<unsigned int> size, unsigned int ind, StdVector<int> &sub) const;

  /** dimension */
  static unsigned int dim_;

  /** degree of splines */
  unsigned int degree_;

  /** number of control points in each dimension */
  Vector<unsigned int> num_cp_;

  /** number of total control points = prod(num_cp) */
  unsigned int total_num_cp_;

  /** vector of current control points */
  StdVector<Point> control_points_;

  /** vector of initial control points */
  StdVector<Point> init_control_points_;

  /** bounding box of spline box */
  Matrix<double> bounding_box_;

  /** if outermost layer of splinebox control points is fixed */
  bool fixed_boundary_;

  /** if forward (FE mesh optimization) or backward (feature mapping) thinking */
  bool forward_;

  /** This are the design variables, which are offset to initial control_points
   *  (i.e. distortion/deformation). Contains offset for all control_points*/
  StdVector<ShapeParamElement> distortion_;

  /** If distortion variable is subject to optimization. same size as distortion_ */
  StdVector<bool> is_opt_;

  /** This are the external distortion variables which means shape_param
   *  w/o fixed and optional symmetry */
  StdVector<ShapeParamElement*> opt_distortion_;

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

  /** interpolation method for density field */
  Interpolation interpolation_;

  /** factor for smooth max function after interpolation */
  double beta_;

  /** 2D interpolator for density field */
  BiCubicInterpolate* bicubicInterpolator_;

  /** 3D interpolator for density field */
  TriCubicInterpolate* tricubicInterpolator_;

  typedef enum { FILE, SUM_OF_SINE, MAX_SINE } AnalyticFunc;

  static Enum<AnalyticFunc> analyticFunc;

  AnalyticFunc analyticFunc_;

  double analyticPeriod_;


  /** reference to optimization as we need it in MapFeatureGradient() to get the functions */
  Optimization* opt_ = NULL; // set in PostInit() if we have optimization and not only external design for sim

  /** this is the design_id for the last MapFeatureToDensity() run */
  int mapped_design_ = -1;

  /** Measure MapFeatureToDensity() */
  shared_ptr<Timer> mapping_timer_;

  /** Measure MapFeatureGradient() */
  shared_ptr<Timer> gradient_timer_;

  /** to conveniently handle the mapping shape param to design */
  struct Item
  {
    /** our Design Element */
    DesignElement* rho;

    /** This is the minimal corner value of all corners. Set by EvalAllCorners() */
    Vector<double> min_corner_value;
    Vector<double> max_corner_value;

    /** Determines based on corner_vals the order of integration for each shape. 0=void, 1=solid, >=2 integration.
     * @param order output
     * @param max_order the maximal order. This may conflict with sensitivity
     * @return the maximal number in order. This is <= the parameter max_order */
    int GetOrder(Vector<int>& order, const NumInt& numInt) const;

    /** Set generalized integration points and weights.
     *  Assumes rectangular finite elements!
     * @param ip result: vector [0...1]^dim
     * @param ip_x, ip_y, ip_y index of integration < max_order
     * @return closed Newton-Cotes weight for the 2D/3D point */
    static double SetIPGetWeight(const SplineBoxDesign* sbd, StdVector<double>& ip, int ip_x, int ip_y, int ip_z, int max_order);

    /** maximal diff max_corner_value and min_corner_value for all shapes */
    double MaxDiffCornerValue() const;
  };

  /** mapping with size of rho from ShapeParamElement pointers to distortion */
  StdVector<Item> map_;

  unsigned int int_order_;

  NumInt numInt_;

};

} //end of namespace

#endif /* OPTIMIZATION_DESIGN_SPLINEBOXDESIGN_HH_ */
