#ifndef OPTIMIZATION_DESIGN_FEATUREMAPPINGDESIGN_HH_
#define OPTIMIZATION_DESIGN_FEATUREMAPPINGDESIGN_HH_

#include "Optimization/Design/FeaturedDesign.hh"
#include "Utils/ToolsFull.hh"
#include <forward_list>

namespace CoupledField
{
class CoefFunctionOpt;
class DesignMaterial;

/** This is a rather traditional feature  implementation. ShapeMappingDesign is Fabian's original shape variant
 * SpaghettiDesign has the continuous spline variant to represent continuous fibers and is based in the python script 
 * for the core functionality. A reduced spaghetti noodle with one segment is same as a feature mapping pill.
 * 
 * For feature mapping we have pills or obrounds which are two parallel lines with half circles at the ends.
 * 
 * Use feature_mapping.py to visualize the geometry information from the .density.xml.
 */
class FeatureMappingDesign : public FeaturedDesign
{
  struct ItemIP; // forward declaration
  struct IntegrationPoint;

public:
  /** @method either SPAGHETTI or SPAGHETTI_PARAM_MAT or FEATURE_MAPPING or FEATURE_MAPPING_PARAM_MAT
   * @param anisotropic the param mat case where we orient by feature and do all in FeatureMappingDesign
   */
  FeatureMappingDesign(StdVector<RegionIdType>& regionIds, PtrParamNode pn, ErsatzMaterial::Method method);

  virtual ~FeatureMappingDesign() {};

  /** late inital MapDesignToDensity() */
  virtual void PostInit(int objectives, int constraints) override;

    /** read the design and do a conditional mapping */
  int ReadDesignFromExtern(const double* space_in, bool doNotSetAndWriteCurrent = false) override;

  void ToInfo(ErsatzMaterial* em) override;

  /** Called from DensityFile::ReadErsatzMaterial() with load ersatz material (-x)
   * See SetupDesign for expected ordering
   * @param set the set from the density.xml
   * @param lower_violation the maximal violation */
  void ReadDensityXml(PtrParamNode set, double& lower_violation, double& upper_violation) override;

  /** add information about transition zone */
  void AddToDensityHeader(PtrParamNode pn) override;

  /** in the anisotropic case DesignSpace::ApplyPhysicalDesign() calls ParamMat::GetMechTensor() which calls this.
   * Will be called indirectly by AnisotropicGradientHelper num_var_by_feature * num_features * num_element!
   * Uses aniso_current_variable structure to know which variable we have to derive to.
   * Needs to be thread safe! */
  bool GetAnisoMechTensor(const Elem* elem, MaterialTensor<double>& mt, DesignElement::Type direction, const CoefFunctionOpt* coef, DesignMaterial* dm) const;

  /** to handle nodal special results like "featureDistance". order needs to be >= 2.
   * Uses node_ip_result_map */
  double GetNodalSpecialResult(unsigned int nodeNumb, const ResultDescription& descr);

  /** Triggered by the 'hessexport' attribute: for each function with curvature information
   * (Function::CalcCurvature()) append the exact shape Hessian and the current function value to
   * <simname>.hessian.xml, accumulating one <iteration> block (design variables + Hessian) for all
   * iterations */
  void WriteHessExportFile() override;

  /** assemble the exact total shape Hessian of func in optimization variable space (n x n).
   * For cfs.evalhessian, used by a second-order python optimizer. The objective gradient of func
   * must have been evaluated before (the aggregation and feature terms need d_func/d_rho_e).
   * @param H_opt resized to n x n, n = opt_shape_param_.GetSize()
   * @return false if func has no curvature information */
  bool CalcShapeHessian(const Function* func, Matrix<double>& H_opt) override;

  /** To be used in GetAnisoMechTensor and set in AnisotropicGradientHelper().
   * In the anisotropic case we compute the gradient much different from the isotropic case. 
   * In the isotropic case all it is a product from feature to shape to element to integration point 
   * In the anisotropic case we go by variable and need to communicate FeatureMappingParamMat::SetElement() the current variable 
   * From the AnisotropicGradientHelper, s * N times GetAnisoMechTensor() is called */
  struct AnisoCurrentVariable
  {
    void Reset(); // set to non-initialized
    BaseDesignElement::Type var_type = BaseDesignElement::NO_DERIVATIVE; // set to FEATURE_MAPPING_PX, ... for shape derivatives in SetElementK()
    int var_idx = -1; // 0 ... 4 max number of a single feature variables
    int feature_idx = -1; // together with aniso_shape_derivative this specifies what we are currently deriving. Only used in the FeatureMappingParamMat case
    Matrix<double> RCdR;  // working array for d_C_0/d_s which is the only feature and variable dependent (without rho_e)
 
    StdVector<int> dJ_ds_res_idx;      // d_J/d_s:      <result value="featureGrad"      detail="compliance"   design="feature_var_Px" generic="1" ../>
    StdVector<int> dmrho_drho_res_idx; // d_mrho/d_rho: <result value="featureGrad"      detail="combine"      design="feature"        generic="1" ../>
    StdVector<int> dmrho_ds_res_idx;   // d_mrho/d_s:   <result value="featureGrad"      detail="combine"      design="feature_var_Px" generic="1" ../>
    StdVector<int> drho_ds_res_idx;    // d_rho/d_s:    <result value="featureGrad"      detail="featureRho"   design="feature_var_Px" generic="1" ../>
  };

  AnisoCurrentVariable aniso_current_variable; // set by AnisotropicGradientHelper() and used in FeatureMappingParamMat::SetElementK()

protected:
  /** A pill is a bar with P and Q end points and a circular cap and a width p. 
   * This are exactly the basic Feature variables. We add mapping to the base Feature */
  struct Pill : public Feature  
  {
  public:

    const char* GetName() override { return "pill"; } 

    /** Updates internal variables for efficient Distance and Grad calculation.
     * Call this after every variable change. */
    void Update() override;

    /** Assumes Update() has been called properly. 
     * As long as we don't need this is SpaghettiDesign or any thing other than Pill, not virtual 
     * @param part if given, the closes part is returned (relevant for gradient)*/
    double Distance(const Point& p, FeatureVariable::Tip* part = nullptr) const override;

    /** for efficiency reason calculate all gradients at once. Is a waste for fixed variables */
    void GradDistance(const Point& X, FeatureVariable::Tip part, Vector<double>& out) const override;

    /** second derivatives of the distance, symmetric. For p always zero.
     * @param out resized to num_var_by_feature^2 and fully set */
    void HessDistance(const Point& X, FeatureVariable::Tip part, Matrix<double>& out) const override;

    /** exact Hessian of the pill length ||Q-P|| (the native 'distance' constraint) w.r.t. the feature
     * variables [Px Py (Pz) Qx Qy (Qz) P]; the profile p does not enter. Purely geometric (no transition
     * function), so independent of the boundary order. @param out resized to num_var_by_feature^2 and fully set. */
    void HessLength(Matrix<double>& out) const;

    /** calculate the gradient of the angle in the anisotropic case, 0 for p
     * @param s index 0 ... 4 */
    inline double GradAngle(int s) const;

    double angle = -1.0; // only for anisotropic

    /** material tensor rotated by angle. Only anisotropic */ 
    Matrix<double> rot_mat;
  protected:
    // this are helper variables for Pill for efficient distance calculation. Call Update() on every change!
    // P and Q are in Feature
    Point U; // normalized vector Q-P
    Point V; // normal to U - 2D only, in 3D the perpendicular direction depends on the point
    double length = 0.0; // || start_ - end_ ||
    double length2 = 0.0; // length^2
    double dx = 0.0; // the x-component of the vector from start to end - 2D only
    double dy = 0.0; // 2D only
    double dp_norm = 0.0; // Q_x*P_y - Q_y*P_x - 2D only
  };

  /** base class for boundary function */
  struct BoundaryFunction
  {
    virtual ~BoundaryFunction() {}
    virtual double Eval(double dist) const = 0; // outside is > 0, inside < 0
    virtual double Grad(double dist) const = 0;
    /** second derivative for exact Hessian computation, requires a C^2 boundary function.
     * Implemented by Quintic and Bezier, the default throws as Poly is only C^1. */
    virtual double Hessian(double dist) const;
    inline void Eval(const Vector<double>& dist, Vector<double>& eval) const;
    bool IsConstant(double dist) const { return dist <= -h || dist >= h; }
    inline bool IsConstant(const ItemIP* item_ip) const;
    bool IsInside(double dist) const { return dist <= -h;}
    bool IsOutside(double dist) const { return dist >= h;}
    double Integrate(const Vector<double>& projected) const { return projected.Sum() / projected.GetSize(); }
    /** for debuging purpose. Write some standard cases to LOG */
    std::string DebugLog();

    double h = -1; // half transition zone!
    double rhomin = -1;
  };

  /** third-order polynomial. See (15) in the review paper
   * prj(x) = 3/4*(1-rho_min) * (phi/h) - phi^3/(3 h^3)) + .5*(1+rho_min) with phi=is phi(x) */
  struct Poly : public BoundaryFunction
  {
    Poly(double transition, double rhomin);

    inline double Eval(double dist) const override;
    inline double Grad(double dist) const override;
    private:
    double fac = -1;
    double bias = -1;
    double h3 = -1;    
  };

  /** symmetric quintic polynomial transition function, the C^2 counterpart of Poly.
   * With h = transition/2 and s = d/h: H(d) = 1/2 - (15s - 10s^3 + 3s^5)/16, hence
   * H'(-h) = H'(h) = H''(-h) = H''(h) = 0 and the second derivative is continuous everywhere. */
  struct Quintic : public BoundaryFunction
  {
    Quintic(double transition, double rhomin);

    inline double Eval(double dist) const override;
    inline double Grad(double dist) const override;
    double Hessian(double dist) const override;
  private:
    double fac = -1;      // (1-rhomin)/16
    double bias = -1;     // (1+rhomin)/2
    double grad_fac = -1; // 15*(1-rhomin)/(16*h)
    double hess_fac = -1; // 15*(1-rhomin)/(4*h*h)
  };

  /** Quintic Bezier based transition function with optional asymmetric transition zone.
    * The curve is parameterized by t in [0,1] with control points W_i. Used via sampling. 
    * See bezier test for paper reference. */
  struct Bezier : public BoundaryFunction
  {
    /** @param transition the full width of the symmetric transition zone as for Poly, a = transition/2
     * @param extension enlarges the outer transition zone b = transition/2 + extension, 0 is the symmetric case */
    Bezier(double transition, double rhomin, double extension = 0.0);

    /** look-up table based, equivalent to Eval(dist, false) */
    inline double Eval(double dist) const override { return Eval(dist, false); }
    inline double Grad(double dist) const override { return Grad(dist, false); }

    /** second derivative H''(dist), continuous as the quintic curve is C^2.
     * look-up table based, equivalent to Hessian(dist, false) */
    inline double Hessian(double dist) const override { return Hessian(dist, false); }

    /** @param exact true solves t(dist) by bisection instead of interpolating the look-up tables */
    double Eval(double dist, bool exact) const;
    double Grad(double dist, bool exact) const;
    double Hessian(double dist, bool exact) const;

    double a = -1;    // inner transition zone, Eval(dist <= -a) = 1
    double b = -1;    // outer transition zone, Eval(dist >= b) = rhomin
    double gamma = 0; // shape parameter for the x-components of the inner control points

  private:
    /** Bernstein sum of degree n: sum_i C[i] * t^i * (1-t)^(n-i) * W[i] with binomials C */
    static double Bernstein(int n, const double* C, const double* W, double t);

    /** quintic Bezier component value at curve parameter t */
    static double EvalByT(const double* W, double t);

    /** first parametric derivative d/dt, a quartic Bezier of the control point differences */
    static double FirstDerivT(const double* W, double t);

    /** second parametric derivative d^2/dt^2 */
    static double SecondDerivT(const double* W, double t);

    /** third parametric derivative d^3/dt^3 */
    static double ThirdDerivT(const double* W, double t);

    /** H'(d) = y'(t)/x'(t) */
    static double BezierDyDx(const double* Wx, const double* Wy, double t);

    /** H''(d) = (y'' x' - y' x'') / x'^3 */
    static double BezierD2yDx2(const double* Wx, const double* Wy, double t);

    /** H'''(d) = (x'(y''' x' - y' x''') - 3 x''(y'' x' - y' x'')) / x'^5, only as Hermite slope data for the Hessian table */
    static double BezierD3yDx3(const double* Wx, const double* Wy, double t);

    /** invert x(t) = dist by bisection, x is strictly monotone in [0,1] by FindGamma()
     * @param tol absolute precision of the curve parameter t. The default keeps the look-up
     * table node error (~tol in H) an order below the Hermite interpolation error */
    double SolveT(double dist, double tol = 1e-12) const;

    /** gamma* = argmin_gamma max_t |H''(t; gamma)| s.t. x'(t) > 0, grid search as in the paper */
    double FindGamma() const;

    /** cubic Hermite interpolation in a look-up table sampled on the uniform grid over [-a, b]
     * @param dtab the exactly sampled derivative of tab, used as Hermite slopes */
    inline double Interpolate(const Vector<double>& tab, const Vector<double>& dtab, double dist) const;

    double Wx[6]; // x-components of the control points, set after FindGamma()
    static constexpr double Wy[6] = {1.0, 1.0, 1.0, 0.0, 0.0, 0.0}; // fixed y-components

    /** number of look-up table intervals, set in the constructor as 500*(b/a), capped at 20000:
     * the curve is sharpest in the inner zone of scale a, so the resolution follows the
     * asymmetry ratio. See the constructor for the accuracy this provides */
    int n_tab = -1;
    Vector<double> H_tab;    // n_tab+1 samples of H(d) in [0,1] (without rhomin scaling)
    Vector<double> dH_tab;   // n_tab+1 samples of H'(d)
    Vector<double> ddH_tab;  // n_tab+1 samples of H''(d)
    Vector<double> dddH_tab; // n_tab+1 samples of H'''(d), only Hermite slope data for ddH_tab
    double delta = -1;     // grid spacing (a+b)/n_tab of the look-up tables
  };

  struct CombineFunction
  {
    virtual ~CombineFunction() {};
    /**@param projected list of projected distances (BoundaryFunction::Eval(dist)) for each feature */
    
    virtual double Eval(const Vector<double>& projected) const = 0;

    /** gradient with respect to feature number */
    virtual double Grad(const Vector<double>& projected, unsigned int num) const = 0;

    /** second derivative with respect to feature numbers f and g, symmetric in f and g */
    virtual double Hessian(const Vector<double>& projected, unsigned int f, unsigned int g) const = 0;
  };

  /** non-differentiable combine function, primarily useful for debugging */
  struct Max : public CombineFunction
  {
    inline double Eval(const Vector<double>& projected) const override;
    inline double Grad(const Vector<double>& projected, unsigned int num)  const override;
    inline double Hessian(const Vector<double>& projected, unsigned int f, unsigned int g) const override { return 0.0; }
  };

  /** combine rho_e by p-max aggregation */
  struct P_Norm : public CombineFunction
  {
    P_Norm(int p) { this->p = p; }
    double Eval(const Vector<double>& projected) const override;
    double Grad(const Vector<double>& projected, unsigned int num)  const override;
    double Hessian(const Vector<double>& projected, unsigned int f, unsigned int g) const override;
    int p; // norm exponent
  };

  /** combine rho_e by Kreiselmeier-Steinhauser aggregation */
  struct KreisselmeierSteinhauser : public CombineFunction
  {
    KreisselmeierSteinhauser(int beta) { this->beta = beta; }
    double Eval(const Vector<double>& projected) const override {
      return SmoothMax(projected, beta, false);  };
    double Grad(const Vector<double>& projected, unsigned int num)  const override {
      return DerivSmoothMax(projected, beta, num); }
    double Hessian(const Vector<double>& projected, unsigned int f, unsigned int g) const override;
    int beta;
  };

  /** combine rho_e by soft max aggregation */
  struct SoftMax : public CombineFunction
  {
    SoftMax(int beta) { this->beta = beta; }
    double Eval(const Vector<double>& projected) const override;
    double Grad(const Vector<double>& projected, unsigned int num)  const override;
    double Hessian(const Vector<double>& projected, unsigned int f, unsigned int g) const override;
    int beta;
  };

  std::unique_ptr<BoundaryFunction> bnd_fnc;
  std::unique_ptr<CombineFunction>  cmb_fnc;
  double cmb_param = -1; // p for P_Norm, beta for KreiselmeierSteinhauser and SoftMax
private:

  /** add our integration framework*/
  void SetupDesign(PtrParamNode pn) override;

  /** add our extension */
  void SetupMapping() override;

  /** for feature mapping we have fixed integration points to be shared over elements (Item) 
   * for efficient determination where we are outside/inside/need integration. After SetupMapping() */
  void SetupFixedIntegrationPoints();

  /** little helper */
  void SetupParsedFeatures(PtrParamNode base) override;

  /** Map structure to rho (DesignSpace::data). Sets DesignSpace::data.
   *  Shall be called by ReadDesignFromExtern(). */
  void MapFeatureToDensity() override;

  /** Takes the density gradients and sums it up on the shape variables.
   *  To be called within WriteGradientToExtern().
   *  @param f the function we add the stuff to the gradient. */
  void MapFeatureGradient(Function* f) override;

  /** helper for the isotropic case for MapFeatureGradient() */
  void IsotropicGradientHelper(const Function* func, Vector<double>& dJ_ds) const;

  /** computes d_J/d_s = sum_e <lmbd_e, dK/ds u_e> with d_K/d_s containing d_rho/d_s and d_phi/d_s
    We precompute stuff for FeatureMappingParamMat::SetElementK(), hence not thread safe! */
  void AnisotropicGradientHelper(Function* func, Vector<double>& dJ_ds);

  /** helper to calculate d_rho/d_s for a feature.  For all 5 variables at once.
   * Not for rho_max!!
   * ddist_ds and elem_ip are provided working data to save allocation, resized and filled by the function.
   * This function is thread save */
  inline void CalcGradRhoByFeature(const Item& item, const Feature& feature, Vector<double>& out, Vector<double>& ddist_ds, StdVector<IntegrationPoint*>& elem_ip) const;

  /** helper to compute the hessians. the Jacobian D = d_mrho/d_s for all elements in the full
   * variable space (features * num_var_by_feature, including fixed variables).
   * Uses the precomputed ItemIP::rho and ItemIP::drho_ds_full.
   * @param D resized to (elements x full variables) */
  void CalcShapeJacobian(Matrix<double>& D) const;

  /** assemble the three Hessian terms (objective D^T diag(curv) D, aggregation, feature) in full
   * feature variable space. Shared by WriteHessExportFile() and CalcShapeHessian().
   * @return false if func has no curvature information */
  bool CalcShapeHessianTerms(const Function* func, Matrix<double>& H_obj, Matrix<double>& H_agg, Matrix<double>& H_feat);

  /** index of each full feature space variable in opt_shape_param_, -1 for fixed variables.
   * Handles linking (map_to). @param opt_idx resized to pills * num_var_by_feature */
  void BuildOptIndexMap(StdVector<int>& opt_idx) const;

  /** helper for CalcShapeHessianTerms() which sets H_agg and H_feat, 
   * the two geometric Hessian terms in full variable space, see eqn. 'hessian_terms' in the tracking paper.
   * H_agg: sum_e w_e * d^2 mrho/(d_rho_f d_rho_g) * d_rho_f/d_s_i * d_rho_g/d_s_j (cross-feature)
   * H_feat: sum_e w_e * d_mrho/d_rho_f * d^2 rho_f/(d_s_i d_s_j) (within-feature blocks)
   * with w_e = d_func/d_mrho_e from the DesignElement gradient of func */
  void AssembleHessianTerms(const Function* func, Matrix<double>& H_agg, Matrix<double>& H_feat) const;

  /** helper for WriteHessExportFile(): build one term element <name> with its non-empty (f,g) blocks
   * as children, sliced from the full-space term matrix H. Term-major: the term is the parent, each
   * <block f g dim1 dim2> its child (data without the <real> wrapper). diagonal_only restricts to the
   * f==g blocks (the feature/geometry term). Returns a self-closing <name/> when H has no content
   * (e.g. the curvature term of a function linear in rho). @param nv num_var_by_feature, block size */
  std::string HessExportTermXML(const std::string& name, unsigned int nv,
      const Matrix<double>& H, bool diagonal_only) const;

  /** traverse all existing features (pills) and searches for first occurrence key.
   * consider to add simple formulas (-key, 1-key, key + 5, ...) when needed
   * @return nullptr if not found */
  FeatureVariable* FindKey(const std::string& key) const;

  /** count number of occurrences of key. > 1 will be an error */
  int CountKey(const std::string& key) const;

  /** An integration point is usually shared by multiple Item (fe mesh cell). 
   * We have corner points and higher order (order > 2) integration points which are dynamically added and removed  */ 
  struct IntegrationPoint
  {
    Point loc; // the coordinates of the integration point (make it pointer for automatically 2D or 3D) in meter space
    Vector<double> dist; // distance to the i-th feature im meter. Size of features_.
    StdVector<FeatureVariable::Tip> part; // to which part (start, end, inner) we have the shortest distance. Speeds up gradient

    /** debug helper. Prints loc */
    static std::string ToString(StdVector<IntegrationPoint*>& ips);
  };

  /** We extend FeaturedDesign::Item to hold our integration points as extension. */
  struct ItemIP : FeaturedDesign::ItemExtension
  {
    typedef enum { CORNER, INNER } Storage;
    StdVector<IntegrationPoint*> corner; // shared integration points at mesh cell corners - constant
    StdVector<IntegrationPoint*> inner;  // higher order integration, also shared at edges
    Point center; // barycenter for debug and some order stuff
    Vector<double> rho; // integrated rho for each feature to be combined. Necessary for gradient calculation
    Vector<double> drho_ds_full; // for any variable of any feature
        
    std::string ToString() override;
    
    /** Get all distances for a given feature num from corner and inner */
    void FillDist(Vector<double>& data, unsigned int num) const;

    /** helper which dynamically combines corner and inner */
    void GetAllIP(StdVector<IntegrationPoint*>& out);
    /** debug variant as it is too expensive (logging) */
    StdVector<IntegrationPoint*> GetAllIP();
  };

  /* helper for SetupDesign() to create an ip and add it to the extension of the given item
   @param ref the to be created integration point is ref + (dx,dy,dz) */
  IntegrationPoint* SetupDesignCreateAddIP(ItemIP::Storage storage, const Point& ref, ItemIP* item_ip, double dx, double dy, double dz = 0.0);

  /** helper for SetupDesign() to add ip at extension of item to corner or inner. */
  void SetupDesignAddIP(ItemIP::Storage storage, Item& item, IntegrationPoint* ip);

  /** when we have nodal results, we set up here the mapping from node to index in corners. -1 for no mapping */
  void SetupNodeIPMap(); 

  /** convenience helper which does the dynamic_cast */
  inline ItemIP* GetItemIP(unsigned int item_idx); 
  inline const ItemIP* GetItemIP(unsigned int item_idx) const; 

  /** helper for special results with a lot of sanity checks and exceptions in case */
  const Pill& GetFeatureForResult(const ResultDescription& rd) const;
  /* helper to read 'detail' 'feature_var_Px', ... 'feature_var_P' 
  @ return zero base < num_var_by_feature (5) */
  int GetVariableForResult(const ResultDescription& rd) const;

  /** for d_J/d_s, d_mrho/d_s or d_mrho/d_feature where mrho is 'combine'
   * 
   * when function is given, we have the case
   *  <result value="featureGrad" detail="compliance" design="feature_var_Px" generic="1" id="optResult_9" />
   * and return of size all variables (#feature * 5) the special result index for DesignElement or -1 if not matching.
   * By this we can have many results for the different features and variables.
   * 
   * This is d_mrho/d_f 
   * <result value="featureGrad" detail="combine" design="feature" generic="1" id="optResult_11" />
   * the return vector then has the size of number of features 
   * @param d_s false for d_mrho/d_f */
  StdVector<int> GetSpecialResultIndices(const Function* f, DesignElement::Detail detail = DesignElement::NONE, bool d_s = true) const;

  /** All integration points at element corners shared by up to 4/8 items. Constant for lifetime
   * Used to compute if we are within the transition zone and might have to do higher order integration.
   * There is no specific ordering guaranteed!! */
  StdVector<IntegrationPoint> corners; 

  /** Cache finding proper node for GetNodalSpecialResults(). -1 if there is no mapping to the 1-base nodeNumber.
   * Has size corners+1 with element 0 = -1 as there is no nodeNumber 0 and nodeNumber 1 belongs to index node 0.
   * The value is the 0-based index of corners which corresponds to the 1-based node number.
   * Requires order >= 2.
   * @see SetupNodeIPMap(); */
  StdVector<int> node_ip_result_map; 

  /** If integration order is 1 or > 2 there are the "inner" ip, 
   * however many of them are also shared on edges by 2 Item for order > 2.
   * We have a list as we don't know in advance the number of ip in the design and we must not Resize, otherwise
   * ItemIP::inner would dangle. 
   * For order > 2 we dynamically delete the data and recreate by MapFeatureToDensity() */ 
  std::forward_list<IntegrationPoint> inners;

  /** Add Pill : public Feature once needed */
  StdVector<Pill> pills;
  /** the number of variables by feature, assume constant for all features.
   * The order matches Feature::GetAllVariables(): 2D [Px Py Qx Qy p], 3D [Px Py Pz Qx Qy Qz p] */
  const unsigned int num_var_by_feature = 2 * dim_ + 1;

  /** for boundary functions linear and poly this is the full transition zone 2*h -> move to FeaturedDesign */
  double transition = -1;

  /** asymmetric outer transition zone extension for the bezier boundary (b = transition/2 + extension), 0 is symmetric */
  double extension = 0;

  /** assemble and write the exact shape Hessian every iteration, see WriteHessExportFile().
   * To actually trigger hessian optimization you need to specify an own Python optimizer */
  bool hessexport_ = false;

  /** persistent root of <simname>.hessian.xml. WriteHessExportFile() appends one <iteration> child
   * per call so the file accumulates all iterations; it is re-written each time. */
  PtrParamNode hessexport_root_;

  /** the rhomin we use, extracted from the first density variable. */
  double rhomin = -1;
  double rhomax = -1;

  /** do we do anisotropic feature mapping assuming the simulation material for the optimization regions to be 
   * transversal isotropic (orthotropic in 2D) with the horizontal axis the main axis */
  bool anisotropic = false;

  /** for the anisotropic case only: we assume a single anisotropic tensor for optimization to be rotated by feature */
  static Matrix<double> aniso_base_tensor;

  PtrParamNode fm_info_; // our own info
};


} // end of name space


#endif /* OPTIMIZATION_DESIGN_FEATUREMAPPINGDESIGN_HH_ */
