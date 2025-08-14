#ifndef OPTIMIZATION_DESIGN_FEATUREMAPPINGDESIGN_HH_
#define OPTIMIZATION_DESIGN_FEATUREMAPPINGDESIGN_HH_

#include "Optimization/Design/FeaturedDesign.hh"
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
   * Needs to be thread safe! */
  bool GetAnisoMechTensor(const Elem* elem, MaterialTensor<double>& mt, DesignElement::Type direction, const CoefFunctionOpt* coef, DesignMaterial* dm) const;

  /** in the anisotropic case we compute the gradient much different from the isotropic case. 
   * In the isotropic case all is a product from feature to shape to element to integration point 
   * In the anisotropic case we go by variable and need to communicate FeatureMappingParamMat::SetElement() the current variable */
  struct AnisoCurrentVariable
  {
    void Reset(); // set to non-initialized
    BaseDesignElement::Type var_type = BaseDesignElement::NO_DERIVATIVE; // set to FEATURE_MAPPING_PX, ... for shape derivatives in SetElementK()
    int var_idx = -1; // 0 ... 4 max number of a single feature variables
    int feature_idx = -1; // together with aniso_shape_derivative this specifies what we are currently deriving. Only used in the FeatureMappingParamMat case
    Matrix<double> RCdR;  // working array for d_C_0/d_s which is the only feature and variable dependent (without rho_e)
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
    Point V; // normal to U
    double length = 0.0; // || start_ - end_ ||
    double length2 = 0.0; // dx*dx + dy*dy w.o. sqrt
    double dx = 0.0; // the x-component of the vector from start to end
    double dy = 0.0; // 
    double dp_norm = 0.0; // Q_x*P_y - Q_y*P_x
  };

  /** base class for boundary function */
  struct BoundaryFunction
  {
    virtual ~BoundaryFunction() {}
    virtual double Eval(double dist) const = 0; // outside is > 0, inside < 0
    virtual double Grad(double dist) const = 0;
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

  struct CombineFunction
  {
    virtual ~CombineFunction() {};
    /**@param projected list of projected distances (BoundaryFunction::Eval(dist)) for each feature */
    
    virtual double Eval(const Vector<double>& projected) const = 0;
    
    /** gradient with respect to feature number */
    virtual double Grad(const Vector<double>& projected, unsigned int num) const = 0;
  };

  /** non-differentiable combine function, primarily useful for debugging */
  struct Max : public CombineFunction
  {
    inline double Eval(const Vector<double>& projected) const override;
    inline double Grad(const Vector<double>& projected, unsigned int num)  const override;
  };

  /** combine rho_e by p-max aggregation */
  struct P_Norm : public CombineFunction
  {
    P_Norm(int p) { this->p = p; }
    double Eval(const Vector<double>& projected) const override;
    double Grad(const Vector<double>& projected, unsigned int num)  const override;
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
    int beta; 
  };

  /** combine rho_e by soft max aggregation */
  struct SoftMax : public CombineFunction
  {
    SoftMax(int beta) { this->beta = beta; }
    double Eval(const Vector<double>& projected) const override;
    double Grad(const Vector<double>& projected, unsigned int num)  const override;
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
   @param ref the to be created integration point is ref + (dx,dy) */
  IntegrationPoint* SetupDesignCreateAddIP(ItemIP::Storage storage, const Point& ref, ItemIP* item_ip, double dx, double dy);

  /** helper for SetupDesign() to add ip at extension of item to corner or inner. */
  void SetupDesignAddIP(ItemIP::Storage storage, Item& item, IntegrationPoint* ip);

  /** convenience helper which does the dynamic_cast */
  inline ItemIP* GetItemIP(unsigned int item_idx); 
  inline const ItemIP* GetItemIP(unsigned int item_idx) const; 

  /** All integration points at element corners shared by up to 4/8 items. Constant for lifetime
   * Used to compute if we are within the transition zone and might have to do higher order integration.
   * There is no specific ordering guaranteed */
  StdVector<IntegrationPoint> corners; 

  /** If integration order is 1 or > 2 there are the "inner" ip, 
   * however many of them are also shared on edges by 2 Item for order > 2.
   * We have a list as we don't know in advance the number of ip in the design and we must not Resize, otherwise
   * ItemIP::inner would dangle. 
   * For order > 2 we dynamically delete the data and recreate by MapFeatureToDensity() */ 
  std::forward_list<IntegrationPoint> inners;

  /** Add Pill : public Feature once needed */
  StdVector<Pill> pills;
  /** the number of variables by feature, assume constant for all features */
  const unsigned int num_var_by_feature = 5;

  /** for boundary functions linear and poly this is the full transition zone 2*h -> move to FeaturedDesign */
  double transition = -1;

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
