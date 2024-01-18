#ifndef OPTIMIZATION_DESIGN_FEATUREDDESIGN_HH_
#define OPTIMIZATION_DESIGN_FEATUREDDESIGN_HH_

#include "Optimization/Design/AuxDesign.hh"

namespace CoupledField {

class Optimization;

class FeaturedDesign : public AuxDesign
{
public:
  FeaturedDesign(StdVector<RegionIdType>& regionIds, PtrParamNode pn, ErsatzMaterial::Method method = ErsatzMaterial::NO_METHOD);

  virtual ~FeaturedDesign() { } ;

  virtual void PostInit(int objectives, int constraints) override;

  /** @see DesignSpace::ReadDesignFromExtern().
   * Already increments design_id in case!
   * This is the default implementation for part of SplineBox and almost all of Spaghetti.
   * It cannot apply for ShapeMap as there opt_shape_param_ is overwritten */
  virtual int ReadDesignFromExtern(const double* space_in, bool setAndWriteCurrent = true) override;

  /** overwrites DesignSpace::CompareDesign() */
  virtual bool CompareDesign(const double* space_in) override;

  /** writes design to the vector, beginning with shape_param_ and then aux_design_ */
  virtual int WriteDesignToExtern(double* space_out, bool scaling = true) const override;

  /** write gradient out to the vector, appending with shape gradient
   * Sparse and dense! */
  virtual void WriteGradientToExtern(StdVector<double>& out, DesignElement::ValueSpecifier vs, DesignElement::Access access, Function* f, bool scaling = true) override;

  /** write a line with all variables and all gradients to a single line. Called for each iteration */
  void WriteGradientFile() override;

  /** same as in DesignSpace, setting elements to zero, but also aux elements */
  virtual void Reset(DesignElement::ValueSpecifier vs, DesignElement::Type design = DesignElement::DEFAULT) override;

  virtual void WriteBoundsToExtern(double* x_l, double* x_u) const override;

  /** total of real variables, which is non-fixed feature variables plus slack variables from AuxDesign */
  unsigned int GetNumberOfVariables() const override {
    return aux_design_.GetSize() + opt_shape_param_.GetSize();
  }

  /** all feature variables including fixed ones and without slack variables */
  int GetNumberOfFeatureMappingVariables() const override { return shape_param_.GetSize(); }

  BaseDesignElement* GetDesignElement(unsigned int idx) override {
    return idx < opt_shape_param_.GetSize() ? opt_shape_param_[idx] : AuxDesign::GetDesignElement(idx);
  }

  /** does not go on opt_shape_param_ but on shape_param_ */
  BaseDesignElement* GetFeaturedDesignElement(unsigned int idx) { return shape_param_[idx];}

  /** handles spaghetti and spline box */
  virtual int FindDesign(DesignElement::Type dt, bool throw_exception) const override;

  virtual void ToInfo(ErsatzMaterial* em) override = 0 ;

  /** Called from DensityFile::ReadErsatzMaterial() with load ersatz material (-x)
   * @param set the set from the density.xml
   * @param lower_violation the maximal violation */
  virtual void ReadDensityXml(PtrParamNode set, double& lower_violation, double& upper_violation) = 0;

  /** This is the variant of Function::Local::SetupVirtualElementMap() for slope constraints on ShapeParamElements.
   * This function is called within Function::Local() constructor, therefore Function::GetLocal() cannot work yet!
   * @param locality just given to assert() it is PREV_AND_NEXT */
  virtual void SetupVirtualShapeElementMap(Function* f, StdVector<Function::Local::Identifier>& virtual_element_map, Function::Local::Locality locality) {
    throw Exception("FeatureMapping does not define SetupVirtualShapeElementMap for " + f->ToString());
  }

  /** For SHAPE_MAP design. Combines NODE and PROFILE. Simple implementation, does not handle symmetry */
  virtual  void SetupVirtualMultiShapeElementMap(Function* f, StdVector<Function::Local::Identifier>& virtual_element_map, Function::Local::Locality locality) {
    throw Exception("FeatureMapping does not define SetupVirtualMultiShapeElementMap for " + f->ToString());
  }

  /** Variant of SetupVirtualShapeElementMap() for the periodic constraint which is the first element of a shape minus the last */
  virtual void SetupCyclicVirtualShapeElementMap(Function* f, StdVector<Function::Local::Identifier>& virtual_element_map, Function::Local::Locality locality) {
    throw Exception("FeatureMapping does not define SetupCyclicVirtualShapeElementMap for " + f->ToString());
  }

  /** The integration strategy applied forItem::GetOrder() */
  typedef enum { CONSTANT_FULL, FULL_OR_NOTHING, TAILORED } IntStrategy;

  static Enum<IntStrategy> intStrategy_;

  /** combines our settings for numerical integration */
  struct NumInt
  {
    /** the constructor */
    void Init(FeaturedDesign* fd, PtrParamNode pn, PtrParamNode info);

    /** standard output, does not print warning here */
    void ToInfo(PtrParamNode info) const;

    IntStrategy strategy_;

    int max_order_;

    /** this is the decision value for CloseEnough() */
    double sensitivity_;

    /** for statistics: cells which need integration */
    unsigned int int_cells_cnt_ = 0;

    /** for statistics: sum of 1D integration on cells we integrate (-> avg. order)  */
    unsigned int int_cells_order_sum_ = 0;

  protected:
    /** number of finite elements in mesh */
    unsigned int cells_;

    /** the element spacing */
    double h = -1.0;
  };

protected:

  /** meant to set up the design variables. Possibly call SetupMapping() next() */
  virtual void SetupDesign(PtrParamNode pn) = 0;

  /** Map (distorted) structure to rho (DesignSpace::data). Sets DesignSpace::data.
   *  Shall be called by ReadDesignFromExtern(). */
  virtual void MapFeatureToDensity() = 0;

  /** Takes the density gradients and sums it up on the shape variables.
   *  To be called within WriteGradientToExtern().
   *  @param f the function we add the stuff to the gradient. */
  virtual void MapFeatureGradient(const Function* f) = 0;

  /** set n_, nx_, ny_, nz_. Shall do verification of the lexicographic ordering later.
   * @see SetupLexicographicMesh() */
  void SetupMeshStructure();

  /** set map_ */
  void SetupMapping();

  /** called in the optimization case in PostInit().
   * E.g. check meaningful constraints line no linear volume */
  virtual void CheckPlausibility();

  /** to conveniently handle the mapping shape param to design */
  struct Item
  {
    /** our Design Element, which is almost always the density (rho) but might also be an angle */
    DesignElement* elemval = nullptr;

    /** this is the index of the element within the rectangular lexicographic n_ domain.
     * usually often the position of Item within map_ but not for complex designs as in solar heater */
    int lexicographic_pos = -1;

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
    static double SetIPGetWeight(const FeaturedDesign* fd, StdVector<double>& ip, int ip_x, int ip_y, int ip_z, int max_order);

    /** maximal diff max_corner_value and min_corner_value for all shapes */
    double MaxDiffCornerValue() const;
  };

  /** MAX means that at each ip we consider only the shape which has the largest rho. Only the gradient of that shape will be considered at that ip.
   * The drawbacks are loss of material at overlaps and doubtful differentiability.
   *
   * TANH_SUM limits via tanh_l( sum(tanh(shape) ) where tanh_l maps to 0..1 with own beta and the beta within sum(tanh(shape)) is halfed.
   *
   * An issue is if the gradients shall be scaled down to match the factor by the cutting of max(sum,1) */
  typedef enum {NO_COMBINE, MAX, TANH_SUM, KS, P_NORM, SOFTMAX} Combine;

  /** this describes the continuation of a structure in 1D. See feature mapping review.
   * Not every class uses all boundary functions, this is handled in the schema file */
  typedef enum { NO_BOUNDARY, TANH, LINEAR, POLY } Boundary;

  /** gives the fiber orientation for anisotropic material in spaghettiParamMat
   * rounded: follows the spaghetti direction but at the endings uses orientations parallel to the boundary
   * straight: follows the spaghetti direction extended to the round endings */
  typedef enum {STRAIGHT, ROUNDED} Orientation;

  /** for a full virtual lexicographic mesh, the position from coordinates */
  unsigned int LexicographicPos(unsigned int x, unsigned int y, unsigned int z) const { return z * (ny_*nx_) + y * nx_ + x; }

  Boundary GetBoundary() const { return boundary_; }

  /** no need for static */
  Enum<Combine> combine;

  Enum<Boundary> boundary;

  Enum<Orientation> orientation;

  /** handles the overlapping of shapes, controls MapShapeToDensity() and has a very strong impact on MapShapeGradietn() */
  Combine combine_ = NO_COMBINE;

  Boundary boundary_ = NO_BOUNDARY;

  unsigned int order = 1;

  Orientation orientation_ = ROUNDED;

  /** shortcut to the dimension (2|3) */
  static unsigned int dim_;

  /** number of elements of elemval in x-direction. +1 for nodes! */
  unsigned int nx_ = 0;
  unsigned int ny_ = 0;
  unsigned int nz_ = 0; // 1 for 2D

  /** discretization spacing (assumed regular!) */
  double dx_ = 0;

  /** repeats nx_, ny_ and nz_ */
  Vector<unsigned int> n_;

  /** this is the design_id for the last MapFeatureToDensity() run */
  int mapped_design_ = -2; // make sure we do not match uninitialized design_id

  /** Measure MapFeatureToDensity() */
  shared_ptr<Timer> mapping_timer_;

  /** Measure MapFeatureGradient() */
  shared_ptr<Timer> gradient_timer_;

  /** This are the design parameters. This includes fixed elements which are not object
   * to optimization (scpip cannot handle lower bound == upper bound). See opt_shape_param_ */
  StdVector<ShapeParamElement*> shape_param_;

  /** This are the external design parameters which means shape_param_ w/o fixed and optional links.
   * Note that in ShapeMapDesign the variable is a vector of another type */
  StdVector<ShapeParamElement*> opt_shape_param_;

  /** mapping with size of elemval to ShapeParamElement pointers to design variables */
  StdVector<Item> map_;

  /** reference to optimization as we need it in MapFeatureGradient() to get the functions */
  Optimization* opt_ = NULL; // set in PostInit() if we have optimization and not only external design for sim

private:
  void SetEnums();
};

/** definition of inline static function has to be in header file */
inline double FeaturedDesign::Item::SetIPGetWeight(const FeaturedDesign* fd, StdVector<double>& ip, int ip_x, int ip_y, int ip_z, int order)
{
  assert(ip.GetSize() == dim_);

  double weight = 1.0;

  if(order == 1) {
    // midpoint integration
    ip[0] = 0.5;
    ip[1] = 0.5;
    if(dim_ == 3)
      ip[2] = 0.5;
  } else {
    assert(ip_x >= 0 && ip_x < order);
    assert(ip_y >= 0 && ip_y < order);
    assert((dim_ == 2 && ip_z == 0) || (dim_ == 3 && ip_z >= 0 && ip_z < order));

    ip[0] = (double) ip_x / (double) (order-1);
    ip[1] = (double) ip_y / (double) (order-1);
    if(dim_ == 3)
      ip[2] = (double) ip_z / (double) (order-1);

    assert(ip[0] >= 0 && ip[0] <= 1);
    assert(ip[1] >= 0 && ip[1] <= 1);
    assert(dim_ == 2 || (ip[2] >= 0 && ip[2] <= 1));

    weight *= 1./(order-1) * (ip_x == 0 || ip_x == order-1 ? 0.5 : 1.0);
    weight *= 1./(order-1) * (ip_y == 0 || ip_y == order-1 ? 0.5 : 1.0);
    if(dim_ == 3)
      weight *= 1./(order-1) * (ip_z == 0 || ip_z == order-1 ? 0.5 : 1.0);
  }

  return weight;
}


} //end of namespace

#endif /*OPTIMIZATION_DESIGN_FEATUREDDESIGN_HH_*/
