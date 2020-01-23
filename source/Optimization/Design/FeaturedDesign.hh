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

  virtual void PostInit(int objectives, int constraints) { AuxDesign::PostInit(objectives, constraints); };

  /** @see DesignSpace::ReadDesignFromExtern() */
  virtual int ReadDesignFromExtern(const double* space_in) = 0;

  /** overwrites DesignSpace::CompareDesign() */
  virtual bool CompareDesign(const double* space_in) = 0;

  /** writes design to the vector, beginning with shape_param_ and then aux_design_ */
  virtual int WriteDesignToExtern(double* space_out, bool scaling = true) const = 0;

  /** write gradient out to the vector, appending with shape gradient
   * Sparse and dense! */
  virtual void WriteGradientToExtern(StdVector<double>& out, DesignElement::ValueSpecifier vs, DesignElement::Access access, Function* f, bool scaling = true) = 0;

  /** same as in DesignSpace, setting elements to zero, but also aux elements */
  virtual void Reset(DesignElement::ValueSpecifier vs, DesignElement::Type design = DesignElement::DEFAULT) = 0;

  virtual void WriteBoundsToExtern(double* x_l, double* x_u) const = 0;

  virtual unsigned int GetNumberOfVariables() const = 0;

  virtual int GetNumberOfFeatureMappingVariables() const { return shape_param_.GetSize(); }

  virtual BaseDesignElement* GetDesignElement(unsigned int idx) = 0;

  /** does not go on opt_shape_param_ but on shape_param_ */
  BaseDesignElement* GetFeaturedDesignElement(unsigned int idx) { return shape_param_[idx];}

  virtual void ToInfo(ErsatzMaterial* em) = 0;

  /** Called from DensityFile::ReadErsatzMaterial() with load ersatz material (-x)
   * @param set the set from the density.xml
   * @param lower_violation the maximal violation */
  virtual void ReadDensityXml(PtrParamNode set, double& lower_violation, double& upper_violation) = 0;

  /** This is the variant of Function::Local::SetupVirtualElementMap() for slope constraints on ShapeParamElements.
   * This function is called within Function::Local() constructor, therefore Function::GetLocal() cannot work yet!
   * @param locality just given to assert() it is PREV_AND_NEXT */
  virtual void SetupVirtualShapeElementMap(Function* f, StdVector<Function::Local::Identifier>& virtual_element_map, Function::Local::Locality locality) = 0;

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

  virtual void SetupDesign(PtrParamNode pn) = 0;

  virtual void SetupOptParam() = 0;

  /** Map (distorted) structure to rho (DesignSpace::data). Sets DesignSpace::data.
   *  Shall be called by ReadDesignFromExtern(). */
  virtual void MapFeatureToDensity() = 0;

  /** Takes the density gradients and sums it up on the shape variables.
   *  To be called within WriteGradientToExtern().
   *  @param f the function we add the stuff to the gradient. */
  virtual void MapFeatureGradient(const Function* f) = 0;

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
    static double SetIPGetWeight(const FeaturedDesign* fd, StdVector<double>& ip, int ip_x, int ip_y, int ip_z, int max_order);

    /** maximal diff max_corner_value and min_corner_value for all shapes */
    double MaxDiffCornerValue() const;
  };

  /** shortcut to the dimension (2|3) */
  static unsigned int dim_;

  /** number of elements of rho in x-direction. +1 for nodes! */
  unsigned int nx_ = 0;
  unsigned int ny_ = 0;
  unsigned int nz_ = 0; // 1 for 2D

  /** repeats nx_, ny_ and nz_ */
  Vector<unsigned int> n_;

  /** this is the design_id for the last MapFeatureToDensity() run */
  int mapped_design_ = -1;

  /** Measure MapFeatureToDensity() */
  shared_ptr<Timer> mapping_timer_;

  /** Measure MapFeatureGradient() */
  shared_ptr<Timer> gradient_timer_;

  /** This are the design parameters. This includes fixed elements which are not object
   * to optimization (scpip cannot handle lower bound == upper bound). See opt_shape_param_ */
  StdVector<ShapeParamElement*> shape_param_;

  /** This are the external design parameters which means shape_param_ w/o fixed and optional links */
  StdVector<ShapeParamElement*> opt_shape_param_;

  /** mapping with size of rho to ShapeParamElement pointers to design variables */
  StdVector<Item> map_;

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
