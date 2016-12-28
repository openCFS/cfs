/*
 * ShapeMatDesign.hh
 *
 *  Created on: Mar 11, 2016
 *      Author: fwein
 */
#ifndef OPTIMIZATION_DESIGN_SHAPEMAPDESIGN_HH_
#define OPTIMIZATION_DESIGN_SHAPEMAPDESIGN_HH_

#include "Optimization/Design/AuxDesign.hh"
#include "Optimization/Function.hh"


namespace CoupledField
{

class Optimization;

/** Holds the data for ShapeMapping. The map from the parameterized shape to the pseudo densities.
 * With respect to external design ordering ShapeMapDesign::shape_param_ takes the role of DesignSpace::data
 * and is before AuxDesign::aux_design_.
 * "Internally" however DesignSpace::data is used for the mapping from  ShapeMapDesign::shape_param_ to DesignSpace::data
 * via MadShapeToDensity()
 * Also standard SIMP gradients are computed and stored in DesignSpace::data and added on shape_param_ via MapShapeGradient().*/
class ShapeMapDesign : public AuxDesign
{
public:
  ShapeMapDesign(StdVector<RegionIdType>& regionIds, PtrParamNode pn, ErsatzMaterial::Method method = ErsatzMaterial::NO_METHOD);

  virtual ~ShapeMapDesign() { } ;

  virtual void PostInit(int objectives, int constraints);

  /** @see DesignSpace::ReadDesignFromExtern() */
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

  virtual int GetNumberOfShapeMappingVariables() const { return shape_param_.GetSize(); }

  /** In case DesignSpace::FindDesign() searches for NODE and PROFILE.
   * @return either DesignSpace::FindDesign() or the index within shape_ */
  virtual int FindDesign(DesignElement::Type dt, bool throw_exception = true) const;

  /** goes on opt_shape_param_ only! */
  virtual BaseDesignElement* GetDesignElement(unsigned int idx);

  /** does not go on opt_shape_param_ but on shape_param_ */
  ShapeParamElement* GetShapeMapDesignElement(unsigned int idx) { return &(shape_param_[idx]);}

  virtual void ToInfo(ErsatzMaterial* em);

  /** creates a gnuplot file for the current iteration with the design value and derivatives */
  virtual void WriteGradientFile();

  /** Called from DensityFile::ReadErsatzMaterial() with load ersatz material (-x)
   * @param set the set from the density.xml
   * @param lower_violation the maximal violation */
  void ReadDensityXml(PtrParamNode set, double& lower_violation, double& upper_violation);

  /** This is the variant of Function::Local::SetupVirtualElementMap() for slope constraints on ShapeParamElements.
   * This function is called within Function::Local() constructor, therefore Function::GetLocal() cannot work yet!
   * @param locality just given to assert() it is PREV_AND_NEXT
   * @param phase just given to assert() it is BOTH  */
  void SetupVirtualShapeElementMap(Function* f, StdVector<Function::Local::Identifier>& virtual_element_map, Function::Local::Locality locality, Function::Local::Phase ph);

  /** Variant of SetupVirtualShapeElementMap() for the periodic constraint which is the first element of a shape minus the last */
  void SetupCyclicVirtualShapeElementMap(Function* f, StdVector<Function::Local::Identifier>& virtual_element_map, Function::Local::Locality locality);

  /** this maps the mesh to a regular lexicographic design representation. For the moment is assumes the mesh to be already
   * lexicographic but this might be extended transparently when required. Used also by LatticeBoltzmannPDE, therefore static!
  @param design_reg extend to vector if necessary!
  @return the size of the mapping as nx, ny, nz with nz = 1 for 2D */
  static StdVector<unsigned int> SetupLexicographicMesh(Grid* grid, RegionIdType design_reg, StdVector<int>& elem_to_idx, StdVector<int>& idx_to_elem);

  /** This types are repeated in BaseDesignElement::Type */
  typedef enum { NODE, PROFILE } Type;

  static Enum<Type> type;

  /** The symmetry is defined globally but stored in ShapeParam */
  typedef enum { NONE, MIRROR } Symmetry;

  static Enum<Symmetry> symmetry;

  /** convert from ShapeMapDesign::NODE to DesignElement::NODE and the same for PROFILE */
  static BaseDesignElement::Type Convert(Type type);

  /** convert "x" to 0 and "y" to 1 */
  static int Dof(const std::string& str);

  /** @see Dof() */
  static std::string Dof(int dof);

  /** store what we read from xml. Will be multiplied to BaseDesignElement in shape_param_ */
  struct ShapeParam
  {
    /** Note that we need a default constructor for StdVector. The symmetry pointers for profile are not set yet!
     * @param pn node for the shape
     * @param node reference for profile only to copy sym and orientation, ...
     * @param flip_orientation reinterpret direction for diagonal mirroring. Only for node, not again for profile */
    void Init(PtrParamNode pn, unsigned int idx, ShapeParam* node, bool flip_orientation = false);
    void ToInfo(PtrParamNode pn);

    /** indicate the symmetry data that an additional shape shall be induced?. Checks the orientation of the shape.
     * This is for orthogonal mirroring, depending on the orientation we map or mirror orthogonal */
    bool ShallInduceOrthogonalSymmetry() const;

    /** diagonal symmetry induces a new mirror shape of switched direction. Used for square symmetric for band-gaps */
    bool ShallInduceDiagonalSymmetry() const;

    /** indicates that only half of the shape is for optimization, the other is mapped. Checks orientation of the shape */
    bool ShallMapHalfShape() const;

    /** for debug purpose */
    std::string ToString() const;

    Type type = NODE; // NODE or PROFILE will also be set in Init
    int idx = -1;
    int dof = -1; // x =0, y = 1, z = 2,
    double lower = -1.0;
    double upper = -1.0;
    double value = -1.0; // initial or fixed

    /** if >= 0 this is the relative bound for the first and last element */
    double clamp = -1.0;

    /** in case we have a symmetry where we induce a shape and mirror it value goes to max - value. Max is the node value */
    double max = 1.0; // fixme an make it smart

    /** subject to optimization or fixed */
    bool fixed = false;

    /** this stores the reference to shape_param_ */
    int start_param = -1;
    /** this stores the reference to the first parameter within opt_shape_param_. -1 if no design.
     * It becomes complex due to fixed and partial symmetry */
    int start_opt = -1;

    /** this end does not reflect symmetry */
    int end_param = -1;
    /** this is the end of the optimization, reflects symmetry and is -1 if no design */
    int end_opt = -1;

    /** the orientiation of the shape is in 2D the complementary */
    int orientation = -1;

    /** the x_symmetry for dof=y means we copy from left to right. x_symmetry for dof=x means we need to induce an
     * additional shape */
    Symmetry x_sym = NONE;
    Symmetry y_sym = NONE;
    /** diag always induces a new shape */
    Symmetry diag  = NONE;

    /** a shape with dof x and x_symmetry mirror or diagonal mirror means that an additional mirror induced shape needs to be inserted.
     * This bool indicates if this shape is such a mirror induced shape */
    bool sym_induced = false;

   /** For a non-sym_induced shape this are the links to the induced shapes.
    * for full symmetry as for bloch modes one structure becomes four structures, here are the three.
    * Note that there is for full symmetry also mapping,  only 1/8 of the data is opt! */
    ShapeParam* sym_ortho = NULL; // same direction
    ShapeParam* sym_diag  = NULL; // changed direction to base shape
    ShapeParam* sym_diag_ortho = NULL; // first diag, then ortho -> changed direction to base shape

  };

protected:

  /** Setup shape design variables from a shape map description
   * @param pn a shapeMap element from problem.xml */
  void SetupShapeDesign(PtrParamNode pn);


  /** map shape design to rho (DesignSpace::data). Sets DesignSpace::data. Shall be called by ReadDesignFromExtern().
   * Sets Item::ip_param_idx within map_ for fast MapShapeGradient() */
  void MapShapeToDensity();

  /** Takes the density gradients and sums it up on the shape variables using map_. To be called within WriteGradientToExtern().
   * All the tanh stuff is repeatedly calculated for each function. However WriteGradientToExtern(f) is not called after all simp function gradients are set,
   * therefore we cannot cache it.
   * Uses Item::ip_param_idx within map_ set by MapShapeToDensity()
   * @param f the function we add the stuff to the gradient. */
  void MapShapeGradient(const Function* f);

  /** Index of rho in DesignSpace::data() by element coordinate */
  unsigned int DensityIdx(int x, int y) const { return nx_ * y + x; }

  /** Index for integration point */
  unsigned int IntPointIdx(unsigned int ip_x, unsigned int ip_y) const { return ip_y*order_+ip_x; }

  /** Search in shape_ */
  StdVector<ShapeMapDesign::ShapeParam*> FindShape(Type type, int dof);

  /** search for the corresponding shape */
  ShapeParam* FindShape(const ShapeParamElement* spe);

  /** helper to fill shape_param_
   * @param free corresponds to the node counter, not element counter as max free is ny_ and not ny_-1
   * @param start_end indicate the first and last element to enable check for clamped */
  void CreateShapeVariable(const ShapeParam* param, int free, bool start_end);

  /** helper for debugging */
  void DumpMap();

  /** Evaluate the function at the given integration point. The integration mapping is cartesian oriented
   * @oaram s1 and s2 are both nodes! Eval finds the profiles by itself!
   * @param coords of the density design element
   * @param beta @see tanh()
   * @param ip_x in range of order_. 0 for the left side of the element within s1/s2, )order_-1) for the right side
   * @param grad_a false for tanh, true for d_tanh_da
   * @param grad_w false for tanh, true for d_tanh_dw. */
  double Eval(const ShapeParamElement* s1, const ShapeParamElement* s2, const Matrix<double>& coords, double beta, unsigned int ip_x, unsigned int ip_y, bool grad_a, bool grad_w) const;

  /** Decides if we element given by the coordinates is close enough to the nodal shapes (profiles found implicitly) such that
   * it is worth to consider them. */
  bool CloseEnough(const ShapeParamElement* s1, const ShapeParamElement* s2, const Matrix<double>& coords) const;

  /** tanh performs the smoothing from the mapping
   * @param beta for overlap = max or open_sum use 2*beta_ for historical reasons
   * @param x is the coordinate (x or y)
   * @param a the shape variable (center of object)
   * @param w half of the thickness/profile of the shape */
  double tanh(double beta, double x, double a, double w) const;

  /** derivative of tanh w.r.t. a */
  double d_tanh_da(double beta, double x, double a, double w) const;

  /** derivative of tanh w.r.t. w which is the half profile */
  double d_tanh_dw(double beta, double x, double a, double w) const;

  /** small helper */
  ShapeParam& GetProfile(const ShapeParam& node) { return shape_[node.idx + num_node_shapes_]; }

  const ShapeParamElement* GetProfile(const ShapeParamElement* node) const { return &shape_param_[node->GetIndex() + num_node_shape_params_]; }
  ShapeParamElement* GetProfile(const ShapeParamElement* node) { return &shape_param_[node->GetIndex() + num_node_shape_params_]; }

  /** do we use a fixed profile? Then opt_shape_param_ smaller shape_param_ */
  bool IsProfileFixed() const;

  /** small helper which gives the start index of the element based on type (default, node or profile) (shape_param_ or opt_shape_param_)
   * @param opt if false is based on shappe_param_ if true is based on opt_shape_param_ which is the same if we have no fixed profile AND if we have no symmetries! */
  unsigned int GetFirstVarIdx(const Function* f, bool opt) const;

  /** small helper which gives the  index *after* the element based on type (node or profile) shape_param_) */
  unsigned int GetEndVarIdx(const Function* f, bool opt) const;

  /** similar to GetFirstVarIdx() but for shape_ instead of shape_param_ */
  unsigned int GetFirstShapeIdx(const Function* f, bool opt) const;

  unsigned int GetEndShapeIdx(const Function* f, bool opt) const;

  /** This are our shape parameters which are blown up to shape_param_. Whend induced, the ortho induces follows the shape, then the diagonal induced
   * First node then profile, therefore always even size */
  StdVector<ShapeParam> shape_;

  /** helper for shape_: number of node which is also the first index of the first profile. equals shape_.GetSize() / 2.
   * This is NOT the size within shape_param_! */
  int num_node_shapes_ = -1;

  /** This are the shape parameters, defined in ersatzMaterial/shapeMap. This includes fixed elements which are not object
   * to optimization (scpip cannot handle lower bound == upper bound). See opt_shape_param_ */
  StdVector<ShapeParamElement> shape_param_;

  /** helper for shape_param_: number of nodes within shape_param_ which is  shape_param_.GetSize() / 2 */
  int num_node_shape_params_ = -1;

  /** same as num_node_shape_params_ but based on opt_shape_param_ */
  int num_node_opt_shape_params_ = -1;

  /** symmetry means that fewer data is in opt_shape_param_ but all is in shape_param_. There are different ways to map
   * from opt_shape_param_ to shape_param_ which is stored in opt_sym_param_ */
  struct ElementSymmetry
  {
    ElementSymmetry(ShapeParamElement* base);

    /** short cut to check if we have any symmetry */
    bool HasSymmetry() const { return !hidden.IsEmpty(); }

    /** apply the value for the opt element to the symmetry elements considering all special cases */
    void ApplyDesign();

    /** returns all plain gradients for the sym elements or 0.0 if no sym */
    double GetPlainSymGradient(const Function* f) const;

    void Reset(DesignElement::ValueSpecifier vs);

    void PostInit(int objectives, int constraints);

    /** @param elem the virtual element
     * @param shape the original shape (when mapping) or the induced one
     * @param map do we map or only mirror
     * @param reciprocal if the value is shape.max - base->value */
    void AddSymmetryReference(ShapeParamElement* elem, ShapeParam* shape, bool map, bool reciprocal);

    /** for logging
     * @param grad add gradient details */
    std::string ToString(bool grad=false) const;

    /** the symmetry information. */
    struct Virtual
    {
      /** just for StdVector */
      Virtual() {};

      ShapeParamElement* elem  = NULL;  // the virtual element from opt_shape_param_
      ShapeParam*        shape = NULL;  // the original shape (when only map) or the induced shape
      bool               map   = false;    // just for debugging, do we stem from debugging information?
      bool               reciprocal = false; // copy the original value or shape.max - val, also for gradient!
    };

    /** This is the mappings we have */
    StdVector<Virtual> hidden;

    /** This is the element we build our symmetries to if we do */
    ShapeParamElement* base = NULL;
  };

  /** Only some of the ShapeParamElements are really for optimization. Skipped are fixed shapes and symmetric shapes/elements. */
  struct OptVar
  {
    /** for StdVector only! */
    OptVar() { };
    ~OptVar();

    void Init(ShapeParamElement* elem);

    /** Link to the shape param element which is actually for optimization. We don't own this! */
    ShapeParamElement* elem = NULL;
    /** in case of symmetry here the links to the hidden ShapeParamElements are contained. We own this! */
    ElementSymmetry*   sym = NULL;
  };

  /** This are the external shape param variables which means shape_param_ w/o fixed and optional symmetry links */
  StdVector<OptVar> opt_shape_param_;

  /** to conveniently handle the mapping shape param to design */
  struct Item
  {
    /** our Design Element */
    DesignElement* rho;

    /** the node variables the mapping is based on.
     * ShapeParamElement is connected to ShapeParam in shape_ (same order).
     * nodes are sufficient as the profile is  */
    StdVector<ShapeParamElement*> nodes;

    /** this is the current subset of nodes which gives the relevant shapes.
     * We wouldn't need it for overlap_ == MAX as ip_param_idx has this information more detailed.
     * Set in MapShapeDesign() when checking for CloseEnough() and used in MapShapeGradient().
     * Has maximal size of nodes */
    StdVector<int> relevant_nodes;

    /** For overlap == max only!
     * for each integration point order_*order_ (x the fastest variable) the index within param for the largest density from tanh.
     * -1 if this value is too small and the gradient shall be 0.
     * Used to compute the gradient which takes the shape_param for each ip where the corresponding rho is max
     * Has size order_ * order_ */
    StdVector<int> ip_param_idx;

    /** Stores da_norm to prevent recomputing. For 3D make sure you have not too much integration points!
     * Size and access is relevant_nodes * order_ * order */
    StdVector<double> da_norm_cache;
    StdVector<double> dw_norm_cache;
  };

  /** mapping with size of rho to ShapeParamElement pointers to shape_param_   */
  StdVector<Item> map_;

  /** this is the design_id for the last MapShapeToDensit() run */
  int mapped_design_ = -1;

  /** controls the boundary. Relates to meter so it also depends on discretication. 30 is a small value (gray) and 70 gives a smaller boundary. */
  double beta_;

  /** this is the decision value for CloseEnough() */
  double sensitivity_;

  /** MAX means that at each ip we consider only the shape which has the largest rho. Only the gradient of that shape will be considered at that ip.
   * The drawbacks are loss of material at overlaps and doubtful differentiability.
   *
   * OPEN_SUM means sum tanh(shape) which will be 2.0 at overlaps > rho_max! However there are no issues with differentiability -> use it for reference only!
   *
   * TANH_SUM limits via tanh_l( sum(tanh(shape) ) where tanh_l maps to 0..1 with own beta and the beta within sum(tanh(shape)) is halfed.
   *
   * An issue is if the gradients shall be scaled down to match the factor by the cutting of max(sum,1) */
  typedef enum { MAX, OPEN_SUM, TANH_SUM } Overlap;

  /** no need for static */
  Enum<Overlap> overlap;

  /** small helper for tanh_sum */
  struct TanhSum
  {
  public:
    TanhSum();

    /** tanh for limiting with exactly 0->0 and >=1 -> 1!!
     * @param x is the sum of tanh(shape1) + tanh(shape2) on a specific integration point */
    double map(double x);

    /** the derivative of map
     * @param x is sum of tanh(shape1) + tanh(shape2)
     * @param dx is d_tanh_a(shape) for the only shape  */
    double d_map(double x, double dx);

  private:
    /** tanh for limiting with approx 0->0 and >=1 -> 1 */
    double tanh(double x);
    double beta = 11;
    double offset = -1.0; // by this we correct tanh(x=0,1) such that map(x=0,1) is 0,1
    double scale = -1.0;  // works with offset
  } tanh_sum_;



  /** handles the overlapping of shapes, controls MapShapeToDensity() and has a very strong impact on MapShapeGradietn() */
  Overlap overlap_;

  /** number of elements of rho in x-direction. +1 for nodes! */
  unsigned int nx_ = 0;
  unsigned int ny_ = 0;
  unsigned int nz_ = 0; // 1 for 2D

  /** this is the order of integration with order^dim evaluations.
   * Smaller 5 has poor numerics, larger 10 might become expensive */
  unsigned int order_;

  /** order_ * order */
  unsigned int order_order_;

  /** shortcut to the dimension (2,3) */
  unsigned int dim_;

  /** checks lower and upper when loading from ersatz material */
  bool enforce_bounds_;

  /** set element upper and lower relative to the initial value from load ersatz material.
   * lower is start - half rel_bound and upper is start + half rel_bound.
   * To overwritten by ShapeParam::clamped. smaller 0 disables */
  double relative_node_bound_;
  double relative_profile_bound_;

  /** reference to optimization as we need it in MapShapeGradient() to get the functions */
  Optimization* opt_ = NULL; // set in PostInit() if we have optimization and not only external design for sim
};

} // end of name space


#endif /* OPTIMIZATION_DESIGN_SHAPEMAPDESIGN_HH_ */
