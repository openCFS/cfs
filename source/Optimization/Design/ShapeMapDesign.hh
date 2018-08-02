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
class ShapeParam;


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

  /** Flip dof, means give the complementary dof. For 2D X->Y and Y->X, for 3D this fails! */
  static inline ShapeParamElement::Dof Flip(ShapeParamElement::Dof dof);
  /** Flip dof for center pairs: X,Y -> Z */
  static inline ShapeParamElement::Dof Flip(ShapeParamElement::Dof first, ShapeParamElement::Dof second);

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
   * @param locality just given to assert() it is PREV_AND_NEXT */
  void SetupVirtualShapeElementMap(Function* f, StdVector<Function::Local::Identifier>& virtual_element_map, Function::Local::Locality locality);

  /** For SHAPE_MAP design. Combines NODE and PROFILE. Simple implementation, does not handle symmetry */
  void SetupVirtualMultiShapeElementMap(Function* f, StdVector<Function::Local::Identifier>& virtual_element_map, Function::Local::Locality locality);

  /** Variant of SetupVirtualShapeElementMap() for the periodic constraint which is the first element of a shape minus the last */
  void SetupCyclicVirtualShapeElementMap(Function* f, StdVector<Function::Local::Identifier>& virtual_element_map, Function::Local::Locality locality);

  /** this maps the mesh to a regular lexicographic design representation. For the moment is assumes the mesh to be already
   * lexicographic but this might be extended transparently when required. Used also by LatticeBoltzmannPDE, therefore static!
  @param design_reg extend to vector if necessary!
  @return the size of the mapping as nx, ny, nz with nz = 1 for 2D */
  static Vector<unsigned int> SetupLexicographicMesh(Grid* grid, RegionIdType design_reg, StdVector<int>& elem_to_idx, StdVector<int>& idx_to_elem);

  /** Translates from rho element idx to the smallest corner. The largest corner is out + coord_step_ */
  void MapIdxToCoords(const Vector<unsigned int>& idx, Vector<double>& out) const;

  /** This are principal shape types from where NODE and PROFILE are repeated in BaseDesignElement::Type.
   * a "center" structure contains two nodes is only 3D and contains  */
  typedef enum { CENTER = -1, NODE = 0, PROFILE = 1} Type;

  static Enum<Type> type;

  /** The integration strategy applied forItem::GetOrder() */
  typedef enum { CONSTANT_FULL, FULL_OR_NOTHING, TAILORED } IntStrategy;

  static Enum<IntStrategy> intStrategy;

  /** this describes the continuation of a strucuture in 1D. Either tanh as in Wein, Stingl; 2017 or piecewise linear for cheaper numerical integration. */
  typedef enum { TANH, LINEAR } ShapeFunc;

  Enum<ShapeFunc> shapeFunc;

  ShapeFunc GetShapeFunc() const { return shapeFunc_; }

  const Vector<unsigned int>& GetDiscretization() const { return n_; }

  double GetBeta() const { return beta_; }

  /** convert from ShapeMapDesign::NODE to DesignElement::NODE and the same for PROFILE */
  static BaseDesignElement::Type Convert(Type type);

  /** @see GetNewtonCotes() */
  static StdVector<Vector<double> > newtonCotes;

  /** A ShapeParam corresponds with the xml entries node and profile and stores the upper and lower bounds.
   *  It corresponds with many design elements of type ShapeParamElement.
   *  The 3D center case is composed by two NODE elements - this connection is encoded in center_idx. */
  struct ShapeParam
  {
    ShapeParam() { induced.Reserve(4);} // should not be necessary!

    /** Note that we need a default constructor for StdVector. The symmetry pointers for profile are not set yet!
     * We assume idx to be set to >=0 for nodes or -1 for center
     * @param pn node for the shape
     * @param base reference for profile or center nodes only to copy sym and orientation, ... */
    void ParseAndInit(PtrParamNode pn, ShapeParam* base = NULL);

    /** Copy properties, similar effect as ParseAndInit()
     * @param copy_master_data is used to set data for a 3D slave node, does not touch dof and slave attribute */
    void CopyProperties(const ShapeParam* ref, bool copy_master_data = false);

    /** flip orientation. flip dof and orientation but also adapt and x_sym, y_sym.
     * This is a service function for diagonal induced shapes.
     * @pram center_node for 3D. This will flip both, first and second center node, together. Such that first->orientation == second->orientation
     *                   0 will swap orientation with first->dof and 1 will swap orientation with second->dof */
    void FlipOrientation(int center_node = -1);

    void ToInfo(PtrParamNode pn);

    /** indicate the symmetry data that an additional shape shall be induced?. Checks the orientation of the shape.
     * This is for parallel mirroring, depending on the orientation we map or mirror. Mirroring means reciprocal data (x - > 1-x)
     * A 2d example is a horizontal bar (dof=y, orientation=x) and bottom_up/y_sym: Then there is a parallel bar induced.
     * A 3d center example is standing shape (dof=x,z, orientation=y)
     *  For a a left_right_sym/x_sym: we induce and mirror the x-node. But the z-node is induced but not mirrored!   */
    bool ShallInduceMirrorSymmetry() const;

    /** only 3D (center nodes). See the description of ShallInduceMirrorSymmetry(). It results in a non reciprocal induce setting */
    bool ShallInduceCloneSymmetry() const;

    /** diagonal symmetry induces a new mirror shape of switched direction. Used for square symmetric for band-gaps.
     * The rotation angle is orthogonal to the dof-orientation plane. For 3D center nodes we do not differentiate but do all there flips together*/
    bool ShallInduceDiagonalSymmetry() const;

    /** indicates that only half of the shape is for optimization, the other is mapped. Checks orientation of the shape */
    bool ShallMapHalfShape() const;

    /** indicates that the other center node is a slave of this one */
    bool ShallBeMasterOfSlave() const;

    /** as long as we have no surface, we are !IsCenterShape() */
    bool Is2DShape() const;

    /** are we first or second 3d center node or its profile? */
    bool IsCenterShape() const;

    /** are we a 3D surface (not yet implemented yet) */
    bool IsSurfaceShape() const { return false; }

    /** for the 3D center node case give back the first center node. This can be called for the first center node, the second
     * center node and the common profile node. Note that the implementation is at the end of this file.
     * @return NULL if none of the three cases above holds */
    inline ShapeParam* GetFirstCenterNode();

    inline bool IsFirstCenterNode() const;

    inline bool IsSecondCenterNode() const;

    /** @see GetFirstCenterNode() */
    inline ShapeParam* GetSecondCenterNode();

    /** The reference for 2d shape is the shape, for 3d second the first and for profiles the first shape.
     * We need to implement here */
    int GetReferenceId() const { return type == ShapeMapDesign::PROFILE ? partner->idx : (other_center != NULL ? std::min(other_center->idx, idx) : idx); }

    /** test if the param is part of the shape */
    bool IsPart(const ShapeParamElement* test) const { return (int) test->GetIndex() >= start_param && (int) test->GetIndex() < end_param; }

    /** are we the master shape or an induced sape? */
    bool IsInduced() const;

    /** for debug purpose */
    std::string ToString() const;

    Type type = NODE; // NODE or PROFILE will also be set in Init
    int idx = -2; // number of shape param. -1 for CENTER

    /** a 3D center has two nodes. They point via other_center to each other. NULL means that we are not part of a center.
     * The shape_param with the lower idx is defined as first, the other as second.
     * The links are chained, hence if other_center != NULL then this shall hold: this->other_center->other_center == this.
     * Only a node has other_center. The corresponding profiles need to be NULL two center nodes share a single partner. */
    ShapeParam* other_center = NULL;

    /** the partner of a node is a profile and vice versa. node idx < profile idx. Note that for 3D center nodes two nodes share a
     * single profile partner. The Profile partner points back to the first center node only. */
    ShapeParam* partner = NULL;

    /** the dof of this shape. Note that this shape might be part of a center shape, where only the nodal shapes have dofs
     * A profile has no dof!
     * @see Flip() */
    ShapeParamElement::Dof dof = ShapeParamElement::NOT_SET;

    /** For a dof x the orientation is y. For 3D center pairs X and Y it is Z. A profile inherits the orientation */
    ShapeParamElement::Dof orientation = ShapeParamElement::NOT_SET;

    std::string lower; // to be math parsed with coordinates xi to be used as xi/nx
    std::string upper; //
    std::string value; // initial or fixed

    /** if >= 0 this is the relative bound for the first and last element */
    double clamp = -1.0;

    /** in case we have a symmetry where we induce a shape and mirror it value goes to max - value. Max is the node value */
    double max = 1.0; // fixme an make it smart

    /** in case of 3D bloch optimization we need cubical symmetry and there the two dofs of a center node need to be the same.
     * In the slave case fixed, initial, lower, upper must not be given. */
    bool slave = false;

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

    /** the x_symmetry for dof=y means we copy from left to right. x_symmetry for dof=x means we need to induce an
     * additional shape */
    bool x_sym = false;
    bool y_sym = false;
    bool z_sym = false;
    /** diag always induces a new shape */
    bool diag  = false;

    /** induced shapes have an identification of their kind
     * In 2D this are three combinations. Mirrored (which means the value is reciprocal), diagonal (which means the value is copied)
     * and both (first diagonal and then mirrored) With the master shape (source for the induced) this are four shapes
     * for square symmetry in 2D.
     * For 3D center nodes we have also the case of a copied value. Hence the same information for virtual as reciprocal.
     * The diagonal flipping is always around the axis othogonal to the place dof-orientation. Hence in 2D around the x-axis.*/
    struct Induce
    {
      /** the master shape is the base for the induced one. Then the other flags need to be off. Not that only
       * the master has the induced vector. An each entry has in Induce the description */
      ShapeParam* master = NULL;

      /** are the values induced reciprocal? This is the case for a mirrored shape - e.g. a horizontal shape mirrored at the
       * x-axis. All other cases are copying the value. On comparing the dof and orientation we can identify which of the symmetry
       * cases we have (see above) */
      bool reciprocal = false;
    };

    /** The master shape has no induced information but his induced child shapes have it. */
    Induce induce;

   /** For the master shape only (induce.master = true) this are the links to the induced shapes.
    * The induced shapes have this vector empty!
    * for full symmetry as for bloch modes one structure becomes four structures, here are the three.
    * Note that there is for full symmetry also mapping,  only 1/8 of the data is opt! */
    StdVector<ShapeParam*> induced;

  private:
    /** little helper which reads only bounds and values */
    static void ParseBounds(ShapeParam* target, const PtrParamNode& pn);
    static void ParseSymmetry(ShapeParam* target, const PtrParamNode& pn);
    /** copy attributes found in the envelope center node (symmetry and orientation only) */
    void CopyBaseCenterProperties(ShapeParam* base);
  };

  /** Give ShapeParam, very fast! */
  ShapeParam* GetShape(const ShapeParamElement* spe) { return shape_param_map_[spe->GetIndex()]; }
  ShapeParam* GetShape(const ShapeParamElement& spe) { return GetShape(&spe); }
  const ShapeParam* GetShape(const ShapeParamElement* spe) const { return shape_param_map_[spe->GetIndex()]; }
  const ShapeParam* GetShape(const ShapeParamElement& spe) const { return GetShape(&spe); }

  /** combines our settings for numerical integration */
  struct NumInt
  {
    /** the constructor with the ShapeMap pn. Calls SetTailored conditionally */
    void Init(ShapeMapDesign* smd, PtrParamNode pn, PtrParamNode info);

    /** standard output, does not print warning here */
    void ToInfo(PtrParamNode info) const;

    /** searches tailored_bounds and returns the appropriate tailored_order content */
    int GetTailoredOrder(double max_min) const;

    /** for Item::GetOrder(), however it makes only sense in debug to reduce it. */
    int max_order;

    IntStrategy strategy;

    /** this is the decision value for CloseEnough() */
    double sensitivity;

    /** for tailored order we set the necessary order for max-min bounds. */
    Vector<double> tailored_bounds; // 1e-10, 1e-9, ..., 0.1, 1
    Vector<int>    tailored_order;

    /** for statistics: cells which need integration */
    unsigned int int_cells_cnt = 0;

    /** for statistics: sum of 1D integration on cells we integrate (-> avg. order)  */
    unsigned int int_cells_order_sum = 0;

  private:
    /** determined tailored one by evaluating necessary order for dtanh(beta) on a 1D min(n) grid
     * based on the tanh(beta) max - min bounds per element
     * @param info will add warnings there when max_order is too small or the max newtwon-cotes order is not sufficient */
    void SetTailoredTanh(PtrParamNode info);

    /** For LINEAR shape function we integrate first order. In the grayness we do full integration. The order is determined
     * here based on sensitivity.
     * @see linear_int_order_ */
    void SetLinearIntOrder(PtrParamNode info);

    /** tanh function */
    double Func(double x, double pos) const;

    /** dtanh function */
    double GradFunc(double x, double pos) const;

    /** returns the error for numerical vs analytical integration of dtanh(!)
     * @param order this is the number of integration points. Newton Cotes order 1 has 2 points (1/2,1/2), hence order=2 is minimum
     * @return abs(analytical - numerical)*/
    double IntGradError(double x1, double x2, double pos, int order) const;

    /** find the necessary order, where order=2 means two integration points which is newton cotes order 1
     * @return newtonCotes.GetSize() when the order was not sufficient   */
    int FindOrder(double x1, double x2, double pos, double accuracy) const;

    double beta = -1.0;

    /** product of n */
    unsigned int cells_;

    ShapeMapDesign::ShapeFunc sf_;

    /** what GetTailoredOrder() returns in the LINEAR case with 1st order integration */
    int linear_int_order_ = -1;

    /** the element spacing */
    double h = -1.0;
  };

private:

  /** Setup shape variables from a shape map description. These are abstract structures
   * @param pn a shapeMap element from problem.xml */
  void SetupShapeDesign(PtrParamNode pn);

  /** from the shapes creatate the rho ShapeParamElement variables to be post processed by SetupOptShapeParam() */
  void SetupShapeParam();

  /** From the rho variables from SetupShapeParam() intialize the optimization variables which is complex due to symmetry. */
  void SetupOptShapeParam();

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
  unsigned int DensityIdx(int x, int y) const { return y * nx_ + x; }
  unsigned int DensityIdx(int x, int y, int z) const { return dim_ == 3 ? z * nx_*ny_ + y * nx_ + x : DensityIdx(x, y); }

  /** This is the inverse to DensityIdx */
  void DensityIdx(unsigned int i, Vector<unsigned int>& idx) const;

  /** Index for integration point */
  unsigned int IntPointIdx(int order, int ip_x, int ip_y) const { return ip_y * order+ip_x; }
  unsigned int IntPointIdx(int order, int ip_x, int ip_y, int ip_z) const { return dim_ == 3 ? ip_z * order * order + ip_y * order + ip_x : IntPointIdx(order, ip_x, ip_y); }

  /** Search in shape_ */
  StdVector<ShapeMapDesign::ShapeParam*> FindShape(Type type, ShapeParamElement::Dof dof);

  /** centers have in xml two nodes as children. They are not stored as center in shape_ but need to be searched. */
  StdVector<std::pair<ShapeMapDesign::ShapeParam*, ShapeMapDesign::ShapeParam*> > FindCenters();


  /** slow version of GetShape() when shape_param_map_ is not yet initialized */
  ShapeParam* FindShape(const ShapeParamElement* spe);

  /** Searches shape matching to the function by order.
   * @param function if function design is not PROFILE the first shape is a candidate.
   * @param opt if set only shapes with opt_idx are considered */
  const ShapeParam* FindShape(const Function* f, bool opt) const;

  /** for shape which is either first or second center node (3D!) and a test which is part of first or second
   * center node, return the second center node param. This might be test
   * @return never NULL */
  ShapeParamElement* GetSecondCenterNodeParam(ShapeParam* shape, ShapeParamElement* test);

  /** helper to fill shape_param_
   * @param free_dof for 2D this is Flip(param.dof) (X or Y) for 3D xy, zy, xz this is the current one (X, Y or Z)
   * @param free_idx corresponds to the node counter, not element counter as max free is ny_ and not ny_-1
   * @param start_end indicate the first and last element to enable check for clamped */
  void CreateShapeVariable(const ShapeParam* param, ShapeParamElement::Dof free_dof, int free_idx, bool start_end);

  /** some sanity checks, e.g. volume shall not be linear */
  void CheckPlausibility();

  /** helper for debugging */
  void DumpMap();

  /** Set all Item::corner_val in map_ */
  void EvalAllCornerValues();

  /** Find the corresponding profile variable - needs to identify the shape first :( */
  inline const ShapeParamElement* GetProfile(const ShapeParamElement* node) const;
  inline ShapeParamElement* GetProfile(const ShapeParamElement* node);

  /** do we use a fixed profile? Then opt_shape_param_ is smaller than shape_param_ */
  bool IsProfileFixed() const;

  /** do we have at least a single node shape which is not fixed?.
   * contains a loop, hence cache! */
  bool IsAllNodeFixed() const;
  /** small helper which gives the start index of the element based on type (default, node or profile) of opt_shape_param_ */
  unsigned int GetFirstOptVarIdx(const Function* f) const;

  /** small helper which gives the  index *after* the element based on type (node or profile) shape_param_) */
  unsigned int GetEndOptVarIdx(const Function* f) const;

  /** similar to GetFirstVarIdx() but for shape_ instead of opt_shape_param_. Checks only shapes which do optimization.
   * No fixed and no symmetry
   * @return 0 or num_node_shape_ */
  unsigned int GetFirstShapeIdx(const Function* f) const;

  /* @return num_node_shape_ or shape_.GetSize() */
  unsigned int GetEndShapeIdx(const Function* f) const;

  /** creates for 2D induced symmetry nodes. */
  void Induce2DSymmetryNodes(ShapeParam& ref_node);

  /** create 3d center node symmetry nodes. Works for x_sym, y_sym and z_sym. diag is always for all possible to flips.
   * In contrast to Induce2dSymmetryNodes() the shapes shall be already parsed and initialized, therefore not paramnode necessary.
   * @param cn are the two center nodes */
  void InduceCenterSymmetryNodes(ShapeParam& first, ShapeParam& second);

  /** create a new node, add it as induced to ref_node and set properties */
  ShapeParam* InduceSymmetryNodeHelper(ShapeParam& ref_node);
  void InduceSymmetryNodeHelper(ShapeParam& first, ShapeParam& second);

  /** helper for Levelset */
  double DensityToLevelSet(int x, int y) const;

  /** export a levelset file in vtk format */
  void ExportLevelSet() const;


  /** This are our shape parameters which are blown up to shape_param_. When induced, the ortho induces follows the shape, then the diagonal induced
   * First node then profile, therefore always even size. */
  StdVector<ShapeParam> shape_;

  /** helper for shape_: number of node which is also the first index of the first profile. Is not shape_.GetSize() / 2 when we have 3D center condes .
   * This is NOT the size within shape_param_! */
  int num_node_shapes_ = -1;

  /** This are the shape parameters, defined in ersatzMaterial/shapeMap. This includes fixed elements which are not object
   * to optimization (scpip cannot handle lower bound == upper bound). See opt_shape_param_ */
  StdVector<ShapeParamElement> shape_param_;

  /** This is a map from shape_param_ to shape of size shape_param_. */
  StdVector<ShapeParam*> shape_param_map_;

  /** helper for shape_param_: number of nodes within shape_param_ which not necessarily 1:1 nodes and profiles as 3d center nodes share a profile */
  int num_node_shape_params_ = -1;

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
     * @param shape the source shape. When mapping the original or the induced one otherwise
     * @param reciprocal if the value is shape.max - base->value */
    void AddSymmetryReference(ShapeParamElement* elem, ShapeParam* shape, bool reciprocal);

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

    /** This is static information.
     * The node variables the mapping is based on. Profiles via the nodes
     * The major order is structures, the SPE* vector are the variables for the structure
     * Note that in 3D with center nodes #nodes < num_node_shapes_ as the a and b shapes are collected in the SPE* vector.
     * The major size in 3D is num_node_shapes_ - FindCenterNodes().GetSize() -> 2 center nodes have 4 shapes but 2 * [4*SPE]
     * The SPE* vector is 2 for 2D and 4 for 3D center nodes.*
     * 2D: The ordering of the SPE is for the two nodes which we interpolate between (1-t)*a0 + t*a1
     * 3D: The ordering is a0,b0,a1,b1
     * Filled in SetupShapeParam() */
    StdVector<StdVector<ShapeParamElement*> > nodes;

    /** This is for each "real" node the minimal corner value of all 4 (8 corners). Set by EvalAllCorners() */
    Vector<double> min_corner_value;
    Vector<double> max_corner_value;

    /** Determines based on corner_vals the order of integration for each shape. 0=void, 1=solid, >=2 integration.
     * @param order output
     * @param sensitivity to identify void and solid and set the proper order. This shall come from numerical studies!
     * @param max_order the maximal order. This may conflict with sensitivity
     * @return the maximal number in order. This is <= the parameter max_order */
    int GetOrder(Vector<int>& order, const NumInt& numInt) const;

    /** Set generalized integration points and weights
     * @param ip result: vector [0...1]^dim
     * @param ip_x, ip_y, ip_y index of integration < max_order
     * @return closed Newton-Cotes weight for the 2D/3D point */
    static double SetIPGetWeight(const ShapeMapDesign* smd, StdVector<double>& ip, int ip_x, int ip_y, int ip_z, int max_order);

    /** maximal diff max_corner_value and min_corner_value for all shapes */
    double MaxDiffCornerValue() const;
  };

  /** mapping with size of rho to ShapeParamElement pointers to shape_param_   */
  StdVector<Item> map_;

  /** This structure evaluates a single integration point. It stores exp calculations to reused for gradient calculations. It is aware of
   * Item::nodes */
  struct EvalAtIp
  {
    EvalAtIp(ShapeMapDesign* smd) { Init(smd); };

    /** Empty constructor for vector reason. Call Init()! */
    EvalAtIp() { };

    /** call Init() first in the object lifetime.
     * @param n is ShapeMapDesign::n_ */
    void Init(ShapeMapDesign* smd);

    /** gives the coordinates for evaluation directly.
     * @param idx element identification with map_ but the element is not read.
     * @param ip only used to set a, b, w, r, t but to for x, y */
    double Setup(const StdVector<ShapeParamElement*>& nodes, const Vector<unsigned int>& idx, const StdVector<double>& ip, double beta);

    /** no gradient */
    double ShapeFunc() const;
    /** use order information! 0 = void, 1 = solid, 2 = eval */
    double SmartShapeFunc(int order) const;

    double GradShapeFunc(bool grad_a, bool grad_b, bool grad_w) const;

    /** user order information. See SmartTanh() */
    double SmartGradShapeFunc(int order, bool grad_a, bool grad_b, bool grad_w) const;

  private:
    double Setup2d(const StdVector<ShapeParamElement*>& nodes, const Vector<unsigned int>& idx, const StdVector<double>& ip, double beta);
    double Setup3d(const StdVector<ShapeParamElement*>& nodes, const Vector<unsigned int>& idx, const StdVector<double>& ip, double beta);

    double EvalTanh2d() const;
    double EvalTanh3d() const;
    double EvalTanhGrad2d(bool grad_a, bool grad_w) const;
    double EvalTanhGrad3d(bool grad_a, bool grad_b, bool grad_w) const;

    double EvalLinear2d() const;
    double EvalLinear3d() const;
    double EvalLinearGrad2d(bool grad_a, bool grad_w) const;
    double EvalLinearGrad3d(bool grad_a, bool grad_b, bool grad_w) const;

    // 2D and 3D common
    double a = -1;
    double w = -1;
    /** normalized. x is dof variable in 2D */
    double x = -1;
    /** this is the interpolation variable a = (1-t) * a1 + t * a2 with t = 0 ... 1 */
    double t = -1;

    // extension for 3D
    /** normalized. x is a dof in 3D center nodes and y is b dof */
    double b = -1; // 3D
    double y = -1;
    double r = -1; // 3D: distance (x,y) -> (a,b)

    // shapeFunc == TANH
    double beta = -1;
    double exapw = -1;
    double examw = -1;
    double erw = -1;

    // shapeFunc == LINEAR
    /** linear spacing = 1/element size. The maximal for all directions */
    double h = -1;

    // technical
    int dim = -1;
    ShapeMapDesign* smd_ = NULL;
    /** helper for Setup */
    Vector<double> coord_;
  };

  /** this is the design_id for the last MapShapeToDensit() run */
  int mapped_design_ = -1;

  ShapeFunc shapeFunc_;

  /** controls the boundary. Relates to meter so it also depends on discretication. 30 is a small value (gray) and 70 gives a smaller boundary. */
  double beta_;

  /** MAX means that at each ip we consider only the shape which has the largest rho. Only the gradient of that shape will be considered at that ip.
   * The drawbacks are loss of material at overlaps and doubtful differentiability.
   *
   * TANH_SUM limits via tanh_l( sum(tanh(shape) ) where tanh_l maps to 0..1 with own beta and the beta within sum(tanh(shape)) is halfed.
   *
   * An issue is if the gradients shall be scaled down to match the factor by the cutting of max(sum,1) */
  typedef enum { MAX, TANH_SUM } Overlap;

  /** no need for static */
  Enum<Overlap> overlap;

  /** handles the overlapping of shapes, controls MapShapeToDensity() and has a very strong impact on MapShapeGradietn() */
  Overlap overlap_;


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


  /** number of elements of rho in x-direction. +1 for nodes! */
  unsigned int nx_ = 0;
  unsigned int ny_ = 0;
  unsigned int nz_ = 0; // 1 for 2D

  /** repeats nx_, ny_ and nz_ */
  Vector<unsigned int> n_;

  /** this is the minimal and maximal coordinate */
  Vector<double> coord_min_;
  Vector<double> coord_max_;
  /** this is the edge size */
  Vector<double> coord_step_;

  /** shortcut to the dimension (2,3) */
  static unsigned int dim_;

  /** checks lower and upper when loading from ersatz material */
  bool enforce_bounds_;

  /** set element upper and lower relative to the initial value from load ersatz material.
   * lower is start - half rel_bound and upper is start + half rel_bound.
   * To overwritten by ShapeParam::clamped. smaller 0 disables */
  double relative_node_bound_;
  double relative_profile_bound_;

  /** reference to optimization as we need it in MapShapeGradient() to get the functions */
  Optimization* opt_ = NULL; // set in PostInit() if we have optimization and not only external design for sim

  /** Measure MapShapeToDensity() */
  shared_ptr<Timer> mapping_timer_;

  /** Measure MapShapeGradient() */
  shared_ptr<Timer> gradient_timer_;

  /** do we export a levelset file? Shall become a struct with more options.
   * In the simplest case writes when we export the density. New files with -d else always overwrite. Make smarter! */
  bool export_leveset_;

  NumInt numInt;
};



/** as the inline function is used in DensityFile it needs to be given in the header */
ShapeMapDesign::ShapeParam* ShapeMapDesign::ShapeParam::GetFirstCenterNode()
{
  assert(!(other_center != NULL && other_center->other_center != this));
  assert(!(other_center != NULL && type == PROFILE)); // due to unsymmetry only nodes have the link
  assert(!(other_center != NULL && idx == other_center->idx));
  assert(idx >= 0 && (other_center == NULL || other_center->idx >= 0));
  assert(!(type == PROFILE && partner != NULL));
  if(other_center != NULL)
    return idx < other_center->idx ? this : other_center;
  if(type == PROFILE && partner->other_center != NULL)
    return partner->idx < partner->other_center->idx ? partner : partner->other_center;
  return NULL;
}

} // end of name space


#endif /* OPTIMIZATION_DESIGN_SHAPEMAPDESIGN_HH_ */
