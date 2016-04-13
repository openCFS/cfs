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

  virtual BaseDesignElement* GetDesignElement(unsigned int idx);

  virtual void ToInfo(PtrParamNode in, ErsatzMaterial* em);

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

  typedef enum { NODE, PROFILE } Type;

  static Enum<Type> type;

  /** convert from ShapeMapDesign::NODE to DesignElement::NODE and the same for PROFILE */
  BaseDesignElement::Type Convert(Type type) const;

  /** store what we read from xml. Will be multiplied to BaseDesignElement in shape_param_ */
  struct ShapeParam
  {
    // note that we would need a default constructor for StdVector
    void Init(PtrParamNode pn, unsigned int idx);
    void ToInfo(PtrParamNode pn);

    Type type = NODE; // NODE or PROFILE will also be set in Init
    int idx = -1;
    int dof = -1; // x =0, y = 1, z = 2,
    double lower = -1.0;
    double upper = -1.0;
    double value = -1.0; // initial or fixed
    bool fixed = false;

    /** this stores the reference to shape_param_ */
    int start_param = -1;
    int end_param = -1;
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

  /** Search in shape_ */
  StdVector<ShapeParam*> FindShape(Type type, int dof);

  /** helper to fill shape_param_
   * @param free corresponds to the node counter, not element counter as max free is ny_ and not ny_-1*/
  void CreateShapeVariable(const ShapeParam* param, int free);

  /** helper for debugging */
  void DumpMap();

  /** Evaluate the function at the given integration point. The integration mapping is cartesian oriented
   * @oaram s1 and s2 are both nodes! Eval finds the profiles by itself!
   * @param coords of the density design element
   * @param ip_x in range of order_. 0 for the left side of the element within s1/s2, )order_-1) for the right side
   * @param grad_a false for tanh, true for d_tanh_da
   * @param grad_w false for tanh, true for d_tanh_dw. */
  double Eval(const ShapeParamElement* s1, const ShapeParamElement* s2, const Matrix<double>& coords, unsigned int ip_x, unsigned int ip_y, bool grad_a, bool grad_w) const;

  /** Aprroximate the maximal rho for the extremal integration points. Note that one can also save by returning the minimal rho
   * for completely within the structure */
  double ApproxMaxRho(const ShapeParamElement* s1, const ShapeParamElement* s2, const Matrix<double>& coords) const;

  /** tanh performs the smoothing from the mapping
   * @param x is the coordinte (x or y)
   * @param a the shape variable (center of object)
   * @param w half of the thickness/profile of the shape */
  double tanh(double x, double a, double w) const;

  /** derivative of tanh w.r.t. a */
  double d_tanh_da(double x, double a, double w) const;

  /** derivative of tanh w.r.t. w which is the half profile */
  double d_tanh_dw(double x, double a, double w) const;

  /** small helper */
  ShapeParam& GetProfile(const ShapeParam& node) { return shape_[node.idx + num_node_shapes_]; }

  const ShapeParamElement* GetProfile(const ShapeParamElement* node) const { return &shape_param_[node->GetIndex() + num_node_shape_params_]; }
  ShapeParamElement* GetProfile(const ShapeParamElement* node) { return &shape_param_[node->GetIndex() + num_node_shape_params_]; }

  /** do we use a fixed profile? Then opt_shape_param_ smaller shape_param_ */
  bool IsProfileFixed() const;

  /** small helper which gives the start index of the element based on type (default, node or profile) (shape_param_ or opt_shape_param_)
   * @param opt if false is based on shappe_param_ if true is based on opt_shape_param_ which is the same if we have no fixed profile */
  unsigned int GetFirstVarIdx(const Function* f, bool opt) const;

  /** small helper which gives the  index *after* the element based on type (node or profile) shape_param_) */
  unsigned int GetEndVarIdx(const Function* f, bool opt) const;

  /** similar to GetFirstVarIdx() but for shape_ instead of shape_param_ */
  unsigned int GetFirstShapeIdx(const Function* f, bool opt) const;

  unsigned int GetEndShapeIdx(const Function* f, bool opt) const;

  /** This are our shape parameters which are blown up to shape_param_.
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

  /** This are the external shape param variables which means shape_param_ w/o fixed */
  StdVector<ShapeParamElement*> opt_shape_param_;

  /** to conveniently handle the mapping shape param to design */
  struct Item
  {
    /** our Design Element */
    DesignElement* rho;
    /** the node variables the mapping is based on.
     * ShapeParamElement is connected to ShapeParam in shape_ (same order).
     * nodes are sufficient as the profile is  */
    StdVector<ShapeParamElement*> nodes;

    /** for each integration point order_*order_ (x the fastest variable) the index within param for the largest density from tanh.
     * -1 if this value is too small and the gradient shall be 0.
     * Used to compute the gradient which takes the shape_param for each ip where the corresponding rho is max */
    Vector<int> ip_param_idx;
  };

  /** mapping with size of rho to ShapeParamElement pointers to shape_param_   */
  StdVector<Item> map_;

  /** this is the design_id for the last MapShapeToDensit() run */
  int mapped_design_ = -1;

  double beta_;

  /** number of elements of rho in x-direction. +1 for nodes! */
  unsigned int nx_ = 0;
  unsigned int ny_ = 0;
  unsigned int nz_ = 0; // 1 for 2D

  /** this is the order of integration with order^dim evaluations.
   * Smaller 5 has poor numerics, larger 10 might become expensive */
  const unsigned int order_;

  /** shortcut to the dimension (2,3) */
  unsigned int dim_;

  /** reference to optimization as we need it in MapShapeGradient() to get the functions */
  Optimization* opt_ = NULL; // set in PostInit() if we have optimization and not only external design for sim
};

} // end of name space


#endif /* OPTIMIZATION_DESIGN_SHAPEMAPDESIGN_HH_ */
