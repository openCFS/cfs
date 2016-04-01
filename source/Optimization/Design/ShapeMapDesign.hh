/*
 * ShapeMatDesign.hh
 *
 *  Created on: Mar 11, 2016
 *      Author: fwein
 */
#ifndef OPTIMIZATION_DESIGN_SHAPEMAPDESIGN_HH_
#define OPTIMIZATION_DESIGN_SHAPEMAPDESIGN_HH_

#include "Optimization/Design/AuxDesign.hh"

namespace CoupledField
{

class Optimization;

/** Holds the data for ShapeMapping. The map from the parameterized shape to the pseudo densities. */
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

  /** writes design to the vector, prepending with shape variables */
  virtual int WriteDesignToExtern(double* space_out, bool scaling = true) const;

  /** write gradient out to the vector, appending with shape gradient
   * Sparse and dense! */
  virtual void WriteGradientToExtern(StdVector<double>& out, DesignElement::ValueSpecifier vs, DesignElement::Access access, Function* f, bool scaling = true);

  /** same as in DesignSpace, setting elements to zero, but also aux elements */
  virtual void Reset(DesignElement::ValueSpecifier vs, DesignElement::Type design = DesignElement::DEFAULT);

  virtual void WriteBoundsToExtern(double* x_l, double* x_u) const;

  virtual unsigned int GetNumberOfVariables() const;

  virtual BaseDesignElement* GetDesignElement(unsigned int idx);

  virtual void ToInfo(PtrParamNode in, ErsatzMaterial* em);

  /** this maps the mesh to a regular lexicographic design representation. For the moment is assumes the mesh to be already
   * lexicographic but this might be extended transparently when required. Used also by LatticeBoltzmannPDE, therefore static!
  @param design_reg extend to vector if necessary!
  @return the size of the mapping as nx, ny, nz with nz = 1 for 2D */
  static StdVector<unsigned int> SetupLexicographicMesh(Grid* grid, RegionIdType design_reg, StdVector<int>& elem_to_idx, StdVector<int>& idx_to_elem);

  typedef enum { NODE, PROFILE } Type;

  static Enum<Type> type;

  /** store what we read from xml. Will be multiplied to BaseDesignElement in shape_param_ */
  struct ShapeParam
  {
    ShapeParam() { type = NODE; }; // default constructor for StdVector()
    void Init(PtrParamNode pn, unsigned int idx);
    void ToInfo(PtrParamNode pn);

    Type type; // NODE or PROFILE
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
   * Uses Item::ip_param_idx within map_ set by MapShapeToDensity()
   * @param obj true for cost functions false for gradients. Sets mapped_obj_gradient_ or mapped_constr_gradient_ */
  void MapShapeGradient(bool obj);

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
   * @param coords of the density design element
   * @param ip_x in range of order_. 0 for the left side of the element within s1/s2, )order_-1) for the right side */
  double Eval(const ShapeParamElement* s1, const ShapeParamElement* s2, const Matrix<double>& coords, unsigned int ip_x, unsigned int ip_y, bool derivative) const;

  /** Aprroximate the maximal rho for the extremal integration points. Note that one can also save by returning the minimal rho
   * for completely within the structure */
  double ApproxMaxRho(const ShapeParamElement* s1, const ShapeParamElement* s2, const Matrix<double>& coords) const;

  /** tanh performs the smoothing from the mapping
   * @param x is the coordinte (x or y)
   * @param a the shape variable (center of object)
   * @param w half of the thickness of the shape */
  double tanh(double x, double a, double w) const;

  /** derivative of tanh w.r.t. a */
  double d_tanh_da(double x, double a, double w) const;

  /** This are our shape parameters which are blown up to shape_param_ */
  StdVector<ShapeParam> shape_;

  /** This are the shape parameters, defined in ersatzMaterial/shapeMap/shapeParam */
  StdVector<ShapeParamElement> shape_param_;

  /** to conveniently handle the mapping shape param to design */
  struct Item
  {
    /** our Design Element */
    DesignElement* rho;
    /** the variable version of shape_ in the same order*/
    StdVector<ShapeParamElement*> param;

    /** for each integration point order_*order_ (x the fastest variable) the index within param for the largest density from tanh.
     * -1 if this value is too small and the gradient shall be 0.
     * Used to compute the gradient which takes the shape_param for each ip where the corresponding rho is max */
    Vector<int> ip_param_idx;
  };

  /** mapping with size of rho to ShapeParamElement pointers to shape_param_   */
  StdVector<Item> map_;

  /** this is the design_id for the last MapShapeToDensit() run */
  int mapped_design_ = -1;

  /** this is the design_id for last MapShapeGradient() run for objective gradient evaluation.
   * We do double work but are on the save side as the the rho gradients might not be evaluated yet. */
  int mapped_obj_gradient_ = -1;
  int mapped_constr_gradient_ = -1;

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
