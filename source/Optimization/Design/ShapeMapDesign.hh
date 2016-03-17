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

/** Holds the data for ShapeMapping. The map from the parameterized shape to the pseudo densities. */
class ShapeMapDesign : public AuxDesign
{
public:
  ShapeMapDesign(StdVector<RegionIdType>& regionIds, PtrParamNode pn, ErsatzMaterial::Method method = ErsatzMaterial::NO_METHOD);

  virtual ~ShapeMapDesign() { } ;

  /** this maps the mesh to a regular lexicographic design representation. For the moment is assumes the mesh to be already
   * lexicographic but this might be extended transparently when required. Used also by LatticeBoltzmannPDE, therefore static!
  @param design_reg extend to vector if necessary!
  @return the size of the mapping as nx, ny, nz with nz = 1 for 2D */
  static StdVector<unsigned int> SetupLexicographicMesh(Grid* grid, RegionIdType design_reg, StdVector<int>& elem_to_idx, StdVector<int>& idx_to_elem);

  /** map shape design to rho (DesignSpace::data). Write DesignSpace::data */
  void MapShapeToDensity();

  typedef enum { NODE, PROFILE } Type;

  static Enum<Type> type;

  /** store what we read from xml. Will be multiplied to BaseDesignElement in shape_param_ */
  struct ShapeParam
  {
    ShapeParam() { type = NODE; }; // default constructor for StdVector()
    void Init(PtrParamNode pn);

    Type type; // NODE or PROFILE
    int nr = -1;
    int dof = -1; // x =0, y = 1, z = 2,
    double lower = -1.0;
    double upper = -1.0;
    double value = -1.0; // initial or fixed
    bool fixed = false;
  };

protected:

  /** Setup shape design variables from a shape map description
   * @param pn a shapeMap element from problem.xml */
  void SetupShapeDesign(PtrParamNode pn);

  /** Index of rho in DesignSpace::data() by element coordinate */
  unsigned int DensityIdx(int x, int y) const { return nx_ * y + x; }

  /** This are the shape parameters, defined in ersatzMaterial/shapeMap/shapeParam */
  StdVector<ShapeParamElement> shape_param_;

  /** Search in shape_ */
  StdVector<ShapeParam> FindShape(Type type, int dof);

  /** helper to fill shape_param_
   * @param free corresponds to the node counter, not element counter as max free is ny_ and not ny_-1*/
  void CreateShapeVariable(ShapeParam& param, int free);

  /** helper for debugging */
  void DumpMap();

  /** tanh performs the smoothing from the mapping
   * @param x is the coordinte (x or y)
   * @param a the shape variable (center of object)
   * @param w half of the thickness of the shape */
  double tanh(double x, double a, double w) const;


  /** This are our shape parameters which are blown up to shape_param_ */
  StdVector<ShapeParam> shape_;

  /** to conveniently handle the mapping shape param to design */
  struct Item
  {
    DesignElement* rho;
    StdVector<ShapeParamElement*> param;
  };

  /** mapping with size of rho */
  StdVector<Item> map_;

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
};

} // end of name space


#endif /* OPTIMIZATION_DESIGN_SHAPEMAPDESIGN_HH_ */
