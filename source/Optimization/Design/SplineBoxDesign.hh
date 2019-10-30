#ifndef OPTIMIZATION_DESIGN_SPLINEBOXDESIGN_HH_
#define OPTIMIZATION_DESIGN_SPLINEBOXDESIGN_HH_

#include "Optimization/Design/AuxDesign.hh"
#include "Utils/BSpline.hh"

namespace CoupledField {

class SplineBox : public AuxDesign
{
public:

  SplineBox(StdVector<RegionIdType>& regionIds, PtrParamNode pn, ErsatzMaterial::Method method = ErsatzMaterial::NO_METHOD);

  /** Add stuff here to keep out of long constructor */
  void SetupDesign(PtrParamNode pn);

  /** conditionally calls UpdateCoordinates()
   *  @see AuxDesign::ReadDesignFromExtern() */
  virtual int ReadDesignFromExtern(const double* space_in);

  /** writes design to the vector, beginning with shape variables (shape_param_) and then aux_design_ */
  virtual int WriteDesignToExtern(double* space_out, bool scaling = true) const;

  /** write gradient out to the vector, appending with shape gradient
   * Sparse and dense! */
  virtual void WriteGradientToExtern(StdVector<double>& out, DesignElement::ValueSpecifier vs, DesignElement::Access access, Function* f, bool scaling = true);

  /** Set or move one control point */
  void SetControlPoint(int idx, Point coords, bool add = false);

  /** Set or move all control points*/
  void SetControlPoints(StdVector<Point> coords, bool add = false);

  /** Returns the control points of the BSpline curve */
  StdVector<Point> GetControlPoints() { return this->control_points; }

  /** Evaluate the splinebox with current control points and set nodes in FE mesh accordingly */
  void UpdateCoordinates();

  void EvalDerivative();

  /** Write to file */
  void Write();

private:

  /** dimension */
  unsigned int dim;

  /** degree of splines */
  unsigned int degree;

  /** number of control points in each dimension */
  StdVector<unsigned int> num_cp;

  /** number of total control points = prod(num_cp) */
  unsigned int total_num_cp;

  /** vector of control points */
  StdVector<Point> control_points;

  /** design is offset to initial control_points */
  StdVector<BaseDesignElement> shape_param;

  /** Multiple subscripts from linear index */
  void sub2ind(StdVector<unsigned int> size, StdVector<int> sub, unsigned int &ind);

  /** Linear index from multiple subscripts */
  void ind2sub(StdVector<unsigned int> size, unsigned int ind, StdVector<int> &sub);

  void GenerateBasis(Matrix<double> bounding_box, StdVector<BSpline*> bsplines);

  /** Generate a matrix from the bspline basis, s.t. deformed_mesh = matrix * control_points
   *  Size of matrix is nodes inside splinebox times total_num_cp */
  Matrix<double> GetBoxMatrix();

  // StdVector<BSpline*> bsplines; needed?
  StdVector<Matrix<double>> basis;
  StdVector<bool> inside;
  StdVector<StdVector<double>*> unique_coords;
  StdVector<StdVector<double>*> unique_inverse;

};

} //end of namespace

#endif /* OPTIMIZATION_DESIGN_SPLINEBOXDESIGN_HH_ */
