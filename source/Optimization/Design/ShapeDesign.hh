#ifndef SHAPEDESIGN_HH_
#define SHAPEDESIGN_HH_

#include "Optimization/Design/DesignSpace.hh"

namespace CoupledField{
  /** This adds the functionality needed for Shape Optimization to the DesignSpace:
   * Storage of the Design Parameters and Derivatives
   * as well as storage of the mesh-deformations
   * it is kept in a way to allow extension to Shape & Material Optimization, so there are some unnecessary parameters at times */
  class ShapeDesign : public DesignSpace {

  public:

    ShapeDesign(StdVector<RegionIdType>& regionIds, ParamNodeList& design, ParamNodeList& transfer, ParamNodeList& result,
        ErsatzMaterial::Method method = ErsatzMaterial::NO_METHOD);

    virtual ~ShapeDesign();

    /** Overwrite the virtual base function.
     * We have always a non-regular grid for SHAPE_OPT and SHAPE_PARAM_MAT */
    virtual bool IsRegular() const
    {
      return false;
    }

    void Configure(PtrParamNode pn, int objectives, int constraints);

    /** Reads design from optimizer, note that first nshapeparams_ entries are the Shape parameters,
     * the rest may be other design parameters */
    int ReadDesignFromExtern(const double* space_in);

    /** applies the calculated design variables to the node offsets */
    void UpdateCoordinates();

    /** writes design to the vector, prepending with shape variables */
    int WriteDesignToExtern(double* space_out, bool scaling = true) const;

    /** write gradient out to the vector, prepending with shape gradient */
    void WriteDenseGradientToExtern(StdVector<double>& out, DesignElement::ValueSpecifier vs,
                               DesignElement::Access access, Condition* constraint = NULL, bool scaling = true) const;

    /** write the shape gradient part */
    void WriteShapeGradientToExtern(StdVector<double>& out, Condition* constraint) const;
    
    /** same as in DesignSpace, setting elements to zero, but also shape elements */
    void Reset(DesignElement::ValueSpecifier vs, DesignElement::Type design = DesignElement::DEFAULT);

    void WriteBoundsToExtern(double* x_l, double* x_u) const;

    virtual unsigned int GetNumberOfVariables() const;

    int GetNumberOfShapeParameters() const {
      return nshapeparams_;
    }
    
    /** return whether this element does depend on any deformations at all */
    bool IsElemDependentAtAll(const StdVector<UInt>& connect);

    /** Get the derivative of the CornerCoords of one element in direction of parameter-th shape param */
    bool GetElemNodesCoordDerivative(Matrix<Double> & coordMat, const StdVector<UInt> & connect, const int parameter);

    void AddShapeDerivatives(Objective* objective, Condition* constraint, StdVector<double>& d, double weight);

    /** return whether also material optimization (SIMP, ParamOpt, ...?) is done */
    bool AlsoMatOpt() const {
      return(alsomatopt_);
    }

  private:
    bool alsomatopt_;

    unsigned int nshapeparams_;

    StdVector<BaseDesignElement> shapeparams_;

    /** deformation-dependency tensor (3rd order)
     * for every node in grid this is a matrix, with dim rows and nshapeparams_ columns */
    StdVector<Matrix<double>* > nodedeformations_;

    UInt dim_;
    
    /** factor for scaling the shapeparameters, this scales the gradient compared to material parameters gradient */
    double scaling;

  };

}

#endif /*SHAPEDESIGN_HH_*/
