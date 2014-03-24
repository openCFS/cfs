#ifndef SHAPEDESIGN_HH_
#define SHAPEDESIGN_HH_

#include <stddef.h>

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "General/defs.hh"
#include "General/environment.hh"
#include "Optimization/Design/DesignElement.hh"
#include "Optimization/Design/AuxDesign.hh"
#include "Optimization/ErsatzMaterial.hh"
#include "Utils/StdVector.hh"

namespace CoupledField{
  /** Extends AuxDesign to have shape parameters:
   * Storage of the Design Parameters and Derivatives
   * as well as storage of the mesh-deformations
   * it is kept in a way to allow extension to Shape & Material Optimization, so there are some unnecessary parameters at times */
template <class TYPE> class Matrix;

  class ShapeDesign : public AuxDesign {

  public:

    ShapeDesign(StdVector<RegionIdType>& regionIds, PtrParamNode pn, ErsatzMaterial::Method method = ErsatzMaterial::NO_METHOD);

    virtual ~ShapeDesign();

    /** Overwrite the virtual base function.
     * We have always a non-regular grid for SHAPE_OPT and SHAPE_PARAM_MAT */
    virtual bool IsRegular() const { return false; }

    /** Read parameters from xml file and initialize our design elements */
    void Configure(PtrParamNode pn, int objectives, int constraints);

    /** change the node positions in the mesh to correspond to the current design,
     * called from ReadDesignFromExtern */
    void UpdateCoordinates();

    /** conditionally calls UpdateCoordinates()
     *  @see AuxDesign::ReadDesignFromExtern() */
    int ReadDesignFromExtern(const double* space_in);
    
    /** return whether this element does depend on any deformations at all 
     * @param connect nodes of this element */
    bool IsElemDependentAtAll(const StdVector<UInt>& connect);

    /** Get the derivative of the CornerCoords of one element in direction of parameter-th shape param 
     * @param coordMat returns the derivative of the CornerCoords 
     * @param connect nodes of the element we want the derivative of 
     * @param parameter index to the shape parameter we calculate the derivative with respect to */
    bool GetElemNodesCoordDerivative(Matrix<Double> & coordMat, const StdVector<UInt> & connect, const int parameter);

    /** return whether also material optimization (SIMP, ParamOpt, ...?) is done */
    bool AlsoMatOpt() const {
      return(alsomatopt_);
    }
    
    /** storage type for shape constraints 
     * a constraint has the type factor[0] * param[0] - factor[1] * param[1] <= 1.0 */
    struct ShapeConstraint {
      int param[2];
      double factor[2];
    };

    /** return the shapeconstraints, this is called from Function.cc */
    StdVector<ShapeConstraint>& GetShapeConstraints();
    
  private:

    /** deformation-dependency tensor (3rd order)
     * for every node in grid this is a matrix, with dim rows and nshapeparams_ columns */
    StdVector<Matrix<double>* > nodedeformations_;

    /** shortcut to dimension */
    UInt dim_;
    
    /** Vector containing all shape constraints */
    StdVector<ShapeConstraint> shapeconstraints_;
  };

}

#endif /*SHAPEDESIGN_HH_*/
