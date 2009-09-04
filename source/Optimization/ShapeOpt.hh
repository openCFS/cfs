#ifndef SHAPEOPT_HH_
#define SHAPEOPT_HH_

#include "Optimization/ParamMat.hh"
#include "Optimization/ShapeDesign.hh"

namespace CoupledField {

  /** This class implements shape optimization, it uses lots of things from Ersatzmaterial
   * but overrides everything concerning derivative calculation
   * it uses ShapeDesign instead of the standard design space 
   * it inherites from ParamMat to be able to use shape & param concurrently */
  class ShapeOpt : public ParamMat {
  public:

    ShapeOpt();

    virtual double CalcVolume(Objective* f, Condition* g, bool derivative, bool normalized = true);

    /** Calculate part of the derivative u1 dK/dShape u2 of compliance and tracking w.r.t. shape */
    void CalcMinusU1dKU2(StdVector<SingleVector*>& u1, StdVector<SingleVector*>& u2, Objective* f, Condition* constraint, double w);
    
    /** Calculate part of the derivative u dF/dShape of compliance and tracking w.r.t. shape */
    void CalcUdF(Excitation& excite, StdVector<SingleVector*>& u, Objective* f, Condition* constraint, double w);

    virtual double CalcCompliance(Excitation& excite, Objective* f, Condition* constraint, bool derivative);

    virtual double CalcTracking(Excitation& excite, Objective* f, Condition* constraint, bool derivative, bool solveproblem = true);

    ShapeDesign* shapedesign;

  private:
    bool alsomatopt_;
  };

}

#endif /*SHAPEOPT_HH_*/
