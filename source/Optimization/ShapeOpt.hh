#ifndef SHAPEOPT_HH_
#define SHAPEOPT_HH_

#include "Optimization/ParamMat.hh"
#include "Optimization/Design/ShapeDesign.hh"

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

    virtual double CalcTracking(Excitation& excite, Objective* f, Condition* constraint, bool derivative);

    ShapeDesign* shapedesign;
    
  protected:
    /** Store the results from the forward/adjoint problem. Handles multiple excitations
     * @param read_sol only if true do s.th. if false just pass on to ersatzmaterial
     * @param read_rhs is only interesting for the forward problem
     * @param save_sol set this in the adjoint problem -> see Solution::Read()
     * @param comment is just to LOG_DBG */
    virtual void StorePDESolution(Solutions& solutions, Excitation &excite, Function* f, UInt timestep, bool read_sol, bool read_rhs, bool save_sol, const std::string& comment);

  private:

    bool alsomatopt_;
  };

}

#endif /*SHAPEOPT_HH_*/
