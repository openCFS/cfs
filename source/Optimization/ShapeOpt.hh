#ifndef SHAPEOPT_HH_
#define SHAPEOPT_HH_

#include <stddef.h>
#include <string>

#include "General/defs.hh"
#include "MatVec/matrix.hh"
#include "Optimization/ParamMat.hh"

namespace CoupledField {

  /** This class implements shape optimization, it uses lots of things from Ersatzmaterial
   * but overrides everything concerning derivative calculation
   * it uses ShapeDesign instead of the standard design space 
   * it inherites from ParamMat to be able to use shape & param concurrently */
class Condition;
class Excitation;
class Function;
class Objective;
class ShapeDesign;
template <class TYPE> class Vector;

  class ShapeOpt : public ParamMat {
  public:

    ShapeOpt();

    virtual double CalcVolume(Objective* f, Condition* g, bool derivative, bool normalized = true);

    /** Calculate part of the derivative u1 dK/dShape u2 of compliance and tracking w.r.t. shape */
    void CalcMinusU1dKU2(Solutions& forward, Solutions& adjoint, Objective* f, Condition* constraint, const Matrix<double>* tensor_diff = NULL);
    
    /** Calculate part of the derivative u dF/dShape of compliance and tracking w.r.t. shape */
    void CalcUdF(Solutions& adjoint, Objective* f, Condition* constraint, double w = 1.0);

    virtual double CalcCompliance(Excitation& excite, Objective* f, Condition* constraint, bool derivative);

    virtual double CalcTracking(Excitation& excite, Objective* f, Condition* constraint, bool derivative);
    
    virtual void CalcHomogenizedTrackingGradient(const Matrix<double>& target, const Matrix<double>& hom, Function* f);
    
    virtual Matrix<double> CalcHomogenizedTensor();

    ShapeDesign* shapedesign;
    
  protected:
    /** Store the results from the forward/adjoint problem. Handles multiple excitations
     * @param read_sol only if true do s.th. if false just pass on to ersatzmaterial
     * @param read_rhs is only interesting for the forward problem
     * @param save_sol set this in the adjoint problem -> see Solution::Read()
     * @param comment is just to LOG_DBG */
    virtual void StorePDESolution(Solutions& solutions, Excitation &excite, Function* f, UInt timestep, bool read_sol, bool read_rhs, bool save_sol, DERIVType derivative, const std::string& comment);

    /** Subtract the current Testdisplacement from a given vector */
    inline void SubtractTestDisplacement(unsigned int idx, Matrix<double>& CornerCoords, Vector<double>& result, Matrix<double>& tmp_strain, Matrix<double>& tmp_displacement);

  private:

    bool alsomatopt_;
  };

}

#endif /*SHAPEOPT_HH_*/
