#ifndef DUMAS_MMA_HH_
#define DUMAS_MMA_HH_

#include "Optimization/Optimizer/BaseOptimizer.hh"

class MMASolver;
class GCMMASolver;

namespace CoupledField
{
  /** this holds the C++ MMA and GCMMA implementation based on Jérémie Dumas' https://github.com/jdumas/mma
      The name is chosen to separate from our native MMA implementation. */
  class DumasMMA : public BaseOptimizer
  {
  public:
    /** @param dumas to be either DUMAS_MMA or DUMAS_GCMMA */
    DumasMMA(Optimization* optimization, PtrParamNode pn, Optimization::Optimizer dumas);
    
    virtual ~DumasMMA();

  protected:

    void PostInit() override;

    void SolveProblem() override;
    
    void ToInfo(PtrParamNode pn) override;

    void LogFileHeader(Optimization::Log& log) override;

    void LogFileLine(std::ofstream* out, PtrParamNode iteration)  override;

    void DescribeProperties(StdVector<std::pair<std::string, std::string> >& map) const override;

    void PythonSetProperty(pyObject* args) override;

  private:

    /** we have either mma or gcmma */
    MMASolver*   mma = nullptr; // -> MMASolver.h from Jérémie Dumas
    GCMMASolver* gcmma = nullptr; // -> GCMMASolver.h from Jérémie Dumas

    /** the number of gcmma search steps, to be logged */
    int gc_search = 0;

    Vector<double> xval;
    StdVector<double> dfdx;
    Vector<double> g;
    StdVector<double> dgdx;
    Vector<double> xmin;
    Vector<double> xmax;

    /** the current move limit, subject to be tuned */
    double move_limit = 0.5;

    double asyminit = 0.5;
    double asymdec = 0.7;
    double asyminc = 1.2;

    /** Svanberg mmasub.m but we have a, c, d as scalar in MMASolver
     %      Minimize  f_0(x) + a_0*z + sum( c_i*y_i + 0.5*d_i*(y_i)^2 )
     %    subject to  f_i(x) - a_i*z - y_i <= 0,  i = 1,...,m
     %                xmin_j <= x_j <= xmax_j,    j = 1,...,n
     %                z >= 0,   y_i >= 0,         i = 1,...,m */
    double ai_lin_aux_z = 0.0; // a0 is fixed as 1
    double c_lin_aux_y = 1000.0;
    double d_quad_aux_y = 0.0;

  };
} // end of name space

#endif /*DUMAS_MMA_HH_*/
