#ifndef _FEAS_SCP_HH_
#define _FEAS_SCP_HH_

#include "Optimization/Optimizer/SCPIP.hh"

namespace CoupledField
{

/** This is based on ScpSolver.h (Sonja Ertel, MichaelStingel) and called ASCPSolver there */
class FeasSCP : public SCPIP
{
public:

  FeasSCP(Optimization* optimization, PtrParamNode pn, Optimization::Optimizer type = Optimization::FEAS_SCP_SOLVER);

	virtual ~FeasSCP();

	/** overwrites SCPIP::ToInfo() to add ipopt parameters */
	void ToInfo(PtrParamNode pn);

	/** implements the communication with the feasible scp from Sonja Lehmann, called here scpip40i
	 *  @see SCPIPBase::SolveProblem() */
	int solve_problem(bool fromWarmstart);

private:
	/** initializes additional parameters.
	 * @see SCPIPBase::AllocateProblem() */
	void AllocateProblem();

	/** identifies linear constraints. 0 for linear, 1 else. Is m_linear in ScpSolver.*/
  StdVector<int> linear;

  /** virtual cfs constraint indices of feasibility constraints */
  StdVector<int> mf_idx;

  /** Number of feasibility constraints (including linear). Is m_linear in ScpSolver. */
  int mf;

  /** Number of constraints that have to be active in each iterate. Is m_lact in ScpSolver. */
  int lact;

  /** Array of constraints which have to be active in each iteration. Is m_setact in ScpSolver. */
  StdVector<int> setact;

  /** Array of ipopt 2.3 subsolver option values */
  StdVector<double> ipopt_args;

  /** This contains the ipopt 2.3 parameter names in a primitive way at every 20th position. */
  char ipopt_cargs[200];
};

} // end of namespace
#endif // _FEAS_SCP_HH_
