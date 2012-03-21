#ifndef _FEAS_SCP_HH_
#define _FEAS_SCP_HH_

#include "Optimization/Optimizer/scpip40i.hh"
#include "Optimization/Optimizer/BaseOptimizer.hh"

namespace CoupledField
{

/** This is based on ScpSolver.h (Sonja Ertel, MichaelStingel) and called ASCPSolver there */
class FeasSCP : public BaseOptimizer, public SCPIP
{
public:

  FeasSCP(Optimization* optimization, PtrParamNode pn);

	virtual ~FeasSCP();

  void SetInitialSolution(double* xinit); /* requires call to SetNumVars first! */
	void SetIntegerOption(int idx, int nVal); 
	void SetRealOption(int idx, double rVal); 
	void SetIntegerOptions(int* nVal); 
	void SetRealOptions(double* rVal); 

	/**  Access methods (for results) */
	inline int GetStatus() {return m_status;};
	double* GetSolution() ;

	/** implements the communication with the feasible scp from Sonja Lehmann, called here scpip40i
	 *  @see SCPIPBase::SolveProblem() */
	int SolveProblem(bool fromWarmstart);

	/**  Mutator methods */
  static bool IsMajor() { return (m_nFEvals % m_nMajor == 0) ;};
  static void Increment() { m_nFEvals++;};

private:

	static double m_UsrGlin (int l, double *y);
  static void m_dUsrGlin (int l, double *y, int* nnz, int* ind, double* val);

	/* solutions status */
	int m_status;
	int m_iemax;

	/* nvars (int)                                                                 */
	/* on entry: number of primal variables                                        */
	/* on exit : unchanged                                                         */
	int m_nvars;

	/* nobj (int)                                                                  */
	/* on entry: number of objectives                                     */
	/* on exit : unchanged                                                         */
	int m_nobj;

	/* nconstr (int)                                                               */
	/* on entry: number of all equality constraints (including linear)             */
	/* on exit : unchanged                                                         */
	int m_neq;

	/* nie (int)                                                                   */
	/* on entry: number of inequality constraints (including linear)               */
	/* on exit : unchanged                                                         */
	int m_nie;

	/* mf (int)                                                                   */
	/* on entry: number of feasibility constraints (including linear)               */
	/* on exit : unchanged                                                         */
	int m_mf;
	
	/* xinit (double array of size nvars)                                          */
	/* on entry: primal initial iterate                                            */
	/* on exit : approximate primal solution vector                                */
	/* comments: set flag ioptions[4] = 1 to ignore initial iterates               */
	double *m_xinit;

	/* lbv (double array of size nconstr)                                          */
	/* on entry: loweer bounds on scalar variables                                 */
	/* on exit : unchanged                                                         */
	/* comments: if constraint should not be bounded above,                        */
	/*           set corresponding value to HUGE_VAL                               */
	double *m_lbv;

	/* ubv (double array of size nconstr)                                          */
	/* on entry: upper bounds on scalar variables                                  */
	/* on exit : unchanged                                                         */
	/* comments: if constraint should not be bounded above,                        */
	/*           set corresponding value to HUGE_VAL                               */
	double *m_ubv;

  /** The design space */
  StdVector<double> x;

  /** lower bound of design space */
  StdVector<double> x_l;

  /** upper bound of design space */
  StdVector<double> x_u;


       
  /* temporary objective function value */
  double m_f_org;
           
  /* temporary function values of inequalities */
  double* m_h_org;    

  /* temporary function values of equalities */
  double* m_g_org;       

  /* temporary derivatives of objective */
  double* m_df;       

  /* temporary Lagrangian multipliers ineq. */
  double* m_y_ie;

  /* temporary Lagrangian multipliers eq. */
  double* m_y_eq;

  /* lower bounds ALL Lagrangian multipliers for lower bounds */
  double* m_y_l;

  /* upper bounds ALL Lagrangian multipliers for upper bounds  */
  double* m_y_u;

  /* doptions (double array of size 6)                                           */
	/* on entry: floating point options                                            */
	/* on exit : unchanged                                                         */
	/* 
	/*   option    | meaning                                       | default       */
	/* --------------------------------------------------------------------------- */

	double *m_doptions;

	/* ioptions (integer array of size 13)                                         */
	/* on entry: integer options                                                   */
	/* on exit : unchanged                                                         */
	/* 
	/*   option    | meaning                                       | default       */
	/* --------------------------------------------------------------------------- */

	int *m_ioptions;

	/* ninfo (integer array of size 23)                                         */
	int *m_ninfo;

	/* rinfo (double array of size 5)                                           */
	double *m_rinfo;
  
  /* output level */
  int m_nout;

  /* temporary double storage for SCPIP */
  double* m_r_scp;
      
  /* length of array above */
  int m_rdim;

  /* temporary double storage for SCPIP sub problem */
  double* m_r_sub;
      
  /* length of array above */
  int m_rsubdim;

  /* temporary integer storage for SCPIP */
  int* m_i_scp;
      
  /* length of array above */
  int m_idim;

  /* temporary integer storage for SCPIP sub problem */
  int* m_i_sub;
      
  /* length of array above */
  int m_isubdim;

  /* integer array indicating which gradients have to be evaluated */
  int* m_active; 
  int  m_mactiv;   /* number of active constraints */

  /* calling modus (fixed to 2=reverse comm.) */
  int m_mode; 

  /* Current Jacobian inequalities (sparse!) */
  int* m_iern;                 /* row indices */
  int* m_iecn;                 /* col indices */
  double* m_iederv;            /* values */
  int m_ielpar;                /* maximal length */ 
  int m_ieleng;                /* and current length */
  //int m_lengmie;                /* and current length */

  /* Current Hessian inequalities (sparse!) */
  //int* m_hcns;                /* row indices */
  //int* m_hrns;                /* col indices */
  //double* m_hval;             /* values */
  //int m_hleng;                /* and current length */
  //int m_hmax;                 /* maximal length */ 

  /* Current Jacobian equalities (sparse!) */
  int* m_eqrn;                 /* row indices */
  int* m_eqcn;                 /* col indices */
  double* m_eqcoef;            /* values */
  int m_eqlpar;                /* maximal length */ 
  int m_eqleng;                /* and current length */

  /* few more working arrays ... */
  int* m_spiw;                    
  int m_spiwdim;
  double* m_spdw;               
  int m_spdwdim;                      

  int m_linsys;
  int m_spstrat;
  int* m_linear;
  int m_lact;
  int* m_setact;


  // Internal scaling parameters
  double m_rScaleF;
  double m_rScaleE;
  bool m_bAutoScale;


  static int m_nFEvals;
  static int m_nMajor;

  const StdVector<int> matrix_sizes;
};

} // end of namespace
#endif // _FEAS_SCP_HH_
