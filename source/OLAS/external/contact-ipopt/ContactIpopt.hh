#ifndef IPOPTINTERFACE_HH_
#define IPOPTINTERFACE_HH_

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "MatVec/basematrix.hh"
#include "MatVec/matrix.hh"
#include "OLAS/solver/basesolver.hh"
#include "Utils/StdVector.hh"

#include <coin/IpTNLP.hpp>
#include <coin/IpIpoptApplication.hpp>

namespace CoupledField
{

class MathParser;
class CoordSystem;
class Grid;
class EqnMap;
class SurfElem;
template <class TYPE> class Vector;

using namespace Ipopt;

class ContactIpopt : public BaseIterativeSolver
{
public:
  /** @param pn here we can have options - might be NULL! */
  ContactIpopt(PtrParamNode param, PtrParamNode olasInfo, BaseMatrix::EntryType type);

  virtual ~ContactIpopt();

  void Setup(BaseMatrix &sysMat, PtrParamNode analysis_id);

  void Solve(const BaseMatrix &base_mat, const BasePrecond &base_precond, 
      const BaseVector &base_rhs,  BaseVector &base_sol, PtrParamNode analysis_id);

  /** Return what solver is used (here CholMod) */
  SolverType GetSolverType() { return CONTACTIPOPT; }

  void InitParameters();
  
  virtual void PrepareForAdjoint(BaseVector& sol);

protected:
  class TNLP : public Ipopt::TNLP {
  public:
    TNLP(ContactIpopt * const contactIpopt);

    virtual ~TNLP();

    /** Method to return some info about the nlp */
    bool get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
        Index& nnz_h_lag, IndexStyleEnum& index_style);

    /** Method to return the bounds for my problem */
    bool get_bounds_info(Index n, Number* x_l, Number* x_u,
        Index m, Number* g_l, Number* g_u);

    /** Method to return the starting point for the algorithm */
    bool get_starting_point(Index n, bool init_x, Number* x,
        bool init_z, Number* z_L, Number* z_U,
        Index m, bool init_lambda,
        Number* lambda);

    /** Method to return the objective value */
    bool eval_f(Index n, const Number* x, bool new_x, Number& obj_value);

    /** Method to return the gradient of the objective */
    bool eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f);

    /** Method to return the constraint residuals */
    bool eval_g(Index n, const Number* x, bool new_x, Index m, Number* g);

    /** Method to return:
     *   1) The structure of the jacobian (if "values" is NULL)
     *   2) The values of the jacobian (if "values" is not NULL)
     */
    bool eval_jac_g(Index n, const Number* x, bool new_x,
        Index m, Index nele_jac, Index* iRow, Index *jCol,
        Number* values);

    /** Method to return:
     *   1) The structure of the hessian of the lagrangian (if "values" is NULL)
     *   2) The values of the hessian of the lagrangian (if "values" is not NULL)
     */
    virtual bool eval_h(Index n, const Number* x, bool new_x,
        Number obj_factor, Index m, const Number* lambda,
        bool new_lambda, Index nele_hess, Index* iRow,
        Index* jCol, Number* values);

    /** This is called for each iteration and we use it to write the
     * current state to the result files */
    bool intermediate_callback(AlgorithmMode mode,
        Index iter, Number obj_value,
        Number inf_pr, Number inf_du,
        Number mu, Number d_norm,
        Number regularization_size,
        Number alpha_du, Number alpha_pr,
        Index ls_trials,
        const IpoptData* ip_data,
        IpoptCalculatedQuantities* ip_cq);

    /** This method is called when the algorithm is complete so the TNLP can store/write the solution */
    void finalize_solution(SolverReturn status,
        Index n, const Number* x, const Number* z_L, const Number* z_U,
        Index m, const Number* g, const Number* lambda,
        Number obj_value, const IpoptData* ip_data, IpoptCalculatedQuantities* ip_cq);

    /** this method allows to do constraint scaling. It is only called by IPOPT if
     *  "nlp_scaling_method" is set to "user-scaling". This is done by this class 
     * automatically if there is a constraint scaling provided. */  
    /*  bool get_scaling_parameters(Number& obj_scaling, 
                                 bool& use_x_scaling, Index n, Number* x_scaling,
                                 bool& use_g_scaling, Index m, Number* g_scaling);*/

  protected:
    ContactIpopt * const contactIpopt_;
    
    // storage for the indexes in the hessian for the corresponding constraints
    StdVector<Integer> validx;

  };

  struct ContactConstraintNode {
    UInt node;
    Integer solIdx[3];
  };

  struct ContactConstraint { // this exists per obstacle
    StdVector<ContactConstraintNode*> nodes;
    UInt mathParserHandleConstraint;
    UInt mathParserHandleGradient[3];
    UInt mathParserHandleHessian[3][3]; // only the lower triangle of this matrix is used
  };
  
  struct NodeElemContact {
    StdVector<UInt> nodeList;
    StdVector<SurfElem*> elemList;
    StdVector<UInt> hessianNodeNodeIdx;  // index for hessian
    StdVector<UInt> hessianEleEleIdx;    // index for hessian
    Matrix<UInt> hessianNodeEleIdx;      // index for hessian
  };
  
  struct BilateralContactConstraint {
    NodeElemContact contactRegions[2];
  };

  // Helper function to get node coordinate updated with solution
  void GetNodeCoordinate(CoupledField::Vector<Double>& Coord, const ContactConstraintNode* node, const Number* solution);
  
  // Helper function to get node coordinate updated with solution
  void GetNodeCoordWithSolution(Vector<Double>& Coord, UInt node, const Number* solution);
  
  // Check whether the a point lies left of the line from l1 to l2
  bool PointLiesLeftOfLine(Vector<Double>& p, Vector<Double>& l1, Vector<Double>& l2);
  
  double Distance(UInt node, StdVector<SurfElem*>& elemList, const Number* solution, Vector<Double>& GlobCoord, Vector<Double>& GlobCoord1, Vector<Double>& GlobCoord2, Vector<Double>& GlobCoordDiff, Vector<Double>& GlobCoordTmp, bool& segment, UInt& jmin);
  
  StdVector<ContactConstraint> contactConstraints;
  
  StdVector<BilateralContactConstraint> bilateralContactConstraints;
  
  /** set in prepareadjoint, so next solve solves an adjoint system */
  bool nextisadjoint;

  /** Base matrix */
  const BaseMatrix * sysMat_;

  /** RHS */
  const BaseVector * rhs_;

  /** Solution 
   * place the solution is stored to 
   * is also used to point to the forward solution corresponding to the next adjoint to be solved (if nextisadjoint) */
  BaseVector * sol_;
  
  /** Warmstart variables */
  StdVector<Double> x;
  StdVector<Double> lambda;
  bool initlambda;

  SmartPtr<IpoptApplication> app;

  SmartPtr<Ipopt::TNLP> tnlp;

  UInt numberofconstraints;
  UInt numberofjacobianentries;
  UInt numberofunilateraljacobianentries;

  // usefull shortcuts
  MathParser* mathParser;
  UInt dim;  
  CoordSystem* coosy;
  Grid* grid;
  EqnMap* eqnMap;
  
};


} // end of namespace
#endif /*IPOPTINTERFACE_HH_*/
