#include "DataInOut/Logging/cfslog.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/programOptions.hh"
#include "Domain/domain.hh"
#include "Domain/grid.hh"
#include "Domain/surfElem.hh"
#include "General/exception.hh"
#include "MatVec/scrs_matrix.hh"
#include "MatVec/vector.hh"
#include "OLAS/external/contact-ipopt/ContactIpopt.hh"
#include "PDE/eqnMap.hh"
#include "PDE/StdPDE.hh"
#include "Utils/tools.hh"

#include "coin/IpSolveStatistics.hpp"

using namespace CoupledField;

// declare class specific logging stream
DECLARE_LOG(ci)
DEFINE_LOG(ci, "contactIpopt")

/*
 * We interface the following problem to IpOpt
 * min 0.5 * u' * K * u
 * s.t. g(u) <= 0
 * where g(u) contains gap functions for obstacles
 * and "gap" functions of the type point to surface elements for bilateral contact
 */


ContactIpopt::ContactIpopt(PtrParamNode param, PtrParamNode olasInfo, BaseMatrix::EntryType type) {
  // we work with out
  xml_ =  param->Get("contact-ipopt", ParamNode::PASS);
  solverInfo_ = olasInfo->Get("contact-ipopt");

  if(type != BaseMatrix::DOUBLE){ // && type != BaseMatrix::COMPLEX && 
    EXCEPTION("unhandled type " << type);
  }
  
  app = new IpoptApplication();
  app->Initialize();
  tnlp = new TNLP(this);
  app->Options()->SetNumericValue("tol", 1e-10);
  // app->Options()->SetStringValue("linear_system_scaling", "none");
  InitParameters();
  mathParser = domain->GetMathParser();
  coosy = domain->GetCoordSystem();
  dim = domain->GetDim();
  numberofconstraints = 0;
  numberofjacobianentries = 0;
  numberofunilateraljacobianentries = 0;
  nextisadjoint = false;
  grid = domain->GetGrid();
  shared_ptr<EqnMap> eqnMapPtr = domain->GetStdPDE("mechanic")->GetEqnMap();
  eqnMap = eqnMapPtr.get();
  if(xml_->Has("obstacles")){
    ParamNodeList obstacles = xml_->Get("obstacles")->GetChildren();
    int nctct = obstacles.GetSize();
    contactConstraints.Resize(nctct);
    for(int ctct = 0; ctct < nctct; ctct++){
      PtrParamNode p = obstacles[ctct];
      ContactConstraint& constr = contactConstraints[ctct];
      if(p->GetName() == "sphere") {
        double obs_r = p->Get("radius")->As<Double>();
        double obs_x = p->Get("x")->As<Double>();
        double obs_y = p->Get("y")->As<Double>();
        double obs_z = p->Get("z")->As<Double>();
        constr.mathParserHandleConstraint = mathParser->GetNewHandle();
        mathParser->SetExpr(constr.mathParserHandleConstraint, "(x - obs_x)^2 + (y - obs_y)^2 + (z - obs_z)^2 - obs_r^2");
        mathParser->SetValue(constr.mathParserHandleConstraint, "obs_r", obs_r);
        mathParser->SetValue(constr.mathParserHandleConstraint, "obs_x", obs_x);
        mathParser->SetValue(constr.mathParserHandleConstraint, "obs_y", obs_y);
        mathParser->SetValue(constr.mathParserHandleConstraint, "obs_z", obs_z);
        mathParser->SetValue(constr.mathParserHandleConstraint, "z", 0.0); // we have to set z in case we have a 2D simulation
        constr.mathParserHandleGradient[0] = mathParser->GetNewHandle();
        mathParser->SetExpr(constr.mathParserHandleGradient[0], "2*(x - obs_x)");
        mathParser->SetValue(constr.mathParserHandleGradient[0], "obs_x", obs_x);
        constr.mathParserHandleGradient[1] = mathParser->GetNewHandle();
        mathParser->SetExpr(constr.mathParserHandleGradient[1], "2*(y - obs_y)");
        mathParser->SetValue(constr.mathParserHandleGradient[1], "obs_y", obs_y);
        constr.mathParserHandleGradient[2] = mathParser->GetNewHandle();
        mathParser->SetExpr(constr.mathParserHandleGradient[2], "2*(z - obs_z)");
        mathParser->SetValue(constr.mathParserHandleGradient[2], "obs_z", obs_z);
        mathParser->SetValue(constr.mathParserHandleGradient[2], "z", 0.0);
        for(int i = 0; i < 3; i++){
          for(int j = 0; j < 3; j++){
            if(i == j){
              constr.mathParserHandleHessian[i][j] = mathParser->GetNewHandle();
              mathParser->SetExpr(constr.mathParserHandleHessian[i][j], "2.0");
            }else{
              constr.mathParserHandleHessian[i][j] = 0;
            }
          }
        }
      }
      constr.nodes.Reserve(grid->GetNumNodes(p->Get("nodes")->As<std::string>()));
      shared_ptr<EntityList> el = grid->GetEntityList( EntityList::NODE_LIST, p->Get("nodes")->As<std::string>(), EntityList::NAMED_NODES); //THIS IS NEEDED: GetEntityList()->GetIterator() does not work, the list is immediately freed and the iterator walks nowhere 
      EntityIterator it = el->GetIterator();
      for(it.Begin(); !it.IsEnd(); it++) {
        ContactConstraintNode* node = new ContactConstraintNode;
        node->node = it.GetNode();
        for(UInt dof = 1; dof <= dim; dof++){
          const Integer eqnr = eqnMap->GetNodeEqn(node->node, dof) - 1;
          node->solIdx[dof-1] = eqnr; // equation numbers are 1 based
          if(eqnr >= 0){
            numberofjacobianentries++;
            numberofunilateraljacobianentries++;
          }
        }
        constr.nodes.Push_back(node);
        numberofconstraints++;
      }
    }
  }
  // bilateral part
  if(xml_->Has("bilateral")){
    ParamNodeList bilateral = xml_->Get("bilateral")->GetChildren();
    int nctct = bilateral.GetSize();
    bilateralContactConstraints.Resize(nctct);
    for(int ctct = 0; ctct < nctct; ctct++){
      PtrParamNode p = bilateral[ctct];
      BilateralContactConstraint& constr = bilateralContactConstraints[ctct];
      StdVector<RegionIdType> regId = StdVector<RegionIdType>();
      regId.Resize(2);
      regId[0] = grid->GetRegion().Parse(p->Get("region1")->As<std::string>());
      regId[1] = grid->GetRegion().Parse(p->Get("region2")->As<std::string>());
      for(int dir = 0; dir < 2; dir++){
        StdVector<UInt>& nodeList = constr.contactRegions[dir].nodeList;
        StdVector<SurfElem*>& elemList = constr.contactRegions[dir].elemList;      
        grid->GetNodesByRegion(nodeList, regId[dir]);
        grid->GetSurfElems(elemList, regId[1-dir]);
        int numnodes = nodeList.GetSize();
        int numelems = elemList.GetSize();
        numberofconstraints += numnodes;
        numberofjacobianentries += numnodes * dim * (1 + 2 * numelems); // see jacobian structure for this formula
      }
    }
  }
  lambda.Resize(numberofconstraints, 0.0);
  initlambda = false;
}

ContactIpopt::~ContactIpopt() {
  //TODO: release ContactConstraints & Nodes
}

void ContactIpopt::Setup(BaseMatrix &sysMat, PtrParamNode analysis_id) {
  // do we really want to create a new entry? Might blast up the output
  ParamNode::ActionType at = progOpts->DoDetailedInfo() ? ParamNode::APPEND : ParamNode::DEFAULT;
  PtrParamNode out = solverInfo_->Get(ParamNode::PROCESS)->Get("setup", at);
  out->Get("analysis_id")->SetValue(analysis_id->Get("analysis_id"));
  
  LOG_TRACE2(ci) <<  "Setup: matrix -> " << sysMat.ToString();
  sysMat_ = &sysMat;
  if(x.GetSize() != sysMat_->GetNumCols()){
    x.Resize(sysMat_->GetNumCols(), 0.0);
  }
}


void ContactIpopt::Solve(const BaseMatrix &base_mat, const BasePrecond &base_precond, 
    const BaseVector &base_rhs,  BaseVector &base_sol, PtrParamNode analysis_id)
{
  ParamNode::ActionType at = progOpts->DoDetailedInfo() ? ParamNode::APPEND : ParamNode::DEFAULT;
  PtrParamNode out = solverInfo_->Get(ParamNode::PROCESS)->Get("solver", at);
  out->Get("analysis_id")->SetValue(analysis_id->Get("analysis_id"));
  
  sysMat_ = &base_mat;
  rhs_ = &base_rhs;
  sol_ = &base_sol;
  
  app->Options()->SetStringValue("warm_start_init_point", initlambda ? "yes" : "no"); // if we have a forward problem, we reuse lambda
  ApplicationReturnStatus status = app->OptimizeTNLP(tnlp);

  nextisadjoint = false;

  if (status == Solve_Succeeded) {
    // Retrieve some statistics about the solve
    Index iter_count = app->Statistics()->IterationCount();
    Number final_obj = app->Statistics()->FinalObjective();
    
    std::cout << std::endl << "Problem solved in " << iter_count 
              << " iterations, final objective value is " << final_obj << std::endl; 
    
    out->Get("converged")->SetValue("yes");
    return;
  }

  switch(status)
  {
    case NonIpopt_Exception_Thrown:
      out->Get("converged")->SetValue("no");
      out->Get("reason/msg")->SetValue("non IPOPT exception occured");
      throw Exception("IPOPT stopped due to non-IPOPT exception. Try again with '-f'.");
         
    case Restoration_Failed:
      out->Get("converged")->SetValue("no");
      out->Get("reason/msg")->SetValue("IPOPT: 'Restoration failed'");
      throw Exception("IPOPT stopped with 'Restoration failed'");
         
    case Insufficient_Memory:
      out->Get("converged")->SetValue("no");
      out->Get("reason/msg")->SetValue("IPOPT: insufficient memory");
     throw Exception("IPOPT reports insufficient memory.");
         
    case Maximum_Iterations_Exceeded:
      out->Get("converged")->SetValue("no");
      out->Get("reason/msg")->SetValue("Maximum iterations exceeded");
      break;

    default:
      // positive is warning
      // Maximum_Iterations_Exceeded == -1
      if(status < Maximum_Iterations_Exceeded) {
        out->Get("converged")->SetValue("no");
        out->Get("reason/msg")->SetValue("IPOPT: error");
        out->Get("reason/error")->SetValue(status);
        EXCEPTION("IPOPT reported error " << status);
      }
      break;
      // else is no bad error and exits this void method :)  
  }
}

void ContactIpopt::InitParameters() {
  LOG_TRACE2(ci) <<  "InitParameters";

  // dump the parameter block and overwrite
  ParamNodeList options = xml_->GetList("option");
  for(UInt i = 0; i < options.GetSize(); i++) {
    PtrParamNode p = options[i];
    if(p->Get("type")->As<std::string>() == "integer"){
      app->Options()->SetIntegerValue(p->Get("name")->As<std::string>(), p->Get("value")->As<Integer>());
    }else if(p->Get("type")->As<std::string>() == "real"){
      app->Options()->SetNumericValue(p->Get("name")->As<std::string>(), p->Get("value")->As<Double>());
    }else{
      app->Options()->SetStringValue(p->Get("name")->As<std::string>(), p->Get("value")->As<std::string>());
    }
  }
}

void ContactIpopt::PrepareForAdjoint(BaseVector& sol){
  nextisadjoint = true;
  initlambda = false;
  sol_ = &sol;
  for(UInt i = 0; i < sol_->GetSize(); i++){
    sol_->GetEntry(i, x[i]);
  }
}


ContactIpopt::TNLP::TNLP(ContactIpopt * const contactIpopt) 
: contactIpopt_(contactIpopt)
{
  // empty method
}

ContactIpopt::TNLP::~TNLP(){  
}

bool ContactIpopt::TNLP::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g, Index& nnz_h_lag, IndexStyleEnum& index_style) {

  const SCRS_Matrix<Double>* scrs = dynamic_cast<const SCRS_Matrix<Double>*>(contactIpopt_->sysMat_);
  n = scrs->GetNumCols();
  
  m = contactIpopt_->numberofconstraints;
  nnz_jac_g = contactIpopt_->numberofjacobianentries;
  
  nnz_h_lag = scrs->GetNumEntries();
  index_style = C_STYLE;
  
  // the bilateral part contains for each bilateralcontactconstraint with each direction
  // #nodesInNodeList * 3 + #nodesInNodeList * #elemsInElemList * 2 * 2 * 2 + #elemsInElemList * (3+4+3) entries 
  for(UInt ctct = 0; ctct < contactIpopt_->bilateralContactConstraints.GetSize(); ctct++){
    BilateralContactConstraint& cct = contactIpopt_->bilateralContactConstraints[ctct];
    for(int dir = 0; dir < 2; dir++){
      int nrNode = cct.contactRegions[dir].nodeList.GetSize();
      int nrElem = cct.contactRegions[dir].elemList.GetSize();
      nnz_h_lag += 3 * nrNode + 8 * nrNode * nrElem + 10 * nrElem;
    }
  }
  
  LOG_TRACE2(ci) << ":get_nlp_info n <- " << n << "; m <- " << m << " nnz_jac_g <- " 
                 << nnz_jac_g << " nnz_h_lag <- " << nnz_h_lag << " index_style <- "
                 << (index_style == C_STYLE ? "C_STYLE" : "FORTRAN_STYLE");
  return true;
}

bool ContactIpopt::TNLP::get_bounds_info(Index n, Number* x_l, Number* x_u, Index m, Number* g_l, Number* g_u) {
  LOG_TRACE2(ci) << "get_bounds_info: n = " << n << "; m = " << m;
  double lbi, ubi;
  contactIpopt_->app->Options()->GetNumericValue("nlp_lower_bound_inf", lbi, "");
  contactIpopt_->app->Options()->GetNumericValue("nlp_upper_bound_inf", ubi, "");
  for(int i = 0; i < n; i++){
    x_l[i] = lbi;
    x_u[i] = ubi;
  }
  if(contactIpopt_->nextisadjoint){ // for the adjoint, active inequality constraints get equality constraints, inactive ones are ignored
    Vector<Double> g;
    g.Resize(m);
    eval_g(n, dynamic_cast<Vector<Double>*>(contactIpopt_->sol_)->GetPointer(), true, m, g.GetPointer() );
    for(int j = 0; j < m; j++){
      if(g[j] > -1e-6){
        g_l[j] = 0.0;
        g_u[j] = 0.0;
      }else{
        g_l[j] = lbi;
        g_u[j] = ubi;
      }
    }
  }else{
    for(int j = 0; j < m; j++){
      g_l[j] = 0;
      g_u[j] = ubi;
    }
  }
  return true;
}

bool ContactIpopt::TNLP::get_starting_point(Index n, bool init_x, Number* x,
                               bool init_z, Number* z_L, Number* z_U,
                               Index m, bool init_lambda,
                               Number* lambda)
{
  LOG_TRACE2(ci) << "get_starting_point: n = " << n << "; m = " << m;

  // Here, we assume we only have starting values for x, if you code
  // your own NLP, you can provide starting values for the others if
  // you wish.
  assert(init_x == true);
  
  assert((UInt)n == contactIpopt_->x.GetSize());
  for(int i = 0; i < n; i++){
    x[i] = contactIpopt_->x[i];
  }
  
  if(init_z){ // we have only unbound variables, however, ipopt, still wants values here
    for(int i = 0; i < n; i++){
      z_L[i] = 0.0;
      z_U[i] = 0.0;
    }
  }
  
  if(init_lambda){
    assert((UInt)m == contactIpopt_->lambda.GetSize());
    for(int i = 0; i < m; i++){
      lambda[i] = contactIpopt_->lambda[i];
    }
  }

  return true;
}

bool ContactIpopt::TNLP::eval_f(Index n, const Number* x, bool new_x, Number& obj_value) {
  CoupledField::Vector<Double> u;
  u.Replace(n, const_cast<Double*>(x), false);
  CoupledField::Vector<Double> Ku(n);
  contactIpopt_->sysMat_->Mult(u, Ku);
  obj_value = 0.5 * u.Inner(Ku) - u.Inner(*(dynamic_cast<const SingleVector*>(contactIpopt_->rhs_)));
  return true;
}

bool ContactIpopt::TNLP::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f) {
  CoupledField::Vector<Double> u;
  u.Replace(n, const_cast<Double*>(x), false);
  CoupledField::Vector<Double> g;
  g.Replace(n, grad_f, false);
  contactIpopt_->sysMat_->Mult(u, g); // g = Ku
  g.Add(-1.0, *(dynamic_cast<const SingleVector*>(contactIpopt_->rhs_))); // g -= right hand side
  return true;
}

void ContactIpopt::GetNodeCoordinate(Vector<Double>& Coord, const ContactConstraintNode * node, const Number* solution){
  grid->GetNodeCoordinate(Coord, node->node);
  for(UInt dof = 0; dof < dim; dof++){
    const Integer solIdx = node->solIdx[dof];
    if(solIdx >= 0){ // this might be wrong for inhomogeneous dirichlet points
      Coord[dof] += solution[solIdx];          
    }
  }
}

void ContactIpopt::GetNodeCoordWithSolution(Vector<Double>& Coord, UInt node, const Number* solution){
  grid->GetNodeCoordinate(Coord, node, true);
  if(solution != NULL){
    for(UInt d = 1; d <= dim; d++){
      int eqnr = eqnMap->GetNodeEqn(node, d) - 1;
      if(eqnr >= 0){
        Coord[d-1] += solution[eqnr];
      }
    }
  }
}

bool ContactIpopt::PointLiesLeftOfLine(Vector<Double>& p, Vector<Double>& l1, Vector<Double>& l2){
  // p lies left <==> the trianlge l1 l2 p is oriented math. positive (ccw) <==> the determinant |l2-l1 p-l1| is positive
  // this works just for 2d
  assert(dim == 2);
  return( (l2[0]-l1[0])*(p[1]-l1[1]) - (l2[1]-l1[1])*(p[0]-l1[0]) > 0);
}

double ContactIpopt::Distance(UInt node, StdVector<SurfElem*>& elemList, const Number* solution, Vector<Double>& GlobCoord, Vector<Double>& GlobCoord1, Vector<Double>& GlobCoord2, Vector<Double>& GlobCoordDiff, Vector<Double>& GlobCoordTmp, bool& segment, UInt& jmin){
  // this computes the squared euclidean norm as a distance!
  Double d = DBL_MAX; Double lambdamin = 0.0;
  GetNodeCoordWithSolution(GlobCoord, node, solution);
  for(UInt j = 0; j < elemList.GetSize(); j++){
    GetNodeCoordWithSolution(GlobCoord1, elemList[j]->connect[0], solution);
    GetNodeCoordWithSolution(GlobCoord2, elemList[j]->connect[1], solution);
    // we calculate lambda = (GlobCoord - GlobCoord1)*(GlobCoord2-GlobCoord1) / (GlobCoord2-GlobCoord1)^2
    // then (1-lambda)*ElemCoord[0] + lambda*ElemCoord[1] is the nearest point on this line
    GlobCoordDiff = GlobCoord2 - GlobCoord1;
    GlobCoordTmp = GlobCoord - GlobCoord1;
    Double lambda = GlobCoordDiff.Inner(GlobCoordTmp) / GlobCoordDiff.NormL2Squared();
    // for a line segment, projecting lambda to [0,1] yields the nearest point 
    Double lambdap = (lambda < 0.0) ? 0.0 : ( (lambda > 1.0) ? 1.0 : lambda);
    GlobCoordDiff = GlobCoord;
    GlobCoordDiff.Add(lambdap - 1.0, GlobCoord1);
    GlobCoordDiff.Add(-lambdap, GlobCoord2);
    Double td = GlobCoordDiff.NormL2Squared();
    if(td < d || (td == d && lambda < lambdamin)){ // the latter part ensures, that if we are closest to a node, we always find the segment starting with this node
      d = td;
      jmin = j;
      lambdamin = lambda;
    }
  }
  // we have found a nearest segment/node, we have to check, whether we lie inside or outside
  segment = lambdamin > 0.0 && lambdamin < 1.0;
  if(segment){ // we have a segment
    GetNodeCoordWithSolution(GlobCoord1, elemList[jmin]->connect[0], solution);
    GetNodeCoordWithSolution(GlobCoord2, elemList[jmin]->connect[1], solution);
    if(PointLiesLeftOfLine(GlobCoord, GlobCoord1, GlobCoord2)){ // point lies left, i.e. inside
      return(-d);
    }else{
      return(d);
    }
  }else{
    // the situation here is more difficult, we have only a node, we have to consider both segments:
    // if both segments together make a left turn at this node, then the point lies inside, if it lies left of both segments
    // if both segments together make a right turn at this node, then the point lies inside, if it lies left of at least one of both segments
    // first, we have to find the other segment
    UInt searchnode = elemList[jmin]->connect[0];
    UInt foundnode = searchnode;
    for(UInt j = 0; j < elemList.GetSize(); j++){
      if(elemList[j]->connect[1] == searchnode){
        foundnode = elemList[j]->connect[0];
        break;
      }
    }
    assert(foundnode != searchnode);
    GetNodeCoordWithSolution(GlobCoordDiff, searchnode, solution); // the middle node          
    GetNodeCoordWithSolution(GlobCoord1, foundnode, solution);
    GetNodeCoordWithSolution(GlobCoord2, elemList[jmin]->connect[1], solution);
    if(PointLiesLeftOfLine(GlobCoordDiff, GlobCoord1, GlobCoord2)) { // we make a right turn
      if(PointLiesLeftOfLine(GlobCoord, GlobCoord1, GlobCoordDiff) || PointLiesLeftOfLine(GlobCoord, GlobCoordDiff, GlobCoord2)){
        return(-d);
      }else{
        return(d);
      }
    }else{
      if(PointLiesLeftOfLine(GlobCoord, GlobCoord1, GlobCoordDiff) && PointLiesLeftOfLine(GlobCoord, GlobCoordDiff, GlobCoord2)){
        return(-d);
      }else{
        return(d);
      }
    }
  }
 
}

bool ContactIpopt::TNLP::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g) {
  Vector<Double> GlobCoord, GlobCoord1, GlobCoord2, GlobCoordDiff, GlobCoordTmp;
  bool dummy1; UInt dummy2;
  
  UInt c = 0;
  for(UInt ctct = 0; ctct < contactIpopt_->contactConstraints.GetSize(); ctct++){
    ContactConstraint& cct = contactIpopt_->contactConstraints[ctct];
    for(UInt i = 0; i < cct.nodes.GetSize(); i++){
      ContactConstraintNode* node = cct.nodes[i];
      contactIpopt_->GetNodeCoordinate(GlobCoord, node, x);
      contactIpopt_->mathParser->SetCoordinates(cct.mathParserHandleConstraint, *(contactIpopt_->coosy), GlobCoord);
      g[c++] = contactIpopt_->mathParser->Eval(cct.mathParserHandleConstraint);
    }    
  }
  // bilateral part
  for(UInt ctct = 0; ctct < contactIpopt_->bilateralContactConstraints.GetSize(); ctct++){
    BilateralContactConstraint& cct = contactIpopt_->bilateralContactConstraints[ctct];
    for(int dir = 0; dir < 2; dir++){ // border1 to 2 and border2 to 1
      StdVector<UInt>& nodeList = cct.contactRegions[dir].nodeList;
      StdVector<SurfElem*>& elemList = cct.contactRegions[dir].elemList;
      // for each node on one side of the border, find the closest segment / node
      for(UInt i = 0; i < nodeList.GetSize(); i++){
        double v = contactIpopt_->Distance(nodeList[i], elemList, x, GlobCoord, GlobCoord1, GlobCoord2, GlobCoordDiff, GlobCoordTmp, dummy1, dummy2); 
        LOG_DBG3(ci) << "constraint " << c << " value: " << v;
        g[c++] = v;
      }
    }    
  }
  assert(c == (UInt)m);
  return true;
}

bool ContactIpopt::TNLP::eval_jac_g(Index n, const Number* x, bool new_x,
                       Index m, Index nele_jac, Index* iRow, Index *jCol,
                       Number* values)
{
  Vector<Double> GlobCoord, GlobCoord1, GlobCoord2, GlobCoordDiff, GlobCoordTmp;
  
  UInt dim = contactIpopt_->dim;
  
  UInt c = 0;
  if (values == NULL) { // only structure
    UInt j = 0;
    for(UInt ctct = 0; ctct < contactIpopt_->contactConstraints.GetSize(); ctct++){
      ContactConstraint cct = contactIpopt_->contactConstraints[ctct];
      for(UInt i = 0; i < cct.nodes.GetSize(); i++){
        ContactConstraintNode* node = cct.nodes[i];
        for(UInt dof = 0; dof < contactIpopt_->dim; dof++){
          Integer solIdx = node->solIdx[dof];
          if(solIdx >= 0){
            LOG_DBG3(ci) << "Jacobian structure: index " << c << " at position " << j << ", " << solIdx;
            iRow[c] = j;
            jCol[c++] = solIdx;
          }
        }
        j++;
      }
    }
    // The jacobian of each constraint has entries corresponding to all nodes of all elements of the element list and to itself.
    for(UInt ctct = 0; ctct < contactIpopt_->bilateralContactConstraints.GetSize(); ctct++){
      BilateralContactConstraint& cct = contactIpopt_->bilateralContactConstraints[ctct];
      for(int dir = 0; dir < 2; dir++){
        StdVector<UInt>& nodeList = cct.contactRegions[dir].nodeList;
        StdVector<SurfElem*>& elemList = cct.contactRegions[dir].elemList;
        for(UInt i = 0; i < nodeList.GetSize(); i++){ 
          // every node here corresponds to one single constraint. it can have entries at itself and at every other node of the corresponding element list
          // we store them per element however, this makes the size double but is much more simple
          // each constraint has then dim * ( 1 + 2 * elemList.GetSize()) possible non-zeros
          for(UInt d = 1; d <= dim; d++){
            int eqnr = contactIpopt_->eqnMap->GetNodeEqn(nodeList[i], d) - 1;
            LOG_DBG3(ci) << "Jacobian structure: index " << c << " at position " << j << ", " << eqnr;
            iRow[c] = j;
            jCol[c++] = eqnr;
          }
          for(UInt ei = 0; ei < elemList.GetSize(); ei++){
            const StdVector<UInt>& nodes = elemList[ei]->connect;
            for(UInt ni = 0; ni < 2; ni++){ // only the first two nodes are used (quadratic elements might not work)
              for(UInt d = 1; d <= dim; d++){
                int eqnr = contactIpopt_->eqnMap->GetNodeEqn(nodes[ni], d) - 1;
                LOG_DBG3(ci) << "Jacobian structure: index " << c << " at position " << j << ", " << eqnr;
                iRow[c] = j;
                jCol[c++] = eqnr;
              }
            }
          }
          j++;
        }
      }
    }
    assert(j == (UInt)m);
  }else{ // do evaluation
    
    // first reset all values as bilateral part does not always set every value
    for(Integer i = 0; i < nele_jac; i++){
      values[i] = 0;
    }
    
    for(UInt ctct = 0; ctct < contactIpopt_->contactConstraints.GetSize(); ctct++){
      ContactConstraint cct = contactIpopt_->contactConstraints[ctct];
      for(UInt i = 0; i < cct.nodes.GetSize(); i++){
        ContactConstraintNode* node = cct.nodes[i];
        contactIpopt_->GetNodeCoordinate(GlobCoord, node, x);
        for(UInt dof = 0; dof < contactIpopt_->dim; dof++){
          contactIpopt_->mathParser->SetCoordinates(cct.mathParserHandleGradient[dof], *(contactIpopt_->coosy), GlobCoord);
        }
        for(UInt dof = 0; dof < contactIpopt_->dim; dof++){
          Integer solIdx = node->solIdx[dof];
          if(solIdx >= 0){
            double v = contactIpopt_->mathParser->Eval(cct.mathParserHandleGradient[dof]);
            LOG_DBG3(ci) << "Jacobian values: index " << c << ", value: " << v;
            values[c++] = v;
          }
        }
      }
    }
    // bilateral part
    for(UInt ctct = 0; ctct < contactIpopt_->bilateralContactConstraints.GetSize(); ctct++){
      BilateralContactConstraint& cct = contactIpopt_->bilateralContactConstraints[ctct];
      for(int dir = 0; dir < 2; dir++){ // border1 to 2 and border2 to 1
        StdVector<UInt>& nodeList = cct.contactRegions[dir].nodeList;
        StdVector<SurfElem*>& elemList = cct.contactRegions[dir].elemList;
        for(UInt i = 0; i < nodeList.GetSize(); i++){
          bool segment;
          UInt jmin;
          double dist = contactIpopt_->Distance(nodeList[i], elemList, x, GlobCoord, GlobCoord1, GlobCoord2, GlobCoordDiff, GlobCoordTmp, segment, jmin);
          bool inside = dist < 0;
          double sig = inside ? -1.0 : 1.0;
          if(segment){ // we currently are nearest to a segment, we use the following distance formula for derivation
            // det(v-p1,p2-p1)^2 / |p2-p1|^2 = n / d
            // we have the coordinates of v = GlobCoord, p1 = GlobCoord1, p2 = GlobCoord2
            double t = (GlobCoord[0] - GlobCoord1[0])*(GlobCoord2[1] - GlobCoord1[1]) - (GlobCoord[1] - GlobCoord1[1])*(GlobCoord2[0] - GlobCoord1[0]); // the determinant
            double tdGlobCoord0 = (GlobCoord2[1] - GlobCoord1[1]);
            double tdGlobCoord1 = (GlobCoord1[0] - GlobCoord2[0]);
            double tdGlobCoord10 = (GlobCoord[1] - GlobCoord2[1]);
            double tdGlobCoord11 = (GlobCoord2[0] - GlobCoord[0]);
            double tdGlobCoord20 = (GlobCoord1[1] - GlobCoord[1]);
            double tdGlobCoord21 = (GlobCoord[0] - GlobCoord1[0]);            
            double n = t*t;
            double ndGlobCoord0 = 2.0 * t * tdGlobCoord0;
            double ndGlobCoord1 = 2.0 * t * tdGlobCoord1;
            double ndGlobCoord10 = 2.0 * t * tdGlobCoord10;
            double ndGlobCoord11 = 2.0 * t * tdGlobCoord11;
            double ndGlobCoord20 = 2.0 * t * tdGlobCoord20;
            double ndGlobCoord21 = 2.0 * t * tdGlobCoord21;
            double d = (GlobCoord1[0] - GlobCoord2[0])*(GlobCoord1[0] - GlobCoord2[0]) + (GlobCoord1[1] - GlobCoord2[1])*(GlobCoord1[1] - GlobCoord2[1]);
            double ddGlobCoord10 = 2.0 * (GlobCoord1[0] - GlobCoord2[0]);
            double ddGlobCoord11 = 2.0 * (GlobCoord1[1] - GlobCoord2[1]);
            double ddGlobCoord20 = 2.0 * (GlobCoord2[0] - GlobCoord1[0]);
            double ddGlobCoord21 = 2.0 * (GlobCoord2[1] - GlobCoord1[1]);
            // sig * n / d
            values[c                   ] = sig * ndGlobCoord0 / d; // dGlobCoord[0]
            values[c+1                 ] = sig * ndGlobCoord1 / d; // dGlobCoord[1]
            values[c+dim*(1+2*jmin)    ] = sig * ( d * ndGlobCoord10 - n * ddGlobCoord10 ) / (d*d); // dGlobCoord1[0]
            values[c+dim*(1+2*jmin)+1  ] = sig * ( d * ndGlobCoord11 - n * ddGlobCoord11 ) / (d*d); // dGlobCoord1[1]
            values[c+dim*(1+2*jmin+1)  ] = sig * ( d * ndGlobCoord20 - n * ddGlobCoord20 ) / (d*d); // dGlobCoord2[0]
            values[c+dim*(1+2*jmin+1)+1] = sig * ( d * ndGlobCoord21 - n * ddGlobCoord21 ) / (d*d); // dGlobCoord2[1]           
            LOG_DBG3(ci) << "Jacobian values: index " << (c                   ) << ", value: " << (sig * ndGlobCoord0 / d); // dGlobCoord[0]
            LOG_DBG3(ci) << "Jacobian values: index " << (c+1                 ) << ", value: " << (sig * ndGlobCoord1 / d); // dGlobCoord[1]
            LOG_DBG3(ci) << "Jacobian values: index " << (c+dim*(1+2*jmin)    ) << ", value: " << (sig * ( d * ndGlobCoord10 - n * ddGlobCoord10 ) / (d*d)); // dGlobCoord1[0]
            LOG_DBG3(ci) << "Jacobian values: index " << (c+dim*(1+2*jmin)+1  ) << ", value: " << (sig * ( d * ndGlobCoord11 - n * ddGlobCoord11 ) / (d*d)); // dGlobCoord1[1]
            LOG_DBG3(ci) << "Jacobian values: index " << (c+dim*(1+2*jmin+1)  ) << ", value: " << (sig * ( d * ndGlobCoord20 - n * ddGlobCoord20 ) / (d*d)); // dGlobCoord2[0]
            LOG_DBG3(ci) << "Jacobian values: index " << (c+dim*(1+2*jmin+1)+1) << ", value: " << (sig * ( d * ndGlobCoord21 - n * ddGlobCoord21 ) / (d*d)); // dGlobCoord2[1]           
          }else{ // we are nearest to a point, the formula is much simpler
            // |v-p|^2, with v = GlobCoord and p = GlobCoordDiff
            // (GlobCoord[0] - GlobCoordDiff[0])*(GlobCoord[0] - GlobCoordDiff[0]) + (GlobCoord[1] - GlobCoordDiff[1])*(GlobCoord[1] - GlobCoordDiff[1]);
            values[c                 ] = sig * 2.0 * (GlobCoord[0] - GlobCoordDiff[0]);     // dGlobCoord[0]
            values[c+1               ] = sig * 2.0 * (GlobCoord[1] - GlobCoordDiff[1]);     // dGlobCoord[1]
            values[c+dim*(1+2*jmin)  ] = sig * 2.0 * (GlobCoordDiff[0] - GlobCoord[0]); // dGlobCoordDiff[0]
            values[c+dim*(1+2*jmin)+1] = sig * 2.0 * (GlobCoordDiff[1] - GlobCoord[1]); // dGlobCoordDiff[1]
            LOG_DBG3(ci) << "Jacobian values: index " << (c                 ) <<  ", value: " << (sig * 2.0 * (GlobCoord[0] - GlobCoordDiff[0]));     // dGlobCoord[0]
            LOG_DBG3(ci) << "Jacobian values: index " << (c+1               ) <<  ", value: " << (sig * 2.0 * (GlobCoord[1] - GlobCoordDiff[1]));     // dGlobCoord[1]
            LOG_DBG3(ci) << "Jacobian values: index " << (c+dim*(1+2*jmin)  ) <<  ", value: " << (sig * 2.0 * (GlobCoordDiff[0] - GlobCoord[0])); // dGlobCoordDiff[0]
            LOG_DBG3(ci) << "Jacobian values: index " << (c+dim*(1+2*jmin)+1) <<  ", value: " << (sig * 2.0 * (GlobCoordDiff[1] - GlobCoord[1])); // dGlobCoordDiff[1]
          }
          c += dim * (1 + 2 * elemList.GetSize()); // advance the width of the jacobian block
        }
      }    
    }
  }
  assert(c == (UInt)nele_jac);
  return true;
}

bool ContactIpopt::TNLP::eval_h(Index n, const Number* x, bool new_x,
    Number obj_factor, Index m, const Number* lambda,
    bool new_lambda, Index nele_hess, Index* iRow,
    Index* jCol, Number* values)
{
  Vector<Double> GlobCoord, GlobCoord1, GlobCoord2, GlobCoordDiff, GlobCoordTmp;
  
  UInt dim = contactIpopt_->dim;

  const SCRS_Matrix<Double>* scrs = dynamic_cast<const SCRS_Matrix<Double>*>(contactIpopt_->sysMat_);
  if(values == NULL){ // we need the structure of the hessian
    const UInt* rows = scrs->GetRowPointer();
    const UInt* cols = scrs->GetColPointer();
    UInt idx = 0;
    UInt nrows = scrs->GetNumRows();
    for(UInt r = 0; r < nrows; r++){
      UInt mcol = (r == nrows-1) ? scrs->GetNumEntries() : rows[r+1];
      for(UInt c = rows[r]; c < mcol; c++){
        LOG_DBG3(ci) << "Hessian structure: index " << idx << " at position " << cols[c] << ", " << r; 
        iRow[idx] = cols[c];
        jCol[idx] = r;
        idx++;
      }
    }
    // the hessian is split into two different parts, which apart from the diagonal is disjoint for (unilateral & system) and bilateral contacts constraints
    // thus, we built the two blocks seperately
    // note that for ipopt, the entries to the matrix do not have to be sorted, just have to be in the same sequence for structure and actual hessian
    // first for the system matrix and unilateral part:
    // we prepare the structure of the constraint hessian, the structure is rather simple, we may have entries on the diagonal and up to 2 cols away from it
    // note that these entries belong to one single node, so they are certainly already found in the system matrix
    validx.Resize(0);
    validx.Reserve(contactIpopt_->numberofunilateraljacobianentries*dim*(dim+1)/2);
    // the array validx stores the indexes to the array values, when walking over the constraints in the following fashion
    for(UInt ctct = 0; ctct < contactIpopt_->contactConstraints.GetSize(); ctct++){
      ContactConstraint cct = contactIpopt_->contactConstraints[ctct];
      for(UInt i = 0; i < cct.nodes.GetSize(); i++){
        ContactConstraintNode* node = cct.nodes[i];
        for(UInt dof = 0; dof < dim; dof++){
          Integer solIdx = node->solIdx[dof];
          if(solIdx >= 0){
            for(UInt dof2 = dof; dof2 < dim; dof2++){
              if(cct.mathParserHandleHessian[dof][dof2] > 0){
                Integer solIdx2 = node->solIdx[dof2];
                if(solIdx2 >= 0){
                  // we found a possible entry at (solIdx,solIdx2), we need to find the right index in scrs for it, note that we only may be (dof2-dof) columns away
                  // we usually are dof2-dof columns away, but if one dof is missing, this might be lower
                  // note that solIdx <= solIdx2 as dof <= dof2
                  UInt r = rows[solIdx];
#ifndef NDEBUG
                  bool found = false;
#endif
                  for(UInt c = r; c <= r+dof2-dof; c++){
                    if(cols[c] == (UInt)solIdx2){
                      validx.Push_back(c);
#ifndef NDEBUG
                      found = true;
#endif
                      break;
                    }
                  }
#ifndef NDEBUG
                  assert(found);
#endif
                }
              }
            }
          }
        }
      }
    }
    // bilateral part
    // each constraint consists actually of three different contributions to the hessian
    // 1: node to itself, 2: node to elem, 3: elem to itself
    // parts 1 & 3 can already be found in the system matrix, however, we add them again, this can be done with ipopt
    // a search comparable to the one above could be used instead, feel free to change this
    // part 1 is disjoint for each constraint
    // part 2 also is disjoint for each constraint
    // part 3 however, is not, every constraint belonging to one BilateralContactConstraint shares the same elements, and thus the same block
    // as we have no object for each individual constraint, we store the indexes in BilateralContactConstraints
    // part 1 in hessianNodeNodeIdx[#nodesInNodeList]
    // part 2 in hessianNodeEleIdx[#nodesInNodeList][#elemsInElemList]
    // part 3 in hessianEleEleIdx[#elemsInElemList]
    // see value part for more information about the single hessian entries 
    // the bilateral part contains for each bilateralcontactconstraint with each direction
    // #nodesInNodeList * 3 + #nodesInNodeList * #elemsInElemList * 2 * 2 * 2 + #elemsInElemList * (3+4+3) entries 
    for(UInt ctct = 0; ctct < contactIpopt_->bilateralContactConstraints.GetSize(); ctct++){
      BilateralContactConstraint& cct = contactIpopt_->bilateralContactConstraints[ctct];
      for(int dir = 0; dir < 2; dir++){
        StdVector<UInt>& nodeList = cct.contactRegions[dir].nodeList;
        StdVector<SurfElem*>& elemList = cct.contactRegions[dir].elemList;
        StdVector<UInt>& hessianNodeNodeIdx = cct.contactRegions[dir].hessianNodeNodeIdx;
        StdVector<UInt>& hessianEleEleIdx = cct.contactRegions[dir].hessianEleEleIdx;
        Matrix<UInt>& hessianNodeEleIdx = cct.contactRegions[dir].hessianNodeEleIdx;
        hessianNodeNodeIdx.Reserve(nodeList.GetSize());
        // part1, for each node in the nodelist, we have to find its equationnumbers
        for(UInt i = 0; i < nodeList.GetSize(); i++){
          hessianNodeNodeIdx.Push_back(idx);
          for(UInt dof1 = 1; dof1 <= dim; dof1++){
            int eqnr1 = contactIpopt_->eqnMap->GetNodeEqn(nodeList[i], dof1) - 1;
            for(UInt dof2 = dof1; dof2 <= dim; dof2++){
              int eqnr2 = contactIpopt_->eqnMap->GetNodeEqn(nodeList[i], dof2) - 1;
              LOG_DBG3(ci) << "Hessian structure: index " << idx << " at position " << eqnr1 << ", " << eqnr2; 
              iRow[idx] = eqnr1;
              jCol[idx] = eqnr2;
              idx++;
            }
          }
        }
        // part2, for each node in the nodeList and each elem in the elemList, we have to find the corresponding equation numbers
        hessianNodeEleIdx.Resize(nodeList.GetSize(), elemList.GetSize());
        for(UInt i = 0; i < nodeList.GetSize(); i++){
          for(UInt j = 0; j < elemList.GetSize(); j++){
            const StdVector<UInt>& nodes = elemList[j]->connect;
            hessianNodeEleIdx[i][j] = idx;
            for(UInt node = 0; node < 2; node++){ // node of the elem
              for(UInt dof1 = 1; dof1 <= dim; dof1++){
                int eqnr1 = contactIpopt_->eqnMap->GetNodeEqn(nodeList[i], dof1) - 1;
                for(UInt dof2 = 1; dof2 <= dim; dof2++){
                  int eqnr2 = contactIpopt_->eqnMap->GetNodeEqn(nodes[node], dof2) - 1;
                  LOG_DBG3(ci) << "Hessian structure: index " << idx << " at position " << eqnr1 << ", " << eqnr2; 
                  iRow[idx] = eqnr1;
                  jCol[idx] = eqnr2;
                  idx++;
                }
              }
            }
          }
        }
        // part3, for each elem in the elemList, we have to find the corresponding equation numbers, this is a little more difficult
        hessianEleEleIdx.Reserve(elemList.GetSize());
        for(UInt i = 0; i < elemList.GetSize(); i++){
          hessianEleEleIdx.Push_back(idx);
          const StdVector<UInt>& nodes = elemList[i]->connect;
          for(UInt node1 = 0; node1 < 2; node1++){
            for(UInt node2 = node1; node2 < 2; node2++){
              for(UInt dof1 = 1; dof1 <= dim; dof1++){
                int eqnr1 = contactIpopt_->eqnMap->GetNodeEqn(nodes[node1], dof1) - 1;
                for(UInt dof2 = (node2 == node1 ? dof1 : 1); dof2 <= dim; dof2++){
                  int eqnr2 = contactIpopt_->eqnMap->GetNodeEqn(nodes[node2], dof2) - 1;
                  LOG_DBG3(ci) << "Hessian structure: index " << idx << " at position " << eqnr1 << ", " << eqnr2; 
                  iRow[idx] = eqnr1;
                  jCol[idx] = eqnr2;
                  idx++;
                }
              }
            }
          }
        }
        // end of all 3 parts of bilateral constraint
      } // direction      
    } // bilateral Contact Constraint
    assert(idx == (UInt)nele_hess);
  }else{ // we need the actual hessian
    
    // first reset all values as bilateral part does not always set every value
    for(Integer i = 0; i < nele_hess; i++){
      values[i] = 0;
    }
    
    Vector<Double> GlobCoord;
    const Double* vals = scrs->GetDataPointer();
    for(UInt i = 0; i < scrs->GetNumEntries(); i++){
      values[i] = obj_factor*vals[i];
    }
    // the hessian of the constraints
    UInt cidx = 0; // index of the constraint
    UInt idx = 0; // index in validx
    for(UInt ctct = 0; ctct < contactIpopt_->contactConstraints.GetSize(); ctct++){
      ContactConstraint cct = contactIpopt_->contactConstraints[ctct];
      for(UInt i = 0; i < cct.nodes.GetSize(); i++){
        ContactConstraintNode* node = cct.nodes[i];
        contactIpopt_->GetNodeCoordinate(GlobCoord, node, x);
        for(UInt dof = 0; dof < contactIpopt_->dim; dof++){
          Integer solIdx = node->solIdx[dof];
          if(solIdx >= 0){
            for(UInt dof2 = dof; dof2 < contactIpopt_->dim; dof2++){
              if(cct.mathParserHandleHessian[dof][dof2] > 0){
                contactIpopt_->mathParser->SetCoordinates(cct.mathParserHandleHessian[dof][dof2], *(contactIpopt_->coosy), GlobCoord);
                Integer solIdx2 = node->solIdx[dof2];
                if(solIdx2 >= 0){
                  Double value = contactIpopt_->mathParser->Eval(cct.mathParserHandleHessian[dof][dof2]);
                  values[validx[idx++]] += lambda[cidx] * value;
                }
              }
            }
          }
        }
        cidx++;
      }
    }
    assert(idx == validx.GetSize());
    // bilateral part
    for(UInt ctct = 0; ctct < contactIpopt_->bilateralContactConstraints.GetSize(); ctct++){
      BilateralContactConstraint& cct = contactIpopt_->bilateralContactConstraints[ctct];
      for(int dir = 0; dir < 2; dir++){ // border1 to 2 and border2 to 1
        StdVector<UInt>& nodeList = cct.contactRegions[dir].nodeList;
        StdVector<SurfElem*>& elemList = cct.contactRegions[dir].elemList;
        Matrix<UInt>& hessianNodeEleIdx = cct.contactRegions[dir].hessianNodeEleIdx;
        StdVector<UInt>& hessianNodeNodeIdx = cct.contactRegions[dir].hessianNodeNodeIdx;
        StdVector<UInt>& hessianEleEleIdx = cct.contactRegions[dir].hessianEleEleIdx;
        for(UInt i = 0; i < nodeList.GetSize(); i++){
          UInt hNN = hessianNodeNodeIdx[i]; // contains only xx xy and yy (3 hessian entries)
          assert(dim == 2);
          bool segment;
          UInt jmin;
          double dist = contactIpopt_->Distance(nodeList[i], elemList, x, GlobCoord, GlobCoord1, GlobCoord2, GlobCoordDiff, GlobCoordTmp, segment, jmin);
          bool inside = dist < 0;
          LOG_DBG3(ci) << "Hessian values for constraint " << cidx << " with factor " << lambda[cidx];
          double fac = lambda[cidx++] * (inside ? -1.0 : 1.0);
          if(segment){ // we currently are nearest to a segment, we use the following distance formula for derivation
            // det(v-p1,p2-p1)^2 / |p2-p1|^2 = n / d
            // we have the coordinates of v = GlobCoord, p1 = GlobCoord1, p2 = GlobCoord2
            UInt hNE1 = hessianNodeEleIdx[i][jmin]; // this contains: node-node1: xx xy yx yy, node-node2: xx xy yx yy (8 hessian entries)
            UInt hNE2 = hNE1 + 4; // shortcut to node-node2
            UInt hE1E1 = hessianEleEleIdx[jmin]; // contains 10 hessian entries first node1-node1 xx xy yy
            UInt hE1E2 = hE1E1 + 3; // then node1-node2 xx xy yx yy
            UInt hE2E2 = hE1E2 + 4; // then node2-node2 xx yy xy
            double t = (GlobCoord[0] - GlobCoord1[0])*(GlobCoord2[1] - GlobCoord1[1]) - (GlobCoord[1] - GlobCoord1[1])*(GlobCoord2[0] - GlobCoord1[0]);
            double tdGlobCoord0 = (GlobCoord2[1] - GlobCoord1[1]);
            double tdGlobCoord1 = (GlobCoord1[0] - GlobCoord2[0]);
            double tdGlobCoord10 = (GlobCoord[1] - GlobCoord2[1]);
            double tdGlobCoord11 = (GlobCoord2[0] - GlobCoord[0]);
            double tdGlobCoord20 = (GlobCoord1[1] - GlobCoord[1]);
            double tdGlobCoord21 = (GlobCoord[0] - GlobCoord1[0]);            
            double n = t*t;
            double ndGlobCoord10 = 2.0 * t * tdGlobCoord10;
            double ndGlobCoord11 = 2.0 * t * tdGlobCoord11;
            double ndGlobCoord20 = 2.0 * t * tdGlobCoord20;
            double ndGlobCoord21 = 2.0 * t * tdGlobCoord21;
            double d = (GlobCoord1[0] - GlobCoord2[0])*(GlobCoord1[0] - GlobCoord2[0]) + (GlobCoord1[1] - GlobCoord2[1])*(GlobCoord1[1] - GlobCoord2[1]);
            double ddGlobCoord10 = 2.0 * (GlobCoord1[0] - GlobCoord2[0]);
            double ddGlobCoord11 = 2.0 * (GlobCoord1[1] - GlobCoord2[1]);
            double ddGlobCoord20 = 2.0 * (GlobCoord2[0] - GlobCoord1[0]);
            double ddGlobCoord21 = 2.0 * (GlobCoord2[1] - GlobCoord1[1]);
            // sig * n / d
            // dGlobCoord[0] -> sig * ndGlobCoord0 / d = sig * 2.0 * t * tdGlobCoord0 / d
            values[hNN]    = fac * 2.0 * tdGlobCoord0 * tdGlobCoord0 / d; // dGlobCoord[0]
            values[hNN+1]  = fac * 2.0 * tdGlobCoord1 * tdGlobCoord0 / d; // dGlobCoord[1]
            values[hNE1]   = fac * 2.0 * tdGlobCoord0 * (d * tdGlobCoord10 - t * ddGlobCoord10) / (d*d); // dGlobCoord1[0]
            values[hNE1+1] = fac * 2.0 * (d * tdGlobCoord11 * tdGlobCoord0 - d * t - t * tdGlobCoord0 * ddGlobCoord11) / (d*d); // dGlobCoord1[1]
            values[hNE2]   = fac * 2.0 * tdGlobCoord0 * (d * tdGlobCoord20 - t * ddGlobCoord20) / (d*d); // dGlobCoord2[0]
            values[hNE2+1] = fac * 2.0 * (d * tdGlobCoord21 * tdGlobCoord0 + d * t - t * tdGlobCoord0 * ddGlobCoord21) / (d*d); // dGlobCoord2[1]
            // dGlobCoord[1] -> sig * ndGlobCoord1 / d = sig * 2.0 * t * tdGlobCoord1 / d
            values[hNN+2]  = fac * 2.0 * tdGlobCoord1 * tdGlobCoord1 / d; // dGlobCoord[1]
            values[hNE1+2] = fac * 2.0 * (d * tdGlobCoord10 * tdGlobCoord1 + d * t - t * tdGlobCoord1 * ddGlobCoord10) / (d*d); // dGlobCoord1[0]
            values[hNE1+3] = fac * 2.0 * tdGlobCoord1 * (d * tdGlobCoord11 - t * ddGlobCoord11) / (d*d); // dGlobCoord1[1]
            values[hNE2+2] = fac * 2.0 * (d * tdGlobCoord20 * tdGlobCoord1 - d * t - t * tdGlobCoord1 * ddGlobCoord20) / (d*d); // dGlobCoord2[0]
            values[hNE2+3] = fac * 2.0 * tdGlobCoord1 * (d * tdGlobCoord21 - t * ddGlobCoord21) / (d*d); // dGlobCoord2[1]
            // dGlobCoord1[0] -> sig * (d * ndGlobCoord10 - n * ddGlobCoord10 ) / (d*d) = sig * (2.0 * t * tdGlobCoord10 / d - n * ddGlobCoord10 / (d*d))
            values[hE1E1]   = fac * ( 2.0 * tdGlobCoord10 * (d * tdGlobCoord10 - t * ddGlobCoord10) / (d*d) - (d*d * (ndGlobCoord10 * ddGlobCoord10 + n * 2.0) - n * ddGlobCoord10 * 2.0 * d * ddGlobCoord10) / (d*d*d*d) ); // dGlobCoord1[0]
            values[hE1E1+1] = fac * ( 2.0 * tdGlobCoord10 * (d * tdGlobCoord11 - t * ddGlobCoord11) / (d*d) - ddGlobCoord10 * (d*d * ndGlobCoord11 - n * 2.0 * d * ddGlobCoord11) / (d*d*d*d) ); // dGlobCoord1[1]
            values[hE1E2]   = fac * ( 2.0 * tdGlobCoord10 * (d * tdGlobCoord20 - t * ddGlobCoord20) / (d*d) - (d*d * (ndGlobCoord20 * ddGlobCoord10 - n * 2.0) - n * ddGlobCoord10 * 2.0 * d * ddGlobCoord20) / (d*d*d*d) ); // dGlobCoord2[0]
            values[hE1E2+1] = fac * ( 2.0 * (d * tdGlobCoord21 * tdGlobCoord10 - d * t - t * tdGlobCoord10 * ddGlobCoord21) / (d*d) - ddGlobCoord10 * (d*d * ndGlobCoord21 - n * 2.0 * d * ddGlobCoord21) / (d*d*d*d) ); // dGlobCoord2[1]
            // dGlobCoord1[1] -> sig * (d * ndGlobCoord11 - n * ddGlobCoord11 ) / (d*d) = sig * (2.0 * t * tdGlobCoord11 / d - n * ddGlobCoord11 / (d*d))
            values[hE1E1+2] = fac * ( 2.0 * tdGlobCoord11 * (d * tdGlobCoord11 - t * ddGlobCoord11) / (d*d) - (d*d * (ndGlobCoord11 * ddGlobCoord11 + n * 2.0) - n * ddGlobCoord11 * 2.0 * d * ddGlobCoord11) / (d*d*d*d) ); // dGlobCoord1[1]
            values[hE1E2+2] = fac * ( 2.0 * (d * tdGlobCoord20 * tdGlobCoord11 + d * t - t * tdGlobCoord11 * ddGlobCoord20) / (d*d) - ddGlobCoord11 * (d*d * ndGlobCoord20 - n * 2.0 * d * ddGlobCoord20) / (d*d*d*d) ); // dGlobCoord2[0]
            values[hE1E2+3] = fac * ( 2.0 * tdGlobCoord11 * (d * tdGlobCoord21 - t * ddGlobCoord21) / (d*d) - (d*d * (ndGlobCoord21 * ddGlobCoord11 - n * 2.0) - n * ddGlobCoord11 * 2.0 * d * ddGlobCoord21) / (d*d*d*d) ); // dGlobCoord2[1]
            // dGlobCoord2[0] -> sig * (d * ndGlobCoord20 - n * ddGlobCoord20 ) / (d*d) = sig * (2.0 * t * tdGlobCoord20 / d - n * ddGlobCoord20 / (d*d))
            values[hE2E2]   = fac * ( 2.0 * tdGlobCoord20 * (d * tdGlobCoord20 - t * ddGlobCoord20) / (d*d) - (d*d * (ndGlobCoord20 * ddGlobCoord20 + n * 2.0) - n * ddGlobCoord20 * 2.0 * d * ddGlobCoord20) / (d*d*d*d) ); // dGlobCoord2[0]
            values[hE2E2+1] = fac * ( 2.0 * tdGlobCoord20 * (d * tdGlobCoord21 - t * ddGlobCoord21) / (d*d) - ddGlobCoord20 * (d*d * ndGlobCoord21 - n * 2.0 * d * ddGlobCoord21) / (d*d*d*d) ); // dGlobCoord2[1]
            // dGlobCoord2[1] -> sig * (d * ndGlobCoord21 - n * ddGlobCoord21 ) / (d*d) = sig * (2.0 * t * tdGlobCoord21 / d - n * ddGlobCoord21 / (d*d))
            values[hE2E2+2] = fac * ( 2.0 * tdGlobCoord21 * (d * tdGlobCoord21 - t * ddGlobCoord21) / (d*d) - (d*d * (ndGlobCoord21 * ddGlobCoord21 + n * 2.0) - n * ddGlobCoord21 * 2.0 * d * ddGlobCoord21) / (d*d*d*d) ); // ddGlobCoord2[1]
#ifndef NDEBUG
            // LOG Output
            LOG_DBG3(ci) << "Hessian values: index " << (hNN) << ", value: " << values[hNN];
            LOG_DBG3(ci) << "Hessian values: index " << (hNN+1) << ", value: " << values[hNN+1];
            LOG_DBG3(ci) << "Hessian values: index " << (hNE1) << ", value: " << values[hNE1];
            LOG_DBG3(ci) << "Hessian values: index " << (hNE1+1) << ", value: " << values[hNE1+1];
            LOG_DBG3(ci) << "Hessian values: index " << (hNE2) << ", value: " << values[hNE2];
            LOG_DBG3(ci) << "Hessian values: index " << (hNE2+1) << ", value: " << values[hNE2+1];
            LOG_DBG3(ci) << "Hessian values: index " << (hNN+2) << ", value: " << values[hNN+2];
            LOG_DBG3(ci) << "Hessian values: index " << (hNE1+2) << ", value: " << values[hNE1+2];
            LOG_DBG3(ci) << "Hessian values: index " << (hNE1+3) << ", value: " << values[hNE1+3];
            LOG_DBG3(ci) << "Hessian values: index " << (hNE2+2) << ", value: " << values[hNE2+2];
            LOG_DBG3(ci) << "Hessian values: index " << (hNE2+3) << ", value: " << values[hNE2+3];
            LOG_DBG3(ci) << "Hessian values: index " << (hE1E1) << ", value: " << values[hE1E1];
            LOG_DBG3(ci) << "Hessian values: index " << (hE1E1+1) << ", value: " << values[hE1E1+1];
            LOG_DBG3(ci) << "Hessian values: index " << (hE1E2) << ", value: " << values[hE1E2];
            LOG_DBG3(ci) << "Hessian values: index " << (hE1E2+1) << ", value: " << values[hE1E2+1];
            LOG_DBG3(ci) << "Hessian values: index " << (hE1E1+2) << ", value: " << values[hE1E1+2];
            LOG_DBG3(ci) << "Hessian values: index " << (hE1E2+2) << ", value: " << values[hE1E2+2];
            LOG_DBG3(ci) << "Hessian values: index " << (hE1E2+3) << ", value: " << values[hE1E2+3];
            LOG_DBG3(ci) << "Hessian values: index " << (hE2E2) << ", value: " << values[hE2E2];
            LOG_DBG3(ci) << "Hessian values: index " << (hE2E2+1) << ", value: " << values[hE2E2+1];
            LOG_DBG3(ci) << "Hessian values: index " << (hE2E2+2) << ", value: " << values[hE2E2+2];
#endif
          }else{ // we are nearest to a point, the formula is much simpler
            // |v-p|^2, with v = GlobCoord and p = GlobCoordDiff
            // (GlobCoord[0] - GlobCoordDiff[0])*(GlobCoord[0] - GlobCoordDiff[0]) + (GlobCoord[1] - GlobCoordDiff[1])*(GlobCoord[1] - GlobCoordDiff[1]);
            // some of the entries are zero and therefore omitted
            UInt hNE1 = hessianNodeEleIdx[i][jmin]; // this contains: node-node1: xx xy yx yy, node-node2: xx xy yx yy (8 hessian entries)
            UInt hE1E1 = hessianEleEleIdx[jmin]; // contains 10 hessian entries first node1-node1 xx xy yy
            values[hNN]     = fac * 2.0;  // dGlobCoord[0] dGlobCoord[0]
            values[hNE1]    = fac * -2.0; // dGlobCoord[0] dGlobCoordDiff[0]
            values[hNN+2]   = fac * 2.0;  // dGlobCoord[1] dGlobCoord[1]
            values[hNE1+1]  = fac * -2.0; // dGlobCoord[1] dGlobCoordDiff[1]
            values[hE1E1]   = fac * 2.0;  // dGlobCoordDiff[0] dGlobCoordDiff[0]
            values[hE1E1+2] = fac * 2.0;  // dGlobCoordDiff[1] dGlobCoordDiff[1]
#ifndef NDEBUG
            // LOG Output
            LOG_DBG3(ci) << "Hessian values: index " << (hNN) << ", value: " << values[hNN];
            LOG_DBG3(ci) << "Hessian values: index " << (hNE1) << ", value: " << values[hNE1];
            LOG_DBG3(ci) << "Hessian values: index " << (hNN+2) << ", value: " << values[hNN+2];
            LOG_DBG3(ci) << "Hessian values: index " << (hNE1+1) << ", value: " << values[hNE1+1];
            LOG_DBG3(ci) << "Hessian values: index " << (hE1E1) << ", value: " << values[hE1E1];
            LOG_DBG3(ci) << "Hessian values: index " << (hE1E1+2) << ", value: " << values[hE1E1+2];
#endif
          } // if(segment)
        } // for node in nodeList
      } // for dir
    } // for bilateralcontactconstraint
    assert(cidx == (UInt)m);
  } // if structure / values
  return true;
}

void ContactIpopt::TNLP::finalize_solution(SolverReturn status,
                              Index n, const Number* x, const Number* z_L, const Number* z_U,
                              Index m, const Number* g, const Number* lambda,
                              Number obj_value, const IpoptData* ip_data, IpoptCalculatedQuantities* ip_cq)
{
  LOG_TRACE(ci) << "finalize_solution: status = " << status << "; n = " << n << "; m = " << m 
                    << " obj_value = " << obj_value << " x_avg = " << Average(x, n) << " x_std_dev = " 
                    << StandardDeviation(x, n) << " z_l_avg = " << Average(z_L, n) << " z_l_std_dev = "
                    << StandardDeviation(z_L, n) << " z_u_avg = " << Average(z_L, n) << " z_u_std_dev = "
                    << StandardDeviation(z_U, n);
  // save back the solution
  for(Index i = 0; i < n; i++){
    contactIpopt_->sol_->SetEntry(i, x[i]);
    contactIpopt_->x[i] = x[i];
  }
  if(!contactIpopt_->nextisadjoint){ // we have solved forward problem
    for(Index i = 0; i < m; i++){
      contactIpopt_->lambda[i] = lambda[i];
    }
    contactIpopt_->initlambda = true;
  }
}

bool ContactIpopt::TNLP::intermediate_callback(AlgorithmMode mode,Index iter, Number obj_value,
     Number inf_pr, Number inf_du, Number mu, Number d_norm,Number regularization_size,
     Number alpha_du, Number alpha_pr, Index ls_trials, const IpoptData* ip_data,
     IpoptCalculatedQuantities* ip_cq)
{
  LOG_TRACE2(ci) << "intermediate_callback: mode = " << mode << "; iter = " << iter 
                    << " obj_value = " << obj_value << "; ls_trials = " << ls_trials;

  return true;
}
