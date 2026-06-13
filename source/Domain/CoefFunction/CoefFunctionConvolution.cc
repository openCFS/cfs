/*
 * CoefFunctionConvolution.cc
 *
 *  Created on: 29 Jan 2016
 *      Author: INTERN\ftoth
 */
#include "CoefFunctionConvolution.hh"
#include <cmath>
#include "DataInOut/Logging/LogConfigurator.hh"

namespace CoupledField {
  DEFINE_LOG(visco, "visco")

  void CoefFunctionConvolution::Init() {

    StdVector<LocPoint> intPoints; // only used for call of GetIntPoints()
    StdVector<Double> weights; // only used for call of GetIntPoints()

    // initialize map for data values
    EntityIterator it = entityList_->GetIterator();

    // this could be simplified
    UInt n = 0;
    StdVector<UInt> nids;
    std::map<UInt, UInt>::iterator it2;
    for ( it.Begin(); !it.IsEnd(); it++)
      {
        const Elem* ptEl = it.GetElem();
        nids = ptEl->connect;
        for (UInt i=0; i<nids.GetSize(); i++)
          {
            it2 = indexMap_.find(nids[i]);
            if (it2 == indexMap_.end())
              {
                indexMap_[nids[i]] = n;
                n = n+vecSize_; // increase by number of components, which are ordered conscutively
              }
          }
      }

    // allocate value vector
    histvals_.Resize(n,0.0);
    u_.Resize(n,0.0);

    // callback if time changes
    // obtain handle from internal variable coefficient function
    mp_ = domain->GetMathParser();
    mHandleStep_ = mp_->GetNewHandle(true);
    mp_->SetExpr(mHandleStep_,"step"); // current step
    mHandleTime_ = mp_->GetNewHandle(true);
    mp_->SetExpr(mHandleTime_,"t"); // current step
    mHandleDt_ = mp_->GetNewHandle(true);
    mp_->SetExpr(mHandleDt_,"dt"); // current step
    //t_ = this->mp_->Eval(mHandle_); // eval time
    // register callback mechanism if the step changes
    // this is called automatically at the beginning of each step
    mp_->AddExpChangeCallBack( boost::bind(&CoefFunctionConvolution::ComputeRecursion, this ), mHandleStep_ );
    //ComputeRecursion(); // call to set f1_ and f2_
  }

  CoefFunctionConvolution::CoefFunctionConvolution(Double tau, shared_ptr<ElemList> actList, shared_ptr<FeSpace> ptFeSpace, BaseBOperator* bOp ) : CoefFunction () {
    dimType_ = VECTOR; // should be determined from input
    dependType_ = SOLUTION;
    isAnalytic_ = false;
    isComplex_  = false;
    vecSize_= ptFeSpace->GetNumDofs();

    tau_ = tau;
    // allocate history variables
    entityList_ = actList;
    // set fespace and inegration scheme
    ptFeSpace_ = ptFeSpace;
    intScheme_ = ptFeSpace_->GetIntScheme();
    ptFeFct_ = shared_ptr<BaseFeFunction>(ptFeSpace_->GetFeFunction());
    // sould be used ptFeSpace->GetSpaceType()
    // dimension should be decided
    bOp_ = bOp;
    Init();

  }

  CoefFunctionConvolution::~CoefFunctionConvolution(){
    ;
  }

  void CoefFunctionConvolution::ComputeRecursion() {

    // build up vector of solution variables
    Vector<Double> u_2;
    u_2.Resize(u_.GetSize());
    // make a node list for all nodes in region
    NodeList nodeList(entityList_->GetGrid());
    nodeList.SetNodesOfRegion(entityList_->GetRegion());
    EntityIterator it = nodeList.GetIterator();
    Vector<Double>  sol; // solution vector
    UInt pos; // position index
    for ( it.Begin(); !it.IsEnd(); it++)
      {
        ptFeFct_->GetEntitySolution(sol, it); // get solution at node
        pos = GetId(it.GetNode()); // find correct position
        for (UInt i=0; i<sol.GetSize();i++) // set all vector components
          {
            u_2[pos+i] = sol[i];
          }
      }

    // compute current factors for update
    Double dt = this->mp_->Eval(mHandleDt_);
    lam_ = dt/tau_;
    f1_ = exp(-lam_);
    if (lam_<1.0e-6) { // lam->0:  (1-exp(-lam))/lam = 0/0 = exp(-lam)
        f2_ = f1_;
    }
    else {
        f2_ = tau_*(1-f1_)/dt;
    }

    Vector<Double> cvals; // current (updated) values
    ComputeCurrentVals(cvals,&histvals_,&u_,&u_2);

    // debug output
    Double t = this->mp_->Eval(mHandleTime_);
    UInt step = this->mp_->Eval(mHandleStep_);
    LOG_DBG(visco) << "step-"<<step<<" "<<t<<" u "<<u_.ToString();
    Vector<Double> du; du.Add(1.0,u_2,-1.0,u_);
    LOG_DBG(visco) << "step-"<<step<<" "<<t<<" du "<<du.ToString();
    LOG_DBG(visco) << "step-"<<step<<" "<<t<<" h "<<cvals.ToString();

    // update history variables
    u_ = u_2;
    histvals_ = cvals;
  }

  void CoefFunctionConvolution::GetVector(Vector<Double>& retvec, const LocPointMapped& lpm ) {

    // get variables for all nodes of the element
    const StdVector<UInt>& nids = lpm.ptEl->connect;
    // history variables at nodes
    Vector<Double> h_n = GetSavedVals(nids,&histvals_);
    // node variables of last (converged) step
    Vector<Double> u_n = GetSavedVals(nids,&u_);
    // current solution
    Vector<Double> u;
    ptFeFct_->GetEntitySolution(u,lpm.ptEl);

    // compute update
    Vector<Double> cvals;
    ComputeCurrentVals(cvals,&h_n,&u_n,&u);

    //get operator matrix
    Matrix<Double> bMat;
    BaseFE* ptFe = ptFeSpace_->GetFe( lpm.ptEl->elemNum ); // get the FE from the feSpace
    bOp_->CalcOpMat(bMat,lpm,ptFe);
    // compute result
    retvec = bMat*(cvals);
  }

  UInt CoefFunctionConvolution::GetVecSize() const {
    return bOp_->GetDimDMat();
  }

  UInt CoefFunctionConvolution::GetId(UInt nodeId) {
    UInt ret;
    std::map<UInt, UInt>::iterator it = indexMap_.find(nodeId);
    if (it != indexMap_.end())
      {
        ret = it->second;
      }
    else
      {
        EXCEPTION("nid not found - this should not happen!");
      }
    return ret;
  }

  Vector<Double> CoefFunctionConvolution::GetSavedVals(const StdVector<UInt> nids, const Vector<Double> * vals) {
    UInt N = nids.GetSize()*vecSize_;
    Vector<Double> ret(N);
    Double val;
    for (UInt i=0; i<nids.GetSize(); i++)
      {
        for (UInt j=0; j<vecSize_;j++)
          {
            vals->GetEntry(GetId(nids[i])+j,val);
            ret[i*vecSize_+j] = val;
          }
      }
    return ret;
  }

  void CoefFunctionConvolution::ComputeCurrentVals(Vector<Double>& res, const Vector<Double> * h0, const Vector<Double> * u0, const Vector<Double> * u) {
    // compute recursion u = f1*h0 + f2*(u-u0)
    res = *h0;
    res.ScalarMult(f1_);
    Vector<Double> du; du.Add(1.0,*u0,-1.0,*u);
    du.ScalarMult(f2_);
    res = res + du;
  }
}
