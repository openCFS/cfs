/*
 * PhistLinearSolver.cc
 *
 *  Created on: Mar 7, 2018
 *      Author: sri
 */

#include "PhistLinearSolver.hh"


namespace CoupledField {

PhistLinearSolver::PhistLinearSolver(PtrParamNode param, PtrParamNode olasInfo, BaseMatrix::EntryType type) {


  A_op_=NULL;
  xml_=param;

  infoNode_=olasInfo->Get("phist");
  tolerance_=xml_->Has("tolerance") ? xml_->Get("tolerance")->As<double>() : 1e-6;
  blockSize_ =xml_->Has("blockSize") ?  xml_->Get("blockSize")->As<int>():4;
  maxIter_= xml_->Has("maxIter") ?xml_->Get("maxIter")->As<int>():10000;

  PtrParamNode hdr = infoNode_->Get(ParamNode::HEADER);
  hdr->Get("maxIter")->SetValue(maxIter_);
  hdr->Get("tolerance")->SetValue(double(tolerance_));
}


void PhistLinearSolver::Setup(BaseMatrix &sysmat){

  shared_ptr<Timer> setup_phistMatrix = infoNode_->Get(ParamNode::SUMMARY)->Get("setup_phistMatrix/timer")->AsTimer(GetSetupTimer());
  setup_phistMatrix->Start();

  phistMat_ = (phist::types<double>::sparseMat_ptr) InitMatrix<double>(sysmat, &A_, 1.0);

  setup_phistMatrix->Stop();
  if (firstSetup_){
    A_op_ = new phist::types<double>::linearOp();
    P_op_ = new phist::types<double>::linearOp();
    P_op_=NULL;
    // we need the domain map of the matrix
    phist_const_map_ptr map = NULL;

    phist::kernels<double>::sparseMat_get_domain_map(phistMat_,&map,&iflag_);
    assert(iflag_ == 0);

    phist::core<double>::linearOp_wrap_sparseMat(A_op_,phistMat_,&iflag_);

    phist::kernels<double>::mvec_create(&sol_,map,1,&iflag_);
    assert(iflag_ == 0);
    phist::kernels<double>::mvec_create(&rhs_,map,1,&iflag_);
    assert(iflag_ == 0);
  }

 }

void PhistLinearSolver::Solve(const BaseMatrix &sysmat, const BaseVector &rhs, BaseVector &sol){

  if (firstSetup_){
    rhsVec_.Resize(rhs.GetSize());
    solVec_.Resize(rhs.GetSize());
  }


  for(unsigned int i=0; i<rhs.GetSize(); i++){
    Double myEnt =0;
    rhs.GetEntry(i,myEnt);
    rhsVec_[i]=myEnt;
  }
  int rowMajor=0;
  int rowSize=1;



  phist::kernels<double>::mvec_set_data(rhs_,rhsVec_.GetPointer(),rowSize,rowMajor,&iflag_);
  assert(iflag_==0);


  // TODO Add all other solvers(uncomment or add as struct) when phist is updated further
  // phist::krylov<double>::restartedBlockedGMRES(A_op_,P_op_,rhs_,sol_,nrhs_,&maxIter_,&tolerance_,blockSize_,maxBlockSize,&iflag_);
  // phist::krylov<double>::restartedGMRES(A_op_,P_op_,rhs_,sol_,&maxIter_,tolerance_,0,&iflag_);
  // phist::krylov<double>::PCG(A_op_,NULL,rhs_,sol_,&maxIter_,tolerance_,&iflag_);

  phist::krylov<double>::BiCGStab(A_op_,P_op_,rhs_,sol_,&maxIter_,tolerance_,&iflag_);
  assert(iflag_==0);

  phist::kernels<double>::mvec_get_data(sol_,solVec_.GetPointer(),rowSize,rowMajor,&iflag_);
  assert(iflag_==0);

  for (unsigned int j=0;j<solVec_.GetSize();j++){
    sol.SetEntry(j,solVec_[j]);
  }


}

PhistLinearSolver::~PhistLinearSolver() {
  delete []  A_op_;
  delete [] P_op_;
  phist::kernels<double>::sparseMat_delete(A_, &iflag_);
  phist::kernels<double>::sparseMat_delete(phistMat_, &iflag_);
  phist::kernels<double>::mvec_delete(sol_, &iflag_);
  phist::kernels<double>::mvec_delete(rhs_, &iflag_);
}

} /* namespace CoupledField */
