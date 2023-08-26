/*
 * PhistLinearSolver.hh
 *
 *  Created on: Mar 7, 2018
 *      Author: sri
 */

#ifndef OLAS_EXTERNAL_PHIST_PHISTLINEARSOLVER_HH_
#define OLAS_EXTERNAL_PHIST_PHISTLINEARSOLVER_HH_

#include "PhistCore.hh"



namespace CoupledField {
class BaseMatrix;
class BaseVector;
class Flags;
class sparseMat_t; // phist matrix type
class StdMatrix;


class PhistLinearSolver  : public BaseIterativeSolver , PhistCore
{

public:

PhistLinearSolver(PtrParamNode param, PtrParamNode olasInfo, BaseMatrix::EntryType type);

void Setup(BaseMatrix &sysmat);


//! Dummy method: Notify the solver that a new matrix pattern has been set
void SetNewMatrixPattern() {EXCEPTION("SetNewMatrixPattern not implemented!");};


void Solve( const BaseMatrix &sysmat, const BaseVector &rhs, BaseVector &sol);


BaseSolver::SolverType GetSolverType() { return BaseSolver::PHIST; }


~PhistLinearSolver();

private:

int iflag_=0;
PtrParamNode xml_;
int maxIter_ = 1000;
double tolerance_ = 1.0;
int blockSize_=4;
int nrhs_=1;

bool firstSetup_=true;
StdVector<double> rhsVec_;
StdVector<double>  solVec_;

sparseMat_t* A_ = nullptr;
phist::types<double>::sparseMat_ptr phistMat_=nullptr;
phist::types<double>::mvec_ptr sol_ = nullptr;
phist::types<double>::mvec_ptr rhs_ = nullptr;
phist::types<double>::linearOp_ptr  A_op_=nullptr;
phist::types<double>::linearOp_ptr  P_op_=nullptr;



};

} /* namespace CoupledField */

#endif /* OLAS_EXTERNAL_PHIST_PHISTLINEARSOLVER_HH_ */
