// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     PETSCSolver.hh
 *       \brief    This implements the interface to the PETSC solver package
 *        
 *        The interface consist of  a solver class (ie master)  and a worker class. The sysmat and rhs
 *        obtained is set into PetscMat only using the master class ,all other commands must be called
 *        by all workers
 *              
 *        
 *       \date     oct 30, 2017
 *       \author   sri
 */
//================================================================================================



#ifndef PETSCSOLVER_HH_
#define PETSCSOLVER_HH_




#include "General/Environment.hh"
#include "OLAS/solver/BaseSolver.hh"
#include "PDE/BasePDE.hh"
#include "Domain/Mesh/Grid.hh"
#include "Domain/Domain.hh"
#include "General/Enum.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "PETScCommon.hh"

// include the original PETSC header
//extern "C"

#include "petsc.h"
//#include "petsc/private/dmdaimpl.h"

#define DIETAG 0
#define SETUP_MATRIX 1
#define ASSEMBLE_MAT 2
#define ASSEMBLE_VEC_RHS  3
#define SETUP_SOLVER_CONTEXT 4
#define SOLVE 5
#define DATA 6
#define GET_GLOBAL_VEC 7
#define SOLVER_STRING 8
#define SETUP_MG 9
#define ISSYMMETRIC 10
#define SET_HOM_DIR_VEC 11
#define HOM_DIR_PENALTY 12
#define ASSEMBLE_VEC_DIR 13
#define GET_SOL 14
#define MAT_ZERO_ENTRIES 15



namespace CoupledField
{
class BaseMatrix;
class BaseVector;
class Flags;


class PETSCSolver : public BaseIterativeSolver ,  PETScCommon
{
public:
PETSCSolver(PtrParamNode param, PtrParamNode olasInfo, BaseMatrix::EntryType type);

virtual ~PETSCSolver();

/** Every call sets up a new preconditionier. */
void Setup(BaseMatrix &sysmat);

//! Dummy method: Notify the solver that a new matrix pattern has been set
void SetNewMatrixPattern() {EXCEPTION("SetNewMatrixPattern not implemented for PETSC. GetRidOfZeros for NCIs will not work.");};

void Solve( const BaseMatrix &sysmat, const BaseVector &rhs, BaseVector &sol);

void SendWorkerCommand(int TAG);


BaseSolver::SolverType GetSolverType() { return BaseSolver::PETSC; }

private:
//PETSC Error Code
PetscErrorCode ierr=0;

PC precond_;//ksp linear precond context

KSP solver_; //ksp linear solver context

//stores the system Matrix
Mat A_;
//stores the current solution
Vec x_;
//stores the current RHS
Vec b_;

Vec dirNodeVec_;

// Create DM which are grid management
DM da_nodes;
//pointer to xml node
PtrParamNode xml_;

//internal status flag for updated matrix computations
bool firstSetup_=true;
/** with throw exception when exceeded */
PetscInt maxIter_ = -1;
/** not that PETSC has three ways to calculate the residuum */
PetscScalar tolerance_ = -1.0;
/** throws only an exception on maxIter exceeded if also minTol_ is not met */
PetscScalar minTol_ = -1.0;

/** Logs extra info when enabled */
bool logging_ = false;


std::string solverstring_;
std::string precondstring_;
bool MG_FLAG =false;
bool symmetric=false;
Vec N_;//dirVector which consist of 0 when a eqnNr corresponds to HomDirBC ,all other place the value is 1
StdVector<unsigned int> cfsEqnMap_;
Vec dirNodeVecGlobal_=nullptr;

//Info dumped to .info.xml
PetscInt totalSolverIter=0; //This sums iterations from all the optimization iteration. Useful for understanding mg.
PetscInt niter=0;



};


class PETSCWorker : public PETScCommon
{
public:

void run();
PETSCWorker(int argc,const char **argv);
~PETSCWorker();

private:

void InitPetscWorker();
void AssembleMatrix();
//PETSC Error Code
PetscErrorCode ierr=0;

PC precond_;//ksp linear precond context

KSP solver_; //ksp linear solver context

//stores the system Matrix
Mat A_;

//stores the current solution
Vec x_;

///stores the current RHS
Vec b_;
// Create DM which are grid management
DM da_nodes;

Vec dirNodeVec_;

/** with throw exception when exceeded */
PetscInt maxIter_ = 10000;

/** note that PETSC has three ways to calculate the residuum */
PetscScalar tolerance_ = 1e-12;

/** throws only an exception on maxIter exceeded if also minTol_ is not met */
PetscScalar minTol_ =  1e-11;

PtrParamNode xml_;

bool MG_FLAG=false;
bool symmetric=true;
//Strings for setting Solver and Preconditioner
std::string solverstring_;
std::string precondstring_;
Vec N_;
Vec dirNodeVecGlobal_=nullptr;
//Number of nodes in x y and z directions respectively in case of structured with hex elem it is always elemNum in one direction +1
int nx=0,ny=0,nz=0,dimension=0;

};
  
}
#endif
