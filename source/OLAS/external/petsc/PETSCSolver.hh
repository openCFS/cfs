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
#include <def_expl_templ_inst.hh>

#include "General/Environment.hh"
#include "OLAS/solver/BaseSolver.hh"
#include "PDE/BasePDE.hh"
#include "Domain/Mesh/Grid.hh"
#include "Domain/Domain.hh"
#include "General/Enum.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "DataInOut/Logging/log.hpp"
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
#define GET_SOL 7
#define SOLVER_STRING 8
#define SETUP_MG 9
#define ISSYMMETRIC 10
#define HOM_DIR_ELIMINATION 11


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

/** To satisfy the compiler
* @param sysmat shall be the one Setup() is called with */
void Solve( const BaseMatrix &sysmat, const BaseVector &rhs, BaseVector &sol);

void SendWorkerCommand(int TAG);

//    void DMDAGetElements_3D(DM dm,PetscInt *nel,PetscInt *nen,const PetscInt *e[]);
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

Vec N_;

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

//To find the rank of the processor its currently in
int rank_=0;
int size_=0;
//    //Coarsegrid solver parameters
//    PetscScalar coarse_rtol = 1.0e-8;
//    PetscScalar coarse_atol = 1.0e-50;
//    PetscScalar coarse_dtol = 1e3;
//    PetscInt coarse_maxits = 30;
//
//    //Number of levels
//    PetscInt nlvls=2;
//
//    // Number of smoothening iterations per up/down smooth_sweeps
//    PetscInt smooth_sweeps = 4;

std::string solverstring_;
std::string precondstring_;
bool MG_FLAG =false;
bool symmetric=false;

};


class PETSCWorker : public PETScCommon
{
public:

void run();
PETSCWorker(int argc,const char **argv);
~PETSCWorker();

private:

void InitPetscWorker();

void HomDirElimination();
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

Vec N_;

/** with throw exception when exceeded */
PetscInt maxIter_ = 10000;

/** note that PETSC has three ways to calculate the residuum */
PetscScalar tolerance_ = 1e-12;

/** throws only an exception on maxIter exceeded if also minTol_ is not met */
PetscScalar minTol_ =  1e-11;


////Coarsegrid solver parameters
//PetscScalar coarse_rtol = 1.0e-8;
//PetscScalar coarse_atol = 1.0e-50;
//PetscScalar coarse_dtol = 1e10;
//PetscInt coarse_maxits = 30;
//
////Number of levels
//PetscInt nlvls=4;
//
//// Number of smoothening iterations per up/down smooth_sweeps
//PetscInt smooth_sweeps = 4;
//

//pointer to xml node
PtrParamNode xml_;

bool MG_FLAG=false;
bool symmetric=false;
//Strings for setting Solver and Preconditioner
std::string solverstring_;
std::string precondstring_;
ParamNodeList regionList_;

};
  
}
#endif
