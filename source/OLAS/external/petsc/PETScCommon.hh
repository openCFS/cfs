/*
 * PETScCommon.hh
 *
 *  Created on: Feb 23, 2018
 *      Author: sri
 */

#ifndef OLAS_EXTERNAL_PETSC_PETSCCOMMON_HH_
#define OLAS_EXTERNAL_PETSC_PETSCCOMMON_HH_


#include <def_expl_templ_inst.hh>

#include "General/Environment.hh"
#include "OLAS/solver/BaseSolver.hh"
#include "PDE/BasePDE.hh"
#include "Domain/Mesh/Grid.hh"
#include "Domain/Domain.hh"
#include "General/Enum.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "DataInOut/Logging/log.hpp"
// include the original PETSC header
#include "petsc.h"

namespace CoupledField {
class BaseMatrix;
class BaseVector;
class Flags;


class PETScCommon
{
public:
PETScCommon();
void SetupMatrix();
void AssembleMatrixMG(Mat &sysMat,DM &daNodes,Vec &solVec,Vec &rhsVec,Vec &dirNodeVec,Vec &dirVec,int nx,int ny,int nz);
void GetHomDirNodes(Vec &dirVec);
void SetLinRhs(Vec &rhsVec);
void GetCFSEqnMapMG(StdVector<unsigned int> &cfsEqnMap);
void GetSolution(bool multigrid);
void CreateDMDA(DM & daNodes,Mat &sysMat,Vec &solVec,Vec &rhsVec,Vec &dirVec,int nx,int ny,int nz,int dim);
void SetupMGSolver(DM &da_nodes,PC &precond_);
void SetupMatrix(UInt dim,UInt* nnzr,bool symmetric, Mat &sysMat,Vec &solVec, Vec &rhsVec);
void SetupSolverContext(Mat &sysMat,KSP &solContext,PC &preContext,string solverstring,string precondstring,PetscScalar tolerance ,PetscScalar minTol,PetscInt maxIter,bool MG_FLAG);
void GetGlobalVec(Vec &x,Vec &xGlobal);
void GetGridInfoMG(int &nx,int &ny,int &nz,int &dimension);
void PenaltyHomDir(Mat &sysMatVec,Vec &rhsVec,Vec &f);
~PETScCommon();
std::string CreateSolverString(PtrParamNode);
std::string CreatePrecondString(PtrParamNode);

private:
PetscErrorCode ierr=0;
//Coarsegrid solver parameters
PetscScalar coarse_rtol = 1.0e-8;
PetscScalar coarse_atol = 1.0e-50;
PetscScalar coarse_dtol = 1e3;
PetscInt coarse_maxits = 30;

//Number of levels
PetscInt nlvls=2;

// Number of smoothening iterations per up/down smooth_sweeps
PetscInt smooth_sweeps = 4;
};

} /* namespace CoupledField */

#endif /* OLAS_EXTERNAL_PETSC_PETSCCOMMON_HH_ */
