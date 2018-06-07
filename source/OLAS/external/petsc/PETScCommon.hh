/*
 * PETScCommon.hh
 *
 *  Created on: Feb 23, 2018
 *      Author: sri
 *      Description: This class contains methods that are used by both petsc worker and the main solver classes. Since
 *      the master and worker needs to call the same functions but cannot be derived from the same base class as the worker
 *      class is required in the main function and PETSCSolver class is not called until solve step.
 *
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
#include "Forms/LinForms/LinearForm.hh"


namespace CoupledField {
class BaseMatrix;
class BaseVector;
class Flags;


class PETScCommon
{

protected:

PETScCommon();
~PETScCommon();
void SetupMatrix();
void CheckLevels(int nx,int ny,int nz);
void AssembleMatrixMG(Mat &sysMat,DM &daNodes,Vec &solVec,Vec &rhsVec,Vec &dirNodeVec,Vec &dirVec,int nx,int ny,int nz);
void GetHomDirNodes(Vec &dirVec);
void SetLinRhs(Vec &rhsVec);
void GetCFSEqnMapMG(StdVector<unsigned int> &cfsEqnMap);
void GetSolution(bool multigrid);
void CreateDMDA(DM & daNodes,Mat &sysMat,Vec &solVec,Vec &rhsVec,Vec &dirVec,int nx,int ny,int nz,int dim);
void SetupMGSolver(DM &da_nodes,PC &precond_);
void SetupMatrix(UInt dim,UInt* nnzr,bool symmetric, Mat &sysMat,Vec &solVec, Vec &rhsVec);
void SetupSolverContext(Mat &sysMat,KSP &solContext,PC &preContext,string solverstring,string precondstring,PetscScalar tolerance ,PetscScalar minTol,PetscInt maxIter,bool MG_FLAG);
void GetGlobalVec(Vec &x,Vec &xGlobal,bool master);
void GetGridInfoMG(int &nx,int &ny,int &nz,int &dimension);
void PenaltyHomDir(Mat &sysMatVec,Vec &rhsVec,Vec &f);
void CacluateElementStiffnessMatrix(EntityIterator &firstEntIt, EntityIterator &secondEntIt,StdVector<int> &globalNodeNum,
    int &nen , BiLinearForm * form ,  int &elem,PetscScalar KE[]);
void CalculateDirNodesAndEdof(Vec &dirVec,PetscInt edof[],const int &nen,const PetscScalar * dirArray,
    const StdVector<int> &globalNodeNum,const int &elem);


std::string CreateSolverString(PtrParamNode);
std::string CreatePrecondString(PtrParamNode);

int rank_=0;
int size_=0;
PetscInt nlvls=2;
//Assemble context
Assemble * assemble_=nullptr;
bool disableParalleAssemby=true;



protected:
PetscErrorCode ierr=0;
//Coarsegrid solver parameters
PetscScalar coarse_rtol = 1.0e-8;
PetscScalar coarse_atol = 1.0e-50;
PetscScalar coarse_dtol = 1e3;
PetscInt coarse_maxits = 30;

//Number of levels

//just for warning in recursion loop
int MGLevels=nlvls;



// Number of smoothening iterations per up/down smooth_sweeps
PetscInt smooth_sweeps = 4;

std::string innerSovler;



};

} /* namespace CoupledField */

#endif /* OLAS_EXTERNAL_PETSC_PETSCCOMMON_HH_ */
