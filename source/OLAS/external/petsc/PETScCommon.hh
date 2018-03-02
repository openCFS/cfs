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
void AssembleMatrixMG(Mat &sysMat,DM &daNodes,Vec &solVec,Vec &rhsVec,Vec &dirVec,StdVector<int> &homDirEquNr,StdVector<int> &dirchletValue);
void GetHomDirMG(Vec &dirVec);
void SetLinRhs(Vec &rhsVec);
void GetCFSEqnMapMG(StdVector<unsigned int> &cfsEqnMap);
void GetSolution(bool multigrid);
void CreateDMDA(DM & daNodes,Mat &sysMat,Vec &solVec,Vec &rhsVec,Vec &dirVec,int nx,int ny,int nz,int dim);
void SetupMGSolver();
void SetupMatrix(UInt dim,UInt* nnzr,bool symmetric, Mat &sysMat,Vec &solVec, Vec &rhsVec);
void SetupSolverContext(Mat &sysMat,KSP &solContext,PC &preContext,string solverstring,string precondstring,PetscScalar tolerance ,PetscScalar minTol,PetscInt maxIter,bool MG_FLAG);
void GetSol(Vec &solVec,Vec &solGlobal);
void GetGridInfoMG(int &nx,int &ny,int &nz,int &dimension);
~PETScCommon();
std::string CreateSolverString(PtrParamNode);
std::string CreatePrecondString(PtrParamNode);

private:
PetscErrorCode ierr=0;

};

} /* namespace CoupledField */

#endif /* OLAS_EXTERNAL_PETSC_PETSCCOMMON_HH_ */
