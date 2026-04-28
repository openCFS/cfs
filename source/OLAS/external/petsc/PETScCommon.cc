/*
 * PETScCommon.cc
 *
 *  Created on: Feb 23, 2018
 *      Author: sri
 */

#include "PETScCommon.hh"

#include <string>
#include <sstream>
#include <stdio.h>
#include <mpi.h>
#include "MatVec/SparseOLASMatrix.hh"
#include "MatVec/SCRS_Matrix.hh"
#include "MatVec/CRS_Matrix.hh"
#include "DataInOut/ProgramOptions.hh"
#include "DataInOut/ParamHandling/XmlReader.hh"
#include "Driver/Assemble.hh"
#include "Forms/LinForms/LinearForm.hh"
#include "PDE/SinglePDE.hh"
#include "Driver/TimeSchemes/BaseTimeScheme.hh"


namespace CoupledField {

#define CHKERRXX(ierr)  do {if (PetscUnlikely(ierr)) {PetscError(PETSC_COMM_SELF,__LINE__,PETSC_FUNCTION_NAME,__FILE__,ierr,PETSC_ERROR_IN_CXX,0);}} while(0)

PETScCommon::PETScCommon() {

  MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
  MPI_Comm_size(MPI_COMM_WORLD,&size_);

}

PETScCommon::~PETScCommon() {


}


void PETScCommon::PenaltyHomDir(Mat &sysMat,Vec &rhsVec,Vec &f){
  MatDiagonalScale(sysMat,f,f);
  Vec NI;
  VecDuplicate(f,&NI);
  VecSet(NI,1.0);
  VecAXPY(NI,-1.0,f);
  MatDiagonalSet(sysMat,NI,ADD_VALUES);
  // Zero out possible loads in the RHS that coincide
  // with Dirichlet conditions
  VecPointwiseMult(rhsVec,rhsVec,f);
  VecDestroy(&NI);
}


void PETScCommon::SetupMatrix(UInt dim,UInt* nnzr,bool symmetric, Mat &sysMat,Vec &solVec,Vec &rhsVec){

  if (symmetric){
    ierr=MatCreateSBAIJ(PETSC_COMM_WORLD,PetscInt(1),PETSC_DECIDE,PETSC_DECIDE,
        PetscInt(dim),PetscInt(dim),0,(PetscInt *)nnzr,0,(PetscInt *)nnzr,&sysMat);
    ierr=MatSetUp(sysMat);
    ierr=MatSetOption (sysMat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);CHKERRXX(ierr);
    ierr=MatCreateVecs(sysMat,&rhsVec,&solVec);CHKERRXX(ierr);
  }
  else {
  //USE THESE WHEN MATRIX IS OF TYPE SPARSE NON SYM
    ierr=MatCreate(PETSC_COMM_WORLD,&sysMat); CHKERRXX(ierr);
    ierr=MatSetSizes(sysMat,PETSC_DECIDE,PETSC_DECIDE,PetscInt(dim),PetscInt(dim));CHKERRXX(ierr);
    ierr=MatSetFromOptions(sysMat);CHKERRXX(ierr);
    ierr=MatSetUp(sysMat);CHKERRXX(ierr);
    ierr=MatSeqAIJSetPreallocation(sysMat,0,(PetscInt *)nnzr);CHKERRXX(ierr);
    ierr=MatMPIAIJSetPreallocation(sysMat,0,(PetscInt *)nnzr,0,(PetscInt *)nnzr);CHKERRXX(ierr);
    ierr=MatSetOption (sysMat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);CHKERRXX(ierr);
    ierr=MatCreateVecs(sysMat,&rhsVec,&solVec);CHKERRXX(ierr);
  }
}

void PETScCommon::SetupSolverContext(Mat &sysMat,KSP &solver,PC &precond,string solverstring,string precondstring,PetscScalar tol ,PetscScalar minTol,PetscInt maxIter,bool MG_FLAG  ){
  //Create Vector and solver objects and set parameters
  //create rhs and sol vector and allocate size in setup step
  //Allocate size for rhs vec assuming the size is equal to dim of sysmat
  //creates rhs and lhs vector from Matrix with appropriate size
  //setup linear solver context with preconditioner for petsc ksp methods //now hard coded
  ierr=KSPCreate(PETSC_COMM_WORLD,&(solver));CHKERRXX(ierr);
  ierr = KSPSetTolerances(solver,minTol,tol,PETSC_DEFAULT,maxIter);CHKERRXX(ierr);
  ierr = KSPSetInitialGuessNonzero(solver,PETSC_TRUE);CHKERRXX(ierr);
  ierr = KSPSetType(solver,solverstring.c_str());CHKERRXX(ierr);
  // The preconditinoer
  ierr=KSPGetPC(solver,&precond);CHKERRXX(ierr);
  // Set the preconditioner
  PCSetType(precond,precondstring.c_str());CHKERRXX(ierr);
  PetscOptionsInsert(NULL,NULL,NULL,NULL);CHKERRXX(ierr);
  // Set solver from options
  KSPSetFromOptions(solver);CHKERRXX(ierr);
  // Get the prec again - check if it has changed
  KSPGetPC(solver,&precond);CHKERRXX(ierr);
  ierr = KSPSetOperators(solver,sysMat,sysMat);CHKERRXX(ierr);

}


void PETScCommon::GetGlobalVec(Vec &x,Vec &xGlobal,bool master){
  VecScatter ctx;
  //Collect the solution vector from different procs
  master?ierr=VecScatterCreateToZero(x,&ctx,&xGlobal):VecScatterCreateToAll(x,&ctx,&xGlobal);CHKERRXX(ierr);
  ierr=VecScatterBegin(ctx,x,xGlobal,INSERT_VALUES,SCATTER_FORWARD);CHKERRXX(ierr);
  ierr=VecScatterEnd(ctx,x,xGlobal,INSERT_VALUES,SCATTER_FORWARD);CHKERRXX(ierr);
  ierr=VecScatterDestroy(&ctx);CHKERRXX(ierr);
}
void PETScCommon::SetupMGSolver(DM &da_nodes,PC &precond_){
  DM  *da_list,*daclist;
  Mat R;
  //create DMs for all levels of the multigrid so that Petsc can construct coarse grids
  PetscMalloc(sizeof(DM)*nlvls,&da_list);
  PetscMalloc(sizeof(DM)*nlvls,&daclist);
  for (PetscInt k=0; k<nlvls; k++){
    da_list[k] = NULL;
    daclist[k] = NULL;
  }

  //Set the finest nodal mesh as the first item in the da list
  daclist[0]=da_nodes;
  // petsc constructs DM structures for all other grid level
  DMCoarsenHierarchy(da_nodes,nlvls-1,&daclist[1]);
  for (PetscInt k=0; k<nlvls; k++) {
    // NOTE: finest grid is nlevels - 1: PCMG MUST USE THIS ORDER ???
    da_list[k] = daclist[nlvls-1-k];
  }

  ierr=PCMGSetLevels(precond_,nlvls,NULL);CHKERRXX(ierr);
  ierr=PCMGSetType(precond_,PC_MG_FULL); CHKERRXX(ierr);// Default
  ierr = PCMGSetCycleType(precond_,PC_MG_CYCLE_V);CHKERRXX(ierr);
  ierr=PCMGSetGalerkin(precond_,PC_MG_GALERKIN_BOTH);CHKERRXX(ierr);

  for (PetscInt k=1; k<nlvls; k++) {
    DMCreateInterpolation(da_list[k-1],da_list[k],&R,NULL);CHKERRXX(ierr);
    PCMGSetInterpolation(precond_,k,R);CHKERRXX(ierr);
    MatDestroy(&R);CHKERRXX(ierr);
  }
  // Now all DM data structure other than the finest level can be destroyed
  for (PetscInt k=1; k<nlvls; k++) { // DO NOT DESTROY LEVEL 0
    DMDestroy(&daclist[k]);CHKERRXX(ierr);
  }

  PetscFree(da_list);CHKERRXX(ierr);
  PetscFree(daclist);CHKERRXX(ierr);

  KSP cksp;
  PCMGGetCoarseSolve(precond_,&cksp);
  // The solver
  ierr = KSPSetType(cksp,KSPGMRES);

  ierr=KSPSetTolerances(cksp,coarse_rtol,coarse_atol,coarse_dtol,coarse_maxits);

  PC cpc;
  KSPGetPC(cksp,&cpc);
  PCSetType(cpc,PCSOR);

  for (PetscInt k=1;k<nlvls;k++){
    KSP dksp;
    PCMGGetSmoother(precond_,k,&dksp);
    PC dpc;
    KSPGetPC(dksp,&dpc);
    ierr = KSPSetType(dksp,innerSovler.c_str()); // KSPCG, KSPGMRES, KSPCHEBYSHEV (VERY GOOD FOR SPD)

    ierr = KSPSetTolerances(dksp,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,smooth_sweeps); // NOTE in the above maxitr=restart;
    PCSetType(dpc,PCJACOBI);// PCJACOBI, PCSOR for KSPCHEBYSHEV very good
  }
}

void PETScCommon::CheckLevels(int nx,int ny,int nz){
 //nx, ny, nz are node numbers in x,y and z direction.



  PetscScalar divisor = PetscPowScalar(2.0,(PetscScalar)nlvls-1.0);
  if ( std::floor((PetscScalar)(nx-1)/divisor) != (nx-1.0)/((PetscInt)divisor)
      || std::floor((PetscScalar)(ny-1)/divisor) != (ny-1.0)/((PetscInt)divisor)
          || std::floor((PetscScalar)(nz-1)/divisor) != (nz-1.0)/((PetscInt)divisor)) {
    nlvls=nlvls-1;
    CheckLevels(nx,ny,nz);
  }
  else
    return;
}


void PETScCommon::GetGridInfoMG(int &nx,int &ny,int &nz,int &dimension){
  //Get a grid pointer so that useful informations can be extracted
  Grid * grid=domain->GetGrid();
  ParamNodeList regionList;
  //Get some param pointers for the domain and region list so that we can iterate through all regions and sum up the number of elements in x,y and z
  PtrParamNode paramNode=domain->GetParamRoot()->Get("domain");


  regionList = paramNode->GetList("regionList");

  StdVector<RegionIdType> regionIds;

  for(unsigned int i=0; i < regionList.GetSize(); i++){
    // we are compatible with the region attribute and unbounded region elements
    string reg = regionList[i]->Get("region")->Has("name") ? regionList[i]->Get("region")->Get("name")->As<std::string>() : regionList[i]->As<std::string>();
    regionIds.Push_back(grid->GetRegionId(reg));
  }//loop over all regions

  //Number of nodes in x y and z directions respectively in case of structured with hex elem it is always elemNum in one direction +1
   nx=1,ny=1,nz=1;
  for (unsigned int j=0;j<regionIds.GetSize();j++){
    StdVector<unsigned int> dim =grid->GetRegularDiscretization(regionIds[j]); //returns number of elements in each direction
    nx+=dim[0];
    ny+=dim[1];
    nz+=dim[2];
  }



  if(domain->HasPerdiodicBC()){
  // Boundary types: DMDA_BOUNDARY_NONE, DMDA_BOUNDARY_GHOSTED, DMDA_BOUNDARY_PERIODIC
    EXCEPTION("Multigrid Doesn't Work for periodic domain yet");
  }

  paramNode=domain->GetParamRoot()->Get("sequenceStep")->Get("linearSystems")->Get("system")
     ->Get("solutionStrategy")->Get("standard");

  if (paramNode->Get("matrix/reordering")->As<std::string>() !="noReordering")
   EXCEPTION("Set matrix reordering to noReordering");


  dimension=grid->GetDim();


}


void PETScCommon::CreateDMDA(DM & daNodes,Mat &sysMat,Vec &solVec,Vec &rhsVec,Vec &dirVec,int nx,int ny,int nz,int dim){

  // number of nodal dofs
  PetscInt numnodaldof = 3;

  //std::cout<<"Assuming the nodal Degrees of freedom is 3 ,if something else not implemented yet for Multigrid"<<std::endl;

  // Stencil width: each node connects to a box around it - linear elements
  PetscInt stencilwidth = 1;


  if (dim==3){
    ierr = DMDACreate3d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_BOX,nx,ny,nz,1,1,size_,
    numnodaldof,stencilwidth,0,0,0,&(daNodes));
  }
  else{
    EXCEPTION("Not yet implemented for 2d or 1d");
  }
  //Setup the DM which is used for Grid Management
  DMSetUp(daNodes);
  DMSetMatrixPreallocateOnly(daNodes,PETSC_TRUE);
  DMCreateMatrix(daNodes,&sysMat);
  ierr=MatSetOption (sysMat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);CHKERRXX(ierr);

  //  DMDASetElementType(daNodes, DMDA_ELEMENT_Q1);
  ierr = DMCreateGlobalVector(daNodes,&(solVec));
  VecDuplicate(solVec,&(rhsVec));
  PetscInt sizeGloabalVec=0;
  VecGetSize(solVec,&sizeGloabalVec);
  VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE,sizeGloabalVec/PetscInt(3),&dirVec);
  VecSet(dirVec,1.0);
//  ierr=MatSetOption (sysMat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);CHKERRXX(ierr);

 }


void PETScCommon::GetCFSEqnMapMG(StdVector<unsigned int> &cfsEqnMap, bool& MG_FLAG){
  //This function also sets the assemble flag to not assemble in cfs.
  // Calling the grid class a lot might affect performance
  //check here for performance issues


  Grid * grid =domain->GetGrid();
  assemble_=domain->GetBasePDE()->GetAssemble();

  // if grid is regular we proceed else we switch off the multigrid assembly and proceed to CFS assembly
  if(!grid->IsGridRegular())
    MG_FLAG = false ;
  else
    assemble_->SkipElemAssembly() ;



  ParamNodeList regionList;
  //Get some param pointers for the domain and region list so that we can iterate through all regions and sum up the number of elements in x,y and z
  PtrParamNode paramNode=domain->GetParamRoot()->Get("domain");
  regionList = paramNode->GetList("regionList");




  for(unsigned int i=0; i < regionList.GetSize(); i++)
  {
  // we are compatible with the region attribute and unbounded region elements
    string reg = regionList[i]->Get("region")->Has("name") ? regionList[i]->Get("region")->Get("name")->As<std::string>() : regionList[i]->As<std::string>();
    shared_ptr<EntityList> nodesList = grid->GetEntityList(EntityList::NODE_LIST, reg);
    EntityIterator nodeIt=nodesList->GetIterator();
    for (nodeIt.Begin(); !nodeIt.IsEnd(); nodeIt++)
    {
     StdVector<int> Eqns;
     StdVector<SinglePDE *> PDEs=domain->GetSinglePDEs();
     for(unsigned int i = 0; i < PDEs.GetSize(); i++)
     {
       SolutionType pde_solution=PDEs[i]->GetNativeSolutionType();
       shared_ptr<BaseFeFunction> feFunction = PDEs[i]->GetFeFunction(pde_solution);
       feFunction->GetFeSpace()->GetEqns(Eqns,nodeIt);
//      std::cout<<nodeIt.GetNode()<<std::endl;
//      std::cout<<Eqns.ToString()<<std::endl;
       for (unsigned int ldof=0;ldof<Eqns.GetSize();ldof++)
       {
         // cfsEqnMap is zero based and will be used to map
         cfsEqnMap.Push_back(Eqns[ldof]) ;
        } //loop over dofs of each node
      }//loop over all pdes
    }//loop over all nodes
  }//loop over all regions
}



void PETScCommon::AssembleMatrixMG(Mat &sysMat,DM &daNodes,Vec &solVec,Vec &rhsVec,Vec &dirNodeVec,Vec &dirVec,
    int nx,int ny,int nz)
{
  //This method gets the local element matrix and assembels into global matrix
  PetscInt edof[24]; // connctinvity of element with eqnNr for assembly.

  //Set the Dirchlet Vector With the homogeneous boundary conditions
  PetscScalar * dirArray;
  VecGetArray(dirNodeVec,&dirArray);

  //Number of element nodes
  int nen=0;


  //One based in CFS , but we will use 0 based in petsc
  StdVector<int> globalNodeNum;

  Vector<double>elemVector;
  PetscScalar KE[24*24]; // Element stiffness matrix

   // Loop over all bilinear forms and calculate the element stiffness matrix for all regions, in this case since its a structured grid one
   // element stiffness matrix should be sufficient.
  int elem=0;
  for (const auto& [_, contexts] : assemble_->GetBiLinForms()) {
    for(BiLinFormContext* context : contexts) {
      BiLinearForm* form = context->GetIntegrator();

      EntityIterator firstEntIt = context->GetFirstEntities()->GetIterator();
      EntityIterator secondEntIt = context->GetSecondEntities()->GetIterator();

      UInt numElems=std::max(firstEntIt.GetSize(), secondEntIt.GetSize());
      firstEntIt.Begin();
      secondEntIt.Begin();
      UInt chunksize;

      disableParalleAssemby ? chunksize= numElems:chunksize= std::floor(numElems/size_);

      UInt start=chunksize*rank_;
      UInt end=(rank_==size_-1)? numElems:(chunksize*(rank_+1));
      firstEntIt+=start;
      secondEntIt+=start;
      for (UInt i=start;i<end;++i)
      {
        CacluateElementStiffnessMatrix(firstEntIt,secondEntIt,globalNodeNum,nen,form,elem,KE);
        CalculateDirNodesAndEdof(dirVec,edof,nen,dirArray,globalNodeNum,elem);
        ierr = MatSetValues(sysMat,24,edof,24,edof,KE,ADD_VALUES); CHKERRXX(ierr);
        firstEntIt++;
        secondEntIt++;
      }
    }
  }
}

void PETScCommon::CalculateDirNodesAndEdof(Vec &dirVec,PetscInt edof[],const int &nen,const PetscScalar * dirArray,
    const StdVector<int> &globalNodeNum,const int &elem){

  // loop over element nodes
  for (PetscInt elem_node=0;elem_node<nen;elem_node++){
    // Get local dofs
    for (PetscInt dof=0;dof <3;dof++){
      // element local idx 0 ... 23
      int idx = elem_node*3+dof;
      // global idx
      int gidx = globalNodeNum[elem*nen+elem_node]-1;
      int eqnNr = 3*(globalNodeNum[elem*nen+elem_node]-1)+dof; //CFS node number is one  based indexing ,here its zero based
      if (dirArray[gidx] == -1.0 ){
        VecSetValue(dirVec,(PetscInt)eqnNr,0,INSERT_VALUES);
      }
      edof[idx]=eqnNr;
  //          edof[elem_node*3+dof] = 3*necon[elem*nen+elem_node]+dof;
    }//dof Loop
  }//Element Node Loop

}

void PETScCommon::CacluateElementStiffnessMatrix(EntityIterator &firstEntIt, EntityIterator &secondEntIt,StdVector<int> &globalNodeNum
    ,int &nen , BiLinearForm  *form ,int &elem , PetscScalar KE[] ){

  StdVector<unsigned int> nodesInElem=firstEntIt.GetElem()->connect;

  //Extract the element node connectivity information which is useful for assembly of Stiffness matrix
  // Everything here is zero based and the globalNodeNum array consist of all the nodes of each element starting
  // from zero to N. It is also assumed that the number of nodes per element is 8

  elem=(firstEntIt.GetElem()->elemNum)-1;
  nen=nodesInElem.GetSize();
  //Matrix to store the element stiffness matrix write a method to convert the Matrix to Petsc Matrix
  Matrix<double>elemMatrix;
  for (unsigned int i =0;i<nodesInElem.GetSize();i++)
    globalNodeNum.Push_back(nodesInElem[i]);

  form->CalcElementMatrix(elemMatrix, firstEntIt, secondEntIt); // use the region part

  //copy the element matrix to the PetscVector
  for( UInt i=0; i <elemMatrix.GetNumRows(); i++)
    for( UInt j=0; j < elemMatrix.GetNumCols(); j++)
      elemMatrix.GetEntry(i,j,KE[i*(elemMatrix.GetNumCols()) + j]);

}


void PETScCommon::SetLinRhs(Vec &rhsVec){

 //Get the RHS

 Vector<Double> elemVec;
 StdVector<LinearFormContext *> linForms_=assemble_->GetLinForms();
 StdVector<LinearFormContext*>::iterator formsIt;
 // iterate over all descriptors
 for(formsIt = linForms_.Begin(); formsIt != linForms_.End(); formsIt++)
 {
  // get integrator
  LinearFormContext& actContext = **formsIt;
  LinearForm* form = actContext.GetIntegrator();
    // get entity iterator
    EntityIterator  entIt = actContext.GetEntities()->GetIterator();

      // That should be STATIC
      Vector<Double> elemVec;
      // iterate over all entities
      for ( entIt.Begin(); !entIt.IsEnd(); entIt++ ) {
        StdVector<Integer> eqnVec;
        // Calculate real valued element vector

        form->CalcElemVector(elemVec, entIt);
        //Find the equNr in Petsc Format for the given element
          for (int dof=0;dof<3;dof++){
            int eqnNr=3* (entIt.GetNode()-1)+dof;
            eqnVec.Push_back(eqnNr);
          }
        assert(!elemVec.ContainsNaN() && !elemVec.ContainsInf());
//        std::cout<<elemVec.ToString()<<std::endl;
//        std::cout<<eqnVec.ToString()<<std::endl;
        VecSetValues(rhsVec,eqnVec.GetSize(),eqnVec.GetPointer(),elemVec.GetPointer(),ADD_VALUES);
    }

  }

}


void  PETScCommon::GetHomDirNodes(Vec &dirVec){

  //Homogeneous Dirchlet Boundary conditions

 //Get Node Number where the homogeneous dirichlet boundary conditions apply
 StdVector<int> DirchletNodeVec;

 //Set HomDirchletValue to -1 for all nodes with Dirchlet Boundary conditions
 StdVector<double> HomDirchletValue;
 //Get all the PDEs to loop over , in this case it mostly should be one
 StdVector<SinglePDE *> PDEs=domain->GetSinglePDEs();

 for(unsigned int i = 0; i < PDEs.GetSize(); i++){
   SolutionType pde_solution=PDEs[i]->GetNativeSolutionType();
   shared_ptr<BaseFeFunction> feFunction = PDEs[i]->GetFeFunction(pde_solution);
   //Get the list of hdbc from the fe_function which returns a list of structs which are of type HomDirichletBc
   HdBcList hbc_list=feFunction->GetHomDirichletBCs();
   //Get a vector of all the nodes for which hdbc is prescribed
   StdVector<unsigned int> hbcNodes;
//      std::cout<<"hbcListSize "<<hbc_list.GetSize()<<std::endl;
   for (unsigned int i=0;i<hbc_list.GetSize();i++){
     HomDirichletBc hdbc= *hbc_list[i];
//        std::cout<<"hdbcNodeSize"<<hdbc.entities.get()->ToString()<<std::endl;
     EntityIterator hdbcEntIt=hdbc.entities.get()->GetIterator();

     // iterate over all entities
     for ( hdbcEntIt.Begin(); !hdbcEntIt.IsEnd(); hdbcEntIt++ ) {
       // loop over all dofs
       DirchletNodeVec.Push_back(hdbcEntIt.GetNode()-1);
       HomDirchletValue.Push_back(-1.0);
//            std::cout<<eqnVec.ToString()<<std::endl;
     } // loop: dofs

   } // loop over hdbcs

 }//loop over pdes

//    std::cout<<"HomDirchlet Nodes"<<DirchletNodeVec.ToString()<<std::endl;
 VecSetValues(dirVec,DirchletNodeVec.GetSize(),DirchletNodeVec.GetPointer(),HomDirchletValue.GetPointer(),INSERT_VALUES);
}

std::string PETScCommon::CreatePrecondString(PtrParamNode xml){
  string output;
  PtrParamNode sNode = xml->Get("precond",ParamNode::PASS);
  ParamNodeList sol = sNode->GetChildren();
  UInt numChilds = sol.GetSize();
  if(numChilds==2){
    output = sol[1]->GetName();
    if (output=="mg" or output=="gamg"){
      nlvls=sol[1]->Get("nlvls")->As<int>();
      coarse_atol=sol[1]->Get("coarse_atol")->As<double>();
      coarse_rtol=sol[1]->Get("coarse_rtol")->As<double>();
      coarse_dtol=sol[1]->Get("coarse_dtol")->As<double>();
      coarse_maxits=sol[1]->Get("coarse_maxits")->As<int>();
      innerSovler=sol[1]->Get("innerSolver")->As<string>();
    }
  }
  else{
    output="jacobi";
  }
  return output;
}


std::string PETScCommon::CreateSolverString(PtrParamNode xml){
  string output;
  PtrParamNode sNode = xml->Get("solver",ParamNode::PASS);
  ParamNodeList sol = sNode->GetChildren();

  UInt numChilds = sol.GetSize();
  if(numChilds==2){
    output = sol[1]->GetName();
  }
  else{
    output="cg";
  }
  return output;
}
} /* namespace CoupledField */
