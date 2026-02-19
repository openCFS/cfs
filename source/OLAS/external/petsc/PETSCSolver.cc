#include "PETSCSolver.hh"

#include <string>
#include <sstream>
#include <stdio.h>
#include <mpi.h>
#include "MatVec/SparseOLASMatrix.hh"
#include "MatVec/SCRS_Matrix.hh"
#include "MatVec/CRS_Matrix.hh"
#include "DataInOut/ProgramOptions.hh"
#include "DataInOut/ParamHandling/XmlReader.hh"
#include "Utils/Timer.hh"



//Checks Petsc error code, if non-zero it calls the C++ error handler which throws an exception
#define CHKERRXX(ierr)  do {if (PetscUnlikely(ierr)) {PetscError(PETSC_COMM_SELF,__LINE__,PETSC_FUNCTION_NAME,__FILE__,ierr,PETSC_ERROR_IN_CXX,0);}} while(0)

using std::string;

namespace CoupledField{

DEFINE_LOG(petsc, "petscSolver")
//initialize PETSCSolver class
PETSCSolver::PETSCSolver(PtrParamNode pn, PtrParamNode olasInfo, BaseMatrix::EntryType type){
  A_ = nullptr;
  da_nodes = nullptr;
  x_ = nullptr;
  b_ = nullptr;
  precond_ = nullptr;
  solver_ = nullptr;
  firstSetup_ = true;
  dirNodeVec_ = nullptr;
  N_ = nullptr;
  xml_ = pn;
  infoNode_ =  olasInfo->Get("petsc");

  maxIter_    = xml_->Get("maxIter")->As<int>();
  tolerance_  =  xml_->Get("tolerance")->As<double>();
  minTol_  =  xml_->Get("minimalTolerance")->As<double>();

  solverstring_ = xml_->Get("solver")->As<std::string>();
  precondstring_ = CreatePrecondString(xml_); //this fn also gets some values to be set for the multigrid if precond=mg


  PtrParamNode hdr = infoNode_->Get(ParamNode::HEADER);
  hdr->Get("maxIter")->SetValue(maxIter_);
  hdr->Get("tolerance")->SetValue(double(tolerance_));
  hdr->Get("solver")->SetValue(solverstring_);
  hdr->Get("precond")->SetValue(precondstring_);
  if (precondstring_=="mg" ||precondstring_=="gamg" ){
     MG_FLAG=true;
     // This fn also sets the assemble class flag skipElemAssembly_ ,
     // also returns MG_FLAG false in case of unstructured meshes
     GetCFSEqnMapMG(cfsEqnMap_,MG_FLAG);
     if (!MG_FLAG && precondstring_=="mg")
       Exception("Geometric Multigrid is implemented only for regular grids");
     hdr->Get("innerSolver")->SetValue(innerSovler);
     hdr->Get("coarse_maxits")->SetValue(coarse_maxits);
   }

//  std::cout<<"MG Flag" << MG_FLAG <<std::endl;


}

	//Destructor for solver
PETSCSolver::~PETSCSolver(){

  VecDestroy(&x_);
  VecDestroy(&b_);
  MatDestroy(&A_);
  KSPDestroy(&solver_);
  VecDestroy(&dirNodeVecGlobal_);
  VecDestroy(&dirNodeVec_);
  VecDestroy(&N_);
  cfsEqnMap_.Clear(false);
  DMDestroy(&da_nodes);

}

	

void PETSCSolver::Setup(BaseMatrix &sysmat){


  //create petsc matrix form the sysmatrix
  const StdMatrix& stdmat = static_cast<const StdMatrix&>(sysmat);

  BaseMatrix::EntryType etype = stdmat.GetEntryType(); //currently only real values can be solved using Petsc should implement for complex values

  BaseMatrix::StorageType stype = stdmat.GetStorageType();

  if (stype==BaseMatrix::SPARSE_SYM)
    symmetric=true;


  const SCRS_Matrix<Double>& crs = static_cast<const SCRS_Matrix<Double>&>(sysmat);

  if(crs.GetNumCols() != crs.GetNumRows())
    EXCEPTION("PETSC solver only tested for quadratic matrices");

  if(!MG_FLAG && precondstring_ =="gamg" && symmetric)
    EXCEPTION("Algebraic Multigrid works only when matrix storage is sparseNonSym with non-regular grids ")


  //gather info
  UInt dim = 0;
  PetscInt* rowPtr = nullptr;
  PetscInt* colPtr = nullptr;
  PetscScalar * dataPtr = nullptr;
  dim = crs.GetNumRows();
  rowPtr = (PetscInt*) crs.GetRowPointer();
  colPtr = (PetscInt*) crs.GetColPointer();
  dataPtr = const_cast<PetscScalar *>(crs.GetDataPointer());

  int nx=0,ny=0,nz=0,dimension=0;
  if (firstSetup_ ){

    SendWorkerCommand(SETUP_MATRIX);

    // The master determines if it is okay for multigrid to be used and communicate to workers
    for (int rank=1;rank<size_;rank++)
      ierr=MPI_Send(&MG_FLAG,sizeof(bool),MPIU_BOOL,rank,DATA,PETSC_COMM_WORLD);CHKERRXX(ierr);


    if (MG_FLAG)
    {
      GetGridInfoMG(nx,ny,nz,dimension);
      for (int rank=1;rank<size_;rank++){
        ierr=MPI_Send(&nx,sizeof(int),MPI_INT,rank,DATA,PETSC_COMM_WORLD);CHKERRXX(ierr);
        ierr=MPI_Send(&ny,sizeof(int),MPI_INT,rank,DATA,PETSC_COMM_WORLD);CHKERRXX(ierr);
        ierr=MPI_Send(&nz,sizeof(int),MPI_INT,rank,DATA,PETSC_COMM_WORLD);CHKERRXX(ierr);
        ierr=MPI_Send(&dimension,sizeof(int),MPI_INT,rank,DATA,PETSC_COMM_WORLD);CHKERRXX(ierr);
      }
      CreateDMDA(da_nodes,A_,x_,b_,dirNodeVec_,nx,ny,nz,dimension);
      CheckLevels(nx,ny,nz);
      GetHomDirNodes(dirNodeVec_);
      SendWorkerCommand(SET_HOM_DIR_VEC);
      VecDuplicate(b_,&N_);
      VecSet(N_,1.0);
      SendWorkerCommand(GET_GLOBAL_VEC);

      VecAssemblyBegin(dirNodeVec_);
      VecAssemblyEnd(dirNodeVec_);
      GetGlobalVec(dirNodeVec_,dirNodeVecGlobal_,true); // dirNodeVecGlobal_ has entire vector in master node

    }
    else
    {

      // We are here when PETSc iterative solvers are used or when the grid is non regular


      //To preallocate without grid info we require info about number of non-zero per row
      UInt * nnzr = new UInt[dim];

      //100 is just a big enough number which ensures there is no lack of memory allocation for petsc matix in parallel
      for (UInt ii=0;ii<dim;ii++)
        nnzr[ii]=crs.GetRowSize(ii)+100;



      //only time where data is sent using mpi(can be implemented better)
      for (int rank=1;rank<size_;rank++){
        ierr=MPI_Send(&nnzr[0], sizeof(int)*dim,MPI_INT,rank,DATA,PETSC_COMM_WORLD);CHKERRXX(ierr);
        ierr=MPI_Send(&symmetric,sizeof(bool),MPIU_BOOL,rank,ISSYMMETRIC,PETSC_COMM_WORLD);CHKERRXX(ierr);
      }
      SetupMatrix(dim,nnzr,symmetric,A_,b_,x_);
    }
  }


  if (MG_FLAG){
    //Zeros all matrix entry without affecting the non-zero structure of sparse matrix
    SendWorkerCommand(MAT_ZERO_ENTRIES);
    MatZeroEntries(A_);

    shared_ptr<Timer> Assemble = infoNode_->Get(ParamNode::SUMMARY)->Get("petsc_assembly/timer")->AsTimer(GetSetupTimer());
    Assemble->Start();


    //Assembly of sys mat done only in master yet. We can switch to multiple node assembly
    SendWorkerCommand(ASSEMBLE_MAT);
    AssembleMatrixMG(A_,da_nodes,x_,b_,dirNodeVecGlobal_,N_,nx,ny,nz);
    ierr=MatAssemblyBegin(A_, MAT_FINAL_ASSEMBLY);CHKERRXX(ierr);
    ierr=MatAssemblyEnd(A_, MAT_FINAL_ASSEMBLY);CHKERRXX(ierr);
    Assemble->Stop();

    SendWorkerCommand(ASSEMBLE_VEC_DIR);
    ierr=VecAssemblyBegin(N_);CHKERRXX(ierr);
    ierr=VecAssemblyEnd(N_);CHKERRXX(ierr);

    // Set the right hand side vector with master node only.
    SetLinRhs(b_);
    SendWorkerCommand(ASSEMBLE_VEC_RHS);
    ierr=VecAssemblyBegin(b_);CHKERRXX(ierr);
    ierr=VecAssemblyEnd(b_);CHKERRXX(ierr);

    // Apply penalty to zero out the nodes.
    SendWorkerCommand(HOM_DIR_PENALTY);
    PenaltyHomDir(A_,b_,N_);

  }

  else {
    /*setting values only using the proc 0 while others wait and then asemble is called in all procs so
    the matix can be distributed.Set each matrix value one by one (highly inefficient there are better ways to do
    this which requires some tweakng to the data structure we receive from the matrix class  )*/
    for (int i=0;i<int(dim);i++)
     for (int j=rowPtr[i];j<rowPtr[i+1];j++)
       ierr=MatSetValue(A_,i,colPtr[j],dataPtr[j],INSERT_VALUES);
    //  Distribute the assembeled matrix accross all process
    SendWorkerCommand(ASSEMBLE_MAT);
    ierr=MatAssemblyBegin(A_,MAT_FINAL_ASSEMBLY);CHKERRXX(ierr);
    ierr=MatAssemblyEnd(A_,MAT_FINAL_ASSEMBLY);CHKERRXX(ierr);
  }
  if (firstSetup_){
    SendWorkerCommand(SETUP_SOLVER_CONTEXT);
    SetupSolverContext(A_,solver_,precond_,solverstring_,precondstring_,tolerance_,minTol_,maxIter_,MG_FLAG);
    if(MG_FLAG){SetupMGSolver(da_nodes,precond_);}
    infoNode_->Get(ParamNode::HEADER)->Get("mgLevels")->SetValue(nlvls);
  }

  firstSetup_=false;
}


void PETSCSolver::Solve( const BaseMatrix &sysmat, const BaseVector &rhs, BaseVector &sol){

  if (!MG_FLAG){
    SendWorkerCommand(ASSEMBLE_VEC_RHS);
    ierr= VecAssemblyBegin(b_);CHKERRXX(ierr);
    ierr= VecAssemblyEnd(b_);CHKERRXX(ierr);
    if(sysmat.GetEntryType() == BaseMatrix::DOUBLE) {
      for(PetscInt i=0, n=(PetscInt)rhs.GetSize(); i<n; i++){
        Double myEnt =0;
        rhs.GetEntry((UInt)i,myEnt);
        ierr=VecSetValue(b_,i,PetscScalar(myEnt),INSERT_VALUES);
      }
    }
    SendWorkerCommand(ASSEMBLE_VEC_RHS);
    ierr= VecAssemblyBegin(b_);CHKERRXX(ierr);
    ierr= VecAssemblyEnd(b_);CHKERRXX(ierr);
  }



  double t1,t2;
  t1 = MPI_Wtime();
  //Solve the linear system
  SendWorkerCommand(SOLVE);
  ierr = KSPSolve(solver_,b_,x_);CHKERRXX(ierr);
  PetscReal RHSnorm;
  VecNorm(b_,NORM_2,&RHSnorm);
  VecSet(b_,0.0);
  t2 = MPI_Wtime();

  ParamNode::ActionType at = progOpts->DoDetailedInfo() ? ParamNode::APPEND : ParamNode::DEFAULT;
  PtrParamNode out = infoNode_->Get(ParamNode::PROCESS)->Get("solver", at);

  PetscScalar rnorm =0.0 ;

  KSPGetIterationNumber(solver_,&niter);CHKERRXX(ierr);
  KSPGetResidualNorm(solver_,&rnorm);CHKERRXX(ierr);
  totalSolverIter += niter;
  rnorm = rnorm/RHSnorm;
  if (progOpts->DoDetailedInfo())
    out->Get("timing")->SetValue(t2-t1);
  out->Get("residualNorm")->SetValue(rnorm);
  out->Get("iterationCurrentStep")->SetValue(niter);
  out->Get("totalIterations")->SetValue(totalSolverIter+niter);


//  if( rnorm > tolerance_ && rnorm > minTol_){
//    EXCEPTION("Linear Solver has not converged, reducing the tolerance or increasing the iterations might work sometime");
//  }

  //Create a global vector and scatter context to collect the solution vector to master process
  SendWorkerCommand(GET_SOL);
  Vec x_global;
  GetGlobalVec(x_,x_global,true);


  //gets the solution vector in bufx array from petsc array
  PetscScalar *bufx;
  ierr=VecGetArray(x_global,&bufx);CHKERRXX(ierr);
  PetscInt Size;
  VecGetSize(x_global,&Size);
//		sol.Resize(Size);
//  std::cout<<cfsEqnMap.ToString()<<std::endl;

  if (MG_FLAG){
    for(unsigned int i=0; i<cfsEqnMap_.GetSize()/3;i++)
        for (unsigned int ldof=0;ldof<3;ldof++)
          if (cfsEqnMap_[(i*3)+ldof]!=0)
            sol.SetEntry(cfsEqnMap_[(i*3)+ldof]-1,bufx[PetscInt(i*3)+ldof]);
  }
  else {
    for(UInt i=0; i<sol.GetSize();i++)
     sol.SetEntry(i,bufx[PetscInt(i)]);
  }

//		std::cout<<"RHS from CFS"<<rhs.ToString()<<std::endl;
//		std::cout<<"Solution Vector in PETSc" <<sol.ToString()<<std::endl;
  ierr=VecRestoreArray(x_global,&bufx);CHKERRXX(ierr);



}

		
void PETSCSolver::SendWorkerCommand(int TAG){

  for (int rank = 1; rank < size_; rank++) {
    ierr=MPI_Send(0, 0, MPI_INT, rank, TAG, MPI_COMM_WORLD);CHKERRXX(ierr);
  }
}

//Methods for PETSCWorker class
PETSCWorker::PETSCWorker(int argc,const char **argv){
  //init petsc objects
  A_ =nullptr;
  b_=nullptr;
  solver_=nullptr;
  x_=nullptr;
  precond_=nullptr;
  dirNodeVec_=nullptr;
  da_nodes=nullptr;
  N_=nullptr;
  progOpts = new ProgramOptions(argc, argv);

  // Parse command line
  progOpts->ParseData();

  string xmlFile = progOpts->GetParamFileStr();
  string schema = progOpts->GetSchemaPathStr() + "/CFS-Simulation/CFS.xsd";

  PtrParamNode root_node = XmlReader::ParseFile(xmlFile, schema, "http://www.cfs++.org/simulation");

  root_node = root_node->Has("sequenceStep") ? root_node->Get("sequenceStep"): nullptr;

  if (root_node != nullptr)
  root_node = root_node->Has("linearSystems") ? root_node->Get("linearSystems") : nullptr;

  if (root_node != nullptr)
    xml_=root_node->Has("system") ? (root_node->Get("system")->Has("solverList") ? root_node->Get("system")->Get("solverList"): nullptr):nullptr;



  if(xml_!= nullptr && xml_->Has("petsc")){

    xml_=xml_->Get("petsc");
    maxIter_    = xml_->Get("maxIter")->As<int>();
    tolerance_  = xml_->Get("tolerance")->As<double>();
    minTol_  = xml_->Get("minimalTolerance")->As<double>();

    solverstring_=xml_->Get("solver")->As<std::string>();
    precondstring_=CreatePrecondString(xml_);
  }


}



PETSCWorker::~PETSCWorker(){

  VecDestroy(&x_);
  VecDestroy(&b_);

  MatDestroy(&A_);

  KSPDestroy(&solver_);

  VecDestroy(&dirNodeVecGlobal_);
  VecDestroy(&dirNodeVec_);
  VecDestroy(&N_);
  DMDestroy(&da_nodes);

  //delete the PtrParamNode  and progOpts of worker
  delete progOpts;
  progOpts = NULL;
  xml_.reset();

}



void PETSCWorker::run(){



  bool done=false;
  MPI_Status stat;
  int work;
  //Worker loop which runs till a kill command is received
  // MPI communication for starting some process in worker
  while (!done){
    MPI_Recv(&work,1,MPI_INT,0,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
    switch(stat.MPI_TAG){
      case DIETAG:
        done=true;
        break;
      case SET_HOM_DIR_VEC:
        VecDuplicate(b_,&N_);
        VecSet(N_,1.0);
        break;
      case SETUP_MATRIX:
        InitPetscWorker();
        break;
      case ASSEMBLE_MAT:
        if (!disableParalleAssemby)
          AssembleMatrixMG(A_,da_nodes,x_,b_,dirNodeVecGlobal_,N_,nx,ny,nz);
        ierr=	MatAssemblyBegin(A_,MAT_FINAL_ASSEMBLY);
        ierr=	MatAssemblyEnd(A_,MAT_FINAL_ASSEMBLY);
        break;
      case SETUP_SOLVER_CONTEXT:
        SetupSolverContext(A_,solver_,precond_,solverstring_,precondstring_,tolerance_,minTol_,maxIter_,MG_FLAG);
        if(MG_FLAG){SetupMGSolver(da_nodes,precond_);}
        break;
      case ASSEMBLE_VEC_RHS:
        ierr=	VecAssemblyBegin(b_);
        ierr=	VecAssemblyEnd(b_);
        break;
      case ASSEMBLE_VEC_DIR:
        ierr= VecAssemblyBegin(N_);
        ierr= VecAssemblyEnd(N_);
        break;
      case SOLVE:
        ierr = KSPSolve(solver_,b_,x_);
        PetscReal RHSnorm;
        VecNorm(b_,NORM_2,&RHSnorm);
        VecSet(b_,0.0);
        break;
      case GET_SOL:
        Vec x_global;
        GetGlobalVec(x_,x_global,true);
        break;
      case GET_GLOBAL_VEC:
        VecAssemblyBegin(dirNodeVec_);
        VecAssemblyEnd(dirNodeVec_);
        GetGlobalVec(dirNodeVec_,dirNodeVecGlobal_,true);
        break;
      case HOM_DIR_PENALTY:
        PenaltyHomDir(A_,b_,N_);
        break;
      case MAT_ZERO_ENTRIES:
        MatZeroEntries(A_);
        break;
    }
  }
}


void PETSCWorker::InitPetscWorker(){

  MPI_Status status;

  MPI_Recv(&MG_FLAG,1,MPIU_BOOL,0,DATA,PETSC_COMM_WORLD,&status);

//  std::cout<<"MG Flag" << MG_FLAG <<std::endl;

  if (MG_FLAG){

    MPI_Recv(&nx,1,MPI_INT,0,DATA,PETSC_COMM_WORLD,&status);
    MPI_Recv(&ny,1,MPI_INT,0,DATA,PETSC_COMM_WORLD,&status);
    MPI_Recv(&nz,1,MPI_INT,0,DATA,PETSC_COMM_WORLD,&status);
    MPI_Recv(&dimension,1,MPI_INT,0,DATA,PETSC_COMM_WORLD,&status);
    CreateDMDA(da_nodes,A_,x_,b_,dirNodeVec_,nx,ny,nz,dimension);
    CheckLevels(nx,ny,nz);
//    assemble_=domain->GetBasePDE()->GetAssemble();
  }
  else {
    int size;
    MPI_Probe(0,DATA,PETSC_COMM_WORLD,&status);
    MPI_Get_count(&status,MPI_INT,&size);
    int dim=size/sizeof(int);
    UInt *nnzr = new UInt[dim];
    //Receive the nnzr array
    MPI_Recv(nnzr,dim,MPI_INT,0,DATA,PETSC_COMM_WORLD,&status);
    MPI_Recv(&symmetric,1,MPIU_BOOL,0,ISSYMMETRIC,PETSC_COMM_WORLD,&status);
    SetupMatrix(dim,nnzr,symmetric,A_,b_, x_);

  }

}




}	


