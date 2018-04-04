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

DECLARE_LOG(petsc)
DEFINE_LOG(petsc, "petscSolver")
//initialize PETSCSolver class
PETSCSolver::PETSCSolver(PtrParamNode pn, PtrParamNode olasInfo, BaseMatrix::EntryType type){
  A_=nullptr;
  da_nodes=nullptr;
  x_=nullptr;
  b_=nullptr;
  precond_=nullptr;
  solver_=nullptr;
  firstSetup_ = true;
  dirNodeVec_=nullptr;
  N_=nullptr;

  infoNode_ =  olasInfo->Get("petsc");
  xml_ = pn;
  maxIter_    = xml_->Has("maxIter") ? pn->Get("maxIter")->As<int>() : 10000;

  tolerance_  = xml_->Has("tolerance") ? pn->Get("tolerance")->As<double>() : 1e-12;
  minTol_     = xml_->Has("minimalTolerance") ? pn->Get("minimalTolerance")->As<double>() : 1e-11;
  logging_    = xml_->Has("logging") ? pn->Get("logging")->As<bool>() : false;

  PtrParamNode hdr = infoNode_->Get(ParamNode::HEADER);
  hdr->Get("maxIter")->SetValue(maxIter_);
  hdr->Get("tolerance")->SetValue(double(tolerance_));
  hdr->Get("minimalTolerance")->SetValue(double(minTol_));

  solverstring_=CreateSolverString(xml_);
  precondstring_=CreatePrecondString(xml_); //this fn also gets some values to be set for the multigrid if precond=mg
  if (precondstring_=="mg"){
    MG_FLAG=true;
    GetCFSEqnMapMG(cfsEqnMap_); //This fn also sets the assemble class flag skipElemAssembly_.
  }




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
  assert(etype==BaseMatrix::DOUBLE);

  BaseMatrix::StorageType stype = stdmat.GetStorageType();

  if (stype==BaseMatrix::SPARSE_SYM){symmetric=true;}

  //Requires the grid information ie dimension of the mesh, global_coordinates, boundary_type,stencil_type,ndof,stencil_width

  const SCRS_Matrix<Double>& crs = static_cast<const SCRS_Matrix<Double>&>(sysmat);

  if(crs.GetNumCols() != crs.GetNumRows()){
  EXCEPTION("PETSC solver only tested for quadratic matrices");
  }
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
    if (MG_FLAG){


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
      GetGlobalVec(dirNodeVec_,dirNodeVecGlobal_); // dirNodeVecGlobal_ has entire vector in master node

    }
    else {

      //To preallocate without grid info we require info about number of non-zero per row
      UInt * nnzr = new UInt[dim];
      for (UInt ii=0;ii<dim;ii++){
      //100 is just a big enough number which ensures there is no lack of memory allocation for petsc matix in parallel
       nnzr[ii]=crs.GetRowSize(ii)+100;
      }
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

    shared_ptr<Timer> Assemble = infoNode_->Get(ParamNode::SUMMARY)->Get("petsc_assembly/timer")->AsTimer();
    Assemble->SetSub();
    Assemble->Start();


    //Assembly of sys mat done only in master yet. We can switch to multiple node assembly
    AssembleMatrixMG(A_,da_nodes,x_,b_,dirNodeVecGlobal_,N_,nx,ny,nz);
    SendWorkerCommand(ASSEMBLE_MAT);
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
    for (int i=0;i<int(dim);i++){
     for (int j=rowPtr[i];j<rowPtr[i+1];j++){
       ierr=MatSetValue(A_,i,colPtr[j],dataPtr[j],INSERT_VALUES);
     }
    }
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
  // Log the the norm and max iter
  PtrParamNode curr = infoNode_->Get(ParamNode::PROCESS)->Get("solve", ParamNode::APPEND);
  PetscInt niter=0;
  PetscScalar rnorm=0.0;

  KSPGetIterationNumber(solver_,&niter);CHKERRXX(ierr);
  KSPGetResidualNorm(solver_,&rnorm);CHKERRXX(ierr);

  rnorm = rnorm/RHSnorm;
  curr->Get("timing")->SetValue(t2-t1);
  curr->Get("residualNorm")->SetValue(rnorm);
  curr->Get("iterations")->SetValue(niter);
  //Create a global vector and scatter context to collect the solution vector to master process
  SendWorkerCommand(GET_SOL);
  Vec x_global;
  GetGlobalVec(x_,x_global);


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
  string schema = progOpts->GetSchemaPathStr();

  schema += "/CFS-Simulation/CFS.xsd";

  xml_ = XmlReader::ParseFile(xmlFile, schema)->Get("sequenceStep")->Get("linearSystems")->Get("system")->Get("solverList");
  if(xml_->GetChild()->GetName()=="petsc"){

    xml_=xml_->Get("petsc");
    maxIter_    = xml_->Has("maxIter") ? xml_->Get("maxIter")->As<int>() : 10000;
    tolerance_  = xml_->Has("tolerance") ? xml_->Get("tolerance")->As<double>() : 1e-12;
    minTol_     = xml_->Has("minimalTolerance") ? xml_->Get("minimalTolerance")->As<double>() : 1e-11;
    solverstring_=CreateSolverString(xml_);
    precondstring_=CreatePrecondString(xml_);
    if (precondstring_=="mg")
      MG_FLAG=true;
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
        GetGlobalVec(x_,x_global);
        break;
      case GET_GLOBAL_VEC:
        VecAssemblyBegin(dirNodeVec_);
        VecAssemblyEnd(dirNodeVec_);
        Vec dirNodeVecGlobal_;
        GetGlobalVec(dirNodeVec_,dirNodeVecGlobal_);
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
  if (MG_FLAG){
    int nx,ny,nz,dimension;
    MPI_Recv(&nx,1,MPI_INT,0,DATA,PETSC_COMM_WORLD,&status);
    MPI_Recv(&ny,1,MPI_INT,0,DATA,PETSC_COMM_WORLD,&status);
    MPI_Recv(&nz,1,MPI_INT,0,DATA,PETSC_COMM_WORLD,&status);
    MPI_Recv(&dimension,1,MPI_INT,0,DATA,PETSC_COMM_WORLD,&status);
    CreateDMDA(da_nodes,A_,x_,b_,dirNodeVec_,nx,ny,nz,dimension);
    if (MG_FLAG){CheckLevels(nx,ny,nz);}
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


