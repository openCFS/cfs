#include "PETSCSolver.hh"

#include <string>
#include "MatVec/SparseOLASMatrix.hh"
#include "MatVec/SCRS_Matrix.hh"
#include "MatVec/CRS_Matrix.hh"
#include "DataInOut/ProgramOptions.hh"
#include <sstream>
#include <stdio.h>
#include<mpi.h>

#define INIT_MAT_STRUCT 0
#define ASSEMBLE_MAT 1
#define ASSEMBLE_VEC_RHS  2
#define SETUP_MATRIX 3
#define SOLVE 4
#define SEND_NNZR 5
#define GET_SOL 6
#define DIETAG 7


using std::string;

namespace CoupledField{
	//initialize PETSCSolver class 	
	PETSCSolver::PETSCSolver(PtrParamNode pn, PtrParamNode olasInfo, BaseMatrix::EntryType type){

		//init petsc objects
		A_ =NULL;
		A0_=NULL;
		b_=NULL;
		solver_=NULL;

		
		infoNode_ =  olasInfo->Get("petsc");
		
		xml_ = pn;
		firstSetup_ = true;
		ownMatrixA_ = false;
		
		maxIter_    = pn->Has("maxIter") ? pn->Get("maxIter")->As<int>() : 10000;
		tolerance_  = pn->Has("tolerance") ? pn->Get("tolerance")->As<double>() : 1e-12;
		minTol_     = pn->Has("minimalTolerance") ? pn->Get("minimalTolerance")->As<double>() : 1e-11;
		//logging_    = pn->Has("logging") ? pn->Get("logging")->As<bool>() : false;
		//resetXZero_ = pn->Has("zeroInitialValue") ? pn->Get("zeroInitialValue")->As<bool>() : false;
		
		PtrParamNode hdr = infoNode_->Get(ParamNode::HEADER);
		hdr->Get("maxIter")->SetValue(maxIter_);
		hdr->Get("tolerance")->SetValue(tolerance_);
		hdr->Get("minimalTolerance")->SetValue(minTol_);
		
		// Solve() also sets solver and precond as attributes to xml_ but w/o arguments and also when the elements here are not given
		if(xml_->Has("solver"))
			hdr->Get("solver")->SetValue(xml_->Get("solver"), false);
		if(xml_->Has("precond"))
			hdr->Get("precond")->SetValue(xml_->Get("precond"),false);
		
	}

	PETSCSolver::~PETSCSolver(){

		VecDestroy(&x_);
		VecDestroy(&b_);
		
		MatDestroy(&A_);
		MatDestroy(&A0_);
		
		KSPDestroy(&solver_);
	}

	void PETSCSolver::Setup(BaseMatrix &sysmat){


		int rank;
		int size;
		//find which is my rank
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		MPI_Comm_size(MPI_COMM_WORLD,&size);
		PetscErrorCode ierr=0;
		
		//create petsc matrix form the sysmatrix
		const StdMatrix& stdmat = dynamic_cast<const StdMatrix&>(sysmat);
		
		BaseMatrix::EntryType etype = stdmat.GetEntryType(); //currently only real values can be solved using Petsc should implement for complex values
		BaseMatrix::StorageType stype = stdmat.GetStorageType();

		const CRS_Matrix<Double>& crs = dynamic_cast<const CRS_Matrix<Double>&>(stdmat);
		if(crs.GetNumCols() != crs.GetNumRows()){
      EXCEPTION("IS solver only tested for quadratic matrices");
    }
		//gather info
    UInt dim = crs.GetNumRows();

		PetscInt * rowPtr = (PetscInt *)crs.GetRowPointer();
    PetscInt * colPtr = (PetscInt *)crs.GetColPointer();
		PetscScalar * dataPtr = const_cast<PetscScalar *>(crs.GetDataPointer());
		

		
		int * nnzr = new int[dim];
		
		for (int ii=0;ii<dim;ii++){
			nnzr[ii]=crs.GetRowSize(ii)+100;
		}
		SendWorkerCommand(INIT_MAT_STRUCT);
	
		for (rank=1;rank<size;++rank){	
			MPI_Send(nnzr, sizeof(int)*dim,MPI_INT,rank,SEND_NNZR,PETSC_COMM_WORLD);
		}

		
		if(firstSetup_) {
			
			ierr=MatCreate(PETSC_COMM_WORLD,&A_); 
			ierr=MatSetSizes(A_,PETSC_DECIDE,PETSC_DECIDE,PetscInt(dim),PetscInt(dim));
			ierr=MatSetFromOptions(A_);
			ierr=MatSetUp(A_);
			ierr=MatSeqAIJSetPreallocation(A_,0,(PetscInt *)nnzr);
			ierr=MatMPIAIJSetPreallocation(A_,0,(PetscInt *)nnzr,0,(PetscInt *)nnzr);
			//setting values only using the proc 0 while others wait and then asemble is called in all procs so the matix can be distributed
		}
		else {
			ierr=MatDestroy(&A0_);
		}

		//set matrix values 
		//TODO add a method in the class for running different code root and other procs 
		
		for (int i=0;i<dim;i++){
			for (int j=rowPtr[i];j<rowPtr[i+1];j++){
				ierr=MatSetValue(A_,i,colPtr[j],dataPtr[j],INSERT_VALUES);
			}
		}
			
			// ierr=MatCreateSeqAIJWithArrays(PETSC_COMM_WORLD,PetscInt(dim),PetscInt(dim),rowPtr,colPtr,dataPtr,&A_);
			
		SendWorkerCommand(ASSEMBLE_MAT);	

		MatAssemblyBegin(A_,MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(A_,MAT_FINAL_ASSEMBLY);
		
		
				
	
			
		
		SendWorkerCommand(SETUP_MATRIX);
		//Using this matrix A0_ for solving the linear system so original matrix is untouched maybe Petsc doesnt require this.
		ierr=MatDuplicate(A_,MAT_COPY_VALUES,&A0_);
		
		//Create Vector and solver objects and set parameters 
		if(firstSetup_) {
			//create rhs and sol vector and allocate size in setup step
			//Allocate size for rhs vec assuming the size is equal to dim of sysmat
			ierr=VecCreate(PETSC_COMM_WORLD,&b_);
			ierr=VecSetSizes(b_,PETSC_DECIDE,PetscInt(dim));
			ierr=VecSetFromOptions(b_);
			
			
			//Allocate sol vec 
			ierr=VecDuplicate(b_,&x_);
			
			//setup linear solver context with preconditioner for petsc ksp methods //now hard coded 
			ierr=KSPCreate(PETSC_COMM_WORLD,&(solver_));
			
			// Set up the solver
			ierr = KSPSetType(solver_,KSPCG); // KSPCG - CG SOLVER 

			ierr = KSPSetTolerances(solver_,minTol_,tolerance_,PETSC_DEFAULT,maxIter_);

			ierr = KSPSetInitialGuessNonzero(solver_,PETSC_TRUE);
		
			ierr = KSPSetOperators(solver_,A0_,A0_);

			// The preconditinoer
			KSPGetPC(solver_,&precond_);
			// Make jacobi the default solver
			PCSetType(precond_,PCJACOBI);

			// Set solver from options
			KSPSetFromOptions(solver_);

			// Get the prec again - check if it has changed
			KSPGetPC(solver_,&precond_);

		}
		

	}


	void PETSCSolver::Solve( const BaseMatrix &sysmat, const BaseVector &rhs, BaseVector &sol){
	
		PetscErrorCode ierr=0;
			
		//Suggested to use vecsetvalues instead of this but for now ignoring since there are no pointer to vec data array found	
	
		SendWorkerCommand(ASSEMBLE_VEC_RHS);
		ierr=	VecAssemblyBegin(b_);
		ierr=	VecAssemblyEnd(b_);         

		
		if(sysmat.GetEntryType() == BaseMatrix::DOUBLE) {		
			for(PetscInt i=0, n=(PetscInt)rhs.GetSize(); i<n; i++){
				Double myEnt =0;
				rhs.GetEntry((UInt)i,myEnt);
				ierr=VecSetValue(b_,i,PetscScalar(myEnt),INSERT_VALUES);
			
				
			}
			
		} 
		SendWorkerCommand(ASSEMBLE_VEC_RHS);
		ierr=	VecAssemblyBegin(b_);
		ierr=	VecAssemblyEnd(b_);
		


	
		SendWorkerCommand(SOLVE);
		ierr = KSPSolve(solver_,b_,x_);
		

		
		// PetscInt niter;
		// PetscScalar rnorm;
		// KSPGetIterationNumber(solver_,&niter);
		// KSPGetResidualNorm(solver_,&rnorm); 	
		
		// PetscPrintf(PETSC_COMM_WORLD,"PETSC_SOLVER:  iter: %i, rerr.: %e, time: %f\n",niter,rnorm);
		
		if(ierr){
		EXCEPTION("Solver returned ierror code: " << ierr << " ...aborting")
		}
	
		SendWorkerCommand(GET_SOL);
		
		Vec x_global;
		VecScatter ctx;
		
		//Collect the solution vector from different procs 
		VecScatterCreateToZero(x_,&ctx,&x_global);
		VecScatterBegin(ctx,x_,x_global,INSERT_VALUES,SCATTER_FORWARD);
		VecScatterEnd(ctx,x_,x_global,INSERT_VALUES,SCATTER_FORWARD);
		VecScatterDestroy(&ctx);
	
		//gets the solution vector in bufx array from petsc array 
		PetscScalar *bufx;
		VecGetArray(x_global,&bufx);
		
		if(sysmat.GetEntryType() == BaseMatrix::DOUBLE) {
			for(UInt i=0; i<sol.GetSize();i++){
				
				sol.SetEntry(i,bufx[PetscInt(i)]);
			}
		}
		ierr=VecRestoreArray(x_global,&bufx);
	}	

	
	void PETSCSolver::SendWorkerCommand(int TAG){
		int rank;
		int size;
		//find which is my rank
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		MPI_Comm_size(MPI_COMM_WORLD,&size);
		for (rank = 1; rank < size; ++rank) {
			MPI_Send(0, 0, MPI_INT, rank, TAG, MPI_COMM_WORLD);
		}		
	}
	//Methods for PETSCWorker class
	PETSCWorker::PETSCWorker(){
		
		//init petsc objects
		A_ =NULL;
		A0_=NULL;
		b_=NULL;
		solver_=NULL;
	}


	PETSCWorker::~PETSCWorker(){
		
		VecDestroy(&x_);
		VecDestroy(&b_);
		
		MatDestroy(&A_);
		MatDestroy(&A0_);
		
		KSPDestroy(&solver_);
	}

 
 
	void PETSCWorker::run(){


		PetscErrorCode ierr;
		bool done=false;
		MPI_Status stat;
		int work;
		int dim;

		//Worker loop which runs till a kill command is received
		
		// MPI communication for starting starting some process in worker 
		// The only case where data is sent using MPI is the nnzr array which is required for PetscMat Preallocation

		while (!done){
			MPI_Recv(&work,1,MPI_INT,0,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
			//Need to change to switch case 
			if (stat.MPI_TAG==DIETAG){
				done=true;
			}
			else if (stat.MPI_TAG==INIT_MAT_STRUCT){
				dim=InitPetscWorker();	
			}
			else if (stat.MPI_TAG==ASSEMBLE_MAT){
	
				ierr=	MatAssemblyBegin(A_,MAT_FINAL_ASSEMBLY);
				ierr=	MatAssemblyEnd(A_,MAT_FINAL_ASSEMBLY);
			}
			else if (stat.MPI_TAG==SETUP_MATRIX){
				SetupPetscWorker(dim);
			}
			else if (stat.MPI_TAG==ASSEMBLE_VEC_RHS){
				
				ierr=	VecAssemblyBegin(b_);
				ierr=	VecAssemblyEnd(b_);
			}
			else if (stat.MPI_TAG==SOLVE){
				ierr = KSPSolve(solver_,b_,x_);
			}
			else if (stat.MPI_TAG==GET_SOL){
				GetSol();
			}
		}
  }		 
 

	int PETSCWorker::InitPetscWorker(){
		PetscErrorCode ierr;
		
		MPI_Status status;
		int size;
				
		
		MPI_Probe(0,SEND_NNZR,PETSC_COMM_WORLD,&status);
		MPI_Get_count(&status,MPI_INT,&size);
		int dim=size/sizeof(int);
		int *nnzr = new int[dim];
		
		MPI_Recv(nnzr,dim,MPI_INT,0,SEND_NNZR,PETSC_COMM_WORLD,&status);

		if(firstSetup_) {
			
			ierr=MatCreate(PETSC_COMM_WORLD,&A_); 
			ierr=MatSetSizes(A_,PETSC_DECIDE,PETSC_DECIDE,PetscInt(dim),PetscInt(dim));
			ierr=MatSetFromOptions(A_);
			ierr=MatSetUp(A_);
			ierr=MatSeqAIJSetPreallocation(A_,0,(PetscInt *)nnzr);
			ierr=MatMPIAIJSetPreallocation(A_,0,(PetscInt *)nnzr,0,(PetscInt *)nnzr);
			//setting values only using the proc 0 while others wait and then asemble is called in all procs so the matix can be distributed
		}
		else {
			ierr=MatDestroy(&A0_);
		}

		return dim;

	}

	void PETSCWorker::SetupPetscWorker(int dim){
		PetscErrorCode ierr;
		//Using this matrix A0_ for solving the linear system so original matrix is untouched maybe Petsc doesnt require this.
		ierr=MatDuplicate(A_,MAT_COPY_VALUES,&A0_);
		
		//Create Vector and solver objects and set parameters 
		if(firstSetup_) {
			//create rhs and sol vector and allocate size in setup step
			//Allocate size for rhs vec assuming the size is equal to dim of sysmat
			ierr=VecCreate(PETSC_COMM_WORLD,&b_);
			ierr=VecSetSizes(b_,PETSC_DECIDE,PetscInt(dim));
			ierr=VecSetFromOptions(b_);
			
			
			//Allocate sol vec 
			ierr=VecDuplicate(b_,&x_);
			
			//setup linear solver context with preconditioner for petsc ksp methods //now hard coded 
			ierr=KSPCreate(PETSC_COMM_WORLD,&(solver_));
			
			// Set up the solver
			ierr = KSPSetType(solver_,KSPCG); // KSPCG - CG SOLVER 

			ierr = KSPSetTolerances(solver_,minTol_,tolerance_,PETSC_DEFAULT,maxIter_);

			ierr = KSPSetInitialGuessNonzero(solver_,PETSC_TRUE);
		
			ierr = KSPSetOperators(solver_,A0_,A0_);

			// The preconditinoer
			KSPGetPC(solver_,&precond_);
			// Make jacobi the default solver
			PCSetType(precond_,PCJACOBI);

			// Set solver from options
			KSPSetFromOptions(solver_);

			// Get the prec again - check if it has changed
			KSPGetPC(solver_,&precond_);

		}

	}

	void PETSCWorker::GetSol(){

		Vec x_global;
		VecScatter ctx;
		VecScatterCreateToZero(x_,&ctx,&x_global);
		VecScatterBegin(ctx,x_,x_global,INSERT_VALUES,SCATTER_FORWARD);
		VecScatterEnd(ctx,x_,x_global,INSERT_VALUES,SCATTER_FORWARD);
		VecScatterDestroy(&ctx);
		

	}
	
	

}	


	