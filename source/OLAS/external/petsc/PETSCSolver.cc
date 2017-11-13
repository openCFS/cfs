#include "PETSCSolver.hh"

#include <string>
#include "MatVec/SparseOLASMatrix.hh"
#include "MatVec/SCRS_Matrix.hh"
#include "MatVec/CRS_Matrix.hh"
#include "DataInOut/ProgramOptions.hh"
#include <sstream>
#include <stdio.h>
#include<mpi.h>

#define DIETAG 0
#define INIT_MAT_STRUCT 1
#define ASSEMBLE_MAT 2
#define ASSEMBLE_VEC_RHS  3
#define SETUP_MATRIX 4
#define SOLVE 5
#define DATA 6
#define GET_SOL 7
#define SOLVER_STRING 8

//Checks Petsc error code, if non-zero it calls the C++ error handler which throws an exception
#define CHKERRXX(ierr)  do {if (PetscUnlikely(ierr)) {PetscError(PETSC_COMM_SELF,__LINE__,PETSC_FUNCTION_NAME,__FILE__,ierr,PETSC_ERROR_IN_CXX,0);}} while(0)

using std::string;

namespace CoupledField{


	std::string PETSCSolver::CreateSolverString(){
		
		string output;
		PtrParamNode sNode = xml_->Get("solver",ParamNode::PASS);
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
	
	
	std::string PETSCSolver::CreatePrecondString(){
		
		string output;
		PtrParamNode sNode = xml_->Get("precond",ParamNode::PASS);
		ParamNodeList sol = sNode->GetChildren();
	
		UInt numChilds = sol.GetSize();
		if(numChilds==2){
			output = sol[1]->GetName();
		}
		else{
			output="jacobi";
		}
		return output;
	}
		
	//initialize PETSCSolver class 	
	PETSCSolver::PETSCSolver(PtrParamNode pn, PtrParamNode olasInfo, BaseMatrix::EntryType type){

		//init petsc objects
		A_ =NULL;
		b_=NULL;
		solver_=NULL;

		
		infoNode_ =  olasInfo->Get("petsc");
		
		xml_ = pn;
		firstSetup_ = true;
	
		
		maxIter_    = pn->Has("maxIter") ? pn->Get("maxIter")->As<int>() : 10000;
		tolerance_  = pn->Has("tolerance") ? pn->Get("tolerance")->As<double>() : 1e-12;
		minTol_     = pn->Has("minimalTolerance") ? pn->Get("minimalTolerance")->As<double>() : 1e-11;
		logging_    = pn->Has("logging") ? pn->Get("logging")->As<bool>() : false;
	
		PtrParamNode hdr = infoNode_->Get(ParamNode::HEADER);
		hdr->Get("maxIter")->SetValue(maxIter_);
		hdr->Get("tolerance")->SetValue(tolerance_);
		hdr->Get("minimalTolerance")->SetValue(minTol_);
		
		// Solve() also sets solver and precond as attributes to xml_ but w/o arguments and also when the elements here are not given
		if(xml_->Has("solver"))
			hdr->Get("solver")->SetValue(xml_->Get("solver"), false);
		if(xml_->Has("precond"))
			hdr->Get("precond")->SetValue(xml_->Get("precond"),false);
		
		MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
		MPI_Comm_size(MPI_COMM_WORLD,&size_);	

		solverstring_=CreateSolverString();
		precondstring_=CreatePrecondString();
	
		
	}

	//Destructor for solver
	PETSCSolver::~PETSCSolver(){

		VecDestroy(&x_);
		VecDestroy(&b_);		
		MatDestroy(&A_);
		KSPDestroy(&solver_);
		
	}

	

	void PETSCSolver::Setup(BaseMatrix &sysmat){

		//create petsc matrix form the sysmatrix
		// const StdMatrix& stdmat = static_cast<const StdMatrix&>(sysmat);
		
		// BaseMatrix::EntryType etype = stdmat.GetEntryType(); //currently only real values can be solved using Petsc should implement for complex values
		// BaseMatrix::StorageType stype = stdmat.GetStorageType();

		

		const CRS_Matrix<Double>& crs = static_cast<const CRS_Matrix<Double>&>(sysmat);
		if(crs.GetNumCols() != crs.GetNumRows()){
      EXCEPTION("PETSC solver only tested for quadratic matrices");
    }
		//gather info
    UInt dim = crs.GetNumRows();

		PetscInt * rowPtr = (PetscInt *)crs.GetRowPointer();
    PetscInt * colPtr = (PetscInt *)crs.GetColPointer();
		PetscScalar * dataPtr = const_cast<PetscScalar *>(crs.GetDataPointer());
		

		
		UInt * nnzr = new UInt[dim];
		for (UInt ii=0;ii<dim;ii++){
			//100 is just a big enough number which ensures there is no lack of memory allocation for petsc matix in parallel 
			nnzr[ii]=crs.GetRowSize(ii)+100;
		}
		//Create Matrix and solver objects only during the first setup 
		if (firstSetup_){

			SendWorkerCommand(INIT_MAT_STRUCT);
			
			//only time where data is sent using mpi(can be implemented better)
			for (int rank=1;rank<size_;rank++){	
				ierr=MPI_Send(&nnzr[0], sizeof(int)*dim,MPI_INT,rank,DATA,PETSC_COMM_WORLD);CHKERRXX(ierr);
				ierr=MPI_Send(&minTol_,sizeof(double),MPI_DOUBLE,rank,DATA,PETSC_COMM_WORLD);CHKERRXX(ierr);
				ierr=MPI_Send(&tolerance_,sizeof(double),MPI_DOUBLE,rank,DATA,PETSC_COMM_WORLD);CHKERRXX(ierr);
				ierr=MPI_Send(&maxIter_,sizeof(int),MPI_INT,rank,DATA,PETSC_COMM_WORLD);CHKERRXX(ierr);
				ierr=MPI_Send(solverstring_.c_str(),solverstring_.size()+1,MPI_CHAR,rank,SOLVER_STRING,PETSC_COMM_WORLD);CHKERRXX(ierr);
				ierr=MPI_Send(precondstring_.c_str(),precondstring_.size()+1,MPI_CHAR,rank,SOLVER_STRING,PETSC_COMM_WORLD);CHKERRXX(ierr);
			}
			
			//Create PETSC matrix structure
			
			ierr=MatCreate(PETSC_COMM_WORLD,&A_); CHKERRXX(ierr);
			ierr=MatSetSizes(A_,PETSC_DECIDE,PETSC_DECIDE,PetscInt(dim),PetscInt(dim));CHKERRXX(ierr);
			ierr=MatSetFromOptions(A_);CHKERRXX(ierr);
			ierr=MatSetUp(A_);CHKERRXX(ierr);
			ierr=MatSeqAIJSetPreallocation(A_,0,(PetscInt *)nnzr);CHKERRXX(ierr);
			ierr=MatMPIAIJSetPreallocation(A_,0,(PetscInt *)nnzr,0,(PetscInt *)nnzr);CHKERRXX(ierr);
			ierr=MatSetOption (A_, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);CHKERRXX(ierr);
		}
			
		//setting values only using the proc 0 while others wait and then asemble is called in all procs so the matix can be distributed
		for (int i=0;i<int(dim);i++){
			for (int j=rowPtr[i];j<rowPtr[i+1];j++){
				ierr=MatSetValue(A_,i,colPtr[j],dataPtr[j],INSERT_VALUES);
			}
		}
			//Deprecated Version of setting petsc matrix for sequential jobs. Not recommended to use it.
			// ierr=MatCreateSeqAIJWithArrays(PETSC_COMM_WORLD,PetscInt(dim),PetscInt(dim),rowPtr,colPtr,dataPtr,&A_);
		
		//Distribute the assembeled matrix accross all process	
		SendWorkerCommand(ASSEMBLE_MAT);	
		ierr=MatAssemblyBegin(A_,MAT_FINAL_ASSEMBLY);CHKERRXX(ierr);
		ierr=MatAssemblyEnd(A_,MAT_FINAL_ASSEMBLY);CHKERRXX(ierr);
		
		
				
	
		if (firstSetup_){

			
			SendWorkerCommand(SETUP_MATRIX);
			
			//Create Vector and solver objects and set parameters 
			
			//create rhs and sol vector and allocate size in setup step
			//Allocate size for rhs vec assuming the size is equal to dim of sysmat
			ierr=VecCreate(PETSC_COMM_WORLD,&b_);CHKERRXX(ierr);
			ierr=VecSetSizes(b_,PETSC_DECIDE,PetscInt(dim));CHKERRXX(ierr);
			ierr=VecSetFromOptions(b_);CHKERRXX(ierr);
			
			
			//Allocate sol vec 
			ierr=VecDuplicate(b_,&x_);CHKERRXX(ierr);
			
			//setup linear solver context with preconditioner for petsc ksp methods //now hard coded 
			ierr=KSPCreate(PETSC_COMM_WORLD,&(solver_));CHKERRXX(ierr);
			
			// Set up the solver
			ierr = KSPSetType(solver_,solverstring_.c_str()); CHKERRXX(ierr);// KSPCG - CG SOLVER 
			
			ierr = KSPSetTolerances(solver_,tolerance_,minTol_,PETSC_DEFAULT,maxIter_);CHKERRXX(ierr);
			
			ierr = KSPSetInitialGuessNonzero(solver_,PETSC_TRUE);CHKERRXX(ierr);
			
			ierr = KSPSetOperators(solver_,A_,A_);CHKERRXX(ierr);
			
		
			
			// The preconditinoer
			KSPGetPC(solver_,&precond_);CHKERRXX(ierr);
			// Make jacobi the default preconditioner
			PCSetType(precond_,precondstring_.c_str());CHKERRXX(ierr);
			
			PetscOptionsInsert(NULL,NULL,NULL,NULL);
		
			// // Set solver from options
			KSPSetFromOptions(solver_);CHKERRXX(ierr);
			
			// Get the prec again - check if it has changed
			KSPGetPC(solver_,&precond_);CHKERRXX(ierr);
			
		}		
		firstSetup_=false;
			
	}


	void PETSCSolver::Solve( const BaseMatrix &sysmat, const BaseVector &rhs, BaseVector &sol){
	
		
			
		//Suggested to use vecsetvalues instead of this but for now ignoring since there are no pointer to vec data array found	
		
		//Assemble the vector before setting values using a single process 
		//doesnt have any effect on performance in single, may be crucial when using more than one node so keeping it
		SendWorkerCommand(ASSEMBLE_VEC_RHS);
		ierr=	VecAssemblyBegin(b_);CHKERRXX(ierr);
		ierr=	VecAssemblyEnd(b_);CHKERRXX(ierr);         
		if(sysmat.GetEntryType() == BaseMatrix::DOUBLE) {		
			for(PetscInt i=0, n=(PetscInt)rhs.GetSize(); i<n; i++){
				Double myEnt =0;
				rhs.GetEntry((UInt)i,myEnt);
				ierr=VecSetValue(b_,i,PetscScalar(myEnt),INSERT_VALUES);
			}
		} 

		//Distribute the vector after values are set
		SendWorkerCommand(ASSEMBLE_VEC_RHS);
		ierr=	VecAssemblyBegin(b_);CHKERRXX(ierr);
		ierr=	VecAssemblyEnd(b_);CHKERRXX(ierr);

		//Time for solve
		double t1,t2;
		t1 = MPI_Wtime();
		//Solve the linear system
		SendWorkerCommand(SOLVE);
		ierr = KSPSolve(solver_,b_,x_);CHKERRXX(ierr);
		t2 = MPI_Wtime();

		// Log the the norm and max iter 
		PtrParamNode curr = infoNode_->Get(ParamNode::PROCESS)->Get("solve", ParamNode::APPEND); 
		PetscInt niter=0;
		PetscScalar rnorm=0.0;
		
		KSPGetIterationNumber(solver_,&niter);CHKERRXX(ierr);
		KSPGetResidualNorm(solver_,&rnorm);CHKERRXX(ierr);
		
		curr->Get("timing")->SetValue(t2-t1);
		curr->Get("residualNorm")->SetValue(rnorm);
		curr->Get("iterations")->SetValue(niter);

		//Create a global vector and scatter context to collect the solution vector to master process
		SendWorkerCommand(GET_SOL);
		Vec x_global;
		VecScatter ctx;
		
		//Collect the solution vector from different procs 
		ierr=VecScatterCreateToZero(x_,&ctx,&x_global);CHKERRXX(ierr);
		ierr=VecScatterBegin(ctx,x_,x_global,INSERT_VALUES,SCATTER_FORWARD);CHKERRXX(ierr);
		ierr=VecScatterEnd(ctx,x_,x_global,INSERT_VALUES,SCATTER_FORWARD);CHKERRXX(ierr);
		ierr=VecScatterDestroy(&ctx);CHKERRXX(ierr);
	
		//gets the solution vector in bufx array from petsc array 
		PetscScalar *bufx;
		ierr=VecGetArray(x_global,&bufx);CHKERRXX(ierr);
		
		if(sysmat.GetEntryType() == BaseMatrix::DOUBLE) {
			for(UInt i=0; i<sol.GetSize();i++){
				sol.SetEntry(i,bufx[PetscInt(i)]);
			}
		}
		ierr=VecRestoreArray(x_global,&bufx);CHKERRXX(ierr);
	}	

		
	void PETSCSolver::SendWorkerCommand(int TAG){
	
		for (int rank = 1; rank < size_; rank++) {
			ierr=MPI_Send(0, 0, MPI_INT, rank, TAG, MPI_COMM_WORLD);CHKERRXX(ierr);
		}		
	}


	//Methods for PETSCWorker class
	PETSCWorker::PETSCWorker(){
		//init petsc objects
		A_ =NULL;
		b_=NULL;
		solver_=NULL;
		x_=NULL;		
	}


	PETSCWorker::~PETSCWorker(){
		
		VecDestroy(&x_);
		VecDestroy(&b_);
		
		MatDestroy(&A_);
	
		
		KSPDestroy(&solver_);

		delete [] solverstring_;
		delete [] precondstring_;
		
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
				case INIT_MAT_STRUCT:
					InitPetscWorker();	
					break;
				case ASSEMBLE_MAT:
					ierr=	MatAssemblyBegin(A_,MAT_FINAL_ASSEMBLY);
					ierr=	MatAssemblyEnd(A_,MAT_FINAL_ASSEMBLY);
					break;
				case SETUP_MATRIX:
					SetupPetscWorker();
					break;
				case ASSEMBLE_VEC_RHS:
					ierr=	VecAssemblyBegin(b_);
					ierr=	VecAssemblyEnd(b_);
					break;
				case SOLVE:
					ierr = KSPSolve(solver_,b_,x_);
					break;
				case GET_SOL:
					GetSol();
					break;
			}
		}
  }		 
 

	void PETSCWorker::InitPetscWorker(){
		
		MPI_Status status;
		int size;		
		
		MPI_Probe(0,DATA,PETSC_COMM_WORLD,&status);
		MPI_Get_count(&status,MPI_INT,&size);
		dim=size/sizeof(int);
		UInt *nnzr = new UInt[dim];
		
		//Receive the nnzr array and all other global parameters for solver
		MPI_Recv(nnzr,dim,MPI_INT,0,DATA,PETSC_COMM_WORLD,&status);
		MPI_Recv(&minTol_,1,MPI_DOUBLE,0,DATA,PETSC_COMM_WORLD,&status);
		MPI_Recv(&tolerance_,1,MPI_DOUBLE,0,DATA,PETSC_COMM_WORLD,&status);
		MPI_Recv(&maxIter_,1,MPI_INT,0,DATA,PETSC_COMM_WORLD,&status);
		
		//Receive the solver type data from master
		int solverStringSize;
		MPI_Probe(0,SOLVER_STRING,PETSC_COMM_WORLD,&status);
		MPI_Get_count(&status,MPI_CHAR,&solverStringSize);
		solverstring_=new char[solverStringSize];
		ierr=MPI_Recv(solverstring_,solverStringSize,MPI_CHAR,0,SOLVER_STRING,PETSC_COMM_WORLD,&status);
	
		

		//Receive the preconditioner type data from  master
		int precondStringSize;
		MPI_Probe(0,SOLVER_STRING,PETSC_COMM_WORLD,&status);
		MPI_Get_count(&status,MPI_CHAR,&precondStringSize);
		precondstring_=new char[precondStringSize];
		ierr=MPI_Recv(precondstring_,precondStringSize,MPI_CHAR,0,SOLVER_STRING,PETSC_COMM_WORLD,&status);
		
		
		
		//Same as master class but all workers also needs to call the petsc commands since most petsc commands are collective
		ierr=MatCreate(PETSC_COMM_WORLD,&A_); 
		ierr=MatSetSizes(A_,PETSC_DECIDE,PETSC_DECIDE,PetscInt(dim),PetscInt(dim));
		ierr=MatSetFromOptions(A_);
		ierr=MatSetUp(A_);
		ierr=MatSeqAIJSetPreallocation(A_,0,(PetscInt *)nnzr);
		ierr=MatMPIAIJSetPreallocation(A_,0,(PetscInt *)nnzr,0,(PetscInt *)nnzr);
		ierr=MatSetOption (A_, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);CHKERRXX(ierr);

	}

	void PETSCWorker::SetupPetscWorker(){
		
		
		//Create Vector and solver objects and set parameters 
		
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
		ierr = KSPSetType(solver_,string(solverstring_).c_str()); // KSPCG - CG SOLVER 

		ierr = KSPSetTolerances(solver_,tolerance_,minTol_,PETSC_DEFAULT,maxIter_);

		ierr = KSPSetInitialGuessNonzero(solver_,PETSC_TRUE);
	
		ierr = KSPSetOperators(solver_,A_,A_);

		// The preconditinoer
		KSPGetPC(solver_,&precond_);CHKERRXX(ierr);
		// Make jacobi the default preconditioner
		PCSetType(precond_,string(precondstring_).c_str());CHKERRXX(ierr);
		
		PetscOptionsInsert(NULL,NULL,NULL,NULL);
		
		// Set solver from options
		KSPSetFromOptions(solver_);

		// Get the prec again - check if it has changed
		KSPGetPC(solver_,&precond_);
	
		
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


	