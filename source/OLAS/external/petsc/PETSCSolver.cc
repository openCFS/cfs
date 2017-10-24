#include "PETSCSolver.hh"

#include <string>
#include "MatVec/SparseOLASMatrix.hh"
#include "MatVec/SCRS_Matrix.hh"
#include "MatVec/CRS_Matrix.hh"
#include "DataInOut/ProgramOptions.hh"
#include <sstream>
#include <stdio.h>

#define TAG_RESULT 0
#define TAG_ASK_FOR_JOB 1
#define SET_MAT 2
#define DIETAG 3


#ifdef _OPENMP
  #include <omp.h>
#endif

using std::string;

namespace CoupledField{
	//initialize PETSCSolver class 	
	PETSCSolver::PETSCSolver(PtrParamNode pn, PtrParamNode olasInfo, BaseMatrix::EntryType type){
		

		// int argc = -1;                   /* dummy arg count */
		// char **argv = NULL;              /* dummy arg */
		// PetscErrorCode ierr;		
		// ierr = PetscInitialize(NULL,NULL,PETSC_NULL,PETSC_NULL);  
    A_ =NULL;
    A0_=NULL;
    b_=NULL;
    solver_=NULL;
	
	  
		//what is olasinfo get for petsc ?
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
		
		int rank;
		int size;
		//find which is my rank
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		MPI_Comm_size(MPI_COMM_WORLD,&size);

		
		if(firstSetup_) {
			
			ierr=MatCreate(PETSC_COMM_WORLD,&A_);CHKERRV(ierr); 
			ierr=MatSetSizes(A_,PETSC_DECIDE,PETSC_DECIDE,PetscInt(dim),PetscInt(dim));CHKERRV(ierr);
			ierr=MatSetFromOptions(A_);CHKERRV(ierr);
			ierr=MatSetUp(A_);CHKERRV(ierr);
			ierr=MatSeqAIJSetPreallocation(A_,0,(PetscInt *)nnzr);CHKERRV(ierr);
			ierr=MatMPIAIJSetPreallocation(A_,0,(PetscInt *)nnzr,0,(PetscInt *)nnzr);
			//setting values only using the proc 0 while others wait and then asemble is called in all procs so the matix can be distributed
		}
		else {
			ierr=MatDestroy(&A0_);CHKERRV(ierr);
		}

		//set matrix values 
		//TODO add a method in the class for running different code root and other procs 
		if (rank==0){
			for (int i=0;i<dim;i++){
				for (int j=rowPtr[i];j<rowPtr[i+1];j++){
					ierr=MatSetValue(A_,i,colPtr[j],dataPtr[j],INSERT_VALUES);CHKERRV(ierr);
				}
			}
			
			// ierr=MatCreateSeqAIJWithArrays(PETSC_COMM_WORLD,PetscInt(dim),PetscInt(dim),rowPtr,colPtr,dataPtr,&A_);
			if (size>1){
				for (rank = 1; rank < size; ++rank) {
					MPI_Send(0, 0, MPI_INT, rank, DIETAG, MPI_COMM_WORLD);
					
					
				}	
			}
			MatAssemblyBegin(A_,MAT_FINAL_ASSEMBLY);CHKERRV(ierr);
			MatAssemblyEnd(A_,MAT_FINAL_ASSEMBLY);CHKERRV(ierr);
		}	
		else{
			bool done=false;
			MPI_Status stat;
			int work;
			while (!done){
				MPI_Recv(&work,1,MPI_INT,0,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
				if (stat.MPI_TAG==DIETAG){
					MatAssemblyBegin(A_,MAT_FINAL_ASSEMBLY);CHKERRV(ierr);
					MatAssemblyEnd(A_,MAT_FINAL_ASSEMBLY);CHKERRV(ierr);
					done=true;
				}  
				
				
	
			}
		}

		//Using this matrix A0_ for solving the linear system so original matrix is untouched maybe Petsc doesnt require this.
		ierr=MatDuplicate(A_,MAT_COPY_VALUES,&A0_);CHKERRV(ierr);
		
		//Create Vector and solver objects and set parameters 
		if(firstSetup_) {
			//create rhs and sol vector and allocate size in setup step
			//Allocate size for rhs vec assuming the size is equal to dim of sysmat
			ierr=VecCreate(PETSC_COMM_WORLD,&b_);CHKERRV(ierr);
			ierr=VecSetSizes(b_,PETSC_DECIDE,PetscInt(dim));CHKERRV(ierr);
			ierr=VecSetFromOptions(b_);CHKERRV(ierr);
			
			
			//Allocate sol vec 
			ierr=VecDuplicate(b_,&x_);CHKERRV(ierr);
			
			//setup linear solver context with preconditioner for petsc ksp methods //now hard coded 
			ierr=KSPCreate(PETSC_COMM_WORLD,&(solver_));CHKERRV(ierr);
			
			// Set up the solver
			ierr = KSPSetType(solver_,KSPCG);CHKERRV(ierr); // KSPCG - CG SOLVER 

			ierr = KSPSetTolerances(solver_,minTol_,tolerance_,PETSC_DEFAULT,maxIter_);CHKERRV(ierr);

			ierr = KSPSetInitialGuessNonzero(solver_,PETSC_TRUE);CHKERRV(ierr);
		
			ierr = KSPSetOperators(solver_,A0_,A0_);CHKERRV(ierr);

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
		int rank;
		int size;

		
		MPI_Comm_size(MPI_COMM_WORLD, &size);         
		MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
		ierr=	VecAssemblyBegin(b_);CHKERRV(ierr);
		ierr=	VecAssemblyEnd(b_);CHKERRV(ierr);         
		
		if (rank==0){
			if(sysmat.GetEntryType() == BaseMatrix::DOUBLE) {		
				for(PetscInt i=0, n=(PetscInt)rhs.GetSize(); i<n; i++){
					Double myEnt =0;
					rhs.GetEntry((UInt)i,myEnt);
					ierr=VecSetValue(b_,i,PetscScalar(myEnt),INSERT_VALUES);CHKERRV(ierr);
				}
				if (size>1){
					for (rank = 1; rank < size; ++rank) {
						MPI_Send(0, 0, MPI_INT, rank, DIETAG, MPI_COMM_WORLD);
					}					
				}
				
				ierr=	VecAssemblyBegin(b_);CHKERRV(ierr);
				ierr=	VecAssemblyEnd(b_);CHKERRV(ierr);
			} 
		}

		else{
			bool done=false;
			MPI_Status stat;
			int work;
			while (!done){
				MPI_Recv(&work,1,MPI_INT,0,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
				
				if (stat.MPI_TAG==DIETAG){
					ierr=VecAssemblyBegin(b_);CHKERRV(ierr);
					ierr=VecAssemblyEnd(b_);CHKERRV(ierr);
					done=true;
					
				}  
			}
		}
	

		

		
		ierr = KSPSolve(solver_,b_,x_);CHKERRV(ierr);

		
		PetscInt niter;
		PetscScalar rnorm;
		KSPGetIterationNumber(solver_,&niter);
		KSPGetResidualNorm(solver_,&rnorm); 	
		
		// PetscPrintf(PETSC_COMM_WORLD,"PETSC_SOLVER:  iter: %i, rerr.: %e, time: %f\n",niter,rnorm);
		
		if(ierr){
		EXCEPTION("Solver returned ierror code: " << ierr << " ...aborting")
		}
		//lis_solve(A_, b_, x_, solver_);


		//copy solution
		//Lazy initialization of Petsc array pointer to get the petsc vec data
	
	// PetscScalar *bufx;
	// PetscInt low;
	// PetscInt high;
	
	Vec x_global;
	VecDuplicate(x_,&x_global);
	VecScatter ctx;
	PetscScalar *bufx;
	
	
	VecScatterCreateToAll(x_,&ctx,&x_global);
	VecScatterBegin(ctx,x_,x_global,INSERT_VALUES,SCATTER_FORWARD);
	VecScatterEnd(ctx,x_,x_global,INSERT_VALUES,SCATTER_FORWARD);
	VecScatterDestroy(&ctx);
	
	VecGetArray(x_global,&bufx);
	
	if(sysmat.GetEntryType() == BaseMatrix::DOUBLE) {
		for(UInt i=0; i<sol.GetSize();i++){
			
			sol.SetEntry(i,bufx[PetscInt(i)]);
		}
	}
	ierr=VecRestoreArray(x_,&bufx);CHKERRV(ierr);
 }	
}	


	