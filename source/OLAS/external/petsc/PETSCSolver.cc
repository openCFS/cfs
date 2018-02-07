#include "PETSCSolver.hh"

#include <string>
#include "MatVec/SparseOLASMatrix.hh"
#include "MatVec/SCRS_Matrix.hh"
#include "MatVec/CRS_Matrix.hh"
#include "DataInOut/ProgramOptions.hh"
#include "DataInOut/ParamHandling/XmlReader.hh"
#include "Driver/Assemble.hh"
#include "Forms/LinForms/LinearForm.hh"
#include <sstream>
#include <stdio.h>
#include <mpi.h>
#include "PDE/SinglePDE.hh"
#include "Driver/TimeSchemes/BaseTimeScheme.hh"



//Checks Petsc error code, if non-zero it calls the C++ error handler which throws an exception
#define CHKERRXX(ierr)  do {if (PetscUnlikely(ierr)) {PetscError(PETSC_COMM_SELF,__LINE__,PETSC_FUNCTION_NAME,__FILE__,ierr,PETSC_ERROR_IN_CXX,0);}} while(0)

using std::string;

namespace CoupledField{

  DECLARE_LOG(petsc)
  DEFINE_LOG(petsc, "petscSolver")

  void PETSCSolver::createDMDA(){

    //Get a grid pointer so that useful informations can be extracted
    Grid * grid=domain->GetGrid();

    //Get some param pointers for the domain and region list so that we can iterate through all regions and sum up the number of elements in x,y and z
    PtrParamNode paramNode=domain->GetParamRoot()->Get("domain");
    ParamNodeList region_list;

    region_list = paramNode->GetList("regionList");

    StdVector<RegionIdType> regionIds;

    for(unsigned int i=0; i < region_list.GetSize(); i++){
      // we are compatible with the region attribute and unbounded region elements
      string reg = region_list[i]->Get("region")->Has("name") ? region_list[i]->Get("region")->Get("name")->As<std::string>() : region_list[i]->As<std::string>();
      regionIds.Push_back(grid->GetRegionId(reg));
    }
    //Number of nodes in x y and z directions respectively in case of structured with hex elem it is always elemNum in one direction +1
    int nx=1,ny=1,nz=1;
    for (unsigned int j=0;j<regionIds.GetSize();j++){
      StdVector<unsigned int> dim=grid->GetBoundaries(regionIds[j]);
      nx+=dim[0];
      ny+=dim[1];
      nz+=dim[2];
    }

    UInt dimension=grid->GetDim();
    DMBoundaryType bx = DM_BOUNDARY_NONE;
    DMBoundaryType by = DM_BOUNDARY_NONE;
    DMBoundaryType bz = DM_BOUNDARY_NONE;

    if(domain->HasPerdiodicBC()){
    // Boundary types: DMDA_BOUNDARY_NONE, DMDA_BOUNDARY_GHOSTED, DMDA_BOUNDARY_PERIODIC
      EXCEPTION("Multigrid Doesn't Work for periodic domain yet");
    }

    // Stencil type - box since this is closest to FEM
    DMDAStencilType  stype = DMDA_STENCIL_BOX;

    // number of nodal dofs
    PetscInt numnodaldof = 3;
//    std::cout<<"Assuming the nodal Degrees of freedom is 3 ,if something else not implemented yet for Multigrid"<<std::endl;

    // Stencil width: each node connects to a box around it - linear elements
    PetscInt stencilwidth = 1;


    if (dimension==3){
    ierr = DMDACreate3d(PETSC_COMM_WORLD,bx,by,bz,stype,nx,ny,nz,PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,
    numnodaldof,stencilwidth,0,0,0,&(da_nodes));
    }
    else{
    EXCEPTION("Not yet implemented for 2d or 1d");

    }
    //Setup the DM which is used for Grid Management
    DMSetUp(da_nodes);
    DMCreateMatrix(da_nodes,&A_);
    DMDASetElementType(da_nodes, DMDA_ELEMENT_Q1);
    ierr = DMCreateGlobalVector(da_nodes,&(x_));
    VecDuplicate(x_,&(bcopy_));
    PetscInt sizeGloabalVec=0;
    VecGetSize(x_,&sizeGloabalVec);

    VecCreateSeq(PETSC_COMM_WORLD,sizeGloabalVec/PetscInt(3),&N_);

    VecSet(x_,0.0);
    VecSet(N_,1.0);


  }



  void PETSCSolver::AssembleStiffnessMatrix(){
  //This method gets the local element matrix and assembels into global matrix

    //Get the assemble context from the domain
    Assemble* assemble=domain->GetBasePDE()->GetAssemble();

    StdVector<int> eqnVec1;

    //One based in CFS , but we will use 0 based in petsc
    StdVector<int> globalNodeNum;


    //Matrix to store the element stiffness matrix write a method to convert the Matrix to Petsc Matrix
    Matrix<double>elemMatrix;


    Vector<double>elemVector;
    PetscScalar KE[24*24]; // Element stiffness matrix

    // Loop over all bilinear forms and calculate the element stiffness matrix for all regions, in this case since its a structured grid one
    // element stiffness matrix should be sufficient.

    // Using the same implemetation as in LocalElementCache and Assemble.cc , need to verify if its right !!
    for(std::set<BiLinFormContext*>::iterator it = assemble->GetBiLinForms().begin();
        it != assemble->GetBiLinForms().end(); it++ ){
      BiLinFormContext& context = **it;
      BiLinearForm*     form    = context.GetIntegrator();
      RegionIdType reg = context.GetFirstEntities()->GetRegion();
      StdVector<Elem*> elems;
      ElemList elemList(domain->GetGrid());
      // region elements
      domain->GetGrid()->GetElems(elems, reg);
      EntityIterator ei = elemList.GetIterator();
      elemList.SetElement(elems[0]); // region's first element
      assert(ei.GetElem()->elemNum == elems[0]->elemNum);

      //Extract the element node connectivity information which is useful for assembly of Stiffness matrix
      // Everything here is zero based and the globalNodeNum array consist of all the nodes of each element starting
      // from zero to N. It is also assumed that the number of nodes per element is 8
      EntityIterator  firstEntIt = context.GetFirstEntities()->GetIterator();
      for ( firstEntIt.Begin(); !firstEntIt.IsEnd(); firstEntIt++ ) {
        StdVector<unsigned int> nodesInElem=firstEntIt.GetElem()->connect;
        for (unsigned int i =0;i<nodesInElem.GetSize();i++){
          globalNodeNum.Push_back(nodesInElem[i]);
        }

//        std::cout<<"NodeNum for Elements " <<elemNumber<<"   "<<nodesInElem.ToString()<<std::endl;
//        std::cout<<"equVec2ForElements"<<eqnVec2.ToString()<<std::endl;

      }

//      std::cout<<"NodeNum for elements indexed from 1 to n   "<<globalNodeNum.ToString()<<std::endl;

      form->CalcElementMatrix(elemMatrix, ei, ei); // use the region part
//      std::cout<<elemMatrix.ToString()<<std::endl;
      }


    //copy the element matrix to the PetscVector
    for( UInt i=0; i <elemMatrix.GetNumRows(); i++){
         for( UInt j=0; j < elemMatrix.GetNumCols(); j++){
           elemMatrix.GetEntry(i,j,KE[i*(elemMatrix.GetNumCols()) + j]);
//           std::cout<<KE[i*(elemMatrix.GetNumCols()) + j]<<std::endl;
         }
    }

    //Get the number of nodes and elements for the local processor
    PetscInt nel, nen;



    //Get the connectivity informantion: size: (num_elem * num_nodes_per_elem) / proc: content: proc global node number, 0 based
    const PetscInt* necon;

    DMDAGetElements(da_nodes,&nel,&nen,&necon);

    // Zero the matrix
    MatZeroEntries(A_);



    // Edof array
    PetscInt edof[24];

    //Set the Dirchlet Vector With the homogeneous boundary conditions
    GetDirchletVector();
    PetscScalar * dirArray;
    VecGetArray(N_,&dirArray);

    // Loop over elements. nel is number of elements within processor
    for (PetscInt elem=0;elem<nel;elem++){
//      std::cout<<"NodeNum for Elements    "<<elem+1;
      // loop over element nodes
      for (PetscInt elem_node=0;elem_node<nen;elem_node++){

//        std::cout<<" "<<necon[elem*nen+elem_node]+1;
        // Get local dofs
        for (PetscInt dof=0;dof <3;dof++){
          // element local idx 0 ... 23
          int idx = elem_node*3+dof;
          // global idx
          int gidx = globalNodeNum[elem*nen+elem_node]-1;
          int eqnNr = 3*(globalNodeNum[elem*nen+elem_node]-1)+dof; //CFS node number is one  based indexing ,here its zero based
          if (dirArray[gidx] == -1.0){
            homDirEquNr_.Push_back(eqnNr);
            dirchletValue_.Push_back(0);
          }
          edof[idx]=eqnNr;


//          edof[elem_node*3+dof] = 3*necon[elem*nen+elem_node]+dof;
//          std::cout<<eqnNr<<std::endl;
        }//dof Loop
      }//Element Node Loop
//      std::cout<<" "<<std::endl;
      //Set the element Stiffness matrix value to the Global Stiffness matrix
      ierr = MatSetValuesLocal(A_,24,edof,24,edof,KE,ADD_VALUES);
    }//Local Element Loopcoe

//    std::cout<<homDirEquNr_.ToString()<<std::endl;
    MatAssemblyBegin(A_, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A_, MAT_FINAL_ASSEMBLY);




//    PetscInt row=0,col=0;
//    MatGetSize(A_,&row,&col);
//    std::cout<<"Matrix Size"<<row<<"   "<<col<<std::endl;
    DMDARestoreElements(da_nodes,&nel,&nen,&necon);
  }


  void PETSCSolver::SetLinRhs(){

    //Get the RHS


    Vector<Double> elemVec;
    Assemble * assemble=domain->GetBasePDE()->GetAssemble();

    StdVector<LinearFormContext *> linForms_=assemble->GetLinForms();
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
//           StdVector<unsigned int> nodeVec=entIt.GetElem()->connect;
//           for (unsigned int elem_node=0;elem_node<nodeVec.GetSize();elem_node++  ){
             for (int dof=0;dof<3;dof++){
               int eqnNr=3* (entIt.GetNode()-1)+dof;
               eqnVec.Push_back(eqnNr);
             }

//           }
           assert(!elemVec.ContainsNaN() && !elemVec.ContainsInf());
//           std::cout<<elemVec.ToString()<<std::endl;
//           std::cout<<eqnVec.ToString()<<std::endl;
           VecSetValues(bcopy_,eqnVec.GetSize(),eqnVec.GetPointer(),elemVec.GetPointer(),ADD_VALUES);
       }

     }
    VecAssemblyBegin(bcopy_);
    VecAssemblyEnd(bcopy_);
  }


  void  PETSCSolver::GetDirchletVector(){
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
    VecSetValues(N_,DirchletNodeVec.GetSize(),DirchletNodeVec.GetPointer(),HomDirchletValue.GetPointer(),INSERT_VALUES);
  }







	//initialize PETSCSolver class 	
	PETSCSolver::PETSCSolver(PtrParamNode pn, PtrParamNode olasInfo, BaseMatrix::EntryType type){

		//init petsc objects
		A_ =NULL;
		b_=NULL;
		solver_=NULL;
		da_nodes=NULL;
		
		// petsc is mandatory
		infoNode_ =  olasInfo->Get("petsc");
		
		xml_ = pn;
		firstSetup_ = true;

		maxIter_    = pn->Has("maxIter") ? pn->Get("maxIter")->As<int>() : 10000;

		tolerance_  = pn->Has("tolerance") ? pn->Get("tolerance")->As<double>() : 1e-12;
		minTol_     = pn->Has("minimalTolerance") ? pn->Get("minimalTolerance")->As<double>() : 1e-11;
		logging_    = pn->Has("logging") ? pn->Get("logging")->As<bool>() : false;
	
		PtrParamNode hdr = infoNode_->Get(ParamNode::HEADER);
		hdr->Get("maxIter")->SetValue(maxIter_);
		hdr->Get("tolerance")->SetValue(double(tolerance_));
		hdr->Get("minimalTolerance")->SetValue(double(minTol_));
		
		// Solve() also sets solver and precond as attributes to xml_ but w/o arguments and also when the elements here are not given
		if(xml_->Has("solver"))
			hdr->Get("solver")->SetValue(xml_->Get("solver"), false);
		if(xml_->Has("precond"))
			hdr->Get("precond")->SetValue(xml_->Get("precond"),false);
		
		MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
		MPI_Comm_size(MPI_COMM_WORLD,&size_);	

		solverstring_=CreateSolverString(xml_);
		precondstring_=CreatePrecondString(xml_);

		if (precondstring_=="mg"){MG_FLAG=true;}
		createDMDA();
		AssembleStiffnessMatrix();
		SetLinRhs();
		//		GetDirchletVector();
	}

	//Destructor for solver
	PETSCSolver::~PETSCSolver(){

		VecDestroy(&x_);
		VecDestroy(&b_);		
		MatDestroy(&A_);
		KSPDestroy(&solver_);
		
	}

	

	void PETSCSolver::Setup(BaseMatrix &sysmat){


	  MatZeroRowsColumns(A_,homDirEquNr_.GetSize(),homDirEquNr_.GetPointer(),PetscScalar(0.0),NULL,NULL);

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
	  UInt dim = crs.GetNumRows();

	  PetscInt* rowPtr = (PetscInt*) crs.GetRowPointer();
	  PetscInt* colPtr = (PetscInt*) crs.GetColPointer();
	  PetscScalar * dataPtr = const_cast<PetscScalar *>(crs.GetDataPointer());

	  //To preallocate without grid info we require info about number of non-zero per row
	  UInt * nnzr = new UInt[dim];
	  for (UInt ii=0;ii<dim;ii++){
		//100 is just a big enough number which ensures there is no lack of memory allocation for petsc matix in parallel
		  nnzr[ii]=crs.GetRowSize(ii)+100;
	  }



	  if (firstSetup_){
		  SendWorkerCommand(INIT_MAT_STRUCT);
		  //only time where data is sent using mpi(can be implemented better)
		  for (int rank=1;rank<size_;rank++){
        ierr=MPI_Send(&nnzr[0], sizeof(int)*dim,MPI_INT,rank,DATA,PETSC_COMM_WORLD);CHKERRXX(ierr);
        ierr=MPI_Send(&symmetric,sizeof(bool),MPIU_BOOL,rank,ISSYMMETRIC,PETSC_COMM_WORLD);
		  }


		  if (MG_FLAG){
		      createDMDA();
		  }
		  else {


			//Create PETSC matrix structure
        if (symmetric){
//          ierr=MatCreateSBAIJ(PETSC_COMM_WORLD,PetscInt(1),PETSC_DECIDE,PETSC_DECIDE,
//              PetscInt(dim),PetscInt(dim),0,(PetscInt *)nnzr,0,(PetscInt *)nnzr,&A_);
//          ierr=MatSetUp(A_);
//          ierr=MatSetOption (A_, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);CHKERRXX(ierr);
        }
        else {
          //USE THESE WHEN MATRIX IS OF TYPE SPARSE NON SYM
            ierr=MatCreate(PETSC_COMM_WORLD,&A_); CHKERRXX(ierr);
            ierr=MatSetSizes(A_,PETSC_DECIDE,PETSC_DECIDE,PetscInt(dim),PetscInt(dim));CHKERRXX(ierr);
            ierr=MatSetFromOptions(A_);CHKERRXX(ierr);
            ierr=MatSetUp(A_);CHKERRXX(ierr);
            ierr=MatSeqAIJSetPreallocation(A_,0,(PetscInt *)nnzr);CHKERRXX(ierr);
            ierr=MatMPIAIJSetPreallocation(A_,0,(PetscInt *)nnzr,0,(PetscInt *)nnzr);CHKERRXX(ierr);
            ierr=MatSetOption (A_, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);CHKERRXX(ierr);
        }
		  }

	  }


	//setting values only using the proc 0 while others wait and then asemble is called in all procs so t
	// the matix can be distributed
	//Set each matrix value one by one (highly inefficient there are better ways to do this which requires some tweakng to
	//the data structure we receive from the matrix class	)
//	for (int i=0;i<int(dim);i++){
//	  for (int j=rowPtr[i];j<rowPtr[i+1];j++){
//	    ierr=MatSetValue(A_,i,colPtr[j],dataPtr[j],INSERT_VALUES);
//	  }
//	}

	//Distribute the assembeled matrix accross all process
//	SendWorkerCommand(ASSEMBLE_MAT);
//	ierr=MatAssemblyBegin(A_,MAT_FINAL_ASSEMBLY);CHKERRXX(ierr);
//	ierr=MatAssemblyEnd(A_,MAT_FINAL_ASSEMBLY);CHKERRXX(ierr);

	if (firstSetup_){


		SendWorkerCommand(SETUP_MATRIX);

		//Create Vector and solver objects and set parameters

		//create rhs and sol vector and allocate size in setup step
		//Allocate size for rhs vec assuming the size is equal to dim of sysmat

		//creates rhs and lhs vector from Matrix with appropriate size
//		ierr=MatCreateVecs(A_,&b_,&x_);CHKERRXX(ierr);

		//setup linear solver context with preconditioner for petsc ksp methods //now hard coded
		ierr=KSPCreate(PETSC_COMM_WORLD,&(solver_));CHKERRXX(ierr);

		ierr = KSPSetTolerances(solver_,tolerance_,minTol_,PETSC_DEFAULT,maxIter_);CHKERRXX(ierr);

		ierr = KSPSetInitialGuessNonzero(solver_,PETSC_TRUE);CHKERRXX(ierr);

		ierr = KSPSetType(solver_,solverstring_.c_str());CHKERRXX(ierr);

		// The preconditinoer
		ierr=KSPGetPC(solver_,&precond_);CHKERRXX(ierr);



		// Set the preconditioner
		PCSetType(precond_,precondstring_.c_str());CHKERRXX(ierr);

//		if(MG_FLAG){

//			DM  *da_list,*daclist;
//			Mat R;
//
//
//			//create DMs for all levels of the multigrid so that Petsc can construct coarse grids
//			PetscMalloc(sizeof(DM)*nlvls,&da_list);
//			PetscMalloc(sizeof(DM)*nlvls,&daclist);
//			for (PetscInt k=0; k<nlvls; k++){
//			  da_list[k] = NULL;
//			  daclist[k] = NULL;
//			}
//
//			//Set the finest nodal mesh as the first item in the da list
//			daclist[0]=da_nodes;
//			// petsc constructs DM structures for all other grid level
//			DMCoarsenHierarchy(da_nodes,nlvls-1,&daclist[1]);
//			for (PetscInt k=0; k<nlvls; k++) {
//			  // NOTE: finest grid is nlevels - 1: PCMG MUST USE THIS ORDER ???
//			  da_list[k] = daclist[nlvls-1-k];
//			  DMDASetUniformCoordinates(da_list[k],0,1,0,1,0,1);CHKERRXX(ierr);
//
//			}
//
//			PCMGSetLevels(precond_,nlvls,NULL);CHKERRXX(ierr);
//			PCMGSetType(precond_,PC_MG_MULTIPLICATIVE); CHKERRXX(ierr);// Default
//			ierr = PCMGSetCycleType(precond_,PC_MG_CYCLE_V);CHKERRXX(ierr);
//			PCMGSetGalerkin(precond_,PC_MG_GALERKIN_BOTH);CHKERRXX(ierr);
//
//			for (PetscInt k=1; k<nlvls; k++) {
//			  DMCreateInterpolation(da_list[k-1],da_list[k],&R,NULL);CHKERRXX(ierr);
//			  PCMGSetInterpolation(precond_,k,R);CHKERRXX(ierr);
//			  MatDestroy(&R);CHKERRXX(ierr);
//			}
//			// Now all DM data structure other than the finest level can be destroyed
//			for (PetscInt k=1; k<nlvls; k++) { // DO NOT DESTROY LEVEL 0
//			  DMDestroy(&daclist[k]);CHKERRXX(ierr);
//			}
//
//			PetscFree(da_list);CHKERRXX(ierr);
//			PetscFree(daclist);CHKERRXX(ierr);
//
//			KSP cksp;
//			PCMGGetCoarseSolve(precond_,&cksp);
//		  // The solver
//			ierr = KSPSetType(cksp,KSPGMRES);
//
//			ierr=KSPSetTolerances(cksp,coarse_rtol,coarse_atol,coarse_dtol,coarse_maxits);
//
//			PC cpc;
//			KSPGetPC(cksp,&cpc);
//			PCSetType(cpc,PCSOR);
//
//			for (PetscInt k=1;k<nlvls;k++){
//			  KSP dksp;
//			  PCMGGetSmoother(precond_,k,&dksp);
//			  PC dpc;
//			  KSPGetPC(dksp,&dpc);
//			  ierr = KSPSetType(dksp,KSPGMRES); // KSPCG, KSPGMRES, KSPCHEBYSHEV (VERY GOOD FOR SPD)
//
//			  ierr = KSPSetTolerances(dksp,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,smooth_sweeps); // NOTE in the above maxitr=restart;
//			  PCSetType(dpc,PCSOR);// PCJACOBI, PCSOR for KSPCHEBYSHEV very good
//			}
//		}

		PetscOptionsInsert(NULL,NULL,NULL,NULL);CHKERRXX(ierr);

		// Set solver from options
		KSPSetFromOptions(solver_);CHKERRXX(ierr);

		// Get the prec again - check if it has changed
		KSPGetPC(solver_,&precond_);CHKERRXX(ierr);

		ierr = KSPSetOperators(solver_,A_,A_);CHKERRXX(ierr);



	}
	firstSetup_=false;
			
	}


	void PETSCSolver::Solve( const BaseMatrix &sysmat, const BaseVector &rhs, BaseVector &sol){
	

			
		//Suggested to use vecsetvalues instead of this but for now ignoring since there are no pointer to vec data array found	
		
		//Assemble the vector before setting values using a single process 
		//doesnt have any effect on performance in single, may be crucial when using more than one node so keeping it

//		SendWorkerCommand(ASSEMBLE_VEC_RHS);
//		ierr=	VecAssemblyBegin(b_);CHKERRXX(ierr);
//		ierr=	VecAssemblyEnd(b_);CHKERRXX(ierr);
//		if(sysmat.GetEntryType() == BaseMatrix::DOUBLE) {
//			for(PetscInt i=0, n=(PetscInt)rhs.GetSize(); i<n; i++){
//				Double myEnt =0;
//				rhs.GetEntry((UInt)i,myEnt);
//				ierr=VecSetValue(b_,i,PetscScalar(myEnt),INSERT_VALUES);
//			}
//		}
//		//Distribute the vector after values are set
//		SendWorkerCommand(ASSEMBLE_VEC_RHS);
//		ierr=	VecAssemblyBegin(b_);CHKERRXX(ierr);
//		ierr=	VecAssemblyEnd(b_);CHKERRXX(ierr);
//		PetscBool  vecequal;
//		PetscInt bsize=0,bcopysize=0;
//		VecEqual(b_,bcopy_,&vecequal);
//		VecGetSize(b_,&bsize);
//		VecGetSize(bcopy_,&bcopysize);
//		std::cout<<bsize<<std::endl;
//		std::cout<<bcopysize<<std::endl;

//	  Vec f;
//	  //Time for solve
//
//	  VecDuplicate(bcopy_,&f);
//	  VecSet(f,1.0);
//	  VecSetValues(f,homDirEquNr_.GetSize(),homDirEquNr_.GetPointer(),(PetscScalar*)(dirchletValue_.GetPointer()),INSERT_VALUES);
//
//	  MatDiagonalScale(A_,f,f);
//	  Vec NI;
//    VecDuplicate(f,&NI);
//    VecSet(NI,1.0);
//    VecAXPY(NI,-1.0,f);
//    MatDiagonalSet(A_,NI,ADD_VALUES);
//
//  // Zero out possible loads in the RHS that coincide
//    // with Dirichlet conditions
//    VecPointwiseMult(bcopy_,bcopy_,f);
//
//    VecDestroy(&NI);


		double t1,t2;
		t1 = MPI_Wtime();
		//Solve the linear system
		SendWorkerCommand(SOLVE);
		ierr = KSPSolve(solver_,bcopy_,x_);CHKERRXX(ierr);
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
		PetscInt Size;
		VecGetSize(x_global,&Size);
		sol.Resize(Size);

    for(int i=0; i<Size;i++){
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
		A_ =NULL;
		b_=NULL;
		solver_=NULL;
		x_=NULL;	

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
			if (precondstring_=="mg"){MG_FLAG=true;}
		}
	
	}


	PETSCWorker::~PETSCWorker(){
		
		VecDestroy(&x_);
		VecDestroy(&b_);
		
		MatDestroy(&A_);
	
		KSPDestroy(&solver_);

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
		//Receive the nnzr array 
		MPI_Recv(nnzr,dim,MPI_INT,0,DATA,PETSC_COMM_WORLD,&status);
		MPI_Recv(&symmetric,1,MPIU_BOOL,0,ISSYMMETRIC,PETSC_COMM_WORLD,&status);
		
	    if (MG_FLAG){
		  // Boundary types: DMDA_BOUNDARY_NONE, DMDA_BOUNDARY_GHOSTED, DMDA_BOUNDARY_PERIODIC
		  DMBoundaryType bx = DM_BOUNDARY_NONE;
		  DMBoundaryType by = DM_BOUNDARY_NONE;
		  DMBoundaryType bz = DM_BOUNDARY_NONE;

		  // Stencil type - box since this is closest to FEM (i.e. STAR is FV/FD)
		  DMDAStencilType  stype = DMDA_STENCIL_BOX;


		  //Set the global dimension from the mesh info
		  PetscInt nx = 11;
		  PetscInt ny = 11;
		  PetscInt nz =	11;

		  // number of nodal dofs
		  PetscInt numnodaldof = 3;

		  // Stencil width: each node connects to a box around it - linear elements
		  PetscInt stencilwidth = 1;




		  ierr = DMDACreate3d(PETSC_COMM_WORLD,bx,by,bz,stype,nx,ny,nz,PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,
			  numnodaldof,stencilwidth,0,0,0,&(da_nodes));

		  DMSetUp(da_nodes);
		  DMCreateMatrix(da_nodes,&A_);
		  ierr=MatSetOption (A_, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);CHKERRXX(ierr);

	  }
	  else {
		  //Create PETSC matrix structure
		if (symmetric){
			ierr=MatCreateSBAIJ(PETSC_COMM_WORLD,PetscInt(1),PETSC_DECIDE,PETSC_DECIDE,
					PetscInt(dim),PetscInt(dim),0,(PetscInt *)nnzr,0,(PetscInt *)nnzr,&A_);
			ierr=MatSetUp(A_);
			ierr=MatSetOption (A_, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);CHKERRXX(ierr);
		}
		else{
			//USE THESE WHEN MATRIX IS OF TYPE SPARSE NON SYM
				ierr=MatCreate(PETSC_COMM_WORLD,&A_); CHKERRXX(ierr);
				ierr=MatSetSizes(A_,PETSC_DECIDE,PETSC_DECIDE,PetscInt(dim),PetscInt(dim));CHKERRXX(ierr);
				ierr=MatSetFromOptions(A_);CHKERRXX(ierr);
				ierr=MatSetUp(A_);CHKERRXX(ierr);
				ierr=MatSeqAIJSetPreallocation(A_,0,(PetscInt *)nnzr);CHKERRXX(ierr);
				ierr=MatMPIAIJSetPreallocation(A_,0,(PetscInt *)nnzr,0,(PetscInt *)nnzr);CHKERRXX(ierr);
				ierr=MatSetOption (A_, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);CHKERRXX(ierr);
		}
	  }

	}





	void PETSCWorker::SetupPetscWorker(){
		
		
	//creates rhs and lhs vector from Matrix with appropriate size
		ierr=MatCreateVecs(A_,&b_,&x_);CHKERRXX(ierr);

		//setup linear solver context with preconditioner for petsc ksp methods //now hard coded
		ierr=KSPCreate(PETSC_COMM_WORLD,&(solver_));CHKERRXX(ierr);

		ierr = KSPSetTolerances(solver_,tolerance_,minTol_,PETSC_DEFAULT,maxIter_);CHKERRXX(ierr);

		ierr = KSPSetInitialGuessNonzero(solver_,PETSC_TRUE);CHKERRXX(ierr);

		ierr = KSPSetType(solver_,solverstring_.c_str());CHKERRXX(ierr);

		// The preconditinoer
		ierr=KSPGetPC(solver_,&precond_);

		// Set the preconditioner
		PCSetType(precond_,precondstring_.c_str());

		PetscOptionsInsert(NULL,NULL,NULL,NULL);

		// Set solver from options
		KSPSetFromOptions(solver_);CHKERRXX(ierr);

		// Get the prec again - check if it has changed
		KSPGetPC(solver_,&precond_);CHKERRXX(ierr);

		ierr = KSPSetOperators(solver_,A_,A_);CHKERRXX(ierr);
	}

	void PETSCWorker::GetSol(){

		Vec x_global;
		VecScatter ctx;
		VecScatterCreateToZero(x_,&ctx,&x_global);
		VecScatterBegin(ctx,x_,x_global,INSERT_VALUES,SCATTER_FORWARD);
		VecScatterEnd(ctx,x_,x_global,INSERT_VALUES,SCATTER_FORWARD);
		VecScatterDestroy(&ctx);	

	}


	std::string PETSCSolver::CreateSolverString(PtrParamNode xml_){
		
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
	
	
	std::string PETSCSolver::CreatePrecondString(PtrParamNode xml_){
		
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


	std::string PETSCWorker::CreateSolverString(PtrParamNode xml_){
		
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
	
	
	std::string PETSCWorker::CreatePrecondString(PtrParamNode xml_){
		
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



}	


