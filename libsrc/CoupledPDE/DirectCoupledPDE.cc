#include "DirectCoupledPDE.hh"

#include "BasePairCoupling.hh"

// include solveStep drivers
#include "Driver/stdSolveStep.hh"
#include "Driver/solveStepMech.hh"

// include PDE classes
#include "PDE/SinglePDE.hh"

#include "PDE/nodeEQN.hh"


namespace CoupledField{


DirectCoupledPDE::DirectCoupledPDE(Grid *aptgrid,  
				   WriteResults *aOutFile, 
				   TimeFunc *aTimeFunc)
  : StdPDE(aptgrid, aOutFile, aTimeFunc)
{
  ENTER_FCN( "DirectCoupledPDE::DirectCoupledPDE" );

}


DirectCoupledPDE::~DirectCoupledPDE()
{
  ENTER_FCN( "DirectCoupledPDE::~DirectCoupledPDE" );

  delete algsys_;
  delete solveStep_;
  
}


void DirectCoupledPDE::SetSinglePDEs( const StdVector<SinglePDE*> &pdes)
{
  ENTER_FCN( "DirectCoupledPDE::SetSinglePDEs" );

  singlePDEs_ = pdes;
}

void DirectCoupledPDE::SetCouplings( const StdVector<BasePairCoupling*> &couplings)
{
  ENTER_FCN( "DirectCoupledPDE::SetCouplings" );

  couplings_ = couplings;
}

void DirectCoupledPDE::Init(Integer bcSequenceStep ,
			    std::string  bcSequenceTag)
{
  ENTER_FCN( "DirectCoupledPDE::Init" );
  
  //std::cerr << "Entering DirectCoupledPDE::Init\n";

  bcSequenceIndex_ = bcSequenceStep;
  bcSequenceTag_ = bcSequenceTag;

  // Create algebraic system and pass it to SinglePDEs
  algsys_ = new StandardSystem();
  
  // Get parameter and report object of OLAS
  olasParams_ = algsys_->GetOLASParams();
  olasReport_ = algsys_->GetOLASReport();
  
  // define solveStep-driver
  DefineSolveStep();
  
  // activate direct coupling information
  // and initialize all single pdes
  for (Integer i=0; i<singlePDEs_.GetSize(); i++) {
    //std::cerr << "set Direct Coupling for pde " 
    //      << singlePDEs_[i]->GetName() << "\n";
    
    singlePDEs_[i]->SetDirectCoupling( algsys_, solveStep_ );
    
    // Initialize all SinglePDEs
    singlePDEs_[i]->Init( bcSequenceStep, bcSequenceTag );
  }

  // HARD CODED:
  // We simply take the assemble object of the first SinglePDE
  assemble_ = singlePDEs_[0]->getPDE_assemble();
 
  // Get information about number of dirichlet values,
  // dofs, constraints and needed matrices
  
  // Iterate over all single PDEs and collect data about
  // included boundary conditions
  numBuildInDirichletBCs_ = 0;
  NodeEQN * eqn;
  for ( Integer i=0; i<singlePDEs_.GetSize(); i++ ) {
    eqn = singlePDEs_[i]->getPDE_eqnData();
    numBuildInDirichletBCs_ += eqn->GetNumBuildInDirichletEQNs();
  }

  // get information about the needed matrix types 
  std::set<FEMatrixType> matTypes;
  for ( Integer i=0; i<singlePDEs_.GetSize(); i++ ) {
    singlePDEs_[i]->GetMatrixTypes( matTypes );
    matrixTypes_.insert( matTypes.begin(), matTypes.end() );
  }
    
  //std::cerr << "Leaving DirectCoupledPDE::Init\n";
}


void DirectCoupledPDE::WriteGeneralPDEdefines()
{

  ENTER_FCN( "DirectCoupledPDE::WriteGeneralPDEdefines" );
}

void DirectCoupledPDE::SetBCs(const Double atimestep)
{
  ENTER_FCN( "DirectCoupledPDE::SetBCs" );
  for (Integer i=0; i<singlePDEs_.GetSize(); i++) 
    singlePDEs_[i]->SetBCs( atimestep );
}


void DirectCoupledPDE::DefineAlgSys()
{
  ENTER_FCN( "DirectCoupledPDE::DefineAlgSys" );
 
  //std::cerr << "Entering DirectCoupledPDE::DefineAlgSys()\n";
 
  std::string pdeName;
  PdeIdType pdeId;
  NodeEQN *eqn = NULL;
  Integer numDirichletEQNs = 0;

  // Begin setup of the matrix graph
  //std::cerr << "algsys_ = " << &(*algsys_) << std::endl;
  algsys_->GraphSetupInit( singlePDEs_.GetSize() );

  // iterate over all singlePDE and register them
  for (Integer i=0; i<singlePDEs_.GetSize(); i++) {

    // obtain PDE identification tag from algebraic system
    // and set number of dirichlet and constraint equations
    pdeName= singlePDEs_[i]->GetName();
    eqn = singlePDEs_[i]->getPDE_eqnData();
    pdeId = algsys_->RegisterPDE( pdeName, eqn->GetNumEQNs() );
    //std::cerr << "Registering PDE " << pdeName << " with id= " 
    //      << pdeId  << "\n";
    singlePDEs_[i]->SetPDEId( pdeId);
    //numDirichletEQNs = singlePDEs_[i]->GetNumRestraints();
    //numDirichletEQNs -= eqn->GetNumBuildInDirichletEQNs();
    //std::cerr << "Setting Number of Dirichlet BCs to " 
    //      << numDirichletEQNs << std::endl;
    //algsys_->SetNumDirichletBCs( pdeId, numDirichletEQNs);
    //algsys_->SetNumDof( pdeId, eqn->GetDofsPerEQN() );

    // HARD CODED 
    // -> call CFSOLASPARAMS for each individual single PDE
    // lateron this has to be done by the directcoupled PDE directly
    singlePDEs_[i]->DefineAlgSys();
  }
  
  
    // iterate over all singlePDE and setup matrix graph
  for (Integer i=0; i<singlePDEs_.GetSize(); i++) {
    //std::cerr << "Setting up Matrix graph for PDE " 
    //<< singlePDEs_[i]->GetName() << std::endl;
    // trigger the creation and assembly of the matrix graph
    singlePDEs_[i]->SetupMatrixGraph();
    //    std::cerr << " DONE\n";
  }
  
  
  // Initialize all Coupling Objects and setup their matrix graph
  for (Integer i=0; i<couplings_.GetSize(); i++) {
    couplings_[i]->SetAlgSys( algsys_ );
    couplings_[i]->Init( bcSequenceIndex_, bcSequenceTag_ );
    //    std::cerr << "Setting up matrix graph of coupling " << i << "\n";
    couplings_[i]->SetupMatrixGraph();
    //    std::cerr << "DONE\n" << std::endl;
  }

  //std::cerr << "Call GraphSetupDon\n";
  // finish the assembly of the matrix graph
  algsys_->GraphSetupDone();
  
  // obtain reordering of the matrix graph for each PDE
  // and give it to the EQN-object
  for (Integer i=0; i<singlePDEs_.GetSize(); i++) {
    pdeId = singlePDEs_[i]->GetPDEId();
    eqn = singlePDEs_[i]->getPDE_eqnData();
    Integer * newOrder = algsys_->GetReordering( pdeId );
    //std::cerr << "Reorder PDE " << pdeId << std::endl;
    eqn->ReorderMapping( newOrder );
    delete[] newOrder;
  }
  
  //std::cerr << "Call 'CreateMatrices_Solver'\n";
  // Allocate the necessary matrices as well as solver and preconditioner
  CreateMatrices_Solver();
  
  //std::cerr << "Leaving DirectCoupledPDE::DefineAlgSys()\n";
}


Integer DirectCoupledPDE::GetNumRestraints() {
  ENTER_FCN( "DirectCoupledPDE::GetNumRestraints" );

  Integer totalNumRestraints = 0;
  for (Integer i=0; i<singlePDEs_.GetSize(); i++) 
    totalNumRestraints += singlePDEs_[i]->GetNumRestraints();
  
  return totalNumRestraints;
}
void DirectCoupledPDE::GetMatrixTypes( std::set<FEMatrixType> &matTypes) {
  ENTER_FCN( "DirectCoupledPDE::GetMatrixTypes" );
  Error ( "Not implemented", __FILE__, __LINE__ );
}


void DirectCoupledPDE::SaveSolution() {
  ENTER_FCN( "DirectCoupledPDE::SaveSolution" );

  Double *ptSol = NULL;
  Integer size = 0;
  BaseNodeStoreSol *ptNodeSol;

  for (Integer i=0; i<singlePDEs_.GetSize(); i++) {
    
    // get associated solution from algsys
    size = algsys_->GetSolutionVal( ptSol, singlePDEs_[i]->GetPDEId() );

    // set pointer to solution object of the PDE
    ptNodeSol = singlePDEs_[i]->getPDESolution();
    ptNodeSol->SetAlgSysDataPointer( size, ptSol );
    
  }
}

// ======================================================
// POSTPROC SECTION
// ======================================================

void DirectCoupledPDE::PostProcess()
{
  ENTER_FCN( "DirectCoupledPDE::PostProcess" );
  
  for (Integer i=0; i<singlePDEs_.GetSize(); i++) {
    singlePDEs_[i]->PostProcess();
  }
}


void DirectCoupledPDE::WriteResultsInFile(const Integer kstep,
					  const Double asteptime,
					  Integer stepOffset,
					  Double timeOffset)
{
  ENTER_FCN( "DirectCoupledPDE::WriteResultsInFile" );
  
  for (Integer i=0; i<singlePDEs_.GetSize(); i++) {
      singlePDEs_[i]->WriteResultsInFile( kstep, asteptime, 
					  stepOffset, timeOffset);
  }

}

// ======================================================
// COUPLING SECTION
// ======================================================

void DirectCoupledPDE::InitCoupling(PDECoupling * Coupling)
{
  ENTER_FCN( "DirectCoupledPDE::InitCoupling" );
  Warning( "DirectCoupledPDE::InitCoupling:Not implemented yet", __FILE__ , __LINE__ );
}
void DirectCoupledPDE::ResetCoupling()
{
  ENTER_FCN( "DirectCoupledPDE::ResetCoupling ");
  Warning( "Not implemented yet", __FILE__ , __LINE__ );
}
	   
void DirectCoupledPDE::CalcInputCoupling()
{
  ENTER_FCN( "DirectCoupledPDE::CalcInputCoupling" );
  Warning( "Not implemented yet", __FILE__ , __LINE__ );
}

void DirectCoupledPDE::CalcOutputCoupling()
{
  ENTER_FCN( "DirectCoupledPDE::CalcOutputCoupling" );
  Warning( "Not implemented yet", __FILE__ , __LINE__ );
}


// ======================================================
// METHODS FOR ASSEMBLING
// ======================================================


void DirectCoupledPDE::AssembleMatrices() {

  // Assembly of diagonal-matrices
  for (Integer i=0; i<singlePDEs_.GetSize(); i++) 
    {
      //std::cerr << "Assembling Matrices for PDE " 
      //<< singlePDEs_[i]->GetName() << std::endl;
      singlePDEs_[i]->AssembleMatrices();
    }

  // Assembly of off-diagonal entries (coupling objcts)
  for (Integer i=0; i<couplings_.GetSize(); i++) 
    {
      //      std::cerr << "Assembling Matrices for Coupling " 
      //<< couplings_[i]->GetName() << std::endl;
      couplings_[i]->AssembleMatrices();
    }
}

void DirectCoupledPDE::AssembleSrcRHS(const Double time) {
  ENTER_FCN( "DirectCoupledPDE::AssembleSrcRHS" );

  for (Integer i=0; i<singlePDEs_.GetSize(); i++) 
    {
      singlePDEs_[i]->AssembleSrcRHS( time );
    }
  
}

void DirectCoupledPDE::AssembleNLRHS(const Double time) {
  ENTER_FCN( "DirectCoupledPDE::AssembleNLRHS" );

  for (Integer i=0; i<singlePDEs_.GetSize(); i++) 
    {
      singlePDEs_[i]->AssembleNLRHS( time );
    }
}

void DirectCoupledPDE::AssembleSprings( const Double time) {
  ENTER_FCN( "DirectCoupledPDE::AssembleSprings" );

  for (Integer i=0; i<singlePDEs_.GetSize(); i++) 
    {
      singlePDEs_[i]->AssembleSprings(time );
    }
}

void DirectCoupledPDE::InitNonLinMatrices() {
  ENTER_FCN( "DirectCoupledPDE::InitNonLinMatrices" );

  for (Integer i=0; i<singlePDEs_.GetSize(); i++) 
    {
      singlePDEs_[i]->InitNonLinMatrices();
    }
}

void DirectCoupledPDE::SetupMatrixGraph() {
  // nothing to do here, since this method gets only called for 
  // SinglePDEs
}


void DirectCoupledPDE::DefineSolveStep() {
  ENTER_FCN( "DirectCoupledPDE::DefineSolveStep" );
  
  // HARD CODED
  solveStep_ = new SolveStepMech(*this);
}

} // end of namesapce


