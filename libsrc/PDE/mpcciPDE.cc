#include <fstream>
#include <iostream>
#include <sstream>
#include <math.h>

#include "DataInOut/Unverg/outUnverg.hh"
#include "DataInOut/GMV/outGMV.hh"
#include "Forms/forms_header.hh"
#include "Estimator/spaceerror.hh"
#include "DataInOut/WriteInfo.hh"
#include "Driver/assemble.hh"
#include "General/defs.hh"

#include <Matrix/matrix.hh>
#include <Utils/vector.hh>
#include <Forms/gradfieldop.hh>
#include "mpcciPDE.hh"
#include <General/defs.hh>
#include "DataInOut/ParamHandling/BaseParamHandler.hh" 
#include <string>
#include "Utils/StdVector.hh"

#ifdef MpCCI
#include <MpCCIcpl/MpCCIexch.hh>
#endif

#include "nodeEQN.hh"
#include "Driver/solveStepMpCCI.hh"

namespace CoupledField {

  // ***************
  //   Constructor
  // ***************
  MpcciPDE::MpcciPDE(Grid * aptgrid, BCs *aptbcs, TimeFunc *aptTimeFunc,
		   FileType *aptFileType, WriteResults *aptOut)
    :SinglePDE(aptgrid, aptbcs, aptFileType, aptOut, aptTimeFunc) 
{

    ENTER_FCN( "MpcciPDE::MpcciPDE" );

    // =====================================================================
    // set solution information
    // =====================================================================
    dofspernode_ = 1;
    solTypes_ = MPCCI;
    solDofs_ = 1;
    pdename_          = "mpcci";
    pdematerialclass_ = "piezo"; 
    flagFirstTimeStep_= TRUE;
    converged_=FALSE;
}

void MpcciPDE::Init(Integer bcSequenceIndex, std::string  bcSequenceTag)
{

  ENTER_FCN( "MpCCI::Init()" );

  bcSequenceIndex_ = bcSequenceIndex;
  bcSequenceTag_ = bcSequenceTag;

  // =====================================================================
  // get regions/subdomains for PDE
  // =====================================================================
  params->GetList( "name", subdoms_, pdename_, "region" );
  Info->PrintF( pdename_, " %s lives on regions:\n", pdename_.c_str());
  for ( Integer k = 0; k < subdoms_.GetSize(); k++ ) 
    {
      Info->PrintF( pdename_, " %s\n", subdoms_[k].c_str() );
    }
 
  eqnData_ = new ScalarNodeEQN( ptgrid_, ptBCs_, subdoms_, actlevel_, dofspernode_ );
  eqnData_->CalcMpcciMapping();
  numPDENodes_ = eqnData_->GetNumLocalNodes();
  numElems_ = eqnData_->GetNumLocalElems();
    
  PreparePDE4Computation();
}



void MpcciPDE::PreparePDE4Computation()
{
  ENTER_FCN( "MpcciPDE::PreparePDE4Computation" );

#ifdef MpCCI

  params->GetList( "name", subdoms_, pdename_, "region" );
    
  StdVector<Elem*> elemssd;

  if (subdoms_.GetSize() != 1)
    {
      std::cerr << "currently only one subdomain can be coupled with mpcciPDE" << std::endl;
      exit(0);
    }
  else
    {
      ptgrid_->GetElemSD(elemssd,subdoms_[0],actlevel_);
      ptgrid_->CalcNumberOfNodesInPatch(elemssd,mapSD_);

      MpCCInodes_= mapSD_.GetSize();
      ptMpCCIexch_ = new MpCCIexch(ptgrid_,MpCCInodes_);
    }

  ptMpCCIexch_->PutExchangeGrid2MpCCI(subdoms_, *eqnData_);
#endif

}


void MpcciPDE::DefineIntegrators(const Integer level)
{
  ENTER_FCN( "MpcciPDE::DefineIntegerators" );
}


void MpcciPDE::DefineSolveStep()
{
  ENTER_FCN( "MpcciPDE::DefineSolveStep" );
  
  solveStep_ = new SolveStepMpCCI(*this);
}


// ======================================================
// POSTPROCESSING SECTION
// ======================================================

void MpcciPDE::WriteResultsInFile(Integer stepOffset,
				 Double timeOffset)
{
  ENTER_FCN( "MpcciPDE::WriteResultsInFile" );
}

void MpcciPDE::PostProcess(const Integer level)
{
  ENTER_FCN( "MpcciPDE::PostProcess" );
}


// ======================================================
// COUPLING SECTION
// ======================================================

void MpcciPDE::InitCoupling(PDECoupling * Coupling)
{
  ENTER_FCN( "MpcciPDE::InitCoupling" );
  
  pdeIsCoupled_ = TRUE;
  ptCoupling_   = Coupling;

  const Integer numCouplings = ptCoupling_->GetNumOutputCouplings();
  

  nonLin_ = FALSE;

  // Initialization of coupling helper arrays
  std::string quantity;
  StdVector<Integer> * couplingnodes = NULL;
  StdVector<Elem*> interface_tmp;
  StdVector<StdVector<ShortInt> > isBoundaryNode_tmp;
  StdVector<std::string> * neighRegions = NULL;
  StdVector<StdVector<Integer> > elemNodeToCouplingNode_tmp;
  F_Interface_.Resize(numCouplings);
  isBoundaryNode_.Resize(numCouplings);
  elemNodeToCouplingNode_.Resize(numCouplings);

  for (Integer actCoupling=0; actCoupling<numCouplings; actCoupling++)
    {
      // Initialize arrays for coupling surface elements
      if (ptCoupling_->GetOutputQuantity(actCoupling) == FLUID_FORCE)
	{

	  
	  ptCoupling_->GetOutputNodes(actCoupling, couplingnodes);
	  if (couplingnodes == 0)
	    std::cerr << "Couplingnodes = 0!!!!" << std::endl;
	  
	  F_Interface_[actCoupling] = interface_tmp;

	  // Intialize the memory of the coupling values
	  ptCoupling_->CreateCouplingVector(actCoupling,isComplex_);

	  isBoundaryNode_tmp.Clear();
	  isBoundaryNode_tmp.Resize(interface_tmp.GetSize());
	  elemNodeToCouplingNode_tmp.Clear();
	  elemNodeToCouplingNode_tmp.Resize(interface_tmp.GetSize());
	 
	  
	  for (Integer ielem=0; ielem<interface_tmp.GetSize(); ielem++)
	    {
	      isBoundaryNode_tmp[ielem].Resize(interface_tmp[ielem]->connect.GetSize());
	      elemNodeToCouplingNode_tmp[ielem].Resize(interface_tmp[ielem]->connect.GetSize());

	      // Determine Boundary Nodes
	      for (Integer ielemnode=0; ielemnode<isBoundaryNode_tmp[ielem].GetSize(); ielemnode++)
		for (Integer inodes=0; inodes<(*couplingnodes).GetSize(); inodes++)
		  if (interface_tmp[ielem]->connect[ielemnode] == (*couplingnodes)[inodes] )
		    {
		      isBoundaryNode_tmp[ielem][ielemnode] = 1;
		      elemNodeToCouplingNode_tmp[ielem][ielemnode] = inodes;
		      break;
		    } // end if

	    } // end for (ielems)

	  isBoundaryNode_[actCoupling] = isBoundaryNode_tmp;
	  elemNodeToCouplingNode_[actCoupling]  = elemNodeToCouplingNode_tmp;
	} // end if
            
    } // end for (actNode)
}
  


void MpcciPDE::CalcInputCoupling()
{

  ENTER_FCN( "MpcciPDE::CalcInputCoupling" );

  std::string errMsg;
  StdVector<Integer> * nodes;
  CFSVector * val;
  Integer pdeNode, eqnNr,eqnDof;
  Integer couplingDof;
  Boolean clearCoords = TRUE;

  // Reset counter for boundary conditions
  couplingBCsCounter_ = 0;
  
  // Outer loop over all INPUT coupling terms
  for (Integer i=0; i<ptCoupling_->GetNumInputCouplings(); i++)
	{

	  //    ptCoupling_ = &ptCoupling_[i];
	  ptCoupling_->GetInputValues(i, val);
	  couplingDof = ptCoupling_->GetInputDof(i);
    
	  // Up to now, Coupling is only possible with
	  // Real valued solutions
	  Vector<Double> const & help = dynamic_cast<Vector<Double>&>(*val);
	  Vector<Double> displ;
	  switch(ptCoupling_->GetInputType(i))
		{
		  // -------------------
		  // COORDINATE COUPLING
		  // -------------------
		case COORD:
	  
		  ptCoupling_->GetInputNodes(i, nodes);

		  displ.Resize(nodes->GetSize() * ptCoupling_->GetInputDof(i) );
		  displ.Init(-1);

		  for (Integer j=0; j<nodes->GetSize(); j++)
		    {
			for (Integer dof=0; dof<ptCoupling_->GetInputDof(i); dof++)
			  {
			    pdeNode = eqnData_->Mesh2PDENode((*nodes)[j]);
			    //std::cerr << "pdeNode " << pdeNode << "=" << (*nodes)[j] << std::endl;
				if (pdeNode==-1) {
				  errMsg =  pdename_;
				  errMsg += "PDE: Coupling node Nr. ";
				  errMsg += Info->GenStr((*nodes)[j]);
				  errMsg += " is not in contained in list of my subdomains!";
				  Error(errMsg.c_str(), __FILE__, __LINE__);
				}
				displ[dof + (pdeNode-1)*dim_] = help[dof + j*dim_];
			  }
		    }
		  if (flagFirstTimeStep_== TRUE)
		    {
		      flagFirstTimeStep_=FALSE;
		    }
#ifdef MpCCI
		  else ptMpCCIexch_->CouplSendPhase(displ,converged_);
#endif
		  break;
		}  // end switch
	} // end for
}


void MpcciPDE::CalcOutputCoupling()
{
  ENTER_FCN( "MpcciPDE::CalcOutputCoupling" );

  SolutionType quantity;
  StdVector<Integer> * couplingNodes     = NULL;
  CFSVector * values = 0;
  Integer forcesCount = 0;

  // loop over all output coupling quantities
  for (Integer actCoupling=0; actCoupling<ptCoupling_->GetNumOutputCouplings(); actCoupling++)
    {
      quantity = ptCoupling_->GetOutputQuantity(actCoupling);
      ptCoupling_->GetOutputValues(actCoupling, values);

      Vector<Double> * temp = dynamic_cast<Vector<Double> *>(values);
      
      switch(ptCoupling_->GetOutputType(actCoupling))
	{
	  
	case NODE:	  
	  ptCoupling_->GetOutputNodes(actCoupling, couplingNodes);
	  
	  if (quantity == FLUID_FORCE)
	    {

#ifdef MpCCI
	      ptMpCCIexch_->CouplRecvPhase(temp[actCoupling]);
#endif

	      forcesCount++;
	    }
	  break;
	  
	case ELEM:
	  Error("No Element coupling output", __FILE__,__LINE__);
	}
    }
}


Boolean MpcciPDE::HasOutput(SolutionType output)
{
  ENTER_FCN( "MpcciPDE::HasOutput" );
  
  switch (output)
    {
    case FLUID_FORCE:
      return TRUE;
      break;
    default:
      return FALSE;
      break;
    }
  return FALSE;
}


void MpcciPDE::ReadStoreResults() {

  ENTER_FCN( "MpcciPDE::ReadStoreResults" );
}

} // end of namespace

