#include <fstream>
#include <iostream>
#include <sstream>
#include <math.h>

#include "DataInOut/Unverg/outUnverg.hh"
#include "DataInOut/GMV/outGMV.hh"
#include "Forms/forms_header.hh"
#include "Forms/linElecInt3D.hh"
#include "Forms/linElecInt2D.hh"
#include "Estimator/spaceerror.hh"
#include "DataInOut/WriteInfo.hh"
#include "Driver/assemble.hh"
#include "General/defs.hh"

#include <Matrix/matrix.hh>
#include <Utils/vector.hh>
#include <Forms/gradfieldop.hh>
#include "elecPDE.hh"
#include <General/defs.hh>
#include "DataInOut/ParamHandling/BaseParamHandler.hh" 
#include <string>
#include "Utils/StdVector.hh"
#include "Driver/solveStepElec.hh"
#include "CoupledPDE/pdecoupling.hh"

namespace CoupledField {

  // ***************
  //   Constructor
  // ***************
  ElecPDE::ElecPDE(Grid * aptgrid, BCs *aptbcs, TimeFunc *aptTimeFunc,
		   FileType *aptFileType, WriteResults *aptOut)
    :SinglePDE(aptgrid, aptbcs, aptFileType, aptOut, aptTimeFunc) {

    ENTER_FCN( "ElecPDE::ElecPDE" );

    // =====================================================================
    // set solution information
    // =====================================================================
    dofspernode_ = 1;
    solTypes_ = ELEC_POTENTIAL;
    solDofs_ = 1;
    pdename_          = "electrostatic";
    pdematerialclass_ = "piezo"; 
 
    geoUpdate_ = FALSE;
    nonLin_    = FALSE;
    isAlwaysStatic_ = TRUE;
    isPiezoCoupled_ = FALSE;

    //check, if problem is axisymmetric
    if ( params->HasValue( "type", "axi", "geometry" ) ) isaxi_ = TRUE;

    
    //\todo Is this variable needed?
    SolverCFS_ = FALSE;
  }
  


void ElecPDE::DefineIntegrators(const Integer level)
{
  ENTER_FCN( "ElecPDE::DefineIntegerators" );

  BaseForm *form;


  // if the pde is piezo-coupled, the electrostatic entries
  // have to multiplied with -1
  Double factor = 1.0;
  if ( isPiezoCoupled_ == TRUE )
    factor *= -1.0;  
  
  for (int actSD = 0; actSD < subdoms_.GetSize(); actSD++)
    {
      //reads eps33 (matrix notation starts with 0)
      Double eps33 = materialData_[actSD].GetPermittivity(2,2);
      
      if (dim_ == 3) {
	form = new linElecInt3D( materialData_[actSD] );
	form->SetFactor( factor );
      }
      else {
	form = new linElecInt2D( materialData_[actSD], isaxi_ );
	form->SetFactor( factor );
      }
	//form = new LaplaceInt( eps33*factor, isaxi_ );
      
      
      assemble_->AddIntegrator(form, subdoms_[actSD], SYSTEM, nonLin_);
    } // for
}

void ElecPDE::DefineSolveStep()
{
  ENTER_FCN( "ElecPDE::DefineSolveStep" );
  
  solveStep_ = new SolveStepElec(*this);
}



// ======================================================
// POSTPROCESSING SECTION
// ======================================================

void ElecPDE::WriteResultsInFile(const Integer kstep,
				 const Double asteptime,
				 Integer stepOffset,
				 Double timeOffset)
{
  ENTER_FCN( "ElecPDE::WriteResultsInFile" );

  Integer actStep = kstep + stepOffset;
  Double actTime = lasttimecalc_ + timeOffset;
  
#ifdef PARALLEL //only one thread should write output
  int commrank;
  MPI_Comm_rank(MPI_COMM_WORLD,&commrank);
  if (!commrank){
#endif
    
    // ATTENTION:
  // The errorMap should be assigned as a StoreSolution, not as a 
  // Vector. This is only temporarely
  ElemStoreSol<Double> error, error_Mesh;
  NodeStoreSol<Double> * solConverted;
  if (analysistype_ == STATIC)
    {
      solConverted = dynamic_cast<NodeStoreSol<Double>*>(sol_);
      
      
      // write electric potential
      if (saveSol_)
	outFile_->WriteNodeSolutionTransient(*solConverted, actStep, actTime);
      
      if (saveSolHist_)
	outFile_->WriteNodeHistoryTransient(*solConverted, actStep, actTime);
      
      if (calcEfield_.GetSize() !=0 )
	{
	  outFile_->WriteElemSolutionTransient(E_, actStep, actTime);
	}
      
      if (calcCharges_.GetSize() !=0 )
	{
	  outFile_->WriteElemSolutionTransient(charges_, actStep, actTime);
	}
      
      if (flags->CalcErrorMap_
	  )
	{
	  // this is only a temporar solution
	  error.SetNumSolutions(1);
	  error.SetNumElems(errorMap_.GetSize());
	  error.SetSolutionType(NO_SOLUTION_TYPE);
	  error.SetNumDofs(dofspernode_);
	  error.Init(0);
      
	  Error("Not implemented. Talk to Andreas and Elena", __FILE__, __LINE__);
	  
	  //Error.SetAlgSysVector(errorMap_);
	  
	  // ATTENTION!!
	  // up to now now transformation of the Error performed,
	  // since the calculation of the error is done on the global element numeration
	  //Error.TransformElemSolution(Error_Mesh,subdoms_,ptgrid_,actlevel_);
	  //OutFile_->WriteElemSolution(errorMap_, laststepcalc_, time, "relERR-E-Potential"); 
	  outFile_->WriteElemSolutionTransient(error_Mesh, actStep, actTime); 
	}
      
      if (calcEnergy_.GetSize() !=0 )
	CalcEnergy();
      
    }
  else
   Error("ElecPDE: Only static results can be written", __FILE__, __LINE__);


  // The following section was used by Gerhard to compute sum of forces over
  // different iteration and time steps. The sum was written into the .data stream.
  // Since this is not available anymore, this is commented out
#ifdef COMMENTET_OUT
  if (isIterCoupled_ == TRUE) {
    //   // TMPORARILY
    SolutionType quantity;
    StdVector<Integer> * couplingNodes     = NULL;
    CFSVector * values = 0;
    Vector<Double> sumForces(dim_);
    sumForces.Init();
    
    
    
    // loop over all output coupling quantities
    for (Integer actCoupling=0; actCoupling<ptCoupling_->GetNumOutputCouplings(); actCoupling++)
      {
	quantity = ptCoupling_->GetOutputQuantity(actCoupling);
	ptCoupling_->GetOutputValues(actCoupling, values);
	
	Vector<Double> const & temp = dynamic_cast<Vector<Double> &>(*values);
	switch(ptCoupling_->GetOutputType(actCoupling))
	  {
	    
	  case NODE:	  
	    ptCoupling_->GetOutputNodes(actCoupling, couplingNodes);
	    if (quantity == ELEC_FORCE_VWP)
	      {
		for (Integer iDof=0; iDof<dim_; iDof++)
		  for (Integer iNode=0; iNode<couplingNodes->GetSize(); iNode++)
		    sumForces[iDof] += temp[iNode*dim_ + iDof];
		
		*data << lasttimecalc_ << "\t";
		for (Integer i=0; i<dim_; i++)
		  *data << sumForces[i]<< "\t";
		
		*data << std::endl;
		
	      }
	    break;

	  case ELEM:
	    Error( "Element input coupling not implemented for elecPDE",
		   __FILE__, __LINE__ );
	  } // switch
      } // for
  }
#endif

#ifdef PARALLEL
    }//!commrank
#endif

}

void ElecPDE::PostProcess(const Integer level)
{
  ENTER_FCN( "ElecPDE::PostProcess" );


  // *** Electric Field Intensity ***
  if (calcEfield_.GetSize() !=0 ){
    CalcElectricField();
  }
  
  // *** Electric Charges ***
  if (calcCharges_.GetSize() !=0 ) {
    CalcCharges();
  }
}

void ElecPDE::CalcElectricField()
{
  ENTER_FCN( "ElecPDE::PostProcess" );
  
  NodeStoreSol<Double> & solhelp = dynamic_cast<NodeStoreSol<Double>&>(*sol_);
  //  GradientFieldOp * FieldOp = new GradientFieldOp(ptgrid_, this, eqnData_,
  //					  solhelp, ELEC_POTENTIAL,
  //					  actlevel_, isaxi_);


  GradientFieldOp<Double> * FieldOp = new GradientFieldOp<Double>(ptgrid_, this, eqnData_,
      						      solhelp, ELEC_POTENTIAL,
      						      actlevel_, isaxi_);
  
  Vector<Double> LCoord;
  
  StdVector<Elem*> elemssd;
  Integer counterElems=0;
  Vector<Double> TempE;
  Integer pdeElem;
  
  // loop over all subdomains
  for (Integer isd=0; isd<calcEfield_.GetSize(); isd++)
    {

      // ------ Calculation of the electric field ------

      //     Vector<Double> LCoord;
      // LCoord.Resize(dim_);
      // LCoord[0] = 0;
      // LCoord[1] = 0;
      // get vector of Elem of subdomain with color: subdoms[isd]
      ptgrid_->GetElemSD(elemssd,calcEfield_[isd],actlevel_);
      
      // loop over elements of subdomain
      for (Integer iel=0; iel< elemssd.GetSize(); iel++,counterElems++)
	{
	  elemssd[iel]->ptElem->GetCoordMidPoint(LCoord);
	  FieldOp->CalcElemGradField( TempE, elemssd[iel], LCoord, 1); 
	  pdeElem = eqnData_->Mesh2PDEElem(elemssd[iel]->elemNum);
	  E_.SetElemResult(pdeElem-1,TempE);
	}
    }
  delete FieldOp;
}

void ElecPDE::CalcCharges()
{
  ENTER_FCN( "ElecPDE::CalcCharges" );

  NodeStoreSol<Double> * solhelp = dynamic_cast<NodeStoreSol<Double>*>(sol_);
  StdVector<Elem*> surfElems, volElems;

  Vector<Double> lCoordSurf, lCoordVol, elemDField, normal;
  BaseFE * ptSurfElem, * ptVolElem;
  Double permittivity = 0.0;
  Double elemNormalD = 0.0;
  Double charge = 0.0;
  Double sumOfCharges = 0.0;
  Integer pdeElemNum = 0;
 
 

  
  // Create vector with interpolation coordinate.
  // For simplicity we only evaluate the integral
  // in coordinate origin
  lCoordSurf.Resize(dim_-1);
  lCoordSurf.Init(0);
  
  // Create operator for electric flux density and charge calculation
  GradientFieldOp<Double> * dFieldOp = new GradientFieldOp<Double>(ptgrid_, this, eqnData_, *solhelp,
				 ELEC_POTENTIAL, actlevel_, isaxi_);
  ElecChargeOp * chargeOp = new ElecChargeOp(ptgrid_, this, eqnData_, actlevel_, isaxi_);
  
  // loop over all subdomains
  for (Integer iSD=0; iSD<calcCharges_.GetSize(); iSD++){
    
    // get surface and acoording volume elements
    if (dim_ == 3)
      surfElems = ptBCs_->getFacesBC(calcCharges_[iSD], actlevel_);
    else if (dim_ == 2)
      surfElems = ptBCs_->getEdgesBC(calcCharges_[iSD], actlevel_);
    
    // get neighbouring volume elements of
    // surface elements
    ptgrid_->GetVolNeighboursForSurf(surfElems,chargeNeighborRegion_,
				     volElems, actlevel_);
    
    // loop over all surface elements
    for (Integer iElem=0; iElem<surfElems.GetSize(); iElem++)
      {
	
	ptSurfElem = surfElems[iElem]->ptElem;
	ptVolElem = volElems[iElem]->ptElem;
	const StdVector<Integer> & surfConnect = surfElems[iElem]->connect;
	const StdVector<Integer> & volConnect = volElems[iElem]->connect;
	
	// calculate volume integration coordinates from
	// surfe integration coordinat for evalauting the 
	// electric flux density on the surface of the volume
	// element
	ptSurfElem->GetCoordMidPoint(lCoordSurf);
	ptVolElem->GetLocalIntPoints4Surface(surfConnect, volConnect,
					     lCoordSurf, lCoordVol);


	//std::cerr << "Volume coordinates are:\n" << lCoordVol << std::endl;
	// Get the right material parameter for actual volume element
	for (Integer i=0; i<subdoms_.GetSize(); i++)
	  {
	    if (subdoms_[i] == volElems[iElem]->namesd)
	      permittivity  = materialData_[i].GetPermittivity(2,2);
	  }
	
	// Calc electric flux density
	dFieldOp->CalcElemGradField(elemDField, volElems[iElem], 
				    lCoordVol,permittivity);
	
	// Calc global normal
	ptgrid_->CalcSurfNormalOutOfVol(normal, *surfElems[iElem], 
					*volElems[iElem]);
	normal *= -1;
	//std::cerr << "new normal =\n " << normal << std::endl << std::endl;
	// Calc normal flux density component
	Vector<Double> elemNoralD;

	// Since the routine CalcLinNormal always computes a normal
	// which points in the OPPOSITE direction of the volume elment,
	// we have to multiply the normal with -1 to get the correct sign for
	// the charges
	elemNormalD =  normal * elemDField;
	
	chargeOp->CalcElemCharge(charge, surfElems[iElem], 
				 lCoordSurf, elemNormalD);

	pdeElemNum = eqnData_->Mesh2PDEElem(volElems[iElem]->elemNum);
	
	// Create temporar vector, since SetElemResult only
	// can handle these
	Vector<Double> chargeVec(1);
	chargeVec[0] = charge;
	sumOfCharges +=charge ;
	charges_.SetElemResult(pdeElemNum-1, chargeVec);
	
      }
  }
  std::string outstring = "Sum of electric charges:\n";
  outstring += Info->GenStr(sumOfCharges);
  Info->PrintF(pdename_, outstring.c_str());

  delete chargeOp;
  delete dFieldOp;
}

  


void ElecPDE::CalcNodeForce(Vector<Double> & force, 
			    StdVector<Integer> & nodes, 
			    StdVector<Elem*> & elems,
			    StdVector<StdVector<ShortInt> > & isBoundaryNode,
			    StdVector<StdVector<Integer> > & elemNodeToCouplingNode)
{
  ENTER_FCN( "ElecPDE::CalcNodeForce" );

  NodeStoreSol<Double> * solhelp = dynamic_cast<NodeStoreSol<Double> *>(sol_);

  
  ElecForceOp * ForceOp; // = new ElecForceOp(ptgrid_, this, eqnData_, *solhelp, actlevel_, isaxi_);
   
  ElemStoreSol<Double> force_temp;
  
  force.Init(0.0);
  
  for (Integer ielem=0; ielem<elems.GetSize(); ielem++)
    {
      // Get Material Parameter
      Double epsilon;
      
      for (Integer i=0; i<subdoms_.GetSize(); i++)
	{
	  if (elems[ielem]->namesd == subdoms_[i]) 
	  epsilon = materialData_[i].GetPermittivity(2,2); 
	}
      
      // Only for testing
      //epsilon = 8.854e-12;
      
      ForceOp->CalcElemElecForce( force_temp, elems[ielem], epsilon, isBoundaryNode[ielem]);
      

      // Add the element force to the according coupling node
      for (Integer ielemnode=0; ielemnode<elems[ielem]->connect.GetSize(); ielemnode++)
	for( Integer idim=0; idim<dim_; idim++) {
	  force[elemNodeToCouplingNode[ielem][ielemnode]*dim_+idim] += 
	    force_temp(ielemnode,idim);
	}
    }


  Vector<Double> sum;
  sum.Resize(dim_);
  
  for (Integer i=0; i<nodes.GetSize(); i++)
    for (Integer dim=0; dim<dim_; dim++)
      sum[dim] += force[i*dim_+dim];

  delete ForceOp;


  //std::cerr << "Sum of E-Force:" << std::endl << sum << std::endl;
  
  // write information in .info-file
  Info->PrintF(pdename_, "Sum of electrostatic force (VWM):\n");
  Info->PrintVec(sum);
      
}





void ElecPDE::CalcEnergy()
{
  ENTER_FCN( "ElecPDE::CalcEnergy" );

  Matrix<Double> elemmat;  
  Matrix<Double> ptCoord;
  BaseFE         * ptElem;
  Double totalE=0;

  StdVector<Integer> connecth;
  Vector<double> help;

  Integer i, j;
  Vector<Double> energy(subdoms_.GetSize());

  for (i=0; i<subdoms_.GetSize(); i++)
    {
      //reads eps33 (matrix notation starts with 0)
      Double eps33 = materialData_[i].GetPermittivity(2,2);

      StdVector<Elem*> elemssd;
      ptgrid_->GetElemSD(elemssd,subdoms_[i],actlevel_);

      energy[i] = 0;
      for (j=0; j < elemssd.GetSize(); j++)
	{  
	  ptElem=elemssd[j]->ptElem;
	  BaseForm * bilinear_stiff = new LaplaceInt(ptElem, eps33, isaxi_);

	  connecth=elemssd[j]->connect;
	  GetElemCoords(connecth, ptCoord, actlevel_);
	  bilinear_stiff->CalcElementMatrix(ptCoord, elemmat);

	  // Mape Mesh to PDE node numbers
	  //Mesh2PDENode(connect_PDE,connecth,mesh2PDENode_);

// 	  EqnData_->Mesh2Eqn(Eqns,connecth);
// 	  (*debug) << "Nodes:" << connecth << std::endl;
// 	  (*debug) << "Eqns :" << Eqns << std::endl;


	  Vector<Double> elpot;
	  //GetSolOfElement(elpot, connect_PDE);	 
	  sol_->GetElemSolution(elpot, connecth);
	  help = elemmat * elpot;
	  energy[i] += help * elpot;

	  delete bilinear_stiff;	  
	}  

      totalE += 0.5*energy[i];
    }

  std::string unit = "(Ws)";
  std::string resulttype = "Electric Energy";
  std::string analysis;
  Double analysisVal;
  if ( analysistype_ == HARMONIC ) {
    analysis    = "Frequency:";
    analysisVal = actFrequency_;
  }
  else {
    analysis    = "Time:";
    analysisVal = lasttimecalc_;
  }

  Info->WriteResult(pdename_,  resulttype, subdoms_, energy, unit,
		    analysis, analysisVal);


  StdVector<std::string> suball(1);
  Vector<Double> tmp(1);
  suball[0] = "Summe";
  tmp[0] = totalE;
  Info->WriteResult(pdename_,  resulttype, suball, tmp, unit,
		    analysis, analysisVal);

}




// ======================================================
// COUPLING SECTION
// ======================================================



void ElecPDE::InitCoupling(PDECoupling * Coupling)
{
  ENTER_FCN( "ElecPDE::InitCoupling" );
  
  isIterCoupled_ = TRUE;
  ptCoupling_   = Coupling;

  const Integer numCouplings = ptCoupling_->GetNumOutputCouplings();
  

  nonLin_ = FALSE;

  // Initialization of coupling helper arrays
  std::string quantity;
  StdVector<Integer> * couplingnodes = NULL;

  for (Integer actCoupling=0; actCoupling<numCouplings; actCoupling++) {
    // Initialize arrays for coupling surface elements
    if (ptCoupling_->GetOutputQuantity(actCoupling) == ELEC_FORCE_VWP
	|| ptCoupling_->GetOutputQuantity(actCoupling) == ELEC_INTERFACE_FORCE) {
      
      ptCoupling_->GetOutputNodes(actCoupling, couplingnodes);
      if (couplingnodes == NULL)
	std::cerr << "Couplingnodes = 0!!!!" << std::endl;
      
      // if quantity is elecFocrceVWP, get volume neighbours lying next to
      // coupling nodes, because these volume elements have to be 
      // moved 'virtually'
      if (ptCoupling_->GetOutputQuantity(actCoupling) == ELEC_FORCE_VWP) {
	NodeStoreSol<Double> * solhelp = dynamic_cast<NodeStoreSol<Double> *>(sol_);
	ForceOp_ = new  ElecForceOp(ptgrid_, this, eqnData_, *solhelp, dim_, materialData_,
				    subdoms_, actlevel_, isaxi_);
	ForceOp_->Setup(subdoms_, *couplingnodes);
      }
      
      else if (ptCoupling_->GetOutputQuantity(actCoupling) == ELEC_INTERFACE_FORCE) {
	Error("Currently ELEC_INTERFACE_FORCE not supported");
      }
      
      else {
	Enum2String(ptCoupling_->GetOutputQuantity(actCoupling), quantity);
	std::string errMsg = "Coupling " + quantity +  " not known! ";	  
	Error(errMsg.c_str(), __FILE__,__LINE__);
      }
      
      // Intialize the memory of the coupling values
      ptCoupling_->CreateCouplingVector(actCoupling,isComplex_);
      
      
    } // end for (actNode)
     
  } 


//   // Initialization of coupling helper arrays
//   std::string quantity;
//   StdVector<Integer> * couplingnodes = NULL;
//   StdVector<Elem*> interface_tmp;
//   StdVector<StdVector<ShortInt> > isBoundaryNode_tmp;
//   StdVector<std::string> * neighRegions = NULL;
//   //StdVector<Integer> numBoundaryNodes_tmp;
//   StdVector<StdVector<Integer> > elemNodeToCouplingNode_tmp;
//   F_Interface_.Resize(numCouplings);
//   isBoundaryNode_.Resize(numCouplings);
//   elemNodeToCouplingNode_.Resize(numCouplings);

//   for (Integer actCoupling=0; actCoupling<numCouplings; actCoupling++)
//     {
//       // Initialize arrays for coupling surface elements
//       if (ptCoupling_->GetOutputQuantity(actCoupling) == ELEC_FORCE_VWP
// 	  || ptCoupling_->GetOutputQuantity(actCoupling) == ELEC_INTERFACE_FORCE)
// 	{

	  
// 	  ptCoupling_->GetOutputNodes(actCoupling, couplingnodes);
// 	  if (couplingnodes == NULL)
// 	    std::cerr << "Couplingnodes = 0!!!!" << std::endl;
	  
// 	  // if quantity is elecFocrceVWP, get volume neighbours lying next to
// 	  // coupling nodes, because these volume elements have to be 
// 	  // moved 'virtually'
// 	  if (ptCoupling_->GetOutputQuantity(actCoupling) == ELEC_FORCE_VWP) {
// 	    ptCoupling_->GetOutputNeighbourRegion(actCoupling, neighRegions);
// 	    //	    ptgrid_->GetInterfaceNeighbours(*couplingnodes, *neighRegions, 
// 	    ptgrid_->GetInterfaceNeighbours(*couplingnodes, subdoms_, 
// 					    interface_tmp, actlevel_);
// 	  }
// 	  else if (ptCoupling_->GetOutputQuantity(actCoupling) == ELEC_INTERFACE_FORCE)
// 	    {
// 	      // help construction for correct assignement of predefined values ... :O(
// 	      StdVector<Elem*>* interface_tmp_Ptr;
// 	      ptCoupling_->GetOutputNeighbourElems(actCoupling, interface_tmp_Ptr);
// 	      interface_tmp = *interface_tmp_Ptr;
// 	    }
// 	    else 
// 	      {
// 		Enum2String(ptCoupling_->GetOutputQuantity(actCoupling), quantity);
// 		std::string errMsg = "Coupling " + quantity +  " not known! ";	  
// 		Error(errMsg.c_str(), __FILE__,__LINE__);
// 	      }
	  
// 	  F_Interface_[actCoupling] = interface_tmp;

// 	  // Intialize the memory of the coupling values
// 	  ptCoupling_->CreateCouplingVector(actCoupling,isComplex_);
	  

// 	  isBoundaryNode_tmp.Clear();
// 	  isBoundaryNode_tmp.Resize(interface_tmp.GetSize());
// 	  elemNodeToCouplingNode_tmp.Clear();
// 	  elemNodeToCouplingNode_tmp.Resize(interface_tmp.GetSize());
	 
	  
// 	  for (Integer ielem=0; ielem<interface_tmp.GetSize(); ielem++)
// 	    {
// 	      isBoundaryNode_tmp[ielem].Resize(interface_tmp[ielem]->connect.GetSize());
// 	      elemNodeToCouplingNode_tmp[ielem].Resize(interface_tmp[ielem]->connect.GetSize());

// 	      // Determine Boundary Nodes
// 	      for (Integer ielemnode=0; ielemnode<isBoundaryNode_tmp[ielem].GetSize(); ielemnode++)
// 		for (Integer inodes=0; inodes<(*couplingnodes).GetSize(); inodes++)
// 		  if (interface_tmp[ielem]->connect[ielemnode] == (*couplingnodes)[inodes] )
// 		    {
// 		      isBoundaryNode_tmp[ielem][ielemnode] = 1;
// 		      elemNodeToCouplingNode_tmp[ielem][ielemnode] = inodes;
// 		      break;
// 		    } // end if

// 	    } // end for (ielems)

// 	  isBoundaryNode_[actCoupling] = isBoundaryNode_tmp;
// 	  elemNodeToCouplingNode_[actCoupling]  = elemNodeToCouplingNode_tmp;
// 	} // end if
            
//     } // end for (actNode)

}
  


void ElecPDE::CalcOutputCoupling()
{
  ENTER_FCN( "ElecPDE::CalcOutputCoupling" );

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
	  
	  if (quantity == ELEC_POTENTIAL)
	    sol_->NodeSolutionToCoupling(*values, *couplingNodes);
	    
	  if (quantity == ELEC_FORCE_VWP)
	    {
	      Vector<Double> totalForce;
	      ForceOp_->CalcNodeForce(*temp, totalForce);

	      // write information in .info-file
	      Info->PrintF(pdename_, "Sum of electrostatic force (VWM):\n");
	      Info->PrintVec(totalForce);

// 	      CalcNodeForce(*temp, 
// 			    *couplingNodes, 
// 			    F_Interface_[forcesCount], 
// 			    isBoundaryNode_[forcesCount], 
// 			    elemNodeToCouplingNode_ [forcesCount]);
	      
	      forcesCount++;
	    }


	  if (quantity == ELEC_INTERFACE_FORCE)
	    CalcInterfaceForces(actCoupling);
	    	  
	  break;
	  
	case ELEM:
	  Error("No Element coupling output", __FILE__,__LINE__);
	}
    }
}


Boolean ElecPDE::HasOutput(SolutionType output)
{
  ENTER_FCN( "ElecPDE::HasOutput" );
  
  switch (output)
    {
    case ELEC_FORCE_VWP:
      return TRUE;
      break;
    case ELEC_POTENTIAL:
      return TRUE;
      break;
    case ELEC_FIELD_INTENSITY:
      return TRUE;
      break;
    case ELEC_INTERFACE_FORCE:
      return TRUE;
      break;
    default:
      return FALSE;
      break;
    }
  return FALSE;
}


void ElecPDE::SetPiezoCoupling()
{
  ENTER_FCN( "ElecPDE::SetPiezoCoupling" );
  
  isPiezoCoupled_ = TRUE;

}


void ElecPDE::CalcInterfaceForces(Integer actCoupling)
{
  ENTER_FCN( "ElecPDE::CalcInterfaceForces" );

  StdVector<Integer> *      couplingNodes          = NULL;
  CFSVector  *              elemCouplingSolsTemp   = NULL;
  StdVector<Elem*> *        couplingElems          = NULL;
  StdVector<Elem*> *        outerInterfaceVolElems = NULL;
  StdVector<Elem*> *        innerInterfaceVolElems = NULL;
  StdVector<MaterialData*>* innerCouplingMaterials = NULL;
  StdVector<MaterialData*>* outerCouplingMaterials = NULL;

  ptCoupling_->GetOutputNodes    (actCoupling, couplingNodes);
  ptCoupling_->GetOutputValues   (actCoupling, elemCouplingSolsTemp);
  ptCoupling_->GetOutputElements (actCoupling, couplingElems);
  ptCoupling_->GetOwnMaterials   (actCoupling, innerCouplingMaterials);
  ptCoupling_->GetOppositeMaterials (actCoupling, outerCouplingMaterials);
  ptCoupling_->GetOutputNeighbourElems(actCoupling, innerInterfaceVolElems);
  ptCoupling_->GetInputNeighbourElems (actCoupling, outerInterfaceVolElems);
 
  Integer couplingDof = ptCoupling_->GetOutputDof(actCoupling);

  Vector<Double> * elemCouplingSols = dynamic_cast<Vector<Double> *>(elemCouplingSolsTemp);
   
  elemCouplingSols->Init(0.0);

  Vector<Double> xPosCoupleNode(couplingNodes->GetSize());
  

  for (Integer actElem=0; actElem < couplingElems->GetSize(); actElem++)
    {
      Elem * actCoupleElem     = (*couplingElems)[actElem];
      Elem * actVolElem        = (*innerInterfaceVolElems)[actElem];
      // BaseFE * ptVolElem       = actVolElem->ptElem;
      BaseFE * ptCoupleElem    = actCoupleElem->ptElem;

      const Vector<Double> & intWeights = ptCoupleElem->GetIntWeights();  
      

      StdVector<Integer> connecth = actCoupleElem->connect;
      const Integer spaceDim   = 2;
      
      if (ptCoupleElem->GetDim() == 3)
	Error("CalcInterfaceForces only implemented for 2D!", __FILE__, __LINE__);


      Matrix<Double> ptCoord; 
      //StdVector<Integer> connect_PDE;

      GetElemCoords(connecth, ptCoord, actlevel_);



      // get correct permittivity belonging to the neighbouring element of the interface
      Double eps1 = (*innerCouplingMaterials)[actElem]->GetPermittivity(2,2);
      Double eps2 = (*outerCouplingMaterials)[actElem]->GetPermittivity(2,2);


      // the normal vector points outwards of the MECHANICAL domain
      // (see. Kaltenbacher, "Num. Sim. of Mechatr. Act. & Sens." chapter 8.2)
      Vector<Double> n;


      // Note: The following line has to be replaced by a call to 
      ptgrid_->CalcSurfNormalOutOfVol(n, *actCoupleElem, *(*innerInterfaceVolElems)[actElem]); 
      //CalcLineNormalVec(n, *actCoupleElem, *(*innerInterfaceVolElems)[actElem]); // points outward own domain



      StdVector<Integer> boundNodesOfVolElem;
  
      // Integer ptCount = 0;
      for (Integer nNode=0; nNode < isBoundaryNode_[actCoupling][actElem].GetSize(); nNode++)
	if (isBoundaryNode_[actCoupling][actElem][nNode] == 1)
	  boundNodesOfVolElem.Push_back(nNode);

      if (boundNodesOfVolElem.GetSize() != 2)
	Error("In CalcInterfaceForces boundary nodes != 2!", __FILE__, __LINE__);
      


      // "interfaceForceVec" holds the absolute value of the forces on every node of an interface vector.
      // To establish the final force vectors, every force in every node has to be multiplied by
      // the normal vector of the interface element
      Vector<Double> interfaceForceOnNodes(connecth.GetSize());   // is automatically initialized by 0
      

      for (Integer actIP=1; actIP <= ptCoupleElem->GetNumIntPoints(); actIP++)
	{
	  Vector<Double> shapeFncAtIP;
	  Vector<Double> coordAtIP;
	  ptCoupleElem->GetShFncAtIp(shapeFncAtIP, actIP);
	  coordAtIP = ptCoord * shapeFncAtIP;
	  double jacDet = ptCoupleElem->CalcJacobianDetAtIp(actIP, ptCoord);
	  

	  Vector<Double> tempE;       // value of a field at lCoord
	  CalcEfieldAtCoupleElemIP(actVolElem, actCoupleElem, coordAtIP, boundNodesOfVolElem, tempE);


	  Double abs_E_normal = 0;
	  
	  for (Integer actSpaceDim=0; actSpaceDim < spaceDim; actSpaceDim++)
	    abs_E_normal += tempE[actSpaceDim] * n[actSpaceDim];
	  
	  Vector<Double> E_tangential(tempE);
	  E_tangential -= n * abs_E_normal;
	  Double abs_E_tangential = E_tangential.NormL2();
	  
	  // D is calculated in region 1 (D = E_1 * eps1), see Kaltenbacher: "Num. Sim ... " p. 136
	  Double interfaceForce = sqr(abs_E_normal * eps1 ) / 2 * (1/eps2 - 1/eps1) 
	    + sqr(abs_E_tangential) * (eps1 - eps2) / 2;

 	  if (isaxi_)
	    jacDet *= 2 * PI * coordAtIP[0];


	  Vector<Double> partInterfaceForceOnNodes;   
	  partInterfaceForceOnNodes =  shapeFncAtIP;
	  partInterfaceForceOnNodes *= interfaceForce * jacDet * intWeights[actIP-1];

	  interfaceForceOnNodes += partInterfaceForceOnNodes;  // force on every node! It is now vector yet!!	  
	}


      


      // copy result into final solution vector
      for (Integer actNode=0; actNode < ptCoord.GetSizeRow(); actNode++)
	{
	  Integer nodePos = 0;
	  
	  while(connecth[actNode] != (*couplingNodes)[nodePos] && 
		nodePos < couplingNodes->GetSize()) 
	    nodePos++;
	  
	  for (Integer actDof=0; actDof < couplingDof ; actDof++)
	    (*elemCouplingSols)[nodePos*couplingDof+actDof]  += 
	      interfaceForceOnNodes[actNode] * n[actDof];

	  xPosCoupleNode[nodePos] = ptCoord[0][actNode];
	}
    }

  Vector<Double> sum;
  sum.Resize(dim_);
  
  Integer k = 0;
  for (Integer i=0; i<elemCouplingSols->GetSize(); i++)
    for (Integer dim=0; dim<dim_; dim++)
      {
	sum[dim] += (*elemCouplingSols)[k];
	k++;
      }

  Info->PrintF(pdename_, "Sum of electrostatic force (Interface):\n");
  Info->PrintVec(sum); 

}




void ElecPDE::CalcEfieldAtCoupleElemIP(Elem * actVolElem,
				       Elem * actCoupleElem,
				       Vector<Double>& coordAtIP, 
				       StdVector<Integer>& boundNodesOfVolElem,
				       Vector<Double>& tempE)
{
  ENTER_FCN( "ElecPDE::CalcEfieldAtCoupleElemIP" );

  BaseFE * ptVolElem    = actVolElem->ptElem;
  BaseFE * ptCoupleElem = actCoupleElem->ptElem;
    
  const Integer spaceDim   = 2;  


  Matrix<Double> boundElLocCoords;
  ptCoupleElem->GetLocalCornerCoords(boundElLocCoords);  //coords: x:number, y:Dim
  
  Matrix<Double> volElLocCoords;
  ptVolElem->GetLocalCornerCoords(volElLocCoords); //coords: x:number, y:Dim
  
  
  
  // relPosIP = normed distance in the range of 0..1 between IP 1 and coord 1 of 1D element
  const Double locElemLen = boundElLocCoords[0][1] - boundElLocCoords[0][0];
  const Double relPosIP = (coordAtIP[0] - boundElLocCoords[0][0]) / locElemLen;
  
  
  Vector<Double> lCoord(spaceDim); // coords, at which the e-field has to be calculated
  
  
  Double volCoord1X = volElLocCoords[0][ boundNodesOfVolElem[0] ];
  Double volCoord2X = volElLocCoords[0][ boundNodesOfVolElem[1] ];
  Double volCoord1Y = volElLocCoords[1][ boundNodesOfVolElem[0] ];
  Double volCoord2Y = volElLocCoords[1][ boundNodesOfVolElem[1] ];
  
  // local x-coord of integration point of boundary element
  lCoord[0] = volCoord1X + relPosIP * (volCoord2X - volCoord1X);
  //y-coord
  lCoord[1] = volCoord1Y + relPosIP * (volCoord2Y - volCoord1Y);
  
  NodeStoreSol<Double> *solTemp = dynamic_cast<NodeStoreSol<Double>*>(sol_);

  GradientFieldOp<Double> * elecFieldOp = new  GradientFieldOp<Double>(ptgrid_, this, eqnData_, *solTemp, ELEC_POTENTIAL,
			      actlevel_, isaxi_);


  elecFieldOp->CalcElemGradField(tempE, actVolElem, lCoord, 1);

}


// ***********************************************************************
//   Obtain information on desired output quantities from parameter file
// ***********************************************************************
void ElecPDE::ReadStoreResults() {

  ENTER_FCN( "ElecPDE::ReadStoreResults" );

  // Construct vectors for restricted parameter search
  StdVector<std::string> keyVec;
  StdVector<std::string> attrVec;
  StdVector<std::string> valVec;
  std::string quantity;

  // *****************************
  // Determine nodal results
  // *****************************

  // --- Electric Potential ---
  StdVector<std::string> nodeValues;
  Enum2String(ELEC_POTENTIAL, quantity);
  keyVec  = pdename_, "storeResults", "nodeResults", "region";
  attrVec = "", "", "type";
  valVec = "", "", quantity;
  params->GetList( keyVec, attrVec, valVec, nodeValues);
  if (nodeValues.GetSize() > 0) {
    saveSol_ = TRUE;
   	hasOutput_ = TRUE;
  }
  
  // *****************************
  // Determine element results
  // *****************************
  keyVec  = pdename_, "storeResults", "elemResults", "region";
  attrVec = "", "", "type";

  // --- Electric Field Intensity ---
  Enum2String(ELEC_FIELD_INTENSITY, quantity);
  valVec  = "", "", quantity;
  params->GetList( keyVec, attrVec, valVec, calcEfield_ );

  // If the the symbolic name is "all" compute electric field for all regions
  if ( calcEfield_.GetSize() == 1 && calcEfield_[0] == "all" ) {
    calcEfield_ = subdoms_;
  }

  if ( calcEfield_.GetSize() > 0 ) {
   	hasOutput_ = TRUE;
    Info->PrintF( pdename_, " Computing electric field for regions:\n" );
    for ( Integer k = 0; k < calcEfield_.GetSize(); k++ ) {
      Info->PrintF( pdename_, " %s", calcEfield_[k].c_str() );
    }
    E_.SetNumSolutions(1);
    E_.SetSolutionType(ELEC_FIELD_INTENSITY);
    E_.SetNumElems(numElems_);
    E_.SetNumDofs(dim_);
    E_.SetPtrEQNData(eqnData_, ptgrid_, actlevel_);
    E_.Init(); 
  }

  // --- Electric Energy ---
  Enum2String(ELEC_ENERGY, quantity);
  valVec  = "", "", quantity;
  params->GetList( keyVec, attrVec, valVec, calcEnergy_ );

  // If the the symbolic name is "all" compute energy for all regions
  if ( calcEnergy_.GetSize() == 1 && calcEnergy_[0] == "all" ) {
    calcEnergy_ = subdoms_;
  }

  if ( calcEnergy_.GetSize() > 0 ) {
   	hasOutput_ = TRUE;
    Info->PrintF( pdename_, " Computing energy for regions:\n" );
    for ( Integer k = 0; k < calcEnergy_.GetSize(); k++ ) {
      Info->PrintF( pdename_, " %s", calcEnergy_[k].c_str() );
    }
  }

  // --- Electric Charges ---
  //check for charge computation
  params->GetList( "region", chargeNeighborRegion_, pdename_, "charge" );
  params->GetList( "element", calcCharges_, pdename_, "charge" );

  if (calcCharges_.GetSize() > 0)
    {
   	hasOutput_ = TRUE;
     Info->PrintF( pdename_,
		   " Computing electric charges for regions:\n");
     for ( Integer k = 0; k < calcCharges_.GetSize(); k++ ) {
       Info->PrintF( pdename_, " %s", calcCharges_[k].c_str() );
    } 
    // Resize solution arrays
    charges_.SetNumSolutions(1);
    charges_.SetSolutionType(ELEC_CHARGE);
    charges_.SetNumElems(numElems_);
    charges_.SetNumDofs(1);
    charges_.SetPtrEQNData(eqnData_, ptgrid_, actlevel_);
    charges_.Init(0);
    } 
  
  // *****************************
  // Determine nodal history
  // *****************************
  StdVector<std::string> saveNodeHist;
  keyVec  = pdename_, "storeResults", "nodeHistory", "saveNodes";
  attrVec = "", "", "type";

  // --- Electric Potential ---
  Enum2String(ELEC_POTENTIAL, quantity);
  valVec  = "", "", quantity;
  params->GetList( keyVec, attrVec, valVec, saveNodeHist );
  
  if (saveNodeHist.GetSize() > 0) {
    saveSolHist_ = TRUE;
   	hasOutput_ = TRUE;
    Info->PrintF( pdename_, " Saving ElecPotential for Nodes:\n" );
    for ( Integer k = 0; k < saveNodeHist.GetSize(); k++ ) {
      Info->PrintF( pdename_, " %s\n", saveNodeHist[k].c_str() );
    }
  }
  
  // *****************************
  // Determine element history
  // *****************************
  StdVector<std::string> saveElemHist;
  keyVec  = pdename_, "storeResults", "elemHistory", "saveElems";
  attrVec = "", "", "";
  valVec = "", "", "";
  params->GetList(keyVec, attrVec, valVec, saveElemHist);
  
  if (saveElemHist.GetSize() > 0) {
    std::string errMsg = pdename_;
    errMsg += ": Saving history elements is not implemented yet!\n";
    errMsg += "Meanwhile you can use 'unvtool' to extract element data.";
    Error( errMsg.c_str(), __FILE__, __LINE__);
  }
  
}


} // end of namespace

