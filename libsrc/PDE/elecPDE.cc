#include <fstream>
#include <iostream>
#include <sstream>
#include <math.h>

#include <DataInOut/Unverg/outUnverg.hh>
#include <DataInOut/GMV/outGMV.hh>
#include <Forms/forms_header.hh>
#include <Forms/elecfieldop.hh>
#include <Forms/elecforceop.hh>
#include <Estimator/spaceerror.hh>
#include <DataInOut/WriteInfo.hh> 
#include <Driver/assemble.hh>
#include <General/defs.hh>

#include <Matrix/matrix.hh>
#include <Utils/vector.hh>
#include "elecPDE.hh"
#include <General/defs.hh>
#include "DataInOut/ParamHandling/BaseParamHandler.hh" 
#include <string>
#include "Utils/StdVector.hh"

namespace CoupledField {

  // ***************
  //   Constructor
  // ***************
  ElecPDE::ElecPDE(Grid * aptgrid, BCs *aptbcs, TimeFunc *aptTimeFunc,
		   FileType *aptFileType, WriteResults *aptOut)
    :BasePDE(aptgrid, aptbcs, aptFileType, aptOut, aptTimeFunc) {

    ENTER_FCN( "ElecPDE::ElecPDE" );

    dofspernode_ = 1;

    pdename_          = "electrostatic";
    pdematerialclass_ = "piezo"; 
 
#ifndef XMLPARAMS
    conf->getsubdompde(subdoms_,pdename_);
#else
    params->GetList( "name", subdoms_, pdename_, "region" );
#endif

    ReadBCs(pdename_);

    //check, if problem is axisymmetric
#ifndef XMLPARAMS

    isaxi_ = FALSE;
    std::string subtype;
    conf->ifget("subtype",subtype,pdename_);
    if (subtype == "axi")
      isaxi_ = TRUE;

    //check for electric field:
    conf->ifgetliststr("calc_EField",calcEfield_,pdename_); 

    //check for electric field energy:
    conf->ifgetliststr("calc_Energy",calcEnergy_,pdename_); 
#else
    // Check whether problem has axial symmetry
    if ( params->HasValue( "type", "axi", "geometry" ) ) isaxi_ = TRUE;
#endif

    

    

  
    SolverCFS_ = FALSE;

    // only static analysis are possible =======================
    delete assemble_;
    assemble_ = new StaticAssemble(algsys_, ptgrid_);
    analysistype_ = STATIC;

  // set analysis parameters
  assemble_->SetGeneralParams(pdename_, dofspernode_, numPDENodes_, subdoms_, surfdoms_);
  assemble_->SetGraphType(NODEGRAPH);
  assemble_->SetMesh2PDENode(&mesh2PDENode_);

#ifdef USE_OLAS
    assemble_->SetMatrixEntryType(OLAS::DOUBLE);
    assemble_->SetMatrixStorageType(OLAS::SPARSE_NONSYM);
    //assemble_->SetMatrixStorageType(HYPRE_MATRIX);
#else
    assemble_->SetMatrixType(RSCALAR);
#endif

    assemble_->SetNumDirichlet(GetNumRestraints(actlevel_));
    assemble_->SetPtrBCs(ptBCs_);
    assemble_->SetPtr2Sol(sol_);
    assemble_->SetPtr2TimeFnc(ptTimeFunc_);

    
    ReadMaterialData();
   
    DefineIntegrators(actlevel_);  

#ifndef XMLPARAMS
    ReadSavings();
#else
    ReadStoreResults();
#endif

    Reset();
  }
  


void ElecPDE::DefineIntegrators(const Integer level)
{
  ENTER_FCN( "ElecPDE::DefineIntegerators" );

  for (int actSD = 0; actSD < subdoms_.GetSize(); actSD++)
    {
      //reads eps33 (matrix notation starts with 0)
      Double eps33 = materialData_[actSD].GetPermittivity(2,2);

      BaseForm * lapl = new LaplaceInt(eps33, isaxi_);

      assemble_->AddIntegrator(lapl, subdoms_[actSD], SYSTEM, nonLin_);
    }
}


// ======================================================
// SOLVING SECTION
// ======================================================


  


void ElecPDE:: PreStepStatic(const Integer kstep, const Double asteptime,
			     const Integer level, const Boolean reset)
{
  ENTER_FCN( "ElecPDE::PreStepStatic" );

  if (pdeIsCoupled_ )
      algsys_->InitSol();

  if (geoUpdate_)
    {
      assemble_->SetNonlinGeo();

      algsys_->InitRHS();
      algsys_->InitSol();
      assemble_->InitMatrices();

      assemble_->SetReassemble();   
    }
}




void ElecPDE::StepStaticNonLin(const Integer kstep, const Double aTime,
			       const Integer level, const Boolean reset)
{
  ENTER_FCN( "ElecPDE::StepStaticNonLin" );

  Integer job = 1;
  Double * ptsol;

  assemble_->AssembleMatrices(level);
  assemble_->AssembleSrcRHS(level);
  
  updateBCs_ = 0;

  SetBCs(level,updateBCs_,aTime);

#ifdef USE_OLAS
  algsys_->BuildInDirichlet();
  algsys_->SetupPrecond(job);
#else
  algsys_->CalcPrecond(job);
#endif

  algsys_->Solve();

  ptsol = algsys_->GetSolutionVal();

  // save solution
  sol_->SetDataPointer(ptsol);
  
}

void ElecPDE::PostStepStatic(const Integer kstep, const Double asteptime,
			     const Integer level)
{
  ENTER_FCN( "ElecPDE::PostStepStatic" );

  if (pdeIsCoupled_)
    iterCoupledCounter_++;


#ifdef ADAPTGRID
  if (flags->CalcErrorMap_)
    {
      Double         totalErr;
      ElemStoreSol<Double>  Sol_Mesh;
      Vector<Double> solVec;
      
      ptError_=new SpaceErrorEstimator();

      ptError_->Init(this);
      
      sol_->TransformNodeSolution(Sol_Mesh,sol_,eqnData_,ptgrid_);

      sol_->GetCompleteVector(solVec);

      //  solVec.Resize(sol_.size());
      //       int i;
      //       for (i=0; i<sol_.size(); i++)
      // 	solVec[i]=sol_[0][i];

      ptError_->CalcErrorMap(solVec,subdoms_,ptgrid_,errorMap_,totalErr,level);
      
      std::cout << " total error of calculation:: " << totalErr << std::endl;
      *data << errorMap_ << std::endl;
    }
#endif
}




// ======================================================
// POSTPROCESSING SECTION
// ======================================================

void ElecPDE::WriteResultsInFile()
{
  ENTER_FCN( "ElecPDE::WriteResultsInFile" );

  Double time = lasttimecalc_;
  
#ifdef PARALLEL //only one thread should write output
  int commrank;
  MPI_Comm_rank(MPI_COMM_WORLD,&commrank);
  if (!commrank){
#endif
    
    
  
  // ATTENTION:
  // The errorMap should be assigned as a StoreSolution, not as a 
  // Vector. This is only temporarely
  ElemStoreSol<Double> Error, Error_Mesh;


  NodeStoreSol<Double> const & solConverted = 
    dynamic_cast<NodeStoreSol<Double>&>(*sol_);
  
   
  // write results
  if (outFile_->IsGMV())
    {
      // write electric potential
      if (savesol_)
	outFile_->WriteNodeSolution(solConverted, laststepcalc_, time, "E-Potential");
      
      if (calcEfield_.GetSize() !=0 )
	{
	  //E_.TransformElemSolution(E_Mesh,ptgrid_,pde2MeshElem_,actlevel_);
	  outFile_->WriteElemSolution(E_,laststepcalc_,time,"E-Field");
	}
    }
  else
    {
      // write electric potential
      if (savesol_)
	outFile_->WriteNodeSolution(solConverted,laststepcalc_,time,"electric potential");

      if (calcEfield_.GetSize() !=0 )
	{
	  outFile_->WriteElemSolution(E_, laststepcalc_, time, "electric field");
	}
    }

    if (flags->CalcErrorMap_)
      {
	// this is only a temporar solution
	Error.SetNumSolutions(1);
	Error.SetNumNodes(errorMap_.GetSize());
	Error.SetSolutionType(NO_SOLUTION_TYPE);
	Error.SetNumDofs(dofspernode_);
	Error.Init(0);
	Error.SetCompleteVector(errorMap_);
	
	// ATTENTION!!
	// up to now now transformation of the Error performed,
	// since the calculation of the error is done on the global element numeration
	//Error.TransformElemSolution(Error_Mesh,subdoms_,ptgrid_,actlevel_);
	//OutFile_->WriteElemSolution(errorMap_, laststepcalc_, time, "relERR-E-Potential"); 
	outFile_->WriteElemSolution(Error_Mesh, laststepcalc_, time, "relERR-E-Potential"); 
      }

    if (calcEnergy_.GetSize() !=0 )
      CalcEnergy();


#ifdef PARALLEL
  }//!commrank
#endif

}

void ElecPDE::PostProcess(const Integer level)
{
  ENTER_FCN( "ElecPDE::PostProcess" );

  NodeStoreSol<Double> & solhelp = dynamic_cast<NodeStoreSol<Double>&>(*sol_);

  if (calcEfield_.GetSize() !=0 )
    {
      ElecFieldOp * FieldOp = new ElecFieldOp(ptgrid_, this, eqnData_,solhelp, level, isaxi_);

      // ------ Calculation of the electric field ------

      Vector<Double> LCoord;
      LCoord.Resize(dim_);
      LCoord[0] = 0;
      LCoord[1] = 0;
      
      StdVector<Elem*> elemssd;
      Integer counterElems=0;
      Vector<Double> TempE;
      
      // loop over all subdomains
      for (Integer isd=0; isd<calcEfield_.GetSize(); isd++)
	{
	  // get vector of Elem of subdomain with color: subdoms[isd]
	  ptgrid_->GetElemSD(elemssd,calcEfield_[isd],level);
	  
	  // loop over elements of subdomain
	  for (Integer iel=0; iel< elemssd.GetSize(); iel++,counterElems++)
	    {
	      FieldOp->CalcElemElecField( TempE, elemssd[iel], LCoord); 
	      E_.SetNodalResult(mesh2PDEElem_[elemssd[iel]->elemNum - 1]-1,TempE);
	    }
	}
      delete FieldOp;
    }
}


void ElecPDE::CalcNodeForce(Vector<Double> & force, 
			    StdVector<Integer> & nodes, 
			    StdVector<Elem*> & elems,
			    StdVector<StdVector<ShortInt> > & isBoundaryNode,
			    StdVector<StdVector<Integer> > & elemNodeToCouplingNode)
{
  ENTER_FCN( "ElecPDE::CalcNodeForce" );

  NodeStoreSol<Double> * solhelp = dynamic_cast<NodeStoreSol<Double> *>(sol_);

  
  ElecForceOp * ForceOp = new ElecForceOp(ptgrid_, this, eqnData_, *solhelp, actlevel_, isaxi_);
   
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
	for( Integer idim=0; idim<dim_; idim++)
	  force[elemNodeToCouplingNode[ielem][ielemnode]*dim_+idim] += 
	    force_temp(ielemnode,idim);
    }


  Vector<Double> sum;
  sum.Resize(dim_);
  
  for (Integer i=0; i<nodes.GetSize(); i++)
    for (Integer dim=0; dim<dim_; dim++)
      sum[dim] += force[i*dim_+dim];

  delete ForceOp;


  //std::cerr << "Sum of E-Force:" << std::endl << sum << std::endl;
  
  // write information in .info-file
  Info->PrintF(pdename_, "Sum of electrostatic force (VWM):");
  Info->PrintVec(sum);
}





void ElecPDE::CalcEnergy()
{
  ENTER_FCN( "ElecPDE::CalcEnergy" );

  Matrix<Double> elemmat;  
  Matrix<Double> ptCoord;
  BaseFE         * ptElem;

  StdVector<Integer> connecth, connect_PDE, Eqns;  
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
	  Mesh2PDENode(connect_PDE,connecth,mesh2PDENode_);

// 	  EqnData_->Mesh2Eqn(Eqns,connecth);
// 	  (*debug) << "Nodes:" << connecth << std::endl;
// 	  (*debug) << "Eqns :" << Eqns << std::endl;


	  Vector<Double> elpot;
	  GetSolOfElement(elpot, connect_PDE);	 
	  help = elemmat * elpot;
	  energy[i] += help * elpot;

	  delete bilinear_stiff;	  
	}  
    }

  std::string resulttype = "Electric Energy";
  Info->WriteResult(pdename_,  resulttype, subdoms_ , energy);
}



void ElecPDE::GetSolOfElement( Vector<Double>& elpot, 
			       StdVector<Integer>& connect_PDE)
{
  ENTER_FCN( "ElecPDE::GetSolOfElement" );

  ElemStoreSol<Double> * solhelp = dynamic_cast<ElemStoreSol<Double> *>(sol_);

  elpot.Resize(connect_PDE.GetSize());
  for(Integer actNode=0; actNode<connect_PDE.GetSize(); actNode++)
    sol_->Get(connect_PDE[actNode]-1,0,elpot[actNode]);
    //elpot[actNode] = (*solhelp)(connect_PDE[actNode]-1,0);
}


// ======================================================
// GENERAL  SECTION
// ======================================================


void ElecPDE::Reset()
{
  ENTER_FCN( "ElecPDE::Reset" );
    
  // Map global numeration of element and nodes to local one
  AssignPDENodeNumbers(mesh2PDENode_, pde2MeshNode_, subdoms_);  
  AssignPDEElemNumbers(mesh2PDEElem_, pde2MeshElem_, subdoms_);
  numPDENodes_ = pde2MeshNode_.GetSize();
  numElems_ = pde2MeshElem_.GetSize();

 eqnData_  = new ScalarNodeEQN(ptgrid_, ptBCs_, subdoms_, actlevel_, dofspernode_);
 eqnData_->SetHomoDirichletBCs(bcs_hd_, homDirichDof_);
 eqnData_->CalcMapping();
 //eqnData_->Print(std::cerr);
 assemble_->SetPtr2EQNData(eqnData_); 
 
  // Initalize solution class
  sol_->SetNumSolutions(1);
  sol_->SetSolutionType(ELEC_POTENTIAL);
  sol_->SetNumNodes(numPDENodes_);
  sol_->SetNumDofs(dofspernode_);
  sol_->SetPtrEQNData(eqnData_);
  sol_->Init(0.0);
  
  E_.SetNumSolutions(1);
  E_.SetSolutionType(ELEC_FIELD);
  E_.SetNumNodes(numElems_);
  E_.SetNumDofs(dim_);
  E_.Init(0.0); 
  E_.SetElemMapping(pde2MeshElem_);
  
}


// ======================================================
// COUPLING SECTION
// ======================================================



void ElecPDE::InitCoupling(PDECoupling * Coupling)
{
  ENTER_FCN( "ElecPDE::InitCoupling" );
  
  pdeIsCoupled_ = TRUE;
  ptCoupling_   = Coupling;

  const Integer numCouplings = ptCoupling_->GetNumOutputCouplings();
  

#ifndef XMLPARAMS
  //check, if geometric nonlinearity is switched of by the user
  geoUpdate_ = TRUE;
  nonLin_    = TRUE;    //general nonlinear switch in basepde!
  
  if (conf->get_optionNo("nonlingeo",  pdename_ ))
    {
      geoUpdate_ = FALSE;
      nonLin_    = FALSE;  
    }
#else

    // ==============================================================
    // NOTE: Currently we can only treat geometric non-linearity and
    //       we assume that for a mechanic PDE all regions either
    //       are linear or non-linear!
    // ==============================================================
    StdVector<std::string> nonLinRegion;
    params->GetList( "nonLinear", nonLinRegion, pdename_, "region" );
    // Should not happen with validating parser, but beware!
    if ( nonLinRegion.GetSize() == 0 ) {
      nonLin_ = FALSE;
    }
    else {
      for ( Integer k = 0; k <= nonLinRegion.GetSize(); k++ ) {
	if ( nonLinRegion[k] != nonLinRegion[0] ) {
	  Info->Error( "Non-linearity should be the same for all regions!",
		       __FILE__, __LINE__ );
	}
      }
      nonLin_ = nonLinRegion[0] == "geo" ? TRUE : FALSE;
    }

#endif

  // Initialization of coupling helper arrays
  std::string quantity;
  StdVector<Integer> * couplingnodes;
  StdVector<Elem*> interface_tmp;
  StdVector<StdVector<ShortInt> > isBoundaryNode_tmp;
  //StdVector<Integer> numBoundaryNodes_tmp;
  StdVector<StdVector<Integer> > elemNodeToCouplingNode_tmp;

  F_Interface_.Resize(numCouplings);
  isBoundaryNode_.Resize(numCouplings);
  elemNodeToCouplingNode_.Resize(numCouplings);
  

  for (Integer actCoupling=0; actCoupling<numCouplings; actCoupling++)
    {
      // Initialize arrays for coupling surface elements
      if (ptCoupling_->GetOutputQuantity(actCoupling) == "elecforce" 
	  || ptCoupling_->GetOutputQuantity(actCoupling) == "interfaceForce")
	{
	  ptCoupling_->GetOutputNodes(actCoupling, couplingnodes);
	  if (couplingnodes == 0)
	    std::cerr << "Couplingnodes = 0!!!!" << std::endl;
	  
	  if (ptCoupling_->GetOutputQuantity(actCoupling) == "elecforce")
	    ptgrid_->GetInterfaceNeighbours(*couplingnodes, subdoms_, interface_tmp, actlevel_);
	  else if (ptCoupling_->GetOutputQuantity(actCoupling) == "interfaceForce")
	    {
	      // help construction for correct assignement of predefined values ... :O(
	      StdVector<Elem*>* interface_tmp_Ptr;
	      ptCoupling_->GetOutputNeighbourElems(actCoupling, interface_tmp_Ptr);
	      interface_tmp = *interface_tmp_Ptr;
	    }
	    else 
	      {
		std::string errMsg = "Coupling " + ptCoupling_->GetOutputQuantity(actCoupling) + 
		  " not known! ";	  
		Error(errMsg.c_str(), __FILE__,__LINE__);
	      }
	  
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


  iterCoupledCounter_ = 0;
}
  


void ElecPDE::CalcOutputCoupling()
{
  ENTER_FCN( "ElecPDE::CalcOutputCoupling" );

  std::string quantity;
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
	  
	  if (quantity == "elecpotential")
	    sol_->NodeSolutionToCoupling(*values, *couplingNodes);
	    
	  if (quantity == "elecforce")
	    {
	      CalcNodeForce(*temp, 
			    *couplingNodes, 
			    F_Interface_[forcesCount], 
			    isBoundaryNode_[forcesCount], 
			    elemNodeToCouplingNode_ [forcesCount]);
	      
	      forcesCount++;
	    }


	  if (quantity == "interfaceForce")
	    CalcInterfaceForces(actCoupling);
	    	  
	  break;
	  
	case ELEM:
	  Error("No Element coupling output", __FILE__,__LINE__);
	}
    }
}


Boolean ElecPDE::HasOutput(std::string output)
{
  ENTER_FCN( "ElecPDE::HasOutput" );
  
  if (output == "elecforce")
    return TRUE;

  if (output == "elecpotential")
    return TRUE;

  if (output == "elecfield")
    return TRUE;

  if (output == "interfaceForce")
    return TRUE;
   
  return FALSE;

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
      BaseFE * ptVolElem       = actVolElem->ptElem;
      BaseFE * ptCoupleElem    = actCoupleElem->ptElem;

      const Vector<Double> & intWeights = ptCoupleElem->GetIntWeights();  
      

      StdVector<Integer> connecth = actCoupleElem->connect;
      const Integer spaceDim   = 2;
      
      if (ptCoupleElem->GetDim() == 3)
	Error("CalcInterfaceForces only implemented for 2D!", __FILE__, __LINE__);


      Matrix<Double> ptCoord; 
      StdVector<Integer> connect_PDE;

      GetElemCoords(connecth, ptCoord, actlevel_);
      Mesh2PDENode(connect_PDE, connecth, mesh2PDENode_);



      // get correct permittivity belonging to the neighbouring element of the interface
      Double eps1 = (*innerCouplingMaterials)[actElem]->GetPermittivity(2,2);
      Double eps2 = (*outerCouplingMaterials)[actElem]->GetPermittivity(2,2);


      // the normal vector points outwards of the MECHANICAL domain
      // (see. Kaltenbacher, "Num. Sim. of Mechatr. Act. & Sens." chapter 8.2)
      Vector<Double> n;
      CalcLineNormalVec(n, *actCoupleElem, *(*innerInterfaceVolElems)[actElem]); // points outward own domain



      StdVector<Integer> boundNodesOfVolElem;
  
      Integer ptCount = 0;
      for (Integer nNode=0; nNode < isBoundaryNode_[actCoupling][actElem].GetSize(); nNode++)
	if (isBoundaryNode_[actCoupling][actElem][nNode] == 1)
	  boundNodesOfVolElem.Push_back(nNode);

      if (boundNodesOfVolElem.GetSize() != 2)
	Error("In CalcInterfaceForces boundary nodes != 2!", __FILE__, __LINE__);
      


      // "interfaceForceVec" holds the absolute value of the forces on every node of an interface vector.
      // To establish the final force vectors, every force in every node has to be multiplied by
      // the normal vector of the interface element
      Vector<Double> interfaceForceOnNodes(connect_PDE.GetSize());   // is automatically initialized by 0
      

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

  Info->PrintF(pdename_, "Sum of electrostatic force (Interface):");
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
  ElecFieldOp elecFieldOp(ptgrid_, this, eqnData_, *solTemp, actlevel_,
			  isaxi_);

  elecFieldOp.CalcElemElecField(tempE, actVolElem, lCoord);
}


// ***********************************************************************
//   Obtain information on desired output quantities from parameter file
// ***********************************************************************
#ifdef XMLPARAMS
void ElecPDE::ReadStoreResults() {

  ENTER_FCN( "ElecPDE::ReadStoreResults" );

  // By default we only save the solution at nodal values
  savesol_    = TRUE;
  savederiv1_ = FALSE;
  savederiv2_ = FALSE;

  // Determine regions for which electric field must be computed
  params->CGetList( "region", calcEfield_, "type", "efield", 0, pdename_,
		    "elemResults" );

  // If the the symbolic name is "all" compute electric field for all regions
  if ( calcEfield_.GetSize() == 1 && calcEfield_[0] == "all" ) {
    calcEfield_ = subdoms_;
  }

  if ( calcEfield_.GetSize() > 0 ) {
    Info->PrintF( pdename_, " Computing electric field for regions:" );
    for ( Integer k = 0; k < calcEfield_.GetSize(); k++ ) {
      Info->PrintF( pdename_, " %s", calcEfield_[k].c_str() );
    }
  }

  // Determine regions for which energy must be computed
  params->CGetList( "region", calcEnergy_, "type", "energy", 0, pdename_,
		    "elemResults" );

  // If the the symbolic name is "all" compute energy for all regions
  if ( calcEnergy_.GetSize() == 1 && calcEnergy_[0] == "all" ) {
    calcEnergy_ = subdoms_;
  }

  if ( calcEnergy_.GetSize() > 0 ) {
    Info->PrintF( pdename_, " Computing energy for regions:" );
    for ( Integer k = 0; k < calcEnergy_.GetSize(); k++ ) {
      Info->PrintF( pdename_, " %s", calcEnergy_[k].c_str() );
    }
  }

}
#endif

} // end of namespace

