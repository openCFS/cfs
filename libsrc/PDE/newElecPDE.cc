#ifdef NEWBASEPDE

#include <fstream>
#include <iostream>
#include <sstream>
#include <math.h>

#include <DataInOut/Unverg/outUnverg.hh>
#include <Forms/forms_header.hh>
#include <Forms/elecfieldop.hh>
#include <Forms/elecforceop.hh>
#include <Estimator/spaceerror.hh>
#include <DataInOut/WriteInfo.hh> 
#include <AlgebraicSystem/LinAlg/linsystem.hh>
#include <Driver/assemble.hh>
#include "ScalarNodeEQN.hh"

#include "newElecPDE.hh"

namespace CoupledField
{


ElecPDE::ElecPDE(Grid * aptgrid, BCs *aptbcs, TimeFunc *aptTimeFunc, FileType *aptFileType, 
		 WriteResults *aptOut)
:BasePDE(aptgrid, aptbcs, aptFileType, aptOut, aptTimeFunc)
{
#ifdef TRACE
  (*trace) << "entering NewElecPDE::ElecPDE " << std::endl;
#endif

  dofspernode_ = 1;  

  pdename_          = "electrostatic";
  pdematerialclass_ = "piezo"; 
 
  conf->getsubdompde(subdoms_,pdename_);
  ReadBCs(pdename_);

 //check, if problem is axisymmetric
  isaxi_ = FALSE;
  std::string subtype;
  conf->ifget("subtype",subtype,pdename_);
  if (subtype == "axi")
    isaxi_ = TRUE;

  //check for electric field:
  conf->ifgetliststr("calc_EField",calcEfield_,pdename_); 

  //check for electric field energy:
  conf->ifgetliststr("calc_Energy",calcEnergy_,pdename_); 

  Reset();
  
  SolverCFS_ = FALSE;


  // only static analysis are possible =======================
  delete assemble_;
  assemble_ = new StaticAssemble(algsys_, ptgrid_);
  analysistype_ = STATIC;


  // set analysis parameters
  assemble_->SetGeneralParams(pdename_, dofspernode_, numPDENodes_, subdoms_, surfdoms_);
  assemble_->SetGraphType(NODEGRAPH);
  assemble_->SetMesh2PDENode(&Mesh2PDENode_);

#ifdef USE_OLAS
  assemble_->SetMatrixEntryType(DOUBLE);
  assemble_->SetMatrixStorageType(SPARSE_NONSYM);
  //assemble_->SetMatrixStorageType(HYPRE_MATRIX);
#else
  assemble_->SetMatrixType(RSCALAR);
#endif


  assemble_->SetNumDirichlet(GetNumRestraints(actlevel_));
  assemble_->SetPtrBCs(ptBCs_);
  assemble_->SetPtr2Sol(&sol_);
  assemble_->SetPtr2TimeFnc(ptTimeFunc_);

  ReadMaterialData();
   
  DefineIntegrators(actlevel_);  

  ReadSavings();
}



void ElecPDE::DefineIntegrators(const Integer level)
{
#ifdef TRACE
  (*trace) << "entering ElecPDE::DefineIntegerators" << std::endl;
#endif

  for (int actSD = 0; actSD < subdoms_.size(); actSD++)
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
#ifdef TRACE
  (*trace) << "entering ElecPDE:: PreStepStatic" << std::endl;
#endif

  if (PDEisCoupled_ )
      algsys_->InitSol();

  if (GeoUpdate_)
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
#ifdef TRACE
  (*trace) << "entering ElecPDE::StepStaticNonLin" << std::endl;
#endif

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
  Integer k=0;

  for (Integer i=0; i<numPDENodes_; i++)
    for (Integer dim=0; dim<dofspernode_; dim++)
      {
	sol_[dim][i] = ptsol[k++];
      }
  

}

void ElecPDE::PostStepStatic(const Integer kstep, const Double asteptime,
			     const Integer level)
{
#ifdef TRACE
  (*trace) << "entering ElecPDE::PostStepStatic" << std::endl;
#endif

  if (PDEisCoupled_)
    iterCoupledCounter_++;


#ifdef ADAPTGRID
  if (flags->CalcErrorMap_)
    {
      Double         totalErr;
      Array<Double>  Sol_Mesh;
      Vector<Double> solVec;
      
      ptError_=new SpaceErrorEstimator();

      ptError_->Init(this);
      
      TransformNodeSolution(Sol_Mesh,sol_,PDE2MeshNode_);

      solVec.Resize(sol_.size());
      int i;
      for (i=0; i<sol_.size(); i++)
	solVec[i]=sol_[0][i];

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
#ifdef TRACE
  (*trace) << "entering ElecPDE::WriteResultsInFile" << std::endl;
#endif

  Double time = lasttimecalc_;

#ifdef PARALLEL //only one thread should write output
  int commrank;
  MPI_Comm_rank(MPI_COMM_WORLD,&commrank);
  if (!commrank){
#endif

  ShortInt Dim = ptgrid_->GetDim();

  Array<Double> E_Mesh, Force_Mesh, Sol_Mesh;
  
  // transform solution vector for electric potential
  TransformNodeSolution(Sol_Mesh,sol_,PDE2MeshNode_);

   
  // write results
  if (OutFile_->IsGMV())
    {
      // write electric potential
      if (savesol_)
	OutFile_->WriteNodeSolution(Sol_Mesh, laststepcalc_, time, "E-Potential");
      
      if (calcEfield_.size() !=0 )
	{
	  TransformElemSolution(E_Mesh,E_,subdoms_);
	  OutFile_->WriteElemSolution(E_Mesh,laststepcalc_,time,"E-Field");
	}
    }
  else
    {
      // write electric potential
      if (savesol_)
	OutFile_->WriteNodeSolution(Sol_Mesh,laststepcalc_,time,"electric potential");

      if (calcEfield_.size() !=0 )
	  //	  TransformElemSolution(E_Mesh,E_,subdoms_);
	  OutFile_->WriteElemSolution(E_, laststepcalc_, time, "electric field");
	
    }

    if (flags->CalcErrorMap_)
      OutFile_->WriteElemSolution(errorMap_, laststepcalc_, time, "relERR-E-Potential"); 
      

    if (calcEnergy_.size() !=0 )
      CalcEnergy();


#ifdef PARALLEL
  }//!commrank
#endif

}

void ElecPDE::PostProcess(const Integer level)
{
  
#ifdef TRACE
  (*trace) << "entering ElecPDE::PostProcess" << std::endl;
#endif  


  if (calcEfield_.size() !=0 )
    {
      ElecFieldOp * FieldOp = new ElecFieldOp(ptgrid_, this, &Mesh2PDENode_, &sol_[0], level, isaxi_);

      // ------ Calculation of the electric field ------

      std::vector<Double> LCoord;
      LCoord.resize(Dim_);
      LCoord[0] = 0;
      LCoord[1] = 0;
      
      std::vector<Elem*> elemssd;
      Integer counterElems=0;
      Vector<Double> TempE;
      
      // Resize solution arrays
      E_.reshape(Dim_, numElems_);
      E_.init();
      
      
      // loop over all subdomains
      for (Integer isd=0; isd<subdoms_.size(); isd++)
	{
	  // get vector of Elem of subdomain with color: subdoms[isd]
	  ptgrid_->GetElemSD(elemssd,subdoms_[isd],level);
	  
	  // loop over elements of subdomain
	  for (Integer iel=0; iel< elemssd.size(); iel++,counterElems++)
	    {
	      FieldOp->CalcElemElecField( TempE, elemssd[iel], LCoord); 
	      E_.setValuesRow(TempE, elemssd[iel]->ElemNum-1);
	    }
	}
      delete FieldOp;
    }
  
}


void ElecPDE::CalcNodeForce(Array<Double> & force, 
			    std::vector<Integer> & nodes, 
			    std::vector<Elem*> & elems,
			    std::vector<std::vector<ShortInt> > & isBoundaryNode,
			    std::vector<std::vector<Integer> > & elemNodeToCouplingNode)
{
#ifdef TRACE
  (*trace) << "entering ElecPDE::CalcNodeForce" << std::endl;
#endif  
  
  
  ElecForceOp * ForceOp = new ElecForceOp(ptgrid_, this, &Mesh2PDENode_, &sol_[0], actlevel_, isaxi_);
   
  Array<Double> force_temp;
  force.reshape(Dim_, nodes.size());
  force.init();
  
  for (Integer ielem=0; ielem<elems.size(); ielem++)
    {
      // Get Material Parameter
      Double epsilon;
      
      for (Integer i=0; i<subdoms_.size(); i++)
	{
	  if (elems[ielem]->namesd == subdoms_[i]) 
	  epsilon = materialData_[i].GetPermittivity(2,2); 
	}
      
      // Only for testing
      //epsilon = 8.854e-12;

      ForceOp->CalcElemElecForce( force_temp, elems[ielem], epsilon, isBoundaryNode[ielem]);
      

      // Add the element force to the according coupling node
      for (Integer ielemnode=0; ielemnode<elems[ielem]->connect.size(); ielemnode++)
	for( Integer idim=0; idim<Dim_; idim++)
	  force[idim][elemNodeToCouplingNode[ielem][ielemnode]] += force_temp[idim][ielemnode];

    }

  Vector<Double> sum;
  sum.Resize(Dim_);
  
  for (Integer i=0; i<nodes.size(); i++) 
    {
      //std::cerr << "Node[" << nodes[i] << "] = ";
    for (Integer dim=0; dim<Dim_; dim++)
      {
	sum[dim] += force[dim][i];
	//std::cerr << force[dim][i] << " , ";
      }
    //std::cerr << std::endl;
    }
  delete ForceOp;


  //std::cerr << "Sum of E-Force:" << std::endl << sum << std::endl;
  
  // write information in .info-file
  //  Info->PrintF(pdename_, "Sum of electrostatic force:");
  //  Info->PrintVec(sum);

}





void ElecPDE::CalcEnergy()
{
#ifdef TRACE
  (*trace) << "entering ElecPDE::CalcEnergy" << std::endl;
#endif

  Matrix<Double> elemmat;  
  Matrix<Double> ptCoord;
  BaseFE         * ptElem;

  Vector<Integer> connecth, connect_PDE, Eqns;  
  Vector<double> help;

  Integer i, j;
  std::vector<Double> energy(subdoms_.size());

  for (i=0; i<subdoms_.size(); i++)
    {
      //reads eps33 (matrix notation starts with 0)
      Double eps33 = materialData_[i].GetPermittivity(2,2);

      std::vector<Elem*> elemssd;
      ptgrid_->GetElemSD(elemssd,subdoms_[i],actlevel_);

      energy[i] = 0;
      for (j=0; j < elemssd.size(); j++)
	{  
	  ptElem=elemssd[j]->ptElem;
	  BaseForm * bilinear_stiff = new LaplaceInt(ptElem, eps33, isaxi_);

	  connecth=elemssd[j]->connect;
	  GetElemCoords(connecth, ptCoord, actlevel_);
	  bilinear_stiff->CalcElementMatrix(ptCoord, elemmat);

	  // Mape Mesh to PDE node numbers
	  Mesh2PDENode(connect_PDE,connecth,Mesh2PDENode_);

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



void ElecPDE::GetSolOfElement( Vector<Double>& elpot, Vector<Integer>& connect_PDE)
{
#ifdef TRACE
    (*trace) << "entering ElecPDE::GetSolOfElement" << std::endl;
#endif

  // displacement of element nodes
  elpot.Resize(connect_PDE.size());

  for(Integer actNode=0; actNode<connect_PDE.size(); actNode++)
    elpot[actNode] = sol_[0][connect_PDE[actNode]-1];
}


// ======================================================
// GENERAL  SECTION
// ======================================================


void ElecPDE::Reset()
{
#ifdef TRACE
  (*trace) << "entering ElecPDE::Reset" << std::endl;
#endif
    
  //just for Testing
  EqnData_ = new ScalarNodeEQN(ptgrid_,ptBCs_,subdoms_, bcs_hd_,actlevel_);
  EqnData_->Print();

  AssignPDENodeNumbers(Mesh2PDENode_,PDE2MeshNode_,subdoms_);
  numPDENodes_=PDE2MeshNode_.size();

  deltCoords_.reshape(Dim_,numPDENodes_);
  sol_.reshape(dofspernode_,numPDENodes_);
  sol_.init();

  numElems_ = ptgrid_->GetMaxnumElem(actlevel_,subdoms_);
}

// ======================================================
// COUPLING SECTION
// ======================================================



void ElecPDE::InitCoupling(PDECoupling * Coupling)
{
#ifdef TRACE
  (*trace) << "entering ElecPDE::InitCoupling" << std::endl;
#endif

  
  PDEisCoupled_ = TRUE;
  ptCoupling_   = Coupling;

  const Integer numCouplings = ptCoupling_->GetNumOutputCouplings();
  

  //check, if geometric nonlinearity is switched of by the user
  GeoUpdate_ = TRUE;
  nonLin_    = TRUE;    //general nonlinear switch in basepde!
  
  if (conf->get_optionNo("nonlingeo",  pdename_ ))
    {
      GeoUpdate_ = FALSE;
      nonLin_    = FALSE;  
    }

  // Initialization of coupling helper arrays
  std::string quantity;
  std::vector<Integer> * couplingnodes;
  std::vector<Elem*> interface_tmp;
  std::vector<std::vector<ShortInt> > isBoundaryNode_tmp;
  //std::vector<Integer> numBoundaryNodes_tmp;
  std::vector<std::vector<Integer> > elemNodeToCouplingNode_tmp;
  Array<Double> * values;

  F_Interface_.resize(numCouplings);
  isBoundaryNode_.resize(numCouplings);
  elemNodeToCouplingNode_.resize(numCouplings);
  

  for (Integer actCoupling=0; actCoupling<numCouplings; actCoupling++)
    {
      // Initialize arrays for coupling surface elements
      if (ptCoupling_->GetOutputQuantity(actCoupling) == "elecforce" 
	  || ptCoupling_->GetOutputQuantity(actCoupling) == "interfaceForce")
	{
	  ptCoupling_->GetOutputNodes(actCoupling, couplingnodes);
	  if (couplingnodes == 0)
	    std::cerr << "Couplingnodes = 0!!!!" << std::endl;
	  
	  ptgrid_->GetInterfaceNeighbours(*couplingnodes, subdoms_, interface_tmp, actlevel_);
	  
	  F_Interface_[actCoupling] = interface_tmp;

	  // Intialize the memory
	  ptCoupling_->SetOutputDof(actCoupling, Dim_);
	 
	  // Assign arrays for element boundary nodes
	  isBoundaryNode_tmp.clear();
	  isBoundaryNode_tmp.resize(interface_tmp.size());
	  elemNodeToCouplingNode_tmp.clear();
	  elemNodeToCouplingNode_tmp.resize(interface_tmp.size());
	 
	  
	  for (Integer ielem=0; ielem<interface_tmp.size(); ielem++)
	    {
	      //std::cerr << "Interface[" << ielem << "]: " << interface_tmp[ielem]->ElemNum << std::endl;

	      isBoundaryNode_tmp[ielem].clear();
	      isBoundaryNode_tmp[ielem].resize(interface_tmp[ielem]->connect.size());
	      elemNodeToCouplingNode_tmp[ielem].clear();
	      elemNodeToCouplingNode_tmp[ielem].resize(interface_tmp[ielem]->connect.size());

	      // Determine Boundary Nodes
	      for (Integer ielemnode=0; ielemnode<isBoundaryNode_tmp[ielem].size(); ielemnode++)
		for (Integer inodes=0; inodes<(*couplingnodes).size(); inodes++)
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
#ifdef TRACE
  (*trace) << "entering ElecPDE::CalcOutputCoupling" << std::endl;
#endif

  std::string quantity;
  std::vector<Integer> * couplingNodes     = NULL;
  Array<Double> * values = NULL;
  Integer forcesCount = 0;

  // loop over all output coupling quantities
  for (Integer actCoupling=0; actCoupling<ptCoupling_->GetNumOutputCouplings(); actCoupling++)
    {
      quantity = ptCoupling_->GetOutputQuantity(actCoupling);
      
      switch(ptCoupling_->GetOutputType(actCoupling))
	{
	  
	case NODE:
	  
	  ptCoupling_->GetOutputNodes(actCoupling, couplingNodes);
	  ptCoupling_->GetOutputValues(actCoupling, values);
	  
	  if (quantity == "elecpotential")
	    NodeSolutionToCoupling(*values, *couplingNodes);

	    
	  if (quantity == "elecforce")
	    {
	      CalcNodeForce(*values, 
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
#ifdef TRACE
  (*trace) << "entering ElecPDE::HasOutput" << std::endl;
#endif
  
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
#ifdef TRACE
  (*trace) << "entering ElecPDE::CalcInterfaceForces" << std::endl;
#endif

  std::vector<Integer> *      couplingNodes     = NULL;
  Array<Double> *             elemCouplingSols  = NULL;
  std::vector<Elem*> *        couplingElems     = NULL;
  std::vector<Elem*> *        interfaceVolElems = NULL;
  std::vector<MaterialData*>* couplingMaterials = NULL;

  ptCoupling_->GetOutputNodes    (actCoupling, couplingNodes);
  ptCoupling_->GetOutputValues   (actCoupling, elemCouplingSols);
  ptCoupling_->GetOutputElements (actCoupling, couplingElems);
  ptCoupling_->GetOutputMaterials(actCoupling, couplingMaterials);
  ptCoupling_->GetOutputNeighbourElems(actCoupling, interfaceVolElems);
 
  Integer couplingDof = ptCoupling_->GetOutputDof(actCoupling);
  

   
  elemCouplingSols->init();

  for (Integer actElem=0; actElem < couplingElems->size(); actElem++)
    {
      Elem * actCoupleElem     = (*couplingElems)[actElem];
      Elem * actVolElem        = (*interfaceVolElems)[actElem];
      BaseFE * ptVolElem       = actVolElem->ptElem;
      BaseFE * ptCoupleElem    = actCoupleElem->ptElem;
      
      Vector<Integer> connecth = actCoupleElem->connect;
      const Integer spaceDim   = 2;
      
      if (ptCoupleElem->GetDim() == 3)
	Error("CalcInterfaceForces only implemented for 2D!", __FILE__, __LINE__);


      Matrix<Double> ptCoord; 
      Vector<Integer> connect_PDE;

      GetElemCoords(connecth, ptCoord, actlevel_);
      Mesh2PDENode(connect_PDE, connecth, Mesh2PDENode_);



      // get correct permittivity belonging to the neighbouring element of the interface
      Double permittivity = (*couplingMaterials)[actElem]->GetPermittivity(2,2);
      mycout << "Permittivity of interface = " << permittivity << myendl;

//       F_Interface_;
//       elemNodeToCouplingNode_; 



      

      std::vector<Integer> boundNodesOfVolElem;
      
      Integer ptCount = 0;
      for (Integer nNode=0; nNode < isBoundaryNode_[actCoupling][actElem].size(); nNode++)
	if (isBoundaryNode_[actCoupling][actElem][nNode] == 1)
	  boundNodesOfVolElem.push_back(nNode);
      
      
      if (boundNodesOfVolElem.size() != 2)
	Error("In CalcInterfaceForces boundary nodes != 2!", __FILE__, __LINE__);


      // the normal vector points outwards of the MECHANICAL domain
      // (see. Kaltenbacher, "Num. Sim. of Mechatr. Act. & Sens." chapter 8.2)
      Vector<Double> n;
      CalcLineNormalVec(n, *actCoupleElem, *(*interfaceVolElems)[actElem]); // points outward own domain


      ElecFieldOp elecFieldOp(ptgrid_, this, &Mesh2PDENode_, &sol_[0], actlevel_, isaxi_);



      for (Integer actIP=1; actIP <= ptCoupleElem->GetNumIntPoints(); actIP++)
	{
	  std::vector<Double> shpFncAtIp;
	  std::vector<Double> coordAtIP;
	  ptCoupleElem->GetShFncAtIp(shpFncAtIp, actIP);
	  coordAtIP = ptCoord * shpFncAtIp;
	  
	  Matrix<Double> boundElLocCoords;
	  ptCoupleElem->GetLocalCornerCoords(boundElLocCoords);  //coords: x:number, y:Dim

	  mycout << "local coords of interface elem" << myendl << boundElLocCoords << myendl;
	  
	  Matrix<Double> volElLocCoords;
	  ptVolElem->GetLocalCornerCoords(volElLocCoords); //coords: x:number, y:Dim
	  

	  std::vector<Double> lCoord(spaceDim); // coords, at which the e-field has to be calculated
	  Vector<Double> tempE;       // value of a field at lCoord

	  // relPosIP = normed distance in the range of 0..1 between IP 1 and coord 1 of 1D element
	  const Double locElemLen = boundElLocCoords[0][1] - boundElLocCoords[0][0];
	  const Double relPosIP = (coordAtIP[0] - boundElLocCoords[0][0]) / locElemLen;
	  
	  
	  // local x-coord of integration point of boundary element
	  Integer coordIndex = 0;  // x-coord	   
	  lCoord[coordIndex] = volElLocCoords [coordIndex][ boundNodesOfVolElem[0] ] + 
	    relPosIP * (volElLocCoords[coordIndex] [boundNodesOfVolElem[0] ]- volElLocCoords[coordIndex] [boundNodesOfVolElem[1] ]);

	  coordIndex = 1;  //y-coord
	  lCoord[coordIndex] = volElLocCoords[coordIndex][boundNodesOfVolElem[0] ] + 
	    relPosIP * (volElLocCoords[coordIndex][boundNodesOfVolElem[0] ] - volElLocCoords[coordIndex][boundNodesOfVolElem[1] ]);

	  mycout << "Pos of local coord: " << lCoord << myendl;
	  
	  elecFieldOp.CalcElemElecField(tempE, actCoupleElem, lCoord);

      


// 	  if (isaxi_)
// 	  {
// 	    std::vector<Double> shpFncAtIp;
// 	    std::vector<Double> coordAtIP;
// 	    ptCoupleElem->GetShFncAtIp(shpFncAtIp, actIP);
// 	    coordAtIP = ptCoord * ShpFncAtIp;
//             jacDet *= 2 * PI * CoordAtIP[0];
// 	  }

	}

      


//       // copy result into final solution vector
//       for (Integer actNode=0; actNode < ptCoord.size_row(); actNode++)
// 	{
// 	  Integer nodePos = 0;      
	  
// 	  while(connecth[actNode] != (*couplingNodes)[nodePos] && 
// 		nodePos < couplingNodes->size()) 
// 	    nodePos++;
	  
// 	  for (Integer actDof=0; actDof < couplingDof ; actDof++)  
// 	    (*elemCouplingSols)[actDof][nodePos] += tempE[actNode] * n[actDof];
// 	}      
    }
}








} // end of namespace


#endif //#ifdef NEWBASEPDE
