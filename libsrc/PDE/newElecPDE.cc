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

  //check for nonlinearity
  nonLin_ = FALSE;
  if (conf->get_option("nonlin",  pdename_ ))
    nonLin_=TRUE;


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
  assemble_->SetMatrixType(RSCALAR);
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

  Boolean nonLin = FALSE;

  for (int actSD = 0; actSD < subdoms_.size(); actSD++)
    {
      //reads eps33 (matrix notation starts with 0)
      Double eps33 = materialData_[actSD].GetPermittivity(2,2);

      BaseForm * lapl = new LaplaceInt(eps33, isaxi_);

      assemble_->AddIntegrator(lapl, subdoms_[actSD], SYSTEM, nonLin);
      //assemble_->AddIntegrator(lapl, subdoms_[actSD], STIFFNESS, nonLin);
    }
}


// ======================================================
// SOLVING SECTION
// ======================================================


  


void ElecPDE:: PreStepStatic(const Integer level)
{
#ifdef TRACE
  (*trace) << "entering ElecPDE:: PreStepStatic" << std::endl;
#endif

  if (PDEisCoupled_ )
      algsys_->InitSol();

  if (nonLin_)
    {
      assemble_->SetNonlinGeo();

      algsys_->InitRHS();
      algsys_->InitSol();
      assemble_->InitMatrices();

      assemble_->SetReassemble();   
    }
}




void ElecPDE::StepStaticNonLin(const Integer level, const Double aTime)
{
#ifdef TRACE
  (*trace) << "entering ElecPDE::StepStaticNonLin" << std::endl;
#endif

  Integer job = 1;
  Double * ptsol;

  assemble_->AssembleMatrices(level);
  assemble_->AssembleSrcRHS(level);
  
  updateBCs_ = 0;
  SetBCs(level,updateBCs_,0);
  algsys_->CalcPrecond(job);

  algsys_->Solve();

  ptsol = algsys_->GetSolutionVal();

  // save solution
  Integer k=0;
  
  for (Integer i=0; i<numPDENodes_; i++)
    for (Integer dim=0; dim<dofspernode_; dim++)
      sol_[dim][i] = ptsol[k++];

}

void ElecPDE::PostStepStatic(const Integer level)
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

  Integer step=0;
  Double time= lasttimecalc_;
  ShortInt Dim = ptgrid_->GetDim();

  Array<Double> E_Mesh, Force_Mesh, Sol_Mesh;
  
  // transform solution vector for electric potential
  TransformNodeSolution(Sol_Mesh,sol_,PDE2MeshNode_);


  // CHANGE F_Interface_
  // TransformElemSolution(Force_Mesh,Force_,F_Interface_[0]);
   
  // write results
  if (OutFile_->IsGMV())
    {
      // write electric potential
      if (savesol_)
	OutFile_->WriteNodeSolution(Sol_Mesh,step,time,"E-Potential");
      
      if (calcEfield_.size() !=0 )
	{
	  TransformElemSolution(E_Mesh,E_,subdoms_);
	  OutFile_->WriteElemSolution(E_Mesh,step,time,"E-Field");
	  //OutFile_->WriteElemSolution(Force_Mesh,step,time,"E-Force");
	}
    }
  else
    {
      // write electric potential
      if (savesol_)
	OutFile_->WriteNodeSolution(Sol_Mesh,step,time,"electric potential");

      if (calcEfield_.size() !=0 )
	{
	  //	  TransformElemSolution(E_Mesh,E_,subdoms_);
	  //	  OutFile_->WriteElemSolution(E_Mesh,step,time,"electric field");
	  OutFile_->WriteElemSolution(E_,step,time,"electric field");
	  //OutFile_->WriteElemSolution(Force_Mesh,step,time,"mag. volume force");
	}
    }

    if (flags->CalcErrorMap_)
      OutFile_->WriteElemSolution(errorMap_,step,time,"relERR-E-Potential"); 
      

    if (calcEnergy_.size() !=0 )
      CalcEnergy();
}

void ElecPDE::PostProcess(const Integer level)
{
  
#ifdef TRACE
  (*trace) << "entering ElecPDE::PostProcess" << std::endl;
#endif  


  if (calcEfield_.size() !=0 )
    {
      ElecFieldOp * FieldOp = new ElecFieldOp(ptgrid_, this, &Mesh2PDENode_, &sol_[0], level);

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
  
  
  ElecForceOp * ForceOp = new ElecForceOp(ptgrid_, this, &Mesh2PDENode_, &sol_[0], actlevel_);
   
  Array<Double> force_temp;
  force.reshape(Dim_, nodes.size());
  force.init();
  
  for (Integer ielem=0; ielem<elems.size(); ielem++)
    {
      // Get Material Parameter
      Double epsilon;
      
      for (Integer i=0; i<subdoms_.size(); i++)
	{
	  if (elems[i]->namesd == subdoms_[i]) 
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


  // Initialization of coupling helper arrays
  std::string quantity;
  std::vector<Integer> * couplingnodes;
  std::vector<Elem*> interface_tmp;
  std::vector<std::vector<ShortInt> > isBoundaryNode_tmp;
  //std::vector<Integer> numBoundaryNodes_tmp;
  std::vector<std::vector<Integer> > elemNodeToCouplingNode_tmp;
  Array<Double> * values;

  for (Integer i=0; i<ptCoupling_->GetNumOutputCouplings(); i++)
    {
      // Initialize arrays for coupling surface elements
      if (ptCoupling_->GetOutputQuantity(i) == "elecforce")
	{
	  ptCoupling_->GetOutputNodes(i, couplingnodes);
	  if (couplingnodes == 0)
	    std::cerr << "Couplingnodes = 0!!!!" << std::endl;
	  
	  ptgrid_->GetInterfaceNeighbours(*couplingnodes, subdoms_, interface_tmp, actlevel_);
	  
	  F_Interface_.push_back(interface_tmp);

	  // Intialize the memory
	  ptCoupling_->SetOutputDof(i, Dim_);
	 
	  // Assign arrays for element boundary nodes
	  isBoundaryNode_tmp.clear();
	  isBoundaryNode_tmp.resize(interface_tmp.size());
	  elemNodeToCouplingNode_tmp.clear();
	  elemNodeToCouplingNode_tmp.resize(interface_tmp.size());
	  //numBoundaryNodes_tmp.clear();
	  //numBoundaryNodes_tmp.resize(interface_tmp.size());
	 
	  
	  for (Integer ielem=0; ielem<interface_tmp.size(); ielem++)
	    {
	      //std::cerr << "Interface[" << ielem << "]: " << interface_tmp[ielem]->ElemNum << std::endl;

	      isBoundaryNode_tmp[ielem].clear();
	      isBoundaryNode_tmp[ielem].resize(interface_tmp[ielem]->connect.size());
	      elemNodeToCouplingNode_tmp[ielem].clear();
	      elemNodeToCouplingNode_tmp[ielem].resize(interface_tmp[ielem]->connect.size());

	      // Determine Boundary Nodes
	      for (Integer ielemnode=0; ielemnode<isBoundaryNode_tmp[ielem].size(); ielemnode++)
		{
		  for (Integer inodes=0; inodes<(*couplingnodes).size(); inodes++)
		    {
		      if (interface_tmp[ielem]->connect[ielemnode] == (*couplingnodes)[inodes] )
			{
			  isBoundaryNode_tmp[ielem][ielemnode] = 1;
			  elemNodeToCouplingNode_tmp[ielem][ielemnode] = inodes;
			  //numBoundaryNodes_tmp[ielem] +=1;
			  break;
			} // end if
		    } // end for (inodes)
		} // end for (ielemnodes)
	    } // end for (ielems)

	  isBoundaryNode_.push_back(isBoundaryNode_tmp);
	  elemNodeToCouplingNode_.push_back(elemNodeToCouplingNode_tmp);
	  //numBoundaryNodes_.push_back(numBoundaryNodes_tmp);	    
	} // end if
            
    } // end for (i)


  iterCoupledCounter_ = 0;
}
  


void ElecPDE::CalcOutputCoupling()
{
#ifdef TRACE
  (*trace) << "entering ElecPDE::CalcOutputCoupling" << std::endl;
#endif

  std::string quantity;
  std::vector<Integer> * couplingnodes;
  Array<Double> * values;
  Integer forcesCount = 0;
  
  // loop over all output coupling quantities
  for (Integer i=0; i<ptCoupling_->GetNumOutputCouplings(); i++)
    {
      quantity = ptCoupling_->GetOutputQuantity(i);
      
      switch(ptCoupling_->GetOutputType(i))
	{
	  
	case NODE:
	  
	  ptCoupling_->GetOutputNodes(i, couplingnodes);
	  ptCoupling_->GetOutputValues(i, values);
	  
	  if (quantity == "elecpotential")
	    NodeSolutionToCoupling(*values, *couplingnodes);
	    
	  if (quantity == "elecforce")
	    {
	      CalcNodeForce(*values, 
			    *couplingnodes, 
			    F_Interface_[forcesCount], 
			    isBoundaryNode_[forcesCount], 
			    elemNodeToCouplingNode_ [forcesCount]);
	      
	      forcesCount++;
	    }
	  
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

  return FALSE;

}




} // end of namespace


#endif //#ifdef NEWBASEPDE
