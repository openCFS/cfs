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
#include <Driver/assemble.hh>
#include "ScalarNodeEQN.hh"
#include "trapezoidal.hh"

#include "magneticPDE.hh"

namespace CoupledField
{

MagPDE::MagPDE(Grid * aptgrid, BCs *aptbcs, TimeFunc *aptTimeFunc, FileType *aptFileType, 
		 WriteResults *aptOut)
:BasePDE(aptgrid, aptbcs, aptFileType, aptOut, aptTimeFunc)
{
  ENTER_FCN( "MagPDE::MagPDE" );

  dofspernode_ = 1;
  
  pdename_          = "magnetic";
  pdematerialclass_ = "magnetic";
  
  conf->getsubdompde(subdoms_,pdename_);
  ReadBCs(pdename_);

 //check, if problem is axisymmetric
  isaxi_ = FALSE;
  std::string subtype;
  conf->ifget("subtype",subtype,pdename_);
  if (subtype == "axi")
    isaxi_ = TRUE;

  //check for coils
  ReadCoils();

  conf->ifget("calc_UI",UIfilename_,pdename_);
  if (UIfilename_.size() > 0)
    {    
      std::string ui_info = "File name for saving voltages/currents is: " + UIfilename_ + "\n";
      Info->PrintF(pdename_,"%s",ui_info.c_str());
      UIfile_ = new std::ofstream(UIfilename_.c_str());

    if (!UIfile_) 
      Error("Can't open UI-file");

    }
  

  //check for permanent magnets
  conf->ifgetliststr("list_magnets", magnetsDomain_, pdename_);

  //check for postprocessing
  conf->ifgetliststr("calc_BField",calcBfield_,pdename_); 
  conf->ifgetliststr("calc_Energy",calcEnergy_,pdename_); 
  conf->ifgetliststr("calc_Eddy",calcEddy_,pdename_); 


  // Map global numeration of element and nodes to local one
  AssignPDENodeNumbers(mesh2PDENode_, pde2MeshNode_, subdoms_);  
  AssignPDEElemNumbers(mesh2PDEElem_, pde2MeshElem_, subdoms_);
  numPDENodes_ = pde2MeshNode_.size();
  numElems_ = pde2MeshElem_.size();

  deltCoords_.Resize(dim_,numPDENodes_);

  // Initalize solution class
  sol_->SetNumSolutions(1);
  sol_->SetSolutionType(MAG_POTENTIAL);
  sol_->SetNumNodes(numPDENodes_);
  sol_->SetNumDofs(dofspernode_);
  sol_->Init(0.0);
  
  


  numElems_ = ptgrid_->GetMaxnumElem(actlevel_,subdoms_);

  // set analysis parameters
  assemble_->SetGeneralParams(pdename_, dofspernode_, numPDENodes_, subdoms_, surfdoms_);
  assemble_->SetGraphType(NODEGRAPH);
  assemble_->SetMesh2PDENode(&mesh2PDENode_);

#ifdef USE_OLAS
  assemble_->SetMatrixEntryType(DOUBLE);
  assemble_->SetMatrixStorageType(SPARSE_NONSYM);
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
}



void MagPDE::DefineIntegrators(const Integer level)
{
  ENTER_FCN( "MagPDE::DefineIntegerators" );

  Boolean nonLin = FALSE;
  Boolean iscoil;
  Integer idxcoil;

  for (int actSD = 0; actSD < subdoms_.size(); actSD++)
    {
      Double reluctivity  = materialData_[actSD].GetPermiability();
      if ( reluctivity == 0)
	Error("Permeability can not be zero!");

      reluctivity = 1/reluctivity;
      BaseForm * curlcurl2D = new CurlCurlNode2DInt(reluctivity, isaxi_);
      assemble_->AddIntegrator(curlcurl2D, subdoms_[actSD], STIFFNESS, nonLin);

      Double conductivity = materialData_[actSD].GetConductivity();      
      BaseForm * bilinear_mass  = new MassInt(conductivity, dofspernode_, isaxi_);
      assemble_->AddIntegrator(bilinear_mass, subdoms_[actSD], MASS, nonLin);

      iscoil = FALSE;
      for (Integer dom=0; dom<coilDef_.size(); dom++)
	if (subdoms_[actSD] == coilDomain_[dom]) 
	  {
	    iscoil = TRUE;
	    idxcoil = dom;
	  }

      if (iscoil)
	{
	  Double factor = coilDef_[idxcoil].current / coilDef_[idxcoil].coilArea;
	  BaseForm * coilSource = new VolumeSrcInt(factor,isaxi_);
	  assemble_->AddRhsSrcIntegrator(coilSource, subdoms_[actSD],coilDef_[idxcoil].timefnc, nonLin_);
	}
	 
    }
}


// ======================================================
// SOLVING SECTION
// ======================================================


void MagPDE :: InitTimeStepping(const Double dt)
{
  ENTER_FCN( "MagPDE::InitTimeStepping" );

  TS_alg_ = new Trapezoidal(pdename_, algsys_, 1, numPDENodes_*dofspernode_);
  TS_alg_->Init(matrix_factor_, dt);

}


void MagPDE:: PreStepStatic(const Integer level)
{
  ENTER_FCN( "MagPDE::PreStepStatic" );

  if (pdeIsCoupled_)
      algsys_->InitSol();

}


void MagPDE::PostStepStatic(const Integer level)
{
  ENTER_FCN( "MagPDE::PostStepStatic" );

  if (pdeIsCoupled_)
    iterCoupledCounter_++;

}


// ======================================================
// POSTPROCESSING SECTION
// ======================================================

void MagPDE::WriteResultsInFile()
{
  ENTER_FCN( "MagPDE::WriteResultsInFile" );

  ShortInt Dim = ptgrid_->GetDim();

  StoreSol<Double> B_Mesh, Jeddy_Mesh, Force_Mesh, Sol_Mesh;

  // transform solution vector for electric potential
  sol_->TransformNodeSolution(Sol_Mesh,pde2MeshNode_,ptgrid_,actlevel_);

  // CHANGE F_Interface_
  // TransformElemSolution(Force_Mesh,Force_,F_Interface_[0]);
   
  // write results
  if (outFile_->IsGMV())
    {
      // write magnetic potential
      outFile_->WriteNodeSolution(Sol_Mesh,laststepcalc_,lasttimecalc_,"Mag-Potential");
      
      if (calcBfield_.size() !=0 )
	{
	  B_.TransformElemSolution(B_Mesh,pde2MeshElem_,ptgrid_,actlevel_);
	  outFile_->WriteElemSolution(B_Mesh,laststepcalc_,lasttimecalc_,"B-Field");
	  //outFile_->WriteElemSolution(Force_Mesh,step,time,"E-Force");
	}

    }
  else
    {
      // write magnetic potential
      outFile_->WriteNodeSolution(Sol_Mesh,laststepcalc_,lasttimecalc_,"mag. vector potential");

      if (calcBfield_.size() !=0 )
	{
	  B_.TransformElemSolution(B_Mesh,pde2MeshElem_,ptgrid_,actlevel_);
	  outFile_->WriteElemSolution(B_Mesh,laststepcalc_,lasttimecalc_,"mag. flux density");
	}

      if (calcEddy_.size() !=0 )
	{
	  Jeddy_.TransformElemSolution(Jeddy_Mesh,pde2MeshElem_,ptgrid_,actlevel_);
	  outFile_->WriteElemSolution(Jeddy_Mesh,laststepcalc_,lasttimecalc_,"eddy current");
	}
    }


}

void MagPDE::PostProcess(const Integer level)
{
  ENTER_FCN( "MagPDE::PostProcess" );

  TRY_CAST
  PTRCAST(sol_,StoreSol<Double>,solhelp);
  CATCH_CAST

  if (calcEnergy_.size() !=0 )
    CalcEnergy();

  if (calcBfield_.size() !=0 )
    {
      CurlNodeOp * FieldOp = new CurlNodeOp(ptgrid_, this, &mesh2PDENode_, *solhelp, level);
      FieldOp->Set2DType(isaxi_);
 
      // ------ Calculation of the electric field ------

      std::vector<Double> LCoord;
      LCoord.resize(dim_);
      LCoord[0] = 0;
      LCoord[1] = 0;
      
      std::vector<Elem*> elemssd;
      Integer counterElems=0;
      Vector<Double> TempE;
      
      // Resize solution arrays
      B_.SetNumSolutions(1);
      B_.SetSolutionType(MAG_FIELD);
      B_.SetNumNodes(numElems_);
      B_.SetNumDofs(dim_);
      B_.Init(0);
      
      // loop over all subdomains
      for (Integer isd=0; isd<subdoms_.size(); isd++)
	{
	  // get vector of Elem of subdomain with color: subdoms[isd]
	  ptgrid_->GetElemSD(elemssd,subdoms_[isd],level);
	  
	  // loop over elements of subdomain
	  for (Integer iel=0; iel< elemssd.size(); iel++,counterElems++)
	    {
	      FieldOp->CalcElemCurlNode( TempE, elemssd[iel], LCoord); 
	      B_.SetNodalResult(mesh2PDEElem_[elemssd[iel]->elemNum - 1]-1, TempE);
	    }
	}
      delete FieldOp;
    }

  if (calcEddy_.size() !=0 )
    {
      std::vector<Double> LCoord;
      LCoord.resize(dim_);
      LCoord[0] = 0;
      LCoord[1] = 0;
      
      std::vector<Elem*> elemssd;
      std::vector<Double> ShpFnc, tmp;
      Vector<Double> magVecDeriv1Elem;
      Vector<Integer> connect, connect_PDE;

      Integer counterElems=0;

      // dimension hard coded for .unv file!
      Vector<Double> JeddyElem(3);
      JeddyElem.Init();
      
     
      // Resize solution arrays
      Jeddy_.SetNumSolutions(1);
      Jeddy_.SetSolutionType(MAG_EDDY_CURRENT);
      Jeddy_.SetNumNodes(numElems_);
      // dimension hard coded for .unv file!
      Jeddy_.SetNumDofs(3);
      Jeddy_.Init(0);

      // loop over all subdomains
      for (Integer actSD=0; actSD<subdoms_.size(); actSD++)
	{
	  // get vector of Elem of subdomain with color: subdoms[isd]
	  ptgrid_->GetElemSD(elemssd,subdoms_[actSD],level);
	  Double conductivity = materialData_[actSD].GetConductivity(); 	  

	  // loop over elements of subdomain
	  for (Integer actEl=0; actEl< elemssd.size(); actEl++,counterElems++)
	    {
	      BaseFE * ptEl = elemssd[actEl]->ptElem;
	      ptEl->GetShFnc(ShpFnc,LCoord);

	      connect = elemssd[actEl]->connect;
	      // Mape Mesh to PDE node numbers
	      Mesh2PDENode(connect_PDE,connect,mesh2PDENode_);

	      GetSolDerivOfElement(magVecDeriv1Elem,connect_PDE);
	      magVecDeriv1Elem.ToStdVector(tmp);
	      JeddyElem[0] = tmp * ShpFnc;
	      JeddyElem[0] *= -conductivity;
	      Jeddy_.SetNodalResult(mesh2PDEElem_[elemssd[actEl]->elemNum - 1]-1,JeddyElem);
	    }
	}
    }

  if (UIfilename_.size() > 0 && analysistype_ == TRANSIENT)
    {
      Vector<Double> uiSD;
      ComputeUI(uiSD);
      WriteUI2File(uiSD);
    }
  
}


void MagPDE::CalcEnergy()
{
  ENTER_FCN( "MagPDE::CalcEnergy" );

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
	  Mesh2PDENode(connect_PDE,connecth,mesh2PDENode_);

// 	  EqnData_->Mesh2Eqn(Eqns,connecth);
// 	  (*debug) << "Nodes:" << connecth << std::endl;
// 	  (*debug) << "Eqns :" << Eqns << std::endl;


	  Vector<Double> magvecpot;
	  GetSolOfElement(magvecpot, connect_PDE);	 
	  help = elemmat * magvecpot;
	  energy[i] += help * magvecpot;

	  delete bilinear_stiff;	  
	}  
    }

  std::string resulttype = "Electric Energy";
  Info->WriteResult(pdename_,  resulttype, subdoms_ , energy);
}


void MagPDE::ComputeUI(Vector<Double>& uiSD)
{
  ENTER_FCN( "MagPDE::ComputeUI" );

  uiSD.Resize(coilDef_.size());
  
  // loop over all subdomains
  for (Integer actSD=0; actSD<subdoms_.size(); actSD++)
    {
      for (Integer dom=0; dom<coilDef_.size(); dom++)
	if (subdoms_[actSD] == coilDomain_[dom]) 
	  {
	   
	    std::vector<Elem*> elemssd;		
	    // get vector of Elem of subdomain with color: subdoms[isd]
	    ptgrid_->GetElemSD(elemssd,subdoms_[actSD],actlevel_);
	    

	    // loop over elements of subdomain	    
	    for (Integer actEl=0; actEl< elemssd.size(); actEl++)
	      {
		BaseFE * ptEl = elemssd[actEl]->ptElem;
		
		const Integer nrIntPts= ptEl->GetNumIntPoints();
		const Integer nrNodes = ptEl->GetNumNodes();
		const std::vector<Double> & intWeights = ptEl->GetIntWeights();  
		Double jacDet;
		
		
		Vector<Integer> connect, connect_PDE;
		connect = elemssd[actEl]->connect;

		Matrix<Double> ptCoord;
		GetElemCoords(connect, ptCoord, actlevel_);

		// Mape Mesh to PDE node numbers
		Mesh2PDENode(connect_PDE,connect,mesh2PDENode_);
		
		
		Vector<Double> magVecDeriv1Elem;
		GetSolDerivOfElement(magVecDeriv1Elem,connect_PDE);
		
		Double uiElem=0;
		
		for (Integer actIntPt=1; actIntPt<=nrIntPts;  actIntPt++)
		  {
		    std::vector<Double> shapeFnc;
		    jacDet = ptEl->CalcJacobianDetAtIp(actIntPt, ptCoord);	
		    ptEl -> GetShFncAtIp(shapeFnc, actIntPt);

		    uiElem += shapeFnc * magVecDeriv1Elem;
		    
		    if (isaxi_)
		      {
			std::vector<Double> coordAtIP = ptCoord * shapeFnc;
			uiElem += shapeFnc * magVecDeriv1Elem * 2 * PI * coordAtIP[0] * intWeights[actIntPt-1];
		      }
		    else
		      uiElem += shapeFnc * magVecDeriv1Elem * intWeights[actIntPt-1];
		  }

		uiSD[dom] += uiElem;
	      }	    
	  }
    }
}



void MagPDE::WriteUI2File(Vector<Double>& uiSD)
{
  ENTER_FCN( "MagPDE::WriteUI2File" );

  Vector<Integer> coilIDs;   // just positive ids
  coilIDs.Push_back(abs(coilDef_[0].ID));

  Integer maxID = coilDef_[0].ID;

  for (Integer dom=1; dom < coilDef_.size(); dom++)
    {
      Boolean isInVec = FALSE;
      
      for (Integer dom2=0; dom2 < coilIDs.GetSize(); dom2++)	
	if (abs(coilDef_[dom].ID) == coilIDs[dom2])
	  isInVec = TRUE;
      
      if (!isInVec)
	coilIDs.Push_back(abs(coilDef_[dom].ID));
      
      if (abs(coilDef_[dom].ID) > maxID)
	maxID = coilDef_[dom].ID;
    }

  Vector<Double> uiID(maxID);
  uiID.Init();


  for (Integer dom=0; dom < coilDef_.size(); dom++)
    {
      Integer actCoilID = coilDef_[dom].ID;
      uiID[abs(actCoilID)-1] += uiSD[dom] * actCoilID/abs(actCoilID);
    }


  *UIfile_ << lasttimecalc_ << " \t";
  for (Integer actID=0; actID < coilIDs.GetSize(); actID++)
    {
      if ( coilDef_[coilIDs[actID]-1].type == MEASUREMENT)   
	*UIfile_ << uiID[coilIDs[actID]-1] *  coilDef_[coilIDs[actID]-1].coilArea << " \t";
    }
  
  *UIfile_ << myEndl;
  
  
}



void MagPDE::GetSolOfElement( Vector<Double>& magvecpot, Vector<Integer>& connect_PDE)
{
  ENTER_FCN( "GetSolOfElement" );
  TRY_CAST
  PTRCAST(sol_,StoreSol<Double>,solhelp)
  CATCH_CAST

  // displacement of element nodes
  magvecpot.Resize(connect_PDE.GetSize());

  for(Integer actNode=0; actNode<connect_PDE.GetSize(); actNode++)
    magvecpot[actNode] = (*solhelp)(connect_PDE[actNode]-1,0);
}

void MagPDE::GetSolDerivOfElement( Vector<Double>& magvecpot_deriv1, Vector<Integer>& connect_PDE)
{
  ENTER_FCN( "MagPDE::GetSolDerivOfElement" );
  
  // displacement of element nodes
  magvecpot_deriv1.Resize(connect_PDE.GetSize());
  const Vector<Double>  & sol_der1 = getS1();

  for(Integer actNode=0; actNode<connect_PDE.GetSize(); actNode++)
    magvecpot_deriv1[actNode] = sol_der1[connect_PDE[actNode]-1];
}


// reads all information in the config file concerning coils 
void MagPDE::ReadCoils()
{
  ENTER_FCN( "MagEdgePDE::ReadCoils" );
    
  conf->ifgetliststr("list_coils", coilDomain_, pdename_);

  Integer nrCoils = coilDomain_.size();

    if (nrCoils)
      {
	coilDef_.resize(nrCoils);
	for (Integer i=0; i < nrCoils; i++)
	  {
	    conf->getCoilData(coilDef_[i],pdename_,coilDomain_[i]+":");
	    Info->PrintCoil(coilDomain_[i], coilDef_[i], analysistype_);
	  }
      }  
  }


// ======================================================
// COUPLING SECTION
// ======================================================





} // end of namespace

