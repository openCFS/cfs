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

#include "magneticPDE.hh"

namespace CoupledField
{

MagPDE::MagPDE(Grid * aptgrid, BCs *aptbcs, TimeFunc *aptTimeFunc, FileType *aptFileType, 
		 WriteResults *aptOut)
:BasePDE(aptgrid, aptbcs, aptFileType, aptOut, aptTimeFunc)
{
#ifdef TRACE
  (*trace) << "entering NewMagPDE::MagPDE " << std::endl;
#endif

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

  //check for permanent magnets
  conf->ifgetliststr("list_magnets", magnetsDomain_, pdename_);

  //check for postprocessing
  conf->ifgetliststr("calc_BField",calcBfield_,pdename_); 
  conf->ifgetliststr("calc_Energy",calcEnergy_,pdename_); 

  AssignPDENodeNumbers(Mesh2PDENode_,PDE2MeshNode_,subdoms_);
  numPDENodes_=PDE2MeshNode_.size();

  deltCoords_.reshape(Dim_,numPDENodes_);
  sol_.reshape(dofspernode_,numPDENodes_);
  sol_.init();

  numElems_ = ptgrid_->GetMaxnumElem(actlevel_,subdoms_);

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



void MagPDE::DefineIntegrators(const Integer level)
{
#ifdef TRACE
  (*trace) << "entering MagPDE::DefineIntegerators" << std::endl;
#endif

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


void MagPDE:: PreStepStatic(const Integer level)
{
#ifdef TRACE
  (*trace) << "entering MagPDE:: PreStepStatic" << std::endl;
#endif

  if (PDEisCoupled_)
      algsys_->InitSol();

}


void MagPDE::PostStepStatic(const Integer level)
{
#ifdef TRACE
  (*trace) << "entering MagPDE::PostStepStatic" << std::endl;
#endif

  if (PDEisCoupled_)
    iterCoupledCounter_++;

}


// ======================================================
// POSTPROCESSING SECTION
// ======================================================

void MagPDE::WriteResultsInFile()
{
#ifdef TRACE
  (*trace) << "entering MagPDE::WriteResultsInFile" << std::endl;
#endif

  Integer step=0;
  Double time=0;
  ShortInt Dim = ptgrid_->GetDim();

  Array<Double> B_Mesh, Force_Mesh, Sol_Mesh;
  
  // transform solution vector for electric potential
  TransformNodeSolution(Sol_Mesh,sol_,PDE2MeshNode_);

  // CHANGE F_Interface_
  // TransformElemSolution(Force_Mesh,Force_,F_Interface_[0]);
   
  // write results
  if (OutFile_->IsGMV())
    {
      // write electric potential
      OutFile_->WriteNodeSolution(Sol_Mesh,step,time,"Mag-Potential");
      
      if (calcBfield_.size() !=0 )
	{
	  TransformElemSolution(B_Mesh,B_,subdoms_);
	  OutFile_->WriteElemSolution(B_Mesh,step,time,"B-Field");
	  //OutFile_->WriteElemSolution(Force_Mesh,step,time,"E-Force");
	}
    }
  else
    {
      // write electric potential
      OutFile_->WriteNodeSolution(Sol_Mesh,step,time,"mag. vector potential");

      if (calcBfield_.size() !=0 )
	{
	  //	  TransformElemSolution(E_Mesh,E_,subdoms_);
	  //	  OutFile_->WriteElemSolution(E_Mesh,step,time,"electric field");
	  OutFile_->WriteElemSolution(B_,step,time,"mag. flux density");
	  //OutFile_->WriteElemSolution(Force_Mesh,step,time,"mag. volume force");
	}
    }


    if (calcEnergy_.size() !=0 )
      CalcEnergy();


}

void MagPDE::PostProcess(const Integer level)
{
  
#ifdef TRACE
  (*trace) << "entering MagPDE::PostProcess" << std::endl;
#endif  


  if (calcBfield_.size() !=0 )
    {
      CurlNodeOp * FieldOp = new CurlNodeOp(ptgrid_, this, &Mesh2PDENode_, &sol_[0], level);
      FieldOp->Set2DType(isaxi_);
 
      // ------ Calculation of the electric field ------

      std::vector<Double> LCoord;
      LCoord.resize(Dim_);
      LCoord[0] = 0;
      LCoord[1] = 0;
      
      std::vector<Elem*> elemssd;
      Integer counterElems=0;
      Vector<Double> TempE;
      
      // Resize solution arrays
      B_.reshape(Dim_, numElems_);
      B_.init();
      
      
      // loop over all subdomains
      for (Integer isd=0; isd<subdoms_.size(); isd++)
	{
	  // get vector of Elem of subdomain with color: subdoms[isd]
	  ptgrid_->GetElemSD(elemssd,subdoms_[isd],level);
	  
	  // loop over elements of subdomain
	  for (Integer iel=0; iel< elemssd.size(); iel++,counterElems++)
	    {
	      FieldOp->CalcElemCurlNode( TempE, elemssd[iel], LCoord); 
	      B_.setValuesRow(TempE, elemssd[iel]->ElemNum-1);
	    }
	}
      delete FieldOp;
    }
  
}


void MagPDE::CalcEnergy()
{
#ifdef TRACE
  (*trace) << "entering MagPDE::CalcEnergy" << std::endl;
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



void MagPDE::GetSolOfElement( Vector<Double>& magvecpot, Vector<Integer>& connect_PDE)
{
#ifdef TRACE
    (*trace) << "entering MagPDE::GetSolOfElement" << std::endl;
#endif

  // displacement of element nodes
  magvecpot.Resize(connect_PDE.size());

  for(Integer actNode=0; actNode<connect_PDE.size(); actNode++)
    magvecpot[actNode] = sol_[0][connect_PDE[actNode]-1];
}


// reads all information in the config file concerning coils 
void MagPDE::ReadCoils()
{
#ifdef TRACE
  (*trace) << "entering MagEdgePDE::ReadCoils" << std::endl;
#endif  
    
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


#endif //#ifdef NEWBASEPDE
