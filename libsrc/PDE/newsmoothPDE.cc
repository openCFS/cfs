#include <fstream>
#include <iostream>
#include <string>
#include <math.h>

#include <DataInOut/Unverg/outUnverg.hh>
#include <DataInOut/GMV/outGMV.hh>
#include <Forms/forms_header.hh>
#include <Estimator/spaceerror.hh>

#include "newsmoothPDE.hh" 

namespace CoupledField
{

SmoothPDE::SmoothPDE(Grid * aptgrid, BCs *aptbcs, TimeFunc *aptTimeFunc, FileType *aptFileType, 
		 WriteResults *aptOut)
:BasePDE(aptgrid, aptbcs, aptFileType, aptOut, aptTimeFunc)
{
#ifdef TRACE
  (*trace) << "entering SmoothPDE::SmoothPDE " << std::endl;
#endif

  firstTurn_ = TRUE;

  pdename_ = "smooth";
  pdematerialclass_ = "piezo";
  
  
  conf->getstr("subtype", subType_, pdename_ );
  
  if (subType_ == "3d")
    dofspernode_ = 3;
  else if (subType_ == "axi")
    {
      isaxi_ = TRUE;
      dofspernode_ = 2;
    }
  else
    dofspernode_ = 2;

  conf->getsubdompde(subdoms_,pdename_);


  //BCs
  ReadBCs(pdename_);
  numDirichletBCs_ = GetNumRestraints(actlevel_);
  
  conf->ifgetliststr("homogenBCDof", homDirichDof_, pdename_);  
  conf->ifgetliststr("inhomogenBCDof", inhomDirichDof_, pdename_);

  // just for consistency with old script
  conf->ifgetliststr("homoBCDof", homDirichDof_, pdename_);
  conf->ifgetliststr("homoBCdof", homDirichDof_, pdename_);
  conf->ifgetliststr("inhomoBCDof", inhomDirichDof_, pdename_);
  conf->ifgetliststr("inhomoBCdof", inhomDirichDof_, pdename_);

  //check for b.c. input data
  if (bcs_hd_.size() != homDirichDof_.size()) 
     Error("Inconsistent definition of homogeneous Dirichlet Boundary Conditions");
  if (bcs_id_.size() != inhomDirichDof_.size()) 
     Error("Inconsistent definition of inhomogeneous Dirichlet Boundary Conditions");

  AssignPDENodeNumbers(Mesh2PDENode_, PDE2MeshNode_, subdoms_);
  numPDENodes_ = PDE2MeshNode_.size();
  numElems_    = ptgrid_->GetMaxnumElem(actlevel_, subdoms_);
  size_ = numPDENodes_*dofspernode_;


  sol_.reshape(dofspernode_, numPDENodes_);
  sol_.init();
  
  method_ = "mechanic";
  conf->ifget("method", method_, pdename_ );

  factor_.Resize(numElems_);
  for (Integer i=0; i<factor_.size(); i++)
    factor_[i] = 1.0;
    
  // set assemble parameters
  assemble_->SetGeneralParams(pdename_, dofspernode_, numPDENodes_, subdoms_, surfdoms_);
  assemble_->SetGraphType(NODEGRAPH);
  assemble_->SetMesh2PDENode(&Mesh2PDENode_);
  assemble_->SetMatrixType(RBLOCK);
  assemble_->SetNumDirichlet(numDirichletBCs_);
  assemble_->SetPtrBCs(ptBCs_);
  assemble_->SetPtr2Sol(&sol_);
  assemble_->SetPtr2TimeFnc(ptTimeFunc_);
  
  ReadMaterialData();
   
  DefineIntegrators(actlevel_);  

  ReadSavings();

  //is a nonlinear PDE, since in each iteration, we have to setup the matrices new!
  nonLin_ = TRUE;
}


void SmoothPDE::DefineIntegrators(const Integer level)
{
#ifdef TRACE
  (*trace) << "entering SmoothPDE::DefineIntegerators" << std::endl;
#endif

  Boolean nonLin = FALSE;

  for (int actSD = 0; actSD < subdoms_.size(); actSD++)
    {
      // ==============  add stiffness ===========================================

      MaterialData actSDMat(materialData_[actSD]);

      // ==============  add "standard" stiffness ===============================
      BaseForm * bilinearStiff;
      if (subType_ == "plainStrain")
	bilinearStiff = new smoothPlainStrainInt(actSDMat);
      else if (subType_ == "3d")
	bilinearStiff = new smooth3DInt(actSDMat);
      else 
	Error("Unknown subtype in smooth PDE! ",__FILE__,__LINE__);

      assemble_->AddIntegrator(bilinearStiff, subdoms_[actSD], STIFFNESS, nonLin);
    }
}


void SmoothPDE::InitCoupling(PDECoupling * coupling)
{
#ifdef TRACE
  (*trace) << "entering SmoothPDE::Initcoupling" << std::endl;
#endif

  PDEisCoupled_ = TRUE;
  ptCoupling_   = coupling; 

  // input couplings
  for (Integer i=0; i<ptCoupling_->GetNumInputCouplings(); i++)
    {
      // check for input mechanic displacement
      if (ptCoupling_->GetInputQuantity(i) == "mechdisplacement")
	numDirichletBCs_ += (dofspernode_ * ptCoupling_->GetInputNumNodes(i));
    }

  // output couplings
  for (Integer i=0; i<ptCoupling_->GetNumOutputCouplings(); i++)
    {
      // check for output displacement
      if (ptCoupling_->GetOutputQuantity(i) == "smoothdisplacement")
	  ptCoupling_->SetOutputDof(i, dofspernode_);
    }

  // now overwrite number of Dirichlet BCs due to coupling 
  assemble_->SetNumDirichlet(numDirichletBCs_);

  iterCoupledCounter_ = 0;
}
  

void SmoothPDE::PreStepStatic(const Integer level)
{
#ifdef TRACE
  (*trace) << "entering SmoothPDE::PreStepStatic" << std::endl;
#endif

  algsys_->InitRHS();
  algsys_->InitSol();
  assemble_->InitMatrices();

  assemble_->SetReassemble();

}


void SmoothPDE::StepStaticNonLin(const Integer level)
{
#ifdef TRACE
  (*trace) << "entering SmoothPDE::StepStaticLin" << std::endl;
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


void SmoothPDE:: PostStepStatic(const Integer level)
{
#ifdef TRACE
  (*trace) << "entering SmoothPDE::PostStepStatic" << std::endl;
#endif
  
  iterCoupledCounter_++;
}


void SmoothPDE::CalcOutputCoupling()
{
#ifdef TRACE
  (*trace) << "entering SmoothPDE::CalcOutputCoupling" << std::endl;
#endif  

  std::string quantity;
  std::vector<Integer> * couplingnodes;
  Array<Double> * values;

  // loop over all output coupling quantities
  for (Integer i=0; i<ptCoupling_->GetNumOutputCouplings(); i++)
    {
      quantity = ptCoupling_->GetOutputQuantity(i);

      switch(ptCoupling_->GetOutputType(i))
	{

	case NODE:
	  
	  ptCoupling_->GetOutputNodes(i, couplingnodes);
	  ptCoupling_->GetOutputValues(i, values);

	  if (quantity == "smoothdisplacement")
	    {
	      NodeSolutionToCoupling(*values, *couplingnodes);
	    }
	  
	  break;

	case ELEM:
	  Error("No Element coupling output", __FILE__,__LINE__);
	}

    }
}

void SmoothPDE::WriteResultsInFile()
{
#ifdef TRACE
  (*trace) << "entering SmoothPDE::WriteResultsInFile" << std::endl;
#endif

  Array<Double> DispMesh;
 
  TransformNodeSolution(DispMesh, sol_, PDE2MeshNode_);
  OutFile_->WriteNodeSolution(DispMesh, laststepcalc_, lasttimecalc_,"displacement");
}



Integer SmoothPDE::GetNrBCDof(const std::string & dofStartString)
  {
#ifdef TRACE
    (*trace) << "entering SmoothPDE::GetNrBCDof" << std::endl;
#endif

   Integer nrActDof;
   
   if (dofStartString == "ux")
      nrActDof = 1;
    else 
      if (dofStartString == "uy")
	nrActDof = 2;
      else
	if (dofspernode_ == 3)
	  if (dofStartString == "uz")
	    nrActDof = 3;
	  else
	    {
	    Error("Unknown dof-type in homog. BC; substring must start with ux, uy or uz!!",__FILE__,__LINE__);
	    std::cerr << dofStartString << std::endl;
	    }
	else
	  {
	    std::cerr << dofStartString << std::endl;
	    Error("Unknown dof-type in homog. BC; substring must start with ux or uy!!",__FILE__,__LINE__);
	  }
    
    return nrActDof;
  }


Integer SmoothPDE::GetBCDof(const std::string dofString)
{
#ifdef TRACE
  (*trace) << "entering SmoothPDE::GetBCDof " << std::endl;
#endif

  if (dofString == "ux")
    return 1;
  else
    if (dofString == "uy")
      return 2;
    else
      if (dofString == "uz")
	return 3;
      else
	Error("The direction mentioned in the config-file is not implemented! ",__FILE__,__LINE__);
}
  

Boolean SmoothPDE::HasOutput(std::string output)
{
#ifdef TRACE
  (*trace) << "entering SmoothPDE::HasOutput" << std::endl;
#endif
  
  if (output == "smoothdisplacement")
    return TRUE;

  return FALSE;
}


} // end namespace CoupledField
