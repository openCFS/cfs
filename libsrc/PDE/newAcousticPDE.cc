#include <fstream>
#include <iostream>
#include <string>
#include <math.h>

#include "newAcousticPDE.hh" 
#include <DataInOut/Unverg/outUnverg.hh>
#include <DataInOut/GMV/outGMV.hh>
#include <Forms/forms_header.hh>
#include <Estimator/spaceerror.hh>
#include "newmark.hh"
#include "newmarkdamp.hh"
#include <DataInOut/WriteInfo.hh>


namespace CoupledField
{

AcousticPDE::AcousticPDE(Grid * aptgrid, BCs *aptbcs, TimeFunc *aptTimeFunc, FileType *aptFileType, 
			 WriteResults *aptOut)
:BasePDE(aptgrid,aptbcs,aptFileType,aptOut,aptTimeFunc)
{
#ifdef TRACE
  (*trace) << "entering AcousticPDE::AcousticPDE " << std::endl;
#endif

  dofspernode_=1;

  pdename_    ="acoustic";
  pdematerialclass_ = "fluid";

  laststepcalc_=0;

  conf->getsubdompde(subdoms_,pdename_);

  AssignPDENodeNumbers(Mesh2PDENode_, PDE2MeshNode_, subdoms_);  
  numPDENodes_ = PDE2MeshNode_.size();

  size_ = numPDENodes_;

  sol_.reshape(dofspernode_, numPDENodes_);
  sol_.init();

  with_fracdamping_=FALSE;
  std::string dampstr;
  conf->ifget("damping",dampstr,pdename_);

  if (dampstr == "fractional")
    {
      Info->PrintF(pdename_,"\n Attenuation according to power law, number of memory is %g\n\n", frac_memory_);
       with_fracdamping_ = TRUE;
       conf->get("frac_memory",frac_memory_,pdename_);
       damping_type_ = FRACTIONAL;
    }

  with_absBCs_=FALSE;
  conf->ifgetliststr("bnd_for_absBCs",bnd_absBCs_,pdename_); 
  if (bnd_absBCs_.size() > 0)
    {
      with_absBCs_ = TRUE;
      if ( damping_type_ == NONE) damping_type_ = ABCDamp;
    }

  ReadBCs(pdename_);


  if (analysistype_==HARMONIC)
    {
      conf->get("frequency", freq_, pdename_);
      solIm_.reshape(dofspernode_, numPDENodes_);
      solIm_.init();
    }

  // set analysis parameters
  assemble_->SetGeneralParams(pdename_, dofspernode_, numPDENodes_, subdoms_, bnd_absBCs_);
  assemble_->SetGraphType(NODEGRAPH);
  assemble_->SetMesh2PDENode(&Mesh2PDENode_);
  assemble_->SetMatrixType(RSCALAR);
  assemble_->SetNumDirichlet(GetNumRestraints(actlevel_));
  assemble_->SetPtrBCs(ptBCs_);
  assemble_->SetPtr2Sol(&sol_);
  assemble_->SetPtr2TimeFnc(ptTimeFunc_);

  if (with_absBCs_ || with_fracdamping_)
    assemble_->NeedDampingMatrix();

  ReadMaterialData();
   
  DefineIntegrators(actlevel_);

  ReadSavings();  
}


  void AcousticPDE::DefineIntegrators(const Integer level)
  {
#ifdef TRACE
  (*trace) << "entering AcousticPDE::DefineIntegerators" << std::endl;
#endif

  Boolean nonLin = FALSE;

  for (Integer actSD = 0; actSD < subdoms_.size(); actSD++)
    {
      Double density = materialData_[actSD].GetDensity();
      Double compressibility = materialData_[actSD].GetCompressibility();

      //stiffness integrator
      BaseForm * bilinearStiff = new LaplaceInt(density,isaxi_);	  
      assemble_->AddIntegrator(bilinearStiff, subdoms_[actSD], STIFFNESS, nonLin);

      //mass integrator
      Double coeffmass  = density*density/compressibility;
      BaseForm * bilinear_mass  = new MassInt(coeffmass, dofspernode_, isaxi_);
      assemble_->AddIntegrator(bilinear_mass, subdoms_[actSD], MASS, nonLin);
    }

  //surface-integration
  // BEGIN DAMPING MATRIX PART: Absorbing boundaries
  if (with_absBCs_ && analysistype_!=HARMONIC) 
    {  
      for (Integer actSD = 0; actSD < bnd_absBCs_.size(); actSD++)
	{
	  //currently hard-coded!!
	  Double density = materialData_[0].GetDensity();
	  Double compressibility = materialData_[0].GetCompressibility();
	  Double coeffdamp = density/((sqrt(compressibility/density)));

	  BaseForm * bilinear_damp = new MassInt(coeffdamp,dofspernode_, isaxi_);
	  assemble_->AddSurfIntegrator(bilinear_damp,  bnd_absBCs_[actSD], DAMPING, nonLin);
	}
    }
  }


// ======================================================
// SOLVING SECTION
// ======================================================

void AcousticPDE :: InitTimeStepping(const Double dt)
{
#ifdef TRACE
  (*trace) << "entering AcousticPDE::InitTimeStepping" << std::endl;
#endif

  if (with_fracdamping_)
    {
      //currently the parameter y is taken from the first subdomain
      //=> currently just one subdomain makes sense
      Double y = materialData_[0].GetDampingBeta();
      //      TS_alg_ = new NewmarkDamp(pdename_, algsys_, dofspernode_, numPDENodes_, damping_type_,
      //			      frac_memory_,y);
    }
  else
    TS_alg_ = new Newmark(pdename_, algsys_, dofspernode_, numPDENodes_, damping_type_);

  TS_alg_->Init(matrix_factor_, dt);

}




// ======================================================
// COUPLING SECTION
// ======================================================


void AcousticPDE::InitCoupling(PDECoupling * Coupling)
{
#ifdef TRACE
  (*trace) << "entering AcousticPDE::InitCoupling" << std::endl;
#endif
  
  PDEisCoupled_ = TRUE;
  ptCoupling_   = Coupling;
  
  for (Integer i=0; i<ptCoupling_->GetNumOutputCouplings(); i++)
    {
      if (ptCoupling_->GetOutputQuantity(i) == "acousticforce")
	{
	  // Nothing to do yet ;-)
	}
    }

  iterCoupledCounter_ = 0;
}




void AcousticPDE::CalcOutputCoupling()
{
#ifdef TRACE
  (*trace) << "entering AcousticPDE::CalcOutputCoupling" << std::endl;
#endif  

  Integer dof;
  std::string quantity;
  std::vector<Elem*> * couplingElems = NULL;
  std::vector<Integer> * couplingNodes = NULL;
  std::vector<MaterialData*> * couplingMaterials = NULL;
  Array<Double> * values = NULL;
  
  // loop over all output coupling quantities
  for (Integer i=0; i<ptCoupling_->GetNumOutputCouplings(); i++)
    {
      quantity = ptCoupling_->GetOutputQuantity(i);

      switch(ptCoupling_->GetOutputType(i))
	{
	case NODE:	  
	  if (quantity == "acousticforce")
	    {
	      ptCoupling_->GetOutputElements(i, couplingElems);
	      ptCoupling_->GetOutputNodes(i, couplingNodes);
	      ptCoupling_->GetOutputMaterials(i, couplingMaterials);
	      ptCoupling_->GetOutputValues(i, values);
	      dof = ptCoupling_->GetOutputDof(i);
	    
	      CalcMechCouplingRHS(couplingElems, *couplingNodes, couplingMaterials, *values, dof);
	    }	  
	  break;

	case ELEM:
	  Error("No Element coupling output", __FILE__,__LINE__);
	}
    }
}

void AcousticPDE::CalcMechCouplingRHS(std::vector<Elem*> * couplingElems, 
				      std::vector<Integer> & couplingNodes,
				      std::vector<MaterialData*> * couplingMaterials,
				      Array<Double>& elemCouplingSols,
				      Integer couplingdof)
{
#ifdef TRACE
  (*trace) << "entering AcousticPDE::CalcMechCouplingRHS" << std::endl;
#endif

   Double density;
   
  for (Integer actElem=0; actElem<couplingElems->size(); actElem++)
    {
      BaseFE * ptElem = (*couplingElems)[actElem]->ptElem;
      Vector<Integer> connecth = (*couplingElems)[actElem]->connect;
      
      Matrix<Double> ptCoord; 
      GetElemCoords(connecth, ptCoord, actlevel_);

      density = (*couplingMaterials)[actElem]->GetDensity();
          
      BaseForm * bilinear_mass = new MassInt(ptElem, density, isaxi_);
      Matrix<Double> elemmat;
      bilinear_mass->CalcElementMatrix(ptCoord, elemmat);
      delete bilinear_mass;	  


      Vector<Integer> connect_PDE;
      Mesh2PDENode(connect_PDE, connecth, Mesh2PDENode_);
      
      Vector<Double> sol;
      GetSolVecOfElement(sol, connect_PDE);

      Vector<Double> forceOnElem = elemmat * sol;

      Vector<Double> n;
      CalcLineNormalVec(n, ptCoord);

     Integer nodePos = 0;
     
     for (Integer actNode=0; actNode<ptCoord.size_row(); actNode++)
       {
	 nodePos = 0;
	 
	 while(connecth[actNode] != couplingNodes[nodePos] && nodePos < couplingNodes.size()) 
	   nodePos++;
	 
	 for (Integer actDof=0; actDof < couplingdof ; actDof++)  
	   elemCouplingSols[actDof][nodePos] += forceOnElem[actNode] * n[actDof];
       }      
    }  
  
}




Boolean AcousticPDE::HasOutput(std::string output)
{
#ifdef TRACE
  (*trace) << "entering AcousticPDE::HasOutput" << std::endl;
#endif

  if (output == "acousticforce")
    return TRUE;

  return FALSE;
}



// ======================================================
// POSTPROCESSING SECTION
// ======================================================

void AcousticPDE::WriteResultsInFile()
{
#ifdef TRACE
  (*trace) << "entering AcousticPDE::WriteResultsInFile" << std::endl;
#endif

  Array<Double> sol_mesh, solder1_mesh, solder2_mesh, solIm_mesh;
  Array<Double> sol_der1Array, sol_der2Array;

  if (analysistype_==HARMONIC)
    {
      TransformNodeSolution(sol_mesh,sol_,PDE2MeshNode_);
      TransformNodeSolution(solIm_mesh,solIm_,PDE2MeshNode_);      
      OutFile_->WriteNodeSolution(sol_mesh,laststepcalc_,lasttimecalc_,"fluid potential, cw realpart,");
      OutFile_->WriteNodeSolution(solIm_mesh,laststepcalc_,lasttimecalc_,"fluid potential, cw imagpart, ");
    }
  else
    {  
      sol_der1Array = getS1();
      sol_der2Array = getS2();

      if (savesol_)
	TransformNodeSolution(sol_mesh,sol_,PDE2MeshNode_);

      if (savederiv1_)    
	TransformNodeSolution(solder1_mesh,sol_der1Array,PDE2MeshNode_);

      if (savederiv2_)
	TransformNodeSolution(solder2_mesh,sol_der2Array,PDE2MeshNode_);
      
      if (OutFile_->IsGMV())
	{
	  if (savesol_)
	    OutFile_->WriteNodeSolution(sol_mesh,laststepcalc_,lasttimecalc_,"vp");
	  if (savederiv1_)
	    OutFile_->WriteNodeSolution(solder1_mesh,laststepcalc_,lasttimecalc_,"vp_1der");
	  if (savederiv2_)
	    OutFile_->WriteNodeSolution(solder2_mesh,laststepcalc_,lasttimecalc_,"vp_2der");
	}
      else
	{
	  if (savesol_)
	    OutFile_->WriteNodeSolution(sol_mesh,laststepcalc_,lasttimecalc_,"fluid potential");

	  if (savederiv1_)
	    OutFile_->WriteNodeSolution(solder1_mesh,laststepcalc_,lasttimecalc_,"fluid potential, 1st deriv.");
	  if (savederiv2_)
	    OutFile_->WriteNodeSolution(solder2_mesh,laststepcalc_,lasttimecalc_,"fluid potential, 2nd deriv.");
	}
    }
}




}


