#include <fstream>
#include <iostream>
#include <string>
#include <math.h>

#include "acousticPDE.hh" 
#include <Estimator/actimeerror.hh>
#include <DataInOut/Unverg/outUnverg.hh>
#include <DataInOut/GMV/outGMV.hh>
#include <Forms/forms_header.hh>
#include <Estimator/spaceerror.hh>
#include "newmark.hh"

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
  NumPDENodes_ = PDE2MeshNode_.size();

  size_ = NumPDENodes_;

  sol_.reshape(dofspernode_, NumPDENodes_);
  sol_.init();

  with_absBCs_=FALSE;
  std::string absBCs="no";
  conf->ifget("absorbingBCs",absBCs,pdename_);
  if (absBCs == "yes")
      with_absBCs_ = TRUE;

  ReadBCs(pdename_);

}


void AcousticPDE::DiscreteParamsPDE()
{
#ifdef TRACE
  (*trace) << "entering AcousticPDE ::DiscreteParamsPDE" << std::endl;
#endif

  MatrixType_ = RSCALAR;
  GraphType_  = NODEGRAPH; 

  SystemMatrix_     = TRUE;
  StiffnessMatrix_  = TRUE;
  MassMatrix_       = TRUE;

  if (with_absBCs_)
    DampingMatrix_  = TRUE;
}


void AcousticPDE::ComputeRHS(const Double atime)
{
#ifdef TRACE
  (*trace) << "entering Acoustic3dPDE::ComputeRHS" << std::endl;
#endif

}

void AcousticPDE::SolveStepTrans(const Integer kstep, const Double asteptime, 
				   const Integer level, const Boolean reset)
{
#ifdef TRACE
  (*trace) << "entering AcousticPDE::SolveStepTrans" << std::endl;
#endif

  lasttimecalc_= asteptime;
  Boolean Recalc=FALSE;

  if (laststepcalc_==kstep && kstep!=0) Recalc=TRUE;
  else laststepcalc_= kstep;

  Double * ptsol;
  Integer update,job;

  //perform predictor step
  TS_alg->Predictor(sol_);

  if (kstep==0)
    {
      update = 0;
      job = 1;
      SetupMatrices(level);
      algsys_->ConstructEffectiveMatrix(matrix_factor_);
      algsys_->InitRHS();
      ComputeRHS(lasttimecalc_);
      TS_alg->UpdateRHS();
    }
  else if (reset)
    {
      update = 1;
      job    = 1;

      algsys_->InitRHS();
      algsys_->InitMatrix(SYSTEM);
      algsys_->ConstructEffectiveMatrix(matrix_factor_);
      ComputeRHS(lasttimecalc_);
      TS_alg->UpdateRHS();
    }
  else
    {
      update = 1;
      job    = 3;
      algsys_->InitRHS();
      ComputeRHS(lasttimecalc_);
      TS_alg->UpdateRHS();
    };

  SetBCs(level,update,lasttimecalc_);
  algsys_->CalcPrecond(job);
  algsys_->Solve();
  ptsol = algsys_->GetSolutionVal();

  // save solution
  Integer k = 0;
  
  for (Integer i=0; i<NumPDENodes_; i++)
    for (Integer dim=0; dim<dofspernode_; dim++)
      sol_[dim][i] = ptsol[k++];

  //perform corrector step  
  TS_alg->Corrector(sol_);


  if (InfoPrint)
    (*infofile) << "maxnode:" <<  ptgrid_->GetMaxnumnodes(level) << std::endl;

}

void AcousticPDE::WriteResultsInFile()
{
#ifdef TRACE
  (*trace) << "entering AcousticPDE::WriteResultsInFile" << std::endl;
#endif

  Array<Double> sol_mesh, solder1_mesh, solder2_mesh;
  Array<Double> solArray_, sol_der1Array_, sol_der2Array_;
  
  solArray_ = sol_;
  sol_der1Array_ = getS1();
  sol_der2Array_ = getS2();

  
  TransformNodeSolution(sol_mesh,sol_,PDE2MeshNode_);
  TransformNodeSolution(solder1_mesh,sol_der1Array_,PDE2MeshNode_);
  TransformNodeSolution(solder2_mesh,sol_der2Array_,PDE2MeshNode_);

  if (OutFile_->IsGMV())
    {
      OutFile_->WriteNodeSolution(sol_mesh,laststepcalc_,lasttimecalc_,"vp");
      //       OutFile_->WriteNodeSolution(solder1_mesh,laststepcalc_,lasttimecalc_,"vp_1der");
      //       OutFile_->WriteNodeSolution(solder2_mesh,laststepcalc_,lasttimecalc_,"vp_2der");
    }
  else
    {
      OutFile_->WriteNodeSolution(sol_mesh,laststepcalc_,lasttimecalc_,"fluid potential");
      OutFile_->WriteNodeSolution(solder1_mesh,laststepcalc_,lasttimecalc_,"fluid potential, 1st deriv., ");
      //      OutFile_->WriteNodeSolution(solder2_mesh,laststepcalc_,lasttimecalc_,"fluid potential, 2nd deriv., ");
    }

}

void AcousticPDE :: InitTimeStepping(const Double dt)
{
#ifdef TRACE
  (*trace) << "entering AcousticPDE::InitTimeStepping" << std::endl;
#endif

  TS_alg = new Newmark(pdename_, algsys_, dofspernode_, NumPDENodes_, DampingMatrix_);
  TS_alg->Init(matrix_factor_, dt);

}


void AcousticPDE::SetupMatrices(const Integer level)
{
#ifdef TRACE
    (*trace) << "entering AcousticPDE::SetupMatrices" << std::endl;
#endif

  Matrix<Double> elemmat;
  Matrix<Double> ptCoord;


  BaseFE * ptEl;
  Double coeffstiff, coeffmass;
  Vector<Integer> connecth, connect_PDE;
  std::vector<Elem*> elemssd;

  Integer i, j;
 
  for (i=0; i<subdoms_.size(); i++)
    {
      //compute material coefficient
      const Double density = materialData_[i].GetDensity();
      const Double compressibility = materialData_[i].GetCompressibility();

      coeffmass  = density*density/compressibility;
      coeffstiff = density;


      ptgrid_->GetElemSD(elemssd,subdoms_[i],level);

      for (j=0; j< elemssd.size(); j++)
	{
	  ptEl = elemssd[j]->ptElem;
    
	  BaseForm * bilinear_mass  = new MassInt(ptEl, coeffmass);
	  BaseForm * bilinear_stiff = new LaplaceInt(ptEl, coeffstiff);

	  connecth=elemssd[j]->connect;
	  GetElemCoords(connecth, ptCoord, level); 

	  // CHANGE connecth
	  Mesh2PDENode(connect_PDE,connecth,Mesh2PDENode_);

	  // stiffness part
	  bilinear_stiff->CalcElementMatrix(ptCoord, elemmat);

#ifdef DEBUG
	  (*debug) << "Connection array  " << std::endl;
	  (*debug)  << connecth  << std::endl;
	  (*debug) << "Stiffnessmatrix, ElementNumber  " <<   i << std::endl;
	  (*debug) << elemmat << std::endl;
#endif     
	  
	  algsys_->SetElementMatrix(elemmat.getinarray(), connect_PDE.get(), connect_PDE.size(), STIFFNESS);

	  // mass part
	  bilinear_mass->CalcElementMatrix(ptCoord, elemmat);

#ifdef DEBUG
	  (*debug) << "Massmatrix, ElementNumber  " << i << std::endl;
	  (*debug) << elemmat << std::endl;
#endif

	  algsys_->SetElementMatrix(elemmat.getinarray(), connect_PDE.get(), connect_PDE.size(), MASS);

	  delete bilinear_stiff;
	  delete bilinear_mass;

	}
    }

  // BEGIN DAMPING MATRIX PART
  if (with_absBCs_) {
    std::vector<Elem*>  DomainBnd;
    Double coeffdamp;
    std::cout<<"Aqui en damping"<<std::endl;
 
    conf->getliststr("bnd_for_absBCs",bnd_absBCs_,"acoustic");
    DomainBnd=ptBCs_->getFacesBC(bnd_absBCs_[0],level);

    for (j=0; j< DomainBnd.size(); j++)
      {
	const Double density = materialData_[0].GetDensity();
	const Double compressibility = materialData_[0].GetCompressibility();
	coeffdamp = density/((sqrt(compressibility/density)));
      
	ptEl=DomainBnd[j]->ptElem;
	// Here MassInt is used to calculate the damping matrix from the surface elements
	BaseForm * linear_damp = new MassInt(ptEl,coeffdamp);

	Integer ii;
	Integer elsize=(DomainBnd[j]->connect).size();
	connecth.Resize(elsize);
	for (ii=0; ii<elsize; ii++)
	  connecth[ii]=(DomainBnd[j]->connect)[ii];

	GetElemCoords(connecth, ptCoord, level);
	// CHANGE connecth
	Mesh2PDENode(connect_PDE,connecth,Mesh2PDENode_);

	linear_damp->CalcElementMatrix(ptCoord, elemmat);

#ifdef DEBUG
	(*debug) << "Connection array  " << std::endl;
	(*debug)  << connecth  << std::endl;

	(*debug) << "Dampingmatrix, ElementNumber  " <<   j << std::endl;
	(*debug) << elemmat << std::endl;
#endif


	algsys_->SetElementMatrix(elemmat.getinarray(), connect_PDE.get(), connect_PDE.size(), DAMPING);

	delete linear_damp;

      }
  }
  // END DAMPING MATRIX PART



#ifdef TRACE
    (*trace) << "Leaving AcousticPDE::SetupMatrices" << std::endl;
#endif
  }


}


