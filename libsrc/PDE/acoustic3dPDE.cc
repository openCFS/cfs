#include <fstream>
#include <iostream>
#include <string>
#include <math.h>

#include "acoustic3dPDE.hh"
#include <Estimator/actimeerror.hh>
#include <DataInOut/Unverg/outUnverg.hh>
#include <DataInOut/GMV/outGMV.hh>
#include <Forms/forms_header.hh>
#include <Estimator/spaceerror.hh>


namespace CoupledField
{

Acoustic3dPDE::Acoustic3dPDE(AbstractAlgebraicSys * ptalgsys, Grid * aptgrid, Material * ptMaterial, TimeFunc * aptTimeFunc, FileType * aptFileType, WriteResults * aptOut)
:BasePDE(ptalgsys,ptMaterial,aptFileType,aptOut,aptTimeFunc)
{
#ifdef TRACE
  (*trace) << "entering Acoustic3dPDE::Acoustic3dPDE " << std::endl;
#endif
  if (!MatFile_) Error("You didn't specified material file. Check your config-file.");

  dofspernode_=1;
  ptgrid_=aptgrid;

  laststepcalc_=0;
  size_=ptgrid_->GetMaxnumnodes(0);

  sol_.Resize(size_);
  sol_.Init(0);

  sol_der1_.Resize(size_);
  sol_der1_.Init(0);

  sol_der1_old_.Resize(size_);
  sol_der1_old_.Init(0);

  sol_der1_.Resize(size_);
  sol_der1_.Init(0);

  sol_der2_.Resize(size_);
  sol_der2_.Init(0);

  sol_old_.Resize(size_);
  sol_old_.Init(0);

  sol_der2_old_.Resize(size_);
  sol_der2_old_.Init(0);

  s_oldold_.Resize(size_);
  s_oldold_.Init(0);

  conf->get("alpha_NM",alpha_,"Acoustic3d");
  conf->get("beta_NM",beta_,"Acoustic3d");
  conf->get("gamma_NM",gamma_,"Acoustic3d");

  conf->getsubdompde(subdoms_,"Acoustic3d");
  ReadBCs("Acoustic3d");

  //  without_absBCs_=conf->get_option("without_absorbingBCs");
  without_absBCs_= TRUE;

  preComputeRHS();
}

void Acoustic3dPDE::SpecifySolver(Integer &solvertype, Integer &precondtype, Double &eps, Double &dampiter, Integer &maxnumit, Integer &numeqcoarse, Double &coarsealpha)
{
#ifdef TRACE
  (*trace) << "entering Acoustic3dPDE::SpecifySolver" << std::endl;
#endif

  conf->get("eps",eps,"Acoustic3d"); // relative accuracy in the precond. energy
  conf->get("dampiter",dampiter,"Acoustic3d"); // damping parameter for Jacobi, SSOR
  conf->get("maxnumit",maxnumit,"Acoustic3d"); // max number of iterations
  conf->get("solvertype",solvertype,"Acoustic3d"); // Richardson or CG
  conf->get("precondtype", precondtype, "Acoustic3d"); //ID or MG
  conf->get("numeqcoarse",numeqcoarse,"Acoustic3d"); // number of equation for coarsing
  conf->get("coarsealpha",coarsealpha,"Acoustic3d"); // coarsing parameter for AMG
}


void Acoustic3dPDE::SetMatrixFactors()
{
#ifdef TRACE
  (*trace) << "entering Acoustic3dPDE::SetMatrixFactors" << std::endl;
#endif

  matrix_factor_[0] = 1.0;
  matrix_factor_[1] = 1.0*a1_;
  matrix_factor_[2] = 0.0;
  matrix_factor_[3] = 1.0*a0_;
}

void Acoustic3dPDE::SpecifyMatrices(Integer &matrixtype, Integer * matrixsystype, Integer &graphtype, Integer &numdofpernode, Integer &numdirichlets,
Integer &numconstraints)
{
#ifdef TRACE
  (*trace) << "entering Acoustic3dPDE::SpecifyMatrices" << std::endl;
#endif

  matrixtype = 1;     // NOCLASS = 0
                      // RSPARSE = 1
                      // CSPARSE = 2
                      // RBLOCK  = 3
                      // CBLOCK  = 4
                      // RFULL   = 5
                      // CFULL   = 6
                      // MIXED   = 7

  /* matrixsystype: NOTYPE     = 0
                    SYSTEM     = 1
                    STIFFNESS  = 2
                    DAMPING    = 3
                    CONVECTION = 4
                    MASS       = 5
  */

  matrixsystype[0] = 1;   // memory for the system matrix
  matrixsystype[1] = 2;   // memory for the stiffness matrix
  matrixsystype[2] = 3;   // memory for the damping matrix
  matrixsystype[3] = 0;   // memory for the convection matrix
  matrixsystype[4] = 5;   // memory for the mass matrix

  graphtype = 1; // NOGRAPH = 0
                 // NODE   = 1
                 // EDGE   = 2
                 // FACE   = 3
                 // VOLUME = 4

  numdofpernode  = 1;
  numdirichlets = 1;
  numconstraints = 0;
}

void Acoustic3dPDE::SetupMatrices(const Integer level, BCs * ptBCs)
{
#ifdef TRACE
  (*trace) << "entering Acoustic3dPDE::SetupMatrices" << std::endl;
#endif

  Matrix<Double> elemmat;
  Matrix<Double> elemmatbnd; // This is a smaller matrix since it is just for linear 1D elements.
  Point<3> * ptCoord;

  BaseElem * ptEl;

  Integer matrix_stiff=2;
  Integer matrix_mass=5;
  Integer matrix_damp=3;


  Vector<Double> coeffm, coeffst, coeffdamp;
  CalcCoeff(coeffm, coeffst, coeffdamp);

  Vector<Integer> connecth;
  std::vector<Elem*> elemssd;

  Integer i, j;
  for (i=0; i<subdoms_.size(); i++)
    {
      ptgrid_->GetElemSD(elemssd,subdoms_[i],level);

      for (j=0; j< elemssd.size(); j++)
	{
	  ptEl=elemssd[j]->ptElem;

	  BaseForm<3> * bilinear_stiff = new LaplaceInt<3>(ptEl,1);
	  BaseForm<3> * bilinear_mass  = new MassInt<3>(ptEl,1);

	  Integer ii;
	  Integer elsize=(elemssd[j]->connect).size();
	  connecth.Resize(elsize);
	  for (ii=0; ii<elsize; ii++)
	    connecth[ii]=(elemssd[j]->connect)[ii];

	  ptCoord=new Point<3>[connecth.size()];
	  ptgrid_->GetCoordNodesElem(connecth,ptCoord,level);

	  // stiffness part
	  bilinear_stiff->CalcElemMatrix(ptCoord, elemmat);

	  elemmat*=coeffst[i];

#ifdef DEBUG
	  (*debug) << "Connection array  " << std::endl;
	  (*debug)  << connecth  << std::endl;

	  (*debug) << "Stiffnessmatrix, ElementNumber  " <<   j << std::endl;
	  (*debug) << elemmat << std::endl;


#endif

	  ptalgsys_->PutElemMatAlgSys(elemmat.getinarray(), connecth.get(), connecth.size(), as_sysid_, as_sysid_, matrix_stiff);

	  // mass part
	  bilinear_mass->CalcElemMatrix(ptCoord, elemmat);

	  elemmat*=coeffm[i];

#ifdef DEBUG
	  (*debug) << "Massmatrix, ElementNumber  " << j << std::endl;

	  (*debug) << elemmat << std::endl;
#endif

	  ptalgsys_->PutElemMatAlgSys(elemmat.getinarray(), connecth.get(), connecth.size(), as_sysid_, as_sysid_,matrix_mass);

	  delete bilinear_stiff;
	  delete bilinear_mass;
	  //	  delete [] ptCoord;
	}
    }

  // BEGIN DAMPING MATRIX PART
  if (!without_absBCs_) {
  std::vector<Elem*>  DomainBnd;
  conf->getliststr("bnd_for_absBCs",bnd_absBCs_,"Acoustic");
  DomainBnd=ptBCs->getFacesBC(bnd_absBCs_[0],level);

  for (j=0; j< DomainBnd.size(); j++)

    {
      ptEl=DomainBnd[j]->ptElem;

      BaseForm<3> * linear_damp = new DampInt<3>(ptEl,1);

      Integer ii;
      Integer elsize=(DomainBnd[j]->connect).size();
      connecth.Resize(elsize);
      for (ii=0; ii<elsize; ii++)
	connecth[ii]=(DomainBnd[j]->connect)[ii];

      ptCoord=new Point<3>[connecth.size()];
      ptgrid_->GetCoordNodesElem(connecth,ptCoord,level);

      linear_damp->CalcElemMatrix(ptCoord, elemmatbnd);
 
      elemmatbnd*=coeffdamp[0];
#ifdef DEBUG
      (*debug) << "Connection array  " << std::endl;
      (*debug)  << connecth  << std::endl;

      (*debug) << "Dampingmatrix, ElementNumber  " <<   j << std::endl;
      (*debug) << elemmatbnd << std::endl;
#endif

      ptalgsys_->PutElemMatAlgSys(elemmatbnd.getinarray(), connecth.get(), connecth.size(), as_sysid_, as_sysid_, matrix_damp);
      delete linear_damp;
      delete [] ptCoord;
    }
}      // END DAMPING MATRIX PART


#ifdef TRACE
  (*trace) << "Leaving Acoustic3dPDE::SetupMatrices" << std::endl;
#endif
}

void Acoustic3dPDE::SetBCs(BCs * ptBCs, const Integer level, const Integer update, const Double atime)
{
#ifdef TRACE
  (*trace) << "entering Acoustic3dPDE::SetBCs" << std::endl;
#endif

  Integer node;
  Double val, valTF;
  Double time=atime;

  //system matrix: id = 1
  Integer matrix_id = 1;

  // set dirichlet boundary conditions
  val=ptTimeFunc_->TimeFuncAtTime(time,level);
  std::cout<<"Time Function: "<<val<<std::endl;
  std::cout<<"atime: "<<time<<std::endl;

  Integer i,j=0;
  //Double val_hd = 0.0;
  std::list<Integer> nodes;
  Integer sizebc=bcs_hd_.size();
  std::list<Integer> nodes_hd; 
  std::cout << " size: " << bcs_hd_.size() << std::endl;

  //! homogeneous dirichlet BCs
  for (i=0; i< bcs_hd_.size(); i++)
    {
#ifdef DEBUG
       (*debug) << " Dirichlet BC" << std::endl;
#endif
     nodes_hd=ptBCs->GetNodesLevel(bcs_hd_[i]);

     for (std::list<Integer>::const_iterator p=nodes_hd.begin(); p!=nodes_hd.end(); p++, j++)

	{
	  node=*p;
	  if (update==1)
	    {
#ifdef DEBUG
       (*debug) << " node: " << node << " val: " << val << std::endl;
#endif
            ptalgsys_->SetBCDirichletUpdate(j+1, node, val, dofspernode_, as_sysid_, as_sysid_, matrix_id);
	    }
	  else
          {
#ifdef DEBUG
      (*debug) << " node: " << node << " val: " << val << std::endl;
#endif
	   ptalgsys_->SetBCDirichlet(j+1, node, val, dofspernode_, as_sysid_, as_sysid_, matrix_id);
	    }
	}
    }

//   //! inhomogeneous dirichlet BCs
//   //valTF=ptTimeFunc_->TimeFuncAtTime(atime,level);
//  for (i=0; i<bcs_id_.size(); i++)
//    {
//      nodes=ptBCs->GetNodesLevel(bcs_id_[i], level);

//      val=val_id_[i]*valTF;

//      for (std::list<Integer>::const_iterator p=nodes.begin(); p!=nodes.end(); p++, j++)
// 	{
// 	  node=*p;
//          if (update==1)
//            {
//              ptalgsys_->SetBCDirichletUpdate(j+1, node, val, dofspernode_, as_sysid_, as_sysid_, matrix_id);
//            }
//          else
//            {
//              ptalgsys_->SetBCDirichlet(j+1, node, val, dofspernode_, as_sysid_,as_sysid_, matrix_id);
//            }
// 	}
//    }


#ifdef TRACE
  (*trace) << "leaving Acoustic3dPDE::SetBCs" << std::endl;
#endif
}


void Acoustic3dPDE::preComputeRHS()
{
#ifdef TRACE
  (*trace) << "entering Acoustic3dPDE::preComputeRHS" << std::endl;
#endif

  SetRHSFnc = FALSE;
   //  if (conf->is_there("rhs_surfaces")) {
   std::string isthererhs;
   conf->get("load_force",isthererhs,"Acoustic3d");
//  }
}

void Acoustic3dPDE::ComputeRHS(const Double atime, BCs * ptBCs)
{
#ifdef TRACE
  (*trace) << "entering Acoustic3dPDE::ComputeRHS" << std::endl;
#endif

  Integer n;
  Integer node;
  Integer matrix_id;
  Vector<Double> coeffMass, coeffDamp;
	Vector<Double> elemvec;
  std::list<Integer> nodes_hd;
  Integer i;
  matrix_id = 5;
  coeffMass = sol_old_*a0_+sol_der1_old_*a2_+sol_der2_old_*a3_;

  ptalgsys_->UpdateRHS(as_sysid_,as_sysid_,matrix_id,coeffMass.get());

  // damping matrix part
  if (!without_absBCs_) {
  matrix_id = 3; // Setting id for Damping matrix	
  coeffDamp = -sol_der1_old_-sol_der2_old_*a6_;
  
  ptalgsys_->UpdateRHS(as_sysid_,as_sysid_,matrix_id,coeffDamp.get());

  matrix_id = 3; // Setting id for Damping matrix	
  coeffDamp = sol_old_*a1_+sol_der1_old_*a2_*a7_+sol_der2_old_*a7_*a3_;  

  ptalgsys_->UpdateRHS(as_sysid_,as_sysid_,matrix_id,coeffDamp.get());
  }
   
}


void Acoustic3dPDE::SolveStepTrans(BCs * ptBCs, const Integer kstep, const Double asteptime, const Integer level, const Boolean reset)
{
#ifdef TRACE
  (*trace) << "entering Acoustic3dPDE::SolveStepTrans" << std::endl;
#endif

  lasttimecalc_= asteptime;
  Boolean Recalc=FALSE;

  if (laststepcalc_==kstep && kstep!=0) Recalc=TRUE;
  else laststepcalc_= kstep;

  Double * ptsol;
  Integer update,job;

  if (kstep==0)
    {
      update = 0;
      job = 1;
      SetupMatrices(level,ptBCs);
      ptalgsys_->ComputeSysMatrix(as_sysid_,as_sysid_,matrix_factor_);    
      ptalgsys_->ResetRHS(as_sysid_);
      ComputeRHS(lasttimecalc_,ptBCs);
    }
  else if (reset)
    {
      update = 1;
      job    = 1;

      ptalgsys_->ResetRHS(as_sysid_);
      ptalgsys_->ResetMatrix(0,0,1);
      ptalgsys_->ComputeSysMatrix(as_sysid_,as_sysid_,matrix_factor_);
      ComputeRHS(lasttimecalc_,ptBCs);
    }
  else
    {
      update = 1;
      job    = 3;
      ptalgsys_->ResetRHS(as_sysid_);
      ComputeRHS(lasttimecalc_,ptBCs);
    };

  SetBCs(ptBCs,level,update,lasttimecalc_);
  ptalgsys_->ComputePrecond(job,as_sysid_);
  ptalgsys_->SolveAlgSys(as_sysid_);
  ptsol = ptalgsys_->GetSolution(as_sysid_);

  // save solution
  Integer i;
  for (i=0; i<ptgrid_->GetMaxnumnodes(level); i++)
    sol_[i]=ptsol[i];

  // calculation of derivatives of solution
  CalcDerSol(); 

}


void Acoustic3dPDE::WriteResultsInFile()
{
#ifdef TRACE
  (*trace) << "entering Acoustic3dPDE::WriteResultsInFile" << std::endl;
#endif

  if (OutFile_->IsGMV())
    {
      OutFile_->WriteSolution(sol_,laststepcalc_,lasttimecalc_,"vp");
      OutFile_->WriteSolution(sol_der1_,laststepcalc_,lasttimecalc_,"vp_1der");
      OutFile_->WriteSolution(sol_der2_,laststepcalc_,lasttimecalc_,"vp_2der");
    }
  else
    {
            OutFile_->WriteSolution(sol_,laststepcalc_,lasttimecalc_,"fluid potential");
	    //OutFile_->WriteSolution(sol_der1_,laststepcalc_,lasttimecalc_,"fluid potential, 1st deriv., ");
      //      OutFile_->WriteSolution(sol_der2_,laststepcalc_,lasttimecalc_,"fluid potential, 2nd deriv., ");
    }
}

void Acoustic3dPDE :: CalcParameters(const Double dt)
{
#ifdef TRACE
  (*trace) << "entering Acoustic3dPDE::CalcParameters" << std::endl;
#endif

  a2_=1.0/(beta_*dt);
  a0_=a2_*(1/dt);
  a1_=gamma_*a2_;
  a3_=0.5/(beta_)-1.0;
  a4_=gamma_/beta_-1.0;
  a5_=0.5*dt*(a4_-1.0);
  a6_=dt*(1-gamma_);
  a7_=gamma_*dt;
}

void Acoustic3dPDE :: CalcDerSol()
{
#ifdef TRACE
  (*trace) << " entering  Acoustic3dPDE :: CalcDerSol() " << std::endl;
#endif

  sol_der2_=(sol_ - sol_old_)*a0_ - (sol_der1_old_)*a2_ - sol_der2_old_*a3_;
  sol_der1_=sol_der1_old_+sol_der2_old_*a6_+sol_der2_*a7_;
}

void Acoustic3dPDE::SaveSolAsPrevStep()
{
#ifdef TRACE
  (*trace) << " entering  Acoustic3dPDE::SaveSolAsPrevStep() " << std::endl;
#endif

  sol_old_=sol_;
  sol_der1_old_=sol_der1_;
  sol_der2_old_=sol_der2_;
}


void Acoustic3dPDE::CalcCoeff(Vector<Double> & coeffmass, Vector<Double> & coeffstiff, Vector<Double> & coeffdamp)
{
#ifdef TRACE
  (*trace) << " entering Acoustic3dPDE::CalcCoeff() " << std::endl;
#endif
  if (!MatFile_) Error("You didn't specialize material file. Check your config-file.");

  coeffmass.Resize(subdoms_.size());
  coeffstiff.Resize(subdoms_.size());
  if (!without_absBCs_)
  coeffdamp.Resize(subdoms_.size());
  
  Integer i,matnum;
  for (i=0; i<subdoms_.size(); i++)
    {
			conf->get(subdoms_[i],matnum,"list_subdomains");
			if (matnum==11111) break; //This may be commented out.
      // read density and compress with material number matnum
      Double density, compress;
      MatFile_->ReadDensityAndCompressity(density,compress,matnum,"fluid");

      coeffmass[i]  = density*density/compress;
      coeffstiff[i] = density;
 if (!without_absBCs_)
      coeffdamp[i]  = density/((sqrt(compress/density)));
    }
}

Acoustic3dPDE::~Acoustic3dPDE()
{
  ;
}



} // end of namespace

