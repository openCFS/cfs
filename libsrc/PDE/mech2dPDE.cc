#include <fstream>
#include <iostream>
#include <string>
#include <math.h>

#include "mech2dPDE.hh"
#include <Estimator/actimeerror.hh>
#include <DataInOut/Unverg/outUnverg.hh>
#include <DataInOut/GMV/outGMV.hh>
#include <Forms/forms_header.hh>
#include <Estimator/spaceerror.hh>



namespace CoupledField
{


Mech2dPDE::Mech2dPDE(AbstractAlgebraicSys * ptalgsys, Grid * aptgrid, Material * ptMaterial, TimeFunc * aptTimeFunc, FileType * aptFileType, WriteResults * aptOut)
:BasePDE(ptalgsys,ptMaterial,aptFileType,aptOut,aptTimeFunc)
{
#ifdef TRACE
  (*trace) << "entering Mech2dPDE::Mech2dPDE " << std::endl;
#endif

  dofspernode_=2;
  ptgrid_=aptgrid; 

  laststepcalc_=0;
  size_=ptgrid_->GetMaxnumnodes(0);

  sol_.Resize(size_*dofspernode_);
  sol_.Init(0);

  sol_der1_.Resize(size_*dofspernode_);
  sol_der1_.Init(0);

  sol_der2_.Resize(size_*dofspernode_);
  sol_der2_.Init(0);

  //  conf->get("alpha_NM",alpha_,"Mech2d"); 
  //  conf->get("beta_NM",beta_,"Mech2d"); 
  //  conf->get("gamma_NM",gamma_,"Mech2d");

  conf->getsubdompde(subdoms_,"Mech2d");
  ReadBCs("Mech2d");

}

void Mech2dPDE::SpecifySolver(Integer &solvertype, Integer &precondtype, Double &eps, Double &dampiter, Integer &maxnumit, Integer &numeqcoarse, Double &coarsealpha)
{
#ifdef TRACE
  (*trace) << "entering Mech2dPDE::SpecifySolver" << std::endl;
#endif

  conf->get("eps",eps,"Mech2d"); // relative accuracy in the precond. energy
  conf->get("dampiter",dampiter,"Mech2d"); // damping parameter for Jacobi, SSOR
  conf->get("maxnumit",maxnumit,"Mech2d"); // max number of iterations
  conf->get("solvertype",solvertype,"Mech2d"); // Richardson or CG
  conf->get("precondtype", precondtype, "Mech2d"); // ID or MG
  conf->get("numeqcoarse",numeqcoarse,"Mech2d"); // number of equation for coarsing
  conf->get("coarsealpha",coarsealpha,"Mech2d"); // coarsing parameter for AMG
}


void Mech2dPDE::SetMatrixFactors()
{
#ifdef TRACE
  (*trace) << "entering Mech2dPDE::SetMatrixFactors" << std::endl;
#endif
 
  matrix_factor_[0] = 1.0;  // stiffness
  matrix_factor_[1] = 0.0;  // damping
  matrix_factor_[2] = 0.0;  //convection
  matrix_factor_[3] = 0.0; //mass
}

void Mech2dPDE::SpecifyMatrices(Integer &matrixtype, Integer * matrixsystype, Integer &graphtype, Integer &numdofpernode, Integer &numdirichlets,
Integer &numconstraints)
{
#ifdef TRACE
  (*trace) << "entering Mech2dPDE::SpecifyMatrices" << std::endl;
#endif

  matrixtype = RBLOCK;

  matrixsystype[0] = 1;   // memory for the system matrix
  matrixsystype[1] = 2;   // memory for the stiffness matrix
  matrixsystype[2] = 0;   // memory for the damping matrix
  matrixsystype[3] = 0;   // memory for the convection matrix
  matrixsystype[4] = 5;   // memory for the mass matrix

  graphtype = NODEGRAPH;

  numdofpernode  = dofspernode_;

  //not used; just for interface
  numdirichlets  = 1;
  numconstraints = 0;
}

void Mech2dPDE::SetupMatrices(const Integer level, BCs * ptBCs)
{
#ifdef TRACE
  (*trace) << "entering Mech2dPDE::SetupMatrices" << std::endl;
#endif

  Matrix<Double> elemmat;
  Point<2> * ptCoord; 

  BaseElem * ptEl;

  //get material tensor 
  Vector<Double> coeffm, coeffst;
  //  CalcCoeff(coeffm, coeffst);

  Vector<Integer> connecth;
  std::vector<Elem*> elemssd;

  Integer i, j;
  for (i=0; i<subdoms_.size(); i++)
    {
      ptgrid_->GetElemSD(elemssd,subdoms_[i],level);

      for (j=0; j< elemssd.size(); j++)
	{
	  ptEl=elemssd[j]->ptElem;

	  BaseForm<2> * bilinear_stiff = new LaplaceInt<2>(ptEl,1);
	  BaseForm<2> * bilinear_mass  = new MassInt<2>(ptEl,dofspernode_);

	  Integer ii;
	  Integer elsize=(elemssd[j]->connect).size();
	  connecth.Resize(elsize);
	  for (ii=0; ii<elsize; ii++)
	    connecth[ii]=(elemssd[j]->connect)[ii];
  
	  ptCoord=new Point<2>[connecth.size()];
	  ptgrid_->GetCoordNodesElem(connecth,ptCoord,level); 

	  // stiffness part
	  bilinear_stiff->CalcElemMatrix(ptCoord, elemmat);

          Integer blocksize = ptEl->GetNumNodes()*dofspernode_;

          Matrix<Double> blockelemmat(blocksize,blocksize);
	  blockelemmat.Init();
	  Integer off = connecth.size();
          for (Integer kr=0;kr<off;kr++)
	    for (Integer kc=0;kc<off;kc++)
	      {
		blockelemmat[kr*2][kc*2] = elemmat[kr][kc];
		blockelemmat[kr*2+1][kc*2+1] = elemmat[kr][kc];
	      }
          
	  // multiply with material
	  //	  elemmat*=coeffst[i];

#ifdef DEBUG
	  (*debug) << "Connection array  " << std::endl;
	  (*debug)  << connecth  << std::endl;

	  (*debug) << "Stiffnessmatrix, ElementNumber  " <<   j << std::endl;
	  (*debug) << blockelemmat << std::endl;


#endif     

	  ptalgsys_->PutElemMatAlgSys(blockelemmat.getinarray(), connecth.get(), connecth.size(), as_sysid_, 
				      as_sysid_, SYSTEM);

	  // mass part
	  bilinear_mass->CalcElemMatrix(ptCoord, elemmat);

	  //	  elemmat*=coeffm[i];

	  blockelemmat.Init();
          for (Integer kr=0;kr<off;kr++)
	    for (Integer kc=0;kc<off;kc++)
	      {
		blockelemmat[kr][kc] = elemmat[kr][kc];
		blockelemmat[kr+off][kc+off] = elemmat[kr][kc];
	      }

#ifdef DEBUG
	  (*debug) << "Massmatrix, ElementNumber  " << j << std::endl;

	  (*debug) <<blockelemmat << std::endl;
#endif
      
	  //	  ptalgsys_->PutElemMatAlgSys(blockelemmat.getinarray(), connecth.get(), connecth.size(), as_sysid_, 
	  //				      as_sysid_,MASS);

  
	  delete bilinear_stiff;
	  delete bilinear_mass;
	  delete [] ptCoord;
	}
    }
#ifdef TRACE
  (*trace) << "Leaving Mech2dPDE::SetupMatrices" << std::endl;
#endif
}

void Mech2dPDE::SetBCs(BCs * ptBCs, const Integer level, const Integer update, const Double atime)
{
#ifdef TRACE
  (*trace) << "entering Mech2dPDE::SetBCs" << std::endl;
#endif

  Integer node,dof;
  Double val;

  std::string dof_ux = "ux";
  std::string dof_uy = "uy";

  // set dirichlet boundary conditions
  Integer i,j=0;
  std::list<Integer> nodes;
  Integer sizebc=bcs_hd_.size();

  // homogeneous boundary conditions
  val = 0;
  for (i=0; i< bcs_hd_.size(); i++)
    {
      std::string doftype = bcs_hd_[i]; 
      if (doftype.substr(0,2)==dof_ux)
	dof = 1;
      else 
	if (doftype.substr(0,2)==dof_uy)
	dof = 2;
      else
	Error("Unknown dof-type in homog. BC; substring is:must start with ux or uy!!",__FILE__,__LINE__);
      
#ifdef DEBUG
      (*debug) << " Homog. Dirichlet BC" << std::endl;
#endif
      nodes=ptBCs->GetNodesLevel(bcs_hd_[i]);
   
      for (std::list<Integer>::const_iterator p=nodes.begin(); p!=nodes.end(); p++, j++)
	{
	  node=*p;
	  if (update==1)
	    {
#ifdef DEBUG
	      (*debug) << " node: " << node << " dof:" << dof << " val: " << val << std::endl;
#endif
      
      ptalgsys_->SetBCDirichletUpdate(i+1, node, val, dof, as_sysid_, as_sysid_, SYSTEM);

	    }
	  else
	    { 
#ifdef DEBUG
	      (*debug) << " node: " << node << " dof:" << dof << " val: " << val << std::endl;
#endif

      ptalgsys_->SetBCDirichlet(j+1, node, val, dof, as_sysid_, as_sysid_, SYSTEM);

	    }
	}  
    }

  //inhomogeneous boundary conditions
  for (i=0; i<bcs_id_.size(); i++)
    {
      std::string doftype = bcs_hd_[i]; 
      if (doftype.substr(0,2)==dof_ux)
	dof = 1;
      else 
	if (doftype.substr(0,2)==dof_uy)
	dof = 2;
	else
	  Error("Unknown dof-type in inhomog. BC; substring is:must start with ux or uy!!",__FILE__,__LINE__);
      
#ifdef DEBUG
      (*debug) << " Inhomog. Dirichlet BC" << std::endl;
#endif

      nodes=ptBCs->GetNodesLevel(bcs_id_[i], level);

      val=val_id_[i];

      for (std::list<Integer>::const_iterator p=nodes.begin(); p!=nodes.end(); p++, j++)
	{
	  node=*p;
          if (update==1)
            {

#ifdef DEBUG
	      (*debug) << " node: " << node << " dof:" << dof << " val: " << val << std::endl;
#endif

              ptalgsys_->SetBCDirichletUpdate(j+1, node, val, dof, as_sysid_, as_sysid_, SYSTEM);
            }
          else
            {

#ifdef DEBUG
	      (*debug) << " node: " << node << " dof:" << dof << " val: " << val << std::endl;
#endif
              ptalgsys_->SetBCDirichlet(j+1, node, val, dof, as_sysid_,as_sysid_, SYSTEM);
            }
	}
    }

#ifdef TRACE
  (*trace) << "leaving Mech2dPDE::SetBCs" << std::endl;
#endif
}


void Mech2dPDE::SolveStepStatic(BCs * ptBCs, Integer level)
{
#ifdef TRACE
  (*trace) << "entering Mech2dPDE::SolveStepStatic" << std::endl;
#endif

  Integer update = 0;
  Integer job = 1;

  SetupMatrices(level);
  SetBCs(ptBCs,level,update,0);
  //  ptalgsys_->ComputeSysMatrix(as_sysid_,as_sysid_,matrix_factor_);
  ptalgsys_->ComputePrecond(job,as_sysid_);
  ptalgsys_->SolveAlgSys(as_sysid_);

  // get solution
  Double * ptsol;
  ptsol = ptalgsys_->GetSolution(as_sysid_);

  // save solution
  Vector<Double> transsol(ptgrid_->GetMaxnumnodes(level)*dofspernode_, ptsol);
  sol_=transsol;
}

void Mech2dPDE::SolveStepTrans(BCs * ptBCs, const Integer kstep, const Double asteptime, const Integer level, const Boolean reset)
{
#ifdef TRACE
  (*trace) << "entering Mech2dPDE::SolveStepTrans" << std::endl;
#endif

}

void Mech2dPDE::WriteResultsInFile()
{
#ifdef TRACE
  (*trace) << "entering Mech2dPDE::WriteResultsInFile" << std::endl;
#endif

  if (OutFile_->IsGMV())
    {
      OutFile_->WriteSolution(sol_,laststepcalc_,lasttimecalc_,"disp");
      //      OutFile_->WriteSolution(sol_der1_,laststepcalc_,lasttimecalc_,"vp_1der");
      //      OutFile_->WriteSolution(sol_der2_,laststepcalc_,lasttimecalc_,"vp_2der");
    }
  else
    {
      OutFile_->WriteSolution(sol_,laststepcalc_,lasttimecalc_,"disp");
      //      OutFile_->WriteSolution(sol_der1_,laststepcalc_,lasttimecalc_,"fluid potential, 1st deriv., ");
      //      OutFile_->WriteSolution(sol_der2_,laststepcalc_,lasttimecalc_,"fluid potential, 2nd deriv., ");
    }
  
}

void Mech2dPDE :: CalcParameters(const Double dt)
{
#ifdef TRACE
  (*trace) << "entering Mech2dPDE::CalcParameters" << std::endl;
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

void Mech2dPDE::CalcCoeff(Vector<Double> & coeffmass, Vector<Double> & coeffstiff)
{
  if (!MatFile_) Error("You didn't specialize material file. Check your config-file.");

  coeffmass.Resize(subdoms_.size());
  coeffstiff.Resize(subdoms_.size());
 
  Integer i,matnum;
  for (i=0; i<subdoms_.size(); i++)
    {
      conf->get(subdoms_[i],matnum,"list_subdomains");

      // read density and compress with material number matnum
      Double density, compress;
      MatFile_->ReadDensityAndCompressity(density,compress,matnum,"fluid");

      coeffmass[i]  = density*density/compress;
      coeffstiff[i] = density;
    }
}
}
