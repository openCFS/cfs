#include <fstream>
#include <iostream>
#include <string>
#include <math.h>

#include "elec3dPDE.hh"
#include <DataInOut/Unverg/outUnverg.hh>
#include <Forms/forms_header.hh>
#include <Estimator/spaceerror.hh>
 
namespace CoupledField
{

Elec3dPDE::Elec3dPDE(AbstractAlgebraicSys * ptalgsys, Grid * aptgrid, Material * ptMaterial, TimeFunc * aptTimeFunc, 
		     FileType * aptFileType, WriteResults * aptOut)
:BasePDE(ptalgsys,ptMaterial,aptFileType,aptOut,aptTimeFunc)
{
#ifdef TRACE
  (*trace) << "entering Elec3dPDE::Elec3dPDE " << std::endl;
#endif

  dofspernode_=1;
  pdename_    = "Electric3d";

  ptgrid_=aptgrid;

  size_=ptgrid_->GetMaxnumnodes(0);
  sol_.Resize(size_);
  sol_.Init(0);

  conf->getsubdompde(subdoms_,pdename_);
  ReadBCs(pdename_);
}


void Elec3dPDE::SetMatrixFactors()
{
#ifdef TRACE
  (*trace) << "entering Elec3dPDE::SetMatrixFactors" << std::endl;
#endif
 
  matrix_factor_[0] = 1.0;  // factor for stiffness matrix
  matrix_factor_[1] = 0.0;  // factor for damping matrix
  matrix_factor_[2] = 0.0;  // factor for convection matrix
  matrix_factor_[3] = 0.0;  // factor for mass matrix
}

void Elec3dPDE::SpecifyMatrices(Integer &matrixtype, Integer * matrixsystype, Integer &graphtype, 
				Integer &numdofpernode, Integer &numdirichlets, Integer &numconstraints)
{
#ifdef TRACE
  (*trace) << "entering Elec3dPDE::SpecifyMatrices" << std::endl;
#endif

  matrixtype = RSCALAR; 
  graphtype  = NODEGRAPH; 

  matrixsystype[0] = SYSTEM;      // memory for the system matrix
  matrixsystype[1] = STIFFNESS;   // memory for the stiffness matrix
  matrixsystype[2] = 0;           // memory for the damping matrix
  matrixsystype[3] = 0;           // memory for the convection matrix
  matrixsystype[4] = 0;           // memory for the mass matrix


  numdofpernode  = dofspernode_;
  numdirichlets = 1;
  numconstraints = 0;
}

void Elec3dPDE::SetupMatrices(const Integer level, BCs * ptBCs)
{
#ifdef TRACE
  (*trace) << "entering Elec3dPDE::SetupMatrices" << std::endl;
#endif
  
  Matrix<Double> elemmat;
  Point<3> * ptCoord;

  BaseElem * ptElem;

  Vector<Double> coeffst;
  CalcCoeff(coeffst);  

  Vector<Integer> connecth;
  std::vector<Elem*> elemssd;

  Integer i, j;
  for (i=0; i<subdoms_.size(); i++)
    {
      ptgrid_->GetElemSD(elemssd,subdoms_[i],level);

      for (j=0; j < elemssd.size(); j++)
	{  
	  ptElem=elemssd[j]->ptElem;

	  BaseForm<3> * bilinear_stiff = new LaplaceInt<3>(ptElem,1);

	  connecth=elemssd[j]->connect;

	  ptCoord=new Point<3>[connecth.size()];
	  ptgrid_->GetCoordNodesElem(connecth,ptCoord,level);

	  // stiffness part
	  bilinear_stiff->CalcElemMatrix(ptCoord, elemmat);
	  elemmat*=coeffst[i];

	  if (InfoPrint)
	    (*infofile) << elemmat << std::endl;

#ifdef DEBUG
	  (*debug) << "Stiffnessmatrix, ElementNumber  " <<   i << std::endl;

	  (*debug) << elemmat << std::endl;
#endif

	  ptalgsys_->PutElemMatAlgSys(elemmat.getinarray(), connecth.get(), connecth.size(), as_sysid_, as_sysid_, 
				      STIFFNESS);

	  delete bilinear_stiff;
	  delete [] ptCoord;
	}     
    }

#ifdef TRACE
  (*trace) << "Leaving Elec3dPDE::SetupMatrices" << std::endl;
#endif
}

void Elec3dPDE::SetBCs(BCs * ptBCs, const Integer level, const Integer update, const Double atime)
{
#ifdef TRACE
  (*trace) << "entering Elec3dPDE::SetBCs" << std::endl;
#endif

  Integer node;
  Double val, valueTF;

  Integer i,j=0;
  std::list<Integer> nodes;

  for (i=0; i<bcs_hd_.size(); i++)
    {  
      nodes=ptBCs->GetNodesLevel(bcs_hd_[i], level);
  
      for (std::list<Integer>::const_iterator p=nodes.begin(); p!=nodes.end(); p++, j++)
	{
	  node=*p;
	  val=0; 
          if (update==1)
            {

	      if (InfoPrint)
		(*infofile) << " node: " << node << " val: " << val << std::endl;

              ptalgsys_->SetBCDirichletUpdate(j+1, node, val, dofspernode_, as_sysid_, 
					      as_sysid_, SYSTEM);
            }
          else
            {

	      if (InfoPrint)
		(*infofile) << " node: " << node << " val: " << val << std::endl;

              ptalgsys_->SetBCDirichlet(j+1, node, val, dofspernode_, as_sysid_,
					as_sysid_, SYSTEM);
            }
	}
    }

  for (i=0; i<bcs_id_.size(); i++)
    {
      nodes=ptBCs->GetNodesLevel(bcs_id_[i], level);

      val=val_id_[i];

      for (std::list<Integer>::const_iterator p=nodes.begin(); p!=nodes.end(); p++, j++)
	{
	  node=*p;
          if (update==1)
            {

	      if (InfoPrint)
		(*infofile) << " node: " << node << " val: " << val << std::endl;

              ptalgsys_->SetBCDirichletUpdate(j+1, node, val, dofspernode_, as_sysid_, 
					      as_sysid_, SYSTEM);
            }
          else
            {

	      if (InfoPrint)
		(*infofile) << " node: " << node << " val: " << val << std::endl;
	      
              ptalgsys_->SetBCDirichlet(j+1, node, val, dofspernode_, as_sysid_,as_sysid_, 
					SYSTEM);
            }
	}
    }
  
#ifdef TRACE
  (*trace) << " leaving Elec3d::ReadBCs " << std::endl;
#endif
 
}

void Elec3dPDE::SolveStepStatic(BCs * ptBCs, Integer level)
{
#ifdef TRACE
  (*trace) << "entering Elec3dPDE::SolveStepStatic" << std::endl;
#endif

  Integer update = 0;
  Integer job = 1;

  Double * ptsol;

  SetupMatrices(level);
  ptalgsys_->ComputeSysMatrix(as_sysid_,as_sysid_,matrix_factor_);

  SetBCs(ptBCs,level,update,0);

  ptalgsys_->ComputePrecond(job,as_sysid_);
  ptalgsys_->SolveAlgSys(as_sysid_);

  ptsol = ptalgsys_->GetSolution(as_sysid_);

  // save solution
  Vector<Double> transsol(ptgrid_->GetMaxnumnodes(level), ptsol);
  sol_=transsol;

}

void Elec3dPDE::CalcCoeff(Vector<Double> & coeff)
{
  if (!MatFile_) Error("You didn't specialize material file. Use option -m");

  coeff.Resize(subdoms_.size());

  Integer i, matnum;
  for (i=0; i<subdoms_.size(); i++)
    {
      conf->get(subdoms_[i],matnum,"list_subdomains");

      Double dielectr;
      MatFile_->ReadDielectricTerms(dielectr,matnum); 

      coeff[i]=dielectr;
    }
}

void Elec3dPDE::WriteResultsInFile()
{
#ifdef TRACE
  (*trace) << "entering Elec3dPDE::WriteResultsInFile" << std::endl;
#endif

  Integer step=0;
  Double time=0;

  // write results
  if (OutFile_->IsGMV())
    OutFile_->WriteSolution(sol_,step,time,"electric_potential");
  else
    OutFile_->WriteSolution(sol_,step,time,"electric potential");

}

Elec3dPDE::~Elec3dPDE()
{
 ;
}

} // end of namespace


