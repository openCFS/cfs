#include <fstream>
#include <iostream>
#include <string>
#include <math.h>

#include "elec2dPDE.hh"
#include <DataInOut/Unverg/outUnverg.hh>
#include <Forms/forms_header.hh>
#include <AlgebraicSystem/LinAlg/linsystem.hh>
#include <Estimator/spaceerror.hh>

namespace CoupledField
{

Elec2dPDE::Elec2dPDE(AbstractAlgebraicSys * ptalgsys, Grid * aptgrid, Material * ptMaterial, 
		     TimeFunc * aptTimeFunc, FileType * aptFileType, WriteResults * aptOut)
:BasePDE(ptalgsys,ptMaterial,aptFileType,aptOut,aptTimeFunc)
{
#ifdef TRACE
  (*trace) << "entering Elec2dPDE::Elec2dPDE " << std::endl;
#endif

  dofspernode_=1;
  pdename_    = "Electric2d";

  ptgrid_=aptgrid;

  size_=ptgrid_->GetMaxnumnodes(0);
  sol_.Resize(size_);
  sol_.Init(0);

  conf->getsubdompde(subdoms_,pdename_);
  ReadBCs(pdename_);
}

void Elec2dPDE::SetMatrixFactors()
{
#ifdef TRACE
  (*trace) << "entering Elec2dPDE::SetMatrixFactors" << std::endl;
#endif
  
  matrix_factor_[0] = 1.0;  // factor for stiffness matrix
  matrix_factor_[1] = 0.0;  // factor for damping matrix
  matrix_factor_[2] = 0.0;  // factor for convection matrix
  matrix_factor_[3] = 0.0;  // factor for mass matrix
}

void Elec2dPDE::SpecifyMatrices(Integer &matrixtype, Integer * matrixsystype, Integer &graphtype, 
				Integer &numdofpernode, Integer &numdirichlets, Integer &numconstraints)
{
#ifdef TRACE
  (*trace) << "entering Elec2dPDE::SpecifyMatrices" << std::endl;
#endif

  matrixtype = RSCALAR; 
  graphtype  = NODEGRAPH; 

  matrixsystype[0] = SYSTEM;      // memory for the system matrix
  matrixsystype[1] = STIFFNESS;   // memory for the stiffness matrix
  matrixsystype[2] = 0;           // memory for the damping matrix
  matrixsystype[3] = 0;           // memory for the convection matrix
  matrixsystype[4] = 0;           // memory for the mass matrix


  numdofpernode  = dofspernode_;
  numdirichlets  = 1;
  numconstraints = 0;
}

void Elec2dPDE::SetupMatrices(const Integer level, BCs * ptBCs)
{
#ifdef TRACE
  (*trace) << "entering Elec2dPDE::SetupMatrices" << std::endl;
#endif
  
  Matrix<Double> elemmat;
  Point<2> * ptCoord;

  BaseElem * ptElem;

  if (InfoPrint)
    (*infofile) << " ------------------------- Element matrices --------------- " << std::endl;

  Vector<Double> coeffst;
  
  CalcCoeff(coeffst);  

  Integer i, j;
  for (i=0; i<subdoms_.size(); i++)
    {
      std::vector<Elem*> elemssd;
   
      ptgrid_->GetElemSD(elemssd,subdoms_[i],level);

      for (j=0; j < elemssd.size(); j++)
	{  
	  Vector<Integer> connecth;
 
	  ptElem=elemssd[j]->ptElem;

	  BaseForm<2> * bilinear_stiff = new LaplaceInt<2>(ptElem,1);

	  connecth=elemssd[j]->connect;

	  ptCoord=new Point<2>[connecth.size()];
	  ptgrid_->GetCoordNodesElem(connecth,ptCoord,level);

	  // stiffness part
	  bilinear_stiff->CalcElemMatrix(ptCoord, elemmat);
	  elemmat*=coeffst[i];
	  
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
  (*trace) << "Leaving Elec2dPDE::SetupMatrices" << std::endl;
#endif
}

void Elec2dPDE::SetBCs(BCs * ptBCs, const Integer level, const Integer update, const Double atime)
{
#ifdef TRACE
  (*trace) << "entering Elec2dPDE::SetBCs" << std::endl;
#endif

  Integer node;
  Double val, valueTF;

  Integer j=0;
  std::list<Integer> nodes;

  if (InfoPrint)
    (*infofile) << " ---------------- Dirichle boundary condition -------------" << std::endl;

  Integer i;
  for (i=0; i<bcs_hd_.size(); i++)
    {  
      nodes=ptBCs->GetNodesLevel(bcs_hd_[i]);
  
      for (std::list<Integer>::const_iterator p=nodes.begin(); p!=nodes.end(); p++, j++)
	{
	  node=*p;
	  val=0; 
          if (update==1)
            {

	      if (InfoPrint)
		(*infofile) << " node: " << node << " val: " << val << std::endl;

              ptalgsys_->SetBCDirichletUpdate(j+1, node, val, dofspernode_, as_sysid_, as_sysid_, SYSTEM);
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
      nodes=ptBCs->GetNodesLevel(bcs_id_[i]);

      val=val_id_[i];

      for (std::list<Integer>::const_iterator p=nodes.begin(); p!=nodes.end(); p++, j++)
	{
	  node=*p;
          if (update==1)
            {	     

	      if (InfoPrint)
		(*infofile) << " node: " << node << " val: " << val << std::endl;

              ptalgsys_->SetBCDirichletUpdate(j+1, node, val, dofspernode_, as_sysid_, as_sysid_, SYSTEM);
            }
          else
            {
	      
	      if (InfoPrint)
		(*infofile) << " node: " << node << " val: " << val << std::endl;
	      
	      ptalgsys_->SetBCDirichlet(j+1, node, val, dofspernode_, as_sysid_,as_sysid_, SYSTEM);
	   
            }
	}
    }
  
#ifdef TRACE
  (*trace) << " leaving Elec2d::SetBCs " << std::endl;
#endif 
}

void Elec2dPDE::SolveStepStatic(BCs * ptBCs, Integer level)
{
#ifdef TRACE
  (*trace) << "entering Elec2dPDE::SolveStepStatic" << std::endl;
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

void Elec2dPDE::WriteResultsInFile()
{
#ifdef TRACE
  (*trace) << "entering Elec2dPDE::WriteResultsInFile" << std::endl;
#endif
 
  Integer step=0;
  Double time=0;
  if (OutFile_->IsGMV()) 
    OutFile_->WriteSolution(sol_,step,time,"electric_potential");
  else 
    OutFile_->WriteSolution(sol_,step,time,"electric potential");

}

void Elec2dPDE::CalcCoeff(Vector<Double> & coeff)
{
#ifdef TRACE
  (*trace) << " entering Elec2dPDE::CalcCoeff " <<std::endl;
#endif

  if (!MatFile_) Error("You didn't specialize material file.");

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

Elec2dPDE::~Elec2dPDE()
{
#ifdef TRACE
  (*trace) << " entering Elec2dPDE::~Elec2dPDE " << std::endl;
#endif

}

} // end of namespace
