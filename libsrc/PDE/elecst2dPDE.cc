#include <fstream>
#include <iostream>
#include <string>
#include <math.h>

#include "elecst2dPDE.hh"
#include "elecst2dErr.hh"
#include "outUnverg.hh"
#include "forms_header.hh"
 
namespace CoupledField
{

Elecst2dPDE::Elecst2dPDE(AbstractAlgebraicSys * ptalgsys, Grid * aptgrid, Material * ptMaterial, TimeFunc * aptTimeFunc, FileType * aptFileType, WriteResults * aptOut)
:BasePDE(ptalgsys,ptMaterial,aptFileType,aptOut,aptTimeFunc)
{
#ifdef TRACE
  (*trace) << "entering Elecst2dPDE::Electst2dPDE " << std::endl;
#endif

  dofspernode_=1;
  ptgrid_=aptgrid;

  size_=ptgrid_->GetMaxnumnodes(0);
  sol_.Resize(size_);
  sol_.Init(0);

  conf->getsubdompde(subdoms_,"Elecst2d");
  ReadBCs("Elecst2d");
}

void Elecst2dPDE::SpecifySolver(Integer &solvertype, Integer &precondtype, Double &eps, Double &dampiter, Integer &maxnumit, Integer &numeqcoarse, Double &coarsealpha )
{
#ifdef TRACE
  (*trace) << "entering Elecst2dPDE::SpecifySolver" << std::endl;
#endif

  conf->get("eps",eps,"Elecst2d"); // relative accuracy in the precond. energy
  conf->get("dampiter",dampiter,"Elecst2d"); // damping parameter for Jacobi, SSOR
  conf->get("maxnumit",maxnumit,"Elecst2d"); // max number of iterations
  conf->get("solvertype",solvertype,"Elecst2d"); // Richardson or CG
  conf->get("precondtype", precondtype, "Elecst2d"); //ID or MG
  conf->get("numeqcoarse",numeqcoarse,"Elecst2d"); // number of equation for coarsing
  conf->get("coarsealpha",coarsealpha,"Elecst2d"); // coarsing parameter for AMG
}

void Elecst2dPDE::SetMatrixFactors()
{
#ifdef TRACE
  (*trace) << "entering Elecst2dPDE::SetMatrixFactors" << std::endl;
#endif
 
  matrix_factor_[0] = 1.0;
  matrix_factor_[1] = 0.0;
  matrix_factor_[2] = 0.0;
  matrix_factor_[3] = 0.0;
}

void Elecst2dPDE::SpecifyMatrices(Integer &matrixtype, Integer * matrixsystype, Integer &graphtype, Integer &numdofpernode, Integer &numdirichlets,
Integer &numconstraints)
{
#ifdef TRACE
  (*trace) << "entering Elecst2dPDE::SpecifyMatrices" << std::endl;
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
  matrixsystype[2] = 0;   // memory for the damping matrix
  matrixsystype[3] = 0;   // memory for the convection matrix
  matrixsystype[4] = 0;   // memory for the mass matrix

  graphtype = 1; // NOGRAPH = 0
                 // NODE   = 1
                 // EDGE   = 2
                 // FACE   = 3
                 // VOLUME = 4

  numdofpernode  = 1;
  numdirichlets = 1;
  numconstraints = 0;
}

void Elecst2dPDE::SetupMatrices(const Integer level)
{
#ifdef TRACE
  (*trace) << "entering Elecst2dPDE::SetupMatrices" << std::endl;
#endif
  
  Matrix<Double> elemmat;
  Point2D * ptCoord;

  BaseElem * ptElem;

  Integer matrix_stiff=2;

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

	  BaseForm<Point2D> * bilinear_stiff = new LaplaceInt<Point2D>(ptElem,1);

	  connecth=elemssd[j]->connect;

	  ptCoord=new Point2D[connecth.size()];
	  ptgrid_->GetCoordNodesElem(connecth,ptCoord,level);

	  // stiffness part
	  bilinear_stiff->CalcElemMatrix(ptCoord, elemmat);
	  elemmat*=coeffst[i];

#ifdef DEBUG
	  (*debug) << "Stiffnessmatrix, ElementNumber  " <<   j << std::endl;

	  (*debug) << elemmat << std::endl;
#endif

	  ptalgsys_->PutElemMatAlgSys(elemmat.getinarray(), connecth.get(), connecth.size(), as_sysid_, as_sysid_, matrix_stiff);

	  delete bilinear_stiff;
	  delete [] ptCoord;
	}     
    }

#ifdef TRACE
  (*trace) << "Leaving Elecst2dPDE::SetupMatrices" << std::endl;
#endif
}

void Elecst2dPDE::SetBCs(BCs * ptBCs, const Integer level, const Integer update, const Double atime)
{
#ifdef TRACE
  (*trace) << "entering Elecst2dPDE::SetBCs" << std::endl;
#endif

  Integer node;
  Double val, valueTF;

  //system matrix: id = 1
  Integer matrix_id = 1;

  val=ptTimeFunc_->TimeFuncAtTime(atime,level);    

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
              ptalgsys_->SetBCDirichletUpdate(j+1, node, val, dofspernode_, as_sysid_, as_sysid_, matrix_id);
            }
          else
            {
              ptalgsys_->SetBCDirichlet(j+1, node, val, dofspernode_, as_sysid_,
					as_sysid_, matrix_id);
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
              ptalgsys_->SetBCDirichletUpdate(j+1, node, val, dofspernode_, as_sysid_, as_sysid_, matrix_id);
            }
          else
            {
              ptalgsys_->SetBCDirichlet(j+1, node, val, dofspernode_, as_sysid_,as_sysid_, matrix_id);
            }
	}
    }
  
#ifdef TRACE
  (*trace) << " leaving Elecst2d::ReadBCs " << std::endl;
#endif
 
}

void Elecst2dPDE::SolveStepStatic(BCs * ptBCs, Integer level)
{
#ifdef TRACE
  (*trace) << "entering Elecst2dPDE::SolveStepStatic" << std::endl;
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

void Elecst2dPDE:: WriteResultsInFile()
{
#ifdef TRACE
  (*trace) << "entering Elecst2dPDE::WriteResultsInFile" << std::endl;
#endif

  Integer step=0;
  Double time=0;
  if (OutFile_->IsGMV())
    OutFile_->WriteSolution(sol_,step,time,"electric_potential");
  else
    OutFile_->WriteSolution(sol_,step,time,"electric potential");
}

Double Elecst2dPDE::CalcEnergyNorm()
{
  Double help1;
  help1=ptalgsys_->CalcEnergyNorm(0,0,2,sol_.get());
 
  return sqrt(help1);
}

void Elecst2dPDE::CalcCoeff(Vector<Double> & coeff)
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

SpaceErrorEstimator * Elecst2dPDE::CreatePtSpaceError()
{
  return ptSpaceError_=new Elecst2dErrEstimator(this,ptgrid_);
}

Elecst2dPDE::~Elecst2dPDE()
{
 ;
}

} // end of namespace
