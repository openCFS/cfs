#include <fstream>
#include <iostream>
#include <string>
#include <math.h>

#include "elecst3dPDE.hh"
#include "outUnverg.hh"
#include "forms_header.hh"
 
namespace CoupledField
{

Elecst3dPDE::Elecst3dPDE(AbstractAlgebraicSys * ptalgsys, Grid * aptgrid, Material * ptMaterial, TimeFunc * aptTimeFunc, FileType * aptFileType, WriteResults * aptOut)
:BasePDE(ptalgsys,ptMaterial,aptFileType,aptOut,aptTimeFunc)
{
#ifdef TRACE
  (*trace) << "entering Elecst3dPDE::Electst3dPDE " << std::endl;
#endif

//  doftype_=4;
  doftype_=ep_restraint;
  dofspernode_=1;
  ptgrid_=aptgrid;

  size_=ptgrid_->GetMaxnumnodes(0);
  sol_.Resize(size_);
  sol_.Init(0);
}

void Elecst3dPDE::SpecifySolver(Integer &solvertype, Integer &precondtype, Double &eps, Double &dampiter, Integer &maxnumit, Integer &numeqcoarse)
{
#ifdef TRACE
  (*trace) << "entering Elecst3dPDE::SpecifySolver" << std::endl;
#endif

  conf->get("eps",eps,"Elecst3d"); // relative accuracy in the precond. energy
  conf->get("dampiter",dampiter,"Elecst3d"); // damping parameter for Jacobi, SSOR
  conf->get("maxnumit",maxnumit,"Elecst3d"); // max number of iterations
  conf->get("solvertype",solvertype,"Elecst3d"); // Richardson or CG
  conf->get("precondtype", precondtype, "Elecst3d"); //ID or MG
  conf->get("numeqcoarse",numeqcoarse,"Elecst3d"); // number of equation for coarsing
}

void Elecst3dPDE::SetMatrixFactors()
{
#ifdef TRACE
  (*trace) << "entering Elecst3dPDE::SetMatrixFactors" << std::endl;
#endif
 
  matrix_factor_[0] = 1.0;
  matrix_factor_[1] = 0.0;
  matrix_factor_[2] = 0.0;
  matrix_factor_[3] = 0.0;
}

void Elecst3dPDE::SpecifyMatrices(Integer &matrixtype, Integer * matrixsystype, Integer &graphtype, Integer &numdofpernode, Integer &numdirichlets,
Integer &numconstraints)
{
#ifdef TRACE
  (*trace) << "entering Elecst3dPDE::SpecifyMatrices" << std::endl;
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

void Elecst3dPDE::SetupMatrices(const Integer type)
{
#ifdef TRACE
  (*trace) << "entering Elecst3dPDE::SetupMatrices" << std::endl;
#endif
  
  Integer level=0;

  Matrix<Double> elemmat;
  Point3D * ptCoord;

  Integer numelems=ptgrid_->GetMaxnumElem(0);

  BaseElem * ptElem;

  Integer matrix_stiff=2;

  Double coeff;
  CalcCoeff(coeff,0);  //// Attention !!!

  Vector<Integer> connecth;

  Integer i;
  for (i=0; i<numelems; i++)
{  
  ptElem=ptgrid_->GetptElem(i);

 BaseForm<Point3D> * bilinear_stiff = new LaplaceInt<Point3D>(ptElem,1);

   ptgrid_->GetConnection(connecth,i,level);

     ptCoord=new Point3D[connecth.size()];
     ptgrid_->GetCoordOfNodesElem(i,0,connecth.size(),ptCoord);

      // stiffness part
      bilinear_stiff->CalcElemMatrix(ptCoord, elemmat);
      elemmat*=coeff;

#ifdef DEBUG
      (*debug) << "Stiffnessmatrix, ElementNumber  " <<   i << std::endl;

      (*debug) << elemmat << std::endl;
#endif

      ptalgsys_->PutElemMatAlgSys(elemmat.getinarray(), connecth.get(), connecth.size(), AS_sysid_, AS_sysid_, matrix_stiff);

      delete bilinear_stiff;
      delete [] ptCoord;
     
}

#ifdef TRACE
  (*trace) << "Leaving Elecst3dPDE::SetupMatrices" << std::endl;
#endif
}

void Elecst3dPDE::SetBCs(BCs * ptBCs, const Integer level, const Integer update, const Double atime)
{
#ifdef TRACE
  (*trace) << "entering Elecst3dPDE::SetBCs" << std::endl;
#endif

  Integer num, node, type, tfunc;
  Double val, valueTF;


  //system matrix: id = 1
  Integer matrix_id = 1;

    std::list<NodeRestraint> restr;
    ptBCs->GetRestraints(restr,level);

    val=ptTimeFunc_->TimeFuncAtTime(atime,level);    

    Integer i=0;
    for (list<NodeRestraint>::const_iterator p=restr.begin(); p!=restr.end(); p++, i++)
   {
          Integer node=p->nodalnum;
 
          if (p->dof==doftype_)
	    {
          if (update==1)
            {
              ptalgsys_->SetBCDirichletUpdate(i+1, node, val, dofspernode_, AS_sysid_, AS_sysid_, matrix_id);
            }
          else
            {
              ptalgsys_->SetBCDirichlet(i+1, node, val, dofspernode_, AS_sysid_,
AS_sysid_, matrix_id);
            }
	    }
    }
}

void Elecst3dPDE::SolveStepStatic(BCs * ptBCs, Integer level)
{
#ifdef TRACE
  (*trace) << "entering Elecst3dPDE::SolveStepStatic" << std::endl;
#endif

  Integer type = 0;
  Integer update = 0;
  Integer job = 1;

  Double * ptsol;

  SetupMatrices(type);
  ptalgsys_->ComputeSysMatrix(AS_sysid_,AS_sysid_,matrix_factor_);
  SetBCs(ptBCs,level,update,0);
  ptalgsys_->ComputePrecond(job,AS_sysid_);
  ptalgsys_->SolveAlgSys(AS_sysid_);

  ptsol = ptalgsys_->GetSolution(AS_sysid_);

    // save solution
  Vector<Double> transsol(ptgrid_->GetMaxnumnodes(level), ptsol);
  sol_=transsol;
}

void Elecst3dPDE:: WriteResultsInFile()
{
#ifdef TRACE
  (*trace) << "entering Elecst3dPDE::WriteResultsInFile" << std::endl;
#endif

  Integer step=0;
  Double time=0;
if (OutFile_->IsGMV())
 OutFile_->WriteSolution(sol_,step,time,"electric_potential");
else
  OutFile_->WriteSolution(sol_,step,time,"electric potential");

}

Double Elecst3dPDE::CalcEnergyNorm()
{
 Double help1;
 help1=ptalgsys_->CalcEnergyNorm(0,0,2,sol_.get());
 
 return sqrt(help1);
}

void Elecst3dPDE::CalcCoeff(Double & coeff, const Integer numsubdom)
{
 // get material number for subdomain with number numsubdom from config-file
 Integer matnum;
// conf->getmatnum(matnum,"elecst3d");
 conf->get("elecst3d",matnum);

 // read density and compress with material number matnum
 Double dielectr;
 if (!MatFile_) Error("You didn't specialize material file. Use option -m");
 MatFile_->ReadDielectricTerms(dielectr,matnum); 

 coeff=dielectr;
}

Elecst3dPDE::~Elecst3dPDE()
{
 ;
}

} // end of namespace
