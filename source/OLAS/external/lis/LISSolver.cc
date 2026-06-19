#include "LISSolver.hh"

#include <string>
#include "MatVec/SparseOLASMatrix.hh"
#include "MatVec/SCRS_Matrix.hh"
#include "MatVec/CRS_Matrix.hh"
#include "DataInOut/ProgramOptions.hh"
#include <sstream>
#include <stdio.h>

#ifdef _OPENMP
  #include <omp.h>
#endif

// Include lis.h AFTER all std/CFS headers so its complex-mode math macros are not
// substituted into them, then undef the macros so the rest of this translation unit uses normal math.
#include "lis_config.h"
#ifdef USE_COMPLEX
// selects LIS's complex scalar (double[2] in C++); 16 bytes
  #define _COMPLEX        
#endif
#include "lis.h"
#include "lis_precon.h"
#ifdef _COMPLEX
  #undef acos
  #undef acosh
  #undef asin
  #undef asinh
  #undef atan
  #undef atanh
  #undef cos
  #undef cosh
  #undef exp
  #undef fabs
  #undef log
  #undef pow
  #undef proj
  #undef sin
  #undef sinh
  #undef sqrt
  #undef tan
  #undef tanh
#endif

using std::string;
using boost::lexical_cast;

namespace CoupledField{

  static EnumTuple lisSolverTypeTuples[] =
  {
    EnumTuple( LISSolver::NOSOLVER, "noSolver" ),
    EnumTuple( LISSolver::CG, "CG" ),
    EnumTuple( LISSolver::BICG, "BiCG"),
    EnumTuple( LISSolver::CGS, "CGS"),
    EnumTuple( LISSolver::BICGSTAB, "BiCGSTAB" ),
    EnumTuple( LISSolver::BICGSTABL, "BiCGSTABL" ),
    EnumTuple( LISSolver::GPBICG, "CPBiCG"),
    EnumTuple( LISSolver::TFQMR, "TFQMR"),
    EnumTuple( LISSolver::ORTHOMIN, "Orthomin" ),
    EnumTuple( LISSolver::GMRES, "GMRES" ),
    EnumTuple( LISSolver::JACOBI, "Jacobi"),
    EnumTuple( LISSolver::GS, "GaussSeidel"),
    EnumTuple( LISSolver::SOR, "SOR" ),
    EnumTuple( LISSolver::BICGSAFE, "BiCGSave" ),
    EnumTuple( LISSolver::CR, "CR"),
    EnumTuple( LISSolver::BICR, "BiCR"),
    EnumTuple( LISSolver::CRS, "CRS"),
    EnumTuple( LISSolver::BICRSTAB, "BiCRSTAB"),
    EnumTuple( LISSolver::GPBICR, "GPBiCR"),
    EnumTuple( LISSolver::BICRSAFE, "BiCRSAFE"),
    EnumTuple( LISSolver::FGMRES, "FGMRES"),
    EnumTuple( LISSolver::IDRS, "IDRs"),
    EnumTuple( LISSolver::IDR1, "IDR1"),
    EnumTuple( LISSolver::MINRES, "minres"),
    EnumTuple( LISSolver::COCG, "COCG"),
    EnumTuple( LISSolver::COCR, "COCR")
  };

  Enum<LISSolver::LISSolverType> LISSolver::lisSolverType = \
  Enum<LISSolver::LISSolverType>("LIS Solver Types",
      sizeof(lisSolverTypeTuples) / sizeof(EnumTuple),
      lisSolverTypeTuples);

  static EnumTuple lisPrecondTypeTuples[] =
  {
    EnumTuple( LISSolver::NONE, "none" ),
    EnumTuple( LISSolver::JACOBI_PRE, "jacobi" ),
    EnumTuple( LISSolver::BJACOBI, "block_jacobi" ),
    EnumTuple( LISSolver::ILU, "iluk"),
    EnumTuple( LISSolver::SSOR, "ssor"),
    EnumTuple( LISSolver::HYBRID, "hybrid" ),
    EnumTuple( LISSolver::IS, "is" ),
    EnumTuple( LISSolver::SAINV, "sainv"),
    EnumTuple( LISSolver::SAAMG, "saamg"),
    EnumTuple( LISSolver::ILUC, "iluc" ),
    EnumTuple( LISSolver::ILUT, "ilut" )
  };

  Enum<LISSolver::LISPrecondType> LISSolver::lisPrecondType = \
  Enum<LISSolver::LISPrecondType>("LIS Preconditioner Types",
      sizeof(lisPrecondTypeTuples) / sizeof(EnumTuple),
      lisPrecondTypeTuples);


LISSolver::LISSolver(PtrParamNode pn, PtrParamNode olasInfo, BaseMatrix::EntryType type){

  int argc = -1;                   /* dummy arg count */
  char **argv = NULL;              /* dummy arg */
  Integer err = 0;
  err = lis_initialize(&argc,&argv); CHKERR(err);

  infoNode_ =  olasInfo->Get("lis");

  xml_ = pn;
  firstSetup_ = true;
  ownMatrixA_ = false;

  maxIter_    = pn->Has("maxIter") ? pn->Get("maxIter")->As<int>() : 10000;
  tolerance_  = pn->Has("tolerance") ? pn->Get("tolerance")->As<double>() : 1e-12;
  minTol_     = pn->Has("minimalTolerance") ? pn->Get("minimalTolerance")->As<double>() : 1e-11;
  logging_    = pn->Has("logging") ? pn->Get("logging")->As<bool>() : false;
  resetXZero_ = pn->Has("zeroInitialValue") ? pn->Get("zeroInitialValue")->As<bool>() : false;

  PtrParamNode hdr = infoNode_->Get(ParamNode::HEADER);
  hdr->Get("maxIter")->SetValue(maxIter_);
  hdr->Get("tolerance")->SetValue(tolerance_);
  hdr->Get("minimalTolerance")->SetValue(minTol_);
  if(minTol_ < tolerance_)
    hdr->SetWarning("minimal tolerance " + lexical_cast<string>(minTol_) + " ignored as below tolerance " + lexical_cast<string>(tolerance_));

  // Solve() also sets solver and precond as attributes to xml_ but w/o arguments and also when the elements here are not given
  if(xml_->Has("solver"))
    hdr->Get("solver")->SetValue(xml_->Get("solver"), false);
  if(xml_->Has("precond"))
    hdr->Get("precond")->SetValue(xml_->Get("precond"),false);
}

LISSolver::~LISSolver(){
  Integer err = 0;
  err = lis_vector_destroy(x_); CHKERR(err);
  err = lis_vector_destroy(b_); CHKERR(err);

  //matrix A_ shares a pointer in case of real valued problems
  //there would be a double free otherwise
  if(ownMatrixA_) {
    err = lis_matrix_destroy(A_); CHKERR(err);
  }
  err = lis_matrix_destroy(A0_); CHKERR(err);

  err = lis_precon_destroy(precond_);CHKERR(err);
  err = lis_solver_destroy(solver_);CHKERR(err);

  lis_finalize();
}

namespace {
  // Zero a LIS vector without a by-value LIS_SCALAR. In C++ LIS_SCALAR is the double[2]
  // array form, so the scalar overload won't bind and would be ABI-incompatible with the
  // C library regardless. Writing zero bytes is layout/ABI-neutral and yields a complex
  // zero the library reads correctly.
  inline void lisZero(LIS_VECTOR v){
    LIS_INT ln = 0;
    LIS_INT gn = 0;
    lis_vector_get_size(v, &ln, &gn);
    std::memset(v->value, 0, sizeof(LIS_SCALAR) * static_cast<size_t>(ln));
  }
}

void LISSolver::Setup(BaseMatrix &sysmat){

  const StdMatrix& stdmat = dynamic_cast<const StdMatrix&>(sysmat);
  
  BaseMatrix::EntryType etype = stdmat.GetEntryType();
  BaseMatrix::StorageType stype = stdmat.GetStorageType();

  Integer err=0;

  // we create the matrix...
  // Somehow after the first call the matrix assembled by CFS isn't updated anymore if we don't call lis_matrix_create again.
  // However w/o destroy we create memory leaks

  if(firstSetup_) {
    err = lis_matrix_create(0,&A_); CHKERR(err);
  } else {
    err = lis_matrix_destroy(A0_); CHKERR(err);
  }

  if(stype == BaseMatrix::SPARSE_SYM)
  {
    // TODO first validate and second, as we anyway work with A0_ we can convert the matrix A0_ only
    EXCEPTION("LIS solver cannot yet handle SCRS matrices. Please set sparseNonSym as storage type");
  }
  
  if(etype == BaseMatrix::DOUBLE)
  {
    // real case. The complex LIS build expects LIS_SCALAR (complex), so embed the real
    // values as complex with zero imaginary part. cbuf_ (member) persists across the solve.
    const CRS_Matrix<Double>& crs = dynamic_cast<const CRS_Matrix<Double>&>(stdmat);

    if(crs.GetNumCols() != crs.GetNumRows())
      EXCEPTION("LIS solver only tested for quadratic matrices");

    UInt nnz = crs.GetNnz();
    UInt dim = crs.GetNumRows();

    Integer * rowPtr = (Integer *)crs.GetRowPointer();
    Integer * colPtr = (Integer *)crs.GetColPointer();
    const Double * dataPtr = crs.GetDataPointer();

    cbuf_.assign(dataPtr, dataPtr + nnz);   // Double -> Complex(re, 0)

    err = lis_matrix_set_size(A_,dim,0); CHKERR(err);
    err = lis_matrix_set_csr(nnz,rowPtr,colPtr, reinterpret_cast<LIS_SCALAR*>(cbuf_.data()),A_); CHKERR(err);
    err = lis_matrix_assemble(A_); CHKERR(err);

    if(firstSetup_){
      err = lis_vector_duplicate(A_,&b_); CHKERR(err);
      lisZero(b_);
    }
    ownMatrixA_ = false;   // cbuf_ owns the data
  }
  else
  {
    // native complex symmetric case (LIS built with --enable-complex).
    // Complex and LIS_SCALAR share {re,im} layout -> share CFS' CRS arrays zero-copy.
    const CRS_Matrix<Complex>& crs = dynamic_cast<const CRS_Matrix<Complex>&>(stdmat);

    if(crs.GetNumCols() != crs.GetNumRows())
      EXCEPTION("LIS solver only tested for quadratic matrices");

    UInt nnz = crs.GetNnz();
    UInt dim = crs.GetNumRows();

    Integer * rowPtr = (Integer *)crs.GetRowPointer();
    Integer * colPtr = (Integer *)crs.GetColPointer();
    Complex * dataPtr = const_cast<Complex*>(crs.GetDataPointer());

    err = lis_matrix_set_size(A_,dim,0); CHKERR(err);
    err = lis_matrix_set_csr(nnz,rowPtr,colPtr, reinterpret_cast<LIS_SCALAR*>(dataPtr),A_); CHKERR(err);
    err = lis_matrix_assemble(A_); CHKERR(err);

    if(firstSetup_){
      err = lis_vector_duplicate(A_,&b_); CHKERR(err);
      lisZero(b_);
    }
    ownMatrixA_ = false;
  }

  if(firstSetup_){
    err = lis_vector_duplicate(b_,&x_); CHKERR(err);
    lisZero(x_);
  }
  if(resetXZero_ || firstSetup_){
    lisZero(x_);
  }

  //copy matrix (needed as a workaround for multiple iterations with different system matrices to solve without memory leak)
  err = lis_matrix_duplicate(A_,&A0_); CHKERR(err);
  lis_matrix_set_type(A0_,LIS_MATRIX_CSR);
  err = lis_matrix_convert(A_,A0_); CHKERR(err);

  //create the solver
  string config;
  if(firstSetup_ ){//|| solver_->A != A0_){
    CreateConfigString(xml_,config);
    err = lis_solver_create(&solver_); CHKERR(err);
    err = lis_solver_set_option(const_cast<char*>(config.c_str()),solver_);CHKERR(err);
  } else {
    err = lis_precon_destroy(precond_); CHKERR(err);
  }
  solver_->A = A0_;

  err = lis_precon_create(solver_, &precond_);
  CHKERR(err);
  if( err ){
    std::cerr << "There was an error creating the preconditioner. Code: " << err << " ...Going to abort" << std::endl;
    EXCEPTION("ERROR");
  }
  firstSetup_ = false;
}

void LISSolver::Solve( const BaseMatrix &sysmat, const BaseVector &rhs, BaseVector &sol)
{
  // LIS vectors hold LIS_SCALAR (complex); Complex shares the same {re,im} layout.
  const bool isReal = (sysmat.GetEntryType() == BaseMatrix::DOUBLE);

  Complex* bval = reinterpret_cast<Complex*>(b_->value);
  
  for(UInt i=0, n=rhs.GetSize(); i<n; i++){
    if(isReal){ 
      Double r = 0;
      rhs.GetEntry(i, r);
      bval[i] = Complex(r, 0.0);
    } else {
      Complex e = 0;
      rhs.GetEntry(i, e);
      bval[i] = e;
    }
  }

  Integer err = lis_solve_kernel(A_, b_, x_, solver_, precond_);
  if(err) EXCEPTION("Solver returned error code: " << err << " ...aborting");

  Complex* xval = reinterpret_cast<Complex*>(x_->value);
  for(UInt i=0, n=sol.GetSize(); i<n; i++){
    if(isReal) {
      sol.SetEntry(i, xval[i].real());
    } else {
      sol.SetEntry(i, xval[i]);
    }
  }

  // the general stuff
  int solverCode;
  lis_solver_get_solver(solver_,&solverCode);
  infoNode_->Get("solver")->SetValue(lisSolverType.ToString((LISSolverType) solverCode));

  int precondCode;
  lis_solver_get_precon(solver_,&precondCode);
  infoNode_->Get("precond")->SetValue(lisPrecondType.ToString((LISPrecondType) precondCode));

  //ParamNode::ActionType at = progOpts->DoDetailedInfo() ? ParamNode::APPEND : ParamNode::DEFAULT;
  PtrParamNode curr = infoNode_->Get(ParamNode::PROCESS)->Get("solve", ParamNode::APPEND); // for an iterative solve each solution should be interesting

  // this is only the solver time, in some rare cases the preconditioner can be much more expensive
  double lastTime = 0.0;
  lis_solver_get_time(solver_,&lastTime);
  curr->Get("timing")->SetValue(lastTime);

  double norm = 0.0;
  lis_solver_get_residualnorm(solver_,&norm);
  curr->Get("residualNorm")->SetValue(norm);

  int iterations = 0;
  lis_solver_get_iter(solver_,&iterations);
  curr->Get("iterations")->SetValue(iterations);

  if(norm > tolerance_ && norm <= minTol_)
    infoNode_->Get(ParamNode::SUMMARY)->SetWarning("residual norm " + lexical_cast<string>(norm) + " exceeds target "
          + lexical_cast<string>(tolerance_) + " but within minimal tolerance " + lexical_cast<string>(minTol_) + " after " +  lexical_cast<string>(norm) + " iterations");

  if(norm > tolerance_ && norm > minTol_ )
    EXCEPTION("after " << iterations << " iterations reached residual " << norm << " with target " << tolerance_ << " exceeding minminal tolerance " << minTol_); // CFS.cc will add it to info.xml

# // case minTol < tolerance give warning in consructor
}


void LISSolver::CreateConfigString(PtrParamNode configNode, string& output){
  string solStr;
  string precondStr;
  CreateSolverString(configNode,solStr);
  CreatePrecondString(configNode,precondStr);

  std::stringstream globStream;

  //lets read the basics
  globStream << " -maxiter " << maxIter_;

  globStream << " -tol " << tolerance_;

  globStream << " -print " << (logging_ ? "out" : "none");

  output = solStr + " " + precondStr + " " + globStream.str() + " -initx_ones false -initx_zeros false";
  infoNode_->Get("config")->SetValue(output);
  return;
}


void LISSolver::CreateSolverString(PtrParamNode solverNode, string& output){
  string solverString;
  string nodeName = solverNode->GetName();
  PtrParamNode sNode = solverNode->Get("solver",ParamNode::PASS);

  std::stringstream solstream;

  if(sNode){
    ParamNodeList sol = sNode->GetChildren();

    UInt numChilds = sol.GetSize();
    //new know that first child is id the other is the actual solver tag...
    //well, not very filsafe but good for now.
    if(numChilds==2){
      solverString = sol[1]->GetName();
      LISSolverType type;
      type = lisSolverType.Parse(solverString);
      UInt typeID = type;

      solstream << " -i "<< typeID;

      PtrParamNode curNode;
      UInt deg = 0;
      UInt restart = 0;
      Double omega = 0;
      switch(type){
      case BICGSTABL:
        deg = 0;
        sol[1]->GetValue("degree",deg,ParamNode::PASS);
        solstream << " -ell "<< deg;
        break;
      case ORTHOMIN:
        restart = 0;
        sol[1]->GetValue("restart",restart,ParamNode::PASS);
        solstream << " -restart "<< restart;
        break;
      case GMRES:
        restart = 0;
        sol[1]->GetValue("restart",restart,ParamNode::PASS);
        solstream << " -restart "<< restart;
        break;
      case SOR:
        omega = 0;
        sol[1]->GetValue("omega",omega,ParamNode::PASS);
        solstream << " -omega "<< omega;
        break;
      case FGMRES:
        restart = 0;
        sol[1]->GetValue("restart",restart,ParamNode::PASS);
        solstream << " -restart "<< restart;
        break;
      case IDRS:
        restart = 0;
        sol[1]->GetValue("irestart",restart,ParamNode::PASS);
        solstream << " -irestart "<< restart;
        break;
      case CG:
      case BICG:
      case CGS:
      case BICGSTAB:
      case GPBICG:
      case TFQMR:
      case JACOBI:
      case GS:
      case BICGSAFE:
      case CR:
      case BICR:
      case CRS:
      case BICRSTAB:
      case GPBICR:
      case BICRSAFE:
      case MINRES:
      case IDR1:
      case COCG:
      case COCR:
        //nothing to do
        break;
      default:
        EXCEPTION("Could not obtain a valid solver");
        break;
      }

    }else if(numChilds > 2){
      EXCEPTION("Expected 2 Childs of element solver got " << numChilds << "!");
    }
  }

  output = solstream.str();
}


void LISSolver::CreatePrecondString(PtrParamNode precondNode, string& output){
  string precondString;

  PtrParamNode sNode = precondNode->Get("precond",ParamNode::PASS);

  std::stringstream solstream;

  if(sNode){
    ParamNodeList sol = sNode->GetChildren();

    UInt numChilds = sol.GetSize();
    //new know that first child is id the other is the actual solver tag...
    //well, not very filsafe but good for now.
    if(numChilds==2){
      precondString = sol[1]->GetName();
      LISPrecondType type;
      type = lisPrecondType.Parse(precondString);
      UInt typeID = type;

      solstream << " -p "<< typeID;

      PtrParamNode curNode;
      Double fill = 0;
      Double prec_omega = 1;
      Double alpha = 1.0;
      UInt m = 3;
      Double drop = 0.05;
      bool unsym = 0;
      Double theta = 0;
      Double fillInRate = 5.0;
      switch(type){
      case ILU:
        fill = 0;
        sol[1]->GetValue("fill",fill,ParamNode::PASS);
        solstream << " -ilu_fill "<< fill;
        break;
      case SSOR:;
        prec_omega = 1;
        sol[1]->GetValue("omega",prec_omega,ParamNode::PASS);
        solstream << " -ssor_w "<< prec_omega;
        break;
      case HYBRID:
        EXCEPTION("Hybrid preconditioner not supported right now...")
        break;
      case IS:
        alpha = 1.0;
        sol[1]->GetValue("alpha",alpha,ParamNode::PASS);
        solstream << " -is_alpha "<< alpha;
        m = 3;
        sol[1]->GetValue("m",m,ParamNode::PASS);
        solstream << " -is_m "<< alpha;
        break;
      case SAINV:
        drop = 0;
        sol[1]->GetValue("drop",drop,ParamNode::PASS);
        solstream << " -sainv_drop "<< drop;
        break;
      case SAAMG:
        unsym = 0;
        theta = 0;
        sol[1]->GetValue("unsym",unsym,ParamNode::PASS);
        sol[1]->GetValue("theta",theta,ParamNode::PASS);
        solstream << " -saamg_unsym "<< unsym;
        solstream << " -saamg_theta "<< theta;
        break;
      case ILUC:
        drop = 0.05;
        fillInRate = 5.0;
        sol[1]->GetValue("drop",drop,ParamNode::PASS);
        sol[1]->GetValue("fillInRate",fillInRate,ParamNode::PASS);
        solstream << " -iluc_drop "<< drop;
        solstream << " -iluc_rate "<< fillInRate;
        break;
      case ILUT:
        drop = 0.05;
        fillInRate = 5.0;
        sol[1]->GetValue("drop",drop,ParamNode::PASS);
        sol[1]->GetValue("fillInRate",fillInRate,ParamNode::PASS);
        solstream << " -ilut_drop "<< drop;
        solstream << " -ilut_rate "<< fillInRate;
        break;
      case NONE:
      case JACOBI_PRE:
      case BJACOBI:
        //nothing to do
        break;
      default:
        EXCEPTION("Could not obtain a valid preconditioner");
        break;
      }

    }else if(numChilds > 2){
      EXCEPTION("Expected 2 Childs of element solver got " << numChilds << "!");
    }
  }
  output = solstream.str();

}

}
