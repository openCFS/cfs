#include <iostream>
#include <string>
#include <algorithm>

#include <def_use_pardiso.hh>
#include <def_use_metis.hh>

#include "General/Exception.hh"
#include "MatVec/SparseOLASMatrix.hh"
#include "MatVec/SCRS_Matrix.hh"
#include "MatVec/CRS_Matrix.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ProgramOptions.hh"
#include "Utils/StdVector.hh"
#include "Ilupack.hh"

extern "C"
{
  #include "reloj.h"
  #include "InputOutput.h"
  #include "ScalarVectors.h"
  #include "SparseVectorsNew.h"
  #include "SparseSymmetricNew.h"
  #include "SparseHarwellBoeingNew.h"
  #include "EliminationTree.h"
  #include "TaskQueue.h"
  #include "ToolsIlupack.h"
  #include "ToolsOPENMP.h"
  #include "SPDfactorOPENMP.h"
  #include "SPDsolverOPENMP.h"

}

#ifdef _OPENMP
  #include <omp.h>
#endif


DECLARE_LOG(ilupack)
DEFINE_LOG(ilupack, "ilupack")

namespace CoupledField
{

template<typename T>
Ilupack<T>::Ilupack(PtrParamNode xml, PtrParamNode olasInfo, BaseMatrix::EntryType type)
{



  // we work with out 
  xml_ = xml;
  infoNode_ = olasInfo->Get("ilupack");
  
  if (type != BaseMatrix::COMPLEX && type != BaseMatrix::DOUBLE)EXCEPTION("unhandled type " << type);
  isComplex_ = type == BaseMatrix::COMPLEX;

  LOG_TRACE(ilupack) <<  "Ilupack(): isComplex=" << isComplex_;

  // set this to zero to signal the destructor if Setup() was called
  mat.a = NULL;
  mat.ia = NULL;
  mat.ja = NULL;

  // init the enums
  SetEnums();


}

template<typename T>
Ilupack<T>::~Ilupack()
{
  LOG_TRACE(ilupack) <<  "~Ilupack()";
  if(isParallel){
    RemoveSparseMatrix (&spr);
  }
  else{
    IlupackAMGDelete();
  }
}

template<typename T>
void Ilupack<T>::SetMatrix(const BaseMatrix &base_mat)
{
  const SparseOLASMatrix<T>& som =  dynamic_cast<const SparseOLASMatrix<T>&> (base_mat);

  LOG_TRACE2(ilupack) <<  "SetMatrix: SPARSE_SYM="  << (som.GetStorageType() == BaseMatrix::SPARSE_SYM)
                      << " mat_.a=" << mat.a << " .ia=" << mat.ia << ".ja=" << mat.ja;

  // delete the old values before we allocate the new space - if could be faster! 
  if(mat.a != NULL)
  { delete[] mat.a; mat.a = NULL;}
  if(mat.ia != NULL)
  { delete[] mat.ia; mat.ia = NULL;}
  if(mat.ja != NULL)
  { delete[] mat.ja; mat.ja = NULL;}

  unsigned int elements = 0;

  // pointers will point to the first actual element    
  const unsigned int* row_ptr = NULL;
  const unsigned int* col_ptr = NULL;
  const T* val_ptr = NULL;

  if(som.GetStorageType() == BaseMatrix::SPARSE_SYM)
  {
    const SCRS_Matrix<T>& scrs = dynamic_cast<const SCRS_Matrix<T>&>(som);
    // Nnz is the total number = Diagonal + 2 * triagonals. But one is skipped!
    elements= scrs.GetNnz() - (unsigned int) (0.5 * (double) (scrs.GetNnz() - scrs.GetNumRows()));

    row_ptr = scrs.GetRowPointer();
    col_ptr = scrs.GetColPointer();
    val_ptr = scrs.GetDataPointer();
  }
  else
  {
    const CRS_Matrix<T>& crs = dynamic_cast<const CRS_Matrix<T>&>(som);

    elements= crs.GetNnz();

    row_ptr = crs.GetRowPointer();
    col_ptr = crs.GetColPointer();
    val_ptr = crs.GetDataPointer();
  }

  // allocate and copy the stuff - note, that we are 0-based in OLAS
  // but 1-based in Ilupack!

  // rows are simple but we have to handle the tailing element
  mat.ia = new int[som.GetNumRows() + 1]; // one plus due to CRS tail
  std::copy(row_ptr, row_ptr + som.GetNumRows() + 1, mat.ia);

  mat.ja = new int[elements];
  std::copy(col_ptr, col_ptr + elements, mat.ja);

  // values is brutal stuff, because binary the std::complex<double>* IS a double*!
  mat.a = reinterpret_cast<double*>(new T[elements]);
  T* val_src = reinterpret_cast<T*>(mat.a);
  std::copy(val_ptr, val_ptr + elements, val_src);

  // no adjust to 1-based!
  for(unsigned int i = 0, numRows = som.GetNumRows(); i < numRows+1; i++)
    mat.ia[i] += 1;

  for(unsigned int i = 0; i < elements; i++)
    mat.ja[i] += 1;

  mat.nr = som.GetNumRows();
  mat.nc = som.GetNumCols();
  mat.nnz = elements;

  mat.issymmetric = 0;
  mat.isdefinite = 0;
  mat.ishermitian = 0;
  mat.isskew = 0;
  mat.isreal = !isComplex_;
  mat.issingle = 0;
  switch(matrix_)
  {
    case GNL:
      break;
    case PD:
      mat.issymmetric = 1;
      mat.isdefinite = 1;
      if (isComplex_)
      {
        mat.ishermitian = 1;
      }
      break;
    case SYM:
      mat.issymmetric = 1;
      break;
    case HER:
      mat.ishermitian = 1;
      if (mat.isreal)
      {
        throw Exception("ilupack matrix is set to hermitian but not complex");
      }
      break;
    default:
      throw Exception("matrix type does not exist (SSM and SHR are not implemented yet)");
  }

  LOG_TRACE2(ilupack) << "SetMatrix: allocate: mat_.a<T>=" << elements << " .ia=" << (som.GetNumRows() + 1)
                      << ".ia=" << elements << "; .nr=" << mat.nr << " .nc=" << mat.nc;
  LOG_DBG2(ilupack) << "mat_.ia: " << StdVector<int>::ToString(mat.nr + 1, mat.ia);
  LOG_DBG2(ilupack) << "mat_.ja: " << StdVector<int>::ToString(elements, mat.ja);
  LOG_DBG2(ilupack) << "mat_.a: " << StdVector<double>::ToString(elements, mat.a);

}

template<typename T>
void Ilupack<T>::Setup(BaseMatrix &sysMat)
{
  // do we really want to create a new entry? Might blast up the output
  ParamNode::ActionType at = progOpts->DoDetailedInfo() ? ParamNode::APPEND : ParamNode::DEFAULT;
  PtrParamNode out = infoNode_->Get(ParamNode::PROCESS)->Get("setup", at);
  // determine the matrix type. Symmetric/nonsymmetric, positive definite, ...
  // it is optional given in the xml file.

  // GNL, SYM, PD, HER, ....
  DetermineMatrixType(sysMat, out);

  LOG_TRACE2(ilupack) <<  "Setup: matrix -> " << matrix.ToString(matrix_);

  // in case we already run release memory - it's save if first run
//  IlupackAMGDelete();

  // set the ilupack matrix mat from the olas matrix - complex can be casted to double
  SetMatrix(sysMat);

  // gain and modify the default settings, does logging uses mat and sets param 
  InitParameters();
}


template<typename T>
void Ilupack<T>::Solve(const BaseMatrix &base_mat, 
    const BaseVector &base_rhs,  BaseVector &base_sol)
{

  ParamNode::ActionType at = progOpts->DoDetailedInfo() ? ParamNode::APPEND : ParamNode::DEFAULT;
  PtrParamNode out = infoNode_->Get(ParamNode::PROCESS)->Get("solver", at);

  auto parameter= reinterpret_cast<DILUPACKparam*>(&param);// for OMP template is not supported in the lib. So casting it manually to
  int ierr;
  ptr_IlupackFactor  vFact=nullptr ;

  isParallel?ierr=IlupackFactorizationOMP(spr,index,*parameter,nleaves,mtmetis,&vFact):ierr=IlupackAMGFactor();

  // factorize the iLU preconditioner
  std::stringstream ss;
  ss << "Error factorizing Ilupack: ";

  switch (ierr)
  {
    case 0:
    // perfect:
    break;

    case -1:
    ss << "Input matrix may be wrong at level " << precond.nlev;
    break;

    case -2:
    ss << "Out of memory. The matrix L overflows the array alu at level " << precond.nlev;
    break;

    case -3:
    ss << "Out of memory. The matrix U overflows the array alu at level " << precond.nlev;
    break;

    case -4:
    ss << "Illegal value for lfil at level " << precond.nlev;
    break;

    case -5:
    ss << "Zero row encountered at level " <<  precond.nlev;
    break;

    case -6:
    ss << "Zero column encountered at level " <<  precond.nlev;
    break;

    case -7:
    ss << "Buffers are too small";
    break;

    default:
      ss << "Zero pivot encountered at step number " << ierr << " of level " << precond.nlev;
  }

  if(ierr != 0) throw Exception(ss.str());

  if (!isParallel){

    out->Get("levels")->SetValue(precond.nlev);
    CalcFillIn(out);
    PtrParamNode timing = out->Get("timing");
    timing->Get("total_time")->SetValue(ILUPACK_secnds[7]);
    timing->Get("initial_preprocessing")->SetValue(ILUPACK_secnds[0]);
    timing->Get("reordering_remaining_levels")->SetValue(ILUPACK_secnds[1]);
  }





  // the preconditioner sets the ilupack matrix
  if (mat.a == NULL)
    throw Exception("Setup() not called before Solve()");



  if (isParallel){
    double *rhs, *sol;
    // Create the vectors
    CreateDoubles (&rhs, spr.dim1); CreateDoubles (&sol, spr.dim1);

    sol = dynamic_cast<Vector<double>&> (base_sol).GetPointer();
    rhs =  const_cast<double*> (dynamic_cast<const Vector<double>&> (base_rhs).GetPointer());
    InitDoubles (sol, spr.dim1, 0.01, 0.001);
    InitRandDoubles (sol, spr.dim1, 0.0, 1.0);
    IlupackSolverOMPG(spr, index, rhs, sol, vFact, *parameter, &parameter->maxit,&parameter->restol);
  }
  else{
    // we gain the solution pointer
     T* sol_ptr = dynamic_cast<Vector<T>&> (base_sol).GetPointer();
     // additionally remove the constness of the rhs, ilupack knows no const.
     T* rhs_ptr =  const_cast<T*> (dynamic_cast<const Vector<T>&> (base_rhs).GetPointer());
     LOG_TRACE2(ilupack) <<  "Solve: sol_ptr=" << sol_ptr << " rhs_ptr=" << rhs_ptr;
     ierr = IlupackAMGSolver(sol_ptr, rhs_ptr);
  }

  ss << "Error Solving Ilupack ";

  // why did the iterative solver stop?
  switch (ierr)
  {
    case  0:  // everything is fine
      break;

    case -1:  // too many iterations
      ss << "Maximum number of iteration steps has been exceeded.";
      break;

    case -2:
      ss << "Not enough work space provided.";
      break;

    case -3:  /* not enough work space */
      ss << "Algorithm breaks down.";
      break;

    default:
      ss << "Solver exited with error code: " << ierr << ".";
  }

  // why did the iterative solver stop?
  switch (ierr)
  {
    case  0:  // everything is fine
      break;

    case -1:  // too many iterations
      WARN(ss.str());
      break;

    default:
      // stop if necessary
      throw Exception(ss.str());
  }

  out->Get("iterations")->SetValue(param.ipar[26]);
  PtrParamNode timing = out->Get("timing");
  timing->Get("total")->SetValue(ILUPACK_secnds[5]);
  timing->Get("maxtrix_vector_mult")->SetValue(ILUPACK_secnds[6]);
  PtrParamNode norms = out->Get("norms");
  norms->Get("target")->SetValue(param.fpar[23]);


}



template<typename T>
void Ilupack<T>::InitParameters()
{
  LOG_TRACE2(ilupack) <<  "InitParameters";

  // initializes the parameter block with ilupacks default stuff
  PtrParamNode out = infoNode_->Get(ParamNode::HEADER)->Get("parameters");
  bool disableParallel=false;


  CheckParameter(out, &disableParallel, "disableParallel");
  isParallel=!disableParallel;
  if(isParallel){
    if(isComplex_){
      Exception("Parallel version of ilupack doesn't support complex matrix. Please set disable parallel bool");
    }
    else{
      DSPDAMGinit(&mat, reinterpret_cast<DILUPACKparam*>(&param));
    }
   ConvertIlupackToMatrix (mat, &spr, index);
  }
  else{
    IlupackAMGInit();
  }

  CheckParameter(out, reinterpret_cast<bool*>(&param.matching), "matching");
  CheckParameter(out, &param.ordering, "ordering");
  CheckParameter(out, &param.droptol, "dropTolLU");
  CheckParameter(out, &param.droptolS, "dropTolSchur");
  CheckParameter(out, &param.condest, "condest");
  CheckParameter(out, &param.solver, "iterativeSolver/solver");
  CheckParameter(out, &param.restol, "iterativeSolver/residualTol");
  CheckParameter(out, &param.maxit, "iterativeSolver/maxIter");
  CheckParameter(out, &param.elbow, "elbowSpace");
  CheckParameter(out, &param.amg, "amg");
  if (xml_ != NULL && xml_->Has("iterativeSolver/nrestart"))
  {
    CheckParameter(out, &param.nrestart, "iterativeSolver/nrestart");
  }

  CheckParameter(out, &nleaves, "nleaves");
  // TODO we currently ignore saddle point structures
  param.ind = NULL;
  param.nthreads=std::atoi(getenv("OMP_NUM_THREADS"));

}


template<typename T>
void Ilupack<T>::DetermineMatrixType(BaseMatrix &sysMat, PtrParamNode out)
{
  // first determine it manually for some checking
  if(sysMat.GetStructureType() != BaseMatrix::SPARSE_MATRIX)
    EXCEPTION("Sorry, excpect StdMatrix " << sysMat.GetStructureType());

  const StdMatrix& stdMat = dynamic_cast<const StdMatrix&> (sysMat);
  BaseMatrix::StorageType mst = stdMat.GetStorageType();
  if(mst != BaseMatrix::SPARSE_SYM && mst != BaseMatrix::SPARSE_NONSYM)
    EXCEPTION("Sorry, expect sparse matrix " << mst);

  if (xml_ != NULL && xml_->Has("matrix"))
  {
    PtrParamNode pn = xml_->Get("matrix");
    matrix_ = matrix.Parse( pn->As<std::string>() );
    // plausibility check -- killme: what is with hermitian?
    if (mst != BaseMatrix::SPARSE_SYM && matrix_ != GNL)
      EXCEPTION("Matrix storrage is unsymmetric, so given ilupack_matrix is invalid: '" << matrix.ToString(matrix_) << "'");
  }
  else
  {
    // ignore PD and HER
    matrix_ = mst == BaseMatrix::SPARSE_SYM ? SYM : GNL;
  }

  out->Get("ilupackMatrix")->SetValue(matrix.ToString(matrix_));
}

template<typename T>
void Ilupack<T>::SetEnums()
{
  matrix.SetName("Ilupack::Matrix");
  matrix.Add(GNL, "gnl");
  matrix.Add(SYM, "sym");
  matrix.Add(PD, "pd");
  matrix.Add(HER, "her");
}

template<typename T>
void Ilupack<T>::IlupackAMGInit()
{
  switch(matrix_)
  {
  case GNL:
    if(isComplex_) ZGNLAMGinit(reinterpret_cast<Zmat*>(&mat), reinterpret_cast<ZILUPACKparam*>(&param));
              else DGNLAMGinit(&mat, reinterpret_cast<DILUPACKparam*>(&param));
    break;

  case SYM:
    if(isComplex_) ZSYMAMGinit(reinterpret_cast<Zmat*>(&mat), reinterpret_cast<ZILUPACKparam*>(&param));
              else DSYMAMGinit(&mat, reinterpret_cast<DILUPACKparam*>(&param));
    break;

  case PD:
    if(isComplex_) ZHPDAMGinit(reinterpret_cast<Zmat*>(&mat), reinterpret_cast<ZILUPACKparam*>(&param));
              else DSPDAMGinit(&mat, reinterpret_cast<DILUPACKparam*>(&param));
    break;

  case HER:
    if(isComplex_) ZHERAMGinit(reinterpret_cast<Zmat*>(&mat), reinterpret_cast<ZILUPACKparam*>(&param));
    else throw Exception("ilupack matrix is set to hermitian but not complex");
    break;
  }
}

template<typename T>
int Ilupack<T>::IlupackAMGFactor()
{
  switch(matrix_)
  {
  case GNL:
    if(isComplex_) return ZGNLAMGfactor(reinterpret_cast<Zmat*>(&mat), 
                                        reinterpret_cast<ZAMGlevelmat*>(&precond), 
                                        reinterpret_cast<ZILUPACKparam*>(&param));
              else return DGNLAMGfactor(&mat, &precond, reinterpret_cast<DILUPACKparam*>(&param));

  case SYM:
    if(isComplex_) return ZSYMAMGfactor(reinterpret_cast<Zmat*>(&mat),
                                        reinterpret_cast<ZAMGlevelmat*>(&precond), 
                                        reinterpret_cast<ZILUPACKparam*>(&param));
              else return DSYMAMGfactor(&mat, &precond,  reinterpret_cast<DILUPACKparam*>(&param));

  case PD :
    if(isComplex_) return ZHPDAMGfactor(reinterpret_cast<Zmat*>(&mat),
                                        reinterpret_cast<ZAMGlevelmat*>(&precond), 
                                        reinterpret_cast<ZILUPACKparam*>(&param));
              else return DSPDAMGfactor(&mat, &precond, reinterpret_cast<DILUPACKparam*>(&param));

  case HER:
    if(isComplex_) return ZHERAMGfactor(reinterpret_cast<Zmat*>(&mat), 
                                        reinterpret_cast<ZAMGlevelmat*>(&precond),
                                        reinterpret_cast<ZILUPACKparam*>(&param));
    else throw Exception("ilupack matrix is set to hermitian but not complex");
  }
  
  throw Exception("not handled");
}

template<typename T>
int Ilupack<T>::IlupackAMGSolver(T* sol_ptr, T* rhs_ptr)
{
  switch(matrix_)
  {
    case GNL:
      if(isComplex_) return ZGNLAMGsolver(reinterpret_cast<Zmat*>(&mat),
                                          reinterpret_cast<ZAMGlevelmat*>(&precond), 
                                          reinterpret_cast<ZILUPACKparam*>(&param), 
                                          reinterpret_cast<doublecomplex*>(rhs_ptr),
                                          reinterpret_cast<doublecomplex*>(sol_ptr));
                else return DGNLAMGsolver(&mat, &precond, reinterpret_cast<DILUPACKparam*>(&param), 
                                          reinterpret_cast<double*>(rhs_ptr),
                                          reinterpret_cast<double*>(sol_ptr));
    case SYM:
      if(isComplex_) return ZSYMAMGsolver(reinterpret_cast<Zmat*>(&mat),
                                          reinterpret_cast<ZAMGlevelmat*>(&precond), 
                                          reinterpret_cast<ZILUPACKparam*>(&param), 
                                          reinterpret_cast<doublecomplex*>(rhs_ptr),
                                          reinterpret_cast<doublecomplex*>(sol_ptr));
                else return DSYMAMGsolver(&mat, &precond, reinterpret_cast<DILUPACKparam*>(&param), 
                                          reinterpret_cast<double*>(rhs_ptr),
                                          reinterpret_cast<double*>(sol_ptr));
    case PD:
      if(isComplex_) return ZHPDAMGsolver(reinterpret_cast<Zmat*>(&mat),
                                          reinterpret_cast<ZAMGlevelmat*>(&precond), 
                                          reinterpret_cast<ZILUPACKparam*>(&param), 
                                          reinterpret_cast<doublecomplex*>(rhs_ptr),
                                          reinterpret_cast<doublecomplex*>(sol_ptr));
                else return DSPDAMGsolver(&mat, &precond, reinterpret_cast<DILUPACKparam*>(&param), 
                                          reinterpret_cast<double*>(rhs_ptr),
                                          reinterpret_cast<double*>(sol_ptr));
    case HER:
      if(isComplex_) return ZHERAMGsolver(reinterpret_cast<Zmat*>(&mat),
                                          reinterpret_cast<ZAMGlevelmat*>(&precond), 
                                          reinterpret_cast<ZILUPACKparam*>(&param), 
                                          reinterpret_cast<doublecomplex*>(rhs_ptr),
                                          reinterpret_cast<doublecomplex*>(sol_ptr));
      break;
  }
  throw Exception("not handled");
}

template<typename T>
void Ilupack<T>::IlupackAMGDelete()
{
  LOG_TRACE2(ilupack) <<  "ReleaseMemory: mat_.a=" << mat.a << " .ia=" << mat.ia << ".ja=" << mat.ja;
  // call the ilupack delete method only if Setup() which sets mat was called
  if(mat.a == NULL) return;

  switch(matrix_)
  {
    case GNL:
    if(isComplex_) ZGNLAMGdelete(reinterpret_cast<Zmat*>(&mat), reinterpret_cast<ZAMGlevelmat*>(&precond), reinterpret_cast<ZILUPACKparam*>(&param));
              else DGNLAMGdelete(&mat, &precond, reinterpret_cast<DILUPACKparam*>(&param));
    break;

    case SYM:
    if(isComplex_) ZSYMAMGdelete(reinterpret_cast<Zmat*>(&mat), reinterpret_cast<ZAMGlevelmat*>(&precond), reinterpret_cast<ZILUPACKparam*>(&param));
              else DSYMAMGdelete(&mat, &precond, reinterpret_cast<DILUPACKparam*>(&param));
    break;

    case PD:
    if(isComplex_) ZHPDAMGdelete(reinterpret_cast<Zmat*>(&mat), reinterpret_cast<ZAMGlevelmat*>(&precond), reinterpret_cast<ZILUPACKparam*>(&param));
              else DSPDAMGdelete(&mat, &precond, reinterpret_cast<DILUPACKparam*>(&param));
    break;

    case HER:
    if(isComplex_) ZHERAMGdelete(reinterpret_cast<Zmat*>(&mat), reinterpret_cast<ZAMGlevelmat*>(&precond), reinterpret_cast<ZILUPACKparam*>(&param));
    break;

    default: EXCEPTION("invalid matrix type " << matrix_);
  }

  if(mat.a != NULL)
  { delete[] mat.a; mat.a = NULL;}
  if(mat.ia != NULL)
  { delete[] mat.ia; mat.ia = NULL;}
  if(mat.ja != NULL)
  { delete[] mat.ja; mat.ja = NULL;}
}



template<typename T>
void Ilupack<T>::CalcFillIn(PtrParamNode out)
{
  // this is an extract for the symprintperformance.c sample from Ilupack 2.2
  // It is reduced to the total fill-in factor

  int nnzU = 0, tmp0 = 0; // original names
  DAMGlevelmat  *next = &precond;

  if(param.ind != NULL) EXCEPTION("saddle point is disabled");
  
  for(int i = 1; i <= precond.nlev; i++) 
  {
    if(!(param.flags & DISCARD_MATRIX)) 
    {
      if(i<precond.nlev)
      {
        tmp0 += next->A.ia[next->n]-1;
      }
      else 
      {
        if(next->LU.ja!=NULL) 
          tmp0 += next->A.ia[next->n]-1;
      }
    }
    if(i < precond.nlev || next->LU.ja != NULL) 
    {
      nnzU += next->LU.ja[next->LU.nr-1] - next->LU.ja[0]+2*next->nB;
    }
    if(i == precond.nlev) 
    {
      if(next->LU.ja == NULL) 
      {
        int j = next->LU.nr;
        nnzU += (j*(j-1))/2;
      }
    }
    if(i < precond.nlev) 
    {
      if(param.flags & COARSE_REDUCE) 
      {
        // fill-in F
        nnzU+=next->F.ia[next->F.nr]-1;
      }
    }
    next=next->next;
  }

  out->Get("totalFillInSum")->SetValue(nnzU + mat.nr + tmp0);
  out->Get("totalFillInFactor")->SetValue((1.0 * nnzU + tmp0) / mat.nnz);
}

// Explicit template instantiation
#ifdef EXPLICIT_TEMPLATE_INSTANTIATION
template class Ilupack<Double> ;
template class Ilupack<Complex> ;
#endif

}
