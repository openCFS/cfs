// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <iostream>
#include <string>
#include <algorithm>

#include <def_use_pardiso.hh>
#include <def_use_metis.hh>

#include "General/exception.hh"
#include "MatVec/sparseolasmatrix.hh"
#include "MatVec/scrs_matrix.hh"
#include "MatVec/crs_matrix.hh"
#include "DataInOut/Logging/cfslog.hh"
#include "DataInOut/ParamHandling/InfoNode.hh"
#include "Utils/StdVector.hh"
#include "Ilupack.hh"

DECLARE_LOG(ilupack)
DEFINE_LOG(ilupack, "ilupack")

using namespace CoupledField;

template<typename T>
Ilupack<T>::Ilupack(ParamNode* xml, InfoNode* olasInfo, BaseMatrix::EntryType type)
{
  // we work with out 
  xml_ = xml != NULL && xml->Has("ilupack") ? xml->Get("ilupack") : NULL;
  solverInfo_ = olasInfo->Get("ilupack");

  if (type != BaseMatrix::COMPLEX && type != BaseMatrix::DOUBLE)EXCEPTION("unhandled type " << type);
  isComplex_ = type == BaseMatrix::COMPLEX;

  LOG_TRACE(ilupack) <<  "Ilupack(): isComplex=" << isComplex_;

  // set this to zero to signal the destructor if Setup() was called
  mat_.a = NULL;
  mat_.ia = NULL;
  mat_.ja = NULL;

  // set the convenience pointers for the D* instances of ilupack structs.
  dmat_ptr_ = &mat_;
  zmat_ptr_ = reinterpret_cast<Zmat*>(&mat_);

  dparam_ptr_ = &param;
  zparam_ptr_ = reinterpret_cast<ZILUPACKparam*>(&param);

  dprecond_ptr_ = &precond_;
  zprecond_ptr_ = reinterpret_cast<ZAMGlevelmat*>(&precond_);

  // init the enums
  SetEnums();
}

template<typename T>
Ilupack<T>::~Ilupack()
{
  LOG_TRACE(ilupack) <<  "~Ilupack()";
  IlupackAMGDelete();
}

template<typename T>
void Ilupack<T>::SetMatrix(const BaseMatrix &base_mat)
{
  const SparseOLASMatrix<T>& som =  dynamic_cast<const SparseOLASMatrix<T>&> (base_mat);

  LOG_TRACE2(ilupack) <<  "SetMatrix: SPARSE_SYM="  << (som.GetStorageType() == BaseMatrix::SPARSE_SYM)
                      << " mat_.a=" << mat_.a << " .ia=" << mat_.ia << ".ja=" << mat_.ja;

  // delete the old values before we allocate the new space - if could be faster! 
  if(mat_.a != NULL)
  { delete[] mat_.a; mat_.a = NULL;}
  if(mat_.ia != NULL)
  { delete[] mat_.ia; mat_.ia = NULL;}
  if(mat_.ja != NULL)
  { delete[] mat_.ja; mat_.ja = NULL;}

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
  mat_.ia = new int[som.GetNumRows() + 1]; // one plus due to CRS tail
  std::copy(row_ptr, row_ptr + som.GetNumRows() + 1, mat_.ia);

  mat_.ja = new int[elements];
  std::copy(col_ptr, col_ptr + elements, mat_.ja);

  // values is brutal stuff, because binary the std::complex<double>* IS a double*!
  mat_.a = reinterpret_cast<double*>(new T[elements]);
  T* val_src = reinterpret_cast<T*>(mat_.a);
  std::copy(val_ptr, val_ptr + elements, val_src);

  // no adjust to 1-based!
  for(unsigned int i = 0, numRows = som.GetNumRows(); i < numRows+1; i++)
    mat_.ia[i] += 1;

  for(unsigned int i = 0; i < elements; i++)
    mat_.ja[i] += 1;

  mat_.nr = som.GetNumRows();
  mat_.nc = som.GetNumCols();
  mat_.nnz = elements;

  LOG_TRACE2(ilupack) << "SetMatrix: allocate: mat_.a<T>=" << elements << " .ia=" << (som.GetNumRows() + 1)
                      << ".ia=" << elements << "; .nr=" << mat_.nr << " .nc=" << mat_.nc;
  LOG_DBG2(ilupack) << "mat_.ia: " << StdVector<int>::ToString(mat_.nr + 1, mat_.ia);
  LOG_DBG2(ilupack) << "mat_.ja: " << StdVector<int>::ToString(elements, mat_.ja);
  LOG_DBG2(ilupack) << "mat_.a: " << StdVector<double>::ToString(elements, mat_.a);
}

template<typename T>
void Ilupack<T>::Setup(BaseMatrix &sysMat, InfoNode* analysis_id)
{
  InfoNode* out = solverInfo_->Get(InfoNode::PROCESS)->Get("setup", InfoNode::APPEND);
  out->Get("analysis_id")->SetValue(analysis_id->Get("analysis_id"));
  
  // determine the matrix type. Symmetric/nonsymmetric, positive definite, ...
  // it is optional given in the xml file.

  // GNL, SYM, PD, HER, ....
  DetermineMatrixType(sysMat, out);

  LOG_TRACE2(ilupack) <<  "Setup: matrix -> " << matrix.ToString(matrix_);

  // in case we already run release memory - it's save if first run
  IlupackAMGDelete();

  // set the ilupack matrix mat_ from the olas matrix - complex can be casted to double
  SetMatrix(sysMat);

  // gain and modify the default settings, does logging uses mat_ and sets param 
  InitParameters();

  // factorize the iLU preconditioner
  std::stringstream ss;
  ss << "Error facorizing Ilupack: ";

  int ierr = IlupackAMGFactor();

  switch (ierr)
  {
    case 0: 
    // perfect:
    break;

    case -1: 
    ss << "Input matrix may be wrong at level " << precond_.nlev;
    break;

    case -2:
    ss << "Out of memory. The matrix L overflows the array alu at level " << precond_.nlev;
    break;

    case -3:
    ss << "Out of memory. The matrix U overflows the array alu at level " << precond_.nlev;
    break;

    case -4:
    ss << "Illegal value for lfil at level " << precond_.nlev;
    break;

    case -5:
    ss << "Zero row encountered at level " <<  precond_.nlev;
    break;

    case -6:
    ss << "Zero column encountered at level " <<  precond_.nlev;
    break;

    case -7:
    ss << "Buffers are too small";

    default:
      ss << "Zero pivot encountered at step number " << ierr << " of level " << precond_.nlev;
  }
  
  if(ierr != 0) throw Exception(ss.str());
 
  out->Get("levels")->SetValue(precond_.nlev);
  CalcFillIn(out);
  InfoNode* timing = out->Get("timing");
  timing->Get("total_time")->SetValue(ILUPACK_secnds[7]);
  timing->Get("initial_preprocessing")->SetValue(ILUPACK_secnds[0]);
  timing->Get("reordering_remaining_levels")->SetValue(ILUPACK_secnds[1]);
}


template<typename T>
void Ilupack<T>::Solve(const BaseMatrix &base_mat, const BasePrecond &base_precond, 
    const BaseVector &base_rhs,  BaseVector &base_sol, InfoNode* analysis_id)
{
  InfoNode* out = solverInfo_->Get(InfoNode::PROCESS)->Get("solver", InfoNode::APPEND);
  out->Get("analysis_id")->SetValue(analysis_id->Get("analysis_id"));
  
  // the preconditioner sets the ilupack matrix
  if (mat_.a == NULL)
    throw Exception("Setup() not called before Solve()");

  // we gain the solution pointer
  T* sol_ptr = dynamic_cast<Vector<T>&> (base_sol).GetPointer();

  // additionally remove the constness of the rhs, ilupack knows no const.
  T* rhs_ptr =  const_cast<T*> (dynamic_cast<const Vector<T>&> (base_rhs).GetPointer());

  LOG_TRACE2(ilupack) <<  "Solve: sol_ptr=" << sol_ptr << " rhs_ptr=" << rhs_ptr;

  std::stringstream ss;
  ss << "Error solving Ilupack: ";
  int ierr = IlupackAMGSolver(sol_ptr, rhs_ptr);

  // why did the iterative solver stop?
  switch (ierr) 
  {
    case  0:  // everything is fine
      break;

    case -1:  // too many iterations
      ss << "Number of iteration steps exceeds its limit";
      break;

    case -2: 
      ss << "Not enough work space provided";
      break;

    case -3:  /* not enough work space */
      ss << "Algorithm breaks down";
      break;

    default: 
      ss << "Solver exited with error code " << ierr;
  } 

  // stop if necessary
  if(ierr != 0) throw Exception(ss.str());

  out->Get("iterations")->SetValue(param.ipar[26]);
  InfoNode* timing = out->Get("timing");
  timing->Get("total")->SetValue(ILUPACK_secnds[5]);
  timing->Get("maxtrix_vector_mult")->SetValue(ILUPACK_secnds[6]);
  InfoNode* norms = out->Get("norms");
  norms->Get("target")->SetValue(param.fpar[23]);
}



template<typename T>
void Ilupack<T>::InitParameters()
{
  LOG_TRACE2(ilupack) <<  "InitParameters";

  // initializes the parameter block with ilupacks default stuff
  IlupackAMGInit();
  
  // dump the parameter block and overwrite
  InfoNode* out = solverInfo_->Get(InfoNode::HEADER)->Get("parameters");

  CheckParameter(out, reinterpret_cast<bool*>(&param.matching), "matching");
  CheckParameter(out, &param.ordering, "ordering");
  CheckParameter(out, &param.droptol, "dropTolLU");
  CheckParameter(out, &param.droptol, "dropTolSchur");
  CheckParameter(out, &param.condest, "condest");
  CheckParameter(out, &param.solver, "iterativeSolver/solver");
  CheckParameter(out, &param.restol, "iterativeSolver/residualTol");
  CheckParameter(out, &param.maxit, "iterativeSolver/maxIter");
  CheckParameter(out, &param.elbow, "elbowSpace");
  CheckParameter(out, &param.amg, "amg");
  
  // TODO we currently ignore saddle point structures
  param.ind = NULL;
}


template<typename T>
void Ilupack<T>::DetermineMatrixType(BaseMatrix &sysMat, InfoNode* out)
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
    matrix_ = matrix.Parse(xml_->Get("matrix"));
    // plausibility check -- killme: what is with hermitian?
    if (mst != BaseMatrix::SPARSE_SYM && matrix_ != GNL)
      EXCEPTION("Matrix storrage is unsymmetric, so given ilupack_matrix is invalid " << matrix.ToString(matrix_));
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
void Ilupack<T>::CheckParameter(InfoNode* out, char** ilupack_string,
    const char* param_name)
{
  InfoNode* tmp = out->Get(param_name);
  tmp->Get("default")->SetValue(*ilupack_string);
  if (xml_ != NULL && xml_->Has(param_name))
  {
    *ilupack_string
        = const_cast<char*> (xml_->Get(param_name)->AsString().c_str());
    tmp->Get("set")->SetValue(*ilupack_string);
  }
}

template<typename T>
void Ilupack<T>::CheckParameter(InfoNode* out, double* ilupack_val, const char* param_name)
{
  InfoNode* tmp = out->Get(param_name);
  tmp->Get("default")->SetValue(*ilupack_val);
  if (xml_ != NULL && xml_->Has(param_name))
  {
    *ilupack_val = xml_->Get(param_name)->AsDouble();
    tmp->Get("set")->SetValue(*ilupack_val);
  }
}

template<typename T>
void Ilupack<T>::CheckParameter(InfoNode* out, int* ilupack_val, const char* param_name)
{
  InfoNode* tmp = out->Get(param_name);
  tmp->Get("default")->SetValue(*ilupack_val);
  if (xml_ != NULL && xml_->Has(param_name))
  {
    *ilupack_val = xml_->Get(param_name)->AsInt();
    tmp->Get("set")->SetValue(*ilupack_val);
  }
}

template<typename T>
void Ilupack<T>::CheckParameter(InfoNode* out, bool* ilupack_val, const char* param_name)
{
  // by convention we interpret this as "integer"
  integer* int_ptr = reinterpret_cast<integer*>(ilupack_val);
  
  InfoNode* tmp = out->Get(param_name);
  tmp->Get("default")->SetValue(*ilupack_val);
  if (xml_ != NULL && xml_->Has(param_name))
  {
    *int_ptr = xml_->Get(param_name)->AsBool() == false ? 0 : 1;
    tmp->Get("set")->SetValue(*int_ptr == 0 ? false : true);
  }
}


template<typename T>
void Ilupack<T>::IlupackAMGInit()
{
  switch(matrix_)
  {
  case GNL:
    if(isComplex_) ZGNLAMGinit(zmat_ptr_, zparam_ptr_);
    else DGNLAMGinit(dmat_ptr_, dparam_ptr_);
    break;

  case SYM:
    if(isComplex_) ZSYMAMGinit(zmat_ptr_, zparam_ptr_);
    else DSYMAMGinit(dmat_ptr_, dparam_ptr_);
    break;

  case PD:
    if(isComplex_) ZHPDAMGinit(zmat_ptr_, zparam_ptr_);
    else DSPDAMGinit(dmat_ptr_, dparam_ptr_);
    break;

  case HER:
    if(isComplex_) ZHERAMGinit(zmat_ptr_, zparam_ptr_);
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
    if(isComplex_) return ZGNLAMGfactor(zmat_ptr_, zprecond_ptr_, zparam_ptr_);
    else return DGNLAMGfactor(dmat_ptr_, dprecond_ptr_, dparam_ptr_);

  case SYM:
    if(isComplex_) return ZSYMAMGfactor(zmat_ptr_, zprecond_ptr_, zparam_ptr_);
    else return DSYMAMGfactor(dmat_ptr_, dprecond_ptr_, dparam_ptr_);

  case PD :
    if(isComplex_) return ZHPDAMGfactor(zmat_ptr_, zprecond_ptr_, zparam_ptr_);
    else return DSPDAMGfactor(dmat_ptr_, dprecond_ptr_, dparam_ptr_);

  case HER:
    if(isComplex_) return ZHERAMGfactor(zmat_ptr_, zprecond_ptr_, zparam_ptr_);
    else throw Exception("ilupack matrix is set to hermitian but not complex");
  }
  
  throw Exception("not handled");
}

template<typename T>
int Ilupack<T>::IlupackAMGSolver(T* sol_ptr, T* rhs_ptr)
{
  double *d_sol, *d_rhs;
  // is a ilupack type
  doublecomplex *z_sol, *z_rhs;

  d_sol = isComplex_ ? NULL : reinterpret_cast<double*>(sol_ptr);
  d_rhs = isComplex_ ? NULL : reinterpret_cast<double*>(rhs_ptr);
  z_sol = isComplex_ ? reinterpret_cast<doublecomplex*>(sol_ptr) : NULL;
  z_rhs = isComplex_ ? reinterpret_cast<doublecomplex*>(rhs_ptr) : NULL;
  

  switch(matrix_)
  {
    case GNL:
      if(isComplex_) return ZGNLAMGsolver(zmat_ptr_, zprecond_ptr_, zparam_ptr_, z_rhs, z_sol);
                else return DGNLAMGsolver(dmat_ptr_, dprecond_ptr_, dparam_ptr_, d_rhs, d_sol);

    case SYM:
      if(isComplex_) return ZSYMAMGsolver(zmat_ptr_, zprecond_ptr_, zparam_ptr_, z_rhs, z_sol);
                else return DSYMAMGsolver(dmat_ptr_, dprecond_ptr_, dparam_ptr_, d_rhs, d_sol);

    case PD:
      if(isComplex_) return ZHPDAMGsolver(zmat_ptr_, zprecond_ptr_, zparam_ptr_, z_rhs, z_sol);
                else return DSPDAMGsolver(dmat_ptr_, dprecond_ptr_, dparam_ptr_, d_rhs, d_sol);
      break;

    case HER:
      if(isComplex_) return ZHERAMGsolver(zmat_ptr_, zprecond_ptr_, zparam_ptr_, z_rhs, z_sol);
      break;
  }
  throw Exception("not handled");
}

template<typename T>
void Ilupack<T>::IlupackAMGDelete()
{
  LOG_TRACE2(ilupack) <<  "ReleaseMemory: mat_.a=" << mat_.a << " .ia=" << mat_.ia << ".ja=" << mat_.ja;
  // call the ilupack delete method only if Setup() which sets mat_ was called
  if(mat_.a == NULL) return;

  switch(matrix_)
  {
    case GNL:
    if(isComplex_) ZGNLAMGdelete(zmat_ptr_, zprecond_ptr_, zparam_ptr_);
    else DGNLAMGdelete(dmat_ptr_, dprecond_ptr_, dparam_ptr_);
    break;

    case SYM:
    if(isComplex_) ZSYMAMGdelete(zmat_ptr_, zprecond_ptr_, zparam_ptr_);
    else DSYMAMGdelete(dmat_ptr_, dprecond_ptr_, dparam_ptr_);
    break;

    case PD:
    if(isComplex_) ZHPDAMGdelete(zmat_ptr_, zprecond_ptr_, zparam_ptr_);
    else DSPDAMGdelete(dmat_ptr_, dprecond_ptr_, dparam_ptr_);
    break;

    case HER:
    if(isComplex_) ZHERAMGdelete(zmat_ptr_, zprecond_ptr_, zparam_ptr_);
    break;

    default: EXCEPTION("invalid matrix type " << matrix_);
  }

  if(mat_.a != NULL)
  { delete[] mat_.a; mat_.a = NULL;}
  if(mat_.ia != NULL)
  { delete[] mat_.ia; mat_.ia = NULL;}
  if(mat_.ja != NULL)
  { delete[] mat_.ja; mat_.ja = NULL;}
}



template<typename T>
void Ilupack<T>::CalcFillIn(InfoNode* out)
{
  // this is an extract for the symprintperformance.c sample from Ilupack 2.2
  // It is reduced to the total fill-in factor

  int nnzU = 0, tmp0 = 0; // original names
  DAMGlevelmat  *next = &precond_;

  if(param.ind != NULL) EXCEPTION("saddle point is disabled");
  
  for(int i = 1; i <= precond_.nlev; i++) 
  {
    if(!(param.flags & DISCARD_MATRIX)) 
    {
      if(i<precond_.nlev)
      {
        tmp0 += next->A.ia[next->n]-1;
      }
      else 
      {
        if(next->LU.ja!=NULL) 
          tmp0 += next->A.ia[next->n]-1;
      }
    }
    if(i < precond_.nlev || next->LU.ja != NULL) 
    {
      nnzU += next->LU.ja[next->LU.nr-1] - next->LU.ja[0]+2*next->nB;
    }
    if(i == precond_.nlev) 
    {
      if(next->LU.ja == NULL) 
      {
        int j = next->LU.nr;
        nnzU += (j*(j-1))/2;
      }
    }
    if(i < precond_.nlev) 
    {
      if(param.flags & COARSE_REDUCE) 
      {
        // fill-in F
        nnzU+=next->F.ia[next->F.nr]-1;
      }
    }
    next=next->next;
  }

  out->Get("totalFillInSum")->SetValue(nnzU + mat_.nr + tmp0);
  out->Get("totalFillInFactor")->SetValue((1.0 * nnzU + tmp0) / mat_.nnz);
}

// Explicit template instantiation
#ifdef EXPLICIT_TEMPLATE_INSTANTIATION
template class Ilupack<Double> ;
template class Ilupack<Complex> ;
#endif
