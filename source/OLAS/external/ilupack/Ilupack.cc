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

  precondAge_ = 0;
  setupRuns_ = 0;

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
  const UInt* row_ptr = NULL;
  const UInt* col_ptr = NULL;
  const T* val_ptr = NULL;

  if(som.GetStorageType() == BaseMatrix::SPARSE_SYM)
  {
    const SCRS_Matrix<T>& scrs = dynamic_cast<const SCRS_Matrix<T>&>(som);
    // Nnz is the total number = Diagonal + 2 * triagonals. But one is skipped!
    elements= scrs.GetNnz() - (UInt) (0.5 * (double) (scrs.GetNnz() - scrs.GetNumRows()));

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
  UInt numRows = som.GetNumRows();
  for(UInt i = 0; i < numRows+1; i++)
  mat_.ia[i] += 1;

  for(UInt i = 0; i < elements; i++)
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
void Ilupack<T>::Setup(BaseMatrix &sysMat)
{
  InfoNode* out = solverInfo_->Get(InfoNode::PROCESS)->Get("setup", InfoNode::APPEND);
  
  // determine the matrix type. Symmetric/nonsymmetric, positive definite, ...
  // it is optional given in the xml file.

  // GNL, SYM, PD, HER, ....
  DetermineMatrixType(sysMat);

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
    
  // reset counter which holds the number of rhs's solved with this precond
  precondAge_ = 0;
  setupRuns_++;
}


template<typename T>
void Ilupack<T>::Solve(const BaseMatrix &base_mat,
    const BasePrecond &base_precond, const BaseVector &base_rhs,  BaseVector &base_sol)
{
  InfoNode* out = solverInfo_->Get(InfoNode::PROCESS)->Get("solver", InfoNode::APPEND);
  
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
 
  // solved another rhs with the current preconditioner
  precondAge_++;
}



template<typename T>
void Ilupack<T>::InitParameters()
{
  LOG_TRACE2(ilupack) <<  "InitParameters";

  // initializes the parameter block with ilupacks default stuff
  IlupackAMGInit();
  
  // dump the parameter block and overwrite
  InfoNode* out = solverInfo_->Get(InfoNode::HEADER)->Get("parameters");

  CheckParameter(out, &param.matching, "matching");
  CheckParameter(out, &param.ordering, "ordering");
  CheckParameter(out, &param.droptol, "dropTolLU");
  CheckParameter(out, &param.droptol, "dropTolSchur");
  CheckParameter(out, &param.condest, "condest");
  CheckParameter(out, &param.solver, "iterativeSolver/solver");
  CheckParameter(out, &param.restol, "iterativeSolver/residualTol");
  CheckParameter(out, &param.maxit, "iterativeSolver/maxIter");
  CheckParameter(out, &param.elbow, "elbowSpace");
  CheckParameter(out, &param.amg, "amg");
  
  // TODO we currently igonre saddle point structures
  param.ind = NULL;
  
  Ilupack_symmessages();
}


template<typename T>
void Ilupack<T>::DetermineMatrixType(BaseMatrix &sysMat)
{
  // first determine it manually for some checking
  if (sysMat.GetStructureType() != BaseMatrix::SPARSE_MATRIX)EXCEPTION("Sorry, excpect StdMatrix " << sysMat.GetStructureType());

  const StdMatrix& stdMat = dynamic_cast<const StdMatrix&> (sysMat);
  BaseMatrix::StorageType mst = stdMat.GetStorageType();
  if (mst != BaseMatrix::SPARSE_SYM && mst != BaseMatrix::SPARSE_NONSYM)EXCEPTION("Sorry, expect sparse matrix " << mst);

  if (xml_ != NULL && xml_->Has("matrix"))
  {
    matrix_ = matrix.Parse(xml_->Get("matrix"));
    // plausibility check -- killme: what is with hermitian?
    if (mst != BaseMatrix::SPARSE_SYM && matrix_ != GNL)
      throw Exception(
          "Matrix storrage is unsymmetric, so given ilupack_matrix is invalid "
              + matrix.ToString(matrix_));
  }
  else
  {
    // ignore PD and HER
    matrix_ = mst == BaseMatrix::SPARSE_SYM ? SYM : GNL;
  }

  if (setupRuns_ == 0)
    (*cla) << " Use Ilupack matrix type " << matrix.ToString(matrix_)
        << std::endl;
}

template<typename T>
void Ilupack<T>::SetEnums()
{
  matrix.SetName("Ilupack::Matrix");
  matrix.Add(GNL, "gnl");
  matrix.Add(SYM, "sym");
  matrix.Add(PD, "pd");
  matrix.Add(HER, "her");

  matching.Add(MUMPS, "mumps");
  matching.Add(MC64, "mc64");
  matching.Add(MWM, "schenk_pardiso");

  flags.SetName("Ilupack::Flags");
  flags.Add(FL_DROP_INVERSE, "DropInverse");
  flags.Add(FL_NO_SHIFT, "NoShift");
  flags.Add(FL_TISMENETSKY_SC, "TismenetskySC");
  flags.Add(FL_REPEAT_FACT, "RepeatFact");
  flags.Add(FL_IMPROVED_ESTIMATE, "ImprovedEstimate");
  flags.Add(FL_DIAGONAL_COMPENSATION, "DiagonalCompensation");
  flags.Add(FL_COARSE_REDUCE, "CoarseReduce");
  flags.Add(FL_FINAL_PIVOTING, "FinalPivoting");
  flags.Add(FL_ENSURE_SPD, "EnsureSPD");
  flags.Add(FL_SIMPLE_SC, "SimpleSC");
//  flags.Add(FL_precond_PROCESS_INITIAL_SYSTEM, "PreprocessInitialSystem");
//  flags.Add(FL_precond_PROCESS_SUBSYSTEMS, "PreprocessSubsystem");
  flags.Add(FL_MULTI_PILUC, "MultiPiluc");
  flags.Add(FL_RE_FACTOR, "ReFactor");
  flags.Add(FL_AGGRESSIVE_DROPPING, "AggressiveDropping");
  flags.Add(FL_DISCARD_MATRIX, "DiscardMatrix");
  flags.Add(FL_SYMMETRIC_STRUCTURE, "SymmetricStructure");
}


template<typename T>
void Ilupack<T>::DumpFlags(int flag, std::ostream& out)
{
  out << "State\tFlag" << std::endl;

  std::map<int, std::string>::iterator iter;
  for (iter = flags.map.begin(); iter != flags.map.end(); iter++)
  {
    int f = iter->first;
    out << ((flag & f) != 0 ? "true" : "false");
    out << "\t" << iter->second;
    out << std::endl;
  }
}

template<typename T>
void Ilupack<T>::ApplyFlagSettings(int& flag)
{
  if (xml_ == NULL)
    return;

  // in xml we have elements like this: <flag name="AggressiveDropping" state="off"/>
  // we loop throug all flag names and check of on and off
  StdVector<ParamNode*> all_flags = xml_->GetList("flag");

  for (unsigned int i = 0; i < all_flags.GetSize(); i++)
  {
    ParamNode* pn = all_flags[i];

    int value = flags.Parse(pn->Get("name"));

    if (pn->Get("state")->AsString() == "on")
      flag |= value; // or this value in
    else
      flag &= ~value; // AND 111110111
  }
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
    if(!(param.flags&DISCARD_MATRIX)) 
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

template<typename T>
void Ilupack<T>::Ilupack_symmessages()
{
  
  FILE* fo = fopen("out_ilupack_symmessages","aw");
  
  if (param.matching) {
     if (!strncmp("yes",param.FCpart,3)) {
        if (strncmp("none",param.typetv,4)) {
     // 1. maximum weight matching
     // 2. fine/coarse grid partitioning
     // 3. ... based on static or dynamic test vector
     if (!strncmp("amd",param.ordering,3))
        fprintf(fo,"mw/ad/FCv|");
     else if (!strncmp("metisn",param.ordering,6))
        fprintf(fo,"mw/mn/FCv|");
     else if (!strncmp("metise",param.ordering,6))
        fprintf(fo,"mw/me/FCv|");
     else if (!strncmp("rcm",param.ordering,3))
        fprintf(fo,"mw/rc/FCv|");
     else if (!strncmp("amf",param.ordering,3))
        fprintf(fo,"mw/af/FCv|");
     else if (!strncmp("mmd",param.ordering,3))
        fprintf(fo,"mw/md/FCv|");
     else 
        fprintf(fo,"  /  /   |");
  }
  // 3. no test vector
  else {
     // 1. maximum weight matching
     // 2. fine/coarse grid partitioning
     // 3. no test vector
     if (!strncmp("amd",param.ordering,3))
        fprintf(fo,"mw/ad/FC |");
     else if (!strncmp("metisn",param.ordering,6))
        fprintf(fo,"mw/mn/FC |");
     else if (!strncmp("metise",param.ordering,6))
        fprintf(fo,"mw/me/FC |");
     else if (!strncmp("rcm",param.ordering,3))
        fprintf(fo,"mw/rc/FC |");
     else if (!strncmp("amf",param.ordering,3))
        fprintf(fo,"mw/af/FC |");
     else if (!strncmp("mmd",param.ordering,3))
        fprintf(fo,"mw/md/FC |");
     else 
        fprintf(fo,"  /  /   |");
  }
     }    
     // 2. no FCpart
     else {
  // 1. maximum weight matching
  // 2. no fine/coarse grid partitioning
  // 3. test vector is not used in this case, anyway
  if (!strncmp("amd",param.ordering,3))
     fprintf(fo,"mw/ad/   |");
  else if (!strncmp("metisn",param.ordering,6))
     fprintf(fo,"mw/mn/   |");
  else if (!strncmp("metise",param.ordering,6))
     fprintf(fo,"mw/me/   |");
  else if (!strncmp("rcm",param.ordering,3))
     fprintf(fo,"mw/rc/   |");
  else if (!strncmp("amf",param.ordering,3))
     fprintf(fo,"mw/af/   |");
  else if (!strncmp("mmd",param.ordering,3))
     fprintf(fo,"mw/md/   |");
  else 
     fprintf(fo,"  /  /   |");
     }
  }
  // 1. no matching
  else {
     if (!strncmp("yes",param.FCpart,3)) {
        if (strncmp("none",param.typetv,4)) {
     // 1. no matching
     // 2. fine/coarse grid partitioning
     // 3. ... based on static or dynamic test vector
     if (!strncmp("amd",param.ordering,3))
        fprintf(fo,"  /ad/FCv|");
     else if (!strncmp("metisn",param.ordering,6))
        fprintf(fo,"  /mn/FCv|");
     else if (!strncmp("metise",param.ordering,6))
        fprintf(fo,"  /me/FCv|");
     else if (!strncmp("rcm",param.ordering,3))
        fprintf(fo,"  /rc/FCv|");
     else if (!strncmp("amf",param.ordering,3))
        fprintf(fo,"  /af/FCv|");
     else if (!strncmp("mmd",param.ordering,3))
        fprintf(fo,"  /md/FCv|");
     else 
        fprintf(fo,"  /  /   |");
  }
  // 3. no test vector
  else {
     // 1. no matching
     // 2. fine/coarse grid partitioning
     // 3. no test vector
     if (!strncmp("amd",param.ordering,3))
        fprintf(fo,"  /ad/FC |");
     else if (!strncmp("metisn",param.ordering,6))
        fprintf(fo,"  /mn/FC |");
     else if (!strncmp("metise",param.ordering,6))
        fprintf(fo,"  /me/FC |");
     else if (!strncmp("rcm",param.ordering,3))
        fprintf(fo,"  /rc/FC |");
     else if (!strncmp("amf",param.ordering,3))
        fprintf(fo,"  /af/FC |");
     else if (!strncmp("mmd",param.ordering,3))
        fprintf(fo,"  /md/FC |");
     else 
        fprintf(fo,"  /  /   |");
  }
     }    
     // 2. no FCpart
     else {
  // 1. no matching
  // 2. no fine/coarse grid partitioning
  // 3. test vector is not used in this case, anyway
  if (!strncmp("amd",param.ordering,3))
     fprintf(fo,"  /ad/   |");
  else if (!strncmp("metisn",param.ordering,6))
     fprintf(fo,"  /mn/   |");
  else if (!strncmp("metise",param.ordering,6))
     fprintf(fo,"  /me/   |");
  else if (!strncmp("rcm",param.ordering,3))
     fprintf(fo,"  /rc/   |");
  else if (!strncmp("amf",param.ordering,3))
     fprintf(fo,"  /af/   |");
  else if (!strncmp("mmd",param.ordering,3))
     fprintf(fo,"  /md/   |");
  else 
     fprintf(fo,"  /  /   |");
     }
  }


  printf("ILUPACK PARAMETERS:\n");
  printf("   droptol=%g\n",          param.droptol);
  printf("   condest=%g\n",          param.condest);
  printf("   elbow space factor=%5.1f\n",param.elbow);
  if (!strncmp("ilu",param.typecoarse,3))
     printf("   simple ILU-type coarse grid system\n");
  else if (!strncmp("amg",param.typecoarse,3))
     printf("   medium AMG-type coarse grid system\n");

  if (!strncmp("static",param.typetv,6))
     printf("   diagonal compensation using a prescribed static test vector\n");
  else if (!strncmp("none",param.typetv,4))
     printf("   no diagonal compensation\n");
  else
     printf("   diagonal compensation using a dynamically generated test vector\n");

  printf("   type of preconditioner: ");
  if (!strncmp("ilu",param.amg,3))
     printf("multilevel ILU\n");
  else if (!strncmp("amli",param.amg,4)) {
     printf("AMLI-like\n");
     printf("   number of recursive calls: %d\n", param.ncoarse);
  }
  else if (!strncmp("mg",param.amg,2)) {
     printf("multigrid\n");
     printf("   number of pre-smoothing steps:  %d\n", param.npresmoothing);
     printf("   number of recursive calls:      %d\n", param.ncoarse);
     printf("   number of post-smoothing steps: %d\n", param.npostsmoothing);
     printf("   type of pre-smoother:  ");
     if (!strncmp("gsf",param.presmoother,3))
  printf("Gauss-Seidel forward\n");
     else if (!strncmp("gsb",param.presmoother,3))
  printf("Gauss-Seidel backward\n");
     else if (!strncmp("j",param.presmoother,1))
  printf("(damped) Jacobi\n");
     else if (!strncmp("ilu",param.presmoother,3))
  printf("ILU\n");
     else
  printf("custom\n");
     printf("   type of post-smoother: ");
     if (!strncmp("gsf",param.postsmoother,3))
  printf("Gauss-Seidel forward\n");
     else if (!strncmp("gsb",param.postsmoother,3))
  printf("Gauss-Seidel backward\n");
     else if (!strncmp("j",param.postsmoother,1))
  printf("(damped) Jacobi\n");
     else if (!strncmp("ilu",param.postsmoother,3))
  printf("ILU\n");
     else
  printf("custom\n");
  }



  if (param.matching)
     printf("   maximum weight matching prior to reorderings\n");
  else
     printf("   NO maximum weight matching\n");

  if (param.ind!=NULL)
     printf("   saddle point structure is assumed\n");
  else
     printf("   NO saddle point structure\n");


  if (!strncmp("amd",param.ordering,3))
     printf("   reorder systems based on AMD\n");
  else if (!strncmp("metisn",param.ordering,6))
     printf("   reorder systems based on METIS nested dissection by nodes\n");
  else if (!strncmp("metise",param.ordering,6))
     printf("   reorder systems based on METIS nested dissection by edges\n");
  else if (!strncmp("rcm",param.ordering,3))
     printf("   reorder systems based on RCM\n");
  else if (!strncmp("amf",param.ordering,3))
     printf("   reorder systems based on HALOAMD\n");
  else if (!strncmp("mmd",param.ordering,3))
     printf("   reorder systems based on MMD\n");
  else
     printf("   custom reordering strategy\n");

  if (!strncmp("yes",param.FCpart,3))
     if (strncmp("none",param.typetv,4))
  printf("   a priori fine/coarse grid partitioning using test vector\n");
     else
  printf("   a priori fine/coarse grid partitioning\n");
  else
     printf("   NO a priori fine/coarse grid partitioning\n");

  if (!strncmp("ilu",param.typecoarse,3))
     fprintf(fo,"ILUPACK(S)|");
  else if (!strncmp("amg",param.typecoarse,3))
     fprintf(fo,"ILUPACK(M)|");
  else
     fprintf(fo,"ILUPACK   |");

  fclose(fo);
}


template<typename T>
void Ilupack<T>::Ilupack_symprintperformance()
{
  FILE* fo = fopen("out_ilupack_symmprintperformance","aw");
  
  integer      i,j = 0,k,l,m,n,nz,tmp0,tmp,tmp2,tmp3,ierr,
                nnzU, *ptr = NULL;
  DAMGlevelmat  *next;
  double secnds = 0.0;
    
  printf("   final elbow space factor=%8.2f\n",param.elbow+0.005);
  printf("   final condest on level 1=%8.2f\n",param.condest+0.005);
  printf("ILUPACK,   multilevel structure\n");

  fprintf(fo,"%3d|",precond_.nlev);

  next=&precond_;
  nnzU=0;
  tmp=0;
  ierr=0;
  tmp0=0;
  tmp2=0;
  tmp3=0;
  
  n   = mat_.nr;
  nz  = mat_.nnz;
  
  if (param.ind!=NULL) {
    EXCEPTION("saddle point is disabled");
  }

  for (i=1; i<=precond_.nlev; i++) {
      // fill-in LU
      printf("level %3d, block size %7d\n",i,next->LU.nr); fflush(stdout);
if (param.ind!=NULL) {
   // permute index vector
   for (j=0; j<next->n; j++) 
       ptr[j+next->n]=ptr[next->p[j]-1];
   for (j=0; j<next->n; j++) 
       ptr[j]=ptr[j+next->n];
   k=0;
   for (j=0; j<next->nB; j++) 
       if (ptr[j]<0) k++;
   printf("           saddle point block structure (%6d,%6d),",next->nB-k,k);
   k=0;
   for (; j<next->n; j++) 
       if (ptr[j]<0) k++;
   printf("(%6d,%6d)\n",next->n-next->nB-k,k);
   ptr+=next->nB;
}
if (!(param.flags&DISCARD_MATRIX)) {
   if (i<precond_.nlev) {
      printf("  system size=%6d, fill-in=%8.1f(av.)\n",
       next->n,(next->A.ia[next->n]-1)/((double)next->n));
      fflush(stdout);
      if (i>1) 
   tmp0+=next->A.ia[next->n]-1;
   }
   else {
      if (next->LU.ja!=NULL) {
         printf("  system size=%6d, fill-in=%8.1f(av.)\n",
    next->n,(next->A.ia[next->n]-1)/((double)next->n));
   fflush(stdout);
   tmp0+=next->A.ia[next->n]-1;
      }
   }
}
l=nnzU;
if (i<precond_.nlev || next->LU.ja!=NULL) {
   nnzU+=next->LU.ja[next->LU.nr-1]-next->LU.ja[0]+2*next->nB;
   k=0;
   for (m=0; m<next->LU.nr; m++) {
       if (next->LU.ja[next->LU.nr+1+m]>0) 
    k+=2;
   } // end for m
   ierr+=k;
   tmp3+=next->LU.nr;
   printf("           2x2 pivots %5.1f%%\n",(100.0*k)/next->LU.nr); 
   fflush(stdout);
}
if (i==precond_.nlev) {
   if (next->LU.ja==NULL) {
      printf("switched to full matrix processing\n");fflush(stdout);
      tmp=-1;
      j=next->LU.nr;
      nnzU+=(j*(j-1))/2;
   }
}
printf("  local fill-in %7d(%fav)\n",
       nnzU-l+next->LU.nr,(1.0*(nnzU-l+next->LU.nr))/next->LU.nr);

if (i<precond_.nlev) {
   if (param.flags&COARSE_REDUCE) {
      // fill-in F
      nnzU+=next->F.ia[next->F.nr]-1;
      printf("level %3d->%3d, block size (%7d,%7d)\n",i,i+1,next->LU.nr,
       next->F.nc);
      printf("  local fill-in F %7d(%fav pr)\n",
       next->F.ia[next->F.nr]-1,
       (1.0*(next->F.ia[next->F.nr]-1))/next->LU.nr);
   }
   else {
      printf("level %3d->%3d, block size (%7d,%7d)\n",i,i+1,next->LU.nr,
       next->F.nc);
      fflush(stdout);
   }
}
next=next->next;
  }
  printf("\ntotal fill-in sum%8d(%fav)\n",
   nnzU+n+tmp0,(1.0*(nnzU+tmp0))/n);
  printf("fill-in factor:      %f\n",(1.0*nnzU+tmp0)/nz);
  printf("total number of sparse  2x2 pivots %5.1f%%\n",(100.0*ierr)/tmp3); 
  fflush(stdout);
  ierr=tmp3=0;

  if (tmp) {
     // nnzU-j*(j+1)/2+n-j memory for sparse data structures
     //                    indices (weight 1/3) and values (weight 2/3)
     // j*(j+1)/2          memory for dense data, no indices (weight 2/3)
     printf("memory usage factor: %f\n",(1.0*(tmp0+nnzU-(j*(j+1))/2-j))/nz
      +(j*(j+1))/(3.0*nz));
     fprintf(fo,"%5.1f|",(1.0*(tmp0+nnzU-(j*(j+1))/2-j))/nz
       +(j*(j+1))/(3.0*nz));
  }
  else {
     printf("memory usage factor: %f\n",(1.0*(nnzU+tmp0))/nz);
     fprintf(fo,"%5.1f|",(1.0*(nnzU+tmp0))/nz);
  }
  printf("total time: %e [sec]\n",  (double)secnds); 
  printf("            %e [sec]\n\n",(double)ILUPACK_secnds[7]); 
  fprintf(fo,"%e|",(double)secnds);

  printf("refined timings for   ILUPACK multilevel factorization\n"); 
  printf("initial preprocessing:         %e [sec]\n",ILUPACK_secnds[0]); 
  printf("reorderings remaining levels:  %e [sec]\n",ILUPACK_secnds[1]); 
  printf("SYMPILUC (sum over all levels):%e [sec]\n",ILUPACK_secnds[2]); 
  printf("SYMILUC (if used):             %e [sec]\n",ILUPACK_secnds[3]); 
  printf("SPTRF, LAPACK (if used):       %e [sec]\n",ILUPACK_secnds[4]); 
  printf("remaining parts:               %e [sec]\n\n",std::max(0.0,(double)secnds
                                                   -ILUPACK_secnds[0]
               -ILUPACK_secnds[1]
               -ILUPACK_secnds[2]
               -ILUPACK_secnds[3]
               -ILUPACK_secnds[4]));

  fflush(stdout);

  /*
  next=&precond_;
  for (i=1; i<=nlev; i++) {
    printf("%d,%d\n",next->n,next->nB);
  }
  fflush(stdout);

  printf("multilevel structure\n");
  fflush(stdout);
  printf("number of levels: %d\n",nlev);
  fflush(stdout);
  next=&precond_;
  for (i=1; i<=nlev; i++) {
    printf("total size %d\n",next->n);
    fflush(stdout);
    printf("leading block %d\n",next->nB);
    fflush(stdout);

    printf("row permutation\n");
    for (j=0; j<next->n; j++)
printf("%4d",next->p[j]);
    printf("\n");
    fflush(stdout);
    printf("inverse column permutation\n");
    for (j=0; j<next->n; j++)
printf("%4d",next->invq[j]);
    printf("\n");
    fflush(stdout);

    if (nlev==1) {
printf("row scaling\n");
for (j=0; j<next->n; j++)
  printf("%12.4e",next->rowscal[j]);
printf("\n");
fflush(stdout);
printf("column scaling\n");
for (j=0; j<next->n; j++)
  printf("%12.4e",next->colscal[j]);
printf("\n");
fflush(stdout);
    }

    printf("(2,1) block E (%d,%d)\n",next->E.nr,next->E.nc);
    fflush(stdout);
    for (k=0; k<next->E.nr; k++) {
printf("%3d: ",k+1);
for (j=next->E.ia[k]-1; j<next->E.ia[k+1]-1;j++)
  fprintf(stdout,"%12d",next->E.ja[j]);
printf("\n");
fflush(stdout);
printf("     ");
for (j=next->E.ia[k]-1; j<next->E.ia[k+1]-1;j++)
  fprintf(stdout,"%12.4e",next->E.a[j]);
printf("\n");
fflush(stdout);
    }
    printf("(1,2) block F (%d,%d)\n",next->F.nr,next->F.nc);
    fflush(stdout);
    for (k=0; k<next->F.nr; k++) {
printf("%3d: ",k+1);
for (j=next->F.ia[k]-1; j<next->F.ia[k+1]-1;j++)
  fprintf(stdout,"%12d",next->F.ja[j]);
printf("\n");
fflush(stdout);
printf("     ");
for (j=next->F.ia[k]-1; j<next->F.ia[k+1]-1;j++)
  fprintf(stdout,"%12.4e",next->F.a[j]);
printf("\n");
fflush(stdout);
    }

    printf("ILU...\n");
    printf("Diagonal part\n");
    for (k=0; k<next->LU.nr; k++) {
fprintf(stdout,"%12.4e",next->LU.a[k]);
    }
    printf("\n");
    printf("L part\n");
    for (k=0; k<next->LU.nr; k++) {
printf("col %3d: ",k+1);
for (j=next->LU.ja[k]-1; j<next->LU.ia[k]-1;j++)
  fprintf(stdout,"%12d",next->LU.ja[j]);
printf("\n");
fflush(stdout);
printf("         ");
for (j=next->LU.ja[k]-1; j<next->LU.ia[k]-1;j++)
  fprintf(stdout,"%12.4e",next->LU.a[j]);
printf("\n");
fflush(stdout);
    }
    printf("U part\n");
    for (k=0; k<next->LU.nr; k++) {
printf("row %3d: ",k+1);
for (j=next->LU.ia[k]-1; j<next->LU.ja[k+1]-1;j++)
  fprintf(stdout,"%12d",next->LU.ja[j]);
printf("\n");
fflush(stdout);
printf("         ");
for (j=next->LU.ia[k]-1; j<next->LU.ja[k+1]-1;j++)
  fprintf(stdout,"%12.4e",next->LU.a[j]);
printf("\n");
fflush(stdout);
    }

    next=next->next;
  }
  */

//#ifdef PRINT_INFO
  next=&precond_;
  for (i=1; i<=precond_.nlev; i++) {
      // fill-in LU
      printf("level %3d, block size %7d\n",i,next->LU.nr); fflush(stdout);
if (i<precond_.nlev || next->LU.ja!=NULL) {
   printf("U-factor");
   printf("\n");fflush(stdout);
   for (l=0; l<next->LU.nr; ) {
       if (next->LU.ja[next->LU.nr+1+l]==0){
    for (j=next->LU.ja[l];j<next->LU.ja[l+1]; j++) {
        printf("%8d",next->LU.ja[j-1]);
    }
    printf("\n");fflush(stdout);
    for (j=next->LU.ja[l];j<next->LU.ja[l+1]; j++) {
        printf("%8.1e",next->LU.a[j-1]);
    }
    l++;
       }
       else {
    for (j=next->LU.ja[l];j<next->LU.ja[l+1]; j++) {
        printf("%8d",next->LU.ja[j-1]);
    }
    printf("\n");fflush(stdout);
    for (j=next->LU.ja[l];j<next->LU.ja[l+1]; j++) {
        printf("%8.1e",
         next->LU.a[next->LU.ja[l]+2*(j-next->LU.ja[l])-1]);
    }
    printf("\n");fflush(stdout);
    for (j=next->LU.ja[l];j<next->LU.ja[l+1]; j++) {
        printf("%8.1e",
         next->LU.a[next->LU.ja[l]+2*(j-next->LU.ja[l])]);
    }
    l+=2;
       }
       printf("\n");fflush(stdout);
   }

   printf("Block diagonal factor\n");
   for (l=0; l<next->LU.nr;) {
       if (next->LU.ja[next->LU.nr+1+l]==0){
    printf("%8.1e",next->LU.a[l]);
    l++;
       }
       else {
    printf("%8.1e%8.1e",next->LU.a[l],
     next->LU.a[next->LU.nr+1+l]);
    l+=2;
       }
   }
   printf("\n");fflush(stdout);
   for (l=0; l<next->LU.nr; ) {
       if (next->LU.ja[next->LU.nr+1+l]==0) {
   printf("        ");
   l++;
       }
       else {
   printf("%e%e",next->LU.a[next->LU.nr+1+l],
    next->LU.a[l+1]);
   l+=2;
       }
   }
   printf("\n");fflush(stdout);
   
}
if (i==precond_.nlev) {
   if (next->LU.ja==NULL) {
      printf("switched to full matrix processing\n");fflush(stdout);
   }
}

if (i<precond_.nlev) {
   // fill-in F
   nnzU+=next->F.ia[next->F.nr]-1;
   printf("level %3d->%3d, block size (%7d,%7d)\n",i,i+1,next->LU.nr,
    next->F.nc);
   printf("  local fill-in F %7d(%fav pr)\n",
    next->F.ia[next->F.nr]-1,
    (1.0*(next->F.ia[next->F.nr]-1))/next->LU.nr);
}
next=next->next;
  }
//#endif

}




/*
template<typename T>
string Ilupack<T>::DumpPrecond(DAMGlevelmat* precond)
{
  std::stringstream ss;
  ss << "level=" << precond->nlev << " n=" << precond->n << " nB=" << precond->nB << std::endl;
  //ss << "A=" << 
}

*/



/*
template<typename T>
void Ilupack<T>::Ilupack_symfinalres(T* sol_ptr, T* rhs_ptr)
{
  integer      i,l,n,nz,tmp0,tmp,tmp2,tmp3,ierr,
                  nnzU,;
    nnzU=0;
    tmp=0;
    ierr=0;
    tmp0=0;
    tmp2=0;
    tmp3=0;
    
    n   = mat_.nr;
    nz  = mat_.nnz;
  
  double* sol = (double*) sol_ptr;
  double* rhs = (double*) rhs_ptr;
  
  // -------   compute final residual   ------
  DSYMmatvec(mat_,sol,param.dbuff);

  for (i=0; i<n; i++) {
    param.dbuff[i]-=rhs[i+mat_.nr*l];
  } // end for i

  i=1;
  val=NRM(&n, param.dbuff, &i);
  // -----------------------------------------
  printf("current: %8.1le\n",val);


  // release part of rhs that may store the uncompressed rhs
  if (nrhs!=0 && (rhstyp[0]=='M' || rhstyp[0]=='m')) {
    rhs+=n;
    rhs+=mat_.nr*l;
  }


  for (i=0; i<n; i++) {

    if (l==0 && param.flags&DIAGONAL_COMPENSATION &&
        param.flags&(STATIC_TESTVECTOR|DYNAMIC_TESTVECTOR))
      sol[i+mat_.nr*l]-=mytestvector[i];
    else
      sol[i+mat_.nr*l]-=1.0+i*l;
  }
  i=1;
  val=NRM(&n,sol+mat_.nr*l,&i);
  for (i=0; i<n; i++) {
    if (l==0 && param.flags&DIAGONAL_COMPENSATION &&
        param.flags&(STATIC_TESTVECTOR|DYNAMIC_TESTVECTOR))
      sol[i+mat_.nr*l]+=mytestvector[i];
    else
      sol[i+mat_.nr*l]+=1.0+i*l;
  }
  i=1;
  vb=NRM(&n,sol+mat_.nr*l,&i);
  // if (nrhs!=0 && ((rhstyp[2]=='X' || rhstyp[2]=='x') || (rhstyp[0]!='M' && rhstyp[0]!='m')))
  printf("rel. error in the solution: %8.1le\n\n",val/vb);
  //else printf("\n");

}


*/


/*
typedef struct  DAMGLM {
  integer nlev;                  
  integer n;                  
  integer nB;
  Dmat A; 
  Dmat LU;
  integer *LUperm;
  Dmat E;
  Dmat F;
  integer *p;
  integer *invq;
  doubleprecision *rowscal;
  doubleprecision *colscal;
  doubleprecision *absdiag;
  struct DAMGLM *prev;
  struct DAMGLM *next;
  integer *nextblock;
  integer issymmetric;
  integer isdefinite;
  integer ishermitian;
  integer isskew;
  integer isreal;
  integer issingle;
} DAMGlevelmat; 

*/

// Explicit template instantiation
#ifdef EXPLICIT_TEMPLATE_INSTANTIATION
template class Ilupack<Double> ;
template class Ilupack<Complex> ;
#endif
