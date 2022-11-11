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
#include "CholMod.hh"

DEFINE_LOG(cholmod, "cholmod")

namespace CoupledField
{

template<typename T>
CholMod<T>::CholMod(PtrParamNode xml, PtrParamNode olasInfo, BaseMatrix::EntryType type)
{
  // we work with out
  xml_ =  xml;
  infoNode_ = olasInfo->Get("cholmod");
  
  if(type != BaseMatrix::COMPLEX && type != BaseMatrix::DOUBLE){
    EXCEPTION("unhandled type " << type);
  }
  isComplex_ = type == BaseMatrix::COMPLEX;

  LOG_TRACE(cholmod) <<  "CholMod(): isComplex=" << isComplex_;
  
  SetEnums();
  
  cholmod_start(&common_);
  
  fact_ = NULL;
  mat_ = new cholmod_sparse;
  mat_->nzmax = 0;

  InitParameters();
}

template<typename T>
CholMod<T>::~CholMod()
{
  CholModFreeFactorization();
  delete mat_; // this should not free pointers contained in mat, just mat itself
  mat_ = NULL;
  cholmod_finish(&common_);
  LOG_TRACE(cholmod) <<  "~CholMod()";
}

template<typename T>
void CholMod<T>::SetMatrix(const BaseMatrix &base_mat)
{
  const SparseOLASMatrix<T>& som =  dynamic_cast<const SparseOLASMatrix<T>&> (base_mat);

  if(som.GetStorageType() != BaseMatrix::SPARSE_SYM){
    EXCEPTION("CholMod can only be used for s.p.d. matrices!!!");
  }
  
  const SCRS_Matrix<T>& scrs = dynamic_cast<const SCRS_Matrix<T>&>(som);
  
  UInt n = som.GetNumRows(); // = Cols
  
  // we already have memory allocated in constructor for the container
  // we just inject the matrix into that container
  mat_->nrow = n;
  mat_->ncol = n;
  mat_->nzmax = scrs.GetNumEntries();
  
  // note that we transpose the matrix here, as storage in CholMod is Transposed compared to OLAS
  mat_->p = const_cast<void*>(static_cast<const void*>(scrs.GetRowPointer()));
  mat_->i = const_cast<void*>(static_cast<const void*>(scrs.GetColPointer()));
  mat_->x = const_cast<void*>(static_cast<const void*>(scrs.GetDataPointer()));
  
  /*std::cout<<"rowPtr = ";
  for (int i = 0; i <= (int) n; i++) {
    std::cout<<scrs.GetRowPointer()[i]<<", ";
  }
  std::cout<<std::endl;

  std::cout<<"colPtr = ";
  for (int i = 0; i < (int) mat_->nzmax; i++) {
    std::cout<<scrs.GetColPointer()[i]<<", ";
  }
  std::cout<<std::endl;

        std::cout<<"dataPtr = ";
        for (int i = 0; i < (int) mat_->nzmax; i++) {
          std::cout<<scrs.GetDataPointer()[i]<<", ";
        }
        std::cout<<std::endl; */

  mat_->nz = NULL;
  mat_->stype = -1; // s.p.d. matrix, using lower triangular part;
  mat_->itype = CHOLMOD_INT;
  mat_->xtype = isComplex_ ? CHOLMOD_COMPLEX : CHOLMOD_REAL;
  mat_->dtype = CHOLMOD_DOUBLE;
  mat_->sorted = true;
  mat_->packed = true;
}

template<typename T>
void CholMod<T>::Setup(BaseMatrix &sysMat)
{
  // do we really want to create a new entry? Might blast up the output
  ParamNode::ActionType at = progOpts->DoDetailedInfo() ? ParamNode::APPEND : ParamNode::DEFAULT;
  PtrParamNode out = infoNode_->Get(ParamNode::PROCESS)->Get("setup", at);
  // out->Get("analysis_id")->SetValue(analysis_id->Get("analysis_id"));
  
  LOG_TRACE2(cholmod) <<  "Setup: matrix -> " << sysMat.ToString();

  // in case we already run release memory - it's save if first run
  CholModFreeFactorization();

  // interface the matrix to CholMod
  SetMatrix(sysMat);

  CholModAnalyze();
  CholModFactorize();

  // output statistics to info.xml, see CholMod manual for details
  out->Get("selected_factorization_method")->SetValue(common_.selected);
  out->Get("ll_flop_count")->SetValue(common_.fl);
  out->Get("nz_in_L")->SetValue(common_.lnz);
  out->Get("nz_in_tri(A)")->SetValue(common_.anz);
  out->Get("memory_usage")->SetValue(common_.memory_usage);
  out->Get("memory_inuse")->SetValue(common_.memory_inuse);
  out->Get("malloc_count")->SetValue(common_.malloc_count);
  out->Get("nrealloc_col")->SetValue(common_.nrealloc_col);
  out->Get("nrealloc_factor")->SetValue(common_.nrealloc_factor);
  out->Get("ndbounds_hit")->SetValue(common_.ndbounds_hit);
  int u = common_.nmethods + 1;
  if(common_.current + 1 > u){
    u = common_.current + 1;
  }
  if(u > CHOLMOD_MAXMETHODS){
    u = CHOLMOD_MAXMETHODS+1;
  }
  for(int i = 0; i < u; ++i){
    PtrParamNode m = out->Get("method");
    m->Get("id")->SetValue(i);
    m->Get("lnz")->SetValue(common_.method[i].lnz);
    m->Get("fl")->SetValue(common_.method[i].fl);
  }
}


template<typename T>
void CholMod<T>::Solve(const BaseMatrix &base_mat, 
    const BaseVector &base_rhs,  BaseVector &base_sol)

{
  ParamNode::ActionType at = progOpts->DoDetailedInfo() ? ParamNode::APPEND : ParamNode::DEFAULT;
  PtrParamNode out = infoNode_->Get(ParamNode::PROCESS)->Get("solver", at);
  // out->Get("analysis_id")->SetValue(analysis_id->Get("analysis_id"));
  
  // the preconditioner sets the matrix
  if (mat_->nzmax == 0)
    throw Exception("Setup() not called before Solve()");

  // we fetch the rhs pointer
  T* rhs_ptr =  const_cast<T*>(dynamic_cast<const Vector<T>&> (base_rhs).GetPointer());

  LOG_TRACE2(cholmod) <<  "Solve: rhs_ptr=" << rhs_ptr;
  
  // build a cholmod_dense container for the rhs
  UInt s = base_rhs.GetSize();  
  cholmod_dense* rhs = new cholmod_dense; // manual allocation (this is just the container)
  rhs->nrow = s;
  rhs->ncol = 1;
  rhs->nzmax = s;
  rhs->d = s;
  rhs->x = rhs_ptr;
  rhs->xtype = isComplex_ ? CHOLMOD_COMPLEX : CHOLMOD_REAL;
  rhs->dtype = CHOLMOD_DOUBLE;
  
  // now solve the problem
  cholmod_dense* sol = cholmod_solve(CHOLMOD_A, fact_, rhs, &common_);
  TestException("solve");

  delete rhs; // free the container again

  // we gain the solution pointer and copy back the result
  T* sol_ptr = dynamic_cast<Vector<T>&> (base_sol).GetPointer();
  std::copy(static_cast<T*>(sol->x), static_cast<T*>(sol->x) + sol->nrow, sol_ptr);
  cholmod_free_dense(&sol, &common_); // free the temp. solution
  LOG_TRACE2(cholmod) <<  "Solve: sol= " << base_sol.ToString();
}



template<typename T>
void CholMod<T>::InitParameters()
{
  LOG_TRACE2(cholmod) <<  "InitParameters";

  // initializes the parameter block with default stuff
  CholModDefaultParams();
  
  // dump the parameter block and overwrite
  PtrParamNode out = infoNode_->Get(ParamNode::HEADER)->Get("parameters");
  
  CheckParameter(out, &common_.nmethods, "factorization/nmethods");
  CheckParameter(out, &common_.postorder, "factorization/postorder");
  CheckParameter(out, &common_.dbound, "factorization/dbound");
  CheckParameter(out, &analysis, &common_.supernodal, "supernodal/method");
  CheckParameter(out, &common_.supernodal_switch, "supernodal/switch");
  CheckParameter(out, &common_.print, "output/level");
  CheckParameter(out, &common_.precise, "output/precise");
  ParamNodeList mList = out->GetList("method");
  for(unsigned int i = 0; i < mList.GetSize(); ++i){
    PtrParamNode m = static_cast<PtrParamNode>(mList[i]);
    PtrParamNode pn = m->Get("id");
    int idx = pn->As<Integer>();
    CheckParameter(m, &common_.method[idx].prune_dense,"prune_dense");
    CheckParameter(m, &common_.method[idx].prune_dense2, "prune_dense2");
    CheckParameter(m, &common_.method[idx].nd_oksep, "nd_oksep");
    CheckParameter(m, &common_.method[idx].nd_small, "nd_small");
    CheckParameter(m, &common_.method[idx].aggressive, "aggressive");
    CheckParameter(m, &common_.method[idx].nd_compress, "nd_compress");
    CheckParameter(m, &common_.method[idx].nd_camd, "nd_camd");
    CheckParameter(m, &common_.method[idx].nd_components, "nd_components");
    CheckParameter(m, &common_.method[idx].ordering, "ordering");
  }
}

template<typename T>
void CholMod<T>::CholModDefaultParams()
{
  cholmod_defaults(&common_);
  // we only want to analyze/factorize/solve, not modify the factorization, so these parameters are recommended
  common_.grow0 = 0.0;
  common_.grow1 = 0.0;
  common_.grow2 = 0;
}

template<typename T>
void CholMod<T>::CholModAnalyze()
{
  fact_ = cholmod_analyze(mat_, &common_);
  TestException("analyze");
}

template<typename T>
void CholMod<T>::CholModFactorize()
{
  cholmod_factorize(mat_, fact_, &common_);
  TestException("factorize");
}

template<typename T>
void CholMod<T>::CholModFreeFactorization()
{
  cholmod_free_factor(&fact_, &common_);
  fact_ = NULL;
}

template<typename T>
void CholMod<T>::TestException(std::string where)
{
  if(common_.status == CHOLMOD_OK){
    return; // no error
  }
  
  // get error message
  std::stringstream ss;
  ss << "Error CholMod in " << where << ": ";  
  switch(common_.status){
  case CHOLMOD_NOT_INSTALLED:
    ss << "method not installed";
    break;
  case CHOLMOD_OUT_OF_MEMORY:
    ss << "out of memory";
    break;
  case CHOLMOD_TOO_LARGE:
    ss << "integer overflow occured";
    break;
  case CHOLMOD_INVALID:
    ss << "invalid input";
    break;
  case CHOLMOD_NOT_POSDEF:
    ss << "matrix not pos. def.";
    break;
  case CHOLMOD_DSMALL:
    ss << "D for LDL' or diag(L) or LL' has tiny absolute value"; // this is usually only a warning however it should never occur in real problems (matrix is close to not p.d.)
    break;
  default:
    ss << "unknown error";
    break;
  }
  throw Exception(ss.str());
}

template<typename T>
void CholMod<T>::SetEnums()
{
  analysis.SetName("CholMod::Analysis");
  analysis.Add(ANALYSIS_SIMPLICIAL, "simplicial");
  analysis.Add(ANALYSIS_AUTO, "auto");
  analysis.Add(ANALYSIS_SUPERNODAL, "supernodal");
}

// Explicit template instantiation
template class CholMod<Double> ;
template class CholMod<Complex> ;

}
