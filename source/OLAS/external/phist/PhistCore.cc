/*
 * PhistCore.cc
 *
 *  Created on: Mar 5, 2018
 *      Author: sri
 */

#include "PhistCore.hh"



namespace CoupledField {

DECLARE_LOG(phico)
DEFINE_LOG(phico, "phistCore")


PhistCore::PhistCore()
{
  int iflag;
  phist_comm_create(&comm_, &iflag);
  assert(iflag == 0);
  // dummy mpi arguments
  int argc = 0;
  char* v = NULL;
  char** argv = &v;
  phist_kernels_init(&argc, &argv, &iflag);
  assert(iflag== 0);

}

PhistCore::~PhistCore() {

}

/** Most used methods from SCRS_Matrix and CRS_Matrix have no common base :(. As templated methods cannot be
  * virtual, this would also be possible for some. Anyway, we do ugly copy&paste here :( */
template<class TYPE>
int PhistCore::SparseMatRowFunc(ghost_gidx row, ghost_lidx* row_nnz, ghost_gidx* row_col, void* values, void* service_void)
{
  const PhistCore::SparseMatRowFuncService* service = (const PhistCore::SparseMatRowFuncService*) service_void;

  assert(service->mat->GetStorageType() == BaseMatrix::SPARSE_SYM || service->mat->GetStorageType() == BaseMatrix::SPARSE_NONSYM);

  assert(!(boost::is_complex<TYPE>::value && service->mat->GetEntryType() != BaseMatrix::COMPLEX));
  const SparseOLASMatrix<TYPE>* mat   = dynamic_cast<const SparseOLASMatrix<TYPE>*>(service->mat);
  assert(mat != NULL);

  assert(row >= 0);
  assert(row < mat->GetNumRows());

  *row_nnz = mat->GetRowSize(row);

  TYPE* data = (TYPE*) values;

  //StdVector<int> cols(*row_nnz); // KILLME only for debug output
  unsigned int base = mat->GetRowPointer()[row];
  for(int i = 0; i < *row_nnz; i++) {
    row_col[i] = mat->GetColPointer()[base + i];
    //cols[i] = row_col[i];
  }

  for(int i = 0; i < *row_nnz; i++)
    data[i] = service->scale * mat->GetDataPointer()[base + i];

  LOG_DBG2(phico) << "SSMRF row=" << row << " row_nnz=" << *row_nnz << " scale=" << service->scale << " values=" << ToString<double>((double*) values, *row_nnz);
  // LOG_DBG2(phico) << "SSMRF row=" << row << "  row_col=" << cols.ToString();
  return 0; // all ok
}


template<class TYPE>
VsparseMat_t* PhistCore::InitMatrix(const BaseMatrix& cfs, VsparseMat_t** phist, double scale)
{
  // create stiffness or mass matrix for phist - which will be ghost matrices
  assert(!(boost::is_complex<TYPE>::value && cfs.GetEntryType() != BaseMatrix::COMPLEX));
  const SparseOLASMatrix<TYPE>* mat =  dynamic_cast<const SparseOLASMatrix<TYPE>*>(&cfs);
  assert(mat != NULL);
  assert(mat->GetStorageType() == StdMatrix::SPARSE_SYM || mat->GetStorageType() == StdMatrix::SPARSE_NONSYM);

  assert(mat != NULL);

  assert(mat->GetNumRows() == mat->GetNumCols());


  SparseMatRowFuncService service;
  service.mat = dynamic_cast<const StdMatrix*>(&cfs);
  service.scale = scale;

  int iflag = 0;

  // see https://stackoverflow.com/questions/610245/where-and-why-do-i-have-to-put-the-template-and-typename-keywords
  typename phist::types<TYPE>::sparseMat_ptr smp = (typename phist::types<TYPE>::sparseMat_ptr) *phist;

  // if already initialized, e.g. because of optimization, we delete forst
  if(smp != NULL)
    phist::kernels<TYPE>::sparseMat_delete(smp, &iflag);
  assert(iflag == 0);

  phist::kernels<TYPE>::sparseMat_create_fromRowFunc(&smp, comm_, mat->GetNumRows(), mat->GetNumCols(), mat->GetMaxRowSize(), SparseMatRowFunc<TYPE>, (void*) &service, &iflag);
  assert(iflag == 0);

  assert(smp != NULL);
  return (VsparseMat_t*) smp;
}

// template instantiation stuff
template VsparseMat_t* PhistCore::InitMatrix<double>(const BaseMatrix& cfs, VsparseMat_t** phist, double scale);
template VsparseMat_t* PhistCore::InitMatrix<Complex>(const BaseMatrix& cfs, VsparseMat_t** phist, double scale);

template int PhistCore::SparseMatRowFunc<double>(ghost_gidx row, ghost_lidx* row_nnz, ghost_gidx* row_col, void* values, void* service_void);
template int PhistCore::SparseMatRowFunc<Complex>(ghost_gidx row, ghost_lidx* row_nnz, ghost_gidx* row_col, void* values, void* service_void);


} /* namespace CoupledField */
