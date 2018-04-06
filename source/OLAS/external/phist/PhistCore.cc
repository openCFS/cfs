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


PhistCore::PhistCore() {

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
  // TODO Auto-generated destructor stub
}

/** Most used methods from SCRS_Matrix and CRS_Matrix have no common base :(. As templated methods cannot be
  * virtual, this would also be possible for some. Anyway, we do ugly copy&paste here :( */
 int SparseMatRowFunc(ghost_gidx row, ghost_lidx* row_nnz, ghost_gidx* row_col, void* values, void* service_void)
 {
   const PhistCore::SparseMatRowFuncService* service = (const PhistCore::SparseMatRowFuncService*) service_void;

   assert(service->mat->GetStorageType() == BaseMatrix::SPARSE_SYM || service->mat->GetStorageType() == BaseMatrix::SPARSE_NONSYM);

   const SparseOLASMatrix<double>* mat   = dynamic_cast<const SparseOLASMatrix<double>*>(service->mat);
   assert(mat != NULL);

   assert(row >= 0);
   assert(row < mat->GetNumRows());

   *row_nnz = mat->GetRowSize(row);

   double* data = (double*) values;

   StdVector<int> cols(*row_nnz); // KILLME only for debug output
   unsigned int base = mat->GetRowPointer()[row];
   for(int i = 0; i < *row_nnz; i++) {
     row_col[i] = mat->GetColPointer()[base + i];
     cols[i] = row_col[i];
   }

   for(int i = 0; i < *row_nnz; i++)
     data[i] = service->scale * mat->GetDataPointer()[base + i];

   LOG_DBG2(phico) << "SSMRF row=" << row << " row_nnz=" << *row_nnz << " scale=" << service->scale << " row_col=" << cols.ToString() << " values=" << ToString<double>((double*) values, *row_nnz);

   return 0; // all ok
 }


sparseMat_t* PhistCore::InitMatrix(const BaseMatrix& cfs, sparseMat_t** phist, double scale)
 {
   // create stiffness or mass matrix for phist - which will be ghost matrices
   const SparseOLASMatrix<double>* mat =  dynamic_cast<const SparseOLASMatrix<double>*>(&cfs);
   assert(mat != NULL);
   assert(mat->GetStorageType() == StdMatrix::SPARSE_SYM || mat->GetStorageType() == StdMatrix::SPARSE_NONSYM);

   assert(mat != NULL);

   assert(mat->GetNumRows() == mat->GetNumCols());


   SparseMatRowFuncService service;
   service.mat = dynamic_cast<const StdMatrix*>(&cfs);
   service.scale = scale;

   phist::types<double>::sparseMat_ptr smp = (phist::types<double>::sparseMat_ptr) *phist;
   assert(smp == NULL); // we are prior initialization
   int iflag = 0;
   phist::kernels<double>::sparseMat_create_fromRowFunc(&smp, comm_, mat->GetNumRows(), mat->GetNumCols(), mat->GetMaxRowSize(), SparseMatRowFunc, (void*) &service, &iflag);
   assert(iflag == 0);

   assert(smp != NULL);
   return (sparseMat_t*) smp;
 }




} /* namespace CoupledField */
