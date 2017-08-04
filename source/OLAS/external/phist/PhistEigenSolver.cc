#include <limits>
#include <string.h>

#include "MatVec/StdMatrix.hh"
#include "MatVec/CRS_Matrix.hh"
#include "MatVec/generatematvec.hh"
#include "Utils/Timer.hh"
#include "Domain/Domain.hh"
#include "Driver/BaseDriver.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "OLAS/algsys/SolStrategy.hh"

#include "OLAS/precond/generateprecond.hh"
#include "OLAS/precond/BasePrecond.hh"
#include "OLAS/solver/generatesolver.hh"
#include "OLAS/solver/BaseSolver.hh"

#include "PhistEigenSolver.hh"

#include <mpi.h>
#include <phist_kernels.h>
#include "phist_subspacejada.h"
#include "phist_schur_decomp.h"
#include <phist_gen_d.h>
#include "phist_driver_utils.h"


DECLARE_LOG(pes)
DEFINE_LOG(pes, "phistEigenSolver")

namespace CoupledField {

  PhistEigenSolver::PhistEigenSolver( shared_ptr<SolStrategy> strat,
                                        PtrParamNode xml,
                                        PtrParamNode solverList,
                                        PtrParamNode precondList,
                                        PtrParamNode eigenInfo)
    : BaseEigenSolver( strat, xml, solverList, precondList, eigenInfo )
  {
    A_   = NULL;
    B_   = NULL;
    xml_ = xml;
  }

  PhistEigenSolver::~PhistEigenSolver()
  {
    // TODO destroy A_ and B_!!
    // ghost_sparsemat_destroy(A_);
    A_ = NULL;
    B_ = NULL;
    //ghost_densemat_destroy(x_);
    //ghost_densemat_destroy(y_);


    ghost_finalize();


  }
  void PhistEigenSolver::Setup(const BaseMatrix & mat, UInt numFreq, Double freqShift, bool sort)
  {
    LOG_DBG(pes) << "PES:S(stiff)";
    // Set flag for indicating a non-quadratic problem
    isQuadratic_ = false;
    isBloch_ = false;
    sort_ = false;

    // Save frequency parameters
    numFreq_ = numFreq;
    freqShift_ = freqShift;

    ToInfo();
  }

  
  /** implementation of phist_sparseMat_rowFunc which is called for every row to copy the sparse matrix data
   * @param row global row index (int64 by default)
   * @param nnz of row (output, int)
   * @param col - column indices of row (output)
   * @param values values of row (output)
   * @param service - cfs internal use: pointer to cfs matrix */
  int SparseMatRowFunc(ghost_gidx row, ghost_lidx* row_nnz, ghost_gidx* row_col, void* values, void* service)
  {
    const CRS_Matrix<double>* mat = (const CRS_Matrix<double>*) service;
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
      data[i] = mat->GetDataPointer()[base + i];

    LOG_DBG2(pes) << "SMRF row=" << row << " row_nnz=" << *row_nnz << " row_col=" << cols.ToString() << " values=" << ToString<double>((double*) values, *row_nnz);

    return 0; // all ok
  }

  void PhistEigenSolver::Setup(const BaseMatrix& stiffMat, const BaseMatrix& massMat, UInt numFreq, Double freqShift, bool sort, bool bloch)
  {
    LOG_DBG(pes) << "PES:S(stiff, mass)";
    // Set flag for indicating a non-quadratic problem
    isQuadratic_ = false;
    isBloch_ = bloch;
    sort_ = sort;

    // Save frequency parameters
    numFreq_ = numFreq;
    freqShift_ = freqShift;

    // MPI handle which is by default MPI_COMM_WORLD
    phist_comm_ptr comm = NULL;
    int err;
    phist_comm_create(&comm, &err);
    assert(err == 0);
    LOG_DBG(pes) << "phist_comm_create -> " << comm;
    MPI_Init(NULL, NULL);
    int argc = 0;
    char* v = NULL;
    char** argv = &v;
    phist_kernels_init(&argc, &argv, &err);
    assert(err == 0);
    LOG_DBG(pes) << "phist_kernels_init -> " << err;

    // phist parameters usually read from opts-standard.txt

    // fill the jadaOpts struct to pass settings to the solver
    phist_jadaOpts opts;
    phist_jadaOpts_setDefaults(&opts);

    opts.symmetry = phist_GENERAL; // phist_HERMITIAN also for real symmetric
    opts.numEigs = numFreq;
    opts.which = phist_SR; // smalles real part
    opts.convTol = 1e-8;
    opts.blockSize = 4; // usually 1,2 or 4.;
    opts.maxIters = 150;
    opts.minBas = 28; //numFreq + 2 * opts.blockSize;
    opts.maxBas = 60; //opts.minBas + 10 * opts.blockSize;
    // parameters for the linear system stuff

    opts.innerSolvBlockSize = opts.blockSize;
    opts.innerSolvType = phist_MINRES;
    opts.innerSolvMaxBas = 10;
    opts.innerSolvMaxIters = 10;
    opts.innerSolvRobust = 1;

    // for include stuff issues, A_ and B_ attributes are not of full type
    phist_DsparseMat* A = (phist_DsparseMat*) A_;
    phist_DsparseMat* B = (phist_DsparseMat*) B_;

    // create stiffness and mass matrix for phist - which will be ghost matrices
    const CRS_Matrix<double>* stiff =  dynamic_cast<const CRS_Matrix<double>*>(&stiffMat);
    assert(stiff != NULL);
    phist_lidx max_nne = stiff->GetMaxRowSize();
    LOG_DBG(pes) << " max row nnz=" << max_nne;
    LOG_DBG2(pes) << " row=" << ToString<unsigned int>(stiff->GetRowPointer(), stiff->GetNumRows()+1);
    LOG_DBG2(pes) << " col=" << ToString<unsigned int>(stiff->GetColPointer(), stiff->GetNnz());
    LOG_DBG2(pes) << " val=" << ToString<double>(stiff->GetDataPointer(), stiff->GetNnz());
    assert(stiff->GetNumRows() == stiff->GetNumCols());
    phist_DsparseMat_create_fromRowFunc(&A, comm, stiff->GetNumRows(), stiff->GetNumCols(), max_nne, SparseMatRowFunc, (void*) stiff, &err);
    assert(err == 0);
    LOG_DBG(pes) << "create A -> " << err;

    const CRS_Matrix<double>* mass =  dynamic_cast<const CRS_Matrix<double>*>(&massMat);
    assert(mass != NULL);
    phist_DsparseMat_create_fromRowFunc(&B, comm, mass->GetNumRows(), mass->GetNumCols(), mass->GetMaxRowSize(), SparseMatRowFunc, (void*) mass, &err);
    assert(err == 0);
    LOG_DBG(pes) << "create B -> " << err;

    // setup eigenvalue problem - taken from subspacejada (jacobi davidson)
    // create an operator from A
    phist_DlinearOp_ptr opA = new phist_DlinearOp();

    // we need the domain map of the matrix
    phist_const_map_ptr map = NULL;
    phist_DsparseMat_get_domain_map(A,&map,&err);
    LOG_DBG(pes) << "phist_DsparseMat_get_domain_map -> " << err;
    assert(err == 0);

    phist_DlinearOp_ptr opB = new phist_DlinearOp();
    phist_DlinearOp_wrap_sparseMat(opB,B,&err);
    LOG_DBG(pes) << "phist_DlinearOp_wrap_sparseMat -> " << err;
    assert(err == 0);

    phist_DlinearOp_wrap_sparseMat_pair(opA,A,B,&err);
    LOG_DBG(pes) << "phist_DlinearOp_wrap_sparseMat_pair -> " << err;
    assert(err == 0);

    int prob_size = numFreq+opts.blockSize-1;
    LOG_DBG(pes) << "prob_size -> " << prob_size;

    // setup necessary vectors and matrices for the schur form
    mvec_ptr Q = NULL;
    phist_Dmvec_create(&Q,map,prob_size,&err);
    LOG_DBG(pes) << "phist_Dmvec_create -> " << err;
    assert(err == 0);

    sdMat_ptr R = NULL;
    phist_DsdMat_create(&R,prob_size,prob_size,comm,&err);
    LOG_DBG(pes) << "phist_DsdMat_create -> " << err;
    assert(err == 0);

    StdVector<double> resNorm(prob_size);
    StdVector<std::complex<double> > ev(prob_size); // always complex,

    // setup start vector (currently to (1 0 1 0 .. ) )
    mvec_ptr v0 = NULL;
    phist_Dmvec_create(&v0,map,1,&err);
    LOG_DBG(pes) << "phist_Dmvec_create -> " << err;
    assert(err == 0);

    phist_Dmvec_put_value(v0,1.0,&err);
    assert(err == 0);

    // skip residual calculation from subspacejada, not necessary for us
    opts.v0=v0;

    // assume phist_NO_PRECON

    int nIter=opts.maxIters;
    int nEig = (int) numFreq;

    // The actual calculation !!! - standard Ritz values for exterior eigenvalues.
    // QR-factorization, Q=orthogonal ev basis , R=upper triangular matrix, ev are on the diagonal,
    // A*Q = B*Q*R holds
    assert(opts.how != phist_HARMONIC);

    phist_Dsubspacejada(opA, opB, opts, Q, R, ev.GetPointer(), resNorm.GetPointer(), &nEig, &nIter, &err);
    LOG_DBG(pes) << "phist_Dsubspacejada -> nEig=" << nEig << " nIter=" << nIter << " ev=" << ev.ToString() << " -> " << err;
    assert(err >= 0);

    // skip calculation of real residual, res = AQ - BQR

    // compute eigenvectors X: A*X=X*D and diagonal matrix D with eigenvalues
    // (for checking the sorting only).
    mvec_ptr X=NULL;
    phist_Dmvec_create(&X,map,nEig,&err);
    assert(err >= 0);

    phist_DComputeEigenvectors(Q,R,X,&err);
    assert(err >= 0);

    // skip calculation of eigenvector residuals
    
    // download X from GPU if applicable
    phist_Dmvec_from_device(X,&err);
    assert(err == 0);

    // get pointer to row/col major block of vectors
    double* xval = NULL; // will be set to memory owned by phist
    phist_lidx lda, nloc;
    phist_Dmvec_my_length(X,&nloc,&err);
    assert(err == 0);
    phist_Dmvec_extract_view(X,&xval,&lda,&err);
    assert(err == 0);
    
    // this is how one can extract the eigenvectors.
    // TODO: move all of this to CalcEigenFrequencies and other functions below
/*
    for (int j=0; j<nEig; j++)
    {
      for (int i=0; i<nloc; i++)
      {
#ifdef PHIST_MVECS_ROW_MAJOR
         mode_j[i]=xval[i*lda+j];
#else
         mode_j[i]=xval[j*lda+i];
#endif
    }
*/
    ToInfo();
  }


  void PhistEigenSolver::Setup(const BaseMatrix & stiffMat, const BaseMatrix & massMat, const BaseMatrix & dampMat,
                                UInt numFreq, double freqShift, bool sort)
  {
    LOG_DBG(pes) << "PES:S(stiff, mass, damp)";

    // Set flag for indicating a non-quadratic problem
    isQuadratic_ = true;
    isBloch_ = false;
    sort_ = false;
  }

  void PhistEigenSolver::ToInfo()
  {
    PtrParamNode setup = info_->Get(ParamNode::HEADER);
    setup->Get("quadratic")->SetValue(isQuadratic_);
    setup->Get("bloch")->SetValue(isBloch_);
  }




  void PhistEigenSolver::CalcEigenFrequencies(BaseVector &sol, BaseVector &err)
  {
    assert(false);
  }

  void PhistEigenSolver::CalcConditionNumber(const BaseMatrix& mat, Double& condNumber, Vector<Double>& evs, Vector<Double>& err )
  {
    // Set flag for indicating a non-quadratic problem
    isQuadratic_ = false;
    // NOTE: Hard coded as true!!!

  }

  void PhistEigenSolver::GetEigenMode(UInt modeNr, Vector<Complex> & mode)
  {

  }

  void PhistEigenSolver::GetComplexEigenMode(UInt modeNr, Vector<Complex>& mode)
  {
    // in bloch mode case the same as GetEigenMode,
    // in quadratic case the modes have internally double size and we want the upper half
  }

}
