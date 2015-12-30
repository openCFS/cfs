#include <string.h>
// you might want to alternatively allow the academic FEAST package instead the implementation from MKL
#include <mkl.h>

#include "MatVec/StdMatrix.hh"
#include "MatVec/SCRS_Matrix.hh"

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/Logging/LogConfigurator.hh"

#include "OLAS/external/feast/FeastEigenSolver.hh"

DECLARE_LOG(fes)
DEFINE_LOG(fes, "feast")

namespace CoupledField {


FeastEigenSolver::FeastEigenSolver(shared_ptr<SolStrategy> strat, PtrParamNode xml, PtrParamNode solverList,
                                   PtrParamNode precondList, PtrParamNode eigenInfo)
  : BaseEigenSolver(strat, xml, solverList, precondList, eigenInfo)
{
  this->generalized_ = false;
  this->bloch_ = false;
  this->a_ = NULL;
  this->b_ = NULL;
  this->numFreq_ = 0;
  this->freqShift_ = 0;
  this->fpm_.Resize(128, 0);
  this->n_ = 0;
  this->m0_ = 0;
  this->m_ = 0;
  this->info_ = -10; // not im MKL code (page 1635)
}

FeastEigenSolver::~FeastEigenSolver()
{

}

void FeastEigenSolver::Setup(const BaseMatrix& stiffMat, unsigned int numFreq, double freqShift)
{
  numFreq_ = numFreq;
  freqShift_ = freqShift;
  generalized_ = false;
  bloch_ = false;

  // initialize feastinit fpm according to "Extended Eigensolver Input Parameters" in the MKL manual
  feastinit(fpm_.GetPointer());
  fpm_[0] = 1; // runtime messages to the screen
  fpm_[1] = 3;
  fpm_[3] = 50; // max refinement loops
  fpm_[4] = 0;  // no initial subspace
  fpm_[6] = 10; // error trace single precission stopping criteria
  fpm_[27] = 1; // check b for positive definiteness

  n_ = stiffMat.GetNumRows();
  m0_ = numFreq_;
  assert(m0_ <= n_);
  x_.Resize(m0_ * 20000 * n_, 0.0); // FIXME

  // store the matrix
  a_ = &dynamic_cast<const StdMatrix&>(stiffMat);
  assert(a_ != NULL);
  // copy the indices and transform to 1-based
  a_->ExportCRSColumns(ja_, 1); // 1 for 1-based
  a_->ExportCRSRows(ia_, 1, true); // tailing nnz
}

void FeastEigenSolver::Setup( const BaseMatrix& stiffMat, const BaseMatrix& massMat,
            unsigned int numFreq, double freqShift, bool bloch)
{
  Setup(stiffMat, numFreq, freqShift);

  generalized_ = true;
  bloch_ = bloch;

  b_ = &dynamic_cast<const StdMatrix&>(massMat);

  // store the matrix
  b_ = &dynamic_cast<const StdMatrix&>(stiffMat);
  assert(b_ != NULL);
  // copy the indices and transform to 1-based
  b_->ExportCRSColumns(jb_, 1); // 1 for 1-based
  b_->ExportCRSRows(ib_, 1, true); // tailing nnz

  assert(a_->GetNumRows() == b_->GetNumRows());
  assert((int) b_->GetNumRows() == n_);
  assert(a_->GetStorageType() == b_->GetStorageType());
  assert(a_->GetStructureType() == b_->GetStructureType());
  assert(a_->GetEntryType() == b_->GetEntryType());
}

void FeastEigenSolver::Setup( const BaseMatrix& stiffMat,
            const BaseMatrix& massMat,
            const BaseMatrix& dampMat,
            unsigned int numFreq, double freqShift )
{
  assert(false);
}
int info;

void FeastEigenSolver::CalcEigenFrequencies(BaseVector& bvs, BaseVector& err)
{
  assert(a_ != NULL);
  assert(a_->GetStructureType() == BaseMatrix::SPARSE_MATRIX);

  switch(a_->GetStorageType())
  {
  case BaseMatrix::SPARSE_SYM:
  {
    assert(a_->GetEntryType() == BaseMatrix::DOUBLE);
    Vector<double>& e = dynamic_cast<Vector<double>&>(bvs);
    Vector<double>& res = dynamic_cast<Vector<double>&>(err);

    // standard problem: void dfeast_scsrev (const char * uplo, const MKL_INT * n, const double * a, const
    // MKL_INT * ia, const MKL_INT * ja, MKL_INT * fpm, double * epsout, MKL_INT * loop,
    // const double * emin, const double * emax, MKL_INT * m0, double * e, double * x,
    // MKL_INT * m, double * res, MKL_INT * info);

    assert(a_->GetNumRows() == a_->GetNumCols());
    const SCRS_Matrix<double>* a = dynamic_cast<const SCRS_Matrix<double>*>(a_);
    // input parameters
    char uplo = 'U'; // according to CFS doku U
    assert((int) ia_.GetSize() == n_ + 1);
    double emin = 0.0;
    double emax = 10000;
    assert(m0_ >= numFreq_);

    // output parameters
    double epsout;
    int loop;
    e.Resize(m0_); // e: m found eigenvalues
    assert((int) x_.Capacity() >= n_ * m0_);
    res.Resize(m0_); // the first m get the relative residuals

    if(!generalized_)
    {
      dfeast_scsrev(&uplo, &n_, a->GetDataPointer(), ia_.GetPointer(), ja_.GetPointer(),
                    fpm_.GetPointer(), &epsout, &loop, &emin, &emax, &m0_,
                    e.GetPointer(), x_.GetPointer(), &m_, res.GetPointer(), &info_);
    }
    else
    {
      // generalized problem:  void dfeast_scsrgv (const char * uplo, const MKL_INT * n, const double * a, const
      // MKL_INT * ia, const MKL_INT * ja, const double * b, const MKL_INT * ib, const MKL_INT *
      // jb, MKL_INT * fpm, double * epsout, MKL_INT * loop, const double * emin, const double
      // * emax, MKL_INT * m0, double * e, double * x, MKL_INT * m, double * res, MKL_INT * info);

      const SCRS_Matrix<double>* b = dynamic_cast<const SCRS_Matrix<double>*>(b_);
      assert((int) ib_.GetSize() == n_ + 1);

      //LOG_DBG3(fes) << "CEF a: " << StdVector<double>::ToString(a->GetNumEntries(), a->GetDataPointer(), 0);
      //LOG_DBG3(fes) << "CEF ia: " << ia_.ToString();
      //LOG_DBG3(fes) << "CEF ja: " << ja_.ToString();

      assert(b->GetNumEntries() < b->GetNnz()); // for symmetric matrices!
      LOG_DBG3(fes) << "CEF n: " << n_;
      LOG_DBG3(fes) << "CEF b_size: " << b->GetNumEntries();
      // LOG_DBG3(fes) << "CEF b: " << StdVector<double>::ToString(b->GetNumEntries(), b->GetDataPointer(), 0);
      // LOG_DBG3(fes) << "CEF ib: " << ib_.ToString();
      // LOG_DBG3(fes) << "CEF jb: " << jb_.ToString();

      dfeast_scsrgv(&uplo, &n_, a->GetDataPointer(), ia_.GetPointer(), ja_.GetPointer(), b->GetDataPointer(), ib_.GetPointer(), jb_.GetPointer(),
                    fpm_.GetPointer(), &epsout, &loop, &emin, &emax, &m0_,
                    e.GetPointer(), x_.GetPointer(), &m_, res.GetPointer(), &info_);

      LOG_DBG(fes) << "CEF info -> " << info_;
      LOG_DBG(fes) << "CEF loop -> " << loop;
      LOG_DBG(fes) << "CEF m -> " << m_;
      LOG_DBG(fes) << "CEF epsout -> " << epsout;
      e.Resize(m_);
      LOG_DBG(fes) << "CEF e -> " << e.ToString();
      res.Resize(m_);
      LOG_DBG(fes) << "CEF res -> " << res.ToString();

    }


    break;
  }
  default:
    assert(false);
  }

}


void FeastEigenSolver::GetEigenMode(unsigned int modeNr, Vector<Complex>& mode)
{
  mode.Resize(n_, 0.0);
  return;

  assert((int) modeNr <= m_); // assume 0-based!!
  assert(x_.GetSize() >= modeNr * n_);
  mode.Resize(n_);
  for(unsigned int i = 0; i < (unsigned int) n_; i++)
    mode[i] = x_[modeNr * n_ + i];
}

void FeastEigenSolver::GetComplexEigenMode(unsigned int modeNr, Vector<Complex>& mode)
{
  assert(false);
}

void FeastEigenSolver::CalcConditionNumber( const BaseMatrix& mat, double& condNumber, Vector<double>& evs, Vector<double>& err)
{
  assert(false);
}

void FeastEigenSolver::ToInfo()
{

}


} // end of namespace
