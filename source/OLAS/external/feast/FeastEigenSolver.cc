#include <string.h>
#include <def_use_feast.hh>
#ifdef USE_FEAST_COMMUNITY
  // include the needed community FEAST libraries (which are C)
  extern "C"{
    #include "feast.h"
    #include "feast_sparse.h"
  }
#else
  #include <mkl_solvers_ee.h>
#endif

#include "MatVec/StdMatrix.hh"
#include "MatVec/SCRS_Matrix.hh"
#include "MatVec/CRS_Matrix.hh"

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
  xml_ = xml;
}

FeastEigenSolver::~FeastEigenSolver()
{

}


void FeastEigenSolver::Setup(const BaseMatrix & A, bool isHermitian){
    // check for symmetry
    this->CheckMatrix(isReal_, isSymmetric_, A);
    this->SetProblemType(A,isHermitian);
    // initialize feastinit fpm according to "Extended Eigensolver Input Parameters" in the MKL manual
    // see: https://software.intel.com/en-us/mkl-developer-reference-c-extended-eigensolver-input-parameters
    //int feastparam[64];
    feastinit(fpm_.GetPointer());
    // now get values from the <feast> section of the XML
    bool logging = false; xml_->GetValue("logging", logging, ParamNode::INSERT); if (logging == true) fpm_[0] = 1; // runtime messages to the screen
    uint Ne = 8; xml_->GetValue("Ne", Ne, ParamNode::INSERT); fpm_[1] = Ne; // Number of contour points
    uint stopCrit = 12; xml_->GetValue("stopCrit", stopCrit, ParamNode::INSERT); fpm_[2] = stopCrit; // error trace double precision stopping criteria eps=10^-fpm[2]
    uint maxRefinementLoops = 20; xml_->GetValue("maxRefinementLoops", maxRefinementLoops, ParamNode::INSERT); fpm_[3] = maxRefinementLoops; // max refinement loops
    fpm_[4] = 0;  // no initial subspace
    fpm_[26] = 1; // Extended Eigensolver routines check input matrices
    fpm_[27] = 1; // check b for positive definiteness

    n_ = A.GetNumRows(); // size of the problem
    uint m0 = 0; xml_->GetValue("m0", m0, ParamNode::INSERT); if (m0>0) m0_ = m0; else m0_=n_; // initial subspace size
    assert(m0_ <= n_);
    vr_.Resize(m0_ * n_, 0.0);

    // store the matrix

    // a
    // Array containing the nonzero elements of either the full matrix A or the upper or lower
    // triangular part of the matrix A, as specified by uplo.
    a_ = &dynamic_cast<const StdMatrix&>(A);
    assert(a_ != NULL);
    // copy the indices and transform to 1-based

    // ia
    // Array of length n + 1, containing indices of elements in the array a ,
    // such that ia[i - 1] is the index in the array a of the first non-zero element from the row i .
    // The value of the last element ia[n] is equal to the number of non-zeros plus one
    a_->ExportCRSRows(ia_, 1, true); // tailing nnz = number of non-zeros
    ia_[ia_.GetSize()-1] = ia_[ia_.GetSize()-1] + 1; // increase by one, see above

    // ja
    // Array containing the column indices for each non-zero element of the matrix A being represented in the array a.
    // Its length is equal to the length of the array a.
    a_->ExportCRSColumns(ja_, 1); // 1 for 1-based
}

void FeastEigenSolver::Setup(const BaseMatrix & A, const BaseMatrix & B, bool isHermitian){
    generalized_=true;
    // first set A as for standard EVP
    Setup(A, isHermitian);
    // check the B matrix
    bool isReal, isSymmetric;
    this->CheckMatrix(isReal, isSymmetric, B);
    this->SetProblemType(B,isHermitian);
    if (!isReal || !isReal_) isReal_ = false;
    if (!isSymmetric || !isSymmetric_) isSymmetric = false;
    // store matrix
    b_ = &dynamic_cast<const StdMatrix&>(B);
    assert(b_ != NULL);
    // copy the indices and transform to 1-based
    b_->ExportCRSColumns(jb_, 1); // 1 for 1-based
    b_->ExportCRSRows(ib_, 1, true); // nzz = tailing number of non-zeros
    ib_[ib_.GetSize()-1] = ib_[ib_.GetSize()-1] + 1; // according to doc this must be nzz+1
    // check if size is equal
    assert(a_->GetNumRows() == b_->GetNumRows());
    assert((int) b_->GetNumRows() == n_);
    assert(a_->GetStorageType() == b_->GetStorageType());
    assert(a_->GetStructureType() == b_->GetStructureType());
    assert(a_->GetEntryType() == b_->GetEntryType());
}

void FeastEigenSolver::CalcEigenValues(BaseVector& sol, BaseVector& err, Double minVal, Double maxVal){
    assert(a_ != NULL);
    assert(a_->GetStructureType() == BaseMatrix::SPARSE_MATRIX);
    if (a_->GetStructureType() != BaseMatrix::SPARSE_MATRIX) {
        EXCEPTION("FEAST is currently only implemented for sparse matrices.")
    }
    assert((int) ia_.GetSize() == n_ + 1);

    // set some common input paramters
    char uplo = 'U'; // according to CFS doku U
    // set List-I parameters depending on problem type
    // I1, I2 = { Emid, r } : complex-general, complex-symmetric, real-general
    // I1, I2 = { Emin, Emax } : complex-hermitian, real-symmetric
    double I1 [2];
    double I2;

    int loop; // # of feast subspace iteration
    switch (eigenProblemType_) {

    case REAL_GENERAL : {
        LOG_DBG3(fes) << "real-general EVP";
        //Complex Emid = Complex()
        I1[0] = 0.5*(maxVal+minVal); I1[1]=0.0; LOG_DBG3(fes) << "CEF Emid: " << I1;
        I2 = 0.5*(maxVal-minVal); LOG_DBG3(fes) << "CEF r: " << I2;
        LOG_DBG3(fes) << "CEF m0: " << m0_; // specifies the initial guess for subspace dimension to be used
        LOG_DBG3(fes) << "CEF x:" << vr_.ToString();
        LOG_DBG3(fes) << "CEF n: " << n_; // size of the problem
        // setup A-matrix
        const CRS_Matrix<double>* a_const = dynamic_cast<const CRS_Matrix<double>*>(a_);
        CRS_Matrix<double>* a = const_cast<CRS_Matrix<double>*>(a_const);
        //LOG_DBG3(fes) << "CEF a_size: " << a->GetNumEntries();
        //LOG_DBG3(fes) << "CEF a: " << StdVector<double>::ToString(a->GetNumEntries(), a->GetDataPointer(), 0);
        LOG_DBG3(fes) << "CEF ia: " << ia_.ToString();
        LOG_DBG3(fes) << "CEF ja: " << ja_.ToString();

        // output parameters
        double epsout;
        Vector<Complex>& E_c = dynamic_cast<Vector<Complex>&>(sol); // eigenvalues
        Vector<double> E; E.Resize(2*m0_,0.0); // factor 2 for complex
        Vector<Complex>& res_c = dynamic_cast<Vector<Complex>&>(err); // Relative residual
        Vector<double> res; res.Resize(2*m0_,0.0); // factor 2 for complex
        Vector<double> X; X.Resize(2*n_*m0_*2);// eigenvectors (if needed) //factor 2 for complex // factor 2 for L and R
        //assert((int) vr_.Capacity() >= 4 * n_ * m0_);
        if(!generalized_) // solve standard EVP
        {
            dfeast_gcsrev(&n_, a->GetDataPointer(), ia_.GetPointer(), ja_.GetPointer(),
                        fpm_.GetPointer(), &epsout, &loop, I1, &I2, &m0_,
                       E.GetPointer(), X.GetPointer(), &m_, res.GetPointer(), &info_);
        }
        else // generalized EVP
        {
            const CRS_Matrix<double>* bc = dynamic_cast<const CRS_Matrix<double>*>(b_);
            CRS_Matrix<double>* b = const_cast<CRS_Matrix<double>*>(bc);
            assert((int) ib_.GetSize() == n_ + 1);
            //assert(b->GetNumEntries() < b->GetNnz()); // for symmetric matrices!
            //LOG_DBG3(fes) << "CEF b_size: " << b->GetNumEntries();
            //LOG_DBG3(fes) << "CEF b: " << StdVector<double>::ToString(b->GetNumEntries(), b->GetDataPointer(), 0);
            LOG_DBG3(fes) << "CEF ib: " << ib_.ToString();
            LOG_DBG3(fes) << "CEF jb: " << jb_.ToString();
            dfeast_gcsrgv(&n_, a->GetDataPointer(), ia_.GetPointer(), ja_.GetPointer(), b->GetDataPointer(), ib_.GetPointer(), jb_.GetPointer(),
                          fpm_.GetPointer(), &epsout, &loop, I1, &I2, &m0_,
                          E.GetPointer(), X.GetPointer(), &m_, res.GetPointer(), &info_);
        }
        // put the solution into the complex vectors
        E_c.Resize(m_);
        res_c.Resize(m_);
        vr_.Resize(m_*n_);
        vl_.Resize(m_*n_);
        LOG_DBG3(fes) << "CEF # of found eigen values =  " << m_;
        LOG_DBG3(fes) << "CEF # elements in X: "<<X.GetSize();
        for (int i=0;i<m_;i++){
            E_c[i] = Complex( E[2*i],E[2*i+1] );
            res_c[i] = Complex( res[2*i], res[2*i+1] );
            for (int j=0;j<n_;j++){
                vr_[i*n_+j] = Complex( X[i*n_*2+2*j], X[i*n_*2+2*j+1]);
                vl_[i*n_+j] = Complex( X[2*n_*m_+i*n_*2+2*j], X[2*n_*m_+i*n_*2+2*j+1]);
            }
            //printf("   %d %.15e+i*%.15e %.15e\n",i,*(E+2*i),*(E+2*i+1),*(res+i));
        }
        LOG_DBG(fes) << "CEF epsout -> " << epsout;
    } break;;
    case COMPLEX_SYMMETRIC: {
        LOG_DBG3(fes) << "complex-symmetric EVP";
        //I1 = 0.5*(maxVal+minVal); LOG_DBG3(fes) << "CEF Emid: " << I1;
        I2 = 0.5*(maxVal-minVal); LOG_DBG3(fes) << "CEF r: " << I2;
        EXCEPTION("not yet implemented")
    } break;;
    case COMPLEX_GENERAL: {
        LOG_DBG3(fes) << "complex-general EVP";
        //I1 = 0.5*(maxVal+minVal); LOG_DBG3(fes) << "CEF Emid: " << I1;
        I2 = 0.5*(maxVal-minVal); LOG_DBG3(fes) << "CEF r: " << I2;
        EXCEPTION("not yet implemented")


    } break;;
    case COMPLEX_HERMITIAN: {
        LOG_DBG3(fes) << "complex-hermitian EVP (with symmetric matrix storage)";
        I1[0] = minVal; LOG_DBG3(fes) << "CEF Emin: " << I1[0];
        I2 = maxVal; LOG_DBG3(fes) << "CEF Emax: " << I2;
        LOG_DBG3(fes) << "CEF m0: " << m0_; // specifies the initial guess for subspace dimension to be used
        LOG_DBG3(fes) << "CEF x:" << vr_.ToString();
        LOG_DBG3(fes) << "CEF n: " << n_; // size of the problem
        // setup A-matrix
        const SCRS_Matrix<Complex>* a_const = dynamic_cast<const SCRS_Matrix<Complex>*>(a_);
        SCRS_Matrix<Complex>* a = const_cast<SCRS_Matrix<Complex>*>(a_const);
        Vector<double> aDouble; aDouble.Resize(2*a->GetNnz(),0.0); // factor 2 for complex
        for (int i=0;i<(int)a->GetNnz();i++) {
            aDouble[2*i] = a->GetDataPointer()[i].real();
            aDouble[2*i+1] = a->GetDataPointer()[i].imag();
        }
        LOG_DBG3(fes) << "CEF a_size: " << a->GetNumEntries();
        LOG_DBG3(fes) << "CEF a: " << StdVector<Complex>::ToString(a->GetNumEntries(), a->GetDataPointer(), 0);
        LOG_DBG3(fes) << "CEF ia: " << ia_.ToString();
        LOG_DBG3(fes) << "CEF ja: " << ja_.ToString();
        // output parameters
        double epsout;
        Vector<double>& E = dynamic_cast<Vector<double>&>(sol); // eigenvalues
        //Vector<double> E; E.Resize(m0_,0.0); // Eigenvalues (are real)
        Vector<double>& res = dynamic_cast<Vector<double>&>(err); // Relative residual
        E.Resize(m0_); // e: m found eigenvalues
        Vector<double> X; X.Resize(n_*m0_*2); // factor 2 for complex
        //assert((int) X.Capacity() >= n_ * m0_);
        res.Resize(m0_); // the first m get the relative residuals


        if(!generalized_) // solve standard EVP
        {
            zfeast_hcsrev(&uplo, &n_, aDouble.GetPointer(), ia_.GetPointer(), ja_.GetPointer(),
                        fpm_.GetPointer(), &epsout, &loop, I1, &I2, &m0_,
                        E.GetPointer(), X.GetPointer(), &m_, res.GetPointer(), &info_);
        }
        else // generalized EVP
        {
            const SCRS_Matrix<Complex>* bc = dynamic_cast<const SCRS_Matrix<Complex>*>(b_);
            SCRS_Matrix<Complex>* b = const_cast<SCRS_Matrix<Complex>*>(bc);
            Vector<double> bDouble; bDouble.Resize(2*b->GetNnz(),0.0); // factor 2 for complex
            for (int i=0;i<(int)b->GetNnz();i++) {
                bDouble[2*i] = b->GetDataPointer()[i].real();
                bDouble[2*i+1] = b->GetDataPointer()[i].imag();
            }
            assert((int) ib_.GetSize() == n_ + 1);
            assert(b->GetNumEntries() < b->GetNnz()); // for symmetric matrices!
            LOG_DBG3(fes) << "CEF b_size: " << b->GetNumEntries();
            LOG_DBG3(fes) << "CEF b: " << StdVector<Complex>::ToString(b->GetNumEntries(), b->GetDataPointer(), 0);
            LOG_DBG3(fes) << "CEF ib: " << ib_.ToString();
            LOG_DBG3(fes) << "CEF jb: " << jb_.ToString();
            zfeast_hcsrgv(&uplo, &n_, aDouble.GetPointer(), ia_.GetPointer(), ja_.GetPointer(), bDouble.GetPointer(), ib_.GetPointer(), jb_.GetPointer(),
                          fpm_.GetPointer(), &epsout, &loop, I1, &I2, &m0_,
                          E.GetPointer(), X.GetPointer(), &m_, res.GetPointer(), &info_);
        }
        LOG_DBG(fes) << "CEF epsout -> " << epsout;
        // save over the eigenvectors
        for (int i=0;i<n_*m_;i++) {
            vr_[i] = Complex( X[2*i], X[2*i+1] );
        }
    } break;;
    case REAL_SYMMETRIC : {
        LOG_DBG3(fes) << "real-symmetric EVP";
        I1[0] = minVal; LOG_DBG3(fes) << "CEF Emin: " << I1;
        I2 = maxVal; LOG_DBG3(fes) << "CEF Emax: " << I2;
        LOG_DBG3(fes) << "CEF m0: " << m0_; // specifies the initial guess for subspace dimension to be used
        LOG_DBG3(fes) << "CEF x:" << vr_.ToString();
        LOG_DBG3(fes) << "CEF n: " << n_; // size of the problem
        // setup A-matrix
        const SCRS_Matrix<double>* a_const = dynamic_cast<const SCRS_Matrix<double>*>(a_);
        SCRS_Matrix<double>* a = const_cast<SCRS_Matrix<double>*>(a_const);
        LOG_DBG3(fes) << "CEF a_size: " << a->GetNumEntries();
        LOG_DBG3(fes) << "CEF a: " << StdVector<double>::ToString(a->GetNumEntries(), a->GetDataPointer(), 0);
        LOG_DBG3(fes) << "CEF ia: " << ia_.ToString();
        LOG_DBG3(fes) << "CEF ja: " << ja_.ToString();

        // output parameters
        double epsout;
        Vector<double>& E = dynamic_cast<Vector<double>&>(sol); // eigenvalues
        Vector<double>& res = dynamic_cast<Vector<double>&>(err); // Relative residual
        E.Resize(m0_); // e: m found eigenvalues
        Vector<double> X; X.Resize(n_*m0_);
        //assert((int) X.Capacity() >= n_ * m0_);
        res.Resize(m0_); // the first m get the relative residuals


        if(!generalized_) // solve standard EVP
        {
            dfeast_scsrev(&uplo, &n_, a->GetDataPointer(), ia_.GetPointer(), ja_.GetPointer(),
                        fpm_.GetPointer(), &epsout, &loop, I1, &I2, &m0_,
                        E.GetPointer(), X.GetPointer(), &m_, res.GetPointer(), &info_);
        }
        else // generalized EVP
        {
            const SCRS_Matrix<double>* bc = dynamic_cast<const SCRS_Matrix<double>*>(b_);
            SCRS_Matrix<double>* b = const_cast<SCRS_Matrix<double>*>(bc);
            assert((int) ib_.GetSize() == n_ + 1);
            assert(b->GetNumEntries() < b->GetNnz()); // for symmetric matrices!
            LOG_DBG3(fes) << "CEF b_size: " << b->GetNumEntries();
            LOG_DBG3(fes) << "CEF b: " << StdVector<double>::ToString(b->GetNumEntries(), b->GetDataPointer(), 0);
            LOG_DBG3(fes) << "CEF ib: " << ib_.ToString();
            LOG_DBG3(fes) << "CEF jb: " << jb_.ToString();
            dfeast_scsrgv(&uplo, &n_, a->GetDataPointer(), ia_.GetPointer(), ja_.GetPointer(), b->GetDataPointer(), ib_.GetPointer(), jb_.GetPointer(),
                          fpm_.GetPointer(), &epsout, &loop, I1, &I2, &m0_,
                          E.GetPointer(), X.GetPointer(), &m_, res.GetPointer(), &info_);
        }
        LOG_DBG(fes) << "CEF epsout -> " << epsout;
        // save over the eigenvector
        for (int i=0;i<n_*m_;i++) {
            vr_[i] = Complex( X[i], 0.0 );
        }
    } break;;
    default: {
        EXCEPTION("Did you set the Problem Type? Use Setup(...)")
    }
    }
    // resize result
    sol.Resize(m_);
    err.Resize(m_);
    // post log
    LOG_DBG(fes) << "CEF info -> " << info_;
    LOG_DBG(fes) << "CEF loop -> " << loop;
    LOG_DBG(fes) << "CEF m -> " << m_;
    LOG_DBG(fes) << "CEF e -> " << sol.ToString();
    LOG_DBG(fes) << "CEF res -> " << err.ToString();
    LOG_DBG3(fes) << "CEF x:" << vr_.ToString();
}

void FeastEigenSolver::Setup(const BaseMatrix& stiffMat, unsigned int numFreq, double freqShift, bool sort)
{
  numFreq_ = numFreq;
  freqShift_ = freqShift;
  generalized_ = false;
  bloch_ = false;

  // initialize feastinit fpm according to "Extended Eigensolver Input Parameters" in the MKL manual
  // see: https://software.intel.com/en-us/mkl-developer-reference-c-extended-eigensolver-input-parameters
  //int feastparam[64];
  feastinit(fpm_.GetPointer());
  // now get values from the <feast> section of the XML
  bool logging = false; xml_->GetValue("logging", logging, ParamNode::INSERT); if (logging == true) fpm_[0] = 1; // runtime messages to the screen
  uint Ne = 8; xml_->GetValue("Ne", Ne, ParamNode::INSERT); fpm_[1] = Ne; // Number of contour points
  uint stopCrit = 12; xml_->GetValue("stopCrit", stopCrit, ParamNode::INSERT); fpm_[2] = stopCrit; // error trace double precision stopping criteria eps=10^-fpm[2]
  uint maxRefinementLoops = 20; xml_->GetValue("maxRefinementLoops", maxRefinementLoops, ParamNode::INSERT); fpm_[3] = maxRefinementLoops; // max refinement loops
  fpm_[4] = 0;  // no initial subspace
  fpm_[26] = 1; // Extended Eigensolver routines check input matrices
  fpm_[27] = 1; // check b for positive definiteness

  n_ = stiffMat.GetNumRows(); // size of the problem
  uint m0 = 0; xml_->GetValue("m0", m0, ParamNode::INSERT); if (m0>0) m0_ = m0; else m0_=n_; // initial subspace size
  assert(m0_ <= n_);
  vr_.Resize(m0_ * n_, 0.0);

  // store the matrix

  // a
  // Array containing the nonzero elements of either the full matrix A or the upper or lower
  // triangular part of the matrix A, as specified by uplo.
  a_ = &dynamic_cast<const StdMatrix&>(stiffMat);
  assert(a_ != NULL);
  // copy the indices and transform to 1-based

  // ia
  // Array of length n + 1, containing indices of elements in the array a ,
  // such that ia[i - 1] is the index in the array a of the first non-zero element from the row i .
  // The value of the last element ia[n] is equal to the number of non-zeros plus one
  a_->ExportCRSRows(ia_, 1, true); // tailing nnz = number of non-zeros
  ia_[ia_.GetSize()-1] = ia_[ia_.GetSize()-1] + 1; // increase by one, see above

  // ja
  // Array containing the column indices for each non-zero element of the matrix A being represented in the array a.
  // Its length is equal to the length of the array a.
  a_->ExportCRSColumns(ja_, 1); // 1 for 1-based
}

void FeastEigenSolver::Setup( const BaseMatrix& stiffMat, const BaseMatrix& massMat,
            unsigned int numFreq, double freqShift, bool sort, bool bloch)
{
  Setup(stiffMat, numFreq, freqShift, sort);

  generalized_ = true;
  bloch_ = bloch;

  b_ = &dynamic_cast<const StdMatrix&>(massMat);
  assert(b_ != NULL);
  // copy the indices and transform to 1-based
  b_->ExportCRSColumns(jb_, 1); // 1 for 1-based
  b_->ExportCRSRows(ib_, 1, true); // nzz = tailing number of non-zeros
  ib_[ib_.GetSize()-1] = ib_[ib_.GetSize()-1] + 1; // according to doc this must be nzz+1

  assert(a_->GetNumRows() == b_->GetNumRows());
  assert((int) b_->GetNumRows() == n_);
  assert(a_->GetStorageType() == b_->GetStorageType());
  assert(a_->GetStructureType() == b_->GetStructureType());
  assert(a_->GetEntryType() == b_->GetEntryType());
}

void FeastEigenSolver::Setup( const BaseMatrix& stiffMat,
            const BaseMatrix& massMat,
            const BaseMatrix& dampMat,
            unsigned int numFreq, double freqShift , bool sort)
{
  assert(false);
}
int info;



void FeastEigenSolver::GetEigenMode(unsigned int modeNr, Vector<Complex>& mode, bool right){
  mode.Resize(n_, 0.0);
  assert((int) modeNr <= m_); // assume 0-based!!
  if (eigenProblemType_==REAL_SYMMETRIC || eigenProblemType_==COMPLEX_SYMMETRIC || eigenProblemType_==COMPLEX_HERMITIAN ){
    right = true; // for these EVPs we have vl_=vr_
  }
  // select what to return
  StdVector<Complex>* v;
  if (right) {
    v = &vr_;
  }
  else {
    v = &vl_;
  }
  // set return vector
  assert(v->GetSize() >= modeNr*n_);
  for(int i = 0; i < n_; i++) {
    mode[i] = (*v)[modeNr * n_ + i];
  }
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
