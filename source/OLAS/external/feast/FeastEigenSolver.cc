#include <string.h>
#include <def_use_feast.hh>
#include <algorithm>
#include <numeric>

// MKL FEAST is not supported any more (used functions are not available for MKL 2021)
// Interestingly we need no name mangling on Linux/macOS?! Usually we include def_cfs_fortran_interface.hh
// For Linux/macOS we could also include feast.h and feast_sparse.h via extern "C" but this fails for Windows.
// Our standard mangling does not work with FEAST, therefore we use our own feast_cpp.
#include "feast_cpp.h" // we maintain this file in cfsdeps/feast

#include "MatVec/StdMatrix.hh"
#include "MatVec/SCRS_Matrix.hh"
#include "MatVec/CRS_Matrix.hh"

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/Logging/LogConfigurator.hh"

#include "OLAS/external/feast/FeastEigenSolver.hh"

using std::string;

DEFINE_LOG(fes, "feast")

namespace CoupledField {


FeastEigenSolver::FeastEigenSolver(shared_ptr<SolStrategy> strat,
                                   PtrParamNode xml,
                                   PtrParamNode solverList,
                                   PtrParamNode precondList,
                                   PtrParamNode eigenInfo)
  : BaseEigenSolver(strat, xml, solverList, precondList, eigenInfo)
{
  this->p_ = 0;
  this->bloch_ = false;
  this->a_ = NULL;
  this->b_ = NULL;
  this->c_ = NULL;
  this->numFreq_ = 0;
  this->freqShift_ = 0;
  this->fpm_.Resize(128, 0);
  this->n_ = 0;
  this->m0_ = 0;
  this->m_ = 0;
  this->info_ = -10; // not im MKL code (page 1635)
  xml_ = xml;
  eigenInfo_ = eigenInfo;
  eigenSolverName_ = EigenSolverType::FEAST;

  BaseEigenSolver::PostInit();
}

FeastEigenSolver::~FeastEigenSolver()
{

}


void FeastEigenSolver::Setup(const BaseMatrix & A, bool isHermitian){
    // some general settings
    // further case specific paramaters are set in the CalcEigenvalues methods
    // get values from the <feast> section of the XML
    bool logging = false; xml_->GetValue("logging", logging, ParamNode::INSERT); // runtime messages to the screen
    UInt stopCrit = 12; xml_->GetValue("stopCrit", stopCrit, ParamNode::INSERT); // error trace double precision stopping criteria eps=10^-fpm[2]
    UInt maxRefinementLoops = 20; xml_->GetValue("maxRefinementLoops", maxRefinementLoops, ParamNode::INSERT); // max refinement loops
    xml_->GetValue("stochasticEstimate", stochasticEstimate_, ParamNode::INSERT); // use stochastic estimate of eigenvalues within search contour
    // initialize feastinit fpm according to "Extended Eigensolver Input Parameters" in the MKL manual
    // see: https://software.intel.com/en-us/mkl-developer-reference-c-extended-eigensolver-input-parameters
    feastinit(fpm_.GetPointer());
    // set some standard feast parameters
    if (logging == true) fpm_[0] = 1; // runtime messages to the screen
    fpm_[2] = stopCrit; // error trace double precision stopping criteria eps=10^-fpm[2]
    fpm_[3] = maxRefinementLoops; // max refinement loops
    fpm_[4] = 0;  // no initial subspace
    if (stochasticEstimate_) fpm_[13] = 2; // use stochastic estimate for number of eigenvalues
    fpm_[26] = 1; // Extended Eigensolver routines check input matrices
    fpm_[27] = 1; // check b for positive definiteness

    if(xml_->Has("precision")) 
      precision_ = xml_->Get("precision")->As<string>() == "double" ? 0 : 1; // 0 = double, 1 = single (legacy default)
    fpm_[41] = precision_; // FEAST fpm(42): 0 = double, 1 = single precision inner solver
    // possible settings for future implementations
    // elliptical search contour => fpm18, fpm19
    // integration type => fpm16

    // check for symmetry
    this->SetProblemType(A,isHermitian);

    // size of the problem
    n_ = A.GetNumRows();
    UInt m0 = 0; xml_->GetValue("m0", m0, ParamNode::INSERT); if (m0>0) m0_ = m0; else m0_=n_; // initial subspace size

    if (p_ == 2)
    {
      assert(m0_ <= 2*n_); // QEVPS yield up to 2n eigenvalues, but FEAST might not support that...
    }
    else
    {
      assert(m0_ <= n_);
    }
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
    if (p_ != 2) // change order only if setup is not called from quadratic case
    {
      p_ = 1;
    }
    // first set A as for standard EVP
    Setup(A, isHermitian);
    // check the B matrix
    //bool isReal, isStoredSymmetric;
    //this->CheckMatrix(isReal, isStoredSymmetric, B);
    this->SetProblemType(B,isHermitian);
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

void FeastEigenSolver::Setup(const BaseMatrix & K, const BaseMatrix & C, const BaseMatrix & M){
    // currently no hermitian quadratic EVP can appear in CFS
    bool isHermitian = false;
    p_ = 2;
    // set K,M as for generalized EVP -> vectors a_,ia_,ja_,b_,ib_,jb_
    Setup(K, M, isHermitian);
    
    // check the C matrix
    this->SetProblemType(C,isHermitian); //check if problem type is same for C and K
    // setup C-matrix
    // store matrix
    c_ = &dynamic_cast<const StdMatrix&>(C);
    assert(c_ != NULL);
    // copy the indices and transform to 1-based
    c_->ExportCRSColumns(jc_, 1); // 1 for 1-based
    c_->ExportCRSRows(ic_, 1, true); // nzz = tailing number of non-zeros
    ic_[ic_.GetSize()-1] = ic_[ic_.GetSize()-1] + 1; // according to doc this must be nzz+1
    // check if size is equal
    assert(a_->GetNumRows() == c_->GetNumRows());
    assert((int) c_->GetNumRows() == n_);
    assert(a_->GetStorageType() == c_->GetStorageType());
    assert(a_->GetStructureType() == c_->GetStructureType());
    assert(a_->GetEntryType() == c_->GetEntryType());

    // Determine maximum of nonzero elements (according to feast)
    nnza_ = ja_.GetSize();
    nnzc_ = jc_.GetSize();
    nnzb_ = jb_.GetSize();
    nnzmax_ = std::max(nnza_,nnzc_);
    nnzmax_ = std::max(nnzmax_,nnzb_);
    
    /*Stack row and column vectors of CSR format
    * the stacked row array has size 3*(n_ + 1) (system size)
    * the stacked column array has size 3*nnzmax_ (maximum number of non-zero elements)
    * indexing by increasing degree: block 1 = stiffness, block 2 = damping, block 3 = mass
    * refer to FEAST Eigenvalue Solver v4.0 User Guide, section 3.4
    */ 
    isa_.Resize(3*(n_+1));
    for (int i=0; i < n_+1; i++)
    {
      isa_[i] = ia_[i];
      isa_[i+n_+1] = ic_[i]; 
      isa_[i+2*(n_+1)] = ib_[i];
    }
    jsa_.Resize(3*nnzmax_);
    for (int i=0;i<nnza_;i++) 
    {
      jsa_[i] = ja_[i];
    }
    for (int i=0;i<nnzc_;i++) 
    {
      jsa_[i+nnzmax_] = jc_[i];
    }
    for (int i=0;i<nnzb_;i++) 
    {
      jsa_[i+2*nnzmax_] = jb_[i];
    }
}

// interface for routine with customized search contour
void FeastEigenSolver::CalcEigenValues(BaseVector& sol, BaseVector& err){
    // step 1: calculate the integration nodes and weights
    // case 1: rectangle
    if ( xml_->Has("rectangularContour") ) {
      EXCEPTION( "Currently not implemented!" );
    }
    // case 2: full custom
    else if ( xml_->Has("eigenValues/customContour") ) {
      // a customized search contour consists of connected pieces (defined via their end points)
      // each piece is divided into integration intervals
      // currently it is setup in a way that the custom contour parameters overrule the normal contour parameters!

      LOG_DBG3(fes) << "custom contour";

      ParamNodeList contour = xml_->Get("eigenValues/customContour")->GetChildren();
      numP_ = contour.size(); // number of segments
      // resize
      z_edge_.Resize(numP_*2); // 2*P since due to complex values
      t_edge_.Resize(numP_);
      n_edge_.Resize(numP_);
      for (int i=0; i<numP_; i++) {
        // endpoints of segment i zedge
        contour[i]->GetValue("endpoint/real", z_edge_[i*2], ParamNode::EX);
        contour[i]->GetValue("endpoint/imag", z_edge_[i*2+1], ParamNode::EX);
        // type of segment i tedge
        // ellipsoid
        if ( contour[i]->Has("ellipsoid") ) {
          double ratio = 0;
          contour[i]->GetValue("ellipsoid/ratio", ratio, ParamNode::EX);
          if ( ratio>1 ) {
            EXCEPTION( "The ratio of elliptical radii a/b has to be less or equal to 1." )
          }
          else{
            // feast requires the ratio to be the integer of a/b*100
            t_edge_[i] = static_cast<int>(ratio*100);
          }
        }
        // line
        else if ( contour[i]->Has("line") ) {
          t_edge_[i] = 0;
        }
        else {
          EXCEPTION( "Unknown segment type!" );
        }
        // amount of nodes of segment i nedge
        contour[i]->GetValue("nodes", n_edge_[i], ParamNode::EX);
      }
      numNc_ = std::accumulate(n_edge_.begin(), n_edge_.end(), 0); // overall number of nodes

      // log all contour parameters
      LOG_DBG3(fes) << "CEF numP_:" << numP_;
      LOG_DBG3(fes) << "CEF z_edge:" << z_edge_.ToString();
      LOG_DBG3(fes) << "CEF t_edge:" << t_edge_.ToString();
      LOG_DBG3(fes) << "CEF n_edge:" << n_edge_.ToString();
      LOG_DBG3(fes) << "CEF numNc_:" << numNc_;

      if ( ( eigenProblemType_==REAL_SYMMETRIC || eigenProblemType_==COMPLEX_SYMMETRIC || eigenProblemType_==COMPLEX_HERMITIAN )
          && ( p_<2 ) ) { 
        fpm_[1] = numNc_; // integration nodes 
        fpm_[15] = 0; // Gauss integration scheme
      }
      else { 
        fpm_[7] = numNc_; // integration nodes 
        fpm_[15] = 1; // Trapezoidal integration scheme
      }

      // set standard contour parameters to zero since custom search contour is used
      I1_[0] = 0; // mid point or min value
      I2_ = 0; // radius or max value

      // calculate integration nodes (integrationNodes_) and weights (integrationWeights_)
      integrationNodes_.Resize(numNc_*2);
      integrationWeights_.Resize(numNc_*2);
      zfeast_customcontour(&numNc_, &numP_, n_edge_.GetPointer(), t_edge_.GetPointer(), z_edge_.GetPointer(), integrationNodes_.GetPointer(), integrationWeights_.GetPointer());
      LOG_DBG3(fes) << "CEF Zne:" << integrationNodes_.ToString(); // integration nodes
      LOG_DBG3(fes) << "CEF Wne:" << integrationWeights_.ToString(); // integration weights
    }
    else {
      EXCEPTION( "Should only be called if customContour is defined in the solver definition of the xml scheme." )
    }

    // step 2: solve the problem
    DoCalculation(sol, err);

    // step 3: save to .info.xml
    ToInfo();
}

// interface for routine with regular search contour
void FeastEigenSolver::CalcEigenValues(BaseVector& sol, BaseVector& err, Double minVal, Double maxVal){
    LOG_DBG3(fes) << "standard contour";
    // step 1: calculate the integration nodes and weights
    // case 1: hermitian
    if ( ( eigenProblemType_==REAL_SYMMETRIC || eigenProblemType_==COMPLEX_HERMITIAN )
        && ( p_<2 ) ) {
      LOG_DBG3(fes) << "hermitian routine";
      // set parameters
      // in the non-polynomial (p<2) hermitian case the min and max real value of the search contour is required
      I1_[0] = minVal; LOG_DBG3(fes) << "CEF Emin: " << I1_[0]; // min value
      I2_ = maxVal; LOG_DBG3(fes) << "CEF Emax: " << I2_; // max value

      UInt Ne = 8; xml_->GetValue("Ne", Ne, ParamNode::INSERT); // Number of contour points
      fpm_[1] = Ne; // integration nodes 
      fpm_[15] = 0; // Gauss scheme
      fpm_[17] = 100; // circle shape

      // log
      LOG_DBG3(fes) << "CEF contour points: " << fpm_[1]; // contour points
      LOG_DBG3(fes) << "CEF integration type: " << fpm_[15]; // contour points
      LOG_DBG3(fes) << "CEF ellipse ratio a/b*100: " << fpm_[17]; // contour points
      LOG_DBG3(fes) << "CEF m0: " << m0_; // specifies the initial guess for subspace dimension to be used
      LOG_DBG3(fes) << "CEF n: " << n_; // size of the problem

      // calculate integration nodes (integrationNodes_) and weights (integrationWeights_)
      integrationNodes_.Resize(fpm_[1]*2);
      integrationWeights_.Resize(fpm_[1]*2);
      zfeast_contour(&I1_[0], &I2_, &fpm_[1], &fpm_[15], &fpm_[17], integrationNodes_.GetPointer(), integrationWeights_.GetPointer());
      LOG_DBG3(fes) << "CEF Zne:" << integrationNodes_.ToString(); // integration nodes
      LOG_DBG3(fes) << "CEF Wne:" << integrationWeights_.ToString(); // integration weights
    }
    // case 2: non-hermitian/polynomial
    else {
      LOG_DBG3(fes) << "non-hermitian/poly. routine";
      // set parameters
      // in this case the (complex) mid point and radius of the circular search contour is required
      I1_[0] = 0.5*(maxVal+minVal); I1_[1]=0.0; LOG_DBG3(fes) << "CEF Emid: " << I1_[0]; // mid point
      I2_ = 0.5*(maxVal-minVal); LOG_DBG3(fes) << "CEF r: " << I2_; // radius
      LOG_DBG3(fes) << "CEF m0: " << m0_; // specifies the initial guess for subspace dimension to be used
      LOG_DBG3(fes) << "CEF n: " << n_; // size of the problem

      UInt Ne = 16; xml_->GetValue("Ne", Ne, ParamNode::INSERT); // Number of contour points
      fpm_[7] = Ne; // integration nodes 
      fpm_[15] = 1; // Trapezoidal scheme
      fpm_[17] = 100; // circle shape
      fpm_[18] = 0; // ellipsis angle

      // calculate integration nodes (integrationNodes_) and weights (integrationWeights_)
      integrationNodes_.Resize(fpm_[7]*2);
      integrationWeights_.Resize(fpm_[7]*2);
      zfeast_gcontour(I1_, &I2_, &fpm_[7], &fpm_[15], &fpm_[17], &fpm_[18], integrationNodes_.GetPointer(), integrationWeights_.GetPointer());
      LOG_DBG3(fes) << "CEF Zne:" << integrationNodes_.ToString(); // integration nodes
      LOG_DBG3(fes) << "CEF Wne:" << integrationWeights_.ToString(); // integration weights
    }

    // step 2: solve the problem
    DoCalculation(sol, err);

    // step 3: save to .info.xml
    ToInfo();
}

// actually solve the problem for a random seed
void FeastEigenSolver::CalculationAttempt(BaseVector& sol, BaseVector& err) {
    assert(a_ != NULL);
    assert(a_->GetStructureType() == BaseMatrix::SPARSE_MATRIX);
    if (a_->GetStructureType() != BaseMatrix::SPARSE_MATRIX) {
        EXCEPTION("FEAST is currently only implemented for sparse matrices.")
    }
    assert((int) ia_.GetSize() == n_ + 1);
 
    switch (eigenProblemType_) {

    case REAL_GENERAL : {
        LOG_DBG3(fes) << "real-general EVP";
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
        
        switch (p_) {
          case 2: { // quadratic EVP
            const CRS_Matrix<double>* c_const = dynamic_cast<const CRS_Matrix<double>*>(c_);
            CRS_Matrix<double>* c = const_cast<CRS_Matrix<double>*>(c_const);
            const CRS_Matrix<double>* b_const = dynamic_cast<const CRS_Matrix<double>*>(b_);
            CRS_Matrix<double>* b = const_cast<CRS_Matrix<double>*>(b_const);
            LOG_DBG3(fes) << "CEF ic: " << ic_.ToString();
            LOG_DBG3(fes) << "CEF jc: " << jc_.ToString();
            LOG_DBG3(fes) << "CEF ib: " << ib_.ToString();
            LOG_DBG3(fes) << "CEF jb: " << jb_.ToString();
            LOG_DBG3(fes) << "Size of stacked vectors: 3x" << nnzmax_;
            
            // store all non-zero values of the system matrices in a stacked vector
            Vector<double> sa;
            sa.Resize(3*nnzmax_);
            for (int i=0;i<nnza_;i++) 
            {
              sa[i] = a->GetDataPointer()[i];
            }
            for (int i=0;i<nnzc_;i++) 
            {
              sa[i+nnzmax_] = c->GetDataPointer()[i];
            }
            for (int i=0;i<nnzb_;i++) 
            {
              sa[i+2*nnzmax_] = b->GetDataPointer()[i];
            }
            dfeast_gcsrpevx(&p_,&n_,sa.GetPointer(),isa_.GetPointer(),jsa_.GetPointer(),fpm_.GetPointer(), &epsout, &loop_, I1_, &I2_, &m0_,
                          E.GetPointer(), X.GetPointer(), &m_, res.GetPointer(), &info_, integrationNodes_.GetPointer(), integrationWeights_.GetPointer());
            break;
          }
          case 1: { // generalized EVP
            const CRS_Matrix<double>* bc = dynamic_cast<const CRS_Matrix<double>*>(b_);
            CRS_Matrix<double>* b = const_cast<CRS_Matrix<double>*>(bc);
            assert((int) ib_.GetSize() == n_ + 1);
            //assert(b->GetNumEntries() < b->GetNnz()); // for symmetric matrices!
            //LOG_DBG3(fes) << "CEF b_size: " << b->GetNumEntries();
            //LOG_DBG3(fes) << "CEF b: " << StdVector<double>::ToString(b->GetNumEntries(), b->GetDataPointer(), 0);
            LOG_DBG3(fes) << "CEF ib: " << ib_.ToString();
            LOG_DBG3(fes) << "CEF jb: " << jb_.ToString();
            dfeast_gcsrgvx(&n_, a->GetDataPointer(), ia_.GetPointer(), ja_.GetPointer(), b->GetDataPointer(), ib_.GetPointer(), jb_.GetPointer(),
                          fpm_.GetPointer(), &epsout, &loop_, I1_, &I2_, &m0_,
                          E.GetPointer(), X.GetPointer(), &m_, res.GetPointer(), &info_, integrationNodes_.GetPointer(), integrationWeights_.GetPointer());
            break;
          }
          case 0: {// solve standard EVP
            dfeast_gcsrevx(&n_, a->GetDataPointer(), ia_.GetPointer(), ja_.GetPointer(),
                        fpm_.GetPointer(), &epsout, &loop_, I1_, &I2_, &m0_,
                        E.GetPointer(), X.GetPointer(), &m_, res.GetPointer(), &info_, integrationNodes_.GetPointer(), integrationWeights_.GetPointer());
            break;
          }
        }
        if (stochasticEstimate_)
        {
          std::cout << "Stochastic estimate for the number of eigenvalues: " << m_ << "\n";
          // return 0 eigenvalue and residual to avoid exception when eigenvalues/ error bounds are expected
          E_c.Resize(1);
          res_c.Resize(1);
          vr_.Resize(n_);
          vl_.Resize(n_);
        }
        else{
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
                vl_[i*n_+j] = Complex( X[2*n_*m0_+i*n_*2+2*j], X[2*n_*m0_+i*n_*2+2*j+1]);
            }
          }
          LOG_DBG(fes) << "CEF epsout -> " << epsout;
        }
        
        
    } break;;
    case COMPLEX_SYMMETRIC: {
        LOG_DBG3(fes) << "complex-symmetric EVP";
        // setup A-matrix
        const SCRS_Matrix<Complex>* a_const = dynamic_cast<const SCRS_Matrix<Complex>*>(a_);
        SCRS_Matrix<Complex>* a = const_cast<SCRS_Matrix<Complex>*>(a_const);
        Vector<double> aDouble; 
        if (p_ != 2)
        {
          aDouble.Resize(2*a->GetNnz(),0.0); // factor 2 for complex
          for (int i=0;i<(int)a->GetNnz();i++) {
            aDouble[2*i] = a->GetDataPointer()[i].real();
            aDouble[2*i+1] = a->GetDataPointer()[i].imag();
          }
        }
        
        LOG_DBG3(fes) << "CEF a_size: " << a->GetNumEntries();
        LOG_DBG3(fes) << "CEF a: " << StdVector<Complex>::ToString(a->GetNumEntries(), a->GetDataPointer());
        LOG_DBG3(fes) << "CEF ia: " << ia_.ToString();
        LOG_DBG3(fes) << "CEF ja: " << ja_.ToString();
        // output parameters
        double epsout;
        Vector<Complex>& E_c = dynamic_cast<Vector<Complex>&>(sol); // eigenvalues
        Vector<double> E; E.Resize(2*m0_,0.0); // factor 2 for complex
        Vector<Complex>& res_c = dynamic_cast<Vector<Complex>&>(err); // Relative residual
        Vector<double> res; res.Resize(2*m0_,0.0); // factor 2 for complex
        Vector<double> X; X.Resize(2*n_*m0_);// eigenvectors (if needed) //factor 2 for complex

        switch (p_) {
          case 2: { // quadratic EVP
            const SCRS_Matrix<Complex>* c_const = dynamic_cast<const SCRS_Matrix<Complex>*>(c_);
            SCRS_Matrix<Complex>* c = const_cast<SCRS_Matrix<Complex>*>(c_const);
            const SCRS_Matrix<Complex>* b_const = dynamic_cast<const SCRS_Matrix<Complex>*>(b_);
            SCRS_Matrix<Complex>* b = const_cast<SCRS_Matrix<Complex>*>(b_const);
            LOG_DBG3(fes) << "CEF mc: " << StdVector<Complex>::ToString(nnzc_, c->GetDataPointer());
            LOG_DBG3(fes) << "CEF ic: " << ic_.ToString();
            LOG_DBG3(fes) << "CEF jc: " << jc_.ToString();
            LOG_DBG3(fes) << "CEF mm: " << StdVector<Complex>::ToString(nnzb_, b->GetDataPointer());
            LOG_DBG3(fes) << "CEF im: " << ib_.ToString();
            LOG_DBG3(fes) << "CEF jm: " << jb_.ToString();
            LOG_DBG3(fes) << "Size of stacked vectors: 3x" << nnzmax_;
            
            // store all non-zero values of the system matrices in a stacked vector
            Vector<double> sa;
            sa.Resize(6*nnzmax_); // 2*3*nnzmax_ for complex
            for (int i=0;i<nnza_;i++) 
            {
              sa[2*i] = a->GetDataPointer()[i].real();
              sa[2*i+1] = a->GetDataPointer()[i].imag();
            }
            for (int i=0;i<nnzc_;i++) 
            {
              sa[2*i+2*nnzmax_] = c->GetDataPointer()[i].real();
              sa[2*i+1+2*nnzmax_] = c->GetDataPointer()[i].imag();
            }
            for (int i=0;i<nnzb_;i++) 
            {
              sa[2*i+4*nnzmax_] = b->GetDataPointer()[i].real();
              sa[2*i+1+4*nnzmax_] = b->GetDataPointer()[i].imag();
            }
            zfeast_scsrpevx(&uplo_,&p_,&n_,sa.GetPointer(),isa_.GetPointer(),jsa_.GetPointer(),fpm_.GetPointer(), &epsout, &loop_, I1_, &I2_, &m0_,
                          E.GetPointer(), X.GetPointer(), &m_, res.GetPointer(), &info_, integrationNodes_.GetPointer(), integrationWeights_.GetPointer());
            break;
          }  
          case 1: { // generalized EVP
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
            LOG_DBG3(fes) << "CEF b: " << StdVector<Complex>::ToString(b->GetNumEntries(), b->GetDataPointer());
            LOG_DBG3(fes) << "CEF ib: " << ib_.ToString();
            LOG_DBG3(fes) << "CEF jb: " << jb_.ToString();
            zfeast_scsrgvx(&uplo_, &n_, aDouble.GetPointer(), ia_.GetPointer(), ja_.GetPointer(), bDouble.GetPointer(), ib_.GetPointer(), jb_.GetPointer(),
                fpm_.GetPointer(), &epsout, &loop_, I1_, &I2_, &m0_,
                E.GetPointer(), X.GetPointer(), &m_, res.GetPointer(), &info_, integrationNodes_.GetPointer(), integrationWeights_.GetPointer());
            break;
          }
          case 0: { // solve standard EVP
            zfeast_scsrevx(&uplo_, &n_, aDouble.GetPointer(), ia_.GetPointer(), ja_.GetPointer(),
                fpm_.GetPointer(), &epsout, &loop_, I1_, &I2_, &m0_,
                E.GetPointer(), X.GetPointer(), &m_, res.GetPointer(), &info_, integrationNodes_.GetPointer(), integrationWeights_.GetPointer());
            break;
          }
        }
        if (stochasticEstimate_)
        {
          std::cout << "Stochastic estimate for the number of eigenvalues: " << m_ << "\n";
          // return 0 eigenvalue and residual to avoid exception when eigenvalues/ error bounds are expected
          E_c.Resize(1);
          res_c.Resize(1);
          vr_.Resize(n_);
          vl_.Resize(n_);
        }
        else
        {
          LOG_DBG(fes) << "CEF epsout -> " << epsout;
          // put the solution into the complex vectors
          E_c.Resize(m_);
          res_c.Resize(m_);
          vr_.Resize(m_*n_);
          LOG_DBG3(fes) << "CEF # of found eigen values =  " << m_;
          LOG_DBG3(fes) << "CEF # elements in X: "<<X.GetSize();
          for (int i=0;i<m_;i++){
            E_c[i] = Complex( E[2*i],E[2*i+1] );
            res_c[i] = Complex( res[2*i], res[2*i+1] );
            for (int j=0;j<n_;j++){
              vr_[i*n_+j] = Complex( X[i*n_*2+2*j], X[i*n_*2+2*j+1]);
            }
          }
        }
        
    } break;;
    case COMPLEX_GENERAL: {
        LOG_DBG3(fes) << "complex-general EVP";
        // setup A-matrix
        const CRS_Matrix<Complex>* a_const = dynamic_cast<const CRS_Matrix<Complex>*>(a_);
        CRS_Matrix<Complex>* a = const_cast<CRS_Matrix<Complex>*>(a_const);
        Vector<double> aDouble;
        if (p_ != 2)
        {
          aDouble.Resize(2*a->GetNnz(),0.0); // factor 2 for complex
          for (int i=0;i<(int)a->GetNnz();i++) {
            aDouble[2*i] = a->GetDataPointer()[i].real();
            aDouble[2*i+1] = a->GetDataPointer()[i].imag();
          }
        }
        //LOG_DBG3(fes) << "CEF a_size: " << a->GetNumEntries();
        //LOG_DBG3(fes) << "CEF a: " << StdVector<double>::ToString(a->GetNumEntries(), a->GetDataPointer(), 0);
        LOG_DBG3(fes) << "CEF ia: " << ia_.ToString();
        LOG_DBG3(fes) << "CEF ja: " << ja_.ToString();

        // output parameters
        double epsout;
        assert(sol.IsComplex()); // otherwise bad cast below ...
        Vector<Complex>& E_c = dynamic_cast<Vector<Complex>&>(sol); // eigenvalues
        Vector<double> E; E.Resize(2*m0_,0.0); // factor 2 for complex
        assert(err.IsComplex()); // otherwise bad cast below ...
        Vector<Complex>& res_c = dynamic_cast<Vector<Complex>&>(err); // Relative residual
        Vector<double> res; res.Resize(2*m0_,0.0); // factor 2 for complex
        Vector<double> X; X.Resize(2*n_*m0_*2);// eigenvectors (if needed) //factor 2 for complex // factor 2 for L and R
        //assert((int) vr_.Capacity() >= 4 * n_ * m0_);

        switch (p_) {
          case 2: { // quadratic EVP
            const CRS_Matrix<Complex>* c_const = dynamic_cast<const CRS_Matrix<Complex>*>(c_);
            CRS_Matrix<Complex>* c = const_cast<CRS_Matrix<Complex>*>(c_const);
            const CRS_Matrix<Complex>* b_const = dynamic_cast<const CRS_Matrix<Complex>*>(b_);
            CRS_Matrix<Complex>* b = const_cast<CRS_Matrix<Complex>*>(b_const);
            LOG_DBG3(fes) << "CEF ic: " << ic_.ToString();
            LOG_DBG3(fes) << "CEF jc: " << jc_.ToString();
            LOG_DBG3(fes) << "CEF im: " << ib_.ToString();
            LOG_DBG3(fes) << "CEF jm: " << jb_.ToString();
            LOG_DBG3(fes) << "Size of stacked vectors: 3x" << nnzmax_;
            
            // store all non-zero values of the system matrices in a stacked vector
            Vector<double> sa;
            sa.Resize(6*nnzmax_); // 2*3*nnzmax_ for complex
            for (int i=0;i<nnza_;i++) 
            {
              sa[2*i] = a->GetDataPointer()[i].real();
              sa[2*i+1] = a->GetDataPointer()[i].imag();
            }
            for (int i=0;i<nnzc_;i++) 
            {
              sa[2*i+2*nnzmax_] = c->GetDataPointer()[i].real();
              sa[2*i+1+2*nnzmax_] = c->GetDataPointer()[i].imag();
            }
            for (int i=0;i<nnzb_;i++) 
            {
              sa[2*i+4*nnzmax_] = b->GetDataPointer()[i].real();
              sa[2*i+1+4*nnzmax_] = b->GetDataPointer()[i].imag();
            }
            zfeast_gcsrpevx(&p_,&n_,sa.GetPointer(),isa_.GetPointer(),jsa_.GetPointer(),fpm_.GetPointer(), &epsout, &loop_, I1_, &I2_, &m0_,
                          E.GetPointer(), X.GetPointer(), &m_, res.GetPointer(), &info_, integrationNodes_.GetPointer(), integrationWeights_.GetPointer());
            break;
          }  
          case 1: { // generalized EVP
            const CRS_Matrix<Complex>* bc = dynamic_cast<const CRS_Matrix<Complex>*>(b_);
            CRS_Matrix<Complex>* b = const_cast<CRS_Matrix<Complex>*>(bc);
            Vector<double> bDouble; bDouble.Resize(2*b->GetNnz(),0.0); // factor 2 for complex
            for (int i=0;i<(int)b->GetNnz();i++) {
              bDouble[2*i] = b->GetDataPointer()[i].real();
              bDouble[2*i+1] = b->GetDataPointer()[i].imag();
            }
            assert((int) ib_.GetSize() == n_ + 1);
            //LOG_DBG3(fes) << "CEF b_size: " << b->GetNumEntries();
            //LOG_DBG3(fes) << "CEF b: " << StdVector<double>::ToString(b->GetNumEntries(), b->GetDataPointer(), 0);
            LOG_DBG3(fes) << "CEF ib: " << ib_.ToString();
            LOG_DBG3(fes) << "CEF jb: " << jb_.ToString();
            zfeast_gcsrgvx(&n_, aDouble.GetPointer(), ia_.GetPointer(), ja_.GetPointer(), bDouble.GetPointer(), ib_.GetPointer(), jb_.GetPointer(),
                fpm_.GetPointer(), &epsout, &loop_, I1_, &I2_, &m0_,
                E.GetPointer(), X.GetPointer(), &m_, res.GetPointer(), &info_, integrationNodes_.GetPointer(), integrationWeights_.GetPointer());
            break;
          }
          case 0: { // solve standard EVP
            zfeast_gcsrevx(&n_, aDouble.GetPointer(), ia_.GetPointer(), ja_.GetPointer(),
                fpm_.GetPointer(), &epsout, &loop_, I1_, &I2_, &m0_,
                E.GetPointer(), X.GetPointer(), &m_, res.GetPointer(), &info_, integrationNodes_.GetPointer(), integrationWeights_.GetPointer());
            break;
          }
        }
        if (stochasticEstimate_)
        {
          std::cout << "Stochastic estimate for the number of eigenvalues: " << m_ << "\n";
          // return 0 eigenvalue and residual to avoid exception when eigenvalues/ error bounds are expected
          E_c.Resize(1);
          res_c.Resize(1);
          vr_.Resize(n_);
          vl_.Resize(n_);
        }
        else 
        {
          // put the solution into the complex vectors
          // structure of X matrix is explained in feast V3.0 manual chapter 3.1.2 table 3 (page 17/37)
          // X matrix is column-wise (eigenvector-wise) casted into array
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
              vl_[i*n_+j] = Complex( X[2*n_*m0_+i*n_*2+2*j], X[2*n_*m0_+i*n_*2+2*j+1]);
            }
          }
          LOG_DBG(fes) << "CEF epsout -> " << epsout;
        }
        
    } break;;
    case COMPLEX_HERMITIAN: {
        LOG_DBG3(fes) << "complex-hermitian EVP (with symmetric matrix storage)";
        // setup A-matrix
        const SCRS_Matrix<Complex>* a_const = dynamic_cast<const SCRS_Matrix<Complex>*>(a_);
        SCRS_Matrix<Complex>* a = const_cast<SCRS_Matrix<Complex>*>(a_const);
        Vector<double> aDouble; aDouble.Resize(2*a->GetNnz(),0.0); // factor 2 for complex
        for (int i=0;i<(int)a->GetNnz();i++) {
            aDouble[2*i] = a->GetDataPointer()[i].real();
            aDouble[2*i+1] = a->GetDataPointer()[i].imag();
        }
        LOG_DBG3(fes) << "CEF a_size: " << a->GetNumEntries();
        LOG_DBG3(fes) << "CEF a: " << StdVector<Complex>::ToString(a->GetNumEntries(), a->GetDataPointer());
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

        switch (p_) {
          case 2: { // quadratic EVP
            EXCEPTION("Quadratic hermitian EVP is not implemented");
            break;
          }
          case 1: // generalized EVP
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
            LOG_DBG3(fes) << "CEF b: " << StdVector<Complex>::ToString(b->GetNumEntries(), b->GetDataPointer());
            LOG_DBG3(fes) << "CEF ib: " << ib_.ToString();
            LOG_DBG3(fes) << "CEF jb: " << jb_.ToString();
            zfeast_hcsrgvx(&uplo_, &n_, aDouble.GetPointer(), ia_.GetPointer(), ja_.GetPointer(), bDouble.GetPointer(), ib_.GetPointer(), jb_.GetPointer(),
                          fpm_.GetPointer(), &epsout, &loop_, I1_, &I2_, &m0_,
                          E.GetPointer(), X.GetPointer(), &m_, res.GetPointer(), &info_, integrationNodes_.GetPointer(), integrationWeights_.GetPointer());
            break;
          }
          case 0: { // solve standard EVP
            zfeast_hcsrevx(&uplo_, &n_, aDouble.GetPointer(), ia_.GetPointer(), ja_.GetPointer(),
                        fpm_.GetPointer(), &epsout, &loop_, I1_, &I2_, &m0_,
                        E.GetPointer(), X.GetPointer(), &m_, res.GetPointer(), &info_, integrationNodes_.GetPointer(), integrationWeights_.GetPointer());
            break;
          }
        }
        if (stochasticEstimate_)
        {
          std::cout << "Stochastic estimate for the number of eigenvalues: " << m_ << "\n";
          // return 0 eigenvalue and residual to avoid exception when eigenvalues/ error bounds are expected
          E.Resize(1);
          res.Resize(1);
          vr_.Resize(n_);
        }
        else 
        {
          LOG_DBG(fes) << "CEF epsout -> " << epsout;
          // save over the eigenvectors
          for (int i=0;i<n_*m_;i++) {
              vr_[i] = Complex( X[2*i], X[2*i+1] );
          }
        }
        
    } break;;
    case REAL_SYMMETRIC : {
        LOG_DBG3(fes) << "real-symmetric EVP";
        // setup A-matrix
        const SCRS_Matrix<double>* a_const = dynamic_cast<const SCRS_Matrix<double>*>(a_);
        SCRS_Matrix<double>* a = const_cast<SCRS_Matrix<double>*>(a_const);
        LOG_DBG3(fes) << "CEF a_size: " << a->GetNumEntries();
        LOG_DBG3(fes) << "CEF a: " << StdVector<double>::ToString(a->GetNumEntries(), a->GetDataPointer());
        LOG_DBG3(fes) << "CEF ia: " << ia_.ToString();
        LOG_DBG3(fes) << "CEF ja: " << ja_.ToString();

        // declare FEAST output variables that are independent of result type
        double epsout;
        Vector<double> X;
        
        switch (p_) {
          case 2: { // quadratic EVP
            // complex output parameters (eigenvalues for quadratic case can become complex)
            assert(sol.IsComplex());
            Vector<Complex>& E_c = dynamic_cast<Vector<Complex>&>(sol); // eigenvalues
            Vector<double> E; E.Resize(2*m0_,0.0); // factor 2 for complex
            assert(err.IsComplex());
            Vector<Complex>& res_c = dynamic_cast<Vector<Complex>&>(err); // Relative residual
            Vector<double> res;res.Resize(2*m0_,0.0); // factor 2 for complex
            X.Resize(2*n_*m0_*2);// eigenvectors (if needed) //factor 2 for complex // factor 2 for L and R
            const SCRS_Matrix<double>* c_const = dynamic_cast<const SCRS_Matrix<double>*>(c_);
            SCRS_Matrix<double>* c = const_cast<SCRS_Matrix<double>*>(c_const);
            const SCRS_Matrix<double>* b_const = dynamic_cast<const SCRS_Matrix<double>*>(b_);
            SCRS_Matrix<double>* b = const_cast<SCRS_Matrix<double>*>(b_const);
            LOG_DBG3(fes) << "CEF cm: " << StdVector<double>::ToString(nnzc_, c->GetDataPointer());
            LOG_DBG3(fes) << "CEF ic: " << ic_.ToString();
            LOG_DBG3(fes) << "CEF jc: " << jc_.ToString();
            LOG_DBG3(fes) << "CEF mm: " << StdVector<double>::ToString(nnzb_, b->GetDataPointer());
            LOG_DBG3(fes) << "CEF im: " << ib_.ToString();
            LOG_DBG3(fes) << "CEF jm: " << jb_.ToString();
            LOG_DBG3(fes) << "Size of stacked vectors: 3x" << nnzmax_;
            
            // store all non-zero values of the system matrices in a stacked vector
            Vector<double> sa;
            sa.Resize(3*nnzmax_);
            for (int i=0;i<nnza_;i++) 
            {
              sa[i] = a->GetDataPointer()[i];
            }
            for (int i=0;i<nnzc_;i++) 
            {
              sa[i+nnzmax_] = c->GetDataPointer()[i];
            }
            for (int i=0;i<nnzb_;i++) 
            {
              sa[i+2*nnzmax_] = b->GetDataPointer()[i];
            }
            dfeast_scsrpevx(&uplo_, &p_,&n_,sa.GetPointer(),isa_.GetPointer(),jsa_.GetPointer(),fpm_.GetPointer(), &epsout, &loop_, I1_, &I2_, &m0_,
                          E.GetPointer(), X.GetPointer(), &m_, res.GetPointer(), &info_, integrationNodes_.GetPointer(), integrationWeights_.GetPointer());
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
                  vl_[i*n_+j] = Complex( X[2*n_*m0_+i*n_*2+2*j], X[2*n_*m0_+i*n_*2+2*j+1]);
              }
            }
            break;
          }
          case 1: { // generalized EVP
            // output parameters
            Vector<double>& E = dynamic_cast<Vector<double>&>(sol); // eigenvalues
            Vector<double>& res = dynamic_cast<Vector<double>&>(err); // Relative residual
            E.Resize(m0_); // e: m found eigenvalues
            X.Resize(n_*m0_);
            res.Resize(m0_); // the first m get the relative residuals
            const SCRS_Matrix<double>* bc = dynamic_cast<const SCRS_Matrix<double>*>(b_);
            SCRS_Matrix<double>* b = const_cast<SCRS_Matrix<double>*>(bc);
            assert((int) ib_.GetSize() == n_ + 1);
            assert(b->GetNumEntries() < b->GetNnz()); // for symmetric matrices!
            LOG_DBG3(fes) << "CEF b_size: " << b->GetNumEntries();
            LOG_DBG3(fes) << "CEF b: " << StdVector<double>::ToString(b->GetNumEntries(), b->GetDataPointer());
            LOG_DBG3(fes) << "CEF ib: " << ib_.ToString();
            LOG_DBG3(fes) << "CEF jb: " << jb_.ToString();
            dfeast_scsrgvx(&uplo_, &n_, a->GetDataPointer(), ia_.GetPointer(), ja_.GetPointer(), b->GetDataPointer(), ib_.GetPointer(), jb_.GetPointer(),
                          fpm_.GetPointer(), &epsout, &loop_, I1_, &I2_, &m0_,
                          E.GetPointer(), X.GetPointer(), &m_, res.GetPointer(), &info_, integrationNodes_.GetPointer(), integrationWeights_.GetPointer());
            if (stochasticEstimate_)
            { 
              std::cout << "Stochastic estimate for the number of eigenvalues: " << m_ << "\n";
              E.Resize(1);
              res.Resize(1);
              vr_.Resize(n_);
            }
            else 
            {
              // save over the eigenvector
              for (int i=0;i<n_*m_;i++) {
                  vr_[i] = Complex( X[i], 0.0 );
              }
            }
            break;
          }
          case 0: { // solve standard EVP
            // output parameters
            Vector<double>& E = dynamic_cast<Vector<double>&>(sol); // eigenvalues
            Vector<double>& res = dynamic_cast<Vector<double>&>(err); // Relative residual
            E.Resize(m0_); // e: m found eigenvalues
            X.Resize(n_*m0_);
            res.Resize(m0_); // the first m get the relative residuals
            dfeast_scsrevx(&uplo_, &n_, a->GetDataPointer(), ia_.GetPointer(), ja_.GetPointer(),
                        fpm_.GetPointer(), &epsout, &loop_, I1_, &I2_, &m0_,
                        E.GetPointer(), X.GetPointer(), &m_, res.GetPointer(), &info_, integrationNodes_.GetPointer(), integrationWeights_.GetPointer());
            if (stochasticEstimate_)
            {
              std::cout << "Stochastic estimate for the number of eigenvalues: " << m_ << "\n";
              // return 0 eigenvalue and residual to avoid exception when eigenvalues/ error bounds are expected
              E.Resize(1);
              res.Resize(1);
              vr_.Resize(n_);
            }
            else 
            {
              // save over the eigenvector
              for (int i=0;i<n_*m_;i++) {
                  vr_[i] = Complex( X[i], 0.0 );
              }
            }
            break; 
          }
        }
        LOG_DBG(fes) << "CEF epsout -> " << epsout;
    } break;;
    default: {
        EXCEPTION("Did you set the Problem Type? Use Setup(...)")
    }
    }
    if (!stochasticEstimate_)
    {
      // resize result
      sol.Resize(m_);
      err.Resize(m_);
      // post log
      LOG_DBG(fes) << "CEF info -> " << info_;
      LOG_DBG(fes) << "CEF loop -> " << loop_;
      LOG_DBG(fes) << "CEF m -> " << m_;
      LOG_DBG(fes) << "CEF e -> " << sol.ToString();
      LOG_DBG(fes) << "CEF res -> " << err.ToString();
      LOG_DBG3(fes) << "CEF x:" << vr_.ToString();
      // check info
      if (info_!=0) {
        WARN( FeastInfo(info_) );
      }
    }
}

// FEAST can fail intermittently with info=-3 (reduced eigenproblem: B not positive
// definite) from an unlucky random initial subspace. Re-running CalculationAttempt
// draws a fresh subspace (if FEAST re-seeds), so retry a few times before giving up.
void FeastEigenSolver::DoCalculation(BaseVector& sol, BaseVector& err) {
  const int feastMaxAttempts = 3;
  for(int feastAttempt = 1; feastAttempt <= feastMaxAttempts; feastAttempt++) 
  {
    std::cout << ">> FEAST: solve attempt " << feastAttempt << " of " << feastMaxAttempts << std::endl;
    CalculationAttempt(sol, err);
    if(info_ == 0) {
      if(feastAttempt > 1)
        std::cout << ">> FEAST: SUCCEEDED on attempt " << feastAttempt << " of " << feastMaxAttempts << std::endl;
      return;
    }
    if(info_ == -3 && feastAttempt < feastMaxAttempts) {
      std::cout << ">> FEAST: info=-3 (reduced system / B not positive definite, retry with another random subspace" << std::endl;
      continue;
    }
    if(info_ == -3)
      std::cout << ">> FEAST: still info=-3 after " << feastMaxAttempts << " attempts -> giving up " << std::endl;
    return;
  }
}

void FeastEigenSolver::Setup(const BaseMatrix& stiffMat, unsigned int numFreq, double freqShift, bool sort)
{
  numFreq_ = numFreq;
  freqShift_ = freqShift;
  p_ = 0;
  bloch_ = false;

  // initialize feastinit fpm according to "Extended Eigensolver Input Parameters" in the MKL manual
  // see: https://software.intel.com/en-us/mkl-developer-reference-c-extended-eigensolver-input-parameters
  //int feastparam[64];
  feastinit(fpm_.GetPointer());
  // now get values from the <feast> section of the XML
  bool logging = false; xml_->GetValue("logging", logging, ParamNode::INSERT); if (logging == true) fpm_[0] = 1; // runtime messages to the screen
  UInt Ne = 8; xml_->GetValue("Ne", Ne, ParamNode::INSERT); fpm_[1] = Ne; // Number of contour points
  UInt stopCrit = 12; xml_->GetValue("stopCrit", stopCrit, ParamNode::INSERT); fpm_[2] = stopCrit; // error trace double precision stopping criteria eps=10^-fpm[2]
  UInt maxRefinementLoops = 20; xml_->GetValue("maxRefinementLoops", maxRefinementLoops, ParamNode::INSERT); fpm_[3] = maxRefinementLoops; // max refinement loops
  fpm_[4] = 0;  // no initial subspace
  fpm_[26] = 1; // Extended Eigensolver routines check input matrices
  fpm_[27] = 1; // check b for positive definiteness

  n_ = stiffMat.GetNumRows(); // size of the problem
  UInt m0 = 0; xml_->GetValue("m0", m0, ParamNode::INSERT); if (m0>0) m0_ = m0; else m0_=n_; // initial subspace size
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

  p_ = 1;
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


void FeastEigenSolver::GetEigenMode(unsigned int modeNr, Vector<Complex>& mode, bool right){
  mode.Resize(n_, 0.0);
  if (!stochasticEstimate_)
  {
    assert((int) modeNr <= m_); // assume 0-based!!
  }
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

std::string FeastEigenSolver::FeastInfo(Integer info) {
  std::string msg;
  switch (info) {
    case 202:
      msg="Problem with size of the system N"; break;
    case 201:
      msg="Problem with size of subspace M0"; break;
    case 200:
      msg="Problem with Emin,Emax or Emid,r"; break;
    case 6:
      msg="FEAST converges but subspace is not bi-orthonormal"; break;
    case 5:
      msg="Only stochastic estimation of #eigenvalues returned fpm(14)=2"; break;
    case 4:
      msg="Only the subspace has been returned using fpm(14)=1"; break;
    case 3:
      msg="Size of the subspace M0 is too small (M0 less than M)"; break;
    case 2:
      msg="No Convergence (#iteration loops greater than fpm(4))"; break;
    case 1:
      msg="No Eigenvalue found in the search interval"; break;
    case 0:
      msg="Successful exit"; break;
    case -1:
      msg="Internal error for allocation memory"; break;
    case -2:
      msg="Internal error of the inner system solver in FEAST predefined interfaces"; break;
    case -3:
      msg="Internal error of the reduced eigenvalue solver"
          "Possible cause for Hermitian problem: matrix B may not be positive definite";
      break;
    default:
      if (info<-100){
        msg="Problem with the "+std::to_string(info+100)+"th argument of the FEAST interface";
      } else if (info>100) {
        msg="Problem with "+std::to_string(info-100)+"th value of the input FEAST parameter (i.e fpm("+std::to_string(info-100)+"))";
      } else {
        msg="Unknown Error Code.";
      }
  }
  return msg;
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

  // stochastic estimate
  if( stochasticEstimate_ ){
    PtrParamNode setup = eigenInfo_->Get(ParamNode::HEADER);
    setup->Get("N_estimated")->SetValue(m_);
  }

  eigenInfo_->Get("innerSolverPrecision")->SetValue(precision_ == 0 ? "double" : "single");
  // contour information
  // setup tags
  PtrParamNode contour = eigenInfo_->Get("contour",ParamNode::APPEND);
  PtrParamNode nodes = contour->Get("integrationNodes",ParamNode::APPEND);
  PtrParamNode weights = contour->Get("integrationWeights",ParamNode::APPEND);
  // amount of nodes
  int numNodes = integrationNodes_.GetSize()/2;
  // setup vectors
  StdVector<Double> integrationNodesReal;
  StdVector<Double> integrationNodesImag;
  StdVector<Double> integrationWeightsReal;
  StdVector<Double> integrationWeightsImag;
  integrationNodesReal.Resize(numNodes);
  integrationNodesImag.Resize(numNodes);
  integrationWeightsReal.Resize(numNodes);
  integrationWeightsImag.Resize(numNodes);
  // fill vectors
  for(int i = 0; i < numNodes; i++) {
    integrationNodesReal[i] = integrationNodes_[i*2];
    integrationNodesImag[i] = integrationNodes_[i*2+1];
    integrationWeightsReal[i] = integrationWeights_[i*2];
    integrationWeightsImag[i] = integrationWeights_[i*2+1];
  }
  // write in info
  nodes->Get("real")->SetValue(integrationNodesReal.ToString());
  nodes->Get("imag")->SetValue(integrationNodesImag.ToString());
  weights->Get("real")->SetValue(integrationWeightsReal.ToString());
  weights->Get("imag")->SetValue(integrationWeightsImag.ToString());
}


} // end of namespace
