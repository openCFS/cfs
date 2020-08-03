#include "PhistEigenSolver.hh"
#include "DataInOut/Logging/LogConfigurator.hh"

using std::complex;

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
    eigenSolverName_ = EigenSolverType::PHIST;

    comm_ = NULL;

    scale_B_ = xml_->Get("scale_B")->As<bool>();

    which.SetName("phist_EeigSort:which");
    which.Add(phist_LM, "lm");
    which.Add(phist_SM, "sm");
    which.Add(phist_LR, "lr");
    which.Add(phist_SR, "sr");

    linSolv.SetName("phist_ElinSolv:linSolv");
    linSolv.Add(phist_GMRES, "gmres");
    linSolv.Add(phist_MINRES, "minres");
    linSolv.Add(phist_QMR, "qmr");
    linSolv.Add(phist_BICGSTAB, "bicgstab");
    linSolv.Add(phist_CARP_CG, "carp_cg");

    BaseEigenSolver::PostInit();
  }

  PhistEigenSolver::~PhistEigenSolver()
  {
    int iflag = 0;

    assert(A_ != NULL);
    if(complex_)
      phist::kernels<Complex>::sparseMat_delete(static_cast<phist_ZsparseMat_ptr>(A_), &iflag);
    else
      phist::kernels<double>::sparseMat_delete(static_cast<phist_DsparseMat_ptr>(A_), &iflag);

    A_ = NULL;

    assert(B_ != NULL);
    if(complex_)
      phist::kernels<Complex>::sparseMat_delete(static_cast<phist_ZsparseMat_ptr>(B_), &iflag);
    else
      phist::kernels<double>::sparseMat_delete(static_cast<phist_DsparseMat_ptr>(B_), &iflag);
    B_ = NULL;

    ghost_finalize();
  }

  void PhistEigenSolver::Setup(const BaseMatrix & mat, UInt numFreq, Double freqShift, bool sort)
  {
    LOG_DBG(pes) << "PES:S(stiff)";

    // Set flag for indicating a non-quadratic problem
    isQuadratic_ = false;

    SetupCommon(numFreq, freqShift); // no bloch

    ToInfo();
  }

  template<class TYPE>
  void PhistEigenSolver::SaveModes(typename phist::types<TYPE>::mvec_ptr X, int nEig)
  {
    // skip calculation of eigenvector residuals
    // download X from GPU if applicable
    int iflag = 0;
    phist::kernels<TYPE>::mvec_from_device(X, &iflag);

    // get pointer to row/col major block of vectors
    TYPE* xval = NULL; // will be set to memory owned by phist
    phist_lidx lda, nloc;
    phist::kernels<TYPE>::mvec_my_length(X, &nloc, &iflag);
    phist::kernels<TYPE>::mvec_extract_view(X, &xval, &lda, &iflag);
    LOG_DBG(pes) << "SM: lda=" << lda << " nloc=" << nloc << " nEig=" << nEig;
    LOG_DBG3(pes) << "SM: xval=" << ToString(xval, nloc * lda);

    if(lda != nEig)
      EXCEPTION("phist could only solve for " << lda << " of " << nEig << " modes");
    if(nloc == 0)
      EXCEPTION("phist could could not compute the eigenmodes");

    // this is how one can extract the eigenvectors.
    mode_.Resize(nEig, nloc);
    // we rescale x by sqrt(scaling) as we have x^T (scale*B) x = 1
    double rescale = std::sqrt(scale_B_val_);
    LOG_DBG(pes) << "SM: mode_.rows=" << mode_.GetNumRows() << " .cols=" << mode_.GetNumCols();
    assert((int) sort_idx_.GetSize() == nEig);
    for (int j = 0; j < nEig; j++)
      for (int i = 0; i < nloc; i++)
        mode_[sort_idx_[j]][i] = xval[i * lda + j] * rescale ; // mode_ is always complex!
    LOG_DBG3(pes) << "SM: mode=" << mode_.ToString(2,true);
  }

  template<class TYPE>
  void PhistEigenSolver::Setup(const BaseMatrix& A, const BaseMatrix& B, bool isHermitian)
  {
    shared_ptr<Timer> timer = info_->Get(ParamNode::SUMMARY)->Get("phist_setup/timer")->AsTimer();
    timer->Start();

    // TODO one of both shall be enough!
    hermitian_ = isHermitian;
    isBloch_ = isHermitian;
    complex_ = boost::is_complex<TYPE>::value;
    assert(!(hermitian_ && !complex_));

    // we scale the B-Matrix as suggested by Jonas:
    // A*x=lambda*sigma*B*x, with z.B. sigma=1.0/B(1,1) and then rescale
    const SparseOLASMatrix<TYPE>* olas = dynamic_cast<const SparseOLASMatrix<TYPE>* >(&B);
    assert(olas != NULL);
    TYPE b11 = olas->GetAvgDiag();

    assert(((Complex) b11).real() > 1e-15);
    scale_B_val_ = scale_B_ ? ((Complex) (1./b11)).real() : 1.0;
    LOG_DBG(pes) << "PES:S b11=" << b11 << " -> B-scale=" << scale_B_val_;

    // for include stuff issues, A_ and B_ attributes are not of full type
    A_ = InitMatrix<TYPE>(A, &A_, 1.0);
    // TODO: in the Bloch-case the B-Matrix does not change and can be reused!
    B_ = InitMatrix<TYPE>(B, &B_, scale_B_val_);

    timer->Stop();
  }


  template<class TYPE>
  void PhistEigenSolver::CalcEigenValues(BaseVector &sol, BaseVector &err, unsigned int numFreq, double freqShift)
  {
    shared_ptr<Timer> solve = info_->Get(ParamNode::SUMMARY)->Get("phist_solve/timer")->AsTimer();
    solve->Start();

    // for include stuff issues, A_ and B_ attributes are not of full type
    assert(A_ != NULL && B_ != NULL); // Setup() has to been called first!
    typename phist::types<TYPE>::sparseMat_ptr A = (typename phist::types<TYPE>::sparseMat_ptr) A_;
    typename phist::types<TYPE>::sparseMat_ptr B = (typename phist::types<TYPE>::sparseMat_ptr) B_;

    SetupCommon(numFreq, freqShift);

    // setup eigenvalue problem - taken from subspacejada (jacobi davidson)
    // create an operator from A
    typename phist::types<TYPE>::linearOp_ptr opA = new typename phist::types<TYPE>::linearOp();

    // we need the domain map of the matrix
    phist_const_map_ptr map = NULL;
    int iflag = 0;
    phist::kernels<TYPE>::sparseMat_get_domain_map(A,&map,&iflag);
    LOG_DBG(pes) << "sparseMat_get_domain_map -> " << iflag;

    typename phist::types<TYPE>::linearOp_ptr opB = new typename phist::types<TYPE>::linearOp();
    phist::core<TYPE>::linearOp_wrap_sparseMat(opB,B,&iflag);
    LOG_DBG(pes) << "linearOp_wrap_sparseMat -> " << iflag;

    phist::core<TYPE>::linearOp_wrap_sparseMat_pair(opA,A,B,&iflag);
    LOG_DBG(pes) << "linearOp_wrap_sparseMat_pair -> " << iflag;

    int prob_size = numFreq+opts_.blockSize-1;
    LOG_DBG(pes) << "prob_size -> " << prob_size;

    // setup necessary vectors and matrices for the schur form
    typename phist::types<TYPE>::mvec_ptr Q = NULL;
    phist::kernels<TYPE>::mvec_create(&Q,map,prob_size,&iflag);
    LOG_DBG(pes) << "mvec_create -> " << iflag;

    typename phist::types<TYPE>::sdMat_ptr R = NULL;
    phist::kernels<TYPE>::sdMat_create(&R,prob_size,prob_size,comm_,&iflag);
    LOG_DBG(pes) << "sdMat_create -> " << iflag;


    resNorm_.Resize(prob_size);
    ev_.Resize(prob_size);

    // setup start vector (currently to (1 0 1 0 .. ) )
    typename phist::types<TYPE>::mvec_ptr v0 = NULL;
    phist::kernels<TYPE>::mvec_create(&v0,map,1,&iflag);
    LOG_DBG(pes) << "mvec_create -> " << iflag;

    phist::kernels<TYPE>::mvec_put_value(v0,1.0,&iflag);

    // skip residual calculation from subspacejada, not necessary for us
    opts_.v0=v0;

    // assume phist_NO_PRECON

    int nIter=opts_.maxIters;
    int nEig = (int) numFreq;

    // The actual calculation !!! - standard Ritz values for exterior eigenvalues.
    // QR-factorization, Q=orthogonal ev basis , R=upper triangular matrix, ev are on the diagonal,
    // A*Q = B*Q*R holds
    assert(opts_.how != phist_HARMONIC);

    phist::jada<TYPE>::subspacejada(opA, opB, opts_, Q, R, ev_.GetPointer(), resNorm_.GetPointer(), &nEig, &nIter, &iflag);
    LOG_DBG(pes) << "subspacejada -> nEig=" << nEig << " nIter=" << nIter << " ev=" << ev_.ToString() << " -> " << iflag;
    assert(iflag >= 0);
    last_iter_ = nIter;
    sum_iter_ += nIter;
    count_iter_++;



    // skip calculation of real residual, res = AQ - BQR

    // compute eigenvectors X: A*X=X*D and diagonal matrix D with eigenvalues
    // (for checking the sorting only).
    typename phist::types<TYPE>::mvec_ptr X=NULL;
    // we want the information back from cuda in SaveModes(
    iflag = PHIST_MVEC_REPLICATE_DEVICE_MEM; // ignored if not cuda
    phist::kernels<TYPE>::mvec_create(&X,map,nEig,&iflag);
    assert(iflag >= 0); // iflag is reset to 0 for success

    // FIXME phist ComputeEigenVectors moved to core but seg faults in CPU
    phist::core<TYPE>::ComputeEigenvectors(Q,R,X,&iflag);
    assert(iflag >= 0);

    // rescale ev_
    complex<double> rescale(scale_B_val_, 1.0);
    LOG_DBG(pes) << "S: rescale ev_ by " << rescale << " org_ev=" << ev_.ToString();
    ev_.ScalarMult(rescale);
    LOG_DBG(pes) << "S: scaled ev_=" << ev_.ToString();


    double pi_pi = 2.0*M_PI;
    assert(close(pi_pi, 8.0*atan(1.0)));

    // we need to sort our stuff
    StdVector<double> freq(numFreq_);
    for(unsigned int i =0; i < numFreq_; i++) // rigid modes are approx zero and can be negative -> abs to avoid NaN
      freq[i] = sqrt(std::abs(ev_[i].real()))/pi_pi; // this is the omega^2 -> f transformation from ArpackEigenSolver :(
    SetupSortIdx(freq);

    // skip calculation of eigenvector residuals
    // save the modes in a sorted fashion
    SaveModes<TYPE>(X, nEig);

    // store solution and error
    sol.Resize(numFreq_);
    err.Resize(numFreq_);

    assert(sol.GetEntryType() == complex_ ? BaseMatrix::COMPLEX : BaseMatrix::DOUBLE);
    assert(err.GetEntryType() == BaseMatrix::DOUBLE);
    // no set results
    assert(!(complex_ && sol.GetEntryType() != BaseMatrix::COMPLEX));
    Vector<TYPE>& solConverted = dynamic_cast<Vector<TYPE>&>(sol);

    for(unsigned int i =0; i < numFreq_; i++) {
      solConverted[i] = freq[sort_idx_[i]];
      err.SetEntry(i, resNorm_[sort_idx_[i]]);
    }

    //Delete all things declared using new or create.
    delete [] opA;
    delete [] opB;

    phist::kernels<TYPE>::mvec_delete(Q, &iflag);
    phist::kernels<TYPE>::mvec_delete(X, &iflag);
    phist::kernels<TYPE>::sdMat_delete(R, &iflag);

    solve->Stop();
    
    ToInfo();
  }


  void PhistEigenSolver::SetupSortIdx(const StdVector<double>& freq)
  {
    // to be called after within each CalcEigenFrequencies()
    sort_idx_.Resize(freq.GetSize());

    // we sort with a std::pair<||ev||, org_idx>
    std::vector<ev_idx> org;
    org.reserve(freq.GetSize());
    for(unsigned int i = 0; i < freq.GetSize(); i++)
      org.push_back(std::make_pair(freq[i], i));

    std::sort(org.begin(), org.end(), PhistEigenSolver::comperator);

    for(unsigned int i = 0; i < freq.GetSize(); i++)
      sort_idx_[i] = org[i].second;

    std::cout << "phist sort: " << sort_idx_.ToString() << std::endl;
  }


  void PhistEigenSolver::SetupCommon(unsigned int numFreq, Double freqShift)
  {
    LOG_DBG(pes) << "PES:SC sym=" << sym_;

    // Save frequency parameters
    numFreq_ = numFreq;
    freqShift_ = freqShift;

    isBloch_ = hermitian_;
    sort_ = true;
    // phist parameters usually read from opts-standard.txt

    // fill the jadaOpts struct to pass settings to the solver
    phist_jadaOpts_setDefaults(&opts_);

    opts_.symmetry = sym_ ? phist_HERMITIAN : phist_GENERAL; // phist_HERMITIAN also for real symmetric
    opts_.numEigs   = numFreq;
    opts_.which     = which.Parse(xml_->Get("which")->As<string>());
    opts_.convTol   = xml_->Get("convTol")->As<double>();
    opts_.blockSize = xml_->Get("blockSize")->As<int>();
    opts_.maxIters  = xml_->Get("maxIter")->As<int>();
    opts_.minBas    = numFreq + 2 * opts_.blockSize; // 28
    opts_.maxBas    = xml_->Get("maxBas")->As<int>(); //opts.minBas + 10 * opts.blockSize; // 60

    // parameters for the linear system stuff. The element is optional
    bool is = xml_->Has("innerSolv");

    opts_.innerSolvBlockSize = is ? xml_->Get("innerSolv/blockSize")->As<int>() : -1;
    opts_.innerSolvType      = is ? linSolv.Parse(xml_->Get("innerSolv/type")->As<string>()) : phist_MINRES;
    opts_.innerSolvMaxBas    = is ? xml_->Get("innerSolv/maxBas")->As<int>() : -1;
    opts_.innerSolvMaxIters  = is ? xml_->Get("innerSolv/maxIter")->As<int>() : 10;
    opts_.innerSolvRobust    = is ? xml_->Get("innerSolv/robust")->As<int>() : 1;
    opts_.innerSolvBaseTol   = is ? xml_->Get("innerSolv/baseTol")->As<double>() : 0.1;
  }


  bool PhistEigenSolver::IsSymmetric(const BaseMatrix& cfs) const
  {
    const StdMatrix* sm = dynamic_cast<const StdMatrix*>(&cfs);
    assert(sm != NULL);
    return sm->GetStorageType() == StdMatrix::SPARSE_SYM;
  }

  void PhistEigenSolver::ToInfo()
  {
    PtrParamNode setup = info_->Get(ParamNode::HEADER);
    setup->Get("quadratic")->SetValue(isQuadratic_);
    setup->Get("bloch")->SetValue(isBloch_);
    setup->Get("numFreq")->SetValue(numFreq_);

    PtrParamNode phist = setup->Get("phist");
    phist->Get("convTol")->SetValue(opts_.convTol);
    phist->Get("blockSize")->SetValue(opts_.blockSize);
    phist->Get("maxIters")->SetValue(opts_.maxIters);
    phist->Get("minBas")->SetValue(opts_.minBas);
    phist->Get("maxBas")->SetValue(opts_.maxBas);
    phist->Get("which")->SetValue(which.ToString(opts_.which));

    PtrParamNode is = phist->Get("innerSolv");
    is->Get("type")->SetValue(linSolv.ToString(opts_.innerSolvType));
    is->Get("baseTol")->SetValue(opts_.innerSolvBaseTol);
    is->Get("blockSize")->SetValue(opts_.innerSolvBlockSize);
    is->Get("maxBas")->SetValue(opts_.innerSolvMaxBas);
    is->Get("maxIter")->SetValue(opts_.innerSolvMaxIters);
    is->Get("robust")->SetValue(opts_.innerSolvRobust);
    is->Get("blockSize")->SetValue(opts_.innerSolvBlockSize);

    PtrParamNode proc = info_->Get(ParamNode::PROCESS);
    PtrParamNode curr = proc->Get("phist", progOpts->DoDetailedInfo() ? ParamNode::APPEND : ParamNode::INSERT);
    curr->Get("analysis_id")->SetValue(domain->GetDriver()->GetAnalysisId().ToString());
    curr->Get("iter")->SetValue(last_iter_);
    curr->Get("avg_iter")->SetValue(sum_iter_ / (double) count_iter_);
    curr->Get("B_scaling")->SetValue(scale_B_val_);
    curr->Get("sorting")->SetValue(sort_idx_.ToString());

  }

  void PhistEigenSolver::Setup(const BaseMatrix& stiffMat, const BaseMatrix& massMat, UInt numFreq, Double freqShift, bool sort, bool bloch)
  {
    LOG_DBG(pes) << "PES:S(stiff, mass)";


    this->numFreq_ = numFreq;
    this->freqShift_ = freqShift;

    // todo: shall be removed with new interface
    this->sym_ = IsSymmetric(stiffMat);

    // Set flag for indicating a non-quadratic problem
    isQuadratic_ = false;

    if(stiffMat.GetEntryType() == BaseMatrix::DOUBLE)
      Setup<double>(stiffMat, massMat, bloch);
    else
      Setup<Complex>(stiffMat, massMat, bloch);
  }


  void PhistEigenSolver::Setup(const BaseMatrix & stiffMat, const BaseMatrix & massMat, const BaseMatrix & dampMat,
                                UInt numFreq, double freqShift, bool sort)
  {
    LOG_DBG(pes) << "PES:S(stiff, mass, damp)";

    assert(false);
  }

  void PhistEigenSolver::CalcEigenFrequencies(BaseVector &sol, BaseVector &err)
  {
    // remove this, this is the old interface!
    if(hermitian_ == true)
      CalcEigenValues<Complex>(sol, err, numFreq_, freqShift_);
    else
      CalcEigenValues<double>(sol, err, numFreq_, freqShift_);
  }

  void PhistEigenSolver::CalcConditionNumber(const BaseMatrix& mat, Double& condNumber, Vector<Double>& evs, Vector<Double>& err)
  {
    // Set flag for indicating a non-quadratic problem
    isQuadratic_ = false;
    // NOTE: Hard coded as true!!!
    assert(false);

  }

  void PhistEigenSolver::GetEigenMode(unsigned int modeNr, Vector<Complex>& mode, bool right)
  {
    assert(right == true); // no idea what this means ?!
    LOG_DBG(pes) << "GEM: modeNr=" << modeNr << " right=";
    assert(mode_.GetNumRows() > 0 && mode_.GetNumCols() > 0);
    if(modeNr >= mode_.GetNumRows())
      EXCEPTION("request mode " << (modeNr+1) << " (1-based) not available");

    mode.Resize(mode_.GetNumCols());
    mode_.GetRow(mode, modeNr);
    LOG_DBG3(pes) << "GEM -> " << mode.ToString();
  }
}
