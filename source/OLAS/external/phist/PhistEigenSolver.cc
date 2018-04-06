


#include "PhistEigenSolver.hh"


//#include "phist_subspacejada.h"

using std::complex;


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

    comm_ = NULL;

    scale_mass_ = xml_->Get("scale_mass")->As<bool>();

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
  }

  PhistEigenSolver::~PhistEigenSolver()
  {
    int iflag = 0;
    phist::kernels<double>::sparseMat_delete(A_, &iflag);
    LOG_DBG(pes) << "~PES: del A -> " << iflag;
    // assert(iflag = 0);
    A_ = NULL;

    phist::kernels<double>::sparseMat_delete(B_, &iflag);
    LOG_DBG(pes) << "~PES: del B -> " << iflag;
    // assert(iflag = 0);
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

    SetupCommon(IsSymmetric(mat), numFreq, freqShift, sort, false); // no bloch

    ToInfo();
  }

  void PhistEigenSolver::SaveModes(phist::types<double>::mvec_ptr X, int nEig)
  {
    // skip calculation of eigenvector residuals
    // download X from GPU if applicable
    int iflag = 0;
    phist::kernels<double>::mvec_from_device(X, &iflag);
    assert(iflag == 0);
    // get pointer to row/col major block of vectors
    double* xval = NULL; // will be set to memory owned by phist
    phist_lidx lda, nloc;
    phist::kernels<double>::mvec_my_length(X, &nloc, &iflag);
    assert(iflag == 0);
    phist::kernels<double>::mvec_extract_view(X, &xval, &lda, &iflag);
    assert(iflag == 0);
    // this is how one can extract the eigenvectors.
    // TODO: move all of this to CalcEigenFrequencies and other functions below
    mode_.Resize(nEig, nloc);
    for (int j = 0; j < nEig; j++)
    {
      for (int i = 0; i < nloc; i++)
      {
        mode_[j][i] = xval[j * lda + i];
      }
    }
  }

  void PhistEigenSolver::Setup(const BaseMatrix& stiffMat, const BaseMatrix& massMat, UInt numFreq, Double freqShift, bool sort, bool bloch)
  {
    LOG_DBG(pes) << "PES:S(stiff, mass)";

    shared_ptr<Timer> setup = info_->Get(ParamNode::SUMMARY)->Get("phist_setup/timer")->AsTimer();
    setup->Start();

    // Set flag for indicating a non-quadratic problem
    isQuadratic_ = false;

    SetupCommon(IsSymmetric(stiffMat), numFreq, freqShift, sort, bloch);

    // we scale the B-Matrix as suggested by Jonas:
    // A*x=lambda*sigma*B*x, with z.B. sigma=1.0/B(1,1) and then rescale
    double b11;
    dynamic_cast<const StdMatrix*>(&massMat)->GetDiagEntry(0, b11);
    assert(b11 > 1e-15);
    double scale = scale_mass_ ? 1./b11 : 1.0;
    LOG_DBG(pes) << "PES:S b11=" << b11 << " -> B-scale=" << scale;

    shared_ptr<Timer> setup_phistMatrix = info_->Get(ParamNode::SUMMARY)->Get("setup_phistMatrix/timer")->AsTimer();
    setup_phistMatrix->SetSub();
    setup_phistMatrix->Start();


    // for include stuff issues, A_ and B_ attributes are not of full type
    phist::types<double>::sparseMat_ptr A = (phist::types<double>::sparseMat_ptr) InitMatrix(stiffMat, &A_, 1.0);
    phist::types<double>::sparseMat_ptr B = (phist::types<double>::sparseMat_ptr) InitMatrix(massMat, &B_, scale);


    // setup eigenvalue problem - taken from subspacejada (jacobi davidson)
    // create an operator from A
    phist::types<double>::linearOp_ptr opA = new phist::types<double>::linearOp(); // TODO needs to be specialized

    // we need the domain map of the matrix
    phist_const_map_ptr map = NULL;
    int iflag = 0;
    phist::kernels<double>::sparseMat_get_domain_map(A,&map,&iflag);
    LOG_DBG(pes) << "sparseMat_get_domain_map -> " << iflag;
    assert(iflag == 0);

    phist::types<double>::linearOp_ptr opB = new phist::types<double>::linearOp();
    phist::core<double>::linearOp_wrap_sparseMat(opB,B,&iflag);
    LOG_DBG(pes) << "linearOp_wrap_sparseMat -> " << iflag;
    assert(iflag == 0);

    phist::core<double>::linearOp_wrap_sparseMat_pair(opA,A,B,&iflag);
    LOG_DBG(pes) << "linearOp_wrap_sparseMat_pair -> " << iflag;
    assert(iflag == 0);

    setup_phistMatrix->Stop();


    int prob_size = numFreq+opts_.blockSize-1;
    LOG_DBG(pes) << "prob_size -> " << prob_size;

    // setup necessary vectors and matrices for the schur form
    phist::types<double>::mvec_ptr Q = NULL;
    phist::kernels<double>::mvec_create(&Q,map,prob_size,&iflag);
    LOG_DBG(pes) << "mvec_create -> " << iflag;
    assert(iflag == 0);

    phist::types<double>::sdMat_ptr R = NULL;
    phist::kernels<double>::sdMat_create(&R,prob_size,prob_size,comm_,&iflag);
    LOG_DBG(pes) << "sdMat_create -> " << iflag;
    assert(iflag == 0);

    resNorm_.Resize(prob_size);
    ev_.Resize(prob_size);

    // setup start vector (currently to (1 0 1 0 .. ) )
    phist::types<double>::mvec_ptr v0 = NULL;
    phist::kernels<double>::mvec_create(&v0,map,1,&iflag);
    LOG_DBG(pes) << "mvec_create -> " << iflag;
    assert(iflag == 0);

    phist::kernels<double>::mvec_put_value(v0,1.0,&iflag);
    assert(iflag == 0);

    // skip residual calculation from subspacejada, not necessary for us
    opts_.v0=v0;

    // assume phist_NO_PRECON

    int nIter=opts_.maxIters;
    int nEig = (int) numFreq;

    // The actual calculation !!! - standard Ritz values for exterior eigenvalues.
    // QR-factorization, Q=orthogonal ev basis , R=upper triangular matrix, ev are on the diagonal,
    // A*Q = B*Q*R holds
    assert(opts_.how != phist_HARMONIC);

    phist::jada<double>::subspacejada(opA, opB, opts_, Q, R, ev_.GetPointer(), resNorm_.GetPointer(), &nEig, &nIter, &iflag);
    LOG_DBG(pes) << "subspacejada -> nEig=" << nEig << " nIter=" << nIter << " ev=" << ev_.ToString() << " -> " << iflag;
    assert(iflag >= 0);

    // skip calculation of real residual, res = AQ - BQR

    // compute eigenvectors X: A*X=X*D and diagonal matrix D with eigenvalues
    // (for checking the sorting only).
    phist::types<double>::mvec_ptr X=NULL;
    // we want the information back from cuda in SaveModes(
    iflag = PHIST_MVEC_REPLICATE_DEVICE_MEM; // ignored if not cuda
    phist::kernels<double>::mvec_create(&X,map,nEig,&iflag);
    assert(iflag >= 0); // iflag is reset to 0 for success

    setup->Stop();

    shared_ptr<Timer> solve = info_->Get(ParamNode::SUMMARY)->Get("phist_solve/timer")->AsTimer();
    solve->Start();

    phist::jada<double>::ComputeEigenvectors(Q,R,X,&iflag);
    assert(iflag >= 0);

    // rescale ev_
    complex<double> rescale(scale, 1.0);
    LOG_DBG(pes) << "S: rescale ev_ by " << rescale << " org_ev=" << ev_.ToString();
    ev_.ScalarMult(rescale);
    LOG_DBG(pes) << "S: scaled ev_=" << ev_.ToString();


    // skip calculation of eigenvector residuals
    SaveModes(X, nEig);
    solve->Stop();

    //Delete all things declared using new or create.
    delete []  opA;
    delete [] opB;

    phist::kernels<double>::sparseMat_delete(A, &iflag);
    phist::kernels<double>::sparseMat_delete(B, &iflag);
    phist::kernels<double>::mvec_delete(Q, &iflag);
    phist::kernels<double>::mvec_delete(X, &iflag);
    phist::kernels<double>::sdMat_delete(R, &iflag);
    
    ToInfo();
  }


  void PhistEigenSolver::Setup(const BaseMatrix & stiffMat, const BaseMatrix & massMat, const BaseMatrix & dampMat,
                                UInt numFreq, double freqShift, bool sort)
  {
    LOG_DBG(pes) << "PES:S(stiff, mass, damp)";

    // Set flag for indicating a non-quadratic problem
    isQuadratic_ = true;

    SetupCommon(IsSymmetric(stiffMat), numFreq, freqShift, sort, false); // no bloch

    assert(false);
  }


  void PhistEigenSolver::SetupCommon(bool sym, unsigned int numFreq, Double freqShift, bool sort, bool bloch)
  {
    LOG_DBG(pes) << "PES:SC sym=" << sym;

    // Save frequency parameters
    numFreq_ = numFreq;
    freqShift_ = freqShift;

    isBloch_ = bloch;
    sort_ = sort;


    // phist parameters usually read from opts-standard.txt

    // fill the jadaOpts struct to pass settings to the solver
    phist_jadaOpts_setDefaults(&opts_);

    opts_.symmetry = sym ? phist_HERMITIAN : phist_GENERAL; // phist_HERMITIAN also for real symmetric
    opts_.numEigs   = numFreq;
    opts_.which     = which.Parse(xml_->Get("which")->As<string>());
    opts_.convTol   = xml_->Get("convTol")->As<double>();
    opts_.blockSize = xml_->Get("blockSize")->As<int>();
    opts_.maxIters  = xml_->Get("maxIter")->As<int>();
    opts_.minBas    = numFreq + 2 * opts_.blockSize; // 28
    opts_.maxBas    = 60; //opts.minBas + 10 * opts.blockSize; // 60

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
  }


  void PhistEigenSolver::CalcEigenFrequencies(BaseVector &sol, BaseVector &err)
  {
    sol.Resize(numFreq_);
    err.Resize(numFreq_);

    assert(sol.GetEntryType() == BaseMatrix::DOUBLE);
    assert(err.GetEntryType() == BaseMatrix::DOUBLE);

    for(unsigned int i =0; i < numFreq_; i++) {
      // rigid modes are approx zero and can be negative -> abs to avoid NaN
      sol.SetEntry(i, sqrt(std::abs(ev_[i].real()))/(8.0*atan(1.0))); // this is the stupid lambda = omega^2 -> f transformation from ArpackEigenSolver :(
      err.SetEntry(i, resNorm_[i]);
    }
  }

  void PhistEigenSolver::CalcConditionNumber(const BaseMatrix& mat, Double& condNumber, Vector<Double>& evs, Vector<Double>& err)
  {
    // Set flag for indicating a non-quadratic problem
    isQuadratic_ = false;
    // NOTE: Hard coded as true!!!
    assert(false);

  }

  void PhistEigenSolver::GetEigenMode(UInt modeNr, Vector<Complex> & mode)
  {
    if(modeNr >= mode_.GetNumRows())
    {
      std::cout << "PhistEigenSolver::GetEigenMode(" << modeNr << ") not in 0 ... " << mode_.GetNumRows() << std::endl;
      mode.Init(0);
      return;
    }
    assert(modeNr < mode_.GetNumRows()); // assume 0-based

    mode.Resize(mode_.GetNumCols());
    mode_.GetCol(mode, modeNr);
  }
}
