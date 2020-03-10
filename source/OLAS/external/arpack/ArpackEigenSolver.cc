#include <limits>
#include <string.h>

#include "MatVec/StdMatrix.hh"
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

#include "ArpackEigenSolver.hh"

DEFINE_LOG(aes, "arpackEigenSolver")


namespace CoupledField {

  ArpackEigenSolver::ArpackEigenSolver( shared_ptr<SolStrategy> strat,
                                        PtrParamNode xml,
                                        PtrParamNode solverList,
                                        PtrParamNode precondList,
                                        PtrParamNode eigenInfo)
    : BaseEigenSolver( strat, xml, solverList, precondList, eigenInfo ){

    interface_    = NULL;
    matrixA_      = NULL;
    matrixB_      = NULL;
    matrixD_      = NULL;
    zStiff_       = NULL;
    zMass_        = NULL;
    zDamp_        = NULL;
    arpackSolver_ = NULL;
    solver_       = NULL;
    precond_      = NULL;
    which_        = NULL;
    xml_          = xml;
    
    isGeneralized_ = false;
    shiftAndInvert_ = false;
  }

  ArpackEigenSolver::~ArpackEigenSolver() {

    delete solver_;
    delete precond_;
    delete interface_;
    delete arpackSolver_;

    delete zStiff_;
    delete zDamp_;
    delete zMass_;

    delete[] which_;

    matrixA_ = NULL;
    matrixB_ = NULL;
    matrixD_ = NULL;
    zStiff_ = NULL;
    zMass_  = NULL;
    zDamp_  = NULL;
  }
  void ArpackEigenSolver::Setup(const BaseMatrix & mat, UInt numFreq, Double freqShift, bool sort)
  {
    // Set flag for indicating a non-quadratic problem
    isQuadratic_ = false;
    isGeneralized_ = false;
    isBloch_ = false;
    sort_ = false;

    // NOTE: Hard coded as true!!!
    shiftAndInvert_ = true;

    // Copy matrix references and determine size of system
    matrixA_ = & dynamic_cast<const StdMatrix&>(mat);

    UInt size = matrixA_->GetNumRows();

    // Save frequency parameters
    numFreq_ = numFreq;
    freqShift_ = freqShift;

    // Create matrix interface for arpack
    interface_ = new ArpackMatInterface( matrixA_, shiftAndInvert_, freqShift_ );

    // Create solver class
    arpackSolver_ = new ArpackSolver(xml_);

    // Check 'which'-settings regarding the type of eigenvalues searched for
    std::string whichString = "LM";
    xml_->GetValue("which", whichString, ParamNode::INSERT);
    
    which_ = new char[whichString.size()+1];
    strncpy(which_, whichString.c_str(), whichString.size()+1 );

    // Create the solver dependent on the problem type ( regular / shiftAndInvert )
    arpackSolver_->Setup(interface_, size, numFreq_, freqShift_, which_ , (char*) "I", shiftAndInvert_, isBloch_);

    // Create solver
    solver_ = GenerateSolverObject( *matrixA_, solStrat_, solverList_, info_ );

    // Perform check, if matrix is std or sbm
    if ( matrixA_->GetStructureType() == BaseMatrix::SPARSE_MATRIX ) {
      const StdMatrix & mat = dynamic_cast< const StdMatrix &>( *matrixA_ );

      // Create preconditioner
      precond_ = GenerateStdPrecondObject( mat, solStrat_->GetPrecondId(), precondList_, info_ );
    } else {
      EXCEPTION( "No preconditioner available for SBM-matrices!" );
    }

    // Setup matrixinterface
    interface_->Setup( solver_, precond_ );

    ToInfo();
  }
  
  
  void ArpackEigenSolver::Setup(const BaseMatrix& stiffMat, const BaseMatrix& massMat, UInt numFreq, Double freqShift, bool sort, bool bloch)
  {
    // Set flag for indicating a non-quadratic problem
    isQuadratic_ = false;
    isGeneralized_ = true;
    isBloch_ = bloch;
    sort_ = sort;

    // Copy matrix references and determine size of system
    matrixA_ = & dynamic_cast<const StdMatrix&>(stiffMat);

    // bloch works only for non-symmetric matrices as the stiffness matrix needs to be Hermitian.
    // At least Pardiso can be used with <pardiso> <hermitean>yes</hermitean> </pardiso>

    matrixB_ = & dynamic_cast<const StdMatrix&>(massMat);

    UInt size = matrixA_->GetNumRows();

    // Save frequency parameters
    numFreq_ = numFreq;
    freqShift_ = freqShift;

    // NOTE: Hard coded as true!!!
    shiftAndInvert_ = true;

    // Create matrix interface for arpack, it might already exist in the bloch mode case
    if(interface_ == NULL)
      interface_ = new ArpackMatInterface( matrixA_, matrixB_, shiftAndInvert_, freqShift_ );

    LOG_DBG(aes) << "S(S,M) shift=" << freqShift_ << " s&i=" << shiftAndInvert_;

    // Create solver class, it might already exist in the bloch mode case
    if(arpackSolver_  == NULL)
      arpackSolver_ = new ArpackSolver(xml_);

    // Check 'which'-settings regarding the type of eigenvalues searched for
    std::string whichString = "LM";
    xml_->GetValue("which", whichString, ParamNode::INSERT );
    
    if(which_ != NULL) // called multiple times
      delete[] which_;

    which_ = new char[whichString.size()+1];
    strncpy(which_, whichString.c_str(), whichString.size()+1 );

    // Create the solver dependent on the problem type ( regular / shiftAndInvert )
    arpackSolver_->Setup( interface_, size, numFreq_, freqShift_, which_ , (char*) "G",
    		              shiftAndInvert_, isBloch_);

    // Create solver - for every bloch wave vector :( Make sure it will be deleted!
    if(solver_ != NULL)
      delete solver_;
    solver_ = GenerateSolverObject( *matrixB_, solStrat_,  solverList_, info_);

    // Perform check, if matrix is std or sbm
    if ( matrixA_->GetStructureType() == BaseMatrix::SPARSE_MATRIX )
    {
      const StdMatrix & mat = dynamic_cast< const StdMatrix &>( *matrixB_ );
      // Create preconditioner
      if(precond_ == NULL)
        precond_ = GenerateStdPrecondObject( mat, solStrat_->GetPrecondId(), precondList_, info_);
    }
    else
    {
      EXCEPTION( "No preconditioner available for SBM-matrices!" );
    }

    // Setup matrixinterface
    interface_->Setup( solver_, precond_ );

    ToInfo();

  }

  void ArpackEigenSolver::ToInfo()
  {
    PtrParamNode setup = info_->Get(ParamNode::HEADER);
    arpackSolver_->ToInfo(setup);

    setup->Get("generalized")->SetValue(isGeneralized_);
    setup->Get("quadratic")->SetValue(isQuadratic_);
    setup->Get("bloch")->SetValue(isBloch_);
  }


  void ArpackEigenSolver::Setup(const BaseMatrix & stiffMat, const BaseMatrix & massMat, const BaseMatrix & dampMat,
                                UInt numFreq, double freqShift, bool sort) {

    // Set flag for indicating a non-quadratic problem
    isQuadratic_ = true;
    isGeneralized_ = true;
    isBloch_ = false;
    sort_ = false;

    // Copy matrix references and convert them to StdMatrices
    matrixA_ = & dynamic_cast<const StdMatrix&>(stiffMat);

    matrixB_ = & dynamic_cast<const StdMatrix&>(massMat);

    matrixD_ = & dynamic_cast<const StdMatrix&>(dampMat);

    // need twice the original system size and twice the number of frequencies
    UInt size = 2*matrixA_->GetNumRows();
    //numFreq_  = 2*numFreq;
    numFreq_  = numFreq;

    // currently no shifts are allowed
    freqShift_ = freqShift;
    
    // NOTE: Hard coded as true!!!
    shiftAndInvert_ = true;

    // Now generate a new complex matrix
    BaseMatrix::EntryType eType = BaseMatrix::COMPLEX;
    BaseMatrix::StorageType sType = matrixA_->GetStorageType();
    Integer nrows = matrixA_->GetNumRows();
    Integer ncols = matrixA_->GetNumCols();
    Integer fill = matrixA_->GetNnz(); // not really needed for

    zStiff_ = GenerateStdMatrixObject( eType, sType,
                                       nrows, ncols, fill );
    zStiff_->SetSparsityPattern( *matrixA_ );
    zStiff_->Add( 1.0, *matrixA_ );

    sType = matrixB_->GetStorageType();
    nrows = matrixB_->GetNumRows();
    ncols = matrixB_->GetNumCols();
    fill = matrixB_->GetNnz();
    zMass_ = GenerateStdMatrixObject( eType, sType,
                                      nrows, ncols, fill );
    zMass_->SetSparsityPattern( *matrixB_ );
    zMass_->Add( 1.0, *matrixB_ );

    sType = matrixD_->GetStorageType();
    nrows = matrixD_->GetNumRows();
    ncols = matrixD_->GetNumCols();
    fill = matrixD_->GetNnz();
    zDamp_ = GenerateStdMatrixObject( eType, sType,
                                      nrows, ncols, fill );
    zDamp_->SetSparsityPattern( *matrixD_ );
    zDamp_->Add( 1.0, *matrixD_ );

    // Create matrix interface for arpack
    interface_ = new ArpackMatInterface( zStiff_, zMass_, zDamp_, shiftAndInvert_, freqShift_ );

    // Create solver class
    arpackSolver_ = new ArpackSolver(xml_);
    
    // Check 'which'-settings regarding the type of eigenvalues searched for
    std::string whichString = "SI";
    xml_->GetValue("which", whichString, ParamNode::INSERT );
    which_ = new char[whichString.size()+1];
    strncpy(which_, whichString.c_str(), whichString.size()+1 );

    // Create the solver dependent on the problem type ( regular / shiftAndInvert )
    arpackSolver_->QuadSetup( interface_, size, numFreq_, freqShift_, which_ , shiftAndInvert_);

    // Set additional scaling parameter for matrix interface
    Double scaleFac = 1.0;
    xml_->GetValue("scaleFac", scaleFac, ParamNode::INSERT);
    if (scaleFac > 0.0)
      interface_->SetDiagScaling(scaleFac);

    if (scaleFac != 1.0)
      std::cout << "       applied scaling factor      :  " << interface_->GetDiagScaling() << "\n";

    std::cout << std::endl;

    // Create solver
    solver_ = GenerateSolverObject(*zStiff_, solStrat_, solverList_, info_);


    // Perform check, if matrix is std or sbm
    if ( matrixA_->GetStructureType() == BaseMatrix::SPARSE_MATRIX ) {
      const StdMatrix & mat = dynamic_cast< const StdMatrix &>( *matrixB_ );

      // Create preconditioner
      precond_ = GenerateStdPrecondObject( mat, solStrat_->GetPrecondId(), precondList_, info_ );
    } else {
      EXCEPTION( "No preconditioner available for SBM-matrices!" );
    }

    // Setup matrixinterface
    interface_->QuadSetup( solver_, precond_ );

    ToInfo();

  }

  void ArpackEigenSolver::SetupIndex(unsigned int numev)
  {
    // to be called after within each CalcEigenFrequencies()
    idx_.Resize(numev);

    if(sort_)
    {
      // we sort with a std::pair<||ev||, org_idx>
      std::vector<ev_idx> org;
      for(unsigned int i = 0; i < numev; i++)
      {
        double v = isBloch_ || isQuadratic_ ? std::abs(arpackSolver_->CmplxEigenvalue(i)) : arpackSolver_->Eigenvalue(i);
        org.push_back(std::make_pair(v, i));
      }

      std::sort(org.begin(), org.end(), ArpackEigenSolver::comperator);

      for(unsigned int i = 0; i < numev; i++)
      {
        LOG_DBG2(aes) << "SI: i=" << i << " (" << org[i].first << ", " << org[i].second << ")";
        idx_[i] = org[i].second;
      }
    }
    else
      for(unsigned int i = 0; i < numev; i++)
        idx_[i] = i;
  }

  void ArpackEigenSolver::CalcEigenFrequencies(BaseVector &sol, BaseVector &err)
  {
    shared_ptr<Timer> timer = info_->Get("arpack/timer")->AsTimer();
    timer->Start();

    assert(!(isBloch_ && isQuadratic_));

    unsigned int numEVs = 0;
    // case1: generalized real problem
    if(!isQuadratic_ && !isBloch_)
    {
      // Find the eigenvalues and calculate the eigenvectors
      numEVs = arpackSolver_->FindEigenvalues<Double>();

      SetupIndex(numEVs); // setup idx_

      // case1: generalized problem
      Vector<Double> & solConverted = dynamic_cast<Vector<Double>&>(sol);
      solConverted.Resize( numEVs );
      for (UInt i = 0; i < numEVs; i++ )
      {
        solConverted[i] = arpackSolver_->Eigenvalue(idx_[i]);
        // if non-negative eigenvalue, convert to eigenfrequency
        if (solConverted[i] >= 0.0)  {
          LOG_DBG(aes) << "CEF: i=" << i << " ev=" << solConverted[i] << " -> f=" << sqrt(solConverted[i])/(8.0*atan(1.0));
          solConverted[i] = sqrt(solConverted[i])/(8.0*atan(1.0));
        }
      }
    }
    // case2: quadratic complex problem
    // case3: bloch modes are generalized complex problems
    if(isQuadratic_ || isBloch_)
    {
      if(isQuadratic_)
        numEVs = arpackSolver_->FindQuadEigenvalues();
      else
        numEVs = arpackSolver_->FindEigenvalues<Complex>();

      SetupIndex(numEVs);

      Vector<Complex> & solConverted = dynamic_cast<Vector<Complex>&>(sol);
      solConverted.Resize( numEVs );
      for (UInt i = 0; i < numEVs; i++ ) {
        solConverted[i] = arpackSolver_->CmplxEigenvalue(idx_[i]);
        if(isBloch_ && solConverted[i].real() >= 0.0)
          solConverted[i] = sqrt(solConverted[i])/(8.0*atan(1.0));
      }
    }

    PtrParamNode in = this->info_->Get("arpack", progOpts->DoDetailedInfo() ? ParamNode::APPEND : ParamNode::INSERT);
    in->Get("analysis_id")->SetValue(domain->GetDriver()->GetAnalysisId().ToString());
    in->Get("rci/solve_x")->SetValue(arpackSolver_->counter_solve_OP_x);
    in->Get("rci/solve_B_x")->SetValue(arpackSolver_->counter_solve_OP_B_x);
    in->Get("rci/matvec_B_x")->SetValue(arpackSolver_->counter_B_x);
    in->Get("rci/total")->SetValue(arpackSolver_->counter_calll_aupd);

    // Save error norms
    Vector<Double> & errVec = dynamic_cast<Vector<Double>&>(err);
    errVec.Resize( numEVs );
    for (UInt i = 0; i < numEVs; i++ ) {
        errVec[i] = arpackSolver_->Tolerance(idx_[i]);
    }

    timer->Stop();
  }

  void ArpackEigenSolver::CalcConditionNumber(const BaseMatrix& mat, Double& condNumber, Vector<Double>& evs, Vector<Double>& err )
  {
    // Set flag for indicating a non-quadratic problem
    isQuadratic_ = false;
    isGeneralized_ = false;
    // NOTE: Hard coded as true!!!
    shiftAndInvert_ = true;

    // Copy matrix references and determine size of system
    matrixA_ = & dynamic_cast<const StdMatrix&>(mat);

    UInt size = matrixA_->GetNumRows();

    // Create matrix interface for arpack
    interface_ = new ArpackMatInterface( matrixA_, shiftAndInvert_, freqShift_ );

    // Create solver class
    arpackSolver_ = new ArpackSolver(xml_);

    // set mode: we look at both ends of the spectrum for eigenvalues
    std::string whichString = "BE";
    which_ = new char[whichString.size()+1];
    strncpy(which_, whichString.c_str(), whichString.size()+1 );

    // Save frequency parameters
    numFreq_ = 10; 
    freqShift_ = 0.0;

    // Create the solver dependent on the problem type ( regular / shiftAndInvert )
    arpackSolver_->Setup( interface_, size, numFreq_, freqShift_, which_ , (char*) "I", shiftAndInvert_, false);

    // Create standard solver
    solver_ = GenerateSolverObject(*matrixA_, solStrat_, solverList_, info_);

    // Perform check, if matrix is std or sbm
    if ( matrixA_->GetStructureType() == BaseMatrix::SPARSE_MATRIX ) {
      const StdMatrix & mat = dynamic_cast< const StdMatrix &>( *matrixA_ );

      // Create preconditioner
      precond_ = GenerateStdPrecondObject( mat,  solStrat_->GetPrecondId(), precondList_, info_ );
    } else {
      EXCEPTION( "No preconditioner available for SBM-matrices!" );
    }

    // Setup matrixinterface
    interface_->Setup( solver_, precond_ );

    //Find the eigenvalues and calculate the eigenvectors
    UInt numEVs = arpackSolver_->FindEigenvalues<Double>();
    SetupIndex(numEVs);

    if( numEVs == 0) {
      EXCEPTION( "No convergence in search for eigenvalues");
    }

    evs.Resize( numEVs );
    err.Resize( numEVs );
    for (UInt i = 0; i < numEVs; i++ ) {
      evs[i] = arpackSolver_->Eigenvalue(idx_[i]);
      err[i] = arpackSolver_->Tolerance(idx_[i]);
    }

    condNumber = evs[evs.GetSize()-1] / evs[0];
  }

  void ArpackEigenSolver::GetEigenMode(UInt modeNr, Vector<Complex> & mode, bool right)
  {
    UInt size = matrixA_->GetNumRows();
    mode.Resize( size );
    mode.Init();

    for(UInt i = 0; i < size; i++)
      mode[i] = Complex((arpackSolver_->GetEigenvector(idx_[modeNr]))[i],0);

    // BLOCH better? mode.Fill(arpackSolver_->GetEigenvector(modeNr), matrixA_->GetNumRows());

  }

  void ArpackEigenSolver::GetComplexEigenMode(UInt modeNr, Vector<Complex>& mode)
  {
    // in bloch mode case the same as GetEigenMode,
    // in quadratic case the modes have internally double size and we want the upper half
    
    UInt size = matrixA_->GetNumRows();

    if(isQuadratic_)
      mode.Fill(arpackSolver_->GetComplexEigenvector(idx_[modeNr]) + size, size);
    else
      mode.Fill(arpackSolver_->GetComplexEigenvector(idx_[modeNr]), size);
  }
}
