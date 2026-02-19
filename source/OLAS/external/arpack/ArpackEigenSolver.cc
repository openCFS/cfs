#include <limits>
#include <string.h>
#include <boost/type_traits/is_complex.hpp>

#include "MatVec/StdMatrix.hh"
#include "MatVec/generatematvec.hh"
#include "Domain/Domain.hh"
#include "Driver/BaseDriver.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "MatVec/SparseOLASMatrix.hh"
#include "OLAS/algsys/SolStrategy.hh"
#include "OLAS/precond/generateprecond.hh"
#include "OLAS/precond/BasePrecond.hh"
#include "OLAS/solver/generatesolver.hh"
#include "OLAS/solver/BaseSolver.hh"
#include "Utils/Timer.hh"

#include "ArpackEigenSolver.hh"

DEFINE_LOG(aes, "arpackEigenSolver")


namespace CoupledField {

  ArpackEigenSolver::ArpackEigenSolver( shared_ptr<SolStrategy> strat,
                                        PtrParamNode xml,
                                        PtrParamNode solverList,
                                        PtrParamNode precondList,
                                        PtrParamNode eigenInfo)
    : BaseEigenSolver( strat, xml, solverList, precondList, eigenInfo )
  {
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
    
    eigenSolverName_ = EigenSolverType::ARPACK;

    isGeneralized_ = false;
    computeMode_ = ArpackMatInterface::ComputeMode::REGULAR;

    scale_B_ = xml_->Has("scale_B") ? xml_->Get("scale_B")->As<bool>() : false;

    BaseEigenSolver::PostInit();
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

  void ArpackEigenSolver::Setup(const BaseMatrix & A, bool isHermitian) {
    // Set flag for indicating a non-quadratic problem
    isGeneralized_ = false;
    sort_ = false;

    this->SetProblemType(A,isHermitian);

    // NOTE: Hard coded!!!
    computeMode_ = ArpackMatInterface::ComputeMode::SHIFT_INVERT;

    // Copy matrix references and determine size of system
    matrixA_ = & dynamic_cast<const StdMatrix&>(A);

    UInt size = matrixA_->GetNumRows();

    // Create matrix interface for arpack
    interface_ = new ArpackMatInterface( matrixA_, computeMode_ );

    // Create solver class
    arpackSolver_ = new ArpackSolver(xml_);

    // Check 'which'-settings regarding the type of eigenvalues searched for
    std::string whichString = "LM";
    xml_->GetValue("which", whichString, ParamNode::INSERT);

    which_ = new char[whichString.size()+1];
    strncpy(which_, whichString.c_str(), whichString.size()+1 );

    // Create the solver dependent on the problem type ( regular | shiftAndInvert | buckling)
    arpackSolver_->Setup(interface_, size, which_ , (char*) "I", computeMode_, isBloch_||this->matrixA_->GetEntryType() == BaseMatrix::COMPLEX);

    // Create solver
    solver_ = GenerateSolverObject( *matrixA_, solStrat_, solverList_, info_ );

    // Perform check, if matrix is std or sbm
    if ( matrixA_->GetStructureType() == BaseMatrix::SPARSE_MATRIX ) {
      //const StdMatrix & mat = dynamic_cast< const StdMatrix &>( *matrixA_ );

      // Create preconditioner
      precond_ = GenerateStdPrecondObject( *matrixA_, solStrat_->GetPrecondId(), precondList_, info_ );
    } else {
      EXCEPTION( "No preconditioner available for SBM-matrices!" );
    }

    ToInfo();
  }

  void ArpackEigenSolver::Setup(const BaseMatrix & A, const BaseMatrix & B, bool isHermitian) {
    Setup(A, B, ArpackMatInterface::ComputeMode::SHIFT_INVERT, isHermitian);
  }

  void ArpackEigenSolver::Setup(const BaseMatrix & A, const BaseMatrix & B, ArpackMatInterface::ComputeMode computeMode, bool isHermitian) {
    isQuadratic_ = false;
    isGeneralized_ = true;
    sort_ = false;

    this->SetProblemType(A,isHermitian);

    computeMode_ = computeMode;
    assert(computeMode_ == ArpackMatInterface::ComputeMode::SHIFT_INVERT || computeMode_ == ArpackMatInterface::ComputeMode::BUCKLING);

    if(scale_B_)
    {
      // we scale the B-Matrix as suggested by Jonas:
      // A*x=lambda*sigma*B*x, with z.B. sigma=1.0/B(1,1) and then rescale
      if(B.GetEntryType() == BaseMatrix::COMPLEX)
      {
        const SparseOLASMatrix<Complex>* olas = dynamic_cast<const SparseOLASMatrix<Complex>* >(&B);
        assert(olas != NULL);
        Complex b11 = olas->GetAvgDiag();
        assert(b11.real() > 1e-15);
        scale_B_val_ = scale_B_ ? (1./b11).real() : 1.0;
      }
      else
      {
        const SparseOLASMatrix<Double>* olas = dynamic_cast<const SparseOLASMatrix<Double>* >(&B);
        assert(olas != NULL);
        Double b11 = olas->GetAvgDiag();
        assert(b11 > 1e-15);
        scale_B_val_ = scale_B_ ? 1./b11 : 1.0;
      }
    }

    LOG_DBG(aes) << "AES:S B-scale=" << scale_B_val_;

    // Copy matrix references and determine size of system
    matrixA_ = & dynamic_cast<const StdMatrix&>(A);
    if(scale_B_)
    {
      StdMatrix* Temp = CopyStdMatrixObject(dynamic_cast<const StdMatrix&>(B));
      Temp->Scale(scale_B_val_);
      matrixB_ = Temp;
      Temp = NULL;
    }
    else
      matrixB_ = & dynamic_cast<const StdMatrix&>(B);

    UInt size = matrixA_->GetNumRows();

    // Create matrix interface for arpack, it might already exist in the bloch mode case
    if(interface_ == NULL)
      interface_ = new ArpackMatInterface( matrixA_, matrixB_, computeMode_ );

    LOG_DBG(aes) << "computational mode: " << computeMode_;

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

    // Create the solver dependent on the problem type ( regular | shiftAndInvert | buckling)
    arpackSolver_->Setup( interface_, size, which_ , (char*) "G", computeMode_, isBloch_||this->matrixA_->GetEntryType() == BaseMatrix::COMPLEX);

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

    ToInfo();
  }

  void ArpackEigenSolver::Setup(const BaseMatrix & mat, UInt numFreq, Double freqShift, bool sort)
  {
    // Set flag for indicating a non-quadratic problem
    isQuadratic_ = false;
    isGeneralized_ = false;
    isBloch_ = false;
    sort_ = false;

    this->SetProblemType(mat);

    // NOTE: Hard coded!!!
    computeMode_ = ArpackMatInterface::ComputeMode::SHIFT_INVERT;

    // Copy matrix references and determine size of system
    matrixA_ = & dynamic_cast<const StdMatrix&>(mat);

    UInt size = matrixA_->GetNumRows();

    // Save frequency parameters
    numFreq_ = numFreq;
    freqShift_ = freqShift;

    // Create matrix interface for arpack
    interface_ = new ArpackMatInterface( matrixA_, computeMode_);

    // Create solver class
    arpackSolver_ = new ArpackSolver(xml_);

    // Check 'which'-settings regarding the type of eigenvalues searched for
    std::string whichString = "LM";
    xml_->GetValue("which", whichString, ParamNode::INSERT);
    
    which_ = new char[whichString.size()+1];
    strncpy(which_, whichString.c_str(), whichString.size()+1 );

    // Create the solver dependent on the problem type ( regular | shiftAndInvert | buckling)
    arpackSolver_->Setup(interface_, size, which_ , (char*) "I", computeMode_, isBloch_ || this->matrixA_->GetEntryType() == BaseMatrix::COMPLEX);

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
    interface_->Setup( solver_, precond_, pow(freqShift*2.0*M_PI,2), GetSetupTimer(), GetSolveTimer());

    ToInfo();
  }
  
  
  void ArpackEigenSolver::Setup(const BaseMatrix& stiffMat, const BaseMatrix& massMat, UInt numFreq, Double freqShift, bool sort, bool bloch)
  {
    // Set flag for indicating a non-quadratic problem
    isQuadratic_ = false;
    isGeneralized_ = true;
    isBloch_ = bloch;
    sort_ = sort;

    this->SetProblemType(stiffMat);

    if(scale_B_)
    {
      // we scale the B-Matrix as suggested by Jonas:
      // A*x=lambda*sigma*B*x, with z.B. sigma=1.0/B(1,1) and then rescale
      if(massMat.GetEntryType() == BaseMatrix::COMPLEX)
      {
        const SparseOLASMatrix<Complex>* olas = dynamic_cast<const SparseOLASMatrix<Complex>* >(&massMat);
        assert(olas != NULL);
        Complex b11 = olas->GetAvgDiag();
        assert(b11.real() > 1e-15);
        scale_B_val_ = scale_B_ ? (1./b11).real() : 1.0;
      }
      else
      {
        const SparseOLASMatrix<Double>* olas = dynamic_cast<const SparseOLASMatrix<Double>* >(&massMat);
        assert(olas != NULL);
        Double b11 = olas->GetAvgDiag();
        assert(b11 > 1e-15);
        scale_B_val_ = scale_B_ ? 1./b11 : 1.0;
      }
    }

    LOG_DBG(aes) << "Setup: B-scale=" << scale_B_val_;

    // bloch works only for non-symmetric matrices as the stiffness matrix needs to be Hermitian.
    // At least Pardiso can be used with <pardiso> <hermitean>yes</hermitean> </pardiso>

    // Copy matrix references and determine size of system
    matrixA_ = & dynamic_cast<const StdMatrix&>(stiffMat);
    if(scale_B_)
    {
      StdMatrix* Temp = CopyStdMatrixObject(dynamic_cast<const StdMatrix&>(massMat));
      Temp->Scale(scale_B_val_);
      matrixB_ = Temp;
    }
    else
      matrixB_ = & dynamic_cast<const StdMatrix&>(massMat);

    UInt size = matrixA_->GetNumRows();

    // Save frequency parameters
    numFreq_ = numFreq;
    freqShift_ = freqShift;

    if (computeMode_ == ArpackMatInterface::ComputeMode::REGULAR) {
      computeMode_ = ArpackMatInterface::ComputeMode::SHIFT_INVERT;
    }

    // Create matrix interface for arpack, it might already exist in the bloch mode case
    if(interface_ == NULL)
      interface_ = new ArpackMatInterface( matrixA_, matrixB_, computeMode_ );

    LOG_DBG(aes) << "computational mode: " << computeMode_;
    LOG_DBG(aes) << "shift-and-invert mode: S(S,M) shift=" << freqShift;

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

    // Create the solver dependent on the problem type ( regular | shiftAndInvert | buckling)
    arpackSolver_->Setup( interface_, size, which_ , (char*) "G", computeMode_, isBloch_ || this->matrixA_->GetEntryType() == BaseMatrix::COMPLEX);

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
    interface_->Setup( solver_, precond_, pow(freqShift*2.0*M_PI,2), GetSetupTimer(), GetSolveTimer());

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


  void ArpackEigenSolver::Setup(const BaseMatrix & A, const BaseMatrix & B, const BaseMatrix & C)
  {
    // Set flag for indicating a quadratic problem
    isQuadratic_ = true;
    isGeneralized_ = true;
    isBloch_ = false;
    sort_ = false;

    // Copy matrix references and convert them to StdMatrices
    matrixA_ = & dynamic_cast<const StdMatrix&>(A);

    matrixB_ = & dynamic_cast<const StdMatrix&>(B);

    matrixD_ = & dynamic_cast<const StdMatrix&>(C);

    // need twice the original system size and twice the number of frequencies
    UInt size = 2*matrixA_->GetNumRows();
    // NOTE: Hard coded!!!
    computeMode_ = ArpackMatInterface::ComputeMode::SHIFT_INVERT;

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
    interface_ = new ArpackMatInterface( zStiff_, zMass_, zDamp_, computeMode_ );
    eigenProblemType_ = COMPLEX_SYMMETRIC;
    // Create solver class
    arpackSolver_ = new ArpackSolver(xml_);
    // Check 'which'-settings regarding the type of eigenvalues searched for
    std::string whichString = "SI";
    xml_->GetValue("which", whichString, ParamNode::INSERT );
    which_ = new char[whichString.size()+1];
    strncpy(which_, whichString.c_str(), whichString.size()+1 );

    // Create the solver dependent on the problem type ( regular | shiftAndInvert | buckling)
    arpackSolver_->QuadSetup( interface_, size, which_ , computeMode_);

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

    ToInfo();
  }

  void ArpackEigenSolver::Setup(const BaseMatrix & stiffMat, const BaseMatrix & massMat, const BaseMatrix & dampMat,
                                UInt numFreq, double freqShift, bool sort) {

    // Set flag for indicating a non-quadratic problem
    isQuadratic_ = true;
    isGeneralized_ = true;
    isBloch_ = false;
    sort_ = false;

    this->SetProblemType(stiffMat);

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

    // NOTE: Hard coded!!!
    computeMode_ = ArpackMatInterface::ComputeMode::SHIFT_INVERT;

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
    interface_ = new ArpackMatInterface( zStiff_, zMass_, zDamp_, computeMode_ );

    // Create solver class
    arpackSolver_ = new ArpackSolver(xml_);

    // Check 'which'-settings regarding the type of eigenvalues searched for
    std::string whichString = "SI";
    xml_->GetValue("which", whichString, ParamNode::INSERT );
    which_ = new char[whichString.size()+1];
    strncpy(which_, whichString.c_str(), whichString.size()+1 );

    // Create the solver dependent on the problem type ( regular / shiftAndInvert )
    arpackSolver_->QuadSetup( interface_, size, which_ , computeMode_);

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
    interface_->QuadSetup( solver_, precond_, freqShift );

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
    if (isQuadratic_)
      CalcEigenValues(sol,err,numFreq_,freqShift_*2.0*M_PI);
    else
      CalcEigenValues(sol,err,numFreq_,pow(freqShift_*2.0*M_PI,2));

    unsigned int numEVs = sol.GetSize();

    // case1: generalized real problem
    if(!isQuadratic_ && !isBloch_)
    {
      Vector<Double> & solConverted = dynamic_cast<Vector<Double>&>(sol);
      for (UInt i = 0; i < numEVs; i++ ) {
        // if non-negative eigenvalue, convert to eigenfrequency
        if (solConverted[i] >= 0.0)  {
          LOG_DBG(aes) << "CEF: i=" << i << " ev=" << solConverted[i] << " -> f=" << sqrt(solConverted[i])/(2.0*M_PI);
          solConverted[i] = sqrt(solConverted[i])/(2.0*M_PI);
        }
      }
    }
    // case2: quadratic complex problem
    // case3: bloch modes are generalized complex problems
    if(isQuadratic_ || isBloch_)
    {
      Vector<Complex> & solConverted = dynamic_cast<Vector<Complex>&>(sol);
      for (UInt i = 0; i < numEVs; i++ ) {
        // if real part of eigenvalue non-negative, convert to eigenfrequency
        if(isBloch_ && solConverted[i].real() >= 0.0) {
          LOG_DBG(aes) << "CEF: i=" << i << " ev=" << solConverted[i] << " -> f=" << sqrt(solConverted[i])/(2.0*M_PI);
          solConverted[i] = sqrt(solConverted[i])/(2.0*M_PI);
        }
      }
    }
  }

  void ArpackEigenSolver::CalcEigenValues(BaseVector& sol, BaseVector& err, UInt N, Double shiftPoint)
  {
    unsigned int numEVs = 0;
    // generalized real problem
    if( !(isQuadratic_ || isBloch_ || this->matrixA_->GetEntryType() == BaseMatrix::COMPLEX) )
    {
      // Setup matrixinterface
      interface_->Setup( solver_, precond_, shiftPoint, GetSetupTimer(), GetSolveTimer());

      // Find the eigenvalues and calculate the eigenvectors
      numEVs = arpackSolver_->FindEigenvalues<Double>(N, shiftPoint);

      SetupIndex(numEVs);
      Vector<Double> & solConverted = dynamic_cast<Vector<Double>&>(sol);
      solConverted.Resize( numEVs );
      for (UInt i = 0; i < numEVs; i++ )
      {
        solConverted[i] = arpackSolver_->Eigenvalue(idx_[i]);
        LOG_DBG(aes) << "CEV: i=" << i << " ev=" << solConverted[i];
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
    }
    else
      CalcEigenValues(sol, err, N, Complex(shiftPoint, 0.0));
  }

  void ArpackEigenSolver::CalcEigenValues(BaseVector& sol, BaseVector& err, UInt N, Complex shiftPoint)
  {
    assert(!(isBloch_ && isQuadratic_));
    assert(isQuadratic_ || isBloch_ || this->matrixA_->GetEntryType() == BaseMatrix::COMPLEX);

    unsigned int numEVs = 0;
    // complex problem or
    // quadratic complex problem or
    // bloch modes (generalized complex problem)
    if(isQuadratic_) {
      interface_->QuadSetup<Complex>( solver_, precond_, shiftPoint);
      numEVs = arpackSolver_->FindQuadEigenvalues<Complex>(N, shiftPoint);
    } else {
      interface_->Setup<Complex>( solver_, precond_, shiftPoint, GetSetupTimer(), GetSolveTimer());
      numEVs = arpackSolver_->FindEigenvalues<Complex>(N, shiftPoint);
    }

    SetupIndex(numEVs);

    Vector<Complex> & solConverted = dynamic_cast<Vector<Complex>&>(sol);
    solConverted.Resize( numEVs );
    for (UInt i = 0; i < numEVs; i++ ) {
      solConverted[i] = arpackSolver_->CmplxEigenvalue(idx_[i]);
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
  }


  void ArpackEigenSolver::CalcConditionNumber(const BaseMatrix& mat, Double& condNumber, Vector<Double>& evs, Vector<Double>& err )
  {
    // Set flag for indicating a non-quadratic problem
    isQuadratic_ = false;
    isGeneralized_ = false;
    // NOTE: Hard coded!!!
    computeMode_ = ArpackMatInterface::ComputeMode::SHIFT_INVERT;

    // Copy matrix references and determine size of system
    matrixA_ = & dynamic_cast<const StdMatrix&>(mat);

    UInt size = matrixA_->GetNumRows();

    // Create matrix interface for arpack
    interface_ = new ArpackMatInterface( matrixA_, computeMode_ );

    // Create solver class
    arpackSolver_ = new ArpackSolver(xml_);

    // set mode: we look at both ends of the spectrum for eigenvalues
    std::string whichString = "BE";
    which_ = new char[whichString.size()+1];
    strncpy(which_, whichString.c_str(), whichString.size()+1 );

    // Save frequency parameters
    numFreq_ = 10; 
    freqShift_ = 0.0;

    // Create the solver dependent on the problem type ( regular | shiftAndInvert | buckling)
    setupTimer_->Start();
    arpackSolver_->Setup( interface_, size, which_ , (char*) "I", computeMode_, false);

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
    interface_->Setup( solver_, precond_, freqShift_, GetSetupTimer(), GetSolveTimer());
    setupTimer_->Stop();

    //Find the eigenvalues and calculate the eigenvectors
    solveTimer_->Start();
    UInt numEVs = arpackSolver_->FindEigenvalues<Double>(numFreq_, freqShift_);
    solveTimer_->Stop();
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

    // we rescale x by sqrt(scale_B_val_) as we have x^T (scale*B) x = 1
    double rescale = std::sqrt(scale_B_val_);

    if(this->matrixA_->GetEntryType() == BaseMatrix::COMPLEX || isQuadratic_)
    {
    	this->GetComplexEigenMode(modeNr, mode);
    	mode.ScalarMult(rescale);
    }
    else
    {
      for(UInt i = 0; i < size; i++)
        mode[i] = Complex((arpackSolver_->GetEigenvector(idx_[modeNr]))[i]*rescale,0);
    }

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
