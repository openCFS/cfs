#include "EigenFrequencyDriver.hh"

#include <iostream>
#include <iomanip>
#include <cmath>
#include <filesystem>
#include "Utils/mathParser/mathParser.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/SimState.hh"
#include "DataInOut/ResultHandler.hh"
#include "DataInOut/ProgramOptions.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "Driver/SolveSteps/StdSolveStep.hh"
#include "Domain/Domain.hh"
#include "MatVec/SBM_Matrix.hh"
#include "OLAS/solver/BaseEigenSolver.hh"
#include "Optimization/Optimization.hh"
#include "PDE/StdPDE.hh"

using std::cout;
using std::setw;
using std::string;

namespace CoupledField {

  DEFINE_LOG(efd, "eigenFrequencyDriver")

  // forward declaration
  class BaseVector;

  // ***************
  //   Constructor
  // ***************
  EigenFrequencyDriver::EigenFrequencyDriver(UInt sequenceStep,
                                             bool isPartOfSequence,
                                             shared_ptr<SimState> state,
                                             Domain* domain,
                                             PtrParamNode paramNode, PtrParamNode infoNode ) 
    : SingleDriver( sequenceStep, isPartOfSequence, state, domain, paramNode, infoNode)
  {
    // set correct analysistype
    analysis_ = BasePDE::EIGENFREQUENCY;

    numFreq_ = 0;
    freqShift_ = 0.0;
    writeModes_ = true;
    isQuadratic_ = false;
    blochSteps_ = 0;
    ibz_     = false;
    eigenFreqs = NULL;
    minVal_ = 0.0;
    maxVal_ = 0.0;
    eigenValuesAreReal_=true;

    // replace with a concrete element
    param_ = param_->Get("eigenFrequency");
    info_ = info_->Get("eigenFrequency");

    isBloch_ = DoBloch(param_);
    sort_ = param_->Get("sort")->As<bool>();

    boundary.SetName("EigenFrequencyDriver::Boundary");
    boundary.Add(SYMMETRIC, "symmetric");
    boundary.Add(QUADRANT, "quadrant");
    boundary.Add(HORIZONZAL, "horizontal");
    boundary.Add(FULL, "full");
    boundary_ = NO_IBZ;
  }

  EigenFrequencyDriver::~EigenFrequencyDriver()
  {
    if(isBloch_)
    {
      bloch_plot_.close();
      if(domain->GetOptimization())
      {
        if(!progOpts->DoDetailedInfo())
        {
          string file = progOpts->GetSimName() + ".bloch.dat"; // See SetupBlochPlot()
          std::rename(string(file + ".tmp").c_str(), file.c_str());
        }
        else
          std::filesystem::copy_file(bloch_name_, progOpts->GetSimName() + ".bloch.dat", std::filesystem::copy_options::overwrite_existing);
      }
    }

    //delete eigenFreqs;
    eigenFreqs = NULL;
  }

  void EigenFrequencyDriver::Init(bool restart)
  {
    // read required parameters from parameter node
    param_->GetValue( "numModes", numFreq_ , ParamNode::INSERT );
    param_->GetValue( "freqShift", freqShift_ , ParamNode::INSERT ); // this should be shift-point and is in general complex valued
    param_->GetValue( "minVal", minVal_, ParamNode::INSERT );
    param_->GetValue( "maxVal", maxVal_, ParamNode::INSERT );
    param_->GetValue( "writeModes", writeModes_, ParamNode::PASS );
    param_->GetValue( "isQuadratic", isQuadratic_, ParamNode::PASS );
    // determine type of mode normalization
    std::string normString = "solver";
    param_->Get("writeModes")->GetValue("normalization",normString, ParamNode::PASS);
    if (normString == "solver") {
      modeNormalization_ = BaseEigenSolver::NONE;
    } else if (normString == "max") {
      modeNormalization_ = BaseEigenSolver::MAX;
    } else if (normString == "norm") {
      modeNormalization_ = BaseEigenSolver::NORM;
    } else {
      EXCEPTION("Specified mode normalization '"+normString+"' not implemented");
    }
    // read flag if all results should get written to database file section
    // to allow e.g. for general postprocessing or result extraction    
    param_->GetValue("allowPostProc", writeAllSteps_, ParamNode::PASS );
    isBloch_ = param_->Has("bloch");
    if(isBloch_)
    {
      FillWaveVectors(param_->Get("bloch"));
      SetupBlochPlot();
    }

    if(isQuadratic_ && isBloch_)
      throw Exception("Bloch mode Eigenfrequency analysis not implemented for quadratic form");

    // the eigenfrequencies are complex in the quadratic case or in bloch mode
    //if(isQuadratic_ || isBloch_) eigenFreqs = new Vector<Complex>(numFreq_);
    //                        else eigenFreqs = new Vector<Double>(numFreq_);

    eigenFreqs = & frequency_;
    InitializePDEs();
  }

  void EigenFrequencyDriver::SetToStepValue(UInt stepNum, Double stepVal )
  {
    // ensure that this method is only called if simState has input
    if( false ) { // ! simState_->HasInput()
      EXCEPTION( "Can only set external time step, if simulation state "
              << "is read from external file" );
    }

    //actFreqStep_ = stepNum;;
    //actFreq_ = stepVal;

    // Set current frequency value in the mathParser
    domain_->GetMathParser()->SetValue( MathParser_GLOB_HANDLER, "f", stepVal );
    domain_->GetMathParser()->SetValue( MathParser_GLOB_HANDLER, "step", stepNum );
  }


  void EigenFrequencyDriver::SetupBlochPlot()
  {
    bloch_plot_.close();

    // non detailed run:
    //    we write in .tmp the ongoing simulation and rename it when we are finished such that we have the last valid full bloch data for
    //    analysis of running optimization via find_band_gap.py
    // detailed run:
    //    the bloch dat file has the iteration number, but we also copy the last one to .bloch.dat
    if(domain->GetOptimization() && progOpts->DoDetailedInfo() && bloch_name_ != "")
      std::filesystem::copy_file(bloch_name_, progOpts->GetSimName() + ".bloch.dat", std::filesystem::copy_options::overwrite_existing);

    bloch_name_ = progOpts->GetSimName(); // on change check destructor!
    if(domain->GetOptimization() && progOpts->DoDetailedInfo())
      bloch_name_ += ".iter_" + boost::lexical_cast<string>(domain->GetOptimization()->GetCurrentIteration());
    bloch_name_ += ".bloch.dat";

    if(domain->GetOptimization() && !progOpts->DoDetailedInfo())
    {
      // dosn't exist in the first wave_vector of iteration 0
      std::rename(string(bloch_name_ + ".tmp").c_str(), bloch_name_.c_str());
      bloch_name_ += ".tmp";
    }

    int dim = domain->GetGrid()->GetDim();
    int edges = boundary_ == HORIZONZAL ? 1 : boundary_ == SYMMETRIC ? dim + 1 : dim + 2; // copy & paste with FillWaveVector()
    bloch_plot_.open(bloch_name_.c_str() , std::ios::out);
    bloch_plot_ << "# <ibz dim=\"" << dim << "\" edges=\"" << edges << "\"/>\n";
    bloch_plot_ << "#step\tk_x\tk_y";
    if(domain->GetGrid()->GetDim() == 3)
      bloch_plot_ << "\tk_z";
    bloch_plot_ << "\t1.mode\t2.mode\t..." << std::endl;
  }

  unsigned int EigenFrequencyDriver::GetNumModes(PtrParamNode node)
  {
    assert(node->Has("numModes"));
    return node->Get("numModes")->As<unsigned int>();
  }


  unsigned int EigenFrequencyDriver::GetNumBlochWave(PtrParamNode node)
  {
    if(!node->Has("bloch"))
      return 0;

    if(node->Has("bloch/ibz/steps"))
      return node->Get("bloch/ibz/steps")->As<unsigned int>();

    return node->Get("bloch")->GetList("waveVector").GetSize();
  }

  void EigenFrequencyDriver::FillWaveVectors(PtrParamNode bloch_pn)
  {
    ParamNodeList lst = bloch_pn->GetList("waveVector");

    unsigned int dim = domain->GetGrid()->GetDim();

    // check if we have ibz
    if(bloch_pn->Has("ibz"))
    {
      ibz_ = true;
      if(!lst.IsEmpty())
        throw Exception("no wave vectors may be given in bloch mode analysis concurrently with ibz");

      boundary_ = boundary.Parse(bloch_pn->Get("ibz/sample")->As<string>());
      assert(boundary_ != NO_IBZ); // no option in schema

      // copy and paste in SetupBlochPlot()
      int edges = boundary_ == HORIZONZAL ? 1 : boundary_ == SYMMETRIC ? dim + 1 : dim + 2; // GAMMA -> X -> M -> [GAMMA | R -> GAMMA]

      blochSteps_ = bloch_pn->Get("ibz/steps")->As<int>();
      int steps = blochSteps_; // shortcut only

      if((boundary_ != FULL && ((steps % edges) != 0 || steps < edges)) || (boundary_ == FULL && (int) sqrt(steps) * (int) sqrt(steps) != steps))
        EXCEPTION("bloch mode ibz/steps need to a multiple of " << edges << " for boundary sampling or a square number for full sampling.");

      // we need the unit cell dimensions to scale the wave vector from 0 .. pi/d where d is unit cell dimension -> Hussein;2009
      Matrix<double>& box = domain->GetGrid()->CalcGridBoundingBox(NULL, true); // force 3D
      double d_x = box[0][1]-box[0][0];
      double d_y = box[1][1]-box[1][0];
      double d_z = box[2][1]-box[2][0]; // 0.0 for 2D
      wave_vectors.Resize(steps);

      LOG_DBG(efd) << "FWV #wv=" << wave_vectors.GetSize();
      LOG_DBG(efd) << "FWV box: " << box.ToString();
      LOG_DBG(efd) << "FWV box d_x=" << d_x << " d_y=" << d_y << " d_z=" << d_z;

      if(boundary_ == FULL)
      {
        if(boundary_ == FULL && dim == 3)
          throw Exception("for 3D no 'full' wave vector sampling.");

        int root = sqrt(steps) + 0.5;
        for(int y = 0; y < root; y++) {
          for(int x = 0; x < root; x++) {
            int idx = y * root + x;
            wave_vectors[idx].Resize(3);
            wave_vectors[idx][0] = (M_PI/d_x) * (x/(root-1));
            wave_vectors[idx][1] = (M_PI/d_y) * (y/(root-1));
           }
         }
      }
      if(boundary_ == HORIZONZAL)
      {
        assert(edges == 1);
        for(int i = 0; i < steps/edges; i++) {
          wave_vectors[i].Resize(3);
          wave_vectors[i][0] = (i*edges*M_PI/(steps-1)) / d_x; // to reach PI,0,0 which would with symmetric and quadrant be job of the vertical line
          wave_vectors[i][1] = 0.0;
          wave_vectors[i][2] = 0.0;
          LOG_DBG2(efd) << "FWV horizontal i=" << i << " -> " << wave_vectors[i][0];
        }
      }
      if(boundary_ == SYMMETRIC || boundary_ == QUADRANT)
      {
        // we don't repeat the corner points, the are calculated at the start of a line
        // start with horizontal 2D and 3D is the same
        for(int i = 0; i < steps/edges; i++) {
          wave_vectors[i].Resize(3);
          wave_vectors[i][0] = (i*edges*M_PI/steps) / d_x;
          wave_vectors[i][1] = 0.0;
          wave_vectors[i][2] = 0.0;
        }
        // vertical on the right side
        for(int i = 0; i < steps/edges; i++) {
          wave_vectors[steps/edges + i].Resize(3);
          wave_vectors[steps/edges + i][0] = M_PI/d_x;
          wave_vectors[steps/edges + i][1] = (i*edges*M_PI/steps) / d_y;
          wave_vectors[steps/edges + i][2] = 0.0;
        }
        if(dim == 2)
        {
          if(boundary_ == SYMMETRIC)
          {
            // diagonal back
            for(int i = 0; i < steps/edges; i++) {
              wave_vectors[2*steps/edges + i].Resize(3);
              wave_vectors[2*steps/edges + i][0] = M_PI/d_x * (1.0-i * (double) edges/steps);
              wave_vectors[2*steps/edges + i][1] = M_PI/d_y * (1.0-i * (double) edges/steps);
              wave_vectors[2*steps/edges + i][2] = 0.0;
            }
          }
          else // boundary_ == QUADRANT
          {
            // top from right to left
            for(int i = 0; i < steps/edges; i++) {
              wave_vectors[2*steps/edges + i].Resize(3);
              wave_vectors[2*steps/edges + i][0] = M_PI/d_x * (1.0-i * (double) edges/steps);
              wave_vectors[2*steps/edges + i][1] = M_PI/d_y;
              wave_vectors[2*steps/edges + i][2] = 0.0;
            }
            // down at the left side
            for(int i = 0; i < steps/edges; i++) {
              wave_vectors[3*steps/edges + i].Resize(3);
              wave_vectors[3*steps/edges + i][0] = 0;
              wave_vectors[3*steps/edges + i][1] = M_PI/d_y * (1.0-i * (double) edges/steps);
              wave_vectors[3*steps/edges + i][2] = 0.0;
            }
          } // end 2D quadrant
        } // end 2D case
        else // 3D case
        {
          assert(boundary_ == SYMMETRIC); // the other stuff ist not yet implemented
          for(int i = 0; i < steps/edges; i++) {
            wave_vectors[2*steps/edges + i].Resize(3);
            wave_vectors[2*steps/edges + i][0] = M_PI/d_x;
            wave_vectors[2*steps/edges + i][1] = M_PI/d_y;
            wave_vectors[2*steps/edges + i][2] = (i*edges*M_PI/steps) / d_z;
          }
          for(int i = 0; i < steps/edges; i++) {
            wave_vectors[3*steps/edges + i].Resize(3);
            wave_vectors[3*steps/edges + i][0] = M_PI/d_x * (1.0-i * (double) edges/steps);
            wave_vectors[3*steps/edges + i][1] = M_PI/d_y * (1.0-i * (double) edges/steps);
            wave_vectors[3*steps/edges + i][2] = M_PI/d_z * (1.0-i * (double) edges/steps);
          }
        } // end 3D case
      } // end non-Full case
    } // end ibz case
    else
    {
      wave_vectors.Resize(lst.GetSize());

      for(unsigned int i = 0; i < lst.GetSize(); i++)
      {
        PtrParamNode wv = lst[i];
        Vector<double>& p = wave_vectors[i];
        p.Resize(dim);
        p[0] = wv->Get("x")->As<double>();
        p[1] = wv->Get("y")->As<double>();
        if(dim == 3)
          p[2] = wv->Get("z")->As<double>();
        LOG_DBG3(efd) << "EDF:FWV i=" << i << " p=" << p.ToString();
      }
    }

    LOG_DBG(efd) << "FWV -> " << ToString<double>(wave_vectors, true);

    current_wave_vector_ = wave_vectors[0]; // copy constructor :)

    if(wave_vectors.IsEmpty())
      throw Exception("Bloch mode Eigenfrequency analysis requires at least one wave vector.");
  }

  double EigenFrequencyDriver::GetFrequency(unsigned int idx) const
  {
    return frequency_[idx];;
  }

  double EigenFrequencyDriver::GetDamping(unsigned int idx) const
  {
    if(isQuadratic_)
      return dampingRatio_[idx];// dynamic_cast<Vector<Complex>&>(*eigenFreqs)[idx].real();
    else
      return 0.0;
  }


  // *****************
  //   Solve problem
  // *****************
  void EigenFrequencyDriver::SolveProblem()
  {
    // options not implemented
    ResultHandler* resHandler = domain_->GetResultHandler();

    // we write the results only if we are not optimization. Optimization writes the results by itself via calling StoreResults().
    // But we call PrintResults() to do the output to the info.xml even in optimization

    // we have to estimate the number of steps as we might loop over bloch modes
    unsigned int n = numFreq_ * (isBloch_ ? blochSteps_ : 1) + 1; // one for safety

    // notify resultHandler about beginning of new sequence step
    // see comments in StaticDriver::SolveProblem() for the interplay with optimization
    if(!domain->GetOptimization())
       resHandler->BeginMultiSequenceStep(sequenceStep_, analysis_, n); // optimization does it by itself

    if(writeAllSteps_ || isPartOfSequence_) // only for allowPostProc
      simState_->BeginMultiSequenceStep(sequenceStep_, analysis_);

    // ------------------------------
    // Phase 1: calculate eigenvalues( generalized problem)
    // ------------------------------
    // Trigger calculation
    ptPDE_->WriteGeneralPDEdefines();
    BaseSolveStep* step = ptPDE_->GetSolveStep();

    // set the mode normalization

    dynamic_cast<StdSolveStep*>(step)->GetAlgSys()->GetEigenSolver()->SetModeNormalization(modeNormalization_);

    // Start of FEAST section and layout for new, more flexible structure which does not rely on the intermediate
    // functions in algSys, StdSolveStep, ...
    BaseEigenSolver::EigenSolverType solverType = dynamic_cast<StdSolveStep*>(step)->GetAlgSys()->GetEigenSolver()->eigenSolverType_;



      if (maxVal_>0 || solverType==BaseEigenSolver::PALM || solverType==BaseEigenSolver::QUADRATIC || solverType==BaseEigenSolver::EXTERNAL) { // use only for feast and PALM -> TODO: remove this if
      // here we should - do the necessary computation depending on the problem type
      StdSolveStep* sstep = dynamic_cast<StdSolveStep*>(step);
      BaseEigenSolver* eigenSolver = sstep->GetAlgSys()->GetEigenSolver();

      // initialize AlgSys
      sstep->GetAlgSys()->InitSol();
      sstep->GetAlgSys()->InitMatrix();
      sstep->GetAssemble()->AssembleMatrices();
      sstep->GetAlgSys()->ExportLinSys(true,false,false); // export the setup
      // determine which EV problem to set up: we make a generalized one
      // We should probably check if we have both matrices, but currently I do not now a case where we do not ...

      // the old stuff should be moved here, after adaption to the new structure of BaseEigenSolver
      if(isQuadratic_)
      {
        // determine which EV problem to set up: we make a generalized one
        // We should probably check if we have both matrices, but currently I do not now a case where we do not ...
        SBM_Matrix* massMat = sstep->GetAlgSys()->GetMatrix(MASS);
        SBM_Matrix* stiffMat = sstep->GetAlgSys()->GetMatrix(STIFFNESS);
        SBM_Matrix* dampMat = sstep->GetAlgSys()->GetMatrix(DAMPING);
        UInt i = massMat->GetNumCols();
        if (i>1) {
            EXCEPTION("only implemented for SBM matrices with a single block")
        }
        // check matrix dimensions
        assert( massMat->GetNumCols()==massMat->GetNumRows() );
        assert( stiffMat->GetNumCols()==stiffMat->GetNumRows() );
        assert( dampMat->GetNumCols()==dampMat->GetNumRows() );
        // * the quadratic EVP should go somewhere else, as it required a completely different handling of the results
        // setup the eigen solver (problem type is determined in Setup based on matrix properties)
        sstep->GetAlgSys()->GetEigenSolver()->Setup(*(stiffMat->GetPointer(0,0)),*(dampMat->GetPointer(0,0)),*(massMat->GetPointer(0,0)));
      }
      else
      {

        SBM_Matrix* massMat = sstep->GetAlgSys()->GetMatrix(MASS);
        SBM_Matrix* stiffMat = sstep->GetAlgSys()->GetMatrix(STIFFNESS);
        UInt i = massMat->GetNumCols();
        if (i>1) {
            EXCEPTION("only implemented for SBM matrices with a single block")
        }
        // check matrix dimensions
        assert( massMat->GetNumCols()==massMat->GetNumRows() );
        assert( stiffMat->GetNumCols()==stiffMat->GetNumRows() );
        // * the quadratic EVP should go somewhere else, as it required a completely different handling of the results
        // setup the eigen solver (problem type is determined in Setup based on matrix properties)
        sstep->GetAlgSys()->GetEigenSolver()->Setup(*(stiffMat->GetPointer(0,0)),*(massMat->GetPointer(0,0)),isBloch_);
        // check if the eigenvalues will be complex

      }
      bool complexEV = eigenSolver->HasComplexEigenvalues();
      if ((minVal_!=0 || maxVal_!=0) || solverType==BaseEigenSolver::EXTERNAL) { // we have an interval
        if (complexEV || ((solverType==BaseEigenSolver::FEAST)&&(isQuadratic_))) {
          Vector<Complex> evals,errs;
          if (numFreq_ > 0 || freqShift_ != 0) // case for the external eigensolver
          {
            sstep->GetAlgSys()->GetEigenSolver()->CalcEigenValues(evals,errs,numFreq_, 2*M_PI*freqShift_);
          }
          else
          {
            sstep->GetAlgSys()->GetEigenSolver()->CalcEigenValues(evals,errs,minVal_,maxVal_);
          }
          eigsRe_.Resize(evals.GetSize());
          eigsIm_.Resize(evals.GetSize());
          for (int i=0;i<(int)evals.GetSize();i++) {
              eigsRe_[i] = evals[i].real();
              eigsIm_[i] = evals[i].imag();
          }
          if (isQuadratic_){
            Vector<double> fd;
            QuadEig2FreqDamp(eigsRe_,eigsIm_,fd,frequency_,dampingRatio_);
          }
          else {
            Eig2FreqDamp(evals,frequency_,dampingRatio_);
          }
        }
        else {
          Vector<Double> evals,errs;
          if (numFreq_ > 0 || freqShift_ != 0) // case for the external eigensolver
          {
            sstep->GetAlgSys()->GetEigenSolver()->CalcEigenValues(evals,errs,numFreq_, 2*M_PI*freqShift_);
          }
          else
          {
            sstep->GetAlgSys()->GetEigenSolver()->CalcEigenValues(evals,errs,minVal_,maxVal_);
          }
          eigsRe_.Resize(evals.GetSize());
          eigsRe_ = evals;
          eigsIm_.Resize(evals.GetSize(),0.0);
          Eig2Freq(evals,frequency_);
          dampingRatio_.Resize(evals.GetSize(),0.0);
        }

        // info output: ToDo: make this pretty
        std::cout << "eigsRe = " << eigsRe_.ToString() << "\n";
        std::cout << "eigsIm = " << eigsIm_.ToString() << "\n";
        std::cout << "Frequency = " << frequency_.ToString() << "\n";
        std::cout << "dampingRatio = " << dampingRatio_.ToString() << "\n";
        SortModes();
        int n = 15; // field width
        cout << "\n";
        cout << " Mode | ";
        cout << setw(n) << "Frequency in Hz" << " | ";
        cout << setw(n) << "Damping Ratio" << " | ";
        cout << setw(n) << "Re(lamda)" << " | ";
        cout << setw(n) << "Im(lambda)" << "\n";
        // plot sorted by frequency
        for(unsigned int i=0; i < modeOrder_.GetSize(); i++) {
          cout << setw(5) << i+1 << " | ";
          cout << setw(n) << frequency_[modeOrder_[i]]<< " | ";
          cout << setw(n) << dampingRatio_[modeOrder_[i]]<< " | ";
          cout << setw(n) << eigsRe_[modeOrder_[i]]<< " | ";
          cout << setw(n) << eigsIm_[modeOrder_[i]] << "\n";
        }

        // export solution
        PtrParamNode els = sstep->GetAlgSys()->GetExportLinSysParam();
        if (els) {
          if(els->Get("solution")->As<bool>()) {
            BaseMatrix::OutputFormat vec_format = MatrixOutputFormatEnum.Parse(els->Get("vecFormat")->As<std::string>());
            std::string base = els->Has("baseName") ? els->Get("baseName")->As<std::string>() : progOpts->GetSimName();
            if(domain->GetDriver()->GetAnalysisId().ToString(true) != ""){
              base += "_" + domain->GetDriver()->GetAnalysisId().ToString(true);
            }
            for (UInt i=0; i< frequency_.GetSize();i++) {
              Vector<Complex> mode;

              // TU Wien Variant with normalized eigenmodes
              sstep->GetAlgSys()->GetEigenSolver()->GetNormalizedEigenMode(modeOrder_[i],mode);
              mode.Export( base + "_mode_" + std::to_string(i+1),vec_format);
              sstep->GetAlgSys()->GetEigenSolver()->GetNormalizedEigenMode(modeOrder_[i],mode,false);

              // sharedopt variant with non-normalized modes
              /*sstep->GetAlgSys()->GetEigenSolver()->GetEigenMode(i,mode);
              mode.Export( base + "_mode_" + std::to_string(i+1),vec_format);
              sstep->GetAlgSys()->GetEigenSolver()->GetEigenMode(i,mode,false);*/

              mode.Export( base + "_mode-left_" + std::to_string(i+1),vec_format);
            }
          }
        }
      }

      else if ( numFreq_ > 0 || freqShift_ != 0  ) {
        // we have num + shift
        // check if the eigenvalues will be complex
        if(isQuadratic_)
        {
          if(solverType==BaseEigenSolver::PALM)
            {
              Vector<Complex> errs,evals;//PALM wants an complex errs vector, while arpack wants a double
              sstep->GetAlgSys()->GetEigenSolver()->CalcEigenValues( evals, errs, numFreq_, 2*M_PI*freqShift_ );
              eigsRe_.Resize(evals.GetSize());
              eigsIm_.Resize(evals.GetSize());
              for (int i=0;i<(int)evals.GetSize();i++) {
                eigsRe_[i] = evals[i].real();
                eigsIm_[i] = evals[i].imag();
            }}
          else if(solverType==BaseEigenSolver::QUADRATIC)
            {
              Vector<Double> errs;
              Vector<Double> evals;
              sstep->GetAlgSys()->GetEigenSolver()->CalcEigenValues( evals, errs, numFreq_, 2*M_PI*freqShift_ );
              eigsRe_.Resize(evals.GetSize());
              eigsIm_.Resize(evals.GetSize());
              for (int i=0;i<(int)evals.GetSize();i++) {
                eigsRe_[i] = evals[i]; //TODO:arpack expects only real valued solution, which does not correspond to the problem

            }}
          else
            {
              Vector<Double> errs;
              Vector<Complex> evals;
              sstep->GetAlgSys()->GetEigenSolver()->CalcEigenValues( evals, errs, numFreq_, 2*M_PI*freqShift_ );
              eigsRe_.Resize(evals.GetSize());
              eigsIm_.Resize(evals.GetSize());
              for (int i=0;i<(int)evals.GetSize();i++) {
                eigsRe_[i] = evals[i].real();
                eigsIm_[i] = evals[i].imag();
            }}



          Vector<Double> frequency_damped_;
          QuadEig2FreqDamp(eigsRe_, eigsIm_ ,frequency_, frequency_damped_,dampingRatio_);

          std::cout << "eigsRe = " << eigsRe_.ToString() << "\n";
          std::cout << "eigsIm = " << eigsIm_.ToString() << "\n";
          std::cout << "Undamped Frequency = " << frequency_.ToString() << "\n";
          std::cout << "Damped Frequency = " << frequency_damped_.ToString() << "\n";
          std::cout << "dampingRatio = " << dampingRatio_.ToString() << "\n";
          SortModes();
          int n = 35; // field width
          cout << "\n";
          cout << " Mode | ";
          cout << setw(n) << "Frequency (undamped) in Hz" << " | ";
          cout << setw(n) << "Frequency (damped) in Hz" << " | ";
          cout << setw(n) << "Damping Ratio" << " | ";
          cout << setw(n) << "Re(lamda)" << " | ";
          cout << setw(n) << "Im(lambda)" << "\n";
          // plot sorted by frequency
          for(unsigned int i=0; i < modeOrder_.GetSize(); i++) {
            cout << setw(5) << i+1 << " | ";
            cout << setw(n) << frequency_[modeOrder_[i]]<< " | ";
            cout << setw(n) << frequency_damped_[modeOrder_[i]]<< " | ";
            cout << setw(n) << dampingRatio_[modeOrder_[i]]<< " | ";
            cout << setw(n) << eigsRe_[modeOrder_[i]]<< " | ";
            cout << setw(n) << eigsIm_[modeOrder_[i]] << "\n";
          }

          // export solution
          PtrParamNode els = sstep->GetAlgSys()->GetExportLinSysParam();
          if (els) {
            if(els->Get("solution")->As<bool>()) {
              BaseMatrix::OutputFormat vec_format = MatrixOutputFormatEnum.Parse(els->Get("vecFormat")->As<std::string>());
              std::string base = els->Has("baseName") ? els->Get("baseName")->As<std::string>() : progOpts->GetSimName();
              if(domain->GetDriver()->GetAnalysisId().ToString(true) != ""){
                base += "_" + domain->GetDriver()->GetAnalysisId().ToString(true);
              }
              for (UInt i=0; i< frequency_.GetSize();i++) {
                Vector<Complex> mode;

                // TU Wien Variant with normalized eigenmodes
                sstep->GetAlgSys()->GetEigenSolver()->GetNormalizedEigenMode(i,mode);
                mode.Export( base + "_mode_" + std::to_string(i+1),vec_format);
                sstep->GetAlgSys()->GetEigenSolver()->GetNormalizedEigenMode(i,mode,false);
                mode.Export( base + "_mode-left_" + std::to_string(i+1),vec_format);
              }
            }
          }
        }
      }
      else {
        EXCEPTION("this case should not be possible in the XML schema")
      }
    }
    else { // the old stuff to (re)move

      if(isBloch_)
      {
        assert(!wave_vectors.IsEmpty());
        LOG_DBG(efd) << "SP: #wv=" << wave_vectors.GetSize();
        for(unsigned int i = 0; i < wave_vectors.GetSize(); i++)
        {
          ComputeBlochWaveVector(i);
        }
      }
      else {
        // the eigenfrequencies are complex in the quadratic case or in bloch mode
        eigenFreqs->Init();
        errBounds_.Resize(numFreq_);

        if(isQuadratic_)
        {
          Vector<Complex> ef = Vector<Complex>();
          ef.Resize(numFreq_);
          step->CalcEigenFrequencies(ef, errBounds_, numFreq_, freqShift_, sort_, isBloch_);
          frequency_.Resize(ef.GetSize());
          dampingRatio_.Resize(ef.GetSize());
          for (int i=0; i<(int)ef.GetSize();i++){
              frequency_[i] = ef[i].imag()/(2*M_PI);
              dampingRatio_[i] = ef[i].real();
          }
          PrintResult();
        }
        else // real generalized
        {
          Vector<Double>& ef = dynamic_cast<Vector<Double>& >(*eigenFreqs);
          step->CalcEigenFrequencies(ef, errBounds_, numFreq_, freqShift_, sort_);
          PrintResult();
        }
      }
    }// end old stuff

    // in optimization we write the results via StoreResults() because
    // we don't necessarily write every forward step.
    if(!domain->GetOptimization()) // in other words: if not optimization
    {
      StoreResults(1, -1.0);
      handler_->FinishMultiSequenceStep();

      if(!isPartOfSequence_)
        handler_->Finalize(); // to be called only once in a HDF5 lifetime!
    }

    if(writeAllSteps_ || isPartOfSequence_)
      simState_->FinishMultiSequenceStep(true);
  }


  unsigned int EigenFrequencyDriver::GetCurrentWaveVectorIndex() const
  {
    for(unsigned int i = 0; i < wave_vectors.GetSize(); i++)
      if(wave_vectors[i] == current_wave_vector_)
        return i;

    assert(false);
    return 0;
  }

  void EigenFrequencyDriver::ComputeBlochWaveVector(int wave_vector_step)
  {
    eigenFreqs->Init();
    errBounds_.Resize(numFreq_);

    ptPDE_->GetSolveStep()->SetActFreq(0.0); // otherwise it is the last freq from the prev EVA

    current_wave_vector_ = wave_vectors[wave_vector_step]; // StrainOperatorBloch2D has a pointer to current_wave_vector

    LOG_DBG(efd) << "CBWV wvs=" << wave_vector_step << " wv=" << current_wave_vector_.ToString();

    Vector<Complex> ef = Vector<Complex>(numFreq_);//dynamic_cast<Vector<Complex>& >(*eigenFreqs);
    ptPDE_->GetSolveStep()->CalcEigenFrequencies(ef , errBounds_, numFreq_, freqShift_, sort_, isBloch_);

    // put the real part into the "frequency" vector -> it should be the eigenvalue actually
    frequency_.Resize(ef.GetSize());
    for (int i=0;i<(int)ef.GetSize();i++){
        frequency_[i] = ef[i].real();
    }

    PrintResult(wave_vector_step);
    // we need to calculate and output results before the displacements are overwritten.
    // only if not optimization or if optimization when it is evaluatate_initial_design
    // not for the last step as we have a separate store result!
    if((domain->GetOptimization() == NULL || domain->GetOptimization()->GetOptimizerType() == Optimization::EVALUATE_INITIAL_DESIGN) && (wave_vector_step <  (int) wave_vectors.GetSize() - 1))
      StoreResults(wave_vector_step*ef.GetSize(), -1.0);
  }


  unsigned int EigenFrequencyDriver::StoreResults(unsigned int stepNum, double step_val)
  {
    LOG_DBG(efd) << "SR step=" << stepNum << " val=" << step_val;

    unsigned int wvs = isBloch_ ? wave_vectors.GetSize() : 1; // save wave vector size
    unsigned int w = isBloch_ ? GetCurrentWaveVectorIndex() : 0;

    // generates a index-array modeOrder_ containing the mode indices sorted by ascending Frequency value
    SortModes();

    int total = eigenFreqs->GetSize() * wvs;
    string s = boost::lexical_cast<string>(step_val);
    int digs =  boost::lexical_cast<string>(total).size() + s.substr(s.find('.')+1).size();
    double sig = std::pow((float) 10.0, -digs); // 1e-2 -> 10 ^ -2 -> a compiler complained with simply 10

    for(unsigned int fi=0; fi < frequency_.GetSize(); fi++)
    {
      // Phase 2: calculate eigenmodes
      if(writeModes_)
      {
        ptPDE_->GetSolveStep()->SetActStep(fi);
        ptPDE_->GetSolveStep()->SetActFreq(GetFrequency(modeOrder_[fi]));
        ptPDE_->GetDomain()->GetMathParser()->SetValue(MathParser::GLOB_HANDLER, "f", GetFrequency(modeOrder_[fi]) );
        ptPDE_->GetSolveStep()->GetEigenMode(modeOrder_[fi]); // this stores the eigen mode result in AlgSys's sol_

        // stupid paraview needs an increasing series of save_value :(

        double save_value = -1.0;

        if(domain->GetOptimization() && domain->GetOptimization()->GetOptimizerType() != Optimization::EVALUATE_INITIAL_DESIGN)
        {
          // time is step.step_val nr
          save_value = step_val + (w * wvs + fi + 1) * sig; // +1 for one based

          LOG_DBG3(efd) << "SR total=" << total << " digs=" << digs << " sig=" << sig << " count=" << (w * wvs + fi + 1);
        }
        else // for bloch case we label <step>.<nr> from the info.xml
          save_value = isBloch_ ? w + (fi+1.0) / (eigenFreqs->GetSize() < 9 ? 10.0 : 100.0) : std::abs(GetFrequency(modeOrder_[fi]));

        LOG_DBG(efd) << "SR w=" << w << " fi=" << fi << " stepNum=" << stepNum << " save_value=" << save_value;

        handler_->BeginStep(stepNum, save_value);
        ptPDE_->WriteResultsInFile(stepNum, save_value);
        handler_->FinishStep();

        if(writeAllSteps_ || isPartOfSequence_)
          simState_->WriteStep(stepNum, save_value);

        if (!GetResultHandler()->streamOnly)
          stepNum++;
      }
    }

    if (!GetResultHandler()->streamOnly)
      return stepNum-1;
    else
      return stepNum;
  }

  void EigenFrequencyDriver::PrintResult(int wave_vector_step)
  {
    unsigned int dim = domain->GetGrid()->GetDim();

    unsigned int numConverged = eigenFreqs->GetSize();

    // If no frequency at all converged, just leave
    if( numConverged == 0) {
      WARN( "No eigenfrequency converged, so no output will be written to the result files." );
      return;
    }

    // Issue warning, if number of converged eigenvalues differs from
    // number of requested ones
    if( numConverged != numFreq_ ) {
      WARN( "Only " << numConverged << " eigenfrequencies of " 
          << numFreq_ << " converged. To improve convergence, either "
          << "reduce the number of eigenfrequencies or the tolerance." );
    }

    // console output (reduced form bloch)
    if(!domain->GetOptimization())
    {
      if(isBloch_)
      {
        // single line printing
        cout << std::endl << " step=" << wave_vector_step << " wave vector=" << current_wave_vector_.ToString() << " frequencies: ";
      }
      else
      {
        cout << std::endl;
        cout << setw(20) << "Frequency in Hz" << " | ";
        if(isQuadratic_)
          cout << setw(20) << "Damping       " << " | ";
        cout << setw(20) << "Errorbound";
        cout << "\n";
        cout << std::setfill('-') << setw(isQuadratic_ ? 70 : 40) << "" << std::setfill(' ');
        cout << "\n";
      }
    }

    // also log via info node.
    // new result only on detailed output and for bloch case only for step=0
    PtrParamNode res;
    bool append = progOpts->DoDetailedInfo() && (!isBloch_ || wave_vector_step == 0);
    if(append)
      res = info_->Get("result", ParamNode::APPEND);
    else  {
      // create or take the last
      ParamNodeList l = info_->GetList("result");
      res = l.IsEmpty() ? info_->Get("result") : l.Last();
    }
    if(domain->GetOptimization())
      res->Get("iteration")->SetValue(domain->GetOptimization()->GetCurrentIteration());

    // the problem might be that we have no detailed output, hence the res might already exist.
    // then, numConverged might be smaller than an older optimization iteration and the old stuff > numConverged remains
    if(!isBloch_ || wave_vector_step == 0) // don't delete old wave_vectors
      res->GetChildren().Resize(0);

    if(isBloch_)
    {
      res = res->GetByVal("wave_vector", "step", wave_vector_step);
      res->Get("k_x")->SetValue(current_wave_vector_[0]);
      res->Get("k_y")->SetValue(current_wave_vector_[1]);
      if(dim == 3)
        res->Get("k_z")->SetValue(current_wave_vector_[2]);

      std::stringstream ss;
      ss << wave_vector_step << "\t" << current_wave_vector_[0] << "\t" << current_wave_vector_[1] << "\t";
      if(dim == 3)
        ss << current_wave_vector_[2] << "\t";

      bloch_plot_ << ss.str();
      if(wave_vector_step == 0)
        first_plot_line_ = ss.str();
    }

    // is not set for pure simulation
    if(analysis_id_.ToString() != "")
      res->Get("id")->SetValue(analysis_id_.ToString());

    for(unsigned int i=0; i < numConverged; i++)
    {
      double freq = GetFrequency(i);
      double damp = GetDamping(i);

      if(isBloch_)
      {
        std::stringstream ss;
        ss << freq << (i < numConverged-1 ? "\t" : "\n");
        bloch_plot_ << ss.str();
        bloch_plot_.flush();
        if(wave_vector_step == 0) // to be repeated as first line!
          first_plot_line_ += ss.str();
      }

      // command line output only when not optimizing
      if(!domain->GetOptimization())
      {
        if(isBloch_)
        {
          cout << freq << (i < numConverged-1 ? ", " : "");
        }
        else
        {
          cout << setw(20) << freq <<" | ";
          if(isQuadratic_)
            cout << setw(20) << damp << " |  ";
          cout << setw(20) << errBounds_[i] <<  "\n";
        }
      }

      // res got cleared above, so we don't overlap with previous in case numConverged is too small
      PtrParamNode mode = res->Get("mode", ParamNode::APPEND);

      mode->Get("nr")->SetValue(i+1); // not the mode but frequency in list
      mode->Get("frequency")->SetValue(freq,15);
      if(isQuadratic_)
        mode->Get("damping")->SetValue(damp);
      mode->Get("errorbound")->SetValue(errBounds_[i]);

    }

    if(!domain->GetOptimization())
      cout << "\n";

    if(isBloch_ && this->ibz_)
    {
      // repeat the first step at the and of bloch.plot for a full ibz for plotting, not when explicit wave vectors are given
      if(wave_vector_step == (int) wave_vectors.GetSize() - 1 && (boundary_ == SYMMETRIC || boundary_ == QUADRANT)) {
        bloch_plot_ << first_plot_line_;
        bloch_plot_.close();
      }
    }
  }


} // end of namespace
