#include "EigenFrequencyDriver.hh"

#include <iostream>
#include <iomanip>

#include "Driver/SolveSteps/StdSolveStep.hh"

#include "Domain/Domain.hh"

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/SimState.hh"
#include "DataInOut/ResultHandler.hh"
#include "DataInOut/ProgramOptions.hh"
#include "DataInOut/Logging/LogConfigurator.hh"

#include "PDE/StdPDE.hh"

using std::cout;
using std::setw;



namespace CoupledField {

  DECLARE_LOG(efd)
  DEFINE_LOG(efd, "eigenFrequencyDriver")

  // ***************
  //   Constructor
  // ***************
  EigenFrequencyDriver::EigenFrequencyDriver(UInt sequenceStep,
                                             bool isPartOfSequence,
                                             shared_ptr<SimState> state,
                                             Domain* domain,
                                             PtrParamNode paramNode, PtrParamNode infoNode ) 
    : SingleDriver( sequenceStep, isPartOfSequence, state, domain,
                    paramNode, infoNode) {

    // set correct analysistype
    analysis_ = BasePDE::EIGENFREQUENCY;

    numFreq_ = 0;
    freqShift_ = 0.0;
    writeModes_ = true;
    isQuadratic_ = false;
    isBloch_ = false;
    ibz_     = false;
    // replace with a concrete element
    param_ = param_->Get("eigenFrequency");
    info_ = info_->Get("eigenFrequency");
  }


  void EigenFrequencyDriver::Init(bool restart)
  {
    // read required parameters from parameter node
    param_->GetValue( "numModes", numFreq_ );
    param_->GetValue( "freqShift", freqShift_ );
    param_->GetValue( "writeModes", writeModes_, ParamNode::PASS );
    param_->GetValue( "isQuadratic", isQuadratic_, ParamNode::PASS );
    // read flag if all results should get written to database file section
    // to allow e.g. for general postprocessing or result extraction    
    param_->GetValue("allowPostProc", writeAllSteps_, ParamNode::PASS );
    isBloch_ = param_->Has("bloch");
    if(isBloch_)
    {
      FillWaveVectors(param_->Get("bloch"));
      std::string file = progOpts->GetSimName() + ".bloch.dat";
      bloch_plot_.open(file.c_str(), std::ios::out);
      bloch_plot_ << "#step\tk_x\tk_y";
      if(domain->GetGrid()->GetDim() == 3)
        bloch_plot_ << "\tk_z";
      bloch_plot_ << "\t1.mode\t2.mode\t..." << std::endl;
    }

    if(isQuadratic_ && isBloch_)
      throw Exception("Bloch mode Eigenfrequency analysis not implemented for quadratic form");

    InitializePDEs();
  }


  // **********************
  //   Default destructor
  // **********************
  EigenFrequencyDriver::~EigenFrequencyDriver()
  {
    if(isBloch_)
      bloch_plot_.close();
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

      bool do_boundary = bloch_pn->Get("ibz/sample")->As<std::string>() == "boundary";

      blochSteps_     = bloch_pn->Get("ibz/steps")->As<int>();
      int steps = blochSteps_;

      int edges = dim == 2 ? 3 : 4; // GAMMA -> X -> M -> [GAMMA | R -> GAMMA]

      if(!do_boundary && dim == 3)
              throw Exception("for 3D only 'do_boundary' sampling of the wave vector.");

      if((do_boundary && ((steps % edges) != 0 || steps < edges)) || (!do_boundary && (int) sqrt(steps) * (int) sqrt(steps) != steps))
        EXCEPTION("bloch mode ibz/steps need to a multiple of " << edges << " for boundary sampling or a square number for full sampling.");

      // we need the unit cell dimensions to scale the wave vector from 0 .. pi/d where d is unit cell dimension -> Hussein;2009
      Matrix<double>& box = domain->GetGrid()->CalcGridBoundingBox(NULL, true); // force 3D
      double d_x = box[0][1]-box[0][0];
      double d_y = box[1][1]-box[1][0];
      double d_z = box[2][1]-box[2][0]; // 0.0 for 2D
      wave_vectors_.Resize(steps);

      LOG_DBG(efd) << "FWV box: " << box.ToString();
      LOG_DBG(efd) << "FWV box d_x=" << d_x << " d_y=" << d_y << " d_z=" << d_z;

      if(do_boundary)
      {
        // we don't repeat the corner points, the are calculated at the start of a line
        for(int i = 0; i < steps/edges; i++) {
          wave_vectors_[i].Resize(3);
          wave_vectors_[i][0] = (i*edges*PI/steps) / d_x;
          wave_vectors_[i][1] = 0.0;
          wave_vectors_[i][2] = 0.0;
        }
        for(int i = 0; i < steps/edges; i++) {
          wave_vectors_[steps/edges + i].Resize(3);
          wave_vectors_[steps/edges + i][0] = PI/d_x;
          wave_vectors_[steps/edges + i][1] = (i*edges*PI/steps) / d_y;
          wave_vectors_[steps/edges + i][2] = 0.0;
        }
        if(dim == 2)
          for(int i = 0; i < steps/edges; i++) {
            wave_vectors_[2*steps/edges + i].Resize(3);
            wave_vectors_[2*steps/edges + i][0] = PI/d_x * (1.0-i * (double) edges/steps);
            wave_vectors_[2*steps/edges + i][1] = PI/d_y * (1.0-i * (double) edges/steps);
            wave_vectors_[2*steps/edges + i][2] = 0.0;
          }
        else
        {
          for(int i = 0; i < steps/edges; i++) {
            wave_vectors_[2*steps/edges + i].Resize(3);
            wave_vectors_[2*steps/edges + i][0] = PI/d_x;
            wave_vectors_[2*steps/edges + i][1] = PI/d_y;
            wave_vectors_[2*steps/edges + i][2] = (i*edges*PI/steps) / d_z;
          }
          for(int i = 0; i < steps/edges; i++) {
            wave_vectors_[3*steps/edges + i].Resize(3);
            wave_vectors_[3*steps/edges + i][0] = PI/d_x * (1.0-i * (double) edges/steps);
            wave_vectors_[3*steps/edges + i][1] = PI/d_y * (1.0-i * (double) edges/steps);
            wave_vectors_[3*steps/edges + i][2] = PI/d_z * (1.0-i * (double) edges/steps);
          }
        }
      }
      else
      {
        int root = sqrt(steps) + 0.5;
        for(int y = 0; y < root; y++) {
          for(int x = 0; x < root; x++) {
            int idx = y * root + x;
            wave_vectors_[idx].Resize(3);
            wave_vectors_[idx][0] = (PI/d_x) * (x/(root-1));
            wave_vectors_[idx][1] = (PI/d_y) * (y/(root-1));
           }
         }
      }
    }
    else
    {
      wave_vectors_.Resize(lst.GetSize());

      for(unsigned int i = 0; i < lst.GetSize(); i++)
      {
        PtrParamNode wv = lst[i];
        Vector<double>& p = wave_vectors_[i];
        p.Resize(dim);
        p[0] = wv->Get("x")->As<double>();
        p[1] = wv->Get("y")->As<double>();
        if(dim == 3)
          p[2] = wv->Get("z")->As<double>();
        LOG_DBG3(efd) << "EDF:FWV i=" << i << " p=" << p.ToString();
      }
    }

    LOG_DBG(efd) << "FWV -> " << ToString<double>(wave_vectors_, true);

    current_wave_vector_ = wave_vectors_[0]; // copy constructor :)

    if(wave_vectors_.IsEmpty())
      throw Exception("Bloch mode Eigenfrequency analysis requires at least one wave vector.");
  }

  template <class TYPE>
  void EigenFrequencyDriver::PrintResult(SingleVector* freq_ptr, Vector<Double>& errBounds,
                                         ResultHandler* resHandler, UInt numConverged, int wave_vector_step)
  {
    Vector<TYPE>& eigenFreqs = dynamic_cast<Vector<TYPE>&>(*freq_ptr);
    unsigned int dim = domain->GetGrid()->GetDim();

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
    if(isBloch_)
    {
      // single line printing
      cout << std::endl << " step=" << wave_vector_step << " wave vector=" << current_wave_vector_.ToString() << " frequencies: ";
      // also plot
      std::stringstream ss;
      ss << wave_vector_step << "\t" << current_wave_vector_[0] << "\t" << current_wave_vector_[1] << "\t";
      if(dim == 3)
        ss << current_wave_vector_[2] << "\t";

      bloch_plot_ << ss.str();
      if(wave_vector_step == 0)
        first_plot_line_ = ss.str();
    }
    else
    {
      cout << std::endl << std::endl;
      cout << setw(20) << "Frequency [Hz]" << " | ";
      if(isQuadratic_)
        cout << setw(20) << "Damping       " << " | ";
      cout << setw(20) << "Errorbound";
      cout << "\n";
      cout << std::setfill('-') << setw(isQuadratic_ ? 70 : 40) << "" << std::setfill(' ');
      cout << "\n";
    }



    for( UInt i=0; i<numConverged; i++ )
    {
      // allways complex for real templated case
      Complex ev = (Complex) eigenFreqs[i];

      Double freq = isQuadratic_ ? ev.imag()/(8.0*atan(1.0)) : ev.real();
      Double damp = ev.real();

      if(isBloch_)
      {
        cout << freq << (i < numConverged-1 ? ", " : "");
        std::stringstream ss;
        ss << freq << (i < numConverged-1 ? "\t" : "\n");
        bloch_plot_ << ss.str();
        if(wave_vector_step == 0) // to be repeated as first line!
          first_plot_line_ += ss.str();
      }
      else
      {
        cout << setw(20) << freq <<" | ";
        if(isQuadratic_)
          cout << setw(20) << damp << " |  ";
        cout << setw(20) << errBounds[i] <<  "\n";
      }

      // also log via info node
      // BLOCH CECK: PtrParamNode result = info_->Get("result",ParamNode::APPEND);
      PtrParamNode mode = info_->Get("result")->Get("mode",ParamNode::APPEND);
      if(isBloch_) {
        mode->Get("step")->SetValue(wave_vector_step);
        mode->Get("k_x")->SetValue(current_wave_vector_[0]);
        mode->Get("k_y")->SetValue(current_wave_vector_[1]);
        if(dim == 3)
          mode->Get("k_z")->SetValue(current_wave_vector_[2]);
      }

      mode->Get("nr")->SetValue(i+1); // not the mode but frequency in list

      mode->Get("frequency")->SetValue(freq,15);
      if(isQuadratic_)
        mode->Get("damping")->SetValue(damp);
      mode->Get("errorbound")->SetValue(errBounds[i]);

      // Phase 2: calculate eigenmodes
      if ( writeModes_ == true )
      {
        ptPDE_->GetSolveStep()->SetActStep(i);
        ptPDE_->GetSolveStep()->SetActFreq(std::abs(freq));
        ptPDE_->GetSolveStep()->CalcEigenMode(i);
        // in bloch mode analysis we solve the numFreq_ eigenmodes many_ times for different wave vectors
        unsigned int save_step = isBloch_ ? wave_vector_step * numFreq_  + i: i + 1;
        // for bloch case we label <step>.<nr> from the info.xml
        double       save_value = isBloch_ ? wave_vector_step + (i+1.0) / 100.0 : std::abs(freq);
        resHandler->BeginStep(save_step, save_value);
        ptPDE_->WriteResultsInFile(save_step, save_value);
        resHandler->FinishStep();
        if( writeAllSteps_ )
          simState_->WriteStep(save_step, std::abs(eigenFreqs[i]) );
      }
    }

    if(isBloch_)
    {
      // repeat the first step at the and of bloch.plot
      if(wave_vector_step == Integer(wave_vectors_.GetSize()) - 1)
         bloch_plot_ << first_plot_line_;
      bloch_plot_.flush();
    }

  }
  
  // *****************
  //   Solve problem
  // *****************
  void EigenFrequencyDriver::SolveProblem() {
    // options not implemented
    
    analysis_id_ = info_->Get(ParamNode::PN_PROCESS);

    ResultHandler* resHandler = domain_->GetResultHandler();

    // we have to estimate the number of steps as we might loop over bloch modes
    unsigned int n = numFreq_ * (isBloch_ ? blochSteps_ : 1);
    // notify resultHandler about beginning of new sequence step
    resHandler->BeginMultiSequenceStep( sequenceStep_, analysis_, n);
    if( writeAllSteps_ )
      simState_->BeginMultiSequenceStep( sequenceStep_, analysis_ );


    // ------------------------------
    // Phase 1: calculate eigenvalues( generalized problem)
    // ------------------------------
    
    // the eigenfrequencies are complex in the quadratic case or in bloch mode
    SingleVector* eigenFreqs = NULL;
    if(isQuadratic_ || isBloch_) eigenFreqs = new Vector<Complex>(numFreq_);
                            else eigenFreqs = new Vector<Double>(numFreq_);

    Vector<Double> errBounds( numFreq_ );

    // Trigger calculation 
    UInt conv = 0;
    ptPDE_->WriteGeneralPDEdefines();
    BaseSolveStep* step = ptPDE_->GetSolveStep();

    if(isBloch_)
    {
      assert(!wave_vectors_.IsEmpty() && current_wave_vector_[0] == wave_vectors_[0][0]);
      for(unsigned int i = 0; i < wave_vectors_.GetSize(); i++)
      {
        ptPDE_->GetSolveStep()->SetActFreq(0.0); // otherwise it is the last freq from the prev EVA

        current_wave_vector_ = wave_vectors_[i]; // StrainOperatorBloch2D has a pointer to current_wave_vector

        LOG_DBG(efd) << "SP i=" << i << " wv=" << current_wave_vector_.ToString();

        Vector<Complex>& ef = dynamic_cast<Vector<Complex>& >(*eigenFreqs);
        conv = step->CalcEigenFrequencies(ef , errBounds,numFreq_, freqShift_, isBloch_);
        PrintResult<Complex>(eigenFreqs, errBounds, resHandler, conv, i);
      }
    }
    if(isQuadratic_ && !isBloch_)
    {
      Vector<Complex>& ef = dynamic_cast<Vector<Complex>& >(*eigenFreqs);
      conv = step->CalcEigenFrequencies(ef, errBounds,numFreq_, freqShift_, isBloch_);
      PrintResult<Complex>(eigenFreqs, errBounds, resHandler, conv);
    }
    if(!isQuadratic_ && !isBloch_) // real generalized
    {
      Vector<Double>& ef = dynamic_cast<Vector<Double>& >(*eigenFreqs);
      conv = step->CalcEigenFrequencies(ef, errBounds,numFreq_, freqShift_ );
      PrintResult<Double>(eigenFreqs, errBounds, resHandler, conv);
    }
    
    // notify resultHandler about finishing of current sequence step
    resHandler->FinishMultiSequenceStep();
    if( !isPartOfSequence_ ) {
      resHandler->Finalize();
    }
    if( writeAllSteps_ )
      simState_->FinishMultiSequenceStep(true);

    delete eigenFreqs; // hopefully this is no problem for optimization
  }
  
  
} // end of namespace
