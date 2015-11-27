#include "EigenFrequencyDriver.hh"

#include <iostream>
#include <iomanip>
#include <math.h>

#include "Driver/SolveSteps/StdSolveStep.hh"
#include "Domain/Domain.hh"
#include "Optimization/Optimization.hh"

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/SimState.hh"
#include "DataInOut/ResultHandler.hh"
#include "DataInOut/ProgramOptions.hh"
#include "DataInOut/Logging/LogConfigurator.hh"

#include "PDE/StdPDE.hh"

using std::cout;
using std::setw;
using std::string;



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
    save_step_ = 1;

    // replace with a concrete element
    param_ = param_->Get("eigenFrequency");
    info_ = info_->Get("eigenFrequency");

    isBloch_ = DoBloch(param_);

  }

  EigenFrequencyDriver::~EigenFrequencyDriver()
  {
    if(isBloch_)
      bloch_plot_.close();

    delete eigenFreqs;
    eigenFreqs = NULL;
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
      SetupBlochPlot();
    }

    if(isQuadratic_ && isBloch_)
      throw Exception("Bloch mode Eigenfrequency analysis not implemented for quadratic form");

    // the eigenfrequencies are complex in the quadratic case or in bloch mode
    if(isQuadratic_ || isBloch_) eigenFreqs = new Vector<Complex>(numFreq_);
                            else eigenFreqs = new Vector<Double>(numFreq_);

    InitializePDEs();
  }

  void EigenFrequencyDriver::SetupBlochPlot()
  {
    string file = progOpts->GetSimName();
    if(progOpts->DoDetailedInfo() && domain->GetOptimization())
      file += "_iter_" + boost::lexical_cast<string>(domain->GetOptimization()->GetCurrentIteration());
    file += ".bloch.dat";

    bloch_plot_.close();
    bloch_plot_.open(file.c_str(), std::ios::out);
    bloch_plot_ << "#step\tk_x\tk_y";
    if(domain->GetGrid()->GetDim() == 3)
      bloch_plot_ << "\tk_z";
    bloch_plot_ << "\t1.mode\t2.mode\t..." << std::endl;
  }

  unsigned int EigenFrequencyDriver::GetNumBlochWave(PtrParamNode node)
  {
    if(!node->Has("bloch"))
      return 0;

    if(node->Has("bloch/waveVector/ibz/steps"))
      return node->Get("bloch/waveVector/ibz/steps")->As<unsigned int>();

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

      bool do_boundary = bloch_pn->Get("ibz/sample")->As<string>() == "boundary";

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
      wave_vectors.Resize(steps);

      LOG_DBG(efd) << "FWV box: " << box.ToString();
      LOG_DBG(efd) << "FWV box d_x=" << d_x << " d_y=" << d_y << " d_z=" << d_z;

      if(do_boundary)
      {
        // we don't repeat the corner points, the are calculated at the start of a line
        for(int i = 0; i < steps/edges; i++) {
          wave_vectors[i].Resize(3);
          wave_vectors[i][0] = (i*edges*M_PI/steps) / d_x;
          wave_vectors[i][1] = 0.0;
          wave_vectors[i][2] = 0.0;
        }
        for(int i = 0; i < steps/edges; i++) {
          wave_vectors[steps/edges + i].Resize(3);
          wave_vectors[steps/edges + i][0] = M_PI/d_x;
          wave_vectors[steps/edges + i][1] = (i*edges*M_PI/steps) / d_y;
          wave_vectors[steps/edges + i][2] = 0.0;
        }
        if(dim == 2)
          for(int i = 0; i < steps/edges; i++) {
            wave_vectors[2*steps/edges + i].Resize(3);
            wave_vectors[2*steps/edges + i][0] = M_PI/d_x * (1.0-i * (double) edges/steps);
            wave_vectors[2*steps/edges + i][1] = M_PI/d_y * (1.0-i * (double) edges/steps);
            wave_vectors[2*steps/edges + i][2] = 0.0;
          }
        else
        {
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
        }
      }
      else
      {
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
    }
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
    if(isQuadratic_)
      return dynamic_cast<Vector<Complex>&>(*eigenFreqs)[idx].imag() / (2.0 * M_PI);
    if(isBloch_)
      return dynamic_cast<Vector<Complex>&>(*eigenFreqs)[idx].real();
    else
      return dynamic_cast<Vector<double>&>(*eigenFreqs)[idx];
  }

  double EigenFrequencyDriver::GetDamping(unsigned int idx) const
  {
    if(isQuadratic_)
      return dynamic_cast<Vector<Complex>&>(*eigenFreqs)[idx].real();
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
    unsigned int n = numFreq_ * (isBloch_ ? blochSteps_ : 1) + 1; // one for savety

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

    if(isBloch_)
    {
      assert(!wave_vectors.IsEmpty());
      for(unsigned int i = 0; i < wave_vectors.GetSize(); i++)
      {
        ComputeBlochWaveVector(i);
      }
    }
    else
    {
      // the eigenfrequencies are complex in the quadratic case or in bloch mode
      eigenFreqs->Init();
      errBounds_.Resize(numFreq_);

      if(isQuadratic_)
      {
        Vector<Complex>& ef = dynamic_cast<Vector<Complex>& >(*eigenFreqs);
        step->CalcEigenFrequencies(ef, errBounds_, numFreq_, freqShift_, isBloch_);
        PrintResult();
      }
      else // real generalized
      {
        Vector<Double>& ef = dynamic_cast<Vector<Double>& >(*eigenFreqs);
        step->CalcEigenFrequencies(ef, errBounds_, numFreq_, freqShift_);
        PrintResult();
      }
    }
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

    LOG_DBG(efd) << "CBWV wv=" << current_wave_vector_.ToString();

    Vector<Complex>& ef = dynamic_cast<Vector<Complex>& >(*eigenFreqs);
    ptPDE_->GetSolveStep()->CalcEigenFrequencies(ef , errBounds_, numFreq_, freqShift_, isBloch_);

    PrintResult(wave_vector_step);
    // we need to calculate and output results before the displacements are overwritten.
    // only if not optimization or if optimization when it is evaluatate_initial_design
    // not for the last step as we have a separate store result!
    if((domain->GetOptimization() == NULL || domain->GetOptimization()->GetOptimizerType() == Optimization::EVALUATE_INITIAL_DESIGN) && (wave_vector_step <  (int) wave_vectors.GetSize() - 1))
      StoreResults(wave_vector_step, -1.0);
  }


  void EigenFrequencyDriver::StoreResults(unsigned int stepNum, double step_val)
  {
    // stepNum and step_val are ignored
    LOG_DBG(efd) << "SR step=" << stepNum << " val=" << step_val;

    unsigned int wvs = isBloch_ ? wave_vectors.GetSize() : 1; // save wave vector size
    unsigned int w = isBloch_ ? GetCurrentWaveVectorIndex() : 0;

    for(unsigned int fi=0; fi < eigenFreqs->GetSize(); fi++)
    {
      // Phase 2: calculate eigenmodes
      if(writeModes_)
      {
        ptPDE_->GetSolveStep()->SetActStep(fi);
        ptPDE_->GetSolveStep()->SetActFreq(std::abs(GetFrequency(fi)));
        ptPDE_->GetSolveStep()->GetEigenMode(fi);

        // stupid paraview needs an increasing series of save_value :(

        double save_value = -1.0;

        if(domain->GetOptimization() && domain->GetOptimization()->GetOptimizerType() != Optimization::EVALUATE_INITIAL_DESIGN)
        {
          // time is step.nr
          int total = eigenFreqs->GetSize() * wvs;
          int digs =  boost::lexical_cast<string>(total).size();
          double sig = std::pow((float) 10.0, -digs); // 1e-2 -> 10 ^ -2 -> a compiler complained with simply 10
          save_value = stepNum + (w * wvs + fi + 1) * sig; // +1 for one based

          LOG_DBG3(efd) << "SR total=" << total << " digs=" << digs << " sig=" << sig << " count=" << (w * wvs + fi + 1);
        }
        else // for bloch case we label <step>.<nr> from the info.xml
          save_value = isBloch_ ? w + (fi+1.0) / (eigenFreqs->GetSize() < 9 ? 10.0 : 100.0) : std::abs(GetFrequency(fi));

        LOG_DBG(efd) << "SR w=" << w << " fi=" << fi << " save_step_=" << save_step_ << " save_value=" << save_value;

        handler_->BeginStep(save_step_, save_value);
        ptPDE_->WriteResultsInFile(save_step_, save_value);
        handler_->FinishStep();

        if(writeAllSteps_ || isPartOfSequence_)
          simState_->WriteStep(save_step_, save_value);

        save_step_++;
      }
    }
  }

 /* void EigenFrequencyDriver::StoreResults(unsigned int stepNum, double step_val)
  {
    // stepNum and step_val are ignored
    LOG_DBG(efd) << "SR step=" << stepNum << " val=" << step_val;

    unsigned int wvs = isBloch_ ? wave_vectors.GetSize() : 1; // save wave vector size

    for(unsigned int w = 0; w < wvs; w++)
    {
      for(unsigned int fi=0; fi < eigenFreqs->GetSize(); fi++)
      {
        // Phase 2: calculate eigenmodes
        if(writeModes_)
        {
          ptPDE_->GetSolveStep()->SetActStep(fi);
          ptPDE_->GetSolveStep()->SetActFreq(std::abs(GetFrequency(fi)));
          ptPDE_->GetSolveStep()->GetEigenMode(fi);

          // stupid paraview needs an increasing series of save_value :(

          double save_value = -1.0;

          if(domain->GetOptimization())
          {
            // time is step.nr
            int total = eigenFreqs->GetSize() * wvs;
            int digs =  boost::lexical_cast<string>(total).size();
            double sig = std::pow(10, -digs); // 1e-2 -> 10 ^ -2
            save_value = stepNum + (w * wvs + fi + 1) * sig; // +1 for one based

            LOG_DBG3(efd) << "SR total=" << total << " digs=" << digs << " sig=" << sig << " count=" << (w * wvs + fi + 1);
          }
          else // for bloch case we label <step>.<nr> from the info.xml
            save_value = isBloch_ ? w + (fi+1.0) / 100.0 : std::abs(GetFrequency(fi));

          LOG_DBG(efd) << "SR w=" << w << " fi=" << fi << " save_step_=" << save_step_ << " save_value=" << save_value;

          handler_->BeginStep(save_step_, save_value);
          ptPDE_->WriteResultsInFile(save_step_, save_value);
          handler_->FinishStep();

          if(writeAllSteps_ || isPartOfSequence_)
            simState_->WriteStep(save_step_, save_value);

          save_step_++;
        }
      }
    }
  }*/

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

    if(isBloch_ && this->ibz_)
    {
      // repeat the first step at the and of bloch.plot for a full ibz for plotting, not when explicit wave vectors are given
      if(wave_vector_step == (int) wave_vectors.GetSize() - 1) {
        bloch_plot_ << first_plot_line_;
        bloch_plot_.close();
      }
    }
  }


} // end of namespace
