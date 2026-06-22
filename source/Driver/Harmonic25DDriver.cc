// ================================================================================================
/*!
 *       \file     Harmonic25DDriver.cc
 *       \brief    Driver for 2.5D frequency domain analysis
 *
 *       \date     November 2024
 *       \author   Likun Luo
 */
//================================================================================================
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>
#include <boost/lexical_cast.hpp>
#include <boost/smart_ptr/shared_ptr.hpp>
// signal handling for catching Ctr-C
#include <signal.h>

#include "Driver/Harmonic25DDriver.hh"
#include "Driver/SolveSteps/StdSolveStep.hh"
#include "Driver/Assemble.hh"
#include "OLAS/algsys/AlgebraicSys.hh"
#include "DataInOut/SimState.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ResultHandler.hh"
#include "DataInOut/ProgramOptions.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "Utils/Timer.hh"
#include "Utils/mathParser/mathParser.hh"
#include "PDE/StdPDE.hh"
#include "Domain/Domain.hh"

using std::cout;
using std::endl;

namespace CoupledField {
    DEFINE_LOG(harmonic25ddriver,"harmonic25dDriver");

    // ***************
    //   Constructor
    // ***************
    Harmonic25DDriver::Harmonic25DDriver( UInt sequenceStep, bool isPartOfSequence,
                                          shared_ptr<SimState> state, Domain* domain,
                                          PtrParamNode paramNode, PtrParamNode infoNode)
    : SingleDriver( sequenceStep, isPartOfSequence, state, domain, paramNode, infoNode ) {
        // Set analysistype
        analysis_ = BasePDE::HARMONIC25D;

        // Specifying parameter node
        param_ = param_->Get("harmonic25d");
        info_ = info_->Get("harmonic25d");
        info_->Get(ParamNode::HEADER)->Get("unit")->SetValue("Hz");

        // Initialize value
        baseFreq_ = param_->Get("excitationFreq")->MathParse<Double>();
        mathParser_->SetValue( MathParser::GLOB_HANDLER, "baseFreqHarmonic25D", baseFreq_);

        // Check if we should store the calculated wavenumber spectrum to results or not
        // param_->Has("storeSpectrum") ? storeSpectrum_ = param_->Get("storeSpectrum")->As<bool>() : storeSpectrum_ = false;

        /* if (!storeSpectrum_) {
          EXCEPTION("Initializing PDE, 2.5D Harmonic Driver inverse fourier transform not yet implemented!");
        } */
      }
    
    // ***************
    //   Destructor
    // ***************
    Harmonic25DDriver::~Harmonic25DDriver() {}
    
    // ***************
    //   Init
    // ***************
    void Harmonic25DDriver::Init(bool restart) {
      
      ReadWaveNumbers();
      // ReadEvalPosList();
      
      // copied from HarmonicDriver.cc
      // PtrParamNode in = info_->Get(ParamNode::HEADER);
      // in->Get("start")->SetValue(startFreq_);
      // in->Get("end")->SetValue(stopFreq_);
      // in->Get("numFreq")->SetValue(numFreq_);

      InitializePDEs();
    }

    // ****************
    //   SolveProblem
    // ****************
    void Harmonic25DDriver::SolveProblem() {
      //EXCEPTION("Solving Problem, 2.5D Harmonic Driver not yet implemented!");
      CalcWavenumberSpectrum();
      //CalcInverseTransform();
    }

    // Compute the whole wavenumber spectrum
    void Harmonic25DDriver::CalcWavenumberSpectrum() {
      // in harmonics one cannot extraxt the result writing to StoreResults() as
      // we have multiple frequencies. (exceptions is optimization)
      ptPDE_->WriteGeneralPDEdefines();
      handler_->BeginMultiSequenceStep( sequenceStep_, analysis_, numFreq_ );
      // storeSpectrum_ ? simState_->BeginMultiSequenceStep(sequenceStep_, analysis_) : void();

      UInt stpIndexStart = 0;
      // Perform one simulation for each desired wavenumber/frequency
      for ( UInt stpIndex = stpIndexStart; stpIndex < numFreq_; stpIndex++ ) {
        actFreqStep_ = waveNum_[stpIndex].step;
        ComputeFrequencyStep(waveNum_[stpIndex]);

        // Log info for this frequency - suppress in Optimization due to search steps
        //if(progOpts->IsQuiet())
        cout << ptPDE_->GetName() << ": 2.5D Harmonic step " << actFreqStep_ << " (" << stpIndex+1-stpIndexStart << "/" << numFreq_ << ")" << " frequency " << std::setprecision(2) << std::fixed << actFreq_ << "\t\r" << std::flush; 
        //else
        //  cout << endl << ptPDE_->GetName() << ": Harmonic step " << actFreqStep_ << " (" << stpIndex+1-stpIndexStart << "/" << numFreq_ << ")" <<" ======================= " << endl;

        analysis_id_.step = actFreqStep_;
        analysis_id_.freq = actFreq_;

        handler_->BeginStep( actFreqStep_, actFreq_);
        ptPDE_->WriteResultsInFile( actFreqStep_, actFreq_);
        handler_->FinishStep();

        // write to step result in hdf5 if storeSpectrum_ is true
        // storeSpectrum_ ? simState_->WriteStep(actFreqStep_, actFreq_) : void();

        // leave loop, if simulation should be aborted
        if ( abortSimulation_ ) {
          break;
        }
      }

      handler_->FinishMultiSequenceStep();
      // storeSpectrum_ ? simState_->FinishMultiSequenceStep(!abortSimulation_) : void();

      // Perform finalization only if not part of sequence
      if(!isPartOfSequence_) 
        handler_->Finalize();
    }

    // void Harmonic25DDriver::CalcInverseTransform() {
    //   EXCEPTION("Inverse Fourier Transform for 2.5D analysis not yet implemented!");
    // }

    // Compute one frequency step per wavenumber
    Double Harmonic25DDriver::ComputeFrequencyStep(Frequency const& freqStp) {
      assert(freqStp.step >= 1);
      assert(freqStp.step <= stopFreqStep_);

      actFreqStep_ = freqStp.step;
      actFreq_ = freqStp.freq;
      
      LOG_DBG(harmonic25ddriver) << "Step: " << actFreqStep_ << ", Frequency: " << actFreq_ << std::endl;

      this->analysis_id_.step = actFreqStep_;
      this->analysis_id_.time = actFreq_;

      // Set current frequency value in the mathParser
      mathParser_->SetValue( MathParser::GLOB_HANDLER, "f", actFreq_ );
      mathParser_->SetValue( MathParser::GLOB_HANDLER, "step", actFreqStep_ );

      // Perform steps for the solution
      ptPDE_->GetSolveStep()->SetActFreq( actFreq_ );
      ptPDE_->GetSolveStep()->SetActStep( actFreqStep_ );

      BaseSolveStep *step = ptPDE_->GetSolveStep();
      StdSolveStep *sstep = dynamic_cast<StdSolveStep*>(step);
      AlgebraicSys *algsys = sstep->GetAlgSys();
      Assemble *assemble = sstep->GetAssemble();
      StdPDE *pde = dynamic_cast<StdPDE*>(ptPDE_);

      // Assemble rhs and system matrices for this wavenumber step
      algsys->InitRHS();
      algsys->InitSol();
      pde->SetRhsValues();
      assemble->AssembleLinRHS();
      assemble->AssembleMatrices();
      pde->SetBCs();
      algsys->InitMatrix(SYSTEM);

      // 2.5D dynamic-stiffness prefactors: excitation omega vs. wavenumber (kz) omega
      Double BaseOmega    = 2*M_PI*baseFreq_;   // excitation angular frequency
      Double WaveNumOmega = 2*M_PI*actFreq_;    // wavenumber angular frequency

      // Mass Matrix prefactor (omega_kz^2 - omega^2)
      Complex FactorM = Complex(WaveNumOmega*WaveNumOmega - BaseOmega*BaseOmega, 0.0);
      // Damping Matrix prefactor (impedance BC): i*omega
      Complex FactorC = Complex(0.0, BaseOmega);
      // Absorbing-BC auxiliary prefactor, sign-dependent on (omega^2 - omega_kz^2)
      Complex SqrtRootResult = std::sqrt(-FactorM);
      Complex FactorABC = Complex(0.0, 0.0);
      if (std::abs(SqrtRootResult.real()) > 1.0E-12) {
        // propagating: i*sqrt(omega^2 - omega_kz^2)
        FactorABC = Complex(0.0,1.0) * SqrtRootResult;
      } else if (std::abs(SqrtRootResult.imag()) > 1.0E-12) {
        // evanescent: sqrt(omega_kz^2 - omega^2)
        FactorABC = Complex(std::abs(SqrtRootResult.imag()), 0.0);
      }

      // Effective matrix: 1*K + FactorC*C + FactorABC*C_abc + FactorM*M
      std::map<FEMatrixType,Complex> dynamicStiffnessMatrixFactors;
      dynamicStiffnessMatrixFactors.insert( std::pair<FEMatrixType,Complex>(STIFFNESS,   Complex(1.0,0.0)) );
      dynamicStiffnessMatrixFactors.insert( std::pair<FEMatrixType,Complex>(DAMPING,     FactorC) );
      dynamicStiffnessMatrixFactors.insert( std::pair<FEMatrixType,Complex>(DAMPING_AUX, FactorABC) );
      dynamicStiffnessMatrixFactors.insert( std::pair<FEMatrixType,Complex>(MASS,        FactorM) );
      algsys->ConstructEffectiveMatrix(NO_FCT_ID, dynamicStiffnessMatrixFactors);

      // Apply Dirichlet, (re)build precond/solver, solve, sync solution back
      algsys->BuildInDirichlet();
      algsys->SetupPrecond();
      algsys->SetupSolver();
      algsys->Solve();
      sstep->StoreSolutionToFeFunctions();

      return actFreq_;
    }

    void Harmonic25DDriver::ReadWaveNumbers() {
      // get upper bound and resolution of the wavenumber spectrum
      freqCutoff_ = param_->Get("cutoffFreq")->MathParse<Double>();
      freqResolution_ = param_->Get("freqResolution")->MathParse<Double>();

      if (param_->Has("startFreq")) {
        startFreq_ = param_->Get("startFreq")->MathParse<Double>();
      }
      
      assert(freqCutoff_ >= startFreq_ && "Cutoff Frequency must be greater than or equal Start Frequency!");

      // calculate how many steps we have to compute, adding 1 to include the Start Frequency
      numFreq_ = static_cast<int>(std::ceil((freqCutoff_-startFreq_)/ freqResolution_)) + 1;

      // Initialize the vector and set the first element to start Frequency 
      waveNum_.Resize(numFreq_);
      waveNum_[0].freq = startFreq_;
      waveNum_[0].step = 1; //1-based

      // calculate the list of wavenumbers
      for (unsigned int i = 1; i < numFreq_; i++) {
        waveNum_[i].freq = waveNum_[i-1].freq + freqResolution_;
        waveNum_[i].step = i + 1;
      }
      
      // get the actual stop frequency
      stopFreq_ = waveNum_.Last().freq;
      stopFreqStep_ = waveNum_.Last().step;
   
      // check if the last frequency step is the same as the cutoff frequency
      if (std::abs(stopFreq_ - freqCutoff_) > 1e-6) {
        WARN("cutoffFreq is not an exact multiple of freqRes. The cutoffFreq is changed to " << stopFreq_);
      }
    }

    // bool Harmonic25DDriver::ReadEvalPosList() {
    //   // check for existence
    //   if(!param_->Has("posList")) return false;

    //   ParamNodeList& list = param_->Get("posList")->GetChildren();
    //   // check empty list
    //   if(list.GetSize() == 0) EXCEPTION("Cannot have empty position list for IFT calculation");

    //   // read from list
    //   numPos_ = list.GetSize();
    //   iftPos_.Resize(numPos_);
    //   for(int fi = 0; fi < (int) numPos_; fi++)
    //   {
    //     PtrParamNode pn = list[fi];
    //     assert(pn->GetName() == "pos");
    //     Position& f = iftPos_[fi];
    //     f.step = fi+1;
    //     f.pos = pn->Get("value")->MathParse<Double>();

    //     // plausibility checking
    //     for(int o = 0; o < fi-1; o++)
    //       if(iftPos_[o].pos == f.pos)
    //         EXCEPTION("Multiple occurence of instances in posList: pos = " << f.pos << " at position " << (o+1) << " and " << (fi+1));
    //   }
    //   // when finished, return true
    //   return true;
    // }

    void Harmonic25DDriver::SetToStepValue(UInt stepNum, Double stepVal )  {
    // ensure that this method is only called if simState has input
    if( ! simState_->HasInput()) {
      EXCEPTION( "Can only set external time step, if simulation state "
              << "is read from external file" );
    }
    
    actFreqStep_ = stepNum;
    actFreq_ = stepVal;

    // Set current frequency value in the mathParser
    domain_->GetMathParser()->SetValue( MathParser::GLOB_HANDLER, "f", actFreq_ );
    domain_->GetMathParser()->SetValue( MathParser::GLOB_HANDLER, "step", actFreqStep_ );
  }
}




