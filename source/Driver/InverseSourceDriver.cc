#include <fstream>
#include <iostream>
#include <string>
#include <math.h>

#include <boost/lexical_cast.hpp>

// signal handling for catching Ctr-C
#include <signal.h>

#include "Driver/InverseSourceDriver.hh"
#include "Driver/SolveSteps/StdSolveStep.hh"
#include "Driver/Assemble.hh"
#include "DataInOut/SimState.hh"
#include "Utils/Timer.hh"
#include "Utils/mathParser/mathParser.hh"

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ResultHandler.hh"
#include "DataInOut/ProgramOptions.hh"

#include "PDE/StdPDE.hh"
#include "FeBasis/BaseFE.hh"
#include "Domain/Domain.hh"


using std::cout;
using std::endl;


// Define pointer to driver instance, needed for the signal handler
// to communicate with
//InverseSourceDriver * instance = NULL;


namespace CoupledField
{
  // ***************
  //   Constructor
  // ***************
  InverseSourceDriver::InverseSourceDriver( UInt sequenceStep, bool isPartOfSequence,
                                  shared_ptr<SimState> state, Domain* domain,
                                  PtrParamNode paramNode, PtrParamNode infoNode)
    : SingleDriver( sequenceStep, isPartOfSequence, state, domain,
                    paramNode, infoNode ), 
      timer_(new Timer())
  {
    // Set correct analysistype
	//set analysis_ to HARMONIC to be able to visualize the results with paraview!!
    analysis_ = BasePDE::HARMONIC; //INVERSESOURCE;

    param_ = param_->Get("inverseSource");

    // replace our info node by a more detailed level
    info_ = info_->Get("inverseSource");
    info_->Get(ParamNode::HEADER)->Get("unit")->SetValue("Hz");

    fileNameMeasdata_ = param_->Get("measDataFilename")->As<std::string>();
    freq_ = param_->Get( "freq")->MathParse<Double>();

    //regularization parameters
    alpha_ = param_->Get( "alphaPar")->MathParse<Double>();
    beta_ = param_->Get( "betaPar")->MathParse<Double>();
    rho_ = param_->Get( "rhoPar")->MathParse<Double>();
    qExp_  = param_->Get( "expPar")->MathParse<Double>();
    resStopCritRel_ = param_->Get( "resStopCritRel")->MathParse<Double>();
    maxReduceParSteps_   = param_->Get( "maxReduceParSteps")->MathParse<UInt>();
    maxGradSteps_   = param_->Get( "maxGradSteps")->MathParse<UInt>();
    maxNumLineSearch_  = param_->Get( "maxNumLineSearch")->MathParse<UInt>();
    minMicDistant_ = param_->Get( "minMicDistant")->MathParse<Double>();
    logLevel_ =  param_->Get("logLevel")->As<std::string>();
    if ( param_->Has("scale2Val") ) {
    	scale2Val_ = param_->Get( "scale2Val")->MathParse<Double>();
    }
    else
    	scale2Val_ = 1.0;

    //check for correct logLevel
    if ( logLevel_ != "1" && logLevel_ != "2" &&logLevel_ != "3" ) {
    	std::cerr << "Log level has to be 1, 2 or 3; set it now to 1";
    	logLevel_ = "1";
    }

    //! we fix to use delta functions for the source identification
    approxSourceWithDeltaFnc_ = true;

    numFreq_ = 1;
    actFreqStep_ = 0;
    actFreq_ = 1.0;
    stopFreqStep_ = 1.0;
    restartStep_ = 0;
    timePerStep_ = 0;
    measL2squared_ = 0;
    isRestarted_ = false;
    
    // Check for presence of restart flag.
    writeRestart_ = true;
    PtrParamNode restartNode = param_->Get("writeRestart", ParamNode::PASS);
    if (restartNode)
      writeRestart_ = restartNode->As<bool>();
    
    // read flag if all results should get written to database file section
    // to allow e.g. for general postprocessing or result extraction
    param_->GetValue("allowPostProc", writeAllSteps_, ParamNode::PASS );
    
  }

  InverseSourceDriver::~InverseSourceDriver()
  { 
    if( !simState_->HasInput() ) {
      // unregister signal handler and use default action
      // register signal handler
      if( signal( SIGINT, SIG_DFL) == SIG_ERR ) {
        std::cerr << "Could not assign default signal action" << std::endl; // no exceptions in destructors!
        domain->GetInfoRoot()->ToFile();
        exit(-1);
      }

      // set global pointer to zero
      //instance = NULL;
    }
  }


  void InverseSourceDriver::Init(bool restart)
  {
    isRestarted_ = restart;
    
//    PtrParamNode in = info_->Get(ParamNode::PN_HEADER);
//    in->Get("freq")->SetValue(freq_);
//    std::cout << "Read freq: " << freq_ << std::endl;

    actFreq_ = freq_;

    // Initialize first multisequence step, as the method "CheckStoreResults" 
    // relies on the result handler to know already about the current
    // sequencestep. However, in case of optimization, the sequence step
    // gets initialized in Optimization::SolveProblem()
    InitializePDEs();
  }


  // ****************
  //   SolveProblem
  // ****************
  void InverseSourceDriver::SolveProblem()
  {
    // in harmonics one cannot extraxt the result writing to StoreResults() as
    // we have multiple frequencies. (exceptions is optimization)

    ptPDE_->WriteGeneralPDEdefines();
    handler_->BeginMultiSequenceStep( sequenceStep_, analysis_, numFreq_ );
    
    if(writeRestart_ || writeAllSteps_ )
      simState_->BeginMultiSequenceStep( sequenceStep_, analysis_ );
    
    // Read restart information
    ReadRestart();
    numFreq_ = numFreq_ - restartStep_;
    stopFreqStep_ = numFreq_ + restartStep_;
    
    // Perform one simulation for each desired frequency
    for ( actFreqStep_ = restartStep_+1; actFreqStep_ <= numFreq_+restartStep_; actFreqStep_++ )
    {
      // register signal handler only, if it is a child driver
      if( actFreqStep_ > restartStep_+1 && !simState_->HasInput() ) {
        if( signal( SIGINT, InverseSourceDriver::SignalHandler) == SIG_ERR ) {
          EXCEPTION( "Could not register Signal Handler");
        }
      }
      
      // Determine next frequency value
      ComputeFrequencyStep(actFreqStep_);

      // Log info for this frequency - suppress in Optimization due to search steps
      if(progOpts->IsQuiet())
        cout << ptPDE_->GetName() << ": InverseSource step " << actFreqStep_ << " frequency " << actFreq_ << endl; 
      else
        cout << endl << ptPDE_->GetName() << ": InverseSource step " 
        << actFreqStep_ <<" ======================= " << endl;

      handler_->BeginStep( actFreqStep_, actFreq_ );
      ptPDE_->WriteResultsInFile( actFreqStep_, actFreq_ );
      handler_->FinishStep( );

      // write out re-start in case of aborted simulation or if all steps should be written
      if(  abortSimulation_  || writeAllSteps_ ) {
        if( writeRestart_ || writeAllSteps_ )
         simState_->WriteStep( actFreqStep_, actFreq_);
      }
      
      // leave loop, if simulation should be aborted
      if ( abortSimulation_ ) {
        break;
      }
        
      // perform runtime estimation
      Double totalTime = timer_->GetWallTime();
      timePerStep_ = totalTime / (Double) actFreqStep_;
      Double remainingTime = (numFreq_ - actFreqStep_) * timePerStep_;
      auto time = std::chrono::system_clock::now();
      time += std::chrono::seconds(static_cast<long int>(remainingTime));
      //analysis_id_->Get("timePerStep")->SetValue( timePerStep_ );
      PtrParamNode envNode = info_->GetRoot()->Get(ParamNode::HEADER)->Get("environment");
      envNode->Get("estimatedEnd")->SetValue(Timer::TimeStamp(time));
      envNode->Get("remainingTime")->SetValue(remainingTime);
    } // loop: frequencies

    handler_->FinishMultiSequenceStep();
    if(writeRestart_ || writeAllSteps_ )
      simState_->FinishMultiSequenceStep( !abortSimulation_ );

    // Perform finalization only if not part of sequence
    if(!isPartOfSequence_) 
      handler_->Finalize();
  }

  Double InverseSourceDriver::ComputeFrequencyStep(UInt actFreqStep)
  {
	assert(actFreqStep >= 1);
	assert(actFreqStep <= numFreq_+restartStep_);

	actFreqStep_ = actFreqStep;

	this->analysis_id_.step = actFreqStep;
	this->analysis_id_.time = actFreq_;

	//get pointers to CoefFncs
    bool isRHSsource = false;
    bool isRHSmeas = false;
    UInt num = 0;
    std::map<SolutionType, shared_ptr<BaseFeFunction> > rhsFeFunctions;
    rhsFeFunctions = ptPDE_->GetRhsFeFunctions();
    std::map<SolutionType, shared_ptr<BaseFeFunction> >::iterator rFncIt= rhsFeFunctions.begin();
    while ( rFncIt != rhsFeFunctions.end() ) {
     	LoadCoefList loadCoefs = rFncIt->second->GetLoadCoefFunctions();
     	LoadCoefList::iterator it = loadCoefs.begin();
     	for ( ; it != loadCoefs.end(); ++it  ) {
     		PtrCoefFct ptCoef = it->first;
     		CoefFunction::CoefInverseType type = ptCoef->GetInverseType();
     		if ( type == CoefFunction::INVSOURCE ) {
     			rhsSource_ = ptCoef;
   		    if ( approxSourceWithDeltaFnc_ ) {
   		    	rhsSource_->SetInverseSourceApproxType(CoefFunction::DELTA);
   		    }
     			isRHSsource = true;
     			num++;
     		}
     		else if ( type == CoefFunction::INVMEASURE ) {
     			rhsMeas_ = ptCoef;
     			isRHSmeas = true;
     			num++;
     		}
     	}
     	rFncIt++;
    }
    // check
    if ( num != 2 || isRHSsource == false || isRHSmeas == false ) {
    	std::string message = "You have to define two 'rhsValues', one for the searched";
    	message += " for sources and one for the measurements";
    	EXCEPTION(message);

    }

    //set the regularization parameters
    rhsMeas_->SetInverseParam(alpha_, beta_, rho_, qExp_, actFreq_,fileNameMeasdata_,
    		                  logLevel_, scale2Val_);
    rhsSource_->SetInverseParam(alpha_, beta_, rho_, qExp_, actFreq_,fileNameMeasdata_,
    		                    logLevel_, scale2Val_);

    // Set current frequency value in the mathParser
    mathParser_->SetValue( MathParser::GLOB_HANDLER, "f", actFreq_ );
    mathParser_->SetValue( MathParser::GLOB_HANDLER, "step", actFreqStep_ );

    bool gradMethod = true;

    if ( logLevel_ == "2" || logLevel_ == "3" )
    	std::cout << "\nInverse Source Localization: Start\n"
		          << "=================================================="
                  << "==========================\n";

    if ( gradMethod ) {
    	Double res2, funcVal, sigma;
    	Double stepLength=1.0;
    	Double optAmp, optPhase;

    	//solve adj: just for initiallization
    	rhsSource_->SetActive(false);
    	rhsMeas_->SetActive(true);

    	ptPDE_->GetSolveStep()->SetActFreq( actFreq_ );
    	ptPDE_->GetSolveStep()->SetActStep( actFreqStep_ );
    	ptPDE_->GetSolveStep()->PreStepHarmonic();
    	ptPDE_->GetSolveStep()->SolveStepHarmonic();
    	ptPDE_->GetSolveStep()->PostStepHarmonic();

    	//compute with RHS as the current identified sources
    	rhsMeas_->SetActive(false);
    	rhsSource_->SetActive(true);

    	ptPDE_->GetSolveStep()->SetActFreq( actFreq_ );
    	ptPDE_->GetSolveStep()->SetActStep( actFreqStep_ );
    	ptPDE_->GetSolveStep()->PreStepHarmonic();
    	ptPDE_->GetSolveStep()->SolveStepHarmonic();
    	ptPDE_->GetSolveStep()->PostStepHarmonic();

    	//save source data
    	rhsSource_->UpdateSource( stepLength, false );

    	//compute residual squared
    	rhsMeas_->ComputeDiff2Meas( res2 );

    	//compute square of L2-norm of measured pressure at mic-positions
    	rhsMeas_->ComputeMeasL2squared( measL2squared_ );

    	//compute Tikhonov functional
    	rhsSource_->ComputeTikh(funcVal,res2);

    	if ( logLevel_ == "3" )
    		std::cout << "\n Starting absolute L2 error: " << sqrt(res2) << std::endl;

    	Double kappa = 0.5;
    	UInt outerIter = 0;
    	while ( std::sqrt( res2/measL2squared_ ) > resStopCritRel_  &&  outerIter < maxReduceParSteps_ ) {
    		//solve adj
    		if ( logLevel_ == "2" || logLevel_ == "3" )
    			std::cout << "\n Outer step " << outerIter+1 << " using: \n alpha: " << alpha_
				          << ";  beta: " << beta_ << ";  rho: " << rho_ <<  ";  Exponent q: "
						  << qExp_ <<"\n\n";

    		rhsSource_->SetActive(false);
    		rhsMeas_->SetActive(true);
    		if ( !approxSourceWithDeltaFnc_ )
    			ptPDE_->GetSolveStep()->SetAdjointSource();

    		ptPDE_->GetSolveStep()->SetActFreq( actFreq_ );
    		ptPDE_->GetSolveStep()->SetActStep( actFreqStep_ );
    		ptPDE_->GetSolveStep()->PreStepHarmonic();
    		ptPDE_->GetSolveStep()->SolveStepHarmonic();
    		ptPDE_->GetSolveStep()->PostStepHarmonic();

    		//compute optimality condition and delta values of amplitude and phase
    		rhsSource_->ComputeOptCondition(optAmp, optPhase);

    		UInt innerIter=0;
    		Double funcValNew;
    		Double normGrad2 = optAmp + optPhase;
    		Double redNormGrad2 = normGrad2*1.0E-05;

    		if ( logLevel_ == "3" )
    			std::cout << " Outer Step: " << outerIter+1 << ";  normGrad2:"
				          << normGrad2 << std::endl;

    		//do the gradient steps
    		while ( normGrad2 >  redNormGrad2 && innerIter < maxGradSteps_ ) {
    			stepLength = 1.0;
        		sigma      = 1.0e-4;

    			//make an increment on source data
    			rhsSource_->UpdateSource( stepLength, true );

    			//compute with RHS as the current identified sources
    			rhsMeas_->SetActive(false);
    			rhsSource_->SetActive(true);
    			ptPDE_->GetSolveStep()->SetActFreq( actFreq_ );
    			ptPDE_->GetSolveStep()->SetActStep( actFreqStep_ );
    			ptPDE_->GetSolveStep()->PreStepHarmonic();
    			ptPDE_->GetSolveStep()->SolveStepHarmonic();
    			ptPDE_->GetSolveStep()->PostStepHarmonic();

    			rhsMeas_->ComputeDiff2Meas( res2 );
    			rhsSource_->ComputeTikh( funcValNew, res2 );

    			if ( logLevel_ == "3" )
    				std::cout << " Gradien step: " << innerIter+1 << "; Tikhonov value: "
					          << funcValNew << std::endl;

    			UInt iterLineSearch = 0;
    			if ( logLevel_ == "3" )
    				std::cout << "\n Start line search " << std::endl;

    			stepLength = 1.0;
    			kappa = 0.25;
    			while ( funcValNew > funcVal - sigma*stepLength*normGrad2
    					&& iterLineSearch < maxNumLineSearch_ ) {

    				stepLength *= kappa;
    				rhsSource_->UpdateSource( stepLength, true );

    				//compute with RHS as the current identified sources
    				rhsMeas_->SetActive(false);
    				rhsSource_->SetActive(true);
    				ptPDE_->GetSolveStep()->SetActFreq( actFreq_ );
    				ptPDE_->GetSolveStep()->SetActStep( actFreqStep_ );
    				ptPDE_->GetSolveStep()->PreStepHarmonic();
    				ptPDE_->GetSolveStep()->SolveStepHarmonic();
    				ptPDE_->GetSolveStep()->PostStepHarmonic();

    				rhsMeas_->ComputeDiff2Meas( res2 );
    				if ( logLevel_ == "3" )
    					std::cout << " LineSearch: " <<  iterLineSearch+1 << ";  StepLenghth: " << stepLength
						          << "; Relative L2 error: " << std::sqrt(res2 / measL2squared_)*100.0
								  << "%\n"<< std::endl;

    				rhsSource_->ComputeTikh(funcValNew,res2);
    				if ( logLevel_ == "3" )
    					std::cout << std::endl;

    				iterLineSearch++;
    			}

    			if ( logLevel_ == "2" || logLevel_ == "3" )
    				std::cout << " Gradient steps " << innerIter+1
					          << ";  L2 error: " << std::sqrt(res2 / measL2squared_) *100.0
    						  << "%;  Used line search steps: " << iterLineSearch
    						  << std::endl;

    			//save source data
    			rhsSource_->UpdateSource( stepLength, false );

    			funcVal = funcValNew;

    			//solve adj
    			rhsSource_->SetActive(false);
    			rhsMeas_->SetActive(true);
    			if ( !approxSourceWithDeltaFnc_ )
    			  ptPDE_->GetSolveStep()->SetAdjointSource();

    			ptPDE_->GetSolveStep()->SetActFreq( actFreq_ );
    			ptPDE_->GetSolveStep()->SetActStep( actFreqStep_ );
    			ptPDE_->GetSolveStep()->PreStepHarmonic();
    			ptPDE_->GetSolveStep()->SolveStepHarmonic();
    			ptPDE_->GetSolveStep()->PostStepHarmonic();

    			//compute optimality condition and delta values of amplitude and phase
    			rhsSource_->ComputeOptCondition(optAmp, optPhase);
    			if ( logLevel_ == "3" )
    				std::cout << " Gradient step: " << innerIter+1 << ";  OptCondAmp: "
    				          << optAmp << "  OptCondPhase: " << optPhase << std::endl;

    			normGrad2 = optAmp + optPhase;
    			innerIter++;
    		}

    		alpha_ *= 0.5;
    		beta_  *= 0.5;
    		rho_   *= 0.5;
    		rhsSource_->ChangeInverseParam(alpha_, beta_, rho_);
    		rhsMeas_->ChangeInverseParam(alpha_, beta_, rho_);
    		outerIter++;
    	}

    	//final computation
    	rhsMeas_->SetActive(false);
    	rhsSource_->SetActive(true);
    	//scale back to the correct values fitting to the measurements
    	rhsSource_->UpdateSource( stepLength, false, true );
    	ptPDE_->GetSolveStep()->SetActFreq( actFreq_ );
    	ptPDE_->GetSolveStep()->SetActStep( actFreqStep_ );
    	ptPDE_->GetSolveStep()->PreStepHarmonic();
    	ptPDE_->GetSolveStep()->SolveStepHarmonic();
    	ptPDE_->GetSolveStep()->PostStepHarmonic();

    	std::cout << "\nInverse Source Localization: Final Result\n" << "=================================================="
    			  << "==========================\n\n"
    			  << "Relative L2 error: " << std::sqrt(res2 / measL2squared_)*100.0 << "%\n"
				  << "Number of reducing parameter steps: " << outerIter << "\n"
				  << "Final regularization parameters:\n  alpha: " << alpha_*2.0
				  << ";  beta: " << beta_*2.0 << ";  rho: " << rho_*2.0 << ";  Exponent q: " << qExp_ << "\n"
				  << "\n=================================================="
				  << "==========================\n\n";
    }
    else {
    	Double optAmp, optPhase;
    	UInt iter = 0;
    	bool notConverged = true;

    	while ( notConverged && iter < maxReduceParSteps_ ) {
    		//compute with RHS being the conjugate difference of computed acoustic pressure
    		//and measured pressure in the microphone points
    		rhsSource_->SetActive(false);
    		rhsMeas_->SetActive(true);
    		if ( !approxSourceWithDeltaFnc_ )
    			ptPDE_->GetSolveStep()->SetAdjointSource();

    		ptPDE_->GetSolveStep()->SetActFreq( actFreq_ );
    		ptPDE_->GetSolveStep()->SetActStep( actFreqStep_ );
    		ptPDE_->GetSolveStep()->PreStepHarmonic();
    		std::cout << "\n SOLVE_MEASURE \n" << std::endl;
    		ptPDE_->GetSolveStep()->SolveStepHarmonic();
    		std::cout << "\n INV_MEASURE solved \n" << std::endl;
    		ptPDE_->GetSolveStep()->PostStepHarmonic();

    		rhsSource_->ComputeOptCondition(optAmp, optPhase);

    		if ( optAmp < 1e-6 && optPhase < 1e-6 )
    			notConverged = false;

    		std::cout << "\n OptCond, Amp: " << optAmp << "   Phase: " << optPhase << std::endl;

    		//compute with RHS as the current identified sources
    		rhsMeas_->SetActive(false);
    		rhsSource_->SetActive(true);
    		ptPDE_->GetSolveStep()->SetActFreq( actFreq_ );
    		ptPDE_->GetSolveStep()->SetActStep( actFreqStep_ );
    		ptPDE_->GetSolveStep()->PreStepHarmonic();
    		std::cout << "\n SOLVE_SOURCE \n" << std::endl;
    		ptPDE_->GetSolveStep()->SolveStepHarmonic();
    		ptPDE_->GetSolveStep()->PostStepHarmonic();
    		std::cout << "\n INV_SOURCE solved \n" << std::endl;

    		Double error;
    		rhsMeas_->ComputeDiff2Meas( error);
    		std::cout << "\n L2 Error: " << error << std::endl;

    		rhsSource_->SetActive(false);

    		iter++;
    	}
    }


    return actFreq_;
  }


  unsigned int InverseSourceDriver::StoreResults(UInt stepNum, double step_val)
  {
    assert(analysis_ == BasePDE::INVERSESOURCE);

    // Write results into output-file(s)
    handler_->BeginStep(stepNum, step_val);
    ptPDE_->WriteResultsInFile(stepNum, step_val);
    handler_->FinishStep( );

    return stepNum;
  }



  void InverseSourceDriver::ReadRestart() {

    if ( isRestarted_ ) {

      // Ensure simState is present
      assert( simState_ );

      // Create input reader from current output reader
      bool hasInput = simState_->SetInputReaderToSameOutput();
      if( !hasInput)  {
        EXCEPTION( "Can not perform restarted simulation, as HDF5 file "
            << "contains no restart information.");
      }

      if( simState_->IsCompleted( sequenceStep_ )) {
        std::cout << "\n\n";
        std::cout << "*******************************************************\n";
        std::cout << " No restart necessary, as the desired number of \n";
        std::cout << " frequency steps are already computed. \n";
        std::cout << "*******************************************************\n\n";
        restartStep_ = stopFreqStep_ +1; 
        return;

      } else{

        // Obtain last step
        UInt lastStepNum;
        Double lastStepVal;
        simState_->GetLastStepNum(sequenceStep_, lastStepNum, lastStepVal );
        restartStep_ = lastStepNum;

        // if lastStep is 0, no restart possibility
        if( lastStepNum == 0 ) {
          EXCEPTION( "Can not perform restarted simulation, as HDF5 file "
              << "contains no restart information.");
        }
        std::cout << "\n\n";
        std::cout << "*******************************************************\n";
        std::cout << " Continuing simulation from step " << restartStep_  << std::endl;
        std::cout << "*******************************************************\n";
      }
    }

  }


  void InverseSourceDriver::SetToStepValue(UInt stepNum, Double stepVal )  {
    // ensure that this method is only called if simState has input
    if( ! simState_->HasInput()) {
      EXCEPTION( "Can only set external time step, if simulation state "
              << "is read from external file" );
    }
    
    actFreqStep_ = stepNum;;
    actFreq_ = stepVal;

    // Set current frequency value in the mathParser
    domain_->GetMathParser()->SetValue( MathParser::GLOB_HANDLER, "f", actFreq_ );
    domain_->GetMathParser()->SetValue( MathParser::GLOB_HANDLER, "step", actFreqStep_ );
  }
  
  
  void InverseSourceDriver::SignalHandler( int sig ) {

//    if( !instance->abortSimulation_) {
//
//      instance->abortSimulation_ = true;
//
//      std::cout << "\n\n";
//      std::cout << "*******************************************************\n";
//      std::cout << " Simulation will be halted after the current frequency\n";
//      std::cout << " step " << instance->actFreqStep_ << " / " << instance->actFreq_
//                << " Hz.\n\n";
//      std::cout << " Estimated time before end of simulation run: " <<
//          int(instance->timePerStep_) << " s" << std::endl;
//      std::cout << "*******************************************************\n\n";
//
//    }

  }
}
