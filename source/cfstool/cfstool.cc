// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <def_cfs_stats.hh>

#include <iostream>
#include <math.h>

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ProgramOptions.hh"
#include "DataInOut/SimInput.hh"
#include "DataInOut/SimOutput.hh"
#include "Domain/Mesh/GridCFS/GridCFS.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "DataInOut/ColoredConsole.hh"
#include "Forms/Operators/BaseBOperator.hh"

#include "ParamsInit.hh"
#include "HelperFuncs.hh"
#include "WVT.hh"

using namespace CoupledField;

//! Global parameter node and info node instance
PtrParamNode param;
PtrParamNode info;

namespace CFSTool {

  void PrintWarning(Exception& ex) {
    
    // Print warning on command line
    std::string msg = ex.GetMsg();
    std::string fileName = ex.GetFileName();
    unsigned int lineNum = ex.GetLineNum();
    
    std::cerr << "\n "
              << fg_blue << "WARNING:" << fg_reset << "\n "
              << msg << std::endl;
    std::cerr << "\n(" << fileName << ", Line " 
              << lineNum  << ")\n\n";
  }

  void Convert( const std::string& inFile, const std::string& outFile ) {

    // obtain input reader for inFile
    shared_ptr<SimInput> input = GetInputReader( inFile, param, info );

    // read in mesh
    input->InitModule();
    UInt dim = input->GetDim();
    Grid * ptGrid = new GridCFS(dim, param, info);
    input->ReadMesh(ptGrid);
    ptGrid->FinishInit();

    // obtain output writer
    shared_ptr<SimOutput> output = GetOutputWriter( outFile, param, info );

    // obtain number of multiSequenceSteps and get analysis types
    std::map<UInt, BasePDE::AnalysisType> types;
    std::map<UInt, UInt> numSteps;
    input->GetNumMultiSequenceSteps( types, numSteps );
    std::cout << "\nFound " << types.size() << " sequence step(s)\n";

    // check if the input reader has results
    bool printGridOnly = false;
    if( types.size() == 0) {
      printGridOnly = true;
      std::cerr << "Printing only grid\n";
    }

     output->Init( ptGrid, printGridOnly);

     // only iterate over results, if not only the mesh is converted
     if( !printGridOnly ) {


       // iterate over all multiSequenceSteps
       std::map<UInt,UInt>::iterator it;
       for( it = numSteps.begin(); it != numSteps.end(); it++ ) {

         UInt actMsStep = it->first;
         std::cout << "\n----------------------------\n"
                   << " Converting sequence step " << actMsStep << std::endl
                   << "----------------------------\n\n";
         // get resulttypes
         StdVector<shared_ptr<ResultInfo> > infos;
         input->GetResultTypes( actMsStep, infos );

         StdVector<shared_ptr<BaseResult> > results;
         std::map<UInt, Double> stepVals;
         std::map<shared_ptr<ResultInfo>, std::map<UInt, Double> > resultSteps;

         if( infos.GetSize() > 0 ){
             std::cout << "Converting the following results:\n";
         }
         // iterate over all result types
         for( UInt iRes = 0; iRes < infos.GetSize(); iRes++) {

           shared_ptr<ResultInfo> actRes = infos[iRes];
           std::cout << "\t" << actRes->resultName << std::endl;

           // get stepvalues
           input->GetStepValues( actMsStep, actRes, resultSteps[actRes] );
           stepVals.insert( resultSteps[actRes].begin(),
                            resultSteps[actRes].end() );

           // iterate over all regions
           StdVector<shared_ptr<EntityList> > resEntities;
           input->GetResultEntities( actMsStep, actRes, resEntities );
           for( UInt iEntity = 0; iEntity < resEntities.GetSize(); iEntity++ ) {
             // generate new result object and add it to output writer
             shared_ptr<BaseResult > result;
             if( types[actMsStep] == BasePDE::STATIC || 
                 types[actMsStep] == BasePDE::TRANSIENT ) {
               result = shared_ptr<BaseResult>( new Result<Double>() );
             } else {
               result = shared_ptr<BaseResult>( new Result<Complex>() );
             }
             result->SetEntityList( resEntities[iEntity] );
             result->SetResultInfo( actRes );
             results.Push_back( result );
             // Note: as the real values of saveBegin, saveInc and saveEnd are almost
             // nevert queried within an output format. we simply set saveBegin = 1,
             // saveInc = 1, saveEnd = number of result steps.
             output->RegisterResult( result, 1, 1, resultSteps[actRes].size(), false );
           }
         }

         // notify writer
         output->BeginMultiSequenceStep( actMsStep, types[actMsStep], numSteps[actMsStep] );

         // iterate over all stepvalues of this multisequence step
         for( UInt iStep = 0; iStep < numSteps[actMsStep]; iStep++ ) {

           // check, if current step contains any results
           if( stepVals.find( iStep+1) == stepVals.end() )
             continue;

           UInt actStepNum = iStep+1;
           Double actStepVal = stepVals[iStep+1];

           output->BeginStep( actStepNum, actStepVal );

           // iterate over all results
           for( UInt iRes = 0; iRes < results.GetSize(); iRes++) {

             // check if current result is defined within this step
             if( resultSteps[results[iRes]->GetResultInfo()].find(actStepNum)
                 == resultSteps[results[iRes]->GetResultInfo()].end() ) {
               continue;
             }

             try {
               input->GetResult( actMsStep, actStepNum, results[iRes] );
               output->AddResult( results[iRes] );
             } catch (Exception& ex ) {
             std::cerr <<  "\nResult '" << results[iRes]->GetResultInfo()->resultName
                       << "' in MsStep" << actMsStep << ", step " << actStepNum
                       << " could not be converted:\n\n";
             std::cerr << ex.what() << std::endl;
             }
           }
           output->FinishStep();
         }
         output->FinishMultiSequenceStep();

       } // loop over multisequence steps
     } // printGridOnly
     output->Finalize();
     delete ptGrid;
     std::cout << "\nOutput successfully written to " << outFile << std::endl;
  } //Convert()

  Double Diff( const std::string& inFile_ref,
               const std::string& inFile_fut,
               const std::string& outFile,
               bool normedtomax,
               bool isHistory,
               std::string& maxDiffResultName) {

    // obtain input reader for inFiles
    StdVector< shared_ptr<SimInput> > inputs(2);
    inputs[0] = GetInputReader( inFile_ref, param, info );
    inputs[1] = GetInputReader( inFile_fut, param, info );

    // check capabilities of input class
    StdVector<bool> readerCaps(2);
    StdVector<std::string> readerDescriptions(2);
    readerDescriptions[0] = "reference results";
    readerDescriptions[1] = "file under test";
    StdVector<std::string> readerCapsResults(2);
    for(UInt i=0; i<inputs.GetSize(); i++) 
    {
      readerCaps[i] = CheckReaderCapabilities(readerDescriptions[i],
                                              inputs[i],
                                              readerCapsResults[i]);
    }

    if( !readerCaps[0] || !readerCaps[1] ) {
      std::cerr << "Some input files are only capable of handling meshes, not results!\n\n";

      for(UInt i=0; i<inputs.GetSize(); i++) 
      {
        std::cerr << readerCapsResults[i] << std::endl;
      }

      exit(EXIT_FAILURE);
    }

    // obtain input reader for inFiles
    shared_ptr<SimInput>& input_ref = inputs[0];
    shared_ptr<SimInput>& input_fut = inputs[1];

    // read in mesh of reference file
    input_ref->InitModule();
    UInt dim = input_ref->GetDim();
    Grid * ptGrid_ref = new GridCFS(dim, param, info);
    input_ref->ReadMesh(ptGrid_ref);
    ptGrid_ref->FinishInit();

    // Obtain settings for node and element limit check
    bool checkLimits = false;
    if (param->Has("checkLimits"))
    {
      checkLimits = param->Get("checkLimits")->As<bool>();         
    }
    else
    {
      if(inFile_ref.find( ".h5ref") != std::string::npos) 
      {
        checkLimits = true;
      }         
    }

    // Do the actual limits check
    if(checkLimits) 
    {
      UInt numNodes = ptGrid_ref->GetNumNodes();
      UInt numElems = ptGrid_ref->GetNumElems();
      std::string msg =
          "\nIn the interests of fast test execution and the ability to postprocess\n" \
          "using certain postprocessors with limitations (e.g. GiD without a license),\n" \
          "please keep the node and element count below the limits.\n";

      if(numNodes > SOFT_NODE_LIMIT) 
      {
        WARN("Number of nodes " << numNodes << 
             " exceeds soft limit " << SOFT_NODE_LIMIT << ".");
      }

      if(numElems > SOFT_ELEM_LIMIT) 
      {
        WARN("Number of elements " << numElems << 
             " exceeds soft limit " << SOFT_ELEM_LIMIT << ".");
      }

      if(numNodes > HARD_NODE_LIMIT) 
      {
        delete ptGrid_ref;
        EXCEPTION("Number of nodes " << numNodes << 
                  " exceeds hard limit " << HARD_NODE_LIMIT << "." <<
                  msg);
      }

      if(numElems > HARD_ELEM_LIMIT) 
      {
        delete ptGrid_ref;
        EXCEPTION("Number of elements " << numElems << 
                  " exceeds hard limit " << HARD_ELEM_LIMIT << "." << 
                  msg);
      }

    }


    // read in mesh of file under test
    input_fut->InitModule();
    Grid * ptGrid_fut = new GridCFS(dim, param, info);
    input_fut->ReadMesh(ptGrid_fut);
    ptGrid_fut->FinishInit();

    // obtain output writer
    bool printGridOnly = false;
    shared_ptr<SimOutput> output;
    if( outFile != "" ) {
      output = GetOutputWriter( outFile, param, info );
      output->Init( ptGrid_fut, printGridOnly);
    }

    // obtain number of Sequence Steps and get analysis types
    std::map<UInt, BasePDE::AnalysisType> types;
    std::map<UInt, UInt> numSteps;
    input_fut->GetNumMultiSequenceSteps( types, numSteps, isHistory );

    std::cout << "\nFound " << types.size() << " sequence step(s) in '" << inFile_fut << "'\n";
    std::map<UInt, BasePDE::AnalysisType> types_ref;
    std::map<UInt, UInt> numSteps_ref;
    input_ref->GetNumMultiSequenceSteps( types_ref, numSteps_ref, isHistory );
    std::cout << "\nFound " << types_ref.size() << " sequence step(s) in '" << inFile_ref << "'\n";

    if(types.size() != types_ref.size()){
      std::cout << "Reference file '" << inFile_ref << "' and file under test '" << inFile_fut
          << "' have different number of sequence steps!\n";
      exit(EXIT_FAILURE);
    }


    // iterate over all Sequence Steps
    Double maxDiff = 0.0;
    std::map<UInt,UInt>::iterator it;
    for( it = numSteps_ref.begin(); it != numSteps_ref.end(); it++ ) {

      UInt actMsStep = it->first;
      std::cout << " Diffing sequence step " << actMsStep << std::endl
          << "-------------------------\n\n";

      // get resulttypes
      StdVector<shared_ptr<ResultInfo> > infos, infos_ref;
      input_ref->GetResultTypes( actMsStep, infos_ref, isHistory );
      input_fut->GetResultTypes( actMsStep, infos, isHistory );

      StdVector<shared_ptr<BaseResult> > inResults_fut, inResults_ref, outResults;
      // stepnumbers, for which at least one result is defined
      std::map<UInt, Double> stepVals, stepVals_ref;
      // contains the stepnumbers/-values in which the particular result is
      // defined in
      std::map<shared_ptr<ResultInfo>, std::map<UInt, Double> > resultSteps;
      std::map<shared_ptr<ResultInfo>, std::map<UInt, Double> > resultSteps_ref;

      if( infos.GetSize() > 0 ){
        std::cout << "Performing diff on the following results:\n";
      }
      // iterate over all result types of input_ref
      for( UInt iRes=0, numRes=infos_ref.GetSize(); iRes < numRes; ++iRes) {

        std::cout << "\t" << infos_ref[iRes]->resultName << "\n";

        // find the corresponding result in the file to be tested
        Integer testRes = -1;
        for ( UInt j=0, n=infos.GetSize(); j<n; ++j ) {
          if ( infos_ref[iRes]->resultType == NO_SOLUTION_TYPE ) {
            if ( infos_ref[iRes]->resultName == infos[j]->resultName ) {
              testRes = j;
              break;
            }
          } else if ( infos_ref[iRes]->resultType == infos[j]->resultType ) {
            testRes = j;
            break;
          }
        }
        if ( testRes == -1 ) {
          EXCEPTION("Result '" << infos_ref[iRes]->resultName
                    << "' is missing in file '" << inFile_fut << "'.");
        }
        
        // get stepvalues of reference file
        shared_ptr<ResultInfo> actRes_ref = infos_ref[iRes];
        input_ref->GetStepValues( actMsStep, actRes_ref,
                                  resultSteps_ref[actRes_ref], isHistory);
        stepVals_ref.insert( resultSteps_ref[actRes_ref].begin(),
                             resultSteps_ref[actRes_ref].end() );

        // get stepvalues of file under test
        shared_ptr<ResultInfo> actRes = infos[testRes];
        input_fut->GetStepValues( actMsStep, actRes,
                                  resultSteps[actRes], isHistory);
        stepVals.insert( resultSteps[actRes].begin(),
                         resultSteps[actRes].end() );

        // Loop over all step values in both sets and compare them. Thus we can see
        // differences e.g. in eigenfrequency analysis
        std::map<UInt, Double>::const_iterator svIt_fut, svIt_ref;
        svIt_fut = stepVals.begin();
        svIt_ref = stepVals_ref.begin();
        for( ; svIt_ref != stepVals_ref.end(); ++svIt_ref, ++svIt_fut ) {
          if( svIt_fut->first != svIt_ref->first ) {
            EXCEPTION( "Encountered different result steps for result " << 
                       infos[testRes]->resultName );
          } else {

            Double val_fut = svIt_fut->second;
            Double val_ref = svIt_ref->second;
            Double relDiff = std::abs(std::abs(val_fut-val_ref))/std::abs(val_fut);
            if ( relDiff > 1e-4 ) {
              EXCEPTION("Time / Frequency values of step " << svIt_ref->first << " differ by " 
                        << relDiff*100.0 << " %:\n"
                        << "\treference value: " << val_ref <<" s / Hz\n\tcompared value:  " << val_fut
                        << " s / Hz" << std::endl );
            }
          }
        }


        // iterate over all regions
        StdVector<shared_ptr<EntityList> > regions;
        input_ref->GetResultEntities( actMsStep, infos_ref[iRes],
                                      regions, isHistory );
        for( UInt iRegion = 0; iRegion < regions.GetSize(); iRegion++ ) {
          // generate new result object and add it to output writer
          shared_ptr<BaseResult > inResult_fut, inResult_ref, outResult;
          if( types[actMsStep] == BasePDE::TRANSIENT ||
              types[actMsStep] == BasePDE::STATIC ) {
            inResult_fut = shared_ptr<BaseResult>( new Result<Double>() );
            inResult_ref = shared_ptr<BaseResult>( new Result<Double>() );
            outResult = shared_ptr<BaseResult>( new Result<Double>() );
          } else {
            inResult_fut = shared_ptr<BaseResult>( new Result<Complex>() );
            inResult_ref = shared_ptr<BaseResult>( new Result<Complex>() );
            outResult = shared_ptr<BaseResult>( new Result<Complex>() );
          }
          inResult_fut->SetEntityList( regions[iRegion] );
          inResult_ref->SetEntityList( regions[iRegion] );
          outResult->SetEntityList( regions[iRegion] );

          inResult_fut->SetResultInfo( infos[testRes] );

          inResult_ref->SetResultInfo( infos_ref[iRes] );
          outResult->SetResultInfo( infos[testRes] );

          inResults_fut.Push_back( inResult_fut );
          inResults_ref.Push_back( inResult_ref );
          outResults.Push_back( outResult );
          if( output) {
            // Hardcoded: set output format to AMPL_PHASE
            //outResult->GetResultInfo()->complexFormat = AMPLITUDE_PHASE;
            outResult->GetResultInfo()->complexFormat = REAL_IMAG;

            // CAUTION: begin, inc, end are hardcoded and noch checked for each result
            output->RegisterResult( outResult, 1, 1,
                                    resultSteps[actRes].size(),
                                    isHistory );
          }
        }
      }

      Vector<Double> maxResVec_ref;
      maxResVec_ref.Resize( inResults_ref.GetSize() );

      // For transient simulation find maximum amplitude over all timesteps
      if( types[actMsStep] == BasePDE::STATIC ||
          types[actMsStep] == BasePDE::TRANSIENT ) {

        // iterate over all results
        for( UInt iRes = 0; iRes < inResults_ref.GetSize(); iRes++) {

          if(numSteps[actMsStep] != numSteps_ref[actMsStep]){
            std::cout << "Reference file '" << inFile_ref << "' has " << numSteps_ref[actMsStep] << " and '" << inFile_fut
                << "' has " << numSteps[actMsStep] << " time steps!\n";
            exit(EXIT_FAILURE);
          }

          maxResVec_ref[iRes] = 0.0;
          // iterate over all time steps
          for( UInt iStep = 0; iStep < numSteps_ref[actMsStep]; iStep++ ) {

            UInt actStepNum = iStep+1;
            // check if current result is defined within this step
            if( resultSteps[inResults_ref[iRes]->GetResultInfo()].find(actStepNum)
                == resultSteps[inResults_ref[iRes]->GetResultInfo()].end() ) {
              continue;
            }

            input_ref->GetResult( actMsStep, actStepNum, inResults_ref[iRes], isHistory );
            Vector<Double> & inVec_ref =
                dynamic_cast<Result<Double>& >(*inResults_ref[iRes]).GetVector();

            for( UInt i = 0; i<inVec_ref.GetSize(); i++ ) {
              if( std::abs(inVec_ref[i]) > maxResVec_ref[iRes] )
                maxResVec_ref[iRes] = std::abs(inVec_ref[i]);
            }
          }
          std::cout << "For result '" << inResults_ref[iRes]->GetResultInfo()->resultName
              << "' maximum amplitude is: " << maxResVec_ref[iRes] << "\n";
        }
      }


      // notify writer
      if( output) {
        output->BeginMultiSequenceStep( actMsStep, types_ref[actMsStep],
                                        numSteps_ref[actMsStep] );
      }

      // iterate over all time/frequency steps
      for( UInt iStep = 0; iStep < numSteps_ref[actMsStep]; iStep++ ) {

        // check, if current step contains any results
        if( stepVals.find( iStep+1) == stepVals.end() )
          continue;
        UInt actStepNum = iStep+1;
        Double actStepVal = stepVals[iStep+1];

        if( output) {
          output->BeginStep( actStepNum, actStepVal );
        }
        std::cout << "\n\t======================================================\n";
        std::cout << "\t  Treating step nr " << actStepNum << " val " << actStepVal << " s / Hz/ iteration\n";
        std::cout << "\t======================================================\n";

        // iterate over all results
        for( UInt iRes = 0; iRes < inResults_fut.GetSize(); iRes++) {
          // check if current result is defined within this step
          if( resultSteps[inResults_fut[iRes]->GetResultInfo()].find(actStepNum)
              == resultSteps[inResults_fut[iRes]->GetResultInfo()].end() ) {
            continue;
          }

          // obtain both result objects for current step
          input_fut->GetResult( actMsStep, actStepNum, inResults_fut[iRes], isHistory );
          input_ref->GetResult( actMsStep, actStepNum, inResults_ref[iRes], isHistory );

          std::cout << "\t-- Comparing result " <<
              inResults_fut[iRes]->GetResultInfo()->resultName << " on " 
              << inResults_fut[iRes]->GetEntityList()->GetName() << " --\n";

          // get number of dofs of result
          UInt numDofs = inResults_fut[iRes]->GetResultInfo()->dofNames.GetSize();

          // cast result objects, get vector and calculate difference vector
          if( types[actMsStep] == BasePDE::STATIC || 
              types[actMsStep] == BasePDE::TRANSIENT ) {
            Vector<Double> & inVec_fut =
                dynamic_cast<Result<Double>& >(*inResults_fut[iRes]).GetVector();
            Vector<Double> & inVec_ref =
                dynamic_cast<Result<Double>& >(*inResults_ref[iRes]).GetVector();
            Vector<Double> & outVec =
                dynamic_cast<Result<Double>& >(*outResults[iRes]).GetVector();
            outVec.Resize( inVec_fut.GetSize() );

            // find maximum amplitude of inResult_ref
            for( UInt i = 0; i<inVec_ref.GetSize(); i++ ) {
              if( std::abs(inVec_ref[i]) > maxResVec_ref[iRes])
                maxResVec_ref[iRes] = std::abs(inVec_ref[i]);
            }

            // calculate difference entrywise
            outVec = inVec_fut - inVec_ref;
            if (normedtomax == true) {
              outVec /= maxResVec_ref[iRes];
            }

            // find maximum entry in difference vector
            for( UInt i = 0; i < outVec.GetSize(); i++ ) {
              if( std::abs(outVec[i]) > maxDiff) {
                maxDiff = std::abs(outVec[i]) ;
                maxDiffResultName = inResults_fut[iRes]->GetResultInfo()->resultName;
              }
            }

          } else if( types[actMsStep] == BasePDE::EIGENFREQUENCY ||  types[actMsStep] == BasePDE::EIGENVALUE ) {
            // in the eigenfrequency case we only compare the absolute value,
            // as the sign is not defined uniquely
            Vector<Complex> & inVec_fut =
                dynamic_cast<Result<Complex>& >(*inResults_fut[iRes]).GetVector();
            Vector<Complex> & inVec_ref =
                dynamic_cast<Result<Complex>& >(*inResults_ref[iRes]).GetVector();
            Vector<Complex> & outVec =
                dynamic_cast<Result<Complex>& >(*outResults[iRes]).GetVector();
            outVec.Resize( inVec_fut.GetSize() );

            // find maximum amplitude of inResult_ref in every frequency step
            Double maxRes_ref = 0.0;
            for( UInt i = 0; i<inVec_ref.GetSize(); i++ ) {
              if( std::abs(inVec_ref[i]) > maxRes_ref)
                maxRes_ref = std::abs(inVec_ref[i]);
            }

            Double aDiff, pDiff;
            // iterate over all dofs
            for (UInt dof = 0; dof<numDofs ; dof++) {
              // iterate over number of entities
              for( UInt i = 0, nEntities=UInt(inVec_ref.GetSize()/numDofs); i<nEntities; i++ ) {

                // index to access entity 'i' of dof 'dof'
                UInt actIndex = i * numDofs + dof;

                // amplitude difference
                if (normedtomax == true)
                  aDiff = std::abs(( std::abs(inVec_fut[actIndex]) - std::abs(inVec_ref[actIndex]) )/maxRes_ref);
                else 
                  aDiff = std::abs(std::abs(inVec_fut[actIndex]) - std::abs(inVec_ref[actIndex]));


                // phase difference in multiples of pi (just used for
                // grid diffs
                pDiff = RadPhase(inVec_fut[actIndex]) - RadPhase(inVec_ref[actIndex]);

                // correct 2*pi-offset if phase angles have different signs
                if ( (std::abs(pDiff)>M_PI) && (pDiff<0) )
                  pDiff+= 2*M_PI;
                if ( (std::abs(pDiff)>M_PI) && (pDiff>0) )
                  pDiff-= 2*M_PI;

                // Dirty hack! Write differences in real_imag format.
                outVec[actIndex] = Complex( aDiff, pDiff*180/M_PI );

                // return maximum of amplitude difference
                if  (aDiff > maxDiff)
                  maxDiff = aDiff;
              }


            } // loop dof
            std::cout << "\n\tMaximum overall rel. difference = " << maxDiff << "\n\n";

          }  else if( types[actMsStep] == BasePDE::HARMONIC ) { 
            Vector<Complex> & inVec_fut =
                dynamic_cast<Result<Complex>& >(*inResults_fut[iRes]).GetVector();
            Vector<Complex> & inVec_ref =
                dynamic_cast<Result<Complex>& >(*inResults_ref[iRes]).GetVector();
            Vector<Complex> & outVec =
                dynamic_cast<Result<Complex>& >(*outResults[iRes]).GetVector();
            outVec.Resize( inVec_fut.GetSize() );

            // find maximum amplitude of inResult_ref in every frequency step
            Double maxRes_ref = 0.0;
            for( UInt i = 0; i<inVec_ref.GetSize(); i++ ) {
              if( std::abs(inVec_ref[i]) > maxRes_ref)
                maxRes_ref = std::abs(inVec_ref[i]);
            }

            Double aDiff, pDiff, aMax=0.0, aMin=0.0, pMax=0.0, pMin=0.0;
            Double rDiff, iDiff, rMax=0.0, iMax=0.0;

            // iterate over all dofs
            for (UInt dof = 0; dof<numDofs ; dof++) {
              // iterate over number of entities
              for( UInt i = 0, nEntities=UInt(inVec_ref.GetSize()/numDofs); i<nEntities; i++ ) {

                // index to access entity 'i' of dof 'dof'
                UInt actIndex = i * numDofs + dof;

                // amplitude difference
                if (normedtomax == true)
                  aDiff = ( std::abs(inVec_fut[actIndex]) - std::abs(inVec_ref[actIndex]) )/maxRes_ref;
                else
                  aDiff = std::abs(inVec_fut[actIndex]) - std::abs(inVec_ref[actIndex]);

                // phase difference in multiples of pi
                pDiff = RadPhase(inVec_fut[actIndex]) - RadPhase(inVec_ref[actIndex]);

                // correct 2*pi-offset if phase angles have different signs
                if ( (std::abs(pDiff)>M_PI) && (pDiff<0) )
                  pDiff+= 2*M_PI;
                if ( (std::abs(pDiff)>M_PI) && (pDiff>0) )
                  pDiff-= 2*M_PI;

                // Dirty hack! Write differences in real_imag format.
                outVec[actIndex] = Complex( aDiff, pDiff*180/M_PI );

                // maximum and minimum values
                if( pDiff > pMax )
                  pMax = pDiff;
                if( pDiff < pMin )
                  pMin = pDiff;
                if( aDiff > aMax )
                  aMax = aDiff;
                if( aDiff < aMin )
                  aMin = aDiff;

                // maximum difference in real and imaginary part
                rDiff = std::abs( inVec_fut[actIndex].real() - inVec_ref[actIndex].real() )/maxRes_ref;
                iDiff = std::abs( inVec_fut[actIndex].imag() - inVec_ref[actIndex].imag() )/maxRes_ref;
                if ( rDiff > rMax )
                  rMax = rDiff;
                if ( iDiff > iMax)
                  iMax = iDiff;
              }
              if( normedtomax == true)
                std::cout << "\n\tMaximum rel. + amplitude difference:  " << aMax*100 << " %\n"
                << "\tMaximum rel. - amplitude difference: " << aMin*100 << " %\n";
              else
                std::cout << "\n\tMaximum + amplitude difference:  " << aMax <<  "\n"
                << "\tMaximum - amplitude difference: " << aMin <<  "\n";

              std::cout << "\tMaximum + phase difference:      " << pMax*180/M_PI <<  " deg\n"
                  << "\tMaximum - phase difference:     " << pMin*180/M_PI <<  " deg\n";

              // return maxDiff for differences in real and imaginary part
              if ( (rMax > iMax) && (rMax > maxDiff) )
                maxDiff = rMax;
              else if ( (iMax > rMax) && (iMax > maxDiff) )
                maxDiff = iMax;
            } // loop dof
            std::cout << "\n\tMaximum overall rel. difference = " << maxDiff << "\n\n";

          } // switch: Analysitype
          // add result to output file
          if (output )
            output->AddResult( outResults[iRes] );
        } // loop over results
        if( output )
          output->FinishStep();
      }
      if( output )
        output->FinishMultiSequenceStep();

    }
    if( output ) {
      output->Finalize();
      std::cout << "\nOutput successfully written to " << outFile << std::endl;
    }
    delete ptGrid_fut;
    delete ptGrid_ref;

    return maxDiff;

  } //Diff

  /**
   * calcAverage calclulates the average of a physical field at each timestep. It
   * stores the results in outFile
   * @param inFile The file which carries the data (e.g. acoustic pressure , mechanical
   * displacement)
   * @param meshFile The mesh file to associate each value with a point in space
   * Import not all points should be included into the averagin.
   * @return void
   */

  void calcAverage( const std::string& inFile,
              const std::string& outFile)
  {
    std::string isFreeCoord = param->Get("freeCoord")->As<std::string>();
    if (isFreeCoord != "")
    {
      std::cerr << "selection of region not implemented!\n";
      exit(EXIT_FAILURE);
      // TODO: this is the call then for setting everything up
      //SetFreeCoord();
    }

    // obtain input reader for inFiles
    shared_ptr<SimInput> input = GetInputReader(inFile, param, info);

    // check capabilities of input class
    bool printGridOnly = false;
    if (std::find( input->GetCapabilities().begin(),
          input->GetCapabilities().end(),
          SimInput::MESH_RESULTS )
        == input->GetCapabilities().end())
    {
      std::cerr << "input file is only capable of handling mesh, not results!\n";
      exit(EXIT_FAILURE);
    }

    // read in mesh of input
    input->InitModule();
    UInt dim = input->GetDim();
    Grid * ptGrid = new GridCFS(dim, param, info);
    input->ReadMesh(ptGrid);
    ptGrid->FinishInit();

    // obtain output writer
    shared_ptr<SimOutput> output;
    output = GetOutputWriter( outFile, param, info );
    output->Init( ptGrid, printGridOnly);

    // obtain number of Sequence Steps and get analysis types
    std::map<UInt, BasePDE::AnalysisType> types;
    std::map<UInt, UInt> numSteps;
    input->GetNumMultiSequenceSteps( types, numSteps, false );

    std::cout << "\nFound " << types.size() << " sequence step(s)\n";

    // iterate over all Sequence Steps
    std::map<UInt,UInt>::iterator it;
    for (it = numSteps.begin(); it != numSteps.end(); it++)
    {
      UInt actMsStep = it->first;
      std::cout << " averaging step " << actMsStep << std::endl
        << "-------------------------\n\n";

      // get resulttypes
      StdVector<shared_ptr<ResultInfo> > infos;
      input->GetResultTypes( actMsStep, infos, false );

      StdVector<shared_ptr<BaseResult> > inResults, outResults;
      // stepnumbers, for which at least one result is defined
      std::map<UInt, Double> stepVals;
      // contains the stepnumbers/-values in which the particular result is
      // defined in
      std::map<shared_ptr<ResultInfo>, std::map<UInt, Double> > resultSteps;

      if (infos.GetSize() > 0)
      {
        std::cout << "Performing average on the following results:\n";
      }
      // iterate over all result types of input
      for (UInt iRes = 0; iRes < infos.GetSize(); iRes++)
      {
        std::cout << "\t" << infos[iRes]->resultName << "\n\n";

        // get stepvalues
        shared_ptr<ResultInfo> actRes = infos[iRes];
        input->GetStepValues( actMsStep, actRes,
            resultSteps[actRes], false);
        stepVals.insert( resultSteps[actRes].begin(),
            resultSteps[actRes].end() );

        // iterate over all regions
        StdVector<shared_ptr<EntityList> > regions;
        input->GetResultEntities( actMsStep, infos[iRes],
            regions, false );
        for (UInt iRegion = 0; iRegion < regions.GetSize(); iRegion++)
        {
          // generate new result object and add it to output writer
          shared_ptr<BaseResult > inResult, outResult;
          if (types[actMsStep] != BasePDE::HARMONIC)
          {
            inResult  = shared_ptr<BaseResult>( new Result<Double>() );
            outResult = shared_ptr<BaseResult>( new Result<Double>() );
          } else {
            std::cerr << "Averaging over harmonic results does not make sense\n";
            exit(EXIT_FAILURE);
          }
          inResult->SetEntityList( regions[iRegion] );
          outResult->SetEntityList( regions[iRegion] );

          inResult->SetResultInfo( infos[iRes] );
          outResult->SetResultInfo( infos[iRes] );

          inResults.Push_back( inResult );
          outResults.Push_back( outResult );
          if (output)
          {
            // Hardcoded: set output format to AMPL_PHASE
            //outResult->GetResultInfo()->complexFormat = AMPLITUDE_PHASE;
            outResult->GetResultInfo()->complexFormat = REAL_IMAG;

            // CAUTION: begin, inc, end are hardcoded and noch checked for each result
            output->RegisterResult( outResult, 1, 1,
                                    resultSteps[actRes].size(),
                                    false );
          }
        }
      }

      // notify writer
      if (output)
      {
        output->BeginMultiSequenceStep( actMsStep, types[actMsStep],
            numSteps[actMsStep] );
      }

      // iterate over all time/frequency steps
      for (UInt iStep = 1; iStep <= numSteps[actMsStep]; ++iStep)
      {
        // check, if current step contains any results
        if (stepVals.find(iStep) == stepVals.end())
        {
          continue;
        }
        Double actStepVal = stepVals[iStep];

        if (output)
        {
          output->BeginStep(iStep, actStepVal);
        }
        // iterate over all results
        for (UInt iRes = 0; iRes < inResults.GetSize(); iRes++)
        {
          // check if current result is defined within this step
          if (resultSteps[inResults[iRes]->GetResultInfo()].find(iStep)
              == resultSteps[inResults[iRes]->GetResultInfo()].end())
          {
            continue;
          }

          // obtain result objects for current step
          input->GetResult(actMsStep, iStep, inResults[iRes], false);

          // get number of dofs of result
#if 0 // TODO check if needed
          UInt numDofs = inResults[iRes]->GetResultInfo()->dofNames.GetSize();
#endif

          // cast result objects, get vector and calculate difference vector
          if (types[actMsStep] != BasePDE::HARMONIC)
          {
            Double meanVal = 0.0;
            Vector<Double> & inVec =
              dynamic_cast<Result<Double>& >(*inResults[iRes]).GetVector();
            Vector<Double> & outVec =
              dynamic_cast<Result<Double>& >(*outResults[iRes]).GetVector();
            UInt inVec_size = inVec.GetSize();

            // sum up and divide by number of entries <- averaging
            for (UInt i = 0; i < inVec_size; ++i)
            {
              meanVal += inVec[i];
            }
            meanVal /= (Double)inVec_size;
            outVec.Resize(inVec_size);
            for (UInt i = 0; i < inVec_size; ++i)
            {
              outVec[i] = meanVal;
            }
          } else {
            std::cerr << "Averaging over harmonic results does not make sense\n";
            exit(EXIT_FAILURE);
          }
          // add result to output file
          if (output)
          {
              output->AddResult(outResults[iRes]);
          }
        }
        if (output)
        {
            output->FinishStep();
        }
      }
      if (output)
      {
          output->FinishMultiSequenceStep();
      }

    }
    if (output)
    {
      output->Finalize();
      std::cout << "\nOutput successfully written to " << outFile << std::endl;
    }
    delete ptGrid;
  } //calcAverage

  

  Double CheckL2( const std::string& inFile_ref, const std::string& inFile_fut, bool isHistory, Double& maxL2, Double& maxL2rel, std::string& maxDiffInfo, Double tolerance = 0.0) {
      // initialize the diff values
      maxL2 = 0.0;
      maxL2rel = 0.0;
      std::string infoStringRel;
      // obtain input reader for inFiles
      StdVector< shared_ptr<SimInput> > inputs(2);
      inputs[0] = GetInputReader( inFile_ref, param, info );
      inputs[1] = GetInputReader( inFile_fut, param, info );

      // check capabilities of input class
      StdVector<bool> readerCaps(2);
      StdVector<std::string> readerDescriptions(2);
      readerDescriptions[0] = "reference results";
      readerDescriptions[1] = "file under test";
      StdVector<std::string> readerCapsResults(2);
      for(UInt i=0; i<inputs.GetSize(); i++)
      {
          readerCaps[i] = CheckReaderCapabilities(readerDescriptions[i], inputs[i], readerCapsResults[i]);
      }

      if( !readerCaps[0] || !readerCaps[1] ) {
          std::cerr << "Some input files are only capable of handling meshes, not results!\n\n";

          for(UInt i=0; i<inputs.GetSize(); i++)
          {
              std::cerr << readerCapsResults[i] << std::endl;
          }

          exit(EXIT_FAILURE);
      }

      // obtain input reader for inFiles
      shared_ptr<SimInput>& input_ref = inputs[0];
      shared_ptr<SimInput>& input_fut = inputs[1];

      // read in mesh of reference file
      input_ref->InitModule();
      UInt dim = input_ref->GetDim();
      Grid * ptGrid_ref = new GridCFS(dim, param, info);
      input_ref->ReadMesh(ptGrid_ref);
      ptGrid_ref->FinishInit();

      // read in mesh of file under test
      input_fut->InitModule();
      Grid * ptGrid_fut = new GridCFS(dim, param, info);
      input_fut->ReadMesh(ptGrid_fut);
      ptGrid_fut->FinishInit();

      // obtain number of Sequence Steps and get analysis types
      std::map<UInt, BasePDE::AnalysisType> types;
      std::map<UInt, UInt> numSteps;
      input_fut->GetNumMultiSequenceSteps( types, numSteps, isHistory );

      std::cout << "Found " << types.size() << " sequence step(s) in '" << inFile_fut << "'\n";
      std::map<UInt, BasePDE::AnalysisType> types_ref;
      std::map<UInt, UInt> numSteps_ref;
      input_ref->GetNumMultiSequenceSteps( types_ref, numSteps_ref, isHistory );
      std::cout << "Found " << types_ref.size() << " sequence step(s) in '" << inFile_ref << "'\n";

      if(types.size() != types_ref.size()){
          std::cout << "Reference file '" << inFile_ref << "' and file under test '" << inFile_fut
                  << "' have different number of sequence steps!\n";
          exit(EXIT_FAILURE);
      }

      // iterate over all Sequence Steps
      auto checkL2Result = [tolerance](auto& inVec_fut, auto& inVec_ref,
          UInt actMsStep, UInt actStepNum, Double actStepVal,
          const std::string& resultName, const std::string& entityName,
          UInt numDofs, const char* kind,
          Double& maxL2_step, Double& maxL2rel_step, std::string& infoStringRel_step)
      {
          using T = std::remove_reference_t<decltype(inVec_fut[0])>;
          Vector<T> diffVec;
          diffVec.Resize(inVec_fut.GetSize());
          diffVec = inVec_fut - inVec_ref;
          T futMin, futMax, refMin, refMax;
          inVec_fut.MinMax(futMin, futMax);
          inVec_ref.MinMax(refMin, refMax);
          Double diffL2 = diffVec.NormL2();
          Double refL2  = inVec_ref.NormL2();
          Double testL2 = inVec_fut.NormL2();
          std::cout << "\n/Mesh/MultiStep_" << actMsStep << "/Step_" << actStepNum
                    << "/" << resultName << "/" << entityName
                    << " step_val=" << actStepVal
                    << ", " << (inVec_fut.GetSize() / numDofs) << "x" << numDofs << " entries (" << kind << ")\n";
          if(inVec_fut.ContainsNaN() || inVec_ref.ContainsNaN()) 
          {
            if(inVec_fut.ContainsNaN() && !inVec_ref.ContainsNaN()) {
              std::cout << "  ERROR: NaN value encountered in test but not reference result\n";
              diffL2 = std::numeric_limits<Double>::max();
            } 
            else  {
              std::cout << "  WARNING: NaN in reference=" <<  inVec_ref.ContainsNaN() << " and in test=" << inVec_fut.ContainsNaN() << " -> cannot compare\n";
              diffL2 = 0;
            }
          } 
          double relL2 = 0.0;
          if(refL2 == 0.0) {
              if(diffL2 == 0.0) {
                  std::cout << "  L2-Norm: test=" << testL2 << " ref=" << refL2 << " -> OK\n";
              } else {
                  std::cout << "  L2-Norm: test=" << testL2 << " ref=" << refL2 << " -> ERROR\n";
                  diffL2 = std::numeric_limits<Double>::max();

              }
          } else {
          Double relL2 = diffL2 / refL2;
          std::cout << "  L2-Norm: test=" << testL2 << " ref=" << refL2 << " rel_diff=" << relL2;
          if (relL2 <= tolerance) {
              std::cout << " <= " << tolerance << " -> OK\n";
          } else {
              std::cout << " > " << tolerance << " -> ERROR\n";
              std::cout << "  Minimum: test=" << futMin << " ref=" << refMin << "\n";
              std::cout << "  Maximum: test=" << futMax << " ref=" << refMax << "\n";

              int worst_idx = 0;
              double worst_rel_pt = 0.0;
              for (UInt i = 0; i < diffVec.GetSize(); i++) {
                  double rel = std::abs(diffVec[i]) / (std::abs(inVec_ref[i]) + 1e-300);
                  if (rel > worst_rel_pt) { worst_rel_pt = rel; worst_idx = (int)i; }
              }

              std::cout << "  entry with largest relative difference:"
                        << " dataset_index=" << worst_idx
                        << " test=" << inVec_fut[worst_idx]
                        << " ref=" << inVec_ref[worst_idx] << "\n";
          }
          if (diffL2 > maxL2_step) maxL2_step = diffL2;
          if (relL2 >= maxL2rel_step) {
              maxL2rel_step = relL2;
              infoStringRel_step = resultName + " on " + entityName;
          }
          } // closes else (refL2 != 0.0)
      };

      std::map<UInt,UInt>::iterator it;
      for( it = numSteps_ref.begin(); it != numSteps_ref.end(); it++ ) {
          Double maxL2_sequence = 0; // maximum L2 difference in sequence step
          Double maxL2rel_sequence = 0;
          std::string infoStringRel_sequence;
          UInt actMsStep = it->first;

          // get resulttypes
          StdVector<shared_ptr<ResultInfo> > infos, infos_ref;
          input_ref->GetResultTypes( actMsStep, infos_ref, isHistory );
          input_fut->GetResultTypes( actMsStep, infos, isHistory );

          StdVector<shared_ptr<BaseResult> > inResults_fut, inResults_ref, outResults;
          // stepnumbers, for which at least one result is defined
          std::map<UInt, Double> stepVals, stepVals_ref;
          // contains the stepnumbers/-values in which the particular result is defined in
          std::map<shared_ptr<ResultInfo>, std::map<UInt, Double> > resultSteps;
          std::map<shared_ptr<ResultInfo>, std::map<UInt, Double> > resultSteps_ref;

          std::cout << std::endl;
          std::cout << "Sequence step " << actMsStep << "\n---------------\n";
          if(infos.GetSize() > 0) 
            std::cout << "Checking the following results:\n";
         
          // iterate over all result types of input_ref
          for( UInt iRes=0, numRes=infos_ref.GetSize(); iRes < numRes; ++iRes) {
              std::cout << infos_ref[iRes]->resultName;
              // find the corresponding result in the file to be tested
              Integer testRes = -1;
              for ( UInt j=0, n=infos.GetSize(); j<n; ++j ) {
                  if ( infos_ref[iRes]->resultType == NO_SOLUTION_TYPE ) {
                      if ( infos_ref[iRes]->resultName == infos[j]->resultName ) {
                          testRes = j;
                          break;
                      }
                  } else if ( infos_ref[iRes]->resultType == infos[j]->resultType ) {
                      testRes = j;
                      break;
                  }
              }
              if ( testRes == -1 ) {
                  EXCEPTION("Result '" << infos_ref[iRes]->resultName << "' is missing in file '" << inFile_fut << "'.");
              }

              // get stepvalues of reference file
              shared_ptr<ResultInfo> actRes_ref = infos_ref[iRes];
              input_ref->GetStepValues( actMsStep, actRes_ref, resultSteps_ref[actRes_ref], isHistory);
              stepVals_ref.insert( resultSteps_ref[actRes_ref].begin(),  resultSteps_ref[actRes_ref].end() );

              // get stepvalues of file under test
              shared_ptr<ResultInfo> actRes = infos[testRes];
              input_fut->GetStepValues( actMsStep, actRes, resultSteps[actRes], isHistory);
              stepVals.insert( resultSteps[actRes].begin(),  resultSteps[actRes].end() );

              // Loop over all step values in both sets and compare them. Thus we can see
              // differences e.g. in eigenfrequency analysis
              std::map<UInt, Double>::const_iterator svIt_fut, svIt_ref;
              svIt_fut = stepVals.begin();
              svIt_ref = stepVals_ref.begin();
              for( ; svIt_ref != stepVals_ref.end(); ++svIt_ref, ++svIt_fut ) {
                  if( svIt_fut->first != svIt_ref->first ) {
                      EXCEPTION( "Encountered different result steps for result " << infos[testRes]->resultName );
                  } else {
                      Double val_fut = svIt_fut->second;
                      Double val_ref = svIt_ref->second;
                      Double relDiff = std::abs(std::abs(val_fut-val_ref))/std::abs(val_fut);
                      if ( relDiff > 1e-4 ) {
                          EXCEPTION("For step " << svIt_ref->first << " the step values differs by " << relDiff*100.0 << "%: "
                                  << "ref=" << val_ref << " test=" << val_fut);
                      }
                  }
              }


              // iterate over all regions
              StdVector<shared_ptr<EntityList> > regions;
              input_ref->GetResultEntities( actMsStep, infos_ref[iRes], regions, isHistory);
              std::cout << " (on regions: ";
              for( UInt iRegion = 0; iRegion < regions.GetSize(); iRegion++ ) {
                  // generate new result object and add it to output writer
                  shared_ptr<BaseResult > inResult_fut, inResult_ref;
                  if( types[actMsStep] == BasePDE::TRANSIENT || types[actMsStep] == BasePDE::STATIC ) {
                      inResult_fut = shared_ptr<BaseResult>( new Result<Double>() );
                      inResult_ref = shared_ptr<BaseResult>( new Result<Double>() );
                      //outResult = shared_ptr<BaseResult>( new Result<Double>() );
                  } else {
                      inResult_fut = shared_ptr<BaseResult>( new Result<Complex>() );
                      inResult_ref = shared_ptr<BaseResult>( new Result<Complex>() );
                      //outResult = shared_ptr<BaseResult>( new Result<Complex>() );
                  }
                  inResult_fut->SetEntityList( regions[iRegion] );
                  inResult_ref->SetEntityList( regions[iRegion] );
                  //outResult->SetEntityList( regions[iRegion] );

                  inResult_fut->SetResultInfo( infos[testRes] );

                  inResult_ref->SetResultInfo( infos_ref[iRes] );
                  //outResult->SetResultInfo( infos[testRes] );

                  inResults_fut.Push_back( inResult_fut );
                  inResults_ref.Push_back( inResult_ref );
                  //outResults.Push_back( outResult );
                  std::cout << inResult_ref->GetEntityList()->GetName();
                  if (iRegion < regions.GetSize()-1) std::cout << ", ";
              }
              std::cout << ")\n";
          }
          // now check mesh results
          //Vector<Double> maxResVec_ref;
          //maxResVec_ref.Resize( inResults_ref.GetSize() );

          // iterate over all time/frequency steps
          for( UInt iStep = 0; iStep < numSteps_ref[actMsStep]; iStep++ ) {
              Double maxL2_step = 0.0;//, maxL2rel
              Double maxL2rel_step = 0.0;
              std::string infoStringRel_step;
              UInt actStepNum = iStep+1;
              Double actStepVal = stepVals[actStepNum];
              // iterate over all results
              for( UInt iRes = 0; iRes < inResults_fut.GetSize(); iRes++) 
              {
                 const shared_ptr<ResultInfo>& resInfo = inResults_fut[iRes]->GetResultInfo();
                 if(resultSteps.count(resInfo) == 0 || resultSteps[resInfo].count(actStepNum) == 0) 
                   continue;
                  // obtain both result objects for current step
                  input_fut->GetResult( actMsStep, actStepNum, inResults_fut[iRes], isHistory );
                  input_ref->GetResult( actMsStep, actStepNum, inResults_ref[iRes], isHistory );
                  // cast result objects, get vector and calculate difference vector
                  const std::string& resultName = inResults_fut[iRes]->GetResultInfo()->resultName;
                  const std::string& entityName = inResults_fut[iRes]->GetEntityList()->GetName();
                  UInt numDofs = inResults_fut[iRes]->GetResultInfo()->dofNames.GetSize();
                  if (numDofs == 0) numDofs = 1;
                  if( types[actMsStep] == BasePDE::STATIC || types[actMsStep] == BasePDE::TRANSIENT ) {
                      Vector<Double> & inVec_fut = dynamic_cast<Result<Double>& >(*inResults_fut[iRes]).GetVector();
                      Vector<Double> & inVec_ref = dynamic_cast<Result<Double>& >(*inResults_ref[iRes]).GetVector();
                      checkL2Result(inVec_fut, inVec_ref, actMsStep, actStepNum, actStepVal,
                                    resultName, entityName, numDofs, "real",
                                    maxL2_step, maxL2rel_step, infoStringRel_step);
                  } else if(types[actMsStep] == BasePDE::EIGENVALUE || types[actMsStep] == BasePDE::BUCKLING ||  types[actMsStep] == BasePDE::EIGENFREQUENCY || types[actMsStep] == BasePDE::HARMONIC || types[actMsStep] == BasePDE::MULTIHARMONIC || types[actMsStep] == BasePDE::HARMONIC25D) {
                      Vector<Complex> & inVec_fut = dynamic_cast<Result<Complex>& >(*inResults_fut[iRes]).GetVector();
                      Vector<Complex> & inVec_ref = dynamic_cast<Result<Complex>& >(*inResults_ref[iRes]).GetVector();
                      // normalize modes
                      if (types[actMsStep] == BasePDE::EIGENFREQUENCY || types[actMsStep] == BasePDE::BUCKLING|| types[actMsStep] == BasePDE::EIGENVALUE) {
                        int i_ref = 0;
                        Complex norm_ref = inVec_ref.MaxAbs(i_ref);
                        Complex norm_fut = inVec_fut[i_ref];
                        std::cout << "  MAC(fut,ref) = " << inVec_fut.MAC(inVec_ref) << "\n";
                        std::cout << "  Normalizing modes by "<< norm_fut << " (fut), "<< norm_ref << " (ref)\n";
                        inVec_fut.ScalarDiv(norm_fut);
                        inVec_ref.ScalarDiv(norm_ref);
                      }
                      checkL2Result(inVec_fut, inVec_ref, actMsStep, actStepNum, actStepVal,
                                    resultName, entityName, numDofs, "complex",
                                    maxL2_step, maxL2rel_step, infoStringRel_step);
                  } else {
                      EXCEPTION("not impl");
                  } // switch: Analysitype
              } // loop over results
              if (maxL2_step > maxL2_sequence ) maxL2_sequence = maxL2_step;
              if (maxL2rel_step >= maxL2rel_sequence ) {
                  maxL2rel_sequence = maxL2rel_step;
                  std::stringstream tmpstr;
                  tmpstr << "step "<< actStepNum << " with step value " << actStepVal << ": " << infoStringRel_step;
                  infoStringRel_sequence = tmpstr.str();
              }
          }// time step
          std::cout << "\nMaximal error withing multi sequence step: L2-Norm: " << maxL2_sequence << " (difference), " << maxL2rel_sequence <<" (relative) @ " <<infoStringRel_sequence<<"\n";
          if ( maxL2_sequence > maxL2 ) maxL2 = maxL2_sequence;
          if ( maxL2rel_sequence >= maxL2rel ) {
              maxL2rel = maxL2rel_sequence;
              std::stringstream tmpstr;
              tmpstr << "MultiSequenceStep "<< actMsStep << ", " << infoStringRel_sequence;
              infoStringRel = tmpstr.str();
          }
      } // sequence step
      delete ptGrid_fut;
      delete ptGrid_ref;
      maxDiffInfo = infoStringRel;
      return std::max(maxL2,maxL2rel);
  } //CheckL2

} //Namespace CFSTool


int main(int argc, char** argv)
{

  // todo: do better once! - Fabian
  CFSTool::InitEnums(); 
  ElemShape::Initialize();
  shared_ptr<LogConfigurator> logConf;

  domain = NULL;
  
  std::string maxDiffResultName; // that's chaos man!
  std::string infoFileName;
  try
  {
    param.reset(new ParamNode( ParamNode::PASS, ParamNode::ELEMENT));

    CFSTool::ParamsInit(argc, argv, param, logConf);

    // Switch this flag tc true for debugging
    Exception::segfault_ =  param->Get("forceSegFault")->As<bool>();

    // Register callback function with exception class for warning
    Exception::SetCallbackWarn(CFSTool::PrintWarning);

    // Print out hard-coded node and element limits and exit.
    if (param->Has("printLimits"))
    {
      CFSTool::PrintLimits(param->Get("printLimits")->As<std::string>());
      return 0;
    }

    std::cout << std::endl
              << "=======================================================================" << std::endl;
    std::cout << " CFSTOOL - File Conversions/Comparisons for openCFS" << std::endl << std::endl
              << " v. " << CFS_VERSION << " - '" << CFS_NAME << "'"
              << " (rev " << CFS_GIT_COMMIT << ")" << std::endl
              << " compiled " << __DATE__
              << " as " << CMAKE_BUILD_TYPE << std::endl;
    std::cout << "=======================================================================\n";
    //  -m relL2diff --checkLimits false --eps 1e-6 mech.h5ref results_hdf5/mech.cfs
    std::cout << "Command line arguments: ";
    for (int i = 1; i < argc; ++i)
       std::cout << argv[i] << " ";
    std::cout << "\n";    
    
    std::string param_mode = param->Get("mode")->As<std::string>();

    // get filenames from parameter
    std::string file1 = param->Get("file1")->As<std::string>();
    std::string file2 = param->Get("file2")->As<std::string>();
    std::string file3 = param->Get("file3")->As<std::string>();
    std::string file4 = param->Get("file4")->As<std::string>();

    // Initialize vector with output fields
    UInt num_files = 0;
    if (file1 != "")
    {
      ++num_files;
      if (file2 != "")
      {
        ++num_files;
        if (file3 != "")
        {
          ++num_files;
          if (file4 != "")
          {
            ++num_files;
          }
        }
      }
    }

    // This is necessary to run, but I do not know what it is for
    info.reset(new ParamNode(ParamNode::INSERT, ParamNode::ELEMENT));
    infoFileName = param_mode + "_" + file1 + ".info.xml";
    info->SetName("cfsInfo"); 

    if (param_mode == "calcAverage")
    {
      if (file3 != "")
      {
        EXCEPTION( "Too many arguments, please only provide two files. (in- and output file)" );
      }
      file3 = file2;
      if (num_files != 2)
      {
        EXCEPTION( "Please provide a reference file and output File" );
      }
      CFSTool::calcAverage(file1, file3);
    } else if (param_mode == "convert") {
      if (file3 != "")
      {
        EXCEPTION( "Too many arguments, please only provide two files. (in- and output file)" );
      }
      file3 = file2;
      if (num_files != 2)
      {
        EXCEPTION( "Please provide 'input_file' and 'output_file'" );
      }
      CFSTool::Convert( file1, file3 );
    } else if (param_mode == "scalardiff") {
      Double tolerance = param->Get("eps")->As<Double>();
      if (num_files != 2)
      {
        EXCEPTION( "Please provide 'reference_file' and 'file_under_test', detected files: " << num_files
                    << " file1='" << file1 << "' file2='" << file2 << "'");
      }
      Double maxDiffMesh = 0.0, maxDiffHist = 0.0;
      std::cout << "Checking for mesh results:\n==========================\n";
      maxDiffMesh = CFSTool::Diff( file1, file2, "", true, false, maxDiffResultName);
      std::cout << "<DartMeasurement name=\"scalardiff (mesh)\" type=\"numeric/double\">"<<maxDiffMesh<<"</DartMeasurement>\n";
      std::cout << "Checking for history results:\n=============================\n";
      maxDiffHist = CFSTool::Diff( file1, file2, "", true, true, maxDiffResultName );
      std::cout << "<DartMeasurement name=\"scalardiff (history)\" type=\"numeric/double\">"<<maxDiffHist<<"</DartMeasurement>\n";
      Double maxDiff = std::max( maxDiffMesh, maxDiffHist );
      std::cout << "<DartMeasurement name=\"scalardiff\" type=\"numeric/double\">"<<maxDiff<<"</DartMeasurement>\n";
      if( maxDiff > tolerance ) {
        std::cout << "error: maximum difference " << maxDiff << " for '" << maxDiffResultName << "' > " << tolerance
                  << " for '" << file1 << "' and '" << file2 << std::endl;
        exit(EXIT_FAILURE);
      } else {
        std::cout << "  No differences larger than tolerance found.\n";
        exit(EXIT_SUCCESS);
      }
    } else if (param_mode == "absL2diff"){
      Double tolerance = param->Get("eps")->As<Double>();
      if (num_files != 2)
      {
          EXCEPTION( "Please provide 'reference_file' and 'file_under_test', detected files: " << num_files
                          << " file1='" << file1 << "' file2='" << file2 << "'");
      }
      std::cout << "#####################################################\n";
      std::cout << "#                   Mode = absL2diff                #\n";
      std::cout << "#####################################################\n";
      std::cout << "Checking for mesh results:\n"
      << "==========================\n";
      Double maxDiffMesh,maxDiffMeshRel;
      CFSTool::CheckL2(file1,file2,false,maxDiffMesh,maxDiffMeshRel,maxDiffResultName,tolerance);
      std::cout << "Maximum L2 norm = " << maxDiffMesh << " @ "<< maxDiffResultName << "\n";
      std::cout << "<DartMeasurement name=\"absL2diff (mesh)\" type=\"numeric/double\">"<<maxDiffMesh<<"</DartMeasurement>\n";
      std::cout << "Checking for history results:\n"
      << "==========================\n";
      Double maxDiffHist,maxDiffHistRel;
      CFSTool::CheckL2(file1,file2,true,maxDiffHist,maxDiffHistRel,maxDiffResultName,tolerance);
      std::cout << "Maximum L2 norm = " << maxDiffHist << " @ "<< maxDiffResultName <<"\n";
      std::cout << "<DartMeasurement name=\"absL2diff (history)\" type=\"numeric/double\">"<<maxDiffHist<<"</DartMeasurement>\n";
      Double err = std::max(maxDiffMesh,maxDiffHist);
      std::cout << "\n";
      std::cout << "<DartMeasurement name=\"absL2diff\" type=\"numeric/double\">"<<err<<"</DartMeasurement>\n";
      if ( err < tolerance ) {
          std::cout << "========================================================\n";
          std::cout << "Maximum L2 norm = " << err << " < " << tolerance << "\n";
          exit(EXIT_SUCCESS);
      } else {
          std::cout << "========================================================\n";
          std::cout << "Maximum L2 norm = " << err << " > " << tolerance << "\n";
          exit(EXIT_FAILURE);
      }
    } else if (param_mode == "relL2diff") {
        Double tolerance = param->Get("eps")->As<Double>();
        if (num_files != 2)
        {
            EXCEPTION( "Please provide 'reference_file' and 'file_under_test', detected files: " << num_files
                            << " file1='" << file1 << "' file2='" << file2 << "'");
        }
        std::cout << "#####################################################\n";
        std::cout << "#           Checking for mesh results ...           #\n";
        std::cout << "#####################################################\n";
        Double absDiff,maxDiffMeshRel;
        CFSTool::CheckL2(file1,file2,false,absDiff,maxDiffMeshRel,maxDiffResultName,tolerance);
        std::cout << "\n";
        if ( maxDiffMeshRel < tolerance ) 
            std::cout << "OK: relative mesh results maximal err: L2-norm=" << maxDiffMeshRel << " < " << tolerance << " -> all fine for mesh";
        else 
            std::cout << "ERROR: relative mesh results maximal err: L2-norm=" << maxDiffMeshRel << " >= " << tolerance << " @ "<< maxDiffResultName << " continue";

        std::cout << "\n\n<DartMeasurement name=\"absL2diff (mesh)\" type=\"numeric/double\">"<<absDiff<<"</DartMeasurement>\n";
        std::cout << "<DartMeasurement name=\"relL2diff (mesh)\" type=\"numeric/double\">"<<maxDiffMeshRel<<"</DartMeasurement>\n";
        std::cout << "\n";
        std::cout << "#####################################################\n";
        std::cout << "#           Checking for history results ...        #\n";
        std::cout << "#####################################################\n";
        Double maxDiffHistRel;
        CFSTool::CheckL2(file1,file2,true,absDiff,maxDiffHistRel,maxDiffResultName,tolerance);
        std::cout << "========================================================\n";
        std::cout << "\n";

        std::cout << "-> Final maximal relative history error: L2-norm= " << maxDiffHistRel << " @ " << maxDiffResultName << "\n\n";
        std::cout << "<DartMeasurement name=\"absL2diff (history)\" type=\"numeric/double\">"<<absDiff<<"</DartMeasurement>\n";
        std::cout << "<DartMeasurement name=\"relL2diff (history)\" type=\"numeric/double\">"<<maxDiffHistRel<<"</DartMeasurement>\n";
        std::cout << "\n";
        Double err = std::max(maxDiffMeshRel,maxDiffHistRel);
        std::cout << "<DartMeasurement name=\"relL2diff\" type=\"numeric/double\">"<<err<<"</DartMeasurement>\n";
        std::cout << "========================================================\n";
        if ( err < tolerance ) {
            std::cout << "OK: History results maximal err=" << err << " < " << tolerance << " -> all fine for history\n";
            exit(EXIT_SUCCESS);
        } else {
            std::cout << "Error: History results maximal err=" << err << " >= " << tolerance << "\n";
            exit(EXIT_FAILURE);
        }
    } else if (param_mode == "meshdiff") {
      if (num_files != 3)
      {
        EXCEPTION( "Please provide 'reference_file', 'file_under_test' and 'out_file'" );
      }
      CFSTool::Diff( file1, file2, file3, false, false, maxDiffResultName);
    } else if (param_mode == "meshdiffnormed") {
      if (num_files != 3)
      {
        EXCEPTION( "Please provide 'reference_file', 'file_under_test' and 'out_file'" );
      }
      CFSTool::Diff( file1, file2, file3, true, false, maxDiffResultName);
    } else if (param_mode == "wvt") {

      CFSTool::WVT wvt( param, info );
      wvt.PostProcess();

    } else {
      EXCEPTION( "No such mode: " << param_mode <<". See help for available modes." );
      return EXIT_FAILURE;
    }
  } catch(std::exception& ex) {
    std::cerr << "The following error occurred during program execution:\n\n" << ex.what();

    if (info != NULL)
    {
      info->Get(ParamNode::FAIL)->SetValue(ex.what());
      info->ToFile(infoFileName);
    }

    return -1;
  }
  info->ToFile(infoFileName);
  return 0;
}
