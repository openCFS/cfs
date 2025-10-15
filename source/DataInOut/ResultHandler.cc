// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
#include <cmath>



#include "ResultHandler.hh"

#include "DataInOut/ProgramOptions.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "Domain/Domain.hh"
#include "Domain/Results/ResultFunctor.hh"
#include "DataInOut/PostProc.hh"
#include "DataInOut/SimInOut/TextOutput/TextSimOutput.hh"
#include "DataInOut/SimInput.hh"
#include "FeBasis/FeSpace.hh"
#include "Utils/Timer.hh"


namespace CoupledField {

  DEFINE_LOG(resHandler, "resultHandler")
    
  ResultHandler::ResultHandler( PtrParamNode paramNode) {
    param_ = paramNode;
    sequenceStep_ = 1;
    actStep_ = 0;
    finalResultExists_ = false;
    actStepVal_ = 0.0;
    numSteps_ = 0;
    analysisType_ = BasePDE::NO_ANALYSIS;
  }


  ResultHandler::~ResultHandler() {

    // not much to do here
  }

  void ResultHandler::
  RegisterResult( shared_ptr<BaseResult> sol,
                  shared_ptr<ResultFunctor> fnc,
                  UInt sequenceStep,
                  UInt saveBegin, UInt saveInc,
                  UInt saveEnd, 
                  const StdVector<std::string> & outDestNames,
                  const std::string& postProcName,
                  bool writeResult,
                  bool isHistory ) {

    ResultInfo & actDof = *(sol->GetResultInfo());

    // From here on, the entry in the ResultInfo.resultName
    // is needed instead of the Resultinfo.resultType
    if( actDof.resultName == "" ) {
      actDof.resultName = SolutionTypeEnum.ToString(actDof.resultType);
    }

    LOG_DBG(resHandler) << "Registering result:";
    LOG_DBG(resHandler) << "-------------------";
    LOG_DBG(resHandler) << "name: " << actDof.resultName;
    LOG_DBG(resHandler) << "dofs: " << actDof.dofNames.Serialize();
    LOG_DBG(resHandler) << "resultList: " << sol->GetEntityList()->GetName();
    LOG_DBG(resHandler) << "saveBegin: " << saveBegin;
    LOG_DBG(resHandler) << "saveEnd: " << saveEnd;
    LOG_DBG(resHandler) << "saveInc: " << saveInc;
    LOG_DBG(resHandler) << "writeResult: " << writeResult;
    LOG_DBG(resHandler) << "isFinal: 0";
    LOG_DBG(resHandler) << "outputDest: " << outDestNames.Serialize();
    LOG_DBG(resHandler) << "postProcId: " << postProcName;
    LOG_DBG(resHandler) << "-------------------";

    std::cout << "Registering result: (ResultHandler.cc)" << std::endl;
    std::cout << "-------------------" << std::endl;
    std::cout << "name: " << actDof.resultName << std::endl;
    std::cout << "dofs: " << actDof.dofNames.Serialize() << std::endl;
    std::cout << "resultList: " << sol->GetEntityList()->GetName() << std::endl;
    std::cout << "saveBegin: " << saveBegin << std::endl;
    std::cout << "saveEnd: " << saveEnd << std::endl;
    std::cout << "saveInc: " << saveInc << std::endl;
    std::cout << "writeResult: " << writeResult << std::endl;
    std::cout << "isFinal: 0" << std::endl;
    std::cout << "outputDest: " << outDestNames.Serialize() << std::endl;
    std::cout << "postProcId: " << postProcName << std::endl;
    std::cout << "-------------------" << std::endl;



    // check, if result context was already created
    shared_ptr<ResultContext> actContext;
  
    if( resultContexts_.find(sol) != resultContexts_.end() ) {
      EXCEPTION( "Context was already found" );
    }

    // create new ResultContext
    actContext = shared_ptr<ResultContext>( new ResultContext() );
    actContext->result = sol;
    actContext->functor =  fnc;
    actContext->sequenceStep = sequenceStep;
    actContext->saveBegin = saveBegin;
    actContext->saveEnd = saveEnd;
    actContext->saveInc = saveInc;
    actContext->writeResult = writeResult;
    actContext->isFinal = false;
    actContext->isHistory = isHistory;
    resultContexts_[sol] = actContext;
    
    // if no destination output was provided, try to find default one
    StdVector<std::string> newDest;
    if (outDestNames.GetSize() == 1 ) {
      if( outDestNames[0] == "" ) {
        GetDefaultOutputs( sol, newDest );
      } else {
        newDest = outDestNames;
      }
    } else {
      newDest = outDestNames;
    }
    
    // add all output destinations to the the context, if
    // the result should be written to file
    if( writeResult ) {
      for( UInt i = 0; i < newDest.GetSize(); i++ ) {
        
        // check, if output is registered
        if( simOutputHandlers_.find( newDest[i] ) == simOutputHandlers_.end() )
          EXCEPTION( "Output writer '" << newDest[i] << "' was not registered yet!" );

        LOG_DBG(resHandler) << "Registering output '" << newDest[i]
                            << "' with result '" << actDof.resultName << "'";

        actContext->outputIds.Push_back( newDest[i] );

        // skip if only streaming is set and this is no streamingHandler
        if (streamOnly && !simOutputHandlers_[newDest[i]]->IsStreaming())
          continue;

        // register results also at the output writer class
        simOutputHandlers_[newDest[i]]->RegisterResult( sol, saveBegin,
                                               saveInc, saveEnd,
                                               actContext->isHistory );
      }
    }

    // insert context to list
    contexts_.insert( actContext );
    
    // Check for postprocessing method(s), which might also
    // produce new result contexts
    if( postProcName != "" ) {
      RegisterResultRec( *actContext, postProcName );
    }
    LOG_DBG(resHandler) << "Finished registering result";
  }
  
  void ResultHandler::BeginMultiSequenceStep( UInt step, BasePDE::AnalysisType type, UInt numSteps )
  {
    LOG_DBG(resHandler) << "BMS: step=" << step << " n=" << numSteps;

    //assert(numSteps >= numSteps_);

    // store current sequencestep
    sequenceStep_ = step;
    analysisType_ = type;
    numSteps_ = numSteps;
    
    // Iterate over all outfiles to notify them about a new
    // multisequence step
    std::map<std::string, shared_ptr<SimOutput> >::iterator it;
    it = simOutputHandlers_.begin();

    for( ; it != simOutputHandlers_.end(); it++ ) {
      // skip if only streaming is set and this is no streamingHandler
      if (streamOnly && !it->second->IsStreaming())
        continue;

      it->second->BeginMultiSequenceStep( step, type, numSteps );
    }
  }

  void ResultHandler::BeginStep( UInt stepNum, Double stepVal ) {

    LOG_DBG(resHandler) << "Begin step " << stepNum;

    // remember current step values
    actStep_ = stepNum;
    actStepVal_ = stepVal;

    // check all result contexts, which one of them has to be updated
    // (determined by saveBegin, saveEnd and saveInc )
    isNeeded_.clear();
    isUpdated_.clear();
    for(auto contextIt = contexts_.begin(); contextIt != contexts_.end(); contextIt++) {
      BaseResult& actResult  = *((*contextIt)->result);
      LOG_DBG(resHandler) << "IsNeeded for '" << actResult.GetResultInfo()->resultName
                          << "' on '"<< actResult.GetEntityList()->GetName() << "' is '"
                          << (IsOutput( **contextIt ) == true ? "true" : "false") << "'";
      if(IsOutput( **contextIt))
        isNeeded_.insert((*contextIt)->result);
    }
    
    // Trigger function for all output files
    for(auto it = simOutputHandlers_.begin(); it != simOutputHandlers_.end(); it++ ) {
      // skip if only streaming is set and this is no streamingHandler
      if (streamOnly && !it->second->IsStreaming())
        continue;

      it->second->BeginStep(stepNum, stepVal);
    }
    
    LOG_DBG(resHandler) << "Finished beginning of new step";
  }

  void ResultHandler::UpdateResult( shared_ptr<BaseResult> sol ) {


    shared_ptr<Timer> timer = domain->GetInfoRoot()->Get(ParamNode::HEADER)->Get("results/timer")->AsTimer();
    timer->Start();

    LOG_DBG(resHandler) << "UR: " << sol->GetResultInfo()->ToString();
    
    // check, if result is to be updated
    if(isNeeded_.find(sol) == isNeeded_.end())
      WARN("Result '" << sol->GetResultInfo()->resultName << "' on entitylist '" << sol->GetEntityList()->GetName() << "' is not needed in step " << actStep_);

    // Set flag for update to true
    isUpdated_.insert( sol );

    LOG_DBG(resHandler) << "Result '" << sol->GetResultInfo()->resultName << "' was provided on '" << sol->GetEntityList()->GetName() << "' in step " << actStep_;

    // Update results recursively (for postprocessing results )
    UpdateResultRec(*(resultContexts_[sol]));
    LOG_DBG(resHandler) << "UR: finished";

    timer->Stop();
  }

  void ResultHandler::UpdateResults() {
    LOG_DBG(resHandler) << "UR: needed=" << isNeeded_.size();

    // iterate over all results, which are needed
    for(auto it = isNeeded_.begin(); it != isNeeded_.end(); it++)
    {
      ResultContext& actContext = *(resultContexts_[*it]);

      if(actContext.functor) {
        LOG_DBG(resHandler) << "Evaluating result '" << SolutionTypeEnum.ToString(actContext.result->GetResultInfo()->resultType )
                            << "' on '" << actContext.result->GetEntityList()->GetName() << "'";

        shared_ptr<Timer> timer = domain->GetInfoRoot()->Get(ParamNode::HEADER)->Get("results/timer")->AsTimer();
        timer->Start();
        actContext.functor->EvalResult(actContext.result);
        timer->Stop();

        UpdateResult(actContext.result);

      }
    }
  }

  void ResultHandler::FinishStep( )
  {
    LOG_DBG(resHandler) << "FinishStep " << actStep_ << " enter";

    // -----------------------
    // First, update results
    // -----------------------
    UpdateResults();
    
    // shared amongst e.g. WriteResults, UpdateResults, FinishStep
    shared_ptr<Timer> timer = domain->GetInfoRoot()->Get(ParamNode::HEADER)->Get("results/timer")->AsTimer();
    timer->Start();

    // === Primary results ===
    // iterate over all results, which are needed
    for(auto it = isNeeded_.begin(); it != isNeeded_.end(); it++ )
    {
      // store context
      ResultContext & actContext = *(resultContexts_[*it]);
      BaseResult & actResult  = *(actContext.result);

      if(actContext.sequenceStep != sequenceStep_)
        continue;

      LOG_DBG(resHandler) << "Checking result '"<< actResult.GetResultInfo()->resultName
                          << "' on '" << actResult.GetEntityList()->GetName() << "' for writing out";

      // Check, if result was updated at all
      if(isUpdated_.find(*it) == isUpdated_.end())
        WARN("Result '" << actResult.GetResultInfo()->resultName << "' on '"
          << (*it)->GetEntityList()->GetName() << "' was not provided in step " << actStep_);

      std::cout << actResult.ToString();


      // check, if result is to be written 
      if( actContext.writeResult && (!actContext.isFinal ) ) {
        // iterate over all outputs
        for(UInt iOut = 0; iOut < actContext.outputIds.GetSize(); iOut++) {
          // skip if only streaming is set and this is no streamingHandler
          if (streamOnly && !simOutputHandlers_[actContext.outputIds[iOut]]->IsStreaming())
            continue;

          std::string outId = actContext.outputIds[iOut];
          // =================================================
          // Check, if output writer "lives" on the same grid.
          // If not, we have to perform mapping
          // =================================================
          ResultInfo & info = *(actContext.result->GetResultInfo());
          std::cout << "OutputIds_[outId]: " << outGridIds_[outId] << std::endl;
          std::cout << "info.definedOn: " << info.definedOn << std::endl;
          if( outGridIds_[outId] != "default" &&
              (info.definedOn == ResultInfo::NODE || 
              info.definedOn == ResultInfo::ELEMENT ||
              info.definedOn == ResultInfo::SURF_ELEM ) ) {
            
            // security check: if no result functor is present, we leave
            if( !actContext.functor )
              continue;

            // obtain destination grid and get the element / node list
            Grid * destGrid = domain->GetGrid( outGridIds_[outId] );
            shared_ptr<BaseResult> res;
            
            // perform interpolation
            InterpolateRes( actContext, destGrid, res );
            simOutputHandlers_[actContext.outputIds[iOut]]->AddResult( res);
         
          } else {
            // Standard case, no interpolation necessary
            
            // Add current result to given output file
            LOG_DBG(resHandler) << "Adding result '" << actResult.GetResultInfo()->resultName
                << "' on '" << actResult.GetEntityList()->GetName() << "' to '"
                << simOutputHandlers_[actContext.outputIds[iOut]]->GetName() << "'";
            simOutputHandlers_[actContext.outputIds[iOut]]->AddResult( actContext.result );
          }
        }
      }
      // In any case: Process also related secondary (postprocessing) results
      FinishStepRec( actContext );
    }

  // Trigger writing of all output writers
  for(auto fileIt = simOutputHandlers_.begin(); fileIt != simOutputHandlers_.end(); fileIt++ ) {
    // skip if only streaming is set and this is no streamingHandler
    if (streamOnly && !fileIt->second->IsStreaming())
      continue;

    LOG_DBG(resHandler) << "Finishing step for output '" << fileIt->second->GetName() << "'";
    fileIt->second->FinishStep();
  }    

  // something is written, so update also info
  domain->GetInfoRoot()->ToFile();

  LOG_DBG(resHandler) << "FinishStep " << actStep_ << " done";
  timer->Stop();
}

  void ResultHandler::FinishMultiSequenceStep() {

    // Before finishing the multisequence step, we have to
    // finalize all postprocessing steps and write all the results
    // which are only defined for the last step
    std::map<std::string, shared_ptr<SimOutput> >::iterator fileIt;
    if( finalResultExists_ ) {

      LOG_DBG(resHandler) << "There exist some results for the final step";
      actStep_++;
      actStepVal_ = 0.0;

      std::string type;
      // Trigger new step at all output writers
      
      for( fileIt = simOutputHandlers_.begin(); 
      fileIt != simOutputHandlers_.end(); fileIt++ ) {
        // skip if only streaming is set and this is no streamingHandler
        if (streamOnly && !fileIt->second->IsStreaming())
          continue;

        fileIt->second->BeginStep( actStep_, actStepVal_);
      }    

      // === Primary results ===
      // iterate over all results, which are needed
      for(auto it = contexts_.begin(); it != contexts_.end(); it++ )
      {
        // store context
        ResultContext & actContext = **it;
        BaseResult & actResult  = *(actContext.result);
        LOG_DBG(resHandler) << "Testing result '" 
        << actResult.GetResultInfo()->resultName << "' on '"
        << actResult.GetEntityList()->GetName()
        << "'";
        // check, if result is to be written 
        if( actContext.writeResult &&  actContext.isFinal )
        {
          // iterate over all outputs
          for( UInt iOut = 0; iOut < actContext.outputIds.GetSize(); iOut++ )
          {
            // skip if only streaming is set and this is no streamingHandler
            if (streamOnly && !simOutputHandlers_[actContext.outputIds[iOut]]->IsStreaming())
              continue;

            // Add current result to given output file
            simOutputHandlers_[actContext.outputIds[iOut]]->AddResult( actContext.result );
            LOG_DBG(resHandler) << "Adding final result '" << type << "' on '"
            << actResult.GetEntityList()->GetName()
            << "' to '" << simOutputHandlers_[actContext.outputIds[iOut]]->GetName()
            << "'";
          }
        }

        // In any case: Process also related secondary (postprocessing) results
        FinishMultiSequenceStepRec( actContext );
      }
      // Finish newly created step
      for(fileIt = simOutputHandlers_.begin(); fileIt != simOutputHandlers_.end(); fileIt++) {
        // skip if only streaming is set and this is no streamingHandler
        if (streamOnly && !fileIt->second->IsStreaming())
          continue;

        fileIt->second->FinishStep();
      }
    }
    
    
    // Iterate over all outfiles to notify them about the end of a
    // multisequence step
    for(auto it = simOutputHandlers_.begin(); it != simOutputHandlers_.end(); it++) {
      // skip if only streaming is set and this is no streamingHandler
      if (streamOnly && !fileIt->second->IsStreaming())
        continue;

      it->second->FinishMultiSequenceStep( );
    }

    // Delete all results, as they are registered newly for each
    // sequence step
    contexts_.clear();
    resultContexts_.clear();
    isNeeded_.clear();
    isUpdated_.clear();
  }

  void ResultHandler::RegisterResultRec( ResultContext& actContext, const std::string& postProcName )
  {
    LOG_DBG(resHandler) << "Registering (recursively) result of postProc '" << postProcName << "'";

    // if no postProcName is set, just leave
    if( postProcName == "" )
      return;

    // Hack: As there might be times for a sequence step with index 0
    // (e.g. grid results etc.) we look for a suitable sequence step for
    // a step 1 in case 
    
    // fetch current postProcNode
    UInt seqStep = 1;
    if( sequenceStep_ > 1 ) {
      seqStep = sequenceStep_;
    }
    PtrParamNode postProcNode, ppListNode;
    PtrParamNode seqNode = param_->GetByVal("sequenceStep", "index", seqStep);
    if(seqNode)
      ppListNode = seqNode->Get("postProcList");
      
    if(ppListNode)
      postProcNode = ppListNode->GetByVal("postProc", "id", postProcName);
    
    if(!postProcNode)
      EXCEPTION( "A Postprocessing section for '" << postProcName << "' does not exist" );

    // Fetch postprocs
    StdVector<shared_ptr<PostProc> > postProcs;
    PostProc::CreatePostProc(postProcNode,  domain->GetGrid(), postProcs);
                              
    // iterate over all postprocs 
    for( UInt i = 0; i < postProcs.GetSize(); i++ )
    {
      // register current solution with new object
      postProcs[i]->SetResult( actContext.result );

      // create new resultcontext and link it to current one
      shared_ptr<ResultContext> nextContext = shared_ptr<ResultContext>(new ResultContext());
      nextContext->result = postProcs[i]->GetOutputResult();
      nextContext->functor = actContext.functor;
      nextContext->sequenceStep = actContext.sequenceStep;
      nextContext->saveBegin = actContext.saveBegin;
      nextContext->saveEnd = actContext.saveEnd;
      nextContext->saveInc = actContext.saveInc;
      if( postProcs[i]->IsHistory() || actContext.isHistory )
        nextContext->isHistory = true;
      else
        nextContext->isHistory = false;

      nextContext->writeResult = postProcs[i]->IsWriteResult();
      if(postProcs[i]->GetReductionType() == PostProc::TIME_FREQ || actContext.isFinal) {
        nextContext->isFinal = true;
        finalResultExists_ = true;
      } else {
        nextContext->isFinal = false;
      }
      
      // fetch suitable output class
      StdVector<std::string> oldDest, outDest;
      postProcs[i]->GetOutDestNames( oldDest );
      LOG_DBG(resHandler) << "outputDest of " << postProcName << " is '" << oldDest.Serialize() << "'";
      LOG_DBG(resHandler) << "size of outputDest is " << outDest.GetSize() << std::endl;
      if( oldDest.GetSize() == 1 ) {
        if( oldDest[0] == "" ) {
          GetDefaultOutputs( nextContext->result, outDest );
        } else {
          outDest = oldDest;
        }
      } else {
        outDest = oldDest;
      }
      
      
      // Write log information
      ResultInfo & nextInfo = *(nextContext->result->GetResultInfo());
      LOG_DBG(resHandler) << "Registering result (RECURSIVELY):";
      LOG_DBG(resHandler) << "-----------------------------------";
      LOG_DBG(resHandler) << "name: " << nextInfo.resultName;
      LOG_DBG(resHandler) << "dofs: " << nextInfo.dofNames.Serialize();
      LOG_DBG(resHandler) << "resultList: " 
                          << nextContext->result->GetEntityList()->GetName();
      LOG_DBG(resHandler) << "saveBegin: " << nextContext->saveBegin;
      LOG_DBG(resHandler) << "saveEnd: " << nextContext->saveEnd;
      LOG_DBG(resHandler) << "saveInc: " << nextContext->saveInc;
      LOG_DBG(resHandler) << "writeResult: " << nextContext->writeResult;
      LOG_DBG(resHandler) << "isFinal: " << nextContext->isFinal;
      LOG_DBG(resHandler) << "isHistory: " << nextContext->isHistory;
      LOG_DBG(resHandler) << "outputDest: " << outDest.Serialize();
      LOG_DBG(resHandler) << "postProcName: " 
                          << postProcs[i]->GetNextPostProcName() << std::endl;

      for(UInt iOut = 0; iOut < outDest.GetSize(); iOut++)
      {
        // skip if only streaming is set and this is no streamingHandler
        if (streamOnly && !simOutputHandlers_[outDest[iOut]]->IsStreaming())
          continue;

        if( simOutputHandlers_.find( outDest[iOut] ) == simOutputHandlers_.end() )
          EXCEPTION( "Output writer '" << outDest[i] << "' was not registered yet!" );
        
        nextContext->outputIds.Push_back( outDest[iOut] );
      
        // register results also at the output writer class
        if( nextContext->writeResult ) {
          simOutputHandlers_[outDest[iOut]]->RegisterResult(postProcs[i]->GetOutputResult(),
                           nextContext->saveBegin, nextContext->saveInc,
                           nextContext->saveEnd, nextContext->isHistory );
        }
      }
      
      // store postproc and result in current context2
      actContext.nextContexts.Push_back( nextContext );
      actContext.postProcs.Push_back( postProcs[i] );
      
      // call method iteratively
      if(postProcs[i]->GetNextPostProcName() != "")
        RegisterResultRec(*nextContext, postProcs[i]->GetNextPostProcName());
    }
    LOG_DBG(resHandler) << "Finished registering result of postProc '" << postProcName << "'";
  }
    
  void ResultHandler::FinishStepRec( ResultContext& actContext ) {
    
    // ensure that we have for each postprocessing method
    // also the related result context
    assert(actContext.postProcs.GetSize() == actContext.nextContexts.GetSize());

    // iterate over all contexts
    for( UInt i = 0; i < actContext.nextContexts.GetSize(); i++ ) {

      ResultContext & next = *(actContext.nextContexts[i]);
      LOG_DBG(resHandler) << "Checking result '" << next.result->GetResultInfo()->resultName
                          << "' on '" << next.result->GetEntityList()->GetName() << "' for writing out";

      // Check, if next result is only to be computed in final stage
      if( next.writeResult != false  && next.isFinal != true) {
        
        // iterate over all outputs
        for( UInt iOut = 0; iOut < next.outputIds.GetSize(); iOut++ ) {

          // skip if only streaming is set and this is no streamingHandler
          if (streamOnly && !simOutputHandlers_[next.outputIds[iOut]]->IsStreaming())
            continue;

          // Add current result to given output file
          simOutputHandlers_[next.outputIds[iOut]]->AddResult( next.result );
          BaseResult & nextResult  = *(next.result);
          LOG_DBG(resHandler) << "Adding result '" << nextResult.GetResultInfo()->resultName
                              << "' on '" << nextResult.GetEntityList()->GetName()
                              << "' to '" << simOutputHandlers_[next.outputIds[iOut]]->GetName() << "'";
        }
      }

      // Call method recursively to write all outputs of next context
      // Note: We only may write the result, if it is not a "final" result,
      // i.e. it is written only once at the end. In this case, the result
      // is written via the method FinalizeRec()
      if(!next.isFinal)
        FinishStepRec(next);
    }
  }
  

  void ResultHandler::UpdateResultRec( ResultContext& actContext ) {

    assert( actContext.postProcs.GetSize() == actContext.nextContexts.GetSize() );

    // iterate over all postprocs 
    for( UInt i = 0; i < actContext.postProcs.GetSize(); i++ ) {
      LOG_DBG(resHandler) << "Applying postProc '" << actContext.postProcs[i]->GetName()
                          << "' on '" << actContext.nextContexts[i]->result->GetEntityList()->GetName();
                
      actContext.postProcs[i]->Apply();

      if( !actContext.isFinal ) {

        // call method recursively for all following contexts,
        // if current result is not reducing in time (=final).

        // If, for example, a sum is applied with reduction over
        // time and the maximum (spatial) is to be calculated of 
        // this sum, then the calculation of the maximum may only
        // be performed in the final stage (isFinal == true);
        UpdateResultRec( *(actContext.nextContexts[i]) );
      }
    }
      
  }


  void ResultHandler::AddOutputDest(shared_ptr<SimOutput> out, const std::string& id, const std::string& gridId )
  {
    LOG_DBG(resHandler) << "Adding output writer with id '" << id << "', on grid with id '" << gridId << "'";
    simOutputHandlers_[id] = out;
    outGridIds_[id] = gridId;
  }

  void ResultHandler::Init(std::map<std::string, Grid* >& gridMap, bool printGridOnly ) {
    
    // Trigger function with all output files
    for(auto it = simOutputHandlers_.begin(); it != simOutputHandlers_.end(); it++ )
    {
      // skip if only streaming is set and this is no streamingHandler
      if (streamOnly && !it->second->IsStreaming())
        continue;

      // get gridId for this output writer
      std::string gridId = outGridIds_[it->first];
      
      // check, if grid is defined
      if( gridMap.find(gridId) == gridMap.end() ) {
        EXCEPTION( "No grid with id '" << gridId << "' defined" );
      }
      it->second->Init( gridMap[gridId], printGridOnly );
    }
  }
  

  void ResultHandler::Finalize() {
        
    LOG_DBG(resHandler) << "Starting to Finalize";
    
    // Trigger writing and finalizing of all output writers
    for(auto fileIt = simOutputHandlers_.begin();
    fileIt != simOutputHandlers_.end(); fileIt++ ) {
      // skip if only streaming is set and this is no streamingHandler
      if (streamOnly && !fileIt->second->IsStreaming())
        continue;

      LOG_DBG(resHandler) << "Finalizing result for output '" << fileIt->second->GetName() << "'";
      fileIt->second->Finalize( );
    }    

    LOG_DBG(resHandler) << "Finished Finalizing" << std::endl;
  }

  void ResultHandler::FinishMultiSequenceStepRec( ResultContext& actContext ) {

    LOG_DBG(resHandler) << "Starting to finish MsStep  recursively";
    assert( actContext.postProcs.GetSize() == actContext.nextContexts.GetSize() );
    
    std::string type;

    // iterate over all postProcs and trigger finalization
    for(UInt i = 0 ; i < actContext.postProcs.GetSize(); i++)
      actContext.postProcs[i]->Finalize();

    // iterate over all contexts
    for( UInt i = 0; i < actContext.nextContexts.GetSize(); i++ ) {
      
      ResultContext & next = *(actContext.nextContexts[i]);
      
      // only print result, if it is defined for final stage
      if( next.isFinal ) {
        
        // Only write result, if flag is set
        if( next.writeResult ) {

          // iterate over all outputs
          for( UInt iOut = 0; iOut < next.outputIds.GetSize(); iOut++ ) {

            // skip if only streaming is set and this is no streamingHandler
            if (streamOnly && !simOutputHandlers_[next.outputIds[iOut]]->IsStreaming())
              continue;

            // Add current result to given output file
            simOutputHandlers_[next.outputIds[iOut]]->AddResult( next.result );
            BaseResult & nextResult  = *(next.result);
            LOG_DBG(resHandler) << "Adding final result '" << nextResult.GetResultInfo()->resultName
                                << "' on '" << nextResult.GetEntityList()->GetName()
                                << "' to '" << simOutputHandlers_[next.outputIds[iOut]]->GetName() << "'";
          }
        }
        
      } // isFinal

      // Call FinalizeRec recursively for each next-context in any case
      FinishMultiSequenceStepRec( next );
    }
    LOG_DBG(resHandler) << "Finished MsStep recursively";
  }
  

  bool ResultHandler::IsOutput(  ResultContext& context )
  {
    if( actStep_ < context.saveBegin|| actStep_ > context.saveEnd )
      return false;

    if( ((actStep_-context.saveBegin) % context.saveInc) == 0 ) {
      return true;
    }
    
    return false;
  }

  void ResultHandler::Dump()
  {
    std::cout << "ResultHandler[sequenceStep=" << sequenceStep_ 
              << "; actStep=" << actStep_ << " actStepVal=" << actStepVal_ 
              << "; needed=" << isNeeded_.size() << "; updated=" << isUpdated_.size()
              << std::endl;
    
    std::set<shared_ptr<BaseResult> >::iterator iter;
    
    for(iter = isNeeded_.begin(); iter !=  isNeeded_.end(); iter++)
      std::cout << " needed: " << (*iter)->ToString() << std::endl;       

    for(iter = isUpdated_.begin(); iter !=  isUpdated_.end(); iter++)
      std::cout << " updated: " << (*iter)->ToString() << std::endl;       
  }

  void ResultHandler::GetDefaultOutputs( shared_ptr<BaseResult> res,
                                         StdVector<std::string>& out ) {
    
    out.Clear();
      
    // Check, on which type of results the input is defined
    SimOutput::Capability neededCap = SimOutput::NONE;
    ResultInfo & actInfo = *( res->GetResultInfo() );
    ResultInfo::EntityUnknownType definedOn = actInfo.definedOn;
    EntityList::DefineType definedBy = res->GetEntityList()->GetDefineType();
    
    
    // a) Mesh result (Nodes/Elems/Surfelems by region)
    if( definedOn == ResultInfo::NODE ||
        definedOn == ResultInfo::ELEMENT ||
        definedOn == ResultInfo::SURF_ELEM ) {
      if( definedBy == EntityList::REGION ) {
        neededCap = SimOutput::MESH_RESULTS;
      } else {
        neededCap = SimOutput::HISTORY;
      }
    } 
    // b) Named-Node/Elem/Region/Coil-Result (nodes/elems/surfelems by name,
    else {
      neededCap = SimOutput::HISTORY;
    }
    
    // Check all defined output writers for the searched capability
    for(auto it = simOutputHandlers_.begin(); it != simOutputHandlers_.end(); it++)
      if(it->second->GetCapabilities().find(neededCap) != it->second->GetCapabilities().end())
        out.Push_back(it->first);
    
    // Check, if any suitable output class could be found
    if( out.GetSize() == 0 )
    {
      // If no writer with capability MESH_RESULTS is found: error
      if( neededCap == SimOutput::MESH_RESULTS ) {
        EXCEPTION( "No output class was specified, which is capable of "
                   << "writing mesh results. Please specify one within the output section!" );
      }

      // If no writer with capability HISTORY_RESULT is found: 
      // -> Create default text output (a.k.a. standard history-writer)
      
      if( neededCap == SimOutput::HISTORY ) {
        std::string simName = progOpts->GetSimName();
        // check for restart
        bool restart = progOpts->GetRestart();
        shared_ptr<SimOutput> textOut 
          = shared_ptr<SimOutput>(new SimOutputText( simName, PtrParamNode(), 
                                                     PtrParamNode(), restart ) );
        textOut->Init(  domain->GetGrid(), false );
        simOutputHandlers_["histDefault"] = textOut;
        outGridIds_["histDefault"] = "default";
        out.Push_back( "histDefault" );
        
        // Register multisequence step also for the newly created output class
        // if it was already begun
        if( analysisType_ != BasePDE::NO_ANALYSIS ) {
          textOut->BeginMultiSequenceStep( sequenceStep_, analysisType_, numSteps_ );
        }
      }
    }
  }

  void ResultHandler::AddInputReader( shared_ptr<SimInput> inClass, 
                                      const std::string& readerId ) {

    LOG_DBG(resHandler) << "Adding input reader with id '" << readerId << "'";
    inFiles_[readerId] = inClass;

  }

  shared_ptr<SimInput> ResultHandler::GetInputReader( const std::string& readerId ) {
    
    // check, if input read with specified id is present
    if( inFiles_.find(readerId) == inFiles_.end() ) {
      EXCEPTION( "Input reader with id '" << readerId << "' does not exist")
    }
    return inFiles_[readerId];
  }



  void ResultHandler::
  GetNumMultiSequenceSteps( const std::string& readerId,
                            std::map<UInt, BasePDE::AnalysisType>& analysis,
                            std::map<UInt, UInt>& numSteps,
                            bool isHistory ) { 
    
    // check, if input reader exists
    if( inFiles_.find(readerId) == inFiles_.end() ) {
      EXCEPTION( "Input reader with id '" << readerId 
                 << "' is not registered yet" );
    }
    
    inFiles_[readerId]->GetNumMultiSequenceSteps( analysis, numSteps, 
                                                  isHistory );
  }

  void ResultHandler::
  GetResultTypes( const std::string& readerId,
                    UInt sequenceStep,
                    StdVector<shared_ptr<ResultInfo> >& infos,
                    bool isHistory ) {

    // check, if input reader exists
    if( inFiles_.find(readerId) == inFiles_.end() ) {
      EXCEPTION( "Input reader with id '" << readerId 
                 << "' is not registered yet" );
    }

    inFiles_[readerId]->GetResultTypes( sequenceStep,
                                        infos, isHistory );
  }
    
  void ResultHandler::GetStepValues(const std::string& readerId, UInt sequenceStep,
                                    shared_ptr<ResultInfo> info, std::map<UInt, Double>& steps, bool isHistory)
  {
    // check, if input reader exists
    if(inFiles_.find(readerId) == inFiles_.end())
      EXCEPTION( "Input reader with id '" << readerId<< "' is not registered yet" );
    
    inFiles_[readerId]->GetStepValues(sequenceStep, info, steps, isHistory);
  }
  
  void ResultHandler::GetResultEntities(const std::string& readerId, UInt sequenceStep, shared_ptr<ResultInfo> info,
                                       StdVector<shared_ptr<EntityList> >& list, bool isHistory)
  {
    // check, if input reader exists
    if(inFiles_.find(readerId) == inFiles_.end())
      EXCEPTION( "Input reader with id '" << readerId << "' is not registered yet" );

    inFiles_[readerId]->GetResultEntities( sequenceStep, info, list, isHistory);
  }

  void ResultHandler::GetResult(const std::string& readerId, UInt sequenceStep, UInt stepValue,
                                shared_ptr<BaseResult> result, bool isHistory)
  {
    // check, if input reader exists
    if(inFiles_.find(readerId) == inFiles_.end())
      EXCEPTION("Input reader with id '" << readerId << "' is not registered yet");
    
    inFiles_[readerId]->GetResult(sequenceStep, stepValue, result, isHistory);
  }
 
  shared_ptr<BaseResult> ResultHandler::GetResult(const std::string& readerId, UInt sequenceStep,
                                                  UInt stepValue, SolutionType solType, const std::string& regionName)
  {
    shared_ptr<BaseResult> result;

    result = inFiles_[readerId]->GetResult(sequenceStep, stepValue, solType, regionName);
    return result;
  }
  
  
  template<typename TYPE>
  shared_ptr<FeFunction<TYPE> > ResultHandler::GetFeFunction(const std::string& readerId, UInt sequenceStep, UInt stepValue,
                                     SolutionType solType, std::set<std::string> & regionNames, PtrParamNode rootNode)
  {
    return inFiles_[readerId]->GetFeFunction<TYPE>(sequenceStep, stepValue, solType, regionNames);
  }
  
  
  
  void ResultHandler::InterpolateRes( ResultContext& actContext, Grid* destGrid, shared_ptr<BaseResult>& res)
  {
    Grid* srcGrid = domain->GetGrid( "default" );
    std::string entListName = actContext.result->GetEntityList()->GetName();
    EntityList::ListType lType = actContext.result->GetEntityList()->GetType();
    shared_ptr<EntityList> destList = destGrid->GetEntityList( lType, entListName );
    
    StdVector<shared_ptr<EntityList> >lists(1);
    lists[0] = srcGrid->GetEntityList(EntityList::ELEM_LIST, entListName); 
 
    // loop over all elements, get the element midpoint and store it in a vector
    StdVector<Vector<Double> > globPoints(destList->GetSize());
    EntityIterator it = destList->GetIterator();
    if( actContext.result->GetResultInfo()->definedOn == ResultInfo::NODE ) {
      // loop over all nodes
      UInt pos = 0;
      for( ; !it.IsEnd(); it++ ) {
        // get global element midpoint
        UInt nodeNum = it.GetNode();
        destGrid->GetNodeCoordinate( globPoints[pos++], nodeNum, false );
      }
    } else if( actContext.result->GetResultInfo()->definedOn == ResultInfo::ELEMENT ) {
      // loop over all elements
      shared_ptr<ElemShapeMap> esm;
      UInt pos = 0;
      for( ; !it.IsEnd(); it++ ) {

        // get global element midpoint
        const Elem* ptEl = it.GetElem();
        esm = destGrid->GetElemShapeMap( ptEl , false );
        esm->GetGlobMidPoint( globPoints[pos++] );
      } 
    }else {
      // In this case we have results not being defined on nodes or elements, 
      // which should not happen
    }

    // now map global locations to our src grid
    StdVector<LocPoint> lps(destList->GetSize());
    StdVector<const Elem*> elems(destList->GetSize());
    srcGrid->GetElemsAtGlobalCoords( globPoints, lps, elems, lists);


    // Get result functor
    shared_ptr<ResultFunctor> fct = actContext.functor;

    // Create result vector (depending on type)
    if(actContext.result->GetEntryType() == BaseMatrix::DOUBLE ) {

      PtrCoefFct coef;
      // try to get the coefficient function
      if( typeid(*fct) == typeid(FieldCoefFunctor<Double>))
        coef = (dynamic_pointer_cast<FieldCoefFunctor<Double> >(fct))->GetCoefFct();
      else
        EXCEPTION( "Can not interpolate values, not defined on a field");

      res.reset(new Result<Double>());
      Vector<Double> & vals = dynamic_cast<Vector<Double>& >(*res->GetSingleVector());
      UInt numDofs = actContext.result->GetResultInfo()->dofNames.GetSize();
      vals.Resize( lps.GetSize() * numDofs );
      vals.Init(0.0);

      // Loop over all local points
      UInt pos = 0;
      Vector<Double> temp;
      shared_ptr<ElemShapeMap> esm;
      LocPointMapped lpm;
      if( coef->GetDimType() == CoefFunction::SCALAR ) {
        temp.Resize(1);
        for( UInt i = 0; i < lps.GetSize(); ++i ) {
          if( elems[i] == NULL) {
            vals[pos++] = NAN;
            continue;
          }
          esm = srcGrid->GetElemShapeMap( elems[i], true );
          lpm.Set( lps[i], esm, 0.0 );
          coef->GetScalar( temp[0], lpm);
          vals[pos++] = temp[0];
        } // loop over elements

      } else if( coef->GetDimType() == CoefFunction::VECTOR ) {
        for( UInt i = 0; i < lps.GetSize(); ++i ) {
          if( elems[i] == NULL) {
            vals[pos++] = NAN;
            continue;
          }
          esm = srcGrid->GetElemShapeMap( elems[i], true );
          lpm.Set( lps[i], esm, 0.0 );

          // Calculate coefficient function and store the result into the vector
          coef->GetVector( temp, lpm);

          // loop over all dofs and save the result
          for( UInt iDof=0; iDof < numDofs; ++iDof ) {
            vals[pos++] = temp[iDof];
          } // loop over dofs
        } // loop over elements
      }
     
    }  else {
      PtrCoefFct coef;
      // try to get the coefficient function
      if( typeid(*fct) == typeid(FieldCoefFunctor<Complex>))
        coef = (dynamic_pointer_cast<FieldCoefFunctor<Complex> >(fct))->GetCoefFct();
      else
        EXCEPTION("Can not interpolate values, not defined on a field");

      res.reset(new Result<Complex>());
      Vector<Complex> & vals = dynamic_cast<Vector<Complex>& >(*res->GetSingleVector());
      UInt numDofs = actContext.result->GetResultInfo()->dofNames.GetSize();
      vals.Resize( lps.GetSize() * numDofs );
      vals.Init(0.0);

      // Loop over all local points
      UInt pos = 0;
      Vector<Complex> temp;
      shared_ptr<ElemShapeMap> esm;
      LocPointMapped lpm;
      if( coef->GetDimType() == CoefFunction::SCALAR ) {
        temp.Resize(1);
        for( UInt i = 0; i < lps.GetSize(); ++i ) {
          if( elems[i] == NULL)
            continue;
          esm = srcGrid->GetElemShapeMap( elems[i], true );
          lpm.Set( lps[i], esm, 0.0 );
          coef->GetScalar( temp[0], lpm);
          vals[pos++] = temp[0];
        } // loop over elements

      } else if( coef->GetDimType() == CoefFunction::VECTOR ) {
        for( UInt i = 0; i < lps.GetSize(); ++i ) {
          if( elems[i] == NULL)
            continue;
          esm = srcGrid->GetElemShapeMap( elems[i], true );
          lpm.Set( lps[i], esm, 0.0 );

          // Calculate coefficient function and store the result into the vector
          coef->GetVector( temp, lpm);

          // loop over all dofs and save the result
          for( UInt iDof=0; iDof < numDofs; ++iDof ) {
            vals[pos++] = temp[iDof];
          } // loop over dofs
        } // loop over elements
      } // vector / scalar case
    } // complex case
    
    res->SetEntityList(destList);
    res->SetResultInfo(actContext.result->GetResultInfo());
  }


// Explicit template instantiation
  template
  shared_ptr<FeFunction<Double> >
  ResultHandler::GetFeFunction<Double>( const std::string& readerId,
      UInt sequenceStep,
      UInt stepValue,
      SolutionType solType,
      std::set<std::string> & regionNames,
      PtrParamNode rootNode );
  template
  shared_ptr<FeFunction<Complex> >
  ResultHandler::GetFeFunction<Complex>( const std::string& readerId,
      UInt sequenceStep,
      UInt stepValue,
      SolutionType solType,
      std::set<std::string> & regionNames,
      PtrParamNode rootNode );
}
