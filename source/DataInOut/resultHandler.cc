// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "resultHandler.hh"

#include "DataInOut/CommandLine/BaseCommandLineHandler.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/Logging/cfslog.hh"
#include "Domain/domain.hh"
#include "DataInOut/postProc.hh"
#include "DataInOut/SimInOut/TextOutput/textSimOutput.hh"

namespace CoupledField {


  // declare logging stream
  DECLARE_LOG(resHandler)
  DEFINE_LOG(resHandler, "resultHandler")
    
    ResultHandler::ResultHandler( OpMode opMode) {
    ENTER_FCN( "ResultHandler::ResultHandler" );
    
    opMode_ = opMode;
    sequenceStep_ = 1;
    actStep_ = 0;
    finalResultExists_ = false;
    actStepVal_ = 0.0;
  }


  ResultHandler::~ResultHandler() {
    ENTER_FCN( "ResultHandler::~ResultHandler" );

    // not much to do here
  }

  void ResultHandler::
  RegisterResult( shared_ptr<BaseResult> sol,
                  UInt saveBegin, UInt saveInc,
                  UInt saveEnd, 
                  const StdVector<std::string> & outDestNames,
                  const std::string& postProcName,
                  bool writeResult ) {
    ENTER_FCN( "ResultHandler::RegisterResult" );

    ResultInfo & actDof = *(sol->GetResultInfo());

    // From here on, the entry in the ResultInfo.resultName
    // is needed instead of the Resultinfo.resultType
    std::string name;
    Enum2String( actDof.resultType, name );
    if( actDof.resultName == "" ) {
      actDof.resultName = name;
    }
    
    LOG_DBG(resHandler) << "Registering result:";
    LOG_DBG(resHandler) << "-------------------";
    LOG_DBG(resHandler) << "name: " << name;
    LOG_DBG(resHandler) << "dofs: " << actDof.dofNames.Serialize();
    LOG_DBG(resHandler) << "resultList: " << sol->GetEntityList()->GetName();
    LOG_DBG(resHandler) << "saveBegin: " << saveBegin;
    LOG_DBG(resHandler) << "saveEnd: " << saveEnd;
    LOG_DBG(resHandler) << "saveInc: " << saveInc;
    LOG_DBG(resHandler) << "writeResult: " << writeResult;
    LOG_DBG(resHandler) << "isFinal: 0";
    LOG_DBG(resHandler) << "outputDest: " << outDestNames.Serialize();
    LOG_DBG(resHandler) << "postProcId: " << postProcName << std::endl;

    // check, if result context was already created
    shared_ptr<ResultContext> actContext;
  
    if( resultContexts_.find( sol) != resultContexts_.end() ) {
      EXCEPTION( "Context was already found" );
    }

    // create new ResultContext
    actContext = 
      shared_ptr<ResultContext>( new ResultContext() );
    actContext->result = sol;
    actContext->saveBegin = saveBegin;
    actContext->saveEnd = saveEnd;
    actContext->saveInc = saveInc;
    actContext->writeResult = writeResult;
    actContext->isFinal = false;
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
        if( outFiles_.find( newDest[i] ) == outFiles_.end() ) {
          EXCEPTION( "Output writer '" << newDest[i] 
                     << "' was not registered yet!" );
        }
        LOG_DBG(resHandler) << "Registering output '" << newDest[i]
                            << "' with result '" << name;
        actContext->outputs.Push_back( outFiles_[newDest[i]] );
        
        // register results also at the output writer class
        outFiles_[newDest[i]]->RegisterResult( sol, saveBegin,
                                               saveInc, saveEnd );
      }
    }

    // insert context to list
    contexts_.insert( actContext );
    
    // Check for postprocessing method(s), which might also
    // produce new result contexts
    if( postProcName != "" ) {
      RegisterResultRec( *actContext, postProcName );
    }
    LOG_TRACE(resHandler) << "Finished registering result" << std::endl;
  }
  
  void ResultHandler::BeginMultiSequenceStep( UInt step,
                                              AnalysisType type,
                                              UInt numSteps ) {
    ENTER_FCN( "ResultHandler::BeginMultiSequenceStep" );
        
    // store current sequencestep
    sequenceStep_ = step;
    
    // Iterate over all outfiles to notify them about a new
    // multisequence step
    std::map<std::string, shared_ptr<SimOutput> >::iterator it;
    it = outFiles_.begin();

    for( ; it != outFiles_.end(); it++ ) {
      it->second->BeginMultiSequenceStep( step, type, numSteps );
    }
  }

  void ResultHandler::BeginStep( UInt stepNum, Double stepVal ) {
    ENTER_FCN( "ResultHandler::BeginStep" );

    LOG_TRACE(resHandler) << "Begin step " << stepNum;

    // remeber current step values 
    actStep_ = stepNum;
    actStepVal_ = stepVal;

    // check all result contexts, which one of them has to be updated
    // (determined by saveBegin, saveEnd and saveInc )
    isNeeded_.clear();
    isUpdated_.clear();
    std::set<shared_ptr<ResultContext> >::iterator contextIt;
    contextIt = contexts_.begin();
    for( ; contextIt != contexts_.end(); contextIt++ ) {
      BaseResult & actResult  = *((*contextIt)->result);
      LOG_DBG(resHandler) << "IsNeeded for '" 
                          << actResult.GetResultInfo()->resultName 
                          << "' on '"
                          << actResult.GetEntityList()->GetName() << "' is '" 
                          << (IsOutput( **contextIt ) == true ? "true" : "false");
      if( IsOutput( **contextIt ) ) {
        isNeeded_.insert( (*contextIt)->result );
      }
    }
    
    // Trigger function for all output files
    std::map<std::string, shared_ptr<SimOutput> >::iterator it;
    for( it = outFiles_.begin(); it != outFiles_.end(); it++ ) {
      it->second->BeginStep( stepNum, stepVal );
    }
    
    LOG_TRACE(resHandler) << "Finished beginning of new step" << std::endl;
  }

  void ResultHandler::UpdateResult( shared_ptr<BaseResult> sol ) {
    ENTER_FCN( "ResultHandler::UpdateResult" );

    LOG_TRACE(resHandler) << "Updating results";
    
    // check, if result is to be updated
    if( isNeeded_.find(sol) == isNeeded_.end() ) {
      *warning << "Result '" 
               << sol->GetResultInfo()->resultName
               << "' on entitylist '"
               << sol->GetEntityList()->GetName() 
               << "' is not needed in step " << actStep_;
      Warning( __FILE__, __LINE__ );
    }

    // Set flag for update to true
    isUpdated_.insert( sol );

    LOG_DBG(resHandler) << "Result '" 
                        << sol->GetResultInfo()->resultName
                        << "' was provided on '"
                        << sol->GetEntityList()->GetName() 
                        << "' in step " << actStep_;

    // Update results recursively (for postprocessing results )
    UpdateResultRec( *(resultContexts_[sol] ) );
    LOG_TRACE(resHandler) << "Finished updating results" << std::endl;
  }


  void ResultHandler::FinishStep( ) {
    ENTER_FCN( "ResultHandler::FinishStep" );
    
    LOG_TRACE(resHandler) << "Starting to finish step " << actStep_;

    // === Primary results ===
    // iterate over all results, which are needed
    std::set<shared_ptr<BaseResult> >::iterator it = isNeeded_.begin();
    for( ; it != isNeeded_.end(); it++ ) {

      // store context
      ResultContext & actContext = *(resultContexts_[*it]);
      BaseResult & actResult  = *(actContext.result);
                   
      LOG_DBG(resHandler) << "Checking result '"
                          << actResult.GetResultInfo()->resultName
                          << "' on '"
                          << actResult.GetEntityList()->GetName()
                          << "' for writing out";

      // Check, if result was updated at all
      if( isUpdated_.find(*it) == isUpdated_.end() ) {
        *warning << "Result '" 
                 << actResult.GetResultInfo()->resultName 
                 << "' on '" 
                 << (*it)->GetEntityList()->GetName()
                 << "' was not provided in step " << actStep_;
        Warning( __FILE__, __LINE__ );
      }


      // check, if result is to be written 
      if( actContext.writeResult && (!actContext.isFinal ) ) {

        // iterate over all outputs
        for( UInt iOut = 0; iOut < actContext.outputs.GetSize(); iOut++ ) {
          
          // Add current result to given output file
          actContext.outputs[iOut]->AddResult( actContext.result );
          LOG_DBG(resHandler) << "Adding result '" 
                              << actResult.GetResultInfo()->resultName 
                              << "' on '"
                              << actResult.GetEntityList()->GetName()
                              << "' to '" 
                              << actContext.outputs[iOut]->GetName()
                              << "'";
        }
      }
      
      // In any case: Process also related secondary (postprocessing) results
      FinishStepRec( actContext );
      
    }
    

    // Trigger writing of all output writers
    std::map<std::string, shared_ptr<SimOutput> >::iterator fileIt;
    for( fileIt = outFiles_.begin(); 
         fileIt != outFiles_.end(); fileIt++ ) {
      LOG_DBG(resHandler) << "Finishing step for output '"
                          << fileIt->second->GetName() << "'";
      fileIt->second->FinishStep( );
    }    

    LOG_TRACE(resHandler) << "Finished step " << actStep_ << std::endl;
  }

  void ResultHandler::FinishMultiSequenceStep() {
    ENTER_FCN( "ResultHandler::FinishMultiSequenceStep" );

    // Iterate over all outfiles to notify them about the end of a
    // multisequence step
    std::map<std::string, shared_ptr<SimOutput> >::iterator it;
    it = outFiles_.begin();

    for( ; it != outFiles_.end(); it++ ) {
      it->second->FinishMultiSequenceStep( );
    }

    // Delete all results, as they are registered newly for each
    // sequence step
    contexts_.clear();
    resultContexts_.clear();
    isNeeded_.clear();
    isUpdated_.clear();
  }

  void ResultHandler::RegisterResultRec( ResultContext& actContext, 
                                         const std::string& postProcName ) {
    ENTER_FCN( "ResultHandler::RegisterResultRec" );

    LOG_DBG(resHandler) << "Registering (recursively) result of postProc '" 
                        << postProcName << "'";

    // if no postProcName is set, just leave
    if( postProcName == "" ) {
      return;
    }

    // fetch current postProcNode
    ParamNode * postProcNode = 
      param->Get("sequenceStep", "index", GenStr(sequenceStep_) )
      ->Get( "postProcList")->Get( "postProc", "id", postProcName );
    
    if( !postProcNode ) {
      EXCEPTION( "A Postprocessing section for '" << postProcName
                 << "' does not exist" );
    }

    // Fetch postprocs
    StdVector<shared_ptr<PostProc> > postProcs;
    PostProc::CreatePostProc( postProcNode,  domain->GetGrid(),
                              postProcs );
                              
    // iterate over all postprocs 
    for( UInt i = 0; i < postProcs.GetSize(); i++ ) {

      // register current solution with new object
      postProcs[i]->SetResult( actContext.result );

      // create new resultcontext and link it to current one
      shared_ptr<ResultContext> nextContext = 
        shared_ptr<ResultContext>( new ResultContext() );
      nextContext->result = postProcs[i]->GetOutputResult();
      nextContext->saveBegin = actContext.saveBegin;
      nextContext->saveEnd = actContext.saveEnd;
      nextContext->saveInc = actContext.saveInc;
      nextContext->writeResult = postProcs[i]->IsWriteResult();
      if( postProcs[i]->GetReductionType() == PostProc::TIME_FREQ ||
          actContext.isFinal ) {
        nextContext->isFinal = true;
        finalResultExists_ = true;
      } else {
        nextContext->isFinal = false;
      }
      
      // fetch suitable output class
      StdVector<std::string> oldDest, outDest;
      postProcs[i]->GetOutDestNames( oldDest );
      LOG_DBG(resHandler) << "outputDest of " << postProcName << " is '" 
                          << oldDest.Serialize() << "'";
      LOG_DBG(resHandler) << "size of outputDest is " << outDest.GetSize() 
                          << std::endl;
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
      LOG_DBG(resHandler) << "outputDest: " << outDest.Serialize();
      LOG_DBG(resHandler) << "postProcName: " 
                          << postProcs[i]->GetNextPostProcName() << std::endl;

      for( UInt iOut = 0; iOut < outDest.GetSize(); iOut++ ) {
        if( outFiles_.find( outDest[iOut] ) == outFiles_.end() ) {
          EXCEPTION( "Output writer '" << outDest[i] 
                     << "' was not registered yet!" );
        }
        
        nextContext->outputs.Push_back( outFiles_[outDest[iOut]] );
      
        // register results also at the output writer class
        outFiles_[outDest[iOut]]->
          RegisterResult(  postProcs[i]->GetOutputResult(),
                           actContext.saveBegin, actContext.saveInc,
                           actContext.saveEnd );
      }
      
      // store postproc and result in current context2
      actContext.nextContexts.Push_back( nextContext );
      actContext.postProcs.Push_back( postProcs[i] );
      
      // call method iteratively
      if( postProcs[i]->GetNextPostProcName() != "" ) {
        RegisterResultRec( *nextContext, 
                           postProcs[i]->GetNextPostProcName() );
      }
    }
    LOG_DBG(resHandler) << "Finished registering result of postProc '"
                        << postProcName << "'";
  }
    
  void ResultHandler::FinishStepRec( ResultContext& actContext ) {
    ENTER_FCN( " ResultHandler::FinishStepRec" );
    
    // ensure that we have for each postprocessing method
    // also the related result context
    assert( actContext.postProcs.GetSize() ==
            actContext.nextContexts.GetSize() );

    // iterate over all contexts
    for( UInt i = 0; i < actContext.nextContexts.GetSize(); i++ ) {

      ResultContext & next = *(actContext.nextContexts[i]);
      LOG_DBG(resHandler) << "Checking result '"
                          << next.result->GetResultInfo()->resultName
                          << "' on '"
                          << next.result->GetEntityList()->GetName()
                          << "' for writing out";

      // Check, if next result is only to be computed in final stage
      if( next.writeResult != false  && next.isFinal != true) {
        
        // iterate over all outputs
        for( UInt iOut = 0; iOut < next.outputs.GetSize(); iOut++ ) {
          
          // Add current result to given output file
          next.outputs[iOut]->AddResult( next.result );
          BaseResult & nextResult  = *(next.result);
          LOG_DBG(resHandler) << "Adding result '" 
                              << nextResult.GetResultInfo()->resultName
                              << "' on '"
                              << nextResult.GetEntityList()->GetName()
                              << "' to '" << next.outputs[iOut]->GetName()
                              << "'";
        }
      }

      // Call method recursively to write all outputs of next context
      // Note: We only may write the result, if it is not a "final" result,
      // i.e. it is written only once at the end. In this case, the result
      // is written via the method FinalizeRec()
      if( !next.isFinal ) {
        FinishStepRec( next );
      }
    }
  }
  

  void ResultHandler::UpdateResultRec( ResultContext& actContext ) {
    ENTER_FCN( "ResultHandler::UpdateResultRec" );

    assert( actContext.postProcs.GetSize() ==
            actContext.nextContexts.GetSize() );

    // iterate over all postprocs 
    for( UInt i = 0; i < actContext.postProcs.GetSize(); i++ ) {
      LOG_DBG(resHandler) << "Applying postProc '" 
                          << actContext.postProcs[i]->GetName() 
                          << "' on '" 
                          << actContext.nextContexts[i]->result
        ->GetEntityList()->GetName();
                
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


  void ResultHandler:: AddOutputDest( shared_ptr<SimOutput> out, 
                                      const std::string& id ) {
    ENTER_FCN( "ResultHandler::AddOutputDest" );

    LOG_DBG(resHandler) << "Adding output writer with id ";
    outFiles_[id] = out;

  }

  void ResultHandler::Init( Grid* ptGrid, bool printGridOnly ) {
    ENTER_FCN( "ResultHandler::Init" );
    
    // Trigger function with all output files
    std::map<std::string, shared_ptr<SimOutput> >::iterator it;
    for( it = outFiles_.begin(); it != outFiles_.end(); it++ ) {
      it->second->Init( ptGrid, printGridOnly );
    }
  }
  

  void ResultHandler::Finalize() {
    ENTER_FCN( "ResultHandler::Finalize" );
        
    LOG_TRACE(resHandler) << "Starting to Finalize";
    actStep_++;
    actStepVal_ = 0.0;

    std::map<std::string, shared_ptr<SimOutput> >::iterator fileIt;
    
    // Note: At the moment, the 'finalization' is done by
    // writing an aditional last step to the file (with stepVal = 0 );
    // However, this is only done, if an additional result is present
    // for the final step.
    if( finalResultExists_ ) {
      
      std::string type;
      // Trigger new step at all output writers
     
      for( fileIt = outFiles_.begin(); 
           fileIt != outFiles_.end(); fileIt++ ) {
        fileIt->second->BeginStep( actStep_, actStepVal_);
      }    
      
      // === Primary results ===
      // iterate over all results, which are needed
      std::set<shared_ptr<ResultContext> >::iterator it = contexts_.begin();
      for( ; it != contexts_.end(); it++ ) {
        
        // store context
        ResultContext & actContext = **it;
        BaseResult & actResult  = *(actContext.result);
        Enum2String( actResult.GetResultInfo()->resultType, type );
        
        LOG_DBG(resHandler) << "Testing result '" << type << "' on '"
                            << actResult.GetEntityList()->GetName()
                            << "'";
        // check, if result is to be written 
        if( actContext.writeResult &&  actContext.isFinal )  {
          
          // iterate over all outputs
          for( UInt iOut = 0; iOut < actContext.outputs.GetSize(); iOut++ ) {
            
            // Add current result to given output file
            actContext.outputs[iOut]->AddResult( actContext.result );
            LOG_DBG(resHandler) << "Adding final result '" << type << "' on '"
                                << actResult.GetEntityList()->GetName()
                                << "' to '" << actContext.outputs[iOut]->GetName()
                                << "'";
          }
        }
        
        // In any case: Process also related secondary (postprocessing) results
        FinalizeRec( actContext );
      }
    }
    
    // Trigger writing and finalizing of all output writers
    for( fileIt = outFiles_.begin(); 
         fileIt != outFiles_.end(); fileIt++ ) {
      LOG_DBG(resHandler) << "Finalizing result for output '"
                          << fileIt->second->GetName() << "'";
      if( finalResultExists_) {
        fileIt->second->FinishStep( );
      }
      fileIt->second->Finalize( );
    }    

    LOG_TRACE(resHandler) << "Finshed Finalizing" << std::endl;
  }

  void ResultHandler::FinalizeRec( ResultContext& actContext ) {
    ENTER_FCN( "ResultHandler::Finalize" );

    assert( actContext.postProcs.GetSize() ==
            actContext.nextContexts.GetSize() );
    
    std::string type;

    // iterate over all postProcs and trigger finalization
    for( UInt i = 0 ; i < actContext.postProcs.GetSize(); i++ ) {
      actContext.postProcs[i]->Finalize();
    }

    // iterate over all contexts
    for( UInt i = 0; i < actContext.nextContexts.GetSize(); i++ ) {
      
      ResultContext & next = *(actContext.nextContexts[i]);
      
      // only print result, if it is defined for final stage
      if( next.isFinal ) {
        
        // Only write result, if flag is set
        if( next.writeResult ) {

          // iterate over all outputs
          for( UInt iOut = 0; iOut < next.outputs.GetSize(); iOut++ ) {
            
            // Add current result to given output file
            next.outputs[iOut]->AddResult( next.result );
            BaseResult & nextResult  = *(next.result);
            LOG_DBG(resHandler) << "Adding final result '" 
                                << nextResult.GetResultInfo()->resultName
                                << "' on '"
                                << nextResult.GetEntityList()->GetName()
                                << "' to '" << next.outputs[iOut]->GetName()
                                << "'";
          }
        }
        
      } // isFinal

      // Call FinalizeRec recursively for each next-context in any case
      FinalizeRec( next );
    }
  }
  

  bool ResultHandler::IsOutput(  ResultContext& context ) {
    if( actStep_ < context.saveBegin||
        actStep_ > context.saveEnd ) 
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
    ENTER_FCN( "ResultHandler::GetDefaultOutput" );
    
    out.Clear();
      
    // Check, on which type of results the input is defined
    SimOutput::Capability neededCap = SimOutput::NONE;
    ResultInfo & actInfo = *( res->GetResultInfo() );
    ResultInfo::EntityUnknownType definedOn = actInfo.definedOn;
    EntityList::DefineType definedBy = res->GetEntityList()->GetDefineType();
    
    
    // a) Mesh result (Nodes/Elems/Surfelems by region)
    if( definedOn == ResultInfo::PFEM || 
        definedOn == ResultInfo::NODE ||
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
    std::map<std::string, shared_ptr<SimOutput> >::iterator it;
    it = outFiles_.begin();
    for( ; it != outFiles_.end(); it++ ) {
      if( it->second->GetCapabilities().find( neededCap ) 
          != it->second->GetCapabilities().end() ) {
        out.Push_back( it->first );
      }
    }
    
    // Check, if any suitable output class could be found
    if( out.GetSize() == 0 ) {

      // If no writer with capability MESH_RESULTS is found: error
      if( neededCap == SimOutput::MESH_RESULTS ) {
        EXCEPTION( "No output class was specified, which is capable of "
                   << "writing mesh results. Please specify one within the "
                   << "<output> section!" );
      }

      // If no writer with capability HISTORY_RESULT is found: 
      // -> Create default text output (a.k.a. standard history-writer)
      
      if( neededCap == SimOutput::HISTORY ) {
        std::string simName = commandLine->GetSimName();
        shared_ptr<SimOutput> textOut 
          = shared_ptr<SimOutput>(new SimOutputText( simName, NULL ) );
        textOut->Init(  domain->GetGrid(), false );
        outFiles_["histDefault"] = textOut;
        out.Push_back( "histDefault" );
      }
    }
  }

  void ResultHandler::AddInputReader( shared_ptr<SimInput> inClass, const std::string& readerId ) {
    ENTER_FCN( "ResultHandler::AddInputReader" );

    LOG_DBG(resHandler) << "Adding input reader with id ";
    inFiles_[readerId] = inClass;
  }
  
}
