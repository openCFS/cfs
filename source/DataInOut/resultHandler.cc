// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <assert.h>
#include <iostream>
#include <utility>

#include "DataInOut/Logging/cfslog.hh"
#include "DataInOut/Logging/log.hpp"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/SimInOut/TextOutput/textSimOutput.hh"
#include "DataInOut/postProc.hh"
#include "DataInOut/programOptions.hh"
#include "DataInOut/simInput.hh"
#include "DataInOut/simOutput.hh"
#include "Domain/domain.hh"
#include "Domain/entityList.hh"
#include "Domain/grid.hh"
#include "Domain/resultInfo.hh"
#include "General/Enum.hh"
#include "General/exception.hh"
#include "MatVec/vector.hh"
#include "PDE/eqnMap.hh"
#include "Utils/nodestoresol.hh"
#include "Utils/result.hh"
#include "resultHandler.hh"

namespace CoupledField {


  // declare logging stream
  DECLARE_LOG(resHandler)
  DEFINE_LOG(resHandler, "resultHandler")
    
    ResultHandler::ResultHandler( OpMode opMode) {
    
    opMode_ = opMode;
    sequenceStep_ = 1;
    actStep_ = 0;
    finalResultExists_ = false;
    actStepVal_ = 0.0;
  }


  ResultHandler::~ResultHandler() {

    // not much to do here
  }

  void ResultHandler::
  RegisterResult( shared_ptr<BaseResult> sol,
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
        if( outFiles_.find( newDest[i] ) == outFiles_.end() ) {
          EXCEPTION( "Output writer '" << newDest[i] 
                     << "' was not registered yet!" );
        }
        LOG_DBG(resHandler) << "Registering output '" << newDest[i]
                            << "' with result '" << actDof.resultName;
        actContext->outputs.Push_back( outFiles_[newDest[i]] );
        
        // register results also at the output writer class
        outFiles_[newDest[i]]->RegisterResult( sol, saveBegin,
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
    LOG_TRACE(resHandler) << "Finished registering result" << std::endl;
  }
  
  void ResultHandler::BeginMultiSequenceStep( UInt step,
                                              BasePDE::AnalysisType type,
                                              UInt numSteps ) {
        
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

    LOG_TRACE(resHandler) << "Updating results";
    
    // check, if result is to be updated
    if( isNeeded_.find(sol) == isNeeded_.end() ) {
      WARN("Result '" 
           << sol->GetResultInfo()->resultName
           << "' on entitylist '"
           << sol->GetEntityList()->GetName() 
           << "' is not needed in step " << actStep_);
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
        WARN("Result '" 
             << actResult.GetResultInfo()->resultName 
             << "' on '" 
             << (*it)->GetEntityList()->GetName()
             << "' was not provided in step " << actStep_);
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
        LOG_DBG(resHandler) << "Testing result '" 
        << actResult.GetResultInfo()->resultName << "' on '"
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
        FinishMultiSequenceStepRec( actContext );
      }
      // Finish newly created step
      for( fileIt = outFiles_.begin(); 
      fileIt != outFiles_.end(); fileIt++ ) {
        fileIt->second->FinishStep();
      }    
    }
    
    
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

    LOG_DBG(resHandler) << "Registering (recursively) result of postProc '" 
                        << postProcName << "'";

    // if no postProcName is set, just leave
    if( postProcName == "" ) {
      return;
    }

    // fetch current postProcNode
    PtrParamNode postProcNode = 
      param->GetByVal("sequenceStep", std::string("index"), actContext.sequenceStep)
      ->Get( "postProcList")->GetByVal( "postProc", "id", postProcName );
    
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
      nextContext->sequenceStep = actContext.sequenceStep;
      nextContext->saveBegin = actContext.saveBegin;
      nextContext->saveEnd = actContext.saveEnd;
      nextContext->saveInc = actContext.saveInc;
      if( postProcs[i]->IsHistory() ||
          actContext.isHistory ) {
        nextContext->isHistory = true;
      } else {
        nextContext->isHistory = false;
      }
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
      LOG_DBG(resHandler) << "isHistory: " << nextContext->isHistory;
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
        if( nextContext->writeResult ) {
          outFiles_[outDest[iOut]]->
          RegisterResult(  postProcs[i]->GetOutputResult(),
                           nextContext->saveBegin, nextContext->saveInc,
                           nextContext->saveEnd,
                           nextContext->isHistory );
        }
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

    LOG_DBG(resHandler) << "Adding output writer with id ";
    outFiles_[id] = out;

  }

  void ResultHandler::Init( Grid* ptGrid, bool printGridOnly ) {
    
    // Trigger function with all output files
    std::map<std::string, shared_ptr<SimOutput> >::iterator it;
    for( it = outFiles_.begin(); it != outFiles_.end(); it++ ) {
      it->second->Init( ptGrid, printGridOnly );
    }
  }
  

  void ResultHandler::Finalize() {
        
    LOG_TRACE(resHandler) << "Starting to Finalize";
    
    // Trigger writing and finalizing of all output writers
    std::map<std::string, shared_ptr<SimOutput> >::iterator fileIt;
    for( fileIt = outFiles_.begin(); 
    fileIt != outFiles_.end(); fileIt++ ) {
      LOG_DBG(resHandler) << "Finalizing result for output '"
      << fileIt->second->GetName() << "'";
      fileIt->second->Finalize( );
    }    

    LOG_TRACE(resHandler) << "Finished Finalizing" << std::endl;
  }

  void ResultHandler::FinishMultiSequenceStepRec( ResultContext& actContext ) {

    LOG_TRACE(resHandler) << "Starting to finish MsStep  recursively";
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
      FinishMultiSequenceStepRec( next );
    }
    LOG_TRACE(resHandler) << "Finished MsStep recursively";
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
                   << "output section!" );
      }

      // If no writer with capability HISTORY_RESULT is found: 
      // -> Create default text output (a.k.a. standard history-writer)
      
      if( neededCap == SimOutput::HISTORY ) {
        std::string simName = progOpts->GetSimName();
        shared_ptr<SimOutput> textOut 
          = shared_ptr<SimOutput>(new SimOutputText( simName, PtrParamNode() ) );
        textOut->Init(  domain->GetGrid(), false );
        outFiles_["histDefault"] = textOut;
        out.Push_back( "histDefault" );
      }
    }
  }

  void ResultHandler::AddInputReader( shared_ptr<SimInput> inClass, 
                                      const std::string& readerId ) {

    LOG_DBG(resHandler) << "Adding input reader with id ";
    inFiles_[readerId] = inClass;

  }

  shared_ptr<SimInput> ResultHandler::GetInputReader(const std::string& readerId, bool silent)
  {
    // check, if input read with specified id is present
    if(inFiles_.find(readerId) == inFiles_.end())
    {
      if(silent)
        return shared_ptr<SimInput>();
      else
        EXCEPTION( "Input reader with id '" << readerId << "' does not exist")
    }
    return inFiles_[readerId];
  }

  shared_ptr<SimOutput> ResultHandler::GetOutputWriter(const std::string& readerId, bool silent)
  {
    // check, if input read with specified id is present
    if(outFiles_.find(readerId) == outFiles_.end())
    {
      if(silent)
        return shared_ptr<SimOutput>();
      else
        EXCEPTION( "Input reader with id '" << readerId << "' does not exist")
    }
    return outFiles_[readerId];
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
    
  void ResultHandler::
  GetStepValues( const std::string& readerId,
                 UInt sequenceStep,
                 shared_ptr<ResultInfo> info,
                 std::map<UInt, Double>& steps,
                 bool isHistory ) {
    // check, if input reader exists
    if( inFiles_.find(readerId) == inFiles_.end() ) {
      EXCEPTION( "Input reader with id '" << readerId 
          << "' is not registered yet" );
    }
    
    inFiles_[readerId]->GetStepValues( sequenceStep, info,
                                       steps, isHistory );
  }
  
  void ResultHandler::
  GetResultEntities( const std::string& readerId,
                     UInt sequenceStep,
                     shared_ptr<ResultInfo> info,
                     StdVector<shared_ptr<EntityList> >& list,
                     bool isHistory ) {
    // check, if input reader exists
    if( inFiles_.find(readerId) == inFiles_.end() ) {
      EXCEPTION( "Input reader with id '" << readerId 
                 << "' is not registered yet" );
    }

    inFiles_[readerId]
      ->GetResultEntities( sequenceStep, info, list, isHistory );
  }

  void ResultHandler::
  GetResult( const std::string& readerId,
             UInt sequenceStep,
             UInt stepValue,
             shared_ptr<BaseResult> result,
             bool isHistory ) {

    // check, if input reader exists
    if( inFiles_.find(readerId) == inFiles_.end() ) {
      EXCEPTION( "Input reader with id '" << readerId 
                 << "' is not registered yet" );
    }
    
    inFiles_[readerId]->GetResult( sequenceStep, stepValue, 
                                   result, isHistory );
  }
 
  shared_ptr<BaseResult> ResultHandler::
  GetResult( const std::string& readerId,
             UInt sequenceStep,
             UInt stepValue,
             SolutionType solType,
             const std::string& regionName ) {

    // check, if input reader exists
    if( inFiles_.find(readerId) == inFiles_.end() ) {
      EXCEPTION( "Input reader with id '" << readerId 
                 << "' is not registered yet" );
    }
    
    // aquire input reader
    shared_ptr<SimInput> actInput = GetInputReader( readerId );
    
    // get all defined result types
    StdVector<shared_ptr<ResultInfo> > infos;
    GetResultTypes( readerId, sequenceStep, infos );
    
    // find correct one; if multiple are present -> Exception
    bool found = false;
    shared_ptr<ResultInfo> actInfo;
    for( UInt i = 0; i < infos.GetSize(); i++ ) {
      if( infos[i]->resultType == solType ) {
        // check, if result was already found
        if(found) { 
          EXCEPTION( "A result of type '" << SolutionTypeEnum.ToString(solType) << "' was already "
                      << "found in input reader '" << readerId
                      << "' in sequence Step " << sequenceStep );
        }
        actInfo = infos[i];
      }
    }
    
    // check if any result at all was found
    if( !actInfo ) {
     EXCEPTION( "Result was not found in input reader '"
                 << readerId << "'" );
    }
    
    // get all regions for given ResultInfo object
    StdVector<shared_ptr<EntityList> > entList;
    GetResultEntities( readerId, sequenceStep, actInfo, entList);
    
    // find correct one; if none is found -> Exception
    shared_ptr<EntityList> actList;
    for( UInt i = 0; i < entList.GetSize(); i++ ) {
      if( entList[i]->GetName() == regionName )
        actList = entList[i];
    }
    
    // check if any region at all was  found
    if( !actList) {
      EXCEPTION( "No entitylist found for result '"
                << solType << "' on region '" << regionName << "'" ); 
      }
      
    // determine analysistype of current multi sequence step
    std::map<UInt, BasePDE::AnalysisType> analysis;
    std::map<UInt, UInt> numSteps;
    GetNumMultiSequenceSteps( readerId, analysis, numSteps, false );
    
    // create new result object, fill it and return it
    shared_ptr<BaseResult> result;
    if( analysis[sequenceStep] != BasePDE::HARMONIC ) {
      result = shared_ptr<BaseResult>(new Result<Double>() );
    } else {
      result = shared_ptr<BaseResult>(new Result<Complex>() );
    }
    result->SetResultInfo( actInfo );
    result->SetEntityList( actList );
    GetResult( readerId, sequenceStep, stepValue, result);
    
    return result;
  }
  
  
  template<typename TYPE>
  shared_ptr<NodeStoreSol<TYPE> >
  ResultHandler::GetStoreSol( const std::string& readerId,
                              UInt sequenceStep,
                              UInt stepValue, 
                              SolutionType solType,
                              StdVector<std::string>& regionNames ) {
    
    
    // get grid
    Grid * ptGrid = static_pointer_cast<Result<TYPE> >
    (GetResult( readerId, sequenceStep,stepValue, solType, regionNames[0] ) )
    ->GetEntityList()->GetGrid();
    
    shared_ptr<ResultInfo> actInfo =  static_pointer_cast<Result<TYPE> >
    (GetResult( readerId, sequenceStep,stepValue, solType, regionNames[0] ) )
    ->GetResultInfo();             
    
    // create new equation map object
    EqnMap *  eqnMap = new EqnMap( ptGrid, 1, true);
    
    // iterate over all regionNames
    StdVector<shared_ptr<Result<TYPE> > > results;
    results.Resize( regionNames.GetSize() );
    for( UInt i = 0; i < regionNames.GetSize(); i++) {

      // obtain result object
      results[i] = static_pointer_cast<Result<TYPE> >(GetResult(readerId,
                   sequenceStep, stepValue, solType, regionNames[i] ) );

      // pass it ot eqnMap
      shared_ptr<EntityList> entList = results[i]->GetEntityList();
      eqnMap->AddResult( *actInfo, entList );
    }
    
    // finalize eqnMap
    eqnMap->Finalize();
    
    if(progOpts->DoListMapping()) 
      eqnMap->ToInfo(info->Get(ParamNode::HEADER)->Get("mappings", ParamNode::APPEND));
    
    // create new vector with coefficients
    Vector<TYPE> solVec;
    solVec.Resize( eqnMap->GetNumEqns() );
    
    // create new storesolution object
    shared_ptr<NodeStoreSol<TYPE> > sol( new NodeStoreSol<TYPE>());
    sol->SetNumSolutions( 1 );
    sol->SetNumNodes( 1 );
    sol->SetSolutionType( results[0]->GetResultInfo()->resultType );
    sol->SetNumDofs( results[0]->GetResultInfo()->dofNames.GetSize() );
    sol->SetResult( actInfo );
    sol->SetPtrEQNData( eqnMap, ptGrid );
    sol->Init();
    sol->SetAlgSysVector( solVec );

    
    // iterate over all regions
    Double max = 0.0;
    for( UInt i = 0; i < regionNames.GetSize(); i++) {

      // get result and entitylist
      shared_ptr<EntityList> regionList = results[i]->GetEntityList();

      // get related nodelist
      shared_ptr<EntityList> nodeList
      = ptGrid->GetEntityList(EntityList::NODE_LIST,
                              regionList->GetName(),
                              EntityList::REGION );
      EntityIterator it= nodeList->GetIterator();
      Vector<TYPE> & resVec = results[i]->GetVector();
      UInt pos = 0;
      StdVector<Integer> eqns;
      for( it.Begin(); !it.IsEnd(); it++ ) {

        // fetch equations
        eqnMap->GetEqns( eqns, *actInfo, it );
        //std::cerr << "equations: " << eqns.Serialize() << std::endl;

        // iterate over all equations
        for( UInt iEqn = 0; iEqn < eqns.GetSize(); iEqn++ ) {
          solVec[eqns[iEqn]-1] = resVec[pos++];
          if( std::abs(solVec[eqns[iEqn]-1]) > max ) {
            max = std::abs(solVec[eqns[iEqn]-1]);
          }
        } 
      }
    }
   
    
//    StdVector<Integer> eqns;
//    for( UInt i = 0; i < eqnMap->GetNumLocalNodes(); i++ ) {
//      StdVector<UInt> node(1);
//      node.Init(i+1);
//      NodeList list(ptGrid);
//      list.SetNodes( node);
//      EntityIterator nodeIt = list.GetIterator();
//      eqnMap->GetEqns( eqns, *actInfo, nodeIt);
//      std::cerr << "node " << i+1 << "\t";
//      std::cerr << solVec[eqns[0]-1] << ", " << solVec[eqns[1]-1] << std::endl;
//    }

//    std::cerr << "Storing algebraic pointer back\n";
    sol->SetAlgSysVector( solVec );
    //std::cerr << "solution vector is " << solVec << std::endl;
     
    return sol; 
    
  }


  //instantiate template methods
  template
  shared_ptr<NodeStoreSol<Double> >
  ResultHandler::GetStoreSol<Double>( const std::string& readerId,
                                      UInt sequenceStep,
                                      UInt stepValue, 
                                      SolutionType solType,
                                      StdVector<std::string>& regionNames );

  template
  shared_ptr<NodeStoreSol<Complex> >
  ResultHandler::GetStoreSol<Complex>( const std::string& readerId,
                                      UInt sequenceStep,
                                      UInt stepValue, 
                                      SolutionType solType,
                                      StdVector<std::string>& regionNames );
 

  
}
