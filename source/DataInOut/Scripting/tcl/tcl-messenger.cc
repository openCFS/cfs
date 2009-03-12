// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "tcl-messenger.hh"
#include <string.h>
#include <sstream>




namespace CoupledField {
  
  // Declare static functions
   StdVector<std::string> TCL_CFSMessenger::curParams_;

  

  TCL_CFSMessenger::TCL_CFSMessenger(  ) {
    
    // create new interpreter object
    tcl_ = Tcl_CreateInterp();

    // register own functions
    RegisterFunctions(); 

    // register (dummy) event procedures
    RegisterEvents();

    

  }
  void TCL_CFSMessenger::ReadScriptFile( const std::string & fileName ) {
    // evaluate script file
    isEvaluating_ = true;
    int code = Tcl_EvalFile( tcl_, fileName.c_str() );
    if ( code != TCL_OK ) {
      std::string error = "Error in TCL file:\n ";
      if ( *tcl_->result != 0 )
        error += tcl_->result;
      Error( error.c_str(), __FILE__, __LINE__ );
    }
    isEvaluating_ = false;
  }
  
  
  TCL_CFSMessenger::~TCL_CFSMessenger() {

    // Delete interpreter object
    Tcl_DeleteInterp( tcl_ );
    tcl_ = NULL;
  }
  
  void TCL_CFSMessenger::Warning( const char * msg, const char * const filename,
                                  const UInt numline) {
    
    std::stringstream warn;
    warn <<  "TCL warning in function '"  << curEvent_ << "':\n";
    for (UInt i = 0; i < curParams_.GetSize(); i++ ) {
      warn << curParams_[i] << " ";
    }
    warn << ": "<< msg;

    // After having generated the correct error string,
    // the bucket is passed back to the global error handler
    isEvaluating_ = false;
    ::Warning( warn.str().c_str(), filename, numline );
    isEvaluating_ = true;
  }


  void TCL_CFSMessenger::Error( const char * msg, const char * const filename,
                                const UInt numline) {
    std::stringstream error;

    // Generate more accurate error message, if error occurs during
    // calling of a event procedure
    if (curEvent_ != std::string() ) {
      error <<  "TCL error in function '"  << curEvent_ << "':\n\n ";
      for (UInt i = 0; i < curParams_.GetSize(); i++ ) {
        error  << curParams_[i] << " ";
      }
      error << ": " << msg;
    } else {
      error << msg;
    }

    // After having generated the correct error string,
    // the bucket is passed back to the global error handler
    isEvaluating_ = false;
    EXCEPTION( error.str().c_str() );
                                              
  } 

  
  bool TCL_CFSMessenger::TriggerEvent( const EventType event, 
                                          const StdVector<std::string> & context) {
    
    // get name of event
    curEvent_ = eventNames_[event];

    std::stringstream procName;
    procName << eventNames_[event];
    for ( UInt i=0; i<context.GetSize(); i++ ) {
      procName << " " << context[i];
    }
    
    isEvaluating_ = true;
    int code = Tcl_Eval( tcl_, procName.str().c_str() );
    isEvaluating_ = false;

    if ( code != TCL_OK ) {
      std::string error = "TCL error in function '";
      error+= eventNames_[event];
      error+= "':\n\n ";
      if ( *tcl_->result != 0 )
        error += tcl_->result;
      EXCEPTION( error.c_str() );
      return false;
    }
    return true;
    
  }

  void TCL_CFSMessenger::RegisterFunctions() {

    // register Eval-function
    Tcl_CreateCommand( tcl_, "cfs", TCL_CFSMessenger::TCL_CFSEval,
                       (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
  }
  
  void TCL_CFSMessenger::RegisterEvents() {

    // First of all, generate mapping from event enumerations
    // to string representation
    eventNames_[CFS_Init] = "CFS_Init";
    eventNames_[CFS_ReadBCs] = "CFS_ReadBCs";
    eventNames_[CFS_SetBCs] = "CFS_SetBCs";
    eventNames_[CFS_CalcResults] = "CFS_CalcResults";
    eventNames_[CFS_Finish] = "CFS_Finish";
    
    // Check, if for all events the number
    // of parameters has been passed correctly
    if ( eventNumParams_.size() != eventNames_.size() ) {
      EXCEPTION( "Size of eventNumParams_ and eventNames_ "
               << "mismatch!\n Please ensure, that for each "
               << "event the number of parameters is defined!");
    }
    
    // Register all possible events
    std::map<EventType,std::string>::iterator it;
    for ( it = eventNames_.begin(); it != eventNames_.end(); it++ ) {
      std::stringstream procName;
      procName << "proc " << (*it).second;

      // Append dummy arguments for each parameter
      procName << " {";
      for ( UInt i = 0; i < eventNumParams_[(*it).first]; i++ ) {
        procName << " p" << i;
      }
      procName << " } {}";

      // register function ( = evaluate dummy proc body)
      Tcl_GlobalEval( tcl_, procName.str().c_str());
    }
    
  }

  int TCL_CFSMessenger::TCL_CFSEval(ClientData clientdata, Tcl_Interp *interp,
                                    int argc, const char *argv[]) {
    
    bool success;
    StdVector<std::string> retVal;

    // count number of arguments
    if (argc < 2 ) { 
      errMsg_ = "cfs needs at least 2 arguments!";
      success = false;
    } else {
      
      // convert arguments into std::strings
      StdVector<std::string> args; 
      args.Resize(argc); 
      for (UInt i=0; i<args.GetSize(); i++ ) {
        args[i] = argv[i];
      }
      
      // Save arguments to curent argument list
      curParams_ = args;
      
      // Call language independend CFSEval function
      success = CFSMessenger::CFSEval( args, retVal);
    }
    
    if ( success == true ) {
      
      // Append result to interpreter result
      Tcl_ResetResult( interp );
      for (UInt i=0; i<retVal.GetSize(); i++) {
        Tcl_AppendElement( interp, retVal[i].c_str());
      }
    } else {
      Tcl_SetResult(interp, strdup(errMsg_.c_str()), TCL_DYNAMIC);
      return TCL_ERROR;
    }
    return TCL_OK;
  }

    
} // end of namespace

