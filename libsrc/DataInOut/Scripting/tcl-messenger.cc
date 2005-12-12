#include "tcl-messenger.hh"
#include <string.h>
#include <sstream>



namespace CoupledField {
  




  TCL_CFSMessenger::TCL_CFSMessenger( const std::string & scriptFileName ) {
    ENTER_FCN( "CFSMessenger::TCL_CFSMessenger" );
    
    // create new interpreter object
    tcl_ = Tcl_CreateInterp();

    // register own functions
    RegisterFunctions(); 

    // register (dummy) event procedures
    RegisterEvents();

    // evaluate script file
    int code = Tcl_EvalFile( tcl_, scriptFileName.c_str() );
    if ( code != TCL_OK ) {
      std::string error = "CFS_Messenger::EvalFil: ";
      if ( *tcl_->result != 0 )
        error += tcl_->result;
      Error( error.c_str(), __FILE__, __LINE__ );
    }

  }
  
  TCL_CFSMessenger::~TCL_CFSMessenger() {
    ENTER_FCN( "TCL_CFSMessenger::~TCL_CFSMessenger()" );

    // Delete interpreter object
    Tcl_DeleteInterp( tcl_ );
    tcl_ = NULL;
  }
  
  Boolean TCL_CFSMessenger::TriggerEvent( const EventType event, 
                                          const StdVector<std::string> & context) {
    ENTER_FCN( "TCL_CFSMessenger::TriggerEvent" );
    
    // call related procedure in TCL interpreter
    std::stringstream procName;
    procName << eventNames_[event];
    for ( UInt i=0; i<context.GetSize(); i++ ) {
      procName << " " << context[i];
    }
    
    int code = Tcl_Eval( tcl_, procName.str().c_str() );
    
    if ( code != TCL_OK ) {
      std::string error = "CFS_Messenger::TriggerEvent: ";
      if ( *tcl_->result != 0 )
        error += tcl_->result;
      Error( error.c_str(), __FILE__, __LINE__ );
      return FALSE;
    }
    return TRUE;
    
  }

  void TCL_CFSMessenger::RegisterFunctions() {
    ENTER_FCN( "TCL_CFSMessenger::RegisterFunctions") ;

    // regier Set-function
    Tcl_CreateCommand( tcl_, "cfs_set", TCL_CFSMessenger::TCL_Set,
                       (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
    
    // register Get-function
    Tcl_CreateCommand( tcl_, "cfs_get", TCL_CFSMessenger::TCL_Get,
                       (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
  }
  
  void TCL_CFSMessenger::RegisterEvents() {
    ENTER_FCN( "TCL_CFSMessenger::RegisterFunctions") ;

    // First of all, generate mapping from event enumerations
    // to string representation
    eventNames_[CFS_Init] = "CFS_Init";
    eventNames_[CFS_ReadBCs] = "CFS_ReadBCs";
    eventNames_[CFS_SetBCs] = "CFS_SetBCs";
    eventNames_[CFS_PostProcess] = "CFS_PostProcess";
    
    // Check, if for all events the number
    // of parameters has been passed correctly
    if ( eventNumParams_.size() != eventNames_.size() ) {
      (*error) << "Size of eventNumParams_ and eventNames_ "
               << "mismatch!\n Please ensure, that for each "
               << "event the number of parameters is defined!";
      Error( __FILE__, __LINE__ );
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

  int TCL_CFSMessenger::TCL_Set(ClientData clientdata, Tcl_Interp *interp,
                                int argc, const char *argv[]) {
    ENTER_FCN( "TCL_CFSMessenger::TCL_Set") ;

    Boolean success;

    // count number of arguments
    if (argc < 2 ) { 
      errMsg_ = "set needs at least 2 arguments!";
      success = FALSE;
    } else {
      
      // convert arguments into std::strings
      StdVector<std::string> args;
      args.Resize(argc-1); 
      for (UInt i=0; i<args.GetSize(); i++ ) {
        args[i] = argv[i+1];
      }
      
      // Call global set function
      success = Set( args );
    }
    
    if ( success == TRUE ) {
      return TCL_OK;
    } else {
      Tcl_SetResult(interp, strdup(errMsg_.c_str()), TCL_DYNAMIC);
      return TCL_ERROR;
    }
  }   

  int TCL_CFSMessenger::TCL_Get(ClientData clientdata, Tcl_Interp *interp,
                                int argc, const char *argv[]) {
    ENTER_FCN( "TCL_CFSMessenger::TCL_Get");

    Boolean success;
    StdVector<std::string> retVal;

    // count number of arguments
    if (argc < 2 ) { 
      errMsg_ = "set needs at least 2 arguments!";
      success = FALSE;
    } else {
      
      // convert arguments into std::strings
      StdVector<std::string> args; 
      args.Resize(argc-1); 
      for (UInt i=0; i<args.GetSize(); i++ ) {
        args[i] = argv[i+1];
      }
      
      // Call global set function
      success = Get( args, retVal);
    }
    
    if ( success == TRUE ) {
      
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

