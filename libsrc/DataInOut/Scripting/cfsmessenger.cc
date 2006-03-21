#include "cfsmessenger.hh"

namespace CoupledField {


  // define static variable for error message
  std::string CFSMessenger::errMsg_ = std::string();
  
  CFSMessenger::CFSMessenger() {
    ENTER_FCN( "CFSMessenger::CFSMessenger" );

    // Initialize number of parameters for event procedures
    eventNumParams_[CFS_Init] = 1;
    eventNumParams_[CFS_ReadBCs] = 1;
    eventNumParams_[CFS_SetBCs] = 3;
    eventNumParams_[CFS_PostProcess] = 3;
  }
  
  CFSMessenger::~CFSMessenger() {

  }
  
  void CFSMessenger::ReadScriptFile( const std::string & fileName ) {
    ENTER_FCN( "CFSMessenger::ReadScriptFile" );
  }

  Boolean CFSMessenger::
  TriggerEvent( const EventType event, 
                const StdVector<std::string> & context) {
    ENTER_FCN( "CFSMessenger::TriggerEvent" );

    // Nothing to do here, since we do not perform any scripting from
    // within c++ itself (...)
    return TRUE;
  }
  
  void CFSMessenger::Error( const Char * msg, const Char * const filename,
                            const UInt numline) {
    // trick: Pretend, that no script is executing and simply
    // pass the error back to the global function
    isEvaluating_ = FALSE;
    Info->Error( msg, filename, numline );
    
  }
  
  void CFSMessenger::Warning( const Char * msg, const Char * const filename,
                            const UInt numline) {
    // trick: Pretend, that no script is executing and simply
    // pass the error back to the global function
    isEvaluating_ = FALSE;
    Info->Warning( msg, filename, numline );
    
  }
    
  Boolean CFSMessenger::CFSEval( const StdVector<std::string> & args,
                                 StdVector<std::string> & retVal ) {
    ENTER_FCN( "CFSMessenger::CFSEval" );

    UInt argOffset = 2;
    Boolean success = FALSE;
    
    if ( args.GetSize() == 1 ) {
      errMsg_ = "No arguments given";
      return FALSE;
    } 

    // Dispatch function call
    // -- PDEs --
    if ( args[1] == "electrostatic" ) {      
      success = domain->GetSinglePDE("electrostatic")->
        Script_Eval(args, argOffset, retVal);

    } else if ( args[1] == "mechanic" ) {
      success = domain->GetSinglePDE("mechanic")->
        Script_Eval(args, argOffset, retVal);

    } else if ( args[1] == "acoustic" ) {
      success = domain->GetSinglePDE("acoustic")->
        Script_Eval(args, argOffset, retVal);

    } else if ( args[1] == "magnetic" ) {
      success = domain->GetSinglePDE("magnetic")->
        Script_Eval(args, argOffset, retVal);

    } else if ( args[1] == "heatConduction" ) {
      success = domain->GetSinglePDE("heatConduction")->
        Script_Eval(args, argOffset, retVal);

      // -- Grid --
    } else if ( args[1] == "grid" ) {
      success = domain->GetGrid()->Script_Eval(args, argOffset, retVal);

    } else {
      success = FALSE;
      std::ostringstream msg;
      msg << "Command '" << args[1] << "' is not known!";
      errMsg_ = msg.str();
    }

    // Check, if operation was successful
    if ( errMsg_ == std::string() ) {
      Scriptable::Script_GetError( errMsg_ );
    }
    
    return success;
  }

} // end of namespace



