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
  
  Boolean CFSMessenger::
  TriggerEvent( const EventType event, 
                const StdVector<std::string> & context) {
    ENTER_FCN( "CFSMessenger::TriggerEvent" );

    // Nothing to do here, since we do not perform any scripting from
    // within c++ itself (...)
    return TRUE;
  }
  
  CFSMessenger::~CFSMessenger() {}
  Boolean CFSMessenger::Set( const StdVector<std::string> & args ) {
    ENTER_FCN( "CFSMessenger::Set" );
    
    StdVector<std::string> ret;
    UInt argOffset = 1;
    Boolean success = FALSE;

    if ( args.GetSize() == 0 ) {
      errMsg_ = "No arguments given";
      return FALSE;
    } 

    // Dispatch command invocation
    // -- PDEs --
    if ( args[0] == "electrostatic" ) {      
      success = domain->GetSinglePDE("electrostatic")->
        Script_Eval(args, argOffset, ret);

    } else if ( args[0] == "mechanic" ) {
      success = domain->GetSinglePDE("mechanic")->
        Script_Eval(args, argOffset, ret);

    } else if ( args[0] == "acoustic" ) {
      success = domain->GetSinglePDE("acoustic")->
        Script_Eval(args, argOffset, ret);

    } else if ( args[0] == "magnetic" ) {
      success = domain->GetSinglePDE("magnetic")->
        Script_Eval(args, argOffset, ret);

    } else if ( args[0] == "heatConduction" ) {
      success = domain->GetSinglePDE("heatConduction")->
        Script_Eval(args, argOffset, ret);

    // -- Grid --
    } else if ( args[0] == "grid" ) {
      success = domain->GetGrid()->Script_Eval(args, argOffset, ret);
      
    } else {
      success = FALSE;
    }

    // Check, if operation was successful
    if ( errMsg_ != std::string() ) {
      Scriptable::Script_GetError( errMsg_ );
    } else {
      std::ostringstream msg;
      msg << "Command '" << args[0] << "' is not known!";
      errMsg_ = msg.str();
    }

    return success;
  }
    
  Boolean CFSMessenger::Get(const StdVector<std::string> & args,
                            StdVector<std::string> & retVal ) {
    ENTER_FCN( "CFSMessenger::Set" );

    UInt argOffset = 1;
    Boolean success = FALSE;
    
    if ( args.GetSize() == 0 ) {
      errMsg_ = "No arguments given";
      return FALSE;
    } 

    // Dispatch function call
    // -- PDEs --
    if ( args[0] == "electrostatic" ) {      
      success = domain->GetSinglePDE("electrostatic")->
        Script_Eval(args, argOffset, retVal);

    } else if ( args[0] == "mechanic" ) {
      success = domain->GetSinglePDE("mechanic")->
        Script_Eval(args, argOffset, retVal);

    } else if ( args[0] == "acoustic" ) {
      success = domain->GetSinglePDE("acoustic")->
        Script_Eval(args, argOffset, retVal);

    } else if ( args[0] == "magnetic" ) {
      success = domain->GetSinglePDE("magnetic")->
        Script_Eval(args, argOffset, retVal);

    } else if ( args[0] == "heatConduction" ) {
      success = domain->GetSinglePDE("heatConduction")->
        Script_Eval(args, argOffset, retVal);

      // -- Grid --
    } else if ( args[0] == "grid" ) {
      success = domain->GetGrid()->Script_Eval(args, argOffset, retVal);

    } else {
      success = FALSE;
    }

    // Check, if operation was successful
    if ( errMsg_ != std::string() ) {
      Scriptable::Script_GetError( errMsg_ );
    } else {
      std::ostringstream msg;
      msg << "Command '" << args[0] << "' is not known!";
      errMsg_ = msg.str();
    }

    return success;
  }

} // end of namespace



