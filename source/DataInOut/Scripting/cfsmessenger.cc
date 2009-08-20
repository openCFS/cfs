// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "cfsmessenger.hh"

#include "PDE/SinglePDE.hh"
#include "DataInOut/Logging/cfslog.hh"
#include "Driver/assemble.hh"

namespace CoupledField {


  // define static variable for error message
  std::string CFSMessenger::errMsg_ = std::string();
  
  CFSMessenger::CFSMessenger() {

    // Initialize number of parameters for event procedures
    eventNumParams_[CFS_Init] = 1;
    eventNumParams_[CFS_PdeInit] = 1;
    eventNumParams_[CFS_ReadBCs] = 1;
    eventNumParams_[CFS_AssembleMat] = 4;
    eventNumParams_[CFS_AssembleRhs] = 3;
    eventNumParams_[CFS_SetBCs] = 3;
    eventNumParams_[CFS_CalcResults] = 3;
    eventNumParams_[CFS_Finish] = 1;

    isEvaluating_ = false;
    
  }
  
  CFSMessenger::~CFSMessenger() {

  }
  
  void CFSMessenger::ReadScriptFile( const std::string & fileName ) {
  }

  bool CFSMessenger::
  TriggerEvent( const EventType event, 
                const StdVector<std::string> & context) {

    // Nothing to do here, since we do not perform any scripting from
    // within c++ itself (...)
    return true;
  }
  
  void CFSMessenger::Error( const char * msg, const char * const filename,
                            const UInt numline) {
    // trick: Pretend, that no script is executing and simply
    // pass the error back to the global function
    isEvaluating_ = false;
    Info->Error( msg, filename, numline );
    
  }
  
  void CFSMessenger::Warning( const char * msg, const char * const filename,
                            const UInt numline) {
    // trick: Pretend, that no script is executing and simply
    // pass the error back to the global function
    isEvaluating_ = false;
    Info->Warning( msg, filename, numline );
    
  }
    
  bool CFSMessenger::CFSEval( const StdVector<std::string> & args,
                                 StdVector<std::string> & retVal ) {

    UInt argOffset = 2;
    bool success = false;
    
    if ( args.GetSize() == 1 ) {
      errMsg_ = "No arguments given";
      return false;
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

      // -- Assemble --
    } else if ( args[1] == "assemble" ) {
      success = domain->GetBasePDE()->getPDE_assemble()->Script_Eval(args, argOffset, retVal);

      // -- Grid --
    } else if ( args[1] == "grid" ) {
      success = domain->GetGrid()->Script_Eval(args, argOffset, retVal);

      // -- logConfigurator --
    } else if ( args[1] == "logConf" ) {
      success = logConf->Script_Eval(args, argOffset, retVal);
    } else {
      success = false;
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



