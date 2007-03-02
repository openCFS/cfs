#include "scriptable.hh"

#include <cstdio>

namespace CoupledField
{

  // =========================================================================
  //    Scriptable Memeber Functions
  // =========================================================================

  // Declaration of static class component
  std::stringstream Scriptable::errMsg_;
  
  Scriptable::Scriptable() {
    ENTER_FCN( "Scriptable::Scriptable" );

    isExecuting_ = false;
    currentArgs_ = NULL;

    // Register function for returnin all availabe functions
    // of this class
    FctPointer<Scriptable> * pt = 
      new FctPointer<Scriptable>(this, &Scriptable::Wrap_GetCommands);
    Script_RegisterFct( "getCommands", pt, ArgList() );
    
  }

  Scriptable::~Scriptable() {
    ENTER_FCN( "Scriptable::~Scriptable" );

    // Delete all fctPointers
    FctPtMap::iterator it;

    for ( it = fctPointers_.begin(); it != fctPointers_.end(); it++ ) {
      delete (*it).second;
    }

    // Reset pointer currentArgs_
    currentArgs_ = NULL;
  }

  void Scriptable::Script_RegisterFct( const std::string name,
                                       BaseFctPointer  * fctPointer,
                                       const ArgList & argList ) {
    // Check if function with this name was already registered
    if ( fctPointers_.find(name) != fctPointers_.end()
         && argLists_.find(name) != argLists_.end() ) {
      (*error) << "Function '" << name << "' was already registered!";
      Error( __FILE__, __LINE__ );
    }
    
    // Insert pointer and argList
    fctPointers_[name] = fctPointer;
    argLists_[name] = argList;
 
  }

  ArgList * Scriptable::Script_GetArgList(  ) {
    ENTER_FCN( "Scriptable::Script_GetArgList" );
    return currentArgs_;
  }
  
  bool Scriptable::Script_Eval( const StdVector<std::string> & args,
                                UInt & argOffset,
                                StdVector<std::string> & retVal) {
    ENTER_FCN( "Scriptable::Script_Eval" );
    
    // Save complete args vector in serial string
    std::stringstream argStr;
    for (UInt i = 0; i < args.GetSize(); i++ ) {
      argStr << args[i] << " ";
    }
    argStr << ": ";
    
    // Check if a function name at all was provided
    if (args.GetSize() - argOffset < 1) {
      errMsg_ << argStr.str()
              <<"Need at least one additional argument!";
      return false;
    }
    
    // Check if function name exists
    FctPtMap::iterator fctIt = fctPointers_.find( args[argOffset] );
    if( fctIt == fctPointers_.end() ) {
      errMsg_ << argStr.str()
              << "Function '" << args[argOffset] << "' does not exist!";
      return false;
    }

    // Get current argument list and save it to current list
    ArgListMap::iterator argListIt = argLists_.find( args[argOffset] );
    currentArgs_ = &((*argListIt).second);
  
    
    // Check, if correct number of arguments was provided
    std::string signature;
    UInt numArguments = args.GetSize()-argOffset - 1;
    if ( numArguments != currentArgs_->GetNumParams() ) {
      currentArgs_->GetSignature( signature );
      errMsg_ << argStr.str()
              << "Wrong number of arguments provided!";
      errMsg_ << "\n \n Usage: \033[31m" << args[argOffset] 
              <<  "\033[0m " << signature;
      return false;
    }
    
    // Copy parameters from message-list to argList
    StdVector<std::string> subArgs(numArguments);
    for ( UInt i = 0; i < numArguments; i++ ) {
      subArgs[i] = args[argOffset+1+i];
    }
    currentArgs_->SetParams( subArgs );
    
    // Set execution flag
    isExecuting_ = true;
    
    // Call Function object
    fctIt->second->Call();
    
    // Save Return value of argList object and delete it afterwars
    retVal = currentArgs_->GetRetVal();
    currentArgs_->GetRetVal().Clear();

    // Empty pointer to argList and reset status flag
    currentArgs_ = NULL;
    isExecuting_ = false;

    // No error occured, so leave without errors
    return true;
  }
  
  
  void Scriptable::Script_GetCommands( StdVector<std::string> & commands ) {
    ENTER_FCN( "Scriptable::Script_GetCommands" );
    ArgListMap::iterator it;
    std::stringstream out;
    std::string signature;
    
    for ( it = argLists_.begin(); it != argLists_.end(); it++ ) {
      (*it).second.GetSignature( signature );
      out << "\033[31m" << (*it).first << "\033[0m "  
          << signature;
      commands.Push_back( out.str() );
      out.str("");
    }
  }

  void Scriptable::Script_GetError( std::string & errMsg ) {
    ENTER_FCN( "Scriptable::Script_GetError" );

    errMsg = errMsg_.str();
    errMsg_.str("");
  }
  
  void Scriptable::Wrap_GetCommands() {
    Script_GetCommands( SCRIPT_RETVAL );
  }

  // =========================================================================
  //    ArgList methods
  // =========================================================================
  
  UInt ArgList::GetNumParams() {
    return orderedParams_.GetSize();
  }
  
  void ArgList::RegisterParam( const std::string name, 
                               ParamType paramType ) {
    ENTER_FCN( "ArgList::RegisterParam" );
    
    // Check if parameter with this name already exists
    if ( orderedParams_.Find(name) != -1 ) {
      *error << "Parameter with name '" << name << "' already exists!";
      Error( __FILE__, __LINE__ );
    }

    // Add paramter
    orderedParams_.Push_back(name);
    paramMap_[name] = paramType;
    
  }

  void ArgList::GetSignature( std::string & signature ) {
    ENTER_FCN( "ArgList::GetSignature" );
    std::stringstream out;
    
    // Iterate over all arguments
    for (UInt i = 0; i < orderedParams_.GetSize(); i++ ) {
      out << "\033[31m" << orderedParams_[i] << "\033[0m (\033[32m" 
          << ParamType2String(paramMap_[orderedParams_[i]])
          << "\033[0m) ";
    }
    
    signature = out.str();
  }

  void ArgList::SetParams( const StdVector<std::string> params ) {
    ENTER_FCN( "ArgList::SetParams" );

    StdVector<std::string> help;
    
    for ( UInt i = 0; i < params.GetSize(); i++ ) {

      // Clear help vector
      help.Clear();

      // Get name of paramter
      std::string name = orderedParams_[i];
      
      // Get type of parameter
      ParamType type = paramMap_[name];

      
      
      // Dependend on the type, perform a type conversion
      switch (type) {
      case STRING:
        STRINGpool_[name] = params[i];
        break;
      case DOUBLE:
        DOUBLEpool_[name] = String2Double( params[i] );
        break;
      case UINT:
        UINTpool_[name] = String2UInt( params[i] );
        break;
      case INT:
        INTpool_[name] = String2Int( params[i] );
        break;
      case STDVEC_STR:
        SplitStringList( params[i], STDVEC_STRpool_[name], ' ');
        break;
      case STDVEC_DOUBLE:
        SplitStringList( params[i], help, ' ' );
        String2Double( STDVEC_DOUBLEpool_[name], help );
        break;
      case STDVEC_UINT:
        SplitStringList( params[i], help, ' ' );
        String2UInt( STDVEC_UINTpool_[name], help );
        break;
      case STDVEC_INT:
        SplitStringList( params[i], help, ' ' );
        String2Int( STDVEC_INTpool_[name], help );
        break;
      case VEC_UINT:
        SplitStringList( params[i], help, ' ' );
        String2UInt( VEC_UINTpool_[name], help );
        break;
      case VEC_INT:
        SplitStringList( params[i], help, ' ' );
        String2Int( VEC_INTpool_[name], help );
        break;
      case VEC_DOUBLE:
        SplitStringList( params[i], help, ' ' );
        String2Double( VEC_DOUBLEpool_[name], help );
        break;
      } // end case
    } // end for
  }
  
  std::string ArgList::ParamType2String( const ParamType type ) {
    ENTER_FCN( "ArgList::ParamType2String" );
    
    std::string out;
    
    switch(type) {
    case STRING:
      out = "string";
      break;
    case DOUBLE:
      out = "double";
      break;
    case UINT:
      out = "unsigned int";
      break;
    case INT:
      out = "integer";
      break;
    case STDVEC_STR:
      out = "vector<string>";
      break;
    case STDVEC_DOUBLE:
      out = "vector<double>";
      break;
    case STDVEC_UINT:
      out = "vector<unsigned int>";
      break;
    case STDVEC_INT:
      out = "vector<int>";
      break;
    case VEC_UINT:
      out = "vector<unsigned int>";
      break;
    case VEC_INT:
      out = "vector<int>";
      break;
    case VEC_DOUBLE:
      out = "vectr<double>";
      break;
    default:
      Error( "No conversion found for given ParamType",
             __FILE__, __LINE__ );
    }
    return out;
  }
  
  StdVector<std::string> & ArgList::GetRetVal() {
    ENTER_FCN( " ArgList::GetRetVals" );
    
    return returnVals_;
  }
  
} // end of namespace

