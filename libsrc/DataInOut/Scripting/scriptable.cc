#include "scriptable.hh"

#include <cstdio>

namespace CoupledField
{

  // declaration of static class component
  std::stringstream Scriptable::errMsg_;
  
  Scriptable::Scriptable() {
    ENTER_FCN( "Scriptable::Scriptable" );

  }

  Scriptable::~Scriptable() {
    ENTER_FCN( "Scriptable::~Scriptable" );
  }

  
  void Scriptable::Script_GetCommands( StdVector<std::string> & commands,
                                     UInt & argOffset) {
    ENTER_FCN( "Scriptable::Script_GetCommands" );
      commands = "No commands available";
    }

  void Scriptable::Script_GetError( std::string & errMsg ) {
    ENTER_FCN( "Scriptable::Script_GetError" );

    errMsg = errMsg_.str();
    errMsg_.str("");
  }
  
} // end of namespace

