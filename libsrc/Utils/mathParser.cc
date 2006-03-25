#include "mathParser.hh"
#include "Utils/vector.hh"
#include "Utils/coordSystem.hh"


namespace CoupledField {
  
  
  // Static variable instantiation 
  std::list<Double> MathParser::dynamicPool_ = std::list<Double>();

  // Call method at muParser object with exception handling
#define MATHPARSER_EXEC(CALL)                                   \
  try { CALL; } catch(mu::Parser::exception_type &e) {          \
    (*error) << "MathParser reports:\n -------------------\n"   \
             << " Message:  " << e.GetMsg() << "\n"             \
             << " Formula:  " << e.GetExpr() << "\n"            \
             << " Token:    " << e.GetToken() << "\n"           \
             << " Position: " << (int)e.GetPos() << "\n"        \
             << " Errc:     " << e.GetCode() << "\n";           \
    Error( __FILE__, __LINE__ );                                \
  }

  
  MathParser::MathParser() {
    ENTER_FCN( "MathParser::MathParser" );

    // Create global math parser and initialize it
    InitParser( parsers_[GLOB_HANDLER], TRUE );

  }


  MathParser::~MathParser() {
    ENTER_FCN( "MathParser::~MathParser" );

    // Clear all variables in Pools
    PoolMap::iterator mapIt;
    for ( mapIt = pools_.begin(); mapIt != pools_.end(); mapIt++ ) {
      mapIt->second.clear();
      
    }

    // Clear all parsers
    parsers_.clear();
    
    // Clear all dynamically allocated memory
    dynamicPool_.clear();
    
  }

  MathParser::HandlerType MathParser::GetNewHandler() {
    ENTER_FCN( "MathParser::GetNewHandler" );
    
    // Create new parser
    HandlerType newHandle = parsers_.size();
    InitParser( parsers_[newHandle], FALSE );

    // Create new entry in variable pool
    pools_[newHandle] = VarPool();
    
    // return handle
    return newHandle;

  }

  void MathParser::SetValue( HandlerType handler,
                             const std::string & varName,
                             Double val) {
    ENTER_FCN( "MathParser::RegisterVariable" );

    // Get parser related to handler
    mu::Parser & myParser  = GetParser( handler );
    
    // Look for related variable pool
    PoolMap::iterator poolsIt = pools_.find( handler );
    if ( poolsIt == pools_.end () ) {
      (*error) << "RegisterVariable: Variable pool for handler '"
               << handler << "' could not be found!";
      Error( __FILE__, __LINE__ );
    }

    // Get map entry for varName
    VarPool::iterator varIt =  (poolsIt->second).find(varName);

    if ( varIt != (poolsIt->second).end() ) {
      // Case1: variable was already defined

      // Simply change value in pool
      varIt->second = val;
    } else {
      // Case: variable was not defined yet

      // create new entry in variable pool
      poolsIt->second[varName] = val;

      /// register function with related parser object
        myParser.DefineVar( varName,&( poolsIt->second[varName])  );

        // If handler is the GLOB_HANDLER, register variables also
        // in all other parsers
        if ( handler == GLOB_HANDLER ) {

          // iterate over all local handler and define variable
          ParserMap::iterator it = parsers_.begin();
          it++;
          for (; it != parsers_.end(); it++ ) {
            it->second.DefineVar(  varName,&( poolsIt->second[varName]) );
          }
        }
    }
  }


  void MathParser::SetCoordinates( HandlerType handler,
                                   const CoordSystem & coosy,
                                   const Vector<Double> & globCoord ) {
    ENTER_FCN( "MathParser::SetCoordinates" );

    // Get local representation of global vector
    Vector<Double> locCoord;
    coosy.Global2LocalCoord( locCoord, globCoord );

    // Get component names and register them within the specified parser
    std::string tempName;
    for ( UInt i = 1; i <= locCoord.GetSize(); i++ ) {

      // Note: Since at the moment the names of the coordinates include mostly 
      // also  a prefix (e.g. "ux"), we simply skip that character.
      // As soon as the splitting of results and dimension is performed, this 
      // skipping can be removed.
      tempName = coosy.GetDofName(i);
      SetValue( handler, std::string( tempName, 1, tempName.size()), 
                locCoord[i-1] );
    }

  }


  Double MathParser::Eval( HandlerType handler,
                           const std::string & expr) {
    ENTER_FCN( "MathParser::Eval" );

    // Get parser related to handler
    mu::Parser & myParser = GetParser( handler );

    // Evaluate expression with error checking
    Double ret = 0.0;
    MATHPARSER_EXEC( myParser.SetExpr(expr) );
    MATHPARSER_EXEC( ret = myParser.Eval() );


    return ret;
  }

  void MathParser::InitParser( mu::Parser &parser, Boolean isGlobal ) {
    ENTER_FCN( "MathParser::InitParser" );
    
    // Register alias functions for >, >=, <=, <
    parser.DefineOprt( "gt", MathParser::Op_gt, 2);
    parser.DefineOprt( "ge", MathParser::Op_ge, 2);
    parser.DefineOprt( "le", MathParser::Op_le, 2);
    parser.DefineOprt( "lt", MathParser::Op_lt, 2);
    

    // Register factory for dynamic variable registering
    //parser.SetVarFactory( AddVariable );

    // Check if parser is non-global one
    if ( isGlobal != TRUE ) {

      // Register each global variable also within the local 
      // parser context
      VarPool & globPool = pools_[GLOB_HANDLER];
      VarPool::iterator varIt;
      for( varIt = globPool.begin(); varIt != globPool.end(); varIt++ ) {
        parser.DefineVar(  varIt->first, &(varIt->second) );
        
      } // for
    } // if
  }
  
  Double* MathParser::AddVariable( const char *varName ) {
    ENTER_FCN( "MathParser::AddVariable" );
    
    dynamicPool_.push_back( 0.0 );
    return &dynamicPool_.back();
  }
  
  

  mu::Parser& MathParser::GetParser( HandlerType handler ) {
    
    ParserMap::iterator it = parsers_.find( handler );

    if ( it == parsers_.end() ) {
      (*error) << "GetParser: MathParser with handler '" 
               << handler << "' not known!";
      Error( __FILE__, __LINE__ );
    }
    return (*it).second;
  }

}
