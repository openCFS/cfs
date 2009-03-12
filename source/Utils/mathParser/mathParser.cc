// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "mathParser.hh"

#include <boost/algorithm/string/replace.hpp>

#include "MatVec/vector.hh"
#include "Utils/coordSystem.hh"
#include "Utils/interpolate.hh"
#include "DataInOut/ResultCache.hh"
#include "Domain/domain.hh"
#include "Utils/coordSystem.hh"


namespace CoupledField {
  
  
  // Static variable instantiation 
  std::list<Double> MathParser::dynamicPool_ = std::list<Double>();

  // Call method at muParser object with exception handling
#define MATHPARSER_EXEC(CALL)                                   \
  try { CALL; } catch(mu::Parser::exception_type &e) {          \
    EXCEPTION( "MathParser reports:\n -------------------\n"    \
               << " Message:  " << e.GetMsg() << "\n"           \
               << " Formula:  " << e.GetExpr() << "\n"          \
               << " Token:    " << e.GetToken() << "\n"         \
               << " Position: " << (int)e.GetPos() << "\n"      \
               << " Errc:     " << e.GetCode() << "\n";)        \
  }
  
  MathParser::MathParser() {

    // Create global math parser and initialize it
    InitParser( parsers_[GLOB_HANDLER], true );

    // Add global handle to activeHandles
    activeHandles_.insert(GLOB_HANDLER);

  }


  MathParser::~MathParser() {

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

  MathParser::HandleType MathParser::GetNewHandle() {
    
    // Obtain new handle
    HandleType newHandle = GLOB_HANDLER;
   
    if( activeHandles_.size() > 0 ) {
      std::set<HandleType>::iterator it = activeHandles_.end();
      newHandle = *(--it);
    }

    while( activeHandles_.find( newHandle) != activeHandles_.end() ) {
      newHandle++;
    }
    
    // Insert new handle
    activeHandles_.insert( newHandle );

    // Create new parser
    InitParser( parsers_[newHandle], false );

    // Create new entry in variable pool
    pools_[newHandle] = VarPool();
    
    // return handle
    return newHandle;

  }

  void MathParser::ReleaseHandle( HandleType handle) {

    // Check if handle is the global one
    if( handle == GLOB_HANDLER ) {
      EXCEPTION( "Can not dynamically release the GLOBAL handle!" );
    }

    // Check if handle exists
    if( activeHandles_.find(handle) == activeHandles_.end() ) {
      EXCEPTION( "MathParser handle '" << handle 
                 << "' is not registered!" );
    }

    // Release variable pool
    pools_.erase( handle );
    parsers_.erase( handle );

    // Remove handle from set
    activeHandles_.erase( handle );
  }

  void MathParser::SetValue( HandleType handler,
                             const std::string & varName,
                             Double val) {

    // Get parser related to handler
    mu::Parser & myParser  = GetParser( handler );
    
    // Look for related variable pool
    PoolMap::iterator poolsIt = pools_.find( handler );
    if ( poolsIt == pools_.end () ) {
      EXCEPTION( "RegisterVariable: Variable pool for handler '"
                 << handler << "' could not be found!" );
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


  void MathParser::SetCoordinates( HandleType handler,
                                   const CoordSystem & coosy,
                                   const Vector<Double> & globCoord ) {

    // Get local representation of global vector
    Vector<Double> locCoord;
    coosy.Global2LocalCoord( locCoord, globCoord );

    // Get component names and register them within the specified parser
    std::string tempName;
    for ( UInt i = 1; i <= locCoord.GetSize(); i++ ) {
      tempName = coosy.GetDofName(i);
      SetValue( handler, tempName, locCoord[i-1] );
    }

  }

  void MathParser::SetExpr( HandleType handler, const std::string & expr) {

    // Get parser related to handler
    mu::Parser & myParser = GetParser( handler );
    
    // Special handling in case the expression is set in a xml-file:
    // We have to replace all occurences of the single quote '
    // by the related double quote " in order to get muParser
    // recognize the expression as a string.
    std::string modExpr = expr;
    boost::algorithm::replace_all( modExpr, "'", "\"" );
    

    MATHPARSER_EXEC( myParser.SetExpr(modExpr) );
  }


  Double MathParser::Eval( HandleType handler ) {

    // Get parser related to handler
    mu::Parser & myParser = GetParser( handler );

    // Evaluate expression with error checking
    Double ret = 0.0;
    MATHPARSER_EXEC( ret = myParser.Eval() );

    return ret;
  }

  void MathParser::InitParser( mu::Parser &parser, bool isGlobal ) {
    
    // Register alias functions for >, >=, <=, <
    parser.DefineOprt( "gt", MathParser::Op_gt, 2);
    parser.DefineOprt( "ge", MathParser::Op_ge, 2);
    parser.DefineOprt( "le", MathParser::Op_le, 2);
    parser.DefineOprt( "lt", MathParser::Op_lt, 2);

    // Register constant variables
    parser.DefineConst("pi", (double) PI);

    // Register functions from within CFS
    parser.DefineFun("sample1D", Interpolate1D::Interpolate, false );
    parser.DefineFun("locCoord2D", MathParser::LocCoord2D, false );
    parser.DefineFun("locCoord3D", MathParser::LocCoord3D, false );
    parser.DefineFun("input", ResultCache::GetResult, false);


    // Register factory for dynamic variable registering
    //parser.SetVarFactory( AddVariable );

    // Check if parser is non-global one
    if ( isGlobal != true ) {

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
    
    dynamicPool_.push_back( 0.0 );
    return &dynamicPool_.back();
  }
  
  

  mu::Parser& MathParser::GetParser( HandleType handler ) {
    
    ParserMap::iterator it = parsers_.find( handler );

    if ( it == parsers_.end() ) {
      EXCEPTION( "GetParser: MathParser with handler '" 
                 << handler << "' not known!" );
    }
    return (*it).second;
  }

  Double MathParser::LocCoord3D( const char * coordSysId, 
                                 Double dof,
                                 Double x, Double y, Double z ) {
    CoordSystem * cosy = domain->GetCoordSystem(coordSysId);
    Vector<Double> loc(3), glob(3);
    loc[0] = x;
    loc[1] = y;
    loc[2] = z;
    cosy->Global2LocalCoord(loc, glob);
    
    return loc[(UInt)dof];
  }

  
  Double MathParser::LocCoord2D( const char * coordSysId, 
                                 Double dof,
                                 Double x, Double y) {

    CoordSystem * cosy = domain->GetCoordSystem(coordSysId);
    Vector<Double> loc(2), glob(2);
    loc[0] = x;
    loc[1] = y;
    cosy->Global2LocalCoord(loc, glob);
    
    return loc[(UInt)dof];
  }
  
  
}
