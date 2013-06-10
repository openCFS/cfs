#include "mathParser.hh"

#include <def_build_type_options.hh>

#include <boost/algorithm/string/replace.hpp>
#include "MatVec/Vector.hh"
#include "Domain/CoordinateSystems/CoordSystem.hh"
#include "Utils/Interpolate1D.hh"
#include "DataInOut/ResultCache.hh"
#include "Domain/Domain.hh"
#include "Utils/mathfunctions.hh"


namespace CoupledField {
  
// Some remarks:
// One new instance of muParser (i.e. when obtaining a new handle), requires 
// about 3-10 kB of memory. Additional memory is of course required, for
// registered variables.


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

    // Create new memory pool for global parser
    pools_[GLOB_HANDLER] = VarPool();
    
    // Create global math parser and initialize it
    InitParser( parsers_[GLOB_HANDLER], pools_[GLOB_HANDLER], 
                true, false );
    
    // Set default expression 0,0
    SetExpr( GLOB_HANDLER, "0.0 ");
    
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
    
    // Clear registered callback functions
    std::map< HandleType, PtSig >::iterator sigIt = exprChangeSignal_.begin(),
                                            itEnd = exprChangeSignal_.end();
    for ( ; sigIt != itEnd; ++sigIt ) {
      sigIt->second->disconnect_all_slots();
    }
    exprChangeSignal_.clear();
  }

  MathParser::HandleType MathParser::GetNewHandle( bool setDefaults ) {
    
    // Obtain new handle
    HandleType newHandle = GLOB_HANDLER;
    if( activeHandles_.size() > 0 ) {
      std ::set<HandleType>::iterator it = activeHandles_.end();
      newHandle = *(--it);
    }

    while( activeHandles_.find( newHandle) != activeHandles_.end() ) {
      newHandle++;
    }
    
    // Insert new handle
    activeHandles_.insert( newHandle );

    // Create new entry in variable pool
    pools_[newHandle] = VarPool();
    
    // Create new parser
    InitParser( parsers_[newHandle], pools_[newHandle], 
                false, setDefaults );
    
    // Initialize expression to 0.0
    SetExpr( newHandle, "0.0 ");

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
    varsInUse_.erase( handle );
    
    // Disconnect all connected callbacks
    if ( exprChangeSignal_.find(handle) != exprChangeSignal_.end() ) {
      exprChangeSignal_[handle]->disconnect_all_slots();
      exprChangeSignal_.erase( handle );
    }

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
      // -------------------------------------
      //  Case1: variable was already defined
      // -------------------------------------

      // Simply change value in pool
      varIt->second = val;
      
      // if we have global parser instance, notify
      // the depending instances as well for callback
      std::set<HandleType> notifyParser;
      
      // always insert current handler
      notifyParser.insert( handler );
      
      // if we have the global handler, add depending parser
      if ( handler == GLOB_HANDLER ) {
        notifyParser.insert( globVarsInUse_[varName].begin(),
                             globVarsInUse_[varName].end() );
      }
      
      //Iterate over all parser instances and check, if
      // variable is in use
      std::set<HandleType>::const_iterator handleIt = notifyParser.begin();
      for( ; handleIt != notifyParser.end(); ++handleIt ) {
        
        HandleType actHandle = *handleIt;
        // check, if variable is in use
        if( varsInUse_[actHandle].find( varName) 
            != varsInUse_[actHandle].end() ) {
        
          // check if signal is defined for this variable 
          // -> fire signal
          if(exprChangeSignal_.find(actHandle) != exprChangeSignal_.end() ) {
            // Fire signal
            (exprChangeSignal_[actHandle]).operator *()();
          }
        }
      }
      
    } else {
      
      // -------------------------------------
      //  Case2: variable was not defined yet
      // -------------------------------------

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
            
            HandleType actHandle = it->first;
            MATHPARSER_EXEC(
                it->second.DefineVar(  varName,&( poolsIt->second[varName]) );
            )

            // We consider the case, that we have in the 
            // child parsers already a "default" value for this variable,
            // i.e. the parser instance assumes that the variable is
            // locally defined. In this case we have to register
            // the variable also in the globVarsInUse_ map
            if( varsInUse_[actHandle].find( varName) !=
                varsInUse_[actHandle].end() ) {
              globVarsInUse_[varName].insert(actHandle);
              pools_[actHandle].erase( varName );
            }
          }
        }
    }
  }
  
  void MathParser::RegisterExternalVar(  HandleType handle,
                                         const std::string& varName,
                                         Double * ptVar ) {
    // Get parser related to handle
    mu::Parser & myParser  =  GetParser( handle );
    /// register function with related parser object
    myParser.DefineVar( varName, ptVar );

    // If handle is the GLOB_HANDLER, register variables also
    // in all other parsers
    if ( handle == GLOB_HANDLER ) {

      // iterate over all local handler and define variable
      ParserMap::iterator it = parsers_.begin();
      it++;
      for (; it != parsers_.end(); it++ ) {

        HandleType actHandle = it->first;
        MATHPARSER_EXEC(
            it->second.DefineVar( varName, ptVar );
        )

        // We consider the case, that we have in the 
        // child parsers already a "default" value for this variable,
        // i.e. the parser instance assumes that the variable is
        // locally defined. In this case we have to register
        // the variable also in the globVarsInUse_ map
        if( varsInUse_[actHandle].find( varName) !=
            varsInUse_[actHandle].end() ) {
          globVarsInUse_[varName].insert(actHandle);
          pools_[actHandle].erase( varName );
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
    
    // Clear all dependencies on variables of current expression
    std::set<std::string> & mySet = varsInUse_[handler];
    std::set<std::string>::iterator varIt;
    for( varIt = mySet.begin(); varIt != mySet.end(); ++varIt ) {
      
     // check if variable is also in use in global handler
      if( globVarsInUse_.find( *varIt ) != globVarsInUse_.end() ) {
        globVarsInUse_[*varIt].erase( handler );
      }
    }
    // in the end, clear the set itself
    mySet.clear();
    
    // Special handling in case the expression is set in a xml-file:
    // We have to replace all occurences of the single quote '
    // by the related double quote " in order to get muParser
    // recognize the expression as a string.
    std::string modExpr = expr;
    boost::algorithm::replace_all( modExpr, "'", "\"" );
    

    MATHPARSER_EXEC( myParser.SetExpr(modExpr) );
    
    // Now get all depending variables of this expression
    StdVector<std::string> varNames;
    GetExprVars( handler, varNames );
    StdVector<std::string>::iterator it = varNames.Begin();
    VarPool & actPool = pools_[handler];
    VarPool & globPool = pools_[GLOB_HANDLER];
    for( ; it != varNames.End(); it++ ) {
      varsInUse_[handler].insert(*it);

      // check, if variable is local or defined in GLOBAL_HANDLER
      if( actPool.find(*it) == actPool.end() ) {
        if( globPool.find(*it) != globPool.end() ) {
          globVarsInUse_[*it].insert( handler );
        }
      }
    }

  }


  Double MathParser::Eval( HandleType handler ) {

    // Get parser related to handler
    mu::Parser & myParser = GetParser( handler );

#ifdef CHECK_INDEX
    // Consistency check: Ensure, that only 1 entry is set
    if( myParser.GetNumResults() > 1 ) {
      WARN("More than one expression set! Returning just last one.");
    }
#endif
    
    // Evaluate expression with error checking
    Double ret = 0.0;
    MATHPARSER_EXEC( ret = myParser.Eval() );

    return ret;
  }
  
  void MathParser::EvalVector( HandleType handle, Vector<Double>& vec ) {
    
    // Get parser related to handler
    mu::Parser & myParser = GetParser( handle );

    // Evaluate list of expression and store it in vector
    Integer nExpr;
    mu::value_type *v = NULL;
    MATHPARSER_EXEC( v= myParser.Eval(nExpr) );
    vec.Resize(nExpr);
    for (Integer i = 0; i < nExpr; ++i )
    {
      vec[i] = v[i];
    }

  }
     
  void MathParser::EvalMatrix( HandleType handle, Matrix<Double>& matrix,
                               UInt numRows , UInt numCols ) {
    
    // Get parser related to handler
    mu::Parser & myParser = GetParser( handle );
    
    // Evaluate list of expressions
    Integer nExpr;
    mu::value_type *v = NULL;
    MATHPARSER_EXEC( v= myParser.Eval(nExpr) );

    // try to be smart about matrix size
    if( nExpr == 1) {
      // 1 entry
      matrix.Resize(1,1);
      matrix[0][0] = v[0];
    } else {
      // > 1 entries

      // If neither rows/col size is given, try to use size of matrix
      if( numRows == 0 && numCols == 0) {
        numRows = matrix.GetNumRows();
        numCols = matrix.GetNumCols();
      }

      // now we have to check, if two of them are zero
      if( numRows ==0 && numCols == 0) {
        EXCEPTION("Neither numRows/numCols was given, nor a "
            "pre- initialized matrix");
      }
      
      if( numRows == 0) {
        if( nExpr % numCols != 0) {
          EXCEPTION("Can not determine unique matrix size");
        } else {
          numRows = nExpr / numCols;
        }
      }
      
      if( numCols == 0) {
        if( nExpr % numRows != 0) {
          EXCEPTION("Can not determine unique matrix size.");
        } else {
          numCols = nExpr / numRows;
        }
      }
      
      // final consistency check
#ifdef CHECK_INDEX
      if( numCols * numRows != UInt(nExpr)) {
        EXCEPTION("Can not store " << nExpr << " entries in a "
                  << numRows << " x " << numCols << " matrix.");
      }
#endif
      // finally, copy values into matrix
      matrix.Resize(numRows,numCols);
      for( UInt i = 0; i < numRows; ++i ) {
        for( UInt j = 0; j < numCols; ++j ) {
          matrix[i][j] = v[i*numCols+j];
        }
      }

    }    
  }

  void MathParser::InitParser( mu::Parser &parser, VarPool& actPool,
                               bool isGlobal, bool setDefaults ) {
    
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

    // Register signal generating functions
    parser.DefineFun("sinBurst", SinBurst, false );
    parser.DefineFun("fadeIn", FadeIn, false );
    parser.DefineFun("spike", Spike, false );
    //parser.DefineFun("cosPulseComb", CosPulseComb, false );
    parser.DefineFun("squareBurst", SquarePulse, false );
    parser.DefineFun("gauss", Gauss, false );
    
    // Register general functions
    parser.DefineFun("mod", Mod, false );
    parser.DefineFun("besselCylJ", BesselCylJ, false );
    parser.DefineFun("besselCylY", BesselCylY, false );
    parser.DefineFun("besselSphJ", BesselSphJ, false );
    parser.DefineFun("besselSphY", BesselSphY, false );
    
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
    
      
      // if default variables should be set, we define missing variables
      if( setDefaults ) {
        StdVector<std::string> defaults;
        defaults = "t", "f", "x", "y", "z";
        StdVector<std::string>::iterator it = defaults.Begin();
        for( ; it != defaults.End(); ++it ) {
          if ( globPool.find( *it) == globPool.end() )  {
            // create new entry in variable pool
            actPool[*it] = 1.0;
            parser.DefineVar(  *it, &(actPool[*it] ) );
          }
        }
      }
    } // if
    
  }
  

  //! Register callback function for change of value of expression
  boost::signals::connection MathParser::
  AddExpChangeCallBack( const MathParserSignal::slot_function_type
                        &subscriber,
                        HandleType handle ) {
    
    if(  exprChangeSignal_.find( handle ) == exprChangeSignal_.end() ) {
      exprChangeSignal_[handle] = shared_ptr<MathParserSignal>(new MathParserSignal());
    }
    return (exprChangeSignal_[handle])->connect( subscriber );
  }

  
  std::string MathParser::GetExpr( HandleType handle ) {
    // Get the map with the variables
    mu::Parser & actParser = GetParser( handle );
    std::string expr;
    MATHPARSER_EXEC(expr = actParser.GetExpr());
    return expr;
  }
  
  bool MathParser::IsExprConstant( HandleType handle ) {
    // Get all depending variables of this expression
    //mu::Parser & actParser = GetParser( handle );
    
    mu::varmap_type variables;
    //MATHPARSER_EXEC( variables = actParser.GetUsedVar() ); 
    bool isConstant = true;
    if( varsInUse_[handle].size() != 0 ) {
      isConstant = false;
    }
    return isConstant;
  }
  
  bool MathParser::IsExprVariable( HandleType handle, const std::string& var ) {
    StdVector<std::string> usedVars;
    //GetExprVars( handle, usedVars );
    bool found = false;
    if( varsInUse_[handle].find( var) != varsInUse_[handle].end() ) {
      found = true;
    }
    return found;
  }
  
  
  void MathParser::GetExprVars( HandleType handle, 
                                StdVector<std::string>& varNames ) {

    // Get the map with the variables
    mu::Parser & actParser = GetParser( handle );
    mu::varmap_type variables;
    MATHPARSER_EXEC( variables = actParser.GetUsedVar() );

    // Get the number of variables 
    mu::varmap_type::const_iterator item = variables.begin();

    // Copy variable names
    varNames.Reserve( variables.size() );
    for (; item!=variables.end(); ++item) {
      varNames.Push_back( item->first );
    }
  }
  
  
  UInt MathParser::GetNumExprs( HandleType handle ) {
    
    // Get parser related to handle
    mu::Parser & myParser = GetParser( handle );
    return UInt(myParser.GetNumResults() );
  }
  
  Double* MathParser::AddVariable( const char *varName ) {
    
    dynamicPool_.push_back( 0.0 );
    return &dynamicPool_.back();
  }
  
  
  void MathParser::Dump( std::ostream& out) {
    
    out << "====================\n"
           " MATH PARSER STATUS \n"
           "====================\n\n";   
    
    // iterate over all active Handles
    out << " 1) Active Handles\n"
           " -----------------\n";
    std::set<HandleType>::const_iterator handleIt = activeHandles_.begin();
    for( ; handleIt != activeHandles_.end(); ++handleIt ) {
      out << "\t" << *handleIt << std::endl;
    }
   out << std::endl;
   
    // iterate over all active parser instances
   out << " 2) Parser instances\n"
          " --------------------\n";
    PoolMap::iterator poolIt = pools_.begin();
    ParserMap::iterator parseIt = parsers_.begin();
    
    
    for( ; parseIt != parsers_.end(); ++parseIt, ++poolIt ) {
      
      mu::Parser & actParser = parseIt->second;
      // - print handle
      out << " Handle: " << parseIt->first << std::endl
          << " ----------" << std::endl
          << "\tExpression: '" << actParser.GetExpr() << "'\n\n";
          
      // Print variables in use
      std::set<std::string> & mySet = varsInUse_[parseIt->first];
      std::set<std::string>::iterator setIt = mySet.begin();
      out << "\tUsed variables: ";
      for(; setIt != mySet.end(); ++setIt ) {
        out << *setIt << ", ";
      }
      out << std::endl;
      
      // Get the map with the variables
      mu::varmap_type variables = actParser.GetVar();
      out << "\n\t#Variables registered in parser: "
          << (int)variables.size() << "\n";

      // Get the number of variables 
      mu::varmap_type::const_iterator item = variables.begin();

      // Query the variables
      for (; item!=variables.end(); ++item) {
        out << "\t\tvar: '" << item->first 
            << "' Address: [0x" << item->second << "]\n";
      }
      
      // Print pool iunformation
      out << "\n\tVariables in pool:\n";
      VarPool & actPool  = poolIt->second;
      VarPool::const_iterator varIt = actPool.begin();
      for( ; varIt != actPool.end(); ++varIt ) {
        out << "\t\tvar: '" << varIt->first << "' \t value: " 
            << varIt->second << std::endl;
      }
      
      out << std::endl;
    }
    out << "\n\n";
    // print global variables information
    out << "3) Map: global variables <-> parser instances:\n"
        << "--------------------------------------------\n";
    
    std::map<std::string, std::set<HandleType> >::const_iterator gIt;
    gIt = globVarsInUse_.begin();
    for( ; gIt != globVarsInUse_.end(); ++gIt) {
      out << "\tvar '"  << gIt->first << "' used in Parser handles: ";
      std::set<HandleType>::const_iterator hIt =  gIt->second.begin();
      for( ; hIt != gIt->second.end(); ++hIt ) {
        out << *hIt << ",";
      }
      out << std::endl;
      
    }
    
    
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

    // check for the correct coordinate components
    if( dof < 1 || dof > 3 ) {
      EXCEPTION( "The coordinate component for a 3D system can just be 1,2 or 3");
    }
    CoordSystem * cosy = domain->GetCoordSystem(coordSysId);
    Vector<Double> loc(3), glob(3);
    loc[0] = x;
    loc[1] = y;
    loc[2] = z;
    cosy->Global2LocalCoord(loc, glob);
    return loc[(UInt)dof-1];
  }

  
  Double MathParser::LocCoord2D( const char * coordSysId, 
                                 Double dof,
                                 Double x, Double y) {
    // check for the correct coordinate components
    if( dof < 1 || dof > 2 ) {
      EXCEPTION( "The coordinate component for a 2D system can just be 1 or 2");
    }
    CoordSystem * cosy = domain->GetCoordSystem(coordSysId);
    Vector<Double> loc(2), glob(2);
    loc[0] = x;
    loc[1] = y;
    cosy->Global2LocalCoord(loc, glob);
    
    return loc[(UInt)dof-1];
  }
  
  
}
