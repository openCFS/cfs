#include "mathParser.hh"

#include <def_build_type_options.hh>

#include <boost/algorithm/string/replace.hpp>
#include <boost/exception/diagnostic_information.hpp>

#include "registerfunc.hh"
#include "MatVec/Vector.hh"
#include "Domain/CoordinateSystems/CoordSystem.hh"
#include "Domain/Domain.hh"
#include "DataInOut/Logging/LogConfigurator.hh"

DEFINE_LOG(math, "mathParser")

namespace CoupledField {
  
// Some remarks:
// One new instance of muParser (i.e. when obtaining a new handle), requires 
// about 3-10 kB of memory. Additional memory is of course required, for
// registered variables.


  // Static variable instantiation 
  std::list<Double> MathParser::dynamicPool_ = std::list<Double>();

  // Call method at muParser object with exception handling
  // the difference between EXCEPTION and throw Exception is that the later does not report file and line
  // which is for this case for the user not of importance
#define MATHPARSER_EXEC(CALL)                                                                      \
  try { CALL; } catch(mu::Parser::exception_type &e) {                                             \
    throw Exception( "MathParser reports '" + e.GetMsg() + "' in formula '" + e.GetExpr() + "'");  \
  }
  // further codes .GetToken() and (int) e.GetPos() already in e.GetMsg().
  // also there is e.GetCode() (e.g. 1, 25, ...)
  
  MathParser::MathParser() {

    // Create new memory pool for global parser
    pools_[GLOB_HANDLER] = VarPool();
    
    // Create global math parser and initialize it
    InitParser( parsers_[GLOB_HANDLER], pools_[GLOB_HANDLER], 
                true, false );
    
    // Set default expression 0,0
    MathParser::SetExpr( GLOB_HANDLER, "0.0 ");
    
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
    std::map< unsigned int, PtSig >::iterator sigIt = exprChangeSignal_.begin(),
                                            itEnd = exprChangeSignal_.end();
    for ( ; sigIt != itEnd; ++sigIt ) {
      sigIt->second->disconnect_all_slots();
    }
    exprChangeSignal_.clear();
  }

  unsigned int MathParser::GetNewHandle( bool setDefaults ) {
    
    // Obtain new handle
    unsigned int newHandle = GLOB_HANDLER;
    if( activeHandles_.size() > 0 ) {
      std ::set<unsigned int>::iterator it = activeHandles_.end();
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
    MathParser::SetExpr( newHandle, "0.0 ");

    // return handle
    return newHandle;

  }

  void MathParser::ReleaseHandle( unsigned int handle) {

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

  void MathParser::SetValue( unsigned int handler,
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
      std::set<unsigned int> notifyParser;
      
      // always insert current handler
      notifyParser.insert( handler );
      
      // if we have the global handler, add depending parser
      if ( handler == GLOB_HANDLER ) {
        notifyParser.insert( globVarsInUse_[varName].begin(),
                             globVarsInUse_[varName].end() );
      }
      
      //Iterate over all parser instances and check, if
      // variable is in use
      std::set<unsigned int>::const_iterator handleIt = notifyParser.begin();
      for( ; handleIt != notifyParser.end(); ++handleIt ) {
        
        unsigned int actHandle = *handleIt;
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
          
          unsigned int actHandle = it->first;
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
  
  void MathParser::RegisterExternalVar(  unsigned int handle,
                                         const std::string& varName,
                                         Double * ptVar ) {
    LOG_DBG(math) << "registering '" << varName << "'\n";
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

        unsigned int actHandle = it->first;
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


  void MathParser::SetCoordinates( unsigned int handler,
                                   const CoordSystem & coosy,
                                   const Vector<Double> & globCoord ) {

    // Get local representation of global vector
    Vector<Double> locCoord;
    coosy.Global2LocalCoord( locCoord, globCoord );

    // Get component names and register them within the specified parser
    std::string tempName;

    UInt maxDim = locCoord.GetSize();
    if(maxDim==3 && coosy.GetDim()==2){
      //WARN("Dimension of coordinate vector is 3 but 2D setup ignoring third component");
      maxDim =2;
    }

    for ( UInt i = 1; i <= maxDim; i++ ) {
      tempName = coosy.GetDofName(i);
      MathParser::SetValue( handler, tempName, locCoord[i-1] );
    }

  }

  CoordPtrs MathParser::GetCoordPtrs( unsigned int handle,
                                      const CoordSystem& coosy ) {
    CoordPtrs ptrs;
    ptrs.numCoords = coosy.GetDim();

    // Get pointers to coordinate variables
    for ( UInt i = 1; i <= ptrs.numCoords; i++ ) {
      std::string varName = coosy.GetDofName(i);
      ptrs.coords[i-1] = GetValuePtr( handle, varName );
    }

    // Initialize remaining pointers to nullptr
    for ( UInt i = ptrs.numCoords; i < 3; i++ ) {
      ptrs.coords[i] = nullptr;
    }

    return ptrs;
  }

  void MathParser::SetCoordinatesDirect( const CoordPtrs& ptrs,
                                         const CoordSystem& coosy,
                                         const Vector<Double>& globCoord ) {
    // Get local representation of global vector
    Vector<Double> locCoord;
    coosy.Global2LocalCoord( locCoord, globCoord );

    UInt maxDim = locCoord.GetSize();
    if ( maxDim == 3 && coosy.GetDim() == 2 ) {
      maxDim = 2;
    }

    // Directly set coordinate values via cached pointers (no map lookups)
    for ( UInt i = 0; i < maxDim && i < ptrs.numCoords; i++ ) {
      if ( ptrs.coords[i] != nullptr ) {
        *(ptrs.coords[i]) = locCoord[i];
      }
    }
  }

  void MathParser::SetExpr( unsigned int handler, const std::string & expr) {

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
    
    LOG_DBG(math) << "SE expr=" << expr;
    LOG_DBG(math) << "SE modExpr=" << modExpr;

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


  Double MathParser::Eval( unsigned int handler )
  {
    // Get parser related to handler
    mu::Parser & myParser = GetParser( handler );

    #ifdef CHECK_INDEX
    // Consistency check: Ensure, that only 1 entry is set
    if( myParser.GetNumResults() > 1 ) {
      WARN("More than one expression set! Returning just last one.");
    }
    #endif
    
    // Evaluate expression with error checking
    try {
      return myParser.Eval();
    } catch(mu::Parser::exception_type &e) {
      std::string msg = "MathParser error '" + e.GetMsg() + "' in formula '" + e.GetExpr() + "' with registered variables " + GetRegisteredVariables(GLOB_HANDLER);
      if(handler != GLOB_HANDLER && GetRegisteredVariables(handler) != "")
        msg += ", " + GetRegisteredVariables(handler);
      throw Exception(msg);
    }
  }
  
  void MathParser::EvalVector( unsigned int handle, Vector<Double>& vec ) {
    
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

  Double MathParser::DiffVectorEntry(unsigned int handle, std::string varName, Integer VecPos){
    //basically a mod of the original diff implementation
    mu::Parser & myParser = GetParser( handle );
    VarPool &  curPool = pools_[handle];
    Double buffer = curPool[varName];
    Double eps = (buffer==0)? 1e-10 : 1e-7*buffer;

    Integer nExpr;
    mu::value_type *v = NULL;
    Double f1,f2,f3,f4;

    MathParser::SetValue( handle, varName, buffer + 2*eps );
    MATHPARSER_EXEC( v = myParser.Eval(nExpr));
    if(nExpr < VecPos)
      Exception("Invalid indices for vector diff");
    f1 = v[VecPos];

    MathParser::SetValue( handle, varName, buffer + 1*eps );
    MATHPARSER_EXEC( v = myParser.Eval(nExpr));
    if(nExpr < VecPos)
      Exception("Invalid indices for vector diff");
    f2 = v[VecPos];

    MathParser::SetValue( handle, varName, buffer - 1*eps );
    MATHPARSER_EXEC( v = myParser.Eval(nExpr));
    if(nExpr < VecPos)
      Exception("Invalid indices for vector diff");
    f3 = v[VecPos];

    MathParser::SetValue( handle, varName, buffer - 2*eps );
    MATHPARSER_EXEC( v = myParser.Eval(nExpr));
    if(nExpr < VecPos)
      Exception("Invalid indices for vector diff");
    f4 = v[VecPos];

    curPool[varName] =  buffer;
    MathParser::SetValue( handle, varName, buffer );
    return (-f1 + 8*f2 - 8*f3 + f4 ) / (12*eps);
  }

  void MathParser::EvalDivVector( unsigned int handle, Double& divergence ){

    //loop over variable pool and compute divergence
    divergence = 0.0;

    if(this->IsExprVariable(handle,"x")){
      divergence += MathParser::DiffVectorEntry(handle,"x",0);
    }
    if(this->IsExprVariable(handle,"y")){
      divergence += MathParser::DiffVectorEntry(handle,"y",1);
    }
    if(this->IsExprVariable(handle,"z")){
      divergence += MathParser::DiffVectorEntry(handle,"z",2);
    }
    if(this->IsExprVariable(handle,"r")){
      divergence += MathParser::DiffVectorEntry(handle,"r",0);
    }
    if(this->IsExprVariable(handle,"phi")){
      divergence += MathParser::DiffVectorEntry(handle,"phi",1);
    }
  }

     
  void MathParser::EvalMatrix( unsigned int handle, Matrix<Double>& matrix,
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
    
    // Register common math functions
    RegisterFunctions(parser);
    
    // Register functions for MathParser only (not for PyMuParser)
    // TODO: muParser needs to be extended to handle locCoord2D and locCoord2D
    // via strfun_type4 and strfun_type5. This has just been added to the current muParser master
    //parser.DefineFun("locCoord2D", MathParser::LocCoord2D, false );
    //parser.DefineFun("locCoord3D", MathParser::LocCoord3D, false );

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
  boost::signals2::connection MathParser::
  AddExpChangeCallBack( const MathParserSignal::slot_function_type
                        &subscriber,
                        unsigned int handle ) {
    
    if(  exprChangeSignal_.find( handle ) == exprChangeSignal_.end() ) {
      exprChangeSignal_[handle] = shared_ptr<MathParserSignal>(new MathParserSignal());
    }
    return (exprChangeSignal_[handle])->connect( subscriber );
  }

  
  std::string MathParser::GetExpr( unsigned int handle ) {
    // Get the map with the variables
    mu::Parser & actParser = GetParser( handle );
    std::string expr;
    MATHPARSER_EXEC(expr = actParser.GetExpr());
    return expr;
  }
  
  bool MathParser::IsExprConstant( unsigned int handle ) {
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
  
  bool MathParser::IsExprVariable( unsigned int handle, const std::string& var ) {
    StdVector<std::string> usedVars;
    //GetExprVars( handle, usedVars );
    bool found = false;
    if( varsInUse_[handle].find( var) != varsInUse_[handle].end() ) {
      found = true;
    }
    return found;
  }
  
  
  void MathParser::GetExprVars( unsigned int handle, 
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

  Double MathParser::GetExprVars( unsigned int handle, 
                                std::string varName ) {

    // Get the map with the variables
    mu::Parser& actParser = GetParser( handle );
    actParser.InitConst();

    // Get the constant variables
    mu::valmap_type valMap = actParser.GetConst();
    for (mu::valmap_type::iterator item = valMap.begin(); item!=valMap.end(); ++item)
    {
      if (item->first.compare(varName) == 0)
      {
        return item->second;
      }
    }
    mu::varmap_type varMap = actParser.GetVar();
    for (mu::varmap_type::iterator item = varMap.begin(); item!=varMap.end(); ++item)
    {
      if (item->first.compare(varName) == 0)
      {
        return *(item->second);
      }
    }
    EXCEPTION("Variable " << varName << " is not registered in mathparser");
  }


  Double* MathParser::GetValuePtr( unsigned int handle,
                                   const std::string& varName ) {
    // Look for related variable pool
    PoolMap::iterator poolsIt = pools_.find( handle );
    if ( poolsIt != pools_.end() ) {
      VarPool::iterator varIt = poolsIt->second.find( varName );
      if ( varIt != poolsIt->second.end() ) {
        return &(varIt->second);
      }
    }

    // If not in local pool and not global handler, check global pool
    if ( handle != GLOB_HANDLER ) {
      PoolMap::iterator globPoolsIt = pools_.find( GLOB_HANDLER );
      if ( globPoolsIt != pools_.end() ) {
        VarPool::iterator globVarIt = globPoolsIt->second.find( varName );
        if ( globVarIt != globPoolsIt->second.end() ) {
          return &(globVarIt->second);
        }
      }
    }

    return nullptr;
  }


  UInt MathParser::GetNumExprs( unsigned int handle ) {
    
    // Get parser related to handle
    mu::Parser & myParser = GetParser( handle );
    return UInt(myParser.GetNumResults() );
  }
  
  Double* MathParser::AddVariable( const char *varName ) {
    
    dynamicPool_.push_back( 0.0 );
    return &dynamicPool_.back();
  }
  
  
  StdVector<std::pair<std::string, double> > MathParser::GetRegisteredValues(unsigned int handle) const
  {
    StdVector<std::pair<std::string, double> > res;

    if(pools_.find(handle) != pools_.end())
    {
      const VarPool& pool = pools_.find(handle)->second;
      for(VarPool::const_iterator it = pool.begin(); it != pool.end(); ++it)
        res.Push_back(std::make_pair(it->first, it->second));
    }

    return res;
  }

  std::string MathParser::ToString(unsigned int handle) const
  {
    std::stringstream ss;

    for(const auto& pair : GetRegisteredValues(handle))
      ss << pair.first << ":" << pair.second << ", ";
    return ss.str();
  }

  std::string MathParser::GetRegisteredVariables(unsigned int handle) const
  {
    std::stringstream ss;

    StdVector<std::pair<std::string, double> > list = GetRegisteredValues(handle);
    for(unsigned int i = 0; i < list.GetSize(); i++)
    {
      ss << list[i].first;
      if(i < list.GetSize() -1)
        ss << ", ";
    }
    return ss.str();
  }


  void MathParser::ToInfo(PtrParamNode pn, unsigned int handle) const
  {
    StdVector<std::pair<std::string, double> > res = GetRegisteredValues(handle);
    for(unsigned int i = 0; i < res.GetSize(); i++)
      pn->Get(res[i].first)->SetValue(res[i].second);
  }


  void MathParser::Dump( std::ostream& out) {
    
    out << "====================\n"
           " MATH PARSER STATUS \n"
           "====================\n\n";   
    
    // iterate over all active Handles
    out << " 1) Active Handles\n"
           " -----------------\n";
    std::set<unsigned int>::const_iterator handleIt = activeHandles_.begin();
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
    
    std::map<std::string, std::set<unsigned int> >::const_iterator gIt;
    gIt = globVarsInUse_.begin();
    for( ; gIt != globVarsInUse_.end(); ++gIt) {
      out << "\tvar '"  << gIt->first << "' used in Parser handles: ";
      std::set<unsigned int>::const_iterator hIt =  gIt->second.begin();
      for( ; hIt != gIt->second.end(); ++hIt ) {
        out << *hIt << ",";
      }
      out << std::endl;
      
    }
    
    
  }
  

  mu::Parser& MathParser::GetParser( unsigned int handler ) {
    
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
    glob[0] = x;
    glob[1] = y;
    glob[2] = z;
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
    glob[0] = x;
    glob[1] = y;
    cosy->Global2LocalCoord(loc, glob);
    
    return loc[(UInt)dof-1];
  }
  
  
}
