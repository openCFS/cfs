#include "postProc.hh"

#include "Utils/result.hh"
#include "Domain/domain.hh"
#include "Domain/resultInfo.hh"
#include "Domain/grid.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"

namespace CoupledField {



  PostProc::PostProc( Grid * ptGrid, ParamNode * postProcNode ) {
    ENTER_FCN( "PostProc::PostProc" );

    ptGrid_ = ptGrid;
    myParam_ = postProcNode;
    reducType_ = NONE;
    writeResult_ = true;

  }

  PostProc::~PostProc() {
    ENTER_FCN( "PostProc::~PostProc" );
  }

  void PostProc::GetOutDestNames( StdVector<std::string> & outNames ) {
    ENTER_IFCN( "PostProc::GetOutDestNames" );
    outNames = outputNames_;
  }

  bool PostProc::IsWriteResult( ) {
    ENTER_IFCN( "PostProc::IsWriteResult" );
    return writeResult_;
  }

  void PostProc:: GetReducedList( shared_ptr<EntityList>& outList,
                                  ResultInfo::EntityUnknownType& outType,
                                  const shared_ptr<EntityList> inList,
                                  const ResultInfo::EntityUnknownType inType,
                                  const ReductionType reducType ) {
    ENTER_FCN( "PostProc::GetReducedList" );

    shared_ptr<EntityList> retList;
    
    if( reducType == NONE ||
        reducType == TIME_FREQ ) {
      outList = inList;
      outType = inType;
    } 

    if( reducType == SPACE ) {

      // Check original type
      // a) NodeList (by Region) -> REGION-List
      // b) ElemList (by Region) -> REGION-List
      // c) NodeList (by Name)   -> NODE_NAME-List
      // d) ElemList (by Name)   -> ELEM_NAME-List
      // e) RegionList           -> Not possible!
      // f) CoilList             -> Not possible!
      // g) FreeList             -> Not possible!
      if( inList->GetType() == EntityList::NODE_LIST ) {
        
        if( inList->GetDefineType() == EntityList::REGION )  {
          // case a)
          outType = ResultInfo::REGION;
          outList = ptGrid_->GetEntityList( EntityList::REGION_LIST,
                                            inList->GetName(),
                                            EntityList::REGION );
        } else {
          // case c)
          EXCEPTION( "Not yet implemented" );
        } 
      } else if( inList->GetType() == EntityList::ELEM_LIST ||
                 inList->GetType() == EntityList::SURF_ELEM_LIST ) {
        if( inList->GetDefineType() == EntityList::REGION ) {
          // case b)
          outType = ResultInfo::REGION;
          outList = ptGrid_->GetEntityList( EntityList::REGION_LIST,
                                            inList->GetName(),
                                            EntityList::REGION );
        } else {
          // case d)
          EXCEPTION( "Not yet implemented" );
        } 
      } else {
        // case e)
        // case f)
        // case g)
        EXCEPTION( "Not yet implemented" );
      }
    }
  }

  void PostProc::CreatePostProc( ParamNode * procNode,
                                 Grid * ptGrid, 
                                 StdVector<shared_ptr<PostProc> > & postProcs ) {
    ENTER_FCN( "PostProc::CreatePostProc" );

    // fetch single postProcs defined in this section
    StdVector<ParamNode*> procNodes = procNode->GetChildren();

    // iterate over all single postprocessing definitions
    for( UInt i = 0; i < procNodes.GetSize(); i++ ) {
      
      shared_ptr<PostProc> actProc;
      bool found = false;
      std::string procName = procNodes[i]->GetName();
      
      // NOTE: As also the attributes of the "postProc"-node
      // are considered as children of the postProcNode,
      // we have to check if the child name is "name"
      if( procName == "id" )
        continue;


      // === SUM ===
      if( procName == "sum" ) {
        found = true;
        std::string typeString;
        procNodes[i]->Get( "reduction", typeString );
        ReductionType type;
        if( typeString == "space" ) {
          type = SPACE;
        } else {
          type = TIME_FREQ;
        }
        actProc = 
          shared_ptr<PostProcSum>(new PostProcSum( ptGrid, type, procNodes[i] ) );
      } 
      
      // === MAXIMUM ===
      if( procName == "max" ) {
        found = true;
        std::string typeString;
        ReductionType type;
        procNodes[i]->Get( "reduction", typeString );
        if( typeString == "space" ) {
          type = SPACE;
        } else {
          type = TIME_FREQ;
        }
        actProc = shared_ptr<PostProcMax>(new PostProcMax( ptGrid, type,NULL ) );
        } 
      
      
      // === FUNCTION ===
      if( procName == "function" ) {      
        found = true;
        shared_ptr<PostProcFunc> funcProc;
        funcProc = shared_ptr<PostProcFunc>(new PostProcFunc( ptGrid, NULL ) );

        // get resultName
        std::string resultName;
        procNodes[i]->Get( "resultName", resultName );

        // get unit
        std::string unit;
        procNodes[i]->Get( "unit", unit );

        // get dof-nodes
        StdVector<ParamNode*> dofNodes = procNodes[i]->GetList("dof");;
        
        // get dofnames, real functions and imaginary functions for each dof
        StdVector<std::string> dofNames, rFunctions, iFunctions;
        for( UInt iDof = 0; iDof < dofNodes.GetSize(); iDof++ ) {
          dofNames.Push_back( dofNodes[iDof]->Get("name")->AsString() );
          rFunctions.Push_back( dofNodes[iDof]->Get("realFunc")->AsString() );
          iFunctions.Push_back( dofNodes[iDof]->Get("imagFunc")->AsString() );
        }
        funcProc->Initialize( resultName, unit, dofNames, rFunctions, iFunctions );
        actProc = funcProc;
      }

      // Check, if postProc method was found
      if (!found) {
        EXCEPTION( "Postprocessing method with name '" << procName
                   << "' is not known." );
      }
      
      // get next postprocessing routine
      procNodes[i]->Get( "postProcId", actProc->next_ );

      // get writeResult
      procNodes[i]->Get( "writeResult", actProc->writeResult_, true );

      // get outputDestNames
      std::string outNames;
      procNodes[i]->Get( "outputIds", outNames );
      SplitStringList( outNames, actProc->outputNames_, ',' );
      postProcs.Push_back( actProc );
    }

  }

  // =================== SUM-POSTPROCEDURE ===================
  
  PostProcSum::PostProcSum( Grid * ptGrid, ReductionType type, 
                            ParamNode * postProcNode )
    : PostProc( ptGrid, postProcNode ) {
    ENTER_FCN( "PostProcSum::PostProcSum" );
    
    name_ = "sum";
    reducType_ = type;
  }

  PostProcSum::~PostProcSum() {
    ENTER_FCN( "PostProcSum::~PostProcSum" );
    
  }


 

  
  void PostProcSum::SetResult( shared_ptr<BaseResult> res ) {
    ENTER_FCN( "PostProcSum::SetResult" );
    
    // Append result to list of results
    input_ =  res;
    
    // Create correct output result
    // a) Fetch related (reduced) entity list and EntityUnknownType
    shared_ptr<EntityList> newList;
    ResultInfo::EntityUnknownType newDefinedOn;
    ResultInfo & inResult  = *(res->GetResultInfo());
    
    GetReducedList( newList, newDefinedOn, res->GetEntityList(), 
                    inResult.definedOn, reducType_ );
    
    // b) Adjust new ResultInfo
    shared_ptr<ResultInfo> newInfo = 
      shared_ptr<ResultInfo>( new ResultInfo( inResult ) );
    newInfo->definedOn = newDefinedOn;
    newInfo->resultName += "-sum";
    
    // c) Create new Result and set is as output
    shared_ptr<BaseResult> newResult;
    if( res->GetEntryType() == EntryType::DOUBLE ) {
      newResult = shared_ptr<BaseResult>(new Result<Double>() );
    } else {
      newResult = shared_ptr<BaseResult>( new Result<Complex>() );
    }
    newResult->SetEntityList( newList );
    newResult->SetResultInfo( newInfo );
    output_ = newResult;

  }

  void PostProcSum::Apply() {
    ENTER_FCN( "PostProcSum::Apply" );
    
    if( output_->GetEntryType() == EntryType::DOUBLE ) {
      CalcSum<Double>();
    } else {
      CalcSum<Complex>();
    }
  }


  template<class TYPE> 
  void PostProcSum::CalcSum() {
    ENTER_FCN( "PostProcSum::CalcSum" );
    
    // Cast output into correct type
    Result<TYPE> & out = dynamic_cast<Result<TYPE>&>(*output_);
    Result<TYPE> & in =  dynamic_cast<Result<TYPE>&>(*input_ );
    Vector<TYPE> & outVec = out.GetVector();
    Vector<TYPE> & inVec = in.GetVector();
      
  
    if( reducType_ == SPACE ) {

      // If ReductionType == SPACE, we have to iterate over all
      // entities of output_ and get the related EntityList in
      // input_..
      
      shared_ptr<EntityList> outList = out.GetEntityList();
      EntityIterator outIt = outList->GetIterator();
      ResultInfo & outResInfo = *(out.GetResultInfo() );
      outVec.Resize( outResInfo.dofNames.GetSize() );
      outVec.Init();
      outIt.Begin();
      
      EntityIterator inIt = in.GetEntityList()->GetIterator();
      UInt numDofs = in.GetResultInfo()->dofNames.GetSize();
      
      // Iterate over all entities of the input entitylist
      for( inIt.Begin(); !inIt.IsEnd(); inIt++ ) {
        // Iterate over all dofs
        for( UInt iDof = 0; iDof < numDofs; iDof++ ) {
          outVec[iDof] += inVec[inIt.GetPos()*numDofs + iDof];
        } // dofs
      } // input entities
    } else {
      
      // Iterate over all vector entries
      if( outVec.GetSize() != inVec.GetSize() ) {
        outVec.Resize( inVec.GetSize() );
        outVec.Init();
      }
      outVec += inVec;

    }

  }
  
  void PostProcSum::Finalize() {
    ENTER_FCN( "PostProcSum::Finalize" );

  }

// =================== MAXIMUM-POSTPROCEDURE ===================

  PostProcMax::PostProcMax( Grid * ptGrid, ReductionType type,
                            ParamNode * postProcNode )
    : PostProc( ptGrid, postProcNode ) {
    ENTER_FCN( "PostProcMax::PostProcMax" );
    
    name_ = "max";
    reducType_ = type;
  }

  PostProcMax::~PostProcMax() {
    ENTER_FCN( "PostProcMax::~PostProcMax" );
    
  }


 

  
  void PostProcMax::SetResult( shared_ptr<BaseResult> res ) {
    ENTER_FCN( "PostProcMax::SetResult" );
    
    // Append result to list of results
    input_ =  res;
    
    // Create correct output result
    // a) Fetch related (reduced) entity list and EntityUnknownType
    shared_ptr<EntityList> newList;
    ResultInfo::EntityUnknownType newDefinedOn;
    ResultInfo & inResult  = *(res->GetResultInfo());
    
    GetReducedList( newList, newDefinedOn, res->GetEntityList(), 
                    inResult.definedOn, reducType_ );
    
    // b) Adjust new ResultInfo
    shared_ptr<ResultInfo> newInfo = 
      shared_ptr<ResultInfo>( new ResultInfo( inResult ) );
    newInfo->definedOn = newDefinedOn;
    newInfo->resultName += "-max";
    
    // c) Create new Result and set is as output
    shared_ptr<BaseResult> newResult;
    if( res->GetEntryType() == EntryType::DOUBLE ) {
      newResult = shared_ptr<BaseResult>(new Result<Double>() );
    } else {
      newResult = shared_ptr<BaseResult>( new Result<Complex>() );
    }
    newResult->SetEntityList( newList );
    newResult->SetResultInfo( newInfo );
    output_ = newResult;

  }

  void PostProcMax::Apply() {
    ENTER_FCN( "PostProcMax::Apply" );
    
    if( output_->GetEntryType() == EntryType::DOUBLE ) {
      CalcMax<Double>();
    } else {
      CalcMax<Complex>();
    }
  }


  template<class TYPE> 
  void PostProcMax::CalcMax() {
    ENTER_FCN( "PostProcMax::CalcMax" );
    
    // Cast output into correct type
    Result<TYPE> & out = dynamic_cast<Result<TYPE>&>(*output_);
    Result<TYPE> & in =  dynamic_cast<Result<TYPE>&>(*input_ );
    Vector<TYPE> & outVec = out.GetVector();
    Vector<TYPE> & inVec = in.GetVector();
      
  
    if( reducType_ == SPACE ) {

      // If ReductionType == SPACE, we have to iterate over all
      // entities of output_ and get the related EntityList in
      // input_..
      
      shared_ptr<EntityList> outList = out.GetEntityList();
      ResultInfo & outResInfo = *(out.GetResultInfo() );
      outVec.Resize( outResInfo.dofNames.GetSize() );
      outVec.Init();


      EntityIterator outIt = outList->GetIterator();
      EntityIterator inIt = in.GetEntityList()->GetIterator();
      outIt.Begin();
      inIt.Begin();
      UInt numDofs = in.GetResultInfo()->dofNames.GetSize();
      
      // Iterate over all entities of the input entitylist
      for( inIt.Begin(); !inIt.IsEnd(); inIt++ ) {
        // Iterate over all dofs
        for( UInt iDof = 0; iDof < numDofs; iDof++ ) {
          outVec[iDof] = std::max( inVec[inIt.GetPos()*numDofs + iDof], 
                                   outVec[iDof], &opLtSingleDof<TYPE> );
        } // dofs
      } // input entities
    } else {
      
      // Iterate over all vector entries
      if( outVec.GetSize() != inVec.GetSize() ) {
        outVec.Resize( inVec.GetSize() );
        outVec.Init();
      }
      
      for( UInt i = 0 ; i< outVec.GetSize(); i++ ) {
        outVec[i] =  std::max(inVec[i], outVec[i], &opLtSingleDof<TYPE> );
      }
    }


  }
  
  void PostProcMax::Finalize() {
    ENTER_FCN( "PostProcMax::Finalize" );

  }

  template<class TYPE>
  bool PostProcMax::opLtSingleDof( TYPE a, TYPE b ) {
    return a < b;
  }

  template<>
  bool PostProcMax::opLtSingleDof( Complex a, Complex b ) {
    return std::abs(a) < std::abs(b);
  }

  

 // =================== FUNC-POSTPROCEDURE ===================
  

  PostProcFunc::PostProcFunc( Grid * ptGrid,
                              ParamNode * postProcNode )
    : PostProc( ptGrid, postProcNode ) {
    ENTER_FCN( "PostProcFunc::PostProcFunc" );

    name_ = "func";
    reducType_ = NONE;

    // fetch math parser object
    mParser_ = domain->GetMathParser();
  }

  PostProcFunc::~PostProcFunc() {
    ENTER_FCN( "PostProcFunc::~PostProcFunc" );
    
  }

  void PostProcFunc::SetResult( shared_ptr<BaseResult> res ) {
    ENTER_FCN( "PostProcFunc::SetResult" );

    // Set result as input
    input_ =  res;

    // intialize rVals and mathParser handles
    // a) double
    rVals_.Resize( dofNames_.GetSize() );
    rVals_.Init();
    rHandles_.Resize( dofNames_.GetSize() );
    // b) complex
    if( input_->GetEntryType() == EntryType::COMPLEX ) {
      iVals_.Resize( dofNames_.GetSize() );
      iVals_.Init();
      iHandles_.Resize( dofNames_.GetSize() );
    }
    
    // register all new functions with mathParser and new result info  object
    StdVector<std::string> & inDofNames = input_->GetResultInfo()->dofNames;
    for( UInt outDof = 0 ; outDof < dofNames_.GetSize(); outDof++ ) {
      
      // obtain new handle
      rHandles_[outDof] = mParser_->GetNewHandle();
      if( input_->GetEntryType() == EntryType::COMPLEX ) {
        iHandles_[outDof] = mParser_->GetNewHandle();
      }
      
      // register all input dofs 
      // a) real case
      if( input_->GetEntryType() == EntryType::DOUBLE ) {
        for( UInt inDof = 0; inDof < inDofNames.GetSize(); inDof++ ) {
          mParser_->SetValue( rHandles_[outDof], "u"+inDofNames[inDof], 0 );      
        }
        // register dofNames with dummy variable and set expression
        mParser_->SetExpr( rHandles_[outDof], rFuncs_[outDof] );
      } else {
        // b) complex case
        for( UInt inDof = 0; inDof < inDofNames.GetSize(); inDof++ ) {
          mParser_->SetValue( rHandles_[outDof], "u"+inDofNames[inDof]+"_real", 0 );      
          mParser_->SetValue( iHandles_[outDof], "u"+inDofNames[inDof]+"_real", 0 );  
          mParser_->SetValue( rHandles_[outDof], "u"+inDofNames[inDof]+"_imag", 0 );      
          mParser_->SetValue( iHandles_[outDof], "u"+inDofNames[inDof]+"_imag", 0 );  
        }
        // register dofNames with dummy variable and set expression
        if(  rFuncs_[outDof] == "" ) 
          rFuncs_[outDof] = "0"; 
        if(  iFuncs_[outDof] == "" ) 
          iFuncs_[outDof] = "0"; 
        mParser_->SetExpr( rHandles_[outDof], rFuncs_[outDof] );
        mParser_->SetExpr( iHandles_[outDof], iFuncs_[outDof] );
      }
    }

    // Create new result object
    // a) create new ResultInfo object for new result
    shared_ptr<ResultInfo> actInfo = 
      shared_ptr<ResultInfo>( new ResultInfo( *input_->GetResultInfo() ) );
    actInfo->resultType = NO_SOLUTION_TYPE;
    actInfo->resultName = resultName_;
    actInfo->unit = unit_;
    actInfo->dofNames = dofNames_;
    if( dofNames_.GetSize() > 1 ) {
      actInfo->entryType = ResultInfo::VECTOR;
    } else {
      actInfo->entryType = ResultInfo::SCALAR;    
    }

    // b) Create new result object and set it as output
    shared_ptr<BaseResult> newResult;
    if( input_->GetEntryType() == EntryType::DOUBLE ) {
      newResult = shared_ptr<BaseResult>(new Result<Double>() );
    } else {
      newResult = shared_ptr<BaseResult>( new Result<Complex>() );
    }
    newResult->SetResultInfo( actInfo );
    newResult->SetEntityList( input_->GetEntityList() );
    
    
    output_ = newResult;

  }
  
  
  void PostProcFunc::Initialize( const std::string& resultName,
                                 const std::string& unit,
                                 const StdVector<std::string>& dofNames,
                                 const StdVector<std::string>& rFunctions,
                                 const StdVector<std::string>& iFunctions ) {
    ENTER_FCN( "PostProcFunc::Initialize" );
    
    // ensure, that dofNames has the same size as functions
    assert( dofNames.GetSize() == rFunctions.GetSize() );
    assert( dofNames.GetSize() == iFunctions.GetSize() );
    
    // store resultName, dofNames and functions
    resultName_ = resultName;
    unit_ = unit;
    dofNames_ = dofNames;
    rFuncs_ = rFunctions;
    iFuncs_ = iFunctions;
  }

  
  void PostProcFunc::Apply( ) {
    ENTER_FCN( "PostProcFunc::Apply" );

    if( output_->GetEntryType() == EntryType::DOUBLE ) {
      ApplyReal();
    } else {
      ApplyComplex();
    }
  }
  
  void PostProcFunc::ApplyReal( ) {
    ENTER_FCN( "PostProcFunc::ApplyReal" );

    // Cast output into correct type
    Result<Double> & out = dynamic_cast<Result<Double>&>(*output_);
    Result<Double> & in =  dynamic_cast<Result<Double>&>(*input_ );
    Vector<Double> & outVec = out.GetVector();
    Vector<Double> & inVec = in.GetVector();
    UInt numOutDofs = out.GetResultInfo()->dofNames.GetSize();
    UInt numInDofs = in.GetResultInfo()->dofNames.GetSize();
    UInt numEnts = out.GetEntityList()->GetSize();
    EntityIterator outIt = out.GetEntityList()->GetIterator();
    StdVector<std::string> & inDofNames = input_->GetResultInfo()->dofNames;
    

    // intialize output vector
    outVec.Resize( numEnts * numOutDofs );
    
    // Iterate over all output entites
    for( ; !outIt.IsEnd(); outIt++ ) {

      // Iterate over all output dofs
      for( UInt outDof = 0; outDof < numOutDofs; outDof++ ) {
        
        // If output entity is a node 
        // -> register coordinate as px/py/pz
        if( outIt.GetType() == EntityList::NODE_LIST) {
          Vector<Double> coords;
          ptGrid_->GetNodeCoordinate( coords, outIt.GetNode() );
          mParser_->SetCoordinates(  rHandles_[outDof],
                                     *(domain->GetCoordSystem()),
                                     coords );
        }
        
        // Register current input values for dofNames
        for( UInt inDof = 0; inDof < numInDofs; inDof++ ) {
          Double actVal = inVec[outIt.GetPos()*numInDofs + inDof];
          mParser_->SetValue( rHandles_[outDof], 
                              "u"+inDofNames[inDof], actVal );
        }
        // Apply function
        outVec[outIt.GetPos()*numOutDofs+outDof] = 
          mParser_->Eval( rHandles_[outDof]  );
        
      } // loop ober output dofs
    } // loop over entities    
  }
  
  void PostProcFunc::ApplyComplex( ) {
    ENTER_FCN( "PostProcFunc::ApplyComplex" );

    // Cast output into correct type
    Result<Complex> & out = dynamic_cast<Result<Complex>&>(*output_);
    Result<Complex> & in =  dynamic_cast<Result<Complex>&>(*input_ );
    Vector<Complex> & outVec = out.GetVector();
    Vector<Complex> & inVec = in.GetVector();
    UInt numOutDofs = out.GetResultInfo()->dofNames.GetSize();
    UInt numInDofs = in.GetResultInfo()->dofNames.GetSize();
    UInt numEnts = out.GetEntityList()->GetSize();
    EntityIterator outIt = out.GetEntityList()->GetIterator();
    StdVector<std::string> & inDofNames = input_->GetResultInfo()->dofNames;
    
    // intialize output vector
    outVec.Resize( numEnts * numOutDofs );
    
    // Iterate over all output entites
    for( ; !outIt.IsEnd(); outIt++ ) {

      // Iterate over all output dofs
      for( UInt outDof = 0; outDof < numOutDofs; outDof++ ) {
        
        // If output entity is a node 
        // -> register coordinate as px/py/pz
        if( outIt.GetType() == EntityList::NODE_LIST) {
          Vector<Double> coords;
          ptGrid_->GetNodeCoordinate( coords, outIt.GetNode() );
          mParser_->SetCoordinates( rHandles_[outDof],
                                    *(domain->GetCoordSystem()),
                                    coords );
          mParser_->SetCoordinates( iHandles_[outDof],
                                    *(domain->GetCoordSystem()),
                                    coords );
        }
        
        // Register current input values for dofNames
        for( UInt inDof = 0; inDof < numInDofs; inDof++ ) {
          Complex actVal = inVec[outIt.GetPos()*numInDofs + inDof];
          mParser_->SetValue( rHandles_[outDof], 
                              "u"+inDofNames[inDof]+"_real", actVal.real() );
          mParser_->SetValue( rHandles_[outDof], 
                              "u"+inDofNames[inDof]+"_imag", actVal.imag() );
          mParser_->SetValue( iHandles_[outDof], 
                              "u"+inDofNames[inDof]+"_real", actVal.real() );
          mParser_->SetValue( iHandles_[outDof], 
                              "u"+inDofNames[inDof]+"_imag", actVal.imag() );
        }
        // Apply function
        Double real, imag;
        real = mParser_->Eval( rHandles_[outDof] );
        imag = mParser_->Eval( iHandles_[outDof] );
        outVec[outIt.GetPos()*numOutDofs+outDof] = Complex(real, imag);
        
      } // loop ober output dofs
    } // loop over entities    

  }

  void PostProcFunc::Finalize( ) {
    ENTER_FCN( "PostProcFunc::Finalize" );
    
  }
  
  
}
