// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "PostProc.hh"

#include "Utils/mathParser/mathParser.hh"
#include "Domain/Results/BaseResults.hh"
#include "Domain/Domain.hh"
#include "Domain/Results/ResultInfo.hh"
#include "Domain/Mesh/Grid.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "FeBasis/FeFunctions.hh"
#include "PDE/SinglePDE.hh"
#include "CoupledPDE/IterCoupledPDE.hh"
#include "ResultHandler.hh"
#include "Utils/ToolsFull.hh"
#include "Domain/CoefFunction/CoefXpr.hh"
#include "Domain/CoefFunction/CoefFunctionCompound.hh"

namespace CoupledField {



  PostProc::PostProc( Grid * ptGrid, PtrParamNode postProcNode ) {

    ptGrid_ = ptGrid;
    myParam_ = postProcNode;
    reducType_ = NONE;
    writeResult_ = true;
    isHistory_ = false;

  }

  PostProc::~PostProc() {
  }

  void PostProc::GetOutDestNames( StdVector<std::string> & outNames ) {
    outNames = outputNames_;
  }

  bool PostProc::IsWriteResult( ) {
    return writeResult_;
  }
  
  bool PostProc::IsHistory( ) {
      return isHistory_;
    }

  void PostProc::GetReducedList( shared_ptr<EntityList>& outList,
                                  ResultInfo::EntityUnknownType& outType,
                                  const shared_ptr<EntityList> inList,
                                  const ResultInfo::EntityUnknownType inType,
                                  const ReductionType reducType ) {

    shared_ptr<EntityList> retList;
    
    if( reducType == NONE ||
        reducType == TIME_FREQ ) {
      outList = inList;
      outType = inType;
    } 

    if( reducType == SPACE ) {

        isHistory_ = true;
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
                                            inList->GetName() );
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
                                            inList->GetName() );

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

  void PostProc::CreatePostProc( PtrParamNode procNode,
                                 Grid * ptGrid, 
                                 StdVector<shared_ptr<PostProc> > & postProcs ) {

    // fetch single postProcs defined in this section
    ParamNodeList procNodes = procNode->GetChildren();

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
        procNodes[i]->GetValue( "reduction", typeString );
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
        procNodes[i]->GetValue( "reduction", typeString );
        if( typeString == "space" ) {
          type = SPACE;
        } else {
          type = TIME_FREQ;
        }
        actProc = shared_ptr<PostProcMax>(new PostProcMax( ptGrid, type,PtrParamNode() ) );
      } 
      
      
      // === FUNCTION ===
      if( procName == "function" ) {      
        found = true;
        shared_ptr<PostProcFunc> funcProc;
        funcProc = shared_ptr<PostProcFunc>(new PostProcFunc( ptGrid, PtrParamNode() ) );

        // get resultName
        std::string resultName;
        procNodes[i]->GetValue( "resultName", resultName );

        // get unit
        std::string unit;
        procNodes[i]->GetValue( "unit", unit );

        // get inputResult node
        ParamNodeList inputNodes = procNodes[i]->GetList("inputResult");
        
        StdVector<std::string> variableNames;
        for( UInt u = 0; u < inputNodes.GetSize(); u++ ) {
          variableNames.Push_back( inputNodes[u]->Get("variableName")->As<std::string>() );
        }

        // get scalar/vector node
        PtrParamNode scalarNode = procNodes[i]->Get("scalar",ParamNode::PASS);
        PtrParamNode vectorNode = procNodes[i]->Get("vector",ParamNode::PASS);

        ParamNodeList dofNodes;
        bool isScalar = false;
        // we can perform this simple check since scalar and vector are exclusive in the XML scheme
        if( scalarNode ) {
          dofNodes = scalarNode;
          isScalar = true;
        } else if( vectorNode ) {
          dofNodes = vectorNode->GetList("comp");
        } else {
          EXCEPTION("No definition of scalar/vectorial postprocessing found!");
        }

        
        // get dofnames, real functions and imaginary functions for each dof
        StdVector<std::string> dofNames, rFunctions, iFunctions;
        for( UInt iDof = 0; iDof < dofNodes.GetSize(); iDof++ ) {
          if( isScalar ) {
            // for scalar values we don't have to (and can't) read the DoF name and just set it to "".
            dofNames.Push_back("");
          } else {
            dofNames.Push_back( dofNodes[iDof]->Get("dof")->As<std::string>() );
          }
          rFunctions.Push_back( dofNodes[iDof]->Get("realFunc")->As<std::string>() );
          iFunctions.Push_back( dofNodes[iDof]->Get("imagFunc")->As<std::string>() );
        }
        funcProc->Initialize( resultName, unit, dofNames, rFunctions, iFunctions, variableNames, procNodes[i] );
        actProc = funcProc;
      }

      // === L2Norm ===
      if( procName == "L2Norm" ) {      
        found = true;
        shared_ptr<PostProcL2Norm> normProc;
        normProc = shared_ptr<PostProcL2Norm>(new PostProcL2Norm( ptGrid, PtrParamNode() ) );

        // get resultName
        std::string resultName;
        procNodes[i]->GetValue( "resultName", resultName );

        // get integration order
        Integer integrationOrder;
        procNodes[i]->GetValue( "integrationOrder", integrationOrder );

        // get mode
        std::string mode;
        procNodes[i]->GetValue( "mode", mode );

        // get computeSolutionSteps
        std::string computeSolutionSteps;
        procNodes[i]->GetValue( "computeSolutionSteps", computeSolutionSteps );

        // get dof-nodes
        ParamNodeList dofNodes = procNodes[i]->GetList("dof");

        // get dofnames, real functions and imaginary functions for each dof
        StdVector<std::string> dofNames, rFunctions, iFunctions;
        for( UInt iDof = 0; iDof < dofNodes.GetSize(); iDof++ ) {
          dofNames.Push_back( dofNodes[iDof]->Get("name")->As<std::string>() );
          rFunctions.Push_back( dofNodes[iDof]->Get("realFunc")->As<std::string>() );
          iFunctions.Push_back( dofNodes[iDof]->Get("imagFunc")->As<std::string>() );
        }
        normProc->Initialize( resultName, integrationOrder, mode, computeSolutionSteps, dofNames, rFunctions, iFunctions );
        actProc = normProc;
      }

      // === LIMIT ===
      if( procName == "limit" ) {
        found = true;

        shared_ptr<PostProcLimit> limitProc;
        limitProc = shared_ptr<PostProcLimit>(new PostProcLimit( ptGrid, PtrParamNode() ) );

        // set upper and lower value, if present
        if( procNodes[i]->Has("lowerLimit") ){
          limitProc->SetLowerLimit( procNodes[i]->Get("lowerLimit")->As<Double>() );
        }
        if( procNodes[i]->Has("upperLimit") ){
          limitProc->SetUpperLimit( procNodes[i]->Get("upperLimit")->As<Double>() );
        }
        actProc = limitProc;

      } 

      // Check, if postProc method was found
      if (!found) {
        EXCEPTION( "Postprocessing method with name '" << procName
                   << "' is not known." );
      }
      
      // get next postprocessing routine
      procNodes[i]->GetValue( "postProcId", actProc->next_ );

      // get writeResult
      procNodes[i]->GetValue( "writeResult", actProc->writeResult_);

      // get outputDestNames
      std::string outNames;
      procNodes[i]->GetValue( "outputIds", outNames );
      SplitStringList( outNames, actProc->outputNames_, ',' );
      postProcs.Push_back( actProc );
    }

  }

  // =================== SUM-POSTPROCEDURE ===================
  
  PostProcSum::PostProcSum( Grid * ptGrid, ReductionType type, 
                            PtrParamNode postProcNode )
    : PostProc( ptGrid, postProcNode ) {
    
    name_ = "sum";
    reducType_ = type;
  }

  PostProcSum::~PostProcSum() { }
  
  template<class TYPE> 
  void PostProcSum::CalcSum() {

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
  
  void PostProcSum::SetResult( shared_ptr<BaseResult> res ) {
    
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
    if( res->GetEntryType() == BaseMatrix::DOUBLE ) {
      newResult = shared_ptr<BaseResult>(new Result<Double>() );
    } else {
      newResult = shared_ptr<BaseResult>( new Result<Complex>() );
    }
    newResult->SetEntityList( newList );
    newResult->SetResultInfo( newInfo );
    output_ = newResult;
  }

  void PostProcSum::Apply() {
    
    if( output_->GetEntryType() == BaseMatrix::DOUBLE ) {
      CalcSum<Double>();
    } else {
      CalcSum<Complex>();
    }
  }


  void PostProcSum::Finalize() {  }

  // =================== MAXIMUM-POSTPROCEDURE ===================

  PostProcMax::PostProcMax( Grid * ptGrid, ReductionType type,
                            PtrParamNode postProcNode )
    : PostProc( ptGrid, postProcNode ) {
    
    name_ = "max";
    reducType_ = type;
  }

  PostProcMax::~PostProcMax() {
    
  }

  void PostProcMax::SetResult( shared_ptr<BaseResult> res ) {
    
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
    if( res->GetEntryType() == BaseMatrix::DOUBLE ) {
      newResult = shared_ptr<BaseResult>(new Result<Double>() );
    } else {
      newResult = shared_ptr<BaseResult>( new Result<Complex>() );
    }
    newResult->SetEntityList( newList );
    newResult->SetResultInfo( newInfo );
    output_ = newResult;

  }

  void PostProcMax::Apply() {
    
    if( output_->GetEntryType() == BaseMatrix::DOUBLE ) {
      CalcMax<Double>();
    } else {
      CalcMax<Complex>();
    }
  }


  template<class TYPE> 
  void PostProcMax::CalcMax() {
    
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
                              PtrParamNode postProcNode )
    : PostProc( ptGrid, postProcNode ) {

    name_ = "func";
    reducType_ = NONE;

    // fetch math parser object
    mParser_ = domain->GetMathParser();
  }

  PostProcFunc::~PostProcFunc() {
    
  }

  void PostProcFunc::SetResult( shared_ptr<BaseResult> res ) {

    // Set result as input
    input_ =  res;

    // check if we only have one variable name (currently only one is supported, the index is hardcoded)
    UInt i_var = 0;
    if( variableNames_.size()!=1 ) {
      EXCEPTION("PostProcFunc: More than one variable name specified, please only specify one variable name for your result evaluation!");
    }

    // intialize rVals and mathParser handles
    // a) double
    rVals_.Resize( dofNames_.GetSize() );
    rVals_.Init();
    rHandles_.Resize( dofNames_.GetSize() );
    // b) complex
    if( input_->GetEntryType() == BaseMatrix::COMPLEX ) {
      iVals_.Resize( dofNames_.GetSize() );
      iVals_.Init();
      iHandles_.Resize( dofNames_.GetSize() );
    }
    
    // register all new functions with mathParser and new result info  object
    StdVector<std::string> & inDofNames = input_->GetResultInfo()->dofNames;
    for( UInt outDof = 0 ; outDof < dofNames_.GetSize(); outDof++ ) {
      
      // obtain new handle
      rHandles_[outDof] = mParser_->GetNewHandle(true);
      if( input_->GetEntryType() == BaseMatrix::COMPLEX ) {
        iHandles_[outDof] = mParser_->GetNewHandle(true);
      }
      
      // register all input dofs 
      // a) real case
      if( input_->GetEntryType() == BaseMatrix::DOUBLE ) {
        for( UInt inDof = 0; inDof < inDofNames.GetSize(); inDof++ ) {
          mParser_->SetValue( rHandles_[outDof], variableNames_[i_var]+inDofNames[inDof], 0 );      
        }
        // register dofNames with dummy variable and set expression
        mParser_->SetExpr( rHandles_[outDof], rFuncs_[outDof] );
      } else {
        // b) complex case
        for( UInt inDof = 0; inDof < inDofNames.GetSize(); inDof++ ) {
          mParser_->SetValue( rHandles_[outDof], variableNames_[i_var]+inDofNames[inDof]+"_real", 0 );      
          mParser_->SetValue( iHandles_[outDof], variableNames_[i_var]+inDofNames[inDof]+"_real", 0 );  
          mParser_->SetValue( rHandles_[outDof], variableNames_[i_var]+inDofNames[inDof]+"_imag", 0 );      
          mParser_->SetValue( iHandles_[outDof], variableNames_[i_var]+inDofNames[inDof]+"_imag", 0 );  
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

    // check if a coefFunction shall be defined. This will be set to true if parsing an enum is successfull
    bool defineCoefFunc = false;

    // Create new result object
    // a) create new ResultInfo object for new result
    shared_ptr<ResultInfo> actInfo = 
      shared_ptr<ResultInfo>( new ResultInfo( *input_->GetResultInfo() ) );
    
    // If we can parse the result it was already set by the initial definition of the back-coupling (done in SinglePDE::InitStage2 SetBCs())
    // Hence, we have already checked that the result is available in the postProcList.
    // If this is not the case, we just have a standard post-processing result and we can skip the definition of a coefFunction
    
    SolutionType solType;
    try 
    {
      solType = SolutionTypeEnum.Parse(resultName_);
      actInfo->resultType = solType;
      defineCoefFunc = true;
    }  catch ( const Exception& e)
    {
      // parsing was not successful, fall back to NO_SOLUTION_TYPE
      solType = NO_SOLUTION_TYPE;
      actInfo->resultType = solType;
    }  

    actInfo->resultName = resultName_;
    actInfo->dofNames = dofNames_;
    actInfo->unit = unit_;
    if( dofNames_.GetSize() > 1 ) {
      actInfo->entryType = ResultInfo::VECTOR;
    } else {
      actInfo->entryType = ResultInfo::SCALAR;    
    }

    // b) Create new result object and set it as output
    shared_ptr<BaseResult> newResult;
    if( input_->GetEntryType() == BaseMatrix::DOUBLE ) {
      newResult = shared_ptr<BaseResult>( new Result<Double>() );
    } else {
      newResult = shared_ptr<BaseResult>( new Result<Complex>() );
    }
    newResult->SetResultInfo( actInfo );
    newResult->SetEntityList( input_->GetEntityList() );
    
    
    output_ = newResult;


    // c) if we want to use this result for backcoupling, we also need to make it available for the singlePDE in GetCoefFct() via fieldCoefs_
    // additionally, we need to have unique identifier, so this will only be triggered when the identifier is different
    if( defineCoefFunc ) {
      // the complex valued case is not implemented yet
      if( input_->GetEntryType() != BaseMatrix::DOUBLE ) {
        EXCEPTION("PostProcFunc: CoefFunction definitions for complex valued input is not implemented yet! Only real valued functions are allowed for back-coupling!");
      }

      // get coefFunctions of primary variable
      shared_ptr<BaseFeFunction> feFct = res->GetResultInfo()->GetFeFunction().lock();
      PtrCoefFct coefPDE = feFct->GetPDE()->GetCoefFct( res->GetResultInfo()->resultType );

      // sort Funcs
      std::string probGeo = feFct->GetPDE()->GetDomain()->GetParamRoot()->Get("domain")->Get("geometryType")->As<std::string>();
      
      shared_ptr<ResultInfo> resInfo = feFct->GetResultInfo();
      //const StdVector<std::string>& dofNamesRes = resInfo->dofNames;
      
      if( dofNames_.GetSize() > 1 ) {
        StdVector<std::string> dofNamesGrid;
        if ( probGeo == "3d" ) {
          dofNamesGrid = "x", "y", "z";
        }
        else if ( probGeo == "plane")
        {
          dofNamesGrid = "x", "y";
        }
        else if ( probGeo == "axi" ) {
          EXCEPTION("DefineCoefFunction for postProcFunc is not implemented for axi-symmetric simulations, please implement me!");
          dofNamesGrid = "r", "z";
        }

        // now loop over DOFs and match rFuncs and iFuncs accordingly
        for( UInt i = 0; i < dofNamesGrid.GetSize(); i++  ) {
          for( UInt o = 0; o < dofNames_.GetSize(); o++  ) {
            if( dofNamesGrid[i]==dofNames_[o] ) {
              rFuncsSorted_[i] = rFuncs_[o];
              //iFuncsSorted_[i] = iFuncs_[o];
            }
          }
        }
      } else {
        // scalar problem
        rFuncsSorted_[0] = rFuncs_[0];
        //iFuncsSorted_[0] = iFuncs_[0];
      }

      // check scalar/vector

      // we use a coefFunctionCompound to replace the value of u in the function with the actual feFunction
      // create vector expression similar to description in CoefFunctionCompound.hh
      //! std::string vecXpr = "(a_0_R * b_R), (a_1_R * b_R);"
      //! c.SetVector( vecXpr, vars);

      // problem: we have in inconsistent naming convetion (u_0_R vs u0_real) --> parse and convert
      // define the strings the user uses
      StdVector<std::string> userStr;
      StdVector<std::string> mpStr;
      StdVector<std::string> checkStr;
      userStr = variableNames_[i_var], variableNames_[i_var]+"x", variableNames_[i_var]+"y", variableNames_[i_var]+"z";
      mpStr = variableNames_[i_var]+"_R", variableNames_[i_var]+"_0_R", variableNames_[i_var]+"_1_R", variableNames_[i_var]+"_2_R"; // this convention is added/used automatically by the mathParser!
      checkStr = "x", "y", "z";
      
      rFuncsSortedMp_ = rFuncsSorted_; // sorted and adapted rFuncs for the mathParser/coefFunctionCompound
      // loop over all dofs (functions associated with one DoF)
      for( UInt i = 0; i < dofNames_.GetSize(); i++  ) {
        // loop over all possible user strigns (feFunction)
        for( UInt u = 0; u < userStr.GetSize(); u++ ) {
          size_t index = 0;
          UInt curLen = userStr[u].length();
          // go through the string and replace every occurence of the string
          while (true) {
            // search for substring and replace it
            index = rFuncsSortedMp_[i].find(userStr[u], index);
            if (index == std::string::npos) break;

            // in the first loop we will find all "u", but it could be a "u0" to "u2" instead of just an "u",
            // which would lead to a wrong replacement. Hence, we check if the character at the next position
            // is different from 0 to 2, ensuring we are only replacing "u".
            if ( u==0 ) {
              if ( rFuncsSortedMp_[i][index+1] != checkStr[0][0] &&
                  rFuncsSortedMp_[i][index+1] != checkStr[1][0] && 
                  rFuncsSortedMp_[i][index+1] != checkStr[2][0] ) {
                rFuncsSortedMp_[i].replace(index, curLen, mpStr[u]);
              }
            } else {
              rFuncsSortedMp_[i].replace(index, curLen, mpStr[u]);
            }

            // advance index according to the length of the replacement string
            index += curLen;
          }
        }
      } 

      // now we have prepared the correct syntax, but in order to use multiple results at the same time, we have to register
      // them in a unique manner. Hence, we add the number of the current result after each "u".
      std::string solString = SolutionTypeEnum.ToString(solType);
      std::vector<std::string> substrings;
      std::stringstream ss(solString);
      std::string token;
      char delimiter = '_'; // it is mandatory to keep the convention of the naming for generic results as it is, otherwise this will fail
      std::string resNr;

      // Split the string by the delimiter
      while (std::getline(ss, token, delimiter)) {
          substrings.push_back(token);
      }

      // Return the last substring (our result number)
      if (!substrings.empty()) {
          resNr = substrings.back();
      } else {
        EXCEPTION("PostProc: Could not identify result number of generic result");
      }

      // insert the result number to have unique identifiers
      for( UInt i = 0; i < dofNames_.GetSize(); i++  ) {
        size_t index = 0;
        while (true) {
          // search for substring and replace it
          index = rFuncsSortedMp_[i].find(userStr[0], index);
          if (index == std::string::npos) break;

          // insert the result number after each occurrence of "u"
          rFuncsSortedMp_[i].insert(index+1, resNr);

          // advance index according to the length of the replacement string
          index += 2;
        }
      } 

      // debug output of the functions which will be evaluated (as a string)
      //for( UInt i = 0; i < dofNames_.GetSize(); i++  ) {
      //  std::cout << rFuncsSortedMp_[i] << std::endl;
      //}

      // define a coefFunction which can handle the combination of a feFunction with mathParser expressions 
      // similar to mechPDE --> e.g. MECHANIC MINIMUM PRINCIPAL STRESS - VECTOR
      shared_ptr<CoefFunctionCompound<Double> > coefFunc(new CoefFunctionCompound<Double>(mParser_));
      std::map<std::string,PtrCoefFct> symbolsFe;
      std::string curSymb = userStr[0] + resNr;
      symbolsFe[curSymb] = coefPDE;

      if( dofNames_.GetSize() > 1 ) {
        // vectorial quantity
        coefFunc->SetVector(rFuncsSortedMp_,symbolsFe);
      } else {
        // scalar quantity
        coefFunc->SetScalar(rFuncsSortedMp_[0],symbolsFe);
      }
      
      // set the coefFunction to the PDE so that SinglePDE::GetCoefFct() can now access our new result
      // fieldCoefs_ is not directly accessible, we have to use an additional set method
      feFct->GetPDE()->SetFieldCoef( coefFunc, actInfo->resultType );
    }
  }
  
  
  void PostProcFunc::Initialize( const std::string& resultName,
                                 const std::string& unit,
                                 const StdVector<std::string>& dofNames,
                                 const StdVector<std::string>& rFunctions,
                                 const StdVector<std::string>& iFunctions,
                                 const StdVector<std::string>& variableNames,
                                 const PtrParamNode& postProcNode ) {
    
    // ensure, that dofNames has the same size as functions
    assert( dofNames.GetSize() == rFunctions.GetSize() );
    assert( dofNames.GetSize() == iFunctions.GetSize() );
    
    // store resultName, dofNames and functions
    resultName_ = resultName;
    unit_ = unit;
    dofNames_ = dofNames;
    rFuncs_ = rFunctions;
    iFuncs_ = iFunctions;
    variableNames_ = variableNames;
    rFuncsSorted_.Resize(rFuncs_.GetSize());
    //iFuncsSorted_.Resize(iFuncs_.GetSize());
    postProcNode_ = postProcNode;
  }

  
  void PostProcFunc::Apply( ) {

    if( output_->GetEntryType() == BaseMatrix::DOUBLE ) {
      ApplyReal();
    } else {
      ApplyComplex();
    }
  }
  
  void PostProcFunc::ApplyReal( ) {

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
    
    // check if we only have one variable name (currently only one is supported, the index is hardcoded)
    UInt i_var = 0;

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
          ptGrid_->GetNodeCoordinate( coords, outIt.GetNode(), true );
          mParser_->SetCoordinates(  rHandles_[outDof],
                                     *(domain->GetCoordSystem()),
                                     coords );
        } else if (outIt.GetType() == EntityList::ELEM_LIST) {
          
          Vector<Double> globMidPoint;
          shared_ptr<ElemShapeMap> esm =
              ptGrid_->GetElemShapeMap(outIt.GetElem(), true);
          esm->GetGlobMidPoint(globMidPoint);
          mParser_->SetCoordinates(  rHandles_[outDof],
                                     *(domain->GetCoordSystem()),
                                     globMidPoint );
        }
        
        // Register current input values for dofNames
        for( UInt inDof = 0; inDof < numInDofs; inDof++ ) {
          Double actVal = inVec[outIt.GetPos()*numInDofs + inDof];
          mParser_->SetValue( rHandles_[outDof], 
                              variableNames_[i_var]+inDofNames[inDof], actVal );
        }
        // Apply function
        outVec[outIt.GetPos()*numOutDofs+outDof] = 
          mParser_->Eval( rHandles_[outDof]  );
        
      } // loop ober output dofs
    } // loop over entities    
  }
  
  void PostProcFunc::ApplyComplex( ) {

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

    // check if we only have one variable name (currently only one is supported, the index is hardcoded)
    UInt i_var = 0;
    
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
          ptGrid_->GetNodeCoordinate( coords, outIt.GetNode(), true );
          mParser_->SetCoordinates( rHandles_[outDof],
                                    *(domain->GetCoordSystem()),
                                    coords );
          mParser_->SetCoordinates( iHandles_[outDof],
                                    *(domain->GetCoordSystem()),
                                    coords );
        } else if (outIt.GetType() == EntityList::ELEM_LIST) {
          Vector<Double> globMidPoint;
          shared_ptr<ElemShapeMap> esm =
              ptGrid_->GetElemShapeMap(outIt.GetElem(), true);
          esm->GetGlobMidPoint(globMidPoint);
          mParser_->SetCoordinates(  rHandles_[outDof],
                                     *(domain->GetCoordSystem()),
                                     globMidPoint );
        }

        // Register current input values for dofNames
        for( UInt inDof = 0; inDof < numInDofs; inDof++ ) {
          Complex actVal = inVec[outIt.GetPos()*numInDofs + inDof];
          mParser_->SetValue( rHandles_[outDof], 
                              variableNames_[i_var]+inDofNames[inDof]+"_real", actVal.real() );
          mParser_->SetValue( rHandles_[outDof], 
                              variableNames_[i_var]+inDofNames[inDof]+"_imag", actVal.imag() );
          mParser_->SetValue( iHandles_[outDof], 
                              variableNames_[i_var]+inDofNames[inDof]+"_real", actVal.real() );
          mParser_->SetValue( iHandles_[outDof], 
                              variableNames_[i_var]+inDofNames[inDof]+"_imag", actVal.imag() );
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
    
  }

// =================== L2Norm-POSTPROCEDURE ===================
  

  PostProcL2Norm::PostProcL2Norm( Grid * ptGrid,
                              PtrParamNode postProcNode )
    : PostProc( ptGrid, postProcNode ) {

    name_ = "L2Norm";
    reducType_ = SPACE;

    // fetch math parser object
    mParser_ = domain->GetMathParser();

  }

  PostProcL2Norm::~PostProcL2Norm() {
    
  }

  void PostProcL2Norm::SetResult( shared_ptr<BaseResult> res ) {

    // Append result to list of results
    input_ =  res;

    
    // intialize rVals and mathParser handles
    // a) double
    rVals_.Resize( dofNames_.GetSize() );
    rVals_.Init();
    rHandles_.Resize( dofNames_.GetSize() );
    // b) complex
    if( input_->GetEntryType() == BaseMatrix::COMPLEX ) {
      iVals_.Resize( dofNames_.GetSize() );
      iVals_.Init();
      iHandles_.Resize( dofNames_.GetSize() );
    }

    if( !(input_->GetResultInfo()->entryType == ResultInfo::SCALAR) && !(input_->GetResultInfo()->entryType == ResultInfo::VECTOR) ) {
      EXCEPTION("Only scalar and vectorial quantities are supported at the moment!");
    }


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
    if( mode_ == "relative" ) {
      newInfo->resultName += "RelL2Norm";
    } else {
      newInfo->resultName += "L2Norm";
    }
    if( mode_ == "relative" ) {
      newInfo->unit = "";
    } else {
      // check grid dim in order to set the correct unit
      UInt dim = res->GetResultInfo()->GetFeFunction().lock()->GetGrid()->GetDim();
      std::string oldUnit = newInfo->unit;
      if( dim == 2 ) {
        // plane
        newInfo->unit = "((" + oldUnit + ")^2 m^2)^(1/2)";
      } else if ( dim == 3) {
        // 3D
        newInfo->unit = "((" + oldUnit + ")^2 m^3)^(1/2)";
      } else {
        EXCEPTION("Unknown grid dimension!");
      }
    }
    // we calcualte the norm, this will always be scalar regardless of the input (scalar/vector)
    newInfo->dofNames = "";
    newInfo->entryType = ResultInfo::SCALAR;
    
    
    // c) Create new Result and set is as output
    shared_ptr<BaseResult> newResult;
    if( input_->GetEntryType() == BaseMatrix::COMPLEX ) {
      newResult = shared_ptr<BaseResult>(new Result<Complex>() );
    } else {
      newResult = shared_ptr<BaseResult>( new Result<Double>() );
    }
    newResult->SetEntityList( newList );
    newResult->SetResultInfo( newInfo );
    
    output_ = newResult;


    // Adapt functor for result evaluation

    shared_ptr<BaseFeFunction> feFct = res->GetResultInfo()->GetFeFunction().lock();
    // get coefFunctions of primary variable
    PtrCoefFct coefPDE = feFct->GetPDE()->GetCoefFct( res->GetResultInfo()->resultType );

    // sort Funcs
    std::string probGeo = feFct->GetPDE()->GetDomain()->GetParamRoot()->Get("domain")->Get("geometryType")->As<std::string>();
    
    shared_ptr<ResultInfo> resInfo = feFct->GetResultInfo();
    const StdVector<std::string>& dofNamesRes = resInfo->dofNames;
    
    if( dofNamesRes.GetSize() > 1 ) {
      StdVector<std::string> dofNamesGrid;
      if ( probGeo == "3d" ) {
        dofNamesGrid = "x", "y", "z";
      }
      else if ( probGeo == "plane")
      {
        dofNamesGrid = "x", "y";
      }
      else if ( probGeo == "axi" ) {
        EXCEPTION("L2Norm is not implemented for axi-symmetric simulations, please implement me!");
        dofNamesGrid = "r", "z";
      }

      // now loop over DOFs and match rFuncs and iFuncs accordingly
      for( UInt i = 0; i < dofNamesGrid.GetSize(); i++  ) {
        for( UInt o = 0; o < dofNames_.GetSize(); o++  ) {
          if( dofNamesGrid[i]==dofNames_[o] ) {
            rFuncsSorted_[i] = rFuncs_[o];
            iFuncsSorted_[i] = iFuncs_[o];
          }
        }
      }
    } else {
      // scalar problem, check if only one dof is given and assign
      if( dofNames_.GetSize()==1 ) {
        rFuncsSorted_[0] = rFuncs_[0];
        iFuncsSorted_[0] = iFuncs_[0];
      } else {
        EXCEPTION("Multiple DOFs defined for a scalar problem, check your xml-file!");
      }
    }

    
    Global::ComplexPart part;
    bool isComplex;
    if( output_->GetEntryType() == BaseMatrix::COMPLEX ) {
      isComplex = true;
      part = Global::COMPLEX;
    } else {
      isComplex = false;
      part = Global::REAL;
    }

    PtrCoefFct coef;
    if(!isComplex) {
      coef = CoefFunction::Generate(mParser_, Global::REAL, rFuncsSorted_ );
    } else {
      coef = CoefFunction::Generate(mParser_, Global::COMPLEX, rFuncsSorted_, iFuncsSorted_ );
    }
    //! Generate coefFunction of the difference (analytical and numerical solution)
    PtrCoefFct diffCoef;
    diffCoef = CoefFunction::Generate( mParser_, part,
                CoefXprBinOp(mParser_,coefPDE,coef,CoefXpr::OP_SUB));

    this->coef_ = coef;
    this->diffCoef_ = diffCoef;

  }
  
  
  void PostProcL2Norm::Initialize( const std::string& resultName,
                                 const UInt& integrationOrder,
                                 const std::string& mode,
                                 const std::string& computeSolutionSteps,
                                 const StdVector<std::string>& dofNames,
                                 const StdVector<std::string>& rFunctions,
                                 const StdVector<std::string>& iFunctions ) {
    
    // ensure, that dofNames has the same size as functions
    assert( dofNames.GetSize() == rFunctions.GetSize() );
    assert( dofNames.GetSize() == iFunctions.GetSize() );
    
    // store resultName, dofNames and functions
    resultName_ = resultName;
    integrationOrder_ = integrationOrder;
    mode_ = mode;
    computeSolutionSteps_ = computeSolutionSteps;
    dofNames_ = dofNames;
    rFuncs_ = rFunctions;
    iFuncs_ = iFunctions;
    rFuncsSorted_.Resize(rFuncs_.GetSize());
    iFuncsSorted_.Resize(iFuncs_.GetSize());
  }

  
  void PostProcL2Norm::Apply( ) {

    if( output_->GetEntryType() == BaseMatrix::DOUBLE ) {
      CalcNorm<Double>();
    } else {
      CalcNorm<Complex>();
    }
  }
  
  template<class TYPE> 
  void PostProcL2Norm::CalcNorm() {
    // Cast output into correct type
    Result<TYPE> & out = dynamic_cast<Result<TYPE>&>(*output_);
    Vector<TYPE> & outVec = out.GetVector();

    shared_ptr<EntityList> outList = out.GetEntityList();
    EntityIterator outIt = outList->GetIterator();
    ResultInfo & outResInfo = *(out.GetResultInfo() );
    outVec.Resize( outResInfo.dofNames.GetSize() );
    outVec.Init();

    shared_ptr<BaseFeFunction> baseFeFct = output_->GetResultInfo()->GetFeFunction().lock();
    shared_ptr<FeFunction<TYPE> > feFct = dynamic_pointer_cast<FeFunction<TYPE> >(baseFeFct);

    // determine if we want to compute a solution for this specific time/frequency step
    ResultHandler* handler = feFct->GetPDE()->GetDomain()->GetResultHandler();
    bool solveThisStep = false;
    if (computeSolutionSteps_ == "all")
      solveThisStep = true;
    else if (computeSolutionSteps_ == "last" && handler->GetCurrStepNum() == handler->GetNumSteps())
      solveThisStep = true;

    // get lpm for higher integration order and loop over all lpms
    Vector<TYPE> elemSol, temp, fac, refFac;
    // Loop over regions
    for(outIt.Begin(); !outIt.IsEnd(); outIt++ ) {
      // if specified in the xml tag, only compute for desired time/frequency steps
      if (solveThisStep == true) {
        shared_ptr<EntityList> actSDList = outIt.GetGrid()->GetEntityList( EntityList::ELEM_LIST, outIt.GetGrid()->GetRegionName(outIt.GetRegion()) );

        // summation variables
        Double tempNorm = 0.0;
        Double refSum = 0.0;

        // we need to do conventional for loop as somehow openMp can not deal with the iterator
        int numElems =  actSDList->GetIterator().GetSize();
        int iElems;
        // loop over elements
        #pragma omp parallel for num_threads(CFS_NUM_THREADS) private(elemSol, temp, fac, refFac)
        for(iElems = 0; iElems < numElems; ++iElems) {
          // define element iterator
          EntityIterator elemIt = actSDList->GetIterator();
          // increment iterator to current position
          elemIt += iElems;
          // get the respective element
          const Elem * el = elemIt.GetElem();

          // we calcualte the individual contributions by coef^2*J*w
          // the feFunction is needed in order to get the grid and integration stuff
          feFct->GetElemSolution( elemSol, el);

          // Obtain FE element from feSpace and integration scheme
          IntegOrder order;
          // Use higher order integration in order to get accurate results
          order.SetIsoOrder(integrationOrder_);
          shared_ptr<FeSpace> feSpace = feFct->GetFeSpace();
          shared_ptr<IntScheme> intScheme = feSpace->GetIntScheme();
          // Get shape map from grid
          shared_ptr<ElemShapeMap> esm = elemIt.GetGrid()->GetElemShapeMap( el, true );

          // Get integration points
          StdVector<LocPoint> intPoints;
          StdVector<Double> weights;
          intScheme->GetIntPoints( Elem::GetShapeType(el->type), IntScheme::GAUSS, order, intPoints, weights );
          // Loop over all integration points
          LocPointMapped lpm;
          for( UInt i = 0; i < intPoints.GetSize(); i++  ) {
            //std::cerr << "i = " << i << ", point = " << intPoints[i] << ", weight = " << weights[i] << std::endl;
            // Calculate for each integration point the LocPointMapped
            lpm.Set( intPoints[i], esm, weights[i] );
            this->diffCoef_->GetVector(fac, lpm);

            #pragma omp critical(tempNorm)
            tempNorm += fac.NormL2_squared()*lpm.jacDet * weights[i];
            if( mode_ == "relative" ) {
              // for the relative norm we also need to evaluate the reference result and divide by it afterwards
              this->coef_->GetVector(refFac, lpm);
              #pragma omp critical(refSum)
              refSum += refFac.NormL2_squared()*lpm.jacDet * weights[i];
            }
          } // loop integration points
        } // loop elements
        if( mode_ == "relative" ) {
          outVec[outIt.GetPos()] = sqrt(tempNorm/refSum);
        } else {
          outVec[outIt.GetPos()] = sqrt(tempNorm);
        }
      } else {
          // if the solution for this time/frequency step should not be computed, replace by nan
          outVec[outIt.GetPos()] = NAN;
      }
    } // loop regions
  }

  void PostProcL2Norm::Finalize( ) {
    
  }

  
 // =================== LIMIT-POSTPROCEDURE ===================

  PostProcLimit::PostProcLimit( Grid * ptGrid, 
                              PtrParamNode postProcNode )
    : PostProc( ptGrid, postProcNode ) {
    
    name_ = "limit";
    reducType_ = NONE;

    loLimSet_ = false;
    upLimSet_ = false;    
    loLim_ = 0;
    upLim_ = 0;

  }

  PostProcLimit::~PostProcLimit() {
    
  }

  void PostProcLimit::SetResult( shared_ptr<BaseResult> res ) {
    
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
    newInfo->resultName += "-limit";
    
    // c) Create new Result and set is as output
    shared_ptr<BaseResult> newResult;
    if( res->GetEntryType() == BaseMatrix::DOUBLE ) {
      newResult = shared_ptr<BaseResult>(new Result<Double>() );
    } else {
      newResult = shared_ptr<BaseResult>( new Result<Complex>() );
    }
    newResult->SetEntityList( newList );
    newResult->SetResultInfo( newInfo );
    output_ = newResult;

  }

  void PostProcLimit::SetLowerLimit( Double loLimit ) {

    loLim_ = loLimit;
    loLimSet_ = true;
  }

  void PostProcLimit::SetUpperLimit( Double upLimit ) {

    upLim_ = upLimit;
    upLimSet_ = true;
  }

  void PostProcLimit::Apply() {
    
    if( output_->GetEntryType() == BaseMatrix::DOUBLE ) {
      CalcLimit<Double>();
    } else {
      CalcLimit<Complex>();
    }
  }


  template<class TYPE> 
  void PostProcLimit::CalcLimit() {
    
    // Cast output into correct type
    Result<TYPE> & out = dynamic_cast<Result<TYPE>&>(*output_);
    Result<TYPE> & in =  dynamic_cast<Result<TYPE>&>(*input_ );
    Vector<TYPE> & outVec = out.GetVector();
    Vector<TYPE> & inVec = in.GetVector();

    // Copy input vector to output vector
    outVec = inVec;
    
    // check for lower limit
    if( loLimSet_ ) {
      for( UInt i = 0 ; i< outVec.GetSize(); i++ ) {
        if( opLtSingleDof(outVec[i], loLim_) )
          outVec[i] = loLim_;
      }
    }
    
    // check for upper limit
    if( upLimSet_ ) {
      for( UInt i = 0 ; i< outVec.GetSize(); i++ ) {
        if( opGtSingleDof(outVec[i], upLim_) )
          outVec[i] = upLim_;
      }
    }


  }
  
  template<class TYPE1, class TYPE2>
  bool PostProcLimit::opLtSingleDof( TYPE1 a, TYPE2 b ) {
    return a < b;
  }

  template<>
  bool PostProcLimit::opLtSingleDof( Complex a, Double b ) {
    return std::abs(a) < b;
  }  

  template<class TYPE1, class TYPE2>
  bool PostProcLimit::opGtSingleDof( TYPE1 a, TYPE2 b ) {
    return a > b;
  }

  template<>
  bool PostProcLimit::opGtSingleDof( Complex a, Double b ) {
    return std::abs(a) > b;
  }  
}
