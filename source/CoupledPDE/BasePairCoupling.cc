// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "BasePairCoupling.hh"

#include "DataInOut/ParamHandling/ParamNode.hh" 
#include "PDE/SinglePDE.hh"
#include "Domain/domain.hh"
#include "DataInOut/MaterialHandler.hh"
#include "Materials/piezoMaterial.hh"
#include "Driver/assemble.hh"
#include "DataInOut/resultHandler.hh"

namespace CoupledField {


  // ***************
  //   Constructor
  // ***************
  BasePairCoupling::BasePairCoupling(SinglePDE *pde1, SinglePDE *pde2,
                                     PtrParamNode paramNode)    
  {

    
    // initialize pointers
    sol_            = NULL;
    solVec_         = NULL;
    pde1_           = NULL;
    pde2_           = NULL;
    ptGrid_         = NULL;
    algsys_         = NULL;
    nonLin_         = false;
    nonLinMaterial_ = false;
    isHysteresis_   = false;
    nonLinHysteresis_ = false;

    pde1_   = pde1;
    pde2_   = pde2;
    myParam_ = paramNode;
    ptGrid_ = domain->GetGrid();

    isaxi_ = false;
    isComplex_ = false;
    
    dim_ = domain->GetGrid()->GetDim();
    infoNode_ = info->Get("PDE")->Get(couplingName_, ParamNode::APPEND);
  }


  // **************
  //   Destructor
  // **************
  BasePairCoupling::~BasePairCoupling() {


    pde1_ = NULL;
    pde2_ = NULL;

    std::map<RegionIdType, BaseMaterial*>::iterator it;
    for ( it = materials_.begin(); it != materials_.end(); it++ ) {
      delete it->second;
    }
    materials_.clear();

  }


  // ********
  //   Init
  // ********
  void BasePairCoupling::Init( UInt sequenceStep ) {
    
    
    results1_ = pde1_->GetResultInfos();
    results2_ = pde2_->GetResultInfos();
    
    sequenceStep_ = sequenceStep;
    
    // Determine analysistype and use of complex values
    analysisType_ = pde1_->GetAnalysisType();
 
    isComplex_ = BasePDE::IsComplex(analysisType_);

    PtrParamNode in = infoNode_->Get(ParamNode::HEADER); 
    in->Get("sequenceStep")->SetValue(sequenceStep); 
    in->Get("pde1")->SetValue(pde1_->GetName()); 
    in->Get("pde2")->SetValue(pde2_->GetName());

    // get "region" list of current coupling object
    PtrParamNode regionListNode = myParam_->Get("regionList", ParamNode::PASS );
    if(regionListNode) 
    {
      ParamNodeList regionList = regionListNode->GetList("region");

      if(regionList.GetSize() > 0) 
      {
        for( UInt i = 0; i < regionList.GetSize(); i++ ) 
        {
          std::string regionName = regionList[i]->Get("name")->As<std::string>(); 
          RegionIdType regionId = ptGrid_->GetRegion().Parse( regionName );

          in->Get("region", ParamNode::APPEND)->Get("name")->SetValue(regionName);
          
          subdoms_.Push_back( regionId );
          bool complexMat = false;
          regionList[i]->GetValue("complexMaterial",  complexMat, ParamNode::PASS );
          complexMatData_[regionId] = complexMat;
          entityLists_.Push_back( ptGrid_->
                                  GetEntityList( EntityList::ELEM_LIST, 
                                                 regionName, EntityList::REGION ) );
                                                          
        }
      }
    }
    
    // get "surfRegion" list of current coupling object
    PtrParamNode surfRegionListNode = myParam_->Get("surfRegionList", ParamNode::PASS );
    if( surfRegionListNode ) 
    {
      ParamNodeList surfRegionList = surfRegionListNode->GetList("surfRegion");
      
      if(surfRegionList.GetSize() > 0) 
      {
        PtrParamNode list = in->Get("surfaceRegions"); 
        
        for( UInt i = 0; i < surfRegionList.GetSize(); i++ ) 
        {
          std::string surfRegionName = surfRegionList[i]->Get("name")->As<std::string>(); 

          list->Get("surfaceRegion", ParamNode::APPEND)->Get("name")->SetValue(surfRegionName); 

          entityLists_.Push_back( ptGrid_->GetEntityList( EntityList::SURF_ELEM_LIST, 
                                                          surfRegionName, EntityList::REGION ) );
        }
      }
    }

    // ===========================================================
    // Non-conforming grid interfaces
    // ===========================================================
    PtrParamNode ncIfListNode = myParam_->Get("ncInterfaceList", ParamNode::PASS);
    if (ncIfListNode) 
    {
      ParamNodeList ncIfList = ncIfListNode->GetList("ncInterface");

      if (ncIfList.GetSize() > 0) 
      {
        PtrParamNode list = in->Get("nonConformingGridInterfaces"); 

        for (UInt i = 0; i < ncIfList.GetSize(); ++i) {
          std::string ncIfName = ncIfList[i]->Get("name")->As<std::string>(); 
          RegionIdType ncIfId = ptGrid_->GetRegion().Parse(ncIfName);

          PtrParamNode e = list->Get("nc_interface", ParamNode::APPEND); 
          e->Get("name")->SetValue(ncIfName); 
          e->Get("id")->SetValue(ncIfId);

          ncIfaces_.Push_back(ncIfId);
        }
      }
    }
    
    // Determine, if axisymmetric geometry is used
    std::string probGeo;
    param->Get("domain")->GetValue( "geometryType", probGeo );
    if( probGeo == "axi" ) {
      isaxi_ = true;
    } else {
      isaxi_ = false;
    }

    // Get type of analysis and create according 
    // assemble object
    // -> copy simply from first pde
    std::string help = BasePDE::analysisType.ToString((*pde1_).analysistype_);
    //std::cerr << "Analysis of PDE is " 
    //<< help << std::endl;

    if ( algsys_ == NULL ) {
      EXCEPTION("BasePairCoupling::Init: The pointer to the algebraic "
                << "system was not set yet! You must call 'SetAlgSys()' "
                << "before!" );
    }


    eqnMap1_ = pde1_->GetEqnMap();
    eqnMap2_ = pde2_->GetEqnMap();
    assert( eqnMap1_ != NULL);
    assert( eqnMap2_ != NULL);


    // Define available results
    DefineAvailResults();

    // Initialize nonlinearities
    InitNonLin();

    // Read in material data
    ReadMaterialData();
 
    // Define the integrators
    DefineIntegrators();

    
    // define which solution types have to be saved
    ReadStoreResults();
    ReadSpecialResults();

  }

  void BasePairCoupling::ReadMaterialData() {

    // get list of parameter nodes for region definitions
    ParamNodeList regionNodes;

    PtrParamNode regionListNode = param->
      Get("domain")->Get("regionList", ParamNode::PASS );
    if( regionListNode)
      regionNodes = regionListNode->GetList("region");

    // obtain pointer to materialHandler
    MaterialHandler * matLoader = NULL;
    matLoader = domain->GetMaterialHandler();

    
    // -------------------
    // NORMAL MATERIALS
    // -------------------
    std::string region, material, composite, refCoordSys;

    // iterate over all regions
    for( UInt i = 0; i < regionNodes.GetSize(); i++ ) {
      try{ 
        // get data from node
        regionNodes[i]->GetValue( "name", region );
        regionNodes[i]->GetValue( "material", material );
        regionNodes[i]->GetValue( "coordSysId", refCoordSys );
        
        // get regionId
        RegionIdType actRegionId = ptGrid_->GetRegion().Parse( region );
        
        // if no material is set, continue with next loop run
        if( material == "" )
          continue;

        // if region is not contained for current pde, simply continue
        // with next loop
        if( subdoms_.Find( actRegionId) < 0 ) 
          continue;
        
        // Read data
        materials_[actRegionId] = matLoader->LoadMaterial( material, materialClass_ );

        // log the just read material. LoadMaterial() so to say initializes the ToInfo()
        PtrParamNode in = infoNode_->GetByVal("material", "name", material);
        // additional regions are automatically appended
        in->Get("regionList")->GetByVal("region", "name", domain->GetGrid()->GetRegion().ToString(actRegionId));
        materials_[actRegionId]->ToInfo(in);

        
        // Check for local coordinate system
        if( refCoordSys != "" ) {
          CoordSystem * actCoosy = 
            domain->GetCoordSystem( refCoordSys);
          materials_[actRegionId]->SetCoordSys( actCoosy );
        }
        
        // Check for material rotation parameters
        PtrParamNode rotNode = regionNodes[i]->Get("matRotation", ParamNode::INSERT);
        
        Vector<Double> rotVec (3);
        rotVec.Init();
          
        // NOTE: If no rotation is specified and the dimension is
        // 2D, -> material is rotated by
        // alpha = -90 and gamma = -90 degree, 
        // so that we pick by default the yz-plane      
        if( !rotNode->HasChildren() ) {
          if( dim_ == 2) {
            rotVec[0] = -90.0;
            rotVec[2] = -90.0;
            materials_[actRegionId]->
              RotateAllTensorsByRotationAngles( rotVec, true );
          }
          continue;
        } else {
          rotVec[0] = rotNode->Get( "alpha" )->MathParse<Double>();
          rotVec[1] = rotNode->Get( "beta" )->MathParse<Double>(); 
          rotVec[2] = rotNode->Get( "gamma" )->MathParse<Double>();
          materials_[actRegionId]->
            RotateAllTensorsByRotationAngles( rotVec, true );
        }

        
      } catch (Exception& ex ) {
        std::string matClassString;
        Enum2String( materialClass_, matClassString );
        RETHROW_EXCEPTION(ex, "Could not assign material '"
                          << material << "' of materialClass '"
                          << matClassString << "' to region '" 
                          << region << "' within coupling '"
                          << couplingName_ << "'");
      }
    }
      
      
    
    // -------------------
    // COMPOSITE MATERIALS
    // -------------------

    // iterate over all regions
    for( UInt i = 0; i < regionNodes.GetSize(); i++ ) {

      try{ 
        // get data from node
        regionNodes[i]->GetValue( "name", region );
        regionNodes[i]->GetValue( "composite", composite );

        // get regionId
        RegionIdType actRegionId = ptGrid_->GetRegion().Parse( region );
        
        // if no composite is set, continue with next loop run
        if( composite == "" )
          continue;
        
        if( subdoms_.Find( actRegionId) < 0 )
          continue;

        // print logging information
        PtrParamNode in = infoNode_->Get(ParamNode::HEADER)->GetByVal("region", "name", region)->Get("composite");

        // get composite node
        PtrParamNode compNode;
        try {
          compNode = param->Get("domain")->GetByVal("composite", "name", composite );
        } catch( Exception& ex ) {
          RETHROW_EXCEPTION( ex, "No composite material defined with name '" << composite << "'" );
        }

        // get laminaNodes
        ParamNodeList laminaNodes = compNode->GetList("lamina");
        
        // Create new lamina and fill ine materials and thicknesses
        Composite & myMat = compositeMaterials_[actRegionId];
        myMat.name = composite;
        // iterate over all single laminas
        for( UInt j = 0; j < laminaNodes.GetSize(); j++ ) {

          // fetch data for lamina
          std::string lamMaterial;
          Double lamThickness, lamOrientation;
          
          laminaNodes[j]->GetValue( "material", lamMaterial);
          laminaNodes[j]->GetValue( "thickness", lamThickness);
          laminaNodes[j]->GetValue( "orientation", lamOrientation);

          // Print information
          PtrParamNode lam = in->Get("lamina", ParamNode::APPEND);
          lam->Get("nr")->SetValue(j+1);
          lam->Get("thickness")->SetValue(lamThickness);
          lam->Get("orientation")->SetValue(lamOrientation);

          myMat.thickness.Push_back( lamThickness );
          myMat.orientation.Push_back( lamOrientation );
          myMat.materials.Push_back(matLoader->LoadMaterial(lamMaterial, materialClass_ ));
          
          PtrParamNode lm_ = lam->Get("material");
          lm_->Get("name")->SetValue(lamMaterial);
          myMat.materials.Last()->ToInfo(lm_);

        } // over single laminae
      } catch (Exception& ex ) {
        RETHROW_EXCEPTION( ex, "Could not create composite material '"
                           << composite << "'");
      }
    } // over composite 
  }

  void BasePairCoupling::ReadStoreResults() {

 StdVector<std::string> regionNames, nodeNames, writeResults, actOutDest;
    StdVector<std::string> postProcNames, outDestNames, neighborRegions;
    UInt saveBegin = 0, saveEnd = 0, saveInc = 0;
    std::string quantity, complexFormatString, listElemName, entityName;
    ComplexFormat complexFormat;
    shared_ptr<EntityList> actList;

    ResultSet::iterator it;
    // initializing entityType to NODE_LIST or receive compiler warning 
    EntityList::ListType entityType = EntityList::NODE_LIST;
        EntityList::DefineType defineType;
    ResultHandler * resHandler = domain->GetResultHandler();
    
    // initialize map for relating EntityUnknownType and name of xml-element
    std::map<ResultInfo::EntityUnknownType, std::string> elemNames;
    std::map<ResultInfo::EntityUnknownType, bool> isHistory;
    elemNames[ResultInfo::NODE] = "nodeResult";
    elemNames[ResultInfo::ELEMENT] = "elemResult";
    elemNames[ResultInfo::SURF_ELEM] = "surfElemResult";
    elemNames[ResultInfo::REGION] = "regionResult";
    elemNames[ResultInfo::SURF_REGION] = "surfRegionResult";

    isHistory[ResultInfo::NODE] = false;
    isHistory[ResultInfo::PFEM] = false;
    isHistory[ResultInfo::ELEMENT] = false;
    isHistory[ResultInfo::SURF_ELEM] = false;
    isHistory[ResultInfo::REGION] = true;
    isHistory[ResultInfo::SURF_REGION] = true;
    
    // fetch result node and leave, if none is present
    PtrParamNode resultNode = myParam_->Get("storeResults", ParamNode::PASS);
    if( !resultNode )
      return;

    // Iterate over all availabe results
    for (it = availResults_.begin(); it != availResults_.end(); it++ ) {
      
      
      // Convert enum
      quantity = SolutionTypeEnum.ToString((*it)->resultType);
      
      // try to catch possible errors 
      try {
        
      // Get type of result
      std::string xmlElemName = elemNames[(*it)->definedOn];
      if( xmlElemName == "" ){
        break;
      }
      
      // Remeber current result node
      PtrParamNode actResultNode = 
        resultNode->GetByVal(xmlElemName, "type", quantity, ParamNode::PASS );
      
      // Check on which entity type the result is defined on 
      if( (*it)->definedOn == ResultInfo::NODE ) {
        entityType = EntityList::NODE_LIST;
      } else if( (*it)->definedOn == ResultInfo::REGION ) {
        entityType = EntityList::REGION_LIST;
      } else if( (*it)->definedOn == ResultInfo::SURF_REGION ) {
        entityType = EntityList::REGION_LIST;
      } else if( (*it)->definedOn == ResultInfo::SURF_ELEM ) {
        entityType = EntityList::SURF_ELEM_LIST;
      } else if( (*it)->definedOn == ResultInfo::ELEMENT ) {
        entityType = EntityList::ELEM_LIST;
      }
      
      // intialize variables
      neighborRegions.Clear();
      regionNames.Clear();
      //outDestNames.Clear();

      // ========== Look for defineType 'REGION' ==========
      // 1a) Look if result is defined on 'allRegions'
      defineType = EntityList::REGION;
      
      

      // if no node was found, continue with next result
      if( !actResultNode) {
        continue;
      }

      // determine complexFormat
      complexFormatString = "amplPhase";
      actResultNode->GetValue("complexFormat", complexFormatString, ParamNode::PASS);
      String2Enum( complexFormatString, complexFormat );
      
      // otherwise check, if result is to be saved on "allRegions"
      if( actResultNode->Has("allRegions" ) ) {
        ptGrid_->GetRegion().ToString( subdoms_, regionNames );
        
        PtrParamNode allRegionsNode = actResultNode->Get("allRegions");
        
        std::string allPostProcName, allOutDestName;
        
        // fetch postProcNames
        allRegionsNode->GetValue("postProcId", allPostProcName );
        postProcNames.Resize( regionNames.GetSize() );
        postProcNames.Init( allPostProcName );
        
        //fetch outDestName
        allRegionsNode->GetValue("outputIds", allOutDestName );
        outDestNames.Resize( regionNames.GetSize() );
        outDestNames.Init( allOutDestName );

        // fetch saveBegin, saveEnd and saveInc
        saveBegin = allRegionsNode->Get("saveBegin")->MathParse<UInt>();
        saveEnd = allRegionsNode->Get("saveEnd")->MathParse<UInt>();
        saveInc = allRegionsNode->Get("saveInc")->MathParse<UInt>();
        
        // fetch writeResult flag
        std::string writeResult;
        allRegionsNode->GetValue("writeResult", writeResult );
        writeResults.Resize( regionNames.GetSize() );
        writeResults.Init( writeResult );

      } else {

        ParamNodeList regionNodes;
        PtrParamNode listNode;
        // 1b) Look for regions the result is defined on
        if( (*it)->definedOn == ResultInfo::NODE ||
            (*it)->definedOn == ResultInfo::ELEMENT ||
            (*it)->definedOn == ResultInfo::REGION ) {
          listNode = actResultNode->Get("regionList", ParamNode::PASS);
          if( listNode )
            regionNodes = listNode->GetList("region");
        } else if( (*it)->definedOn == ResultInfo::SURF_ELEM ||
                   (*it)->definedOn == ResultInfo::SURF_REGION ) {
          listNode = actResultNode->Get("surfRegionList", ParamNode::PASS);
          if( listNode )
            regionNodes = listNode->GetList("surfRegion");

          // fetch entry with neighboring regions
          for( UInt i = 0; i < regionNodes.GetSize(); i++ ) {
            neighborRegions.Push_back( regionNodes[i]->
                                       Get("neighborRegion")->As<std::string>() );
          }
        } 
        
        // only enter, at least one region is present
        if( listNode->HasChildren() ) {
          // fetch saveBegin, saveEnd and saveInc
          saveBegin = listNode->Get("saveBegin")->MathParse<UInt>();
          saveEnd = listNode->Get("saveEnd")->MathParse<UInt>();
          saveInc = listNode->Get("saveInc")->MathParse<UInt>();
          
          // iterate over all regions 
          regionNames.Clear();
          postProcNames.Clear();
          outDestNames.Clear();
          writeResults.Clear();
          for( UInt i = 0; i < regionNodes.GetSize(); i++ ) {
            regionNames.Push_back( regionNodes[i]->Get("name")->As<std::string>() );
            postProcNames.Push_back( regionNodes[i]->Get("postProcId")->As<std::string>() );
            outDestNames.Push_back( regionNodes[i]->Get("outputIds")->As<std::string>() );
            writeResults.Push_back( regionNodes[i]->Get("writeResult")->As<std::string>() );
          }
        }
      }

      // Check, if any region was found for this result type
      if( regionNames.GetSize() != 0 ) 
      {
        (*it)->complexFormat = complexFormat;
        
        // iterate over all regions
        for( UInt iRegion = 0; iRegion < regionNames.GetSize(); iRegion++ ) 
        {
          actList = ptGrid_->GetEntityList( entityType, regionNames[iRegion], 
                                            defineType );
          shared_ptr<BaseResult> actSol;
          if( isComplex_ ) {
            actSol = shared_ptr<BaseResult>(new Result<Complex>());
          } else {
            actSol = shared_ptr<BaseResult>(new Result<Double>());
          }

          // intialize result object
          actSol->SetResultInfo( *it );
          actSol->SetEntityList( actList );
          resultLists_[*it].Push_back( actSol );

          // extract all output destinations and determine bool flag for writeResult
          SplitStringList( outDestNames[iRegion], actOutDest, ',' );
          bool writeResult = writeResults[iRegion] == "yes"  ? true : false ;

          // pass result to resulthandler
          resHandler->RegisterResult( actSol, sequenceStep_,saveBegin, saveInc, 
                                      saveEnd, actOutDest, 
                                      postProcNames[iRegion], writeResult,
                                      isHistory[(*it)->definedOn] );
          
          // if neighboring region is present, store related volume 
          // neighbor region
          if( neighborRegions.GetSize() > 0 ) {
            if( neighborRegions[iRegion] != "" ) {
              surfNeighborRegions_[actSol] = 
                ptGrid_->GetRegion().Parse( neighborRegions[iRegion] );
            }
          }

        }
      }
      
      
      // ========== Look for defineType  node/elemList (history)  ==========

      std::string entityTypeName;
      StdVector<std::string> histNames;
      neighborRegions.Clear();
      
      PtrParamNode histNode;
      ParamNodeList histEntities;

      if( (*it)->definedOn == ResultInfo::NODE ) {
        defineType = EntityList::NAMED_NODES;  
        histNode = actResultNode->Get("nodeList",ParamNode::PASS);
        if( histNode )
          histEntities = histNode->GetList("nodes");
        entityTypeName = "nodes";

      } else if( (*it)->definedOn == ResultInfo::ELEMENT ) {
        defineType = EntityList::NAMED_ELEMS; 
        histNode = actResultNode->Get("elemList",ParamNode::PASS); 
        if( histNode )
          histEntities = histNode->GetList("elems");
        entityTypeName = "elements";
      
      } else if( (*it)->definedOn == ResultInfo::SURF_ELEM ) {
        defineType = EntityList::NAMED_ELEMS;  
        histNode = actResultNode->Get("surfEelemList",ParamNode::PASS); 
        if( histNode) 
          histEntities = histNode->GetList("surfElems");
        entityTypeName = "surfElems";

        // fetch entry with neighboring regions
        // fetch entry with neighboring regions
        for( UInt i = 0; i < histEntities.GetSize(); i++ ) {
          neighborRegions.Push_back( histEntities[i]->
                                     Get("neighborRegion")->As<std::string>() );
        }
      }

      // only proceed, if any history result is defined+
      if( histNode && histNode->HasChildren() ) {
        
        // fetch saveBegin, saveEnd and saveInc
        saveBegin = histNode->Get("saveBegin")->MathParse<UInt>();
        saveEnd = histNode->Get("saveEnd")->MathParse<UInt>();
        saveInc = histNode->Get("saveInc")->MathParse<UInt>();

        // iterate over all regions
        histNames.Clear();
        postProcNames.Clear();
        outDestNames.Clear();
        writeResults.Clear();
        for( UInt i = 0; i < histEntities.GetSize(); i++ ) {
          histNames.Push_back( histEntities[i]->Get("name")->As<std::string>() );
          postProcNames.Push_back( histEntities[i]->Get("postProcId")->As<std::string>() );
          outDestNames.Push_back( histEntities[i]->Get("outputIds")->As<std::string>() );
          writeResults.Push_back( histEntities[i]->Get("writeResult")->As<std::string>() );
        }
      }
      
      if( histNames.GetSize() > 0 ) {
        
        (*it)->complexFormat = complexFormat;

        // iterate over all entityNames
        for( UInt i = 0; i < histNames.GetSize(); i++ ) 
        {
          actList = ptGrid_->GetEntityList(entityType, histNames[i], defineType);
          
          shared_ptr<BaseResult> actSol;
          if( isComplex_ ) {
            actSol = shared_ptr<BaseResult>(new Result<Complex>());
          } else {
            actSol = shared_ptr<BaseResult>(new Result<Double>());
          }
            
          // Set result info and entitylist at the result object
          actSol->SetResultInfo( *it );
          actSol->SetEntityList( actList );
          resultLists_[*it].Push_back( actSol );
            
          // extract all output destinations and determine bool flag for writeResult
          SplitStringList( outDestNames[i], actOutDest, ',' );
          bool writeResult = (writeResults[i] == "yes"  ? true : false );
            
          resHandler->RegisterResult( actSol, sequenceStep_, 
                                      saveBegin, saveInc, saveEnd,
                                      actOutDest, 
                                      postProcNames[i], writeResult, true );
            
          // if neighboring region is present, store related volume 
          // neighbor region
          if( neighborRegions.GetSize() > 0 ) {
            if( neighborRegions[i] != "" ) {
              surfNeighborRegions_[actSol] = 
                ptGrid_->GetRegion().Parse( neighborRegions[i] );
            }
          }
            
        }
      }
      } catch( Exception &ex ) {
        RETHROW_EXCEPTION(ex, "Could not determine storeResults for quantity '"
                          << quantity << "' within coupling '" 
                          << couplingName_ << "'" );
      }
    }
  }

  void BasePairCoupling::WriteResultsInFile( const UInt kstep,
                                             const Double asteptime ) {

    ResultMap::iterator it = resultLists_.begin();
    ResultHandler * resHandler = domain->GetResultHandler();
    
    // iterate over all results
    for( ; it != resultLists_.end(); it++ ) {
      ResultList & actList = it->second;
      
      // iterate over all solutions for each result type
      for( UInt i = 0; i < actList.GetSize(); i++ ) {

        // get string representation of quantity and entity list
        std::string listName;
        std::string quantity = SolutionTypeEnum.ToString(actList[i]->GetResultInfo()->resultType);
        listName = actList[i]->GetEntityList()->GetName();

        // Only calculate result, if needed
        if( resHandler->IsResultNeeded( actList[i] ) ) {
          try {
            CalcResults( actList[i] );
          } catch (Exception &ex ) {
            RETHROW_EXCEPTION( ex, "Could not calculate result '" << quantity
                               << "' on '" << listName << "' in pde '" 
                               << couplingName_ << "'");
          }
          try {
            resHandler->UpdateResult( actList[i] );
          } catch (Exception &ex ) {
            RETHROW_EXCEPTION( ex, "Could not write result '" << quantity
                               << "' on '" << listName 
                               << "' to output file(s) in coupling '" 
                               << couplingName_ << "'");
          }  
        }
      }
    }
    
  }


  PdeIdType BasePairCoupling::GetPdeId1() {
    return pde1_->GetPDEId();
  }

  PdeIdType BasePairCoupling::GetPdeId2() {
    return pde2_->GetPDEId();
  }

} // end of namespace
