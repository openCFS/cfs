// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "BasePairCoupling.hh"


#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ParamHandling/InfoNode.hh" 
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
                                     ParamNode * paramNode)    
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

    infoNode_ = info->Get("PDEs")->Get(couplingName_, InfoNode::APPEND); 
    InfoNode* in = infoNode_->Get(InfoNode::HEADER); 
    in->Get("sequenceStep")->SetValue(sequenceStep); 
    in->Get("pde1")->SetValue(pde1_->GetName()); 
    in->Get("pde2")->SetValue(pde2_->GetName());

    // get "region" list of current coupling object
    ParamNode * regionListNode = myParam_->Get("regionList", false );
    if(regionListNode) 
    {
      StdVector<ParamNode*> regionList = regionListNode->GetList("region");

      if(regionList.GetSize() > 0) 
      {
        for( UInt i = 0; i < regionList.GetSize(); i++ ) 
        {
          std::string regionName = regionList[i]->Get("name")->AsString(); 
          RegionIdType regionId = ptGrid_->RegionNameToId( regionName ); 

          in->Get("region", InfoNode::APPEND)->Get("name")->SetValue(regionName);

          subdoms_.Push_back( regionId );
          entityLists_.Push_back( ptGrid_->
                                  GetEntityList( EntityList::ELEM_LIST, 
                                                 regionName, EntityList::REGION ) );
                                                          
        }
      }
    }
    
    // get "surfRegion" list of current coupling object
    ParamNode * surfRegionListNode = myParam_->Get("surfRegionList", false );
    if( surfRegionListNode ) 
    {
      StdVector<ParamNode*> surfRegionList = surfRegionListNode->GetList("surfRegion");
      
      if(surfRegionList.GetSize() > 0) 
      {
        InfoNode* list = in->Get("surfaceRegions"); 
        
        for( UInt i = 0; i < surfRegionList.GetSize(); i++ ) 
        {
          std::string surfRegionName = surfRegionList[i]->Get("name")->AsString(); 

          list->Get("surfaceRegion", InfoNode::APPEND)->Get("name")->SetValue(surfRegionName); 

          entityLists_.Push_back( ptGrid_->GetEntityList( EntityList::SURF_ELEM_LIST, 
                                                          surfRegionName, EntityList::REGION ) );
        }
      }
    }

    // ===========================================================
    // Non-conforming grid interfaces
    // ===========================================================
    ParamNode *ncIfListNode = myParam_->Get("ncInterfaceList", false);
    if (ncIfListNode) 
    {
      StdVector<ParamNode*> ncIfList = ncIfListNode->GetList("ncInterface");

      if (ncIfList.GetSize() > 0) 
      {
        InfoNode* list = in->Get("nonConformingGridInterfaces"); 

        for (UInt i = 0; i < ncIfList.GetSize(); ++i) {
          std::string ncIfName = ncIfList[i]->Get("name")->AsString(); 
          RegionIdType ncIfId = ptGrid_->RegionNameToId(ncIfName); 

          InfoNode* e = list->Get("nc_interface", InfoNode::APPEND); 
          e->Get("name")->SetValue(ncIfName); 
          e->Get("id")->SetValue(ncIfId);

          ncIfaces_.Push_back(ncIfId);
        }
      }
    }
    
    // Determine, if axisymmetric geometry is used
    std::string probGeo;
    param->Get("domain")->Get( "geometryType", probGeo );
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
    StdVector<ParamNode*> regionNodes;

    ParamNode * regionListNode = param->
      Get("domain")->Get("regionList", false );
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
        regionNodes[i]->Get( "name", region );
        regionNodes[i]->Get( "material", material );
        regionNodes[i]->Get( "coordSysId", refCoordSys );
        
        // get regionId
        RegionIdType actRegionId = ptGrid_->RegionNameToId( region );
        
        // if no material is set, continue with next loop run
        if( material == "" )
          continue;

        // if region is not contained for current pde, simply continue
        // with next loop
        if( subdoms_.Find( actRegionId) < 0 ) 
          continue;
        
        // print logging information
        Info->PrintF( couplingName_, "Material '%s' for region '%s' (ID = %d) "
                      "follows\n", material.c_str(), region.c_str(),
                      actRegionId );
        // Read data
        materials_[actRegionId] = matLoader->
          LoadMaterial( material, materialClass_ );
        
        // Check for local coordinate system
        if( refCoordSys != "" ) {
          CoordSystem * actCoosy = 
            domain->GetCoordSystem( refCoordSys);
          materials_[actRegionId]->SetCoordSys( actCoosy );
        }
        
        // Check for material rotation parameters
        ParamNode * rotNode = regionNodes[i]->Get("matRotation", false);
        
        Vector<Double> rotVec (3);
        rotVec.Init();
          
        // NOTE: If no rotation is specified and the dimension is
        // 2D, -> material is rotated by
        // alpha = -90 and gamma = -90 degree, 
        // so that we pick by default the yz-plane      
        if( !rotNode ) {
          if( dim_ == 2) {
            rotVec[0] = -90.0;
            rotVec[2] = -90.0;
            materials_[actRegionId]->
              RotateAllTensorsByRotationAngles( rotVec, true );
          }
          continue;
        } else {
          rotNode->Get( "alpha", rotVec[0] );
          rotNode->Get( "beta", rotVec[1] );
          rotNode->Get( "gamma", rotVec[2] );
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
        regionNodes[i]->Get( "name", region );
        regionNodes[i]->Get( "composite", composite );

        // get regionId
        RegionIdType actRegionId = ptGrid_->RegionNameToId( region );
        
        // if no composite is set, continue with next loop run
        if( composite == "" )
          continue;

        // print logging information
        std::ostringstream out;
        out << "Composite material '" << composite << "' for "
            << "region '" << region << "' (ID = " << actRegionId
            << ") follows:\n";
        Info->PrintF( couplingName_, out.str().c_str());
        out.str("");

        // get composite node
        ParamNode * compNode = NULL;
        try {
          compNode = param->Get("domain")
            ->Get("composite", "name", composite );
        } catch( Exception& ex ) {
          RETHROW_EXCEPTION( ex, "No composite material defined with name '"
                             << composite << "'" );
        }

        // get laminaNodes
        StdVector<ParamNode*> laminaNodes = compNode->GetList("lamina");
        
        // Create new lamina and fill ine materials and thicknesses
        Composite & myMat = compositeMaterials_[actRegionId];
        myMat.name = composite;
        // iterate over all single laminas
        for( UInt j = 0; j < laminaNodes.GetSize(); j++ ) {

          // fetch data for lamina
          std::string lamMaterial;
          Double lamThickness, lamOrientation;
          
          laminaNodes[j]->Get( "material", lamMaterial);
          laminaNodes[j]->Get( "thickness", lamThickness);
          laminaNodes[j]->Get( "orientation", lamOrientation);

          // Print information
          out << " --- Lamina " << j+1 << ": thickness = " 
              << lamThickness
              << " m, orientation = " << lamOrientation << "° ---\n";
          Info->PrintF( couplingName_, out.str().c_str());
          out.str("");

          myMat.thickness.Push_back( lamThickness );
          myMat.orientation.Push_back( lamOrientation );
          myMat.materials.Push_back( matLoader->
                                     LoadMaterial( lamMaterial,
                                                   materialClass_ ) );
          
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
    UInt saveBegin, saveEnd, saveInc;
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
    ParamNode * resultNode = myParam_->Get("storeResults", false);
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
      ParamNode* actResultNode = 
        resultNode->Get(xmlElemName, "type", quantity, false );
      
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
      complexFormatString = actResultNode->Get("complexFormat")->AsString();
      String2Enum( complexFormatString, complexFormat );
      
      // otherwise check, if result is to be saved on "allRegions"
      if( actResultNode->Has("allRegions" ) ) {
        ptGrid_->RegionIdToName( regionNames, subdoms_ );
        
        ParamNode * allRegionsNode = actResultNode->Get("allRegions");
        
        std::string allPostProcName, allOutDestName;
        
        // fetch postProcNames
        allRegionsNode->Get("postProcId", allPostProcName );
        postProcNames.Resize( regionNames.GetSize() );
        postProcNames.Init( allPostProcName );
        
        //fetch outDestName
        allRegionsNode->Get("outputIds", allOutDestName );
        outDestNames.Resize( regionNames.GetSize() );
        outDestNames.Init( allOutDestName );

        // fetch saveBegin, saveEnd and saveInc
        allRegionsNode->Get("saveBegin", saveBegin );
        allRegionsNode->Get("saveEnd", saveEnd );
        allRegionsNode->Get("saveInc", saveInc );
        
        // fetch writeResult flag
        std::string writeResult;
        allRegionsNode->Get("writeResult", writeResult );
        writeResults.Resize( regionNames.GetSize() );
        writeResults.Init( writeResult );

      } else {

        StdVector<ParamNode*> regionNodes;
        ParamNode * listNode = NULL;
        // 1b) Look for regions the result is defined on
        if( (*it)->definedOn == ResultInfo::NODE ||
            (*it)->definedOn == ResultInfo::ELEMENT ||
            (*it)->definedOn == ResultInfo::REGION ) {
          listNode = actResultNode->Get("regionList", false);
          if( listNode )
            regionNodes = listNode->GetList("region");
        } else if( (*it)->definedOn == ResultInfo::SURF_ELEM ||
                   (*it)->definedOn == ResultInfo::SURF_REGION ) {
          listNode = actResultNode->Get("surfRegionList", false);
          if( listNode )
            regionNodes = listNode->GetList("surfRegion");

          // fetch entry with neighboring regions
          for( UInt i = 0; i < regionNodes.GetSize(); i++ ) {
            neighborRegions.Push_back( regionNodes[i]->
                                       Get("neighborRegion")->AsString() );
          }
        } 
        
        // only enter, at least one region is present
        if( listNode ) {
          // fetch saveBegin, saveEnd and saveInc
          listNode->Get( "saveBegin", saveBegin );
          listNode->Get( "saveEnd", saveEnd );
          listNode->Get( "saveInc", saveInc );
          
          // iterate over all regions 
          regionNames.Clear();
          postProcNames.Clear();
          outDestNames.Clear();
          writeResults.Clear();
          for( UInt i = 0; i < regionNodes.GetSize(); i++ ) {
            regionNames.Push_back( regionNodes[i]->Get("name")->AsString() );
            postProcNames.Push_back( regionNodes[i]->Get("postProcId")->AsString() );
            outDestNames.Push_back( regionNodes[i]->Get("outputIds")->AsString() );
            writeResults.Push_back( regionNodes[i]->Get("writeResult")->AsString() );
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
          resHandler->RegisterResult( actSol, saveBegin, saveInc, saveEnd,
                                      actOutDest, 
                                      postProcNames[iRegion], writeResult,
                                      isHistory[(*it)->definedOn] );
          
          // if neighboring region is present, store related volume 
          // neighbor region
          if( neighborRegions.GetSize() > 0 ) {
            if( neighborRegions[iRegion] != "" ) {
              surfNeighborRegions_[actSol] = 
                ptGrid_->RegionNameToId( neighborRegions[iRegion] );
            }
          }

        }
        Info->PrintF( couplingName_, "\n");
      }
      
      
      // ========== Look for defineType  node/elemList (history)  ==========

      std::string entityTypeName;
      StdVector<std::string> histNames;
      neighborRegions.Clear();
      
      ParamNode * histNode = NULL;
      StdVector<ParamNode*> histEntities;

      if( (*it)->definedOn == ResultInfo::NODE ) {
        defineType = EntityList::NAMED_NODES;  
        histNode = actResultNode->Get("nodeList",false);
        if( histNode )
          histEntities = histNode->GetList("nodes");
        entityTypeName = "nodes";

      } else if( (*it)->definedOn == ResultInfo::ELEMENT ) {
        defineType = EntityList::NAMED_ELEMS; 
        histNode = actResultNode->Get("elemList",false); 
        if( histNode )
          histEntities = histNode->GetList("elems");
        entityTypeName = "elements";
      
      } else if( (*it)->definedOn == ResultInfo::SURF_ELEM ) {
        defineType = EntityList::NAMED_ELEMS;  
        histNode = actResultNode->Get("surfEelemList",false); 
        if( histNode) 
          histEntities = histNode->GetList("surfElems");
        entityTypeName = "surfElems";

        // fetch entry with neighboring regions
        // fetch entry with neighboring regions
        for( UInt i = 0; i < histEntities.GetSize(); i++ ) {
          neighborRegions.Push_back( histEntities[i]->
                                     Get("neighborRegion")->AsString() );
        }
      }

      // only proceed, if any history result is defined+
      if( histNode ) {
        
        // fetch saveBegin, saveEnd and saveInc
        histNode->Get("saveBegin", saveBegin );
        histNode->Get("saveEnd", saveEnd );
        histNode->Get("saveInc", saveInc );
        
        // iterate over all regions
        histNames.Clear();
        postProcNames.Clear();
        outDestNames.Clear();
        writeResults.Clear();
        for( UInt i = 0; i < histEntities.GetSize(); i++ ) {
          histNames.Push_back( histEntities[i]->Get("name")->AsString() );
          postProcNames.Push_back( histEntities[i]->Get("postProcId")->AsString() );
          outDestNames.Push_back( histEntities[i]->Get("outputIds")->AsString() );
          writeResults.Push_back( histEntities[i]->Get("writeResult")->AsString() );
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
            
          resHandler->RegisterResult( actSol, saveBegin, saveInc, saveEnd,
                                      actOutDest, 
                                      postProcNames[i], writeResult, true );
            
          // if neighboring region is present, store related volume 
          // neighbor region
          if( neighborRegions.GetSize() > 0 ) {
            if( neighborRegions[i] != "" ) {
              surfNeighborRegions_[actSol] = 
                ptGrid_->RegionNameToId( neighborRegions[i] );
            }
          }
            
        }
        Info->PrintF( couplingName_, "\n");
      }
      } catch( Exception &ex ) {
        RETHROW_EXCEPTION(ex, "Could not determine storeResults for quantity '"
                          << quantity << "' within coupling '" 
                          << couplingName_ << "'" );
      }
    }
    Info->PrintF( couplingName_, "\n");
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


  FeFctIdType BasePairCoupling::GetPdeId1() {
    return pde1_->GetPDEId();
  }

  FeFctIdType BasePairCoupling::GetPdeId2() {
    return pde2_->GetPDEId();
  }

} // end of namespace
