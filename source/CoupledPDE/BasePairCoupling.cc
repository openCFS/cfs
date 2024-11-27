#include "BasePairCoupling.hh"

#include "DataInOut/ParamHandling/ParamNode.hh" 
#include "PDE/SinglePDE.hh"
#include "Domain/Domain.hh"
#include "DataInOut/ParamHandling/MaterialHandler.hh"
#include "Domain/Results/ResultFunctor.hh"

#include "Driver/Assemble.hh"
#include "DataInOut/ResultHandler.hh"

// header for saving / retrieving the simulation state
#include "DataInOut/SimState.hh"

// header for Drivers
#include "Driver/TransientDriver.hh"
#include "Driver/TimeSchemes/BaseTimeScheme.hh"

#include <boost/filesystem.hpp>
namespace CoupledField {


  // ***************
  //   Constructor
  // ***************
  BasePairCoupling::BasePairCoupling(SinglePDE *pde1, SinglePDE *pde2,
                                     PtrParamNode paramNode ,
                                     PtrParamNode infoNode,
                                     shared_ptr<SimState> simState,
                                     Domain* domain)    
  {

    
    // initialize pointers
    pde1_           = NULL;
    pde2_           = NULL;
    ptGrid_         = NULL;
    algsys_         = NULL;
    nonLin_         = false;
    simState_       = simState;

    pde1_   = pde1;
    pde2_   = pde2;
    myParam_ = paramNode;
    domain_ = domain;
    ptGrid_ = domain_->GetGrid();
    dim_    = ptGrid_->GetDim(); 

    isaxi_ = false;
    isComplex_ = false;
    
    
    infoNode_ = infoNode;
		considerCounterpart_ = false;
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
    
    infoNode_ = infoNode_->Get(couplingName_, ParamNode::APPEND);
    
    results1_ = pde1_->GetResultInfos();
    results2_ = pde2_->GetResultInfos();
    
    sequenceStep_ = sequenceStep;
    
    // Determine analysistype and use of complex values
    analysisType_ = pde1_->GetAnalysisType();
 
    isComplex_ = domain_->GetDriver()->IsComplex();

    PtrParamNode in = infoNode_->Get(ParamNode::HEADER); 
    in->Get("sequenceStep")->SetValue(sequenceStep); 
    in->Get("pde1")->SetValue(pde1_->GetName()); 
    in->Get("pde2")->SetValue(pde2_->GetName());

    // get "region" list of current coupling object (volumetric coupling)
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
          // Check, if region was already defined an issue a warning otherwise
          if( std::find(subdoms_.Begin(), subdoms_.End(), regionId ) 
          != subdoms_.End() )  {
            WARN( "The region '" << regionName
                  << "' was already defined for coupling '" << couplingName_ 
                  << "'. Please remove duplicate entries." );
          } else {
            subdoms_.Push_back( regionId );
            bool complexMat = false;
            regionList[i]->GetValue("complexMaterial",  complexMat, ParamNode::PASS );
            complexMatData_[regionId] = complexMat;
            entityLists_.Push_back( ptGrid_->
                                    GetEntityList( EntityList::ELEM_LIST, regionName ) ); 
          }
                                                          
        }
      }
    }
    
    // get "surfRegion" list of current coupling object (surface coupling)
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
                                                          surfRegionName ) );
        }
      }
    }

    // ===========================================================
    // Non-conforming grid interfaces
    // ===========================================================
    PtrParamNode ncIfListNode = myParam_->Get("ncInterfaceList", ParamNode::PASS);
    if (ncIfListNode) {
      ParamNodeList ncIfList = ncIfListNode->GetList("ncInterface");
      if (ncIfList.GetSize() > 0) {
        PtrParamNode list = in->Get("nonConformingGridInterfaces"); 
        for (UInt i = 0; i < ncIfList.GetSize(); ++i) {
          std::string ncIfName = ncIfList[i]->Get("name")->As<std::string>();
          RegionIdType ncIfId = ptGrid_->GetRegion().Parse(ncIfName);
          // define non-conforming mortar interface
          StdPDE::NcInterfaceInfo newIface;
          Enum<StdPDE::NcCouplingType> ncCouplingType;
          newIface.interfaceId = ptGrid_->GetNcInterfaceId( ncIfList[i]->Get("name")->As<std::string>() );
          newIface.type = StdPDE::NC_MORTAR;
          newIface.lagrangeMultType = StdPDE::LM_STANDARD;
          ncInterfaceIds_.Push_back(ncIfId);
          ncInterfaces_.Push_back(newIface);
        }
      }
    }
    // Determine, if axisymmetric geometry is used
    isaxi_ = ptGrid_->IsAxi();
    // check if the alg sys has been set
    if ( algsys_ == NULL ) {
      EXCEPTION("BasePairCoupling::Init: The pointer to the algebraic "
                << "system was not set yet! You must call 'SetAlgSys()' "
                << "before!" );
    }
    // Define the FE spaces and their functions
    DefineFeFunctions();
    // Register all fe functions with the algebraic system
    std::map<SolutionType, shared_ptr<BaseFeFunction> >::iterator fncIt= feFunctions_.begin();
    while(fncIt != feFunctions_.end()){
      shared_ptr<BaseFeFunction> actFct = fncIt->second;
      shared_ptr<FeSpace> actSpace = fncIt->second->GetFeSpace();
      std::string fctName = SolutionTypeEnum.ToString(fncIt->first);
      FeFctIdType fctId = algsys_->ObtainFctId( fctName );
      actFct->SetFctId(fctId);
      fncIt++;
    }
    //define primary results (possible new unknowns!)
    DefinePrimaryResults();
    // Define available results
    DefineAvailResults();
    // Initialize nonlinearities
    InitNonLin();
    // Read in material data
    ReadMaterialData();
    // Define the integrators
    DefineIntegrators();
  }
  
  void BasePairCoupling::FinalizeInit() {
    std::map<SolutionType, shared_ptr<BaseFeFunction> >::iterator fncIt= feFunctions_.begin();
    fncIt= feFunctions_.begin();
    // Finalize spaces and fefunctions
    fncIt= feFunctions_.begin();
    while(fncIt != feFunctions_.end()){
      shared_ptr<BaseFeFunction> actFct = fncIt->second;
      shared_ptr<FeSpace> actSpace = fncIt->second->GetFeSpace();
      actSpace->Finalize();
      actSpace->PreCalcShapeFncs();

      // finalize feFunctions
      actFct->Finalize();
      
      // register FeFunctions with SimState class
      simState_->RegisterFeFct( actFct );

      // Pass feFctId of primary result also to RHS result
      rhsFeFunctions_[fncIt->first]->SetFctId(actFct->GetFctId());
      rhsFeFunctions_[fncIt->first]->Finalize();
      fncIt++;
    }

    // Init Time Stepping
    if ( analysisType_ == BasePDE::TRANSIENT ) {
      InitTimeStepping();
    }

    if ( analysisType_ == BasePDE::TRANSIENT ) {
      Double dt;
      dt = dynamic_cast<TransientDriver*>(domain_->GetSingleDriver())
             ->GetDeltaT();
      //WARN("Note: The initialization of the timestepping class is currently wrong: "
      //    "The 2nd argument must be the complete SBM-vector of the algebraic system in "
      //    "order to correclty initialize the internal vectors of the timestepping method. "
      //    "In the old implementation it was sufficient to know the number of unknowns. "
      //    "In the current implementation, the SBM-vectors are just defined within the "
      //    "SolveStep classed. Thus maybe the right thing to do is to shift the creation and "
      //    "initialization of the timestepping scheme to the solveStep classes.")


      // Call the init function of timescheme of each fefunction
      fncIt= feFunctions_.begin();
      while(fncIt != feFunctions_.end()){
        shared_ptr<BaseFeFunction> actFct = fncIt->second;
        actFct->GetTimeScheme()->Init(actFct->GetSingleVector(),dt);
        fncIt++;
      }
    }

    // Define postprocessing results
    DefinePostProcResults();

    // define which solution types have to be saved
    ReadStoreResults();
    ReadSpecialResults();
    ReadSensorArrayResults();
  }
  void BasePairCoupling::ReadSensorArrayResults() {
    // check, if calculation of field variables is requested at all

    ParamNodeList sensorNodes;
    if( !myParam_->Has("storeResults")) 
      return;
    sensorNodes = myParam_->Get("storeResults")->GetList("sensorArray");
    std::string solTypeString;
    static std::map< SolutionType, bool> warningPrinted;

    sensors_.Resize(sensorNodes.GetSize());
    // loop over all parts
    for( UInt iPart = 0; iPart <sensorNodes.GetSize(); ++iPart ) {
      PtrParamNode  actNode = sensorNodes[iPart];
      
      FieldAtPoints & actField = sensors_[iPart];
      actField.fileName = actNode->Get("fileName")->As<std::string>();
      actField.csv = actNode->Get("csv")->As<bool>();
      std::string coordSysId = actNode->Get("coordSysId")->As<std::string>();
      actField.coordSys = domain_->GetCoordSystem(coordSysId);
      
      std::string delim = actNode->Get("delimiter")->As<std::string>();
      if(actField.csv && delim.length() == 0) 
      {
        actField.delim = ',';
      }
      else 
      {        
        actField.delim = delim[0];
      }
      
      // check for solution type
      solTypeString = actNode->Get("type")->As<std::string>();
      
      SolutionType solType = SolutionTypeEnum.Parse(solTypeString);

      // find related result resultinfo
      ResultSet::const_iterator it = availResults_.begin();
      for( ; it != availResults_.end(); ++it ) {
        if( (*it)->resultType == solType ) {
          actField.resultInfo = *it;
          break;
        }
      }
      
      //array of global sensor coordinates
      StdVector< Vector<Double> > globPoints;

      //get entity list of current pde
      StdVector<shared_ptr<EntityList> > lists;
      std::map<RegionIdType, BaseMaterial*> ::iterator regIt = materials_.begin();
      for(; regIt != materials_.end(); regIt++ ) {
        shared_ptr<ElemList> newList(new ElemList(ptGrid_));
        newList->SetRegion(regIt->first);
        lists.Push_back(newList);
      }

      // create list
      // generate new vector
      if(isComplex_) {
        actField.field = new Vector<Complex>();
      } else {
        actField.field = new Vector<Double>();
      }

      if(actNode->Has("parametric")){
        // define sensors according to parametric line definitions
        ParamNodeList listNodes = actNode->Get("parametric")->GetList("list");
        // loop over all components
        StdVector<Double> start(3), stop(3), inc(3);
        StdVector<UInt> numSamples(3);
        start.Init(0);
        stop.Init(0);
        inc.Init(1);
        numSamples.Init(1);
        std::string comp;
        UInt compIndex;
        for( UInt iComp = 0; iComp < listNodes.GetSize(); iComp++ ) {
          PtrParamNode actCompNode = listNodes[iComp];
          actCompNode->GetValue("comp", comp);
          compIndex = actField.coordSys->GetVecComponent(comp)-1;
          start[compIndex]=  actCompNode->Get("start")->MathParse<Double>();
          stop[compIndex]=  actCompNode->Get("stop")->MathParse<Double>();
          inc[compIndex] = actCompNode->Get("inc")->MathParse<Double>();
          numSamples[compIndex] = UInt(floor( (stop[compIndex]-start[compIndex]) / inc[compIndex] ) )+1;
        }

        globPoints.Resize( numSamples[0] *
                           numSamples[1] *
                           numSamples[2] );
        UInt pIdx = 0;

        for( UInt xSample = 0; xSample < numSamples[0]; xSample++ ) {
          Double actX = start[0] + xSample * inc[0];
          for( UInt ySample = 0; ySample < numSamples[1]; ySample++ ) {
            Double actY = start[1] + ySample * inc[1];
            for( UInt zSample = 0; zSample < numSamples[2]; zSample++ ) {
              Double actZ = start[2] + zSample * inc[2];

              // transform global point w.r.t. to coordinate system
              // to global point w.r.t. to global cartesian system
              Vector<Double> globPointcSys;
              globPointcSys.Resize(dim_);

              globPointcSys[0] = actX;
              globPointcSys[1] = actY;
              if( dim_ > 2) {
                globPointcSys[2] = actZ;
              }
              actField.coordSys->Local2GlobalCoord(globPoints[pIdx++],
                                                   globPointcSys);
            } // z
          } // y
        } // x
      }else if(actNode->Has("coordinateFile")){
        globPoints.Reserve(200);

        PtrParamNode coordFileNode = actNode->Get("coordinateFile");
        std::string inFileName = coordFileNode->Get("fileName")->As<std::string>();
        std::string delim = coordFileNode->Get("delimiter")->As<std::string>();
        std::string comment = coordFileNode->Get("commentCharacter")->As<std::string>();
        UInt xCol = coordFileNode->Get("xCoordColumn",ParamNode::PASS)->As<UInt>();
        UInt yCol = coordFileNode->Get("yCoordColumn",ParamNode::PASS)->As<UInt>();
        UInt zCol = coordFileNode->Get("zCoordColumn",ParamNode::PASS)->As<UInt>();

        if(xCol == 0 || yCol ==0 || zCol == 0){
          EXCEPTION("Read coordinate file for sensor array: column indices need to be one based.");
        }

        if(comment.size()>1 || delim.size() > 1){
          WARN("Read coordinate file for sensor array: Comment and delimiter strings need to be single characters!");
        }


        if(!boost::filesystem::exists(inFileName)){
          EXCEPTION("Read coordinate file for sensor array: Could not find coordinate file \"" + inFileName + "\" to read sensor positions!");
          continue;
        }

        std::fstream coordFile(inFileName.c_str(),std::ios::in);

        std::string curLine;
        coordFile >> std::ws;
        UInt lineCounter = 0;
        while(std::getline (coordFile,curLine)){
          lineCounter++;
          //ignore leading whitespace
          std::string::size_type pos = 0;
          while (pos < curLine.size() && std::isspace(curLine[pos], std::locale()))
            pos++;

          curLine.erase(0, pos);

          //check for comment character
          if(curLine.at(0) == comment.at(0)){
            continue;
          }

          //tokenize line with tokenizer
          typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
          boost::char_separator<char> sep(delim.c_str());
          tokenizer tokens(curLine, sep);

          UInt numNumbers = std::distance(tokens.begin(),tokens.end());

          //ignore empty lines
          if(numNumbers==0){
            continue;
          }

          //ignore invalid lines and print a warning
          //here we check for dimension
          if( (dim_ == 2 && numNumbers < 2) ||
              (numNumbers < xCol) ||
              (numNumbers < yCol) ||
              (numNumbers < zCol) ){
            WARN("Read coordinate file for sensor array: Invalid coordinate definition at line: " << lineCounter << " in file : " << inFileName );
          }

          //finally read in the tokens
          Vector<Double> curCoord(dim_);
          tokenizer::iterator tokIter=tokens.begin();
          for(UInt i=0;i<dim_ && tokIter!=tokens.end();i++,tokIter++){
            try{
              curCoord[i] = boost::lexical_cast<Double>(*tokIter);
            }catch(const boost::bad_lexical_cast &e){
              EXCEPTION("Read coordinate file for sensor array: Error reading coordinates in line: " << lineCounter << ". " << *tokIter << " The line was:\n" << curLine << e.what());
            }
          }
          Vector<Double> globPointcSys;
          actField.coordSys->Local2GlobalCoord(globPointcSys,
                                               curCoord);
          globPoints.Push_back(globPointcSys);
        }
        globPoints.Trim();
      }else{
        EXCEPTION("Could not find valid sensor coordinate definition in xml file although tag was given.");
      }
      StdVector< LocPoint > locPoints;
      StdVector< const Elem* > elems;

      // now, map global points to local points restricted to regions of this PDE.
      ptGrid_->GetElemsAtGlobalCoords( globPoints,
                                       locPoints,
                                       elems,
                                       lists,
                                       0.0, 1.0e-2,
                                       false );

      for(UInt i=0, n=globPoints.GetSize(); i<n; i++) {
        const Elem* ptElem = elems[i];
        
        if( !ptElem ) {
          bool wP = !warningPrinted[actField.resultInfo->resultType];
          if( wP ) {
            std::stringstream sstr;
            sstr << "Could not find element at position " 
                 << globPoints[i].ToString()
                 << " for evaluation of field values for "
                 << solTypeString << ".";
            WARN( sstr.str() );
            warningPrinted[actField.resultInfo->resultType] = true;
          }
        } else {
          //               std::cerr << "locPoint for globPoint " << globPoint.ToString() 
          //                                    << " is " << locPoint.ToString() 
          //                                    << " in Elem " << ptElem->elemNum << std::endl;
               
               // check again mapping by performing loc->glob mapping
          // shared_ptr<ElemShapeMap> esm = ptGrid_->GetElemShapeMap(ptElem);
          // shared_ptr<ElemShapeMap> esm(ptElem->ptrShapeMap);
          shared_ptr<ElemShapeMap> esm = ((ptElem))->GetElemShapeMap(ptGrid_);
          //esm->Local2Global(globPoint, locPoint);
          //               std::cerr << "\tAdditional check loc->glob delivers global point " 
          //                   << globPoint.ToString() << std::endl << std::endl;
          
          actField.elems.Push_back(ptElem);
          actField.locPoints.Push_back(locPoints[i]);
        }
      }

      if(warningPrinted[actField.resultInfo->resultType]) {
        std::stringstream sstr;
        sstr << "Could not find " << (globPoints.GetSize()-actField.locPoints.GetSize())
             << " locations for evaluation of field values for "
             << solTypeString << ".";
        WARN( sstr.str() );
      }
    } // loop over <field> entries
  }

  void BasePairCoupling::DefineFeFunctions(){
    //This is the default creation of spaces
    //idee: die PDE gibt zum attribute formulation die passenden space zurueck
    //DOGMA: PRO UNBEKANNTE EINE FUNCTION UND EIN SPACE
    std::string formulation;
    myParam_->GetValue("feSpaceFormulation",formulation,ParamNode::PASS);
    PtrParamNode feSpaceNode = infoNode_->Get("feSpaces");
    std::map<SolutionType, shared_ptr<FeSpace> > spaces;
    CreateFeSpaces(formulation, feSpaceNode, spaces);
    
    //loop over all spaces and set an FeFunction
    std::map<SolutionType, shared_ptr<FeSpace> >::iterator spIt = spaces.begin();
    while(spIt != spaces.end()){
      
      if(feFunctions_.find(spIt->first) != feFunctions_.end()){
        EXCEPTION("It seems that the PDE has created multiple spaces for one result: " << \
                  spIt->first << " This is not how its ought to be!");
      }
      
      if(isComplex_){
        feFunctions_[spIt->first].reset(new FeFunction<Complex>(domain_->GetMathParser()));
        rhsFeFunctions_[spIt->first].reset(new FeFunction<Complex>(domain_->GetMathParser()));
      }
      else{
        feFunctions_[spIt->first].reset(new FeFunction<Double>(domain_->GetMathParser()));
        rhsFeFunctions_[spIt->first].reset(new FeFunction<Double>(domain_->GetMathParser()));
        
        // Note: in the transient case, we also need fefunctions for the time derivatives
        // ... todo: add initialization
      }
      //let the objects know about each other
      spIt->second->AddFeFunction(feFunctions_[spIt->first]);
      feFunctions_[spIt->first]->SetFeSpace(spIt->second);
      //set just PDE1: WORKING?????
      feFunctions_[spIt->first]->SetPDE(pde1_);
      feFunctions_[spIt->first]->SetGrid(ptGrid_);
      
      rhsFeFunctions_[spIt->first]->SetFeSpace(spIt->second);
      //set just PDE1: WORKING?????
      rhsFeFunctions_[spIt->first]->SetPDE(pde1_);
      rhsFeFunctions_[spIt->first]->SetGrid(ptGrid_);
      spIt++;
    }
  }


  void BasePairCoupling::ReadMaterialData() {

    // get list of parameter nodes for region definitions
    ParamNodeList regionNodes;

    PtrParamNode regionListNode = domain_->GetParamRoot()->
      Get("domain")->Get("regionList", ParamNode::PASS );
    if( regionListNode)
      regionNodes = regionListNode->GetList("region");

    // obtain pointer to materialHandler
    MaterialHandler * matLoader = NULL;
    matLoader = domain_->GetMaterialHandler();

    
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
        in->Get("regionList")->GetByVal("region", "name", domain_->GetGrid()->GetRegion().ToString(actRegionId));
        materials_[actRegionId]->ToInfo(in);

        
        // Check for local coordinate system
        if( refCoordSys != "" ) {
          CoordSystem * actCoosy = 
            domain_->GetCoordSystem( refCoordSys);
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
          compNode = domain_->GetParamRoot()->Get("domain")->GetByVal("composite", "name", composite );
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

  void BasePairCoupling::DefineFieldResult( PtrCoefFct coef, shared_ptr<ResultInfo> res ) {

    // create new result functor based on coefficient
    shared_ptr<ResultFunctor> func;
    if( isComplex_ ) {
      func.reset(new FieldCoefFunctor<Complex>(coef, res));
    } else {
      func.reset(new FieldCoefFunctor<Double>(coef, res));
    }
    
    // insert result to list of available results and field functors
    resultFunctors_[res->resultType] = func;
    fieldCoefs_[res->resultType] = coef;
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
    ResultHandler * resHandler = domain_->GetResultHandler();
    
    // initialize map for relating EntityUnknownType and name of xml-element
    std::map<ResultInfo::EntityUnknownType, std::string> elemNames;
    std::map<ResultInfo::EntityUnknownType, bool> isHistory;
    elemNames[ResultInfo::NODE] = "nodeResult";
    elemNames[ResultInfo::ELEMENT] = "elemResult";
    elemNames[ResultInfo::SURF_ELEM] = "surfElemResult";
    elemNames[ResultInfo::REGION] = "regionResult";
    elemNames[ResultInfo::SURF_REGION] = "surfRegionResult";

    isHistory[ResultInfo::NODE] = false;
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
      
      SolutionType actSolType = (*it)->resultType;
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
        entityType = EntityList::NAME_LIST;
      } else if( (*it)->definedOn == ResultInfo::SURF_REGION ) {
        entityType = EntityList::NAME_LIST;
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
        if( listNode != NULL && listNode->HasChildren() ) {
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
          actList = ptGrid_->GetEntityList( entityType, regionNames[iRegion] ); 
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

          // try to get result functor
          shared_ptr<ResultFunctor> fnc;
          if( resultFunctors_.find(actSolType) == 
              resultFunctors_.end() ) {
            EXCEPTION( "No result functor defined for results of type '"
                << quantity << "'");
          }
          fnc = resultFunctors_[actSolType];
          
          // pass result to resulthandler
          resHandler->RegisterResult( actSol, fnc, sequenceStep_,saveBegin, saveInc, 
                                      saveEnd, actOutDest, 
                                      postProcNames[iRegion], writeResult,
                                      isHistory[(*it)->definedOn] );
          
//          // if neighboring region is present, store related volume 
//          // neighbor region
//          if( neighborRegions.GetSize() > 0 ) {
//            if( neighborRegions[iRegion] != "" ) {
//              surfNeighborRegions_[actSol] = 
//                ptGrid_->GetRegion().Parse( neighborRegions[iRegion] );
//            }
//          }

        }
      }
      
      
      // ========== Look for defineType  node/elemList (history)  ==========

      std::string entityTypeName;
      StdVector<std::string> histNames;
      neighborRegions.Clear();
      
      PtrParamNode histNode;
      ParamNodeList histEntities;

      if( (*it)->definedOn == ResultInfo::NODE ) {
        histNode = actResultNode->Get("nodeList",ParamNode::PASS);
        if( histNode )
          histEntities = histNode->GetList("nodes");
        entityTypeName = "nodes";

      } else if( (*it)->definedOn == ResultInfo::ELEMENT ) {
        histNode = actResultNode->Get("elemList",ParamNode::PASS); 
        if( histNode )
          histEntities = histNode->GetList("elems");
        entityTypeName = "elements";
      
      } else if( (*it)->definedOn == ResultInfo::SURF_ELEM ) {
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
          actList = ptGrid_->GetEntityList(entityType, histNames[i] );
          
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
          
          // try to get result functor
          shared_ptr<ResultFunctor> fnc;
          if( resultFunctors_.find(actSolType) == 
              resultFunctors_.end() ) {
            EXCEPTION( "No result functor defined for results of type '"
                << quantity << "'");
          }
          fnc = resultFunctors_[actSolType];
          resHandler->RegisterResult( actSol, fnc, sequenceStep_, 
                                      saveBegin, saveInc, saveEnd,
                                      actOutDest, 
                                      postProcNames[i], writeResult, true );
            
//          // if neighboring region is present, store related volume 
//          // neighbor region
//          if( neighborRegions.GetSize() > 0 ) {
//            if( neighborRegions[i] != "" ) {
//              surfNeighborRegions_[actSol] = 
//                ptGrid_->GetRegion().Parse( neighborRegions[i] );
//            }
//          }
            
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

      // Add calculation of fields, e.g. for the piezo stress
    // ===================================================
       //  Trigger calculation of interpolated field results 
       // ===================================================
       
       // Check for additional field variable
       UInt numFields = sensors_.GetSize();
       
       // loop over all fields variables
       for( UInt i = 0; i < numFields; ++i ) {
       
         // call specialized calculation method in sub-class
         FieldAtPoints& fap = sensors_[i];
         
         
         // Obtain field resultFunctor object
         SolutionType solType = fap.resultInfo->resultType;
         StdVector<std::string> dofNames = fap.resultInfo->dofNames;
         UInt numDofs = dofNames.GetSize();
         std::string solTypeString;
         solTypeString = SolutionTypeEnum.ToString(solType);
         std::map<SolutionType, PtrCoefFct >::iterator fctIt;
         fctIt = fieldCoefs_.find(solType);
         if( fctIt == fieldCoefs_.end() )  {
           EXCEPTION( "Could not find field functor for result '" 
               << SolutionTypeEnum.ToString(solType) << "'");
         }
         PtrCoefFct fct =  fctIt->second;
         
         // calculate vector entries
         if( isComplex_) {
           Vector<Complex> temp;
           Vector<Complex>& vec = dynamic_cast<Vector<Complex> &>(*fap.field);
           vec.Resize(fap.elems.GetSize() * numDofs);
           UInt pos = 0;
           LocPointMapped lpm;
           for ( UInt iElem = 0; iElem < fap.elems.GetSize(); ++iElem ) {
            //  shared_ptr<ElemShapeMap> esm =
            //      ptGrid_->GetElemShapeMap( fap.elems[iElem], true );
            // shared_ptr<ElemShapeMap> esm(fap.elems[iElem]->ptrShapeMap);
            shared_ptr<ElemShapeMap> esm = ((fap.elems[iElem]))->GetElemShapeMap(ptGrid_, true);
             lpm.Set(fap.locPoints[iElem], esm, 0.0);
             fct->GetVector(temp, lpm );
             
             for( UInt i = 0; i < numDofs; ++i ) {
               vec[pos++] = temp[i];
             }
           }
         } else {
           Vector<Double> temp;
           Vector<Double>& vec = dynamic_cast<Vector<Double> &>(*fap.field);
           vec.Resize(fap.elems.GetSize() * numDofs);
           UInt pos = 0;
           LocPointMapped lpm;
           for ( UInt iElem = 0; iElem < fap.elems.GetSize(); ++iElem ) {
            //  shared_ptr<ElemShapeMap> esm =
            //      ptGrid_->GetElemShapeMap( fap.elems[iElem], true );
            // shared_ptr<ElemShapeMap> esm(fap.elems[iElem]->ptrShapeMap);
            shared_ptr<ElemShapeMap> esm = ((fap.elems[iElem]))->GetElemShapeMap(ptGrid_, true);
             lpm.Set(fap.locPoints[iElem], esm, 0.0);
             fct->GetVector(temp, lpm );

             for( UInt i = 0; i < numDofs; ++i ) {
               vec[pos++] = temp[i];
             }
           }
         }


         std::ofstream  out((fap.fileName+"-"+lexical_cast<std::string>(kstep)).c_str(),
                             std::ios::out );
         // Ensure that no precision is lost
         out.precision(15);
               
         Vector<Double> globPoint, globPointcSys;
         
         StdVector<std::string> globCoordNames;
         StdVector<std::string> locCoordNames;
         for(UInt i = 0; i < dim_; ++i ) {
           globCoordNames.Push_back(fap.coordSys->GetDofName(i+1));
         }
         locCoordNames.Push_back("xi");
         locCoordNames.Push_back("eta");
         locCoordNames.Push_back("zeta");      
         std::string delim = "\t";
         if(fap.csv) 
         {
           delim = std::string(1,fap.delim);
         }
         
         
         // print out information
         if( isComplex_ ){
           // cast solution vector
           Vector<Complex>& vec = dynamic_cast<Vector<Complex> &>(*(fap.field));

           // Write header line with descriptions of columns
           if(fap.csv) 
           {
             out << "origElemNum" << delim;        
             for(UInt j = 0; j < dim_; ++j ) {
               out << "globCoord-" << globCoordNames[j] << delim;
             }
             for(UInt j = 0; j < numDofs; ++j ) {
               out << solTypeString << "-real" << "-" << dofNames[j] << delim;
             }
             for(UInt j = 0; j < numDofs; ++j ) {
               out << solTypeString << "-imag" << "-" << dofNames[j] << delim;
             }
             for(UInt j = 0; j < dim_-1; ++j ) {
               out << "locCoord-" << locCoordNames[j] << delim;
             }
             out << "locCoord-" << locCoordNames[dim_-1] << std::endl;
           }
         
           // Loop over all points
           for( UInt iPoint = 0; iPoint < fap.locPoints.GetSize(); iPoint++) { 
             
            //  shared_ptr<ElemShapeMap> esm =
            //      ptGrid_->GetElemShapeMap(fap.elems[iPoint], true);
            // shared_ptr<ElemShapeMap> esm(fap.elems[iPoint]->ptrShapeMap);
            shared_ptr<ElemShapeMap> esm = ((fap.elems[iPoint]))->GetElemShapeMap(ptGrid_, true);
             esm->Local2Global(globPoint, fap.locPoints[iPoint]);
             
             fap.coordSys->Global2LocalCoord(globPointcSys, globPoint);
             
             // write to file
             out << fap.elems[iPoint]->elemNum << delim;
             out << globPointcSys.ToString(TS_PLAIN, delim) << delim;
             for(UInt j = 0; j < numDofs; ++j ) {
               out << vec[iPoint*numDofs + j].real() << delim;
             }
             for(UInt j = 0; j < numDofs; ++j ) {
               out << vec[iPoint*numDofs + j].imag() << delim;
             }
             for(UInt j = 0; j < dim_-1; ++j ) {
               out << fap.locPoints[iPoint][j] << delim;
             }
             out << fap.locPoints[iPoint][dim_-1] << std::endl;
           }

         } else {
           // cast solution vector
           Vector<Double>& vec = dynamic_cast<Vector<Double> &>(*(fap.field));

           // Write header line with descriptions of columns
           if(fap.csv) 
           {
             out << "origElemNum" << delim;        
             for(UInt j = 0; j < dim_; ++j ) {
               out << "globCoord-" << globCoordNames[j] << delim;
             }
             for(UInt j = 0; j < numDofs; ++j ) {
               out << solTypeString << "-" << dofNames[j] << delim;
             }
             for(UInt j = 0; j < dim_-1; ++j ) {
               out << "locCoord-" << locCoordNames[j] << delim;
             }
             out << "locCoord-" << locCoordNames[dim_-1] << std::endl;
           }
           
           // Loop over all points
           for( UInt iPoint = 0; iPoint < fap.locPoints.GetSize(); iPoint++) { 
             
            //  shared_ptr<ElemShapeMap> esm =
            //      ptGrid_->GetElemShapeMap(fap.elems[iPoint], true);
            // shared_ptr<ElemShapeMap> esm(fap.elems[iPoint]->ptrShapeMap);
            shared_ptr<ElemShapeMap> esm = ((fap.elems[iPoint]))->GetElemShapeMap(ptGrid_, true);
             esm->Local2Global(globPoint, fap.locPoints[iPoint]);
             
             fap.coordSys->Global2LocalCoord(globPointcSys, globPoint);
             // write to file
             out << fap.elems[iPoint]->elemNum << delim;
             out << globPointcSys.ToString(TS_PLAIN, delim) << delim;
             for(UInt j = 0; j < numDofs; ++j ) {
               out << vec[iPoint*numDofs + j] << delim;
             }
             for(UInt j = 0; j < dim_-1; ++j ) {
               out << fap.locPoints[iPoint][j] << delim;
             }
             out << fap.locPoints[iPoint][dim_-1] << std::endl;
           }
         }

         out.close();
         
       }
  }

} // end of namespace
