
#include "IPreciceAdapter.hh"
#include "PreciceAdapter.hh"
#include "MinimalXmlParser.hh"
#include "NodeResult.hh"
#include "ElementResult.hh"

#include "DataInOut/Logging/LogConfigurator.hh"

#include "Domain/Domain.hh"
#include "Domain/Mesh/GridCFS/GridCFS.hh"
#include "TransientDriverPrecice.hh"
#include "Driver/SolveSteps/StdSolveStep.hh"
#include "PDE/SinglePDE.hh"
#include "Utils/StdVector.hh"
#include "MatVec/Vector.hh"
#include <numeric>
#include <unordered_map>

namespace CoupledField
{

    // Define / declare logging stream
    DEFINE_LOG(preciceAdapter, "preciceAdapter")

    PreciceAdapter::PreciceAdapter(boost::shared_ptr<ParamNode> paramNode)
   
     : paramNode_(paramNode),
#ifdef USE_PRECICE
      participant_(nullptr),
#endif
      isInit_(false),
      configFileName_(""),
      participantName_(""),
      participantMeshName_(""),
      participantElemMeshName_(""),
      sequenceStep_(),
      nodeNumCoordMap_(),
      cfsNodeNumsVec_(),
      flatCoords_(),
      preciceNodeNumsVec_(),
      elemNumCoordMap_(),
      cfsElemNumsVec_(),
      flatElemCoords_(),
      preciceElemNumsVec_(),
      needElemMesh_(false),
      activeParticipantConfig_(),
      runtimeReadResults_(),
      runtimeWriteResults_(),
      rank_(0),
      size_(1),
      domain_(nullptr),
      solveStep_(nullptr),
      singlePDE_(nullptr)
    {
    }

    PreciceAdapter::~PreciceAdapter()
    {
        finalize();
    }


    // Helper function: add a unique node to the node container and map
    static void addUniqueNode(GridCFS* gridCFS, int nodenum,
                            std::vector<int>& nodeVec,
                            std::map<unsigned int, Vector<double>>& nodeMap) {
        if (nodeMap.find(nodenum) == nodeMap.end()) {
            nodeVec.push_back(nodenum);
            Vector<double> coord;
            gridCFS->GetNodeCoordinate(coord, nodenum, true);
            nodeMap[nodenum] = coord;
        }
    }

    // Helper function: add a unique element to the element container and map
    static void addUniqueElem(GridCFS* gridCFS, int elemNum,
                            std::vector<int>& elemVec,
                            std::map<unsigned int, Vector<double>>& elemMap) {
        if (elemMap.find(elemNum) == elemMap.end()) {
            elemVec.push_back(elemNum);
            Vector<double> center;
            gridCFS->GetElemCentroid(center, elemNum, true);
            elemMap[elemNum] = center;
        }
    }


    void PreciceAdapter::readConfigurationParameters() {
        // Retrieve parameters from paramNode_ (participantName_, participantMeshName_, etc.)
        paramNode_->Get("fileFormats/preciceCoupling/participantName")->GetValue("name", participantName_);
        paramNode_->Get("fileFormats/preciceCoupling/participantMeshName")->GetValue("name", participantMeshName_);
        paramNode_->Get("fileFormats/preciceCoupling/precice_configFile")->GetValue("file", configFileName_);
        sequenceStep_ = paramNode_->Get("fileFormats/preciceCoupling")->Get("sequenceStep")->As<int>();

        std::cout << "PreciceAdapter: Configuration parameters:" << "\n";
        std::cout << "  Participant Name: " << participantName_ << "\n";
        std::cout << "  Participant Mesh Name: " << participantMeshName_ << "\n";
        std::cout << "  Precice Config File: " << configFileName_ << "\n";
        std::cout << "  Sequence Step: " << sequenceStep_ << "\n";
    }

    void PreciceAdapter::readPreciceConfiguration() {
        std::cout << "PreciceAdapter: reading precice configuration from file: " << configFileName_ << "\n";
        // Use PreciceConfigReader to read additional configuration.
        try {
            PreciceConfigReader configReader(configFileName_);
            const auto &participants = configReader.getParticipants();
            bool found = false;
            for (const auto &pc : participants) {
                if (pc.name == participantName_) {
                    activeParticipantConfig_ = pc; // Store the relevant configuration
                    // now let's check if there is a _elem mesh in the provide meshes
                    for (const auto &mesh : pc.provideMeshes) {
                        // Check if mesh.name ends with "_elem"
                        if (mesh.name.size() >= 5 &&
                            mesh.name.substr(mesh.name.size() - 5) == "_elem") {
                            needElemMesh_ = true;
                            break;
                        }
                    }
                    LOG_DBG(preciceAdapter) << "Loaded configuration for participant: " << pc.name;
                    found = true;
                    break;
                }
            }
            if (!found) {
                EXCEPTION("Participant " << participantName_ << " not found in precice-config.xml");
            }
        } catch (const std::exception &e) {
            EXCEPTION("Error reading precice-config.xml: " << e.what());
        }
    }


    std::tuple<std::string, SolutionType> PreciceAdapter::convertResultNamesToCFS(const std::string& precicename) {
        static const std::unordered_map<std::string, std::tuple<std::string, SolutionType>> conversionMap = {
            {"Temperature", {"heatTemperature", SolutionType::HEAT_TEMPERATURE}},
            {"Displacement", {"mechDisplacement", SolutionType::MECH_DISPLACEMENT}},
            {"MagneticFluxDensity", {"magFluxDensity", SolutionType::MAG_FLUX_DENSITY}},
            {"MagneticFieldIntensity", {"magFieldIntensity", SolutionType::MAG_FIELD_INTENSITY}},
            {"JouleLossDensity", {"magJouleLossPowerDensity", SolutionType::MAG_JOULE_LOSS_POWER_DENSITY}},
            {"Pressure", {"acouRhsLoad", SolutionType::ACOU_RHS_LOAD}},
            {"PressureTemporalDerivative", {"acouRhsLoad", SolutionType::ACOU_RHS_LOAD}}
        };

        auto it = conversionMap.find(precicename);
        if (it != conversionMap.end())
            return it->second;
        else
            EXCEPTION("Invalid quantity: " << precicename
                    << ". Currently the adapter only works for one of [\"Temperature\", \"Displacement\", "
                    << "\"MagneticFluxDensity\", \"MagneticFieldIntensity\", \"JouleLossDensity\"]");
    }


    void PreciceAdapter::createPreciceParticipant() {
        std::cout << "PreciceAdapter: creating PreCICE participant..." << "\n";
        if (domain_->GetGridMap().size() > 1)
        {
            EXCEPTION("PreciceAdapter: works only with one grid per simulation - until now");
        }

        // add a check that we dont have a multidriver and that it's transient
        if( !dynamic_cast<TransientDriverPrecice*>(domain_->GetSingleDriver()) )
        {
            EXCEPTION("PreciceAdapter: this only works when using a transient simulation")
        }

        // Create the PreCICE participant instance.
        participant_ = std::make_unique<precice::Participant>(participantName_, configFileName_, rank_, size_);

        std::cout << "PreciceAdapter: created PreCICE participant '" << participantName_
                  << "' with config file '" << configFileName_ << "'" << "\n";

        // Create the relevant data structures for storing data of the participant (read-data and write-data)
        for (const auto &entry : activeParticipantConfig_.writeData) {
            ResultConfig config;
            config.precicename = entry.name;
            auto [tempName, tempSolutionType] = convertResultNamesToCFS(config.precicename);
            config.cfsname = tempName;
            config.solutiontype = tempSolutionType;
            config.meshName = entry.mesh;
            config.quantitydim = (needElemMesh_) ?
                participant_->getDataDimensions(participantMeshName_ + "_elem", config.precicename) :
                participant_->getDataDimensions(participantMeshName_, config.precicename);

            std::cout << "PreciceAdapter: Preparing to register write-data entry: " << config.precicename
                      << " on mesh: " << config.meshName
                      << " with dimension: " << config.quantitydim << "\n";

            std::cout << "CFSName is: " << config.cfsname << "\n";

            // Decide whether to create a node- or element-based result.
            // For example, if you decide based on the mesh name:
            if (config.meshName.find("_elem") != std::string::npos) {
                runtimeWriteResults_.push_back(std::make_unique<ElementResult>(config));
                // Allocate flat data based on the number of CFS elements.
                runtimeWriteResults_.back()->allocateData(cfsElemNumsVec_.size());
            } else {
                runtimeWriteResults_.push_back(std::make_unique<NodeResult>(config));
                runtimeWriteResults_.back()->allocateData(cfsNodeNumsVec_.size());
            }
        }

        std::cout << "PreciceAdapter: Registered " << runtimeWriteResults_.size() << " write-data entries." << "\n";

        for (const auto &entry : activeParticipantConfig_.readData) {
            ResultConfig config;
            config.precicename = entry.name;
            auto [tempName, tempSolutionType] = convertResultNamesToCFS(config.precicename);
            config.cfsname = tempName;
            config.solutiontype = tempSolutionType;
            config.meshName = entry.mesh;
            config.quantitydim = participant_->getDataDimensions(participantMeshName_, config.precicename);

            std::cout << "PreciceAdapter: Preparing to register read-data entry: " << config.precicename
                      << " on mesh: " << config.meshName
                      << " with dimension: " << config.quantitydim << "\n";

            std::cout << "CFSName is: " << config.cfsname << "\n";

            if (config.meshName.find("_elem") != std::string::npos) {
                runtimeReadResults_.push_back(std::make_unique<ElementResult>(config));
                runtimeReadResults_.back()->allocateData(cfsElemNumsVec_.size());
            } else {
                // EXCEPTION("PreciceAdapter: reading nodal values is not tested!");
                runtimeReadResults_.push_back(std::make_unique<NodeResult>(config));
                runtimeReadResults_.back()->allocateData(cfsNodeNumsVec_.size());
                std::cout << "PreciceAdapter: Registered " << runtimeReadResults_.size() << " read-data entries." << "\n";
            }
        }

        // this is not a joke! we need this because precice's setMeshVertices method changes
        // the last argument! by doing this, we have our original node number vector from cfs
        // and the adapted one from precice
        preciceNodeNumsVec_ = cfsNodeNumsVec_;
        participant_->setMeshVertices(participantMeshName_, flatCoords_, preciceNodeNumsVec_);

        // Similarly, register the element mesh with preCICE.
        // Derive the element mesh name by appending "_elem" to the participant mesh name.
        // but only if the "node mesh name"+"_elem" is specified in the precice config
        if(needElemMesh_)
        {
            participantElemMeshName_ = participantMeshName_ + "_elem";
            preciceElemNumsVec_ = cfsElemNumsVec_;
            participant_->setMeshVertices(participantElemMeshName_, flatElemCoords_, preciceElemNumsVec_);
        }

        std::cout << "PreciceAdapter: Registered mesh '" << participantMeshName_ << "' with "
                  << cfsNodeNumsVec_.size() << " nodes." << "\n";

    }

    void PreciceAdapter::setupMeshAndCoordinates() {
        // Get grid from the cfs domain and flatten node coordinates.
        std::cout << "PreciceAdapter: setting up mesh and coordinates..." << "\n";
        GridCFS* gridCFS = dynamic_cast<GridCFS*>(domain_->GetGrid());
        if (!gridCFS) {
            EXCEPTION("PreciceAdapter: Domain's grid is not of type GridCFS as expected.");
        }

        // get the region list from the domain tag in the openCFS simulation xml file
        ParamNodeList regionListNode = paramNode_->Get("domain")->Get("regionList")->GetChildren();
        
        if(regionListNode.GetSize() > 0) 
        {
            for( UInt i = 0; i < regionListNode.GetSize(); i++ ) 
            {
            std::string regionName = regionListNode[i]->Get("name")->As<std::string>(); 
            RegionIdType regionId = gridCFS->GetRegion().Parse( regionName );
            StdVector<UInt> elems, nodes;
            gridCFS->GetElementsByRegion(elems, regionId);
            gridCFS->GetNodesByRegion(nodes, regionId);
            shared_ptr<EntityList> elemlist = gridCFS->GetEntityList(EntityList::ListType::ELEM_LIST, regionName);
            shared_ptr<EntityList> nodelist = gridCFS->GetEntityList(EntityList::ListType::NODE_LIST, regionName);

            // Nodes
            EntityIterator it = nodelist->GetIterator();
            for (it.Begin(); !it.IsEnd(); it++) {
                int nodenum = it.GetNode();
                addUniqueNode(gridCFS, nodenum, cfsNodeNumsVec_, nodeNumCoordMap_);
            }

            // Elements
            it = elemlist->GetIterator();
            for (it.Begin(); !it.IsEnd(); it++) {
                const Elem * el = it.GetElem();
                addUniqueElem(gridCFS, el->GetElemNum(), cfsElemNumsVec_, elemNumCoordMap_);
            }

            }
        }
                

        // Flatten node coordinates into a contiguous vector.
        int cfsGridDim = gridCFS->GetDim();
        flatCoords_.resize(cfsNodeNumsVec_.size() * cfsGridDim);
        for (size_t i = 0; i < cfsNodeNumsVec_.size(); ++i) {
            Vector<double> point = nodeNumCoordMap_[cfsNodeNumsVec_[i]];
            for (int j = 0; j < cfsGridDim; ++j) {
                flatCoords_[i * cfsGridDim + j] = point[j];
            }
        }

        // Flatten element midpoint coordinates.
        flatElemCoords_.resize(cfsElemNumsVec_.size() * cfsGridDim);
        for (size_t i = 0; i < cfsElemNumsVec_.size(); ++i) {
            Vector<double> midPoint = elemNumCoordMap_[cfsElemNumsVec_[i]];
            for (int j = 0; j < cfsGridDim; ++j) {
                flatElemCoords_[i * cfsGridDim + j] = midPoint[j];
            }
        }

    }


    void PreciceAdapter::initialize(Domain *domain, SinglePDE* pde) {
#ifdef USE_PRECICE
    domain_ = domain;

    // Step 1: Retrieve configuration parameters from the openCFS ParamNode.
    readConfigurationParameters();

    // Step 2: Parse the precice-config.xml file 
    readPreciceConfiguration();

    // Step 3: Extract and flatten node coordinates.
    setupMeshAndCoordinates();
    
    // maybe we are not yet in the correct sequence step, therefore,
    // continue
    TransientDriverPrecice* tp = dynamic_cast<TransientDriverPrecice*>(domain_->GetSingleDriver());

    if(!tp){
        return;
    }else{
        if(!(tp->GetActSequenceStep() == sequenceStep_)){
            std::cout << "PreciceAdapter: Not the correct sequence step, skipping precice initialization" << "\n";
            return;
        }
    }

    // Step 4: Create the PreCICE participant.
    if(!isInit_){
        std::cout << "Not yet initialized" << "\n";

        createPreciceParticipant();

        std::cout << "PreciceAdapter: Precice participant created." << "\n";
        
        // This is such a DIRTY HACK!!!! This should actually be set by precice!!!!!!!!!
        // This ONLY WORKS for initializing the "Temperature"
        if(participant_->requiresInitialData()){
            std::cout << "PreciceAdapter: Participant requires initial data - initializing 'Temperature' to 293.15K" << "\n";
            for (auto &result : runtimeWriteResults_) {
                if(result->getConfig().precicename == "Temperature"){
                    if (result->getResultType() == ResultType::NODE) {
                        NodeResult* nr = dynamic_cast<NodeResult*>(result.get());
                        if (nr) {
                            std::vector<double> initVec(nr->getFlatData().size());
                            std::fill(initVec.begin(), initVec.end(), 293.15);
                            participant_->writeData(nr->getConfig().meshName, nr->getConfig().precicename,
                                                    preciceNodeNumsVec_, initVec);
                        }
                    } else { // ELEMENT-based result.
                        ElementResult* er = dynamic_cast<ElementResult*>(result.get());
                        if (er) {
                            std::vector<double> initVec(er->getFlatData().size());
                            std::fill(initVec.begin(), initVec.end(), 293.15);
                            participant_->writeData(er->getConfig().meshName, er->getConfig().precicename,
                                                    preciceElemNumsVec_, initVec);
                        }
                    }
                }
            }
        }

        participant_->initialize();
        isInit_ = true;
    }

    RegisterSinglePDE(pde);

    RegisterExternalResults();
    
#endif
    }



    void PreciceAdapter::RegisterSolveStep(BaseSolveStep *solveStep){
        // Try to cast to StdSolveStep*
        StdSolveStep* stdSolveStep = dynamic_cast<StdSolveStep*>(solveStep);
        if (!stdSolveStep) {
            // The cast failed: solveStep is not a StdSolveStep
            EXCEPTION("PreciceAdapter::RegisterSolveStep expected a StdSolveStep pointer.");
        }
        solveStep_ = stdSolveStep;
    }

    void PreciceAdapter::RegisterSinglePDE(SinglePDE* pde){
        singlePDE_ = pde;
        std::cout << "PreciceAdapter: Registered SinglePDE." << "\n";
    }


    void PreciceAdapter::RegisterTimeStepReadData(){
        if(!singlePDE_->IsInitialized()){
            // initial state case
            return;
        }
        // could happen if we are not in the right sequencestep
        TransientDriverPrecice* tp = dynamic_cast<TransientDriverPrecice*>(domain_->GetSingleDriver());
        if(!tp){
            return;
        }else{
            if(!(tp->GetActSequenceStep() == sequenceStep_)){
                return;
            }
        }

       
        //TODO handle nodal and element cases. Maybe handle it based on the particular result
        // like expectResultType nodes or elements
        // loop over the required results that we need - defined in the precice config
        for (auto &result : runtimeReadResults_) {
            LOG_DBG(preciceAdapter) << "Participant " << activeParticipantConfig_.name
                                    << " reads data: " << result->getConfig().cfsname
                                    << " on mesh: " << result->getConfig().meshName;

            // Depending on whether the result is node- or element-based, call participant_->readData
            if (result->getResultType() == ResultType::ELEMENT) {
                // Resize flat data based on number of elements
                result->allocateData(cfsElemNumsVec_.size());
                participant_->readData(result->getConfig().meshName,
                                        result->getConfig().precicename,
                                        preciceElemNumsVec_,
                                        dynamic_cast<TransientDriverPrecice*>(domain_->GetSingleDriver())->GetDeltaT(),
                                        const_cast<std::vector<double>&>(result->getFlatData()));
                // Let the result object convert the flat data into its internal map.
                result->readData(cfsElemNumsVec_,
                                dynamic_cast<TransientDriverPrecice*>(domain_->GetSingleDriver())->GetDeltaT());
            }
            else { // Node-based result
                result->allocateData(cfsNodeNumsVec_.size());
                participant_->readData(result->getConfig().meshName,
                                        result->getConfig().precicename,
                                        preciceNodeNumsVec_,
                                        dynamic_cast<TransientDriverPrecice*>(domain_->GetSingleDriver())->GetDeltaT(),
                                        const_cast<std::vector<double>&>(result->getFlatData()));
                result->readData(cfsNodeNumsVec_,
                                dynamic_cast<TransientDriverPrecice*>(domain_->GetSingleDriver())->GetDeltaT());
            }
        }


        // ================================== START: Update OpenCFS result contexts

        ResultHandler* resHandler = domain_->GetResultHandler();
        auto* resultContextsPtr = resHandler->GetResultContexts();

        // Iterate over each result context from OpenCFS
        for (auto &entry : *resultContextsPtr) {
            shared_ptr<BaseResult> baseResult = entry.first;
            shared_ptr<ResultHandler::ResultContext> resultContext = entry.second;
            // Compare the result name with the configuration of our result objects.
            std::string cfsResultName = baseResult->GetResultInfo()->resultName;
            
            // Loop over our runtime read results to see if one matches.
            for (const auto &result : runtimeReadResults_) {
                if (result->getConfig().cfsname == cfsResultName) {
                    // We assume here that we are dealing with an element-based result.
                    // (A similar branch would be added for node-based results.)
                    if (result->getResultType() == ResultType::ELEMENT) {
                        Result<Double>& actSol = static_cast<Result<Double>&>(*baseResult);
                        EntityIterator it = actSol.GetEntityList()->GetIterator();
                        Vector<Double>& vec = actSol.GetVector();
                        vec.Resize(it.GetSize() * result->getConfig().quantitydim);
                        ElementResult* elemResult = dynamic_cast<ElementResult*>(result.get());
                        if (!elemResult) {
                            EXCEPTION("Expected an ElementResult for element-based data.");
                        }
                        const auto& elemMap = elemResult->getElementResultMap();
                        for (it.Begin(); !it.IsEnd(); it++) {
                            const Elem* el = it.GetElem();
                            auto search = elemMap.find(el->GetElemNum());
                            if (search != elemMap.end()) {
                                const Vector<double>& tempField = search->second;
                                for (UInt iDim = 0; iDim < result->getConfig().quantitydim; iDim++) {
                                    vec[it.GetPos() * result->getConfig().quantitydim + iDim] = tempField[iDim];
                                }
                            }
                        }
                    }
                }
            }
        }


    }



    void PreciceAdapter::RegisterTimeStepWriteData(){
        TransientDriverPrecice* tp = dynamic_cast<TransientDriverPrecice*>(domain_->GetSingleDriver());
        if(!tp){
            return;
        }else{
            if(!(tp->GetActSequenceStep() == sequenceStep_)){
                return;
            }
        }

        // Get the result contexts from the OpenCFS result handler.
        ResultHandler* resHandler = domain_->GetResultHandler();
        auto* resultContextsPtr = resHandler->GetResultContexts();

        // Build a vector of contexts using boost::shared_ptr.
        std::vector<std::pair<boost::shared_ptr<BaseResult>, boost::shared_ptr<ResultHandler::ResultContext>>> contexts;
        for (auto &entry : *resultContextsPtr) {
            contexts.push_back(std::make_pair(entry.first, entry.second));
        }

        // Process each write result.
        for (auto &result : runtimeWriteResults_) {
            if (result->getResultType() == ResultType::NODE) {
                NodeResult* nr = dynamic_cast<NodeResult*>(result.get());
                std::cout << "WE ACTUALLY HAVE NODE\n"; 
                if (nr) {
                    nr->updateFromOpenCFS(contexts, cfsNodeNumsVec_);
                    nr->writeData(cfsNodeNumsVec_);
                    participant_->writeData(nr->getConfig().meshName, nr->getConfig().precicename,
                                            preciceNodeNumsVec_, nr->getFlatData());
                }
            } else { // ELEMENT-based result.
                ElementResult* er = dynamic_cast<ElementResult*>(result.get());
                std::cout << "WE ACTUALLY HAVE ELEMENT\n"; 
                if (er) {
                    er->updateFromOpenCFS(contexts, cfsElemNumsVec_);
                    er->writeData(cfsElemNumsVec_);
                    participant_->writeData(er->getConfig().meshName, er->getConfig().precicename,
                                            preciceElemNumsVec_, er->getFlatData());
                }
            }
        }

    }





    void PreciceAdapter::MarkReadResultsUpdated()
    {
#ifdef USE_PRECICE
        if (runtimeReadResults_.empty()) return;

        ResultHandler* resHandler = domain_->GetResultHandler();
        auto* resultContextsPtr = resHandler->GetResultContexts();

        for (auto &entry : *resultContextsPtr) {
            shared_ptr<BaseResult> baseResult = entry.first;
            std::string cfsResultName = baseResult->GetResultInfo()->resultName;

            for (const auto &result : runtimeReadResults_) {
                if (result->getConfig().cfsname == cfsResultName) {
                    std::cout << "PreciceAdapter::MarkReadResultsUpdated: marking '"
                              << cfsResultName << "' as updated\n";
                    resHandler->UpdateResult(baseResult);
                    break;
                }
            }
        }
#endif
    }

    void PreciceAdapter::RegisterExternalResults()
    {
    std::cout << "PreciceAdapter: Registering external results with OpenCFS result handler..." << "\n";
    #ifdef USE_PRECICE
    // Get the result handler from the domain.
    ResultHandler* resHandler = domain_->GetResultHandler();

    std::cout << "ResultHandler initialised" << "\n";

    if(!resHandler)
        EXCEPTION("PreciceAdapter::RegisterElementResults: No result handler available.");

    std::cout << "runtimeReadResults.size() = " << runtimeReadResults_.size() << "\n";

    if(runtimeReadResults_.size() == 0) return;

    GridCFS* gridCFS = dynamic_cast<GridCFS*>(domain_->GetGrid());
    //solveStep_

    // get the entity list of the entities on which the pde is defined on
    StdVector<RegionIdType> regions = singlePDE_->GetRegions();

    std::cout << "regions.ToString() (PreciceAdapter.cc): " << regions.ToString() << "\n";
    
    shared_ptr<EntityList> list = gridCFS->GetEntityList(EntityList::ListType::ELEM_LIST, gridCFS->GetRegionName(regions[0]));
    //shared_ptr<EntityList> list = gridCFS->GetEntityList(EntityList::ListType::NODE_LIST, gridCFS->GetRegionName(regions[0]));

   std::cout << list << "\n";

   if(list){
        // For each runtime quantity that is element-based (i.e. does not have nodal results)
        for(auto &result : runtimeReadResults_) {
            // If nodal results are empty but we have element results, then we treat it as an element result.
                // Create a new shallow result object.
                shared_ptr<BaseResult> sol(new Result<Double>());
                sol->SetEntityList(list);
                
                // Create a new ResultInfo to describe this element result.
                shared_ptr<ResultInfo> ri(new ResultInfo());
                ri->resultName = result->getConfig().cfsname;
                ri->resultType = result->getConfig().solutiontype;

                ri->definedOn  = ResultInfo::ELEMENT;
                //ri->definedOn  = ResultInfo::NODE;
                // If the quantity has dimension 1 we use SCALAR; otherwise, we use VECTOR.
                ri->entryType  = (result->getConfig().quantitydim == 1) ? ResultInfo::SCALAR : ResultInfo::VECTOR;
                ri->dofNames = "";
                sol->SetResultInfo(ri);
                
                // Prepare additional parameters for registration.
                StdVector<std::string> outDest;
                outDest.Push_back("");
                
                // Register the result with the result handler.
                // The parameters here (0,0,1,1,outDest,"",true,false) mimic the style in GridCFS::CreateGridInformation.

                // resHandler->RegisterResult(sol, /*functor*/ nullptr, sequenceStep_, 0, 1,
                //                          domain_->GetSingleDriver()->GetNumSteps(),
                //                          outDest, "", true, false);
            
        }
   }


    #endif
    }

    Vector<Double> PreciceAdapter::GetElemResult(SolutionType solType, int elemNum)
    {
        for (const auto &result : runtimeReadResults_) {
            std::cout << "solutiontype: " << result->getConfig().solutiontype << "\n";
            std::cout <<  "resulttype: " << static_cast<int>(result->getResultType()) << "\n";
            if (result->getConfig().solutiontype == solType) {
                if (result->getResultType() == ResultType::ELEMENT) {
                    ElementResult* er = dynamic_cast<ElementResult*>(result.get());
                    if (er) {
                        const auto& mapRef = er->getElementResultMap();
                        auto it = mapRef.find(elemNum);
                        if (it != mapRef.end()) {
                            return it->second;
                        }
                    }
                }
            }
        }
        EXCEPTION("PreciceAdapter: no result with solution type " << solType << " found");
    }

    Vector<Double> PreciceAdapter::GetNodeResult(SolutionType solType, int nodeNum)
    {
        for (const auto &result : runtimeReadResults_) {
            if (result->getConfig().solutiontype == solType) {
                if (result->getResultType() == ResultType::NODE) {
                    NodeResult* er = dynamic_cast<NodeResult*>(result.get());
                    if (er) {
                        const auto& mapRef = er->getNodeResultMap();
                        auto it = mapRef.find(nodeNum);
                        if (it != mapRef.end()) {
                            return it->second;
                        }
                    }
                }
            }
        }
        EXCEPTION("PreciceAdapter: no result with solution type " << solType << " found");
    }

    void PreciceAdapter::finalize()
    {
#ifdef USE_PRECICE
        if (participant_)
        {
            participant_->finalize();
            participant_.reset();
            LOG_DBG(preciceAdapter) << "PreCICE participant finalized.";
        }
#endif
    }



} // end of namespace
