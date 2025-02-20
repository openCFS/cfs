
#include "IPreciceAdapter.hh"
#include "PreciceAdapter.hh"
#include "MinimalXmlParser.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "DataInOut/ResultHandler.hh"
#include "Domain/Domain.hh"
#include "Domain/Results/ResultFunctor.hh"
#include "Domain/Mesh/GridCFS/GridCFS.hh"
#include "Driver/TransientDriver.hh"
#include "Driver/SolveSteps/StdSolveStep.hh"
#include "PDE/SinglePDE.hh"
#include "Utils/StdVector.hh"
#include "MatVec/Vector.hh"
#include <numeric>

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
      runtimeReadQuantities_(),
      runtimeWriteQuantities_(),
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
    }

    void PreciceAdapter::readPreciceConfiguration() {
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


    std::tuple<std::string, SolutionType> PreciceAdapter::convertResultNamesToCFS(std::string precicename){
        std::string convertedName;
        if (precicename == "Temperature") {
            return std::make_tuple("heatTemperature", SolutionType::HEAT_TEMPERATURE) ;
        } else if (precicename == "Displacement") {
            return std::make_tuple("mechDisplacement", SolutionType::MECH_DISPLACEMENT);
        } else if (precicename == "MagneticFluxDensity") {
            return std::make_tuple("magFluxDensity", SolutionType::MAG_FLUX_DENSITY);
        } else if (precicename == "MagneticFieldIntensity") {
            return std::make_tuple("magFieldIntensity", SolutionType::MAG_FIELD_INTENSITY);
        } else if (precicename == "JouleLossDensity") {
            return std::make_tuple("magJouleLossPowerDensity", SolutionType::MAG_JOULE_LOSS_POWER_DENSITY);
        } else {
            EXCEPTION("Invalid quantity: " << precicename
                    << ". Currently the adapter only works for one of"
                    << "[\"Temperature\", \"Displacement\", \"MagneticFluxDensity\", \"MagneticFieldIntensity\",\"JouleLossDensity\",]");
        }
    }


    void PreciceAdapter::createPreciceParticipant() {
        if (domain_->GetGridMap().size() > 1)
        {
            EXCEPTION("PreciceAdapter: works only with one grid per simulation - until now");
        }

        // add a check that we dont have a multidriver and that it's transient
        if( !dynamic_cast<TransientDriver*>(domain_->GetSingleDriver()) )
        {
            EXCEPTION("PreciceAdapter: this only works when using a transient simulation")
        }

        // Create the PreCICE participant instance.
        participant_ = std::make_unique<precice::Participant>(participantName_, configFileName_, rank_, size_);

        // Create the relevant data structures for storing data of the participant (read-data and write-data)
        for (const auto &entry : activeParticipantConfig_.writeData) {
            PreciceRuntimeQuantity rq;
            rq.precicename = entry.name;
        
            // Validate and convert the precicename to an openCFS result name
            auto [tempName, tempSolutionType] = convertResultNamesToCFS(rq.precicename);
            rq.cfsname = tempName;
            rq.solutiontype = tempSolutionType;
            // Determine data dimensions based on whether an element mesh is needed
            rq.quantitydim = (needElemMesh_) ?
                participant_->getDataDimensions(participantMeshName_ + "_elem", rq.precicename) :
                participant_->getDataDimensions(participantMeshName_, rq.precicename);
            rq.meshName = entry.mesh;
            //rq.nodeNumbers = cfsNodeNumsVec_;
            // Reserve space for the actual exchanged data:
            rq.data.resize(cfsNodeNumsVec_.size() * rq.quantitydim);
            rq.type = Exchangetype::WRITE;        
            runtimeWriteQuantities_.push_back(rq);
        }


        for (const auto &entry : activeParticipantConfig_.readData) {
            PreciceRuntimeQuantity rq;
            rq.precicename = entry.name;
        
            // Validate and convert the precicename to an openCFS result name
            auto [tempName, tempSolutionType] = convertResultNamesToCFS(rq.precicename);
            rq.cfsname = tempName;
            rq.solutiontype = tempSolutionType;
            //rq.griddim = 0;
            rq.quantitydim = participant_->getDataDimensions(participantMeshName_, rq.precicename);
            rq.meshName = entry.mesh;
            //rq.nodeNumbers = {0};
            // Reserve space for the actual exchanged data:
            //rq.data.resize(cfsNodeNumsVec_.size() * rq.quantitydim);
            rq.type = Exchangetype::READ;        
            runtimeReadQuantities_.push_back(rq);
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

    }

    void PreciceAdapter::setupMeshAndCoordinates() {
        // Get grid from the cfs domain and flatten node coordinates.
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

            // Nodes
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
    if( !(dynamic_cast<TransientDriver*>(domain_->GetSingleDriver())->GetActSequenceStep() == sequenceStep_) ){ return;}

    // Step 4: Create the PreCICE participant.
    if(!isInit_){
        createPreciceParticipant();
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
    }


    void PreciceAdapter::RegisterTimeStepReadData(){
        // could happen if we are not in the right sequencestep
        if( !(dynamic_cast<TransientDriver*>(domain_->GetSingleDriver())->GetActSequenceStep() == sequenceStep_) ) return;

       
        //TODO handle nodal and element cases. Maybe handle it based on the particular result
        // like expectResultType nodes or elements
        // loop over the required results that we need - defined in the precice config
        for (auto &quantity : runtimeReadQuantities_) {
            LOG_DBG(preciceAdapter) << "Participant " << activeParticipantConfig_.name
                                    << " reads data: " << quantity.cfsname
                                    << " on mesh: " << quantity.meshName;

            // we need to resize the data vector
            quantity.data.resize(cfsElemNumsVec_.size() * quantity.quantitydim);
            participant_->readData(quantity.meshName, quantity.precicename, preciceElemNumsVec_,
                                dynamic_cast<TransientDriver*>(domain_->GetSingleDriver())->GetDeltaT(),
                                quantity.data);
            
            // for easier handling, we also fill the nodeResultMap
            // but keep in mind that the nodes in nodeResultMap are a point cloud and
            // do not have to be real node-values! They can also be element values
            for (std::size_t i = 0; i < preciceElemNumsVec_.size(); ++i) {
                int cfselem = cfsElemNumsVec_[i];
                Vector<double> vals;
                vals.Resize(quantity.quantitydim);
                for (size_t k = 0; k < static_cast<size_t>(quantity.quantitydim); ++k) {
                    vals[k] = quantity.data[i * quantity.quantitydim + k];
                }
                quantity.elemResultMap[cfselem] = vals;
            }



            // ================================== START TESTING
 
            //Get the result contexts from the opencfs result handler.
            ResultHandler * resHandler = domain_->GetResultHandler();
            auto* resultContextsPtr = resHandler->GetResultContexts();

            
            //loop over the required results that we need - defined in the precice config      
            
            //Process each matching result context to fill result maps
            for (auto &entry : *resultContextsPtr) {
                shared_ptr<BaseResult> baseResult = entry.first;
                shared_ptr<ResultHandler::ResultContext> resultContext = entry.second;
                Result<Double>& actSol = static_cast<Result<Double>& >(*baseResult);
                EntityIterator it = actSol.GetEntityList()->GetIterator();
                Vector<Double>& vec = actSol.GetVector();
                Vector<Double> tempField;
                vec.Resize( it.GetSize() * quantity.quantitydim );
                
                // here we just assume the results are EntityList::ELEM_LIST and SCALAR
                // TODO: we need to define somewhere in the xml file if the read
                // quantities are element or node based results

                // similar to have a look at ResultFunctor.cc lines 63ff
                std::string cfsResultName = baseResult->GetResultInfo()->resultName;
                if (cfsResultName == quantity.cfsname) {
                    for(it.Begin(); !it.IsEnd(); it++)
                        {
                        const Elem * el = it.GetElem();
                        // get result of element
                        tempField = quantity.elemResultMap[el->GetElemNum()];
                        // loop over dofs
                        for(UInt iDim = 0; iDim < quantity.quantitydim; iDim++ )
                            vec[it.GetPos()*quantity.quantitydim + iDim] = tempField[iDim];
                        }
                }
            }
            

            // ==================================  END TESTING

        }

    }


    void PreciceAdapter::RegisterTimeStepWriteData(){
        // could happen if we are not in the right sequencestep
        if( !(dynamic_cast<TransientDriver*>(domain_->GetSingleDriver())->GetActSequenceStep() == sequenceStep_) ) return;
        
        // Get the result contexts from the opencfs result handler.
        ResultHandler * resHandler = domain_->GetResultHandler();
        auto* resultContextsPtr = resHandler->GetResultContexts();


        // loop over the required results that we need - defined in the precice config
        for (auto &quantity : runtimeWriteQuantities_) {
            // loop over the provided results from openCFS
            // Process each matching result context to fill result maps
            for (auto &entry : *resultContextsPtr) {
                shared_ptr<BaseResult> baseResult = entry.first;
                shared_ptr<ResultHandler::ResultContext> resultContext = entry.second;

                BaseResult &actResult = *(resultContext->result);
                std::string cfsResultName = actResult.GetResultInfo()->resultName;
                if (cfsResultName == quantity.cfsname) {
                    quantity.available = true;
                    // This is more or less a c&p from FeFunction<T>::ExtractResult( shared_ptr<BaseResult> res )
                    shared_ptr<EntityList> list = actResult.GetEntityList();
                    EntityIterator it = list->GetIterator();
                    EntityList::ListType entityListType = actResult.GetEntityList()->GetType();
                    int quantityDim = quantity.quantitydim;
                    // quick sanity check
                    if (quantityDim != static_cast<int>(actResult.GetResultInfo()->GetDimDof()))
                    {
                        EXCEPTION("PreciceAdapter: the requested quantity should have dimension "
                                << quantityDim << " but openCFS provided a result with dimension "
                                << actResult.GetResultInfo()->GetDimDof());
                    }
                    switch (entityListType) {
                        case EntityList::ELEM_LIST:
                        case EntityList::SURF_ELEM_LIST:
                        {
                            // in FieldCoefFunctor<TYPE>::EvalResult this is handled in a separate
                            // way compared to nodal evaluation. The reason for that is that element
                            // values, e.g., mech stress tensor, is a derived quantity and it's evaluated
                            // in the cell/element center

                            FieldCoefFunctor<double>* fcfunctor = dynamic_cast<FieldCoefFunctor<double>*>(resultContext->functor.get());
                            if (!fcfunctor)
                            {
                                EXCEPTION("PreciceAdapter: no downcast to FieldCoefFunctor possible");
                            }
                            resultContext->functor->EvalResult(resultContext->result);
                            for (it.Begin(); !it.IsEnd(); it++)
                            {
                                const Elem *el = it.GetElem();
                                LocPoint lp = Elem::shapes[el->type].midPointCoord;
                                LocPointMapped lpm;
                                Vector<double> tempField;
                                shared_ptr<ElemShapeMap> esm = it.GetGrid()->GetElemShapeMap(el, true);
                                lpm.Set(lp, esm, 0.0);
                                fcfunctor->GetVector(tempField, lpm);
                                // LOG_DBG(preciceAdapter) << "element "<<el->elemNum
                                //                         <<" has center " << lpm.GetGlobal().ToString()
                                //                         <<" and result "<<tempField.ToString();
                                quantity.elemResultMap[el->GetElemNum()] = tempField;
                            }
                            break;
                        }
                        case EntityList::NODE_LIST: {
                            Vector<double>* res_vec = dynamic_cast<Vector<Double>*>(actResult.GetSingleVector());
                            if (!res_vec)
                                EXCEPTION("PreciceAdapter: Unable to retrieve node result vector for opencfs result " << quantity.cfsname);
                            int pos = 0;
                            for (it.Begin(); !it.IsEnd(); it++) {
                                // we could either go via the FeSpace and get the equation number
                                // for the specific entity but since we know how the actSol
                                // was generated we can skip this re-organization step
                                // Have a look at void FeFunction<T>::ExtractResult( shared_ptr<BaseResult> res )
                                // there the connection to solution (res_vec) and eqnNumbers is set
                                int node = it.GetNode();
                                Vector<double> resultForNode(quantityDim);
                                for (int k = 0; k < quantityDim; ++k) {
                                    resultForNode[k] = (*res_vec)[pos * quantityDim + k];
                                }
                                // LOG_DBG(preciceAdapter)
                                //     << "Reading result for node " << node
                                //     << " : " << resultForNode.ToString();
                                quantity.nodeResultMap[node] = resultForNode;
                                ++pos;
                            }
                            break;
                        }
                        default:
                            EXCEPTION("PreciceAdapter::RegisterTimeStep: Unhandled entity list type.");
                            break;
                    }
                }
            }
            if (!quantity.available)
                EXCEPTION("PreciceAdapter: result " << quantity.cfsname << " was requested from openCFS but not provided");

            LOG_DBG(preciceAdapter) << "Participant " << activeParticipantConfig_.name
                                    << " writes data: " << quantity.precicename
                                    << " on mesh: " << quantity.meshName;

            // Process nodal vs. element results:
            if (!quantity.nodeResultMap.empty()) {
                quantity.data.resize(cfsNodeNumsVec_.size() * quantity.quantitydim);
                for (size_t i = 0; i < cfsNodeNumsVec_.size(); ++i) {
                    int node = cfsNodeNumsVec_[i];
                    Vector<double> nodeResult = quantity.nodeResultMap[node];
                    for (int k = 0; k < quantity.quantitydim; ++k) {
                        quantity.data[i * quantity.quantitydim + k] = nodeResult[k];
                    }
                }
                participant_->writeData(quantity.meshName, quantity.precicename,
                                        preciceNodeNumsVec_, quantity.data);
            } else {
                if (!needElemMesh_)
                    EXCEPTION("PreciceAdapter: you want to write an element result but did not specify the xxx_elem mesh in provided-mesh");

                LOG_DBG(preciceAdapter) << "Element results available: " << quantity.elemResultMap.size();
                LOG_DBG(preciceAdapter) << "Size of preciceElemNumsVec_: " << preciceElemNumsVec_.size();

                quantity.data.resize(quantity.elemResultMap.size() * quantity.quantitydim);
                
                for(size_t i = 0; i < cfsElemNumsVec_.size(); ++i){
                    int cfselemnum = cfsElemNumsVec_[i];
                    int preciceelemnum = preciceElemNumsVec_[i];
                    Vector<double>& value = quantity.elemResultMap[cfselemnum];
                    for (int k = 0; k < quantity.quantitydim; ++k) {
                        quantity.data[preciceelemnum * quantity.quantitydim + k] = value[k];
                    }
                }

                participant_->writeData(quantity.meshName, quantity.precicename,
                                        preciceElemNumsVec_, quantity.data);

            //     // logging
            //     for (size_t i = 0; i < preciceElemNumsVec_.size(); ++i) {
            //         int cfsElem = cfsElemNumsVec_[i];
            //         int preciceElem = preciceElemNumsVec_[i];
            //         std::ostringstream oss;
            //         oss << "CFS Element: " << cfsElem 
            //             << ", PreCICE Element: " << preciceElem 
            //             << ", Result: ";
            //         // For each element, output all components if vector result.
            //         for (int d = 0; d < quantity.quantitydim; ++d) {
            //             oss << quantity.data[i * quantity.quantitydim + d] << " ";
            //         }
            //     LOG_DBG(preciceAdapter) << oss.str();
            // }
            }
        }

        double cfs_dt = dynamic_cast<TransientDriver*>(domain_->GetSingleDriver())->GetDeltaT();
        double precice_dt = this->participant_->getMaxTimeStepSize();
        double dt = cfs_dt <= precice_dt? cfs_dt : precice_dt;
        this->participant_->advance(dt);
    }



    void PreciceAdapter::RegisterExternalResults()
    {
    #ifdef USE_PRECICE
    // Get the result handler from the domain.
    ResultHandler* resHandler = domain_->GetResultHandler();
    if(!resHandler)
        EXCEPTION("PreciceAdapter::RegisterElementResults: No result handler available.");

    if(runtimeReadQuantities_.size() == 0) return;

    GridCFS* gridCFS = dynamic_cast<GridCFS*>(domain_->GetGrid());
    //solveStep_

    // get the entity list of the entities on which the pde is defined on
    StdVector<RegionIdType> regions = singlePDE_->GetRegions();
    

   shared_ptr<EntityList> list = gridCFS->GetEntityList(EntityList::ListType::ELEM_LIST, gridCFS->GetRegionName(regions[0]));
   if(!list)
       EXCEPTION("PreciceAdapter::RegisterElementResults: Could not create element entity list.");

    // For each runtime quantity that is element-based (i.e. does not have nodal results)
    for(auto &quantity : runtimeReadQuantities_) {
        // If nodal results are empty but we have element results, then we treat it as an element result.
            // Create a new shallow result object.
            shared_ptr<BaseResult> sol(new Result<Double>());
            sol->SetEntityList(list);
            
            // Create a new ResultInfo to describe this element result.
            shared_ptr<ResultInfo> ri(new ResultInfo());
            ri->resultName = quantity.cfsname;
            ri->resultType = quantity.solutiontype;
            ri->definedOn  = ResultInfo::ELEMENT;
            // If the quantity has dimension 1 we use SCALAR; otherwise, we use VECTOR.
            ri->entryType  = (quantity.quantitydim == 1) ? ResultInfo::SCALAR : ResultInfo::VECTOR;
            ri->dofNames = "";
            sol->SetResultInfo(ri);
            
            // Prepare additional parameters for registration.
            StdVector<std::string> outDest;
            outDest.Push_back("");
            
            // Register the result with the result handler.
            // The parameters here (0,0,1,1,outDest,"",true,false) mimic the style in GridCFS::CreateGridInformation.
            resHandler->RegisterResult(sol, /*functor*/ nullptr, sequenceStep_, 0, 1,
                                       domain_->GetSingleDriver()->GetNumSteps(),
                                       outDest, "", true, false);
        
    }
    #endif
    }

    Vector<Double> PreciceAdapter::GetElemResult(SolutionType solType, int elemNum)
    {
        bool found = false;
        for(auto &quantity : runtimeReadQuantities_)
        {
            if(quantity.solutiontype == solType)
            {
                found = true;
                return quantity.elemResultMap[elemNum];
            }            

        }
        if(!found)
            EXCEPTION("PreciceAdapter: no result with solution type "<<solType<<" found");

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
