
#include "IPreciceAdapter.hh"
#include "PreciceAdapter.hh"
#include "MinimalXmlParser.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "DataInOut/ResultHandler.hh"
#include "Domain/Domain.hh"
#include "Domain/Mesh/GridCFS/GridCFS.hh"
#include "Driver/TransientDriver.hh"
#include "Driver/SolveSteps/StdSolveStep.hh"
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
          rank_(0),
          size_(1),
          configFileName_(""),
          participantName_(""),
          participantMeshName_(""),
          cfsNodeNumsVec_(0)
    {
    }

    PreciceAdapter::~PreciceAdapter()
    {
        finalize();
    }


    void PreciceAdapter::readConfigurationParameters() {
        // Retrieve parameters from paramNode_ (participantName_, participantMeshName_, etc.)
        paramNode_->Get("fileFormats/preciceCoupling/participantName")->GetValue("name", participantName_);
        paramNode_->Get("fileFormats/preciceCoupling/participantMeshName")->GetValue("name", participantMeshName_);
        paramNode_->Get("fileFormats/preciceCoupling/precice_configFile")->GetValue("file", configFileName_);
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

        /*
            In case I need to have a look at the participant configuration, I can use
            for (const auto &pm : activeParticipantConfig_.provideMeshes) {
                LOG_DBG(preciceAdapter) << "Participant " << activeParticipantConfig_.name
                                        << " provides mesh: " << pm.name;
            }
            or others like    activeParticipantConfig_.{provideMeshes, receiveMeshes, readData, writeData}
            }
        */
    }


    std::string PreciceAdapter::convertResultNamesToCFS(std::string precicename){
        std::string convertedName;
        if (precicename == "Temperature") {
            convertedName = "heatTemperature";
        } else if (precicename == "Displacement") {
            convertedName = "mechDisplacement";
        } else {
            EXCEPTION("Invalid quantity: " << precicename
                    << ". Currently the adapter only works for one of [\"Temperature\", \"Displacement\"]");
        }
        return convertedName;
    }


    void PreciceAdapter::createPreciceParticipant() {
        if (domain_->GetGridMap().size() > 1)
        {
            EXCEPTION("PreciceAdapter: works only with one grid per simulation - until now");
        }

        // add a check that we dont have a multidriver and that it's transient
        if( !dynamic_cast<TransientDriver*>(domain_->GetSingleDriver()) ){
            EXCEPTION("PreciceAdapter: this only works when using a transient simulation")
        }

        // Create the PreCICE participant instance.
        participant_ = std::make_unique<precice::Participant>(participantName_, configFileName_, rank_, size_);

        // Create the relevant data structures for storing data of the participant (read-data and write-data)
        for (const auto &entry : activeParticipantConfig_.writeData) {
            PreciceRuntimeQuantity rq;
            rq.precicename = entry.name;
        
            // Validate and convert the precicename to an openCFS result name
            rq.cfsname = convertResultNamesToCFS(rq.precicename);
            // specified grid dimension in precice config
            //rq.griddim = participant_->getMeshDimensions(participantMeshName_);
            rq.quantitydim = participant_->getDataDimensions(participantMeshName_, rq.precicename);
            rq.meshName = entry.mesh;
            //rq.nodeNumbers = cfsNodeNumsVec_;
            // Reserve space for the actual exchanged data:
            rq.data.resize(cfsNodeNumsVec_.size() * rq.quantitydim);
            rq.type = Exchangetype::WRITE;        
            runtimeWriteQuantities_.push_back(rq);
        }

        // TODO: Reading is not finished yet
        for (const auto &entry : activeParticipantConfig_.readData) {
            PreciceRuntimeQuantity rq;
            rq.precicename = entry.name;
        
            // Validate and convert the precicename to an openCFS result name
            rq.cfsname = convertResultNamesToCFS(rq.precicename);
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
    }

    void PreciceAdapter::setupMeshAndCoordinates() {
        // Get grid from the cfs domain and flatten node coordinates.
        GridCFS* gridCFS = dynamic_cast<GridCFS*>(domain_->GetGrid());
        if (!gridCFS) {
            EXCEPTION("PreciceAdapter: Domain's grid is not of type GridCFS as expected.");
        }

        // Build nodeNumCoordMap_ and cfsNodeNumsVec_.
        ResultHandler* resHandler = domain_->GetResultHandler();
        auto* resultContextsPtr = resHandler->GetResultContexts();
        if (resultContextsPtr->empty()) {
            EXCEPTION("No result contexts available");
        } else {
            bool foundNodeList = false;
            // Loop over all result contexts as long as we find a node result
            for (auto entry = resultContextsPtr->begin(); entry != resultContextsPtr->end(); ++entry) {
                shared_ptr<ResultHandler::ResultContext> resultContext = entry->second;
                BaseResult & actResult = *(resultContext->result);
                // Check if the current result context's entity list is of type NODE_LIST.
                if (actResult.GetEntityList()->GetType() == EntityList::NODE_LIST) {
                    // Process the node-based result.
                    shared_ptr<EntityList> list = actResult.GetEntityList();
                    EntityIterator it = list->GetIterator();
                    for (it.Begin(); !it.IsEnd(); it++) {
                        int nodenum = it.GetNode();
                        cfsNodeNumsVec_.push_back(nodenum);
                        Vector<double> coord;
                        gridCFS->GetNodeCoordinate(coord, nodenum, true);
                        nodeNumCoordMap_[nodenum] = coord;
                    }
                    foundNodeList = true;
                    break; // Stop processing further result contexts.
                }
            }
            if (!foundNodeList) {
                EXCEPTION("No NODE_LIST result context found");
            }
        }

        // Flatten node coordinates.
        int cfsGridDim = gridCFS->GetDim();
        flatCoords_.resize(cfsNodeNumsVec_.size() * cfsGridDim);
        for (size_t i = 0; i < cfsNodeNumsVec_.size(); ++i) {
            Vector<double> point = nodeNumCoordMap_[cfsNodeNumsVec_[i]];
            for (int j = 0; j < cfsGridDim; ++j) {
                flatCoords_[i * cfsGridDim + j] = point[j];
            }
        }
        
    }


    void PreciceAdapter::initialize(Domain *domain) {
#ifdef USE_PRECICE
    domain_ = domain;

    // Step 1: Retrieve configuration parameters from the openCFS ParamNode.
    readConfigurationParameters();

    // Step 2: Parse the precice-config.xml file 
    readPreciceConfiguration();

    // Step 3: Extract and flatten node coordinates.
    setupMeshAndCoordinates();
    
    // Step 4: Create the PreCICE participant.
    createPreciceParticipant();

    // Step 5: Initialize the PreCICE participant.
    participant_->initialize();
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



    void PreciceAdapter::RegisterTimeStepReadData(){
        // loop over the required results that we need - defined in the precice config
        for (auto &quantity : runtimeReadQuantities_) {
            LOG_DBG(preciceAdapter) << "Participant " << activeParticipantConfig_.name
                                    << " reads data: " << quantity.cfsname
                                    << " on mesh: " << quantity.meshName;

            // we need to resize the data vector
            quantity.data.resize(cfsNodeNumsVec_.size() * quantity.quantitydim);
            participant_->readData(quantity.meshName, quantity.precicename, preciceNodeNumsVec_,
                                dynamic_cast<TransientDriver*>(domain_->GetSingleDriver())->GetDeltaT(),
                                quantity.data);
            
            // for easier handling, we also fill the nodeResultMap
            // but keep in mind that the nodes in nodeResultMap are a point cloud and
            // do not have to be real node-values!
            for (std::size_t i = 0; i < preciceNodeNumsVec_.size(); ++i) {
                int cfsnode = cfsNodeNumsVec_[i];
                Vector<double> vals;
                vals.Resize(quantity.quantitydim);
                for(std::size_t k = 0; k < quantity.quantitydim; ++k){
                    vals[k] = quantity.data[i * quantity.quantitydim + k];
                }
                quantity.nodeResultMap[cfsnode] = vals;
            }
        }

    }


    void PreciceAdapter::RegisterTimeStepWriteData(){
        // Get the result contexts from the opencfs result handler.
        ResultHandler * resHandler = domain_->GetResultHandler();
        auto* resultContextsPtr = resHandler->GetResultContexts();

        // loop over the required results that we need - defined in the precice config
        for (auto &quantity : runtimeWriteQuantities_) {
            // loop over the provided results from openCFS
            for (auto& entry : *resultContextsPtr) {   
                // entry.first is a shared_ptr<BaseResult>
                // entry.second is a shared_ptr<ResultContext>
                shared_ptr<BaseResult> baseResult = entry.first;
                shared_ptr<ResultHandler::ResultContext> resultContext = entry.second;

                BaseResult & actResult  = *(resultContext->result);
                std::string cfsResultName = actResult.GetResultInfo()->resultName;
                // check if we have the correct result
                if(cfsResultName == quantity.cfsname){
                    quantity.available = true;
                    // This is more or less a c&p from FeFunction<T>::ExtractResult( shared_ptr<BaseResult> res )
                    shared_ptr<EntityList> list = actResult.GetEntityList();
                    EntityIterator it = list->GetIterator();
                    EntityList::ListType entityListType = actResult.GetEntityList()->GetType();

                    // Get results and location information for precice
                    switch( entityListType ) 
                    {
                        case EntityList::ELEM_LIST:
                        case EntityList::SURF_ELEM_LIST:
                        {
                            // in FieldCoefFunctor<TYPE>::EvalResult this is handled in a separate
                            // way compared to nodal evaluation. The reason for that is that element
                            // values, e.g., mech stress tensor, is a derived quantity and it's evaluated
                            // in the cell/element center
                            // loop over elements
                            // for(it.Begin(); !it.IsEnd(); it++)
                            // {
                            //     const Elem * el = it.GetElem();
                            //     LocPoint lp = Elem::shapes[el->type].midPointCoord;
                            //     LocPointMapped lpm;
                            //     shared_ptr<ElemShapeMap> esm = it.GetGrid()->GetElemShapeMap(el, true);
                            //     lpm.Set(lp, esm, 0.0);
                            //     EXCEPTION("No element value handling via precice yet!");
                                // resultContext->functor->EvalResult(actContext.result);
                                // this->GetVector(tempField, lpm );
                                // loop over dofs
                                // for(UInt iDim = 0; iDim < dim_; iDim++ )
                                //     vec[it.GetPos()*dim_ + iDim] = tempField[iDim];
                                // }
                            // }
                            WARN("PreciceAdapter: Element-based result processing is not yet implemented.");
                            break;
                        }
                        case EntityList::NODE_LIST:
                        {
                            // Determine the dimension for this quantity from the runtime container.
                            int quantityDim = quantity.quantitydim;
                            // quick sanity check
                            if(quantityDim != actResult.GetResultInfo()->GetDimDof()){
                                EXCEPTION("PreciceAdapter: the requested quantity should have dimension "<<
                                         quantityDim<<" but openCFS provided a result with dimension "<<
                                         actResult.GetResultInfo()->GetDimDof());
                            }
                            
                            // Retrieve the result vector from actResult.
                            Vector<double>* res_vec = dynamic_cast<Vector<Double>*>(actResult.GetSingleVector());
                            if (!res_vec) {
                                EXCEPTION("PreciceAdapter: Unable to retrieve node result vector for opencfs result " << quantity.cfsname);
                            }

                            int pos = 0;
                            // Loop over nodes.
                            for(it.Begin(); !it.IsEnd(); it++)
                            {
                                // we could either go via the FeSpace and get the equation number
                                // for the specific entity but since we know how the actSol
                                // was generated we can skip this re-organization step
                                // Have a look at void FeFunction<T>::ExtractResult( shared_ptr<BaseResult> res )
                                // there the connection to solution (res_vec) and eqnNumbers is set
                                int node = it.GetNode();
                                // Create a vector of length 'dim' for the result of this node.
                                Vector<double> resultForNode(quantityDim);
                                for (int k = 0; k < quantityDim; ++k)
                                {
                                    resultForNode[k] = (*res_vec)[pos * quantityDim + k];
                                    LOG_DBG(preciceAdapter)
                                        << "Reading result for node " << node
                                        << " component " << k
                                        << ": " << resultForNode[k];
                                }
                                quantity.nodeResultMap[node] = resultForNode;
                                ++pos;
                            }


                            // reate a flat vector from the node results using the ordering from cfsNodeNumsVec_.
                            quantity.data.resize(cfsNodeNumsVec_.size() * quantityDim);
                            for (std::size_t i = 0; i < cfsNodeNumsVec_.size(); ++i)
                            {
                                int node = cfsNodeNumsVec_[i]; // Get node number in the registered order.
                                Vector<double> nodeResult = quantity.nodeResultMap[node];
                                for (int k = 0; k < quantityDim; ++k)
                                {
                                    quantity.data[i * quantityDim + k] = nodeResult[k];
                                }
                            }
                            break;
                        }
                        default:
                            EXCEPTION("PreciceAdapter::RegisterTimeStep: Unhandled entity list type.");
                            break;
                    }
                }
            }
            if(!quantity.available){
                EXCEPTION("PreciceAdapter: result "<<quantity.cfsname<<" was requested from openCFS but not provided");
            }

            LOG_DBG(preciceAdapter) << "Participant " << activeParticipantConfig_.name
                                    << " writes data: " << quantity.precicename
                                    << " on mesh: " << quantity.meshName;
            participant_->writeData( quantity.meshName, quantity.precicename,
                                    preciceNodeNumsVec_, quantity.data);
        }

        double cfs_dt = dynamic_cast<TransientDriver*>(domain_->GetSingleDriver())->GetDeltaT();
        double precice_dt = this->participant_->getMaxTimeStepSize();
        double dt = cfs_dt <= precice_dt? cfs_dt : precice_dt;
        this->participant_->advance(dt);
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
