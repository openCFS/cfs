
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
          participantExchangeQuantityName_(""),               
          cfsExchangeQuantityName_(""),
          flatResults_(0.0),
          nodenumsvec_(0),
          exchangeQuantity_(0.0)
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
        paramNode_->Get("fileFormats/preciceCoupling/exchangeQuantity")->GetValue("name", participantExchangeQuantityName_);
        paramNode_->Get("fileFormats/preciceCoupling/precice_configFile")->GetValue("file", configFileName_);


        // Validate and convert the participantExchangeQuantityName_ to an openCFS result name
        std::string convertedName;
        if (participantExchangeQuantityName_ == "Temperature") {
            convertedName = "heatTemperature";
        } else if (participantExchangeQuantityName_ == "Displacement") {
            convertedName = "mechDisplacement";
        } else {
            EXCEPTION("Invalid participantExchangeQuantityName_: " << participantExchangeQuantityName_
                    << ". Currently the adapter works for one of [\"Temperature\", \"Displacement\"]");
        }
        cfsExchangeQuantityName_ = convertedName;

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
        participant_ = std::make_unique<precice::Participant>(
            participantName_, configFileName_, rank_, size_);
    }

    void PreciceAdapter::setupMeshAndCoordinates() {
        // Retrieve mesh dimension from participant.
        int meshDim = participant_->getMeshDimensions(participantMeshName_);
        // Get grid from the domain and flatten node coordinates.
        GridCFS* gridCFS = dynamic_cast<GridCFS*>(domain_->GetGrid());
        if (!gridCFS) {
            EXCEPTION("PreciceAdapter: Domain's grid is not of type GridCFS as expected.");
        }

        // Build nodeNumCoordMap_ and nodenumsvec_.
        ResultHandler* resHandler = domain_->GetResultHandler();
        auto* resultContextsPtr = resHandler->GetResultContexts();
        for (auto &entry : *resultContextsPtr) {
            shared_ptr<ResultHandler::ResultContext> resultContext = entry.second;
            BaseResult & actResult = *(resultContext->result);
            shared_ptr<EntityList> list = actResult.GetEntityList();
            EntityIterator it = list->GetIterator();
            if (actResult.GetEntityList()->GetType() == EntityList::NODE_LIST) {
                for (it.Begin(); !it.IsEnd(); it++) {
                    int nodenum = it.GetNode();
                    nodenumsvec_.push_back(nodenum);
                    Vector<double> coord;
                    gridCFS->GetNodeCoordinate(coord, nodenum, true);
                    nodeNumCoordMap_[nodenum] = coord;
                }
            }
        }
        // Flatten node coordinates.
        std::vector<double> flatCoords(nodenumsvec_.size() * meshDim);
        for (size_t i = 0; i < nodenumsvec_.size(); ++i) {
            Vector<double> point = nodeNumCoordMap_[nodenumsvec_[i]];
            for (int j = 0; j < meshDim; ++j) {
                flatCoords[i * meshDim + j] = point[j];
            }
        }
        precicenodenumsvec_ = nodenumsvec_; // if needed.
        participant_->setMeshVertices(participantMeshName_, flatCoords, precicenodenumsvec_);
    }

    void PreciceAdapter::reserveExchangeQuantities() {
        // Reserve space for exchange data based on the dimension.
        int exchangeQuantityDim = participant_->getDataDimensions(participantMeshName_, participantExchangeQuantityName_);
        exchangeQuantity_.resize(nodenumsvec_.size() * exchangeQuantityDim);
    }




    void PreciceAdapter::initialize(Domain *domain) {
#ifdef USE_PRECICE
    domain_ = domain;

    // Step 1: Retrieve configuration parameters from the ParamNode.
    readConfigurationParameters();

    // Step 2: Parse the precice-config.xml file 
    readPreciceConfiguration();

    // Step 3: Create the PreCICE participant.
    createPreciceParticipant();

    // Step 4: Extract and flatten node coordinates.
    setupMeshAndCoordinates();

    // Step 5: Reserve space for exchange quantities.
    reserveExchangeQuantities();

    // Step 6: Initialize the PreCICE participant.
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

    void PreciceAdapter::RegisterTimeStep(){
        ResultHandler * resHandler = domain_->GetResultHandler();

        auto* resultContextsPtr = resHandler->GetResultContexts();
        Vector<Double>*  res_vec = NULL;
        for (auto& entry : *resultContextsPtr) {
            // entry.first is a shared_ptr<BaseResult>
            // entry.second is a shared_ptr<ResultContext>
            shared_ptr<BaseResult> baseResult = entry.first;
            shared_ptr<ResultHandler::ResultContext> resultContext = entry.second;

            BaseResult & actResult  = *(resultContext->result);
            // This is more or less a c&p from FeFunction<T>::ExtractResult( shared_ptr<BaseResult> res )
            shared_ptr<EntityList> list = actResult.GetEntityList();
            EntityIterator it = list->GetIterator();

            // Get results and location information for precice
            Vector<double> tempField;
            EntityList::ListType entityListType = actResult.GetEntityList()->GetType();
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
                    break;
                }
                case EntityList::NODE_LIST:
                {
                    
                    res_vec = dynamic_cast<Vector<Double>*>(actResult.GetSingleVector());
                    
                    int pos = 0;
                    int resdim = participant_->getDataDimensions(participantMeshName_, participantExchangeQuantityName_);
                    for(it.Begin(); !it.IsEnd(); it++)
                    {
                        // we could either go via the FeSpace and get the equation number
                        // for the specific entity but since we know how the actSol
                        // was generated we can skip this re-organization step

                        // Have a look at void FeFunction<T>::ExtractResult( shared_ptr<BaseResult> res )
                        // there the connection to solution (res_vec) and eqnNumbers is set

                        int node = it.GetNode();
                        // Create a vector of length 'dim' for the result of this node.
                        Vector<double> resultForNode(resdim);
                        for (int k = 0; k < resdim; ++k)
                        {
                            resultForNode[k] = (*res_vec)[pos * resdim + k];
                            LOG_DBG(preciceAdapter)
                                << "Reading result for node " << node
                                << " component " << k
                                << ": " << resultForNode[k];
                        }
                        nodeResultMap_[node] = resultForNode;
                        ++pos;
                    }
                    
                    // reate a flat vector from the node results using the ordering from nodenumsvec_.
                    flatResults_.resize(nodenumsvec_.size() * resdim);
                    for (std::size_t i = 0; i < nodenumsvec_.size(); ++i)
                    {
                        int node = nodenumsvec_[i]; // Get node number in the registered order.
                        Vector<double> nodeResult = nodeResultMap_[node];
                        for (int k = 0; k < resdim; ++k)
                        {
                            flatResults_[i * resdim + k] = nodeResult[k];
                        }
                    }
                    break;
                }
            }
        }

        participant_->writeData(participantMeshName_, participantExchangeQuantityName_,
                                precicenodenumsvec_, flatResults_);
        // participant_->readData(participantMeshName_, participantExchangeQuantityName_, precicenodenumsvec_,
        //                        dynamic_cast<TransientDriver*>(domain_->GetSingleDriver())->GetDeltaT(),
        //                        flatResults_);

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
