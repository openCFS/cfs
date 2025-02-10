
#include "IPreciceAdapter.hh"
#include "PreciceAdapter.hh"
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

    void PreciceAdapter::initialize(Domain *domain)
    {
        // Steps
        // 1. Retrieve necessary configuration parameters from the parameter node.
        // 2. Create the preCICE participant using the retrieved configuration.
        // 3. Flatten the cfs grid’s node coordinates into a flat vector format required by preCICE.
        // 4. Generate a sequential list of vertex IDs. == THIS MUST BE ADAPTED FOR REGION HANDLING
        // 5. Register the mesh (with node coordinates and vertex IDs) with preCICE.
        // 6. Determine and reserve the space for the exchange quantity.
        // 7. Initialize the preCICE participant.
#ifdef USE_PRECICE
        if (!paramNode_)
        {
            EXCEPTION("PreciceAdapter: Parameter node is not set.");
        }

        // 1. Retrieve necessary configuration parameters from the parameter node.
        paramNode_->Get("fileFormats/preciceCoupling/participantName")->GetValue("name", participantName_);
        paramNode_->Get("fileFormats/preciceCoupling/participantMeshName")->GetValue("name", participantMeshName_);
        paramNode_->Get("fileFormats/preciceCoupling/exchangeQuantity")->GetValue("name", participantExchangeQuantityName_);
        paramNode_->Get("fileFormats/preciceCoupling/precice_configFile")->GetValue("file", configFileName_);

        LOG_DBG(preciceAdapter) << "Initializing PreCICE participant: " << participantName_
                                << " with config file: " << configFileName_
                                << " at rank: " << rank_ << " out of size: " << size_;

        LOG_DBG(preciceAdapter) << "Provided mesh for participant " << participantName_
                                << " is: " << participantMeshName_;


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


        domain_ = domain;
        if (domain_->GetGridMap().size() > 1)
        {
            EXCEPTION("PreciceAdapter: works only with one grid per simulation - until now");
        }

        // add a check that we dont have a multidriver and that it's transient
        if( !dynamic_cast<TransientDriver*>(domain_->GetSingleDriver()) ){
            EXCEPTION("PreciceAdapter: this only works when using a transient simulation")
        }

        // 2. Create the PreCICE participant using the retrieved configuration.
        participant_ = std::make_unique<precice::Participant>(
            participantName_,
            configFileName_,
            rank_,
            size_);


        // Retrieve the mesh dimension from preCICE using the mesh name.
        int meshDim = participant_->getMeshDimensions(participantMeshName_);

        /* 
        Currently we will not be able to couple higher order results with precice because there is no
        infrastructure that maps virtual nodes (imagine a first order grid and second order Lagrangian shape
        functions, then there are also dof's on edges and faces and interior - depending on the order). But
        the ResultHandler knows that we have a "definedOn" variable which gets set in the actual PDE and let's assume
        it says NODES. Then the ResultHandler evaluates the FeFunction (which includes the virtual nodes from higher
        order polynomials) at exactly the wanted node (now it's a vertex node). Actually not really because we get
        the node number where we want to get the result at and simply take the value there - we are not including
        the other polynomial nodes at this location - which on the other hand might be correct if Lagrangian elements
        are used, then the shape function of other (polynomial) nodes are zero at the evaluation node anyway...
        If we ever wanted to include higher order coupling, we would need to be able to write higher order results.
        So, being able to write higher order results would be the first step towards coupling higher order results
        with precice.
        In the current version, when evaluating the results, the final vector already includes DBC's, so we don't
        need to worry about that at this place.
        */
        /* Alternative end */

        // Downcast the grid pointer to GridCFS to access the no-argument version of GetNodeCoordinates().
        GridCFS* gridCFS = dynamic_cast<GridCFS*>(domain_->GetGrid());
        if (!gridCFS) {
            EXCEPTION("PreciceAdapter: Domain's grid is not of type GridCFS as expected.");
        }

        // node number (1 based, does not have to start at 1)
        ResultHandler * resHandler = domain_->GetResultHandler();
        auto* resultContextsPtr = resHandler->GetResultContexts();
        
        // Make sure the exchanged result via precice is defined in the openCFS xml file
        // this means, we need to check if it there is a Result for that quantity
        bool found = false;
        for (auto& entry : *resultContextsPtr) {
            shared_ptr<ResultHandler::ResultContext> resultContext = entry.second;

            BaseResult & actResult  = *(resultContext->result);
            shared_ptr<EntityList> list = actResult.GetEntityList();
            EntityIterator it = list->GetIterator();

            LOG_DBG(preciceAdapter) << "Registering result '" << actResult.GetResultInfo()->resultName
                                << "' on '" << actResult.GetEntityList()->GetName()<< "' in precice adapter";
            if(cfsExchangeQuantityName_ == actResult.GetResultInfo()->resultName){
                found = true;
            }

            EntityList::ListType entityListType = actResult.GetEntityList()->GetType();
            switch( entityListType ) 
            {
                case EntityList::ELEM_LIST:
                case EntityList::SURF_ELEM_LIST:
                {
                    WARN("PreciceAdapter: result exchange via elements not yet possible");
                    break;
                }
                case EntityList::NODE_LIST:
                {
                    for(it.Begin(); !it.IsEnd(); it++)
                    {
                        int nodenum = it.GetNode();
                        nodenumsvec_.push_back(nodenum);
                        Vector<double> coord;
                        gridCFS->GetNodeCoordinate(coord,nodenum,true);
                        nodeNumCoordMap_[nodenum] = coord;
                    }
                    break;
                }
                default:
                EXCEPTION("PreciceAdapter: Unknown entityListType "<<entityListType);
            }

        }
        if(!found){
            EXCEPTION("PreciceAdapter: I cannot find the result [precicename:"<<participantExchangeQuantityName_<<", cfsname:"<<cfsExchangeQuantityName_<<
                        "] in the available results."<< " Please check the precice and openCFS xml files.");
        }

        // 3. Create a flat vector with the expected size.
        std::vector<double> flatCoords(nodenumsvec_.size() * meshDim);
        // Flatten the node coordinate container into a flat array required by preCICE.
        // The flatCoords vector will have layout: [x0, y0, z0, x1, y1, z1, ...] (for 3D)
        // or [x0, y0, x1, y1, ...] (for 2D), depending on meshDim.
        for (int i = 0; i < nodenumsvec_.size(); ++i)
        {
            // ATTENTION 1 based!!!
            Vector<double> point = nodeNumCoordMap_[i+1];
            for (int j = 0; j < meshDim; ++j)
            {
                flatCoords[i * meshDim + j] = point[j];
            }
        }

        // take care, precicenodenumsvec_ gets adapted by setMeshVertices!!!
        precicenodenumsvec_ = nodenumsvec_;
        participant_->setMeshVertices(participantMeshName_, flatCoords, precicenodenumsvec_);


        // Retrieve the exchange quantity dimension for the provided mesh and exchange data name.
        // This dimension determines how many components each exchange data entry has.
        int exchangeQuantityDim = participant_->getDataDimensions(participantMeshName_, participantExchangeQuantityName_);
        exchangeQuantity_.resize(nodenumsvec_.size() * exchangeQuantityDim);


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
