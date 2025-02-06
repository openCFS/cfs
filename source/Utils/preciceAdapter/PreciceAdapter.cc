
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
          size_(1)
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
        Double dt;
        dt = dynamic_cast<TransientDriver*>(domain_->GetSingleDriver())->GetDeltaT();;

        // 2. Create the PreCICE participant using the retrieved configuration.
        participant_ = std::make_unique<precice::Participant>(
            participantName_,
            configFileName_,
            rank_,
            size_);

        // Retrieve the mesh dimension from preCICE using the mesh name.
        int meshDim = participant_->getMeshDimensions(participantMeshName_);
        // Determine number of vertices
        // TODO: this has to be adjusted to the relevant regions as well
        int vertexSize = domain_->GetGrid()->GetNumNodes();

        /* Alternative start */
        int numNodes = domain_->GetGrid()->GetNumElems(CoupledField::ALL_REGIONS);
        int numElements = domain_->GetGrid()->GetNumNodes(CoupledField::ALL_REGIONS);
        int maxNumEqns = numNodes;
        std::vector<int> entityEqns(maxNumEqns);
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

        // 3. Create a flat vector with the expected size.
        const StdVector<Vector<double>>& coordsContainer = gridCFS->GetNodeCoordinates();

        std::vector<double> flatCoords(vertexSize * meshDim);
        // Flatten the node coordinate container into a flat array required by preCICE.
        // The flatCoords vector will have layout: [x0, y0, z0, x1, y1, z1, ...] (for 3D)
        // or [x0, y0, x1, y1, ...] (for 2D), depending on meshDim.
        for (int i = 0; i < vertexSize; ++i)
        {
            const Vector<double>& point = coordsContainer[i];
            for (int j = 0; j < meshDim; ++j)
            {
                flatCoords[i * meshDim + j] = point[j];
            }
        }

        // asking openCFS for the the vertexIDs is not optimal because in cfs they are based on and connected to elements.
        // for now, let us just hack this
        std::vector<int> vertexIDs(vertexSize);
        std::iota(vertexIDs.begin(), vertexIDs.end(), 0); // fills the vector with 0,1,2,...vertexSize-1
        participant_->setMeshVertices(participantMeshName_, flatCoords, vertexIDs);

        if (coordsContainer.GetSize() != static_cast<size_t>(vertexSize)) {
            EXCEPTION("PreciceAdapter: Mismatch between vertex size (" << vertexSize << ") and coordinate container size (" 
                    << coordsContainer.GetSize() << ")");
        }

        // Retrieve the exchange quantity dimension for the provided mesh and exchange data name.
        // This dimension determines how many components each exchange data entry has.
        int exchangeQuantityDim = participant_->getDataDimensions(participantMeshName_, participantExchangeQuantityName_);
        exchangeQuantity_.resize(vertexSize * exchangeQuantityDim);


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

        // TODO We should probably check this at another place
        // Make sure the exchanged result via precice is defined in the openCFS xml file
        // this means, we need to check if it there is a Result for that quantity
        bool found = false;
        auto* resultContextsPtr = resHandler->GetResultContexts();
        for (auto& entry : *resultContextsPtr) {
            // entry.first is a shared_ptr<BaseResult>
            // entry.second is a shared_ptr<ResultContext>
            shared_ptr<BaseResult> baseResult = entry.first;
            shared_ptr<ResultHandler::ResultContext> resultContext = entry.second;

            BaseResult & actResult  = *(resultContext->result);

            LOG_DBG(preciceAdapter) << "Registering result '"
                             << actResult.GetResultInfo()->resultName
                             << "' on '"
                             << actResult.GetEntityList()->GetName()
                             << "' in precice adapter";
            

            if(cfsExchangeQuantityName_ == actResult.GetResultInfo()->resultName){
                found = true;
                // Get results and nodnumbers (or at least location information)


            }
        }
        if(!found){
            EXCEPTION("PreciceAdapter: I cannot find the result [precicename:"<<participantExchangeQuantityName_<<", cfsname:"<<cfsExchangeQuantityName_<<
                      "] in the available results."<< " Please check the precice and openCFS xml files.");
        }


        //participant_->readData("thermSolver-Mesh", "Temperature", vertexIDs, dt, displacements);

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
