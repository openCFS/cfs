
#include "IPreciceAdapter.hh"
#include "PreciceAdapter.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
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
