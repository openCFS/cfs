#ifndef FILE_CFS_PRECICEADAPTER_HH
#define FILE_CFS_PRECICEADAPTER_HH

#include "IPreciceAdapter.hh"
#include "PreciceConfigReader.hh"
#include <string>
#include <vector>
#include "MatVec/Vector.hh"
#include <boost/shared_ptr.hpp>
#include "def_use_precice.hh"
#ifdef USE_PRECICE
#include <precice.hpp>
#endif

namespace CoupledField
{
        class ParamNode;
        class Domain;
        class BaseSolveStep;
        class GridCFS;

        /**
         * The PreciceAdapter class manages the coupling between openCFS and preCICE.
         *
         * This adapter retrieves configuration parameters, creates the preCICE participant,
         * registers mesh information and exchange quantities, and handles data mapping.
         */
        class PreciceAdapter : public IPreciceAdapter
        {
        public:
                PreciceAdapter(boost::shared_ptr<ParamNode> paramNode);
                ~PreciceAdapter() override;

                // Main adapter methods
                void initialize(Domain *domain) override;
                void RegisterSolveStep(BaseSolveStep *solveStep) override;
                void RegisterTimeStep() override;
                void finalize() override;



        private:

        // --- Helper functions for initialize() ---
        /**
         * Reads configuration parameters (participant name, mesh name, etc.) from paramNode_.
         */
        void readConfigurationParameters();

        /**
         * Reads and processes the preCICE configuration file (e.g. receive-mesh, provide-mesh,
         *        read-data, and write-data entries).
         */
        void readPreciceConfiguration();

        /**
         * Creates the preCICE participant using the configuration parameters.
         */
        void createPreciceParticipant();

        /**
         * Extracts the mesh node coordinates from the domain and registers the mesh with preCICE.
         */
        void setupMeshAndCoordinates();

        /**
         * Reserves the exchange quantity vector based on the data dimensions and number of nodes.
         */
        void reserveExchangeQuantities();
        // --- End of helper functions ---


        // --- Data members ---
        boost::shared_ptr<ParamNode> paramNode_;
#ifdef USE_PRECICE
        std::unique_ptr<precice::Participant> participant_; ///< preCICE participant instance
#endif
        std::string configFileName_;
        std::string participantName_;
        std::string participantMeshName_;
        std::string participantExchangeQuantityName_;
        std::string cfsExchangeQuantityName_;

        // Mapping from node numbers to coordinates
        std::map<unsigned int, Vector<double>> nodeNumCoordMap_;
        // Mapping for results (if needed)
        std::map<int, Vector<double>> nodeResultMap_;
        std::vector<double> flatResults_;
        std::vector<int> nodenumsvec_;
        std::vector<int> precicenodenumsvec_;
        std::vector<double> exchangeQuantity_;
        ParticipantConfig activeParticipantConfig_;

        int rank_;
        int size_;

        Domain *domain_;
        BaseSolveStep *solveStep_;
        // --- End of data members ---

        // Disable copy and assignment.
        PreciceAdapter(const PreciceAdapter &) = delete;
        PreciceAdapter &operator=(const PreciceAdapter &) = delete;

        };
} // end of namespace
#endif
