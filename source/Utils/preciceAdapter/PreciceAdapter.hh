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
                void RegisterTimeStepWriteData() override;
                void RegisterTimeStepReadData() override;
                void finalize() override;

                enum Exchangetype {READ=0, WRITE=1};

                /**
                 * Contains runtime data for an exchange quantity.
                 *
                 * This container links the configuration (the quantity name and the mesh on which it is defined)
                 * with the runtime data: the node numbers at which the quantity is defined,
                 * and the actual exchanged data.
                 */
                struct PreciceRuntimeQuantity {
                        std::string precicename;                // Quantity name (e.g., "Temperature")
                        std::string cfsname;                // Quantity name in cfs (e.g., "heatTemperature")
                        std::string meshName;            // Mesh name on which this quantity is defined.
                        //int griddim;                    // spatial dimension of the grid
                        int quantitydim;                // dimension of the quantity (scalar, vector)
                        //GridCFS* meshPtr;                // Pointer to the corresponding mesh object.
                        //std::vector<int> nodeNumbers;    // Node numbers where the quantity is defined.
                        std::vector<double> data;        // The actual data (flattened) for exchange.
                        std::map<int, Vector<double>> nodeResultMap; // different representation of the data. key is the cfs node number
                        bool available;                  // flag if the requested quantity is ready
                        Exchangetype type;               // flag if this is read or write data
                };

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
        
        /**
         * Convert the result name specified in precice's config to openCFS result name
         */
        std::string convertResultNamesToCFS(std::string precicename);
        // --- End of helper functions ---

        


        // --- Data members ---
        boost::shared_ptr<ParamNode> paramNode_;
#ifdef USE_PRECICE
        std::unique_ptr<precice::Participant> participant_; ///< preCICE participant instance
#endif
        std::string configFileName_;
        std::string participantName_;
        std::string participantMeshName_;

        // Mapping from node numbers to coordinates
        std::map<unsigned int, Vector<double>> nodeNumCoordMap_;
        std::vector<int> cfsNodeNumsVec_;
        std::vector<double> flatCoords_;
        std::vector<int> preciceNodeNumsVec_;
        ParticipantConfig activeParticipantConfig_;

        // runtime containers
        std::vector<PreciceRuntimeQuantity> runtimeReadQuantities_;
        std::vector<PreciceRuntimeQuantity> runtimeWriteQuantities_;

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
