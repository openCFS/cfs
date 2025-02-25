#ifndef FILE_CFS_PRECICEADAPTER_HH
#define FILE_CFS_PRECICEADAPTER_HH

#include "IPreciceAdapter.hh"
#include "PreciceConfigReader.hh"

#include <string>
#include <vector>
#include <boost/shared_ptr.hpp>
#include "def_use_precice.hh"
#ifdef USE_PRECICE
#include <precice.hpp>
#endif

namespace CoupledField
{
        class ParamNode;
        class Domain;
        class StdSolveStep;
        class GridCFS;
        class SinglePDE;
        class ResultBase;
        
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
                void initialize(Domain *domain, SinglePDE* pde) override;
                void RegisterSolveStep(BaseSolveStep *solveStep) override;
                void RegisterTimeStepWriteData() override;
                void RegisterTimeStepReadData() override;
                void finalize() override;
                Vector<Double> GetElemResult(SolutionType solType, int elemNum) override;
                Vector<Double> GetNodeResult(SolutionType solType, int nodeNum) override;

                enum Exchangetype {READ=0, WRITE=1};

                virtual bool IsCouplingOngoing() override{ return this->participant_->isCouplingOngoing();}
                virtual bool RequiresWritingCheckpoint() override{ return this->participant_->requiresWritingCheckpoint();}
                virtual bool RequiresReadingCheckpoint()  override{ return this->participant_->requiresReadingCheckpoint();}
                virtual Double GetMaxTimeStepSize() override{ return this->participant_->getMaxTimeStepSize();}
                virtual void Advance(Double dt) override{this->participant_->advance(dt);}

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
        

        void RegisterSinglePDE(SinglePDE* pde);

        void RegisterExternalResults();


        /**
         * Convert the result name specified in precice's config to openCFS result name
         */
        std::tuple<std::string, SolutionType> convertResultNamesToCFS(const std::string& precicename);
        // --- End of helper functions ---

        


        // --- Data members ---
        bool isInit_;
        boost::shared_ptr<ParamNode> paramNode_;
#ifdef USE_PRECICE
        std::unique_ptr<precice::Participant> participant_; ///< preCICE participant instance
#endif
        std::string configFileName_;
        std::string participantName_;
        std::string participantMeshName_;
        std::string participantElemMeshName_;
        int sequenceStep_;

        // Mapping from node numbers to coordinates
        std::map<unsigned int, Vector<double>> nodeNumCoordMap_;
        std::vector<int> cfsNodeNumsVec_;
        std::vector<double> flatCoords_;
        std::vector<int> preciceNodeNumsVec_;

        // Mapping from element numbers to midpoint coordinates
        std::map<unsigned int, Vector<double>> elemNumCoordMap_;
        std::vector<int> cfsElemNumsVec_;
        std::vector<double> flatElemCoords_;
        std::vector<int> preciceElemNumsVec_;
        bool needElemMesh_;

        ParticipantConfig activeParticipantConfig_;

        // runtime containers
        std::vector<std::unique_ptr<ResultBase>> runtimeReadResults_;
        std::vector<std::unique_ptr<ResultBase>> runtimeWriteResults_;

        int rank_;
        int size_;

        Domain *domain_;
        StdSolveStep *solveStep_;
        SinglePDE *singlePDE_;
        // --- End of data members ---

        // Disable copy and assignment.
        PreciceAdapter(const PreciceAdapter &) = delete;
        PreciceAdapter &operator=(const PreciceAdapter &) = delete;

        };
} // end of namespace
#endif
