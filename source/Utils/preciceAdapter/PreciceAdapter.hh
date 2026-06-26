#ifndef FILE_CFS_PRECICEADAPTER_HH
#define FILE_CFS_PRECICEADAPTER_HH

#include "IPreciceAdapter.hh"
#include "PreciceConfigReader.hh"

#include <string>
#include <vector>
#include <map>
#include <set>
#include <boost/shared_ptr.hpp>
#include "def_use_precice.hh"
#ifdef USE_PRECICE
#include <precice.hpp>
#endif

namespace CoupledField
{
        class ParamNode;
        class Domain;
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
                void RegisterTimeStepWriteData() override;
                void RegisterTimeStepReadData() override;
                void finalize() override;
                Vector<Double> GetElemResult(SolutionType solType, int elemNum) override;
                Vector<Double> GetNodeResult(SolutionType solType, int nodeNum) override;

                enum Exchangetype {READ=0, WRITE=1};
#ifdef USE_PRECICE
                virtual bool IsCouplingOngoing() override{ 
                         if (this->participant_ == nullptr) {
                            EXCEPTION("Error: Precice participant is null in IsCouplingOngoing()."); 
                         }
                        return this->participant_->isCouplingOngoing();}
                virtual bool RequiresWritingCheckpoint() override{ return this->participant_->requiresWritingCheckpoint();}
                virtual bool RequiresReadingCheckpoint()  override{ return this->participant_->requiresReadingCheckpoint();}
                virtual Double GetMaxTimeStepSize() override{ return this->participant_->getMaxTimeStepSize();}
                virtual void Advance(Double dt) override{this->participant_->advance(dt);}
#endif

                virtual bool IsPreciceDummy() override {
                        return false;
                }

                virtual UInt GetSequenceStep() override{
                        return this->sequenceStep_;
                }
                
                virtual void UpdateDomain(Domain* domain) override { this->domain_ = domain;}

                /**
                 * Marks all read results (already filled by RegisterTimeStepReadData) as updated
                 * in the ResultHandler. Must be called after ResultHandler::BeginStep().
                 */
                virtual void MarkReadResultsUpdated() override;
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

                                struct CouplingMeshData {
                                        std::vector<int> cfsNodeNums;
                                        std::vector<double> flatNodeCoords;
                                        std::vector<int> preciceNodeNums;

                                        std::vector<int> cfsElemNums;
                                        std::vector<double> flatElemCoords;
                                        std::vector<int> preciceElemNums;

                                        bool needsNodeData = false;
                                        bool needsElemData = false;
                                        std::vector<std::string> regionNames;
                                };

                                const CouplingMeshData& getMeshData(const std::string& meshName) const;
                                CouplingMeshData& getMeshData(const std::string& meshName);


        /**
         * Convert the result name specified in precice's config to openCFS result name.
         * @param isWrite true for write-data, false for read-data. Needed for the
         *        characteristic (impedance-matched) coupling, where the openCFS quantity
         *        is role-dependent (write -> outgoing invariant acouCharacteristic,
         *        read -> incoming invariant acouCharacteristicCoupling) so two openCFS
         *        participants can couple symmetrically. For all other data names the
         *        mapping is the same for read and write and isWrite is ignored.
         */
        std::tuple<std::string, SolutionType> convertResultNamesToCFS(const std::string& precicename, bool isWrite);
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
        bool elementMeshUsesSurfaceElems_;

        ParticipantConfig activeParticipantConfig_;
        std::map<std::string, CouplingMeshData> meshDataByName_;
        std::vector<std::string> defaultCouplingRegions_;

        // runtime containers
        std::vector<std::unique_ptr<ResultBase>> runtimeReadResults_;
        std::vector<std::unique_ptr<ResultBase>> runtimeWriteResults_;

        int rank_;
        int size_;

        Domain *domain_;
        //! All coupled PDEs registered with this participant. initialize() is
        //! called once per PDE; the participant may drive several (iteratively
        //! coupled) PDEs that exchange data with the same external participant.
        std::vector<SinglePDE*> pdes_;
        //! Guard so external (read) results are registered with the ResultHandler
        //! only once, even though initialize() runs per coupled PDE.
        bool externalResultsRegistered_;
        // --- End of data members ---

        // Disable copy and assignment.
        PreciceAdapter(const PreciceAdapter &) = delete;
        PreciceAdapter &operator=(const PreciceAdapter &) = delete;

        };
} // end of namespace
#endif
