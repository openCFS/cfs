#ifndef FILE_CFS_PRECICEADAPTER_HH
#define FILE_CFS_PRECICEADAPTER_HH

#include "IPreciceAdapter.hh"
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
        class BaseSolveStep;

        class PreciceAdapter : public IPreciceAdapter
        {
        public:
                PreciceAdapter(boost::shared_ptr<ParamNode> paramNode);

                ~PreciceAdapter() override;

                void initialize(Domain *domain) override;
                void RegisterSolveStep(BaseSolveStep *solveStep) override;
                void RegisterTimeStep() override;
                void finalize() override;

        private:
                boost::shared_ptr<ParamNode> paramNode_;
#ifdef USE_PRECICE
                std::unique_ptr<precice::Participant> participant_; ///< PreCICE participant instance
#endif
                std::string configFileName_;
                std::string participantName_;
                std::string participantMeshName_;
                std::string participantExchangeQuantityName_;               
                std::string cfsExchangeQuantityName_;
                std::vector<double> exchangeQuantity_;
                int rank_;
                int size_;

                Domain *domain_;
                BaseSolveStep *solveStep_;

                // Disable copy and assignment
                PreciceAdapter(const PreciceAdapter &) = delete;
                PreciceAdapter &operator=(const PreciceAdapter &) = delete;
        };

} // end of namespace
#endif
