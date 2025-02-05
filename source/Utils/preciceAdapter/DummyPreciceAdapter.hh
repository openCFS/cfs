#ifndef FILE_CFS_DUMMYPRECICEADAPTER_HH
#define FILE_CFS_DUMMYPRECICEADAPTER_HH

#include "IPreciceAdapter.hh"

namespace CoupledField
{
    /**
     * @brief Null implementation of IPreciceAdapter.
     *
     * Performs no operations. Useful when PreCICE is disabled.
     */
    class DummyPreciceAdapter : public IPreciceAdapter {
    public:
        DummyPreciceAdapter() = default;
        ~DummyPreciceAdapter() override = default;

        void initialize(Domain* domain) override {
            // No-op
        }

        void RegisterSolveStep(BaseSolveStep* solveStep) override {
            // No-op
        }

        void finalize() override {
            // No-op
        }

        // Define a non-inline dummy function to force an object file to be generated.
        // void DummyPreciceAdapter_ForceBuild();
    };
}

#endif // FILE_CFS_DUMMYPRECICEADAPTER_HH
