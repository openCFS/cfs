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

        void initialize(Domain* domain, SinglePDE* pde) override {
            // No-op
        }

        void RegisterSolveStep(BaseSolveStep* solveStep) override {
            // No-op
        }

        void RegisterTimeStepWriteData() override {
            // No-op
        }

        void RegisterTimeStepReadData() override {
            // No-op
        }

        void finalize() override {
            // No-op
        }

        Vector<Double> GetElemResult(SolutionType solType, int elemNum) override{
            // No-op
            return Vector<Double>();
        }
        


        virtual bool IsCouplingOngoing() override {
            return true;
        }
        virtual bool RequiresWritingCheckpoint() override {
            return true;
        }
        virtual bool RequiresReadingCheckpoint()  override {
            return false;
        }

        // Define a non-inline dummy function to force an object file to be generated.
        // void DummyPreciceAdapter_ForceBuild();
    };
}

#endif // FILE_CFS_DUMMYPRECICEADAPTER_HH
