#ifndef FILE_CFS_DUMMYPRECICEADAPTER_HH
#define FILE_CFS_DUMMYPRECICEADAPTER_HH

#include <limits>

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

        void Initialize(Domain* domain, SinglePDE* pde) override {
            // No-op
        }

        void RegisterTimeStepWriteData() override {
            // No-op
        }

        void RegisterTimeStepReadData(bool atWindowStart = false) override {
            // No-op
        }

        void Finalize() override {
            // No-op
        }

        Vector<Double> GetElemResult(SolutionType solType, int elemNum) override{
            // No-op
            return Vector<Double>();
        }
        Vector<Double> GetNodeResult(SolutionType solType, int nodeNum) override{
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


        virtual Double GetMaxTimeStepSize() override {
            // Dummy: never actually invoked (the dummy adapter routes to the plain
            // TransientDriver, not TransientDriverPrecice). Return "no constraint" so
            // that, if ever called, dt = min(cfs_dt, precice_dt) is not shrunk to 0.
            return std::numeric_limits<Double>::max();
        }
        virtual void Advance(Double dt) override {
            return;
        }

        virtual bool IsPreciceDummy() override {
            return true;
        }

        virtual UInt GetSequenceStep() override {
            return 0;
        }

        virtual void UpdateDomain(Domain* domain) override{
            return;
        }

        virtual void MarkReadResultsUpdated() override {
            // No-op
        }
    };
}

#endif // FILE_CFS_DUMMYPRECICEADAPTER_HH
