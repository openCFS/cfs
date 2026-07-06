#ifndef FILE_CFS_IPRECICEADAPTER_HH
#define FILE_CFS_IPRECICEADAPTER_HH

#include "MatVec/Vector.hh"

// Forward declaration
namespace CoupledField {
    class ParamNode;
    class Domain;
    class BaseSolveStep;
    class SinglePDE;

}

namespace CoupledField
{
    class IPreciceAdapter;

    //! Define global instance of this class
    extern IPreciceAdapter* gPreciceAdapter;


    /**
     * Abstract interface for PreciceAdapter.
     *
     * Defines the contract for PreCICE interaction. Concrete implementations
     * can either use PreCICE or perform no operations.
     */
    class IPreciceAdapter {
    public:
        virtual ~IPreciceAdapter() = default;

        /**
         * Initializes the PreCICE participant.
         */
        virtual void initialize(Domain* domain, SinglePDE* pde) = 0;

        /**
         * Writes data of the finished time step
         */
        virtual void RegisterTimeStepWriteData() = 0;

        /**
         * Reads data needed for the current time step.
         *
         * atWindowStart = false (default): read at relativeReadTime = deltaT (end of the
         * current time window) - the normal per-iteration read before solving.
         * atWindowStart = true: read at relativeReadTime = 0. Meant for AFTER Advance()
         * completed a CONVERGED window: relative time 0 of the new window is exactly the
         * converged end of the completed one. The per-iteration reads leave the buffers
         * holding the last-used (pre-convergence) iterate; for passive channels (data
         * without a convergence measure, e.g. DisplacementVolume) that iterate is not
         * reproducible across runs - the converged sample is. See TransientDriverPrecice.
         */
        virtual void RegisterTimeStepReadData(bool atWindowStart = false) = 0;


        virtual bool IsCouplingOngoing() = 0;
        virtual bool RequiresWritingCheckpoint() = 0;
        virtual bool RequiresReadingCheckpoint() = 0;
        virtual Double GetMaxTimeStepSize() = 0;
        virtual void Advance(Double dt) = 0;


        virtual Vector<Double> GetElemResult(SolutionType solType, int elemNum) = 0;
        virtual Vector<Double> GetNodeResult(SolutionType solType, int nodeNum) = 0;

        virtual bool IsPreciceDummy() = 0;
        virtual UInt GetSequenceStep() = 0;

        virtual void UpdateDomain(Domain* domain) = 0;

        /**
         * Marks all read results (filled by RegisterTimeStepReadData) as updated
         * in the ResultHandler so they can be written to output even without a functor.
         * Must be called after ResultHandler::BeginStep() and before FinishStep().
         */
        virtual void MarkReadResultsUpdated() = 0;

        /**
         * Finalizes the PreCICE participant.
         *
         * Called automatically upon destruction.
         */
        virtual void finalize() = 0;
    };
}

#endif // FILE_CFS_IPRECICEADAPTER_HH