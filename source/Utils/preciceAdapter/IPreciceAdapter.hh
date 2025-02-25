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
         * Initializes the PreCICE participant.
         */
        virtual void RegisterSolveStep(BaseSolveStep* solvestep) = 0;

        /**
         * Writes data of the finished time step
         */
        virtual void RegisterTimeStepWriteData() = 0;

        /**
         * Reads data needed for the current time step
         */
        virtual void RegisterTimeStepReadData() = 0;


        virtual bool IsCouplingOngoing() = 0;
        virtual bool RequiresWritingCheckpoint() = 0;
        virtual bool RequiresReadingCheckpoint() = 0;
        virtual Double GetMaxTimeStepSize() = 0;
        virtual void Advance(Double dt) = 0;


        virtual Vector<Double> GetElemResult(SolutionType solType, int elemNum) = 0;
        virtual Vector<Double> GetNodeResult(SolutionType solType, int nodeNum) = 0;

        /**
         * Finalizes the PreCICE participant.
         *
         * Called automatically upon destruction.
         */
        virtual void finalize() = 0;
    };
}

#endif // FILE_CFS_IPRECICEADAPTER_HH