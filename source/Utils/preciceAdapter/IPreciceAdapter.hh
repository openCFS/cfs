#ifndef FILE_CFS_IPRECICEADAPTER_HH
#define FILE_CFS_IPRECICEADAPTER_HH


// Forward declaration
namespace CoupledField {
    class ParamNode;
    class Domain;
    class BaseSolveStep;
}

namespace CoupledField
{
    /**
     * @brief Abstract interface for PreciceAdapter.
     *
     * Defines the contract for PreCICE interaction. Concrete implementations
     * can either use PreCICE or perform no operations.
     */
    class IPreciceAdapter {
    public:
        virtual ~IPreciceAdapter() = default;

        /**
         * @brief Initializes the PreCICE participant.
         */
        virtual void initialize(Domain* domain) = 0;

        /**
         * @brief Initializes the PreCICE participant.
         */
        virtual void RegisterSolveStep(BaseSolveStep* solvestep) = 0;

        virtual void RegisterTimeStepWriteData() = 0;

        /**
         * @brief Finalizes the PreCICE participant.
         *
         * Called automatically upon destruction.
         */
        virtual void finalize() = 0;
    };
}

#endif // FILE_CFS_IPRECICEADAPTER_HH