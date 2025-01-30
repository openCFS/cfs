#ifndef FILE_CFS_IPRECICEADAPTER_HH
#define FILE_CFS_IPRECICEADAPTER_HH


// Forward declaration
namespace CoupledField {
    class ParamNode;
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
        virtual void initialize() = 0;


        /**
         * @brief Finalizes the PreCICE participant.
         *
         * Called automatically upon destruction.
         */
        virtual void finalize() = 0;
    };
}

#endif // FILE_CFS_IPRECICEADAPTER_HH