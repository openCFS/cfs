// source/Utils/preciceAdapter/NullPreciceAdapter.hh

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

        void initialize() override {
            // No-op
        }

        void finalize() override {
            // No-op
        }
    };
}

#endif // FILE_CFS_DUMMYPRECICEADAPTER_HH
