// PreciceAdapterFactory.hh

#ifndef FILE_CFS_PRECICEADAPTERFACTORY_HH
#define FILE_CFS_PRECICEADAPTERFACTORY_HH

#include "IPreciceAdapter.hh"
#include "PreciceAdapter.hh"
#include "DummyPreciceAdapter.hh"
#include <memory>
#include "DataInOut/ParamHandling/ParamNode.hh"
#include <boost/shared_ptr.hpp>

namespace CoupledField
{
    /**
     * @brief Factory function to create the appropriate PreciceAdapter.
     *
     * @param paramNode Shared pointer to ParamNode containing configuration.
     * @return Unique pointer to IPreciceAdapter.
     */
    inline std::unique_ptr<IPreciceAdapter> CreatePreciceAdapter(boost::shared_ptr<ParamNode> paramNode) {
#ifdef USE_PRECICE
        // Determine at runtime whether to use PreCICE based on configuration
        return std::make_unique<PreciceAdapter>(paramNode);
        
#else
        // PreCICE is not compiled in; return DummyPreciceAdapter
        return std::make_unique<DummyPreciceAdapter>();
#endif
    }
}

#endif // FILE_CFS_PRECICEADAPTERFACTORY_HH
