// PreciceAdapterFactory.hh

#ifndef FILE_CFS_PRECICEADAPTERFACTORY_HH
#define FILE_CFS_PRECICEADAPTERFACTORY_HH

#include "IPreciceAdapter.hh"
#include "PreciceAdapter.hh"
#include "DummyPreciceAdapter.hh"
#include <memory>
#include "DataInOut/ParamHandling/ParamNode.hh"

namespace CoupledField
{
  /**
   * @brief Factory function to create the appropriate PreciceAdapter.
   *
   * @param paramNode Shared pointer to ParamNode containing configuration.
   * @return Unique pointer to IPreciceAdapter.
   */
  inline std::unique_ptr<IPreciceAdapter> CreatePreciceAdapter(shared_ptr<ParamNode> paramNode) {
#ifdef USE_PRECICE
    // Determine at runtime whether to use PreCICE based on configuration:
    // only build the real adapter when a <fileFormats><preciceCoupling> block is
    // present. Otherwise fall back to the dummy so standalone (uncoupled) transient
    // runs work even in a preCICE-enabled build.
    if (paramNode) {
      auto ff = paramNode->Get("fileFormats", ParamNode::PASS);
      if (ff && ff->Has("preciceCoupling"))
          return std::make_unique<PreciceAdapter>(paramNode);
    }
    return std::make_unique<DummyPreciceAdapter>();
#else
    // PreCICE is not compiled in; return DummyPreciceAdapter
    return std::make_unique<DummyPreciceAdapter>();
#endif
  }
}

#endif // FILE_CFS_PRECICEADAPTERFACTORY_HH
