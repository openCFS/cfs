// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef WVT_HEADER_2013
#define WVT_HEADER_2013

#include <string>

#include "General/Environment.hh"
#include "DataInOut/SimInput.hh"

namespace CFSTool
{

  class WVT 
  {
  public:

    WVT( const std::string& lateral_mode_file,
         const std::string& coriolis_mode_file,
         const std::string& mean_flow_file,
         const std::string& outFile,
         const PtrParamNode& param,
         const PtrParamNode& info);

    ~WVT() {};

  private:

    const PtrParamNode& param_;
    const PtrParamNode& info_;

    bool writeOutputFile_;
    
    StdVector< shared_ptr<SimInput> > inputs_;

  };

}

#endif
