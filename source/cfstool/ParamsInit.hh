#ifndef FILE_PARAMS_2009
#define FILE_PARAMS_2009

#include <string>

namespace CFSTool 
{
  void ParamsInit(int argc, char* argv[],
                  PtrParamNode& param,
                  shared_ptr<LogConfigurator>& logConf);
}

#endif
