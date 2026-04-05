#include "LogConfigurator.hh"

// in release we have no LogConfigurator.cc - have dummy to prevent build issues
namespace CoupledField 
{
  void LogConfigurator::ParseLogConfFile(const std::string& confFile) { }
} 
