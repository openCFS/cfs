#ifndef CFS_LOG_HH
#define CFS_LOG_HH

#include <string>

#ifndef NDEBUG
  // include files of BOOST's log library
  // disable warnings within boost files
  #if !defined(_MSC_VER)
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wall"
  #endif
  #include <boost/log/common.hpp>
  #include <boost/log/expressions.hpp>
  #include <boost/log/sources/severity_logger.hpp>
  #if !defined(_MSC_VER)
    #pragma GCC diagnostic pop
  #endif
#endif 

namespace CoupledField
{
  const int level_disable_all = -1;
  const int level_default_ = 1000;
  const int level_enable_all = 0;

  const int level_fatal  = 2000;
  const int level_err    = 1600;
  const int level_trace  = 1500;
  const int level_trace2 = 1450;
  const int level_dbg    = 1400;
  const int level_dbg2   = 1300;
  const int level_dbg3   = 1200;
  const int level_warn   = 1200;
  const int level_info   = 1000;

#ifdef NDEBUG
  // in release build the logs will be optimized away
  #define LOG_DBG3(logger) if(true) {} else std::cout
  #define LOG_DBG2(logger) if(true) {} else std::cout
  #define LOG_DBG(logger) if(true) {} else std::cout
  #define LOG_TRACE2(logger) if(true) {} else std::cout
  #define LOG_TRACE(logger) if(true) {} else std::cout

  #define DEFINE_LOG(logger, loggerName)   // no-op: logger var never referenced in release
  #define EXTERN_LOG(logger)               // no-op
  #define IS_LOG_ENABLED(a, b) false  

  class LogConfigurator  
  {
  public:
    //! Read logging configuration file (log.xml)
    static void ParseLogConfFile(const std::string& confFile);
  };
#else  
  // in debug we have the logging implemented 
  // enables "named" sinks
  BOOST_LOG_ATTRIBUTE_KEYWORD(bAttrClassName, "bAttrClassName", std::string)
  
  // try to avoid TRACE - the are always written to the log/console
  #define LOG_TRACE(logger) BOOST_LOG_SEV(logger, ::CoupledField::level_trace) << __LINE__ << "] "
  #define LOG_TRACE2(logger) BOOST_LOG_SEV(logger, ::CoupledField::level_trace2) << __LINE__ << "] "
  // use e.g. for basic information about major functions
  #define LOG_DBG(logger) BOOST_LOG_SEV(logger, ::CoupledField::level_dbg) << __LINE__ << "] "
  // detailed info to understand more how a function works 
  #define LOG_DBG3(logger) BOOST_LOG_SEV(logger, ::CoupledField::level_dbg3) << __LINE__ << "] "
  // deepest level for debugging the core of functions
  #define LOG_DBG2(logger) BOOST_LOG_SEV(logger, ::CoupledField::level_dbg2) << __LINE__ << "] "


  // this creates a new logger. E.g. after DEFINE_LOG(mylog,"mylog") use LOG_DBG(mylog) << "hans"; see log.xml
  #define DEFINE_LOG(logger, loggerName) auto logger = ::CoupledField::LogConfigurator::getLogger(loggerName);

  // this allows to use a logger which was defined already elsewhere
  #define EXTERN_LOG(logger) extern boost::log::sources::severity_logger<int> logger;

  #define IS_LOG_ENABLED(asdf, jkl) false /// TODO: check every occurrence

  //! This class manages the output for the different logging streams (only per script)
  class LogConfigurator  
  {
  public:

    //! Read logging configuration file (log.xml)
    static void ParseLogConfFile(const std::string& confFile);
    
    virtual ~LogConfigurator() {}

    static boost::log::sources::severity_logger<int> getLogger(const std::string& loggerName);
  };
#endif // end of the debug case
} // end of namespace CoupledField
#endif // end of CFS_LOG_HH
