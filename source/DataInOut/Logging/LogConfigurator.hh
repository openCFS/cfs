#ifndef CFS_LOG_HH
#define CFS_LOG_HH

// include files of BOOST's log library
#include <boost/log/common.hpp>
#include <boost/log/sinks.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/sources/severity_logger.hpp>
#include <boost/core/null_deleter.hpp>
#include <boost/limits.hpp>

#include <memory>

// enables "named" sinks
BOOST_LOG_ATTRIBUTE_KEYWORD(bAttrClassName, "bAttrClassName", std::string)

namespace CoupledField { namespace logging {
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
}
}

#ifndef NDEBUG
#define LOG_DBG3(logger) BOOST_LOG_SEV(logger, ::CoupledField::logging::level_dbg3) << __LINE__ << "] "
#define LOG_DBG2(logger) BOOST_LOG_SEV(logger, ::CoupledField::logging::level_dbg2) << __LINE__ << "] "
#define LOG_DBG(logger) BOOST_LOG_SEV(logger, ::CoupledField::logging::level_dbg) << __LINE__ << "] "
#define LOG_TRACE2(logger) BOOST_LOG_SEV(logger, ::CoupledField::logging::level_trace2) << __LINE__ << "] "
#define LOG_TRACE(logger) BOOST_LOG_SEV(logger, ::CoupledField::logging::level_trace) << __LINE__ << "] "
#else
#define LOG_DBG3(logger) if(true) {} else std::cout
#define LOG_DBG2(logger) if(true) {} else std::cout
#define LOG_DBG(logger) if(true) {} else std::cout
#define LOG_TRACE2(logger) if(true) {} else std::cout
#define LOG_TRACE(logger) if(true) {} else std::cout
#endif

#define DEFINE_LOG(logger, loggerName) auto logger = ::CoupledField::LogConfigurator::getLogger(loggerName);

#define EXTERN_LOG(logger) extern boost::log::sources::severity_logger<int> logger;

#define IS_LOG_ENABLED(asdf, jkl) false /// TODO: check every occourence

namespace CoupledField {

  //! This class manages the output for the different logging streams (only per script)
  class LogConfigurator  {
    
  public:

    //! Read logging configuration file
    static void ParseLogConfFile(const std::string& confFile);
    
    //! Destructor
    virtual ~LogConfigurator() {}
    
    //! Set logging level
    /*static void SetLogLevel( const std::string& logStream,
                      const std::string& level );*/

    static boost::log::sources::severity_logger<int> getLogger(const std::string& loggerName);

  private:

    //std::atomic<int> loggingLevel;

    //! Constructor - not needed, so private
    LogConfigurator(void);
  };

}

#endif
