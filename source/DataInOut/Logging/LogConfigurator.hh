#ifndef CFS_LOG_HH
#define CFS_LOG_HH


// include files of BOOST's (unofficial) log library
#include "DataInOut/Logging/log.hpp"
#include "DataInOut/Logging/functions.hpp"
#include <boost/limits.hpp>

// this is a dirty hack to prevent issues with parallel loggin.
// We only log for thread 0. It would probably better to move to the current
// boost logging were a multithreaded version is provided: _mt
#include "def_use_openmp.hh"
#ifdef USE_OPENMP
  #include <omp.h>
#endif

using namespace boost::logging;

// Create customized levels
BOOST_LOG_DEFINE_LEVEL(dbg2, 1300)
BOOST_LOG_DEFINE_LEVEL(dbg3, 1200)
BOOST_LOG_DEFINE_LEVEL(trace, 1500)
BOOST_LOG_DEFINE_LEVEL(trace2, 1450)

// Create own logging statements, in order to ease the
// use of logging
#ifndef NDEBUG
#define DECLARE_LOG(log_name) BOOST_DECLARE_LOG(log_name)
#else
#define DECLARE_LOG(log_name) BOOST_DECLARE_LOG_COMPILE_TIME(log_name,0)
#endif

#define LOG_DEFINE_LEVEL(lvl, value) BOOST_LOG_DEFINE_LEVEL(lvl,value)
#define DEFINE_LOG(log_name,log_str) BOOST_DEFINE_LOG(log_name,log_str) 
#define IS_LOG_ENABLED(log_name,lvl) BOOST_IS_LOG_ENABLED(log_name,lvl)


#define LOGL(log_name,lvl) BOOST_LOGL(log_name,lvl)
#define LOG(log_name) BOOST_LOG(log_name)


// see the comments at #include "def_use_openmp.hh"
// we have to set -Wno-dangling-else. Remove it in compiler.cmake if this construct is removed
#define THREAD_0_LOGL(log_name, lvl) if(LogConfigurator::GetThreadNum() == 0) BOOST_LOGL(log_name,lvl)


#ifndef NDEBUG
#define LOG_DBG3(log_name) THREAD_0_LOGL(log_name,dbg3)
#define LOG_DBG2(log_name) THREAD_0_LOGL(log_name,dbg2)
#define LOG_DBG(log_name) THREAD_0_LOGL(log_name,dbg)
#define LOG_TRACE2(log_name) THREAD_0_LOGL(log_name,trace2)
#define LOG_TRACE(log_name) THREAD_0_LOGL(log_name,trace)
#else
#define LOG_DBG3(log_name) if(false) std::cout
#define LOG_DBG2(log_name) if(false) std::cout
#define LOG_DBG(log_name) if(false) std::cout
#define LOG_TRACE2(log_name) if(false) std::cout
#define LOG_TRACE(log_name) if(false) std::cout
#endif

namespace CoupledField {
  
  //! This class manages the output for the different logging streams (only per script)
  class LogConfigurator  {
    
  public:
    
    //! Constructor 
    LogConfigurator(const std::string& confFile);
    
    //! Destructor
    virtual ~LogConfigurator() {}
    
    //! Read logging configuration file
    void ParseLogConfFile();
    
    //! Add output destination for given stream
    void AddAppender( const std::string& logStream,
                      const std::string& destination,
                      const std::string& destDetail );
    
    //! Set logging level
    void SetLogLevel( const std::string& logStream,
                      const std::string& level );

    static int GetThreadNum()
    {
    #ifdef USE_OPENMP
      return 2; //omp_get_thread_num();
    #else
      return 0;
    #endif
    }
    
  protected:
    
    //! Name of XML file containing the settings for logging.
    std::string confFile_;
  };

}

#endif
