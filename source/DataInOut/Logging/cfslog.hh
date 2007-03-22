// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef CFS_LOG_HH
#define CFS_LOG_HH

#include "DataInOut/Scripting/scriptable.hh"

// include files of BOOST's (unofficial) log library
#include "DataInOut/Logging/log.hpp"
#include "DataInOut/Logging/functions.hpp"
#include <boost/limits.hpp>

using namespace boost::logging;

// Create customized levels
BOOST_LOG_DEFINE_LEVEL(dbg2, 1300)
BOOST_LOG_DEFINE_LEVEL(dbg3, 1200)
BOOST_LOG_DEFINE_LEVEL(trace, 1500)
BOOST_LOG_DEFINE_LEVEL(trace2, 1450)

// Create own logging statements, in order to ease the
// use of logging
#ifdef DEBUG
#define DECLARE_LOG(log_name) BOOST_DECLARE_LOG(log_name)
#else
#define DECLARE_LOG(log_name) BOOST_DECLARE_LOG_COMPILE_TIME(log_name,0)
#endif

#define LOG_DEFINE_LEVEL(lvl, value) BOOST_LOG_DEFINE_LEVEL(lvl,value)
#define DEFINE_LOG(log_name,log_str) BOOST_DEFINE_LOG(log_name,log_str) 
#define IS_LOG_ENABLED(log_name,lvl) BOOST_IS_LOG_ENABLED(log_name,lvl)


#define LOGL(log_name,lvl) BOOST_LOGL(log_name,lvl)
#define LOG(log_name) BOOST_LOG(log_name)


#define LOG_DBG3(log_name) BOOST_LOGL(log_name,dbg3)
#define LOG_DBG2(log_name) BOOST_LOGL(log_name,dbg2)
#define LOG_DBG(log_name) BOOST_LOGL(log_name,dbg)
#define LOG_TRACE2(log_name) BOOST_LOGL(log_name,trace2)
#define LOG_TRACE(log_name) BOOST_LOGL(log_name,trace)

namespace CoupledField {
  
  //! This class manages the output for the different logging streams (only per script)
  class LogConfigurator : public Scriptable {
    
  public:
    
    //! Constructor 
    LogConfigurator();
    
    //! Add output destination for given stream
    void AddAppender();
    
    //! Add modifier
    void AddModifier();

    //! Set logging level
    void SetLogLevel();
    
  protected:
    
    //! Register functions with scripting object
    void RegisterFunctions();
  };

}

#endif
