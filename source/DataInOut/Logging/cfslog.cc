// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "cfslog.hh"



namespace CoupledField {

  
  LogConfigurator::LogConfigurator() {
    
    RegisterFunctions();

    // Set some default output modifiers for the log streams
    add_modifier("*", prepend_prefix);
    //add_modifier("*", prepend_time("$hh:$mm:$ss "), "", INT_MAX );
    add_modifier("*", append_enter);
    flush_log_cache();
  }


  void LogConfigurator::AddAppender() {
    SCRIPT_GET( std::string, log ); 
    SCRIPT_GET( std::string, appender ); 
    SCRIPT_GET( std::string, detail ); 
    
    if( appender == "cout" ) {
      add_appender( log, write_to_cout);
    } else if ( appender == "file" ) {
      add_appender( log, write_to_file(detail));
    } else {
      std::cerr << "Not working"  << std::endl;
    }
    flush_log_cache();
  }

  void LogConfigurator::AddModifier() {

  }

  void LogConfigurator::SetLogLevel() {
    SCRIPT_GET( std::string, log ); 
    SCRIPT_GET( std::string, level ); 

    if( level == "off" ) {
      manipulate_logs(log).disable();
    } else if( level == "dbg3" ) {
      manipulate_logs(log).enable(::boost::logging::level::dbg3);
    } else if( level == "dbg2" ) {
      manipulate_logs(log).enable(::boost::logging::level::dbg2);
    } else if( level == "dbg" ) {
      manipulate_logs(log).enable(::boost::logging::level::dbg);
    } else if( level == "trace2" ) {
      manipulate_logs(log).enable(::boost::logging::level::trace2);
    } else if( level == "trace" ) {
      manipulate_logs(log).enable(::boost::logging::level::trace);
    } else if( level == "all" ) {
      manipulate_logs(log).enable(::boost::logging::level::enable_all);
    }
    flush_log_cache();
  }

  void LogConfigurator::RegisterFunctions() {

    typedef FctPointer<LogConfigurator> FCPT;
    StdVector<ArgList> a;
    StdVector<FCPT*> pt;
    StdVector<std::string> name;

    // --- AddAppender ---
    a.Push_back();
    a.Last().RegisterParam("log", ArgList::STRING );
    a.Last().RegisterParam("appender", ArgList::STRING );
    a.Last().RegisterParam("detail", ArgList::STRING );
    pt.Push_back( new FCPT( this, &LogConfigurator::AddAppender) );
    name.Push_back( "addAppender" );

    // --- SetLogLevel ---
    a.Push_back();
    a.Last().RegisterParam("log", ArgList::STRING );
    a.Last().RegisterParam("level", ArgList::STRING );
    pt.Push_back( new FCPT( this, &LogConfigurator::SetLogLevel) );
    name.Push_back( "setLogLevel" );

    // Now register all functions with scripting 
    for (unsigned int i = 0; i < pt.GetSize(); i++ ) {
      Script_RegisterFct(name[i], pt[i], a[i] );
    }

  }



}
