#include "LogConfigurator.hh"
#include "General/Exception.hh"
#include "DataInOut/ProgramOptions.hh"

#include "DataInOut/ParamHandling/Xerces.hh"

namespace CoupledField {

  
  LogConfigurator::LogConfigurator(const std::string& confFile) :
    confFile_(confFile)
  {    
    // Set some default output modifiers for the log streams
    add_modifier("*", prepend_prefix);
    //add_modifier("*", prepend_time("$hh:$mm:$ss "), "", INT_MAX );
    add_modifier("*", append_enter);
    flush_log_cache();
  }


  void LogConfigurator::AddAppender(const std::string& logStream,
                                    const std::string& appender,
                                    const std::string& detail ) {
    if( appender == "cout" ) {
      add_appender( logStream, write_to_cout);
    } else if ( appender == "file" ) {
      add_appender( logStream, write_to_file(detail));
    } else {
      EXCEPTION("LogConfigurator knows only appender of type 'cout' "
                << "and 'file', but not '" << appender << "'");
    }
    flush_log_cache();
  }

  void LogConfigurator::SetLogLevel( const std::string& logStream,
                                     const std::string& level ) {
    if( level == "off" ) {
      manipulate_logs(logStream).disable();
    } else if( level == "dbg3" ) {
      manipulate_logs(logStream).enable(::boost::logging::level::dbg3);
    } else if( level == "dbg2" ) {
      manipulate_logs(logStream).enable(::boost::logging::level::dbg2);
    } else if( level == "dbg" ) {
      manipulate_logs(logStream).enable(::boost::logging::level::dbg);
    } else if( level == "trace2" ) {
      manipulate_logs(logStream).enable(::boost::logging::level::trace2);
    } else if( level == "trace" ) {
      manipulate_logs(logStream).enable(::boost::logging::level::trace);
    } else if( level == "all" ) {
      manipulate_logs(logStream).enable(::boost::logging::level::enable_all);
    } else {
      EXCEPTION("Log level '" << level << "' not defined!");
    }
    flush_log_cache();
  }

  void LogConfigurator::ParseLogConfFile() {
   
    // Query filename from option parser
    if( confFile_ == "") return;
    
    // If file is defined, create new Xerces instance and parse file
    Xerces xmlParser;
    xmlParser.SetFile(confFile_);
    
    // Obtain ParamNode
    PtrParamNode root = xmlParser.CreateParamNodeInstance();
    
    // Get list of classes
    
    if( root->GetName() != "logging" ) {
      EXCEPTION("Root-Node for logging definition must have name <logging>");
    }
    
    ParamNodeList classes  = root->GetChildren();
    
    // Loop over all classes
    for( UInt iClass = 0; iClass < classes.GetSize(); ++iClass ) {

      PtrParamNode classNode = classes[iClass];
      std::string className = classNode->Get("name")->As<std::string>();
      std::string level = classNode->Get("level")->As<std::string>();
      
      // Call for each class the SetLogLevel-method
      this->SetLogLevel(className, level );

      // Loop over all children (= appender definitions)
      ParamNodeList appenderNodes = classNode->GetChildren();

      // Call AddAppender for every childnode
      for( UInt iApp = 0; iApp < appenderNodes.GetSize(); ++iApp ) {
        PtrParamNode appNode = appenderNodes[iApp];
        
        // distinguish type of appender
        if( appNode->GetName() == "cout") {
          this->AddAppender(className, "cout", "");
        }
        if( appNode->GetName() == "file") {
          std::string fileName = appNode->Get("name")->As<std::string>();
          this->AddAppender(className, "file", fileName);
        }
      }
    }
  }



}
