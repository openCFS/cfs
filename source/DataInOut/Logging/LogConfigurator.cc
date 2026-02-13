#include "LogConfigurator.hh"
#include "General/Exception.hh"
#include "DataInOut/ProgramOptions.hh"

#include "DataInOut/ParamHandling/XmlReader.hh"

#include <fstream>
#include <boost/smart_ptr/shared_ptr.hpp>
#include <boost/smart_ptr/make_shared_object.hpp>
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/sinks/sync_frontend.hpp>
#include <boost/log/sinks/text_ostream_backend.hpp>
#include <boost/log/sources/severity_logger.hpp>
#include <boost/log/sources/record_ostream.hpp>
#include <boost/log/utility/setup/common_attributes.hpp>

using namespace boost::log;
namespace expr = boost::log::expressions;

namespace CoupledField {

// an alias to simplify the template stuff
typedef sinks::synchronous_sink<sinks::text_ostream_backend> text_sink;

// needed to extract the correct level from xml file
const std::map<std::string, int> logLevelMap = {
    {"disable_all", logging::level_disable_all},
    {"default_"   , logging::level_default_},
    {"enable_all" , logging::level_enable_all},
    {"all"        , logging::level_enable_all},
    {"fatal"      , logging::level_fatal},
    {"err"        , logging::level_err},
    {"trace"      , logging::level_trace},
    {"trace2"     , logging::level_trace2},
    {"dbg"        , logging::level_dbg},
    {"dbg2"       , logging::level_dbg2},
    {"dbg3"       , logging::level_dbg3},
    {"warn"       , logging::level_warn},
    {"info"       , logging::level_info}
};

void LogConfigurator::ParseLogConfFile(const std::string& confFile) {

  boost::log::add_common_attributes();

  // Query filename from option parser
  if( confFile == "") {
    // Add default logger to prevent default logging from kicking in

    boost::shared_ptr<text_sink> sink = boost::make_shared<text_sink>();
    sink->set_filter(expressions::attr<int>("Severity") > 10000); // set to always fail

    boost::shared_ptr<std::ostream> stream { &std::cout, boost::null_deleter { } };
    sink->locked_backend()->add_stream(stream);

    core::get()->add_sink(sink);

    return;
  }

  // If file is defined, create new xml instance and parse file
  PtrParamNode root = XmlReader::ParseFile(confFile);

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

    // Remember log level
    int logLevel = logLevelMap.at(level);

    // Loop over all children (= appender definitions)
    ParamNodeList appenderNodes = classNode->GetChildren();

    boost::shared_ptr<text_sink> sink = boost::make_shared<text_sink>();
    
    //Check for wildcard 
    if(className == "*") {
      //Logging all registert loggers
      sink->set_filter(expr::attr<int>("Severity") >= logLevel);
    } else {
      sink->set_filter(expr::attr<int>("Severity") >= logLevel && expressions::attr<std::string>("bAttrClassName") == className);
    }
    // line number and "]" added in header file
    sink->set_formatter(expr::stream << " [" << bAttrClassName << ":" << expr::smessage);

    // add stream for every childnode
    for (UInt iApp = 0; iApp < appenderNodes.GetSize(); ++iApp) {
      PtrParamNode appNode = appenderNodes[iApp];

      // distinguish type of appender
      if (appNode->GetName() == "cout") {

        boost::shared_ptr<std::ostream> stream { &std::cout, boost::null_deleter { } };
        sink->locked_backend()->add_stream(stream);

      } else if (appNode->GetName() == "file") {
        std::string fileName = appNode->Get("name")->As<std::string>();
        sink->locked_backend()->add_stream(boost::make_shared< std::ofstream >(fileName));
      }
    }

    core::get()->add_sink(sink);
  }
}

boost::log::sources::severity_logger<int> LogConfigurator::getLogger(const std::string& loggerName) {
    boost::log::sources::severity_logger<int> logger = sources::severity_logger<int>();
    logger.add_attribute("bAttrClassName", boost::log::attributes::constant<std::string>(loggerName));
    return logger;
}

}
