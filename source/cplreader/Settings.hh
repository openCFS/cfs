// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_SETTINGS_2006
#define FILE_SETTINGS_2006

#include <map>
#include <string>
#include <sstream>
#include <memory>

#include <boost/cstdint.hpp>

#include <def_cplreader.hh>

#include "DataInOut/ParamHandling/Xerces.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"

namespace CoupledField
{

  //! \class Settings settings.hh "settings.hh"
  //! \brief A class to store and retrieve settings.
  //!
  //! The Settings class can store strings, ints and doubles
  //! associated with a name. The settings can be accessed from
  //! everywhere by taking advantage of the singleton concept.

  class Settings
  {
  public:
    // ========================================================================
    // TYPEDEFS
    // ========================================================================

    //! Map type for storing settings
    typedef std::map<std::string, std::string> Map;

  public:

    //! Get the instance
    //!
    //! Due to the fact that this class uses private
    //! constructors this method is the only means of
    //! creating/getting the only instance of the
    //! Settings class.
    static Settings& Instance()
    {
      if(!settingsInstance_.get())
      {
        settingsInstance_.reset( new Settings() );
      }

      return *settingsInstance_;
    }

    //! Get string
    //!
    //! This method gets the string with the name
    //! key from the Settings instance.
    //! \return string value belonging to key
    std::string GetString(const std::string& key)
    {
      if(settingsMap_.find(key) == settingsMap_.end())
      {
        std::stringstream sstr;

        EXCEPTION("Key '" << key << "' not found in map.");
      }

      return settingsMap_[key];
    }

    //! Get int
    //!
    //! This method gets the integer with the name
    //! key from the Settings instance.
    //! \return int value belonging to key
    boost::int32_t GetInt(const std::string& key)
    {
      if(settingsMap_.find(key) == settingsMap_.end())
      {
        EXCEPTION("Key '" << key << "' not found in map.");
      }

      std::stringstream sstr;
      int32_t value;
      bool failed;

      sstr << settingsMap_[key];
      failed = !(( sstr >> value ) && ( sstr.eof() ));

      if(failed)
      {
        EXCEPTION("Conversion to int32_t failed.");
      }

      return value;
    }

    //! Get double
    //!
    //! This method gets the double with the name
    //! key from the Settings instance.
    //! \return double value belonging to key
    double GetDouble(const std::string& key)
    {
      if(settingsMap_.find(key) == settingsMap_.end())
      {
        EXCEPTION("Key '" << key << "' not found in map.");
      }

      std::stringstream sstr;
      double value;
      bool failed;

      sstr << settingsMap_[key];
      failed = !(( sstr >> value ) && ( sstr.eof() ));

      if(failed)
      {
        EXCEPTION("Conversion to double failed.");
      }

      return value;
    }

    //! Set string
    //!
    //! This method associates the string val with
    //! the name key.
    //! \param key (input) key
    //! \param value (input) string value
    void SetString(const std::string& key, const std::string& value)
    {
      settingsMap_[key] = value;
    }

    //! Set int
    //!
    //! This method associates the integer val with
    //! the name key.
    //! \param key (input) key
    //! \param value (input) int value
    void SetInt(const std::string& key, boost::int32_t value)
    {
      std::stringstream sstr;

      sstr << value;

      settingsMap_[key] = sstr.str();
    }

    //! Set double
    //!
    //! This method associates the double val with
    //! the name key.
    //! \param key (input) key
    //! \param value (input) double value
    void SetDouble(const std::string& key, double value)
    {
      std::stringstream sstr;

      sstr << value;

      settingsMap_[key] = sstr.str();
    }

    void DumpXML(std::string& dump)
    {
      std::stringstream sstr;
      Map::iterator it, end;
       
      sstr << "<?xml version=\"1.0\"?>" << std::endl;

      sstr << "<cplreaderSettings>" << std::endl;

      it = settingsMap_.begin();
      end = settingsMap_.end();
      
      for( ; it != end; it++ )
      {
        sstr << "  <" << it->first << "> "
             << it->second
             << " </" << it->first << ">"
             << std::endl;
      }
      
      sstr << "</cplreaderSettings>" << std::endl;

      dump = sstr.str();
    }

    void InitParamNode(const std::string& xmlFile)
    {
      // Initialize our xerces dom parser to handle the cfs xml file
      Xerces* xerces = new Xerces(xmlFile, "");
      
      // set the global ParamNode tree pointer
      param_ = xerces->CreateParamNodeInstance();
        // save us in the info stuff, with defaults but no comments
      // todo: info->Get(ParamNode::HEADER)->Get("cfsSimulation")->SetValue(param);
      // release the xerces ressources, param is not affected
      delete xerces;
    }

    PtrParamNode GetParamNode()
    {
      return param_;
    }

  private:

    //! Private constructor.
    //!
    //! This one can never be called from outside.
    //! The only way to instantiate a Settings object is
    //! by using Instance().
    Settings() {param_ = NULL; };

    //! Private copy constructor.
    //!
    //! This one can never be called from outside.
    //! The only way to instantiate a Settings object is
    //! by using Instance().
    Settings(const Settings& s) {};

    //! Here the actual values will be stored
    Map settingsMap_;

    PtrParamNode param_;

    //! This one holds the only Settings instance
    static std::auto_ptr<Settings> settingsInstance_;
  };
}


#endif // FILE_SETTINGS_2006

/// Local Variables:
/// mode: C++
/// c-basic-offset: 2
/// End:
