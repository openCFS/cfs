#include  "Settings.hh"

#include  <cstring>

#ifdef WIN32
#define NO_HOME_PATH_EXPANSION
#endif

#ifndef NO_HOME_PATH_EXPANSION
#include  "pwd.h"
#endif

/**********************************************************
 *
 * Settings.cpp                     
 * Uwe Fabricius                      code@uwefabricius.de
 * Version 2.08.00                              02.06.2004
 *
 **********************************************************/

// granularity for the buffer, into that a line of a
// configuration file is read.
#ifndef  SETTING_DB_LINE_GRANULARITY
#define  SETTING_DB_LINE_GRANULARITY  1024
#endif // SETTING_DB_LINE_GRANULARITY

// definition of directive strings
#ifndef  SETTING_DB_DIRECTIVES
#define  SETTING_DB_DIRECTIVES
#define  SETTING_DB_DIRECTIVE_SECTION_IN   ">>section_enter<<"
#define  SETTING_DB_DIRECTIVE_SECTION_OUT  ">>section_leave<<"
#define  SETTING_DB_DIRECTIVE_END          ">>end<<"     
#endif // SETTING_DB_DIRECTIVES

/**********************************************************/

#ifndef  NEW
#define  NEW  new(nothrow)
#endif

#ifdef  SUPPRESS_COLORED_OUTPUT
#define  RED(s)     s
#define  BLACK_(s)  s
#define  CYAN(s)    s
#else
#define  RED(s)    "\033[31m" << s << "\033[0m"
#define  BLACK_(s) "\033[1;m" << s << "\033[0m"
#define  CYAN(s)   "\033[36m" << s << "\033[0m"
#endif // SUPPRESS_COLORED_OUTPUT

#ifndef  FUNCTION_ERROR
#define  FUNCTION_ERROR(fns)  RED(fns) << "\n >> "
#endif // FUNCTION_ERROR
#ifndef  CALLED_FROM
#define  CALLED_FROM(fns)     "    called from " << CYAN(fns) << "\n"
#endif // CALLED_FROM
#ifndef  INTERPRET_ERROR
#define  INTERPRET_ERROR(param, setting, type) \
         "    parameter " << BLACK_("\""<<param<<"\"") << " of setting " \
         << BLACK_("\""<<setting<<"\"") << "\n    cannot be interpreted "\
         "as " << type << " value -> ignoring this setting\n"
#endif // INTERPRET_ERROR

/**********************************************************/

namespace CoupledField {

/**********************************************************/

// definition of possible directives
const char* const SettingDataBase::
m_directives[SettingDataBase::N_DIRECTIVES] =
  { SETTING_DB_DIRECTIVE_SECTION_IN,
    SETTING_DB_DIRECTIVE_SECTION_OUT,
    SETTING_DB_DIRECTIVE_END };

/**********************************************************
 *
 * class Setting
 *
 **********************************************************/

  Setting::Setting( const SettingDef& settingdef,
                    bool        verbose     )
    : m_nParameters(0),
      m_minNumParams(settingdef.m_minNumParams),
      m_maxNumParams(settingdef.m_maxNumParams),
      m_Attributes(settingdef.m_Attributes),
      m_Type(settingdef.m_Type),
      m_IDString(NULL),
      m_verbose(verbose),
      m_bool(NULL),
      m_int(NULL),
      m_float(NULL),
      m_string(NULL),
      m_next(NULL)
  {
    // copy the ID-string
    if( !(m_IDString = NEW char [strlen(settingdef.m_IDString)+1]) ) {
      cerr << FUNCTION_ERROR( "CFSetting::CFSetting" )
           << " >> could not get, " << (strlen(settingdef.m_IDString)+1)
           << " bytes\n";
    } else {
      strcpy( m_IDString, settingdef.m_IDString );
    }
  }

  /**********************************************************/

  Setting::~Setting()
  {
    delete [] m_IDString;
    delete [] m_bool;
    delete [] m_int;
    if( m_string ) {
      for( int i = 0; i < m_nParameters; i++ ) {
        delete [] m_string[i];
      }
    }
    delete [] m_string;
    delete [] m_float;
    delete    m_next;
  }

  /**********************************************************
   * direct access to the parameters
   **********************************************************/

  bool Setting::getBool( const int i ) const {
    return checkIndexRange(i,"getBoolean",m_bool) ? m_bool[i] : false;
  }
  int Setting::getInt( const int i ) const {
    return checkIndexRange(i,"getInt",m_int) ? m_int[i] : -1;
  }
  char* Setting::getString( const int i ) const {
    return checkIndexRange(i,"getString",m_string) ? m_string[i] : NULL;
  }
  float Setting::getFloat( const int i ) const {
    return checkIndexRange(i,"getFloat",m_float) ? m_float[i] : -1;
  }

  /**********************************************************
   * Assignment operator
   * Assigns all member data of the Setting class, except m_next,
   * because this operator should only change the real data, not
   * their organisation.
   **********************************************************/

  Setting& Setting::operator = ( const Setting& setting )
  {
    // first clean this instance
    switch( m_Type ) {
    case BOOL :
      delete [] m_bool;
      m_bool = NULL;
      break;
    case INT :
      delete [] m_int;
      m_int = NULL;
      break;
    case FLOAT :
      delete [] m_float;
      m_float = NULL;
      break;
    case STRING :
      for( int i = 0; i < m_nParameters; i++ )  delete [] m_string[i];
      delete [] m_string;
      m_string = NULL;
      break;
    case FLAG :
    default :
      break;
    }
	
    // copy the member variables
    // DO NOT COPY : m_next, m_IDString
    m_nParameters  = setting.m_nParameters;
    m_minNumParams = setting.m_minNumParams;
    m_maxNumParams = setting.m_maxNumParams;
    m_Attributes   = setting.m_Attributes;
    m_Type         = setting.m_Type;
    m_verbose      = setting.m_verbose;
	
    // get m_IDString
    delete [] m_IDString;
    m_IDString = NEW char [strlen(setting.m_IDString)+1];
    strcpy( m_IDString, setting.m_IDString );
	
    // get data
    if( m_nParameters ) {
      switch( m_Type ) {
      case BOOL :
        if( NULL == (m_bool = NEW bool [m_nParameters]) ) {
          return *this; }
        memcpy( m_bool, setting.m_bool, m_nParameters*sizeof(bool) );
        break;
      case INT :
        if( NULL == (m_int = NEW int [m_nParameters]) ) {
          return *this; }
        memcpy( m_int, setting.m_int, m_nParameters*sizeof(int) );
        break;
      case FLOAT :
        if( NULL == (m_float = NEW float [m_nParameters]) ) {
          return *this; }
        memcpy( m_float, setting.m_float, m_nParameters*sizeof(float) );
        break;
      case STRING :
        if( NULL == (m_string = NEW char* [m_nParameters]) ) {
          return *this; }
        for( int i = 0; i < m_nParameters; i++ ) {
          m_string[i] = NEW char [strlen(setting.m_string[i])+1];
          strcpy( m_string[i], setting.m_string[i] );
        }
        break;
      case FLAG :
        break;
      }
    }
	
    return *this;
  }

  /**********************************************************
   *
   **********************************************************/

  bool Setting::expandPathString( const int i )
  {
    const char fn[] = "Setting::expandPathString";
	
    if( getType() != STRING ) {
      cerr << FUNCTION_ERROR(fn) << "setting " << BLACK_(getIDString())
           << " is not of type STRING -> cannot be expanded\n";
      return false;
    }
    if( getNumParameters() <= i ) {
      cerr << FUNCTION_ERROR(fn) << "parameter string ["<<i<<"] does "
           << "not exist in setting " << BLACK_(getIDString()) << "\n";
      return false;
    }
    if( false == expandHomePath(getStringPtr()+i) ) {
      cerr << CALLED_FROM(fn);
      return false;
    }

    return true;
  }


  bool Setting::expandAllPathStrings()
  {
    const char fn[] = "Setting::expandAllPathStrings";
	
    if( getType() != STRING ) {
      cerr << FUNCTION_ERROR(fn) << "setting " << BLACK_(getIDString())
           << " is not of type STRING -> cannot be expanded\n";
      return false;
    }

    bool result = true;
    for( int i = 0; i < getNumParameters(); i++ ) {
      result = (result && expandHomePath(getStringPtr()+i));
    }
    if( result == false )  cerr << CALLED_FROM(fn);

    return result;
  }

  /**********************************************************
   * Reads the parameters from the raw, unprepared configuration line.
   *
   * line         : Here the configuration line, terminated with a zero
   *                character is located.
   * FileName     : Contains the name of the file, the line is from.
   * LineNumber   : Line number of the line in file "FileName".
   * pickySuccess : if pickySuccess != NULL, *pickySuccess will be set to
   *                false, if there occured any, even the slightest, error.
   *              
   * "FileName" and "LineNumber" are only used for error messages. Their
   * default values, NULL and -1, cause all error messages to refere to
   * the command line instead of a file.
   *
   * Return value: true,  if this instance contains correct parameters
   *                      after the call
   *               false, otherwise
   **********************************************************/

  bool Setting::readParametersFromRawLine(       char *const line,
                                                 const char *const FileName,
                                                 const int         LineNumber,
                                                 bool *const pickySuccess )
  {
    int nStrings = 0;
	
    // separate the parameters with 0-characters
    nStrings = separateParametersInRawLine( line,       FileName,
                                            LineNumber, pickySuccess );
    // extract the parameters
    if( nStrings > 0 ) {
      if( false == readParametersFromPreparedLine(nStrings, line,
                                                  FileName, LineNumber,
                                                  pickySuccess) ) {
        return false;
      }
    }
	
    return true;
  }

  /**********************************************************
   * Append "pSetting" at this instance.
   *
   * WARNING: "pSetting" is appended exactly at THIS instance, not
   * at the last Setting of the list, this instance is part of.
   * To append a "pSetting" at the end of the list after this instance
   * call Setting::getLast()->append( pSetting ).
   *
   * Return value: Pointer to the previous next setting at m_next.
   *               This makes it easy to insert a setting S1 after a
   *               a setting S2, just calling
   *               S1->append( S2->append( S1 );
   **********************************************************/

  Setting* Setting::append( const Setting* const pSetting )
  {
    Setting *old = m_next;
    m_next = (Setting*)pSetting;
    return old;
  }

  /**********************************************************
   * Returns a pointer to the last element in the list starting
   * with this instance.
   **********************************************************/

  Setting* Setting::getLast() const
  {
    const Setting *last = this;
    while( last->m_next != NULL )  last = last->m_next;
    return (Setting*)last;
  }

  /**********************************************************
   * Separates the parameters in a raw setting line with 0 characters.
   * 
   * line         : Here the configuration line, terminated with a zero
   *                character is located.
   * FileName     : Contains the name of the file, the line is from.
   * LineNumber   : Line number of the line in file "FileName".
   * pickySuccess : if pickySuccess != NULL, *pickySuccess will be set to
   *                false, if there occured any, even the slightest, error.
   *              
   * "FileName" and "LineNumber" are only used for error messages. Their
   * default values, NULL and -1, cause all error messages to refere to
   * the command line instead of a file.
   *
   * Return value: The number of strings that have been separated in the
   *               line. -1, if an error occured.
   **********************************************************/

  int Setting::separateParametersInRawLine(       char *const line,
                                                  const char *const FileName,
                                                  const int         LineNumber,
                                                  bool *const pickySuccess )
  {
    const char FunctionName[] = "Setting::separateParametersInRawLine";
	
    int nStrings = 0,
      iread    = 0,
      iwrite   = 0;
	
    // check IDString and explicit assignment
    while( line[iread] == ' ' )  iread++;  // ignore leading spaces
    if(    line[iread] ==  0  )  return 0; // -> line contains only spaces
    // write IDString at the top of the line
    while( line[iread] != ' ' && line[iread] != '=' && line[iread] != 0 ) {
      line[iwrite++] = line[iread++];
    }
    // search next non-space character
    if( line[iread] != '=' )  while( line[iread] == ' ' )  iread++;
    // check explicit assignment
    if( line[iread] == '=' ) {
      line[iwrite++] = 0;
      iread++;
    } else {
      if( m_Attributes & EXPLICIT_ASSIGNMENT ) {
        SettingDataBase::printLineError( FunctionName, FileName, LineNumber );
        cerr << "    parameters for setting " << BLACK_('\"'<<m_IDString<<'\"')
             << " must be assigned explicitly with =\n";
        if( pickySuccess )  *pickySuccess = false;
        return 0;
      } else {
        line[iwrite++] = 0;
      }
    }
    nStrings++;
	
    // separate the parameters of this setting
    int insideQuoting   = 0,
      insideParameter = 0;
	
    while( line[iread] != 0 )
      {
        switch( line[iread] ) {
          // \ : just copy the next character
        case '\\' :
          line[iwrite++] = line[++iread];
          break;
          // " : enter or leave a string (or simply ignore)
        case '"' :
          if( m_Type == STRING ) {
            insideParameter = (insideQuoting ^= 1);
            if( insideParameter == 0 ) {
              nStrings++;
              line[iwrite++] = 0;
            }
          }
          // if type is not string, simply ignore the "
          break;
          // space : only copy, if it belongs to a string
        case ' ' :
          if( m_Type == STRING && insideQuoting ) {
            line[iwrite++] = line[iread];
          } else if( insideParameter ) {
            insideParameter = insideQuoting = 0;
            line[iwrite++] = 0;
            nStrings++;
          }
          break;
          // comma : eventually separating parameter
        case ',' :
          if( m_Attributes & COMMA_SEPARATED ) {
            if( insideQuoting ) {
              line[iwrite++] = line[iread];
            } else if( insideParameter ) {
              insideParameter = 0;
              line[iwrite++] = 0;
              nStrings++;
            }
          }
          break;
          // copy all other characters
        default :
          insideParameter = 1;
          line[iwrite++] = line[iread];
          break;
        }
        iread++;
      }
    if( insideParameter )  nStrings++;
    line[iwrite] = 0;

    return nStrings;
  }

  /**********************************************************
   * Extracts parameters from a prepared line in a passed buffer.
   *
   * line         : contains all parameter strings, separated by
   *                '\0'-characters, ALSO containing the IDString.
   * nParams      : number of strings in "line", also including the ID-string
   * FileName     : contains a file name, if a configuaration file is the
   *                source of the content of "line"
   * pickySuccess : if pickySuccess != NULL, *pickySuccess will be set to
   *                false, if there occured any, even the slightest, error.
   *
   * Return value: true,  if all parameters could be read properly
   *               false, otherwise. This could be caused by one of the
   *               following reasons.
   *   -> the setting has been read from the command line, but is only
   *      permitted in a configuration file ore vice versa
   *   -> the number of parameters does not reach the required minimum
   *   -> there is not enough memory to store the parameters (no error
   *      message is printed in this case)
   *   -> one of the parameters could not be read as one of the expected
   *      type
   **********************************************************/

  bool Setting::readParametersFromPreparedLine(       int         nStrings,
                                                      const char* const line,
                                                      const char* const FileName,
                                                      const int         LineNumber,
                                                      bool* const pickySuccess )
  {
    const char FunctionName[] = "Setting::readParametersFromPreparedLine";
    int nParameters = nStrings - 1;

    // check, if the source of this setting is valid
    // This might seem not to be the best place for this check, but the
    // programm behaves best, if it is checked here
    if( ((m_Attributes & Setting::COMMAND_LINE_ONLY) &&  FileName) ||
        ((m_Attributes & Setting::CONFIG_FILE_ONLY ) && !FileName)    )
      {
        SettingDataBase::printLineError( FunctionName, FileName, LineNumber );
        cerr << "    setting " << BLACK_('\"'<<m_IDString<<'\"') << " is "
          "only allowed in " << (FileName ? "the command line" :
                                 "a configuration file") << "\n -> setting ignored\n";
        if( pickySuccess )  *pickySuccess = false;
        return false;
      }
	
    // check number of parameters
    if( nParameters < m_minNumParams ) {
      SettingDataBase::printLineError( FunctionName, FileName, LineNumber );
      cerr << "    setting " << BLACK_('\"'<<line<<'\"') << " takes at "
        "least " << m_minNumParams << " parameter"
           << (m_minNumParams > 1 ? "s" : "") << " -> setting ignored\n";
      if( pickySuccess )  *pickySuccess = false;
      return false;
    }
    if( nParameters > m_maxNumParams && m_maxNumParams >= 0 ) {
      SettingDataBase::printLineError( FunctionName, FileName, LineNumber );
      cerr << "    setting " << BLACK_('\"'<<line<<'\"') << " takes at "
        "most " << m_maxNumParams << " parameter"
           << (m_minNumParams > 1 ? "s" : "")
           << "\n    -> overlapping parameters ignored\n";
      if( pickySuccess )  *pickySuccess = false;
      // truncate to maximal number of parameters
      nParameters = m_maxNumParams;
    }

    m_nParameters = 0;

    // create memory for parameters
    switch( m_Type ) {
    case BOOL :
      delete [] m_bool;
      if( !(m_bool = NEW bool [nParameters]) ) {
        if( pickySuccess )  *pickySuccess = false;
        return false;
      }
      break;
    case INT :
      delete [] m_int;
      if( !(m_int = NEW int [nParameters]) ) {
        if( pickySuccess )  *pickySuccess = false;
        return false;
      }
      break;
    case FLOAT :
      delete [] m_float;
      if( !(m_float = NEW float [nParameters]) ) {
        if( pickySuccess )  *pickySuccess = false;
        return false;
      }
      break;
    case STRING :
      delete [] m_string;
      if( !(m_string = NEW char* [nParameters]) ) {
        if( pickySuccess )  *pickySuccess = false;
        return false;
      }
      break;
    case FLAG :
      break;
    }
	
    // get parameters
    const char *pparam = line;
    for( int i = 0; i < nParameters; i++ ) {
      // seek next parameter
      while( *pparam++ != 0 );
      // read parameter
      switch( m_Type )
        {
        case BOOL :
          if( strcmp(pparam, "true") == 0 ||
              strcmp(pparam, "yes" ) == 0 ||
              strcmp(pparam, "1"   ) == 0    )
            {
              m_bool[i] = true;
            } else if( strcmp(pparam, "false") == 0 ||
                       strcmp(pparam, "no"   ) == 0 ||
                       strcmp(pparam, "0"    ) == 0    )
              {
                m_bool[i] = false;
              } else {
                SettingDataBase::printLineError( FunctionName, FileName, LineNumber );
                cerr << INTERPRET_ERROR(pparam, line, "boolean" );
                if( pickySuccess )  *pickySuccess = false;
                return false;
              }
          break;
				
        case INT :
          if( 1 != sscanf(pparam, "%d", &m_int[i]) ) {
            SettingDataBase::printLineError( FunctionName, FileName, LineNumber );
            cerr << INTERPRET_ERROR(pparam, line, "integer" );
            if( pickySuccess )  *pickySuccess = false;
            return false;
          }
          break;
				
        case FLOAT :
          if( 1 != sscanf(pparam, "%f", &m_float[i]) ) {
            SettingDataBase::printLineError( FunctionName, FileName, LineNumber );
            cerr << INTERPRET_ERROR(pparam, line, "float" );
            if( pickySuccess )  *pickySuccess = false;
            return false;
          }
          break;
			
        case STRING :
          if( NULL == (m_string[i] = NEW char [strlen(pparam)+1]) ) {
            if( pickySuccess )  *pickySuccess = false;
            return false;
          }
          strcpy( m_string[i], pparam );
          break;
				
        case FLAG :
          break;
        }
		
      m_nParameters++;
    }
	
    return true;
  }

  /**********************************************************
   * 
   **********************************************************/

  bool Setting::expandHomePath( char** path )
  {
    const char fn[] = "Setting::expandHomePath";
	
    if( *path == NULL || (*path)[0] == 0 ) {
      cerr << FUNCTION_ERROR(fn) << "invalid parameter string\n";
      return false;
    }
    // nothing to expand, nothing to do
    if( (*path)[0] != '~' )  return true;
    // only "~" is not a valid path
    if( (*path)[1] == 0 ) {
      cerr << FUNCTION_ERROR(fn) << "cannot expand string \"~\"\n";
      return false;
    }

#ifdef NO_HOME_PATH_EXPANSION

    cerr << FUNCTION_ERROR(fn) << "the home path expansion has been "
      "deactivated by compiling with NO_HOME_PATH_EXPANSION\n"
      "    cannot expand path " << *path << "\n";
    return false;
	
#else

    const char *expansion       = NULL;
    char       *pslash          = (*path + 1);
    int         expansionLength = 0;
    // only expand HOME path of the current user
    if( *pslash == '/' ) {
      // get home path from environment variable
      expansion = getenv( "HOME" );
      if( NULL ==  expansion ||
          0    >= (expansionLength = strlen(expansion)) ) {
        cerr << FUNCTION_ERROR(fn) << "could not expand HOME "
          "path of current user\n";
        return false;
      }
      // '~' is not followed by a '/', so we need the full path
      // of a different user's home directory
    } else {
      // temporarily crop to the user's name
      while( *pslash != '/' && *pslash )  pslash++;
      const bool cropped = (*pslash == '/');
      *pslash = 0;
      // get the user's full home path from /etc/pwd
      struct passwd *pwdentry = getpwnam( *path+1 );
      if( pwdentry == NULL ) {
        cerr << FUNCTION_ERROR(fn) << "could not read /etc/passwd entry"
          "for user " << BLACK_('\"'<<(*path+1)<<'\"') << "\n";
        // restore path
        if( cropped ) *pslash = '/';
        return false;
      }
      if( NULL ==  pwdentry->pw_dir ||
          0    >= (expansionLength = strlen(pwdentry->pw_dir)) ) {
        cerr << FUNCTION_ERROR(fn) << "could not expand HOME path "
          "of user " << BLACK_('\"'<<(*path+1)<<'\"') << "\n";
        // restore path
        if( cropped ) *pslash = '/';
        return false;
      }
      expansion = pwdentry->pw_dir;
      // restore path
      if( cropped ) *pslash = '/';
    }

    // create new string including expansion
    char *fullPath = NEW char[expansionLength + strlen(pslash) + 1];
    if( fullPath == NULL ) {
      cerr << FUNCTION_ERROR(fn) << "could not get memory for "
        "expanded path\n";
      return false;
    }
    // compose full path
    strcpy( fullPath, expansion );
    strcpy( fullPath + expansionLength, pslash );

    delete [] *path;
    *path = fullPath;
	
    return true;

#endif // NO_HOME_PATH_EXPANSION
  }

  /**********************************************************
   * Simply checks the validity of the index used to access a
   * parameter in the array containing the read parameters.
   **********************************************************/

  bool Setting::checkIndexRange( const int         i,
                                 const char* const FunctionName,
                                 const void* const pointer ) const
  {
    if( i < 0 || i >= m_nParameters ) {
#ifdef _WARNINGS_
      cerr << ERROR_FUNCTION( "Setting::" << FunctionName )
           << "i out of range [0," << m_nParameters-1 << "]\n";
#endif
      return false;
    } else if( pointer == NULL ) {
#ifdef _WARNINGS_
      cerr << ERROR_FUNCTION( "Setting::" << FunctionName )
           << "the setting is not of the guessed type\n";
#endif
      return false;
    } else {
      return true;
    }
  }

  /**********************************************************
   *
   * class SettingDef
   *
   **********************************************************/

  SettingDef::SettingDef( const char *const    IDString,
                          const char *const    altIDString,
                          const Setting::type  Type,
                          const int            minNumParams,
                          const int            maxNumParams,
                          const int            Attributes,
                          const char *const    HelpString    )
    : m_IDString(IDString),
      m_altIDString(altIDString),
      m_Type(Type),
      m_Attributes(Attributes),
      m_minNumParams(minNumParams),
      m_maxNumParams(maxNumParams),
      m_HelpString(HelpString)
  {
    // prohibit less than zero parameter
    if( m_minNumParams < 0 )  m_minNumParams = 0;
    // avoid sensless values (m_maxNumParams < 0 means, that there
    // doesn't exist an upper limit for the number of parameters)
    if( m_minNumParams > m_maxNumParams && m_maxNumParams > 0 )
      m_maxNumParams = m_minNumParams;
    // assure that a flag has zero parameters
    if( m_Type == Setting::FLAG )  m_minNumParams = m_maxNumParams = 0;
  }

  /**********************************************************/

  SettingDef::~SettingDef()
  {}

  /**********************************************************
   * Print operator for the setting definition.
   **********************************************************/

  ostream& operator << ( ostream& os, const SettingDef& sd )
  {
    if( sd.m_IDString == NULL )  return os;
	
    const char* TypeStrings[Setting::STRING+1] =
      { "flag", "boolean", "integer", "float", "string" };
	
#ifndef  SUPPRESS_COLORED_OUTPUT
    os << "\033[1m" << sd.m_IDString;
#else
    os << sd.m_IDString;
#endif

    if( sd.m_altIDString )  os << ", " << sd.m_altIDString;

#ifndef  SUPPRESS_COLORED_OUTPUT
    os << "\033[0m : [";
#else
    os << " : [";
#endif

    if( sd.m_Type != Setting::FLAG ) {
      os << sd.m_minNumParams;
      if( sd.m_minNumParams != sd.m_maxNumParams ) {
        os << "-";
        if( sd.m_maxNumParams > 0 )  os << sd.m_maxNumParams;
        else                         os << "?";
      }
    }
    os << " " << TypeStrings[sd.m_Type]
       << (sd.m_maxNumParams > 1 || sd.m_maxNumParams < 0 ? "s" : "") << "]\n";
	
    if( sd.m_HelpString )  cout << "    " << sd.m_HelpString << endl;
	
    return os;
  }

  /**********************************************************
   *
   * class SettingDataBase
   *
   **********************************************************/
		
  SettingDataBase::SettingDataBase( const bool verbose )
    : m_verbose(verbose),
      m_Settings(NULL),
      m_activeSection(NULL),
      m_insideSection(false),
      m_currentLineLength(0)
  {
  }

  /**********************************************************/

  SettingDataBase::~SettingDataBase()
  {
    delete    m_Settings;
    delete [] m_activeSection;
  }

  /**********************************************************/

  Setting* SettingDataBase::getSetting( const char* const IDString ) const
  {
    Setting* setting = m_Settings;
    while( setting != NULL ) {
      if( 0 == strcmp(IDString, setting->getIDString()) )  break;
      setting = setting->getNext();
    }
    return setting;
  }

  /**********************************************************/

  Setting* SettingDataBase::getSetting( int index ) const
  {
    Setting *setting = m_Settings;
    while( index-- > 0 && setting )  setting = setting->getNext();
    return setting;
  }

  /**********************************************************/

  int SettingDataBase::getNumParameters( const char* const IDString ) const
  {
    Setting* setting = getSetting( IDString );
    return setting != NULL ? setting->getNumParameters() : -1;
  }

  /**********************************************************/

  bool SettingDataBase::expandPaths()
  {
    const char fn[] = "SettingDataBase::expandPaths";

    Setting* setting = m_Settings;
    bool     result  = true;
    while( setting ) {
      if( setting->getType() == Setting::STRING ) {
        result = (result && setting->expandAllPathStrings());
      }
      setting = setting->getNext();
    }
    if( result == false )  cerr << CALLED_FROM(fn);

    return result;
  }

  /**********************************************************
   * Searches the string getSetting(ID)->getString(0) in the passed
   * array SettingStrings. If it can be found, the index in
   * SettingStrings is returned.
   *
   * ID             : ID of the setting, the string is defined in,
   *                  that should be searched for.
   * SettingStrings : Array of strings, the string getSetting(ID)->
   *                  getString(0) should be searched in.
   * maxID          : Either SettingStrings is terminated with a NULL
   *                  string or the length of SettingStrings is passed
   *                  as maxID.
   **********************************************************/
 
  int SettingDataBase::mapStringToIndex( const char*  const ID,
                                         const char** const SettingStrings,
                                         const int          stringIndex,
                                         const bool         verbose ) const
  {
    const char fn[] = "SettingDataBase::mapStringToIndex";
    Setting *setting = getSetting( ID );
	
    // get pointer to matching setting
    if( setting == NULL ) {
      if( verbose ) {
        cerr << FUNCTION_ERROR( fn ) << "setting "
             << BLACK_(ID) << " has not yet been set\n";
      }
      return -1;
    }
    // check type of the setting
    if( setting->getType() != Setting::STRING ) {
      if( verbose ) {
        cerr << FUNCTION_ERROR( fn ) << "setting "
             << BLACK_(ID) << " is not of type STRING\n";
      }
      return -1;
    }
    //
    if( setting->getNumParameters() <= stringIndex ) {
      if( verbose ) {
        cerr << FUNCTION_ERROR( fn ) << "for setting " << BLACK_(ID)
             << " only " << setting->getNumParameters() << " parameter"
             << (setting->getNumParameters() > 1 ? "s are" : " is")
             << " known\n     -> [" << stringIndex << "] is not a valid"
          " parameter index\n";
      }
      return -1;

    }

    int index = 0;
    const char* const pParamString = setting->getString( stringIndex );

    // seek matching string in array
    while( SettingStrings[index] != NULL &&
           strcmp(SettingStrings[index], pParamString) )  index++;
    // if it could not be found -> index < 0
    if( SettingStrings[index] == NULL ) {
      if( verbose ) {
        cerr << FUNCTION_ERROR( fn ) << "none of the passed "
          "strings matched in setting " << BLACK_(ID) << "\n";
      }
      index = -1;
    }

    return index;	
  }

  /**********************************************************
   * Reads Settings from a file.
   *
   * settingdef : contains the definition of possible settings.
   * FileName   : name of the setting file.
   * overwrite  : if true, all already present settings are overwritten,
   *              by new ones with the same ID-String.
   * anyErrors  : an additional monitor for the success of the function.
   *              The return value of readSettings is only false, if there
   *              occured a really critical error, that ommited a sensible
   *              proceding of the configuration file. But even if the
   *              return value is true, there might have occured some smaller
   *              errors, which is indicated by setting pickySuccess to false,
   *              if the pointer pickySuccess is passed as != NULL.
   *
   * Return value: true,  if the settings could be read
   *               false, otherwise
   **********************************************************/

  bool SettingDataBase::readSettings( const SettingDef* const settingdef,
                                      const char*       const FileName,
                                      const bool              overwrite,
                                      bool*       const pickySuccess )
  {
    const char FunctionName[] = "SettingDataBase::readSettings";
	
    // preparations
    // initialize pickySuccess
    if( pickySuccess )  *pickySuccess = true;
    // open configuration file
    FILE *file = fopen( FileName, "rt" );
    if( file == NULL ) {
      cerr << FUNCTION_ERROR(FunctionName) << "could not open "
           << BLACK_('\"'<<FileName<<'\"') << "\n";
      if( pickySuccess )  *pickySuccess = false;
      return false;
    }

    char *line   = NULL;  // pointer to the line buffer
    int   length = 0,     // length of the read line
      number = 0;     // line number (for error messages)
	
    // read line after line until file end
    while( 0 <= (length = readLine(file, &line, FileName, ++number)) ) {
      // find matching definition in "settingdef"
      int idef = findMatchingSetting( line,       // pointer to the line
                                      length,     // line lenght
                                      settingdef, // setting definitions
                                      FileName,   // name of the config file
                                      number,     // line number in the config file
                                      false,      // NOT quiet, print errors
                                      pickySuccess );
      if( idef >= 0 ) {
        // create new setting
        Setting *newSetting = NEW Setting( settingdef[idef] );
        // try to read the parameters and to integrate the new setting
        if( false == (newSetting->readParametersFromRawLine(line,   FileName,
                                                            number, pickySuccess)) ||
            false == integrateSetting(newSetting, overwrite, pickySuccess) )
          {
            delete newSetting;
          }
      }
    }	


    // clean up
    fclose( file );
    delete [] line;
    m_currentLineLength = 0;
	
    return true;
  }

  /**********************************************************
   * Read settings from the command line.
   *
   * settingdef  : array with definition of possible settings.
   * CommandLine : command line as passed to main()
   * nargs       : number of arguments at "CommandLine"
   * expectDashesInCommandLine : if true, all settings defined in settingdef
   *               must be affected by a "-"
   * owerwrite   : if true, all already present settings are overwritten,
   *               by new ones with the same ID-String.
   *
   * Return value: returns true, errors are printed out
   **********************************************************/

  bool SettingDataBase::readSettings( const SettingDef* const settingdef,
                                      const char**      const CommandLine,
                                      const int               nargs,
                                      const bool              overwrite,
                                      const bool              expectDashesInCommandLine,
                                      bool*       const pickySuccess )
  {
    // assume, that CommandLine[0] contains the program name
    int   istart =    1,
      iend   =    1,
      idef_1 =   -1,
      idef_2 =   -1;
    char *line   = NULL;

    // initialize pickySuccess
    if( pickySuccess )  *pickySuccess = true;
    // check, if the first program parameter calls help
    if( nargs > 1 && !strcmp(CommandLine[1], "-help") ) {
      cout << "\nsettings for " << BLACK_(CommandLine[0]) << "\n";
      printHelp( settingdef );
      cout << endl;
      return false;
    }
	
    // scan command line
    while( istart < nargs && iend < nargs ) {
      // end of the previous setting is now start of the next one
      istart = iend;
      // check, if the argument can be identified as a setting (with output)
      idef_1 = identifySettingInCommandLine( settingdef, CommandLine,
                                             istart, expectDashesInCommandLine,
                                             false, pickySuccess );
      if( 0 <= idef_1 ) {
        // walk with "iend" to the next setting
        while( ++iend < nargs ) {
          idef_2 = identifySettingInCommandLine( settingdef, CommandLine,
                                                 iend, expectDashesInCommandLine,
                                                 false, pickySuccess );
          // leave loop, if next setting could be found OR
          // setting should not be used in command line
          if( idef_2 >= 0 || idef_2 == -2 )  break;
        }
        // build a configuration file line from the command line arguments
        buildConfigLineFromCommandLine( CommandLine, istart, iend,
                                        &line, expectDashesInCommandLine );
        // now process the line like it was a line from a config file
        Setting *newSetting = NEW Setting( settingdef[idef_1] );
        if( !newSetting->readParametersFromRawLine(line, NULL, -1,
                                                   pickySuccess  ) ||
            !integrateSetting(newSetting, overwrite, pickySuccess)    )
          {
            delete newSetting;
          }
      } else {
        iend = ++istart;
      }
    }	
	
    delete [] line;
    m_currentLineLength = 0;
	
    return true;
  }

  /**********************************************************
   * Prints help, including string defined for this purpose in settingdef.
   **********************************************************/

  void SettingDataBase::printHelp( const SettingDef* const settingdef ) const
  {
    for( int i = 0; settingdef[i].m_IDString; i++ ) {
      cout << settingdef[i];
    }
  }

  /**********************************************************
   * Restricts the reading of the file to a section
   **********************************************************/

  bool SettingDataBase::restrictToSection( const char* const section  )
  {
    if( section == NULL ) {
      delete [] m_activeSection;
      m_insideSection = false;
      m_activeSection = NULL;
    } else {
      // calculate old and new lengths
      int oldLength = m_activeSection ? strlen(m_activeSection) : 0,
        newLength = strlen(section),
        // allocate buffer in chunks of size 20
        oldBufferLength = 20 * ((int)((oldLength+1) / 20) +
                                ((oldLength+1) % 20 ? 1 : 0)),
        newBufferLength = 20 * ((int)((newLength+1) / 20) +
                                ((newLength+1) % 20 ? 1 : 0));
      // eventually get new buffer for the restriction string
      if( oldBufferLength != newBufferLength || NULL == m_activeSection ) {
        delete [] m_activeSection;
        if( NULL == (m_activeSection = NEW char [newBufferLength]) ) {
          cerr << FUNCTION_ERROR("SettingDataBase::restrictToSection")
               << "could not get enough memory for storing "
            "restriction name\n";
          return false;
        }
      }
      // copy new restriction name
      memcpy( m_activeSection, section, sizeof(char)*(newLength+1) );
    }
	
    return true;
  }

  /**********************************************************
   * Finally check the completeness of the settings, using the definitions
   * in "settingdef".
   **********************************************************/

  bool SettingDataBase::finalCheck( const SettingDef* const settingdef ) const
  {
    bool     result  = true;
    const char FunctionName[] = "SettingDataBase::finalCheck";
	
    for( int i = 0; settingdef[i].m_IDString; i++ ) {
      // check, if all mandatory settings are present
      if( (settingdef[i].m_Attributes & Setting::MANDATORY) &&
          NULL == getSetting(settingdef[i].m_IDString) ) {
        cerr << FUNCTION_ERROR( FunctionName );
        cerr << BLACK_('\"'<<settingdef[i].m_IDString<<'\"')
             << " must be specified\n";
        result = false;
      }
    }
	
    if( false == result ) {
      cerr << "\n try command line option " << BLACK_("\"-help\"\n");
    }
	
    return result;
  }

  /**********************************************************
   * Small function to print out an error with file and line information.
   **********************************************************/

  void SettingDataBase::printLineError( const char* const FunctionName,
                                        const char* const FileName,
                                        const int         LineNumber )
  {
    cerr << FUNCTION_ERROR( FunctionName );
    if( FileName ) {
      cerr << "error in file " << BLACK_('\"'<<FileName<<'\"')
           << ", line " << BLACK_(LineNumber) << "\n";
    } else {
      cerr << "error in command line option\n";
    }
  }

  /**********************************************************
   * Reads on line from the file "file".
   *
   * Reads characters from "file" into the memory block starting at line.
   * Reading stops, if the line ends with a '\n', or a comment starts
   * with '#'. Here it is assumed that the buffer at "*line" has got the
   * length m_currentLineLength. If the lenght of the line read from the
   * file exceedes this lenght, the buffer is enlarged successively by
   * SETTING_DB_LINE_GRANULARITY.
   *
   * Return value: Size of the read line
   **********************************************************/

  int SettingDataBase::readLine(       FILE*       file,
                                       char**      line,
                                       const char* const FileName,
                                       const int         LineNumber )
  {
    int   LineLength     = -1,
      nchars         =  0;
    char *character      = *line - 1,
      FunctionName[] = "SettingDataBase::readLine";
    bool  insidequoting  = false;
	
    // read line until comment or line end ...
    do {// if the line buffer is full, try to enlarge it
      if( (int)((LineLength + 1)*sizeof(char)) >= m_currentLineLength ) {
        char *oldline = *line;
        // == "if( could not enlarge line buffer )"
        if( oldline == (*line = enlargeLineBuffer(*line)) ) {
          cerr << CALLED_FROM( FunctionName );
          return -1;
        } else {
          // position character in the new buffer
          character = *line + m_currentLineLength - 1
            - SETTING_DB_LINE_GRANULARITY;
        }
      }
      LineLength += (nchars = fread( ++character, sizeof(char), 1, file ));
      // check quoting
      if( *character == '"'  )  insidequoting = !insidequoting;
      // check line break by a backslash
      if( *character == '\\' && !insidequoting ) {
        // remove backslash from line buffer
        *character = '\n';
        int lineChar = 0;
        // search end of line
        while( '\n' != (lineChar = fgetc(file)) && EOF != lineChar );
        if( lineChar == EOF )  break;
        // continue reading the logical line
        else { character--;  LineLength--; }
      }
    }
    while(    (*character != '#'  || insidequoting) // comment starts
              &&  *character != '\n'                   // line ends
              &&   nchars );                           // file ends
	
    // ... and search start of the next line
    if( *character == '#' ) {
      do  fread( character, 1, 1, file );
      while( *character != '\n' ); }
    // terminate line
    *character = 0;	

    // prescan line
    return prescanLine( *line, LineLength, FileName, LineNumber );
  }

  /**********************************************************
   * Prescans the raw line in particular for direcives.
   *
   * line        : buffer, containing the raw line
   * FileName    : name of the file, if the line is from a file
   * LineNumber  : line number, if the line is from a file
   *
   * Return value: line length to be returned by
   *               SettingDataBase::readLine
   **********************************************************/
 
  int SettingDataBase::prescanLine( const char* const line,
                                    const int         length,
                                    const char* const FileName,
                                    const int         LineNumber )
  {
    if( length <= 4 )  return length;

    const char FunctionName[] = "SettingDataBase::prescanLine";
	
    int dir   = 0,
      iline = 0;
    // test key ">>" for directives
    if( line[0] == '>' && line[1] == '>' ) {
      // identify directive
      while( dir < N_DIRECTIVES ) {
        int idir  = 0;
        iline = 0;
        while( iline < length && m_directives[dir][idir] ) {
          if( line[iline] != m_directives[dir][idir] )  break;
          iline++;
          idir++;
        }
        if( !m_directives[dir][idir] )  break;
        dir++;
      }
      // unknown directive
      if( dir >= N_DIRECTIVES ) {
        printLineError( FunctionName, FileName, LineNumber );
        cerr << "    unknown directive -> ignoring this line\n";
        return 0;
      }
      // no directive in this line
    } else {
      // there is an active section, but this line is not inside
      if( isRestricted() && !m_insideSection )  return 0;
      return length;                   
    }

    // select the action for the found directive
    switch( dir ) {
      //////////////////////
      // section directive
    case DIRECTIVE_SECTION_ENTER:
    case DIRECTIVE_SECTION_LEAVE: {
      // ignore sections, if no restriction is active
      if( !isRestricted() )  return 0;
      // seek start of section name
      while( iline < length &&
             (line[iline] == ' ' || line[iline] == '\t') )  iline++;
      // section directives need a section name
      if( iline == length ) {
        printLineError( FunctionName, FileName, LineNumber );
        cerr << "    directive " << BLACK_('\"'<<m_directives[dir]<<'\"')
             << " needs a section name\n    -> ignored directive\n";
        return 0;
      }
      // compare section name
      int isec = 0;
      while( m_activeSection[isec] &&
             m_activeSection[isec] == line[iline] ) {
        isec++;  iline++;
      }
      // section name matches
      if( !m_activeSection[isec] &&
          (!line[iline] || line[iline] == ' ' || line[iline] == '\t') ) {
        // adapt position flag
        m_insideSection = (dir == DIRECTIVE_SECTION_ENTER);
      }
      break;
    }
      /////////////////////////////
                                    // stop reading immidiately
                                    case DIRECTIVE_END: {
                                      if( isRestricted() && m_insideSection )  return -1;
                                      break;
                                    }
    }

    return 0;
  }

  /**********************************************************
   * Gets a raw line read from the file by readLine and tries to find a
   * matching setting in settingdefs, i.e. find the first string in the
   * line at "line" and try to identify a setting at "settingdefs".
   * NOTE: "line" is changed, but restored.
   *
   * line        : buffer, containing the raw line
   * linelength  : length of the line
   * settingdefs : setting definitions
   * FileName    : name of the file, if the line is from a file
   * LineNumber  : line number, if the line is from a file
   * quiet       : if true, no error message is printed out (we need that
   *               for the testing of the parameters in the command line)
   *
   * Return value: The index of the matching setting at settingdefs.
   *               -1, if no matching setting could be found.
   **********************************************************/

  int SettingDataBase::findMatchingSetting(       char*             line,
                                                  const int               linelength,
                                                  const SettingDef *const settingdefs,
                                                  const char       *const FileName,
                                                  const int               LineNumber,
                                                  const bool              quiet,
                                                  bool       *const pickySuccess )
  {
    const char FunctionName[] = "SettingDataBase::findMatchingSetting";
	
    int  istart = 0,
      iend   = 0,
      idef   = 0;
    char ctemp  = 0;
	
    // seek start of the IDString
    while( line[istart] == ' ' )  istart++;  // skip leading spaces
    if( istart >= linelength-1 )  return -1; // -> line contains only spaces
    // seek end of the IDString
    iend = istart;
    while( line[iend] != ' ' && line[iend] != '=' && line[iend] != 0 )  iend++;
    // terminate the first string
    ctemp      = line[iend];
    line[iend] = 0;
	
    // search for a matching setting definition
    for( idef = 0; settingdefs[idef].m_IDString; idef++ ) {
      if(  0 == strcmp(line+istart, settingdefs[idef].m_IDString) ||
           ( settingdefs[idef].m_altIDString &&
             !strcmp(line+istart, settingdefs[idef].m_altIDString)) ) {
        break;
      }
    }
    // check, if a matching setting definition exists
    if( settingdefs[idef].m_IDString == NULL ) {
      if( !quiet ) {
        printLineError( FunctionName, FileName, LineNumber );
        cerr << "    " << BLACK_('\"'<<(line+istart)<<'\"')
             << " is not a valid setting -> ignored\n";
      }
      if( pickySuccess )  *pickySuccess = false;
      return -1;
    }
	
    // restore the line
    line[iend] = ctemp;
	
    return idef;
  }

  /**********************************************************
   * Identify a setting in the command line.
   * A little augmented version of "findMatchingSetting"
   **********************************************************/

  int SettingDataBase::
  identifySettingInCommandLine( const SettingDef* const settingdef,
                                const char**      const CommandLine,
                                const int               i,
                                const bool              expectDashesInCommandLine,
                                const bool              quiet,
                                bool*       const pickySuccess )
  {
    int jstart =  0;
    const char FunctionName[] = "SettingDataBase::identifySettingInCommandLine";
	
    if( expectDashesInCommandLine ) {
      if( CommandLine[i][0] == '-' ) {
        jstart = 1;
      } else {
        if( quiet == false ) {
          cerr << FUNCTION_ERROR( FunctionName )
               << "expected " << BLACK_('\"'<<CommandLine[i]<<'\"')
               << " being a setting but a \"-\" must precede\n"
            "    -> skipping to next command line argument\n";
        }
        if( pickySuccess )  *pickySuccess = false;
        return -1;
      }
    }
	
    return findMatchingSetting( (char*)(CommandLine[i]+jstart),
                                strlen( CommandLine[i]+jstart ),
                                settingdef, NULL, -1, quiet,
                                pickySuccess );
  }

  /**********************************************************
   * Builds a configuration file line from one command line setting.
   *
   * CommandLine : command line as passed to main()
   * istart      : analyse the command line starting at CommandLine[istart] ...
   * iend        : ... ending at CommandLine[iend]
   * line        : buffer, where the line should be stored to
   * expectDashesInCommandLine : if true, all settings defined in settingdef
   *               must be affected by a "-"
   *
   * Return value: length of the built configuration line
   *               -1, if an error occured
   **********************************************************/

  int SettingDataBase::
  buildConfigLineFromCommandLine( const char** const CommandLine,
                                  const int          istart,
                                  const int          iend,
                                  char**       line,
                                  const bool         expectDashesInCommandLine )
  {
    int  LineLength =  0,
      iread      = expectDashesInCommandLine ? 1 : 0;
    bool quoting    = false;
    const char FunctionName[] = "SettingDataBase::buildConfigLineFromCommandLine";
	
    for( int iargument = istart; iargument < iend; iargument++ ) {
      // if the argument string contains a space, it must be quoted in the line
      if( NULL != strchr(CommandLine[iargument]+iread, ' ') )  quoting = true;
      // test, if the buffer can take this argument, eventually including quotes
      // (2 characters) and one space and a terminating zero
      if( (int)((  LineLength
                   + strlen(CommandLine[iargument]+iread)
                   + (quoting ? 4 : 2) ) * sizeof(char)) > m_currentLineLength )
        {   // try to enlarge buffer
          char *oldline = *line;
          if( oldline == (*line = enlargeLineBuffer(*line)) ) {
            cerr << CALLED_FROM(FunctionName);
            return -1;
          }
        }
		
      // opening quote, if marked so
      if( quoting )  (*line)[LineLength++] = '"';
      // copy argument
      strcpy( (*line)+LineLength, CommandLine[iargument]+iread );
      // Although strlen returns the string length NOT including the terminating
      // zero, but strcpy DOES copy it, we only enlarge the line by this return
      // value, because we do not want this terminating zero in the line.
      LineLength += strlen( CommandLine[iargument]+iread );
      // closing quote, if marked so
      if( quoting )  (*line)[LineLength++] = '"';
      (*line)[LineLength++] = ' ';
		
      iread = 0;
    }
    // terminate line (in particular if rests of an old,
    // longer line is still present in the buffer)
    (*line)[LineLength++] = 0;
	
    return LineLength;
  }			

  /**********************************************************
   * Integrates a new setting into the list at m_Settings.
   * If a setting with the same IDString already exists, the flag
   * "overwrite" is respected.
   *
   * Return value : true,  if "setting" has been used. That means, that
   *                       it must not be simply deleted in future
   *                false, if "setting" has not been used. It could be
   *                       deleted without harming the list at m_Settings
   **********************************************************/

  bool SettingDataBase::integrateSetting(       Setting *const setting,
                                                const bool           overwrite,
                                                bool    *const pickySuccess )
  {
    // search for an already present setting with the same IDString
    Setting* oldSetting = getSetting( setting->getIDString() );

    // if the setting is already present
    if( oldSetting ) {
      // if the present setting should be overwritten, i.e. replaced
      if( overwrite ) {
        // remove and delete present setting
        Setting *pSetting = m_Settings;
        while( pSetting->getNext() != oldSetting ) {
          pSetting = pSetting->getNext();
        }
        pSetting->append( oldSetting->append(NULL) );
        delete oldSetting;
        // insert the new setting at the head of the list. Note
        // that this code is also correct, if there is no setting
        // present before (m_Settings == NULL).
        setting->append( m_Settings );
        m_Settings = setting;
        // the setting is present and should not be overwritten
      } else {
        return false;
      }
      // setting not yet present
    } else {
      // insert the new setting at the head of the list
      setting->append( m_Settings );
      m_Settings = setting;
    }

    return true;
  }

  /**********************************************************
   * Helping function to enlarge the line buffer at "line".
   * A new buffer with a length augmented by SETTING_DB_LINE_GRANULARITY
   * is allocated, and the content of "line" is copied into it. The new
   * part of the buffer is initialised with 0. If "line" == NULL , this
   * function simply works like a class intern realloc, enlarging just
   * by SETTING_DB_LINE_GRANULARITY chars.
   * return value: Pointer to the new line buffer. If the buffer could
   *               not be enlarged, the old pointer is returned.
   **********************************************************/

  char* SettingDataBase::enlargeLineBuffer( char* line )
  {
    // allocate memory with enlarged length
    char* newLine = NEW char [m_currentLineLength + SETTING_DB_LINE_GRANULARITY];
	
    if( newLine == NULL ) {
      cerr << FUNCTION_ERROR("SettingDataBase::enlargeLineBuffer")
           << "could to enlarge line buffer up to "
           << m_currentLineLength + SETTING_DB_LINE_GRANULARITY
           << " bytes\n";
      newLine = line;
    } else {
      // copy old line into new one
      if( line != NULL ) {
        memcpy( newLine, line, sizeof(char)*m_currentLineLength );
      }
      // initialize new part of the line with 0
      memset( newLine + m_currentLineLength, 0, SETTING_DB_LINE_GRANULARITY );
      m_currentLineLength += SETTING_DB_LINE_GRANULARITY;
      delete [] line;
    }
	
    return newLine;
  }

  /**********************************************************/

} // end of namespace
