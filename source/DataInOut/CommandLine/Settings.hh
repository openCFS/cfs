#ifndef FILE_SETTINGS_HH
#define FILE_SETTINGS_HH

#include  <cstdio>
#include  <iostream>

namespace CoupledField {

  /**********************************************************
   *
   * Settings.hpp                     
   * Uwe Fabricius                      code@uwefabricius.de
   * Version 2.08.00                              02.06.2004
   *
   **********************************************************
   *
   * Application example:
   *
   * --8<---<snip>----8<---<snip>----8<---<snip>---
   * SettingDef Defs[] =
   * {
   *     SettingDef("family",      NULL, Setting::STRING, 1,  4, Setting::COMMA_SEPARATED, "help for \"family\"" ),
   *     SettingDef("age",         NULL, Setting::INT,    0, -1, Setting::COMMA_SEPARATED, NULL ),
   *     SettingDef("intelligent",  "i", Setting::BOOL,   1,  4, 0, NULL ),
   *     SettingDef("athome",      NULL, Setting::STRING, 0, -1, 0, NULL ),
   *     SettingDef()
   * };
   *
   * int main( int nargs, const char **pargs )
   * {
   *     SettingDataBase settings;
   * 
   *     settings.readSettings( Defs, "config.test", true );
   *     settings.readSettings( Defs, pargs, nargs, true, true );
   *     settings.finalCheck( Defs );
   *
   *     for( int i = 0; Defs[i].m_IDString; i++ ) {
   *         Setting* s = settings.getSetting( Defs[i].m_IDString );
   *         if( s ) {
   *             cout << "\033[1m" << s->getIDString() << "\033[0m\n";
   *             for( int j = 0; j < s->getNumParameters(); j++ ) {
   *                 switch( s->getType() ) {
   *                     case Setting::BOOL :
   *                         cout << "   [" << j << "]: " << s->getBool(j) << endl;
   *                         break;
   *                     case Setting::INT :
   *                         cout << "   [" << j << "]: " << s->getInt(j) << endl;
   *                         break;
   *                     case Setting::FLOAT :
   *                         cout << "   [" << j << "]: " << s->getFloat(j) << endl;
   *                         break;
   *                     case Setting::STRING :
   *                         cout << "   [" << j << "]: " << s->getString(j) << endl;
   *                         break;
   *                     case Setting::FLAG :
   *                     default :
   *                         break;
   *                 }
   *             }
   *         } else {
   *             cout << "\033[1m" << Defs[i].m_IDString << "\033[0m : not specified\n";
   *         }
   *     }
   *     return true;
   * }
   * --8<---<snip>----8<---<snip>----8<---<snip>---
   *
   * sample configuration file:
   *
   * --8<---<snip>----8<---<snip>----8<---<snip>---
   * ## part of the family
   * family = "homer, with #" lisa,bart # marge, maggy
   *  age    45  40 12
   * intelligent= false yes 0
   * athome LISA MARGE
   * --8<---<snip>----8<---<snip>----8<---<snip>---
   *
   **********************************************************/

  using namespace std;

  /**********************************************************/

  class SettingDef;
  std::ostream& operator << ( std::ostream&, const SettingDef& );

  /**********************************************************/

  //! class to store <b>one</b> setting including its parameters

  /*! This class represents one single setting beloning to class
   *  SettingDataBase, consisting of an ID-string and the setting's
   *  parameters. The parameters may be of type int, bool, float,
   *  string, or simply used as a flag, but only one type per setting.
   *  A setting is created from an object of type class SettingDef.
   *
   *  <br><br>
   *  <b>PATH EXPANSION</b> <br>
   *  A setting of type STRING might represent a certain path in the file
   *  system. Therefore it may contain start with the abreviation '~' for
   *  the home path of the current user or for example "~homer" the home
   *  path of user "homer". Using such a path inside your program to open
   *  a file will not work, since the '~' or "~user" is expanded by the
   *  shell, if you use it in the terminal, but there is no shell inside
   *  your programm. You have to expand it manually in your program. class
   *  Setting provides this expansion with the methods expandPathString
   *  and expandAllPathStrings. Setting::expandPathString expands a single
   *  specified string in the setting, Setting::expandAllPathStrings expands
   *  all strings in the setting.
   */
  class Setting
  {
  public:

    //! constants for the types of the setting's parameters
    enum type { FLAG = 0, //!< a simple flag
                BOOL,     //!< boolean value, possible parameter strings are:
                //!  "true", "false", "yes", "no", "1", "0"
                INT,      //!< integer value
                FLOAT,    //!< float value
                STRING }; //!< a string
                
    //! some flags for the treatment of the setting

    /*! some flags for the treatment of the setting
     *   - COMMA_SEPARATED: The parameters <b>can</b> be separated
     *     with commas, i.e. "homer,lisa,bart" will be recognised
     *     as the three parameters "homer","lisa","bart". The
     *     specification "homer lisa bart" is understood nevertheless.
     *   - EXPLICIT_ASSIGNMENT:  The parameters must be assigned
     *     explicitly with a "=", where by spaces may separate the
     *     "=" from the setting identifier and the parameters.
     *   - MANDATORY: This setting is marked as mandatory. If all
     *     settings have been loaded, the function
     *     SettingDataBase::finalCheck() will print an error message
     *     for each mandatory setting that has not been set, and
     *     return false.
     *   - COMMAND_LINE_ONLY, CONFIG_FILE_ONLY: Allows restrictions
     *     concerning the source of a setting. (For example it does
     *     not make sense to allow the setting of a configuration file
     *     to specify the name or path of exactly that config file).
     */
    enum { COMMA_SEPARATED     = 1<<0,
           EXPLICIT_ASSIGNMENT = 1<<1,
           MANDATORY           = 1<<2,
           COMMAND_LINE_ONLY   = 1<<3,
           CONFIG_FILE_ONLY    = 1<<4 };

    //! creates a Setting object from a object of type SettingDef
    Setting( const SettingDef& settingdef,
             bool        verbose = true );
    ~Setting();

    //! \name data access
    //@{
    //! returns the number of present parameters
    int  getNumParameters() const { return m_nParameters; }
    //! returns the type of the setting (see the enumeration types)
    type getType() const { return m_Type; }
    //! returns the boolean nr. i in the setting
    bool  getBool( const int i = 0 ) const;
    //! returns the integer nr. i in the setting
    int   getInt( const int i = 0 ) const;
    //! returns the string nr. i in the setting
    char* getString( const int i = 0 ) const;
    //! returns the float nr. i in the setting
    float getFloat( const int i = 0 ) const;
    //! direct access to the array
    char*  getIDString()  const { return m_IDString; }
    //! direct access to the array
    bool*  getBoolPtr()   const { return m_bool; }
    //! direct access to the array
    int*   getIntPtr()    const { return m_int; }
    //! direct access to the array
    char** getStringPtr() const { return (char **const)m_string; }
    //! direct access to the array
    float* getFloatPtr()  const { return m_float; }
    //@}

    //! \name data manipulation
    //@{
    //! expands home path from '~'

    /*! String settings often represent a file or directory path.
     *  This path may be passed relatively to the home directory of
     *  a user, abreviated by '~', more precise by "~/" for the
     *  current user's home directory or "~homer/" for the home
     *  directory of the user "homer". Setting::expandPathString
     *  checks, if a string parameters starts with a '~' and tries
     *  to expand the string using the full path of the according
     *  home directory (so the string is changed irreversibly).
     *  Therefore the method uses the content of the environment
     *  variable "HOME" or the corresponding entry in "/etc/passwd".
     *  The necessary functions are declared in the standard library
     *  header "pwd.h". You can switch of this feature by
     *  <i>#defining</i> NO_HOME_PATH_EXPANSION. Then the method can
     *  still be called, but with no effect, and will return false
     *  whenever the string starts with a '~'.
     *  \param i index of the string in this setting, that should be
     *         expanded
     *  \return true, if the expansion was successfull. This also
     *          includes the case that the string did not start with
     *          a '~'. <br> false, if the string started with a '~'
     *          but could not be expanded, including due to a defined
     *          NO_HOME_PATH_EXPANSION.
     */
    bool expandPathString( const int i = 0 );

    //! expands home path form '~' in all strings

    /*! This function expands all present strings in this setting.
     *  So it is simply a loop of expandPathString(const int) over
     *  all strings.
     *  \return true, if all calls of expandPathString(const int)
     *          returned true, otherwise false
     */
    bool expandAllPathStrings();
    //@}
                
    // operators
    Setting& operator = ( const Setting& setting );
                

    // reading parameters
    // read parameters from one configuration line
    bool readParametersFromRawLine(       char* const line,
                                          const char* const FileName     = NULL,
                                          const int         LineNumber   =   -1,
                                          bool* const pickySuccess = NULL );

                
    //! \name administration routines for the linked list of Settings
    //@{
    Setting* append( const Setting* const pSetting );
    Setting* getLast() const;
    Setting* getNext() const { return m_next; }
    //@}

  protected:
                        
    // separate the parameters in a line by 0 characters
    int separateParametersInRawLine(       char *const line,
                                           const char *const FileName,
                                           const int         LineNumber,
                                           bool *const pickySuccess = NULL );
    // read parameters of a setting in a prepared, i.e zero-separated line
    bool readParametersFromPreparedLine(       int         nStrings,
                                               const char *const line,
                                               const char *const FileName,
                                               const int         LineNumber,
                                               bool *const pickySuccess = NULL );
    // expand a home directory path from '~'
    bool expandHomePath( char** path );
                                        
    // controlling data
    int         m_nParameters,  // # parameters
      m_minNumParams, // min # parameters
      m_maxNumParams, // max # parameters
      m_Attributes;
    type        m_Type;
    char       *m_IDString;
    bool        m_verbose;
    // parameters of the setting
    bool       *m_bool;
    int        *m_int;
    float      *m_float;
    char      **m_string;
    // pointer to the next setting
    Setting    *m_next;
                
  private:
                        
    bool checkIndexRange( const int i,
                          const char* const FunctionName,
                          const void* const pointer ) const;
  };

  /**********************************************************
   * Class only needed as parameter for Setting::Setting(..)
   **********************************************************/

  //! class representing the definition of a setting

  /*! This class is used in combination with class SettingDataBase.
   *  In principal an object of this class only exists to get passed to
   *  to the constructor of class Setting. The user of the settings
   *  framework only will use this class to create an array of SettingDef
   *  objects, that will be passed to the two versions of function
   *  SettingDataBase::readSettings.
   */
  class SettingDef
  {
    friend std::ostream& operator << ( std::ostream&, const SettingDef& );

  public:

    //! constructor

    /*! Since all parameters have a default value, it is possible
     *  to create a SettingDef object with the standard constructor.
     *  This should only be used to create an empty object to
     *  terminate the definition array, that is passed to the
     *  functions SettingDataBase::readSettings.
     *  \param IDString first identification string of a setting
     *  \param alternative IDString first identification string of
     *         a setting (optional)
     *  \param Type the type of the setting (see documentation of
     *         enumeration type Setting::type).
     *  \param minNumParams minimal number of parameters, that should
     *         be passed with this setting
     *  \param minNumParams maximal number of parameters, that should
     *         be passed with this setting. Default value -1 means
     *         <i>unlimited</i>
     *  \param attributes for the setting (combination of anonymous
     *         enumeration constants in class Setting)
     *  \param HelpString a string containing the text, that should
     *         be printed out as documentation of the accordingly
     *         generated Setting (also see SettingDataBase::printHelp)
     */
    SettingDef( const char *const    IDString     = NULL,
                const char *const    altIDString  = NULL,
                const Setting::type  Type         = Setting::FLAG,
                const int            minNumParams = -1,
                const int            maxNumParams = -1,
                const int            Attributes   =  0,
                const char *const    HelpString   = NULL );
    ~SettingDef();
                
                
    const char *const    m_IDString,      //! ID-string
      *const    m_altIDString;   //! alternative ID-String
    const Setting::type  m_Type;          //! storage type
    const int            m_Attributes;    //! attributes
    int                  m_minNumParams,  //! minimal number of parameters
      m_maxNumParams;  //! maximal number of parameters
    const char *const    m_HelpString;    //! string for help
  };

  /**********************************************************
   * Small class for the mapping of string parameters to other values.
   **********************************************************/

  //! class for the mapping of a string to another type
  template <class T>
  class SettingStringMap
  {
  public :

    SettingStringMap( const char* const string = 0,
                      T                 value  = 0  )
      : m_String(string),
        m_Value(value) {}
        
    const char* const  m_String;
    T                  m_Value;
  };

  /**********************************************************
   * ProgramSettings
   **********************************************************/

  //! class for reading settings from command line and configuration files

  /*! This is the class, you work with, if you want read the settings from
   *  the command line or a configuration file. <br>
   *  Create an instance of SettingDataBase with the constructor (optionally
   *  specifying the verbosity flag) and call readSettings, either the
   *  version for data files or the one for the command line options. You
   *  can use the same array of SettingDef objects for both versions of
   *  readSettings. This allows to read settings either from the command
   *  line or from a configuration file or from both. If you read the same
   *  settings from the command line <b>and</b> you can specify by the
   *  boolean parameter <code>overwrite</code>, if already set settings
   *  should be overwritten. <br>
   *  Both versions of readSettings are quite fault-tolerant. Even if there
   *  are some not readable settings in the source, the functions print
   *  out an error message, but continue reading setting as far as possible.
   *  Only in case of really heavy errors the return value will reflect
   *  this in a false. So by only checking the return value the calling
   *  routine cannot decide if really <b>everything</b> was fine. But
   *  readSettings take an optional parameter bool* pickySuccess. This
   *  pointer is passed all through the whole process of reading and
   *  interpreting settings and the variable refered by this pointer is
   *  initialized with true and set to false, whenever a slight error
   *  occures, even if this error does not harm the processing of the
   *  functions and does not cause a false as return value. So by passing
   *  this optional parameter one can still react on small faults.
   *  <br><br>
   *  <b>SYNTAX</b> of the configuration file:
   *   - one setting per line
   *   - everything after a '#' in a line is comment (including the '#')
   *   - ONLY and ALL directive lines start with ">>"
   *   - spaces are ignored, exept, if the separate parameters
   *   - doublequotes to define a string containing spaces, commas or
   *     hashes are recognised (see example)
   *   - overlapping parameters for a setting are ignored (remove comment
   *     in "family" line of the sample file
   *   - a backslash, that is not in included in quoting, causes a line
   *     break. There is no character inserted to separate the items
   *     befor and after this backslash. So you can continue a string
   *     over more that one line. For example (| symbolizes the left
   *     boundary of the text file):
   *     <br><br>
   * \code
   * |age = 10 20 \ this text is not recognized
   * |30
   * |height = 120 190\
   * |230
   * \endcode
   *     <br>
   *     This configuration file will cause: <br>
   *       age[0] = 10, age[1] = 20, age[2] = 30
   *       heigth[0] = 120, height[1] = 190230 
   *
   *  <br><br>
   *  <b>DIRECTIVES</b> <br>
   *  The syntax of the configuration files provides directives, that
   *  means key words, that influence the reading of the file dynamically.
   *  All directives lines must start with the first two characters ">>".
   *  Possible directives
   *   - >>section_enter<< sectionname        <br>
   *     >>section_leave<< sectionname        <br>
   *       Enter and leave a section the reading is restricted to. For
   *       setting and resetting use functions "restrictToSection" and
   *       and "removeSection". If the SettingDataBase class is restricted
   *       to a section "sectionname", only settings braced by these two
   *       directives are read.
   *   - >>end<<    <br>
   *       Reading stops immidiately. This directive might be usefull in
   *       developing states to hide developing rubbish without the need
   *       of repeated (un)commenting. Section restriction is stronger
   *       than ">>end<<", i.e. if the class is restricted to a section
   *       ">>end<<" is only respected, if placed inside this section.
   *
   *  <br><br>
   *  <b>PARAMETER MAPPING</b> <br>
   *  Often one wants a string parameter of a setting beeing mapped to
   *  another type, e.g. a enum constant. This means alway the same code
   *  of a strcmp inside a loop ... The functions mapStringParameter and
   *  mapStringToIndex do exactly that.
   *
   *  <br><br>
   *  <b>PATH EXPANSION</b> <br>
   *  A setting of type STRING might represent a certain path in the file
   *  system. Therefore it may contain start with the abreviation '~' for
   *  the home path of the current user or for example "~homer" the home
   *  path of user "homer". Using such a path inside your program to open
   *  a file will not work, since the '~' or "~user" is expanded by the
   *  shell, if you use it in the terminal, but there is no shell inside
   *  your programm. You have to expand it manually in your program. class
   *  SettingDataBase provides this expansion with its method expandPaths.
   */
  class SettingDataBase
  {
  public:

    //! constructor
    SettingDataBase( const bool verbose = true  );
    //! destructor
    ~SettingDataBase();

    // member access

    //! returns the first member of the settings list
    Setting* getSettings() const { return m_Settings; }
    //! returns a pointer to the setting with the passed ID string
    Setting* getSetting( const char* const IDString ) const;
    //! returns setting [i] in the list (who needs that?)
    Setting* getSetting( int index ) const;
    //! returns the number of parameters of setting "IDString"
    int      getNumParameters( const char* const IDString ) const;
    //! returns a pointer to the name of the current restriction
    char*    getSection() const { return m_activeSection; }
    //! returns true, if a restriction is currently active
    bool     isRestricted() const { return m_activeSection != NULL; }
    //! expands all paths from '~'

    /*! Inspects all settings of type Setting::STRING, whether they
     *  start with a '~'. Such a string is recognized as a path and
     *  is tried to be expanded to a full path, no longer containing
     *  the abreviation '~'. For the expansion of selected strings
     *  use methods Setting::expandPathString and Setting::expandAllPathStrings
     *  after a call of SettingDataBase::getSetting.
     *  \return true, if all paths could be expanded, otherwise false
     */
    bool expandPaths();

    //! translate a setting string into another parameter

    /*! Often one wants a string parameter of a setting beeing mapped
     *  to another type, e.g. a enum constant. This means alway the
     *  same code of a strcmp inside a loop ... The function
     *  SettingDataBase::mapStringSetting does exactly that. Code
     *  example:
     *
     *  \code
     *  SettingStringMap<int> map[] =
     *  {
     *      SettingStringMap<int>( "homer",  0 ),
     *      SettingStringMap<int>( "lisa",   6 ),
     *      SettingStringMap<int>( "bart",  -6 ),
     *      SettingStringMap<int>( )
     *  };
     *
     *  int familyvalue;
     *  settings.mapStringParameter( "family", map,
     *                               familyvalue, 1, true );
     *  \endcode
     *
     * Such a mapping often could be reduced to simply returning
     * the index of an array, the possible strings are contained
     * in (see example). To prevent the creation of SettingStringMap
     * arrays, that contain in the m_Value entry simply the index
     * in the array, this mapping is supported explicitly by the
     * function SettingDataBase::mapStringToIndex.
     *
     *  \param ID the ID string of the setting you want to map
     *  \param map an array of Type SettingStringMap<T>, that contains
     *         the mapping of type STRING to type T. Like the arrays
     *         of SettingDef this array must end with one element
     *         created with the standard constructor (to get
     *         m_String == NULL).
     *  \param value address of a variable of type T, the mapped value
     *         should be stored in
     *  \param iparameter index of the parameter string of the setting
     *         that should be mapped. In most cases these settings
     *         just take one parameter (e.g. runningmode = server),
     *         where iparameter == 0. Therefore the default value of
     *         iparameter is 0.
     *  \param verbose print some informations, if mapping did not
     *         succeed.
     *  \return true, if "value" constains valid data after the call
     */
    template <class T>
    bool mapStringParameter( const char*                const ID,
                             const SettingStringMap<T>* const map,
                             T&   value,
                             const int  iparameter = 0,
                             const bool verbose    = false ) const;
                
    //! maps a string to its index in an array of strings

    /*! maps a string to its index in an array of strings <br>
     *  Code example:
     *
     *  \code
     *  enum { HOMER = 0, MARGE, MAGGY, BART,
     *         LISA, FAMILY_SIZE );
     *  const char* enumStrings[FAMILY_SIZE + 1] = {
     *      "HOMER", "MARGE", "MAGGY",
     *      "BART",  "LISA",   NULL
     *  };
     *  settings.mapStringToIndex( "athome", enumStrings, 0, true );
     *  \endcode
     *
     *  \param ID the ID string of the setting you want to map
     *  \param SettingStrings array containing the possible strings.
     *         This array <b>must</b> be terminated with a NULL
     *         string.
     *  \param stringIndex index of the parameter string of the
     *         setting that should be mapped. In most cases these
     *         settings just take one parameter (e.g. runningmode
     *         = server), where iparameter == 0.
     *  \param verbose print some informations, if mapping did not
     *         succeed.
     *  \param stringIndex of SettingStrings, where the string could
     *         be found
     */
    int  mapStringToIndex( const char*  const ID,
                           const char** const SettingStrings,
                           const int          stringIndex = 0,
                           const bool         verbose     = false ) const;
                
    //! read settings from a file

    /*! Read settings from an configuration file
     *  \param settingdef array of SettingDef defining the accepted
     *         settings. The last element in this array <b>must</b> be
     *         one with IDString == NULL, i.e. simply created with
     *         the constructor using only default arguments.
     *  \param FileName name of the configuration file
     *  \param overwrite if true, all in the class already present
     *         settings are overwritten, if one with the same
     *         ID-String is read. This flag allows for example a
     *         hierarchy of configuration sources, that are read
     *         one after the other.
     *  \param pickySuccess: pointer to a boolean, where the success
     *         is indicated in, in a pickier way
     *  \return true, if all settings could be read
     */
    bool readSettings( const SettingDef* const settingdef,
                       const char*       const FileName,
                       const bool              overwrite    = true,
                       bool*       const pickySuccess = NULL );

    //! read settings from the command line passed as int and char**
    /*!
     *  \param settingdef array of SettingDef defining the accepted
     *         settings
     *  \param CommandLine pointer to the command line, just the
     *         first parameter passed to the main function
     *  \param nargs number of arguments in the command line, just
     *         the second parameter passed to the main function
     *  \param overwrite if true, all in the class already present
     *         settings are overwritten, if one with the same
     *         ID-String is read. This flag allows for example a
     *         hierarchy of configuration sources, that are read
     *         one after the other.
     *  \param pickySuccess pointer to a boolean, where the success
     *         is indicated in, in a pickier way
     *  \param expectDashesInCommandLine : if true, additionally,
     *         not in m_IDString present, a "-" is demanded for a
     *         setting in the command line.
     */           
    bool readSettings( const SettingDef* const settingdef,
                       const char**      const CommandLine,
                       const int               nargs,
                       const bool              overwrite = true,
                       const bool              expectDashesInCommandLine = true,
                       bool*       const pickySuccess = NULL );
                
    //! print out help using the definitions in settingdef

    /*! This function generates a standard help output for all
     *  settings in the array at settingdef. This array must be
     *  terminated by an empty SettingDef object (generated with
     *  with the standard constructor). For each setting the ID
     *  strings are printed, the storage type, limits for number of
     *  parameters, and the help string, if specified
     *  (-> SettingDef::SettingDef).
     */
    void printHelp( const SettingDef* const settingdef ) const;
                
    //! finally check the results of readSettings
    bool finalCheck( const SettingDef* const settingdef ) const;

    // set read restriction

    //! restricst reading of the config file to the passed section

    /*! \param section Zero-terminated name of the section the
     *         reading of the configuration file should be restricted
     *         to (see also <b>DIRECTIVES</b> in documentation of
     *         class SettingDataBase). If NULL is passed, the current
     *         restriction will be removed (-> removeSection)
     *  \return true, if the registration of the restriction was
     *          successfull, false, otherwise
     */
    bool restrictToSection( const char* const section = NULL );
                
    //! Removes the restriction to a section.

    /*! In the implementation this is simply a call of
     *  restrictToSection(NULL).
     */
    void removeSection() { restrictToSection(NULL); }
                
    // static function
    static void printLineError( const char* const FunctionName,
                                const char* const FileName = NULL,
                                const int         LineNumber = -1 );

  protected:
        
    // reads one line of a file into the buffer
    int readLine(       FILE*       file,
                        char**      line,
                        const char* const FileName,
                        const int         LineNumber );
    // prescan line
    int prescanLine( const char* const line,
                     const int         length,
                     const char* const FileName,
                     const int         LineNumber );
    // process the raw line read from a file
    int findMatchingSetting(       char*             line,
                                   const int               linelength,
                                   const SettingDef *const settingdefs,
                                   const char       *const FileName     = NULL,
                                   const int               LineNumber   = -1,
                                   const bool              quiet        = false,
                                   bool       *const pickySuccess = NULL );
    // find the next candidate for a setting in the command line
    int identifySettingInCommandLine( const SettingDef* const settingdef,
                                      const char**      const CommandLine,
                                      const int               i,
                                      const bool              expectDashesInCommandLine,
                                      const bool              quiet = true,
                                      bool*       const pickySuccess = NULL );
    // build a line like in a config file from the command line
    int buildConfigLineFromCommandLine( const char** const CommandLine,
                                        const int          istart,
                                        const int          iend,
                                        char**       line,
                                        const bool         expectDashesInCommandLine );
    // integrate the new setting in the list at "m_Settings"
    bool integrateSetting(       Setting *const setting,
                                 const bool           override,
                                 bool    *const pickySuccess = NULL );
    // enlarge buffer by one granularity
    char* enlargeLineBuffer( char* line );

    /////////////////
    // data members
    bool     m_verbose;
    Setting *m_Settings;
    char    *m_activeSection;
    bool     m_insideSection;             
    // current Length of the line buffer
    int      m_currentLineLength;

    // constants for direcitves in config file
    enum { DIRECTIVE_SECTION_ENTER = 0,
           DIRECTIVE_SECTION_LEAVE,
           DIRECTIVE_END,
           N_DIRECTIVES };
    // static member keeping the possible directives
    static const char* const m_directives[N_DIRECTIVES];
  };

  /**********************************************************
   * Template implementations
   **********************************************************/

  /**********************************************************
   * Automises the mapping from a string parameter to another type T
   *
   * ID      : IDString of the setting
   * map     : array of type SettingStringMap, of which the last element
   *           must be created without arguments (compare class Settingdef).
   * value   : address of variable the mapped value should be written to
   * verbose : if true, "value" contains a valid value
   **********************************************************/
 
#ifndef  FUNCTION_ERROR
#ifdef  SUPPRESS_COLORED_OUTPUT
#define  FUNCTION_ERROR(fns)  fns << "\n >> "
#else
#define  FUNCTION_ERROR(fns)  "\033[31m" << fns << ":\033[0m\n >> "
#endif
#define  FUNCTION_ERROR_WAS_NOT_DEFINED
#endif

#ifndef  CALLED_FROM
#ifdef  SUPPRESS_COLORED_OUTPUT
#define  CALLED_FROM(fns)     "    called from " << fns << "\n"
#else
#define  CALLED_FROM(fns)     "    called from \033[36m" << fns << "\033[0m\n"
#endif
#define  CALLED_FROM_WAS_NOT_DEFINED
#endif

#ifndef  NEW
#define  NEW  new(std::nothrow)
#define  NEW_WAS_NOT_DEFINED
#endif

  template <class T>
  bool SettingDataBase::
  mapStringParameter( const char*                const ID,
                      const SettingStringMap<T>* const map,
                      T&                               value,
                      const int                        iparameter,
                      bool                             verbose ) const
  {
    using namespace std;
    const char fn[]   = "SettingDataBase::translateStringSetting";
    int  length = 0;

    // create an array, that contains only the strings
    while( map[length++].m_String != NULL );
    char** const stringmap = NEW char* [length];
    if( stringmap == NULL ) {
      cerr << FUNCTION_ERROR( fn ) << "could not get enough memory\n";
      return false;
    }
    for( int i = 0; i < length; i++ )  stringmap[i] = (char*)map[i].m_String;

    int index = mapStringToIndex( ID, (const char** const)stringmap,
                                  iparameter, verbose );
    if( index < 0 ) {
      if( verbose ) cerr << CALLED_FROM( fn );
    } else {
      value = map[index].m_Value;
    }

    delete [] stringmap;
    return index >= 0;
  }

#ifdef  FUNCTION_ERROR_WAS_NOT_DEFINED
#undef  FUNCTION_ERROR
#endif

#ifdef  CALLED_FROM_WAS_NOT_DEFINED
#undef  CALLED_FROM
#endif

#ifdef  NEW_WAS_NOT_DEFINED
#undef  NEW
#endif

  /**********************************************************/

}
#endif
