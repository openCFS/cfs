#ifndef FILE_CONFIGFILE_2002
#define FILE_CONFIGFILE_2002

#include "Utils/tools.hh"
#include <fstream>
#include <vector>

namespace CoupledField
{

//! class for handling a config-file
class ConfFile
{
public:
  //! constructor with name of conf-file
   ConfFile(const Char* const afilename);

  //! deconstructor
  ~ConfFile(); 


  //! get from conf-file value(integer, double or string).
  /*!
	\param keyword keyword for the value in config-file
	\param val return value from config-file
	\param section name of a section in which keyword there is. can be omitted.
    	\param subsection name of a subsection of the section in which keyword there is. can be omitted.
	\param subsubsection etc. can be omitted.
   */
  template<class TypeVal>
  void get(const std::string keyword, TypeVal & val, const std::string section="", 
	   const std::string subsection="", const std::string subsubsection="");

  //! get from conf-file value and a string.
  /*!
	\param keyword keyword for the value in config-file
	\param val return value from config-file
	\param fncname name of time function file
	\param section name of a section in which keyword there is. can be omitted.
    	\param subsection name of a subsection of the section in which keyword there is. can be omitted.
	\param subsubsection etc. can be omitted.
   */
  void get2(const std::string keyword, Double & val,std::string & fncname,
	   const std::string section="", const std::string subsection="", 
	   const std::string subsubsection="");


  //! get coil data from conf-file 
  /*!
	\param acoil coild data (see environmanr.hh)
	\param section name of a section in which keyword there is. can be omitted.
    	\param subsection name of a subsection of the section in which keyword there is. can be omitted.
	\param subsubsection etc. can be omitted.
   */
  void getCoilData(coilDefStruct &acoil, const std::string section, const std::string subsection);


   //! get from conf-file value(integer, double or string), if this value is specified there
  /*!
        \param keyword keyword for the value in config-file
        \param val return value from config-file
        \param section name of a section in which keyword there is
        \param subsection name of a subsection of the section in which keyword there is
        \param subsubsection etc.
   */
  template<class TypeVal>
  Boolean ifget(const std::string keyword, TypeVal & val, const std::string section="", const std::string subsection="", const std::string subsubsection="");

  //! use only in cases when value of keyword is yes/no. if the keyword is absent in file, then return FALSE
  /*!
        \param keyword keyword for the value in config-file
        \param section name of a section in which keyword there is
        \param subsection name of a subsection of the section in which keyword t
here is
        \param subsubsection etc.
   */
   Boolean get_option(const std::string keyword, const std::string section="", const std::string subsection = "", const std::string subsubsection = "");

   //! get list of subdomains
  void getsubdom(std::vector<std::string> & subdoms);

   //! get list of integer-values for a keyword(for ex. history nodes)
   /*!
	\param hist list of returned values
	\param seekexp keyword
    */
  void getlist(std::vector<Integer> & hist,const std::string seekexp);

   //! get list of integer-values for a keyword(for ex. history nodes)
   /*!
	\param hist list of returned values
	\param seekexp keyword
	\param section name of a section in which keyword is there. can be omitted
	\param subsection name of a subsection of the section in which keyword is there. can be omitted
    */
  void getlist(const std::string seekexp, std::vector<Double> & hist, const std::string section="", const std::string subsection="");

  //! get list of string-values(for ex. pdes) for a keyword
  /*!
        \param seekexp keyword
        \param pdes vector with returned list of string
	\param section name of a section in which keyword is there. can be omitted
	\param subsection name of a subsection of the section in which keyword is there. can be omitted
  */
  void getliststr(const std::string seekexp, std::vector<std::string> & pdes, const std::string section="", const std::string subsection="");

  //! get string-value for a keyword
  /*!
    \param seekexp keyword
    \param str string
    \param section name of a section in which keyword is there. can be omitted
    \param subsection name of a subsection of the section in which keyword is there. can be omitted
  */
  void getstr(const std::string seekexp, std::string &str, const std::string section="", const std::string subsection="");

   //! get list of string-values(for ex. pdes) for a keyword, if this keyword is in config-file.
  /*!
        \param seekexp keyword
        \param pdes vector with returned list of string
        \param section name of a section in which keyword is there. can be omitted
        \param subsection name of a subsection of the section in which keyword is there. can be omitted
  */
  Boolean ifgetliststr(const std::string seekexp, std::vector<std::string> & pdes, const std::string section="", const std::string subsection="");

  //! get list of subdomains for an equation, which is defined in a section 
    /*!
	\param subdoms vector with returned values
	\param nameSection name of section in which the equation is specified
	*/
  void getsubdompde(std::vector<std::string> & subdoms, const std::string nameSection);


  /// return  name of all subdomains
  std::vector<std::string> & GetAllSubDomains() {return allSubDomains_;};
  

protected:

  //! name of file
  Char * filename;

  //! input file
  std::ifstream infile;  

  //! write an error and stop an execution of the program
  void error(const std::string keyword) const;
  
  //! get position in conf-file
  std::string::size_type getpos(const std::string keyword, 
				const std::string::size_type startpos=0, 
				Boolean inSection=FALSE,
				Boolean inSubSection=FALSE,
				Boolean writeErr=TRUE);
  
  //! get position of section in conf-file
  std::string::size_type getsectionpos(const std::string keyword, const std::string::size_type startpos=0, Boolean writeErr=TRUE);

  //! get poistion of subsection in conf-file
  std::string::size_type getsubsectionpos(const std::string keyword, const std::string::size_type startpos=0, Boolean writeErr=TRUE);
  


private:

  //! final position on file
  std::string::size_type pos_end;

  //! store list of all subdomains. we need this for checking that subdomains for PDE are valid
  std::vector<std::string> allSubDomains_;

  //! checking that the value is from predefined list
  void check(const std::string value, std::vector<std::string> data);

  //!
  void open_file();

};

} // end of namespace

#endif
