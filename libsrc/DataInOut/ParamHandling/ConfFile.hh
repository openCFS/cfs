#ifndef FILE_CONFIGFILE_2002
#define FILE_CONFIGFILE_2002

#include "tools.hh"
#include <vector>

namespace CoupledField
{

class ConfFile
{
public:
  //! constructor with name of conf-file
   ConfFile(const Char* const afilename);

  //! deconstructor
  ~ConfFile(); 

  //! get from conf-file TypeVal(string, integer, double) value
  template<class TypeVal>
  void get(const std::string keyword, TypeVal & val, const std::string section="", const std::string subsection="", const std::string subsubsection="");

   //! get from conf-file TypeVal(string, integer, double) value, if this value is specified there
  template<class TypeVal>
  Boolean ifget(const std::string keyword, TypeVal & val, const std::string section="", const std::string subsection="", const std::string subsubsection="");

  //! return yes/no of keyword. if this keyword is absent in file, then FALSE
   Boolean get_option(const std::string keyword, const std::string section="", const std::string subsection = "", const std::string subsubsection = "");

   //! get type of equation for subdomain with number numsd
  void getsubdom(std::vector<std::string> & subdoms);

   //! get list of integer-values(for ex. history nodes)
   void getlist(std::vector<Integer> & hist,const std::string seekexp);

  //! get list of string-values(for ex. pdes)
  void getliststr(const std::string seekexp, std::vector<std::string> & pdes, const std::string section="", const std::string subsection="");

  //! get list of subdomains for equation, which is described in section nameSection
  void getsubdompde(std::vector<std::string> & subdoms, const std::string nameSection);

protected:

  //! name of file
  Char * filename;

  //! input file
  std::ifstream infile;  

   //! write an error and stop an execution of the program
   void error(const std::string keyword) const;

   //! get position in conf-file
   std::string::size_type getpos(const std::string keyword, const std::string::size_type startpos=0, Boolean writeErr=TRUE);

private:

  //! final position on file
  std::string::size_type pos_end;

  //! store list of all subdomains for subdomain. we need this for checking that subdomains for PDE are valid
  std::vector<std::string> allSubDomains_;

  //! this auxialary function for checking that this parameter is from predefined list
  void check(const std::string value, std::vector<std::string> data);

};

} // end of namespace

#endif
