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

   //! get type of equation for subdomain with number numsd
  void getsubdom(std::vector<std::string> & subdoms);

   //! get list of value like history nodes
   void getlist(std::vector<Integer> & hist,const std::string seekexp);

  //! get arrays of pdes
  void getliststr(const std::string seekexp, std::vector<std::string> & pdes, const std::string section="", const std::string subsection="");

  //!
void getsubdompde(std::vector<std::string> & subdoms, const std::string section);

protected:

  //! name of file
  Char * filename;

  //! input file
  std::ifstream infile;  

   //! write warning
   void error(const std::string keyword) const;

   //! get position in conf-file
   std::string::size_type getpos(const std::string keyword, const std::string::size_type startpos=0);

private:

  //! final position on file
  std::string::size_type pos_end;

};

} // end of namespace

#endif
