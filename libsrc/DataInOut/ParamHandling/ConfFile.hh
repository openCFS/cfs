#ifndef FILE_CONFIGFILE_2002
#define FILE_CONFIGFILE_2002

#include "tools.hh"

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
  void get(const std::string keyword, TypeVal & val);

protected:

  //! name of file
  Char * filename;

  //! input file
  std::ifstream infile;  

   //! write warning
   void error(const std::string keyword) const;

private:

  //! final position on file
  std::string::size_type pos_end;

};

} // end of namespace

#endif
