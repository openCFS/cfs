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

  //! get from conf-file Double value  
  void get(const std::string keyword, Double & val);
 
  //! get from conf-file Integer value
  void get(const std::string keyword, Integer & val);

  //! get from conf-file string value
  void get(const std::string keyword, std::string & val);

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
