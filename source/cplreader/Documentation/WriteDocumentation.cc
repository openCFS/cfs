#include <fstream>

#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/convenience.hpp>
#include <boost/filesystem/exception.hpp>

namespace fs=boost::filesystem;

#include "General/exception.hh"
#include "GetDocumentation.hh"
#include "WriteDocumentation.hh"

namespace CoupledField
{
  void WriteDocumentation()
  {
    std::map<std::string, std::vector<unsigned char> > doc;
    std::map<std::string, std::vector<unsigned char> >::const_iterator it, end;
    std::stringstream sstr;
    std::string pathsep;
    std::string fn;
    static const std::string docDir = "cplreader_docs";
    
    // create directory for HTML help
    try {
      fs::create_directory( docDir );
      pathsep = fs::path("/").native_directory_string();
    } catch (std::exception &ex) {
      EXCEPTION(ex.what());
    }

    GetDocumentation(doc);
    
    for( it=doc.begin(), end = doc.end(); it != end; it++) 
    {
      sstr.clear(); sstr.str("");
      sstr << docDir << pathsep << it->first;
      fn = sstr.str();
      
      std::ofstream docFile(fn.c_str(), std::ios::binary | std::ios::trunc);
      docFile.write((const char*) &it->second[0], it->second.size());
      docFile.close();
    }

  }
  
} // end of namespace

