#include <string>
#include <fstream>
#include <iostream>

#include "conffile.hh"

namespace CoupledField
{

ConfFile::ConfFile(const Char* const afilename)
{
#ifdef TRACE
  (*trace) << "entering ConfFile::ConfFile" << std::endl;
#endif

 filename = new Char[100];
 strcpy(filename,afilename);

 infile.open(strcat(filename,".conf"));

 if (!infile) infile.open("general.conf");
 if (!infile) {std::cerr << "ERROR(" << __FILE__ << " " << __LINE__ <<
                         ") Can't open " << filename << std::endl;
                exit(1);}

 infile.seekg(0, std::ios::end);
 pos_end=infile.tellg();
}

ConfFile::~ConfFile()
{
 infile.close();
}

template<class TypeVal>
void ConfFile::get(const std::string keyword, TypeVal & val, const std::string section)
{
 std::string::size_type pos,pos1=0;

 if (section != "") pos1=getpos(section);

 pos=getpos(keyword,pos1);

 infile.seekg(pos, std::ios::beg);
 infile >> val;
}

std::string::size_type ConfFile::getpos(const std::string keyword,const std::string::size_type startpos)
{
  std::string::size_type help,pos=std::string::npos;
  std::string buf;

  infile.seekg(startpos, std::ios::beg);
   while ( pos == std::string::npos && !infile.eof() )
  {
    help=infile.tellg();
    std::getline(infile, buf, '\n');
    pos=buf.find(keyword);
  }

  if (pos>=pos_end) error(keyword);

  pos=buf.find("=");
  return pos+help+1;
}

#ifdef __GNUC__
template void ConfFile::get(const std::string , std::string &);
template void ConfFile::get(const std::string , Integer &);
template void ConfFile::get(const std::string , Double &);
#endif

void ConfFile::getmatnum(Integer & matnum, const Integer numsubdom)
{
 std::string::size_type pos;
 pos=getpos("numsubdomain");
 infile.seekg(pos,std::ios::beg);
 infile.ignore(100,'\n');

 Integer i=0, numsd=0;
 //
 Integer numsubdomain=numsubdom+1;
 std::string buffer;

 while (numsd != numsubdomain && i!=1000)
 { infile >> numsd >> buffer >> matnum;
   i++;
 }
 
 if (i==1000) error("Error in conf->getmatnum: this numsubdom is absent in config-file");
}

void ConfFile::getequation(std::string & eq, const Integer numsubdom)
{
 std::string::size_type pos;
 pos=getpos("numsubdomain");
 infile.seekg(pos,std::ios::beg);
 infile.ignore(100,'\n');

 Integer numsd=0, i=0;
 Integer numsubdomain=numsubdom+1; 

  while (numsd != numsubdomain && i!=1000)
 {
   infile >> numsd >> eq;
   i++;
 }

 if (i==1000) error("Error in conf->getequation: this numsubdom is absent in config-file");
}

void ConfFile::error(const std::string keyword) const
{
std::cerr << "\033[32m ERROR: \033[0m (" << __FILE__ <<" "<< __LINE__
               << ") Cannot find string: " << keyword ;
          std::cerr << " in your conf-file.\n\t\t Please, check conf-file."<< std::endl;
                      exit(1);
}

} // end of namespace
