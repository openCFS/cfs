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

void ConfFile::get(const std::string keyword, std::string & val)
{
  infile.seekg(0, std::ios::beg);
   std::string buf;
  std::string::size_type pos=std::string::npos,pos1;

  while ( pos == std::string::npos && !infile.eof() )
  {
    pos1=infile.tellg();
    std::getline(infile, buf, '\n');
    pos=buf.find(keyword);
  }

   pos=buf.find("=");
   infile.seekg(pos1+pos+1, std::ios::beg);
   infile >> val;

  if (pos>=pos_end) error(keyword); 
}

void ConfFile::get(const std::string keyword, Double & val)
{
  infile.seekg(0, std::ios::beg);
  std::string buf;
  std::string::size_type pos=std::string::npos,pos1;

  while ( pos == std::string::npos && !infile.eof() )
  {
    pos1=infile.tellg();
    std::getline(infile, buf, '\n');
    pos=buf.find(keyword);
  }

   pos=buf.find("=");
   infile.seekg(pos1+pos+1, std::ios::beg);
   infile >> val;

  if (pos>=pos_end) error(keyword);
}

void ConfFile::get(const std::string keyword, Integer & val)
{
  infile.seekg(0, std::ios::beg);
  std::string buf;
  std::string::size_type pos=std::string::npos,pos1;

  while ( pos == std::string::npos && !infile.eof() )
  {
    pos1=infile.tellg();
    std::getline(infile, buf, '\n');
    pos=buf.find(keyword);
  }

   pos=buf.find("=");
   infile.seekg(pos1+pos+1, std::ios::beg);
   infile >> val;

  if (pos>=pos_end) error(keyword);
}

void ConfFile::error(const std::string keyword) const
{
std::cerr << "\033[32m ERROR: \033[0m (" << __FILE__ <<" "<< __LINE__
               << ") Cannot find string: " << keyword ;
          std::cerr << " in your conf-file.\n\t\t Please, check conf-file."<< std::endl;
                      exit(1);
}

} // end of namespace
