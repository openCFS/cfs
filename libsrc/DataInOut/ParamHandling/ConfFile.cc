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
#ifdef TRACE
  (*trace) << " entering ConfFile::~ConfFile " << std::endl;
#endif

  delete [] filename;

 infile.close();
}

template<class TypeVal>
void ConfFile::get(const std::string keyword, TypeVal & val, const std::string section, const std::string subsection, const std::string subsubsection)
{
 std::string::size_type pos,pos1=0;

 if (section != "") pos1=getpos(section); 
 if (subsection !="") pos1=getpos(subsection,pos1);
 if (subsubsection != "") { pos1=getpos(subsubsection,pos1);
                            infile.seekg(pos1, std::ios::beg);
                            infile.ignore(100,'\n');                      
                            pos1=infile.tellg();
                          }

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

void ConfFile::getsubdom(std::vector<std::string> & subdoms)
{
 std::string::size_type pos;

 Integer nsubds;
 get("subdomains",nsubds);

 pos=getpos("list_subdomains");
 infile.seekg(pos,std::ios::beg);
 infile.ignore(100,'\n');

 subdoms.resize(nsubds);

 Integer i;
 for (i=0; i < nsubds; i++)
 {
   infile >> subdoms[i];
   infile.ignore(100,'\n');
 }
}

void ConfFile::getsubdompde(std::vector<std::string> & subdoms, const std::string section)
{
  std::string::size_type pos=0;

  pos=getpos(section,pos);
  pos=getpos("subdomains",pos);
  infile.seekg(pos,std::ios::beg);
  
  std::string help;
  do
  {
   infile >> help;
   if (help!="non") subdoms.push_back(help);   
  } while(help!="non");
}

void ConfFile::getlist(std::vector<Integer> & hist, const std::string seekexp)
{
 std::string::size_type pos;
 pos=getpos(seekexp);
 infile.seekg(pos,std::ios::beg);

 Integer node;
 do
{
 infile >> node;
 if (node!=-1) hist.push_back(node);
}
 while (node!=-1);
}

void ConfFile :: getliststr( const std::string seekexp, std::vector<std::string> & stlist, const std::string section, const std::string subsection)
{
 std::string::size_type pos=0;

 if (section != "") pos=getpos(section);
 if (subsection !="") pos=getpos(subsection,pos);

 pos=getpos(seekexp,pos);

 infile.seekg(pos,std::ios::beg);

 std::string help;
 do
 {
  infile >> help;

  if (help != "non") stlist.push_back(help);
 } while  (help != "non");

}

void ConfFile::error(const std::string keyword) const
{
std::cerr << "\033[32m ERROR: \033[0m (" << __FILE__ <<" "<< __LINE__
               << ") Cannot find string: " << keyword ;
          std::cerr << " in your conf-file.\n\t\t Please, check conf-file."<< std::endl;
                      exit(1);
}

} // end of namespace
