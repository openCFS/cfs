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
 strcat(filename,".conf");

 open_file();

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
 Boolean inSection = FALSE;
 Boolean inSubSection = FALSE;
 if (section != "") 
   { 
     pos1=getsectionpos(section); 
     inSection=TRUE;
   }
 if (subsection !="") 
   {
     pos1=getsubsectionpos(subsection,pos1);
     inSubSection = TRUE;
   }
 if (subsubsection != "") { pos1=getpos(subsubsection,pos1);
                            infile.seekg(pos1, std::ios::beg);
                            infile.ignore(100,'\n');                      
                            pos1=infile.tellg();
                          }

 pos=getpos(keyword,pos1,inSection,inSubSection);

 infile.seekg(pos, std::ios::beg);
 infile >> val;
}

template<class TypeVal>
Boolean ConfFile::ifget(const std::string keyword, TypeVal & val, const std::string section, const std::string subsection, const std::string subsubsection)
{
 std::string::size_type pos,pos1=0;
 Boolean inSection = FALSE;
 Boolean inSubSection = FALSE;

 if (section != "") 
   { 
     pos1=getsectionpos(section,0, FALSE); 
     inSection=TRUE;
   }
 if (subsection !="") 
   {
     pos1=getsubsectionpos(subsection,pos1,FALSE);
     inSubSection = TRUE;
   }

 if (subsubsection != "") { pos1=getpos(subsubsection,pos1);
                            infile.seekg(pos1, std::ios::beg);
                            infile.ignore(100,'\n');                      
                            pos1=infile.tellg();
                          }

 pos=getpos(keyword,pos1,inSection,inSubSection,FALSE);

 if (pos==std::string::npos) return FALSE;

 infile.seekg(pos, std::ios::beg);
 infile >> val;
 
 return TRUE;  
}

Boolean ConfFile::get_option(const std::string keyword, const std::string section, const std::string subsection, const std::string subsubsection)
{
#ifdef TRACE
  (*trace) << " entering ConfFile::get_option " << std::endl;
#endif

 std::string::size_type pos,pos1=0;
 Boolean inSection = FALSE;
 Boolean inSubSection = FALSE;

 if (section != "") 
   { 
     pos1=getsectionpos(section); 
     inSection=TRUE;
   }
 if (subsection !="") 
   {
     pos1=getsubsectionpos(subsection,pos1);
     inSubSection = TRUE;
   }


 if (subsubsection != "") { pos1=getpos(subsubsection,pos1);
                            infile.seekg(pos1, std::ios::beg);
                            infile.ignore(100,'\n');                      
                            pos1=infile.tellg();
                          }

 pos=getpos(keyword,pos1,inSection,inSubSection,FALSE);

 Boolean val=FALSE;

 if (pos==std::string::npos) return FALSE;

 infile.seekg(pos, std::ios::beg);
 std::string option;
 infile >> option;
 
 if (option == "yes") return TRUE;

 return FALSE;
}

std::string::size_type ConfFile::getpos(const std::string keyword,
					const std::string::size_type startpos,
					Boolean inSection,
					Boolean inSubSection,
					Boolean writeErr)
{
  std::string::size_type help,pos=std::string::npos;
  std::string buf;
  bool nextSectionReached = FALSE;

  if (infile.eof()) infile.clear();
  infile.seekg(startpos, std::ios::beg);
 
     while ( pos == std::string::npos && !infile.eof() )
  {
    help=infile.tellg();
    std::getline(infile, buf, '\n');

    if (buf[0]!='#' && buf.size() != 0) 
      {
	if ((buf[0] != '\t' && inSection) && buf[0] != '\n' && buf[0] != ' ')
	  {
	    nextSectionReached = TRUE;
	    break;
	  }
	if (buf[0] != '\t' && buf[1] != '\t' && inSubSection && buf[0] != '\n' && buf[0] != ' ')
	  {
	    nextSectionReached = TRUE;
	    break;
	  }
	pos=buf.find(keyword);
      }
    
  }
   
   if (pos>=pos_end || nextSectionReached)  {
     if (writeErr) error(keyword);
     return std::string::npos;
   }
   
   pos=buf.find("=");
   return pos+help+1;
      
}

std::string::size_type ConfFile::getsectionpos(const std::string keyword, const std::string::size_type startpos, Boolean writeErr)
{
  std::string::size_type help,pos=std::string::npos;
  std::string buf;
  std::string::size_type pos_helper;

  if (infile.eof()) infile.clear();
  infile.seekg(startpos, std::ios::beg);
    
   while ( pos == std::string::npos && !infile.eof() )
  {
    help=infile.tellg();
    std::getline(infile, buf, '\n');
 
    if (buf[0] !='#' && buf[0] != ' ' && buf[0] != '\t')
      {
	pos_helper=buf.find(keyword);
	if (buf.find(":") < 100)
	  pos = pos_helper;
      }
  }
   
   if (pos>=pos_end)  {
     if (writeErr) error(keyword);
        return std::string::npos;
   }

  return infile.tellg();

}

std::string::size_type ConfFile::getsubsectionpos(const std::string keyword, const std::string::size_type startpos, Boolean writeErr)
{
  std::string::size_type help,pos=std::string::npos;
  std::string buf;
  std::string::size_type pos_helper;
  bool nextSectionReached = FALSE;

  if (infile.eof()) infile.clear();
  infile.seekg(startpos, std::ios::beg);
    
   while ( pos == std::string::npos && !infile.eof() )
  {
    help=infile.tellg();
    std::getline(infile, buf, '\n');
 
    if (buf[0]!='#' && buf.size() != 0) 
      {
	if (buf[0] != '\t' && buf[0] != '\n' && buf[0] != ' ')
	  {
	    nextSectionReached = TRUE;
	    break;
	  }
	
	pos_helper=buf.find(keyword);
	if (buf.find(":") < 100)
	  pos = pos_helper;
      }
  }
   
   if (pos>=pos_end && nextSectionReached)  {
     if (writeErr) error(keyword);
        return std::string::npos;
   }

  
   //return pos+help+1;
   return infile.tellg();

}


#ifdef __GNUC__
template void ConfFile::get(const std::string , std::string &, const std::string, const std::string, const std::string);
template void ConfFile::get(const std::string , Integer &, const std::string, const std::string, const std::string);
template void ConfFile::get(const std::string , Double &,const std::string, const std::string, const std::string);

template Boolean ConfFile::ifget(const std::string , std::string &, const std::string, const std::string, const std::string);
template Boolean ConfFile::ifget(const std::string , Integer &,const std::string, const std::string, const std::string);
template Boolean ConfFile::ifget(const std::string , Double &, const std::string, const std::string, const std::string);
#endif

void ConfFile::getsubdom(std::vector<std::string> & subdoms)
{
 std::string::size_type pos;

 Integer nsubds;
 get("subdomains",nsubds,"","","");

 pos=getpos("list_subdomains");
 infile.seekg(pos,std::ios::beg);
 infile.ignore(100,'\n');

 subdoms.resize(nsubds);

 Integer i;
 for (i=0; i < nsubds; i++)
 {
   infile >> subdoms[i];
   allSubDomains_.push_back(subdoms[i]);
   infile.ignore(100,'\n');
 }
}

void ConfFile::getsubdompde(std::vector<std::string> & subdoms, const std::string section)
{
  std::string::size_type pos=0;

  pos=getsectionpos(section,pos);
  pos=getpos("subdomains",pos);
  infile.seekg(pos,std::ios::beg);
  
  std::string help;
  do
  {
   infile >> help;
   if (help!="non") {
     check(help,allSubDomains_);
     subdoms.push_back(help);   
   }
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



void ConfFile::getlist( const std::string seekexp, std::vector<Double> & hist,
			const std::string section, const std::string subsection)
{
 std::string::size_type pos=0;

 Boolean inSection = FALSE;
 Boolean inSubSection = FALSE;

 if (section != "")
   {
     pos=getsectionpos(section);
     inSection = TRUE;
   }
 if (subsection !="") 
   {   
     pos=getsubsectionpos(subsection,pos);
     inSubSection = TRUE;
   }

 pos=getpos(seekexp,pos);
 infile.seekg(pos,std::ios::beg);

 
 Double help;
 do
   {
     infile >> help;
     
     if (help != 1e100) hist.push_back(help);
   }
 while (help != 1e100);
}





void ConfFile::getliststr( const std::string seekexp, std::vector<std::string> & stlist, const std::string section, const std::string subsection)
{
 std::string::size_type pos=0;

 Boolean inSection = FALSE;
 Boolean inSubSection = FALSE;

 if (section != "")
   {
     pos=getsectionpos(section);
     inSection = TRUE;
   }
 if (subsection !="") 
   {   
     pos=getsubsectionpos(subsection,pos);
     inSubSection = TRUE;
   }

 pos=getpos(seekexp,pos,inSection,inSubSection);

 infile.seekg(pos,std::ios::beg);

 std::string help;
 do
 {
  infile >> help;

  if (help != "non") stlist.push_back(help);
 } while  (help != "non");
}

void ConfFile::getstr( const std::string seekexp, std::string &str, const std::string section, const std::string subsection)
{
 std::string::size_type pos=0;
 Boolean inSection = FALSE;
 Boolean inSubSection = FALSE;

 if (section != "")
   {
     pos=getsectionpos(section,0,TRUE);
     inSection = TRUE;
   }
 if (subsection !="") 
   {   
     pos=getsubsectionpos(subsection,pos);
     inSubSection = TRUE;
   }

 pos=getpos(seekexp,pos,inSection,inSubSection);

 infile.seekg(pos,std::ios::beg);

 std::string help;
 infile >> str;

}

Boolean ConfFile::ifgetliststr( const std::string seekexp, std::vector<std::string> & stlist, const std::string section, const std::string subsection)
{
 std::string::size_type pos=0;
 Boolean inSection = FALSE;
 Boolean inSubSection = FALSE;

 if (section != "")
   {
     pos=getsectionpos(section,0,FALSE);
     inSection = TRUE;
   }
 if (subsection !="") 
   {   
     pos=getsubsectionpos(subsection,pos,FALSE);
     inSubSection = TRUE;
   }

 pos=getpos(seekexp,pos,inSection,inSubSection,FALSE);

 if (pos==std::string::npos) return FALSE;

 infile.seekg(pos,std::ios::beg);

 std::string help;
 do
 {
  infile >> help;

  if (help != "non") stlist.push_back(help);
 } while  (help != "non");

 return TRUE;
}

void ConfFile::open_file()
{

  if (infile.is_open()) infile.close(); 
  infile.clear();
  infile.open(filename);

 if (!infile) infile.open("general.conf");
 if (!infile) {std::cerr << "ERROR(" << __FILE__ << " " << __LINE__ <<
                         ") Can't open " << filename << std::endl;
                exit(1);}

 infile.seekg(0, std::ios::end);
}

void ConfFile::error(const std::string keyword) const
{
std::cerr << "\033[32m ERROR: \033[0m (" << __FILE__ <<" "<< __LINE__
               << ") Cannot find string: " << keyword ;
          std::cerr << " in your conf-file.\n\t\t Please, check conf-file."<< std::endl;
                      exit(1);
}

void ConfFile::check(const std::string value, const std::vector<std::string> data)
{
  Boolean Find=FALSE;
  Integer id;
  for (id=0; id<data.size();id++) {
    if (value==data[id]) {
      Find=TRUE;
      break;
    }
  }

  if (!Find) {
    std::cerr << " \033[32m ERROR: \033[0m (" << __FILE__ <<" "<< __LINE__
               << ") This subdomain: " << value ;
          std::cerr << " is not specified in list of all subdomains.\n\t\t Please, check conf-file."<< std::endl;
                      exit(1);
  }		      

}

} // end of namespace
