#include <string>
#include <iostream>
#include <fstream>
#include <stdio.h>

#include "material.hh"

namespace CoupledField
{

 Material::Material(const Char * const aname)
{
  Char * name = new Char[100];
  strcpy(name,aname);

  infile.open(strcat(name,".dat"));
  if (!infile) {std::cerr << "ERROR(" << __FILE__ << " " << __LINE__ <<
                         ") Can't open " << name << "/ We working only with dat format / " << std::endl;
                exit(1);}
}

  Material :: ~Material()
{
 infile.close();
}

void Material::ReadDensityAndCompressity(Double & density, Double & compress, const Integer matnum, const std::string keyword)
{
#ifdef TRACE
  if (trace) (*trace) << "Entering  ReadMaterial::ReadDensityAndCompress" << std::endl;
#endif

  infile.seekg(0, std::ios::beg);  
  std::string buf;
  std::string::size_type pos=std::string::npos, pos1=std::string::npos;

 // transfer matnum in string
  Char s[20];
  sprintf(s,"%i",matnum);

 while ( (pos == std::string::npos || pos1== std::string::npos) & !infile.eof() )
  { std::getline(infile, buf, '\n');
    pos=buf.find(keyword);
    pos1=buf.find(s);
  }
  
  pos=std::string::npos; 

  while ( pos == std::string::npos & !infile.eof() )
  { std::getline(infile, buf, '\n');
    pos=buf.find("density");
  }

  if (infile.eof()) Error("Can't find data in material dat-file. Check your dat-file",__FILE__,__LINE__);

  infile >> density >> compress;
}

void Material::ReadDielectricTerms(Double & dielectr,const Integer matnum)
{
#ifdef TRACE
  if (trace) (*trace) << "Entering  Material::ReadDielectricTerms" << std::endl;
#endif

  infile.seekg(0, std::ios::beg);
  std::string buf;
  std::string::size_type pos=std::string::npos, pos1=std::string::npos;

 // transfer matnum in string
  Char s[20];
  sprintf(s,"%i",matnum);

  while ( (pos == std::string::npos || pos1 == std::string::npos) & !infile.eof())
  { std::getline(infile, buf, '\n');
    pos=buf.find("piezo");
    pos1=buf.find(s);
  }

  pos=std::string::npos;

  while ( pos == std::string::npos & !infile.eof() )
  { std::getline(infile, buf, '\n');
    pos=buf.find("dielectric");
  }

  std::getline(infile, buf, '\n');
  std::getline(infile, buf, '\n');

  if (infile.eof()) Error("Can't find data in material dat-file. Check your dat-file",__FILE__,__LINE__);

  Double dummy;
  infile >> dummy >> dummy >> dielectr;
}

} // end of namespace
