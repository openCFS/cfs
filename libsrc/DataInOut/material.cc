#include <string>
#include <iostream>
#include <fstream>

#include <general_head.hh> 
#include "tools.hh"

#include "material.hh"
using namespace std;
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
 void Material::ReadDensityAndCompress(Double & density, Double & compress)
{
#ifdef TRACE
  (*trace) << "entering  ReadMaterial::ReadDensityAndCompress" << endl;
#endif
  infile.seekg(0, std::ios::beg);  
  string buf;
  string::size_type pos=string::npos;

  while ( pos == string::npos & !infile.eof() )
  { std::getline(infile, buf, '\n');
    pos=buf.find("density");
  }

  infile >> density >> compress;
}
} // end of namespace
