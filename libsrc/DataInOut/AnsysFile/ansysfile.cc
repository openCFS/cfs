#include <fstream>
#include <iostream>
#include <string>
#include <stdarg.h>

#include "ansysfile.hh"
#include "bcs.hh"

namespace CoupledField
{

AnsysFile :: AnsysFile(const Char * const afilename)
 :FileType(afilename)
{
#ifdef TRACE
  (*trace) << "entering AnsysFile::AnsysFile" << std::endl;
#endif

  infile.open(strcat(filename,".mesh"));
  if (!infile.good()) {std::cerr << "ERROR(" << __FILE__ << " " << __LINE__ <<
                         ") Can't open " << filename << std::endl;
                exit(1);}

  infile.seekg(0, std::ios::end);
  pos_end=infile.tellg();

  ReadDim();
}
  
AnsysFile :: ~AnsysFile()
{
#ifdef TRACE
  (*trace) << "entering AnsysFile::~AnsysFile" << std::endl;
#endif

 infile.close() ;
}

void AnsysFile::ReadDim()
{
#ifdef TRACE
  (*trace) << "entering Ansys::ReadDim" << std::endl;
#endif

  Integer i;
  std::string::size_type pos=0;
  takePosition("Dimension", pos);
  infile.seekg(pos,std::ios::beg);

  infile >> dim_;
}

void AnsysFile::ReadCoordinate(Point2D * const NodesCoord, const Integer maxnumnodes)
{
#ifdef TRACE
  (*trace) << "entering Ansys::ReadCoordinate 2D" << std::endl;
#endif

  Integer i;
  std::string::size_type pos=0;
  takePosition("Nodes", pos);
  infile.seekg(pos,std::ios::beg);

  Integer ibuf, dbuf;
  for (i=0; i < maxnumnodes; i++)
  {
    infile >> ibuf >> dbuf >> NodesCoord[i].x >> NodesCoord[i].y ;
  }
}

void AnsysFile::ReadCoordinate(Point3D * const NodesCoord, const Integer maxnumnodes)
{
#ifdef TRACE
  (*trace) << "entering Ansys::ReadCoordinate 3D" << std::endl;
#endif

  Integer i;
  std::string::size_type pos=0;
  takePosition("Nodes", pos);
  infile.seekg(pos,std::ios::beg);

  Integer ibuf;
  for (i=0; i < maxnumnodes; i++)
  {
    infile >> ibuf >> NodesCoord[i].x >> NodesCoord[i].y >>  NodesCoord[i].z;
  }
}

void AnsysFile::ReadMaxnumnodes(Integer & nnodes)
{
#ifdef TRACE
  (*trace) << "entering Ansys::ReadMaxnumnodes" << std::endl;
#endif

  std::string::size_type pos=0;
  takePosition("NumNodes", pos);
  infile.seekg(pos,std::ios::beg);

  infile >> nnodes;
}

void AnsysFile::ReadMaxnumelem(Integer & nelem)
{
#ifdef TRACE
  (*trace) << "entering Ansys::ReadMaxnumelem" << std::endl;
#endif

std::string keyword;
switch (dim_)
{
case 2:
  keyword="Num3DElements";
  break;
case 3:
  keyword="Num2DElements";
  break;
}

 std::string::size_type pos=0;
 takePosition(keyword, pos);
 infile.seekg(pos,std::ios::beg);

 infile >> nelem;
}

void AnsysFile::ReadMaxnumnodesbc(Integer & nbc)
{
 std::string::size_type pos=0;
 takePosition("NumNodeBC", pos);
 infile.seekg(pos,std::ios::beg);

 infile >> nbc;
}

void AnsysFile::ReadBoundRestr(std::list<NodeRestraint> & restr, Integer & numberRestr)
{
#ifdef TRACE
  (*trace) << "entering Ansys::ReadBoundRestr" << std::endl;
#endif
  
 std::string::size_type pos=0;
 takePosition("Node BC", pos);
 infile.seekg(pos,std::ios::beg);

 Integer numbc;
 ReadMaxnumnodesbc(numbc);

 Integer nrestr=0;
 std::string str; 

 NodeRestraint A;
 Integer i;
 for (i=0; i < numbc; i++)
   {
    infile >> A.nodalnum >> str;
    A.dof=TransformInDof(str.c_str());
    infile.ignore(100,'\n');

    if (A.dof==5)
    { restr.push_back(A); nrestr++;}
   }

  numberRestr=nrestr;

}

void AnsysFile::takePosition(const std::string seekexp, std::string::size_type & pos)
{
  infile.seekg(pos, std::ios::beg);
  std::string buf;
  pos=std::string::npos;

  while (pos == std::string::npos && !infile.eof())
{
  std::getline(infile, buf, '\n');
  pos=buf.find(seekexp);
}
  pos=infile.tellg();

 do
{
  std::getline(infile, buf, '\n');
  if (buf[0]=='#') pos=infile.tellg();
}
 while (buf[0] != '#'); 

 if (pos>=pos_end) {std::cerr << "ERROR: (" << __FILE__ <<" "<< __LINE__
               << ") Cannot find string: " << seekexp ;
          std::cerr << " in your dat file.\n\t\t Please, change your dat file."<< std::endl;
                      exit(1);}
}

Integer AnsysFile::TransformInDof(const Char * el)
{
  Integer result;
  if (!strcmp("vp-restrataint", el)) result=5;
  if (!strcmp("ep-restrataint",el)) result=4;
  else Error(" This type of level for boundary condition is unknown");

  return result; 
}

}
