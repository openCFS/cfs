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
  std::cout << " dim " << dim_ << std::endl;
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
  getPosition("Dimension", pos);
  infile.seekg(pos,std::ios::beg);

  std::string auxname;
  infile >> dim_;
}

void AnsysFile::ReadCoordinate(Point2D * const NodesCoord, const Integer maxnumnodes)
{
#ifdef TRACE
  (*trace) << "entering Ansys::ReadCoordinate 2D" << std::endl;
#endif

  Integer i;
  std::string::size_type pos=0;
  getPosLine("Nodes", pos);
  infile.seekg(pos,std::ios::beg);

  Integer ibuf;
  Double dbuf;

  for (i=0; i < maxnumnodes; i++)
  {
    infile >> ibuf >> dbuf >> NodesCoord[i].x >> NodesCoord[i].y ;
    infile.ignore(100,'\n');
  }
}

void AnsysFile::ReadCoordinate(Point3D * const NodesCoord, const Integer maxnumnodes)
{
#ifdef TRACE
  (*trace) << "entering Ansys::ReadCoordinate 3D" << std::endl;
#endif

  Integer i;
  std::string::size_type pos=0;
  getPosLine("Nodes", pos);
  infile.seekg(pos,std::ios::beg);

  Integer ibuf;
  for (i=0; i < maxnumnodes; i++)
  {
    infile >> ibuf >> NodesCoord[i].x >> NodesCoord[i].y >>  NodesCoord[i].z;
    infile.ignore(100,'\n');
  }
}

void AnsysFile::ReadMaxnumnodes(Integer & nnodes)
{
#ifdef TRACE
  (*trace) << "entering Ansys::ReadMaxnumnodes" << std::endl;
#endif

  std::string::size_type pos=0;
  getPosition("NumNodes", pos);
  infile.seekg(pos,std::ios::beg);

  infile >> nnodes;
  std::cout << nnodes << std::endl;
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
  keyword="Num2DElements";
  break;
case 3:
  keyword="Num3DElements";
  break;
default:
  Error("Wrong dimension in function ReadMaxnumelem");
}

 std::string::size_type pos=0;
 getPosition(keyword,pos);
 infile.seekg(pos,std::ios::beg);

 infile >> nelem;
 std::cout << nelem << std::endl;
}

void AnsysFile::ReadMaxnumnodesbc(Integer & nbc)
{
 std::string::size_type pos=0;
 getPosition("NumNodeBC", pos);
 infile.seekg(pos,std::ios::beg);

 infile >> nbc;
}

void AnsysFile::ReadBoundRestr(std::list<NodeRestraint> & restr, Integer & numberRestr)
{
#ifdef TRACE
  (*trace) << "entering Ansys::ReadBoundRestr" << std::endl;
#endif
  
 Integer numbc;
 ReadMaxnumnodesbc(numbc);
 std::cout << numbc << std::endl;

 std::string::size_type pos=0;
 getPosLine("Node BC", pos);
 infile.seekg(pos,std::ios::beg);

 Integer nrestr=0;
 std::string str; 

 NodeRestraint A;
 Integer i;
 for (i=0; i < numbc; i++)
   {
    infile >> A.nodalnum >> str;
    std::cout << str << std::endl;

    A.dof=TransformInDof(str);
    infile.ignore(100,'\n');

    if (A.dof==5)
    { restr.push_back(A); nrestr++;}
   }
  numberRestr=nrestr;
}

void AnsysFile::getPosLine(const std::string seekexp, std::string::size_type & pos)
{
  infile.seekg(pos, std::ios::beg);
  std::string buf;
  pos=std::string::npos;

  std::string::size_type hpos;

  std::cout << seekexp << std::endl;

  while (pos == std::string::npos && !infile.eof())
{
  hpos=infile.tellg();
  std::getline(infile, buf, '\n');
  std::cout << hpos << buf << std::endl;
  pos=buf.find(seekexp);
}

 pos=infile.tellg();

// check, if there are comments lines
 do
{
  std::getline(infile, buf, '\n');
  if (buf[0]=='#') pos=infile.tellg();
  std::cout << buf << std::endl;
}
 while (buf[0] != '#'); 

 if (pos>=pos_end) {std::cerr << "ERROR: (" << __FILE__ <<" "<< __LINE__
               << ") Cannot find string: " << seekexp ;
          std::cerr << " in your dat file.\n\t\t Please, change your dat file."<< std::endl;
                      exit(1);}
}

void AnsysFile::getPosition(const std::string seekexp, std::string::size_type & pos)
{
  infile.seekg(pos, std::ios::beg);
  std::string buf;
  pos=std::string::npos;

  std::string::size_type hpos;

  std::cout << seekexp << std::endl;

  while (pos == std::string::npos && !infile.eof())
{
  hpos=infile.tellg();
  std::getline(infile, buf, '\n');
  std::cout << hpos << buf << std::endl;
  pos=buf.find(seekexp);
}
  pos+=hpos+seekexp.length();

 if (pos>=pos_end) {std::cerr << "ERROR: (" << __FILE__ <<" "<< __LINE__
               << ") Cannot find string: " << seekexp ;
          std::cerr << " in your dat file.\n\t\t Please, change your dat file."<< std::endl;
                      exit(1);}
}

Integer AnsysFile::TransformInDof(const std::string type_bc)
{
  Integer result;
  std::cout << " tb " << type_bc << std::endl;

  if (type_bc=="vp-restrataint") result=5;
  else
  if (type_bc=="ep-restrataint") result=4;
  else Error(" This type of level for boundary condition in mesh-file (section [Node BC]) is unknown");

  return result;
}

void AnsysFile::ReadNumberNodesPerElem(Integer & nnodesperelem)
{
 std::string::size_type pos=0;
 getPosLine("2D Elements", pos);
 infile.seekg(pos,std::ios::beg);

 std::string buf;
 std::getline(infile, buf, '\n');

 Integer ibuf;
 infile >> ibuf >> ibuf >> nnodesperelem;  
}

  //!
void AnsysFile::ReadElemConnectionGH(const Integer maxelem, Integer * connect, const Integer maxnode, const Integer numelemgr, const Integer startposinarrayconn)
{
 std::string::size_type pos=0;
 getPosLine("2D Elements", pos);
 infile.seekg(pos,std::ios::beg);

 std::string buf;
 std::getline(infile, buf, '\n');
 std::getline(infile, buf, '\n');

 Integer counter=0;

 Integer j, jj;

 std::cout << " me " << maxelem << std::endl;
 for (j=0; j < maxelem; j++)
{
  for (jj=0; jj < maxnode; jj++, counter++)
   infile >> connect[counter];

  infile.ignore(100,'\n');
  std::getline(infile, buf, '\n');
}

 Integer ii;
 for(ii=0; ii < counter; ii++)
  std::cout << connect[ii] << " " << std::endl;

}

}
