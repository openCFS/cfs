#include <fstream>
#include <iostream>
#include <string>
#include <stdarg.h>
#include <vector>

#include "ansysfile.hh"
#include "bcs.hh"
#include "grid_cfs.hh"

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
 getPosLine("[Nodes]", pos);
 infile.seekg(pos,std::ios::beg);

 Integer ibuf;

 for (i=0; i < maxnumnodes; i++)
  {
    infile >> ibuf >> NodesCoord[i].x >> NodesCoord[i].y ;
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
 std::cout << " number nodes" << nnodes << std::endl;
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
 std::cout << " EELEM " << nelem << std::endl;
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
    infile.ignore(100,'\n');

    A.dof=TransformInDof(str);

    restr.push_back(A);
    nrestr++;
   }
  numberRestr=nrestr;
}

void AnsysFile::getPosLine(const std::string seekexp, std::string::size_type & pos)
{
  infile.seekg(pos, std::ios::beg);
  std::string buf;
  pos=std::string::npos;

  std::string::size_type hpos;

  while (pos == std::string::npos && !infile.eof())
{
  hpos=infile.tellg();
  std::getline(infile, buf, '\n');
  pos=buf.find(seekexp);
}

 pos=infile.tellg();

 if (pos>=pos_end) {std::cerr << "ERROR: (" << __FILE__ <<" "<< __LINE__
               << ") Cannot find string: " << seekexp ;
          std::cerr << " in your dat file.\n\t\t Please, change your dat file."<< std::endl;
                      exit(1);}

// check, if there are comments lines
 do
{
  std::getline(infile, buf, '\n');
  if (buf[0] =='#') pos=infile.tellg();
}
 while (buf[0] == '#'); 

}

void AnsysFile::getPosition(const std::string seekexp, std::string::size_type & pos)
{
  infile.seekg(pos, std::ios::beg);
  std::string buf;
  pos=std::string::npos;

  std::cout << seekexp << std::endl;
  std::string::size_type hpos;

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

enum TypeBCs AnsysFile::TransformInDof(const std::string type_bc)
{
  enum TypeBCs result;

  if (type_bc=="vp-restr") result=vp_restraint;
  else
  if (type_bc=="ep-restr") result=ep_restraint;
  else Error(" This type of level for boundary condition in mesh-file (section [Node BC]) is unknown");

  return result;
}

//!
void AnsysFile::ReadElems(std::vector<Elem> & allelems)
{
#ifdef TRACE
 (*trace) << " entering AnsysFile::ReadElems " << std::endl;
#endif

 Integer maxnelems;
 ReadMaxnumelem(maxnelems); 

 allelems.resize(maxnelems);

 std::string::size_type pos=0;

 switch(dim_)
{
 case 2:
 getPosLine("2D Elements", pos);
 break;
 case 3:
 getPosLine("3D Elements", pos);
 break;
}
 infile.seekg(pos,std::ios::beg);

 if (ptTr || ptQ || ptTet)
  Error(" Pointers to BaseElem is not initialized",__FILE__,__LINE__);

 Integer i, ii, ibuf, itype, innodes;
 for (i=0; i<maxnelems; i++)
{
 infile >> ibuf >> itype >> innodes >> allelems[i].namesd;
 infile.ignore(100,'\n');

 allelems[i].ptElem=Type2ptElem(itype);   
 allelems[i].connect.Resize(innodes);
 for (ii=0; ii<innodes; ii++)
  infile >> allelems[i].connect[ii];

 infile.ignore(100,'\n'); 
}

}

//!
BaseElem * AnsysFile::Type2ptElem(const Integer itype)
{
  switch(itype)
{
case 4:
  return ptTr;
case 6:
  return ptQ;
case 8:
  return ptTet;
default:
  Error(" This type is unknown. ", __FILE__,__LINE__);
}
}

}
