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

 dim_=ReadDim();
}
  
AnsysFile :: ~AnsysFile()
{
#ifdef TRACE
  (*trace) << "entering AnsysFile::~AnsysFile" << std::endl;
#endif

 infile.close() ;
}

Integer AnsysFile::ReadDim()
{
#ifdef TRACE
  (*trace) << "entering Ansys::ReadDim" << std::endl;
#endif

 Integer dim;

 Integer i;
 std::string::size_type pos=0;
 getPosition("Dimension", pos);
 infile.seekg(pos,std::ios::beg);

 std::string auxname;
 infile >> dim;

 return dim;
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
 getPosLine("[Nodes]", pos);
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
}

void AnsysFile::ReadMaxnumelem(Integer & nelem,const std::string keyword)
{
 std::string::size_type pos=0;
 getPosition(keyword,pos);
 infile.seekg(pos,std::ios::beg);

 infile >> nelem;
}

void AnsysFile::ReadMaxnumnodesbc(Integer & nbc)
{
 std::string::size_type pos=0;
 getPosition("NumNodeBC", pos);
 infile.seekg(pos,std::ios::beg);

 infile >> nbc;
}

void AnsysFile::ReadBCs(std::list<Integer> * bcs, const std::vector<std::string> levels)
{
#ifdef TRACE
  (*trace) << "entering Ansys::ReadBCs" << std::endl;
#endif
 
 Integer numbc;
 ReadMaxnumnodesbc(numbc);

 std::string::size_type pos=0;
 getPosLine("Node BC", pos);
 infile.seekg(pos,std::ios::beg);
 
 std::string str;

 Integer nodalnum;
 Integer i,j;
 for (i=0; i < numbc; i++)
   {
    infile >> nodalnum >> str;
    infile.ignore(100,'\n');

    Boolean Find=FALSE;
    for (j=0; j<levels.size(); j++)
       { if (str==levels[j]) { Find=TRUE; break;}
       }         

    std::string msg=str+"-this level of BCs from .mesh file is not mentioned in .config file. Please, check .config-file";
    if (!Find) Error(msg.c_str(),__FILE__,__LINE__);

    bcs[j].push_back(nodalnum);
   } 
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

  std::string::size_type hpos;

  while (pos == std::string::npos && !infile.eof())
{
  hpos=infile.tellg();
  std::getline(infile, buf, '\n');
  pos=buf.find(seekexp);
}
  pos+=hpos+seekexp.length();

 if (pos>=pos_end) {std::cerr << "ERROR: (" << __FILE__ <<" "<< __LINE__
               << ") Cannot find string: " << seekexp ;
          std::cerr << " in your dat file.\n\t\t Please, change your dat file."<< std::endl;
                      exit(1);}
}


void AnsysFile::ReadEl(std::vector<Elem> * allelems, const std::vector<std::string> sd)
{
#ifdef TRACE
 (*trace) << " entering AnsysFile::ReadEl " << std::endl;
#endif

 switch(dim_)
{
 case 2:
 ReadEl2d(allelems,sd);
 break;
 case 3:
 ReadEl2d(allelems,sd);
 ReadEl3d(allelems,sd);
 break;
}

}

void AnsysFile::ReadEl2d(std::vector<Elem> * allelems, const std::vector<std::string> sd)
{
#ifdef TRACE
 (*trace) << " entering AnsysFile::ReadEl2D " << std::endl;
#endif

 Integer maxnelems;
 ReadMaxnumelem(maxnelems,"Num2DElements");

if (maxnelems)
{
 std::string::size_type pos=0;

 getPosLine("2D Elements", pos);
 infile.seekg(pos,std::ios::beg);

  if (!ptTr || !ptQ || !ptTet)
  Error(" Pointers to BaseElem is not initialized",__FILE__,__LINE__);

 Integer i, ii, j, ibuf, itype, innodes;
 std::string namesd;
 Elem el;
 for (i=0; i<maxnelems; i++)
{
 infile >> ibuf >> itype >> innodes >> namesd;
 infile.ignore(100,'\n');

 el.ptElem=Type2ptElem(itype);
 el.connect.Resize(innodes);
 for (ii=0; ii<innodes; ii++)
  infile >> el.connect[ii];

 infile.ignore(100,'\n');

 Boolean Find;
 for (j=0; j<sd.size(); j++)
  if (namesd == sd[j]) { allelems[j].push_back(el);
                         Find=TRUE;
                       }
 if (!Find) { std::string msg=namesd + "- this level of element is not mentioned in .conf-file. Please, check .config-file";
              Error(msg.c_str(),__FILE__,__LINE__);
            }
}
}
}

void AnsysFile::ReadEl3d(std::vector<Elem> * allelems, const std::vector<std::string> sd)
{
#ifdef TRACE
 (*trace) << " entering AnsysFile::ReadEl3d " << std::endl;
#endif

 Integer maxnelems;
 ReadMaxnumelem(maxnelems,"Num3DElements");

 std::string::size_type pos=0;
 getPosLine("3D Elements", pos);
 infile.seekg(pos,std::ios::beg);

  if (!ptTr || !ptQ || !ptTet)
  Error(" Pointers to BaseElem is not initialized",__FILE__,__LINE__);

 Integer i, ii, j, ibuf, itype, innodes;
 std::string namesd;
 Elem el;
 for (i=0; i<maxnelems; i++)
{
 infile >> ibuf >> itype >> innodes >> namesd;
 infile.ignore(100,'\n');

 el.ptElem=Type2ptElem(itype);
 el.connect.Resize(innodes);
 for (ii=0; ii<innodes; ii++)
  infile >> el.connect[ii];

 infile.ignore(100,'\n');

 Boolean Find;
 for (j=0; j<sd.size(); j++)
  if (namesd == sd[j]) { allelems[j].push_back(el);
                         Find=TRUE;
                       }
 if (!Find) { std::string msg=namesd + "- this level of element is not mentioned in .conf-file. Please, check .config-file";
              Error(msg.c_str(),__FILE__,__LINE__);
            }
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
