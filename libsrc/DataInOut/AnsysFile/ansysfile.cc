#include <fstream>
#include <iostream>
#include <string>
#include <stdarg.h>
#include <vector>

#include "ansysfile.hh"
#include "Domain/bcs.hh"
#include "Domain/GridCFS/grid_cfs.hh"

#ifdef ADAPTGRID
#include "Triangle.h"
#include "Tetrahedron.h"
#include "Hexahedron.h"
#endif

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
    (*trace) << "entering AnsysFile::ReadDim" << std::endl;
#endif

    Integer dim;

    std::string::size_type pos=0;
    getPosition("Dimension", pos);
    infile.seekg(pos,std::ios::beg);

    // std::string auxname;
    infile >> dim;

    return dim;
  }

  void AnsysFile::ReadCoordinate(Point<2> * const NodesCoord, const Integer maxnumnodes)
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
	infile >> ibuf >> NodesCoord[i][0] >> NodesCoord[i][1] ;
	infile.ignore(100,'\n');
      }
  }

  void AnsysFile::ReadCoordinate(Point<3> * const NodesCoord, const Integer maxnumnodes)
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
	infile >> ibuf >> NodesCoord[i][0] >> NodesCoord[i][1] >>  NodesCoord[i][2];
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

	std::string msg=str+": This level of BCs from .mesh file is not mentioned in .config file. Please, check .config-file";
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


  void AnsysFile::ReadEl(std::vector<Elem*> * allelems, const std::vector<std::string> sd)
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

#ifdef ADAPTGRID
  void AnsysFile::ReadBCs_GridRG(std::vector<Integer> & idBCs,std::vector<Integer> &colorBCs)
  {
#ifdef TRACE
    (*trace) << " entering AnsysFile::ReadGrid_RG " << std::endl;
#endif
  
    std::vector<std::string> levels;
    conf->getliststr("list_nodes",levels);

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

	idBCs.push_back(nodalnum);
	colorBCs.push_back(j);
      } 
  
  }

  void AnsysFile::ReadGrid_RG(std::vector<grd::Element*> & elems, std::vector<grd::Vertex*> * vertex, const std::vector<std::string> sd)
  {
#ifdef TRACE
    (*trace) << " entering AnsysFile::ReadGrid_RG " << std::endl;
#endif

    switch(dim_)
      {
      case 2:
	ReadEl4AdaptGrid2d(elems, vertex, sd);
	break;
      case 3:
	ReadEl4AdaptGrid3d(elems, vertex, sd);
	break;
      } 

#ifdef TRACE
    (*trace) << " leaving AnsysFile::ReadGrid_RG " << std::endl;
#endif
  }

  void AnsysFile::ReadEl4AdaptGrid2d(std::vector<grd::Element*> & elems, std::vector<grd::Vertex*> * vertex,  const std::vector<std::string> sd)
  {
#ifdef TRACE
    (*trace) << " entering AnsysFile::ReadElems4AdaptGrid " << std::endl;
#endif

    Integer maxnelems;
    ReadMaxnumelem(maxnelems,"Num2DElements");

    if (maxnelems)
      {
	std::string::size_type pos=0;

	getPosLine("2D Elements", pos);
	infile.seekg(pos,std::ios::beg);

	Integer i, ii, ibuf, itype, innodes;
	std::string namesd;
	Integer connect[4]; 

	elems.resize(maxnelems);
	for (i=0; i<maxnelems; i++)
	  {
	    infile >> ibuf >> itype >> innodes >> namesd;
	    infile.ignore(100,'\n');

	    for (ii=0; ii<innodes; ii++)
	      infile >> connect[ii];

	    grd::Triangle * tmpTri;
	    grd::Quadrangle * tmpQuad;
	    switch(itype)
	      {
	      case 4:
		tmpTri=new grd::Triangle;

		tmpTri->setVertex(0,(*vertex)[connect[0]-1]);
		tmpTri->setVertex(1,(*vertex)[connect[1]-1]);
		tmpTri->setVertex(2,(*vertex)[connect[2]-1]); 
   
		SetNumSD(tmpTri,namesd,sd);
    
		elems[i]=tmpTri; 
		break;
	      case 6:
		tmpQuad=new grd::Quadrangle;
   
		tmpQuad->setVertex(0,(*vertex)[connect[0]-1]);
		tmpQuad->setVertex(1,(*vertex)[connect[1]-1]);
		tmpQuad->setVertex(2,(*vertex)[connect[2]-1]);
		tmpQuad->setVertex(3,(*vertex)[connect[3]-1]);

		SetNumSD(tmpQuad,namesd,sd);

		elems[i]=tmpQuad;
		break;
 
	      default:
		Error(" This type of elems in mesh file is not implemented yet ");
	      }
	    infile.ignore(100,'\n'); 
	  } 
      } // end of if
  }

  void AnsysFile::ReadEl4AdaptGrid3d(std::vector<grd::Element*> & elems, std::vector<grd::Vertex*> * vertex,  const std::vector<std::string> sd)
  {
#ifdef TRACE
    (*trace) << " entering AnsysFile::ReadElems4AdaptGrid3d " << std::endl;
#endif

    Integer maxnelems;
    ReadMaxnumelem(maxnelems,"Num3DElements");

    if (maxnelems)
      {
	std::string::size_type pos=0;

	getPosLine("3D Elements", pos);
	infile.seekg(pos,std::ios::beg);

	Integer i, ii, ibuf, itype, innodes;
	std::string namesd;
	Integer connect[4]; 

	elems.resize(maxnelems);
	for (i=0; i<maxnelems; i++)
	  {
	    infile >> ibuf >> itype >> innodes >> namesd;
	    infile.ignore(100,'\n');

	    for (ii=0; ii<innodes; ii++)
	      infile >> connect[ii];

	    grd::Tetrahedron * tmpTetra;
	    grd::Hexahedron * tmpHexa;
	    switch(itype)
	      {
	      case 8:
		tmpTetra=new grd::Tetrahedron;

		tmpTetra->setVertex(0,(*vertex)[connect[0]-1]);
		tmpTetra->setVertex(1,(*vertex)[connect[1]-1]);
		tmpTetra->setVertex(2,(*vertex)[connect[2]-1]); 
		tmpTetra->setVertex(3,(*vertex)[connect[3]-1]);
   
		SetNumSD(tmpTetra,namesd,sd);
    
		elems[i]=tmpTetra; 
		break;

	      case 10:
		tmpHexa=new grd::Hexahedron;

		tmpTetra->setVertex(0,(*vertex)[connect[0]-1]);
		tmpTetra->setVertex(1,(*vertex)[connect[1]-1]);
		tmpTetra->setVertex(2,(*vertex)[connect[2]-1]); 
		tmpTetra->setVertex(3,(*vertex)[connect[3]-1]);
		tmpTetra->setVertex(4,(*vertex)[connect[4]-1]);
		tmpTetra->setVertex(5,(*vertex)[connect[5]-1]);
		tmpTetra->setVertex(6,(*vertex)[connect[6]-1]); 
		tmpTetra->setVertex(7,(*vertex)[connect[7]-1]);

		SetNumSD(tmpHexa,namesd,sd);

		elems[i]=tmpHexa;

		break;
 
	      default:
		Error(" This type of elems in mesh file is not implemented yet ");
	      }
	    infile.ignore(100,'\n'); 
	  }

      } // end of if
  }

  void AnsysFile::SetNumSD(grd::Element * ptEl, const std::string namesd, const std::vector<std::string> sd)
  {
    Boolean Find;
    Integer j;
    for (j=0; j<sd.size(); j++)
      if (namesd == sd[j]) { ptEl->setValue(j);
      Find=TRUE;
      }
    if (!Find) { std::string msg=namesd + "- this level of element is not mentioned in .conf-file. Please, check .config-file";
    Error(msg.c_str(),__FILE__,__LINE__);
    }
  }

#endif

  void AnsysFile::ReadEl1d(std::vector<Elem*> * allelems, const std::vector<std::string> sd)
  {
#ifdef TRACE
    (*trace) << " entering AnsysFile::ReadEl1D " << std::endl;
#endif
    
    Integer maxnelems;
    ReadMaxnumelem(maxnelems,"Num1DElements");
    
    if (maxnelems)
      {
	std::string::size_type pos=0;

	getPosLine("1D Elements", pos);
	infile.seekg(pos,std::ios::beg);

	if (!ptL1)
	  Error(" Pointers to BaseElem is not initialized",__FILE__,__LINE__);

	Integer i, ii, j, inum, itype, innodes;
	std::string namesd;

	for (i=0; i<maxnelems; i++)
	  {
	    Elem * el=new Elem();
	    infile >> inum >> itype >> innodes >> namesd;
	    infile.ignore(100,'\n');

	    el->ElemNum=inum;
	    el->ptElem=Type2ptElem(itype);
	    el->namesd = namesd;
	   	
	    el->connect.Resize(innodes);
	    for (ii=0; ii<innodes; ii++)
	      infile >> el->connect[ii];

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

  void AnsysFile::ReadEl2d(std::vector<Elem*> * allelems, const std::vector<std::string> sd)
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

	//	if (!ptTr || !ptQ || !ptTet)
	if (!ptQ || !ptTr1)
	  Error(" Pointers to BaseElem is not initialized",__FILE__,__LINE__);

	Integer i, ii, j, inum, itype, innodes;
	std::string namesd;

	for (i=0; i<maxnelems; i++)
	  {
	    Elem * el=new Elem();
	    infile >> inum >> itype >> innodes >> namesd;
	    infile.ignore(100,'\n');
	    
	    el->ElemNum=inum;
	    el->ptElem=Type2ptElem(itype);
	    el->namesd=namesd;
	    el->connect.Resize(innodes);

	    for (ii=0; ii<innodes; ii++)
	      infile >> el->connect[ii];

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


  void AnsysFile::ReadEl3d(std::vector<Elem*> * allelems, const std::vector<std::string> sd)
  {
#ifdef TRACE
    (*trace) << " entering AnsysFile::ReadEl3d " << std::endl;
#endif

    Integer maxnelems;
    ReadMaxnumelem(maxnelems,"Num3DElements");

    std::string::size_type pos=0;
    getPosLine("3D Elements", pos);
    infile.seekg(pos,std::ios::beg);

    //    if (!ptTr || !ptQ || !ptTet)
    if (!ptTet)
      Error(" Pointers to BaseElem is not initialized",__FILE__,__LINE__);

    Integer i, ii, j, inum, itype, innodes;
    std::string namesd;

    for (i=0; i<maxnelems; i++)
      {
	Elem * el=new Elem();
	infile >> inum >> itype >> innodes >> namesd;
	infile.ignore(100,'\n');

	el->ElemNum=inum;
	el->ptElem=Type2ptElem(itype);
	el->namesd = namesd;
	el->connect.Resize(innodes);
	for (ii=0; ii<innodes; ii++)
	  infile >> el->connect[ii];

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
  BaseFE * AnsysFile::Type2ptElem(const Integer itype)
  {
    switch(itype)
      {
      case 100:
	return ptL1;
      case 4:
	return ptTr1;
      case 6:
	return ptQ;
      case 8:
	return ptTet;
	//  case 10:
	//    return ptHexa;
      default:
	Error(" This type is unknown. ", __FILE__,__LINE__);
      }
  }
}
