#include <fstream>
#include <iostream>
#include <string>
#include <stdarg.h>

#include "ansysfile.hh"
#include "Domain/bcs.hh"
#include "Domain/GridCFS/grid_cfs.hh"
#include "DataInOut/WriteInfo.hh"
#include "Utils/StdVector.hh"

#ifdef ADAPTGRID
#include "DataInOut/ParamHandling/ConfFile.hh"
#include "Triangle.h"
#include "Tetrahedron.h"
#include "Hexahedron.h"
#endif

namespace CoupledField
{

AnsysFile::AnsysFile(const Char * const afilename)
:FileType(afilename)
{
  ENTER_FCN( "AnsysFile::AnsysFile" );

  maxNumElems_ = 0;
  maxNumNodes_ = 0;

    infile.open(strcat(filename,".mesh"));
    if (!infile.good()) {std::cerr << "ERROR(" << __FILE__ << " " << __LINE__ <<
			   ") Can't open " << filename << std::endl;
    exit(1);}

    infile.seekg(0, std::ios::end);
    pos_end=infile.tellg();

    dim_=ReadDim();
}
  
AnsysFile::~AnsysFile()
{
  ENTER_FCN( "AnsysFile::~AnsysFile" );
  
  infile.close() ;
}

Integer AnsysFile::ReadDim()
{
  ENTER_FCN( "AnsysFile::ReadDim" );
  
    Integer dim;
    
    std::string::size_type pos=0;
    getPosition("Dimension", pos);
    infile.seekg(pos,std::ios::beg);
    
    // std::string auxname;
    infile >> dim;
    
    return dim;
}

Integer AnsysFile::GetNum3DElems()
{
  ENTER_FCN( "AnsysFile::GetNum3DElems" );
  
    Integer num;
    
    std::string::size_type pos=0;
    getPosition("Num3DElements", pos);
    infile.seekg(pos,std::ios::beg);
    infile >> num;
    return num;
}

Integer AnsysFile::GetNum2DElems()
{
  ENTER_FCN( "AnsysFile::GetNum2DElems" );
  
    Integer num;
    
    std::string::size_type pos=0;
    getPosition("Num2DElements", pos);
    infile.seekg(pos,std::ios::beg);
    infile >> num;
    return num;
}

Integer AnsysFile::GetNum1DElems()
{
  ENTER_FCN( "AnsysFile::GetNum1DElems" );
  
    Integer num;
    
    std::string::size_type pos=0;
    getPosition("Num1DElements", pos);
    infile.seekg(pos,std::ios::beg);
    infile >> num;
    return num;
}

Integer AnsysFile::GetNumBCs()
{
  ENTER_FCN( "AnsysFile::GetNumBCs" );
  
    Integer num;
    
    std::string::size_type pos=0;
    getPosition("NumNodeBC", pos);
    infile.seekg(pos,std::ios::beg);
    infile >> num;
    return num;
}



Integer AnsysFile::GetNumSaveNodes()
{
  ENTER_FCN( "AnsysFile::GetNumSaveNodes" );

    Integer nrSaveNodes;
    ReadNumSaveNodes(nrSaveNodes);
    return nrSaveNodes;
}



void AnsysFile::ReadCoordinate(Point<2> * const NodesCoord, const Integer maxnumnodes)
{
  ENTER_FCN( "Ansys::ReadCoordinate 2D" );

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
  ENTER_FCN( "Ansys::ReadCoordinate 3D" );
  
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
  ENTER_FCN( "Ansys::ReadMaxnumnodes" );

    std::string::size_type pos=0;
    getPosition("NumNodes", pos);
    infile.seekg(pos,std::ios::beg);

    infile >> nnodes;
    maxNumNodes_ = nnodes;
}

void AnsysFile::ReadMaxnumelem(Integer & nelem,const std::string keyword)
{
  ENTER_FCN( "AnsysFile::ReadMaxnumelem" );

  std::string::size_type pos=0;
  getPosition(keyword,pos);
  infile.seekg(pos,std::ios::beg);
  
  infile >> nelem;
  maxNumElems_ += nelem;
}

void AnsysFile::ReadMaxnumnodesbc(Integer & nbc)
{
  ENTER_FCN( "AnsysFile::ReadMaxnumnodesbc" );
  std::string::size_type pos=0;
  getPosition("NumNodeBC", pos);
  infile.seekg(pos,std::ios::beg);
  
  infile >> nbc;
}


void AnsysFile::ReadNumSaveNodes(Integer & nrSaveNodes)
{
  ENTER_FCN( "AnsysFile::ReadNumSaveNodes" );
  std::string::size_type pos=0;
  getPosition("NumSaveNodes", pos);
  infile.seekg(pos,std::ios::beg);
  
  infile >> nrSaveNodes;
}


void AnsysFile::ReadBCs(std::list<Integer> * bcs, const StdVector<std::string> levels)
{
  ENTER_FCN( "AnsysFile::ReadBCs" );
    
    Integer numbc;
    ReadMaxnumnodesbc(numbc);

    // Vector which contains Booleans for all boundary nodes for all levels.
    // It is used to determine, wether a node occurs several times, which could
    // lead to errors with the number of boundary conditions later on
    StdVector<StdVector<Integer> > nodesPerLevel;
    nodesPerLevel.Resize(levels.GetSize());
    

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

	// Find the according level
	for (j=0; j<levels.GetSize(); j++)
	  if (str==levels[j]) 
	    {
	      Find=TRUE;
	      break;
	    }
   
	if (!Find) 
	  {
	    std::string msg=str+": This level of BCs from .mesh file is not mentioned in xml file. Please, check .xml-file";
	    Error(msg.c_str(),__FILE__,__LINE__);
	  }
	
	// Check, if node was already read in
	Boolean onlyNode = TRUE;
	for (Integer iNode=0; iNode<nodesPerLevel[j].GetSize(); iNode++) 
	  if (nodesPerLevel[j][iNode] == nodalnum) {
	    onlyNode = FALSE;
	    break;
	  }

	if (onlyNode == FALSE)
	  {
	    std::string warnMsg = "ReadBCs: The node with Nr. ";
	    warnMsg += Info->GenStr(nodalnum);
	    warnMsg += " of BC level '";
	    warnMsg += str;
	    warnMsg += "' occured at least two times in the .mesh file.\n";
	    warnMsg += "Please make sure that each node occurs only on time in ";
	    warnMsg += ".mesh file, otherwise some undefined errors may occur!";
	    Warning(warnMsg.c_str(), __FILE__, __LINE__);
	  }
	else
	  {
	    nodesPerLevel[j].Push_back(nodalnum);
	    bcs[j].push_back(nodalnum);
	  }
      } 
}





void AnsysFile::ReadSaveNodes(StdVector<Integer> & saveNodes , const std::string level)
{
  ENTER_FCN( "Ansys::ReadSaveNodes" );

    Integer nrSaveNodes;
    ReadNumSaveNodes(nrSaveNodes);
    saveNodes.Clear();

    std::string::size_type pos=0;
    getPosLine("Save Nodes", pos);
    infile.seekg(pos,std::ios::beg);
    
    std::string str;

    Integer nodalnum;
    Integer i,j;
    for (i=0; i < nrSaveNodes; i++)
      {
	infile >> nodalnum >> str;
	infile.ignore(100,'\n');
	
	if (str == level) 
	  	  saveNodes.Push_back(nodalnum);
      } 
 
    if (saveNodes.GetSize() == 0) {
      std::string msg="The level \"" + str + "\" of 'saveNodes'";
      msg += "is not mentioned in the .mesh file.\n";
      msg += "History nodes are written with the command 'wsavnod'.";
      Error(msg.c_str());
    }
}




void AnsysFile::ReadLevelOfSaveNodes(StdVector<std::string>& levels)
{
  ENTER_FCN( "Ansys::ReadLevelOfSaveNodes" );
    
    Integer nrSaveNodes;
    ReadNumSaveNodes(nrSaveNodes);

    std::string::size_type pos=0;
    getPosLine("Save Nodes", pos);
    infile.seekg(pos,std::ios::beg);
    
    std::string str;

    Integer nodalnum;
    Integer i,j;
    for (i=0; i < nrSaveNodes; i++)
      {
	infile >> nodalnum >> str;
	infile.ignore(100,'\n');
	
	Boolean found=FALSE;

	for (j=0; j<levels.GetSize(); j++) 
	  if (str==levels[j])
	    found=TRUE;

	if (!found) 
	  levels.Push_back(str);
      } 
}




void AnsysFile::ReadBCsConf(StdVector<std::string> &levels)
{
  ENTER_FCN( "Ansys::ReadBCsConf" );
    
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
	
	if (i==0) 
	    levels.Push_back(str);
	else
	  {
	    Integer find = 0;
	    for (j=0; j<levels.GetSize(); j++)
		if (str == levels[j]) find = 1;

	    if (!find) levels.Push_back(str);	      
	  }
      } 
}

void AnsysFile::getPosLine(const std::string seekexp, std::string::size_type & pos)
{
  ENTER_IFCN( "AnsysFile::getPosLine" );
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
  ENTER_IFCN( "AnsysFile::getPosition" );
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


void AnsysFile::ReadEl(StdVector<Elem*> * allelems, const StdVector<std::string> sd)
{
  ENTER_FCN( "AnsysFile::ReadEl" );
  
    switch(dim_)
      {
      case 2:
	ReadEl2d(allelems,sd);
	break;
      case 3:
	ReadEl3d(allelems,sd);
	break;
      }
  }

#ifdef ADAPTGRID
void AnsysFile::ReadBCs_GridRG(std::list<Integer> & idBCs,
			       StdVector<Integer> &colorBCs)
{
  ENTER_FCN( "AnsysFile::ReadGrid_RG" );
  
  StdVector<std::string> levels;
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
	for (j=0; j<levelsGetSize(); j++)
	  { if (str==levels[j]) { Find=TRUE; break;}
	  }         
	
	std::string msg=str+"-this level of BCs from .mesh file is not mentioned in .xml file. Please, check .xml-file";
	if (!Find) Error(msg.c_str(),__FILE__,__LINE__);
	
	idBCs.Push_back(nodalnum);
	colorBCs.Push_back(j);
      } 
  
}

void AnsysFile::ReadGrid_RG(StdVector<grd::Element*> & elems, StdVector<grd::Vertex*> * vertex, const StdVector<std::string> sd)
{
  ENTER_FCN( "AnsysFile::ReadGrid_RG ");

    switch(dim_)
      {
      case 2:
	ReadEl4AdaptGrid2d(elems, vertex, sd);
	break;
      case 3:
	ReadEl4AdaptGrid3d(elems, vertex, sd);
	break;
      } 
  }

void AnsysFile::ReadEl4AdaptGrid2d(StdVector<grd::Element*> & elems, 
				   StdVector<grd::Vertex*> * vertex,  
				   const StdVector<std::string> sd)
  {
    ENER_FCN( "AnsysFile::ReadElems4AdaptGrid" );

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

void AnsysFile::ReadEl4AdaptGrid3d(StdVector<grd::Element*> & elems, 
				   StdVector<grd::Vertex*> * vertex,  
				   const StdVector<std::string> sd)
{
  ENTER_FCN( "AnsysFile::ReadElems4AdaptGrid3d" );

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

void AnsysFile::SetNumSD(grd::Element * ptEl, const std::string namesd, const StdVector<std::string> sd)
{
  ENTER_FCN( "AnsysFile::SetNumSD" );
  Boolean Find;
  Integer j;
  for (j=0; j<sdGetSize(); j++)
    if (namesd == sd[j]) { ptEl->setValue(j);
    Find=TRUE;
    }
  if (!Find) { std::string msg=namesd + "- this level of element is not mentioned in .xml-file. Please, check .xml-file";
  Error(msg.c_str(),__FILE__,__LINE__);
  }
  }

#endif

void AnsysFile::ReadEl1d(StdVector<Elem*> * allelems, const StdVector<std::string> sd) 
{
  ENTER_FCN( "AnsysFile::ReadEl1D" );
    
    Integer maxnelems;
    ReadMaxnumelem(maxnelems,"Num1DElements");
    
    if (maxnelems)
      {
	std::string::size_type pos=0;

	getPosLine("1D Elements", pos);
	infile.seekg(pos,std::ios::beg);

	if (!ptL1 || !ptL2)
	  Error(" Pointers to BaseElem is not initialized",__FILE__,__LINE__);

	Integer i, ii, j, inum, itype, innodes;
	std::string namesd;

	for (i=0; i<maxnelems; i++)
	  {
	    Elem * el=new Elem();
	    infile >> inum >> itype >> innodes >> namesd;
	    infile.ignore(100,'\n');

	    el->elemNum=inum;
	    if (inum > maxNumElems_) 
	      {
		std::string errMsg = "The current element number is higher than the ";
		errMsg += "maximum number of elements in your .mesh-file. Something might ";
		errMsg += "have gone wrong in the meshing process.";
		Error(errMsg.c_str(), __FILE__, __LINE__);
	      }
	    el->ptElem=Type2ptElem(itype);
	    el->namesd = namesd;
	   	
	    el->connect.Resize(innodes);
	    for (ii=0; ii<innodes; ii++)
	      infile >> el->connect[ii];

	    infile.ignore(100,'\n');

	    Boolean Find=FALSE;
	    for (j=0; j<sd.GetSize(); j++)
	      if (namesd == sd[j]) { allelems[j].Push_back(el);
	      Find=TRUE;
	      }
	    if (!Find) { std::string msg=namesd + "- this level of element is not mentioned in .xml-file. Please, check .xml-file";
	    Error(msg.c_str(),__FILE__,__LINE__);
	    }
	  }
      }

 
  }

void AnsysFile::ReadEl2d(StdVector<Elem*> * allelems, const StdVector<std::string> sd)
{
  ENTER_FCN( "AnsysFile::ReadEl2D" );

    Integer actElemNr = 0;
    Integer maxnelems;

    ReadMaxnumelem(maxnelems,"Num2DElements");
    
    if (maxnelems)
      {
	std::string::size_type pos=0;

	getPosLine("2D Elements", pos);
	infile.seekg(pos,std::ios::beg);

	if (!ptQ || !ptQ2 || !ptTr1 || !ptTr2)
	  Error(" Pointers to BaseElem is not initialized",__FILE__,__LINE__);

	Integer i, ii, j, inum, itype, innodes;
	std::string namesd;

	for (i=0; i<maxnelems; i++)
	  {
	    Elem * el=new Elem();
	    infile >> inum >> itype >> innodes >> namesd;
	    infile.ignore(100,'\n');
	    
	    el->elemNum=inum;
	    if (inum > maxNumElems_) 
	      {
		std::string errMsg = "The current element number is higher than the ";
		errMsg += "maximum number of elements in your .mesh-file.\n  Something might ";
		errMsg += "have gone wrong in the meshing process.";
		Error(errMsg.c_str(), __FILE__, __LINE__);
	      }
	    el->ptElem=Type2ptElem(itype);
	    el->namesd=namesd;
	    el->connect.Resize(innodes);

	    for (ii=0; ii<innodes; ii++)
	      infile >> el->connect[ii];

	    infile.ignore(100,'\n');

	    Boolean Find=FALSE;
	    for (j=0; j<sd.GetSize(); j++)
	      if (namesd == sd[j]) 
		{ 
		  allelems[j].Push_back(el);
		  actElemNr++;
		  Find=TRUE;
		}
	    if (!Find) 
	      { 
		std::string msg=namesd + "- this level of element is not mentioned in .xml-file. Please, check .xml-file";
		Error(msg.c_str(),__FILE__,__LINE__);
	      }
	  }
      }
  }


void AnsysFile::ReadEl3d(StdVector<Elem*> * allelems, const StdVector<std::string> sd)
  {
    ENTER_FCN( "AnsysFile::ReadEl3d" );

    Integer maxnelems;
    ReadMaxnumelem(maxnelems,"Num3DElements");

    std::string::size_type pos=0;
    getPosLine("3D Elements", pos);
    infile.seekg(pos,std::ios::beg);

    //    if (!ptTr || !ptQ || !ptTet)
    if (!ptTet || !ptHexa || !ptPyra || !ptWedge)
      Error(" Pointers to BaseElem is not initialized",__FILE__,__LINE__);

    Integer i, ii, j, inum, itype, innodes;
    std::string namesd;

    for (i=0; i<maxnelems; i++)
      {
	Elem * el=new Elem();
	infile >> inum >> itype >> innodes >> namesd;
	infile.ignore(100,'\n');

	el->elemNum=inum;
	if (inum > maxNumElems_) 
	  {
	    std::string errMsg = "The current element number is higher than the ";
	    errMsg += "maximum number of elements in your .mesh-file. Something might ";
	    errMsg += "have gone wrong in the meshing process.";
	    Error(errMsg.c_str(), __FILE__, __LINE__);
	  }
	el->ptElem=Type2ptElem(itype);
	el->namesd = namesd;
	el->connect.Resize(innodes);
	for (ii=0; ii<innodes; ii++)
	  infile >> el->connect[ii];

	infile.ignore(100,'\n');

	Boolean Find=FALSE;
	for (j=0; j<sd.GetSize(); j++)
	  if (namesd == sd[j]) { allelems[j].Push_back(el);
	  Find=TRUE;
	  }
	if (!Find) { std::string msg=namesd + "- this level of element is not mentioned in .xml-file. Please, check .xml-file";
	Error(msg.c_str(),__FILE__,__LINE__);
	}
      }

  }

void AnsysFile::ReadEl3dConf(StdVector<std::string> &sd)
{
  ENTER_FCN( "AnsysFile::ReadEl3dConf" );

    Integer maxnelems;
    ReadMaxnumelem(maxnelems,"Num3DElements");

    std::string::size_type pos=0;
    getPosLine("3D Elements", pos);
    infile.seekg(pos,std::ios::beg);

    if (!ptTet || !ptHexa  || !ptPyra || !ptWedge)
      Error(" Pointers to BaseElem is not initialized",__FILE__,__LINE__);

    Integer i, ii, j, inum, itype, innodes;
    std::string namesd;

    for (i=0; i<maxnelems; i++)
      {
	Elem * el=new Elem();
	infile >> inum >> itype >> innodes >> namesd;
	infile.ignore(100,'\n');

	el->elemNum=inum;
	if (inum > maxNumElems_) 
	  {
	    std::string errMsg = "The current element number is higher than the ";
	    errMsg += "maximum number of elements in your .mesh-file. Something might ";
	    errMsg += "have gone wrong in the meshing process.";
	    Error(errMsg.c_str(), __FILE__, __LINE__);
	  }
	el->ptElem=Type2ptElem(itype);
	el->namesd = namesd;
	el->connect.Resize(innodes);
	for (ii=0; ii<innodes; ii++)
	  infile >> el->connect[ii];

	infile.ignore(100,'\n');
	if (i==0) 
	    sd.Push_back(namesd);
	else
	  {
	    Integer find = 0;
	    for (j=0; j<sd.GetSize(); j++)
		if (namesd == sd[j]) find = 1;

	    if (!find) sd.Push_back(namesd);
	      
	  }
      }

  }


void AnsysFile::ReadEl2dConf(StdVector<std::string> &sd)
{
  ENTER_FCN( "AnsysFile::ReadEl2dConf" );

    Integer maxnelems;
    ReadMaxnumelem(maxnelems,"Num2DElements");

    std::string::size_type pos=0;
    getPosLine("2D Elements", pos);
    infile.seekg(pos,std::ios::beg);

    if (!ptQ || !ptTr1 || !ptTr2)
      Error(" Pointers to BaseElem is not initialized",__FILE__,__LINE__);

    Integer i, ii, j, inum, itype, innodes;
    std::string namesd;

    for (i=0; i<maxnelems; i++)
      {
	Elem * el=new Elem();
	infile >> inum >> itype >> innodes >> namesd;
	infile.ignore(100,'\n');

	el->elemNum=inum;
	if (inum > maxNumElems_) 
	  {
	    std::string errMsg = "The current element number is higher than the ";
	    errMsg += "maximum number of elements in your .mesh-file. Something might ";
	    errMsg += "have gone wrong in the meshing process.";
	    Error(errMsg.c_str(), __FILE__, __LINE__);
	  }
	el->ptElem=Type2ptElem(itype);
	el->namesd = namesd;
	el->connect.Resize(innodes);
	for (ii=0; ii<innodes; ii++)
	  infile >> el->connect[ii];

	infile.ignore(100,'\n');
	//	allelems.Push_back(el);
	if (i==0) 
	    sd.Push_back(namesd);
	else
	  {
	    Integer find = 0;
	    for (j=0; j<sd.GetSize(); j++)
		if (namesd == sd[j]) find = 1;

	    if (!find) sd.Push_back(namesd);
	      
	  }
      }

}

void AnsysFile::ReadEl1dConf(StdVector<std::string> &sd)
{
  ENTER_FCN( "AnsysFile::ReadEl1dConf" );

    Integer maxnelems;
    ReadMaxnumelem(maxnelems,"Num1DElements");

    std::string::size_type pos=0;
    getPosLine("1D Elements", pos);
    infile.seekg(pos,std::ios::beg);

    if (!ptL1 || !ptL2)
      Error(" Pointers to BaseElem is not initialized",__FILE__,__LINE__);

    Integer i, ii, j, inum, itype, innodes;
    std::string namesd;

    for (i=0; i<maxnelems; i++)
      {
	Elem * el=new Elem();
	infile >> inum >> itype >> innodes >> namesd;
	infile.ignore(100,'\n');

	el->elemNum=inum;
	if (inum > maxNumElems_) 
	  {
	    std::string errMsg = "The current element number is higher than the ";
	    errMsg += "maximum number of elements in your .mesh-file. Something might ";
	    errMsg += "have gone wrong in the meshing process.";
	    Error(errMsg.c_str(), __FILE__, __LINE__);
	  }
	el->ptElem=Type2ptElem(itype);
	el->namesd = namesd;
	el->connect.Resize(innodes);
	for (ii=0; ii<innodes; ii++)
	  infile >> el->connect[ii];

	infile.ignore(100,'\n');
	//	allelems.Push_back(el);
	if (i==0) 
	    sd.Push_back(namesd);
	else
	  {
	    Integer find = 0;
	    for (j=0; j<sd.GetSize(); j++)
		if (namesd == sd[j]) find = 1;

	    if (!find) sd.Push_back(namesd);
	      
	  }
      }

}

//!
BaseFE * AnsysFile::Type2ptElem(const Integer itype)
{
  ENTER_IFCN( "AnsysFile::Type2ptElem" );
    switch(itype)
      {
      case 1:
	return ptL2;
      case 100:
	return ptL1;
      case 4:
	return ptTr1;
      case 5:
	return ptTr2;
      case 6:
	return ptQ;
      case 7:
	return ptQ2;
      case 8:
	return ptTet;
      case 10:
	return ptHexa;
      case 12:
	return ptPyra;
      case 13:
	return ptWedge;
      default:
	{
	  std::cout << "Used Element Type: " << itype << std::endl;  
	  Error(" This element type is unknown. ", __FILE__,__LINE__);
	}
	
      }
}
}
