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


namespace CoupledField {


  // ***************
  //   Constructor
  // ***************
  AnsysFile::AnsysFile( const Char *const afilename ) : FileType(afilename) {

    ENTER_FCN( "AnsysFile::AnsysFile" );

    maxNumNodes_   = 0;
    actMaxElemNum_ = 0;

    infile.open(strcat(filename,".mesh"));
    
    if (!infile.good()) {
      std::cerr << "ERROR(" << __FILE__ << " " << __LINE__ <<
        ") Can't open " << filename << std::endl;
      exit(1);
    }
    
    infile.seekg(0, std::ios::end);
    pos_end = infile.tellg();
    
    dim_ = ReadDim();

    // Determine overall number of elements. This is required
    // for reserving space for a vector that stores pointers
    // to the elements in ascending order with respect to the
    // element numbers.
    maxNumElems_  = 0;
    maxNumElems_ += GetNum1DElems();
    maxNumElems_ += GetNum2DElems();
    maxNumElems_ += GetNum3DElems();
  }


  // **************
  //   Destructor
  // **************
  AnsysFile::~AnsysFile() {
    ENTER_FCN( "AnsysFile::~AnsysFile" );
    infile.close() ;
  }


  // *******************************
  //   ReadCoordinate (2D version)
  // *******************************
  void AnsysFile::ReadCoordinate( Point<2> *const NodesCoord,
                                  const Integer maxnumnodes ) {

    ENTER_FCN( "Ansys::ReadCoordinate 2D" );

    Integer i;
    std::string::size_type pos=0;
    getPosLine("[Nodes]", pos);
    infile.seekg(pos,std::ios::beg);
  
    Integer ibuf;
  
    for ( i = 0; i < maxnumnodes; i++ ) {
      infile >> ibuf >> NodesCoord[i][0] >> NodesCoord[i][1] ;
      infile.ignore(100,'\n');
    }
  }


  // *******************************
  //   ReadCoordinate (3D version)
  // *******************************
  void AnsysFile::ReadCoordinate( Point<3> *const NodesCoord,
                                  const Integer maxnumnodes ) {

    ENTER_FCN( "Ansys::ReadCoordinate 3D" );
  
    Integer i;
    std::string::size_type pos=0;
    getPosLine("[Nodes]", pos);
    infile.seekg(pos,std::ios::beg);
  
    Integer ibuf;
    for ( i = 0; i < maxnumnodes; i++ ) {
      infile >> ibuf >> NodesCoord[i][0] >> NodesCoord[i][1]
             >> NodesCoord[i][2];
      infile.ignore(100,'\n');
    }
  }


  // ***********
  //   ReadBCs
  // ***********
  void AnsysFile::ReadBCs( std::list<Integer> *bcs,
                           const StdVector<std::string> levels ) {

    ENTER_FCN( "AnsysFile::ReadBCs" );
    
    Integer numbc;
    ReadMaxnumnodesbc(numbc);

    // Vector which contains Booleans for all boundary nodes for all levels.
    // It is used to determine, wether a node occurs several times, which
    // could lead to errors with the number of boundary conditions later on
    StdVector<StdVector<Integer> > nodesPerLevel;
    nodesPerLevel.Resize(levels.GetSize());
    

    std::string::size_type pos=0;
    std::string::size_type lineEndPos =0;
    getPosLine("Node BC", pos);
    infile.seekg(pos,std::ios::beg);
    
    std::string str, buf, errMsg;

    Integer nodalnum;
    Integer i,j;
    for ( i = 0; i < numbc; i++ ) {
        
      // remember current position and get the position of endline
      pos = infile.tellg();
      std::getline(infile,buf,'\n');
      lineEndPos=infile.tellg();
      infile.seekg(pos,std::ios::beg);
        
      // try to read in the data
      infile >> nodalnum >> str;

      // if read in was successfull, enline position and current
      // position are the same
      infile.ignore(100,'\n');
      pos = infile.tellg();
      if (pos != lineEndPos) {
        (*error) << "AnsysFile:ReadBCs: The node list for the boundary "
                 << "conditions has wrong size or format. Please correct it!";
        Error( __FILE__, __LINE__ );
      }
        
      Boolean Find=FALSE;

      // Find the associated level
      for (j=0; j<levels.GetSize(); j++) {
        if (str==levels[j]) {
          Find=TRUE;
          break;
        }
      }
      if (!Find) {
        (*error) << " The region '" << str << "' for BCs from .mesh file "
                 << "is not mentioned in xml file. Please, check .xml-file";
        Error( __FILE__,__LINE__ );
      }
        
      // Check, if node was already read in
      Boolean onlyNode = TRUE;
      for (Integer iNode=0; iNode<nodesPerLevel[j].GetSize(); iNode++) 
        if (nodesPerLevel[j][iNode] == nodalnum) {
          onlyNode = FALSE;
          break;
        }

      if (onlyNode == FALSE) {
        std::string warnMsg = "ReadBCs: The node with Nr. ";
        warnMsg += Info->GenStr(nodalnum);
        warnMsg += " of BC level '";
        warnMsg += str;
        warnMsg += "' occured at least two times in the .mesh file.\n";
        warnMsg += "Please make sure that each node occurs only on time in ";
        warnMsg += ".mesh file, otherwise some undefined errors may occur!";
        Warning(warnMsg.c_str(), __FILE__, __LINE__);
      }
      else {
        nodesPerLevel[j].Push_back(nodalnum);
        bcs[j].push_back(nodalnum);
      }
    } 

    if (! IsNextLineEmpty(pos)) {
      (*warning) << "AnsysFile::ReadBCs: The line after the last BC node "
                 << "no. " << nodalnum << " in region '" << str
                 << "' seems to contain nodes too. Please check if the "
                 << "number of BC nodes specified in the header of the "
                 << "mesh-file matches the real number of BC  nodes!";
      Warning( __FILE__, __LINE__ );
    }
  }


  // *****************
  //   ReadSaveNodes
  // *****************
  void AnsysFile::ReadSaveNodes( StdVector<Integer> &saveNodes,
                                 const std::string level ) {

    ENTER_FCN( "Ansys::ReadSaveNodes" );

    Integer nrSaveNodes;
    ReadNumSaveNodes(nrSaveNodes);
    saveNodes.Clear();

    std::string::size_type pos=0;
    std::string::size_type lineEndPos=0;
    
    getPosLine("Save Nodes", pos);
    infile.seekg(pos,std::ios::beg);
    
    std::string str, buf;

    Integer nodalnum;
    Integer i;
    for ( i = 0; i < nrSaveNodes; i++ ) {

      // remember current position and get the position of endline
      pos = infile.tellg();
      std::getline(infile,buf,'\n');
      lineEndPos=infile.tellg();
      infile.seekg(pos,std::ios::beg);
        
      // Read in data
      infile >> nodalnum >> str;

      // if read in was successfull, enline position and current
      // position are the same
      infile.ignore(100,'\n');
      pos = infile.tellg();
      if (pos != lineEndPos){
        (*error) << "AnsysFile:ReadSaveNodes: The node list for the "
                 << "saveNodes has wrong size or format. Please correct it!";
        Error( __FILE__, __LINE__ );
      }

      if (str == level) {
        saveNodes.Push_back(nodalnum);
      }
    } 
 
    if (saveNodes.GetSize() == 0) {
      (*error) << " The region '" << str << "' for saveNodes from .mesh file "
               << "is not mentioned in xml file. Please, check .xml-file. "
               << "History nodes are written with the command 'wsavnod'.";
      Error( __FILE__, __LINE__ );
    }

    if (! IsNextLineEmpty(pos)) {
      (*warning) << "AnsysFile::ReadSaveNodes: The line after the last "
                 << "save node no. " << nodalnum << " in region '" << str
                 << "' seems to contain nodes too. Please check if the "
                 << "number of BC nodes specified in the header of the "
                 << "mesh-file matches the real number of BC  nodes!";
      Warning( __FILE__, __LINE__ );
    }
  }


  // ************************
  //   ReadSaveLevelOfNodes
  // ************************
  void AnsysFile::ReadLevelOfSaveNodes( StdVector<std::string> &levels ) {

    ENTER_FCN( "Ansys::ReadLevelOfSaveNodes" );
    
    Integer nrSaveNodes;
    ReadNumSaveNodes(nrSaveNodes);

    std::string::size_type pos=0;
    getPosLine("Save Nodes", pos);
    infile.seekg(pos,std::ios::beg);
    
    std::string str;

    Integer nodalnum;
    Integer i,j;
    for ( i = 0; i < nrSaveNodes; i++ ) {
      infile >> nodalnum >> str;
      infile.ignore(100,'\n');

      Boolean found = FALSE;

      for ( j = 0; j < levels.GetSize(); j++ ) {
        if ( str == levels[j] ) {
          found = TRUE;
        }
      }
      if ( !found ) {
        levels.Push_back(str);
      }
    }
  }


  // ***************
  //   ReadBCsConf
  // ***************
  void AnsysFile::ReadBCsConf( StdVector<std::string> &levels ) {

    ENTER_FCN( "Ansys::ReadBCsConf" );
    
    Integer numbc;
    ReadMaxnumnodesbc(numbc);

    std::string::size_type pos=0;
    getPosLine("Node BC", pos);
    infile.seekg(pos,std::ios::beg);
    
    std::string str;

    Integer nodalnum;
    Integer i,j;
    for ( i = 0; i < numbc; i++ ) {

      infile >> nodalnum >> str;
      infile.ignore(100,'\n');
        
      if ( i == 0 ) {
        levels.Push_back(str);
      }
      else {
        Integer find = 0;
        for (j=0; j<levels.GetSize(); j++)
          if (str == levels[j]) find = 1;
        if (!find) levels.Push_back(str);         
      }
    } 
  }


  // **********
  //   ReadEl
  // **********
  void AnsysFile::ReadEl( StdVector<Elem*> *allElems, 
                          StdVector<Elem*> &orderedElems,
                          const StdVector<std::string> sd ) {

    ENTER_FCN( "AnsysFile::ReadEl" );

    switch( dim_ ) {
    case 2:
      ReadElementInfoFromMeshFile( "2D", allElems, orderedElems, sd );
      break;
    case 3:
      ReadElementInfoFromMeshFile( "3D", allElems, orderedElems, sd );
      break;
    }
  }


  // *******************************
  //   ReadElementInfoFromMeshFile
  // *******************************
  void AnsysFile::
  ReadElementInfoFromMeshFile( std::string elemType,
                               StdVector<Elem*> *elemVec, 
                               StdVector<Elem*> &elemVecSeq,
                               const StdVector<std::string> sd ) {

    ENTER_FCN( "AnsysFile::ReadElementInfoFromMeshFile" );

    // Check that task specifier is correct
    if ( elemType != "1D" && elemType != "2D" && elemType != "3D" ) {
      (*error) << "AnsysFile::ReadElementInfoFromMeshFile: Task specifier "
               << "is '" << elemType << "', but should be one of '1D', '2D' "
               << "or '3D'";
      Error( __FILE__, __LINE__ );
    }

    // Check that pointers to base elements are initialised
    Boolean baseElemInit = TRUE;
    if ( (elemType == "1D" && ( !ptTet1 || !ptHexa1 || !ptHexa2 || !ptPyra1 ||
                                !ptWedge1 || !ptWedge2 ) ) ||
         (elemType == "2D" && ( !ptQ1 || !ptQ2 || !ptTr1 || !ptTr2 ) ) ||
         (elemType == "3D" && ( !ptTet1 || !ptHexa1 || !ptHexa2 || !ptPyra1 ||
                                !ptWedge1 || !ptWedge2 ) ) ) {
      (*error) << "Pointers to " << elemType << " base elements are not "
               << "completely initialized";
        Error( __FILE__, __LINE__ );
    }

    // This string is used for assembling keywords that contain the
    // task specifier elemType
    std::string searchString;

    // Determine the number of elements of respective dimension from
    // the header of the mesh-file
    Integer numElems;
    searchString = "Num" + elemType + "Elements";
    ReadMaxnumelem( numElems, searchString );

    // If there are no elements, we assume that this is fine an
    // simply return
    if ( numElems <= 0 ) {
      return;
    }

    // Check whether for the vector containing all elements (independent of
    // dimension) space has already been reserved
    if ( elemVecSeq.Capacity() < maxNumElems_ ) {
      elemVecSeq.Resize( maxNumElems_ );
    }

    // We need some strings for navigating the mesh-file
    std::string::size_type pos = 0;
    std::string::size_type lineEndPos = 0;
    std::string buf;

    // Position ourselves in the correct setion
    searchString.clear();
    searchString = elemType + " Elements";
    getPosLine( searchString, pos );
    infile.seekg( pos, std::ios::beg );

    // Some additional variables
    Integer i, j, k, eNum, eType, eNodes;
    std::string region;
    Elem *el = NULL;

    // Loop over all elements
    for ( i = 0; i < numElems; i++ ) {

      // Remember current position and get the position of endline
      pos = infile.tellg();
      std::getline( infile, buf, '\n' );
      lineEndPos = infile.tellg();
      infile.seekg( pos, std::ios::beg );

      // try to read data
      infile >> eNum >> eType >> eNodes >> region;       

      // if read in was successfull, enline position and current
      // position are the same
      infile.ignore( 100, '\n' );
      pos = infile.tellg();
      if ( pos != lineEndPos ) {
        (*error) << "An error occured while reading the " << i << "-th "
                 << elemType << " element";
        Error( __FILE__, __LINE__ );
      }

      // Check number of element
      if ( eNum > actMaxElemNum_ ) {
        (*error) << "Current element number = " << eNum << " > "
                 << actMaxElemNum_ << " = actMaxElemNum_. Something might "
                 << "have gone wrong in the meshing process.";
        Error( __FILE__, __LINE__ );
      }

      // Generate new element and insert basic information
      el = new Elem();
      el->elemNum = eNum;
      el->ptElem  = Type2ptElem( eType );
      el->namesd  = region;
      el->connect.Resize( eNodes );

      // Read node numbers and insert them into the element
      for ( k = 0; k < eNodes; k++ ) {
        infile >> el->connect[k];
      }

      // Proceed in mesh-file
      infile.ignore( 100, '\n' );
      pos = infile.tellg();

      // Insert new element into pointer vectors
      elemVecSeq[eNum-1] = el;
        
      Boolean Find = FALSE;
      for ( j = 0; j < sd.GetSize(); j++ ) {
        if ( region == sd[j] ) {
          elemVec[j].Push_back( el );
          Find = TRUE;
        }
      }
      if ( Find == FALSE ) {
        (*error) << elemType << " element " << eNum << " belongs to the "
                 << "region '" << region << "'. However, that regions is "
                 << "not mentioned in the xml-file.";
        Error( __FILE__, __LINE__ );
      }
    }

    // Check that there are no more elements
    if ( !IsNextLineEmpty(pos) ) {
      (*warning) << "AnsysFile::ReadElementInfoFromMeshFile: "
                 << "The line after the last " << elemType
                 << " element no. " << eNum << " in region '" << region
                 << "' seems to contain elements too. Please check if the "
                 << "number of " << elemType << " elements specified in "
                 << "the header of the mesh-file matches the real number of "
                 << elemType << " elements!";
      Warning( __FILE__, __LINE__ );
    }
  }


  // ************
  //   ReadEl1D
  // ************
  void AnsysFile::ReadEl1d( StdVector<Elem*> *allElems, 
                            StdVector<Elem*> &orderedElems,
                            const StdVector<std::string> sd ) {
    ENTER_FCN( "AnsysFile::ReadEl1D" );
    ReadElementInfoFromMeshFile( "1D", allElems, orderedElems, sd );
  }


  // ************
  //   ReadEl2D
  // ************
  void AnsysFile::ReadEl2d( StdVector<Elem*> *allElems, 
                            StdVector<Elem*> &orderedElems,
                            const StdVector<std::string> sd ) {
    ENTER_FCN( "AnsysFile::ReadEl2D" );
    ReadElementInfoFromMeshFile( "2D", allElems, orderedElems, sd );
  }


  // ************
  //   ReadEl3D
  // ************
  void AnsysFile::ReadEl3d( StdVector<Elem*> *allElems, 
                            StdVector<Elem*> &orderedElems,
                            const StdVector<std::string> sd ) {
    ENTER_FCN( "AnsysFile::ReadEl3d" );
    ReadElementInfoFromMeshFile( "3D", allElems, orderedElems, sd );
  }


  // =========================================================================
  // AUXILLIARY METHODS
  // =========================================================================


  // **************
  //   getPosLine
  // **************
  void AnsysFile::getPosLine( const std::string seekexp,
                              std::string::size_type &pos ) {

    ENTER_IFCN( "AnsysFile::getPosLine" );
    infile.seekg(pos, std::ios::beg);
    std::string buf;
    pos=std::string::npos;
    
    std::string::size_type hpos;
    
    while (pos == std::string::npos && !infile.eof()) {
      hpos=infile.tellg();
      std::getline(infile, buf, '\n');
      pos=buf.find(seekexp);
    }

    pos=infile.tellg();

    if (pos>=pos_end) {
      (*error) << "Cannot find string: " << seekexp << " in your mesh-file.";
      Error( __FILE__, __LINE__ );
    }

    // check, if there are comments lines
    do {
      std::getline(infile, buf, '\n');
      if (buf[0] =='#') pos=infile.tellg();
    }
    while (buf[0] == '#'); 
  }


  // ***************
  //   getPosition
  // ***************
  void AnsysFile::getPosition( const std::string seekexp,
                               std::string::size_type &pos ) {

    ENTER_IFCN( "AnsysFile::getPosition" );

    infile.seekg(pos, std::ios::beg);
    std::string buf;
    pos=std::string::npos;
    
    std::string::size_type hpos;

    while (pos == std::string::npos && !infile.eof()) {
      hpos=infile.tellg();
      std::getline(infile, buf, '\n');
      pos=buf.find(seekexp);
    }
    pos+=hpos+seekexp.length();

    if (pos>=pos_end) {
      (*error) << "Cannot find string: " << seekexp << " in your mesh-file.";
      Error( __FILE__, __LINE__ );
    }
  }


  // **************
  //   GetInteger
  // **************
  Integer AnsysFile::GetInteger( std::string seekexp ) {

    ENTER_IFCN( "AnsysFile::GetInteger" );

    std::string::size_type pos = 0;
    std::string::size_type lineEndPos = 0;
    Integer val;
    std::string buf;

    getPosition(seekexp, pos);
    infile.seekg(pos,std::ios::beg);

    // remember current position and get the position of endline
    std::getline(infile,buf,'\n');
    lineEndPos=infile.tellg();
    infile.seekg(pos,std::ios::beg);
  
    // try to read data
    infile >> val;
  
    // if read in was successfull, endline position and current
    // position are the same
    infile.ignore(100,'\n');
    pos = infile.tellg();
    if ( pos != lineEndPos ) {
      (*error) << "AnsysFile::GetInteger: The value for '" << seekexp
               << "' could not be read. Please check your mesh-file";
      Error( __FILE__, __LINE__ );
    }
    return val;
  }


  // *******************
  //   IsNextLineEmpty
  // *******************
  Boolean AnsysFile::IsNextLineEmpty( std::string::size_type actPos ) {

    ENTER_IFCN( "AnsysFile::IsNextLineEmpty" );

    std::string buf = "------";
  
    infile.seekg(actPos,std::ios::beg);  
    std::getline(infile,buf,'\n');
    Boolean retVal;
    if ( buf == "" ) {
      retVal = TRUE;
    }
    else {
      retVal = FALSE;
    }
    return retVal;
  }


  // =========================================================================
  // MISCELLANEOUS METHODS
  // =========================================================================


  // ***************
  //   Type2ptElem
  // ***************
  BaseFE * AnsysFile::Type2ptElem( const Integer itype ) {

    ENTER_IFCN( "AnsysFile::Type2ptElem" );

    BaseFE *retVal = NULL;

    switch( itype ) {

    case 101:
      retVal = ptL2;
      break;
    case 100:
      retVal = ptL1;
      break;
    case 4:
      retVal = ptTr1;
      break;
    case 5:
      retVal = ptTr2;
      break;
    case 6:
      retVal = ptQ1;
      break;
    case 7:
      retVal = ptQ2;
      break;
    case 8:
      retVal = ptTet1;
      break;
    case 10:
      retVal = ptHexa1;
      break;
    case 11:
      retVal = ptHexa2;
      break;
    case 12:
      retVal = ptPyra1;
      break;
      // case 13:
      // retVal = ptPyra2;
    case 14:
      retVal = ptWedge1;
      break;
    case 15:
      retVal = ptWedge2;
      break;
    default:
      (*error) << "AnsysFile::Type2ptElem: No conversion defined for "
               << "element type no. " << itype;
      Error( __FILE__,__LINE__ );
    }

    // Return what we found
    return retVal;
  }


  // ****************
  //   ReadEl3dConf
  // ****************
  void AnsysFile::ReadEl3dConf( StdVector<std::string> &sd ) {

    ENTER_FCN( "AnsysFile::ReadEl3dConf" );

    sd.Clear();
    Integer maxnelems;
    ReadMaxnumelem(maxnelems,"Num3DElements");

    std::string::size_type pos=0;
    getPosLine("3D Elements", pos);
    infile.seekg(pos,std::ios::beg);

    if (!ptTet1 || !ptHexa1 || !ptHexa2 || !ptPyra1 || !ptWedge1)
      Error(" Pointers to BaseElem is not initialized",__FILE__,__LINE__);

    Integer i, ii, j, inum, itype, innodes;
    std::string namesd;

    for (i=0; i<maxnelems; i++)
      {
        Elem * el=new Elem();
        infile >> inum >> itype >> innodes >> namesd;
        infile.ignore(100,'\n');

        el->elemNum=inum;
        if (inum > actMaxElemNum_) {
          (*error) << "The current element number is higher than the "
                   << "maximum number of elements in your .mesh-file. "
                   << "Something might have gone wrong in the meshing "
                   << "process.";
          Error(  __FILE__, __LINE__ );
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


  // ****************
  //   ReadEl2dConf
  // ****************
  void AnsysFile::ReadEl2dConf( StdVector<std::string> &sd ) {

    ENTER_FCN( "AnsysFile::ReadEl2dConf" );

    sd.Clear();
    Integer maxnelems;
    ReadMaxnumelem(maxnelems,"Num2DElements");

    std::string::size_type pos=0;
    getPosLine("2D Elements", pos);
    infile.seekg(pos,std::ios::beg);

    if (!ptQ1 || !ptQ2 || !ptTr1 || !ptTr2)
      Error(" Pointers to BaseElem is not initialized",__FILE__,__LINE__);

    Integer i, ii, j, inum, itype, innodes;
    std::string namesd;

    for (i=0; i<maxnelems; i++)
      {
        Elem * el=new Elem();
        infile >> inum >> itype >> innodes >> namesd;
        infile.ignore(100,'\n');

        el->elemNum=inum;
        if (inum > actMaxElemNum_) {
          (*error) << "The current element number is higher than the "
                   << "maximum number of elements in your .mesh-file. "
                   << "Something might have gone wrong in the meshing "
                   << "process.";
          Error(  __FILE__, __LINE__ );
        }
        el->ptElem=Type2ptElem(itype);
        el->namesd = namesd;
        el->connect.Resize(innodes);
        for (ii=0; ii<innodes; ii++)
          infile >> el->connect[ii];

        infile.ignore(100,'\n');
        //      allelems.Push_back(el);
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


  // ****************
  //   ReadEl1dConf
  // ****************
  void AnsysFile::ReadEl1dConf( StdVector<std::string> &sd ) {

    ENTER_FCN( "AnsysFile::ReadEl1dConf" );

    sd.Clear();
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
        if (inum > actMaxElemNum_) {
          (*error) << "The current element number is higher than the "
                   << "maximum number of elements in your .mesh-file. "
                   << "Something might have gone wrong in the meshing "
                   << "process.";
          Error(  __FILE__, __LINE__ );
        }
        el->ptElem=Type2ptElem(itype);
        el->namesd = namesd;
        el->connect.Resize(innodes);
        for (ii=0; ii<innodes; ii++)
          infile >> el->connect[ii];

        infile.ignore(100,'\n');
        //      allelems.Push_back(el);
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



  // =========================================================================
  // The following methods are concerned with grid adaptation
  // =========================================================================

#ifdef ADAPTGRID

  void AnsysFile::ReadBCs_GridRG( std::list<Integer> &idBCs,
                                  StdVector<Integer> &colorBCs ) {

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
        for (j=0; j<levelsGetSize(); j++) {
          if (str==levels[j]) {
            Find=TRUE;
            break;
          }
        }
        std::string msg=str+"-this level of BCs from .mesh file is not ";
        msg += "mentioned in .xml file. Please, check .xml-file";
        if (!Find) Error(msg.c_str(),__FILE__,__LINE__);
        
        idBCs.Push_back(nodalnum);
        colorBCs.Push_back(j);
      } 
  
  }

  void AnsysFile::ReadGrid_RG(StdVector<grd::Element*> & elems,
                              StdVector<grd::Vertex*> * vertex,
                              const StdVector<std::string> sd)
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
                {
                  (*error) << "This type of elems in mesh file is not "
                           << "implemented yet";
                  Error( __FILE__, __LINE__ );
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
                {
                  (*error) << "This type of elems in mesh file is not "
                           << "implemented yet";
                  Error( __FILE__, __LINE__ );
                }
              }
            infile.ignore(100,'\n'); 
          }

      } // end of if
  }


  // ************
  //   SetNumSD
  // ************
  void AnsysFile::SetNumSD( grd::Element *ptEl, const std::string namesd,
                            const StdVector<std::string> sd ) {

    ENTER_FCN( "AnsysFile::SetNumSD" );

    Boolean Find;
    Integer j;
    for (j=0; j<sdGetSize(); j++)
      if (namesd == sd[j]) { ptEl->setValue(j);
      Find=TRUE;
      }
    if (!Find) {
      (*error) << namesd << "- this level of element is not mentioned in "
               << ".xml-file. Please, check .xml-file";
      Error( __FILE__, __LINE__ );
    }
  }

#endif

}
