#include <fstream>
#include <iostream>
#include <string>
#include <stdarg.h>

#include "ansysfile.hh"
#include "DataInOut/WriteInfo.hh"
#include "Utils/StdVector.hh"
#include "DataInOut/CommandLine/BaseCommandLineHandler.hh"

#ifdef ADAPTGRID
#include "DataInOut/ParamHandling/ConfFile.hh"
#include "Triangle.h"
#include "Tetrahedron.h"
#include "Hexahedron.h"
#endif


namespace CoupledField {


  AnsysFile::AnsysFile( const Char *const afilename ) : FileType(afilename) {

    ENTER_FCN( "AnsysFile::AnsysFile" );

    inFile_.open( afilename );
    if ( !inFile_.good() ) {
      (*error) << "AnsysFile: I am unable to open mesh file '"
               << afilename << "!";
      Error( __FILE__, __LINE__ );
    }

    inFile_.seekg(0, std::ios::end);
    pos_end = inFile_.tellg();

    dim_ = GetDim();
    elemDimReadIn_.Resize(dim_);
    elemDimReadIn_.Init(FALSE);
    maxNumElems_ = GetNumElems();
    maxNumNodes_ = GetNumNodes();
  }


  AnsysFile::~AnsysFile() {
    ENTER_FCN( "AnsysFile::~AnsysFile" );
    inFile_.close() ;
  }


  // ======================================================
  // GENERAL MESH INFORMATION
  // ======================================================
  UInt AnsysFile::GetDim() {
    ENTER_FCN( "AndsysFile::GetDim ");
    return GetInteger("Dimension");
  }
  
  UInt AnsysFile::GetNumNodes(){
    ENTER_FCN( "AndsysFile::GetNumNodes ");
    return GetInteger("NumNodes");
  }
    
  UInt AnsysFile::GetNumElems(const UInt dim){
    ENTER_FCN( "AndsysFile::GetNumElems ");
    
    UInt numElems = 0;
    std::string search;
    

    // 1.) return number of all elements
    if ( dim == 0) {
      if( GetDim() == 3)
        numElems += GetNumElems(3);
      numElems += GetNumElems(2);
      numElems += GetNumElems(1);
    }  
    else if ( dim >=1 && dim <= 3 ) {
      search = "Num" + Info->GenStr(dim) + "DElements";
      numElems = GetInteger(search);
    }
    else {
      (*error) << "AnsysFile::GetNumElems: Dimension " << dim 
               << " is out of range!";
      Error( __FILE__, __LINE__ );
    }
    return numElems;
  }
  
  UInt AnsysFile::GetNumRegions(){
    ENTER_FCN( "AndsysFile::GetNumRegions ");
    Error( "Not implemented", __FILE__, __LINE__ );
    return 0;
  }

  UInt AnsysFile::GetNumNamedNodes(){
    ENTER_FCN( "AndsysFile::GetNumNamedNodes");
    UInt numNamedNodes = 0;
    
    numNamedNodes += GetInteger("NumNodeBC");
    numNamedNodes += GetInteger("NumSaveNodes");

    return numNamedNodes;
  }

  UInt AnsysFile::GetNumNamedElems(){
    ENTER_FCN( "AndsysFile::GetNumNamedElems");
    return GetInteger("NumSaveElements");
  }
  
 
  // ======================================================
  // ENTITY NAME ACCESS
  // ======================================================

  
  void AnsysFile::GetAllRegionNames( StdVector<std::string> & regionNames ){
    ENTER_FCN( "AnsysFile::GetAllRegionNames" );
    
    StdVector<std::string>  names;
    
    regionNames.Clear();
    
    for ( UInt iDim=dim_; iDim>0; iDim-- ) {
      names.Clear();
      GetRegionNamesOfDim(names,iDim);
      for ( UInt iName=0; iName<names.GetSize(); iName++ )
        regionNames.Push_back(names[iName]);

    }

  }
    
    
  void AnsysFile::GetRegionNamesOfDim( StdVector<std::string> & regionNames,
                                       const UInt dim ) {
    ENTER_FCN( "AnsysFile::GetRegionNamesOfDim" );
    
    regionNames.Clear();

    // Check if elements of desired dimension were read in. If not,
    // read them in into dummy variables
    if ( elemDimReadIn_[dim-1] == FALSE ) {
      StdVector<StdVector<Elem*> > dummyElems;
      StdVector<RegionIdType> dummyId;
      GetElements(dummyElems,dummyId,dim);
    }
    
    // Look for region names of desired dimension
    for ( UInt i=0; i<regionDim_.GetSize(); i++ ) 
      if ( regionDim_[i] == dim )
        regionNames.Push_back( regionNames_[i] );
    
  }

  void AnsysFile::GetNodeNames( StdVector<std::string> &nodeNames ) {

    ENTER_FCN( "AnsysFile::GetNodeNames" );
    

    std::string::size_type pos=0;
    std::string str;
    UInt nodalnum;
    UInt i,j;    
    StdVector<std::string> sections;
    StdVector<UInt> numNamedNodes;
    
    nodeNames.Clear();
    sections = "Node BC", "Save Nodes";
    numNamedNodes.Resize(2);
    numNamedNodes[0] = GetInteger("NumNodeBC");
    numNamedNodes[1] = GetInteger("NumSaveNodes");
    
    for ( UInt iSect=0; iSect<sections.GetSize(); iSect++ ) {
      
      GetPosLine(sections[iSect], pos);
      inFile_.seekg(pos,std::ios::beg);
      
      
      for ( i = 0; i < numNamedNodes[iSect]; i++ ) {
        inFile_ >> nodalnum >> str;
        inFile_.ignore(100,'\n');
        
        Boolean found = FALSE;
        
        for ( j = 0; j < nodeNames.GetSize(); j++ ) {
          if ( str == nodeNames[j] ) {
            found = TRUE;
          }
        }
        if ( !found ) {
          nodeNames.Push_back(str);
        } 
      }
    }
  }
  void AnsysFile::GetElemNames( StdVector<std::string> & elemNames ) {
    ENTER_FCN( "AnsysFile::GetElemNames" );

    std::string::size_type pos=0;
    std::string str;
    UInt elemNum, numNamedElems;
    UInt i,j;    
    
    elemNames.Clear();
    numNamedElems = GetInteger("NumSaveElements");
    
    GetPosLine("[Save Elements]", pos);
    inFile_.seekg(pos,std::ios::beg);
    
    
    for ( i = 0; i < numNamedElems; i++ ) {
      inFile_ >> elemNum >> str;
      inFile_.ignore(100,'\n');
      
      Boolean found = FALSE;
      
      for ( j = 0; j < elemNames.GetSize(); j++ ) {
        if ( str == elemNames[j] ) {
          found = TRUE;
        }
      }
      if ( !found ) {
        elemNames.Push_back(str);
      } 
    }
  }
  

  // ======================================================
  // ENTITY ACCESS
  // ======================================================
  
  void AnsysFile::GetCoordinates( StdVector<Point<2> > & nodeCoords ) {
    ENTER_FCN( "AnsysFile::GetCoordinates 2D" );

    UInt i, ibuf;

    std::string::size_type pos=0;
    
    UInt numNodes = GetNumNodes();
    nodeCoords.Resize(numNodes);


    GetPosLine("[Nodes]", pos);
    inFile_.seekg(pos,std::ios::beg);
  
    for ( i = 0; i < numNodes; i++ ) {
      inFile_ >> ibuf >> nodeCoords[i][0] >> nodeCoords[i][1] ;
      inFile_.ignore(100,'\n');
    }
  }


  void AnsysFile::GetCoordinates( StdVector<Point<3> > & nodeCoords ) {
    ENTER_FCN( "AnsysFile::GetCoordinates 3D" );
    
    UInt i, ibuf;
    std::string::size_type pos=0;
    
    UInt numNodes = GetNumNodes();
    nodeCoords.Resize(numNodes);
    
    
    GetPosLine("[Nodes]", pos);
    inFile_.seekg(pos,std::ios::beg);
    
    for ( i = 0; i < numNodes; i++ ) {
      inFile_ >> ibuf >> nodeCoords[i][0] >> nodeCoords[i][1]
              >> nodeCoords[i][2];
      inFile_.ignore(100,'\n');
    }
  }

  void AnsysFile::GetNodesOfRegions( StdVector<StdVector<UInt> > &nodes,
                                     const StdVector<RegionIdType> & regionId ) {

    ENTER_FCN( "AnsysFile::GetNodesOfRegions" );

    std::set<UInt>::iterator it;
    UInt iRegion, index, iNode;
    
    
    nodes.Resize(regionId.GetSize());

    for ( iRegion = 0; iRegion < regionId.GetSize(); iRegion++ ) {
      
      iNode = 0;
      index = regionId[iRegion];
      nodes[iRegion].Resize(regionNodes_[index].size());

      for (it = regionNodes_[index].begin();it != regionNodes_[index].end();
           it++, iNode++ ) {
        nodes[iRegion][iNode] = *it;
      }
        
    }
    
    
  }
    

  void AnsysFile::GetElements( StdVector< StdVector<Elem*> > & elems, 
                               StdVector<RegionIdType> & regionIds,
                               const UInt dim ) {
    ENTER_FCN( "AnsysFile::GetElements" );
    
    // Check that dimension is correct
    if ( dim < 1 || dim > 3 ) {
      (*error) << "AnsysFile::ReadElements: The dimension of elements to "
               << "be read in was specified with " << dim << ", but is only "
               << "allowed to have a value between 1 and 3!";
      Error( __FILE__, __LINE__ );
    }
    
    // Check that pointers to base elements are initialised

    if ( (dim == 1 && ( !ptTet1 || !ptTet2 || !ptHexa1 || !ptHexa2 ||!ptPyra1
			|| !ptPyra2 || !ptWedge1 || !ptWedge2 ) ) ||
         (dim == 2 && ( !ptQ1 || !ptQ2 || !ptTr1 || !ptTr2 ) ) ||
         (dim == 3 && ( !ptTet1 || !ptTet2 || !ptHexa1 || !ptHexa2 ||!ptPyra1
			|| !ptPyra2 || !ptWedge1 || !ptWedge2 ) ) ) {
      (*error) << "Pointers to " << dim << "D base elements are not "

               << "completely initialized";
      Error( __FILE__, __LINE__ );
    }
    
    // This string is used for assembling keywords that contain the
    // task specifier elemType
    std::string searchString;

    // Determine the number of elements of respective dimension from
    // the header of the mesh-file
    UInt numElems = GetNumElems(dim);

    // If there are no elements, we assume that this is fine an
    // simply return
    if ( numElems <= 0 ) {
      return;
    }

    // We need some strings for navigating the mesh-file
    std::string::size_type pos = 0;
    std::string::size_type lineEndPos = 0;
    std::string buf;

    // Position ourselves in the correct setion
    searchString.clear();
    searchString = Info->GenStr(dim) + "D Elements";
    GetPosLine( searchString, pos );
    inFile_.seekg( pos, std::ios::beg );

    // Some additional variables
    UInt i, k, eNum, eType, eNodes;
    std::string region, lastRegion;
    RegionIdType regionId = NO_REGION_ID;
    Integer regionIndex = 0;
    Elem *el = NULL;
    

    // Loop over all elements
    for ( i = 0; i < numElems; i++ ) {

      // Remember current position and get the position of endline
      pos = inFile_.tellg();
      std::getline( inFile_, buf, '\n' );
      lineEndPos = inFile_.tellg();
      inFile_.seekg( pos, std::ios::beg );

      // try to read data
      inFile_ >> eNum >> eType >> eNodes >> region;       

      // if read in was successfull, enline position and current
      // position are the same
      inFile_.ignore( 100, '\n' );
      pos = inFile_.tellg();
      if ( pos != lineEndPos ) {
        (*error) << "An error occured while reading the " << i << "-th "
                 << dim << "D element";
        Error( __FILE__, __LINE__ );
      }


      // Check number of element
      if ( eNum > maxNumElems_ ) {
        (*error) << "Current element number = " << eNum << " > "
                 << maxNumElems_ << " = actMaxElemNum_. Something might "
                 << "have gone wrong in the meshing process.";
        Error( __FILE__, __LINE__ );
      }


      // Check if previous element had the same id. 
      // If not, obtain new region identifier
      if( lastRegion != region ) {
        lastRegion = region;
        regionId = ObtainRegionId(region, dim);
        
        // Check if region of this type already exists, and if not
        // add new vector
        regionIndex = regionIds.Find(regionId);
        if ( regionIndex == -1 ) {
          regionIds.Push_back(regionId);
          elems.Push_back( StdVector<Elem*>() );
          regionNodes_.Push_back(std::set<UInt>());
          regionIndex = regionIds.GetSize() - 1;
        }
      }

      // Generate new element and insert basic information
      el = new Elem();
      el->elemNum = eNum;
      el->ptElem  = Type2ptElem( eType );
      el->connect.Resize( eNodes );
      el->regionId = regionId;

      // Read node numbers and insert them into the element and
      // into the vector with all node-numbers per region
      for ( k = 0; k < eNodes; k++ ) {
        inFile_ >> el->connect[k];
        regionNodes_[regionId].insert(el->connect[k]);
      }

      // Proceed in mesh-file
      inFile_.ignore( 100, '\n' );
      pos = inFile_.tellg();

      elems[regionIndex].Push_back( el );
    }

    // Check that there are no more elements
    if ( !IsNextLineEmpty(pos) ) {
      (*warning) << "AnsysFile::ReadElements: "
                 << "The line after the last " << dim
                 << "D element no. " << eNum << " in region '" << region
                 << "' seems to contain elements too. Please check if the "
                 << "number of " << dim << "D elements specified in "
                 << "the header of the mesh-file matches the real number of "
                 << dim << "D elements!";
      Warning( __FILE__, __LINE__ );
    }

    // Set flag which indicates, that elements of given dimension
    // were read in
    elemDimReadIn_[dim-1] = TRUE;
  }


 
  void AnsysFile::GetNamedNodes( StdVector<StdVector<UInt> > & nodes,
                                 StdVector<std::string> & nodeNames ) {

    ENTER_FCN( "AnsysFile::GetNamedNodes" );
    


    std::string::size_type pos=0;
    std::string::size_type lineEndPos =0;
    std::string lastName = "";
    Integer lastIndex = 0;
    std::string str, buf, errMsg;
    UInt nodalnum;
    UInt i;

    StdVector<std::string> sections;
    StdVector<UInt> numNamedNodes;
    sections = "Node BC", "Save Nodes";
    numNamedNodes.Resize(2);
    numNamedNodes[0] = GetInteger("NumNodeBC");
    numNamedNodes[1] = GetInteger("NumSaveNodes");
    
    for ( UInt iSect=0; iSect<sections.GetSize(); iSect++) {

      GetPosLine(sections[iSect], pos);
      inFile_.seekg(pos,std::ios::beg);

      for ( i = 0; i < numNamedNodes[iSect]; i++ ) {
        
        // remember current position and get the position of endline
        pos = inFile_.tellg();
        std::getline(inFile_,buf,'\n');
        lineEndPos=inFile_.tellg();
        inFile_.seekg(pos,std::ios::beg);
        
        // try to read in the data
        inFile_ >> nodalnum >> str;
        
        // if read in was successfull, enline position and current
        // position are the same
        inFile_.ignore(100,'\n');
        pos = inFile_.tellg();
        if (pos != lineEndPos) {
          (*error) << "AnsysFile:GetNamedNodes: The node list for the boundary "
                   << "conditions has wrong size or format. Please correct it!";
          Error( __FILE__, __LINE__ );
        }
        
        // get according vector index
        if (str != lastName) {
          lastName = str;
          
          // find the associated level
          lastIndex = nodeNames.Find(str);
          
          if( lastIndex == -1 ) {
            nodeNames.Push_back(str);
            nodes.Push_back( StdVector<UInt>() );
            lastIndex = nodes.GetSize()-1; 
          }
        }
        
        nodes[lastIndex].Push_back(nodalnum);
      } 
      
      if (! IsNextLineEmpty(pos)) {
        (*warning) << "AnsysFile::GetNamedNodes: The line after the last BC"
                   << "node "
                   << "no. " << nodalnum << " in region '" << str
                   << "' seems to contain nodes too. Please check if the "
                   << "number of named nodes specified in the header of the "
                   << "mesh-file matches the real number of BC  nodes!";
        Warning( __FILE__, __LINE__ );
      } // end if 
    } // end for
  }


  void AnsysFile::GetNamedElems( StdVector<StdVector<UInt> > & elems,
                                 StdVector<std::string> & elemNames ) {
    ENTER_FCN( "AnsysFile::GetNamedElems" );
    
    std::string::size_type pos=0;
    std::string::size_type lineEndPos =0;
    std::string lastName = "";
    Integer lastIndex = 0;
    std::string str, buf, errMsg;
    UInt elemNum;
    UInt i;

    UInt numNamedElems = GetInteger("NumSaveElements");

    GetPosLine("Save Elements", pos);
    inFile_.seekg(pos,std::ios::beg);
    
    for ( i = 0; i < numNamedElems; i++ ) {
      
      // remember current position and get the position of endline
      pos = inFile_.tellg();
      std::getline(inFile_,buf,'\n');
      lineEndPos=inFile_.tellg();
      inFile_.seekg(pos,std::ios::beg);
      
      // try to read in the data
      inFile_ >> elemNum >> str;
      
      // if read in was successfull, enline position and current
      // position are the same
      inFile_.ignore(100,'\n');
      pos = inFile_.tellg();
      if (pos != lineEndPos) {
        (*error) << "AnsysFile:ReadNamedElems: The node list for the "
                 << "boundary "
                 << "conditions has wrong size or format. Please correct it!";
        Error( __FILE__, __LINE__ );
      }
      
      // get according vector index
      if (str != lastName) {
        lastName = str;
        
        // find the associated level
        lastIndex = elemNames.Find(str);
        
        if( lastIndex == -1 ) {
          elemNames.Push_back(str);
          elems.Push_back( StdVector<UInt>() );
          lastIndex = elems.GetSize()-1; 
        }
      }
      
      elems[lastIndex].Push_back(elemNum);
    } 
    
    if (! IsNextLineEmpty(pos)) {
      (*warning) << "AnsysFile::ReadNamedEles: The line after the last "
                 << "named element "
                 << "no. " << elemNum << " in region '" << str
                 << "' seems to contain nodes too. Please check if the "
                 << "number of BC nodes specified in the header of the "
                 << "mesh-file matches the real number of named elems!";
      Warning( __FILE__, __LINE__ );
    } // end if 
      
  }
  
 


  // =========================================================================
  // AUXILLIARY METHODS
  // =========================================================================


  // **************
  //   GetPosLine
  // **************
  void AnsysFile::GetPosLine( const std::string seekexp,
                              std::string::size_type &pos ) {

    ENTER_IFCN( "AnsysFile::GetPosLine" );
    inFile_.seekg(pos, std::ios::beg);
    std::string buf;
    pos=std::string::npos;
    Boolean found = FALSE;
    
    std::string::size_type hpos;
    
    while (found == FALSE && !inFile_.eof()) {
      hpos=inFile_.tellg();
      std::getline(inFile_, buf, '\n');
      pos=buf.find(seekexp);

      if ( pos != std::string::npos) {
        found = TRUE;
      }
    }

    pos=inFile_.tellg();

    if (pos>=pos_end && found == FALSE) {
      (*error) << "Cannot find string: " << seekexp << " in your mesh-file.";
      Error( __FILE__, __LINE__ );
    }

    // check, if there are comments lines
    do {
      std::getline(inFile_, buf, '\n');
      if (buf[0] =='#') pos=inFile_.tellg();
    }
    while (buf[0] == '#'); 
    
    // reset file pointer
    inFile_.seekg(0, std::ios::beg);
  } 


  // ***************
  //   GetPosition
  // ***************
  void AnsysFile::GetPosition( const std::string seekexp,
                               std::string::size_type &pos ) {

    ENTER_IFCN( "AnsysFile::GetPosition" );

    inFile_.seekg(pos, std::ios::beg);
    std::string buf;
    pos=std::string::npos;
    Boolean found = FALSE;
    std::string::size_type hpos;

    while ( found == FALSE && !inFile_.eof() ) {
      hpos=inFile_.tellg();
      std::getline(inFile_, buf, '\n');
      
      pos=buf.find(seekexp);


      if ( pos != std::string::npos ) {
        found = TRUE;
      }
    }


    pos+=hpos+seekexp.length();

    if ( pos>=pos_end && found == FALSE ) {
      (*error) << "Cannot find string: " << seekexp << " in your mesh-file.";
      Error( __FILE__, __LINE__ );
    }

    // set file pointer to beginning
    inFile_.seekg(0, std::ios::beg);
   
  }


  // **************
  //   GetInteger
  // **************
  UInt AnsysFile::GetInteger( std::string seekexp ) {

    ENTER_IFCN( "AnsysFile::GetInteger" );

    std::string::size_type pos = 0;
    std::string::size_type lineEndPos = 0;
    UInt val;
    std::string buf;

    GetPosition(seekexp, pos);
    inFile_.seekg(pos,std::ios::beg);

    // remember current position and get the position of endline
    std::getline(inFile_,buf,'\n');
    lineEndPos=inFile_.tellg();
    inFile_.seekg(pos,std::ios::beg);
  
    // try to read data
    inFile_ >> val;
  
    // if read in was successfull, endline position and current
    // position are the same
    inFile_.ignore(100,'\n');
    pos = inFile_.tellg();
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
  
    inFile_.seekg(actPos,std::ios::beg);  
    std::getline(inFile_,buf,'\n');
    inFile_.seekg(actPos,std::ios::beg);  
 
    Boolean retVal;
    if ( buf == "" ) {
      retVal = TRUE;
    }
    else {
      retVal = FALSE;
    }
    return retVal;
  }

  // *******************
  //   ObtainRegionId
  // *******************
  
  RegionIdType AnsysFile::ObtainRegionId( const std::string & regionName, 
                                          const UInt dim ) {
    ENTER_FCN( "AnsysFile::ObtainRegionId" );

    Integer index = regionNames_.Find(regionName);
    if( index == -1 ) {
      regionNames_.Push_back(regionName);
      regionDim_.Push_back(dim);
      index = regionNames_.GetSize() - 1;
    }
    
    // remember, of what dimension this region is
    return index;
      
  }
  // =========================================================================
  // MISCELLANEOUS METHODS
  // =========================================================================


  // ***************
  //   Type2ptElem
  // ***************
  BaseFE * AnsysFile::Type2ptElem( const UInt itype ) {

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
    case 9:
      retVal = ptTet2;
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
    case 13:
      (*warning) << "Pyram. with quadratic shape functions" << 
	"do not work well for some cases "<< 
	"(i.e. electric field intensity). Please verify your results."; 
      Warning(__FILE__, __LINE__);
      retVal = ptPyra2;
      break;
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
  
  // =========================================================================
  // The following methods are concerned with grid adaptation
  // =========================================================================

#ifdef ADAPTGRID

  void AnsysFile::ReadBCs_GridRG( std::list<Integer> &idBCs,
                                  StdVector<Integer> &colorBCs ) {

    ENTER_FCN( "AnsysFile::ReadGrid_RG" );
  
    StdVector<std::string> levels;
    conf->getliststr("list_nodes",levels);
    
    UInt numbc;
    ReadMaxnumnodesbc(numbc);

    std::string::size_type pos=0;
    GetPosLine("Node BC", pos);
    inFile_.seekg(pos,std::ios::beg);
    
    std::string str;
  
    UInt nodalnum;
    UInt i,j;
    for (i=0; i < numbc; i++)
      {
        inFile_ >> nodalnum >> str;
        inFile_.ignore(100,'\n');
        
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

    UInt maxnelems;
    ReadMaxnumelem(maxnelems,"Num2DElements");

    if (maxnelems)
      {
        std::string::size_type pos=0;

        GetPosLine("2D Elements", pos);
        inFile_.seekg(pos,std::ios::beg);

        UInt i, ii, ibuf, itype, innodes;
        std::string namesd;
        UInt connect[4]; 

        elems.resize(maxnelems);
        for (i=0; i<maxnelems; i++)
          {
            inFile_ >> ibuf >> itype >> innodes >> namesd;
            inFile_.ignore(100,'\n');

            for (ii=0; ii<innodes; ii++)
              inFile_ >> connect[ii];

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
                inFile_.ignore(100,'\n'); 
              } 
          } // end of if
      }

    void AnsysFile::ReadEl4AdaptGrid3d(StdVector<grd::Element*> & elems, 
                                       StdVector<grd::Vertex*> * vertex,  
                                       const StdVector<std::string> sd)
    {
      ENTER_FCN( "AnsysFile::ReadElems4AdaptGrid3d" );

      UInt maxnelems;
      ReadMaxnumelem(maxnelems,"Num3DElements");

      if (maxnelems)
        {
          std::string::size_type pos=0;

          GetPosLine("3D Elements", pos);
          inFile_.seekg(pos,std::ios::beg);

          UInt i, ii, ibuf, itype, innodes;
          std::string namesd;
          UInt connect[4]; 

          elems.resize(maxnelems);
          for (i=0; i<maxnelems; i++)
            {
              inFile_ >> ibuf >> itype >> innodes >> namesd;
              inFile_.ignore(100,'\n');

              for (ii=0; ii<innodes; ii++)
                inFile_ >> connect[ii];

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
              inFile_.ignore(100,'\n'); 
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
      UInt  j;
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
