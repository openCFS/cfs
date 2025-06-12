// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <algorithm>
#include <cstdarg>

#include "SimInputMESH.hh"
#include "DataInOut/Logging/LogConfigurator.hh"

// declare class specific logging stream
DEFINE_LOG(mesh, "mesh")

namespace CoupledField {


  SimInputMESH::SimInputMESH(std::string fileName, PtrParamNode inputNode, PtrParamNode infoNode) :
    SimInput(fileName, inputNode, infoNode )
  {
    capabilities_.insert( SimInput::MESH);
  }


  SimInputMESH::~SimInputMESH()
  {
    inFile_.close();
  }


  void SimInputMESH::InitModule()
  {
    
    inFile_.open( fileName_.c_str(), std::ios::binary );
    if ( !inFile_.good() ) {
      EXCEPTION("Cannot open mesh file " << fileName_);
    }

    inFile_.seekg(0, std::ios::end);
    pos_end = inFile_.tellg();

    dim_ = GetDim();
    elemDimReadIn_.resize(dim_);
    maxNumElems_ = GetNumElems();
    maxNumNodes_ = GetNumNodes();
  }
    
  void SimInputMESH::ReadMesh(Grid *mi)
  {

    mi_ = mi;
    
    // Get Nodes
    UInt numNodes = GetNumNodes();
    mi_->AddNodes(numNodes);

    std::vector< double > coords;
    GetCoordinates( coords );

    Vector<Double> loc(3);
    for(UInt i=0; i<numNodes; i++) {
      loc[0] = coords[i*3+0];
      loc[1] = coords[i*3+1];
      loc[2] = coords[i*3+2];
      mi_->SetNodeCoordinate(i+1, loc);
    }
      

    // Get Regions
    StdVector<std::string> regionNames;
    StdVector<std::string> names;
    StdVector<Integer> ids;

    GetRegionNamesOfDim(names, 1);
    for(UInt i = 0; i < names.GetSize(); i++)
      regionNames.Push_back(names[i]);
    
    GetRegionNamesOfDim(names, 2);
    for(UInt i = 0; i < names.GetSize(); i++)
      regionNames.Push_back(names[i]);

    if(GetDim() == 3)
    {
      GetRegionNamesOfDim(names, 3);
      for(UInt i = 0; i < names.GetSize(); i++)
        regionNames.Push_back(names[i]);
    }

    mi_->AddRegions(regionNames, ids);

    // Get Elements
    UInt numElems = 0;
    UInt numElems1D = GetNumElems(1);
    UInt numElems2D = GetNumElems(2);
    UInt numElems3D = GetNumElems(3);
    numElems = numElems1D + numElems2D + numElems3D;

    mi_->AddElems(numElems);

    std::vector< std::vector<UInt> > elems;
    std::vector< std::vector<UInt> > elemNums;
    std::vector< std::vector<Elem::FEType> > elemTypes;
    std::vector<RegionIdType> regionId;

    GetElements(elems,elemTypes,elemNums,regionId,1);
    
    UInt n;

    for(UInt j=0; j<elems.size(); j++)
    {
      n=0;
      for(UInt i=0; i<elemTypes[j].size(); i++)
      {
        mi_->SetElemData(elemNums[j][i], elemTypes[j][i], 
                         ids[regionId[j]], &elems[j][n]);
        n += Elem::shapes[elemTypes[j][i]].numNodes;
      }
    }

    elems.clear();
    elemTypes.clear();
    elemNums.clear();
    regionId.clear();
    GetElements(elems,elemTypes,elemNums,regionId,2);
    for(UInt j=0; j<elems.size(); j++)
    {
      n=0;
      for(UInt i=0; i<elemTypes[j].size(); i++)
      {
        mi_->SetElemData(elemNums[j][i], elemTypes[j][i], 
                         ids[regionId[j]], &elems[j][n]);
        n += Elem::shapes[elemTypes[j][i]].numNodes;
      }
    }

    elems.clear();
    elemTypes.clear();
    elemNums.clear();
    regionId.clear();
    GetElements(elems,elemTypes,elemNums,regionId,3);
    n=0;
    for(UInt j=0; j<elems.size(); j++)
    {
      n=0;
      for(UInt i=0; i<elemTypes[j].size(); i++)
      {
        mi_->SetElemData(elemNums[j][i], elemTypes[j][i], 
                         ids[regionId[j]], &elems[j][n]);
        n += Elem::shapes[elemTypes[j][i]].numNodes;
      }
    }


    // Get Named Nodes
    StdVector<StdVector<UInt> > indices;

    names.Clear();
    GetNamedNodes(indices, names);

    for(UInt i = 0; i < names.GetSize(); ++i)
      mi_->AddNamedNodes(names[i], indices[i]);


    // Get Named Elements
    names.Clear();
    indices.Clear();
    
    GetNamedElems(indices, names);

    for(UInt i = 0; i < names.GetSize(); ++i)
      mi_->AddNamedElems(names[i], indices[i]);

  }

  // ======================================================
  // GENERAL MESH INFORMATION
  // ======================================================
  UInt SimInputMESH::GetDim() {
    return GetInteger("Dimension");
  }
  
  UInt SimInputMESH::GetNumNodes(){
    return GetInteger("NumNodes");
  }
    
  UInt SimInputMESH::GetNumElems(const Integer dim){
    
    UInt numElems = 0;
    std::stringstream search;
    

    // 1.) return number of all elements
    if ( dim == 0) {
      if( GetDim() == 3)
        numElems += GetNumElems(3);
      numElems += GetNumElems(2);
      numElems += GetNumElems(1);
    }  
    else if ( dim >=1 && dim <= 3 ) {
      search << "Num" << dim;
      search << "DElements";
      numElems = GetInteger(search.str());
    }
    else {
      EXCEPTION("Dimension " << dim << " is out of range!");
    }
    return numElems;
  }
  
  UInt SimInputMESH::GetNumRegions(){
    if(regionNames_.size() == 0)
    {
      StdVector<std::string> names;

      GetAllRegionNames(names);
    }
    return regionNames_.size();
  }

  UInt SimInputMESH::GetNumNamedNodes(){
    UInt numNamedNodes = 0;
    
    numNamedNodes += GetInteger("NumNodeBC");
    numNamedNodes += GetInteger("NumSaveNodes");

    return numNamedNodes;
  }

  UInt SimInputMESH::GetNumNamedElems(){
    return GetInteger("NumSaveElements");
  }
  
  // ======================================================
  // ENTITY NAME ACCESS
  // ======================================================

  void SimInputMESH::GetAllRegionNames( StdVector<std::string> & regionNames ){
    
    if(regionNames_.size() == 0)
    {
      StdVector<std::string>  names;

      regionNames.Clear();

      for ( UInt iDim=dim_; iDim>0; iDim-- ) {
        names.Clear();
        GetRegionNamesOfDim(names,iDim);
        for ( UInt iName=0; iName<names.GetSize(); iName++ )
          regionNames.Push_back(names[iName]);

      }
    }
    else
      regionNames = regionNames_;


    ///////////////////////////////////////////////////////
    /*
    std::vector<std::string> strs;
    std::vector<std::string>::iterator it, end;
      GetNodeNames( strs );
      std::cout << "Named Nodes: " << std::endl;
      it = strs.begin();
      end = strs.end();
      for(; it != end; it++)
      {
      std::cout << "  " << (*it) << std::endl;
      }

      GetElemNames( strs );
      std::cout << "Named Elems: " << std::endl;
      it = strs.begin();
      end = strs.end();
      for(; it != end; it++)
      {
      std::cout << "  " << (*it) << std::endl;
      }

      std::vector< double > coords;
      std::vector<double>::iterator itc, endc;
      GetCoordinates( coords );
      std::cout << "Coordinates: " << std::endl;
      itc = coords.begin();
      endc = coords.end();
      for(; itc != endc; )
      {
      std::cout << "  " << (*itc) << " "; itc++;
      std::cout << (*itc) << " "; itc++;
      std::cout << (*itc) << std::endl; itc++;
      }
    */

    ///////////////////////////////////////////////////////

  }
    
  void SimInputMESH::GetRegionNamesOfDim( StdVector<std::string> & regionNames,
                                       const UInt dim ) {
    
    regionNames.Clear();

    // Check if elements of desired dimension were read in. If not,
    // read them in into dummy variables
    if ( elemDimReadIn_[dim-1] == false ) {
      std::vector< std::vector<UInt> > elems, elemNums;
      std::vector< std::vector<Elem::FEType> > elemTypes;
      std::vector<RegionIdType> dummyId;
      GetElements(elems,elemTypes,elemNums,dummyId,dim);
    }
    
    // Look for region names of desired dimension
    for ( UInt i=0; i<regionDim_.size(); i++ ) 
      if ( regionDim_[i] == dim )
        regionNames.Push_back( regionNames_[i] );
    
  }

  void SimInputMESH::GetNodeNames( StdVector<std::string> &nodeNames ) {

    
    std::string::size_type pos=0;
    std::string str;
    UInt nodalnum;
    UInt i;    
    std::vector<std::string> sections;
    std::vector<UInt> numNamedNodes;
    
    nodeNames.Clear();
    sections.push_back("Node BC");
    sections.push_back("Save Nodes");
    numNamedNodes.resize(2);
    numNamedNodes[0] = GetInteger("NumNodeBC");
    numNamedNodes[1] = GetInteger("NumSaveNodes");
    
    for ( UInt iSect=0; iSect<sections.size(); iSect++ ) {
      
      GetPosLine(sections[iSect], pos);
      inFile_.seekg(pos,std::ios::beg);
      
      
      for ( i = 0; i < numNamedNodes[iSect]; i++ ) {
        inFile_ >> nodalnum >> str;
        inFile_.ignore(100,'\n');


        if ( nodeNames.Find(str) == -1 ) {
          nodeNames.Push_back(str);
        } 
      }
    }
  }

  void SimInputMESH::GetElemNames( StdVector<std::string> & elemNames ) {

    std::string::size_type pos=0;
    std::string str;
    UInt elemNum, numNamedElems;
    UInt i;
    
    elemNames.Clear();
    numNamedElems = GetInteger("NumSaveElements");
    
    GetPosLine("[Save Elements]", pos);
    inFile_.seekg(pos,std::ios::beg);
    
    
    for ( i = 0; i < numNamedElems; i++ ) {
      inFile_ >> elemNum >> str;
      inFile_.ignore(100,'\n');
      
      if ( elemNames.Find(str) == -1 ) {
        elemNames.Push_back(str);
      } 
    }
  }

  // ======================================================
  // ENTITY ACCESS
  // ======================================================
  
  void SimInputMESH::GetCoordinates( std::vector< double > & nodeCoords ) {

    UInt i, ibuf;

    std::string::size_type pos=0;
    
    UInt numNodes = GetNumNodes();
    nodeCoords.resize(numNodes*3);


    GetPosLine("[Nodes]", pos);
    inFile_.seekg(pos,std::ios::beg);
  
    for ( i = 0; i < numNodes; i++ ) {
      inFile_ >> ibuf
              >> nodeCoords[i*3+0]
              >> nodeCoords[i*3+1]
              >> nodeCoords[i*3+2];
      inFile_.ignore(100,'\n');
    }
  }


  void SimInputMESH::GetNodesOfRegions( std::vector<std::vector<UInt> > &nodes,
                                     const std::vector<RegionIdType> & regionId ) {


    std::set<UInt>::iterator it;
    UInt iRegion, index, iNode;
    
    
    nodes.resize(regionId.size());

    for ( iRegion = 0; iRegion < regionId.size(); iRegion++ ) {
      
      iNode = 0;
      index = regionId[iRegion];
      nodes[iRegion].resize(regionNodes_[index].size());

      for (it = regionNodes_[index].begin();it != regionNodes_[index].end();
           it++, iNode++ ) {
        nodes[iRegion][iNode] = *it;
      }
    }
  }
    
  void SimInputMESH::GetElements( std::vector< std::vector<UInt> > & elems,
                                  std::vector< std::vector<Elem::FEType> > & elemTypes,
                                  std::vector< std::vector<UInt> > & elemNums,                                
                                  std::vector<RegionIdType> & regionIds,
                                  const UInt dim ) {
    
    // Check that dimension is correct
    if ( dim < 1 || dim > 3 ) {
      EXCEPTION("The dimension of elements to be read in was specified with "
                << dim << "but is only allowed to have a value between 1 and 3!");
    }
    
    // This string is used for assembling keywords that contain the
    // task specifier elemType
    std::stringstream searchString;

    // Determine the number of elements of respective dimension from
    // the header of the mesh-file
    UInt numElems = GetNumElems(dim);

    // If there are no elements, we assume that this is fine and
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
    searchString << dim;
    searchString << "D Elements";
    GetPosLine( searchString.str(), pos );
    inFile_.seekg( pos, std::ios::beg );

    // Some additional variables
    UInt i, k, eNum, eType, eNodes;
    std::string region, lastRegion;
    RegionIdType regionId = 0;
    Integer regionIndex = 0;
    
    // Loop over all elements
    for ( i = 0; i < numElems; i++ ) {

      // Remember current position and get the position of endline
      pos = inFile_.tellg();
      std::getline( inFile_, buf, '\n' );
      lineEndPos = inFile_.tellg();
      inFile_.seekg( pos, std::ios::beg );

      // try to read data
      inFile_ >> eNum >> eType >> eNodes >> region;       

      // if read in was successful, enline position and current
      // position are the same
      inFile_.ignore( 100, '\n' );
      pos = inFile_.tellg();
      LOG_DBG3(mesh) << "GE dim=" << dim << " numElems=" << numElems << " i=" << i << " pos=" << pos << " eNum=" << eNum << " eType=" << eType << " eNodes=" << eNodes << " region=" << region;
      if(pos != lineEndPos)
        EXCEPTION("An error occurred while reading the " << i << "-th " << dim << "D element");


      // Check number of element
      if(eNum > maxNumElems_)
        EXCEPTION("Current element number = " << eNum << " > " << maxNumElems_
               << " = actMaxElemNum_. Something might have gone wrong in the meshing process.");

      // Check if previous element had the same id. 
      // If not, obtain new region identifier
      if( lastRegion != region ) {
        lastRegion = region;
        regionId = ObtainRegionId(region, dim);
        
        // Check if region of this type already exists, and if not
        // add new vector
        std::vector<RegionIdType>::iterator it, end;
                
        end = regionIds.end();

        it = std::find(regionIds.begin(), end, regionId);
        
        if ( it == end ) {
          regionIds.push_back(regionId);
          elems.push_back( std::vector<UInt>() );
          elemNums.push_back( std::vector<UInt>() );
          elemTypes.push_back( std::vector<Elem::FEType>() );
          regionNodes_.push_back(std::set<UInt>());
          regionIndex = regionIds.size() - 1;
        } else {
          regionIndex = std::distance(regionIds.begin(), it );
        }
      }

      // Generate new element and insert basic information
      //            el = new Elem();
      //            el->elemNum = eNum;
      //            el->ptElem  = Type2ptElem( eType );
      //            el->connect.resize( eNodes );
      //            el->regionId = regionId;

            

      // Read node numbers and insert them into the element and
      // into the vector with all node-numbers per region
      UInt dummy;
      for ( k = 0; k < eNodes; k++ ) {
        inFile_ >> dummy;
        elems[regionIndex].push_back(dummy);
        regionNodes_[regionId].insert(dummy);
      }

      elemTypes[regionIndex].push_back(AnsysType2ElemType(eType));
      elemNums[regionIndex].push_back( eNum );
            
      // Proceed in mesh-file
      inFile_.ignore( 100, '\n' );
      pos = inFile_.tellg();

      //            elems[regionIndex].push_back( el );
    }

    // Check that there are no more elements
    if ( !IsNextLineEmpty(pos) ) {
      EXCEPTION("The line after the last " << dim
                << "D element no. " << eNum << " in region '" << region
                << "' seems to contain elements too. Please check if the "
                << "number of " << dim << "D elements specified in "
                << "the header of the mesh-file matches the real number of "
                << dim << "D elements!");
    }

    // Set flag which indicates, that elements of given dimension
    // were read in
    elemDimReadIn_[dim-1] = true;
  }

  void SimInputMESH::GetNamedNodes(StdVector<StdVector<UInt> > &nodes,
                                   StdVector<std::string> &nodeNames )
  {
    std::string::size_type pos=0;
    std::string::size_type lineEndPos =0;
    std::string lastName = "";
    Integer lastIndex = 0;
    std::string str, buf, errMsg;
    UInt nodalnum;
    UInt i;

    std::vector<std::string> sections;
    std::vector<UInt> numNamedNodes;
    sections.push_back("Node BC");
    sections.push_back("Save Nodes");
    numNamedNodes.resize(2);
    numNamedNodes[0] = GetInteger("NumNodeBC");
    numNamedNodes[1] = GetInteger("NumSaveNodes");
    
    for ( UInt iSect=0; iSect<sections.size(); iSect++) {

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
          EXCEPTION("The node list for the boundary "
                    << "conditions has wrong size or format. Please correct it!");
        }
        
        // get according vector index
        if (str != lastName) {
          lastName = str;
          
          // find the associated level

          StdVector<std::string>::iterator it, end;
                
          end = nodeNames.End();

          it = std::find(nodeNames.Begin(), end, str);
                
          if ( it == end ) {
            nodeNames.Push_back(str);
            nodes.Push_back( StdVector<UInt>() );
            lastIndex = nodes.GetSize()-1; 
          }
          else
          {
            lastIndex = std::distance(nodeNames.Begin(), it);
          }
        }
        
        nodes[lastIndex].Push_back(nodalnum);
      } 
      
      if (! IsNextLineEmpty(pos)) {
        EXCEPTION("The line after the last BC"
                  << "node "
                  << "no. " << nodalnum << " in region '" << str
                  << "' seems to contain nodes too. Please check if the "
                  << "number of named nodes specified in the header of the "
                  << "mesh-file matches the real number of BC  nodes!");
      } // end if 
    } // end for
  }

  void SimInputMESH::GetNamedElems(StdVector<StdVector<UInt> > & elems,
                                   StdVector<std::string> & elemNames)
  {  
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
      
      // if read in was successful, endline position and current
      // position are the same
      inFile_.ignore(100,'\n');
      pos = inFile_.tellg();
      if (pos != lineEndPos) {
        EXCEPTION("The node list for the boundary conditions has wrong size or format. Please correct it! Pos: " << pos << " lineEndPos " << lineEndPos );
      }
      
      // get according vector index
      if (str != lastName) {
        lastName = str;
        
        // find the associated level
        StdVector<std::string>::iterator it, end;
                
        end = elemNames.End();

        it = std::find(elemNames.Begin(), end, str);
                
        if ( it == end ) {
          elemNames.Push_back(str);
          elems.Push_back( StdVector<UInt>() );
          lastIndex = elems.GetSize()-1; 
        }
        else {
          lastIndex = std::distance(elemNames.Begin(), it);
        }
      }
      
      elems[lastIndex].Push_back(elemNum);
    } 
    
    if (! IsNextLineEmpty(pos)) {
      EXCEPTION("The line after the last "
                << "named element "
                << "no. " << elemNum << " in region '" << str
                << "' seems to contain nodes too. Please check if the "
                << "number of BC nodes specified in the header of the "
                << "mesh-file matches the real number of named elems!");
    } // end if 
      
  }
  

  // =========================================================================
  // AUXILLIARY METHODS
  // =========================================================================


  // **************
  //   GetPosLine
  // **************
  void SimInputMESH::GetPosLine( const std::string seekexp,
                              std::string::size_type &pos ) {

    inFile_.seekg(pos, std::ios::beg);
    std::string buf;
    pos=std::string::npos;
    bool found = false;
    
    // std::string::size_type hpos;
    
    while (found == false && !inFile_.eof()) {
      // hpos=inFile_.tellg();
      std::getline(inFile_, buf, '\n');
      pos=buf.find(seekexp);

      if ( pos != std::string::npos) {
        found = true;
      }
    }

    pos=inFile_.tellg();

    if (pos>=pos_end && found == false) {
      EXCEPTION("Cannot find string: "
                << seekexp << " in your mesh-file.");
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
  void SimInputMESH::GetPosition( const std::string seekexp,
                               std::string::size_type &pos ) {


    inFile_.seekg(pos, std::ios::beg);
    std::string buf;
    pos=std::string::npos;
    bool found = false;
    std::string::size_type hpos = 0;

    while ( found == false && !inFile_.eof() ) {
      hpos=inFile_.tellg();
      std::getline(inFile_, buf, '\n');
      
      pos=buf.find(seekexp);


      if ( pos != std::string::npos ) {
        found = true;
      }
    }


    pos+=hpos+seekexp.length();

    if ( pos>=pos_end && found == false ) {
      EXCEPTION("Cannot find string: " << seekexp << " in your mesh-file.");
    }

    // set file pointer to beginning
    inFile_.seekg(0, std::ios::beg);
   
  }

  // **************
  //   GetInteger
  // **************
  UInt SimInputMESH::GetInteger( std::string seekexp ) {


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
      EXCEPTION("The value for " << seekexp
                << " could not be read. Please check your mesh-file");
    }
    return val;
  }

  // *******************
  //   IsNextLineEmpty
  // *******************
  bool SimInputMESH::IsNextLineEmpty( std::string::size_type actPos ) {


    std::string buf = "------";
  
    inFile_.seekg(actPos,std::ios::beg);  
    std::getline(inFile_,buf,'\n');
    inFile_.seekg(actPos,std::ios::beg);  
 
    bool retVal;
    if ( buf == "" || buf == "\r") {
      retVal = true;
    }
    else {
      retVal = false;
    }
    return retVal;
  }

  // *******************
  //   ObtainRegionId
  // *******************
  
  RegionIdType SimInputMESH::ObtainRegionId( const std::string & regionName, 
                                          const UInt dim ) {

    RegionIdType index;
    std::vector< std::string >::iterator it, end;

    end = regionNames_.end();

    it = std::find(regionNames_.begin(),
                   end,
                   regionName);
        
    if( it == end ) {
      regionNames_.push_back(regionName);
      // remember, of what dimension this region is
      regionDim_.push_back(dim);
      index = regionNames_.size() - 1;
    }
    else
      index = std::distance(regionNames_.begin(), it);
    
    return index;
  }

  // =========================================================================
  // MISCELLANEOUS METHODS
  // =========================================================================


  // ***************
  //   Type2ptElem
  // ***************
  Elem::FEType SimInputMESH::AnsysType2ElemType( const UInt itype ) {


    switch( itype ) {

    case 101:
      return Elem::ET_LINE3;
    case 100:
      return Elem::ET_LINE2;
    case 4:
      return Elem::ET_TRIA3;
    case 5:
      return Elem::ET_TRIA6;
    case 6:
      return Elem::ET_QUAD4;
    case 7:
      return Elem::ET_QUAD8;
    case 107:
      return Elem::ET_QUAD9;
    case 8:
      return Elem::ET_TET4;
    case 9:
      return Elem::ET_TET10;
    case 10:
      return Elem::ET_HEXA8;
    case 11:
      return Elem::ET_HEXA20;
    case 111:
      return Elem::ET_HEXA27;
    case 12:
      return Elem::ET_PYRA5;
    case 13:
      //    (*warning) << "Pyram. with quadratic shape functions" << 
      //    "do not work well for some cases "<< 
      //    "(i.e. electric field intensity). Please verify your results."; 
      return Elem::ET_PYRA13;
    case 14:
      return Elem::ET_WEDGE6;
    case 15:
      return Elem::ET_WEDGE15;
    default:
      break;
    }
    
    // This place should never be reached!
    return Elem::ET_UNDEF;
  }

}
