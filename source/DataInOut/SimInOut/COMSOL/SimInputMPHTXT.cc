// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

/*  
   ElmerGrid - A simple mesh generation and manipulation utility  
   Copyright (C) 1995- , CSC - IT Center for Science Ltd.   

   Author: Peter Råback
   Email: Peter.Raback@csc.fi
   Address: CSC - IT Center for Science Ltd.
            Keilaranta 14
            02101 Espoo, Finland

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.
*/

#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <algorithm>
#include <cstdarg>

#include <unzip.h>

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ParamHandling/Xerces.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "SimInputMPHTXT.hh"

namespace CoupledField {

  // declare logging stream
  DECLARE_LOG(simInputMPHTXT)
  DEFINE_LOG(simInputMPHTXT, "simInput.MPHTXT")


#define MAXLINESIZE 200     /* maximum length of line to be read */
#define MAXNODESD3 64       /* maximum number of 3D nodes */ 
#define MAXNODESD2 27       /* maximum number of 2D nodes */ 
#define MAXNODESD1 9        /* maximum number of 1D nodes */
#define MAXFILESIZE 600     /* maximum filenamesize for i/o files */

  int SimInputMPHTXT::next_int(char **start)
  {
    int i;
    char *end;
    
    i = strtol(*start,&end,10);
    *start = end;
    return(i);
  }
  
  Double SimInputMPHTXT::next_real(char **start)
  {
    Double r;
    char *end;
    
    r = strtod(*start,&end);
    
    *start = end;
    return(r);
  }
  
  int SimInputMPHTXT::Comsolrow(char *line1,FILE *io) 
  {
    int i,isend;
    char line0[MAXLINESIZE],*charend;
    
    for(i=0;i<MAXLINESIZE;i++) 
      line0[i] = ' ';
    
    // newline:
    
    charend = fgets(line0,MAXLINESIZE,io);
    isend = (charend == NULL);
    
    if(isend) return(1);
    
    for(i=0;i<MAXLINESIZE;i++) line1[i] = line0[i];    
    
    return(0);
  }
  
  void SimInputMPHTXT::ReorderComsolNodes(int elementtype,UInt *topo) 
  {
    int i,tmptopo[MAXNODESD2];
    int order306[]={1,2,3,4,6,5};
    int order409[]={1,2,4,3,5,8,9,6,7};
    int order510[]={1,2,3,4,5,7,6,8,9,10};
    int order614[]={1,2,4,3,5,6,9,10,7,11,12,14,23,8};
    int order718[]={1,2,3,4,5,6,7,9,8,16,18,17,10,12,15,11,14,13};
    int order827[]={1,2,4,3,5,6,8,7,9,12,13,10,23,26,27,24,14,16,22,20,15,19,21,17,11,25,18};
    
    for(i=0;i<elementtype%100;i++) 
      tmptopo[i] = topo[i];
    
    switch (elementtype) {
      
    case 306:
      for(i=0;i<elementtype%100;i++) 
        topo[i] = tmptopo[order306[i]-1];
      break;
      
    case 404:
    case 409:
      for(i=0;i<elementtype%100;i++) 
        topo[i] = tmptopo[order409[i]-1];
      break;
      
    case 510:
      for(i=0;i<elementtype%100;i++) 
        topo[i] = tmptopo[order510[i]-1];
      break;
      
    case 718:
      for(i=0;i<elementtype%100;i++) 
        topo[i] = tmptopo[order718[i]-1];
      break;
      
    case 605:
    case 614:
      for(i=0;i<elementtype%100;i++) 
        topo[i] = tmptopo[order614[i]-1];
      break;
      
    case 808:
    case 827:
      for(i=0;i<elementtype%100;i++) 
        topo[i] = tmptopo[order827[i]-1];
      break;
      
    case 408:
    case 820:
    case 715:
      EXCEPTION("Reordering for Elmer elemtype " << elementtype << " not yet implemented.");    
      break;
      
    default:
      break;
    }
    
  }
  
  int SimInputMPHTXT::LoadComsolMesh(Grid *data,
                                     const char *prefix,
                                     bool readInfos)
    /* Load the grid in Comsol Multiphysics mesh format */
  {
    int noknots,noelements,noelements_dom,maxnodes,material;
    // int elemcode;
    // int foundsame;
    // int mode,nvalue,maxknot,nosides,sideelemtype;
    int allocated;
    //  int boundarytype,materialtype,boundarynodes,side,parent,elemsides;
    int elemnodes, elembasis, elemtype;
    elembasis = 0;
    elemnodes = 0;
    dim_ = 0;
    // int bulkdone, usedmax,hits;
    // int label;
    int offset,domains,mindom,minbc,elemdim = 0;
    char line[MAXLINESIZE],*cp;
    // int l,n,ind,inds[MAXNODESD2],sideind[MAXNODESD1];
    int i,j,k;
    FILE *in;
    //  Double x,y,z;
    Elem::FEType fetype = Elem::ET_UNDEF;
    StdVector<RegionIdType> ids;
    std::set<UInt> surfDoms, volDoms;
    std::vector<UInt> newSurfRegions, newVolRegions;
    Vector<UInt> globTopo;
    UInt maxNumNodesPerElem = 0;
    
    std::stringstream sstr;
    sstr << prefix; //  << ".mphtxt";
    if ((in = fopen(sstr.str().c_str(),"r")) == NULL) {
      EXCEPTION("LoadComsolMesh: opening of the Comsol mesh file '" <<
                sstr.str() << "' wasn't succesfull!");
    }
    
    LOG_TRACE(simInputMPHTXT) << "Reading mesh from Comsol mesh file "
                              << sstr.str() << std::endl;
    // InitializeKnots(data);
    
    allocated = false;
    
    mindom = 1000;
    minbc = 1000;
    offset = 1;
    
  omstart:
    
    maxnodes = 0;
    noknots = 0;
    noelements = 0;
    noelements_dom = 0;
    // elemcode = 0;
    material = 0;
    domains = 0;
    
    
    for(;;) {
      
      if(Comsolrow(line,in))
        goto end;

      //if(line == NULL) goto end;  clang is right here:  error: comparison of array 'line' equal to a null pointer is always false
      
      if(strstr(line,"# sdim")) {
        cp = line;
        dim_ = next_int(&cp);
        LOG_DBG(simInputMPHTXT) << "dim=" << dim_ << std::endl;
      }
      
      if(strstr(line,"# class")) {
        std::string class_line = line;

        if(class_line.find( "Mesh") == std::string::npos) {
          EXCEPTION("Cannot handle .mphtxt classes other than 'Mesh'.\n"
                    << "Offending class line --> " << class_line);
        }
          
        LOG_DBG(simInputMPHTXT) << "class line=" << class_line << std::endl;
      }

      else if(strstr(line,"# number of mesh points")) {
        cp = line;
        noknots = next_int(&cp);
        LOG_DBG(simInputMPHTXT) << "noknots=" << noknots << std::endl;
      }
      
      else if(strstr(line,"# lowest mesh point index")) {
        cp = line;
        offset = 1 - next_int(&cp);
        LOG_DBG(simInputMPHTXT) << "offset=" << offset << std::endl;
      }
      
      else if(strstr(line,"# type name")) {
        if(strstr(line,"vtx")) elembasis = 100;
        else if(strstr(line,"edg")) elembasis = 200;
        else if(strstr(line,"tri")) elembasis = 300;
        else if(strstr(line,"quad")) elembasis = 400;
        else if(strstr(line,"tet")) elembasis = 500;
        else if(strstr(line,"pyr")) elembasis = 600;
        else if(strstr(line,"prism")) elembasis = 700;
        else if(strstr(line,"hex")) elembasis = 800;
        else printf("unknown element type = %s",line);
      }
      
      else if(strstr(line,"# number of nodes per element")) {
        cp = line;
        elemnodes = next_int(&cp);
        if(elemnodes > maxnodes) maxnodes = elemnodes;      
        LOG_DBG(simInputMPHTXT) << "elemnodes=" << elemnodes << std::endl;
        
        fetype = ElmerType2ElemType(elembasis + elemnodes);
        
      }
      
      else if(strstr(line,"# Mesh point coordinates")) {
        LOG_TRACE(simInputMPHTXT) << "Loading " << noknots <<
          " coordinates" << std::endl;
        
        for(i=1;i<=noknots;i++) {
          Comsolrow(line,in);	
          
          if(allocated) {
            cp = line;
            //	  data->x[i] = next_real(&cp);
            //	  data->y[i] = next_real(&cp);
            //	  if(dim_ == 3) data->z[i] = next_real(&cp);
            
            UInt nodeNum = i;
            Vector<Double> coord(3);
            coord[2] = 0;
            for(UInt j=0; j<dim_; j++) 
            {
              coord[j] = next_real(&cp);
            }
            
            data->SetNodeCoordinate( nodeNum, coord );
          }
        }
      }
      
      else if(strstr(line,"# number of elements")) {
        
        cp = line;
        k = next_int(&cp);
        
        Comsolrow(line,in);	            
        elemtype = elemnodes + elembasis;
        // elemdim = GetElementDimension(elemtype);
        elemdim = Elem::shapes[fetype].dim;
        
        LOG_DBG(simInputMPHTXT) << "Loading " << k << " elements of type " <<
          elemtype << "." << std::endl;

        domains = noelements;
        
        for(i=1;i<=k;i++) {
          Comsolrow(line,in);	
          
          if(dim_ == 3 && elembasis < 300) continue;
          if(dim_ == 2 && elembasis < 200) continue;
          
          noelements = noelements + 1;
          if(allocated) {
            //	  data->elementtypes[noelements] = elemtype;
            //	  data->material[noelements] = 1;
            
            cp = line;
            for(j=0;j<elemnodes;j++)
              globTopo[(noelements-1)*maxNumNodesPerElem + j] = next_int(&cp) + offset;
            //	    data->topology[noelements][j] = next_int(&cp) + offset;
            
            if(1) ReorderComsolNodes(elemtype, &globTopo[(noelements-1)*maxNumNodesPerElem] );
            
          }
        }
      }
      
      else if(strstr(line,"# number of domains") ||
              strstr(line,"# number of geometric entity indices")) {
        
        cp = line;
        k = next_int(&cp);
        
        Comsolrow(line,in);	            
        LOG_DBG(simInputMPHTXT) << "Loading " << k <<
          " domains for the elements." << std::endl;
        
        for(i=1;i<=k;i++) {
          Comsolrow(line,in);	
          
          if(dim_ == 3 && elembasis < 300) continue;
          if(dim_ == 2 && elembasis < 200) continue;
          
          domains = domains + 1;
          cp = line;
          material = next_int(&cp);
          
          noelements_dom = noelements_dom + 1;
          if(allocated) {
            // UInt topo[128];
            RegionIdType id;
            
            // data->GetElemData(noelements_dom, fetype, id, topo);
            
            if(elemdim < (int) dim_) 
            {
              UInt distance = std::distance(surfDoms.begin(), surfDoms.find(material));
              id = ids[ newSurfRegions[distance] ];
              
              material = material - minbc + 1;
              
              data->SetElemData(noelements_dom, fetype, id, &globTopo[(noelements_dom-1)*maxNumNodesPerElem]);
            }
            else
            {
              UInt distance = std::distance(volDoms.begin(), volDoms.find(material));
              id = ids[ newVolRegions[distance] ];
              
              material = material - mindom + 1;
              //	  data->material[domains] = material;
              
              data->SetElemData(noelements_dom, fetype, id, &globTopo[(noelements_dom-1)*maxNumNodesPerElem]);
            }
            
          }
          else {
            if(elemdim < (int) dim_) {
              if(minbc > material) minbc = material;
              surfDoms.insert(material);
            }
            else { 
              if(mindom > material) mindom = material;	  
              volDoms.insert(material);
            }
          }
          
        }
      }      
      else if(strstr(line,"#")) {
        LOG_DBG(simInputMPHTXT) << "Unused command: " << line << "." << std::endl;
      }      
    }
    
  end:
    
    if(!allocated) {
      
      if(noknots == 0 || noelements == 0 || maxnodes == 0) {
        printf("Invalid mesh consits of %d knots and %d %d-node elements.\n",
               noknots,noelements,maxnodes);     
        fclose(in);
        return(2);
      }
      
      rewind(in);
      //    data->noknots = noknots;
      //    data->noelements = noelements;
      //    data->maxnodes = maxnodes;
      //    data->dim = dim_;
      
      LOG_TRACE(simInputMPHTXT) << "Allocating for " << noknots << 
        " knots and " << noelements << " " << maxnodes <<
        "-node elements." << std::endl;
      
      if(readInfos) 
      {
        maxNumNodes_ = noknots;
        maxNumElems_ = noelements;
      }
      else 
      {
        data->AddElems(noelements);
        data->AddNodes(noknots);
        
        // Since the .mphtxt file only contains informations about
        // geometrical entities like points, edges, areas and volumes
        // a large number of these may be mapped to surf_domain_* or
        // vol_domain_* regions without physical meaning in CFS++.
        // Therefore, these entities can be grouped either manually in
        // the XML file or programmatically by reading the model.xml
        // from the COMSOL .mph file.
        RegionMapType volDomMap, surfDomMap;
        BuildRegionMaps(volDomMap, surfDomMap);
        
        StdVector<std::string> names;
        
        std::set<UInt>::iterator it, end;
        std::map<std::string, std::vector<UInt> >::iterator mapit, mapend;
        
        it = volDoms.begin();
        end = volDoms.end();
        UInt idx = 0;
        
        newVolRegions.resize(volDoms.size());
        
        for( ; it != end; it++ ) 
        {
          std::stringstream sstr;
          
          mapit = volDomMap.begin();
          mapend = volDomMap.end();
          
          bool foundDomain = false;
          std::string name;
          for( ; mapit != mapend; mapit++) 
          {
            std::vector<UInt>::iterator foundIt;
            
            foundIt = std::find(mapit->second.begin(), mapit->second.end(), (*it));
            
            if(foundIt != mapit->second.end()) {
              foundDomain = true;
              name = mapit->first;
            }
          }
          
          if(!foundDomain) {
            sstr << "vol_domain_" << (*it);
            name = sstr.str();
          }
          
          StdVector<std::string>::iterator nameIt;
          nameIt = std::find(names.Begin(), names.End(), name);
          
          if(nameIt == names.End()) {
            names.Push_back(name);
          }
          
          nameIt = std::find(names.Begin(), names.End(), name);
          
          newVolRegions[idx++] = std::distance(names.Begin(), nameIt);
        }
        
        it = surfDoms.begin();
        end = surfDoms.end();
        idx = 0;
        
        newSurfRegions.resize(surfDoms.size());
        
        for( ; it != end; it++ ) 
        {
          std::stringstream sstr;
          
          mapit = surfDomMap.begin();
          mapend = surfDomMap.end();
          
          bool foundDomain = false;
          std::string name;
          for( ; mapit != mapend; mapit++) 
          {
            std::vector<UInt>::iterator foundIt;
            
            foundIt = std::find(mapit->second.begin(), mapit->second.end(), (*it));
            
            if(foundIt != mapit->second.end()) {
              foundDomain = true;
              name = mapit->first;
            }
          }
          
          if(!foundDomain) {
            sstr << "surf_domain_" << (*it);
            name = sstr.str();
          }
          
          StdVector<std::string>::iterator nameIt;
          nameIt = std::find(names.Begin(), names.End(), name);
          
          if(nameIt == names.End()) {
            names.Push_back(name);
          }
          
          nameIt = std::find(names.Begin(), names.End(), name);
          
          newSurfRegions[idx++] = std::distance(names.Begin(), nameIt);
        }
        
        data->AddRegions(names, ids);
        maxNumNodesPerElem = maxnodes;
        
        globTopo.Resize(noelements*maxNumNodesPerElem);
        //   AllocateKnots(data);
        allocated = true;
      }
      
      if(!readInfos) 
        goto omstart;    
    }
    fclose(in);
    
    LOG_TRACE(simInputMPHTXT) << "The Comsol mesh was loaded from file " <<
      sstr.str() << std::endl;

    return(0);
  }
  
  
  
  SimInputMPHTXT::SimInputMPHTXT(std::string fileName, PtrParamNode inputNode,
                                 PtrParamNode infoNode ) :
    SimInput(fileName, inputNode, infoNode)
  {
    capabilities_.insert( SimInput::MESH);

    std::string mph = ""; 
    myParam_->GetValue("mph", mph, ParamNode::PASS );
    mph_ = mph;
  }
  
  
  SimInputMPHTXT::~SimInputMPHTXT()
  {
    inFile_.close();
  }
  
  
  void SimInputMPHTXT::InitModule()
  {
    
    inFile_.open( fileName_.c_str(), std::ios::binary );
    if ( !inFile_.good() ) {
      EXCEPTION("I am unable to open COMSOL .mphtxt file '" << fileName_ << "'");
    }
    
    inFile_.seekg(0, std::ios::end);
    pos_end = inFile_.tellg();
    
    LoadComsolMesh(0, fileName_.c_str(), true);
    
  }
  
  void SimInputMPHTXT::ReadMesh(Grid *mi)
  {
    
    mi_ = mi;
    
    LoadComsolMesh(mi, fileName_.c_str(), false);
  }
  
  void SimInputMPHTXT::BuildRegionMaps(RegionMapType& volDomMap,
                                       RegionMapType& surfDomMap)
  {
    typedef boost::tokenizer<boost::char_separator<char> > Tok;
    boost::char_separator<char> sep(",;| ");

    if(mph_ != "") 
    {
      LOG_TRACE(simInputMPHTXT) << "Multiphysics file: " << mph_ << std::endl;

      // Read model.xml from COMSOL .mph multiphysics file.
      unzFile uf = unzOpen64(mph_.c_str());
      if(!uf)
      {
        EXCEPTION("Could not open COMSOL multiphysics file '" << mph_ << "'!");
      }      

      int res = unzLocateFile( uf, "model.xml", 1 ); // case sensitive
      if( res != UNZ_OK ) 
      {
        unzClose(uf);
        EXCEPTION("Could not locate model.xml inside '" << mph_ << "'!");
      }      

      unz_file_info fi;
      unzGetCurrentFileInfo( uf, &fi, NULL, 0, NULL, 0, NULL, 0 );
      if( fi.uncompressed_size == 0 )
      {
        unzClose(uf);
        EXCEPTION("Could not get file info for model.xml!");
      }      
      int outLength = fi.uncompressed_size;
      
      StdVector<char> outData(outLength+1);
      
      res = unzOpenCurrentFile( uf );
      if( res != UNZ_OK )
      {
        unzClose(uf);
        EXCEPTION("Could not open model.xml!");
      }      
      
      int bytesRead = unzReadCurrentFile( uf, &outData[0], outLength );
      if( bytesRead != outLength ) 
      {
        unzCloseCurrentFile( uf );
        unzClose(uf);
        EXCEPTION("Number of bytes read " << bytesRead << " is different " << 
                  "from size " << outLength << " of model.xml!");
      }      
        
      res = unzCloseCurrentFile( uf );
      if( res != UNZ_OK )
      {
        unzClose(uf);
        EXCEPTION("Problem while closing model.xml!");
      }      

      unzClose(uf);
      
      outData[outLength] = 0;
      std::string modelXmlStr = &outData[0];
      
      // Parse model.xml into a ParamNode using our Xerces parser.
      CoupledField::Xerces xerces;
      xerces.SetString(modelXmlStr);
      PtrParamNode modelNode = xerces.CreateParamNodeInstance();

      std::cout << modelNode->GetName() << std::endl;
      PtrParamNode selectionNode = modelNode->Get("selection");
      std::cout << selectionNode->GetName() << std::endl;

      ParamNodeList regions = selectionNode->GetList("selections");
      ParamNodeList::iterator regit = regions.Begin();
      std::string name, domains, domtype;
      for( ; regit != regions.End(); regit++ ) {
        if((*regit)->Has( "entityName" )) 
        {  
          // For COMSOL 4.2
          (*regit)->GetValue( "entityName", name );
          (*regit)->GetValue( "domains", domains );
          (*regit)->GetValue( "hDim", domtype );
        }
        else 
        {
          // For COMSOL 4.4
          PtrParamNode valNode = (*regit)->Get("value");
          valNode->GetValue( "entityName", name );
          valNode->GetValue( "domains", domains );
          valNode->GetValue( "hDim", domtype );
        }        

        LOG_TRACE(simInputMPHTXT) << "name: " << name << " domains: " << domains
                                  << " domType " << domtype << std::endl;

        std::stringstream sstr;
        UInt hDim;
        sstr << domtype;
        sstr >> hDim;
        
        Tok t(domains, sep);
        Tok::iterator it, end;
        it = t.begin(); it++;
        end = t.end();
        
        for( ; it != end; it++) {
          sstr.clear(); sstr.str(*it);
          UInt dom;
          sstr >> dom;
          
          if(hDim == dim_) {
            volDomMap[name].push_back(dom);
          } else {
            surfDomMap[name].push_back(dom-1);
          }
        }
      }      
    }
    else 
    {
      // Aquire infos about physical meaning of geometric entities inside
      // .mphtxt from .xml file.
      ParamNodeList regions = myParam_->GetList("region");
      ParamNodeList::iterator regit = regions.Begin();
      std::string name, domains, domtype;
      
      for( ; regit != regions.End(); regit++ ) {
        
        (*regit)->GetValue( "name", name );
        (*regit)->GetValue( "domains", domains );
        (*regit)->GetValue( "domtype", domtype );
        
        Tok t(domains, sep);
        Tok::iterator it, end;
        it = t.begin();
        end = t.end();
        
        for( ; it != end; it++) {
          std::stringstream sstr;
          sstr << (*it);
          UInt dom;
          sstr >> dom;
          
          if(domtype == "volume") {
            volDomMap[name].push_back(dom);
          } else {
            surfDomMap[name].push_back(dom);
          }
        }
      }
    }
  }

  // ======================================================
  // GENERAL MESH INFORMATION
  // ======================================================
  UInt SimInputMPHTXT::GetDim() {
    return dim_;
  }
  
  UInt SimInputMPHTXT::GetNumNodes(){
    return maxNumNodes_;
  }
    
  UInt SimInputMPHTXT::GetNumElems(const Integer dim){    
    return maxNumElems_;
  }
  
  UInt SimInputMPHTXT::GetNumRegions(){
    return regionNames_.size();
  }

  UInt SimInputMPHTXT::GetNumNamedNodes(){
    return 0;
  }

  UInt SimInputMPHTXT::GetNumNamedElems(){
    return 0;
  }
  
  // ======================================================
  // ENTITY NAME ACCESS
  // ======================================================

  void SimInputMPHTXT::GetAllRegionNames( StdVector<std::string> & regionNames ){
  }
    
  void SimInputMPHTXT::GetRegionNamesOfDim( StdVector<std::string> & regionNames,
                                       const UInt dim ) {
  }

  void SimInputMPHTXT::GetNodeNames( StdVector<std::string> &nodeNames ) {
  }

  void SimInputMPHTXT::GetElemNames( StdVector<std::string> & elemNames ) {
  }

  // ======================================================
  // ENTITY ACCESS
  // ======================================================
  
  void SimInputMPHTXT::GetCoordinates( std::vector< double > & nodeCoords ) {
  }


  void SimInputMPHTXT::GetNodesOfRegions( std::vector<std::vector<UInt> > &nodes,
                                     const std::vector<RegionIdType> & regionId ) {
  }
    
  void SimInputMPHTXT::GetElements( std::vector< std::vector<UInt> > & elems,
                                  std::vector< std::vector<Elem::FEType> > & elemTypes,
                                  std::vector< std::vector<UInt> > & elemNums,                                
                                  std::vector<RegionIdType> & regionIds,
                                  const UInt dim ) {
  }

  void SimInputMPHTXT::GetNamedNodes(StdVector<StdVector<UInt> > &nodes,
                                   StdVector<std::string> &nodeNames )
  {
  }

  void SimInputMPHTXT::GetNamedElems(StdVector<StdVector<UInt> > & elems,
                                   StdVector<std::string> & elemNames)
  {  
  }
  
  // ***************
  //   Type2ptElem
  // ***************
  Elem::FEType SimInputMPHTXT::ElmerType2ElemType( const UInt itype ) {


    switch( itype ) {
    case 101:
      return Elem::ET_POINT;
    case 202:
      return Elem::ET_LINE2;
    case 203:
      return Elem::ET_LINE3;
    case 303:
      return Elem::ET_TRIA3;
    case 306:
      return Elem::ET_TRIA6;
    case 404:
      return Elem::ET_QUAD4;
    case 408:
      return Elem::ET_QUAD8;
    case 409:
      return Elem::ET_QUAD9;
    case 504:
      return Elem::ET_TET4;
    case 510:
      return Elem::ET_TET10;
    case 808:
      return Elem::ET_HEXA8;
    case 820:
      return Elem::ET_HEXA20;
    case 827:
      return Elem::ET_HEXA27;
    case 605:
      return Elem::ET_PYRA5;
    case 613:
      return Elem::ET_PYRA13;
    case 614:
      return Elem::ET_PYRA14;
    case 706:
      return Elem::ET_WEDGE6;
    case 715:
      return Elem::ET_WEDGE15;
    case 718:
      return Elem::ET_WEDGE18;
    default:
      break;
    }
    
    // This place should never be reached!
    return Elem::ET_UNDEF;
  }

}
