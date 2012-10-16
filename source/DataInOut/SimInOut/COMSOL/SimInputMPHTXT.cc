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

#include "SimInputMPHTXT.hh"

namespace CoupledField {

#define MAXLINESIZE 200     /* maximum length of line to be read */
#define MAXNODESD3 64       /* maximum number of 3D nodes */ 
#define MAXNODESD2 27       /* maximum number of 2D nodes */ 
#define MAXNODESD1 9        /* maximum number of 1D nodes */
#define MAXFILESIZE 600     /* maximum filenamesize for i/o files */

int next_int(char **start)
{
  int i;
  char *end;

  i = strtol(*start,&end,10);
  *start = end;
  return(i);
}

Double next_real(char **start)
{
  Double r;
  char *end;

  r = strtod(*start,&end);

  *start = end;
  return(r);
}

static int Comsolrow(char *line1,FILE *io) 
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

static void ReorderComsolNodes(int elementtype,UInt *topo) 
{
  int i,tmptopo[MAXNODESD2];
  int order306[]={1,2,3,4,6,5};
  int order404[]={1,2,4,3};
  int order409[]={1,2,4,3,5,8,9,6,7};
  int order510[]={1,2,3,4,5,7,6,8,9,10};
  int order718[]={1,2,3,4,5,6,7,9,8,16,18,17,10,12,15,11,14,13};
  int order808[]={1,2,4,3,5,6,8,7};
 
  for(i=0;i<elementtype%100;i++) 
    tmptopo[i] = topo[i];
   
  switch (elementtype) {
 
  case 306:
    for(i=0;i<elementtype%100;i++) 
      topo[i] = tmptopo[order306[i]-1];
    break;

  case 404:
    for(i=0;i<elementtype%100;i++) 
      topo[i] = tmptopo[order404[i]-1];
    break;

  case 409:
    for(i=0;i<elementtype%100;i++) 
      topo[i] = tmptopo[order409[i]-1];
    break;

  case 510:
    for(i=0;i<elementtype%100;i++) 
      topo[i] = tmptopo[order510[i]-1];
    break;
    
  case 808:
    for(i=0;i<elementtype%100;i++) 
      topo[i] = tmptopo[order808[i]-1];
    break;

  case 718:
    for(i=0;i<elementtype%100;i++) 
      topo[i] = tmptopo[order718[i]-1];
    break;

  case 408:
  case 820:
  case 827:
    // case 510:
  case 715:
    // case 718:
    EXCEPTION("Reordering for Elmer elemtype " << elementtype << " not yet implemented.");    
    break;

  default:
    break;
  }

}

  int LoadComsolMesh(Grid *data, const char *prefix, int info, bool readInfos, UInt& dim, UInt& numNodes, UInt& numElems, PtrParamNode myParam)
/* Load the grid in Comsol Multiphysics mesh format */
{
  int noknots,noelements,noelements_dom,elemcode,maxnodes,material;
  // int foundsame;
  // int mode,nvalue,maxknot,nosides,sideelemtype;
  int allocated;
  //  int boundarytype,materialtype,boundarynodes,side,parent,elemsides;
  int elemnodes, elembasis, elemtype;
  elembasis = 0;
  elemnodes = 0;
  dim = 0;
  // int bulkdone, usedmax,hits;
  // int label;
  int debug,offset,domains,mindom,minbc,elemdim = 0;
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
  UInt maxNumNodes = 0;

  std::stringstream sstr;
  sstr << prefix; //  << ".mphtxt";
  if ((in = fopen(sstr.str().c_str(),"r")) == NULL) {
      printf("LoadComsolMesh: opening of the Comsol mesh file '%s' wasn't succesfull !\n",
	     sstr.str().c_str());
      return(1);
  }

  printf("Reading mesh from Comsol mesh file %s.\n",sstr.str().c_str());
  // InitializeKnots(data);

  debug = true;
  allocated = false;

  mindom = 1000;
  minbc = 1000;
  offset = 1;

omstart:

  maxnodes = 0;
  noknots = 0;
  noelements = 0;
  noelements_dom = 0;
  elemcode = 0;
  material = 0;
  domains = 0;


  for(;;) {

    if(Comsolrow(line,in)) goto end;
    if(line == NULL) goto end;

    if(strstr(line,"# sdim")) {
      cp = line;
      dim = next_int(&cp);
      if(debug) printf("dim=%d\n",dim);
    }

    else if(strstr(line,"# number of mesh points")) {
      cp = line;
      noknots = next_int(&cp);
      if(debug) printf("noknots=%d\n",noknots);
    }

    else if(strstr(line,"# lowest mesh point index")) {
      cp = line;
      offset = 1 - next_int(&cp);
      if(debug) printf("offset=%d\n",offset);
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
      if(debug) printf("elemnodes=%d\n",elemnodes);

      switch(elembasis) 
      {
      case 100:
        fetype = Elem::ET_POINT;
        break;
      case 200:
        switch(elemnodes) 
        {
        case 3:
          fetype = Elem::ET_LINE3;
          break;
        default:
          fetype = Elem::ET_LINE2;
          break;
        }
        break;
      case 300:
        switch(elemnodes) 
        {
        case 6:
          fetype = Elem::ET_TRIA6;
          break;
        default:
          fetype = Elem::ET_TRIA3;
          break;
        }
        break;
      case 400:
        switch(elemnodes) 
        {
        case 8:
          fetype = Elem::ET_QUAD8;
          break;
        case 9:
          fetype = Elem::ET_QUAD9;
          break;
        default:
          fetype = Elem::ET_QUAD4;
          break;
        }
        break;
      case 500:
        switch(elemnodes) 
        {
        case 10:
          fetype = Elem::ET_TET10;
          break;
        default:
          fetype = Elem::ET_TET4;
          break;
        }
        break;
      case 600:
        switch(elemnodes) 
        {
        case 13:
          fetype = Elem::ET_PYRA13;
          break;
        case 14:
          fetype = Elem::ET_PYRA14;
          break;
        default:
          fetype = Elem::ET_PYRA5;
          break;
        }
        break;
      case 700:
        switch(elemnodes) 
        {
        case 15:
          fetype = Elem::ET_WEDGE15;
          break;
        case 18:
          fetype = Elem::ET_WEDGE18;
          break;
        default:
          fetype = Elem::ET_WEDGE6;
          break;
        }
        break;
      case 800:
        switch(elemnodes) 
        {
        case 20:
          fetype = Elem::ET_HEXA20;
          break;
        case 27:
          fetype = Elem::ET_HEXA27;
          break;
        default:
          fetype = Elem::ET_HEXA8;
          break;
        }
        break;
      }      
    }

    else if(strstr(line,"# Mesh point coordinates")) {
      printf("Loading %d coordinates\n",noknots);

      for(i=1;i<=noknots;i++) {
	Comsolrow(line,in);	

	if(allocated) {
	  cp = line;
          //	  data->x[i] = next_real(&cp);
          //	  data->y[i] = next_real(&cp);
          //	  if(dim == 3) data->z[i] = next_real(&cp);

          UInt nodeNum = i;
          Vector<Double> coord(3);
          coord[2] = 0;
          for(UInt j=0; j<dim; j++) 
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

      if(debug) printf("Loading %d elements of type %d\n",k,elemtype);
      domains = noelements;

      for(i=1;i<=k;i++) {
	Comsolrow(line,in);	

	if(dim == 3 && elembasis < 300) continue;
	if(dim == 2 && elembasis < 200) continue;

	noelements = noelements + 1;
	if(allocated) {
          //	  data->elementtypes[noelements] = elemtype;
          //	  data->material[noelements] = 1;
          
          //          UInt topo[128];
	  cp = line;
	  for(j=0;j<elemnodes;j++)
            globTopo[(noelements-1)*maxNumNodes + j] = next_int(&cp) + offset;
            //	    data->topology[noelements][j] = next_int(&cp) + offset;

          if(1) ReorderComsolNodes(elemtype, &globTopo[(noelements-1)*maxNumNodes] );

          /*
          if(fetype == Elem::ET_POINT) 
          {
            topo[1] = topo[0];
            fetype = Elem::ET_LINE2;
          }
          */

          // data->SetElemData(noelements, fetype, 0, topo);

	}
      }
    }

    else if(strstr(line,"# number of domains") || strstr(line,"# number of geometric entity indices")) {

      cp = line;
      k = next_int(&cp);

      Comsolrow(line,in);	            
      if(debug) printf("Loading %d domains for the elements\n",k);

      for(i=1;i<=k;i++) {
	Comsolrow(line,in);	

	if(dim == 3 && elembasis < 300) continue;
	if(dim == 2 && elembasis < 200) continue;

	domains = domains + 1;
	cp = line;
	material = next_int(&cp);

	noelements_dom = noelements_dom + 1;
	if(allocated) {
          // UInt topo[128];
          RegionIdType id;
          
          // data->GetElemData(noelements_dom, fetype, id, topo);

	  if(elemdim < (int) dim) 
          {
            UInt distance = std::distance(surfDoms.begin(), surfDoms.find(material));
            id = ids[ newSurfRegions[distance] ];

	    material = material - minbc + 1;

            data->SetElemData(noelements_dom, fetype, id, &globTopo[(noelements_dom-1)*maxNumNodes]);
          }
	  else
          {
            UInt distance = std::distance(volDoms.begin(), volDoms.find(material));
            id = ids[ newVolRegions[distance] ];

	    material = material - mindom + 1;
          //	  data->material[domains] = material;

            data->SetElemData(noelements_dom, fetype, id, &globTopo[(noelements_dom-1)*maxNumNodes]);
          }
          
	}
	else {
	  if(elemdim < (int) dim) {
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
      if(debug) printf("Unused command:  %s",line);
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
    //    data->dim = dim;
    
    if(info) {
      printf("Allocating for %d knots and %d %d-node elements.\n",
	     noknots,noelements,maxnodes);
    }
    
    if(readInfos) 
    {
      numNodes = noknots;
      numElems = noelements;
    }
    else 
    {
      data->AddElems(noelements);
      data->AddNodes(noknots);

      ParamNodeList regions = myParam->GetList("region");
      ParamNodeList::iterator regit = regions.Begin();
      std::string name, domains, domtype;

      std::map<std::string, std::vector<UInt> > volDomMap;
      std::map<std::string, std::vector<UInt> > surfDomMap;

      for( ; regit != regions.End(); regit++ ) {
 
        (*regit)->GetValue( "name", name );
        (*regit)->GetValue( "domains", domains );
        (*regit)->GetValue( "domtype", domtype );

        typedef boost::tokenizer<boost::char_separator<char> > Tok;
        boost::char_separator<char> sep(";| ");
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
       
/*      
      volDomMap["wasser"].push_back(3);
      volDomMap["wasser"].push_back(4);
*/      
      
/*      
      surfDomMap["surfaces"].push_back(0);
      surfDomMap["surfaces"].push_back(1);
      surfDomMap["surfaces"].push_back(2);
      surfDomMap["surfaces"].push_back(3);
      surfDomMap["surfaces"].push_back(4);
      surfDomMap["surfaces"].push_back(5);
      surfDomMap["surfaces"].push_back(6);
      surfDomMap["surfaces"].push_back(7);
      surfDomMap["surfaces"].push_back(8);
      surfDomMap["surfaces"].push_back(9);
      surfDomMap["surfaces"].push_back(10);
      surfDomMap["surfaces"].push_back(11);
      surfDomMap["surfaces"].push_back(12);
      surfDomMap["surfaces"].push_back(13);
      surfDomMap["surfaces"].push_back(14);
      surfDomMap["surfaces"].push_back(15);
      surfDomMap["surfaces"].push_back(16);
      surfDomMap["surfaces"].push_back(17);
      surfDomMap["surfaces"].push_back(18);
      surfDomMap["surfaces"].push_back(19);
      surfDomMap["surfaces"].push_back(20);
      surfDomMap["surfaces"].push_back(21);
*/
      
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
      maxNumNodes = maxnodes;

      globTopo.Resize(noelements*maxNumNodes);
      //   AllocateKnots(data);
      allocated = true;
    }

    if(!readInfos) 
      goto omstart;    
  }
  fclose(in);

  if(info) printf("The Comsol mesh was loaded from file %s.\n\n",sstr.str().c_str());
  return(0);
}



  SimInputMPHTXT::SimInputMPHTXT(std::string fileName, PtrParamNode inputNode) :
    SimInput(fileName, inputNode)
  {
    capabilities_.insert( SimInput::MESH);
  }


  SimInputMPHTXT::~SimInputMPHTXT()
  {
    inFile_.close();
  }


  void SimInputMPHTXT::InitModule()
  {
    
    inFile_.open( fileName_.c_str(), std::ios::binary );
    if ( !inFile_.good() ) {
      EXCEPTION("I am unable to open Comsol .mphtxt file " << fileName_);
    }

    inFile_.seekg(0, std::ios::end);
    pos_end = inFile_.tellg();

    // dim_ = GetDim();
    // elemDimReadIn_.resize(dim_);
    // maxNumElems_ = GetNumElems();
    // maxNumNodes_ = GetNumNodes();

    LoadComsolMesh(0, fileName_.c_str(), 1, true, dim_, maxNumNodes_, maxNumElems_, myParam_);

    std::cout << "DIM: " << GetDim() << std::endl;
    std::cout << "NUM NODES: " << GetNumNodes() << std::endl;
    std::cout << "NUM ELEMS: " << GetNumElems() << std::endl;

    //    EXCEPTION("HALT");
  }
    
  void SimInputMPHTXT::ReadMesh(Grid *mi)
  {

    mi_ = mi;

    std::cout << "DIM: " << GetDim() << std::endl;
    std::cout << "NUM NODES: " << GetNumNodes() << std::endl;
    std::cout << "NUM ELEMS: " << GetNumElems() << std::endl;
    LoadComsolMesh(mi, fileName_.c_str(), 1, false, dim_, maxNumNodes_, maxNumElems_, myParam_);

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
#if 0
    case 12:
      return Elem::ET_PYRA5;
    case 13:
      //    (*warning) << "Pyram. with quadratic shape functions" << 
      //    "do not work well for some cases "<< 
      //    "(i.e. electric field intensity). Please verify your results."; 
      return Elem::ET_PYRA13;
#endif
    case 706:
      return Elem::ET_WEDGE6;
    case 715:
      return Elem::ET_WEDGE15;
    case 718:
      return Elem::ET_WEDGE15;
    default:
      break;
    }
    
    // This place should never be reached!
    return Elem::ET_UNDEF;
  }

}
