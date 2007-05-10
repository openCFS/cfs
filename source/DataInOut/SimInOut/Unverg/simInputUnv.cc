
#include <algorithm>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <set>
#include <map>

#include "simInputUnv.hh"
#include "unv_if.hh"

extern const char *nodeDataTypesStr[30];
extern const char *elemDataTypesStr[10];

namespace CoupledField {

  // declare logging stream
  DEFINE_LOG(simInputUNV, "SimInputUNV")

  void SimInputUnv::InitModule(Grid *mi)
  {
    mi_ = mi;
  }

  // ======================================================
  // GENERAL MESH INFORMATION
  // ======================================================
  UInt SimInputUnv::GetDim() {
    CoupledField::Warning("SimInputUnv::ReadMesh() not implemented");
    return 0;
  }
  
  UInt SimInputUnv::GetNumNodes(){
    CoupledField::Warning("SimInputUnv::ReadMesh() not implemented");
    return 0;
  }
    
  UInt SimInputUnv::GetNumElems(const Integer dim){
    CoupledField::Warning("SimInputUnv::ReadMesh() not implemented");
    return 0;
  }
  
  UInt SimInputUnv::GetNumRegions(){
    CoupledField::Warning("SimInputUnv::ReadMesh() not implemented");
    return 0;
  }

  UInt SimInputUnv::GetNumNamedNodes(){
    CoupledField::Warning("SimInputUnv::ReadMesh() not implemented");
    return 0;
  }

  UInt SimInputUnv::GetNumNamedElems(){
    CoupledField::Warning("SimInputUnv::ReadMesh() not implemented");
    return 0;
  }
  
  // ======================================================
  // ENTITY NAME ACCESS
  // ======================================================

  void SimInputUnv::GetAllRegionNames( std::vector<std::string> & regionNames ){
    CoupledField::Warning("SimInputUnv::ReadMesh() not implemented");
  }

  void SimInputUnv::GetRegionNamesOfDim( StdVector<std::string> & regionNames,
                                   const UInt dim )
  {
    EXCEPTION("SimInputUnv::GetRegionNamesOfDim() not implemented");
  }
    
  void SimInputUnv::GetNodeNames( StdVector<std::string> & nodeNames )
  {
    EXCEPTION("SimInputUnv::GetNodeNames() not implemented");
  }
  
  void SimInputUnv::GetElemNames( StdVector<std::string> & elemNames )
  {
    EXCEPTION("SimInputUnv::GetElemNames() not implemented");
  }

  void
  SimInputUnv::ReadMesh()
  {
    //!! TODO also extract surface meshes

    LOG_TRACE(simInputUNV) << "reading base mesh from file " << fileName_ << ":";
    
    long n,set;
    CapaInterfaceC capaInterface;
    GDataInfo      datainfo;
    

    if(capaInterface.ReadUniversalfile(fileName_.c_str()) < 0)
    {
      EXCEPTION("Could not read universal file " << fileName_);
      return;
    }
    
    LOG_TRACE(simInputUNV) << "Finished reading";
    
    int nnumtimesteps=0;
    int nnodes=capaInterface.GetNumOfNodes();
    int nelems=capaInterface.GetNumOfElements();
    int dim=capaInterface.GetDimension();
    int nset55=capaInterface.GetNumOfNodeDataSets();
    int nset56=capaInterface.GetNumOfElemDataSets();

    capaInterface.GetDataInfo(datainfo);
 #if 0
    for (set=0; set<nset55; set++) {
      std::cout << "*** node data set: " << set << std::endl;
      std::cout << " datainfo.types55[set] " << datainfo.types55[set] << " => " << nodeDataTypesStr[datainfo.types55[set]] << std::endl;
      int setid = datainfo.Nsetinfo[0][set].idx;
      int edata = datainfo.Nsetinfo[0][set].ndata;
      std::cout << " idx   = " << setid <<std::endl;
      std::cout << " ndata = " << edata <<std::endl;

    }
    for (set=0; set<nset56; set++) {
      std::cout << "*** element data set: " << set << std::endl;
      std::cout << " datainfo.types56[set] " << datainfo.types56[set] << " => " << elemDataTypesStr[datainfo.types56[set]] << std::endl;
      int setid = datainfo.Esetinfo[0][set].idx;
      int edata = datainfo.Esetinfo[0][set].ndata;
      std::cout << " idx   = " << setid <<std::endl;
      std::cout << " ndata = " << edata <<std::endl;
    }
 #endif

    LOG_TRACE(simInputUNV) << "Number of nodes            : "<<nnodes;
    LOG_TRACE(simInputUNV) << "Number of elements         : "<<nelems;
    LOG_TRACE(simInputUNV) << "Number of node data sets   : "<<nset55;
    LOG_TRACE(simInputUNV) << "Number of element data sets: "<<nset56;
    

    // get coordinates
    LOG_TRACE(simInputUNV) << "reading vertex coordinates";

    Point  p;
    long   node,elem;
    bool   isDimPlausible = false;
    //    ocs::Vec3d oldNormal;

    // Check if dim is really sane!
    for (node=1; node<=nnodes; node++) {
      capaInterface.GetPos(node, &p[0]);

      if(p[0] != 0)
      {
        dim = 3;
        break;
      }
    }

    mi_->AddNodes(nnodes);

    for (node=1; node<=nnodes; node++) {
      capaInterface.GetPos(node, &p[0]);

      // In the 2D case CFS++ writes the 3D coordinates in the
      // following form: 0 x y
      if(dim == 2)
      {
        p[0] = p[1];
        p[1] = p[2];
        p[2] = 0;
      }

      mi_->SetNodeCoordinate(node, p);
    }

    

    // get elements
    long elemColor, unvElemType, numElemNodes, elemNodes[32], dummy;
    std::vector< FEType >  element_types(nelems);
    std::vector< UInt > num_vertices_of_elements(nelems);
    std::vector< UInt > element_partition(nelems);
    std::set<UInt> partitions;
    FEType eType;
    std::stringstream strBuffer;


    mi_->AddElems(nelems);

    long totalidxs = 0;
    for (n=0; n<nelems; ++n) {
      capaInterface.GetElemColor(n+1, elemColor);
      capaInterface.GetElemType(n+1, unvElemType);
      capaInterface.GetElemNodes(n+1, numElemNodes, elemNodes);

      eType = UnvType2ElemType(unvElemType);

      if(dim == 2)
      {
        switch(eType)
        {
        case ET_TRIA6:
          dummy = elemNodes[1];
          elemNodes[1] = elemNodes[2];
          elemNodes[2] = elemNodes[4];
          elemNodes[4] = elemNodes[3];
          elemNodes[3] = dummy;
          break;
        case ET_QUAD8:
          dummy = elemNodes[1];
          elemNodes[1] = elemNodes[2];
          elemNodes[2] = elemNodes[4];
          elemNodes[4] = elemNodes[3];
          elemNodes[3] = elemNodes[6];
          elemNodes[6] = elemNodes[5];
          elemNodes[5] = elemNodes[4];
          elemNodes[4] = dummy;
          break;
        }
      }        

      mi_->SetElemData(elem, eType, elemColor-1,
                       (UInt*)elemNodes);

      partitions.insert(elemColor);
    }

    std::vector<std::string> regionNames;
    std::vector<UInt> regionDims;
    std::vector<RegionIdType> regionIds;
    std::set<UInt>::iterator pit, pend;

    pit = partitions.begin();
    pend = partitions.end();

    for(; pit != pend ; pit++)
    {
      strBuffer.str("");
      strBuffer.clear();
      strBuffer << "Partition" << (*pit);
      regionNames.push_back(strBuffer.str());
      regionDims.push_back(3);
    }
    mi_->AddRegions(regionNames, regionIds);

//    for (n=0; n<nelems; ++n) {
//      mi_->SetElemRegion(eit, regionIds[mi_->GetElemRegion(eit)]);
//    }
    /*
      std::vector< UInt > element_vertex_ids(totalidxs);
      std::vector< UInt >::iterator idxptr = element_vertex_ids.begin();
      long *clbuf;
      int i;
      for (n=1; n<=nelems; ++n ) {
      capaInterface.GetElemNodes(n, numElemNodes, elemNodes);

      if(noQuadElems)
      {
      if(dim == 3) 
      {
      switch (numElemNodes) {
      case 10:
      numElemNodes = 4;
      break; // quad. tetrahedron
      case 13:
      numElemNodes = 5;
      break; // quad. prism
      case 15:
      numElemNodes = 6;
      break; // quad. pyramid
      case 20:
      numElemNodes = 8;
      break; // quad. hexahedron
      }
      }
      else // dim == 2
      {
      switch (numElemNodes) {
      case 6:
      numElemNodes = 3;
      quadElem = true;
      break; //quad. triangle
      case 8:
      numElemNodes = 4;
      quadElem = true;
      break; //quad. quad
      }
      }

      // for 2D meshes UNV defines a different ordering of quadratic nodes
      if((dim == 2) && quadElem)
      {
      for(i=0,clbuf=elemNodes;i<numElemNodes;i++)
      {
      *idxptr++ = *(clbuf + (i*2));
      }
      }        
      else
      {
                
      for(i=0,clbuf=elemNodes;i<numElemNodes;i++)
      {
      *idxptr++ = *clbuf++;
      }
      }
            
      }
      else
      {
      // for 2D meshes UNV defines a different ordering of quadratic nodes
      if((dim == 2) && ((element_types[n-1] == ocs::OC_QUAD_TRIANGLE) ||
      (element_types[n-1] == ocs::OC_QUAD_QUADRILATERAL)))
      {
      int offset = numElemNodes/2;
            
      for(i=0,clbuf=elemNodes;i<numElemNodes;i++)
      {
      if(i < offset)
      {
      *idxptr++ = *(clbuf + (i*2));
      }
      else
      {
      *idxptr++ = *(clbuf + ((i-offset)*2+1));
      }
      }
      }        
      else
      {     
      for(i=0,clbuf=elemNodes;i<numElemNodes;i++)
      {
      *idxptr++ = *clbuf++;
      }
      }
      }
        
            
      }


    

      // --------------------------------------------------------------------------------
      // ORCAN SPECIFIC
      // --------------------------------------------------------------------------------

      // setup mesh

      ocs::VolMesh::Interfaces::CreatePtrType         mc;
      ocs::VolMesh::Interfaces::ElementAttribPtrType  ma;
      ocs::VolMesh::Interfaces::AttribInfoPtrType     mi;
      ocs::VolMesh::Interfaces::VertexAttribPtrType   mv;


      if( !(mc=mMesh.I.CreatePtr) ) {
      MESHIO_ERROR( "could not get IVolMeshCreate interface" );
      return false;
      }
      if( !(ma=mMesh.I.ElementAttribPtr) ) {
      MESHIO_ERROR( "could not get IVolMeshElementAttrib interface" );
      return false;
      }
      if( !(mi=mMesh.I.AttribInfoPtr) ) {
      MESHIO_ERROR( "could not get IVolMeshAttribInfo interface" );
      return false;
      }
      if( !(mv=mMesh.I.VertexAttribPtr) ) {
      MESHIO_ERROR( "could not get IVolMeshVertexAttrib interface" );
      return false;
      }

      // create vertex-array
      std::vector<UInt> vtx_ids;
      mc->CreateVertices(vv.size(),vtx_ids);

      // change vertex-IDs (convert FORTRAN -> C numbering)
      std::vector<UInt>::iterator it,eit;
      for( it=element_vertex_ids.begin(),eit=element_vertex_ids.end(); it!=eit; ++it) {
      *it = vtx_ids[*it-1]; 
      }
    
      // add coordinates
      MESHIO_DEBUG( ".unv adding coordinates" );
    
      UInt attribid;
      std::vector< UInt >::iterator it0;
      std::vector< UInt >::const_iterator it2,eit2;
      std::vector< ocs::Vec3d >::iterator it3;
      real64 min, max, length;
      std::vector<real64> coord_maximums;
      std::vector<real64> coord_minimums;

      if( mv->CreateVertexAttrib("Coordinate",3,"real64",ocs::Vec3d(0.,0.,0.),attribid) ) {        
      it3 = vv.begin();
      length = sqrt((*it3)[0]*(*it3)[0] + (*it3)[1]*(*it3)[1] + (*it3)[2] * (*it3)[2]);
      min = length;
      max = length;
      for( it2=vtx_ids.begin(), eit2=vtx_ids.end(); it2!=eit2; ++it2 ) {
      mv->SetVertexAttrib( attribid, *it2, *it3++ );
      length = sqrt((*it3)[0]*(*it3)[0] + (*it3)[1]*(*it3)[1] + (*it3)[2] * (*it3)[2]);
      min = length < min ? length : min;
      max = length > max ? length : max;
      }
      }
      coord_maximums.push_back(max);
      coord_minimums.push_back(min);

      // create elment-array
      MESHIO_DEBUG( ".unv creating element array" );
      MESHIO_DEBUG( ".unv adding " << num_vertices_of_elements.size() << " elements" );
      std::vector< UInt > element_ids;
      mc->CreateElements( num_vertices_of_elements, 
      element_vertex_ids,
      element_types,
      element_ids );

      std::set<UInt> regionset;

      // create element partition attrib
      MESHIO_DEBUG( ".unv creating element-partition attribute" );

      if( ma->CreateElementAttrib("Partition",1,"UInt",(UInt)0,attribid) ) {
        
      std::vector< UInt >::iterator it4,eit4;
      it0 = element_ids.begin();
      for( it4=element_partition.begin(), eit4=element_partition.end(); it4!=eit4; ++it4 ) {
      ma->SetElementAttrib(attribid,*it0++,*it4);
      regionset.insert(*it4);
      }
      }


      // --------------------------------------------------------------------------------
      // UNV SPECIFIC
      // --------------------------------------------------------------------------------

      std::string attrname;
      UInt k;
      std::vector<UInt> stepnrs;
      std::vector<real64> timesteps;
      std::vector<std::string> gridlabels;
      std::vector<UInt> regions;
      std::vector<std::string> attributes;

      for(std::set<UInt>::iterator regit = regionset.begin(); regit != regionset.end(); regit++)
      {
      regions.push_back(*regit);
      }
    

      // get TimeRegionMetaInfo Interface
      LSETimeRegionMeshInterface* trmi=NULL;

      try { 
      if( mMesh.GetProperties().HasProperty("LSETimeRegionMeshInterface") ) {
      oc::CPtr<LSETimeRegionMeshInterface> cptr = mMesh.GetProperties()["LSETimeRegionMeshInterface"](FIELD_ID(oc::CPtr<LSETimeRegionMeshInterface>));
      trmi = cptr.GetCPtr();
      }
      }
      catch(...) {
      CoupledField::Warning("Could not get LSETimeRegionMeshInterface");
      return false;
      }

      trmi->SetRegionAttrName("Partition");
      trmi->SetDisplacementAttrName("");
      attributes.push_back("Partition");
      attributes.push_back("Coordinate");
      trmi->SetAttrList(attributes);
      trmi->SetRegions(regions);

      // get the max/min of the region
      std::vector<UInt>::iterator region, regionend;
      std::vector<real64> reg_maximums;
      std::vector<real64> reg_minimums;
      region = regions.begin();
      regionend = regions.end();
      max = *region;
      min = *region;
      region++;

      for( ; region != regionend; region++)
      {
      max = *region > max ? *region : max;
      min = *region < min ? *region : min;
      }
      reg_maximums.push_back(max);
      reg_minimums.push_back(min);

      for(k=0; k<datainfo.numtype55; k++)
      {
      attrname = nodeDataTypesStr[datainfo.types55[k]];
      attributes.push_back(attrname);
      }

      for(k=0; k<datainfo.numtype56; k++)
      {
      attrname = elemDataTypesStr[datainfo.types56[k]];
      attributes.push_back(attrname);
      MESHIO_INFO("attribute " << k << ": " << attrname);
      MESHIO_INFO("attribute " << 1 << ": " << elemDataTypesStr[datainfo.types56[1]]);
      }

      trmi->SetAttrList(attributes);
    
      double data1[10],data2[10]; 
      std::string name;
      double tensor[6];

      std::map< std::string, std::vector<double> >::iterator find_it;

      // get nodal results and set up min/max etc.
      std::map< std::string, std::vector<double> >     minimums55;
      std::map< std::string, std::vector<double> >     maximums55;

      for(k=0; k<datainfo.numtype55; k++)
      {
      attrname = nodeDataTypesStr[datainfo.types55[k]];
      trmi->SetAttrRegions(attrname, regions);
      if(attrname == "Displacements")
      trmi->SetDisplacementAttrName("Displacements");

      stepnrs.clear();
      timesteps.clear();
      gridlabels.clear();

      std::vector<double>     dv_scalars;
      std::vector<ocs::Vec3d> dv_vectors;
      std::vector<double>     dv_tensors;

      for (set=0; set<datainfo.n55[k]; set++)
      {
      stepnrs.push_back(datainfo.Nsetinfo[k][set].idx);
      timesteps.push_back(datainfo.Nsetinfo[k][set].time);
            
      // generate string for gridlabel
      name = attrname;
      char            buf[128];
      sprintf(buf, "%d", datainfo.Nsetinfo[k][set].idx);
      name += "_step_";
      name += buf;
      name += "_time_";
      sprintf(buf, "%f", datainfo.Nsetinfo[k][set].time);
      name += buf;
            
      gridlabels.push_back(name);

      find_it = minimums55.find(attrname);
      if(find_it == minimums55.end())
      {
      minimums55[attrname] = std::vector<double>(0);
      }

      find_it = maximums55.find(attrname);
      if(find_it == maximums55.end())
      {
      maximums55[attrname] = std::vector<double>(0);
      }


      long ndata=datainfo.Nsetinfo[k][set].ndata;
      switch( ndata ) {
      case 1:
      {
      dv_scalars = std::vector<double>(nnodes);
      for (node=1; node<=nnodes; node++) {
      capaInterface.GetNodeData(node, datainfo.Nsetinfo[k][set].idx, data1, data2);
      dv_scalars[node-1] = data1[0];
      if(node == 1)
      {
      max = data1[0];
      min = data1[0];
      }
      else
      {
      max = data1[0] > max ? data1[0] : max;
      min = data1[0] < min ? data1[0] : min;
      }
      }
      minimums55[attrname].push_back(min);
      maximums55[attrname].push_back(max);
                    
      // add scalar vertex attributes
      if( mv->CreateVertexAttrib(name,1,"real64",real64(0.),attribid) ) {
      std::vector< double >::iterator it3,eit3;
      it0 = vtx_ids.begin();
      for( it3=dv_scalars.begin(), eit3=dv_scalars.end(); it3!=eit3; ++it3 ) {
      mv->SetVertexAttrib(attribid,*it0++,real64(*it3));
      }
      }
      }
      break;
      case 3:
      {
      dv_vectors = std::vector<ocs::Vec3d>(nnodes);
      for (node=1; node<=nnodes; node++) {
      capaInterface.GetNodeData(node, datainfo.Nsetinfo[k][set].idx, data1, data2);

      if(dim == 2)
      {
      dv_vectors[node-1] = ocs::Vec3d(data1[1],data1[2], 0);
      }
      else
      {
      dv_vectors[node-1] = ocs::Vec3d(data1[0],data1[1],data1[2]);
      }
                        
      length = sqrt(data1[0]*data1[0]+
      data1[1]*data1[1]+
      data1[2]*data1[2]);

      if(node == 1)
      {
      max = length;
      min = length;
      }
      else
      {
      max = length > max ? length : max;
      min = length < min ? length : min;
      }
      }
      minimums55[attrname].push_back(min);
      maximums55[attrname].push_back(max);

      // add vector vertex attributes
      if( mv->CreateVertexAttrib(name,3,"real64",ocs::Vec3d(0.,0.,0.),attribid) ) {
                        
      std::vector< ocs::Vec3d >::iterator it3,eit3;
      it0 = vtx_ids.begin();
      for( it3=dv_vectors.begin(), eit3=dv_vectors.end(); it3!=eit3; ++it3 ) {
      mv->SetVertexAttrib(attribid,*it0++,*it3);
      }
      }

      }
      break;
      case 6:
      {
      dv_tensors = std::vector<double>(nnodes*6);
      for (node=1; node<=nnodes; node++) {
      capaInterface.GetNodeData(node, datainfo.Nsetinfo[k][set].idx, data1, data2);
      dv_tensors[(node-1)*6+0] = data1[0];
      dv_tensors[(node-1)*6+1] = data1[1];
      dv_tensors[(node-1)*6+2] = data1[2];
      dv_tensors[(node-1)*6+3] = data1[3];
      dv_tensors[(node-1)*6+4] = data1[4];
      dv_tensors[(node-1)*6+5] = data1[5];
      length = sqrt(data1[0]*data1[0]+
      data1[1]*data1[1]+
      data1[2]*data1[2]+
      data1[3]*data1[3]+
      data1[4]*data1[4]+
      data1[5]*data1[5]);
      // I don't know whether this is the right tensor norm???
      if(node == 1)
      {
      max = length;
      min = length;
      }
      else
      {
      max = length > max ? length : max;
      min = length < min ? length : min;
      }
      }
      minimums55[attrname].push_back(min);
      maximums55[attrname].push_back(max);

      // add scalar vertex attributes
                    
      for(int i = 0; i<6; i++)
      {
      tensor[i] = 0;
      }
                    
      if( mv->CreateVertexAttrib(name,6,"real64",tensor,sizeof(tensor),attribid) ) {
                        
      it0 = vtx_ids.begin();
      for( int i=0; i<dv_tensors.size(); i+=6 ) {
      tensor[0] = dv_tensors[i+0];
      tensor[1] = dv_tensors[i+1];
      tensor[2] = dv_tensors[i+2];
      tensor[3] = dv_tensors[i+3];
      tensor[4] = dv_tensors[i+4];
      tensor[5] = dv_tensors[i+5];
                            
      mv->SetVertexAttrib(attribid,*it0++,tensor, sizeof(tensor));
      }
      }
      }
      break;
      default:
      CoupledField::Warning("Element " << node << " has unusable data");
      break;
      }
            
            
      }

      regionend = regions.end();
      for(region = regions.begin(); region != regionend; region++)
      {
      trmi->SetAttrTimeStepInfos(attrname, *region, stepnrs, timesteps, gridlabels, minimums55[attrname], maximums55[attrname]);
      }
      }


      // get results per element and set up min/max etc.    
      std::map< std::string, std::vector<double> >     minimums56;
      std::map< std::string, std::vector<double> >     maximums56;

      for(k=0; k<datainfo.numtype56; k++)
      {
      attrname = elemDataTypesStr[datainfo.types56[k]];
      trmi->SetAttrRegions(attrname, regions);

      stepnrs.clear();
      timesteps.clear();
      gridlabels.clear();

      std::vector<double>     de_scalars;
      std::vector<ocs::Vec3d> de_vectors;
      std::vector<double>     de_tensors;

      for (set=0; set<datainfo.e56[k]; set++)
      {
      stepnrs.push_back(datainfo.Esetinfo[k][set].idx);
      timesteps.push_back(datainfo.Esetinfo[k][set].time);
            
      name = attrname;
      char            buf[128];
      sprintf(buf, "%d", datainfo.Esetinfo[k][set].idx);
      name += "_step_";
      name += buf;
      name += "_time_";
      sprintf(buf, "%f", datainfo.Esetinfo[k][set].time);
      name += buf;
            
      gridlabels.push_back(name);

      find_it = minimums56.find(attrname);
      if(find_it == minimums56.end())
      {
      minimums56[attrname] = std::vector<double>(0);
      }

      find_it = maximums56.find(attrname);
      if(find_it == maximums56.end())
      {
      maximums56[attrname] = std::vector<double>(0);
      }

      long edata=datainfo.Esetinfo[k][set].ndata;
      switch( edata ) {
      case 1:
      {
      de_scalars = std::vector<double>(nelems);
      for (elem=1; elem<=nelems; elem++) {
      capaInterface.GetElemData(elem, datainfo.Esetinfo[k][set].idx, data1, data2);
      de_scalars[elem-1] = data1[0];
      if(elem == 1)
      {
      max = data1[0];
      min = data1[0];
      }
      else
      {
      max = data1[0] > max ? data1[0] : max;
      min = data1[0] < min ? data1[0] : min;
      }
      }
      minimums56[attrname].push_back(min);
      maximums56[attrname].push_back(max);

      if( ma->CreateElementAttrib(name,1,"real64",real64(0.),attribid) )
      {
      std::vector< double >::iterator it3,eit3;
      it0 = element_ids.begin();
      for( it3=de_scalars.begin(), eit3=de_scalars.end(); it3!=eit3; ++it3 ) {
      ma->SetElementAttrib(attribid,*it0++,(real64)*it3);
      }
      }
      }
      break;
      case 3:
      {
      de_vectors = std::vector<ocs::Vec3d>(nelems);
      for (elem=1; elem<=nelems; elem++) {
      capaInterface.GetElemData(elem, datainfo.Esetinfo[k][set].idx, data1, data2);

      if(dim == 2)
      {
      de_vectors[elem-1] = ocs::Vec3d(data1[1],data1[2], 0);
      }
      else
      {
      de_vectors[elem-1] = ocs::Vec3d(data1[0],data1[1],data1[2]);
      }

      length = sqrt(data1[0]*data1[0]+
      data1[1]*data1[1]+
      data1[1]*data1[1]);

      if(elem == 1)
      {
      max = length;
      min = length;
      }
      else
      {
      max = length > max ? length : max;
      min = length < min ? length : min;
      }
      }
      minimums56[attrname].push_back(min);
      maximums56[attrname].push_back(max);

      if( ma->CreateElementAttrib(name,3,"real64",ocs::Vec3d(0.,0.,0.),attribid) )
      {      
      std::vector< ocs::Vec3d >::iterator it3,eit3;
      it0 = element_ids.begin();
      for( it3=de_vectors.begin(), eit3=de_vectors.end(); it3!=eit3; ++it3 ) {
      ma->SetElementAttrib(attribid,*it0++,*it3);
      }
      }

      }
      break;
      case 6:
      {
      de_tensors = std::vector<double>(nelems*6);
      for (elem=1; elem<=nelems; elem++) {
      capaInterface.GetElemData(elem, datainfo.Esetinfo[k][set].idx, data1, data2);
      de_tensors[(elem-1)*6+0] = data1[0];
      de_tensors[(elem-1)*6+1] = data1[1];
      de_tensors[(elem-1)*6+2] = data1[2];
      de_tensors[(elem-1)*6+3] = data1[3];
      de_tensors[(elem-1)*6+4] = data1[4];
      de_tensors[(elem-1)*6+5] = data1[5];
      length = sqrt(data1[0]*data1[0]+
      data1[1]*data1[1]+
      data1[2]*data1[2]+
      data1[3]*data1[3]+
      data1[4]*data1[4]+
      data1[5]*data1[5]);
      // I don't know wether this is the right tensor norm???
      if(elem == 1)
      {
      max = length;
      min = length;
      }
      else
      {
      max = length > max ? length : max;
      min = length < min ? length : min;
      }
      }
      minimums56[attrname].push_back(min);
      maximums56[attrname].push_back(max);

      for(int i = 0; i<6; i++)
      {
      tensor[i] = 0;
      }

      std::cout << name << std::endl;
      if( ma->CreateElementAttrib(name,6,"real64",tensor,sizeof(tensor),attribid) ) 
      {
      it0 = element_ids.begin();
      int i=0;
                        
      for( ; it0 != element_ids.end(); it0++ ) {
      tensor[0] = de_tensors[i+0];
      tensor[1] = de_tensors[i+1];
      tensor[2] = de_tensors[i+2];
      tensor[3] = de_tensors[i+3];
      tensor[4] = de_tensors[i+4];
      tensor[5] = de_tensors[i+5];
                            
      ma->SetElementAttrib(attribid,*it0,tensor, sizeof(tensor));
      i+=6;
      }
      }
      }
      break;
      default:
      CoupledField::Warning("Element "<<elem <<" has unusable data");
      break;
      }
            
      }
        
      regionend = regions.end();
      for(region = regions.begin(); region != regionend; region++)
      {
      trmi->SetAttrTimeStepInfos(attrname, *region, stepnrs, timesteps, gridlabels, minimums56[attrname], maximums56[attrname]);
      }
      }


      // add Partition and Coordinate attribute to mesh
      stepnrs.clear();
      timesteps.clear();
      gridlabels.clear();
        
      stepnrs.push_back(0);
      timesteps.push_back(0);

      trmi->SetAttrRegions("Partition", regions);
      gridlabels.push_back("Partition");

      regionend = regions.end();
      for(region = regions.begin(); region != regionend; region++)
      {
      trmi->SetAttrTimeStepInfos("Partition", *region, stepnrs, timesteps, gridlabels, reg_minimums, reg_maximums);
      }


      trmi->SetAttrRegions("Coordinate", regions);
      gridlabels.clear();
      gridlabels.push_back("Coordinate");

      for(region = regions.begin(); region != regionend; region++)
      {
      trmi->SetAttrTimeStepInfos("Coordinate", *region, stepnrs, timesteps, gridlabels, coord_minimums, coord_maximums);
      }
    */  
    return;
  }

  FEType SimInputUnv::UnvType2ElemType( const UInt elemType ) {

    switch (elemType) {
    case 91:  return ET_TRIA3;
    case 94:  return ET_QUAD4;
    case 92:  return ET_TRIA6;
    case 95:  return ET_QUAD8;
    case 111: return ET_TET4;  // tetraeder 1.ord
    case 112: return ET_WEDGE6;  // prism     1.ord
    case 115: return ET_HEXA8;  // hexaeder  1.ord
    case 113: return ET_WEDGE15;  // prism     2.ord
    case 116: return ET_HEXA20;  // hexaeder  2.ord
    }

    // This place should never be reached!
    return ET_UNDEF;
  }

  bool SimInputUnv::ReadData()
  {
  }
  
}