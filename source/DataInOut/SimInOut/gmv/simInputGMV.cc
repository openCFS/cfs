// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

// This code is also based on the gmvread library from the official GMV
// website at: http://www-xdiv.lanl.gov/XCM/gmv/GMVHome.html

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <set>

#include "simInputGMV.hh"

#undef RDATA_INIT
#include "gmvread.hh"

namespace CoupledField {

  std::vector< std::string > SimInputGMV::mPossibleAttribs;
  gmv_data_struct gmv_data;
  gmv_meshdata_struct gmv_meshdata;
  gmvray_data_struct gmvray_data;

  SimInputGMV::SimInputGMV(std::string fileName, ParamNode * inputNode) :
      SimInput(fileName, inputNode)
  {
    if(mPossibleAttribs.empty())
    {
      mPossibleAttribs.push_back("mechAcceleration");
      mPossibleAttribs.push_back("mechVelocity");
      mPossibleAttribs.push_back("mechForce");
      mPossibleAttribs.push_back("mechStress");
      mPossibleAttribs.push_back("mechStrain");
      mPossibleAttribs.push_back("mechEnergy");
      mPossibleAttribs.push_back("mechDisplacement");
      mPossibleAttribs.push_back("elecPotential");
      mPossibleAttribs.push_back("elecFieldIntensity");
      mPossibleAttribs.push_back("elecForceVWP");
      mPossibleAttribs.push_back("elecInterfaceForce");
      mPossibleAttribs.push_back("elecCharge");
      mPossibleAttribs.push_back("elecFluxDensity");
      mPossibleAttribs.push_back("elecEnergy");
      mPossibleAttribs.push_back("smoothDisplacement");
      mPossibleAttribs.push_back("acouPotential");
      mPossibleAttribs.push_back("acouPressure");
      mPossibleAttribs.push_back("acouForce");
      mPossibleAttribs.push_back("acouPotentialD1");
      mPossibleAttribs.push_back("acouPotentialD2");
      mPossibleAttribs.push_back("acouRHSval");
      mPossibleAttribs.push_back("magPotential");
      mPossibleAttribs.push_back("magFluxDensity");
      mPossibleAttribs.push_back("magEddyCurrent");
      mPossibleAttribs.push_back("magForceVWP");
      mPossibleAttribs.push_back("magForceLorentz");
      mPossibleAttribs.push_back("magEnergy");
      mPossibleAttribs.push_back("fluidForce");
      mPossibleAttribs.push_back("elecField");
    }

    mCycleNos.resize(1);
    mProbTimes.resize(1);
  }

  SimInputGMV::~SimInputGMV()
  {
  }

  void SimInputGMV::InitModule()
  {
  }

  // ======================================================
  // GENERAL MESH INFORMATION
  // ======================================================
  UInt SimInputGMV::GetDim() {
    return 3;
  }
  
  void SimInputGMV::GetNumMultiSequenceSteps( std::map<UInt, BasePDE::AnalysisType>& analysis,
                                              std::map<UInt, UInt>& numSteps,
                                              bool isHistory) {
    // At the moment only the grid may be read.
    analysis.clear();
    numSteps.clear();
  }

  UInt SimInputGMV::GetNumNodes(){
    EXCEPTION("SimInputGMV::GetNumNodes() not implemented");
    return 0;
  }
    
  UInt SimInputGMV::GetNumElems(const Integer dim){
    EXCEPTION("SimInputGMV::GetNumElems() not implemented");
    return 0;
  }
  
  UInt SimInputGMV::GetNumRegions(){
    EXCEPTION("SimInputGMV::GetNumRegions() not implemented");
    return 0;
  }

  UInt SimInputGMV::GetNumNamedNodes(){
    EXCEPTION("SimInputGMV::GetNumNamedNodes() not implemented");
    return 0;
  }

  UInt SimInputGMV::GetNumNamedElems(){
    EXCEPTION("SimInputGMV::GetNumNamedElems() not implemented");
    return 0;
  }
  
  // ======================================================
  // ENTITY NAME ACCESS
  // ======================================================

  void SimInputGMV::GetAllRegionNames( StdVector<std::string> & regionNames ){
    EXCEPTION("SimInputGMV::GetAllRegionNames() not implemented");
  }

  void SimInputGMV::GetRegionNamesOfDim( StdVector<std::string> & regionNames,
                                   const UInt dim )
  {
    EXCEPTION("SimInputGMV::GetRegionNamesOfDim() not implemented");
  }
  

  void SimInputGMV::GetNodeNames( StdVector<std::string> & nodeNames )
  {
    EXCEPTION("SimInputGMV::GetNodeNames() not implemented");
  }
  
  void SimInputGMV::GetElemNames( StdVector<std::string> & elemNames )
  {
    EXCEPTION("SimInputGMV::GetElemNames() not implemented");
  }


  void SimInputGMV::ReadMesh(Grid *mi)
  {
    Integer retcode;

    mi_ = mi;
    
    gmvread_printon();
    
    //    gmv_base_dir = mFile.GetParent().GetAbsolutePath();

    try {
      retcode = gmvread_open(fileName_.c_str());
    } catch (GMVReadException& e)
    { 
      EXCEPTION("GMVReadException!")
    }

    if(retcode > 0)
      return;
    
    UInt polys = 0;
    UInt iend = 0;
    while (iend == 0)
    {
      try {
        gmvread_data();
      } catch (GMVReadException& e)
      {
        Cleanup();
        EXCEPTION("GMVReadException!!");
      }

      if(gmv_data.keyword == GMVEND)
        iend = 1;

      if(gmv_data.keyword == GMVERROR)
      {
        Cleanup();
  
        EXCEPTION("The GMV reader encountered an error!");
      }
        
      switch (gmv_data.keyword)
      {
      case(NODES):
        // MESHIO_DEBUG("NODES");

        try {
          gmvread_mesh();
        } catch (GMVReadException& e)
        {
          Cleanup();
          EXCEPTION("GMVReadException!!");
        }

        if(gmv_data.keyword == GMVERROR)
          continue;

        if(!ProcessMesh())
        {
          Cleanup();
          return;
        }
        break;
      case(MATERIAL):
//        MESHIO_DEBUG("MATERIAL");
        if(!ProcessMaterials())
        {
          mMatNames.Clear();
          mMatNames.Push_back("Material");
          mMatNums.Push_back(0);
          mElementMaterials.resize(gmv_meshdata.ncells);
        }
        break;
      case CELLS:
        std::cout << "CELLS" << std::endl;
        break;
      case FACES:
        std::cout << "FACES" << std::endl;
        
        break;
      case VFACES:
        std::cout << "VFACES" << std::endl;
        
        break;
      case XFACES:
        std::cout << "XFACES" << std::endl;
        
        break;
      case VELOCITY:
        std::cout << "VELOCITY" << std::endl;
        
        if(!ProcessVelocities())
        {
          return;
        }
        break;
      case VARIABLE:
        std::cout << "VARIABLE" << std::endl;
        
        if(!ProcessVariables())
        {
          return;
        }
        break;
      case FLAGS:
        std::cout << "FLAGS" << std::endl;
        
        break;
      case POLYGONS:
        //            MESHIO_DEBUG("POLYGONS");
        //            MESHIO_DEBUG(polys);
        polys++;
        break;
      case TRACERS:
        std::cout << "TRACERS" << std::endl;
        
        break;
      case PROBTIME:
        std::cout << "PROBTIME" << std::endl;
        
        if(!ProcessProbtime())
        {
          mProbTimes.resize(1);
          mProbTimes[0]=0;
        }
        break;
      case CYCLENO:
        std::cout << "CYCLENO" << std::endl;
        
        if(!ProcessCycleNo())
        {
          mCycleNos.resize(1);
          mCycleNos[0]=0;
        }
        break;
      case NODEIDS:
        std::cout << "NODEIDS" << std::endl;
        
        break;
      case CELLIDS:
        std::cout << "CELLIDS" << std::endl;
        
        break;
      case SURFACE:
        std::cout << "SURFACE" << std::endl;
        
        break;
      case SURFMATS:
        std::cout << "SURFMATS" << std::endl;
        
        break;
      case SURFVEL:
        std::cout << "SURFVEL" << std::endl;
        
        break;
      case SURFVARS:
        std::cout << "SURFVARS" << std::endl;
        
        break;
      case SURFFLAG:
        std::cout << "SURFFLAG" << std::endl;
        
        break;
      case UNITS:
        std::cout << "UNITS" << std::endl;
        
        break;
      case VINFO:
        std::cout << "VINFO" << std::endl;
        
        break;
      case TRACEIDS:
        std::cout << "TRACEIDS" << std::endl;
        
        break;
      case GROUPS:
        std::cout << "GROUPS" << std::endl;
        
        if(!ProcessGroups())
        {
          return;
        }
        break;
      case FACEIDS:
        std::cout << "FACEIDS" << std::endl;
        
        break;
      case SURFIDS:
        std::cout << "SURFIDS" << std::endl;
        
        break;
      case CELLPES:
        std::cout << "CELLPES" << std::endl;
        break;
      case SUBVARS:
        std::cout << "SUBVARS" << std::endl;
        break;
      case GHOSTS:
        std::cout << "GHOSTS" << std::endl;
        break;
      case CODENAME:
        std::cout << "CODENAME" << std::endl;
        break;
      case CODEVER:
        std::cout << "CODEVER" << std::endl;
        break;
      case SIMDATE:
        std::cout << "SIMDATE" << std::endl;
        break;
      case VECTOR:
        std::cout << "VECTOR" << std::endl;
        if(!ProcessVectors())
        {
          return;
        }            
        break;
      case INVALIDKEYWORD:
        std::cout << "INVALIDKEYWORD" << std::endl;
        break;
      }
    }  
    
    SetupGridAndRegions();

    Cleanup();
  }

  void 
  SimInputGMV::Cleanup()
  {
    gmvread_close();
    
    mMatNames.Clear();
    mMatNums.Clear();
    mElementMaterials.clear();
    mElementIds.clear();
    gmv_cell_types.clear();
    mCycleNos.resize(1);
    mCycleNos[0]=0;
    mProbTimes.resize(1);
    mProbTimes[0]=0;
    mGridLabels.clear();
    mPossibleAttribs.clear();
  }


  bool 
  SimInputGMV::ProcessMesh()
  {
    std::vector< Double > vv;
    long   node,elem;
    UInt iz, iy, ix;
    UInt elemsx, elemsy, elemsz;
    UInt offset_z, offset_z1;
    UInt offset_yz, offset_y1z;
    UInt offset_yz1, offset_y1z1;
    const UInt LSTRUCT = STRUCT + LOGICALLY_STRUCT;
    UInt intype;

    if ((gmv_meshdata.intype == STRUCT) || (gmv_meshdata.intype == LOGICALLY_STRUCT))
      intype = LSTRUCT;
    else
      intype = gmv_meshdata.intype;


    std::vector< UInt > element_vertex_ids(100);

    // --------------------------------------------------------------------------------
    //  GET cells / nodes from GMVs weired cell->face datastructure thingies!
    // --------------------------------------------------------------------------------

    mElementMaterials.resize(gmv_meshdata.ncells);

    UInt numNodes, numNodesOffset;
    numNodesOffset = 0;
    UInt cellNodesPtr[100];

    switch(intype)
    {
    case CELLS:
      mi_->AddNodes(gmv_meshdata.nnodes);

      for (node=0; node<gmv_meshdata.nnodes; node++) {
        Point p;
        
        p[0] = gmv_meshdata.x[node];
        p[1] = gmv_meshdata.y[node];
        p[2] = gmv_meshdata.z[node];
        
        mi_->SetNodeCoordinate(node+1, p);
      }

      mi_->AddElems(gmv_meshdata.ncells);

      for (elem=0; elem<gmv_meshdata.ncells; elem++) {
        //            OCINFO("faces for elem " << elem << ": " << gmv_meshdata.celltoface[elem]);
        // for the meaning of numverts look at regcell() in gmvread.cpp
        //            long *face1, *face2, *face3, *face4, *face5, *face6;            

        numNodes = gmv_meshdata.cellnnode[elem];
        long* pt = gmv_meshdata.cellnodes + numNodesOffset;
        long* ptend = pt + numNodes;
        
        for(UInt i=0; pt!=ptend; pt++, i++)
        {          
          cellNodesPtr[i] = (UInt)(*pt);
        }
        
        if(numNodes == 0)
        {
          mi_->SetElemData(elem+1, Elem::UNDEF, 0, cellNodesPtr);
          continue;
        }
        //        cellNodesPtr = (long*)(gmv_meshdata.cellnodes + numNodesOffset);
            
        switch(gmv_cell_types[elem])
        {

        case LINE:
          mi_->SetElemData(elem+1, Elem::LINE2, 0, cellNodesPtr);
          break;

        case LINE3:
          element_vertex_ids[0] = cellNodesPtr[0];
          element_vertex_ids[1] = cellNodesPtr[2];
          element_vertex_ids[2] = cellNodesPtr[1];
          mi_->SetElemData(elem+1, Elem::LINE3, 0, &element_vertex_ids[0]);
          break;

        case TRI:
          mi_->SetElemData(elem+1, Elem::TRIA3, 0, cellNodesPtr);
          break;

        case TRI6:
          mi_->SetElemData(elem+1, Elem::TRIA6, 0, cellNodesPtr);
          break;

        case QUAD:
          mi_->SetElemData(elem+1, Elem::QUAD4, 0, cellNodesPtr);
          break;

        case QUAD8:               
          mi_->SetElemData(elem+1, Elem::QUAD8, 0, cellNodesPtr);
          break;

        case TET:
          element_vertex_ids[0] = cellNodesPtr[1];
          element_vertex_ids[1] = cellNodesPtr[2];
          element_vertex_ids[2] = cellNodesPtr[3];
          element_vertex_ids[3] = cellNodesPtr[0];
          mi_->SetElemData(elem+1, Elem::TET4, 0, &element_vertex_ids[0]);
          break;

        case PTET4:
          mi_->SetElemData(elem+1, Elem::TET4, 0, cellNodesPtr);
          break;

        case PTET10:
          mi_->SetElemData(elem+1, Elem::TET10, 0, cellNodesPtr);
          break;

        case PYRAMID:
          element_vertex_ids[0] = cellNodesPtr[1];
          element_vertex_ids[1] = cellNodesPtr[2];
          element_vertex_ids[2] = cellNodesPtr[3];
          element_vertex_ids[3] = cellNodesPtr[4];
          element_vertex_ids[4] = cellNodesPtr[0];
          mi_->SetElemData(elem+1, Elem::PYRA5, 0, &element_vertex_ids[0]);
          break;

        case PPYRMD5:                    
          mi_->SetElemData(elem+1, Elem::PYRA5, 0, cellNodesPtr);
          break;

        case PPYRMD13:               
          mi_->SetElemData(elem+1, Elem::PYRA13, 0, cellNodesPtr);
          break;

        case PRISM:
          element_vertex_ids[0] = cellNodesPtr[3];
          element_vertex_ids[1] = cellNodesPtr[4];
          element_vertex_ids[2] = cellNodesPtr[5];
          element_vertex_ids[3] = cellNodesPtr[0];
          element_vertex_ids[4] = cellNodesPtr[1];
          element_vertex_ids[5] = cellNodesPtr[2];
          mi_->SetElemData(elem+1, Elem::WEDGE6, 0, &element_vertex_ids[0]);
          break;

        case PPRISM6:
          mi_->SetElemData(elem+1, Elem::WEDGE6, 0, cellNodesPtr);
          break;

        case PPRISM15:
          mi_->SetElemData(elem+1, Elem::WEDGE15, 0, cellNodesPtr);
          break;

        case HEX:
          element_vertex_ids[0] = cellNodesPtr[4];
          element_vertex_ids[1] = cellNodesPtr[5];
          element_vertex_ids[2] = cellNodesPtr[6];
          element_vertex_ids[3] = cellNodesPtr[7];
          element_vertex_ids[4] = cellNodesPtr[0];
          element_vertex_ids[5] = cellNodesPtr[1];
          element_vertex_ids[6] = cellNodesPtr[2];
          element_vertex_ids[7] = cellNodesPtr[3];
          mi_->SetElemData(elem+1, Elem::HEXA8, 0, &element_vertex_ids[0]);
          break;

        case PHEX8:               
          mi_->SetElemData(elem+1, Elem::HEXA8, 0, cellNodesPtr);
          break;

        case PHEX20:               
          mi_->SetElemData(elem+1, Elem::HEXA20, 0, cellNodesPtr);
          break;

        case PHEX27:  
          mi_->SetElemData(elem+1, Elem::HEXA27, 0, cellNodesPtr);
          break;
        }

        numNodesOffset += numNodes;
      }

      break;

    case LSTRUCT: // STRUCT || LOGICALLY_STRUCT
      mi_->AddNodes(gmv_meshdata.nnodes);

      for (node=0; node<gmv_meshdata.nnodes; node++) {
        Point p;
        
        p[0] = gmv_meshdata.x[node];
        p[1] = gmv_meshdata.y[node];
        p[2] = gmv_meshdata.z[node];
        
        mi_->SetNodeCoordinate(node+1, p);
      }

      mi_->AddElems(gmv_meshdata.ncells);
      elem = 0;

      elemsx = gmv_meshdata.nxv-1;
      elemsy = gmv_meshdata.nyv-1;
      elemsz = gmv_meshdata.nzv-1;
        
      for (iz=0; iz<elemsz; iz++) {
        offset_z = iz*gmv_meshdata.nyv;
        offset_z1 = offset_z + gmv_meshdata.nyv;
        for (iy=0; iy<elemsy; iy++) {
          offset_yz = (offset_z + iy) * gmv_meshdata.nxv;
          offset_y1z = (offset_z + iy + 1) * gmv_meshdata.nxv;
          offset_yz1 = (offset_z1 + iy) * gmv_meshdata.nxv;
          offset_y1z1 = (offset_z1 + iy + 1) * gmv_meshdata.nxv;
          for (ix=0; ix<elemsx; ix++) {

            element_vertex_ids[0] = offset_yz + ix +1;
            element_vertex_ids[1] = offset_yz + ix+1 +1;
            element_vertex_ids[2] = offset_y1z + ix+1 +1;
            element_vertex_ids[3] = offset_y1z + ix +1;
            element_vertex_ids[4] = offset_yz1 + ix +1;
            element_vertex_ids[5] = offset_yz1 + ix+1 +1;
            element_vertex_ids[6] = offset_y1z1 + ix+1 +1;
            element_vertex_ids[7] = offset_y1z1 + ix +1;

            mi_->SetElemData(elem+1, Elem::HEXA8, 0, &element_vertex_ids[0]);
                    
            elem++;
          }
        }
      }

      break;
        
    case AMR:
      EXCEPTION("Loading of AMR type grids not supported at this time!");
      return false;
    }

    return true;
  }

  bool
  SimInputGMV::ProcessMaterials()
  {
    if(gmv_data.datatype != CELL)
    {
      EXCEPTION("Can only process material data on cells!");
      return false;
    }
    
    UInt numElems;
    
    mMatNames.Clear();
    mMatNums.Clear();

    for(Integer i=0; i<gmv_data.num; i++)
    {
      mMatNames.Push_back(gmv_data.chardata1+i*33);
      mMatNums.Push_back(i+1);
//      MESHIO_DEBUG(mMatNames[i]);
    }


    mi_->AddRegions(mMatNames, mMatNums);

    StdVector<Elem*> elems;
    mi_->GetElems(elems, ALL_REGIONS);
    numElems = elems.GetSize();
    
    for( UInt i=0; i<numElems; i++ ) {
      elems[i]->regionId = mMatNums[gmv_data.longdata1[i]-1];
    }

    return true;
  }

  bool
  SimInputGMV::ProcessGroups()
  {
    if((gmv_data.datatype == FACE) ||
       (gmv_data.datatype == SURF))
    {
      EXCEPTION("Can only process node or elem groups!");
      return false;
    }
    
    if(gmv_data.datatype == ENDKEYWORD)
    {
//      MESHIO_DEBUG("All groups have been read!");
      return true;
    }

    std::string entityName;
    StdVector <UInt> entities;
    entities.Resize(gmv_data.nlongdata1);

    entityName = gmv_data.name1;

    long *entpt, *entptend;
    
    entpt = gmv_data.longdata1;
    entptend = entpt + gmv_data.nlongdata1;
    
    for(UInt i = 0; entpt != entptend; entpt++, i++)
    {
      entities[i] = (UInt)(*entpt);
    }

    if(gmv_data.datatype == NODE)
    {
      mi_->AddNamedNodes(entityName, entities);
    }
    else
    {
      mi_->AddNamedElems(entityName, entities);
    }

    return true;
  }


  bool
  SimInputGMV::ProcessVariables()
  {
    /*
      ocs::VolMesh::Interfaces::CreatePtrType         mc;
      ocs::VolMesh::Interfaces::ElementAttribPtrType  ma;
      ocs::VolMesh::Interfaces::VertexAttribPtrType   mv;
      UInt attribid;

      if(gmv_data.datatype == FACE)
      {
      MESHIO_WARN("Can only process variables on nodes or cells!");
      return false;
      }

      if(gmv_data.datatype == ENDKEYWORD)
      {
      MESHIO_DEBUG("Done reading variables.");
      return true;
      }

      if( !(mc=mMesh.I.CreatePtr) ) {
      MESHIO_ERROR( "could not get IVolMeshCreate interface" );
      return false;
      }
      if( !(ma=mMesh.I.ElementAttribPtr) ) {
      MESHIO_ERROR( "could not get IVolMeshElementAttrib interface" );
      return false;
      }
      if( !(mv=mMesh.I.VertexAttribPtr) ) {
      MESHIO_ERROR( "could not get IVolMeshVertexAttrib interface" );
      return false;
      }

      std::map< std::string, GridLabelInfo >::iterator f_label;
      std::string gridlabel;
      GridLabelInfo glinfo;

      // remove trailing spaces from data name
      // PROBLEM: if spaces are at the end of the gmv_data.name1 string
      // the AdvXMLParser in ORCAN hangs the whole application when
      // switching to the attribute in VisScalarMapping

      int32 len=std::strlen(gmv_data.name1)-1;

      for( ; len >= 0; len--)
      {
      if(gmv_data.name1[len] == ' ')
      gmv_data.name1[len] = 0;
      else
      break;
      }
    
      f_label = mGridLabels.find(gmv_data.name1);
      gridlabel = gmv_data.name1;

      if(f_label != mGridLabels.end())
      {
      oc::GUID guid = oc::GUID::Generate();
      gridlabel += "_";
      gridlabel += guid.GetString();
      }

      glinfo.dim = 1;

      switch(gmv_data.datatype)
      {
      case NODE:
      if( mv->CreateVertexAttrib(gridlabel,1,"Double",(Double)0,attribid) ) {
      UInt i=0;
      std::vector< UInt >::iterator it,eit;

      glinfo.elemAttrib = false;
      glinfo.min = gmv_data.doubledata1[0];
      glinfo.max = glinfo.min;
            
      for( it=mVertexIds.begin(), eit=mVertexIds.end(); it!=eit; ++it ) {
      mv->SetVertexAttrib(attribid,*it,gmv_data.doubledata1[i]);

      glinfo.min = glinfo.min < gmv_data.doubledata1[i] ? glinfo.min : gmv_data.doubledata1[i];
      glinfo.max = glinfo.max > gmv_data.doubledata1[i] ? glinfo.max : gmv_data.doubledata1[i];
                
      i++;
      }
      }
      break;

      case CELL:
      if( ma->CreateElementAttrib(gridlabel,1,"Double",(Double)0,attribid) ) {
      UInt i=0;
      std::vector< UInt >::iterator it,eit;

      glinfo.elemAttrib = true;
      glinfo.min = gmv_data.doubledata1[0];
      glinfo.max = glinfo.min;

      for( it=mElementIds.begin(), eit=mElementIds.end(); it!=eit; ++it ) {
      ma->SetElementAttrib(attribid,*it,gmv_data.doubledata1[i]);

      glinfo.min = glinfo.min < gmv_data.doubledata1[i] ? glinfo.min : gmv_data.doubledata1[i];
      glinfo.max = glinfo.max > gmv_data.doubledata1[i] ? glinfo.max : gmv_data.doubledata1[i];
                
      i++;
      }
      }
      break;
      }
    
      MESHIO_DEBUG("variable name: " << gridlabel);
      std::pair< std::string, GridLabelInfo > value( gridlabel, glinfo );
      mGridLabels.insert(value);
    */
    return true;

  }


  bool
  SimInputGMV::ProcessProbtime()
  {
    mProbTimes[0] = gmv_data.doubledata1[0];
    return true;
  }

  bool
  SimInputGMV::ProcessCycleNo()
  {
    mCycleNos[0] = gmv_data.num;
    return true;
  }

  bool
  SimInputGMV::ProcessVelocities()
  {
    /*
      ocs::VolMesh::Interfaces::CreatePtrType         mc;
      ocs::VolMesh::Interfaces::ElementAttribPtrType  ma;
      ocs::VolMesh::Interfaces::VertexAttribPtrType   mv;
      UInt attribid;

      if(gmv_data.datatype == FACE)
      {
      MESHIO_WARN("Can only process velocities on nodes or cells!");
      return false;
      }

      if(gmv_data.datatype == ENDKEYWORD)
      {
      MESHIO_DEBUG("Done reading velocities.");
      return true;
      }

      if( !(mc=mMesh.I.CreatePtr) ) {
      MESHIO_ERROR( "could not get IVolMeshCreate interface" );
      return false;
      }
      if( !(ma=mMesh.I.ElementAttribPtr) ) {
      MESHIO_ERROR( "could not get IVolMeshElementAttrib interface" );
      return false;
      }
      if( !(mv=mMesh.I.VertexAttribPtr) ) {
      MESHIO_ERROR( "could not get IVolMeshVertexAttrib interface" );
      return false;
      }

      std::map< std::string, GridLabelInfo >::iterator f_label;
      std::string gridlabel;
      GridLabelInfo glinfo;

      gridlabel = "Velocity";
      f_label = mGridLabels.find(gridlabel);

      if(f_label != mGridLabels.end())
      {
      oc::GUID guid = oc::GUID::Generate();
      gridlabel += "_";
      gridlabel += guid.GetString();
      }

      glinfo.dim = 3;

      switch(gmv_data.datatype)
      {
      case NODE:
      if( mv->CreateVertexAttrib(gridlabel,3,"Double",ocs::Vec3d(0.,0.,0.),attribid) ) {
      UInt i=0;
      ocs::Vec3d vec;
      Double length;
      std::vector< UInt >::iterator it,eit;

      glinfo.elemAttrib = false;
            
      for( it=mVertexIds.begin(), eit=mVertexIds.end(); it!=eit; ++it ) {
      mv->SetVertexAttrib(attribid,*it,vec);

      vec[0] = gmv_data.doubledata1[i];
      vec[1] = gmv_data.doubledata2[i];
      vec[2] = gmv_data.doubledata3[i];
      length = std::sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);

      if(i == 0)
      {
      glinfo.min = length;
      glinfo.max = glinfo.min;
      }
      else
      {
      glinfo.min = glinfo.min < length ? glinfo.min : length;
      glinfo.max = glinfo.max > length ? glinfo.max : length;
      }
                
      i++;
      }
      }
      break;

      case CELL:
      if( ma->CreateElementAttrib(gridlabel,3,"Double",ocs::Vec3d(0.,0.,0.),attribid) ) {
      UInt i=0;
      ocs::Vec3d vec;
      Double length;
      std::vector< UInt >::iterator it,eit;

      glinfo.elemAttrib = true;
            
      for( it=mElementIds.begin(), eit=mElementIds.end(); it!=eit; ++it ) {
      ma->SetElementAttrib(attribid,*it,vec);

      vec[0] = gmv_data.doubledata1[i];
      vec[1] = gmv_data.doubledata2[i];
      vec[2] = gmv_data.doubledata3[i];
      length = std::sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);

      if(i == 0)
      {
      glinfo.min = length;
      glinfo.max = glinfo.min;
      }
      else
      {
      glinfo.min = glinfo.min < length ? glinfo.min : length;
      glinfo.max = glinfo.max > length ? glinfo.max : length;
      }
                
      i++;
      }
      }
      break;
      }
    
      MESHIO_DEBUG("Read velocity: " << gridlabel);
      std::pair< std::string, GridLabelInfo > value( gridlabel, glinfo );
      mGridLabels.insert(value);
    */    
    return true;
  }

  bool
  SimInputGMV::ProcessVectors()
  {
    /*
      ocs::VolMesh::Interfaces::CreatePtrType         mc;
      ocs::VolMesh::Interfaces::ElementAttribPtrType  ma;
      ocs::VolMesh::Interfaces::VertexAttribPtrType   mv;
      UInt attribid;

      if(gmv_data.datatype == FACE)
      {
      MESHIO_WARN("Can only process vectors on nodes or cells!");
      return false;
      }

      if(gmv_data.datatype == ENDKEYWORD)
      {
      MESHIO_DEBUG("Done reading vectors.");
      return true;
      }

      if( !(mc=mMesh.I.CreatePtr) ) {
      MESHIO_ERROR( "could not get IVolMeshCreate interface" );
      return false;
      }
      if( !(ma=mMesh.I.ElementAttribPtr) ) {
      MESHIO_ERROR( "could not get IVolMeshElementAttrib interface" );
      return false;
      }
      if( !(mv=mMesh.I.VertexAttribPtr) ) {
      MESHIO_ERROR( "could not get IVolMeshVertexAttrib interface" );
      return false;
      }

      std::map< std::string, GridLabelInfo >::iterator f_label;
      std::string gridlabel;
      GridLabelInfo glinfo;

      // remove trailing spaces from data name
      // PROBLEM: if spaces are at the end of the gmv_data.name1 string
      // the AdvXMLParser in ORCAN hangs the whole application when
      // switching to the attribute in VisScalarMapping

      int32 len=std::strlen(gmv_data.name1)-1;

      for( ; len >= 0; len--)
      {
      if(gmv_data.name1[len] == ' ')
      gmv_data.name1[len] = 0;
      else
      break;
      }
    
      f_label = mGridLabels.find(gmv_data.name1);
      gridlabel = gmv_data.name1;

      if(f_label != mGridLabels.end())
      {
      oc::GUID guid = oc::GUID::Generate();
      gridlabel += "_";
      gridlabel += guid.GetString();
      }

      glinfo.dim = gmv_data.num;

      Double* vec = new Double[glinfo.dim];

      switch(gmv_data.datatype)
      {
      case NODE:
      if( mv->CreateVertexAttrib(gridlabel,glinfo.dim,"Double",vec,glinfo.dim*sizeof(double),attribid) ) {
      UInt i=0;
      Double length;
      std::vector< UInt >::iterator it,eit;

      glinfo.elemAttrib = false;
            
      for( it=mVertexIds.begin(), eit=mVertexIds.end(); it!=eit; ++it ) {

      length = 0.0;
      for(UInt n=0; n<gmv_data.num; n++)
      {
      vec[n] = gmv_data.doubledata1[n*gmv_data.ndoubledata1+i];
      length += vec[n] * vec[n];
      }
      length = std::sqrt(length);

      mv->SetVertexAttrib(attribid,*it,*vec);

      if(i == 0)
      {
      glinfo.min = length;
      glinfo.max = glinfo.min;
      }
      else
      {
      glinfo.min = glinfo.min < length ? glinfo.min : length;
      glinfo.max = glinfo.max > length ? glinfo.max : length;
      }
                
      i++;
      }
      }
      break;

      case CELL:
      if( ma->CreateElementAttrib(gridlabel,glinfo.dim,"Double",vec,glinfo.dim*sizeof(double),attribid) ) {
      UInt i=0;
      Double length;
      std::vector< UInt >::iterator it,eit;

      glinfo.elemAttrib = true;
            
      for( it=mElementIds.begin(), eit=mElementIds.end(); it!=eit; ++it ) {

      length = 0.0;
      for(UInt n=0; n<gmv_data.num; n++)
      {
      vec[n] = gmv_data.doubledata1[n*gmv_data.ndoubledata1+i];
      length += vec[n] * vec[n];
      }
      length = std::sqrt(length);

      ma->SetElementAttrib(attribid,*it,*vec);

      if(i == 0)
      {
      glinfo.min = length;
      glinfo.max = glinfo.min;
      }
      else
      {
      glinfo.min = glinfo.min < length ? glinfo.min : length;
      glinfo.max = glinfo.max > length ? glinfo.max : length;
      }
                
      i++;
      }

      }
      break;
      }
    
      MESHIO_DEBUG("Read vector: " << gridlabel);
      std::pair< std::string, GridLabelInfo > value( gridlabel, glinfo );
      mGridLabels.insert(value);
    */    
    return true;
  }


  bool
  SimInputGMV::SetupGridAndRegions()
  {
    std::string regionName = "Material";
    /*
      ocs::VolMesh::Interfaces::CreatePtrType         mc;
      ocs::VolMesh::Interfaces::ElementAttribPtrType  ma;
      ocs::VolMesh::Interfaces::AttribInfoPtrType     mi;
      ocs::VolMesh::Interfaces::VertexAttribPtrType   mv;
      UInt attribid;

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

      if( ma->CreateElementAttrib(regionName,1,"UInt",(uint32)0,attribid) ) {
        
      std::vector< UInt >::iterator it4,eit4,it0;
      it0 = mElementIds.begin();
      for( it4=mElementMaterials.begin(), eit4=mElementMaterials.end(); it4!=eit4; ++it4 ) {
      ma->SetElementAttrib(attribid,*it0++,*it4);
      }
      }

      mGridLabels[regionName].dim = 1;
      mGridLabels[regionName].min = 1;
      mGridLabels[regionName].max = mMatNums.size();
      mGridLabels[regionName].elemAttrib = true;

      // build up list of attributes
      std::map< std::string, GridLabelInfo >::iterator it1, eit1;
      std::vector< std::string >::iterator it2, eit2;
      UInt length, i;
      std::map< std::string, GridLabelInfo > cfsAttribs;
      std::string dispAttribName;
      std::string finalAttribName;

      // find attributes in GMV files written with CFS++
      for(it1=mGridLabels.begin(), eit1 = mGridLabels.end(); it1 != eit1; it1++)
      {
      for(it2=mPossibleAttribs.begin(), eit2 = mPossibleAttribs.end(); it2 != eit2; it2++)
      {
      length = it2->length();

      if( it1->first.compare(0, length, *it2) == 0)
      {
      if((*it2) == "mechDisplacement")
      dispAttribName = "mechDisplacement";


      finalAttribName = (*it2) + "-Ampl";
      if( it1->first.compare(0, finalAttribName.length(), finalAttribName) != 0)
      {
      finalAttribName = (*it2) + "-Phase";
      if( it1->first.compare(0, finalAttribName.length(), finalAttribName) != 0)
      {
      finalAttribName = (*it2) + "-Real";
      if( it1->first.compare(0, finalAttribName.length(), finalAttribName) != 0)
      {
      finalAttribName = (*it2) + "-Imag";
      if( it1->first.compare(0, finalAttribName.length(), finalAttribName) != 0)
      {
      finalAttribName = *it2;
      }
      }
      }
      }

      if(it1->second.dim == 1)
      {
      cfsAttribs[finalAttribName].dim++;
      cfsAttribs[finalAttribName].elemAttrib = it1->second.elemAttrib;
                    
      MESHIO_DEBUG("Found 1D CFS++ attribute: " << finalAttribName << " dim: " << cfsAttribs[finalAttribName].dim);

      break;
      }
      else
      {
      MESHIO_DEBUG("Found CFS++ attribute " << finalAttribName << " with dim: " << it1->second.dim);
      break;
      }
                
      }
      }
      }

      // now remove the "primary" CFS++ attributes from the grid and
      // mGridLabels map and build up the the real attributes in the
      // case of 3-dimensional attributes
      UInt attribid1, attribid2, attribid3;
    
      for(it1=cfsAttribs.begin(), eit1 = cfsAttribs.end(); it1 != eit1; it1++)
      {
      if(cfsAttribs[it1->first].dim == 1)
      break;
        
      attribid1 = mi->GetAttributeID(it1->first + "1");
      attribid2 = mi->GetAttributeID(it1->first + "2");
      if(it1->second.dim == 3)
      attribid3 = mi->GetAttributeID(it1->first + "3");

      if(cfsAttribs[it1->first].elemAttrib)
      {
      //            cfsAttribs[it1->first].elemAttrib = true;
            
      if( ma->CreateElementAttrib(it1->first,3,"Double",ocs::Vec3d(0.,0.,0.),attribid) ) {
      std::vector< UInt >::iterator it,eit;
      ocs::Vec3d vec;
      Double length;
                
      i=0;

      for( it=mElementIds.begin(), eit=mElementIds.end(); it!=eit; ++it ) {
      ma->GetElementAttrib(attribid1,*it,vec[0]);
      ma->GetElementAttrib(attribid2,*it,vec[1]);
      if(it1->second.dim == 3)
      ma->GetElementAttrib(attribid3,*it,vec[2]);
      else
      vec[2] = 0.0;

      ma->SetElementAttrib(attribid,*it,vec);

      length = std::sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
                    
      if(i == 0)
      {
      cfsAttribs[it1->first].min = length;
      cfsAttribs[it1->first].max = length;
      }
                    
      cfsAttribs[it1->first].min = cfsAttribs[it1->first].min < length ? cfsAttribs[it1->first].min : length;
      cfsAttribs[it1->first].max = cfsAttribs[it1->first].max > length ? cfsAttribs[it1->first].max : length;
                
      i++;
      }

      ma->DestroyElementAttrib(attribid1);
      ma->DestroyElementAttrib(attribid2);
      if(it1->second.dim == 3)
      ma->DestroyElementAttrib(attribid3);
      }
      }
      else
      {
      // cfsAttribs[it1->first].elemAttrib = false;
            
      if( mv->CreateVertexAttrib(it1->first,3,"Double",ocs::Vec3d(0.,0.,0.),attribid) ) {
      std::vector< UInt >::iterator it,eit;
      ocs::Vec3d vec;
      Double length;
                
      i=0;

      for( it=mVertexIds.begin(), eit=mVertexIds.end(); it!=eit; ++it ) {
      mv->GetVertexAttrib(attribid1,*it,vec[0]);
      mv->GetVertexAttrib(attribid2,*it,vec[1]);
      if(it1->second.dim == 3)
      mv->GetVertexAttrib(attribid3,*it,vec[2]);
      else
      vec[2] = 0.0;

      mv->SetVertexAttrib(attribid,*it,vec);

      length = std::sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
                    
      if(i == 0)
      {
      cfsAttribs[it1->first].min = length;
      cfsAttribs[it1->first].max = length;
      }
                    
      cfsAttribs[it1->first].min = cfsAttribs[it1->first].min < length ? cfsAttribs[it1->first].min : length;
      cfsAttribs[it1->first].max = cfsAttribs[it1->first].max > length ? cfsAttribs[it1->first].max : length;
                
      i++;
      }

      mv->DestroyVertexAttrib(attribid1);
      mv->DestroyVertexAttrib(attribid2);
      if(it1->second.dim == 3)
      mv->DestroyVertexAttrib(attribid3);
      }
      }

      // Erase the 3 primary attribute's infos from mGridLabels
      mGridLabels.erase(it1->first + "1");
      mGridLabels.erase(it1->first + "2");
      if(it1->second.dim == 3)
      mGridLabels.erase(it1->first + "3");        
      }

      for(it1=cfsAttribs.begin(), eit1 = cfsAttribs.end(); it1 != eit1; it1++)
      {
      mGridLabels.insert(*it1);
      }
            
      std::vector<UInt> stepnrs;
      std::vector<Double> timesteps;
      std::vector<std::string> gridlabels;
      std::vector<std::string> attributes;
      std::vector<Double> minimums;
      std::vector<Double> maximums;

      // get TimeRegionMetaInfo Interface
      LSETimeRegionMeshInterface* trmi=NULL;

      try { 
      if( mMesh.GetProperties().HasProperty("LSETimeRegionMeshInterface") ) {
      oc::CPtr<LSETimeRegionMeshInterface> cptr = mMesh.GetProperties()["LSETimeRegionMeshInterface"](FIELD_ID(oc::CPtr<LSETimeRegionMeshInterface>));
      trmi = cptr.GetCPtr();
      }
      }
      catch(...) {
      MESHIO_WARN("Could not get LSETimeRegionMeshInterface");
      return false;
      }

      trmi->SetRegionAttrName(regionName);
      trmi->SetDisplacementAttrName(dispAttribName);
      trmi->SetRegions(mMatNums);
      trmi->SetRegionNames(mMatNames);

      for(it1=mGridLabels.begin(), eit1 = mGridLabels.end(); it1 != eit1; it1++)
      {
      MESHIO_DEBUG("gridlabels: " << it1->first);

      attributes.push_back(it1->first);
      }

      trmi->SetAttrList(attributes);
    
      for(it1=mGridLabels.begin(), eit1 = mGridLabels.end(); it1 != eit1; it1++)
      {
      trmi->SetAttrRegions(it1->first, mMatNums);
      stepnrs.push_back(mCycleNos[0]);
      timesteps.push_back(mProbTimes[0]);
      gridlabels.push_back(it1->first);
      minimums.push_back(it1->second.min);
      maximums.push_back(it1->second.max);
      for(i=0; i<mMatNums.size(); i++)
      {
      trmi->SetAttrTimeStepInfos(it1->first, mMatNums[i], stepnrs, timesteps, gridlabels, minimums, maximums);
      }

      stepnrs.clear();
      timesteps.clear();
      gridlabels.clear();
      minimums.clear();
      maximums.clear();
      }
    */    
    return true;
  }
}
