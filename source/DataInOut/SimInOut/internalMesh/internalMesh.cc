// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <string>
#include <algorithm>

#include "internalMesh.hh"
#include "DataInOut/ParamHandling/Xerces.hh"


namespace CoupledField
{
using std::string;

InternalMesh::InternalMesh(std::string fileName, PtrParamNode inputNode) :
        SimInput(fileName, inputNode),
        dim_(0),
        minimal_(Point(0.0, 0.0, 0.0)),
        maximal_(Point(1.0, 1.0, 0.0)), // maxz is set to 0.0 for 2D
        maxNumElems_(0),
        maxNumNodes_(0)
{
    // in base class
    capabilities_.insert(SimInput::MESH);
    
    // create temporary xerces file parser
    Xerces* xerces = new Xerces(fileName_);
    // parse the file and create corresponding paramnode
    xml_ = xerces->CreateParamNodeInstance();
    // delete the xerces object, does not affect the paramnode!
    delete xerces;
    
    if(xml_->GetName() != string("cfsInternalMesh"))
      EXCEPTION("file " << fileName_ << " is not a cfsInternalMesh-file!");
    
    xml_ = xml_->Get("grid"); 

    xml_->GetValue("dim", dim_);
		assert(dim_ == 2 || dim_ == 3);

		xml_->GetValue("nx", nelems_[0]);
		xml_->GetValue("ny", nelems_[1]);
		maxNumElems_ = nelems_[0] * nelems_[1];
		maxNumNodes_ = (nelems_[0]+1) * (nelems_[1]+1);

		assert(maxNumElems_ > 0);
		if(dim_ == 3)
		{
			assert(xml_->Has("nz"));
			xml_->GetValue("nz", nelems_[2]);
			maxNumElems_ *= nelems_[2];
			maxNumNodes_ *= (nelems_[2]+1);
			assert(maxNumElems_ > 0);
			
			maximal_[2] = 1.0; // set maximal z to 1.0
		}
    else
      nelems_[2] = 0;

    // set the region name
    if(xml_->Has("region"))
		  regionNames_.Push_back(xml_->Get("region")->Get("name")->As<string>());
    else
		  regionNames_.Push_back("mech"); // default region name
    
    // FIXME put this in the log output
    std::cout << "dim = " << dim_
							<< ", nx = " << nelems_[0]
							<< ", ny = " << nelems_[1]
							<< ", nz = " << nelems_[2]
							<< ", numElems = " << maxNumElems_
							<< ", numNodes = " << maxNumNodes_
							<< ", regionname = " << string(regionNames_[0])
							<< std::endl;

    elemDimReadIn_.Resize(dim_);
  }

  void InternalMesh::InitModule()
  {
    if(xml_->Has("dimensions"))
    {
      PtrParamNode dims = xml_->Get("dimensions");
      if(dims->Has("minimal"))
      {
        dims->Get("minimal")->GetValue("x", minimal_[0]);
        dims->Get("minimal")->GetValue("y", minimal_[1]);
        if(dim_ == 3)
          dims->Get("minimal")->GetValue("z", minimal_[2]);
      }

      if(dims->Has("maximal"))
      {
        dims->Get("maximal")->GetValue("x", maximal_[0]);
        dims->Get("maximal")->GetValue("y", maximal_[1]);
        if(dim_ == 3)
          dims->Get("maximal")->GetValue("z", maximal_[2]);
      }
    }

    if(xml_->Has("boundary"))
    {
      ParseBoundary(xml_->Get("boundary"));
    }
  }
    
  void InternalMesh::ReadMesh(Grid *mi)
  {
    mi_ = mi;
    
    // Get Nodes
    mi_->AddNodes(maxNumNodes_);

		const double xinc((maximal_[0] - minimal_[0])/nelems_[0]);
		const double yinc((maximal_[1] - minimal_[1])/nelems_[1]);

		// write node coordinates directly into the grid
		for(UInt y = 0; y <= nelems_[0]; ++y)
    {
			for(UInt x = 0; x <= nelems_[0]; ++x)
			{
				// FIXME for 2D only atm
				assert(dim_ == 2);
				const UInt nodeNum((nelems_[0]+1)*y + x + 1);
				assert(nodeNum > 0 && nodeNum <= GetNumNodes());

				mi_->SetNodeCoordinate(nodeNum, Point(xinc*x, yinc*y, 0.0));
			}
    }


    // Get Regions
		RegionIdType rId;
    // add the one region and directly set regular to true
		rId = mi_->AddRegion(regionNames_[0], true);





    
    // Get Elements
    UInt numElems = 0;
    UInt numElems1D = 0; //GetNumElems(1);
    UInt numElems2D = maxNumElems_; // FIXME
    UInt numElems3D = 0; //GetNumElems(3);
    numElems = numElems1D + numElems2D + numElems3D;

    mi_->AddElems(numElems);
    
    // TODO: set the element neighborhood here!
    /*
    for(UInt i = 0; i < maxNumElems_; ++i)
    {
      // get pointer to element
      const Elem *el = mi_->GetElem(i+1);
      
    }
    */

    StdVector<StdVector<UInt> > elems;
    StdVector<StdVector<UInt> > elemNums;
    StdVector<StdVector<Elem::FEType> > elemTypes;
    StdVector<RegionIdType> regionId;

    GetElements(elems,elemTypes,elemNums,regionId,2);
    
    UInt n(0);
    for(UInt j = 0; j < elems.GetSize(); ++j)
    {
      n=0;
      for(UInt i = 0; i < elemTypes[j].GetSize(); ++i)
      {
        mi_->SetElemData(elemNums[j][i], elemTypes[j][i], 
                         rId, &elems[j][n]);
        n += Elem::GetNumElemNodes(elemTypes[j][i]);
      }
    }





    // Get Named Nodes
    StdVector<StdVector<UInt> > indices;
    StdVector<string> nodeNames;

    GetNamedNodes(indices, nodeNames);
    
    // FIXME put this in the log output
    // for(UInt i = 0; i < indices.GetSize(); ++i)
      // std::cout << std::endl
                // << nodeNames[i] << ": " 
                // << indices[i].ToString() 
                // << std::endl;
    
    for(UInt i = 0, ss = nodeNames.GetSize(); i < ss; ++i)
      mi_->AddNamedNodes(nodeNames[i], indices[i]);







    // Get Named Elements
		/* no named elements atm
    std::vector< std::string > elemNames;

    GetElemNames( names );
    indices.clear();
    GetNamedElems( indices, elemNames );

    names.Clear();
    for(UInt i = 0; i<elemNames.size(); i++)
      names.Push_back(elemNames[i]);

    for(UInt i = 0; i<elemNames.size(); i++) {
      mi_->AddNamedElems(names[i], indices[i]);
    }
		*/
  }
    
  UInt InternalMesh::GetNumElems(const Integer dim)
  {
    return maxNumElems_;
  }
  
  UInt InternalMesh::GetNumRegions()
  {
		// we only have one region
    return 1;
  }

  UInt InternalMesh::GetNumNamedNodes()
  {
    // FIXME
    return 2*(nelems_[0] + 1) + 2*(nelems_[1] + 1);
  }

  UInt InternalMesh::GetNumNamedElems()
  {
    //return GetInteger("NumSaveElements");
    return 0;
  }
  
  // ======================================================
  // ENTITY NAME ACCESS
  // ======================================================

  void InternalMesh::GetAllRegionNames(StdVector<string> &regionNames)
  {
    regionNames = regionNames_;
  }
    
  void InternalMesh::GetRegionNamesOfDim(StdVector<string> &regionNames,
                                         const UInt dim) 
  {
    regionNames.Clear();

    // Check if elements of desired dimension were read in. If not,
    // read them in into dummy variables
    if(elemDimReadIn_[dim-1] == false)
    {
      StdVector<StdVector<UInt> > elems, elemNums;
      StdVector<StdVector<Elem::FEType> > elemTypes;
      StdVector<RegionIdType> dummyId;
      GetElements(elems,elemTypes,elemNums,dummyId,dim);
    }
    
    // Look for region names of desired dimension
    for(UInt i=0; i<regionDim_.GetSize(); i++) 
      if(regionDim_[i] == dim)
        regionNames.Push_back(regionNames_[i]);
    
  }

  void InternalMesh::GetNodeNames(StdVector<string> &nodeNames) 
	{
    nodeNames.Reserve(dim_ == 2 ? 4 : 6);
    // FIXME hard coded node names
    nodeNames.Push_back("bottom");
    nodeNames.Push_back("top");
    nodeNames.Push_back("left");
    nodeNames.Push_back("right");

    if(dim_ == 3)
    {
      nodeNames.Push_back("front");
      nodeNames.Push_back("back");
    }
	}

  void InternalMesh::GetElemNames(StdVector<string> &elemNames) 
	{
		// no named element for now
  }

  // ======================================================
  // ENTITY ACCESS
  // ======================================================
  
  void InternalMesh::GetNodesOfRegions(StdVector<StdVector<UInt> > &nodes,
                                     const StdVector<RegionIdType> &regionId) 
	{
    std::set<UInt>::iterator it;
    UInt index, iNode;
    
    nodes.Resize(regionId.GetSize());

    for(UInt iRegion = 0, ss = regionId.GetSize(); iRegion < ss; ++iRegion) 
		{
      iNode = 0;
      index = regionId[iRegion];
      nodes[iRegion].Resize(regionNodes_[index].size());

      for(it = regionNodes_[index].begin(); it != regionNodes_[index].end();
          ++it, ++iNode)
      {
        nodes[iRegion][iNode] = *it;
      }
    }
  }
    
  void InternalMesh::GetElements(StdVector<StdVector<UInt> > & elems,
                                 StdVector<StdVector<Elem::FEType> > & elemTypes,
                                 StdVector<StdVector<UInt> > & elemNums,                                
                                 StdVector<RegionIdType> & regionIds,
                                 const UInt dim)
  {  
    // Check that dimension is correct
    if(dim < 1 || dim > 3)
      EXCEPTION("The dimension of elements to be read in was specified with "
                << dim << "but is only allowed to have a value between 1 and 3!");
    
    // If there are no elements, we assume that this is fine and
    // simply return
    if(maxNumElems_ <= 0)
      return;

    // Some additional variables
    UInt i(0), eType(6), eNodes(4);
    string lastRegion;
    RegionIdType regionId(-1);
    Integer regionIndex(0);
    
    // Loop over all elements
    if(dim_ == 3)
    {
      eType = 10; // FIXME
      eNodes = 8;
    }
    
    UInt eNum(0);
    for(i = 0; i < maxNumElems_; ++i)
    {
      eNum += 1;
      
      // Check if previous element had the same id. 
      // If not, obtain new region identifier
      if(lastRegion != regionNames_[0])
      {
        lastRegion = regionNames_[0];
        regionId = 0;
        
        // Check if region of this type already exists, and if not
        // add new vector
        StdVector<RegionIdType>::iterator it, end;
                
        end = regionIds.End();

        it = std::find(regionIds.Begin(), end, regionId);
        
        if(it == end)
        {
          regionIds.Push_back(regionId);
          elems.Push_back( StdVector<UInt>() );
          elems.Last().Reserve(maxNumElems_);
          elemNums.Push_back( StdVector<UInt>() );
          elemNums.Last().Reserve(maxNumElems_);
          elemTypes.Push_back( StdVector<Elem::FEType>() );
          elemTypes.Last().Reserve(maxNumElems_);
          regionNodes_.Push_back(std::set<UInt>());
          regionIndex = regionIds.GetSize() - 1;
        }
        else
        {
          regionIndex = std::distance(regionIds.Begin(), it );
        }
      }
      
      if(eNum % (nelems_[0]+1) == 0)
        eNum += 1;
      
      // Read node numbers and insert them into the element and
      // into the vector with all node-numbers per region
      elems[regionIndex].Push_back(eNum);
      regionNodes_[regionId].insert(eNum);
      elems[regionIndex].Push_back(eNum + 1);
      regionNodes_[regionId].insert(eNum + 1);
      elems[regionIndex].Push_back(eNum + (nelems_[0]+1) + 1);
      regionNodes_[regionId].insert(eNum + (nelems_[0]+1) + 1);
      elems[regionIndex].Push_back(eNum + (nelems_[0]+1));
      regionNodes_[regionId].insert(eNum + (nelems_[0]+1));

//      std::cout << "elem " << i+1
//          << ": n1 = " << eNum
//          << ", n2 = " << eNum+1
//          << ", n3 = " << eNum + (nx_+1) + 1
//          << ", n4 = " << eNum + (nx_+1)
//          << std::endl;
          
      elemTypes[regionIndex].Push_back(AnsysType2ElemType(eType));
      elemNums[regionIndex].Push_back(i+1);
    }

    // Set flag which indicates, that elements of given dimension
    // were read in
    elemDimReadIn_[dim-1] = true;
  }

  void InternalMesh::GetNamedNodes(StdVector<StdVector<UInt> > &nodes,
                                   StdVector<string> &nodeNames) 
	{
		nodeNames.Reserve(dim_ == 2 ? 4 : 6);
		
		// bottom
		nodeNames.Push_back("bottom");
		nodes.Push_back( StdVector<UInt>(nelems_[0]+1) );
		for(UInt i = 0; i <= nelems_[0]; ++i)
		  nodes[0][i] = i + 1;

		// top
		nodeNames.Push_back("top");
		nodes.Push_back( StdVector<UInt>(nelems_[0]+1) );
		for(UInt i = 0; i <= nelems_[0]; ++i)
		  nodes[1][i] = nelems_[1] * (nelems_[0]+1) + i + 1;

		// left
		nodeNames.Push_back("left");
		nodes.Push_back( StdVector<UInt>(nelems_[1]+1) );
		for(UInt i = 0; i <= nelems_[1]; ++i)
		  nodes[2][i] = i*(nelems_[0]+1) + 1;

		// right
		nodeNames.Push_back("right");
		nodes.Push_back( StdVector<UInt>(nelems_[1]+1) );
		for(UInt i = 0; i <= nelems_[1]; ++i)
		  nodes[3][i] = (i+1)*(nelems_[0]+1);
  }

  void InternalMesh::GetNamedElems(StdVector<StdVector<UInt> > & elems,
                                   StdVector<string> & elemNames) 
	{
		// no named elements atm
  }
 
  Elem::FEType InternalMesh::AnsysType2ElemType(const UInt itype) 
	{
    switch(itype)
		{
    case 101:
      return Elem::LINE3;
    case 100:
      return Elem::LINE2;
    case 4:
      return Elem::TRIA3;
    case 5:
      return Elem::TRIA6;
    case 6:
      return Elem::QUAD4;
    case 7:
      return Elem::QUAD8;
    case 107:
      return Elem::QUAD9;
    case 8:
      return Elem::TET4;
    case 9:
      return Elem::TET10;
    case 10:
      return Elem::HEXA8;
    case 11:
      return Elem::HEXA20;
    case 111:
      return Elem::HEXA27;
    case 12:
      return Elem::PYRA5;
    case 13:
      return Elem::PYRA13;
    case 14:
      return Elem::WEDGE6;
    case 15:
      return Elem::WEDGE15;
    default:
			// This place should never be reached!
			return Elem::UNDEF;
    }
  }

  void InternalMesh::ParseBoundary(PtrParamNode bdr)
  {
    ParamNodeList tmp;

    if(bdr->Has("nodes"))
    {
      tmp = bdr->Get("nodes")->GetChildren();
      std::cout << "nodes.size = " << tmp.GetSize() << std::endl;
    }

    if(bdr->Has("edges"))
    {
      tmp = bdr->Get("edges")->GetChildren();
      std::cout << "edges.size = " << tmp.GetSize() << std::endl;
    }

    if(bdr->Has("surfaces")) // FIXME && dim_ == 2) // ignore surfaces in 2D
    {
      tmp = bdr->Get("surfaces")->GetChildren();
      std::cout << "surfaces.size = " << tmp.GetSize() << std::endl;
    }
  }
}
