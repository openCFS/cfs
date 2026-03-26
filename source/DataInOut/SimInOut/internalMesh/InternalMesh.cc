// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <string>
#include <algorithm>

#include "InternalMesh.hh"
#include "DataInOut/ParamHandling/XmlReader.hh"
#include "DataInOut/SimInOut/AnsysFile/SimInputMESH.hh"

namespace CoupledField
{
using std::string;
using std::cout;
using std::endl;

InternalMesh::InternalMesh(std::string fileName, PtrParamNode inputNode,
                           PtrParamNode infoNode) :
            SimInput(fileName, inputNode, infoNode),
            dim_(0),
            maxNumElems_(0),
            maxNumNodes_(0),
            info_(infoNode->Get("header")->Get("domain")->Get("internal"))
{
  
  // initialize minimal / maximal point
  minimal_.Resize(3);
  minimal_.Init(0.0);
  maximal_.Resize(3);
  maximal_[0] = 1.0;  // maxz is set to 0.0 for 2D
  maximal_[1] = 1.0;   // maxz is set to 0.0 for 2D
  
  // in base class
  capabilities_.insert(SimInput::MESH);

  // parse the file and create corresponding paramnode
  xml_ = XmlReader::ParseFile(fileName_);

  
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

    maximal_[2] = 1.0; // set default maximal z to 1.0
  }
  else
    nelems_[2] = 0;

  // set the region name
  if(xml_->Has("region"))
    regionNames_.Push_back(xml_->Get("region")->Get("name")->As<string>());
  else
    regionNames_.Push_back("mech"); // default region name

  // FIXME put this in the log output
  info_->Get("basic")->Get("nx")->SetValue(nelems_[0]);
  info_->Get("basic")->Get("ny")->SetValue(nelems_[1]);
  info_->Get("basic")->Get("nz")->SetValue(nelems_[2]);
  
  cout << "dim = " << dim_
       << ", nx = " << nelems_[0]
       << ", ny = " << nelems_[1]
       << ", nz = " << nelems_[2]
       << ", numElems = " << maxNumElems_
       << ", numNodes = " << maxNumNodes_
       << ", regionname = " << string(regionNames_[0])
       << endl;

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

    // sanity checks
    assert(minimal_[0] < maximal_[0]);
    assert(minimal_[1] < maximal_[1]);
    assert(dim_ == 3 ? (minimal_[2] < maximal_[2]) : true);
  }

  if(xml_->Has("boundary"))
    ParseBoundary(xml_->Get("boundary"));
    
  PtrParamNode dims_ = info_->Get("coordinate_dimensions");
  dims_->Get("minimal")->Get("x")->SetValue(minimal_[0]);
  dims_->Get("minimal")->Get("y")->SetValue(minimal_[1]);
  dims_->Get("minimal")->Get("z")->SetValue(minimal_[2]);
  dims_->Get("maximal")->Get("x")->SetValue(maximal_[0]);
  dims_->Get("maximal")->Get("y")->SetValue(maximal_[1]);
  dims_->Get("maximal")->Get("z")->SetValue(maximal_[2]);
}

void InternalMesh::ReadMesh(Grid *mi)
{
  mi_ = mi;

  // Get Nodes
  mi_->AddNodes(maxNumNodes_);

  const double xinc((maximal_[0] - minimal_[0])/nelems_[0]);
  const double yinc((maximal_[1] - minimal_[1])/nelems_[1]);
  const double zinc(dim_ == 3 ? ((maximal_[2] - minimal_[2])/nelems_[2]) : 0.0);

  // write node coordinates directly into the grid
  // regionNodes_ must contain all the node numbers and that can be done here, too
  // because we only have one region!
  regionNodes_.Push_back(std::set<UInt>());
  UInt nodeNum(0);
  for(UInt z = 0; z <= nelems_[2]; ++z)
    for(UInt y = 0; y <= nelems_[1]; ++y)
      for(UInt x = 0; x <= nelems_[0]; ++x)
      {
        Vector<Double> loc(3);
        loc[0] = xinc*x;
        loc[1] = yinc*y;
        loc[2] = zinc*z;
        mi_->SetNodeCoordinate(++nodeNum, loc);
        regionNodes_[0].insert(nodeNum);
      }

  // Get Regions
  RegionIdType rId;
  // add the one region and directly set regular to true
  rId = mi_->AddRegion(regionNames_[0], true);


  // Get Elements
  mi_->AddElems(maxNumElems_);

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

  GetElements(elems, elemTypes, elemNums, regionId, dim_);

  UInt n(0);
  for(UInt j = 0; j < elems.GetSize(); ++j)
  {
    n=0;
    for(UInt i = 0; i < elemTypes[j].GetSize(); ++i)
    {
      mi_->SetElemData(elemNums[j][i], elemTypes[j][i], rId, &elems[j][n]);
      n += Elem::shapes[elemTypes[j][i]].numNodes;
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
  // FIXME this must be adjusted to what named nodes are defined
  // in the xml file
  // for now we use the ones for homogenization boundaries
  UInt num(0);
  switch(dim_)
  {
  case 2:
    num = 2*(nelems_[0] + 1) + 2*(nelems_[1] + 1);
    break;
  case 3:
    num =  2*(nelems_[0] + 1)*(nelems_[1] + 1);
    num += 2*(nelems_[0] + 1)*(nelems_[2] + 1);
    num += 2*(nelems_[1] + 1)*(nelems_[2] + 1);
    break;
  default:
    assert(false);
  }

  return num;
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


// ======================================================
// ENTITY ACCESS
// ======================================================


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
  UInt eType(6);
  // UInt eNodes(4);
  RegionIdType regionId(0);
  Integer regionIndex(0);

  if(dim_ == 3)
  {
    eType = 10;
    // eNodes = 8;
  }

  // prepare the vectors
  // we only need to do this once because we only have one region!
  regionIds.Push_back(regionId);
  elems.Push_back( StdVector<UInt>() );
  elems.Last().Reserve(maxNumElems_);

  elemNums.Push_back( StdVector<UInt>() );
  elemNums.Last().Reserve(maxNumElems_);

  elemTypes.Push_back( StdVector<Elem::FEType>() );
  elemTypes.Last().Reserve(maxNumElems_);


  UInt eNum(0);
  UInt num(0);
  // std::cout << std::endl; // FIXME

  for(UInt z = 0; z < (dim_ == 2 ? 1 : nelems_[2]); ++z)
    for(UInt y = 0; y < nelems_[1]; ++y)
      for(UInt x = 0; x < nelems_[0]; ++x)
      {
        ++eNum;

        switch(dim)
        {
        case 2:
        {
          num = y*(nelems_[0]+1) + x + 1;
          elems[regionIndex].Push_back(num);
          ++num;
          elems[regionIndex].Push_back(num);

          num = (y+1)*(nelems_[0]+1) + x + 2;
          elems[regionIndex].Push_back(num);
          --num;
          elems[regionIndex].Push_back(num);

          /* FIXME put this into the log file
						std::cout << "elem " << eNum
											<< ": n1 = " << y*(nelems_[0]+1) + x + 1
											<< ", n2 = " << y*(nelems_[0]+1) + x + 2
											<< ", n3 = " << (y+1)*(nelems_[0]+1) + x + 2
											<< ", n4 = " << (y+1)*(nelems_[0]+1) + x + 1
											<< std::endl;
           */

          break;
        }

        case 3:
        {
          num = ((z+1)*(nelems_[1]+1) + y)*(nelems_[0]+1) + x + 1;
          elems[regionIndex].Push_back(num);
          num = ((z+1)*(nelems_[1]+1) + y+1)*(nelems_[0]+1) + x + 1;
          elems[regionIndex].Push_back(num);
          ++num;
          elems[regionIndex].Push_back(num);
          num = ((z+1)*(nelems_[1]+1) + y)*(nelems_[0]+1) + x + 2;
          elems[regionIndex].Push_back(num);

          num = (z*(nelems_[1]+1) + y)*(nelems_[0]+1) + x + 1;
          elems[regionIndex].Push_back(num);
          num = (z*(nelems_[1]+1) + y+1)*(nelems_[0]+1) + x + 1;
          elems[regionIndex].Push_back(num);
          ++num;
          elems[regionIndex].Push_back(num);
          num = (z*(nelems_[1]+1) + y)*(nelems_[0]+1) + x + 2;
          elems[regionIndex].Push_back(num);

          break;
        }

        default:
          EXCEPTION("not supported dimension for internal mesh");
        }

        elemTypes[regionIndex].Push_back(SimInputMESH::AnsysType2ElemType(eType));
        elemNums[regionIndex].Push_back(eNum);
      }

  // Set flag which indicates, that elements of given dimension
  // were read in
  elemDimReadIn_[dim-1] = true;
}

void InternalMesh::GetNamedNodes(StdVector<StdVector<UInt> > &nodes,
    StdVector<string> &nodeNames) 
{
  nodeNames.Reserve(dim_ == 2 ? 4 : 6);

  UInt nel(0);
  switch(dim_)
  {
  case 2:
  {
    nel = nelems_[0]+1;
    // bottom
    nodeNames.Push_back("bottom");
    nodes.Push_back( StdVector<UInt>(nel) );
    for(UInt i = 0; i < nel; ++i)
      nodes[0][i] = i + 1;

    // top
    nodeNames.Push_back("top");
    nodes.Push_back( StdVector<UInt>(nel) );
    for(UInt i = 0; i < nel; ++i)
      nodes[1][i] = nelems_[1]*nel + i + 1;

    nel = nelems_[1]+1;

    // left
    nodeNames.Push_back("left");
    nodes.Push_back( StdVector<UInt>(nel) );
    for(UInt i = 0; i < nel; ++i)
      nodes[2][i] = i*(nelems_[0]+1) + 1;

    // right
    nodeNames.Push_back("right");
    nodes.Push_back( StdVector<UInt>(nel) );
    for(UInt i = 0; i < nel; ++i)
      nodes[3][i] = (i+1)*(nelems_[0]+1);

    break;
  }
  case 3:
  {
    UInt num(0);

    // note that the node number is x+(nx+1)*(y+(ny+1)*z)

    // left (x=minimal_)
    nodeNames.Push_back("left");
    nel = (nelems_[1]+1) * (nelems_[2]+1);
    nodes.Push_back( StdVector<UInt>(nel) );
    UInt c(0);
    for(UInt z = 0; z < nelems_[2]+1; ++z)
      for(UInt y = 0; y < nelems_[1]+1; ++y)
      {
        num = (z*(nelems_[1]+1)+y)*(nelems_[0]+1) + 1;
        assert(num < maxNumNodes_);
        nodes[0][c++] = num;
      }

    // right (x=maximal_)
    nodeNames.Push_back("right");
    nodes.Push_back( StdVector<UInt>(nel) );
    c = 0;
    for(UInt z = 0; z < nelems_[2]+1; ++z)
      for(UInt y = 0; y < nelems_[1]+1; ++y)
      {
        num = (z*(nelems_[1]+1)+y)*(nelems_[0]+1) + nelems_[0] + 1;
        assert(num <= maxNumNodes_);
        nodes[1][c++] = num;
      }

    // bottom (y=minimal_)
    nodeNames.Push_back("bottom");
    nel = (nelems_[0]+1) * (nelems_[2]+1);
    nodes.Push_back( StdVector<UInt>(nel) );
    c = 0;
    for(UInt z = 0; z < nelems_[2]+1; ++z)
      for(UInt x = 0; x < nelems_[0]+1; ++x)
      {
        num = z*(nelems_[1]+1)*(nelems_[0]+1) + x + 1;
        assert(num <= maxNumNodes_);
        nodes[2][c++] = num;
      }

    // top (y=maximal_)
    nodeNames.Push_back("top");
    nodes.Push_back( StdVector<UInt>(nel) );
    c = 0;
    for(UInt z = 0; z < nelems_[2]+1; ++z)
      for(UInt x = 0; x < nelems_[0]+1; ++x)
      {
        num = (nelems_[0]+1)*nelems_[1] + z*(nelems_[1]+1)*(nelems_[0]+1) + x + 1;
        assert(num <= maxNumNodes_);
        nodes[3][c++] = num;
      }

    // note that the node number is x+(nx+1)*y
    // front (z=minimal_)
    nodeNames.Push_back("front");
    nel = (nelems_[0]+1) * (nelems_[1]+1);
    nodes.Push_back( StdVector<UInt>(nel) );
    c = 0;
    for(UInt y = 0; y < nelems_[1]+1; ++y)
      for(UInt x = 0; x < nelems_[0]+1; ++x)
      {
        num = y*(nelems_[0]+1) + x + 1;
        assert(num <= maxNumNodes_);
        nodes[4][c++] = num;
      }

    // back (z=maximal_)
    nodeNames.Push_back("back");
    nodes.Push_back( StdVector<UInt>(nel) );
    c = 0;
    for(UInt y = 0; y < nelems_[1]+1; ++y)
      for(UInt x = 0; x < nelems_[0]+1; ++x)
      {
        num = (nelems_[0]+1)*(nelems_[1]+1)*nelems_[2] + y*(nelems_[0]+1) + x + 1;
        assert(num <= maxNumNodes_);
        nodes[5][c++] = num;
      }

    break;
  }
  default:
    EXCEPTION("not supported dimension for internal mesh");
  }
}

void InternalMesh::ParseBoundary(PtrParamNode bdr)
{
  ParamNodeList tmp;

  if(bdr->Has("nodes"))
  {
    tmp = bdr->Get("nodes")->GetChildren();
    cout << "nodes.size = " << tmp.GetSize() << endl;
  }

  if(bdr->Has("edges"))
  {
    tmp = bdr->Get("edges")->GetChildren();
    cout << "edges.size = " << tmp.GetSize() << endl;
  }

  if(bdr->Has("surfaces")) // FIXME && dim_ == 2) // ignore surfaces in 2D
  {
    tmp = bdr->Get("surfaces")->GetChildren();
    cout << "surfaces.size = " << tmp.GetSize() << endl;
  }
}
}
