#include <cmath>
#include <string>

#include "SimInputRegular.hh"
#include "Domain/Mesh/GridCFS/GridCFS.hh"
#include "General/EnvironmentTypes.hh"
#include "General/Exception.hh"
#include "Utils/ToolsFull.hh" 
#include "DataInOut/ProgramOptions.hh"
#include "DataInOut/ParamHandling/XmlReader.hh"
#include "DataInOut/Logging/LogConfigurator.hh"

// declare class specific logging stream
DEFINE_LOG(regular, "regular")

namespace CoupledField
{
using std::string;

SimInputRegular::SimInputRegular(std::string fileName, PtrParamNode pn, PtrParamNode info) :
            SimInput(fileName, pn, info)
{
  capabilities_.insert(SimInput::MESH);
  info_ = info->Get("header/input/regular");

  n_.Resize(3, 0);
  h_.Resize(3, 0.0);
  min_.Resize(3, 0.0);

  assert(!(fileName.empty() && pn == nullptr)); // we need either a file or a node or both
  if(fileName_.empty())
  {
    // keep the xml; the actual parsing happens in ReadMesh() where the global domain
    // math parser (used by MathParse and the GridCFS Select* functions) is initialized
    xml_ = pn;
  }
  else // we have given file name for an external file
  {
    // check that we have only attributes but no child elements
    for(const auto& child : pn->GetChildren())
      if(child->GetChildren().GetSize() > 0)
        throw Exception("regular: fileName is given, but child " + child->GetName() + " is given.");
  
    // as sanity check of the input xml, we parse first without schema (shall be very fast)
    auto root = XmlReader::ParseFile(fileName_);
    if(root->GetName() != "cfsMeshDescription")
      throw Exception("regular: file " + fileName_ + " does not have 'cfsMeshDescription' as root element.");

    root = XmlReader::ParseFile(fileName_, progOpts->GetSchemaPathStr() + "/CFS-Simulation/CFS.xsd", "http://www.cfs++.org/simulation");
    xml_ = root->Get("regular");
    if(xml_->Get("fileName")->As<std::string>() != "default")
      throw Exception("regular: fileName in file " + fileName_ + " is not allowed.");
  }

  // pn can be null if we are loaded via -m, but then we have a fileName
  assert(xml_->Get("geometryType")->As<std::string>() == "plane" || xml_->Get("geometryType")->As<std::string>() == "3d");
  dim_ = (xml_->Get("geometryType")->As<std::string>() == "plane") ? 2 : 3;
  
  // whenever region is given in the problem.xml, it overwrites the one in an external file.
  region_ = (pn != nullptr ? pn : xml_)->Get("region")->As<std::string>();
}

void SimInputRegular::ParseXML(const PtrParamNode& pn)
{
  const PtrParamNode& box = pn->Get("box"); // all values have defaults
  Vector<double> max(3); // max needs not to be stored (min_ + n_*h_ is enough), but we do sanity checks.
  min_[0] = box->Get("xMin")->MathParse<double>();
  max[0]  = box->Get("xMax")->MathParse<double>();
  min_[1] = box->Get("yMin")->MathParse<double>();
  max[1]  = box->Get("yMax")->MathParse<double>();
  min_[2] = box->Get("zMin")->MathParse<double>(); // we have defaults - ignored in 2D
  max[2]  = box->Get("zMax")->MathParse<double>();
  if(min_[0] >= max[0] || min_[1] >= max[1] || (dim_ == 3 && min_[2] >= max[2]))
    throw Exception("regular mesh box maximal coordinates must be larger than the minimals!");
  
  LOG_DBG(regular) << "PX: min=" << min_.ToString() << " max=" << max.ToString();

  // --- exactly one of <elements> / <spacing>, guaranteed by the schema choice ---
  // the given quantity is taken as-is, the other is derived; the first axis (x) drives
  // the default (cubic elements) for any y/z the user omitted.
  n_.Init(0);   // ensures size 3 and zeroes the z-component that stays unused in 2D
  h_.Init(0.0);
  if(pn->Has("elements"))
  {
    n_[0] = pn->Get("elements/nx")->MathParse<int>(); // mandatory, positive per schema
    h_[0] = (max[0] - min_[0]) / n_[0];

    // parse or construct
    n_[1] = pn->Has("elements/ny") ? pn->Get("elements/ny")->MathParse<int>() : std::lround((max[1] - min_[1]) / h_[0]);
    h_[1] = (max[1] - min_[1]) / n_[1];

    if(dim_ == 3) {
      n_[2] = pn->Has("elements/nz") ? pn->Get("elements/nz")->MathParse<int>() : std::lround((max[2] - min_[2]) / h_[0]);
      h_[2] = (max[2] - min_[2]) / n_[2];
    }
    LOG_DBG(regular) << "PX: elements n=" << n_.ToString() << " -> h=" << h_.ToString();
  }
  else // <spacing>
  {
    const PtrParamNode& s = pn->Get("spacing");
    h_[0] = s->Get("hx")->MathParse<double>(); // mandatory, positive per schema
    n_[0] = std::lround((max[0] - min_[0]) / h_[0]);

    h_[1] = s->Has("hy") ? s->Get("hy")->MathParse<double>() : (max[1] - min_[1]) / n_[0];
    n_[1] = std::lround((max[1] - min_[1]) / h_[1]);

    if(dim_ == 3) {
      h_[2] = s->Has("hz") ? s->Get("hz")->MathParse<double>() : (max[2] - min_[2]) / n_[0];
      n_[2] = std::lround((max[2] - min_[2]) / h_[2]);
    }
    LOG_DBG(regular) << "PX: spacing h=" << h_.ToString() << " -> n=" << n_.ToString();
  }

  info_->Get("dim")->SetValue(dim_);
  info_->Get("region")->SetValue(region_);
  info_->Get("elements")->SetValue(n_.ToString());
  info_->Get("spacing")->SetValue(h_.ToString());
  info_->Get("min")->SetValue(min_.ToString());
}

void SimInputRegular::ReadMesh(Grid *basegrid)
{
  GridCFS* grid = dynamic_cast<GridCFS*>(basegrid);
  assert(grid != nullptr);

  if(dim_ != 0 && (dim_ != grid->GetDim()))
    throw Exception("Regular Mesh dimension does not match grid dimension");
  dim_ = grid->GetDim(); // in case it was auto, we set it to the grid dimension

  // parse box/elements/spacing now: only here is the global domain math parser ready
  // (RegularMesh is constructed in SetupIO, before 'domain' exists). Sets n_, h_, min_.
  ParseXML(xml_);

  base_ = grid->GetNumNodes(); // almost always = 0

  unsigned int nx = n_[0];
  unsigned int nxx = nx+1;
  unsigned int ny = n_[1];  
  unsigned int nyy = ny+1;
  unsigned int nz = dim_ == 3 ? n_[2] : 1;
  unsigned int nzz = nz+1;

  // reserve empty node space, overwritten by SetNodeCoordinate()
  unsigned int nnodes = nxx * nyy * (dim_ == 3 ? nzz : 1);
  grid->AddNodes(nnodes);
  Vector<double> loc(3);

  for(unsigned int z = 0; z < (dim_ == 3 ? nzz : 1); z++)
    for(unsigned int y = 0; y < nyy; y++)
      for(unsigned int x = 0; x < nxx; x++) // x is fastest 
      {
        unsigned int nodeNr = base_ + z*nyy*nxx + y*nxx + x + 1; // 1-based numbering
        loc[0] = min_[0] + x * h_[0];
        loc[1] = min_[1] + y * h_[1];
        loc[2] = min_[2] + z * h_[2];
        
        grid->SetNodeCoordinate(nodeNr, loc); // for GridCFS simply sets coords_
      }

  // master region - probably partly overwritten afterwards  - true for regular
  RegionIdType regid = grid->AddRegion(region_, true); // populate just the name
  Vector<unsigned int> conn(dim_ == 3 ? 8 : 4);  // connectivity
  unsigned int elemBase = grid->GetNumElems(); // almost always = 0
  grid->AddElems(nx*ny*nz); // reserve empty element space, overwritten by SetElemData()
  Elem* elem = nullptr;
  for(unsigned int z = 0; z < nz; z++)
  {
    for(unsigned int y = 0; y < ny; y++)
    {
      for(unsigned int x = 0; x < nx; x++) 
      {
        if(dim_ == 2) {
          conn[0] = base_ + y*nxx + x + 1;     // n1 = N(x,   y)
          conn[1] = base_ + y*nxx + x + 2;     // n2 = N(x+1, y)
          conn[2] = base_ + (y+1)*nxx + x + 2; // n3 = N(x+1, y+1)
          conn[3] = base_ + (y+1)*nxx + x + 1; // n4 = N(x,   y+1)
          elem =grid->SetElemData(elemBase + z*ny*nx + y*nx + x + 1, Elem::ET_QUAD4, regid, conn.GetPointer());
        }
        else {
          // left-bottom-back node number
          unsigned int ll = base_ + z*nyy*nxx + y*nxx + x + 1; // 1-based node number of left-bottom-back corner
          // offsets cf. Python: bottom face (z) CCW n1..n4, then top face (z+1) n5..n8
          conn[0] = ll;                     // n1 = N(x,   y,   z)
          conn[1] = ll + 1;                 // n2 = N(x+1, y,   z)
          conn[2] = ll + nxx + 1;           // n3 = N(x+1, y+1, z)
          conn[3] = ll + nxx;               // n4 = N(x,   y+1, z)
          conn[4] = ll + nxx*nyy;           // n5 = N(x,   y,   z+1)
          conn[5] = ll + nxx*nyy + 1;       // n6 = N(x+1, y,   z+1)
          conn[6] = ll + nxx*nyy + nxx + 1; // n7 = N(x+1, y+1, z+1)
          conn[7] = ll + nxx*nyy + nxx;     // n8 = N(x,   y+1, z+1)
          elem = grid->SetElemData(elemBase + z*ny*nx + y*nx + x + 1, Elem::ET_HEXA8, regid, conn.GetPointer());
        }
        assert(elem->extended != nullptr);
        elem->extended->barycenter.Set(min_[0] + (x+0.5) * h_[0], min_[1] + (y+0.5) * h_[1], min_[2] + (z+0.5) * h_[2]); // z is zero for 2D

      } // end x
    } // end y 
  } // end z

  // we define our own named nodes -> see mesh_tools.py
  if(dim_ == 2)
    Set2DNamedNodes(grid);
  else
    Set3DNamedNodes(grid);

  // read the optional region definitions which overwrite the master region - easy before FinalizeInit()
  for(auto pn : xml_->GetList("subRegion"))
    SetRegion(grid, pn);

  // we cannot call grid->CreateUserDefinedNodesElems() before FinishInit(). 
  // store the optional elemList and nodeList to be read by GridCFS::CreateUserDefinedNodesElems()
  grid->pendingDefinedNodesElemsLists.Push_back(xml_); // will search for elemList and nodeList
}

void SimInputRegular::SetRegion(GridCFS* grid, const PtrParamNode& pn)
{
  // the body is like defining named elements with any choice: bounds, list, condition ...
  // <region name="region_1">
  //   <bounds eps="0.001">
  //     <range comp="x" min=".8" max=".8" />
  //     <range comp="y" min=".4" max=".6" />
  //   </bounds>
  // </region>
  string name = pn->Get("name")->As<string>();
  RegionIdType regid = grid->GetRegionId(name, true); // silent for new region  
  if(regid == NO_REGION_ID) 
    regid = grid->AddRegion(name); // populate just the name

  StdVector<unsigned int> found_entities;

  // traverse all entities and either test by condition or check against given bounds
  if(pn->Has("test") || pn->Has("bounds"))
    grid->SelectByTraversal(name, false, pn->GetOneOf("test", "bounds"), found_entities);

  // check if node is defined by virtual coordinate
  if(pn->Has("coord"))
    grid->SelectByCoord(name, false, pn->Get("coord"), found_entities);

  // check if entity is defined by parametric virtual coordinates
  if(pn->Has("list") || pn->Has("expression"))
    grid->SelectBySearch(name, false, pn->GetOneOf("list", "expression"), found_entities);

  // overwrite the region for the found elements -> easy as Grid::FinishInit() is called only after ReadMesh()
  for(unsigned int elemNr : found_entities) {
    const Elem* old = grid->GetElem(elemNr);
    StdVector<unsigned int> connect = old->connect; // deep copy of the connectione
    grid->SetElemData(elemNr, old->type, regid, connect.GetPointer()); // overwrites existing element data and registers elemDim_
  }
    

  if(found_entities.GetSize() == 0)
    throw Exception("Region '" + name + "' in defined in regular mesh input, but no elements were found for it.");  

  auto spn = info_->Get("subRegion", ParamNode::APPEND);
  spn->Get("name")->SetValue(name);
  spn->Get("numElems")->SetValue(found_entities.GetSize());

  grid->regionData[regid].regular = true;
  grid->regionData[regid].barycenters = true; // set in ReadMesh()
}

void SimInputRegular::AddNamedNodes(Grid* grid, const std::string& name, unsigned int start, unsigned int end, unsigned int step)
{
  StdVector<unsigned int> range;
  if(end == 0) // TODO handle base
    range = Range<unsigned int>(base_ + start+1, base_ + start+2);
  else
    range = Range<unsigned int>(base_ + start+1, base_ + end+1, step);

  grid->AddNamedNodes(name, range);
}

void SimInputRegular::Set2DNamedNodes(Grid* grid)
{
  assert(dim_ == 2);

  unsigned int nx = n_[0];
  unsigned int ny = n_[1];  
  
  // based on create_2d_mesh() from mesh_tools.py 
  AddNamedNodes(grid, "bottom", 0, (nx+1));
  AddNamedNodes(grid, "top", (nx+1)*ny, (nx+1)*(ny+1));
  AddNamedNodes(grid, "left", 0, (nx+1)*ny+1, nx+1);
  AddNamedNodes(grid, "right", nx, (nx+1)*(ny+1), nx+1);
  AddNamedNodes(grid, "bottom_left", 0);
  AddNamedNodes(grid, "bottom_right", nx);
  AddNamedNodes(grid, "top_left", (nx+1)*ny);
  AddNamedNodes(grid, "top_right", (nx+1)*(ny+1)-1);
}

void SimInputRegular::Set3DNamedNodes(Grid* grid)
{
  assert(dim_ == 3);

  unsigned int nx = n_[0];
  unsigned int nxx = nx+1;
  unsigned int ny = n_[1];  
  unsigned int nyy = ny+1;
  unsigned int nz = n_[2];
  unsigned int nzz = nz+1;
  
  // based on create_3d_mesh() from mesh_tools.py 
  // faces
  // left (x=0) and right (x=nx): every nxx-th node
  AddNamedNodes(grid, "left",  0,  nxx*nyy*nz + nxx*ny + 1, nxx);
  AddNamedNodes(grid, "right", nx, nxx*nyy*nzz + 1,         nxx);
  // back (z=0) and front (z=nz) as they appear in ParaView: contiguous node blocks
  AddNamedNodes(grid, "back",  0,          nxx*nyy);
  AddNamedNodes(grid, "front", nz*nxx*nyy, nzz*nxx*nyy);
  // bottom (y=0) is a comb over z and x -> not a single range. top (y=ny) is the
  // same set shifted by the constant ny*nxx, so build bottom once and offset for top.
  StdVector<unsigned int> data; // bottom
  data.Reserve(nzz * nxx);
  for(unsigned int z = 0; z < nzz; z++)
    for(unsigned int x = 0; x < nxx; x++)
      data.Push_back(base_ +z*nyy*nxx + x + 1); // y=0, 1-based
  grid->AddNamedNodes("bottom", data);
  Vector<unsigned int> top(data.GetPointer(), data.GetSize(), WRAP);
  top.Add(ny*nxx); // shift the y=0 plane to y=ny directly on the data vector data (-> wrap!)
  grid->AddNamedNodes("top", data);

  // corners
  AddNamedNodes(grid, "left_bottom_back",   0);
  AddNamedNodes(grid, "right_bottom_back",  nx);
  AddNamedNodes(grid, "left_top_back",      nxx*ny);
  AddNamedNodes(grid, "right_top_back",     nxx*nyy - 1);
  AddNamedNodes(grid, "left_bottom_front",  nxx*nyy*nz);
  AddNamedNodes(grid, "right_bottom_front", nxx*nyy*nz + nx);
  AddNamedNodes(grid, "left_top_front",     nxx*nyy*nz + nxx*ny);
  AddNamedNodes(grid, "right_top_front",    nxx*nyy*nzz - 1);

  // edges
  AddNamedNodes(grid, "bottom_back",  0,                 nxx);
  AddNamedNodes(grid, "bottom_front", nxx*nyy*(nzz-1),   nxx*nyy*(nzz-1) + nxx);
  AddNamedNodes(grid, "bottom_left",  0,                 nxx*nyy*nzz - nxx - 1,   nxx*nyy);
  AddNamedNodes(grid, "bottom_right", nxx-1,             nxx*nyy*nzz - 2,         nxx*nyy);
  AddNamedNodes(grid, "top_back",     nxx*nyy - nxx,     nxx*nyy);
  AddNamedNodes(grid, "top_front",    nxx*nyy*nzz - nxx, nxx*nyy*nzz);
  AddNamedNodes(grid, "top_left",     nxx*nyy - nxx,     nxx*nyy*nzz,             nxx*nyy);
  AddNamedNodes(grid, "top_right",    nxx*nyy - 1,       nxx*nyy*nzz,             nxx*nyy);
  AddNamedNodes(grid, "back_left",    0,                 nxx*nyy - nxx + 1,       nxx);
  AddNamedNodes(grid, "back_right",   nxx-1,             nxx*nyy,                 nxx);
  AddNamedNodes(grid, "front_left",   nxx*nyy*nz,        nxx*nyy*nz + nxx*ny + 1, nxx);
  AddNamedNodes(grid, "front_right",  nxx*nyy*nz + nx,   nxx*nyy*nzz,             nxx);
}

} // namespace CoupledField
