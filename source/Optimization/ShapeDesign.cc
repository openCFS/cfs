#include "Optimization/ShapeDesign.hh"
#include "Utils/tools.hh"
#include "Domain/domain.hh"
#include "Domain/grid.hh"
#include "DataInOut/programOptions.hh"
#include "DataInOut/ParamHandling/Xerces.hh"

using namespace CoupledField;

ShapeDesign::ShapeDesign(StdVector<RegionIdType>& regionIds, ParamNodeList& design, ParamNodeList& transfer, ParamNodeList& result, ErsatzMaterial::Method method)
  : DesignSpace(regionIds, design, transfer, result, method)
{
  dim_ = domain->GetGrid()->GetDim();
  alsomatopt_ = method == ErsatzMaterial::SHAPE_PARAM_MAT;
}

void ShapeDesign::Configure(PtrParamNode pn, int objectives, int constraints){
  // done as in domain.cc Domain::ReadErsatzMaterial
  /* it would be nicer to use the schema to check, but this takes over 8 hours on a 3D example 
  std::string schema = progOpts->GetSchemaPathStr();
  schema += "/CFS-Simulation/Schemas/CFS_MeshDeformations.xsd";
  Xerces* xerces = new Xerces(pn->Get("meshdeformationfile")->As<std::string>(), schema); */
  Xerces* xerces = new Xerces(pn->Get("meshdeformationfile")->As<std::string>());
  PtrParamNode xml = xerces->CreateParamNodeInstance();
  delete xerces;
  PtrParamNode meshdefs = xml->Get("meshdeformations");
  nshapeparams_ = meshdefs->Get("parameters")->As<Integer>();
  ParamNodeList deformations = meshdefs->GetList("deformation");
  unsigned int numnodes = domain->GetGrid()->GetNumNodes();
  if((unsigned int)meshdefs->Get("nodes")->As<Integer>() != numnodes){
    EXCEPTION("The number of nodes given in the mesh deformation file for Shape Optimization does not correspond to the number of nodes in the grid!");
  }
  if((unsigned int)meshdefs->Get("dim")->As<Integer>() > dim_){ // one can specify fewer dimensions as in 3D usually only 2D movements are needed
    EXCEPTION("The dimension of the mesh deformation file is higher than the current problem!");
  }
  nodedeformations_.Reserve(numnodes);
  Matrix<double> empty; // an empty matrix
  for(unsigned int i = 0; i < numnodes; i++){
    nodedeformations_.Push_back(empty);
  }
  for(unsigned int i = 0; i < deformations.GetSize(); i++){
    PtrParamNode deformation = deformations[i];
    unsigned int node = deformation->Get("node")->As<Integer>();
    Matrix<double>& m = nodedeformations_[node-1]; // Nodes are numbered in the xml corresponding to the mesh, but we have 0-based index here
    if(m.GetNumRows() == 0){ // resize if first param for this node
      m.Resize(dim_, nshapeparams_);
      m.Init(); // this initializes with 0.0
    }
    unsigned int param = deformation->Get("param")->As<Integer>();
    unsigned int direction = deformation->Get("direction")->As<Integer>();
    double value = deformation->Get("value")->As<Double>();
    m[direction][param] = value;
  }
  // note that there exist empty matrices in the nodedeformations_
  shapeparams_.Reserve(nshapeparams_);
  double l = -1.0;
  double u = 1.0;
  scaling = 1.0;
  if(pn->Has("allShapeParams")){
    l = pn->Get("allShapeParams")->Get("lower")->As<Double>();
    u = pn->Get("allShapeParams")->Get("upper")->As<Double>();
    scaling = pn->Get("allShapeParams")->Get("scaling")->As<Double>();
  }
  for(unsigned int i = 0; i < nshapeparams_; i++){
    BaseDesignElement de;
    de.SetLowerBound(l);
    de.SetUpperBound(u);
    de.PostInit(objectives, constraints);
    shapeparams_.Push_back(de);
  }
  ParamNodeList shapeparams = pn->GetList("shapeParam");
  for(unsigned int i = 0; i < shapeparams.GetSize(); i++){
    unsigned int p = shapeparams[i]->Get("param")->As<Integer>();
    if(p >= nshapeparams_){
      throw Exception("ShapeParameter does  not exist!");
    }
    if(shapeparams[i]->Has("lower")){
      shapeparams_[p].SetLowerBound(shapeparams[i]->Get("lower")->As<Double>());
    }
    if(shapeparams[i]->Has("upper")){
      shapeparams_[p].SetUpperBound(shapeparams[i]->Get("upper")->As<Double>());
    }
  }
}

int ShapeDesign::ReadDesignFromExtern(const double* space_in){
  int old_design = design_id;
  bool new_design = false;
  for(unsigned int i=0; i < nshapeparams_; i++){
    double v = space_in[i] * scaling;
    if(!new_design && v != shapeparams_[i].GetDesign()){
      new_design = true;
    }
    shapeparams_[i].SetDesign(v);
  }
  if(new_design){
    UpdateCoordinates();
  }
  if(alsomatopt_){
    DesignSpace::ReadDesignFromExtern(space_in + nshapeparams_);
  }
  if(new_design && design_id <= old_design) design_id++; // if new design and not already changed by DesignSpace
  return(design_id);
}

void ShapeDesign::UpdateCoordinates(){
  Grid* grd = domain->GetGrid();
  int n = grd->GetNumNodes();
  for(int i = 0; i < n; i++){
    Matrix<double>& nodedefs = nodedeformations_[i]; // dim_ x nshapeparams_
    if(nodedefs.GetNumRows() > 0){ // if this node is dependent on any parameters at all
      Point p; // p = nodedefs * shapeparams
      for(unsigned int j = 0; j < dim_; j++) {
        double v = nodedefs[j][0] * shapeparams_[0].GetDesign();
        for(unsigned int k = 1; k < nshapeparams_; k++) {
          v += nodedefs[j][k] * shapeparams_[k].GetDesign();
        }
        p[j] = v;
      }
      grd->SetNodeOffset(i, p);
    }
  }
}

int ShapeDesign::WriteDesignToExtern(double* space_out, bool scaling) const {
  for(unsigned int i=0; i < nshapeparams_; i++){
    space_out[i] = shapeparams_[i].GetDesign() / scaling;
  }
  if(alsomatopt_){
    DesignSpace::WriteDesignToExtern(space_out + nshapeparams_, scaling);
  }
  return(design_id);
}

void ShapeDesign::WriteDenseGradientToExtern(StdVector<double>& out, DesignElement::ValueSpecifier vs, DesignElement::Access access, Condition* g, bool scaling) const
{
  // makes use of the window within out even  if only a part of the window is used in the alsomatopt_ case
  WriteShapeGradientToExtern(out, g);

  if(alsomatopt_)
  {
    // we call DesignSpace::WriteDenseGradientToExtern() for the ersatz material part.
    // therefore the window needs to be adjusted to a "subwindow" as WriteDenseGradientToExtern()
    // orientates itself on the given window.
    // The original window needs to be restored afterwards for assert() in BaseOptimization which
    // dies not know about separation into shape and mat variables
    StdVector<double>::Window org_window = out.window; // I like standard constructors :)
    // shift the window to do the rest
    assert(nshapeparams_ <= out.window.GetSize());
    // replace the original window by a subwindow excluding the shape stuff
    out.window.Set(out.window.GetStart() + nshapeparams_, out.window.GetSize() - nshapeparams_);
    DesignSpace::WriteDenseGradientToExtern(out, vs, access, g, scaling);

    // restore original window
    out.window = org_window;
  }
}

void ShapeDesign::WriteShapeGradientToExtern(StdVector<double>& out, Condition* g) const {
  unsigned int base = out.window.GetStart();
  assert(nshapeparams_ <= out.window.GetSize());
  for(unsigned int i=0; i < nshapeparams_; i++){
    out[base + i] = shapeparams_[i].GetGradient(NULL, g) * scaling;
  }
}

void ShapeDesign::Reset(DesignElement::ValueSpecifier vs, DesignElement::Type design){
  for(unsigned int i=0; i < nshapeparams_; i++)
    shapeparams_[i].Reset(vs);

  DesignSpace::Reset(vs, design);
}

void ShapeDesign::WriteBoundsToExtern(double* x_l, double* x_u) const {
  for(unsigned int i=0; i < nshapeparams_; i++){
    x_l[i] = shapeparams_[i].GetLowerBound() / scaling;
    x_u[i] = shapeparams_[i].GetUpperBound() / scaling;
  }
  if(alsomatopt_){
    DesignSpace::WriteBoundsToExtern(x_l + nshapeparams_, x_u + nshapeparams_);
  }
}

unsigned int ShapeDesign::GetNumberOfVariables() const {
  if(alsomatopt_){
    return(nshapeparams_ + DesignSpace::GetNumberOfVariables());
  }else{
    return(nshapeparams_);
  }
}

bool ShapeDesign::IsElemDependentAtAll(const StdVector<UInt>& connect){
  for(UInt k = 0; k < connect.GetSize(); ++k){
    if(nodedeformations_[connect[k]-1].GetNumRows() > 0){
      return(true);
    }
  }
  return(false);
}

bool ShapeDesign::GetElemNodesCoordDerivative(Matrix<Double> & coordMat, const StdVector<UInt> & connect, const int parameter){
  bool allIsZero = true;
  coordMat.Resize(dim_, connect.GetSize());
  for (UInt k=0; k < connect.GetSize(); k++) {
    Matrix<double>& nodedefs = nodedeformations_[connect[k]-1]; // complete deformation matrix of this node
    if(nodedefs.GetNumRows() > 0){ // if this node is dependent on any parameters at all
      // we have to extract the parameter-th column
      for (UInt actDim=0; actDim < dim_; actDim++) {
        double d = nodedefs[actDim][parameter];
        coordMat[actDim][k] = d;
        if(allIsZero && d != 0.0){ // this comparison is ok, not used nodes do really have the value 0.0 (it is never modified)
          allIsZero = false;
        }
      }
    }else{
      for(UInt actDim = 0; actDim < dim_; ++actDim){
        coordMat[actDim][k] = 0.0;
      }
    }
  }
  return(!allIsZero);
}

void ShapeDesign::AddShapeDerivatives(Objective* f, Condition* g, StdVector<double>& d, double weight){
  assert(d.GetSize() == nshapeparams_);
  for(unsigned int i = 0; i < nshapeparams_; i++){
    shapeparams_[i].AddGradient(f, g, d[i]*weight);
  }
}

