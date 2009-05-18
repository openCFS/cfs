#include "Optimization/ShapeDesign.hh"
#include "Utils/tools.hh"
#include "Domain/domain.hh"
#include "Domain/grid.hh"
#include "DataInOut/programOptions.hh"
#include "DataInOut/ParamHandling/Xerces.hh"

using namespace CoupledField;

ShapeDesign::ShapeDesign(StdVector<RegionIdType>& regionIds, StdVector<ParamNode*>& design, StdVector<ParamNode*>& transfer, StdVector<ParamNode*>& result, ErsatzMaterial::Method method)
  : DesignSpace(regionIds, design, transfer, result, method)
{
  dim_ = domain->GetGrid()->GetDim();
  alsomatopt_ = method == ErsatzMaterial::SHAPE_PARAM_MAT;
}

void ShapeDesign::Configure(ParamNode* pn, int constraints){
  // done as in domain.cc Domain::ReadErsatzMaterial
  std::string schema = progOpts->GetSchemaPathStr();
  schema += "/CFS-Simulation/Schemas/CFS_MeshDeformations.xsd";
  Xerces* xerces = new Xerces(pn->Get("meshdeformationfile")->AsString(), schema);
  ParamNode* xml = xerces->CreateParamNodeInstance();
  delete xerces;
  ParamNode* meshdefs = xml->Get("meshdeformations");
  nshapeparams_ = meshdefs->Get("numberofparameters")->AsInt();
  StdVector<ParamNode*> deformations = meshdefs->GetList("deformation");
  unsigned int numnodes = domain->GetGrid()->GetNumNodes();
  nodedeformations_.Reserve(numnodes);
  for(unsigned int i = 0; i < numnodes; i++){
    Matrix<double> d(dim_, nshapeparams_);
    nodedeformations_.Push_back(d);
  }
  for(unsigned int i = 0; i < deformations.GetSize(); i++){
    ParamNode* deformation = deformations[i];
    unsigned int node = deformation->Get("node")->AsInt();
    unsigned int param = deformation->Get("param")->AsInt();
    unsigned int direction = deformation->Get("direction")->AsInt();
    double value = deformation->Get("value")->AsDouble();
    nodedeformations_[node-1][direction][param] = value; // Nodes are numbered in the xml corresponding to the mesh, but we have 0-based index here
  }
  shapeparams_.Reserve(nshapeparams_);
  double l = -1.0;
  double u = 1.0;
  scaling = 1.0;
  if(pn->Has("allShapeParams")){
    l = pn->Get("allShapeParams")->Get("lower")->AsDouble();
    u = pn->Get("allShapeParams")->Get("upper")->AsDouble();
    scaling = pn->Get("allShapeParams")->Get("scaling")->AsDouble();
  }
  for(unsigned int i = 0; i < nshapeparams_; i++){
    BaseDesignElement de;
    de.SetLowerBound(l);
    de.SetUpperBound(u);
    de.PostInit(constraints);
    shapeparams_.Push_back(de);
  }
  StdVector<ParamNode*> shapeparams = pn->GetList("shapeParam");
  for(unsigned int i = 0; i < shapeparams.GetSize(); i++){
    unsigned int p = shapeparams[i]->Get("param")->AsInt();
    if(p >= nshapeparams_){
      throw Exception("ShapeParameter does  not exist!");
    }
    if(shapeparams[i]->Has("lower")){
      shapeparams_[p].SetLowerBound(shapeparams[i]->Get("lower")->AsDouble());
    }
    if(shapeparams[i]->Has("upper")){
      shapeparams_[p].SetUpperBound(shapeparams[i]->Get("upper")->AsDouble());
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
    Point p; // p = nodedefs * shapeparams
    for(unsigned int j = 0; j < dim_; j++) {
      double *s = nodedefs[j];
      double v = s[0] * shapeparams_[0].GetDesign();
      for(unsigned int k = 1; k < nshapeparams_; k++) {
        v += s[k] * shapeparams_[k].GetDesign();
      }
      p[j] = v;
    }
    grd->SetNodeOffset(i, p);
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

void ShapeDesign::WriteGradientToExtern(double* out, DesignElement::ValueSpecifier vs, DesignElement::Access access, Condition* constraint, bool scaling) const {
  WriteShapeGradientToExtern(out, constraint);
  if(alsomatopt_){
    DesignSpace::WriteGradientToExtern(out + nshapeparams_, vs, access, constraint, scaling);
  }
}

void ShapeDesign::WriteShapeGradientToExtern(double* out, Condition* constraint) const {
  for(unsigned int i=0; i < nshapeparams_; i++){
    out[i] = shapeparams_[i].GetConstraintGradient(constraint) * scaling;
  }
}

void ShapeDesign::Reset(DesignElement::ValueSpecifier vs, DesignElement::Type design){
  if(design == DesignElement::DEFAULT){
    for(unsigned int i=0; i < nshapeparams_; i++){
      switch(vs){
      case DesignElement::DESIGN:
           shapeparams_[i].SetDesign(0.0);
           break;

      case DesignElement::COST_GRADIENT:
           shapeparams_[i].SetObjectiveGradient(0.0);
           break;

      default: throw Exception("value specifier not handled");
      }
    }
  }
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

void ShapeDesign::GetElemNodesCoordDerivative(Matrix<Double> & coordMat, const StdVector<UInt> & connect, const int parameter){
  coordMat.Resize(dim_, connect.GetSize());
  for (UInt k=0; k < connect.GetSize(); k++) {
    Matrix<double>& nodedefs = nodedeformations_[connect[k]-1]; // complete deformation matrix of this node
    // we have to extract the parameter-th column
    for (UInt actDim=0; actDim < dim_; actDim++) {
      coordMat[actDim][k] = nodedefs[actDim][parameter];
    }
  }
}

void ShapeDesign::SetShapeDerivatives(Condition* constraint, StdVector<double>& d){
  assert(d.GetSize() == nshapeparams_);
  for(unsigned int i = 0; i < nshapeparams_; i++){
    shapeparams_[i].SetConstraintGradient(constraint, d[i]);
  }
}

void ShapeDesign::AddShapeDerivatives(Condition* constraint, StdVector<double>& d, double weight){
  assert(d.GetSize() == nshapeparams_);
  for(unsigned int i = 0; i < nshapeparams_; i++){
    const double ov = constraint ? 0.0 : shapeparams_[i].GetConstraintGradient(constraint); // currently constraints are not calculated per excitation
    const double v = ov + d[i]*weight; 
    shapeparams_[i].SetConstraintGradient(constraint, v);
  }
}

