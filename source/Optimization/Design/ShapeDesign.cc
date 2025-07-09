#include <assert.h>
#include <ostream>
#include <string>

#include "DataInOut/Logging/LogConfigurator.hh"
#include "DataInOut/ParamHandling/XmlReader.hh"
#include "DataInOut/ProgramOptions.hh"
#include "Domain/Domain.hh"
#include "Domain/Mesh/Grid.hh"
#include "General/Exception.hh"
#include "MatVec/Matrix.hh"
#include "Optimization/Design/ShapeDesign.hh"
#include "Utils/Point.hh"

namespace CoupledField {

class Condition;
class Objective;

DEFINE_LOG(ShDes, "ShapeDesign")

ShapeDesign::ShapeDesign(StdVector<RegionIdType>& regions,  PtrParamNode pn, ErsatzMaterial::Method method)
  : AuxDesign(regions, pn, method)
{
  dim_ = domain->GetDim();
  export_fe_design_ = method == ErsatzMaterial::SHAPE_PARAM_MAT;
}

ShapeDesign::~ShapeDesign(){
  unsigned int n = nodedeformations_.GetSize(); 
  for(unsigned int i = 0; i < n; ++i){
    delete nodedeformations_[i];
  }
}

void ShapeDesign::Configure(PtrParamNode pn, int objectives, int constraints){
  // done as in domain.cc Domain::ReadErsatzMaterial
  // it would be nicer to use the schema to check, but this takes over 8 hours on a 3D example
  std::string schema = progOpts->GetSchemaPathStr();
  schema += "/CFS-Simulation/Schemas/CFS_MeshDeformations.xsd";
  PtrParamNode xml = XmlReader::ParseFile(pn->Get("meshdeformationfile")->As<std::string>());
  PtrParamNode meshdefs = xml->Get("meshdeformations");
  unsigned int nshapeparams = meshdefs->Get("parameters")->As<Integer>();
  ParamNodeList deformations = meshdefs->GetList("deformation");
  ParamNodeList constrs = meshdefs->GetList("constraint");
  unsigned int numnodes = domain->GetGrid()->GetNumNodes();
  if((unsigned int)meshdefs->Get("nodes")->As<Integer>() != numnodes){
    EXCEPTION("The number of nodes given in the mesh deformation file for Shape Optimization does not correspond to the number of nodes in the grid!");
  }
  if((unsigned int)meshdefs->Get("dim")->As<Integer>() > dim_){ // one can specify fewer dimensions as in 3D usually only 2D movements are needed
    EXCEPTION("The dimension of the mesh deformation file is higher than the current problem!");
  }
  nodedeformations_.Resize(numnodes, NULL);
  for(unsigned int i = 0; i < deformations.GetSize(); i++){
    PtrParamNode deformation = deformations[i];
    unsigned int node = deformation->Get("node")->As<Integer>();
    Matrix<double>* pm = nodedeformations_[node-1]; // Nodes are numbered in the xml corresponding to the mesh, but we have 0-based index here
    if(pm == NULL){ // resize if first param for this node
      pm = new Matrix<double>(dim_, nshapeparams);
      pm->Init(); // this initializes with 0.0
      nodedeformations_[node-1] = pm;
    }
    Matrix<double>& m = *pm;
    unsigned int param = deformation->Get("param")->As<Integer>();
    unsigned int direction = deformation->Get("direction")->As<Integer>();
    double value = deformation->Get("value")->As<Double>();
    m[direction][param] = value;
  }
  shapeconstraints_.Resize(constrs.GetSize());
  for(unsigned int i = 0; i < constrs.GetSize(); i++){
    ParamNodeList constr = constrs[i]->GetList("pair");
    ShapeConstraint& c = shapeconstraints_[i];
    assert(constr.GetSize() == 1); // more complex constraints could have several pairs, i.e. if not taking the inf-norm. This is not yet implemented
    PtrParamNode con = constr[0];
    c.param[0] = con->Get("p")->As<Integer>();
    c.factor[0] = con->Get("fp")->As<Double>();
    c.param[1] = con->Get("q")->As<Integer>();
    c.factor[1] = con->Get("fq")->As<Double>();
    assert(c.factor[1] == 0.0 || c.param[1] >= 0);
  }
  // note that there exist empty matrices in the nodedeformations_
  aux_design_.Reserve(nshapeparams);
  double l = -1.0;
  double u = 1.0;
  double v = 0.0;
  scaling_ = 1.0;
  if(pn->Has("allShapeParams")){
    l = pn->Get("allShapeParams")->Get("lower")->As<Double>();
    u = pn->Get("allShapeParams")->Get("upper")->As<Double>();
    v = pn->Get("allShapeParams")->Get("initial")->As<Double>();
    scaling_ = pn->Get("allShapeParams")->Get("scaling")->As<Double>();
  }
  int numberOfIndexesBefore = DesignSpace::GetNumberOfVariables();
  for(unsigned int i = 0; i < nshapeparams; i++){
    ShapeDesignElement de(numberOfIndexesBefore + i);
    de.SetLowerBound(l);
    de.SetUpperBound(u);
    de.SetDesign(v);
    de.PostInit(objectives, constraints);
    aux_design_.Push_back(de);
  }
  ParamNodeList shapeparams = pn->GetList("shapeParam");
  for(unsigned int i = 0; i < shapeparams.GetSize(); i++){
    unsigned int p = shapeparams[i]->Get("param")->As<Integer>();
    if(p >= nshapeparams){
      throw Exception("ShapeParameter does  not exist!");
    }
    if(shapeparams[i]->Has("lower")){
      aux_design_[p].SetLowerBound(shapeparams[i]->Get("lower")->As<Double>());
    }
    if(shapeparams[i]->Has("upper")){
      aux_design_[p].SetUpperBound(shapeparams[i]->Get("upper")->As<Double>());
    }
    if(shapeparams[i]->Has("initial")){
      aux_design_[p].SetDesign(shapeparams[i]->Get("initial")->As<Double>());
    }
  }
  UpdateCoordinates();
  
  assert(aux_design_.GetSize() == nshapeparams);
}

int ShapeDesign::ReadDesignFromExtern(const double* space_in, bool setAndWriteCurrent)
{
  int old_design = design_id;

  AuxDesign::ReadDesignFromExtern(space_in, setAndWriteCurrent); // sets design_id

  if(design_id != old_design)
    UpdateCoordinates();

  return design_id;
}

void ShapeDesign::UpdateCoordinates(){
  Grid* grd = domain->GetGrid();
  int n = grd->GetNumNodes();
  for(int i = 0; i < n; i++){
    Matrix<double>* pnodedefs = nodedeformations_[i];
    if(pnodedefs != NULL){ // if this node is dependent on any parameters at all
      Matrix<double>& nodedefs = *pnodedefs; // dim_ x nshapeparams_
      Point p; // p = nodedefs * shapeparams
      for(unsigned int j = 0; j < dim_; j++) {
        double v = nodedefs[j][0] * aux_design_[0].GetDesign();
        for(unsigned int k = 1; k < aux_design_.GetSize(); k++) {
          v += nodedefs[j][k] * aux_design_[k].GetDesign();
        }
        p[j] = v;
      }
      assert(false);
      // FIXME comment with fe-space merge grd->SetNodeOffset(i, p);
    }
  }
}



bool ShapeDesign::IsElemDependentAtAll(const StdVector<UInt>& connect){
  for(UInt k = 0; k < connect.GetSize(); ++k){
    if(nodedeformations_[connect[k]-1] != NULL){
      return(true);
    }
  }
  return(false);
}

bool ShapeDesign::GetElemNodesCoordDerivative(Matrix<Double>& coordMat, const StdVector<UInt>& connect, const int parameter){
  bool allIsZero = true;
  coordMat.Resize(dim_, connect.GetSize());
  for (UInt k=0; k < connect.GetSize(); k++) {
    Matrix<double>* pnodedefs = nodedeformations_[connect[k]-1];
    if(pnodedefs != NULL){ // if this node is dependent on any parameters at all
      Matrix<double>& nodedefs = *pnodedefs; // complete deformation matrix of this node
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

StdVector<ShapeDesign::ShapeConstraint>& ShapeDesign::GetShapeConstraints() {
  return(shapeconstraints_);
}


} // end of namespace
