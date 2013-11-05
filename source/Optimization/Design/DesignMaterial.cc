#include <math.h>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>

#include "DataInOut/Logging/cfslog.hh"
#include "DataInOut/Logging/log.hpp"
#include "DesignMaterial.hh"
#include "Domain/domain.hh"
#include "Domain/grid.hh"
#include "General/defs.hh"
#include "General/exception.hh"
#include "Optimization/Design/DesignElement.hh"
#include "PDE/SinglePDE.hh"
#include "Utils/StdVector.hh"
#include "Elements/2D/quad9fe.hh"
#include "DataInOut/Logging/cfslog.hh"
#include "DataInOut/Logging/log.hpp"

#include "DataInOut/ParamHandling/ParamTools.hh"
#include "DataInOut/ParamHandling/Xerces.hh"
#include "MatVec/matrix.hh"

DECLARE_LOG(dm)
DEFINE_LOG(dm, "designMaterial")

using namespace CoupledField;
using std::string;

Enum<DesignMaterial::Type>         DesignMaterial::type;
Enum<DesignMaterial::TransIsoType> DesignMaterial::transIsoType;
Enum<DesignMaterial::Notation>     DesignMaterial::notation;

DesignMaterial::DesignMaterial(PtrParamNode pn, OptimizationMaterial::System material, StdVector<DesignID>& design)
{
  type_ = type.Parse(pn->Get("type")->As<string>());
  
  dim = domain->GetGrid()->GetDim();
  
  transIsoType_ = transIsoType.Parse(pn->Get("isoplane")->As<string>());
  
  massIsDesign_ = pn->Get("optimizeMass")->As<bool>();
  
  massFactor_ = pn->Get("massFactor")->As<double>();
  
  penalty_ = pn->Get("penalty")->As<double>();
  
  trace_ = pn->Get("trace")->As<double>();
  
  dampingIsDesign_ = pn->Get("optimizeDamping")->As<bool>();

  deb_ = true;

  // collect all designs here, to check whether all are given
  unsigned int r = RequiredParameters(material);
  StdVector<DesignElement::Type> d;
  d.Reserve(r);
  // copy the ones from DesignSpace
  for(unsigned int i=0; i < design.GetSize(); ++i){
    d.Push_back(design[i].design);
  }
  // read non-design parameters
  ParamNodeList params = pn->GetList("param");
  for(unsigned int i=0; i < params.GetSize(); i++){
    DesignElement::Type dt = DesignElement::type.Parse(params[i]->Get("name")->As<string>());
    SetParameter(dt, params[i]->Get("value")->As<Double>());
    if(d.Find(dt) < 0){
      d.Push_back(dt);
    }
  }
  if(!CheckRequiredDesigns(d)){
    throw Exception("Not all Parameters for chosen DesignMaterial given");
  }else if(design.GetSize() > r){ // design.GetSize() < r is impossible as CheckRequiredDesigns passed
    info->Get("optimization/header/designSpace")->Get(ParamNode::WARNING)->SetValue("There are designs specified that are not used!");
  }

  if(type_ == HOM_RECT)
  {
    PtrParamNode hr = pn->Get("homRect");
    hom_rect_samples_.Resize(9, 6);
    FillHomRectSamples(hr, 0, "0.0", "0.0");
    FillHomRectSamples(hr, 1, "0.5", "0.0");
    FillHomRectSamples(hr, 2, "0.5", "0.5");
    FillHomRectSamples(hr, 3, "0.0", "0.5");
    FillHomRectSamples(hr, 4, "0.25", "0.0");
    FillHomRectSamples(hr, 5, "0.5", "0.25");
    FillHomRectSamples(hr, 6, "0.25", "0.5");
    FillHomRectSamples(hr, 7, "0.0", "0.25");
    FillHomRectSamples(hr, 8, "0.25", "0.25");
  }
  if(type_ == HOM_RECT_C1) {
    if (dim == 2) {
      PtrParamNode hr = pn->Get("homRectC1");
      std::string file = hr->Get("file")->As<std::string>();
      Xerces xerces(file);
      PtrParamNode root = xerces.CreateParamNodeInstance();
      int dim1 = root->Get("coeff11/matrix/dim1")->As<int>();
      int dim2 = root->Get("coeff11/matrix/dim2")->As<int>();
      int dim3 = root->Get("a/matrix/dim1")->As<int>();
      int dim4 = root->Get("b/matrix/dim1")->As<int>();
      ParamTools::AsTensor<double>(root->Get("coeff11/matrix/real"),dim1, dim2, hom_rect_coeff11_);
      ParamTools::AsTensor<double>(root->Get("coeff12/matrix/real"),dim1, dim2, hom_rect_coeff12_);
      ParamTools::AsTensor<double>(root->Get("coeff22/matrix/real"),dim1, dim2, hom_rect_coeff22_);
      ParamTools::AsTensor<double>(root->Get("coeff33/matrix/real"),dim1, dim2, hom_rect_coeff33_);
      ParamTools::AsTensor<double>(root->Get("a/matrix/real"),dim3, 1, hom_rect_a_);
      ParamTools::AsTensor<double>(root->Get("b/matrix/real"),dim4, 1, hom_rect_b_);
      Notation notation = root->Get("notation")->As<string>() == "voigt" ? VOIGT : HILL_MANDEL;
      hom_rect_coeff33_ = hom_rect_coeff33_ * (notation == VOIGT ? 2.0 : 1.0);
    } else if (dim == 3) {
      PtrParamNode hr = pn->Get("homRectC1");
      std::string file = hr->Get("file")->As<std::string>();
      Xerces xerces(file);
      PtrParamNode root = xerces.CreateParamNodeInstance();
      int dim1 = root->Get("coeff11/matrix/dim1")->As<int>();
      int dim2 = root->Get("coeff11/matrix/dim2")->As<int>();
      int dim3 = root->Get("a/matrix/dim1")->As<int>();
      int dim4 = root->Get("b/matrix/dim1")->As<int>();
      int dim5 = root->Get("c/matrix/dim1")->As<int>();
      ParamTools::AsTensor<double>(root->Get("coeff11/matrix/real"),dim1, dim2, hom_rect_coeff11_);
      ParamTools::AsTensor<double>(root->Get("coeff12/matrix/real"),dim1, dim2, hom_rect_coeff12_);
      ParamTools::AsTensor<double>(root->Get("coeff22/matrix/real"),dim1, dim2, hom_rect_coeff22_);
      ParamTools::AsTensor<double>(root->Get("coeff33/matrix/real"),dim1, dim2, hom_rect_coeff33_);
      ParamTools::AsTensor<double>(root->Get("coeff13/matrix/real"),dim1, dim2, hom_rect_coeff13_);
      ParamTools::AsTensor<double>(root->Get("coeff23/matrix/real"),dim1, dim2, hom_rect_coeff23_);
      ParamTools::AsTensor<double>(root->Get("coeff44/matrix/real"),dim1, dim2, hom_rect_coeff44_);
      ParamTools::AsTensor<double>(root->Get("coeff55/matrix/real"),dim1, dim2, hom_rect_coeff55_);
      ParamTools::AsTensor<double>(root->Get("coeff66/matrix/real"),dim1, dim2, hom_rect_coeff66_);

      ParamTools::AsTensor<double>(root->Get("a/matrix/real"),dim3, 1, hom_rect_a_);
      ParamTools::AsTensor<double>(root->Get("b/matrix/real"),dim4, 1, hom_rect_b_);
      ParamTools::AsTensor<double>(root->Get("c/matrix/real"),dim5, 1, hom_rect_c_);
      // the internal tensor representation in hom_rect_samples_ is HILL-MANDEL!
      Notation notation = root->Get("notation")->As<string>() == "voigt" ? VOIGT : HILL_MANDEL;
      hom_rect_coeff44_ = hom_rect_coeff44_ * (notation == VOIGT ? 2.0 : 1.0);
      hom_rect_coeff55_ = hom_rect_coeff55_ * (notation == VOIGT ? 2.0 : 1.0);
      hom_rect_coeff66_ = hom_rect_coeff66_ * (notation == VOIGT ? 2.0 : 1.0);
      // the tensor is orthotropic
    }
  }
}

void DesignMaterial::FillHomRectSamples(PtrParamNode homRect, unsigned int idx, const string& a, const string& b)
{
  // the internal tensor representation in hom_rect_samples_ is HILL-MANDEL!
  Notation notation = homRect->Get("notation")->As<string>() == "voigt" ? VOIGT : HILL_MANDEL;

  PtrParamNode data = homRect->GetByVal("data", "a", a, "b", b);
  hom_rect_samples_[idx][DesignElement::TENSOR11 - DesignElement::TENSOR11] = data->Get("e11")->As<double>();
  hom_rect_samples_[idx][DesignElement::TENSOR12 - DesignElement::TENSOR11] = data->Get("e12")->As<double>();
  hom_rect_samples_[idx][DesignElement::TENSOR22 - DesignElement::TENSOR11] = data->Get("e22")->As<double>();
  hom_rect_samples_[idx][DesignElement::TENSOR33 - DesignElement::TENSOR11] = data->Get("e33")->As<double>() * (notation == VOIGT ? 2.0 : 1.0);
  hom_rect_samples_[idx][DesignElement::TENSOR13 - DesignElement::TENSOR11] = 0.0;
  hom_rect_samples_[idx][DesignElement::TENSOR23 - DesignElement::TENSOR11] = 0.0;
}

unsigned int DesignMaterial::RequiredParameters(OptimizationMaterial::System material)
{
  unsigned int r = MassIsDesign() ? 1 : 0;
  if(DampingIsDesign()){
    r += 2;
  }
  switch(type_){
  case FMO:
    assert(material == OptimizationMaterial::MECH || material == OptimizationMaterial::PIEZOCOUPLING);
    return r + (material == OptimizationMaterial::MECH ? 6 : 15);
  case ISOTROPIC:
  case LAME_ISOTROPIC:
    return r+2;
  case TRANSVERSAL_ISOTROPIC:
  case TRANSVERSAL_ISOTROPIC_BOXED:
    return domain->GetSinglePDE("mechanic")->GetSubTensorType() == PLANE_STRESS ? r+4 : r+5;
  case DENSITY_TIMES_TRANSVERSAL_ISOTROPIC:
  case DENSITY_TIMES_TRANSVERSAL_ISOTROPIC_BOXED:
  case LAMINATES:
    return r+5;
  case HOM_RECT:
    return r+3;
  case HOM_RECT_C1:
    return r+3;
  case DENSITY_TIMES_2D_TENSOR_CONSTANT_TRACE:
  case DENSITY_TIMES_ROT_TRANSVERSAL_ISOTROPIC_BOXED:
    return r+6;
  case DENSITY_TIMES_2D_TENSOR:
  case DENSITY_TIMES_ROTATED_2D_TENSOR:
    return r+7;
  }

  assert(false);
  return 0;
}

bool DesignMaterial::CheckRequiredDesigns(StdVector<DesignElement::Type>& design){
  if(MassIsDesign() && design.Find(DesignElement::MASS) < 0){
    return(false);
  }
  if(DampingIsDesign() 
      && (design.Find(DesignElement::DAMPINGALPHA) < 0 || design.Find(DesignElement::DAMPINGBETA) < 0) ){ 
    return(false);
  }
  switch(type_){
  case FMO:
    return(design.Find(DesignElement::TENSOR11) >= 0
        && design.Find(DesignElement::TENSOR22) >= 0
        && design.Find(DesignElement::TENSOR33) >= 0
        && design.Find(DesignElement::TENSOR23) >= 0
        && design.Find(DesignElement::TENSOR13) >= 0
        && design.Find(DesignElement::TENSOR12) >= 0);
  case ISOTROPIC:
    return(design.Find(DesignElement::EMODUL) >=0
        && design.Find(DesignElement::POISSON) >= 0);
  case LAME_ISOTROPIC:
    return(design.Find(DesignElement::LAMELAMBDA) >= 0 
        && design.Find(DesignElement::LAMEMU) >= 0);
  case TRANSVERSAL_ISOTROPIC:
  case TRANSVERSAL_ISOTROPIC_BOXED:
    return(design.Find(DesignElement::EMODULISO) >= 0
        && design.Find(DesignElement::EMODUL) >= 0 
        && design.Find(DesignElement::POISSON) >= 0 
        && design.Find(DesignElement::GMODUL) >= 0
        && domain->GetSinglePDE("mechanic")->GetSubTensorType() == PLANE_STRESS ? true : design.Find(DesignElement::POISSONISO) >= 0);
  case DENSITY_TIMES_TRANSVERSAL_ISOTROPIC:
  case DENSITY_TIMES_TRANSVERSAL_ISOTROPIC_BOXED:
    return(design.Find(DesignElement::DENSITY) >= 0
        && design.Find(DesignElement::EMODULISO) >= 0
        && design.Find(DesignElement::EMODUL) >= 0
        && design.Find(DesignElement::POISSON) >= 0
        && design.Find(DesignElement::GMODUL) >= 0);
  case DENSITY_TIMES_ROT_TRANSVERSAL_ISOTROPIC_BOXED:
    return(design.Find(DesignElement::DENSITY) >= 0
        && design.Find(DesignElement::EMODULISO) >= 0
        && design.Find(DesignElement::EMODUL) >= 0
        && design.Find(DesignElement::POISSON) >= 0
        && design.Find(DesignElement::GMODUL) >= 0
        && design.Find(DesignElement::ROTANGLE));
  case DENSITY_TIMES_2D_TENSOR:
    return(design.Find(DesignElement::DENSITY) >= 0
        && design.Find(DesignElement::TENSOR11) >= 0
        && design.Find(DesignElement::TENSOR22) >= 0
        && design.Find(DesignElement::TENSOR33) >= 0
        && design.Find(DesignElement::TENSOR23) >= 0
        && design.Find(DesignElement::TENSOR13) >= 0
        && design.Find(DesignElement::TENSOR12) >= 0);
  case DENSITY_TIMES_2D_TENSOR_CONSTANT_TRACE:
    return(design.Find(DesignElement::DENSITY) >= 0
        && design.Find(DesignElement::TENSOR11) >= 0
        && design.Find(DesignElement::TENSOR22) >= 0
        && design.Find(DesignElement::TENSOR23) >= 0
        && design.Find(DesignElement::TENSOR13) >= 0
        && design.Find(DesignElement::TENSOR12) >= 0);
  case DENSITY_TIMES_ROTATED_2D_TENSOR:
    return(design.Find(DesignElement::DENSITY) >= 0
        && design.Find(DesignElement::TENSOR11) >= 0
        && design.Find(DesignElement::TENSOR33) >= 0
        && design.Find(DesignElement::TENSOR23) >= 0
        && design.Find(DesignElement::TENSOR13) >= 0
        && design.Find(DesignElement::TENSOR12) >= 0
        && design.Find(DesignElement::ROTANGLE) >= 0);
  case LAMINATES:
    return(design.Find(DesignElement::STIFF1) >= 0
        && design.Find(DesignElement::STIFF2) >= 0
        && design.Find(DesignElement::ROTANGLE) >= 0
        && design.Find(DesignElement::EMODUL) >= 0
        && design.Find(DesignElement::POISSON) >= 0);
  case HOM_RECT:
    return(design.Find(DesignElement::STIFF1) >= 0
        && design.Find(DesignElement::STIFF2) >= 0
        && design.Find(DesignElement::ROTANGLE) >= 0);
  case HOM_RECT_C1:
    if (dim == 3) {
    return(design.Find(DesignElement::STIFF1) >= 0
           && design.Find(DesignElement::STIFF2) >= 0
           && design.Find(DesignElement::STIFF3) >= 0
           && design.Find(DesignElement::ROTX) >= 0
           && design.Find(DesignElement::ROTY) >= 0
           && design.Find(DesignElement::ROTZ) >= 0);
    } else {
      return(design.Find(DesignElement::STIFF1) >= 0
             && design.Find(DesignElement::STIFF2) >= 0
             && design.Find(DesignElement::ROTANGLE) >= 0);
    }
  }
  assert(false);
  return false;
}

void DesignMaterial::SetParameter(const DesignElement::Type p, const double value){
  params_[p] = value;
  // LOG_DBG3(dm) << "SP p=" << DesignElement::type.ToString(p) << " v=" << value;
}


void DesignMaterial::GetIsoMaterialTensor(Matrix<double>& t, SubTensorType subTensor, DesignElement::Type direction){
  double E = params_[DesignElement::EMODUL];
  double nu = params_[DesignElement::POISSON];
  switch(direction){
  case DesignElement::NO_DERIVATIVE:
  {
    double lambda = nu * E / ((1.0+nu)*(1.0-2.0*nu));
    double mu = E / (2.0*(1.0+nu));
    double diag = lambda + 2.0*mu;
    SetIsoTensor(t, subTensor, diag, lambda, mu);
    break;
  }
  case DesignElement::EMODUL:
  {
    double dlambda_dE = nu/((1.0+nu)*(1.0-2.0*nu));
    double dmu_dE = 1.0/(2.0*(1.0+nu));
    double ddiag_dE = dlambda_dE + 2.0*dmu_dE;
    SetIsoTensor(t, subTensor, ddiag_dE, dlambda_dE, dmu_dE);
    break;
  }
  case DesignElement::POISSON:
  {
    double dlambda_dnu = (1.0+2.0*nu*nu)*E / ((1.0+nu)*(1.0+nu)*(1.0-2.0*nu)*(1.0-2.0*nu));
    double dmu_dnu = E / (-2.0*(1.0+nu)*(1.0+nu));
    double ddiag_dnu = dlambda_dnu + 2.0*dmu_dnu;
    SetIsoTensor(t, subTensor, ddiag_dnu, dlambda_dnu, dmu_dnu);
    break;
  }
  default:
    ZeroTensor(t, subTensor); // any derivative in any direction other than EMODUL or POISSON is zero
    break;
  }  
}

double DesignMaterial::GetIsoMaterialMass(DesignElement::Type direction){
  double E = params_[DesignElement::EMODUL];
  double nu = params_[DesignElement::POISSON];
  switch(direction){
  case DesignElement::NO_DERIVATIVE:
  {
    double lambda = nu * E / ((1.0+nu)*(1.0-2.0*nu));
    double mu = E / (2.0*(1.0+nu));
    double diag = lambda + 2.0*mu;
    return(GetIsoMass(diag, mu));
  }
  case DesignElement::EMODUL:
  {
    double dlambda_dE = nu/((1.0+nu)*(1.0-2.0*nu));
    double dmu_dE = 1.0/(2.0*(1.0+nu));
    double ddiag_dE = dlambda_dE + 2.0*dmu_dE;
    return(GetIsoMass(ddiag_dE, dmu_dE));
  }
  case DesignElement::POISSON:
  {
    double dlambda_dnu = (1.0+2.0*nu*nu)*E / ((1.0+nu)*(1.0+nu)*(1.0-2.0*nu)*(1.0-2.0*nu));
    double dmu_dnu = E / (-2.0*(1.0+nu)*(1.0+nu));
    double ddiag_dnu = dlambda_dnu + 2.0*dmu_dnu;
    return(GetIsoMass(ddiag_dnu, dmu_dnu));
  }
  default:
    return(0.0); // any derivative in any direction other than EMODUL or POISSON is zero
  }  
}

void DesignMaterial::GetLameMaterialTensor(Matrix<double>& t, SubTensorType subTensor, DesignElement::Type direction){
  switch(direction){
  case DesignElement::NO_DERIVATIVE:
  {
    double lambda = params_[DesignElement::LAMELAMBDA];
    double mu = params_[DesignElement::LAMEMU];
    double diag = lambda + 2.0*mu;
    SetIsoTensor(t, subTensor, diag, lambda, mu);
    break;
  }
  case DesignElement::LAMELAMBDA:
    SetIsoTensor(t, subTensor, 1.0, 1.0, 0.0);
    break;
  case DesignElement::LAMEMU:
    SetIsoTensor(t, subTensor, 2.0, 0.0, 1.0);
    break;
  default:
    ZeroTensor(t, subTensor); // any derivative in any direction other than EMODUL or POISSON is zero
    break;
  }
}

double DesignMaterial::GetLameMaterialMass(DesignElement::Type direction){
  switch(direction){
  case DesignElement::NO_DERIVATIVE:
  {
    double lambda = params_[DesignElement::LAMELAMBDA];
    double mu = params_[DesignElement::LAMEMU];
    double diag = lambda + 2.0*mu;
    return(GetIsoMass(diag, mu));
  }
  case DesignElement::LAMELAMBDA:
    return(GetIsoMass(1.0, 0.0));
  case DesignElement::LAMEMU:
    return(GetIsoMass(2.0, 1.0));
  default:
    return(0.0); // any derivative in any direction other than EMODUL or POISSON is zero
  }
}

void DesignMaterial::GetTransIsoMaterialTensor(Matrix<double>& t, SubTensorType subTensor, DesignElement::Type direction){
  double E = params_[DesignElement::EMODULISO];
  double E3 = params_[DesignElement::EMODUL];
  double G3 = params_[DesignElement::GMODUL];
  if(type_ == DENSITY_TIMES_ROT_TRANSVERSAL_ISOTROPIC_BOXED)
    G3 = 2*G3;
  double nu13 = params_[DesignElement::POISSON]; //used as theta in the boxed formulations

  if(subTensor == PLANE_STRESS){ //This is only implemented for density times tensor formulations yet
    double dens(1.0), ninv2(0.0), D(0.0), D3(0.0), nD3(0.0);
    if(type_ == DENSITY_TIMES_TRANSVERSAL_ISOTROPIC || type_ == DENSITY_TIMES_TRANSVERSAL_ISOTROPIC_BOXED
        || type_ == DENSITY_TIMES_ROT_TRANSVERSAL_ISOTROPIC_BOXED)
    {
      dens = params_[DesignElement::DENSITY];
      dens = std::pow(dens, penalty_);
    }
    if(type_ == TRANSVERSAL_ISOTROPIC || type_ == DENSITY_TIMES_TRANSVERSAL_ISOTROPIC)
    {
      ninv2 = E3-nu13*nu13*E;
      //assert(ninv2>0.0); //positivity of the elasticity tensor is violated. Use constraint "parametrized-plane-stress-pos-def" > 0
      ninv2 = 1/(ninv2*ninv2);
    }
    switch(direction){
    case DesignElement::NO_DERIVATIVE:
    case DesignElement::ROTANGLE:
    {
      if(type_ == TRANSVERSAL_ISOTROPIC || type_ == DENSITY_TIMES_TRANSVERSAL_ISOTROPIC)
      {
        D = E*E3/(E3-nu13*nu13*E);
        D3 = E3*E3/(E3-nu13*nu13*E);
        nD3 = nu13*E*E3/(E3-nu13*nu13*E);
      } else
      {
        D = E/(1-nu13);
        D3 = E3/(1-nu13);
        nD3 = sqrt(E*E3*nu13)/(1-nu13);
      }
      SetTransIsoTensor(t, subTensor, dens*D, 0, 0, dens*D3, dens*nD3, dens*G3);
      break;
    }
    case DesignElement::DENSITY:
    {
      if(type_ == DENSITY_TIMES_TRANSVERSAL_ISOTROPIC || type_ == DENSITY_TIMES_TRANSVERSAL_ISOTROPIC_BOXED
          || type_ == DENSITY_TIMES_ROT_TRANSVERSAL_ISOTROPIC_BOXED)
      {
        dens = params_[DesignElement::DENSITY];
        if(penalty_ == 1.0){
          dens = 1.0;
        }else dens = penalty_*std::pow(dens, penalty_-1);
      }
      if(type_ == TRANSVERSAL_ISOTROPIC || type_ == DENSITY_TIMES_TRANSVERSAL_ISOTROPIC)
      {
        D = E*E3/(E3-nu13*nu13*E);
        D3 = E3*E3/(E3-nu13*nu13*E);
        nD3 = nu13*E*E3/(E3-nu13*nu13*E);
      }else
      {
        D = E/(1-nu13);
        D3 = E3/(1-nu13);
        nD3 = sqrt(E*E3*nu13)/(1-nu13);
      }
      SetTransIsoTensor(t, subTensor, dens*D, 0, 0, dens*D3, dens*nD3, dens*G3);
      break;
    }
    case DesignElement::EMODULISO:
    {
      if(type_ == TRANSVERSAL_ISOTROPIC || type_ == DENSITY_TIMES_TRANSVERSAL_ISOTROPIC)
      {
        D = E3*E3*ninv2;
        D3 = nu13*nu13*E3*E3*ninv2;
        nD3 = nu13*E3*E3*ninv2;
      } else
      {
        D = 1/(1-nu13);
        nD3 = sqrt(E3*nu13/E)/(2-2*nu13);
      }
      SetTransIsoTensor(t, subTensor, dens*D, 0, 0, dens*D3, dens*nD3, 0);
      break;
    }
    case DesignElement::EMODUL:
    {
      if(type_ == TRANSVERSAL_ISOTROPIC || type_ == DENSITY_TIMES_TRANSVERSAL_ISOTROPIC)
      {
        D = -E*E*nu13*nu13*ninv2;
        D3 = -E3*(-E3+2*nu13*nu13*E)*ninv2;
        nD3 = -nu13*nu13*nu13*E*E*ninv2;
      }else
      {
        D3 = 1/(1-nu13);
        nD3 = sqrt(E*nu13/E3)/(2-2*nu13);
      }
      SetTransIsoTensor(t, subTensor, dens*D, 0, 0, dens*D3, dens*nD3, 0);
      break;
    }
    case DesignElement::POISSON:
    {
      if(type_ == TRANSVERSAL_ISOTROPIC || type_ == DENSITY_TIMES_TRANSVERSAL_ISOTROPIC)
      {
        D = 2*nu13*E*E*E3*ninv2;
        D3 = 2*nu13*E*E3*E3*ninv2;
        nD3 = E*E3*(nu13*nu13*E+E3)*ninv2;
      }else
      {
        D = 1/((1-nu13)*(1-nu13));
        D3 = E3*D;
        nD3 = sqrt(E*E3/nu13)*(nu13+1)*D*0.5;
        D = E*D;
      }
      SetTransIsoTensor(t, subTensor, dens*D, 0, 0, dens*D3, dens*nD3, 0);
      break;
    }
    case DesignElement::GMODUL:
    {
      SetTransIsoTensor(t, subTensor, 0, 0, 0, 0, 0, dens);
      break;
    }
    default:
      ZeroTensor(t, subTensor);
      return;
    } // switch direction
    if(type_ == DENSITY_TIMES_ROTATED_2D_TENSOR || type_ == DENSITY_TIMES_ROT_TRANSVERSAL_ISOTROPIC_BOXED){
      double rotAngle = params_[DesignElement::ROTANGLE];
      RotateHMStiffnessTensor(t, subTensor, direction, rotAngle);
      //    static int count(0);
      //    if (count % 10 == 0 && count/100 % 10 == 0){
      ////      std::cout << "(" << (count/100 % 10)*(count % 10)+1 << ")" << t.ToString() << std::endl;
      //      std::cout << t(0,0) << " " << t(0,1) << " " << t(1,1) << " " << t(0,2) << " " << t(1,2) << " " << t(2,2) << std::endl;
      //    }
      //    count++;
    }
    return;
  } // PLANE_STRESS

  double nu = params_[DesignElement::POISSONISO];
  double nu3;
  double n3;
  double c;
  if(type_ == TRANSVERSAL_ISOTROPIC){
    nu3 = nu13 * E3/E;
    n3 = nu3*nu3*E/E3;
    c = (1.0-nu-2.0*n3); // this is the interesting thing, this must not get 0, however this would imply a volume (trace of tensor) of infinity, so it is hopefully not occuring
    if(c < 1e-8) {
      c = 1e-8;
    }
  }else{
    nu3 = sqrt(0.5*(1.0-nu)*E3/E)*nu13;
    n3 = nu3*nu3*E/E3;
    c = (1.0-nu-2.0*n3);
  }
  double f = E/((1.0+nu)*c);
  double dE = 0.0, dE3 = 0.0, dnu = 0.0, dnu3 = 0.0, dn3 = 0.0, dG3 = 0.0;
  
  switch(direction){
  case DesignElement::NO_DERIVATIVE:
  {
    double D = (1.0-n3)*f;
    double D3 = (1.0-nu)*E3/c;
    double nD = (nu+n3)*f;
    double nD3 = (1.0+nu)*nu3*f;
    double G = 0.5*E/(1.0+nu);
    SetTransIsoTensor(t, subTensor, D, nD, G, D3, nD3, G3);
    return;
  }
  case DesignElement::DENSITY:
  {
    throw Exception("direction DENSITY not implemented yet");
  }
  case DesignElement::EMODULISO:
    dE = 1.0;
    if(type_ == TRANSVERSAL_ISOTROPIC){
      dnu3 = -E3*nu13/(E*E);
      dn3 = nu3/E3 * (2.0*E*dnu3 + nu3);
    }else{
      dnu3 = -sqrt(0.125*(1.0-nu)*E3/E)*nu13/E;
    }
    break;
  case DesignElement::EMODUL:
    dE3 = 1.0;
    if(type_ == TRANSVERSAL_ISOTROPIC){
      dnu3 = nu13/E;
      dn3 = nu3*E/E3 * (2.0*dnu3 - nu3/E3);
    }else{
      dnu3 = sqrt(0.125*(1.0-nu)/(E*E3))*nu13;
    }
    break;
  case DesignElement::POISSONISO:
    dnu = 1.0;
    if(type_ == TRANSVERSAL_ISOTROPIC_BOXED){ // else = 0
      dnu3 = -sqrt(0.125*E3/(E*(1.0-nu)))*nu13;
      dn3 = -0.5*nu13*nu13;
    }
  break;
  case DesignElement::POISSON:
    if(type_ == TRANSVERSAL_ISOTROPIC){
      dnu3 = 1.0;
      dn3 = 2.0*nu3*E/E3*dnu3;
    }else{
      dnu3 = sqrt(0.5*(1.0-nu)*E3/E);
      dn3 = (1.0-nu)*nu13;
    }
    break;
  case DesignElement::GMODUL:
    dG3 = 1.0;
    break;
  default:
    ZeroTensor(t, subTensor);
    return;
  } // switch(direction)
  double dc = -dnu-2.0*dn3;
  double df = ( dE - E*dnu/(1.0+nu) - E*dc/c ) / ((1.0+nu)*c);
  double dD = (1.0-n3)*df - dn3*f;
  double dnD = (nu+n3)*df + (dnu+dn3)*f;
  double dD3 = ( (1.0-nu)*dE3 - dnu*E3 - (1.0-nu)*E3*dc/c ) / c;
  double dnD3 = (1.0+nu)*nu3*df + (1.0+nu)*dnu3*f + dnu*nu3*f;
  double dG = 0.5 * ( (1.0+nu)*dE - E*dnu ) / ( (1.0+nu)*(1.0+nu) );
  SetTransIsoTensor(t, subTensor, dD, dnD, dG, dD3, dnD3, dG3);
}
 
double DesignMaterial::GetTransIsoMaterialMass(DesignElement::Type direction){
  double E = params_[DesignElement::EMODULISO];
  double E3 = params_[DesignElement::EMODUL];
  double nu = params_[DesignElement::POISSONISO];
  double nu13 = params_[DesignElement::POISSON];
  double nu3 = nu13 * E3/E;
  double n3 = nu3*nu3*E/E3;
  double c = (1.0-nu-2.0*n3); // this is the interesting thing, this must not get 0, however this would imply a volume (trace of tensor) of infinity, so it is hopefully not occuring
  if(type_ == TRANSVERSAL_ISOTROPIC){
    if(c < 1e-8) {
      c = 1e-8;
    }
  }else{
    nu3 = sqrt(0.5*(1.0-nu)*E3/E)*nu13;
    n3 = nu3*nu3*E/E3;
    c = (1.0-nu-2.0*n3);
  }
  double f = E/((1.0+nu)*c);
  double dE = 0.0, dE3 = 0.0, dnu = 0.0, dnu3 = 0.0, dn3 = 0.0, dG3 = 0.0;
  
  switch(direction){
  case DesignElement::NO_DERIVATIVE:
  {
    double D = (1.0-n3)*f;
    double D3 = (1.0-nu)*E3/c;
    double G3 = params_[DesignElement::GMODUL];
    double G = 0.5*E/(1.0+nu);
    return(GetTransIsoMass(D, G, D3, G3));
  }
  case DesignElement::EMODULISO:
    dE = 1.0;
    if(type_ == TRANSVERSAL_ISOTROPIC){
      dnu3 = -E3*nu13/(E*E);
      dn3 = nu3/E3 * (2.0*E*dnu3 + nu3);
    }else{
      dnu3 = -sqrt(0.125*(1.0-nu)*E3/E)*nu13/E;
    }
    break;
  case DesignElement::EMODUL:
    dE3 = 1.0;
    if(type_ == TRANSVERSAL_ISOTROPIC){
      dnu3 = nu13/E;
      dn3 = nu3*E/E3 * (2.0*dnu3 - nu3/E3);
    }else{
      dnu3 = sqrt(0.125*(1.0-nu)/(E*E3))*nu13;
    }
    break;
  case DesignElement::POISSONISO:
    dnu = 1.0;
    if(type_ == TRANSVERSAL_ISOTROPIC_BOXED){ // else = 0
      dnu3 = -sqrt(0.125*E3/(E*(1.0-nu)))*nu13;
      dn3 = -0.5*nu13*nu13;
    }
  break;
  case DesignElement::POISSON:
    if(type_ == TRANSVERSAL_ISOTROPIC){
      dnu3 = 1.0;
      dn3 = 2.0*nu3*E/E3*dnu3;
    }else{
      dnu3 = sqrt(0.5*(1.0-nu)*E3/E);
      dn3 = (1.0-nu)*nu13;
    }
    break;
  case DesignElement::GMODUL:
    dG3 = 1.0;
    break;
  default:
    return(0.0);
  } // switch(direction)
  double dc = -dnu-2.0*dn3;
  double df = ( dE - E*dnu/(1.0+nu) - E*dc/c ) / ((1.0+nu)*c);
  double dD = (1.0-n3)*df - dn3*f;
  double dD3 = ( (1.0-nu)*dE3 - dnu*E3 - (1.0-nu)*E3*dc/c ) / c;
  double dG = 0.5 * ( (1.0+nu)*dE - E*dnu ) / ( (1.0+nu)*(1.0+nu) );
  return(GetTransIsoMass(dD, dG, dD3, dG3));
}

void DesignMaterial::GetDensityTimes2dTensorTensor(Matrix<double>& t, SubTensorType subTensor, DesignElement::Type direction)
{
  // DumpParams();
  double e11 = 0;
  double e22 = 0;
  double e33 = 0;
  double e23 = 0;
  double e13 = 0;
  double e12 = 0;
  if(direction == DesignElement::NO_DERIVATIVE || direction == DesignElement::DENSITY || direction == DesignElement::ROTANGLE){
   e11 = params_[DesignElement::TENSOR11];
   if(type_ == DENSITY_TIMES_2D_TENSOR_CONSTANT_TRACE){
     e22 = params_[DesignElement::TENSOR22];
     e33 = 0.5 * (trace_ - e11 - e22); 
   }else if(type_ == DENSITY_TIMES_ROTATED_2D_TENSOR){
     e22 = 15 - e11;
     e11 += 1.0;
     e33 = params_[DesignElement::TENSOR33];
   }else{
     e22 = params_[DesignElement::TENSOR22];
     e33 = params_[DesignElement::TENSOR33];
   }
   e23 = params_[DesignElement::TENSOR23];
   e13 = params_[DesignElement::TENSOR13];
   e12 = params_[DesignElement::TENSOR12];
  }
  double d = params_[DesignElement::DENSITY];
  switch(direction){
  case DesignElement::NO_DERIVATIVE:
  case DesignElement::ROTANGLE:
    d = pow(d, penalty_);
    Set2dVoigtTensor(t, d*e11, d*e22, d*e33, d*e23, d*e13, d*e12);
    break;
  case DesignElement::DENSITY:
    if(penalty_ == 1.0){
      d = 1.0;
    }else{
      d = penalty_*pow(d, penalty_-1);
    }
    Set2dVoigtTensor(t, d*e11, d*e22, d*e33, d*e23, d*e13, d*e12);
    break;
  case DesignElement::TENSOR11:
    d = pow(d, penalty_);
    Set2dVoigtTensor(t, d, type_ == DENSITY_TIMES_ROTATED_2D_TENSOR ? -d : 0.0, type_ == DENSITY_TIMES_2D_TENSOR_CONSTANT_TRACE ? -0.5*d : 0.0, 0.0, 0.0, 0.0);
    break;
  case DesignElement::TENSOR22:
    d = pow(d, penalty_);
    Set2dVoigtTensor(t, 0.0, d, type_ == DENSITY_TIMES_2D_TENSOR_CONSTANT_TRACE ? -0.5*d : 0.0, 0.0, 0.0, 0.0);
    break;
  case DesignElement::TENSOR33:
    d = pow(d, penalty_);
    Set2dVoigtTensor(t, 0.0, 0.0, d, 0.0, 0.0, 0.0);
    break;
  case DesignElement::TENSOR23:
    d = pow(d, penalty_);
    Set2dVoigtTensor(t, 0.0, 0.0, 0.0, d, 0.0, 0.0);
    break;
  case DesignElement::TENSOR13:
    d = pow(d, penalty_);
    Set2dVoigtTensor(t, 0.0, 0.0, 0.0, 0.0, d, 0.0);
    break;
  case DesignElement::TENSOR12:
    d = pow(d, penalty_);
    Set2dVoigtTensor(t, 0.0, 0.0, 0.0, 0.0, 0.0, d);
    break;
  default:
    ZeroTensor(t, subTensor);
    break;
  }
  if(type_ == DENSITY_TIMES_ROTATED_2D_TENSOR){
    double rotAngle = params_[DesignElement::ROTANGLE];
    RotateHMStiffnessTensor(t, subTensor, direction, rotAngle);
//    static int count(0);
//    if (count % 10 == 0 && count/100 % 10 == 0){
////      std::cout << "(" << (count/100 % 10)*(count % 10)+1 << ")" << t.ToString() << std::endl;
//      std::cout << t(0,0) << " " << t(0,1) << " " << t(1,1) << " " << t(0,2) << " " << t(1,2) << " " << t(2,2) << std::endl;
//    }
//    count++;
    return;
  }
}

void DesignMaterial::GetAnisotropicTensor(Matrix<double>& t, DesignElement::Type direction, Notation notation)
{
  // We use the anisotropic tensor only for solving FMO problems. Then we assume the design to be in Hill-Mandel
  // notation and therefore we need to transform it for using it in CFS
  double e11 = 0;
  double e22 = 0;
  double e33 = 0;
  double e23 = 0;
  double e13 = 0;
  double e12 = 0;
  assert(direction != DesignElement::DENSITY);
  if(direction == DesignElement::NO_DERIVATIVE)
  {
    e11 = params_[DesignElement::TENSOR11];
    e22 = params_[DesignElement::TENSOR22];
    e33 = params_[DesignElement::TENSOR33] * (notation == VOIGT ? 0.5 : 1.0);
    e23 = params_[DesignElement::TENSOR23] * (notation == VOIGT ? 1.0/sqrt(2.0) : 1.0);
    e13 = params_[DesignElement::TENSOR13] * (notation == VOIGT ? 1.0/sqrt(2.0) : 1.0);
    e12 = params_[DesignElement::TENSOR12];
  }
  switch(direction)
  {
  case DesignElement::NO_DERIVATIVE:
    Set2dVoigtTensor(t, e11, e22, e33, e23, e13, e12);
    break;
  case DesignElement::TENSOR11:
    Set2dVoigtTensor(t, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    break;
  case DesignElement::TENSOR22:
    Set2dVoigtTensor(t, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0);
    break;
  case DesignElement::TENSOR33:
    Set2dVoigtTensor(t, 0.0, 0.0, notation == VOIGT ? 0.5 : 1.0, 0.0, 0.0, 0.0);
    break;
  case DesignElement::TENSOR23:
    Set2dVoigtTensor(t, 0.0, 0.0, 0.0, notation == VOIGT ? 1.0/sqrt(2.0) : 1.0, 0.0, 0.0);
    break;
  case DesignElement::TENSOR13:
    Set2dVoigtTensor(t, 0.0, 0.0, 0.0, 0.0, notation == VOIGT ? 1.0/sqrt(2.0) : 1.0, 0.0);
    break;
  case DesignElement::TENSOR12:
    Set2dVoigtTensor(t, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0);
    break;
  default:
    // for piezo FMO the derivative w.r.t. dielec_11, ... is zero
    ZeroTensor(t, PLANE_STRAIN);
    return;
  }
}

void DesignMaterial::GetHomRectTensor(Matrix<double>& E, SubTensorType subTensor, DesignElement::Type direction, Notation notation)
{
   // only relevant for hom_rect
   Quad9FE fe;

   double a = params_[DesignElement::STIFF1];
   double b = params_[DesignElement::STIFF2];
   double c = subTensor == FULL ? params_[DesignElement::STIFF3] : 0.0;
   double alpha = 0.;
   double beta = 0.;
   double gamma = 0.;
   double rotAngle = 0.;
   if (subTensor == FULL) {
     alpha = params_[DesignElement::ROTX];
     beta = params_[DesignElement::ROTY];
     gamma = params_[DesignElement::ROTZ];
   } else {
     rotAngle = params_[DesignElement::ROTANGLE];
   }

   Vector<double> p(subTensor == FULL ? 3 : 2);



   if (type_== HOM_RECT) {
     p[0] = -1.0 + 4 * a; // assume max 0.5
     p[1] = -1.0 + 4 * b; // assume max 0.5
   }
   if (type_ == HOM_RECT_C1) {
     p[0] = a;
     p[1] = b;
     if (subTensor == FULL) {
       p[2] = c;
     }
   }
   #ifndef NDEBUG
     Vector<double> peps(p);
     double eps = 1e-8;
     Matrix<double> Eeps(E);
     Matrix<double> Etmp(E);
     Vector<double> Diff(subTensor == FULL ? 3 : 2);
   #endif


   LOG_DBG2(dm) << "GHRT: dir=" << (direction == DesignElement::NO_DERIVATIVE ? "no_derivative" : DesignElement::type.ToString(direction))
                << " not=" << notation << " rotAngle=" << rotAngle << " a=" << a << " b=" << b <<" c="<<(subTensor == FULL ? c : 0.0)<< " -> " << p.ToString();

   switch(direction)
   {
   case DesignElement::NO_DERIVATIVE:
   case DesignElement::ROTANGLE:
   {
     if(type_ == HOM_RECT) {
       Vector<double> shape;
       fe.GetShFnc(shape, p, NULL);
       ApplyHomRectTensor(E, shape);
       LOG_DBG2(dm) << "GHRT: shape=" << shape.ToString();
     }

     if(type_ == HOM_RECT_C1) {
       ApplyHomRectC1Tensor(E,p,direction, subTensor);
     }


     break;
   }

   case DesignElement::STIFF1:
   case DesignElement::STIFF2:
   case DesignElement::STIFF3:
   {
     if(type_ == HOM_RECT) {
     Matrix<double> jac;
     Matrix<double> dummy; // not used -> strange function ?! :(
     fe.GetLocDerivShFnc(jac, p, dummy, NULL);
     LOG_DBG3(dm) << "GHRT: jac=" << jac.ToString(2);

     Vector<double> d_shape;
     jac.GetCol(d_shape, direction == DesignElement::STIFF1 ? 0 : 1); // a or by

     ApplyHomRectTensor(E, d_shape);

     // correct scaling to local FE coordinates
     E *= 4;

     LOG_DBG2(dm) << "GHRT: d_shape=" << d_shape.ToString();
     }
     if(type_ == HOM_RECT_C1) {
       ApplyHomRectC1Tensor(E, p,direction,subTensor);
       if (subTensor == FULL) {
               #ifndef NDEBUG
                 if (direction == DesignElement::STIFF1) {
                   peps[0] += eps;
                 } else if (direction == DesignElement::STIFF2) {
                   peps[1] += eps;
                 } else if (direction == DesignElement::STIFF3) {
                   peps[2] += eps;
                 }
                 deb_ = false;
                 ApplyHomRectC1Tensor(Eeps,peps,DesignElement::NO_DERIVATIVE,subTensor);
                 ApplyHomRectC1Tensor(Etmp,p,DesignElement::NO_DERIVATIVE,subTensor);
                 deb_ = true;
                 LOG_DBG(dm)<<"Eeps11: "<<std::setprecision(10)<<Eeps[0][0]<<", E11: "<<Etmp[0][0]<<" Diff: "<<Eeps[0][0]-Etmp[0][0];
                 LOG_DBG(dm)<<"Eeps13: "<<std::setprecision(10)<<Eeps[0][2]<<", E13: "<<Etmp[0][2]<<" Diff: "<<Eeps[0][2]-Etmp[0][2];
                 double e11 = (Eeps[0][0]-Etmp[0][0])/eps;
                 double e12 = (Eeps[0][1]-Etmp[0][1])/eps;
                 double e13 = (Eeps[0][2]-Etmp[0][2])/eps;
                 double e22 = (Eeps[1][1]-Etmp[1][1])/eps;
                 double e23 = (Eeps[1][2]-Etmp[1][2])/eps;
                 double e33 = (Eeps[2][2]-Etmp[2][2])/eps;
                 double e44 = (Eeps[3][3]-Etmp[3][3])/eps;
                 double e55 = (Eeps[4][4]-Etmp[4][4])/eps;
                 double e66 = (Eeps[5][5]-Etmp[5][5])/eps;
                 LOG_DBG(dm)<<"FD Derivative "<<((direction == DesignElement::STIFF1)?"1":((direction == DesignElement::STIFF2) ? "2":"3"))<<" E11= "<<e11<<" E12= "<<e12<<" E22= "<< e22<<
                 " E33= "<<e33<<" E23= "<<e23<<" E13= "<<e13<<" E44= "<<e44<<" E55= "<<e55<<" E66= "<<e66;
                 LOG_DBG(dm)<<"deriv p= "<<p[0]<<", "<<p[1]<<", "<<p[2];
                 LOG_DBG(dm)<<"FD Derivative - Derivative: "<<((direction == DesignElement::STIFF1)?"1":((direction == DesignElement::STIFF2) ? "2":"3"))<<" diff E11= "<<E[0][0]-e11<<" diff E12= "<<E[0][1]-e12<<" diff E22= "<< E[1][1]-e22<<
                                  " diff E33= "<<E[2][2]-e33<<" diff E23= "<<E[1][2]-e23<<" diff E13= "<<E[0][2]-e13<<" diff E44= "<<E[3][3]-e44<<" diff E55= "<<E[4][4]-e55<<" diff E66= "<<E[5][5]-e66;
               #endif
             }
     }
     break;
   }

   default:
     ZeroTensor(E, subTensor == FULL ? FULL : PLANE);
   }
   LOG_DBG2(dm) << "GHRT: E before rotation = " << E.ToString(2);
   if (subTensor == FULL) {
     RotateHMStiffnessTensor(E, FULL, direction, alpha,notation,beta, gamma);
   }else {
     RotateHMStiffnessTensor(E, PLANE, direction, rotAngle, notation);
   }
   LOG_DBG2(dm) << "GHRT: E after rotation =  " << E.ToString(2);


/*   for(double y = 0; y <= 0.5; y += 0.25)
   {
     for(double x = 0; x <= 0.5; x += 0.25 )
     {
       p[0] = -1.0 + 4 * x; //
       p[1] = -1.0 + 4 * y;
       fe.GetShFnc(shape, p, NULL);
       hom_rect_samples_.GetCol(data, DesignElement::TENSOR11 - DesignElement::TENSOR11);
       assert(shape.GetSize() == data.GetSize());

       double val = shape * data;
       std::cout << "x=" << x << " y=" << y << " xi=" << p[0] << " eta=" << p[1] << " -> " << val << std::endl; //" s=" << shape.ToString() << " d=" << data.ToString() << std::endl;
     }
   }
*/
}

void DesignMaterial::ApplyHomRectTensor(Matrix<double>& E, const Vector<double>& shape) const
{
  E.Resize(3,3);
  E.Init(); // for off-diagonal

  Vector<double> data;

  hom_rect_samples_.GetCol(data, DesignElement::TENSOR11 - DesignElement::TENSOR11);
  E[1-1][1-1] = shape * data;
  LOG_DBG(dm) << "AHRT 11=" << E[1-1][1-1] << " data=" << data.ToString();

  hom_rect_samples_.GetCol(data, DesignElement::TENSOR12 - DesignElement::TENSOR11);
  E[1-1][2-1] = shape * data;
  E[2-1][1-1] = E[1-1][2-1];
  LOG_DBG(dm) << "AHRT 12=" << E[1-1][2-1] << " data=" << data.ToString();

  hom_rect_samples_.GetCol(data, DesignElement::TENSOR22 - DesignElement::TENSOR11);
  E[2-1][2-1] = shape * data;
  LOG_DBG(dm) << "AHRT 22=" << E[2-1][2-1] << " data=" << data.ToString();

  hom_rect_samples_.GetCol(data, DesignElement::TENSOR33 - DesignElement::TENSOR11);
  E[3-1][3-1] = shape * data;
  LOG_DBG(dm) << "AHRT 33=" << E[3-1][3-1] << " data=" << data.ToString();


}
void DesignMaterial::ApplyHomRectC1Tensor(Matrix<double>& E, Vector<double>& p,DesignElement::Type direction, SubTensorType subTensor) const
{
  PtrParamNode inf_warn = info->Get("optimization/header/designSpace");
  int m = hom_rect_a_.GetNumRows();
  int n = hom_rect_b_.GetNumRows();
  int o = (subTensor == FULL) ? hom_rect_c_.GetNumRows() : 0;
  double da = hom_rect_a_[1][0] - hom_rect_a_[0][0];
  double db = hom_rect_b_[1][0] - hom_rect_b_[0][0];
  double dc = (subTensor == FULL) ? hom_rect_c_[1][0]-hom_rect_c_[0][0] : 1;
  int j = -1;
  for (int i=0;i<m-1;i++) {
    if (hom_rect_a_[i][0] <= p[0] && p[0] < hom_rect_a_[i+1][0]) {
      j=i;
      break;
    } else if (p[0] == hom_rect_a_[m-1][0]) {
      j=m-2;
      break;
    }else if (p[0] > hom_rect_a_[m-1][0]){
      j=m-2;
      p[0] = 1.;
      if (p[0]>1.01) {
        inf_warn->Get(ParamNode::WARNING)->SetValue("Interpolation of Hom_RectC1 tensor failed. Design Variable p[0]" +lexical_cast<string>(p[0])+ " out of bounds ");
      }
      break;
    } else if (p[0] < 0.) {
      j=0;
      p[0] = 0.;
      if (p[0]<-0.01) {
        inf_warn->Get(ParamNode::WARNING)->SetValue("Interpolation of Hom_RectC1 tensor failed. Design Variable p[1]" +lexical_cast<string>(p[0])+ " out of bounds ");
      }
      break;
    }
  }
  assert(j!=-1);
  int k = -1;
  for (int i=0;i<n-1;i++) {
    if (hom_rect_b_[i][0] <= p[1] && p[1] < hom_rect_b_[i+1][0]) {
      k=i;
      break;
    } else if (p[1] == hom_rect_b_[n-1][0]) {
      k=n-2;
      break;
    } else if (p[1] > hom_rect_b_[n-1][0]){
      k=n-2;
      p[1] = 1.;
      if (p[1]>1.01) {
        inf_warn->Get(ParamNode::WARNING)->SetValue("Interpolation of Hom_RectC1 tensor failed. Design Variable p[1]" +lexical_cast<string>(p[1])+ " out of bounds ");
      }
      break;
    } else if (p[1] < 0.) {
      k=0;
      p[1] = 0.;
      if (p[1]<-0.01) {
        inf_warn->Get(ParamNode::WARNING)->SetValue("Interpolation of Hom_RectC1 tensor failed. Design Variable p[1]" +lexical_cast<string>(p[1])+ " out of bounds ");
      }
      break;
    }
  }
  assert(k != -1);
  if (subTensor == FULL) {
    int l = -1;
    E.Resize(6,6);
    E.Init(); // for off-diagonal
     for (int i=0;i<o-1;i++) {
       if (hom_rect_c_[i][0] <= p[2] && p[2] < hom_rect_c_[i+1][0]) {
         l=i;
         break;
       } else if (p[2] == hom_rect_c_[o-1][0]) {
         l=o-2;
         break;
       } else if (p[2] > hom_rect_c_[o-1][0]){
         l=o-2;
         p[2] = 1.;
         if (p[2]>1.01) {
           inf_warn->Get(ParamNode::WARNING)->SetValue("Interpolation of Hom_RectC1 tensor failed. Design Variable p[2]" +lexical_cast<string>(p[2])+ " out of bounds ");
         }
         break;
       } else if (p[2] < 0.) {
         l=0;
         p[2] = 0.;
         if (p[2]<-0.01) {
           inf_warn->Get(ParamNode::WARNING)->SetValue("Interpolation of Hom_RectC1 tensor failed. Design Variable p[2]" +lexical_cast<string>(p[2])+ " out of bounds ");
         }
         break;
       }

     }
     assert(l != -1);
     if (direction == DesignElement::NO_DERIVATIVE || direction == DesignElement::ROTANGLE) {
       E[1-1][1-1] = EvaluateC1Interpolation_3D(E, p, hom_rect_coeff11_, da,db,dc,j,k,l,m,n,o);
       E[1-1][2-1] = EvaluateC1Interpolation_3D(E, p, hom_rect_coeff12_, da,db,dc,j,k,l,m,n,o);
       E[1-1][3-1] = EvaluateC1Interpolation_3D(E, p, hom_rect_coeff13_, da,db,dc,j,k,l,m,n,o);
       E[2-1][3-1] = EvaluateC1Interpolation_3D(E, p, hom_rect_coeff23_, da,db,dc,j,k,l,m,n,o);
       E[2-1][2-1] = EvaluateC1Interpolation_3D(E, p, hom_rect_coeff22_, da,db,dc,j,k,l,m,n,o);
       E[3-1][3-1] = EvaluateC1Interpolation_3D(E, p, hom_rect_coeff33_, da,db,dc,j,k,l,m,n,o);
       E[4-1][4-1] = EvaluateC1Interpolation_3D(E, p, hom_rect_coeff44_, da,db,dc,j,k,l,m,n,o);
       E[5-1][5-1] = EvaluateC1Interpolation_3D(E, p, hom_rect_coeff55_, da,db,dc,j,k,l,m,n,o);
       E[6-1][6-1] = EvaluateC1Interpolation_3D(E, p, hom_rect_coeff66_, da,db,dc,j,k,l,m,n,o);
       E[2-1][1-1] = E[1-1][2-1];
       E[3-1][1-1] = E[1-1][3-1];
       E[3-1][2-1] = E[2-1][3-1];
       if (deb_) {
         LOG_DBG(dm)<<"E11= "<<E[0][0]<<" E12= "<<E[0][1]<<" E22= "<< E[1][1]<<" E33= "<<E[2][2]<<" E23= "<<E[1][2]<<" E13= "<<E[0][2]<<" E44= "<<E[3][3]<<" E55= "<<E[4][4]<<" E66= "<<E[5][5];
       }
     } else {
         E[1-1][1-1] = EvaluateC1Interpolation_Deriv_3D(E, p, hom_rect_coeff11_, da,db,dc,j,k,l,m,n,o,direction);
         E[1-1][2-1] = EvaluateC1Interpolation_Deriv_3D(E, p, hom_rect_coeff12_, da,db,dc,j,k,l,m,n,o,direction);
         E[1-1][3-1] = EvaluateC1Interpolation_Deriv_3D(E, p, hom_rect_coeff13_, da,db,dc,j,k,l,m,n,o,direction);
         E[2-1][3-1] = EvaluateC1Interpolation_Deriv_3D(E, p, hom_rect_coeff23_, da,db,dc,j,k,l,m,n,o,direction);
         E[2-1][2-1] = EvaluateC1Interpolation_Deriv_3D(E, p, hom_rect_coeff22_, da,db,dc,j,k,l,m,n,o,direction);
         E[3-1][3-1] = EvaluateC1Interpolation_Deriv_3D(E, p, hom_rect_coeff33_, da,db,dc,j,k,l,m,n,o,direction);
         E[4-1][4-1] = EvaluateC1Interpolation_Deriv_3D(E, p, hom_rect_coeff44_, da,db,dc,j,k,l,m,n,o,direction);
         E[5-1][5-1] = EvaluateC1Interpolation_Deriv_3D(E, p, hom_rect_coeff55_, da,db,dc,j,k,l,m,n,o,direction);
         E[6-1][6-1] = EvaluateC1Interpolation_Deriv_3D(E, p, hom_rect_coeff66_, da,db,dc,j,k,l,m,n,o,direction);
         E[2-1][1-1] = E[1-1][2-1];
         E[3-1][1-1] = E[1-1][3-1];
         E[3-1][2-1] = E[2-1][3-1];
         if (deb_) {
           LOG_DBG(dm)<<"Derivative "<<((direction == DesignElement::STIFF1)?"1":((direction == DesignElement::STIFF2) ? "2":"3"))<<" E11= "<<E[0][0]<<" E12= "<<E[0][1]<<" E22= "<< E[1][1]<<" E33= "<<E[2][2]<<" E23= "<<E[1][2]<<" E13= "<<E[0][2]<<" E44= "<<E[3][3]<<" E55= "<<E[4][4]<<" E66= "<<E[5][5];
         }
     }
  } else {
    E.Resize(3,3);
    E.Init(); // for off-diagonal
    if (direction == DesignElement::NO_DERIVATIVE || direction == DesignElement::ROTANGLE) {
      E[1-1][1-1] = EvaluateC1Interpolation(E, p, hom_rect_coeff11_,da,db,j,k,m,n);
      E[1-1][2-1] = EvaluateC1Interpolation(E, p, hom_rect_coeff12_,da,db ,j,k,m,n);
      E[2-1][1-1] = E[1-1][2-1];
      E[2-1][2-1] = EvaluateC1Interpolation(E, p, hom_rect_coeff22_, da,db,j,k,m,n);
      E[3-1][3-1] = EvaluateC1Interpolation(E, p, hom_rect_coeff33_, da,db,j,k,m,n);
      LOG_DBG(dm)<<"E11= "<<E[0][0]<<" E12= "<<E[0][1]<<" E22= "<< E[1][1]<<" E33= "<<E[2][2];
      LOG_DBG(dm)<<"hom_rect_coeff11 = "<<hom_rect_coeff11_[0][0];
    } else {
        E[1-1][1-1] = EvaluateC1Interpolation_Deriv(E, p, hom_rect_coeff11_, da,db,j,k,m,n,direction);
        E[1-1][2-1] = EvaluateC1Interpolation_Deriv(E, p, hom_rect_coeff12_, da,db,j,k,m,n,direction);
        E[2-1][1-1] = E[1-1][2-1];
        E[2-1][2-1] = EvaluateC1Interpolation_Deriv(E, p, hom_rect_coeff22_, da,db,j,k,m,n,direction);
        E[3-1][3-1] = EvaluateC1Interpolation_Deriv(E, p, hom_rect_coeff33_, da,db,j,k,m,n,direction);
        LOG_DBG(dm)<<"Derivative "<<((direction == DesignElement::STIFF1)?"1":"2")<<" E11= "<<E[0][0]<<" E12= "<<E[0][1]<<" E22= "<< E[1][1]<<" E33= "<<E[2][2];
    }
  }

}

double DesignMaterial::EvaluateC1Interpolation_3D(Matrix<double>& E, Vector<double>& p,const Matrix<double> & coeff, double & da,double & db, double & dc, int & j, int & k, int & l,int & m,int & n, int &o) const{
    LOG_DBG(dm) <<"p=["<<p[0]<<","<<p[1]<<", "<<p[2]<<"]";
    double t=(p[0]-hom_rect_a_[j][0])/da;
    double u =(p[1]-hom_rect_b_[k][0])/db;
    double v=(p[2]-hom_rect_c_[l][0])/dc;
    if (deb_) {
      LOG_DBG(dm)<<"u = "<<u<<" t= "<<t<<" v= "<<v;
    }
    double res = 0;
    for (int ii = 0;ii<4;ii++) {
      for (int jj=0;jj<4;jj++) {
        for (int kk=0;kk<4;kk++) {
          res += coeff[(n-1)*(o-1)*j+(o-1)*k+l][ii+4*jj+16*kk]*pow(t,ii)*pow(u,jj)*pow(v,kk);
        }
      }
    }
    if (deb_) {
      LOG_DBG(dm) << "Result =" << res;
    }
    return res;
}

double DesignMaterial::EvaluateC1Interpolation_Deriv_3D(Matrix<double>& E,  Vector<double>& p,const Matrix<double> & coeff, double & da,double & db,double & dc,int & j, int & k, int & l,int & m,int & n, int & o, DesignElement::Type direction) const{
    double u =(p[0]-hom_rect_a_[j][0])/(da);
    double t=(p[1]-hom_rect_b_[k][0])/(db);
    double v = (p[2]-hom_rect_c_[l][0])/dc;
    if (deb_) {
      LOG_DBG(dm)<<"Deriv: u = "<<u<<" t= "<<t<<" v= "<<v<<" j= "<<j<<" k= "<<k<<" l= "<<l;
      LOG_DBG(dm)<<"p_deriv: ["<<p[0]<<", "<<", "<<p[1]<<", "<<p[2];
    }
    double deriv = 0;
    if (direction == DesignElement::STIFF1){
      for (int ii = 1;ii<4;ii++) {
        for (int jj=0;jj<4;jj++) {
          for (int kk=0;kk<4;kk++) {
            deriv += coeff[(n-1)*(o-1)*j+(o-1)*k+l][ii+4*jj+16*kk]*ii*pow(u,ii-1)*pow(t,jj)*pow(v,kk);
          }
        }
      }
      deriv /= da;
    }
    if (direction == DesignElement::STIFF2) {
      for (int ii = 0;ii<4;ii++) {
        for (int jj=1;jj<4;jj++) {
          for (int kk=0;kk<4;kk++) {
            deriv += coeff[(n-1)*(o-1)*j+(o-1)*k+l][ii+4*jj+16*kk]*jj*pow(u,ii)*pow(t,jj-1)*pow(v,kk);
          }
        }
      }
      deriv /= db;
    }
    if (direction == DesignElement::STIFF3) {
      for (int ii = 0;ii<4;ii++) {
        for (int jj=0;jj<4;jj++) {
          for (int kk=1;kk<4;kk++) {
            deriv += coeff[(n-1)*(o-1)*j+(o-1)*k+l][ii+4*jj+16*kk]*kk*pow(u,ii)*pow(t,jj)*pow(v,kk-1);
          }
        }
      }
      deriv /= dc;
    }
    if (deb_) {
      LOG_DBG(dm) << "Deriv Result =" << deriv;
    }
    return deriv;
}
double DesignMaterial::EvaluateC1Interpolation(Matrix<double>& E,  Vector<double>& p,const Matrix<double> & coeff, double & da,double & db,int & j, int & k,int & m,int & n) const{
    LOG_DBG(dm) <<"p=["<<p[0]<<","<<p[1]<<"]";
    double u =(p[1]-hom_rect_b_[k][0])/(db);
    double t=(p[0]-hom_rect_a_[j][0])/(da);
    LOG_DBG(dm)<<"u = "<<u<<" t= "<<t<<"\n";
    LOG_DBG(dm)<<"u = "<<u<<" t= "<<t;
    LOG_DBG(dm)<<"j = "<<j<<" k= "<<k;
    double res = 0;
    for (int i = 0;i<4;i++) {
      for (int l=0;l<4;l++) {
        res += coeff[(n-1)*j+k][(i)*4+l]*pow(t,i)*pow(u,l);
      }
    }
    LOG_DBG(dm) << "Result =" << res;
    return res;
}

double DesignMaterial::EvaluateC1Interpolation_Deriv(Matrix<double>& E,  Vector<double>& p,const Matrix<double> & coeff, double & da,double & db,int & j, int & k,int & m,int & n, DesignElement::Type direction) const{
    double u =(p[1]-hom_rect_b_[k][0])/(db);
    double t=(p[0]-hom_rect_a_[j][0])/(da);
    LOG_DBG(dm)<<"Deriv: u = "<<u<<" t= "<<t<<"\n";

    double deriv = 0;
    if (direction == DesignElement::STIFF1){
      for (int i = 1;i<4;i++) {
        for (int l=0;l<4;l++) {
          deriv += coeff[(n-1)*j+k][(i)*4+l] *i*pow(t,i-1)*pow(u,l);
        }
      }
      deriv /= da;
    }
    if (direction == DesignElement::STIFF2) {
      for (int i = 0;i<4;i++) {
        for (int l=1;l<4;l++) {
          deriv += coeff[(n-1)*j+k][(i)*4+l] *l*pow(t,i)*pow(u,l-1);
        }
      }
      deriv /= db;
    }
    LOG_DBG(dm) << "Deriv Result =" << deriv;
    return deriv;
}

void DesignMaterial::GetLaminatesTensor(Matrix<double>& t, SubTensorType subTensor, DesignElement::Type direction, Notation notation){
  switch(subTensor){
  case PLANE_STRAIN:
  {
    t.Resize(3,3);
    t.Init();
    double eps = 5.0625e-4;
    double density = params_[DesignElement::DENSITY];
    double stiff1 = density*params_[DesignElement::STIFF1];
    double stiff2 = density*params_[DesignElement::STIFF2];
    double E = params_[DesignElement::EMODUL];
    double nu = params_[DesignElement::POISSON];
    double lambda = E*nu/((1+nu)*(1-2*nu));
    double mu = E/(2*(1+nu));
    Matrix<double> D(3,3);
    Matrix<double> Dinv(3,3);
    D.SetEntry(0,0,1/(4*(mu+lambda))+1/(4*mu));
    D.SetEntry(0,1,1/(4*(mu+lambda))-1/(4*mu));
    D.SetEntry(1,0,D(0,1));
    D.SetEntry(1,1,D(0,0));
    D.SetEntry(2,2,1/(2*mu));
    D *= 1/(eps-1);
    D.AddToEntry(0,0,stiff2/(2*mu+lambda));
    D.AddToEntry(2,2,stiff2/(2*mu));
    D.AddToEntry(1,1,stiff1*(1-stiff2)/(2*mu+lambda));
    D.AddToEntry(2,2,stiff1*(1-stiff2)/(2*mu));
    D.Invert(Dinv);
    switch(direction){
    case DesignElement::NO_DERIVATIVE:
    case DesignElement::ROTANGLE:
      D.Init();
//      mu /= eps;
//      lambda /= eps;
      D.SetEntry(0,0,2*mu+lambda);
      D.SetEntry(1,1,D(0,0));
      D.SetEntry(2,2,2*mu);
      D.SetEntry(0,1,lambda);
      D.SetEntry(1,0,lambda);
      t.Add((1-stiff2)*(1-stiff1), Dinv);
      t.Add(1.0, D);
      break;
    case DesignElement::STIFF1:
      t.SetEntry(1,1,1/(2*mu+lambda));
      t.SetEntry(2,2,1/(2*mu));
      t *= (stiff2-1)*(1-stiff2)*(1-stiff1);
      Dinv.Mult(t, D);
      D.Mult(Dinv, t);
      t.Add(stiff2-1, Dinv);
      t*=density;
      break;
    case DesignElement::STIFF2:
      t.SetEntry(0,0,1/(2*mu+lambda));
      t.SetEntry(1,1,-stiff1/(2*mu+lambda));
      t.SetEntry(2,2,(1-stiff1)/(2*mu));
      t *= (stiff2-1)*(1-stiff1);
      Dinv.Mult(t, D);
      D.Mult(Dinv, t);
      t.Add(stiff1-1, Dinv);
      t*=density;
      break;
    default:
      ZeroTensor(t, subTensor);
      return;
    }
    break;
  }
  case PLANE_STRESS:
  {
    double E33 = 1e-5;
    double density = params_[DesignElement::DENSITY];
    double stiff1 = density*params_[DesignElement::STIFF1];
    double stiff2 = density*params_[DesignElement::STIFF2];
    double E = params_[DesignElement::EMODUL];
    double nu = params_[DesignElement::POISSON];
    double n = (stiff2 + stiff1*stiff2*(nu*nu - 1) - 1);
    switch(direction)
    {
    case DesignElement::NO_DERIVATIVE:
    case DesignElement::ROTANGLE:
    {
      double E11 = -(E*stiff1)/n;
      double E22 = stiff2*E+stiff2*stiff2*nu*nu*E11;
      double E12 = stiff2*nu*E11;
      Set2dVoigtTensor(t, E11, E22, E33, 0.0, 0.0, E12);
      break;
    }
    case DesignElement::STIFF1:
    {
      double E11 = -(E*(stiff2 - 1))/(n*n);
      double E22 = stiff2*stiff2*nu*nu*E11;
      double E12 = stiff2*nu*E11;
      Set2dVoigtTensor(t, E11, E22, 0.0, 0.0, 0.0, E12);
      t*=density;
      break;
    }
    case DesignElement::STIFF2:
    {
      double E11 = (E*stiff1*(stiff1*(nu*nu - 1) + 1))/(n*n);
      double E22 = E - 2*stiff2*nu*nu*E*stiff1/n+stiff2*stiff2*nu*nu*E11;
      double E12 = -nu*E*stiff1/n+stiff2*nu*E11;
      Set2dVoigtTensor(t, E11, E22, 0.0, 0.0, 0.0, E12);
      t*=density;
    break;
    }
    default:
      ZeroTensor(t, subTensor);
        return;
      }
      break;
    }
  default:
    throw Exception("subTensor not implemented yet");
  }

  double rotAngle = params_[DesignElement::ROTANGLE];
  RotateHMStiffnessTensor(t, subTensor, direction, rotAngle, notation);
  return;
}


void DesignMaterial::ZeroTensor(Matrix<double>& t, SubTensorType subTensor){
  switch(subTensor){
  case FULL:
    t.Resize(6, 6);
    LOG_DBG(dm)<<"Zero Tensor: "<<t.ToString(2);
    break;
  case PLANE_STRAIN:
  case PLANE_STRESS:
  case PLANE:
    t.Resize(3, 3);
    break;
  default:
    throw Exception("subTensor not implemented yet");
  }
  t.Init();
}

void DesignMaterial::Set2dVoigtTensor(Matrix<double>& t, double t11, double t22, double t33, double t23, double t13, double t12)
{
  t.Resize(3,3);
  t.Init();
  t[0][0] = t11; t[0][1] = t12; t[0][2] = t13;
  t[1][0] = t12; t[1][1] = t22; t[1][2] = t23;
  t[2][0] = t13; t[2][1] = t23; t[2][2] = t33;
}


void DesignMaterial::SetTransIsoTensor(Matrix<double>& t, SubTensorType subTensor, double iD, double inD, double iG, double oD, double onD, double oG){
  switch(subTensor){
  case FULL:
    t.Resize(6, 6);
    t.Init();
    switch(transIsoType_){
    case TRANSISO_XY:
      t[0][0] = iD;  t[0][1] = inD; t[0][2] = onD;
      t[1][0] = inD; t[1][1] = iD;  t[1][2] = onD;
      t[2][0] = onD; t[2][1] = onD; t[2][2] = oD;
      t[3][3] = oG;  t[4][4] = oG;  t[5][5] = iG;
      break;
    case TRANSISO_YZ:
      t[0][0] = oD;  t[0][1] = onD; t[0][2] = onD;
      t[1][0] = onD; t[1][1] = iD;  t[1][2] = inD;
      t[2][0] = onD; t[2][1] = inD; t[2][2] = iD;
      t[3][3] = iG;  t[4][4] = oG;  t[5][5] = oG;
      break;
    case TRANSISO_XZ:
      t[0][0] = iD;  t[0][1] = onD; t[0][2] = inD;
      t[1][0] = onD; t[1][1] = oD;  t[1][2] = onD;
      t[2][0] = inD; t[2][1] = onD; t[2][2] = iD;
      t[3][3] = oG;  t[4][4] = iG;  t[5][5] = oG;
      break;
    }
    break;
  case PLANE_STRAIN:
  case PLANE_STRESS:
    t.Resize(3, 3);
    t.Init();
    switch(transIsoType_){
    case TRANSISO_XY:
      t[0][0] = iD;  t[0][1] = inD;
      t[1][0] = inD; t[1][1] = iD;
      t[2][2] = iG;
      break;
    case TRANSISO_YZ:
      t[0][0] = oD;  t[0][1] = onD;
      t[1][0] = onD; t[1][1] = iD;
      t[2][2] = oG;
      break;
    case TRANSISO_XZ:
      t[0][0] = iD;  t[0][1] = onD;
      t[1][0] = onD; t[1][1] = oD;
      t[2][2] = oG;
      break;
    }
    break;
  default:
    throw Exception("subTensor not implemented yet");
  }
}

void DesignMaterial::SetIsoTensor(Matrix<double>& t, SubTensorType subTensor, double D, double nd, double G){
  SetTransIsoTensor(t, subTensor, D, nd, G, D, nd, G);
}

void DesignMaterial::RotateHMStiffnessTensor(Matrix<double>& t, SubTensorType subTensor, DesignElement::Type direction,
    double alpha, Notation notation,double beta, double gamma){
  const double sq2inv = 1/sqrt(2);
  switch(subTensor){
  case PLANE_STRAIN:
  case PLANE_STRESS:
  case PLANE:
  {
     Matrix<double> theta(3,3);
     Matrix<double> help(3,3);
     theta.SetEntry(0,0, pow(cos(alpha),2));
     theta.SetEntry(0,1, pow(sin(alpha),2));
     theta.SetEntry(0,2, -sqrt(2)/2*sin(2*alpha));
     theta.SetEntry(1,0, theta(0,1));
     theta.SetEntry(1,1, theta(0,0));
     theta.SetEntry(1,2, -theta(0,2));
     theta.SetEntry(2,0, theta(1,2));
     theta.SetEntry(2,1, theta(0,2));
     theta.SetEntry(2,2, cos(2*alpha));
     LOG_DBG3(dm) << "RTHM: Theta =  " << theta.ToString(2) << " d=" << DesignElement::type.ToString(direction);

     t.Mult(theta, help);
     if(direction == DesignElement::ROTANGLE){
       Matrix<double> dtheta(3,3);
       dtheta.SetEntry(0,0, -sin(2*alpha));
       dtheta.SetEntry(0,1, -dtheta(0,0));
       dtheta.SetEntry(0,2, -sqrt(2)*cos(2*alpha));
       dtheta.SetEntry(1,0, dtheta(0,1));
       dtheta.SetEntry(1,1, dtheta(0,0));
       dtheta.SetEntry(1,2, -dtheta(0,2));
       dtheta.SetEntry(2,0, dtheta(1,2));
       dtheta.SetEntry(2,1, dtheta(0,2));
       dtheta.SetEntry(2,2, -2*sin(2*alpha));
       Matrix<double> dthetaTttheta(3,3);
       dtheta.MultT(help, dthetaTttheta);
       t.Mult(dtheta, help);
       theta.MultT(help, dtheta);
       t = dthetaTttheta + dtheta;
       if(notation != HILL_MANDEL)
       {
         t(0,2)*=sq2inv;
         t(1,2)*=sq2inv;
         t(2,2)/=2;
         t(2,0)*=sq2inv;
         t(2,1)*=sq2inv;
       }
       return;
     }
     theta.MultT(help, t);
     if(notation != HILL_MANDEL)
     {
       t(0,2)*=sq2inv;
       t(1,2)*=sq2inv;
       t(2,2)/=2;
       t(2,0)*=sq2inv;
       t(2,1)*=sq2inv;
     }
     return;
  }
  case FULL:
  {
    Matrix<double> Q(6, 6);
    Matrix<double> help(6, 6);
    Matrix<double> R(3, 3);
    // source: share/scripts/paraview_fmo.py
    R[0][0] = cos(beta) * cos(gamma);
    R[0][1] = -cos(beta) * sin(gamma);
    R[0][2] = sin(beta);
    R[1][0] = cos(alpha) * sin(gamma) + sin(alpha) * sin(beta) * cos(gamma);
    R[1][1] = cos(alpha) * cos(gamma) - sin(alpha) * sin(beta) * sin(gamma);
    R[1][2] = -sin(alpha) * cos(beta);
    R[2][0] = sin(alpha) * sin(gamma) - cos(alpha) * sin(beta) * cos(gamma);
    R[2][1] = sin(alpha) * cos(gamma) + cos(alpha) * sin(beta) * sin(gamma);
    R[2][2] = cos(alpha) * cos(beta);

    Q[0][0] = R[0][0] * R[0][0];
    Q[0][1] = R[0][1] * R[0][1];
    Q[0][2] = R[0][2] * R[0][2];
    Q[0][3] = 2.0 * R[0][1] * R[0][2];
    Q[0][4] = 2.0 * R[0][0] * R[0][2];
    Q[0][5] = 2.0 * R[0][0] * R[0][1];

    Q[1][0] = R[1][0] * R[1][0];
    Q[1][1] = R[1][1] * R[1][1];
    Q[1][2] = R[1][2] * R[1][2];
    Q[1][3] = 2.0 * R[1][1] * R[1][2];
    Q[1][4] = 2.0 * R[1][0] * R[1][2];
    Q[1][5] = 2.0 * R[1][0] * R[1][1];

    Q[2][0] = R[2][0] * R[2][0];
    Q[2][1] = R[2][1] * R[2][1];
    Q[2][2] = R[2][2] * R[2][2];
    Q[2][3] = 2.0 * R[2][1] * R[2][2];
    Q[2][4] = 2.0 * R[2][0] * R[2][2];
    Q[2][5] = 2.0 * R[2][0] * R[2][1];

    Q[3][0] = R[1][0] * R[2][0];
    Q[3][1] = R[1][1] * R[2][1];
    Q[3][2] = R[1][2] * R[2][2];
    Q[3][3] = R[1][1] * R[2][2] + R[1][2] * R[2][1];
    Q[3][4] = R[1][0] * R[2][2] + R[1][2] * R[2][0];
    Q[3][5] = R[1][0] * R[2][1] + R[1][1] * R[2][0];

    Q[4][0] = R[0][0] * R[2][0];
    Q[4][1] = R[0][1] * R[2][1];
    Q[4][2] = R[0][2] * R[2][2];
    Q[4][3] = R[0][1] * R[2][2] + R[0][2] * R[2][1];
    Q[4][4] = R[0][0] * R[2][2] + R[0][2] * R[2][0];
    Q[4][5] = R[0][0] * R[2][1] + R[0][1] * R[2][0];

    Q[5][0] = R[0][0] * R[1][0];
    Q[5][1] = R[0][1] * R[1][1];
    Q[5][2] = R[0][2] * R[1][2];
    Q[5][3] = R[0][1] * R[1][2] + R[0][2] * R[1][1];
    Q[5][4] = R[0][0] * R[1][2] + R[0][2] * R[1][0];
    Q[5][5] = R[0][0] * R[1][1] + R[0][1] * R[1][0];
    t.Mult(Q, help);
    if (direction == DesignElement::ROTX || direction == DesignElement::ROTY || direction == DesignElement::ROTZ) {
      Matrix<double> dQ(6, 6);
      Matrix<double> dR(3,3);
      switch (direction) {
      case DesignElement::ROTX:
        dR[0][0] =  0.;
        dR[0][1] = 0.;
        dR[0][2] =  0.;
        dR[1][0] =  -cos(alpha)*sin(gamma) + cos(alpha)*sin(beta)*cos(gamma);
        dR[1][1] =  -sin(alpha)*cos(gamma) - cos(alpha)*sin(beta)*sin(gamma);
        dR[1][2] = -cos(alpha)*cos(beta);
        dR[2][0] =  cos(alpha)*sin(gamma) + sin(alpha)*sin(beta)*cos(gamma);
        dR[2][1] =  cos(alpha)*cos(gamma) - sin(alpha)*sin(beta)*sin(gamma);
        dR[2][2] =  -sin(alpha)*cos(beta);
        break;
      case DesignElement::ROTY:
        dR[0][0] =  -sin(beta) * cos(gamma);
        dR[0][1] = sin(beta) * sin(gamma);
        dR[0][2] =  cos(beta);
        dR[1][0] =  sin(alpha)*cos(beta)*cos(gamma);
        dR[1][1] =  -sin(alpha)*cos(beta)*sin(gamma);
        dR[1][2] = sin(alpha)*cos(beta);
        dR[2][0] =  - cos(alpha)*cos(beta)*cos(gamma);
        dR[2][1] =  cos(alpha)*cos(beta)*sin(gamma);
        dR[2][2] =  -cos(alpha)*sin(beta);
        break;
      case DesignElement::ROTZ:
        dR[0][0] =  -cos(beta) * sin(gamma);
        dR[0][1] = -cos(beta) * cos(gamma);
        dR[0][2] =  0.;
        dR[1][0] =  cos(alpha)*cos(gamma) - sin(alpha)*sin(beta)*sin(gamma);
        dR[1][1] =  -cos(alpha)*sin(gamma) - sin(alpha)*sin(beta)*cos(gamma);
        dR[1][2] = 0.;
        dR[2][0] =  sin(alpha)*cos(gamma) + cos(alpha)*sin(beta)*sin(gamma);
        dR[2][1] =  -sin(alpha)*sin(gamma) + cos(alpha)*sin(beta)*cos(gamma);
        dR[2][2] =  0.;
        break;
      default:
          throw Exception("Rotation desing variable not implemented yet");
      }
      dQ[0][0] = 2.*R[0][0]*dR[0][0];
      dQ[0][1] = 2.*R[0][1]*dR[0][1];
      dQ[0][2] = 2.*R[0][2]*dR[0][2];
      dQ[0][3] = 2.0*R[0][1]*dR[0][2]+2.*R[0][2]*dR[0][1];
      dQ[0][4] = 2.0*R[0][0]*dR[0][2]+2.*R[0][2]*dR[0][0];
      dQ[0][5] = 2.0*R[0][1]*dR[0][0]+2.*R[0][0]*dR[0][1];

      dQ[1][0] = 2.*R[1][0]*dR[1][0];
      dQ[1][1] = 2.* R[1][1]*dR[1][1];
      dQ[1][2] = 2.*R[1][2]*dR[1][2];
      dQ[1][3] = 2.0*R[1][1]*dR[1][2]+2.0*R[1][2]*dR[1][1];
      dQ[1][4] = 2.0*R[1][0]*dR[1][2]+2.0*dR[1][0]*R[1][2];
      dQ[1][5] = 2.0*R[1][0]*dR[1][1]+2.0*dR[1][0]*R[1][1];

      dQ[2][0] = 2.*R[2][0]*dR[2][0];
      dQ[2][1] = 2.*R[2][1]*dR[2][1];
      dQ[2][2] = 2.*R[2][2]*dR[2][2];
      dQ[2][3] = 2.0*R[2][1]*dR[2][2]+ 2.0*dR[2][1]*R[2][2];
      dQ[2][4] = 2.0*R[2][0]*dR[2][2]+2.0*dR[2][0]*R[2][2];
      dQ[2][5] = 2.0*R[2][0]*dR[2][1]+2.0*dR[2][0]*R[2][1];

      dQ[3][0] = R[1][0]*dR[2][0] +  dR[1][0]*R[2][0];
      dQ[3][1] = R[1][1]*dR[2][1] + dR[1][1]*R[2][1];
      dQ[3][2] = R[1][2]*dR[2][2] + dR[1][2]*R[2][2];
      dQ[3][3] = R[1][1]*dR[2][2] + R[1][2]*dR[2][1] + dR[1][1]*R[2][2] + dR[1][2]*R[2][1];
      dQ[3][4] = R[1][0]*dR[2][2] + R[1][2]*dR[2][0] + dR[1][0]*R[2][2] + dR[1][2]*R[2][0];
      dQ[3][5] = R[1][0]*dR[2][1] + R[1][1]*dR[2][0] + dR[1][0]*R[2][1] + dR[1][1]*R[2][0];

      dQ[4][0] = R[0][0]*dR[2][0] + dR[0][0]*R[2][0];
      dQ[4][1] = R[0][1]*dR[2][1] + dR[0][1]*R[2][1];
      dQ[4][2] = R[0][2]*dR[2][2] + dR[0][2]*R[2][2];
      dQ[4][3] = R[0][1]*dR[2][2] + R[0][2]*dR[2][1] + dR[0][1]*R[2][2] + dR[0][2]*R[2][1];
      dQ[4][4] = R[0][0]*dR[2][2] + R[0][2]*dR[2][0] + dR[0][0]*R[2][2] + dR[0][2]*R[2][0];
      dQ[4][5] = R[0][0]*dR[2][1] + R[0][1]*dR[2][0] + dR[0][0]*R[2][1] + dR[0][1]*R[2][0];

      dQ[5][0] = R[0][0]*dR[1][0] + dR[0][0]*R[1][0];
      dQ[5][1] = R[0][1]*dR[1][1] + dR[0][1]*R[1][1];
      dQ[5][2] = R[0][2]*dR[1][2] + dR[0][2]*R[1][2];
      dQ[5][3] = R[0][1]*dR[1][2] + R[0][2]*dR[1][1] + dR[0][1]*R[1][2] + dR[0][2]*R[1][1];
      dQ[5][4] = R[0][0]*dR[1][2] + R[0][2]*dR[1][0] + dR[0][0]*R[1][2] + dR[0][2]*R[1][0];
      dQ[5][5] = R[0][0]*dR[1][1] + R[0][1]*dR[1][0] + dR[0][0]*R[1][1] + dR[0][1]*R[1][0];

      Matrix<double> dQTtQ(6, 6);
      dQ.MultT(help, dQTtQ);
      t.Mult(dQ, help);
      Q.MultT(help, dQ);
      t = dQTtQ + dQ;
      if (notation != HILL_MANDEL) {
        for (int i=0;i<6;i++) {
          for (int j=0;j<6;j++) {
            if ((i<3 && j>2) || (i>2 && j<3)) {
              t(i,j) *= sq2inv;
            } else if (i>2 && j>2) {
              t(i,j) /= 2.;
            }
          }
        }
      }
      return;
    }
    Q.MultT(help, t);
    if (notation != HILL_MANDEL) {
      for (int i=0;i<6;i++) {
        for (int j=0;j<6;j++) {
          if ((i<3 && j>2) || (i>2 && j<3)) {
            t(i,j) *= sq2inv;
          } else if (i>2 && j>2) {
            t(i,j) /= 2.;
          }
        }
      }
    }
    return;
  }
  default:
    throw Exception("subTensor not implemented yet");
  }
}

double DesignMaterial::GetTransIsoMass(double iD, double iG, double oD, double oG){
  switch(dim){
  case 2:
    switch(transIsoType_){
    case TRANSISO_XY:
      return(2.0*iD + iG);
    case TRANSISO_YZ:
    case TRANSISO_XZ:
      return(iD + oD + oG);
    default:
      throw Exception("transIsoType not implemented yet");
    }
    break;
  case 3:
    return(2.0*iD + oD + iG + 2.0*oG);
  default:
    throw Exception("strange dimension");
  }
}

double DesignMaterial::GetIsoMass(double D, double G){
  return(GetTransIsoMass(D, G, D, G));
}

void DesignMaterial::GetMaterialTensor(Matrix<double>& t, SubTensorType subTensor, DesignElement::Type direction, Notation notation)
{
  assert(!(notation == HILL_MANDEL && type_ != FMO && type_ != LAMINATES && type_ != HOM_RECT && type_ != HOM_RECT_C1));

  switch(type_){
  case FMO:
    GetAnisotropicTensor(t, direction, notation);
    break;
  case ISOTROPIC:
    GetIsoMaterialTensor(t, subTensor, direction);
    break;
  case LAME_ISOTROPIC: // LAME_ISOTROPIC
    GetLameMaterialTensor(t, subTensor, direction);
    break;
  case TRANSVERSAL_ISOTROPIC:
  case TRANSVERSAL_ISOTROPIC_BOXED:
  case DENSITY_TIMES_TRANSVERSAL_ISOTROPIC:
  case DENSITY_TIMES_TRANSVERSAL_ISOTROPIC_BOXED:
  case DENSITY_TIMES_ROT_TRANSVERSAL_ISOTROPIC_BOXED:
    GetTransIsoMaterialTensor(t, subTensor, direction);
    break;
  case DENSITY_TIMES_2D_TENSOR:
  case DENSITY_TIMES_2D_TENSOR_CONSTANT_TRACE:
  case DENSITY_TIMES_ROTATED_2D_TENSOR:
    GetDensityTimes2dTensorTensor(t, subTensor, direction);
    break;
  case LAMINATES:
    GetLaminatesTensor(t, subTensor, direction, notation);
    break;
  case HOM_RECT:
    GetHomRectTensor(t,subTensor, direction, notation);
    break;
  case HOM_RECT_C1:
    GetHomRectTensor(t, subTensor,direction, notation);
    break;
  default: // case default
    throw Exception("DesignMaterial Type not implemented yet");
  }  
}

void DesignMaterial::GetDielecTensor(Matrix<double>& t, DesignElement::Type direction)
{
  // only 2D!
  double e11 = 0;
  double e22 = 0;
  double e12 = 0;
  if(direction == DesignElement::NO_DERIVATIVE)
  {
    e11 = params_[DesignElement::DIELEC_11];
    e22 = params_[DesignElement::DIELEC_22];
    e12 = params_[DesignElement::DIELEC_12];
  }
  t.Resize(2,2);
  t.Init();

  switch(direction)
  {
  case DesignElement::NO_DERIVATIVE:
    // negative for the piezo case
    t[0][0] = -e11; t[0][1] = -e12;
    t[1][0] = -e12; t[1][1] = -e22;
    break;
  case DesignElement::DIELEC_11:
    t[0][0] = -1.0;
    break;
  case DesignElement::DIELEC_22:
    t[1][1] = -1.0;
    break;
  case DesignElement::DIELEC_12:
    t[0][1] = -1.0;
    break;
  default:
    // sensitivity is zero!
    break;
  }
}

void DesignMaterial::GetPiezoCouplingTensor(Matrix<double>& t, DesignElement::Type direction)
{
  // only 2D!
  double e11 = 0;
  double e12 = 0;
  double e13 = 0;
  double e21 = 0;
  double e22 = 0;
  double e23 = 0;
  if(direction == DesignElement::NO_DERIVATIVE)
  {
    e11 = params_[DesignElement::PIEZO_11];
    e12 = params_[DesignElement::PIEZO_12];
    e13 = params_[DesignElement::PIEZO_13];
    e21 = params_[DesignElement::PIEZO_21];
    e22 = params_[DesignElement::PIEZO_22];
    e23 = params_[DesignElement::PIEZO_23];
  }
  t.Resize(2,3);
  t.Init();

  switch(direction)
  {
  case DesignElement::NO_DERIVATIVE:
    t[0][0] = e11; t[0][1] = e12; t[0][2] = e13;
    t[1][0] = e21; t[1][1] = e22; t[1][2] = e23;
    break;
  case DesignElement::PIEZO_11:
    t[0][0] = 1.0;
    break;
  case DesignElement::PIEZO_12:
    t[0][1] = 1.0;
    break;
  case DesignElement::PIEZO_13:
    t[0][2] = 1.0;
    break;
  case DesignElement::PIEZO_21:
    t[1][0] = 1.0;
    break;
  case DesignElement::PIEZO_22:
    t[1][1] = 1.0;
    break;
  case DesignElement::PIEZO_23:
    t[1][2] = 1.0;
    break;

  default:
    // empty, sensitivity is zero
    break;
  }
}


double DesignMaterial::GetMaterialMass(DesignElement::Type direction){
  if(massIsDesign_){
    switch(direction){
    case DesignElement::MASS:
      return(massFactor_);
    case DesignElement::NO_DERIVATIVE:
      return(params_[DesignElement::MASS]*massFactor_);
    default:
      return(0.0);
    }
  }else{
    switch(type_){
    case ISOTROPIC:
      return(GetIsoMaterialMass(direction)*massFactor_);
    case LAME_ISOTROPIC: // LAME_ISOTROPIC
      return(GetLameMaterialMass(direction)*massFactor_);
    case TRANSVERSAL_ISOTROPIC:
    case TRANSVERSAL_ISOTROPIC_BOXED:
      return(GetTransIsoMaterialMass(direction)*massFactor_);
    default: // case default
      throw Exception("DesignMaterial Type not implemented yet");
    }
  }
}

bool DesignMaterial::GetMaterialDamping(double& alpha, double& beta, DesignElement::Type direction){
  if(DampingIsDesign()){
    switch(direction){
    case DesignElement::DAMPINGALPHA:
      alpha = 1.0;
      beta = 0.0;
      break;
    case DesignElement::DAMPINGBETA:
      alpha = 0.0;
      beta = 1.0;
      break;
    case DesignElement::NO_DERIVATIVE:
      alpha = params_[DesignElement::DAMPINGALPHA];
      beta = params_[DesignElement::DAMPINGBETA];
      break;
    default:
      alpha = 0.0;
      beta = 0.0;
      break;
    }
    return(true);
  }else{
    return(false);
  }
}

void DesignMaterial::DumpParams()
{
  std::map<DesignElement::Type, double>::iterator iter;

  for(iter = params_.begin(); iter != params_.end(); ++iter)
    std::cout << "params[" << DesignElement::type.ToString(iter->first) << "] = " << iter->second << std::endl;
}

void DesignMaterial::SetEnums()
{
  type.SetName("DesignMaterial::Type");
  type.Add(FMO, "fmo");
  type.Add(ISOTROPIC, "isotropic");
  type.Add(LAME_ISOTROPIC, "lame-isotropic");
  type.Add(TRANSVERSAL_ISOTROPIC, "transversal-isotropic");
  type.Add(TRANSVERSAL_ISOTROPIC_BOXED, "transversal-isotropic-boxed");
  type.Add(DENSITY_TIMES_TRANSVERSAL_ISOTROPIC, "density-times-transversal-isotropic");
  type.Add(DENSITY_TIMES_TRANSVERSAL_ISOTROPIC_BOXED, "density-times-transversal-isotropic-boxed");
  type.Add(DENSITY_TIMES_ROT_TRANSVERSAL_ISOTROPIC_BOXED, "density-times-rotated-transversal-isotropic-boxed");
  type.Add(DENSITY_TIMES_2D_TENSOR, "density-times-2dtensor");
  type.Add(DENSITY_TIMES_2D_TENSOR_CONSTANT_TRACE, "density-times-2dtensor-constant-trace");
  type.Add(DENSITY_TIMES_ROTATED_2D_TENSOR, "density-times-rotated-2dtensor");
  type.Add(LAMINATES, "laminates");
  type.Add(HOM_RECT, "hom-rect");
  type.Add(HOM_RECT_C1,"hom-rect-C1");

  transIsoType.SetName("DesignMaterial::TransIsoType");
  transIsoType.Add(TRANSISO_XY, "xy");
  transIsoType.Add(TRANSISO_YZ, "yz");
  transIsoType.Add(TRANSISO_XZ, "xz");

  notation.SetName("DesignMaterial::Notation");
  notation.Add(VOIGT, "voigt");
  notation.Add(HILL_MANDEL, "hill_mandel");
}
