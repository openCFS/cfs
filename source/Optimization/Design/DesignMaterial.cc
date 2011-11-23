#include "DesignMaterial.hh"
#include "Domain/domain.hh"

using namespace CoupledField;

Enum<DesignMaterial::Type>      DesignMaterial::type;
Enum<DesignMaterial::TransIsoType> DesignMaterial::transIsoType;

DesignMaterial::DesignMaterial(PtrParamNode pn, StdVector<DesignElement::Type>& design){
  type_ = type.Parse(pn->Get("type")->As<std::string>());
  
  dim = domain->GetGrid()->GetDim();
  
  transIsoType_ = transIsoType.Parse(pn->Get("isoplane")->As<std::string>());
  
  massIsDesign_ = pn->Get("optimizeMass")->As<bool>();
  
  massFactor_ = pn->Get("massFactor")->As<double>();
  
  penalty_ = pn->Get("penalty")->As<double>();
  
  trace_ = pn->Get("trace")->As<double>();
  
  dampingIsDesign_ = pn->Get("optimizeDamping")->As<bool>();

  // collect all designs here, to check whether all are given
  unsigned int r = RequiredParameters();
  StdVector<DesignElement::Type> d;
  d.Reserve(r);
  // copy the ones from DesignSpace
  for(unsigned int i=0; i < design.GetSize(); ++i){
    d.Push_back(design[i]);
  }
  // read non-design parameters
  ParamNodeList params = pn->GetList("param");
  for(unsigned int i=0; i < params.GetSize(); i++){
    DesignElement::Type dt = DesignElement::type.Parse(params[i]->Get("name")->As<std::string>());
    SetParameter(dt, params[i]->Get("value")->As<Double>());
    if(d.Find(dt) < 0){
      d.Push_back(dt);
    }
  }
  if(!CheckRequiredDesigns(d)){
    throw Exception("Not all Parameters for chosen DesignMaterial given");
  }else if(design.GetSize() > r){ // design.GetSize() < r is impossible as CheckRequiredDesigns passed  
    WARN("There are designs specified that are not used!");
  }
}

unsigned int DesignMaterial::RequiredParameters(){
  unsigned int r = MassIsDesign() ? 1 : 0;
  if(DampingIsDesign()){
    r += 2;
  }
  switch(type_){
  case ISOTROPIC:
  case LAME_ISOTROPIC:
    return r+2;
  case TRANSVERSAL_ISOTROPIC:
  case TRANSVERSAL_ISOTROPIC_BOXED:
    return r+5;
  case DENSITY_TIMES_2D_TENSOR_CONSTANT_TRACE:
    return r+6;
  case DENSITY_TIMES_2D_TENSOR:
    return r+7;
  default:
    throw Exception("DesignMaterial Type not implemented yet (RequiredParameters)");
  }
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
  case ISOTROPIC:
    return(design.Find(DesignElement::EMODUL) >=0
        && design.Find(DesignElement::POISSON) >= 0);
  case LAME_ISOTROPIC:
    return(design.Find(DesignElement::LAMELAMBDA) >= 0 
        && design.Find(DesignElement::LAMEMU) >= 0);
  case TRANSVERSAL_ISOTROPIC:
  case TRANSVERSAL_ISOTROPIC_BOXED:
    return(design.Find(DesignElement::EMODULISO) >= 0 
        && design.Find(DesignElement::POISSONISO) >= 0 
        && design.Find(DesignElement::EMODUL) >= 0 
        && design.Find(DesignElement::POISSON) >= 0 
        && design.Find(DesignElement::GMODUL) >= 0);
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
  default:
    throw Exception("DesignMaterial Type not implemented yet (CheckRequiredDesigns)");
  }
}

void DesignMaterial::SetParameter(const DesignElement::Type p, const double value){
  params_[p] = value;
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
  double nu = params_[DesignElement::POISSONISO];
  double nu13 = params_[DesignElement::POISSON];
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
    double G3 = params_[DesignElement::GMODUL];
    double G = 0.5*E/(1.0+nu);
    SetTransIsoTensor(t, subTensor, D, nD, G, D3, nD3, G3);
    return;
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

void DesignMaterial::GetDensityTimes2dTensorTensor(Matrix<double>& t, SubTensorType subTensor, DesignElement::Type direction){
  double e11 = 0;
  double e22 = 0;
  double e33 = 0;
  double e23 = 0;
  double e13 = 0;
  double e12 = 0;
  if(direction == DesignElement::NO_DERIVATIVE || direction == DesignElement::DENSITY){
   e11 = params_[DesignElement::TENSOR11];
   e22 = params_[DesignElement::TENSOR22];
   if(type_ == DENSITY_TIMES_2D_TENSOR_CONSTANT_TRACE){
     e33 = 0.5 * (trace_ - e11 - e22); 
   }else{
     e33 = params_[DesignElement::TENSOR33];
   }
   e23 = params_[DesignElement::TENSOR23];
   e13 = params_[DesignElement::TENSOR13];
   e12 = params_[DesignElement::TENSOR12];
  }
  double d = params_[DesignElement::DENSITY];
  switch(direction){
  case DesignElement::NO_DERIVATIVE:
    d = exp(log(d)*penalty_);
    Set2dVoigtTensor(t, subTensor, d*e11, d*e22, d*e33, d*e23, d*e13, d*e12);
    break;
  case DesignElement::DENSITY:
    if(penalty_ == 1.0){
      d = 1.0;
    }else{
      d = penalty_*exp(log(d)*(penalty_-1));
    }
    Set2dVoigtTensor(t, subTensor, d*e11, d*e22, d*e33, d*e23, d*e13, d*e12);
    break;
  case DesignElement::TENSOR11:
    d = exp(log(d)*penalty_);
    Set2dVoigtTensor(t, subTensor, d, 0.0, type_ == DENSITY_TIMES_2D_TENSOR_CONSTANT_TRACE ? -0.5*d : 0.0, 0.0, 0.0, 0.0);
    break;
  case DesignElement::TENSOR22:
    d = exp(log(d)*penalty_);
    Set2dVoigtTensor(t, subTensor, 0.0, d, type_ == DENSITY_TIMES_2D_TENSOR_CONSTANT_TRACE ? -0.5*d : 0.0, 0.0, 0.0, 0.0);
    break;
  case DesignElement::TENSOR33:
    d = exp(log(d)*penalty_);
    Set2dVoigtTensor(t, subTensor, 0.0, 0.0, d, 0.0, 0.0, 0.0);
    break;
  case DesignElement::TENSOR23:
    d = exp(log(d)*penalty_);
    Set2dVoigtTensor(t, subTensor, 0.0, 0.0, 0.0, d, 0.0, 0.0);
    break;
  case DesignElement::TENSOR13:
    d = exp(log(d)*penalty_);
    Set2dVoigtTensor(t, subTensor, 0.0, 0.0, 0.0, 0.0, d, 0.0);
    break;
  case DesignElement::TENSOR12:
    d = exp(log(d)*penalty_);
    Set2dVoigtTensor(t, subTensor, 0.0, 0.0, 0.0, 0.0, 0.0, d);
    break;
  default:
    ZeroTensor(t, subTensor);
    return;
  }
}
 
void DesignMaterial::ZeroTensor(Matrix<double>& t, SubTensorType subTensor){
  switch(subTensor){
  case FULL:
    t.Resize(6, 6);
    break;
  case PLANE_STRAIN:
    t.Resize(3, 3);
    break;
  default:
    throw Exception("subTensor not implemented yet");
  }
  t.Init();
}

void DesignMaterial::Set2dVoigtTensor(Matrix<double>& t, SubTensorType subTensor, double t11, double t22, double t33, double t23, double t13, double t12){
  switch(subTensor){
  case PLANE_STRAIN:
    t.Resize(3,3);
    t.Init();
    t[0][0] = t11; t[0][1] = t12; t[0][2] = t13;
    t[1][0] = t12; t[1][1] = t22; t[1][2] = t23;
    t[2][0] = t13; t[2][1] = t23; t[2][2] = t33;
    break;
  default:
    throw Exception("subTensor not implemented yet");
  }
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

void DesignMaterial::GetMaterialTensor(Matrix<double>& t, SubTensorType subTensor, DesignElement::Type direction){
  switch(type_){
  case ISOTROPIC:
    GetIsoMaterialTensor(t, subTensor, direction);
    break;
  case LAME_ISOTROPIC: // LAME_ISOTROPIC
    GetLameMaterialTensor(t, subTensor, direction);
    break;
  case TRANSVERSAL_ISOTROPIC:
  case TRANSVERSAL_ISOTROPIC_BOXED:
    GetTransIsoMaterialTensor(t, subTensor, direction);
    break;
  case DENSITY_TIMES_2D_TENSOR:
  case DENSITY_TIMES_2D_TENSOR_CONSTANT_TRACE:
    GetDensityTimes2dTensorTensor(t, subTensor, direction);
    break;
  default: // case default
    throw Exception("DesignMaterial Type not implemented yet");
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

void DesignMaterial::SetEnums(){
  type.SetName("DesignMaterial::Type");
  type.Add(ISOTROPIC, "isotropic");
  type.Add(LAME_ISOTROPIC, "lame-isotropic");
  type.Add(TRANSVERSAL_ISOTROPIC, "transversal-isotropic");
  type.Add(TRANSVERSAL_ISOTROPIC_BOXED, "transversal-isotropic-boxed");
  type.Add(DENSITY_TIMES_2D_TENSOR, "density-times-2dtensor");
  type.Add(DENSITY_TIMES_2D_TENSOR_CONSTANT_TRACE, "density-times-2dtensor-constant-trace");
  transIsoType.SetName("DesignMaterial::TransIsoType");
  transIsoType.Add(TRANSISO_XY, "xy");
  transIsoType.Add(TRANSISO_YZ, "yz");
  transIsoType.Add(TRANSISO_XZ, "xz");
}

