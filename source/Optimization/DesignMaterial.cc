#include "DesignMaterial.hh"

using namespace CoupledField;

Enum<DesignMaterial::Type>      DesignMaterial::type;
Enum<DesignMaterial::TransIsoType> DesignMaterial::transIsoType;

DesignMaterial::DesignMaterial(PtrParamNode pn){
  type_ = type.Parse(pn->Get("type")->As<std::string>());
  
  // read non-design parameters
  ParamNodeList params = pn->GetList("param");
  nparams_ = 0;
  for(unsigned int i=0; i < params.GetSize(); i++){
    SetParameter(DesignElement::type.Parse(params[i]->Get("name")->As<std::string>()), params[i]->Get("value")->As<Double>());
    nparams_++;
  }

  transIsoType_ = transIsoType.Parse(pn->Get("isoplane")->As<std::string>());
}

unsigned int DesignMaterial::RequiredParameters(){
  switch(type_){
  case ISOTROPIC:
    return 2;
  case LAME_ISOTROPIC:
    return 2;
  case TRANSVERSAL_ISOTROPIC:
  case TRANSVERSAL_ISOTROPIC_BOXED:
    return 5;
  default:
    throw Exception("DesignMaterial Type not implemented yet");
  }
}

void DesignMaterial::SetParameter(const DesignElement::Type p, const double value){
  params_[p] = value;
}
 
void DesignMaterial::GetMaterialTensor(Matrix<double>& t, SubTensorType subTensor){
  switch(type_){
  case ISOTROPIC:
  {
    double E = params_[DesignElement::EMODUL];
    double nu = params_[DesignElement::POISSON];
    double lambda = nu * E / ((1+nu)*(1-2*nu));
    double mu = E / (2*(1+nu));
    double diag = lambda + 2*mu;
    SetIsoTensor(t, subTensor, diag, lambda, mu);
      break;
    }
  case LAME_ISOTROPIC:
  {
    double lambda = params_[DesignElement::LAMELAMBDA];
    double mu = params_[DesignElement::LAMEMU];
    double diag = lambda + 2*mu;
    SetIsoTensor(t, subTensor, diag, lambda, mu);
      break;
    }
  case TRANSVERSAL_ISOTROPIC:
  case TRANSVERSAL_ISOTROPIC_BOXED:
  {
    double E = params_[DesignElement::EMODULISO];
    double E3 = params_[DesignElement::EMODUL];
    double nu = params_[DesignElement::POISSONISO];
    double G3 = params_[DesignElement::GMODUL];
    double G = 0.5*E/(1+nu);
    double nu13 = params_[DesignElement::POISSON]; // nu_io
    double nu3 = nu13 * E3/E; // nu_oi
    double n3 = nu3*nu3*E/E3;
    double c = (1-nu-2*n3); // this is the interesting thing, this must not get 0, however this would imply a volume (trace of tensor) of infinity, so it is hopefully not occuring
    if(type_ == TRANSVERSAL_ISOTROPIC){
      if(c < 1e-8){
        c = 1e-8;
      }
    }else{ // here the parameter POISSON is really: nu_io*sqrt(2*E3/(1-nu_iso)/E) and valid values for (poisson,poissoniso) are [-1,1]²
      nu3 = sqrt(0.5*(1-nu)*E3/E)*nu13; // this is nu_oi
      n3 = 0.5*(1-nu)*nu13*nu13; // this is nu_oi*nu_io
      c = (1-nu-2*n3);
    }
    double f = E/((1+nu)*c);
    double D = (1-n3)*f;
    double D3 = (1-nu)*E3/c;
    double nD = (nu+n3)*f;
    double nD3 = (1+nu)*nu3*f;
    SetTransIsoTensor(t, subTensor, D, nD, G, D3, nD3, G3);
    break;
  }  
  default:
    throw Exception("DesignMaterial Type not implemented yet");
  }
}
 
void DesignMaterial::GetMaterialTensorDerivative(Matrix<double>& t, SubTensorType subTensor, DesignElement::Type direction){
  switch(type_){
  case ISOTROPIC:
  {
    double E = params_[DesignElement::EMODUL];
    double nu = params_[DesignElement::POISSON];
    switch(direction){
    case DesignElement::EMODUL: // ISOTROPIC, EMODUL
    {
      double dlambda_dE = nu/((1+nu)*(1-2*nu));
      double dmu_dE = 1/(2*(1+nu));
      double ddiag_dE = dlambda_dE + 2*dmu_dE;
      SetIsoTensor(t, subTensor, ddiag_dE, dlambda_dE, dmu_dE);
    } // case ISOTROPIC, EMODUL
    break;
    case DesignElement::POISSON: // ISOTROPIC, POISSON
    {
      double dlambda_dnu = (1+2*nu*nu)*E / ((1+nu)*(1+nu)*(1-2*nu)*(1-2*nu));
      double dmu_dnu = E / (-2*(1+nu)*(1+nu));
      double ddiag_dnu = dlambda_dnu + 2*dmu_dnu;
      SetIsoTensor(t, subTensor, ddiag_dnu, dlambda_dnu, dmu_dnu);
    } // case ISOTROPIC, POISSON
        break;
      default:
      ZeroTensor(t, subTensor); // any derivative in any direction other than EMODUL or POISSON is zero
    } // switch(direction)
  } // case ISOTROPIC
  break;
  case LAME_ISOTROPIC: // LAME_ISOTROPIC
  {
    switch(direction){
    case DesignElement::LAMELAMBDA: // LAME_ISOTROPIC, LAMELAMBDA
      SetIsoTensor(t, subTensor, 1, 1, 0);
      break;
    case DesignElement::LAMEMU: // LAME_ISOTROPIC, LAMEMU
      SetIsoTensor(t, subTensor, 2, 0, 1);
      break;
      default:
      ZeroTensor(t, subTensor); // any derivative in any direction other than EMODUL or POISSON is zero
    } // switch(direction)
  } // case LAME_ISOTROPIC
    break;
  case TRANSVERSAL_ISOTROPIC:
  case TRANSVERSAL_ISOTROPIC_BOXED:
    {
    double E = params_[DesignElement::EMODULISO];
    double E3 = params_[DesignElement::EMODUL];
    double nu = params_[DesignElement::POISSONISO];
    double nu13 = params_[DesignElement::POISSON];
    double nu3 = nu13 * E3/E;
    double n3 = nu3*nu3*E/E3;
    double c = (1-nu-2*n3); // this is the interesting thing, this must not get 0, however this would imply a volume (trace of tensor) of infinity, so it is hopefully not occuring
    if(type_ == TRANSVERSAL_ISOTROPIC){
      if(c < 1e-8) {
        c = 1e-8;
      }
    }else{
      nu3 = sqrt(0.5*(1-nu)*E3/E)*nu13;
      n3 = nu3*nu3*E/E3;
      c = (1-nu-2*n3);
    }
    double f = E/((1+nu)*c);
    double dE = 0, dE3 = 0, dnu = 0, dnu3 = 0, dn3 = 0, dG3 = 0;
    switch(direction){
    case DesignElement::EMODULISO:
      dE = 1;
      if(type_ == TRANSVERSAL_ISOTROPIC){
        dnu3 = -E3*nu13/(E*E);
        dn3 = nu3/E3 * (2*E*dnu3 + nu3);
      }else{
        dnu3 = -sqrt(0.125*(1-nu)*E3/E)*nu13/E;
      }
      break;
    case DesignElement::EMODUL:
      dE3 = 1;
      if(type_ == TRANSVERSAL_ISOTROPIC){
        dnu3 = nu13/E;
        dn3 = nu3*E/E3 * (2*dnu3 - nu3/E3);
      }else{
        dnu3 = sqrt(0.125*(1-nu)/(E*E3))*nu13;
      }
      break;
    case DesignElement::POISSONISO:
      dnu = 1;
      if(type_ == TRANSVERSAL_ISOTROPIC_BOXED){ // else = 0
        dnu3 = -sqrt(0.125*E3/(E*(1-nu)))*nu13;
        dn3 = -0.5*nu13*nu13;
      }
    break;
    case DesignElement::POISSON:
      if(type_ == TRANSVERSAL_ISOTROPIC){
        dnu3 = 1;
        dn3 = 2*nu3*E/E3*dnu3;
      }else{
        dnu3 = sqrt(0.5*(1-nu)*E3/E);
        dn3 = (1-nu)*nu13;
      }
      break;
    case DesignElement::GMODUL:
      dG3 = 1;
      break;
    default:;
    } // switch(direction)
    double dc = -dnu-2*dn3;
    double df = ( dE - E*dnu/(1+nu) - E*dc/c ) / ((1+nu)*c);
    double dD = (1-n3)*df - dn3*f;
    double dnD = (nu+n3)*df + (dnu+dn3)*f;
    double dD3 = ( (1-nu)*dE3 - dnu*E3 - (1-nu)*E3*dc/c ) / c;
    double dnD3 = (1+nu)*nu3*df + (1+nu)*dnu3*f + dnu*nu3*f;
    double dG = 0.5 * ( (1+nu)*dE - E*dnu ) / ( (1+nu)*(1+nu) );
    SetTransIsoTensor(t, subTensor, dD, dnD, dG, dD3, dnD3, dG3);
  } // case TRANSVERSAL_ISOMETRIC
  break;
  default: // case default
    throw Exception("DesignMaterial Type not implemented yet");
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

void DesignMaterial::SetEnums(){
  type.SetName("DesignMaterial::Type");
  type.Add(ISOTROPIC, "isotropic");
  type.Add(LAME_ISOTROPIC, "lame-isotropic");
  type.Add(TRANSVERSAL_ISOTROPIC, "transversal-isotropic");
  type.Add(TRANSVERSAL_ISOTROPIC_BOXED, "transversal-isotropic-boxed");
  transIsoType.SetName("DesignMaterial::TransIsoType");
  transIsoType.Add(TRANSISO_XY, "xy");
  transIsoType.Add(TRANSISO_YZ, "yz");
  transIsoType.Add(TRANSISO_XZ, "xz");
}

