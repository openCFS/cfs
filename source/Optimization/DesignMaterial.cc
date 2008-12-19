#include "DesignMaterial.hh"

using namespace CoupledField;

Enum<DesignMaterial::Type>      DesignMaterial::type;

DesignMaterial::DesignMaterial(ParamNode* pn){
  type_ = type.Parse(pn->Get("type")->AsString());
  
  // read non-design parameters
  StdVector<ParamNode*> params = pn->GetList("param");
  nparams_ = 0;
  for(unsigned int i=0; i < params.GetSize(); i++){
    SetParameter(DesignElement::type.Parse(params[i]->Get("name")->AsString()), params[i]->Get("value")->AsDouble());
    nparams_++;
  }
}

unsigned int DesignMaterial::RequiredParameters(){
  switch(type_){
  case ISOTROPIC:
    return 2;
  case LAME_ISOTROPIC:
    return 2;
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
    switch(subTensor){
    case FULL:
      t.Resize(6, 6);
      t.Init();
      t[0][0] = diag;
      t[0][1] = lambda;
      t[0][2] = lambda;
      t[1][0] = lambda;
      t[1][1] = diag;
      t[1][2] = lambda;
      t[2][0] = lambda;
      t[2][1] = lambda;
      t[2][2] = diag;
      t[3][3] = mu;
      t[4][4] = mu;
      t[5][5] = mu;
      break;
    case PLANE_STRAIN:
      t.Resize(3, 3);
      t.Init();
      t[0][0] = diag;
      t[0][1] = lambda;
      t[1][0] = lambda;
      t[1][1] = diag;
      t[2][2] = mu;
      break;
    default:
      throw Exception("SubTensor not implemented yet");
    }
    break;
  }
  case LAME_ISOTROPIC:
  {
    double lambda = params_[DesignElement::LAMELAMBDA];
    double mu = params_[DesignElement::LAMEMU];
    double diag = lambda + 2*mu;
    switch(subTensor){
    case FULL:
      t.Resize(6,6);
      t.Init();
      t[0][0] = diag;
      t[0][1] = lambda;
      t[0][2] = lambda;
      t[1][0] = lambda;
      t[1][1] = diag;
      t[1][2] = lambda;
      t[2][0] = lambda;
      t[2][1] = lambda;
      t[2][2] = diag;
      t[3][3] = mu;
      t[4][4] = mu;
      t[5][5] = mu;
      break;
    case PLANE_STRAIN:
      t.Resize(3, 3);
      t.Init();
      t[0][0] = diag;
      t[0][1] = lambda;
      t[1][0] = lambda;
      t[1][1] = diag;
      t[2][2] = mu;
      break;
    default:
      throw Exception("SubTensor not implemented yet");
    }
    break;
  }  
  default:
    throw Exception("DesignMaterial Type not implemented yet");
  }
}
 
void DesignMaterial::GetMaterialTensorDerivative(Matrix<double>& t, SubTensorType subTensor, DesignElement::Type direction){
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
      switch(subTensor){
      case FULL: // ISOTROPIC, EMODUL, FULL
        t[0][0] = ddiag_dE;
        t[0][1] = dlambda_dE;
        t[0][2] = dlambda_dE;
        t[1][0] = dlambda_dE;
        t[1][1] = ddiag_dE;
        t[1][2] = dlambda_dE;
        t[2][0] = dlambda_dE;
        t[2][1] = dlambda_dE;
        t[2][2] = ddiag_dE;
        t[3][3] = dmu_dE;
        t[4][4] = dmu_dE;
        t[5][5] = dmu_dE;
        break;
      case PLANE_STRAIN: // ISOTROPIC, EMODUL, PLANE_STRAIN
        t[0][0] = ddiag_dE;
        t[0][1] = dlambda_dE;
        t[1][0] = dlambda_dE;
        t[1][1] = ddiag_dE;
        t[2][2] = dmu_dE;
        break;
      default:
        throw Exception("SubTensor not implemented yet");
      } // switch(subTensor)
    } // case ISOTROPIC, EMODUL
    break;
    case DesignElement::POISSON: // ISOTROPIC, POISSON
    {
      double dlambda_dnu = (1+2*nu*nu)*E / ((1+nu)*(1+nu)*(1-2*nu)*(1-2*nu));
      double dmu_dnu = E / (-2*(1+nu)*(1+nu));
      double ddiag_dnu = dlambda_dnu + 2*dmu_dnu;
      switch(subTensor){
      case FULL: // ISOTROPIC, POISSON, FULL
        t[0][0] = ddiag_dnu;
        t[0][1] = dlambda_dnu;
        t[0][2] = dlambda_dnu;
        t[1][0] = dlambda_dnu;
        t[1][1] = ddiag_dnu;
        t[1][2] = dlambda_dnu;
        t[2][0] = dlambda_dnu;
        t[2][1] = dlambda_dnu;
        t[2][2] = ddiag_dnu;
        t[3][3] = dmu_dnu;
        t[4][4] = dmu_dnu;
        t[5][5] = dmu_dnu;
        break;
      case PLANE_STRAIN: // ISOTROPIC, POISSON, PLANE_STRAIN
        t[0][0] = ddiag_dnu;
        t[0][1] = dlambda_dnu;
        t[1][0] = dlambda_dnu;
        t[1][1] = ddiag_dnu;
        t[2][2] = dmu_dnu;
        break;
      default:
        throw Exception("SubTensor not implemented yet");
      } // switch(subTensor)
    } // case ISOTROPIC, POISSON
    break;
    default:; // any derivative in any direction other than EMODUL or POISSON is zero
    } // switch(direction)
  } // case ISOTROPIC
  break;
  case LAME_ISOTROPIC: // LAME_ISOTROPIC
  {
    switch(direction){
    case DesignElement::LAMELAMBDA: // LAME_ISOTROPIC, LAMELAMBDA
    {
      switch(subTensor){
      case FULL: // LAME_ISOTROPIC, LAMELAMBDA, FULL
        t[0][0] = 1; t[0][1] = 1; t[0][2] = 1;
        t[1][0] = 1; t[1][1] = 1; t[1][2] = 1;
        t[2][0] = 1; t[2][1] = 1; t[2][2] = 1;
      break;
      case PLANE_STRAIN: // LAME_ISOTROPIC, LAMELAMBDA, PLANE_STRAIN
        t[0][0] = 1; t[0][1] = 1;
        t[1][0] = 1; t[1][1] = 1;
      break;
      default:
        throw Exception("SubTensor not implemented yet");
      } // switch(subTensor)      
    } // case LAME_ISOTROPIC, LAMELAMBDA
    break;
    case DesignElement::LAMEMU: // LAME_ISOTROPIC, LAMEMU
    {
      switch(subTensor){
      case FULL: // LAME_ISOTROPIC, LAMEMU, FULL
        t[0][0] = 2; t[1][1] = 2; t[2][2] = 2; t[3][3] = 1; t[4][4] = 1; t[5][5] = 1;
      break;
      case PLANE_STRAIN: // LAME_ISOTROPIC, LAMEMU, PLANE_STRAIN
        t[0][0] = 2; t[1][1] = 2; t[2][2] = 1;
      break;
      default:
        throw Exception("SubTensor not implemented yet");
      }
    } // case LAME_ISOTROPIC, LAMEMU
    break;
    default:; // any derivative in any direction other than EMODUL or POISSON is zero
    } // switch(direction)
  } // case LAME_ISOTROPIC
  break;
  default: // case default
    throw Exception("DesignMaterial Type not implemented yet");
  }  
}

void DesignMaterial::SetEnums(){
  type.SetName("DesignMaterial::Type");
  type.Add(ISOTROPIC, "isotropic");
  type.Add(LAME_ISOTROPIC, "lame-isotropic");
  type.Add(TRANSVERSAL_ISOTROPIC, "transversal-isotropic");
}

