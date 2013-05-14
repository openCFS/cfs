#include <math.h>
#include <string>

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

  if((type_ == REDBAS_PARAM) | (type_ == REDBAS_FREE))
  {
    std::string file = pn->Get("modRed/file")->As<std::string>();
    Xerces xerces(file);
    PtrParamNode root = xerces.CreateParamNodeInstance();
    PtrParamNode hr = root->Get("modRed");
    dimension_ = hr->Get("dimension")->As<UInt>();
    all_param_ = hr->Get("allParameters")->As<bool>();
    PtrParamNode mean = hr->Get("meantensor");
    PtrParamNode tensor = mean->Get("tensor");
    Matrix<double> mat(3,1);
    ParamTools::AsTensor<double>(tensor->Get("real"),3, 1, mat);
    mean_tensor_ = mat;

    StdVector<Matrix<double> > matrices(30);
    StdVector<Vector<double> > vectors(12);
    mod_red_matrices_ = matrices;
    mod_red_vectors_ = vectors;
//    PtrParamNode matnode = hr->Get("matrices");
//    PtrParamNode vecnode = hr->Get("vectors");

    StdVector<std::string> tensor_comp(3);// =["11", "12", "33"];
    tensor_comp[0]="11";
    tensor_comp[1]="12";
    tensor_comp[2]="33";

    for (int tensor_int =0; tensor_int <3; tensor_int ++)
    {

        FillModRedMatrices(hr, tensor_comp, tensor_int, dimension_);
        FillModRedVectors(hr, tensor_comp,tensor_int, dimension_);
    }
  }


  if((type_ == GREEDY_FREE) | (type_ == GREEDY_PARAM))
  {
    std::string file = pn->Get("modRed/file")->As<std::string>();
    Xerces xerces(file);
    PtrParamNode root = xerces.CreateParamNodeInstance();
    PtrParamNode hr = root->Get("modRed");

    //Here dimension is the number of terms in the tensor-product expansion
    dimension_ = hr->Get("dimension")->As<UInt>();
   // dimension_tot_ = hr->Get("dimensionTotale")->As<UInt>();
    Nl_ = hr->Get("Nl")->As<UInt>();
    Na_ = hr->Get("Na")->As<UInt>();
    lmin_ = hr->Get("lmin")->As<double>();
    lmax_ = hr->Get("lmax")->As<double>();

    all_param_ = hr->Get("allParameters")->As<bool>();

    //If type_ == GREEDY_FREE, necessarily, we must have all_param_
    assert((type_==GREEDY_PARAM) || (all_param_));

    PtrParamNode mean = hr->Get("meantensor");
    PtrParamNode tensor = mean->Get("tensor");
    Matrix<double> mat(3,1);
    ParamTools::AsTensor<double>(tensor->Get("real"),3, 1, mat);
    mean_tensor_ = mat;

    int num_param=3;
    if (all_param_) num_param = 4;

    StdVector<Matrix<double> > mat_param(num_param);
    matrices_param_ = mat_param;

///////////////Initialisation of matl2////////////////////////

    mat.Resize(Nl_+1,3*dimension_);


    PtrParamNode matl2 = hr->Get("matl2");
    tensor = matl2->Get("tensor");
    Matrix<double> mattrans(Nl_+1,3*dimension_);
    ParamTools::AsTensor<double>(tensor->Get("real"),Nl_+1, 3*dimension_, mat);
    matrices_param_[num_param-1] = mat;

  /*  PtrParamNode matl2 = hr->Get("matl2");
    tensor = matl2->Get("tensor");
    Matrix<double> mattrans(Nl_+1,3*dimension_tot_);
    ParamTools::AsTensor<double>(tensor->Get("real"),Nl_+1, 3*dimension_tot_, mattrans);


    mat.Resize(Nl_+1,3*dimension_);
    mat.Init();

    Matrix<double> mathelp(Nl_+1,dimension_);
    mattrans.GetSubMatrix(mathelp,0,0);
    mat.SetSubMatrix(mathelp,0,0);
    mattrans.GetSubMatrix(mathelp,0,dimension_tot_);
    mat.SetSubMatrix(mathelp,0,dimension_);
    mattrans.GetSubMatrix(mathelp,0,2*dimension_tot_);
    mat.SetSubMatrix(mathelp,0,2*dimension_);

    matrices_param_[num_param-1] = mat;*/

/////////Initialisation of matl1 /////////////////////////////

    mat.Resize(Nl_+1,3*dimension_);
    PtrParamNode matl1 = hr->Get("matl1");
    tensor = matl1->Get("tensor");
    ParamTools::AsTensor<double>(tensor->Get("real"),Nl_+1, 3*dimension_, mat);
    matrices_param_[num_param-2] = mat;

    /*mattrans.Resize(Nl_+1,3*dimension_tot_);
    mattrans.Init();

    PtrParamNode matl1 = hr->Get("matl1");
    tensor = matl1->Get("tensor");
    ParamTools::AsTensor<double>(tensor->Get("real"),Nl_+1, 3*dimension_tot_, mattrans);

    mat.Resize(Nl_+1,3*dimension_);
    mat.Init();

    mathelp.Resize(Nl_+1,dimension_);
    mathelp.Init();
    mattrans.GetSubMatrix(mathelp,0,0);
    mat.SetSubMatrix(mathelp,0,0);
    mattrans.GetSubMatrix(mathelp,0,dimension_tot_);
    mat.SetSubMatrix(mathelp,0,dimension_);
    mattrans.GetSubMatrix(mathelp,0,2*dimension_tot_);
    mat.SetSubMatrix(mathelp,0,2*dimension_);

    matrices_param_[num_param-2] = mat;*/



///////////////Inintialisation of matphi//////////////////////////
  mat.Resize(2*Na_+1,3*dimension_);
  PtrParamNode matphi = hr->Get("matphi");
  tensor = matphi->Get("tensor");
  ParamTools::AsTensor<double>(tensor->Get("real"),2*Na_+1, 3*dimension_, mat);

  matrices_param_[num_param-3] = mat;


   /* mattrans.Resize(2*Na_+1,3*dimension_tot_);
    mattrans.Init();

    PtrParamNode matphi = hr->Get("matphi");
    tensor = matphi->Get("tensor");
    ParamTools::AsTensor<double>(tensor->Get("real"),2*Na_+1, 3*dimension_tot_, mattrans);


    mat.Resize(2*Na_+1,3*dimension_);
    mat.Init();
    matrices_param_[num_param-3] = mat;*/

   // std::cout << matrices_param_[num_param-3](0,0);

    if (all_param_)
    {

          PtrParamNode mattheta = hr->Get("mattheta");
          tensor = mattheta->Get("tensor");
          mat.Resize(2*Na_+1,3*dimension_);
          ParamTools::AsTensor<double>(tensor->Get("real"),2*Na_+1, 3*dimension_, mat);
          matrices_param_[num_param-4] = mat;
    }


    StdVector<Matrix<double> > matrices(30);
    StdVector<Vector<double> > vectors(12);
    mod_red_matrices_ = matrices;
    mod_red_vectors_ = vectors;

        StdVector<std::string> tensor_comp(3);// =["11", "12", "33"];
           tensor_comp[0]="11";
           tensor_comp[1]="12";
           tensor_comp[2]="33";

           for (int tensor_int =0; tensor_int <3; tensor_int ++)
           {

               FillModRedMatrices(hr, tensor_comp, tensor_int, 3*dimension_);
               FillModRedVectors(hr, tensor_comp,tensor_int, 3*dimension_);
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




void DesignMaterial::FillModRedMatrices(PtrParamNode matnode, const StdVector<std::string>& tensor_comp, const int& tensor_int, const UInt& dimbas)
{

      PtrParamNode mat1 = matnode->Get("matuxux" + tensor_comp[tensor_int]);
      PtrParamNode tensor = mat1->Get("tensor");
      Matrix<double> mat(dimbas,dimbas);
      ParamTools::AsTensor<double>(tensor->Get("real"),dimbas, dimbas, mat);
      mod_red_matrices_[10*tensor_int +0] = mat;

      PtrParamNode mat2= matnode->Get("matuxvx" + tensor_comp[tensor_int] );
      tensor = mat2->Get("tensor");
      ParamTools::AsTensor<double>(tensor->Get("real"),dimbas, dimbas, mat);
       mod_red_matrices_[10*tensor_int + 1] = mat;

      PtrParamNode mat3= matnode->Get("matvxvx" + tensor_comp[tensor_int]);
      tensor = mat3->Get("tensor");
      ParamTools::AsTensor<double>(tensor->Get("real"),dimbas, dimbas, mat);
      mod_red_matrices_[10*tensor_int + 2] = mat;


      PtrParamNode mat4 = matnode->Get("matuxuy" + tensor_comp[tensor_int]);
      tensor = mat4->Get("tensor");
      ParamTools::AsTensor<double>(tensor->Get("real"),dimbas, dimbas, mat);
      mod_red_matrices_[ 10*tensor_int +3] = mat;

      PtrParamNode mat5 = matnode->Get("matuxvy" + tensor_comp[tensor_int]);
      tensor = mat5->Get("tensor");
      ParamTools::AsTensor<double>(tensor->Get("real"),dimbas, dimbas, mat);
       mod_red_matrices_[10*tensor_int + 4] = mat;

       PtrParamNode mat6 = matnode->Get("matuyvx" + tensor_comp[tensor_int]);
       tensor = mat6->Get("tensor");
       ParamTools::AsTensor<double>(tensor->Get("real"),dimbas, dimbas, mat);
       mod_red_matrices_[ 10*tensor_int + 5] = mat;

      PtrParamNode mat7 = matnode->Get("matvxvy" + tensor_comp[tensor_int]);
      tensor = mat7->Get("tensor");
      ParamTools::AsTensor<double>(tensor->Get("real"),dimbas, dimbas, mat);
      mod_red_matrices_[10*tensor_int + 6] = mat;


      PtrParamNode mat8 = matnode->Get("matuyuy" + tensor_comp[tensor_int]);
      tensor = mat8->Get("tensor");
      ParamTools::AsTensor<double>(tensor->Get("real"),dimbas, dimbas, mat);
       mod_red_matrices_[10*tensor_int + 7] = mat;

       PtrParamNode mat9 = matnode->Get("matuyvy" + tensor_comp[tensor_int]);
       tensor = mat9->Get("tensor");
       ParamTools::AsTensor<double>(tensor->Get("real"),dimbas, dimbas, mat);
       mod_red_matrices_[10*tensor_int + 8] = mat;

      PtrParamNode mat10 = matnode->Get("matvyvy" + tensor_comp[tensor_int]);
      tensor = mat10->Get("tensor");
      ParamTools::AsTensor<double>(tensor->Get("real"),dimbas, dimbas, mat);
      mod_red_matrices_[10*tensor_int + 9] = mat;

}

void DesignMaterial::FillModRedVectors(PtrParamNode vecnode, const StdVector<std::string>& tensor_comp, const int& tensor_int, const UInt& dimbas)
{

      PtrParamNode vecux = vecnode->Get("vecux" + tensor_comp[tensor_int] );
      PtrParamNode tensor = vecux->Get("tensor");
      Matrix<double> mat(dimbas,1);
      Vector<double> vec(dimbas);
      ParamTools::AsTensor<double>(tensor->Get("real"),dimbas, 1, mat);
      mat.GetCol(vec,0);
      mod_red_vectors_[4*tensor_int +0] = vec;

      PtrParamNode vecuy = vecnode->Get("vecuy" + tensor_comp[tensor_int]);
      tensor = vecuy->Get("tensor");
      ParamTools::AsTensor<double>(tensor->Get("real"),dimbas, 1, mat);
      mat.GetCol(vec,0);
      mod_red_vectors_[4*tensor_int +1] = vec;

      PtrParamNode vecvx = vecnode->Get("vecvx" + tensor_comp[tensor_int]);
      tensor = vecvx->Get("tensor");
      ParamTools::AsTensor<double>(tensor->Get("real"),dimbas, 1, mat);
      mat.GetCol(vec,0);
      mod_red_vectors_[ 4*tensor_int +2] = vec;

      PtrParamNode vecvy = vecnode->Get("vecvy" + tensor_comp[tensor_int]);
      tensor = vecvy->Get("tensor");
      ParamTools::AsTensor<double>(tensor->Get("real"),dimbas, 1, mat);
      mat.GetCol(vec,0);
      mod_red_vectors_[ 4*tensor_int +3] = vec;

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
  case DENSITY_TIMES_2D_TENSOR_CONSTANT_TRACE:
  case DENSITY_TIMES_ROT_TRANSVERSAL_ISOTROPIC_BOXED:
    return r+6;
  case DENSITY_TIMES_2D_TENSOR:
  case DENSITY_TIMES_ROTATED_2D_TENSOR:
    return r+7;
  case REDBAS_PARAM:
  case REDBAS_FREE:
  case GREEDY_PARAM:
  case GREEDY_FREE:
	  return r+4;
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
  case REDBAS_PARAM:
  case GREEDY_PARAM:
    return(design.Find(DesignElement::ROTANGLE) >= 0
        && design.Find(DesignElement::ROTANGLE2) >= 0
        && design.Find(DesignElement::SCALING1) >= 0
        && design.Find(DesignElement::SCALING2) >= 0);
  case REDBAS_FREE:
  case GREEDY_FREE:
      return(design.Find(DesignElement::G11) >= 0
          && design.Find(DesignElement::G12) >= 0
          && design.Find(DesignElement::G21) >= 0
          && design.Find(DesignElement::G22) >= 0);
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

void DesignMaterial::GetHomRectTensor(Matrix<double>& E, DesignElement::Type direction, Notation notation)
{
   Quad9FE fe;

   double a = params_[DesignElement::STIFF1];
   double b = params_[DesignElement::STIFF2];
   double rotAngle = params_[DesignElement::ROTANGLE];

   Vector<double> p(2);
   p[0] = -1.0 + 4 * a; // assume max 0.5
   p[1] = -1.0 + 4 * b; // assume max 0.5

   LOG_DBG2(dm) << "GHRT: dir=" << (direction == DesignElement::NO_DERIVATIVE ? "no_derivative" : DesignElement::type.ToString(direction))
                << " not=" << notation << " rotAngle=" << rotAngle << " a=" << a << " b=" << b << " -> " << p.ToString();

   switch(direction)
   {
   case DesignElement::NO_DERIVATIVE:
   case DesignElement::ROTANGLE:
   {
     Vector<double> shape;
     fe.GetShFnc(shape, p, NULL);
     ApplyHomRectTensor(E, shape);

     LOG_DBG2(dm) << "GHRT: shape=" << shape.ToString();
     break;
   }

   case DesignElement::STIFF1:
   case DesignElement::STIFF2:
   {
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

     break;
   }

   default:
     ZeroTensor(E, PLANE);
   }

   LOG_DBG2(dm) << "GHRT: E before rotation = " << E.ToString(2);
   RotateHMStiffnessTensor(E, PLANE, direction, rotAngle, notation);

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

void DesignMaterial::GetRedBasCorrector(StdVector<Vector<double> >& corrector_, const Matrix<double>& G)
{
   corrector_ = StdVector<Vector<double> >(3);
  //std::cout << "corrector"<< std::endl;
    double G11 = G[1-1][1-1];
    double G12 = G[1-1][2-1];
    double G21 = G[2-1][1-1];
    double G22 = G[2-1][2-1];

   //std::cout << "G11=" << G11 << " G12=" << G12 << " G21=" << G21 << " G22=" << G22;

  Matrix<double> mat(dimension_, dimension_);
  mat.Init();

  //mat = G11*G11*matuxux11 + G11*G12*(matuxuy11 + matuxuy11t) + G12*G12*matuyuy11
  //    + G21*G21*matvxvx11 + G21*G22*(matvxvy11 + matvxvy11t) + G22*G22*matvyvy11
  //    + G11*G21*(matuxvx12 + matuxvx12t) + G12*G21*(matuyvx12 + matuyvx12t) + G11*G22*(matuxvy12 + matuxvy12t) + G12*G22*(matuyvy12 + matuyvy12t)
  //    +0.25*(G11*G11*matvxvx33 + G11*G12*(matvxvy33 + matvxvy33t) + G12*G12*matvyvy33)
  //    +0.25*(G21*G21*matuxux33 + G21*G22*(matuxuy33 + matuxuy33t) + G22*G22*matuyuy33)
  //    +0.25*( G11*G21*(matuxvx33 + matuxvx33t) + G12*G21*(matuxvy33 + matuxvy33t) + G11*G22*(matuyvx33 + matuyvx33t) + G12*G22*(matuyvy33 + matuyvy33t))

  /*mat = mod_red_matrices_[10*0 + 0]*G11*G11
      +(mod_red_matrices_[10*0 + 3] + Transpose(mod_red_matrices_[10*0 + 3]))*(G11*G12)
      + mod_red_matrices_[10*0 + 7]*G12*G12

      + mod_red_matrices_[10*0 + 2]*G21*G21
      +(mod_red_matrices_[10*0 + 6]+ Transpose(mod_red_matrices_[10*0 + 6]))*G21*G22
      + mod_red_matrices_[10*0 + 9]*G22*G22


      + (mod_red_matrices_[10*1 + 1] + Transpose(mod_red_matrices_[10*1 + 1]))*G11*G21
      + (mod_red_matrices_[10*1 + 5] + Transpose(mod_red_matrices_[10*1 + 5]))*G12*G21
      + (mod_red_matrices_[10*1 + 4] + Transpose(mod_red_matrices_[10*1 + 4]))*G11*G22
      + (mod_red_matrices_[10*1 + 8] + Transpose(mod_red_matrices_[10*1 + 8]))*G12*G22



      + mod_red_matrices_[10*2 + 0]*G21*G21
      +(mod_red_matrices_[10*2 + 3] + Transpose(mod_red_matrices_[10*2 + 3]))*G21*G22
      + mod_red_matrices_[10*2 + 7]*G22*G22

      + mod_red_matrices_[10*2 + 2]*G11*G11
      +(mod_red_matrices_[10*2 + 6]+ Transpose(mod_red_matrices_[10*2 + 6]))*G11*G12
      + mod_red_matrices_[10*2 + 9]*G12*G12

      + (mod_red_matrices_[10*2 + 1] + Transpose(mod_red_matrices_[10*2 + 1]))*G11*G21
      + (mod_red_matrices_[10*2 + 5] + Transpose(mod_red_matrices_[10*2 + 5]))*G11*G22
      + (mod_red_matrices_[10*2 + 4] + Transpose(mod_red_matrices_[10*2 + 4]))*G12*G21
      + (mod_red_matrices_[10*2 + 8] + Transpose(mod_red_matrices_[10*2 + 8]))*G12*G22;*/


  mat.Add(G11*G11, mod_red_matrices_[10*0 + 0]);
   mat.Add(G11*G12, mod_red_matrices_[10*0 + 3]);
   mat.AddT(G11*G12, mod_red_matrices_[10*0 + 3]);
   mat.Add(G12*G12, mod_red_matrices_[10*0 + 7]);

   mat.Add(G21*G21, mod_red_matrices_[10*0 + 2]);
   mat.Add(G21*G22 ,mod_red_matrices_[10*0 + 6]);
   mat.AddT(G21*G22,mod_red_matrices_[10*0 + 6]);
   mat.Add(G22*G22, mod_red_matrices_[10*0 + 9]);

   mat.Add(G11*G21, mod_red_matrices_[10*1 + 1] );
   mat.AddT(G11*G21, mod_red_matrices_[10*1 + 1] );
   mat.Add(G12*G21, mod_red_matrices_[10*1 + 5]);
   mat.AddT(G12*G21, mod_red_matrices_[10*1 + 5]);
   mat.Add(G11*G22, mod_red_matrices_[10*1 + 4]);
   mat.AddT(G11*G22, mod_red_matrices_[10*1 + 4]);
   mat.Add(G12*G22, mod_red_matrices_[10*1 + 8]);
   mat.AddT(G12*G22, mod_red_matrices_[10*1 + 8]);

   mat.Add(G21*G21, mod_red_matrices_[10*2 + 0]);
   mat.Add(G21*G22, mod_red_matrices_[10*2 + 3] );
   mat.AddT(G21*G22, mod_red_matrices_[10*2 + 3] );
   mat.Add(G22*G22,mod_red_matrices_[10*2 + 7] );

   mat.Add(G11*G11,mod_red_matrices_[10*2 + 2] );
   mat.Add(G11*G12,mod_red_matrices_[10*2 + 6] );
   mat.AddT(G11*G12,mod_red_matrices_[10*2 + 6] );
   mat.Add(G12*G12, mod_red_matrices_[10*2 + 9]);

   mat.Add(G11*G21, mod_red_matrices_[10*2 + 1]);
   mat.AddT(G11*G21, mod_red_matrices_[10*2 + 1]);
   mat.Add(G11*G22, mod_red_matrices_[10*2 + 5]);
   mat.AddT(G11*G22, mod_red_matrices_[10*2 + 5]);
   mat.Add(G12*G21, mod_red_matrices_[10*2 + 4]);
   mat.AddT(G12*G21, mod_red_matrices_[10*2 + 4]);
   mat.Add(G12*G22, mod_red_matrices_[10*2 + 8]);
   mat.AddT(G12*G22, mod_red_matrices_[10*2 + 8]);


  for (int corrector_int =0; corrector_int<3; corrector_int++)
  {

    Vector<double> vec(dimension_);

        if (corrector_int == 0)
        {
          //vec = G11*vecux11 + G12vecuy11 + G21*vecvx12 + G22*vecvy12
          vec = -(mod_red_vectors_[4*0 +0]*G11 + mod_red_vectors_[4*0 +1]*G12
              +mod_red_vectors_[4*1 +2]*G21 + mod_red_vectors_[4*1 +3]*G22);
        }
        else if (corrector_int == 1)
        {
          //vec = G21*vecvx11 + G22*vecvy11 + G11*vecux12 + G12*vecuy12
          vec = -(mod_red_vectors_[4*0 +2]*G21 + mod_red_vectors_[4*0 +3]*G22
               +mod_red_vectors_[4*1 +0]*G11 + mod_red_vectors_[4*1 +1]*G12);
        }
        else if (corrector_int == 2)
        {
          //vec=0.5*(G21*vecux33 + G22*vecuy33 + G11*vecvx33 + G12*vecvy33)
          vec = -(mod_red_vectors_[ 4*2 +0]*G21 + mod_red_vectors_[4*2 +1]*G22
                +mod_red_vectors_[4*2 +2]*G11 + mod_red_vectors_[4*2 +3]*G12);
        }

        //lapackSysMatType lmst = ZHESV;
        //lapackSysMatType lmst =ZSYSV;
        //Matrix<Complex> matcomp(dimension_, dimension_);
        //Matrix<Complex> veccomp(dimension_,1);
        //matcomp.SetPart(Global::REAL,mat);
        //veccomp.SetPart(Global::REAL,vec);
        //matcomp.solveWithLapack(veccomp,lmst);
        //corrector_[corrector_int] = veccomp.GetPart(Global::REAL);

        Vector<double> sol(dimension_);
        mat.DirectSolve(sol,vec);
        corrector_[corrector_int] = sol;
  }
}

void DesignMaterial::GetModRedHomTensor(Matrix<double>& E, const Matrix<double>& G, const StdVector<Vector<double> >& corrector_, Notation notation)
{

  E.Resize(3,3);
  E.Init();

 // std::cout << "GetModRedTensorHom" << std::endl;
  double G11 = G(0,0);
  double G12 = G(0,1);
  double G21 = G(1,0);
  double G22 = G(1,1);

  UInt dimbas =0;
  if ((type_ == REDBAS_PARAM) | (type_ == REDBAS_FREE)) dimbas = dimension_;
  if ((type_ == GREEDY_PARAM) | (type_ == GREEDY_FREE)) dimbas = 3*dimension_;

  Matrix<double> mat(dimbas, dimbas);
  mat.Init();

  /*mat = mod_red_matrices_[10*0 + 0]*G11*G11
      +(mod_red_matrices_[10*0 + 3] + Transpose(mod_red_matrices_[10*0 + 3]))*(G11*G12)
      + mod_red_matrices_[10*0 + 7]*G12*G12

      + mod_red_matrices_[10*0 + 2]*G21*G21
      +(mod_red_matrices_[10*0 + 6]+ Transpose(mod_red_matrices_[10*0 + 6]))*G21*G22
      + mod_red_matrices_[10*0 + 9]*G22*G22


      + (mod_red_matrices_[10*1 + 1] + Transpose(mod_red_matrices_[10*1 + 1]))*G11*G21
      + (mod_red_matrices_[10*1 + 5] + Transpose(mod_red_matrices_[10*1 + 5]))*G12*G21
      + (mod_red_matrices_[10*1 + 4] + Transpose(mod_red_matrices_[10*1 + 4]))*G11*G22
      + (mod_red_matrices_[10*1 + 8] + Transpose(mod_red_matrices_[10*1 + 8]))*G12*G22



      + mod_red_matrices_[10*2 + 0]*G21*G21
      +(mod_red_matrices_[10*2 + 3] + Transpose(mod_red_matrices_[10*2 + 3]))*G21*G22
      + mod_red_matrices_[10*2 + 7]*G22*G22

      + mod_red_matrices_[10*2 + 2]*G11*G11
      +(mod_red_matrices_[10*2 + 6]+ Transpose(mod_red_matrices_[10*2 + 6]))*G11*G12
      + mod_red_matrices_[10*2 + 9]*G12*G12

      + (mod_red_matrices_[10*2 + 1] + Transpose(mod_red_matrices_[10*2 + 1]))*G11*G21
      + (mod_red_matrices_[10*2 + 5] + Transpose(mod_red_matrices_[10*2 + 5]))*G11*G22
      + (mod_red_matrices_[10*2 + 4] + Transpose(mod_red_matrices_[10*2 + 4]))*G12*G21
      + (mod_red_matrices_[10*2 + 8] + Transpose(mod_red_matrices_[10*2 + 8]))*G12*G22;*/


  mat.Add(G11*G11, mod_red_matrices_[10*0 + 0]);
  mat.Add(G11*G12, mod_red_matrices_[10*0 + 3]);
  mat.AddT(G11*G12, mod_red_matrices_[10*0 + 3]);
  mat.Add(G12*G12, mod_red_matrices_[10*0 + 7]);

  mat.Add(G21*G21, mod_red_matrices_[10*0 + 2]);
  mat.Add(G21*G22 ,mod_red_matrices_[10*0 + 6]);
  mat.AddT(G21*G22,mod_red_matrices_[10*0 + 6]);
  mat.Add(G22*G22, mod_red_matrices_[10*0 + 9]);

  mat.Add(G11*G21, mod_red_matrices_[10*1 + 1] );
  mat.AddT(G11*G21, mod_red_matrices_[10*1 + 1] );
  mat.Add(G12*G21, mod_red_matrices_[10*1 + 5]);
  mat.AddT(G12*G21, mod_red_matrices_[10*1 + 5]);
  mat.Add(G11*G22, mod_red_matrices_[10*1 + 4]);
  mat.AddT(G11*G22, mod_red_matrices_[10*1 + 4]);
  mat.Add(G12*G22, mod_red_matrices_[10*1 + 8]);
  mat.AddT(G12*G22, mod_red_matrices_[10*1 + 8]);

  mat.Add(G21*G21, mod_red_matrices_[10*2 + 0]);
  mat.Add(G21*G22, mod_red_matrices_[10*2 + 3] );
  mat.AddT(G21*G22, mod_red_matrices_[10*2 + 3] );
  mat.Add(G22*G22,mod_red_matrices_[10*2 + 7] );

  mat.Add(G11*G11,mod_red_matrices_[10*2 + 2] );
  mat.Add(G11*G12,mod_red_matrices_[10*2 + 6] );
  mat.AddT(G11*G12,mod_red_matrices_[10*2 + 6] );
  mat.Add(G12*G12, mod_red_matrices_[10*2 + 9]);

  mat.Add(G11*G21, mod_red_matrices_[10*2 + 1]);
  mat.AddT(G11*G21, mod_red_matrices_[10*2 + 1]);
  mat.Add(G11*G22, mod_red_matrices_[10*2 + 5]);
  mat.AddT(G11*G22, mod_red_matrices_[10*2 + 5]);
  mat.Add(G12*G21, mod_red_matrices_[10*2 + 4]);
  mat.AddT(G12*G21, mod_red_matrices_[10*2 + 4]);
  mat.Add(G12*G22, mod_red_matrices_[10*2 + 8]);
  mat.AddT(G12*G22, mod_red_matrices_[10*2 + 8]);


  Vector<double> mat0(dimbas);
  Vector<double> mat1(dimbas);
  Vector<double> mat2(dimbas);

  mat.Mult(corrector_[0], mat0);
  mat.Mult(corrector_[1], mat1);
  mat.Mult(corrector_[2], mat2);

  Vector<double> vec0(dimbas);
  Vector<double> vec1(dimbas);
  Vector<double> vec2(dimbas);

  vec0 = mod_red_vectors_[4*0+0]*G11 + mod_red_vectors_[4*0+1]*G12 + mod_red_vectors_[4*1+2]*G21 + mod_red_vectors_[4*1+3]*G22;
  vec1 = mod_red_vectors_[4*1+0]*G11 + mod_red_vectors_[4*1+1]*G12 + mod_red_vectors_[4*0+2]*G21 + mod_red_vectors_[4*0+3]*G22;
  vec2 = mod_red_vectors_[4*2+0]*G21 + mod_red_vectors_[4*2+1]*G22 + mod_red_vectors_[4*2+2]*G11 + mod_red_vectors_[4*2+3]*G12;


 // double p11 = corrector_[0]*vec1;
 // double p12 = corrector_[0]*vec2;
 // double p13 = corrector_[0]*vec3;

 // double p21 = corrector_[1]*vec1;
 // double p22 = corrector_[1]*vec2;
 // double p23 = corrector_[1]*vec3;

 // double p31 = corrector_[2]*vec1;
 // double p32 = corrector_[2]*vec2;
 // double p33 = corrector_[2]*vec3;

   //E11 = E11mean + E11*dxg(u1) + E12*dyg(v1)
   E[1-1][1-1] = mean_tensor_[1-1][0]

               + corrector_[0].Inner(mat0) + 2.0*(corrector_[0].Inner(vec0));

                /*+ G11*(corrector_[0]*mod_red_vectors_[4*0+0])
                + G12*(corrector_[0]*mod_red_vectors_[4*0+1])
                + G21*(corrector_[0]*mod_red_vectors_[4*1+2])
                + G22*(corrector_[0]*mod_red_vectors_[4*1+3])

                + G11*(corrector_[0]*mod_red_vectors_[4*0+0])
                + G12*(corrector_[0]*mod_red_vectors_[4*0+1])
                + G21*(corrector_[0]*mod_red_vectors_[4*1+2])
                + G22*(corrector_[0]*mod_red_vectors_[4*1+3])

                +corrector_[0]*mat0;*/



  // E12 = E12mean - E11*dxg(u2) - E12*dyg(v2)
   E[1-1][2-1] = mean_tensor_[2-1][0]

               + corrector_[1].Inner(mat0)  + corrector_[1].Inner(vec0) + corrector_[0].Inner(vec1);



                /*+ G11*(corrector_[1]*mod_red_vectors_[4*0+0])
                + G12*(corrector_[1]*mod_red_vectors_[4*0+1])
                + G21*(corrector_[1]*mod_red_vectors_[4*1+2])
                + G22*(corrector_[1]*mod_red_vectors_[4*1+3])

                + G11*(corrector_[0]*mod_red_vectors_[4*1+0])
                + G12*(corrector_[0]*mod_red_vectors_[4*1+1])
                + G21*(corrector_[0]*mod_red_vectors_[4*0+2])
                + G22*(corrector_[0]*mod_red_vectors_[4*0+3])

                +corrector_[1]*mat0;*/

   E[2-1][1-1] =  E[1-1][2-1] ;


  //E13 = - E11*dxg(u3) - E12*dyg(v3)
   E[1-1][3-1] = 0

             + corrector_[2].Inner(mat0) + corrector_[2].Inner(vec0) + corrector_[0].Inner(vec2);

             /*+ G11*(corrector_[2]*mod_red_vectors_[4*0+0])
             + G12*(corrector_[2]*mod_red_vectors_[4*0+1])
             + G21*(corrector_[2]*mod_red_vectors_[4*1+2])
             + G22*(corrector_[2]*mod_red_vectors_[4*1+3])

             + G21*(corrector_[0]*mod_red_vectors_[4*2+0])
             + G22*(corrector_[0]*mod_red_vectors_[4*2+1])
             + G11*(corrector_[0]*mod_red_vectors_[4*2+2])
             + G12*(corrector_[0]*mod_red_vectors_[4*2+3])

             +corrector_[2]*mat0;*/

   E[3-1][1-1] = E[1-1][3-1];

  //E22 = E11mean - E12*dxg(u2) - E11*dyg(v2)
   E[2-1][2-1] = mean_tensor_[1-1][0]

              + corrector_[1].Inner(mat1) + 2.0*(corrector_[1].Inner(vec1));

               /* + G11*(corrector_[1]*mod_red_vectors_[4*1+0])
                + G12*(corrector_[1]*mod_red_vectors_[4*1+1])
                + G21*(corrector_[1]*mod_red_vectors_[4*0+2])
                + G22*(corrector_[1]*mod_red_vectors_[4*0+3])

                + G11*(corrector_[1]*mod_red_vectors_[4*1+0])
                + G12*(corrector_[1]*mod_red_vectors_[4*1+1])
                + G21*(corrector_[1]*mod_red_vectors_[4*0+2])
                + G22*(corrector_[1]*mod_red_vectors_[4*0+3])

                +corrector_[1]*mat1;*/


  //E23 = - E12*dxg(u3) - E11*dyg(v3)
   E[2-1][3-1] = 0

             + corrector_[2].Inner(mat1) + corrector_[2].Inner(vec1) + corrector_[1].Inner(vec2);

               /*+ G11*(corrector_[2]*mod_red_vectors_[4*1+0])
               + G12*(corrector_[2]*mod_red_vectors_[4*1+1])
               + G21*(corrector_[2]*mod_red_vectors_[4*0+2])
               + G22*(corrector_[2]*mod_red_vectors_[4*0+3])

               + G21*(corrector_[1]*mod_red_vectors_[4*2+0])
               + G22*(corrector_[1]*mod_red_vectors_[4*2+1])
               + G11*(corrector_[1]*mod_red_vectors_[4*2+2])
               + G12*(corrector_[1]*mod_red_vectors_[4*2+3])

               +corrector_[2]*mat1;*/

   E[3-1][2-1] = E[2-1][3-1];

  //E33 = E33mean - E33*0.5*(dyg(u3) + dxg(v3))
   E[3-1][3-1] = mean_tensor_[3-1][0]

               + corrector_[2].Inner(mat2) + 2.0*(corrector_[2].Inner(vec2));

                /*+ G21*(corrector_[2]*mod_red_vectors_[4*2+0])
                + G22*(corrector_[2]*mod_red_vectors_[4*2+1])
                + G11*(corrector_[2]*mod_red_vectors_[4*2+2])
                + G12*(corrector_[2]*mod_red_vectors_[4*2+3])

                + G21*(corrector_[2]*mod_red_vectors_[4*2+0])
                + G22*(corrector_[2]*mod_red_vectors_[4*2+1])
                + G11*(corrector_[2]*mod_red_vectors_[4*2+2])
                + G12*(corrector_[2]*mod_red_vectors_[4*2+3])

               +corrector_[2]*mat2;*/

   E[1-1][3-1] = sqrt(2)*E[1-1][3-1];
   E[3-1][1-1] = E[1-1][3-1];
   E[2-1][3-1] = sqrt(2)*E[2-1][3-1];
   E[3-1][2-1] = E[2-1][3-1];
   E[3-1][3-1] = 2*E[3-1][3-1];

   if(notation != HILL_MANDEL)
   {
     E[1-1][3-1] = 1.0/sqrt(2)*E[1-1][3-1];
     E[3-1][1-1] = E[1-1][3-1];
     E[2-1][3-1] = 1.0/sqrt(2)*E[2-1][3-1];
     E[3-1][2-1] = E[2-1][3-1];
     E[3-1][3-1] = 0.5*E[3-1][3-1];
   }
}

void DesignMaterial::GetModRedHomTensor(Matrix<double>& E, const Matrix<double>& G, const Matrix<double>& Gderiv, const StdVector<Vector<double> >& corrector_, Notation notation)
{

   //Should be done with Hill-Mandel notation
  E.Resize(3,3);
  E.Init();

  double G11 = G[1-1][1-1];
  double G12 = G[1-1][2-1];
  double G21 = G[2-1][1-1];
  double G22 = G[2-1][2-1];

  double G11d = Gderiv[1-1][1-1];
  double G12d = Gderiv[1-1][2-1];
  double G21d = Gderiv[2-1][1-1];
  double G22d = Gderiv[2-1][2-1];

  UInt dimbas =0;
  if ((type_ == REDBAS_PARAM) | (type_ == REDBAS_FREE)) dimbas = dimension_;
  if ((type_ == GREEDY_PARAM) | (type_ == GREEDY_FREE)) dimbas = 3*dimension_;

  Matrix<double> matd(dimbas, dimbas);
  matd.Init();

  /*matd= mod_red_matrices_[10*0 + 0]*(G11d*G11 + G11*G11d)
        +(mod_red_matrices_[10*0 + 3] + Transpose(mod_red_matrices_[10*0 + 3]))*(G11d*G12 + G11*G12d)
        + mod_red_matrices_[10*0 + 7]*(G12d*G12 + G12*G12d)

        + mod_red_matrices_[10*0 + 2]*(G21d*G21 + G21*G21d)
        +(mod_red_matrices_[10*0 + 6]+ Transpose(mod_red_matrices_[10*0 + 6]))*(G21d*G22 + G21*G22d)
        + mod_red_matrices_[10*0 + 9]*(G22d*G22 + G22*G22d)


        + (mod_red_matrices_[10*1 + 1] + Transpose(mod_red_matrices_[10*1 + 1]))*(G11d*G21 + G11*G21d)
        + (mod_red_matrices_[10*1 + 5] + Transpose(mod_red_matrices_[10*1 + 5]))*(G12d*G21 + G12*G21d)
        + (mod_red_matrices_[10*1 + 4] + Transpose(mod_red_matrices_[10*1 + 4]))*(G11d*G22 + G11*G22d)
        + (mod_red_matrices_[10*1 + 8] + Transpose(mod_red_matrices_[10*1 + 8]))*(G12d*G22 + G12*G22d)



        + mod_red_matrices_[10*2 + 0]*(G21d*G21 + G21*G21d)
        +(mod_red_matrices_[10*2 + 3] + Transpose(mod_red_matrices_[10*2 + 3]))*(G21d*G22 + G21*G22d)
        + mod_red_matrices_[10*2 + 7]*(G22d*G22 + G22*G22d)

        + mod_red_matrices_[10*2 + 2]*(G11d*G11 + G11*G11d)
        +(mod_red_matrices_[10*2 + 6]+ Transpose(mod_red_matrices_[10*2 + 6]))*(G11d*G12 + G11*G12d)
        + mod_red_matrices_[10*2 + 9]*(G12d*G12 + G12*G12d)

        + (mod_red_matrices_[10*2 + 1] + Transpose(mod_red_matrices_[10*2 + 1]))*(G11d*G21 + G11*G21d)
        + (mod_red_matrices_[10*2 + 5] + Transpose(mod_red_matrices_[10*2 + 5]))*(G11d*G22 + G11*G22d)
        + (mod_red_matrices_[10*2 + 4] + Transpose(mod_red_matrices_[10*2 + 4]))*(G12d*G21 + G12*G21d)
        + (mod_red_matrices_[10*2 + 8] + Transpose(mod_red_matrices_[10*2 + 8]))*(G12d*G22 + G12*G22d);*/

  matd.Add(G11d*G11 + G11*G11d, mod_red_matrices_[10*0 + 0]);
  matd.Add(G11d*G12 + G11*G12d, mod_red_matrices_[10*0 + 3]);
  matd.AddT(G11d*G12 + G11*G12d, mod_red_matrices_[10*0 + 3]);
  matd.Add(G12d*G12 + G12*G12d, mod_red_matrices_[10*0 + 7]);

  matd.Add(G21d*G21 + G21*G21d, mod_red_matrices_[10*0 + 2]);
  matd.Add(G21d*G22 + G21*G22d,mod_red_matrices_[10*0 + 6]);
  matd.AddT(G21d*G22 + G21*G22d,mod_red_matrices_[10*0 + 6]);
  matd.Add(G22d*G22 + G22*G22d, mod_red_matrices_[10*0 + 9]);

  matd.Add(G11d*G21 + G11*G21d, mod_red_matrices_[10*1 + 1] );
  matd.AddT(G11d*G21 + G11*G21d, mod_red_matrices_[10*1 + 1] );
  matd.Add(G12d*G21 + G12*G21d, mod_red_matrices_[10*1 + 5]);
  matd.AddT(G12d*G21 + G12*G21d, mod_red_matrices_[10*1 + 5]);
  matd.Add(G11d*G22 + G11*G22d, mod_red_matrices_[10*1 + 4]);
  matd.AddT(G11d*G22 + G11*G22d, mod_red_matrices_[10*1 + 4]);
  matd.Add(G12d*G22 + G12*G22d, mod_red_matrices_[10*1 + 8]);
  matd.AddT(G12d*G22 + G12*G22d, mod_red_matrices_[10*1 + 8]);

  matd.Add(G21d*G21 + G21*G21d, mod_red_matrices_[10*2 + 0]);
  matd.Add(G21d*G22 + G21*G22d, mod_red_matrices_[10*2 + 3] );
  matd.AddT(G21d*G22 + G21*G22d, mod_red_matrices_[10*2 + 3] );
  matd.Add(G22d*G22 + G22*G22d,mod_red_matrices_[10*2 + 7] );

  matd.Add(G11d*G11 + G11*G11d,mod_red_matrices_[10*2 + 2] );
  matd.Add(G11d*G12 + G11*G12d,mod_red_matrices_[10*2 + 6] );
  matd.AddT(G11d*G12 + G11*G12d,mod_red_matrices_[10*2 + 6] );
  matd.Add(G12d*G12 + G12*G12d, mod_red_matrices_[10*2 + 9]);

  matd.Add(G11d*G21 + G11*G21d, mod_red_matrices_[10*2 + 1]);
  matd.AddT(G11d*G21 + G11*G21d, mod_red_matrices_[10*2 + 1]);
  matd.Add(G11d*G22 + G11*G22d, mod_red_matrices_[10*2 + 5]);
  matd.AddT(G11d*G22 + G11*G22d, mod_red_matrices_[10*2 + 5]);
  matd.Add(G12d*G21 + G12*G21d, mod_red_matrices_[10*2 + 4]);
  matd.AddT(G12d*G21 + G12*G21d, mod_red_matrices_[10*2 + 4]);
  matd.Add(G12d*G22 + G12*G22d, mod_red_matrices_[10*2 + 8]);
  matd.AddT(G12d*G22 + G12*G22d, mod_red_matrices_[10*2 + 8]);



   Vector<double> mat0(dimbas);
   Vector<double> mat1(dimbas);
   Vector<double> mat2(dimbas);

   matd.Mult(corrector_[0], mat0);
   matd.Mult(corrector_[1], mat1);
   matd.Mult(corrector_[2], mat2);


   Vector<double> vec0(dimbas);
    Vector<double> vec1(dimbas);
    Vector<double> vec2(dimbas);

    vec0 = mod_red_vectors_[4*0+0]*G11d + mod_red_vectors_[4*0+1]*G12d + mod_red_vectors_[4*1+2]*G21d + mod_red_vectors_[4*1+3]*G22d;
    vec1 = mod_red_vectors_[4*1+0]*G11d + mod_red_vectors_[4*1+1]*G12d + mod_red_vectors_[4*0+2]*G21d + mod_red_vectors_[4*0+3]*G22d;
    vec2 = mod_red_vectors_[4*2+0]*G21d + mod_red_vectors_[4*2+1]*G22d + mod_red_vectors_[4*2+2]*G11d + mod_red_vectors_[4*2+3]*G12d;



  //E11 =- E11*dxgd(u1) - E12*dygd(v1) - E11*dxgd(u1) - E12*dygd(v1) + corr0'*mat*corr0;
  E[1-1][1-1] = 0

                /*+ G11d*(corrector_[0]*mod_red_vectors_[4*0+0])
                + G12d*(corrector_[0]*mod_red_vectors_[4*0+1])
                + G21d*(corrector_[0]*mod_red_vectors_[4*1+2])
                + G22d*(corrector_[0]*mod_red_vectors_[4*1+3])

                + G11d*(corrector_[0]*mod_red_vectors_[4*0+0])
                + G12d*(corrector_[0]*mod_red_vectors_[4*0+1])
                + G21d*(corrector_[0]*mod_red_vectors_[4*1+2])
                + G22d*(corrector_[0]*mod_red_vectors_[4*1+3])

                +corrector_[0]*mat0; */

      +corrector_[0].Inner(mat0) + 2.0*(corrector_[0].Inner(vec0));

 // E12 = - E11*dxgd(u2) - E12*dygd(v2) - E12*dxgd(u1) - E22*dygd(v1) + corr1'*mat*corr0;
  E[1-1][2-1] =  0

                 /*+ G11d*(corrector_[1]*mod_red_vectors_[4*0+0])
                 + G12d*(corrector_[1]*mod_red_vectors_[4*0+1])
                 + G21d*(corrector_[1]*mod_red_vectors_[4*1+2])
                 + G22d*(corrector_[1]*mod_red_vectors_[4*1+3])

                 + G11d*(corrector_[0]*mod_red_vectors_[4*1+0])
                 + G12d*(corrector_[0]*mod_red_vectors_[4*1+1])
                 + G21d*(corrector_[0]*mod_red_vectors_[4*0+2])
                 + G22d*(corrector_[0]*mod_red_vectors_[4*0+3])

                 +corrector_[1]*mat0;*/

               +corrector_[1].Inner(mat0) + corrector_[1].Inner(vec0) + corrector_[0].Inner(vec1);

  E[2-1][1-1] =  E[1-1][2-1] ;


 //E13 = - E11*dxgd(u3) - E12*dygd(v3) - E33*0.5*(dygd(u1)+dxgd(v1)) + corr2'*mat*corr0;
  E[1-1][3-1] =  0

                /* + G11d*(corrector_[2]*mod_red_vectors_[4*0+0])
                 + G12d*(corrector_[2]*mod_red_vectors_[4*0+1])
                 + G21d*(corrector_[2]*mod_red_vectors_[4*1+2])
                 + G22d*(corrector_[2]*mod_red_vectors_[4*1+3])

                 + G21d*(corrector_[0]*mod_red_vectors_[4*2+0])
                 + G22d*(corrector_[0]*mod_red_vectors_[4*2+1])
                 + G11d*(corrector_[0]*mod_red_vectors_[4*2+2])
                 + G12d*(corrector_[0]*mod_red_vectors_[4*2+3])

                 +corrector_[2]*mat0;*/

                + corrector_[2].Inner(mat0) + corrector_[2].Inner(vec0) + corrector_[0].Inner(vec2);
               // + corrector_[2].Inner(mat0)  + corrector_[2].Inner(vec0) + corrector_[0]*vec2;


  E[3-1][1-1] = E[1-1][3-1];

 //E22 = - E12*dxgd(u2) - E11*dygd(v2)- E12*dxgd(u2) - E11*dygd(v2) + corr1'*mat*corr1;
  E[2-1][2-1] =  0

                 /*+ G11d*(corrector_[1]*mod_red_vectors_[4*1+0])
                 + G12d*(corrector_[1]*mod_red_vectors_[4*1+1])
                 + G21d*(corrector_[1]*mod_red_vectors_[4*0+2])
                 + G22d*(corrector_[1]*mod_red_vectors_[4*0+3])

                 + G11d*(corrector_[1]*mod_red_vectors_[4*1+0])
                 + G12d*(corrector_[1]*mod_red_vectors_[4*1+1])
                 + G21d*(corrector_[1]*mod_red_vectors_[4*0+2])
                 + G22d*(corrector_[1]*mod_red_vectors_[4*0+3])

                 +corrector_[1]*mat1;*/

                + corrector_[1].Inner(mat1) + 2.0*(corrector_[1].Inner(vec1));

 //E23 = - E12*dxgd(u3) - E11*dygd(v3) - E33*05*(dygd(u2) + dxgd(v2)) + corr1'*mat*corr2
  E[2-1][3-1] =  0

               /*+ G11d*(corrector_[2]*mod_red_vectors_[4*1+0])
               + G12d*(corrector_[2]*mod_red_vectors_[4*1+1])
               + G21d*(corrector_[2]*mod_red_vectors_[4*0+2])
               + G22d*(corrector_[2]*mod_red_vectors_[4*0+3])

               + G21d*(corrector_[1]*mod_red_vectors_[4*2+0])
               + G22d*(corrector_[1]*mod_red_vectors_[4*2+1])
               + G11d*(corrector_[1]*mod_red_vectors_[4*2+2])
               + G12d*(corrector_[1]*mod_red_vectors_[4*2+3])

               +corrector_[2]*mat1;*/

              + corrector_[2].Inner(mat1) + corrector_[2].Inner(vec1) + corrector_[1].Inner(vec2);


  E[3-1][2-1] = E[2-1][3-1];

 //E33 = - E33*0.5*(dygd(u3) + dxgd(v3))- E33*0.5*(dygd(u3) + dxgd(v3)) + corr2'*mat*corr2;
  E[3-1][3-1] = 0

            /*  + G21d*corrector_[2].FrobeniusProduct(mod_red_vectors_[4*2+0])
              + G22d*corrector_[2].FrobeniusProduct(mod_red_vectors_[4*2+1])
              + G11d*corrector_[2].FrobeniusProduct(mod_red_vectors_[4*2+2])
              + G12d*corrector_[2].FrobeniusProduct(mod_red_vectors_[4*2+3])

              + G21d*corrector_[2].FrobeniusProduct(mod_red_vectors_[4*2+0])
              + G22d*corrector_[2].FrobeniusProduct(mod_red_vectors_[4*2+1])
              + G11d*corrector_[2].FrobeniusProduct(mod_red_vectors_[4*2+2])
              + G12d*corrector_[2].FrobeniusProduct(mod_red_vectors_[4*2+3])

              +corrector_[2].FrobeniusProduct(mat2);*/

             +corrector_[2].Inner(mat2) + 2.0*(corrector_[2].Inner(vec2));

                E[1-1][3-1] = sqrt(2)*E[1-1][3-1];
                E[3-1][1-1] = E[1-1][3-1];
                E[2-1][3-1] = sqrt(2)*E[2-1][3-1];
                E[3-1][2-1] = E[2-1][3-1];
                E[3-1][3-1] = 2*E[3-1][3-1];

          if (notation != HILL_MANDEL)
          {
              E[1-1][3-1] = 1.0/sqrt(2)*E[1-1][3-1];
              E[3-1][1-1] = E[1-1][3-1];
              E[2-1][3-1] = 1.0/sqrt(2)*E[2-1][3-1];
              E[3-1][2-1] = E[2-1][3-1];
              E[3-1][3-1] = 0.5*E[3-1][3-1];
          }



}

void DesignMaterial::GetModRedGTensor(Matrix<double>& G, DesignElement::Type direction, const bool& all_param)
{

  G.Resize(2,2);
  G.Init();

  assert((type_== REDBAS_PARAM) | (type_== REDBAS_FREE) | (type_== GREEDY_PARAM) | (type_== GREEDY_FREE));

  if((type_ == REDBAS_PARAM) | (type_ == GREEDY_PARAM))
  {
  //Quad9FE fe;
  // std::cout << "GetModRedTensor" << std::endl;

  double theta = params_[DesignElement::ROTANGLE];
  double phi = params_[DesignElement::ROTANGLE2];
  double l1 = params_[DesignElement::SCALING1];
  double l2 = params_[DesignElement::SCALING2];

  switch(direction)
  {
     case DesignElement::NO_DERIVATIVE:
     {
       if (!all_param)
       {
          G(0,0) = l1*cos(phi); //matrix element G11
          G(0,1) = l1*sin(phi); //matrix element G12
          G(1,0) = -l2*sin(phi); //matrix element G21
          G(1,1) = l2*cos(phi); //matrix element G22
       }
       else
       {
           G(0,0) = l1*cos(theta)*cos(phi) -l2*sin(theta)*sin(phi); //matrix element G11
           G(0,1) = l1*cos(theta)*sin(phi) +l2*sin(theta)*cos(phi); //matrix element G12
           G(1,0) = -l1*sin(theta)*cos(phi) -l2*cos(theta)*sin(phi); //matrix element G21
           G(1,1) = -l1*sin(theta)*sin(phi) + l2*cos(theta)*cos(phi); //matrix element G22
       }
       break;
     }

     case DesignElement::ROTANGLE:
     {
         if (!all_param)
         {
           throw Exception("GetModRedGTensor should not be called when all_param_ = false and Design::Element direction = ROTANGLE");
         }
         else
         {
           G(0,0) = -l1*sin(theta)*cos(phi) -l2*cos(theta)*sin(phi); //matrix element G11
           G(0,1) = -l1*sin(theta)*sin(phi) +l2*cos(theta)*cos(phi); //matrix element G12
           G(1,0) = -l1*cos(theta)*cos(phi) +l2*sin(theta)*sin(phi); //matrix element G21
           G(1,1) = -l1*cos(theta)*sin(phi) - l2*sin(theta)*cos(phi); //matrix element G22
         }
       break;
     }

     case DesignElement::ROTANGLE2:
     {
       if (!all_param)
       {
         G(0,0) = -l1*sin(phi); //matrix element G11
         G(0,1) = l1*cos(phi); //matrix element G12
         G(1,0) = - l2*cos(phi); //matrix element G21
         G(1,1) = - l2*sin(phi); //matrix element G22
       }
       else
       {
         G(0,0) = -l1*cos(theta)*sin(phi) - l2*sin(theta)*cos(phi); //matrix element G11
         G(0,1) = l1*cos(theta)*cos(phi) - l2*sin(theta)*sin(phi); //matrix element G12
         G(1,0) = l1*sin(theta)*sin(phi) - l2*cos(theta)*cos(phi); //matrix element G21
         G(1,1) = -l1*sin(theta)*cos(phi) - l2*cos(theta)*sin(phi); //matrix element G22
       }
       break;
     }

     case DesignElement::SCALING1:
     {
       if (!all_param)
       {
         G(0,0) = cos(phi); //matrix element G11
         G(0,1) = sin(phi); //matrix element G12
         G(1,0) = 0.0; //matrix element G21
         G(1,1) = 0.0; //matrix element G22
       }
       else
       {
         G(0,0) = cos(theta)*cos(phi); //matrix element G11
         G(0,1) = cos(theta)*sin(phi); //matrix element G12
         G(1,0) = -sin(theta)*cos(phi); //matrix element G21
         G(1,1) = -sin(theta)*sin(phi); //matrix element G22
       }
       break;
     }

     case DesignElement::SCALING2:
     {
       if (!all_param)
       {
         G(0,0) =0.0;
         G(0,1) =0.0;
         G(1,0) =-sin(phi);
         G(1,1) = cos(phi);
       }
       else
       {
         G(0,0) =-sin(theta)*sin(phi);
         G(0,1) =sin(theta)*cos(phi);
         G(1,0) =-cos(theta)*sin(phi);
         G(1,1) = cos(theta)*cos(phi);
       }
       break;
     }

     default: //Zero matrix
       G(0,0) =0.0;
       G(0,1) =0.0;
       G(1,0) =0.0;
       G(1,1) =0.0;
     break;
     }

  }
  else if ((type_ == REDBAS_FREE) | (type_ == GREEDY_FREE))
  {
     double G11 = params_[DesignElement::G11];
     double G12 = params_[DesignElement::G12];
     double G21 = params_[DesignElement::G21];
     double G22 = params_[DesignElement::G22];

     switch(direction)
       {
          case DesignElement::NO_DERIVATIVE:
          {
             G(0,0) = G11; //matrix element G11
             G(0,1) = G12; //matrix element G12
             G(1,0) = G21; //matrix element G21
             G(1,1) = G22; //matrix element G22

            break;
          }
          case DesignElement::G11:
          {
              G(0,0) = 1.0; //matrix element G11
              G(0,1) = 0.0; //matrix element G12
              G(1,0) = 0.0; //matrix element G21
              G(1,1) = 0.0; //matrix element G22
            break;
          }
          case DesignElement::G12:
          {
             G(0,0) = 0.0; //matrix element G11
             G(0,1) = 1.0; //matrix element G12
             G(1,0) = 0.0; //matrix element G21
             G(1,1) = 0.0; //matrix element G22
            break;
          }
          case DesignElement::G21:
          {

            G(0,0) = 0.0; //matrix element G11
            G(0,1) = 0.0; //matrix element G12
            G(1,0) = 1.0; //matrix element G21
            G(1,1) = 0.0; //matrix element G22

            break;
          }
          case DesignElement::G22:
          {
            G(0,0) =0.0;
            G(0,1) =0.0;
            G(1,0) =0.0;
            G(1,1) =1.0;
            break;
          }
          default: //Identity matrix
            G(0,0) =0.0;
            G(0,1) =0.0;
            G(1,0) =0.0;
            G(1,1) =0.0;
            break;

          }
  }


}

void DesignMaterial::GetModRedTensor(Matrix<double>& E, DesignElement::Type direction, Notation notation)
{

  //if type_ == REDBAS_FREE or type_ == GREEDY_FREE then we must have all_param_
  assert((all_param_) | (type_ == REDBAS_PARAM) | (type_ == GREEDY_PARAM));

  E.Resize(3,3);
  E.Init();

  Matrix<double> G(2,2);
  GetModRedGTensor(G,DesignElement::NO_DERIVATIVE, all_param_);

  Vector<double> paramvec(4);
  GetModRedParamVector(paramvec);

  double theta = paramvec[0];

  StdVector<Vector<double> > corrector_(3);
  //The corrector calculation is always done with Voigt notation
  if ((type_ == REDBAS_PARAM) | (type_== REDBAS_FREE)) GetRedBasCorrector(corrector_, G);
  else if ((type_ == GREEDY_PARAM) | (type_== GREEDY_FREE)) GetGreedyCorrector(corrector_, paramvec, all_param_);

  if ((direction == DesignElement::NO_DERIVATIVE) | ((!all_param_) & (direction == DesignElement::ROTANGLE)))
  {
    GetModRedHomTensor(E, G,corrector_, notation);
  }
  else
  {
    Matrix<double> Gderiv(2,2);
    GetModRedGTensor(Gderiv,direction, all_param_);
    GetModRedHomTensor(E, G, Gderiv, corrector_, notation);
  }
  if (!all_param_) RotateHMStiffnessTensor(E, PLANE, direction, theta, notation);


}


void DesignMaterial::GetModRedParamVector(Vector<double>& params)
{

  params.Resize(4,1);
  params.Init();

  if ((type_==GREEDY_PARAM) | (type_ == REDBAS_PARAM))
  {
    //Quad9FE fe;
    // std::cout << "GetModRedTensor" << std::endl;

    double theta = params_[DesignElement::ROTANGLE];
    double phi = params_[DesignElement::ROTANGLE2];
    double l1 = params_[DesignElement::SCALING1];
    double l2 = params_[DesignElement::SCALING2];

    params[0]=theta;
    params[1]=phi;
    params[2]=l1;
    params[3]=l2;

  }
  else if ((type_==GREEDY_FREE) | (type_ == REDBAS_FREE))
  {
      Matrix<double> G(2,2);
      GetModRedGTensor(G,DesignElement::NO_DERIVATIVE, true);

      GetSVDGTensorParameters(G, params);

  }
}


/*void DesignMaterial::GetGreedyTensor(Matrix<double>& E, DesignElement::Type direction, Notation notation)
{

  E.Resize(3,3);
  E.Init();

  Matrix<double> params(4,1);
  GetGreedyParams(params);

  GetGreedyTensor(E,params,all_param_, direction,notation);

}*/



void DesignMaterial::GetGreedyCorrector(StdVector<Vector<double> >& corrector_, const Vector<double>& params, const bool& all_param)
{

  assert((type_== GREEDY_PARAM) | (type_ == GREEDY_FREE));

  double theta = params[0];
  double phi = params[1];
  double l1 = params[2];
  double l2 = params[3];

  if (l1 < lmin_) l1 = lmin_;
  else if (l1 > lmax_) l1 = lmax_;

  if (l2 < lmin_) l2 = lmin_;
  else if (l2 > lmax_) l2 = lmax_;

  corrector_ = StdVector<Vector<double> >(3);
  corrector_[0] = Vector<double>(3*dimension_);
  corrector_[1] = Vector<double>(3*dimension_);
  corrector_[2] = Vector<double>(3*dimension_);


  // corrector_[0] = Matrix<double>(dimension_, 1);
  // corrector_[1] = Matrix<double>(dimension_, 1);
  // corrector_[2] = Matrix<double>(dimension_, 1);

  //We begin by computing the coefficients of the three correctors in the reduced basis associated to the value of these parameters
  if (all_param)
  {

      for (UInt k=0; k<dimension_; k++)
      {
        Vector<double> coord1theta(2*Na_+1);
        Vector<double> coord2theta(2*Na_+1);
        Vector<double> coord3theta(2*Na_+1);
        Vector<double> coord1phi(2*Na_+1);
        Vector<double> coord2phi(2*Na_+1);
        Vector<double> coord3phi(2*Na_+1);
        Vector<double> coord1l1(Nl_+1);
        Vector<double> coord2l1(Nl_+1);
        Vector<double> coord3l1(Nl_+1);
        Vector<double> coord1l2(Nl_+1);
        Vector<double> coord2l2(Nl_+1);
        Vector<double> coord3l2(Nl_+1);

        matrices_param_[0].GetCol(coord1theta, k);
        matrices_param_[0].GetCol(coord2theta, dimension_+ k);
        matrices_param_[0].GetCol(coord3theta, 2*dimension_+ k);

        matrices_param_[1].GetCol(coord1phi, k);
        matrices_param_[1].GetCol(coord2phi, dimension_+ k);
        matrices_param_[1].GetCol(coord3phi, 2*dimension_+ k);

        matrices_param_[2].GetCol(coord1l1, k);
        matrices_param_[2].GetCol(coord2l1, dimension_+ k);
        matrices_param_[2].GetCol(coord3l1, 2*dimension_+ k);

        matrices_param_[3].GetCol(coord1l2, k);
        matrices_param_[3].GetCol(coord2l2, dimension_+ k);
        matrices_param_[3].GetCol(coord3l2, 2*dimension_+ k);

        //corrector_[0][k] = -AngleGreedyCalculus(coord1theta, theta)*AngleGreedyCalculus(coord1phi, phi)*ScalingGreedyCalculus(coord1l1,l1)*ScalingGreedyCalculus(coord1l2,l2);
        //corrector_[1][k] = -AngleGreedyCalculus(coord2theta, theta)*AngleGreedyCalculus(coord2phi, phi)*ScalingGreedyCalculus(coord2l1,l1)*ScalingGreedyCalculus(coord2l2,l2);
        //corrector_[2][k] = -AngleGreedyCalculus(coord3theta, theta)*AngleGreedyCalculus(coord3phi, phi)*ScalingGreedyCalculus(coord3l1,l1)*ScalingGreedyCalculus(coord3l2,l2);

        corrector_[0][k] = -AngleGreedyCalculus(coord1theta, theta)*AngleGreedyCalculus(coord1phi, phi)*ScalingGreedyCalculus(coord1l1,l1)*ScalingGreedyCalculus(coord1l2,l2);
        corrector_[1][dimension_ +k] = -AngleGreedyCalculus(coord2theta, theta)*AngleGreedyCalculus(coord2phi, phi)*ScalingGreedyCalculus(coord2l1,l1)*ScalingGreedyCalculus(coord2l2,l2);
        corrector_[2][2*dimension_+k] = -AngleGreedyCalculus(coord3theta, theta)*AngleGreedyCalculus(coord3phi, phi)*ScalingGreedyCalculus(coord3l1,l1)*ScalingGreedyCalculus(coord3l2,l2);


      }


  }
  else
  {
          for (UInt k=0; k<dimension_; k++)
          {
            Vector<double> coord1phi(2*Na_+1);
            Vector<double> coord2phi(2*Na_+1);
            Vector<double> coord3phi(2*Na_+1);
            Vector<double> coord1l1(Nl_+1);
            Vector<double> coord2l1(Nl_+1);
            Vector<double> coord3l1(Nl_+1);
            Vector<double> coord1l2(Nl_+1);
            Vector<double> coord2l2(Nl_+1);
            Vector<double> coord3l2(Nl_+1);

            matrices_param_[0].GetCol(coord1phi, k);
            matrices_param_[0].GetCol(coord2phi, dimension_+ k);
            matrices_param_[0].GetCol(coord3phi, 2*dimension_+ k);

            matrices_param_[1].GetCol(coord1l1, k);
            matrices_param_[1].GetCol(coord2l1, dimension_+ k);
            matrices_param_[1].GetCol(coord3l1, 2*dimension_+ k);

            matrices_param_[2].GetCol(coord1l2, k);
            matrices_param_[2].GetCol(coord2l2, dimension_+ k);
            matrices_param_[2].GetCol(coord3l2, 2*dimension_+ k);

            //corrector_[0][k] = -AngleGreedyCalculus(coord1phi, phi)*ScalingGreedyCalculus(coord1l1,l1)*ScalingGreedyCalculus(coord1l2,l2);
            //corrector_[1][k] = -AngleGreedyCalculus(coord2phi, phi)*ScalingGreedyCalculus(coord2l1,l1)*ScalingGreedyCalculus(coord2l2,l2);
            //corrector_[2][k] = -AngleGreedyCalculus(coord3phi, phi)*ScalingGreedyCalculus(coord3l1,l1)*ScalingGreedyCalculus(coord3l2,l2);

            corrector_[0][k] = -AngleGreedyCalculus(coord1phi, phi)*ScalingGreedyCalculus(coord1l1,l1)*ScalingGreedyCalculus(coord1l2,l2);
            corrector_[1][dimension_ +k] = -AngleGreedyCalculus(coord2phi, phi)*ScalingGreedyCalculus(coord2l1,l1)*ScalingGreedyCalculus(coord2l2,l2);
            corrector_[2][2*dimension_+k] = -AngleGreedyCalculus(coord3phi, phi)*ScalingGreedyCalculus(coord3l1,l1)*ScalingGreedyCalculus(coord3l2,l2);

          }

  }
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

void DesignMaterial::RotateHMStiffnessTensor(Matrix<double>& t, SubTensorType subTensor, DesignElement::Type direction, double a, Notation notation){
  switch(subTensor){
  case PLANE_STRAIN:
  case PLANE_STRESS:
  case PLANE:
  {
    Matrix<double> theta(3,3);
     Matrix<double> help(3,3);
     const double sq2inv = 1/sqrt(2);
     theta.SetEntry(0,0, pow(cos(a),2));
     theta.SetEntry(0,1, pow(sin(a),2));
     theta.SetEntry(0,2, -sqrt(2)/2*sin(2*a));
     theta.SetEntry(1,0, theta(0,1));
     theta.SetEntry(1,1, theta(0,0));
     theta.SetEntry(1,2, -theta(0,2));
     theta.SetEntry(2,0, theta(1,2));
     theta.SetEntry(2,1, theta(0,2));
     theta.SetEntry(2,2, cos(2*a));
     t.Mult(theta, help);
     if(direction == DesignElement::ROTANGLE){
       Matrix<double> dtheta(3,3);
       dtheta.SetEntry(0,0, -sin(2*a));
       dtheta.SetEntry(0,1, -dtheta(0,0));
       dtheta.SetEntry(0,2, -sqrt(2)*cos(2*a));
       dtheta.SetEntry(1,0, dtheta(0,1));
       dtheta.SetEntry(1,1, dtheta(0,0));
       dtheta.SetEntry(1,2, -dtheta(0,2));
       dtheta.SetEntry(2,0, dtheta(1,2));
       dtheta.SetEntry(2,1, dtheta(0,2));
       dtheta.SetEntry(2,2, -2*sin(2*a));
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
  assert(!(notation == HILL_MANDEL && type_ != FMO && type_ != LAMINATES && type_ != HOM_RECT && type_ != REDBAS_PARAM && type_ != REDBAS_FREE && type_ != GREEDY_PARAM && type_ != GREEDY_FREE));
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
    GetHomRectTensor(t, direction, notation);
    break;
  case GREEDY_PARAM:
  case GREEDY_FREE:
  case REDBAS_PARAM:
  case REDBAS_FREE:
    GetModRedTensor(t, direction, notation);
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
  type.Add(REDBAS_PARAM, "reducedBasis-param");
  type.Add(REDBAS_FREE, "reducedBasis-free");
  type.Add(GREEDY_PARAM, "greedy-param");
  type.Add(GREEDY_FREE, "greedy-free");

  transIsoType.SetName("DesignMaterial::TransIsoType");
  transIsoType.Add(TRANSISO_XY, "xy");
  transIsoType.Add(TRANSISO_YZ, "yz");
  transIsoType.Add(TRANSISO_XZ, "xz");

  notation.SetName("DesignMaterial::Notation");
  notation.Add(VOIGT, "voigt");
  notation.Add(HILL_MANDEL, "hill_mandel");

}

void DesignMaterial::GetSVDGTensorParameters(const Matrix<double>& G, Vector<double>& paramvec)
{
   /* G= R(theta)*(l1 0 \\ 0l2) R(phi) paramvec = (theta, phi, l1, l2) */
  //Check first that G is a 2 by 2 matrix
  assert((G.GetNumRows() == 2) & (G.GetNumCols() == 2));

  if((G.GetNumRows() != 2) || (G.GetNumCols() != 2)){
      throw Exception("GetSVDGTensor only works for 2*2 matrices");
  }


  paramvec.Resize(4);
  paramvec.Init();

  //help = G'*G
  //help2 = G*G'
  Matrix<double> help(2,2);
  Matrix<double> transpose(2,2);
  G.Transpose(transpose);
  Matrix<double> help2(2,2);
  G.MultT(G,help);
  G.Mult(transpose,help2);

  double theta = 0.5*std::atan2(help2(0,1)+help2(1,0), help2(0,0)-help2(1,1));

  double phi = -0.5*std::atan2(help(0,1)+help(1,0), help(0,0)-help(1,1));

  Matrix<double> Rthetainv(2,2);
  Matrix<double> Rphiinv(2,2);
  Rthetainv(0,0) = cos(theta);
  Rthetainv(1,1) = cos(theta);
  Rthetainv(0,1) = -sin(theta);
  Rthetainv(1,0) = sin(theta);

  Rphiinv(0,0) = cos(phi);
  Rphiinv(1,1) = cos(phi);
  Rphiinv(0,1) = -sin(phi);
  Rphiinv(1,0) = sin(phi);

  Matrix<double> t1(2,2);
  Matrix<double> t2(2,2);

  G.Mult(Rphiinv, t1);
  Rthetainv.Mult(t1,t2);

  double l1 = t2(0,0);
  double l2 = t2(1,1);

  if (l1 <0)
  {
    l1 = -l1;
    l2 = -l2;
    if (phi< 0) phi = PI +phi;
    else phi = phi - PI;
  }

  if (phi <0) phi = phi +2*PI;
  if (theta<0 ) theta = theta + 2*PI;

  paramvec[0] = theta;
  paramvec[1] = phi;
  paramvec[2] = l1;
  paramvec[3] = l2;

}

double DesignMaterial::AngleGreedyCalculus(const Vector<double>& coeffs, const double& angle)
{
   assert(coeffs.size() == 2*Na_+1);

   double s=0;
   for (UInt k= 0; k<2*Na_+1; k++ )
   {
     int indk = k - Na_;
     if (indk<0)
     {
       s += coeffs[k]*sin(-indk*angle);
     }
     else if (indk>0)
     {
       s += coeffs[k]*cos(indk*angle);
     }
     else
     {
       s += coeffs[k];
     }
   }
   return s;

}

double DesignMaterial::ScalingGreedyCalculus(const Vector<double>& coeffs, const double& l)
{
   assert(coeffs.size() == Nl_+1);
   assert((l>=lmin_) & (l<=lmax_));
   double dl = (lmax_ -lmin_)/Nl_;

   int ind = static_cast <int> (std::floor((l-lmin_)/dl));

   if (ind < 0) ind=0;
   if (ind > int(Nl_)-1) ind = Nl_-1;
   double r = l-lmin_-ind*dl;

   return (1.0-r)*coeffs[ind] + r*coeffs[ind+1];

}

