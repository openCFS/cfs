// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <iostream>
#include <fstream>

#include "baseForm.hh"
#include "Domain/domain.hh"
#include "Domain/grid.hh"

namespace CoupledField
{

 
BaseForm::MaterialDescriptor::MaterialDescriptor()
{
  Init(NOT_SET, NO_CLASS, NO_MATERIAL, NO_MATERIAL, Global::INTEGER);
}


BaseForm::MaterialDescriptor::MaterialDescriptor(Type type, MaterialClass mat_class, MaterialType mat_1, Global::ComplexPart data_type)
{
  Init(type, mat_class, mat_1, NO_MATERIAL, data_type);
}


BaseForm::MaterialDescriptor::MaterialDescriptor(Type type, MaterialClass mat_class, MaterialType mat_1, MaterialType mat_2, Global::ComplexPart data_type)
{
  Init(type, mat_class, mat_1, mat_2, data_type);
}

void BaseForm::MaterialDescriptor::Init(Type type, MaterialClass mat_class, MaterialType mat_1, MaterialType mat_2, Global::ComplexPart data_type)
{
  this->type = type;
  this->mat_class = mat_class;
  this->mat_1 = mat_1;
  this->mat_2 = mat_2;
  this->data_type = data_type;
}


double BaseForm::MaterialDescriptor::GetScalar(BaseMaterial* bm)
{
  assert(data_type == Global::REAL);

  double result = -1.0;
  switch(type)
  {
  case SCALAR: // mechanic mass, acoustic stiffness
  case MINUS_SCALAR: // -1 acoustic in mech coupling case
    assert(mat_1 != NO_MATERIAL);
    bm->GetScalar(result,mat_1,data_type);

    if(type == MINUS_SCALAR)
      result *= -1.0;
    break;

  case MAT_1_MAT_1_BY_MAT_2: // acoustic mass
  case MINUS_MAT_1_MAT_1_BY_MAT_2: // -1 acoustic mass in coupling case
  {
    assert(mat_1 != NO_MATERIAL);
    bm->GetScalar(result,mat_1,data_type);
    double tmp;
    assert(mat_2 != NO_MATERIAL);
    bm->GetScalar(tmp,mat_2,data_type);
    result = (result * result) / tmp;

    if(type == MINUS_MAT_1_MAT_1_BY_MAT_2)
      result *= -1.0;
    break;
  }

  default:
    assert(false);
  }

  return result;
}


double BaseForm::MaterialDescriptor::GetErsatzMaterial(BaseForm* form, const Elem* elem, double default_mat_value)
{
  // always obtain the material data if we have a MaterialDescriptor to be able to
  // do bimaterial topology optimization. For matrices we also do it every time!
  assert((Enabled() && form->GetMaterial() != NULL) || !Enabled());
  double density = Enabled() ? GetScalar(form->GetMaterial()) : default_mat_value;

  double top_opt = form->GetErsatzMaterialFactor(elem); // returns 1 of not relevant

  // check for bimaterial optimization
  BaseMaterial* bi_mat = domain->GetErsatzBiMaterial(elem, mat_class);
  if(bi_mat != NULL && top_opt != 1.0)
  {
    double density2 = GetScalar(bi_mat);

    density = top_opt * density + (1.0 - top_opt) *  density2;
  }
  if(bi_mat == NULL && top_opt != 1.0)
  {
    density = top_opt * density;
  }

  return density;
}

  BaseForm::BaseForm( BaseMaterial* matData, SubTensorType type, 
                      bool coordUpdate )
    : coordUpdate_( coordUpdate),
      isaxi_(false), 
      subTensorType_(type)
  {
    isSetIntPoint_ = false;
    isSymmetric_ = true;
    ptMaterial = matData;
    ptelem = NULL;
    sol_ = NULL;
    sol1_ = NULL; 
    sol2_ = NULL; 
    solDeriv1_ = NULL;
    solDeriv2_ = NULL;
    TS_alg_ = NULL;

    // We generate the object, so we will delete it
    matDataType_ = Global::REAL;
    delMatDataAtEnd_ = false;

    baseType_ = NOTYPE;
    isSolDependent_ = false;
    isComplex_ = false;
    softeningModel_ = "no";   

    // Initialize form with standard Lagrange ansatz fct object
    shared_ptr<AnsatzFct> fct( new LagrangeFct() );
    ansatzFct1_ = fct;
    ansatzFct2_ = fct;

#ifndef INTEGLIB

    // Get grip of a new math parser object.
    // This object has per default the variable
    // f (harmonic) or t (transient) registered with the
    // current value
    if(domain)
    {
      mHandle_ =  domain->GetMathParser()->GetNewHandle();
      mParser_ = domain->GetMathParser();
      mParser_->SetExpr( mHandle_, "1.0" );
    }
    else
    {
      mHandle_ = 0;
    }
    
#endif
  }


  // **************
  //   Destructor
  // **************
  BaseForm::~BaseForm() {

#ifndef INTEGLIB
    // Release math parser object
    if(mHandle_)
      mParser_->ReleaseHandle( mHandle_ );
#endif

//     if ( delMatDataAtEnd_ == true ) {
//       delete ptMaterial;
//     }
  }


  // ***************
  //   SetMaterial
  // ***************
  void BaseForm::SetMaterial( BaseMaterial *matPtr ) {


    // Avoid leaks
    if ( ptMaterial != NULL && delMatDataAtEnd_ == true ) {
      delete ptMaterial;
    }

    // Set new material pointer and avoid deleting it
    ptMaterial = matPtr;
    delMatDataAtEnd_ = false;
  }

 
  void BaseForm::SetAnsatzFct( shared_ptr<AnsatzFct> actFct1,
                               shared_ptr<AnsatzFct> actFct2 ) {

    assert( actFct1 != NULL );
    ansatzFct1_ = actFct1;

    if( actFct2 != NULL ) {
      ansatzFct2_ = actFct2;
    } else {
      ansatzFct2_ = ansatzFct1_;
    } 

  }

   /** for ersatz material w and w/o SIMP. */
  Double BaseForm::GetErsatzMaterialFactor(const Elem* elem)
  {
    Double factor;
    bool ok = domain->GetErsatzMaterial(elem, this, factor);
    return ok ? factor : 1.0;
  }
  
#ifndef INTEGLIB
  void BaseForm::ExtractElemInfo( EntityIterator& it ) {
    ptelem = it.GetElem()->ptElem;
    it1_ = it;

    domain->GetGrid()->GetElemNodesCoord( ptCoord_, 
                                          it.GetElem()->connect,
                                          coordUpdate_ );
  }
#endif




  // ------------- SURFACE BILINEAR FORMS -------------

  SurfForm::SurfForm() 
    : BaseForm( NULL ){

  }


  SurfForm::~SurfForm() {
  }
  
  void SurfForm::SetSurfElem( SurfElem * ptSurfElem) {
    actElem_ = ptSurfElem;
  }
  
  void SurfForm::SetFirstVoluNormal( Vector<Double> & n ) {
    normal_ = n;
  }

  void SurfForm::SetFirstVoluInfo( const std::string & name,
                                   const std::map<RegionIdType, BaseMaterial*>& materials ) {
    
    firstPDEName_ = name;
    firstMaterials_ = materials;
  }

  void SurfForm::SetSecondVoluInfo( const std::string & name,
                                    const std::map<RegionIdType, BaseMaterial*>& materials ) {                         
    
    secondPDEName_ = name;
    secondMaterials_ = materials;
  }

  
#ifndef INTEGLIB
  void SurfForm::ExtractElemInfo( EntityIterator& it ) {
    ptelem = it.GetElem()->ptElem;
    
    if( it.GetType() == EntityList::SURF_ELEM_LIST ) {
      actElem_ = it.GetSurfElem();
      domain->GetGrid()->CalcSurfNormal(normal_,*actElem_ );
      normal_ *= (Double) actElem_->normalSign;
    }

    domain->GetGrid()->GetElemNodesCoord( ptCoord_, 
                                          it.GetElem()->connect,
                                          coordUpdate_ );
  }
  
  
#endif

}
