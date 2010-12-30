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
  this->type = NOT_SET;
}

BaseForm::MaterialDescriptor::MaterialDescriptor(Type type, MaterialType matType, Global::ComplexPart dataType)
{
  this->type = type;
  this->matType = matType;
  this->dataType = dataType;
}

double BaseForm::MaterialDescriptor::GetScalar(BaseMaterial* bm)
{
  assert(type == SCALAR);
  assert(dataType == Global::REAL);
  // set the density
  Double val;
  bm->GetScalar(val,matType,dataType);
  return val;
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
