#include <iostream>
#include <fstream>

#include "baseForm.hh"
#include "Domain/domain.hh"
#include "Domain/grid.hh"

namespace CoupledField
{

 


  BaseForm::BaseForm(BaseMaterial* matData, SubTensorType type, 
                     bool coordUpdate )
    : coordUpdate_( coordUpdate),
      isaxi_(false), 
      subTensorType_(type), 
      isFracDamping_(false)
  {
    ENTER_FCN( "BaseForm::BaseForm" );
    isSetIntPoint_ = false;
    ptMaterial = matData;
    ptelem = NULL;
    sol_ = NULL;
    solDeriv1_ = NULL;
    solDeriv2_ = NULL;

    // We generate the object, so we will delete it
    matDataType_ = REAL;
    delMatDataAtEnd_ = true;

    baseType_ = NOTYPE;
    isSolDependent_ = false;
    softeningModel_ = "no";   

    // Get grip of a new math parser object.
    // This object has per default the variable
    // f (harmonic) or t (transient) registered with the
    // current value
    mHandle_ =  domain->GetMathParser()->GetNewHandle();
    mParser_ = domain->GetMathParser();
  }

 
 


  // **************
  //   Destructor
  // **************
  BaseForm::~BaseForm() {
    ENTER_FCN( "BaseForm::~BaseForm" );

    // Release math parser object
    mParser_->ReleaseHandle( mHandle_ );

//     if ( delMatDataAtEnd_ == true ) {
//       delete ptMaterial;
//     }
  }


  // ***************
  //   SetMaterial
  // ***************
  void BaseForm::SetMaterial( BaseMaterial *matPtr ) {

    ENTER_FCN( "BaseForm::SetMaterial" );

    // Avoid leaks
    if ( ptMaterial != NULL && delMatDataAtEnd_ == true ) {
      delete ptMaterial;
    }

    // Set new material pointer and avoid deleting it
    ptMaterial = matPtr;
    delMatDataAtEnd_ = false;
  }

 
  
  void BaseForm::ExtractElemInfo( EntityIterator& it ) {
    ptelem = it.GetElem()->ptElem;
    
    domain->GetGrid()->GetElemNodesCoord( ptCoord_, 
                                          it.GetElem()->connect,
                                          coordUpdate_ );
  }

  // ------------- SURFACE BILINEAR FORMS -------------

  SurfForm::SurfForm() 
    : BaseForm( NULL ){
    ENTER_FCN( "SurfForm::SurfForm" );

    factor_ = 1.0;
    formulation_ = NO_SOLUTION_TYPE;
  }


  SurfForm::~SurfForm() {
    ENTER_FCN( "SurfForm::~SurfForm" );
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

  void SurfForm::SetFactor( Double factor ) {
    factor_ = factor;
  } 
  
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
}
