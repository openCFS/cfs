#include <iostream>
#include <fstream>

#include "baseForm.hh"

namespace CoupledField
{

  BaseForm::BaseForm(BaseFE * aptelem, BaseMaterial* matData,
		     SubTensorType type) 
    : ptelem(aptelem), isaxi_(FALSE), subTensorType_(type),
      isFracDamping_(FALSE), isRaylDamping_(FALSE), dofzero_(0)
  {
    ENTER_FCN( "BaseForm::BaseForm" );
    isSetIntPoint_ = FALSE;
    piezoMatType_ = REALMATERIALPARAMETER;

    // We generate the object, so we will delete it
    //    Error("Copy constructor not implemented",__FILE__,__LINE__);
    ptMaterial = matData;
    delMatDataAtEnd_ = true;

    baseType_ = NOTYPE;
    materialArray_ = NULL;
    softeningModel_ = "no";
  }

  BaseForm::BaseForm(BaseMaterial* matData, SubTensorType type)
    :isaxi_(FALSE), subTensorType_(type), isFracDamping_(FALSE), 
     isRaylDamping_(FALSE), dofzero_(0)
  {
    ENTER_FCN( "BaseForm::BaseForm" );
    isSetIntPoint_ = FALSE;
    ptMaterial = matData;

    // We generate the object, so we will delete it
    piezoMatType_ = REALMATERIALPARAMETER;
    delMatDataAtEnd_ = true;

    baseType_ = NOTYPE;
    materialArray_ = NULL;
    softeningModel_ = "no";   
  }

  BaseForm::BaseForm(BaseFE * aptelem)
    : ptelem(aptelem), isaxi_(FALSE), isFracDamping_(FALSE), 
      isRaylDamping_(FALSE), dofzero_(0)
  {
    ENTER_FCN( "BaseForm::BaseForm" );
    isSetIntPoint_ = FALSE;
    piezoMatType_ = REALMATERIALPARAMETER;

    // We do not generate the object, so we will not delete it
    ptMaterial = NULL;
    delMatDataAtEnd_ = false;

    baseType_ = NOTYPE;
    materialArray_ = NULL;
    softeningModel_ = "no";
  }
 
  BaseForm::BaseForm()
    : isaxi_(FALSE), isFracDamping_(FALSE), 
      isRaylDamping_(FALSE), dofzero_(0)
  {
    ENTER_FCN( "BaseForm::BaseForm" );
    isSetIntPoint_ = FALSE;
    piezoMatType_ = REALMATERIALPARAMETER;

    // We do not generate the object, so we will not delete it
    ptMaterial = NULL;
    delMatDataAtEnd_ = false;

    baseType_ = NOTYPE;
    materialArray_ = NULL;
    softeningModel_ = "no";
  }
 


  // **************
  //   Destructor
  // **************
  BaseForm::~BaseForm() {
    ENTER_FCN( "BaseForm::~BaseForm" );

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


  void BaseForm::CalcElementMatrix(Matrix<Double>& ptCoord, Matrix<Double>& elemMat) 
  {
    ENTER_FCN( "BaseForm::CalcElemMatrix" );
    Error(" Function BaseForm::CalcElementMatrix is virtual. You can use it for derived classes.",__FILE__,__LINE__);
  }


  void BaseForm::CalcComplexElementMatrix(Matrix<Double>& ptCoord, Matrix<Complex>& elemMat, Double & beta, Double & omega) 
  {
    ENTER_FCN( "BaseForm::CalcComplexElemMatrix" );
    Error(" Function BaseForm::CalcComplexElementMatrix is virtual. You can use it for derived classes.",__FILE__,__LINE__);
  }


  void BaseForm::Print(std::ostream * out, const Matrix<Double> Result) const
  {
    ENTER_FCN( "BaseForm::Print" );
    Error(" Function BaseForm::Print is virtual. You can use it for derived classes.",__FILE__,__LINE__);
  }


  // ------------- SURFACE BILINEAR FORMS -------------

  SurfForm::SurfForm() {
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
                                   const StdVector<RegionIdType> & regionIds,
                                   const StdVector<BaseMaterial*>& materials ) {
    
    firstPDEName_ = name;
    firstRegionIds_ = regionIds;
    firstMaterials_.Resize(materials.GetSize());
    for ( UInt k=0; k<materials.GetSize(); k++ ) {
      firstMaterials_[k] = materials[k];
    }
  }

  void SurfForm::SetSecondVoluInfo( const std::string & name,
                                    const StdVector<RegionIdType> & regionIds,
                                    const StdVector<BaseMaterial*>& materials ) {
    
    secondPDEName_ = name;
    secondRegionIds_ = regionIds;
    secondMaterials_.Resize(materials.GetSize());
    for ( UInt k=0; k<materials.GetSize(); k++ ) {
      secondMaterials_[k] = materials[k];
    }    
  }

  void SurfForm::SetFactor( Double factor ) {
    factor_ = factor;
  } 
  
}
