#include "linSurfForm.hh"

namespace CoupledField 
{

  
  LinearSurfForm::LinearSurfForm() {
    ENTER_FCN( "LinearSurfForm::LinearSurfForm" );

    factor_ = 0.0;
    actElem_ = NULL;   
  }

  
  LinearSurfForm::~LinearSurfForm() {
    ENTER_FCN( "LinearSurfForm::~LinearSurfForm" );
  }

  
  void LinearSurfForm::SetSurfElem( SurfElem * ptSurfElem) {
    ENTER_FCN( "LinearSurfForm::SetSurfElem" );

    actElem_ = ptSurfElem;
  }

  
  void LinearSurfForm::SetVoluNormal( Vector<Double> & n ) {
    ENTER_FCN( "LinearSurfForm::SetVoluNormal" );

    normal_ = n;
  }
    
  
  void LinearSurfForm::SetVoluInfo( const StdVector<RegionIdType> & regionIds,
                                    const StdVector<BaseMaterial*>& materials ) {
    ENTER_FCN( "LinearSurfForm::SetVoluInfo" );

    regionIds_ = regionIds;
    materials_.Resize(materials.GetSize());
    for ( UInt k=0; k<materials.GetSize(); k++ ) {
      materials_[k] = materials[k];
    }    
  }
    
  
  void LinearSurfForm::SetFactor(Double factor) {
    ENTER_FCN( "LinearSurfForm::SetFactor" );

    factor_ = factor;
  }
   

} // end of namespace

