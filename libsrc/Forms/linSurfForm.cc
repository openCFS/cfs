#include "linSurfForm.hh"
#include "Domain/domain.hh"
#include "Domain/grid.hh"

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
    
  
  void LinearSurfForm::SetVoluInfo( std::map<RegionIdType, BaseMaterial*> 
                                    & materials ) {
    ENTER_FCN( "LinearSurfForm::SetVoluInfo" );

    materials_ = materials;
  }
    
  
  void LinearSurfForm::SetFactor(Double factor) {
    ENTER_FCN( "LinearSurfForm::SetFactor" );

    factor_ = factor;
  }

  void LinearSurfForm::ExtractElemInfo( EntityIterator& it) {
    ENTER_FCN( "LinearSurfForm::ExtractElemInfo" );
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
  

} // end of namespace

