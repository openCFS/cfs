#include "linSurfForm.hh"
#include "Domain/domain.hh"
#include "Domain/grid.hh"

namespace CoupledField 
{

  
  LinearSurfForm::LinearSurfForm() {
    ENTER_FCN( "LinearSurfForm::LinearSurfForm" );

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

  void LinearSurfForm:: 
  RegisterSurfElemMidPoint( MathParser::HandleType handle,
                            const SurfElem * ptSurfElem,
                            const Elem * ptVolElem ) {
    ENTER_FCN( "LinearSurfForm::RegisterSurfElemMidPoint" );

    Vector<Double> locMidPoint, locMidPointVol, globMidPointVol;

    // fetch global coordinate system
    CoordSystem * coosy = domain->GetCoordSystem();
    
    // get connectivity 
    const StdVector<UInt> & surfConnect = ptSurfElem->connect;
    const StdVector<UInt> & volConnect = ptVolElem->connect;    

    // get local surface element midpoint
    ptSurfElem->ptElem->GetCoordMidPoint(locMidPoint);

    // map local surfac element midpoint to local volume element midpoint
    ptVolElem->ptElem->GetLocalIntPoints4Surface(surfConnect, volConnect,
                                                 locMidPoint, locMidPointVol);
    
    // Map to global coordinate system
    Matrix<Double> volCoordMat;
    domain->GetGrid()->GetElemNodesCoord( volCoordMat, volConnect, coordUpdate_ ); 
    ptVolElem->ptElem->Local2GlobalCoord( globMidPointVol, locMidPointVol,
                                           volCoordMat, ptVolElem );
    
    // Update variables for mathParser
    mParser_->SetCoordinates( handle, *coosy, globMidPointVol );
  }
  

} // end of namespace

