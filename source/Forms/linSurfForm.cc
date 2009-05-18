// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "linSurfForm.hh"
#include "Domain/domain.hh"
#include "Domain/grid.hh"

namespace CoupledField 
{

  
  LinearSurfForm::LinearSurfForm() {

    actElem_ = NULL;   
  }

  
  LinearSurfForm::~LinearSurfForm() {
  }

  
  void LinearSurfForm::SetSurfElem( SurfElem * ptSurfElem) {

    actElem_ = ptSurfElem;
  }

  
  void LinearSurfForm::SetVoluNormal( Vector<Double> & n ) {

    normal_ = n;
  }
    
  
  void LinearSurfForm::SetVoluInfo( std::map<RegionIdType, BaseMaterial*> 
                                    & materials ) {

    materials_ = materials;
  }
    
  
  void LinearSurfForm::ExtractElemInfo( EntityIterator& it) {
    ptelem = it.GetElem()->ptElem;

 
    
    if( it.GetType() == EntityList::SURF_ELEM_LIST ) {
      actElem_ = it.GetSurfElem();
      domain->GetGrid()->CalcSurfNormal(normal_,*actElem_, coordUpdate_ );
      if( std::abs(actElem_->normalSign) > EPS ) {
        normal_ *= (Double) actElem_->normalSign;
      }
    }
    
    domain->GetGrid()->GetElemNodesCoord( ptCoord_, 
                                          it.GetElem()->connect,
                                          coordUpdate_ );
    
  }

  void LinearSurfForm::RegisterSurfElemMidPoint( MathParser::HandleType handle,
                                                 const SurfElem * ptSurfElem,
                                                 const Elem * ptVolElem ) {

    Vector<Double> locMidPoint, locMidPointVol, globMidPointVol;

    // fetch global coordinate system
    CoordSystem * coosy = domain->GetCoordSystem();
    
    // get connectivity 
    const StdVector<UInt> & surfConnect = ptSurfElem->connect;
    
    // get local surface element midpoint
    ptSurfElem->ptElem->GetCoordMidPoint(locMidPoint);
    
    //check, if volume element is zero
    if( ptVolElem != NULL ) {
      const StdVector<UInt> & volConnect = ptVolElem->connect;    

      // map local surfac element midpoint to local volume element midpoint
      ptVolElem->ptElem->GetLocalIntPoints4Surface(surfConnect, volConnect,
                                                   locMidPoint, locMidPointVol);

      // Map to global coordinate system
      Matrix<Double> volCoordMat;
      domain->GetGrid()->GetElemNodesCoord( volCoordMat, volConnect, coordUpdate_ ); 
      ptVolElem->ptElem->Local2GlobalCoord( globMidPointVol, locMidPointVol,
                                            volCoordMat, ptVolElem );
    } else {
      Matrix<Double> surfCoordMat;
      domain->GetGrid()->GetElemNodesCoord( surfCoordMat, surfConnect, coordUpdate_ );
      ptSurfElem->ptElem->Local2GlobalCoord( globMidPointVol, locMidPoint,
                                             surfCoordMat, ptSurfElem );
    }
    
    // Update variables for mathParser
    mParser_->SetCoordinates( handle, *coosy, globMidPointVol );
  }
  

} // end of namespace

