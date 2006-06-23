#include "linNeumannInt.hh"

namespace CoupledField 
{

  
  LinNeumannInt::LinNeumannInt( Double amplitude,
                                MaterialType materialParam,
                                bool isaxi ) 
    : LinearSurfForm() {
    ENTER_FCN( "LinNeumannInt::LinNeumannInt" );

    isaxi_ = isaxi;
    amplitude_ = amplitude;
    materialParam_ = materialParam;
    
  }

  
  LinNeumannInt::~LinNeumannInt() {
    ENTER_FCN( "LinNeumannInt::~LinNeumannInt" );

    actElem_ = NULL;
    
  }

  void LinNeumannInt::CalcElemVector( Vector<Double> & elemVec,
                                      EntityIterator& ent) {
    ENTER_FCN( "LinNeumannInt::CalcElemVector" );
    
    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent );
    
    const UInt nrIntPts = ptelem->GetNumIntPoints();
    const UInt nrNodes  = ptelem->GetNumNodes();
    const Vector<Double> & intWeights = ptelem->GetIntWeights();  
    Vector<Double> shapeFnc;
  
    elemVec.Resize(nrNodes);
    elemVec.Init();
  
    // Get density of relating volume neighbour region
    std::map<RegionIdType, BaseMaterial*>::iterator it = 
      materials_.find(actElem_->ptVolElem1->regionId );

    if ( it == materials_.end() && actElem_->ptVolElem2 != NULL ) {
      it = materials_.find( actElem_->ptVolElem2->regionId );
    } 
    
    if( it == materials_.end() ) {
      (*error) << "LinNeumannInt: Surface element number " << actElem_->elemNum
               << " has no acoustic volume element neighbour!";
      Error( __FILE__, __LINE__ );
    }

    Double factor;
    if (materialParam_ == DENSITY) {
      it->second->GetScalar(factor,DENSITY,REAL);
    }
    else if (materialParam_ == HEAT_CONDUCTIVITY) {
      Double k;
      it->second->GetScalar(k,HEAT_CONDUCTIVITY,REAL);
      factor = 1.0 / k;
    }


    Double value = 0.0;
    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++) {
      
      ptelem->GetShFncAtIp(shapeFnc,actIntPt);
      value = ptelem->CalcJacobianDetAtIp(actIntPt, ptCoord_) * 
        intWeights[actIntPt-1] * amplitude_ * factor;
      
      if (isaxi_) {
        
        Vector<Double> CoordAtIP;
        CoordAtIP = ptCoord_ * shapeFnc;
        value *=  2 * PI * CoordAtIP[0];
      }
      
      shapeFnc *= value;
      elemVec += shapeFnc;
    }
  }
} // end of namespace
