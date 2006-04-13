#include "heatNeumannInt.hh"

namespace CoupledField 
{

  
  HeatNeumannInt::HeatNeumannInt( Double factor, 
                                  Boolean isaxi ) 
    : LinearSurfForm() {
    ENTER_FCN( "HeatNeumannInt::HeatNeumannInt" );

    isaxi_ = isaxi;
    factor_ = factor;
    
  }

  
  HeatNeumannInt::~HeatNeumannInt() {
    ENTER_FCN( "HeatNeumannInt::~HeatNeumannInt" );

    actElem_ = NULL;
    
  }

  void HeatNeumannInt::CalcElemVector( Matrix<Double>& ptCoord, 
                                       Vector<Double> & elemVec ) {
    ENTER_FCN( "HeatNeumannInt::CalcElemVector" );
    
    const UInt nrIntPts = ptelem->GetNumIntPoints();
    const UInt nrNodes  = ptelem->GetNumNodes();
    const Vector<Double> & intWeights = ptelem->GetIntWeights();  
    Vector<Double> shapeFnc;
  
    elemVec.Resize(nrNodes);
    elemVec.Init();
  
    // Get thermal conductivity of relating volume neighbour region
    std::map<RegionIdType, BaseMaterial*>::iterator it = 
      materials_.find(actElem_->ptVolElem1->regionId );

    if ( it == materials_.end() && actElem_->ptVolElem2 != NULL ) {
      it = materials_.find( actElem_->ptVolElem2->regionId );
    } 
    
    if( it == materials_.end() ) {
      (*error) << "HeatNeumannInt: Surface element No." << actElem_->elemNum
               << " has no acoustic volume element neighbour!";
      Error( __FILE__, __LINE__ );
    }

    Double k;
    it->second->GetScalar(k,HEAT_CONDUCTIVITY,REAL);

    Double value = 0.0;
    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++)
      {  
        ptelem->GetShFncAtIp(shapeFnc,actIntPt);
        value = ptelem->CalcJacobianDetAtIp(actIntPt, ptCoord) * 
          intWeights[actIntPt-1] * factor_ / k;
      
        if (isaxi_)
          {
            Vector<Double> CoordAtIP;
            CoordAtIP = ptCoord * shapeFnc;
            value *=  2 * PI * CoordAtIP[0];
          }
      
        shapeFnc *= value;
        elemVec += shapeFnc;
      }

  }
    

} // end of namespace
