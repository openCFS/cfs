#include "acouNeumannInt.hh"

namespace CoupledField 
{

  
  AcouNeumannInt::AcouNeumannInt( Double factor, 
                                  Boolean isaxi ) 
    : LinearSurfForm() {
    ENTER_FCN( "AcouNeumannInt::AcouNeumannInt" );

    isaxi_ = isaxi;
    factor_ = factor;
    
  }

  
  AcouNeumannInt::~AcouNeumannInt() {
    ENTER_FCN( "AcouNeumannInt::~AcouNeumannInt" );

    actElem_ = NULL;
    materials_ = NULL;
    
  }

  void AcouNeumannInt::CalcElemVector( Matrix<Double>& ptCoord, 
                                       Vector<Double> & elemVec ) {
    ENTER_FCN( "AcouNeumannInt::CalcElemVector" );
    
    const UInt nrIntPts = ptelem->GetNumIntPoints();
    const UInt nrNodes  = ptelem->GetNumNodes();
    const Vector<Double> & intWeights = ptelem->GetIntWeights();  
    Vector<Double> shapeFnc;
  
    elemVec.Resize(nrNodes);
    elemVec.Init();
  
    // Get density of relating volume neighbour region
    Integer index = regionIds_.Find( actElem_->ptVolElem1->regionId );
    if ( index == -1 && actElem_->ptVolElem2 != NULL ) {
      index = regionIds_.Find( actElem_->ptVolElem2->regionId );
    } 
    
    if( index == -1 ) {
      (*error) << "AcouNeumannInt: Surface element number " << actElem_->elemNum
               << " has no acoustic volume element neighbour!";
      Error( __FILE__, __LINE__ );
    }

    Double density = materials_[index].GetDensity();

    Double value = 0.0;
    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++)
      {  
        ptelem->GetShFncAtIp(shapeFnc,actIntPt);
        value = ptelem->CalcJacobianDetAtIp(actIntPt, ptCoord) * 
          intWeights[actIntPt-1] * density * factor_;
      
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
