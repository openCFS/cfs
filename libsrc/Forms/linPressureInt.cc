#include "linPressureInt.hh"

namespace CoupledField {

  // =================================================================
  // pressureLinForm 
  // =================================================================


  PressureLinForm::PressureLinForm(Double aVal, bool isaxi)
    : LinearSurfForm(), multiplier_(aVal)
  {
    ENTER_FCN( "PressureLinForm::PressureLinForm" );

    name_ = "PressureLinForm";
    isaxi_ = isaxi;
  }



  PressureLinForm::~PressureLinForm()
  {
    ENTER_FCN( "PressureLinForm::~PressureLinForm" );
  }


  void PressureLinForm::CalcElemVector(Vector<Double> & elemVec,  
                                       EntityIterator& ent ) 
  {
    ENTER_FCN( "PressureLinForm::CalcElemVector" );
    
    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent );

    ptelem->SetAnsatzFct( ansatzFct1_ );
    UInt numFncs = ptelem->GetNumFncs( ansatzFct1_ );
    const UInt nrIntPts = ptelem->GetNumIntPoints();
    const UInt dim      = ptCoord_.GetSizeRow();
    const Vector<Double> & intWeights = ptelem->GetIntWeights();
    Vector<Double> shapeFnc;
    

 //    // calculcate correct orientation of normal
//     std::map<RegionIdType, BaseMaterial*>::iterator it = 
//       materials_.find( actElem_->ptVolElem1->regionId );

//     if ( it == materials_.end() && actElem_->ptVolElem2 != NULL ) {
//       it = materials_.find( actElem_->ptVolElem2->regionId );
//     } else {
//       normal_ *= -1;
//     }
    
//     if( it == materials_.end() ) {
//       (*error) << "PressureLinForm: Surface element number " << actElem_->elemNum
//                << " has no mechanic volume element neighbour!";
//       Error( __FILE__, __LINE__ );
//     };

    elemVec.Resize(numFncs*dim);
    elemVec.Init(0);

    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++)
      {
        ptelem->GetShFncAtIp(shapeFnc, actIntPt, ent.GetElem() );
        Double jacDet = ptelem->CalcJacobianDetAtIp(actIntPt, ptCoord_,
                                                    ent.GetElem() );
        Double factor = multiplier_ * intWeights[actIntPt-1] * jacDet;

        if (isaxi_)
          {
            Vector<Double> CoordAtIP;
            CoordAtIP = ptCoord_ * shapeFnc;
            factor *=  2 * PI * CoordAtIP[0];
          }
      
        shapeFnc *= factor;      
     
        Vector<Double> helpVec;
        UInt idx;
        for (UInt i=0; i<dim; i++) 
          {
            //multiply with corresponding component of normal vector
            //to get the x-,y-,z-component
            helpVec = shapeFnc * normal_[i];
            for (UInt j=0; j<helpVec.GetSize(); j++) {
              idx = j*(dim) +i;
              elemVec[idx] -= helpVec[j];
            }
          }
      }
  } // end of method



} // end of namespace
