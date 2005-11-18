#include "linPressureInt.hh"

namespace CoupledField {

  // =================================================================
  // pressureLinForm 
  // =================================================================


  PressureLinForm::PressureLinForm(Double aVal, Boolean isaxi)
    : LinearSurfForm(), multiplier_(aVal)
  {
    ENTER_FCN( "PressureLinForm::PressureLinForm" );
    isaxi_ = isaxi;
  }



  PressureLinForm::~PressureLinForm()
  {
    ENTER_FCN( "PressureLinForm::~PressureLinForm" );
  }


  void PressureLinForm::CalcElemVector(Matrix<Double>& ptCoord, 
                                       Vector<Double> & elemVec)
  {
    ENTER_FCN( "PressureLinForm::CalcElemVector" );

    const UInt nrIntPts = ptelem->GetNumIntPoints();
    const UInt nrNodes  = ptelem->GetNumNodes();
    const UInt dim      = ptCoord. GetSizeRow();
    const Vector<Double> & intWeights = ptelem->GetIntWeights();
    Vector<Double> shapeFnc;
    

    // calculcate correct orientation of normal
    Integer index = regionIds_.Find( actElem_->ptVolElem1->regionId );
    if ( index == -1 && actElem_->ptVolElem2 != NULL ) {
      index = regionIds_.Find( actElem_->ptVolElem2->regionId );
    } else {
      normal_ *= -1;
    }
    
    if( index == -1 ) {
      (*error) << "PressureLinForm: Surface element number " << actElem_->elemNum
               << " has no mechanic volume element neighbour!";
      Error( __FILE__, __LINE__ );
    };

    UInt zerodim = 0;
    if (dofzero_ >0) {
      //now ndDofs = dim + 1 for electric potential
      elemVec.Resize(nrNodes*(dim+1));
      zerodim = 1;
    }
    else
      elemVec.Resize(nrNodes*dim);

    elemVec.Init(0);    // set elems to 0

    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++)
      {
        ptelem->GetShFncAtIp(shapeFnc, actIntPt);
        Double jacDet = ptelem->CalcJacobianDetAtIp(actIntPt, ptCoord);
        Double factor = multiplier_ * intWeights[actIntPt-1] * jacDet;

        if (isaxi_)
          {
            Vector<Double> CoordAtIP;
            CoordAtIP = ptCoord * shapeFnc;
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
              idx = j*(dim+zerodim) +i;
              elemVec[idx] -= helpVec[j];
            }
          }
      }
  } // end of method



} // end of namespace
