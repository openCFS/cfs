// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "linPressureInt.hh"



#include "DataInOut/Logging/cfslog.hh"

DECLARE_LOG(forms)

namespace CoupledField {

  // =================================================================
  // pressureLinForm 
  // =================================================================


  PressureLinForm::PressureLinForm( const std::string& value,
                                    const std::string& phase,
                                    bool isAxi )
    : LinearSurfForm() {

    name_ = "PressureLinForm";
    isaxi_ = isAxi;

    // store value and phase string
    value_ = value;
    phase_ = phase;
      
  }



  PressureLinForm::~PressureLinForm() {
  }



  void PressureLinForm::CalcElemVector( Vector<Double> & elemVec,
                                        EntityIterator& ent ) {

    // compute element vector
    PrepareElemVec( elemVec, ent );
    
    // register global coordinates of element midpoint
    Elem * ptVolElem = actElem_->ptVolElem1;
    RegisterSurfElemMidPoint( mHandle_, ent.GetSurfElem(), ptVolElem );
    
    // evaluate value for current element
    mParser_->SetExpr( mHandle_, value_ );
    Double factor = mParser_->Eval( mHandle_ );

    // When we do SIMP of an Piezo we might have pressure and charge density
    // on surface elements. Then scale our element contribution by the corresponding
    // volume element which has the scaling factor.
    double density = GetErsatzMaterialFactor(ptVolElem);
    factor *= density;
    
    LOG_DBG3(forms) << "PressureLinForm::CalcElemVector<double> elem=" 
                    << actElem_->elemNum << " peseudo density=" << density
                    << " factor=" << factor; 
    
    // multiply element vector with factor
    elemVec *= factor;
  }

  void PressureLinForm::CalcElemVector( Vector<Complex> & elemVec,
                                        EntityIterator& ent ) {

    // compute element vector
    Vector<Double> helpVec;
    PrepareElemVec( helpVec, ent );
    
    // register global coordinates of element midpoint
    Elem * ptVolElem = actElem_->ptVolElem1;
    RegisterSurfElemMidPoint( mHandle_, ent.GetSurfElem(), ptVolElem );


    // evaluate value and phase for current element
    mParser_->SetExpr( mHandle_, value_ );
    Double value = mParser_->Eval( mHandle_ );

    mParser_->SetExpr( mHandle_, phase_ );
    Double phase = mParser_->Eval( mHandle_ );


    // Note: Since phase is in (grad), we have to transform it into
    //        rad-value
    Double realPart = value * cos( phase / 180 * PI );
    Double imagPart = value * sin( phase / 180 * PI );

    // multiply element vector with complex factor
    Complex factor(realPart, imagPart);
    elemVec = helpVec * factor;
  }

  void PressureLinForm::PrepareElemVec( Vector<Double>& elemVec, 
                                        EntityIterator& ent ) {
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
        Double factor = intWeights[actIntPt-1] * jacDet;

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
  }


} // end of namespace
