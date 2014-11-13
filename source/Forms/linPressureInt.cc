// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <math.h>
#include <complex>
#include <ostream>

#include "DataInOut/Logging/cfslog.hh"
#include "DataInOut/Logging/log.hpp"
#include "DataInOut/ResultCache.hh"
#include "Domain/entityList.hh"
#include "Domain/surfElem.hh"
#include "Elements/basefe.hh"
#include "General/environment.hh"
#include "MatVec/exprt/xpr1.hh"
#include "MatVec/matrix.hh"
#include "MatVec/vector.hh"
#include "Utils/mathParser/mathParser.hh"
#include "linPressureInt.hh"

namespace CoupledField {
struct Elem;
}  // namespace CoupledField

DECLARE_LOG(forms)

namespace CoupledField {

  // =================================================================
  // pressureLinForm 
  // =================================================================


  PressureLinForm::PressureLinForm( const std::string& value,
                                    const std::string& phase,
                                    const std::string& subType,
                                    bool isAxi )
    : LinearSurfForm() {

    name_ = "PressureLinForm";
    isaxi_ = isAxi;

    // store value and phase string
    value_   = value;
    phase_   = phase;
    subType_ = subType;
      
  }



  PressureLinForm::~PressureLinForm() {
  }

  Double PressureLinForm::GetPressureFactor(EntityIterator& ent){
    // register global coordinates of element midpoint
    const SurfElem* elem = ent.GetSurfElem();
    Elem* ptVolElem = elem->ptVolElem1;
    RegisterSurfElemMidPoint( mHandle_, elem, ptVolElem );
    
    // provide info to ResultCache, so we can use the input-function
    // in math parser expressions
    ResultCache::SetInfo( ResultCache::OUT_REAL, 1, ent.GetName(),
                          MECH_PRESSURE );
    ResultCache::SetIndex(elem->elemNum);
    
    // evaluate value for current element
    mParser_->SetExpr( mHandle_, value_ );
    Double factor = mParser_->Eval( mHandle_ );

    // When we do SIMP of an Piezo we might have pressure and charge density
    // on surface elements. Then scale our element contribution by the corresponding
    // volume element which has the scaling factor.
    Double density = GetErsatzMaterialFactor(ptVolElem);
    return(factor * density);
  }
    
  void PressureLinForm::CalcElemVector( Vector<Double> & elemVec,
                                        EntityIterator& ent ) {

    // compute element vector
    PrepareElemVec( elemVec, ent );

    Double factor = GetPressureFactor(ent);

    LOG_DBG3(forms) << "PressureLinForm::CalcElemVector<double> elem=" 
                    << actElem_->elemNum << " (peseudo density * factor) =" << factor;
    
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

    // provide info to ResultCache, so we can use the input-function
    // in math parser expressions
    ResultCache::SetInfo( ResultCache::OUT_AMPL, 1, ent.GetName(),
                          MECH_PRESSURE );
    ResultCache::SetIndex(actElem_->elemNum);

    // evaluate value and phase for current element
    mParser_->SetExpr( mHandle_, value_ );
    Double value = mParser_->Eval( mHandle_ );

    ResultCache::SetOutputType(ResultCache::OUT_PHASE);
    mParser_->SetExpr( mHandle_, phase_ );
    Double phase = mParser_->Eval( mHandle_ );


    // Note: Since phase is in degrees, we have to transform it into radians
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
    const UInt dim      = ptCoord_.GetNumRows();
    const Vector<Double> & intWeights = ptelem->GetIntWeights();
    const Vector<Double> * intPoints = ptelem->GetIntPoints();
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
    //       EXCEPTION("PressureLinForm: Surface element number " << actElem_->elemNum
    //                << " has no mechanic volume element neighbour!");
    //     };

    UInt numDofs;
    if( subType_ == "flatShell")
      numDofs = 5;
    else
       numDofs = dim;
    elemVec.Resize(numFncs*numDofs);
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
          ptelem->Local2GlobalCoord( CoordAtIP, intPoints[actIntPt-1],
                                     ptCoord_, ent.GetElem() );
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
              idx = j*(numDofs) +i;
              elemVec[idx] -= helpVec[j];
            }
          }
      }
  }


} // end of namespace
