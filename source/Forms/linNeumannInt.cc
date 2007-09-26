// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "linNeumannInt.hh"

#include "Domain/domain.hh"
#include "Domain/grid.hh"

namespace CoupledField 
{

  
  LinNeumannInt::LinNeumannInt( std::string amplitudeStr,
                                std::string phaseStr,
                                MaterialType materialParam,
                                bool isaxi ) 
    : LinearSurfForm() {

    name_ = "LinNeumannInt";
    isaxi_ = isaxi;

    // store value and phase string
    amplitude_ = amplitudeStr;
    phase_ = phaseStr;
    materialParam_ = materialParam;

#ifdef DEBUG
    (*debug) << "In LinNeumannInt::CalcElemVector output " << std::endl
             << " coordinate     amplitude    phase " << std::endl;
#endif

  }


  LinNeumannInt::~LinNeumannInt() {

    actElem_ = NULL;

  }

  void LinNeumannInt::PrepareElemVector( Vector<Double> & elemVec,
                                      EntityIterator& ent) {

    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent );

    ptelem->SetAnsatzFct( ansatzFct1_ );
    const UInt nrIntPts = ptelem->GetNumIntPoints();
    UInt numFncs = ptelem->GetNumFncs( ansatzFct1_ );

    const Vector<Double> & intWeights = ptelem->GetIntWeights();  
    Vector<Double> shapeFnc, CoordAtIP;


    // Find material of region of relating volume element neighbour
    std::map<RegionIdType, BaseMaterial*>::iterator it = 
      materials_.find(actElem_->ptVolElem1->regionId );
    ptVolElem_ = actElem_->ptVolElem1;

    if ( it == materials_.end() && actElem_->ptVolElem2 != NULL ) {
      it = materials_.find( actElem_->ptVolElem2->regionId );
      ptVolElem_ = actElem_->ptVolElem2;
    } 
    
    if( it == materials_.end() ) {
      (*error) << "LinNeumannInt: Surface element number " << actElem_->elemNum
               << " has no corresponding volume element neighbour!";
      Error( __FILE__, __LINE__ );
    }

    Double factor;
    switch(materialParam_)
    {
       case DENSITY: 
           it->second->GetScalar(factor,DENSITY,REAL);
           break;
           
       case HEAT_CONDUCTIVITY:
            {
              Double k;
              it->second->GetScalar(k,HEAT_CONDUCTIVITY,REAL);
              factor = 1.0 / k;
            }
            break;
            
       case NO_MATERIAL:
            factor = 1.0;
            break;
            
       default: EXCEPTION("material parameter " << materialParam_ << " not implemented");     
    }            

    // Calculate element vector  
    elemVec.Resize(numFncs);
    elemVec.Init(0.0);
    Double value = 0.0;
    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++) {
      
      ptelem->GetShFncAtIp(shapeFnc,actIntPt,ent.GetElem() );
      value = ptelem->CalcJacobianDetAtIp(actIntPt, ptCoord_, ent.GetElem()) * 
        intWeights[actIntPt-1] * factor;
      
      if (isaxi_) {
        
        CoordAtIP = ptCoord_ * shapeFnc;
        value *=  2 * PI * CoordAtIP[0];
      }
      
      shapeFnc *= value;
      elemVec += shapeFnc;
    }
  }



  void LinNeumannInt::CalcElemVector( Vector<Double> & elemVec,
                                      EntityIterator& ent) {
    
    // compute element vector
    PrepareElemVector( elemVec, ent );

    // register global coordinates of element midpoint
    Elem * ptVolElem = actElem_->ptVolElem1;
    RegisterSurfElemMidPoint( mHandle_, ent.GetSurfElem(), ptVolElem );

    // evaluate value for current element
    mParser_->SetExpr( mHandle_, amplitude_ );
    Double factor = mParser_->Eval( mHandle_ );

    // When we do SIMP of an Piezo we might have pressure and charge density
    // on surface elements. Then scale our element contribution by the corresponding
    // volume element which has the scaling factor.
    double density = GetErsatzMaterialFactor(ptVolElem);
    factor *= density;

    // we assume the surface normal points out of domain,
    //  but we want to take deflections into the domain positive
    factor *= -1.0;

    // multiply element vector with factor
    elemVec *= factor;
  }

  void LinNeumannInt::CalcElemVector( Vector<Complex> & elemVec,
                                      EntityIterator& ent) {

    // compute element vector
    Vector<Double> helpVec;
    PrepareElemVector( helpVec, ent );
    
    // register global coordinates of element midpoint
    Elem * ptVolElem = actElem_->ptVolElem1;
    RegisterSurfElemMidPoint( mHandle_, ent.GetSurfElem(), ptVolElem );


    // evaluate value and phase for current element
    mParser_->SetExpr( mHandle_, amplitude_ );
    Double value = mParser_->Eval( mHandle_ );
    // we assume the surface normal points out of domain,
    //  but we want to take deflections into the domain positive
    value *= - 1.0;

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


} // end of namespace
