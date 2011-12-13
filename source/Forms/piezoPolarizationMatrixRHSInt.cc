// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <assert.h>
#include <stddef.h>
#include <map>
#include <ostream>
#include <string>

#include "DataInOut/Logging/cfslog.hh"
#include "DataInOut/Logging/log.hpp"
#include "Domain/domain.hh"
#include "Domain/elem.hh"
#include "Domain/entityList.hh"
#include "Domain/grid.hh"
#include "Domain/surfElem.hh"
#include "Elements/basefe.hh"
#include "General/environment.hh"
#include "General/exception.hh"
#include "MatVec/exprt/xpr1.hh"
#include "piezoPolarizationMatrixRHSInt.hh"

namespace CoupledField {
class BaseMaterial;
}  // namespace CoupledField

DECLARE_LOG(forms)

namespace CoupledField
{

  ////////////////////////////
  // MECH
  ////////////////////////////
  PiezoPolarizationMatrixMechRHSInt::PiezoPolarizationMatrixMechRHSInt(const Vector<Double> &vals,
      const int num) :
    LinearSurfForm(), number(num), rhsvals(vals)
  {
    // we need exactly 3 (2D) / 6 (3D) values
    assert(vals.GetSize() == (domain->GetGrid()->GetDim() == 2 ? 3 : 6));
    name_ = "PiezoPolarizationMatrixMechRHSInt";
    
    // ptMaterial = mat;
  }


  PiezoPolarizationMatrixMechRHSInt::~PiezoPolarizationMatrixMechRHSInt()
  {
    actElem_ = NULL;
  }

  void PiezoPolarizationMatrixMechRHSInt::PrepareElemVector(Vector<Double> & elemVec, EntityIterator& ent)
  {
    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent );
    
    ptelem->SetAnsatzFct( ansatzFct1_ );
    const UInt nrIntPts = ptelem->GetNumIntPoints();
    const UInt numFncs = ptelem->GetNumFncs( ansatzFct1_ );

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
      EXCEPTION( "PiezoPolarizationMatrixRHSInt: Surface element number " << actElem_->elemNum
               << " has no corresponding volume element neighbour!" );
    }

    /*
    Double factor(1.0);
    switch(materialParam_){
       case DENSITY:
           it->second->GetScalar(factor,DENSITY,Global::REAL); break;
       case HEAT_CONDUCTIVITY:
            { Double k;
              it->second->GetScalar(k,HEAT_CONDUCTIVITY,Global::REAL);
              factor = 1.0 / k; } break;
       case NO_MATERIAL:
            factor = 1.0; break;
       default: EXCEPTION("material parameter " << materialParam_ << " not implemented"); }
    */
    
    // Calculate element vector
    elemVec.Resize(numFncs);
    elemVec.Init();
    Double value(0.0);
    for(UInt actIntPt = 1; actIntPt <= nrIntPts; ++actIntPt)
    {
      ptelem->GetShFncAtIp(shapeFnc, actIntPt, ent.GetElem());
      value = ptelem->CalcJacobianDetAtIp(actIntPt, ptCoord_, ent.GetElem()) * intWeights[actIntPt-1];
      // * factor;

      assert(!isaxi_); // axisymmetric case cannot be handled at the moment!
      
      shapeFnc *= value;
      elemVec += shapeFnc;
    }
  }

  void PiezoPolarizationMatrixMechRHSInt::CalcElemVector( Vector<Double> & elemVec, EntityIterator& ent)
  {
    // compute element vector
    PrepareElemVector( elemVec, ent );

    // register global coordinates of element midpoint
    Elem * ptVolElem = actElem_->ptVolElem1;
    RegisterSurfElemMidPoint( mHandle_, ent.GetSurfElem(), ptVolElem );

    /*
    
    // evaluate value for current element
    mParser_->SetExpr( mHandle_, amplitude_ );
    Double factor = mParser_->Eval( mHandle_ );

    // When we do SIMP of an Piezo we might have pressure and charge density
    // on surface elements. Then scale our element contribution by the corresponding
    // volume element which has the scaling factor.
    double density = GetErsatzMaterialFactor(ptVolElem);
    factor *= density;
    */
    
    // we assume the surface normal points out of domain,
    //  but we want to take deflections into the domain positive
    Double factor(-1.0);
    
    LOG_DBG3(forms) << "LinNeumannInt::CalcElemVector<double> elem = "
                    //<< actElem_->elemNum << " pseudo density=" << density
                    << ", factor = " << factor 
                    << ", elem vec = " << elemVec.ToString();

    // multiply element vector with factor
    elemVec *= factor;
  }
  
  
  ////////////////////////////
  // ELEC
  ////////////////////////////

  PiezoPolarizationMatrixElecRHSInt::PiezoPolarizationMatrixElecRHSInt(const Vector<Double> &vals,
      const int num) :
    LinearSurfForm(), number(num), rhsvals(vals)
  {
    // we need exactly 2 (2D) / 3 (3D) values
    assert(vals.GetSize() == domain->GetGrid()->GetDim());
    name_ = "PiezoPolarizationMatrixMechRHSInt";
    // ptMaterial = mat;
  }


  PiezoPolarizationMatrixElecRHSInt::~PiezoPolarizationMatrixElecRHSInt()
  {
    actElem_ = NULL;
  }

  void PiezoPolarizationMatrixElecRHSInt::PrepareElemVector(Vector<Double> & elemVec, EntityIterator& ent)
  {
    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent );
    
    ptelem->SetAnsatzFct( ansatzFct1_ );
    const UInt nrIntPts = ptelem->GetNumIntPoints();
    const UInt numFncs = ptelem->GetNumFncs( ansatzFct1_ );

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
      EXCEPTION( "PiezoPolarizationMatrixRHSInt: Surface element number " << actElem_->elemNum
               << " has no corresponding volume element neighbour!" );
    }

    /*
    Double factor(1.0);
    switch(materialParam_){
       case DENSITY:
           it->second->GetScalar(factor,DENSITY,Global::REAL); break;
       case HEAT_CONDUCTIVITY:
            { Double k;
              it->second->GetScalar(k,HEAT_CONDUCTIVITY,Global::REAL);
              factor = 1.0 / k; } break;
       case NO_MATERIAL:
            factor = 1.0; break;
       default: EXCEPTION("material parameter " << materialParam_ << " not implemented"); }
    */
    
    // Calculate element vector
    elemVec.Resize(numFncs);
    elemVec.Init();
    Double value(0.0);
    for(UInt actIntPt = 1; actIntPt <= nrIntPts; ++actIntPt)
    {
      ptelem->GetShFncAtIp(shapeFnc, actIntPt, ent.GetElem());
      value = ptelem->CalcJacobianDetAtIp(actIntPt, ptCoord_, ent.GetElem()) * intWeights[actIntPt-1];
      // * factor;

      assert(!isaxi_); // axisymmetric case cannot be handled at the moment!
      
      shapeFnc *= value;
      elemVec += shapeFnc;
    }
  }

  void PiezoPolarizationMatrixElecRHSInt::CalcElemVector( Vector<Double> & elemVec, EntityIterator& ent)
  {
    // compute element vector
    PrepareElemVector( elemVec, ent );

    // register global coordinates of element midpoint
    Elem * ptVolElem = actElem_->ptVolElem1;
    RegisterSurfElemMidPoint( mHandle_, ent.GetSurfElem(), ptVolElem );

    /*
    
    // evaluate value for current element
    mParser_->SetExpr( mHandle_, amplitude_ );
    Double factor = mParser_->Eval( mHandle_ );

    // When we do SIMP of an Piezo we might have pressure and charge density
    // on surface elements. Then scale our element contribution by the corresponding
    // volume element which has the scaling factor.
    double density = GetErsatzMaterialFactor(ptVolElem);
    factor *= density;
    */
    
    // we assume the surface normal points out of domain,
    //  but we want to take deflections into the domain positive
    Double factor(-1.0);
    
    LOG_DBG3(forms) << "LinNeumannInt::CalcElemVector<double> elem = "
                    //<< actElem_->elemNum << " pseudo density=" << density
                    << ", factor = " << factor 
                    << ", elem vec = " << elemVec.ToString();

    // multiply element vector with factor
    elemVec *= factor;
  }

} // end of namespace
