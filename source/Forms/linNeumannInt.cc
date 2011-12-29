// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <math.h>
#include <stddef.h>
#include <complex>
#include <map>
#include <ostream>
#include <utility>

#include "DataInOut/Logging/cfslog.hh"
#include "DataInOut/Logging/log.hpp"
#include "DataInOut/ResultCache.hh"
#include "Domain/domain.hh"
#include "Domain/elem.hh"
#include "Domain/entityList.hh"
#include "Domain/grid.hh"
#include "Domain/surfElem.hh"
#include "Elements/basefe.hh"
#include "General/exception.hh"
#include "MatVec/exprt/xpr1.hh"
#include "MatVec/exprt/xpr2.hh"
#include "MatVec/matrix.hh"
#include "MatVec/vector.hh"
#include "Materials/baseMaterial.hh"
#include "Utils/StdVector.hh"
#include "Utils/mathParser/mathParser.hh"
#include "linNeumannInt.hh"

namespace CoupledField {
class CoordSystem;
}  // namespace CoupledField

DECLARE_LOG(forms)

namespace CoupledField
{


  LinNeumannInt::LinNeumannInt( std::string amplitudeStr,
                                std::string phaseStr,
                                SolutionType quantity,
                                MaterialType materialParam,
                                bool isaxi )
    : LinearSurfForm() {

    name_ = "LinNeumannInt";
    isaxi_ = isaxi;

    // store value and phase string
    amplitude_ = amplitudeStr;
    phase_ = phaseStr;
    quantity_ = quantity;
    materialParam_ = materialParam;
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
      EXCEPTION( "LinNeumannInt: Surface element number " << actElem_->elemNum
               << " has no corresponding volume element neighbour!" );
    }

    Double factor;
    switch(materialParam_)
    {
       case DENSITY:
           it->second->GetScalar(factor,DENSITY,Global::REAL);
           break;

       case HEAT_CONDUCTIVITY:
            {
              Double k;
              it->second->GetScalar(k,HEAT_CONDUCTIVITY,Global::REAL);
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
    elemVec.Init();
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
    Double density, factor;
    
    // compute element vector
    PrepareElemVector( elemVec, ent );

    Elem * ptVolElem = actElem_->ptVolElem1;
    //RegisterSurfElemMidPoint( mHandle_, ent.GetSurfElem(), ptVolElem );

    // register global coordinates of nodes
    const StdVector<UInt> & surfConnect = ent.GetSurfElem()->connect;
    Matrix<Double> surfCoordMat;
    domain->GetGrid()->GetElemNodesCoord( surfCoordMat, surfConnect, coordUpdate_ );
    Vector<Double> nodeCoord;
    nodeCoord.Resize(surfCoordMat.GetNumRows());

    // fetch global coordinate system
    CoordSystem * coosy = domain->GetCoordSystem();

    // set info for ResultCache, which handles the input function
    // (reads data from input files, i.e. hdf5, unv, etc.)
    // in math parser expressions
    ResultCache::SetInfo( ResultCache::OUT_REAL,
                          0,
                          ent.GetName(),
                          quantity_ );

    mParser_->SetExpr( mHandle_, amplitude_ );

    for ( UInt i=0, n=elemVec.GetSize(); i<n; ++i )
    {
      for ( UInt j=0, m=nodeCoord.GetSize(); j<m; ++j )
      {
        nodeCoord[j] = surfCoordMat[j][i];
      }
      mParser_->SetCoordinates( mHandle_, *coosy, nodeCoord );

      ResultCache::SetIndex(surfConnect[i]);
      
      // evaluate value for current element
      factor = mParser_->Eval( mHandle_ );

      // When we do SIMP of a Piezo we might have pressure and charge density
      // on surface elements. Then scale our element contribution by the
      // corresponding volume element which has the scaling factor.
      density = GetErsatzMaterialFactor(ptVolElem);
      factor *= density;

      // we assume the surface normal points out of domain,
      //  but we want to take deflections into the domain positive
      factor *= -1.0;

      // multiply element vector with factor
      elemVec[i] = elemVec[i] * factor;
    }
    
    LOG_DBG3(forms) << "LinNeumannInt::CalcElemVector<double> elem="
                    << actElem_->elemNum << " pseudo density=" << density
                    << " factor=" << factor 
                    << " elem vec = " << elemVec.ToString();
  }

  void LinNeumannInt::CalcElemVector( Vector<Complex> & elemVec,
                                      EntityIterator& ent) {
    Double value, phase, realPart, imagPart;
    Vector<Double> helpVec;
    
    // compute element vector
    PrepareElemVector( helpVec, ent );
    elemVec.Resize(helpVec.GetSize());

    //Elem * ptVolElem = actElem_->ptVolElem1;
    //RegisterSurfElemMidPoint( mHandle_, ent.GetSurfElem(), ptVolElem );

    // register global coordinates of nodes
    const StdVector<UInt> & surfConnect = ent.GetSurfElem()->connect;
    Matrix<Double> surfCoordMat;
    domain->GetGrid()->GetElemNodesCoord( surfCoordMat, surfConnect, coordUpdate_ );
    Vector<Double> nodeCoord;
    nodeCoord.Resize(surfCoordMat.GetNumRows());

    // fetch global coordinate system
    CoordSystem * coosy = domain->GetCoordSystem();

    // set info for ResultCache, which handles the input function
    // (reads data from input files, i.e. hdf5, unv, etc.)
    // in math parser expressions
    ResultCache::SetInfo( ResultCache::OUT_REAL,
                          0,
                          ent.GetName(),
                          quantity_ );

    for ( UInt i=0, n=helpVec.GetSize(); i<n; ++i )
    {
      for (UInt j=0, m=nodeCoord.GetSize(); j<m; ++j )
      {
        nodeCoord[j] = surfCoordMat[j][i];
      }
      mParser_->SetCoordinates( mHandle_, *coosy, nodeCoord );
      
      ResultCache::SetIndex(surfConnect[i]);
      
      // evaluate value and phase for current element
      ResultCache::SetOutputType(ResultCache::OUT_AMPL);
      mParser_->SetExpr( mHandle_, amplitude_ );
      value = mParser_->Eval( mHandle_ );
      // we assume the surface normal points out of domain,
      //  but we want to take deflections into the domain positive
      value *= - 1.0;

      ResultCache::SetOutputType(ResultCache::OUT_PHASE);
      mParser_->SetExpr( mHandle_, phase_ );
      phase = mParser_->Eval( mHandle_ );

      // Note: Since phase is in degrees, we have to transform it into radians
      realPart = value * cos( phase / 180 * PI );
      imagPart = value * sin( phase / 180 * PI );

      // multiply element vector with complex factor
      Complex factor(realPart, imagPart);
      elemVec[i] = helpVec[i] * factor;
    }
  }


} // end of namespace
