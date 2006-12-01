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
    ENTER_FCN( "LinNeumannInt::LinNeumannInt" );

    isaxi_ = isaxi;
    amplitudeStr_ = amplitudeStr;
    phaseStr_ = phaseStr;
    materialParam_ = materialParam;

    std::cerr << "In LinNeumannInt::CalcElemVector output " << std::endl
              << " coordinate     amplitude    phase " << std::endl;
    
  }

  
  LinNeumannInt::~LinNeumannInt() {
    ENTER_FCN( "LinNeumannInt::~LinNeumannInt" );

    actElem_ = NULL;
    
  }

  void LinNeumannInt::PrepareElemVector( Vector<Double> & elemVec,
                                      EntityIterator& ent) {
    ENTER_FCN( "LinNeumannInt::PrepareElemVector" );

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
    if (materialParam_ == DENSITY) {
      it->second->GetScalar(factor,DENSITY,REAL);
    }
    else if (materialParam_ == HEAT_CONDUCTIVITY) {
      Double k;
      it->second->GetScalar(k,HEAT_CONDUCTIVITY,REAL);
      factor = 1.0 / k;
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
    ENTER_FCN( "LinNeumannInt::CalcElemVector" );

    // we assume the surface normal points out of domain,
    //  but we want to take deflections into the domain positive
    Double amplitude = -1.0 * String2Double(amplitudeStr_);

    //Double phase = String2Double( actBc.phase);

    PrepareElemVector( elemVec, ent);
    elemVec *= amplitude;
  }

  void LinNeumannInt::CalcElemVector( Vector<Complex> & elemVec,
                                      EntityIterator& ent) {
    ENTER_FCN( "LinNeumannInt::CalcElemVector" );

    // In this order, because ExtractElemInfo(ent)  has to be executed at first
    Vector<Double> helpVec;
    PrepareElemVector(helpVec, ent);

    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent );

    // get global coordinate system and math parser
    CoordSystem * coosy = domain->GetCoordSystem();
    MathParser * parser = domain->GetMathParser();

    // Get midpoint of element
    Vector<Double> locMidPoint, locMidPointVol, globMidPointVol;

    const StdVector<UInt> & surfConnect = actElem_->connect;
    const StdVector<UInt> & volConnect = ptVolElem_->connect;

    actElem_->ptElem->GetCoordMidPoint(locMidPoint);
    ptVolElem_->ptElem->GetLocalIntPoints4Surface(surfConnect, volConnect,
                                                  locMidPoint, locMidPointVol);

    // Map to global coordinate system
    Matrix<Double> volCoordMat;
    domain->GetGrid()->GetElemNodesCoord( volCoordMat, volConnect, coordUpdate_ ); 
    ptVolElem_->ptElem->Local2GlobalCoord( globMidPointVol, locMidPointVol,
                                           volCoordMat, ptVolElem_ );

    // Update variables for mathParser
    parser->SetCoordinates( mHandle_, *coosy, globMidPointVol );

    // Get amplitude and phase
    Double amplitude;
    parser->SetExpr( mHandle_, amplitudeStr_ );
    amplitude = parser->Eval( mHandle_ );
    // we assume the surface normal points out of domain,
    //  but we want to take deflections into the domain positive
    amplitude *= - 1.0;

    Double phase;
    parser->SetExpr( mHandle_, phaseStr_ );
    phase = parser->Eval( mHandle_ );

    // Note: Since phase is in °(grad), we have to transform it into
    //        rad-value
    Double realPart = amplitude * cos(phase/180*PI);
    Double imagPart = amplitude * sin(phase/180*PI);

    Complex val = Complex( realPart, imagPart); 

    elemVec =  helpVec * val;

    std::cerr << globMidPointVol[0] 
              << "   " << amplitude 
              << "   " << phase << std::endl;

  }


} // end of namespace
