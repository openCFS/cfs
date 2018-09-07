#include "CoefFunctionSurf.hh"


#include "Domain/Results/ResultInfo.hh"
namespace CoupledField {


CoefFunctionSurf::CoefFunctionSurf( bool mapNormal, 
                                    Double factor,
                                    shared_ptr<ResultInfo> surfInfo ) 
: CoefFunction() {


  // not sure about the following one
  dependType_ = GENERAL;
  isAnalytic_ = false;
  isComplex_ =  false;
  mapNormal_ = mapNormal;
  factor_ = factor;
  
  if( this->mapNormal_) {
    if( !surfInfo) {
      EXCEPTION( "The resultinfo object must be set in case of normal mapping!")
    }
    switch(surfInfo->entryType ) {
      case ResultInfo::SCALAR:
        this->dimType_ = SCALAR;
        break;
      case ResultInfo::VECTOR:
        this->dimType_ = VECTOR;
        break;
      case ResultInfo::TENSOR:
        EXCEPTION( "A normal-mapped result can not be of type TENSOR")
        break;
      default:
        EXCEPTION("Unhandled case for normal mapping")
    }
  } else {
    // In this case, we simply take the dimensionality of the first
    // volume coef in the method ::AddVolumeCoef
    this->dimType_ = NO_DIM;
  }
}



void CoefFunctionSurf::AddVolumeCoef( RegionIdType region, PtrCoefFct coef ) {
  
  // check, if there are already volume coeffiencts set
  if (coefs_.size() ) {
    std::map<RegionIdType, PtrCoefFct>::iterator it = coefs_.begin();
    for(;  it != coefs_.end(); ++it )  {
      if( it->second->GetDimType() != coef->GetDimType() ) {
        EXCEPTION( "CoefFunctions for surface CoefFunction have " <<
                   "inconsistent dimension type" );
      }
    }
  } else {
    isComplex_ = coef->IsComplex();
    if( !this->mapNormal_)
      dimType_ = coef->GetDimType();
  }
  
  // adjust dependency type
  dependType_ = std::max(this->GetDependency(), 
                         coef->GetDependency());
  regions_.insert(region);
  coefs_[region] = coef;
}


void CoefFunctionSurf::SetVolumeCoefs( std::map<RegionIdType, PtrCoefFct> coefs ) {
  // ensure that at least one entry is present
  assert( coefs.size() != 0);
  std::map<RegionIdType, PtrCoefFct>::iterator it = coefs.begin();

  // loop over all entries and set following variables
  for( ; it != coefs.end(); ++it ) {
    AddVolumeCoef( it->first, it->second);
  }
}

CoefFunctionSurf::~CoefFunctionSurf() {

}


void CoefFunctionSurf::GetTensor(Matrix<Double>& coefMat, 
                            const LocPointMapped& lpm ){
  assert(this->dimType_ == TENSOR);
}

void CoefFunctionSurf::GetTensor(Matrix<Complex>& coefMat, 
                            const LocPointMapped& lpm ){
  assert(this->dimType_ == TENSOR);
}


void CoefFunctionSurf::GetVector(Vector<Double>& coefVec, 
                                 const LocPointMapped& lpm ) {
  assert(this->dimType_ == VECTOR);

  // create local point for surface
  LocPointMapped surfLpm(lpm);
  surfLpm.SetSurfInfo( regions_);
  RegionIdType region = surfLpm.lpmVol->ptEl->regionId;
  if( mapNormal_ ) {
    Vector<Double> coefTens;
    coefs_[region]->GetVector(coefTens, *surfLpm.lpmVol );
    MapTensorNormal<Double>( coefVec, coefTens, surfLpm.normal);
  } else {
    coefs_[region]->GetVector(coefVec, *surfLpm.lpmVol );
  }
  coefVec *= factor_;
}

void CoefFunctionSurf::GetVector(Vector<Complex>& coefVec, 
                               const LocPointMapped& lpm ) {
  assert(this->dimType_ == VECTOR);

  // create local point for surface
  LocPointMapped surfLpm(lpm);
  surfLpm.SetSurfInfo( regions_);
  RegionIdType region = surfLpm.lpmVol->ptEl->regionId;
  if( mapNormal_ ) {
    Vector<Complex> coefTens;
    coefs_[region]->GetVector(coefTens, *surfLpm.lpmVol );
    MapTensorNormal<Complex>( coefVec, coefTens, surfLpm.normal);
  } else {
    coefs_[region]->GetVector(coefVec, *surfLpm.lpmVol );
  }
  coefVec *= factor_;
}


void CoefFunctionSurf::GetScalar(Double& coefScalar, 
                                 const LocPointMapped& lpm ){
  assert(this->dimType_ == SCALAR);

  // create local point for surface
  LocPointMapped surfLpm(lpm);
  surfLpm.SetSurfInfo( regions_);
  RegionIdType region = surfLpm.lpmVol->ptEl->regionId;
  if( mapNormal_ ) {
    Vector<Double> coefVec;
    coefs_[region]->GetVector(coefVec, *surfLpm.lpmVol );
    if (coefVec.GetSize() == surfLpm.normal.GetSize())
      coefScalar = coefVec * surfLpm.normal;
    else // this happens in 2.5D case where 'coefVec' and 'normal' have 3 and 2 components, respectively
      coefScalar = coefVec[0]*surfLpm.normal[0] + coefVec[1]*surfLpm.normal[1];
  } else {
    coefs_[region]->GetScalar(coefScalar, *surfLpm.lpmVol );
  }
  coefScalar *= factor_;
  //std::cout << "Value: " << coefScalar << std::endl;
}

void CoefFunctionSurf::GetScalar(Complex& coefScalar, 
                                 const LocPointMapped& lpm ){
  assert(this->dimType_ == SCALAR);

  // create local point for surface
  LocPointMapped surfLpm(lpm);
  surfLpm.SetSurfInfo( regions_);
  RegionIdType region = surfLpm.lpmVol->ptEl->regionId;
  if( mapNormal_ ) {
    Vector<Complex> coefVec;
    coefs_[region]->GetVector(coefVec, *surfLpm.lpmVol );
    if (coefVec.GetSize() == surfLpm.normal.GetSize())
      coefScalar = coefVec * surfLpm.normal;
    else // this happens in 2.5D case where 'coefVec' and 'normal' have 3 and 2 components, respectively
      coefScalar = coefVec[0]*surfLpm.normal[0] + coefVec[1]*surfLpm.normal[1];
  } else {
    coefs_[region]->GetScalar(coefScalar, *surfLpm.lpmVol );
  }
  coefScalar *= factor_;
}

UInt CoefFunctionSurf::GetVecSize() const {
  return 0;
}


void CoefFunctionSurf::GetTensorSize( UInt& numRows, UInt& numCols ) const {

}

std::string  CoefFunctionSurf::ToString() const {
  std::stringstream ret;
    ret << "CoefFunctionSurf defined on:\n";
    std::map<RegionIdType,PtrCoefFct >::const_iterator it = coefs_.begin();
    for( ; it != coefs_.end(); ++it ) {
      ret << "regionId " << it->first << ", value:" << it->second->ToString() << std::endl;
    }
  return ret.str();
}


template<typename TYPE>
void CoefFunctionSurf::MapVecNormal( TYPE& ret, const Vector<TYPE>& vec,
                                     const Vector<Double>& normal ) {
  ret = vec * normal;
}

//! Mapping operation for tensors in Voigt notation
template<typename TYPE>
void CoefFunctionSurf::MapTensorNormal( Vector<TYPE>& ret, const Vector<TYPE>& tensor,
                                        const Vector<Double>& n ) {
  if( tensor.GetSize() == 6 ) {
    // FULL 3D case
    ret.Resize(3);
    ret[0] = tensor[0]*n[0] + tensor[4]*n[2] + tensor[5]*n[1];
    ret[1] = tensor[1]*n[1] + tensor[3]*n[2] + tensor[5]*n[0];
    ret[2] = tensor[2]*n[2] + tensor[3]*n[1] + tensor[4]*n[0]; 
  } else if( tensor.GetSize() == 4 || tensor.GetSize() == 3 ) {
    // 2D / AXI-Symmetric case
    ret.Resize(2);
    ret[0] = tensor[0]*n[0] + tensor[2]*n[1];
    ret[1] = tensor[1]*n[1] + tensor[2]*n[0];
  } else {
    EXCEPTION( "Case not implemented" );
  }
}
  
//============= Maxwell's stress tensor =======================================================
CoefFunctionSurfMaxwell::CoefFunctionSurfMaxwell( bool mapNormal,
		                            std::map<SolutionType, shared_ptr<CoefFunctionMulti> > matCoefs,
									Grid* ptGrid,
									Double factor,
									shared_ptr<ResultInfo> surfInfo)
: CoefFunctionSurf(mapNormal, factor, surfInfo) {


  // not sure about the following one
  matCoef_ = matCoefs;
  ptGrid_ = ptGrid;

}

void CoefFunctionSurfMaxwell::GetVector(Vector<Double>& coefVec,
                                 const LocPointMapped& lpm ) {
  assert(this->dimType_ == VECTOR);// Richtungsvektor

  // create local point for surface
  LocPointMapped surfLpm(lpm);//Permeabilitätsänderung

  //get regionId from surfRegion
  RegionIdType surRegId = surfLpm.ptEl->regionId;
  RegionIdType volNeighborRegionId = neighborRegionId_[surRegId];
  surfLpm.SetSurfInfo( regions_, volNeighborRegionId);

  //RegionIdType region = surfLpm.lpmVol->ptEl->regionId;
  //std::string regionName = ptGrid_->GetRegion().ToString(region);
  //std::cout << "RegName: " << regionName << std::endl;

  //get magnetic flux density
  Vector<Double> Bvec;
  coefs_[volNeighborRegionId]->GetVector(Bvec, *surfLpm.lpmVol );
  //std::cout << "Bvec:  " << Bvec << std::endl;

  //get permeability
  Double permeability, matFactor;
  std::map<RegionIdType,PtrCoefFct > permFncs = matCoef_[MAG_ELEM_PERMEABILITY]->GetRegionCoefs();
  permFncs[volNeighborRegionId]->GetScalar(permeability, *surfLpm.lpmVol );
  matFactor = 1.0 / permeability;

  //compute: factors * ( (B*n)*B - 1/2*B^2*n )
  Double Bn = Bvec * surfLpm.normal; // normal component of B
  Double B2 = Bvec.Inner(); // square of B
  Vector<Double> BnB = Bvec; BnB.ScalarMult(Bn); // B * normal component
  Vector<Double> B2n = surfLpm.normal; B2n.ScalarMult(0.5*B2); // n * B^2/2
  coefVec = BnB - B2n;
  coefVec.ScalarMult(factor_ * matFactor); // don't forget factors
}

void CoefFunctionSurfMaxwell::GetVector(Vector<Complex>& coefVec,
                               const LocPointMapped& lpm ) {
  assert(this->dimType_ == VECTOR);

  EXCEPTION("SurfMaxwell for Harmonic Analysis not implemented");

}

CoefFunctionSurfMaxwell::~CoefFunctionSurfMaxwell() {

}


//===================Virtual Work Principle =================================================
CoefFunctionSurfVWP::CoefFunctionSurfVWP( bool mapNormal,
		                            std::map<SolutionType, shared_ptr<CoefFunctionMulti> > matCoefs,
									Double factor,
									shared_ptr<ResultInfo> surfInfo)
: CoefFunctionSurf(mapNormal, factor, surfInfo) {


  // not sure about the following one
  matCoef_ = matCoefs;
}

void CoefFunctionSurfVWP::GetVector(Vector<Double>& coefVec,
                                    const LocPointMapped& lpm ) {
  assert(this->dimType_ == VECTOR);

  //get magnetic flux density
  Vector<Double> Bvec;
  RegionIdType surRegId = lpm.ptEl->regionId;
  //RegionIdType volNeighborRegionId = neighborRegionId_[surRegId];
  coefs_[surRegId]->GetVector(Bvec, lpm );

  //get permeability
  Double permeability;
  std::map<RegionIdType,PtrCoefFct > permFncs = matCoef_[MAG_ELEM_PERMEABILITY]->GetRegionCoefs();
  permFncs[surRegId]->GetScalar(permeability, lpm );

  coefVec = Bvec / std::sqrt( permeability );
}

void CoefFunctionSurfVWP::GetVector(Vector<Complex>& coefVec,
                               const LocPointMapped& lpm ) {
  assert(this->dimType_ == VECTOR);

  EXCEPTION("CoefFunctionSurfVWP for Harmonic Analysis not implemented");

}

CoefFunctionSurfVWP::~CoefFunctionSurfVWP() {
}

}
