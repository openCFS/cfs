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
  // set the dimension based on the result (if possible)
  if( surfInfo ) {
    dim_= surfInfo->GetDimDof();
  } else {
    dim_ = 0;
  }
  
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
  if( surfInfo) {
    resultType_ = surfInfo->resultType;
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
  // create local point for surface
  LocPointMapped surfLpm(lpm);
  surfLpm.SetSurfInfo( regions_);
  RegionIdType region = surfLpm.lpmVol->ptEl->regionId;

  //Confirm that output variable is defined as a vector
  assert(this->dimType_ == VECTOR);

  //For tensor or vector volume coefficient function
  if (coefs_[region]->GetDimType() == TENSOR || coefs_[region]->GetDimType() == VECTOR ){
    //normal mapping (from Tensor)
    if( mapNormal_ ) {
      Vector<Double> coefTens;
      coefs_[region]->GetVector(coefTens, *surfLpm.lpmVol );
      MapTensorNormal<Double>( coefVec, coefTens, surfLpm.normal);
    //get vector without normal mapping
    } else {
      coefs_[region]->GetVector(coefVec, *surfLpm.lpmVol );
    }
  }

  //normal mapping from scalar CoefFunction
  else if (mapNormal_ && coefs_[region]->GetDimType() == SCALAR){
    Double coefScal;
    coefs_[region]->GetScalar(coefScal, *surfLpm.lpmVol);
    MapScalNormal<Double>( coefVec, coefScal, surfLpm.normal);
  }

  else{
    EXCEPTION("Case not implemented.");
  }

  coefVec *= factor_;
}


void CoefFunctionSurf::GetVector(Vector<Complex>& coefVec, 
                               const LocPointMapped& lpm ) {
  // create local point for surface
  LocPointMapped surfLpm(lpm);
  surfLpm.SetSurfInfo( regions_);
  RegionIdType region = surfLpm.lpmVol->ptEl->regionId;

  //Confirm that output variable is defined as a vector
  assert(this->dimType_ == VECTOR);

  //For tensor or vector volume coefficient function
  if (coefs_[region]->GetDimType() == TENSOR || coefs_[region]->GetDimType() == VECTOR ){
    //normal mapping (from Tensor)
    if( mapNormal_ ) {
      Vector<Complex> coefTens;
      coefs_[region]->GetVector(coefTens, *surfLpm.lpmVol );
      MapTensorNormal<Complex>( coefVec, coefTens, surfLpm.normal);
    //get vector without normal mapping
    } else {
      coefs_[region]->GetVector(coefVec, *surfLpm.lpmVol );
    }
  }

  //normal mapping from scalar CoefFunction
  else if (mapNormal_ && coefs_[region]->GetDimType() == SCALAR){
    Complex coefScal;
    coefs_[region]->GetScalar(coefScal, *surfLpm.lpmVol);
    MapScalNormal<Complex>( coefVec, coefScal, surfLpm.normal);
  }

  else{
    EXCEPTION("Case not implemented.");
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
  assert(this->dimType_ == VECTOR);
  // determine size from secondary CoefFunction: only possible if there is no surface mapping
  //UInt dim;
  //for (auto const& it : coefs_)
  //{
  //  dim = (it.second)->GetVecSize();
  //}
  return dim_;
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
void CoefFunctionSurf::MapScalNormal( Vector<TYPE>& ret, TYPE& scal,
                                     const Vector<Double>& normal ) {
//  ret = normal * scal;
  ret.Resize(normal.GetSize());
  for(UInt i=0; i<normal.GetSize();++i)
    ret[i] = normal[i]*scal;
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
  } else if (tensor.GetSize() == 9) {
    // 3D non-symmetric tensor (defined row-wise, currenty only used in WaterWavePDE)
    ret.Resize(3);
    ret[0] = tensor[0]*n[0] + tensor[1]*n[1] + tensor[2]*n[2];
    ret[1] = tensor[3]*n[0] + tensor[4]*n[1] + tensor[5]*n[2];
    ret[2] = tensor[6]*n[0] + tensor[7]*n[1] + tensor[8]*n[2];
  } else if (tensor.GetSize() == 2 ) {
    // 2D cross product case (currenty only used in WaterWavePDE)
    ret.Resize(2); // needs to be 2 since the result functor wants to evaluate a 2d vector in a 2d analysis
    ret[0] = tensor[0]*n[0] + tensor[1]*n[1]; // this is the z-component (only physical component)
    ret[1] = 0.0; // dummy component = always zero
  }
  else {
    EXCEPTION( "Case not implemented" << tensor.GetSize());
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
  LocPointMapped surfLpm(lpm); //Permeabilitätsänderung

  //get regionId from surfRegion
  RegionIdType surRegId = surfLpm.ptEl->regionId;
  RegionIdType volNeighborRegionId = neighborRegionId_[surRegId];
  surfLpm.SetSurfInfo( regions_, volNeighborRegionId);

  //get magnetic flux density
  Vector<Double> Bvec;
  coefs_[volNeighborRegionId]->GetVector(Bvec, *surfLpm.lpmVol );

  //get permeability
  Double permeability, matFactor;
  std::map<RegionIdType,PtrCoefFct > permFncs = matCoef_[MAG_ELEM_PERMEABILITY]->GetRegionCoefs();
  permFncs[volNeighborRegionId]->GetScalar(permeability, *surfLpm.lpmVol );
  matFactor = 1.0 / permeability;

  //compute: factors * ( (B*n)*B - 1/2*B^2*n )
  Vector<Double> normalVec = surfLpm.normal;
  Double Bn = Bvec * normalVec; // normal component of B
  Double B2 = Bvec.Inner(); // square of B
  Vector<Double> BnB = Bvec; BnB.ScalarMult(Bn); // B * normal component
  Vector<Double> B2n = normalVec; B2n.ScalarMult(0.5*B2); // n * B^2/2
  coefVec = BnB - B2n;
  coefVec.ScalarMult(factor_ * matFactor); // don't forget factors

  if ( resultType_ == MAG_NORMALFORCE_MAXWELL_DENSITY
       || resultType_ == MAG_TANGENTIALFORCE_MAXWELL_DENSITY) {
    //save the total surface force density vector
    Vector<Double> fVecSurface(coefVec);

    Double fNormal = coefVec * normalVec;
    coefVec = normalVec; coefVec.ScalarMult(fNormal); //vector in normal direction

    if (resultType_ == MAG_TANGENTIALFORCE_MAXWELL_DENSITY) {
      //coefVec contains already the force density vector in normal direction
      coefVec *= -1.0;
      coefVec += fVecSurface;
    }
  }
}


void CoefFunctionSurfMaxwell::GetVector(Vector<Complex>& coefVec,
                               const LocPointMapped& lpm ) {

  assert(this->dimType_ == VECTOR);// surface force vector

  // create local point for surface
  LocPointMapped surfLpm(lpm);//Permeabilitätsänderung

  //get regionId from surfRegion
  RegionIdType surRegId = surfLpm.ptEl->regionId;
  RegionIdType volNeighborRegionId = neighborRegionId_[surRegId];
  surfLpm.SetSurfInfo( regions_, volNeighborRegionId);

  // normal unit vector
  Vector<Double> normalVec = surfLpm.normal;
  Vector<Complex> normalVecC(normalVec.GetSize());
  for ( UInt i=0; i<normalVec.GetSize(); i++) {
    normalVecC[i] = Complex(normalVec[i], 0.0);
  }

  //get magnetic flux density
  Vector<Complex> Bvec;
  Vector<Complex> conjBvec;
  coefs_[volNeighborRegionId]->GetVector(Bvec, *surfLpm.lpmVol );
  conjBvec = Bvec.Conj();

  //get permeability
  Double permeability, matFactor;
  std::map<RegionIdType,PtrCoefFct > permFncs = matCoef_[MAG_ELEM_PERMEABILITY]->GetRegionCoefs();
  permFncs[volNeighborRegionId]->GetScalar(permeability, *surfLpm.lpmVol );
  matFactor = 1.0 / permeability;

  //compute B* \cdot n; Bn = B \cdot n; compute B^2
  Complex Bn1 = conjBvec * normalVecC; // normal component of B
  Complex Bn2 = Bvec * normalVecC;
  Complex B2 = Bvec.Inner(); // square of B

  //compute BnB1 = (B* \cdot n)*B; BnB2 = (B \cdot n)*B*
  Vector<Complex> BnB1;
  Vector<Complex> BnB2;
  BnB1 = Bvec; BnB1.ScalarMult(0.25*Bn1);
  BnB2 = conjBvec; BnB2.ScalarMult(0.25*Bn2);

  Vector<Complex> B2n(normalVecC);
  B2n.ScalarMult(0.25*B2); // n * B^2/2

  //compute force vector on the surface: (1/mu) * ( 1/4 ((B* \cdot n) B + (B \cdot n) B*) - 1/4*(B \cdot B*) n )
  coefVec = BnB1 + BnB2 - B2n;
  coefVec.ScalarMult(factor_ * matFactor);

  if ( resultType_ == MAG_NORMALFORCE_MAXWELL_DENSITY
       || resultType_ == MAG_TANGENTIALFORCE_MAXWELL_DENSITY) {
    //save the total surface force density vector
    Vector<Complex> fVecSurface(coefVec);

    //compute surface force in normal direction
    Complex fNormal = coefVec * normalVecC;
    coefVec = normalVecC;
    coefVec.ScalarMult(fNormal);

    if (resultType_ == MAG_TANGENTIALFORCE_MAXWELL_DENSITY) {
      //coefVec contains already the force density vector in normal direction
      coefVec.ScalarMult(-1.0);
      coefVec += fVecSurface;
    }
  }
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
