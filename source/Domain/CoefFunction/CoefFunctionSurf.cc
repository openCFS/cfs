#include "CoefFunctionSurf.hh"
#include "FeBasis/FeSpace.hh"
#include "FeBasis/H1/H1Elems.hh"
#include "FeBasis/HCurl/HCurlElems.hh"
#include "Driver/SolveSteps/BaseSolveStep.hh"
#include "PDE/SinglePDE.hh"
#include "Domain/Results/ResultInfo.hh"
#include <type_traits>
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
    numDofs_ = surfInfo->dofNames.size();
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

void CoefFunctionSurf::AddEntity(shared_ptr<EntityList> entity) {  
  switch (entity->GetType()) {
    case EntityList::SURF_ELEM_LIST:
      if (entities_.find(entity->GetName()) ==entities_.end()) {
          entities_[entity->GetName()] = entity;
      }
      break;

    case EntityList::NODE_LIST:
      if (entities_.find(entity->GetName()) == entities_.end()) {
          entities_[entity->GetName()] =
            entity->GetGrid()->GetEntityList(EntityList::SURF_ELEM_LIST, entity->GetName());
      }
      break;

    case EntityList::NAME_LIST:
      {
        EntityIterator nameIt = entity->GetIterator();
        for (nameIt.Begin(); !nameIt.IsEnd(); nameIt++) {
          if (entities_.find(nameIt.GetName()) == entities_.end()) {
            entities_[nameIt.GetName()] =
                nameIt.GetGrid()->GetEntityList(EntityList::SURF_ELEM_LIST, nameIt.GetName());      
          }
        }
      }
      break;

    default:
      EXCEPTION("Unsupported list type " << entity->GetType());
      break;
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
template<class FE, class DATA_TYPE>
CoefFunctionSurfVWP<FE,DATA_TYPE>::CoefFunctionSurfVWP(PtrCoefFct matCoef,
                                                   shared_ptr<ResultInfo> surfInfo,
                                                   Grid* ptGrid,
                                                   shared_ptr<BaseFeFunction> feFnc)
: CoefFunctionSurf(false, 1.0, surfInfo)
{
  matCoef_ = matCoef;
  ptGrid_ = ptGrid;
  FeFunction_ = feFnc; 
  cacheStep_ = 0;
  isComplex_ = false;
  if (std::is_same<DATA_TYPE, Complex>::value)
    isComplex_ = true;
}

template<class FE, class DATA_TYPE>
CoefFunctionSurfVWP<FE,DATA_TYPE>::~CoefFunctionSurfVWP() {
}

// Returns the total force summing up all nodal forces over a group
template<class FE, class DATA_TYPE>
void CoefFunctionSurfVWP<FE,DATA_TYPE>::GetTotalForce(const std::string & entityName,
                                                         Vector<DATA_TYPE> & totalForce)
{
    UpdateCache();

    if (totalForces_.find(entityName) == totalForces_.end()) {
      EXCEPTION("No data available for entity '" << entityName << "'");
    }

    totalForce = totalForces_[entityName];
}

template<class FE, class DATA_TYPE>
void CoefFunctionSurfVWP<FE,DATA_TYPE>::GetVector(Vector<DATA_TYPE>& coefVec,
                                                     const LocPointMapped& lpm ) {
  assert(this->dimType_ == VECTOR);

  UpdateCache();
  
  coefVec.Resize(numDofs_);
  coefVec.Init();

  BaseFE * fe;
  if (std::is_same<FE, FeH1>::value) {
    //nodal H1 elements
    fe =  static_cast<FeH1*>(FeFunction_->GetFeSpace()->GetFe(lpm.ptEl->elemNum));
    
    Vector<Double> shapeFncs;
    fe->GetShFnc(shapeFncs, lpm.lp, lpm.ptEl);

    const StdVector<UInt> & connect = lpm.ptEl->connect;
    UInt numElemNodes = connect.size();
    for (UInt node = 0; node < numElemNodes; ++node) {
      if (nodalForces_.find(connect[node]) == nodalForces_.end()) {
        EXCEPTION("Nodal force of node " << connect[node] << " not found in cache.");
      }
      coefVec += nodalForces_[connect[node]] * shapeFncs[node];
    }
    coefVec.ScalarDiv(lpm.shapeMap->CalcVolume(true));
  }
  else if (std::is_same<FE, FeHCurl>::value) {
    //edge H(curl) elements    
    fe =  static_cast<FeHCurl*>(FeFunction_->GetFeSpace()->GetFe(lpm.ptEl->elemNum));
    
    Matrix<Double> shapeFncs;
    fe->GetCurlShFnc(shapeFncs, lpm, lpm.ptEl);

    const StdVector<UInt> & connect = lpm.ptEl->connect;
    UInt numElemNodes = connect.size();
    for (UInt node = 0; node < numElemNodes; ++node) {
      if (nodalForces_.find(connect[node]) == nodalForces_.end()) {
        EXCEPTION("Nodal force of node " << connect[node] << " not found in cache.");
      }

      const Vector<DATA_TYPE> & nodeForce = nodalForces_[connect[node]];
      for (UInt dof = 0; dof < numDofs_; ++dof) {
        coefVec[dof] += nodeForce[dof] * shapeFncs[dof][node];
     }
    }
    coefVec.ScalarDiv(lpm.shapeMap->CalcVolume(true));
  }
  else  
    EXCEPTION("CoefFunctionSurfVWP<FE,DATA_TYPE>::GetVector: FE-Type unknown");
}

// Update the cache with nodal forces 
template<class FE, class DATA_TYPE> 
void CoefFunctionSurfVWP<FE,DATA_TYPE>::UpdateCache() {
  if (FeFunction_->GetPDE()->GetSolveStep()->GetActStep() == cacheStep_) {
    return;
  }

  std::map< std::string, shared_ptr<EntityList> >::iterator srcIt = entities_.begin();
  std::map< std::string, shared_ptr<EntityList> >::iterator srcEnd = entities_.end();

  StdVector<UInt> surfNodeList;
  StdVector<const Elem*> elemList;
  StdVector< StdVector<UInt> > elemNodeToCouplingNode;
  StdVector< std::vector<bool> > isBoundaryNode;
  Matrix<DATA_TYPE> elemForce;
  StdVector<RegionIdType> neighborIds(1);

  //make sure, that when array has already be used, zero it!
  for ( UInt i = 0; i< nodalForces_.size(); i++ ) {
    Vector<DATA_TYPE> & nodeForce = nodalForces_[i];
    if ( nodeForce.size() != 0) 
      nodeForce.Init();
  }

  // Loop over groups
  for (srcIt = entities_.begin(); srcIt != srcEnd; ++srcIt)  {
    totalForces_[srcIt->first].Resize(numDofs_);
    totalForces_[srcIt->first].Init(0.0);

    // get nodes belonging to the surface elements
    ptGrid_->GetNodesByRegion(surfNodeList, ptGrid_->GetRegionId(srcIt->first) );

    // get volume elements next to nodes !!!HIER 
    neighborIds[0] = this->GetVolNeighborRegionId(srcIt->second->GetRegion()); 
    ptGrid_->GetElemsNextToNodes( elemList, surfNodeList, neighborIds); 

    // get memory
    isBoundaryNode.Resize(elemList.GetSize());
    isBoundaryNode.Init();
    elemNodeToCouplingNode.Resize(elemList.GetSize());
    elemNodeToCouplingNode.Init();

    //search for boundary nodes, which are allowed to move
    for (UInt iElem = 0; iElem < elemList.GetSize(); ++iElem) {
      isBoundaryNode[iElem].resize(elemList[iElem]->connect.GetSize(), false);
      elemNodeToCouplingNode[iElem].Resize(elemList[iElem]->connect.GetSize());
      elemNodeToCouplingNode[iElem].Init();

      // Determine Boundary Nodes
      for (UInt iElemNode = 0; iElemNode < isBoundaryNode[iElem].size(); ++iElemNode) {
        for (UInt iNodes = 0; iNodes < surfNodeList.GetSize(); ++iNodes) {
          if (elemList[iElem]->connect[iElemNode] == surfNodeList[iNodes] ) {
            isBoundaryNode[iElem][iElemNode] = true;
            elemNodeToCouplingNode[iElem][iElemNode] = iNodes;
            break;
          }
        }
      }
    } // end elements

    for (UInt iElem = 0; iElem < elemList.GetSize(); ++iElem) {
      CalcElemForce(elemForce, elemList[iElem], isBoundaryNode[iElem]);

      const StdVector<UInt> & connect = elemList[iElem]->connect;
      UInt numElemNodes = connect.size();

      // Add the element force to the according coupling node
      for (UInt iElemNode = 0; iElemNode < numElemNodes; ++iElemNode) {
        Vector<DATA_TYPE> & nodeForce = nodalForces_[connect[iElemNode]];
        if (nodeForce.size() == 0) {
          nodeForce.Resize(numDofs_);
          nodeForce.Init();
        }

        for (UInt dof = 0; dof < numDofs_; ++dof) {
          nodeForce[dof] += elemForce(iElemNode, dof);
          totalForces_[srcIt->first][dof] += elemForce(iElemNode, dof);
        }
      }
    } // end elements!
  } // end groups

  // //make check
  // std::cout << "Num nodal forces: " << nodalForces_.size() << std::endl;
  // Vector<DATA_TYPE> totForce(numDofs_); totForce.Init(); 
  // for ( UInt i = 0; i< nodalForces_.size(); i++ ) {
  //   Vector<DATA_TYPE> & nodeForce = nodalForces_[i];
  //   if ( nodeForce.size() != 0) {
  //     for ( UInt j = 0; j < numDofs_; j++ ) {
  //       totForce[j] +=  nodeForce[j];
  //     }
  //   }
  // }
  // std::cout << "Total Force:  " << totForce << std::endl;
  // for (srcIt = entities_.begin(); srcIt != srcEnd; ++srcIt)  {
  //   std::cout << "Real Total Force: " << totalForces_[srcIt->first] << std::endl;
  // }

  //node forces computed!
  cacheStep_ = FeFunction_->GetPDE()->GetSolveStep()->GetActStep();
}

template<class FE, class DATA_TYPE>
void CoefFunctionSurfVWP<FE,DATA_TYPE>::CalcElemForce(Matrix<DATA_TYPE>& force,
                                            const Elem* ptElement,
                                            const std::vector<bool> & isBoundaryNode)
{
    UInt numNodes = ptElement->connect.GetSize();

    Double permeability, jDet, DetdJ_dr;
    DATA_TYPE fieldSqr;
    Vector<DATA_TYPE> field, conjField;
    Vector<DATA_TYPE> JInvTimesdJdr(numDofs_), JInvTimesdJdrConj(numDofs_);
    Matrix<Double> J, JinvT, J_r_Trans;
    Matrix<Double> SpecCornerCoords(numDofs_, numNodes);

    //if complex, we have to multiply by 1/4, since we are interested in the time averaged force!
    Double preFactor = 1.0;
    if ( isComplex_ )
      preFactor = 0.25;

    // Obtain FE element from feSpace and integration scheme
    IntegOrder order;
    order.SetIsoOrder(2);

    shared_ptr<FeSpace> feSpace = FeFunction_->GetFeSpace(); 
    // Get shape map from grid
    shared_ptr<ElemShapeMap> esm = FeFunction_->GetGrid()->GetElemShapeMap(ptElement, true);

    // Get integration points
    StdVector<LocPoint> intPoints;
    StdVector<Double> weights;
    shared_ptr<IntScheme> intScheme = feSpace->GetIntScheme();
    intScheme->GetIntPoints( Elem::GetShapeType(ptElement->type),
                             IntScheme::GAUSS, order, intPoints, weights );

    force.Resize(isBoundaryNode.size(), numDofs_);
    force.Init();

    // Loop over all integration points
    LocPointMapped lpm;
    lpm.SetCheckJacobi(false); //Jacobian dterminat may be zero or negative!!!
    for (UInt i = 0; i < intPoints.GetSize(); ++i) {
      // Calculate for each integration point the LocPointMapped
      lpm.Set((intPoints)[i], esm, (weights)[i]);

      // Get field vector
      coefs_[lpm.ptEl->regionId]->GetVector(field, lpm);

      // Get permeability
      matCoef_->GetScalar(permeability, lpm);

      // Scale field by square root of material parameter
      field /= std::sqrt(permeability);
      fieldSqr = field.Inner();

      //in complex case, we also need the conjugate complex
      if ( isComplex_ ) 
        conjField = field.Conj();

      jDet = lpm.jacDet;
      J = lpm.jac;
      lpm.jacInv.Transpose(JinvT);

      for (UInt nNode = 0; nNode < isBoundaryNode.size(); ++nNode) {
        // loop over all dimension
        for (UInt idim = 0; idim < numDofs_; ++idim) {
          // form SpecCornerCoords-Array
          SpecCornerCoords.Init();
          if (isBoundaryNode[nNode]) {
            SpecCornerCoords[idim][nNode] = 1.0;
          }
          else {
            break;
          }

          lpm.Set(intPoints[i], esm, weights[i], SpecCornerCoords);

          // calculate dJ_dr
          DetdJ_dr = CalcDetJDr(J, lpm.jac, lpm);
          lpm.jac.Transpose(J_r_Trans);
          JInvTimesdJdr = JinvT * J_r_Trans * field;
          if ( isComplex_ )
            JInvTimesdJdrConj = JInvTimesdJdr.Conj(); //JinvT * J_r_Trans * conjField;       

          if ( isComplex_ ) {
            force(nNode,idim) -= preFactor * ( ( conjField.Inner(JInvTimesdJdr)
                                                 + field.Inner(JInvTimesdJdrConj)) * jDet
                                  - fieldSqr * DetdJ_dr) * weights[i];
          }
          else {
            force(nNode,idim) -= (field.Inner(JInvTimesdJdr) * jDet
                                  - 0.5 * fieldSqr * DetdJ_dr) * weights[i];
          }

        } // loop over dimension
      } // loop over boundary nodes
    } // loop over integration points

}

template<class FE, class DATA_TYPE>
Double CoefFunctionSurfVWP<FE,DATA_TYPE>::CalcDetJDr(const Matrix<Double> &J,
                                              const Matrix<Double> &dJ_dr,
                                              const LocPointMapped &lpm) {
  Double det;

  if (J.GetNumRows() == 2) {
    det = dJ_dr[0][0] * J[1][1] + dJ_dr[1][1] * J[0][0]
          - dJ_dr[0][1] * J[1][0] - dJ_dr[1][0] * J[0][1];
    // has to be done ????????
    // if (lpm.shapeMap->IsAxi()) {
    //   det *= 2 * M_PI * lpm.globPoint[0];
    // }
  }
  else {
    det =  dJ_dr[0][0] * (J[1][1]*J[2][2] - J[1][2]*J[2][1])
          -  dJ_dr[0][1] * (J[1][0]*J[2][2] - J[1][2]*J[2][0])
          +  dJ_dr[0][2] * (J[1][0]*J[2][1] - J[1][1]*J[2][0])
          -  dJ_dr[1][0] * (J[0][1]*J[2][2] - J[0][2]*J[2][1])
          +  dJ_dr[1][1] * (J[0][0]*J[2][2] - J[0][2]*J[2][0])
          -  dJ_dr[1][2] * (J[0][0]*J[2][1] - J[0][1]*J[2][0])
          +  dJ_dr[2][0] * (J[0][1]*J[1][2] - J[0][2]*J[1][1])
          -  dJ_dr[2][1] * (J[0][0]*J[1][2] - J[0][2]*J[1][0])
          +  dJ_dr[2][2] * (J[0][0]*J[1][1] - J[0][1]*J[1][0]);
  }
  return det;
}

// explicit template instantiation
template class CoefFunctionSurfVWP<FeH1, Double>;
template class CoefFunctionSurfVWP<FeHCurl, Double>;
template class CoefFunctionSurfVWP<FeH1, Complex>;
template class CoefFunctionSurfVWP<FeHCurl, Complex>;


}
