#include <boost/assign/list_of.hpp>

#include "Utils/tools.hh"
#include "ElemShapeMap.hh"
#include "Domain/Domain.hh"
#include "Domain/Mesh/Grid.hh"
#include "FeBasis/H1/H1ElemsLagExpl.hh"

#include "DataInOut/Logging/LogConfigurator.hh"

#ifdef OPENMP
#include <omp.h>
#endif

namespace CoupledField {
DEFINE_LOG(locPoint, "ElemShapeMap.LocPoint");
DEFINE_LOG(locPointMapped, "ElemShapeMap.LocPointMapped");
DEFINE_LOG(elemShapeMap, "ElemShapeMap.ElemShapeMap");
DEFINE_LOG(lagrangeElemShapeMap, "ElemShapeMap.LagrangeElemShapeMap");

// ===========================================================================
//  C L A S S  LocPoint
// ===========================================================================

LocPoint::LocPoint(const Vector<Double>& vec) {
  coord = vec;
}

std::ostream& operator <<(std::ostream& out, const LocPoint& lp) {
  out << "number = " << lp.number << ", ";
  out << "vector: " << lp.coord.ToString();
  return out;
}

// ===========================================================================
//  C L A S S   LocPointMapped
// ===========================================================================
LocPointMapped::LocPointMapped() :
    ptEl(nullptr), weight(0.0), jacDet(0.0), isSurface(false), checkJacobi_(true) {

}


// ---------- Pyramid rational map helpers ----------
static inline void PyraLamAndGrads(double x, double y, double z,
                                   double lam[5], double dlam[5][3])
{
  // rho = 1 - z
  const double rho = 1.0 - z;
  const double invrho = (std::fabs(rho) > 1e-12) ? 1.0/rho : 1.0/1e-12;
  const double invrho2 = invrho*invrho;

  // f = x*y*z / rho
  const double f     = x*y*z * invrho;
  const double dfx   = y*z * invrho;
  const double dfy   = x*z * invrho;
  const double dfz   = (x*y)*invrho + x*y*z*invrho2;

  // lambda_1 ... lambda_4
  lam[0] = 0.25*((1.0+x)*(1.0+y) - z + f);   // lambda_1
  lam[1] = 0.25*((1.0-x)*(1.0+y) - z - f);   // lambda_2
  lam[2] = 0.25*((1.0-x)*(1.0-y) - z + f);   // lambda_3
  lam[3] = 0.25*((1.0+x)*(1.0-y) - z - f);   // lambda_4
  lam[4] = z;                                // lambda_5

  const double dP1dx =  (1.0+y);
  const double dP1dy =  (1.0+x);
  const double dP1dz = -1.0;

  const double dP2dx = -(1.0+y);
  const double dP2dy =  (1.0-x);
  const double dP2dz = -1.0;

  const double dP3dx = -(1.0-y);
  const double dP3dy = -(1.0-x);
  const double dP3dz = -1.0;

  const double dP4dx =  (1.0-y);
  const double dP4dy = -(1.0+x);
  const double dP4dz = -1.0;

  dlam[0][0] = 0.25*(dP1dx + dfx);
  dlam[0][1] = 0.25*(dP1dy + dfy);
  dlam[0][2] = 0.25*(dP1dz + dfz);

  dlam[1][0] = 0.25*(dP2dx - dfx);
  dlam[1][1] = 0.25*(dP2dy - dfy);
  dlam[1][2] = 0.25*(dP2dz - dfz);

  dlam[2][0] = 0.25*(dP3dx + dfx);
  dlam[2][1] = 0.25*(dP3dy + dfy);
  dlam[2][2] = 0.25*(dP3dz + dfz);

  dlam[3][0] = 0.25*(dP4dx - dfx);
  dlam[3][1] = 0.25*(dP4dy - dfy);
  dlam[3][2] = 0.25*(dP4dz - dfz);

  dlam[4][0] = 0.0;
  dlam[4][1] = 0.0;
  dlam[4][2] = 1.0;
}


/*
 * NACS ported version; main difference: compute JacobianDeterminant via CalcJDet function,
 * taking the depth of the setup into account; 
 * note: CalcJDet is also changed to the NACS version, i.e., it now contains the calculation
 * of 1d Elements in 3d but no longer the scaling for the axi-symmetric case
 */
void LocPointMapped::Set(const LocPoint& lp,
                         shared_ptr<ElemShapeMap> esm,
                         Double weight)
{
  this->shapeMap = esm;
  this->lp = lp;
  this->weight = weight;
  this->ptEl = shapeMap->GetElem();
  this->isSurface = false;

  // Calculate Jacobian, its inverse as well as determinant for local point
  jacDet = esm->CalcJDet(jac, lp, true);
  
  // The inversion can only be performed in case we have a quadratic Jacobian
  // i.e. the dimension of the element is the dimension of the grid
  if (jac.GetNumCols() == jac.GetNumRows()) {
    // == normal volume element case (2D elems in 2D, 3D elems in 3D) ===
    jac.Invert(jacInv);
  }
  else {
    // In case of surface elements, we calculate the "pseudo" inverse
    Matrix<Double> prod, tmpInv, jacTrans;
    jacTrans = Transpose(jac);
    prod = jacTrans * jac;
    prod.Invert(tmpInv);
    jacInv = (tmpInv * jacTrans);
  }

  // safety check for negative Jacobian determinant
  if (jacDet <= 0.0) {
    EXCEPTION(
        "Jacobian determinant of element " << ptEl->elemNum << " with connectivity " << ptEl->connect.ToString() 
        << " in region '" << shapeMap->GetGrid()->GetRegion().ToString(ptEl->regionId) 
        << "' is negative in local point "
        << lp.coord.ToString() << "!");
  }

  // Check, if geometry is axi-symmetric. In this case we adjust the "volume" of
  // the Jacobian and store the global point
  Vector<Double> globPoint;
  if (esm->IsAxi()) {
    esm->Local2Global(globPoint, lp);
    jacDet *= 2 * M_PI * globPoint[0];
  } else {
    globPoint.Clear();
  }
}

// same for NACS and CFS (but why does NACS not use CalcJDet here
// and why do we not scale with depth_ here?
// > it just has nevern been tested according to Jens)
// Note: added depth_ scaling (which does not do any harm as long as depth_ = default = 1.0)
// can/should be commented out if not working correctly for depth_ != 1.0
void LocPointMapped::Set(const LocPoint& lp, shared_ptr<ElemShapeMap> esm,
                         Double weight, Matrix<Double>& cornerCoord) {
  this->shapeMap = esm;
  this->lp = lp;
  this->weight = weight;
  this->ptEl = shapeMap->GetElem();
  this->isSurface = false;

  // Calculate Jacobian, its inverse as well as determinant for local point
  esm->CalcJ(this->jac, lp, cornerCoord);

  // The inversion can only be performed in case we have a quadratic Jacobian
  // i.e. the dimension of the element is the dimension of the grid
  if (jac.GetNumCols() == jac.GetNumRows()) {
    // == normal volume element case (2D elemens in 2D, 3D elems in 3D) ===
    jac.Invert(jacInv);
    jac.Determinant(jacDet);
    
  } else if (jac.GetNumRows() == 3 && jac.GetNumCols() == 2) {
    // === 2D elements in 3D === 
    Vector<Double> normal;
    normal.Resize(3);
    normal[0] = jac[1][0] * jac[2][1] - jac[2][0] * jac[1][1];
    normal[1] = jac[2][0] * jac[0][1] - jac[0][0] * jac[2][1];
    normal[2] = jac[0][0] * jac[1][1] - jac[1][0] * jac[0][1];
    jacDet = sqrt(
        normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]);
    
  } else if (jac.GetNumRows() == 3 && jac.GetNumCols() == 1) {
    // === 1D elements in 3D ===
    jacDet = sqrt( jac[0][0] * jac[0][0]
                 + jac[1][0] * jac[1][0]
                 + jac[2][0] * jac[2][0]);
    
  } else if (jac.GetNumRows() == 2) {
    // === 1D elements in 2D ===
    //see kaltenbacher, p.23, eq.(2.122)
    jacDet = sqrt(jac[0][0] * jac[0][0] + jac[1][0] * jac[1][0]);
  };
  
  // safety check for negative Jacobian determinant
  if ( checkJacobi_ ) {
	  if ( jacDet <= 0.0) {
		  EXCEPTION("Jacobian determinant of element " << ptEl->elemNum << " with connectivity " << ptEl->connect.ToString() << " in region '" << shapeMap->GetGrid()->GetRegion().ToString(ptEl->regionId) << "' is negative! The Jacobian was:\n " << jac << " Coordinates were: \n" << shapeMap->CalcVolume());
	  }
  }

  // Check, if geometry is axi-symmetric. In this case scale the
  // Jacobian determinant with 2*pi*r
  Vector<Double> globPoint;
  if (esm->IsAxi()) {
    esm->Local2Global(globPoint, lp);
    jacDet *= 2 * M_PI * globPoint[0];
  } else {
    // scale with depth; only relevant for 2d-plane
    // note that variable depth_ will be 1.0 in 3d case as it will only be set in 2d case
    jacDet *= shapeMap->GetModelDepth();
  }
}


void LocPointMapped::SetWithSurface(const LocPoint& lp, shared_ptr<ElemShapeMap> esm, const std::set<RegionIdType>& myRegions, Double weight) {
  // ------------------------------------------------------
  //  1) Set "normal" element information (Jacobian, etc.)
  // ------------------------------------------------------
  // first, call "normal" set method, also valid for volume elements
  this->Set(lp, esm, weight);
  // ------------------------------------------------------
  //  2) Perform surface-element specific tasks
  // ------------------------------------------------------
  this->SetSurfInfo(myRegions);
}


void LocPointMapped::SetSurfInfo(const std::set<RegionIdType>& myRegions, const RegionIdType volRegId) {
  // if we have not previously selected the correct volume neighbor, we have to do it now
  // check, if previously set element is the same as this one. Only assign for the first time.
  bool isSameElem = (shapeMap->GetElem() == this->ptEl && isSurface == true);
  // set flag for surface mapped element
  this->isSurface = true;
  // get surface element
  const SurfElem& surfElem = *(shapeMap->GetSurfElem());
  // shape map for assigned volume element
  shared_ptr<ElemShapeMap> esmVol;
  // pointer to volume element that gets assigned
  const Elem * ptVolElem = nullptr;
  if (!isSameElem) {
    // loop over volume element neighbors of the surface element and check,
    // if the regionId of the element is in the map "myRegions"
    for (auto it = surfElem.ptVolElems.begin(); it != surfElem.ptVolElems.end(); it++) {
      // check if element is set at all
      if (*it) {
        // check if regionId is in the "allowed" list
        if (myRegions.find((*it)->regionId) != myRegions.end()) {
          if ( volRegId == NO_REGION_ID ) {
            // assign the first vol elem if no region is specified
            ptVolElem = *it;
            break;
          }
          // assign vol elem of the specific region if specified
          else if ( (*it)->regionId == volRegId ) {
              ptVolElem = *it;
              break;
          }
        }
      }
    } // loop over volume element neighbors
    // check, if element could be found
    if (ptVolElem == nullptr)
      EXCEPTION("Could not find a suitable volume neighbor for surface element #" << surfElem.elemNum << ". ");

    // create new local point for volume element
    lpmVol.reset(new LocPointMapped());
    esmVol = shapeMap->GetGrid()->GetElemShapeMap(ptVolElem, shapeMap->IsUpdated());
  } else {
    esmVol = lpmVol->shapeMap;
  } // !isSameElem

  LocPoint lpVol;
  Vector<Double> locNormal;
  // assign the mapped integration point to the volume lpm for computing the global normal vector
  esmVol->GetLocalIntPoints4Surface(ptEl->connect, lp, lpVol, locNormal);
  // assign integration points to the (new) volume element
  lpmVol->Set(lpVol, esmVol, weight);
  // calculate global normal pointing out of new volume element
  normal = Transpose(lpmVol->jacInv) * locNormal;
  normal /= normal.NormL2();
#ifdef NDEBUG
#else
  // sanity check if the surf->vol mapping is correct
  Vector<Double> globSurfIntPoint;
  Vector<Double> globVolIntPoint;
  this->shapeMap->Local2Global(globSurfIntPoint, lp);
  esmVol->Local2Global(globVolIntPoint, lpVol);
  LOG_DBG2(locPointMapped) << "------------------------------------------------------------------------------------------------------------" << std::endl;
  LOG_DBG2(locPointMapped) "globCoordSurf: " << globSurfIntPoint << std::endl;
  LOG_DBG2(locPointMapped) "globCoordVol: " << globVolIntPoint << std::endl;
  LOG_DBG2(locPointMapped) "locCoordSurf: " << lp.coord << std::endl;
  LOG_DBG2(locPointMapped) "locCoordVol: " << lpVol.coord << std::endl << std::endl;
  LOG_DBG2(locPointMapped) "surfElemNum---:" << surfElem.elemNum << std::endl << std::endl;
  LOG_DBG2(locPointMapped) "normal loc-----------:" << locNormal << std::endl;
  LOG_DBG2(locPointMapped) "normal glob-------------------:" << normal << std::endl << std::endl;
  LOG_DBG2(locPointMapped) << "------------------------------------------------------------------------------------------------------------" << std::endl;
  for (UInt iDim = 0; iDim < globSurfIntPoint.GetSize(); ++iDim)
    assert(abs(globSurfIntPoint[iDim] - globVolIntPoint[iDim]) < EPS);
  assert(esmVol->CoordIsInsideElem(lpVol.coord, NORM_EPS));
#endif
}

void LocPointMapped::SetWithNitscheSurface(const LocPoint& lp, shared_ptr<ElemShapeMap> esm, const bool usePrimary, Double weight) {
  // ------------------------------------------------------
  //  1) Set "normal" element information (Jacobian, etc.)
  // ------------------------------------------------------
  // first, call "normal" set method, also valid for volume elements
  this->Set(lp, esm, weight);
  // ------------------------------------------------------
  //  2) Perform surface-element specific tasks
  // ------------------------------------------------------
  // set flag for surface mapped element (Set() sets it to false before)
  this->isSurface = true;
  // get surface element and cast down to mortar element
  const MortarNcSurfElem* mortarElem = dynamic_cast<const MortarNcSurfElem*>(this->shapeMap->GetSurfElem());
  // get correct surface element
  const SurfElem* surfElem = (usePrimary) ? mortarElem->ptPrimary : mortarElem->ptSecondary;
  // get correct volume element
  const Elem* volElem = surfElem->ptVolElems[0];
  // shape map for assigned volume element
  const shared_ptr<ElemShapeMap> esmVol =  this->shapeMap->GetGrid()->GetElemShapeMap(volElem, this->shapeMap->IsUpdated());
  // get shape map of the corresponding surface element
  const shared_ptr<ElemShapeMap> esmSurf = this->shapeMap->GetGrid()->GetElemShapeMap(surfElem, this->shapeMap->IsUpdated());

  //in this special case, we need to transform local2Global wrt mortarelem
  //then global2local wrt volume elem...
  // global integration point, which is the lp transformed via th intersection element
  Vector<Double> globIntPoint;
  // local integration point determined by mapping the global point via the surfElem
  // stored in a LocPoint object for lpmVol that is used to determine the normal vector
  LocPoint lpVol, lpSurf;
  // local normal vector of lpmVol that points out of the volElem in direction of the surfElem
  Vector<Double> locNormal;
  // Calculate global coordinates of integration point
  this->shapeMap->Local2Global(globIntPoint, lp);

  // -------kirill (periodic boundary conditions)-------------------------------------------------------------------------------------
  // For periodic boundary conditions, the intersection element may not coincide with the primary/secondary element.
  // Therefore, the corresponding volume element will not lie adjacent to it. For this reason
  // we must do the calculations with respect to the chosen 'surfElem' which is either primary or secondary
  // and definitely connected with 'volElem'.
  // We cannot simply map the global point lying on our surface element to the volume element: in case of mortar element
  // for p.b.c., the global point can lie out of the volume element. That's why we first map the global point form the NC-element to
  // the primary (or secondary) surface element in order to calculate the corresponding local point in the adjacent volume element.
  // The mapping is performed in accordance with the transformation between the primary and the secondary surfaces.
  if (mortarElem->transVect.GetSize() != 0 && !mortarElem->rotationCenter.GetSize()) {
    esmSurf->TranslatePointOntoSurface(mortarElem->transVect, globIntPoint);
  }
  else if (mortarElem->rotationCenter.GetSize()) {
    // in case of cake-piece projection we also need to rotate...
    esmSurf->TransferPointOntoSurface(mortarElem->transVect, mortarElem->rotationCenter, mortarElem->rotationAngle, globIntPoint);
  }
  // -------kirill (periodic boundary conditions)-------------------------------------------------------------------------------------

  // ATTENTION! Global2Local on the volume element will cause integration points outside the element when the interface is not coplanar!
  esmVol->Global2Local(lpVol.coord, globIntPoint);
  // create new local point mapped for volume element
  this->lpmVol.reset(new LocPointMapped());
  esmVol->CalcNormalOutOfVolume(locNormal,lpVol,volElem,surfElem);
  // assign integration points to the (new) volume element
  this->lpmVol->Set(lpVol, esmVol, this->weight);
  // calculate global normal pointing out of new volume element
  this->normal = Transpose(this->lpmVol->jacInv) * locNormal;
  this->normal /= this->normal.NormL2();
  if (!usePrimary)
    this->normal *= -1.0;

#ifndef NDEBUG
  // sanity check if the surf->vol mapping is correct
  Vector<Double> globVolIntPoint;
  esmVol->Local2Global(globVolIntPoint, lpVol);
  LOG_DBG2(locPointMapped) << "------------------------------------------------------------------------------------------------------------" << std::endl;
  LOG_DBG2(locPointMapped) "globCoordVol: " << globVolIntPoint << std::endl;
  LOG_DBG2(locPointMapped) "locCoordVol: " << lpVol.coord << std::endl << std::endl;
  LOG_DBG2(locPointMapped) "globCoordMortarElem: " << globIntPoint << std::endl;
  LOG_DBG2(locPointMapped) "locCoordMortarElem: " << this->lp.coord << std::endl;
  LOG_DBG2(locPointMapped) "usePrimary:" << usePrimary << std::endl;
  LOG_DBG2(locPointMapped) "surfElemNums---:" << mortarElem->ptPrimary->elemNum << ", " << mortarElem->ptSecondary->elemNum << std::endl << std::endl;
  LOG_DBG2(locPointMapped) "volElemNums---:" << mortarElem->ptPrimary->ptVolElems[0]->elemNum << ", " << mortarElem->ptSecondary->ptVolElems[0]->elemNum << std::endl << std::endl;
  LOG_DBG2(locPointMapped) "normal loc-----------:" << locNormal << std::endl;
  LOG_DBG2(locPointMapped) "normal glob-------------------:" << this->normal << std::endl << std::endl;
  LOG_DBG2(locPointMapped) << "------------------------------------------------------------------------------------------------------------" << std::endl;
  for (UInt iDim = 0; iDim < globIntPoint.GetSize(); ++iDim)
  {
    assert(abs(globIntPoint[iDim] - globVolIntPoint[iDim]) < EPS);
    assert(esmVol->CoordIsInsideElem(lpVol.coord, NORM_EPS*1e3));
  }
#endif
}

Vector<double>& LocPointMapped::GetGlobal(Vector<double>& coord, const LocPoint* loc, bool fallback, bool update) const
{
  loc = loc != NULL ? loc : &(this->lp); // use own lp as default

  if(shapeMap)
    shapeMap->Local2Global(coord,*loc);
  else
  {
    assert(fallback);
    assert(loc->number != LocPoint::NOT_SET);
    assert(loc->number > 0);
    domain->GetGrid()->GetNodeCoordinate(coord, (unsigned int) loc->number, update);
  }
  return coord;
}

// ========================================================================
//  ElemShapeMap
// ========================================================================

// Definition of finite element type mappings
static EnumTuple elemShapeTuples[] = { EnumTuple(ElemShapeMap::NO_TYPE,
    "NO_TYPE"), EnumTuple(ElemShapeMap::LAGRANGE, "LAGRANGE"), EnumTuple(
    ElemShapeMap::LAGRANGE_BLENDED, "LAGRANGE_BLENDED"), EnumTuple(
    ElemShapeMap::ANALYTICAL, "ANALYTICAL") };

Enum<ElemShapeMap::ShapeMapType> ElemShapeMap::shapeMapType = Enum<
    ElemShapeMap::ShapeMapType>("Finite Element Shape Mapping Types",
    sizeof(elemShapeTuples) / sizeof(EnumTuple), elemShapeTuples);

ElemShapeMap::ElemShapeMap(Grid* ptGrid) {
  ptGrid_ = ptGrid;
}

ElemShapeMap::~ElemShapeMap() {

}

std::string ElemShapeMap::ToString() const
{
  std::stringstream ss;
  ss << "e=" << ptElem_->elemNum << " c=" << ptElem_->connect.ToString() << " bc=" << ptElem_->extended->barycenter.ToString();
  return ss.str();
}

void ElemShapeMap::SetElem(const Elem* ptElem, bool isUpdated) {
  ptElem_ = ptElem;
  isUpdated_ = isUpdated;
  isAxi_ = ptGrid_->IsAxi();
  SetModelDepth(ptGrid_->GetDepth2dPlane());

  // Check, if the element is a surface element
  if (Elem::shapes[ptElem->type].dim == ptGrid_->GetDim() - 1) {
    try {
      ptSurfElem_ = dynamic_cast<const SurfElem*>(ptElem);
    } catch (...) {
      EXCEPTION("Could not convert element #" << ptElem->elemNum << " to a surface element");
    }
  } else {
    ptSurfElem_ = nullptr;
  }
}

// ========================================================================
//  Lagrangian Element Shape Map
// ========================================================================

LagrangeElemShapeMap::LagrangeMapSingleton::LagrangeMapSingleton() {
  //obtain thread local copy this method will only be called once in the program!
  for(UInt aT = 0; aT<feMap_.GetNumSlots();++aT){
    std::map<Elem::FEType, FeH1LagrangeExpl* >& tMap = feMap_.Mine(aT);
    tMap[Elem::ET_LINE2]   = new FeH1LagrangeLine1();
    tMap[Elem::ET_LINE3]   = new FeH1LagrangeLine2();
    tMap[Elem::ET_TRIA3]   = new FeH1LagrangeTria1();
    tMap[Elem::ET_TRIA6]   = new FeH1LagrangeTria2();
    tMap[Elem::ET_QUAD4]   = new FeH1LagrangeQuad1();
    tMap[Elem::ET_QUAD8]   = new FeH1LagrangeQuad2();
    tMap[Elem::ET_QUAD9]   = new FeH1LagrangeQuad9();
    tMap[Elem::ET_TET4]    = new FeH1LagrangeTet1();
    tMap[Elem::ET_TET10]   = new FeH1LagrangeTet2();
    tMap[Elem::ET_HEXA8]   = new FeH1LagrangeHex1();
    tMap[Elem::ET_HEXA20]  = new FeH1LagrangeHex2();
    tMap[Elem::ET_HEXA27]  = new FeH1LagrangeHex27();
    tMap[Elem::ET_WEDGE6]  = new FeH1LagrangeWedge1();
    tMap[Elem::ET_WEDGE15] = new FeH1LagrangeWedge2();
    tMap[Elem::ET_WEDGE18] = new FeH1LagrangeWedge18();
    tMap[Elem::ET_PYRA5]   = new FeH1LagrangePyra1();
    tMap[Elem::ET_PYRA13]  = new FeH1LagrangePyra2();
    tMap[Elem::ET_PYRA14]  = new FeH1LagrangePyra14();
  }
}

LagrangeElemShapeMap::LagrangeMapSingleton::LagrangeMapSingleton(const LagrangeMapSingleton&) {
  return;
}            

LagrangeElemShapeMap::LagrangeMapSingleton& LagrangeElemShapeMap::LagrangeMapSingleton::operator=(const LagrangeMapSingleton&) {
  return *this;
}

LagrangeElemShapeMap::LagrangeMapSingleton::~LagrangeMapSingleton() {
  // delete reference elements
  feMap_.Clear();
//  for(UInt aT = 0; aT<feMap_.GetNumSlots();++aT){
//    std::map<Elem::FEType, FeH1LagrangeExpl* >& tMap = feMap_.Mine(aT);
//    std::map<Elem::FEType, FeH1LagrangeExpl*>::iterator it = tMap.begin();
//    for (; it != tMap.end(); ++it) {
//      delete it->second;
//    }
//    tMap.clear();
//  }
}

LagrangeElemShapeMap::LagrangeMapSingleton&
LagrangeElemShapeMap::LagrangeMapSingleton::getInstance() {
  static LagrangeMapSingleton instance;
  return instance;
}

LagrangeElemShapeMap::LagrangeElemShapeMap(Grid* ptGrid) :
    ElemShapeMap(ptGrid), elems_(LagrangeMapSingleton::getInstance()) {
  type_ = LAGRANGE;
  intScheme_ = ptGrid_->GetIntegrationScheme();
  ptFe_ = nullptr;
  shape_ = NULL;
}

LagrangeElemShapeMap::~LagrangeElemShapeMap() {

//  // Remove pointers to all reference elements
//  std::map<Elem::FEType, FeH1LagrangeExpl *>::iterator it = feMap_.begin();
//  for( ; it != feMap_.end(); it++ ) {
//    delete  it->second;
//  }
//  feMap_.clear();
}

void LagrangeElemShapeMap::Local2Global(Vector<Double>& globPoint,
    const LocPoint& lp) {

  //step 1: evaluate shape fncs. at local coordinate
  ptFe_->GetShFnc(shFnc_, lp, nullptr, 0);

  // step2: multiply shape fncs for each dimension with according matrix entries
  globPoint = coords_ * shFnc_;
}

void LagrangeElemShapeMap::Global2Local(Vector<Double>& locPoint,
    const Vector<Double>& globalPoint) {

  // first of all check if the coordinate coincides with one node
  bool isNode = Global2LocalOnNode(locPoint, globalPoint);
  if (isNode) {
    return;
  } else {
    locPoint.Init();
  }

  Elem::FEType curType = ptFe_->FeType();
  switch (curType) {
  case Elem::ET_LINE3:
    //WARN("Using Global2LocalLine2 for ET_LINE3!");
  case Elem::ET_LINE2:
    Global2LocalLine2(locPoint, globalPoint);
    break;
  case Elem::ET_QUAD4:
  case Elem::ET_QUAD8:
  case Elem::ET_QUAD9:
    // Use specialization of Global2Local for quadrilaterals
    // that works only in 2D!
    if (ptGrid_->GetDim() == 2) {
      Global2LocalQuad4(locPoint, globalPoint);
    } else {
      Global2LocalDuester(locPoint,globalPoint);
    }
    break;
  case Elem::ET_TRIA6:
  case Elem::ET_HEXA20:
  case Elem::ET_HEXA27:
  case Elem::ET_HEXA8:
  case Elem::ET_PYRA5:
  case Elem::ET_TET10:
  case Elem::ET_PYRA13:
  case Elem::ET_WEDGE6:
  case Elem::ET_WEDGE15:
  case Elem::ET_WEDGE18:
    Global2LocalDuester(locPoint, globalPoint);
    break;
  case Elem::ET_TRIA3:
  case Elem::ET_TET4:
    Global2LocalBarycentric(locPoint, globalPoint);
    break;
  default:
    Global2LocalGeneral(locPoint, globalPoint);
    break;
  }
}

void LagrangeElemShapeMap::TranslatePointOntoSurface(const Vector<Double>& direction, Vector<Double>& point)
{
  const ElemShape& shape = *shape_;

  // make sure that we are not trying to project the point onto a volume element
  if (shape.dim == ptGrid_->GetDim())
  {
    EXCEPTION("It is possible to project onto a surface element only.");
    return;
  }

  if (ptGrid_->GetDim() == 2)
  {
    // TODO: think how to do in 2D
    // projection onto a segment
    Vector<Double> p1, p2, norm(2);
    Double lambda;
    Local2Global(p1, shape.nodeCoords[0]);
    Local2Global(p2, shape.nodeCoords[1]);
    p2 -= p1;

    norm[0] = p2[1];
    norm[1] = -p2[0];
    norm.Normalize();
    p2 = point - p1;
    lambda = norm.Inner(p2)/norm.Inner(direction);
    point -= direction*lambda;
  }
  else
  {
    // projection onto a polygon
    Vector<Double> p0, p1, p2, norm;
    Double lambda;
    Local2Global(p0, shape.nodeCoords[0]);
    Local2Global(p1, shape.nodeCoords[1]);
    Local2Global(p2, shape.nodeCoords[2]);

    p1 -= p0;
    p2 -= p0;
    p1.CrossProduct(p2, norm);
    norm.Normalize();
    p1 = point - p0;
    lambda = norm.Inner(p1)/norm.Inner(direction);
    point -= direction*lambda;
  }
}

void LagrangeElemShapeMap::TransferPointOntoSurface(const Vector<Double>& translationVec, 
                                                    const Vector<Double>& rotationCenter, 
                                                    const Double& rotationAngle,
                                                    Vector<Double>& point) {
  const ElemShape& shape = *shape_;
  // disable for 3D case
  if (ptGrid_->GetDim() != 2)
    EXCEPTION("Cake-piece projection only implemented in 2D so far!");

  // make sure that we are not trying to project the point onto a volume element
  if (shape.dim == ptGrid_->GetDim())
    EXCEPTION("It is possible to project onto a surface element only.");

  Vector<Double> p1, p2, p3;
  // get coordinates of 2 arbitrary nodes on current surface
  Local2Global(p1, shape.nodeCoords[0]);
  Local2Global(p2, shape.nodeCoords[1]);
  // construct tangential direction vector of currenct surface in p2
  p2 -= p1;
  p2.Normalize();
  // check wether the integration point is on the current surface
  // in this case, no transfer is necessary
  p3 = point - p1;
  p3.Normalize();
  Double cross = p2[0]*p3[1] - p2[1]*p3[0];
  if (fabs(cross) < NORM_EPS*1e3)
    return;

  // check for the direction of translation to decide in which direction we need to transfer
  if (translationVec.Inner(p3) < 0) {
    // 1. translate first
    point += translationVec;
    // 1. translate point to center of rotation
    point -= rotationCenter;
    // 2. rotate
    Vector<Double> tmp(2);
    tmp[0] = point[0]*cos(-rotationAngle) - point[1]*sin(-rotationAngle);
    tmp[1] = point[0]*sin(-rotationAngle) + point[1]*cos(-rotationAngle);
    point = tmp;
    // 3. translate back and onto the other surface
    point += rotationCenter;
  } else {
    // 1. translate point to center of rotation
    point -= rotationCenter;
    // 2. rotate
    Vector<Double> tmp(2);
    tmp[0] = point[0]*cos(rotationAngle) - point[1]*sin(rotationAngle);
    tmp[1] = point[0]*sin(rotationAngle) + point[1]*cos(rotationAngle);
    point = tmp;
    // 3. translate back and onto the other surface
    point += rotationCenter - translationVec;
  }

  LOG_DBG2(lagrangeElemShapeMap) << "------------------------------------------------------------------------------------------------------------" << std::endl;
  LOG_DBG2(lagrangeElemShapeMap) << "p2: " <<  p2 << std::endl;
  LOG_DBG2(lagrangeElemShapeMap) << "p3: " <<  p3 << std::endl;
  LOG_DBG2(lagrangeElemShapeMap) << "cross: " <<  cross << std::endl;
  LOG_DBG2(lagrangeElemShapeMap) << "translationVec: " << translationVec << std::endl;
  LOG_DBG2(lagrangeElemShapeMap) << "rotationCenter: " << rotationCenter << std::endl;
  LOG_DBG2(lagrangeElemShapeMap) << "rotationAngle: " << rotationAngle << std::endl;
  LOG_DBG2(lagrangeElemShapeMap) << "point: " << point << std::endl;
  LOG_DBG2(lagrangeElemShapeMap) << "point-final: " << point << std::endl;
  LOG_DBG2(lagrangeElemShapeMap) << "------------------------------------------------------------------------------------------------------------" << std::endl;
}

void LagrangeElemShapeMap::Global2LocalBarycentric(Vector<Double>& locPoint,
    const Vector<Double>& globalPoint) {

  const ElemShape & shape = *shape_;
  
  if (ptFe_->FeType() == Elem::ET_TRIA3) {

    Double lamb1, lamb2;

    if (ptGrid_->GetDim() == 3)
    {
      // kirill:
      // form here "http://gamedev.stackexchange.com/questions/23743/whats-the-most-efficient-way-to-find-barycentric-coordinates"
      // in 3D the global point has to be projected to the element plane
      Vector<Double> p0, p1, p2, norm;
      Local2Global(p0, shape.nodeCoords[0]);
      Local2Global(p1, shape.nodeCoords[1]);
      Local2Global(p2, shape.nodeCoords[2]);

      p1 -= p0; // vector AB
      p2 -= p0; // vector AC
      p0 -= globalPoint; // vector PA
      p1.CrossProduct(p2, norm);
      // we calculate the doubled square of the initial triangle ABC
      Double sqrABC_2 = norm.NormL2();
      // and the doubled squares of the triangles PBC and PAC, respectively
      p2.CrossProduct(p0, norm);
      Double sqrPAC_2 = norm.NormL2();
      p1.CrossProduct(p0, norm);
      Double sqrPBC_2 = sqrABC_2 - sqrPAC_2 - norm.NormL2();

      lamb1 = sqrPBC_2/sqrABC_2;
      lamb2 = sqrPAC_2/sqrABC_2;
    }
    else
    {
      lamb1 = (coords_[1][1] - coords_[1][2])
            * (globalPoint[0] - coords_[0][2])
            + (coords_[0][2] - coords_[0][1]) * (globalPoint[1] - coords_[1][2]);
      lamb1 /= ((coords_[1][1] - coords_[1][2]) * (coords_[0][0] - coords_[0][2])
            + (coords_[0][2] - coords_[0][1]) * (coords_[1][0] - coords_[1][2]));

      lamb2 = (coords_[1][2] - coords_[1][0])
            * (globalPoint[0] - coords_[0][2])
            + (coords_[0][0] - coords_[0][2]) * (globalPoint[1] - coords_[1][2]);
      lamb2 /= ((coords_[1][1] - coords_[1][2]) * (coords_[0][0] - coords_[0][2])
            + (coords_[0][2] - coords_[0][1]) * (coords_[1][0] - coords_[1][2]));
    }

    Double lamb3 = 1 - lamb1 - lamb2;
    locPoint.Resize(2);

    locPoint[0] = lamb1 * shape.nodeCoords[0][0]
        + lamb2 * shape.nodeCoords[1][0] + lamb3 * shape.nodeCoords[2][0];
    locPoint[1] = lamb1 * shape.nodeCoords[0][1]
        + lamb2 * shape.nodeCoords[1][1] + lamb3 * shape.nodeCoords[2][1];
  } else if (ptFe_->FeType() == Elem::ET_TET4) {
    Matrix<Double> baryMat(3, 3);
    Matrix<Double> baryMatInv(3, 3);
    Vector<Double> b(3);
    Vector<Double> lamb(3);
    locPoint.Resize(3);

    baryMat.Init();
    baryMat[0][0] = coords_[0][0] - coords_[0][3];
    baryMat[0][1] = coords_[0][1] - coords_[0][3];
    baryMat[0][2] = coords_[0][2] - coords_[0][3];
    baryMat[1][0] = coords_[1][0] - coords_[1][3];
    baryMat[1][1] = coords_[1][1] - coords_[1][3];
    baryMat[1][2] = coords_[1][2] - coords_[1][3];
    baryMat[2][0] = coords_[2][0] - coords_[2][3];
    baryMat[2][1] = coords_[2][1] - coords_[2][3];
    baryMat[2][2] = coords_[2][2] - coords_[2][3];

    baryMat.Invert(baryMatInv);

    b[0] = globalPoint[0] - coords_[0][3];
    b[1] = globalPoint[1] - coords_[1][3];
    b[2] = globalPoint[2] - coords_[2][3];

    lamb = baryMatInv * b;
    Double lamb4 = 1 - lamb[0] - lamb[1] - lamb[2];

    locPoint[0] = lamb[0] * shape.nodeCoords[0][0]
        + lamb[1] * shape.nodeCoords[1][0] + lamb[2] * shape.nodeCoords[2][0]
        + lamb4 * shape.nodeCoords[3][0];
    locPoint[1] = lamb[0] * shape.nodeCoords[0][1]
        + lamb[1] * shape.nodeCoords[1][1] + lamb[2] * shape.nodeCoords[2][1]
        + lamb4 * shape.nodeCoords[3][1];
    locPoint[2] = lamb[0] * shape.nodeCoords[0][2]
        + lamb[1] * shape.nodeCoords[1][2] + lamb[2] * shape.nodeCoords[2][2]
        + lamb4 * shape.nodeCoords[3][2];
  }
}

bool LagrangeElemShapeMap::Global2LocalOnNode(Vector<Double>& locPoint,
    const Vector<Double>& glob) {

  const ElemShape & shape = *shape_;
  
  //first get the local node coordinates of the given element
  Vector<Double> curNodeGlobCoord;
  Vector<Double> diffVec;
  LocPoint curNodePoint;
  curNodePoint.coord.Resize(coords_.GetNumRows(), 0.0);
  locPoint.Resize(coords_.GetNumRows());
  bool retVal = false;
  for (UInt i = 0; i < coords_.GetNumCols(); i++) {
    for (UInt d = 0; d < shape.dim; ++d) {
      curNodePoint.coord[d] = shape.nodeCoords[i][d];
    }
    Local2Global(curNodeGlobCoord, curNodePoint);
    diffVec = glob - curNodeGlobCoord;
    if (diffVec.NormL2() < 1e-10) {
      locPoint = curNodePoint.coord;
      retVal = true;
      break;
    }
  }
  return retVal;
}

void LagrangeElemShapeMap::Global2LocalQuad4(Vector<Double>& locPoint,
    const Vector<Double>& globalPoint) {
  //based on FEM skript found by simon at
  //http://www.colorado.edu/engineering/cas/courses.d/IFEM.d/IFEM.Ch23.d/IFEM.Ch23.pdf

  //corner coords
  Double x1 = 0.0, x2 = 0.0, x3 = 0.0, x4 = 0.0;
  Double y1 = 0.0, y2 = 0.0, y3 = 0.0, y4 = 0.0;
  //global point
  Double x_p, y_p;
  Double x_b, y_b;
  Double x_cx, y_cx, x_ce, y_ce, A, J_1, J_2, x_0, y_0, x_p0, y_p0;
  Double b_xi, b_eta, c_xi, c_eta;

  locPoint.Resize(2);
  locPoint.Init();

  x1 = coords_[0][0];
  x2 = coords_[0][1];
  x3 = coords_[0][2];
  x4 = coords_[0][3];

  y1 = coords_[1][0];
  y2 = coords_[1][1];
  y3 = coords_[1][2];
  y4 = coords_[1][3];

  x_p = globalPoint[0];
  y_p = globalPoint[1];

  //begin algorithm
  x_b = x1 - x2 + x3 - x4;
  y_b = y1 - y2 + y3 - y4;

  x_cx = x1 + x2 - x3 - x4;
  y_cx = y1 + y2 - y3 - y4;
  x_ce = x1 - x2 - x3 + x4;
  y_ce = y1 - y2 - y3 + y4;

  A = 0.5 * ((x3 - x1) * (y4 - y2) - (x4 - x2) * (y3 - y1));
  J_1 = (x3 - x4) * (y1 - y2) - (x1 - x2) * (y3 - y4);
  J_2 = (x2 - x3) * (y1 - y4) - (x1 - x4) * (y2 - y3);

  x_0 = 0.25 * (x1 + x2 + x3 + x4);
  y_0 = 0.25 * (y1 + y2 + y3 + y4);

  x_p0 = x_p - x_0;
  y_p0 = y_p - y_0;

  b_xi = A - (x_p0 * y_b) + (y_p0 * x_b);
  b_eta = (-1.0 * A) - (x_p0 * y_b) + (y_p0 * x_b);

  c_xi = (x_p0 * y_cx) - (y_p0 * x_cx);
  c_eta = (x_p0 * y_ce) - (y_p0 * x_ce);

  locPoint[0] = 2.0 * c_xi
      / ((-1.0 * std::sqrt(b_xi * b_xi - 2.0 * J_1 * c_xi)) - b_xi);
  locPoint[1] = 2.0 * c_eta
      / (std::sqrt(b_eta * b_eta + (2.0 * J_2 * c_eta)) - b_eta);
  /*
   //BaseFE::Global2LocalCoords(localCoords,globalCoords,coordMat);
   //return;
   //Version 0 algorithm from duester script
   Vector<Double> delta_xi; // update for Newton-Raphson method
   Vector<Double> xi_start; // local start point for Newton-Raphson method
   Vector<Double> xi_k; // local point at iteration k
   Vector<Double> f; //current right hand side
   Vector<Double> f_start; //current right hand side
   Double f_old; //storing the absolute value of search direction
   Double f_test; //storing the absolute value of search direction
   UInt k = 0; //iteration counter
   Integer l = 0; //stepping value
   Double jacDet = 0;
   Matrix<Double> J; // Jacobian at local point xi_k
   locPoint.Resize(2);
   locPoint.Init();

   Double tolerance = 1e-10;

   //initialize everything
   xi_start.Resize(2);
   delta_xi.Resize(2);
   xi_k.Resize(2);
   J.Resize(2, 2);

   //first initial guess is zero
   f.Resize(2);

   // Perform Newton-Raphson method on the list of global points
   // to find local coordinates within this element.

   f.Init();
   J.Init();
   xi_start.Init();
   xi_k.Init();
   k = 0;
   f_old = 222e20;
   f_test = 0;

   Local2Global(f, xi_k);

   f = f - globalPoint;
   f_old = f.NormL2();
   f_start = f;
   for(Double w=0.5; w <= 1.0; w+=0.5){
   for(UInt j=1;j<3;j++){
   for(UInt m=1;m<3;m++){
   xi_k[0] = pow(-1,j)*w;
   xi_k[1] = pow(-1,m)*w;
   Local2Global(f, xi_k);
   f = f - globalPoint;
   f_test = f.NormL2();
   if(f_old > f_test){
   xi_start  = xi_k;
   f_start = f;
   f_old = f_test;
   }
   }
   }
   }
   xi_k = xi_start;
   f = f_start;
   while(f_test > tolerance && k < 20){
   delta_xi.Init();
   xi_start.Init();
   // Calculate Jacobian at iteration point xi_k
   this->CalcJ(J,xi_k);
   jacDet = + J[0][0]*J[1][1] - J[1][0]*J[0][1];
   //  ( f_0  J_01 )
   //  ( f_1  J_11 )
   delta_xi[0] =  (-J[1][1]*f[0] + J[0][1]*f[1])/jacDet;

   //  ( J_00  f_0 )
   //  ( J_10  f_1 )
   delta_xi[1] =  (-J[0][0]*f[1] + J[1][0]*f[0])/jacDet;
   //xi_start = xi_k + (delta_xi * 2.0);
   //Local2GlobalCoord(f, xi_start, coordMat, NULL);
   //f = f - globalPoint;
   //f_test = f.NormL2();
   l = 1;
   //perform damping
   while(l<30 && f_test >= f_old){
   Double dampFac = 1.0/std::pow(2.0,l-1.0);
   xi_start = xi_k + (delta_xi * dampFac);
   Local2Global(f, xi_start);
   f = f - globalPoint;
   f_test = f.NormL2();
   l++;
   }
   f_old = f_test;
   //if(l >= 30){
   //  std::cout << "Daming iteration was not convergend" << std::endl;
   //}
   xi_k = xi_start;
   k++;
   }
   //if( f_test > tolerance){
   //  std::cout << "performed " << k << " iterations to reach the point" << std::endl<< xi_k << std::endl;
   //  //ONLY for DEBUGGING
   //  Local2Global(f, xi_k);
   //  std::cout << "Calculated global point :" << std::endl << f << std::endl;
   //  std::cout << "Original global coord " << std::endl << globalPoint << std::endl;
   //  std::cout << "The error was: " << f_test << std::endl<< std::endl;
   //  //std::cout << std::endl;
   //}

   // Put local coordinate of point into matrix.
   for(UInt j = 0; j < 2; j++) {
   if(xi_k.GetSize() == 0){
   std::cerr << "local2global messed up setting everything to zero" << std::endl;
   std::cerr << globalPoint << std::endl;
   xi_k.Resize(2);
   xi_k.Init();
   }
   locPoint[j] = xi_k[j];
   }*/
}

Double LagrangeElemShapeMap::GetLocDirJac( Vector<Double>& delta_xi, 
                                           const Vector<Double>& f,
                                           const Matrix<Double>& J) {
  
  // note: in NACS jacDet is computed via CalcJDet function without setting the flag useDepth, so using the default
  // (false?); 
  // here the determinant is computed explicitly (currently); apply no depth scaling here
  // BUT: we cannot use the CalcJDet function here as in NACS as this function computes J and here J is passed
  // in NACS xi is passed instead of J
  Double jacDet = 0.0;
  UInt globDim = J.GetNumRows(); // dimension of the grid
  UInt locDim = J.GetNumCols();  // dimension of the element
  
  // ===============================================================
  //  2 DIMENSIONAL MESH
  // ===============================================================
  if (globDim == 2) {
    if( locDim == 2) {
      // Find new local search direction for 2D -> 2D mapping using
      // Cramer's rule.

      // Jacobian Matrix:
      //
      //  J = ( J_00  J_01 )
      //      ( J_10  J_11 )
      jacDet = +J[0][0] * J[1][1] - J[1][0] * J[0][1];

      //  ( f_0  J_01 )
      //  ( f_1  J_11 )
      delta_xi[0] = -J[1][1] * f[0] + J[0][1] * f[1];

      //  ( J_00  f_0 )
      //  ( J_10  f_1 )
      delta_xi[1] = -J[0][0] * f[1] + J[1][0] * f[0];
    } else if( locDim == 1 ) {
      EXCEPTION("Not implemented yet")
    }

  } else {
    // ===============================================================
    //  3 DIMENSIONAL MESH
    // ===============================================================
    if (locDim == 3) {
      // Find new local search direction for 3D -> 3D mapping using
       // Cramer's rule.

       //      ( J_00  J_01  J_02 )
       //  J = ( J_10  J_11  J_12 )
       //      ( J_20  J_21  J_22 )

       jacDet = +J[0][0] * J[1][1] * J[2][2] - J[0][0] * J[2][1] * J[1][2]
                - J[1][0] * J[0][1] * J[2][2] + J[2][1] * J[1][0] * J[0][2]
                - J[1][1] * J[2][0] * J[0][2] + J[2][0] * J[0][1] * J[1][2];

       //  ( f_0  J_01  J_02 )
       //  ( f_1  J_11  J_12 )
       //  ( f_2  J_21  J_22 )
       delta_xi[0] = -J[1][1] * J[2][2] * f[0] + J[2][1] * J[1][2] * f[0]
                     - J[2][1] * J[0][2] * f[1] + J[0][1] * J[2][2] * f[1]
                     - J[0][1] * J[1][2] * f[2] + J[1][1] * J[0][2] * f[2];

       //  ( J_00  f_0  J_02 )
       //  ( J_10  f_1  J_12 )
       //  ( J_20  f_2  J_22 )
       delta_xi[1] = -J[1][2] * J[2][0] * f[0] + J[1][0] * J[2][2] * f[0]
                     - J[0][0] * J[2][2] * f[1] + J[2][0] * J[0][2] * f[1]
                     - J[1][0] * J[0][2] * f[2] + J[1][2] * J[0][0] * f[2];

       //  ( J_00  J_01  f_0 )
       //  ( J_10  J_11  f_1 )
       //  ( J_20  J_21  f_2 )
       delta_xi[2] = -J[1][0] * J[2][1] * f[0] + J[2][0] * J[1][1] * f[0]
                    - J[0][1] * J[2][0] * f[1] + J[2][1] * J[0][0] * f[1]
                    - J[0][0] * J[1][1] * f[2] + J[1][0] * J[0][1] * f[2];
    } else if( locDim == 2 ) {
      // Find new local search direction for 2D -> 3D mapping.
      // Project 3D difference vector onto 2D basis given by
      // the Jacobian to find the new 2D search direction.
      // we use the cross product to add a third orthogonal base vector for 
      // the Jacobi matrix: e3  = J1 x J2
      // this should mimic a 3D element of unit depth in the 3rd dimension
      // Jacobian Matrix:
      //
      //      ( J_00  J_01  e3_1 )
      //  J = ( J_10  J_11  e3_2 )
      //      ( J_20  J_21  e3_2 )
      Vector<Double> e3(3);
      e3[0] = J[1][0] * J[2][1] - J[2][0] * J[1][1];
      e3[1] = J[2][0] * J[0][1] - J[0][0] * J[2][1];
      e3[2] = J[0][0] * J[1][1] - J[1][0] * J[0][1];
      Double e3norm;
      e3norm = e3.NormL2();
      e3 /= e3norm;

      // determinant
      jacDet = +J[0][0] * J[1][1] * e3[2] - J[0][0] * J[2][1] * e3[1]
               -J[1][0] * J[0][1] * e3[2] + J[2][1] * J[1][0] * e3[0]
               -J[1][1] * J[2][0] * e3[0] + J[2][0] * J[0][1] * e3[1];

      // first local offset
      //  ( f_0  J_01  J_02 )
      //  ( f_1  J_11  J_12 )
      //  ( f_2  J_21  J_22 )
      delta_xi[0] = -J[1][1] * e3[2] * f[0] + J[2][1] * e3[1] * f[0]
                    -J[2][1] * e3[0] * f[1] + J[0][1] * e3[2] * f[1]
                    -J[0][1] * e3[1] * f[2] + J[1][1] * e3[0] * f[2];

      // second local offset
      //  ( J_00  f_0  J_02 )
      //  ( J_10  f_1  J_12 )
      //  ( J_20  f_2  J_22 )
      delta_xi[1] = -e3[1]   * J[2][0] * f[0] + J[1][0] * e3[2]   * f[0]
                    -J[0][0] * e3[2]   * f[1] + J[2][0] * e3[0]   * f[1]
                    -J[1][0] * e3[0]   * f[2] + e3[1]   * J[0][0] * f[2];

    } else {
      // Find new local search direction for 1D -> 3D mapping.
      EXCEPTION("Not implemented yet")

    }
  }
  
  // Normalize the local search direction
  delta_xi /= jacDet;  
  return jacDet;
}

void LagrangeElemShapeMap::Global2LocalDuester(
    Vector<Double> &locPoint, const Vector<Double> &globalPoint) {
  // Perform Newton-Raphson method on the list of global points to find local
  // coordinates within this element. Version 0 algorithm from Düster script
  Vector<Double> delta_xi;     // update for Newton-Raphson method
  Vector<Double> xi_start;     // local start point for Newton-Raphson method
  Vector<Double> f;            // current right hand side
  Vector<Double> f_start;      // current right hand side
  Double f_old = 0;            // storing the absolute value of search direction
  Double f_test = 0;           // storing the absolute value of search direction
  UInt k = 0;                  // iteration counter (outer loop)
  Integer l = 0;               // exponent for step size (inner loop)
  Double tolerance = EPS;      // tolerance for convergence
  Matrix<Double> J;            // Jacobian at local point xi_k
  UInt elemDim = shape_->dim;  // dimension of the element

  locPoint.Resize(elemDim); // we use locPoint directly to iteratively find the
                            // local coordinates (xi_k)
  locPoint.Init();
  xi_start.Resize(elemDim);
  xi_start.Init();
  delta_xi.Resize(elemDim);
  J.Resize(elemDim, elemDim);
  J.Init();
  f.Resize(elemDim);
  f.Init();

  // get global coordinates to initialize f
  Local2Global(f, locPoint);
  f = f - globalPoint;
  f_old = f.NormL2();
  f_start = f; // store the initial value of f
  // try to find better initial f values
  // loop through possible starting positions: [-0.5, 0.5, -1.0, 1.0] and keep
  // the best one
  if (elemDim == 1) {
    for (Double w = 0.5; w <= 1.0; w += 0.5) {
      for (UInt j = 1; j < 3; j++) {
        locPoint[0] = pow(-1, j) * w;
        Local2Global(f, locPoint);
        f = f - globalPoint;
        f_test = f.NormL2();
        // if the new f value is smaller than the old one, keep the new one and
        // set the new start positions
        if (f_old > f_test) {
          xi_start = locPoint;
          f_start = f;
          f_old = f_test;
        }
      }
    }
  } else if (elemDim == 2) {
    for (Double w = 0.5; w <= 1.0; w += 0.5) {
      for (UInt j = 1; j < 3; j++) {
        for (UInt m = 1; m < 3; m++) {
          locPoint[0] = pow(-1, j) * w;
          locPoint[1] = pow(-1, m) * w;
          Local2Global(f, locPoint);
          f = f - globalPoint;
          f_test = f.NormL2();
          if (f_old > f_test) {
            xi_start = locPoint;
            f_start = f;
            f_old = f_test;
          }
        }
      }
    }
  } else if (elemDim == 3) {
    for (Double w = 0.5; w <= 1.0; w += 0.5) {
      for (UInt j = 1; j < 3; j++) {
        for (UInt m = 1; m < 3; m++) {
          for (UInt n = 1; n < 3; n++) {
            locPoint[0] = pow(-1, j) * w;
            locPoint[1] = pow(-1, m) * w;
            locPoint[2] = pow(-1, n) * w;
            Local2Global(f, locPoint);
            f = f - globalPoint;
            f_test = f.NormL2();
            if (f_old > f_test) {
              xi_start = locPoint;
              f_start = f;
              f_old = f_test;
            }
          }
        }
      }
    }
  }
  // assnign the determined start values
  locPoint = xi_start;
  f = f_start;
  while (f_test > tolerance && k < 60) {
    delta_xi.Init();
    xi_start.Init();
    // Calculate Jacobian (J) at iteration point xi_k (locPoint)
    this->CalcJ(J, locPoint);
    // Calculate local search direction delta_xi
    GetLocDirJac(delta_xi, f, J);

    l = 0;
    // perform damping to obtain coordinate iteratively
    while (l < 60 && f_test >= f_old) {
      Double dampFac = 1.0 / std::pow(2.0, (Double)l);
      xi_start = locPoint + (delta_xi * dampFac);
      Local2Global(f, xi_start);
      f = f - globalPoint;
      f_test = f.NormL2();
      l++;
    }
    f_old = f_test;
    locPoint = xi_start;
    k++;
    // raise exception if the algorithm does not converge
    if (k == 60) {
      WARN("Newton-Raphson algorithm did not converge for global to local "
           "conversion in point"
           << globalPoint << ". The determined local coordinates are "
           << locPoint << ". The remaining f_test is " << f_test << ".");
    }
  }
}

void LagrangeElemShapeMap::Global2LocalLine2(Vector<Double>& locPoint,
    const Vector<Double>& globalPoint) {
  Vector<Double> c0, c1; // endpoint-coordinates
  Vector<Double> xi_k(1); // local point at iteration k
  Vector<Double> diffVecToPoint, diffVecWholeElem;
  Double s;
  Double lengthOfWholeElem, lengthToPoint, fac;
  UInt globDim = globalPoint.GetSize();
  const ElemShape & shape = *shape_;

  // x-----------o-------------x
  // c0          point         coordMat[i][1]
  //
  //  ------------------------>
  //  diffVecWholeElem
  //
  //  ---------->
  //  diffVecToPoint
  //
  //  lengthOfWholeElem = length(diffVecWholeElem)
  //  lengthToPoint = length(diffVecToPoint)
  //  fac = 1 / lengthOfWholeElem
  //  s = lengthToPoint / lengthWholeElem

  // Get coordinates of the endpoints
  c0.Resize(globDim);
  diffVecToPoint.Resize(globDim);
  diffVecWholeElem.Resize(globDim);

  xi_k[0] = shape.nodeCoords[0][0];
  Local2Global(c0, xi_k);

  xi_k[0] = shape.nodeCoords[1][0];
  Local2Global(c1, xi_k);

  diffVecWholeElem = c1 - c0;

  lengthOfWholeElem = diffVecWholeElem.NormL2();
  fac = 1.0 / lengthOfWholeElem;

  for (UInt j = 0; j < globDim; j++) {
    diffVecToPoint[j] = globalPoint[j] - c0[j];
  }

  lengthToPoint = diffVecToPoint.NormL2();
  diffVecToPoint.Inner(diffVecWholeElem, s);
  s = fac * lengthToPoint;

  locPoint[0] = shape.nodeCoords[0][0]
      + (shape.nodeCoords[1][0] - shape.nodeCoords[0][0]) * s;
}

void LagrangeElemShapeMap::Global2LocalGeneral(Vector<Double>& locPoint,
    const Vector<Double>& globalPoint) {
  
  const ElemShape & shape = *shape_;
  //UInt globDim = globalPoint.GetSize(); // determine global dimension
  UInt locDim = shape.dim; // dimension of current element

  Vector<Double> xi_start; // local start point for Newton-Raphson method
  Vector<Double> xi_k; // local point at iteration k
  Vector<Double> xi_min; // local point with minimal global distance
  Vector<Double> delta_xi; // local search direction
  Vector<Double> f, f1, f2; // global (negative) search direction
  Vector<Double> f_start; // start value for global search direction
  Matrix<Double> J; // Jacobian at local point xi_k
  Double f_min; // minimal distance between global corners and global point
  Double distance; // global distance!
  Double distance_l; // global distance at iteration l
  Double distMeasure; // global distance or local distance
  // depending on sqrt(|jacDet|)
  //Double TOL = 0.0000001; // tolerance for global distance
  Double TOL = 1e-5;
  Double jacDet; // denominator for Cramer's rule.
  Double distNormalizer; // denominator for Cramer's rule.
  // bool divergence; // does the Newton-Raphson algorithm diverge?
  bool converged; // have we found the local point?
  UInt iter = 0;
  const Double golden_ratio = (3 - sqrt(5)) / 2;
  Double minDist = 999999;

  // Initialize variables
  xi_start.Resize(3);
  xi_start.Init();
  xi_k.Resize(3);
  xi_k.Init();
  delta_xi.Resize(3);
  delta_xi.Init();
  J.Resize(3, 3);
  J.Init();

  // Perform Newton-Raphson method on global point
  // Find good startpoint xi_k among local node coordinates
  f_min = 999e5; // really big value!
  for (UInt k = 0; k < shape.numNodes; k++) {
    for (UInt l = 0; l < locDim; l++) {
      xi_k[l] = shape.nodeCoords[k][l];
    }

    Local2Global(f, xi_k);
    f = f - globalPoint;
    distance = f.NormL2();
    if (distance < f_min) {
      f_min = distance;
      xi_start = xi_k;
      f_start = f;
      minDist = distance;
    }
  }
  xi_k = xi_start;
  xi_min = xi_start;
  f = f_start;
  distance = f_min;
  // divergence = false;
  converged = false;

  do {
    // Calculate Jacobian at iteration point xi_k
    CalcJ(J, xi_k);

    // locDim should never be 1 since line elements will
    // be handled differently
    if (locDim == 1) {
      EXCEPTION("Line elements should not use the Newton-Raphson method!");
      return;
    }
    
    jacDet = GetLocDirJac(delta_xi, f, J);
    distNormalizer = sqrt(fabs(jacDet));
//
//    if (globDim == 2) {
//      // Find new local search direction for 2D -> 2D mapping using
//      // Cramer's rule.
//
//      // Jacobian Matrix:
//      //
//      //  J = ( J_00  J_01 )
//      //      ( J_10  J_11 )
//
//      jacDet = +J[0][0] * J[1][1] - J[1][0] * J[0][1];
//
//      //  ( f_0  J_01 )
//      //  ( f_1  J_11 )
//      delta_xi[0] = -J[1][1] * f[0] + J[0][1] * f[1];
//
//      //  ( J_00  f_0 )
//      //  ( J_10  f_1 )
//      delta_xi[1] = -J[0][0] * f[1] + J[1][0] * f[0];
//
//      distNormalizer = sqrt(fabs(jacDet));
//    } else {
//      if (locDim == 2) {
//        // Find new local search direction for 2D -> 3D mapping.
//        // Project 3D difference vector onto 2D basis given by
//        // the Jacobian to find the new 2D search direction.
//
//        // Jacobian Matrix:
//        //
//        //      ( J_00  J_01  1 )
//        //  J = ( J_10  J_11  1 )
//        //      ( J_20  J_21  1 )
//
//        jacDet = +J[1][0] * J[2][1] - J[2][0] * J[1][1] - J[0][0] * J[2][1]
//            + J[2][0] * J[0][1] + J[0][0] * J[1][1] - J[1][0] * J[0][1];
//
//        // Calculate negative Jacobian determinants of the following matrices
//        // to find the local search direction which has to point in the opposite
//        // direction of the backprojected global error vector f.
//        //
//        //  ( f_0  J_01  1 )
//        //  ( f_1  J_11  1 )
//        //  ( f_2  J_21  1 )
//
//        delta_xi[0] = -f[0] * J[1][1] - f[1] * J[2][1] - f[2] * J[0][1]
//            + f[2] * J[1][1] + f[1] * J[0][1] + f[0] * J[2][1];
//
//        //  ( J_00  f_0  1 )
//        //  ( J_10  f_1  1 )
//        //  ( J_20  f_2  1 )
//
//        delta_xi[1] = -f[0] * J[2][0] - f[1] * J[0][0] - f[2] * J[1][0]
//            + f[2] * J[0][0] + f[1] * J[2][0] + f[0] * J[1][0];
//
//        distNormalizer = sqrt(fabs(jacDet));
//      } else {
//        // Find new local search direction for 3D -> 3D mapping using
//        // Cramer's rule.
//
//        //      ( J_00  J_01  J_02 )
//        //  J = ( J_10  J_11  J_12 )
//        //      ( J_20  J_21  J_22 )
//
//        jacDet = +J[0][0] * J[1][1] * J[2][2] - J[0][0] * J[2][1] * J[1][2]
//            - J[1][0] * J[0][1] * J[2][2] + J[2][1] * J[1][0] * J[0][2]
//            - J[1][1] * J[2][0] * J[0][2] + J[2][0] * J[0][1] * J[1][2];
//
//        //  ( f_0  J_01  J_02 )
//        //  ( f_1  J_11  J_12 )
//        //  ( f_2  J_21  J_22 )
//        delta_xi[0] = -J[1][1] * J[2][2] * f[0] + J[2][1] * J[1][2] * f[0]
//            - J[2][1] * J[0][2] * f[1] + J[0][1] * J[2][2] * f[1]
//            - J[0][1] * J[1][2] * f[2] + J[1][1] * J[0][2] * f[2];
//
//        //  ( J_00  f_0  J_02 )
//        //  ( J_10  f_1  J_12 )
//        //  ( J_20  f_2  J_22 )
//        delta_xi[1] = -J[1][2] * J[2][0] * f[0] + J[1][0] * J[2][2] * f[0]
//            - J[0][0] * J[2][2] * f[1] + J[2][0] * J[0][2] * f[1]
//            - J[1][0] * J[0][2] * f[2] + J[1][2] * J[0][0] * f[2];
//
//        //  ( J_00  J_01  f_0 )
//        //  ( J_10  J_11  f_1 )
//        //  ( J_20  J_21  f_2 )
//        delta_xi[2] = -J[1][0] * J[2][1] * f[0] + J[2][0] * J[1][1] * f[0]
//            - J[0][1] * J[2][0] * f[1] + J[2][1] * J[0][0] * f[1]
//            - J[0][0] * J[1][1] * f[2] + J[1][0] * J[0][1] * f[2];
//
//        distNormalizer = std::pow(fabs(jacDet), 1.0 / 3.0);
//      }
//    }

    // Here is the new local search direction. We normalize it so we
    // can be sure, that we have a local search vector of a length
    // comparable to the local element diameter.
    Double len;

//    delta_xi[0] /= jacDet;
//    delta_xi[1] /= jacDet;
//    delta_xi[2] /= jacDet;

    len = delta_xi.NormL2();
    delta_xi[0] /= len;
    delta_xi[1] /= len;
    delta_xi[2] /= len;

    // If global element is smaller use local distance as a measure.
    // If global element is bigger use global distance as a measure.
    distMeasure = distNormalizer < 1 ? distance / distNormalizer : distance;
    if (distMeasure < TOL) {
      converged = true;
      break;
    }

    // Perform damping iterations to find good damping coefficient.
    // That means we search for a factor along the local search direction
    // which minimizes the global error vector. We do it by braketing the
    // the interval for the factor until we encounter a case were the
    // global error vectors face into opposite directions. This means
    // that the point we search lies somewhere between the mapped
    // interval limits.
    xi_start = xi_k;
    Double interval[2];
    Double dist[2];
    interval[1] = 1.0;
    interval[0] = 0.0;
    dist[0] = distance;

    // Braket the location of the local point
    UInt l = 0;

    for (l = 0; l < 20; l++) {
      xi_k = xi_start + delta_xi * interval[1];
      Local2Global(f2, xi_k);
      f2 = f2 - globalPoint;
      dist[1] = f2.NormL2();

      distMeasure = distNormalizer < 1 ? dist[1] / distNormalizer : dist[1];
      if (distMeasure < TOL) {
        converged = true;
        break;
      }

      if (f.Inner(f2) < 0) {
        break;
      } else {
        interval[0] = interval[1];
        interval[1] *= 4.0;
        dist[0] = dist[1];
        f = f2;
      }

      if (dist[0] < minDist) {
        minDist = dist[0];
        xi_min = xi_start + delta_xi * interval[0];
      }

      if (dist[1] < minDist) {
        minDist = dist[1];
        xi_min = xi_start + delta_xi * interval[1];
      }
    }

    if (converged)
      break;

    // Try to narrow down the extents of the interval by searching points
    // ever closer from below and above the minumum in the current search
    // direction.
    for (l = 0; l < 20; l++) {
      Double x3 = interval[0] + (interval[1] - interval[0]) * golden_ratio;
      Double x4 = interval[1] - (interval[1] - interval[0]) * golden_ratio;

      distMeasure = distNormalizer < 1 ? dist[0] / distNormalizer : distance;
      if (distMeasure < TOL) {
        converged = true;
        break;
      }

      distMeasure = distNormalizer < 1 ? dist[1] / distNormalizer : distance;
      if (distMeasure < TOL) {
        converged = true;
        break;
      }

      xi_k = xi_start + delta_xi * x3;
      Local2Global(f1, xi_k);

      f1 = f1 - globalPoint;

      xi_k = xi_start + delta_xi * x4;
      Local2Global(f2, xi_k);

      f2 = f2 - globalPoint;

      // If both points lie on the same side of the searched point try to find
      // points on opposite sides. Note that this check is not completely clean
      // since a line in the local coordinate system might get mapped to an
      // arbitrary curve in the global coordinate system and thus both vectors
      // might point into slightly different directions. We assume however that
      // we are close to a minimum and that the vectors point into nearly the
      // same direction.
      //
      // Situation in global space:
      //
      //          error vector 1      error vector 2
      //        o---------------> X <----------------o
      // 
      //   interval[0]        searched        interval[1]
      //                       point

      if (f1.Inner(f2) > 0) {
        // Either x3 switched to the side of x4 in respect to the searched
        // point or x4 switched to the side of x3.
        UInt m;
        const UInt n = 16;

        for (m = 0; m < n && f1.Inner(f2) > 0; m++) {
          x3 -= (interval[1] - interval[0]) * (golden_ratio / n);

          xi_k = xi_start + delta_xi * x3;
          Local2Global(f1, xi_k);

          f1 = f1 - globalPoint;
        }

        if (m == n) {
          for (m = 0; m < n && f1.Inner(f2) > 0; m++) {
            x4 += (interval[1] - interval[0]) * (golden_ratio / n);

            xi_k = xi_start + delta_xi * x4;
            Local2Global(f2, xi_k);

            f2 = f2 - globalPoint;
          }
        } else {
          for (m = 0; m < n && f1.Inner(f2) < 0; m++) {
            x4 -= (interval[1] - interval[0]) * (golden_ratio / n);

            xi_k = xi_start + delta_xi * x4;
            Local2Global(f2, xi_k);

            f2 = f2 - globalPoint;
          }
        }

        if (x4 < x3)
          x4 += (interval[1] - interval[0]) * (golden_ratio / n);
      }

      // Update interval and distance array with new values.
      distance_l = f1.NormL2();

      if (distance_l < dist[0]) {
        interval[0] = x3;
        dist[0] = distance_l;
      }

      distance_l = f2.NormL2();

      if (distance_l < dist[1]) {
        interval[1] = x4;
        dist[1] = distance_l;
      }

      if (dist[0] < minDist) {
        minDist = dist[0];
        xi_min = xi_start + delta_xi * interval[0];
      }

      if (dist[1] < minDist) {
        minDist = dist[1];
        xi_min = xi_start + delta_xi * interval[1];
      }
    }

    // If distances increase again, break
    if (dist[0] > minDist && dist[1] > minDist) {
      xi_k = xi_min;
      break;
    }

    if (dist[0] < dist[1])
      xi_k = xi_start + delta_xi * interval[0];
    else
      xi_k = xi_start + delta_xi * interval[1];

    distance = minDist;
    xi_start = xi_k;
    Local2Global(f, xi_k);
    f = f - globalPoint;

    if (converged)
      break;

    iter++;
  } while (iter < 6);

  if (xi_k.GetSize() == 0) {
    WARN("Global2Local lost its memory, Resetting the point to 99");
    xi_k.Resize(locDim);
    xi_k.Init(99);
  }

//  Vector<Double> globMapped, errVec;
//  Local2Global(globMapped, xi_k);
//  errVec = globMapped - globalPoint;
  locPoint = xi_k;
}

void LagrangeElemShapeMap::GetGlobMidPoint(Vector<Double>& midPoint) {
  Local2Global(midPoint, shape_->midPointCoord);
}

Double LagrangeElemShapeMap::CalcVolume(bool useDepth) {

  // Get integration points
  StdVector<LocPoint> intPoints;
  StdVector<Double> weights;

  // Order: use element order and add 2 to be sure for curved elements
  UInt order = shape_->order + 2;
  intScheme_->GetIntPoints(Elem::GetShapeType(ptElem_->type), IntScheme::GAUSS,
      order, intPoints, weights);
  Double vol = 0.0;
  Double jacDet = 0.0;
  Matrix<Double> jac;
  // Loop over all integration points
  LocPointMapped lpm;
  for (UInt i = 0; i < intPoints.GetSize(); i++) {
    jacDet = CalcJDet(jac, intPoints[i],useDepth);
    vol += jacDet * weights[i];
  }
  return vol;
}

void LagrangeElemShapeMap::CalcNormal(Vector<Double>& normal,
    const LocPoint& lp) {

  // check, that element is a surface element at all
  if (shape_->dim != ptGrid_->GetDim() - 1) {
    EXCEPTION(
        "Can not calculate normal of element #" << ptElem_->elemNum << " which is of dimension " 
        << shape_->dim << " in a " << ptGrid_->GetDim() << "-dimensional grid!");
  }

  // Get neighboring volume element. 
  // Here we always use the first neighbor, so the resulting
  // normal will point OUT of the first volume neighor
  const SurfElem & surfEl = *ptSurfElem_;
  Elem * ptVolEl = surfEl.ptVolElems[0];
  assert(ptVolEl);

  // Obtain shape map of neighboring volume element
  LagrangeElemShapeMap sm(ptGrid_);
  sm.SetElem(ptVolEl, isUpdated_);

  // Map local point of surface to global point and obtain
  // the normal vector (pointing OUT of the volume element)
  // in local coordinates (w.r.t to the volume element)
  Vector<Double> locNormal;
  LocPoint volPoint;
  sm.ptFe_->GetLocalIntPoints4Surface(ptElem_->connect, ptVolEl->connect, lp,
      volPoint, locNormal);
//  std::cerr << "\nlp = " << lp.coord.ToString() << std::endl;
//  std::cerr << "volPoint = " << volPoint.coord.ToString() << std::endl;
//  std::cerr << "locNormal = " << locNormal.ToString() << std::endl;

  // Calculate Jacobian matrix of volume element in that point
  Matrix<Double> jac;
  sm.CalcJ(jac, volPoint);

  // Calculate global normal
  Matrix<Double> jInv;
  jac.Invert(jInv);
  normal = Transpose(jInv) * locNormal;

  // normalize normal
  Double norm = normal.NormL2();
  normal /= norm;
}

bool LagrangeElemShapeMap::CalcNormalOutOfVolume(Vector<Double> & normal,
    const LocPoint & lp,
    const Elem* volElem,
    const SurfElem* edgeFaceElem) {
  //now we check for the element dimension
  //in 2D its quite simple, loop over the edges
  //and check if the local point is contained in the edge
  //in 3D it gets harder but we will come to this later
  const ElemShape & shape = *shape_;
  bool success = false;
  bool edgeFaceFound = false;
  
  if(shape.dim == 2){
    //for the 2D case, we just determine the correct side of the
    //volume element by checking the edge number
    UInt locENum = 0;
    UInt numEdges = volElem->extended->edges.GetSize();

    for(locENum = 0; locENum < numEdges; locENum++){
      if(abs(edgeFaceElem->extended->edges[0]) == abs(volElem->extended->edges[locENum])){
        edgeFaceFound = true;
        break;
      }
    }
    if(!edgeFaceFound){
      WARN("cannot determine corresponding edge to volume element for normal computation");
    }
    //now we have the local edge number, lets get the connectivity and compute the normal
    //get Vertices of current edge
    StdVector<UInt> eVert = shape.edgeVertices[locENum];
    Vector<Double> c1 = shape.nodeCoords[eVert[0]-1];
    Vector<Double> c2 = shape.nodeCoords[eVert[1]-1];
    Vector<Double> diff;
    diff = c2 - c1;
    normal.Resize(2,0.0);

    normal[0] = diff[1];
    normal[1] = -diff[0];

    success = true;
  }else{
    //loop over faces
    UInt numFaces = shape.numFaces;
    UInt locFaceNum = 0;
    //determine on which face the local point is located and the compute the normal for this face
    for(locFaceNum = 0; locFaceNum < numFaces ; locFaceNum++){
      if(abs(edgeFaceElem->extended->faces[0]) == abs(volElem->extended->faces[locFaceNum])){
        edgeFaceFound = true;
        break;
      }
    }
    if(!edgeFaceFound){
      WARN("cannot determine corresponding face to volume element for normal computation");
    }
    StdVector<UInt> fVert = shape.faceVertices[locFaceNum];
    //take the first three vertices to span our surface (we will always have at least three vertices
    Vector<Double> v1(3),v2(3),v3(3);
    v1 = shape.nodeCoords[fVert[1]-1];
    v2 = shape.nodeCoords[fVert[0]-1];
    v3 = shape.nodeCoords[fVert[2]-1];
    Vector<Double> c1(3);
    Vector<Double> c2(3);
    c1 = v1 - v2;
    c2 = v1 - v3;
    //compute cross product
    normal.Resize(3,0.0);
    c1.CrossProduct(c2,normal);
    //normalize
    normal /= normal.NormL2();

    success = true;
  }

  // make sure the normal vector points outside of the element
  // the inner product of normal vector and center_to_surface vector is negative if they point into different directions
  // multiplying with it will flop the direction of the normal vector
  Vector<Double> center_to_surface;
  center_to_surface = lp.coord - shape.midPointCoord;
  normal *= normal * center_to_surface;

  if(!success){
    WARN("could not determine surface normal.. to be checked!")
  }
  return success;
}

// The follogin method is not needed anymore. 
//void  LagrangeElemShapeMap::
//CalcNormalOutOfVol( Vector<Double> & normal,
//                    const LocPoint& lp, 
//                    const Elem & volElem ) {
//
//  // Obtain shape map of neighboring volume element
//  LagrangeElemShapeMap sm(ptGrid_);
//  sm.SetElem(&volElem, isUpdated_);
//  
//  Matrix<Double> & volCoords = sm.coords_; 
//
//  // Calculate surface normal without defined sign
//  CalcNormal( normal, lp);
//  
//  std::cerr << "====== In CalcNormal out of Vol for surfElem " << ptElem_->elemNum << "======\n";  
//  std::cerr << "unoriented normal: " << normal.ToString() << std::endl;
//  UInt volCorners = sm.shape_.numVertices;
//
//  /* Idea:
//       - find a common surface / volume element node
//       - find additional node, which is not common (i.e. lies on the "volume side"
//       - calculate difference vector from  "additional" node to first common
//         one ( = pointing out of  the volume element)
//       - the scalar product of the normal and the difference vector determines
//         the normalSign
//   */
//
//  // find first common vertex index
//  Integer firstCommonIndex = -1;
//  for (UInt i=0; i<volCorners; i++)
//    if (volElem.connect[i] == ptElem_->connect[0]){
//      firstCommonIndex = i;
//      break;
//    }
//
//  // find additional node of the volume node, which is not contained in the
//  // surface element
//  std::set<UInt> volSet, surfSet, diffSet;
//  volSet = std::set<UInt>(sm.ptElem_->connect.Begin(), 
//                          sm.ptElem_->connect.End());
//  
//  surfSet = std::set<UInt>(ptElem_->connect.Begin(), 
//                           ptElem_->connect.End());
//  std::set_difference( volSet.begin(), volSet.end(),
//                       surfSet.begin(), surfSet.end(),
//                       std::inserter(diffSet,diffSet.begin()) );
//
//  Vector<double> diffVec( sm.shape_.dim );
//  Double scalarProd = 0.0;
//  std::set<UInt>::const_iterator it = diffSet.begin();
//  for( ; it != diffSet.end(); it++ ) {
//
//    UInt node = *it;
//
//    // search for volume 
//    UInt index = volElem.connect.Find(node);
//
//    // calculate difference vector (pointing out of the volume)
//    for( UInt iDim = 0; iDim < sm.shape_.dim; ++iDim ) {
//      diffVec[iDim] =  volCoords[iDim][firstCommonIndex]
//                     - volCoords[iDim][index];
//    }
//    // normalize difference vector to 1.0
//    diffVec /= diffVec.NormL2();
//
//    // calculate scalar product
//    scalarProd = diffVec * normal;
//
//    // check if scalar product is != 0 (otherwise we may have a degenerated
//    // element, where two nodes lie on the same location)
//    if( std::abs(scalarProd) < EPS ) {
//      // we have to find another node
//      continue;
//    } else {
//      if( scalarProd < 0.0 ) {
//        // original normal vector points into the volume -> reorient
//        normal *= -1.0;
//      }
//      break;
//    }
//  }
//  std::cerr << "oriented normal: " << normal.ToString() << "\n\n";
//}

void LagrangeElemShapeMap::GetLocalIntPoints4Surface( const StdVector<UInt> & surfConnect, const LocPoint & surfIntPoint,
  LocPoint & volIntPoint, Vector<Double>& locNormal)
{
  ptFe_->GetLocalIntPoints4Surface(surfConnect, ptElem_->connect, surfIntPoint, volIntPoint, locNormal);
}

void LagrangeElemShapeMap::MapSurfLocDirs(const Elem* ptSurfElem,
  StdVector<UInt>& surfLocDirs)
{
  const ElemShape & shape = *shape_;
  // determine dimension of element
  // 1: look for edges
  // 2: look for faces

  UInt surfDim = Elem::shapes[ptSurfElem->type].dim;
  assert( surfDim < shape.dim);
  surfLocDirs.Resize(surfDim);

  if (surfDim == 1) {
    // -------------
    //  Common Edge
    // -------------
    // look for common edge
    // A 1D surface element only has 1 edge
    UInt edgeNum = std::abs(ptSurfElem->extended->edges[0]);
    Integer index = ptElem_->extended->edges.Find(edgeNum);
    if (index < 0) {
      index = ptElem_->extended->edges.Find(-Integer(edgeNum));
      if (index < 0) {
        EXCEPTION("edge not found");
      }
    }
    surfLocDirs[0] = shape.edgeLocDirs[index];

  } else if (surfDim == 2) {
    // -------------
    //  Common Face
    // -------------
    // a 2D surface element only has one face
    UInt faceNum = std::abs(ptSurfElem->extended->faces[0]);
    Integer index = ptElem_->extended->faces.Find(faceNum);
    if (index < 0) {
      EXCEPTION("face not found");
    }
    surfLocDirs[0] = shape.faceLocDirs[index][0];
    surfLocDirs[1] = shape.faceLocDirs[index][1];

    // as we have identified the correct face, we have already the two 
    // involved global directions. In order to determine the correct ordering
    // of the directions, we compare the orientation of the first edge of each
    // face.
    
    // get local volume nodes of surface edge #1, pointing in local xi/0 direction  
    Integer volNode1, volNode2;
    volNode1 = ptElem_->connect.Find(ptSurfElem->connect[0])+1;
    volNode2 = ptElem_->connect.Find(ptSurfElem->connect[1])+1;

    Integer edgeIndex = -1;
    for( UInt i = 0; i < shape.numEdges; ++i ) {
      const StdVector<UInt> & edgeNodes = shape.edgeNodes[i];
      if( (edgeNodes[0] == UInt(volNode1) && edgeNodes[1] == UInt(volNode2) ) ||
          (edgeNodes[1] == UInt(volNode1) && edgeNodes[0] == UInt(volNode2) ) ) {
        edgeIndex = i;
        break;
      }
    }
    if(UInt(shape.edgeLocDirs[UInt(edgeIndex)]) != surfLocDirs[0]) {
      std::swap(surfLocDirs[0], surfLocDirs[1]);
    }
  } else {
    EXCEPTION("Can only handle 1D or 2D elements.")
  }

}

bool LagrangeElemShapeMap::CoordIsInsideElem(const Vector<Double>& point,
  Double tolerance) {

return ptFe_->CoordIsInsideElem(point, tolerance);
}

void LagrangeElemShapeMap::CalcDiameter(Vector<Double>& diameter) {
  Vector<Double> mins(shape_->dim), maxs(shape_->dim);
  mins.Init( std::numeric_limits<double>::max());
  maxs.Init(-std::numeric_limits<double>::max());
  diameter.Resize(shape_->dim);
  diameter.Init();
  for (UInt dim = 0; dim < shape_->dim; dim++) {
    for (UInt k = 0, n_elems = coords_.GetNumCols(); k < n_elems; k++) {
      Double test = coords_[dim][k];
      Double& min = mins[dim];
      Double& max = maxs[dim];
      min = std::min(min, test);
      max = std::max(max, test);
    }

    diameter[dim] = maxs[dim] - mins[dim];
  }
}

void LagrangeElemShapeMap::CalcBarycenter(Point& barycenter)
{
  CalcBarycenter(barycenter, coords_, domain);
}

void LagrangeElemShapeMap::CalcBarycenter(Point& barycenter, const Matrix<double>& coords, const Domain* domain)
{
  UInt n_dims  = coords.GetNumRows();
  UInt n_elems = coords.GetNumCols();

  //sometimes there is no domain pointer
  if(domain){
    assert(n_dims == (unsigned int) domain->GetDim());
  }

  // init barycenter for safty reason
  barycenter.SetZero();

  // std::cout << "calc a new barycenter" << std::endl;
  // a barycenter is simply the average of all coordinates
  for (UInt dim=0; dim < n_dims; dim++)
  {
    // std::cout << "dim = " << dim << "  ";
    for (UInt k=0; k < n_elems; k++)
    {
      // the constructor of Point initializes
      barycenter[dim] += coords[dim][k];
      // std::cout << coords[dim][k] << "->" << barycenter[dim] << "\t";
    }

    barycenter[dim] /= (double) n_elems;
    // std::cout << " average: " << (barycenter[dim]) << std::endl;
  }
}

void LagrangeElemShapeMap::GetMaxMinEdgeLength(Double& max, Double& min) {
  const ElemShape & shape = *shape_;
  max = 1e-100;
  min = 1e+100;
  Double length, dl;
  //check for surface element
  UInt dim = std::max(this->ptGrid_->GetDim(), shape.dim);
  for (UInt i = 0; i < shape.numEdges; ++i) {
    length = 0.0;
    for (UInt iDim = 0; iDim < dim; ++iDim) {
      dl = coords_[iDim][shape.edgeVertices[i][1] - 1] - coords_[iDim][shape.edgeVertices[i][0] - 1];
      length += dl * dl;
    }
    length = sqrt(length);
    max = max > length ? max : length;
    min = min < length ? min : length;
  }
}

void LagrangeElemShapeMap::GetEdgeLength(StdVector<Double>& edges_out)
{
  const ElemShape & shape = *shape_;

  assert(Elem::GetShapeType(ptElem_->type) == Elem::ST_QUAD || Elem::GetShapeType(ptElem_->type) == Elem::ST_HEXA);

  edges_out.Resize(shape.dim, 0.0);

//  for (UInt i = 0; i < shape.numEdges; ++i)
//    for (UInt iDim = 0; iDim < shape.dim; ++iDim){
//      std::cout << "nodes of edge " << i << ": " << shape.edgeVertices[i][0] << "," << shape.edgeVertices[i][1] << std::endl;
//      edges_out[iDim] += abs(coords_[iDim][shape.edgeVertices[i][1] - 1] - coords_[iDim][shape.edgeVertices[i][0] - 1]) / shape.dim;
//      std::cout <<"idim " << iDim << ":" <<  edges_out[iDim] << "+=" <<  abs(coords_[iDim][shape.edgeVertices[i][1] - 1] - coords_[iDim][shape.edgeVertices[i][0] - 1])  << "/" << shape.dim << std::endl;
//
//      std::cout << "edge number: " << i << " edges_out[" << iDim <<"]=" << edges_out[iDim] << std::endl;
//    }
//
//   std::cout << "edges=" << edges_out.ToString() << std::endl;
   // maybe computation is not correct but I do not see where we need a loop?
   // see RectangleFE::GetEdgeLength()

   for(UInt i = 0; i < shape.dim; ++i)
   {
     // for all dimensions, add offset; only in one direction this is not 0
     edges_out[i]  = abs(coords_[i][0] - coords_[i][1]);
     edges_out[i] += abs(coords_[i][0] - coords_[i][3]);
     if (shape.dim == 3)
       edges_out[i] += abs(coords_[i][0] - coords_[i][4]);
   }
}

void LagrangeElemShapeMap::GetExtensionLocalDir(Vector<Double>& extension) {
  const Double MIN = 1e100;
  const ElemShape & shape = *shape_;
  Vector<Double> min(shape.dim);
  min.Init(MIN);
  Double length, dl;
  for (UInt i = 0; i < shape.numEdges; ++i) {
    length = 0.0;
    for (UInt iDim = 0; iDim < shape.dim; ++iDim) {
      dl = coords_[iDim][shape.edgeVertices[i][1] - 1]
                         - coords_[iDim][shape.edgeVertices[i][0] - 1];
      length += dl * dl;
    }
    length = sqrt(length);
    Integer locDir = shape.edgeLocDirs[i];
    if (locDir > -1) {
      if (length < min[locDir]) {
        min[locDir] = length;
      }
    }
  }

  extension.Resize(shape.dim);
  extension.Init(0.0);
  for (UInt i = 0; i < shape.dim; ++i) {
    if (min[i] != MIN)
      extension[i] = min[i];
  }
}

/*
 * 21.04.2020
 * > imported from NACS
 */
void LagrangeElemShapeMap::SetModelDepth(Double depth) {
  if (depth < 0.0) {
    EXCEPTION("Depth of model must be > 0.")
  }
//  std::cout << "Model depth set to value " << depth << std::endl;
  this->depth_ = depth;
}

Double LagrangeElemShapeMap::GetModelDepth(){
  return this->depth_;
}

void LagrangeElemShapeMap::CalcJ(Matrix<Double>& jac, const LocPoint& lp) {
  // --- SPECIAL CASE: rational pyramid --- //
  if (ptElem_->type == Elem::ET_PYRA5 ||
      ptElem_->type == Elem::ET_PYRA13 ||
      ptElem_->type == Elem::ET_PYRA14)
  {
    double x = lp.coord[0];
    double y = lp.coord[1];
    double z = lp.coord[2];

    double lam[5];
    double dlam[5][3];
    PyraLamAndGrads(x,y,z, lam, dlam);

    jac.Resize(3,3);
    jac.Init();

    const UInt nNodes = std::min<UInt>(coords_.GetNumCols(), 5);
    for (UInt a=0; a<nNodes; ++a)
    {
      jac[0][0] += coords_[0][a] * dlam[a][0];
      jac[0][1] += coords_[0][a] * dlam[a][1];
      jac[0][2] += coords_[0][a] * dlam[a][2];

      jac[1][0] += coords_[1][a] * dlam[a][0];
      jac[1][1] += coords_[1][a] * dlam[a][1];
      jac[1][2] += coords_[1][a] * dlam[a][2];

      jac[2][0] += coords_[2][a] * dlam[a][0];
      jac[2][1] += coords_[2][a] * dlam[a][1];
      jac[2][2] += coords_[2][a] * dlam[a][2];
    }
  }else{


    jac = coords_ * ptFe_->GetLocDerivShFnc( lp, ptElem_);
  //  jac *= depth_; // explicitly include depth_ of setup
    /*
    * Note: Scaling of jacobian not correct; only scale its determinant
    * Explanation (thanks Jens!):  The Jacobian determinant is part of every integral form
    * and thus causes the integrals to consider the "correct" depth of the 2d plane setup.
    * In case of differential operations the inverse Jacobian is required in addition for the
    * transform of the differential operator. Applying a scaling with depth_ to the Jacobian 
    * would therewith affect these terms in a stronger way than terms without differential operators,
    * like e.g., terms for the mass matrix. 
    */
  }

}

//! Calculation of Jacobian with given coordinates
void LagrangeElemShapeMap::CalcJ( Matrix<Double>& jac,
       		                     const LocPoint& lp,
   				                 Matrix<Double>& cornerCoords) {
	jac = cornerCoords * ptFe_->GetLocDerivShFnc( lp, ptElem_);
//	jac *= depth_; // explicitly include depth_ of setup
}

/*
 * Ported from NACS
 * difference to CFS:
 * 1. depth of setup optional (makes no difference if default depth of 1m is used
 * 2. scaling for axi-symmetric case NOT included
 * 3. 1d elems in 3d included here
 */
Double LagrangeElemShapeMap::CalcJDet(Matrix<Double>& jac,
                                      const LocPoint& lp,
                                      bool useDepth)
{
  // --- SPECIAL CASE: rational pyramid --- //
  if (ptElem_->type == Elem::ET_PYRA5 ||
      ptElem_->type == Elem::ET_PYRA13 ||
      ptElem_->type == Elem::ET_PYRA14)
  {
    // use rational lambdas instead of polynomial Lagrange derivs
    double x = lp.coord[0];
    double y = lp.coord[1];
    double z = lp.coord[2];

    double lam[5];
    double dlam[5][3];
    PyraLamAndGrads(x,y,z, lam, dlam);

    // coords_: dim x nNodes (dim=3), nNodes>=5
    jac.Resize(3,3);
    jac.Init();

    const UInt nNodes = std::min<UInt>(coords_.GetNumCols(), 5);
    for (UInt a=0; a<nNodes; ++a)
    {
      jac[0][0] += coords_[0][a] * dlam[a][0];
      jac[0][1] += coords_[0][a] * dlam[a][1];
      jac[0][2] += coords_[0][a] * dlam[a][2];

      jac[1][0] += coords_[1][a] * dlam[a][0];
      jac[1][1] += coords_[1][a] * dlam[a][1];
      jac[1][2] += coords_[1][a] * dlam[a][2];

      jac[2][0] += coords_[2][a] * dlam[a][0];
      jac[2][1] += coords_[2][a] * dlam[a][1];
      jac[2][2] += coords_[2][a] * dlam[a][2];
    }

    double jacDet = 0.0;
    jac.Determinant(jacDet);

    if (useDepth) jacDet *= depth_;

    return jacDet;
  }else{


    deriv_ = ptFe_->GetLocDerivShFnc( lp, ptElem_);
    jac = coords_ * deriv_;
    Double jacDet = 0.0;

    if (jac.GetNumCols() == jac.GetNumRows()) {
      jac.Determinant(jacDet);
    } else if (jac.GetNumRows() == 3 && jac.GetNumCols() == 2) {
      // 2D elements in 3D
      Vector<Double> normal;
      normal.Resize(3);
      normal[0] = jac[1][0] * jac[2][1] - jac[2][0] * jac[1][1];
      normal[1] = jac[2][0] * jac[0][1] - jac[0][0] * jac[2][1];
      normal[2] = jac[0][0] * jac[1][1] - jac[1][0] * jac[0][1];
      jacDet = sqrt(
          normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]);
    } else if (jac.GetNumRows() == 3 && jac.GetNumCols() == 1) {
      // === 1D elements in 3D ===
      jacDet = sqrt( jac[0][0] * jac[0][0]
                    + jac[1][0] * jac[1][0]
                    + jac[2][0] * jac[2][0]);
    } else if (jac.GetNumRows() == 2) {
      // 1D elements in 2D
      //see kaltenbacher, p.23, eq.(2.122)
      jacDet = sqrt(jac[0][0] * jac[0][0] + jac[1][0] * jac[1][0]);
    }
    if (useDepth) {
      jacDet *= depth_; // explicitly include depth_ of setup Note: but ONLY to determinant!
    }
    
    return jacDet;
  }
}

BaseFE* LagrangeElemShapeMap::GetBaseFE() {
  return ptFe_;
}

std::string LagrangeElemShapeMap::ToString() const
{
  std::stringstream ss;
  ss << ElemShapeMap::ToString();
  ss << " co=" << coords_.ToString();
  return ss.str();
}

void LagrangeElemShapeMap::SetElem(const Elem* ptElem, bool isUpdated)
{
  ElemShapeMap::SetElem(ptElem, isUpdated);
  //  ptElem_ = ptElem;
  //  isUpdated_ = isUpdated;
  //  isAxi_ = ptGrid_->IsAxi();
  //
  // call setElem at base class

  // get coordinates from grid
  ptGrid_->GetElemNodesCoord(coords_, ptElem->connect, isUpdated_);
  //  std::cerr << "**** Coordinates for element " << ptElem->elemNum
  //      << " with connect " << ptElem->connect.ToString()
  //
  //      << "(" << (isUpdated ? "updated)" : "original)") << std::endl
  //      << coords_ << std::endl;

  // set reference element
  #ifndef NDEBUG
    if( elems_.feMap_.find(ptElem->type) == elems_.feMap_.end())
      EXCEPTION("Element of type '" << Elem::feType.ToString(ptElem->type) << "' not defined for Lagrangian Shape Map!");
  #endif

  ptFe_ = elems_.feMap_[ptElem->type];
  shape_ = &Elem::shapes[ptElem_->type];
}

void LagrangeElemShapeMap::SetElem(const Elem* ptElem, const Matrix<double>& coords)
{
  this->coords_ = coords;
  assert(elems_.feMap_.find(ptElem->type) != elems_.feMap_.end());
  ptFe_ = elems_.feMap_[ptElem->type];
  shape_ = &Elem::shapes[ptElem_->type];
}

} // namespace CoupledField
