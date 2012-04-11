#include <boost/assign/list_of.hpp>

#include "ElemShapeMap.hh"
#include "Domain/Mesh/Grid.hh"
#include "FeBasis/H1/H1ElemsLagExpl.hh"

namespace CoupledField {

// ===========================================================================
//  C L A S S  LocPoint
// ===========================================================================

  LocPoint::LocPoint() 
  : number( NOT_SET ) {
  }
  
  LocPoint::LocPoint(const Vector<Double>& vec ) {
    number = NOT_SET;
    coord = vec;
  }
  
  std::ostream& operator << ( std::ostream& out, 
                              const LocPoint& lp ) {
    out << "number = " << lp.number << std::endl;
    out << "vector: " << lp.coord.ToString();
    return out;
  }

// ===========================================================================
//  C L A S S   LocPointMapped
// ===========================================================================
  LocPointMapped::LocPointMapped()
  : ptEl( NULL ),
    jacDet( 0.0 ),
    isSurface( false ),
    normalFactor( 0.0 )
    
  {

  }
  
  void LocPointMapped::Set( const LocPoint& lp, shared_ptr<ElemShapeMap> esm ) {

    this->shapeMap = esm;
    this->lp = lp;
    this->ptEl = shapeMap->GetElem();
    this->isSurface = false;

    // Calculate Jacobian, its inverse as well as determinant for local point
    esm->CalcJ( jac, lp);

    // The inversion can only be performed in case we have a quadratic Jacobian
    // i.e. the dimension of the element is the dimension of the grid
    if( jac.GetNumCols() == jac.GetNumRows() ) {
      jac.Invert( jacInv);
      jac.Determinant(jacDet);
    } else if ( jac.GetNumRows() == 3) {
      // 2D elements in 3D
      Vector<Double> normal; 
      normal.Resize(3);
      normal[0]= jac[1][0]* jac[2][1]- jac[2][0]* jac[1][1];
      normal[1]=jac[2][0]*jac[0][1]- jac[0][0]* jac[2][1];
      normal[2]= jac[0][0]* jac[1][1]- jac[1][0]*jac[0][1];
      jacDet = sqrt(normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]);

    } else if ( jac.GetNumRows() == 2) {
      // 1D elements in 2D
      //see kaltenbacher, p.23, eq.(2.122)
      jacDet = sqrt(jac[0][0]*jac[0][0] + jac[1][0]*jac[1][0]);
    };

    // Check, if geometry is axi-symmetric. In this case scale the 
    // Jacobian determinant with 2*pi*r
    Vector<Double> globPoint;
    if( esm->IsAxi() ) {
      esm->Local2Global( globPoint, lp);
      jacDet *= 2 * PI * globPoint[0];
    }
  }
  
  void LocPointMapped::Set( const LocPoint& lp, shared_ptr<ElemShapeMap> esm,
                            const std::set<RegionIdType>& myRegions) {


    // ------------------------------------------------------
    //  1) Set "normal" element information (Jacobian, etc.)
    // ------------------------------------------------------

    // first, call "normal" set method, also valid for volume elements
    this->Set(lp, esm);

    // ------------------------------------------------------
    //  2) Perform surface-element specific tasks
    // ------------------------------------------------------
    SetSurfInfo( myRegions );
  }
  
  void LocPointMapped::SetSurfInfo( const std::set<RegionIdType>& myRegions ) {
    // check, if previously set element is the same as this one
    bool isSameElem = (shapeMap->GetElem() == this->ptEl 
        && isSurface == true ); 

    // set flag for surface mapped element
    this->isSurface = true;

    // get surface element
    const SurfElem& surfElem = *(shapeMap->GetSurfElem()); 

    // if we have not previously selected the correct volume neighbor, we
    // have to do it now 
    shared_ptr<ElemShapeMap> esmVol;
    if( !isSameElem) {
      const Elem * ptVolElem = NULL;
      normalFactor = Double(surfElem.normalSign) * -1.0;  
      
      // loop over volume element neighbors of the surface element and check,
      // if the regionId of the element is in the map "myRegions"
      boost::array<Elem*,2>::const_iterator it = surfElem.ptVolElems.begin();
      for( ; it != surfElem.ptVolElems.end(); it++ ) {
        normalFactor *= -1;
        // check, if element is set at all
        if( *it) {
          if(myRegions.find( (*it)->regionId) != myRegions.end()) {
            ptVolElem = *it;
            break;
          }
        }
      } // loop over volume element neighbors

      // check, if element could be found
      if( !ptVolElem ) {
        EXCEPTION("Could not find a suitable volume neighbor for surface element #"
            << surfElem.elemNum << ". " );
      }
      // create new local point for volume element 
      lpmVol.reset(new LocPointMapped());
      esmVol = 
          shapeMap->GetGrid()->GetElemShapeMap( ptVolElem, 
                                                shapeMap->IsUpadet() );
    } else {
      esmVol = lpmVol->shapeMap;
    } // elemIsSame

    // calculate local point in volume element, corresponding to
    // the surface element point
    //perhaps here we can distinguish between the NC_SURF_ELEM case
    //in which the local coordinates are already available

    LocPoint lpVol;
    Vector<Double> locNormal;
    esmVol->GetLocalIntPoints4Surface(  ptEl->connect, lp, 
                                        lpVol,locNormal);
    lpmVol->Set( lpVol, esmVol);

    // calculate global normal pointing into current volume element
    normal = Transpose(lpmVol->jacInv) * locNormal;
    normal /= normal.NormL2();

    // multiply normal by normal sign, so it points out of the first
    // volume element
    normal *= normalFactor;
  }

// ========================================================================
//  ElemShapeMap
// ========================================================================

// Definition of finite element type mappings
static EnumTuple elemShapeTuples[] = 
{
 EnumTuple(ElemShapeMap::NO_TYPE,  "NO_TYPE"), 
 EnumTuple(ElemShapeMap::LAGRANGE,  "LAGRANGE"),
 EnumTuple(ElemShapeMap::LAGRANGE_BLENDED,  "LAGRANGE_BLENDED"),
 EnumTuple(ElemShapeMap::ANALYTICAL,  "ANALYTICAL")
};

Enum<ElemShapeMap::ShapeMapType> ElemShapeMap::shapeMapType = \
    Enum<ElemShapeMap::ShapeMapType>("Finite Element Shape Mapping Types",
        sizeof(elemShapeTuples) / sizeof(EnumTuple),
        elemShapeTuples);
   
ElemShapeMap::ElemShapeMap( Grid* ptGrid ) {
  type_ = NO_TYPE;
  ptGrid_ = ptGrid;
  depth_ = 1.0;
  isAxi_ = false;
  ptElem_ = NULL;
  ptSurfElem_ = NULL;
}

ElemShapeMap::~ElemShapeMap() {
  
}

void ElemShapeMap::SetElem( const Elem* ptElem, bool isUpdated ) {
  ptElem_ = ptElem;
  isUpdated_ = isUpdated;
  isAxi_ = ptGrid_->IsAxi();
  
  // Check, if the element is a surface element
  if( Elem::shapes[ptElem->type].dim == ptGrid_->GetDim()-1 ) {
    try {
      ptSurfElem_ = dynamic_cast<const SurfElem*>(ptElem);
    } catch(...) {
      EXCEPTION( "Could not convert element #" << ptElem->elemNum
                 << " to a surface element" );
    }
  } else {
    ptSurfElem_ = NULL;
  }
}


// ========================================================================
//  Lagrangian Element Shape Map
// ========================================================================

// declaration of static member data for geometric reference elements
std::map<Elem::FEType, FeH1LagrangeExpl* > LagrangeElemShapeMap::feMap_;

void LagrangeElemShapeMap::InitStaticMembers() {
   
  feMap_ = 
      boost::assign::list_of< std::pair<Elem::FEType, FeH1LagrangeExpl* > > 
  (Elem::ET_LINE2, new FeH1LagrangeLine1()) 
  (Elem::ET_LINE3, new FeH1LagrangeLine2()) 
  (Elem::ET_TRIA3, new FeH1LagrangeTria1()) 
  (Elem::ET_TRIA6, new FeH1LagrangeTria2()) 
  (Elem::ET_QUAD4, new FeH1LagrangeQuad1()) 
  (Elem::ET_QUAD8, new FeH1LagrangeQuad2()) 
  (Elem::ET_QUAD9, new FeH1LagrangeQuad9())
  (Elem::ET_TET4,  new FeH1LagrangeTet1()) 
  (Elem::ET_TET10, new FeH1LagrangeTet2()) 
  (Elem::ET_HEXA8, new FeH1LagrangeHex1()) 
  (Elem::ET_HEXA20, new FeH1LagrangeHex2()) 
  (Elem::ET_HEXA27, new FeH1LagrangeHex27()) 
  (Elem::ET_WEDGE6, new FeH1LagrangeWedge1()) 
  (Elem::ET_WEDGE15, new FeH1LagrangeWedge2())
  (Elem::ET_WEDGE18, new FeH1LagrangeWedge18()) 
  (Elem::ET_PYRA5, new FeH1LagrangePyra1()) 
  (Elem::ET_PYRA13, new FeH1LagrangePyra2())
  (Elem::ET_PYRA14, new FeH1LagrangePyra14());
}



LagrangeElemShapeMap::LagrangeElemShapeMap( Grid* ptGrid  ) 
: ElemShapeMap(ptGrid ) {
  type_ = LAGRANGE; 
  intScheme_ = ptGrid_->GetIntegrationScheme();
  
//  // Fill map with reference elements
//  feMap_[Elem::ET_LINE2] = new FeH1LagrangeLine1(); 
//  feMap_[Elem::ET_LINE3] = new FeH1LagrangeLine2();
//  feMap_[Elem::ET_TRIA3] = new FeH1LagrangeTria1();
//  feMap_[Elem::ET_TRIA6] = new FeH1LagrangeTria2();
//  feMap_[Elem::ET_QUAD4] = new FeH1LagrangeQuad1();
//  feMap_[Elem::ET_QUAD8] = new FeH1LagrangeQuad2();
//  feMap_[Elem::ET_QUAD9] = new FeH1LagrangeQuad9();
//  feMap_[Elem::ET_TET4]  = new FeH1LagrangeTet1();
//  feMap_[Elem::ET_TET10] = new FeH1LagrangeTet2();
//  feMap_[Elem::ET_HEXA8] = new FeH1LagrangeHex1();
//  feMap_[Elem::ET_HEXA20] = new FeH1LagrangeHex2();
//  feMap_[Elem::ET_HEXA27] = new FeH1LagrangeHex27();
//  feMap_[Elem::ET_WEDGE6] = new FeH1LagrangeWedge1();
//  feMap_[Elem::ET_WEDGE15] = new FeH1LagrangeWedge2();
//  feMap_[Elem::ET_WEDGE18] = new FeH1LagrangeWedge18();
//  feMap_[Elem::ET_PYRA5] = new FeH1LagrangePyra1();
//  feMap_[Elem::ET_PYRA13] = new FeH1LagrangePyra2();
//  feMap_[Elem::ET_PYRA14] = new FeH1LagrangePyra14();
}

LagrangeElemShapeMap::~LagrangeElemShapeMap() {
  
//  // Remove pointers to all reference elements
//  std::map<Elem::FEType, FeH1LagrangeExpl *>::iterator it = feMap_.begin();
//  for( ; it != feMap_.end(); it++ ) {
//    delete  it->second;
//  }
//  feMap_.clear();
}

void LagrangeElemShapeMap::Local2Global( Vector<Double>& globPoint, 
                                         const LocPoint& lp ) {
  
  //step 1: evaluate shape fncs. at local coordinate
  Vector<Double> shFnc;
  ptFe_->GetShFnc( shFnc, lp, NULL, 0 );

  // step2: multiply shape fncs for each dimension with according matrix entries
  globPoint = coords_ * shFnc;
}

   
void LagrangeElemShapeMap::Global2Local( Vector<Double>& locPoint, 
                                         const Vector<Double>& globalPoint ) {


    Elem::FEType curType = ptFe_->FeType();
    switch(curType){
    case Elem::ET_QUAD4:
      //Global2LocalGeneral(locPoint,globalPoint);
      Global2LocalQuad4(locPoint,globalPoint);
      break;
    default:
      Global2LocalGeneral(locPoint,globalPoint);
      break;
    }
  }

void LagrangeElemShapeMap::Global2LocalQuad4( Vector<Double>& locPoint,
                                         const Vector<Double>& globalPoint ){
  //based on FEM skript found by simon at
  //http://www.colorado.edu/engineering/cas/courses.d/IFEM.d/IFEM.Ch23.d/IFEM.Ch23.pdf

  //corner coords
  Double x1=0.0,x2=0.0,x3=0.0,x4=0.0;
  Double y1=0.0,y2=0.0,y3=0.0,y4=0.0;
  //global point
  Double x_p,y_p;
  Double x_b,y_b;
  Double x_cx,y_cx,x_ce,y_ce,A,J_1,J_2,x_0,y_0,x_p0,y_p0;
  Double b_xi,b_eta,c_xi,c_eta;

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
  x_b = x1-x2+x3-x4;
  y_b = y1-y2+y3-y4;

  x_cx = x1+x2-x3-x4;
  y_cx = y1+y2-y3-y4;
  x_ce = x1-x2-x3+x4;
  y_ce = y1-y2-y3+y4;

  A = 0.5 * ( (x3-x1)*(y4-y2) - (x4-x2)*(y3-y1) );
  J_1 = (x3-x4)*(y1-y2) - (x1-x2)*(y3-y4);
  J_2 = (x2-x3)*(y1-y4) - (x1-x4)*(y2-y3);

  x_0 = 0.25 * (x1+x2+x3+x4);
  y_0 = 0.25 * (y1+y2+y3+y4);

  x_p0 = x_p - x_0;
  y_p0 = y_p - y_0;

  b_xi =  A - (x_p0*y_b) + (y_p0*x_b); b_eta = (-1.0 * A) - (x_p0*y_b) + (y_p0*x_b);

  c_xi = (x_p0*y_cx) - (y_p0*x_cx); c_eta = (x_p0*y_ce) - (y_p0*x_ce);

  locPoint[0] = 2.0*c_xi / ((-1.0 * std::sqrt(b_xi*b_xi - 2.0 * J_1 * c_xi))  - b_xi);
  locPoint[1] = 2.0*c_eta / ( std::sqrt(b_eta*b_eta + (2.0 * J_2 * c_eta)) - b_eta);
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


void LagrangeElemShapeMap::Global2LocalGeneral( Vector<Double>& locPoint,
                                         const Vector<Double>& globalPoint ){

  // increate counter
  UInt globDim = globalPoint.GetSize(); // determine global dimension
  UInt locDim = shape_.dim; // dimension of current element

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
  xi_start.Resize(3); xi_start.Init();
  xi_k.Resize(3); xi_k.Init();
  delta_xi.Resize(3);  delta_xi.Init();
  J.Resize(3, 3); J.Init();

  // Perform Newton-Raphson method on global point
  // Find good startpoint xi_k among local node coordinates
  f_min = 999e5; // really big value!
  for(UInt k = 0; k < shape_.numNodes; k++)
  {
    for(UInt l = 0; l < locDim; l++) {
      xi_k[l] = shape_.nodeCoords[k][l];
    }

    Local2Global(f, xi_k);
    f = f - globalPoint;
    distance = f.NormL2();
    if( distance < f_min ) {
      f_min = distance;
      xi_start = xi_k;
      f_start = f;
      minDist = distance;
    }
  }
  xi_k = xi_start;
  f = f_start;
  distance = f_min;
  // divergence = false;
  converged = false;

  do {
    // Calculate Jacobian at iteration point xi_k
    CalcJ(J, xi_k);

    // locDim should never be 1 since line elements will
    // be handled differently
    if(locDim == 1)
    {
      EXCEPTION("Line elements should not use the Newton-Raphson method!");
      return;
    }

    if(globDim == 2)
    {
      // Find new local search direction for 2D -> 2D mapping using
      // Cramer's rule.

      // Jacobian Matrix:
      //
      //  J = ( J_00  J_01 )
      //      ( J_10  J_11 )

      jacDet = + J[0][0]*J[1][1] - J[1][0]*J[0][1];

      //  ( f_0  J_01 )
      //  ( f_1  J_11 )
      delta_xi[0] = - J[1][1]*f[0] + J[0][1]*f[1];

      //  ( J_00  f_0 )
      //  ( J_10  f_1 )
      delta_xi[1] = - J[0][0]*f[1] + J[1][0]*f[0];

      distNormalizer = sqrt(fabs(jacDet));
    }
    else
    {
      if(locDim == 2)
      {
        // Find new local search direction for 2D -> 3D mapping.
        // Project 3D difference vector onto 2D basis given by
        // the Jacobian to find the new 2D search direction.

        // Jacobian Matrix:
        //
        //      ( J_00  J_01  1 )
        //  J = ( J_10  J_11  1 )
        //      ( J_20  J_21  1 )

        jacDet = + J[1][0] * J[2][1] - J[2][0] * J[1][1]
                 - J[0][0] * J[2][1] + J[2][0] * J[0][1]
                 + J[0][0] * J[1][1] - J[1][0] * J[0][1];

        // Calculate negative Jacobian determinants of the following matrices
        // to find the local search direction which has to point in the opposite
        // direction of the backprojected global error vector f.
        //
        //  ( f_0  J_01  1 )
        //  ( f_1  J_11  1 )
        //  ( f_2  J_21  1 )

        delta_xi[0] = - f[0] * J[1][1] - f[1] * J[2][1] - f[2] * J[0][1]
                      + f[2] * J[1][1] + f[1] * J[0][1] + f[0] * J[2][1];

        //  ( J_00  f_0  1 )
        //  ( J_10  f_1  1 )
        //  ( J_20  f_2  1 )

        delta_xi[1] = - f[0] * J[2][0] - f[1] * J[0][0] - f[2] * J[1][0]
                      + f[2] * J[0][0] + f[1] * J[2][0] + f[0] * J[1][0];

        distNormalizer = sqrt(fabs(jacDet));
      }
      else
      {
        // Find new local search direction for 3D -> 3D mapping using
        // Cramer's rule.

        //      ( J_00  J_01  J_02 )
        //  J = ( J_10  J_11  J_12 )
        //      ( J_20  J_21  J_22 )

        jacDet = + J[0][0]*J[1][1]*J[2][2] - J[0][0]*J[2][1]*J[1][2]
                 - J[1][0]*J[0][1]*J[2][2] + J[2][1]*J[1][0]*J[0][2]
                                                                                                                   - J[1][1]*J[2][0]*J[0][2] + J[2][0]*J[0][1]*J[1][2];

        //  ( f_0  J_01  J_02 )
        //  ( f_1  J_11  J_12 )
        //  ( f_2  J_21  J_22 )
        delta_xi[0] = - J[1][1]*J[2][2]*f[0] + J[2][1]*J[1][2]*f[0]
                      - J[2][1]*J[0][2]*f[1] + J[0][1]*J[2][2]*f[1]
                      - J[0][1]*J[1][2]*f[2] + J[1][1]*J[0][2]*f[2];

        //  ( J_00  f_0  J_02 )
        //  ( J_10  f_1  J_12 )
        //  ( J_20  f_2  J_22 )
        delta_xi[1] = - J[1][2]*J[2][0]*f[0] + J[1][0]*J[2][2]*f[0]
                      - J[0][0]*J[2][2]*f[1] + J[2][0]*J[0][2]*f[1]
                      - J[1][0]*J[0][2]*f[2] + J[1][2]*J[0][0]*f[2];

        //  ( J_00  J_01  f_0 )
        //  ( J_10  J_11  f_1 )
        //  ( J_20  J_21  f_2 )
        delta_xi[2] = - J[1][0]*J[2][1]*f[0] + J[2][0]*J[1][1]*f[0]
                      - J[0][1]*J[2][0]*f[1] + J[2][1]*J[0][0]*f[1]
                      - J[0][0]*J[1][1]*f[2] + J[1][0]*J[0][1]*f[2];

        distNormalizer = std::pow(fabs(jacDet), 1.0 / 3.0);
      }
    }

    // Here is the new local search direction. We normalize it so we
    // can be sure, that we have a local search vector of a length
    // comparable to the local element diameter.
    Double len;

    delta_xi[0] /= jacDet;
    delta_xi[1] /= jacDet;
    delta_xi[2] /= jacDet;

    len = delta_xi.NormL2();
    delta_xi[0] /= len;
    delta_xi[1] /= len;
    delta_xi[2] /= len;

    // If global element is smaller use local distance as a measure.
    // If global element is bigger use global distance as a measure.
    distMeasure = distNormalizer < 1 ? distance/distNormalizer : distance;
    if(distMeasure < TOL)
    {
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

    for(l=0; l<20; l++)
    {
      xi_k = xi_start + delta_xi * interval[1];
      Local2Global(f2, xi_k);
      f2 = f2 - globalPoint;
      dist[1] = f2.NormL2();

      distMeasure = distNormalizer < 1 ? dist[1]/distNormalizer : dist[1];
      if(distMeasure < TOL)
      {
        converged = true;
        break;
      }

      if(f.Inner(f2) < 0) 
      {
        break;
      } else 
      {
        interval[0] = interval[1];
        interval[1] *= 4.0;
        dist[0] = dist[1];
        f = f2;
      }

      if(dist[0] < minDist) {
        minDist = dist[0];
        xi_min = xi_start + delta_xi * interval[0];
      }

      if(dist[1] < minDist) {
        minDist = dist[1];
        xi_min = xi_start + delta_xi * interval[1];
      }
    }

    if(converged)
      break;

    // Try to narrow down the extents of the interval by searching points
    // ever closer from below and above the minumum in the current search
    // direction.
    for(l=0; l<20; l++)
    {
      Double x3 = interval[0] + (interval[1] - interval[0]) * golden_ratio;
      Double x4 = interval[1] - (interval[1] - interval[0]) * golden_ratio;

      distMeasure = distNormalizer < 1 ? dist[0]/distNormalizer : distance;
      if(distMeasure < TOL)
      {
        converged = true;
        break;
      }

      distMeasure = distNormalizer < 1 ? dist[1]/distNormalizer : distance;
      if(distMeasure < TOL)
      {
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

      if(f1.Inner(f2) > 0) 
      {
        // Either x3 switched to the side of x4 in respect to the searched
        // point or x4 switched to the side of x3.
        UInt m;
        const UInt n=16;

        for(m = 0; m < n && f1.Inner(f2) > 0; m++)
        {
          x3 -= (interval[1] - interval[0]) * (golden_ratio/n);

          xi_k = xi_start + delta_xi * x3;
          Local2Global(f1, xi_k);

          f1 = f1 - globalPoint;
        }

        if(m==n) 
        {              
          for(m = 0; m < n && f1.Inner(f2) > 0; m++)
          {
            x4 += (interval[1] - interval[0]) * (golden_ratio/n);

            xi_k = xi_start + delta_xi * x4;
            Local2Global(f2, xi_k);

            f2 = f2 - globalPoint;
          }
        } else 
        {
          for(m = 0; m < n && f1.Inner(f2) < 0; m++)
          {
            x4 -= (interval[1] - interval[0]) * (golden_ratio/n);

            xi_k = xi_start + delta_xi * x4;
            Local2Global(f2, xi_k);

            f2 = f2 - globalPoint;
          }
        }

        if(x4 < x3)
          x4 += (interval[1] - interval[0]) * (golden_ratio/n);
      }


      // Update interval and distance array with new values.
      distance_l = f1.NormL2();

      if(distance_l < dist[0]) 
      {
        interval[0] = x3;
        dist[0] = distance_l;
      }


      distance_l = f2.NormL2();

      if(distance_l < dist[1]) 
      {
        interval[1] = x4;
        dist[1] = distance_l;
      }

      if(dist[0] < minDist) {
        minDist = dist[0];
        xi_min = xi_start + delta_xi * interval[0];
      }

      if(dist[1] < minDist) {
        minDist = dist[1];
        xi_min = xi_start + delta_xi * interval[1];
      }
    }

    // If distances increase again, break
    if(dist[0] > minDist && dist[1] > minDist)
    {
      xi_k = xi_min;
      break;
    }

    if(dist[0] < dist[1])
      xi_k = xi_start + delta_xi * interval[0];
    else
      xi_k = xi_start + delta_xi * interval[1];

    distance = minDist;
    xi_start = xi_k;
    Local2Global(f, xi_k);
    f = f - globalPoint;

    if(converged)
      break;

    iter++;
  }
  while(iter < 6);

  Vector<Double> globMapped, errVec;
  Local2Global(globMapped, xi_k);
  errVec = globMapped - globalPoint;
  locPoint = xi_k;
}

   
void LagrangeElemShapeMap::GetGlobMidPoint( Vector<Double>& midPoint ) {
  Local2Global(midPoint, shape_.midPointCoord);
}

Double LagrangeElemShapeMap::CalcVolume( ) {
  
  // Get integration points
  StdVector<LocPoint> intPoints;
  StdVector<Double> weights;
  
  // Order: use element order and add 2 to be sure for curved elements
  UInt order = shape_.order + 2;
  intScheme_->GetIntPoints( Elem::GetShapeType(ptElem_->type), 
                            IntScheme::GAUSS, order, 
                            intPoints, weights );
  Double vol = 0.0;
  Double jacDet = 0.0;
  Matrix<Double> jac;
  // Loop over all integration points
  LocPointMapped lpm;
  for( UInt i = 0; i < intPoints.GetSize(); i++  ) {
    jacDet = CalcJDet(jac, intPoints[i] );
    vol +=  jacDet * weights[i];
  }
  return  vol;
}


void LagrangeElemShapeMap::CalcNormal( Vector<Double>& normal, 
                                       const LocPoint& lp ) {
  
  
  // check, that element is a surface element at all
  if( shape_.dim != ptGrid_->GetDim() - 1 ) {
    EXCEPTION("Can not calculate normal of element #"
              << ptElem_->elemNum << " which is of dimension "
              <<  shape_.dim << " in a " 
              << ptGrid_->GetDim() << "-dimensional grid!");
  }
  
  // get neighboring volume element
  const SurfElem & surfEl = *ptSurfElem_;
  Elem * ptVolEl = surfEl.ptVolElems[0];
  
  
  // Obtain shape map of neighboring volume element
  LagrangeElemShapeMap sm(ptGrid_);
  sm.SetElem(ptVolEl, isUpdated_);
  
  // Map local point of surface to global point 
  Vector<Double> locNormal;
  LocPoint volPoint;
  sm.ptFe_->GetLocalIntPoints4Surface( ptElem_->connect, 
                                       ptVolEl->connect,
                                       lp, volPoint, locNormal);
//  std::cerr << "\nlp = " << lp.coord.ToString() << std::endl;
//  std::cerr << "volPoint = " << volPoint.coord.ToString() << std::endl;
//  std::cerr << "locNormal = " << locNormal.ToString() << std::endl;
  
  // Calculate Jacobian matrix of volume element in that point
  Matrix<Double> jac;
  sm.CalcJ( jac, volPoint);
  
  // Calculate global normal
  Matrix<Double> jInv;
  jac.Invert(jInv);
  normal = Transpose(jInv) * locNormal;

  // normalize normal
  Double norm = normal.NormL2();
  normal /= norm;
}

void  LagrangeElemShapeMap::
CalcNormalOutOfVol( Vector<Double> & normal,
                    const LocPoint& lp, 
                    const Elem & volElem ) {

  // Obtain shape map of neighboring volume element
  LagrangeElemShapeMap sm(ptGrid_);
  sm.SetElem(&volElem, isUpdated_);
  
  Matrix<Double> & volCoords = sm.coords_; 

  // Calculate surface normal without defined sign
  CalcNormal( normal, lp);
  
  UInt volCorners = sm.shape_.numVertices;

  /* Idea:
       - find a common surface / volume element node
       - find additional node, which is not common (i.e. lies on the "volume side"
       - calculate difference vector from  "additional" node to first common
         one ( = pointing out of  the volume element)
       - the scalar product of the normal and the difference vector determines
         the normalSign
   */

  // find first common vertex index
  Integer firstCommonIndex = -1;
  for (UInt i=0; i<volCorners; i++)
    if (volElem.connect[i] == ptElem_->connect[0]){
      firstCommonIndex = i;
      break;
    }

  // find additional node of the volume node, which is not contained in the
  // surface element
  std::set<UInt> volSet, surfSet, diffSet;
  volSet = std::set<UInt>(sm.ptElem_->connect.Begin(), 
                          sm.ptElem_->connect.End());
  
  surfSet = std::set<UInt>(ptElem_->connect.Begin(), 
                           ptElem_->connect.End());
  std::set_difference( volSet.begin(), volSet.end(),
                       surfSet.begin(), surfSet.end(),
                       std::inserter(diffSet,diffSet.begin()) );

  Vector<double> diffVec( sm.shape_.dim );
  Double scalarProd = 0.0;
  std::set<UInt>::const_iterator it = diffSet.begin();
  for( ; it != diffSet.end(); it++ ) {

    UInt node = *it;

    // search for volume 
    UInt index = volElem.connect.Find(node);

    // calculate difference vector (pointing out of the volume)
    for( UInt iDim = 0; iDim < sm.shape_.dim; ++iDim ) {
      diffVec[iDim] =  volCoords[iDim][firstCommonIndex]
                     - volCoords[iDim][index];
    }
    // normalize difference vector to 1.0
    diffVec /= diffVec.NormL2();

    // calculate scalar product
    scalarProd = diffVec * normal;

    // check if scalar product is != 0 (otherwise we may have a degenerated
    // element, where two nodes lie on the same location)
    if( std::abs(scalarProd) < EPS ) {
      // we have to find another node
      continue;
    } else {
      if( scalarProd < 0.0 ) {
        // original normal vector points into the volume -> reorient
        normal *= -1.0;
      }
      break;
    }
  }
}

void LagrangeElemShapeMap::
GetLocalIntPoints4Surface( const StdVector<UInt> & surfConnect,
                           const LocPoint & surfIntPoint,
                           LocPoint & volIntPoint,
                           Vector<Double>& locNormal ) {
  ptFe_->GetLocalIntPoints4Surface(  surfConnect, ptElem_->connect,
                                    surfIntPoint, volIntPoint, locNormal );
}

bool LagrangeElemShapeMap::
 CoordIsInsideElem( const Vector<Double>& point, Double tolerance ) {
  
  return ptFe_->CoordIsInsideElem(point, tolerance);
}

void LagrangeElemShapeMap::CalcDiameter( Vector<Double>& diameter ) {
  Vector<Double> mins(shape_.dim), maxs(shape_.dim);
  mins.Init(std::numeric_limits<double>::max());
  maxs.Init();
  diameter.Resize(shape_.dim);
  diameter.Init();
  for (UInt dim = 0; dim < shape_.dim; dim++)
  {
    for (UInt k = 0, n_elems = coords_.GetNumCols(); k < n_elems; k++)
    {
      Double test = coords_[dim][k];
      Double& min = mins[dim];
      Double& max = maxs[dim];
      min = std::min(min, test);
      max = std::max(max, test);
    }

    diameter[dim] = maxs[dim] - mins[dim];
  }
}

void LagrangeElemShapeMap::CalcBarycenter( Vector<Double>& baryCenter ) {
  EXCEPTION("Not implemented");
  //    UInt n_dims  = coords.GetNumRows();
  //    UInt n_elems = coords.GetNumCols();
  //
  //    // init barycenter for safty reason
  //    barycenter.SetZero();
  //
  //    // std::cout << "calc a new barycenter" << std::endl;
  //    // a barycenter is simply the average of all coordinates
  //    for (UInt dim=0; dim < n_dims; dim++)
  //    {
  //      // std::cout << "dim = " << dim << "  ";
  //      for (UInt k=0; k < n_elems; k++)
  //      {
  //        // the constructor of Point initializes
  //        barycenter[dim] += coords[dim][k];
  //        // std::cout << coords[dim][k] << "->" << barycenter[dim] << "\t";
  //      }
  //
  //      barycenter[dim] /= (double) n_elems;
  //      // std::cout << " average: " << (barycenter[dim]) << std::endl;
  //    }
}

void LagrangeElemShapeMap::GetMaxMinEdgeLength(Double& max, Double& min ) {
  max = 1e-100;
  min = 1e+100;
  Double length, dl;
  for( UInt i = 0; i < shape_.numEdges; ++i  ) {
    length = 0.0;
    for( UInt iDim = 0; iDim < shape_.dim; ++iDim ) {
      dl = coords_[iDim][shape_.edgeVertices[i][1]-1]
          -coords_[iDim][shape_.edgeVertices[i][0]-1];
     length += dl*dl;
    }
    length = sqrt(length);
    max = max > length ? max : length;
    min = min < length ? min : length;
  }
}

void LagrangeElemShapeMap::GetEdgeLength( StdVector<Double>& edges_out) {
  EXCEPTION("Not implemented");
}

void LagrangeElemShapeMap::CalcJ( Matrix<Double>& jac, 
                                  const LocPoint& lp ) {
  Matrix<Double> deriv;
  ptFe_->GetLocDerivShFnc( deriv, lp, ptElem_ );
  jac = coords_ * deriv;
  jac *= depth_; // explicitly include depth_ of setup
}

Double LagrangeElemShapeMap::CalcJDet( Matrix<Double>& jac, 
                                       const LocPoint& lp ) {
  Matrix<Double> deriv;
  
  ptFe_->GetLocDerivShFnc( deriv, lp, ptElem_ );
  jac = coords_ * deriv;
  jac *= depth_; // explicitly include depth_ of setup
  
  Double jacDet = 0.0;
  // redundant code, but necessary at some points
  if( jac.GetNumCols() == jac.GetNumRows() ) {
    jac.Determinant(jacDet);
  } else if ( jac.GetNumRows() == 3) {
    // 2D elements in 3D
    Vector<Double> normal; 
    normal.Resize(3);
    normal[0]= jac[1][0]* jac[2][1]- jac[2][0]* jac[1][1];
    normal[1]=jac[2][0]*jac[0][1]- jac[0][0]* jac[2][1];
    normal[2]= jac[0][0]* jac[1][1]- jac[1][0]*jac[0][1];
    jacDet = sqrt(normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]);

  } else if ( jac.GetNumRows() == 2) {
    // 1D elements in 2D
    //see kaltenbacher, p.23, eq.(2.122)
    jacDet = sqrt(jac[0][0]*jac[0][0] + jac[1][0]*jac[1][0]);
  }
  
  return jacDet;
}
   
   
void LagrangeElemShapeMap::SetElem( const Elem* ptElem, bool isUpdated ) {
  
  ElemShapeMap::SetElem( ptElem, isUpdated);
//  ptElem_ = ptElem;
//  isUpdated_ = isUpdated;
//  isAxi_ = ptGrid_->IsAxi();
//  
  // call setElem at base class
  
  
  // get coordinates from grid
  ptGrid_->GetElemNodesCoord( coords_,ptElem->connect, isUpdated_ );

  // set reference element
#ifndef NDEBUG
  if( feMap_.find(ptElem->type) == feMap_.end()) {
    EXCEPTION("Element of type '" << Elem::feType.ToString(ptElem->type)
              << "' not defined for Lagrangian Shape Map!");
  }
#endif
  ptFe_ = feMap_[ptElem->type];
  shape_ = Elem::shapes[ptElem_->type];
}

} // namespace CoupledField
