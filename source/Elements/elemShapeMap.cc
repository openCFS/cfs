#include "elemShapeMap.hh"
#include "Domain/grid.hh"
#include "H1ElemsLagExpl.hh"

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

// ===========================================================================
//  C L A S S   LocPointMapped
// ===========================================================================
  LocPointMapped::LocPointMapped()
  : jacDet( 0.0 ) {
    
  }
  
  void LocPointMapped::Set( const LocPoint& lp, shared_ptr<ElemShapeMap> esm ) {

    shapeMap = esm;
    this->lp = lp;

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
    }

    // Check, if geometry is axi-symmetric. In this case scale the 
    // Jacobian determinant with 2*pi*r
    Vector<Double> globPoint;
    if( esm->IsAxi() ) {
      esm->Local2Global( globPoint, lp);
      jacDet *= 2 * PI * globPoint[0];
    }
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
   
ElemShapeMap::ElemShapeMap( Grid* ptGrid, bool isUpdated ) {
  type_ = NO_TYPE;
  ptGrid_ = ptGrid;
  isAxi_ = ptGrid->IsAxi();
  isUpdated_ = isUpdated;
  depth_ = 1.0;
}

ElemShapeMap::~ElemShapeMap() {
  
}

void ElemShapeMap::SetElem( const Elem* ptElem ) {
  ptElem_ = ptElem;
}


// ========================================================================
//  Lagrangian Element Shape Map
// ========================================================================

LagrangeElemShapeMap::LagrangeElemShapeMap( Grid* ptGrid, bool isUpdated ) 
: ElemShapeMap(ptGrid, isUpdated ) {
  type_ = LAGRANGE; 
  
  
  // Fill map with reference elements
  feMap_[Elem::ET_LINE2] = new FeH1LagrangeLine1(); 
  feMap_[Elem::ET_LINE3] = new FeH1LagrangeLine2();
//  feMap_[ET_TRIA3] = new FeH1Lagrange();
//  feMap_[ET_TRIA6] = new FeH1LagrangeLine1();
  feMap_[Elem::ET_QUAD4] = new FeH1LagrangeQuad1();
  feMap_[Elem::ET_QUAD8] = new FeH1LagrangeQuad2();
//  feMap_[ET_QUAD9] = new FeH1LagrangeQuad1();
//  feMap_[ET_TET4] = new FeH1LagrangeLine1();
//  feMap_[ET_TET10] = new FeH1LagrangeLine1();
  feMap_[Elem::ET_HEXA8] = new FeH1LagrangeHex1();
//  feMap_[ET_HEXA20] = new FeH1LagrangeHexa2();
//  feMap_[ET_HEXA27] = new FeH1LagrangeLine1();
//  feMap_[ET_PYRA5] = new FeH1LagrangeLine1();
//  feMap_[ ET_PYRA13] = new FeH1LagrangeLine1();
//  feMap_[ET_WEDGE6] = new FeH1LagrangeLine1();
//  feMap_[ET_WEDGE15] = new FeH1LagrangeLine1();
}

LagrangeElemShapeMap::~LagrangeElemShapeMap() {
  
  // Remove pointers to all refrence elements
  std::map<Elem::FEType, FeH1 *>::iterator it = feMap_.begin();
  for( ; it != feMap_.end(); it++ ) {
    delete  it->second;
  }
  feMap_.clear();
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
                                               const Vector<Double>& glob ) {
  EXCEPTION("Not implemented");
  // Following: Old code from BaseFE
    // ===============================
    //    Vector<Double> globalPoint; // global point coordinates
    //    UInt globDim = globalCoords.GetNumRows(); // determine global dimension
    //    UInt numPoints = globalCoords.GetNumCols(); // number of global points
    //    UInt locDim = LCornerCoords_.GetNumRows(); // dimension of current element
    //    //    UInt numCorners = LCornerCoords_.GetNumCols(); // number of element corners
    //
    //    Vector<Double> xi_start; // local start point for Newton-Raphson method
    //    Vector<Double> xi_k; // local point at iteration k
    //    Vector<Double> delta_xi; // local search direction
    //    Vector<Double> f; // global (negative) search direction
    //    Vector<Double> f_start; // start value for global search direction
    //    Matrix<Double> J; // Jacobian at local point xi_k
    //    Double f_min; // minimal distance between global corners and global point
    //    Double distance; // global distance!
    //    Double distance_l; // global distance at iteration l
    //    Double distMeasure; // global distance or local distance
    //                        // depending on sqrt(|jacDet|)
    //    Double TOL = 0.0000001; // tolerance for global distance
    //    Double jacDet; // denominator for Cramer's rule.
    //    Double distNormalizer; // denominator for Cramer's rule.
    //    Double coeff; // damping coefficient
    //    bool divergence; // does the Newton-Raphson algorithm diverge?
    //    bool converged; // have we found the local point?
    //    UInt iter = 0;
    //
    //    // Initialize variables
    //    globalPoint.Resize(globDim);
    //    xi_start.Resize(3); xi_start.Init();
    //    xi_k.Resize(3); xi_k.Init();
    //    delta_xi.Resize(3); ; delta_xi.Init();
    //    localCoords.Resize(locDim, numPoints);
    //    J.Resize(3, 3); J.Init();
    //
    //    // Perform Newton-Raphson method on the list of global points
    //    // to find local coordinates within this element.
    //    for(UInt i = 0; i < numPoints; i++)
    //    {
    //      // Copy global point coordinates into a Vector for algebraic ops.
    //      for(UInt j = 0; j < globDim; j++)
    //      {
    //        globalPoint[j] = globalCoords[j][i];
    //      }
    //
    //      // Find good startpoint xi_k among local corner coordinates
    //      f_min = 999e5; // really big value!
    //      for(UInt k = 0; k < NumCorners_; k++)
    //      {
    //        for(UInt l = 0; l < locDim; l++)
    //        {
    //          xi_k[l] = LCornerCoords_[l][k];
    //        }
    //
    //        Local2GlobalCoord(f, xi_k, coordMat, NULL);
    //        f = f - globalPoint;
    //        distance = f.NormL2();
    //        if( distance < f_min )
    //        {
    //          f_min = distance;
    //          xi_start = xi_k;
    //          f_start = f;
    //        }
    //      }
    //
    //      /*
    //      std::cout << "Beginning to find local coords..." << std::endl;
    //      std::cout << "Global coords..." << globalPoint << std::endl;
    //      */
    //      xi_k = xi_start;
    //      f = f_start;
    //      distance = f_min;
    //      divergence = false;
    //      converged = false;
    //
    //      do
    //      {
    //        // Calculate Jacobian at iteration point xi_k
    //        CalcJacobian(J, xi_k, coordMat, NULL );
    //
    //        /*
    //        std::cout << std::endl << "++++++ xi_k (" << xi_k[0] << ", "
    //                  << xi_k[1] << ", " << xi_k[2] << ") ++++++" << std::endl;
    //        std::cout << "f (" << f[0] << ", " << f[1]
    //                  << ", " << f[2] << ")" << std::endl;
    //        std::cout << "distance " << distance << std::endl;
    //        */
    //
    //        // locDim should never be 1 since line elements will
    //        // be handled differently
    //        if(locDim == 1)
    //        {
    //          EXCEPTION("Line elements should not use the Newton-Raphson method!");
    //          return;
    //        }
    //
    //        if(globDim == 2)
    //        {
    //          // Find new local search direction for 2D -> 2D mapping using
    //          // Cramer's rule.
    //          jacDet = + J[0][0]*J[1][1] - J[1][0]*J[0][1];
    //
    //          delta_xi[0] = - J[1][1]*f[0] + J[0][1]*f[1];
    //
    //          delta_xi[1] = - J[0][0]*f[1] + J[1][0]*f[0];
    //
    //          distNormalizer = sqrt(fabs(jacDet));
    //        }
    //        else
    //        {
    //          if(locDim == 2)
    //          {
    //            // Find new local search direction for 2D -> 3D mapping.
    //            // Project 3D difference vector onto 2D basis given by
    //            // the Jacobian to find the new 2D search direction.
    //            //            Vector<Double> normal;
    //
    //            jacDet = + J[1][0] * J[2][1] - J[2][0] * J[1][1]
    //                     - J[0][0] * J[2][1] + J[2][0] * J[0][1]
    //                     + J[0][0] * J[1][1] - J[1][0] * J[0][1];
    //            jacDet = -fabs(jacDet);
    //            
    //            delta_xi[0] = (f[0] * J[0][0] +
    //                           f[1] * J[1][0] +
    //                           f[2] * J[2][0]);
    //
    //            delta_xi[1] = (f[0] * J[0][1] +
    //                           f[1] * J[1][1] +
    //                           f[2] * J[2][1]);
    //
    //            distNormalizer = sqrt(fabs(jacDet));
    //
    //            /*
    //            std::cout << "J1 (" << J[0][0] << ", "
    //                      << J[1][0] << ", " << J[2][0] << ")" << std::endl;
    //            std::cout << "J2 (" << J[0][1] << ", "
    //                      << J[1][1] << ", " << J[2][1] << ")" << std::endl;
    //            std::cout << "jacDet " << jacDet << std::endl;
    //            */
    //          }
    //          else
    //          {
    //            // Find new local search direction for 3D -> 3D mapping using
    //            // Cramer's rule.
    //
    //            jacDet = + J[0][0]*J[1][1]*J[2][2] - J[0][0]*J[2][1]*J[1][2]
    //                    - J[1][0]*J[0][1]*J[2][2] + J[2][1]*J[1][0]*J[0][2]
    //                    - J[1][1]*J[2][0]*J[0][2] + J[2][0]*J[0][1]*J[1][2];
    //
    //            delta_xi[0] += + J[0][1]*J[2][2]*f[1] - J[0][2]*J[2][1]*f[1]
    //                           - J[1][1]*J[2][2]*f[0] + J[2][1]*J[1][2]*f[0]
    //                           - J[0][1]*J[1][2]*f[2] + J[0][2]*J[1][1]*f[2];
    //
    //            delta_xi[1] += - J[0][0]*J[2][2]*f[1] + J[2][0]*J[0][2]*f[1]
    //                           - J[1][2]*J[2][0]*f[0] + J[1][0]*J[2][2]*f[0]
    //                           + J[0][0]*J[1][2]*f[2] - J[1][0]*J[0][2]*f[2];
    //
    //            delta_xi[2] = - J[2][1]*J[1][0]*f[0] + J[0][0]*J[2][1]*f[1]
    //                          - J[0][0]*J[1][1]*f[2] + J[1][1]*J[2][0]*f[0]
    //                          + J[1][0]*J[0][1]*f[2] - J[2][0]*J[0][1]*f[1];
    //
    //            distNormalizer = std::pow(fabs(jacDet), 0.3333333);
    //          }
    //        }
    //
    //        // Here is the new local search direction.
    //        delta_xi[0] /= jacDet;
    //        delta_xi[1] /= jacDet;
    //        delta_xi[2] /= jacDet;
    //
    //        // If global element is smaller use local distance as a measure.
    //        // If global element is bigger use global distance as a measure.
    //        distMeasure = distNormalizer < 1 ? distance/distNormalizer : distance;
    //        if(distMeasure < TOL)
    //        {
    //          converged = true;
    //          break;
    //        }
    //
    //        /*
    //        std::cout << "delta_xi (" << delta_xi[0] << ", "
    //                  << delta_xi[1] << ")" << std::endl;
    //        */
    //
    //        // Perform damping iterations to find good damping coefficient
    //        xi_start = xi_k;
    //        coeff = 1;
    //        for(UInt l = 0; l < 8; l++)
    //        {
    //          xi_k = xi_start + delta_xi * coeff;
    //          Local2GlobalCoord(f, xi_k, coordMat, NULL);
    //
    //          f = f - globalPoint;
    //          distance_l = f.NormL2();
    //
    //          //          std::cout << "coeff " << coeff << " distance_l "
    //          //                    << distance_l << std::endl;
    //
    //          if(distance_l < distance)
    //          {
    //            distance = distance_l;
    //
    //            /*
    //            std::cout << "coeff " << coeff << " distance "
    //                      << distance << std::endl;
    //            */
    //
    //            // If global element is smaller use local distance as a measure.
    //            // If global element is bigger use global distance as a measure.
    //            distMeasure = distNormalizer < 1 ? distance/distNormalizer : distance;
    //            if(distMeasure < TOL)
    //            {
    //              converged = true;
    //              break;
    //            }
    //
    //            coeff /= 2;
    //          }
    //          else
    //          {
    //            divergence = true;
    //            break;
    //          }
    //        }
    //
    //        if(divergence)
    //          break;
    //
    //        if(converged)
    //          break;
    //
    //        iter++;
    //      }
    //      while(iter < 10);
    //
    //      // Put local coordinate of point into matrix.
    //      for(UInt l = 0; l < locDim; l++)
    //      {
    //        localCoords[l][i] = xi_k[l];
    //      }
    //
    //      /*
    //      std::cout << std::endl << "++++++ xi_k final (" << xi_k[0]
    //                << ", " << xi_k[1] << ", " << xi_k[2] << ") ++++++"
    //                << std::endl;
    //      std::cout << "f (" << f[0] << ", " << f[1] << ", " << f[2] << ")" << std::endl;
    //      //      std::cout << "distance " << distance << std::endl;
    //      std::cout << "divergence " << divergence << std::endl;
    //      std::cout << "distance " << distance << std::endl;
    //      std::cout << "dist_measure " << distMeasure << std::endl;
    //      */
    //    }
}

   
void LagrangeElemShapeMap::GetGlodMidPoint( Vector<Double>& locPoint ) {
  EXCEPTION("Not implemented");
}

Double LagrangeElemShapeMap::CalcVolume( ) {
  EXCEPTION("Not implemented");
  //    Double elemVol = 0;
  //    Double  jacDet, partVol;
  //    for (UInt actIntPt=1; actIntPt <= NumIntPoints_; actIntPt++) {
  //
  //      jacDet = CalcJacobianDetAtIp(actIntPt, CornerCoords, NULL);
  //
  //      if (isaxi) {
  //        Vector<Double> shapeFncAtIp;
  //        Vector<Double> CoordAtIP;
  //        GetShFncAtIp(shapeFncAtIp, actIntPt, NULL);
  //        CoordAtIP = CornerCoords * shapeFncAtIp;
  //        partVol = 2 * PI * IntWeights_[actIntPt-1] * jacDet * CoordAtIP[0];
  //      }
  //      else
  //        partVol = IntWeights_[actIntPt-1] * jacDet;
  //
  //      elemVol += partVol;
  //    }
  //
  //    return elemVol;
}

void LagrangeElemShapeMap::CalcNormal( Vector<Double>& normal, 
                                       const LocPoint& lp ) {
  EXCEPTION("Not implemented");
}


bool LagrangeElemShapeMap::
CoordIsInsideElem( const StdVector<LocPoint>& lps,
                     const Double tolerance,
                     StdVector<bool>& coordsInside ) {
  EXCEPTION("Not implemented");
}

void LagrangeElemShapeMap::CalcDiameter( Vector<Double>& diameter ) {
  EXCEPTION("Not implemented");
  // This method has to be implemented in the related finite element
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

void LagrangeElemShapeMap::GetMaxMinEdgeLength( ) {
  EXCEPTION("Not implemented");
}

void LagrangeElemShapeMap::GetEdgeLength( StdVector<Double>& edges_out) {
  EXCEPTION("Not implemented");
}

void LagrangeElemShapeMap::CalcJ( Matrix<Double>& jac, 
                                       const LocPoint& lp ) {
  Matrix<Double> deriv;
  ptFe_->GetDerivShFnc( deriv, lp, ptElem_ );
  jac = coords_ * deriv;
  jac *= depth_; // explicitly include depth_ of setup
}

Double LagrangeElemShapeMap::CalcJDet( Matrix<Double>& jac, 
                                       const LocPoint& lp ) {
  Matrix<Double> deriv;
  
  ptFe_->GetDerivShFnc( deriv, lp, ptElem_ );
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
   
   
void LagrangeElemShapeMap::SetElem( const Elem* ptElem ) {
  
  ptElem_ = ptElem;
  
  // get coordinates from grid
  ptGrid_->GetElemNodesCoord( coords_,ptElem->connect, isUpdated_ );

  // set reference element
  ptFe_ = feMap_[ptElem->type];
}

} // namespace CoupledField
