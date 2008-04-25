// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <iostream>
#include <fstream>
#include <string>

#include <string.h>
#include "basefe.hh"
#include "1D/line1fe.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#ifndef INTEGLIB
#include "olas.hh"
#endif
#include <Domain/elem.hh> 
#include <Utils/tools.hh> 
#include <Matrix/matrix.hh>
#include <DataInOut/WriteInfo.hh> 
 
namespace CoupledField
{

  BaseFE::BaseFE()
  {
  
    // initializing dynamic objects
    ShFncAtIp_             = NULL;
    ShFncDerivAtIp_        = NULL; 
    ShFncICModesDerivAtIp_ = NULL; 
    IntPoints_             = NULL; 

    // Create dummy ansatzFct object for lagrange functions
    shared_ptr<AnsatzFct> fct( new LagrangeFct() );
    actFct_ = fct;
    actNumFcns_ = 0;

    ICModes_     = false;
    CalcICModes_ = false;
    edgeIndices_ = NULL; 
    faceIndices_ = NULL; 
  }
 
  BaseFE::~BaseFE()
  {
 
    if( ShFncAtIp_ ) delete[] ShFncAtIp_;
    if( ShFncDerivAtIp_ ) delete[] ShFncDerivAtIp_;
    if( IntPoints_ ) delete[] IntPoints_;

    if(edgeIndices_ != NULL) { delete[] edgeIndices_; edgeIndices_ = NULL; } 
    if(faceIndices_ != NULL) { delete[] faceIndices_; faceIndices_ = NULL; } 
    if(ShFncICModesDerivAtIp_) { delete[]  ShFncICModesDerivAtIp_ ; ShFncICModesDerivAtIp_ = NULL; }    
    
    // delete our map. The content of normal integration rules are static array
    // the cartesian rules are complete dynamic
    std::map<const std::string, StdVector<Double*>*>::iterator iter; 
    for(iter = IntegrationPointsMap_.begin(); 
        iter != IntegrationPointsMap_.end(); iter++)  {
      StdVector<Double*>* data = iter->second; 
      if(data != NULL) {
        if(iter->first.find("Cartesian", 0 ) != std::string::npos)  { 
          for(UInt i = 0; i < data->GetSize(); i++)  {
            delete (*data)[i];
          } 
        }
        
        delete data;
      }    
    }
    
  }
  
  void BaseFE::GetShFnc(Vector<double> & S, 
                          const Vector<Double> & LCoord,
                          const Elem* elem,
                          UInt dof  )
  {

    CalcShapeFnc(S, LCoord, elem, dof );

  }
  
  void BaseFE::Local2GlobalCoord(Vector<Double> & globCoord,
                                   const Vector<Double> & locCoord, 
                                   const Matrix<Double> & coordMat,
                                   const Elem* elem ) {
    
    
    // step 1: evaluate shape fncs. at local coordinate
    Vector<Double> shFnc;
    CalcShapeFnc(shFnc, locCoord, elem, 1, AnsatzFct::NODE);

    // step2: multiply shape fncs for each dimension with according matrix entries
    globCoord.Resize(coordMat.GetSizeRow());
    coordMat.Mult(shFnc,globCoord);    
  }

  void BaseFE::Global2LocalCoords(Matrix<Double> & localCoords,
                                    const Matrix<Double> & globalCoords,
                                    const Matrix<Double> & coordMat ) {
    
    
    Vector<Double> globalPoint; // global point coordinates
    UInt globDim = globalCoords.GetSizeRow(); // determine global dimension
    UInt numPoints = globalCoords.GetSizeCol(); // number of global points
    UInt locDim = LCornerCoords_.GetSizeRow(); // dimension of current element
    //    UInt numCorners = LCornerCoords_.GetSizeCol(); // number of element corners

    Vector<Double> xi_start; // local start point for Newton-Raphson method
    Vector<Double> xi_k; // local point at iteration k
    Vector<Double> delta_xi; // local search direction
    Vector<Double> f; // global (negative) search direction
    Vector<Double> f_start; // start value for global search direction 
    Matrix<Double> J; // Jacobian at local point xi_k
    Double f_min; // minimal distance between global corners and global point
    Double distance; // global distance!
    Double distance_l; // global distance at iteration l
    Double distMeasure; // global distance or local distance
                        // depending on sqrt(|jacDet|)
    Double TOL = 0.0000001; // tolerance for global distance
    Double jacDet; // denominator for Cramer's rule.
    Double distNormalizer; // denominator for Cramer's rule.
    Double coeff; // damping coefficient
    bool divergence; // does the Newton-Raphson algorithm diverge?
    bool converged; // have we found the local point?
    UInt iter = 0;
   
    // Initialize variables
    globalPoint.Resize(globDim);
    xi_start.Resize(3); xi_start.Init();
    xi_k.Resize(3); xi_k.Init();
    delta_xi.Resize(3); ; delta_xi.Init();
    localCoords.Resize(locDim, numPoints);
    J.Resize(3, 3); J.Init();

    // Perform Newton-Raphson method on the list of global points
    // to find local coordinates within this element.
    for(UInt i = 0; i < numPoints; i++)
    {
      // Copy global point coordinates into a Vector for algebraic ops.
      for(UInt j = 0; j < globDim; j++)
      {
        globalPoint[j] = globalCoords[j][i];
      }

      // Find good startpoint xi_k among local corner coordinates
      f_min = 999e5; // really big value!
      for(UInt k = 0; k < NumCorners_; k++)
      {
        for(UInt l = 0; l < locDim; l++)
        {
          xi_k[l] = LCornerCoords_[l][k];
        }
        
        Local2GlobalCoord(f, xi_k, coordMat, NULL);
        f = f - globalPoint;
        distance = f.NormL2();
        if( distance < f_min )
        {
          f_min = distance;
          xi_start = xi_k;
          f_start = f;
        }
      }

      /*
      std::cout << "Beginning to find local coords..." << std::endl;
      std::cout << "Global coords..." << globalPoint << std::endl;
      */
      xi_k = xi_start;
      f = f_start;
      distance = f_min;
      divergence = false;
      converged = false;
      
      do
      {
        // Calculate Jacobian at iteration point xi_k
        CalcJacobian(J, xi_k, coordMat, NULL );

        /*
        std::cout << std::endl << "++++++ xi_k (" << xi_k[0] << ", "
                  << xi_k[1] << ", " << xi_k[2] << ") ++++++" << std::endl;
        std::cout << "f (" << f[0] << ", " << f[1]
                  << ", " << f[2] << ")" << std::endl;
        std::cout << "distance " << distance << std::endl;
        */

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
          jacDet = + J[0][0]*J[1][1] - J[1][0]*J[0][1];

          delta_xi[0] = - J[1][1]*f[0] + J[0][1]*f[1];
        
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
            //            Vector<Double> normal;

            jacDet = + J[1][0] * J[2][1] - J[2][0] * J[1][1]
                    - J[0][0] * J[2][1] + J[2][0] * J[0][1]
                    + J[0][0] * J[1][1] - J[1][0] * J[0][1];

            delta_xi[0] = (f[0] * J[0][0] +
                           f[1] * J[1][0] +
                           f[2] * J[2][0]);            

            delta_xi[1] = (f[0] * J[0][1] +
                           f[1] * J[1][1] +
                           f[2] * J[2][1]);            

            distNormalizer = sqrt(fabs(jacDet));

            /*
            std::cout << "J1 (" << J[0][0] << ", "
                      << J[1][0] << ", " << J[2][0] << ")" << std::endl;
            std::cout << "J2 (" << J[0][1] << ", "
                      << J[1][1] << ", " << J[2][1] << ")" << std::endl;
            std::cout << "jacDet " << jacDet << std::endl;
            */
          }
          else
          {
            // Find new local search direction for 3D -> 3D mapping using
            // Cramer's rule.
            
            jacDet = + J[0][0]*J[1][1]*J[2][2] - J[0][0]*J[2][1]*J[1][2]
                    - J[1][0]*J[0][1]*J[2][2] + J[2][1]*J[1][0]*J[0][2]
                    - J[1][1]*J[2][0]*J[0][2] + J[2][0]*J[0][1]*J[1][2];

            delta_xi[0] += + J[0][1]*J[2][2]*f[1] - J[0][2]*J[2][1]*f[1]
                           - J[1][1]*J[2][2]*f[0] + J[2][1]*J[1][2]*f[0]
                           - J[0][1]*J[1][2]*f[2] + J[0][2]*J[1][1]*f[2];
          
            delta_xi[1] += - J[0][0]*J[2][2]*f[1] + J[2][0]*J[0][2]*f[1]
                           - J[1][2]*J[2][0]*f[0] + J[1][0]*J[2][2]*f[0]
                           + J[0][0]*J[1][2]*f[2] - J[1][0]*J[0][2]*f[2];
          
            delta_xi[2] = - J[2][1]*J[1][0]*f[0] + J[0][0]*J[2][1]*f[1]
                          - J[0][0]*J[1][1]*f[2] + J[1][1]*J[2][0]*f[0]
                          + J[1][0]*J[0][1]*f[2] - J[2][0]*J[0][1]*f[1];

            distNormalizer = std::pow(fabs(jacDet), 0.3333333);
          }
        }

        // Here is the new local search direction.
        delta_xi[0] /= jacDet;
        delta_xi[1] /= jacDet;
        delta_xi[2] /= jacDet;

        // If global element is smaller use local distance as a measure.
        // If global element is bigger use global distance as a measure.
        distMeasure = distNormalizer < 1 ? distance/distNormalizer : distance;
        if(distMeasure < TOL)
        {
          converged = true;
          break;
        }

        /*        
        std::cout << "delta_xi (" << delta_xi[0] << ", "
                  << delta_xi[1] << ")" << std::endl;
        */

        // Perform damping iterations to find good damping coefficient
        xi_start = xi_k;
        coeff = 1;
        for(UInt l = 0; l < 8; l++)
        {
          xi_k = xi_start + delta_xi * coeff;
          Local2GlobalCoord(f, xi_k, coordMat, NULL);

          f = f - globalPoint;
          distance_l = f.NormL2();

          //          std::cout << "coeff " << coeff << " distance_l "
          //                    << distance_l << std::endl;

          if(distance_l < distance)
          {
            distance = distance_l;

            /*
            std::cout << "coeff " << coeff << " distance "
                      << distance << std::endl;
            */

            // If global element is smaller use local distance as a measure.
            // If global element is bigger use global distance as a measure.
            distMeasure = distNormalizer < 1 ? distance/distNormalizer : distance;
            if(distMeasure < TOL)
            {
              converged = true;
              break;
            }

            coeff /= 2;
          }
          else
          {
            divergence = true;
            break;
          }
        }

        if(divergence)
          break;

        if(converged)
          break;

        iter++;
      }
      while(iter < 10);

      // Put local coordinate of point into matrix.
      for(UInt l = 0; l < locDim; l++)
      {
        localCoords[l][i] = xi_k[l];
      }

      /*
      std::cout << std::endl << "++++++ xi_k final (" << xi_k[0]
                << ", " << xi_k[1] << ", " << xi_k[2] << ") ++++++"
                << std::endl;
      std::cout << "f (" << f[0] << ", " << f[1] << ", " << f[2] << ")" << std::endl;
      //      std::cout << "distance " << distance << std::endl;
      std::cout << "divergence " << divergence << std::endl;
      std::cout << "distance " << distance << std::endl;
      std::cout << "dist_measure " << distMeasure << std::endl;
      */
    }
  }


  void BaseFE::GetShFncAtIp(Vector<Double> & S, const UInt ip,
                              const Elem * elem, UInt dof )
  {
    
    if( actFct_->GetType() == AnsatzFct::LAGRANGE ) {
      S = ShFncAtIp_[ip-1];
    } else {
      CalcShapeFnc( S, IntPoints_[ip-1], elem, dof, 
                    AnsatzFct::ALL );
    }
  }

  void BaseFE::GetGlobDerivShFnc(Matrix<Double>  &Deriv, 
                                   const Vector<Double> & LCoord,
                                   const Matrix<Double> & CornerCoords,
                                   const Elem * elem, 
                                   UInt dof )
  {

    
    CalcInvJacobian(JInv, LCoord, CornerCoords, elem );
    CalcLocalDerivShapeFnc(LDeriv, LCoord, elem, dof, AnsatzFct::ALL);
    Deriv = LDeriv * JInv;
  }

  void BaseFE::GetGlobDerivShFncAtIp(Matrix<Double> & Deriv, 
                                       const UInt ip,
                                       const Matrix<Double> & CornerCoords,
                                       Double & jacDet,
                                       const Elem * elem, 
                                       UInt dof )
  {

    std::string errMsg;

    //  Deriv.Resize(NumNodes_,Dim_);
    Double JInvDet;
  
    CalcInvJacobianAtIp(JInv, ip, CornerCoords, elem);

    if( actFct_->GetType() == AnsatzFct::LAGRANGE ) {
      Deriv = ShFncDerivAtIp_[ip-1] * JInv;
    } else {
      CalcLocalDerivShapeFnc( LDeriv, IntPoints_[ip-1],
                              elem, dof, AnsatzFct::ALL );
      Deriv = LDeriv * JInv;
    }
        
    // det(A) = 1 / det(A^(-1))
    JInv.Determinant(JInvDet);
    jacDet = 1.0 / JInvDet;

    if ( jacDet < 0.0 ){
      EXCEPTION("BaseFE:GetGlobDerivShFncAtIp: Negative Jacobian Determinant "
                << "at element " << elem->elemNum 
                << " with connectivity " << elem->connect.Serialize());
    }
    
  }


  void BaseFE::GetGlobDerivShFncAtIp(Matrix<Double> & Deriv, 
                                       const UInt ip,
                                       const Matrix<Double> & CornerCoords,
                                       const Elem* elem,
                                       UInt dof )
  {
  
    //  Deriv.Resize(NumNodes_,Dim_);

    if ( CalcICModes_ && ICModes_ ) {
      //incompatible modes
      Vector<Double> LCoord( CornerCoords.GetSizeRow());
      LCoord.Init();
      CalcInvJacobian(JInv, LCoord, CornerCoords, elem);
      Deriv = ShFncICModesDerivAtIp_[ip-1] * JInv;
    }
    else {
      CalcInvJacobianAtIp(JInv, ip, CornerCoords, elem);

      if( actFct_->GetType() == AnsatzFct::LAGRANGE ) {
        Deriv = ShFncDerivAtIp_[ip-1] * JInv;
      } 
      else {
        CalcLocalDerivShapeFnc( LDeriv, IntPoints_[ip-1],
                                elem, dof, AnsatzFct::ALL );
        Deriv = LDeriv * JInv;
      }
      //std::cerr << "Deriv = \n" << Deriv << std::endl;
    }
  }




  void BaseFE::CalcJacobian(Matrix<Double> & J, 
                              const Vector<Double> & LCoord, 
                              const Matrix<Double> & CornerCoords,
                              const Elem* elem )
  {

    //  J.Resize(Dim_,Dim_);

    CalcLocalDerivShapeFnc(LDeriv, LCoord, elem, 1,AnsatzFct::NODE );
    J = CornerCoords * LDeriv;
  }


  void BaseFE::CalcJacobianAtIp(Matrix<Double> & J, 
                                  const UInt ip, 
                                  const Matrix<Double> & CornerCoords,
                                  const Elem* elem)
  {

    if (CornerCoords.GetSizeRow()==3 && Dim_==2) // Surface element in 3D
      J.Resize(CornerCoords.GetSizeRow(),Dim_);

    J = CornerCoords * ShFncDerivAtIp_[ip-1];
  }

  Double BaseFE::CalcJacobianDet(const Vector<Double> & LCoord,
                                   const Matrix<Double> & CornerCoords,
                                   const Elem* elem )
  {

    std::string errMsg;
    Double jacDet;

    CalcJacobian( J, LCoord, CornerCoords, elem  );
    //    J.Determinant(jacDet);

    if (CornerCoords.GetSizeRow()==3 && Dim_==2)
    {
      Vector<Double> normal;
      normal.Resize(CornerCoords.GetSizeRow());  
      normal[0]= J[1][0]* J[2][1]- J[2][0]* J[1][1];
      normal[1]=J[2][0]*J[0][1]- J[0][0]* J[2][1];  
      normal[2]= J[0][0]* J[1][1]- J[1][0]*J[0][1];

      jacDet = sqrt(normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]);
    }
    else  {
      J.Determinant(jacDet);
    }

    if ( jacDet < 0.0 ){
      EXCEPTION("BaseFE:CalcJacobianDet: Negative Jacobian Determinant "
                << "at element " << elem->elemNum 
                << " with connectivity " << elem->connect.Serialize());
    }

    return jacDet;
  }


  Double BaseFE::CalcJacobianDetAtIp(const UInt ip, 
                                       const Matrix<Double> & CornerCoords,
                                       const Elem* elem)
  {

    std::string errMsg;

    CalcJacobianAtIp( J, ip, CornerCoords, elem);

    if (CornerCoords.GetSizeRow()==3 && Dim_==2)
    {
      Vector<Double> normal;
      normal.Resize(CornerCoords.GetSizeRow());  
      normal[0]= J[1][0]* J[2][1]- J[2][0]* J[1][1];
      normal[1]=J[2][0]*J[0][1]- J[0][0]* J[2][1];  
      normal[2]= J[0][0]* J[1][1]- J[1][0]*J[0][1];

      Double detJ = sqrt(normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]);

      if ( detJ < 0.0 ){
        EXCEPTION("BaseFE::CalcJacobianDetAtIp: Negative Jacobian Determinant "
                  << "at element " << elem->elemNum 
                  << " with connectivity " << elem->connect.Serialize());
      }
      return detJ;
    }

    else  {
      Double jacDet ;
      J.Determinant(jacDet);
      if ( jacDet < 0.0 ){
        EXCEPTION("BaseFE::CalcJacobianDetAtIp: Negative Jacobian Determinant "
                  << "at element " << elem->elemNum 
                  << " with connectivity " << elem->connect.Serialize());
      }
      return jacDet;
    }  
  }




  void BaseFE::CalcInvJacobian(Matrix<Double> & JInv,
                                 const Vector<Double> & LCoord,
                                 const Matrix<Double> & CornerCoords,
                                 const Elem* elem )
  {
  
    Matrix<Double> J, LDeriv;
    Double det;

    CalcLocalDerivShapeFnc(LDeriv, LCoord, elem, 1,AnsatzFct::NODE);

    J = CornerCoords * LDeriv;

    if ( CalcICModes_ && ICModes_ ) {
      // modified inverse: do not divide by the determinant!!
      J.Determinant( det );
      J.Invert(JInv);
      JInv *= det;
    }
    else 
      J.Invert(JInv);
  }



 
  void BaseFE::CalcInvJacobianAtIp(Matrix<Double> & JInv,
                                     const UInt ip,
                                     const Matrix<Double> & CornerCoords,
                                     const Elem* elem)
  {

    JInv.Resize(Dim_,Dim_);

    //  J.Resize(Dim_,Dim_);

    J = CornerCoords * ShFncDerivAtIp_[ip-1];

    J.Invert(JInv);
  }


  void BaseFE::SetShapeFncAtIp()
  {
  
    if (!ShFncAtIp_) {
      ShFncAtIp_ = new Vector<Double>[NumIntPoints_];
    } else{ 
      delete[] ShFncAtIp_ ;
      ShFncAtIp_ = new Vector<Double>[NumIntPoints_];
    }


    for( UInt i=0; i<NumIntPoints_; i++ )
    {
      CalcShapeFnc( ShFncAtIp_[i], IntPoints_[i], 
                    NULL, 1, AnsatzFct::NODE );
    }
  }
  
  void BaseFE::SetShapeFncDerivAtIp()
  {

    if( !ShFncDerivAtIp_) {
      ShFncDerivAtIp_ = new Matrix<Double>[NumIntPoints_]; 
    } 
    else{ 
      delete[] ShFncDerivAtIp_ ;
      ShFncDerivAtIp_ = new Matrix<Double>[NumIntPoints_];
    }

    for( UInt i=0; i<NumIntPoints_; i++ )
      CalcLocalDerivShapeFnc( ShFncDerivAtIp_[i], IntPoints_[i], 
                              NULL, 1, AnsatzFct::NODE );

    // check, if incompatible modes are used
    if ( ICModes_ ) {
      if( !ShFncICModesDerivAtIp_ ) {
        ShFncICModesDerivAtIp_ = new Matrix<Double>[NumIntPoints_]; 
      } 
      else{ 
        delete[]  ShFncICModesDerivAtIp_ ;
        ShFncICModesDerivAtIp_ = new Matrix<Double>[NumIntPoints_];
      }
      for( UInt i=0; i<NumIntPoints_; i++ ) {
	CalcLocalICModesDerivShapeFnc(  ShFncICModesDerivAtIp_[i], IntPoints_[i], 
					NULL, 1, AnsatzFct::NODE );
      }
    }
  }



  void BaseFE::GetGlobalEdgeIndicesAtIP( Vector<Double> & globCoord, UInt ip,
                                         const Matrix<Double> & cornerCoords)
  {
    // cornerCoords: nrCorners x dim

    globCoord.Resize(Dim_);
  
    globCoord =  cornerCoords * ShFncAtIp_[ip-1];
  }




  // calculate global derivates of edge shape functions in integration point ip
  void BaseFE::GetEdgeGlobDerivShFncAtIp(StdVector< Matrix<Double>* > & deriv, 
                                         const UInt ip,
                                         const Matrix<Double> & cornerCoords)
  {

    // vector of coordinates of the desired integration point
    Vector<Double> lCoord;
    lCoord.Resize(Dim_);

    for (UInt i=0; i<Dim_; i++)
      lCoord[i] = IntPoints_[ip-1][i];
  
    GetEdgeGlobalDerivShapeFnc(deriv, lCoord, cornerCoords);
  }

  void BaseFE::CalcEdgeShapeFncAtIp(Matrix<Double> & shape, 
                                    const UInt ip,
                                    const Matrix<Double> & cornerCoords)
  {

    Vector<Double> lCoord;
    lCoord.Resize(Dim_);

    for (UInt i=0; i<Dim_; i++)
      lCoord[i] = IntPoints_[ip-1][i];
  
    CalcEdgeShapeFnc(shape, lCoord, cornerCoords);
  }


  void BaseFE::GetGlobalEdgeIndices(StdVector<UInt> & globEdgeIndex,
                                    UInt * pDENodes, 
                                    BaseSystem * algsys)
  {

    EXCEPTION( "Edge functions are currently not supported!");

    // define the global edge number
    //for(UInt actEdge=0; actEdge < globEdgeIndex.GetSize(); actEdge++)
    //  globEdgeIndex[actEdge] = 
    //    algsys->GetNode2Edge(pDENodes[ edgeVertices_[actEdge][0]],
    // pDENodes[ edgeVertices_[actEdge][1]]);
  }

  std::string BaseFE::CoordMatrix2String(const Matrix<Double> & coordMat)
  {
    std::string ret;
    for (UInt j=0; j<coordMat.GetSizeCol(); j++) {
      ret += "(";
      for (UInt i=0; i<coordMat.GetSizeRow()-1; i++) {
        ret += GenStr(coordMat[i][j]);
        ret += ", ";
      }
      ret += GenStr(coordMat[coordMat.GetSizeRow()-1][j]);
      ret +=")\n";
    }
    return ret;
  }


  Double BaseFE::CalcVolume(const Matrix<Double> & CornerCoords, 
                              const bool isaxi)
  {

    Double elemVol = 0;
    Double  jacDet, partVol;
    for (UInt actIntPt=1; actIntPt <= NumIntPoints_; actIntPt++) {

      jacDet = CalcJacobianDetAtIp(actIntPt, CornerCoords, NULL);
        
      if (isaxi) {
        Vector<Double> shapeFncAtIp;
        Vector<Double> CoordAtIP;
        GetShFncAtIp(shapeFncAtIp, actIntPt, NULL);
        CoordAtIP = CornerCoords * shapeFncAtIp;
        partVol = 2 * PI * IntWeights_[actIntPt-1] * jacDet * CoordAtIP[0];
      }
      else 
        partVol = IntWeights_[actIntPt-1] * jacDet;
    
      elemVol += partVol;
    }
  
    return elemVol;
  }

  void BaseFE::CalcBarycenter(const Matrix<Double>& coords, Point& barycenter) 
  {
    UInt n_dims  = coords.GetSizeRow();
    UInt n_elems = coords.GetSizeCol();
       
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


  /////////////////////////////////////////////////////////////////////
  /////// For geometrical transformation using direction cosines //////

  void BaseFE::CoordTrans( const Matrix<Double> &ptCoord, 
                           Matrix<Double> &TransMat, 
                           Matrix<Double> &ShellCoord )
  {

    // Method based on book XXX and Master thesis from Xiong //
    const Integer spaceDim = 3;
    std::vector<Double> Vx, Vy, Vz;
    const UInt row = ptCoord.GetSizeRow();
    const UInt col = ptCoord.GetSizeCol();

    TransMat.Resize(spaceDim,true);
    TransMat.Init();

    Matrix<Double> NewCoord, temp;
    NewCoord.Resize( row, col );
    NewCoord.Init();
    temp.Resize( row, col );
    temp.Init();

    ShellCoord.Resize( row - 1, col );
    ShellCoord.Init();

    Double length;

    Vx.push_back( ptCoord[0][1] - ptCoord[0][0] );
    Vx.push_back( ptCoord[1][1] - ptCoord[1][0] );
    Vx.push_back( ptCoord[2][1] - ptCoord[2][0] );

    length = 1.0 / sqrt( Vx[0] * Vx[0] + Vx[1] * Vx[1] + Vx[2] * Vx[2] );
 
    //   cos(x_,x), cos(x_,y), cos(x_,z)
    TransMat[0][0] = Vx[0] * length;
    TransMat[0][1] = Vx[1] * length;
    TransMat[0][2] = Vx[2] * length;

    Vy.push_back( ptCoord[0][2] - ptCoord[0][0] );
    Vy.push_back( ptCoord[1][2] - ptCoord[1][0] );
    Vy.push_back( ptCoord[2][2] - ptCoord[2][0] );

    Vz.push_back( Vx[1] * Vy[2] - Vx[2] * Vy[1] );
    Vz.push_back( Vx[2] * Vy[0] - Vx[0] * Vy[2] );
    Vz.push_back( Vx[0] * Vy[1] - Vx[1] * Vy[0] );

    length = 1.0 / sqrt( Vz[0] * Vz[0] + Vz[1] * Vz[1] + Vz[2] * Vz[2] );

    TransMat[2][0] = Vz[0] * length;
    TransMat[2][1] = Vz[1] * length;
    TransMat[2][2] = Vz[2] * length;

    TransMat[1][0] = TransMat[0][2] * TransMat[2][1] - TransMat[0][1] 
                     * TransMat[2][2];
    TransMat[1][1] = TransMat[0][0] * TransMat[2][2] - TransMat[0][2] 
                     * TransMat[2][0];
    TransMat[1][2] = TransMat[0][1] * TransMat[2][0] - TransMat[0][0] 
                     * TransMat[2][1];

    // transform geometry from real(global) coordinate to standard(local) 
    // coordinate

    for( int i = row - 1; i >= 0; i-- )
      for( int j = col - 1; j >= 0; j-- )
	//transform
	temp[i][j] = ptCoord[i][j] - ptCoord[i][0];

    //     std::cout << "The new coordinate is\n" << temp << std::endl;

    //     std::cout << "The base TransMatrix is\n" << TransMat << std::endl;

    NewCoord = TransMat * temp;

    //     std::cout << "The 3D LocalCoord matrix is\n" << NewCoord << std::endl;

    //rotate
    for( UInt i = 0; i < row - 1; i++)
      for( UInt j = 0; j < col; j++)
	ShellCoord[i][j] = NewCoord[i][j];

  }

  void BaseFE::CoordTrans2D( const Matrix<Double> &ptCoord, 
			     Matrix<Double> &TransMat, 
			     Matrix<Double> &ShellCoord )
  {

    // Method based on book XXX and Master thesis from Xiong //
    const Integer spaceDim = 2;
    std::vector<Double> Vx, Vy;
    const UInt row = ptCoord.GetSizeRow();
    const UInt col = ptCoord.GetSizeCol();

    TransMat.Resize(spaceDim);
    TransMat.Init();

    Matrix<Double> NewCoord, temp;
    NewCoord.Resize( row, col );
    NewCoord.Init();
    temp.Resize( row, col );
    temp.Init();

    ShellCoord.Resize( row - 1, col );
    ShellCoord.Init();

    Double length;

    Vx.push_back( ptCoord[0][1] - ptCoord[0][0] );
    Vx.push_back( ptCoord[1][1] - ptCoord[1][0] );
 

    length = 1.0 / sqrt( Vx[0] * Vx[0] + Vx[1] * Vx[1] );
 
    //   cos(x_,x), cos(x_,y), cos(y_,x), cos(y_,y)
    TransMat[0][0] = Vx[0] * length;
    TransMat[0][1] = Vx[1] * length;
    TransMat[1][0] = Vx[1] * length;
    TransMat[1][1] = Vx[0] * length;


    // transform geometry from real(global) coordinate to standard(local) 
    // coordinate

    for( UInt i = row - 1; i >= 0; i-- )
      for( UInt j = col - 1; j >= 0; j-- )
	//transform
	temp[i][j] = ptCoord[i][j] - ptCoord[i][0];

    std::cout << "The new coordinate is\n" << temp << std::endl;

    std::cout << "The base TransMatrix is\n" << TransMat << std::endl;

    NewCoord = TransMat * temp;

    std::cout << "The 2D LocalCoord matrix is\n" << NewCoord << std::endl;

    //rotate
    for( UInt i = 0; i < row - 1; i++)
      for( UInt j = 0; j < col; j++)
	ShellCoord[i][j] = NewCoord[i][j];

  }
  
   
  void BaseFE::MakeKey(IntegrationMethod type, int order, std::string &out)
  {
    Enum2String(type, out);       
    std::stringstream ss;
    ss << GetShapeName() << " (" << out << ") order " << order;
    out.assign(ss.str());
  }   

  /** private order encoder. Makes a check and exits on error */
  int BaseFE::EncodeCartesianOrder(int order1, int order2, int order3)
  {
    if(order1 > 99 || order2 > 99 || order3 > 99) 
      EXCEPTION("Cartesian product numerical integration only with orders < 99");
        
    return order1 + 100 * order2 + (order3 > 0 && Dim_ == 3 ? 10000 : 0) * order3;
  }

  /** private order decoder */
  void BaseFE::DecodeCartesianOrder(int encoded_order, int* order1, int* order2, int* order3)
  {
    if(encoded_order <= 11 || encoded_order > 999999) 
      EXCEPTION("Invalid encoded cartesian integration order");
           
           
    *order3 = (encoded_order >= 10000) ? encoded_order/10000 : 0;   
    encoded_order -= 10000 * (*order3);
        
    *order2 = encoded_order/100;
    encoded_order -= 100 * (*order2);        
        
    *order1 = encoded_order;
        
    // we cannot set it to 0 before as order2 needs to be < 19
    if(Dim_ != 3) *order3 = 0;
  }


  void BaseFE::MakeKey(int order1, int order2, int order3, std::string &out)
  {
    // this is only for debugging, we can encode the orders as defined in environment.hh
    int order = EncodeCartesianOrder(order1, order2, order3);
       
    MakeKey(CARTESIAN, order, out);
  }   



  /** check CreateCartesian... for similar implementation! */   
  void BaseFE::AddIntegrationPoints(IntegrationMethod type, int order, 
                                    int numberOfRows, Double* data)
  {

    std::string key;
    MakeKey(type, order, key);
        
    // check the key for uniqueness
    if(IntegrationPointsMap_.find(key) != IntegrationPointsMap_.end())
    {
      key.append(" is already in BaseFE::IntegrationPointsMap_");
      EXCEPTION(key.c_str()); 
    }

    // create the payload
    StdVector<Double*>*  value = new StdVector<Double*>(numberOfRows);
    for(int i = 0; i < numberOfRows; i++) (*value)[i]=(data + i * (Dim_+1)); // the data is Dim_ * coordintates + weight

    IntegrationPointsMap_[key] = value;  

    // sum up and print the weights for debugging         
    //Double sum = 0;
    //for(int i = 0; i < numberOfRows; i++) sum += *((*value)[i]+Dim_);
    //std::cout << key << " has weight sum " << sum << std::endl;
       
  }

  StdVector<Double*>* BaseFE::GetIntegrationPoints(IntegrationMethod type, int order, bool search_upwards, 
                                                   bool search_downwards, bool fallback)     
  {
    std::string key;
    int org_order = order;

    do
    {
      MakeKey(type, order, key);
          
      // store current value so the state is set to what we have
      IntegMethod = type;
      IntegOrder  = order;
   	      
      if(IntegrationPointsMap_.find(key) != IntegrationPointsMap_.end())
      {
        if(order != org_order)
        {
          *warning << "Use Integration " << key << " instead of the wanted order " << org_order << std::endl;
          Warning(__FILE__, __LINE__);
        }    
        return IntegrationPointsMap_[key];
      }   
   	      
      // nothing found
      if(search_upwards) order++;
      if(search_downwards) order--;   
      if(search_upwards && search_downwards)
        EXCEPTION("searching up- and downward concurrently!! Stupid!");
    }    	
    while(order > 0 && order < 42 && (search_upwards || search_downwards));
    
    // this is a call from the fallback-part
    if(!fallback) return NULL;  

    // print msg;
    MakeKey(type, org_order, key);
    *warning << "No integration points defined for " << key << "! Available are: ";
    std::map<const std::string, StdVector<Double*>*>::iterator iter;
   
    for(iter = IntegrationPointsMap_.begin(); iter != IntegrationPointsMap_.end(); iter++) {
      *warning << iter->first << "\n";
    }                    
    Warning(__FILE__, __LINE__);
           
    // use default
    // restrict to order 2 but stay on 1 if desired
    int alt_o = IntegOrder >= 2 ? 2 : 1;

    // when no economical already use it - if economical doesn't work use classic
    IntegrationMethod alt_m = IntegOrder == ECONOMICAL ? CLASSICAL : ECONOMICAL;
    StdVector<Double*>* data = GetIntegrationPoints(alt_m, alt_o, false, false, false);
      
    // more tries
    if(data == NULL)
    {
      // try again with flipped method
      alt_m = alt_m == ECONOMICAL ? CLASSICAL : ECONOMICAL;          
      data = GetIntegrationPoints(alt_m, alt_o, false, false, false);
    }

    // special pyramid handling :(
    if(data == NULL && alt_o == 1)
    {
      alt_o = 2;
      alt_m = CLASSICAL;
      data = GetIntegrationPoints(alt_m, alt_o, false, false, false);
    } 

    // no more tries   
    if(data != NULL)
    {
      MakeKey(alt_m, alt_o, key);
      *warning << "Use fallback " << key << " for integration\n";
      Warning(__FILE__, __LINE__);
      return data;
    }
          
    EXCEPTION("No default integration found");  
    exit(-1);
  } 

  void BaseFE::CommonInit(IntegrationMethod method, int order)
  {
    
 #ifndef INTEGLIB
    // if undefined we have an empty constructor and search :)
    if(method == UNDEFINED)
    {  
      // when we habe a XML setting use it
      // check defaults in case of legacy XML files

      // Check, if element "integRules" was defined
      if(!param) {
        SetDefaultIntegration();
      } else {
        ParamNode * integNode = param->Get( "integRules", false );
        if(integNode)    
        { 
          std::string type = integNode->Get("type")->AsString();
          // we have values in the XML file
          String2Enum(type,IntegMethod );    

          std::string str;
          integNode->Get( "order", str );
          IntegOrder = String2Int(str);    
        }    
        else
        {
          // The default Integration rules
          SetDefaultIntegration();
        }
      }
    } 
    else
    {
      IntegMethod = method;
      IntegOrder  = order;
    }
 #else
    // The default Integration rules
    SetDefaultIntegration();
 #endif        

    // first set integration points and corner coords ...
    // load all into the map
    FillIntegrationPoints();    
    // get the values by IntegMethod and IntegOrder
    SetIntPoints();
    SetCornerCoords();

    // ... then calc shape function values at integration points
    SetShapeFncAtIp();
    SetShapeFncDerivAtIp();      
  }

  void BaseFE::DumpIntegrationPoints()
  {
    std::string str;
    Enum2String(IntegMethod, str);
    
    std::cout << "Current integration points for " << GetShapeName() << " method=" << str << " order=" << IntegOrder << std::endl;

    for(UInt ip=0; ip<NumIntPoints_; ip++)
    {         
      for(UInt dim = 0; dim < Dim_; dim++)
      {
        std::cout << ip << ":" << dim << "=" << IntPoints_[ip][dim] << "  ";
      }
          
      std::cout << "weight=" << IntWeights_[ip] << std::endl;
    }
  }

  void BaseFE::SetCartesianInteg(int order1, int order2, int order3, bool create_only)
  {
    // check if we already have
    std::string key;
    MakeKey(order1, order2, order3, key);

    if(IntegrationPointsMap_.find(key) == IntegrationPointsMap_.end())
    {
      // doesn't exist -> create
      Line1FE* line1 = new Line1FE(ECONOMICAL, order1);
      Line1FE* line2 = new Line1FE(ECONOMICAL, order2);
      Line1FE* line3 = Dim_ == 3 ? new Line1FE(ECONOMICAL, order3) : NULL;
          
      // store in map 
      CreateCartesianIntegration(line1, line2, line3);
          
      // destroy temp objects
      if(line3) delete line3;
      delete line2;
      delete line1;
    }   
      
    // set the state variables so SetIntPoints can load the data from the map
    IntegMethod = CARTESIAN;
    IntegOrder  = EncodeCartesianOrder(order1, order2, order3);
      
    // false only for the call from SetIntPoints()
    if(!create_only) SetIntPoints(); // this calls this method again but the map entry will already the set
  }
  

  /** Creates the integration points via cartesian product for 2D and 3D elements. The result is stored in the map.
   *  @param element1 IntegMethod and IntegOrder shall be set
   *  @param element2 IntegMethod and IntegOrder shall be set
   *  @param element3 the only optional element */
  void BaseFE::CreateCartesianIntegration(BaseFE* line1, BaseFE* line2, BaseFE* line3)
  {
    // we store in map format, as such we have dim corrdinates + weight
    int rows = line1->NumIntPoints_ * line2->NumIntPoints_ * (Dim_ == 3 ? line3->NumIntPoints_ : 1);

    std::string key;
    MakeKey(line1->IntegOrder, line2->IntegOrder, line3 != NULL ? line3->IntegOrder : 0, key);
        
    // check the key for uniqueness
    if(IntegrationPointsMap_.find(key) != IntegrationPointsMap_.end())
    {
      key.append(" is already in BaseFE::IntegrationPointsMap_");
      EXCEPTION(key.c_str()); 
    }
     
    StdVector<Double*>* data = new StdVector<Double*>(rows);
    Double* row = NULL;
     
    int counter = 0;
    for(UInt i = 0; i < line1->NumIntPoints_; i++)
    {
      for(UInt j = 0; j < line2->NumIntPoints_; j++)
      {
        // in 2D run the k loop one time for every i,j iteration
        for(UInt k = 0; k < (Dim_ == 3 ? line3->NumIntPoints_ : 1); k++)
        {
          // create and store our data
          row = new Double[Dim_ + 1];
          (*data)[counter] = row;
          counter++;
                
          // store the coordinates and the weight
          *(row + 0)= line1->IntPoints_[i][0];
          *(row + 1)= line2->IntPoints_[j][0];                
          if(Dim_ == 3) *(row + 2)= line3->IntPoints_[k][0];                
          *(row + Dim_)= line1->IntWeights_[i] * line2->IntWeights_[j] * (Dim_  == 3 ? line3->IntWeights_[k] : 1.);
        }
      }
    }

    IntegrationPointsMap_[key] = data;  
  }

  /** Expicitly set and load the integration type  */
  void BaseFE::SetIntPoints(IntegrationMethod method, int order)
  {
    IntegMethod = method;
    IntegOrder  = order;

    SetIntPoints();

    SetShapeFncAtIp();
    SetShapeFncDerivAtIp();  
  }


  
  // reads from the map but generates cartesian integration points on the 
  // fly when set in XML for the proper elements!
  void BaseFE::SetIntPoints()
  {
      
    // if we are not a valid element for CARTESIAN (assuming reading from XML)
    // the fallback will work
    if(IntegMethod == CARTESIAN 
       && (strcmp(GetShapeName(), ELEM_TYPE_NAMES[ET_HEXA8].c_str()) == 0  
           || strcmp(GetShapeName(), ELEM_TYPE_NAMES[ET_QUAD4].c_str()) == 0))
    {        
      // create the map entry on the fly
      int order1, order2, order3;
      DecodeCartesianOrder(IntegOrder, &order1, &order2, &order3);
         
      // this is the cached version, true to omit calling this method recursively!
      SetCartesianInteg(order1, order2, order3, true);
    }
      
    // searched upwards, e.g. a quadrilateral for order 2,3 should be 
    // stored with higher order 3, so we find for 2
    StdVector<Double*>* data = GetIntegrationPoints(IntegMethod, IntegOrder, true, false, true);
      
    NumIntPoints_= data->GetSize();

    if(IntPoints_ == NULL) {
      IntPoints_ = new Vector<Double>[NumIntPoints_];
    } else {
      delete[] IntPoints_;
      IntPoints_ = new Vector<Double>[NumIntPoints_];
    }
      
    IntWeights_.Resize(NumIntPoints_);
      
    for(UInt ip=0; ip<NumIntPoints_; ip++)
    {         
      IntPoints_[ip].Resize(Dim_);      	
          
      for(UInt dim = 0; dim < Dim_; dim++)
      {
        IntPoints_[ip][dim] = *((*data)[ip]+dim);
      }
          
      IntWeights_[ip]= *((*data)[ip]+Dim_);
    }
      
    // DumpIntegrationPoints();
  }
  
  void BaseFE::DumpIntegrationPointsMap()
  {
    std::string key;
    MakeKey(IntegMethod, IntegOrder, key);
        
    std::cout << "Current Integration method with " << key << " and " 
              << NumIntPoints_ << " integration points\n";
      
    std::map<const std::string, StdVector<Double*>*>::iterator iter;

    for(iter = IntegrationPointsMap_.begin(); iter != IntegrationPointsMap_.end(); iter++) {
      std::cout << iter->first << std::endl;
    }                    
    std::cout << std::flush;
  }  
 
  const char* BaseFE::GetShapeName() 
  {
    // the string representation of the FEType (environment.cc) does not
    // always match the names used on the XML file for integration type :(

    return ELEM_TYPE_NAMES[feType()].c_str();
  }
  
  // =======================================================================
  // L E G E N D R E    P A R T
  // =======================================================================
  void  BaseFE::GetNumFncs(Vector<UInt>& numFcns,  
                           const shared_ptr<AnsatzFct>& fcnType, 
                           AnsatzFct::FctEntityType fctEntityType, 
                           UInt dof) {

    
    // Check ansatzFctType
    if( fcnType->GetType() == AnsatzFct::LAGRANGE ) {
      numFcns.Resize( NumNodes_ );
      numFcns.Init(1);
    } else {
      EXCEPTION("In base class only implemented for Lagrange functions!");
    }
  }
  
  
  UInt BaseFE::GetNumFncs( const shared_ptr<AnsatzFct>& fcnType ) {

    // Check ansatzFctType
    if( fcnType->GetType() == AnsatzFct::LAGRANGE ) {
      return GetNumNodes();
    } else {
      EXCEPTION("In base class only implemented for Lagrange functions!");
    }

    return 0;
    
  }
  
  
  void BaseFE::EvalPolynom( Double& value, Double& deriv,
                            const UInt order, const Double* coeff, 
                            const Double xVal ) {

    // Consider the following expression
    // f(xVal) = a0 * (a1*x^order + a2*x^(order-1) + .. + a(order+1))
    // The coefficients a0..a(order+1) are stored in the coeff-array

    value = coeff[1];
    deriv = 0.0;
    for( UInt i = 2; i < order+2; i++ ) {
      deriv = deriv * xVal + value;
      value = value * xVal + coeff[i];
    }
    // Multiply by pre-factor
    deriv *= coeff[0];
    value *= coeff[0];
  }
  
  
  // Define coefficients for legendre ansatz functions up to order 8
  Double  BaseFE::lCoeff_[9][10] = {
    {0.5                  ,   -1, 1,    0, 0,   0, 0,    0, 0, 0 },
    {0.5                  ,    1, 1,    0, 0,   0, 0,    0, 0, 0 },
    {0.25*sqrt(6.0)       ,    1, 0,   -1, 0,   0, 0,    0, 0, 0 },
    {0.25*sqrt(10.0)      ,    1, 0,   -1, 0,   0, 0,    0, 0, 0 },
    {1.0/16.0*sqrt(14.0)  ,    5, 0,   -6, 0,   1, 0,    0, 0, 0 },
    {3.0/16.0*sqrt(2.0)   ,    7, 0,  -10, 0,   3, 0,    0, 0, 0 },
    {1.0/32.0*sqrt(22.0)  ,   21, 0,  -35, 0,  15, 0,   -1, 0, 0 },
    {1.0/32.0*sqrt(26.0)  ,   33, 0,  -63, 0,  35, 0,   -5, 0, 0 },
    {1.0/256.0*sqrt(30.0) ,  429, 0, -924, 0, 630, 0, -140, 0, 5 }
  };
  
} // end namespace CoupledField
