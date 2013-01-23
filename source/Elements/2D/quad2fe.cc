// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <fstream>

#include "Elements/2D/rectanglefe.hh"
#include "General/exception.hh"
#include "MatVec/matrix.hh"
#include "quad2fe.hh"

namespace CoupledField
{
  const UInt Quad2FE::nodeNumbering_[8] = { 0, 4, 1, 5, 2, 6, 3, 7};


  Quad2FE::Quad2FE() : RectangleFE()
  {

    Init();
  }

  Quad2FE::~Quad2FE()
  {
  }

  void Quad2FE::Init()
  {
    NumNodes_ = 8;

    CommonInit();
    SetShapeFnc2ndDerivAtIp();
  }
  // Should be called SetNodalCoords!!
  void Quad2FE::SetCornerCoords()
  {

    LCornerCoords_.Resize(Dim_,NumNodes_);

    LCornerCoords_[0][0] = -1;
    LCornerCoords_[1][0] = -1;
    LCornerCoords_[0][1] =  1;
    LCornerCoords_[1][1] = -1;
    LCornerCoords_[0][2] =  1;
    LCornerCoords_[1][2] =  1;
    LCornerCoords_[0][3] = -1;
    LCornerCoords_[1][3] =  1;

    LCornerCoords_[0][4] =  0;
    LCornerCoords_[1][4] = -1;
    LCornerCoords_[0][5] =  1;
    LCornerCoords_[1][5] =  0;
    LCornerCoords_[0][6] =  0;
    LCornerCoords_[1][6] =  1;
    LCornerCoords_[0][7] = -1;
    LCornerCoords_[1][7] =  0;

  }

  void Quad2FE :: CalcShapeFnc(Vector<Double> & Shape,
                               const Vector<Double> & LCoord,
                               const Elem*, UInt dof,
                               AnsatzFct::FctEntityType )
  {

    Shape.Resize(NumNodes_);
    // From Zienkiewicz, The Finite Element Method. Vol 1, page 122.
    // corner nodes
    for( UInt i=0; i<(NumNodes_/2); i++)
      Shape[i] = 0.25 * (1 + LCornerCoords_[0][i] * LCoord[0])
        * (1 + LCornerCoords_[1][i] * LCoord[1])
        * (LCornerCoords_[0][i] * LCoord[0] +
           LCornerCoords_[1][i] * LCoord[1] - 1);
    // midside node
    for( UInt i=(NumNodes_/2); i<NumNodes_; i=i+2)
      {
        Shape[i]   = 0.5 * (1 - LCoord[0]*LCoord[0])*
          (1 + LCornerCoords_[1][i] * LCoord[1]);
        Shape[i+1] = 0.5 * (1 - LCoord[1]*LCoord[1])*
          (1 + LCornerCoords_[0][i+1] * LCoord[0]);
      }
  }


  void Quad2FE::CalcLocalDerivShapeFnc(Matrix<Double> & LDeriv,
                                         const Vector<Double> & LCoord,
                                         const Elem*, UInt dof,
                                         AnsatzFct::FctEntityType )
  {

    LDeriv.Resize(NumNodes_,Dim_);

    // corner nodes
    for( UInt i=0; i<(NumNodes_/2); i++)
      {
        LDeriv[i][0] = 0.25 * LCornerCoords_[0][i] *
          (1 + LCornerCoords_[1][i] * LCoord[1]) *
          (2 * LCornerCoords_[0][i] * LCoord[0] +
           LCornerCoords_[1][i] * LCoord[1]);

        LDeriv[i][1] = 0.25 * LCornerCoords_[1][i] *
          (1 + LCornerCoords_[0][i] * LCoord[0]) *
          (2 * LCornerCoords_[1][i] * LCoord[1] +
           LCornerCoords_[0][i] * LCoord[0]);
      }

    // midside node
    for( UInt i=(NumNodes_/2); i<NumNodes_; i=i+2)
      {
        LDeriv[i][0]   = - LCoord[0]*
          (1 + LCornerCoords_[1][i] * LCoord[1]);

        LDeriv[i][1]   = 0.5 * LCornerCoords_[1][i] *
          (1 - LCoord[0]*LCoord[0]);

        LDeriv[i+1][0] = 0.5 *  LCornerCoords_[0][i+1] *
          (1 - LCoord[1]*LCoord[1]);

        LDeriv[i+1][1] = -LCoord[1]*
          (1 + LCornerCoords_[0][i+1] * LCoord[0]);
      }

  }

  void Quad2FE :: CalcLocal2ndDerivShapeFnc(Matrix<Double> & L2ndDeriv, 
                                            const Vector<Double> & LCoord,
                                            const Elem* elem,
                                            UInt dof,
                                            AnsatzFct::FctEntityType)
  {

    if (Dim_==2) {
      L2ndDeriv.Resize(NumNodes_,Dim_+1);

      // corner nodes
      for( UInt i=0; i<(NumNodes_/2); i++)
        {
          L2ndDeriv[i][0] = 0.5 * LCornerCoords_[0][i] * LCornerCoords_[0][i] * 
            ( 1 + LCornerCoords_[1][i] * LCoord[1] );

          L2ndDeriv[i][1] = 0.5 * LCornerCoords_[1][i] * LCornerCoords_[1][i] * 
            ( 1 + LCornerCoords_[0][i] * LCoord[0] );

          L2ndDeriv[i][2] = 0.25 * LCornerCoords_[0][i] * LCornerCoords_[1][i] * 
            ( (2 * LCornerCoords_[0][i] * LCoord[0]) + 
              (2 * LCornerCoords_[1][i] * LCoord[1]) + 1);

        }
  
      // midside node
      for( UInt i=(NumNodes_/2); i<NumNodes_; i+=2)
        {


          L2ndDeriv[i][0]   = - (1 + LCornerCoords_[1][i] * LCoord[1]);

          L2ndDeriv[i][1]   = 0.0;

          L2ndDeriv[i][2]   = - LCornerCoords_[1][i] * LCoord[0];


          L2ndDeriv[i+1][0] = 0.0;

          L2ndDeriv[i+1][1] = - (1 + LCornerCoords_[0][i+1] * LCoord[0]);

          L2ndDeriv[i+1][2] = - LCornerCoords_[0][i+1] * LCoord[1];
        }
    }
  }


  Double Quad2FE::CalcMeanStrain(Matrix<Double> &cornerCoords, Matrix<Double> &displacements)
  {

    EXCEPTION( "Quad2FE::CalcDistortion. This function has not yet been implemented in "
             << " quadratic quadrilaterals" );
    return -1;

    //   Double factor;
    //   Double eps1, eps2, eps4, eps5, eps11, eps12, eps21, eps22, eps41, eps42, eps51, eps52;
    //   Double length1, length2, length11, length12, length21, length22;


    //   // calculate inital size of element
    //   length11 = abs(cornerCoords[0][0] - cornerCoords[0][1]);
    //   length12 = abs(cornerCoords[0][3] - cornerCoords[0][2]);
    //   length1 = (length11+length12) * 0.5;

    //   length21 = abs(cornerCoords[1][0] - cornerCoords[1][3]);
    //   length22 = abs(cornerCoords[1][1] - cornerCoords[1][2]);
    //   length2 = (length21+length22) * 0.5;

    //    // calculate absolute change of size
    //   eps11 = displacements[0][0] - displacements[0][1];
    //   eps12 = displacements[0][3] - displacements[0][2];
    //   eps1 = (eps11+eps12) * 0.5;

    //   eps21 = displacements[1][0] - displacements[1][3];
    //   eps22 = displacements[1][1] - displacements[1][2];
    //   eps2 = (eps21+eps22) * 0.5;

    //   eps41 = displacements[0][2] - displacements[0][1];
    //   eps42 = displacements[0][3] - displacements[0][0];
    //   eps4 = (eps41+eps42)*0.5;

    //   eps51 = displacements[1][1] - displacements[1][0];
    //   eps52 = displacements[1][3] - displacements[1][2];
    //   eps5= (eps51+eps52)*0.5;

    //   factor =  0.2 * ((eps1*eps1/(length1*length1))
    //               + (eps2*eps2/(length2*length2))
    //               + (eps5*eps5/(length1*length2))
    //               + (eps4*eps4/(length1*length1)));

    //   return factor;


  }

  void Quad2FE::fastGlobal2LocalCoords(Matrix<Double> & localCoords,
                                   const Matrix<Double> & globalCoords,
                                   const Matrix<Double> & coordMat ){
    //BaseFE::Global2LocalCoords(localCoords,globalCoords,coordMat);
    //return;
    //Version 0 algorithm from duester script
    Vector<Double> delta_xi; // update for Newton-Raphson method
    Vector<Double> xi_start; // local start point for Newton-Raphson method
    Vector<Double> xi_k; // local point at iteration k
    Vector<Double> f; //current right hand side
    Vector<Double> f_start; //current right hand side
    Vector<Double> globalPoint; // global point coordinates
    UInt numPoints = globalCoords.GetNumCols(); // number of global points
    Double f_old; //storing the absolute value of search direction
    Double f_test; //storing the absolute value of search direction
    UInt k = 0; //iteration counter
    Integer l = 0; //stepping value
    Double jacDet = 0;
    Matrix<Double> J; // Jacobian at local point xi_k
    localCoords.Resize(2, numPoints);
    localCoords.Init();

    Double tolerance = 1e-14;

    //initialize everything
    xi_start.Resize(2);
    delta_xi.Resize(2);
    xi_k.Resize(2);
    J.Resize(2, 2);
    globalPoint.Resize(2);

    //first initial guess is zero
    f.Resize(2);

   // Perform Newton-Raphson method on the list of global points
   // to find local coordinates within this element.
   for(UInt i = 0; i < numPoints; i++)
   {
     f.Init();
     globalPoint.Init();
     J.Init();
     xi_start.Init();
     xi_k.Init();
     k = 0;
     f_old = 222e20;
     f_test = 0;
     for(UInt j = 0; j < 2; j++) {
       globalPoint[j] = globalCoords[j][i];
     }
     Local2GlobalCoord(f, xi_k, coordMat, NULL);
     f = f - globalPoint;
     f_old = f.NormL2();
     f_start = f;
     for(Double w=0.5; w <= 1.0; w+=0.5){
       for(UInt j=1;j<3;j++){
         for(UInt m=1;m<3;m++){
           xi_k[0] = pow(-1,j)*w;
           xi_k[1] = pow(-1,m)*w;
           Local2GlobalCoord(f, xi_k, coordMat, NULL);
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
     while(f_test > tolerance && k < 1000){
        delta_xi.Init();
        xi_start.Init();
        // Calculate Jacobian at iteration point xi_k
        CalcJacobian(J, xi_k, coordMat, NULL );
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
           Local2GlobalCoord(f, xi_start, coordMat, NULL);
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
     if( f_test > tolerance){
       std::cout << "performed " << k << " iterations to reach the point" << std::endl<< xi_k << std::endl;
       //ONLY for DEBUGGING
       Local2GlobalCoord(f, xi_k, coordMat, NULL);
       std::cout << "Calculated global point :" << std::endl << f << std::endl;
       std::cout << "Original global coord " << std::endl << globalPoint << std::endl;
       std::cout << "The error was: " << f_test << std::endl<< std::endl;
       //std::cout << std::endl;
     }

     // Put local coordinate of point into matrix.
     for(UInt j = 0; j < 2; j++) {
       if(xi_k.GetSize() == 0){
         std::cerr << "local2global messed up setting everything to zero" << std::endl;
         std::cerr << globalCoords << std::endl;
         xi_k.Resize(2);
         xi_k.Init();
       }
       localCoords[j][i] = xi_k[j];
     }
   }
 }

} // end of namespace


