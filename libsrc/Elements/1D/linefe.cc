#include <iostream>
#include <fstream>
#include <string>

#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "DataInOut/WriteInfo.hh"
#include "General/environment.hh"
#include "linefe.hh"

namespace CoupledField
{


  LineFE :: LineFE()
  {
    ENTER_FCN( "LineFE::LineFE" );

    Dim_ = 1;
    NumEdges_   = 1;
    NumFaces_   = 1;
    NumCorners_ = 2;
    numChilds_  = 1;
    MidPoint_.Resize(2);
    MidPoint_[0] = 0.0;
    MidPoint_[0] = 0.0;
  }


  LineFE :: ~LineFE()
  {
    ENTER_FCN( "LineFE::~LineFE" );
  }


  void LineFE::FillIntegrationPoints()
  {
      ENTER_IFCN("LineFE::FillIntegrationPoints");  	
      
   	  //  Gauss  quadrature  points  and  weights  order  1 - Fabian
      static Double a1[][2] = { 
      	{0, 1}
      };
   	  AddIntegrationPoints(ECONOMICAL, 1, 1, (Double*) a1);
  	
   	
   	  //  Gauss  quadrature  points  and  weights  order  2n-1=3
      static Double a3[][2] = { 
      	{-0.5773502691896257645091488,  1.0000000000000000000000000},
      	{0.5773502691896257645091488,  1.0000000000000000000000000}
      };
   	  AddIntegrationPoints(ECONOMICAL, 3, 2, (Double*) a3);

 	  //  Gauss  quadrature  points  and  weights  order  2n-1=5
      static Double a5[][2] = { 
        {-0.7745966692414833770358531,  0.5555555555555555555555556},
        {0.0000000000000000000000000,  0.8888888888888888888888889},
        {0.7745966692414833770358531,  0.5555555555555555555555556}
      };
   	  AddIntegrationPoints(ECONOMICAL, 5, 3, (Double*) a5);
   	
 	  //  Gauss  quadrature  points  and  weights  order 2n-1=7
      static Double a7[][2] = { 
        {-0.8611363115940525752239465,  0.3478548451374538573730639},
        {-0.3399810435848562648026658,  0.6521451548625461426269361},
        {0.3399810435848562648026658, 0.6521451548625461426269361},
        {0.8611363115940525752239465,  0.3478548451374538573730639}
      };
   	  AddIntegrationPoints(ECONOMICAL, 7, 4, (Double*) a7);

 	  //  Gauss  quadrature  points  and  weights  order 2n-1=9
      static Double a9[][2] = { 
        {-0.9061798459386639927976269,  0.2369268850561890875142640},
        {-0.5384693101056830910363144,  0.4786286704993664680412915},
        {0.0000000000000000000000000, 0.5688888888888888888888889},
        {0.5384693101056830910363144,  0.4786286704993664680412915},
        {0.9061798459386639927976269,  0.2369268850561890875142640}
      };
   	  AddIntegrationPoints(ECONOMICAL, 9, 5, (Double*) a9);

 	  //  Gauss  quadrature  points  and  weights  order 2n-1=11
      static Double a11[][2] = { 
        {-0.9324695142031520278123016,  0.1713244923791703450402961},
        {-0.6612093864662645136613996,  0.3607615730481386075698335},
        {-0.2386191860831969086305017,  0.4679139345726910473898703},
        {0.2386191860831969086305017, 0.4679139345726910473898703},
        {0.6612093864662645136613996,  0.3607615730481386075698335},
        {0.9324695142031520278123016,  0.1713244923791703450402961}
      };
   	  AddIntegrationPoints(ECONOMICAL, 11, 6, (Double*) a11);

        
      // The original values from SetIntPoints() just have less significatn digits and another order count
      static Double c2[][2] = { 
        {-0.57735026919,  1},
        {0.57735026919,   1}
      };
      AddIntegrationPoints(CLASSICAL, 2, 2, (Double*) c2);


      // Lobatto  (Radau)  quadrature  constants  order  2n-3=3
      // Solin, Segeth, Dolezel, Higher-Order Finit Element Methods, p. 224
      static Double l3[][2] = { 
        {1.000000000000000000000000,  0.3333333333333333333333333},
        {0.000000000000000000000000,  1.3333333333333333333333333},
        {-1.000000000000000000000000,  0.3333333333333333333333333}
      };
      AddIntegrationPoints(LOBATTO, 3, 3, (Double*) l3);

      // Lobatto  (Radau)  quadrature  constants  order  2n-3=5
      // Solin, Segeth, Dolezel, Higher-Order Finit Element Methods, p. 224
      static Double l5[][2] = { 
        {1.000000000000000000000000,  0.1666666666666666666666667},
        {0.4472135954999579392818347,  0.8333333333333333333333333},        
        {-0.4472135954999579392818347,  0.8333333333333333333333333},                
        {-1.000000000000000000000000,  0.1666666666666666666666667}
      };
      AddIntegrationPoints(LOBATTO, 5, 4, (Double*) l5);
      
      // Lobatto  (Radau)  quadrature  constants  order  2n-3=7
      // Solin, Segeth, Dolezel, Higher-Order Finit Element Methods, p. 224
      static Double l7[][2] = { 
        {1.000000000000000000000000,  0.1000000000000000000000000},
        {0.6546536707079771437982925,  0.5444444444444444444444444},        
        {0.0000000000000000000000000,  0.7111111111111111111111111},                
        {-0.6546536707079771437982925,  0.5444444444444444444444444},
        {-1.000000000000000000000000,  0.1000000000000000000000000}
      };
      AddIntegrationPoints(LOBATTO, 7, 5, (Double*) l7);

      // Lobatto  (Radau)  quadrature  constants  order  2n-3=9
      // Solin, Segeth, Dolezel, Higher-Order Finit Element Methods, p. 224
      static Double l9[][2] = { 
        {1.000000000000000000000000,  0.0666666666666666666666667},
        {0.7650553239294646928510030,  0.3784749562978469803166128},        
        {0.2852315164806450963141510,  0.5548583770354863530167205},                
        {-0.2852315164806450963141510,  0.5548583770354863530167205},
        {-1.7650553239294646928510030,  0.3784749562978469803166128},
        {-1.000000000000000000000000,  0.0666666666666666666666667}        
      };
      AddIntegrationPoints(LOBATTO, 9, 6, (Double*) l9);

      // Lobatto  (Radau)  quadrature  constants  order  2n-3=11
      // Solin, Segeth, Dolezel, Higher-Order Finit Element Methods, p. 224
      static Double l11[][2] = { 
        {1.000000000000000000000000,  0.0476190476190476190476190},
        {0.8302238962785669298720322,  0.2768260473615659480107004},        
        {0.4688487934707142138037719,  0.4317453812098626234178710},                
        {0.0000000000000000000000000,  0.4876190476190476190476190},                        
        {-0.4876190476190476190476190,  0.4317453812098626234178710},
        {-0.8302238962785669298720322,  0.2768260473615659480107004},
        {-1.000000000000000000000000,  0.0476190476190476190476190}        
      };
      AddIntegrationPoints(LOBATTO, 11, 7, (Double*) l11);

      // Lobatto  (Radau)  quadrature  constants  order  2n-3=13
      // Solin, Segeth, Dolezel, Higher-Order Finit Element Methods, p. 224
      static Double l13[][2] = { 
        {1.000000000000000000000000,  0.0357142857142857142857143},
        {0.8717401485096066153374457,  0.2107042271435060393829921},        
        {0.5917001814331423021445107,  0.3411226924835043647642407},                
        {0.2092992179024788687686573,  0.4124587946587038815670530},                        
        {-0.2092992179024788687686573,  0.4124587946587038815670530},
        {-0.5917001814331423021445107,  0.3411226924835043647642407},
        {-0.8717401485096066153374457,  0.2107042271435060393829921},        
        {-1.000000000000000000000000,  0.0357142857142857142857143}        
      };
      AddIntegrationPoints(LOBATTO, 13, 8, (Double*) l13);
      
      // Chebyshev  quadrature  constants  order  n+1=3
      // Solin, Segeth, Dolezel, Higher-Order Finit Element Methods, p. 222
      static Double ch3[][2] = { 
        {0.5773502691896257645091488,  1.},
        {-0.5773502691896257645091488,  1.}
      };
      AddIntegrationPoints(CHEBYSHEV, 3, 2, (Double*) ch3);
      
      // Chebyshev  quadrature  constants  order  n+1=4
      // Solin, Segeth, Dolezel, Higher-Order Finit Element Methods, p. 222
      static Double ch4[][2] = { 
        {0.7071067811865475244008444,  2./3.},
        {0.0000000000000000000000000,  2./3.},
        {-0.7071067811865475244008444,  2./3.}      
      };
      AddIntegrationPoints(CHEBYSHEV, 4, 3, (Double*) ch4);


      // Chebyshev  quadrature  constants  order  n+1=5
      // Solin, Segeth, Dolezel, Higher-Order Finit Element Methods, p. 222
      static Double ch5[][2] = { 
        {0.7946544722917661229555309,  1./2.},
        {0.1875924740850798998601393,  1./2.},
        {-0.1875924740850798998601393,  1./2.},
        {-0.7946544722917661229555309,  1./2.}        
      };
      AddIntegrationPoints(CHEBYSHEV, 5, 4, (Double*) ch5);


      // Chebyshev  quadrature  constants  order  n+1=6
      // Solin, Segeth, Dolezel, Higher-Order Finit Element Methods, p. 222
      static Double ch6[][2] = { 
        {0.8324974870009818758925836,  2./5.},
        {0.3745414095535810655860444,  2./5.},
        {0.0000000000000000000000000,  2./5.},
        {-0.3745414095535810655860444,  2./5.},
        {-0.8324974870009818758925836,  2./5.}        
      };
      AddIntegrationPoints(CHEBYSHEV, 6, 5, (Double*) ch6);

      // Chebyshev  quadrature  constants  order  n+1=7
      // Solin, Segeth, Dolezel, Higher-Order Finit Element Methods, p. 222
      static Double ch7[][2] = { 
        {0.8662468181078205913835981,  1./3.},
        {0.4225186537611115291185464,  1./3.},        
        {0.2666354015167047203315346,  1./3.},        
        {-0.2666354015167047203315346,  1./3.},       
        {-0.4225186537611115291185464,  1./3.},                         
        {-0.8662468181078205913835981,  1./3.}        
      };
      AddIntegrationPoints(CHEBYSHEV, 7, 6, (Double*) ch7);

      // Chebyshev  quadrature  constants  order  n+1=8
      // Solin, Segeth, Dolezel, Higher-Order Finit Element Methods, p. 222
      static Double ch8[][2] = { 
        {0.8838617007580490357042241,  2./7.},
        {0.5296567752851568113850475,  2./7.},        
        {0.3239118105199076375196731,  2./7.},                
        {0.0000000000000000000000000,  2./7.},                        
        {-0.3239118105199076375196731,  2./7.},       
        {-0.5296567752851568113850475,  2./7.},                                         
        {-0.8838617007580490357042241,  2./7.}        
      };
      AddIntegrationPoints(CHEBYSHEV, 8, 7, (Double*) ch8);

      // Chebyshev  quadrature  constants  order  n+1=8
      // Solin, Segeth, Dolezel, Higher-Order Finit Element Methods, p. 222
      static Double ch10[][2] = { 
        {0.9115893077284344736649486,  2./9.},
        {0.6010186553802380714281279,  2./9.},
        {0.5287617830578799932601816,  2./9.},
        {0.1679061842148039430680319,  2./9.},
        {0.0000000000000000000000000,  2./9.},
        {-0.1679061842148039430680319,  2./9.},
        {-0.5287617830578799932601816,  2./9.},        
        {-0.6010186553802380714281279,  2./9.},        
        {-0.9115893077284344736649486,  2./9.}        
      };
      AddIntegrationPoints(CHEBYSHEV, 10, 9, (Double*) ch10);
      
   }   



  void LineFE::SetShapeFncAtIp()
  {
    ENTER_FCN( "LineFE::SetShapeFncAtIp" );

    if (!ShFncAtIp_)
      ShFncAtIp_ = new Vector<Double>[NumIntPoints_];
  
    for( UInt i=0; i<NumIntPoints_; i++ )
      CalcShapeFnc( ShFncAtIp_[i], IntPoints_[i]);
  
  }


  void LineFE::SetShapeFncDerivAtIp()
  {
    ENTER_FCN( "LineFE::SetShapeFncDerivAtIp" );

    if( !ShFncDerivAtIp_)
      ShFncDerivAtIp_ = new Matrix<Double>[NumIntPoints_];
  
    for( UInt i=0; i<NumIntPoints_; i++ )
      CalcLocalDerivShapeFnc( ShFncDerivAtIp_[i], IntPoints_[i]);
  }


  Double LineFE::CalcJacobianDet(const Vector<Double> & LCoord,
                                 const Matrix<Double> & CornerCoords)
  {
    ENTER_FCN( "LineFE::CalcJacobianDet" );

    Matrix<Double> J;
  
    CalcJacobian( J, LCoord, CornerCoords );
    return J[0][0];
  
  }


  Double LineFE::CalcJacobianDetAtIp(const UInt ip, 
                                     const Matrix<Double> & CornerCoords)
  {
    ENTER_FCN( "LineFE::CalcJacobianDetAtIp" );

    //   Matrix<Double> J;

    //   CalcJacobianAtIp( J, ip, CornerCoords);
    //   return J[0][0];

    Double length;
    if (CornerCoords.GetSizeRow()==2) 
      {
        //see kaltenbacher, p.23, eq.(2.122)
        Matrix<Double> J;
        CalcJacobianAtIp( J, ip, CornerCoords);
        length = sqrt(J[0][0]*J[0][0] + J[1][0]*J[1][0]);
      }
    else 
      {
        // Length/2 is simply the jacDet for a line
        length=dist_Mat(CornerCoords)/2;
      }

    return length;

  
  }

  void LineFE:: CalcJacobian(Matrix<Double> & J, 
                             const Vector<Double> & LCoord, 
                             const Matrix<Double> & CornerCoords)
  {
    ENTER_FCN( "LineFE::CalcJacobian" );

    J.Resize(1,1);

    Matrix<Double> LDeriv;
  
    CalcLocalDerivShapeFnc(LDeriv, LCoord);
    J = CornerCoords * LDeriv;
  }


  void LineFE::CalcJacobianAtIp(Matrix<Double> & J, 
                                const UInt ip, 
                                const Matrix<Double> & CornerCoords)
  {
    ENTER_FCN( "LineFE::CalcJacobianAtIp" );

    if (CornerCoords.GetSizeRow()==2) {
      // Surface element in 2D 
      J.Resize(CornerCoords.GetSizeRow(),Dim_);
    } 
    else
      J.Resize(1,1);

    J = CornerCoords * ShFncDerivAtIp_[ip-1];
  }


  void LineFE::CalcInvJacobian(Matrix<Double> & JInv,
                               const Vector<Double> & LCoord,
                               const Matrix<Double> & CornerCoords)
  {
    ENTER_FCN( "LineFE::CalcInvJacobian" );
  
    JInv.Resize(1,1);

    Matrix<Double> J, LDeriv;

    J.Resize(1,1);

    CalcLocalDerivShapeFnc(LDeriv, LCoord);
    J = CornerCoords * LDeriv;

    JInv[0][0] = 1 / J[0][0];
  }

  void LineFE::CalcInvJacobianAtIp(Matrix<Double> & JInv,
                                   const UInt ip,
                                   const Matrix<Double> & CornerCoords)
  {
    ENTER_FCN( "LineFE::CalcInvJacobianAtIp" );
  
    JInv.Resize(1,1);

    Matrix<Double> J;

    J.Resize(1,1);

    J = CornerCoords * ShFncDerivAtIp_[ip-1];

    JInv[0][0] = 1 / J[0][0];
  }

} // end of namespace
