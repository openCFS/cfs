#include <fstream>
#include <iostream>
#include <string>
#include <stdio.h>

#include <cmath>     


#include "Gilmoredimlos.hh"


namespace CoupledField
{


  Gilmoredimlos::Gilmoredimlos(Double  RadiusInit, 
			       Double  density,
			       Double  sonicVel,
			       Double  pStatic,
			       Double  pVapour,
			       Double  surfacTen,
			       Double  polytrop,
			       Double  viscosity) {
    
    ENTER_FCN( "Gilmoredimlos::Gilmoredimlos" );

    RadiusInit_ = RadiusInit;
    density_    = density;
    sonicVel_   = sonicVel;
    pStatic_    = pStatic;
    pVapour_    = pVapour;
    surfacTen_  = surfacTen;
    polytrop_   = polytrop;
    viscosity_  = viscosity;

    n_ = 7.0;
    A_ = 3.001e8; 
    B_ = 3e8;

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //Compute dimensionless constants needed in CompDeriv

    sigma_ = 2 * surfacTen_ / (RadiusInit_ * pStatic_);
    c0_    = sonicVel_ / (sqrt( pStatic_ / density_ ));
    mu_    = viscosity_ / ( pStatic_ * RadiusInit_) * sqrt( pStatic_ / density_); 
    a_     = A_ / pStatic_;
    b_     = B_ / pStatic_;
    pV_    = pVapour_ / pStatic_;


  }




  void  Gilmoredimlos::CompDeriv(const Double &t,
				 const StdVector<Double> &y,
				 StdVector<Double> &dydt){
    ENTER_FCN( "Gilmoredimlos::CompDeriv" );


    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Dimensionless computation of the right hand sight of the 
    // Gilmore equation y=dydt(y);
    // Input and output values of y und dydt are dimensionless,
    // therefore choose right call in bubbledriver.cc and 
    // solveStepAcousticBubble.cc
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    // For testing only:-------------------------------------------

//     Double tdim;
//     StdVector<Double> ydim(2);
//     Double ampl, freq;


//     ydim[0] = y[0] * RadiusInit_ ;
//     ydim[1] = y[1] * sqrt(pStatic_ / density_); 
//     tdim = t * RadiusInit_ / (sqrt( pStatic_ / density_));
//     printf("%10.6e \t %10.6e \t %10.6e \n",tdim, ydim[0],ydim[1]);

//     ampl  = 24.0e5;            //Anregungsamplitude
//     freq  = 27000.0;           //Anregungsfrequenz
//     ampl  = ampl / pStatic_;
//     freq  = freq *  RadiusInit_ / (sqrt( pStatic_ / density_));
//     p_    = -ampl *sin(2 * (PI) *freq * (t));                        
//     dpdt_ = -ampl *(2 *(PI) * freq)* cos(2 * (PI) *freq * (t));
    // -------------------------------------------------------------


    // Dimensionless Gilmore model 




    Double temp1, temp2, temp3, temp4, temp5, temp6, temp7;

    dydt[0] = y[1];
    
    temp1  = (1 - pV_ + sigma_) / std::pow((y[0]),(3.0*polytrop_));
    temp1 += (- sigma_) / y[0] - 4.0 * mu_ * y[1] / y[0] + pV_ + b_;
    temp2  = std::pow(temp1,((n_-1.0)/n_));
    temp3  = std::pow((1 + p_ + b_),((n_-1.0)/n_));
    H_      = (temp2 - temp3) * std::pow(a_,(1.0/n_)) * n_/(n_-1.0);
    C_      = sqrt(c0_ * c0_ + (n_-1.0) * H_);


    temp4  = (1 - pV_ + sigma_ ) /std::pow((y[0]),(3.0*polytrop_))
      *(-3.0) * polytrop_ * y[1] / y[0];
    temp4 += sigma_ / (y[0] * y[0]) *y[1] 
      + 4.0* mu_ *(y[1] * y[1]) / (y[0] * y[0]);
    temp5  = temp4 * std::pow(temp1,(-1.0/n_));
    temp6  = std::pow((1 + p_ + b_),(-1.0/n_)) * dpdt_;

    dydt[1]= (temp5 - temp6) * std::pow(a_,(1.0/n_))
      * (1.0-(y[1]/C_)) * y[0]/ C_;
    dydt[1]+= (1.0+ (y[1]/C_)) * H_ 
      - 3.0 / 2.0 * y[1] * y[1]*(1.0 - (y[1]/(3.0*C_)));
    
    temp7  = (1.0 - (y[1] / C_)) * y[0]* 
      (1.0+( 4.0 * mu_ / (C_*y[0]) *std::pow(temp1,(-1.0/n_)) 
	     * std::pow(a_,(1.0/n_))));
    dydt[1]= dydt[1] / temp7;

  }


  // +++++++++++++++++++++++++++++++++++++++++++++++++++++
  // Choose only one of the following methods
  // to obtain the Jacobian Matrix

  // Function to compute the Jacobian
//   void  Gilmoredimlos::Jacobi(StdVector<Double> &y,
// 			      Matrix<Double> &dfdy,
// 			      Double &t){
    
//     ENTER_FCN( "Gilmoredimlos::Jacobi" );

//     dfdy(0,0) = 0;
//     dfdy(0,1) = 1;
//     dfdy(1,0) = df2dy1(y,t);
//     dfdy(1,1) = df2dy2(y,t);
// //    //    std::cout << dfdy << std::endl;
//   }

  // Second order function to approximate the Jacobian
  void  Gilmoredimlos::Jacobi(StdVector<Double> &y,
   			      Matrix<Double> &dfdy,
   			      Double &t){
  
    ENTER_FCN( "Gilmoredimlos::Jacobi" );
  
    Integer n = y.GetSize();
    StdVector<Double> yPlus(n);
    StdVector<Double> yMinus(n);
    StdVector<Double> dydtPlus(n);
    StdVector<Double> dydtMinus(n);
    Double Tiny=1e-30;
  
    yPlus  = y;
    yMinus = y;
    for (Integer i=0; i<n; i++){
      yPlus[i]  = 1.001 *( y[i]+Tiny);
      yMinus[i] = 0.999 *( y[i]+Tiny);
      CompDeriv(t,yPlus,dydtPlus);
      CompDeriv(t,yMinus,dydtMinus);
      for (Integer k=0; k<n; k++){
      	dfdy(k,i) = 0.5*(dydtPlus[k]-dydtMinus[k])/(yPlus[i]-y[i]);
      }
      yPlus[i]  = y[i];
      yMinus[i] = y[i];
    }
// //Nur für Vergleichszwecke der beiden Jacobiverfahren
// //      Matrix<Double> dfdy2(2,2);
// //      dfdy2(0,0) = 0;
// //      dfdy2(0,1) = 1;
// //      dfdy2(1,0) = df2dy1(y,t);
// //      dfdy2(1,1) = df2dy2(y,t);
// //      std::cout << dfdy2 << "    "<< dfdy<<std::endl;
  }

  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  // Help-function to determine the Jacobian,
  // constructed with the help of maple,
  // only needed in case of computation
  Double Gilmoredimlos::df2dy1(StdVector<Double> &y,
			       Double &t){


    ENTER_FCN( "Gilmoredimlos::df2dy1" );

    Double y1 = y[0];
    Double y2 = y[1];


    // For testing only:-------------------------------------------
    // Double ampl, freq;

    //ampl  = 24.0e5;            //Anregungsamplitude
    //freq  = 27000.0;           //Anregungsfrequenz
    //ampl  = ampl / pStatic_;
    //freq  = freq *  RadiusInit_ / (sqrt( pStatic_ / density_));
    //p_    = -ampl *sin(2 * (PI) *freq * (t));                        
    //dpdt_ = -ampl *(2 *(PI) * freq)* cos(2 * (PI) *freq * (t));
    // -------------------------------------------------------------

    Double t39;
    Double t78;
    Double t53;
    Double t26;
    Double t20;
    Double t67;
    Double t34;
    Double t100;
    Double t10;
    Double t59;
    Double t41;
    Double t88;
    Double t89;
    Double t17;
    Double t65;
    Double t80;
    Double t129;
    Double t21;
    Double t15;
    Double t99;
    Double t37;
    Double t38;
    Double t22;
    Double t23;
    Double t101;
    Double t28;
    Double t18;
    Double t19;
    Double t104;
    Double t11;
    Double t60;
    Double t6;
    Double t146;
    Double t147;
    Double t76;
    Double t77;
    Double t137;
    Double t136;
    Double t12;
    Double t13;
    Double t42;
    Double t43;
    Double t44;
    Double t8;
    Double t110;
    Double t3;
    Double t127;
    Double t25;
    Double t55;
    Double t62;
    Double t92;
    Double t35;
    Double t36;
    Double t158;
    Double t4;
    Double t5;
    Double t131;
    Double t132;
    Double t1;
    Double t48;




    t1 = y2 * y2;
    t3 = c0_ * c0_;
    t4 = 0.1e1 / n_;
    t5 = 0.10e1 * t4;
    t6 = std::pow(a_, t5);
    t8 = 0.1e1 - pV_ + sigma_;
    t10 = std::pow(y1, -0.3e1 * polytrop_);
    t11 = t8 * t10;
    t12 = 0.1e1 / y1;
    t13 = sigma_ * t12;
    t15 = mu_ * t12 * y2;
    t17 = t11 - t13 - 0.4e1 * t15 + pV_ + b_;
    t18 = n_ - 0.1e1;
    t19 = t18 * t4;
    t20 = std::pow(t17, t19);
    t21 = 0.1e1 + p_ + b_;
    t22 = std::pow(t21, t19);
    t23 = t20 - t22;
    t25 = t3 + t6 * n_ * t23;
    t26 = sqrt(t25);
    t28 = 0.1e1 / t26 / t25;
    t34 = 0.3e1 * t11 * polytrop_ * t12;
    t35 = y1 * y1;
    t36 = 0.1e1 / t35;
    t37 = sigma_ * t36;
    t38 = mu_ * t36;
    t39 = t38 * y2;
    t41 = -t34 + t37 + 0.4e1 * t39;
    t42 = 0.1e1 / t17;
    t43 = t41 * t42;
    t44 = t20 * t18 * t43;
    t48 = t6 * t6;
    t53 = 0.1e1 / (n_ - 0.10e1);
    t55 = n_ * t53 * t23;
    t59 = 0.1e1 / t26;
    t60 = y2 * t59;
    t62 = (0.1e1 + t60) * t6;
    t65 = 0.1e1 - t60;
    t67 = t8 * polytrop_;
    t76 = std::pow(t17, -t5);
    t77 = (-0.3e1 * t67 * t10 * y2 * t12 + t37 * y2 + 0.40e1 * t38 * t1) * t76;
    t78 = std::pow(t21, -t5);
    t80 = t77 - dpdt_ * t78;
    t88 = t18 * t41 * t42;
    t89 = t80 * t20 * t88;
    t92 = t25 * t25;
    t99 = y1 * t59;
    t100 = t65 * t6;
    t101 = polytrop_ * polytrop_;
    t104 = t10 * t36 * y2;
    t110 = 0.1e1 / t35 / y1;
    t127 = 0.1e1 / t65;
    t129 = t6 * mu_;
    t131 = t11 - t13 - 0.40e1 * t15 + pV_ + b_;
    t132 = std::pow(t131, -t5);
    t136 = y1 + 0.4e1 * t129 * t59 * t132;
    t137 = 0.1e1 / t136;
    t146 = -0.1500000000e1 * (0.1e1 - t60 / 0.3e1) * t1 + t62 * t55 + t99 * t100 * t80;
    t147 = t65 * t65;
    t158 = t136 * t136;
    return((-0.2500000000e0 * t1 * y2 * t28 * t6 * t44 - y2 * t28 * t48 * t20 * t18 * t43 * t55 / 0.2e1 + t62 * t53 * t44 + t59 * t65 * t6 * t80 - y1 * t28 * t65 * t48 * t89 / 0.2e1 + y1 / t92 * y2 * t48 * t89 / 0.2e1 + t99 * t100 * ((0.9e1 * t8 * t101 * t104 + 0.3e1 * t67 * t104 - 0.2e1 * sigma_ * t110 * y2 - 0.80e1 * mu_ * t110 * t1) * t76 - 0.10e1 * t77 * t4 * t41 * t42)) * t127 * t137 - t146 / t147 * t137 * y2 * t28 * t6 * t20 * t88 / 0.2e1 - t146 * t127 / t158 * (0.1e1 - 0.2e1 * t48 * mu_ * t28 * t132 * t44 - 0.40e1 * t129 * t59 * t132 * t4 * (-t34 + t37 + 0.40e1 * t39) / t131));

  }


  // Help-function to determine the Jacobian,
  // constructed with the help of maple 
  // only needed in case of computation
  Double Gilmoredimlos::df2dy2(StdVector<Double> &y,
			       Double &t){

    ENTER_FCN( "Gilmoredimlos::df2dy2" );

    Double y1 = y[0];
    Double y2 = y[1];

    // For testing only:-------------------------------------------
    //Double ampl, freq;

    //ampl  = 24.0e5;            //Anregungsamplitude
    //freq  = 27000.0;           //Anregungsfrequenz
    //ampl  = ampl / pStatic_;
    //freq  = freq *  RadiusInit_ / (sqrt( pStatic_ / density_));
    //p_    = -ampl *sin(2 * (PI) *freq * (t));                        
    //dpdt_ = -ampl *(2 *(PI) * freq)* cos(2 * (PI) *freq * (t));
    // -------------------------------------------------------------


    Double t48;
    Double t51;
    Double t23;
    Double t24;
    Double t25;
    Double t63;
    Double t17;
    Double t18;
    Double t19;
    Double t20;
    Double t21;
    Double t77;
    Double t78;
    Double t84;
    Double t42;
    Double t69;
    Double t70;
    Double t71;
    Double t6;
    Double t35;
    Double t36;
    Double t16;
    Double t3;
    Double t4;
    Double t2;
    Double t13;
    Double t15;
    Double t89;
    Double t90;
    Double t28;
    Double t61;
    Double t110;
    Double t114;
    Double t115;
    Double t119;
    Double t120;
    Double t73;
    Double t9;
    Double t10;
    Double t11;
    Double t8;
    Double t128;
    Double t134;
    Double t136;
    Double t33;
    Double t34;
    Double t39;
    Double t94;
    Double t56;
    Double t79;
    Double t81;
    Double t44;
    Double t1;
    Double t127;
    Double t53;
    Double t12;



    t1 = c0_ * c0_;
    t2 = 0.1e1 / n_;
    t3 = 0.10e1 * t2;
    t4 = std::pow(a_, t3);
    t6 = 0.1e1 - pV_ + sigma_;
    t8 = std::pow(y1, -0.3e1 * polytrop_);
    t9 = t6 * t8;
    t10 = 0.1e1 / y1;
    t11 = sigma_ * t10;
    t12 = mu_ * t10;
    t13 = t12 * y2;
    t15 = t9 - t11 - 0.4e1 * t13 + pV_ + b_;
    t16 = n_ - 0.1e1;
    t17 = t16 * t2;
    t18 = std::pow(t15, t17);
    t19 = 0.1e1 + p_ + b_;
    t20 = std::pow(t19, t17);
    t21 = t18 - t20;
    t23 = t1 + t4 * n_ * t21;
    t24 = sqrt(t23);
    t25 = 0.1e1 / t24;
    t28 = 0.1e1 / t24 / t23;
    t33 = 0.1e1 / t15;
    t34 = t10 * t33;
    t35 = t16 * mu_ * t34;
    t36 = y2 * t28 * t4 * t18 * t35;
    t39 = y2 * y2;
    t42 = y2 * t25;
    t44 = 0.1e1 - t42 / 0.3e1;
    t48 = t25 + 0.2e1 * t36;
    t51 = 0.1e1 / (n_ - 0.10e1);
    t53 = n_ * t51 * t21;
    t56 = (0.1e1 + t42) * t4;
    t61 = 0.1e1 - t42;
    t63 = t4 * t4;
    t69 = y1 * y1;
    t70 = 0.1e1 / t69;
    t71 = sigma_ * t70;
    t73 = mu_ * t70;
    t77 = std::pow(t15, -t3);
    t78 = (-0.3e1 * t6 * polytrop_ * t8 * y2 * t10 + t71 * y2 + 0.40e1 * t73 * t39) * t77;
    t79 = std::pow(t19, -t3);
    t81 = t78 - dpdt_ * t79;
    t84 = t18 * t16;
    t89 = y1 * t25;
    t90 = -t48;
    t94 = t61 * t4;
    t110 = 0.1e1 / t61;
    t114 = t9 - t11 - 0.40e1 * t13 + pV_ + b_;
    t115 = std::pow(t114, -t3);
    t119 = y1 + 0.4e1 * t4 * mu_ * t25 * t115;
    t120 = 0.1e1 / t119;
    t127 = -0.1500000000e1 * t44 * t39 + t56 * t53 + t89 * t94 * t81;
    t128 = t61 * t61;
    t134 = t119 * t119;
    t136 = mu_ * mu_;
    return((-0.1500000000e1 * (-t25 / 0.3e1 - 0.2e1 / 0.3e1 * t36) * t39 - 0.3000000000e1 * t44 * y2 + t48 * t4 * t53 - 0.4e1 * t56 * t51 * t18 * t35 + 0.2e1 * t28 * t61 * t63 * t81 * t84 * mu_ * t33 + t89 * t90 * t4 * t81 + t89 * t94 * ((-0.3e1 * t9 * polytrop_ * t10 + t71 + 0.80e1 * t73 * y2) * t77 + 0.40e1 * t78 * t2 * t12 * t33)) * t110 * t120 - t127 / t128 * t120 * t90 - t127 * t110 / t134 * (0.8e1 * t63 * t136 * t28 * t115 * t84 * t34 + 0.1600e2 * t4 * t136 * t25 * t115 * t2 * t10 / t114));

  }




 
} // end of nanespace

