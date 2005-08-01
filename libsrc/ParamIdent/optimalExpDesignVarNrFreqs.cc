#include <iostream>
#include <fstream>
#include <string>
#include "DataInOut/GMV/outGMV.hh"
#include "General/environment.hh"
#include "PDE/SinglePDE.hh" 

#include "piezoParamIdent.hh"
#include "Forms/baseForm.hh"
#include "Utils/vector.hh"
#include "Utils/nodestoresol.hh"
#include "Utils/elemstoresol.hh"
#include "DataInOut/MaterialData.hh"
#include "PDE/timestepping.hh"
#include "Utils/baseelemstoresol.hh"
#include "Driver/singleDriver.hh"
#include "PDE/nodeEQN.hh"
#include <Domain/elem.hh>
#include "Forms/forms_header.hh"
#include "Utils/mathfunctions.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"


#ifdef __sgi
#include <stdarg.h>
#include <stdio.h>
#include <math.h>
#define POW pow
#else
#include <cstdarg>
#include <cstdio>
#include <cmath>
#define POW std::pow
#endif

#include <stdlib.h>
#include <sstream>
#include <iomanip>


#include "Utils/tools.hh"
#include <PDE/pdes_header.hh>


namespace CoupledField
{

  // xxxxxxxxxxxxxxxxxxxxxxxxxx nuMethods  - ITERATION xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  void piezoParamIdent::optimalExpDesignVarNrFreqs(){

    ENTER_FCN("piezoParamIdent::optimalDesignVarNrFreqs");
    std::cout<<"call for optimum experiment desing with variable number of frequencies "<<std::endl;


    Double hOmega = (stopfreq-startfreq)/nrfreq;
    freqs.Resize(nrfreq);
    nrMeasuredData=nrfreq;
    
    for (UInt actFreq=0;actFreq<nrfreq;actFreq++)
      freqs[actFreq] = startfreq+actFreq*hOmega;

    std::cout<<freqs<<std::endl;

    UInt maxNrBreakpoints = nrfreq;

    Complex J;

    Double dOmega = 0.0001;
    rhos.Resize(nrfreq);

    for(UInt i=0;i<nrfreq;i++)
      rhos[i]=0.1;


    descentMethodRho(J);

    std::cout<<"Final rhos ..." <<std::endl;
    std::cout<<rhos <<std::endl;
    //  std::cout<< " J after descent - method = " << J <<std::endl;
 


  }

  void piezoParamIdent::createJRho(Complex &J, Boolean writeOutCov){
    ENTER_FCN("piezoParamIdent::createJRho");

    Double  hOmega = (stopfreq-startfreq)/nrfreq;

    //    Complex J;
    J=Complex(0.0,0.0);

    Vector<Double> newFreqs;
    Vector<Complex> jacobi;
    Vector<Complex> jacobiH;

    Double penaltyFactor1=0.1;
    Double penaltyFactor2=0.1;


    readInMeasurement(newFreqs);
    calc_measuredCharge(freqs, real, imag, y_hat); // out of new measurements

    fr_=4000/(2*thickness);
    // minus 5 percent
    fr_=0.94*fr_;
    // fa = 1/thickness * sqrt(c_33/rho)
   
    Double rho;
    rho = ptMaterial->GetDensity();

    fa_ = 1.0/(2*thickness) * std::sqrt(1.3e+11/rho);
    fa_=1.2*fa_;

//     std::cout<<"++ Bounds for resonance and antiresonace frequency"<<std::endl;
//     std::cout<<"fr:"<<std::endl;
//     std::cout<<fr_<<std::endl;

//     std::cout<<"fa:"<<std::endl;
//     std::cout<<fa_<<std::endl;

    //    MaterialData * ptMaterial=ptMyPDE_->getPDEMaterialData();   // Pointer to MaterialData
    Complex F_hat_incr1,F_hat_incr2;

    Matrix<Complex> cov(actNrParameter+actNrParameterC, actNrParameter+actNrParameterC);
    Matrix<Complex> covTemp(actNrParameter+actNrParameterC, actNrParameter+actNrParameterC);

     for (UInt actFreq=0;actFreq<nrfreq;actFreq++){
    
       createJacobian(jacobi, freqs[actFreq]);

       //       std::cout<<"Frequency = " << actFreq<<") "<< freqs[actFreq]<<std::endl;
//        std::cout<<"jacobi:"<<std::endl;
//        std::cout<<jacobi<<std::endl;

       jacobiH.Resize(actNrParameter+actNrParameterC);
       
       //       std::cout<<jacobi<<std::endl;
       for (UInt i=0;i<jacobi.GetSize();i++)
         jacobiH[i]=Complex(jacobi[i].real(),-jacobi[i].imag());
       
       for(UInt i=0;i<actNrParameter+actNrParameterC;i++)
         for(UInt j=0;j<actNrParameter+actNrParameterC;j++){
	   //           if (freqs[actFreq]<=fr_||freqs[actFreq]>=fa_){
         //           covTemp[i][j]=jacobiH[i]*jacobi[j]/(delta*delta*std::abs(y_hat[actFreq])*std::abs(y_hat[actFreq]));
             if(actFreq==0||actFreq==nrfreq)
               covTemp[i][j]=rhos[actFreq]*hOmega*0.5*jacobiH[i]*jacobi[j]/((1.0+delta)*(1.0+delta)*std::abs(y_hat[actFreq])*std::abs(y_hat[actFreq]));
             //covTemp[i][j]=rhos[actFreq]*hOmega*0.5*jacobiH[i]*jacobi[j]/((1.0+delta)*(1.0+delta)*std::abs(y_hat[actFreq]));
             else 
               covTemp[i][j]=rhos[actFreq]*hOmega*jacobiH[i]*jacobi[j]/((1.0+delta)*(1.0+delta)*std::abs(y_hat[actFreq])*std::abs(y_hat[actFreq]));
             //covTemp[i][j]=rhos[actFreq]*hOmega*jacobiH[i]*jacobi[j]/((1.0+delta)*(1.0+delta)*std::abs(y_hat[actFreq]));
             if (i==j)
               covTemp[i][i]+=covTemp[i][i];// *1.0e-10;
           }
       //}
       cov+=covTemp;
//        std::cout<<"cov"<<std::endl;
//        std::cout<<cov<<std::endl;
     } // end for actFreq ...

    
     //     invertWithLapack(cov);

     Matrix<Complex> data;
     data.Resize(actNrParameter+actNrParameterC,actNrParameter+actNrParameterC);
     data=cov;     
     invert(cov);
     //     std::cout<<data*cov<<std::endl;

     J=Complex(0.0,0.0);


     for(UInt parInd=0;parInd<actNrParameter+actNrParameterC;parInd++)
       J+=cov[parInd][parInd];
     J/=(actNrParameter+actNrParameterC);

     Double sum1=0.0;
     Double sum2=0.0;

     for (UInt actFreq=0;actFreq<nrfreq;actFreq++){
       sum1+=std::min(0.0,rhos[actFreq])*std::min(0.0,rhos[actFreq]);
       sum2+= (1-std::max(1.0,rhos[actFreq]))*(1-std::max(1.0,rhos[actFreq]));
       if (sum1!=0.0||sum2!=0){
         std::cout<<" sum1 ... sum2 " <<std::endl;
         std::cout<< sum1 << " ... " << sum2 <<std::endl;
       }
     }

     J+=Complex((1.0/penaltyFactor1)*sum1+(1.0/penaltyFactor2)*sum2,0.0);
    //  penaltyFactor1*=10;
//      penaltyFactor2*=10;


  }// end createJRHo

  void piezoParamIdent::createGradientRho(Vector<Complex> & grad, Double dOmega){
    ENTER_FCN("piezoParamIdent::createGradientRho");

    Complex J1, J2;
    createJRho(J1,FALSE);
    std::cout<<"Value of J1 = "<<J1 <<std::endl;
    Vector<Double> rhosTemp;
    rhosTemp.Resize(nrfreq);
    grad.Resize(nrfreq);
    rhosTemp=rhos;

    for(UInt actFreq=0; actFreq<nrMeasuredData; actFreq ++){
      rhos[actFreq] +=0.001;//*rhos[actFreq];
      createJRho(J2,FALSE);
      grad[actFreq]=(J2-J1)/(rhos[actFreq]-rhosTemp[actFreq]);
      rhos=rhosTemp;
    }
    std::cout<<"gradient von rhoJ ..."<<std::endl;
    std::cout<<grad<<std::endl;


  }


  void piezoParamIdent::descentMethodRho(Complex & functional){
    ENTER_FCN("piezoParamIdent::descentMethodRho");

    UInt maxNumberDescentIterations=20;
    Vector<Complex> grad;
    Vector<Double> rhosOld;
    rhosOld.Resize(nrMeasuredData);
    rhosOld=rhos;
    Double lambda=1.0e+8;

    Complex J_old,J;

   
    Double dOmega=0.0001;
    for (UInt descIter=0;descIter<maxNumberDescentIterations;descIter++){

      createJRho(J_old,FALSE);
      
      J=J_old;
      
      lambda=10*lambda;
      createGradientRho(grad, dOmega);

      while (J_old.real()<=J.real()){

	for(UInt fr=0;fr<nrMeasuredData;fr++)
          rhos[fr]=rhos[fr]-lambda*grad[fr].real();
	createJRho(J,FALSE);
	std::cout<<"NEW Functional J = "<<J.real()<<std::endl;
	std::cout<<"OLD Functional J = "<<J_old.real()<<std::endl;

	if (J.real()>=J_old.real()){
	  lambda=0.1*lambda;
	  rhos=rhosOld;
	  std::cout<<"lambda = " << lambda << std::endl;
	}
	rhosOld=rhos;
      }
      std::cout<<"rhos:"<<std::endl;
      std::cout<<rhos<<std::endl;
      Double rhoInt=0.0;
      *optimalFreqs<<J.real()<<"  ";
      
      for (UInt actFreq=0;actFreq<nrMeasuredData;actFreq++)
        rhoInt+=rhos[actFreq];
      rhoInt/=nrMeasuredData;

      std::cout<<"rhoInt = "<< rhoInt<<std::endl;
      *optimalFreqs<<rhoInt;

      for (UInt actFreq=0;actFreq<nrMeasuredData;actFreq++)
        *optimalFreqs<<rhos[actFreq]<<"  " ;
      *optimalFreqs<<std::endl;


      if (rhoInt>=0.85){
        std::cout<<" Sum rho_i exceeds upper bound ... -> break min J(w) " <<std::endl;
        break;
      }
    }


  }// end descentMethodRho

}
