#include "PDE/SinglePDE.hh" 
#include "piezoParamIdent.hh"

namespace CoupledField
{

  // xxxxxxxxxxxxxxxxxxxxxxxxxx nuMethods  - ITERATION xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  void piezoParamIdent::optimalExpDesignVarNrFreqs(){

    ENTER_FCN("piezoParamIdent::optimalDesignVarNrFreqs");
    std::cout<<"call for optimum experiment desing with variable number of frequencies "<<std::endl;

    Boolean firstTime=TRUE;
    parameterInitial = parameter;
    residuumParIdentOld = residuumParIdent=0.5;


    for(UInt optExpVarFreqCout=0; optExpVarFreqCout<15; optExpVarFreqCout++){

      Double hOmega = (stopfreq-startfreq)/nrfreq;
      freqs.Resize(nrfreq);
      nrMeasuredData=nrfreq;
    
      for (UInt actFreq=0;actFreq<nrfreq;actFreq++)
        freqs[actFreq] = startfreq+actFreq*hOmega;

      Vector<Double> newFreqs;
      readInMeasurement(newFreqs);
      calc_measuredCharge(freqs, real, imag, y_hat); // out of new measurements
    
      if (firstTime==TRUE){
        std::cout<<freqs<<std::endl;
      firstTime=FALSE;
    }

      UInt maxNrBreakpoints = nrfreq;

      Complex J;

      Double dOmega = 0.0001;
      rhos.Resize(nrfreq);

      for(UInt i=0;i<nrfreq;i++)
        rhos[i]=0.1;

      descentMethodRho(J);
      std::cout<<"rhos: "<<std::endl;
      std::cout<<rhos<<std::endl;

      UInt nrNuMethods=0;

      // choose 8 highest values ...
      Vector<Double> highestFreqs;
      Vector<Double> highestRhos;
      UInt maxIndex;
      Double sumRhos=0.0;
      highestFreqs.Resize(nrfreq);
      highestRhos.Resize(nrfreq);

      for (UInt i=0;i<nrfreq;i++){
        maxIndex=0;
        highestFreqs[i]=rhos[0];

        for (UInt j=0;j<rhos.GetSize();j++){
          if (highestRhos[i]<=rhos[j]){
            highestFreqs[i]=freqs[j];
            highestRhos[i]=rhos[j];
            maxIndex=j;
          }
        }
        rhos[maxIndex]=0.0;

        Double omegaDiff=highestFreqs[i]-highestFreqs[0];

        // Get Neighbor of \omega_i
        for(UInt nI=0;nI<highestFreqs.GetSize();nI++){
          if ((highestFreqs[i]-highestFreqs[nI])<omegaDiff)
            omegaDiff=highestFreqs[i]-highestFreqs[nI];
          std::cout<<" omegaDiff = " <<omegaDiff<<std::endl;
        }
        sumRhos+=highestRhos[i]/omegaDiff;

        if (sumRhos>1.15 && i>3){
          std::cout<<"sumRhos = "<<sumRhos <<std::endl;
          //          getchar();
          break;
        }
      }
             
      Integer howManyFreqs=0;
      for (UInt i=0;i<nrfreq;i++)
        if (highestFreqs[i]!=0.0)
          howManyFreqs++;

      std::cout<<"howManyFreqs = "<< howManyFreqs<<std::endl;

      // sort highestFreqs
      Double dv,dr;
      Integer k, ii;

      for (k=0; k<howManyFreqs-1; ++k)
        for (ii=k+1; ii<howManyFreqs; ++ii)
          if (highestFreqs[ii] < highestFreqs[k])
            {
              dv = highestFreqs[k];
              dr = highestRhos[k];
              highestFreqs[k] = highestFreqs[ii];
              highestRhos[k] = highestRhos[ii];
              highestFreqs[ii] = dv;
              highestRhos[ii] = dr;
            }

 //      std::cout<<"highestRhos"<<std::endl;
//       std::cout<<highestRhos<<std::endl;
      std::cout<<"highestFreqs"<<std::endl;
      std::cout<<highestFreqs<<std::endl;

      for(UInt i=0; i<highestFreqs.GetSize();i++)
        *optimalFreqs<<highestFreqs[i]<<"  ";
      *optimalFreqs<<std::endl;
      

      nrMeasuredData=howManyFreqs;
      rhos.Resize(howManyFreqs);
      freqs.Resize(howManyFreqs);
      real.Resize(howManyFreqs);
      imag.Resize(howManyFreqs);

      rhos=highestRhos;
      freqs=highestFreqs;

      readInMeasurement(newFreqs);

      calc_measuredCharge(freqs, real, imag, y_hat); // out of new measurements
      newtonCounter=0;
    
      while (nrNuMethods<maxNumberNewtonLoops){
        //    nuMethodsC2();
        nuMethods();
        nrNuMethods++;
        for (UInt i=0; i<parameter.GetSize(); i++)
          *parFinal<<parameter[i]<<", ";
        *parFinal<<"/"<<std::endl;

        for (UInt i=0; i<parameterC.GetSize(); i++)
          *parFinal<<parameterC[i]<<", ";
        *parFinal<<"/"<<std::endl;

        newtonCounter++;
      }

      //  nrfreq++;

    }

  }

  void piezoParamIdent::createJRho(Complex &J, Boolean writeOutCov){
    ENTER_FCN("piezoParamIdent::createJRho");

    Double  hOmega = (stopfreq-startfreq)/nrfreq;

    //    Complex J;
    J=Complex(0.0,0.0);


    Vector<Complex> jacobi;
    Vector<Complex> jacobiH;

    Double penaltyFactor1=0.000000001;
    Double penaltyFactor2=0.000000001;

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

       if (freqs[actFreq]<=fr_||freqs[actFreq]>=fa_){   

         for(UInt i=0;i<actNrParameter+actNrParameterC;i++)
           for(UInt j=0;j<actNrParameter+actNrParameterC;j++){
             //           covTemp[i][j]=jacobiH[i]*jacobi[j]/(delta*delta*std::abs(y_hat[actFreq])*std::abs(y_hat[actFreq]));
             if(actFreq==0||actFreq==nrfreq)
               //  covTemp[i][j]=rhos[actFreq]*hOmega*0.5*jacobiH[i]*jacobi[j]/
               // ((1.0+delta)*(1.0+delta)*std::abs(y_hat[actFreq])*std::abs(y_hat[actFreq]));

               covTemp[i][j]=rhos[actFreq]*hOmega*0.5*jacobiH[i]*jacobi[j]/
                 (0.8*freqs[actFreq]*(1.0+delta)*(1.0+delta)*std::abs(y_hat[actFreq])*std::abs(y_hat[actFreq]));

//                 covTemp[i][j]=rhos[actFreq]*hOmega*0.5*jacobiH[i]*jacobi[j]/
//                   (0.8*freqs[actFreq]*(1+delta)*(1+delta)*std::abs(y_hat[actFreq])*std::abs(y_hat[actFreq]));

             else 
               covTemp[i][j]=rhos[actFreq]*hOmega*jacobiH[i]*jacobi[j]/
                 (0.8*freqs[actFreq]*(1.0+delta)*(1.0+delta)*std::abs(y_hat[actFreq])*std::abs(y_hat[actFreq]));

//                covTemp[i][j]=rhos[actFreq]*hOmega*jacobiH[i]*jacobi[j]/
//                  (freqs[actFreq]*(1+delta)*(1+delta)*std::abs(y_hat[actFreq])*std::abs(y_hat[actFreq]));

             if (i==j)
               covTemp[i][i]+=covTemp[i][i];// *1.0e-10;
           }
       }
       cov+=covTemp;
//        std::cout<<"cov"<<std::endl;
//        std::cout<<cov<<std::endl;
     } // end for actFreq ...

    
     //     invertWithLapack(cov);

     Matrix<Complex> data;
     data.Resize(actNrParameter+actNrParameterC,actNrParameter+actNrParameterC);
     data=cov;     
     invert(cov);
     Vector<Double> covDiag(actNrParameter+actNrParameterC);
     for (UInt ii=0;ii<actNrParameter+actNrParameterC;ii++)
       covDiag[ii]=cov[ii][ii].real();
//      std::cout<<"covDiag"<<std::endl;
//      std::cout<<covDiag<<std::endl;
     //     std::cout<<data*cov<<std::endl;

     J=Complex(0.0,0.0);

     if (writeOutCov==TRUE){

       if (actNrParameter==3){      
         *confInterval<<parameter[1]+parameterInitial[1]*std::sqrt(cov[0][0].real()*0.115*0.115)<<"  ";
         //        std::cout<<parameter[1]+parameter[1]*std::sqrt(cov[1][1].real()*0.115*0.115)<<std::endl;
         *confInterval<<parameter[1]<<"  ";
         *confInterval<<parameter[1]-parameterInitial[1]*std::sqrt(cov[0][0].real()*0.115*0.115)<<"  ";
         //        std::cout<<parameter[1]-1000.0*std::sqrt(cov[0][0].real()*parameter[1]*0.115*0.115)<<std::endl;
         
         
         *confInterval<<parameter[7]+parameterInitial[7]*std::sqrt(cov[1][1].real()*0.115*0.115)<<"  ";
         *confInterval<<parameter[7]<<"  ";
         *confInterval<<parameter[7]-parameterInitial[7]*std::sqrt(cov[1][1].real()*0.115*0.115)<<"  ";
         
         
         *confInterval<<parameter[9]+parameterInitial[9]*std::sqrt(cov[2][2].real()*0.115*0.115)<<"  ";
         *confInterval<<parameter[9]<<"  ";
         *confInterval<<parameter[9]-parameterInitial[9]*std::sqrt(cov[2][2].real()*0.115*0.115)<<"  ";
       }

       if (actNrParameter==10)
         for (UInt i=0;i<10;i++){
           *confInterval<<parameter[i]+parameterInitial[i]*std::sqrt(cov[i][i].real()*0.155*0.155)<<"  ";
           *confInterval<<parameter[i]<<"  ";
           *confInterval<<parameter[i]-parameterInitial[i]*std::sqrt(cov[i][i].real()*0.155*0.155)<<"  ";
         }


       *confInterval<<std::endl;
     }

      




     for(UInt parInd=0;parInd<actNrParameter+actNrParameterC;parInd++)
       J+=cov[parInd][parInd];
     J/=(actNrParameter+actNrParameterC);

     Double sum1=0.0;
     Double sum2=0.0;

     for (UInt actFreq=0;actFreq<nrfreq;actFreq++){
       sum1+=std::min(0.0,rhos[actFreq])*std::min(0.0,rhos[actFreq]);
       sum2+= (1-std::max(1.0,rhos[actFreq]))*(1-std::max(1.0,rhos[actFreq]));
//        }
     }

     if (sum1!=0.0||sum2!=0){
       std::cout<<" sum1 ... sum2 " <<std::endl;
       std::cout<<1.0/penaltyFactor1 * sum1 << " ... " << 1.0/penaltyFactor1*sum2 <<std::endl;
     }

     J+=Complex((1.0/penaltyFactor1)*sum1+(1.0/penaltyFactor2)*sum2,0.0);
    //  penaltyFactor1*=10;
//      penaltyFactor2*=10;


  }// end createJRHo

  void piezoParamIdent::createGradientRho(Vector<Complex> & grad, Double dOmega){
    ENTER_FCN("piezoParamIdent::createGradientRho");

    Complex J1, J2;
    createJRho(J1,TRUE);
    std::cout<<"Value of J1 = "<<J1 <<std::endl;
    Vector<Double> rhosTemp;
    rhosTemp.Resize(nrfreq);
    grad.Resize(nrfreq);
    rhosTemp=rhos;

    for(UInt actFreq=0; actFreq<nrMeasuredData; actFreq ++){
      rhos[actFreq] = 1.0001*rhos[actFreq];
      createJRho(J2,FALSE);
      grad[actFreq]=(J2-J1)/(rhos[actFreq]-rhosTemp[actFreq]);
      rhos=rhosTemp;
    }
//     std::cout<<"gradient von rhoJ ..."<<std::endl;
//     std::cout<<grad<<std::endl;

  }


  void piezoParamIdent::descentMethodRho(Complex & functional){
    ENTER_FCN("piezoParamIdent::descentMethodRho");

    UInt maxNumberDescentIterations=1;
    Vector<Complex> grad;
    Vector<Double> rhosOld;
    rhosOld.Resize(nrMeasuredData);
    rhosOld=rhos;
    Double lambda=1.0e-4; // for 10 parameters
    //    Double lambda=1.0e-1; // for 3 parameters
    //    Double lambda=1.0e-2; // for 4 parameters

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
//       std::cout<<"rhos:"<<std::endl;
//       std::cout<<rhos<<std::endl;
      Double rhoInt=0.0;
      *rhosOut<<J.real()<<"  ";
      
      for (UInt actFreq=0;actFreq<nrMeasuredData;actFreq++)
        rhoInt+=rhos[actFreq];
      rhoInt/=nrMeasuredData;

      std::cout<<"rhoInt = "<< rhoInt<<std::endl;
      *rhosOut<<rhoInt;

      for (UInt actFreq=0;actFreq<nrMeasuredData;actFreq++)
        *rhosOut<<rhos[actFreq]<<"  " ;
      *rhosOut<<std::endl;

      if (rhoInt>=0.85){
        std::cout<<" Sum rho_i exceeds upper bound ... -> break min J(w) " <<std::endl;
        break;
      }
    }


  }// end descentMethodRho

}
