#include "PDE/SinglePDE.hh" 
#include "PDE/piezoPDE.hh"
#include "piezoParamIdent.hh"



namespace CoupledField
{

  // xxxxxxxxxxxxxxxxxxxxxxxxxx nuMethods  - ITERATION xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  void piezoParamIdent::optimalExpDesign(){

    ENTER_FCN("piezoParamIdent::optimalDesign");
    std::cout<<"call for optimum experiment desing .."<<std::endl;

#ifndef USE_LAPACK
    std::cout<<"\n! Optimum experiment design works with LAPACK Routines"<<std::endl;
    std::cout<<"! Please set LAPACK = yes in your Makefile.option (CFS & OLAS)" <<std::endl;
    std::cout<<"! Please, leave cfs with Ctrl + c ... " <<std::endl;
    getchar();
#endif
#ifdef USE_LAPACK
    std::cout<<"! LAPACK routines will be called throughout the program run" <<std::endl;
#endif



    // optimalExpDesignDiffNumberFreqs();
      nrMeasuredData=10;
      Vector<Double> freqs5;
      freqs5.Resize(10);
      freqs5[0]=3.0e+06;
      freqs5[1]=3.5e+06;
      freqs5[2]=3.75e+06;
      freqs5[3]=4.0e+06;
      freqs5[4]=4.1e+06;
      freqs5[5]=4.2e+06;
      freqs5[6]=5.6e+06;
      freqs5[7]=6.0e+06;
      freqs5[8]=6.5e+06;
      freqs5[9]=7.0e+06;
      //      freqs5[10]=7.75e+06;
      freqs=freqs5;

//      nrMeasuredData=5;
//      Vector<Double> freqs5;
//      freqs5.Resize(5);

//      freqs5[0]=2.0e+06;
//      freqs5[1]=3.5e+06;
//      freqs5[2]=4.5e+06;
//      freqs5[3]=5.8e+06;
//      freqs5[4]=6.9e+06;
//      freqs=freqs5;

    Vector<Double> newFreqs;
   
    for(UInt fr=0;fr<nrMeasuredData;fr++)
      *impedCurve<< freqs[fr]<<"  ";
    *impedCurve<<std::endl;

    readInMeasurement(newFreqs);

    calc_measuredCharge(freqs, real, imag, y_hat); // out of new measurements

    //  actNrParameterC=0;
    MaterialData * ptMaterial=ptMyPDE_->getPDEMaterialData();   // Pointer to MaterialData

    // Frequency bounds

    // fr = 4000m/s / 2*thickness

    fr_=4000/(2*thickness);
    // minus 5 percent
    fr_=0.94*fr_;
    // fa = 1/thickness * sqrt(c_33/rho)
   
    Double rho;
    rho = ptMaterial->GetDensity();

    fa_ = 1.0/(2*thickness) * std::sqrt(1.3e+11/rho);
    fa_=1.2*fa_;

    std::cout<<"++ Bounds for resonance and antiresonace frequency"<<std::endl;
    std::cout<<"fr:"<<std::endl;
    std::cout<<fr_<<std::endl;

    std::cout<<"fa:"<<std::endl;
    std::cout<<fa_<<std::endl;
    Vector<Complex> jacobi;


    Complex functional;

    for (UInt nrOptExpSteps=0; nrOptExpSteps<20; nrOptExpSteps++){ 
      UInt nrNuMethods=0;

      //      descentMethod(functional);
      readInMeasurement(newFreqs);

      for(UInt fr=0;fr<nrMeasuredData;fr++)
        *impedCurve<< freqs[fr]<<"  ";
      *impedCurve<<std::endl;

      calc_measuredCharge(freqs, real, imag, y_hat); // out of new measurements
      newtonCounter=0;

      while (nrNuMethods<maxNumberNewtonLoops && residuum< 1.0e-2*tau*(1+delta)){
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
      
    } // end nrOptExpSteps
    
    std::cout<<"\n END OF OPTIMAL EXPERIMENT DESIGN" <<std::endl;
     
  } // end optimalExpDesign

  void piezoParamIdent::optimalExpDesignDiffNumberFreqs(){

    ENTER_FCN("piezoParamIdent::optimalExpDesignDiffNumberFreqs");

    nrMeasuredData=1;
    Vector<Double> newFreqs;

    Vector<Double> freqs2;
    freqs2.Resize(2);
    freqs2[0]=3.5e+06;
    freqs2[1]=6.0e+06;
    Vector<Double> freqs3;
    freqs3.Resize(3);
    freqs3[0]=3.5e+06;
    freqs3[1]=4.0e+06;
    freqs3[2]=6.0e+06;
    Vector<Double> freqs4;
    freqs4.Resize(4);
    freqs4[0]=3.5e+06;
    freqs4[1]=4.0e+06;
    freqs4[2]=5.9e+06;
    freqs4[3]=6.5e+06;
    Vector<Double> freqs5;
    freqs5.Resize(5);
    freqs5[0]= 3.0e+06;
    freqs5[1]= 4.0e+06;
    freqs5[2]= 5.9e+06;
    freqs5[3]= 6.5e+06;
    freqs5[4]= 7.0e+06;
    Vector<Double> freqs6;
    freqs6.Resize(6);
    freqs6[0]= 3.5e+06;
    freqs6[1]= 3.9e+06;
    freqs6[2]= 4.3e+06;
    freqs6[3]= 5.8e+06;
    freqs6[4]= 6.2e+06;
    freqs6[5]= 6.5e+06;
    Vector<Double> freqs7;
    freqs7.Resize(7);
    freqs7[0]= 3.3e+06;
    freqs7[1]= 3.8e+06;
    freqs7[2]= 4.0e+06;
    freqs7[3]= 4.5e+06;
    freqs7[4]= 5.8e+06;
    freqs7[5]= 6.4e+06;
    freqs7[6]= 7.0e+06;

    newtonCounter=0;

    for (UInt nrOptExpSteps=0; nrOptExpSteps<10; nrOptExpSteps++){ 

      nrMeasuredData++;
      
      freqs.Resize(nrMeasuredData);

      if (nrMeasuredData==2)
        freqs=freqs2;
      else  if (nrMeasuredData==3)
        freqs=freqs3;
      if (nrMeasuredData==4)
        freqs=freqs4;
      else  if (nrMeasuredData==5)
        freqs=freqs5;
      if (nrMeasuredData==6)
        freqs=freqs6;
      else  if (nrMeasuredData==7)
        freqs=freqs7;
      else if(nrMeasuredData>7)
        break;

      readInMeasurement(newFreqs);
      
      calc_measuredCharge(freqs, real, imag, y_hat); // out of new measurements 

      UInt nrNuMethods=0;
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

    }


  } // end optimalExpDesignDiffNumberFreqs



  void piezoParamIdent::descentMethod(Complex & functional){
    ENTER_FCN("piezoParamIdent::descentMethod");

    Vector<Complex> gradient;
    //    Double lambda = 0.00002; // for three parameter
    //Double lambda = 1.0e+8; // for three parameter
    // Double lambda = 1.0e-19; // for nine parameter
    //Double lambda = 1.0e+3; // for nine parameter
    Double lambda = 1.0e+5; // for nine parameter
    //Double lambda = 1.0e-7; // for nine parameter
    //    Double lambda = 1.0; // for nine parameter using M - confidence crit.
    //    Double lambda = 1.0e+7; // for nine parameter using M - confidence crit.
    //Double lambda = 40; // for three parameters using M - confidence crit.

    UInt maxNumberDescentIterations=2;
   
    for (UInt descIter=0;descIter<maxNumberDescentIterations;descIter++){

      Complex J, JOld;
      Vector<Double> freqsOld;
      freqsOld.Resize(nrMeasuredData);
      freqsOld=freqs;
      createCovA(J, TRUE);
      std::cout<<" Value of functional J(w) = "<< J <<std::endl;
      JOld=J;
      //      Boolean scaleFlag=FALSE;

      std::cout<<J<<std::endl;

      //      createGradient(gradient,10.0);
      createGradient(gradient,1.0e-5);
      lambda = 2*lambda;

      std::cout<<"gradient"<<std::endl;
      std::cout<<gradient<<std::endl;
      //      getchar();
      
      *optimalFreqs<< descIter <<"  "<< J.real() <<"  ";
      //     for(UInt fr=0;fr<nrMeasuredData;fr++){
      // *impedCurve<< freqs[fr]<<"  ";
         //         std::cout<<"freqs written out, "<< fr<<") = "<<freqs[fr]<<std::endl;
      //       }
      *optimalFreqs<<std::endl;    

      UInt innerCounter=0;
      //      while(innerCounter<1){
        while(JOld.real()<=J.real()){
        innerCounter ++;
      
        for(UInt fr=0;fr<nrMeasuredData;fr++){
          freqs[fr]=freqs[fr]-lambda*gradient[fr].real();
          
          std::cout<<"freqs:( " << fr <<" ) = "<< freqs[fr]<<std::endl;

      //     if (freqs[fr]<=1.11e+06||freqs[fr]>=8.0e+06)
//             {
//               std::cout<<"! There are no measurements for this frequency ...  " << std::endl;
//               scaleFlag = TRUE;
//             }
        } // end for fr ...

      //   if (lambda<1.0e-9)
//           break;
        if (innerCounter>40)
          break;

        createCovA(J,FALSE);

        //        if (J.real()>=JOld.real()||scaleFlag){
        if (J.real()>=JOld.real()){
//           if(J.real()/JOld.real()<1.1){
//             std::cout<<"J/JOld = " <<J.real()/JOld.real()<<std::endl;
//             break;
//           }

          lambda=0.5*lambda;
          freqs=freqsOld;
          std::cout<<"! Reduced lambda in descent method ... lambda =  " << lambda<<std::endl;
          if (J.real()>=JOld.real())
            std::cout<<"  - because JOld<=J " << JOld.real()<<" <= " <<J.real() << std::endl;
//           if (scaleFlag)
//             std::cout<<"  - because scaleFlag==TRUE " << std::endl;

 //          scaleFlag=FALSE;
        }
        }
        for(UInt fr=0;fr<nrMeasuredData;fr++){
          if ((freqs[fr]<=fa_)&&(freqs[fr]>=fr_)){
            //            std::cout<<"! Frequency in dangerous intervall [" <<fr_<< ", "<< fa_<<"] "<<std::endl;
            if(freqs[fr]<(fr_+fa_)/2){
              freqs[fr]=fr_;
              std::cout<<"! Frequency "<< fr <<" to close to fr, fixed it at f= " << fr_<< std::endl;
            }
            else{
              freqs[fr]=fa_;
              std::cout<<"! Frequency "<< fr <<" to close to fa, fixed it at f= " << fa_<< std::endl;
            }
          } 
          if (freqs[fr]<=2.0e+06){
            std::cout<<"! There are no measurements for frequency "<< fr  <<std::endl;
            std::cout<<"! Fixed frequency at 2.0e+06 " <<std::endl;
            freqs[fr]=2.0e+06;
          }
          if(freqs[fr]>=8.0e+06){
            std::cout<<"! There are no measurements for frequency"<<fr <<std::endl;
            std::cout<<"! Fixed frequency at 9.0e+06 " <<std::endl;
            freqs[fr]=8.0e+06;
          }
            
            
        }// end in freqs[fr]<=fa_ ...
          
          
        // end while ....
      
        for(UInt fr=0;fr<nrMeasuredData;fr++)
          *optimalFreqs<< freqs[fr]<<"  ";
        *optimalFreqs<<std::endl;    
          
    }      

  }// descentMethod


  void piezoParamIdent::createCovA(Complex &J, Boolean writeOutCov){
    ENTER_FCN("piezoParamIdent::createCovA");

    //    Complex J;
    J=Complex(0.0,0.0);

    Vector<Complex> jacobi;
    Vector<Complex> jacobiH;
    //    MaterialData * ptMaterial=ptMyPDE_->getPDEMaterialData();   // Pointer to MaterialData
    Complex F_hat_incr1,F_hat_incr2;

    Matrix<Complex> cov(actNrParameter+actNrParameterC, actNrParameter+actNrParameterC);
    Matrix<Complex> covTemp(actNrParameter+actNrParameterC, actNrParameter+actNrParameterC);

     for (UInt actFreq=0;actFreq<nrMeasuredData;actFreq++){

       createJacobian(jacobi, freqs[actFreq]);
       jacobiH.Resize(actNrParameter+actNrParameterC);

       for (UInt i=0;i<jacobi.GetSize();i++)
         jacobiH[i]=Complex(jacobi[i].real(),-jacobi[i].imag());
       
       for(UInt i=0;i<actNrParameter+actNrParameterC;i++)
         for(UInt j=0;j<actNrParameter+actNrParameterC;j++){
           //           covTemp[i][j]=jacobiH[i]*jacobi[j]/(delta*delta*std::abs(y_hat[actFreq])*std::abs(y_hat[actFreq]));
           covTemp[i][j]=jacobiH[i]*jacobi[j]/((1.0+delta)*(1.0+delta)*std::abs(y_hat[actFreq])*std::abs(y_hat[actFreq]));
             
           if (i==j)
             covTemp[i][i]=covTemp[i][i]+covTemp[i][i];// *1.0e-10;
         }
       cov+=covTemp;
     } // end for actFreq ...

     //     invertWithLapack(cov);

//      for(UInt i=0;i<cov.GetSizeRow();i++)
//        for (UInt j=0; j<cov.GetSizeCol();j++){
//          std::cout<<cov[i][j].real()<<"+"<<cov[i][j].imag()<<"i"<< ", ";
//           //std::cout<<"F'("<<i<<")("<<j<<")= "<< JacobiMatrix[i][j]<<"; \t";
//           if (j==cov.GetSizeCol()-1)
//              std::cout<<";\n";
//        }

     Matrix<Complex> data;
     data.Resize(actNrParameter+actNrParameterC,actNrParameter+actNrParameterC);
     data=cov;     
     invert(cov);

     //     std::cout<<"inv(cov):"<<std::endl;
//        std::cout<<cov<<std::endl;
    //    std::cout<<"inv(cov)*cov :"<<std::endl;

     data=data*cov;
     Complex diagC=Complex(0.0,0.0);
     for(UInt i=0;i<actNrParameter+actNrParameterC;i++)
       diagC+=data[i][i];
     diagC/=Complex(actNrParameter+actNrParameterC,0.0);
     //     std::cout<<"Trace of inv(cov)*cov = "<< diagC <<std::endl;
     if (std::abs(diagC)>1.5)
       std::cout<<" ! Inversion of Cov failed!!"<<std::endl;

  //    cov[0][0]/=parameter[1];
//      cov[1][1]/=parameter[7];
//      cov[2][2]/=parameter[9];


   //   std::cout<<"cov"<<std::endl;
//      std::cout<<cov<<std::endl;

// //       getchar();
//       cov=data;


     // confidence - intervall criterion
     if (FALSE){
       Double JJ=0.0;
       for (UInt i=0;i<actNrParameter+actNrParameterC;i++){
         if (std::sqrt(cov[i][i].real()) >= JJ)
           JJ=std::sqrt(cov[i][i].real());
         if (writeOutCov==TRUE)
           *optimalFreqs<<cov[i][i].real()<<"  ";
       }
       J=Complex(JJ,0.0);
     }

     // A - Criterion
    if(TRUE){
      Vector<Double> covDiag(actNrParameter+actNrParameterC);
      for(UInt i=0; i<actNrParameter+actNrParameterC; i++){
        J+=cov[i][i];
        covDiag[i]=cov[i][i].real();
      }
      J=J/Complex(actNrParameter+actNrParameterC,0.0);
      
      if (writeOutCov==TRUE){
        std::cout<<"covDiag"<<std::endl;
        std::cout<<covDiag<<std::endl;
      
        for(UInt i=0; i<actNrParameter+actNrParameterC; i++){
          *confInterval<<cov[i][i].real()<<"  ";
        }

        if (actNrParameter==3){      
          *confInterval<<parameter[1]+parameter[1]*std::sqrt(cov[0][0].real()*0.115*0.115)<<"  ";
        //        std::cout<<parameter[1]+parameter[1]*std::sqrt(cov[1][1].real()*0.115*0.115)<<std::endl;
          *confInterval<<parameter[1]<<"  ";
          *confInterval<<parameter[1]-parameter[1]*std::sqrt(cov[0][0].real()*0.115*0.115)<<"  ";
        //        std::cout<<parameter[1]-1000.0*std::sqrt(cov[0][0].real()*parameter[1]*0.115*0.115)<<std::endl;


          *confInterval<<parameter[7]+parameter[7]*std::sqrt(cov[1][1].real()*0.115*0.115)<<"  ";
          *confInterval<<parameter[7]<<"  ";
          *confInterval<<parameter[7]-parameter[7]*std::sqrt(cov[1][1].real()*0.115*0.115)<<"  ";


          *confInterval<<parameter[9]+parameter[9]*std::sqrt(cov[2][2].real()*0.115*0.115)<<"  ";
          *confInterval<<parameter[9]<<"  ";
          *confInterval<<parameter[9]-parameter[9]*std::sqrt(cov[2][2].real()*0.115*0.115)<<"  ";
        }

        else if (actNrParameter==4){
          for (UInt ii=0;ii<=4;ii++){
            if(ii<=1){
              *confInterval<<parameter[ii]+parameter[ii]*std::sqrt(cov[ii][ii].real()*0.115*0.115)<<"  ";
              *confInterval<<parameter[ii]<<"  ";
              *confInterval<<parameter[ii]-parameter[ii]*std::sqrt(cov[ii][ii].real()*0.115*0.115)<<"  ";
            }
            else if(ii>2)
              {
                *confInterval<<parameter[ii]+parameter[ii]*std::sqrt(cov[ii-1][ii-1].real()*0.115*0.115)<<"  ";
                *confInterval<<parameter[ii]<<"  ";
                *confInterval<<parameter[ii]-parameter[ii]*std::sqrt(cov[ii-1][ii-1].real()*0.115*0.115)<<"  ";
              }
          }
        }
        if (actNrParameter==10){
          for (UInt ii=0;ii<10;ii++){
            *confInterval<<parameter[ii]+parameter[ii]*std::sqrt(cov[ii][ii].real()*2.5*2.5)<<"  ";
            *confInterval<<parameter[ii]<<"  ";
            *confInterval<<parameter[ii]-parameter[ii]*std::sqrt(cov[ii][ii].real()*2.5*2.5)<<"  ";
          }
        }
        if (actNrParameter==5){
          for (UInt ii=0;ii<5;ii++){
            *confInterval<<parameter[ii]+parameter[ii]*std::sqrt(cov[ii][ii].real()*0.55*0.55)<<"  ";
            *confInterval<<parameter[ii]<<"  ";
            *confInterval<<parameter[ii]-parameter[ii]*std::sqrt(cov[ii][ii].real()*0.55*0.55)<<"  ";
          }
        }   

      }
    }

    
    // smallest eigenvalue criterion:
    if(FALSE){

      Vector<Double> eig;
      eig.Resize(actNrParameter+actNrParameterC);
#ifdef USE_LAPACK
      cov.eigenvaluesWithLapack(eig);
#endif 

      std::cout<<"eigenvalues"<<std::endl;
      std::cout<<eig<<std::endl;
      
      Double eigMin=std::abs(eig[0]);
      for(UInt eigInd=0; eigInd<actNrParameter+actNrParameterC; eigInd++){
        if(std::abs(eig[eigInd])>eigMin)
          eigMin=std::abs(eig[eigInd]);
      }
      J=std::log(eigMin);


// #ifndef USE_LAPACK
//     std::cout<<"Optimum experiment design works with LAPACK Routines"<<std::endl;
//     std::cout<<"Please set LAPACK = yes in your Makefile.option (CFS & OLAS)" <<std::endl;
// #endif

// #ifdef USE_LAPACK
//       char lp_jobz='V';
//       char lp_uplo='L';
//       Integer lp_N=actNrParameter+actNrParameterC;                
//       Integer lp_lda=actNrParameter+actNrParameterC;
      
//       // array contains ev in ascending order
//       Vector<Double> lp_w;
//       lp_w.Resize(actNrParameter+actNrParameterC);
      
//       Integer lp_lworkf77=99;
      
//       // workspace array - complex 16 array
//       Vector<Complex> lp_work;
//       lp_work.Resize(lp_lworkf77);
      
//       // workspace array - double precission
//       Vector<Double> lp_rwork;
//       lp_rwork.Resize(3*actNrParameter+actNrParameterC-2);
      
//       Integer lp_infof77;
//       F77complex16 auxValC;
//       F77real8 auxValR;
      
      
//       // Convert CFS++ Vector<Complex> to Vector<F77complex16>
//       for ( UInt count = 0; count < actNrParameter+actNrParameterC; count++ ) 
//         for ( UInt countC = 0; countC < actNrParameter+actNrParameterC; countC++ ) {
//           CC2F77( cov[count][countC], auxValC );
//           lp_af77[countC*actNrParameter+actNrParameterC+count] = auxValC;
//         }

//       for ( UInt count = 0; count < actNrParameter+actNrParameterC; count++ ) {
//         CC2F77( lp_w[count], auxValR );
//         lp_wf77[count] = auxValR;
//       }
      
//       for (UInt count=0; count < lp_rwork.GetSize();count++){
//         CC2F77(lp_rwork[count], auxValR);
//         lp_rworkf77[count] = auxValR;
//       }
      
//       //     void LP_ZHEEV( char*, char*, int*, F77complex16*, int *, F77real8*, F77complex16*, int*, F77real8* ,int* ); 
//       zheev_( &lp_jobz, &lp_uplo, &lp_N, lp_af77, &lp_lda, lp_wf77, lp_workf77, &lp_lworkf77, lp_rworkf77 ,&lp_infof77); 
      
//       // reconvert f772C++
//       for (UInt count=0; count < lp_work.GetSize();count++)
//         F772CC(lp_workf77[count], lp_work[count]);
      
//       for ( UInt count = 0; count < actNrParameter+actNrParameterC; count++ ) 
//         F772CC( lp_wf77[count], lp_w[count] );
      
// //        std::cout<<"eigenvalues"<<std::endl;
// //        std::cout<<lp_w<<std::endl;
      
//       Double eigMin=std::abs(lp_w[0]);
//       for(UInt eigInd=0; eigInd<actNrParameter+actNrParameterC; eigInd++){
//         if(std::abs(lp_w[eigInd])>eigMin)
//           eigMin=std::abs(lp_w[eigInd]);
//       }
//       J=std::log(eigMin);
// #endif

    } // end minim EV criterion

    //      std::cout<<" A-Kriterium variance - covariance Matrix J= " << J <<std::endl;  

  }// end createCovA


  void piezoParamIdent::createGradient(Vector<Complex> & grad, Double dOmega){

    ENTER_FCN("piezoParamIdent::createGradient");

    Vector<Complex> jacobi1, jacobi2;
    Vector<Double> freqsOld;
    freqsOld.Resize(nrMeasuredData);
    freqsOld=freqs;

    grad.Resize(nrMeasuredData);
    Complex J1,J2;
    createCovA(J1,FALSE);

    for (UInt j=0;j<nrMeasuredData;j++){
      freqs[j]=freqs[j]+dOmega*freqs[j];
      createCovA(J2,FALSE);
      for(UInt i=0; i<actNrParameter+actNrParameterC;i++)
      grad[j]= (J2-J1)* voltage*voltage / (dOmega*freqs[j]) ;
      freqs[j]=freqsOld[j];
    }
    //    grad=grad*Complex(1.0e+05,0);

//      std::cout<<"Gradient of J" <<std::endl;
//      std::cout<<grad <<std::endl;
//      getchar();


  } // end create Gradient

  void piezoParamIdent::invertWithLapack(Matrix<Complex> & data){
    ENTER_FCN("piezoParamIdent::invertWithLapack");

#ifndef USE_LAPACK
    std::cout<<"Optimum experiment design works with LAPACK Routines"<<std::endl;
    std::cout<<"Please set LAPACK = yes in your Makefile.option (CFS & OLAS)" <<std::endl;
#endif

#ifdef USE_LAPACK

    Matrix<Complex> cov;
    cov.Resize(actNrParameter+actNrParameterC,actNrParameter+actNrParameterC);
    cov=data;

    // data.Resize(actNrParameter+actNrParameterC, actNrParameter+actNrParameterC);
//     data[0][0]=Complex(1.0,0.0);
//     data[0][1]=Complex(1.0,1.0);
//     data[0][2]=Complex(2.0,0.0);
//     data[1][0]=Complex(1.0,1.0);
//     data[1][1]=Complex(1.0,1.0);
//     data[1][2]=Complex(0.0,0.0);
//     data[2][0]=Complex(1.0,0.0);
//     data[2][1]=Complex(0.0,0.0);
//     data[2][2]=Complex(1.0,1.0); 



   //  Vector<Complex> lp_covVec;
//     lp_covVec.Resize(actNrParameter+actNrParameterC*actNrParameter+actNrParameterC);
//     //  lp_covVec.Resize(6);

//     // Complex divisor = cov[0][0];
//     // std::cout<<"\n wo geht es raus ? 0 " <<std::endl;
//     for (UInt i=0;i<actNrParameter+actNrParameterC;i++)
//       for (UInt j=0;j<actNrParameter+actNrParameterC;j++){
//         lp_covVec[i+j*actNrParameter+actNrParameterC]=cov[i][j];// /divisor;
//        }
 
//     char  lp_matType;
//     lp_matType='L';

//     Integer lp_nrRHS, lp_loadDim, lp_loadDimRHS, lp_info, lp_lwork,lp_dim;
//     lp_nrRHS = actNrParameter+actNrParameterC;

//     lp_lwork=192;
//     lp_loadDim =actNrParameter+actNrParameterC;
//     lp_loadDimRHS = actNrParameter+actNrParameterC;
//     lp_dim=actNrParameter+actNrParameterC;

//     Integer *lp_interchanges;
//     Vector<Complex> lp_rhsVec, lp_rhsVec1, lp_rhsVec2, lp_work;

//     lp_rhsVec.Resize(actNrParameter+actNrParameterC*actNrParameter+actNrParameterC);
//     lp_work.Resize(192);

//     F77complex16 *lp_rhsVecf77, *lp_covVecf77, *lp_workf77;

    Matrix<Complex> rhsMat;
    rhsMat.Resize(actNrParameter+actNrParameterC, actNrParameter+actNrParameterC);
    for(UInt i=0;i<actNrParameter+actNrParameterC;i++)
      rhsMat[i][i]=Complex(1.0,0.0);

// //     std::cout<<"rhsMat:"<<std::endl;
// //     std::cout<<rhsMat<<std::endl;

//     Matrix<Complex> x;

#ifdef USE_LAPACK
    lapackSysMatType LAPACK_SYS_MAT_TYPE = ZGESV;
    data.solveWithLapack(rhsMat,LAPACK_SYS_MAT_TYPE);
#endif 

    data=rhsMat;


    // for (UInt i=0;i<actNrParameter+actNrParameterC;i++)
//       for (UInt j=0;j<actNrParameter+actNrParameterC;j++){
//         lp_rhsVec[i+j*actNrParameter+actNrParameterC]=rhsMat[i][j];
//       }
// //     std::cout<<"lp_rhsVec:"<<std::endl;
// //     std::cout<<lp_rhsVec<<std::endl;


//     NewArray( lp_rhsVecf77, F77complex16, actNrParameter+actNrParameterC*actNrParameter+actNrParameterC );
//     NewArray( lp_interchanges, int, actNrParameter+actNrParameterC );
//     NewArray( lp_covVecf77, F77complex16, actNrParameter+actNrParameterC*actNrParameter+actNrParameterC );
//     NewArray( lp_workf77, F77complex16, 192 );
//     F77complex16 auxVal2;

//     // Convert CFS++ Vector<Complex> to Vector<F77complex16>
//     for ( UInt count = 0; count < actNrParameter+actNrParameterC*actNrParameter+actNrParameterC; count++ ) {
//       CC2F77( lp_rhsVec[count], auxVal2 );
//       lp_rhsVecf77[count] = auxVal2;
//     }

//     for (UInt count = 0; count < actNrParameter+actNrParameterC*actNrParameter+actNrParameterC; count++ ) {
//       CC2F77( lp_covVec[count], auxVal2 );
//       lp_covVecf77[count] = auxVal2;
//     }

//     for (UInt count=0; count < lp_work.GetSize();count++){
//       CC2F77(lp_work[count], auxVal2);
//       lp_workf77[count] = auxVal2;
//     }

//     Integer lp_lda, lp_ldb;
//     lp_lda=actNrParameter+actNrParameterC;
//     lp_ldb=actNrParameter+actNrParameterC;


//     // Solve symmetric linear system
//     //        zsysv_(&lp_matType, &lp_dim , &lp_nrRHS, lp_covVecf77, &lp_loadDim, &lp_interchanges,
//     //               lp_rhsVecf77, &lp_loadDimRHS, lp_workf77, &lp_lwork, &lp_info);
    
//     // Solves system with hermitian matrix
//     //     zhesv_(&lp_matType, &lp_dim , &lp_nrRHS, lp_covVecf77, &lp_loadDim, &lp_interchanges,
//     //            lp_rhsVecf77, &lp_loadDimRHS, lp_workf77, &lp_lwork, &lp_info);

//     // Solves general kind of linear system with pivotin
//     zgesv_(&lp_dim , &lp_nrRHS, lp_covVecf77, &lp_lda, lp_interchanges, lp_rhsVecf77, &lp_ldb, &lp_info);


//     //    std::cout<<"Finished solving ... the following lines got interchanged"<<std::endl;
//     for (UInt count =0 ;count < actNrParameter+actNrParameterC;count++)
//       std::cout<<(int)lp_interchanges[count]<<std::endl;

//     Vector<Complex> covVec2;
//     covVec2.Resize(actNrParameter+actNrParameterC*actNrParameter+actNrParameterC);


//     // reconvert: fortran -> cfs++
//     for ( UInt count = 0; count < actNrParameter+actNrParameterC*actNrParameter+actNrParameterC; count++ ) 
//       F772CC( lp_rhsVecf77[count], lp_rhsVec[count] );

//     for ( UInt count = 0; count < actNrParameter+actNrParameterC*actNrParameter+actNrParameterC; count++ ) 
//       F772CC( lp_covVecf77[count], lp_covVec[count]);

//     for (UInt count=0; count < lp_work.GetSize();count++)
//       F772CC(lp_workf77[count], lp_work[count]);

//     //    std::cout<<"Conversion f772cc Solved ZSYSV successfull ...?? " <<std::endl;

//     Matrix<Complex> inverse(actNrParameter+actNrParameterC, actNrParameter+actNrParameterC);
//     for (UInt i=0;i<actNrParameter+actNrParameterC;i++)
//       for (UInt j=0;j<actNrParameter+actNrParameterC;j++){
//         inverse[i][j]=lp_rhsVec[i+j*actNrParameter+actNrParameterC];
//       }
//     // std::cout<<"inverse:"<<std::endl;
// //     std::cout<<inverse<<std::endl;

// //    std::cout<<" Covariance  - Inverse : " <<std::endl;

// //     for(UInt i=0;i<inverse.GetSizeRow();i++)
// //       for (UInt j=0; j<inverse.GetSizeCol();j++){
// //         std::cout<<inverse[i][j].real()<<"+"<<inverse[i][j].imag()<<"i"<< ", ";
// //         //std::cout<<"F'("<<i<<")("<<j<<")= "<< JacobiMatrix[i][j]<<"; \t";
// //         if (j==inverse.GetSizeCol()-1)
// //           std::cout<<";\n";
// //       }
    
    data=rhsMat;

    //Matrix<Complex> covInv;
//     covInv=cov*inverse;
//     std::cout<<"Inv(cov)*cov:"<<std::endl;
//     std::cout<<inverse*cov<<std::endl;

//     for(UInt i=0;i<covInv.GetSizeRow();i++)
//       for (UInt j=0; j<covInv.GetSizeCol();j++){
//         std::cout<<covInv[i][j].real()<<"+"<<covInv[i][j].imag()<<"i"<< ", ";
//         if (j==covInv.GetSizeCol()-1)
//           std::cout<<";\n";
//       }

 //    std::cout<<"deleting variables in inversion of covariance Matrix ... " <<std::endl;

//     DeleteArray(lp_rhsVecf77);
//     DeleteArray(lp_interchanges);
//     DeleteArray(lp_covVecf77);
//     DeleteArray(lp_workf77 );

//     std::cout<<"leaving inversion of covariance Matrix ... " <<std::endl;

#endif


  } // end invertWithLapack

  void piezoParamIdent::invert(Matrix<Complex> & data)  {
//     data[0][0]=1;
//     data[0][1]=2;
//     data[0][2]=3;
//     data[1][0]=2;
//     data[1][1]=3;
//     data[1][2]=1;
//     data[2][0]=1;
//     data[2][1]=1;
//     data[2][2]=1;
    Matrix<Complex> dataTemp;
    dataTemp.Resize(actNrParameter+actNrParameterC,actNrParameter+actNrParameterC);
    dataTemp=data;

    Vector<Complex> e;
    Vector<Complex> x;
    e.Resize(actNrParameter+actNrParameterC);
    x.Resize(actNrParameter+actNrParameterC);
    Matrix<Complex> inverse;
    inverse.Resize(actNrParameter+actNrParameterC, actNrParameter+actNrParameterC);

    for (UInt ind=0;ind<data.GetSizeRow();ind++){
      e[ind]=1.0;
//       e[ind+1]=1.0;
//       e[ind+2]=1.0;
//      std::cout<<e<<std::endl;
      //      std::cout<<data.GetSizeRow()<<"x"<<data.GetSizeCol()<<std::endl;
      dataTemp.DirectSolve(x,e);
//       std::cout<<dataTmep<<std::endl;
//       std::cout<<data<<std::endl;
      dataTemp=data;
      //std::cout<<x<<std::endl;
      //      getchar();
      for (UInt indC=0;indC<data.GetSizeCol();indC++)
        inverse[indC][ind]=x[indC];
      x.Resize(actNrParameter+actNrParameterC);
      e[ind]=0.0;
    }
//     std::cout<<" inverse in invert: " <<std::endl;
//     std::cout<<inverse<<std::endl;
//     std::cout<<" data * inverse " <<std::endl;
//     td::cout<<data*inverse<<std::endl;
//      getchar();
     data=inverse;

  }
  void piezoParamIdent::createJacobian(Vector<Complex> & jacobi, Double omega){

    ENTER_FCN("piezoParamIdent::createJacobian");
    //    std::cout<<"create Jacobian .."<<std::endl;

    Vector<Double>parIncr1, parIncr2;
    parIncr1.Resize(nrParameter);
    parIncr2.Resize(nrParameter);
    Vector<Double>parIncrC1, parIncrC2;
    parIncrC1.Resize(nrParameter);
    parIncrC2.Resize(nrParameter);
    parIncr1=parameter;
    parIncr2=parameter;
    parIncrC1=parameterC;
    parIncrC2=parameterC;
    Complex F_hat_incr1,F_hat_incr2;
    Complex F_hat_incrC1,F_hat_incrC2;
    Integer parInd=0;
    jacobi.Resize(actNrParameter+actNrParameterC);

    MaterialData * ptMaterial=ptMyPDE_->getPDEMaterialData();   // Pointer to MaterialData

     for (UInt ind_param=0;ind_param<nrParameter;ind_param++){ 
      if (whichParameterToUpdate[ind_param]==1){
        
        parIncr1[ind_param]=1.001*parameter[ind_param];
        //        std::cout<<parIncr1<<std::endl;
        updateMaterialData(parIncr1,ptMaterial);
        createFVec(F_hat_incr1,FALSE,omega);
        // std::cout<<"par="<<ind_param<<" ="<<F_hat_incr1<<std::endl;
        
        parIncr2[ind_param]=0.999*parameter[ind_param];  

        updateMaterialData(parIncr2,ptMaterial);
        createFVec(F_hat_incr2,FALSE,omega);

//         std::cout<<"scaling"<<std::endl;
//         std::cout<<scaling<<std::endl;
        jacobi[parInd]=0.5*(F_hat_incr1-F_hat_incr2)/
          ((parIncr1[ind_param]-parameter[ind_param])*scaling[ind_param]);
        
        parIncr1[ind_param]=parameter[ind_param];
        parIncr2[ind_param]=parameter[ind_param];
        parInd++;

      }
    }
     parInd=0;
     for (UInt ind_param=0;ind_param<nrParameter;ind_param++){ 
      if (whichParameterToUpdateC[ind_param]==1){
        
        parIncrC1[ind_param]=1.000001*parameterC[ind_param];
        //      std::cout<<parameter_incr<<std::endl
        updateMaterialData(parIncrC1,ptMaterial);
        createFVec(F_hat_incrC1,FALSE,omega);
        
        parIncrC2[ind_param]=0.999999*parameterC[ind_param];  

        updateMaterialData(parIncrC2,ptMaterial);
        createFVec(F_hat_incrC2,FALSE,omega);

        jacobi[actNrParameter+parInd]=0.5*(F_hat_incrC1-F_hat_incrC2)/
          ((parIncrC1[ind_param]-parIncrC2[ind_param])*scalingC[ind_param]);
        
        parIncrC1[ind_param]=parameterC[ind_param];
        parIncrC2[ind_param]=parameterC[ind_param];
        parInd++;

      }
    }
     // std::cout<<"jacobi"<<std::endl;
//      std::cout<<jacobi<<std::endl;
//      getchar();
  } // end create Jacobian

  void piezoParamIdent::createFVec(Complex & F_hat, Boolean typeOut,
                                   Double frequency){
    ENTER_FCN("PiezoParamIdent:createFVec");
    //    std::cout<<"createFVec ...."<<std::endl;
        
    ptMaterial=ptMyPDE_->getPDEMaterialData();   // Pointer to MaterialData
    ptAssemble = ptMyPDE_->getPDE_assemble();
    ptAlgsys = ptMyPDE_->getPDE_algsys();

    Boolean reset = TRUE;
      
    ptAlgsys->InitMatrix();
    ptAssemble->SetReassemble();
    ptAlgsys->InitRHS();

  
   
      ////////////////////////////////////////////////////////
      //                   SOLVES PDE                      //
      ///////////////////////////////////////////////////////  

        //      ptMyPDE_->WriteGeneralPDEdefines();   // should not be used, overwrites to much!!    

      reset=TRUE;
      ptMyPDE_->GetSolveStep()->SetActFreq(frequency); 
      //      std::cout<<"\n piezoParam:createF PreStepHarmonic"<<std::endl;
      ptMyPDE_->GetSolveStep()->SetActStep(0);       
      //      std::cout<<"\n piezoParam:createF SetActStepHarmonic"<<std::endl;
      ptMyPDE_->GetSolveStep()->PreStepHarmonic(0); 
      //      std::cout<<"\n piezoParam:createF PreStepHarmonic"<<std::endl;        
      ptMyPDE_->GetSolveStep()->SolveStepHarmonic(reset);
      //      std::cout<<"\n piezoParam:createF SolveStepHarmonic"<<std::endl;

      ptMyPDE_->GetSolveStep()->PostStepHarmonic(reset);
      //      std::cout<<"\n piezoParam:createF PostStepHarmonic"<<std::endl;

      ptMyPDE_->PostProcess();


      //////////////////////////////////////////////////////////
      //Retrieves & stores Solution for further calculations  //
      /////////////////////////////////////////////////////////

        BaseNodeStoreSol * ptSol = ptMyPDE_->getPDESolution();
        NodeStoreSol<Complex> * ptNodeStoreSol;
        ptNodeStoreSol = dynamic_cast<NodeStoreSol<Complex>*>(ptSol);     

        Vector<Complex> chargeVec =   ptMyPDE_->getPDE_complexValuedCharge(); // Vector wich contains charges for each element !

        Complex charge=Complex(0.0,0.0);

         
        for (UInt i=0;i<chargeVec.GetSize();i++){
          charge+=chargeVec[i];
        }
        Integer fstep=0;
        Double x=real[fstep]*cos(PI/180*imag[fstep]);
        Double y=real[fstep]*sin(PI/180*imag[fstep]);
        Complex Z=Complex(x,y);

	// Logarithmic value of F
        F_hat=(sign*charge*Z)/std::log(Z); // without minus --- classical way ...
        //F_hat=charge;
        //	F_hat[fstep]=-(sign*charge*Z)/std::log(Z); // the "-" became important within the stack


        if (typeOut==TRUE){
          //      std::cout<<"\nFinished to create F ... here it is:"<<std::endl;

            std::cout<<"F(p)="<<F_hat<<"; \t";
          std::cout<<"\n ------------------------------- " <<std::endl;
      
          
    }

  } // end createF
  void piezoParamIdent::readInMeasurement(Vector<Double> & frequencies){

    ENTER_FCN("piezoParamIdent::readInMeasurement");

    Char* measuremets="mess.dat";
    mess = new std::ifstream(measuremets, std::basic_ios<char>::in);
    if (!mess){
      std::cerr << "\n File measuredData.dat does not exist!" << std::endl;
    }

    frequencies.Resize(nrMeasuredData);

    real.Resize(nrMeasuredData);
    imag.Resize(nrMeasuredData);

    char mDataRow[256], helpChar[64];
    UInt i=0, j=0;
    Double newFreq, amplitude,phase;

    while(mess->getline(mDataRow, 265)){
      //std::cout<<mDataRow<<std::endl;
      i=0;j=0;
      while (mDataRow[i]!='\t'){
        helpChar[j]=mDataRow[i];
        i++;j++;
      }// end while madataRow
      newFreq=atof(helpChar);
//       std::cout<<"newFreq"<<std::endl;
//       std::cout<<newFreq<<std::endl;
      i++;
      j=0;
      for(UInt k=0;k<i;k++) // Delete content of helpChar
        helpChar[k]='0';
      while (mDataRow[i]!='\t'){
        helpChar[j]=mDataRow[i];
        i++;j++;
      }// end while madataRow
      amplitude=atof(helpChar);
      //      std::cout<<amplitude<<std::endl;
      i++;
      j=0;
      for(UInt k=0;k<i;k++) // Delete content of helpChar
        helpChar[k]='0';
      while (mDataRow[i]!='\t'){
        helpChar[j]=mDataRow[i];
        i++;j++;
      }// end while madataRow
      phase=atof(helpChar);
      //      std::cout<<phase<<std::endl;
//       std::cout<<"freqs:"<<std::endl;
//       std::cout<<freqs<<std::endl;

      freqs[0]=startfreq;
      for(UInt mInd=0;mInd<nrMeasuredData;mInd++){
        if(std::abs(freqs[mInd]-newFreq)<std::abs(freqs[mInd]-frequencies[mInd])){
          frequencies[mInd]=newFreq;
          real[mInd]=amplitude;
          imag[mInd]=phase;
        }
      }
    }// end while mess
    //  std::cout<<frequencies<<std::endl;

    freqs.Resize(nrMeasuredData);
    freqs=frequencies;   
//     std::cout<<" New Measurements"<<std::endl;
//     for(UInt m=0;m<nrMeasuredData;m++)
//       std::cout<<freqs[m]<<",  " << real[m] << ",  " << imag[m] <<std::endl;

    mess->close();

  }// end readInMeasurements

      
} // end namespace
