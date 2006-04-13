#include "PDE/SinglePDE.hh" 
#include "PDE/piezoPDE.hh"
#include "piezoParamIdent.hh"



namespace CoupledField
{

  // xxxxxxxxxxxxxxxxxxxxxxxxxx nuMethods  - ITERATION xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  void piezoParamIdent::optimalExpDesign(){

    ENTER_FCN("piezoParamIdent::optimalDesign");

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
    nrMeasuredData=15;
    Vector<Double> freqs5;
    freqs5.Resize(15);
    freqs.Resize(15);
   
  
   
    //initial guesses suggested by optExpDesignVarNr ...for thickness mode
    // Thicknes = 0.1mm
    if (FALSE){
      freqs5[0]=1.1e+06;
      freqs5[1]=1.3e+06;
      freqs5[2]=1.4e+06;
      freqs5[3]=1.5e+06;
      freqs5[4]=1.6e+06;
      freqs5[5]=1.7e+06;
      freqs5[6]=1.8e+06;
      freqs5[7]=2.5e+06;
      freqs5[8]=2.7e+06;
      freqs5[9]=3.0e+06;
      freqs5[10]=3.2e+06;
      freqs5[11]=3.4e+06;
      freqs5[12]=3.6e+06;
      freqs5[13]=3.8e+06;
    }

    //initial guesses suggested by optExpDesignVarNr ...for thickness mode
    // Thicknes = 0.05mm
    if (FALSE){
      freqs5[0]=2550000;
      freqs5[1]=2733330;
      freqs5[2]=2916670;
      freqs5[3]=3100000;
      freqs5[4]=3283330;
      freqs5[5]=3466670;
      freqs5[6]=3650000;
      freqs5[7]=3690500;
      freqs5[8]=3705000;
      freqs5[9]=5200000;
      freqs5[10]=5236880;
      freqs5[11]=5300000;
      freqs5[12]=5483330;
      freqs5[13]=6.0e+06;
    }


    // Fitting an CeramTecMaterial, 5V Vorspannung, Radius 11mm, Dicke 1mm
    if (FALSE){
      freqs5[0]=0.2e+05;
      freqs5[1]=0.35e+05;
      freqs5[2]=0.45e+05;
      freqs5[3]=0.55e+05;
      freqs5[4]=0.6e+05;
      freqs5[5]=0.7e+05;
      freqs5[6]=0.8e+05;
      freqs5[7]=1.15e+05;
      freqs5[8]=1.25e+05;
      freqs5[9]=1.3e+05;
      freqs5[10]=1.4e+05;
    }

    // Fitting an CeramTecMaterial, 5V Vorspannung, Radius 11mm, Dicke 1mm
    if (TRUE){
      freqs5[0]=0.2e+05;
      freqs5[1]=0.35e+05;
      freqs5[2]=0.45e+05;
      freqs5[3]=0.55e+05;
      freqs5[4]=0.6e+05;
      freqs5[5]=0.69e+05;
      freqs5[6]=0.72e+05;
      freqs5[7]=0.82e+05;
      freqs5[8]=1.15e+05;
      freqs5[9]=1.25e+05;
      freqs5[10]=1.3e+05;
      freqs5[11]=1.4e+05;
      freqs5[12]=1.6e+05;
      freqs5[13]=1.8e+05;
      freqs5[14]=2.0e+05;
    }


    //initial guesses for planar, radial mode
    if (FALSE){
      freqs5[0]=0.31e+5;
      freqs5[1]=0.3e+5;
      freqs5[2]=0.45e+5;
      freqs5[3]=0.55e+5;
      freqs5[4]=0.65e+5;
      freqs5[5]=0.7e+5;
      freqs5[6]=0.8e+5;
      freqs5[7]=0.85e+5;
      freqs5[8]=1.45e+5;
      freqs5[9]=1.5e+5;
      freqs5[10]=1.6e+5;
      freqs5[11]=1.65e+5;
      freqs5[12]=1.7e+5;
      freqs5[13]=1.75e+5;
      freqs5[14]=1.8e+5;
    }
      

    //   freqs5[0]=1.1e+06;
    //     freqs5[1]=1.3e+06;
    //     freqs5[2]=1.5e+06;
    //     freqs5[3]=1.6e+06;
    //     freqs5[4]=1.7e+06;
    //     freqs5[5]=1.8e+06;
    //     freqs5[6]=2.5e+06;
    //     freqs5[7]=2.7e+06;
    //     freqs5[8]=3.0e+06;
    //     freqs5[9]=3.1e+06;
    //     freqs5[10]=3.2e+06;
    //     freqs5[11]=3.4e+06;
    //     freqs5[12]=3.6e+06;
    //     freqs5[13]=3.8e+06;
    //     freqs5[14]=3.95e+06;

    freqs=freqs5;

    //Writes out initial guesses of optimalFreqs

    *optimalFreqs<< 0 <<"  "<< 0 <<"  ";      
    
    for(UInt fr=0;fr<nrMeasuredData;fr++)
      *optimalFreqs<< freqs[fr]<<"  ";
    *optimalFreqs<<std::endl;         
   

    projGradientFlags.Resize(nrMeasuredData);
    parameterInitial = parameter;
    residuumParIdentOld = residuumParIdent=5.0e-3;
       

    std::cout<<"++ Optimum experiment desing with "<<nrMeasuredData <<
      " frequencies and "<<actNrParameter+actNrParameterC <<
      " parameter ..." << std::endl;
    Vector<Double> newFreqs;
   
    readInMeasurement(newFreqs);
    //     newFreqs[0]=newFreqs[1];
    //     real[0]=real[1];
    //     imag[0]=imag[1];
    //     freqs[0]=freqs[1];
    std::cout<<"++ Frequencies taken from file mess.dat:"<<std::endl;
    std::cout<<newFreqs<<std::endl;

//     for(UInt fr=0;fr<nrMeasuredData;fr++)
//       *impedCurve<< freqs[fr]<<"  ";
//     *impedCurve<<std::endl;

    calc_measuredCharge(freqs, real, imag, y_hat); // out of new measurements

    //  actNrParameterC=0;
    MaterialMap ptMaterial;

    if(directCoupling==TRUE)
      ptMaterial=ptPDE1_->getPDEMaterialData();   // Pointer to MaterialData
    else
      ptMaterial=ptMyPDE_->getPDEMaterialData();   // Pointer to MaterialData

    // Frequency bounds

    // fr = 4000m/s / 2*thickness

  
    Double rho;
    Vector<Complex> jacobi;
    Complex functional;

    ptMaterial.begin()->second->GetScalar(rho,DENSITY,REAL);

    // fa_ = 1.0/(2*thickness) * std::sqrt(1.3e+11/rho);

    fa_=antiResonanceFrequency_;
    fa_=1.1*fa_;

    fr_=resonanceFrequency_;
    fr_=0.9*fr_;
    // fa = 1/thickness * sqrt(c_33/rho)

    std::cout<<"++ Bounds for resonance and antiresonace frequency: "<<std::endl;
    std::cout<<"fr= "<< fr_ <<std::endl;
    std::cout<<"fa= "<< fa_ <<std::endl;

    for (UInt nrOptExpSteps=0; nrOptExpSteps<20; nrOptExpSteps++){ 
      UInt nrNuMethods=0;
      // That at least we perform one descent step
      normGradient=2.0;
      normGradientOld=1.0;

      //      descentMethod(functional);

      readInMeasurement(newFreqs);

      std::cout<<"++ Frequencies, just taken from mess.dat - file ... " <<std::endl;
      std::cout<<newFreqs<<std::endl;

      for(UInt fr=0;fr<nrMeasuredData;fr++)
        *impedCurve<< freqs[fr]<<"  ";
      *impedCurve<<std::endl;

      calc_measuredCharge(freqs, real, imag, y_hat); // out of new measurements
      newtonCounter=0;

      while (nrNuMethods<maxNumberNewtonLoops && 
             residuumParIdent > 0.05*residuumParIdentOld){
        //while (nrNuMethods<maxNumberNewtonLoops){
        ///        nuMethodsC2();
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

      Vector<Double> saveFreqs=freqs;

      Double startFreqTemp;
      startFreqTemp=startfreq;
      Double freqincr=(stopfreq-startfreq)/50;
      freqs.Resize(50);
      for(UInt i=0;i<50;i++){
        startFreqTemp+=freqincr;
        freqs[i]=startFreqTemp;
      }
      calcImpedanceCurve();
      freqs=saveFreqs;

      fa_=antiResonanceFrequency_;
      fa_=1.1*fa_;
      
      fr_=resonanceFrequency_;
      fr_=0.9*fr_;

      std::cout<<"++ New bounds for resonance and antiresonace frequency: "<<std::endl;
      std::cout<<"fr= "<< fr_ <<std::endl;
      std::cout<<"fa= "<< fa_ <<std::endl;


      residuumParIdentOld = residuumParIdent;
      
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
    Double lambda;
    UInt descIter=0;
    Double gradientNormUpdate;
    // fitting am Dicken Mode!
    if (actNrParameter==3)
      lambda = 0.00002; // for three parameter
    else if (actNrParameter==10||actNrParameter==9)
      lambda=0.35;
    //      lambda= 0.002; // for 10 parameters with StepWidth
    // lambda = 0.000002; // for nine parameter
    else if (actNrParameter==6)
      lambda=0.025;

    // Fitting mit thickness = 0.05 mm
    lambda = 3.500e-03;

    // Fitting am radial, planar Mode
    //       lambda = 1.0e-14;
    //        Double lambda = 1.0e-5; // for ten parameter
    //    Double lambda = 1.0; // for nine parameter using M - confidence crit.
    //    Double lambda = 1.0e+7; // for nine parameter using M - confidence crit.
    //Double lambda = 40; // for three parameters using M - confidence crit.


    UInt maxNumberDescentIterations=5;
    gradientNormUpdate = std::abs(normGradient - normGradientOld)/normGradientOld;
    std::cout<<"gradientNormUpdate before descent iteration"<<std::endl;
    std::cout<<gradientNormUpdate<<std::endl;

   
    while((descIter<maxNumberDescentIterations && gradientNormUpdate>0.1)|| descIter<2){
      descIter++;
      //    projGradientFlags.Resize(nrMeasuredData);
      Complex J, JOld;
      Vector<Double> freqsOld;
      freqsOld.Resize(nrMeasuredData);
      freqsOld=freqs;
      createCovA(J,FALSE);
      JOld=J;

      if (JOld.real()<=0.0){
        std::cout<< "! ! ! we should be in trouble here, since Re(J)<0 "<<std::endl;
        std::cout<<JOld.real()<<std::endl;
        //        getchar();
      }

      normGradientOld=normGradient;
      //overwrites normGradient!
      createGradient(gradient,1.0e-5);
      gradientNormUpdate = std::abs(normGradient - normGradientOld)/normGradientOld;

      //       std::cout<<"gradientNormOld"<<std::endl;
      //       std::cout<<normGradientOld<<std::endl;
      //       std::cout<<"gradientNorm"<<std::endl;
      //       std::cout<<normGradient<<std::endl;

      //       std::cout<<"gradientNormUpdate"<<std::endl;
      //       std::cout<<gradientNormUpdate<<std::endl;

      lambda = 1*lambda;

      std::cout<<"gradient"<<std::endl;
      std::cout<<gradient<<std::endl;
    
      UInt innerCounter=0;
      //      while(innerCounter<1){
      while(JOld.real()<=J.real()){      
        for(UInt fr=0;fr<nrMeasuredData;fr++)
          freqs[fr]=freqs[fr]-lambda*gradient[fr].real();
    
        for(UInt fr=0;fr<nrMeasuredData;fr++)
          if (freqs[fr]<=0.0){
            std::cout<<"! Frequency "<< fr << " is negative " <<std::endl;
            std::cout<<"! Frequency "<< fr << " will be reset to previous value " <<std::endl;
            freqs[fr]=freqs[fr]+lambda*gradient[fr].real();
          }
        
        std::cout<<"freqs after update by -lambda*GradJ:"<<std::endl;
        std::cout<<freqs<<std::endl;
      
        createCovA(J,FALSE);

        
        if (J.real()>JOld.real()){
          lambda=0.1*lambda;
          freqs=freqsOld;
          std::cout<<"! Reduced lambda in descent method ... lambda =  " << lambda<<std::endl;
          if (J.real()>=JOld.real())
            std::cout<<"  -> because JOld = " << JOld.real()<<" <= JNew = " <<J.real() << std::endl;
          innerCounter ++;
          if (innerCounter>10)
            break;
        
        }
        else 
          std::cout<< "JNew is smaller than JOld ... continue with lambda = " << lambda <<std::endl;
        
        for(UInt fr=0;fr<nrMeasuredData;fr++){
          
          if(TRUE){
            if ((freqs[fr]<=fa_)&&(freqs[fr]>=fr_)){
              if(freqs[fr]<(fr_+fa_)/2){
                freqs[fr]=fr_;
                projGradientFlags[fr]=freqs[fr];
                std::cout<<"! Frequency "<< fr <<" to close to fr, fixed it at f= " << freqs[fr]<< std::endl;
                if (fr>=1)
                  if (freqs[fr]==freqs[fr-1]){
                    freqs[fr]=0.95*freqs[fr];
                    std::cout<<"Frequency "<< fr <<" (fr) coincides with previos one( " << freqs[fr-1] <<" ), set freq= " << freqs[fr]<<std::endl;
                  }
                projGradientFlags[fr]=freqs[fr];
              }
              else{
                freqs[fr]=fa_;
                projGradientFlags[fr]=freqs[fr];
                std::cout<<"! Frequency "<< fr <<" to close to fa, fixed it at f= " << freqs[fr]<< std::endl;
                if (fr>=1)
                  if (freqs[fr]==freqs[fr-1]){
                    freqs[fr]=1.05*fa_;
                    std::cout<<"Frequency "<< fr <<" (fr) coincides with previos one, set freq= " << freqs[fr]<<std::endl;
                  }
              }
            } 
          } // end if FALSE/TRUE
          
          if (freqs[fr]< startfreq){
            std::cout<<"! There are no measurements for frequency ("<<fr<<") = " << freqs[fr]  <<std::endl;
            std::cout<<"! Fixed frequency at " <<startfreq<< std::endl;
            freqs[fr]=startfreq;
            if (fr>=1)
              if (freqs[fr]==freqs[fr-1]){
                freqs[fr]=startfreq+fr*0.1*startfreq;
                std::cout<<"! Frequency "<< fr <<" to close to startfreq and coincides with previous one, fixed it at f= " << freqs[fr]<< std::endl;
              }
            projGradientFlags[fr]=freqs[fr];
          }
          if(freqs[fr]>stopfreq){
            std::cout<<"! There are no measurements for frequency (" << fr<< ") = " <<freqs[fr] <<std::endl;
            std::cout<<"! Fixed frequency at  " <<stopfreq<<std::endl;
            freqs[fr]=stopfreq;
            if (fr>=1)
              if (freqs[fr]==freqs[fr-1]){
                freqs[fr]=stopfreq-fr*0.1*startfreq;
                std::cout<<"! Frequency "<< fr <<" to close to stopfreq and coincides with previous one, fixed it at f= " << freqs[fr]<< std::endl;
              }

            projGradientFlags[fr]=freqs[fr];
          }
            
            
        }// end in freqs[fr]<=fa_ ...

        //   std::cout<<" frequencies before final check of coincidences:" <<std::endl;
        //       std::cout<<freqs<<std::endl;


        // check if two frequencies coincide:
        // if they do so reduce them by 5 per cent

        Boolean flag=TRUE;

        while(flag){
          for(UInt i=0;i<freqs.GetSize();i++){
            for(UInt j=0;j<freqs.GetSize();j++){
              if((freqs[i]==freqs[j]) && (i!=j))
                if (i<freqs.GetSize()-1)
                  freqs[j]=1.025*freqs[j];
                else
                  freqs[j]=0.975*freqs[j];
              else
                flag=FALSE;
            }
            // check if frequencies are out of range        
            while(freqs[i]>stopfreq)
              freqs[i]=0.98*freqs[i];
            while(freqs[i]<startfreq)
              freqs[i]=1.02*freqs[i];             
          }
        }
         
        createCovA(J,FALSE);         
        // end while ....         
      }      

      // will write also confidence intervals if we found a smaller  J
      // if JNew>Jold after 10 inner loops, we might assume, that there is no improvement for 
      // the location of the frequencies and therefore write nothing into confidence intervals
      if (innerCounter<10)
        createCovA(J,TRUE);
      std::cout<<" Value J(w) calculated with projected gradient "<< J <<std::endl;
      std::cout<<" frequencies at end of descent Method:" <<std::endl;
      std::cout<<freqs<<std::endl;

      *optimalFreqs<< descIter <<"  "<< J.real() <<"  ";      

      for(UInt fr=0;fr<nrMeasuredData;fr++)
        *optimalFreqs<< freqs[fr]<<"  ";
      *optimalFreqs<<std::endl;               
    } // while descIter<max Number


  }// descentMethod


  void piezoParamIdent::createCovA(Complex &J, Boolean writeOutCov){
    ENTER_FCN("piezoParamIdent::createCovA");

    //    Complex J;
    J=Complex(0.0,0.0);
    Double deltaLocal=1.01;

    Vector<Complex> jacobi;
    Vector<Complex> jacobiH;
    //    MaterialData * ptMaterial=ptMyPDE_->getPDEMaterialData();   // Pointer to MaterialData
    Complex F_hat_incr1,F_hat_incr2;

    Matrix<Complex> cov(actNrParameter+actNrParameterC, actNrParameter+actNrParameterC);
    Matrix<Complex> covTemp(actNrParameter+actNrParameterC, actNrParameter+actNrParameterC);

    for (UInt actFreq=0;actFreq<nrMeasuredData;actFreq++){

      Double stepWidth=0.0;
      stepWidth=0.5*std::abs((freqs[std::min(actFreq+1,nrfreq-1)]-freqs[std::max((Integer)actFreq-1,0)]));
       
      createJacobian(jacobi, freqs[actFreq]);
      jacobiH.Resize(actNrParameter+actNrParameterC);
       
      for (UInt i=0;i<jacobi.GetSize();i++)
        jacobiH[i]=Complex(jacobi[i].real(),-jacobi[i].imag());
       
      for(UInt i=0;i<actNrParameter+actNrParameterC;i++)
        for(UInt j=0;j<actNrParameter+actNrParameterC;j++){
          //           ovTemp[i][j]=jacobiH[i]*jacobi[j]/(delta*delta*std::abs(y_hat[actFreq])*std::abs(y_hat[actFreq]));
          covTemp[i][j]= stepWidth*jacobiH[i]*jacobi[j]/
            ((deltaLocal)*(deltaLocal)*std::abs(y_hat[actFreq])*std::abs(y_hat[actFreq]));
           
          //    if (i==j)
          //    covTemp[i][i]=covTemp[i][i]+covTemp[i][i];// *1.0e-10;
        }
      cov+=covTemp;
         
    } // end for actFreq ...

    for(UInt i=0; i<actNrParameter+actNrParameterC; i++)
      if(cov[i][i].real()<0.0){
        std::cout<<"WARNING: Entries "<< i  <<" on diagonal in (F'*F) are negative"<<std::endl;
        std::cout<<cov<<std::endl;
        std::cout<<jacobi<<std::endl;
        std::cout<<jacobiH<<std::endl;
      }

       
    //     invertWithLapack(cov);

    Matrix<Complex> data, data1;
    //     data.Resize(actNrParameter+actNrParameterC,actNrParameter+actNrParameterC);
    data1=data=cov;     

    //      std::cout<<"cov:"<<std::endl;
    //      std::cout<<cov<<std::endl;
    //      getchar();

    invert(cov);

    // std::cout<<"cov:"<<std::endl;
    //     std::cout<<cov<<std::endl;
    //     getchar();

   
    data1=data*cov;
    Complex diagC=Complex(0.0,0.0);
    for(UInt i=0;i<actNrParameter+actNrParameterC;i++)
      diagC+=data1[i][i];
    diagC/=Complex(actNrParameter+actNrParameterC,0.0);
    //     std::cout<<"Trace of inv(cov)*cov = "<< diagC <<std::endl;
    if (std::abs(diagC)>1.5)
      std::cout<<" ! Inversion of Cov failed!!"<<std::endl;

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

    // A - Criterion for all parameters
    if(TRUE){
      Vector<Double> covDiag(actNrParameter+actNrParameterC);
      for(UInt i=0; i<actNrParameter+actNrParameterC; i++){
        J+=cov[i][i];
        covDiag[i]=cov[i][i].real();
        if(cov[i][i].real()<0.0)
          std::cout<<"WARNING: Entries "<< i <<" on diagonal in cov =inv(F'*F) are negative"<<std::endl;
      }

      Double JReal = J.real();
      JReal = JReal/actNrParameter+actNrParameterC;
      J=J/Complex(actNrParameter+actNrParameterC,0.0);

      //optimality criterion just for one parameter, namely eps_11!
      if(FALSE){
        J=cov[7][7].real();
      }
      //optimality criterion just for four parameter, namely eps_11, e15, e13, c13!
      if(FALSE){
        J=cov[7][7].real()+cov[2][2].real()+cov[4][4].real()+cov[5][5].real();
        J=J/Complex(4.0,0.0);
      }


         
      //      writeOutCov=TRUE;
      if (writeOutCov==TRUE){
   
        if (actNrParameter==3){      
          //          *confInterval<<parameter[1]+parameterInitial[1]*std::sqrt(cov[0][0].real()*0.115*0.115)<<"  ";
          *confInterval<<parameter[1]+parameterInitial[1]*(cov[0][0].real()*0.7)<<"  ";
          //        std::cout<<parameter[1]+parameter[1]*std::sqrt(cov[1][1].real()*0.115*0.115)<<std::endl;
          *confInterval<<parameter[1]<<"  ";
          //          *confInterval<<parameter[1]-parameterInitial[1]*std::sqrt(cov[0][0].real()*0.115*0.115)<<"  ";
          *confInterval<<parameter[1]-parameterInitial[1]*(cov[0][0].real()*0.7)<<"  ";
          //        std::cout<<parameter[1]-1000.0*std::sqrt(cov[0][0].real()*parameter[1]*0.115*0.115)<<std::endl;


          //           *confInterval<<parameter[7]+parameterInitial[7]*std::sqrt(cov[1][1].real()*0.115*0.115)<<"  ";
          //           *confInterval<<parameter[7]<<"  ";
          //           *confInterval<<parameter[7]-parameterInitial[7]*std::sqrt(cov[1][1].real()*0.115*0.115)<<"  ";

          *confInterval<<parameter[7]+parameterInitial[7]*(cov[1][1].real()*0.7)<<"  ";
          *confInterval<<parameter[7]<<"  ";
          *confInterval<<parameter[7]-parameterInitial[7]*(cov[1][1].real()*0.7)<<"  ";


          //           *confInterval<<parameter[9]+parameterInitial[9]*std::sqrt(cov[2][2].real()*0.115*0.115)<<"  ";
          //           *confInterval<<parameter[9]<<"  ";
          //           *confInterval<<parameter[9]-parameterInitial[9]*std::sqrt(cov[2][2].real()*0.115*0.115)<<"  ";

          *confInterval<<parameter[9]+parameterInitial[9]*(cov[2][2].real()*0.7)<<"  ";
          *confInterval<<parameter[9]<<"  ";
          *confInterval<<parameter[9]-parameterInitial[9]*(cov[2][2].real()*0.7)<<"  ";

        }

        else if (actNrParameter==4){
          for (UInt ii=0;ii<=4;ii++){
            if(ii<=1){
              //               *confInterval<<parameter[ii]+parameter[ii]*std::sqrt(cov[ii][ii].real()*0.115*0.115)<<"  ";
              //               *confInterval<<parameter[ii]<<"  ";
              //               *confInterval<<parameter[ii]-parameter[ii]*std::sqrt(cov[ii][ii].real()*0.115*0.115)<<"  ";

              *confInterval<<parameter[ii]+parameterInitial[ii]*(cov[ii][ii].real()*0.7)<<"  ";
              *confInterval<<parameter[ii]<<"  ";
              *confInterval<<parameter[ii]-parameterInitial[ii]*(cov[ii][ii].real()*0.7)<<"  ";

            }
            else if(ii>2)
              {
                //                 *confInterval<<parameter[ii]+parameterInitial[ii]*std::sqrt(cov[ii-1][ii-1].real()*0.115*0.115)<<"  ";
                //                 *confInterval<<parameter[ii]<<"  ";
                //                 *confInterval<<parameter[ii]-parameterInitial[ii]*std::sqrt(cov[ii-1][ii-1].real()*0.115*0.115)<<"  ";

                *confInterval<<parameter[ii]+parameterInitial[ii]*(cov[ii-1][ii-1].real()*0.7)<<"  ";
                *confInterval<<parameter[ii]<<"  ";
                *confInterval<<parameter[ii]-parameterInitial[ii]*(cov[ii-1][ii-1].real()*0.7)<<"  ";

              }
          }
        }
        if (actNrParameter==10){
          for (UInt ii=0;ii<10;ii++){
            *confInterval<<parameter[ii]+parameterInitial[ii]*sqrt(cov[ii][ii].real()*2.5*2.5)<<"  ";
            *confInterval<<parameter[ii]<<"  ";
            *confInterval<<parameter[ii]-parameterInitial[ii]*sqrt(cov[ii][ii].real()*2.5*2.5)<<"  ";
          }
        }

        if (actNrParameter==9){
          for (UInt ii=0;ii<10;ii++){
            if(ii<=1){
              *confInterval<<parameter[ii]+parameterInitial[ii]*sqrt(cov[ii][ii].real()*2.5*2.5)<<"  ";
              //              std::cout<<parameter[ii]+parameterInitial[ii]*sqrt(cov[ii][ii].real()*2.5*2.5)<<"  ";
              *confInterval<<parameter[ii]<<"  ";
              //              std::cout<<parameter[ii]<<std::endl;
              *confInterval<<parameter[ii]-parameterInitial[ii]*sqrt(cov[ii][ii].real()*2.5*2.5)<<"  ";
              //              std::cout<<parameter[ii]-parameterInitial[ii]*sqrt(cov[ii][ii].real()*2.5*2.5)<<"  ";
              //              getchar();
            }
            else if (ii>2){
              *confInterval<<parameter[ii]+parameterInitial[ii]*sqrt(cov[ii-1][ii-1].real()*2.5*2.5)<<"  ";
              *confInterval<<parameter[ii]<<"  ";
              *confInterval<<parameter[ii]-parameterInitial[ii]*sqrt(cov[ii-1][ii-1].real()*2.5*2.5)<<"  ";
            }

          }
        }

        if (actNrParameter==5){
          for (UInt ii=0;ii<5;ii++){
            *confInterval<<parameter[ii]+parameterInitial[ii]*std::sqrt(cov[ii][ii].real()*0.55*0.55)<<"  ";
            *confInterval<<parameter[ii]<<"  ";
            *confInterval<<parameter[ii]-parameterInitial[ii]*std::sqrt(cov[ii][ii].real()*0.55*0.55)<<"  ";
          }
        }   

        if (actNrParameter==6){
          
          *confInterval<<parameter[0]+parameterInitial[0]*(cov[0][0].real()*0.03)<<"  ";
          *confInterval<<parameter[0]<<"  ";
          *confInterval<<parameter[0]-parameterInitial[0]*(cov[0][0].real()*0.03)<<"  ";

          *confInterval<<parameter[1]+parameterInitial[1]*(cov[1][1].real()*1.0)<<"  ";
          *confInterval<<parameter[1]<<"  ";
          *confInterval<<parameter[1]-parameterInitial[1]*(cov[1][1].real()*1.0)<<"  ";

          *confInterval<<parameter[3]+parameterInitial[3]*(cov[2][2].real()*0.03)<<"  ";
          *confInterval<<parameter[3]<<"  ";
          *confInterval<<parameter[3]-parameterInitial[3]*(cov[2][2].real()*0.03)<<"  ";

          *confInterval<<parameter[4]+parameterInitial[4]*(cov[3][3].real()*0.03)<<"  ";
          *confInterval<<parameter[4]<<"  ";
          *confInterval<<parameter[4]-parameterInitial[4]*(cov[3][3].real()*0.03)<<"  ";

          *confInterval<<parameter[7]+parameterInitial[7]*(cov[4][4].real()*0.75)<<"  ";
          *confInterval<<parameter[7]<<"  ";
          *confInterval<<parameter[7]-parameterInitial[7]*(cov[4][4].real()*0.75)<<"  ";

          *confInterval<<parameter[9]+parameterInitial[9]*(cov[5][5].real()*0.75)<<"  ";
          *confInterval<<parameter[9]<<"  ";
          *confInterval<<parameter[9]-parameterInitial[9]*(cov[5][5].real()*0.75)<<"  ";

        }
        *confInterval<<std::endl;
      }
    }
  

    
    // smallest eigenvalue criterion:
    if(FALSE){

      Vector<Double> eig;
      eig.Resize(actNrParameter+actNrParameterC);
#ifdef USE_LAPACK
      cov.eigenvaluesWithLapack(eig);
#endif 

      //  std::cout<<"eigenvalues"<<std::endl;
      //       std::cout<<eig<<std::endl;
      
      Double eigMin=std::abs(eig[0]);
      for(UInt eigInd=0; eigInd<actNrParameter+actNrParameterC; eigInd++){
        if(std::abs(eig[eigInd])>eigMin)
          eigMin=std::abs(eig[eigInd]);
      }
      J=std::log(eigMin);


    } // end minim EV criterion

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
      if (projGradientFlags[j]==0.0){
        freqs[j]=freqs[j]+dOmega*freqs[j];
        createCovA(J2,FALSE);
        for(UInt i=0; i<actNrParameter+actNrParameterC;i++)
          grad[j]= (J2-J1)* voltage*voltage / (dOmega*freqs[j]) ;
        //grad[j]= (J2-J1) / (dOmega*freqs[j]) ;
        freqs[j]=freqsOld[j];
      }

    }

    for (UInt i=0; i<nrMeasuredData;i++)
      normGradient+=grad[i].real()*grad[i].real();
    normGradient=std::sqrt(normGradient);


  } // end create Gradient

#ifdef USE_LAPACK
  void piezoParamIdent::invertWithLapack(Matrix<Complex> & data){
    ENTER_FCN("piezoParamIdent::invertWithLapack");
    
    std::cout<<"Optimum experiment design works with LAPACK Routines"<<std::endl;
    std::cout<<"Please set LAPACK = yes in your Makefile.option (CFS & OLAS)" <<std::endl;
    
    Matrix<Complex> cov;
    cov.Resize(actNrParameter+actNrParameterC,actNrParameter+actNrParameterC);
    cov=data;
     
    Matrix<Complex> rhsMat;
    rhsMat.Resize(actNrParameter+actNrParameterC, actNrParameter+actNrParameterC);
    for(UInt i=0;i<actNrParameter+actNrParameterC;i++)
      rhsMat[i][i]=Complex(1.0,0.0);
    lapackSysMatType LAPACK_SYS_MAT_TYPE = ZGESV;
    data.solveWithLapack(rhsMat,LAPACK_SYS_MAT_TYPE);
    data=rhsMat;
     
  } // end invertWithLapack
  
#endif 

  void piezoParamIdent::invert(Matrix<Complex> & data)  {

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

      dataTemp.DirectSolve(x,e);
      dataTemp=data;
      for (UInt indC=0;indC<data.GetSizeCol();indC++)
        inverse[indC][ind]=x[indC];
      x.Resize(actNrParameter+actNrParameterC);
      e[ind]=0.0;
    }
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

    MaterialMap ptMaterial;

    if(directCoupling==TRUE)
      ptMaterial=ptPDE1_->getPDEMaterialData();   // Pointer to MaterialData
    else
      ptMaterial=ptMyPDE_->getPDEMaterialData();   // Pointer to MaterialData
    
    for (UInt ind_param=0;ind_param<nrParameter;ind_param++){ 
      if (whichParameterToUpdate[ind_param]==1){
        
        parIncr1[ind_param]=1.001*parameter[ind_param];
        updateMaterialData(parIncr1,ptMaterial);
        createFVec(F_hat_incr1,FALSE,omega);
        
        parIncr2[ind_param]=0.999*parameter[ind_param];  
        updateMaterialData(parIncr2,ptMaterial);
        createFVec(F_hat_incr2,FALSE,omega);

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

  } // end create Jacobian

  void piezoParamIdent::createFVec(Complex & F_hat, Boolean typeOut,
                                   Double frequency){
    ENTER_FCN("PiezoParamIdent:createFVec");
    //    std::cout<<"createFVec ...."<<std::endl;
      

    if(directCoupling==TRUE){
      ptMaterial=ptPDE1_->getPDEMaterialData();   // Pointer to MaterialData
      ptAssemble = ptPDE1_->getPDE_assemble();
      ptAlgsys = ptPDE1_->getPDE_algsys();    

    }
    else{
      ptMaterial=ptMyPDE_->getPDEMaterialData();   // Pointer to MaterialData
      ptAssemble = ptMyPDE_->getPDE_assemble();
      ptAlgsys = ptMyPDE_->getPDE_algsys();    
      ptAlgsys->InitMatrix();
      ptAlgsys->InitRHS();
    }

    ptAssemble->SetReassemble();
   

    Boolean reset = TRUE;
      

    ////////////////////////////////////////////////////////
    //                   SOLVES PDE                      //
    ///////////////////////////////////////////////////////  
      
      reset=TRUE;
      ptPDE_->GetSolveStep()->SetActFreq(frequency); 
      ptPDE_->GetSolveStep()->SetActStep(0);       
      ptPDE_->GetSolveStep()->PreStepHarmonic(0); 
      ptPDE_->GetSolveStep()->SolveStepHarmonic(reset);
      ptPDE_->GetSolveStep()->PostStepHarmonic(reset);
      ptPDE_->PostProcess();

      //////////////////////////////////////////////////////////
      //Retrieves & stores Solution for further calculations  //
      /////////////////////////////////////////////////////////

        Vector<Complex> chargeVec;
        if(directCoupling==TRUE)      
          chargeVec = ptPDE1_->getPDE_complexValuedCharge(); 
        else
          chargeVec = ptMyPDE_->getPDE_complexValuedCharge();

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
    
        if (typeOut==TRUE){
          std::cout<<"F(p)="<<F_hat<<"; \t";
          std::cout<<"\n ------------------------------- " <<std::endl;
               
        }

  } // end createF
  void piezoParamIdent::readInMeasurement(Vector<Double> & frequencies){

    ENTER_FCN("piezoParamIdent::readInMeasurement");

    std::cout<<"++ open and read file mess.dat ... " <<std::endl;    

    frequencies.Resize(nrMeasuredData);
    real.Resize(nrMeasuredData);
    imag.Resize(nrMeasuredData);

    Char* measuremets="mess.dat";
    mess = new std::ifstream(measuremets, std::basic_ios<char>::in);

    if (!mess){
      std::cerr << "\n File measuredData.dat does not exist!" << std::endl;
    }
   
   
    char mDataRow[256], helpChar[64];
    UInt i=0, j=0;
    Double newFreq, amplitude,phase;

    while(mess->getline(mDataRow, 265)){
      //  std::cout<<mDataRow<<std::endl;
      i=0;j=0;
      while (mDataRow[i]!='\t'){
        helpChar[j]=mDataRow[i];
        i++;j++;
      }// end while madataRow
      newFreq=atof(helpChar);
      i++;
      j=0;
      for(UInt k=0;k<i;k++) // Delete content of helpChar
        helpChar[k]='0';

      while (mDataRow[i]!='\t'){
        helpChar[j]=mDataRow[i];
        i++;j++;
      }// end while mdataRow
      amplitude=atof(helpChar);
      i++;
      j=0;
      for(UInt k=0;k<i;k++) // Delete content of helpChar
        helpChar[k]='0';
      while (mDataRow[i]!='\n'){
        helpChar[j]=mDataRow[i];
        i++;j++;
      }// end while mdataRow
      phase=atof(helpChar);

      for(UInt mInd=0;mInd<nrMeasuredData;mInd++){
        if(std::abs(freqs[mInd]-newFreq)<std::abs(freqs[mInd]-frequencies[mInd])){
          frequencies[mInd]=newFreq;
          real[mInd]=amplitude;
          imag[mInd]=phase;
        }
      }
    }// end while mess

    freqs.Resize(nrMeasuredData);
    freqs=frequencies;   
    mess->close();
    std::cout<<"++ open and read file finished ... " <<std::endl;    
  }// end readInMeasurements
      
} // end namespace
