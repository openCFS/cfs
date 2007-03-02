
#include "PDE/SinglePDE.hh" 
#include "Driver/stdSolveStep.hh"
#include "piezoParamIdent.hh"



namespace CoupledField
{

  // xxxxxxxxxxxxxxxxxxxxxxxxxx optimalExpDesign  - ITERATION xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

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
    nrMeasuredData=14;
    Vector<Double> freqs5;
    freqs5.Resize(14);
    freqs_.Resize(14);
   
  
   
    //initial guesses suggested by optExpDesignVarNr ...for thickness mode
    // Thicknes = 0.1mm
    if (true){
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
    if (false){
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
    if (false){
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
    if (false){
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
    if (false){
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

    freqs_=freqs5;

    //Writes out initial guesses of optimalFreqs

    *optimalFreqs<< 0 <<"  "<< 0 <<"  ";      
    
    for(UInt fr=0;fr<nrMeasuredData;fr++)
      *optimalFreqs<< freqs_[fr]<<"  ";
    *optimalFreqs<<std::endl;         
   

    projGradientFlags_.Resize(nrMeasuredData);
    projGradientFlags_.Init();
    parameterInitial_ = parameter_;
    residuumParIdentOld_ = residuumParIdent_=5.0e-3;
       

    std::cout<<"++ Optimum experiment desing with "<<nrMeasuredData <<
      " frequencies and "<<actNrParameter+actNrParameterC <<
      " parameter_ ..." << std::endl;
    Vector<Double> newFreqs;
    Vector<Double> newFreqsMech;
   
    readInMeasurement(newFreqs,newFreqsMech);
    //     newFreqs[0]=newFreqs[1];
    //     real_[0]=real_[1];
    //     imag_[0]=imag_[1];
    //     freqs_[0]=freqs_[1];
    std::cout<<"++ Frequencies taken from file mess.dat:"<<std::endl;
    std::cout<<newFreqs<<std::endl;

//     for(UInt fr=0;fr<nrMeasuredData;fr++)
//       *impedCurve<< freqs_[fr]<<"  ";
//     *impedCurve<<std::endl;

    calc_measuredCharge(freqs_, real_, imag_, y_hat_); // out of new measurements

    //  actNrParameterC=0;

//     if(directCoupling_==true)
//       ptMaterial_=ptPDE1_->getPDEMaterialData();   // Pointer to MaterialData
//     else
//       ptMaterial_=ptMyPDE_->getPDEMaterialData();   // Pointer to MaterialData

    // Frequency bounds

    // fr = 4000m/s / 2*thickness_

  
    Double rho;
    Vector<Complex> jacobi;
    Complex functional;

    //    rho = ptMaterial_->GetDensity();

    // fa_ = 1.0/(2*thickness_) * std::sqrt(1.3e+11/rho);

    fa_=antiResonanceFrequency_;
    fa_=1.1*fa_;

    fr_=resonanceFrequency_;
    fr_=0.9*fr_;
    // fa = 1/thickness_ * sqrt(c_33/rho)

    std::cout<<"++ Bounds for resonance and antiresonace frequency: "<<std::endl;
    std::cout<<"fr= "<< fr_ <<std::endl;
    std::cout<<"fa= "<< fa_ <<std::endl;

    for (UInt nrOptExpSteps=0; nrOptExpSteps<20; nrOptExpSteps++){ 
      UInt nrNuMethods=0;
      // That at least we perform one descent step
      normGradient_=2.0;
      normGradientOld_=1.0;

      descentMethod(functional);

      readInMeasurement(newFreqs,newFreqsMech);

      std::cout<<"++ Frequencies, just taken from mess.dat - file ... " <<std::endl;
      std::cout<<newFreqs<<std::endl;

      for(UInt fr=0;fr<nrMeasuredData;fr++)
        *impedCurve<< freqs_[fr]<<"  ";
      *impedCurve<<std::endl;

      calc_measuredCharge(freqs_, real_, imag_, y_hat_); // out of new measurements
      newtonCounter_=0;

      while (nrNuMethods<maxNumberNewtonLoops_ && 
             residuumParIdent_ > 0.05*residuumParIdentOld_){
        //while (nrNuMethods<maxNumberNewtonLoops_){
        ///        nuMethodsC2();
        nuMethods();
        nrNuMethods++;
        for (UInt i=0; i<parameter_.GetSize(); i++)
          *parFinal<<parameter_[i]<<", ";
        *parFinal<<"/"<<std::endl;

        for (UInt i=0; i<parameterC_.GetSize(); i++)
          *parFinal<<parameterC_[i]<<", ";
        *parFinal<<"/"<<std::endl;

        newtonCounter_++;
      }

      Vector<Double> saveFreqs=freqs_;

      Double startFreqTemp;
      startFreqTemp=startfreq_;
      Double freqincr=(stopfreq_-startfreq_)/50;
      freqs_.Resize(50);
      for(UInt i=0;i<50;i++){
        startFreqTemp+=freqincr;
        freqs_[i]=startFreqTemp;
      }
      calcImpedanceCurve();
      freqs_=saveFreqs;

      fa_=antiResonanceFrequency_;
      fa_=1.1*fa_;
      
      fr_=resonanceFrequency_;
      fr_=0.9*fr_;

      std::cout<<"++ New bounds for resonance and antiresonace frequency: "<<std::endl;
      std::cout<<"fr= "<< fr_ <<std::endl;
      std::cout<<"fa= "<< fa_ <<std::endl;


      residuumParIdentOld_ = residuumParIdent_;
      
    } // end nrOptExpSteps
    
    std::cout<<"\n END OF OPTIMAL EXPERIMENT DESIGN" <<std::endl;
     
  } // end optimalExpDesign

  void piezoParamIdent::optimalExpDesignDiffNumberFreqs(){

    ENTER_FCN("piezoParamIdent::optimalExpDesignDiffNumberFreqs");

    nrMeasuredData=1;
    Vector<Double> newFreqs;
    Vector<Double> newFreqsMech;

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

    newtonCounter_=0;

    for (UInt nrOptExpSteps=0; nrOptExpSteps<10; nrOptExpSteps++){ 

      nrMeasuredData++;
      
      freqs_.Resize(nrMeasuredData);

      if (nrMeasuredData==2)
        freqs_=freqs2;
      else  if (nrMeasuredData==3)
        freqs_=freqs3;
      if (nrMeasuredData==4)
        freqs_=freqs4;
      else  if (nrMeasuredData==5)
        freqs_=freqs5;
      if (nrMeasuredData==6)
        freqs_=freqs6;
      else  if (nrMeasuredData==7)
        freqs_=freqs7;
      else if(nrMeasuredData>7)
        break;

      readInMeasurement(newFreqs,newFreqsMech);
      
      calc_measuredCharge(freqs_, real_, imag_, y_hat_); // out of new measurements 

      UInt nrNuMethods=0;
      while (nrNuMethods<maxNumberNewtonLoops_){
        //    nuMethodsC2();
        nuMethods();
        nrNuMethods++;
        for (UInt i=0; i<parameter_.GetSize(); i++)
          *parFinal<<parameter_[i]<<", ";
        *parFinal<<"/"<<std::endl;
        
        for (UInt i=0; i<parameterC_.GetSize(); i++)
          *parFinal<<parameterC_[i]<<", ";
        *parFinal<<"/"<<std::endl;
        
        newtonCounter_++;
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
    //    lambda = 3.500e-03;

    // Fitting am radial, planar Mode
    //       lambda = 1.0e-14;
    //        Double lambda = 1.0e-5; // for ten parameter
    //    Double lambda = 1.0; // for nine parameter using M - confidence crit.
    //    Double lambda = 1.0e+7; // for nine parameter using M - confidence crit.
    //Double lambda = 40; // for three parameters using M - confidence crit.


    UInt maxNumberDescentIterations=5;
    gradientNormUpdate = std::abs(normGradient_ - normGradientOld_)/normGradientOld_;
    std::cout<<"gradientNormUpdate before descent iteration"<<std::endl;
    std::cout<<gradientNormUpdate<<std::endl;

    projGradientFlags_.Resize(nrMeasuredData);
    projGradientFlags_.Init();
   
    while((descIter<maxNumberDescentIterations && gradientNormUpdate>0.1)|| descIter<2){
      descIter++;
      Complex J, JOld;
      Vector<Double> freqsOld;
      freqsOld.Resize(nrMeasuredData);
      freqsOld=freqs_;
      createCovA(J,false);
      JOld=J;

      if (JOld.real()<=0.0){
        std::cout<< "! ! ! we should be in trouble here, since Re(J)<0 "<<std::endl;
        std::cout<<JOld.real()<<std::endl;
        //        getchar();
      }

      normGradientOld_=normGradient_;
      //overwrites normGradient!
      createGradient(gradient,1.0e-5);
      gradientNormUpdate = std::abs(normGradient_ - normGradientOld_)/normGradientOld_;

      //       std::cout<<"gradientNormOld"<<std::endl;
      //       std::cout<<normGradientOld_<<std::endl;
      //       std::cout<<"gradientNorm_"<<std::endl;
      //       std::cout<<normGradient_<<std::endl;

      //       std::cout<<"gradientNormUpdate"<<std::endl;
      //       std::cout<<gradientNormUpdate<<std::endl;

      lambda = 1.2*lambda;

      std::cout<<"gradient"<<std::endl;
      std::cout<<gradient<<std::endl;
    
      UInt innerCounter=0;
      //      while(innerCounter<1){
      while(JOld.real()<=J.real()){      
        for(UInt fr=0;fr<nrMeasuredData;fr++)
          freqs_[fr]=freqs_[fr]-lambda*gradient[fr].real();
    
        for(UInt fr=0;fr<nrMeasuredData;fr++)
          if (freqs_[fr]<=0.0){
            std::cout<<"! Frequency "<< fr << " is negative " <<std::endl;
            std::cout<<"! Frequency "<< fr << " will be reset to previous value " <<std::endl;
            freqs_[fr]=freqs_[fr]+lambda*gradient[fr].real();
          }
        
        std::cout<<"freqs after update by -lambda*GradJ:"<<std::endl;
        std::cout<<freqs_<<std::endl;
      
        createCovA(J,false);

        
        if (J.real()>JOld.real()){
          lambda=0.1*lambda;
          freqs_=freqsOld;
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
          
          if(true){
            if ((freqs_[fr]<=fa_)&&(freqs_[fr]>=fr_)){
              if(freqs_[fr]<(fr_+fa_)/2){
                freqs_[fr]=fr_;
                projGradientFlags_[fr]=freqs_[fr];
                std::cout<<"! Frequency "<< fr <<" to close to fr, fixed it at f= " << freqs_[fr]<< std::endl;
                if (fr>=1)
                  if (freqs_[fr]==freqs_[fr-1]){
                    freqs_[fr]=0.95*freqs_[fr];
                    std::cout<<"Frequency "<< fr <<" (fr) coincides with previos one( " << freqs_[fr-1] <<" ), set freq= " << freqs_[fr]<<std::endl;
                  }
                projGradientFlags_[fr]=freqs_[fr];
              }
              else{
                freqs_[fr]=fa_;
                projGradientFlags_[fr]=freqs_[fr];
                std::cout<<"! Frequency "<< fr <<" to close to fa, fixed it at f= " << freqs_[fr]<< std::endl;
                if (fr>=1)
                  if (freqs_[fr]==freqs_[fr-1]){
                    freqs_[fr]=1.05*fa_;
                    std::cout<<"Frequency "<< fr <<" (fr) coincides with previos one, set freq= " << freqs_[fr]<<std::endl;
                  }
              }
            } 
          } // end if false/true
          
          if (freqs_[fr]< startfreq_){
            std::cout<<"! There are no measurements for frequency ("<<fr<<") = " << freqs_[fr]  <<std::endl;
            std::cout<<"! Fixed frequency at " <<startfreq_<< std::endl;
            freqs_[fr]=startfreq_;
            if (fr>=1)
              if (freqs_[fr]==freqs_[fr-1]){
                freqs_[fr]=startfreq_+fr*0.1*startfreq_;
                std::cout<<"! Frequency "<< fr <<" to close to startfreq_ and coincides with previous one, fixed it at f= " << freqs_[fr]<< std::endl;
              }
            projGradientFlags_[fr]=freqs_[fr];
          }
          if(freqs_[fr]>stopfreq_){
            std::cout<<"! There are no measurements for frequency (" << fr<< ") = " <<freqs_[fr] <<std::endl;
            std::cout<<"! Fixed frequency at  " <<stopfreq_<<std::endl;
            freqs_[fr]=stopfreq_;
            if (fr>=1)
              if (freqs_[fr]==freqs_[fr-1]){
                freqs_[fr]=stopfreq_-fr*0.1*startfreq_;
                std::cout<<"! Frequency "<< fr <<" to close to stopfreq_ and coincides with previous one, fixed it at f= " << freqs_[fr]<< std::endl;
              }

            projGradientFlags_[fr]=freqs_[fr];
          }
            
            
        }// end in freqs_[fr]<=fa_ ...

        //   std::cout<<" frequencies before final check of coincidences:" <<std::endl;
        //       std::cout<<freqs_<<std::endl;


        // check if two frequencies coincide:
        // if they do so reduce them by 5 per cent

        bool flag=true;

        while(flag){
          for(UInt i=0;i<freqs_.GetSize();i++){
            for(UInt j=0;j<freqs_.GetSize();j++){
              if((freqs_[i]==freqs_[j]) && (i!=j))
                if (i<freqs_.GetSize()-1)
                  freqs_[j]=1.025*freqs_[j];
                else
                  freqs_[j]=0.975*freqs_[j];
              else
                flag=false;
            }
            // check if frequencies are out of range        
            while(freqs_[i]>stopfreq_)
              freqs_[i]=0.98*freqs_[i];
            while(freqs_[i]<startfreq_)
              freqs_[i]=1.02*freqs_[i];             
          }
        }
         
        createCovA(J,false);         
        // end while ....         
      }      

      // will write also confidence intervals if we found a smaller  J
      // if JNew>Jold after 10 inner loops, we might assume, that there is no improvement for 
      // the location of the frequencies and therefore write nothing into confidence intervals
      if (innerCounter<10)
        createCovA(J,true);
      std::cout<<" Value J(w) calculated with projected gradient "<< J <<std::endl;
      std::cout<<" frequencies at end of descent Method:" <<std::endl;
      std::cout<<freqs_<<std::endl;

      *optimalFreqs<< descIter <<"  "<< J.real() <<"  ";      

      for(UInt fr=0;fr<nrMeasuredData;fr++)
        *optimalFreqs<< freqs_[fr]<<"  ";
      *optimalFreqs<<std::endl;               
    } // while descIter<max Number


  }// descentMethod


  void piezoParamIdent::createCovA(Complex &J, bool writeOutCov){
    ENTER_FCN("piezoParamIdent::createCovA");

    //    Complex J;
    J=Complex(0.0,0.0);
    Double deltaLocal=1.01;

    Vector<Complex> jacobi;
    Vector<Complex> jacobiH;
    Complex F_hat__incr1,F_hat__incr2;

    Matrix<Complex> cov(actNrParameter+actNrParameterC, actNrParameter+actNrParameterC);
    Matrix<Complex> covTemp(actNrParameter+actNrParameterC, actNrParameter+actNrParameterC);

    for (UInt actFreq=0;actFreq<nrMeasuredData;actFreq++){

      Double stepWidth=0.0;
      stepWidth=0.5*std::abs((freqs_[std::min(actFreq+1,nrfreq_-1)]-freqs_[std::max((Integer)actFreq-1,0)]));
       
      createJacobian(jacobi, freqs_[actFreq]);
      jacobiH.Resize(actNrParameter+actNrParameterC);
      jacobiH.Init();
       
      for (UInt i=0;i<jacobi.GetSize();i++)
        jacobiH[i]=Complex(jacobi[i].real(),-jacobi[i].imag());
       
      for(UInt i=0;i<actNrParameter+actNrParameterC;i++)
        for(UInt j=0;j<actNrParameter+actNrParameterC;j++){
          //           ovTemp[i][j]=jacobiH[i]*jacobi[j]/(delta_*delta_*std::abs(y_hat_[actFreq])*std::abs(y_hat_[actFreq]));
          covTemp[i][j]= stepWidth*jacobiH[i]*jacobi[j]/
            ((deltaLocal)*(deltaLocal)*std::abs(y_hat_[actFreq])*std::abs(y_hat_[actFreq]));
           
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
    if (false){
      Double JJ=0.0;
      for (UInt i=0;i<actNrParameter+actNrParameterC;i++){
        if (std::sqrt(cov[i][i].real()) >= JJ)
          JJ=std::sqrt(cov[i][i].real());
        if (writeOutCov==true)
          *optimalFreqs<<cov[i][i].real()<<"  ";
      }
      J=Complex(JJ,0.0);
    }

    // A - Criterion for all parameters
    if(true){
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
      if(false){
        J=cov[7][7].real();
      }
      //optimality criterion just for four parameter, namely eps_11, e15, e13, c13!
      if(true){
        J=cov[7][7].real()+cov[2][2].real()+cov[4][4].real()+cov[0][0].real();
        J=J/Complex(4.0,0.0);
      }


         
      //      writeOutCov=true;
      if (writeOutCov==true){
   
        if (actNrParameter==3){      
          //          *confInterval<<parameter_[1]+parameterInitial_[1]*std::sqrt(cov[0][0].real()*0.115*0.115)<<"  ";
          *confInterval<<parameter_[1]+parameterInitial_[1]*(cov[0][0].real()*0.7)<<"  ";
          //        std::cout<<parameter_[1]+parameter[1]*std::sqrt(cov[1][1].real()*0.115*0.115)<<std::endl;
          *confInterval<<parameter_[1]<<"  ";
          //          *confInterval<<parameter_[1]-parameterInitial_[1]*std::sqrt(cov[0][0].real()*0.115*0.115)<<"  ";
          *confInterval<<parameter_[1]-parameterInitial_[1]*(cov[0][0].real()*0.7)<<"  ";
          //        std::cout<<parameter_[1]-1000.0*std::sqrt(cov[0][0].real()*parameter_[1]*0.115*0.115)<<std::endl;


          //           *confInterval<<parameter[7]+parameterInitial_[7]*std::sqrt(cov[1][1].real()*0.115*0.115)<<"  ";
          //           *confInterval<<parameter[7]<<"  ";
          //           *confInterval<<parameter[7]-parameterInitial_[7]*std::sqrt(cov[1][1].real()*0.115*0.115)<<"  ";

          *confInterval<<parameter_[7]+parameterInitial_[7]*(cov[1][1].real()*0.7)<<"  ";
          *confInterval<<parameter_[7]<<"  ";
          *confInterval<<parameter_[7]-parameterInitial_[7]*(cov[1][1].real()*0.7)<<"  ";


          //           *confInterval<<parameter_[9]+parameterInitial_[9]*std::sqrt(cov[2][2].real()*0.115*0.115)<<"  ";
          //           *confInterval<<parameter_[9]<<"  ";
          //           *confInterval<<parameter_[9]-parameterInitial_[9]*std::sqrt(cov[2][2].real()*0.115*0.115)<<"  ";

          *confInterval<<parameter_[9]+parameterInitial_[9]*(cov[2][2].real()*0.7)<<"  ";
          *confInterval<<parameter_[9]<<"  ";
          *confInterval<<parameter_[9]-parameterInitial_[9]*(cov[2][2].real()*0.7)<<"  ";

        }

        else if (actNrParameter==4){
          for (UInt ii=0;ii<=4;ii++){
            if(ii<=1){
              //               *confInterval<<parameter_[ii]+parameter[ii]*std::sqrt(cov[ii][ii].real()*0.115*0.115)<<"  ";
              //               *confInterval<<parameter_[ii]<<"  ";
              //               *confInterval<<parameter_[ii]-parameter[ii]*std::sqrt(cov[ii][ii].real()*0.115*0.115)<<"  ";

              *confInterval<<parameter_[ii]+parameterInitial_[ii]*(cov[ii][ii].real()*0.7)<<"  ";
              *confInterval<<parameter_[ii]<<"  ";
              *confInterval<<parameter_[ii]-parameterInitial_[ii]*(cov[ii][ii].real()*0.7)<<"  ";

            }
            else if(ii>2)
              {
                //                 *confInterval<<parameter_[ii]+parameterInitial_[ii]*std::sqrt(cov[ii-1][ii-1].real()*0.115*0.115)<<"  ";
                //                 *confInterval<<parameter_[ii]<<"  ";
                //                 *confInterval<<parameter_[ii]-parameterInitial_[ii]*std::sqrt(cov[ii-1][ii-1].real()*0.115*0.115)<<"  ";

                *confInterval<<parameter_[ii]+parameterInitial_[ii]*(cov[ii-1][ii-1].real()*0.7)<<"  ";
                *confInterval<<parameter_[ii]<<"  ";
                *confInterval<<parameter_[ii]-parameterInitial_[ii]*(cov[ii-1][ii-1].real()*0.7)<<"  ";

              }
          }
        }
        if (actNrParameter==10){
          for (UInt ii=0;ii<10;ii++){
            *confInterval<<parameter_[ii]+parameterInitial_[ii]*sqrt(cov[ii][ii].real()*2.5*2.5)<<"  ";
            *confInterval<<parameter_[ii]<<"  ";
            *confInterval<<parameter_[ii]-parameterInitial_[ii]*sqrt(cov[ii][ii].real()*2.5*2.5)<<"  ";
          }
        }

        if (actNrParameter==9){
          for (UInt ii=0;ii<10;ii++){
            if(ii<=1){
              *confInterval<<parameter_[ii]+parameterInitial_[ii]*sqrt(cov[ii][ii].real()*2.5*2.5)<<"  ";
              //              std::cout<<parameter[ii]+parameterInitial_[ii]*sqrt(cov[ii][ii].real()*2.5*2.5)<<"  ";
              *confInterval<<parameter_[ii]<<"  ";
              //              std::cout<<parameter[ii]<<std::endl;
              *confInterval<<parameter_[ii]-parameterInitial_[ii]*sqrt(cov[ii][ii].real()*2.5*2.5)<<"  ";
              //              std::cout<<parameter[ii]-parameterInitial_[ii]*sqrt(cov[ii][ii].real()*2.5*2.5)<<"  ";
              //              getchar();
            }
            else if (ii>2){
              *confInterval<<parameter_[ii]+parameterInitial_[ii]*sqrt(cov[ii-1][ii-1].real()*2.5*2.5)<<"  ";
              *confInterval<<parameter_[ii]<<"  ";
              *confInterval<<parameter_[ii]-parameterInitial_[ii]*sqrt(cov[ii-1][ii-1].real()*2.5*2.5)<<"  ";
            }

          }
        }

        if (actNrParameter==5){
          for (UInt ii=0;ii<5;ii++){
            *confInterval<<parameter_[ii]+parameterInitial_[ii]*std::sqrt(cov[ii][ii].real()*0.55*0.55)<<"  ";
            *confInterval<<parameter_[ii]<<"  ";
            *confInterval<<parameter_[ii]-parameterInitial_[ii]*std::sqrt(cov[ii][ii].real()*0.55*0.55)<<"  ";
          }
        }   

        if (actNrParameter==6){
          
          *confInterval<<parameter_[0]+parameterInitial_[0]*(cov[0][0].real()*0.03)<<"  ";
          *confInterval<<parameter_[0]<<"  ";
          *confInterval<<parameter_[0]-parameterInitial_[0]*(cov[0][0].real()*0.03)<<"  ";

          *confInterval<<parameter_[1]+parameterInitial_[1]*(cov[1][1].real()*1.0)<<"  ";
          *confInterval<<parameter_[1]<<"  ";
          *confInterval<<parameter_[1]-parameterInitial_[1]*(cov[1][1].real()*1.0)<<"  ";

          *confInterval<<parameter_[3]+parameterInitial_[3]*(cov[2][2].real()*0.03)<<"  ";
          *confInterval<<parameter_[3]<<"  ";
          *confInterval<<parameter_[3]-parameterInitial_[3]*(cov[2][2].real()*0.03)<<"  ";

          *confInterval<<parameter_[4]+parameterInitial_[4]*(cov[3][3].real()*0.03)<<"  ";
          *confInterval<<parameter_[4]<<"  ";
          *confInterval<<parameter_[4]-parameterInitial_[4]*(cov[3][3].real()*0.03)<<"  ";

          *confInterval<<parameter_[7]+parameterInitial_[7]*(cov[4][4].real()*0.75)<<"  ";
          *confInterval<<parameter_[7]<<"  ";
          *confInterval<<parameter_[7]-parameterInitial_[7]*(cov[4][4].real()*0.75)<<"  ";

          *confInterval<<parameter_[9]+parameterInitial_[9]*(cov[5][5].real()*0.75)<<"  ";
          *confInterval<<parameter_[9]<<"  ";
          *confInterval<<parameter_[9]-parameterInitial_[9]*(cov[5][5].real()*0.75)<<"  ";

        }
        *confInterval<<std::endl;
      }
    }
  

    
    // smallest eigenvalue criterion:
    if(false){

      Vector<Double> eig;
      eig.Resize(actNrParameter+actNrParameterC);
      eig.Init();
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
    freqsOld=freqs_;

    grad.Resize(nrMeasuredData);
    grad.Init();
    Complex J1,J2;
    createCovA(J1,false);

    for (UInt j=0;j<nrMeasuredData;j++){
      if (projGradientFlags_[j]==0.0){
        freqs_[j]=freqs_[j]+dOmega*freqs_[j];
        createCovA(J2,false);
        for(UInt i=0; i<actNrParameter+actNrParameterC;i++)
          grad[j]= (J2-J1)* voltage_*voltage_ / (dOmega*freqs_[j]) ;
        //grad[j]= (J2-J1) / (dOmega*freqs_[j]) ;
        freqs_[j]=freqsOld[j];
      }

    }

    for (UInt i=0; i<nrMeasuredData;i++)
      normGradient_+=grad[i].real()*grad[i].real();
    normGradient_=std::sqrt(normGradient_);


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
    e.Init();
    x.Resize(actNrParameter+actNrParameterC);
    x.Init();
    Matrix<Complex> inverse;
    inverse.Resize(actNrParameter+actNrParameterC, actNrParameter+actNrParameterC);
    inverse.Init();

    for (UInt ind=0;ind<data.GetSizeRow();ind++){
      e[ind]=1.0;

      dataTemp.DirectSolve(x,e);
      dataTemp=data;
      for (UInt indC=0;indC<data.GetSizeCol();indC++)
        inverse[indC][ind]=x[indC];
      x.Resize(actNrParameter+actNrParameterC);
      x.Init();
      e[ind]=0.0;
    }
    data=inverse;

  }
  void piezoParamIdent::createJacobian(Vector<Complex> & jacobi, Double omega){

    ENTER_FCN("piezoParamIdent::createJacobian");
    //    std::cout<<"create Jacobian .."<<std::endl;

    Vector<Double>parIncr1, parIncr2;
    parIncr1.Resize(nrParameter_);
    parIncr1.Init();
    parIncr2.Resize(nrParameter_);
    parIncr2.Init();
    Vector<Double>parIncrC1, parIncrC2;
    parIncrC1.Resize(nrParameter_);
    parIncrC1.Init();
    parIncrC2.Resize(nrParameter_);
    parIncrC2.Init();
    parIncr1=parameter_;
    parIncr2=parameter_;
    parIncrC1=parameterC_;
    parIncrC2=parameterC_;
    Complex F_hat__incr1,F_hat__incr2;
    Complex F_hat__incrC1,F_hat__incrC2;
    Integer parInd=0;
    jacobi.Resize(actNrParameter+actNrParameterC);
    jacobi.Init();


//     if(directCoupling_==true)
//       ptMaterial_=ptPDE1_->getPDEMaterialData();   // Pointer to MaterialData
//     else
//       ptMaterial_=ptMyPDE_->getPDEMaterialData();   // Pointer to MaterialData
    
    for (UInt ind_param=0;ind_param<nrParameter_;ind_param++){ 
      if (whichParameterToUpdate_[ind_param]==1){
        
        parIncr1[ind_param]=1.001*parameter_[ind_param];
        updateMaterialData(parIncr1);
        createFVec(F_hat__incr1,false,omega);
        
        parIncr2[ind_param]=0.999*parameter_[ind_param];  
        updateMaterialData(parIncr2);
        createFVec(F_hat__incr2,false,omega);

        jacobi[parInd]=0.5*(F_hat__incr1-F_hat__incr2)/
          ((parIncr1[ind_param]-parameter_[ind_param])*scaling_[ind_param]);
        
        parIncr1[ind_param]=parameter_[ind_param];
        parIncr2[ind_param]=parameter_[ind_param];
        parInd++;

      }
    }
    parInd=0;
    for (UInt ind_param=0;ind_param<nrParameter_;ind_param++){ 
      if (whichParameterToUpdateC_[ind_param]==1){
        
        parIncrC1[ind_param]=1.000001*parameterC_[ind_param];
        //      std::cout<<parameter_incr<<std::endl
        updateMaterialData(parIncrC1);
        createFVec(F_hat__incrC1,false,omega);
        
        parIncrC2[ind_param]=0.999999*parameterC_[ind_param];  

        updateMaterialData(parIncrC2);
        createFVec(F_hat__incrC2,false,omega);

        jacobi[actNrParameter+parInd]=0.5*(F_hat__incrC1-F_hat__incrC2)/
          ((parIncrC1[ind_param]-parIncrC2[ind_param])*scalingC_[ind_param]);
        
        parIncrC1[ind_param]=parameterC_[ind_param];
        parIncrC2[ind_param]=parameterC_[ind_param];
        parInd++;

      }
    }

  } // end create Jacobian

  void piezoParamIdent::createFVec(Complex & F_hat_, bool typeOut,
                                   Double frequency){
    ENTER_FCN("PiezoParamIdent:createFVec");
    //    std::cout<<"createFVec ...."<<std::endl;
      

    if(directCoupling_==true){
      //      ptMaterial_=ptPDE1_->getPDEMaterialData();   // Pointer to MaterialData
      ptAssemble_ = ptPDE1_->getPDE_assemble();
      ptAlgsys_ = ptPDE1_->getPDE_algsys();    

    }
    else{
      //      ptMaterial_=ptMyPDE_->getPDEMaterialData();   // Pointer to MaterialData
      ptAssemble_ = ptMyPDE_->getPDE_assemble();
      ptAlgsys_ = ptMyPDE_->getPDE_algsys();    
      ptAlgsys_->InitMatrix();
      ptAlgsys_->InitRHS();
    }

    Error( "SetReassemble not present anymore -> contact Andreas!",
           __FILE__, __LINE__ );

    //ptAssemble_->SetReassemble();
   

      

    ////////////////////////////////////////////////////////
    //                   SOLVES PDE                      //
    ///////////////////////////////////////////////////////  
      
      ptPDE_->GetSolveStep()->SetActFreq(frequency); 
      ptPDE_->GetSolveStep()->SetActStep(0);       
      ptPDE_->GetSolveStep()->PreStepHarmonic(); 
      ptPDE_->GetSolveStep()->SolveStepHarmonic();
      ptPDE_->GetSolveStep()->PostStepHarmonic();

      //////////////////////////////////////////////////////////
      //Retrieves & stores Solution for further calculations  //
      /////////////////////////////////////////////////////////

        Vector<Complex> chargeVec;
        if(directCoupling_==true)      
          chargeVec = ptPDE1_->getPDE_complexValuedCharge(); 
        else
          chargeVec = ptMyPDE_->getPDE_complexValuedCharge();

        Complex charge=Complex(0.0,0.0);
         
        for (UInt i=0;i<chargeVec.GetSize();i++){
          charge+=chargeVec[i];
        }

        Integer fstep=0;
        Double x=real_[fstep]*cos(PI/180*imag_[fstep]);
        Double y=real_[fstep]*sin(PI/180*imag_[fstep]);
        Complex Z=Complex(x,y);

        // Logarithmic value of F
        F_hat_=(sign_*charge*Z)/std::log(Z); // without minus --- classical way ...
    
        if (typeOut==true){
          std::cout<<"F(p)="<<F_hat_<<"; \t";
          std::cout<<"\n ------------------------------- " <<std::endl;
               
        }

  } // end createF
 
      
} // end namespace






































































