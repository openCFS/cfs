// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "PDE/SinglePDE.hh" 
#include "piezoParamIdent.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include <iostream>

namespace CoupledField
{

  // xxxxxxxxxxxxxxxxxxxxxxxxxx optimalExDesignVarNrFreqs xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  void piezoParamIdent::optimalExpDesignVarNrFreqs(){


    bool firstTime=true;
    parameterInitial_ = parameter_;
    //    std::cout<<"parameter_"<<std::endl;
    //     std::cout<<parameter_<<std::endl;

    residuumParIdentOld_ = residuumParIdent_=2000;

    for(UInt optExpVarFreqCout=0; optExpVarFreqCout<35; optExpVarFreqCout++){

      myParam_->Get( "numFreq", nrfreq_ );

      Vector<Double> newFreqs;
      Vector<Double> newFreqsMech;

      //      if (optExpVarFreqCout==0){

        Double hOmega = (stopfreq_-startfreq_)/nrfreq_;
        freqs_.Resize(nrfreq_);
        freqs_.Init();
        nrMeasuredData=nrfreq_;

        // writes out confidence intervals for maximal number of frequencies
        if (false){
          Double dOmega = 0.0001;
          rhos.Resize(nrfreq_);      
          rhos.Init();
          for(UInt i=0;i<nrfreq_;i++)
            rhos[i]=1.0;
          for (UInt actFreq=0;actFreq<nrfreq_;actFreq++)
            freqs_[actFreq] = startfreq_+actFreq*hOmega;
          std::cout<<"message 0 " <<std::endl;
          readInMeasurement(newFreqs,newFreqsMech);  // get frequencies
          std::cout<<"message 1 " <<std::endl;
          calc_measuredCharge(freqs_, real_, imag_, y_hat_); // indispensable for calculation of conf intervals
          std::cout<<"message 2 " <<std::endl;
          writeOutConfInterval(); 
          std::cout<<" Wrote out confidence interval with " << nrfreq_ << " frequencies " <<std::endl;
          std::cout<<" get any key to continue ... " <<std::endl;
          getchar();
        }


        // to compare write out confidence interval for twelve arbitrarily,
        // but equidistantly chosen frequencies        
        if (false){
          std::cout<<"message 00 " <<std::endl;
          UInt nrfreqOld = nrfreq_;
          hOmega = (stopfreq_-startfreq_)/12;
          freqs_.Resize(nrfreq_);
          freqs_.Init();
          nrMeasuredData=12;
          nrfreq_=12;
          std::cout<<"message 0 " <<std::endl;
          for (UInt actFreq=0;actFreq<12;actFreq++)
            freqs_[actFreq] = startfreq_+actFreq*hOmega;
          readInMeasurement(newFreqs,newFreqsMech);
          std::cout<<"message 1 " <<std::endl;
          calc_measuredCharge(freqs_, real_, imag_, y_hat_); // out of new measurements
          std::cout<<"message 2 " <<std::endl;
          writeOutConfInterval();
          std::cout<<"message 3 " <<std::endl;

          // reset everything to avoid erroneous calculations in the following
          nrfreq_=nrfreqOld;
          hOmega = (stopfreq_-startfreq_)/nrfreq_;
          freqs_.Resize(nrfreq_);
          freqs_.Init();
          nrMeasuredData=nrfreq_;
        }

        // the normal procedure
        for (UInt actFreq=0;actFreq<nrfreq_;actFreq++)
          freqs_[actFreq] = startfreq_+actFreq*hOmega;
        
        std::cout<<"Until here " <<std::endl;
        readInMeasurement(newFreqs,newFreqsMech);
        std::cout<<"newFreqs:"<<std::endl;
        //        std::cout<<newFreqs<<std::endl;
        calc_measuredCharge(freqs_, real_, imag_, y_hat_); // out of new measurements
        
        if (firstTime==true){
          std::cout<<"++ Calculate weighting function rho(omega) for the following frequencies"<<std::endl;
          std::cout<<newFreqs<<std::endl;
          firstTime=false;
        }
        
        UInt maxNrBreakpoints = nrfreq_;
        Complex J;
        
        Double dOmega = 0.0001;
        rhos.Resize(nrfreq_);
        rhos.Init();
        
        for(UInt i=0;i<nrfreq_;i++)
          rhos[i]=0.01;

        writeOutConfInterval();
        
        std::cout<<"++ Minimize the functional J(omega) ... "<<std::endl;

        descentMethodRho(J);

        // with all frequencies:
        globalIndexSet_.Resize(nrfreq_);
        for (UInt k=0;k<nrfreq_;k++)
          globalIndexSet_[k]=1;

        //        writeOutConfInterval();
        std::cout<<"++ Calculated rhos: ... "<<std::endl;        
        std::cout<<rhos<<std::endl;
        
        // choose 8 highest values ...
        Vector<Double> highestFreqs;
        Vector<Double> highestRhos;
        UInt maxIndex;
        Double sumRhos=0.0;
        highestFreqs.Resize(nrfreq_);
        highestFreqs.Init();
        highestRhos.Resize(nrfreq_);
        highestRhos.Init();
        
        omegaDiffVec_.Resize(nrfreq_);
        omegaDiffVec_.Init();
        
        for (UInt i=0;i<nrfreq_;i++){
          maxIndex=0;
          highestFreqs[i]=rhos[0];
          
          for (UInt j=0;j<rhos.GetSize();j++){
            if (highestRhos[i]<=rhos[j]){ // hier muss " <= " hin!!!
              highestFreqs[i]=freqs_[j];
              highestRhos[i]=rhos[j];
              maxIndex=j;
            }
          }
          rhos[maxIndex]=0.0;
          
          Double omegaDiff=nrfreq_*hOmega;
          
          // Get Neighbor of \omega_i
          for(UInt nI=0;nI<highestFreqs.GetSize();nI++){
            if (std::abs((highestFreqs[i]-highestFreqs[nI]))<omegaDiff&&(nI!=i))
              omegaDiff=std::abs(highestFreqs[i]-highestFreqs[nI]);
            
          }
          //  std::cout<<" omegaDiff = " <<omegaDiff<<std::endl;
          omegaDiffVec_[i]=omegaDiff;
          sumRhos+=highestRhos[i]*omegaDiff;
          std::cout<<" - sumRhos = " << sumRhos<<std::endl;
          
          // fits for thickness 0.1mm and thickness mode
          //  if ((sumRhos>4.9e+06 && i>10) || i>16){

          // fits for thickness 0.05mm and thickness mode
          if ((sumRhos>8.4e+06 && i>10) || i>16){

            //        if (sumRhos>6*hOmega && i>3){
            std::cout<<"sumRhos = "<<sumRhos <<" ecxeeds 33.5*hOmega = " 
                     << 33.5*hOmega << "; hOmega = " << hOmega<< std::endl;
            //          getchar();
            break;
          }
        } // end for i ...
            
        globalIndexSet_.Resize(nrfreq_);
        globalIndexSet_.Init();
        
        Integer howManyFreqs=0;
        for (UInt i=0;i<nrfreq_;i++)
          if (highestFreqs[i]!=0.0){
            howManyFreqs++;
          }
                    
        std::cout<<"howManyFreqs = "<< howManyFreqs<<std::endl;
        nrfreq_=howManyFreqs;
        *nrOfFreqs<< howManyFreqs<<std::endl;
        
        // sort highestFreqs
        Double dv,dr,dO;
        Integer k, ii;
        
        for (k=0; k<howManyFreqs-1; ++k)
          for (ii=k+1; ii<howManyFreqs; ++ii)
            if (highestFreqs[ii] < highestFreqs[k])
              {
                dv = highestFreqs[k];
                dr = highestRhos[k];
                dO = omegaDiffVec_[k];
                highestFreqs[k] = highestFreqs[ii];
                highestRhos[k] = highestRhos[ii];
                omegaDiffVec_[k]=omegaDiffVec_[ii];
                highestFreqs[ii] = dv;
                highestRhos[ii] = dr;
                omegaDiffVec_[ii]=dO;
              }

        for (UInt i=0; i<highestFreqs.GetSize(); i++)
          for (UInt j=0;j<howManyFreqs;j++)
            if (highestFreqs[j]==freqs_[i])
              globalIndexSet_[i]=1;

        std::cout<<" ++ finished to solve optimal experiment design " <<std::endl;
        std::cout<<"    the following results might be used as initial guesses "<<std::endl;
        std::cout<<"    for optimal experiment design with fixed number of frequencies:" <<std::endl;
        
        std::cout<<"highestRhos"<<std::endl;
        std::cout<<highestRhos<<std::endl;
        
        std::cout<<"highestFreqs"<<std::endl;
        std::cout<<highestFreqs<<std::endl;

        nrMeasuredData=howManyFreqs;
        //      rhos.Resize(howManyFreqs);
        freqs_.Resize(howManyFreqs);
        freqs_.Init();
        real_.Resize(howManyFreqs);
        real_.Init();
        imag_.Resize(howManyFreqs);
        imag_.Init();
        

        for (UInt j=0;j<nrfreq_;j++){
          freqs_[j]=highestFreqs[j];
          //rhos[j]=highestRhos[j];
        }

        // Project frequencies into the feasible set ...
        fa_=antiResonanceFrequency_;
        fa_=1.05*fa_;
        
        fr_=resonanceFrequency_;
        fr_=0.95*fr_;
        
        for(UInt fr=0;fr<freqs_.GetSize();fr++){

          if(true){
            if ((freqs_[fr]<=fa_)&&(freqs_[fr]>=fr_)){
              if(freqs_[fr]<(fr_+fa_)/2){
                freqs_[fr]=fr_;
                //                projGradientFlags_[fr]=freqs_[fr];
                std::cout<<"! Frequency "<< fr <<" to close to fr, fixed it at f= " << freqs_[fr]<< std::endl;
                if (fr>=1)
                  if (freqs_[fr]==freqs_[fr-1]){
                    freqs_[fr]=0.95*freqs_[fr];
                    std::cout<<"Frequency "<< fr <<" (fr) coincides with previos one( " << freqs_[fr-1] <<" ), set freq= " << freqs_[fr]<<std::endl;
                  }
                //                projGradientFlags_[fr]=freqs_[fr];
              }
              else{
                freqs_[fr]=fa_;
                //              projGradientFlags_[fr]=freqs_[fr];
                std::cout<<"! Frequency "<< fr <<" to close to fa, fixed it at f= " << freqs_[fr]<< std::endl;
                if (fr>=1)
                  if (freqs_[fr]==freqs_[fr-1]){
                    freqs_[fr]=1.05*fa_;
                    std::cout<<"Frequency "<< fr <<" (fr) coincides with previos one, set freq= " << freqs_[fr]<<std::endl;
                  }
              }
            } 
          }
          
          if (freqs_[fr]< startfreq_){
            std::cout<<"! There are no measurements for frequency ("<<fr<<") = " << freqs_[fr]  <<std::endl;
            std::cout<<"! Fixed frequency at " <<startfreq_<< std::endl;
            freqs_[fr]=startfreq_;
            if (fr>=1)
              if (freqs_[fr]==freqs_[fr-1]){
                freqs_[fr]=startfreq_+fr*0.1*startfreq_;
                std::cout<<"! Frequency "<< fr <<" to close to startfreq_ and coincides with previous one, fixed it at f= " << freqs_[fr]<< std::endl;
              }
            //          projGradientFlags_[fr]=freqs_[fr];
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

            //          projGradientFlags_[fr]=freqs_[fr];
          }
            
            
        }// end in freqs_[fr]<=fa_ ...

        std::cout<<" frequencies before final check of coincidences:" <<std::endl;
        std::cout<<freqs_<<std::endl;

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
             
        std::cout<<" frequencies at end of descent Method:" <<std::endl;
        std::cout<<freqs_<<std::endl;
        
        // end while ....

        Integer nrfreqAll;
        myParam_->Get( "numFreq", nrfreqAll );
                
        for(UInt i=0; i<nrfreqAll;i++)
          if(i<freqs_.GetSize())
            *optimalFreqs<<freqs_[i]<<"  ";
          else
            *optimalFreqs<<" 0.0 ";
        
        *optimalFreqs<<std::endl;

        std::cout<<"Press any key to continue with paramter identificaton"<<std::endl;        
        //        getchar();
        
        projGradientFlags_.Resize(nrMeasuredData);
        projGradientFlags_.Init();
        parameterInitial_ = parameter_;
        //        residuumParIdentOld_ = residuumParIdent_=0.5;
       
        std::cout<<"optimum experiment desing with "<<nrMeasuredData 
                 <<" frequencies and "<<actNrParameter+actNrParameterC
                 << " parameter! " << std::endl;
        
        readInMeasurement(newFreqs,newFreqsMech);
        std::cout<<newFreqs<<std::endl;
        
//         for(UInt fr=0;fr<nrMeasuredData;fr++)
//           *impedCurve<< freqs_[fr]<<"  ";
//         *impedCurve<<std::endl;
        
        calc_measuredCharge(freqs_, real_, imag_, y_hat_); // out of new measurements
        
        //  actNrParameterC=0;
        
        Double rho;
        Vector<Complex> jacobi;
        Complex functional;
        
        std::cout<<"++ Bounds for resonance and antiresonace frequency: "<<std::endl;
        std::cout<<"fr= "<< fr_ <<std::endl;
        std::cout<<"fa= "<< fa_ <<std::endl;
                
        //      for (UInt nrOptExpSteps=0; nrOptExpSteps<20; nrOptExpSteps++){ 
        //           UInt nrNuMethods=0;
        // That at least we perform one descent step
        normGradient_=2.0;
        normGradientOld_=1.0;
                                  
        newtonCounter_=0;
        UInt nrNuMethods=0;
          
        //        while (nrNuMethods<maxNumberNewtonLoops_ && residuumParIdent_ > 0.05*residuumParIdentOld_){
        while (nrNuMethods<maxNumberNewtonLoops_){
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
        residuumParIdentOld_ = residuumParIdent_;
        
        //} // end nrOptExpSteps
    }
    
    std::cout<<"\n END OF OPTIMAL EXPERIMENT DESIGN" <<std::endl;
  }

  void piezoParamIdent::writeOutConfInterval(){

    std::cout<<"++ Write out confidence intervals into file confInterval.dat " <<std::endl;

    //     std::cout<<"globalIndexSet:"<<std::endl;
    //     std::cout<<globalIndexSet<<std::endl;

    std::cout.precision( 8 );


    Double deltaLocal=1.01;
    Double  hOmega = (stopfreq_-startfreq_)/nrfreq_;

    Vector<Complex> jacobi;
    Vector<Complex> jacobiH;
  
    Matrix<Complex>covTempWithOutWeight(actNrParameter+actNrParameterC, actNrParameter+actNrParameterC);
    Matrix<Complex>covWithOutWeight(actNrParameter+actNrParameterC, actNrParameter+actNrParameterC);

    for (UInt actFreq=0;actFreq<nrfreq_;actFreq++){

      // std::cout<<freqs_<<std::endl;
      createJacobian(jacobi, freqs_[actFreq]);
      jacobiH.Resize(actNrParameter+actNrParameterC);
      jacobiH.Init();
      if(actFreq==0){
        std::cout<<"jacobi:"<<std::endl;
        std::cout<<jacobi<<std::endl;
      }

      for (UInt i=0;i<jacobi.GetSize();i++)
        jacobiH[i]=Complex(jacobi[i].real(),-jacobi[i].imag());

      globalJacobi_=jacobi;
      globalJacobiH_=jacobiH;

      //    for(UInt i=0; i<globalIndexSet.GetSize(); i++)
      // if (globalIndexSet[actFreq]==1){
      //     if (freqs_[actFreq]<=fr_||freqs_[actFreq]>=fa_){

      Double stepWidth=0.0;
      stepWidth=0.5*(freqs_[std::min(actFreq+1,nrfreq_-1)]-freqs_[std::max((Integer)actFreq-1,0)]);
      //    std::cout<<"right Freq"<<std::endl;
      //           std::cout<<freqs_[std::min(actFreq+1,nrfreq_-1)]<<std::endl;
      //           std::cout<<"left Freq"<<std::endl;
      //           std::cout<<freqs_[std::max((Integer)actFreq-1,0)]<<std::endl;
      //           //          stepWidth/=1.0e+05;
      //           //   std::cout<<0.5*Double(nrfreq_)<<std::endl;

      //           std::cout<<"stepWidth: = "<<stepWidth<<std::endl;

      //     std::cout<<" _________________"<<std::endl;
      for(UInt i=0;i<actNrParameter+actNrParameterC;i++)
        for(UInt j=0;j<actNrParameter+actNrParameterC;j++){
              
          //                if(actFreq==0||actFreq==nrfreq_)
          //                  covTempWithOutWeight[i][j]=0.5*stepWidth*globalJacobiH_[i]*globalJacobi_[j]/
          //                    (freqs_[actFreq]*(deltaLocal)*(deltaLocal)*std::abs(y_hat_[actFreq])*std::abs(y_hat_[actFreq]));
          //                else
          covTempWithOutWeight[i][j]= rhos[actFreq]*stepWidth*globalJacobiH_[i]*globalJacobi_[j]/
            ((deltaLocal)*(deltaLocal)*std::abs(y_hat_[actFreq])*std::abs(y_hat_[actFreq]));

          //   covTempWithOutWeight[i][j]=globalJacobiH_[i]*globalJacobi_[j]/
          //                 ((deltaLocal)*(deltaLocal)*std::abs(y_hat_[actFreq])*std::abs(y_hat_[actFreq]));


          //    if (i==j)
          //    covWithOutWeight[i][i]+=covWithOutWeight[i][i];
        }
      // }
      covWithOutWeight+=covTempWithOutWeight;

      //            for(UInt i=0;i<actNrParameter+actNrParameterC;i++)
      //              for(UInt j=0;j<actNrParameter+actNrParameterC;j++)
      //                covWithOutWeight[i][j]=covWithOutWeight[i][j]/Complex(Double(nrfreq_),0);
      //       /
      //        std::cout<<"cov"<<std::endl;
      //        std::cout<<cov<<std::endl;
      // } // end if global IndexSet .
    } // end for actFreq ..

    
    globalCov_=covWithOutWeight;
    //     std::cout<<"globalCov_ = cov wihtout Weight :"<<std::endl;
    //     std::cout<<covWithOutWeight<<std::endl;

    invert(covWithOutWeight);

    Matrix<Complex> testId;
    testId=globalCov_*covWithOutWeight;
    //    std::cout<<testId<<std::endl;
  
    globalCov_=covWithOutWeight;

    //     std::cout<<"globalCov_ = cov wihtout Weight inverted:"<<std::endl;
    //     std::cout<<globalCov_<<std::endl;
   
    //    getchar();

    //   for(UInt i=0;i<actNrParameter+actNrParameterC;i++)
    //       for(UInt j=0;j<actNrParameter+actNrParameterC;j++)
    //         covWithOutWeight[i][j]=1.0e-9*covWithOutWeight[i][j];

    
    if (actNrParameter==9){
      std::cout<<"++ Confidence intervals will be written into confInterval.dat " <<std::endl;
      for (UInt ii=0;ii<10;ii++){
        if(ii<=1){
          std::cout<<parameter_[ii]+parameterInitial_[ii]*sqrt(globalCov_[ii][ii].real()*2.5*2.5)<<"  ";
          *confInterval<<parameter_[ii]+parameterInitial_[ii]*sqrt(globalCov_[ii][ii].real()*2.5*2.5)<<"  ";
          *confInterval<<parameter_[ii]<<"  ";
          std::cout<<parameter_[ii]-parameterInitial_[ii]*sqrt(globalCov_[ii][ii].real()*2.5*2.5)<<"  ";
          *confInterval<<parameter_[ii]-parameterInitial_[ii]*sqrt(globalCov_[ii][ii].real()*2.5*2.5)<<"  ";
          // getchar();
        }
        else if (ii>2){
          *confInterval<<parameter_[ii]+parameterInitial_[ii]*sqrt(globalCov_[ii-1][ii-1].real()*2.5*2.5)<<"  ";
          *confInterval<<parameter_[ii]<<"  ";
          *confInterval<<parameter_[ii]-parameterInitial_[ii]*sqrt(globalCov_[ii-1][ii-1].real()*2.5*2.5)<<"  ";
        }
        
      }
      *confInterval<<std::endl;
    }
    std::cout.precision(6);
    
  } // end writeOutConfInterval

  void piezoParamIdent::createJRho(Complex &J, bool writeOutCov){

    Double  hOmega = (stopfreq_-startfreq_)/nrfreq_;
    delta_=0.00;
    Double deltaLocal=0.01;

    //    Complex J;
    J=Complex(0.0,0.0);

    Vector<Complex> jacobi;
    Vector<Complex> jacobiH;

    Double penaltyFactor1=1.0;
    Double penaltyFactor2=1.0;
    // no penalization

    fr_=4000/(2*thickness_);
    // minus 5 percent
    fr_=0.94*fr_;
    // fa = 1/thickness_ * sqrt(c_33/rho)
   
    Double rho;
    //    rho = ptMaterial_->GetDensity();

    fa_ = 1.0/(2*thickness_) * std::sqrt(1.3e+11/rho);
    fa_=1.2*fa_;

    // for thicker probes:
    //     fa_=1.2e+06;
    //     fr_=0.85e+06;

    //     std::cout<<"++ Bounds for resonance and antiresonace frequency"<<std::endl;
    //     std::cout<<"fr:"<<std::endl;
    //     std::cout<<fr_<<std::endl;

    //     std::cout<<"fa:"<<std::endl;
    //     std::cout<<fa_<<std::endl;

    //    ptMaterial_=ptMyPDE_->getPDEMaterialData();   // Pointer to MaterialData
    Complex F_hat__incr1,F_hat__incr2;

    Matrix<Complex> cov(actNrParameter+actNrParameterC, actNrParameter+actNrParameterC);
    Matrix<Complex> covTemp(actNrParameter+actNrParameterC, actNrParameter+actNrParameterC);
    Matrix<Complex>covTempWithOutWeight(actNrParameter+actNrParameterC, actNrParameter+actNrParameterC);
    Matrix<Complex>covWithOutWeight(actNrParameter+actNrParameterC, actNrParameter+actNrParameterC);

    for (UInt actFreq=0;actFreq<nrfreq_;actFreq++){
    
      createJacobian(jacobi, freqs_[actFreq]);

      jacobiH.Resize(actNrParameter+actNrParameterC);
      jacobiH.Init();
       
      //       std::cout<<jacobi<<std::endl;
      for (UInt i=0;i<jacobi.GetSize();i++)
        jacobiH[i]=Complex(jacobi[i].real(),-jacobi[i].imag());

      globalJacobi_=jacobi;
      globalJacobiH_=jacobiH;

      //    if (freqs_[actFreq]<=fr_||freqs_[actFreq]>=fa_){   

      for(UInt i=0;i<actNrParameter+actNrParameterC;i++)
        for(UInt j=0;j<actNrParameter+actNrParameterC;j++){

          covTempWithOutWeight[i][j]=jacobiH[i]*jacobi[j]/(std::abs(y_hat_[actFreq])*std::abs(y_hat_[actFreq]));

          //           covTemp[i][j]=jacobiH[i]*jacobi[j]/(delta*delta*std::abs(y_hat_[actFreq])*std::abs(y_hat_[actFreq]));
          if(actFreq==0||actFreq==nrfreq_)
            covTemp[i][j]=rhos[actFreq]*hOmega*0.5*jacobiH[i]*jacobi[j]/
              ((deltaLocal)*(deltaLocal)*std::abs(y_hat_[actFreq])*std::abs(y_hat_[actFreq]));

          //              covTemp[i][j]=rhos[actFreq]*hOmega*0.5*jacobiH[i]*jacobi[j]/
          //  (freqs_[actFreq]*(deltaLocal)*(deltaLocal)*std::abs(y_hat_[actFreq])*std::abs(y_hat_[actFreq]));

          //                 covTemp[i][j]=rhos[actFreq]*hOmega*0.5*jacobiH[i]*jacobi[j]/
          //                   (0.8*freqs_[actFreq]*(1+delta_)*(1+delta_)*std::abs(y_hat_[actFreq])*std::abs(y_hat_[actFreq]));

          else 
            covTemp[i][j]=rhos[actFreq]*hOmega*jacobiH[i]*jacobi[j]/
              ((deltaLocal)*(deltaLocal)*std::abs(y_hat_[actFreq])*std::abs(y_hat_[actFreq]));

          //              covTemp[i][j]=rhos[actFreq]*hOmega*jacobiH[i]*jacobi[j]/
          //  (freqs_[actFreq]*(deltaLocal)*(deltaLocal)*std::abs(y_hat_[actFreq])*std::abs(y_hat_[actFreq]));

          //                covTemp[i][j]=rhos[actFreq]*hOmega*jacobiH[i]*jacobi[j]/
          //                  (freqs_[actFreq]*(1+delta_)*(1+delta_)*std::abs(y_hat_[actFreq])*std::abs(y_hat_[actFreq]));

          //  if (i==j){
          //               covTemp[i][i]+=covTemp[i][i];// *1.0e-10;
          //               covWithOutWeight[i][i]+=covWithOutWeight[i][i];
          //             }
        }
      //}
      cov+=covTemp;
      // covWithOutWeight+=covTempWithOutWeight;
      //        std::cout<<"cov"<<std::endl;
      //        std::cout<<cov<<std::endl;
    } // end for actFreq ...

    
    //     invertWithLapack(cov);

    Matrix<Complex> data;
    data.Resize(actNrParameter+actNrParameterC,actNrParameter+actNrParameterC);
    data=cov; 
   
    //     std::cout<<"cov:"<<std::endl; 
    //     std::cout<<cov<<std::endl;
    //     getchar();
    invert(cov);

    Matrix<Complex> Identity;
    Identity = data * cov;
    for (UInt k=0;k<Identity.GetSizeRow();k++)
      if (Identity[k][k].real()>1.5){
        std::cout<<"Inversion of covMatrix failed"<<std::endl;
        std::cout<<cov<<std::endl;
        getchar();
      }
      else
        if(k==0&&writeOutCov==true)
          std::cout<<"inversion of cov-Matrix went fine" <<std::endl;
    // data=cov;
    //   invert(covWithOutWeight);
    //     globalCov=covWithOutWeight;

    J=Complex(0.0,0.0);

    if (writeOutCov==true){
      // cov=covWithOutWeight;

      if (actNrParameter==3){      
        *confInterval<<parameter_[1]+parameterInitial_[1]*std::sqrt(cov[0][0].real()*0.115*0.115)<<"  ";
        //        std::cout<<parameter[1]+parameter[1]*std::sqrt(cov[1][1].real()*0.115*0.115)<<std::endl;
        *confInterval<<parameter_[1]<<"  ";
        *confInterval<<parameter_[1]-parameterInitial_[1]*std::sqrt(cov[0][0].real()*0.115*0.115)<<"  ";
        //        std::cout<<parameter[1]-1000.0*std::sqrt(cov[0][0].real()*parameter[1]*0.115*0.115)<<std::endl;
         
         
        *confInterval<<parameter_[7]+parameterInitial_[7]*std::sqrt(cov[1][1].real()*0.115*0.115)<<"  ";
        *confInterval<<parameter_[7]<<"  ";
        *confInterval<<parameter_[7]-parameterInitial_[7]*std::sqrt(cov[1][1].real()*0.115*0.115)<<"  ";
         
         
        *confInterval<<parameter_[9]+parameterInitial_[9]*std::sqrt(cov[2][2].real()*0.115*0.115)<<"  ";
        *confInterval<<parameter_[9]<<"  ";
        *confInterval<<parameter_[9]-parameterInitial_[9]*std::sqrt(cov[2][2].real()*0.115*0.115)<<"  ";
      }

      if (actNrParameter==10)
        for (UInt i=0;i<10;i++){
          *confInterval<<parameter_[i]+parameterInitial_[i]*std::sqrt(cov[i][i].real()*0.155*0.155)<<"  ";
          *confInterval<<parameter_[i]<<"  ";
          *confInterval<<parameter_[i]-parameterInitial_[i]*std::sqrt(cov[i][i].real()*0.155*0.155)<<"  ";
        }

      if (actNrParameter==9){
        for (UInt ii=0;ii<10;ii++){
          if(ii<=1){
            *confInterval<<parameter_[ii]+parameterInitial_[ii]*sqrt(cov[ii][ii].real()*2.5*2.5)<<"  ";
            *confInterval<<parameter_[ii]<<"  ";
            *confInterval<<parameter_[ii]-parameterInitial_[ii]*sqrt(cov[ii][ii].real()*2.5*2.5)<<"  ";
          }
          else if (ii>2){
            *confInterval<<parameter_[ii]+parameterInitial_[ii]*sqrt(cov[ii-1][ii-1].real()*2.5*2.5)<<"  ";
            *confInterval<<parameter_[ii]<<"  ";
            *confInterval<<parameter_[ii]-parameterInitial_[ii]*sqrt(cov[ii-1][ii-1].real()*2.5*2.5)<<"  ";
          }
           
//           if(ii<=1){
//             *confInterval<<parameter_[ii]+parameterInitial_[ii]*sqrt(cov[ii][ii].real()*0.01)<<"  ";
//             *confInterval<<parameter_[ii]<<"  ";
//             *confInterval<<parameter_[ii]-parameterInitial_[ii]*sqrt(cov[ii][ii].real()*0.01)<<"  ";
//           }
//           else if (ii>2){
//             *confInterval<<parameter_[ii]+parameterInitial_[ii]*sqrt(cov[ii-1][ii-1].real()*0.01)<<"  ";
//             *confInterval<<parameter_[ii]<<"  ";
//             *confInterval<<parameter_[ii]-parameterInitial_[ii]*sqrt(cov[ii-1][ii-1].real()*0.01)<<"  ";
             
//           }
           
        }
      }
       
      if (actNrParameter==10)
        for (UInt i=0;i<10;i++){
          *confInterval<<parameter_[i]+parameterInitial_[i]*std::sqrt(cov[i][i].real()*0.155*0.155)<<"  ";
          *confInterval<<parameter_[i]<<"  ";
          *confInterval<<parameter_[i]-parameterInitial_[i]*std::sqrt(cov[i][i].real()*0.155*0.155)<<"  ";
        }

      *confInterval<<std::endl;
      //cov=data;
    }


    // Averaged criterion for all parameters
    if(false){
      for(UInt parInd=0;parInd<actNrParameter+actNrParameterC;parInd++)
        J+=cov[parInd][parInd];
      J/=(actNrParameter+actNrParameterC);
    }
    // Criterion just for one parameter, here eps_11
    if(true)
      J=cov[7][7];

    //   Double sum1=0.0;
    //     Double sum2=0.0;

    //  for (UInt actFreq=0;actFreq<nrfreq_;actFreq++){
    //       sum1+=std::min(0.0,rhos[actFreq])*std::min(0.0,rhos[actFreq]);
    //       sum2+= (1-std::max(1.0,rhos[actFreq]))*(1-std::max(1.0,rhos[actFreq]));
    //       //        }
    //     }

    //     if (sum1!=0.0||sum2!=0){
    //       std::cout<<" sum1 ... sum2 " <<std::endl;
    //       std::cout<<1.0/penaltyFactor1 * sum1 << " ... " << 1.0/penaltyFactor1*sum2 <<std::endl;
    //     }

    //     J+=Complex((1.0/penaltyFactor1)*sum1+(1.0/penaltyFactor2)*sum2,0.0);
    //  penaltyFactor1*=10;
    //      penaltyFactor2*=10;


  }// end createJRHo

  void piezoParamIdent::createGradientRho(Vector<Complex> & grad, Double dOmega){

    Complex J1, J2;
    createJRho(J1,false);
    std::cout<<"Value of J1 = "<<J1 <<std::endl;
    Vector<Double> rhosTemp;
    rhosTemp.Resize(nrfreq_);
    rhosTemp.Init();
    grad.Resize(nrfreq_);
    grad.Init();
    rhosTemp=rhos;

    for(UInt actFreq=0; actFreq<nrMeasuredData; actFreq ++){
      rhos[actFreq] = 1.000001*rhos[actFreq];
      createJRho(J2,false);
      grad[actFreq]=(J2-J1)/(rhos[actFreq]-rhosTemp[actFreq]);
      rhos=rhosTemp;
    }
    //     std::cout<<"gradient von rhoJ ..."<<std::endl;
    //     std::cout<<grad<<std::endl;

  }


  void piezoParamIdent::descentMethodRho(Complex & functional){

    UInt maxNumberDescentIterations=1;
    Vector<Complex> grad;
    Vector<Double> rhosOld;
    rhosOld.Resize(nrMeasuredData);
    rhosOld=rhos;
    // suits for thickness mode and 0.1mm radius
    //    Double lambda=1.25e+4; // for 10 parameters withoutStepWidth

    // suits for thickness mode and 0.05mm radius
    Double lambda=0.9e+2; // for 10 parameters withoutStepWidth

    // for radial, planar mode
    //    Double lambda=-1.0e-20; // for 10 parameters withoutStepWidth


    
    //    Double lambda=1.0e-1; // for 3 parameters
    //    Double lambda=1.0e-2; // for 4 parameters
    //    Double lambda=1.0e-2; // for 6 parameters

    Complex J_old,J;

    Double dOmega=0.0001;
    for (UInt descIter=0;descIter<maxNumberDescentIterations;descIter++){

      createJRho(J_old,false);
      
      J=J_old;
      
      //      lambda=2*lambda;
      createGradientRho(grad, dOmega);

      std::cout<<"gradient J(omega)"<<std::endl;
      std::cout<<grad<<std::endl;

      while (J_old.real()<=J.real()){

        for(UInt fr=0;fr<nrMeasuredData;fr++){
          rhos[fr]=rhos[fr]-lambda*grad[fr].real();
          //  if (rhos[fr]>1.0)
          //             rhos[fr]=1.0;
          //           else if (rhos[fr]<0.0)
          //             rhos[fr]=0.0;
        }
        createJRho(J,false);
        std::cout<<"NEW Functional J = "<<J.real()<<std::endl;
        std::cout<<"OLD Functional J = "<<J_old.real()<<std::endl;

        if (J.real()>=J_old.real()){
          lambda=0.25*lambda;
          rhos=rhosOld;
          std::cout<<"lambda = " << lambda << std::endl;
        }
        rhosOld=rhos;
      }
      std::cout<<" rhos vor Projektion " <<std::endl;
      std::cout<<rhos<<std::endl;

      for(UInt fr=0;fr<nrMeasuredData;fr++){
        if (rhos[fr]>1.0)
          rhos[fr]=1.0;
        else if (rhos[fr]<0.0)
          rhos[fr]=0.0;
      }
      std::cout<<" rhos nach Projektion " <<std::endl;
      std::cout<<rhos<<std::endl;
      
      Double rhoInt=0.0;
      *rhosOut<<J.real()<<"  ";
      
      for (UInt actFreq=0;actFreq<nrMeasuredData;actFreq++)
        rhoInt+=rhos[actFreq];
      rhoInt/=nrMeasuredData;

    //   std::cout<<"rhoInt = "<< rhoInt<<std::endl;
//       *rhosOut<<rhoInt;

      Integer nrfreqAll;
      myParam_->Get( "numFreq", nrfreqAll);

      for (UInt actFreq=0;actFreq<nrfreqAll;actFreq++)
        if (actFreq<nrMeasuredData)
          *rhosOut<<rhos[actFreq]<<"  " ;
        else
          *rhosOut<<rhos[actFreq]<<"  0.0  " ;
      *rhosOut<<std::endl;

      if (rhoInt>=0.85){
        std::cout<<" Sum rho_i exceeds upper bound ... -> break min J(w) " <<std::endl;
        break;
      }
    }


  }// end descentMethodRho

}

