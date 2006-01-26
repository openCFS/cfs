#include "PDE/SinglePDE.hh" 
#include "piezoParamIdent.hh"
#include <iostream>

namespace CoupledField
{

  // xxxxxxxxxxxxxxxxxxxxxxxxxx nuMethods  - ITERATION xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  void piezoParamIdent::optimalExpDesignVarNrFreqs(){

    ENTER_FCN("piezoParamIdent::optimalDesignVarNrFreqs");

    Boolean firstTime=TRUE;
    parameterInitial = parameter;
 //    std::cout<<"parameter"<<std::endl;
//     std::cout<<parameter<<std::endl;

    residuumParIdentOld = residuumParIdent=2000;


    for(UInt optExpVarFreqCout=0; optExpVarFreqCout<35; optExpVarFreqCout++){

      Vector<Double> newFreqs;

      if (optExpVarFreqCout==0){

        Double hOmega = (stopfreq-startfreq)/nrfreq;
        freqs.Resize(nrfreq);
        nrMeasuredData=nrfreq;

        // writes out confidence intervals for maximal number of frequencies
        if (FALSE){
          for (UInt actFreq=0;actFreq<nrfreq;actFreq++)
            freqs[actFreq] = startfreq+actFreq*hOmega;
          readInMeasurement(newFreqs);  // get frequencies
          calc_measuredCharge(freqs, real, imag, y_hat); // indispensable for calculation of conf intervals
          writeOutConfInterval(); 
        }


        // to compare write out confidence interval for twelfe arbitrarily,
        // but equidistantly chosen frequencies        
        if (FALSE){
          UInt nrfreqOld = nrfreq;
          hOmega = (stopfreq-startfreq)/12;
          freqs.Resize(nrfreq);
          nrMeasuredData=12;
          nrfreq=12;
          for (UInt actFreq=0;actFreq<12;actFreq++)
            freqs[actFreq] = startfreq+actFreq*hOmega;
          readInMeasurement(newFreqs);
          calc_measuredCharge(freqs, real, imag, y_hat); // out of new measurements
          writeOutConfInterval();

          // reset everything to avoid erroneous calculations in the following
          nrfreq=nrfreqOld;
          hOmega = (stopfreq-startfreq)/nrfreq;
          freqs.Resize(nrfreq);
          nrMeasuredData=nrfreq;
        }

        for (UInt actFreq=0;actFreq<nrfreq;actFreq++)
          freqs[actFreq] = startfreq+actFreq*hOmega;
        
        readInMeasurement(newFreqs);
        calc_measuredCharge(freqs, real, imag, y_hat); // out of new measurements
        
        if (firstTime==TRUE){
          std::cout<<"++ Calculate weighting function rho(omega) for the following frequencies"<<std::endl;
          std::cout<<newFreqs<<std::endl;
          firstTime=FALSE;
        }
        
        UInt maxNrBreakpoints = nrfreq;
        
        Complex J;
        
        Double dOmega = 0.0001;
        rhos.Resize(nrfreq);
        
        for(UInt i=0;i<nrfreq;i++)
          rhos[i]=0.01;
        
        std::cout<<"++ Minimize the functional J(omega) ... "<<std::endl;
        descentMethodRho(J);

        // with all frequencies:
        globalIndexSet.Resize(nrfreq);
        for (UInt k=0;k<nrfreq;k++)
          globalIndexSet[k]=1;

        //        writeOutConfInterval();
        std::cout<<"++ Calculated rhos: ... "<<std::endl;        
        std::cout<<rhos<<std::endl;
        
        UInt nrNuMethods=0;

        // choose 8 highest values ...
        Vector<Double> highestFreqs;
        Vector<Double> highestRhos;
        UInt maxIndex;
        Double sumRhos=0.0;
        highestFreqs.Resize(nrfreq);
        highestRhos.Resize(nrfreq);
        
        omegaDiffVec.Resize(nrfreq);
        
        for (UInt i=0;i<nrfreq;i++){
          maxIndex=0;
          highestFreqs[i]=rhos[0];
          
          for (UInt j=0;j<rhos.GetSize();j++){
            if (highestRhos[i]<=rhos[j]){ // hier muss " <= " hin!!!
              highestFreqs[i]=freqs[j];
              highestRhos[i]=rhos[j];
              maxIndex=j;
            }
          }
          rhos[maxIndex]=0.0;
          
          Double omegaDiff=nrfreq*hOmega;
          
          // Get Neighbor of \omega_i
          for(UInt nI=0;nI<highestFreqs.GetSize();nI++){
            if (std::abs((highestFreqs[i]-highestFreqs[nI]))<omegaDiff&&(nI!=i))
              omegaDiff=std::abs(highestFreqs[i]-highestFreqs[nI]);
            
          }
          //  std::cout<<" omegaDiff = " <<omegaDiff<<std::endl;
          omegaDiffVec[i]=omegaDiff;
          sumRhos+=highestRhos[i]*omegaDiff;
          std::cout<<" - sumRhos = " << sumRhos<<std::endl;
          
          if ((sumRhos>33.5*hOmega && i>10) || i>=11){
            //        if (sumRhos>6*hOmega && i>3){
            std::cout<<"sumRhos = "<<sumRhos <<" ecxeeds 33.5*hOmega = " 
                     << 33.5*hOmega << "; hOmega = " << hOmega<< std::endl;
            //          getchar();
            break;
          }
        }
            
        globalIndexSet.Resize(nrfreq);
        
        Integer howManyFreqs=0;
        for (UInt i=0;i<nrfreq;i++)
          if (highestFreqs[i]!=0.0){
            howManyFreqs++;
          }
                    
        std::cout<<"howManyFreqs = "<< howManyFreqs<<std::endl;
        nrfreq=howManyFreqs;
        
        // sort highestFreqs
        Double dv,dr,dO;
        Integer k, ii;
        
        
        for (k=0; k<howManyFreqs-1; ++k)
          for (ii=k+1; ii<howManyFreqs; ++ii)
            if (highestFreqs[ii] < highestFreqs[k])
              {
                dv = highestFreqs[k];
                dr = highestRhos[k];
                dO = omegaDiffVec[k];
                highestFreqs[k] = highestFreqs[ii];
                highestRhos[k] = highestRhos[ii];
                omegaDiffVec[k]=omegaDiffVec[ii];
                highestFreqs[ii] = dv;
                highestRhos[ii] = dr;
                omegaDiffVec[ii]=dO;
              }

        for (UInt i=0; i<highestFreqs.GetSize(); i++)
          for (UInt j=0;j<howManyFreqs;j++)
            if (highestFreqs[j]==freqs[i])
              globalIndexSet[i]=1;

        std::cout<<" ++ finished to solve optimal experiment design " <<std::endl;
        std::cout<<"    the following results might be used as initial guesses "<<std::endl;
        std::cout<<"    for optimal experiment design with fixed number of frequencies:" <<std::endl;

        std::cout<<"highestRhos"<<std::endl;
        std::cout<<highestRhos<<std::endl;
      
        std::cout<<"highestFreqs"<<std::endl;
        std::cout<<highestFreqs<<std::endl;

      //       std::cout<<"omegaDiffVec"<<std::endl;
      //       std::cout<<omegaDiffVec<<std::endl;


        for(UInt i=0; i<highestFreqs.GetSize();i++)
          *optimalFreqs<<highestFreqs[i]<<"  ";
        *optimalFreqs<<std::endl;

        std::cout<<"Press any key to continue with paramter identificaton"<<std::endl;        
        getchar();

        if (optExpVarFreqCout==0){

          nrMeasuredData=howManyFreqs;
        
          rhos.Resize(howManyFreqs);
          freqs.Resize(howManyFreqs);
          real.Resize(howManyFreqs);
          imag.Resize(howManyFreqs);
          
          rhos=highestRhos;
          freqs=highestFreqs;
                   
          readInMeasurement(newFreqs);

          calc_measuredCharge(freqs, real, imag, y_hat); // out of new measurements

          writeOutConfInterval();

          getchar();

          newtonCounter=0;
          
          if (FALSE){
            while (nrNuMethods<maxNumberNewtonLoops && residuumParIdent > 0.1*residuumParIdentOld){
              //      while (nrNuMethods<maxNumberNewtonLoops){
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

          residuumParIdentOld = residuumParIdent;
        //  nrfreq++;
        } // end if if (optExpVarFreqCout==2){

      //   Complex JTemp;
//         createJRho(JTemp,FALSE);
//         writeOutConfInterval();

      }      
      else { // do optimization of fixed number of frequencies!
        projGradientFlags.Resize(nrMeasuredData);
        parameterInitial = parameter;
        residuumParIdentOld = residuumParIdent=0.5;
       

        std::cout<<"optimum experiment desing with "<<nrMeasuredData 
                 <<" frequencies and "<<actNrParameter+actNrParameterC
                 << " parameter! " << std::endl;
        Vector<Double> newFreqs;
        
        readInMeasurement(newFreqs);
        std::cout<<newFreqs<<std::endl;
        
        for(UInt fr=0;fr<nrMeasuredData;fr++)
          *impedCurve<< freqs[fr]<<"  ";
        *impedCurve<<std::endl;
        
        calc_measuredCharge(freqs, real, imag, y_hat); // out of new measurements
        
        //  actNrParameterC=0;
        MaterialData * ptMaterial;
        
        if(directCoupling==TRUE)
          ptMaterial=ptPDE1_->getPDEMaterialData();   // Pointer to MaterialData
        else
          ptMaterial=ptMyPDE_->getPDEMaterialData();   // Pointer to MaterialData
        
        // Frequency bounds
        
        // fr = 4000m/s / 2*thickness
        
        fr_=4000/(2*thickness);
        // minus 5 percent
        fr_=0.94*fr_;
        // fa = 1/thickness * sqrt(c_33/rho)
        
        Double rho;
        Vector<Complex> jacobi;
        Complex functional;
        
        rho = ptMaterial->GetDensity();
        
        fa_ = 1.0/(2*thickness) * std::sqrt(1.3e+11/rho);
        fa_=1.2*fa_;
        
        std::cout<<"++ Bounds for resonance and antiresonace frequency: "<<std::endl;
        std::cout<<"fr= "<< fr_ <<std::endl;
        std::cout<<"fa= "<< fa_ <<std::endl;
        
        
        for (UInt nrOptExpSteps=0; nrOptExpSteps<20; nrOptExpSteps++){ 
          UInt nrNuMethods=0;
      // That at least we perform one descent step
          normGradient=2.0;
          normGradientOld=1.0;

          descentMethod(functional);
          
          readInMeasurement(newFreqs);
          
          std::cout<<newFreqs<<std::endl;
          
          for(UInt fr=0;fr<nrMeasuredData;fr++)
            *impedCurve<< freqs[fr]<<"  ";
          *impedCurve<<std::endl;
          
          calc_measuredCharge(freqs, real, imag, y_hat); // out of new measurements
          newtonCounter=0;
          
          while (nrNuMethods<maxNumberNewtonLoops && residuumParIdent > 0.2*residuumParIdentOld){
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
          residuumParIdentOld = residuumParIdent;
      
        } // end nrOptExpSteps
    
        std::cout<<"\n END OF OPTIMAL EXPERIMENT DESIGN" <<std::endl;

      } //end of else (fixed nr of frequencies!)
      }
      
    }

 void piezoParamIdent::writeOutConfInterval(){
    ENTER_FCN("piezoParamIdent::writeOutConfInterval");

    std::cout<<"++ Write out confidence intervals into file confInterval.dat " <<std::endl;

//     std::cout<<"globalIndexSet:"<<std::endl;
//     std::cout<<globalIndexSet<<std::endl;

    std::cout.precision( 2 );


    Double deltaLocal=1.01;

    Vector<Complex> jacobi;
    Vector<Complex> jacobiH;
  
    Matrix<Complex>covTempWithOutWeight(actNrParameter+actNrParameterC, actNrParameter+actNrParameterC);
    Matrix<Complex>covWithOutWeight(actNrParameter+actNrParameterC, actNrParameter+actNrParameterC);

    for (UInt actFreq=0;actFreq<nrfreq;actFreq++){

      // std::cout<<freqs<<std::endl;
      createJacobian(jacobi, freqs[actFreq]);
      jacobiH.Resize(actNrParameter+actNrParameterC);
      if(actFreq==0){
        std::cout<<"jacobi:"<<std::endl;
        std::cout<<jacobi<<std::endl;
      }

      for (UInt i=0;i<jacobi.GetSize();i++)
        jacobiH[i]=Complex(jacobi[i].real(),-jacobi[i].imag());

      globalJacobi=jacobi;
      globalJacobiH=jacobiH;

      //    for(UInt i=0; i<globalIndexSet.GetSize(); i++)
      // if (globalIndexSet[actFreq]==1){
      //     if (freqs[actFreq]<=fr_||freqs[actFreq]>=fa_){

          Double stepWidth=0.0;
          stepWidth=0.5*(freqs[std::min(actFreq+1,nrfreq-1)]-freqs[std::max((Integer)actFreq-1,0)]);
       //    std::cout<<"right Freq"<<std::endl;
//           std::cout<<freqs[std::min(actFreq+1,nrfreq-1)]<<std::endl;
//           std::cout<<"left Freq"<<std::endl;
//           std::cout<<freqs[std::max((Integer)actFreq-1,0)]<<std::endl;
//           //          stepWidth/=1.0e+05;
//           //   std::cout<<0.5*Double(nrfreq)<<std::endl;

//           std::cout<<"stepWidth: = "<<stepWidth<<std::endl;

          //     std::cout<<" _________________"<<std::endl;
          for(UInt i=0;i<actNrParameter+actNrParameterC;i++)
            for(UInt j=0;j<actNrParameter+actNrParameterC;j++){
              
//                if(actFreq==0||actFreq==nrfreq)
//                  covTempWithOutWeight[i][j]=0.5*stepWidth*globalJacobiH[i]*globalJacobi[j]/
//                    (freqs[actFreq]*(deltaLocal)*(deltaLocal)*std::abs(y_hat[actFreq])*std::abs(y_hat[actFreq]));
//                else
                 covTempWithOutWeight[i][j]= stepWidth*globalJacobiH[i]*globalJacobi[j]/
                   ((deltaLocal)*(deltaLocal)*std::abs(y_hat[actFreq])*std::abs(y_hat[actFreq]));

            //   covTempWithOutWeight[i][j]=globalJacobiH[i]*globalJacobi[j]/
//                 ((deltaLocal)*(deltaLocal)*std::abs(y_hat[actFreq])*std::abs(y_hat[actFreq]));


              //    if (i==j)
              //    covWithOutWeight[i][i]+=covWithOutWeight[i][i];
            }
          // }
          covWithOutWeight+=covTempWithOutWeight;

//            for(UInt i=0;i<actNrParameter+actNrParameterC;i++)
//              for(UInt j=0;j<actNrParameter+actNrParameterC;j++)
//                covWithOutWeight[i][j]=covWithOutWeight[i][j]/Complex(Double(nrfreq),0);
//       /
            //        std::cout<<"cov"<<std::endl;
      //        std::cout<<cov<<std::endl;
        // } // end if global IndexSet .
    } // end for actFreq ..

    
    globalCov=covWithOutWeight;
//     std::cout<<"globalCov = cov wihtout Weight :"<<std::endl;
//     std::cout<<covWithOutWeight<<std::endl;

    invert(covWithOutWeight);
    Matrix<Complex> testId;
    testId=globalCov*covWithOutWeight;
    std::cout<<testId<<std::endl;
  
    globalCov=covWithOutWeight;

    // std::cout<<"globalCov = cov wihtout Weight inverted:"<<std::endl;
//     std::cout<<globalCov<<std::endl;
   
    // getchar();

  //   for(UInt i=0;i<actNrParameter+actNrParameterC;i++)
//       for(UInt j=0;j<actNrParameter+actNrParameterC;j++)
//         covWithOutWeight[i][j]=1.0e-9*covWithOutWeight[i][j];

    
    if (actNrParameter==9){
      std::cout<<"++ Confidence intervals will be written into confInterval.dat " <<std::endl;
      for (UInt ii=0;ii<10;ii++){
        if(ii<=1){
          *confInterval<<parameter[ii]+parameterInitial[ii]*sqrt(globalCov[ii][ii].real()*2.5*2.5)<<"  ";
          *confInterval<<parameter[ii]<<"  ";
          *confInterval<<parameter[ii]-parameterInitial[ii]*sqrt(globalCov[ii][ii].real()*2.5*2.5)<<"  ";
        }
        else if (ii>2){
          *confInterval<<parameter[ii]+parameterInitial[ii]*sqrt(globalCov[ii-1][ii-1].real()*2.5*2.5)<<"  ";
          *confInterval<<parameter[ii]<<"  ";
          *confInterval<<parameter[ii]-parameterInitial[ii]*sqrt(globalCov[ii-1][ii-1].real()*2.5*2.5)<<"  ";
        }
        
      }
      *confInterval<<std::endl;
    }
    std::cout.precision(6);
    
 } // end writeOutConfInterval

  void piezoParamIdent::createJRho(Complex &J, Boolean writeOutCov){
    ENTER_FCN("piezoParamIdent::createJRho");

    Double  hOmega = (stopfreq-startfreq)/nrfreq;
    delta=0.00;
    Double deltaLocal=0.01;

    //    Complex J;
    J=Complex(0.0,0.0);

    Vector<Complex> jacobi;
    Vector<Complex> jacobiH;

    Double penaltyFactor1=1.0;
    Double penaltyFactor2=1.0;
    // no penalization

    fr_=4000/(2*thickness);
    // minus 5 percent
    fr_=0.94*fr_;
    // fa = 1/thickness * sqrt(c_33/rho)
   
    Double rho;
    rho = ptMaterial->GetDensity();

    fa_ = 1.0/(2*thickness) * std::sqrt(1.3e+11/rho);
    fa_=1.2*fa_;

    // for thicker probes:
//     fa_=1.2e+06;
//     fr_=0.85e+06;

    //     std::cout<<"++ Bounds for resonance and antiresonace frequency"<<std::endl;
    //     std::cout<<"fr:"<<std::endl;
    //     std::cout<<fr_<<std::endl;

    //     std::cout<<"fa:"<<std::endl;
    //     std::cout<<fa_<<std::endl;

    //    MaterialData * ptMaterial=ptMyPDE_->getPDEMaterialData();   // Pointer to MaterialData
    Complex F_hat_incr1,F_hat_incr2;

    Matrix<Complex> cov(actNrParameter+actNrParameterC, actNrParameter+actNrParameterC);
    Matrix<Complex> covTemp(actNrParameter+actNrParameterC, actNrParameter+actNrParameterC);
    Matrix<Complex>covTempWithOutWeight(actNrParameter+actNrParameterC, actNrParameter+actNrParameterC);
    Matrix<Complex>covWithOutWeight(actNrParameter+actNrParameterC, actNrParameter+actNrParameterC);

    for (UInt actFreq=0;actFreq<nrfreq;actFreq++){
    
      createJacobian(jacobi, freqs[actFreq]);

      jacobiH.Resize(actNrParameter+actNrParameterC);
       
      //       std::cout<<jacobi<<std::endl;
      for (UInt i=0;i<jacobi.GetSize();i++)
        jacobiH[i]=Complex(jacobi[i].real(),-jacobi[i].imag());

      globalJacobi=jacobi;
      globalJacobiH=jacobiH;

      //    if (freqs[actFreq]<=fr_||freqs[actFreq]>=fa_){   

        for(UInt i=0;i<actNrParameter+actNrParameterC;i++)
          for(UInt j=0;j<actNrParameter+actNrParameterC;j++){

            covTempWithOutWeight[i][j]=jacobiH[i]*jacobi[j]/(std::abs(y_hat[actFreq])*std::abs(y_hat[actFreq]));

            //           covTemp[i][j]=jacobiH[i]*jacobi[j]/(delta*delta*std::abs(y_hat[actFreq])*std::abs(y_hat[actFreq]));
            if(actFreq==0||actFreq==nrfreq)
              covTemp[i][j]=rhos[actFreq]*hOmega*0.5*jacobiH[i]*jacobi[j]/
                ((deltaLocal)*(deltaLocal)*std::abs(y_hat[actFreq])*std::abs(y_hat[actFreq]));

            //              covTemp[i][j]=rhos[actFreq]*hOmega*0.5*jacobiH[i]*jacobi[j]/
            //  (freqs[actFreq]*(deltaLocal)*(deltaLocal)*std::abs(y_hat[actFreq])*std::abs(y_hat[actFreq]));

            //                 covTemp[i][j]=rhos[actFreq]*hOmega*0.5*jacobiH[i]*jacobi[j]/
            //                   (0.8*freqs[actFreq]*(1+delta)*(1+delta)*std::abs(y_hat[actFreq])*std::abs(y_hat[actFreq]));

            else 
              covTemp[i][j]=rhos[actFreq]*hOmega*jacobiH[i]*jacobi[j]/
                ((deltaLocal)*(deltaLocal)*std::abs(y_hat[actFreq])*std::abs(y_hat[actFreq]));

            //              covTemp[i][j]=rhos[actFreq]*hOmega*jacobiH[i]*jacobi[j]/
            //  (freqs[actFreq]*(deltaLocal)*(deltaLocal)*std::abs(y_hat[actFreq])*std::abs(y_hat[actFreq]));

            //                covTemp[i][j]=rhos[actFreq]*hOmega*jacobiH[i]*jacobi[j]/
            //                  (freqs[actFreq]*(1+delta)*(1+delta)*std::abs(y_hat[actFreq])*std::abs(y_hat[actFreq]));

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
        if(k==0&&writeOutCov==TRUE)
          std::cout<<"inversion of cov-Matrix went fine" <<std::endl;
  // data=cov;
  //   invert(covWithOutWeight);
//     globalCov=covWithOutWeight;

    J=Complex(0.0,0.0);

    if (writeOutCov==TRUE){
      // cov=covWithOutWeight;

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

      if (actNrParameter==9){
        for (UInt ii=0;ii<10;ii++){
          if(ii<=1){
            *confInterval<<parameter[ii]+parameterInitial[ii]*sqrt(cov[ii][ii].real()*2.5*2.5)<<"  ";
            *confInterval<<parameter[ii]<<"  ";
            *confInterval<<parameter[ii]-parameterInitial[ii]*sqrt(cov[ii][ii].real()*2.5*2.5)<<"  ";
          }
          else if (ii>2){
            *confInterval<<parameter[ii]+parameterInitial[ii]*sqrt(cov[ii-1][ii-1].real()*2.5*2.5)<<"  ";
            *confInterval<<parameter[ii]<<"  ";
            *confInterval<<parameter[ii]-parameterInitial[ii]*sqrt(cov[ii-1][ii-1].real()*2.5*2.5)<<"  ";
          }
           
          //  if(ii<=1){
          //              *confInterval<<parameter[ii]+parameterInitial[ii]*sqrt(cov[ii][ii].real()*0.01)<<"  ";
          //              *confInterval<<parameter[ii]<<"  ";
          //              *confInterval<<parameter[ii]-parameterInitial[ii]*sqrt(cov[ii][ii].real()*0.01)<<"  ";
          //            }
          //            else if (ii>2){
          //              *confInterval<<parameter[ii]+parameterInitial[ii]*sqrt(cov[ii-1][ii-1].real()*0.01)<<"  ";
          //              *confInterval<<parameter[ii]<<"  ";
          //              *confInterval<<parameter[ii]-parameterInitial[ii]*sqrt(cov[ii-1][ii-1].real()*0.01)<<"  ";
             
          //            }
           
        }
      }
       
      if (actNrParameter==10)
        for (UInt i=0;i<10;i++){
          *confInterval<<parameter[i]+parameterInitial[i]*std::sqrt(cov[i][i].real()*0.155*0.155)<<"  ";
          *confInterval<<parameter[i]<<"  ";
          *confInterval<<parameter[i]-parameterInitial[i]*std::sqrt(cov[i][i].real()*0.155*0.155)<<"  ";
        }

      *confInterval<<std::endl;
      //cov=data;
    }



    for(UInt parInd=0;parInd<actNrParameter+actNrParameterC;parInd++)
      J+=cov[parInd][parInd];
    J/=(actNrParameter+actNrParameterC);

  //   Double sum1=0.0;
//     Double sum2=0.0;

   //  for (UInt actFreq=0;actFreq<nrfreq;actFreq++){
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
    ENTER_FCN("piezoParamIdent::createGradientRho");

    Complex J1, J2;
    createJRho(J1,FALSE);
    std::cout<<"Value of J1 = "<<J1 <<std::endl;
    Vector<Double> rhosTemp;
    rhosTemp.Resize(nrfreq);
    grad.Resize(nrfreq);
    rhosTemp=rhos;

    for(UInt actFreq=0; actFreq<nrMeasuredData; actFreq ++){
      rhos[actFreq] = 1.000001*rhos[actFreq];
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
    Double lambda=1.25e+4; // for 10 parameters withoutStepWidth
    
    //    Double lambda=1.0e-1; // for 3 parameters
    //    Double lambda=1.0e-2; // for 4 parameters
    //    Double lambda=1.0e-2; // for 6 parameters

    Complex J_old,J;

    Double dOmega=0.0001;
    for (UInt descIter=0;descIter<maxNumberDescentIterations;descIter++){

      createJRho(J_old,FALSE);
      
      J=J_old;
      
      lambda=2*lambda;
      createGradientRho(grad, dOmega);

      while (J_old.real()<=J.real()){

        for(UInt fr=0;fr<nrMeasuredData;fr++){
          rhos[fr]=rhos[fr]-lambda*grad[fr].real();
         //  if (rhos[fr]>1.0)
//             rhos[fr]=1.0;
//           else if (rhos[fr]<0.0)
//             rhos[fr]=0.0;
        }
        createJRho(J,FALSE);
        std::cout<<"NEW Functional J = "<<J.real()<<std::endl;
        std::cout<<"OLD Functional J = "<<J_old.real()<<std::endl;

        if (J.real()>=J_old.real()){
          lambda=0.5*lambda;
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
