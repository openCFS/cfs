// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;




#include "PDE/SinglePDE.hh" 
#include "piezoParamIdent.hh"


namespace CoupledField
{

  
  // xxxxxxxxxxxxxxxxxxxxxxxxxxxxx least square xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  void piezoParamIdent::leastSquare(){ 

    //ptMaterial_=ptMyPDE_->getPDEMaterialData();   // Pointer to MaterialData


    Vector<Double>freqs5;
    if(false){
      Integer numMeasPoints=25;
      Double fmin=0.2e+5;
      Double fmax=2.2e+5;
      Double fincr=(fmax-fmin)/numMeasPoints;
      freqs5.Resize(numMeasPoints);
      nrMeasuredData=numMeasPoints;
      for(UInt i=0;i<numMeasPoints;i++)
        freqs5[i]=fmin+i*fincr;
    }
      
    //    freqs_=freqs5;


    Vector<Complex> F_y(nrMeasuredData);
    Vector <Double> parameter_old(10);
    Vector <Double> parameter_oldC(10);
    Vector <Double> parameter_initial(10);
    Vector <Double> parameter_initialC(10);
    Vector <Double> gradP(actNrParameter+actNrParameterC);
    Vector<Double> newFreqs;
    parameter_old=parameter_;
    parameter_oldC=parameterC_;
    parameter_initial=parameter_;
    parameter_initialC=parameterC_;
    Double normFy, maxres, normFy0,normFy1;
    Integer indPar=0;
    Integer indParC=0;
    Double lambda=50.0;
    bool negFlag=false;      


 //    std::cout<<"actNrParameter"<<actNrParameter <<std::endl;
//     std::cout<<"actNrParameterC"<<actNrParameterC <<std::endl;

    std::cout<<"++ The following parameters will be updated:"<<std::endl;
    std::cout<<"Real parts\n "<<whichParameterToUpdate_ <<std::endl;
    std::cout<<"Imaginary parts\n "<<whichParameterToUpdateC_ <<std::endl;

    calc_measuredCharge(freqs_, real_, imag_, y_hat_); // out of new measurements

    for (UInt iterIndex=0; iterIndex<maxNumberNewtonLoops_;iterIndex++){

      if( imagMaterialParam_ ) {
        updateComplexMaterialData(parameterC_);
      }
      updateMaterialData(parameter_);

      createF(F_hat_, false);

//       std::cout<<"y_hat"<<std::endl;
//       std::cout<<y_hat_<<std::endl;
     
//       std::cout<<"F_hat"<<std::endl;
//       std::cout<<F_hat_<<std::endl;

      if ((iterIndex+1)%10==0){
        Integer nrfreqTemp=100;
        Vector<Double> freqsTemp = freqs_;
        freqs_.Resize(nrfreqTemp);
        Double startFreqTemp;
        startFreqTemp=startfreq_;
        Double freqincr=(stopfreq_-startfreq_)/nrfreqTemp;
        for(UInt i=0;i<nrfreqTemp;i++){
          startFreqTemp+=freqincr;
          freqs_[i]=startFreqTemp;
        }
        if (impedCurve)
          impedCurve->close();
        if(mechDispl)
          mechDispl->close();

        std::ostringstream num;
        num<<iterIndex;
        std::string filename= "imped";
        std::string numString = num.str();
        //        filename.append(numString);
        filename.append(".dat");
        impedCurve = new std::ofstream(filename.c_str(),std::basic_ios<char>::out);
        filename="mechDispl.dat";
        mechDispl = new std::ofstream(filename.c_str(),std::basic_ios<char>::out);

        calcImpedanceCurve();
        freqs_ = freqsTemp;
      }

      
      for (UInt i=0;i<nrMeasuredData;i++)
        F_y[i]=F_hat_[i]-y_hat_[i];
      

      //       std::cout<<"F_hat_:"<<std::endl;
      //       std::cout<<F_hat_<<std::endl;

      //       std::cout<<"y_hat_:"<<std::endl;
      //       std::cout<<y_hat_<<std::endl;

      norm(F_y,normFy0,maxres,y_hat_);
      
      std::cout<<"\n ||F(p)-y|| = " << normFy0<<std::endl;
      std::cout<< "--------------------------------\n"<<std::endl;

      *parLog<<normFy0<<"  ";
      for (UInt par=0;par<nrParameter_;par++)
        if (whichParameterToUpdate_[par]==1){
          *parLog<<parameter_[par]<<"  ";
        }
      for (UInt par=0;par<nrParameter_;par++)
        if (whichParameterToUpdateC_[par]==1){
          *parLog<<parameterC_[par]<<"  ";
        }

      *parLog<<std::endl;

      
      indPar=0;
      indParC=0;

      for (UInt par=0;par<nrParameter_;par++){

        if (whichParameterToUpdate_[par]==1){

          parameter_[par]=1.000001*parameter_[par];
          
          updateMaterialData(parameter_);
          
          createF(F_hat_, false);
          
          for (UInt i=0;i<nrMeasuredData;i++)
            F_y[i]=F_hat_[i]-y_hat_[i];
          
          norm(F_y,normFy,maxres,y_hat_);
          
          gradP[indPar]=normFy0-normFy;
          
          parameter_[par]=parameter_old[par];

          indPar++;
        }

        if (whichParameterToUpdateC_[par]==1){
          
          parameterC_[par]=1.1*parameterC_[par];
          
         
          if( imagMaterialParam_ ) {
            updateComplexMaterialData(parameterC_);
          }
          updateMaterialData(parameter_);
          
          createF(F_hat_, false);
          
          for (UInt i=0;i<nrMeasuredData;i++)
            F_y[i]=F_hat_[i]-y_hat_[i];
          
          norm(F_y,normFy,maxres,y_hat_);
          
          //            gradP[indParC+actNrParameter]=1.0*(normFy0-normFy);
          gradP[indParC+actNrParameter]=0.000001*(normFy0-normFy);
          
          parameterC_[par]=parameter_oldC[par];
          
          indParC++;
        }

      }
      
      computeScaling();     
          
      
      //     for (UInt par=0;par<nrParameter;par++)
      //       gradP[par]=gradP[par]/scaling_[par];
      
      std::cout<<"Gradient of parameters: "<<std::endl;
      std::cout<<gradP<<std::endl;

//       std::cout<<"scalingC " <<std::endl;
//       std::cout<<scalingC_<<std::endl;

      indPar=0;
      lambda=5*lambda;
      for (UInt par=0;par<nrParameter_;par++)
        if (whichParameterToUpdate_[par]==1){
          Integer innerLoopCount2=0;
          Double lambdaLocal=lambda;
          parameter_[par]+=lambda*gradP[indPar]/scaling_[par];
          indPar++;
        }

      indParC=0;
      for (UInt par=0;par<parameterC_.GetSize();par++)
        if (whichParameterToUpdateC_[par]==1){
          parameterC_[par]+=lambda*gradP[indParC+actNrParameter]/scalingC_[par];
          std::cout<<"in Iteration: parameterC_[" << par <<"] = " << parameterC_[par]<<std::endl;
          indParC++;
        }
      
      if( imagMaterialParam_ ) {
        updateComplexMaterialData(parameterC_);
      }
      updateMaterialData(parameter_);

      createF(F_hat_, false);
      for (UInt i=0;i<nrMeasuredData;i++)
        F_y[i]=F_hat_[i]-y_hat_[i];
      norm(F_y,normFy1,maxres,y_hat_);

      
      for (UInt par=0;par<nrParameter_;par++)
        if (parameter_[par]<0.0&&par!=6){
          std::cout<<"parameter( " << par << " ) was negative " <<std::endl;
          negFlag=true;
        }
      
      for (UInt par=0;par<nrParameter_;par++)
        if (parameterC_[par]<0.0){
          std::cout<<"imag part parameter( " << par << " ) was negative " <<std::endl;
          negFlag=true;
        }


      Integer lineSearchCount=0;

      std::cout<<"++ Performing line search ... " <<std::endl;

      while (normFy1>normFy0||negFlag==true){
        if (negFlag==true)
          std::cout<< "reduce lambda due to negative parameters or non reduction of norm ... "<<std::endl;
        parameter_=parameter_old;
        parameterC_=parameter_oldC;
        lambda=0.8*lambda;
        std::cout<<"lambda = " << lambda <<std::endl;
        negFlag=false;

        //avoids that programm stagnates!
        lineSearchCount++;
        if (lineSearchCount>15){
          std::cout<<"! TERMINATION of line search ... suspicious event ...\n " <<std::endl;
          std::cout<<"either we are close to a solution or something went completely wrong.\n" <<std::endl;
          std::cout<<"Check your results (convergence?) and input data!"<<std::endl;
          break;}
        
        indPar=0;
        for (UInt par=0;par<nrParameter_;par++)
          if (whichParameterToUpdate_[par]==1){
            Integer innerLoopCount=0;
            Double lambdaLocal=lambda;
            parameter_[par]+=lambda*gradP[indPar]/scaling_[par];
            indPar++;
          }
        indParC=0;
        
        for (UInt par=0;par<parameterC_.GetSize();par++)
          if (whichParameterToUpdateC_[par]==1){
            parameterC_[par]+=lambda*gradP[indParC+actNrParameter]/scalingC_[par];
            indParC++;
          }      

        if( imagMaterialParam_ ) 
          updateComplexMaterialData(parameterC_);
        
        updateMaterialData(parameter_);
        createF(F_hat_, false);
        for (UInt i=0;i<nrMeasuredData;i++)
          F_y[i]=F_hat_[i]-y_hat_[i];
        norm(F_y,normFy1,maxres,y_hat_);
        
        for (UInt par=0;par<nrParameter_;par++)
          if (parameter_[par]<0.0&&par!=6){
            std::cout<<"parameter( " << par << " ) is negative " <<std::endl;
            negFlag=true;
          }
        for (UInt par=0;par<nrParameter_;par++)
          if (parameterC_[par]<0.0){
            std::cout<<"imaginary part of parameter( " << par << " ) is negative " <<std::endl;
            negFlag=true;
          }

      
      }
      normFy0=normFy1;

      
      std::cout<<"New parameter:"<<std::endl;

      //std::cout<<parameter_<<std::endl;
      for (UInt i=0; i<parameter_.GetSize(); i++){
        std::cout<<i<<") "<< parameter_[i] << " + "<< parameterC_[i] <<"i "<<std::endl;
        *parFinal<<parameter_[i]<<", ";
      }

      *parFinal<<std::endl;
      
      for (UInt i=0; i<parameter_.GetSize(); i++)
        *parFinal<<parameterC_[i]<<", ";

      *parFinal<<std::endl;
            
      //       std::cout<<"parameter_old:"<<std::endl;
      //       std::cout<<parameter_old<<std::endl;
      parameter_old=parameter_;
      parameter_oldC=parameterC_;
    } // end iter index
    



          
  } // end least square

} // end namespace

