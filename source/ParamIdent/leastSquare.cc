// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;




#include "PDE/SinglePDE.hh" 
#include "piezoParamIdent.hh"


namespace CoupledField
{

  
  // xxxxxxxxxxxxxxxxxxxxxxxxxxxxx least square xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  void piezoParamIdent::leastSquare(){ 
    ENTER_FCN("piezoParamIdent::leastSquare");    

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


    calc_measuredCharge(freqs_, real_, imag_, y_hat_); // out of new measurements

    for (UInt iterIndex=0; iterIndex<maxNumberNewtonLoops_;iterIndex++){

      if( imagMaterialParam_ ) {
        updateComplexMaterialData(parameterC_);
      }
      updateMaterialData(parameter_);

      createF(F_hat_, false);

      std::cout<<"y_hat"<<std::endl;
      std::cout<<y_hat_<<std::endl;
     
      std::cout<<"F_hat"<<std::endl;
      std::cout<<F_hat_<<std::endl;

      if ((iterIndex+1)%5==0){
        Integer nrfreqTemp=200;
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
        std::string filename= "imped.dat";
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
      
      std::cout<<"||F(p)-y|| = " << normFy0<<std::endl;
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
          
          parameterC_[par]=1.000001*parameterC_[par];
          
         
          if( imagMaterialParam_ ) {
            updateComplexMaterialData(parameterC_);
          }
          updateMaterialData(parameter_);
          
          createF(F_hat_, false);
          
          for (UInt i=0;i<nrMeasuredData;i++)
            F_y[i]=F_hat_[i]-y_hat_[i];
          
          norm(F_y,normFy,maxres,y_hat_);
          
          //            gradP[indParC+actNrParameter]=1.0*(normFy0-normFy);
          gradP[indParC+actNrParameter]=100.0*(normFy0-normFy);
          
          parameterC_[par]=parameter_oldC[par];
          
          indParC++;
        }

      }
      
      computeScaling();     
          
      
      //     for (UInt par=0;par<nrParameter;par++)
      //       gradP[par]=gradP[par]/scaling_[par];
      
      std::cout<<"Gradient of parameters: "<<std::endl;
      std::cout<<gradP<<std::endl;

      indPar=0;
      lambda=100*lambda;
      for (UInt par=0;par<nrParameter_;par++)
        if (whichParameterToUpdate_[par]==1){
          Integer innerLoopCount2=0;
          Double lambdaLocal=lambda;
          parameter_[par]+=lambda*gradP[indPar]/scaling_[par];

          // Now we do not allow non piezo parameters to be more than 2 percent away from the initial guess

//            while((par>=10&&std::abs((parameter_[par]-parameter_initial[par])/parameter_initial[par])>=0.02)
//                  ||(par==7&&std::abs((parameter_[par]-parameter_initial[par])/parameter_initial[par])>=0.01)){

          while(par>=15&&std::abs((parameter_[par]-parameter_initial[par])/parameter_initial[par])>=0.02){                        
            lambdaLocal=0.5*lambdaLocal;

            parameter_[par] = parameter_old[par] + lambdaLocal*gradP[indPar]/scaling_[par];

            innerLoopCount2++;

            if(innerLoopCount2>1000){
              std::cout<<"Help, we are hang up in least-square, penalization of parameter " << par <<std::endl;
              std::cout<<"lambdaLocal = "<< lambdaLocal<<std::endl;
              std::cout<<"(p-p0)/p = " 
                       << std::abs((parameter_[par]-parameter_initial[par])/parameter_initial[par])<<std::endl;
              std::cout<<"reset parameter [" << par << " ] from " << parameter_[par] << " to " << parameter_old[par]<<std::endl;
              parameter_[par]=parameter_old[par];
              break;
            }
            
          }
          indPar++;
        }
      indParC=0;
      for (UInt par=0;par<nrParameter_;par++)
        if (whichParameterToUpdateC_[par]==1){
          parameterC_[par]+=lambda*gradP[indParC+actNrParameter]/scalingC_[par];
          indParC++;
        }
      

      for (UInt par=0;par<nrParameter_;par++)
        if (parameter_[par]<0.0&&par!=6){
          //  parameter[par]=parameter_old[par];
          std::cout<<"parameter( " << par << " ) was negative " <<std::endl;
          negFlag=true;
        }  

      if( imagMaterialParam_ ) {
        updateComplexMaterialData(parameterC_);
      }
      updateMaterialData(parameter_);
      createF(F_hat_, false);
      for (UInt i=0;i<nrMeasuredData;i++)
        F_y[i]=F_hat_[i]-y_hat_[i];
      norm(F_y,normFy1,maxres,y_hat_);


      Integer lineSearchCount=0;

      while (normFy1>normFy0||negFlag==true){
        if (negFlag==true)
          std::cout<< "reduce lambda due to negative parameters or non reduction of norm ... "<<std::endl;
        parameter_=parameter_old;
        parameterC_=parameter_oldC;
        lambda=0.5*lambda;
        std::cout<<"lambda = " << lambda <<std::endl;
        negFlag=false;

        //avoids that programm stagnates!
        lineSearchCount++;
        if (lineSearchCount>1000)
          break;
        
        indPar=0;
        for (UInt par=0;par<nrParameter_;par++)
          if (whichParameterToUpdate_[par]==1){
            Integer innerLoopCount=0;
            Double lambdaLocal=lambda;
            parameter_[par]+=lambda*gradP[indPar]/scaling_[par];

            // Now we do not allow the single parameter to be more than 3 percent away from the initial guess

   //           while((par>=10&&std::abs((parameter_[par]-parameter_initial[par])/parameter_initial[par])>=0.02)   
//                    ||(par==7&&std::abs((parameter_[par]-parameter_initial[par])/parameter_initial[par])>=0.01)){

            while(par>=17 && std::abs((parameter_[par]-parameter_initial[par])/parameter_initial[par])>=0.02){

              //               std::cout<<std::abs((parameter_[par]-parameter_initial[par])/parameter_initial[par])<<std::endl;
//               std::cout<<"parameter["<<par<<"] out of 3 percent circle " <<std::endl;
              lambdaLocal=0.5*lambdaLocal;
              //              std::cout<<"lambdaLocal = "<< lambdaLocal<<std::endl;
              parameter_[par] = parameter_old[par] + lambdaLocal*gradP[indPar]/scaling_[par];
              innerLoopCount++;
              if(innerLoopCount>200){
                std::cout<<"Help, we are hang up in least-square, penalization of parameter " << par <<std::endl;
                std::cout<<"lambdaLocal = "<< lambdaLocal<<std::endl;
                std::cout<<"(p-p0)/p = " 
                         << std::abs((parameter_[par]-parameter_initial[par])/parameter_initial[par])<<std::endl;
                std::cout<<"reset parameter [" << par << " ] from " << parameter_[par] << " to " << parameter_old[par]<<std::endl;
                parameter_[par]=parameter_old[par];
                break;
              }
                             
            }
            indPar++;
          }
        indParC=0;
        for (UInt par=0;par<nrParameter_;par++)
          if (whichParameterToUpdateC_[par]==1){
            parameterC_[par]+=lambda*gradP[indParC+actNrParameter]/scalingC_[par];
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
            //  parameter[par]=parameter_old[par];
            std::cout<<"parameter( " << par << " ) was negative " <<std::endl;
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

