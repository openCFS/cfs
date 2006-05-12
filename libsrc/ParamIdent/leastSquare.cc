


#include "PDE/SinglePDE.hh" 
#include "piezoParamIdent.hh"


namespace CoupledField
{

  
  // xxxxxxxxxxxxxxxxxxxxxxxxxxxxx least square xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  void piezoParamIdent::leastSquare(){ 
    ENTER_FCN("piezoParamIdent::leastSquare");    
    std::cout<<" Entering least square ..."<<std::endl;

    //ptMaterial_=ptMyPDE_->getPDEMaterialData();   // Pointer to MaterialData


    Vector<Double>freqs5;
    if(FALSE){
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
    Vector <Double> gradP(actNrParameter+actNrParameterC);
    Vector<Double> newFreqs;
    parameter_old=parameter_;
    parameter_oldC=parameterC_;
    Double normFy, maxres, normFy0,normFy1;
    Integer indPar=0;
    Integer indParC=0;
    Double lambda=10.0;
    Boolean negFlag=FALSE;      

    //    readInMeasurement(newFreqs);

    std::cout<<"newFreqs:"<<std::endl;
    std::cout<<newFreqs<<std::endl;

    calc_measuredCharge(freqs_, real_, imag_, y_hat_); // out of new measurements

    for (UInt iterIndex=0; iterIndex<1000;iterIndex++){

      updateMaterialData(parameter_);

    if( params->HasValue( "type", "imagMaterialParameter", "materialDataType" ) )
      updateComplexMaterialData(parameterC_);

      createF(F_hat_, FALSE);

      if ((iterIndex+1)%20==0){
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
        std::string filename= "imped.dat";
        impedCurve = new std::ofstream(filename.c_str(),std::basic_ios<char>::out);
        if (!impedCurve){
          std::cerr << "\n ImpedanceCurve.dat could not be initialized" << std::endl;
        }
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
      
      std::cout<<"normF_y0 = " << normFy0<<std::endl;
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

      for (UInt par=0;par<nrParameter_;par++)
        if (whichParameterToUpdate_[par]==1){
          parameter_[par]=1.000001*parameter_[par];
          
          updateMaterialData(parameter_);
          
          createF(F_hat_, FALSE);
          
          for (UInt i=0;i<nrMeasuredData;i++)
            F_y[i]=F_hat_[i]-y_hat_[i];
          
          norm(F_y,normFy,maxres,y_hat_);
          
          gradP[indPar]=normFy0-normFy;
          
          parameter_[par]=parameter_old[par];

          indPar++;

          if (whichParameterToUpdateC_[par]==1){

            parameterC_[par]=1.000001*parameterC_[par];

            updateMaterialData(parameter_);

            if( params->HasValue( "type", "imagMaterialParameter", "materialDataType" ) )
              updateComplexMaterialData(parameterC_);
         
            createF(F_hat_, FALSE);
          
            for (UInt i=0;i<nrMeasuredData;i++)
              F_y[i]=F_hat_[i]-y_hat_[i];
          
            norm(F_y,normFy,maxres,y_hat_);

            gradP[indParC+actNrParameter]=normFy0-normFy;
          
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
          parameter_[par]+=lambda*gradP[indPar]/scaling_[par];
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
          negFlag=TRUE;
        }  

      updateMaterialData(parameter_);
      if( params->HasValue( "type", "imagMaterialParameter", "materialDataType" ) )
        updateComplexMaterialData(parameterC_);
      createF(F_hat_, FALSE);
      for (UInt i=0;i<nrMeasuredData;i++)
        F_y[i]=F_hat_[i]-y_hat_[i];
      norm(F_y,normFy1,maxres,y_hat_);


      while (normFy1>normFy0||negFlag==TRUE){
        if (negFlag==TRUE)
          std::cout<< "reduce lambda due to negative parameters"<<std::endl;
        parameter_=parameter_old;
        parameterC_=parameter_oldC;
        lambda=0.5*lambda;
        std::cout<<"lambda = " << lambda <<std::endl;
        negFlag=FALSE;
        
        indPar=0;
        for (UInt par=0;par<nrParameter_;par++)
          if (whichParameterToUpdate_[par]==1){
            parameter_[par]+=lambda*gradP[indPar]/scaling_[par];
            indPar++;
          }
        indParC=0;
        for (UInt par=0;par<nrParameter_;par++)
          if (whichParameterToUpdateC_[par]==1){
            parameterC_[par]+=lambda*gradP[indParC+actNrParameter]/scalingC_[par];
            indParC++;
          }      

        updateMaterialData(parameter_);
        if( params->HasValue( "type", "imagMaterialParameter", "materialDataType" ) )
          updateComplexMaterialData(parameterC_);
        createF(F_hat_, FALSE);
        for (UInt i=0;i<nrMeasuredData;i++)
          F_y[i]=F_hat_[i]-y_hat_[i];
        norm(F_y,normFy1,maxres,y_hat_);
        
        for (UInt par=0;par<nrParameter_;par++)
          if (parameter_[par]<0.0&&par!=6){
            //  parameter[par]=parameter_old[par];
            std::cout<<"parameter( " << par << " ) was negative " <<std::endl;
            negFlag=TRUE;
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

