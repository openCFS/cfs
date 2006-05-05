


#include "PDE/SinglePDE.hh" 
#include "piezoParamIdent.hh"


namespace CoupledField
{

  
  // xxxxxxxxxxxxxxxxxxxxxxxxxxxxx least square xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  void piezoParamIdent::leastSquare(){ 
    ENTER_FCN("piezoParamIdent::leastSquare");    
    std::cout<<" Entering least square ..."<<std::endl;

    //ptMaterial_=ptMyPDE_->getPDEMaterialData();   // Pointer to MaterialData

    nrMeasuredData=31;
    Vector<Double> freqs5;
    freqs5.Resize(31);

    if(TRUE){
      freqs5[0]=0.5e+06;
      freqs5[1]=0.7e+06;
      freqs5[2]=1.0e+06;
      freqs5[3]=1.2e+06;
      freqs5[4]=1.3e+06;
      freqs5[5]=1.4e+06;
    freqs5[6]=1.5e+06;
    freqs5[7]=2.5e+06;
    freqs5[8]=2.7e+06;
    freqs5[9]=3.0e+06;
    freqs5[10]=3.2e+06;
    freqs5[11]=3.4e+06;
    freqs5[12]=3.8e+06;
    freqs5[13]=4.1e+06;
    freqs5[14]=4.5e+06;
    freqs5[15]=4.8e+06;
    freqs5[16]=5.0e+06;
    freqs5[17]=5.4e+06;
    freqs5[18]=5.8e+06;
    freqs5[19]=6.1e+06;
    freqs5[20]=6.2e+06;
    freqs5[21]=6.4e+06;
    freqs5[22]=7.5e+06;
    freqs5[23]=8.0e+06;
    freqs5[24]=8.1e+06;
    freqs5[25]=8.2e+06;
    freqs5[26]=8.3e+06;
    freqs5[27]=8.4e+06;
    freqs5[28]=8.5e+06;
    freqs5[29]=8.75e+06;
    freqs5[30]=8.9e+06;
    }
    if(FALSE){

    freqs5[0]=5.1e+06;
    freqs5[1]=5.3e+06;
    freqs5[2]=5.5e+06;
    freqs5[3]=5.7e+06;
    freqs5[4]=5.8e+06;
    freqs5[5]=5.9e+06;
    freqs5[6]=6.0e+06;
    freqs5[7]=6.1e+06;
    freqs5[8]=6.2e+06;
    freqs5[9]=6.3e+06;
    freqs5[10]=6.4e+06;
    freqs5[11]=6.5e+06;
    freqs5[12]=6.6e+06;
    freqs5[13]=6.7e+06;
    freqs5[14]=6.8e+06;
    freqs5[15]=6.9e+06;
    freqs5[16]=7.5e+06;
    freqs5[17]=7.6e+06;
    freqs5[18]=7.7e+06;
    freqs5[19]=7.75e+06;
    freqs5[20]=7.8e+06;
    freqs5[21]=7.85e+06;
    freqs5[22]=7.9e+06;
    freqs5[23]=8.0e+06;
    freqs5[24]=8.1e+06;
    freqs5[25]=8.2e+06;
    freqs5[26]=8.3e+06;
    freqs5[27]=8.4e+06;
    freqs5[28]=8.5e+06;
    freqs5[29]=8.75e+06;
    freqs5[30]=8.9e+06;
    }

    if (FALSE){
      freqs5.Resize(11);
      nrMeasuredData=11;
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
    }
    // 10 V
    if (FALSE){
      freqs5[0]=0.2e+05;
      freqs5[1]=0.35e+05;
      freqs5[2]=0.45e+05;
      freqs5[3]=0.55e+05;
      freqs5[4]=0.6e+05;
      freqs5[5]=0.69e+05;
      freqs5[6]=0.72e+05;
      freqs5[7]=0.75e+05;
      freqs5[8]=0.82e+05;
      freqs5[9]=0.78e+05;
      freqs5[10]=0.82e+05;
      freqs5[11]=0.83e+05;
      freqs5[12]=1.1e+05;
      freqs5[13]=1.2e+05;
      freqs5[14]=1.3e+05;
      freqs5[15]=1.4e+05;
      freqs5[16]=1.5e+05;
      freqs5[17]=1.6e+05;
      freqs5[18]=1.8e+05;
      freqs5[19]=2.0e+05;
    } 
    // 100V
    if (FALSE){
      nrMeasuredData=22;
      freqs5.Resize(22);

      freqs5[0]=0.2e+05;
      freqs5[1]=0.35e+05;
      freqs5[2]=0.45e+05;
      freqs5[3]=0.55e+05;
      freqs5[4]=0.6e+05;
      freqs5[5]=0.69e+05;
      freqs5[6]=0.72e+05;
      freqs5[7]=0.75e+05;
      freqs5[8]=0.77e+05;
      freqs5[9]=0.79e+05;
      freqs5[10]=0.81e+05;
//       freqs5[11]=1.085e+05;
//       freqs5[12]=1.1e+05;
      freqs5[11]=1.11e+05;
      freqs5[12]=1.15e+05;
      freqs5[13]=1.2e+05;
      freqs5[14]=1.3e+05;
      freqs5[15]=1.4e+05;
      freqs5[16]=1.6e+05;
      freqs5[17]=1.8e+05;
      freqs5[18]=1.9e+05;
      freqs5[19]=2.0e+05;
      freqs5[20]=2.2e+05;
      freqs5[21]=2.5e+05;

    }

    // for 50V
    if (FALSE){
      nrMeasuredData=26;
      freqs5.Resize(26);
      freqs5[0]=0.2e+05;
      freqs5[1]=0.3e+05;
      freqs5[2]=0.4e+05;
      freqs5[3]=0.5e+05;
      freqs5[4]=0.6e+05;
      freqs5[5]=0.65e+05;
      freqs5[6]=0.7e+05;
      freqs5[7]=0.75e+05;
      freqs5[8]=0.77e+05;
      freqs5[9]=0.82e+05;
      freqs5[10]=0.85e+05;
      freqs5[11]=1.095e+05;
      freqs5[12]=1.1e+05;
      freqs5[13]=1.15e+05;
      freqs5[14]=1.18e+05;
      freqs5[15]=1.22e+05;
      freqs5[16]=1.28e+05;
      freqs5[17]=1.32e+05;
      freqs5[18]=1.4e+05;
      freqs5[19]=1.5e+05;
      freqs5[20]=1.6e+05;
      freqs5[21]=1.7e+05;
      freqs5[22]=1.8e+05;
      freqs5[23]=2.0e+05;
      freqs5[24]=2.2e+05;
      freqs5[25]=2.4e+05;

    }

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
      
    freqs_=freqs5;


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

    readInMeasurement(newFreqs);

    std::cout<<"newFreqs:"<<std::endl;
    std::cout<<newFreqs<<std::endl;

    calc_measuredCharge(freqs_, real_, imag_, y_hat_); // out of new measurements

    for (UInt iterIndex=0; iterIndex<1000;iterIndex++){

      updateMaterialData(parameter_);

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
      
      Matrix<Double> piezoMat,stiffMat, permMat;
      Matrix<Complex> piezoMatC, stiffMatC, permMatC;

      ptMaterialPiezo_[0]->GetTensor(piezoMat,PIEZO_TENSOR,REAL,FULL);
      ptMaterialMech_[0]->GetTensor(stiffMat,MECH_STIFFNESS_TENSOR, REAL,FULL);
      ptMaterialElec_[0]->GetTensor(permMat,ELEC_PERMITTIVITY,REAL,FULL);

      ptMaterialPiezo_[0]->GetTensor(piezoMatC,PIEZO_TENSOR,IMAG,FULL);
      ptMaterialMech_[0]->GetTensor(stiffMatC,MECH_STIFFNESS_TENSOR,IMAG,FULL);
      ptMaterialElec_[0]->GetTensor(permMatC,ELEC_PERMITTIVITY,IMAG,FULL);

      scaling_[0]=1.0/(stiffMat[0][0]); 
      scaling_[1]=1.0/(stiffMat[2][2]);
      scaling_[2]=1.0/(stiffMat[1][0]);
      scaling_[3]=1.0/(stiffMat[0][2]);
      scaling_[4]=1.0/(stiffMat[3][3]); 
      scaling_[5]=1.0/(piezoMat[1][3]);
      scaling_[6]=std::abs(1.0/((piezoMat)[2][0]));
      scaling_[7]=1.0/((piezoMat)[2][2]);
      scaling_[8]=1.0/((permMat)[0][0]); 
      scaling_[9]=1.0/((permMat)[2][2]);
      
      scalingC_[0]=1.0/(stiffMatC[0][0].imag()); 
      scalingC_[1]=1.0/(stiffMatC[2][2].imag());
      scalingC_[2]=1.0/(stiffMatC[1][0].imag());
      scalingC_[3]=1.0/(stiffMatC[0][2].imag());
      scalingC_[4]=1.0/(stiffMatC[3][3].imag()); 
      scalingC_[5]=1.0/(piezoMatC[1][3].imag());
      scalingC_[6]=std::abs(1.0/((piezoMatC)[2][0].imag()));
      scalingC_[7]=1.0/((piezoMatC)[2][2].imag());
      scalingC_[8]=1.0/((permMatC)[0][0].imag()); 
      scalingC_[9]=1.0/((permMatC)[2][2].imag());
      
          
      
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

