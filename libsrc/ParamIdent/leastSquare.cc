#include "PDE/SinglePDE.hh" 
#include "piezoParamIdent.hh"


namespace CoupledField
{

  
  // xxxxxxxxxxxxxxxxxxxxxxxxxxxxx NEWTON CG 3 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  void piezoParamIdent::leastSquare(){ 
    ENTER_FCN("piezoParamIdent::leastSquare");    

    MaterialData * ptMaterial=ptMyPDE_->getPDEMaterialData();   // Pointer to MaterialData

    nrMeasuredData=20;
    Vector<Double> freqs5;
    freqs5.Resize(20);

   
    if (TRUE){

      nrMeasuredData=15;
      freqs5.Resize(15);

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
      freqs5[10]=0.85e+05;
      freqs5[11]=0.9e+05;
      freqs5[12]=1.0e+05;
      freqs5[13]=1.1e+05;
      freqs5[14]=1.2e+05;
      freqs5[15]=1.3e+05;
      freqs5[16]=1.4e+05;
      freqs5[17]=1.6e+05;
      freqs5[18]=1.8e+05;
      freqs5[19]=2.0e+05;
    } 
    // 100V
    if (FALSE){
      nrMeasuredData=23;
      freqs5.Resize(23);

      nrMeasuredData=23;
      freqs5.Resize(23);
      freqs5[0]=0.2e+05;
      freqs5[1]=0.35e+05;
      freqs5[2]=0.45e+05;
      freqs5[3]=0.55e+05;
      freqs5[4]=0.6e+05;
      freqs5[5]=0.65e+05;
      freqs5[6]=0.70e+05;
      freqs5[7]=0.72e+05;
      freqs5[9]=0.75e+05;
      freqs5[10]=0.76e+05;
      freqs5[11]=0.78e+05;
      freqs5[12]=0.79e+05;
      freqs5[13]=0.8e+05;
      freqs5[14]=1.09e+05;
      freqs5[15]=1.1e+05;
      freqs5[16]=1.2e+05;
      freqs5[17]=1.3e+05;
      freqs5[18]=1.4e+05;
      freqs5[19]=1.5e+05;
      freqs5[20]=1.6e+05;
      freqs5[21]=1.8e+05;
      freqs5[22]=2.0e+05;
    }


  // 150V
    if (FALSE){
      nrMeasuredData=22;
      freqs5.Resize(22);


      freqs5[0]=0.2e+05;
      freqs5[1]=0.3e+05;
      freqs5[2]=0.35e+05;
      freqs5[3]=0.4e+05;
      freqs5[4]=0.45e+05;
      freqs5[5]=0.5e+05;
      freqs5[6]=0.55e+05;
      freqs5[7]=0.60e+05;
      freqs5[8]=0.65e+05;
      freqs5[9]=0.7e+05;
      freqs5[10]=0.72e+05;
  //     freqs5[11]=0.74e+05;
//       freqs5[12]=0.75e+05;
      //    freqs5[13]=1.1e+05;
      // freqs5[14]=1.15e+05;
      freqs5[11]=1.2e+05;
      freqs5[12]=1.25e+05;
      freqs5[13]=1.3e+05;
      freqs5[14]=1.35e+05;
      freqs5[15]=1.4e+05;
      freqs5[16]=1.45e+05;
      freqs5[17]=1.5e+05;
      freqs5[18]=1.6e+05;
      freqs5[19]=1.7e+05;
      freqs5[20]=1.8e+05;
      freqs5[21]=2.0e+05;

    }
   
    if(FALSE){
      Integer numMeasPoints=30;
      Double fmin=0.2e+5;
      Double fmax=2.0e+5;
      Double fincr=(fmax-fmin)/numMeasPoints;
      freqs5.Resize(numMeasPoints);
      nrMeasuredData=numMeasPoints;
      for(UInt i=0;i<numMeasPoints;i++)
        freqs5[i]=fmin+i*fincr;
    }
      
    freqs=freqs5;


    Vector<Complex> F_y(nrMeasuredData);
    Vector <Double> parameter_old(10);
    Vector <Double> parameter_oldC(10);
    Vector <Double> gradP(actNrParameter+actNrParameterC);
    Vector<Double> newFreqs;
    parameter_old=parameter;
    parameter_oldC=parameterC;
    Double normFy, maxres, normFy0,normFy1;
    Integer indPar=0;
    Integer indParC=0;
    Double lambda=1.0;
    Boolean negFlag=FALSE;      

    //    std::cout<<"Before read in Measuerment"<<std::endl;
  
    readInMeasurement(newFreqs);

    std::cout<<"++ Selected frequencies for fitting:"<<std::endl;
    std::cout<<newFreqs<<std::endl;

    calc_measuredCharge(freqs, real, imag, y_hat); // out of new measurements

    for (UInt iterIndex=0; iterIndex<1000;iterIndex++){

      updateMaterialData(parameter,ptMaterial);
      updateComplexMaterialData(parameterC,ptMaterial);

      createF(ptMaterial, F_hat, FALSE);

      if ((iterIndex+1)%10==0){
        Integer nrfreqTemp=100;
        Vector<Double> freqsTemp = freqs;
        freqs.Resize(nrfreqTemp);
        Double startFreqTemp;
        startFreqTemp=startfreq;
        Double freqincr=(stopfreq-startfreq)/nrfreqTemp;
        for(UInt i=0;i<nrfreqTemp;i++){
          startFreqTemp+=freqincr;
          freqs[i]=startFreqTemp;
        }
        if (impedCurve)
          impedCurve->close();
        std::string filename= "imped.dat";
        impedCurve = new std::ofstream(filename.c_str(),std::basic_ios<char>::out);
        if (!impedCurve){
          std::cerr << "\n ImpedanceCurve.dat could not be initialized" << std::endl;
        }
        calcImpedanceCurve();
        freqs = freqsTemp;
      }
      
      for (UInt i=0;i<nrMeasuredData;i++)
        F_y[i]=F_hat[i]-y_hat[i];
      

      //       std::cout<<"F_hat:"<<std::endl;
      //       std::cout<<F_hat<<std::endl;

      //       std::cout<<"y_hat:"<<std::endl;
      //       std::cout<<y_hat<<std::endl;

      norm(F_y,normFy0,maxres,y_hat);
      
      std::cout<<"Residue (||F(p)-y||) = " << normFy0<<std::endl;
      std::cout<< "--------------------------------\n"<<std::endl;

      *parLog<<normFy0<<"  ";
      for (UInt par=0;par<nrParameter;par++)
        if (whichParameterToUpdate[par]==1){
          *parLog<<parameter[par]<<"  ";
        }
      for (UInt par=0;par<nrParameter;par++)
        if (whichParameterToUpdateC[par]==1){
          *parLog<<parameterC[par]<<"  ";
        }

      *parLog<<std::endl;

      
      indPar=0;
      indParC=0;

      for (UInt par=0;par<nrParameter;par++)
        if (whichParameterToUpdate[par]==1){
          parameter[par]=1.000001*parameter[par];
          
          updateMaterialData(parameter,ptMaterial);
          
          createF(ptMaterial, F_hat, FALSE);
          
          for (UInt i=0;i<nrMeasuredData;i++)
            F_y[i]=F_hat[i]-y_hat[i];
          
          norm(F_y,normFy,maxres,y_hat);
          
          gradP[indPar]=normFy0-normFy;
          
          parameter[par]=parameter_old[par];

          indPar++;

          if (whichParameterToUpdateC[par]==1){

            parameterC[par]=1.000001*parameterC[par];

            updateMaterialData(parameter,ptMaterial);

            updateComplexMaterialData(parameterC,ptMaterial);
         
            createF(ptMaterial, F_hat, FALSE);
          
            for (UInt i=0;i<nrMeasuredData;i++)
              F_y[i]=F_hat[i]-y_hat[i];
          
            norm(F_y,normFy,maxres,y_hat);

            gradP[indParC+actNrParameter]=normFy0-normFy;
          
            parameterC[par]=parameter_oldC[par];

            indParC++;
          
          }
        }
      
      Matrix<Double> *matMat = ptMaterial->GetMatrix();
      scaling[0]=1.0/((*matMat)[0][0]); 
      scaling[1]=1.0/((*matMat)[2][2]);
      scaling[2]=1.0/((*matMat)[1][0]);
      scaling[3]=1.0/((*matMat)[0][2]);
      scaling[4]=1.0/((*matMat)[3][3]); 
      scaling[5]=1.0/((*matMat)[6][4]);
      scaling[6]=std::abs(1.0/((*matMat)[8][0]));
      scaling[7]=1.0/((*matMat)[8][2]);
      scaling[8]=1.0/((*matMat)[6][6]); 
      scaling[9]=1.0/((*matMat)[8][8]);

      Matrix<Double> *matMatC = ptMaterial->GetMatrixC();
      scalingC[0]=1.0/((*matMatC)[0][0]); 
      scalingC[1]=1.0/((*matMatC)[2][2]);
      scalingC[2]=1.0/((*matMatC)[1][0]);
      scalingC[3]=1.0/((*matMatC)[0][2]);
      scalingC[4]=1.0/((*matMatC)[3][3]); 
      scalingC[5]=1.0/((*matMatC)[6][4]);
      scalingC[6]=std::abs(1.0/((*matMatC)[8][0]));
      scalingC[7]=1.0/((*matMatC)[8][2]);
      scalingC[8]=1.0/((*matMatC)[6][6]); 
      scalingC[9]=1.0/((*matMatC)[8][8]);

      //     for (UInt par=0;par<nrParameter;par++)
      //       gradP[par]=gradP[par]/scaling[par];
      
      //      std::cout<<"gradP: "<<std::endl;
      // std::cout<<gradP<<std::endl;

      indPar=0;
      lambda=100*lambda;
      for (UInt par=0;par<nrParameter;par++)
        if (whichParameterToUpdate[par]==1){
          parameter[par]+=lambda*gradP[indPar]/scaling[par];
          indPar++;
        }
      indParC=0;
      for (UInt par=0;par<nrParameter;par++)
        if (whichParameterToUpdateC[par]==1){
          parameterC[par]+=lambda*gradP[indParC+actNrParameter]/scalingC[par];
          indParC++;
        }


      for (UInt par=0;par<nrParameter;par++)
        if (parameter[par]<0.0&&par!=6){
          //  parameter[par]=parameter_old[par];
          std::cout<<"parameter( " << par << " ) is negative " <<std::endl;
          negFlag=TRUE;
        }  

      updateMaterialData(parameter,ptMaterial);
      updateComplexMaterialData(parameterC,ptMaterial);
      createF(ptMaterial, F_hat, FALSE);
      for (UInt i=0;i<nrMeasuredData;i++)
        F_y[i]=F_hat[i]-y_hat[i];
      norm(F_y,normFy1,maxres,y_hat);


      while (normFy1>normFy0||negFlag==TRUE){
        if (negFlag==TRUE)
          std::cout<< "reduce lambda due to a negative parameter or a nondecreasing residue"<<std::endl;
        parameter=parameter_old;
        parameterC=parameter_oldC;
        lambda=0.1*lambda;
        std::cout<<"lambda = " << lambda <<std::endl;
        negFlag=FALSE;
        
        indPar=0;
        for (UInt par=0;par<nrParameter;par++)
          if (whichParameterToUpdate[par]==1){
            parameter[par]+=lambda*gradP[indPar]/scaling[par];
            indPar++;
          }
        indParC=0;
        for (UInt par=0;par<nrParameter;par++)
          if (whichParameterToUpdateC[par]==1){
            parameterC[par]+=lambda*gradP[indParC+actNrParameter]/scalingC[par];
            indParC++;
          }      

        updateMaterialData(parameter,ptMaterial);
        updateComplexMaterialData(parameterC,ptMaterial);
        createF(ptMaterial, F_hat, FALSE);
        for (UInt i=0;i<nrMeasuredData;i++)
          F_y[i]=F_hat[i]-y_hat[i];
        norm(F_y,normFy1,maxres,y_hat);
        
        for (UInt par=0;par<nrParameter;par++)
          if (parameter[par]<0.0&&par!=6){
            //  parameter[par]=parameter_old[par];
            std::cout<<"parameter( " << par << " ) was negative " <<std::endl;
            negFlag=TRUE;
          }
      
      }
      normFy0=normFy1;
      

      std::cout<<"New parameters:"<<std::endl;

      //std::cout<<parameter<<std::endl;
      for (UInt i=0; i<parameter.GetSize(); i++){
        std::cout<<parameter[i] << " + "<< parameterC[i] <<"i"<<std::endl;
        *parFinal<<parameter[i]<<", ";
      }

      *parFinal<<std::endl;
      
      for (UInt i=0; i<parameter.GetSize(); i++)
        *parFinal<<parameterC[i]<<", ";

      *parFinal<<std::endl;
            
      //       std::cout<<"parameter_old:"<<std::endl;
      //       std::cout<<parameter_old<<std::endl;
      parameter_old=parameter;
      parameter_oldC=parameterC;
    } // end iter index
    



          
  } // end least square

} // end namespace
