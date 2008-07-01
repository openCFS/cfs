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


    Vector<Complex> F_y(nrMeasuredData_);
    Vector <Double> parameter_old(nrParameter_);
    Vector <Double> parameter_oldC(nrParameter_);
    Vector <Double> parameter_initial(nrParameter_);
    Vector <Double> parameter_initialC(nrParameter_);
    Vector <Double> gradP(actNrParameter_+actNrParameterC_);
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

    evaluateMeasuredData(freqs_, real_, imag_, y_hat_); // out of new measurements

    for (UInt iterIndex=0; iterIndex<maxNumberNewtonLoops_;iterIndex++){
      
      std::cout<<"\nLeast-Squares step " << iterIndex << std::endl;

      if( imagMaterialParam_ ) {
        updateComplexMaterialData(parameterC_);
      }
      updateMaterialData(parameter_);

      // forward simulation, computes F(p) ...
      createF(F_hat_, false);
      
      // std::cout<<F_hat_<<std::endl;
      // std::cout<<y_hat_<<std::endl;

      // compute impedance or mechanical displacement curve after every computeImpedanceCurveAfterStep_ step  
      if ((iterIndex+1)%computeImpedanceCurveAfterStep_==0){
        Integer nrfreqTemp=nrfreq_;
        Vector<Double> freqsTemp = freqs_;
        freqs_.Resize(nrfreqTemp);
        Double startFreqTemp;
        startFreqTemp=startfreq_;
        Double freqincr=(stopfreq_-startfreq_)/nrfreqTemp;
        for(UInt i=0;i<nrfreqTemp;i++){
          startFreqTemp+=freqincr;
          freqs_[i]=startFreqTemp;
        }
        if (impedCurve_)
          impedCurve_->close();
        if(mechDispl_)
          mechDispl_->close();

        std::ostringstream num;
        num<<iterIndex;
        std::string filename= "imped";
        std::string numString = num.str();
        filename.append(".dat");
        impedCurve_ = new std::ofstream(filename.c_str(),std::basic_ios<char>::out);
        filename="mechDispl.dat";
        mechDispl_ = new std::ofstream(filename.c_str(),std::basic_ios<char>::out);

        calcImpedanceCurve();
        freqs_ = freqsTemp;
        nrfreq_=nrfreqTemp;
      }

      
      for (UInt i=0;i<nrMeasuredData_;i++)
        F_y[i]=F_hat_[i]-y_hat_[i];
     
      norm(F_y,normFy0,maxres,y_hat_);
      
      if (normFy0 <= stopRes_){
        std::cout<<"Terminate iteration, since the norm of the residual is smaller than " << stopRes_ <<std::endl;
        break;
      }
      
      std::cout<<"\n||F(p)-y|| = " << normFy0<<std::endl;
      std::cout<< "--------------------------------\n"<<std::endl;

      *parLog_<<normFy0<<"  ";
      for (UInt par=0;par<nrParameter_;par++)
        if (whichParameterToUpdate_[par]==1){
          *parLog_<<parameter_[par]<<"  ";
        }
      for (UInt par=0;par<nrParameter_;par++)
        if (whichParameterToUpdateC_[par]==1){
          *parLog_<<parameterC_[par]<<"  ";
        }

      *parLog_<<std::endl;

      
      indPar=0;
      indParC=0;

      // compute gradient numerically
      for (UInt par=0;par<nrParameter_;par++){

        if (whichParameterToUpdate_[par]==1){

          parameter_[par]=1.00001*parameter_[par];
          
          updateMaterialData(parameter_);
          
          createF(F_hat_, false);
          
          for (UInt i=0;i<nrMeasuredData_;i++)
            F_y[i]=F_hat_[i]-y_hat_[i];
          
          norm(F_y,normFy,maxres,y_hat_);
          
          gradP[indPar]=normFy0-normFy;
          
          parameter_[par]=parameter_old[par];

          indPar++;
        }
      }
      
      // compute part of the gradient for the imaginary parts
      for (UInt par=0;par<nrParameter_;par++){

        if (whichParameterToUpdateC_[par]==1){
          
          parameterC_[par]=1.00001*parameterC_[par];          
         
          if( imagMaterialParam_ ) {
            updateComplexMaterialData(parameterC_);
          }
          updateMaterialData(parameter_);
          
          createF(F_hat_, false);
          
          for (UInt i=0;i<nrMeasuredData_;i++)
            F_y[i]=F_hat_[i]-y_hat_[i];
          
          norm(F_y,normFy,maxres,y_hat_);       
          
          if (whichNormCriteria_=="phaseMech" || whichNormCriteria_=="logAmplitudeMech" || whichNormCriteria_=="displacementMech" )
            gradP[indParC+actNrParameter_]=1.0*(normFy0-normFy);
          else
            gradP[indParC+actNrParameter_]=-1.0*(normFy0-normFy);
          
          parameterC_[par]=parameter_oldC[par];
          
          indParC++;
        }

      }
            
      // we always need proper scaling
      computeScaling();     
               
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
          parameterC_[par]+=lambda*10*gradP[indParC+actNrParameter_]/scalingC_[par];
          indParC++;
        }
      
      if( imagMaterialParam_ ) {
        updateComplexMaterialData(parameterC_);
      }
      updateMaterialData(parameter_);

      createF(F_hat_, false);
      for (UInt i=0;i<nrMeasuredData_;i++)
        F_y[i]=F_hat_[i]-y_hat_[i];
      norm(F_y,normFy1,maxres,y_hat_);

      for (UInt par=0;par<nrParameter_;par++)
        if (parameter_[par]<0.0&&par!=6){
          negFlag=true;
        }
      
      Integer lineSearchCount=0;

      std::cout<<"++ Performing line search ... \n " <<std::endl;
      
      while (normFy1>normFy0||negFlag==true){
  
        parameter_=parameter_old;
        parameterC_=parameter_oldC;
        lambda=0.7*lambda;
        negFlag=false;

        //avoids that programm stagnates!
        lineSearchCount++;
        if (lineSearchCount>20){
          std::cout<<"! TERMINATION of line search after 20 backtracking steps (with stepwidth = " << lambda <<") ... which is a suspicious event ... " <<std::endl;
          std::cout<<"Either we are close to a solution or something went completely wrong." <<std::endl;
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
            parameterC_[par]+=lambda*gradP[indParC+actNrParameter_]/scalingC_[par];
            indParC++;
          }      

        if( imagMaterialParam_ ) 
          updateComplexMaterialData(parameterC_);
        
        updateMaterialData(parameter_);
        createF(F_hat_, false);
        for (UInt i=0;i<nrMeasuredData_;i++)
          F_y[i]=F_hat_[i]-y_hat_[i];
        norm(F_y,normFy1,maxres,y_hat_);
        
        for (UInt par=0;par<nrParameter_;par++)
          if (parameter_[par]<0.0&&par!=6){
            std::cout<<"parameter( " << par << " ) is negative " <<std::endl;
            negFlag=true;
          }
             
      }
      normFy0=normFy1;

      
      std::cout<<"New parameter:"<<std::endl;

      //std::cout<<parameter_<<std::endl;
      for (UInt i=0; i<parameter_.GetSize(); i++){
        std::cout<<i<<") "<< parameter_[i] << " + "<< parameterC_[i] <<"i "<<std::endl;
        *parFinal_<<parameter_[i]<<", ";
      }

      *parFinal_<<std::endl;
      
      for (UInt i=0; i<parameter_.GetSize(); i++)
        *parFinal_<<parameterC_[i]<<", ";

      *parFinal_<<std::endl;
            
      parameter_old=parameter_;
      parameter_oldC=parameterC_;
    } // end iter index
          
  } // end least square

} // end namespace

