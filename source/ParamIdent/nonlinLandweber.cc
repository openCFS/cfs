// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "PDE/SinglePDE.hh" 
#include "piezoParamIdent.hh"

namespace CoupledField {

// xxxxxxxxxxxxxxxxxxxxxxxxxxxxx least square xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

void piezoParamIdent::nonlinLandweber() {

  evaluateMeasuredData(freqs_, real_, imag_, y_hat_); // out of new measurements

  Vector<Complex> act_res(nrMeasuredData_);
  Vector<Complex> parUpdate(actNrParameter_+actNrParameterC_);
  Matrix<Complex> ImgSpaceScaling_Mat(nrMeasuredData_, nrMeasuredData_);
  Vector<Double> parameterNew(nrParameter_);
  Vector<Double> parameterCNew(nrParameter_);
  Vector<Complex> F_y(nrMeasuredData_);
  Matrix<Complex> adjJacobiTemp;
  bool negFlag=false;

  act_res.Init();
  parUpdate.Init();
  ImgSpaceScaling_Mat.Init();
  parameterNew.Init();
  parameterCNew.Init();
  F_y.Init();
  adjJacobiTemp.Init();

  
  Double omega = 0.1; //7.812500e-03;
  Double normFy, maxres, normFy0, normFy1;
  Integer indPar=0;
  Integer indParC=0;
  parameterNew = parameter_;
  parameterCNew = parameterC_;

  
  for (UInt iterIndex=0; iterIndex<maxNumberNewtonLoops_; iterIndex++) {

    updateMaterialData(parameter_);
    if (imagMaterialParam_ ) {
      updateComplexMaterialData(parameterC_);
    }
    
    createF(F_hat_, FALSE);
    
    act_res = y_hat_-F_hat_;
    norm(act_res, normFy0, maxres, y_hat_);
    
    std::cout<<"\n||F(p)-y|| = "<< normFy0<<std::endl;
    std::cout<< "--------------------------------\n"<<std::endl;

    *parLog_<<normFy0<<"  ";

    for (UInt par=0; par<nrParameter_; par++)
      if (whichParameterToUpdate_[par]==1) {
        *parLog_<<parameter_[par]<<"  ";
      }

    for (UInt par=0; par<nrParameter_; par++)
      if (whichParameterToUpdateC_[par]==1) {
        *parLog_<<parameterC_[par]<<"  ";
      }

    *parLog_<<std::endl;

    if (whichMethod_=="Landweber")
     std::cout<<"\n"<< iterIndex << ". Landweber step ... "<<std::endl;
    else if (whichMethod_=="minimalError")
      std::cout<<"\n"<< iterIndex << ". minimal error step ... "<<std::endl;
    else if (whichMethod_=="steepestDescent")
      std::cout<<"\n"<< iterIndex << ". steepest descent step ..." <<std::endl;   
            

    if (imagMaterialParam_ )
      computeJacobiMatrixC();
    else
      computeJacobiMatrix();

    JacobiMatrix_=approxJacobiMatrix_;

    for (UInt i=0; i<nrMeasuredData_; i++)
      for (UInt j=0; j<nrMeasuredData_; j++)
        if (i==j)
          if (whichNormCriteria_=="logAmplitude")
            ImgSpaceScaling_Mat[i][j] = 1.0/std::log(std::abs(Complex(real_[i],
                imag_[i])));
          else if (whichNormCriteria_=="logImpedance")
            ImgSpaceScaling_Mat[i][j] = 1.0/std::log(std::abs(Complex(real_[i],
                imag_[i]))); ///(180.0/PI*std::atan2(imag_[i],real_[i]));
          
          else if (whichNormCriteria_=="amplitude")
            ImgSpaceScaling_Mat[i][j] = 1.0/(std::abs(Complex(real_[i], imag_[i]))); ///(180.0/PI*std::atan2(imag_[i],real_[i]));
       
          else if (whichNormCriteria_=="phase")
            ImgSpaceScaling_Mat[i][j] = 1.0/(180.0/PI*std::atan2(imag_[i],real_[i]));
    
          else if (whichNormCriteria_=="amplitudeMech")
            ImgSpaceScaling_Mat[i][j]=1000.0/(std::abs(Complex(realMech_[i], imagMech_[i])));

          else if (whichNormCriteria_=="logAmplitudeMech")
            ImgSpaceScaling_Mat[i][j]=-1.0/(std::log(std::abs(Complex(realMech_[i],imagMech_[i]))));

          else if (whichNormCriteria_=="displacementMech")
            ImgSpaceScaling_Mat[i][j]=1.0/(0.0001*std::abs((realMech_[i],imagMech_[i])));

          else if (whichNormCriteria_=="phaseMech")
            ImgSpaceScaling_Mat[i][j]=1.0/(180.0/PI * std::atan2(realMech_[i], imagMech_[i]));
          else 
            std::cerr<<"Your choice of the fitting quantity seems to be invalid" <<std::endl;


    createAdjointJacobiMatrix();

    adjJacobiTemp = adjJacobiMatrix_*ImgSpaceScaling_Mat;
    adjJacobiMatrix_ = adjJacobiTemp;

    parUpdate=adjJacobiMatrix_*act_res;
    
    // make Landwebers iteration a steepest descent method
    if (whichMethod_=="steepestDescent") {
      Vector<Complex> adjResidual(actNrParameter_+actNrParameterC_);
      Vector<Complex> normalResidual(nrMeasuredData_);

      adjResidual.Init();
      normalResidual.Init();

      adjResidual = adjJacobiMatrix_*act_res;
      normalResidual=JacobiMatrix_*adjResidual;
      
      Double normAdjResidual = a2norm(adjResidual);

      // relative Norm:
      Double normNormalResidual;
      
     
      if (whichNormCriteria_=="phase")
        normNormalResidual = a2norm(normalResidual);
      else
        norm(normalResidual, normNormalResidual, maxres, y_hat_);
                
      omega = normAdjResidual/normNormalResidual;
                
      if(whichNormCriteria_=="phase")
        omega=0.5*omega;
         
    }

    // make Landweber a minimal error method:
    else if (whichMethod_=="minimalError") {

      Vector<Complex> adjResidual(actNrParameter_+actNrParameterC_);

      adjResidual.Init();
      adjResidual = adjJacobiMatrix_*act_res;

      Double normAdjResidual = a2norm(adjResidual);

      omega =normFy0/normAdjResidual;
    
    }

    computeScaling();
    
    indPar=0;
    for (UInt par=0; par<nrParameter_; par++)
      if (whichParameterToUpdate_[par]==1) {
        parameterNew[par]=parameter_[par]+omega*parUpdate[indPar].real()/scaling_[par];
        indPar++;
      }
    indParC=0;
    for (UInt par=0; par<nrParameter_; par++)
      if (whichParameterToUpdateC_[par]==1) {
        parameterCNew[par]=parameterC_[par]-10.0*omega*parUpdate[indParC+actNrParameter_].real()/scalingC_[par];
        indParC++;
      }

    updateMaterialData(parameterNew);
    if (imagMaterialParam_ ) {
      updateComplexMaterialData(parameterCNew);
    }

    createF(F_hat_, FALSE);
    act_res = y_hat_-F_hat_;
    norm(act_res, normFy1, maxres, y_hat_);
    
    if (normFy0 <= stopRes_){
      std::cout<<"Terminate iteration, since the norm of the residual is smaller than " << stopRes_ <<std::endl;
      break;
    }

    Integer lineSearchCount=0;

    for (UInt par=0; par<nrParameter_; par++)
      if (parameterNew[par]<0.0&&par!=6) {
        std::cout<<"parameter( "<< par << " ) was negative "<<std::endl;
        negFlag=true;
      }

    std::cout<<"++ Performing line search ... " <<std::endl;
    while (normFy1>normFy0||negFlag==true) {
    
      omega=0.8*omega;
      //std::cout<<"Linesearchstep = "<< lineSearchCount << "\t omega "
       //   << omega << "\t ||F(p)-y|| = "<< normFy1 << std::endl;
      negFlag=false;

      //avoids that programm stagnates!
      lineSearchCount++;
      if (lineSearchCount>=30){
        std::cout<<"! TERMINATION of line search after 30 backtracking steps (with step width = " << omega <<") ... which is a suspicious event ... " <<std::endl;
                  std::cout<<"Either we are close to a solution or something went completely wrong." <<std::endl;
                  std::cout<<"Check your results (convergence?) and input data!"<<std::endl;
        break;
      }
        

      indPar=0;

      for (UInt par=0; par<nrParameter_; par++)
        if (whichParameterToUpdate_[par]==1) {
          parameterNew[par]=parameter_[par]+omega*parUpdate[indPar].real()
              /scaling_[par];
          indPar++;
        }
      indParC=0;
      for (UInt par=0; par<nrParameter_; par++)
        if (whichParameterToUpdateC_[par]==1) {
          parameterCNew[par]=parameterC_[par]-10.0*omega
              *parUpdate[indParC+actNrParameter_].real()/scalingC_[par];
          indParC++;
        }

      updateMaterialData(parameterNew);

      if (imagMaterialParam_ ) {
        updateComplexMaterialData(parameterCNew);
      }

      createF(F_hat_, false);
      for (UInt i=0; i<nrMeasuredData_; i++)
        F_y[i]=F_hat_[i]-y_hat_[i];

      norm(F_y, normFy1, maxres, y_hat_);

      for (UInt par=0; par<nrParameter_; par++)
        if (parameterNew[par]<0.0&&par!=6) {
          std::cout<<"parameter( "<< par << " ) was negative "<<std::endl;
          negFlag=true;
        }
    }

    normFy0=normFy1;

    parameter_=parameterNew;
    parameterC_=parameterCNew;

    std::cout<<"\n New parameter:"<<std::endl;

    for (UInt i=0; i<parameter_.GetSize(); i++) {
      std::cout<<i<<") "<< parameter_[i]<< " + "<< parameterC_[i]<<"i "
          <<std::endl;
      *parFinal_<<parameter_[i]<<", ";
    }

    *parFinal_<<std::endl;

    for (UInt i=0; i<parameter_.GetSize(); i++)
      *parFinal_<<parameterC_[i]<<", ";

    *parFinal_<<std::endl;

    // computes impedance curve after each 10 iteration steps
    if ((iterIndex+1)%computeImpedanceCurveAfterStep_==0) {
      
            
      Vector<Double> freqsTemp = freqs_;
      
      freqs_.Resize(nrfreq_);
      freqs_.Init();
      Double startFreqTemp=startfreq_;
      Double freqincr=(stopfreq_-startfreq_)/nrfreq_;
      
      for (UInt i=0; i<nrfreq_; i++) {
        startFreqTemp+=freqincr;
       freqs_[i]=startFreqTemp;
      }
      
      if (impedCurve_)
        impedCurve_->close();
      std::string filename= "imped.dat";
      impedCurve_ = new std::ofstream(filename.c_str(),std::basic_ios<char>::out);
          
      if(mechDispl_)
        mechDispl_->close();
      filename="mechDispl.dat";
      mechDispl_ = new std::ofstream(filename.c_str(),std::basic_ios<char>::out);
    
     calcImpedanceCurve();
     
     // std::cout<<freqsTemp<<std::endl;
     
     freqs_ = freqsTemp;
     
    }

  } // end iter index


} // end nonlinlandweber
} // end namespace

