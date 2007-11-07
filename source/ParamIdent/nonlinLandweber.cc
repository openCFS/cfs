// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "PDE/SinglePDE.hh" 
#include "piezoParamIdent.hh"

namespace CoupledField
{

  
  // xxxxxxxxxxxxxxxxxxxxxxxxxxxxx least square xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  void piezoParamIdent::nonlinLandweber(){ 

    //ptMaterial_=ptMyPDE_->getPDEMaterialData();   // Pointer to MaterialData



    calc_measuredCharge(freqs_, real_, imag_, y_hat_); // out of new measurements

    Vector<Complex> act_res(nrMeasuredData);
    Vector<Complex> parUpdate(actNrParameter+actNrParameterC);
    Matrix<Complex> ImgSpaceScaling_Mat(nrMeasuredData,nrMeasuredData);
    Vector<Double> parameterNew(nrParameter_);
    Vector<Double> parameterCNew(nrParameter_);
   //  Vector <Double> parameter_old(10);
//     Vector <Double> parameter_oldC(10);
    Vector<Complex> F_y(nrMeasuredData);
    Matrix<Complex> adjJacobiTemp;
    bool negFlag=false;      

    act_res.Init();
    parUpdate.Init();
    ImgSpaceScaling_Mat.Init();
    parameterNew.Init();
    parameterCNew.Init();
    F_y.Init();
    adjJacobiTemp.Init();

    //    Double omega = 0.015; // geht gut!
    //    Double omega = 0.05;
    // omega für piezo:
    //    Double omega = 0.1;
    // omega für ninas transducer
    Double omega = 0.0005;
    Double normFy, maxres, normFy0,normFy1;
    Integer indPar=0;
    Integer indParC=0;
    parameterNew = parameter_;
    parameterCNew = parameterC_;

    std::cout<<"y_hat"<<std::endl;
    std::cout<<y_hat_<<std::endl;


    for (UInt iterIndex=0; iterIndex<maxNumberNewtonLoops_;iterIndex++){
     
      updateMaterialData(parameter_);      
      if( imagMaterialParam_ ) {
        updateComplexMaterialData(parameterC_);
      }
      
      createF(F_hat_, FALSE);

      std::cout<<"F_hat_"<<std::endl;
      std::cout<<F_hat_<<std::endl;

      act_res = y_hat_-F_hat_;

      norm(act_res,normFy0,maxres,y_hat_);
      
      std::cout<<"||F(p)-y|| = " << normFy0<<std::endl;
      std::cout<< "--------------------------------\n"<<std::endl;

      std::cout<<"\n"<< iterIndex << " - Landweber step " <<std::endl;

      if( imagMaterialParam_ ) 
        testJacobiMatrixC(F_hat_, JacobiMatrix_, parameter_);
      else
        testJacobiMatrix2(F_hat_, JacobiMatrix_, 
                          parameter_, parameterIncrement_,
                          solElecPot_, solMechDispl_);
   
      JacobiMatrix_=approxJacobiMatrix_;

//       std::cout<<"JacobiMatrix_"<<std::endl;
//       std::cout<<JacobiMatrix_<<std::endl;

      for (UInt i=0;i<nrMeasuredData;i++)
        for (UInt j=0;j<nrMeasuredData;j++)
          if (i==j)
            ImgSpaceScaling_Mat[i][j] = 1.0/std::log(std::abs(Complex(real_[i],imag_[i])));
      //ImgSpaceScaling_Mat[i][j] = 1.0/std::log(Complex(real_[i],imag_[i]));
      

      createAdjointJacobiMatrix(JacobiMatrix_,adjJacobiMatrix_);

      adjJacobiTemp = adjJacobiMatrix_*ImgSpaceScaling_Mat;
      adjJacobiMatrix_ = adjJacobiTemp;

      parUpdate=adjJacobiMatrix_*act_res;

//       std::cout<<"parUpdate"<<std::endl;
//       std::cout<<parUpdate<<std::endl;

      // make Landwebers iteration a steepest descent method
       if (true){
         Vector<Complex> adjResidual(actNrParameter+actNrParameterC);
         Vector<Complex> normalResidual(nrMeasuredData);
       
         adjResidual.Init();
         normalResidual.Init();

         adjResidual = adjJacobiMatrix_*act_res;
         normalResidual=JacobiMatrix_*adjResidual;

         Double normAdjResidual = a2norm(adjResidual);

         // relative Norm:
         Double normNormalResidual;

         norm(normalResidual,normNormalResidual,maxres,y_hat_);

//          std::cout<<"adjResidual Norm = " << normAdjResidual<< std::endl;
//          std::cout<<"normalResidual Norm (relative) = " << normNormalResidual<< std::endl;
//          std::cout<<"omega = " << normAdjResidual/normNormalResidual<< std::endl;
         omega = 50*normAdjResidual/normNormalResidual;
         *piezoLog<<omega;
      }

       // make Landweber a minimal error method:
       if (false){

         Vector<Complex> adjResidual(actNrParameter+actNrParameterC);
        
         adjResidual.Init();
         adjResidual = adjJacobiMatrix_*act_res;

         Double normAdjResidual = a2norm(adjResidual);

         std::cout<<"adjResidual Norm = " << normAdjResidual<< std::endl;
         std::cout<<"normFy0= " << normFy0 << std::endl;
         std::cout<<"omega = " << normFy0/normAdjResidual<< std::endl;
         omega =0.0001*normFy0/normAdjResidual;
         *piezoLog<<omega;
      }


      computeScaling();     

      indPar=0;
      for (UInt par=0;par<nrParameter_;par++)
        if (whichParameterToUpdate_[par]==1){
          parameterNew[par]=parameter_[par]+omega*parUpdate[indPar].real()/scaling_[par];
          indPar++;
        }
      indParC=0;
      for (UInt par=0;par<nrParameter_;par++)
        if (whichParameterToUpdateC_[par]==1){
          parameterCNew[par]=parameterC_[par]+1.0e+7*omega*parUpdate[indParC+actNrParameter].real()/scalingC_[par];
          indParC++;
        }

      updateMaterialData(parameterNew);      
      if( imagMaterialParam_ ) {
        updateComplexMaterialData(parameterCNew);
      }

      createF(F_hat_, FALSE);
      act_res = y_hat_-F_hat_;
      norm(act_res,normFy1,maxres,y_hat_);

      std::cout<<"||F(p)-y||_old = " << normFy0<<std::endl;
      std::cout<<"||F(p)-y||_new = " << normFy1<<std::endl;
      std::cout<< "--------------------------------\n"<<std::endl;


      Integer lineSearchCount=0;

      while (normFy1>normFy0||negFlag==true){
        if (negFlag==true)
          std::cout<< "reduce omega due to negative parameters  "<<std::endl;

        omega=0.5*omega;
        std::cout<<"Linesearch - Step = " << lineSearchCount << "\t omega = " << omega << "\t norm = " << normFy1 << std::endl;
        negFlag=false;

        //avoids that programm stagnates!
        lineSearchCount++;
        if (lineSearchCount>20)
          break;
        
        indPar=0;
     
        for (UInt par=0;par<nrParameter_;par++)
          if (whichParameterToUpdate_[par]==1){
            parameterNew[par]=parameter_[par]+omega*parUpdate[indPar].real()/scaling_[par];
            indPar++;
          }
        indParC=0;
        for (UInt par=0;par<nrParameter_;par++)
          if (whichParameterToUpdateC_[par]==1){
            parameterCNew[par]=parameterC_[par]+1.0e+7*omega*parUpdate[indParC+actNrParameter].real()/scalingC_[par];
            indParC++;
          }
        
        
        updateMaterialData(parameterNew);
        
        if( imagMaterialParam_ ) {
          updateComplexMaterialData(parameterCNew);
        }
        
        createF(F_hat_, false);
        for (UInt i=0;i<nrMeasuredData;i++)
          F_y[i]=F_hat_[i]-y_hat_[i];
        
        norm(F_y,normFy1,maxres,y_hat_);
        
        for (UInt par=0;par<nrParameter_;par++)
          if (parameterNew[par]<0.0&&par!=6){
            //  parameter[par]=parameter_old[par];
            std::cout<<"parameter( " << par << " ) was negative " <<std::endl;
            negFlag=true;
          }
      }

      std::cout<<"omega after linesearch = " << omega <<std::endl;
      *piezoLog<<",   "<<omega<<std::endl;


      *parLog<<normFy0<<"  ";
      
      normFy0=normFy1;
      
      parameter_=parameterNew;
      parameterC_=parameterCNew;
      
      for (UInt par=0;par<nrParameter_;par++)
        if (whichParameterToUpdate_[par]==1){
          *parLog<<parameter_[par]<<"  ";
        }
      for (UInt par=0;par<nrParameter_;par++)
        if (whichParameterToUpdateC_[par]==1){
          *parLog<<parameterC_[par]<<"  ";
        }
      
      *parLog<<std::endl;
      
      
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

      // computes impedance curve after each 10 iteration steps
      if ((iterIndex+1)%11==0){
        Integer nrfreqTemp=100;
        Vector<Double> freqsTemp = freqs_;
        freqs_.Resize(nrfreqTemp);
        freqs_.Init();
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


    } // end iter index
    
          
  } // end nonlinlandweber
} // end namespace

