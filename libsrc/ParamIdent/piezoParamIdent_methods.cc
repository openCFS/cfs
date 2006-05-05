#include <iostream>
#include <fstream>
#include <string>

#include <PDE/SinglePDE.hh>
#include "piezoParamIdent.hh"
#include "PDE/piezoPDE.hh"
#include "Domain/domain.hh"
#include "Forms/forms_header.hh"



namespace CoupledField
{

  // ========================================================================
  // ========================= piezoParamIdent - Part ===========================
  // ========================================================================

  
  void piezoParamIdent::calc_measuredCharge(Vector<Double> freqs_, 
                                            Vector<Double> & absZ, 
                                            Vector<Double> & phi, 
                                            Vector<Complex> & y_hat){
    ENTER_FCN("piezoParamIdent::calc_measuredCharge");
    Complex Z,j;
    Double x, y;
    j=Complex(0,1);
    Double phase;
    Double randFactor=0.0;

    y_hat_.Resize(nrMeasuredData);

    Vector<Complex> rand(nrMeasuredData);
    for (UInt i=0; i<nrMeasuredData; i++){
      x=absZ[i]*cos(PI/180*phi[i]);
      y=absZ[i]*sin(PI/180*phi[i]);
      Z=Complex(x,y);
      phase = 180.0/PI*std::arg(Z);
    
      y_hat_[i]=sign_*voltage_/(2.0*PI*std::log(Z)*freqs_[i]*j);

    }
   
    if (TRUE){ // generates synthetically created radom noise
      for (UInt i=0;i<nrMeasuredData;i++){
        rand[i] = Complex(2.0*Double(std::rand())/RAND_MAX-1);
        if (randFactor<std::abs(rand[i]))
          randFactor = delta_/std::abs(rand[i]);
      }
      for (UInt i=0;i<nrMeasuredData;i++)
        rand[i]=Complex(randFactor*rand[i])*y_hat_[i];
      if(delta_!=0.0)
        std::cout<<"\n Random noise with data error delta_ = "<< delta_<<std::endl;
      //    std::cout<<rand<<std::endl;
      //       std::cout<<y_hat_<<std::endl;
     
      for (UInt i=0;i<nrMeasuredData;i++)
        y_hat_[i]=y_hat_[i]+rand[i];
    
      Double average_error=0.0;
      for (UInt i=0;i<nrMeasuredData;i++)
        average_error+=std::abs((y_hat_[i]-rand[i])/y_hat_[i]);
      average_error/=nrMeasuredData;
      if(delta_!=0.0){
        std::cout<<"\n The average data error is about ~ " << 
          std::abs(average_error-1.0)*100<<" % " << std::endl;
      }
      
    }

  }// end calc_measuredCharge()


  void piezoParamIdent::calcSyntheticData(Vector<Complex> & y_hat_){
    ENTER_FCN("piezoParamIdent::calcSyntheticData");
    //   std::cout<<"\n We are generating synthetic data, i.e. 
    //   we solve the piezo-equation with exact material - parameters"<<std::endl;
    //   std::cout<<"and alienate the results by alternating +-10Percent"<<std::endl;

//     if(directCoupling_==TRUE)
//       ptMaterial_= ptPDE1_->getPDEMaterialData();   // Pointer to MaterialData
//     else
//       ptMaterial_= ptMyPDE_->getPDEMaterialData();   // Pointer to MaterialData
     
    createF(y_hat_,TRUE); // calculates only forward problems over all omegas_

    for (UInt i=0;i<y_hat_.GetSize();i++){
      if (i%2==0)
        y_hat_[i]*=1.01;
      else
        y_hat_[i]*=0.99;
    }
  }// end calcSyntheticData


  void piezoParamIdent::calcImpedanceCurve(){
    ENTER_FCN("PiezoParamIdent::caclImpedanceCurve");
    
    Boolean reset = TRUE;
     
    updateMaterialData(parameter_);
    updateComplexMaterialData(parameterC_);

    Double maxImpedance=0.0;
    Double minImpedance=1.0e+10;
    Integer freqAtMaxImpedance=0;
    Integer freqAtMinImpedance=0;


    for (UInt fstep = 0; fstep < freqs_.GetSize(); fstep++) { 

          
      ////////////////////////////////////////////////////////
      //                   SOLVES PDE                      //
      ///////////////////////////////////////////////////////  

      //    ptAssemble_->SetReassemble();

      // Set curent frequency value in the mathParser
      domain->GetMathParser()->SetValue( MathParser::GLOB_HANDLER,
                                         "f", actFreq_ );

      ptPDE_->GetSolveStep()->SetActFreq(freqs_[fstep]); 
      ptPDE_->GetSolveStep()->SetActStep(fstep); 
      ptPDE_->GetSolveStep()->PreStepHarmonic(reset); 
      ptPDE_->GetSolveStep()->SolveStepHarmonic(reset);
      ptPDE_->GetSolveStep()->PostStepHarmonic( reset);
      ptPDE_->PostProcess();   
      
      //    ptPDE_->WriteResultsInFile(fstep,freqs_[fstep]);

      Vector<Complex> chargeVec;
      // Vector wich contains charges for each element !
      if(directCoupling_==TRUE)      
        chargeVec = ptPDE1_->getPDE_complexValuedCharge(); 
      else
        chargeVec = ptMyPDE_->getPDE_complexValuedCharge();
      
      Complex charge=Complex(0.0,0.0);
	
      for (UInt i=0;i<chargeVec.GetSize();i++){
       charge+=chargeVec[i];
      }

      if (std::abs(charge)==0){
        std::cout<<"! WARNING: Calculated charge equals zero."<<std::endl;
        std::cout<<" Impedance cannot be calculated"<<std::endl;
        std::cout<<" please check your xml and mesh file for surface charges"<<std::endl;
        std::cout<<" press any key to continue"<<std::endl;
        getchar();
      }
    
      Double imped, phase;
      Complex impedC;
	
      Complex im=Complex(0.0,1);
      impedC=voltage_/(charge*2.0*PI*freqs_[fstep]*im);
      imped = std::abs(voltage_/(charge*2.0*PI*freqs_[fstep]*im)); 
      //    phase = 180.0/PI*(std::arg(charge));
      phase = 180.0/PI*(std::arg(impedC));

      std::cout << fstep <<");\t Frequenz: " << freqs_[fstep] << ";\t Impedanz: "<< imped
                << ";\t Phase: " << phase <<";\t Volt = "<<voltage_<<";\t Ladung = "<< charge<< std::endl;
    
      *impedCurve <<"\n" << freqs_[fstep] << " " << imped << "  " << phase << "  " 
                  << impedC.real()<<"  " << impedC.imag() << "  " << charge.real()
                  << "  " << charge.imag()<< std::endl;

      synMess->setf(std::ios::fixed); // output should by in fixed point numbers,
                                      // i.e. 100.23 instead of 1.0023e+2

      *synMess <<freqs_[fstep]<<"\t"<<imped<<"\t"<<phase<<"\t"<<fstep<<std::endl;

      // determine resonance and antiresonance frequency
      if (imped>maxImpedance){
        maxImpedance=imped;
        antiResonanceFrequency_=freqs_[fstep];
      }

      if (imped<minImpedance){
        minImpedance=imped;
        resonanceFrequency_=freqs_[fstep];
      }

    } //  end loop over freqs

    std::cout<<"++ Resonance lies at " << resonanceFrequency_ 
             <<" Hz with |Z| = " << minImpedance<<std::endl;
    std::cout<<"++ Antiresonance lies at " << antiResonanceFrequency_ 
             <<" Hz with |Z| = " << maxImpedance<<std::endl;

  } // end calcImpedance Curve

  void piezoParamIdent::calcMechDisplCurve(){
    ENTER_FCN("PiezoParamIdent::calcMechDisplCurve");

    
      Boolean reset = TRUE;
      
      if(directCoupling_==TRUE)      
        ptAssemble_ = ptPDE1_->getPDE_assemble(); // Vector wich contains charges for each element !
      else
        ptAssemble_ = ptMyPDE_->getPDE_assemble();
      
      
      updateMaterialData(parameter_);
      updateComplexMaterialData(parameterC_);
      
      // Boolean  adjustDamping = params->IsSet("adjustDamping",  "harmonic");

      for (UInt fstep = 0; fstep < freqs_.GetSize(); fstep++) {       
        
        ////////////////////////////////////////////////////////
        //                   SOLVES PDE                      //
        ///////////////////////////////////////////////////////  
        
        ptPDE_->GetSolveStep()->SetActFreq(freqs_[fstep]); 
        ptPDE_->GetSolveStep()->SetActStep(fstep); 
        ptPDE_->GetSolveStep()->PreStepHarmonic(reset);         
        ptPDE_->GetSolveStep()->SolveStepHarmonic(reset);        
        ptPDE_->GetSolveStep()->PostStepHarmonic(reset);
        ptPDE_->PostProcess();   
        
        //ptMyPDE_->WriteResultsInFile();
        
        BaseNodeStoreSol * ptSol;
        if(directCoupling_==TRUE)      
          ptSol = ptPDE1_->getPDESolution(); // Vector wich contains charges for each element !
        else
          ptSol = ptMyPDE_->getPDESolution();

        NodeStoreSol<Complex> * ptNodeStoreSol;
        ptNodeStoreSol = dynamic_cast<NodeStoreSol<Complex>*>(ptSol);     
        
                //      ptNodeStoreSol->GetGlobalSolVector(ELEC_POTENTIAL, solElecPot_);
        ptNodeStoreSol->GetGlobalSolVector(MECH_DISPLACEMENT, solMechDispl_);
        
        Vector<Complex> nodeSol;
        Vector<Complex> nodeResult(12);
        
        StdVector<std::string> keyVec, attrVec, valVec;
        attrVec = "tag";
        valVec = "paramIdent";
        UInt node1, node2, node3, node4;
        
        keyVec = "paramIdent", "mechDisplAtNode1";
        params->Get( keyVec, attrVec, valVec, node1 );
        keyVec = "paramIdent", "mechDisplAtNode2";
        params->Get( keyVec, attrVec, valVec, node2 );
        keyVec = "paramIdent", "mechDisplAtNode3";
        params->Get( keyVec, attrVec, valVec, node3 );
        keyVec = "paramIdent", "mechDisplAtNode4";
        params->Get( keyVec, attrVec, valVec, node4 );
        
        ptNodeStoreSol->Get(node1,1,nodeResult[0]);	
        ptNodeStoreSol->Get(node1,2,nodeResult[1]);	
        ptNodeStoreSol->Get(node1,3,nodeResult[2]);	
        ptNodeStoreSol->Get(node2,1,nodeResult[3]);	
        ptNodeStoreSol->Get(node2,2,nodeResult[4]);	
        ptNodeStoreSol->Get(node2,3,nodeResult[5]);	
        ptNodeStoreSol->Get(node3,1,nodeResult[6]);	
        ptNodeStoreSol->Get(node3,2,nodeResult[7]);	
        ptNodeStoreSol->Get(node3,3,nodeResult[8]);	
        ptNodeStoreSol->Get(node4,1,nodeResult[9]);	
        ptNodeStoreSol->Get(node4,2,nodeResult[10]);	
        ptNodeStoreSol->Get(node4,3,nodeResult[11]);	

        std::cout<<"Mechanical displacement al selected nodes at frequency "
                 <<freqs_[fstep]<< std::endl;   

        *mechDispl<<freqs_[fstep];
        for(UInt imech=0;imech<nodeResult.GetSize();imech++)
          *mechDispl<<"  "<<nodeResult[imech].real();
        *mechDispl<<std::endl;

        
        
      } //  end loop over freqs
  }


  void piezoParamIdent::createF(Vector<Complex> & F_hat_, Boolean typeOut){
    ENTER_FCN("PiezoParamIdent:createF");
    //   std::cout<<"\nF wil be created ..."<<std::endl;


    ptAssemble_ = ptPDE1_->getPDE_assemble();
    ptAlgsys_ = ptPDE1_->getPDE_algsys();


    F_hat_.Resize(nrMeasuredData);

    Boolean reset = TRUE;
        
    // updateMaterialData(parameter_);
    // updateComplexMaterialData(parameterC_);
    
   
    for (UInt fstep = 0; fstep < nrMeasuredData; fstep++) { // harmonic solver for different frequency - values

          ////////////////////////////////////////////////////////
      //                   SOLVES PDE                      //
      ///////////////////////////////////////////////////////  

      //      ptPDE_->WriteGeneralPDEdefines();   // should not be used, overwrites to much!!    

      //      std::cout<<"\n piezoParam:createF PreStepHarmonic"<<std::endl;

      reset=TRUE;
      ptPDE_->GetSolveStep()->SetActFreq(freqs_[fstep]); 
      ptPDE_->GetSolveStep()->SetActStep(fstep);       
      ptPDE_->GetSolveStep()->PreStepHarmonic(reset); 
         
      //  std::cout<<"\n piezoParam:createF SolveStepHarmonic"<<std::endl;
      ptPDE_->GetSolveStep()->SolveStepHarmonic(reset);
      //std::cout<<"\n after SolveStepHarm " <<std::endl;

      //        std::cout<<"\n piezoParam:createF PostStepHarmonic"<<std::endl;
      ptPDE_->GetSolveStep()->PostStepHarmonic(reset);

      //std::cout<<"\n piezoParam:createF PostProcess at step  "<< fstep << std::endl;
      ptPDE_->PostProcess();


      //////////////////////////////////////////////////////////
      //Retrieves & stores Solution for further calculations  //
      /////////////////////////////////////////////////////////
        
        Vector<Complex> chargeVec;
      
        chargeVec = ptPDE1_->getPDE_complexValuedCharge(); // Vector wich contains charges for each element !
             
        Complex charge=Complex(0.0,0.0);
         
        for (UInt i=0;i<chargeVec.GetSize();i++){
          charge+=chargeVec[i];
        }

        Double x=real_[fstep]*cos(PI/180*imag_[fstep]);
        Double y=real_[fstep]*sin(PI/180*imag_[fstep]);
        Complex Z=Complex(x,y);

        //F_hat_[fstep]=charge;

        // Logarithmic value of F
        F_hat_[fstep]=(sign_*charge*Z)/std::log(Z); // without minus --- classical way ...
        //	F_hat_[fstep]=-(sign_*charge*Z)/std::log(Z); // the "-" became important within the stack


//      } // end loop over elements
        calcAbsImped(charge, freqs_[fstep], fstep,typeOut);   // calculates |Z| and writes results in File

    } // end of loop over all frequencies


    if (typeOut==TRUE){
      for (UInt i=0;i<F_hat_.GetSize();i++)
        std::cout<<"F("<<i<<")="<<F_hat_[i]<<"; \t";
      std::cout<<"\n ------------------------------- " <<std::endl;
      
    }

  
  } // end createF
  // ___________________________________________________________________________________________
  //
  /////////////////////// JacobiMatrix for complex - valued parameter ///////////////////////////////
  // ___________________________________________________________________________________________


  void piezoParamIdent::testJacobiMatrix(Vector<Complex> & F_hat_, 
                                         Matrix<Complex> & JacobiMatrix_, 
                                         Vector<Double> & parameter_,
                                         Vector<Double> & parameterIncrement_, 
                                         Vector<Complex>& solElecPot_,
                                         Vector<Complex> &solMechDispl_){
    ENTER_FCN("piezoParamIdent::testJacobiMatrix");

    Vector<Complex> F_hat__incr(F_hat_.GetSize());
    approxJacobiMatrix_.Resize(JacobiMatrix_.GetSizeRow(), JacobiMatrix_.GetSizeCol());
    Vector<Double> parameter_incr(parameter_.GetSize());
    parameter_incr=parameter_;

    updateMaterialData(parameter_);

    if(directCoupling_==TRUE){
      updateMaterialData(parameter_);
      updateMaterialData(parameter_);
    }
    
    createF(F_hat_, FALSE);

    for (UInt ind_param=0;ind_param<nrParameter_;ind_param++){

      parameter_incr[ind_param]=1.0001*parameter_[ind_param];
      //  std::cout<<parameter_incr[ind_param]<<std::endl;
      updateMaterialData(parameter_incr);

      createF(F_hat__incr,FALSE);

      for (UInt j=0;j<nrMeasuredData;j++)
        approxJacobiMatrix_[j][ind_param]=-(F_hat_[j]-F_hat__incr[j])/((parameter_incr[ind_param]-parameter_[ind_param])*scaling_[ind_param]);

      parameter_incr[ind_param]=parameter_[ind_param];
    }
    //     std::cout<<"\n Here we see the approx. Jacobian and the created Jacobian Matrix:"<<std::endl;
    //     std::cout<<approxJacobiMatrix_<<std::endl;
    //     std::cout<<JacobiMatrix_<<std::endl;
    // getchar();   

  }// end testJacobiMatrix

  void piezoParamIdent::testJacobiMatrix2(Vector<Complex> & F_hat_, Matrix<Complex> & JacobiMatrix_,
                                          Vector<Double> & parameter_,
                                          Vector<Double> & parameterIncrement_, Vector<Complex>& solElecPot_,
                                          Vector<Complex> &solMechDispl_){
    ENTER_FCN("piezoParamIdent::testJacobiMatrix");

    Vector<Complex> F_hat__incr(F_hat_.GetSize());
    Vector<Complex> F_hat__incr2(F_hat_.GetSize());
    Vector<Complex> F_hat__incr3(F_hat_.GetSize());
    Vector<Complex> F_hat__incr4(F_hat_.GetSize());
    approxJacobiMatrix_.Resize(nrMeasuredData,actNrParameter);
    Vector<Double> parameter_incr(nrParameter_);
    Vector<Double> parameter_incr2(nrParameter_);
    Vector<Double> parameter_incr3(nrParameter_);
    Vector<Double> parameter_incr4(nrParameter_);

    parameter_incr=parameter_;
    parameter_incr2=parameter_;
    UInt parInd=0;

    updateMaterialData(parameter_);
    createF(F_hat_, FALSE);

    for (UInt ind_param=0;ind_param<nrParameter_;ind_param++){ 
      if (whichParameterToUpdate_[ind_param]==1){

        parameter_incr[ind_param]=1.00001*parameter_[ind_param];
        //      std::cout<<parameter_incr<<std::endl
        updateMaterialData(parameter_incr);
        createF(F_hat__incr,FALSE);


        parameter_incr2[ind_param]=0.99999*parameter_[ind_param];  
        //      std::cout<<parameter_incr2<<std::endl;
        updateMaterialData(parameter_incr2);
        createF(F_hat__incr2,FALSE);

        // second order FD approximation
        for (UInt j=0;j<nrMeasuredData;j++)
          approxJacobiMatrix_[j][parInd]=0.5*(F_hat__incr[j]-F_hat__incr2[j])/
            ((parameter_incr[ind_param]-parameter_[ind_param])*scaling_[ind_param]);

        parInd++;
        //  std::cout<<"\n Performed second order FD Approx of Jacobian"<<std::endl;
        parameter_incr[ind_param]=parameter_[ind_param];
        parameter_incr2[ind_param]=parameter_[ind_param];
        parameter_incr3[ind_param]=parameter_[ind_param];
        parameter_incr4[ind_param]=parameter_[ind_param];

      }
    }

  }// end testJacobiMatrix


  void piezoParamIdent::testJacobiMatrixC(Vector<Complex> & F_hat_, Matrix<Complex> & JacobiMatrix_, 
                                          Vector<Double> & parameter_){
    ENTER_FCN("piezoParamIdent::testJacobiMatrix");

    Vector<Complex> F_hat__incr(F_hat_.GetSize());
    Vector<Complex> F_hat__incr2(F_hat_.GetSize());
    Vector<Complex> F_hat__incr3(F_hat_.GetSize());
    Vector<Complex> F_hat__incr4(F_hat_.GetSize());
    approxJacobiMatrix_.Resize(nrMeasuredData,actNrParameter+actNrParameterC);
    Vector<Double> parameter_incr(nrParameter_);
    Vector<Double> parameter_incr2(nrParameter_);
    Vector<Double> parameter_incr3(nrParameter_);
    Vector<Double> parameter_incr4(nrParameter_);

    parameter_incr=parameter_;
    parameter_incr2=parameter_;
    UInt parInd=0;

    updateMaterialData(parameter_);
    createF(F_hat_, FALSE);

    for (UInt ind_param=0;ind_param<nrParameter_;ind_param++){ 
      if (whichParameterToUpdate_[ind_param]==1){

        parameter_incr[ind_param]=1.001*parameter_[ind_param];
        //      std::cout<<parameter_incr<<std::endl
        updateMaterialData(parameter_incr);
        createF(F_hat__incr,FALSE);

        parameter_incr2[ind_param]=0.999*parameter_[ind_param];  
        //      std::cout<<parameter_incr2<<std::endl;
        updateMaterialData(parameter_incr2);
        createF(F_hat__incr2,FALSE);

        // second order FD approximation

        for (UInt j=0;j<nrMeasuredData;j++)
          approxJacobiMatrix_[j][parInd]=0.5*(F_hat__incr[j]-F_hat__incr2[j])/
            ((parameter_incr[ind_param]-parameter_[ind_param])*scaling_[ind_param]);

        parInd++;
        //  std::cout<<"\n Performed second order FD Approx of Jacobian"<<std::endl;
        parameter_incr[ind_param]=parameter_[ind_param];
        parameter_incr2[ind_param]=parameter_[ind_param];
        parameter_incr3[ind_param]=parameter_[ind_param];
        parameter_incr4[ind_param]=parameter_[ind_param];

      }
    }

    parInd=0;
    for (UInt ind_param=0;ind_param<nrParameter_;ind_param++){ 
      if (whichParameterToUpdateC_[ind_param]==1){

        parameter_incr[ind_param]=1.001*parameterC_[ind_param];
        //      std::cout<<parameter_incr<<std::endl

        updateComplexMaterialData(parameter_incr);
        createF(F_hat__incr,FALSE);

        parameter_incr2[ind_param]=0.999*parameterC_[ind_param]; 
        //      std::cout<<parameter_incr2<<std::endl;
        updateComplexMaterialData(parameter_incr2);
        createF(F_hat__incr2,FALSE);

        // second order FD approximation
        for (UInt j=0;j<nrMeasuredData;j++)
          approxJacobiMatrix_[j][actNrParameter+parInd]=0.5*(F_hat__incr[j]-F_hat__incr2[j])/
            ((parameter_incr[ind_param]-parameter_[ind_param])*scaling_[ind_param]);

        parInd++;
        //  std::cout<<"\n Performed second order FD Approx of Jacobian"<<std::endl;
        parameter_incr[ind_param]=parameterC_[ind_param];
        parameter_incr2[ind_param]=parameterC_[ind_param];
        parameter_incr3[ind_param]=parameterC_[ind_param];
        parameter_incr4[ind_param]=parameterC_[ind_param];

      }
    }
        
  } // end testJacobiMatrix



  void piezoParamIdent::createAdjointJacobiMatrix(Matrix<Complex> & JacobiMatrix_,
                                                  Matrix<Complex> & adjJacobiMatrix_){
    ENTER_FCN("piezoParamIdent::createAdjointJacobiMatrix");
    //    std::cout<<"\n Adjoint Jacobian will be created ... "<<std::endl;
    adjJacobiMatrix_.Resize(JacobiMatrix_.GetSizeCol(),JacobiMatrix_.GetSizeRow());
    for (UInt i=0;i<JacobiMatrix_.GetSizeCol();i++)
      for (UInt j=0;j<JacobiMatrix_.GetSizeRow();j++){
        //adjJacobiMatrix_[i][j] = JacobiMatrix_[j][i];
        adjJacobiMatrix_[i][j] = std::conj(JacobiMatrix_[j][i]);
        //std::cout<<"F*("<<i<<")("<<j<<")= "<< adjJacobiMatrix_[i][j]<<";\t ";
      }
  } // end createAdjointJacobiMatrix




} // end namespace
