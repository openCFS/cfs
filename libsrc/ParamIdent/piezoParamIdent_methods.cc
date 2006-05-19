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
    
      // two criteria .. either charge
      if(whichNormCriteria_==1)
        y_hat_[i]=sign_*voltage_/(2.0*PI*std::log(Z)*freqs_[i]*j);
      
      else if(whichNormCriteria_==2)
        y_hat_[i]=std::log(std::abs(Z));      // or log imedance ...



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

    std::cout<<"++ Starting to compute impedance curve with " << freqs_.GetSize()  <<" steps " <<std::endl;
    std::cout<<"++ Results are written in file imped.dat " <<std::endl;
     
    updateMaterialData(parameter_);
    if( params->HasValue( "type", "imagMaterialParameter", "materialDataType" ) )
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
      //    domain->GetMathParser()->SetValue( MathParser::GLOB_HANDLER,
      //                                 "f", actFreq_ );

      ptPDE_->GetSolveStep()->SetActFreq(freqs_[fstep]); 
      ptPDE_->GetSolveStep()->SetActStep(fstep); 
      ptPDE_->GetSolveStep()->PreStepHarmonic(reset); 
      ptPDE_->GetSolveStep()->SolveStepHarmonic(reset);
      ptPDE_->GetSolveStep()->PostStepHarmonic( reset);
      ptPDE_->PostProcess();   
      
      Vector<Complex> chargeVec;
  
      chargeVec = ptPDE1_->getPDE_complexValuedCharge(); 
        
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
      phase = 180.0/PI*(std::arg(impedC));

      std::cout << fstep <<");\t Frequency: " << freqs_[fstep] << ";\t Impedance: "<< imped
                << ";\t Phase: " << phase <<";\t Current = "<<voltage_<<";\t Charge = "<< charge<< std::endl;
    
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
      
      ptAssemble_ = ptPDE1_->getPDE_assemble(); // Vector wich contains charges for each element !
      
      updateMaterialData(parameter_);
      updateComplexMaterialData(parameterC_);
      
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
       
        ptSol = ptPDE1_->getPDESolution(); // Vector wich contains charges for each element !
       

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

    F_hat_.Resize(nrMeasuredData);

    Boolean reset = TRUE;
           
    for (UInt fstep = 0; fstep < nrMeasuredData; fstep++) { 

          ////////////////////////////////////////////////////////
      //                   SOLVES PDE                      //
      ///////////////////////////////////////////////////////  

      ptPDE_->GetSolveStep()->SetActFreq(freqs_[fstep]); 
      ptPDE_->GetSolveStep()->SetActStep(fstep);       
      ptPDE_->GetSolveStep()->PreStepHarmonic(reset); 
      ptPDE_->GetSolveStep()->SolveStepHarmonic(reset);
      ptPDE_->GetSolveStep()->PostStepHarmonic(reset);
      ptPDE_->PostProcess();


      //////////////////////////////////////////////////////////
      //Retrieves & stores Solution for further calculations  //
      /////////////////////////////////////////////////////////
        
        Vector<Complex> chargeVec;
      
        chargeVec = ptPDE1_->getPDE_complexValuedCharge(); // Vector wich contains charges for each element !
             
        Complex charge=Complex(0.0,0.0);
        Complex im=Complex(0.0,1.0);
         
        for (UInt i=0;i<chargeVec.GetSize();i++){
          charge+=chargeVec[i];
        }

        Double x=real_[fstep]*cos(PI/180*imag_[fstep]);
        Double y=real_[fstep]*sin(PI/180*imag_[fstep]);
        Complex Z=Complex(x,y);

        //F_hat_[fstep]=charge;


        // Logarithmic value of charge

        if (whichNormCriteria_==1)
          F_hat_[fstep]=(sign_*charge*Z)/std::log(Z); // without minus --- classical way ...     
        
        else  if (whichNormCriteria_==2)         // logarithmic value of impedance
          F_hat_[fstep]=std::log(std::abs(voltage_/(charge*2.0*PI*freqs_[fstep]*im)));

        //calcAbsImped(charge, freqs_[fstep], fstep,typeOut);   // calculates |Z| and writes results in File

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
        approxJacobiMatrix_[j][ind_param]=
          -(F_hat_[j]-F_hat__incr[j])/
          ((parameter_incr[ind_param]-parameter_[ind_param])*scaling_[ind_param]);

      parameter_incr[ind_param]=parameter_[ind_param];
    }
    //     std::cout<<"\n Here we see the approx. Jacobian 
    //     and the created Jacobian Matrix:"<<std::endl;
    //     std::cout<<approxJacobiMatrix_<<std::endl;
    //     std::cout<<JacobiMatrix_<<std::endl;
    // getchar();   

  }// end testJacobiMatrix

  void piezoParamIdent::testJacobiMatrix2(Vector<Complex> & F_hat_,
                                          Matrix<Complex> & JacobiMatrix_,
                                          Vector<Double> & parameter_,
                                          Vector<Double> & parameterIncrement_, 
                                          Vector<Complex>& solElecPot_,
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




  void piezoParamIdent::calcAbsImped(Complex & charge, Double & freq, UInt & fstep, Boolean typeOut){
    Double imped, phase;
    Complex impedC;

    if (!impedCurve)
      std::cerr<<"Error opening 'ImpedCurve.dat' "<<std::endl;

    Complex im=Complex(0.0,1);
    impedC=voltage_/(charge*2.0*PI*freq*im);
    // We need the following line for a comparison with CAPA
    //    imped = std::abs(voltage_/(charge*2.0*PI*freq*im)); phase = -90 - 180.0/PI*(std::arg(charge));
    // This line makes sense in this routine!

    imped = std::abs(voltage_/(charge*2.0*PI*freq*im));
    phase = 180/PI*(std::arg(impedC));


    if(typeOut==TRUE){
      std::cout<<std::setprecision(10);
      //    std::cout<<"\n Frequency - Impendace - Phase: "<<std::endl;
      std::cout <<"\n Frequency: "<< freq << ", |Z|: "<< std::abs(impedC) << "; Phase: "<< phase << std::endl;
      *impedCurve <<"\n" << freq << "  " << impedC.real()<<"  " << impedC.imag() << imped << "   " << phase << std::endl;
    }



  }  // end calcAbsImped 

  void piezoParamIdent::norm(Vector<Complex> &  vec, Double & norm, Double & norm2,Vector<Complex> & q_meas){
    ENTER_FCN("piezoParamIdent::norm");

    Vector<Complex> y_comp(nrParameter_);
    Vector<Complex> y_temp(nrParameter_);

    switch (whichNorm_){

    case 1:
      norm = a2norm(vec);
      //       std::cout<<"\n l2-Norm = "<<norm<<std::endl;
      break;
    case 2:
      maxAndWeightedResNorm(vec,norm2,norm, q_meas); // for real -  valued driver suitable
      //       std::cout<<"\n weighted-Norm = "<<norm<<std::endl;
      break;
    case 3:
      maxAndEuclNorm(vec,norm,norm2);
      //       std::cout<<"\n max-Norm = "<<norm<<std::endl;
      break;
    case 4:
      //       std::cout<<"\n weighted - logarithmic Norm will be determined ..."<< std::endl;
      for(UInt i=0;i<nrParameter_;i++){
        y_comp[i]=q_meas[i] - vec[i];
        y_comp[i]=std::log(y_comp[i]);
        y_temp[i]=std::log(q_meas[i]);
        vec[i]= std::abs(y_comp[i]-y_temp[i]);
      }
      //       norm=std::sqrt(a2norm(vec));
      maxAndWeightedResNorm(vec,norm2,norm,y_temp);
      // std::cout<<"\n weighted - logarithmic Norm = "<< norm <<std::endl;
      break;

    case 5:
      maxAndWeightedResNorm(vec,norm2,norm, q_meas);  // for complex valued problem
      //       std::cout<<"\n weighted-Norm = "<<norm<<std::endl;
      break;
    case 6:
      maxAndWeightedResNorm(vec,norm2,norm, q_meas); // for real -  valued driver suitable
      //       std::cout<<"\n weighted-Norm = "<<norm<<std::endl;

    case 7:
      logNorm(vec, norm);
      break;

    default:
      norm=a2norm(vec);

    }
  } // end norm
  
  

  Double piezoParamIdent::calcEuclidianMatrixNorm(Matrix<Complex> & mat){
    ENTER_FCN("piezoParamIdent::calcEuclidianMatrixNorm");

    Double norm=0.0;
    for (UInt i=0;i<mat.GetSizeRow();i++)
      for (UInt j=0;j<mat.GetSizeCol();j++)
        norm+=std::abs(mat[i][j])*std::abs(mat[i][j]);
    norm=sqrt(norm);
    return norm;

  } // end calcEuclidianMatrixNorm

  void piezoParamIdent::maxAndEuclNorm(Vector<Complex> & vec,
                                       Double & maxNorm, Double & euclNorm){
    ENTER_FCN("piezoParamIdent::maxAndEuclNorm");
    Double maxNormTemp = 0.0;
    maxNorm=0.0;
    euclNorm=0.0;

    for (UInt i=0;i<vec.GetSize();i++){
      maxNormTemp=std::abs(vec[i]);
      euclNorm += maxNormTemp*maxNormTemp;
      if (maxNormTemp>maxNorm)
        maxNorm=maxNormTemp;
    }
    //    euclNorm=std::sqrt(euclNorm);

  } // end maxAndEuclNorm

  void piezoParamIdent::logNorm(Vector<Complex> & vec, Double & logNorm){
    ENTER_FCN("piezoParamIdent::logNorm");
    logNorm=0.0;
    for (UInt i=0;i<vec.GetSize();i++){
      logNorm = logNorm + std::pow(std::log(std::abs(std::abs(vec[i]))),2);
    }
    //    euclNorm=std::sqrt(euclNorm);
  } // end logNorm



  void piezoParamIdent::maxAndWeightedResNorm(Vector<Complex> & vec,
                                              Double & maxNorm, Double & wNorm,
                                              Vector<Complex> & q_meas){
    ENTER_FCN("piezoParamIdent::maxAndWeightedResNorm");
    Double maxNormTemp = 0.0;
    maxNorm=0.0;
    wNorm=0.0;
    Double Denominator=0.0;

    for (UInt i=0;i<vec.GetSize();i++){
      maxNormTemp=std::abs(vec[i]);
      Denominator = std::abs(q_meas[i])*std::abs(q_meas[i]);
      if (whichNorm_==2){
        //      wNorm = wNorm+((1.0/Denominator)*vec[i]*vec[i]).real(); 
        // this is a good running version!
        wNorm = wNorm+((1.0/Denominator)*std::abs(vec[i])*std::abs(vec[i]));
        //std::cout<<"\n WeightedResNorm = " << std::abs(vec[i])*std::abs(vec[i])<< std::endl;
        //      std::cout<<wNorm<<std::endl;
      }
      else if (whichNorm_==5)
        wNorm = wNorm+((1.0/Denominator)*std::abs(vec[i])*std::abs(vec[i]));
      else if (whichNorm_==6){
        if (whichNewtonCG_!=11){
          std::cout<<"This choice of norm is not valid for your case "<<std::endl;
          std::cout<<"Set Norm in measuredData.dat = 2"<<std::endl;
          getchar();
        }
        Double stepWidth=0.0;
        stepWidth=0.5*std::abs(freqs_[std::min(i+1,freqs_.GetSize())]-
                               freqs_[std::max((Integer)i-1,0)]);
        stepWidth/=1.0e+06;
        //wNorm = wNorm+stepWidth*rhos[i]*((1.0/Denominator)*std::abs(vec[i])*std::abs(vec[i]));
        wNorm = wNorm + 
          omegaDiffVec_[i]*rhos[i]*((1.0/Denominator)*std::abs(vec[i])*std::abs(vec[i]));
        //        wNorm = wNorm+rhos[i]*((1.0/Denominator)*std::abs(vec[i])*std::abs(vec[i]));
        //      getchar();
      }
      
      if (maxNormTemp>maxNorm)
        maxNorm=maxNormTemp;
    }

    //wNorm=std::sqrt(wNorm);
  } // end maxAndWeightedNorm


  void piezoParamIdent::calcNorm2Resid(Vector<Complex> &res,
                                       Double & anorm_, 
                                       UInt nrMeasuredData){
    ENTER_FCN("piezoParamIdent::calcNorm2Resid");
    anorm_=0.0;
    for (UInt i=0;i<2*nrMeasuredData;i++){
      anorm_+=res[i].real()*res[i].real()+ res[i].imag()*res[i].imag();
      anorm_=sqrt(anorm_);
    }
  } // end calcNorm2Resid

  Double piezoParamIdent::norm2Real(Vector<Complex> &vec){
    ENTER_FCN("piezoParamIdent::realA2norm");
    Double result=0.0; 
    //      Double real_result;
    for(UInt i=0;i<vec.GetSize();i++)
      result+=vec[i].real()*vec[i].real();
    result=sqrt(result);
    return result;
  }

  Double piezoParamIdent::realA2norm(Vector<Complex> &vec){
    ENTER_FCN("piezoParamIdent::realA2norm");
    Double result=0.0; 
    Complex resultC = Complex(0.0,0.0);
    //      Double real_result;
    for(UInt i=0;i<vec.GetSize();i++)
      resultC+=vec[i]*vec[i];
    result=resultC.real();
    //    std::cout <<" \n Result in realA2norm = " << result << std::endl;
    //    result=sqrt(result);
    return result;
  }

  Double piezoParamIdent::a2norm(Vector<Complex> &vec){
    ENTER_FCN("piezoParamIdent::a2norm");
    Double result=0.0; //Complex(0.0,0.0);
    //      Double real_result;
    for(UInt i=0;i<vec.GetSize();i++)
      result+=std::abs(vec[i])*std::abs(vec[i]);
    result=sqrt(result);
    return result;
  }

  Double piezoParamIdent::a2norm(Vector<Double> &vec){
    ENTER_FCN("piezoParamIdent::a2norm");
    Double result=0.0; //Complex(0.0,0.0);
    //      Double real_result;
    for(UInt i=0;i<vec.GetSize();i++)
      result+=std::abs(vec[i])*std::abs(vec[i]);
    result=sqrt(result);
    return result;
  }



  void piezoParamIdent::readMeasuredData(Vector<Double> & freqs_, 
                                         Vector<Double> & real_,
                                         Vector<Double> & imag_ ,
                                         Vector<Double> & parameter_, 
                                         Double & voltage_,
                                         UInt & nrMeasuredData, 
                                         Double & thickness_, Double & radius_,
                                         Double & delta_){
    ENTER_FCN( "piezoParamIdent::readMeasuredData" );
    char mDataRow[256], helpChar[64];
    UInt i=0, j=0, k=0;
    while(allMeasuredData->getline(mDataRow, 256)){
      if (mDataRow[0]=='1')
        {i=2;
          while(mDataRow){
            if (mDataRow[i]=='/')
              break;
            if(mDataRow[i]!=','){
              helpChar[k]=mDataRow[i];
              k++; i++;
            }
            else{
              freqs_[j]=atof(helpChar);
              for(UInt l=0;l<=k;l++)
                helpChar[l]=0;
              j++; i++; k=0;
            }
          }
          nrMeasuredData = j;
        }
      else if (mDataRow[0]=='2'){
        i=2;j=0;k=0;
        while(mDataRow){
          if (mDataRow[i]=='/')
            break;
          if(mDataRow[i]!=','){
            helpChar[k]=mDataRow[i];
            k++; i++;
          }
          else{
            real_[j]=atof(helpChar);
            for(UInt l=0;l<=k;l++)
              helpChar[l]=0;
            j++; i++; k=0;
          }
        }
      }
      else if (mDataRow[0]=='3'){
        i=2; k=0; j=0;
        while(mDataRow){
          if (mDataRow[i]=='/')
            break;
          if(mDataRow[i]!=','){
            helpChar[k]=mDataRow[i];
            k++; i++;
          }
          else{
            imag_[j]=atof(helpChar);
            for(UInt l=0;l<=k;l++)
              helpChar[l]=0;
            j++; i++; k=0;
          }
        }
      }
      else if (mDataRow[0]=='4'){
        i=2; k=0; j=0;
        while(mDataRow){
          if (mDataRow[i]=='/')
            break;
          if(mDataRow[i]!=','){
            helpChar[k]=mDataRow[i];
            k++; i++;
          }
          else{
            parameter_[j]=atof(helpChar);
            for(UInt l=0;l<=k;l++)
              helpChar[l]=0;
            j++; i++; k=0;
          }
        }
      }
      else if (mDataRow[0]=='i'){
        i=2; k=0; j=0;
        while(mDataRow){
          if (mDataRow[i]=='/')
            break;
          if(mDataRow[i]!=','){
            helpChar[k]=mDataRow[i];
            k++; i++;
          }
          else{
            parameterC_[j]=atof(helpChar);
            for(UInt l=0;l<=k;l++)
              helpChar[l]=0;
            j++; i++; k=0;
          }
        }
      }
      else if (mDataRow[0]=='P'){
        i=2; k=0; j=0;
        while(mDataRow){
          if (mDataRow[i]=='/')
            break;
          if(mDataRow[i]!=','){
            helpChar[k]=mDataRow[i];
            k++; i++;
          }
          else{
            whichParameterToUpdate_[j]=atoi(helpChar);
            for(UInt l=0;l<=k;l++)
              helpChar[l]=0;
            j++; i++; k=0;
          }
        }
      }
      else if (mDataRow[0]=='Q'){
        i=2; k=0; j=0;
        while(mDataRow){
          if (mDataRow[i]=='/')
            break;
          if(mDataRow[i]!=','){
            helpChar[k]=mDataRow[i];
            k++; i++;
          }
          else{
            whichParameterToUpdateC_[j]=atoi(helpChar);
            for(UInt l=0;l<=k;l++)
              helpChar[l]=0;
            j++; i++; k=0;
          }
        }
      }
      else if (mDataRow[0]=='5'){
        i=2; k=0;
        while(mDataRow){
          if (mDataRow[i]=='/')
            break;
          helpChar[k]=mDataRow[i];
          k++; i++;
        }
        voltage_=atof(helpChar);
        for (UInt l=0;l<=k;l++)
          helpChar[l]=0;
      }
      else if (mDataRow[0]=='6'){
        i=2; k=0;
        while(mDataRow){
          if (mDataRow[i]=='/')
            break;
          helpChar[k]=mDataRow[i];
          k++; i++;
        }
        whichNormCriteria_=atoi(helpChar);
        for (UInt l=0;l<=k;l++)
          helpChar[l]=0;
      }
      else if (mDataRow[0]=='I'){
        i=2; k=0;
        while(mDataRow){
          if (mDataRow[i]=='/')
            break;
          helpChar[k]=mDataRow[i];
          k++; i++;
        }
        CalcImpedanceCurve_=atoi(helpChar);
        for (UInt l=0;l<=k;l++)
          helpChar[l]=0;
      }
     
      else if (mDataRow[0]=='M'){
        i=2; k=0;
        while(mDataRow){
          if (mDataRow[i]=='/')
            break;
          helpChar[k]=mDataRow[i];
          k++; i++;
        }
        whichNewtonCG_=atoi(helpChar); 
        for (UInt l=0;l<=k;l++)
          helpChar[l]=0;
      }
      else if (mDataRow[0]=='N'){
        i=2; k=0;
        while(mDataRow){
          if (mDataRow[i]=='/')
            break;
          helpChar[k]=mDataRow[i];
          k++; i++;
        }
        whichNorm_=atoi(helpChar); 
        for (UInt l=0;l<=k;l++)
          helpChar[l]=0;
      }
      else if (mDataRow[0]=='r'){
        i=2; k=0;
        while(mDataRow){
          if (mDataRow[i]=='/')
            break;
          helpChar[k]=mDataRow[i];
          k++; i++;
        }
        relaxParameter=atof(helpChar); 
        for (UInt l=0;l<=k;l++)
          helpChar[l]=0;
      }
      else if (mDataRow[0]=='7'){
        i=2; k=0;
        while(mDataRow){
          if (mDataRow[i]=='/')
            break;
          helpChar[k]=mDataRow[i];
          k++; i++;
        }
        radius_=atof(helpChar);  
        for (UInt l=0;l<=k;l++)
          helpChar[l]=0;
      }
      else if (mDataRow[0]=='8'){
        i=2; k=0;
        while(mDataRow){
          if (mDataRow[i]=='/')
            break;
          helpChar[k]=mDataRow[i];
          k++; i++;
        }
        delta_=atof(helpChar); // delta - Fehlerniveau
        for (UInt l=0;l<=k;l++)
          helpChar[l]=0;
      }
    }

    if (whichNormCriteria_!=1 && whichNormCriteria_!=2){
      std::cout<<" Check your norm criteria selected in measuredData.dat!!" <<std::endl;
      std::cout<<" It is at point 6) and should be 1 for log(q) and 2 for log(Z)" <<std::endl;
      getchar();
    }
  } // end read MeasuredData


 void piezoParamIdent::readInMeasurement(Vector<Double> & frequencies){

    ENTER_FCN("piezoParamIdent::readInMeasurement");

    std::cout<<"++ Open and read file mess.dat ... " <<std::endl;    

    frequencies.Resize(nrMeasuredData);
    real_.Resize(nrMeasuredData);
    imag_.Resize(nrMeasuredData);

    //! input file, reads set of measurements, frequencies, re(Z), im(z)
    std::ifstream * mess; 
    Char* messChar="mess.dat";
    mess = new std::ifstream(messChar, std::basic_ios<char>::in);

    if (!mess){
      std::cerr << "\n File mess.dat does not exist!" << std::endl;
    }

    std::cout<<"++ Found file mess.dat  ... " <<std::endl;  
      
    char mDataRow[512], helpChar[64];
    UInt i=0, j=0;
    Double newFreq, amplitude, phase;

    while(mess->getline(mDataRow, 512)){
      //      std::cout<<mDataRow<<std::endl;
      i=0;j=0;
      while (mDataRow[i]!='\t'){
        helpChar[j]=mDataRow[i];
        i++;j++;
      }// end while madataRow
      newFreq=atof(helpChar);
      i++;
      j=0;
      for(UInt k=0;k<i;k++) // Delete content of helpChar
        helpChar[k]='0';

      while (mDataRow[i]!='\t'){
        helpChar[j]=mDataRow[i];
        i++;j++;
      }// end while mdataRow
      amplitude=atof(helpChar);
      i++;
      j=0;
      for(UInt k=0;k<i;k++) // Delete content of helpChar
        helpChar[k]='0';
      while (mDataRow[i]!='\n'){
        helpChar[j]=mDataRow[i];
        i++;j++;
      }// end while mdataRow
      phase=atof(helpChar);
    
      for(UInt mInd=0;mInd<nrMeasuredData;mInd++){
        if(std::abs(freqs_[mInd]-newFreq)<std::abs(freqs_[mInd]-frequencies[mInd])){
          frequencies[mInd]=newFreq;
          real_[mInd]=amplitude;
          imag_[mInd]=phase;
        }
      }
    }// end while mess

    freqs_=frequencies;   
    mess->close();
    std::cout<<"++ Open and read of file mess.dat finished ... " <<std::endl;    
  }// end readInMeasurements

  
} // end namespace
