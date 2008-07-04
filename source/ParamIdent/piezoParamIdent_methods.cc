// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <iostream>
#include <fstream>
#include <string>

#include <PDE/SinglePDE.hh>
#include "piezoParamIdent.hh"
#include "Domain/domain.hh"
#include "Driver/stdSolveStep.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "CoupledPDE/PiezoCoupling.hh"
#include "DataInOut/resultHandler.hh"

namespace CoupledField {

  // ========================================================================
  // ========================= piezoParamIdent - Part ===========================
  // ========================================================================


  void piezoParamIdent::evaluateMeasuredData(Vector<Double> freqs_,
      Vector<Double> & absZ, Vector<Double> & phi, Vector<Complex> & y_hat) {
    Complex Z, j;
    Double x, y;
    j=Complex(0, 1);
    Double phase;
    Double randFactor=0.0;

    y_hat_.Resize(nrMeasuredData_);
    y_hat_.Init();

    for (UInt i=0; i<nrMeasuredData_; i++) {
      x=absZ[i]*cos(PI/180*phi[i]);
      y=absZ[i]*sin(PI/180*phi[i]);
      Z=Complex(x, y);
      phase = 180.0/PI * std::atan2(Z.imag(), Z.real());

      if (whichNormCriteria_=="logAmplitude")
      y_hat_[i]=std::log(std::abs(Z)); // or log abs |impedance ...|

      else if (whichNormCriteria_=="logImpedance")
      y_hat_[i]=std::log(Z); // or  imedance ...

      else if (whichNormCriteria_=="amplitude")
      y_hat_[i]=std::abs(Z); // or log impedance ...

      // consider just the phase
      else if (whichNormCriteria_=="phase")
      y_hat_[i]=phase;

      else if (whichNormCriteria_=="zeros")
      y_hat_[i]=0.0;

      Complex u;
      u= Complex(realMech_[i],imagMech_[i]);

      if (whichNormCriteria_=="amplitudeMech")
      y_hat_ [i]=Complex(std::abs(u), 0.0);

      else if (whichNormCriteria_=="logAmplitudeMech")
      y_hat_[i]=Complex(std::log(std::abs(u)), 0.0);

      else if (whichNormCriteria_=="displacementMech")
      y_hat_[i]=u;

      else if (whichNormCriteria_=="phaseMech")
      y_hat_[i]=Complex(180.0/PI * std::atan2(u.imag(), u.real()), 0.0);

    }

    myParam_->Get("artDataNoise", delta_ );

    if (delta_!=0.0) { // generates synthetically created radom noise
      std::cout<<"\n Create random noise with data error delta_ = "<< delta_
      <<std::endl;
      getchar();
      Vector<Complex> rand(nrMeasuredData_);

      for (UInt i=0; i<nrMeasuredData_; i++) {
        rand[i] = Complex(2.0*Double(std::rand())/RAND_MAX-1);
        if (randFactor<std::abs(rand[i]))
        randFactor = delta_/std::abs(rand[i]);
      }
      for (UInt i=0; i<nrMeasuredData_; i++)
      rand[i]=Complex(randFactor*rand[i])*y_hat_[i];
      if (delta_!=0.0)
      std::cout<<"\n Random noise with data error delta_ = "<< delta_
      <<std::endl;

      for (UInt i=0; i<nrMeasuredData_; i++)
      y_hat_[i]=y_hat_[i]+rand[i];

      Double average_error=0.0;
      for (UInt i=0; i<nrMeasuredData_; i++)
      average_error+=std::abs((y_hat_[i]-rand[i])/y_hat_[i]);
      average_error/=nrMeasuredData_;
      if (delta_!=0.0) {
        std::cout<<"\n The average data error is about ~ "
        <<std::abs(average_error-1.0)*100<<" % "<< std::endl;
      }

    }

  }

  void piezoParamIdent::calcImpedanceCurve() {

    if( CalcImpedanceCurve_ == true) {
      std::cout<<"++ Compute impedance curve with "<< freqs_.GetSize()
      <<" steps ... "<<std::endl;
      std::cout<<"++ Results are written in file imped.dat ...\n"<<std::endl;
    }
    if(CalcMechDisplCurve_ == true) {
      std::cout<<"++ Compute mechanical displacement curve with "<< freqs_.GetSize()
      <<" steps ... "<<std::endl;
      std::cout<<"++ Results are written in file mechDispl.dat ...\n"<<std::endl;
    }

    ResultHandler * resHandler = domain->GetResultHandler();
    
    // we only write results during the final computation of the impedance/mechDispl curve
    if (writeResults_==true)
      resHandler->BeginMultiSequenceStep( 1, analysis_, freqs_.GetSize() );

    if (imagMaterialParam_ ) {
      updateComplexMaterialData(parameterC_);
    }
    updateMaterialData(parameter_);

    Double maxImpedance=0.0;
    Double minImpedance=1.0e+10;

    for (UInt fstep = 0; fstep < freqs_.GetSize(); fstep++) {

      ////////////////////////////////////////////////////////
      //                   SOLVES PDE                      //
      ///////////////////////////////////////////////////////

      // Set curent frequency value in the mathParser
      domain->GetMathParser()->SetValue(MathParser::GLOB_HANDLER, "f", freqs_[fstep]);
      domain->GetMathParser()->SetValue(MathParser::GLOB_HANDLER, "step", fstep );

      ptPDE_->GetSolveStep()->SetActFreq(freqs_[fstep]);
      ptPDE_->GetSolveStep()->SetActStep(fstep);
      ptPDE_->GetSolveStep()->PreStepHarmonic();
      ptPDE_->GetSolveStep()->SolveStepHarmonic();
     
      if(writeResults_==true){
        resHandler->BeginStep(fstep+1, freqs_[fstep] );
        ptPDE_->WriteResultsInFile(fstep, freqs_[fstep]);
      }
      resHandler->FinishStep();

      if (CalcImpedanceCurve_==true) {
        
        piezoCpl_->CalcCharges<Complex>( charges_ , chargeNeighborRegion_ );
        Vector<Complex> & chargeVec = charges_->GetVector();
        Complex charge=Complex(0.0, 0.0);

        for (UInt i=0; i<chargeVec.GetSize(); i++) {
          charge+=chargeVec[i];
        }

        if (std::abs(charge)==0) {
          std::cout<<"! WARNING: Calculated charge equals zero."<<std::endl;
          std::cout<<" Impedance cannot be calculated"<<std::endl;
          std::cout<<" please check your xml and mesh file for surface charges"
          <<std::endl;
          std::cout<<" press any key to continue"<<std::endl;
        }

        Double imped, phase;
        Complex impedC;

        Complex im=Complex(0.0, 1);
        impedC=voltage_/(2.0*PI*charge*freqs_[fstep]*im);
        imped = std::abs(voltage_/(2.0*PI*charge*freqs_[fstep]*im));

        phase = 180.0/PI * (std::atan2(impedC.imag(), impedC.real()));

        std::cout << fstep <<"):\t Frequency: "<< freqs_[fstep]
        << ";\t Impedance: "<< imped<< ";\t Phase: "<< phase<< std::endl;

        *impedCurve_ <<"\n"<< freqs_[fstep]<< " "<< imped << "  "<< phase << "  "
        << impedC.real()<<"  "<< impedC.imag() << "  "<< charge.real()<< "  "
        << charge.imag()<< std::endl;


        // determine resonance and antiresonance frequency
        if (imped>maxImpedance) {
          maxImpedance=imped;
          antiResonanceFrequency_=freqs_[fstep];
        }

        if (imped<minImpedance) {
          minImpedance=imped;
          resonanceFrequency_=freqs_[fstep];
        }

      }

      //compute mechanical displacement at nodes specified in xml - Scheme

      if (CalcMechDisplCurve_==true) {

        BaseNodeStoreSol * ptSol;
        ptSol = ptPDE1_->getPDESolution(); // Vector wich contains charges for each element !
        NodeStoreSol<Complex> * ptNodeStoreSol;
        ptNodeStoreSol = dynamic_cast<NodeStoreSol<Complex>*>(ptSol);
        ptNodeStoreSol->GetGlobalSolVector(MECH_DISPLACEMENT, solMechDispl_);
        Vector<Complex> nodeSol;
        Vector<Complex> nodeResult(6);

        UInt node0, node1, node2, node3, node4, node5;
        Vector<UInt> dof(6);
        dof.Init();

        if (myParam_->Has("mechDisplAtNode0")) {
          myParam_->Get("mechDisplAtNode0", node0 );
          myParam_->Get("dofOfMechDispl0", dof[0]);
          ptNodeStoreSol->Get(node0, dof[0], nodeResult[0]);
        }

        if (myParam_->Has("mechDisplAtOppositeNode0")) {
          myParam_->Get("mechDisplAtOppositeNode0", node1 );
          ptNodeStoreSol->Get(node1, dof[0], nodeResult[1]);
        }

        if (myParam_->Has("mechDisplAtNode1")) {
          myParam_->Get("mechDisplAtNode1", node2 );
          myParam_->Get("dofOfMechDispl1", dof[1]);
          ptNodeStoreSol->Get(node2, dof[1], nodeResult[2]);
        }

        if (myParam_->Has("mechDisplAtOppositeNode1")) {
          myParam_->Get("mechDisplAtOppositeNode1", node3 );
          ptNodeStoreSol->Get(node3, dof[1], nodeResult[3]);
        }

        if (myParam_->Has("mechDisplAtNode2")) {
          myParam_->Get("mechDisplAtNode2", node4 );
          myParam_->Get("dofOfMechDispl2", dof[2]);
          ptNodeStoreSol->Get(node4, dof[2], nodeResult[4]);
        }

        if (myParam_->Has("mechDisplAtOppositeNode2")) {
          myParam_->Get("mechDisplAtOppositeNode2", node5 );
          ptNodeStoreSol->Get(node5, dof[2], nodeResult[5]);
        }

        Complex u;

        u=(nodeResult[0]-nodeResult[1]) + (nodeResult[2]-nodeResult[3])
        + (nodeResult[4]-nodeResult[5]);

        std::cout << fstep <<");\t Frequency: "<< freqs_[fstep]
        << ";\t Mechanical displacement" << u << std::endl;

        *mechDispl_<<freqs_[fstep];

        *mechDispl_<<"\t"<< u.real() << "\t"<< u.imag() << "\t"<< std::abs(u)
        << "\t"<< std::atan2(u.imag(), u.real());

        *mechDispl_<<std::endl;
      }

    } //  end loop over freqs

    if( CalcImpedanceCurve_ == true) {
      std::cout<<"\nResonance lies at "<< resonanceFrequency_<<" Hz with |Z| = "
      << minImpedance<<std::endl;
      std::cout<<"Antiresonance lies at "<< antiResonanceFrequency_
      <<" Hz with |Z| = "<< maxImpedance<<std::endl;
      std::cout<<"Piezoelectric coupling k_t ~(fa-fr)/fa = "
      << (antiResonanceFrequency_ - resonanceFrequency_)
      /antiResonanceFrequency_ <<std::endl;
    }

  } // end calcImpedance Curve


  void piezoParamIdent::createF(Vector<Complex> & F_hat_, bool typeOut) {

    // In case that we only consider electrical measurements the solution (F_hat) contains as
    // many entries as we consider frequencies to evaluate.
    // If additionally mechanical measurements are to be considered F_hat is
    // twice the number of measurements

    F_hat_.Resize(nrMeasuredData_);
    F_hat_.Init();

    ResultHandler * resHandler = domain->GetResultHandler();

    for (UInt fstep = 0; fstep < nrMeasuredData_; fstep++) {

      ////////////////////////////////////////////////////////
      //                   SOLVES PDE                      //
      ///////////////////////////////////////////////////////

      domain->GetMathParser()->SetValue(MathParser::GLOB_HANDLER, "f", freqs_[fstep]);
      domain->GetMathParser()->SetValue(MathParser::GLOB_HANDLER, "step", fstep );

      ptPDE_->GetSolveStep()->SetActFreq(freqs_[fstep]);
      ptPDE_->GetSolveStep()->SetActStep(fstep);
      ptPDE_->GetSolveStep()->PreStepHarmonic();
      ptPDE_->GetSolveStep()->SolveStepHarmonic();
      resHandler->FinishStep();
  
      //////////////////////////////////////////////////////////
      //Retrieves & stores Solution for further calculations  //
      /////////////////////////////////////////////////////////
      
      piezoCpl_->CalcCharges<Complex>( charges_ , chargeNeighborRegion_ );
      Vector<Complex> & chargeVec = charges_->GetVector();
      Complex charge=Complex(0.0, 0.0);
      Complex im=Complex(0.0, 1.0);

      for (UInt i=0; i<chargeVec.GetSize(); i++) {
        charge+=chargeVec[i];
      }

      Complex Z=voltage_/(2*PI*charge*freqs_[fstep]*im);

      if (whichNormCriteria_=="logAmplitude") // logarithmic value of impedance
      F_hat_[fstep]=std::log(std::abs(voltage_/(2*PI*charge*freqs_[fstep]*im)));

      else if (whichNormCriteria_=="logImpedance")
      F_hat_[fstep]=std::log(voltage_/(2*PI*charge*freqs_[fstep]*im));

      else if (whichNormCriteria_=="amplitude")
      F_hat_[fstep]=std::abs(voltage_/(2*PI*charge*freqs_[fstep]*im));

      // consider just the phase
      else if (whichNormCriteria_=="phase")
      F_hat_[fstep] = 180.0/PI * std::atan2(Z.imag(), Z.real());

      else if (whichNormCriteria_=="zeros")
      F_hat_[fstep] = 0.0;

      //} // end of loop over all frequencies


      if (whichNormCriteria_=="amplitudeMech"|| whichNormCriteria_
          =="logAmplitudeMech"|| whichNormCriteria_=="phaseMech"
          ||whichNormCriteria_=="displacementMech") {

        //////////////////////////////////////////////////////////
        //Retrieves & stores mechanical solution for further calculations  //
        /////////////////////////////////////////////////////////


        // write mechanical deformation at second part of F_hat:
        BaseNodeStoreSol * ptSol;
        ptSol = ptPDE1_->getPDESolution(); // Vector wich contains charges for each element !
        NodeStoreSol<Complex> * ptNodeStoreSol;
        ptNodeStoreSol = dynamic_cast<NodeStoreSol<Complex>*>(ptSol);
        ptNodeStoreSol->GetGlobalSolVector(MECH_DISPLACEMENT, solMechDispl_);
        Vector<Complex> nodeSol;
        Vector<Complex> nodeResult(6);
        nodeResult.Init();

        UInt node0, node1, node2, node3, node4, node5;
        Vector<UInt> dof(6);
        dof.Init();

        if (myParam_->Has("mechDisplAtNode0")) {
          myParam_->Get("mechDisplAtNode0", node0 );
          myParam_->Get("dofOfMechDispl0", dof[0]);
          ptNodeStoreSol->Get(node0, dof[0], nodeResult[0]);
        }

        if (myParam_->Has("mechDisplAtOppositeNode0")) {
          myParam_->Get("mechDisplAtOppositeNode0", node1 );
          ptNodeStoreSol->Get(node1, dof[0], nodeResult[1]);
        }

        if (myParam_->Has("mechDisplAtNode1")) {
          myParam_->Get("mechDisplAtNode1", node2 );
          myParam_->Get("dofOfMechDispl1", dof[1]);
          ptNodeStoreSol->Get(node2, dof[1], nodeResult[2]);
        }

        if (myParam_->Has("mechDisplAtOppositeNode1")) {
          myParam_->Get("mechDisplAtOppositeNode1", node3 );
          ptNodeStoreSol->Get(node3, dof[1], nodeResult[3]);
        }

        if (myParam_->Has("mechDisplAtNode2")) {
          myParam_->Get("mechDisplAtNode2", node4 );
          myParam_->Get("dofOfMechDispl2", dof[2]);
          ptNodeStoreSol->Get(node4, dof[2], nodeResult[4]);
        }

        if (myParam_->Has("mechDisplAtOppositeNode2")) {
          myParam_->Get("mechDisplAtOppositeNode2", node5 );
          ptNodeStoreSol->Get(node5, dof[2], nodeResult[5]);
        }

        Complex u;

        u=(nodeResult[0]-nodeResult[1]) + (nodeResult[2]-nodeResult[3])
        + (nodeResult[4]-nodeResult[5]);

        if (whichNormCriteria_=="amplitudeMech")
        F_hat_[fstep]=Complex(std::abs(u), 0.0);

        else if (whichNormCriteria_=="logAmplitudeMech")
        F_hat_[fstep]=Complex(std::log(std::abs(u)), 0.0);

        else if (whichNormCriteria_=="displacementMech")
        F_hat_[fstep]=u;

        else if (whichNormCriteria_=="phaseMech")
        F_hat_[fstep]=Complex(180.0/PI * std::atan2(u.imag(), u.real()), 0.0);

        else {
          std::cout
          <<"! your selected fitting quanitity is not supperted -> see documentation"
          <<std::endl;
          getchar();
        }

      }

      if (typeOut==true) {
        for (UInt i=0; i<F_hat_.GetSize(); i++)
        std::cout<<"F("<<i<<")="<<F_hat_[i]<<"; \t";
        std::cout<<"\n ------------------------------- "<<std::endl;
      }

    }

  } // end createF

  void piezoParamIdent::createFVec(Complex & F_hat_, bool typeOut,
      Double frequency) {

    ResultHandler * resHandler = domain->GetResultHandler();

    domain->GetMathParser()->SetValue( MathParser::GLOB_HANDLER,
        "f", UInt(frequency));

    domain->GetMathParser()->SetValue( MathParser::GLOB_HANDLER,
        "step", 0 );

    domain->GetMathParser()->SetValue(MathParser::GLOB_HANDLER, "f", UInt(frequency));
    ptPDE_->GetSolveStep()->SetActFreq(frequency);
    ptPDE_->GetSolveStep()->SetActStep(0);
    ptPDE_->GetSolveStep()->PreStepHarmonic();
    ptPDE_->GetSolveStep()->SolveStepHarmonic();
    ptPDE_->GetSolveStep()->PostStepHarmonic();

    resHandler->FinishStep();

    //////////////////////////////////////////////////////////
    //Retrieves & stores Solution for further calculations  //
    /////////////////////////////////////////////////////////
    
    piezoCpl_->CalcCharges<Complex>( charges_ , chargeNeighborRegion_ );
    Vector<Complex> & chargeVec = charges_->GetVector();
    Complex charge=Complex(0.0, 0.0);
    Complex im=Complex(0.0, 1.0);     
    
    for (UInt i=0;i<chargeVec.GetSize();i++) {
      charge+=chargeVec[i];
    }

    Complex Z=voltage_/(2*PI*charge*frequency*im);

    if (whichNormCriteria_=="logAmplitude") // logarithmic value of impedance
    F_hat_=std::log(std::abs(voltage_/(2*PI*charge*frequency*im)));

    else if (whichNormCriteria_=="logImpedance")
    F_hat_=std::log(voltage_/(2*PI*charge*frequency*im));

    else if (whichNormCriteria_=="amplitude")
    F_hat_=std::abs(voltage_/(2*PI*charge*frequency*im));

    // consider just the phase
    else if (whichNormCriteria_=="phase")
    F_hat_= 180.0/PI * std::atan2(Z.imag(), Z.real());

    else {
      std::cout<<"The fitting quantity "<< whichNormCriteria_ << "is not supported for the computation of confidence intervals" << std::endl;
    }

  } // end createFVec

  void piezoParamIdent::invert(Matrix<Complex> & data) {

    Matrix<Complex> dataTemp;
    dataTemp.Resize(actNrParameter_+actNrParameterC_,actNrParameter_+actNrParameterC_);
    dataTemp=data;

    Vector<Complex> e;
    Vector<Complex> x;
    e.Resize(actNrParameter_+actNrParameterC_);
    e.Init();
    x.Resize(actNrParameter_+actNrParameterC_);
    x.Init();
    Matrix<Complex> inverse;
    inverse.Resize(actNrParameter_+actNrParameterC_, actNrParameter_+actNrParameterC_);
    inverse.Init();

    for (UInt ind=0;ind<data.GetSizeRow();ind++) {
      e[ind]=1.0;

      dataTemp.DirectSolve(x,e);
      dataTemp=data;
      for (UInt indC=0;indC<data.GetSizeCol();indC++)
      inverse[indC][ind]=x[indC];
      x.Resize(actNrParameter_+actNrParameterC_);
      x.Init();
      e[ind]=0.0;
    }
    data=inverse;

  }

  void piezoParamIdent::computeJacobiMatrix() {

    Vector<Complex> F_hat__incr(F_hat_.GetSize());
    Vector<Complex> F_hat__incr2(F_hat_.GetSize());
    Vector<Complex> F_hat__incr3(F_hat_.GetSize());
    Vector<Complex> F_hat__incr4(F_hat_.GetSize());
    approxJacobiMatrix_.Resize(nrMeasuredData_, actNrParameter_);
    approxJacobiMatrix_.Init();
    Vector<Double> parameter_incr(nrParameter_);
    Vector<Double> parameter_incr2(nrParameter_);
    Vector<Double> parameter_incr3(nrParameter_);
    Vector<Double> parameter_incr4(nrParameter_);

    parameter_incr=parameter_;
    parameter_incr2=parameter_;
    UInt parInd=0;

    updateMaterialData(parameter_);
    createF(F_hat_, false);

    for (UInt ind_param=0; ind_param<nrParameter_; ind_param++) {
      if (whichParameterToUpdate_[ind_param]==1) {

        parameter_incr[ind_param]=1.00001*parameter_[ind_param];
        //      std::cout<<parameter_incr<<std::endl
        updateMaterialData(parameter_incr);
        createF(F_hat__incr, false);

        parameter_incr2[ind_param]=0.99999*parameter_[ind_param];
        //      std::cout<<parameter_incr2<<std::endl;
        updateMaterialData(parameter_incr2);
        createF(F_hat__incr2, false);

        // second order FD approximation
        for (UInt j=0; j<nrMeasuredData_; j++)
        approxJacobiMatrix_[j][parInd]=0.5*(F_hat__incr[j]-F_hat__incr2[j])
        /((parameter_incr[ind_param]-parameter_[ind_param])
            *scaling_[ind_param]);

        parInd++;
        //  std::cout<<"\n Performed second order FD Approx of Jacobian"<<std::endl;
        parameter_incr[ind_param]=parameter_[ind_param];
        parameter_incr2[ind_param]=parameter_[ind_param];
        parameter_incr3[ind_param]=parameter_[ind_param];
        parameter_incr4[ind_param]=parameter_[ind_param];

      }
    }

  }// end testJacobiMatrix


  void piezoParamIdent::computeJacobiMatrixC() {

    Vector<Complex> F_hat__incr(F_hat_.GetSize());
    Vector<Complex> F_hat__incr2(F_hat_.GetSize());
    Vector<Complex> F_hat__incr3(F_hat_.GetSize());
    Vector<Complex> F_hat__incr4(F_hat_.GetSize());
    approxJacobiMatrix_.Resize(nrMeasuredData_, actNrParameter_+actNrParameterC_);
    approxJacobiMatrix_.Init();
    Vector<Double> parameter_incr(nrParameter_);
    Vector<Double> parameter_incr2(nrParameter_);
    Vector<Double> parameter_incr3(nrParameter_);
    Vector<Double> parameter_incr4(nrParameter_);

    parameter_incr=parameter_;
    parameter_incr2=parameter_;
    UInt parInd=0;

    updateComplexMaterialData(parameterC_);
    updateMaterialData(parameter_);

    createF(F_hat_, false);

    for (UInt ind_param=0; ind_param<nrParameter_; ind_param++) {
      if (whichParameterToUpdate_[ind_param]==1) {

        parameter_incr[ind_param]=1.001*parameter_[ind_param];
        //      std::cout<<parameter_incr<<std::endl
        updateComplexMaterialData(parameterC_);
        updateMaterialData(parameter_incr);
        createF(F_hat__incr, false);

        parameter_incr2[ind_param]=0.999*parameter_[ind_param];
        //      std::cout<<parameter_incr2<<std::endl;
        updateComplexMaterialData(parameterC_);
        updateMaterialData(parameter_incr2);
        createF(F_hat__incr2, false);

        // second order FD approximation

        for (UInt j=0; j<nrMeasuredData_; j++)
        approxJacobiMatrix_[j][parInd]=0.5*(F_hat__incr[j]-F_hat__incr2[j])
        /((parameter_incr[ind_param]-parameter_[ind_param])
            *scaling_[ind_param]);

        parInd++;
        //  std::cout<<"\n Performed second order FD Approx of Jacobian"<<std::endl;
        parameter_incr[ind_param]=parameter_[ind_param];
        parameter_incr2[ind_param]=parameter_[ind_param];
        parameter_incr3[ind_param]=parameter_[ind_param];
        parameter_incr4[ind_param]=parameter_[ind_param];

      }
    }

    parInd=0;
    for (UInt ind_param=0; ind_param<nrParameter_; ind_param++) {
      if (whichParameterToUpdateC_[ind_param]==1) {

        parameter_incr[ind_param]=1.001*parameterC_[ind_param];
        //      std::cout<<parameter_incr<<std::endl

        updateComplexMaterialData(parameter_incr);
        updateMaterialData(parameter_);
        createF(F_hat__incr, false);

        parameter_incr2[ind_param]=0.999*parameterC_[ind_param];
        //      std::cout<<parameter_incr2<<std::endl;
        updateComplexMaterialData(parameter_incr2);
        updateMaterialData(parameter_);
        createF(F_hat__incr2, false);

        // second order FD approximation
        for (UInt j=0; j<nrMeasuredData_; j++)
        approxJacobiMatrix_[j][actNrParameter_+parInd]=0.5*(F_hat__incr[j]
            -F_hat__incr2[j])
        /((parameter_incr[ind_param]-parameter_[ind_param])
            *scaling_[ind_param]);

        parInd++;
        //  std::cout<<"\n Performed second order FD Approx of Jacobian"<<std::endl;
        parameter_incr[ind_param]=parameterC_[ind_param];
        parameter_incr2[ind_param]=parameterC_[ind_param];
        parameter_incr3[ind_param]=parameterC_[ind_param];
        parameter_incr4[ind_param]=parameterC_[ind_param];

      }
    }

  } // end testJacobiMatrix


  void piezoParamIdent::createAdjointJacobiMatrix() {
    //    std::cout<<"\n Adjoint Jacobian will be created ... "<<std::endl;
    adjJacobiMatrix_.Resize(JacobiMatrix_.GetSizeCol(),
        JacobiMatrix_.GetSizeRow());
    adjJacobiMatrix_.Init();
    for (UInt i=0; i<JacobiMatrix_.GetSizeCol(); i++)
    for (UInt j=0; j<JacobiMatrix_.GetSizeRow(); j++) {
      //adjJacobiMatrix_[i][j] = JacobiMatrix_[j][i];
      adjJacobiMatrix_[i][j] = std::conj(JacobiMatrix_[j][i]);
      //std::cout<<"F*("<<i<<")("<<j<<")= "<< adjJacobiMatrix_[i][j]<<";\t ";
    }
  } // end createAdjointJacobiMatrix


  void piezoParamIdent::computeVecOfSingleDeriv(Vector<Complex> & jacobi, Double omega) {

    Vector<Double>parIncr1, parIncr2;
    parIncr1.Resize(nrParameter_);
    parIncr1.Init();
    parIncr2.Resize(nrParameter_);
    parIncr2.Init();
    Vector<Double>parIncrC1, parIncrC2;
    parIncrC1.Resize(nrParameter_);
    parIncrC1.Init();
    parIncrC2.Resize(nrParameter_);
    parIncrC2.Init();
    parIncr1=parameter_;
    parIncr2=parameter_;
    parIncrC1=parameterC_;
    parIncrC2=parameterC_;
    Complex F_hat_incr1,F_hat_incr2;
    Complex F_hat_incrC1,F_hat_incrC2;
    Integer parInd=0;
    jacobi.Resize(actNrParameter_+actNrParameterC_);
    jacobi.Init();

    for (UInt ind_param=0;ind_param<nrParameter_;ind_param++) {
      if (whichParameterToUpdate_[ind_param]==1) {

        parIncr1[ind_param]=1.1*parameter_[ind_param];
        updateMaterialData(parIncr1);
        createFVec(F_hat_incr1,false,omega);

        parIncr2[ind_param]=0.9*parameter_[ind_param];
        updateMaterialData(parIncr2);
        createFVec(F_hat_incr2,false,omega);

        jacobi[parInd]=0.5*(F_hat_incr1-F_hat_incr2)/
        ((parIncr1[ind_param]-parameter_[ind_param])*scaling_[ind_param]);

        parIncr1[ind_param]=parameter_[ind_param];
        parIncr2[ind_param]=parameter_[ind_param];
        parInd++;

      }
      updateMaterialData(parameter_);
    }
    parInd=0;
    for (UInt ind_param=0;ind_param<nrParameter_;ind_param++) {
      if (whichParameterToUpdateC_[ind_param]==1) {

        parIncrC1[ind_param]=1.000001*parameterC_[ind_param];
        //      std::cout<<parameter_incr<<std::endl
        updateComplexMaterialData(parIncrC1);
        createFVec(F_hat_incrC1,false,omega);

        parIncrC2[ind_param]=0.999999*parameterC_[ind_param];

        updateComplexMaterialData(parIncrC2);
        createFVec(F_hat_incrC2,false,omega);

        jacobi[actNrParameter_+parInd]=0.5*(F_hat_incrC1-F_hat_incrC2)/
        ((parIncrC1[ind_param]-parIncrC2[ind_param])*scalingC_[ind_param]);

        parIncrC1[ind_param]=parameterC_[ind_param];
        parIncrC2[ind_param]=parameterC_[ind_param];
        parInd++;

      }
      updateComplexMaterialData(parameterC_);
    }
  } // end create Jacobian


  void piezoParamIdent::norm(Vector<Complex> & vec, Double & norm,
      Double & norm2, Vector<Complex> & q_meas) {

    maxAndWeightedResNorm(vec, norm2, norm, q_meas); // for real -  valued driver suitable

  } // end norm


  Double piezoParamIdent::calcEuclidianMatrixNorm(Matrix<Complex> & mat) {

    Double norm=0.0;
    for (UInt i=0; i<mat.GetSizeRow(); i++)
    for (UInt j=0; j<mat.GetSizeCol(); j++)
    norm+=std::abs(mat[i][j])*std::abs(mat[i][j]);
    norm=sqrt(norm);
    return norm;

  } // end calcEuclidianMatrixNorm


  void piezoParamIdent::maxAndWeightedResNorm(Vector<Complex> & vec,
      Double & maxNorm, Double & wNorm, Vector<Complex> & q_meas) {
    Double maxNormTemp = 0.0;
    maxNorm=0.0;
    wNorm=0.0;
    Double Denominator=0.0;

    for (UInt i=0; i<vec.GetSize(); i++) {
      maxNormTemp=std::abs(vec[i]);
      Denominator = std::abs(q_meas[i])*std::abs(q_meas[i]);
      wNorm = wNorm+((1.0/Denominator)*std::abs(vec[i])*std::abs(vec[i]));
    }
    wNorm=std::sqrt(wNorm);
  } // end maxAndWeightedNorm


  Double piezoParamIdent::a2norm(Vector<Complex> &vec) {
    Double result=0.0; //Complex(0.0,0.0);
    //      Double real_result;
    for (UInt i=0; i<vec.GetSize(); i++)
    result+=std::abs(vec[i])*std::abs(vec[i]);
    result=sqrt(result);
    return result;
  }

  Double piezoParamIdent::a2norm(Vector<Double> &vec) {
    Double result=0.0; //Complex(0.0,0.0);
    //      Double real_result;
    for (UInt i=0; i<vec.GetSize(); i++)
    result+=std::abs(vec[i])*std::abs(vec[i]);
    result=sqrt(result);
    return result;
  }

  void piezoParamIdent::readMeasuredData() {
    char mDataRow[1024], helpChar[128];
    UInt i=0, j=0, k=0;
    for (UInt ii=0; ii<64; ii++)
    helpChar[ii]=0;

    UInt additionalParameters1=0;
    UInt additionalParametersC1=0;

    UInt additionalParameters2=0;
    UInt additionalParametersC2=0;

    Vector<Double> parameterAddMaterialReal1(2);
    Vector<Double> parameterAddMaterialImag1(2);

    Vector<Double> parameterAddMaterialReal2(2);
    Vector<Double> parameterAddMaterialImag2(2);

    Vector<UInt> whichParameterToUpdateAdd1(2);
    Vector<UInt> whichParameterToUpdateAddC1(2);

    Vector<UInt> whichParameterToUpdateAdd2(2);
    Vector<UInt> whichParameterToUpdateAddC2(2);

    while (allMeasuredData_->getline(mDataRow, 1024)) {
      if (mDataRow[0]=='f') {
        i=2;
        j=0;
        k=0;
        while (mDataRow[i]!='/') {
          if (mDataRow[i]!=',') {
            helpChar[k]=mDataRow[i];
            k++;
            i++;
          } else {
            freqs_[j]=std::atof(helpChar);
            for (UInt l=0; l<=k; l++)
            helpChar[l]=0;
            j++;
            i++;
            k=0;
          }
        }
        nrMeasuredData_ = j;
      }
      else if (mDataRow[0]=='r') {
        i=2;
        k=0;
        j=0;
        while (mDataRow[i]!='/') {
          if (mDataRow[i]!=',') {
            helpChar[k]=mDataRow[i];
            k++;
            i++;
          } else {
            parameter_[j]=std::atof(helpChar);
            for (UInt l=0; l<=k; l++)
            helpChar[l]=0;
            j++;
            i++;
            k=0;
          }
        }
      } else if (mDataRow[0]=='i') {
        i=2;
        k=0;
        j=0;
        while (mDataRow[i]!='/') {
          if (mDataRow[i]!=',') {
            helpChar[k]=mDataRow[i];
            k++;
            i++;
          } else {
            parameterC_[j]=std::atof(helpChar);
            for (UInt l=0; l<=k; l++)
            helpChar[l]=0;
            j++;
            i++;
            k=0;
          }
        }
      }
      else if (mDataRow[0]=='a') {
        i=2;
        k=0;
        j=0;
        while (mDataRow[i]!='/') {
          if (mDataRow[i]!=',') {
            helpChar[k]=mDataRow[i];
            k++;
            i++;
          } else {
            parameterAddMaterialReal1[j]=std::atof(helpChar);
            for (UInt l=0; l<=k; l++)
            helpChar[l]=0;
            j++;
            i++;
            k=0;
          }
        }
        additionalParameters1+=j;
      }
      else if (mDataRow[0]=='b') {
        i=2;
        k=0;
        j=0;
        while (mDataRow[i]!='/') {
          if (mDataRow[i]!=',') {
            helpChar[k]=mDataRow[i];
            k++;
            i++;
          } else {
            parameterAddMaterialImag1[j]=std::atof(helpChar);
            for (UInt l=0; l<=k; l++)
            helpChar[l]=0;
            j++;
            i++;
            k=0;
          }
        }
        additionalParametersC1+=j;
      }
      else if (mDataRow[0]=='c') {
        i=2;
        k=0;
        j=0;
        while (mDataRow[i]!='/') {
          if (mDataRow[i]!=',') {
            helpChar[k]=mDataRow[i];
            k++;
            i++;
          } else {
            parameterAddMaterialReal2[j]=std::atof(helpChar);
            for (UInt l=0; l<=k; l++)
            helpChar[l]=0;
            j++;
            i++;
            k=0;
          }
        }
        additionalParameters2+=j;
      }
      else if (mDataRow[0]=='d') {
        i=2;
        k=0;
        j=0;
        while (mDataRow[i]!='/') {
          if (mDataRow[i]!=',') {
            helpChar[k]=mDataRow[i];
            k++;
            i++;
          } else {
            parameterAddMaterialImag2[j]=std::atof(helpChar);
            for (UInt l=0; l<=k; l++)
            helpChar[l]=0;
            j++;
            i++;
            k=0;
          }
        }
        additionalParametersC2 += j;
      }

      else if (mDataRow[0]=='R') {
        i=2;
        k=0;
        j=0;
        while (mDataRow[i]!='/') {
          if (mDataRow[i]!=',') {
            helpChar[k]=mDataRow[i];
            k++;
            i++;
          } else {
            whichParameterToUpdate_[j]=std::atoi(helpChar);
            for (UInt l=0; l<=k; l++)
            helpChar[l]=0;
            j++;
            i++;
            k=0;
          }
        }
      }
      else if (mDataRow[0]=='A') {
        i=2;
        k=0;
        j=0;
        while (mDataRow[i]!='/') {
          if (mDataRow[i]!=',') {
            helpChar[k]=mDataRow[i];
            k++;
            i++;
          } else {
            whichParameterToUpdateAdd1[j]=std::atoi(helpChar);
            for (UInt l=0; l<=k; l++)
            helpChar[l]=0;
            j++;
            i++;
            k=0;
          }
        }
      }
      else if (mDataRow[0]=='B') {
        i=2;
        k=0;
        j=0;
        while (mDataRow[i]!='/') {
          if (mDataRow[i]!=',') {
            helpChar[k]=mDataRow[i];
            k++;
            i++;
          } else {
            whichParameterToUpdateAddC1[j]=std::atoi(helpChar);
            for (UInt l=0; l<=k; l++)
            helpChar[l]=0;
            j++;
            i++;
            k=0;
          }
        }
      } else if (mDataRow[0]=='I') {
        i=2;
        k=0;
        j=0;
        while (mDataRow[i]!='/') {
          if (mDataRow[i]!=',') {
            helpChar[k]=mDataRow[i];
            k++;
            i++;
          } else {
            whichParameterToUpdateC_[j]=std::atoi(helpChar);
            for (UInt l=0; l<=k; l++)
            helpChar[l]=0;
            j++;
            i++;
            k=0;
          }
        }
      }
      else if (mDataRow[0]=='C') {
        i=2;
        k=0;
        j=0;
        while (mDataRow[i]!='/') {
          if (mDataRow[i]!=',') {
            helpChar[k]=mDataRow[i];
            k++;
            i++;
          } else {
            whichParameterToUpdateAdd2[j]=std::atoi(helpChar);
            for (UInt l=0; l<=k; l++)
            helpChar[l]=0;
            j++;
            i++;
            k=0;
          }
        }
      }
      else if (mDataRow[0]=='D') {
        i=2;
        k=0;
        j=0;
        while (mDataRow[i]!='/') {
          if (mDataRow[i]!=',') {
            helpChar[k]=mDataRow[i];
            k++;
            i++;
          } else {
            whichParameterToUpdateAddC2[j]=std::atoi(helpChar);
            for (UInt l=0; l<=k; l++)
            helpChar[l]=0;
            j++;
            i++;
            k=0;
          }
        }
      }
    }

      // Append additional material parameters (Lame coefficients of adaption layers)
      // to parameter vector

      Vector <Double> parTemp = parameter_;
      Vector <Double> parTempC = parameterC_;
      Vector <UInt> whichParTemp = whichParameterToUpdate_;
      Vector <UInt> whichParTempC = whichParameterToUpdateC_;

      parameter_.Resize(10+additionalParameters1 + additionalParameters2);
      parameterC_.Resize(10+additionalParametersC1 + additionalParametersC2);
      whichParameterToUpdate_.Resize(10+additionalParameters1 + additionalParameters2);
      whichParameterToUpdateC_.Resize(10+additionalParametersC1 + additionalParametersC2);

      for (UInt i=0;i<10;i++) {
        parameter_[i] = parTemp[i];
        parameterC_[i] = parTempC[i];
        whichParameterToUpdate_[i] = whichParTemp[i];
        whichParameterToUpdateC_[i] = whichParTempC[i];
      }

      for (UInt i=0;i<additionalParameters1;i++){
        parameter_[i+10] = parameterAddMaterialReal1[i];
        whichParameterToUpdate_[i+10] = whichParameterToUpdateAdd1[i];
      }
      for (UInt i=0;i<additionalParameters2;i++){
        parameter_[i+10+additionalParameters1] = parameterAddMaterialReal2[i];
        whichParameterToUpdate_[i+10+additionalParameters1] = whichParameterToUpdateAdd2[i];
      }


      for (UInt i=0;i<additionalParametersC1;i++){
        parameterC_[i+10] = parameterAddMaterialImag1[i];
        whichParameterToUpdateC_[i+10] = whichParameterToUpdateAddC1[i];
      }
      for (UInt i=0;i<additionalParametersC2;i++){
        parameterC_[i+10+additionalParametersC1] = parameterAddMaterialImag2[i];
        whichParameterToUpdateC_[i+10+additionalParametersC1] = whichParameterToUpdateAddC2[i];
      }

      std::cout<<"\nParameter vector, real parts:" <<std::endl;
      std::cout<<parameter_<<std::endl;
      std::cout<<"Parameter vector, imaginary parts:" <<std::endl;
      std::cout<<parameterC_<<std::endl;


    } // end read MeasuredData


    void piezoParamIdent::readInMeasurement(Vector<Double> & frequenciesElec) {

      if (whichNormCriteria_=="logAmplitudeMech"|| whichNormCriteria_
          =="amplitudeMech"|| whichNormCriteria_=="phaseMech"||whichNormCriteria_
          =="displacementMech")
      std::cout<<"++ Reading file messMech.dat ... "<<std::endl;
      else
      std::cout<<"++ Reading file mess.dat ... "<<std::endl;

      frequenciesElec.Resize(nrMeasuredData_);
      frequenciesElec.Init();
      real_.Resize(nrMeasuredData_);
      real_.Init();
      imag_.Resize(nrMeasuredData_);
      imag_.Init();

      realMech_.Resize(nrMeasuredData_);
      realMech_.Init();
      imagMech_.Resize(nrMeasuredData_);
      imagMech_.Init();

      Double newFreq, amplitude, phase;
      
      // read only file mess.dat if we have one of the following fitting quantities
      if (whichNormCriteria_=="logImpedance"|| whichNormCriteria_=="logAmplitude"
        || whichNormCriteria_=="amplitude"||whichNormCriteria_=="phase") {

        inMess_.open("mess.dat");
        std::string line;

        while (!inMess_.eof()) {
          std::getline(inMess_, line, '\n');

          inMess_ >> newFreq >> amplitude >> phase;

          for (UInt mInd=0; mInd<nrMeasuredData_; mInd++) {
            if (std::abs(freqs_[mInd]-newFreq)<std::abs(freqs_[mInd]
                                                               -frequenciesElec[mInd])) {
              frequenciesElec[mInd]=newFreq;
              real_[mInd]=amplitude;
              imag_[mInd]=phase;
            }
          } 
        }// end while mess
      }

      // read only file messMech.dat if we have one of the following fitting quantities
      if (whichNormCriteria_=="logAmplitudeMech"|| whichNormCriteria_
          =="amplitudeMech"|| whichNormCriteria_=="phaseMech"||whichNormCriteria_
          =="displacementMech") {

        frequenciesElec.Resize(nrMeasuredData_);
        frequenciesElec.Init();
        realMech_.Resize(nrMeasuredData_);
        realMech_.Init();
        imagMech_.Resize(nrMeasuredData_);
        imagMech_.Init();

        Double newFreq, amplitude, phase;

        inMechMess_.open("messMech.dat");
        std::string line;

        while (!inMechMess_.eof()) {
          std::getline(inMechMess_, line, '\n');

          inMechMess_ >> newFreq >> amplitude >> phase;

          for (UInt mInd=0; mInd<nrMeasuredData_; mInd++) {
            if (std::abs(freqs_[mInd]-newFreq)<std::abs(freqs_[mInd]
                    -frequenciesElec[mInd])) {
              frequenciesElec[mInd]=newFreq;
              realMech_[mInd]=amplitude;
              imagMech_[mInd]=phase;
            }
          }
        }// end while mess


      }// end readInMeasurements
    }

  } // end namespace
