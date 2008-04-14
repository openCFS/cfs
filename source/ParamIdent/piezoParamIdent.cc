// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <iomanip>
#include <PDE/SinglePDE.hh>
#include "Domain/domain.hh"
#include "piezoParamIdent.hh"
#include "CoupledPDE/DirectCoupledPDE.hh"
#include "CoupledPDE/BasePairCoupling.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/resultHandler.hh"
#include "PDE/mechPDE.hh"
#include "Forms/linPressureInt.hh"
#include "Materials/baseMaterial.hh"



namespace CoupledField
{

  // ========================================================================
  // ========================= piezoParamIdent - Part ===========================
  // ========================================================================

  //constructor
  // opens datafiles: measuredData.dat for input, imedCurve.dat and piezoLog.dat for output

  piezoParamIdent :: piezoParamIdent( UInt sequenceStep,
                                      bool isPartOfSequence )
    :SingleDriver( sequenceStep, isPartOfSequence ){


    // Set analysistype
    analysis_ = BasePDE::HARMONIC;
    ptMyPDE_ = NULL;
    residuumParIdent_=1.0;
    resonanceFrequency_=0;
    antiResonanceFrequency_=0;
    tau_=1.0;
    sign_=1.0;
    whichNormCriteria_=1;
    nrMeasuredDataElec_=0;
    nrMeasuredDataMech_=0;
    isPartOfSequence_ = isPartOfSequence;

    // get parameter node
    std::string name = "sequenceStep";
    std::string idx = "index";
    std::string one = "1";
    
    myParam_ = param->Get(name, idx, one )
               ->Get("analysis")->Get("paramIdent");

    // check, if "imagMaterialData" is used
    imagMaterialParam_ = false;
    ParamNode * matNode = 
      param->Get("sequenceStep")->Get("couplingList")->Get("direct")
      ->Get("piezoDirect")->Get("materialDataType", false );

    if( matNode ) {
      imagMaterialParam_ = 
        matNode->Get("type")->AsString() == "imagMaterialParameter";
    }

  }
  
  void piezoParamIdent :: Init() {
    
    // Note: directCoupling_ is always true, as we de no
    // longer have the old structure with a SinglePDE as PiezoPDE
    directCoupling_ = true;

    //  StdVector<std::string> pdeList;
    //     params->GetPDEList( pdeList );
    
    //     numMechMeasurements_=1;
    
    //     if (pdeList[0]=="piezo")
    //       directCoupling_=false;
    //     else 
    //       directCoupling_=true;


    Info->StartProgress( "\n Opening in and output files ... " );

    std::string  filenameMeasuredData="measuredData.dat";
    allMeasuredData = new std::ifstream(filenameMeasuredData.c_str(), 
                                        std::basic_ios<char>::in);
    if (!allMeasuredData){
      std::cerr << "\n File measuredData.dat does not exist!" << std::endl;
    }

//     filenameMeasuredData="mess.dat";
//     mess = new std::ifstream(filenameMeasuredData.c_str(), std::basic_ios<char>::in);
//     if (!mess){
//       std::cerr << "\n File mess.dat does not exist!" << std::endl;
//     }

//     filenameMeasuredData="messMech.dat";
//     messMech = new std::ifstream(filenameMeasuredData.c_str(), std::basic_ios<char>::in);
//     if (!messMech){
//       std::cerr << "\n File mess.dat does not exist!" << std::endl;
//     }

    std::string filename= "imped.dat";
    impedCurve = new std::ofstream(filename.c_str(),std::basic_ios<char>::out);
    if (!impedCurve){
      std::cerr << "\n ImpedanceCurve.dat could not be initialized" << std::endl;
    }

    std::string filenameSynMess= "synMess.dat";
    synMess = new std::ofstream(filenameSynMess.c_str(),std::basic_ios<char>::out);
    if (!synMess){
      std::cerr << "\n syMess.dat could not be initialized" << std::endl;
    }

    std::string filenameLog= "piezoLog.dat";
    piezoLog = new std::ofstream(filenameLog.c_str(),std::basic_ios<char>::out);
    if (!piezoLog){
      std::cerr << "\n piezoLog.dat could not be initialized" << std::endl;
    }

    std::string filenameMechDispl= "mechDispl.dat";
    mechDispl = new std::ofstream(filenameMechDispl.c_str(),std::basic_ios<char>::out);
    if (!mechDispl){
      std::cerr << "\n mechDispl.dat could not be initialized" << std::endl;
    }

    std::string filenameParLog= "parLog.dat";
    parLog = new std::ofstream(filenameParLog.c_str(),std::basic_ios<char>::out);
    if (!parLog){
      std::cerr << "\n parLog.dat could not be initialized" << std::endl;
    }

    std::string filenameParFinal= "parFinal.dat";
    parFinal = new std::ofstream(filenameParFinal.c_str(),std::basic_ios<char>::out);
    if (!parFinal){
      std::cerr << "\n parFinal.dat could not be initialized" << std::endl;
    }

    std::string filenameAllTensors = "allTensors.dat";
    allTensors = new std::ofstream(filenameAllTensors.c_str(),std::basic_ios<char>::out);
    if (!allTensors){
      std::cerr << "\n allTensors.dat could not be initialized" << std::endl;
    }

    std::string filenameOptimalFreqs= "optimalFreqs.dat";
    optimalFreqs = new std::ofstream(filenameOptimalFreqs.c_str(),std::basic_ios<char>::out);
    if (!optimalFreqs){
      std::cerr << "\n optimalFreqs.dat could not be initialized" << std::endl;
    }

    std::string filenameConfInterval = "confInterval.dat";
    confInterval = new std::ofstream(filenameConfInterval.c_str(),std::basic_ios<char>::out);
    if (!confInterval){
      std::cerr << "\n confInterval.dat could not be initialized" << std::endl;
    }

    std::string filenameRhosOut = "rhosOut.dat";
    rhosOut = new std::ofstream(filenameRhosOut.c_str(),std::basic_ios<char>::out);
    if (!rhosOut){
      std::cerr << "\n rhosOut.dat could not be initialized" << std::endl;
    }

    std::string filenamenrOfFreqs = "nrOfFreqs.dat";
    nrOfFreqs = new std::ofstream(filenamenrOfFreqs.c_str(),std::basic_ios<char>::out);
    if (!rhosOut){
      std::cerr << "\n nrOfFreqs.dat could not be initialized" << std::endl;
    }

    Info->FinishProgress();

    // query parameters from "paramIdent" node (myParam_)
    myParam_->Get( "startFreq", startfreq_ );
    myParam_->Get( "stopFreq",  stopfreq_ );
    myParam_->Get( "numFreq", nrfreq_ );
    myParam_->Get( "numMechMeasurements", numMechMeasurements_ );  

    
    // should we calculate the impedance curve?
    myParam_->Get( "calcImpedanceCurve", CalcImpedanceCurve_, false );

    // should we calculate the mechanical displacement curve at selected node?
    myParam_->Get( "calcMechDisplCurve", CalcMechDisplCurve_, false );

    // at which node should the mechDisplCurve be calculated?
    //     keyVec="paramIdent", "mechDisplAtNode_";
    //     params->Get(keyVec, attrVec, valVec, mechDisplAtNode_);

    // further global constants for truncated Newton methods
    myParam_->Get( "maxNrInnerIterations", maxNumberInnerLoops_ );
    myParam_->Get( "maxNrOuterIterations", maxNumberNewtonLoops_ );
    myParam_->Get( "artDataNoise", delta_);
    

  } // end of constructor

  // destructor
  piezoParamIdent :: ~piezoParamIdent()
  {
    if (allMeasuredData)
      allMeasuredData->close();
    if (impedCurve)
      impedCurve->close();
    if(piezoLog)
      piezoLog->close();
    if(parLog)
      parLog->close();
    if(synMess)
      synMess->close();
    if(piezoLog)
      piezoLog->close();
    if(mechDispl)
      mechDispl->close();
    if(optimalFreqs)
      optimalFreqs->close();
    if(confInterval)
      confInterval->close();
    if(rhosOut)
      rhosOut->close();

    if(inMess_)
      inMess_.close();
    if(inMechMess_)
      inMechMess_.close();
  }

  void piezoParamIdent::SolveProblem(bool write_results, const std::string& comment) {
    assert(write_results = true);
    assert(comment == "");
    
    UInt highestAssumableNrOfMeasData=100;

    // notify resultHandler about beginning of new sequence step 
    ResultHandler * resHandler = domain->GetResultHandler();
    resHandler->BeginMultiSequenceStep( 1, analysis_, 1 );

    InitializePDEs();
    //domain->PrintGrid();
    DirectCoupledPDE* ptCoupledPDE =  domain->GetDirectCoupledPDE();
      
    ptPDE1_=domain-> GetSinglePDE("mechanic");
    ptPDE2_=domain-> GetSinglePDE("electrostatic");

    subdomsMech_ = ptPDE1_->getPDE_subdoms();
    subdomsElec_ = ptPDE1_->getPDE_subdoms();

    if(subdomsMech_.GetSize()==1){
      nrParameter_ = 10;
      actNrParameter = 10;
      actNrParameterC = 10;
    }
    else if (subdomsMech_.GetSize()>1){
    nrParameter_ = 16;
    actNrParameter = 16;
    actNrParameterC = 16;
    }

    parameter_.Resize(nrParameter_);
    parameterC_.Resize(nrParameter_);
    parameter_.Init();
    parameterC_.Init();

    whichParameterToUpdate_.Resize(nrParameter_);
    whichParameterToUpdateC_.Resize(nrParameter_);
    whichParameterToUpdate_.Init();
    whichParameterToUpdateC_.Init();

    whichParToUpInd_.Init();
    whichParToUpIndC_.Init();

    //     parameterIncrement.Resize(nrParameter_);
    //parameterIncrement = parameter_;
    omegas_.Resize(highestAssumableNrOfMeasData);
    freqs_.Resize(highestAssumableNrOfMeasData);
    freqsMech_.Resize(highestAssumableNrOfMeasData);
    real_.Resize(highestAssumableNrOfMeasData);
    imag_.Resize(highestAssumableNrOfMeasData);
    amplitude_phase.Resize(highestAssumableNrOfMeasData);
    F_hat_.Resize(highestAssumableNrOfMeasData);

    omegas_.Init();
    freqs_.Init();
    freqsMech_.Init();
    real_.Init();
    imag_.Init();
    amplitude_phase.Init();
    F_hat_.Init();  

    //will be overwritten in readMeasuredData ...
    nrMeasuredDataElec_=0;
    nrMeasuredDataMech_=0;

   
    // the following passage reads Data from file measuredData.dat
    // The rows are containing the values of the given frequencies, such as phase and amplitude!
    readMeasuredData(freqs_, real_, imag_, parameter_, voltage_, 
                     nrMeasuredData, thickness_, radius_, delta_);

    Vector<Double> freqsTemp;
    freqsTemp = freqs_;
    freqs_.Resize(nrMeasuredDataElec_);
    freqs_.Init();

    for(UInt i=0;i<nrMeasuredDataElec_;i++)
      freqs_[i]=freqsTemp[i];

    Vector<Double> freqsTempMech = freqsMech_;
    freqsMech_.Resize(nrMeasuredDataMech_);
    freqsMech_.Init();

    for(UInt i=0;i<nrMeasuredDataMech_;i++)
      freqsMech_[i]=freqsTempMech[i];

    nrMeasuredData = nrMeasuredDataElec_ + nrMeasuredDataMech_;
    
     std::cout<<"\nSampling points for impedance:" <<std::endl;
     std::cout<<freqs_<<std::endl;
     std::cout<<"Sampling points for mechanical measurements:" <<std::endl;
     std::cout<<freqsMech_<<std::endl;


    y_hat_.Resize(nrMeasuredData);
    //    bas.Resize(nrParameter_);
    res_NE_new.Resize(actNrParameter+actNrParameterC);
    res_NE.Resize(actNrParameter+actNrParameterC);
    lin_res.Resize(nrMeasuredData);
    res.Resize(nrMeasuredData);
    bas_bar.Resize(nrMeasuredData);
    s.Resize(actNrParameter+actNrParameterC);
    scaling_.Resize(nrParameter_);
    scalingC_.Resize(nrParameter_);
    F_hat_.Resize(nrMeasuredData);
    overall_res0.Resize(nrMeasuredData);
    parameter_new_.Resize(nrParameter_);

    actNrParameter=0;
    actNrParameterC=0;

    y_hat_.Init();
    res_NE_new.Init();
    res_NE.Init();
    lin_res.Init();
    res.Init();
    bas_bar.Init();
    s.Init();
    scaling_.Init(0.0);
    scalingC_.Init(0.0);
    F_hat_.Init();
    overall_res0.Init();
    parameter_new_.Init();


    ptPDE_->WriteGeneralPDEdefines();

    // how many parameters are there actually to set?
    for (UInt i=0;i<whichParameterToUpdate_.GetSize();i++)
      if (whichParameterToUpdate_[i]==1)
        actNrParameter++;

    for (UInt i=0;i<whichParameterToUpdateC_.GetSize();i++)
      if (whichParameterToUpdateC_[i]==1)
        actNrParameterC++;

    if (actNrParameter!=0){
      whichParToUpInd_.Resize(actNrParameter);
      whichParToUpInd_.Init();
    }


    if (actNrParameterC!=0) {
      whichParToUpIndC_.Resize(actNrParameterC);
      whichParToUpIndC_.Init();
    }

//      std::cout<<"actNrParameter:"<<std::endl;
//      std::cout<<"actNrParameterC:"<<std::endl;
//      std::cout<<actNrParameter<<std::endl;
//      std::cout<<actNrParameterC<<std::endl;

//      std::cout<<"WhichParameter to update " <<std::endl;
//      std::cout<<whichParameterToUpdate_<<std::endl;
//      std::cout<<whichParameterToUpdateC_<<std::endl;

    UInt intTemp=0;

    for (UInt i=0;i<whichParameterToUpdate_.GetSize();i++)
      if (whichParameterToUpdate_[i]==1){        
        whichParToUpInd_[intTemp]=i;
        intTemp++;
      }


    intTemp=0;
    for (UInt i=0;i<whichParameterToUpdateC_.GetSize();i++)
      if (whichParameterToUpdateC_[i]==1) {
        whichParToUpIndC_[intTemp]=i;
        intTemp++;
      }


    // ************************************************************************
    // Communication with BasePDE ... gets i.G. pointers to objects involved  *
    // ************************************************************************


    StdVector<BasePairCoupling*> ptCoupling = ptCoupledPDE->GetCouplingsObject();
    
    // if driver is not part of multiSequence Driver, get list
    // of pdes which have to be solved and intialize them
    
    ptMaterialMech_  = ptPDE1_[0].getPDEMaterialData();   // Pointer to mech. MaterialData
    ptMaterialElec_  = ptPDE2_[0].getPDEMaterialData();   // Pointer to elec. MaterialData
    ptMaterialPiezo_ = ptCoupling[0]->getPDEMaterialData();   // Pointer to piezo MaterialData

    
    //get the material tensors
    Matrix<Double> piezoMat,stiffMat, stiffMatAlu, stiffMatSteel, permMat;
    Matrix<Complex> piezoMatC, stiffMatC, permMatC;
    
    // Get first regionId of mechanic pde
    RegionIdType actId = subdomsMech_[0];

     if (subdomsMech_.GetSize()==1){
       ptMaterialPiezo_[actId]->GetTensor(piezoMat,PIEZO_TENSOR,REAL,FULL);
       ptMaterialMech_[actId]->GetTensor(stiffMat,MECH_STIFFNESS_TENSOR, REAL,FULL);
       ptMaterialElec_[actId]->GetTensor(permMat,ELEC_PERMITTIVITY,REAL,FULL);

       //       getchar();
     }


     if( imagMaterialParam_ ){
      updateComplexMaterialData(parameterC_);
    }

    updateMaterialData(parameter_);

    std::cout<<"Real parts of material parameters (initial values):"<<std::endl;
    std::cout<<parameter_<<std::endl;
    std::cout<<"Imaginary parts of material parameters (initial values):"<<std::endl;
    std::cout<<parameterC_<<std::endl;

    parameterIncrement_=parameter_;

    
    // <<<<<<<<<<<<<< calc mechanical displacement curve <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    if (CalcMechDisplCurve_ == true){
      Vector<Double> freqsTemp = freqs_;
      freqs_.Resize(nrfreq_);
      Double startFreqTemp;
      startFreqTemp=startfreq_;
      Double freqincr=(stopfreq_-startfreq_)/nrfreq_;
      for(UInt i=0;i<nrfreq_;i++){
        startFreqTemp+=freqincr;
        freqs_[i]=startFreqTemp;
      }
      calcMechDisplCurve();
      freqs_ = freqsTemp;
      std::cout<<"\n Press any key to continue ... "<<std::endl;
      //      getchar();
    }


    // <<<<<<<<<<<<<< calc impedance curve <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
     
    if (CalcImpedanceCurve_ == true){
      //      std::cout<<"++ Computing impedance curve ..."<<std::endl;
      Vector<Double> freqsTemp = freqs_;
      freqs_.Resize(nrfreq_);
      Double startFreqTemp;
      startFreqTemp=startfreq_;
      Double freqincr=(stopfreq_-startfreq_)/nrfreq_;
      for(UInt i=0;i<nrfreq_;i++){
        startFreqTemp+=freqincr;
        freqs_[i]=startFreqTemp;
      }
      calcImpedanceCurve();
      freqs_ = freqsTemp;
      std::cout<<"\n Press any key to continue ... "<<std::endl;
      getchar();
    }
   

    computeScaling();
  

    // xxxxxxxxxxxxxxxxxxxxxxx Choose different regularizing solvers here xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx //
    

    if (whichNewtonCG_==1){
      std::cout<<"++ Nonlinear Landweber's iteration ... " <<std::endl;
      Vector<Double> newFreqs;
      Vector<Double> newFreqsMech;
      readInMeasurement(newFreqs,newFreqsMech);
      std::cout<<"++ Fitting will be performed with the following " 
               << nrMeasuredData << " frequencies:" <<std::endl;
      std::cout<<freqs_ <<std::endl;
      std::cout<<freqsMech_ <<std::endl;
      nonlinLandweber();
      CalcImpedanceCurve_ = false;
    }

    if (whichNewtonCG_==7){
      if( imagMaterialParam_ ) {
        std::cerr<<" This version of nuMethods only works for real valued paraeters " <<std::endl;
        std::cerr<<" Choose either nuMethods=8 in your input file or " <<std::endl;
        std::cerr<<" use only real - valued parameters " <<std::endl;
        //        getchar();
      }
      std::cout<<"++ Newton - nu Methods " <<std::endl;
      UInt nrNuMethods=0;
      newtonCounter_=0;
      inner_eta_=1.0;

      Vector<Double> newFreqs;
      Vector<Double> newFreqsMech;
      readInMeasurement(newFreqs,newFreqsMech);
      calc_measuredCharge(freqs_, real_, imag_, y_hat_); // out of new measurements

      while (nrNuMethods<maxNumberNewtonLoops_){
 
        nuMethods();
        nrNuMethods++;

        if (nrNuMethods%10==0){
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
          if(mechDispl)
            mechDispl->close();
          
          std::string filename= "imped";
          filename.append(".dat");
          impedCurve = new std::ofstream(filename.c_str(),std::basic_ios<char>::out);
          filename="mechDispl.dat";
          mechDispl = new std::ofstream(filename.c_str(),std::basic_ios<char>::out);
          
          calcImpedanceCurve();
          freqs_ = freqsTemp;
        }
        

        for (UInt i=0; i<parameter_.GetSize(); i++)
          *parFinal<<parameter_[i]<<", ";
        *parFinal<<"/"<<std::endl;
        newtonCounter_++;
      }

    }

    else if (whichNewtonCG_==8){
      std::cout<<"++ Newton - nu MethodsC " <<std::endl;

      Vector<Double> newFreqs;
      Vector<Double> newFreqsMech;
      readInMeasurement(newFreqs,newFreqsMech);
      calc_measuredCharge(freqs_, real_, imag_, y_hat_); // out of new measurements

      UInt nrNuMethodsC=0;
      newtonCounter_=0;
      inner_eta_=1.0;
      while (nrNuMethodsC<maxNumberNewtonLoops_){

        nuMethodsC2();
        
        *parLog<<nrNuMethodsC << parameter_ << parameterC_<<std::endl;
        for (UInt i=0; i<parameter_.GetSize(); i++)
          *parFinal<<parameter_[i]<<", ";
        *parFinal<<"/"<<std::endl;
        for (UInt i=0; i<parameterC_.GetSize(); i++)
          *parFinal<<parameterC_[i]<<", ";
        *parFinal<<"/"<<std::endl;

        nrNuMethodsC++;
        newtonCounter_++;

      }
    }

    else if (whichNewtonCG_==9){
      std::cout<<"++ Opimtal experiment Design"<<std::endl;
      optimalExpDesign();
    }

    else if (whichNewtonCG_==10){
      std::cout<<"++ Optimal experiment Design - Complex parameters"<<std::endl;
      optimalExpDesign();
    }

    else if (whichNewtonCG_==11){
      std::cout<<"++ Optimal experiment Design - variable number of frequencies"<<std::endl;
      optimalExpDesignVarNrFreqs();
    }

    else if (whichNewtonCG_==12){

      std::cout<<"++ Least squares fitting started with maximal " << maxNumberNewtonLoops_ << " descent steps ... "<<std::endl;
      Vector<Double> newFreqs;
      Vector<Double> newFreqsMech;
      readInMeasurement(newFreqs,newFreqsMech);
      calc_measuredCharge(freqs_, real_, imag_, y_hat_); // out of new measurements
      std::cout<<"++ Fitting will be performed with the following " 
               << nrMeasuredData << " frequencies:" <<std::endl;

      std::cout<<freqs_ <<std::endl;
      std::cout<<freqsMech_ <<std::endl;

      leastSquare();
      CalcImpedanceCurve_ = false;
    }



    else
      std::cout<<"\n There was no valid NewtonCG method specified - see in your measuredData.dat -file "<<std::endl;
    //    tichonov();
    //  NewtonLandweber();
    //  createF(ptMaterial_, F_hat_,false); // calculates only forward problems over all omegas_

    // xxxxxxxxxxxxxxxxxxxxxxx End of choice xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx //

    if ( maxNumberNewtonLoops_!=0){
      std::cout<<"\n\n *** FINALLY CALCULATED PARAMETERS *** " <<std::endl;
      
      *parFinal<<"4) "<<std::endl;
      for (UInt i=0;i<parameter_.GetSize();i++){
        std::cout<<"par[" << i<<"]="<< parameter_[i]<<" + " << parameterC_[i]<<"i"<<std::endl;
        *parFinal<<parameter_[i]<<", ";
      }
      *parFinal<<"/"<<std::endl;
      std::cout<<"\n ****************************************** \n " <<std::endl;
    }


    // <<<<<<<<<<<<<< for a hopefully nice imped curve after identification !! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

     
    if (CalcImpedanceCurve_ == true && maxNumberNewtonLoops_!=0){
      Vector<Double> freqsTemp = freqs_;
      freqs_.Resize(nrfreq_);
      Double freqincr=(stopfreq_-startfreq_)/nrfreq_;
      for(UInt i=0;i<nrfreq_;i++){
        startfreq_+=freqincr;
        freqs_[i]=startfreq_;
      }

      if (impedCurve)
        impedCurve->close();
      if(mechDispl)
        mechDispl->close();
      
      std::string filename= "imped";
      filename.append(".dat");
      impedCurve = new std::ofstream(filename.c_str(),std::basic_ios<char>::out);
      filename="mechDispl.dat";
      mechDispl = new std::ofstream(filename.c_str(),std::basic_ios<char>::out);
      
      calcImpedanceCurve();
      
      freqs_ = freqsTemp;
    }

    // notify resultHandler about finishing of current sequence step
      resHandler->FinishMultiSequenceStep(); 
  }// End solveProblem

  
  // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXxx - now some methods are following ...

  //! Updates material data & updates system matrices!!
  void piezoParamIdent::updateMaterialData(Vector<Double> & parameter_){
    //    std::cout<<"++ updateMaterialData " <<std::endl;

    UInt dim=ptPDE1_->getPDE_spaceDim();
    
    Matrix<Double> stiffTensor;
    Matrix<Double> piezoTensor;
    Matrix<Double> permTensor;

    stiffTensor.Resize(6,6);
    piezoTensor.Resize(3,6);
    permTensor.Resize(3,3);

    stiffTensor.Init(0);
    piezoTensor.Init(0);
    permTensor.Init(0);


    subdomsMech_ = ptPDE1_->getPDE_subdoms();

    RegionIdType actId = subdomsMech_[0];

    
    // Materialparameteridentifizierung fuer nur eine Scheibe:
    if (subdomsMech_.GetSize()==1){         
      stiffTensor[0][0] = parameter_[0]; //c_11
      stiffTensor[1][1] = parameter_[0]; //c_11
      stiffTensor[2][2] = parameter_[1]; //c_33
      stiffTensor[0][1] = parameter_[2]; //c_12
      stiffTensor[1][0] = parameter_[2]; //c_12
      stiffTensor[0][2] = parameter_[3]; //c_13
      stiffTensor[2][0] = parameter_[3]; //c_13
      stiffTensor[1][2] = parameter_[3]; //c_13
      stiffTensor[2][1] = parameter_[3]; //c_13
      stiffTensor[3][3] = parameter_[4]; //c_44
      stiffTensor[4][4] = parameter_[4]; //c_44
      stiffTensor[5][5] = 0.5*(parameter_[0]-parameter_[2]); //c_66

      if (dim==2){
        stiffTensor[0][0] = parameter_[0]; //c_11
        stiffTensor[2][2] = parameter_[0]; //c_11
        stiffTensor[1][1] = parameter_[1]; //c_33
        stiffTensor[0][2] = parameter_[2]; //c_12
        stiffTensor[2][0] = parameter_[2]; //c_12
        stiffTensor[0][1] = parameter_[3]; //c_13
        stiffTensor[1][0] = parameter_[3]; //c_13
        stiffTensor[1][2] = parameter_[3]; //c_13
        stiffTensor[2][1] = parameter_[3]; //c_13
        stiffTensor[3][3] = parameter_[4]; //c_44
        stiffTensor[5][5] = parameter_[4]; //c_44
        stiffTensor[4][4] = 0.5*(parameter_[0]-parameter_[2]); //c_66
      }
 
      piezoTensor[1][3]=parameter_[5];  //e_15
      piezoTensor[0][4]=parameter_[5];  //e_15
      piezoTensor[2][0]=parameter_[6]; //e_31
      piezoTensor[2][1]=parameter_[6]; //e_31
      piezoTensor[2][2]=parameter_[7]; // e_33

      if (dim==2){
        piezoTensor[2][3]=parameter_[5];  //e_15
        piezoTensor[0][5]=parameter_[5];  //e_15
        piezoTensor[1][0]=parameter_[6]; //e_31
        piezoTensor[1][2]=parameter_[6]; //e_31
        piezoTensor[1][1]=parameter_[7]; // e_33
      }

      permTensor[0][0] = parameter_[8]; //eps_11
      permTensor[1][1] = parameter_[8]; //eps_11
      permTensor[2][2] = parameter_[9]; //eps_33

      if(dim==2){
        permTensor[0][0] = parameter_[8]; //eps_11
        permTensor[2][2] = parameter_[8]; //eps_11
        permTensor[1][1] = parameter_[9]; //eps_33
      }

   
      ptMaterialPiezo_[actId]->SetTensor(piezoTensor,PIEZO_TENSOR,REAL);
      ptMaterialMech_[actId]->SetTensor(stiffTensor,MECH_STIFFNESS_TENSOR, REAL);
      ptMaterialElec_[actId]->SetTensor(permTensor,ELEC_PERMITTIVITY,REAL);

    }
   
  } // end updateMaterialData

  void piezoParamIdent::updateComplexMaterialData(Vector<Double> & parameterC_){
    //    std::cout<<"++ updateComplexMaterialData " <<std::endl;

    Matrix<Double> stiffTensorC;
    Matrix<Double> piezoTensorC;
    Matrix<Double> permTensorC;

    stiffTensorC.Resize(6,6);
    piezoTensorC.Resize(3,6);
    permTensorC.Resize(3,3);

    stiffTensorC.Init(0);
    piezoTensorC.Init(0); 
    permTensorC.Init(0);


    UInt dim=ptPDE1_->getPDE_spaceDim();
    subdomsMech_ = ptPDE1_->getPDE_subdoms();
    
    if (subdomsMech_.GetSize()==1){

      stiffTensorC[0][0] = parameterC_[0]; //c_11
      stiffTensorC[1][1] = parameterC_[0]; //c_11
      stiffTensorC[2][2] = parameterC_[1]; //c_33
      stiffTensorC[0][1] = parameterC_[2]; //c_12
      stiffTensorC[1][0] = parameterC_[2]; //c_12
      stiffTensorC[0][2] = parameterC_[3]; //c_13
      stiffTensorC[2][0] = parameterC_[3]; //c_13
      stiffTensorC[1][2] = parameterC_[3]; //c_13
      stiffTensorC[2][1] = parameterC_[3]; //c_13
      stiffTensorC[3][3] = parameterC_[4]; //c_44
      stiffTensorC[4][4] = parameterC_[4]; //c_44
      stiffTensorC[5][5] = 0.5*(parameterC_[0]-parameterC_[2]); //c_66

      if (dim==2){
        stiffTensorC[0][0] = parameterC_[0]; //c_11
        stiffTensorC[2][2] = parameterC_[0]; //c_11
        stiffTensorC[1][1] = parameterC_[1]; //c_33
        stiffTensorC[0][2] = parameterC_[2]; //c_12
        stiffTensorC[2][0] = parameterC_[2]; //c_12
        stiffTensorC[0][1] = parameterC_[3]; //c_13
        stiffTensorC[1][0] = parameterC_[3]; //c_13
        stiffTensorC[1][2] = parameterC_[3]; //c_13
        stiffTensorC[2][1] = parameterC_[3]; //c_13
        stiffTensorC[3][3] = parameterC_[4]; //c_44
        stiffTensorC[5][5] = parameterC_[4]; //c_44
        stiffTensorC[4][4] = 0.5*(parameterC_[0]-parameterC_[2]); //c_66
      }
      
      piezoTensorC[1][3]=parameterC_[5];  //e_15
      piezoTensorC[0][4]=parameterC_[5];  //e_15
      piezoTensorC[2][0]=parameterC_[6]; //e_31
      piezoTensorC[2][1]=parameterC_[6]; //e_31
      piezoTensorC[2][2]=parameterC_[7]; // e_33

      if (dim==2){
        piezoTensorC[2][3]=parameterC_[5];  //e_15
        piezoTensorC[0][5]=parameterC_[5];  //e_15
        piezoTensorC[1][0]=parameterC_[6]; //e_31
        piezoTensorC[1][2]=parameterC_[6]; //e_31
        piezoTensorC[1][1]=parameterC_[7]; // e_33
      }
      
      permTensorC[0][0] = parameterC_[8]; //eps_11
      permTensorC[1][1] = parameterC_[8]; //eps_11
      permTensorC[2][2] = parameterC_[9]; //eps_33

      if(dim==2){
        permTensorC[0][0] = parameterC_[8]; //eps_11
        permTensorC[2][2] = parameterC_[8]; //eps_11
        permTensorC[1][1] = parameterC_[9]; //eps_33
      }
  
      // Get first regionId of mechanic pde
      RegionIdType actId = subdomsMech_[0];
      
      ptMaterialPiezo_[actId]->SetTensor(piezoTensorC,PIEZO_TENSOR,IMAG);
      ptMaterialMech_[actId]->SetTensor(stiffTensorC,MECH_STIFFNESS_TENSOR, IMAG);
      ptMaterialElec_[actId]->SetTensor(permTensorC,ELEC_PERMITTIVITY,IMAG);

    }

    
  } // end updateComplexMaterialData

  void piezoParamIdent::setNewParameterSet(Vector<Double> & par,
                                           Vector<Double> &  par_new,
                                           Vector<Double> & scaling_,
                                           Double & theta,Vector<Double> & uStep,
                                           Vector<UInt> & whichParameterToUpdate_){
    UInt helpInd=0;
    for (UInt i=0;i<nrParameter_;i++){
      //      std::cout<<"\n setNewParameterSet " << whichParameterToUpdate_[i]<<", ";
      if (whichParameterToUpdate_[i]==1){
        par_new[i]=par[i]+(1.0/scaling_[i])*theta*uStep[helpInd];
        //      std::cout<<"\n parNew = " << par_new[i]<<", step = " << uStep[helpInd] << std::endl;
        helpInd++;
      }
    }
    std::cout<<"\n-----------------------------"<<std::endl;
    //    getchar();


  } // end setNewParameterSet


 void piezoParamIdent::computeScaling(){

//     Matrix<Double> piezoMat,stiffMat, stiffMatAlu, stiffMatSteel, permMat, stiffMatSchraube;
//     Matrix<Double> piezoMatC, stiffMatAluC, stiffMatSteelC, stiffMatC, permMatC, stiffMatSchraubeC;

//     stiffMat.Resize(6,6);
//     stiffMat.Init();
    
//     piezoMat.Resize(3,6);
//     piezoMat.Init();
    
//     permMat.Resize(3,3);
//     permMat.Init();

//     stiffMatC.Resize(6,6);
//     stiffMatC.Init();
    
//     piezoMatC.Resize(3,6);
//     piezoMatC.Init();
    
//     permMatC.Resize(3,3);
//     permMatC.Init();

//     // Get first regionId of mechanic pde
//     RegionIdType actId = subdomsMech_[0];
    
//     if (subdomsMech_.GetSize()==1){
//       ptMaterialPiezo_[actId]->GetTensor(piezoMat,PIEZO_TENSOR,REAL,FULL);
//       ptMaterialMech_[actId]->GetTensor(stiffMat,MECH_STIFFNESS_TENSOR, REAL,FULL);
//       ptMaterialElec_[actId]->GetTensor(permMat,ELEC_PERMITTIVITY,REAL,FULL);

//       scaling_.Resize(nrParameter_);
//       scaling_.Init();

//       scaling_[0]=1.0/(stiffMat[0][0]); 
//       scaling_[1]=1.0/(stiffMat[1][1]);
//       scaling_[2]=1.0/(stiffMat[1][0]);
//       scaling_[3]=1.0/(stiffMat[0][2]);
//       scaling_[4]=1.0/(stiffMat[3][3]); 
//       scaling_[5]=1.0/(piezoMat[1][3]);
//       scaling_[6]=std::abs(1.0/((piezoMat)[2][0]));
//       scaling_[7]=1.0/((piezoMat)[2][2]);
//       scaling_[8]=1.0/((permMat)[0][0]); 
//       scaling_[9]=1.0/((permMat)[2][2]);


//       if( imagMaterialParam_ ) {
//         ptMaterialPiezo_[actId]->GetTensor(piezoMatC,PIEZO_TENSOR,IMAG,FULL);
//         ptMaterialMech_[actId]->GetTensor(stiffMatC,MECH_STIFFNESS_TENSOR,IMAG,FULL);
//         ptMaterialElec_[actId]->GetTensor(permMatC,ELEC_PERMITTIVITY,IMAG,FULL);

//  //        std::cout<<"piezoMatC"<<std::endl;
// //         std::cout<<piezoMatC<<std::endl;

// //         std::cout<<"stiffMatC"<<std::endl;
// //         std::cout<<stiffMatC<<std::endl;

// //         std::cout<<"permMatC"<<std::endl;
// //         std::cout<<permMatC<<std::endl;

//         scalingC_.Resize(nrParameter_);
//         scalingC_.Init();

//         scalingC_[0]=1.0/(stiffMatC[0][0]); 
//         scalingC_[1]=1.0/(stiffMatC[1][1]);
//         scalingC_[2]=1.0/(stiffMatC[1][0]);
//         scalingC_[3]=1.0/(stiffMatC[0][2]);
//         scalingC_[4]=1.0/(stiffMatC[3][3]); 
//         scalingC_[5]=1.0/(piezoMatC[1][3]);
//         scalingC_[6]=std::abs(1.0/((piezoMatC)[2][0]));
//         scalingC_[7]=1.0/((piezoMatC)[2][2]);
//         scalingC_[8]=1.0/((permMatC)[0][0]); 
//         scalingC_[9]=1.0/((permMatC)[2][2]);
//       }
//     }

//     else{
  
//       ptMaterialMech_[1]->GetTensor(stiffMatAlu,MECH_STIFFNESS_TENSOR, REAL,FULL);
//       ptMaterialPiezo_[2]->GetTensor(piezoMat,PIEZO_TENSOR,REAL,FULL);
//       ptMaterialMech_[2]->GetTensor(stiffMat,MECH_STIFFNESS_TENSOR, REAL,FULL);
//       ptMaterialElec_[2]->GetTensor(permMat,ELEC_PERMITTIVITY,REAL,FULL);
//       // ptMaterialMech_[3]->GetTensor(stiffMat,MECH_STIFFNESS_TENSOR, REAL,FULL);
//       ptMaterialMech_[4]->GetTensor(stiffMatSteel,MECH_STIFFNESS_TENSOR, REAL,FULL);
//       ptMaterialMech_[5]->GetTensor(stiffMatSchraube,MECH_STIFFNESS_TENSOR, REAL,FULL);

//       scaling_[0]=1.0/(stiffMat[0][0]); 
//       scaling_[1]=1.0/(stiffMat[2][2]);
//       scaling_[2]=1.0/(stiffMat[1][0]);
//       scaling_[3]=1.0/(stiffMat[0][2]);
//       scaling_[4]=1.0/(stiffMat[3][3]); 
//       scaling_[5]=1.0/(piezoMat[1][3]);
//       scaling_[6]=std::abs(1.0/((piezoMat)[2][0]));
//       scaling_[7]=1.0/((piezoMat)[2][2]);
//       scaling_[8]=1.0/((permMat)[0][0]); 
//       scaling_[9]=1.0/((permMat)[2][2]);

//       scaling_[10]=1.0/(stiffMatSteel[3][3]); 
//       scaling_[11]=1.0/(stiffMatSteel[0][1]); 

//       scaling_[12]=1.0/(stiffMatAlu[3][3]); 
//       scaling_[13]=1.0/(stiffMatAlu[0][1]);

//       scaling_[14]=1.0/(stiffMatSchraube[3][3]); 
//       scaling_[15]=1.0/(stiffMatSchraube[0][1]);

        
//       if( imagMaterialParam_ ) {

//         ptMaterialMech_[1]->GetTensor(stiffMatAluC,MECH_STIFFNESS_TENSOR,IMAG,FULL);
//         ptMaterialPiezo_[2]->GetTensor(piezoMatC,PIEZO_TENSOR,IMAG,FULL);
//         ptMaterialMech_[2]->GetTensor(stiffMatC,MECH_STIFFNESS_TENSOR,IMAG,FULL);
//         ptMaterialElec_[2]->GetTensor(permMatC,ELEC_PERMITTIVITY,IMAG,FULL);
//         ptMaterialMech_[4]->GetTensor(stiffMatSteelC,MECH_STIFFNESS_TENSOR, IMAG,FULL);
//         ptMaterialMech_[5]->GetTensor(stiffMatSchraubeC,MECH_STIFFNESS_TENSOR, REAL,FULL);

//         scalingC_[0]=1.0/(stiffMatC[0][0]); 
//         scalingC_[1]=1.0/(stiffMatC[2][2]);
//         scalingC_[2]=1.0/(stiffMatC[1][0]);
//         scalingC_[3]=1.0/(stiffMatC[0][2]);
//         scalingC_[4]=1.0/(stiffMatC[3][3]); 
//         scalingC_[5]=1.0/(piezoMatC[1][3]);
//         scalingC_[6]=std::abs(1.0/((piezoMatC)[2][0]));
//         scalingC_[7]=1.0/((piezoMatC)[2][2]);
//         scalingC_[8]=1.0/((permMatC)[0][0]); 
//         scalingC_[9]=1.0/((permMatC)[2][2]);
          
//         scalingC_[10]=1.0/(stiffMatSteelC[3][3]); 
//         scalingC_[11]=1.0/(stiffMatSteelC[0][1]); 
          
//         scalingC_[12]=1.0/(stiffMatAluC[3][3]); 
//         scalingC_[13]=1.0/(stiffMatAluC[0][1]);                

//         scalingC_[14]=1.0/(stiffMatSchraubeC[3][3]); 
//         scalingC_[15]=1.0/(stiffMatSchraubeC[0][1]);        


//       }

      for (UInt i=0;i<parameter_.GetSize(); i++)
        if (parameter_[i]!=0.0)
          scaling_[i]= std::abs(1.0/parameter_[i]);
        else scaling_[i]=1.0;

      for (UInt i=0;i<parameterC_.GetSize(); i++)
        if (parameterC_[i]!=0.0)
          scalingC_[i]= std::abs(1.0/parameterC_[i]);
        else scalingC_[i]=1.0;

      // }

 }//end compute scaling


  void piezoParamIdent::writeTensorsInFile(){

    Matrix<Complex> piezoMatC, stiffMatC, permMatC;

    // Get first regionId of mechanic pde
    RegionIdType actId = subdomsMech_[0];

    ptMaterialPiezo_[actId]->GetTensor(piezoMatC,PIEZO_TENSOR,COMPLEX,FULL);
    ptMaterialMech_[actId]->GetTensor(stiffMatC,MECH_STIFFNESS_TENSOR, COMPLEX,FULL);
    ptMaterialElec_[actId]->GetTensor(permMatC,ELEC_PERMITTIVITY,COMPLEX,FULL);
    
    Matrix<Complex> sE(6,6);
    sE.Init();
    for(UInt i=0;i<6;i++)
      sE[i][i]=Complex(1.0,0.0);
        
#ifdef USE_LAPACK
    
    lapackSysMatType LAPACK_SYS_MAT_TYPE = ZGESV;
    stiffMatC.solveWithLapack(sE,LAPACK_SYS_MAT_TYPE);

#endif

    *allTensors<<"Mechanical Modulus Tensor cE = " <<std::endl;
    *allTensors<<stiffMatC<<std::endl;

    *allTensors<<"Mechanical Elasticity Tensor sE = " <<std::endl;
    *allTensors<<sE<<std::endl;


    *allTensors<<"Piezoelectric Coupling Tensor e = " <<std::endl;
    *allTensors<<piezoMatC<<std::endl;


    Matrix<Complex> d(piezoMatC.GetSizeRow(),piezoMatC.GetSizeCol());
    d.Init();
    d=piezoMatC*sE;  
    *allTensors<<"Piezoelectric Coupling Tensor d = " <<std::endl;
    *allTensors<<d<<std::endl;


    *allTensors<<"Permittivity Tensor epsS = " <<std::endl;
    *allTensors<<permMatC<<std::endl;

    Matrix<Complex> epsT (3,3);
    Matrix<Complex> epsTrans;
    piezoMatC.Transpose(epsTrans);

    epsT.Init();
   
    epsT = permMatC + d*epsTrans;

    *allTensors<<"Permittivity Tensor epsT = " <<std::endl;
    *allTensors<<epsT<<std::endl;

    Complex influenceConstant = Complex(1.0/8.8542e-12,0.0);
    permMatC *= influenceConstant;
    epsT *= influenceConstant;


    *allTensors<<"Relative Permittivity Tensor epsS = " <<std::endl;
    *allTensors<<permMatC<<std::endl;

    *allTensors<<"Relative Permittivity Tensor epsT = " <<std::endl;
    *allTensors<<epsT<<std::endl;


    std::cout<<"++ All Material Tesnors are written into File 'allTensors.dat' " <<std::endl;
   

  
  }


} // end namespace CoupledField

 


































































