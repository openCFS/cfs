#include <iomanip>
#include <PDE/SinglePDE.hh>
#include "Domain/domain.hh"
#include "piezoParamIdent.hh"
#include "CoupledPDE/DirectCoupledPDE.hh"
#include "CoupledPDE/BasePairCoupling.hh"



namespace CoupledField
{

  // ========================================================================
  // ========================= piezoParamIdent - Part ===========================
  // ========================================================================

  //constructor
  // opens datafiles: measuredData.dat for input, imedCurve.dat and piezoLog.dat for output

  piezoParamIdent :: piezoParamIdent(Domain * adomain,
                                     Integer stepOffset,
                                     Double timeOffset,
                                     std::string driverTag,
                                     Boolean isPartOfSequence)
    :SingleDriver(adomain, stepOffset, timeOffset, 
                  driverTag, isPartOfSequence){

    ENTER_FCN( "piezoParamIdent::piezoParamIdent" );

    ptDomain_ = adomain;
    ptMyPDE_ = NULL;
    residuumParIdent_=1.0;
    resonanceFrequency_=0;
    antiResonanceFrequency_=0;
    tau_=1.0;
    sign_=1.0;
    StdVector<std::string> pdeList;
    params->GetPDEList( pdeList );

    if (pdeList[0]=="piezo")
      directCoupling_=FALSE;
    else 
      directCoupling_=TRUE;

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

    Info->StartProgress( "\n Opening output files \n imped.dat, piezoLog.dat, parLog.dat, mechDispl.dat, parFinal.dat ... " );

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


    // in future, several parameters wwill be taken from the xml - file ...
    StdVector<std::string> keyVec, attrVec, valVec;
    
    attrVec = "tag";
    valVec = driverTag_;
    
    // Get time stepping information from parameter object
    keyVec = "paramIdent", "startFreq";
    params->Get(keyVec, attrVec, valVec, startfreq_);
    
    keyVec = "paramIdent", "stopFreq";
    params->Get(keyVec, attrVec, valVec, stopfreq_);
    
    keyVec = "paramIdent", "numFreq";
    params->Get(keyVec, attrVec, valVec, nrfreq_);

    
    // should we calculate the impedance curve?
    CalcImpedanceCurve_ = params->IsSet("calcImpedanceCurve",  "paramIdent");

    // should we calculate the mechanical displacement curve at selected node?
    CalcMechDisplCurve_ = params->IsSet("calcMechDisplCurve",  "paramIdent");

    // at which node should the mechDisplCurve be calculated?
    //     keyVec="paramIdent", "mechDisplAtNode_";
    //     params->Get(keyVec, attrVec, valVec, mechDisplAtNode_);

    // DOF of mechDispl.
    keyVec="paramIdent", "dofOfMechDispl";
    params->Get(keyVec, attrVec, valVec, dofOfMechDispl_);


    // further important constants for truncated Newton methods
    keyVec="paramIdent", "maxNrInnerIterations";
    params->Get(keyVec, attrVec, valVec, maxNumberInnerLoops_);

    keyVec="paramIdent", "maxNrOuterIterations";
    params->Get(keyVec, attrVec, valVec, maxNumberNewtonLoops_);

    keyVec="paramIdent", "artDataNoise";
    params->Get(keyVec, attrVec, valVec, delta_);
    

  } // end of constructor

  // destructor
  piezoParamIdent :: ~piezoParamIdent()
  {
    ENTER_FCN( "piezoParamIdent::~piezoParamIdent" );
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
    //     if(mess)
    //       mess->close();


#ifdef USE_LAPACK
    if(lp_af77)
      DeleteArray(lp_af77);
    std::cout<<"delete wf77 in destructor ..."<<std::endl;
    if(&lp_wf77)
      DeleteArray(lp_wf77);
    std::cout<<"delete workf77"<<std::endl;
    if(lp_workf77)
      DeleteArray(lp_workf77);
    std::cout<<"delete rwork77"<<std::endl;
    if(lp_rworkf77)
      DeleteArray(lp_rworkf77);
#endif
  }

  void piezoParamIdent :: SolveProblem() {
    ENTER_FCN( "piezoParamIdent::SolveProblem" );


    UInt highestAssumableNrOfMeasData=25;
    
    ptDomain_->PrintGrid();
    GetMyPDEs();
    DirectCoupledPDE* ptCoupledPDE =  ptDomain_->GetDirectCoupledPDE();
    Info->StartProgress ("Starting to solve problem", FALSE);
  
    ptPDE1_=ptDomain_-> GetSinglePDE("mechanic");
    ptPDE2_=ptDomain_-> GetSinglePDE("electrostatic");

    subdomsMech_ = ptPDE1_->getPDE_subdoms();
    subdomsElec_ = ptPDE1_->getPDE_subdoms();

    if(subdomsMech_.GetSize()==1){
      nrParameter_ = 10;
      actNrParameter = 10;
      actNrParameterC = 10;
    }
    else if (subdomsMech_.GetSize()>1){
    nrParameter_ = 14;
    actNrParameter = 14;
    actNrParameterC = 14;
    }

    parameter_.Resize(nrParameter_);
    parameterC_.Resize(nrParameter_);
    whichParameterToUpdate_.Resize(nrParameter_);
    whichParameterToUpdateC_.Resize(nrParameter_);

    //     parameterIncrement.Resize(nrParameter_);
    //parameterIncrement = parameter_;
    omegas_.Resize(highestAssumableNrOfMeasData);
    freqs_.Resize(highestAssumableNrOfMeasData);
    real_.Resize(highestAssumableNrOfMeasData);
    imag_.Resize(highestAssumableNrOfMeasData);
    amplitude_phase.Resize(highestAssumableNrOfMeasData);
    F_hat_.Resize(highestAssumableNrOfMeasData);
    
    // the following passage reads Data from file measuredData.dat
    // The rows are containing the values of the given frequencies, such as phase and amplitude!
    readMeasuredData(freqs_, real_, imag_, parameter_, voltage_, 
                     nrMeasuredData, thickness_, radius_, delta_);

    Vector<Double> freqsTemp = freqs_;
    freqs_.Resize(nrMeasuredData);

    for(UInt i=0;i<nrMeasuredData;i++)
      freqs_[i]=freqsTemp[i];

    y_hat_.Resize(2*nrMeasuredData);
    //    bas.Resize(nrParameter_);
    res_NE_new.Resize(actNrParameter+actNrParameterC);
    res_NE.Resize(actNrParameter+actNrParameterC);
    lin_res.Resize(2*nrMeasuredData);
    res.Resize(2*nrMeasuredData);
    bas_bar.Resize(2*nrMeasuredData);
    s.Resize(actNrParameter+actNrParameterC);
    scaling_.Resize(nrParameter_);
    scalingC_.Resize(nrParameter_);
    F_hat_.Resize(2*nrMeasuredData);
    overall_res0.Resize(2*nrMeasuredData);
    parameter_new_.Resize(nrParameter_);

    actNrParameter=0;
    actNrParameterC=0;


 

    ptPDE_->WriteGeneralPDEdefines();

    // how many parameters are there actually to set?
    for (UInt i=0;i<whichParameterToUpdate_.GetSize();i++)
      if (whichParameterToUpdate_[i]==1)
        actNrParameter++;

    for (UInt i=0;i<whichParameterToUpdateC_.GetSize();i++)
      if (whichParameterToUpdateC_[i]==1)
        actNrParameterC++;
    if (actNrParameter!=0)   
      whichParToUpInd_.Resize(actNrParameter);

    if (whichNewtonCG_==4||whichNewtonCG_==6||whichNewtonCG_==8
        ||whichNewtonCG_==10||whichNewtonCG_==12)
      if (actNrParameterC!=0)   
        whichParToUpIndC_.Resize(actNrParameterC);

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
    
  
    //     if( params->HasValue( "type", "imagMaterialParameter", "materialDataType" ) ){
    //       std::cerr << "should not be here " <<std::endl;
    //       ptMaterialPiezo_[0]->GetTensor(piezoMatC,PIEZO_TENSOR,IMAG,FULL);
    //       ptMaterialMech_[0]->GetTensor(stiffMatC,MECH_STIFFNESS_TENSOR,IMAG,FULL);
    //       ptMaterialElec_[0]->GetTensor(permMatC,ELEC_PERMITTIVITY,IMAG,FULL);
    //       std::cout<< " Mech-Tensor (Imag) \n" << stiffMatC  << std::endl;
    //       std::cout<< " Elec-Tensor (Imag) \n" << permMatC   << std::endl;
    //       std::cout<< " Piezo-Tensor (Imag) \n" << piezoMatC << std::endl;
    //     }


    if (subdomsMech_.GetSize()==1){
      ptMaterialPiezo_[0]->GetTensor(piezoMat,PIEZO_TENSOR,REAL,FULL);
      ptMaterialMech_[0]->GetTensor(stiffMat,MECH_STIFFNESS_TENSOR, REAL,FULL);
      ptMaterialElec_[0]->GetTensor(permMat,ELEC_PERMITTIVITY,REAL,FULL);
    }
    else {      
      ptMaterialMech_[1]->GetTensor(stiffMatAlu,MECH_STIFFNESS_TENSOR, REAL,FULL);

      ptMaterialPiezo_[2]->GetTensor(piezoMat,PIEZO_TENSOR,REAL,FULL);
      ptMaterialMech_[2]->GetTensor(stiffMat,MECH_STIFFNESS_TENSOR, REAL,FULL);
      ptMaterialElec_[2]->GetTensor(permMat,ELEC_PERMITTIVITY,REAL,FULL);
      
      ptMaterialMech_[3]->GetTensor(stiffMat,MECH_STIFFNESS_TENSOR, REAL,FULL);

      //ptMaterialMech_[4]->GetTensor(stiffMatSteel,MECH_STIFFNESS_TENSOR, REAL,FULL);
      ptMaterialMech_[5]->GetTensor(stiffMatSteel,MECH_STIFFNESS_TENSOR, REAL,FULL);

      std::cout<< " Mech-Tensor Steel (Real) \n" << stiffMatSteel   << std::endl;
      std::cout<< " Mech-Tensor Alu (Real) \n" << stiffMatAlu << std::endl;

    }
    std::cout<< " Mech-Tensor (Real) \n" << stiffMat  << std::endl;
    std::cout<< " Elec-Tensor (Real) \n" << permMat   << std::endl;
    std::cout<< " Piezo-Tensor (Real) \n" << piezoMat << std::endl;

    
    updateMaterialData(parameter_);

    if( params->HasValue( "type", "imagMaterialParameter", "materialDataType" ) )
      updateComplexMaterialData(parameterC_);

    parameterIncrement_=parameter_;

    
    // <<<<<<<<<<<<<< calc mechanical displacement curve <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    if (CalcMechDisplCurve_ == TRUE){
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
      getchar();
    }


    // <<<<<<<<<<<<<< calc impedance curve <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
     
    if (CalcImpedanceCurve_ == TRUE){
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
    

    if (whichNewtonCG_==7){
      if( params->HasValue( "type", "imagMaterialParameter", "materialDataType" ) ){
        std::cerr<<" This version of nuMethods only works for real valued paraeters " <<std::endl;
        std::cerr<<" Choose either nuMethods=8 in your input file or " <<std::endl;
        std::cerr<<" use only real - valued parameters " <<std::endl;
        getchar();
      }
      std::cout<<"++ Newton - nu Methods " <<std::endl;
      UInt nrNuMethods=0;
      newtonCounter_=0;
      inner_eta_=1.0;

      Vector<Double> newFreqs;
      readInMeasurement(newFreqs);
      calc_measuredCharge(freqs_, real_, imag_, y_hat_); // out of new measurements

      while (nrNuMethods<maxNumberNewtonLoops_){
 
        nuMethods();
        nrNuMethods++;

        for (UInt i=0; i<parameter_.GetSize(); i++)
          *parFinal<<parameter_[i]<<", ";
        *parFinal<<"/"<<std::endl;
        newtonCounter_++;
      }

    }

    else if (whichNewtonCG_==8){
      std::cout<<"++ Newton - nu MethodsC " <<std::endl;

      Vector<Double> newFreqs;
      readInMeasurement(newFreqs);
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
      Vector<Double> newFreqs;
      readInMeasurement(newFreqs);

      std::cout<<" we red in the following frequencies " <<std::endl;
      std::cout<<freqs_ <<std::endl;

      calc_measuredCharge(freqs_, real_, imag_, y_hat_); // out of new measurements
      std::cout<<"++ Least square fitting"<<std::endl;
      leastSquare();
    }



    else
      std::cout<<"\n There was no valid NewtonCG method specified - see in your measuredData.dat -file "<<std::endl;
    //    tichonov();
    //  NewtonLandweber();
    //  createF(ptMaterial_, F_hat_,FALSE); // calculates only forward problems over all omegas_

    // xxxxxxxxxxxxxxxxxxxxxxx End of choice xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx //

    if ( maxNumberNewtonLoops_!=0){
      std::cout<<"\n\n *** FINALLY CALCULATED PARAMETERS *** ... here they are: " <<std::endl;
      
      *parFinal<<"4) "<<std::endl;
      for (UInt i=0;i<parameter_.GetSize();i++){
        std::cout<<"par[" << i<<"]="<< parameter_[i]<<" + " << parameterC_[i]<<"i"<<std::endl;
        *parFinal<<parameter_[i]<<", ";
      }
      *parFinal<<"/"<<std::endl;
    }


    // <<<<<<<<<<<<<< for a hopefully nice imped curve after identification !! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

     
    if (CalcImpedanceCurve_ == TRUE && maxNumberNewtonLoops_!=0){
      Vector<Double> freqsTemp = freqs_;
      freqs_.Resize(nrfreq_);
      Double freqincr=(stopfreq_-startfreq_)/nrfreq_;
      for(UInt i=0;i<nrfreq_;i++){
        startfreq_+=freqincr;
        freqs_[i]=startfreq_;
      }
      calcImpedanceCurve();
      std::cout<<"\n Press any key to continue ... "<<std::endl;
      freqs_ = freqsTemp;
    }
    
  }// End solveProblem

  
  // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXxx - now some methods are following ...

  //! Updates material data & updates system matrices!!
  void piezoParamIdent::updateMaterialData(Vector<Double> & parameter_){
    ENTER_FCN("piezoParamIdent::updateMaterialData");    
    
    Matrix<Double> stiffTensor,stiffTensorSteel,stiffTensorAlu;
    Matrix<Double> piezoTensor;
    Matrix<Double> permTensor;

    stiffTensor.Resize(6,6);
    stiffTensorSteel.Resize(6,6);
    stiffTensorAlu.Resize(6,6);

    piezoTensor.Resize(3,6);
    permTensor.Resize(3,3);

    subdomsMech_ = ptPDE1_->getPDE_subdoms();

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
      
      piezoTensor[1][3]=parameter_[5];  //e_15
      piezoTensor[0][4]=parameter_[5];  //e_15
      piezoTensor[2][0]=parameter_[6]; //e_31
      piezoTensor[2][1]=parameter_[6]; //e_31
      piezoTensor[2][2]=parameter_[7]; // e_33
      
      permTensor[0][0] = parameter_[8]; //eps_11
      permTensor[1][1] = parameter_[8]; //eps_11
      permTensor[2][2] = parameter_[9]; //eps_33

      Double a1, a2, a3;
      a1=a2=a3=0;
      
      if( params->HasValue( "x", "1", "piezo", "polingDirection" ) )
        a1=1;
      
      if( params->HasValue( "y", "1", "piezo", "polingDirection" ) )
        a2=1;
      
      if( params->HasValue( "z", "1", "piezo", "polingDirection" ) )
        a3=1;
      
      if (a1==0&&a2==0&&a3==0)
        a3=1.0;    // if no poling direction is specified, the z-direction is chosen by default 
      
      //     ptMaterial_[2].RotateMaterialMatrix(a1,a2,a3); //?

      ptMaterialPiezo_[0]->SetTensor(piezoTensor,PIEZO_TENSOR,REAL);
      ptMaterialMech_[0]->SetTensor(stiffTensor,MECH_STIFFNESS_TENSOR, REAL);
      ptMaterialElec_[0]->SetTensor(permTensor,ELEC_PERMITTIVITY,REAL);
    }
    else if (subdomsMech_.GetSize()>1){

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
      
      piezoTensor[1][3]=parameter_[5];  //e_15
      piezoTensor[0][4]=parameter_[5];  //e_15
      piezoTensor[2][0]=parameter_[6]; //e_31
      piezoTensor[2][1]=parameter_[6]; //e_31
      piezoTensor[2][2]=parameter_[7]; // e_33
      
      permTensor[0][0] = parameter_[8]; //eps_11
      permTensor[1][1] = parameter_[8]; //eps_11
      permTensor[2][2] = parameter_[9]; //eps_33

      // par[10] - nu
      // par[11] - lambda
     
      // steel
      Double lambda2nu=parameter_[11]+2*parameter_[10];

      stiffTensorSteel[0][0] = lambda2nu; //c_11 
      stiffTensorSteel[1][1] = lambda2nu; //c_11
      stiffTensorSteel[2][2] = lambda2nu; //c_33
      stiffTensorSteel[0][1] = parameter_[11]; //c_12
      stiffTensorSteel[1][0] = parameter_[11]; //c_12
      stiffTensorSteel[0][2] = parameter_[11]; //c_13
      stiffTensorSteel[2][0] = parameter_[11]; //c_13
      stiffTensorSteel[1][2] = parameter_[11]; //c_13
      stiffTensorSteel[2][1] = parameter_[11]; //c_13
      stiffTensorSteel[3][3] = parameter_[10]; //c_44
      stiffTensorSteel[4][4] = parameter_[10]; //c_44
      stiffTensorSteel[5][5] = parameter_[10]; //c_66

      // Aluminium
      lambda2nu=parameter_[13]+2*parameter_[12];

      stiffTensorAlu[0][0] = lambda2nu; //c_11 
      stiffTensorAlu[1][1] = lambda2nu; //c_11
      stiffTensorAlu[2][2] = lambda2nu; //c_33
      stiffTensorAlu[0][1] = parameter_[13]; //c_12
      stiffTensorAlu[1][0] = parameter_[13]; //c_12
      stiffTensorAlu[0][2] = parameter_[13]; //c_13
      stiffTensorAlu[2][0] = parameter_[13]; //c_13
      stiffTensorAlu[1][2] = parameter_[13]; //c_13
      stiffTensorAlu[2][1] = parameter_[13]; //c_13
      stiffTensorAlu[3][3] = parameter_[12]; //c_44
      stiffTensorAlu[4][4] = parameter_[12]; //c_44
      stiffTensorAlu[5][5] = parameter_[12]; //c_66

    
      // alu
      ptMaterialMech_[1]->SetTensor(stiffTensorAlu,MECH_STIFFNESS_TENSOR, REAL);
      //piezo
      ptMaterialPiezo_[2]->SetTensor(piezoTensor,PIEZO_TENSOR,REAL);
      ptMaterialMech_[2]->SetTensor(stiffTensor,MECH_STIFFNESS_TENSOR, REAL);
      ptMaterialElec_[2]->SetTensor(permTensor,ELEC_PERMITTIVITY,REAL);
      ptMaterialPiezo_[3]->SetTensor(piezoTensor,PIEZO_TENSOR,REAL);
      ptMaterialMech_[3]->SetTensor(stiffTensor,MECH_STIFFNESS_TENSOR, REAL);
      ptMaterialElec_[3]->SetTensor(permTensor,ELEC_PERMITTIVITY,REAL);
      //Steel
      ptMaterialMech_[4]->SetTensor(stiffTensorSteel,MECH_STIFFNESS_TENSOR, REAL);
      ptMaterialMech_[5]->SetTensor(stiffTensorSteel,MECH_STIFFNESS_TENSOR, REAL);

    }
    
   
  } // end updateMaterialData

  void piezoParamIdent::updateComplexMaterialData(Vector<Double> & parameterC_){
    ENTER_FCN("piezoParamIdent::updateComplexMaterialData");    
    
    Matrix<Double> stiffTensorC, stiffTensorAluC, stiffTensorSteelC;
    Matrix<Double> piezoTensorC;
    Matrix<Double> permTensorC;

    stiffTensorAluC.Resize(6,6);
    stiffTensorSteelC.Resize(6,6);
    stiffTensorC.Resize(6,6);

    piezoTensorC.Resize(3,6);
    permTensorC.Resize(3,3);

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
      
      piezoTensorC[1][3]=parameterC_[5];  //e_15
      piezoTensorC[0][4]=parameterC_[5];  //e_15
      piezoTensorC[2][0]=parameterC_[6]; //e_31
      piezoTensorC[2][1]=parameterC_[6]; //e_31
      piezoTensorC[2][2]=parameterC_[7]; // e_33
      
      permTensorC[0][0] = parameterC_[8]; //eps_11
      permTensorC[1][1] = parameterC_[8]; //eps_11
      permTensorC[2][2] = parameterC_[9]; //eps_33
  

      ptMaterialPiezo_[0]->SetTensor(piezoTensorC,PIEZO_TENSOR,IMAG);
      ptMaterialMech_[0]->SetTensor(stiffTensorC,MECH_STIFFNESS_TENSOR, IMAG);
      ptMaterialElec_[0]->SetTensor(permTensorC,ELEC_PERMITTIVITY,IMAG);

    }
    else if (subdomsMech_.GetSize()>1){

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
      
      piezoTensorC[1][3]=parameterC_[5];  //e_15
      piezoTensorC[0][4]=parameterC_[5];  //e_15
      piezoTensorC[2][0]=parameterC_[6]; //e_31
      piezoTensorC[2][1]=parameterC_[6]; //e_31
      piezoTensorC[2][2]=parameterC_[7]; // e_33
      
      permTensorC[0][0] = parameterC_[8]; //eps_11
      permTensorC[1][1] = parameterC_[8]; //eps_11
      permTensorC[2][2] = parameterC_[9]; //eps_33
      //       // par[10] - nu
      //       // par[11] - lambda
     
      //       // steel
      Double lambda2nuC=parameterC_[11]+2*parameterC_[10];

      stiffTensorSteelC[0][0] = lambda2nuC; //c_11 
      stiffTensorSteelC[1][1] = lambda2nuC; //c_11
      stiffTensorSteelC[2][2] = lambda2nuC; //c_33
      stiffTensorSteelC[0][1] = parameterC_[11]; //c_12
      stiffTensorSteelC[1][0] = parameterC_[11]; //c_12
      stiffTensorSteelC[0][2] = parameterC_[11]; //c_13
      stiffTensorSteelC[2][0] = parameterC_[11]; //c_13
      stiffTensorSteelC[1][2] = parameterC_[11]; //c_13
      stiffTensorSteelC[2][1] = parameterC_[11]; //c_13
      stiffTensorSteelC[3][3] = parameterC_[10]; //c_44
      stiffTensorSteelC[4][4] = parameterC_[10]; //c_44
      stiffTensorSteelC[5][5] = parameterC_[10]; //c_66

      ///       // Aluminium
        lambda2nuC=parameterC_[13]+2*parameterC_[12];
        
        stiffTensorAluC[0][0] = lambda2nuC; //c_11 
        stiffTensorAluC[1][1] = lambda2nuC; //c_11
        stiffTensorAluC[2][2] = lambda2nuC; //c_33
        stiffTensorAluC[0][1] = parameterC_[13]; //c_12
        stiffTensorAluC[1][0] = parameterC_[13]; //c_12
        stiffTensorAluC[0][2] = parameterC_[13]; //c_13
        stiffTensorAluC[2][0] = parameterC_[13]; //c_13
        stiffTensorAluC[1][2] = parameterC_[13]; //c_13
        stiffTensorAluC[2][1] = parameterC_[13]; //c_13
        stiffTensorAluC[3][3] = parameterC_[12]; //c_44
        stiffTensorAluC[4][4] = parameterC_[12]; //c_44
        stiffTensorAluC[5][5] = parameterC_[12]; //c_66

  //       // alu
        ptMaterialMech_[1]->SetTensor(stiffTensorAluC,MECH_STIFFNESS_TENSOR,IMAG);
        ptMaterialMech_[1]->GetTensor(stiffTensorAluC,MECH_STIFFNESS_TENSOR,IMAG,FULL);
        ptMaterialMech_[1]->GetTensor(stiffTensorAluC,MECH_STIFFNESS_TENSOR,REAL,FULL);

        //       //piezo
        ptMaterialPiezo_[2]->SetTensor(piezoTensorC,PIEZO_TENSOR,IMAG);
        ptMaterialMech_[2]->SetTensor(stiffTensorC,MECH_STIFFNESS_TENSOR,IMAG);
        ptMaterialElec_[2]->SetTensor(permTensorC,ELEC_PERMITTIVITY,IMAG);
        ptMaterialPiezo_[3]->SetTensor(piezoTensorC,PIEZO_TENSOR,IMAG);
        ptMaterialMech_[3]->SetTensor(stiffTensorC,MECH_STIFFNESS_TENSOR, IMAG);
        ptMaterialElec_[3]->SetTensor(permTensorC,ELEC_PERMITTIVITY,IMAG);
        //       //Steel
        ptMaterialMech_[4]->SetTensor(stiffTensorSteelC,MECH_STIFFNESS_TENSOR, IMAG);
        ptMaterialMech_[5]->SetTensor(stiffTensorSteelC,MECH_STIFFNESS_TENSOR, IMAG);
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
    ENTER_FCN("piezoParamIdent::coumputeScaling");

    Matrix<Double> piezoMat,stiffMat, stiffMatAlu, stiffMatSteel, permMat;
    Matrix<Double> piezoMatC, stiffMatAluC, stiffMatSteelC, stiffMatC, permMatC;
    
    if (subdomsMech_.GetSize()==1){
      ptMaterialPiezo_[0]->GetTensor(piezoMat,PIEZO_TENSOR,REAL,FULL);
      ptMaterialMech_[0]->GetTensor(stiffMat,MECH_STIFFNESS_TENSOR, REAL,FULL);
      ptMaterialElec_[0]->GetTensor(permMat,ELEC_PERMITTIVITY,REAL,FULL);

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
        

      if( params->HasValue( "type", "imagMaterialParameter", "materialDataType" ) ){
        ptMaterialPiezo_[0]->GetTensor(piezoMatC,PIEZO_TENSOR,IMAG,FULL);
        ptMaterialMech_[0]->GetTensor(stiffMatC,MECH_STIFFNESS_TENSOR,IMAG,FULL);
        ptMaterialElec_[0]->GetTensor(permMatC,ELEC_PERMITTIVITY,IMAG,FULL);
        scalingC_[0]=1.0/(stiffMatC[0][0]); 
        scalingC_[1]=1.0/(stiffMatC[2][2]);
        scalingC_[2]=1.0/(stiffMatC[1][0]);
        scalingC_[3]=1.0/(stiffMatC[0][2]);
        scalingC_[4]=1.0/(stiffMatC[3][3]); 
        scalingC_[5]=1.0/(piezoMatC[1][3]);
        scalingC_[6]=std::abs(1.0/((piezoMatC)[2][0]));
        scalingC_[7]=1.0/((piezoMatC)[2][2]);
        scalingC_[8]=1.0/((permMatC)[0][0]); 
        scalingC_[9]=1.0/((permMatC)[2][2]);
      }
    }

    else{
  
      ptMaterialMech_[1]->GetTensor(stiffMatAlu,MECH_STIFFNESS_TENSOR, REAL,FULL);
      ptMaterialPiezo_[2]->GetTensor(piezoMat,PIEZO_TENSOR,REAL,FULL);
      ptMaterialMech_[2]->GetTensor(stiffMat,MECH_STIFFNESS_TENSOR, REAL,FULL);
      ptMaterialElec_[2]->GetTensor(permMat,ELEC_PERMITTIVITY,REAL,FULL);
      // ptMaterialMech_[3]->GetTensor(stiffMat,MECH_STIFFNESS_TENSOR, REAL,FULL);
      ptMaterialMech_[5]->GetTensor(stiffMatSteel,MECH_STIFFNESS_TENSOR, REAL,FULL);

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

      scaling_[10]=1.0/(stiffMatSteel[3][3]); 
      scaling_[11]=1.0/(stiffMatSteel[0][1]); 

      scaling_[12]=1.0/(stiffMatAlu[3][3]); 
      scaling_[13]=1.0/(stiffMatAlu[0][1]);
        
      if( params->HasValue( "type", "imagMaterialParameter", "materialDataType" ) ){

        ptMaterialMech_[1]->GetTensor(stiffMatAluC,MECH_STIFFNESS_TENSOR,IMAG,FULL);
        ptMaterialPiezo_[2]->GetTensor(piezoMatC,PIEZO_TENSOR,IMAG,FULL);
        ptMaterialMech_[2]->GetTensor(stiffMatC,MECH_STIFFNESS_TENSOR,IMAG,FULL);
        ptMaterialElec_[2]->GetTensor(permMatC,ELEC_PERMITTIVITY,IMAG,FULL);
        ptMaterialMech_[5]->GetTensor(stiffMatSteelC,MECH_STIFFNESS_TENSOR, IMAG,FULL);

        scalingC_[0]=1.0/(stiffMatC[0][0]); 
        scalingC_[1]=1.0/(stiffMatC[2][2]);
        scalingC_[2]=1.0/(stiffMatC[1][0]);
        scalingC_[3]=1.0/(stiffMatC[0][2]);
        scalingC_[4]=1.0/(stiffMatC[3][3]); 
        scalingC_[5]=1.0/(piezoMatC[1][3]);
        scalingC_[6]=std::abs(1.0/((piezoMatC)[2][0]));
        scalingC_[7]=1.0/((piezoMatC)[2][2]);
        scalingC_[8]=1.0/((permMatC)[0][0]); 
        scalingC_[9]=1.0/((permMatC)[2][2]);
          
        scalingC_[10]=1.0/(stiffMatSteelC[3][3]); 
        scalingC_[11]=1.0/(stiffMatSteelC[0][1]); 
          
        scalingC_[12]=1.0/(stiffMatAluC[3][3]); 
        scalingC_[13]=1.0/(stiffMatAluC[0][1]);                
      }
    }

}//end compute scaling





} // end namespace CoupledField

 


































































