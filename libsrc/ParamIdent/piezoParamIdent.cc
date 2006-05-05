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
    //  std::string measurements="mess.dat";
    //     mess = new std::ifstream(measurements.c_str(), std::basic_ios<char>::in);
    //     if (!mess){
    //       std::cerr << "\n File measuredData.dat does not exist!" << std::endl;
    //     }

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
    if(mess)
      mess->close();


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
    nrParameter_ = 10;

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

     
    ptDomain_->PrintGrid();
    GetMyPDEs();
    DirectCoupledPDE* ptCoupledPDE =  ptDomain_->GetDirectCoupledPDE();
    std::cerr << "ptCoupledPDE = " << ptCoupledPDE << std::endl;

    Info->StartProgress ("Starting to solve problem", FALSE);

    ptPDE1_=ptDomain_-> GetSinglePDE("mechanic");
    ptPDE2_=ptDomain_-> GetSinglePDE("electrostatic");

    ptPDE_->WriteGeneralPDEdefines();

    actNrParameter=0;
    actNrParameterC=0;

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

    // if driver is not part of multiSequence Driver, get list
    // of pdes which have to be solved and intialize them

    
    StdVector<BasePairCoupling*> ptCoupling = ptCoupledPDE->GetCouplingsObject();
    std::cerr << "Name of coupled PDE:" << ptCoupledPDE->GetName()
              << std::endl;

    ptMaterialMech_  = ptPDE1_->getPDEMaterialData();   // Pointer to mech. MaterialData
    ptMaterialElec_  = ptPDE2_->getPDEMaterialData();   // Pointer to elec. MaterialData
    ptMaterialPiezo_ = ptCoupling[0]->getPDEMaterialData();   // Pointer to piezo MaterialData

    //get the material tensors
    Matrix<Double> piezoMat,stiffMat, permMat;
    Matrix<Complex> piezoMatC, stiffMatC, permMatC;
      
    //    if( params->HasValue( "type", "imagMaterialParameter", "materialDataType" ) ){
      ptMaterialPiezo_[0]->GetTensor(piezoMatC,PIEZO_TENSOR,IMAG,FULL);
      ptMaterialMech_[0]->GetTensor(stiffMatC,MECH_STIFFNESS_TENSOR,IMAG,FULL);
      ptMaterialElec_[0]->GetTensor(permMatC,ELEC_PERMITTIVITY,IMAG,FULL);
      std::cout<< " Mech-Tensor\n" << stiffMatC  << std::endl;
      std::cout<< " Elec-Tensor\n" << permMatC   << std::endl;
      std::cout<< " Piezo-Tensor\n" << piezoMatC << std::endl;
      //    }
      //    else{
      ptMaterialPiezo_[0]->GetTensor(piezoMat,PIEZO_TENSOR,REAL,FULL);
      ptMaterialMech_[0]->GetTensor(stiffMat,MECH_STIFFNESS_TENSOR, REAL,FULL);
      ptMaterialElec_[0]->GetTensor(permMat,ELEC_PERMITTIVITY,REAL,FULL);
      std::cout<< " Mech-Tensor\n" << stiffMat  << std::endl;
      std::cout<< " Elec-Tensor\n" << permMat   << std::endl;
      std::cout<< " Piezo-Tensor\n" << piezoMat << std::endl;
      // }

    updateMaterialData(parameter_);
    updateComplexMaterialData(parameter_);
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
   

    // <<<<<<<<<<<<<<<<<<<<<<<< now we have it ... <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  
    //     matMat = ptMaterial_->GetMatrix();
    //     matMatC = ptMaterial_->GetMatrixC();

    //     std::cout<<"We start the calculation with the following material!"<<std::endl;
    //     std::cout<<*matMat<<std::endl;
    //    std::cout<<matMatC<<std::endl;

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



    Vector<Double> c33history(500);
    Vector<Double> e33history(500);
    Vector<Double> eps33history(500);
    Vector<Double> c33historyC(500);
    Vector<Double> e33historyC(500);
    Vector<Double> eps33historyC(500);

    // if we do not wanna scale ..
    //  for (UInt i=0;i<nrParameter_;i++)
    //scaling_[i]=1.0;
  

    // xxxxxxxxxxxxxxxxxxxxxxx Choose different regularizing solvers here xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx //
    

    if (whichNewtonCG_==7){
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
      std::cout<<"++ Optimal experiment Design"<<std::endl;
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

    //    std::cout<<matMatStart<<std::endl;
    //std::cout<<matMatCStart<<std::endl;

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


  void piezoParamIdent::typeOutSolutionOnScreen(Vector<Complex> & solElecPot_,
                                                Vector<Complex> & solMechDispl_){
    ENTER_FCN("piezoParamIdent::typeOutSolutionOnScreen");
    Double sol_real, sol_imag;
    //    std::cout<<"\nElecPot: Amplitude & Phase:"<<std::endl;
    for(UInt i=0;i<solElecPot_.GetSize();i++){
      //      sol_real=solElecPot_[i].real();
      //      sol_imag=solElecPot_[i].imag();
      //   std::cout << "solElecPot_("<< i<< ")=" << sol_real << " + " << sol_imag <<" i " <<std::endl;
      std::cout<<"ElecPot: Amplitude ("<< i <<") = "<< std::abs(solElecPot_[i])
               << ";  Phase ("<< i <<") = "<< std::arg(solElecPot_[i])*180/PI<<std::endl;
    }
    for(UInt i=0;i<solMechDispl_.GetSize();i++){
      sol_real=solMechDispl_[i].real();
      sol_imag=solMechDispl_[i].imag();
      std::cout<<"\nMechDispl: Real & Imag :"<<std::endl;
      std::cout << "solMechDispl_( " << i<< ")=" << sol_real << " + " << sol_imag <<" i " <<std::endl;
    }
  }// end typeOutSolutionOnSreen



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
    while(allMeasuredData->getline(mDataRow, 265)){
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
        thickness_=atof(helpChar);
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
  } // end read MeasuredData

  //! Updates material data & updates system matrices!!
  void piezoParamIdent::updateMaterialData(Vector<Double> & parameter_){
    ENTER_FCN("piezoParamIdent::updateMaterialData");    
    
    Matrix<Double> stiffTensor;
    Matrix<Double> piezoTensor;
    Matrix<Double> permTensor;

    stiffTensor.Resize(6,6);
    piezoTensor.Resize(3,6);
    permTensor.Resize(3,3);


    subdoms_ = ptPDE1_->getPDE_subdoms();
    
    if (subdoms_.GetSize()==1){
        
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

    }

    ptAssemble_ = ptPDE1_->getPDE_assemble();
    ptAssemble2_ = ptPDE2_->getPDE_assemble();

//     ptAssemble_->SetAlternatingMaterial(TRUE);
//     ptAssemble2_->SetAlternatingMaterial(TRUE);
      
    
    // Consider poling of piezoelectric body
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
   
   
  } // end updateMaterialData

  void piezoParamIdent::updateComplexMaterialData(Vector<Double> & parameterC_){
    ENTER_FCN("piezoParamIdent::updateComplexMaterialData");    
    
    Matrix<Double> stiffTensorC;
    Matrix<Double> piezoTensorC;
    Matrix<Double> permTensorC;

    stiffTensorC.Resize(6,6);
    piezoTensorC.Resize(3,6);
    permTensorC.Resize(3,3);


    subdoms_ = ptPDE1_->getPDE_subdoms();
    
    if (subdoms_.GetSize()==1){

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

    }

    ptMaterialPiezo_[0]->SetTensor(piezoTensorC,PIEZO_TENSOR,IMAG);
    ptMaterialMech_[0]->SetTensor(stiffTensorC,MECH_STIFFNESS_TENSOR, IMAG);
    ptMaterialElec_[0]->SetTensor(permTensorC,ELEC_PERMITTIVITY,IMAG);
   


  } // end updateMaterialData

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

} // end namespace CoupledField

 


































































