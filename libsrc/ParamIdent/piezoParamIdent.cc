#include <iomanip>
#include <PDE/SinglePDE.hh>
#include "Domain/domain.hh"
#include "piezoParamIdent.hh"



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

    ptDomain = adomain;
    ptMyPDE_ = NULL;
    residuumParIdent=1.0;
    tau=1.0;
    StdVector<std::string> pdeList;
      params->GetPDEList( pdeList );

      if (pdeList[0]=="piezo")
        directCoupling=FALSE;
      else 
        directCoupling=TRUE;

    std::string  filenameMeasuredData="measuredData.dat";
    allMeasuredData = new std::ifstream(filenameMeasuredData.c_str(), std::basic_ios<char>::in);
    if (!allMeasuredData){
      std::cerr << "\n File measuredData.dat does not exist!" << std::endl;
    }

    Info->StartProgress( "\n Opening output files \n imped.dat, piezoLog.dat, parLog.dat, mechDispl.dat, parFinal.dat ... " );
    //     //    std::cout<<"\n Opening output files (imped.dat, piezoLog.dat, parLog.dat, mechDispl.dat, parFinal.dat) ... "<<std::endl;
    Info->FinishProgress();

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


    Char* measuremets="mess.dat";
    mess = new std::ifstream(measuremets, std::basic_ios<char>::in);
    if (!mess){
      std::cerr << "\n File measuredData.dat does not exist!" << std::endl;
    }

    // in future, several parameters wwill be taken from the xml - file ...
    StdVector<std::string> keyVec, attrVec, valVec;
    
    attrVec = "tag";
    valVec = driverTag_;
    
    // Get time stepping information from parameter object
    keyVec = "paramIdent", "startFreq";
    params->Get(keyVec, attrVec, valVec, startfreq);
    
    keyVec = "paramIdent", "stopFreq";
    params->Get(keyVec, attrVec, valVec, stopfreq);
    
    keyVec = "paramIdent", "numFreq";
    params->Get(keyVec, attrVec, valVec, nrfreq);

    
    // should we calculate the impedance curve?
    CalcImpedanceCurve = params->IsSet("calcImpedanceCurve",  "paramIdent");

    // should we calculate the mechanical displacement curve at selected node?
    CalcMechDisplCurve = params->IsSet("calcMechDisplCurve",  "paramIdent");

    // at which node should the mechDisplCurve be calculated?
    keyVec="paramIdent", "mechDisplAtNode";
    params->Get(keyVec, attrVec, valVec, mechDisplAtNode);

    // DOF of mechDispl.
    keyVec="paramIdent", "dofOfMechDispl";
    params->Get(keyVec, attrVec, valVec, dofOfMechDispl);


    // further important constants for truncated Newton methods
    keyVec="paramIdent", "maxNrInnerIterations";
    params->Get(keyVec, attrVec, valVec, maxNumberInnerLoops);

    keyVec="paramIdent", "maxNrOuterIterations";
    params->Get(keyVec, attrVec, valVec, maxNumberNewtonLoops);

    keyVec="paramIdent", "artDataNoise";
    params->Get(keyVec, attrVec, valVec, delta);
    

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
    nrParameter = 10;

    parameter.Resize(nrParameter);
    parameterC.Resize(nrParameter);
    whichParameterToUpdate.Resize(nrParameter);
    whichParameterToUpdateC.Resize(nrParameter);
    whichParameterToUpdateRC.Resize(1);

    //     parameterIncrement.Resize(nrParameter);
    //parameterIncrement = parameter;
    omegas.Resize(highestAssumableNrOfMeasData);
    freqs.Resize(highestAssumableNrOfMeasData);
    real.Resize(highestAssumableNrOfMeasData);
    imag.Resize(highestAssumableNrOfMeasData);
    amplitude_phase.Resize(highestAssumableNrOfMeasData);
    F_hat.Resize(highestAssumableNrOfMeasData);

    // the following passage reads Data from file measuredData.dat
    // The rows are containing the values of the given frequencies, such as phase and amplitude!
    readMeasuredData(freqs, real, imag, parameter, voltage, nrMeasuredData, thickness, radius, delta);
     
    if (! isPartOfSequence_)
      ptdomain_->PrintGrid();
    if (isPartOfSequence_ == FALSE){     
      GetMyPDEs();
      //      StdVector<DirectCoupledPDE*> ptDirectCoupledPde_ = ptDomain->GetBasePDE();
      //! cast pointer to BasePDE * to pointer of SinglePDE *
      ptMyPDE_ = dynamic_cast<SinglePDE*>(ptPDE_);
      Info->StartProgress ("Starting to solve problem", FALSE);

      //      ptPDE1_=ptPDE_->pde1_;

    }
    if(directCoupling==TRUE){
      ptPDE1_=ptDomain-> GetSinglePDE("mechanic");
      ptPDE2_=ptDomain-> GetSinglePDE("electrostatic");
      ptPDE_->WriteGeneralPDEdefines();
    }
    else
      ptMyPDE_->WriteGeneralPDEdefines();


    //ptMyPDE_->WriteGeneralPDEdefines();
    //    pdeId_ = ptMyPDE_->GetPDEId();
    //    pdeId_ = ptPDE_->GetPDEId();

    actNrParameter=0;
    actNrParameterC=0;

    // how many parameters are there actually to set?
    for (UInt i=0;i<whichParameterToUpdate.GetSize();i++)
      if (whichParameterToUpdate[i]==1)
        actNrParameter++;

    for (UInt i=0;i<whichParameterToUpdateC.GetSize();i++)
      if (whichParameterToUpdateC[i]==1)
        actNrParameterC++;
    if (actNrParameter!=0)   
      whichParToUpInd.Resize(actNrParameter);

    if (whichNewtonCG==4||whichNewtonCG==6||whichNewtonCG==8||whichNewtonCG==10||whichNewtonCG==12)
      if (actNrParameterC!=0)   
        whichParToUpIndC.Resize(actNrParameterC);

    UInt intTemp=0;

    for (UInt i=0;i<whichParameterToUpdate.GetSize();i++)
      if (whichParameterToUpdate[i]==1){        
        whichParToUpInd[intTemp]=i;
        intTemp++;
      }

    intTemp=0;
    for (UInt i=0;i<whichParameterToUpdateC.GetSize();i++)
      if (whichParameterToUpdateC[i]==1) {
        whichParToUpIndC[intTemp]=i;
        intTemp++;
      }

    whichParameterToUpdateRC=whichParameterToUpdate;
    whichParameterToUpdateRC.InsertVector(whichParameterToUpdateC,10);


    // ************************************************************************
    // Communication with BasePDE ... gets i.G. pointers to objects involved  *
    // ************************************************************************

    // if driver is not part of multiSequence Driver, get list
    // of pdes which have to be solved and intialize them


    if(directCoupling==TRUE)
      ptMaterial=ptPDE1_->getPDEMaterialData();   // Pointer to MaterialData
    else
      ptMaterial=ptMyPDE_->getPDEMaterialData();   // Pointer to MaterialData

    
    Matrix<Double> *matMat = ptMaterial->GetMatrix();	
    Matrix<Double> *matMatC = ptMaterial->GetMatrixC();	
    parameterIncrement=parameter;

     
   //  //Generates ansys code for stack actuator
//     if (FALSE){
//       UInt numOfLayers=200;

//       for(UInt i=148;i<numOfLayers;i++){
// 	*parLog<<"asel,s,loc,z,inactiveZ+"<<i<<"*th_layer-eps,inactiveZ+"<<i<<"*th_layer+eps"<<std::endl;
// 	*parLog<<"cm,layer"<<i+1<<",area"<<std::endl;
// 	*parLog<<"vext,layer"<<i+1<<",,,zero,zero,th_layer"<<std::endl;
// 	*parLog<<"allsel"<<std::endl;
//       }

//       for(UInt i=148;i<numOfLayers;i=i+2){
// 	*parLog<<"vsel,a,loc,z,inactiveZ+"<<i+3.5<<"*th_layer-eps,inactiveZ+"<<i+3.5<<"*th_layer+eps"<<std::endl;
//       }

//       for(UInt i=148;i<numOfLayers;i=i+2){
// 	*parLog<<"vsel,a,loc,z,inactiveZ+"<<i+2.5<<"*th_layer-eps,inactiveZ+"<<i+2.5<<"*th_layer+eps"<<std::endl;
//       }


//       for(UInt i=148;i<numOfLayers;i=i+2){
// 	*parLog<<"asel,a,loc,z,inactiveZ+"<<i+3<<"*th_layer-eps,inactiveZ+"<<i+3<<"*th_layer+eps"<<std::endl;
//       }

//       for(UInt i=148;i<numOfLayers;i=i+2){
// 	*parLog<<"asel,a,loc,z,inactiveZ+"<<i+2<<"*th_layer-eps,inactiveZ+"<<i+2<<"*th_layer+eps"<<std::endl;
//       }

//       std::cout<< " Stack Commands written in parLog.dat ..." <<std::endl;
//       getchar();
//     }
    

    updateMaterialData(parameter,ptMaterial);
    updateComplexMaterialData(parameterC,ptMaterial);
  

    Matrix<Double> matMatStart(10,10); // = ptMaterial->GetMatrix();
    Matrix<Double> matMatCStart(10,10); // = ptMaterial->GetMatrixC();

    for(UInt i=0;i<9;i++)
      for(UInt j=0;j<9;j++){
        matMatStart[i][j]=(*matMat)[i][j];
        matMatCStart[i][j]=(*matMatC)[i][j];
      }


    //xxxxxxxxxxxxxxxx Initialize and resize all matrices and vectors involved xxxxxxxxxx


    freqs.Part(0,nrMeasuredData);
    real.Part(0,nrMeasuredData);
    imag.Part(0,nrMeasuredData);
    amplitude_phase.Part(0,nrMeasuredData);
 
    y_hat.Resize(2*nrMeasuredData);
    s_0.Resize(actNrParameter+actNrParameterC);    
    //    bas.Resize(nrParameter);
    res_NE_new.Resize(actNrParameter+actNrParameterC);
    res_NE.Resize(actNrParameter+actNrParameterC);
    lin_res.Resize(2*nrMeasuredData);
    res.Resize(2*nrMeasuredData);
    bas_bar.Resize(2*nrMeasuredData);
    s.Resize(actNrParameter+actNrParameterC);
    scaling.Resize(nrParameter);
    scalingC.Resize(nrParameter);
    F_hat.Resize(2*nrMeasuredData);
    overall_res0.Resize(2*nrMeasuredData);
    parameter_new.Resize(nrParameter);

    //     updateMaterialData(parameter, ptMaterial);
    //     updateComplexMaterialData(parameterC, ptMaterial);


    // ~~~~~~~~~~~~~ modificate the algorithm ~~~~~~~~~~~~~~~

    // If we donnot want to consider the mechanical deformation ...
    considerMechDeformation=FALSE;
    // calculates and determines the ImpedanceCurve before and after identification
    sign=1.0;
    // ~~~~~~~~~~ end of modification part  ~~~~~~~~~~~~~~

    if (considerMechDeformation==FALSE){
      y_hat.Resize(nrMeasuredData);
      lin_res.Resize(nrMeasuredData);
      res.Resize(nrMeasuredData);
      bas_bar.Resize(nrMeasuredData);
      s.Resize(nrParameter);
      scaling.Resize(nrParameter);
      F_hat.Resize(nrMeasuredData);
      overall_res0.Resize(nrMeasuredData);      //    nrMeasuredData=1.0/2.0*nrMeasuredData;
      std::cout<<"\n NRMEASURED DATA = " <<nrMeasuredData<<std::endl;
    }

    // <<<<<<<<<<<<<< for a hopefully nice imped curve <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    if (CalcMechDisplCurve == TRUE){
      Vector<Double> freqsTemp = freqs;
      freqs.Resize(nrfreq);
      Double startFreqTemp;
      startFreqTemp=startfreq;
      Double freqincr=(stopfreq-startfreq)/nrfreq;
      for(UInt i=0;i<nrfreq;i++){
        startFreqTemp+=freqincr;
        freqs[i]=startFreqTemp;
      }
      calcMechDisplCurve();
      freqs = freqsTemp;
      std::cout<<"\n Press any key to continue ... "<<std::endl;
      getchar();
    }

     
    if (CalcImpedanceCurve == TRUE){
      Vector<Double> freqsTemp = freqs;
      freqs.Resize(nrfreq);
      Double startFreqTemp;
      startFreqTemp=startfreq;
      Double freqincr=(stopfreq-startfreq)/nrfreq;
      for(UInt i=0;i<nrfreq;i++){
        startFreqTemp+=freqincr;
        freqs[i]=startFreqTemp;
      }
      calcImpedanceCurve();
      freqs = freqsTemp;
      std::cout<<"\n Press any key to continue ... "<<std::endl;
      getchar();
    }
   

    // <<<<<<<<<<<<<<<<<<<<<<<< now we have it ... <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  
    matMat = ptMaterial->GetMatrix();
    matMatC = ptMaterial->GetMatrixC();

    //     std::cout<<"We start the calculation with the following material!"<<std::endl;
    //     std::cout<<*matMat<<std::endl;
    //    std::cout<<matMatC<<std::endl;

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

    scalingC[0]=1.0/((*matMatC)[0][0]); 
    scalingC[1]=1.0/((*matMatC)[2][2]);
    scalingC[2]=10.0/((*matMatC)[1][0]);
    scalingC[3]=1.0/((*matMatC)[0][2]);
    scalingC[4]=1.0/((*matMatC)[3][3]); 
    scalingC[5]=1.0/((*matMatC)[6][4]);
    scalingC[6]=std::abs(1.0/((*matMatC)[8][0]));
    scalingC[7]=1.0/((*matMatC)[8][2]);
    scalingC[8]=1.0/((*matMatC)[6][6]); 
    scalingC[9]=1.0/((*matMatC)[8][8]);

    Vector<Double> c33history(500);
    Vector<Double> e33history(500);
    Vector<Double> eps33history(500);
    Vector<Double> c33historyC(500);
    Vector<Double> e33historyC(500);
    Vector<Double> eps33historyC(500);

    // if we do not wanna scale ..
    //  for (UInt i=0;i<nrParameter;i++)
    //scaling[i]=1.0;
  

    // xxxxxxxxxxxxxxxxxxxxxxx Choose different regularizing solvers here xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx //
    
    //  NewtonCG();
    //         NewtonCG2();


    if (whichNewtonCG==3){


     nrMeasuredData=5;
     Vector<Double> freqs5;
     freqs5.Resize(5);

     freqs5[0]=2.0e+06;
     freqs5[1]=3.5e+06;
     freqs5[2]=4.5e+06;
     freqs5[3]=5.8e+06;
     freqs5[4]=6.9e+06;
     freqs=freqs5;

     Vector<Double> newFreqs;
     readInMeasurement(newFreqs);

     std::cout<<"newFreqs:"<<std::endl;
     std::cout<<newFreqs<<std::endl;
     calc_measuredCharge(freqs, real, imag, y_hat); // out of new measurements

      UInt nNewtonCG =0;
      //    while (nNewtonCG<100){
      c33history[nNewtonCG] = parameter[1];
      e33history[nNewtonCG] = parameter[7];
      eps33history[nNewtonCG] = parameter[8];
      std::cout<<"\n Nr: "<< nNewtonCG << ", start next NewtonCG Iteration?"<<std::endl;
      //arLog <<nNewtonCG <<"  "<< c33history[nNewtonCG]<<"  " <<e33history[nNewtonCG]<<"   " <<eps33history[nNewtonCG]<<std::endl;
      NewtonCG3();
      nNewtonCG++;
      //  }
    }
    else if (whichNewtonCG==4){
      calc_measuredCharge(freqs, real, imag, y_hat); // out of new measurements
      UInt nNewtonCG =0;
      while (nNewtonCG<200){
        c33history[nNewtonCG] = parameter[1];
        e33history[nNewtonCG] = parameter[7];
        eps33history[nNewtonCG] = parameter[9];
        c33historyC[nNewtonCG] = parameterC[1];
        e33historyC[nNewtonCG] = parameterC[7];
        eps33historyC[nNewtonCG] = parameterC[9];
        std::cout<<"\n Nr: "<< nNewtonCG << ", start next NewtonCG Iteration?"<<std::endl;
        //      getchar();
        *parLog <<nNewtonCG <<"  "<< c33history[nNewtonCG]<<"  " << 
          e33history[nNewtonCG]<<"   " <<eps33history[nNewtonCG] 
                <<"  "<< c33historyC[nNewtonCG]<<"  " <<
          e33historyC[nNewtonCG]<<"   " <<eps33historyC[nNewtonCG]<<std::endl;

        NewtonCG4();
        nNewtonCG++;
      }
    }
    else if (whichNewtonCG==5){
      calc_measuredCharge(freqs, real, imag, y_hat); // out of new measurements
      UInt nrNewtonLandweber=0;
      while (nrNewtonLandweber<maxNumberNewtonLoops){
        std::cout<<"\n Nr: "<< nrNewtonLandweber << ", start next Newton Landweber?"<<std::endl;
        c33history[nrNewtonLandweber] = parameter[1];
        e33history[nrNewtonLandweber] = parameter[7];
        eps33history[nrNewtonLandweber] = parameter[9];
        //getchar();
        *parLog <<nrNewtonLandweber <<"  "<< c33history[nrNewtonLandweber]<<"  " <<e33history[nrNewtonLandweber]<<"   " <<eps33history[nrNewtonLandweber]<<std::endl;
        NewtonLandweber();
        nrNewtonLandweber++;

      }
    }

    else if (whichNewtonCG==6){
      calc_measuredCharge(freqs, real, imag, y_hat); // out of new measurements
      UInt nrNewtonLandweber=0;
      while (nrNewtonLandweber<maxNumberNewtonLoops){
        std::cout<<"\n Nr: "<< nrNewtonLandweber << ", start next Newton Landweber?"<<std::endl;
        c33history[nrNewtonLandweber] = parameter[1];
        e33history[nrNewtonLandweber] = parameter[7];
        eps33history[nrNewtonLandweber] = parameter[9];
        c33historyC[nrNewtonLandweber] = parameterC[1];
        e33historyC[nrNewtonLandweber] = parameterC[7];
        eps33historyC[nrNewtonLandweber] = parameterC[9];
        //      //getchar();
        *parLog <<nrNewtonLandweber <<"  "<< c33history[nrNewtonLandweber]<<"  " <<
          e33history[nrNewtonLandweber]<<"   " <<eps33history[nrNewtonLandweber] 
                <<"  "<< c33historyC[nrNewtonLandweber]<<"  " <<
          e33historyC[nrNewtonLandweber]<<"   " <<eps33historyC[nrNewtonLandweber]<<std::endl;
        NewtonLandweberC();
        nrNewtonLandweber++;


      }
    }

    else if (whichNewtonCG==7){
      calc_measuredCharge(freqs, real, imag, y_hat); // out of new measurements
      UInt nrNuMethods=0;
      newtonCounter=0;
      inner_eta=1.0;
      while (nrNuMethods<maxNumberNewtonLoops){
 
        nuMethods();
        nrNuMethods++;

        for (UInt i=0; i<parameter.GetSize(); i++)
          *parFinal<<parameter[i]<<", ";
        *parFinal<<"/"<<std::endl;
        newtonCounter++;
      }

    }

    else if (whichNewtonCG==8){
      calc_measuredCharge(freqs, real, imag, y_hat); // out of new measurements
      UInt nrNuMethodsC=0;
      newtonCounter=0;
      inner_eta=1.0;
      while (nrNuMethodsC<maxNumberNewtonLoops){

        nuMethodsC2();
        //        *parLog <<nrNuMethodsC <<"  "<< parameter[1]<<"  " <<parameterC[1]<<"   " <<parameter[9]<<"  "<<
        //c33history[nrNuMethodsC]<<"  " <<e33history[nrNuMethodsC]<<"   " <<eps33history[nrNuMethodsC]<<std::endl;
        *parLog<<nrNuMethodsC << parameter << parameterC<<std::endl;
        for (UInt i=0; i<parameter.GetSize(); i++)
          *parFinal<<parameter[i]<<", ";
        *parFinal<<"/"<<std::endl;
        for (UInt i=0; i<parameterC.GetSize(); i++)
          *parFinal<<parameterC[i]<<", ";
        *parFinal<<"/"<<std::endl;

        nrNuMethodsC++;
        newtonCounter++;

      }
    }

    else if (whichNewtonCG==9){
      std::cout<<"++ Optimal experiment Design"<<std::endl;
      optimalExpDesign();
    }

    else if (whichNewtonCG==10){
      std::cout<<"++ Optimal experiment Design - Complex parameters"<<std::endl;
      optimalExpDesign();
    }

    else if (whichNewtonCG==11){
      std::cout<<"++ Optimal experiment Design - variable number of frequencies"<<std::endl;
      optimalExpDesignVarNrFreqs();
    }

    else if (whichNewtonCG==12){
      std::cout<<"++ Optimal experiment Design - variable number of frequencies"<<std::endl;
      leastSquare();
    }



    else
      std::cout<<"\n There was no valid NewtonCG method specified - see in your measuredData.dat -file "<<std::endl;
    //    tichonov();
    //  NewtonLandweber();
    //  createF(ptMaterial, F_hat,FALSE); // calculates only forward problems over all omegas

    // xxxxxxxxxxxxxxxxxxxxxxx End of choice xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx //

    if ( maxNumberNewtonLoops!=0){
      std::cout<<"\n\n *** FINALLY CALCULATED PARAMETERS *** ... here they are: " <<std::endl;
      
      *parFinal<<"4) "<<std::endl;
      for (UInt i=0;i<parameter.GetSize();i++){
        std::cout<<"par[" << i<<"]="<< parameter[i]<<" + " << parameterC[i]<<"i"<<std::endl;
        *parFinal<<parameter[i]<<", ";
      }
      *parFinal<<"/"<<std::endl;
    }

    //    std::cout<<matMatStart<<std::endl;
    //std::cout<<matMatCStart<<std::endl;

    // <<<<<<<<<<<<<< for a hopefully nice imped curve after identification !! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

     
    if (CalcImpedanceCurve == TRUE && maxNumberNewtonLoops!=0){
      Vector<Double> freqsTemp = freqs;
      freqs.Resize(nrfreq);
      Double freqincr=(stopfreq-startfreq)/nrfreq;
      for(UInt i=0;i<nrfreq;i++){
        startfreq+=freqincr;
        freqs[i]=startfreq;
      }
      calcImpedanceCurve();
      std::cout<<"\n Press any key to continue ... "<<std::endl;
      freqs = freqsTemp;
    }
    
  }// End solveProblem



  // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXxx - now some methods are following ...

  void piezoParamIdent::calcAbsImped(Complex & charge, Double & freq, UInt & fstep, Boolean typeOut){
    Double imped, phase;
    Complex impedC;

    if (!impedCurve)
      std::cerr<<"Error opening 'ImpedCurve.dat' "<<std::endl;

    Complex im=Complex(0.0,1);
    impedC=voltage/(charge*2.0*PI*freq*im);
    // We need the following line for a comparison with CAPA
    //    imped = std::abs(voltage/(charge*2.0*PI*freq*im)); phase = -90 - 180.0/PI*(std::arg(charge));
    // This line makes sense in this routine!

    imped = std::abs(voltage/(charge*2.0*PI*freq*im));
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

    Vector<Complex> y_comp(nrParameter);
    Vector<Complex> y_temp(nrParameter);

    switch (whichNorm){

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
      for(UInt i=0;i<nrParameter;i++){
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

  void piezoParamIdent::maxAndEuclNorm(Vector<Complex> & vec, Double & maxNorm, Double & euclNorm){
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
      logNorm = logNorm + std::abs(std::log(vec[i]*vec[i]));
    }
    //    euclNorm=std::sqrt(euclNorm);
  } // end logNorm



  void piezoParamIdent::maxAndWeightedResNorm(Vector<Complex> & vec, Double & maxNorm, Double & wNorm, Vector<Complex> & q_meas){
    ENTER_FCN("piezoParamIdent::maxAndWeightedResNorm");
    Double maxNormTemp = 0.0;
    maxNorm=0.0;
    wNorm=0.0;
    Double Denominator=0.0;

    for (UInt i=0;i<vec.GetSize();i++){
      maxNormTemp=std::abs(vec[i]);
      Denominator = std::abs(q_meas[i])*std::abs(q_meas[i]);
      if (whichNorm==2){
        //      wNorm = wNorm+((1.0/Denominator)*vec[i]*vec[i]).real(); // this is a good running version!
        wNorm = wNorm+((1.0/Denominator)*std::abs(vec[i])*std::abs(vec[i]));
        //std::cout<<"\n WeightedResNorm = " << std::abs(vec[i])*std::abs(vec[i])<< std::endl;
        //      std::cout<<wNorm<<std::endl;
      }
      else if (whichNorm==5)
        wNorm = wNorm+((1.0/Denominator)*std::abs(vec[i])*std::abs(vec[i]));
      else if (whichNorm==6){
        if (whichNewtonCG!=11){
          std::cout<<"This choice of norm is not valid for your case. Set Norm in measuredData.dat = 2"<<std::endl;
          getchar();
        }
        Double stepWidth=0.0;
        stepWidth=std::abs(0.5*freqs[std::min(i+1,actNrParameter)]-freqs[std::max((Integer)i-1,1)]);
        stepWidth/=1.0e+06;
        //wNorm = wNorm+stepWidth*rhos[i]*((1.0/Denominator)*std::abs(vec[i])*std::abs(vec[i]));
        wNorm = wNorm + omegaDiffVec[i]*rhos[i]*((1.0/Denominator)*std::abs(vec[i])*std::abs(vec[i]));
        //        wNorm = wNorm+rhos[i]*((1.0/Denominator)*std::abs(vec[i])*std::abs(vec[i]));
        //      getchar();
      }
      
      if (maxNormTemp>maxNorm)
        maxNorm=maxNormTemp;
    }

    //wNorm=std::sqrt(wNorm);
  } // end maxAndWeightedNorm


  void piezoParamIdent::calcNorm2Resid(Vector<Complex> &res, Double & anorm, UInt nrMeasuredData){
    ENTER_FCN("piezoParamIdent::calcNorm2Resid");
    anorm=0.0;
    for (UInt i=0;i<2*nrMeasuredData;i++){
      anorm+=res[i].real()*res[i].real()+ res[i].imag()*res[i].imag();
      anorm=sqrt(anorm);
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


  void piezoParamIdent::measureMechDeformationInZ_Direction(Vector<Complex> & mechDisplacement, Double & Radius, Double & meanValueMechDeformation, UInt dof){
    ENTER_FCN("piezoParamIdent::measureMechDeformationInZ_Direction");
    meanValueMechDeformation=0.0;

    UInt spacedim = ptMyPDE_->getPDE_spaceDim();

    StdVector<UInt> bcs_list;
    std::string BCName="ep-top";

    
    ptdomain_->GetGrid()->GetNodesByName(bcs_list,BCName);
    if (spacedim==3){
      for (UInt iNode=0; iNode<bcs_list.GetSize(); iNode++ ) 
        {
          //    std::cout<<"\n MECHDISPL "<< mechDisplacement[(*it)*(dof-1)-1]<< " & it " << (*it)*(dof-1)-1 << std::endl;
          //meanValueMechDeformation+=std::abs(mechDisplacement[(*it)*(dofs-1)-1]);
          meanValueMechDeformation+=mechDisplacement[bcs_list[iNode]*(dof-1)-1].real();
        }
      meanValueMechDeformation=meanValueMechDeformation/(PI*Radius*Radius*mechDisplacement.GetSize()/(dof-1));
    }

    else if (spacedim==2){
      for (UInt iNode=0; iNode<bcs_list.GetSize(); iNode++ ) 
        {
          //std::cout<<"\n MECHDISPL "<< mechDisplacement[(*it)*(dof-1)-1]<< " & it " << (*it)*(dof-1)-1 << std::endl;
          meanValueMechDeformation+=mechDisplacement[bcs_list[iNode]*(dof-1)-1].real();
        }
      meanValueMechDeformation=meanValueMechDeformation/(Radius*mechDisplacement.GetSize()/(dof-1));
    }
  } // end measureMechDeformationInZDirection


  void piezoParamIdent::typeOutSolutionOnScreen(Vector<Complex> & solElecPot,Vector<Complex> & solMechDispl){
    ENTER_FCN("piezoParamIdent::typeOutSolutionOnScreen");
    Double sol_real, sol_imag;
    //    std::cout<<"\nElecPot: Amplitude & Phase:"<<std::endl;
    for(UInt i=0;i<solElecPot.GetSize();i++){
      //      sol_real=solElecPot[i].real();
      //      sol_imag=solElecPot[i].imag();
      //   std::cout << "solElecPot("<< i<< ")=" << sol_real << " + " << sol_imag <<" i " <<std::endl;
      std::cout<<"ElecPot: Amplitude ("<< i <<") = "<< std::abs(solElecPot[i])<< ";  Phase ("<< i <<") = "<< std::arg(solElecPot[i])*180/PI<<std::endl;
    }
    for(UInt i=0;i<solMechDispl.GetSize();i++){
      sol_real=solMechDispl[i].real();
      sol_imag=solMechDispl[i].imag();
      std::cout<<"\nMechDispl: Real & Imag :"<<std::endl;
      std::cout << "solMechDispl( " << i<< ")=" << sol_real << " + " << sol_imag <<" i " <<std::endl;
    }
  }// end typeOutSolutionOnSreen

  void piezoParamIdent::calcInitialResidual(Vector<Complex> & res, Vector<Complex> & y_hat, Vector<Complex> & PHI_p, UInt fstep, Vector<Complex> & solElecPot, Double & meanValueMechDeformation){
    StdVector<UInt> bcs_list;
    ptdomain_->GetGrid()->GetNodesByName(bcs_list,"ep-top"); // for cube3dharmonic; zero, because level=0
    //bcs_list=ptBCs->GetNodesLevel("pot", 0); // for cubexi, zero, because level=0
    PHI_p[fstep]=solElecPot[bcs_list[0]];
    res[fstep]=y_hat[fstep]-PHI_p[fstep];
    res[y_hat.GetSize()+fstep]=meanValueMechDeformation;
    std::cout << "residual ( " << fstep << ")=" << res[fstep].real() << " + " << res[fstep].imag()<<" i " <<std::endl;
  }


  void piezoParamIdent::createMaterialTensorMatrices(Vector<Double> & parameter, Matrix<Double> & couplingMatrix, Matrix<Double> & dielectricMatrix, UInt spaceDim){
    ENTER_FCN("piezoParamIdent::createMaterialTensorMatrices");
    if (spaceDim==2){ // the rotational symmetric case;  couplingMatrix = e
      couplingMatrix[1][0] = couplingMatrix[1][3] = parameter[6]; //e_31
      couplingMatrix[1][1] = parameter[7];      // e_33
      couplingMatrix[0][2] = parameter[5];      //e_15
      dielectricMatrix[0][0] = parameter[8]; // \eps_11
      dielectricMatrix[1][1] = parameter[9]; // \eps_33
    }
    else if (spaceDim==3){
      couplingMatrix[2][0] = couplingMatrix[2][1] = parameter[6]; // e_31
      couplingMatrix[2][2] = parameter[7]; // e_33
      couplingMatrix[1][4] = couplingMatrix[0][5] = parameter[5]; // e_15
      dielectricMatrix[0][0] = dielectricMatrix[1][1] = parameter[8]; // \eps_11
      dielectricMatrix[2][2] = parameter[9]; // \eps_33
    }
  } // end createMaterialTensorMatrix


  void piezoParamIdent::readMeasuredData(Vector<Double> & freqs, Vector<Double> & real, Vector<Double> & imag ,Vector<Double> & parameter, Double & voltage, UInt & nrMeasuredData, Double & thickness, Double & radius,  Double & delta){
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
            freqs[j]=atof(helpChar);
            for(UInt l=0;l<=k;l++)
              helpChar[l]=0;
            j++; i++; k=0;
          }
        }
        nrMeasuredData = j;
        //        std::cout<<"Nr_of_measured_Data in readMeasuredData = "<< nrMeasuredData<<std::endl;
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
            real[j]=atof(helpChar);
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
            imag[j]=atof(helpChar);
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
            parameter[j]=atof(helpChar);
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
            parameterC[j]=atof(helpChar);
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
            whichParameterToUpdate[j]=atoi(helpChar);
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
            whichParameterToUpdateC[j]=atoi(helpChar);
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
        voltage=atof(helpChar);
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
        thickness=atof(helpChar);
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
        CalcImpedanceCurve=atoi(helpChar);
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
        whichNewtonCG=atoi(helpChar); 
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
        whichNorm=atoi(helpChar); 
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
        radius=atof(helpChar);  
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
        delta=atof(helpChar); // delta - Fehlerniveau
        for (UInt l=0;l<=k;l++)
          helpChar[l]=0;
      }
    }
  } // end read MeasuredData

  //! Updates material data & updates system matrices!!
  void piezoParamIdent::updateMaterialData(Vector<Double> & parameter, MaterialData * ptMaterial){
    ENTER_FCN("piezoParamIdent::updateMaterialData");    
    // std::cout<<"updateMaterialData"<<std::endl;
    // std::cout<<parameter<<std::endl;


    for(UInt i=0;i<9;i++)
      for(UInt j=0;j<9;j++)
        ptMaterial->SetPiezoMatrixData(i,j,0.0);
    
    StdVector<Integer>subdoms;
    if(directCoupling==TRUE)
      subdoms = ptPDE1_->getPDE_subdoms();
    else 
      subdoms=ptMyPDE_->getPDE_subdoms();
    
    if (subdoms.GetSize()==1){
         
      ptMaterial[0].SetPiezoMatrixData(0,0, parameter[0]);
      ptMaterial[0].SetPiezoMatrixData(1,1, parameter[0]);
      ptMaterial[0].SetPiezoMatrixData(2,2, parameter[1]);
      ptMaterial[0].SetPiezoMatrixData(0,1, parameter[2]);
      ptMaterial[0].SetPiezoMatrixData(1,0, parameter[2]);
      ptMaterial[0].SetPiezoMatrixData(0,2, parameter[3]);
      ptMaterial[0].SetPiezoMatrixData(2,0, parameter[3]);
      ptMaterial[0].SetPiezoMatrixData(1,2, parameter[3]);
      ptMaterial[0].SetPiezoMatrixData(2,1, parameter[3]);
      ptMaterial[0].SetPiezoMatrixData(3,3, parameter[4]);
      ptMaterial[0].SetPiezoMatrixData(4,4, parameter[4]);
      // std::cout<<"updateMaterialData Set Data 44"<<std::endl;
      ptMaterial[0].SetPiezoMatrixData(5,5, 0.5*(parameter[0]-parameter[2]));
      ptMaterial[0].SetPiezoMatrixData(6,4, parameter[5]);
      ptMaterial[0].SetPiezoMatrixData(7,3, parameter[5]);
      ptMaterial[0].SetPiezoMatrixData(4,6, parameter[5]);
      ptMaterial[0].SetPiezoMatrixData(3,7, parameter[5]);
      ptMaterial[0].SetPiezoMatrixData(8,0, parameter[6]);
      ptMaterial[0].SetPiezoMatrixData(8,1, parameter[6]);
      ptMaterial[0].SetPiezoMatrixData(0,8, parameter[6]);
      ptMaterial[0].SetPiezoMatrixData(1,8, parameter[6]);
      ptMaterial[0].SetPiezoMatrixData(8,2, parameter[7]);
      ptMaterial[0].SetPiezoMatrixData(2,8, parameter[7]);
      ptMaterial[0].SetPiezoMatrixData(6,6, parameter[8]);
      ptMaterial[0].SetPiezoMatrixData(7,7, parameter[8]);
      ptMaterial[0].SetPiezoMatrixData(8,8, parameter[9]);
    }
    else{

      ptMaterial[0].SetPiezoMatrixData(0,0, parameter[0]);
      ptMaterial[0].SetPiezoMatrixData(1,1, parameter[0]);
      ptMaterial[0].SetPiezoMatrixData(2,2, parameter[1]);
      ptMaterial[0].SetPiezoMatrixData(0,1, parameter[2]);
      ptMaterial[0].SetPiezoMatrixData(1,0, parameter[2]);
      ptMaterial[0].SetPiezoMatrixData(0,2, parameter[3]);
      ptMaterial[0].SetPiezoMatrixData(2,0, parameter[3]);
      ptMaterial[0].SetPiezoMatrixData(1,2, parameter[3]);
      ptMaterial[0].SetPiezoMatrixData(2,1, parameter[3]);
      ptMaterial[0].SetPiezoMatrixData(3,3, parameter[4]);
      ptMaterial[0].SetPiezoMatrixData(4,4, parameter[4]);
      // std::cout<<"updateMaterialData Set Data 44"<<std::endl;
      ptMaterial[0].SetPiezoMatrixData(5,5, 0.5*(parameter[0]-parameter[2]));
      ptMaterial[0].SetPiezoMatrixData(6,4, parameter[5]);
      ptMaterial[0].SetPiezoMatrixData(7,3, parameter[5]);
      ptMaterial[0].SetPiezoMatrixData(4,6, parameter[5]);
      ptMaterial[0].SetPiezoMatrixData(3,7, parameter[5]);
      ptMaterial[0].SetPiezoMatrixData(8,0, parameter[6]);
      ptMaterial[0].SetPiezoMatrixData(8,1, parameter[6]);
      ptMaterial[0].SetPiezoMatrixData(0,8, parameter[6]);
      ptMaterial[0].SetPiezoMatrixData(1,8, parameter[6]);
      ptMaterial[0].SetPiezoMatrixData(8,2, parameter[7]);
      ptMaterial[0].SetPiezoMatrixData(2,8, parameter[7]);
      ptMaterial[0].SetPiezoMatrixData(6,6, parameter[8]);
      ptMaterial[0].SetPiezoMatrixData(7,7, parameter[8]);
      ptMaterial[0].SetPiezoMatrixData(8,8, parameter[9]);
      ptMaterial[1].SetPiezoMatrixData(0,0, parameter[0]);
      ptMaterial[1].SetPiezoMatrixData(1,1, parameter[0]);
      ptMaterial[1].SetPiezoMatrixData(2,2, parameter[1]);
      ptMaterial[1].SetPiezoMatrixData(0,1, parameter[2]);
      ptMaterial[1].SetPiezoMatrixData(1,0, parameter[2]);
      ptMaterial[1].SetPiezoMatrixData(0,2, parameter[3]);
      ptMaterial[1].SetPiezoMatrixData(2,0, parameter[3]);
      ptMaterial[1].SetPiezoMatrixData(1,2, parameter[3]);
      ptMaterial[1].SetPiezoMatrixData(2,1, parameter[3]);
      ptMaterial[1].SetPiezoMatrixData(3,3, parameter[4]);
      ptMaterial[1].SetPiezoMatrixData(4,4, parameter[4]);
      // std::cout<<"updateMaterialData Set Data 44"<<std::endl;
      ptMaterial[1].SetPiezoMatrixData(5,5, 0.5*(parameter[0]-parameter[2]));
      ptMaterial[1].SetPiezoMatrixData(6,4, -parameter[5]);
      ptMaterial[1].SetPiezoMatrixData(7,3, -parameter[5]);
      ptMaterial[1].SetPiezoMatrixData(4,6, -parameter[5]);
      ptMaterial[1].SetPiezoMatrixData(3,7, -parameter[5]);
      ptMaterial[1].SetPiezoMatrixData(8,0, -parameter[6]);
      ptMaterial[1].SetPiezoMatrixData(8,1, -parameter[6]);
      ptMaterial[1].SetPiezoMatrixData(0,8, -parameter[6]);
      ptMaterial[1].SetPiezoMatrixData(1,8, -parameter[6]);
      ptMaterial[1].SetPiezoMatrixData(8,2, -parameter[7]);
      ptMaterial[1].SetPiezoMatrixData(2,8, -parameter[7]);
      ptMaterial[1].SetPiezoMatrixData(6,6, parameter[8]);
      ptMaterial[1].SetPiezoMatrixData(7,7, parameter[8]);
      ptMaterial[1].SetPiezoMatrixData(8,8, parameter[9]);

      ptMaterial[2].SetPiezoMatrixData(0,0, parameter[0]);
      ptMaterial[2].SetPiezoMatrixData(1,1, parameter[0]);
      ptMaterial[2].SetPiezoMatrixData(2,2, parameter[1]);
      ptMaterial[2].SetPiezoMatrixData(0,1, parameter[2]);
      ptMaterial[2].SetPiezoMatrixData(1,0, parameter[2]);
      ptMaterial[2].SetPiezoMatrixData(0,2, parameter[3]);
      ptMaterial[2].SetPiezoMatrixData(2,0, parameter[3]);
      ptMaterial[2].SetPiezoMatrixData(1,2, parameter[3]);
      ptMaterial[2].SetPiezoMatrixData(2,1, parameter[3]);
      ptMaterial[2].SetPiezoMatrixData(3,3, parameter[4]);
      ptMaterial[2].SetPiezoMatrixData(4,4, parameter[4]);
      ptMaterial[2].SetPiezoMatrixData(5,5, 0.5*(parameter[0]-parameter[2]));
      ptMaterial[2].SetPiezoMatrixData(6,4, 0.0);
      ptMaterial[2].SetPiezoMatrixData(7,3, 0.0);
      ptMaterial[2].SetPiezoMatrixData(4,6, 0.0);
      ptMaterial[2].SetPiezoMatrixData(3,7, 0.0);
      ptMaterial[2].SetPiezoMatrixData(8,0, 0.0);
      ptMaterial[2].SetPiezoMatrixData(8,1, 0.0);
      ptMaterial[2].SetPiezoMatrixData(0,8, 0.0);
      ptMaterial[2].SetPiezoMatrixData(1,8, 0.0);
      ptMaterial[2].SetPiezoMatrixData(8,2, 0.0);
      ptMaterial[2].SetPiezoMatrixData(2,8, 0.0);
      ptMaterial[2].SetPiezoMatrixData(6,6, parameter[8]);
      ptMaterial[2].SetPiezoMatrixData(7,7, parameter[8]);
      ptMaterial[2].SetPiezoMatrixData(8,8, parameter[9]);
    }

    if(directCoupling==TRUE)
      ptAssemble = ptPDE1_->getPDE_assemble();
    else
      ptAssemble = ptMyPDE_->getPDE_assemble();
    
    ptAssemble->SetAlternatingMaterial(TRUE);
    ptAssemble->SetMaterialPointer(ptMaterial);
    //   std::cout<< parameter <<std::endl;

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
    
    ptMaterial->RotateMaterialMatrix(a1,a2,a3);
    //     ptMaterial[1].RotateMaterialMatrix(a1,a2,a3);
    //     ptMaterial[2].RotateMaterialMatrix(a1,a2,a3);
   
   
  } // end updateMaterialData

  void piezoParamIdent::updateComplexMaterialData(Vector<Double> & parameterC, MaterialData * ptMaterial){
    ENTER_FCN("piezoParamIdent::updateComplexMaterialData");    
    //    std::cout<<"updateComplexMaterialData"<<std::endl;

    StdVector<Integer>subdoms;
    if(directCoupling==TRUE)
      subdoms = ptPDE1_->getPDE_subdoms();
    else 
      subdoms=ptMyPDE_->getPDE_subdoms();
    
    if (subdoms.GetSize()==1){

      ptMaterial[0].SetPiezoMatrixDataC(0,0, parameterC[0]);
      ptMaterial[0].SetPiezoMatrixDataC(1,1, parameterC[0]);
      ptMaterial[0].SetPiezoMatrixDataC(2,2, parameterC[1]);
      ptMaterial[0].SetPiezoMatrixDataC(0,1, parameterC[2]);
      ptMaterial[0].SetPiezoMatrixDataC(1,0, parameterC[2]);
      ptMaterial[0].SetPiezoMatrixDataC(0,2, parameterC[3]);
      ptMaterial[0].SetPiezoMatrixDataC(2,0, parameterC[3]);
      ptMaterial[0].SetPiezoMatrixDataC(1,2, parameterC[3]);
      ptMaterial[0].SetPiezoMatrixDataC(2,1, parameterC[3]);
      ptMaterial[0].SetPiezoMatrixDataC(3,3, parameterC[4]);
      ptMaterial[0].SetPiezoMatrixDataC(4,4, parameterC[4]);
      // std::cout<<"updateMaterialData Set Data 44"<<std::endl;
      ptMaterial[0].SetPiezoMatrixDataC(5,5, 0.5*(parameterC[0]-parameterC[2]));
      ptMaterial[0].SetPiezoMatrixDataC(6,4, parameterC[5]);
      ptMaterial[0].SetPiezoMatrixDataC(7,3, parameterC[5]);
      ptMaterial[0].SetPiezoMatrixDataC(4,6, parameterC[5]);
      ptMaterial[0].SetPiezoMatrixDataC(3,7, parameterC[5]);
      ptMaterial[0].SetPiezoMatrixDataC(8,0, parameterC[6]);
      ptMaterial[0].SetPiezoMatrixDataC(8,1, parameterC[6]);
      ptMaterial[0].SetPiezoMatrixDataC(0,8, parameterC[6]);
      ptMaterial[0].SetPiezoMatrixDataC(1,8, parameterC[6]);
      ptMaterial[0].SetPiezoMatrixDataC(8,2, parameterC[7]);
      ptMaterial[0].SetPiezoMatrixDataC(2,8, parameterC[7]);
      ptMaterial[0].SetPiezoMatrixDataC(6,6, parameterC[8]);
      ptMaterial[0].SetPiezoMatrixDataC(7,7, parameterC[8]);
      ptMaterial[0].SetPiezoMatrixDataC(8,8, parameterC[9]);
    } else
      {
        ptMaterial[0].SetPiezoMatrixDataC(0,0, parameterC[0]);
        ptMaterial[0].SetPiezoMatrixDataC(1,1, parameterC[0]);
        ptMaterial[0].SetPiezoMatrixDataC(2,2, parameterC[1]);
        ptMaterial[0].SetPiezoMatrixDataC(0,1, parameterC[2]);
        ptMaterial[0].SetPiezoMatrixDataC(1,0, parameterC[2]);
        ptMaterial[0].SetPiezoMatrixDataC(0,2, parameterC[3]);
        ptMaterial[0].SetPiezoMatrixDataC(2,0, parameterC[3]);
        ptMaterial[0].SetPiezoMatrixDataC(1,2, parameterC[3]);
        ptMaterial[0].SetPiezoMatrixDataC(2,1, parameterC[3]);
        ptMaterial[0].SetPiezoMatrixDataC(3,3, parameterC[4]);
        ptMaterial[0].SetPiezoMatrixDataC(4,4, parameterC[4]);
        // std::cout<<"updateMaterialData Set Data 44"<<std::endl;
        ptMaterial[0].SetPiezoMatrixDataC(5,5, 0.5*(parameterC[0]-parameterC[2]));
        ptMaterial[0].SetPiezoMatrixDataC(6,4, parameterC[5]);
        ptMaterial[0].SetPiezoMatrixDataC(7,3, parameterC[5]);
        ptMaterial[0].SetPiezoMatrixDataC(4,6, parameterC[5]);
        ptMaterial[0].SetPiezoMatrixDataC(3,7, parameterC[5]);
        ptMaterial[0].SetPiezoMatrixDataC(8,0, parameterC[6]);
        ptMaterial[0].SetPiezoMatrixDataC(8,1, parameterC[6]);
        ptMaterial[0].SetPiezoMatrixDataC(0,8, parameterC[6]);
        ptMaterial[0].SetPiezoMatrixDataC(1,8, parameterC[6]);
        ptMaterial[0].SetPiezoMatrixDataC(8,2, parameterC[7]);
        ptMaterial[0].SetPiezoMatrixDataC(2,8, parameterC[7]);
        ptMaterial[0].SetPiezoMatrixDataC(6,6, parameterC[8]);
        ptMaterial[0].SetPiezoMatrixDataC(7,7, parameterC[8]);
        ptMaterial[0].SetPiezoMatrixDataC(8,8, parameterC[9]);


        ptMaterial[1].SetPiezoMatrixDataC(0,0, parameterC[0]);
        ptMaterial[1].SetPiezoMatrixDataC(1,1, parameterC[0]);
        ptMaterial[1].SetPiezoMatrixDataC(2,2, parameterC[1]);
        ptMaterial[1].SetPiezoMatrixDataC(0,1, parameterC[2]);
        ptMaterial[1].SetPiezoMatrixDataC(1,0, parameterC[2]);
        ptMaterial[1].SetPiezoMatrixDataC(0,2, parameterC[3]);
        ptMaterial[1].SetPiezoMatrixDataC(2,0, parameterC[3]);
        ptMaterial[1].SetPiezoMatrixDataC(1,2, parameterC[3]);
        ptMaterial[1].SetPiezoMatrixDataC(2,1, parameterC[3]);
        ptMaterial[1].SetPiezoMatrixDataC(3,3, parameterC[4]);
        ptMaterial[1].SetPiezoMatrixDataC(4,4, parameterC[4]);
        // std::cout<<"updateMaterialData Set Data 44"<<std::endl;
        ptMaterial[1].SetPiezoMatrixDataC(5,5, 0.5*(parameterC[0]-parameterC[2]));
        ptMaterial[1].SetPiezoMatrixDataC(6,4, -parameterC[5]);
        ptMaterial[1].SetPiezoMatrixDataC(7,3, -parameterC[5]);
        ptMaterial[1].SetPiezoMatrixDataC(4,6, -parameterC[5]);
        ptMaterial[1].SetPiezoMatrixDataC(3,7, -parameterC[5]);
        ptMaterial[1].SetPiezoMatrixDataC(8,0, -parameterC[6]);
        ptMaterial[1].SetPiezoMatrixDataC(8,1, -parameterC[6]);
        ptMaterial[1].SetPiezoMatrixDataC(0,8, -parameterC[6]);
        ptMaterial[1].SetPiezoMatrixDataC(1,8, -parameterC[6]);
        ptMaterial[1].SetPiezoMatrixDataC(8,2, -parameterC[7]);
        ptMaterial[1].SetPiezoMatrixDataC(2,8, -parameterC[7]);
        ptMaterial[1].SetPiezoMatrixDataC(6,6, parameterC[8]);
        ptMaterial[1].SetPiezoMatrixDataC(7,7, parameterC[8]);
        ptMaterial[1].SetPiezoMatrixDataC(8,8, parameterC[9]);

        ptMaterial[2].SetPiezoMatrixDataC(0,0, parameterC[0]);
        ptMaterial[2].SetPiezoMatrixDataC(1,1, parameterC[0]);
        ptMaterial[2].SetPiezoMatrixDataC(2,2, parameterC[1]);
        ptMaterial[2].SetPiezoMatrixDataC(0,1, parameterC[2]);
        ptMaterial[2].SetPiezoMatrixDataC(1,0, parameterC[2]);
        ptMaterial[2].SetPiezoMatrixDataC(0,2, parameterC[3]);
        ptMaterial[2].SetPiezoMatrixDataC(2,0, parameterC[3]);
        ptMaterial[2].SetPiezoMatrixDataC(1,2, parameterC[3]);
        ptMaterial[2].SetPiezoMatrixDataC(2,1, parameterC[3]);
        ptMaterial[2].SetPiezoMatrixDataC(3,3, parameterC[4]);
        ptMaterial[2].SetPiezoMatrixDataC(4,4, parameterC[4]);
        ptMaterial[2].SetPiezoMatrixDataC(5,5, 0.5*(parameterC[0]-parameterC[2]));
        ptMaterial[2].SetPiezoMatrixDataC(6,4, 0.0);
        ptMaterial[2].SetPiezoMatrixDataC(7,3, 0.0);
        ptMaterial[2].SetPiezoMatrixDataC(4,6, 0.0);
        ptMaterial[2].SetPiezoMatrixDataC(3,7, 0.0);
        ptMaterial[2].SetPiezoMatrixDataC(8,0, 0.0);
        ptMaterial[2].SetPiezoMatrixDataC(8,1, 0.0);
        ptMaterial[2].SetPiezoMatrixDataC(0,8, 0.0);
        ptMaterial[2].SetPiezoMatrixDataC(1,8, 0.0);
        ptMaterial[2].SetPiezoMatrixDataC(8,2, 0.0);
        ptMaterial[2].SetPiezoMatrixDataC(2,8, 0.0);
        ptMaterial[2].SetPiezoMatrixDataC(6,6, parameterC[8]);
        ptMaterial[2].SetPiezoMatrixDataC(7,7, parameterC[8]);
        ptMaterial[2].SetPiezoMatrixDataC(8,8, parameterC[9]);

        //     ptMaterial->SetPiezoMatrixDataC(0,0, parameterC[0]);
        //     ptMaterial->SetPiezoMatrixDataC(1,1, parameterC[0]);
        //     ptMaterial->SetPiezoMatrixDataC(2,2, parameterC[1]);
        //     ptMaterial->SetPiezoMatrixDataC(0,1, parameterC[2]);
        //     ptMaterial->SetPiezoMatrixDataC(1,0, parameterC[2]);
        //     ptMaterial->SetPiezoMatrixDataC(0,2, parameterC[3]);
        //     ptMaterial->SetPiezoMatrixDataC(2,0, parameterC[3]);
        //     ptMaterial->SetPiezoMatrixDataC(1,2, parameterC[3]);
        //     ptMaterial->SetPiezoMatrixDataC(2,1, parameterC[3]);
        //     ptMaterial->SetPiezoMatrixDataC(3,3, parameterC[4]);
        //     ptMaterial->SetPiezoMatrixDataC(4,4, parameterC[4]);
        //     ptMaterial->SetPiezoMatrixDataC(5,5, 0.5*(parameterC[0]-parameterC[2]));
        //     ptMaterial->SetPiezoMatrixDataC(6,4, parameterC[5]);
        //     ptMaterial->SetPiezoMatrixDataC(7,3, parameterC[5]);
        //     ptMaterial->SetPiezoMatrixDataC(4,6, parameterC[5]);
        //     ptMaterial->SetPiezoMatrixDataC(3,7, parameterC[5]);
        //     ptMaterial->SetPiezoMatrixDataC(8,0, parameterC[6]);
        //     ptMaterial->SetPiezoMatrixDataC(8,1, parameterC[6]);
        //     ptMaterial->SetPiezoMatrixDataC(0,8, parameterC[6]);
        //     ptMaterial->SetPiezoMatrixDataC(1,8, parameterC[6]);
        //     ptMaterial->SetPiezoMatrixDataC(8,2, parameterC[7]);
        //     ptMaterial->SetPiezoMatrixDataC(2,8, parameterC[7]);
        //     ptMaterial->SetPiezoMatrixDataC(6,6, parameterC[8]);
        //     ptMaterial->SetPiezoMatrixDataC(7,7, parameterC[8]);
        //     ptMaterial->SetPiezoMatrixDataC(8,8, parameterC[9]);
      }

    if(directCoupling==TRUE)
      ptAssemble = ptPDE1_->getPDE_assemble();
    else
      ptAssemble = ptMyPDE_->getPDE_assemble();
    ptAssemble->SetAlternatingMaterial(TRUE);
    ptAssemble->SetMaterialPointer(ptMaterial);

  } // end updateMaterialData

  void piezoParamIdent::setNewParameterSet(Vector<Double> & par,Vector<Double> &  par_new,Vector<Double> & scaling,Double & theta,Vector<Double> & uStep, Vector<UInt> & whichParameterToUpdate){
    UInt helpInd=0;
    for (UInt i=0;i<nrParameter;i++){
      //      std::cout<<"\n setNewParameterSet " << whichParameterToUpdate[i]<<", ";
      if (whichParameterToUpdate[i]==1){
        par_new[i]=par[i]+(1.0/scaling[i])*theta*uStep[helpInd];
        //      std::cout<<"\n parNew = " << par_new[i]<<", step = " << uStep[helpInd] << std::endl;
        helpInd++;
      }
    }
    std::cout<<"\n-----------------------------"<<std::endl;
    //    getchar();


  } // end setNewParameterSet

} // end namespace CoupledField

 
