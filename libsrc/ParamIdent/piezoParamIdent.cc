#include <iostream>
#include <fstream>
#include <string>
//include "staticdriver.hh"
#include "DataInOut/GMV/outGMV.hh"
#include "CoupledPDE/basecoupledpde.hh"
#include "General/environment.hh"
#include <PDE/basePDE.hh>

#include "piezoParamIdent.hh"
#include "Forms/baseForm.hh"
#include "Utils/vector.hh"
#include "Utils/nodestoresol.hh"
#include "Utils/elemstoresol.hh"
#include "DataInOut/MaterialData.hh"
#include "PDE/timestepping.hh"
#include "Utils/baseelemstoresol.hh"
#include "Driver/singleDriver.hh"
#include "PDE/nodeEQN.hh"

#include "DataInOut/WriteInfo.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"

#include <Domain/elem.hh>

#include <sstream>
#include <iomanip>


#ifdef __sgi
#include <stdarg.h>
#include <stdio.h>
#include <math.h>
#define POW pow
#else
#include <cstdarg>
#include <cstdio>
#include <cmath>
#define POW std::pow
#endif

#include "Utils/tools.hh"
#include <PDE/pdes_header.hh>



//#include "/../OLAS/algsys/basesystem.hh"
//#include "DataInOut/piezoParameterData.hh"



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
 
    Char* measuredData="measuredData.dat";
    allMeasuredData = new std::ifstream(measuredData, std::basic_ios<char>::in);

    if (!allMeasuredData)
      {
	std::cerr << "\n File measuredData.dat does not exist!" << std::endl;
	exit(1);
      }

    std::cout<<"\n Opens impedCurve.dat and piezoLog.dat ... "<<std::endl;

    std::string filename= "imped.dat";

    impedCurve = new std::ofstream(filename.c_str(),std::basic_ios<char>::out);

    if (!impedCurve)
      {
	std::cerr << "\n ImpedanceCurve.dat could not be initialized" << std::endl;
      }

    std::string filenameLog= "piezoLog.dat";

    piezoLog = new std::ofstream(filenameLog.c_str(),std::basic_ios<char>::out);

    if (!piezoLog)
      {
	std::cerr << "\n piezoLog.dat could not be initialized" << std::endl;
      }

    std::string filenameParLog= "parLog.dat";
    parLog = new std::ofstream(filenameParLog.c_str(),std::basic_ios<char>::out);

    if (!parLog)
      {
	std::cerr << "\n piezoLog.dat could not be initialized" << std::endl;
      }

  } // end of constructor

  // destructor
  piezoParamIdent :: ~piezoParamIdent()
  {
    ENTER_FCN( "piezoParamIdent::~piezoParamIdent" );
    allMeasuredData->close();
    impedCurve->close();
    piezoLog->close();
    parLog->close();
  }

  void piezoParamIdent :: SolveProblem() {
    ENTER_FCN( "piezoParamIdent::SolveProblem" );

    Integer highestAssumableNrOfMeasData=25;
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

    //    Double pi = 3.14159265358979;
    Double tau=1.5;
   

    // the following passage reads Data from file measuredData.dat
    // The rows are containing the values of the given frequencies, such as phase and amplitude!
    readMeasuredData(freqs, real, imag, parameter, voltage, nrMeasuredData, thickness, radius, delta);
     

    // std::cout<<whichParameterToUpdate<<std::endl;

    // std::cout<<"\n oben wichParToUp ... unten whichParToUpC"<<std::endl;
    std::cout<<whichParameterToUpdateC<<std::endl;
    // real - entspricht |Z|, Betrag der Impedanz
    // imag - entspricht \phi, gemessener Phasenwinkel


    //Kreisfrequenz oder nicht?!   
    //for (int i=0; i<freqs.GetSize();i++)
    // freqs[i]=2*pi*freqs[i];

    //    std::cout<<"\n Size of piezoElectric Body:"<< thickness << " x " << radius <<std::endl;
    //    std::cout<<"\n Number of measure points: " << nrMeasuredData << " with DataError: " << delta <<  std::endl;

    //Settings for harmonic PDE - Driver
    Integer level=0;
    Boolean reset = TRUE;
    Integer pdenumber  = 0;


    if (! isPartOfSequence_)
      ptdomain_->PrintGrid(level);

    if (PrintGridOnly)
      exit(0);

    if (isPartOfSequence_ == FALSE){     
      GetMyPDEs();
      Info->StartProgress ("Starting to solve problem", FALSE);
    }
    pdes_[0]->WriteGeneralPDEdefines();

    actNrParameter=0;
    actNrParameterC=0;

    // how many parameters are there actually to set?
    for (Integer i=0;i<whichParameterToUpdate.GetSize();i++)
      if (whichParameterToUpdate[i]==1)
	actNrParameter++;

     for (Integer i=0;i<whichParameterToUpdateC.GetSize();i++)
       if (whichParameterToUpdateC[i]==1)
 	actNrParameterC++;
   
    whichParToUpInd.Resize(actNrParameter);
    if (whichNewtonCG==4||whichNewtonCG==6||whichNewtonCG==8)
      whichParToUpIndC.Resize(actNrParameterC);

    Integer intTemp=0;

    for (Integer i=0;i<whichParameterToUpdate.GetSize();i++)
      if (whichParameterToUpdate[i]==1){	
	  whichParToUpInd[intTemp]=i;
	  intTemp++;
      }


    intTemp=0;
    for (Integer i=0;i<whichParameterToUpdateC.GetSize();i++)
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

    MaterialData * ptMaterial=pdes_[0]->getPDEMaterialData();   // Pointer to MaterialData

    Matrix<Double> *matMat = ptMaterial->GetMatrix();
    Matrix<Double> *matMatC = ptMaterial->GetMatrixC();
    parameterIncrement=parameter;

    //Alternate Materialparameter by x percent
    if (TRUE)
      for (Integer i=0;i<nrParameter;i++){
	if(i!=2)
	  if (i==1||i==7||i==9)
	    parameter[i]=parameter[i]-0.05*parameter[i];
 	  else if (i==0||i==5||i==6)
 	    parameter[i]=parameter[i]-0.05*parameter[i];	     
	if (i==1)
	  parameter[i]=parameter[i]+0.1*parameter[i];
	if (i==0)
	  parameter[i]=parameter[i]+0.1*parameter[i];	  
      }

    if (FALSE){
      parameter[0]=1.458491485e+11;
      parameter[1]=1.132373842e+11;
      parameter[2]=1.05e+11;
      parameter[3]=9.320007465e+10;
      parameter[4]=2.322576951e+10;
      parameter[5]=11.14101571;
      parameter[6]=-2.959244238;
      parameter[7]=16.008769;
      parameter[8]=9.479219002e-09;
      parameter[9]=8.098761552e-09;
  }
    

    updateMaterialData(parameter,ptMaterial);
    updateComplexMaterialData(parameterC,ptMaterial);


    //    std::cout<<*matMat<<std::endl;
    //getchar();

    Matrix<Double> matMatStart(9,9); // = ptMaterial->GetMatrix();
    Matrix<Double> matMatCStart(9,9); // = ptMaterial->GetMatrixC();

    for(Integer i=0;i<9;i++)
      for(Integer j=0;j<9;j++){
	matMatStart[i][j]=(*matMat)[i][j];
	matMatCStart[i][j]=(*matMatC)[i][j];
      }


 //    ptBCs = pdes_[0]->getPDE_BCs();                             // Pointer to BCs
//     ptAlgsys = pdes_[0]->getPDE_algsys();                       //Pointer to AlgebraicSystem
//     Integer numElems = pdes_[0]->getPDE_numElems();
//     dofs=pdes_[0]->getPDE_dofspernode();
//     numNodes= pdes_[0]->getPDE_numPDENodes();

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

    //    for(Integer i=0;i<nrParameter;i++)
    // parameterC[i]=0.0;
    //    parameterC[7]=1.0;

//     updateMaterialData(parameter, ptMaterial);
//     updateComplexMaterialData(parameterC, ptMaterial);

    //	pdes_[0]->DefineIntegratorsWithMatInfo(level,ptMaterial); // deletes all Integrators and creates new ones with Material in ptMaterial


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


     
    if (CalcImpedanceCurve == 1){
      Vector<Double> freqsTemp = freqs;
      freqs.Resize(nrfreq);
      Double startFreqTemp;
      startFreqTemp=startfreq;
      Double freqincr=(stopfreq-startfreq)/nrfreq;
      for(Integer i=0;i<nrfreq;i++){
	startFreqTemp+=freqincr;
	freqs[i]=startFreqTemp;
      }
      calcImpedanceCurve();
      freqs = freqsTemp;
      getchar();
    }
   

    // <<<<<<<<<<<<<<<<<<<<<<<< now we have it ... <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


    //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx - MeasuredData - xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    // Which data for the right had side of F(p)=y_hat should be taken? ... either .. or ... !!!!
    
      calc_measuredCharge(freqs, real, imag, y_hat); // out of measurements, values are taken from measuredData.dat
      // calcSyntheticData(y_hat); // Generates synthetic data, i.e. one forward simulation will be performed.

     std::cout<<"\n Measured Data - Set: "<<std::endl;
//     for(int i=0;i<y_hat.GetSize();i++)
//       std::cout<<"y("<<i<<")= "<< y_hat[i]<<"; ";
//     std::cout<<"\n"<<std::endl;

    // some values for typical mechanical displacements:
    if (considerMechDeformation==TRUE){

      Integer spacedim = pdes_[0]->getPDE_spaceDim();

      if (spacedim==2){

	y_hat[11]=Complex(-4.475871e-11,0); 
	y_hat[12]=Complex(-5.002088e-11,0); 
	y_hat[13]=Complex(-5.879827e-11,0); 
	y_hat[14]=Complex(-7.289186e-11,0); 
	y_hat[15]=Complex(-1.067941e-10,0); 
	y_hat[16]=Complex(-2.390268e-10,0);
	y_hat[17]=Complex(-1.460878e-09,0); 
	y_hat[18]=Complex(-1.946208e-10,0); 
	y_hat[19]=Complex(-6.589395e-1,0);
	y_hat[20]=Complex(-5.830401e-11,0); 
	y_hat[21]=Complex(-4.383149e-11,0); 

      }

      else if (spacedim==3){
	y_hat[11]=Complex(3.875446e-07,0.000000e+00);
	y_hat[12]=Complex(4.739544e-07,0.000000e+00); 
	y_hat[13]=Complex(3.883185e-07,0.000000e+00); 
	y_hat[14]=Complex(4.742736e-07,0.000000e+00); 
	y_hat[15]=Complex(3.886104e-07,0.000000e+00); 
	y_hat[16]=Complex(4.753440e-07,0.000000e+00); 
	y_hat[17]=Complex(3.892408e-07,0.000000e+00); 
	y_hat[18]=Complex(4.765859e-07,0.000000e+00); 
	y_hat[19]=Complex(3.906905e-07,0.000000e+00); 
	y_hat[20]=Complex(4.785139e-07,0.000000e+00); 
	y_hat[21]=Complex(3.923972e-07,0.000000e+00);
      }

    }

    //    std::cout<<"\n\n yhat_0"<<y_hat[0];


    // xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx Driver Part xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    // since this driver normally will not be used in a coupled 
    // field simulation, we simply can access the first PDE,
    // as it will be the single one
 
    pdes_[0]->WriteGeneralPDEdefines(); 
   
    Complex misfit;

    alpha_m; 
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
    scalingC[2]=1.0/((*matMatC)[1][0]);
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
    //  for (Integer i=0;i<nrParameter;i++)
    //scaling[i]=1.0;
  

    // xxxxxxxxxxxxxxxxxxxxxxx Choose different regularizing solvers here xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx //
    
    //  NewtonCG();
    //         NewtonCG2();
    if (whichNewtonCG==3){
      Integer nNewtonCG =0;
      //    while (nNewtonCG<100){
	 c33history[nNewtonCG] = parameter[1];
	 e33history[nNewtonCG] = parameter[7];
	 eps33history[nNewtonCG] = parameter[8];
	 std::cout<<"\n Nr: "<< nNewtonCG << ", start next NewtonCG Iteration?"<<std::endl;
	 //	getchar();
	 //arLog <<nNewtonCG <<"  "<< c33history[nNewtonCG]<<"  " <<e33history[nNewtonCG]<<"   " <<eps33history[nNewtonCG]<<std::endl;
         NewtonCG3();
	 nNewtonCG++;
	 //  }
    }
    else if (whichNewtonCG==4){
      Integer nNewtonCG =0;
      while (nNewtonCG<200){
	c33history[nNewtonCG] = parameter[1];
 	e33history[nNewtonCG] = parameter[7];
 	eps33history[nNewtonCG] = parameter[9];
	c33historyC[nNewtonCG] = parameterC[1];
	e33historyC[nNewtonCG] = parameterC[7];
	eps33historyC[nNewtonCG] = parameterC[9];
	 std::cout<<"\n Nr: "<< nNewtonCG << ", start next NewtonCG Iteration?"<<std::endl;
	//	getchar();
	  *parLog <<nNewtonCG <<"  "<< c33history[nNewtonCG]<<"  " << 
	    e33history[nNewtonCG]<<"   " <<eps33history[nNewtonCG] 
		<<"  "<< c33historyC[nNewtonCG]<<"  " <<
	  e33historyC[nNewtonCG]<<"   " <<eps33historyC[nNewtonCG]<<std::endl;

         NewtonCG4();
	 nNewtonCG++;
      }
    }
    else if (whichNewtonCG==5){
      Integer nrNewtonLandweber=0;
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
      Integer nrNewtonLandweber=0;
      while (nrNewtonLandweber<maxNumberNewtonLoops){
	std::cout<<"\n Nr: "<< nrNewtonLandweber << ", start next Newton Landweber?"<<std::endl;
 	c33history[nrNewtonLandweber] = parameter[1];
 	e33history[nrNewtonLandweber] = parameter[7];
 	eps33history[nrNewtonLandweber] = parameter[9];
 	c33historyC[nrNewtonLandweber] = parameterC[1];
 	e33historyC[nrNewtonLandweber] = parameterC[7];
 	eps33historyC[nrNewtonLandweber] = parameterC[9];
// 	//getchar();
 	*parLog <<nrNewtonLandweber <<"  "<< c33history[nrNewtonLandweber]<<"  " <<
	  e33history[nrNewtonLandweber]<<"   " <<eps33history[nrNewtonLandweber] 
		<<"  "<< c33historyC[nrNewtonLandweber]<<"  " <<
	  e33historyC[nrNewtonLandweber]<<"   " <<eps33historyC[nrNewtonLandweber]<<std::endl;
	NewtonLandweberC();
	nrNewtonLandweber++;

    }
  }

    else if (whichNewtonCG==7){
      Integer nrNuMethods=0;
      newtonCounter=0;
      while (nrNuMethods<maxNumberNewtonLoops){
      //      while (nrNuMethods<2){
	std::cout<<"\n Nr: "<< nrNuMethods << ", start next Newton NuMethod?"<<std::endl;
// 	c33history[nrNuMethods] = parameter[1];
// 	e33history[nrNuMethods] = parameter[7];
// 	eps33history[nrNuMethods] = parameter[9];
// 	//getchar();
// 	*parLog <<nrNuMethods <<"  "<< c33history[nrNuMethods]<<"  " <<e33history[nrNuMethods]<<"   " <<eps33history[nrNuMethods]<<std::endl;
	nuMethods();
	nrNuMethods++;
	newtonCounter++;
    }
  }

    else if (whichNewtonCG==8){
      Integer nrNuMethodsC=0;
      while (nrNuMethodsC<maxNumberNewtonLoops){
	std::cout<<"\n Nr: "<< nrNuMethodsC << ", start next Newton NuMethodC?"<<std::endl;
	c33history[nrNuMethodsC] = parameterC[1];
	e33history[nrNuMethodsC] = parameterC[7];
	eps33history[nrNuMethodsC] = parameterC[9];
	//getchar();
	*parLog <<nrNuMethodsC <<"  "<< c33history[nrNuMethodsC]<<"  " <<e33history[nrNuMethodsC]<<"   " <<eps33history[nrNuMethodsC]<<std::endl;
	nuMethodsC();
	nrNuMethodsC++;

    }
  }

    else if (whichNewtonCG==9)
      tichonov();
    else
      std::cout<<"\n There was no valid NewtonCG method specified - see in your measuredData.dat -file "<<std::endl;
    //    tichonov();
    //  NewtonLandweber();
    //  createF(ptMaterial, ptBCs, F_hat,FALSE); // calculates only forward problems over all omegas

    // xxxxxxxxxxxxxxxxxxxxxxx End of choice xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx //


    std::cout<<"\n\n *** FINALLY CALCULATED PARAMETERS *** ... here they are: " <<std::endl;

    for (int i=0;i<parameter.GetSize();i++)
      std::cout<<"par[" << i<<"]="<< parameter[i]<<" + " << parameterC[i]<<"i"<<std::endl;


    //    std::cout<<matMatStart<<std::endl;
    //std::cout<<matMatCStart<<std::endl;

// <<<<<<<<<<<<<< for a hopefully nice imped curve after identification !! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

     
    if (CalcImpedanceCurve == 1){
      getchar();
      Vector<Double> freqsTemp = freqs;
      freqs.Resize(nrfreq);
      Double freqincr=(stopfreq-startfreq)/nrfreq;
      for(Integer i=0;i<nrfreq;i++){
	startfreq+=freqincr;
	freqs[i]=startfreq;
      }
      calcImpedanceCurve();
      freqs = freqsTemp;
    }
    
  }// End solveProblem



  // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXxx - now some methods are following ...

  void piezoParamIdent::calcAbsImped(Complex & charge, Double & freq, Integer & fstep, Boolean typeOut){
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
	  //	std::cout<<"\n Frequency - Impendace - Phase: "<<std::endl;
	  std::cout <<"\n Frequency: "<< freq << ", |Z|: "<< std::abs(impedC) << "; Phase: "<< phase << std::endl;
	}

	*impedCurve <<"\n" << freq << "  " << impedC.real()<<"  " << impedC.imag() << imped << "   " << phase << std::endl;

  }  // end calcAbsImped 

   void piezoParamIdent::norm(Vector<Complex> &  vec, Double & norm, Double & norm2,Vector<Complex> & q_meas){
     ENTER_FCN("piezoParamIdent::norm");

     Vector<Complex> y_comp(nrParameter);
     Vector<Complex> y_temp(nrParameter);

     switch (whichNorm){

     case 1:
       norm = a2norm(vec);
       std::cout<<"\n l2-Norm = "<<norm<<std::endl;
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
       for(Integer i=0;i<nrParameter;i++){
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

     default:
       norm=a2norm(vec);

     }
   } // end norm
  
  

  Double piezoParamIdent::calcEuclidianMatrixNorm(Matrix<Complex> & mat){
    ENTER_FCN("piezoParamIdent::calcEuclidianMatrixNorm");

    Double norm=0.0;
    for (Integer i=0;i<mat.GetSizeRow();i++)
      for (Integer j=0;j<mat.GetSizeCol();j++)
	norm+=std::abs(mat[i][j])*std::abs(mat[i][j]);
    norm=sqrt(norm);
    return norm;

  } // end calcEuclidianMatrixNorm

  void piezoParamIdent::maxAndEuclNorm(Vector<Complex> & vec, Double & maxNorm, Double & euclNorm){
    ENTER_FCN("piezoParamIdent::maxAndEuclNorm");
    Double maxNormTemp = 0.0;
    maxNorm=0.0;
    euclNorm=0.0;

    for (Integer i=0;i<vec.GetSize();i++){
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
    for (Integer i=0;i<vec.GetSize();i++){
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

    for (Integer i=0;i<vec.GetSize();i++){
      maxNormTemp=std::abs(vec[i]);
      Denominator = std::abs(q_meas[i])*std::abs(q_meas[i]);
      if (whichNorm==2){
	//	wNorm = wNorm+((1.0/Denominator)*vec[i]*vec[i]).real(); // this is a good running version!
	wNorm = wNorm+((1.0/Denominator)*std::abs(vec[i])*std::abs(vec[i]));
	//        std::cout<<"\n WeightedResNorm = " << std::abs(vec[i])*std::abs(vec[i])<< std::endl;
      }
      else if (whichNorm==5)
	wNorm = wNorm+((1.0/Denominator)*std::abs(vec[i])*std::abs(vec[i]));

      if (maxNormTemp>maxNorm)
	maxNorm=maxNormTemp;
    }

    //wNorm=std::sqrt(wNorm);
  } // end maxAndWeightedNorm


  void piezoParamIdent::calcNorm2Resid(Vector<Complex> &res, Double & anorm, Integer nrMeasuredData){
    ENTER_FCN("piezoParamIdent::calcNorm2Resid");
    anorm=0.0;
    for (int i=0;i<2*nrMeasuredData;i++){
      anorm+=res[i].real()*res[i].real()+ res[i].imag()*res[i].imag();
      anorm=sqrt(anorm);
    }
  } // end calcNorm2Resid

 Double piezoParamIdent::norm2Real(Vector<Complex> &vec){
    ENTER_FCN("piezoParamIdent::realA2norm");
    Double result=0.0; 
    //      Double real_result;
    for(int i=0;i<vec.GetSize();i++)
      result+=vec[i].real()*vec[i].real();
    result=sqrt(result);
    return result;
  }

  Double piezoParamIdent::realA2norm(Vector<Complex> &vec){
    ENTER_FCN("piezoParamIdent::realA2norm");
    Double result=0.0; 
    Complex resultC = Complex(0.0,0.0);
    //      Double real_result;
    for(int i=0;i<vec.GetSize();i++)
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
    for(int i=0;i<vec.GetSize();i++)
      result+=std::abs(vec[i])*std::abs(vec[i]);
    result=sqrt(result);
    return result;
  }

  Double piezoParamIdent::a2norm(Vector<Double> &vec){
    ENTER_FCN("piezoParamIdent::a2norm");
    Double result=0.0; //Complex(0.0,0.0);
    //      Double real_result;
    for(int i=0;i<vec.GetSize();i++)
      result+=std::abs(vec[i])*std::abs(vec[i]);
    result=sqrt(result);
    return result;
  }


  void piezoParamIdent::measureMechDeformationInZ_Direction(Vector<Complex> & mechDisplacement, Double & Radius, Double & meanValueMechDeformation, int dof){
    ENTER_FCN("piezoParamIdent::measureMechDeformationInZ_Direction");
    meanValueMechDeformation=0.0;

    Integer spacedim = pdes_[0]->getPDE_spaceDim();

    std::list<Integer> bcs_list;
    std::string BCName="ep-top";

    bcs_list=ptBCs->GetNodesLevel(BCName, 0);
    std::list<Integer>::const_iterator it;
    it=bcs_list.begin();
    if (spacedim==3){
      while (it!=bcs_list.end())
	{
	  //	std::cout<<"\n MECHDISPL "<< mechDisplacement[(*it)*(dof-1)-1]<< " & it " << (*it)*(dof-1)-1 << std::endl;
	  //meanValueMechDeformation+=std::abs(mechDisplacement[(*it)*(dofs-1)-1]);
	  meanValueMechDeformation+=mechDisplacement[(*it)*(dof-1)-1].real();
	  it++;
	}
      meanValueMechDeformation=meanValueMechDeformation/(PI*Radius*Radius*mechDisplacement.GetSize()/(dof-1));
    }

    else if (spacedim==2){
      while (it!=bcs_list.end())
	{
	  //std::cout<<"\n MECHDISPL "<< mechDisplacement[(*it)*(dof-1)-1]<< " & it " << (*it)*(dof-1)-1 << std::endl;
	  meanValueMechDeformation+=mechDisplacement[(*it)*(dof-1)-1].real();
	  it++;
	}
      meanValueMechDeformation=meanValueMechDeformation/(Radius*mechDisplacement.GetSize()/(dof-1));
    }
  } // end measureMechDeformationInZDirection


  void piezoParamIdent::typeOutSolutionOnScreen(Vector<Complex> & solElecPot,Vector<Complex> & solMechDispl){
    ENTER_FCN("piezoParamIdent::typeOutSolutionOnScreen");
    Double sol_real, sol_imag;
    std::cout<<"\nElecPot: Amplitude & Phase:"<<std::endl;
    for(int i=0;i<solElecPot.GetSize();i++){
      //      sol_real=solElecPot[i].real();
      //      sol_imag=solElecPot[i].imag();
      //   std::cout << "solElecPot("<< i<< ")=" << sol_real << " + " << sol_imag <<" i " <<std::endl;
            std::cout<<"ElecPot: Amplitude ("<< i <<") = "<< std::abs(solElecPot[i])<< ";  Phase ("<< i <<") = "<< std::arg(solElecPot[i])*180/PI<<std::endl;
    }
    for(int i=0;i<solMechDispl.GetSize();i++){
      sol_real=solMechDispl[i].real();
      sol_imag=solMechDispl[i].imag();
      std::cout<<"\nMechDispl: Real & Imag :"<<std::endl;
      std::cout << "solMechDispl( " << i<< ")=" << sol_real << " + " << sol_imag <<" i " <<std::endl;
    }
  }// end typeOutSolutionOnSreen

  void piezoParamIdent::calcInitialResidual(Vector<Complex> & res, Vector<Complex> & y_hat, Vector<Complex> & PHI_p, Integer fstep, Vector<Complex> & solElecPot, Double & meanValueMechDeformation){
    std::list<Integer> bcs_list;
    bcs_list=ptBCs->GetNodesLevel("ep-top", 0); // for cube3dharmonic; zero, because level=0
    //bcs_list=ptBCs->GetNodesLevel("pot", 0); // for cubexi, zero, because level=0
    std::list<Integer>::const_iterator it;
    it=bcs_list.begin();
    PHI_p[fstep]=solElecPot[*it];
    res[fstep]=y_hat[fstep]-PHI_p[fstep];
    res[y_hat.GetSize()+fstep]=meanValueMechDeformation;
    std::cout << "residual ( " << fstep << ")=" << res[fstep].real() << " + " << res[fstep].imag()<<" i " <<std::endl;
  }


  void piezoParamIdent::createMaterialTensorMatrices(Vector<Double> & parameter, Matrix<Double> & couplingMatrix, Matrix<Double> & dielectricMatrix, Integer spaceDim){
    ENTER_FCN("piezoParamIdent::createMaterialTensorMatrices");
    if (spaceDim==2){ // the rotational symmetric case;  couplingMatrix = e
      couplingMatrix[1][0] = couplingMatrix[1][3] = parameter[6]; //e_31
      couplingMatrix[1][1] = parameter[7];	// e_33
      couplingMatrix[0][2] = parameter[5];	//e_15
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


  void piezoParamIdent::readMeasuredData(Vector<Double> & freqs, Vector<Double> & real, Vector<Double> & imag ,Vector<Double> & parameter, Double & voltage, Integer & nrMeasuredData, Double & thickness, Double & radius,  Double & delta){
    ENTER_FCN( "piezoParamIdent::readMeasuredData" );
    char mDataRow[256], helpChar[64];
    Integer i=0, j=0, k=0;
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
	    for(int l=0;l<=k;l++)
	      helpChar[l]=0;
	    j++; i++; k=0;
	  }
	}
	nrMeasuredData = j;
	//	  std::cout<<"Nr_of_measured_Data in readMeasuredData = "<< nrMeasuredData<<std::endl;
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
	    for(int l=0;l<=k;l++)
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
	    for(int l=0;l<=k;l++)
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
	    for(int l=0;l<=k;l++)
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
	    for(int l=0;l<=k;l++)
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
	    for(int l=0;l<=k;l++)
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
	    for(int l=0;l<=k;l++)
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
	for (int l=0;l<=k;l++)
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
	for (int l=0;l<=k;l++)
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
	for (int l=0;l<=k;l++)
	  helpChar[l]=0;
      }
      else if (mDataRow[0]=='S'){
	i=2; k=0;
	while(mDataRow){
	  if (mDataRow[i]=='/')
	    break;
	  helpChar[k]=mDataRow[i];
	  k++; i++;
	}
	nrfreq=atoi(helpChar);
	for (int l=0;l<=k;l++)
	  helpChar[l]=0;
      }
      else if (mDataRow[0]=='L'){
	i=2; k=0;
	while(mDataRow){
	  if (mDataRow[i]=='/')
	    break;
	  helpChar[k]=mDataRow[i];
	  k++; i++;
	}
	startfreq=atof(helpChar);
	for (int l=0;l<=k;l++)
	  helpChar[l]=0;
      }
      else if (mDataRow[0]=='R'){
	i=2; k=0;
	while(mDataRow){
	  if (mDataRow[i]=='/')
	    break;
	  helpChar[k]=mDataRow[i];
	  k++; i++;
	}
	stopfreq=atof(helpChar);
	for (int l=0;l<=k;l++)
	  helpChar[l]=0;
      }
      else if (mDataRow[0]=='C'){
	i=2; k=0;
	while(mDataRow){
	  if (mDataRow[i]=='/')
	    break;
	  helpChar[k]=mDataRow[i];
	  k++; i++;
	}
	maxNumberInnerLoops=atoi(helpChar); 
	for (int l=0;l<=k;l++)
	  helpChar[l]=0;
      }
      else if (mDataRow[0]=='O'){
	i=2; k=0;
	while(mDataRow){
	  if (mDataRow[i]=='/')
	    break;
	  helpChar[k]=mDataRow[i];
	  k++; i++;
	}
	maxNumberNewtonLoops=atoi(helpChar); 
	for (int l=0;l<=k;l++)
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
	for (int l=0;l<=k;l++)
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
	for (int l=0;l<=k;l++)
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
	for (int l=0;l<=k;l++)
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
	for (int l=0;l<=k;l++)
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
	for (int l=0;l<=k;l++)
	  helpChar[l]=0;
      }
    }
  } // end read MeasuredData

  //! Updates material data & updates system matrices!!
  void piezoParamIdent::updateMaterialData(Vector<Double> & parameter, MaterialData * ptMaterial){
    ENTER_FCN("piezoParamIdent::updateMaterialData");    
    // std::cout<<"updateMaterialData"<<std::endl;
      // std::cout<<parameter<<std::endl;

    for(Integer i=0;i<9;i++)
       for(Integer j=0;j<9;j++)
	 ptMaterial->SetPiezoMatrixData(i,j,0.0);

    ptMaterial->SetPiezoMatrixData(0,0, parameter[0]);
    ptMaterial->SetPiezoMatrixData(1,1, parameter[0]);
    ptMaterial->SetPiezoMatrixData(2,2, parameter[1]);
    ptMaterial->SetPiezoMatrixData(0,1, parameter[2]);
    ptMaterial->SetPiezoMatrixData(1,0, parameter[2]);
    ptMaterial->SetPiezoMatrixData(0,2, parameter[3]);
    ptMaterial->SetPiezoMatrixData(2,0, parameter[3]);
    ptMaterial->SetPiezoMatrixData(1,2, parameter[3]);
    ptMaterial->SetPiezoMatrixData(2,1, parameter[3]);
    ptMaterial->SetPiezoMatrixData(3,3, parameter[4]);
    ptMaterial->SetPiezoMatrixData(4,4, parameter[4]);
    // std::cout<<"updateMaterialData Set Data 44"<<std::endl;
     ptMaterial->SetPiezoMatrixData(5,5, 0.5*(parameter[0]-parameter[2]));
    ptMaterial->SetPiezoMatrixData(6,4, parameter[5]);
    ptMaterial->SetPiezoMatrixData(7,3, parameter[5]);
    ptMaterial->SetPiezoMatrixData(4,6, parameter[5]);
    ptMaterial->SetPiezoMatrixData(3,7, parameter[5]);
    ptMaterial->SetPiezoMatrixData(8,0, parameter[6]);
    ptMaterial->SetPiezoMatrixData(8,1, parameter[6]);
    ptMaterial->SetPiezoMatrixData(0,8, parameter[6]);
    ptMaterial->SetPiezoMatrixData(1,8, parameter[6]);
    ptMaterial->SetPiezoMatrixData(8,2, parameter[7]);
    ptMaterial->SetPiezoMatrixData(2,8, parameter[7]);
    ptMaterial->SetPiezoMatrixData(6,6, parameter[8]);
    ptMaterial->SetPiezoMatrixData(7,7, parameter[8]);
    ptMaterial->SetPiezoMatrixData(8,8, parameter[9]);

    ptAssemble = pdes_[0]->getPDE_assemble();
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
   
   
  } // end updateMaterialData

 void piezoParamIdent::updateComplexMaterialData(Vector<Double> & parameterC, MaterialData * ptMaterial){
    ENTER_FCN("piezoParamIdent::updateComplexMaterialData");    
    //    std::cout<<"updateComplexMaterialData"<<std::endl;

    ptMaterial->SetPiezoMatrixDataC(0,0, parameterC[0]);
    ptMaterial->SetPiezoMatrixDataC(1,1, parameterC[0]);
    ptMaterial->SetPiezoMatrixDataC(2,2, parameterC[1]);
    ptMaterial->SetPiezoMatrixDataC(0,1, parameterC[2]);
    ptMaterial->SetPiezoMatrixDataC(1,0, parameterC[2]);
    ptMaterial->SetPiezoMatrixDataC(0,2, parameterC[3]);
    ptMaterial->SetPiezoMatrixDataC(2,0, parameterC[3]);
    ptMaterial->SetPiezoMatrixDataC(1,2, parameterC[3]);
    ptMaterial->SetPiezoMatrixDataC(2,1, parameterC[3]);
    ptMaterial->SetPiezoMatrixDataC(3,3, parameterC[4]);
    ptMaterial->SetPiezoMatrixDataC(4,4, parameterC[4]);
    ptMaterial->SetPiezoMatrixDataC(5,5, 0.5*(parameterC[0]-parameterC[2]));
    ptMaterial->SetPiezoMatrixDataC(6,4, parameterC[5]);
    ptMaterial->SetPiezoMatrixDataC(7,3, parameterC[5]);
    ptMaterial->SetPiezoMatrixDataC(4,6, parameterC[5]);
    ptMaterial->SetPiezoMatrixDataC(3,7, parameterC[5]);
    ptMaterial->SetPiezoMatrixDataC(8,0, parameterC[6]);
    ptMaterial->SetPiezoMatrixDataC(8,1, parameterC[6]);
    ptMaterial->SetPiezoMatrixDataC(0,8, parameterC[6]);
    ptMaterial->SetPiezoMatrixDataC(1,8, parameterC[6]);
    ptMaterial->SetPiezoMatrixDataC(8,2, parameterC[7]);
    ptMaterial->SetPiezoMatrixDataC(2,8, parameterC[7]);
    ptMaterial->SetPiezoMatrixDataC(6,6, parameterC[8]);
    ptMaterial->SetPiezoMatrixDataC(7,7, parameterC[8]);
    ptMaterial->SetPiezoMatrixDataC(8,8, parameterC[9]);

    ptAssemble = pdes_[0]->getPDE_assemble();
    ptAssemble->SetAlternatingMaterial(TRUE);
    ptAssemble->SetMaterialPointer(ptMaterial);

  } // end updateMaterialData

  void piezoParamIdent::setNewParameterSet(Vector<Double> & par,Vector<Double> &  par_new,Vector<Double> & scaling,Double & theta,Vector<Double> & uStep, Vector<Integer> & whichParameterToUpdate){
    Integer helpInd=0;
    for (Integer i=0;i<nrParameter;i++){
      //      std::cout<<"\n setNewParameterSet " << whichParameterToUpdate[i]<<", ";
      if (whichParameterToUpdate[i]==1){
	par_new[i]=par[i]+(1.0/scaling[i])*theta*uStep[helpInd];
	//	std::cout<<"\n parNew = " << par_new[i]<<", step = " << uStep[helpInd] << std::endl;
	helpInd++;
      }
    }
    std::cout<<"\n-----------------------------"<<std::endl;
    //    getchar();


  } // end setNewParameterSet

} // end namespace CoupledField

 
