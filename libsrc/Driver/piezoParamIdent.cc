#include <iostream>
#include <fstream>
#include <string>
#include "staticdriver.hh"
#include "DataInOut/GMV/outGMV.hh"
#include "CoupledPDE/basecoupledpde.hh"
#include "General/environment.hh"
#include "PDE/basePDE.hh"

#include "piezoParamIdent.hh"
#include "Forms/baseForm.hh"
#include "Utils/vector.hh"
#include "Utils/nodestoresol.hh"
#include "Utils/elemstoresol.hh"
#include "DataInOut/MaterialData.hh"
#include "PDE/timestepping.hh"
#include "Utils/baseelemstoresol.hh"
#include "singleDriver.hh"
#include "PDE/nodeEQN.hh"

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
  piezoParamIdent::piezoParamIdent(Domain * adomain,
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
    else
      std::cerr <<"\n File measuredData is opened to be read" << std::endl;

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
	std::cerr << "piezoLog.dat could not be initialized" << std::endl;
      }



  } // end of constructor


  // destructor
  piezoParamIdent :: ~piezoParamIdent()
  {
    ENTER_FCN( "piezoParamIdent::~piezoParamIdent" );
    allMeasuredData->close();
    std::cout<<"File measuredData.dat is closed" << std::endl;
  }

  void piezoParamIdent :: SolveProblem() {
    ENTER_FCN( "piezoParamIdent::SolveProblem" );

    Integer highestAssumableNrOfMeasData=20;

    parameter.Resize(10);
    whichParameterToUpdate.Resize(10);
    parameterIncrement.Resize(10);
    parameterIncrement = parameter;
    omegas.Resize(highestAssumableNrOfMeasData);
    freqs.Resize(highestAssumableNrOfMeasData);
    real.Resize(highestAssumableNrOfMeasData);
    imag.Resize(highestAssumableNrOfMeasData);
    amplitude_phase.Resize(highestAssumableNrOfMeasData);
    F_hat.Resize(highestAssumableNrOfMeasData);
    nrParameter = 10;
    Double pi = 3.14159265358979;
    Double tau=1.0;
   

    // the following passage reads Data from file measuredData.dat
    // The rows are containing the values of the given frequencies, such as phase and amplitude!
    readMeasuredData(freqs, real, imag, parameter, voltage, nrMeasuredData, thickness, radius, delta);
    // real - entspricht |Z|, Betrag der Impedanz
    // imag - entspricht \phi, gemessener Phasenwinkel


    //Kreisfrequenz oder nicht?!   
    //for (int i=0; i<freqs.GetSize();i++)
    // freqs[i]=2*pi*freqs[i];

    std::cout<<"\n Size of piezoElectric Body:"<< thickness << " x " << radius <<std::endl;
    std::cout<<"Number of measure points: " << nrMeasuredData << " with DataError: " << delta <<  std::endl;
    std::cout<<"Number of parameters: " << nrParameter<< std::endl;

    //Settings for harmonic PDE - Driver
    Integer level=0;
    Boolean reset = TRUE;
    Integer pdenumber  = 0;

    if (!isPartOfSequence_)
      ptdomain_->PrintGrid(level);

    if (PrintGridOnly)
      exit(0);

    if (isPartOfSequence_ == FALSE){     
      GetMyPDEs();
      Info->StartProgress ("Starting to solve problem", FALSE);
    }


    //    std::cout<<"Here begins communication with base PDE" << std::endl;

    // ************************************************************************
    // Communication with BasePDE ... gets i.G. pointers to objects involved  *
    // ************************************************************************

    // if driver is not part of multiSequence Driver, get list
    // of pdes which have to be solved and intialize them


    MaterialData * ptMaterial=pdes_[0]->getPDEMaterialData();   // Pointer to MaterialData
    //   Matrix<Double> * matMatrix =  ptMaterial->GetMatrix();
    ptBCs = pdes_[0]->getPDE_BCs();                             // Pointer to BCs
    ptAlgsys = pdes_[0]->getPDE_algsys();                       //Pointer to AlgebraicSystem
    Integer numElems = pdes_[0]->getPDE_numElems();
    dofs=pdes_[0]->getPDE_dofspernode();
    numNodes= pdes_[0]->getPDE_numPDENodes();

    //xxxxxxxxxxxxxxxx Initialize and resize all matrices and vectors involved xxxxxxxxxx

    freqs.Part(0,nrMeasuredData);
    real.Part(0,nrMeasuredData);
    imag.Part(0,nrMeasuredData);
    amplitude_phase.Part(0,nrMeasuredData);
 
    y_hat.Resize(2*nrMeasuredData);
    s_0.Resize(nrParameter);    
    //    bas.Resize(nrParameter);
    res_NE_new.Resize(nrParameter);
    res_NE.Resize(nrParameter);
    lin_res.Resize(2*nrMeasuredData);
    res.Resize(2*nrMeasuredData);
    bas_bar.Resize(2*nrMeasuredData);
    s.Resize(nrParameter);
    scaling.Resize(nrParameter);
    F_hat.Resize(2*nrMeasuredData);
    overall_res0.Resize(2*nrMeasuredData);
    parameter_new.Resize(nrParameter);


    // ~~~~~~~~~~~~~ modificate the algorithm ~~~~~~~~~~~~~~~

    // If we donnot want to consider the mechanical deformation ...
    considerMechDeformation=FALSE;
    // calculates and determines the ImpedanceCurve before and after identification

    sign=-1.0;

    // ~~~~~~~~~~ end of modification part  ~~~~~~~~~~~~~~


    if (considerMechDeformation==FALSE){
      y_hat.Resize(nrMeasuredData);
      lin_res.Resize(nrMeasuredData);
      res.Resize(nrMeasuredData);
      bas_bar.Resize(nrMeasuredData);
      s.Resize(nrParameter);
      scaling.Resize(nrParameter);
      F_hat.Resize(nrMeasuredData);
      overall_res0.Resize(nrMeasuredData);
      //    nrMeasuredData=1.0/2.0*nrMeasuredData;
      std::cout<<"\n NRMEASURED DATA = " <<nrMeasuredData<<std::endl;
    }

    // <<<<<<<<<<<<<< for a hopefully nice imped curve <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


     
    if (CalcImpedanceCurve == 1){
      Vector<Double> freqsTemp = freqs;
      Integer nrfreq=100;
      freqs.Resize(nrfreq);
      Double startfreq=2.0e+06;
      Double stopfreq=6.0e+06;
      Double freqincr=(stopfreq-startfreq)/nrfreq;
      for(Integer i=0;i<nrfreq;i++){
	startfreq+=freqincr;
	freqs[i]=startfreq;
      }
      calcImpedanceCurve();
      freqs = freqsTemp;
    }
   

    // <<<<<<<<<<<<<<<<<<<<<<<< now we have it ... <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


    //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx - MeasuredData - xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    // Which data for the right had side of F(p)=y_hat should be taken? ... either .. or ... !!!!
    
      calc_measuredCharge(freqs, real, imag, y_hat); // out of measurements, values are taken from measuredData.dat
      // calcSyntheticData(y_hat); // Generates synthetic data, i.e. one forward simulation will be performed.

    std::cout<<"\n"<<std::endl;
    for(int i=0;i<y_hat.GetSize();i++)
      std::cout<<"y("<<i<<")= "<< y_hat[i]<<"; ";
    std::cout<<"\n"<<std::endl;

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

    updateMaterialData(parameter, ptMaterial);         //Writes initial guesses of parameters (read from MeasuredData.dat) to system
    //    std::cout<<parameter<<std::endl;
 
    pdes_[0]->WriteGeneralPDEdefines(); 
   
    Complex misfit;

    alpha_m; 
    Matrix<Double> *matMat = ptMaterial->GetMatrix();
    std::cout<<"\n before scaling"<<std::endl;

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




//     scaling[0]=1.0e-11;    scaling[1]=1.0e-11;    scaling[2]=1.0e-11;    scaling[3]=1.0e-11;    scaling[4]=1.0e-10;
//     scaling[5]=0.1;     scaling[6]=1.0;     scaling[7]=0.1;
//     scaling[8]=1.0e+8; scaling[9]=1.0e+9;

    // if we do not wanna scale ..
    //  for (Integer i=0;i<nrParameter;i++)
    //scaling[i]=1.0;
  

    // xxxxxxxxxxxxxxxxxxxxxxx Choose different regularizing solvers here xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx //
    
    //  NewtonCG();
    //         NewtonCG2();
      NewtonCG3();
    //    tichonov();
    //   NewtonLandweber();
    //  createF(ptMaterial, ptBCs, F_hat); // calculates only forward problems over all omegas

    // xxxxxxxxxxxxxxxxxxxxxxx End of choice xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx //


    std::cout<<"\n\n *** FINALLY CALCULATED PARAMETERS *** ... here they are:"<<std::endl;

    for (int i=0;i<parameter.GetSize();i++)
      std::cout<<"par[" << i<<"]="<< parameter[i]<<";\n";


// <<<<<<<<<<<<<< for a hopefully nice imped curve after identification !! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

     
    if (CalcImpedanceCurve == 1){
      updateMaterialData(parameter,ptMaterial);
      Integer nrfreq=100;
      freqs.Resize(nrfreq);
      Double startfreq=2.0e+06;
      Double stopfreq=6.0e+06;
      Double freqincr=(stopfreq-startfreq)/nrfreq;
      for(Integer i=0;i<nrfreq;i++){
	startfreq+=freqincr;
	freqs[i]=startfreq;
      }
      calcImpedanceCurve();
    }

   
  }// End solveProblem



  // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXxx - now some methods are following ...



//   void piezoParamIdent::calcAbsImped(Complex & charge, Double & freq, Integer & fstep){
//     Double imped, phase;
//     Complex impedC;

//     if (!impedCurve)
//       std::cerr<<"Error opening 'ImpedCurve.dat' "<<std::endl;
//     Complex im=Complex(0.0,1);
//     impedC=voltage/(charge*2.0*PI*freq*im);
//     imped = std::abs(voltage/(charge*2.0*PI*freq*im)); phase = -90 - 180.0/PI*(std::arg(charge));
//     std::cout <<"Impedanz: "<< impedC << "; Frequenz: " << freq <<"; Phase: "<< phase << std::endl;
//     *impedCurve <<"\n" << freq << "  " << impedC.real()<<"  " << impedC.imag() << imped << "   " << phase << std::endl;

//   }  // end calcAbsImped 


  void piezoParamIdent::calcAbsImped(Complex & charge, Double & freq, Integer & fstep){
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

	std::cout <<"Frequency: "<< freq << ", |Z|: "<< std::abs(impedC) << "; Phase: "<< phase << std::endl;

	*impedCurve <<"\n" << freq << "  " << impedC.real()<<"  " << impedC.imag() << imped << "   " << phase << std::endl;

  }  // end calcAbsImped 

  /*  void piezoParamIdent::norm(Vector<Complex> &  vec, Double & norm, Double & 2ndNorm,Double & q_meas){

    switch (whichNorm)

     case l:
  norm = a2norm(vec);
    break;
  case w:
    maxAndWeightedResNorm(vec,maxNorm,wNorm, q_meas)
  */
  

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


  void piezoParamIdent::maxAndWeightedResNorm(Vector<Complex> & vec, Double & maxNorm, Double & wNorm, Vector<Complex> & q_meas){
    ENTER_FCN("piezoParamIdent::maxAndWeightedResNorm");
    Double maxNormTemp = 0.0;
    maxNorm=0.0;
    wNorm=0.0;
    Double Denominator=0.0;

    for (Integer i=0;i<vec.GetSize();i++){
      maxNormTemp=std::abs(vec[i]);
      Denominator = std::abs(q_meas[i])*std::abs(q_meas[i]);
      wNorm = wNorm+((1.0/Denominator)*vec[i]*vec[i]).real();

      if (maxNormTemp>maxNorm)
	maxNorm=maxNormTemp;
    }
    std::cout<<"\n WeightedResNorm = " << wNorm<< std::endl;
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

    /*  Double piezoParamIdent::a2norm(Vector<Complex> &vec){
	ENTER_FCN("piezoParamIdent::a2norm");
	Complex result=Complex(0.0,0.0);
	Double real_result;
	for(int i=0;i<vec.GetSize();i++)
	result+=vec[i]*vec[i];
	result=sqrt(result);
	real_result=result.real();
	///    real_result=sqrt(real_result);
	std::cout<<" \n real_result"<<real_result<<std::endl;
	return real_result;
	}*/
   
    /*Double piezoParamIdent::a2norm(Vector<Complex> &vec){
      ENTER_FCN("piezoParamIdent::a2norm");
      Complex result=Complex(0.0,0.0);
      Double real_result;
      for(int i=0;i<vec.GetSize();i++)
      real_result+=std::abs(vec[i].real())+std::abs(vec[i].imag());
      //    real_result=sqrt(real_result);
      //    real_result=result.real();
      ///    real_result=sqrt(real_result);
      //std::cout<<" \n real_result"<<real_result<<std::endl;
      return real_result;
      }*/

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

  //   void piezoParamIdent::updateRHS(Vector<Complex> & solElecPot, Vector<Complex> & solMechDispl, Double omega){
  //     ENTER_FCN (" piezoParamIdent::updateRHS");
  //     Vector<Double> new_RHS;
  //     new_RHS.Resize(numNodes * dofs);
  //     for(int i=0;i<numNodes;i++){
  //       new_RHS[i]=solElecPot[i].real();
  //       for(int j=0;j<dofs-1;j++){	
  // 	new_RHS[(j+1)*numNodes+i]=solMechDispl[j*numNodes+i].real();
  // 	//	std::cout<<new_RHS[(j+1)*numNodes+i]<<"; ";
  // 	// matVecRHS(IncrementedRHSMatrix,solMechDispl, solElecPot,newRHS);	
  //       }
  //     }

  //     ptAlgsys->UpdateRHS(1,new_RHS.GetPointer());
  //     std::cout<<"RHS was updated"<< std::endl;
  //   } // end update RHS

  //   void piezoParamIdent::updateRHS(Vector<Complex> & RHSsol){
  //     ENTER_FCN("piezoParamIdent::updateRHS");
  //     Vector<Double> temp;
  //     temp.Resize(RHSsol.GetSize());
  //     for(int i=0;i<RHSsol.GetSize();i++){
  //       temp[i]=RHSsol[i].real();
  //       std::cout<<RHSsol[i]<<"; ";
  //     }
  //     //    ptAlgsys->InitRHS();
  //     ptAlgsys->UpdateRHS(1,temp.GetPointer());
  //     std::cout<<"RHS was updated"<< std::endl;

  //   }

  //   void piezoParamIdent::updateRHS2(Vector<Complex> & RHSsol){
  //     ENTER_FCN("piezoParamIdent::updateRHS2");
  //     Vector<Double> temp;
  //     temp.Resize(RHSsol.GetSize());
  //     for(int i=0;i<RHSsol.GetSize();i++){
  //       temp[i]=RHSsol[i].real();
  //       std::cout<<RHSsol[i]<<"; ";
  //     }
  //     //    ptAlgsys->InitRHS();
  //     ptAlgsys->UpdateRHS(1,temp.GetPointer());
  //     std::cout<<"RHS2 was updated"<< std::endl;

  //   }


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
      else if (mDataRow[0]=='N'){
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
    //  std::cout<<"updateMaterialData"<<std::endl;

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

//     StdVector<std::string> pdeNames;
//     params->GetPDEList( pdeNames );
//     StdVector<std::string> tags;
//     tags.Resize(pdes_.GetSize());
//     tags.Init("anyTag");
    //ptDomain->InitPDEs(pdeNames,1,tags);

  } // end updateMaterialData


  void piezoParamIdent::setNewParameterSet(Vector<Double> & parameter,Vector<Double> &  parameter_new,Vector<Double> & scaling,Double & theta,Vector<Complex> & step, Vector<Integer> & whichParameterToUpdate){

    for (Integer i=0;i<nrParameter;i++){
      std::cout<<whichParameterToUpdate[i]<<", ";
      if (whichParameterToUpdate[i]==1)
	parameter_new[i]=parameter[i]+(1.0/scaling[i])*theta*step[i].real(); 
    }
    std::cout<<"\n-----------------------------"<<std::endl;


  } // end setNewParameterSet


} // end namespace CoupledField

 
