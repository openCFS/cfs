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

#include <stdlib.h>
#include <Domain/elem.hh>

#include <sstream>
#include <iomanip>

#include <stdio.h>

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
 
    //Char* measuredData="measuredData.dat";
    allMeasuredData.open("measuredData.dat");

    if (!allMeasuredData.good())
      {
	std::cerr << "File measuredData.dat does not exist!" << std::endl;
	exit(1);
      }
    else
      std::cerr <<"File measuredData is opened to be read" << std::endl;

    std::cout<<"Opens impedCurve.dat ... "<<std::endl;
    std::string filename= "imped.dat";
    impedCurve = new std::ofstream(filename.c_str(),std::basic_ios<char>::out);

    //      Char* char_impedCurve="impedCurve.dat";
    //      impedCurve.open('impedCurve.dat', basic_ios::out);
    //      std::ofstream impedCurve = new std::ofstream("impedCurve.dat",std::basic_ios::out);
    //    std::string filename= "impedCurve.dat";

    //       impedCurve = new std::ofstream(filename.c_str());
    //    impedCurve.open(char_impedCurve,std::out);
    //      std::ofstream impedCurve('impedCurve.dat', basic_ios::out);

    if (!impedCurve)
      {
	std::cerr << "ImpedanceCurve.dat could not be initialized" << std::endl;
      }
    std::cout<<"Opend impedCurve.dat ... "<<std::endl;


  } // end of constructor


  // destructor
  piezoParamIdent :: ~piezoParamIdent()
  {
    ENTER_FCN( "piezoParamIdent::~piezoParamIdent" );
    allMeasuredData.close();
    std::cout<<"File measuredData.dat is closed" << std::endl;
  }

  void piezoParamIdent :: SolveProblem() {
    ENTER_FCN( "piezoParamIdent::SolveProblem" );

    Integer highestAssumableNrOfMeasData=20;

    parameter.Resize(10);
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

    std::cout<<"Size of piezoElectric Body:"<< thickness << " x " << radius <<std::endl;
    std::cout<<"Number of measure points: " << nrMeasuredData << " with DataError: " << delta <<  std::endl;
    std::cout<<"Number of parameters: " << nrParameter<< std::endl;

    freqs.Part(0,nrMeasuredData);
    real.Part(0,nrMeasuredData);
    imag.Part(0,nrMeasuredData);
    amplitude_phase.Part(0,nrMeasuredData);
    y_hat.Resize(2*nrMeasuredData);

    calc_measuredCharge(freqs, real, imag, y_hat);

    for(int i=0;i<nrMeasuredData;i++)
      std::cout<<"; q("<<i<<")= "<< y_hat[i];

    //Ladungen mit exakten Parametern berechnet:
    //   y_hat[0]=(-8.16274e-08,2.91866e-06); // (-6.98873e-08,2.49888e-06); 
    //     y_hat[1]=(-1.07233e-07,2.91787e-06); // (-9.181e-08,2.4982e-06);
    //     y_hat[2]=(-1.94767e-07,2.91337e-06); //(-1.66754e-07,2.49434e-06);
    //     y_hat[3]=(-2.66432e-06,-1.19465e-06); //(-2.28111e-06,-1.02283e-06);
    //     y_hat[4]=(-3.54947e-07,-2.89826e-06); // (-3.03895e-07,-2.4814e-06); 
    //     y_hat[5]=(-1.20177e-06,2.66115e-06); // (-1.02892e-06,2.27839e-06); 
    //     y_hat[6]=(-1.06335e-07,2.918e-06); // (-9.10406e-08,2.4983e-06); 
    //     y_hat[7]=(-7.99403e-08,2.91885e-06); // (-6.84422e-08,2.49902e-06); 
    //     y_hat[8]=(1.03852e-09,0); 
    y_hat[9]=(8.24177e-10,0); 
    y_hat[10]=(6.6726e-10,0); 
    y_hat[11]=(5.51238e-10,0); 
    y_hat[12]=(4.63042e-10,0); 
    y_hat[13]=(3.94437e-10,0); 
    y_hat[14]=(3.40021e-10,0); 
    y_hat[15]=(2.96137e-10,0);

    /* for (int i=0; i<nrMeasuredData;i++){
       y_hat[i]=std::polar(-real[i], imag[i]);
       y_hat[i]=voltage*y_hat[i];
       //y_hat[i+nrMeasuredData]=       
       std::cout <<"Phase&Amplitude (polar): " <<  y_hat[i].real() <<" + "<< y_hat[i].imag() << "i" <<std::endl;
       }*/

    // Settings for Newton-CG - routine
    Integer nrIterations=0;
    Double tol=0.001, tolCG=0.001, err=1.0, res_norm;

    //Settings for harmonic PDE - Driver
    Integer level=0;
    Boolean reset = TRUE;
    Integer pdenumber  = 0;
    if (!isPartOfSequence_)
      ptdomain_->PrintGrid(level);

    if (PrintGridOnly)
      exit(0);

    //    std::cout<<"Here begins communication with base PDE" << std::endl;

    // ************************************************************************
    // Communication with BasePDE ... gets i.G. pointers to objects involved  *
    // ************************************************************************


    // if driver is not part of multiSequence Driver, get list
    // of pdes which have to be solved and intialize them

    if (isPartOfSequence_ == FALSE){     
      GetMyPDEs();
      Info->StartProgress ("Starting to solve problem", FALSE);
    }

    // since this driver normally will not be used in a coupled 
    // field simulation, we simply can access the first PDE,
    // as it will be the single one
 
    pdes_[0]->WriteGeneralPDEdefines(); 
 
    MaterialData * ptMaterial=pdes_[0]->getPDEMaterialData();   // Pointer to MaterialData
    Matrix<Double> * matMatrix =  ptMaterial->GetMatrix();

    updateMaterialData(parameter, ptMaterial);         //Writes initial guesses of parameters (read from MeasuredData.dat) to system

    ptBCs = pdes_[0]->getPDE_BCs();                             // Pointer to BCs
    ptAlgsys = pdes_[0]->getPDE_algsys();                       //Pointer to AlgebraicSystem
    Integer numElems = pdes_[0]->getPDE_numElems();
    dofs=pdes_[0]->getPDE_dofspernode();
    numNodes= pdes_[0]->getPDE_numPDENodes();
   

    y_hat.Resize(2*nrMeasuredData);
    s_0.Resize(nrParameter);    
    bas.Resize(nrParameter);
    res_NE_new.Resize(nrParameter);
    res_NE.Resize(nrParameter);
    res.Resize(2*nrMeasuredData);
    bas_bar.Resize(2*nrMeasuredData);
    s.Resize(nrParameter);
    res.Resize(2*nrMeasuredData);
    scaling.Resize(nrParameter);
    F_hat.Resize(2*nrMeasuredData);
    overall_res0.Resize(2*nrMeasuredData);

    Complex misfit;

    alpha_m; 
     scaling[0]=1.0e-11;    scaling[1]=1.0e-11;    scaling[2]=1.0e-11;    scaling[3]=1.0e-11;    scaling[4]=1.0e-11;
   scaling[5]=0.1;     scaling[6]=0.1;     scaling[7]=0.1;
   scaling[8]=1.0e+8; scaling[9]=1.0e+8;
   // for(int i=0;i<scaling.GetSize();i++)
   //  scaling[i]=1.0;

    //here begins the Newton-Routine:   

   // we do a bit of scaling here ...
   //   for(int i=0;i<scaling.GetSize();i++)
   //    parameter[i]=parameter[i]*scaling[i];

   // we do a bit of rescaling here ...
   //  for(int i=0;i<scaling.GetSize();i++)
   //    parameter[i]=parameter[i]/scaling[i];


    while (nrIterations<3) {
      nrIterations++;
      std::cout<<"\n NewtonCG ... Newton-Iteration-Nr = "<<nrIterations<<std::endl;

      // Create the Matrices F, F', F*
      createF(ptMaterial, ptBCs, F_hat);
      std::cout<<"Parameter-to-solution-map Matrix is built up!"<<std::endl;
      createJacobiMatrix(ptMaterial, ptBCs, F_hat, parameterIncrement,JacobiMatrix, solElecPot, solMechDispl);
      std::cout<<"JacobiMatrix was created ..."<<std::endl;
      createAdjointJacobiMatrix(parameterIncrement, parameter, JacobiMatrix,solElecPot,solMechDispl,freqs,adjJacobiMatrix);
      std::cout<<"Adjoint JacMatrix was created ..."<<std::endl;

      Double real_misfit, real_misfit_new,eta_max, eta_new, t, a,b,c, theta, theta_min, theta_max, gamma, al;
      Complex misfit;
      Vector<Complex> y_hat_F_hat=y_hat-F_hat;

      real_misfit=a2norm(y_hat_F_hat);
      std::cout<<"\n Vor CG real(misfit) = real(y_hat-F_hat) = "<< real_misfit;

 // we do a bit of scaling here ...
     for(int i=0;i<scaling.GetSize();i++)
       parameter[i]=parameter[i]*scaling[i];

      CG(); // calls CG Routine ... uses F, F', F'*    

      y_hat_F_hat=y_hat-F_hat;
      real_misfit=a2norm(y_hat_F_hat);
  
      //here begins attempt to implement the backtracking Algo .... ;-)

      t=1.0e-4;
      eta_max=0.9;
      eta=0.8;
      theta_min=0.1;
      theta_max=0.5;
      gamma=0.5;
      al=1.5;
      std::cout<<"\n NACH CG real(misfit) = real(y_hat-F_hat) = "<< real_misfit;
      int backtrackIterator=0;
  
      while((real_misfit>(1-t*(1-eta))*real_misfit)&&(eta<eta_max)){
	std::cout<<"backtracking ... "<< backtrackIterator<< std::endl;
	backtrackIterator++;

	// choose theta:
	a=real_misfit*real_misfit;
	b=2*real_misfit;
	c=real_misfit*real_misfit-b-a;
	std::cout<<"c"<<c<<std::endl;
	if (c==0.0){
	  if(b<=0)
	    theta=theta_max;
	  else 
	    theta=theta_min;
	}
	else if (c>0){
	  if (b<-2*c*theta_max)
	    theta = theta_max;
	  else if (b>-2*c*theta_min)
	    theta=theta_min;
	  else
	    theta=-b/(2*c);
	}
	else{
	  if (b<-2*c*theta_min)
	    theta=theta_max;
	  else if (b>-2*c*theta_max)
	    theta=theta_min;
	  else
	    {
	      if ((a+b*theta_min+c*theta_min*theta_min)>=(a+b*theta_max+c*theta_max*theta_max))
		theta=theta_max;
	      else theta=theta_min;
	    }
	} // end else if c>0

	std::cout<<"Choice of theta = " << theta;
	//end choose theta
	for (int i=0;i<s.GetSize();i++){
	 s[i]=theta*s[i];
	// std::cout<<"; s= "<< s[i];
	}
	eta=1-theta*(1-eta);

	for (int i=0;i<parameter.GetSize();i++)
	  parameter[i]=parameter[i]+s[i].real();

	createF(ptMaterial, ptBCs, F_hat);       
	y_hat_F_hat=y_hat-F_hat;
	real_misfit=a2norm(y_hat_F_hat);
	//	real_misfit_new=misfit.real();
	//	std::cout<<"\n real_misfit_new = "<<real_misfit_new;

      } //end while, end backtracking

      //choose eta
      eta_new=gamma*POW((real_misfit_new/real_misfit),al);
      if (gamma*POW(eta,al) > 0.1)
        eta_new=std::max(eta_new,gamma*POW(eta,al));
      if (eta_new>eta_max)
        eta_new=eta_max;
      if(eta_new<=2*(tau*delta)/real_misfit_new)
	eta_new=0.8*(tau*delta)/real_misfit_new;
      //end choose eta
      std::cout<<"Choice of eta ..." << eta_new; 
      eta=eta_new;
      
      for (int i=0;i<parameter.GetSize();i++){
	parameter[i]=parameter[i]+s[i].real();
	std::cout<<"par[" << i<<"]="<< parameter[i]<<"; ";
      }

   // we do a bit of rescaling here ...
     for(int i=0;i<scaling.GetSize();i++)
       parameter[i]=parameter[i]/scaling[i];

      updateMaterialData(parameter,ptMaterial);

    } // end while-Newton

    std::cout<<"\n\n *** FINALLY CALCULATED PARAMETERS *** ... here they are:"<<std::endl;
    for (int i=0;i<parameter.GetSize();i++)
      std::cout<<"par[" << i<<"]="<< parameter[i]<<"; ";

  }// End solveProblem

  // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXxx - now some methods are following ...

  void piezoParamIdent::calcAbsImped(Complex & charge, Double & freq, Integer & fstep){
    Double imped;

    if (!impedCurve)
      std::cerr<<"Error opening 'ImpedCurve.dat' "<<std::endl;
    Complex temp=Complex(0.0,1);
    imped = std::abs(voltage/charge*freq*temp);
    std::cout <<"Betrag der Impedanz:"<< imped << std::endl;
    *impedCurve <<"\n" << imped << "   " << freq<< std::endl;

  }  // end calcAbsImped 


  void piezoParamIdent::calcNorm2Resid(Vector<Complex> &res, Double & anorm, Integer nrMeasuredData){
    anorm=0.0;
    for (int i=0;i<2*nrMeasuredData;i++){
      anorm+=std::abs(res[i])*std::abs(res[i]);
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

  Double piezoParamIdent::a2norm(Vector<Complex> &vec){
      ENTER_FCN("piezoParamIdent::a2norm");
      Complex result=Complex(0.0,0.0);
      Double real_result;
      for(int i=0;i<vec.GetSize();i++)
      result+=vec[i]*vec[i];
      //result=sqrt(result);
      real_result=result.real();
      // real_result=sqrt(real_result);
       std::cout<<" \n real_result"<<real_result<<std::endl;
       return real_result;
      }

  void piezoParamIdent::updateRHS(Vector<Complex> & solElecPot, Vector<Complex> & solMechDispl, Double omega){
    ENTER_FCN (" piezoParamIdent::updateRHS");
    Vector<Double> new_RHS;
    new_RHS.Resize(numNodes * dofs);
    for(int i=0;i<numNodes;i++){
      new_RHS[i]=solElecPot[i].real();
      for(int j=0;j<dofs-1;j++){	
	new_RHS[(j+1)*numNodes+i]=solMechDispl[j*numNodes+i].real();
	//	std::cout<<new_RHS[(j+1)*numNodes+i]<<"; ";
	// matVecRHS(IncrementedRHSMatrix,solMechDispl, solElecPot,newRHS);	
      }
    }

    ptAlgsys->UpdateRHS(1,new_RHS.GetPointer());
    std::cout<<"RHS was updated"<< std::endl;
  } // end update RHS

  void piezoParamIdent::updateRHS(Vector<Complex> & RHSsol){
    Vector<Double> temp;
    temp.Resize(RHSsol.GetSize());
    for(int i=0;i<RHSsol.GetSize();i++){
      temp[i]=RHSsol[i].real();
      std::cout<<RHSsol[i]<<"; ";
    }
    ENTER_FCN("piezoParamIdent::updateRHS2");
    ptAlgsys->InitRHS();
    ptAlgsys->UpdateRHS(1,temp.GetPointer());
    std::cout<<"RHS was updated"<< std::endl;

  }


  void piezoParamIdent::measureMechDeformationInZ_Direction(Vector<Complex> & mechDisplacement, Double & Radius, Double & meanValueMechDeformation, int dof){
    meanValueMechDeformation=0.0;
    for (int i=dof-2;i<mechDisplacement.GetSize();i+=(dof-1)) // only in Z-Direction
      meanValueMechDeformation+=std::abs(mechDisplacement[i]);
    // std::cout<<"meanValueMechDef: "<< meanValueMechDeformation <<std::endl;
    meanValueMechDeformation=meanValueMechDeformation/(PI*Radius*Radius*mechDisplacement.GetSize()*(dof-1));
    //    std::cout <<"\nMeanValue of MechDeformation: " << meanValueMechDeformation<< std::endl;
    // std::cout<<"AREA "<< PI*Radius*Radius*mechDisplacement.GetSize()*(dof-1) <<std::endl;
  }


  void piezoParamIdent::typeOutSolutionOnScreen(Vector<Complex> & solElecPot,Vector<Complex> & solMechDispl){
    Double sol_real, sol_imag;
    for(int i=0;i<solElecPot.GetSize();i++){
      sol_real=solElecPot[i].real();
      sol_imag=solElecPot[i].imag();
      std::cout << "solElecPot( " << i<< ")=" << sol_real << " + " << sol_imag <<" i " <<std::endl;
    }
    for(int i=0;i<solMechDispl.GetSize();i++){
      sol_real=solMechDispl[i].real();
      sol_imag=solMechDispl[i].imag();
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
    while(allMeasuredData.getline(mDataRow, 265)){
      if (mDataRow[0]=='1')
	{ //std::cout<<"1";
	  i=2;
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


  void piezoParamIdent::updateMaterialData(Vector<Double> & parameter, MaterialData * ptMaterial){
    ENTER_FCN("piezoParamIdent::updateMaterialData");    //std::cout<<"updateMaterialData"<<std::endl;

    ptMaterial->SetPiezoMatrixData(0,0, parameter[0]);
    ptMaterial->SetPiezoMatrixData(2,2, parameter[1]);
    ptMaterial->SetPiezoMatrixData(0,1, parameter[2]);
    ptMaterial->SetPiezoMatrixData(1,0, parameter[2]);
    ptMaterial->SetPiezoMatrixData(0,2, parameter[3]);
    ptMaterial->SetPiezoMatrixData(2,0, parameter[3]);
    ptMaterial->SetPiezoMatrixData(1,2, parameter[3]);
    ptMaterial->SetPiezoMatrixData(1,2, parameter[3]);
    ptMaterial->SetPiezoMatrixData(3,3, parameter[4]);
    ptMaterial->SetPiezoMatrixData(4,4, parameter[4]);
    ptMaterial->SetPiezoMatrixData(5,5, 0.5*(parameter[0]-parameter[2]));
    ptMaterial->SetPiezoMatrixData(6,4, parameter[5]);
    ptMaterial->SetPiezoMatrixData(7,3, parameter[5]);
    ptMaterial->SetPiezoMatrixData(8,0, parameter[6]);
    ptMaterial->SetPiezoMatrixData(8,1, parameter[6]);
    ptMaterial->SetPiezoMatrixData(8,2, parameter[7]);
    ptMaterial->SetPiezoMatrixData(6,6, parameter[8]);
    ptMaterial->SetPiezoMatrixData(7,7, parameter[8]);
    ptMaterial->SetPiezoMatrixData(8,8, parameter[9]);

  } // end updateMaterialData


} // end namespace CoupledField

/*	for(int i=0;i<nrMeasuredData;i++){
	amplitude_phase[i]=Complex(real[i],imag[i]);
	std::cout <<"Phase&Amplitude: " <<  amplitude_phase[i].real() <<" + "<< amplitude_phase[i].imag() << "i" <<std::endl;
	amplitude_phase[i]=std::polar(real[i],imag[i]);
	std::cout <<"Phase&Amplitude (polar): " <<  amplitude_phase[i].real() <<" + "<< amplitude_phase[i].imag() << "i"    <<std::endl; }





	std::cout << " SizeOf amplitude_phase" <<amplitude_phase.GetSize()<< std::endl;
	std::cout << " SizeOf real" <<real.GetSize()<< "nrmeasuredData= " << nrMeasuredData << std::endl;
*/
