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
#include "Forms/forms_header.hh"


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

#include <stdlib.h>
#include <sstream>
#include <iomanip>



#include "Utils/tools.hh"
#include <PDE/pdes_header.hh>

#include <Driver/piezoParamIdent.hh>



//#include "/../OLAS/algsys/basesystem.hh"
//#include "DataInOut/piezoParameterData.hh"



namespace CoupledField
{

  // ========================================================================
  // ========================= piezoParamIdent - Part ===========================
  // ========================================================================

//   void piezoParamIdent::calc_measuredCharge(Vector<Double> freqs, Vector<Double> & absZ, Vector<Double> & phi, Vector<Complex> & y_hat){
//     ENTER_FCN("piezoParamIdent::calc_measuredCharge");
//     Complex Z,j;
//     Double x, y;
//     j=Complex(0,1);
//     for (int i=0; i<nrMeasuredData; i++){
//       x=absZ[i]*cos(phi[i]);
//       y=absZ[i]*sin(phi[i]);
//       Z=Complex(x,y);
//       y_hat[i]=voltage/(2.0*PI*Z*freqs[i]*j);
//       //      y_hat[i]=voltage/(Z*freqs[i]*j);
//       std::cout<<i<<") = " << phi[i] << ",\t "<< freqs[i] <<",\t q = y_hat = " << y_hat[i]<<",\t Z= " << Z << std::endl;
//     }
//   }// end calc_measuredCharge()


  void piezoParamIdent::calc_measuredCharge(Vector<Double> freqs, Vector<Double> & absZ, Vector<Double> & phi, Vector<Complex> & y_hat){
    ENTER_FCN("piezoParamIdent::calc_measuredCharge");
    Complex Z,j;
    Double x, y;
    j=Complex(0,1);
    Double phase;
    for (int i=0; i<nrMeasuredData; i++){
      x=absZ[i]*cos(PI/180*phi[i]);
      y=absZ[i]*sin(PI/180*phi[i]);
      Z=Complex(x,y);
      phase = 180.0/PI*std::arg(Z);
    
	y_hat[i]=sign*voltage/(2.0*PI*Z*freqs[i]*j);
      //      y_hat[i]=voltage/(Z*freqs[i]*j);

	std::cout<<"\n Frequenz; " << freqs[i] << ", messZ: " << absZ[i] << ", phase: " << phi[i] << std::endl;
	std::cout<<" Frequenz; " << freqs[i] << ", calcZ: " << std::abs(Z) << ", phase: " << phase << std::endl << std::endl;
     
      //      std::cout<<i<<") = " << phi[i] << ",\t "<< freqs[i] <<",\t q = y_hat = " << y_hat[i]<<",\t Z= " << Z << " phase " << phase << std::endl;
      }
    

  }// end calc_measuredCharge()


  void piezoParamIdent::calcSyntheticData(Vector<Complex> & y_hat){
    ENTER_FCN("piezoParamIdent::calcSyntheticData");
    //   std::cout<<"\n We are generating synthetic data, i.e. we solve the piezo-equation with exact material - parameters"<<std::endl;
    //    std::cout<<"and alienate the results by alternating +-10Percent"<<std::endl;

    MaterialData * ptMaterial=pdes_[0]->getPDEMaterialData();   // Pointer to MaterialData
    Matrix<Double> * matMatrix =  ptMaterial->GetMatrix();
    //    std::cout<<*matMatrix<<std::endl;
    ptBCs = pdes_[0]->getPDE_BCs();     
    
    createF(ptMaterial, ptBCs, y_hat); // calculates only forward problems over all omegas

    for (Integer i=0;i<y_hat.GetSize();i++){
      if (i%2==0)
	y_hat[i]*=1.01;
      else
	y_hat[i]*=0.99;
    }
  }// end calcSyntheticData


  void piezoParamIdent::calcImpedanceCurve(){
    ENTER_FCN("PiezoParamIdent::caclImpedanceCurve");

    ptMaterial=pdes_[0]->getPDEMaterialData();   // Pointer to MaterialData
    ptBCs = pdes_[0]->getPDE_BCs();                             // Pointer to BCs
    Integer level=0;
    Boolean reset = TRUE;
    Integer pdenumber = 0;

    reset = TRUE;
    for (Integer fstep = 0; fstep < freqs.GetSize(); fstep++) { // harmonic solver for different frequency - values
    
      ////////////////////////////////////////////////////////
      //                   SOLVES PDE                      //
      ///////////////////////////////////////////////////////  
      pdes_[0]->PreStepHarmonic(fstep, freqs[fstep], level, reset); 
      pdes_[0]->SolveStepHarmonic(fstep, freqs[fstep], level, reset);
      pdes_[0]->PostStepHarmonic(fstep, freqs[fstep], level, reset);
      pdes_[0]->PostProcess(level);
	
      /////////////////////////////////////////////////////////

	Vector<Complex> chargeVec = pdes_[0]->getPDE_complexValuedCharge(); // Vector wich contains charges for each element !

	Complex charge=Complex(0.0,0.0);

	for (int i=0;i<chargeVec.GetSize();i++){
	  charge+=chargeVec[i];
	}
	//	calcAbsImped(charge, freqs[fstep], fstep);   // calculates |Z| and writes results in File
	//	charge=-charge/Complex(chargeVec.GetSize());

	Double imped, phase;
	Complex impedC;

	if (!impedCurve)
	  std::cerr<<"Error opening 'ImpedCurve.dat' "<<std::endl;

	Complex im=Complex(0.0,1);
	impedC=voltage/(charge*2.0*PI*freqs[fstep]*im);
	imped = std::abs(voltage/(charge*2.0*PI*freqs[fstep]*im)); 
	//	phase = 180.0/PI*(std::arg(charge));
	phase = 180.0/PI*(std::arg(impedC));
	std::cout << fstep <<");\t Frequenz: " << freqs[fstep] << ";\t Impedanz: "<< imped << ";\t Phase: " << phase << std::endl;
	*impedCurve <<"\n" << freqs[fstep] << " " << imped << "  " << phase << "  " << impedC.real()<<"  " << impedC.imag() << std::endl;

    } //  end loop over freqs

  } // end calcImpedance Curve

  void piezoParamIdent::createF(MaterialData * ptMaterial, BCs * ptBCs, Vector<Complex> & F_hat){
    ENTER_FCN("PiezoParamIdent:createF");
    //  std::cout<<"\nF wil be created ..."<<std::endl;

    ptMaterial=pdes_[0]->getPDEMaterialData();   // Pointer to MaterialData
    ptBCs = pdes_[0]->getPDE_BCs();                             // Pointer to BCs
    Grid * ptGrid =   ptdomain_->GetGrid();

    Integer level=0;
    Boolean reset = TRUE;
    Integer pdenumber = 0;

    Integer numElems = pdes_[0]->getPDE_numElems();
    Integer dofs=pdes_[0]->getPDE_dofspernode();  
    Integer numNodes= pdes_[0]->getPDE_numPDENodes();  

      
    // Explicit calculation of tensor matrices

    //  Matrix<Double> couplingMatrix; // \bf e
    //  Matrix<Double> dielectricMatrix; // \eps^S
    //  couplingMatrix.Resize(spaceDim, 2*spaceDim);
    //  dielectricMatrix.Resize(spaceDim, spaceDim);
    //  createMaterialTensorMatrices(parameter, couplingMatrix, dielectricMatrix, spaceDim); 

   
    Matrix<Double> * matMatrix =  ptMaterial->GetMatrix();
    //    std::cout<<*matMatrix<<std::endl;

       for (Integer fstep = 0; fstep < nrMeasuredData; fstep++) { // harmonic solver for different frequency - values

      //      Info->WriteHarmonicStep(pdes_[0]->GetName(), fstep, freqs[fstep]);
    

      ////////////////////////////////////////////////////////
      //                   SOLVES PDE                      //
      ///////////////////////////////////////////////////////  

	//	pdes_[0]->WriteGeneralPDEdefines();   // should not be used, overwrites to much!!    

	pdes_[0]->PreStepHarmonic(fstep, freqs[fstep], level, reset); 
	
	//	updateMaterialData(parameter,ptMaterial);
   
	//std::cout<<"\n piezoParam:createF SolveStepHarmonic"<<std::endl;
	pdes_[0]->SolveStepHarmonic(fstep, freqs[fstep], level, reset);

	//std::cout<<"\n piezoParam:createF setBCs"<<std::endl;
	//	pdes_[0]-> setBCs_id_phase_(0, imag[fstep]);      

	//std::cout<<"\n piezoParam:createF PostStepHarmonic"<<std::endl;
	pdes_[0]->PostStepHarmonic(fstep, freqs[fstep], level, reset);

	//	std::cout<<"\n piezoParam:createF PostProcess at step  "<< fstep << std::endl;
	pdes_[0]->PostProcess(level);


	//std::cout<<"\n piezoParam:createF PostStepHarmonic"<<std::endl;

        //////////////////////////////////////////////////////////
	//Retrieves & stores Solution for further calculations  //
	/////////////////////////////////////////////////////////

	  BaseNodeStoreSol * ptSol = pdes_[0]->getPDESolution();
	  NodeStoreSol<Complex> * ptNodeStoreSol;
	  ptNodeStoreSol = dynamic_cast<NodeStoreSol<Complex>*>(ptSol);	  

	  if (considerMechDeformation==TRUE){
	    ptNodeStoreSol->GetGlobalSolVector(ELEC_POTENTIAL, solElecPot);
	    ptNodeStoreSol->GetGlobalSolVector(MECH_DISPLACEMENT, solMechDispl);
	    //typeOutSolutionOnScreen(solElecPot, solMechDispl);              //member function of piezoParamIdent
	    meanValueMechDeformation=0.0;
	    measureMechDeformationInZ_Direction(solMechDispl,radius,meanValueMechDeformation,dofs); 
	    // std::cout<<"meanValueMechDef: "<< meanValueMechDeformation <<std::endl; 
	    F_hat[fstep + nrMeasuredData]=meanValueMechDeformation;
	  }

	  Vector<Complex> chargeVec =   pdes_[0]->getPDE_complexValuedCharge(); // Vector wich contains charges for each element !

	  Complex charge=Complex(0.0,0.0);

	   //      std::cout<<"\n Mean - Value Charge: "<< mean_value_charge << " for frequency " << freqs[fstep]<< std::endl;
	 
	  for (int i=0;i<chargeVec.GetSize();i++){
	      //	std::cout<<"\n|charge("<<i<<")|= "<<std::abs(chargeVec[i])<<";\t charge("<<i<<")= "<<(chargeVec[i]);
	    charge+=chargeVec[i];
	  }
	  F_hat[fstep]=sign*charge; 
	 
	  	  calcAbsImped(charge, freqs[fstep], fstep);   // calculates |Z| and writes results in File
     
	  //      Info->PrintPiezoMat(*ptMaterial);
	  pdes_[0]->WriteResultsInFile();     
      
	  //Gathers solution (d_k^l, \phi_k^l) to provide it for RHS of PDE which solution gives F'
	  //      Vector<Complex> algSysSolVector;
	  // algSysSolVector=ptNodeStoreSol->GetAlgSysVector();
	  //if (fstep==0)
	  //	completeSolOf_F.Resize(nrMeasuredData,algSysSolVector.GetSize());
	  //for (int i=0;i<algSysSolVector.GetSize();i++)
	  //	completeSolOf_F[fstep][i]=algSysSolVector[i];

	  // BIG CHANGE!!!!
	  // Instead of completeSolof_f, we try to save all element results ...
      
	  StdVector<Elem*> elemssd;
	  subdoms = pdes_[0]->getPDE_subdoms();
	  //	std::cout<<"PDE_SUBDOM[0] = "<< subdoms[0]; 
	  ptGrid->GetElemSD(elemssd,subdoms[0], level);
     
	  if (fstep==0)
	    allElemsVec.Resize(nrMeasuredData,elemssd.GetSize()*dofs*numNodes);
	  for (int actEl=0;actEl<elemssd.GetSize();actEl++){
	    BaseFE * ptEl = elemssd[actEl]->ptElem;
	    StdVector<Integer> connecth = elemssd[actEl]->connect;

	    Vector<Complex> elSolVec; 
	    ptNodeStoreSol->GetElemSolution(elSolVec,connecth);
	    Matrix<Double> coordinateMatrix;

	    //	pdes_[0]->GetElemCoords(connecth, coordinateMatrix,0);
	    //	std::cout<<"\n coordinateMatrix:" <<std::endl;
	    //	std::cout<<coordinateMatrix<<std::cout;

	    for (int k=0;k<elSolVec.GetSize();k++)
	      allElemsVec[fstep][actEl*elSolVec.GetSize()+k] = elSolVec[k];
	  } // end loop over elements

    } // end of loop over all frequencies

        std::cout<<"\nFinished to create F ... here it is:"<<std::endl;
           for (int i=0;i<F_hat.GetSize();i++)
             std::cout<<"F("<<i<<")="<<F_hat[i]<<"; \t";
    std::cout<<"\n ------------------------------- " <<std::endl;

    // std::cout<<" \n all ElemsVec: "<<std::endl;
//     Complex nullC = Complex(0.0,0.0);
//  for(int i=0;i<allElemsVec.GetSizeRow();i++)
//       for (int j=0; j<allElemsVec.GetSizeCol();j++){
// 	//	std::cout<<JacobiMatrix[i][j].real()<<"+"<<JacobiMatrix[i][j].imag()<<"i"<< ", ";
// 	//	if (allElemsVec[i][j]!=nullC)
// 	std::cout<<"aEv'("<<i<<")("<<j<<")= "<< allElemsVec[i][j]<<"; \t";
//     	if (j==allElemsVec.GetSizeCol()-1)
//     	  std::cout<<";\n";
//       }



  } // end createF



  void piezoParamIdent::createJacobiMatrix(MaterialData * ptMaterial, BCs * ptBCs, Vector<Complex> & F_hat, Vector<Double> & parameterIncrement, Matrix<Complex> & JacobiMatrix, Vector<Complex> & solElecPot,Vector<Complex> & solMechDispl){
    ENTER_FCN("piezoParamIdent::createJacobiMatrix");
    std::cout<<"JacobiMatrix will be created"<<std::endl;
    Vector<Double> IncrementedRHSMatrix;   
  
    //    Matrix<Double> * matMatrix =  ptMaterial->GetMatrix();
    ptBCs = pdes_[0]->getPDE_BCs();                             // Pointer to BCs
    ptAlgsys = pdes_[0]->getPDE_algsys();
    Integer level=0;
    Boolean reset = TRUE;
    ptGrid = pdes_[0]->getPDE_grid();
    ptNodeEqn = pdes_[0]->getPDE_eqnData();
    ptAssemble = pdes_[0]->getPDE_assemble();
    //Integer nrMeasuredData = freqs.GetSize();
    Integer numElems_ = pdes_[0]->getPDE_numElems();
    nrParameter=parameter.GetSize();

    Integer job;

    JacobiMatrix.Resize(2*nrMeasuredData,nrParameter);
    if (considerMechDeformation==FALSE)
      JacobiMatrix.Resize(nrMeasuredData,nrParameter);

    Integer pdenumber = 0;

    Vector<Double> dparameter(nrParameter);

    for(int ind_param=0; ind_param<nrParameter;ind_param++){ // loop over different parameter increments
      
      //-------------       first strategy --------------------------------------
      parameter[ind_param]+= 1.0/(scaling[ind_param]);       // we are incrementing one parameter after another with 1.0/scale
      if (ind_param>0)
      	parameter[ind_param-1]-=1.0/(scaling[ind_param-1]);
      updateMaterialData(parameter, ptMaterial);         //member function of piezoParamIdent, recalculates stiffness, etc. ...

      if (ind_param==nrParameter-1)
 	parameter[ind_param]-=1.0/scaling[ind_param];

      //std::cout<<"\n"<<std::endl;
      //for(Integer i=0;i<parameter.GetSize();i++)
      //   std::cout<<parameter[i]<<"; ";
      
      //------------------------------------------------------------

      // ~~~~~~~~~~~~~~~~~~~~~  second strategy ~~~~~~~~~~~~~~~~~~~~~~~~~
      //        dparameter[ind_param]=1.0/scaling[ind_param];
      //        if (ind_param>0)
      //  	dparameter[ind_param-1]=0.0; 
      //        updateMaterialData(dparameter, ptMaterial);    
      //        if (ind_param==nrParameter-1)
      //  	dparameter[ind_param-1]=0.0; 

      //       std::cout<<"\n"<<std::endl;
      //       for(Integer i=0;i<parameter.GetSize();i++)
      //         	std::cout<<dparameter[i]<<"; ";

      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



      //      pdes_[0]->DefineIntegrators(0);

      for (Integer fstep = 0; fstep < nrMeasuredData; fstep++) { // harmonic solver for different frequency - values
	reset = TRUE;

        pdes_[0]-> setBCs_id_phase_(0, imag[fstep]);

	//	Info->WriteHarmonicStep(pdes_[0]->GetName(), fstep, freqs[fstep]);
	
	// 	pdes_[0]->WriteGeneralPDEdefines();
 	pdes_[0]->PreStepHarmonic(fstep, freqs[fstep], level, reset);
        
	//	Cannot use SolveStepHarmonic, since it overwrites RHS ...
	//      for this, I have copied the method StepHarmonicLin to this place 
	//      void BasePDE::StepHarmonicLin(const Integer freqStep, const Double frequency, const Integer level, const Boolean reset)
	
	ptAssemble->AssembleMatrices(level);
       
	// The folowing method creates and calculates the RHS for harmonic Problems in CreateJacobian ...
	createAndSetRHSforJacobian(fstep);

	//this has to be done each time!
	ptAssemble->AssembleSrcRHS(level, freqs[fstep]);

        Double * ptsol;
 	BaseNodeStoreSol * ptSol = pdes_[0]->getPDESolution();
 	NodeStoreSol<Complex> * ptNodeStoreSol;
 	ptNodeStoreSol = dynamic_cast<NodeStoreSol<Complex>*>(ptSol);
		
	pdes_[0]-> setBCs_id_phase_(0, imag[fstep]);

	if (reset)
	  {
	    //account for bcs
	    pdes_[0]->SetBCs(level, freqs[fstep]);
	    job = 1; // calc new preconditioner
	  }
	else
	  job = 3;
	//	  std::cout<<"SetBcs ..."<<std::endl;
#ifdef USE_OLAS
	ptAlgsys->BuildInDirichlet();
	ptAlgsys->SetupPrecond(job);
#else
	ptAlgsys->CalcPrecond(job);
#endif
	ptAlgsys->Solve();
	ptsol = ptAlgsys->GetSolutionVal();
	ptSol->CopyFromAlgSysDataPointer(ptsol);



	if (considerMechDeformation==TRUE){	  
	  //	ptNodeStoreSol->GetGlobalSolVector(ELEC_POTENTIAL, solElecPot);
	  ptNodeStoreSol->GetGlobalSolVector(MECH_DISPLACEMENT, solMechDispl);
	  //	  typeOutSolutionOnScreen(solElecPot, solMechDispl);              //member function of piezoParamIdent
	  measureMechDeformationInZ_Direction(solMechDispl,radius,meanValueMechDeformation,dofs); // Braucht üÜberarbeitung!!
	  JacobiMatrix[fstep + nrMeasuredData][ind_param]=meanValueMechDeformation;    
	}

        pdes_[0]->PostStepHarmonic(fstep, freqs[fstep], level, reset);
	pdes_[0]->PostProcess(level);
        pdes_[0]->PostStepHarmonic(fstep, freqs[fstep], level, reset);
	Vector<Complex> chargeVec =  pdes_[0]->getPDE_complexValuedCharge();

	Complex charge=Complex(0.0,0.0);


	for (int i=0; i<chargeVec.GetSize();i++){
	  charge=charge+chargeVec[i];

	}
       
	//	mean_value_charge = mean_value_charge/(Double(chargeVec.GetSize()));
	//	std::cout<<"\nCHARGE VEC SIZE = "<< chargeVec.GetSize()<<"mean-value-charge = " << mean_value_charge << std::endl;

	  JacobiMatrix[fstep][ind_param]=sign*charge; //pdes_[0]->getPDE_complexValuedCharge();    


	//   Info->PrintPiezoMat(*ptMaterial);
	//     pdes_[0]->WriteResultsInFile();     
     
      }		// end for over all frequencies
      std::cout << " \nCreateJacobian Line " << ind_param  <<std::endl;
      //for (int i=0;i<JacobiMatrix.GetSizeRow();i++)
      //	std::cout<<"F'("<<i<<")("<<ind_param<<")= "<< JacobiMatrix[i][ind_param]<<"; ";


    } //end loop over paramters
    
    std::cout<<"JACOBI - MATRIX 1: " <<std::endl;
    for(int i=0;i<JacobiMatrix.GetSizeRow();i++)
      for (int j=0; j<JacobiMatrix.GetSizeCol();j++){
	//	std::cout<<JacobiMatrix[i][j].real()<<"+"<<JacobiMatrix[i][j].imag()<<"i"<< ", ";
	std::cout<<"F'("<<i<<")("<<j<<")= "<< JacobiMatrix[i][j]<<"; \t";
    	if (j==JacobiMatrix.GetSizeCol()-1)
    	  std::cout<<";\n";
      }

  }            //end CreateJacobiMatrix


  void piezoParamIdent::createJacobiMatrix2(Matrix<Complex> & JacobiMatrix){
    ENTER_FCN("piezoParamIdent::createJacobiMatrix2");
    std::cout<<"JacobiMatrix2 will be created"<<std::endl;

  
    //    Matrix<Double> * matMatrix =  ptMaterial->GetMatrix();
    ptBCs = pdes_[0]->getPDE_BCs();                             // Pointer to BCs
    ptAlgsys = pdes_[0]->getPDE_algsys();
    Integer level=0;
    Boolean reset = TRUE;
    ptGrid = pdes_[0]->getPDE_grid();
    ptNodeEqn = pdes_[0]->getPDE_eqnData();
    ptAssemble = pdes_[0]->getPDE_assemble();
    //Integer nrMeasuredData = freqs.GetSize();
    Integer numElems_ = pdes_[0]->getPDE_numElems();
    nrParameter=parameter.GetSize();
	Integer job;
    
    Integer spaceDim = pdes_[0]->getPDE_spaceDim();
    Double * ptsol;
    StdVector<Elem*> elemssd;
    subdoms = pdes_[0]->getPDE_subdoms();
    ptMaterial=pdes_[0]->getPDEMaterialData();   // Pointer to MaterialData
    ptGrid->GetElemSD(elemssd,subdoms[0], level); // gets element list elemssd

    BaseNodeStoreSol * ptSol = pdes_[0]->getPDESolution();
    NodeStoreSol<Complex> * ptNodeStoreSol;
    ptNodeStoreSol = dynamic_cast<NodeStoreSol<Complex>*>(ptSol);
    Vector<Complex> algSysSolVector;
    algSysSolVector=ptNodeStoreSol->GetAlgSysVector();
    Vector<Complex> RHSVec(algSysSolVector.GetSize());

    JacobiMatrix.Resize(2*nrMeasuredData,nrParameter);
    if (considerMechDeformation==FALSE)
      JacobiMatrix.Resize(nrMeasuredData,nrParameter);

    Integer pdenumber = 0;

    Vector<Double> dparameter(nrParameter);

    Vector<Double> tempHarm;
    StdVector<Integer> connect_PDE;	   	    

    Matrix<Double> *matMat = ptMaterial->GetMatrix();

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

    for(int ind_param=0; ind_param<nrParameter;ind_param++){ // loop over different parameter increments
      

      // ~~~~~~~~~~~~~~~~~~~~~  second strategy ~~~~~~~~~~~~~~~~~~~~~~~~~
      //              dparameter[ind_param]=1.2/scaling[ind_param];
      // diese Zeile funzt!!      
      //     dparameter[ind_param]=1.3/scaling[ind_param]*bas[ind_param]; // 1.0/scaling[ind_param];

      //andere Versuche:
      dparameter[ind_param]=1.0/scaling[ind_param];//*bas[ind_param]; // 1.0/scaling[ind_param];
      
     if (ind_param>0)
      	dparameter[ind_param-1]=0.0; 

      updateMaterialData(dparameter, ptMaterial);    

      //       std::cout<<"\n"<<std::endl;
      // for(Integer i=0;i<parameter.GetSize();i++)
      //	std::cout<<dparameter[i]<<"; ";
      //  std::cout<<"\n"<<std::endl;


      if (ind_param==nrParameter-1)
	dparameter[ind_param-1]=0.0; 

      
      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      //       pdes_[0]->DefineIntegrators(0);

      for (Integer fstep = 0; fstep < nrMeasuredData; fstep++) { // harmonic solver for different frequency - values
	reset = TRUE;

	ptAlgsys->InitRHS();

      //   pdes_[0]-> setBCs_id_phase_(0, imag[fstep]);
  	
	//loop over elements
	for (int actEl=0; actEl< elemssd.GetSize(); actEl++) {
	  BaseFE * ptEl = elemssd[actEl]->ptElem;
	  StdVector<Integer> connecth = elemssd[actEl]->connect;
		             
	  Matrix<Double> ptCoord;
	  ptGrid->GetCoordNodesElemMat(connecth, ptCoord, 0);
		    
	  // map connect to PDE node numbers

	  ptNodeEqn->Node2EQN(connecth, connect_PDE);

	  //Vector<Complex> elSolVec; 
	  // ptNodeStoreSol->GetElemSolution(elSolVec,connecth);

	  MaterialData actSDMat(*ptMaterial);
	  Boolean isdamping=TRUE;
	     
	  // Create new Integrator, with this calculate elemMat for sysmat in RHS 
	  BaseForm * bilinearStiff;
	  //	   std::cout<<"create bilinearStiff for element "<< actEl << " and for fstep = " <<fstep <<std::endl;

	  if (spaceDim==3)
	    bilinearStiff =  new linPiezo3DInt(actSDMat,isdamping);
	  else if (spaceDim==2)
	    bilinearStiff = new piezoAxiInt(actSDMat, isdamping);

	  updateMaterialData(dparameter,ptMaterial);

	  IntegratorDescriptor *actIntDescrStiff = new IntegratorDescriptor(bilinearStiff, STIFFNESS);

	  bilinearStiff->SetElemPtr(ptEl);
	  Matrix<Complex> elemMat;
	  //      Double damp_beta =1.0e-9; // in future, beta will be dependent of omega_l ...
	  Double damp_beta = ptMaterial->GetDampingBeta();

	  Double omega = 2.0*PI*freqs[fstep];
	  bilinearStiff->CalcComplexElementMatrix(ptCoord,elemMat,damp_beta,omega);

	  Complex sum=Complex(0.0,0.0);
	  Vector<Complex> temp;
	  temp.Resize(elemMat.GetSizeRow());
	  for (int i=0; i<elemMat.GetSizeRow();i++){
	    for (int j=0; j<elemMat.GetSizeCol();j++)
	      sum=sum+elemMat[i][j]*allElemsVec[fstep][fstep*elemMat.GetSizeRow()+j];  
	    temp[i]=sum;
	    sum=Complex(0.0,0.0);
	  }
	  // save values like it was done in transformElemMat2harmonic (see assemle.cc) ...
	  tempHarm.Resize(2*temp.GetSize());
	  for (int i=0; i<temp.GetSize();i++){
	    tempHarm[i]=-temp[i].real();
	    tempHarm[i+temp.GetSize()]=-temp[i].imag();
	  } 

	  ptAlgsys->SetElementRHS(&tempHarm[0], connect_PDE.GetPointer(), 
	  	  connect_PDE.GetSize());


	} // end for over all Elems

	updateMaterialData(parameter,ptMaterial);

	  // the following lines are the content of: PreStepHarmonic:
	  // pdes_[0]->PreStepHarmonic(fstep, freqs[fstep], level, reset);

	pdes_[0]->setPDE_actFrequency(freqs[fstep]);
	pdes_[0]->setPDE_actFreqStep(fstep);
	ptAssemble->SetFrequency(freqs[fstep]);


	  if (reset)
	   ptAssemble->InitMatrices();
        
	//	Cannot use SolveStepHarmonic, since it overwrites RHS ...
	//      for this, I have copied the method StepHarmonicLin to this place 
	//      void BasePDE::StepHarmonicLin(const Integer freqStep, const Double frequency, const Integer level, const Boolean reset)
	
	ptAssemble->AssembleMatrices(level);

	       
	// The folowing method creates and calculates the RHS for harmonic Problems in CreateJacobian ...
	 	
	ptAssemble->AssembleSrcRHS(level, freqs[fstep]);

 	ptNodeStoreSol = dynamic_cast<NodeStoreSol<Complex>*>(ptSol);
		
	if (reset)
	  {
	    //account for bcs
	     pdes_[0]->SetBCs(level, freqs[fstep]);
	    job = 1; // calc new preconditioner
	  }
	else
	  job = 3;

#ifdef USE_OLAS
	ptAlgsys->BuildInDirichlet();
	ptAlgsys->SetupPrecond(job);
#else
	ptAlgsys->CalcPrecond(job);
#endif
	ptAlgsys->Solve();

	ptsol = ptAlgsys->GetSolutionVal();
	ptSol->CopyFromAlgSysDataPointer(ptsol);

	if (considerMechDeformation==TRUE){
	  //      ptNodeStoreSol->GetGlobalSolVector(ELEC_POTENTIAL, solElecPot);
	  ptNodeStoreSol->GetGlobalSolVector(MECH_DISPLACEMENT, solMechDispl);
	  //	  typeOutSolutionOnScreen(solElecPot, solMechDispl);              //member function of piezoParamIdent
	  measureMechDeformationInZ_Direction(solMechDispl,radius,meanValueMechDeformation,dofs); // Braucht üÜberarbeitung!!
	  JacobiMatrix[fstep + nrMeasuredData][ind_param]=meanValueMechDeformation;    
	}

        pdes_[0]->PostStepHarmonic(fstep, freqs[fstep], level, reset);
	pdes_[0]->PostProcess(level);
        pdes_[0]->PostStepHarmonic(fstep, freqs[fstep], level, reset);

	Vector<Complex> chargeVec =  pdes_[0]->getPDE_complexValuedCharge();

	Complex charge=Complex(0.0,0.0);

	for (int i=0; i<chargeVec.GetSize();i++)
	  charge=charge+chargeVec[i];
    
	JacobiMatrix[fstep][ind_param]=sign*charge; //pdes_[0]->getPDE_complexValuedCharge();    
     
      }		// end for over all frequencies

    } //end loop over paramters

//     for(int i=0;i<JacobiMatrix.GetSizeRow();i++)
//       for (int j=0; j<JacobiMatrix.GetSizeCol();j++){
// 	//	std::cout<<JacobiMatrix[i][j].real()<<"+"<<JacobiMatrix[i][j].imag()<<"i"<< ", ";
// 		std::cout<<"F'("<<i<<")("<<j<<")= "<< JacobiMatrix[i][j]<<"; \t";
//     	if (j==JacobiMatrix.GetSizeCol()-1)
// 	 std::cout<<";\n";
//       }

  }            //end CreateJacobiMatrix 2


  void piezoParamIdent::testJacobiMatrix(Vector<Complex> & F_hat, Matrix<Complex> & JacobiMatrix, Vector<Double> & parameter,BCs * ptBCs,MaterialData * ptMaterial,Vector<Double> & parameterIncrement, Vector<Complex>& solElecPot,Vector<Complex> &solMechDispl){
    ENTER_FCN("piezoParamIdent::testJacobiMatrix");

    Vector<Complex> F_hat_incr(F_hat.GetSize());
    Matrix<Complex> approxJacobiMatrix(JacobiMatrix.GetSizeRow(), JacobiMatrix.GetSizeCol());
    Vector<Double> parameter_incr(parameter.GetSize());

    for(Integer i=0;i<parameter.GetSize();i++){
      parameter_incr[i]=1.1*parameter[i];
      //      std::cout<<"Denominator("<<i<<") = " << (Complex(parameter[i])-Complex(parameter_incr[i]))<<std::endl;
    }

    updateMaterialData(parameter_incr, ptMaterial);
    createF(ptMaterial, ptBCs, F_hat_incr);
    std::cout<<"\n"<<std::endl;
    //    for (int i=0;i<F_hat.GetSize();i++)
      //      std::cout<<"F("<<i+1<<")="<<F_hat[i]<< " <-> " << F_hat_incr[i]<<"; ";
    // createJacobiMatrix(ptMaterial, ptBCs, F_hat_incr, parameterIncrement,JacobiMatrix, solElecPot, solMechDispl);
      //    std::cout<<"\n  - - - - - - - - - - - - - - - \n JacobiMatrix <-> approxJacobiMatrix"<< std::endl;

    for (Integer i=0; i< JacobiMatrix.GetSizeRow();i++)
      for (Integer j=0; j< JacobiMatrix.GetSizeCol();j++){
	if (j<parameter.GetSize())
	  approxJacobiMatrix[i][j]=(F_hat[i]-F_hat_incr[i])/(Complex(scaling[j]*parameter[j])-Complex(scaling[j]*parameter_incr[j]));
	else {
	  Integer jj=i-parameter.GetSize();
	  approxJacobiMatrix[i][j]=(F_hat[i]-F_hat_incr[i])/(Complex(scaling[jj]*parameter[jj])-Complex(scaling[jj]*parameter_incr[jj]));
	}
	std::cout<<"F~("<<i<<")("<<j<<")= "<< approxJacobiMatrix[i][j]<<";\t "; // <<" <-> " <<approxJacobiMatrix[i][j]<< "; ";
	//	if (j==JacobiMatrix.GetSizeCol()-1)
	//std::cout<<"\n " ;
      }
    //       JacobiMatrix = approxJacobiMatrix;
  }// end testJacobiMatrix



  void piezoParamIdent::createAdjointJacobiMatrix(Matrix<Complex> & JacobiMatrix, Matrix<Complex> & adjJacobiMatrix){
    ENTER_FCN("piezoParamIdent::createAdjointJacobiMatrix");
      std::cout<<"\n Adjoint Jacobian will be created ... "<<std::endl;
    adjJacobiMatrix.Resize(JacobiMatrix.GetSizeCol(),JacobiMatrix.GetSizeRow());
    for (int i=0;i<JacobiMatrix.GetSizeCol();i++)
      for (int j=0;j<JacobiMatrix.GetSizeRow();j++){
	//adjJacobiMatrix[i][j] = JacobiMatrix[j][i];
		adjJacobiMatrix[i][j] = std::conj(JacobiMatrix[j][i]);
		//std::cout<<"F*("<<i<<")("<<j<<")= "<< adjJacobiMatrix[i][j]<<";\t ";
      }
  } // end createAdjointJacobiMatrix



  void piezoParamIdent::createAndSetRHSforJacobian(Integer & fstep)
  { 
    ENTER_FCN("piezoParamIdent::createAndSetRHSforJacobian");
    //    std::cout<<"piezoParamIdent::createAndSetRHSforJacobian 1 "<< std::endl; 
    Integer level =0;
    Integer spaceDim = pdes_[0]->getPDE_spaceDim();
    Double * ptsol;
    StdVector<Elem*> elemssd;
    subdoms = pdes_[0]->getPDE_subdoms();
    ptMaterial=pdes_[0]->getPDEMaterialData();   // Pointer to MaterialData
    ptGrid->GetElemSD(elemssd,subdoms[0], level); // gets element list elemssd

    BaseNodeStoreSol * ptSol = pdes_[0]->getPDESolution();
    NodeStoreSol<Complex> * ptNodeStoreSol;
    ptNodeStoreSol = dynamic_cast<NodeStoreSol<Complex>*>(ptSol);
    Vector<Complex> algSysSolVector;
    algSysSolVector=ptNodeStoreSol->GetAlgSysVector();
    Vector<Complex> RHSVec(algSysSolVector.GetSize());
    StdVector<Integer> connect_PDE;	   	    
	
    //loop over elements
    for (int actEl=0; actEl< elemssd.GetSize(); actEl++) {
      BaseFE * ptEl = elemssd[actEl]->ptElem;
      StdVector<Integer> connecth = elemssd[actEl]->connect;
		             
      Matrix<Double> ptCoord;
      ptGrid->GetCoordNodesElemMat(connecth, ptCoord, 0);
		    
      // map connect to PDE node numbers

      ptNodeEqn->Node2EQN(connecth, connect_PDE);

      Vector<Complex> elSolVec; 
      ptNodeStoreSol->GetElemSolution(elSolVec,connecth);

      MaterialData actSDMat(*ptMaterial);
      Boolean isdamping=TRUE;
	     
      // Create new Integrator, with this calculate elemMat for sysmat in RHS 
      BaseForm * bilinearStiff;

      if (spaceDim==3)
	bilinearStiff =  new linPiezo3DInt(actSDMat,isdamping);
      else if (spaceDim==2)
	bilinearStiff = new piezoAxiInt(actSDMat, isdamping);

      IntegratorDescriptor *actIntDescrStiff = new IntegratorDescriptor(bilinearStiff, STIFFNESS);
      bilinearStiff->SetElemPtr(ptEl);
      Matrix<Complex> elemMat;
      //      Double damp_beta =1.0e-9; // in future, beta will be dependent of omega_l ...
      Double damp_beta = ptMaterial->GetDampingBeta();
      //     std::cout<<"\n DampingBeta = " << damp_beta << std::endl;
      Double omega_temp = 2*PI*freqs[fstep];
      bilinearStiff->CalcComplexElementMatrix(ptCoord,elemMat,damp_beta,omega_temp);

      /*      std::cout <<"\n ELEMENT - MATRIX " << std::endl;
	      for (int i=0;i<elemMat.GetSizeRow();i++)
	      for (int j=0;j<elemMat.GetSizeCol();j++){
	      std::cout<<elemMat[i][j]<<"; ";
	      if (j==elemMat.GetSizeCol()-1)
	      std::cout<<"\n";
	      }*/

      //hardcoded temp = elemmat*elemvec;
      Complex sum=0.0;
      Vector<Complex> temp;
      temp.Resize(elemMat.GetSizeRow());
      for (int i=0; i<elemMat.GetSizeRow();i++){
	for (int j=0; j<elemMat.GetSizeCol();j++)
	  sum=sum+elemMat[i][j]*allElemsVec[fstep][fstep*elemMat.GetSizeRow()+j];  
	temp[i]=sum;
	sum=Complex(0,0);
      }
      // save values like it was done in transformElemMat2harmonic (see assemle.cc) ...
      Vector<Double> tempHarm(2*temp.GetSize());
      for (int i=0; i<temp.GetSize();i++){
	tempHarm[i]=temp[i].real();
	tempHarm[i+temp.GetSize()]=RHSVec[i].imag();
      } 

      ptAlgsys->SetElementRHS(&tempHarm[0], connect_PDE.GetPointer(), 
			      connect_PDE.GetSize());

    } // end for elemssd 

  } // end createAndSetRHSforJacobian();

} // end namespace

