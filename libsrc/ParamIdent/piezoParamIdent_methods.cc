#include <iostream>
#include <fstream>
#include <string>
//#include "staticdriver.hh"
#include "DataInOut/GMV/outGMV.hh"
#include "General/environment.hh"
#include <PDE/SinglePDE.hh>

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
#include "Domain/elem.hh"
#include "Domain/domain.hh"
#include "Forms/forms_header.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"


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


//#include "/../OLAS/algsys/basesystem.hh"
//#include "DataInOut/piezoParameterData.hh"



namespace CoupledField
{

  // ========================================================================
  // ========================= piezoParamIdent - Part ===========================
  // ========================================================================

  
  void piezoParamIdent::calc_measuredCharge(Vector<Double> freqs, Vector<Double> & absZ, Vector<Double> & phi, Vector<Complex> & y_hat){
    ENTER_FCN("piezoParamIdent::calc_measuredCharge");
    Complex Z,j;
    Double x, y;
    j=Complex(0,1);
    Double phase;
    Double randFactor=0.0;
    Vector<Complex> rand(nrMeasuredData);
    for (int i=0; i<nrMeasuredData; i++){
      x=absZ[i]*cos(PI/180*phi[i]);
      y=absZ[i]*sin(PI/180*phi[i]);
      Z=Complex(x,y);
      phase = 180.0/PI*std::arg(Z);
    
      y_hat[i]=sign*voltage/(2.0*PI*std::log(Z)*freqs[i]*j);
//       std::cout<<y_hat<<std::endl;
//       std::cout<<"This is y_hat!!"<<std::endl;
      //      getchar();

      //      std::cout<<"\n Frequenz; " << freqs[i] << ", messZ: " << absZ[i] << ", phase: " << phi[i] << std::endl;
      //std::cout<<" Frequenz; " << freqs[i] << ", calcZ: " << std::abs(Z) << ", phase: " << phase << std::endl << std::endl;
     
      //      std::cout<<i<<") = " << phi[i] << ",\t "<< freqs[i] <<",\t q = y_hat = " << y_hat[i]<<",\t Z= " << Z << " phase " << phase << std::endl;
    }
   
    if (TRUE){
      for (Integer i=0;i<nrMeasuredData;i++){
	rand[i] = Complex(2.0*Double(std::rand())/RAND_MAX-1);
	if (randFactor<std::abs(rand[i]))
	  randFactor = delta/std::abs(rand[i]);
    }
      for (Integer i=0;i<nrMeasuredData;i++)
	rand[i]=Complex(randFactor*rand[i])*y_hat[i];
      if(delta!=0.0)
	std::cout<<"\n Random noise with data error delta = "<< delta<<std::endl;
      //	std::cout<<rand<<std::endl;
      //	std::cout<<y_hat<<std::endl;
     
      for (Integer i=0;i<nrMeasuredData;i++)
	y_hat[i]=y_hat[i]+rand[i];
    
      Double average_error=0.0;
      for (Integer i=0;i<nrMeasuredData;i++)
	average_error+=std::abs((y_hat[i]-rand[i])/y_hat[i]);
      average_error/=nrMeasuredData;
      if(delta!=0.0){
	//std::cout<<"\n average_error = " <<average_error<<std::endl;
	std::cout<<"\n The average data error is about ~ " << std::abs(average_error-1.0)*100<<" % " << std::endl;
	std::cout<<"\n Press any key to continue ... " <<std::endl;
	getchar();
      }
      
    }

  }// end calc_measuredCharge()


  void piezoParamIdent::calcSyntheticData(Vector<Complex> & y_hat){
    ENTER_FCN("piezoParamIdent::calcSyntheticData");
    //   std::cout<<"\n We are generating synthetic data, i.e. we solve the piezo-equation with exact material - parameters"<<std::endl;
    //    std::cout<<"and alienate the results by alternating +-10Percent"<<std::endl;

    ptMaterial= ptMyPDE_->getPDEMaterialData();   // Pointer to MaterialData
    Matrix<Double> * matMatrix =  ptMaterial->GetMatrix();
    //    std::cout<<*matMatrix<<std::endl;
    ptBCs = ptMyPDE_->getPDE_BCs();     
    
    createF(ptMaterial, ptBCs, y_hat,TRUE); // calculates only forward problems over all omegas

    for (Integer i=0;i<y_hat.GetSize();i++){
      if (i%2==0)
        y_hat[i]*=1.01;
      else
        y_hat[i]*=0.99;
    }
  }// end calcSyntheticData


  void piezoParamIdent::calcImpedanceCurve(){
    ENTER_FCN("PiezoParamIdent::caclImpedanceCurve");

    ptMaterial=ptMyPDE_->getPDEMaterialData();   // Pointer to MaterialData
    ptBCs = ptMyPDE_->getPDE_BCs();                             // Pointer to BCs
    Integer level=0;
    Boolean reset = TRUE;
    Integer pdenumber = 0;
    ptAssemble = ptMyPDE_->getPDE_assemble();

    // ptAlgsys->InitMatrix();
    // ptAssemble->SetReassemble();
   
    updateMaterialData(parameter,ptMaterial);
    updateComplexMaterialData(parameterC,ptMaterial);

    //    ptMaterial->RotateMaterialMatrix(1,0,1);

    Boolean  adjustDamping = params->IsSet("adjustDamping",  "harmonic");
    if(adjustDamping)
      ptPDE_->getPDE_assemble()->SetStartFrequency(freqs[0]);

    for (Integer fstep = 0; fstep < freqs.GetSize(); fstep++) { 

//       if (reset){
//         ptAlgsys->InitMatrix();
//         ptAssemble->SetReassemble();
//         }
//       ptAssemble->CreateMatrices();
      //  reset=TRUE;

      // harmonic solver for different frequency - values
      
      //	Info->WriteHarmonicStep(ptMyPDE_->GetName(), fstep, freqs[fstep]);
      
      ////////////////////////////////////////////////////////
      //                   SOLVES PDE                      //
      ///////////////////////////////////////////////////////  
        ptMyPDE_->GetSolveStep()->PreStepHarmonic(fstep, freqs[fstep], level, reset); 

        ptMyPDE_->GetSolveStep()->SolveStepHarmonic(fstep, freqs[fstep], level, reset);

        ptMyPDE_->GetSolveStep()->PostStepHarmonic(fstep, freqs[fstep], level, reset);

        ptMyPDE_->PostProcess(level);   
        /////////////////////////////////////////////////////////

	  ptMyPDE_->WriteResultsInFile();

        
          Vector<Complex> chargeVec = ptMyPDE_->getPDE_complexValuedCharge(); // Vector wich contains charges for each element !

          Complex charge=Complex(0.0,0.0);

          for (int i=0;i<chargeVec.GetSize();i++){
            charge+=chargeVec[i];
          }
          //    calcAbsImped(charge, freqs[fstep], fstep);   // calculates |Z| and writes results in File
          //    charge=-charge/Complex(chargeVec.GetSize());

          Double imped, phase;
          Complex impedC;

          if (!impedCurve)
            std::cerr<<"Error opening 'ImpedCurve.dat' "<<std::endl;


          Complex im=Complex(0.0,1);
          impedC=voltage/(charge*2.0*PI*freqs[fstep]*im);
          imped = std::abs(voltage/(charge*2.0*PI*freqs[fstep]*im)); 
          //    phase = 180.0/PI*(std::arg(charge));
          phase = 180.0/PI*(std::arg(impedC));
          std::cout << fstep <<");\t Frequenz: " << freqs[fstep] << ";\t Impedanz: "<< std::log(imped) 
		    << ";\t Phase: " << phase <<";\t Volt = "<<voltage<<";\t Charge = "<< charge<< std::endl;
          *impedCurve <<"\n" << freqs[fstep] << " " << std::log(imped) << "  " << phase << "  " 
		      << impedC.real()<<"  " << impedC.imag() << "  " << charge.real()<< "  " << charge.imag()<< std::endl;

    } //  end loop over freqs

  } // end calcImpedance Curve

  void piezoParamIdent::createF(MaterialData * ptMaterial, BCs * ptBCs, Vector<Complex> & F_hat, Boolean typeOut){
    ENTER_FCN("PiezoParamIdent:createF");
    //   std::cout<<"\nF wil be created ..."<<std::endl;

    ptMaterial=ptMyPDE_->getPDEMaterialData();   // Pointer to MaterialData
    ptBCs = ptMyPDE_->getPDE_BCs();                             // Pointer to BCs
    Grid * ptGrid =   ptdomain_->GetGrid();
    ptAssemble = ptMyPDE_->getPDE_assemble();
    ptAlgsys = ptMyPDE_->getPDE_algsys();

    Integer level=0;
    Boolean reset = TRUE;
    Integer pdenumber = 0;

    Integer numElems = ptMyPDE_->getPDE_numElems();
    Integer dofs=ptMyPDE_->getPDE_dofspernode();  
    Integer numNodes= ptMyPDE_->getPDE_numPDENodes();  

    //    updateMaterialData(parameter,ptMaterial);
    //updateComplexMaterialData(parameterC,ptMaterial);

    //Matrix<Double> * matMatrix =  ptMaterial->GetMatrix();

    //std::cout<<*matMatrix<<std::endl;

    //    ptMyPDE_->Init();
    //    ptMyPDE->DeleteAlgSys();
    // ptAlgsys->InitMatrix();
    // ptAssemble->SetReassemble();
    //    ptAssemble->InitRHS();
        
    ptAlgsys->InitMatrix();
    ptAssemble->SetReassemble();
    ptAlgsys->InitRHS();

    ptAlgsys->InitMatrix();
    ptAssemble->SetReassemble();
    ptAlgsys->InitRHS();
      
    // updateMaterialData(parameter,ptMaterial);
//     updateComplexMaterialData(parameterC,ptMaterial);


    Boolean aTime=TRUE;
    //  ptAssemble->setPDE_readInAnotherTime(aTime);
    
    Boolean  adjustDamping = params->IsSet("adjustDamping",  "harmonic");
    if(adjustDamping)
      ptPDE_->getPDE_assemble()->SetStartFrequency(freqs[0]);
    
    for (Integer fstep = 0; fstep < nrMeasuredData; fstep++) { // harmonic solver for different frequency - values

    
      ////////////////////////////////////////////////////////
      //                   SOLVES PDE                      //
      ///////////////////////////////////////////////////////  

      ptMyPDE_->WriteGeneralPDEdefines();   // should not be used, overwrites to much!!    

      //      std::cout<<"\n piezoParam:createF PreStepHarmonic"<<std::endl;

      reset=TRUE;
      ptMyPDE_->GetSolveStep()->PreStepHarmonic(fstep, freqs[fstep], level, reset); 
         
        
      //         updateMaterialData(parameter,ptMaterial);
      //updateComplexMaterialData(parameterC,ptMaterial);

      //  std::cout<<"\n piezoParam:createF SolveStepHarmonic"<<std::endl;
      ptMyPDE_->GetSolveStep()->SolveStepHarmonic(fstep, freqs[fstep], level, reset);
      //std::cout<<"\n after SolveStepHarm " <<std::endl;

      //        std::cout<<"\n piezoParam:createF PostStepHarmonic"<<std::endl;
      ptMyPDE_->GetSolveStep()->PostStepHarmonic(fstep, freqs[fstep], level, reset);

      //std::cout<<"\n piezoParam:createF PostProcess at step  "<< fstep << std::endl;
      ptMyPDE_->PostProcess(level);


      //////////////////////////////////////////////////////////
      //Retrieves & stores Solution for further calculations  //
      /////////////////////////////////////////////////////////

        BaseNodeStoreSol * ptSol = ptMyPDE_->getPDESolution();
        NodeStoreSol<Complex> * ptNodeStoreSol;
        ptNodeStoreSol = dynamic_cast<NodeStoreSol<Complex>*>(ptSol);     

        if (considerMechDeformation==TRUE){
          ptNodeStoreSol->GetGlobalSolVector(ELEC_POTENTIAL, solElecPot);
          ptNodeStoreSol->GetGlobalSolVector(MECH_DISPLACEMENT, solMechDispl);
          typeOutSolutionOnScreen(solElecPot, solMechDispl);              //member function of piezoParamIdent
          meanValueMechDeformation=0.0;
          measureMechDeformationInZ_Direction(solMechDispl,radius,meanValueMechDeformation,dofs); 
	  //          std::cout<<"meanValueMechDef: "<< meanValueMechDeformation <<std::endl; 
          F_hat[fstep + nrMeasuredData]=meanValueMechDeformation;
        }

        Vector<Complex> chargeVec =   ptMyPDE_->getPDE_complexValuedCharge(); // Vector wich contains charges for each element !

        Complex charge=Complex(0.0,0.0);

         
        for (int i=0;i<chargeVec.GetSize();i++){
          charge+=chargeVec[i];
        }

	Double x=real[fstep]*cos(PI/180*imag[fstep]);
	Double y=real[fstep]*sin(PI/180*imag[fstep]);
	Complex	Z=Complex(x,y);

        //F_hat[fstep]=charge;

	// Logarithmic value of F
	F_hat[fstep]=(sign*charge*Z)/std::log(Z);

// 	std::cout<<F_hat<<std::endl;
// 	std::cout<<"This is F_hat!!"<<std::endl;
	//	getchar(); 
         
        calcAbsImped(charge, freqs[fstep], fstep, typeOut);   // calculates |Z| and writes results in File
     
        ptMyPDE_->WriteResultsInFile();     
           
        StdVector<Elem*> elemssd;
        subdoms = ptMyPDE_->getPDE_subdoms();
        //      std::cout<<"PDE_SUBDOM[0] = "<< subdoms[0]; 
        ptGrid->GetElemSD(elemssd,subdoms[0], level);
     
        if (fstep==0)
          allElemsVec.Resize(nrMeasuredData,elemssd.GetSize()*dofs*numNodes);
        for (int actEl=0;actEl<elemssd.GetSize();actEl++){
          BaseFE * ptEl = elemssd[actEl]->ptElem;
          StdVector<Integer> connecth = elemssd[actEl]->connect;

          Vector<Complex> elSolVec; 
          ptNodeStoreSol->GetElemSolution(elSolVec,connecth);
          Matrix<Double> coordinateMatrix;
          //      std::cout<<elSolVec<<std::endl;

          //    ptMyPDE_->GetElemCoords(connecth, coordinateMatrix,0);
          //    std::cout<<"\n coordinateMatrix:" <<std::endl;
          //    std::cout<<coordinateMatrix<<std::cout;

          for (int k=0;k<elSolVec.GetSize();k++)
            allElemsVec[fstep][actEl*elSolVec.GetSize()+k] = elSolVec[k];
        } // end loop over elements

    } // end of loop over all frequencies
    //std::cout<<"\n Number of Integrators in CreateF: " << ptAssemble->integrators_[0]->GetSize()<< std::endl;

    if (typeOut==true){
      //      std::cout<<"\nFinished to create F ... here it is:"<<std::endl;
      for (int i=0;i<F_hat.GetSize();i++)
        std::cout<<"F("<<i<<")="<<F_hat[i]<<"; \t";
      std::cout<<"\n ------------------------------- " <<std::endl;
    }

    //     std::cout<<" \n all ElemsVec: "<<std::endl;
    //          Complex nullC = Complex(0.0,0.0);
    //       for(int i=0;i<allElemsVec.GetSizeRow();i++)
    //            for (int j=0; j<allElemsVec.GetSizeCol();j++){
    //       //                 std::cout<<JacobiMatrix[i][j].real()<<"+"<<JacobiMatrix[i][j].imag()<<"i"<< ", ";
    //                  if (allElemsVec[i][j]!=nullC)
    //          std::cout<<"aEv'("<<i<<")("<<j<<")= "<< allElemsVec[i][j]<<"; \t";
    //          if (j==allElemsVec.GetSizeCol()-1)
    //                    std::cout<<";\n";
    //            }

  } // end createF
  // ___________________________________________________________________________________________
  //
  /////////////////////// JacobiMatrix for complex - valued parameter ///////////////////////////////
  // ___________________________________________________________________________________________


  void piezoParamIdent::createJacobiMatrixC(Matrix<Complex> & JacobiMatrix){
    ENTER_FCN("piezoParamIdent::createJacobiMatrixC");
    std::cout<<"JacobiMatrixC will be created"<<std::endl;
  
    //    Matrix<Double> * matMatrix =  ptMaterial->GetMatrix();
    ptBCs = ptMyPDE_->getPDE_BCs();                             // Pointer to BCs
    ptAlgsys = ptMyPDE_->getPDE_algsys();
    Integer level=0;
    Boolean reset = TRUE;
    ptGrid = ptMyPDE_->getPDE_grid();
    ptNodeEqn = ptMyPDE_->getPDE_eqnData();
    ptAssemble = ptMyPDE_->getPDE_assemble();
    //Integer nrMeasuredData = freqs.GetSize();
    Integer numElems_ = ptMyPDE_->getPDE_numElems();
    //    nrParameter=parameter.GetSize();
    Integer job;

    IntegratorDescriptor *actIntDescrStiff;
    IntegratorDescriptor *actIntDescrStiffC;
    
    
    Integer spaceDim = ptMyPDE_->getPDE_spaceDim();
    Double * ptsol;
    StdVector<Elem*> elemssd;
    subdoms = ptMyPDE_->getPDE_subdoms();
    ptMaterial=ptMyPDE_->getPDEMaterialData();   // Pointer to MaterialData
    ptGrid->GetElemSD(elemssd,subdoms[0], level); // gets element list elemssd

   
//     std::cout<<whichParToUpInd<<std::endl;
//     std::cout<<whichParToUpIndC<<std::endl;
//     std::cout<<"\nactNrParameters: "<<actNrParameter <<std::endl;
//     std::cout<<"\nactNrParametersC: "<<actNrParameterC <<std::endl;
    BaseNodeStoreSol * ptSol = ptMyPDE_->getPDESolution();
    NodeStoreSol<Complex> * ptNodeStoreSol;
    ptNodeStoreSol = dynamic_cast<NodeStoreSol<Complex>*>(ptSol);
   
    Vector<Complex> algSysSolVector;

    algSysSolVector=ptNodeStoreSol->GetAlgSysVector();

    Vector<Complex> RHSVec(algSysSolVector.GetSize());
   

    JacobiMatrix.Resize(2*nrMeasuredData,actNrParameter+actNrParameterC);
    if (considerMechDeformation==FALSE)
      JacobiMatrix.Resize(nrMeasuredData,actNrParameter+actNrParameterC);
   

    Integer pdenumber = 0;

    Vector<Double> dparameter(nrParameter);
    Vector<Double> dparameterC(nrParameter);

    piezoMaterialType realMatPar = REALMATERIALPARAMETER; 
    piezoMaterialType imagMatPar = IMAGMATERIALPARAMETER; 

    Vector<Double> tempHarm;
    StdVector<Integer> connect_PDE;                 

    Matrix<Double> *matMat = ptMaterial->GetMatrix();
    Matrix<Double> *matMatC = ptMaterial->GetMatrixC();

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

    //     scalingC[0]=0.001/((*matMatC)[0][0]); 
    //     scalingC[1]=0.001/((*matMatC)[2][2]);
    //     scalingC[2]=0.001/((*matMatC)[1][0]);
    //     scalingC[3]=0.001/((*matMatC)[0][2]);
    //     scalingC[4]=0.001/((*matMatC)[3][3]); 
    //     scalingC[5]=0.001/((*matMatC)[6][4]);
    //     scalingC[6]=std::abs(0.001/((*matMatC)[8][0]));
    //     scalingC[7]=0.001/((*matMatC)[8][2]);
    //     scalingC[8]=0.001/((*matMatC)[6][6]); 
    //     scalingC[9]=0.001/((*matMatC)[8][8]);

    Integer parIndex=0;


    for(int ind_param=0; ind_param<2*nrParameter;ind_param++){ // loop over different parameter increments

      // ~~~~~~~~~~~~~~~~~~~~~  second strategy ~~~~~~~~~~~~~~~~~~~~~~~~~
      if (whichParameterToUpdateRC[ind_param]==1){   
        if (ind_param<nrParameter){
          dparameter[ind_param]=relaxParameter/scaling[ind_param]*basC[parIndex].real(); // 1.0/scaling[ind_param];
        }
        else if (ind_param>=nrParameter){
	  //  std::cout<<"indParam-nrParamerer "<<ind_param-nrParameter<<std::endl;
          dparameterC[ind_param-nrParameter]=100.0*relaxParameter/scalingC[ind_param-nrParameter]*basC[parIndex-actNrParameter].real(); // 1.0/scaling[ind_param];
          //      dparameterC[ind_param-nrParameter]=1.1/scalingC[ind_param-nrParameter]*basC[parIndex-actNrParameter].imag(); // 1.0/scaling[ind_param];

        }
        //      if (ind_param>0&&ind_param<nrParameter)
        //        dparameter[whichParToUpInd[parIndex-1]]=0.0;//parameter[ind_param]; 
        //      else if(ind_param>=nrParameter){
        //                std::cout<<"Wo steige ich aus  = "<< parIndex-1-actNrParameter<<"  " << whichParToUpIndC[parIndex-1-actNrParameter]<<std::endl;
        //                dparameterC[whichParToUpIndC[parIndex-1-actNrParameter]]=1.0;// parameterC[ind_param-1-nrParameter]; 
        //      }
        //    if (ind_param>0){
        //              dparameter[ind_param-1]=0.0; 
        //      dparameter[ind_param-1]=0.0; 
        //       }      

        updateMaterialData(dparameter, ptMaterial);    

        updateComplexMaterialData(dparameterC, ptMaterial);

        //      ptMyPDE_->DefineIntegratorsWithMatInfo(level,ptMaterial);  


        if(FALSE){
          std::cout<<"\n"<<std::endl;
          for(Integer i=0;i<nrParameter;i++)
            std::cout<<dparameter[i]<<"; ";
          std::cout<<"\n"<<std::endl;
          for(Integer i=0;i<nrParameter;i++)
            std::cout<<dparameterC[i]<<"; ";
          std::cout<<"\n"<<std::endl;}     


        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        //       ptMyPDE_->DefineIntegrators(0);

        for (Integer fstep = 0; fstep < nrMeasuredData; fstep++) { // harmonic solver for different frequency - values
          reset = TRUE;

          ptAlgsys->InitRHS();

          //   ptMyPDE_-> setBCs_id_phase_(0, imag[fstep]);
        
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
            BaseForm * bilinearStiffC;

            if (spaceDim==3){
              bilinearStiff =  new linPiezo3DInt(actSDMat,isdamping);
              bilinearStiffC =  new linPiezo3DInt(actSDMat,isdamping);
            }
            else if (spaceDim==2){
              bilinearStiff = new piezoAxiInt(actSDMat, isdamping);
              bilinearStiffC = new piezoAxiInt(actSDMat, isdamping);
            }

            updateMaterialData(dparameter,ptMaterial);
            updateComplexMaterialData(dparameterC,ptMaterial);
            //      ptMyPDE_->DefineIntegratorsWithMatInfo(level,ptMaterial);  

	  
	    actIntDescrStiff = new IntegratorDescriptor(bilinearStiff, STIFFNESS);
	    actIntDescrStiffC = new IntegratorDescriptor(bilinearStiffC, STIFFNESS);


            bilinearStiff->SetPiezoMaterialType(realMatPar);
            actIntDescrStiff->SetPiezoMaterialType(realMatPar);

            bilinearStiffC->SetPiezoMaterialType(imagMatPar);
            actIntDescrStiffC->SetPiezoMaterialType(imagMatPar);

            bilinearStiff->SetElemPtr(ptEl);
            bilinearStiffC->SetElemPtr(ptEl);

            Matrix<Complex> elemMat;
            Matrix<Complex>elemMatC;
            //      Double damp_beta =1.0e-9; // in future, beta will be dependent of omega_l ...
            Double damp_beta = ptMaterial->GetDampingBeta();

            Double omega = 2.0*PI*freqs[fstep];

            bilinearStiff->CalcComplexElementMatrix(ptCoord,elemMat,damp_beta,omega);
            bilinearStiffC->CalcComplexElementMatrix(ptCoord,elemMatC,damp_beta,omega);

            if (actEl==-1){
              std::cout<<"\n"<<std::endl;
              std::cout<<elemMat<<std::endl;
              std::cout<<"\n"<<std::endl;        
              std::cout<<elemMatC<<std::endl;
              std::cout<<"\n"<<std::endl;        
            }

            elemMat+=elemMatC;

            Complex sum=Complex(0.0,0.0);
            Vector<Complex> temp;
            temp.Resize(elemMat.GetSizeRow());

            for (int i=0; i<elemMat.GetSizeRow();i++){
              for (int j=0; j<elemMat.GetSizeCol();j++)
                sum=sum+elemMat[i][j]*allElemsVec[fstep][actEl*elemMat.GetSizeRow()+j];  
              temp[i]=sum;
              sum=Complex(0.0,0.0);
            }

            // save values like it was done in transformElemMat2harmonic (see assemle.cc) ...
            tempHarm.Resize(2*temp.GetSize());
            for (int i=0; i<temp.GetSize();i++){
              tempHarm[i]=-temp[i].real();
              tempHarm[i+temp.GetSize()]=-temp[i].imag();
            } 

            ptAlgsys->SetElementRHS(&tempHarm[0], pdeId_, connect_PDE.GetPointer(), 
                                    connect_PDE.GetSize());

	    //  if (spaceDim==3){
	                delete  bilinearStiff;
	                delete bilinearStiffC;
          } // end for over all Elems

	  
          updateMaterialData(parameter,ptMaterial);
          updateComplexMaterialData(parameterC,ptMaterial);
          //      ptMyPDE_->DefineIntegratorsWithMatInfo(level,ptMaterial);  

          // the following lines are the content of: PreStepHarmonic:
          // ptMyPDE_->PreStepHarmonic(fstep, freqs[fstep], level, reset);

          ptMyPDE_->setPDE_actFrequency(freqs[fstep]);
          ptMyPDE_->setPDE_actFreqStep(fstep);
          ptAssemble->SetFrequency(freqs[fstep]);


          if (reset){
	    ptAlgsys->InitMatrix();
	    ptAssemble->SetReassemble();
	  }
	  
          //      ptMyPDE_->DefineIntegratorsWithMatInfo(level,ptMaterial);  
        
          //    Cannot use SolveStepHarmonic, since it overwrites RHS ...
          //      for this, I have copied the method StepHarmonicLin to this place 
          //      void BasePDE::StepHarmonicLin(const Integer freqStep, const Double frequency, const Integer level, const Boolean reset)
        
          ptAssemble->AssembleMatrices(level);

               
          // The folowing method creates and calculates the RHS for harmonic Problems in CreateJacobian ...
                
          ptAssemble->AssembleSrcRHS(level, freqs[fstep]);

          ptNodeStoreSol = dynamic_cast<NodeStoreSol<Complex>*>(ptSol);
                
          if (reset)
            {
              //account for bcs
              ptMyPDE_->SetBCs(level, freqs[fstep]);
              job = 1; // calc new preconditioner
            }
          else
            job = 3;

          ptAlgsys->BuildInDirichlet();
          if ( job == 1 ) {
            ptAlgsys->SetupSolver(job);
            ptAlgsys->SetupPrecond(job);
          }

          ptAlgsys->Solve();

          ptAlgsys->GetSolutionVal(ptsol);
          ptSol->CopyFromAlgSysDataPointer(ptsol);

          if (considerMechDeformation==TRUE){
            //      ptNodeStoreSol->GetGlobalSolVector(ELEC_POTENTIAL, solElecPot);
            ptNodeStoreSol->GetGlobalSolVector(MECH_DISPLACEMENT, solMechDispl);
            //    typeOutSolutionOnScreen(solElecPot, solMechDispl);              //member function of piezoParamIdent
            measureMechDeformationInZ_Direction(solMechDispl,radius,meanValueMechDeformation,dofs); // Braucht üÜberarbeitung!!
            JacobiMatrix[fstep + nrMeasuredData][parIndex]=meanValueMechDeformation;    
          }

          ptMyPDE_->GetSolveStep()->PostStepHarmonic(fstep, freqs[fstep], level, reset);
          ptMyPDE_->PostProcess(level);

          Vector<Complex> chargeVec =  ptMyPDE_->getPDE_complexValuedCharge();

          Complex charge=Complex(0.0,0.0);

          for (int i=0; i<chargeVec.GetSize();i++)
            charge=charge+chargeVec[i];
    
          JacobiMatrix[fstep][parIndex]=sign*charge; //ptMyPDE_->getPDE_complexValuedCharge();    
     
        }               // end for over all frequencies
        if (ind_param<nrParameter)
          dparameter[ind_param]=0.0;// parameter[ind_param]; 

        if (ind_param>=nrParameter)
          dparameterC[ind_param-nrParameter]=0.0;// parameterC[ind_param-1-nrParameter];

        parIndex++;
      } // end if whichParameterToUpdateRC==1

    } //end loop over paramters

    //     for(int i=0;i<JacobiMatrix.GetSizeRow();i++)
    //       for (int j=0; j<JacobiMatrix.GetSizeCol();j++){
    //  //      std::cout<<JacobiMatrix[i][j].real()<<"+"<<JacobiMatrix[i][j].imag()<<"i"<< ", ";
    //          std::cout<<"F'("<<i<<")("<<j<<")= "<< JacobiMatrix[i][j]<<"; \t";
    //          if (j==JacobiMatrix.GetSizeCol()-1)
    //   std::cout<<";\n";
    //       }

    //     std::cout<<JacobiMatrix<<std::endl;

    //        std::cout<<"\n end CreateJacobiMatrix C"<<std::endl;

  }            //end CreateJacobiMatrix 2

  // ___________________________________________________________________________________________
  //
  /////////////////////// JacobiMatrix for real-valued parameter ///////////////////////////////
  // ___________________________________________________________________________________________


  void piezoParamIdent::createJacobiMatrix2(Matrix<Complex> & JacobiMatrix){
    ENTER_FCN("piezoParamIdent::createJacobiMatrix2");
    //    std::cout<<"JacobiMatrix2 will be created"<<std::endl;
  
    //    Matrix<Double> * matMatrix =  ptMaterial->GetMatrix();
    ptBCs = ptMyPDE_->getPDE_BCs();                             // Pointer to BCs
    ptAlgsys = ptMyPDE_->getPDE_algsys();
    Integer level=0;
    Boolean reset = TRUE;
    ptGrid = ptMyPDE_->getPDE_grid();
    ptNodeEqn = ptMyPDE_->getPDE_eqnData();
    ptAssemble = ptMyPDE_->getPDE_assemble();
    //Integer nrMeasuredData = freqs.GetSize();
    Integer numElems_ = ptMyPDE_->getPDE_numElems();
    //    nrParameter=parameter.GetSize();
    Integer job;
    
    Integer spaceDim = ptMyPDE_->getPDE_spaceDim();
    Double * ptsol;
    StdVector<Elem*> elemssd;
    subdoms = ptMyPDE_->getPDE_subdoms();
    ptMaterial=ptMyPDE_->getPDEMaterialData();   // Pointer to MaterialData
    ptGrid->GetElemSD(elemssd,subdoms[0], level); // gets element list elemssd

    BaseNodeStoreSol * ptSol = ptMyPDE_->getPDESolution();
    NodeStoreSol<Complex> * ptNodeStoreSol;
    ptNodeStoreSol = dynamic_cast<NodeStoreSol<Complex>*>(ptSol);
    Vector<Complex> algSysSolVector;
    algSysSolVector=ptNodeStoreSol->GetAlgSysVector();
    Vector<Complex> RHSVec(algSysSolVector.GetSize());

    //    JacobiMatrix.Resize(2*nrMeasuredData,actNrParameter);
    //   if (considerMechDeformation==FALSE)
    JacobiMatrix.Resize(nrMeasuredData,actNrParameter);

    Integer pdenumber = 0;

    Vector<Double> dparameter(nrParameter);
  
    StdVector<Integer> connect_PDE;                 
    ptMyPDE_->WriteGeneralPDEdefines();

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
    Integer parIndex=0;

    for(int ind_param=0; ind_param<nrParameter;ind_param++){ // loop over different parameter increments
      if (whichParameterToUpdate[ind_param]==1){
      
        // ~~~~~~~~~~~~~~~~~~~~~  second strategy ~~~~~~~~~~~~~~~~~~~~~~~~~
        //              dparameter[ind_param]=1.2/scaling[ind_param];
        //     dparameter[ind_param]=1.3/scaling[ind_param]*bas[ind_param]; // 1.0/scaling[ind_param];
 	if (whichNewtonCG==3)
           dparameter[ind_param]=relaxParameter/scaling[ind_param]*bas[parIndex]; // 1.0/scaling[ind_param];
         else 
          dparameter[ind_param]=relaxParameter/scaling[ind_param];
	
	//        std::cout<<whichParToUpInd<<std::endl;
	//        std::cout<<"\n parIndex "<<parIndex<<", ind_param: "<<ind_param <<", whichParToUpInd[parIndex] " <<whichParToUpInd[parIndex] <<std::endl;

      
	if (parIndex>0){
	  dparameter[whichParToUpInd[parIndex-1]]=0.0; 
	  //dparameter[ind_param-1]=0.0;
	}

	updateMaterialData(dparameter, ptMaterial);  
	// ptMyPDE_->DefineIntegratorsWithMatInfo(level,ptMaterial);  

	if(FALSE){
	  std::cout<<"\n"<<std::endl;
	  for(Integer i=0;i<parameter.GetSize();i++)
	    std::cout<<dparameter[i]<<"; ";
	  std::cout<<"\n"<<std::endl;
	}

      
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	//       ptMyPDE_->DefineIntegrators(0);

	for (Integer fstep = 0; fstep < nrMeasuredData; fstep++) { // harmonic solver for different frequency - values
	  reset = TRUE;

	  ptAlgsys->InitRHS();
	  ptAlgsys->InitMatrix();
	  ptAssemble->SetReassemble();
        

	  //   ptMyPDE_-> setBCs_id_phase_(0, imag[fstep]);
        
	  //loop over elements
	  for (int actEl=0; actEl< elemssd.GetSize(); actEl++) {
	    BaseFE * ptEl = elemssd[actEl]->ptElem;
	    StdVector<Integer> connecth = elemssd[actEl]->connect;
                             
	    Matrix<Double> ptCoord;
	    ptGrid->GetCoordNodesElemMat(connecth, ptCoord, 0);
	    //        std::cout<<ptCoord<<std::endl;
	    //  std::cout<<connecth<<std::endl;
                    
	    // map connect to PDE node numbers

	    ptNodeEqn->Node2EQN(connecth, connect_PDE);
	    //        std::cout<<"nodeToEquation:"<<std::endl;
	    //        std::cout<<connect_PDE<<std::endl;

	    //Vector<Complex> elSolVec; 
	    // ptNodeStoreSol->GetElemSolution(elSolVec,connecth);

	    Boolean isdamping=TRUE;
             
	    // Create new Integrator, with this calculate elemMat for sysmat in RHS 
	    BaseForm * bilinearStiff;

	    //        updateMaterialData(dparameter,ptMaterial);      

	    if (spaceDim==3)
	      bilinearStiff =  new linPiezo3DInt(*ptMaterial,isdamping);
	    else if (spaceDim==2)     
	      bilinearStiff = new piezoAxiInt(*ptMaterial, isdamping);

	    //                std::cout<<"matMat 2"<<std::endl;
	    //        matMat = ptMaterial->GetMatrix();
	    //        std::cout<<*matMat<<std::endl;

	    updateMaterialData(dparameter,ptMaterial);

	    //        std::cout<<"matMat 3"<<std::endl;
	    //        matMat = ptMaterial->GetMatrix();
	    //        std::cout<<*matMat<<std::endl;

          
	    //        IntegratorDescriptor *actIntDescrStiff = new IntegratorDescriptor(bilinearStiff, STIFFNESS);

	    bilinearStiff->SetElemPtr(ptEl);
	    Matrix<Complex> elemMat;
	    //      Double damp_beta =1.0e-9; // in future, beta will be dependent of omega_l ...
	    Double damp_beta = ptMaterial->GetDampingBeta();

	    if (params->IsSet("adjustDamping",  "harmonic"))
		if (freqs[0] > 0 && fstep > 0 ) {
		  damp_beta = damp_beta*freqs[0] / freqs[fstep];
		  //		  std::cout<<"damp-beta = " << damp_beta << std::endl;
		}
		

	    Double omega = 2.0*PI*freqs[fstep];
	    bilinearStiff->CalcComplexElementMatrix(ptCoord,elemMat,damp_beta,omega);
	    //        std::cout<<"\n dparam = " << ind_param<<"; actEl = "<<actEl<< std::endl;
	    //                std::cout<<elemMat<<std::endl;
	    //        std::cout<<"\n Size of allElemsVec = "<<allElemsVec.GetSizeCol()<<std::endl;

	    Complex sum=Complex(0.0,0.0);
	    Vector<Complex> tempSolution;
	    tempSolution.Resize(elemMat.GetSizeRow());
	    for (int i=0; i<elemMat.GetSizeRow();i++){
	      for (int j=0; j<elemMat.GetSizeCol();j++)
		sum=sum+elemMat[i][j]*allElemsVec[fstep][actEl*elemMat.GetSizeRow()+j];  
	      tempSolution[i]=sum;
	      sum=Complex(0.0,0.0);
	    }
	    // save values like it was done in transformElemMat2harmonic (see assemle.cc) ...
	    Vector<Double> tempHarm;
	    tempHarm.Resize(2*tempSolution.GetSize());
	    for (int i=0; i<tempSolution.GetSize();i++){
	      tempHarm[i]=-tempSolution[i].real();
	      tempHarm[i+tempSolution.GetSize()]=-tempSolution[i].imag();
	    } 
	    //        std::cout<<tempHarm<<std::endl;

        
	    ptAlgsys->SetElementRHS(&tempHarm[0], pdeId_, connect_PDE.GetPointer(), 
				    connect_PDE.GetSize());
	    delete bilinearStiff;  

	  } // end for over all Elems

	  updateMaterialData(parameter,ptMaterial);
	  //          ptMyPDE_->DefineIntegratorsWithMatInfo(level,ptMaterial);

	  // the following lines are the content of: PreStepHarmonic:
	  // ptMyPDE_->PreStepHarmonic(fstep, freqs[fstep], level, reset);

	  ptMyPDE_->setPDE_actFrequency(freqs[fstep]);
	  ptMyPDE_->setPDE_actFreqStep(fstep);
	  ptAssemble->SetFrequency(freqs[fstep]);
	  //        updateMaterialData(parameter, ptMaterial);    // is neccessary, since otherwise we wouldd solve PDE with sparse Mat-Data
	  //        ptMyPDE_->DefineIntegratorsWithMatInfo(level,ptMaterial);

	  if (reset) {
	    ptAlgsys->InitMatrix();
	    ptAssemble->SetReassemble();
	  }
	  //        ptMyPDE_->DefineIntegratorsWithMatInfo(level,ptMaterial);
        
	  //        Cannot use SolveStepHarmonic, since it overwrites RHS ...
	  //      for this, I have copied the method StepHarmonicLin to this place 
	  //      void BasePDE::StepHarmonicLin(const Integer freqStep, const Double frequency, const Integer level, const Boolean reset)
        
	  ptAssemble->AssembleMatrices(level);
               
	  // The folowing method creates and calculates the RHS for harmonic Problems in CreateJacobian ...
                
	  ptAssemble->AssembleSrcRHS(level, freqs[fstep]);

	  ptNodeStoreSol = dynamic_cast<NodeStoreSol<Complex>*>(ptSol);
                
	  if (reset)
	    {
	      //account for bcs
	      ptMyPDE_->SetBCs(level, freqs[fstep]);
	      job = 1; // calc new preconditioner
	    }
	  else
	    job = 3;

	  ptAlgsys->BuildInDirichlet();
	  if ( job == 1 ) {
	    ptAlgsys->SetupPrecond(job);
	    ptAlgsys->SetupSolver(job);
	  }

	  ptAlgsys->Solve();
	  ptAlgsys->GetSolutionVal( ptsol );
	  ptSol->CopyFromAlgSysDataPointer(ptsol);

	  if (considerMechDeformation==TRUE){
	    //      ptNodeStoreSol->GetGlobalSolVector(ELEC_POTENTIAL, solElecPot);
	    ptNodeStoreSol->GetGlobalSolVector(MECH_DISPLACEMENT, solMechDispl);
	    //        typeOutSolutionOnScreen(solElecPot, solMechDispl);              //member function of piezoParamIdent
	    measureMechDeformationInZ_Direction(solMechDispl,radius,meanValueMechDeformation,dofs); // Braucht üÜberarbeitung!!
	    JacobiMatrix[fstep + nrMeasuredData][ind_param]=meanValueMechDeformation;    
	  }

	  ptMyPDE_->GetSolveStep()->PostStepHarmonic(fstep, freqs[fstep], level, reset);
	  ptMyPDE_->PostProcess(level);
        
	  Vector<Complex> chargeVec =  ptMyPDE_->getPDE_complexValuedCharge();

	  Complex charge=Complex(0.0,0.0);

	  for (int i=0; i<chargeVec.GetSize();i++)
	    charge=charge+chargeVec[i];
    
	  JacobiMatrix[fstep][parIndex]=sign*charge; //ptMyPDE_->getPDE_complexValuedCharge();    
     
	}           // end for over all frequencies

	if (ind_param==nrParameter-1){
	  dparameter[ind_param-1]=0.0;
	}
	parIndex++;
      } // end if whichParameterToUpdate==1

    } //end loop over paramters

    if(FALSE){
      for(int i=0;i<JacobiMatrix.GetSizeRow();i++)
        for (int j=0; j<JacobiMatrix.GetSizeCol();j++){
          //std::cout<<JacobiMatrix[i][j].real()<<"+"<<JacobiMatrix[i][j].imag()<<"i"<< ", ";
          std::cout<<"F'("<<i<<")("<<j<<")= "<< JacobiMatrix[i][j]<<"; \t";
          //    std::cout<<std::setprecision(10);
          if (j==JacobiMatrix.GetSizeCol()-1)
            std::cout<<";\n";
        }
    }

  }            //end CreateJacobiMatrix 2


  void piezoParamIdent::testJacobiMatrix(Vector<Complex> & F_hat, Matrix<Complex> & JacobiMatrix, Vector<Double> & parameter,BCs * ptBCs,MaterialData * ptMaterial,Vector<Double> & parameterIncrement, Vector<Complex>& solElecPot,Vector<Complex> &solMechDispl){
    ENTER_FCN("piezoParamIdent::testJacobiMatrix");

    Vector<Complex> F_hat_incr(F_hat.GetSize());
    approxJacobiMatrix.Resize(JacobiMatrix.GetSizeRow(), JacobiMatrix.GetSizeCol());
    Vector<Double> parameter_incr(parameter.GetSize());
    parameter_incr=parameter;

    //     for(Integer i=0;i<parameter.GetSize();i++){
    //       parameter_incr[i]=1.1*parameter[i];
    //       //      std::cout<<"Denominator("<<i<<") = " << (Complex(parameter[i])-Complex(parameter_incr[i]))<<std::endl;
    //     }

    //     updateMaterialData(parameter_incr, ptMaterial);
    updateMaterialData(parameter, ptMaterial);
    createF(ptMaterial, ptBCs, F_hat, FALSE);
    //     std::cout<<"\n"<<std::endl;
    //     //    for (int i=0;i<F_hat.GetSize();i++)
    //     //      std::cout<<"F("<<i+1<<")="<<F_hat[i]<< " <-> " << F_hat_incr[i]<<"; ";
    //     // createJacobiMatrix(ptMaterial, ptBCs, F_hat_incr, parameterIncrement,JacobiMatrix, solElecPot, solMechDispl);
    //     //    std::cout<<"\n  - - - - - - - - - - - - - - - \n JacobiMatrix <-> approxJacobiMatrix"<< std::endl;

    //     std::cout<<F_hat<<std::endl;
    //     std::cout<<F_hat_incr<<std::endl;
    //     for (Integer i=0; i< JacobiMatrix.GetSizeRow();i++)
    //       for (Integer j=0; j< JacobiMatrix.GetSizeCol();j++){
    //  if (j<parameter.GetSize())
    //    approxJacobiMatrix[i][j]=(F_hat_incr[i]-F_hat[i])/(scaling[i]*(parameter_incr[j]-parameter[i]));
    //  else {
    //    Integer jj=i-parameter.GetSize();
    //    approxJacobiMatrix[i][j]=(F_hat_incr[i]-F_hat[i])/(Complex(scaling[jj]*parameter[jj])-Complex(scaling[jj]*parameter_incr[jj]));
    //  }
    //  std::cout<<"F~("<<i<<")("<<j<<")= "<< approxJacobiMatrix[i][j]<<";\t "; // <<" <-> " <<approxJacobiMatrix[i][j]<< "; ";
    //  //      if (j==JacobiMatrix.GetSizeCol()-1)
    //  //std::cout<<"\n " ;
    //       }
    //     //       JacobiMatrix = approxJacobiMatrix;
    //     updateMaterialData(parameter,ptMaterial);

    for (Integer ind_param=0;ind_param<nrParameter;ind_param++){

      parameter_incr[ind_param]=1.0001*parameter[ind_param];
      //      std::cout<<parameter_incr[ind_param]<<std::endl;
      updateMaterialData(parameter_incr,ptMaterial);
      createF(ptMaterial,ptBCs,F_hat_incr,FALSE);

      for (Integer j=0;j<nrMeasuredData;j++)
        approxJacobiMatrix[j][ind_param]=-(F_hat[j]-F_hat_incr[j])/((parameter_incr[ind_param]-parameter[ind_param])*scaling[ind_param]);

      parameter_incr[ind_param]=parameter[ind_param];
    }
//     std::cout<<"\n Here we see the approx. Jacobian and the created Jacobian Matrix:"<<std::endl;
//     std::cout<<approxJacobiMatrix<<std::endl;
//     std::cout<<JacobiMatrix<<std::endl;
    // getchar();   

  }// end testJacobiMatrix

void piezoParamIdent::testJacobiMatrix2(Vector<Complex> & F_hat, Matrix<Complex> & JacobiMatrix, Vector<Double> & parameter,BCs * ptBCs,MaterialData * ptMaterial,Vector<Double> & parameterIncrement, Vector<Complex>& solElecPot,Vector<Complex> &solMechDispl){
    ENTER_FCN("piezoParamIdent::testJacobiMatrix");

    Vector<Complex> F_hat_incr(F_hat.GetSize());
    Vector<Complex> F_hat_incr2(F_hat.GetSize());
    Vector<Complex> F_hat_incr3(F_hat.GetSize());
    Vector<Complex> F_hat_incr4(F_hat.GetSize());
    approxJacobiMatrix.Resize(nrMeasuredData,actNrParameter);
    Vector<Double> parameter_incr(nrParameter);
    Vector<Double> parameter_incr2(nrParameter);
    Vector<Double> parameter_incr3(nrParameter);
    Vector<Double> parameter_incr4(nrParameter);

    parameter_incr=parameter;
    parameter_incr2=parameter;
    Integer parInd=0;

    updateMaterialData(parameter, ptMaterial);
    createF(ptMaterial, ptBCs, F_hat, FALSE);

    for (Integer ind_param=0;ind_param<nrParameter;ind_param++){ 
      if (whichParameterToUpdate[ind_param]==1){

	parameter_incr[ind_param]=1.001*parameter[ind_param];
	//	std::cout<<parameter_incr<<std::endl
	updateMaterialData(parameter_incr,ptMaterial);
	createF(ptMaterial,ptBCs,F_hat_incr,FALSE);

	parameter_incr2[ind_param]=0.999*parameter[ind_param];	
	//	std::cout<<parameter_incr2<<std::endl;
	updateMaterialData(parameter_incr2,ptMaterial);
	createF(ptMaterial,ptBCs,F_hat_incr2,FALSE);

	//	parameter_incr3[ind_param]=1.005*parameter[ind_param];
	//	std::cout<<parameter_incr<<std::endl
	//updateMaterialData(parameter_incr3,ptMaterial);
	//	createF(ptMaterial,ptBCs,F_hat_incr3,FALSE);

	//	parameter_incr4[ind_param]=0.995*parameter[ind_param];	
	//	std::cout<<parameter_incr2<<std::endl;
	//	updateMaterialData(parameter_incr4,ptMaterial);
	//	createF(ptMaterial,ptBCs,F_hat_incr4,FALSE);



      // First order approximation
//       for (Integer j=0;j<nrMeasuredData;j++)
//         approxJacobiMatrix[j][ind_param]=-(F_hat[j]-F_hat_incr[j])/((parameter_incr[ind_param]-parameter[ind_param])*scaling[ind_param]);

      // second order FD approximation
  	for (Integer j=0;j<nrMeasuredData;j++)
  	  approxJacobiMatrix[j][parInd]=0.5*(F_hat_incr[j]-F_hat_incr2[j])/((parameter_incr[ind_param]-parameter[ind_param])*scaling[ind_param]);


      // forth order FD approximation
//  	for (Integer j=0;j<nrMeasuredData;j++)
//  	  approxJacobiMatrix[j][parInd]=1.0/6.0*(8.0*F_hat_incr3[j]-8.0*F_hat_incr4[j]-F_hat_incr[j]+F_hat_incr2[j])/((parameter_incr[ind_param]-parameter[ind_param])*scaling[ind_param]);


	    parInd++;
	    //	std::cout<<"\n Performed second order FD Approx of Jacobian"<<std::endl;
	parameter_incr[ind_param]=parameter[ind_param];
	parameter_incr2[ind_param]=parameter[ind_param];
	parameter_incr3[ind_param]=parameter[ind_param];
	parameter_incr4[ind_param]=parameter[ind_param];

      }
    }
    //    std::cout<<"\n Here we see the approx. Jacobian and the created Jacobian Matrix:"<<std::endl;
  
    // JacobiMatrix=approxJacobiMatrix;

  }// end testJacobiMatrix


  void piezoParamIdent::testJacobiMatrixC(Vector<Complex> & F_hat, Matrix<Complex> & JacobiMatrix, 
					  Vector<Double> & parameter,BCs * ptBCs,MaterialData * ptMaterial){
    ENTER_FCN("piezoParamIdent::testJacobiMatrix");

    Vector<Complex> F_hat_incr(F_hat.GetSize());
    Vector<Complex> F_hat_incr2(F_hat.GetSize());
    Vector<Complex> F_hat_incr3(F_hat.GetSize());
    Vector<Complex> F_hat_incr4(F_hat.GetSize());
    approxJacobiMatrix.Resize(nrMeasuredData,actNrParameter+actNrParameterC);
    Vector<Double> parameter_incr(nrParameter);
    Vector<Double> parameter_incr2(nrParameter);
    Vector<Double> parameter_incr3(nrParameter);
    Vector<Double> parameter_incr4(nrParameter);

    parameter_incr=parameter;
    parameter_incr2=parameter;
    Integer parInd=0;

    updateMaterialData(parameter, ptMaterial);
    createF(ptMaterial, ptBCs, F_hat, FALSE);

    for (Integer ind_param=0;ind_param<nrParameter;ind_param++){ 
      if (whichParameterToUpdate[ind_param]==1){

	parameter_incr[ind_param]=1.001*parameter[ind_param];
	//	std::cout<<parameter_incr<<std::endl
	updateMaterialData(parameter_incr,ptMaterial);
	createF(ptMaterial,ptBCs,F_hat_incr,FALSE);

	parameter_incr2[ind_param]=0.999*parameter[ind_param];	
	//	std::cout<<parameter_incr2<<std::endl;
	updateMaterialData(parameter_incr2,ptMaterial);
	createF(ptMaterial,ptBCs,F_hat_incr2,FALSE);

      // second order FD approximation
  	for (Integer j=0;j<nrMeasuredData;j++)
  	  approxJacobiMatrix[j][parInd]=0.5*(F_hat_incr[j]-F_hat_incr2[j])/
	    ((parameter_incr[ind_param]-parameter[ind_param])*scaling[ind_param]);

	    parInd++;
	    //	std::cout<<"\n Performed second order FD Approx of Jacobian"<<std::endl;
	parameter_incr[ind_param]=parameter[ind_param];
	parameter_incr2[ind_param]=parameter[ind_param];
	parameter_incr3[ind_param]=parameter[ind_param];
	parameter_incr4[ind_param]=parameter[ind_param];

      }
    }

    parInd=0;
    for (Integer ind_param=0;ind_param<nrParameter;ind_param++){ 
      if (whichParameterToUpdateC[ind_param]==1){

	parameter_incr[ind_param]=1.001*parameterC[ind_param];
	//	std::cout<<parameter_incr<<std::endl

	updateComplexMaterialData(parameter_incr,ptMaterial);
	createF(ptMaterial,ptBCs,F_hat_incr,FALSE);

	parameter_incr2[ind_param]=0.999*parameterC[ind_param];	
	//	std::cout<<parameter_incr2<<std::endl;
	updateComplexMaterialData(parameter_incr2,ptMaterial);
	createF(ptMaterial,ptBCs,F_hat_incr2,FALSE);

      // second order FD approximation
  	for (Integer j=0;j<nrMeasuredData;j++)
  	  approxJacobiMatrix[j][actNrParameter+parInd]=0.5*(F_hat_incr[j]-F_hat_incr2[j])/
	    ((parameter_incr[ind_param]-parameter[ind_param])*scaling[ind_param]);

	    parInd++;
	    //	std::cout<<"\n Performed second order FD Approx of Jacobian"<<std::endl;
	parameter_incr[ind_param]=parameterC[ind_param];
	parameter_incr2[ind_param]=parameterC[ind_param];
	parameter_incr3[ind_param]=parameterC[ind_param];
	parameter_incr4[ind_param]=parameterC[ind_param];

      }
    }
        
  } // end testJacobiMatrix



  void piezoParamIdent::createAdjointJacobiMatrix(Matrix<Complex> & JacobiMatrix, Matrix<Complex> & adjJacobiMatrix){
    ENTER_FCN("piezoParamIdent::createAdjointJacobiMatrix");
    //    std::cout<<"\n Adjoint Jacobian will be created ... "<<std::endl;
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
    Integer spaceDim = ptMyPDE_->getPDE_spaceDim();
    Double * ptsol;
    StdVector<Elem*> elemssd;
    subdoms = ptMyPDE_->getPDE_subdoms();
    ptMaterial=ptMyPDE_->getPDEMaterialData();   // Pointer to MaterialData
    ptGrid->GetElemSD(elemssd,subdoms[0], level); // gets element list elemssd

    BaseNodeStoreSol * ptSol = ptMyPDE_->getPDESolution();
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

      ptAlgsys->SetElementRHS(&tempHarm[0], pdeId_, connect_PDE.GetPointer(), 
                              connect_PDE.GetSize());

    } // end for elemssd 

  } // end createAndSetRHSforJacobian();




  void piezoParamIdent::createJacobiMatrix(MaterialData * ptMaterial, BCs * ptBCs, Vector<Complex> & F_hat, Vector<Double> & parameterIncrement, Matrix<Complex> & JacobiMatrix, Vector<Complex> & solElecPot,Vector<Complex> & solMechDispl){
    ENTER_FCN("piezoParamIdent::createJacobiMatrix");
    std::cout<<"JacobiMatrix will be created"<<std::endl;
    Vector<Double> IncrementedRHSMatrix;   
  
    //    Matrix<Double> * matMatrix =  ptMaterial->GetMatrix();
    ptBCs = ptMyPDE_->getPDE_BCs();                             // Pointer to BCs
    ptAlgsys = ptMyPDE_->getPDE_algsys();
    Integer level=0;
    Boolean reset = TRUE;
    ptGrid = ptMyPDE_->getPDE_grid();
    ptNodeEqn = ptMyPDE_->getPDE_eqnData();
    ptAssemble = ptMyPDE_->getPDE_assemble();
    //Integer nrMeasuredData = freqs.GetSize();
    Integer numElems_ = ptMyPDE_->getPDE_numElems();
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
      //        dparameter[ind_param-1]=0.0; 
      //        updateMaterialData(dparameter, ptMaterial);    
      //        if (ind_param==nrParameter-1)
      //        dparameter[ind_param-1]=0.0; 

      //       std::cout<<"\n"<<std::endl;
      //       for(Integer i=0;i<parameter.GetSize();i++)
      //                std::cout<<dparameter[i]<<"; ";

      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



      //      ptMyPDE_->DefineIntegrators(0);

      for (Integer fstep = 0; fstep < nrMeasuredData; fstep++) { // harmonic solver for different frequency - values
        reset = TRUE;

        ptMyPDE_-> setBCs_id_phase_(0, imag[fstep]);

        //      Info->WriteHarmonicStep(ptMyPDE_->GetName(), fstep, freqs[fstep]);
        
        ptMyPDE_->WriteGeneralPDEdefines();
        ptMyPDE_->GetSolveStep()->PreStepHarmonic(fstep, freqs[fstep], level, reset);
        
        //      Cannot use SolveStepHarmonic, since it overwrites RHS ...
        //      for this, I have copied the method StepHarmonicLin to this place 
        //      void BasePDE::StepHarmonicLin(const Integer freqStep, const Double frequency, const Integer level, const Boolean reset)
        
        ptAssemble->AssembleMatrices(level);
       
        // The folowing method creates and calculates the RHS for harmonic Problems in CreateJacobian ...
        createAndSetRHSforJacobian(fstep);

        //this has to be done each time!
        ptAssemble->AssembleSrcRHS(level, freqs[fstep]);

        Double * ptsol;
        BaseNodeStoreSol * ptSol = ptMyPDE_->getPDESolution();
        NodeStoreSol<Complex> * ptNodeStoreSol;
        ptNodeStoreSol = dynamic_cast<NodeStoreSol<Complex>*>(ptSol);
                
        ptMyPDE_-> setBCs_id_phase_(0, imag[fstep]);

        if (reset)
          {
            //account for bcs
            ptMyPDE_->SetBCs(level, freqs[fstep]);
            job = 1; // calc new preconditioner
          }
        else
          job = 3;
        //        std::cout<<"SetBcs ..."<<std::endl;

        ptAlgsys->BuildInDirichlet();
        if ( job == 1 ) {
          ptAlgsys->SetupPrecond(job);
          ptAlgsys->SetupSolver(job);
        }

        ptAlgsys->Solve();
        ptAlgsys->GetSolutionVal( ptsol );
        ptSol->CopyFromAlgSysDataPointer(ptsol);



        if (considerMechDeformation==TRUE){       
          //    ptNodeStoreSol->GetGlobalSolVector(ELEC_POTENTIAL, solElecPot);
          ptNodeStoreSol->GetGlobalSolVector(MECH_DISPLACEMENT, solMechDispl);
          //      typeOutSolutionOnScreen(solElecPot, solMechDispl);              //member function of piezoParamIdent
          measureMechDeformationInZ_Direction(solMechDispl,radius,meanValueMechDeformation,dofs); // Braucht üÜberarbeitung!!
          JacobiMatrix[fstep + nrMeasuredData][ind_param]=meanValueMechDeformation;    
        }

        ptMyPDE_->GetSolveStep()->PostStepHarmonic(fstep, freqs[fstep], level, reset);
        ptMyPDE_->PostProcess(level);
        ptMyPDE_->GetSolveStep()->PostStepHarmonic(fstep, freqs[fstep], level, reset);
        Vector<Complex> chargeVec =  ptMyPDE_->getPDE_complexValuedCharge();

        Complex charge=Complex(0.0,0.0);


        for (int i=0; i<chargeVec.GetSize();i++){
          charge=charge+chargeVec[i];

        }
       
        //      mean_value_charge = mean_value_charge/(Double(chargeVec.GetSize()));
        //      std::cout<<"\nCHARGE VEC SIZE = "<< chargeVec.GetSize()<<"mean-value-charge = " << mean_value_charge << std::endl;

        JacobiMatrix[fstep][ind_param]=sign*charge; //ptMyPDE_->getPDE_complexValuedCharge();    


        //   Info->PrintPiezoMat(*ptMaterial);
        //     ptMyPDE_->WriteResultsInFile();     
     
      }         // end for over all frequencies
      std::cout << " \nCreateJacobian Line " << ind_param  <<std::endl;
      //for (int i=0;i<JacobiMatrix.GetSizeRow();i++)
      //        std::cout<<"F'("<<i<<")("<<ind_param<<")= "<< JacobiMatrix[i][ind_param]<<"; ";


    } //end loop over paramters
    
    std::cout<<"JACOBI - MATRIX 1: " <<std::endl;
    for(int i=0;i<JacobiMatrix.GetSizeRow();i++)
      for (int j=0; j<JacobiMatrix.GetSizeCol();j++){
        //      std::cout<<JacobiMatrix[i][j].real()<<"+"<<JacobiMatrix[i][j].imag()<<"i"<< ", ";
        std::cout<<"F'("<<i<<")("<<j<<")= "<< JacobiMatrix[i][j]<<"; \t";
        if (j==JacobiMatrix.GetSizeCol()-1)
          std::cout<<";\n";
      }

  }            //end CreateJacobiMatrix



} // end namespace
