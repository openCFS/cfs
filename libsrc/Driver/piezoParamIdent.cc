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
#include "DataInOut/MaterialData.hh"
#include "PDE/timestepping.hh"
#include "Utils/baseelemstoresol.hh"
#include "singleDriver.hh"
//#include "/../OLAS/algsys/basesystem.hh"
//#include "DataInOut/piezoParameterData.hh"


#include <fstream>
#include <iostream>
#include <string>


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

    Char* measuredData="measuredData.dat";
    allMeasuredData.open("measuredData.dat");

    if (!allMeasuredData.good())
      {
	std::cerr << "File measuredData.dat does not exist!" << std::endl;
	exit(1);
      }
    else
      std::cerr <<"File measuredData is opened to be read" << std::endl;
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
    std::cout<<"\n Calc parameter() in piezoParamIdent ... \n";

    parameter.Resize(10);
    parameterIncrement.Resize(10);
    parameterIncrement = parameter;
    omegas.Resize(15);
    freqs.Resize(15);
    real.Resize(15);
    imag.Resize(15);
    amplitude_phase.Resize(15);
    F_hat.Resize(15);

    // the following passage reads Data from file measuredData.dat
    // The rows are containing the values of the given frequencies, such as phase and amplitude!
    readMeasuredData(freqs, real, imag, parameter, voltage, nrMeasuredData, thickness, radius, delta);
    std::cout<<"Size of piezoElectric Body:"<< thickness << " x " << radius <<std::endl;
    tau=1.0;

    Vector<Complex> y_hat(nrMeasuredData);
    Vector<Complex> PHI_p(nrMeasuredData);		// /Phi(p^k)
    Vector<Complex> res(2*nrMeasuredData);
    std::cout<<"Number of measure points: " << nrMeasuredData << " with DataError: " << delta <<  std::endl;

    freqs.Part(0,nrMeasuredData);
    real.Part(0,nrMeasuredData);
    imag.Part(0,nrMeasuredData);
    amplitude_phase.Part(0,nrMeasuredData);


    // Settings for Newton-CG - routine
    int nrIterations=0, nrCGIterations=0;
    Double tol=0.001, tolCG=0.001, err=1.0, res_norm;

    //Settings for harmonic PDE - Driver
    Integer level=0;
    Boolean reset = TRUE;
    Integer pdenumber  = 0;
    if (!isPartOfSequence_)
      ptdomain_->PrintGrid(level);

    if (PrintGridOnly)
      exit(0);

    std::cout<<"Here begins communication with base PDE" << std::endl;

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

 
    std::cout<<" WriteGeneralPDEdefines() ... "<<std::endl;
    pdes_[0]->WriteGeneralPDEdefines(); 
 
    MaterialData * ptMaterial=pdes_[0]->getPDEMaterialData();   // Pointer to MaterialData
    Matrix<Double> * matMatrix =  ptMaterial->GetMatrix();
 
    ptBCs = pdes_[0]->getPDE_BCs();                             // Pointer to BCs
 
     ptAlgsys = pdes_[0]->getPDE_algsys();                       //Pointer to AlgebraicSystem
   
 
    Integer numElems = pdes_[0]->getPDE_numElems();
 
    dofs=pdes_[0]->getPDE_dofspernode();
   
    numNodes= pdes_[0]->getPDE_numPDENodes();
 
   std::cout<<"TEStMark:4"<<std::endl;
   
    while (nrIterations<2){
      nrIterations++;
      std::cout<<"\n NewtonCG ... Newton-Iteration-Nr = "<<nrIterations<<std::endl;
      calcNorm2Resid(res,anorm,nrMeasuredData);
      std::cout<<"A-Norm:" << anorm << std::endl;

      createF(ptMaterial, ptBCs, F_hat);
       createJacobiMatrix(ptMaterial, ptBCs, F_hat, parameterIncrement,JacobiMatrix, solElecPot, solMechDispl);
       createAdjointJacobiMatrix(parameterIncrement, parameter, JacobiMatrix,solElecPot,solMechDispl,freqs,adjJacobiMatrix);
      updateMaterialData(parameter, ptMaterial);         //member function of piezoParamIdent

      while(nrCGIterations<5&&res_norm>tolCG){ // CG
	nrCGIterations++;

      } // end while - CG

    } // end while-Newton

  }// End solveProblem

  void piezoParamIdent::calcNorm2Resid(Vector<Complex> &res, Double & anorm, Integer nrMeasuredData){
    anorm=0.0;
    for (int i=0;i<2*nrMeasuredData;i++){
      anorm+=abs(res[i])*abs(res[i]);
      anorm=sqrt(anorm);
    }
  }

  void piezoParamIdent::createF(MaterialData * ptMaterial, BCs * ptBCs, Vector<Complex> & F_hat){
    ENTER_FCN("PiezoParamIdent:createF");
    std::cout<<"creates F ..."<<std::endl;

    ptMaterial=pdes_[0]->getPDEMaterialData();   // Pointer to MaterialData
    Matrix<Double> * matMatrix =  ptMaterial->GetMatrix();
    ptBCs = pdes_[0]->getPDE_BCs();                             // Pointer to BCs
    ptAlgsys = pdes_[0]->getPDE_algsys();

    Integer level=0;
    Boolean reset = TRUE;
    Integer pdenumber = 0;

    Integer numElems = pdes_[0]->getPDE_numElems();

    Integer dofs=pdes_[0]->getPDE_dofspernode();
    std::cout<<"Create F, 2"<<std::endl;

    Integer numNodes= pdes_[0]->getPDE_numPDENodes();
    std::cout<<"Create F, 3"<<std::endl;

    Integer spaceDim = pdes_[0]->getPDE_spaceDim();
    std::cout<<"Create F, 4"<<std::endl;

    tau=1.0;
    Matrix<Double> couplingMatrix; // \bf e
    Matrix<Double> dielectricMatrix; // \eps^S
    couplingMatrix.Resize(spaceDim, 2*spaceDim);
    dielectricMatrix.Resize(spaceDim, spaceDim);
    Integer numElems_ = pdes_[0]->getPDE_numElems();


    createMaterialTensorMatrices(parameter, couplingMatrix, dielectricMatrix, spaceDim);
    // std::cout<<"Create F,3"<<std::endl;


    for (Integer fstep = 0; fstep < nrMeasuredData; fstep++) { // harmonic solver for different frequency - values

      Info->WriteHarmonicStep(pdes_[0]->GetName(), fstep, freqs[fstep]);


      Info->WriteHarmonicStep(pdes_[0]->GetName(), fstep, freqs[fstep]);
      std::cout<<"Create F, 5"<<std::endl;

      pdes_[0]->SolveStepHarmonic(fstep, freqs[fstep], level, reset);

      BaseNodeStoreSol * ptSol = pdes_[0]->getPDESolution();

      NodeStoreSol<Complex> * ptNodeStoreSol;
      ptNodeStoreSol = dynamic_cast<NodeStoreSol<Complex>*>(ptSol);
      std::cout<<"Create F, 7"<<std::endl;

      ptNodeStoreSol->GetGlobalSolVector(ELEC_POTENTIAL, solElecPot);
      ptNodeStoreSol->GetGlobalSolVector(MECH_DISPLACEMENT, solMechDispl);
      // Now we have \hat{d_k} and {\phi_k} for each omega ...
      //ptElem->GetGlobalDerivShFnc();

      // typeOutSolutionOnScreen(solElecPot, solMechDispl);              //member function of piezoParamIdent
      measureMechDeformationInZ_Direction(solMechDispl,radius,meanValueMechDeformation,dofs); // Braucht üÜberarbeitung!!



      pdes_[0]->PostStepHarmonic(fstep, freqs[fstep], level, reset);
      std::cout<<"Create F, 9"<<std::endl;
      //      actPDE->PostProcess(level);
      pdes_[0]->PostProcess(level);
      // pdes_[0]->PostStepHarmonic(fstep, freqs[fstep], level, reset);

      std::cout<<"Create F, 10"<<std::endl;
      Info->PrintPiezoMat(*ptMaterial);

      pdes_[0]->WriteResultsInFile();     
      std::cout<<"Finished to create F ..."<<std::endl;

    } // end of loop over all frequencies
  } // end createF



  void piezoParamIdent::createJacobiMatrix(MaterialData * ptMaterial, BCs * ptBCs, Vector<Complex> & F_hat, Vector<Double> & parameterIncrement,Matrix<Complex> & JacobiMatrix, Vector<Complex> & solElecPot,Vector<Complex> & solMechDispl){
    ENTER_FCN("piezoParamIdent::createJacobiMatrix");
    std::cout<<"JacobiMatrix will be created"<<std::endl;
    Vector<Double> IncrementedRHSMatrix;
   
    //    pdes_[0]->WriteGeneralPDEdefines();

    Matrix<Double> * matMatrix =  ptMaterial->GetMatrix();
    ptBCs = pdes_[0]->getPDE_BCs();                             // Pointer to BCs
    ptAlgsys = pdes_[0]->getPDE_algsys();
    Integer level=0;
    Boolean reset = TRUE;
    //Integer nrMeasuredData = freqs.GetSize();
    Integer numElems_ = pdes_[0]->getPDE_numElems();
    JacobiMatrix.Resize(2*nrMeasuredData);
    Integer pdenumber = 0;
    

    //     for(int i=0;i<parameter.GetSize();i++){
    //       parameter[i]=1; // noch reelle Parameterwerte !!
    //       parameter[i-1]-=1;
    //for (int i=0;i<parameter.GetSize();i++)
    // parameter[i]+=parameterIncrement[i];  //parameterIncrement noch ohne Inhalt!!
    //   } // should include the following lines!!! Closes after building up entries of F'(p^k) ....

    // updateMaterialData(parameter, ptMaterial);         //member function of piezoParamIdent

    for (Integer fstep = 0; fstep < nrMeasuredData; fstep++) { // harmonic solver for different frequency - values

      pdes_[0]->CreateIncrementedRHSMatrix(IncrementedRHSMatrix, freqs[fstep], level);

      // Now we calculate the MatVecProduct in the Right Hand Site, only real parts!!:
      BaseNodeStoreSol * ptSol = pdes_[0]->getPDESolution();

      NodeStoreSol<Complex> * ptNodeStoreSol;
      ptNodeStoreSol = dynamic_cast<NodeStoreSol<Complex>*>(ptSol);

      Vector<Complex> algSysSolVector;
      algSysSolVector=ptNodeStoreSol->GetAlgSysVector();


      Complex sum=0.0;
      Vector<Complex> RHSsol(algSysSolVector.GetSize());
      Integer j=0;
      Integer k=0;
      for(int i=0;i<IncrementedRHSMatrix.GetSize()/2; i++){
	sum=sum+IncrementedRHSMatrix[i]*algSysSolVector[k];
	k++;
	if ((i%algSysSolVector.GetSize())==0){
	  RHSsol[j]=sum;
          sum=0.0;
          j++;
 	  k=0;
	}
      }
      updateRHS(RHSsol);	//member function of piezoParamIdent
      Info->WriteHarmonicStep(pdes_[0]->GetName(), fstep, freqs[fstep]);
      pdes_[0]->PreStepHarmonic(fstep, freqs[fstep], level, reset);
      pdes_[0]->SolveStepHarmonic(fstep, freqs[fstep], level, reset);
      // ptNodeStoreSol->GetGlobalSolVector(ELEC_POTENTIAL, solElecPot);
      //ptNodeStoreSol->GetGlobalSolVector(MECH_DISPLACEMENT, solMechDispl);
      

      /*  for(int i=0; i<algSysSolVector.GetSize();i++){
	  std::cout<<"Algsysvector("<<i<<") " << algSysSolVector[i].real() << "+" << algSysSolVector[i].imag()<< std::endl;
	  std::cout<<"RHS(" <<i<<") " << IncrementedRHSMatrix[i] << std::endl;
	  }*/

      //  pdes_[0]->PostStepHarmonic(fstep, freqs[fstep], level, reset);

      pdes_[0]->PostProcess(level);

      Info->PrintPiezoMat(*ptMaterial);
      pdes_[0]->WriteResultsInFile();
     
     
    }		// end for over all frequencies
  }            //end CreateJacobiMatrix


  void piezoParamIdent::createAdjointJacobiMatrix(Vector<Double> & parameterIncrement,Vector<Double> &  parameter, Matrix<Complex> & JacobiMatrix, Vector<Complex> & solElecPot,Vector<Complex> & solMechDispl, Vector<Double> & freqs, Matrix<Complex> & adjJacobiMatrix){
    ENTER_FCN("piezoParamIdent::createAdjointJacobiMatrix");
    std::cout<<"Jacobi Matrix will be created"<<std::endl;
    Integer size = JacobiMatrix.GetSizeRow();
    //std::<<"SIZE OF JACOBIMATRIX"<<size<<std::endl;
    adjJacobiMatrix.Resize(size,size);
    for (int i=0;i<JacobiMatrix.GetSizeRow();i++)
      for (int j=0;j<JacobiMatrix.GetSizeCol();j++)
	adjJacobiMatrix[j][i] = std::conj(JacobiMatrix[i][j]);
  } // end createAdjointJacobiMatrix


  void piezoParamIdent::updateRHS(Vector<Complex> & solElecPot, Vector<Complex> & solMechDispl, Double omega){
    ENTER_FCN (" piezoParamIdent::updateRHS");
    Vector<Double> new_RHS;
    new_RHS.Resize(numNodes * dofs);
    for(int i=0;i<numNodes;i++){
      new_RHS[i]=solElecPot[i].real();
      for(int j=0;j<dofs-1;j++){	
	new_RHS[(j+1)*numNodes+i]=solMechDispl[j*numNodes+i].real();
	// matVecRHS(IncrementedRHSMatrix,solMechDispl, solElecPot,newRHS);	
      }
    }

    ptAlgsys->UpdateRHS(1,new_RHS.GetPointer());
    std::cout<<"RHS was updated"<< std::endl;
  } // end update RHS

  void piezoParamIdent::updateRHS(Vector<Complex> & RHSsol){
    Vector<Double> temp;
    temp.Resize(RHSsol.GetSize());
    for(int i=0;i<RHSsol.GetSize();i++)
      temp[i]=RHSsol[i].real();
    ENTER_FCN("piezoParamIdent::updateRHS2");
    ptAlgsys->InitRHS();
    ptAlgsys->UpdateRHS(1,temp.GetPointer());
    std::cout<<"RHS was updated"<< std::endl;
  }



  void piezoParamIdent::measureMechDeformationInZ_Direction(Vector<Complex> & mechDisplacement, Double & Radius, Double meanValueMechDeformation, int dof){
    meanValueMechDeformation=0.0;
    for (int i=dof-2;i<mechDisplacement.GetSize();i+=(dof-1)) // only in Z-Direction
      meanValueMechDeformation+=abs(mechDisplacement[i]);
    meanValueMechDeformation=meanValueMechDeformation/(PI*Radius*Radius*mechDisplacement.GetSize()*(dof-1));
    std::cout <<"MeanValue of MechDeformation: " << meanValueMechDeformation<< std::endl;
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
	  std::cout<<"Nr_of_measured_Data in readMeasuredData = "<< nrMeasuredData<<std::endl;
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
    ptMaterial->SetMatrixData(0,0, parameter[0]);
    ptMaterial->SetMatrixData(2,2, parameter[1]);
    ptMaterial->SetMatrixData(0,1, parameter[2]);
    ptMaterial->SetMatrixData(1,0, parameter[2]);
    ptMaterial->SetMatrixData(0,2, parameter[3]);
    ptMaterial->SetMatrixData(2,0, parameter[3]);
    ptMaterial->SetMatrixData(1,2, parameter[3]);
    ptMaterial->SetMatrixData(1,2, parameter[3]);
    ptMaterial->SetMatrixData(3,3, parameter[4]);
    ptMaterial->SetMatrixData(4,4, parameter[4]);
    ptMaterial->SetMatrixData(5,5, 0.5*(parameter[0]-parameter[2]));
    ptMaterial->SetMatrixData(6,4, parameter[5]);
    ptMaterial->SetMatrixData(7,3, parameter[5]);
    ptMaterial->SetMatrixData(8,0, parameter[6]);
    ptMaterial->SetMatrixData(8,1, parameter[6]);
    ptMaterial->SetMatrixData(8,2, parameter[7]);
    ptMaterial->SetMatrixData(6,6, parameter[8]);
    ptMaterial->SetMatrixData(7,7, parameter[8]);
    ptMaterial->SetMatrixData(8,8, parameter[9]);
  } // end updateMaterialData


} // end namespace CoupledField

/*	for(int i=0;i<nrMeasuredData;i++){
	amplitude_phase[i]=Complex(real[i],imag[i]);
	std::cout <<"Phase&Amplitude: " <<  amplitude_phase[i].real() <<" + "<< amplitude_phase[i].imag() << "i" <<std::endl;
	amplitude_phase[i]=std::polar(real[i],imag[i]);
	std::cout <<"Phase&Amplitude (polar): " <<  amplitude_phase[i].real() <<" + "<< amplitude_phase[i].imag() << "i"    <<std::endl; }


	for (int i=0; i<nrMeasuredData;i++){
	y_hat[i]=std::polar(-real[i], imag[i]);
	y_hat[i]=voltage*y_hat[i];
	std::cout <<"Phase&Amplitude (polar): " <<  y_hat[i].real() <<" + "<< y_hat[i].imag() << "i" <<std::endl;
	}


	std::cout << " SizeOf amplitude_phase" <<amplitude_phase.GetSize()<< std::endl;
	std::cout << " SizeOf real" <<real.GetSize()<< "nrmeasuredData= " << nrMeasuredData << std::endl;
*/
