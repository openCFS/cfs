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

#include <stdlib.h>

#include <sstream>
#include <iomanip>

#include <stdio.h>
#include <cstdarg>
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

  void piezoParamIdent::calc_measuredCharge(Vector<Double> freqs, Vector<Double> & absZ, Vector<Double> & phi, Vector<Complex> & y_hat){
    Complex Z,j;
    Double x, y;
    j=Complex(0,1);
    Double pi = 3.14159265358979;
    for (int i=0; i<nrMeasuredData; i++){
      //      x=absZ[i]*cos(2*pi*phi[i]/360);
      //      y=absZ[i]*sin(2*pi*phi[i]/360);
      x=absZ[i]*cos(phi[i]);
      //      std::cout<<"\n cos (phi)=" << cos(phi[i]) << ", phi="<< phi[i];
      y=absZ[i]*sin(phi[i]);
      //      std::cout<<"\n sin (phi)=" << sin(phi[i]) << ", phi="<< phi[i];
      Z=Complex(x,y);
      y_hat[i]=voltage/(2.0*pi*Z*freqs[i]*j);
      std::cout<<i<<") = " << phi[i] << ", "<< freqs[i] <<", q= " << y_hat[i]<<", Z= " << Z << std::endl;
    }
  }// end calc_measuresCharge()


 void piezoParamIdent::CG(){
    ENTER_FCN("piezoParamIdent::CG");
    std::cout<<"HERE STARTS CG ..."<<std::endl;
    Complex Jfs;
    Integer nrCGIt=0;
    s.Resize(nrParameter);    //sets all s_i to zero

      // res = y_hat-F_hat-JacobiMatrix*s_0
      for (int i=0; i< res.GetSize();i++){
	for (int j=0;j<s_0.GetSize();j++)
	  Jfs+=JacobiMatrix[i][j]*s[j];
	res[i]=y_hat[i]-F_hat[i]-Jfs;
	overall_res0[i]=y_hat[i]-F_hat[i];
	//	std::cout<<" res["<<i<<"]="<<res[i];
	Jfs=Complex(0,0);
      }
      //   std::cout<<"res = y_hat-F_hat-JacobiMatrix*s_0"<<std::endl;

      // bas = res_NE = adjJacMat*res
      for (int i=0; i< bas.GetSize();i++){
	for (int j=0;j<res.GetSize();j++)
	  Jfs+=adjJacobiMatrix[i][j]*res[j];
	//		bas[i]=res_NE[i]=1/scaling[i]*Jfs;
		 bas[i]=res_NE[i]=Jfs;
		//	std::cout<<" bas["<<i<<"]="<<bas[i];
	Jfs=Complex(0,0);
      }
      //      std::cout<<"  bas = res_NE = adjJacMat*res"<<std::endl;
      //      std::cout<<"\n";
      //      eta_k=eta_k*0.9;
      Vector<Complex> res_temp(y_hat.GetSize());
      res_temp=y_hat-F_hat;

      norm_res=a2norm(res);
      overall_res=a2norm(overall_res0);
      std::cout<<"\nOVERALL_RES: " << overall_res << "; norm_res" << norm_res;

      while((nrCGIt<15&&norm_res<eta*overall_res)||nrCGIt<3){ // CG
	nrCGIt++;
	std::cout<<"\nCG-Iteration -Nr: "<< nrCGIt;
        //bas_bar=JacMatr*bas
        for (int i=0; i< bas_bar.GetSize();i++){
	  for (int j=0;j<bas.GetSize();j++)
	    Jfs+=JacobiMatrix[i][j]*bas[j];
	  bas_bar[i]=Jfs;
	  //	std::cout<<" bas_bar["<<i<<"]="<<bas_bar[i];
	  Jfs=Complex(0,0);
	  // std::cout<<"res= "<< res[i];
	}
	//	std::cout<<"bas_bar=JacMatr*bas"<<std::endl;

	//scaling
	for(int i=0; i<res_NE.GetSize(); i++){
	  res_NE[i]=res_NE[i]*scaling[i];
	  bas_bar[i]=bas_bar[i]*scaling[i];
	  //	  std::cout<<" bas_bar["<<i<<"]="<<bas_bar[i];
	}

	//		alpha_m=pow(a2norm(res_NE),2)/pow(a2norm(bas_bar),2);
	alpha_m=a2norm(res_NE)/a2norm(bas_bar);

	std::cout<<"\n alpha = " <<alpha_m<<std::endl;

	//      std::cout<<"\n";

	// s_m = s_{m-1} + \alpha*bas
	for(int i=0; i< s.GetSize();i++){
	  s[i]=s[i]+alpha_m*bas[i];
      	  std::cout<<"s("<<i<<")="<<s[i]<<"; ";
	}
	//res=res-alpha_m*bas_bar
	for(int i=0; i<res.GetSize();i++)
	  res[i]-=alpha_m*bas_bar[i];

       //res_NE = adjJacDet*res
       for (int i=0; i< res_NE.GetSize();i++){
	  for (int j=0;j<res.GetSize();j++)
	    Jfs+=adjJacobiMatrix[i][j]*res[j];
	   res_NE_new[i]=1/(scaling[i]*scaling[i])*Jfs;
	  Jfs=Complex(0,0);
       }

       // scaling again:
       for(int i=0;i<res_NE.GetSize();i++){
	  res_NE_new[i]=res_NE_new[i]*scaling[i];
	  res_NE[i]=res_NE[i]*scaling[i];
       }

       // beta_m=pow(a2norm(res_NE_new),2)/pow(a2norm(res_NE),2);
        beta_m=a2norm(res_NE_new)/a2norm(res_NE);
       res_NE=res_NE_new;
       //       std::cout<<"\n beta_m= "<< beta_m;

	       //bas=res_NE+beta_m*bas
      for (int i=0;i<bas.GetSize();i++)
	       bas[i]=res_NE[i]+beta_m*bas[i];

      for (int i=0; i< res.GetSize();i++){
	for (int j=0;j<s_0.GetSize();j++)
	  Jfs+=JacobiMatrix[i][j]*s[j];
	res[i]=y_hat[i]-F_hat[i]-Jfs;
      }
	norm_res=a2norm(res);

      } // end while - CG

  }// end CG


void piezoParamIdent::createF(MaterialData * ptMaterial, BCs * ptBCs, Vector<Complex> & F_hat){
    ENTER_FCN("PiezoParamIdent:createF");
    //    std::cout<<"creates F ..."<<std::endl;

    ptMaterial=pdes_[0]->getPDEMaterialData();   // Pointer to MaterialData
    Matrix<Double> * matMatrix =  ptMaterial->GetMatrix();
    ptBCs = pdes_[0]->getPDE_BCs();                             // Pointer to BCs
    ptAlgsys = pdes_[0]->getPDE_algsys();

    Integer level=0;
    Boolean reset = TRUE;
    Integer pdenumber = 0;

    Integer numElems = pdes_[0]->getPDE_numElems();

    Integer dofs=pdes_[0]->getPDE_dofspernode();  

    Integer numNodes= pdes_[0]->getPDE_numPDENodes();  

    //   Integer spaceDim = pdes_[0]->getPDE_spaceDim();  

    tau=1.0;
    // Matrix<Double> couplingMatrix; // \bf e
    //  Matrix<Double> dielectricMatrix; // \eps^S
    //  couplingMatrix.Resize(spaceDim, 2*spaceDim);
    //  dielectricMatrix.Resize(spaceDim, spaceDim);
    Integer numElems_ = pdes_[0]->getPDE_numElems();
    Grid * ptGrid =   ptdomain_->GetGrid();

    //completeSolOf_F.Resize(1,numNodes*dofs);

    //    createMaterialTensorMatrices(parameter, couplingMatrix, dielectricMatrix, spaceDim); 

    for (Integer fstep = 0; fstep < nrMeasuredData; fstep++) { // harmonic solver for different frequency - values
      reset = TRUE;

      Info->WriteHarmonicStep(pdes_[0]->GetName(), fstep, freqs[fstep]);
     
      pdes_[0]-> setBCs_id_phase_(0, imag[fstep]);
      

      // !!!! Boundary  Conditions are set, but do not go over into algebraic system!!

      // ptdomain_->Update(level);      
      // ptBCs->Update(ptGrid);
      // ptdomain_->UpdateAlgSys(level);
    
      pdes_[0]->WriteGeneralPDEdefines(); 
      pdes_[0]->PreStepHarmonic(fstep, freqs[fstep], level, reset);
      pdes_[0]->SolveStepHarmonic(fstep, freqs[fstep], level, reset);

      BaseNodeStoreSol * ptSol = pdes_[0]->getPDESolution();

      NodeStoreSol<Complex> * ptNodeStoreSol;
      ptNodeStoreSol = dynamic_cast<NodeStoreSol<Complex>*>(ptSol);

      ptNodeStoreSol->GetGlobalSolVector(ELEC_POTENTIAL, solElecPot);
      ptNodeStoreSol->GetGlobalSolVector(MECH_DISPLACEMENT, solMechDispl);
    
      //      typeOutSolutionOnScreen(solElecPot, solMechDispl);              //member function of piezoParamIdent
      measureMechDeformationInZ_Direction(solMechDispl,radius,meanValueMechDeformation,dofs); // Braucht üÜberarbeitung!!

      // std::cout<<"meanValueMechDef: "<< meanValueMechDeformation <<std::endl;
    
      F_hat[fstep + nrMeasuredData]=meanValueMechDeformation;

      pdes_[0]->PostStepHarmonic(fstep, freqs[fstep], level, reset);
      pdes_[0]->PostProcess(level);
      pdes_[0]->PostStepHarmonic(fstep, freqs[fstep], level, reset);

      Vector<Complex> chargeVec =   pdes_[0]->getPDE_complexValuedCharge();
      Complex mean_value_charge=Complex(0.0,0.0);
      for (int i=0;i<chargeVec.GetSize();i++)
	mean_value_charge+=chargeVec[i];
      mean_value_charge = (1/Double(chargeVec.GetSize())) * mean_value_charge;
      //      std::cout<<"Integral over Charge: "<< mean_value_charge << std::endl;

      F_hat[fstep]=mean_value_charge; 

      calcAbsImped(mean_value_charge, freqs[fstep], fstep);   // calculates |Z| and writes results in File
     
      Info->PrintPiezoMat(*ptMaterial);

      pdes_[0]->WriteResultsInFile();     
      
      //Gathers solution (d_k^l, \phi_k^l) to provide it for RHS of PDE which solution gives F'
      Vector<Complex> algSysSolVector;
      algSysSolVector=ptNodeStoreSol->GetAlgSysVector();
      if (fstep==0)
	completeSolOf_F.Resize(nrMeasuredData,algSysSolVector.GetSize());
      for (int i=0;i<algSysSolVector.GetSize();i++)
	completeSolOf_F[fstep][i]=algSysSolVector[i];

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

	      for (int k=0;k<elSolVec.GetSize();k++)
		allElemsVec[fstep][actEl*elSolVec.GetSize()+k] = elSolVec[k];
	  } // end loop over elements

    } // end of loop over all frequencies

    //    std::cout<<"Finished to create F ... here it is:"<<std::endl;

    // for (int i=0;i<F_hat.GetSize();i++)
    // std::cout<<"F("<<i+1<<")="<<F_hat[i]<<"; ";

  } // end createF



  void piezoParamIdent::createJacobiMatrix(MaterialData * ptMaterial, BCs * ptBCs, Vector<Complex> & F_hat, Vector<Double> & parameterIncrement, Matrix<Complex> & JacobiMatrix, Vector<Complex> & solElecPot,Vector<Complex> & solMechDispl){
    ENTER_FCN("piezoParamIdent::createJacobiMatrix");
    std::cout<<"JacobiMatrix will be created"<<std::endl;
    Vector<Double> IncrementedRHSMatrix;
   
    //    pdes_[0]->WriteGeneralPDEdefines();

    Matrix<Double> * matMatrix =  ptMaterial->GetMatrix();
    ptBCs = pdes_[0]->getPDE_BCs();                             // Pointer to BCs
    ptAlgsys = pdes_[0]->getPDE_algsys();
    Integer level=0;
    Boolean reset = TRUE;
    ptGrid = pdes_[0]->getPDE_grid();
    ptNodeEqn = pdes_[0]->getPDE_eqnData();
    //Integer nrMeasuredData = freqs.GetSize();
    Integer numElems_ = pdes_[0]->getPDE_numElems();
    nrParameter=parameter.GetSize();
    JacobiMatrix.Resize(2*nrMeasuredData,nrParameter);
    //    std::
    Integer pdenumber = 0;
    
    for(int ind_param=0; ind_param<nrParameter;ind_param++){

      //         parameter[ind_param]+=std::abs(parameter[ind_param]);
      // if(ind_param>0)
      //	parameter[ind_param-1]-=0.5*std::abs(parameter[ind_param-1]);
       
       parameter[ind_param]+= 1/(scaling[ind_param]);
       if (ind_param>0)
       	parameter[ind_param-1]-=1/(scaling[ind_param-1]);
       //       std::cout<<"Value of parameter["<< ind_param <<"] is now = "<<parameter[ind_param]<<std::endl;
       //std::cout<<"Value of parameter["<< ind_param-1 <<"] was reset to = "<<parameter[ind_param-1]<<std::endl;

       updateMaterialData(parameter, ptMaterial);         //member function of piezoParamIdent

       if (ind_param==nrMeasuredData-1)
       	parameter[ind_param]-=10/scaling[ind_param];


       //      pdes_[0]->DefineIntegrators(0);

      for (Integer fstep = 0; fstep < nrMeasuredData; fstep++) { // harmonic solver for different frequency - values
	reset = TRUE;

	//	Matrix<Complex> sysmat = ptAlgsys->GetSysMat();

	// try to generate Sysmat which is used in RHS ...
       
        StdVector<Elem*> elemssd;
	subdoms = pdes_[0]->getPDE_subdoms();
	//	std::cout<<"PDE_SUBDOM[0] = "<< subdoms[0]; 
	 ptGrid->GetElemSD(elemssd,subdoms[0], level);

	BaseNodeStoreSol * ptSol = pdes_[0]->getPDESolution();
	NodeStoreSol<Complex> * ptNodeStoreSol;
	ptNodeStoreSol = dynamic_cast<NodeStoreSol<Complex>*>(ptSol);

	Vector<Complex> algSysSolVector;
	algSysSolVector=ptNodeStoreSol->GetAlgSysVector();
	std::cout<<"\n algSysSolVec : "<< algSysSolVector.GetSize()<<std::endl;
	Vector<Complex> RHSVec(algSysSolVector.GetSize());
	
	
	for (int actEl=0; actEl< elemssd.GetSize(); actEl++) {
	    BaseFE * ptEl = elemssd[actEl]->ptElem;
	    StdVector<Integer> connecth = elemssd[actEl]->connect;
		             
	    Matrix<Double> ptCoord;
	    ptGrid->GetCoordNodesElemMat(connecth, ptCoord, 0);
		    
	    // map connect to PDE node numbers
	    StdVector<Integer> connect_PDE;
	    
	   	    
	    ptNodeEqn->Node2EQN(connecth, connect_PDE);

	    std::cout<<"CONNECTh"<<std::endl;
	     for(int i=0;i<connecth.GetSize();i++)
	      std::cout<< connecth[i] << "; "<<std::endl;

	     std::cout<<"CONNECT_ PDE"<<std::endl;
	     for(int i=0;i<connect_PDE.GetSize();i++)
	      std::cout<< connect_PDE[i] << "; "<<std::endl;


	     Matrix<Complex> elSolMat;
	     //  ptSol = pdes_[0]->getPDESolution();

	     NodeStoreSol<Complex> * ptNodeStoreSol;
	     ptNodeStoreSol = dynamic_cast<NodeStoreSol<Complex>*>(ptSol);
	    
	     ptSol->GetElemSolutionAsMatrix(elSolMat, connecth);
	     //	      std::cout<<"\n\n - ELSOLMAT - :"<<std::endl;	    

	      Vector<Complex> elSolVec; 
              ptNodeStoreSol->GetElemSolution(elSolVec,connecth);

	      //	      std::cout<<"\n\n - ELSOLVEC - :"<<std::endl;
	      //for (int i=0;i<elSolVec.GetSize();i++)
	      //	std::cout<<elSolVec[i]<<"; ";

	      //std::cout.flush();
	     //	     BaseForm * bilinearStiff = GetStiffIntegrator(0);
	     //    else if (subType_ == "3d")

	     MaterialData actSDMat(*ptMaterial);
	     Boolean isdamping=FALSE;

	     //             BaseForm * bilinearStiff = GetStiffIntegrator(actSDMat);

	     BaseForm * bilinearStiff;  
	     bilinearStiff = new linPiezo3DInt(actSDMat,isdamping);

	      IntegratorDescriptor *actIntDescrStiff = new IntegratorDescriptor(bilinearStiff, STIFFNESS);
	      // ptAssemble->AddIntegrator(actIntDescrStiff, subdoms[0]);

	      //bilinearStiff->SetActElemSol(elSolMat);
	      Matrix<Double> elemmat;
	      bilinearStiff->SetElemPtr(ptEl);
	      bilinearStiff->CalcElementMatrix(ptCoord, elemmat);
	      std::cout<<" \n \n ELEMMAT: "<<elemmat.GetSizeRow() << " x " <<elemmat.GetSizeCol() << std::endl;     
	      //	     for (int i=0;i<elemmat.GetSizeRow();i++)
	      // for(int j=0;j<elemmat.GetSizeCol();j++){
	      //	 std::cout<<elemmat[i][j]<<"; ";
	      //	 if (j==elemmat.GetSizeRow()-1)
	      //	   std::cout<<"\n";
	      // }
	
	      //temp = elemmat*elemvec;
	      Complex sum=0.0;
	      Vector<Complex> temp;
	      temp.Resize(elemmat.GetSizeRow());

	      for (int i=0; i<elemmat.GetSizeRow();i++){
		for (int j=0; j<elemmat.GetSizeCol();j++)
		  sum=sum+elemmat[i][j]*allElemsVec[fstep][fstep*elemmat.GetSizeRow()+j];  
		temp[i]=sum;
		std::cout<< temp [i]<< "; ";
		sum=Complex(0,0);
	      }
	      std::cout<<"First for for runs well ... "<<std::endl;
	      Integer iEQN, iEQNDof;
	      std::cout<<" Size of connect_PDE: " << connect_PDE.GetSize()<< std::endl;
	      std::cout<<" Size of connecth: "<< connecth.GetSize()<< std::endl;
	      std::cout<<" Size of dofs: " << dofs << std::endl;
	     
	      for(int iNode=0;iNode<connecth.GetSize()*dofs;iNode++)
	      	for(int iDof=0;iDof<dofs;iDof++){
	      	  ptNodeEqn->Node2EQN(iNode,iDof,iEQN,iEQNDof);
		  //				  RHSVec[iEQN*dofs+iEQNDof]=temp[iNode*dofs+iDof];
	      	}

	 } // end for elemssd
		
 
        pdes_[0]-> setBCs_id_phase_(0, imag[fstep]);


	//	pdes_[0]->CreateIncrementedRHSMatrix(IncrementedRHSMatrix, freqs[fstep], level);
	// Now we calculate the MatVecProduct in the Right Hand Site, only real parts!!:

	//	Vector<Complex> algSysSolVector;
	algSysSolVector=ptNodeStoreSol->GetAlgSysVector();
	//	std::cout<<"\n algSysSolVec : "<< algSysSolVector.GetSize()<<std::endl;

	Complex sum=0.0;
	Vector<Complex> RHSsol(algSysSolVector.GetSize());
	Integer j=0;
	Integer k=0;
	for(int i=0;i<IncrementedRHSMatrix.GetSize()/2; i++){
	  //	  sum=sum+IncrementedRHSMatrix[i]*algSysSolVector[k];
	  sum=sum+IncrementedRHSMatrix[i]*completeSolOf_F[fstep][k];
	  //	  std::cout<<IncrementedRHSMatrix[i]<<"; ";
	  k++;
	  if ((i%algSysSolVector.GetSize())==0){
	    RHSsol[j]=sum;
	    sum=0.0;
	    j++;
	    k=0;
	  }
	  //	  std::cout<<RHSsol[j]<<"; ";
	}
	//	std::cout<<"\n --- \n";
       	updateRHS(RHSsol);	//member function of piezoParamIdent
	Info->WriteHarmonicStep(pdes_[0]->GetName(), fstep, freqs[fstep]);
 	pdes_[0]->WriteGeneralPDEdefines();
 	pdes_[0]->PreStepHarmonic(fstep, freqs[fstep], level, reset);
 	pdes_[0]->SolveStepHarmonic(fstep, freqs[fstep], level, reset); 

	//      ptNodeStoreSol->GetGlobalSolVector(ELEC_POTENTIAL, solElecPot);
	ptNodeStoreSol->GetGlobalSolVector(MECH_DISPLACEMENT, solMechDispl);
    
	// typeOutSolutionOnScreen(solElecPot, solMechDispl);              //member function of piezoParamIdent
	measureMechDeformationInZ_Direction(solMechDispl,radius,meanValueMechDeformation,dofs); // Braucht üÜberarbeitung!!
    
	JacobiMatrix[fstep + nrMeasuredData][ind_param]=meanValueMechDeformation;    

        pdes_[0]->PostStepHarmonic(fstep, freqs[fstep], level, reset);
	pdes_[0]->PostProcess(level);
        pdes_[0]->PostStepHarmonic(fstep, freqs[fstep], level, reset);
	Vector<Complex> chargeVec =  pdes_[0]->getPDE_complexValuedCharge();
	Complex<Double> mean_value_charge;

	for (int i=0; i<chargeVec.GetSize();i++)
	  mean_value_charge=mean_value_charge+chargeVec[i];
	mean_value_charge = mean_value_charge/(Double(chargeVec.GetSize()));

	JacobiMatrix[fstep][ind_param]=mean_value_charge; //pdes_[0]->getPDE_complexValuedCharge();    

	//   Info->PrintPiezoMat(*ptMaterial);
	//     pdes_[0]->WriteResultsInFile();     
     
      }		// end for over all frequencies

    } //end loop over paramters

    // for(int i=0;i<JacobiMatrix.GetSizeRow();i++)
      // for (int j=0; j<JacobiMatrix.GetSizeCol();j++){
	//	std::cout<<"F'("<<i<<")("<<j<<")= "<< JacobiMatrix[i][j]<<"; ";
	//	if (j==JacobiMatrix.GetSizeCol()-1)
	//  std::cout<<"\n";
    //      }

  }            //end CreateJacobiMatrix


  void piezoParamIdent::createAdjointJacobiMatrix(Vector<Double> & parameterIncrement,Vector<Double> &  parameter, Matrix<Complex> & JacobiMatrix, Vector<Complex> & solElecPot,Vector<Complex> & solMechDispl, Vector<Double> & freqs, Matrix<Complex> & adjJacobiMatrix){
    ENTER_FCN("piezoParamIdent::createAdjointJacobiMatrix");
    std::cout<<"Adjoint Jacobian will be created ... "<<std::endl;
    adjJacobiMatrix.Resize(JacobiMatrix.GetSizeRow(),JacobiMatrix.GetSizeCol());
    for (int i=0;i<JacobiMatrix.GetSizeRow();i++)
      for (int j=0;j<JacobiMatrix.GetSizeCol();j++)
	adjJacobiMatrix[j][i] = std::conj(JacobiMatrix[i][j]);
  } // end createAdjointJacobiMatrix



} // end namespace

