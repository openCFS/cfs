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

  void piezoParamIdent::calc_measuredCharge(Vector<Double> freqs, Vector<Double> & absZ, Vector<Double> & phi, Vector<Complex> & y_hat){
    Complex Z,j;
    Double x, y;
    j=Complex(0,1);
    for (int i=0; i<nrMeasuredData; i++){
      x=absZ[i]*cos(phi[i]);
      y=absZ[i]*sin(phi[i]);
      Z=Complex(x,y);
      y_hat[i]=voltage/(2.0*PI*Z*freqs[i]*j);
      //      y_hat[i]=voltage/(Z*freqs[i]*j);
      std::cout<<i<<") = " << phi[i] << ",\t "<< freqs[i] <<",\t q = y_hat = " << y_hat[i]<<",\t Z= " << Z << std::endl;
    }
  }// end calc_measuresCharge()


  void piezoParamIdent::CG(){
    ENTER_FCN("piezoParamIdent::CG");
    std::cout<<"\nHERE STARTS CG ..."<<std::endl;
    Complex Jfs=Complex(0,0);
    Integer nrCGIt=0;
    s.Resize(nrParameter);    //sets all s_i to zero
    bas.Resize(nrParameter); 

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
      //		bas[i]=res_NE[i]=1.0/scaling[i]*Jfs;
      //bas[i]=res_NE[i]=1.0/(scaling[i]*scaling[i])*Jfs;

      bas[i]=res_NE[i]=Jfs/(scaling[i]*scaling[i]);
      //	std::cout<<" bas["<<i<<"]="<<bas[i];
      Jfs=Complex(0,0);
    }

    Vector<Complex> res_temp(y_hat.GetSize());
    overall_res0=y_hat-F_hat;
    norm_res=a2norm(res);
    overall_res=a2norm(overall_res0);
    //    std::cout<<"\nOVERALL_RES: " << overall_res << "; norm_res" << norm_res;

    while((nrCGIt<10&&norm_res<=eta*eta*overall_res)||nrCGIt<2){ // CG
      nrCGIt++;
      std::cout<<"\nCG-Iteration -Nr: "<< nrCGIt << std::endl;
      std::cout<< "\n norm_res = " << norm_res << "; overall_res = "<< overall_res << "; eta^2*overall-norm= "<< eta*eta*overall_res-norm_res << std::endl;

       //bas_bar=JacMatr*bas
      JacobiMatrix.Mult(bas,bas_bar);
//       for (int i=0; i< bas_bar.GetSize();i++){
// 	for (int j=0;j<bas.GetSize();j++)
// 	  Jfs+=JacobiMatrix[i][j]*bas[j];
// 	bas_bar[i]=Jfs;
// 	Jfs=Complex(0,0);
// 	// std::cout<<"res= "<< res[i];
//       }
      //	std::cout<<"bas_bar=JacMatr*bas"<<std::endl;

      //scaling
      for(int i=0; i<res_NE.GetSize(); i++){
     	res_NE[i]=res_NE[i]*scaling[i];
     	bas_bar[i]=bas_bar[i]*scaling[i];
	//	bas[i]=bas[i]*scaling[i];
	//	  std::cout<<" bas_bar["<<i<<"]="<<bas_bar[i];
      }

      //alpha_m=POW(a2norm(res_NE),2)/POW(a2norm(bas_bar),2);
      alpha_m=a2norm(res_NE)/a2norm(bas_bar);

      std::cout<<"\n alpha = " <<alpha_m<<std::endl;

      // s_m = s_{m-1} + \alpha*bas
      for(int i=0; i< s.GetSize();i++){
	s[i]=s[i]+alpha_m*bas[i];
	std::cout<<"s("<<i<<")="<<s[i]<<"; ";
      }
      //res=res-alpha_m*bas_bar
      for(int i=0; i<res.GetSize();i++)
	res[i]=res[i]-alpha_m*bas_bar[i];
    
      //res_NE = adjJacDet*res
      for (int i=0; i< res_NE.GetSize();i++){
	for (int j=0;j<res.GetSize();j++)
	  Jfs+=adjJacobiMatrix[i][j]*res[j];
	res_NE_new[i]=1.0/(scaling[i]*scaling[i])*Jfs;
       	//res_NE_new[i]=1.0/(scaling[i])*Jfs;
	res_NE_new[i]=Jfs;
	Jfs=Complex(0,0);
      }

      // scaling again:
           for(int i=0;i<res_NE.GetSize();i++){
	     //	res_NE_new[i]=res_NE_new[i]*scaling[i];
	     res_NE[i]=res_NE[i]*scaling[i];
       }

      // beta_m=POW(a2norm(res_NE_new),2)/POW(a2norm(res_NE),2);
      beta_m=a2norm(res_NE_new)/a2norm(res_NE);
      res_NE=res_NE_new;
      std::cout<<"\n beta_m= "<< beta_m;

      //bas=res_NE+beta_m*bas
      for (int i=0;i<bas.GetSize();i++)
	bas[i]=res_NE[i]+beta_m*bas[i];

      for (int i=0; i< res.GetSize();i++){
	for (int j=0;j<s.GetSize();j++)
	  Jfs+=JacobiMatrix[i][j]*s[j];
	lin_res[i]=y_hat[i]-F_hat[i]-Jfs;
	Jfs=Complex(0,0);
      }
      norm_res=a2norm(lin_res);
     
    } // end while - CG

  }// end CG

  void piezoParamIdent::backtracking(Double & eta, Double & theta, Vector<Complex> & s, Double & norm_res, Double & norm_res_new){
   ENTER_FCN("piezoParamIdent::backtracking");



  Double real_misfit, real_misfit_new,eta_max, eta_new, t, a,b,c, theta_min, theta_max, gamma, al;

  real_misfit_new=norm_res_new;
  real_misfit = norm_res;

      t=1.0e-4;
      eta_max=0.9;
      eta=0.5;
      theta_min=0.1;
      theta_max=0.5;
      gamma=0.5;
      al=1.5;
      tau=0.3;
      Vector<Complex> y_hat_F_hat (2*nrMeasuredData);

      int backtrackIterator=0;
      Double temp = (1-t*(1-eta))*real_misfit;
      std::cout<<"\n backtracking: real(misfit) = " <<real_misfit_new<< "; 1-t(1-eta)*real_misfit= "<< temp<<std::endl;
  
      while((real_misfit_new>(1-t*(1-eta))*real_misfit)&&(eta<eta_max)){
	std::cout<<"backtracking ... "<< backtrackIterator<< std::endl;
	backtrackIterator++;

	// choose theta:
	a=real_misfit*real_misfit;
	b=2*real_misfit;
	c=real_misfit*real_misfit-b-a;
	std::cout<<"c"<<c<<std::endl;
	if (c==0.0){
	  if(b<=0.0)
	    theta=theta_max;
	  else 
	    theta=theta_min;
	}
	else if (c>0.0){
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

	std::cout<<"\n Choice of theta = " << theta<<std::endl;
	//end choose theta

	for (int i=0;i<s.GetSize();i++){
	  s[i]=theta*s[i];
	  std::cout<<"; s= "<< s[i];
	}
	eta=1-theta*(1-eta);
	std::cout<<"\n ETA = " << eta << std::endl;

	// we do a bit of scaling here ...
	for(int i=0;i<scaling.GetSize();i++)
	  parameter_new[i]=parameter_new[i]*scaling[i];

	for (int i=0;i<parameter.GetSize();i++)
	  parameter_new[i]=parameter[i]+s[i].real();

	// we do a bit of rescaling here ...
	for(int i=0;i<scaling.GetSize();i++)
	  parameter_new[i]=parameter_new[i]*1/scaling[i];

	updateMaterialData(parameter_new, ptMaterial); 
	createF(ptMaterial, ptBCs, F_hat);       

	// we do a bit of scaling here ...
	//	for(int i=0;i<scaling.GetSize();i++)
	// parameter_new[i]=parameter_new[i]*scaling[i];

	y_hat_F_hat=y_hat-F_hat;
	real_misfit_new=a2norm(y_hat_F_hat);

       	std::cout<<"\n real_misfit_new = "<<real_misfit_new;
      } // end Liniensuche


	//choose eta
	eta_new=gamma*POW((real_misfit_new/real_misfit),al);
	if (gamma*POW(eta,al) > 0.1)
	  eta_new=std::max(eta_new,gamma*POW(eta,al));
	//       std::cout<< "\n eta_new " << eta_new <<std::endl;
	if (eta_new>eta_max)
	  eta_new=eta_max;
	// std::cout<< "\n eta_new " << eta_new <<std::endl;
	if(eta_new<=2*(tau*delta)/real_misfit_new)
	  eta_new=0.8*(tau*delta)/real_misfit_new;
	//end choose eta
	std::cout<< "\nChoice of eta ..." << eta_new<<std::endl; 
	eta=eta_new;
	real_misfit=real_misfit_new;
	parameter=parameter_new;      


  }// end backtracking

  void piezoParamIdent::createF(MaterialData * ptMaterial, BCs * ptBCs, Vector<Complex> & F_hat){
    ENTER_FCN("PiezoParamIdent:createF");
    //    std::cout<<"creates F ..."<<std::endl;

    ptMaterial=pdes_[0]->getPDEMaterialData();   // Pointer to MaterialData
    //    Matrix<Double> * matMatrix =  ptMaterial->GetMatrix();
    ptBCs = pdes_[0]->getPDE_BCs();                             // Pointer to BCs
    ptAlgsys = pdes_[0]->getPDE_algsys();
    Grid * ptGrid =   ptdomain_->GetGrid();

    Integer level=0;
    Boolean reset = TRUE;
    Integer pdenumber = 0;

    Integer numElems = pdes_[0]->getPDE_numElems();
    Integer dofs=pdes_[0]->getPDE_dofspernode();  
    Integer numNodes= pdes_[0]->getPDE_numPDENodes();  

    tau=1.0;
    // Matrix<Double> couplingMatrix; // \bf e
    //  Matrix<Double> dielectricMatrix; // \eps^S
    //  couplingMatrix.Resize(spaceDim, 2*spaceDim);
    //  dielectricMatrix.Resize(spaceDim, spaceDim);
    //    createMaterialTensorMatrices(parameter, couplingMatrix, dielectricMatrix, spaceDim); 

       updateMaterialData(parameter, ptMaterial);         //member function of piezoParamIdent

    for (Integer fstep = 0; fstep < nrMeasuredData; fstep++) { // harmonic solver for different frequency - values
      reset = TRUE;

      //      Info->WriteHarmonicStep(pdes_[0]->GetName(), fstep, freqs[fstep]);
     
            pdes_[0]-> setBCs_id_phase_(0, imag[fstep]);      

      // ptdomain_->Update(level);      
      // ptBCs->Update(ptGrid);
      // ptdomain_->UpdateAlgSys(level);

	////////////////////////////////////////////////////////
	//                   SOLVES PDE                      //
	///////////////////////////////////////////////////////  

      pdes_[0]->WriteGeneralPDEdefines();       
      pdes_[0]->PreStepHarmonic(fstep, freqs[fstep], level, reset);      
      pdes_[0]->SolveStepHarmonic(fstep,freqs[fstep], level, reset);
      pdes_[0]-> setBCs_id_phase_(0, imag[fstep]);      

      pdes_[0]->PostStepHarmonic(fstep, freqs[fstep], level, reset);
      pdes_[0]->PostProcess(level);
      pdes_[0]->PostStepHarmonic(fstep, freqs[fstep], level, reset);

        //////////////////////////////////////////////////////////
	//Retrieves & stores Solution for further calculations  //
	/////////////////////////////////////////////////////////

      BaseNodeStoreSol * ptSol = pdes_[0]->getPDESolution();
      NodeStoreSol<Complex> * ptNodeStoreSol;
      ptNodeStoreSol = dynamic_cast<NodeStoreSol<Complex>*>(ptSol);

      if (considerMechDeformation==TRUE){
      // ptNodeStoreSol->GetGlobalSolVector(ELEC_POTENTIAL, solElecPot);
      ptNodeStoreSol->GetGlobalSolVector(MECH_DISPLACEMENT, solMechDispl);
      //typeOutSolutionOnScreen(solElecPot, solMechDispl);              //member function of piezoParamIdent
      meanValueMechDeformation=0.0;
      measureMechDeformationInZ_Direction(solMechDispl,radius,meanValueMechDeformation,dofs); // Braucht üÜberarbeitung!!
      // std::cout<<"meanValueMechDef: "<< meanValueMechDeformation <<std::endl; 
      F_hat[fstep + nrMeasuredData]=meanValueMechDeformation;
      }

      Vector<Complex> chargeVec =   pdes_[0]->getPDE_complexValuedCharge();

      Complex mean_value_charge=Complex(0.0,0.0);

      for (int i=0;i<chargeVec.GetSize();i++){
	//	std::cout<<"\n|charge("<<i<<")|= "<<std::abs(chargeVec[i])<<";\t charge("<<i<<")= "<<(chargeVec[i]);
	mean_value_charge+=chargeVec[i];
      }

      //        mean_value_charge = (1.0/Double(chargeVec.GetSize())) * mean_value_charge;
        std::cout<<"Integral over Charge: "<< mean_value_charge << std::endl;

      F_hat[fstep]=mean_value_charge; 

      // calcAbsImped(mean_value_charge, freqs[fstep], fstep);   // calculates |Z| and writes results in File
     
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

     std::cout<<"\n Finished to create F ... here it is:"<<std::endl;
         for (int i=0;i<F_hat.GetSize();i++)
    std::cout<<"F("<<i<<")="<<F_hat[i]<<"; ";

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

    JacobiMatrix.Resize(2*nrMeasuredData,nrParameter);
    if (considerMechDeformation==FALSE)
      JacobiMatrix.Resize(nrMeasuredData,nrParameter);

    Integer pdenumber = 0;

    for(int ind_param=0; ind_param<nrParameter;ind_param++){ // loop over different parameter increments
    // std::cout<<"\n"<<std::endl;
      //      for(Integer i=0;i<parameter.GetSize();i++)
      //	std::cout<<parameter[i]<<"; ";
      //        parameter[ind_param]+=std::abs(parameter[ind_param]);
      // if(ind_param>0)
      //	parameter[ind_param-1]-=0.5*std::abs(parameter[ind_param-1]);
       
      parameter[ind_param]+= 1.0/(scaling[ind_param]);       // we are incrementing one parameter after another
       if (ind_param>0)
           	parameter[ind_param-1]-=1.0/(scaling[ind_param-1]);
       updateMaterialData(parameter, ptMaterial);         //member function of piezoParamIdent
      
 if (ind_param==nrParameter-1)
       	 parameter[ind_param]-=1.0/scaling[ind_param];
       
 //   std::cout<<"\n"<<std::endl;
 //   for(Integer i=0;i<parameter.GetSize();i++)
 //   	std::cout<<parameter[i]<<"; ";

      //      pdes_[0]->DefineIntegrators(0);

      for (Integer fstep = 0; fstep < nrMeasuredData; fstep++) { // harmonic solver for different frequency - values
	reset = TRUE;

        pdes_[0]-> setBCs_id_phase_(0, imag[fstep]);

	//	Info->WriteHarmonicStep(pdes_[0]->GetName(), fstep, freqs[fstep]);
	
 	pdes_[0]->WriteGeneralPDEdefines();
 	pdes_[0]->PreStepHarmonic(fstep, freqs[fstep], level, reset);
        
	//	Cannot use SolveStepHarmonic, since it overwrites RHS ...
	//      for this, I have copied the method StepHarmonicLin to this place 
	//      void BasePDE::StepHarmonicLin(const Integer freqStep, const Double frequency, const Integer level, const Boolean reset)
	
	Integer job;

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
	  //      ptNodeStoreSol->GetGlobalSolVector(ELEC_POTENTIAL, solElecPot);
	  ptNodeStoreSol->GetGlobalSolVector(MECH_DISPLACEMENT, solMechDispl);
	  // typeOutSolutionOnScreen(solElecPot, solMechDispl);              //member function of piezoParamIdent
	  measureMechDeformationInZ_Direction(solMechDispl,radius,meanValueMechDeformation,dofs); // Braucht üÜberarbeitung!!
	  JacobiMatrix[fstep + nrMeasuredData][ind_param]=meanValueMechDeformation;    
	}

        pdes_[0]->PostStepHarmonic(fstep, freqs[fstep], level, reset);
	pdes_[0]->PostProcess(level);
        pdes_[0]->PostStepHarmonic(fstep, freqs[fstep], level, reset);
	Vector<Complex> chargeVec =  pdes_[0]->getPDE_complexValuedCharge();
	Complex mean_value_charge=Complex(0,0);

	for (int i=0; i<chargeVec.GetSize();i++)
	  mean_value_charge=mean_value_charge+chargeVec[i];
       
	mean_value_charge = mean_value_charge/(Double(chargeVec.GetSize()));
	//	std::cout<<"\nCHARGE VEC SIZE = "<< chargeVec.GetSize()<<"mean-value-charge = " << mean_value_charge << std::endl;

	JacobiMatrix[fstep][ind_param]=2.0*mean_value_charge; //pdes_[0]->getPDE_complexValuedCharge();    

	//   Info->PrintPiezoMat(*ptMaterial);
	//     pdes_[0]->WriteResultsInFile();     
     
      }		// end for over all frequencies
    std::cout << " \nCreateJacobian Line " << ind_param  <<std::endl;
      //for (int i=0;i<JacobiMatrix.GetSizeRow();i++)
      //	std::cout<<"F'("<<i<<")("<<ind_param<<")= "<< JacobiMatrix[i][ind_param]<<"; ";
    } //end loop over paramters
    
    for(int i=0;i<JacobiMatrix.GetSizeRow();i++)
      for (int j=0; j<JacobiMatrix.GetSizeCol();j++){
    	std::cout<<"F'("<<i<<")("<<j<<")= "<< JacobiMatrix[i][j]<<"; ";
    	if (j==JacobiMatrix.GetSizeCol()-1)
    	  std::cout<<"\n";
      }

  }            //end CreateJacobiMatrix


  void piezoParamIdent::testJacobiMatrix(Vector<Complex> & F_hat, Matrix<Complex> & JacobiMatrix, Vector<Double> & parameter,BCs * ptBCs,MaterialData * ptMaterial,Vector<Double> & parameterIncrement, Vector<Complex>& solElecPot,Vector<Complex> &solMechDispl){
    ENTER_FCN("piezoParamIdent::testJacobiMatrix");

    Vector<Complex> F_hat_incr(F_hat.GetSize());
    Matrix<Complex> approxJacobiMatrix(JacobiMatrix.GetSizeRow(), JacobiMatrix.GetSizeCol());
    Vector<Double> parameter_incr(parameter.GetSize());

    for(Integer i=0;i<parameter.GetSize();i++)
      parameter_incr[i]=10.1*parameter[i];

    updateMaterialData(parameter_incr, ptMaterial);
    createF(ptMaterial, ptBCs, F_hat_incr);
     for (int i=0;i<F_hat.GetSize();i++)
       std::cout<<"F("<<i+1<<")="<<F_hat[i]<< " <-> " << F_hat_incr[i]<<"; ";
    // createJacobiMatrix(ptMaterial, ptBCs, F_hat_incr, parameterIncrement,JacobiMatrix, solElecPot, solMechDispl);
     std::cout<<"\n  - - - - - - - - - - - - - - - \n JacobiMatrix <-> approxJacobiMatrix"<< std::endl;

    for (Integer i=0; i< JacobiMatrix.GetSizeRow();i++)
      for (Integer j=0; j< JacobiMatrix.GetSizeCol();j++){
	if (j<=parameter.GetSize())
	  approxJacobiMatrix[i][j]=0.5*(F_hat[i]-F_hat_incr[i])/(Complex(scaling[j]*parameter[j])-Complex(scaling[j]*parameter_incr[j]));
	else {
	  int jj=i-parameter.GetSize();
	  approxJacobiMatrix[i][j]=0.5*(F_hat[i]-F_hat_incr[i]);// /(Complex(scaling[jj]*parameter[jj])-Complex(scaling[jj]*parameter_incr[jj]));
								 }
	std::cout<<"F'("<<i<<")("<<j<<")= "<< JacobiMatrix[i][j]<<" <-> " <<approxJacobiMatrix[i][j]<< "; " ;
	if (j==JacobiMatrix.GetSizeCol()-1)
	  std::cout<<"\n " ;
      }
  }// end testJacobiMatrix

  void piezoParamIdent::createAdjointJacobiMatrix(Vector<Double> & parameterIncrement,Vector<Double> &  parameter, Matrix<Complex> & JacobiMatrix, Vector<Complex> & solElecPot,Vector<Complex> & solMechDispl, Vector<Double> & freqs, Matrix<Complex> & adjJacobiMatrix){
    ENTER_FCN("piezoParamIdent::createAdjointJacobiMatrix");
    //    std::cout<<"Adjoint Jacobian will be created ... "<<std::endl;
    adjJacobiMatrix.Resize(JacobiMatrix.GetSizeCol(),JacobiMatrix.GetSizeRow());
    for (int i=0;i<JacobiMatrix.GetSizeRow();i++)
      for (int j=0;j<JacobiMatrix.GetSizeCol();j++)
	adjJacobiMatrix[j][i] = std::conj(JacobiMatrix[i][j]);
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
	
    //loop over elements
    for (int actEl=0; actEl< elemssd.GetSize(); actEl++) {
      BaseFE * ptEl = elemssd[actEl]->ptElem;
      StdVector<Integer> connecth = elemssd[actEl]->connect;
		             
      Matrix<Double> ptCoord;
      ptGrid->GetCoordNodesElemMat(connecth, ptCoord, 0);
		    
      // map connect to PDE node numbers
      StdVector<Integer> connect_PDE;	   	    
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
	
      bilinearStiff->CalcComplexElementMatrix(ptCoord,elemMat,damp_beta,freqs[fstep]);

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

