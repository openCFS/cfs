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

  // xxxxxxxxxxxxxxxxxxxxxxxxxx LANDWEBER  - ITERATION xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  void piezoParamIdent::NewtonLandweber(){

    ENTER_FCN("piezoParamIdent::NewtonLandweber()");
    std::cout<<"\n Entering piezoParamIdent::NewtonLandweber()"<<std::endl;
    std::cout<<"\n Landweber 1"<<std::endl;

    // Settings for Newton-CG - routine
    Integer nrIterations=0;
    Integer nLandweber=0;

    Double tol=0.001, tolCG=0.001, err=1.0, res_norm, theta;
    Double gamma = 0.8;
    tau=1.0;


    MaterialData * ptMaterial=pdes_[0]->getPDEMaterialData();   // Pointer to MaterialData
  
    updateMaterialData(parameter, ptMaterial);         //Writes initial guesses of parameters (read from MeasuredData.dat) to system

    ptBCs = pdes_[0]->getPDE_BCs();                             // Pointer to BCs
    ptAlgsys = pdes_[0]->getPDE_algsys();                       //Pointer to AlgebraicSystem
    Integer numElems = pdes_[0]->getPDE_numElems();
    dofs=pdes_[0]->getPDE_dofspernode();
    numNodes= pdes_[0]->getPDE_numPDENodes();

    //we do a bit of scaling here ...
    //   for(int i=0;i<scaling.GetSize();i++)
    //    parameter[i]=parameter[i]*scaling[i];

    eta=0.9;
    Double eta_acc=0.5;

      Double normJacMat, old_res_outer, new_res_inner, old_res_inner, new_res_outer;
      Matrix<Complex> Identity(actNrParameter, actNrParameter);
      //we need the Identity of size nrParameter x nrParameter and some other temporary matrices
      for (Integer i=0;i<actNrParameter;i++)
	Identity[i][i]=1;
      //  std::cout<<"\n Landweber 2"<<std::endl;

      Matrix<Complex> adjFJacF (actNrParameter, 2*nrMeasuredData);
      Vector<Complex> adjFRes (actNrParameter);
      Vector<Complex> act_res(2*nrMeasuredData);
      Matrix<Complex> I_adjFF(actNrParameter,actNrParameter);
      Vector<Complex> adjFF_res(actNrParameter);
      Vector<Complex> JacFs(2*nrMeasuredData);
      Vector<Complex> JacFs_res(2*nrMeasuredData);
      Matrix<Complex> s_old(maxNumberInnerLoops,actNrParameter);
      Vector<Complex> s_temp(actNrParameter);
      Matrix<Complex> FF;
      Vector<Double> stepR(actNrParameter);

      //  std::cout<<"\n Landweber 3"<<std::endl;
      if (considerMechDeformation==FALSE){
	adjFJacF.Resize(actNrParameter, actNrParameter);
	act_res.Resize(nrMeasuredData);
	JacFs.Resize(nrMeasuredData);
	JacFs_res.Resize(nrMeasuredData);
      }

      //  std::cout<<"\n Landweber 4"<<std::endl;
      createF(ptMaterial, ptBCs, F_hat,TRUE);
      act_res = y_hat-F_hat;
      new_res_outer=old_res_outer=a2norm(act_res);


      //    std::cout<<"maxNumberNewtonLoops = " << maxNumberNewtonLoops <<std::endl;
      //      getchar();
    while (new_res_outer<=old_res_outer && nrIterations<maxNumberNewtonLoops) {

 for (Integer i=0;i<parameter.GetSize();i++)
	scaling[i]=1.0;
    
      nrIterations++;
      std::cout<<"\n NewtonLandweber ... Newton-Iteration-Nr = "<<nrIterations<<std::endl;
      s.Resize(actNrParameter);
      s_old.Resize(maxNumberInnerLoops,actNrParameter);
      s_temp.Resize(actNrParameter);


      // Create the Matrices F, F', F*
      createF(ptMaterial, ptBCs, F_hat,TRUE);
      // std::cout<<"Parameter-to-solution-map Matrix is built up!"<<std::endl;
      //      createJacobiMatrix(ptMaterial, ptBCs, F_hat, parameterIncrement,JacobiMatrix, solElecPot, solMechDispl);
      createJacobiMatrix2(JacobiMatrix);

       std::cout<<"\n Landweber 2"<<std::endl;
      // testJacobiMatrix(F_hat, JacobiMatrix, parameter, ptBCs, ptMaterial,parameterIncrement, solElecPot, solMechDispl);
      //JacobiMatrix = approxJacobiMatrix;
      std::cout<<"JacobiMatrix was created ..."<<std::endl;
      createAdjointJacobiMatrix(JacobiMatrix,adjJacobiMatrix);
      std::cout<<"Adjoint JacMatrix was created ..."<<std::endl;
      // std::cout<<"\n\ntest Jacobi Matrix ..."<<std::endl;
      //  std::cout<<"\n Landweber 6"<<std::endl;
      //testJacobiMatrix(F_hat, JacobiMatrix, parameter, ptBCs, ptMaterial,parameterIncrement, solElecPot, solMechDispl);

      // TEST MAT_MULT 

      // FF = F' F'*
      FF.Resize(adjJacobiMatrix.GetSizeRow(),JacobiMatrix.GetSizeCol());
      adjJacobiMatrix.Mult(JacobiMatrix,FF);

      Double symsysmat = calcEuclidianMatrixNorm(FF);
      //    std::cout<<"Norm of symmetric system matrix = " << symsysmat <<std::endl;


      //     for(int i=0;i<FF.GetSizeRow();i++)
      //       for (int j=0; j<FF.GetSizeCol();j++){
      //     	std::cout<<"FF("<<i<<")("<<j<<")= "<< FF[i][j]<<"; \t";
      //     	if (j==FF.GetSizeCol()-1)
      //     	  std::cout<<"\n";
      //       }

       std::cout<<"\n Landweber 3"<<std::endl;
      Complex Jfs = Complex(0.0,0.0);
      Double w=0.9;

      //  std::cout<<"\n Landweber 7"<<std::endl;
 
      normJacMat=calcEuclidianMatrixNorm(JacobiMatrix);
      while (w>=1/(normJacMat*normJacMat))
	w=0.9*w;
      // w=0.9;
      //    std::cout<<"Choice of w = "<< w << "; Norm JacMat*AdjJacMat = " << normJacMat << std::endl;
      //      getchar();

      // y_hat-F_hat
      for (Integer i=0;i<y_hat.GetSize();i++)
	act_res[i]=y_hat[i]-F_hat[i];
       std::cout<<"\n Landweber 4"<<std::endl;
	  
      old_res_outer=a2norm(act_res);
      new_res_outer = old_res_outer;

      new_res_inner=old_res_outer;
      old_res_inner=old_res_outer;

      //    std::cout<<"Before Landweber, nLandweber = " << nLandweber << ", new_res = " << new_res_outer << ", eta * old_resid_outer = " << eta*old_res_outer <<std::endl;
    
      // LANDWEBER ITERATION   s^k+1 = s^k - w F'*(F's^k - (y-F))
      while(new_res_inner<=old_res_inner&&nLandweber<maxNumberInnerLoops){
	nLandweber++;
	//	std::cout <<"\n Here starts the Landweber Iteration - Nr " << nLandweber<< std::endl;
	old_res_inner=new_res_inner;

	// F' * s^k
	JacobiMatrix.Mult(s,JacFs);
       std::cout<<"\n Landweber 5"<<std::endl;

	// F'sk - (y_hat-F_hat)
	for (Integer i=0;i<nrMeasuredData;i++){
	  act_res[i]=y_hat[i]-F_hat[i];
	  JacFs_res[i]=JacFs[i]-act_res[i];
	}
	
	// F'*(F's^k - (y-F))
	adjJacobiMatrix.Mult(JacFs_res,adjFF_res);
	
	//	std::cout<< "\n Parameter increments after Landweber Step " << nLandweber <<std::endl;
	//  s^k+1 = s^k - w F'*(F's^k - (y-F))
	//std::cout<<"adjFF_res.GetSize() = " << adjFF_res.GetSize()<<std::endl;
       std::cout<<"\n Landweber 6"<<std::endl;

       std::cout<<s_old<<std::endl;
	
       // if accelerated
       if (FALSE){
	 for(Integer i=0;i<actNrParameter;i++){
	   if (nLandweber>2)
	     s[i]=s[i]+eta_acc*(s[i]-s_old[nLandweber-2][i])-100.0*w*adjFF_res[i];
	   else
	     s[i]=s[i]-100.0*w*adjFF_res[i];
// 	   std::cout<<"s("<<i<<")= "<<s[i]<<"; "<<std::endl;
	 }
	 for(Integer i=0;i<actNrParameter;i++)
	   s_old[nLandweber][i]=s[i];
       }
       

	if (TRUE){
	  for(Integer i=0;i<actNrParameter;i++){
	    s[i]=s[i]-100.0*w*adjFF_res[i];
	    //std::cout<<"s("<<i<<")= "<<s[i]<<"; "<<std::endl;
	  }
	}
       std::cout<<"\n Landweber 6.4"<<std::endl;
	// JacFs = F'*s
	JacobiMatrix.Mult(s,JacFs);
       std::cout<<"\n Landweber 6.5"<<std::endl;
	//F'(p^k)(s^k)-(y-F(p^k))
	for (Integer i=0;i<nrMeasuredData;i++)
	  act_res[i]-=JacFs[i];
       std::cout<<"\n Landweber 6.6"<<std::endl;
	
	//	for (Integer i=0;i<act_res.GetSize();i++)
	// std::cout<<"act_res= "<< act_res[i]<<std::endl;
	
	new_res_inner=a2norm(act_res);
       std::cout<<"\n Landweber 7"<<std::endl;

	if (new_res_inner>old_res_inner){
	  std::cout << " \n !! New_res_inner is worse than old_res_inner!! "<< std::endl;
	  break;
	  getchar();
	}
	
	//	std::cout<<"\n new_res_inner = " << new_res_inner << ";  old_res_outer = " << old_res_outer << "; tau*delta= " << tau*delta << std::endl;
		  	
	//	std::cout<<"\n end of inner Landweber Iter ..."<<std::endl;
      } // end while Landweber ...


      nLandweber=0;

      parameter_new = parameter;

      Double old_resid2=old_res_outer;
      Double new_resid2=new_res_outer;

      // backtracking(et , theta, s, old_resid2, new_resid2); 

      theta = 0.5;
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

    for (Integer i=0;i<actNrParameter;i++)
      stepR[i]=s[i].real();
    
    parameter_new=parameter;
      
      setNewParameterSet(parameter, parameter_new, scaling, theta, stepR, whichParameterToUpdate);

      // if no backtracking is specified, please include the following lines!
    for (Integer i=0;i<nrParameter;i++){
//  	parameter_new[i]=scaling[i]*parameter[i];
//  	parameter_new[i]+=s[i].real();
//  	parameter[i]=1/scaling[i]*parameter_new[i];
 	std::cout<<" paramter("<<i<<") = " << parameter_new[i]<<std::endl;
       }
    parameter=parameter_new;
      
      updateMaterialData(parameter, ptMaterial);
      createF(ptMaterial, ptBCs, F_hat,FALSE);


      for (Integer i=0;i<y_hat.GetSize();i++)
      	act_res[i]=y_hat[i]-F_hat[i];
      new_res_outer=(a2norm(act_res));
      std::cout<<"\n new_res_outer = " << new_res_outer <<std::endl;

      if (new_res_outer>=old_res_outer){
	std::cout<<"\n Warning: residual norm gets worse!" <<std::endl;
	//	getchar();
	//parameter=parameter_new;
	//	NewtonLandweber();
      }
      
    } // end while nrIterations

    //    delete adjFJacF;
   
  } // end NewtonLandweber
} // end namespace
