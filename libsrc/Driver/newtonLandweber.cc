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

    ENTER_FCN("piezoParamIdent::Newton()");

    // Settings for Newton-CG - routine
    Integer nrIterations=0;
    Integer nLandweber=0;

    Double tol=0.001, tolCG=0.001, err=1.0, res_norm, theta;


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

    while (nrIterations<3) {
      nrIterations++;
      std::cout<<"\n NewtonLandweber ... Newton-Iteration-Nr = "<<nrIterations<<std::endl;
      s.Resize(nrParameter);

      // Create the Matrices F, F', F*
      createF(ptMaterial, ptBCs, F_hat);
      // std::cout<<"Parameter-to-solution-map Matrix is built up!"<<std::endl;
      createJacobiMatrix(ptMaterial, ptBCs, F_hat, parameterIncrement,JacobiMatrix, solElecPot, solMechDispl);
      //std::cout<<"JacobiMatrix was created ..."<<std::endl;
      createAdjointJacobiMatrix(JacobiMatrix,adjJacobiMatrix);
      //    std::cout<<"Adjoint JacMatrix was created ..."<<std::endl;
      // std::cout<<"\n\ntest Jacobi Matrix ..."<<std::endl;

      //testJacobiMatrix(F_hat, JacobiMatrix, parameter, ptBCs, ptMaterial,parameterIncrement, solElecPot, solMechDispl);


      // TEST MAT_MULT 

      // FF = F' F'*
      Matrix<Complex> FF(adjJacobiMatrix.GetSizeRow(),JacobiMatrix.GetSizeCol());
      adjJacobiMatrix.Mult(JacobiMatrix,FF);

      Double symsysmat = calcEuclidianMatrixNorm(FF);
      std::cout<<"Norm of symmetric system matrix = " << symsysmat <<std::endl;


      //     for(int i=0;i<FF.GetSizeRow();i++)
      //       for (int j=0; j<FF.GetSizeCol();j++){
      //     	std::cout<<"FF("<<i<<")("<<j<<")= "<< FF[i][j]<<"; \t";
      //     	if (j==FF.GetSizeCol()-1)
      //     	  std::cout<<"\n";
      //       }



      //we need the Identity of size nrParameter x nrParameter and some other temporary matrices
      Matrix<Complex> Identity(nrParameter,nrParameter);
      for (Integer i=0;i<nrParameter;i++)
	Identity[i][i]=1;
      Matrix<Complex> adjFJacF (nrParameter, 2*nrMeasuredData);
      Vector<Complex> adjFRes (nrParameter);
      Vector<Complex> act_res(2*nrMeasuredData);
      Matrix<Complex> I_adjFF(nrParameter,nrParameter);
      Vector<Complex> adjFF_res(nrParameter);
      Vector<Complex> JacFs(2*nrMeasuredData);
      Vector<Complex> JacFs_res(2*nrMeasuredData);
      Vector<Complex> s_old(nrParameter);

      if (considerMechDeformation==FALSE){
	adjFJacF.Resize(nrParameter, nrMeasuredData);
	act_res.Resize(nrMeasuredData);
	JacFs.Resize(nrMeasuredData);
	JacFs_res.Resize(nrMeasuredData);
      }

 
      Complex Jfs = Complex(0.0,0.0);
      Double w=0.9;
      Double normJacMat, old_resid, new_resid;
 
      normJacMat=calcEuclidianMatrixNorm(JacobiMatrix);
      while (w>=1/normJacMat)
	w=0.9*w;
      std::cout<<"Choice of w = "<< w << "; Norm JacMat = " << normJacMat << std::endl;

      // y_hat-F_hat
      for (Integer i=0;i<y_hat.GetSize();i++)
	act_res[i]=y_hat[i]-F_hat[i];
	  
      old_resid=a2norm(act_res);
      new_resid = old_resid;

      std::cout<<"Before Landweber, nLandweber = " << nLandweber << ", new_resid = " << new_resid << ", eta * old_resid = " << eta*old_resid <<std::endl;
    
      // LANDWEBER ITERATION   s^k+1 = s^k - w F'*(F's^k - (y-F))
      while(new_resid>eta*old_resid&&nLandweber<4){
	nLandweber++;
	std::cout <<"\n Here starts the Landweber Iteration - Nr " << nLandweber<< std::endl;

	// F' * s^k
	JacobiMatrix.Mult(s,JacFs);

	// F'sk - (y_hat-F_hat)
	for (Integer i=0;i<y_hat.GetSize();i++){
	  act_res[i]=y_hat[i]-F_hat[i];
	  JacFs_res[i]=JacFs[i]-act_res[i];
	}
	
	// F'*(F's^k - (y-F))
	adjJacobiMatrix.Mult(JacFs_res,adjFF_res);


	std::cout<< "\n Parameter increments after Landweber Step " << nLandweber <<std::endl;
	//  s^k+1 = s^k - w F'*(F's^k - (y-F))
	//std::cout<<"adjFF_res.GetSize() = " << adjFF_res.GetSize()<<std::endl;


	for(Integer i=0;i<nrParameter;i++){
	  s[i]=s[i]-w*adjFF_res[i];
	  std::cout<<"s("<<i<<")= "<<s[i]<<"; "<<std::endl;
	}

	// JacFs = F'*s
	JacobiMatrix.Mult(s,JacFs);

	//F'(p^k)(s^k)-(y-F(p^k))
	for (Integer i=0;i<adjJacobiMatrix.GetSizeRow();i++)
	  act_res[i]-=JacFs[i];
	
	//	for (Integer i=0;i<act_res.GetSize();i++)
	// std::cout<<"act_res= "<< act_res[i]<<std::endl;
	
	new_resid=a2norm(act_res);
	std::cout<<"\n new_resid = " << new_resid << ";  old_resid = " << old_resid << "; tau*delta= " << tau*delta << std::endl;
	s_old=s;
	  	
      } // end while Landweber ...

      nLandweber=0;

      parameter_new = parameter;

      Double old_resid2=old_resid;
      Double new_resid2=new_resid;

      // backtracking(eta, theta, s, old_resid2, new_resid2); 
      
      // if no backtracking is specified, please include the following lines!
      for (Integer i=0;i<nrParameter;i++){
 	parameter_new[i]=scaling[i]*parameter[i];
 	parameter_new[i]+=s[i].real();
 	parameter_new[i]=1/scaling[i]*parameter[i];
      }
      
      updateMaterialData(parameter, ptMaterial);
      createF(ptMaterial, ptBCs, F_hat);

      for (Integer i=0;i<y_hat.GetSize();i++)
	act_res[i]=y_hat[i]-F_hat[i];
      new_resid=(a2norm(act_res));


    } // end while nrIterations
   
  } // end NewtonLandweber
} // end namespace
