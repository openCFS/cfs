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

#include <sstream>
#include <iomanip>



#include "Utils/tools.hh"
#include <PDE/pdes_header.hh>

#include <Driver/piezoParamIdent.hh>



//#include "/../OLAS/algsys/basesystem.hh"
//#include "DataInOut/piezoParameterData.hh"



namespace CoupledField
{

  // xxxxxxxxxxxxxxxxxxxxxxxxxxxxx NEWTON CG 3 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  void piezoParamIdent::tichonov(){   // here we scale in advance, just once before NewtonCG and afterwards
    ENTER_FCN("piezoParamIdent::tichonov");
    std::cout<<"\n Entering piezoParamIdent::tichonov ... "<< std::endl;

    Integer nrNewtonIterations=0;
    Integer backtrackIterator=0;
    Integer i;
    Double misfit, misfit_new, lin_misfit;

    MaterialData * ptMaterial=pdes_[0]->getPDEMaterialData();   // Pointer to MaterialData

    ptBCs = pdes_[0]->getPDE_BCs();
    
    Vector<Complex> y_hat_F_hat (nrMeasuredData);
    Vector<Complex> JacFS(nrMeasuredData);
    Vector<Complex> linres(nrMeasuredData);
    Vector<Complex> basC;
    bas.Resize(nrParameter);
    basC.Resize(nrParameter);

    for(Integer i=0;i<nrParameter;i++)
      bas[i]=1.0;

    Vector<Complex> step(nrParameter);

    createF(ptMaterial, ptBCs, F_hat);

    for (i=0; i<nrMeasuredData;i++)
      y_hat_F_hat[i]=y_hat[i]-F_hat[i];

    misfit=sqrt(realA2norm(y_hat_F_hat)); 
    misfit_new=misfit;
    lin_misfit;
    Double eta = 0.9;

    std::cout << "\n misfit = |y-F| =  " <<  misfit <<std::endl;

    while((misfit>tau*delta||nrNewtonIterations<25)&&nrNewtonIterations<9){ // Newton

      nrNewtonIterations++;

      createJacobiMatrix2(JacobiMatrix);

      createAdjointJacobiMatrix(JacobiMatrix,adjJacobiMatrix);

      adjJacobiMatrix.Mult(y_hat_F_hat,basC);
      
      std::cout<<"\n update bas"<<std::endl;
      for(Integer i=0;i<nrParameter;i++)
	bas[i]=basC[i].real();


      Double reg_alpha=5.0e-16;

      Matrix<Complex> JacobiMatrixNE(JacobiMatrix.GetSizeCol(), JacobiMatrix.GetSizeCol());
      Vector<Complex> y_hat_F_hatNE(JacobiMatrix.GetSizeCol());

       adjJacobiMatrix.Mult(JacobiMatrix,JacobiMatrixNE);
       adjJacobiMatrix.Mult(y_hat_F_hat, y_hat_F_hatNE);

      // Choice of regularisation parameter alpha!
       Integer innerIter=0;

       while (lin_misfit>eta*misfit_new&&innerIter<2||innerIter<1){
	 innerIter++;

	 std::cout<<"\n --> Value of reg_alpha = " << reg_alpha << std::endl;
	 std::cout<<"\n lin_misfit = " << lin_misfit << "; misfit_new = "<<misfit_new<<std::endl;

       std::cout<<"\n Build up normal equation"<<std::endl;

       for (Integer i=0;i<JacobiMatrixNE.GetSizeRow();i++)
	 for (Integer j=0;j<JacobiMatrixNE.GetSizeCol();j++){
	   JacobiMatrixNE[i][j]=JacobiMatrixNE[i][j].real();
	   if(i==j)
	     JacobiMatrixNE[i][i]=JacobiMatrixNE[i][i]+reg_alpha;
	   // std::cout << std::setprecision(15);
	    std::cout<<JacobiMatrixNE[i][j]<<" ";
	   if(j==JacobiMatrix.GetSizeCol()-1)
	    std::cout<<"; \n";
	 }

       std::cout<<"\n"<<std::endl;

       for (Integer j=0;j<y_hat_F_hatNE.GetSize();j++){
	   y_hat_F_hatNE[j]=y_hat_F_hatNE[j].real();
	   //std::cout<<y_hat_F_hatNE[j].real()<<" ";
       }

       JacobiMatrixNE.DirectSolve(step,y_hat_F_hatNE);
      
        std::cout<< "\n s determined with direct solver ...\n " <<std::endl;
              for (int i=0;i<parameter.GetSize();i++){ 
        	std::cout<<"s("<<i<<")="<<step[i]<<"; \t";
		parameter_new[i]=parameter[i]; //+(1.0/scaling[i])*step[i].real();        
		//		parameter_new=parameter;
		//	std::cout<<"parameter("<<i<<")="<<parameter[i]<<"; \t";
              }	     

	parameter_new[9]=parameter[9]+(1.0/scaling[9])*step[9].real();        // eps33
	parameter_new[1]=parameter[1]+(1.0/scaling[1])*step[1].real(); // c33
	parameter_new[7]=parameter[7]+(1.0/scaling[7])*step[7].real(); //e33
	      
        std::cout<<"\n"<<std::endl;
	updateMaterialData(parameter_new,ptMaterial);

	createF(ptMaterial, ptBCs, F_hat);
	y_hat_F_hat=y_hat-F_hat;
	//	misfit = realA2norm(y_hat_F_hat);
	JacobiMatrix.Mult(step,JacFS);
	linres = JacFS-y_hat_F_hat;
	lin_misfit = sqrt(realA2norm(linres));
	misfit_new = sqrt(realA2norm(y_hat_F_hat));
	reg_alpha = 1.5 * reg_alpha;
	
	for (Integer i=0;i<nrParameter;i++){
	  std::cout<<"Paramter_new["<<i<<"]= " << parameter_new[i]<<std::endl; 
      }      

 }


      

      std::cout<<"\n"<<std::endl;

      // precautionary measure, if paramters tend to far away from initial values ...
//       for(Integer i=0;i<nrParameter;i++)
// 	if (std::abs(parameter_new[i]-parstart[i])>0.5*std::abs(parstart[i])){
// 	  parameter_new[i]=parameter[i];
// 	  std::cout<<"\n parameter("<<i<<") was reset to value " << parameter[i] <<std::endl;
// 	}

      updateMaterialData(parameter_new, ptMaterial);

      createF(ptMaterial, ptBCs, F_hat);
      y_hat_F_hat=y_hat-F_hat;
      misfit = realA2norm(y_hat_F_hat);
      parameter=parameter_new;
	

    } // end while
  }
}// end namepsace
