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

  // xxxxxxxxxxxxxxxxxxxxxxxxxxxxx NEWTON CG 3 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  void piezoParamIdent::tichonov(){   // here we scale in advance, just once before NewtonCG and afterwards
    ENTER_FCN("piezoParamIdent::NewtonCG3");
    std::cout<<"\n Entering piezoParamIdent::NewtonCG3 ... "<< std::endl;

    Integer nrNewtonIterations=0;
    Integer backtrackIterator=0;
    Integer i;

    MaterialData * ptMaterial=pdes_[0]->getPDEMaterialData();   // Pointer to MaterialData

    ptBCs = pdes_[0]->getPDE_BCs();
    
    Vector<Complex> y_hat_F_hat (nrMeasuredData);

    //    real_misfit=norm_res;
    //real_misfit_new=norm_res_new;

    Vector<Complex> step(nrParameter);
    Vector<Double> parstart (nrParameter);

    Vector<Complex> res0(nrMeasuredData);

    //    Vector<Complex> y_hat_F_hat(nrMeasuredData);

    Double misfit, misfit0,misfitnew, normres02, normres2, normresNEold2, normresNEold02, normresNE2;

     //    updateMaterialData(parameter,ptMaterial);

    // Create the Matrices F, F', F*
    createF(ptMaterial, ptBCs, F_hat);

    for (i=0; i<nrMeasuredData;i++)
      y_hat_F_hat[i]=y_hat[i]-F_hat[i];

    res=y_hat-F_hat;
      
    misfit=sqrt(realA2norm(y_hat_F_hat)); 
    misfit0=misfit;

    std::cout << "\n misfit = |y-F| =  " <<  misfit <<std::endl;

    while((misfit>tau*delta||nrNewtonIterations<25)&&nrNewtonIterations<25){ // Newton
      res0=res;
      normres02=norm2Real(res0);

      nrNewtonIterations++;

      createJacobiMatrix2(JacobiMatrix);
      createAdjointJacobiMatrix(JacobiMatrix,adjJacobiMatrix);

    } // end while
  }
}// end namepsace
