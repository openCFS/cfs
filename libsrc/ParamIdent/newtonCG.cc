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

  void piezoParamIdent::NewtonCG3(){   // here we scale in advance, just once before NewtonCG and afterwards
    ENTER_FCN("piezoParamIdent::NewtonCG3");
    std::cout<<"\n Entering piezoParamIdent::NewtonCG3 ... "<< std::endl;

    Integer nrNewtonIterations=0;
    Integer backtrackIterator=0;

    MaterialData * ptMaterial=pdes_[0]->getPDEMaterialData();   // Pointer to MaterialData

    ptBCs = pdes_[0]->getPDE_BCs();

    Double eta_max, eta_new, t, aa,b,c, theta_min, theta_max, gamma, al;
    Double alpha, beta, tau;
    tau = 1.1;
    Integer i,j;
    Double theta=0.5;
    delta = 0.1;

    t=1.0e-4;
    eta_max=0.9;
    eta=0.9;
    theta_min=0.1;
    theta_max=0.5;
    gamma=0.5;
    al=1.5;
    
    
    Vector<Complex> y_hat_F_hat (nrMeasuredData);

    //    real_misfit=norm_res;
    //real_misfit_new=norm_res_new;

    Vector<Complex> step(nrParameter);
    Vector<Double> parstart (nrParameter);
    Vector<Double> res_NE(nrParameter);
    Vector<Double> res_NE_rescaled(nrParameter);
    Vector<Complex> res_NE_C(nrParameter);
    //    Vector<Double> bas(nrParameter);
    bas.Resize(nrParameter);
    Vector<Complex> basbar(nrMeasuredData);
    Vector<Complex> basbar_old(nrMeasuredData);
    // Vector<Complex> res_old(nrMeasuredData);
    Vector<Complex> res(nrMeasuredData);
    Vector<Complex> JMats(nrMeasuredData);
    Vector<Complex> res0(nrMeasuredData);

    //    Vector<Complex> y_hat_F_hat(nrMeasuredData);
    Complex new_res_res0;
    Double misfit, misfit0,misfitnew, normres02, normres2, normresNEold2, normresNEold02, normresNE2;


//     scaling[0]=1.0e-11;    scaling[1]=1.0e-11;    scaling[2]=1.0e-11;    scaling[3]=1.0e-11;    scaling[4]=1.0e-10;
//       scaling[5]=0.1;     scaling[6]=1.0;     scaling[7]=0.1;
//       scaling[8]=1.0e+8; scaling[9]=1.0e+9;

// NO scaling in NewtonCG!! 
    for (Integer i=0;i<parameter.GetSize();i++)
      scaling[i]=1.0;
    //    parameter[i]=scaling[i]*parameter[i];

     parstart=parameter;

     //    updateMaterialData(parameter,ptMaterial);

    // Create the Matrices F, F', F*
    createF(ptMaterial, ptBCs, F_hat);

    for (i=0; i<nrMeasuredData;i++)
      y_hat_F_hat[i]=y_hat[i]-F_hat[i];

    for (i=0;i<nrParameter;i++)
      bas[i]=1.0;

    

    res=y_hat-F_hat;
      
    misfit=norm2Real(y_hat_F_hat); 
    misfit0=misfit;

    std::cout << "\n misfit = |y-F| =  " <<  misfit <<std::endl;

    *piezoLog << "\n Initial misfit ||y-F|| = "<<misfit<<std::endl;

    nrNewtonIterations=0;

    // NEWTON ITERATION -- outer Loop!!
    while((misfit>tau*delta||nrNewtonIterations<5)&&nrNewtonIterations<25){ // Newton
      *piezoLog << "\n Newton-Iteration: " << nrNewtonIterations <<std::endl;
      *piezoLog <<"------------------------"<<std::endl;
      nrNewtonIterations++;

      res0=res;
      normres02=norm2Real(res0);

      std::cout<<"\n Newton-Iteration "<< nrNewtonIterations <<std::endl;
      step.Resize(nrParameter);
    
      res=y_hat-F_hat;

      // Create the matrices ...
      createJacobiMatrix2(JacobiMatrix);
      createAdjointJacobiMatrix(JacobiMatrix,adjJacobiMatrix);
      // testJacobiMatrix(F_hat, JacobiMatrix, parameter, ptBCs, ptMaterial,parameterIncrement, solElecPot, solMechDispl);

    for (Integer i=0;i<parameter.GetSize();i++)
      scaling[i]=1.0;

      //  JacobiMatrix.MatVecMult_CD(step,JMats); // is not neccessary, since we start with s_i=0;
      //res = res0-JMats;
      
	std::cout<<"\n \nres = res0 - JMats: "<<std::endl;
	for (Integer i=0;i<basbar.GetSize();i++){
	  std::cout<<"res="<<res[i]<<"; ";
	  //  std::cout<<"JMats="<<JMats[i]<<"; "<<std::endl;
	}
	std::cout<<"\n  " <<std::endl;	  

      normres2 = norm2Real(res);

      adjJacobiMatrix.Mult(res, res_NE_C);

      std::cout<<"\nbas = res NE = (F* res).real(): "<<std::endl;
      for (Integer i=0;i<res_NE_C.GetSize();i++){
	res_NE[i]=(1.0/scaling[i]*scaling[i])*res_NE_C[i].real();
	res_NE_rescaled[i]=scaling[i]*res_NE[i];
	 std::cout<<"res_NE="<<res_NE[i]<<"; ";
      }

      bas = res_NE;
      
      std::cout<<"\n  " <<std::endl;	  

      normresNEold2 = POW(a2norm(res_NE_rescaled),2);

      basbar_old.Resize(nrParameter);


      // ----- INNER LOOP --- CG-ITERATION ---
      Integer nrCGIter=0;
      while((normres2>eta*eta*normres02)&&(nrCGIter<10)){ // CG
	std::cout<<"\n CG Iteration " << nrCGIter << std::endl;
	nrCGIter++;

	JacobiMatrix.MatVecMult_CD(bas,basbar);
	
	std::cout<<"\nbasbar = F' bas:  "<<std::endl;
	for (Integer i=0;i<basbar.GetSize();i++){
	  basbar[i]=basbar[i].real();
	  std::cout<<"basbar="<<basbar[i]<<"; ";
	}
	std::cout<<"\n "<<std::endl;	  

	basbar_old=basbar;

	Double norm_bas_bar = POW(a2norm(basbar),2);

	alpha=normresNEold2 / norm_bas_bar;

	std::cout<<"\n alpha = " << alpha << ",\t || normresNEold || = " << normresNEold2 << ",\t ||bas_bar|| = " << norm_bas_bar << std::endl;

	std::cout<<"\nstep = step + alpha *bas: "<<std::endl;
	for (Integer i=0;i<nrParameter;i++){
	  step[i]+=alpha*bas[i];
	  std::cout<<"step("<<i<<")= "<< step[i] << "; \t ";
	}
	

	std::cout<<"\n\nres = res - alpha*bas_bar: "<<std::endl;
	for (Integer i=0;i<res.GetSize();i++){
	  res[i]=res[i]-alpha*basbar[i];
	  res[i]=res[i].real();
	  std::cout<<"res = "<<res[i]<<"; ";
	}
	std::cout<<"\n"<<std::endl;

	//	std::cout<<"\n res.GetSize() " << res.GetSize() << "; res_NE_NEW_COMP.GetSize()= "<< res_NE_new_compl.GetSize()<<std::endl;

	adjJacobiMatrix.Mult(res,res_NE_C);

	std::cout<<"\nres_NE = (F*' res).real():"<<std::endl;
	for(Integer i=0;i<res_NE_new.GetSize();i++){
	  res_NE[i]=(1.0/(scaling[i]*scaling[i]))*res_NE_C[i].real();
	   res_NE_rescaled[i]=res_NE[i]*scaling[i];
	//	  res_NE_new[i]=1.0/(scaling[i]*scaling[i])*res_NE_new_compl[i].real();
	  std::cout<<"resNE = " << res_NE[i];
	}

	normresNE2=POW(a2norm(res_NE_rescaled),2);

	std::cout<<"\n\n normresNE2 =  " << normresNE2<<std::endl; 

	beta = normresNE2/normresNEold2;
	
	std::cout<<"\n beta = " << beta << " = " << normresNE2 << " / " << normresNEold2<< std::endl; 

	normresNEold2 = normresNE2;

	std::cout<<"\nbas = resNE + beta * bas: "<<std::endl;	
	for(Integer i=0;i<nrParameter;i++){
	  bas[i]= res_NE[i]+beta*bas[i];
	  std::cout<<"bas = "<<bas[i]<<"; ";
	}

	normres2=norm2Real(res);

      }// end while CG
	*piezoLog << "\t Number of CG - Iterations performed " << nrCGIter <<std::endl;

      parameter_new=parameter; // no update for other parameters

      Matrix<Double> *matMat = ptMaterial->GetMatrix();

      scaling[1]=1.0/((*matMat)[2][2]);
      scaling[7]=1.0/((*matMat)[8][2]);
      scaling[9]=1.0/((*matMat)[8][8]);
      
     	parameter_new[9]=parameter[9]+(1.0/scaling[9])*step[9].real();        // eps33
	parameter_new[1]=parameter[1]+(1.0/scaling[1])*step[1].real(); // c33
	parameter_new[7]=parameter[7]+(1.0/scaling[7])*step[7].real(); //e33
	      
        std::cout<<"\n"<<std::endl;


      std::cout<<"\n"<<std::endl;

      for (Integer i=0;i<nrParameter;i++){
// 	std::cout<<"step("<<i<<")= " << step[i]<< "; \t";
// 	//	parameter_new[i]=parameter[i]; //*scaling[i];
// 	//	parameter_new[i]=parameter[i]+step[i];
// 	//	parameter_new[i]=(1.0/scaling[i])*parameter_new[i];
	  std::cout<<"Paramter_new["<<i<<"]= " << parameter_new[i]<<std::endl; 
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

      for(Integer i=0;i<F_hat.GetSize();i++)
	y_hat_F_hat[i]=y_hat[i]-F_hat[i];


      //	backtracking(eta, theta, s, a, a_lin_new);
      misfitnew=norm2Real(y_hat_F_hat);
      *piezoLog << "\t misfit before Linesearch " << misfitnew <<std::endl;

      backtrackIterator=0;

      while ((misfitnew>(1-t*(1-eta))*misfit) && backtrackIterator<15 && eta<eta_max || backtrackIterator<2){ // Liniensuche
	std::cout<<"\n backtracking ... "<< backtrackIterator<< std::endl;
	backtrackIterator++;

	// choose theta:
	aa=misfit*misfit;
	new_res_res0=Complex(0.0,0.0);

	for(Integer i=0;i<nrMeasuredData;i++)
	  new_res_res0=new_res_res0+y_hat_F_hat[i]*(res[i].real()-y_hat_F_hat[i]);

	b=2.0*new_res_res0.real();
	std::cout<<"\n b = " << b << std::endl;
	c=misfitnew*misfitnew-b-aa;
	std::cout<<"\n c = "<<c<<std::endl;
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
	      if ((aa+b*theta_min+c*theta_min*theta_min)>=(aa+b*theta_max+c*theta_max*theta_max))
		theta=theta_max;
	      else theta=theta_min;
	    }
	} // end else if c>0
	//	theta = 0.35;

	std::cout<<"\n Choice of theta = " << theta<<std::endl;
	*piezoLog<<"\t Choice of theta = " << theta<<std::endl;

	std::cout<<"\n"<<std::endl;

	std::cout<<"\n Parameter after backtracking ... " <<std::endl;
	
	parameter_new[9]=parameter[9]+(1.0/scaling[9])*theta*step[9].real();        // eps33
	parameter_new[1]=parameter[1]+(1.0/scaling[1])*theta*step[1].real(); // c33
	parameter_new[7]=parameter[7]+(1.0/scaling[7])*theta*step[7].real(); //e33
	
	for(Integer i=0;i<nrParameter;i++){
	  step[i]=theta*step[i]; 
	
	  //	  parameter_new[i]=parameter[i]; //*scaling[i];
	  //	  parameter_new[i]=parameter[i]+(1.0/scaling[i])*step[i].real();
	  //	  parameter_new[i]=(1.0/scaling[i])*parameter_new[i];
	  std::cout<<"paramter_new["<<i<<"]= " << parameter_new[i]<<std::endl;
	  eta = 1- theta*(1-eta);
	}

	updateMaterialData(parameter_new,ptMaterial);
	createF(ptMaterial, ptBCs, F_hat);

	for(Integer i=0;i<F_hat.GetSize();i++)
	  y_hat_F_hat[i]=y_hat[i].real()-F_hat[i].real();
	misfitnew=norm2Real(y_hat_F_hat);
	*piezoLog << "\t misfit after Linesearch " << misfitnew <<std::endl;

	std::cout<<"\n misfitnew = " << misfitnew <<std::endl;
	

      }// end linesearch ...

     
      //choose eta
      eta_new=gamma*POW((misfitnew/misfit),al);
      std::cout<<"\n eta_new =  " << eta_new << std::endl;

      if (gamma*POW(eta,al) > 0.1) // safeguard choice 2, see Pernice, Walker
	eta_new=std::max(eta_new,gamma*POW(eta,al));
      //       std::cout<< "\n eta_new " << eta_new <<std::endl;

      if (eta_new>eta_max)
	eta_new=eta_max;

      if(eta_new<=2*(tau*delta)/misfitnew) // final safeguard
	eta_new=0.8*(tau*delta)/misfitnew;
      //end choose eta
      eta=eta_new;

      if (eta>=eta_max)
	eta=eta_max;

      std::cout<< "\n\n *** Choice of eta = " << eta<<std::endl; 
      *piezoLog<<"\t Choice of eta = " << eta<<std::endl;

      // end choose eta
      parameter=parameter_new;
      misfit=misfitnew;

    }// end while Newton
    std::cout<<"\n leaving Newton CG 2 " <<std::endl;
    std::cout<<"\n\n *** FINALLY CALCULATED PARAMETERS *** ... here they are:"<<std::endl;
    *piezoLog<<"\n\n *** FINALLY CALCULATED PARAMETERS *** ... here they are:"<<std::endl;

    for (Integer  i=0;i<parameter.GetSize();i++){
      std::cout<<"par[" << i<<"]="<< parameter[i]<<";\n";
      *piezoLog<<"parameter["<<i<<"] = "<< parameter[i]<<std::endl;
    }

  }// end NewtonCG 3

  // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


  void piezoParamIdent::NewtonCG2(){
    ENTER_FCN("piezoParamIdent::NewtonCG2");
    std::cout<<"\n Entering piezoParamIdent::NewtonCG2 ... "<< std::endl;

    Integer nrNewtonIterations=0;
    Integer backtrackIterator=0;

    MaterialData * ptMaterial=pdes_[0]->getPDEMaterialData();   // Pointer to MaterialData

    ptBCs = pdes_[0]->getPDE_BCs();

    Double real_misfit, real_misfit_new,eta_max, eta_new, t, aa,b,c, theta_min, theta_max, gamma, al;
    Double res_0, a, a_lin_new, a_new, alpha, beta, tau;
    tau = 1.5;
    Integer i,j;
    Double theta=0.5;

    t=1.0e-4;
    eta_max=0.9;
    eta=0.9;
    theta_min=0.1;
    theta_max=0.5;
    gamma=0.5;
    al=1.5;
    
    
    Vector<Complex> y_hat_F_hat (y_hat.GetSize());

    //    real_misfit=norm_res;
    //real_misfit_new=norm_res_new;

    Vector<Double> step (nrParameter);
    Vector<Double> parstart (nrParameter);
    Vector<Complex> temp_res_NE(nrParameter);
    Vector<Double> temp1(nrParameter);
    Vector<Double> temp2(nrParameter);
    Vector<Double> res_NE_old(nrParameter);
    Vector<Double> res_NE_new(nrParameter);
    Vector<Double> bas_old(nrParameter);
    Vector<Double> bas_new(nrParameter);
    // Vector<Complex> res_old(nrMeasuredData);
    Vector<Complex> res(nrMeasuredData);
    Vector<Complex> res_NE_new_compl(nrParameter);
    Vector<Complex> new_res(nrMeasuredData);
    Complex new_res_res0;
    Vector<Complex> bas_old_compl(nrParameter);

    parstart=parameter;

    // Create the Matrices F, F', F*
    createF(ptMaterial, ptBCs, F_hat);

    for (i=0; i<nrMeasuredData;i++)
      y_hat_F_hat[i]=y_hat[i].real()-F_hat[i].real();

    res=y_hat-F_hat;
      
    a=sqrt(realA2norm(y_hat_F_hat)); 

    std::cout << "\n a = |y-F| =  " <<a <<std::endl;

    while((a>tau*delta*delta||nrNewtonIterations<5)&&nrNewtonIterations<10){ // Newton

      nrNewtonIterations++;
      std::cout<<"\n Newton-Iteration "<< nrNewtonIterations <<std::endl;
      step.Resize(nrParameter);

      //            createJacobiMatrix(ptMaterial, ptBCs, F_hat, parameterIncrement,JacobiMatrix, solElecPot, solMechDispl);
      createJacobiMatrix2(JacobiMatrix);

      createAdjointJacobiMatrix(JacobiMatrix,adjJacobiMatrix);
      //      testJacobiMatrix(F_hat, JacobiMatrix, parameter, ptBCs, ptMaterial,parameterIncrement, solElecPot, solMechDispl);

      adjJacobiMatrix.Mult(res, temp_res_NE);
      
      std::cout<<"\n bas_old = resNE= scal^-2 * F'* res "<<std::endl;

      for(i=0;i<nrParameter;i++){
	//bas_old[i] = res_NE_old[i] = temp_res_NE[i].real();
       bas_old[i] = res_NE_old[i] = 1.0/(scaling[i]*scaling[i])*temp_res_NE[i].real();
       std::cout<<"bas_old = " << bas_old[i]<<std::endl;
      }

      a_lin_new=a;

       std::cout<<"\n a_lin_new  = " << a <<std::endl;

      Integer nrCGIter=0;

      while((a_lin_new>eta*eta*a)&&(nrCGIter<10)){ // CG

	std::cout<<"\n CG Iteration " << nrCGIter << std::endl;
	nrCGIter++;

	JacobiMatrix.MatVecMult_CD(bas_old,bas_bar);

	for (i=0;i<bas_bar.GetSize();i++)
	 std::cout<<"bas_bar="<<bas_bar[i]<<"; ";
	 std::cout<<"\n temp_res_NE:"<<std::endl;
	  
	temp_res_NE.Resize(res_NE_old.GetSize());

	for (i=0;i<res_NE_old.GetSize();i++){
	  temp_res_NE[i]=res_NE_old[i]*scaling[i];
	  std::cout<<"temp_res_NE="<<temp_res_NE[i]<<"; ";
	}

	std::cout<<"\n"<<std::endl;

	Double norm_temp_res_NE = POW(a2norm(temp_res_NE),2);

	Double norm_bas_bar = realA2norm(bas_bar);

	alpha=norm_temp_res_NE / norm_bas_bar;

	std::cout<<"\n alpha = " << alpha << ",   || temp_res_NE || = " << norm_temp_res_NE << ",    || bas_bar || = " << norm_bas_bar << std::endl;

	std::cout<<"\n step = step + alpha *bas: "<<std::endl;

	for (i=0;i<nrParameter;i++){
	  step[i]+=alpha*bas_old[i];
	  std::cout<<"step("<<i<<")= "<< step[i] << "; \t ";
	}
	std::cout<<"\n res = res - alpha*bas_bar "<<std::endl;

	for (i=0;i<res.GetSize();i++){
	  res[i]=res[i]-alpha*bas_bar[i];
	  std::cout<<"res = "<<res[i]<<"; ";
	}
	std::cout<<"\n"<<std::endl;

	//	std::cout<<"\n res.GetSize() " << res.GetSize() << "; res_NE_NEW_COMP.GetSize()= "<< res_NE_new_compl.GetSize()<<std::endl;

	adjJacobiMatrix.Mult(res,res_NE_new_compl);

	for(i=0;i<res_NE_new.GetSize();i++)
	  res_NE_new[i]=1.0/(scaling[i]*scaling[i])*res_NE_new_compl[i].real();

	a_lin_new=POW(realA2norm(res),0.5);

	std::cout<<"\n a_lin_new =  " << a_lin_new<<std::endl; 

	temp1.Resize(nrParameter);
	temp2.Resize(nrParameter);
	for (i=0;i<nrParameter;i++){
	  temp1[i]=scaling[i]*res_NE_new[i];
	  temp2[i]=scaling[i]*res_NE_old[i];
	}
	beta = (POW(a2norm(temp1),2))/(POW(a2norm(temp2),2));
	std::cout<<"\n beta = " << beta <<std::endl; 

	std::cout<<"\n bas = res + beta * bas"<<std::endl;	

	for(i=0;i<nrParameter;i++){
	  bas_new[i]= res_NE_new[i]+beta*bas_old[i];
	  std::cout<<"bas = "<<bas_new [i]<<"; ";
	}


	bas_old=bas_new;
	res_NE_old=res_NE_new;

      }// end while CG


      std::cout<<"\n"<<std::endl;

      for (i=0;i<nrParameter;i++){
	std::cout<<"step("<<i<<")= " << step[i]<< "; \t";
	parameter_new[i]=parameter[i]; //*scaling[i];
	parameter_new[i]=parameter_new[i]+step[i];
	//	parameter_new[i]=(1.0/scaling[i])*parameter_new[i];
	  std::cout<<"Paramter_new["<<i<<"]= " << parameter_new[i]<<std::endl;
      }
      std::cout<<"\n"<<std::endl;

      // precautionary measure, if paramters tend to far away from initial values ...
      for(i=0;i<nrParameter;i++)
	if (std::abs(parameter_new[i]-parstart[i])>0.5*std::abs(parstart[i])){
	  parameter_new[i]=parameter[i];
	  std::cout<<"\n parameter("<<i<<") was reset to value " << parameter[i] <<std::endl;
	}

      updateMaterialData(parameter_new, ptMaterial);
      createF(ptMaterial, ptBCs, F_hat);

      for(i=0;i<F_hat.GetSize();i++)
	new_res[i]=y_hat[i]-F_hat[i];


      //	backtracking(eta, theta, s, a, a_lin_new);
      a_new=POW(realA2norm(new_res),0.5);

      while ((a_new>=(1-t*(1-eta))*a) && backtrackIterator<20 && eta<eta_max){ // Liniensuche
	std::cout<<"\n backtracking ... "<< backtrackIterator<< std::endl;
	backtrackIterator++;

	// choose theta:
	aa=a*a;
	new_res_res0=Complex(0.0,0.0);

	for(i=0;i<nrMeasuredData;i++)
	  new_res_res0=new_res_res0+new_res[i]*(res[i].real()-new_res[i]);

	b=2.0*new_res_res0.real();
	std::cout<<"\n b = " << b << std::endl;
	c=a_new*a_new-b-aa;
	std::cout<<"\n c = "<<c<<std::endl;
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
	      if ((aa+b*theta_min+c*theta_min*theta_min)>=(aa+b*theta_max+c*theta_max*theta_max))
		theta=theta_max;
	      else theta=theta_min;
	    }
	} // end else if c>0

	std::cout<<"\n Choice of theta = " << theta<<std::endl;
	std::cout<<"\n"<<std::endl;

	for(i=0;i<nrParameter;i++){
	  step[i]=theta*step[i];
	  parameter_new[i]=parameter[i]; //*scaling[i];
	  parameter_new[i]+=step[i];
	  //	  parameter_new[i]=(1.0/scaling[i])*parameter_new[i];
	  std::cout<<"Paramter_new["<<i<<"]= " << parameter_new[i]<<std::endl;
	  eta = 1- theta*(1-eta);
	}

	updateMaterialData(parameter_new,ptMaterial);
	createF(ptMaterial, ptBCs, F_hat);

	for(i=0;i<F_hat.GetSize();i++)
	  y_hat_F_hat[i]=y_hat[i]-F_hat[i];
	a_new=POW(realA2norm(res),0.5);

	std::cout<<"\n a_new = " << a_new <<std::endl;
	

      }// end linesearch ...

     
      //choose eta
      eta_new=gamma*POW((a_new/a),al);
      std::cout<<"\n eta_new =  " << eta_new << std::endl;

      if (gamma*POW(eta,al) > 0.1) // safeguard choice 2, see Pernice, Walker
	eta_new=std::max(eta_new,gamma*POW(eta,al));
      //       std::cout<< "\n eta_new " << eta_new <<std::endl;

      if (eta_new>eta_max)
	eta_new=eta_max;

      if(eta_new<=2*(tau*delta)/a_new) // final safeguard
	eta_new=0.8*(tau*delta)/a_new;
      //end choose eta
      eta=eta_new;
      if (eta>=eta_max)
	eta=eta_max;

      std::cout<< "\n\n *** Choice of eta = " << eta<<std::endl; 

      // end choose eta
      parameter=parameter_new;

    }// end while Newton
    std::cout<<"\n leaving Newton CG 2 " <<std::endl;
 std::cout<<"\n\n *** FINALLY CALCULATED PARAMETERS *** ... here they are:"<<std::endl;

    for (int i=0;i<parameter.GetSize();i++)
      std::cout<<"par[" << i<<"]="<< parameter[i]<<";\n";

  }// end NewtonCG 2




  void piezoParamIdent::NewtonCG(){

    ENTER_FCN("piezoParamIdent::Newton()");

    // Settings for Newton-CG - routine
    Integer nrIterations=0;
    Double tol=0.001, tolCG=0.001, err=1.0, res_norm;


    MaterialData * ptMaterial=pdes_[0]->getPDEMaterialData();   // Pointer to MaterialData
    Matrix<Double> * matMatrix =  ptMaterial->GetMatrix();

    updateMaterialData(parameter, ptMaterial);         //Writes initial guesses of parameters (read from MeasuredData.dat) to system

    ptBCs = pdes_[0]->getPDE_BCs();                             // Pointer to BCs
    ptAlgsys = pdes_[0]->getPDE_algsys();                       //Pointer to AlgebraicSystem
    Integer numElems = pdes_[0]->getPDE_numElems();
    dofs=pdes_[0]->getPDE_dofspernode();
    numNodes= pdes_[0]->getPDE_numPDENodes();

    while (nrIterations<5) {
      nrIterations++;
      std::cout<<"\n ******************** \n NewtonCG ... Newton-Iteration-Nr = "<<nrIterations<<std::endl;

      // Create the Matrices F, F', F*
      createF(ptMaterial, ptBCs, F_hat);
      //std::cout<<"Parameter-to-solution-map Matrix is built up!"<<std::endl;
      createJacobiMatrix(ptMaterial, ptBCs, F_hat, parameterIncrement,JacobiMatrix, solElecPot, solMechDispl);
      // std::cout<<"JacobiMatrix was created ..."<<std::endl;
      createAdjointJacobiMatrix(JacobiMatrix,adjJacobiMatrix);
      // std::cout<<"Adjoint JacMatrix was created ..."<<std::endl;
       
      //std::cout<<"\n\ntest Jacobi Matrix ..."<<std::endl;
      //  testJacobiMatrix(F_hat, JacobiMatrix, parameter, ptBCs, ptMaterial,parameterIncrement, solElecPot, solMechDispl);

      Double real_misfit, real_misfit_new,eta_max, eta, theta;
      Complex misfit;
      Vector<Double> y_hat_F_hat;
      for (Integer i=0;i<y_hat.GetSize();i++)
	y_hat[i].real()-F_hat[i].real();

      real_misfit=POW(a2norm(y_hat_F_hat),2);
      std::cout<<"\n Vor CG real(misfit) = real(y_hat-F_hat) = "<< real_misfit<<std::endl;

      // we do a bit of scaling here ...
      //      for(int i=0;i<scaling.GetSize();i++)
      //	parameter[i]=parameter[i]*scaling[i];

      // xxxxxxxxxxxxxxxx how should we solve the inner problem? xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

      CG(); // calls CG Routine ... uses F, F', F'*    

      std::cout<< "\n s determined with CG ... " <<std::endl;
      for (int i=0;i<parameter.GetSize();i++){ 
	std::cout<<"s("<<i<<")="<<s[i]<<"; \t";
       	parameter_new[i]=parameter[i]+s[i].real();        
      }
      std::cout<<"\n"<<std::endl;

      // Matrix<Complex> JacobiMatrixNE(JacobiMatrix.GetSizeCol(), JacobiMatrix.GetSizeCol());
//       Vector<Complex> y_hat_F_hatNE(JacobiMatrix.GetSizeCol());

//       adjJacobiMatrix.Mult(JacobiMatrix,JacobiMatrixNE);
//       adjJacobiMatrix.Mult(y_hat_F_hat, y_hat_F_hatNE);


//         directSolve(JacobiMatrixNE, s , y_hat_F_hatNE);
      
//        std::cout<< "\n s determined with direct solver ... " <<std::endl;
//              for (int i=0;i<parameter.GetSize();i++){ 
//        	std::cout<<"s("<<i<<")="<<s[i]<<"; \t";
//        	parameter[i]=parameter[i]+s[i].real();        
//        	std::cout<<"par_new("<<i<<")="<<parameter_new[i]<<"; \t";
//              }
//              std::cout<<"\n"<<std::endl;

      // xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


      // we do a bit of rescaling here ...
      //      for(int i=0;i<scaling.GetSize();i++)
      //	parameter[i]=parameter[i]*1.0/scaling[i];

      updateMaterialData(parameter_new, ptMaterial); 

      createF(ptMaterial, ptBCs, F_hat);

      // we do a bit of scaling here ...
      //        for(int i=0;i<scaling.GetSize();i++)
      //parameter_new[i]=parameter_new[i]*scaling[i];

      for(Integer i=0;i<y_hat.GetSize();i++)
	y_hat_F_hat=y_hat[i].real()-F_hat[i].real();
      real_misfit_new=a2norm(y_hat_F_hat);
  
      //here begins attempt to implement the backtracking Algo .... ;-)

      backtracking(eta, theta, s, real_misfit, real_misfit_new);
     
      std::cout<<"\n\n currently calculated parameters: " << std::endl;
      
      // we do a bit of rescaling here ...
      //     for(int i=0;i<scaling.GetSize();i++)
      //	parameter[i]=parameter[i]/scaling[i];
      //parameter=parameter_new;

      for (int i=0;i<parameter.GetSize();i++){
	std::cout<<"par[" << i<<"]="<< parameter[i]<<"; \t";
	std::cout<<"\n"<<std::endl;
      }
      ///	parameter[i]=parameter[i]+s[i].real();
      	
      updateMaterialData(parameter,ptMaterial);

    } // end while-Newton


  } // end Newton


}// end namespace
