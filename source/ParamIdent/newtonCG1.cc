// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <iostream>
#include <fstream>
#include <string>
//#include "staticdriver.hh"
#include "DataInOut/GMV/outGMV.hh"
#include "CoupledPDE/basecoupledpde.hh"
#include "General/environment.hh"
#include "PDE/basePDE.hh" 

#include "piezoParamIdent.hh"
#include "Forms/baseForm.hh"
#include "MatVec/vector.hh"
#include "Utils/nodestoresol.hh"
#include "Utils/elemstoresol.hh"
#include "DataInOut/MaterialData.hh"
#include "PDE/timestepping.hh"
#include "Utils/baseelemstoresol.hh"
#include "Driver/singleDriver.hh"
#include "PDE/nodeEQN.hh"
#include <Domain/elem.hh>


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



//#include "/../OLAS/algsys/algebraicSys.hh"
//#include "DataInOut/piezoParameterData.hh"



namespace CoupledField
{

  // xxxxxxxxxxxxxxxxxxxxxxxxxxxxx NEWTON CG 4 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  void piezoParamIdent::NewtonCG1(){   // here we scale in advance, just once before NewtonCG and afterwards

    Integer nrNewtonIterations=0;
    Integer backtrackIterator=0;

    MaterialData * ptMaterial=pdes_[0]->getPDEMaterialData();   // Pointer to MaterialData

    Double eta_max, eta_new, t, aa,b,c, theta_min, theta_max, gamma, al;
    Double alpha, beta, tau;
    tau = 1.1;
    Integer i,j;
    Double theta=0.5;
 
    t=1.0e-4;
    eta_max=0.9;
    eta=0.9;
    theta_min=0.1;
    theta_max=0.5;
    gamma=0.5;
    al=1.5;
    
    Vector<Complex> res(nrMeasuredData);
    Vector<Complex> step(actNrParameter+actNrParameterC);
    Vector<Double> stepR(actNrParameter);    
    res_NE.Resize(actNrParameter);
    res_NE.Init();
    Vector<Complex> res_NE_rescaled(actNrParameter);
    bas.Resize(actNrParameter);
    bas.Init();
    Vector<Complex> basbar(nrMeasuredData);
    Vector<Complex> Res_outer(nrMeasuredData);
    Vector<Complex> Res_inner(nrMeasuredData);
    Vector<Complex> JacFs(nrMeasuredData);
    Vector<Complex> res_JacFs(nrMeasuredData);
    Vector<Complex> basC(actNrParameter);
    res_NE.Resize(actNrParameter);
    res_NE.Init();

    Double res_inner, res_inner_old, res_outer, res_outer_old, lin_res, lin_res_max, max_res, max_res_outer, max_res_inner;

    Double norm_bas_bar, norm_resNE, norm_resNEOld;

    // NO scaling in NewtonCG!! 
    for (Integer i=0;i<parameter.GetSize();i++){
      scaling[i]=1.0;
    }

    for (i=0;i<actNrParameter;i++){
      bas[i]=1.0;
    }

    
    // NEWTON ITERATION -- outer Loop!!
    while(nrNewtonIterations<maxNumberNewtonLoops){ // Newton
      *piezoLog << "\n Newton-Iteration: " << nrNewtonIterations <<std::endl;
      *piezoLog <<"------------------------"<<std::endl;
    
      nrNewtonIterations++;

      step.Resize(actNrParameter);
      step.Init();

      updateMaterialData(parameter, ptMaterial);

      createF(ptMaterial, ptBCs, F_hat,false);

      res=y_hat-F_hat;
      Res_inner=res;

      for(Integer i=0;i<actNrParameter;i++){
        bas[i]= basC[i].real();
      }
       
      // Create the matrices ...
      createJacobiMatrix2(JacobiMatrix);

      createAdjointJacobiMatrix(JacobiMatrix,adjJacobiMatrix);
  
      for (Integer i=0;i<parameter.GetSize();i++){
        scaling[i]=1.0;
      }

      adjJacobiMatrix.Mult(res, res_NE);

      //       for (Integer i=0;i<actNrParameter;i++){
      //      res_NE[i]=(1.0/scaling[whichParToUpInd[i]]*scaling[whichParToUpInd[i]])*res_NE[i];
      //      res_NE_rescaled[i]=scaling[whichParToUpInd[i]]*res_NE[i];
      //        }

      basC = res_NE;

      // ----- INNER LOOP --- CG-ITERATION ---
      Integer nrCGIter=0;

      res_outer_old=res_outer;
      res_inner=res_outer;
      
      while(nrCGIter<maxNumberInnerLoops){ // CG
        nrCGIter++;

        JacobiMatrix.Mult(basC,basbar);

        norm_bas_bar = POW(a2norm(basbar),2);

        norm_resNE = POW(a2norm(res_NE),2);

        res_inner_old=res_inner;

        norm_resNEOld = norm_resNE;

        alpha=norm_resNE / norm_bas_bar;
        
        //      std::cout<<"\n alpha = " << alpha << ",\t || normresNEold || = " << norm_resNEOld <<
        // ",\t ||bas_bar|| = " << norm_bas_bar << std::endl;


        //      std::cout<<"\nstep = step + alpha *bas: "<<std::endl;
        for (Integer i=0;i<actNrParameter;i++){
          step[i]+=Complex(alpha,alpha)*basC[i];
        }
        std::cout<<"\n basC:"<<std::endl;
        std::cout<<basC<<std::endl;

        //      std::cout<<"\n\nres = res - alpha*bas_bar: "<<std::endl;
        for (Integer i=0;i<nrMeasuredData;i++){
          Res_inner[i]=Res_inner[i]-Complex(alpha,alpha)*basbar[i];
        }

        norm(Res_inner,res_inner,max_res_inner,y_hat);

        adjJacobiMatrix.Mult(res,res_NE);

        norm_resNE=POW(a2norm(res_NE),2);

        beta = norm_resNE/norm_resNEOld;

        std::cout<<"\n beta = " << beta << std::endl;


        //      std::cout<<"\nbas = resNE + beta * bas: "<<std::endl;   
        for(Integer i=0;i<actNrParameter;i++){
          basC[i]= res_NE[i]+beta*basC[i];
        }

        JacobiMatrix.Mult(step,JacFs);

        res_JacFs = res-JacFs;

        std::cout<<res_JacFs<<std::endl;

        norm(res_JacFs, lin_res, lin_res_max ,y_hat);

        std::cout<<"\n CG-Iter: " << nrCGIter << " || res_JacFs || = " << lin_res_max << ", max_res_outer = " << max_res_outer << std::endl;

        *piezoLog<<"\n CG-Iter: " << nrCGIter << " || res_JacFs || = " << lin_res_max << ", Res_outer = " << res_outer << std::endl;

        eta=1.0;
        if (lin_res_max <= 0.1*max_res_outer){
          *piezoLog<<"\n Terminated CG-Iteration after " << nrCGIter << " steps with res_JacFs = " << lin_res << std::endl;
          std::cout<<"\n Terminated CG-Iteration after " << nrCGIter << " steps with res_JacFs = " << lin_res << std::endl;
          getchar();
          break;
        }
        
      }// end while CG

      parameter_new=parameter; 
     
      Matrix<Double> *matMat = ptMaterial->GetMatrix();
     
      scaling[0]=1.0/((*matMat)[0][0]); 
      scaling[1]=7.5/((*matMat)[2][2]);
      scaling[2]=1.0/((*matMat)[1][0]);
      scaling[3]=1.0/((*matMat)[0][2]);
      scaling[4]=1.0/((*matMat)[3][3]); 
      scaling[5]=1.0/((*matMat)[6][4]);
      scaling[6]=std::abs(1.0/((*matMat)[8][0]));
      scaling[7]=1.0/((*matMat)[8][2]);
      scaling[8]=1.0/((*matMat)[6][6]); 
      scaling[9]=1.0/((*matMat)[8][8]);

      //      parameter_new[9]=parameter[9]+(1.0/scaling[9])*step[9].real();        // eps33      
      for (Integer i=0;i<actNrParameter;i++)
        stepR[i]=step[i].real();
      std::cout<<step<<std::endl;

      setNewParameterSet(parameter, parameter_new, scaling,theta,stepR,whichParameterToUpdate);

      updateMaterialData(parameter_new, ptMaterial);

      createF(ptMaterial, ptBCs, F_hat,false);
      
      res=y_hat-F_hat;

      norm(res,res_outer,max_res_outer,y_hat);

      parameter=parameter_new;

      std::cout<<parameter<<std::endl;
      *parLog <<nrNewtonIterations <<"  "<< parameter[1]<<"  " <<parameter[7]<<"   " <<parameter[9]<<std::endl;

      if(max_res_outer<=5.e-05)
        getchar();

      if (res_outer>=res_outer_old){
        *piezoLog<<"\n Newton Iteration terminated after "<< nrNewtonIterations << " steps with res = " << res_outer << std::endl;
        std::cout<<"\n Newton Iteration terminated after "<< nrNewtonIterations << " steps with res = " << res_outer << std::endl;
        //      getchar();
        //break;
        //      std::cout<<parameter<<std::endl;
        //getchar();
      }
   
      //        backtracking(eta, theta, s, a, a_lin_new);

      backtrackIterator=0;

      while ((res_outer>POW((1.0-t*(1.0-eta)),2)*res_outer_old) && backtrackIterator<0){ // Liniensuche
        //      std::cout<<"\n backtracking ... "<< backtrackIterator<< std::endl;
        backtrackIterator++;

        b=0.0;
        for(Integer i=0;i<nrMeasuredData;i++)
          b+= 2.0*(res[i]*(res_JacFs[i]-res[i])).real();

        // choose theta:
        aa=res_outer_old;

        c=res_outer-b-res_outer_old;

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
        //      theta = 0.35;

        *piezoLog<<"\t Choice of theta = " << theta<<std::endl;

        // for (Integer i=0;i<actNrParameter;i++)
        //        step[i]=theta*step[i];

        for (Integer i=0;i<actNrParameter;i++)
          stepR[i]=theta*stepR[i];

        setNewParameterSet(parameter, parameter_new, scaling,theta,stepR,whichParameterToUpdate);

        for(Integer i=0;i<nrParameter;i++)
          std::cout<<"paramter_new["<<i<<"]= " << parameter_new[i]<<" + " << parameter_newC[i]<<std::endl;

        //      b = theta*b;
        eta = 1- theta*(1-eta);

        updateMaterialData(parameter_new,ptMaterial);

        createF(ptMaterial, ptBCs, F_hat,false);

        res=y_hat-F_hat;

        norm(res,res_outer,max_res_outer,y_hat);

        std::cout<<"\n New res_outer = " << res_outer <<std::endl;
        //      getchar();

        //      *piezoLog << "\t aval_new after Linesearch " << aval_new <<std::endl;

      }// end linesearch ...

      //choose eta

      eta_new = gamma*POW((res_outer/res_outer_old),al);
      // std::cout<<"\n eta_new =  " << eta_new << std::endl;
      Double gammaEtaAl = gamma*POW(eta,al);
      if (gammaEtaAl > 0.1) // safeguard choice 2, see Pernice, Walker
        if (gammaEtaAl>eta)
          eta_new=gammaEtaAl;

      if(eta_new<=2*(tau*delta)/sqrt(res_outer)); // final safeguard
      eta_new=0.8*(tau*delta)/sqrt(res_outer);
      //end choose eta

      if (eta_new>eta_max)
        eta=eta_max;
      else
        eta=eta_new;

      //    std::cout<< "\n\n *** Choice of eta = " << eta<<std::endl; 
      //      *piezoLog<<"\t Choice of eta = " << eta<<std::endl;

      // end choose eta
      parameter=parameter_new;
   
      res_outer_old=res_outer;

    }// end while Newton
    //   std::cout<<"\n leaving Newton CG 4 " <<std::endl;
    //std::cout<<"\n\n *** FINALLY CALCULATED PARAMETERS *** ... here they are:"<<std::endl;
    //    *piezoLog<<"\n\n *** FINALLY CALCULATED PARAMETERS *** ... here they are:"<<std::endl;

    //        for (Integer  i=0;i<parameter.GetSize();i++){
    //std::cout<<"par[" << i<<"]="<< parameter[i]<<";\n";
    //    *piezoLog<<"parameter["<<i<<"] = "<< parameter[i]<<std::endl;
    //}

    //    *piezoLog <<"MaxResiduum: " << max_res_outer <<", weighted norm = " <<aval_new<< std::endl;
    finalnorm=max_res_outer;

  }// end NewtonCG 4
  

  // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

} // end namespace ...
