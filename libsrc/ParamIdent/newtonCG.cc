#include "PDE/SinglePDE.hh" 
#include "piezoParamIdent.hh"


namespace CoupledField
{

  
  // xxxxxxxxxxxxxxxxxxxxxxxxxxxxx NEWTON CG 3 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  void piezoParamIdent::NewtonCG3(){   // here we scale in advance, just once before NewtonCG and afterwards
    ENTER_FCN("piezoParamIdent::NewtonCG3");
    std::cout<<"\n Entering piezoParamIdent::NewtonCG3 ... "<< std::endl;

    UInt nrNewtonIterations=0;
    UInt backtrackIterator=0;

    MaterialData * ptMaterial=ptMyPDE_->getPDEMaterialData();   // Pointer to MaterialData

    Double eta_max, eta_new, t, theta_min, theta_max, gamma, al, eta_min;
    Double alpha, beta, tau;
    tau = 1.1;
    UInt i;
    Double theta=1.0;
    //delta = 0.1;

    t=1.0e-4;
    eta_max=0.9;
    eta=0.9;
    eta_min=0.1;
    theta_min=0.1;
    theta_max=1.0;
    gamma=0.5;
    al=1.5;
    

    Vector<Complex> y_hat_F_hat (nrMeasuredData);

    Vector<Complex> step(actNrParameter);
    Vector<Double> stepR(actNrParameter);
    Vector<Double> parstart (nrParameter);
    
    Vector<Double> res_NE(actNrParameter);
    Vector<Double> res_NE_rescaled(actNrParameter);
    Vector<Complex> res_NE_C(actNrParameter);
    bas.Resize(actNrParameter);
    Vector<Complex> basbar(nrMeasuredData);
    Vector<Complex> basbar_old(nrMeasuredData);
    Vector<Complex> res(nrMeasuredData);
    Vector<Complex> JMats(nrMeasuredData);
    Vector<Complex> res0(nrMeasuredData);
    Vector<Complex> Res_outer(nrMeasuredData);
    Vector<Complex> Res_inner(nrMeasuredData);
    Vector<Complex> Res_linear(nrMeasuredData);
    Vector<Complex> JacFs_res(nrMeasuredData);
    Vector<Complex> JacFs(nrMeasuredData);


    Complex new_res_res0;
    Double misfit, misfit0, misfitnew, normres02, normresNEold2, normresNE2;

    Double new_res_inner;

    Double aval, aval_new, alin_new, maxres_outer, maxres_inner, aval_old, res_linear, res_linear_max;


    // NO scaling in NewtonCG!! 
    for (UInt i=0;i<parameter.GetSize();i++)
      scaling[i]=1.0;
    //    parameter[i]=scaling[i]*parameter[i];

    parstart=parameter;

    updateMaterialData(parameter,ptMaterial);
    createF(ptMaterial, F_hat, FALSE);

    std::cout<<"F_hat:"<<std::endl;
    std::cout<<F_hat<<std::endl;

    std::cout<<"y_hat:"<<std::endl;
    std::cout<<y_hat<<std::endl;

    for (i=0; i<nrMeasuredData;i++)
      y_hat_F_hat[i]=y_hat[i]-F_hat[i];

    for (i=0;i<actNrParameter;i++)
      bas[i]=1.0;

    Res_outer=res=y_hat-F_hat;
      
    //misfit=norm2Real(y_hat_F_hat); 

    norm(Res_outer,aval,maxres_outer,y_hat);

    std::cout<<"\n Newton Step 0: Maxres = " << maxres_outer << ", aval = "<< aval<<std::endl;

    aval_new=aval;
    aval_old=aval;

    misfit0=misfit;


    *piezoLog << "\n Initial misfit ||y-F|| = "<<aval<< " and maxres_outer = max |y-F| = " << maxres_outer << std::endl;

    nrNewtonIterations=0;

    // NEWTON ITERATION -- outer Loop!!
    while(aval>tau*delta&&nrNewtonIterations<maxNumberNewtonLoops){ // Newton
      *piezoLog << "\n Newton-Iteration: " << nrNewtonIterations <<std::endl;
      *piezoLog <<"------------------------"<<std::endl;
      std::cout<<"\n Newton-Iteration "<< nrNewtonIterations <<std::endl;

      nrNewtonIterations++;

      res0=res;
      normres02=norm2Real(res0);
      step.Resize(actNrParameter);    
    
      // Create the matrices ...
      updateMaterialData(parameter,ptMaterial);      
      createF(ptMaterial, F_hat,TRUE);
      //     getchar();
      //     getchar();
      //     getchar();
      //     getchar();

      std::cout<<parameter<<std::endl;
      // getchar();
      //      createJacobiMatrix2(JacobiMatrix);
      testJacobiMatrix2(F_hat, JacobiMatrix, parameter, ptMaterial,parameterIncrement, solElecPot, solMechDispl);

      //std::cout<<JacobiMatrix<<std::endl;

      //      std::cout<<approxJacobiMatrix<<std::endl;
      //getchar();

      JacobiMatrix = approxJacobiMatrix;
      createAdjointJacobiMatrix(JacobiMatrix,adjJacobiMatrix);
           
      for (UInt i=0;i<parameter.GetSize();i++)
        scaling[i]=1.0;

      //JacobiMatrix.MatVecMult_CD(step,JMats); // is not neccessary, since we start with s_i=0;
      //res = res0-JMats;
      
      //       normres2 = norm2Real(res);

      adjJacobiMatrix.Mult(res, res_NE_C);

      for (UInt i=0;i<res_NE_C.GetSize();i++){
        res_NE[i]=(1.0/scaling[whichParToUpInd[i]]*scaling[whichParToUpInd[i]])*res_NE_C[i].real();
        res_NE_rescaled[i]=scaling[whichParToUpInd[i]]*res_NE[i];
        // std::cout<<"res_NE="<<res_NE[i]<<"; ";
      }

      //  for (UInt i=0;i<res_NE_C.GetSize();i++){
      //      res_NE[i]=res_NE_C[i].real();
      //      res_NE_rescaled[i]=res_NE[i];
      //      // std::cout<<"res_NE="<<res_NE[i]<<"; ";
      //       }

      bas = res_NE;
      
      alin_new = aval;
      *piezoLog << "\n \t alin_new = aval = " << alin_new <<"(e_33, eps_33 = "<<parameter[7]<<", " <<parameter[9]<<")";

      normresNEold2 = POW(a2norm(res_NE_rescaled),2);
      basbar_old.Resize(actNrParameter);

      Res_inner = Res_outer;


      // ----- INNER LOOP --- CG-ITERATION ---
      UInt nrCGIter=0;
      //      eta=1.0;
      ///      while((alin_new<=eta*eta*aval)&&(nrCGIter<maxNumberInnerLoops)||nrCGIter<1){ // CG
        while(nrCGIter<maxNumberInnerLoops){ // CG
          //      while((normres2>eta*eta*normres02)&&(nrCGIter<9)){ // CG
          //      std::cout<<"\n CG Iteration " << nrCGIter << std::endl;
          nrCGIter++;

        
          JacobiMatrix.MatVecMult_CD(bas,basbar);
        
          //      std::cout<<"\nbasbar = F' bas:  "<<std::endl;
          //      for (UInt i=0;i<basbar.GetSize();i++){
          //        basbar[i]=basbar[i].real();
          //std::cout<<"basbar="<<basbar[i]<<"; ";
          //      }
          //std::cout<<"\n "<<std::endl;    

          basbar_old=basbar;

          Double norm_bas_bar = POW(a2norm(basbar),2);

          alpha=normresNEold2 / norm_bas_bar;

          std::cout<<"\n alpha = " << alpha << ",\t || normresNEold || = " << normresNEold2 << ",\t ||bas_bar|| = " << norm_bas_bar << std::endl;

          //      getchar();
          //      std::cout<<"\nstep = step + alpha *bas: "<<std::endl;
          for (UInt i=0;i<actNrParameter;i++){
            step[i]+=1.0*alpha*bas[i];
            //      std::cout<<"step("<<i<<")= "<< step[i] << "; \t ";
          }
          std::cout<<step<<std::endl;

          JacobiMatrix.Mult(step,JacFs);
        
          //F'(p^k)(s^k)-(y-F(p^k))
          for (UInt i=0;i<nrMeasuredData;i++){
            //res[i]=y_hat[i]-F_hat[i];
            JacFs_res[i]=JacFs[i]-res[i];
          }

          //      for (UInt i=0;i<nrMeasuredData;i++)
          //act_res[i]-=JacFs[i];
        
          new_res_inner=a2norm(JacFs_res); 
          //      std::cout<<"\n new_res_inner = "<< new_res_inner <<", old_res_inner = " << old_res_inner<< std::endl;
          std::cout<<"\n new_res_inner = "<< new_res_inner <<", normres02 = " << normres02 << std::endl;

        
//           if (new_res_inner>0.00010*normres02){
//             std::cout << " \n !! New_res_inner is worse than old_res_inner -> break of inner Loop! "<< std::endl;
//             //std::cout<<"\n Nr of CG Iterations = " << nrCGIter <<std::endl;
//             //            getchar();
//             //      break;
//           }

          //      std::cout<<"\n\nres = res - alpha*bas_bar: "<<std::endl;
          for (UInt i=0;i<res.GetSize();i++){
            Res_inner[i]=Res_inner[i]-alpha*basbar[i];
            res[i]=res[i]-alpha*basbar[i];
            //      res[i]=res[i].real();
            //      std::cout<<"Res_inner = "<<res[i]<<"; ";
          }

          //std::cout<<"\n"<<std::endl;

          norm(Res_inner,alin_new,maxres_inner,y_hat);
          //      maxAndEuclNorm(Res_inner,maxres_inner,alin_new);
          //              maxAndWeightedResNorm(Res_inner, maxres_inner, alin_new, y_hat);

          //      std::cout<<"\n res.GetSize() " << res.GetSize() << "; res_NE_NEW_COMP.GetSize()= "<< res_NE_new_compl.GetSize()<<std::endl;

          adjJacobiMatrix.Mult(res,res_NE_C);

          //      std::cout<<"\nres_NE = (F*' res).real():"<<std::endl;
          for(UInt i=0;i<res_NE_new.GetSize();i++){
            res_NE[i]=(1.0/(scaling[whichParToUpInd[i]]*scaling[whichParToUpInd[i]]))*res_NE_C[i].real();
            res_NE_rescaled[i]=1.0/(scaling[whichParToUpInd[i]])*res_NE_C[i].real();
            //res_NE[i]*scaling[whichParToUpInd[i]];
            //      res_NE_new[i]=1.0/(scaling[i]*scaling[i])*res_NE_new_compl[i].real();
            //      std::cout<<"resNE = " << res_NE[i]<<std::endl;
          }

          normresNE2=POW(a2norm(res_NE_rescaled),2);

          //      std::cout<<"\n\n normresNE2 =  " << normresNE2<<std::endl; 

          beta = normresNE2/normresNEold2;
        
          JacobiMatrix.Mult(s,JacFs);
          Res_linear=res-JacFs;
          norm(Res_linear,res_linear,res_linear_max,y_hat);
        

//           if (res_linear<alin_new){
//             std::cout<<"\n res_linear = " << res_linear << " < = alin_new = " << alin_new << std::endl;
//             getchar();
//             break;
//           }

          normresNEold2 = normresNE2;

          //      std::cout<<"\nbas = resNE + beta * bas: "<<std::endl;   
          for(UInt i=0;i<actNrParameter;i++){
            bas[i]= res_NE[i]+beta*bas[i];
            //std::cout<<"bas = "<<bas[i]<<"; ";
          }

          //      normres2=norm2Real(res);
          aval=alin_new;
          norm(res,alin_new,maxres_inner,y_hat);
          //maxAndEuclNorm(res,maxres_inner, alin_new);
          //maxAndWeightedResNorm(res, maxres_inner, alin_new, y_hat);

          *piezoLog << "\n \t \t alin_new  = " << alin_new << "(e_33, eps_33 = "<<parameter[7]<<", " <<parameter[9]<<")";

          //      *piezoLog <<"\t\tnormres2 = ||res|| = "<<normres2 << " normres02 = "<< normres02 << " \nand maxres_inner = " << maxres_inner<< std::endl;

        }// end while CG

        *piezoLog << "\n \t Number of CG - Iterations performed " << nrCGIter <<std::endl;

        parameter_new=parameter; // no update for other parameters

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

        for(UInt i=0;i<actNrParameter;i++)
          stepR[i]=step[i].real();  
        std::cout<<stepR<<std::endl;
        getchar();
     
        //      parameter_new[9]=parameter[9]+(1.0/scaling[9])*step[9].real();        // eps33      
        theta=1.0;
        setNewParameterSet(parameter, parameter_new, scaling, theta, stepR, whichParameterToUpdate);

        //       for (UInt i=0;i<nrParameter;i++){
        //        //      std::cout<<"step("<<i<<")= " << step[i]<< "; \t";
        //        //      //      parameter_new[i]=parameter[i]; //*scaling[i];
        //        //      //      parameter_new[i]=parameter[i]+step[i];
        //        //      //      parameter_new[i]=(1.0/scaling[i])*parameter_new[i];
        //        //      std::cout<<"Paramter_new["<<i<<"]= " << parameter_new[i]<<std::endl; 
        //       }
        //       //      std::cout<<"\n"<<std::endl;

        //       // precautionary measure, if paramters tend to far away from initial values ...
        //       //       for(UInt i=0;i<nrParameter;i++)
        //       //       if (std::abs(parameter_new[i]-parstart[i])>0.5*std::abs(parstart[i])){
        //       //         parameter_new[i]=parameter[i];
        //       //         std::cout<<"\n parameter("<<i<<") was reset to value " << parameter[i] <<std::endl;
        //       //       }

        std::cout<<"\n Parameter directly after CG ..." <<std::endl;
        std::cout<<parameter_new<<std::endl;
        parameter=parameter_new;

        updateMaterialData(parameter_new, ptMaterial);
        createF(ptMaterial, F_hat,FALSE);

        for(UInt i=0;i<F_hat.GetSize();i++)
          y_hat_F_hat[i]=y_hat[i]-F_hat[i];

        Res_outer = y_hat-F_hat;


        norm(Res_outer,aval_new,maxres_outer,y_hat);

        std::cout<<"\n ||res|| = " << Res_outer << std::endl;
        //      maxAndWeightedResNorm(Res_outer, maxres_outer, aval_new, y_hat);
        //      maxAndEuclNorm(Res_outer, maxres_outer, aval_new);

        std::cout<< parameter<<std::endl;

//         if(maxres_outer<=1.0e-8)
//           getchar();


        //        backtracking(eta, theta, s, a, a_lin_new);
        misfitnew=norm2Real(y_hat_F_hat);
        *piezoLog << "\t misfit before Linesearch " << aval_new << "(e_33, eps_33 = "<<parameter[7]<<", "<<parameter_new[9]<<")"<<std::endl;

        backtrackIterator=0;

     
        //       //      while ((misfitnew>(1-t*(1-eta))*misfit) && backtrackIterator<15 && eta<eta_max || backtrackIterator<2){ // Liniensuche
        //       while ((aval_new>POW((1.0-t*(1.0-eta)),2)*aval) && backtrackIterator<0 && eta<eta_max || backtrackIterator<0){ // Liniensuche
        //      std::cout<<"\n backtracking ... "<< backtrackIterator<< std::endl;
        //      backtrackIterator++;

        //      b=0.0;
        //      for(UInt i=0;i<nrMeasuredData;i++)
        //        b+= 2.0*(Res_outer[i]*(Res_inner[i]-Res_outer[i])).real();


        //      // choose theta:
        //      aa=aval;
        //      new_res_res0=Complex(0.0,0.0);

        //      //        new_res_res0=new_res_res0+y_hat_F_hat[i]*(res[i]-y_hat_F_hat[i]);
        //      //      b=2.0*new_res_res0.real();

        //      //      std::cout<<"\n b = " << b << std::endl;

        //      //      c=misfitnew*misfitnew-b-aa;
        //      c=aval_new-b-aval;
        //      //      std::cout<<"\n c = "<<c<<std::endl;

        //      if (c==0.0){
        //        if(b<=0.0)
        //          theta=theta_max;
        //        else 
        //          theta=theta_min;
        //      }
        //      else if (c>0.0){
        //        if (b<-2*c*theta_max)
        //          theta = theta_max;
        //        else if (b>-2*c*theta_min)
        //          theta=theta_min;
        //        else
        //          theta=-b/(2*c);
        //      }
        //      else{
        //        if (b<-2*c*theta_min)
        //          theta=theta_max;
        //        else if (b>-2*c*theta_max)
        //          theta=theta_min;
        //        else
        //          {
        //            if ((aa+b*theta_min+c*theta_min*theta_min)>=(aa+b*theta_max+c*theta_max*theta_max))
        //              theta=theta_max;
        //            else theta=theta_min;
        //          }
        //      } // end else if c>0
        //      //      theta = 0.35;

        //      std::cout<<"\n Choice of theta = " << theta<<std::endl;
        //      *piezoLog<<"\t Choice of theta = " << theta<<std::endl;

        //      //      std::cout<<"\n"<<std::endl;

        //      //      for(UInt i=0;i<actNrParameter;i++)
        //      //        step[i]=theta*step[i]; 

        //      for (UInt i=0;i<actNrParameter;i++)
        //        stepR[i]=step[i].real();

        //      setNewParameterSet(parameter, parameter_new, scaling, theta, stepR, whichParameterToUpdate);

        //      std::cout<<"\n New parameter after backtracking ... " <<std::endl;      

        
        //      for(UInt i=0;i<nrParameter;i++)
        //        std::cout<<"paramter_new["<<i<<"]= " << parameter_new[i] << " ~ " << 
        //          100-(1-parameter[i]/parameterIncrement[i])*100 << " Prozent"<<std::endl;
        
        //      //              getchar();

        //      //      b = theta*b;
        //      eta = 1- theta*(1-eta);

        //      updateMaterialData(parameter_new,ptMaterial);
        //      createF(ptMaterial, F_hat,FALSE);

        //      //for(UInt i=0;i<F_hat.GetSize();i++)
        //      //y_hat_F_hat[i]=y_hat[i].real()-F_hat[i].real();
        //      //misfitnew=norm2Real(y_hat_F_hat);

        //      Res_outer=y_hat-F_hat;

        //      //      maxAndWeightedResNorm(Res_outer, maxres_outer, aval_new, y_hat);
        //      //              maxAndEuclNorm(Res_outer,maxres_outer,aval_new); 
        //      norm(Res_outer,aval_new,maxres_outer,y_hat);
        //      std::cout<<"\n || res || = " << aval_new;
        //      std::cout<<"\n || max_res || = " << maxres_outer;

        //      *piezoLog << "\t misfit after Linesearch " << aval_new << "(e_33, eps_33 = "<<parameter[7]<<", " <<parameter_new[9]<<")"<<std::endl;

        //      //      std::cout<<"\n misfitnew = " << misfitnew <<std::endl;


        //       }// end linesearch ...

     
        //choose eta
        //      eta_new=gamma*POW((misfitnew/misfit),al);
     
        eta_new = gamma*POW((aval/aval_old),al);
 
        //  std::cout<<"\n eta_new =  " << eta_new << std::endl;
        Double gammaEtaAl = gamma*POW(eta,al);
        if (gammaEtaAl > 0.1) // safeguard choice 2, see Pernice, Walker
          if (gammaEtaAl>eta)
            eta_new=gammaEtaAl;

        //eta_new=std::max(eta_new,gamma*POW(eta,al));
        //       std::cout<< "\n eta_new " << eta_new <<std::endl;
      

        if(eta_new<=2*(tau*delta)/sqrt(aval)); // final safeguard
        eta_new=0.8*(tau*delta)/sqrt(aval);
        //end choose eta

        if (eta_new>eta_max)
          eta=eta_max;
        else if (eta_new<eta_min)
          eta = eta_min;
        else
          eta=eta_new;


        std::cout<< "\n\n *** Choice of eta = " << eta<<std::endl; 
        *piezoLog<<"\t Choice of eta = " << eta<<std::endl;

        // end choose eta
        parameter=parameter_new;
        misfit=misfitnew;

        aval_old=aval;
        aval=aval_new;
        Res_inner = Res_outer;

        //      std::cout<<"\n before end newtoin ..."<<std::endl;

    }// end while Newton
    //   std::cout<<"\n leaving Newton CG 2 " <<std::endl;
    //std::cout<<"\n\n *** FINALLY CALCULATED PARAMETERS *** ... here they are:"<<std::endl;
    *piezoLog<<"\n\n *** FINALLY CALCULATED PARAMETERS *** ... here they are:"<<std::endl;

    for (UInt  i=0;i<parameter.GetSize();i++){
      //std::cout<<"par[" << i<<"]="<< parameter[i]<<";\n";
      *piezoLog<<"parameter["<<i<<"] = "<< parameter[i]<<std::endl;
   
    }

  }// end NewtonCG 3

  // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


  void piezoParamIdent::NewtonCG4(){   // here we scale in advance, just once before NewtonCG and afterwards
    ENTER_FCN("piezoParamIdent::NewtonCG4");
    //   std::cout<<"\n Entering piezoParamIdent::NewtonCG4 ... "<< std::endl;

    UInt nrNewtonIterations=0;
    UInt backtrackIterator=0;

    MaterialData * ptMaterial=ptMyPDE_->getPDEMaterialData();   // Pointer to MaterialData

    Double eta_max, eta_new, t, aa,b,c, theta_min, theta_max, gamma, al;
    Double alpha, beta, tau;
    tau = 1.1;
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
    Vector<Complex> step(actNrParameter+actNrParameterC);
    Vector<Double> stepR(actNrParameter);    
    Vector<Double> stepC(actNrParameterC);
    Vector<Complex> parstart (2*nrParameter);
    res_NE.Resize(actNrParameter+actNrParameterC);
    Vector<Complex> res_NE_rescaled(actNrParameter+actNrParameterC);
    bas.Resize(actNrParameter+actNrParameterC);
    basC.Resize(actNrParameter+actNrParameterC);
    Vector<Complex> basbar(nrMeasuredData);
    Vector<Complex> basbar_old(nrMeasuredData);
    Vector<Complex> res(nrMeasuredData);
    Vector<Complex> res0(nrMeasuredData);
    Vector<Complex> Res_outer(nrMeasuredData);
    Vector<Complex> Res_inner(nrMeasuredData);

    Complex new_res_res0;
    Double normres02, normres2, normresNEold2, normresNE2;

    Double aval, aval_new, alin_new, maxres_outer, maxres_inner, aval_old;

    // NO scaling in NewtonCG!! 
    for (UInt i=0;i<parameter.GetSize();i++){
      scaling[i]=1.0;
      scalingC[i]=1.0;
    }
    updateMaterialData(parameter, ptMaterial);
    updateComplexMaterialData(parameterC, ptMaterial);

    createF(ptMaterial,F_hat,TRUE);

    for (UInt i=0;i<actNrParameter+actNrParameterC;i++){
      bas[i]=1.0;
      basC[i]=Complex(1.0,1.0);
    }

    Res_outer=res=y_hat-F_hat;
    norm(Res_outer,aval,maxres_outer,y_hat);

    //    std::cout<<"\n Newton Step 0: Maxres = " << maxres_outer << ", aval = "<< aval<<std::endl;
    aval_new=aval;
    aval_old=aval;

    *piezoLog << "\n Initial misfit ||y-F|| = "<<aval<< " and maxres_outer = max |y-F| = " << maxres_outer << std::endl;

  
    // NEWTON ITERATION -- outer Loop!!
    while((aval>tau*delta||nrNewtonIterations<maxNumberNewtonLoops)&&nrNewtonIterations<maxNumberNewtonLoops){ // Newton
      //while(nrNewtonIterations<2){ // Newton
      *piezoLog << "\n Newton-Iteration: " << nrNewtonIterations <<std::endl;
      *piezoLog <<"------------------------"<<std::endl;
      //      std::cout<<"\n Newton-Iteration "<< nrNewtonIterations <<std::endl;

      nrNewtonIterations++;

      res0=res;
      normres02=norm2Real(res0);
      step.Resize(actNrParameter+actNrParameterC);
      res=y_hat-F_hat;

      updateMaterialData(parameter, ptMaterial);
      updateComplexMaterialData(parameterC, ptMaterial);

      // Create the matrices ...
      //      createJacobiMatrixC(JacobiMatrix);

      //      std::cout <<"\n NewtonCG4 ... 2 "<<std::endl;
      //       std::cout<<"Jacobi created ..." <<std::endl;

      createAdjointJacobiMatrix(JacobiMatrix,adjJacobiMatrix);
      //      std::cout<<"Adjoint Jacobi created ..." <<std::endl;
      //std::cout<<JacobiMatrix<<std::endl;
  
      for (UInt i=0;i<parameter.GetSize();i++){
        scaling[i]=1.0;
        scalingC[i]=1.0;
      }


      normres2 = a2norm(res);

      adjJacobiMatrix.Mult(res, res_NE);

      for (UInt i=0;i<actNrParameter;i++){
        res_NE[i]=(1.0/scaling[whichParToUpInd[i]]*scaling[whichParToUpInd[i]])*res_NE[i];
        res_NE_rescaled[i]=scaling[whichParToUpInd[i]]*res_NE[i];
      }


      for (UInt i=actNrParameter;i<actNrParameter+actNrParameterC;i++){
        res_NE[i]=(1.0/scalingC[whichParToUpIndC[i-actNrParameter]]*scalingC[whichParToUpIndC[i-actNrParameter]])*res_NE[i];
        res_NE_rescaled[i]=scaling[whichParToUpIndC[i-actNrParameter]]*res_NE[i];
      }


      basC = res_NE;

    
      alin_new = aval;
      *piezoLog << "\n \t alin_new = aval = " << alin_new<<std::endl;

      normresNEold2 = POW(a2norm(res_NE_rescaled),2);
      basbar_old.Resize(nrMeasuredData);

    
      Res_inner = Res_outer;

      // ----- INNER LOOP --- CG-ITERATION ---
      UInt nrCGIter=0;

      while((alin_new>eta*eta*aval)&&(nrCGIter<maxNumberInnerLoops)||nrCGIter<1){ // CG
        //      std::cout<<"\n CG Iteration " << nrCGIter << std::endl;
        nrCGIter++;

        JacobiMatrix.Mult(basC,basbar);

        //      for (UInt i=0;i<basC.GetSize();i++)
        //        std::cout<<"; basC = "<< basC[i]<<std::endl;

        //      for (UInt i=0;i<basbar.GetSize();i++)
        //        std::cout<<"; basbar = "<< basbar[i]<<std::endl;

        basbar_old=basbar;

        Double norm_bas_bar = POW(a2norm(basbar),2);

        alpha=normresNEold2 / norm_bas_bar;

        // std::cout<<"\n alpha = " << alpha << ",\t || normresNEold || = " << normresNEold2 << ",\t ||bas_bar|| = " << norm_bas_bar << std::endl;

        //      std::cout<<"\nstep = step + alpha *bas: "<<std::endl;
        for (UInt i=0;i<actNrParameter+actNrParameterC;i++){
          step[i]+=alpha*basC[i];
          //std::cout<<"step("<<i<<")= "<< step[i] << "; \t ";
        }


        //      std::cout<<"\n\nres = res - alpha*bas_bar: "<<std::endl;
        for (UInt i=0;i<res.GetSize();i++){
          Res_inner[i]=Res_inner[i]-alpha*basbar[i];
          res[i]=res[i]-alpha*basbar[i];
          //      std::cout<<"; res = "<< res[i]<<std::endl;
        }


        norm(Res_inner,alin_new,maxres_inner,y_hat);

        adjJacobiMatrix.Mult(res,res_NE);

        for(UInt i=0;i<actNrParameter;i++){
          res_NE[i]=(1.0/scaling[whichParToUpInd[i]]*scaling[whichParToUpInd[i]])*res_NE[i];
          res_NE_rescaled[i]=res_NE[i]*scaling[whichParToUpInd[i]];
          //       std::cout<<"; resNE_rescaled = "<< res_NE_rescaled[i]<<std::endl;

        }
        for(UInt i=actNrParameter;i<actNrParameter+actNrParameterC;i++){
          res_NE[i]=(1.0/scalingC[whichParToUpIndC[i-actNrParameter]]*scalingC[whichParToUpIndC[i-actNrParameter]])*res_NE[i];
          res_NE_rescaled[i]=scaling[whichParToUpIndC[i-actNrParameter]]*res_NE[i];
          //      res_NE[i]=(1.0/(scalingC[i-nrParameter]*scalingC[i-nrParameter]))*res_NE[i];
          //res_NE_rescaled[i]=res_NE[i]*scalingC[i-nrParameter];
          //       std::cout<<"; resNE_rescaledC = "<< res_NE_rescaled[i]<<std::endl;
        }


        //       std::cout<<"\n Successfully updated resNE and so on second time... "<<std::endl;

        //      std::cout<<"\n NewtonCG4 ... 9 "<<std::endl;
        normresNE2=POW(a2norm(res_NE_rescaled),2);
        beta = normresNE2/normresNEold2;
        normresNEold2 = normresNE2;

        //      std::cout<<"\nbas = resNE + beta * bas: "<<std::endl;   
        for(UInt i=0;i<actNrParameter;i++){
          basC[i]= res_NE[i]+beta*basC[i];
        }

        //      normres2=a2norm(res);
        norm(res,alin_new,maxres_inner,y_hat);
        *piezoLog << "\n \t \t alin_new  = " << alin_new;

      }// end while CG
      *piezoLog << "\n \t Number of CG - Iterations performed " << nrCGIter <<std::endl;

      parameter_new=parameter; 
      parameter_newC=parameterC;

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


      //      parameter_new[9]=parameter[9]+(1.0/scaling[9])*step[9].real();        // eps33      
      for (UInt i=0;i<actNrParameter;i++)
        stepR[i]=step[i].real();

      for (UInt i=actNrParameter;i<actNrParameter+actNrParameterC;i++)
        stepC[i-actNrParameter]=step[i].imag();
      //      std::cout<<" stepR = " <<stepR[i]<<"; stepC[i] = " << stepC[i]; 
  
      setNewParameterSet(parameter, parameter_new, scaling,theta,stepR,whichParameterToUpdate);
      setNewParameterSet(parameterC, parameter_newC, scalingC,theta,stepC,whichParameterToUpdateC);


      for (UInt i=0;i<nrParameter;i++){
        //      std::cout<<"Paramter_new["<<i<<"]= " << parameter_new[i]<< " + " << parameter_newC[i]<<"i"<<std::endl; 
      }
      //      std::cout<<"\n"<<std::endl;

      updateMaterialData(parameter_new, ptMaterial);
      updateComplexMaterialData(parameter_newC, ptMaterial);

      createF(ptMaterial, F_hat,FALSE);


      for(UInt i=0;i<F_hat.GetSize();i++)
        y_hat_F_hat[i]=y_hat[i]-F_hat[i];

      Res_outer = y_hat-F_hat;

      norm(Res_outer,aval_new,maxres_outer,y_hat);

      //        backtracking(eta, theta, s, a, a_lin_new);
      *piezoLog << "\t aval_new before Linesearch " << aval_new <<std::endl;

      backtrackIterator=0;

      while ((aval_new>POW((1.0-t*(1.0-eta)),2)*aval) && backtrackIterator<15 && eta<eta_max || backtrackIterator<1){ // Liniensuche
        //      std::cout<<"\n backtracking ... "<< backtrackIterator<< std::endl;
        backtrackIterator++;

        b=0.0;
        for(UInt i=0;i<nrMeasuredData;i++)
          b+= 2.0*(Res_outer[i]*(Res_inner[i]-Res_outer[i])).real();

        // choose theta:
        aa=aval;
        new_res_res0=Complex(0.0,0.0);

        c=aval_new-b-aval;

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

        for (UInt i=0;i<actNrParameter+actNrParameterC;i++)
          step[i]=theta*step[i];

        for (UInt i=0;i<actNrParameter;i++)
          stepR[i]=step[i].real();

        for (UInt i=actNrParameter;i<actNrParameter+actNrParameterC;i++)
          stepC[i-actNrParameter]=step[i].imag();
        

        setNewParameterSet(parameter, parameter_new, scaling,theta,stepR,whichParameterToUpdate);
        setNewParameterSet(parameterC, parameter_newC, scaling,theta,stepC,whichParameterToUpdateC);

        for(UInt i=0;i<actNrParameter+actNrParameterC;i++)
          step[i]=theta*step[i]; 

        for(UInt i=0;i<nrParameter;i++)
          std::cout<<"paramter_new["<<i<<"]= " << parameter_new[i]<<" + " << parameter_newC[i]<<std::endl;

        //      b = theta*b;
        eta = 1- theta*(1-eta);

        updateMaterialData(parameter_new,ptMaterial);
        updateComplexMaterialData(parameter_newC, ptMaterial);
        createF(ptMaterial, F_hat,FALSE);

        Res_outer=y_hat-F_hat;

        norm(Res_outer,aval_new,maxres_outer,y_hat);

        std::cout<<"\n New res = " <<  aval_new <<std::endl;
        //      getchar();

        //      *piezoLog << "\t aval_new after Linesearch " << aval_new <<std::endl;

      }// end linesearch ...

      //choose eta

      eta_new = gamma*POW((aval/aval_old),al);
      // std::cout<<"\n eta_new =  " << eta_new << std::endl;
      Double gammaEtaAl = gamma*POW(eta,al);
      if (gammaEtaAl > 0.1) // safeguard choice 2, see Pernice, Walker
        if (gammaEtaAl>eta)
          eta_new=gammaEtaAl;

      if(eta_new<=2*(tau*delta)/sqrt(aval)); // final safeguard
      eta_new=0.8*(tau*delta)/sqrt(aval);
      //end choose eta

      if (eta_new>eta_max)
        eta=eta_max;
      else
        eta=eta_new;

      //    std::cout<< "\n\n *** Choice of eta = " << eta<<std::endl; 
      //      *piezoLog<<"\t Choice of eta = " << eta<<std::endl;

      // end choose eta
      parameter=parameter_new;
      parameterC=parameter_newC;
   
      aval_old=aval;
      aval=aval_new;
      Res_inner = Res_outer;

    }// end while Newton
    //   std::cout<<"\n leaving Newton CG 4 " <<std::endl;
    //std::cout<<"\n\n *** FINALLY CALCULATED PARAMETERS *** ... here they are:"<<std::endl;
    //    *piezoLog<<"\n\n *** FINALLY CALCULATED PARAMETERS *** ... here they are:"<<std::endl;

    //        for (UInt  i=0;i<parameter.GetSize();i++){
    //std::cout<<"par[" << i<<"]="<< parameter[i]<<";\n";
    //    *piezoLog<<"parameter["<<i<<"] = "<< parameter[i]<<std::endl;
    //}

    //    *piezoLog <<"MaxResiduum: " << maxres_outer <<", weighted norm = " <<aval_new<< std::endl;
    finalnorm=maxres_outer;

  }// end NewtonCG 4


 //  void piezoParamIdent::NewtonCG2(){
//     ENTER_FCN("piezoParamIdent::NewtonCG2");
//     std::cout<<"\n Entering piezoParamIdent::NewtonCG2 ... "<< std::endl;

//     UInt nrNewtonIterations=0;
//     UInt backtrackIterator=0;

//     MaterialData * ptMaterial=ptMyPDE_->getPDEMaterialData();   // Pointer to MaterialData


//     Double eta_max, eta_new, t, aa,b,c, theta_min, theta_max, gamma, al;
//     Double a, a_lin_new, a_new, alpha, beta, tau;
//     tau = 1.5;
//     UInt i;
//     Double theta=0.5;

//     t=1.0e-4;
//     eta_max=0.5;
//     eta=0.5;
//     theta_min=0.1;
//     theta_max=0.5;
//     gamma=0.5;
//     al=1.5;
    
    
//     Vector<Complex> y_hat_F_hat (y_hat.GetSize());

//     //    real_misfit=norm_res;
//     //real_misfit_new=norm_res_new;

//     Vector<Double> step (nrParameter);
//     Vector<Double> parstart (nrParameter);
//     Vector<Complex> temp_res_NE(nrParameter);
//     Vector<Double> temp1(nrParameter);
//     Vector<Double> temp2(nrParameter);
//     Vector<Double> res_NE_old(nrParameter);
//     Vector<Double> res_NE_new(nrParameter);
//     Vector<Double> bas_old(nrParameter);
//     Vector<Double> bas_new(nrParameter);
//     // Vector<Complex> res_old(nrMeasuredData);
//     Vector<Complex> res(nrMeasuredData);
//     Vector<Complex> res_NE_new_compl(nrParameter);
//     Vector<Complex> new_res(nrMeasuredData);
//     Complex new_res_res0;
//     Vector<Complex> bas_old_compl(nrParameter);

//     parstart=parameter;

//     // Create the Matrices F, F', F*
//     updateMaterialData(parameter,ptMaterial);
//     createF(ptMaterial, F_hat,FALSE);

//     for (i=0; i<nrMeasuredData;i++)
//       y_hat_F_hat[i]=y_hat[i].real()-F_hat[i].real();

//     res=y_hat-F_hat;
      
//     a=sqrt(realA2norm(y_hat_F_hat)); 

//     std::cout << "\n a = |y-F| =  " <<a <<std::endl;

//     //    while((a>tau*delta*delta&&nrNewtonIterations<20)&&nrNewtonIterations<20){ // Newton
//     while((nrNewtonIterations<20)&&nrNewtonIterations<20){ // Newton

//       nrNewtonIterations++;
//       std::cout<<"\n Newton-Iteration "<< nrNewtonIterations <<std::endl;
//       step.Resize(nrParameter);

//       //            createJacobiMatrix(ptMaterial, F_hat, parameterIncrement,JacobiMatrix, solElecPot, solMechDispl);
//       createJacobiMatrix2(JacobiMatrix);

//       createAdjointJacobiMatrix(JacobiMatrix,adjJacobiMatrix);
//       //      testJacobiMatrix(F_hat, JacobiMatrix, parameter,  ptMaterial,parameterIncrement, solElecPot, solMechDispl);

//       adjJacobiMatrix.Mult(res, temp_res_NE);
      
//       std::cout<<"\n bas_old = resNE= scal^-2 * F'* res "<<std::endl;

//       for(i=0;i<nrParameter;i++){
//         //bas_old[i] = res_NE_old[i] = temp_res_NE[i].real();
//         bas_old[i] = res_NE_old[i] = 1.0/(scaling[i]*scaling[i])*temp_res_NE[i].real();
//         std::cout<<"bas_old = " << bas_old[i]<<std::endl;
//       }

//       a_lin_new=a;

//       std::cout<<"\n a_lin_new  = " << a <<std::endl;

//       UInt nrCGIter=0;

//       while((a_lin_new>eta*eta*a)&&(nrCGIter<10)){ // CG

//         std::cout<<"\n CG Iteration " << nrCGIter << std::endl;
//         nrCGIter++;

//         JacobiMatrix.MatVecMult_CD(bas_old,bas_bar);

//         for (i=0;i<bas_bar.GetSize();i++)
//           std::cout<<"bas_bar="<<bas_bar[i]<<"; ";
//         std::cout<<"\n temp_res_NE:"<<std::endl;
          
//         temp_res_NE.Resize(res_NE_old.GetSize());

//         for (i=0;i<res_NE_old.GetSize();i++){
//           temp_res_NE[i]=res_NE_old[i]*scaling[i];
//           std::cout<<"temp_res_NE="<<temp_res_NE[i]<<"; ";
//         }

//         std::cout<<"\n"<<std::endl;

//         Double norm_temp_res_NE = POW(a2norm(temp_res_NE),2);

//         Double norm_bas_bar = realA2norm(bas_bar);

//         alpha=norm_temp_res_NE / norm_bas_bar;

//         std::cout<<"\n alpha = " << alpha << ",   || temp_res_NE || = " << norm_temp_res_NE << ",    || bas_bar || = " << norm_bas_bar << std::endl;

//         std::cout<<"\n step = step + alpha *bas: "<<std::endl;

//         for (i=0;i<nrParameter;i++){
//           step[i]+=alpha*bas_old[i];
//           std::cout<<"step("<<i<<")= "<< step[i] << "; \t ";
//         }
//         std::cout<<"\n res = res - alpha*bas_bar "<<std::endl;

//         for (i=0;i<res.GetSize();i++){
//           res[i]=res[i]-alpha*bas_bar[i];
//           std::cout<<"res = "<<res[i]<<"; ";
//         }
//         std::cout<<"\n"<<std::endl;

//         //      std::cout<<"\n res.GetSize() " << res.GetSize() << "; res_NE_NEW_COMP.GetSize()= "<< res_NE_new_compl.GetSize()<<std::endl;

//         adjJacobiMatrix.Mult(res,res_NE_new_compl);

//         for(i=0;i<res_NE_new.GetSize();i++)
//           res_NE_new[i]=1.0/(scaling[i]*scaling[i])*res_NE_new_compl[i].real();

//         a_lin_new=POW(realA2norm(res),0.5);

//         std::cout<<"\n a_lin_new =  " << a_lin_new<<std::endl; 

//         temp1.Resize(nrParameter);
//         temp2.Resize(nrParameter);
//         for (i=0;i<nrParameter;i++){
//           temp1[i]=scaling[i]*res_NE_new[i];
//           temp2[i]=scaling[i]*res_NE_old[i];
//         }
//         beta = (POW(a2norm(temp1),2))/(POW(a2norm(temp2),2));
//         std::cout<<"\n beta = " << beta <<std::endl; 

//         std::cout<<"\n bas = res + beta * bas"<<std::endl;      

//         for(i=0;i<nrParameter;i++){
//           bas_new[i]= res_NE_new[i]+beta*bas_old[i];
//           std::cout<<"bas = "<<bas_new [i]<<"; ";
//         }


//         bas_old=bas_new;
//         res_NE_old=res_NE_new;

//       }// end while CG


//       std::cout<<"\n"<<std::endl;

//       for (i=0;i<nrParameter;i++){
//         std::cout<<"step("<<i<<")= " << step[i]<< "; \t";
//         parameter_new[i]=parameter[i]; //*scaling[i];
//         parameter_new[i]=parameter_new[i]+step[i];
//         //      parameter_new[i]=(1.0/scaling[i])*parameter_new[i];
//         std::cout<<"Paramter_new["<<i<<"]= " << parameter_new[i]<<std::endl;
//       }
//       std::cout<<"\n"<<std::endl;

//       // precautionary measure, if paramters tend to far away from initial values ...
//       for(i=0;i<nrParameter;i++)
//         if (std::abs(parameter_new[i]-parstart[i])>0.5*std::abs(parstart[i])){
//           parameter_new[i]=parameter[i];
//           std::cout<<"\n parameter("<<i<<") was reset to value " << parameter[i] <<std::endl;
//         }

//       updateMaterialData(parameter_new, ptMaterial);
//       createF(ptMaterial, F_hat,FALSE);

//       for(i=0;i<F_hat.GetSize();i++)
//         new_res[i]=y_hat[i]-F_hat[i];


//       //        backtracking(eta, theta, s, a, a_lin_new);
//       a_new=POW(realA2norm(new_res),0.5);

//       while ((a_new>=(1-t*(1-eta))*a) && backtrackIterator<20 && eta<eta_max){ // Liniensuche
//         std::cout<<"\n backtracking ... "<< backtrackIterator<< std::endl;
//         backtrackIterator++;

//         // choose theta:
//         aa=a*a;
//         new_res_res0=Complex(0.0,0.0);

//         for(i=0;i<nrMeasuredData;i++)
//           new_res_res0=new_res_res0+new_res[i]*(res[i].real()-new_res[i]);

//         b=2.0*new_res_res0.real();
//         std::cout<<"\n b = " << b << std::endl;
//         c=a_new*a_new-b-aa;
//         std::cout<<"\n c = "<<c<<std::endl;
//         if (c==0.0){
//           if(b<=0.0)
//             theta=theta_max;
//           else 
//             theta=theta_min;
//         }
//         else if (c>0.0){
//           if (b<-2*c*theta_max)
//             theta = theta_max;
//           else if (b>-2*c*theta_min)
//             theta=theta_min;
//           else
//             theta=-b/(2*c);
//         }
//         else{
//           if (b<-2*c*theta_min)
//             theta=theta_max;
//           else if (b>-2*c*theta_max)
//             theta=theta_min;
//           else
//             {
//               if ((aa+b*theta_min+c*theta_min*theta_min)>=(aa+b*theta_max+c*theta_max*theta_max))
//                 theta=theta_max;
//               else theta=theta_min;
//             }
//         } // end else if c>0

//         std::cout<<"\n Choice of theta = " << theta<<std::endl;
//         std::cout<<"\n"<<std::endl;

//         for(i=0;i<nrParameter;i++){
//           step[i]=theta*step[i];
//           parameter_new[i]=parameter[i]; //*scaling[i];
//           parameter_new[i]+=step[i];
//           //      parameter_new[i]=(1.0/scaling[i])*parameter_new[i];
//           std::cout<<"Paramter_new["<<i<<"]= " << parameter_new[i]<<std::endl;
//           eta = 1- theta*(1-eta);
//         }

//         updateMaterialData(parameter_new,ptMaterial);
//         createF(ptMaterial, F_hat,FALSE);

//         for(i=0;i<F_hat.GetSize();i++)
//           y_hat_F_hat[i]=y_hat[i]-F_hat[i];
//         a_new=POW(realA2norm(res),0.5);

//         std::cout<<"\n a_new = " << a_new <<std::endl;
        

//       }// end linesearch ...

     
//       //choose eta
//       eta_new=gamma*POW((a_new/a),al);
//       std::cout<<"\n eta_new =  " << eta_new << std::endl;

//       if (gamma*POW(eta,al) > 0.1) // safeguard choice 2, see Pernice, Walker
//         eta_new=std::max(eta_new,gamma*POW(eta,al));
//       //       std::cout<< "\n eta_new " << eta_new <<std::endl;

//       if (eta_new>eta_max)
//         eta_new=eta_max;

//       if(eta_new<=2*(tau*delta)/a_new) // final safeguard
//         eta_new=0.8*(tau*delta)/a_new;
//       //end choose eta
//       eta=eta_new;
//       if (eta>=eta_max)
//         eta=eta_max;

//       std::cout<< "\n\n *** Choice of eta = " << eta<<std::endl; 

//       // end choose eta
//       parameter=parameter_new;

//     }// end while Newton
//     std::cout<<"\n leaving Newton CG 2 " <<std::endl;
//     std::cout<<"\n\n *** FINALLY CALCULATED PARAMETERS *** ... here they are:"<<std::endl;

//     for (UInt i=0;i<parameter.GetSize();i++)
//       std::cout<<"par[" << i<<"]="<< parameter[i]<<";\n";

//   }// end NewtonCG 2




//   void piezoParamIdent::NewtonCG(){

//     ENTER_FCN("piezoParamIdent::Newton()");

//     // Settings for Newton-CG - routine
//     UInt nrIterations=0;


//     MaterialData * ptMaterial=ptMyPDE_->getPDEMaterialData();   // Pointer to MaterialData

//     updateMaterialData(parameter, ptMaterial);         //Writes initial guesses of parameters (read from MeasuredData.dat) to system

//     ptAlgsys = ptMyPDE_->getPDE_algsys();                       //Pointer to AlgebraicSystem
//     dofs=ptMyPDE_->getPDE_dofspernode();
//     numNodes= ptMyPDE_->getPDE_numPDENodes();

//     while (nrIterations<5) {
//       nrIterations++;
//       std::cout<<"\n ******************** \n NewtonCG ... Newton-Iteration-Nr = "<<nrIterations<<std::endl;

//       // Create the Matrices F, F', F*
//       createF(ptMaterial,  F_hat,FALSE);
//       //std::cout<<"Parameter-to-solution-map Matrix is built up!"<<std::endl;
//       createJacobiMatrix2(JacobiMatrix);
//       // std::cout<<"JacobiMatrix was created ..."<<std::endl;
//       createAdjointJacobiMatrix(JacobiMatrix,adjJacobiMatrix);
//       // std::cout<<"Adjoint JacMatrix was created ..."<<std::endl;
       
//       //std::cout<<"\n\ntest Jacobi Matrix ..."<<std::endl;
//       //  testJacobiMatrix(F_hat, JacobiMatrix, parameter, ptMaterial,parameterIncrement, solElecPot, solMechDispl);

//       Double real_misfit, real_misfit_new;
//       Complex misfit;
//       Vector<Double> y_hat_F_hat;
//       for (UInt i=0;i<y_hat.GetSize();i++)
//         y_hat[i].real()-F_hat[i].real();

//       real_misfit=POW(a2norm(y_hat_F_hat),2);
//       std::cout<<"\n Vor CG real(misfit) = real(y_hat-F_hat) = "<< real_misfit<<std::endl;

//       // we do a bit of scaling here ...
//       //      for(int i=0;i<scaling.GetSize();i++)
//       //        parameter[i]=parameter[i]*scaling[i];

//       // xxxxxxxxxxxxxxxx how should we solve the inner problem? xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

//       //      CG(); // calls CG Routine ... uses F, F', F'*    

//       std::cout<< "\n s determined with CG ... " <<std::endl;
//       for (UInt i=0;i<parameter.GetSize();i++){ 
//         std::cout<<"s("<<i<<")="<<s[i]<<"; \t";
//         parameter_new[i]=parameter[i]+s[i].real();        
//       }
//       std::cout<<"\n"<<std::endl;

//       // Matrix<Complex> JacobiMatrixNE(JacobiMatrix.GetSizeCol(), JacobiMatrix.GetSizeCol());
//       //       Vector<Complex> y_hat_F_hatNE(JacobiMatrix.GetSizeCol());

//       //       adjJacobiMatrix.Mult(JacobiMatrix,JacobiMatrixNE);
//       //       adjJacobiMatrix.Mult(y_hat_F_hat, y_hat_F_hatNE);


//       //         directSolve(JacobiMatrixNE, s , y_hat_F_hatNE);
      
//       //        std::cout<< "\n s determined with direct solver ... " <<std::endl;
//       //              for (int i=0;i<parameter.GetSize();i++){ 
//       //                std::cout<<"s("<<i<<")="<<s[i]<<"; \t";
//       //                parameter[i]=parameter[i]+s[i].real();        
//       //                std::cout<<"par_new("<<i<<")="<<parameter_new[i]<<"; \t";
//       //              }
//       //              std::cout<<"\n"<<std::endl;

//       // xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


//       // we do a bit of rescaling here ...
//       //      for(UInt i=0;i<scaling.GetSize();i++)
//       //        parameter[i]=parameter[i]*1.0/scaling[i];

//       updateMaterialData(parameter_new, ptMaterial); 

//       createF(ptMaterial, F_hat,FALSE);

//       // we do a bit of scaling here ...
//       //        for(UInt i=0;i<scaling.GetSize();i++)
//       //parameter_new[i]=parameter_new[i]*scaling[i];

//       for(UInt i=0;i<y_hat.GetSize();i++)
//         y_hat_F_hat=y_hat[i].real()-F_hat[i].real();
//       real_misfit_new=a2norm(y_hat_F_hat);
  
//       //here begins attempt to implement the backtracking Algo .... ;-)

//       //    backtracking(eta, theta, s, real_misfit, real_misfit_new);
     
//       std::cout<<"\n\n currently calculated parameters: " << std::endl;
      
//       // we do a bit of rescaling here ...
//       //     for(UInt i=0;i<scaling.GetSize();i++)
//       //        parameter[i]=parameter[i]/scaling[i];
//       //parameter=parameter_new;

//       for (UInt i=0;i<parameter.GetSize();i++){
//         std::cout<<"par[" << i<<"]="<< parameter[i]<<"; \t";
//         std::cout<<"\n"<<std::endl;
//       }
//       ///       parameter[i]=parameter[i]+s[i].real();
        
//                                                updateMaterialData(parameter,ptMaterial);

//     } // end while-Newton


//   } // end Newton


}// end namespace
