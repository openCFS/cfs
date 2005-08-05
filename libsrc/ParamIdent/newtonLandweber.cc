
#include "PDE/SinglePDE.hh" 
#include "piezoParamIdent.hh"


namespace CoupledField
{

  // xxxxxxxxxxxxxxxxxxxxxxxxxx LANDWEBER  - ITERATION xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  void piezoParamIdent::NewtonLandweber(){

    ENTER_FCN("piezoParamIdent::NewtonLandweber()");
    std::cout<<"\n Entering piezoParamIdent::NewtonLandweber()"<<std::endl;
    //    std::cout<<"\n Landweber 1"<<std::endl;

    // Settings for Newton-CG - routine
    UInt nrIterations=0;
    UInt nLandweber=0;
    Double theta, eta_acc;

    MaterialData * ptMaterial=ptMyPDE_->getPDEMaterialData();   // Pointer to MaterialData
    updateMaterialData(parameter, ptMaterial);         //Writes initial guesses of parameters (read from MeasuredData.dat) to system

    ptAlgsys = ptMyPDE_->getPDE_algsys();                       //Pointer to AlgebraicSystem

    Double normJacMat, old_res_outer, new_res_inner, old_res_inner, new_res_outer;
    Matrix<Complex> Identity(actNrParameter, actNrParameter);
    //we need the Identity of size nrParameter x nrParameter and some other temporary matrices
    for (UInt i=0;i<actNrParameter;i++)
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
    createF(ptMaterial, F_hat,FALSE);
    act_res = y_hat-F_hat;
    new_res_outer=old_res_outer=a2norm(act_res);

    //    std::cout<<"maxNumberNewtonLoops = " << maxNumberNewtonLoops <<std::endl;
    //      getchar();
    //    while (new_res_outer<=old_res_outer && nrIterations<maxNumberNewtonLoops) {
    while (nrIterations<1) {

      for (UInt i=0;i<parameter.GetSize();i++)
        scaling[i]=1.0;
    
      nrIterations++;
      std::cout<<"\n NewtonLandweber ... Newton-Iteration-Nr = "<<nrIterations<<std::endl;
      s.Resize(actNrParameter);
      s_old.Resize(maxNumberInnerLoops,actNrParameter);
   

      // Create the Matrices F, F', F*
      createF(ptMaterial, F_hat,FALSE);
      createJacobiMatrix2(JacobiMatrix);
      createAdjointJacobiMatrix(JacobiMatrix,adjJacobiMatrix);

      // TEST MAT_MULT 

      // FF = F' F'*
      FF.Resize(adjJacobiMatrix.GetSizeRow(),JacobiMatrix.GetSizeCol());
      adjJacobiMatrix.Mult(JacobiMatrix,FF);

      //      Double symsysmat = calcEuclidianMatrixNorm(FF);
      //    std::cout<<"Norm of symmetric system matrix = " << symsysmat <<std::endl;


      //      Complex Jfs = Complex(0.0,0.0);

      Double w=0.9;

      normJacMat=calcEuclidianMatrixNorm(JacobiMatrix);

      while (w>=1/(normJacMat*normJacMat))
        w=0.9*w;

      // if nu-method ..
      JacobiMatrix = JacobiMatrix*Complex(w,w);
      adjJacobiMatrix = adjJacobiMatrix*Complex(w,w);
      
      // w=0.9;
      //    std::cout<<"Choice of w = "<< w << "; Norm JacMat*AdjJacMat = " << normJacMat << std::endl;
      //      getchar();

      // y_hat-F_hat
      //       for (UInt i=0;i<y_hat.GetSize();i++)
      //      act_res[i]=y_hat[i]-F_hat[i];
      //       std::cout<<"\n Landweber 4"<<std::endl;
          
      old_res_outer=a2norm(act_res);
      new_res_outer = old_res_outer;

      new_res_inner=old_res_outer;
      old_res_inner=old_res_outer;

      //    std::cout<<"Before Landweber, nLandweber = " << nLandweber << ", new_res = " << new_res_outer << ", eta * old_resid_outer = " << eta*old_res_outer <<std::endl;
    
      // LANDWEBER ITERATION   s^k+1 = s^k - w F'*(F's^k - (y-F))
      // while(new_res_inner<=old_res_inner&&nLandweber<maxNumberInnerLoops){
      while(nLandweber<maxNumberInnerLoops){
        nLandweber++;
        //      std::cout <<"\n Here starts the Landweber Iteration - Nr " << nLandweber<< std::endl;
        //      old_res_inner=new_res_inner;
        
        // F' * s^k
        JacobiMatrix.Mult(s,JacFs);
        //      //       std::cout<<"\n Landweber 5"<<std::endl;

        //      // F'sk - (y_hat-F_hat)
        //      for (UInt i=0;i<nrMeasuredData;i++){
        //        act_res[i]=y_hat[i]-F_hat[i];
        //        JacFs_res[i]=JacFs[i]-act_res[i];
        //      }
        
        // F'*(F's^k - (y-F))
        adjJacobiMatrix.Mult(JacFs_res,adjFF_res);
        
        //      std::cout<< "\n Parameter increments after Landweber Step " << nLandweber <<std::endl;
        //  s^k+1 = s^k - w F'*(F's^k - (y-F))
        //std::cout<<"adjFF_res.GetSize() = " << adjFF_res.GetSize()<<std::endl;
        //       std::cout<<"\n Landweber 6"<<std::endl;

        //       std::cout<<s_old<<std::endl;
        
        // if accelerated
        if (FALSE){
          for(UInt i=0;i<actNrParameter;i++){
            if (nLandweber>2)
              s[i]=s[i]+eta_acc*(s[i]-s_old[nLandweber-1][i])-100.0*w*adjFF_res[i];
            else
              s[i]=s[i]-100.0*w*adjFF_res[i];
            //         std::cout<<"s("<<i<<")= "<<s[i]<<"; "<<std::endl;
          }
          for(UInt i=0;i<actNrParameter;i++)
            s_old[nLandweber][i]=s[i];
        }
       
        // classical version of Landweber ....
        if (TRUE){
          for(UInt i=0;i<actNrParameter;i++){
            s[i]=s[i]-100.0*w*adjFF_res[i];
            //std::cout<<"s("<<i<<")= "<<s[i]<<"; "<<std::endl;
          }
        }
        //       std::cout<<"\n Landweber 6.4"<<std::endl;
        // JacFs = F'*s
        JacobiMatrix.Mult(s,JacFs);
        //       std::cout<<"\n Landweber 6.5"<<std::endl;

        //F'(p^k)(s^k)-(y-F(p^k))
        for (UInt i=0;i<nrMeasuredData;i++){
          act_res[i]=y_hat[i]-F_hat[i];
          JacFs_res[i]=JacFs[i]-act_res[i];
        }

        //      for (UInt i=0;i<nrMeasuredData;i++)
        //act_res[i]-=JacFs[i];
        
        new_res_inner=a2norm(JacFs_res);
        std::cout<<"\n new_res_inner = "<< new_res_inner <<", old_res_inner = " << old_res_inner<< std::endl;

        
        if (new_res_inner>4.5*old_res_inner){
          std::cout << " \n !! New_res_inner is worse than old_res_inner -> break of inner Loop! "<< std::endl;
          std::cout<<"\n Nr of Landweber = " << nLandweber <<std::endl;
          //      getchar();
          break;
        }



        //      //F'(p^k)(s^k)-(y-F(p^k))
        //      for (UInt i=0;i<nrMeasuredData;i++)
        //        act_res[i]-=JacFs[i];
        //      //       std::cout<<"\n Landweber 6.6"<<std::endl;
        
        //      //      for (UInt i=0;i<act_res.GetSize();i++)
        //      // std::cout<<"act_res= "<< act_res[i]<<std::endl;
        
        //      new_res_inner=a2norm(act_res);
        //      //             std::cout<<"\n Landweber 7"<<std::endl;


        //      if (new_res_inner>old_res_inner){
        //        std::cout << " \n !! New_res_inner is worse than old_res_inner!! "<< std::endl;
        //        //              break;
        //        //      getchar();
        //      }
        
        //      std::cout<<"\n new_res_inner = " << new_res_inner << ";  old_res_outer = " << old_res_outer << "; tau*delta= " << tau*delta << std::endl;
                        
        std::cout<<"\n end of inner Landweber Iter ..."<<std::endl;
      } // end while Landweber ...

      std::cout<<"\n Leaving Landwebers itaration ..."<<std::endl;

      std::cout<<"\n Number of inner iterations = " << nLandweber<<std::endl;

      nLandweber=0;


      // backtracking(et , theta, s, old_resid2, new_resid2); 

      theta = 1.0;
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

      for (UInt i=0;i<actNrParameter;i++)
        stepR[i]=s[i].real();
    
      //    parameter_new=parameter;
      theta=1.0;
      setNewParameterSet(parameter, parameter, scaling, theta, stepR, whichParameterToUpdate);

      // if no backtracking is specified, please include the following lines!
      for (UInt i=0;i<nrParameter;i++){
        //      parameter_new[i]=scaling[i]*parameter[i];
        //      parameter_new[i]+=s[i].real();
        //      parameter[i]=1/scaling[i]*parameter_new[i];
        std::cout<<" paramter("<<i<<") = " << parameter[i]<<std::endl;
      }
      // parameter=parameter_new;
      
      updateMaterialData(parameter, ptMaterial);
      createF(ptMaterial, F_hat,FALSE);


      for (UInt i=0;i<y_hat.GetSize();i++)
        act_res[i]=y_hat[i]-F_hat[i];
      new_res_outer=(a2norm(act_res));
      std::cout<<"\n new_res_outer = " << new_res_outer <<std::endl;

      if(new_res_outer<=1.0000e-08)
        getchar();

      if (new_res_outer>=old_res_outer){
        std::cout<<"\n Warning: residual norm gets worse!" <<std::endl;
        //      getchar();
        //parameter=parameter_new;
        //      NewtonLandweber();
      }
      

    } // end while nrIterations

    //    delete adjFJacF;
   
  } // end NewtonLandweber

  // --------------------------------------------------------------------------------------------------
  //
  // ---------  piezoParamIdent::NewtonLandweberC() ---------------------------------------------------
  //
  // --------------------------------------------------------------------------------------------------

  void piezoParamIdent::NewtonLandweberC(){

    ENTER_FCN("piezoParamIdent::NewtonLandweberC()");
    std::cout<<"\n Entering piezoParamIdent::NewtonLandweberC()"<<std::endl;

    UInt nrIterations=0;
    UInt nLandweber=0;
    Double theta, eta_acc;

    MaterialData * ptMaterial=ptMyPDE_->getPDEMaterialData();   // Pointer to MaterialData
    updateMaterialData(parameter, ptMaterial);         //Writes initial guesses of parameters (read from MeasuredData.dat) to system
    updateMaterialData(parameterC, ptMaterial);         //Writes initial guesses of parameters (read from MeasuredData.dat) to system

    ptAlgsys = ptMyPDE_->getPDE_algsys();                       //Pointer to AlgebraicSystem

    Double normJacMat, old_res_outer, new_res_inner, old_res_inner, new_res_outer;
    Matrix<Complex> Identity(actNrParameter+actNrParameterC, actNrParameter+actNrParameterC);

    //we need the Identity of size nrParameter x nrParameter and some other temporary matrices
    for (UInt i=0;i<actNrParameter+actNrParameterC;i++)
      Identity[i][i]=1;

    std::cout<<"\n actNrParameter " << actNrParameter << " actNrParameterC = " << actNrParameterC << std::endl;
    Matrix<Complex> adjFJacF(actNrParameter+actNrParameterC, 2*nrMeasuredData);
    Vector<Complex> adjFRes(actNrParameter+actNrParameterC);
    Vector<Complex> act_res(2*nrMeasuredData);
    Matrix<Complex> I_adjFF(actNrParameter,actNrParameter+actNrParameterC);
    Vector<Complex> adjFF_res(actNrParameter+actNrParameterC);
    Vector<Complex> JacFs(2*nrMeasuredData);
    Vector<Complex> JacFs_res(2*nrMeasuredData);
    Matrix<Complex> s_old(maxNumberInnerLoops,actNrParameter+actNrParameterC);
    Vector<Complex> s_temp(actNrParameter+actNrParameterC);
    Matrix<Complex> FF;
    Vector<Double> stepR(actNrParameter);
    Vector<Double> stepC(actNrParameterC);
    bas.Resize(actNrParameter+actNrParameterC);
    basC.Resize(actNrParameter+actNrParameterC);
    Matrix<Double> *matMat;
    Matrix<Double> *matMatC;

    if (considerMechDeformation==FALSE){
      adjFJacF.Resize(actNrParameter+actNrParameterC, actNrParameter+actNrParameterC);
      act_res.Resize(nrMeasuredData);
      JacFs.Resize(nrMeasuredData);
      JacFs_res.Resize(nrMeasuredData);
    }


    createF(ptMaterial, F_hat,TRUE);
    act_res = y_hat-F_hat;
    new_res_outer=old_res_outer=a2norm(act_res);

    for (UInt i=0;i<actNrParameter+actNrParameterC;i++){
      bas[i]=1.0;
      basC[i]=Complex(1.0,1.0);
    }

    std::cout<<"maxNumberNewtonLoops = " << maxNumberNewtonLoops <<std::endl;
    //      getchar();
    //    while (new_res_outer<=old_res_outer && nrIterations<maxNumberNewtonLoops) {
    while (nrIterations<1) {

      //  for (UInt i=0;i<parameter.GetSize();i++){
      //      scaling[i]=1.0;
      //      scalingC[i]=1.0;
      //       }
    
      nrIterations++;
      std::cout<<"\n NewtonLandweber ... Newton-Iteration-Nr = "<<nrIterations<<std::endl;
      s.Resize(actNrParameter+actNrParameterC);
      s_old.Resize(maxNumberInnerLoops,actNrParameter+actNrParameterC);
      // s_temp.Resize(actNrParameter);

      updateMaterialData(parameter, ptMaterial);
      updateComplexMaterialData(parameterC, ptMaterial);

      // Create the Matrices F, F', F*
      createF (ptMaterial, F_hat,FALSE);
      createJacobiMatrixC(JacobiMatrix);
      createAdjointJacobiMatrix(JacobiMatrix,adjJacobiMatrix);
     

      // TEST MAT_MULT 

      // FF = F' F'*
      FF.Resize(adjJacobiMatrix.GetSizeRow(),JacobiMatrix.GetSizeCol());
      adjJacobiMatrix.Mult(JacobiMatrix,FF);

      //      Double symsysmat = calcEuclidianMatrixNorm(FF);
      //          std::cout<<"Norm of symmetric system matrix = " << symsysmat <<std::endl;


      //      Complex Jfs = Complex(0.0,0.0);

      Double w=0.9;

      normJacMat=calcEuclidianMatrixNorm(JacobiMatrix);

      while (w>=1/(normJacMat*normJacMat))
        w=0.9*w;
      // w=0.9;
      //    std::cout<<"Choice of w = "<< w << "; Norm JacMat*AdjJacMat = " << normJacMat << std::endl;
      //      getchar();

      // y_hat-F_hat
      //       for (UInt i=0;i<y_hat.GetSize();i++)
      //      act_res[i]=y_hat[i]-F_hat[i];
      //       std::cout<<"\n Landweber 4"<<std::endl;
          
      old_res_outer=a2norm(act_res);
      new_res_outer = old_res_outer;

      new_res_inner=old_res_outer;
      old_res_inner=old_res_outer;

      //    std::cout<<"Before Landweber, nLandweber = " << nLandweber << ", new_res = " << new_res_outer << ", eta * old_resid_outer = " << eta*old_res_outer <<std::endl;
    
      // LANDWEBER ITERATION   s^k+1 = s^k - w F'*(F's^k - (y-F))
      //      while(new_res_inner<=old_res_inner&&nLandweber<maxNumberInnerLoops){
      while(nLandweber<maxNumberInnerLoops){

        nLandweber++;
        //      std::cout <<"\n Here starts the Landweber Iteration - Nr " << nLandweber<< std::endl;
        //      old_res_inner=new_res_inner;

        // F' * s^k
        JacobiMatrix.Mult(s,JacFs);
        //       std::cout<<"\n Landweber 5"<<std::endl;

        // F'sk - (y_hat-F_hat)
        for (UInt i=0;i<nrMeasuredData;i++){
          act_res[i]=y_hat[i]-F_hat[i];
          JacFs_res[i]=JacFs[i]-act_res[i];
        }
        
        // F'*(F's^k - (y-F))
        adjJacobiMatrix.Mult(JacFs_res,adjFF_res);
        
        //      std::cout<< "\n Parameter increments after Landweber Step " << nLandweber <<std::endl;
        //  s^k+1 = s^k - w F'*(F's^k - (y-F))
        //std::cout<<"adjFF_res.GetSize() = " << adjFF_res.GetSize()<<std::endl;
        //       std::cout<<"\n Landweber 6"<<std::endl;

        //       std::cout<<s_old<<std::endl;
        
        // if accelerated
        if (FALSE){
          for(UInt i=0;i<actNrParameter+actNrParameterC;i++){
            if (nLandweber>2)
              s[i]=s[i]+eta_acc*(s[i]-s_old[nLandweber-2][i])-100.0*w*adjFF_res[i];
            else
              s[i]=s[i]-100.0*w*adjFF_res[i];
            //         std::cout<<"s("<<i<<")= "<<s[i]<<"; "<<std::endl;
          }
          for(UInt i=0;i<actNrParameter;i++)
            s_old[nLandweber][i]=s[i];
        }
       
        //       std::cout<<"\n omega = " << w << std::endl;

        if (TRUE){
          for(UInt i=0;i<actNrParameter+actNrParameterC;i++){
            s[i]=s[i]-10.0*w*adjFF_res[i];
            //std::cout<<"s("<<i<<")= "<<s[i]<<"; "<<std::endl;
          }
        }
        //       std::cout<<"\n Landweber 6.4"<<std::endl;
        // JacFs = F'*s
        JacobiMatrix.Mult(s,JacFs);
        //       std::cout<<"\n Landweber 6.5"<<std::endl;
        //F'(p^k)(s^k)-(y-F(p^k))
        for (UInt i=0; i<nrMeasuredData; i++)
          act_res[i]-=JacFs[i];
        //       std::cout<<"\n Landweber 6.6"<<std::endl;
        
        //for (UInt i=0;i<act_res.GetSize();i++)
        //std::cout<<"act_res= "<< act_res[i]<<std::endl;
        
        //      new_res_inner=a2norm(JacFs_res);
        new_res_inner=a2norm(act_res);
        //       std::cout<<"\n Landweber 7"<<std::endl;

        if (new_res_inner>1.25*old_res_inner){
          std::cout << " \n !! New_res_inner is worse than old_res_inner!! "<< std::endl;
          break;
          //getchar();
        }
        
        //      std::cout<<"\n new_res_inner = " << new_res_inner << ";  old_res_outer = " << old_res_outer << "; tau*delta= " << tau*delta << std::endl;
                        
        std::cout<<"\n end of inner Landweber Iter ..."<<nLandweber << std::endl;
      } // end while Landweber ...


      nLandweber=0;



      // backtracking(et , theta, s, old_resid2, new_resid2); 

      theta = 1.0;
      matMat = ptMaterial->GetMatrix();
      matMatC = ptMaterial->GetMatrixC();
      
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
      scalingC[1]=0.1/((*matMatC)[2][2]);
      scalingC[2]=1.0/((*matMatC)[1][0]);
      scalingC[3]=1.0/((*matMatC)[0][2]);
      scalingC[4]=1.0/((*matMatC)[3][3]); 
      scalingC[5]=1.0/((*matMatC)[6][4]);
      scalingC[6]=std::abs(1.0/((*matMatC)[8][0]));
      scalingC[7]=0.01/((*matMatC)[8][2]);
      scalingC[8]=1.0/((*matMatC)[6][6]); 
      scalingC[9]=-0.001/((*matMatC)[8][8]);

      //       std::cout<<"\n s:" <<std::endl;
      //       std::cout<<s<<std::endl;

      for (UInt i=0;i<actNrParameter;i++)
        stepR[i]=s[i].real();

      //    std::cout<<stepR<<std::endl;

      for (UInt i=actNrParameter;i<actNrParameter+actNrParameterC;i++)
        stepC[i-actNrParameter]=s[i].imag();
      //       std::cout<<"\n stepC:"<<std::endl;
      //       std::cout<<stepC*scalingC[1]<<std::endl;
      // getchar();

    
      //    parameter_new=parameter;
      theta=1.0;
      setNewParameterSet(parameter, parameter, scaling, theta, stepR, whichParameterToUpdate);
      std::cout<<parameterC<<std::endl;
      setNewParameterSet(parameterC, parameterC, scalingC, theta, stepC, whichParameterToUpdateC);
  
      // if no backtracking is specified, please include the following lines!
      for (UInt i=0;i<nrParameter;i++){
        //      parameter_new[i]=scaling[i]*parameter[i];
        //      parameter_new[i]+=s[i].real();
        //      parameter[i]=1/scaling[i]*parameter_new[i];

        std::cout<<" parameter("<<i<<") = " << parameter[i]<< " + " << parameterC[i]<< " i " << std::endl;
      }
      // parameter=parameter_new;
      
      updateMaterialData(parameter, ptMaterial);
      updateComplexMaterialData(parameterC, ptMaterial);
      createF(ptMaterial, F_hat,FALSE);


      for (UInt i=0;i<y_hat.GetSize();i++)
        act_res[i]=y_hat[i]-F_hat[i];
      new_res_outer=(a2norm(act_res));
      std::cout<<"\n new_res_outer = " << new_res_outer <<std::endl;

      if(new_res_outer<=1.0e-07)
        getchar();

      if (new_res_outer>=old_res_outer){
        std::cout<<"\n Warning: residual norm gets worse!" <<std::endl;
        //      getchar();
        //parameter=parameter_new;
        //      NewtonLandweber();
      }
      
    } // end while nrIterations

    //    delete adjFJacF;
  

  } // end NewtonLandweberC

} // end namespace
