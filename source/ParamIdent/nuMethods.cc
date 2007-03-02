#include "PDE/SinglePDE.hh" 
#include "piezoParamIdent.hh"
#include "Utils/mathfunctions.hh"


namespace CoupledField
{

  // xxxxxxxxxxxxxxxxxxxxxxxxxx nuMethods  - ITERATION xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  void piezoParamIdent::nuMethods(){

    ENTER_FCN("piezoParamIdent::nuMethod()");

    UInt nrIterations=0;
    UInt nNuMethods=0;
    Double theta, eta_acc, nu, omega;

    updateMaterialData(parameter_);         //Writes initial guesses of parameters (read from MeasuredData.dat) to system

    Double normJacMat, old_res_outer, res_outer, new_res_inner, old_res_inner, new_res_outer, maxres_inner;
    Double relax;
    relax=10.0;

    Vector<Complex> act_res(nrMeasuredData);
    Vector<Complex> JacFs(nrMeasuredData);
    Vector<Complex> z(nrMeasuredData);
    Matrix<Complex> z_old(maxNumberInnerLoops_+2,nrMeasuredData);
    Vector<Double> stepR(actNrParameter);
    Vector<Complex> s_old(actNrParameter);
    Vector<Complex> JacFs_res(nrMeasuredData);
    Matrix<Complex> ImgSpaceScaling_Mat(nrMeasuredData,nrMeasuredData);
    Vector<Double> parameter_old(nrParameter_);

    updateMaterialData(parameter_);        
    createF(F_hat_,false);

    act_res = y_hat_-F_hat_;
    norm(act_res,new_res_outer,maxres_inner,y_hat_);

    std::cout<<"Norm of residual " << new_res_outer <<std::endl;
    res_outer=new_res_outer;

    *parLog<< new_res_outer; 

    for (UInt i=0;i<nrParameter_;i++)
      if (whichParameterToUpdate_[i]==1)
        *parLog<<"  " << parameter_[i];
    *parLog<<"  " << newtonCounter_;
    *parLog<<std::endl;
          
          
    nrIterations++;
    s.Resize(actNrParameter);
    s.Init();
    z.Resize(nrMeasuredData);
    z.Init();
    z_old.Resize(maxNumberInnerLoops_+2,nrMeasuredData);
    z_old.Init();
    parameter_old=parameter_;

      
    // Create the Matrices F, F', F*
    testJacobiMatrix2(F_hat_, JacobiMatrix_, parameter_, parameterIncrement_, solElecPot_, solMechDispl_);
   
    JacobiMatrix_=approxJacobiMatrix_;

    for (UInt i=0;i<nrMeasuredData;i++)
      for (UInt j=0;j<nrMeasuredData;j++)
        if (i==j)
          ImgSpaceScaling_Mat[i][j] = 1.0/std::log(Complex(real_[i],imag_[i]));


    createAdjointJacobiMatrix(JacobiMatrix_,adjJacobiMatrix_);
    Matrix<Complex> adjJacobiTemp;
    adjJacobiTemp = adjJacobiMatrix_*ImgSpaceScaling_Mat;
    adjJacobiMatrix_ = adjJacobiTemp;


    // XXXXXXXXXXXXXXX SPECTRUM OF F'*F XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX    
    if (false){

      Matrix<Complex> JacobiMatrixNE(JacobiMatrix_.GetSizeCol(), JacobiMatrix_.GetSizeCol());
      Matrix<Double> JacobiMatrixR(JacobiMatrix_.GetSizeCol(), JacobiMatrix_.GetSizeCol());
      Matrix<Double> JacobiMatrixNE_R(JacobiMatrix_.GetSizeCol(), JacobiMatrix_.GetSizeCol());
      Vector<Complex> y_hat_F_hat_NE(JacobiMatrix_.GetSizeCol());

      adjJacobiMatrix_.Mult(JacobiMatrix_,JacobiMatrixNE);
      adjJacobiMatrix_.Mult(act_res, y_hat_F_hat_NE);
      //std::cout<<JacobiMatrixNE<<std::endl;
      
      for (UInt i=0;i<JacobiMatrixNE.GetSizeRow();i++)
        for (UInt j=0;j<JacobiMatrixNE.GetSizeCol();j++){
          JacobiMatrixNE[i][j]=JacobiMatrixNE[i][j].real();
          JacobiMatrixNE_R[i][j]=JacobiMatrixNE[i][j].real();
        }
      std::cout<<"JacobiMatrixNE_R"<<std::endl;
      std::cout<<JacobiMatrixNE_R<<std::endl;

      for (UInt i=0;i<JacobiMatrix_.GetSizeRow();i++)
        for (UInt j=0;j<JacobiMatrix_.GetSizeCol();j++)
          JacobiMatrixR[i][j]=JacobiMatrix_[i][j].real();
      std::cout<<"JacobiMatrixR"<<std::endl;
      std::cout<<JacobiMatrixR<<std::endl;

      Vector<Double> eigenvalues(JacobiMatrixNE_R.GetSizeRow());
      eigenValues(JacobiMatrixNE_R,0.000001,eigenvalues);
      std::cout<<"\n Eigenvalues of F'*F:"<<std::endl;
      std::cout<<eigenvalues<<std::endl;

      for(UInt i=0;i<eigenvalues.GetSize();i++)
        std::cout <<" eig ("<<i<<") = "<< eigenvalues[i] << ", "<<std::endl; 

#ifdef USE_LAPACK
      Vector<Double> eigenvaluesWL;
      JacobiMatrixNE.eigenvaluesWithLapack(eigenvaluesWL);
      std::cout<<"eigenvaluesWL:"<<std::endl;
      std::cout<<eigenvaluesWL<<std::endl;
#endif
      getchar();

    } // end if true/false

    // XXXXXXXXXXXXXXX END SPECTRUM OF F'*F XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX    

    
    // TEST MAT_MULT 
          
    Double w=0.9;
    normJacMat=calcEuclidianMatrixNorm(JacobiMatrix_);

    while (w>=1/(normJacMat*normJacMat))
      w=0.9*w;
 
    JacobiMatrix_ *= Complex(w,w);
     adjJacobiMatrix_ *= Complex(w,w);
      
    //Norm ersetzt:
    norm(act_res,old_res_outer,maxres_inner,y_hat_);
     new_res_outer = old_res_outer;
    
    new_res_inner=old_res_outer;
    old_res_inner=old_res_outer;

    while(nNuMethods<maxNumberInnerLoops_){
      s_old=s;
        
      nNuMethods++;
      s.Resize(actNrParameter);
      s.Init();
       old_res_inner=new_res_inner;

       nu=1.5;

      // fitting am dicken mode
      relax=10.0;
     
      if (newtonCounter_>=25)
        relax=5.0;
      if (newtonCounter_>=30)
        relax=5.0;

      // fitting am radial mode
      //      relax=100.0;

      eta_acc = ((nNuMethods-1)*(2*nNuMethods-3)*(2*nNuMethods+2*nu-1))/
        ((nNuMethods+2*nu-1)*(2*nNuMethods+4*nu-1)*(2*nNuMethods+2*nu-3));
      omega = 4*((2*nNuMethods + 2*nu -1)*(nNuMethods+nu-1))/((nNuMethods+2*nu-1)*(2*nNuMethods+4*nu-1));

      for(UInt i=0;i<nrMeasuredData;i++){
        if (nNuMethods>1){
          z[i]=z[i]+relax*(eta_acc*(z[i]-z_old[nNuMethods-2][i])+omega*(act_res[i]-JacFs[i]));
        }
        else{
          z[i]=relax*(4*nu+2)/(4*nu+1)*act_res[i];
          //     std::cout<<"s("<<i<<")= "<<s[i]<<"; "<<std::endl;
        }
      }
   
      adjJacobiMatrix_.Mult(z,s);
         
      for(UInt i=0;i<nrMeasuredData;i++)
        z_old[nNuMethods][i]=z[i];
     
      JacobiMatrix_.Mult(s,JacFs);

      //F'(p^k)(s^k)-(y-F(p^k))
      for (UInt i=0;i<nrMeasuredData;i++){
        act_res[i]=y_hat_[i]-F_hat_[i];
        JacFs_res[i]=JacFs[i]-act_res[i];
      }

      norm(JacFs_res,new_res_inner,maxres_inner,y_hat_);
           

      if (new_res_inner>1.15*old_res_inner){
              std::cout<<"\n " <<  nNuMethods <<" inner iterations during Newton-step  " 
                 << newtonCounter_ <<" with " << nrMeasuredData<< " nrMeasuredData"<<std::endl;
              
        // The following part prevents stagnation of algorithm in case that
        // no better update was found during the first iteration step.
        
        Complex sum=Complex(0.0,0.0);
         for (UInt i=0;i<s_old.GetSize();i++)
           sum+=s_old[i];
         if(std::abs(sum)!=0.0)
           s=s_old;
         else
           for (UInt i=0;i<s_old.GetSize();i++)
             s[i]=0.1*s[i];

        break;
      }
        
    } // end while nuMethod ...


    nNuMethods=0;
    computeScaling();

    for (UInt i=0;i<actNrParameter;i++)
      stepR[i]=s[i].real();

    theta=5.0;

    parameter_old=parameter_;
    std::cout<<"stepR"<<std::endl;
    std::cout<<stepR<<std::endl;
    setNewParameterSet(parameter_, parameter_, scaling_, theta, stepR, whichParameterToUpdate_);
    updateMaterialData(parameter_);
    createF(F_hat_,false);
    
    for (UInt i=0;i<nrMeasuredData;i++)
      act_res[i]=y_hat_[i]-F_hat_[i];

    norm(act_res,new_res_outer,maxres_inner,y_hat_);

    std::cout<<"new_res_outer:"<<std::endl;
    std::cout<<new_res_outer<<std::endl;

    Integer lineSearchCount=0;

    while (new_res_outer>old_res_outer){
      theta = 0.5*theta;
      std::cout<<"theta = "<<theta<<std::endl;
      parameter_=parameter_old;
      setNewParameterSet(parameter_, parameter_, scaling_, theta, stepR, whichParameterToUpdate_);
      updateMaterialData(parameter_);
      createF(F_hat_,false);
      
      for (UInt i=0;i<nrMeasuredData;i++)
        act_res[i]=y_hat_[i]-F_hat_[i];
 
      norm(act_res,new_res_outer,maxres_inner,y_hat_);
      std::cout<<"new_res_outer = " << new_res_outer <<std::endl;
      std::cout<<"res_outer = " << res_outer <<std::endl;
      
      lineSearchCount++;
      if (lineSearchCount>12)
        break;

    }


    // if no backtracking is specified, please include the following lines!
    for (UInt i=0;i<nrParameter_;i++){
      //      parameter_new_[i]=scaling_[i]*parameter_[i];
      //      parameter_new_[i]+=s[i].real();
      //      parameter_[i]=1/scaling_[i]*parameter_new_[i];
      if (whichParameterToUpdate_[i]==1)
        std::cout<<" paramter("<<i<<") = " << parameter_[i]<<
          " ( ~ "<< 100-(1-parameter_[i]/parameterIncrement_[i])*100<<" Prozent) "<< std::endl;
    }
    // parameter_=parameter_new_;
      
    updateMaterialData(parameter_);
    createF(F_hat_,false);

    for (UInt i=0;i<nrMeasuredData;i++)
      act_res[i]=y_hat_[i]-F_hat_[i];
    //Norm ersetzt:
    //      std::cout<<act_res<<std::endl;
    norm(act_res,new_res_outer,maxres_inner,y_hat_);

    if(new_res_outer<=1.0e-12)
      getchar();

    //       while (new_res_outer>old_res_outer){
    //      std::cout<<"\n Warning: residual norm gets worse!" <<std::endl;

    //      parameter_=parameter_old;
    //      theta = 0.1*theta;
    //      inner_eta_=0.99*inner_eta_;
    //      setNewParameterSet(parameter_, parameter_, scaling_, theta, stepR, whichParameterToUpdate_);
    //      updateMaterialData(parameter_, ptMaterial_);
    //      createF(ptMaterial_, F_hat_,false);
    //      for (UInt i=0;i<y_hat_.GetSize();i++)
    //        act_res[i]=y_hat_[i]-F_hat_[i];
    //      //Norm ersetzt:
    //      norm(act_res,new_res_outer,maxres_inner,y_hat_);
    //      std::cout<<"\n Drifted away! So theta = " << theta << std::endl;
    //      std::cout<<"\n NewresOuter = " << new_res_outer << "old_res_outer = " << old_res_outer<< std::endl;
    //      std::cout<<"\n inner_eta_= " << inner_eta_ <<std::endl;
    //      //      getchar();
    //      //parameter_=parameter_new_;
    //      //      NewtonNuMethods();
    //       }

    theta = 1.0;
    //       Double t=1.0e-04;
    //       Double theta_min=0.1;
    //       Double theta_max=0.5;
    eta=1.15;
    //      UInt Backtrackiterator;

    //        while(new_res_outer >(1-t*(1-eta)*old_res_outer)&& Backtrackiterator<5){
    //       Backtrackiterator++;

    //       Double b,c;
    //       b=0.0; c=0.0;
    //       for(UInt i=0;i<nrMeasuredData;i++)
    //         b+= 2.0*(act_res[i]*(JacFs_res[i]-act_res[i])).real();


    //      // choose theta:
    //       Double aa=old_res_outer;
    //       Double aval=old_res_outer;
    //       Double aval_new = new_res_outer;

    //       //      new_res_res0=Complex(0.0,0.0);
         
    //       c=aval_new-b-aval;
         
    //       if (c==0.0){
    //         if(b<=0.0)
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
    //       std::cout<<"\n choice of theta = " << theta <<std::endl;

    //       for (UInt i=0;i<stepR.GetSize();i++)
    //         stepR[i]=theta*stepR[i];
    //       Double theta_temp=1.0;
    //       //      parameter_=parameter_old;
    //       //      setNewParameterSet(parameter_, parameter_, scaling_,theta_temp , stepR, whichParameterToUpdate_);
    //       //std::cout<<parameter_<<std::endl;
                 
    //       //updateMaterialData(parameter_, ptMaterial_);
    //       //createF(ptMaterial_, F_hat_,false);
         
    //       for (UInt i=0;i<y_hat_.GetSize();i++)
    //         act_res[i]=y_hat_[i]-F_hat_[i];
    //       //Norm ersetzt:
    //       norm(act_res,new_res_outer,maxres_inner,y_hat_);
    //        }

    parameter_old=parameter_;
    residuumParIdent_=new_res_outer;

  } // end NewtonNuMethods


 
  // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  //
  // ------------------------ nu METHODS-C II --------------------------------------------------------------------
  //
  // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  void piezoParamIdent::nuMethodsC2(){

    ENTER_FCN("piezoParamIdent::nuMethodsC2");

    UInt nrIterations=0;
    UInt nNuMethods=0;
    Double theta, eta_acc, nu, omega;

    updateMaterialData(parameter_);         //Writes initial guesses of parameters (read from MeasuredData.dat) to system
    updateComplexMaterialData(parameterC_);         //Writes initial guesses of parameterC_

    Double normJacMat, max_res_inner, maxres_inner, old_res_outer, new_res_inner, old_res_inner, new_res_outer;

    Vector<Complex> act_res(nrMeasuredData);
    Vector<Complex> JacFs(nrMeasuredData);
    Vector<Complex> z(nrMeasuredData);
    Matrix<Complex> z_old(maxNumberInnerLoops_+2,nrMeasuredData);
    Vector<Double> stepR(actNrParameter);
    Vector<Complex> s_old(actNrParameter+actNrParameterC);
    Vector<Double> stepC (actNrParameterC);
    Vector<Complex> JacFs_res(nrMeasuredData);
    Matrix<Complex> ImgSpaceScaling_Mat(nrMeasuredData, nrMeasuredData);
    Vector<Double> parameter_old(nrParameter_);
    Vector<Double> parameter_oldC(nrParameter_);

    bas.Resize(actNrParameter+actNrParameterC);
    bas.Init();
    basC.Resize(actNrParameter+actNrParameterC);
    basC.Init();

    for (UInt i=0;i<actNrParameter+actNrParameterC;i++){
      bas[i]=1.0;
      basC[i]=Complex(1.0,1.0);
    }

    createF(F_hat_,false);
    act_res = y_hat_-F_hat_;
    //    new_res_outer=old_res_outer=a2norm(act_res);
    norm(act_res, new_res_outer, maxres_inner,y_hat_);
    std::cout<<"\n Norm of residual = " << new_res_outer <<std::endl;
     
    nrIterations++;
    std::cout<<"\n Newton NuMethodsC ... Newton-Iteration-Nr = "<<nrIterations<<std::endl;
    s.Resize(actNrParameter+actNrParameterC);
    s.Init();
    z.Resize(nrMeasuredData);
    z.Init();
    z_old.Resize(maxNumberInnerLoops_+2,nrMeasuredData);
    z_old.Init();
    parameter_old = parameter_;
    parameter_oldC = parameterC_;
      
    // Create the Matrices F, F', F*
    createF(F_hat_,false);
    testJacobiMatrixC(F_hat_, JacobiMatrix_, parameter_);
    JacobiMatrix_=approxJacobiMatrix_;


    for (UInt i=0; i<nrMeasuredData;i++)
      for (UInt j=0; j<nrMeasuredData;j++)
        if(i==j)
          ImgSpaceScaling_Mat[i][j]=1.0/std::log(Complex(real_[i],imag_[i]));

    createAdjointJacobiMatrix(JacobiMatrix_,adjJacobiMatrix_);
    
    Matrix<Complex> adjJacobiTemp;
    adjJacobiTemp = adjJacobiMatrix_*ImgSpaceScaling_Mat;
    adjJacobiMatrix_ = adjJacobiTemp;
    
  
    // XXXXXXXXXXXXXXX SPECTRUM OF F'*F XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    if (false){

      Matrix<Complex> JacobiMatrixNE(JacobiMatrix_.GetSizeCol(), JacobiMatrix_.GetSizeCol());
      //      Matrix<Double> JacobiMatrixR(JacobiMatrix.GetSizeCol(), JacobiMatrix_.GetSizeCol());
      Matrix<Double> JacobiMatrixNE_R(JacobiMatrix_.GetSizeCol(), JacobiMatrix_.GetSizeCol());
      //       Vector<Complex> y_hat_F_hat_NE(JacobiMatrix.GetSizeCol());

      adjJacobiMatrix_.Mult(JacobiMatrix_,JacobiMatrixNE);
      //       adjJacobiMatrix_.Mult(act_res, y_hat_F_hat_NE);
      //       std::cout<<JacobiMatrixNE<<std::endl;
      
      for (UInt i=0;i<JacobiMatrixNE.GetSizeRow();i++)
        for (UInt j=0;j<JacobiMatrixNE.GetSizeCol();j++){
          JacobiMatrixNE[i][j]=JacobiMatrixNE[i][j].real();
          JacobiMatrixNE_R[i][j]=JacobiMatrixNE[i][j].real();
        }
 
      std::cout<<JacobiMatrixNE_R<<std::endl;
      Vector<Double> eigenvalues(JacobiMatrixNE_R.GetSizeRow());
      eigenValues(JacobiMatrixNE_R,0.000001,eigenvalues);
      std::cout<<"\n Eigenvalues of F'*F:"<<std::endl;
      std::cout<<eigenvalues<<std::endl;

#ifdef USE_LAPACK
      Vector<Double> eigenvaluesWL;
      JacobiMatrixNE.eigenvaluesWithLapack(eigenvaluesWL);
      std::cout<<"eigenvaluesWL:"<<std::endl;
      std::cout<<eigenvaluesWL<<std::endl;
#endif
      getchar();

    } // end if true/false

    // XXXXXXXXXXXXXXX SPECTRUM OF F'*F XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX




    // TEST MAT_MULT 
    
    Double w=0.9;

    normJacMat=calcEuclidianMatrixNorm(JacobiMatrix_);

    while (w>=1/(normJacMat*normJacMat))
      w=0.9*w;

    // if nu-method ..
    //JacobiMatrix = JacobiMatrix*Complex(w,w);
    JacobiMatrix_ *= Complex(w,w);
    //adjJacobiMatrix = adjJacobiMatrix*Complex(w,w);
    adjJacobiMatrix_ *= Complex(w,w);
      
    //    std::cout<<adjJacobiMatrix<<std::endl;
      
    // old_res_outer=a2norm(act_res);
    norm(act_res, old_res_outer, maxres_inner,y_hat_);
    new_res_outer = old_res_outer;

    new_res_inner=old_res_outer;
    old_res_inner=old_res_outer;

    while(nNuMethods<maxNumberInnerLoops_){
      //while(nnuMethods<maxNumberInnerLoops_){
      s_old=s;
      s.Resize(actNrParameter+actNrParameterC);
      s.Init();
      //  z.Resize(nrMeasuredData);

      nNuMethods++;
      //        std::cout <<"\n Here starts the nuMethods Iteration - Nr " << nnuMethods<< std::endl;
      old_res_inner=new_res_inner;

      //   F' * s^k
      JacobiMatrix_.Mult(s,JacFs);
        
      nu=1.5;
      Double relax=10.0;
      eta_acc = ((nNuMethods-1)*(2*nNuMethods-3)*(2*nNuMethods+2*nu-1))/
        ((nNuMethods+2*nu-1)*(2*nNuMethods+4*nu-1)*(2*nNuMethods+2*nu-3));
      omega = 4*((2*nNuMethods + 2*nu -1)*(nNuMethods+nu-1))/((nNuMethods+2*nu-1)*(2*nNuMethods+4*nu-1));
      
      for(UInt i=0;i<nrMeasuredData;i++){
        if (nNuMethods>1){
          z[i]=z[i]+relax*(eta_acc*(z[i]-z_old[nNuMethods-2][i])+omega*(act_res[i]-JacFs[i]));
        }
        else{
          z[i]=relax*(4*nu+2)/(4*nu+1)*act_res[i];
          //     std::cout<<"s("<<i<<")= "<<s[i]<<"; "<<std::endl;
        }
      }
     
      adjJacobiMatrix_.Mult(z,s);
      //  std::cout<<s<<std::endl;
              
      for(UInt i=0;i<nrMeasuredData;i++)
        z_old[nNuMethods][i]=z[i];     
         
      JacobiMatrix_.Mult(s,JacFs);
        
      //F'(p^k)(s^k)-(y-F(p^k))
      for (UInt i=0;i<nrMeasuredData;i++){
        act_res[i]=y_hat_[i]-F_hat_[i];
        JacFs_res[i]=JacFs[i]-act_res[i];
      }

      norm(JacFs_res, new_res_inner, max_res_inner,y_hat_);

      if (new_res_inner>1.15*old_res_inner){
        //      std::cout << " \n !! New_res_inner is worse than old_res_inner -> break of inner Loop! "<< std::endl;
        std::cout<<"\n " <<  nNuMethods <<" inner iterations during Newton-step  " 
                 << newtonCounter_ <<" with " << nrMeasuredData<< " nrMeasuredData"<<std::endl;
        //      getchar();
        
        // The following part prevents stagnation of algorithm in case that
        // no better update was found during the first iteration step.
        
        Complex sum=Complex(0.0,0.0);
         for (UInt i=0;i<s_old.GetSize();i++)
           sum+=s_old[i];
         if(std::abs(sum)!=0.0)
           s=s_old;
         else
           for (UInt i=0;i<s_old.GetSize();i++)
             s[i]=0.1*s[i];

        break;
      }
      //       }
        
    } // end while nuMethod ...


    nNuMethods=0;

    //     Double old_resid2=old_res_outer;
    //     Double new_resid2=new_res_outer;

    // backtracking(et , theta, s, old_resid2, new_resid2); 

    computeScaling();

    Double mult=1.0;
    
    for (UInt i=0;i<actNrParameter;i++)
      stepR[i]=s[i].real();

    std::cout<<actNrParameter<<std::endl;
    std::cout<<actNrParameterC<<std::endl;

    for (UInt i=actNrParameter;i<actNrParameter+actNrParameterC;i++)
      stepC[i-actNrParameter]=-1000.0*s[i].imag();

    std::cout<<"stepR:" <<std::endl;
    std::cout<<stepR<<std::endl;

    std::cout<<"stepC:" <<std::endl;
    std::cout<<stepC<<std::endl;

    
    //    parameter_new_=parameter;
    theta=1.0;
    Double thetaC=1.0;
    setNewParameterSet(parameter_, parameter_, scaling_, theta, stepR, whichParameterToUpdate_);
    setNewParameterSet(parameterC_, parameterC_, scalingC_, thetaC, stepC, whichParameterToUpdateC_);

    // if no backtracking is specified, please include the following lines!
    for (UInt i=0;i<nrParameter_;i++){
      std::cout<<" paramter("<<i<<") = " << parameter_[i]<<" + " << parameterC_[i] << " i " << std::endl;
    }
    // parameter=parameter_new_;
      
    updateMaterialData(parameter_);
    updateComplexMaterialData(parameterC_);

    createF(F_hat_,false);

    for (UInt i=0;i<nrMeasuredData;i++)
      act_res[i]=y_hat_[i]-F_hat_[i];
   
    norm(act_res, new_res_outer, max_res_inner,y_hat_);
    //      new_res_outer=(a2norm(act_res));

    std::cout<<"new_res_outer:"<<std::endl;
    std::cout<<new_res_outer<<std::endl;

    Integer lineSearchCount=0;

    while (new_res_outer>old_res_outer){
      theta = 0.5*theta;
      thetaC = 0.5*thetaC;
      std::cout<<"theta = "<<theta<<std::endl;
      parameter_=parameter_old;
      parameterC_=parameter_oldC;

      setNewParameterSet(parameter_, parameter_, scaling_, theta, stepR, whichParameterToUpdate_);
      setNewParameterSet(parameterC_, parameterC_, scalingC_, thetaC, stepC, whichParameterToUpdateC_);

      updateMaterialData(parameter_);
      updateComplexMaterialData(parameterC_);

      createF(F_hat_,false);
      
      for (UInt i=0;i<nrMeasuredData;i++)
        act_res[i]=y_hat_[i]-F_hat_[i];
 
      norm(act_res,new_res_outer,maxres_inner,y_hat_);
      std::cout<<"new_res_outer = " << new_res_outer <<std::endl;
      
      lineSearchCount++;
      if (lineSearchCount>10)
        break;

    }


    *parLog<<" "<< new_res_outer<<"  "; 

    for (UInt i=0; i<10;i++)
      if (whichParameterToUpdateC_[i]==1)
        *parLog<<"  " << parameterC_[i];

    *parLog<<std::endl;

    parameter_old=parameter_;
    parameter_oldC=parameterC_;
    residuumParIdent_=new_res_outer;
      
   
  } // end NewtonNuMethods

}






