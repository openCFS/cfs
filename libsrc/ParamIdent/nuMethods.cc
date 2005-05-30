#include <iostream>
#include <fstream>
#include <string>
#include "DataInOut/GMV/outGMV.hh"
#include "General/environment.hh"
#include "PDE/SinglePDE.hh" 

#include "piezoParamIdent.hh"
#include "Forms/baseForm.hh"
#include "Utils/vector.hh"
#include "Utils/nodestoresol.hh"
#include "Utils/elemstoresol.hh"
#include "DataInOut/MaterialData.hh"
#include "PDE/timestepping.hh"
#include "Utils/baseelemstoresol.hh"
#include "Driver/singleDriver.hh"
#include "PDE/nodeEQN.hh"
#include <Domain/elem.hh>
#include "Forms/forms_header.hh"
#include "Utils/mathfunctions.hh"


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

namespace CoupledField
{

  // xxxxxxxxxxxxxxxxxxxxxxxxxx LANDWEBER  - ITERATION xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  void piezoParamIdent::nuMethods(){

    ENTER_FCN("piezoParamIdent::nuMethod()");
    //    std::cout<<"\n Entering piezoParamIdent::nuMethod()"<<std::endl;

    UInt nrIterations=0;
    UInt nNuMethods=0;
    Double theta, eta_acc, nu, omega;

    MaterialData * ptMaterial=ptMyPDE_->getPDEMaterialData();   // Pointer to MaterialData
    updateMaterialData(parameter, ptMaterial);         //Writes initial guesses of parameters (read from MeasuredData.dat) to system

    Double normJacMat, old_res_outer, new_res_inner, old_res_inner, new_res_outer, maxres_inner;
    Double relax;
    relax=1.0;

    Vector<Complex> act_res(nrMeasuredData);
    Vector<Complex> JacFs(nrMeasuredData);
    Vector<Complex> z(nrMeasuredData);
    Matrix<Complex> z_old(maxNumberInnerLoops+2,nrMeasuredData);
    Vector<Double> stepR(actNrParameter);
    Vector<Complex> s_old(actNrParameter);
    Vector<Complex> JacFs_res(nrMeasuredData);
    Matrix<Complex> ImgSpaceScalingMat(nrMeasuredData,nrMeasuredData);
    Vector<Double> parameter_old(nrParameter);



    //  while (nrIterations<maxNumberNewtonLoops){
    //      std::cout<<"\n Nr: "<< nrIterations << ", start next Newton NuMethod?"<<std::endl;
    

      updateMaterialData(parameter, ptMaterial);         //Writes initial guesses of parameters (read from MeasuredData.dat) to system
      createF(ptMaterial, F_hat,FALSE);
      act_res = y_hat-F_hat;
      std::cout<<"act_res = " <<std::endl;
//       std::cout<<act_res<<std::endl;
      std::cout<<y_hat<<std::endl;
      std::cout<<F_hat<<std::endl;
      //   getchar();

    //       std::cout<<"act_res = " <<std::endl;
    //       std::cout<<act_res<<std::endl;
    //       std::cout<<y_hat<<std::endl;
    //       std::cout<<F_hat<<std::endl;
    //   getchar();
      
    // Norm ersetzt:
    //      new_res_outer=old_res_outer=a2norm(act_res);      
    norm(act_res,new_res_outer,maxres_inner,y_hat);

    *parLog<< new_res_outer; 

    for (UInt i=0;i<nrParameter;i++)
      if (whichParameterToUpdate[i]==1)
        *parLog<<"  "<< parameter[i]/parameterIncrement[i];
    *parLog<<"  " << newtonCounter;
    *parLog<<std::endl;
          
          
    nrIterations++;
    //      std::cout<<"\n Newton NuMethods ... Newton-Iteration-Nr = "<<nrIterations<<std::endl;
    s.Resize(actNrParameter);
    z.Resize(nrMeasuredData);
    z_old.Resize(maxNumberInnerLoops+2,nrMeasuredData);
    parameter_old=parameter;

      
    // Create the Matrices F, F', F*
    //    createF(ptMaterial, F_hat,FALSE);
    //  createJacobiMatrix2(JacobiMatrix);
    //std::cout<<JacobiMatrix<<std::endl;
    testJacobiMatrix2(F_hat, JacobiMatrix, parameter, ptMaterial,parameterIncrement, solElecPot, solMechDispl);
    //std::cout<<approxJacobiMatrix<<std::endl;
    //std::cout<<JacobiMatrix<<std::endl;
    //      std::cout<<approxJacobiMatrix<<std::endl;
    JacobiMatrix=approxJacobiMatrix;
      

    for (UInt i=0;i<nrMeasuredData;i++)
      for (UInt j=0;j<nrMeasuredData;j++)
        if (i==j)
          ImgSpaceScalingMat[i][j] = 1.0/std::log(real[i]);
    //            ImgSpaceScalingMat[i][j] = Complex(1.0/real[i],1.0/imag[i]);

    //       std::cout<<"\n ImgspaceScalingMat"<<std::endl;

    //       std::cout<<ImgSpaceScalingMat<<std::endl;

    //       std::cout<<"\n Adjoint Matrix"<<std::endl;

    createAdjointJacobiMatrix(JacobiMatrix,adjJacobiMatrix);
    //std::cout<<adjJacobiMatrix<<std::endl;
    //      std::cout<<"\n Adjoint Matrix * ImgSpaceScalingMat"<<std::endl;
    adjJacobiMatrix = adjJacobiMatrix*ImgSpaceScalingMat;
      
    std::cout<<adjJacobiMatrix<<std::endl;
    
  
    // XXXXXXXXXXXXXXX SPECTRUM OF F'*F XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    //         Matrix<Complex> JacobiMatrixNE(JacobiMatrix.GetSizeCol(), JacobiMatrix.GetSizeCol());
    //      Matrix<Double> JacobiMatrixR(JacobiMatrix.GetSizeCol(), JacobiMatrix.GetSizeCol());
    //    Matrix<Double> JacobiMatrixNE_R(JacobiMatrix.GetSizeCol(), JacobiMatrix.GetSizeCol());
    //       Vector<Complex> y_hat_F_hatNE(JacobiMatrix.GetSizeCol());

    //    adjJacobiMatrix.Mult(JacobiMatrix,JacobiMatrixNE);
    //       adjJacobiMatrix.Mult(act_res, y_hat_F_hatNE);
    //       std::cout<<JacobiMatrixNE<<std::endl;
      
    //   for (UInt i=0;i<JacobiMatrixNE.GetSizeRow();i++)
    //      for (UInt j=0;j<JacobiMatrixNE.GetSizeCol();j++){
    //        JacobiMatrixNE[i][j]=JacobiMatrixNE[i][j].real();
    //        JacobiMatrixNE_R[i][j]=JacobiMatrixNE[i][j].real();
    //      }
    //       for (UInt i=0;i<JacobiMatrix.GetSizeRow();i++)
    //      for (UInt j=0;j<JacobiMatrix.GetSizeCol();j++)
    //        JacobiMatrixR[i][j]=JacobiMatrix[i][j].real();
    //       std::cout<<"JacobiMatrixR"<<std::endl;
    //       std::cout<<JacobiMatrixR<<std::endl;

    // //       std::cout<<JacobiMatrixNE_R<<std::endl;
    //    Vector<Double> eigenvalues(JacobiMatrixNE_R.GetSizeRow());
    //        eigenValues(JacobiMatrixNE_R,0.000001,eigenvalues);
    //        std::cout<<"\n Eigenvalues of F'*F:"<<std::endl;
    //        std::cout<<eigenvalues<<std::endl;
    //       getchar();

    //      for(UInt i=0;i<eigenvalues.GetSize();i++)
    // std::cout <<" eig ("<<i<<") = "<< eigenvalues[i] << ", "<<std::endl; 
    
    // TEST MAT_MULT 
          
    Double w=0.9;
    normJacMat=calcEuclidianMatrixNorm(JacobiMatrix);

    while (w>=1/(normJacMat*normJacMat))
      w=0.9*w;
    // if nu-method ..
    JacobiMatrix = JacobiMatrix*Complex(w,w);
    adjJacobiMatrix = adjJacobiMatrix*Complex(w,w);
      
    //Norm ersetzt:
    norm(act_res,old_res_outer,maxres_inner,y_hat);
    //      old_res_outer=a2norm(act_res);
    new_res_outer = old_res_outer;
    
    new_res_inner=old_res_outer;
    old_res_inner=old_res_outer;

    while(nNuMethods<maxNumberInnerLoops){
      s_old=s;
        
      nNuMethods++;
      s.Resize(actNrParameter);
      //      std::cout <<"\n Here starts the nuMethods Iteration - Nr " << nnuMethods<< std::endl;
      old_res_inner=new_res_inner;

      //      // F' * s^k
      //      JacobiMatrix.Mult(s,JacFs);
        
      //      // F'sk - (y_hat-F_hat)
      //      for (UInt i=0;i<nrMeasuredData;i++){
      //        act_res[i]=y_hat[i]-F_hat[i];
      //        JacFs_res[i]=JacFs[i]-act_res[i];
      //      }
        
      nu=1.5;

      relax=10.0;

      if (newtonCounter>=25)
        relax=5.0;
      //      if (newtonCounter>=30)
      //        relax=5.0;

      relax=10.0;
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
      //      std::cout<<"\n mu = " << eta_acc << ", omega = " << omega << ", nNuMehtods= " << nNuMethods << std::endl;
      //       std::cout<<z<<std::endl;
      adjJacobiMatrix.Mult(z,s);
      //      std::cout<<s<<std::endl;
         
         
      for(UInt i=0;i<nrMeasuredData;i++)
        z_old[nNuMethods][i]=z[i];

     
      JacobiMatrix.Mult(s,JacFs);
        
      //F'(p^k)(s^k)-(y-F(p^k))
      for (UInt i=0;i<nrMeasuredData;i++){
        act_res[i]=y_hat[i]-F_hat[i];
        JacFs_res[i]=JacFs[i]-act_res[i];
      }

      //      new_res_inner=a2norm(JacFs_res);        
      norm(JacFs_res,new_res_inner,maxres_inner,y_hat);

      //      std::cout<<"\n new_res_inner = "<< new_res_inner <<", old_res_inner = " << old_res_inner<< std::endl;

        
      //      if (new_res_inner>1.15*inner_eta*old_res_inner){
      if (new_res_inner>1.15*old_res_inner){
        //      std::cout << " \n !! New_res_inner is worse than old_res_inner -> break of inner Loop! "<< std::endl;
        std::cout<<"\n Nr of nuMethods during "<< newtonCounter << " Newton-step = " << nNuMethods <<std::endl;
        //      getchar();
        s=s_old;
        break;
      }
        
      //      std::cout<<"\n end of inner nuMethods Iter after "<<  nNuMethods << " iterations" <<std::endl;
    } // end while nuMethod ...


    nNuMethods=0;


    Matrix<Double> *matMat = ptMaterial->GetMatrix();
      
    scaling[0]=1.0/((*matMat)[0][0]); 
    scaling[1]=1.0/((*matMat)[2][2]);
    scaling[2]=1.0/((*matMat)[1][0]);
    scaling[3]=1.0/((*matMat)[0][2]);
    scaling[4]=1.0/((*matMat)[3][3]); 
    scaling[5]=1.0/((*matMat)[6][4]);
    scaling[6]=1.0/((*matMat)[8][0]);
    scaling[7]=1.0/((*matMat)[8][2]);
    scaling[8]=1.0/((*matMat)[6][6]); 
    scaling[9]=1.0/((*matMat)[8][8]);

    //     if(100-(1-parameter[6]/parameterIncrement[6])*100>=100)
    //      scaling[6]=scaling[6]*100;
    //       if(100-(1-parameter[4]/parameterIncrement[4])*100>=100)
    //      scaling[4]=scaling[4]*10;
    //       if(100-(1-parameter[8]/parameterIncrement[8])*100>=100)
    //      scaling[8]=scaling[8]*70;


    //   if (newtonCounter<=1){
    //       stepR[1]=s[1].real();
    //       stepR[6]=s[5].real();
    //       stepR[8]=s[7].real();
    //              theta=1.0;
    //        }
    //        else
    for (UInt i=0;i<actNrParameter;i++)
      stepR[i]=s[i].real();
    //      stepR[7]=s[0].real();
    //stepR[8]=s[1].real();       

    theta=1.0;

    ///std::cout<<stepR<<std::endl;
      //std::cout<<s<<std::endl;
      //    getchar();      

      setNewParameterSet(parameter, parameter, scaling, theta, stepR, whichParameterToUpdate);

      // if no backtracking is specified, please include the following lines!
      for (UInt i=0;i<nrParameter;i++){
        //      parameter_new[i]=scaling[i]*parameter[i];
        //      parameter_new[i]+=s[i].real();
        //      parameter[i]=1/scaling[i]*parameter_new[i];
        std::cout<<" paramter("<<i<<") = " << parameter[i]<< " ( ~ "<< 100-(1-parameter[i]/parameterIncrement[i])*100<<" Prozent) "<< std::endl;
      }
      // parameter=parameter_new;
      
      updateMaterialData(parameter, ptMaterial);
      createF(ptMaterial, F_hat,FALSE);

      for (UInt i=0;i<y_hat.GetSize();i++)
        act_res[i]=y_hat[i]-F_hat[i];
      //Norm ersetzt:
      //      std::cout<<act_res<<std::endl;
      norm(act_res,new_res_outer,maxres_inner,y_hat);
      //      new_res_outer=(a2norm(act_res));
      std::cout<<"\n Norm of residual = " << new_res_outer <<std::endl;
      //getchar();
      
      //  *parLog<< new_res_outer; 
      //       for (UInt i=0;i<nrParameter;i++)
      //      if (whichParameterToUpdate[i]==1)
      //        *parLog<<"  "<< parameter[i]/parameterIncrement[i];
      //       *parLog<<std::endl;

      if(new_res_outer<=1.0e-5)
        getchar();

      //       while (new_res_outer>old_res_outer){
      //      std::cout<<"\n Warning: residual norm gets worse!" <<std::endl;

      //      parameter=parameter_old;
      //      theta = 0.1*theta;
      //      inner_eta=0.99*inner_eta;
      //      setNewParameterSet(parameter, parameter, scaling, theta, stepR, whichParameterToUpdate);
      //      updateMaterialData(parameter, ptMaterial);
      //      createF(ptMaterial, F_hat,FALSE);
      //      for (UInt i=0;i<y_hat.GetSize();i++)
      //        act_res[i]=y_hat[i]-F_hat[i];
      //      //Norm ersetzt:
      //      norm(act_res,new_res_outer,maxres_inner,y_hat);
      //      std::cout<<"\n Drifted away! So theta = " << theta << std::endl;
      //      std::cout<<"\n NewresOuter = " << new_res_outer << "old_res_outer = " << old_res_outer<< std::endl;
      //      std::cout<<"\n inner_eta= " << inner_eta <<std::endl;
      //      //      getchar();
      //      //parameter=parameter_new;
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
      //       //      parameter=parameter_old;
      //       //      setNewParameterSet(parameter, parameter, scaling,theta_temp , stepR, whichParameterToUpdate);
      //       //std::cout<<parameter<<std::endl;
                 
      //       //updateMaterialData(parameter, ptMaterial);
      //       //createF(ptMaterial, F_hat,FALSE);
         
      //       for (UInt i=0;i<y_hat.GetSize();i++)
      //         act_res[i]=y_hat[i]-F_hat[i];
      //       //Norm ersetzt:
      //       norm(act_res,new_res_outer,maxres_inner,y_hat);
      //        }

      parameter_old=parameter;

  } // end NewtonNuMethodsy


  // xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  //
  //     xxxxx       nu MethodsC     xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  //
  // xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  void piezoParamIdent::nuMethodsC(){

    ENTER_FCN("piezoParamIdent::nuMethods(C)");
    std::cout<<"\n Entering piezoParamIdent::nuMethodsC()"<<std::endl;

    UInt nrIterations=0;
    UInt nNuMethods=0;
    Double theta, eta_acc, nu, omega;

    MaterialData * ptMaterial=ptMyPDE_->getPDEMaterialData();   // Pointer to MaterialData
    updateMaterialData(parameter, ptMaterial);         //Writes initial guesses of parameters (read from MeasuredData.dat) to system
    updateComplexMaterialData(parameterC, ptMaterial);         //Writes initial guesses of parameterC

    Double normJacMat, old_res_outer, new_res_inner, old_res_inner, new_res_outer, eta;
    eta =0.9;

    Vector<Complex> act_res(nrMeasuredData);
    Vector<Complex> JacFs(nrMeasuredData);
    Vector<Complex> z(nrMeasuredData);
    Matrix<Complex> z_old(maxNumberInnerLoops+2,nrMeasuredData);
    Vector<Double> stepR(actNrParameter);
    Vector<Complex> s_old(actNrParameter+actNrParameterC);
    Vector<Double> stepC (actNrParameterC);
    Vector<Complex> JacFs_res(nrMeasuredData);

    bas.Resize(actNrParameter+actNrParameterC);
    basC.Resize(actNrParameter+actNrParameterC);


    for (UInt i=0;i<actNrParameter+actNrParameterC;i++){
      bas[i]=1.0;
      basC[i]=Complex(1.0,1.0);
    }


    //  std::cout<<"\n Landweber 4"<<std::endl;
    createF(ptMaterial, F_hat,FALSE);
    act_res = y_hat-F_hat;
    new_res_outer=old_res_outer=a2norm(act_res);

    //    std::cout<<"maxNumberNewtonLoops = " << maxNumberNewtonLoops <<std::endl;
    //      getchar();
    //   while (new_res_outer<=old_res_outer && nrIterations<maxNumberNewtonLoops) {
    
    for (UInt i=0;i<parameter.GetSize();i++){
      scaling[i]=1.0;
      scalingC[i]=1.0;
    }
    
    nrIterations++;
    std::cout<<"\n Newton NuMethodsC ... Newton-Iteration-Nr = "<<nrIterations<<std::endl;
    s.Resize(actNrParameter+actNrParameterC);
    z.Resize(nrMeasuredData);
    z_old.Resize(maxNumberInnerLoops+2,nrMeasuredData);
      
    // Create the Matrices F, F', F*
    createF(ptMaterial, F_hat,FALSE);
    createJacobiMatrixC(JacobiMatrix);
    createAdjointJacobiMatrix(JacobiMatrix,adjJacobiMatrix);

    // XXXXXXXXXXXXXXX SPECTRUM OF F'*F XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    Matrix<Complex> JacobiMatrixNE(JacobiMatrix.GetSizeCol(), JacobiMatrix.GetSizeCol());
    //      Matrix<Double> JacobiMatrixR(JacobiMatrix.GetSizeCol(), JacobiMatrix.GetSizeCol());
    Matrix<Double> JacobiMatrixNE_R(JacobiMatrix.GetSizeCol(), JacobiMatrix.GetSizeCol());
    //       Vector<Complex> y_hat_F_hatNE(JacobiMatrix.GetSizeCol());

    adjJacobiMatrix.Mult(JacobiMatrix,JacobiMatrixNE);
    //       adjJacobiMatrix.Mult(act_res, y_hat_F_hatNE);
    //       std::cout<<JacobiMatrixNE<<std::endl;
      
    for (UInt i=0;i<JacobiMatrixNE.GetSizeRow();i++)
      for (UInt j=0;j<JacobiMatrixNE.GetSizeCol();j++){
        JacobiMatrixNE[i][j]=JacobiMatrixNE[i][j].real();
        JacobiMatrixNE_R[i][j]=JacobiMatrixNE[i][j].real();
      }
    //       for (UInt i=0;i<JacobiMatrix.GetSizeRow();i++)
    //      for (UInt j=0;j<JacobiMatrix.GetSizeCol();j++)
    //        JacobiMatrixR[i][j]=JacobiMatrix[i][j].real();
    //       std::cout<<"JacobiMatrixR"<<std::endl;
    //       std::cout<<JacobiMatrixR<<std::endl;

    // //       std::cout<<JacobiMatrixNE_R<<std::endl;
    //        Vector<Double> eigenvalues(JacobiMatrixNE_R.GetSizeRow());
    //        eigenValues(JacobiMatrixNE_R,0.000001,eigenvalues);
    //        std::cout<<"\n Eigenvalues of F'*F:"<<std::endl;
    //        std::cout<<eigenvalues<<std::endl;
    //       getchar();

    // TEST MAT_MULT 

      
    Double w=0.9;

    normJacMat=calcEuclidianMatrixNorm(JacobiMatrix);

    while (w>=1/(normJacMat*normJacMat))
      w=0.9*w;

    // if nu-method ..
    JacobiMatrix = JacobiMatrix*Complex(w,w);
    adjJacobiMatrix = adjJacobiMatrix*Complex(w,w);
      
      
    old_res_outer=a2norm(act_res);
    new_res_outer = old_res_outer;

    new_res_inner=old_res_outer;
    old_res_inner=old_res_outer;

    while(nNuMethods<maxNumberInnerLoops){
      //while(nnuMethods<maxNumberInnerLoops){
      s_old=s;
      s.Resize(actNrParameter+actNrParameterC);
      z.Resize(nrMeasuredData);

      nNuMethods++;
      //        std::cout <<"\n Here starts the nuMethods Iteration - Nr " << nnuMethods<< std::endl;
      old_res_inner=new_res_inner;

      // F' * s^k
      JacobiMatrix.Mult(s,JacFs);
        
      // F'sk - (y_hat-F_hat)
      for (UInt i=0;i<nrMeasuredData;i++){
        act_res[i]=y_hat[i]-F_hat[i];
        JacFs_res[i]=JacFs[i]-act_res[i];
      }
        
       
      nu=0.004515;

      eta_acc = ((nNuMethods-1)*(2*nNuMethods-3)*(2*nNuMethods+2*nu-1))/
        ((nNuMethods+2*nu-1)*(2*nNuMethods+4*nu-1)*(2*nNuMethods+2*nu-3));
      omega = 4*((2*nNuMethods + 2*nu -1)*(nNuMethods+nu-1))/((nNuMethods+2*nu-1)*(2*nNuMethods+4*nu-1));

      for(UInt i=0;i<nrMeasuredData;i++){
        if (nNuMethods>1){
          z[i]=z[i]+5.0*(eta_acc*(z[i]-z_old[nNuMethods-2][i])+omega*(act_res[i]-JacFs[i]));
          //         z[i]=z[i]+10.0*(Complex(eta_acc,eta_acc)*(z[i]-z_old[nNuMethods-2][i])+Complex(omega,omega)*(act_res[i]-JacFs[i]));
        }
        else{
          //         z[i]=10.0*Complex((4*nu+2)/(4*nu+1),(4*nu+2)/(4*nu+1))*act_res[i];
          z[i]=5.0*((4*nu+2)/(4*nu+1))*act_res[i];
          //       std::cout<<"s("<<i<<")= "<<s[i]<<"; "<<std::endl;
        }
      }
      std::cout<<"\n mu = " << eta_acc << ", omega = " << omega << ", nNuMehtods= " << nNuMethods << std::endl;
      std::cout<<z<<std::endl;
      adjJacobiMatrix.Mult(z,s);
      //  std::cout<<s<<std::endl;
         
         
      for(UInt i=0;i<nrMeasuredData;i++)
        z_old[nNuMethods][i]=z[i];
     
      JacobiMatrix.Mult(s,JacFs);
        
      //F'(p^k)(s^k)-(y-F(p^k))
      for (UInt i=0;i<nrMeasuredData;i++){
        act_res[i]=y_hat[i]-F_hat[i];
        JacFs_res[i]=JacFs[i]-act_res[i];
      }

      //        for (UInt i=0;i<nrMeasuredData;i++)
      //act_res[i]-=JacFs[i];
        
      new_res_inner=a2norm(JacFs_res);
      std::cout<<"\n new_res_inner = "<< new_res_inner <<", old_res_inner = " << old_res_inner<< std::endl;


        
      //       if (new_res_inner>1.1*old_res_inner){
      //      std::cout << " \n !! New_res_inner is worse than old_res_inner!! "<< std::endl;
      //      std::cout<<"\n Nr of nuMethodsC = " << nNuMethods <<std::endl;
      //      //        getchar();
      //      s=s_old;
      //      break;
      //       }
        
      std::cout<<"\n end of inner nuMethodsC Iter ..."<<std::endl;
    } // end while nuMethod ...


    nNuMethods=0;


    // backtracking(et , theta, s, old_resid2, new_resid2); 

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

    Double mult;
    //     if (nNuMethods<5)
    //       mult=1.5;
    //     else
    //       mult=1.0;

    scalingC[0]=1.0/((*matMatC)[0][0]); 
    scalingC[1]=1.0/((*matMatC)[2][2]);
    scalingC[2]=1.0/((*matMatC)[1][0]);
    scalingC[3]=1.0/((*matMatC)[0][2]);
    scalingC[4]=1.0/((*matMatC)[3][3]); 
    scalingC[5]=1.0/((*matMatC)[6][4]);
    scalingC[6]=std::abs(1.0/((*matMatC)[8][0]));
    scalingC[7]=mult/((*matMatC)[8][2]);
    scalingC[8]=1.0/((*matMatC)[6][6]); 
    scalingC[9]=mult/((*matMatC)[8][8]);

    for (UInt i=0;i<actNrParameter;i++)
      stepR[i]=s[i].real();

    for (UInt i=actNrParameter;i<actNrParameter+actNrParameterC;i++)
      stepC[i-actNrParameter]=-s[i].imag();
    
    //    parameter_new=parameter;
    theta=2.0;
    setNewParameterSet(parameter, parameter, scaling, theta, stepR, whichParameterToUpdate);
    setNewParameterSet(parameterC, parameterC, scalingC, theta, stepC, whichParameterToUpdateC);

    // if no backtracking is specified, please include the following lines!
    for (UInt i=0;i<nrParameter;i++){
      //        parameter_new[i]=scaling[i]*parameter[i];
      //        parameter_new[i]+=s[i].real();
      //        parameter[i]=1/scaling[i]*parameter_new[i];
      std::cout<<" paramter("<<i<<") = " << parameter[i]<<" + " << parameterC[i] << " i " << std::endl;
    }
    // parameter=parameter_new;
      
    updateMaterialData(parameter, ptMaterial);
    updateComplexMaterialData(parameterC, ptMaterial);
    createF(ptMaterial, F_hat,FALSE);


    for (UInt i=0;i<y_hat.GetSize();i++)
      act_res[i]=y_hat[i]-F_hat[i];
    new_res_outer=(a2norm(act_res));
    std::cout<<"\n new_res_outer = " << new_res_outer <<std::endl;

    *parLog<<" "<< new_res_outer<<"  "; 


    //     if(new_res_outer<=1.0e-08)
    //       getchar();


    if (new_res_outer>=old_res_outer){
      std::cout<<"\n Warning: residual norm gets worse!" <<std::endl;
      //        getchar();
      //parameter=parameter_new;
      //        NewtonNuMethods();
    }
      
  } // end while nrIterations

    //    delete adjFJacF;
   
  //} // end NewtonNuMethodsy

  // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  //
  // ------------------------ nu METHODS-C II --------------------------------------------------------------------
  //
  // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  void piezoParamIdent::nuMethodsC2(){

    ENTER_FCN("piezoParamIdent::nuMethodsC2");
    std::cout<<"\n MUMETHODS C 2 "<<std::endl;

    UInt nrIterations=0;
    UInt nNuMethods=0;
    Double theta, eta_acc, nu, omega;

    MaterialData * ptMaterial=ptMyPDE_->getPDEMaterialData();   // Pointer to MaterialData
    updateMaterialData(parameter, ptMaterial);         //Writes initial guesses of parameters (read from MeasuredData.dat) to system
    updateComplexMaterialData(parameterC, ptMaterial);         //Writes initial guesses of parameterC

    Double normJacMat, max_res_inner, maxres_inner, old_res_outer, new_res_inner, old_res_inner, new_res_outer;

    Vector<Complex> act_res(nrMeasuredData);
    Vector<Complex> JacFs(nrMeasuredData);
    Vector<Complex> z(nrMeasuredData);
    Matrix<Complex> z_old(maxNumberInnerLoops+2,nrMeasuredData);
    Vector<Double> stepR(actNrParameter);
    Vector<Complex> s_old(actNrParameter+actNrParameterC);
    Vector<Double> stepC (actNrParameterC);
    Vector<Complex> JacFs_res(nrMeasuredData);
    Matrix<Complex> ImgSpaceScalingMat(nrMeasuredData, nrMeasuredData);
    Vector<Double> parameter_old(nrParameter);
    Vector<Double> parameter_oldC(nrParameter);

    bas.Resize(actNrParameter+actNrParameterC);
    basC.Resize(actNrParameter+actNrParameterC);

    for (UInt i=0;i<actNrParameter+actNrParameterC;i++){
      bas[i]=1.0;
      basC[i]=Complex(1.0,1.0);
    }

    createF(ptMaterial, F_hat,FALSE);
    act_res = y_hat-F_hat;
    //    new_res_outer=old_res_outer=a2norm(act_res);
    norm(act_res, new_res_outer, maxres_inner,y_hat);
    *parLog<<new_res_outer;

    //    std::cout<<"maxNumberNewtonLoops = " << maxNumberNewtonLoops <<std::endl;
    //      getchar();
    //   while (new_res_outer<=old_res_outer && nrIterations<maxNumberNewtonLoops) {
    
    //  for (UInt i=0;i<parameter.GetSize();i++){
    //       scaling[i]=1.0;
    //       scalingC[i]=1.0;
    //     }
    
    nrIterations++;
    std::cout<<"\n Newton NuMethodsC ... Newton-Iteration-Nr = "<<nrIterations<<std::endl;
    s.Resize(actNrParameter+actNrParameterC);
    z.Resize(nrMeasuredData);
    z_old.Resize(maxNumberInnerLoops+2,nrMeasuredData);
    parameter_old = parameter;
    parameter_oldC = parameterC;
      
    // Create the Matrices F, F', F*
    createF(ptMaterial, F_hat,FALSE);
    //    createJacobiMatrixC(JacobiMatrix);
    testJacobiMatrixC(F_hat, JacobiMatrix, parameter, ptMaterial);
    JacobiMatrix=approxJacobiMatrix;
    

    for (UInt i=0; i<nrMeasuredData;i++)
      for (UInt j=0; j<nrMeasuredData;j++)
        if(i==j)
          ImgSpaceScalingMat[i][j]=1.0/std::log(real[i]);
    
    createAdjointJacobiMatrix(JacobiMatrix,adjJacobiMatrix);
    
    adjJacobiMatrix = adjJacobiMatrix*ImgSpaceScalingMat;

    //    std::cout<<approxJacobiMatrix<<std::endl;
    // std::cout<<adjJacobiMatrix<<std::endl;
    //getchar();



    // XXXXXXXXXXXXXXX SPECTRUM OF F'*F XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    //  Matrix<Complex> JacobiMatrixNE(JacobiMatrix.GetSizeCol(), JacobiMatrix.GetSizeCol());
    //      Matrix<Double> JacobiMatrixR(JacobiMatrix.GetSizeCol(), JacobiMatrix.GetSizeCol());
    // Matrix<Double> JacobiMatrixNE_R(JacobiMatrix.GetSizeCol(), JacobiMatrix.GetSizeCol());
    //       Vector<Complex> y_hat_F_hatNE(JacobiMatrix.GetSizeCol());
    
    //adjJacobiMatrix.Mult(JacobiMatrix,JacobiMatrixNE);
    //       adjJacobiMatrix.Mult(act_res, y_hat_F_hatNE);
    //       std::cout<<JacobiMatrixNE<<std::endl;
    
    //   for (UInt i=0;i<JacobiMatrixNE.GetSizeRow();i++)
    //       for (UInt j=0;j<JacobiMatrixNE.GetSizeCol();j++){
    //      JacobiMatrixNE[i][j]=JacobiMatrixNE[i][j].real();
    //      JacobiMatrixNE_R[i][j]=JacobiMatrixNE[i][j].real();
    //      }
    //       for (UInt i=0;i<JacobiMatrix.GetSizeRow();i++)
    //      for (UInt j=0;j<JacobiMatrix.GetSizeCol();j++)
    //        JacobiMatrixR[i][j]=JacobiMatrix[i][j].real();
    //       std::cout<<"JacobiMatrixR"<<std::endl;
    //       std::cout<<JacobiMatrixR<<std::endl;

    // //       std::cout<<JacobiMatrixNE_R<<std::endl;
    //        Vector<Double> eigenvalues(JacobiMatrixNE_R.GetSizeRow());
    //        eigenValues(JacobiMatrixNE_R,0.000001,eigenvalues);
    //        std::cout<<"\n Eigenvalues of F'*F:"<<std::endl;
    //        std::cout<<eigenvalues<<std::endl;
    //       getchar();

    // TEST MAT_MULT 
    
    Double w=0.9;

    normJacMat=calcEuclidianMatrixNorm(JacobiMatrix);

    while (w>=1/(normJacMat*normJacMat))
      w=0.9*w;

    // if nu-method ..
    JacobiMatrix = JacobiMatrix*Complex(w,w);
    adjJacobiMatrix = adjJacobiMatrix*Complex(w,w);
      
    std::cout<<adjJacobiMatrix<<std::endl;
      
    // old_res_outer=a2norm(act_res);
    norm(act_res, old_res_outer, maxres_inner,y_hat);
    new_res_outer = old_res_outer;

    new_res_inner=old_res_outer;
    old_res_inner=old_res_outer;

    while(nNuMethods<maxNumberInnerLoops){
      //while(nnuMethods<maxNumberInnerLoops){
      s_old=s;
      s.Resize(actNrParameter+actNrParameterC);
      //  z.Resize(nrMeasuredData);

      nNuMethods++;
      //        std::cout <<"\n Here starts the nuMethods Iteration - Nr " << nnuMethods<< std::endl;
      old_res_inner=new_res_inner;

      //   F' * s^k
      JacobiMatrix.Mult(s,JacFs);
        
      //       // F'sk - (y_hat-F_hat)
      //        for (UInt i=0;i<nrMeasuredData;i++){
      //      act_res[i]=y_hat[i]-F_hat[i];
      //      JacFs_res[i]=JacFs[i]-act_res[i];
      //        }
        
       
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
      //        std::cout<<"\n mu = " << eta_acc << ", omega = " << omega << ", nNuMehtods= " << nNuMethods << std::endl;
      //         std::cout<<z<<std::endl;
      adjJacobiMatrix.Mult(z,s);
      //        std::cout<<s<<std::endl;
      
         
      for(UInt i=0;i<nrMeasuredData;i++)
        z_old[nNuMethods][i]=z[i];     
         
      JacobiMatrix.Mult(s,JacFs);
        
      //F'(p^k)(s^k)-(y-F(p^k))
      for (UInt i=0;i<nrMeasuredData;i++){
        act_res[i]=y_hat[i]-F_hat[i];
        JacFs_res[i]=JacFs[i]-act_res[i];
      }

      //        for (UInt i=0;i<nrMeasuredData;i++)
      //act_res[i]-=JacFs[i];
      norm(JacFs_res, new_res_inner, max_res_inner,y_hat);
      //     new_res_inner=a2norm(JacFs_res);
      std::cout<<"\n new_res_inner = "<< new_res_inner <<", old_res_inner = " << old_res_inner<< std::endl;


        
      if (new_res_inner>1.15*old_res_inner){
        std::cout << " \n !! New_res_inner is worse than old_res_inner!! "<< std::endl;
        std::cout<<"\n Nr of nuMethodsC = " << nNuMethods <<std::endl;
        //      getchar();
        s=s_old;
        break;
      }
        
      std::cout<<"\n end of inner nuMethodsC Iter ..."<<std::endl;
    } // end while nuMethod ...


    nNuMethods=0;

//     Double old_resid2=old_res_outer;
//     Double new_resid2=new_res_outer;

    // backtracking(et , theta, s, old_resid2, new_resid2); 

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

    Double mult=1.0;
    //     if (nNuMethods<5)
    //       mult=1.5;
    //     else
    //       mult=1.0;

    scalingC[0]=1.0/((*matMatC)[0][0]); 
    scalingC[1]=1.0/((*matMatC)[2][2]);
    scalingC[2]=1.0/((*matMatC)[1][0]);
    scalingC[3]=1.0/((*matMatC)[0][2]);
    scalingC[4]=1.0/((*matMatC)[3][3]); 
    scalingC[5]=1.0/((*matMatC)[6][4]);
    scalingC[6]=std::abs(1.0/((*matMatC)[8][0]));
    scalingC[7]=mult/((*matMatC)[8][2]);
    scalingC[8]=1.0/((*matMatC)[6][6]); 
    scalingC[9]=mult/((*matMatC)[8][8]);

    for (UInt i=0;i<actNrParameter;i++)
      stepR[i]=s[i].real();

    for (UInt i=actNrParameter;i<actNrParameter+actNrParameterC;i++)
      //  stepC[i-actNrParameter]=s[i].real();
      stepC[i-actNrParameter]=s[i].imag();

    std::cout<<"stepR:" <<std::endl;
    std::cout<<stepR<<std::endl;


    std::cout<<"stepC:" <<std::endl;
    std::cout<<stepC<<std::endl;

    
    //    parameter_new=parameter;
    theta=1.0;
    Double thetaC=1000.0;
    setNewParameterSet(parameter, parameter, scaling, theta, stepR, whichParameterToUpdate);
    setNewParameterSet(parameterC, parameterC, scalingC, thetaC, stepC, whichParameterToUpdateC);

    // if no backtracking is specified, please include the following lines!
    for (UInt i=0;i<nrParameter;i++){
      //        parameter_new[i]=scaling[i]*parameter[i];
      //        parameter_new[i]+=s[i].real();
      //        parameter[i]=1/scaling[i]*parameter_new[i];
      std::cout<<" paramter("<<i<<") = " << parameter[i]<<" + " << parameterC[i] << " i " << std::endl;
    }
    // parameter=parameter_new;
      
    updateMaterialData(parameter, ptMaterial);
    updateComplexMaterialData(parameterC, ptMaterial);

    createF(ptMaterial, F_hat,FALSE);

    for (UInt i=0;i<y_hat.GetSize();i++)
      act_res[i]=y_hat[i]-F_hat[i];
   
    norm(act_res, new_res_outer, max_res_inner,y_hat);
    //      new_res_outer=(a2norm(act_res));

    std::cout<<"\n Norm of residual = " << new_res_outer <<std::endl;

    *parLog<<" "<< new_res_outer<<"  "; 


    //     if(new_res_outer<=1.0e-08)
    //       getchar();


    if (new_res_outer>old_res_outer){
      std::cout<<"\n Warning: residual norm gets worse!" <<std::endl;
      //        getchar();
      //parameter=parameter_new;
      //        NewtonNuMethods();
    }
    parameter_old=parameter;
      
   
  } // end NewtonNuMethodsy

}


