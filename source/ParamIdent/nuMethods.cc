// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "PDE/SinglePDE.hh"
#include "piezoParamIdent.hh"
#include "Utils/mathfunctions.hh"

namespace CoupledField {

// xxxxxxxxxxxxxxxxxxxxxxxxxx nuMethods  - ITERATION xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

void piezoParamIdent::nuMethods() {

  UInt nrIterations=0;
  UInt nNuMethods=0;
  Double theta, eta_acc, nu, omega;
  nu=1.5;

  updateMaterialData(parameter_); //Writes initial guesses of parameters (read from MeasuredData.dat) to system

  Double normJacMat, old_res_outer, res_outer, new_res_inner, old_res_inner,
      new_res_outer, maxres_inner;
  Double relax;
  relax=10.0;

  Vector<Complex> act_res(nrMeasuredData_);
  Vector<Complex> JacFs(nrMeasuredData_);
  Vector<Complex> z(nrMeasuredData_);
  Matrix<Complex> z_old(maxNumberInnerLoops_+2, nrMeasuredData_);
  Vector<Double> stepR(actNrParameter_);
  Vector<Complex> s_old(actNrParameter_);
  Vector<Complex> JacFs_res(nrMeasuredData_);
  Matrix<Complex> ImgSpaceScaling_Mat(nrMeasuredData_, nrMeasuredData_);
  Vector<Double> parameter_old(nrParameter_);

  updateMaterialData(parameter_);

  if (imagMaterialParam_ ) {
    updateComplexMaterialData(parameterC_);
  }

  createF(F_hat_, false);

  act_res = y_hat_-F_hat_;
  norm(act_res, new_res_outer, maxres_inner, y_hat_);

  std::cout<<"\n Norm of residual ||F(p)-y|| = "<< new_res_outer <<std::endl;
  res_outer=new_res_outer;

  *parLog_<<" "<< new_res_outer<<"  ";

  // *parLog<< new_res_outer;

  /*for (UInt i=0;i<nrParameter_;i++)
    if (whichParameterToUpdate_[i]==1)
      *parLog_<<"  " << parameter_[i];
    *parLog_<<"  " << newtonCounter_;
    *parLog_<<std::endl; */


  nrIterations++;
  s_.Resize(actNrParameter_);
  s_.Init();
  z.Resize(nrMeasuredData_);
  z.Init();
  z_old.Resize(maxNumberInnerLoops_+2, nrMeasuredData_);
  z_old.Init();
  parameter_old=parameter_;

  // Create the Matrices F, F', F*
  computeJacobiMatrix();

  JacobiMatrix_=approxJacobiMatrix_;

  for (UInt i=0; i<nrMeasuredData_; i++)
    for (UInt j=0; j<nrMeasuredData_; j++)
      if (i==j)
      {
        if (whichNormCriteria_=="logAmplitude")
          ImgSpaceScaling_Mat[i][j] = 1.0/std::log(std::abs(Complex(real_[i],
              imag_[i])));

        else if (whichNormCriteria_=="logImpedance")
          ImgSpaceScaling_Mat[i][j] = 1.0/std::log(std::abs(Complex(real_[i],
              imag_[i]))); ///(180.0/PI*std::atan2(imag_[i],real_[i]));

        else if (whichNormCriteria_=="amplitude")
          ImgSpaceScaling_Mat[i][j] = 1.0
              /(std::abs(Complex(real_[i], imag_[i]))); ///(180.0/PI*std::atan2(imag_[i],real_[i]));

        else if (whichNormCriteria_=="phase")
          ImgSpaceScaling_Mat[i][j] = 1.0/(180.0/PI*std::atan2(imag_[i],
              real_[i]));

        else if (whichNormCriteria_=="amplitudeMech")
          ImgSpaceScaling_Mat[i][j]=1000.0/(std::abs(Complex(realMech_[i],
              imagMech_[i])));

        else if (whichNormCriteria_=="logAmplitudeMech")
          ImgSpaceScaling_Mat[i][j]=-1.0/(std::log(std::abs(Complex(
              realMech_[i], imagMech_[i]))));

        else if (whichNormCriteria_=="displacementMech")
          ImgSpaceScaling_Mat[i][j]=1.0/(0.0001*std::abs((realMech_[i],
              imagMech_[i])));

        else if (whichNormCriteria_=="phaseMech")
          ImgSpaceScaling_Mat[i][j]=1.0/(180.0/PI * std::atan2(realMech_[i],
              imagMech_[i]));

        else
          std::cerr<<"Your choice of the fitting quantity seems to be invalid"
              <<std::endl;
      }

  createAdjointJacobiMatrix();
  Matrix<Complex> adjJacobiTemp;
  adjJacobiTemp = adjJacobiMatrix_*ImgSpaceScaling_Mat;
  adjJacobiMatrix_ = adjJacobiTemp;

  Double w=0.9;
  normJacMat=calcEuclidianMatrixNorm(JacobiMatrix_);

  while (w>=1/(normJacMat*normJacMat))
    w=0.9*w;

  JacobiMatrix_ *= Complex(w, w);
  adjJacobiMatrix_ *= Complex(w, w);

  //Norm ersetzt:
  norm(act_res, old_res_outer, maxres_inner, y_hat_);
  new_res_outer = old_res_outer;

  new_res_inner=old_res_outer;
  old_res_inner=old_res_outer;

  while (nNuMethods<maxNumberInnerLoops_) {
    s_old=s_;

    nNuMethods++;
    s_.Resize(actNrParameter_);
    s_.Init();
    old_res_inner=new_res_inner;

    // fitting am dicken mode
    relax=10.0;

    if (newtonCounter_>=25)
      relax=5.0;
    if (newtonCounter_>=30)
      relax=5.0;

    // fitting am radial mode
    //      relax=100.0;

    eta_acc = ((nNuMethods-1)*(2*nNuMethods-3)*(2*nNuMethods+2*nu_-1))
        /((nNuMethods+2*nu_-1)*(2*nNuMethods+4*nu_-1)*(2*nNuMethods+2*nu_-3));
    omega = 4*((2*nNuMethods + 2*nu_ -1)*(nNuMethods+nu_-1))/((nNuMethods+2*nu_-1)
        *(2*nNuMethods+4*nu_-1));

    for (UInt i=0; i<nrMeasuredData_; i++) {
      if (nNuMethods>1) {
        z[i]=z[i]+relax*(eta_acc*(z[i]-z_old[nNuMethods-2][i])+omega
            *(act_res[i]-JacFs[i]));
      } else {
        z[i]=relax*(4*nu_+2)/(4*nu_+1)*act_res[i];
        //     std::cout<<"s("<<i<<")= "<<s[i]<<"; "<<std::endl;
      }
    }

    adjJacobiMatrix_.Mult(z, s_);

    for (UInt i=0; i<nrMeasuredData_; i++)
      z_old[nNuMethods][i]=z[i];

    JacobiMatrix_.Mult(s_, JacFs);

    //F'(p^k)(s^k)-(y-F(p^k))
    for (UInt i=0; i<nrMeasuredData_; i++) {
      act_res[i]=y_hat_[i]-F_hat_[i];
      JacFs_res[i]=JacFs[i]-act_res[i];
    }

    norm(JacFs_res, new_res_inner, maxres_inner, y_hat_);

    if (new_res_inner>1.15*old_res_inner) {
      std::cout<<"\n "<< nNuMethods <<" inner iterations during Newton-step  "
          << newtonCounter_ << std::endl;

      // The following part prevents stagnation of algorithm in case that
      // no better update was found during the first iteration step.

      Complex sum=Complex(0.0, 0.0);
      for (UInt i=0; i<s_old.GetSize(); i++)
        sum+=s_old[i];
      if (std::abs(sum)!=0.0)
        s_=s_old;
      else
        for (UInt i=0; i<s_old.GetSize(); i++)
          s_[i]=0.1*s_[i];

      break;
    }

  } // end while nuMethod ...


  nNuMethods=0;
  computeScaling();

  for (UInt i=0; i<actNrParameter_; i++)
    stepR[i]=s_[i].real();

  theta=5.0;

  parameter_old=parameter_;
  setNewParameterSet(parameter_, parameter_, scaling_, theta, stepR,
      whichParameterToUpdate_);
  updateMaterialData(parameter_);
  createF(F_hat_, false);

  for (UInt i=0; i<nrMeasuredData_; i++)
    act_res[i]=y_hat_[i]-F_hat_[i];

  norm(act_res, new_res_outer, maxres_inner, y_hat_);

  Integer lineSearchCount=0;

  std::cout<<"++ Performing line search ... "<<std::endl;

  bool negFlag=false;

   for (UInt par=0; par<nrParameter_; par++)
        if (parameter_[par]<0.0&&par!=6) {
          negFlag=true;
        }

   while (new_res_outer>old_res_outer || negFlag==true) {
     theta = 0.5*theta;
     negFlag=false;
    //std::cout<<"theta = "<<theta<<std::endl;
    parameter_=parameter_old;
    setNewParameterSet(parameter_, parameter_, scaling_, theta, stepR,
        whichParameterToUpdate_);
    updateMaterialData(parameter_);
    createF(F_hat_, false);

    for (UInt i=0; i<nrMeasuredData_; i++)
      act_res[i]=y_hat_[i]-F_hat_[i];

    for (UInt par=0; par<nrParameter_; par++)
          if (parameter_[par]<0.0&&par!=6)
            negFlag=true;

    norm(act_res, new_res_outer, maxres_inner, y_hat_);
    //std::cout<<"new_res_outer = " << new_res_outer <<std::endl;
    //std::cout<<"res_outer = " << res_outer <<std::endl;

    lineSearchCount++;
    if (lineSearchCount>35) {
      std::cout
          <<"! TERMINATION of line search after 35 backtracking steps (with step width = "
          << theta <<")\n ... which is a suspicious event ... "<<std::endl;
      std::cout
          <<"Either we are close to a solution or something went completely wrong."
          <<std::endl;
      std::cout<<"Check your results (convergence?) and input data!"<<std::endl;
      nu_ = 1.0 + (nu_-1.0)*0.5;
      std::cout<<"nu_: "<<nu_<<std::endl;
      break;
    }

  }

  *parLog_<<" "<< new_res_outer<<"  ";

  for (UInt i=0; i<nrParameter_; i++)
    if (whichParameterToUpdate_[i]==1)
      *parLog_<<"  "<< parameter_[i];

  for (UInt i=0; i<nrParameter_; i++)
    if (whichParameterToUpdateC_[i]==1)
      *parLog_<<"  "<< parameterC_[i];

  *parLog_<<std::endl;

  std::cout<<"New parameter:"<< std::endl;
  // if no backtracking is specified, please include the following lines!
  for (UInt i=0; i<nrParameter_; i++) {
    //      parameter_new_[i]=scaling_[i]*parameter_[i];
    //      parameter_new_[i]+=s[i].real();
    //      parameter_[i]=1/scaling_[i]*parameter_new_[i];
    if (whichParameterToUpdate_[i]==1)
      std::cout<<"Parameter("<<i<<") = "<< parameter_[i]<<" ( ~ "<< 100-(1
          -parameter_[i]/parameterIncrement_[i])*100
          <<" percent of initial guess) "<< std::endl;
  }
  // parameter_=parameter_new_;

  updateMaterialData(parameter_);
  createF(F_hat_, false);

  for (UInt i=0; i<nrMeasuredData_; i++)
    act_res[i]=y_hat_[i]-F_hat_[i];
  //Norm ersetzt:
  //      std::cout<<act_res<<std::endl;
  norm(act_res, new_res_outer, maxres_inner, y_hat_);
  parameter_old=parameter_;
  residuumParIdent_=new_res_outer;

} // end NewtonNuMethods


// XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//
// ------------------------ nu METHODS-C II --------------------------------------------------------------------
//
// XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

void piezoParamIdent::nuMethodsC2() {

  UInt nrIterations=0;
  UInt nNuMethods=0;
  Double theta, eta_acc, omega;

  updateMaterialData(parameter_); //Writes initial guesses of parameters (read from MeasuredData.dat) to system
  updateComplexMaterialData(parameterC_); //Writes initial guesses of parameterC_

  Double normJacMat, max_res_inner, maxres_inner, old_res_outer, new_res_inner,
      old_res_inner, new_res_outer;

  Vector<Complex> act_res(nrMeasuredData_);
  Vector<Complex> JacFs(nrMeasuredData_);
  Vector<Complex> z(nrMeasuredData_);
  Matrix<Complex> z_old(maxNumberInnerLoops_+2, nrMeasuredData_);
  Vector<Double> stepR(actNrParameter_);
  Vector<Complex> s_old(actNrParameter_+actNrParameterC_);
  Vector<Double> stepC(actNrParameterC_);
  Vector<Complex> JacFs_res(nrMeasuredData_);
  Matrix<Complex> ImgSpaceScaling_Mat(nrMeasuredData_, nrMeasuredData_);
  Vector<Double> parameter_old(nrParameter_);
  Vector<Double> parameter_oldC(nrParameter_);

  createF(F_hat_, false);
  act_res = y_hat_-F_hat_;
  //    new_res_outer=old_res_outer=a2norm(act_res);
  norm(act_res, new_res_outer, maxres_inner, y_hat_);
  std::cout<<"\n||F(p)-y|| = "<< new_res_outer <<std::endl;

  nrIterations++;
  // std::cout<<"\n Newton NuMethodsC ... Newton-Iteration-Nr = "<<nrIterations<<std::endl;
  s_.Resize(actNrParameter_+actNrParameterC_);
  s_.Init();
  z.Resize(nrMeasuredData_);
  z.Init();
  z_old.Resize(maxNumberInnerLoops_+2, nrMeasuredData_);
  z_old.Init();
  parameter_old = parameter_;
  parameter_oldC = parameterC_;

  // Create the Matrices F, F', F*
  createF(F_hat_, false);
  computeJacobiMatrixC();
  JacobiMatrix_=approxJacobiMatrix_;

  for (UInt i=0; i<nrMeasuredData_; i++)
    for (UInt j=0; j<nrMeasuredData_; j++)
      if (i==j) {
        if (whichNormCriteria_=="logAmplitude")
          ImgSpaceScaling_Mat[i][j] = 1.0/std::log(std::abs(Complex(real_[i],
              imag_[i])));
        else if (whichNormCriteria_=="logImpedance")
          ImgSpaceScaling_Mat[i][j] = 1.0/std::log(std::abs(Complex(real_[i],
              imag_[i])));
        else if (whichNormCriteria_=="amplitude")
          ImgSpaceScaling_Mat[i][j] = 1.0
              /(std::abs(Complex(real_[i], imag_[i])));
        else if (whichNormCriteria_=="phase")
          ImgSpaceScaling_Mat[i][j] = 1.0/(180.0/PI*std::atan2(imag_[i],
              real_[i]));
        else if (whichNormCriteria_=="amplitudeMech")
          ImgSpaceScaling_Mat[i][j]=1000.0/(std::abs(Complex(realMech_[i],
              imagMech_[i])));

        else if (whichNormCriteria_=="logAmplitudeMech")
          ImgSpaceScaling_Mat[i][j]=-1.0/(std::log(std::abs(Complex(
              realMech_[i], imagMech_[i]))));

        else if (whichNormCriteria_=="displacementMech")
          ImgSpaceScaling_Mat[i][j]=1.0/(0.0001*std::abs((realMech_[i],
              imagMech_[i])));

        else if (whichNormCriteria_=="phaseMech")
          ImgSpaceScaling_Mat[i][j]=1.0/(180.0/PI * std::atan2(realMech_[i],
              imagMech_[i]));

        else
          std::cerr<<"Your choice of the fitting quantity seems to be invalid"
              <<std::endl;
        // ImgSpaceScaling_Mat[i][j]=1.0/std::log(Complex(real_[i],imag_[i]));
      }

  createAdjointJacobiMatrix();

  Matrix<Complex> adjJacobiTemp;
  adjJacobiTemp = adjJacobiMatrix_*ImgSpaceScaling_Mat;
  adjJacobiMatrix_ = adjJacobiTemp;

  Double w=0.9;

  normJacMat=calcEuclidianMatrixNorm(JacobiMatrix_);

  while (w>=1/(normJacMat*normJacMat))
    w=0.9*w;

  // if nu-method ..
  //JacobiMatrix = JacobiMatrix*Complex(w,w);
  JacobiMatrix_ *= Complex(w, w);
  //adjJacobiMatrix = adjJacobiMatrix*Complex(w,w);
  adjJacobiMatrix_ *= Complex(w, w);

  norm(act_res, old_res_outer, maxres_inner, y_hat_);
  new_res_outer = old_res_outer;

  new_res_inner=old_res_outer;
  old_res_inner=old_res_outer;

  *parLog_<<" "<< new_res_outer<<"  ";

  for (UInt i=0; i<nrParameter_; i++)
    if (whichParameterToUpdate_[i]==1)
      *parLog_<<"  "<< parameter_[i];

  for (UInt i=0; i<nrParameter_; i++)
    if (whichParameterToUpdateC_[i]==1)
      *parLog_<<"  "<< parameterC_[i];

  *parLog_<<std::endl;

  while (nNuMethods<maxNumberInnerLoops_) {
    //while(nnuMethods<maxNumberInnerLoops_){
    s_old=s_;
    s_.Resize(actNrParameter_+actNrParameterC_);
    s_.Init();
    //  z.Resize(nrMeasuredData_);

    nNuMethods++;
    //        std::cout <<"\n Here starts the nuMethods Iteration - Nr " << nnuMethods<< std::endl;
    old_res_inner=new_res_inner;

    //   F' * s^k
    JacobiMatrix_.Mult(s_, JacFs);

    Double relax=10.0;
    if (newtonCounter_>=25)
      relax=5.0;
    else if (newtonCounter_>=30)
      relax=5.0;

    eta_acc = ((nNuMethods-1)*(2*nNuMethods-3)*(2*nNuMethods+2*nu_-1))
        /((nNuMethods+2*nu_-1)*(2*nNuMethods+4*nu_-1)*(2*nNuMethods+2*nu_-3));
    omega = 4*((2*nNuMethods + 2*nu_ -1)*(nNuMethods+nu_-1))/((nNuMethods+2*nu_-1)
        *(2*nNuMethods+4*nu_-1));

    for (UInt i=0; i<nrMeasuredData_; i++) {
      if (nNuMethods>1) {
        z[i]=z[i]+relax*(eta_acc*(z[i]-z_old[nNuMethods-2][i])+omega
            *(act_res[i]-JacFs[i]));
      } else {
        z[i]=relax*(4*nu_+2)/(4*nu_+1)*act_res[i];
        //     std::cout<<"s("<<i<<")= "<<s[i]<<"; "<<std::endl;
      }
    }

    adjJacobiMatrix_.Mult(z, s_);

    for (UInt i=0; i<nrMeasuredData_; i++)
      z_old[nNuMethods][i]=z[i];

    JacobiMatrix_.Mult(s_, JacFs);

    //F'(p^k)(s^k)-(y-F(p^k))
    for (UInt i=0; i<nrMeasuredData_; i++) {
      act_res[i]=y_hat_[i]-F_hat_[i];
      JacFs_res[i]=JacFs[i]-act_res[i];
    }

    norm(JacFs_res, new_res_inner, max_res_inner, y_hat_);

    if (new_res_inner>1.05*old_res_inner) {
    //if (new_res_inner>1.15*old_res_inner) {
      //      std::cout << " \n !! New_res_inner is worse than old_res_inner -> break of inner Loop! "<< std::endl;
      std::cout<<"\n Computed "<< nNuMethods
          <<" inner iterations during Newton-step  "<< newtonCounter_ <<"\n"
          << std::endl;
      //      getchar();

      // The following part prevents stagnation of algorithm in case that
      // no better update was found during the first iteration step.

      Complex sum=Complex(0.0, 0.0);
      for (UInt i=0; i<s_old.GetSize(); i++)
        sum+=s_old[i];
      if (std::abs(sum)!=0.0)
        s_=s_old;
      else
        for (UInt i=0; i<s_old.GetSize(); i++)
          s_[i]=0.1*s_[i];

      break;
    }
    //       }

  } // end while nuMethod ...


  nNuMethods=0;

  //     Double old_resid2=old_res_outer;
  //     Double new_resid2=new_res_outer;

  // backtracking(et , theta, s, old_resid2, new_resid2);

  computeScaling();

  // Double mult=1.0; // TODO: Check if this is still needed

  for (UInt i=0; i<actNrParameter_; i++)
    stepR[i]=s_[i].real();

  //  std::cout<<actNrParameter_<<std::endl;
  //  std::cout<<actNrParameterC_<<std::endl;
  //  std::cout<<s<<std::endl;

  for (UInt i=actNrParameter_; i<actNrParameter_+actNrParameterC_; i++)
    stepC[i-actNrParameter_]= -1000.0*s_[i].real();

    //    parameter_new_=parameter;
  theta=1.0;
  Double thetaC=1.0;
  setNewParameterSet(parameter_, parameter_, scaling_, theta, stepR,
      whichParameterToUpdate_);
  setNewParameterSet(parameterC_, parameterC_, scalingC_, thetaC, stepC,
      whichParameterToUpdateC_);

  // if no backtracking is specified, please include the following lines!

  // parameter=parameter_new_;

  updateMaterialData(parameter_);
  updateComplexMaterialData(parameterC_);

  createF(F_hat_, false);

  for (UInt i=0; i<nrMeasuredData_; i++)
    act_res[i]=y_hat_[i]-F_hat_[i];

  norm(act_res, new_res_outer, max_res_inner, y_hat_);
  //      new_res_outer=(a2norm(act_res));

  // std::cout<<"new_res_outer:"<<std::endl;
  //std::cout<<new_res_outer<<std::endl;

  Integer lineSearchCount=0;

  bool negFlag=false;

  for (UInt par=0; par<nrParameter_; par++)
       if (parameter_[par]<0.0&&par!=6) {
         negFlag=true;
       }


  std::cout<<"Performing line search ... "<< std::endl;

  while (new_res_outer>old_res_outer || negFlag==true) {
    negFlag=false;
    theta = 0.5*theta;
    thetaC = 0.5*thetaC;
    //std::cout<<"theta = "<<theta<<std::endl;
    parameter_=parameter_old;
    parameterC_=parameter_oldC;

    setNewParameterSet(parameter_, parameter_, scaling_, theta, stepR,
        whichParameterToUpdate_);
    setNewParameterSet(parameterC_, parameterC_, scalingC_, thetaC, stepC,
        whichParameterToUpdateC_);

    updateMaterialData(parameter_);
    updateComplexMaterialData(parameterC_);

    createF(F_hat_, false);

    for (UInt i=0; i<nrMeasuredData_; i++)
      act_res[i]=y_hat_[i]-F_hat_[i];

    for (UInt par=0; par<nrParameter_; par++)
      if (parameter_[par]<0.0&&par!=6)
        negFlag=true;

    norm(act_res, new_res_outer, maxres_inner, y_hat_);
    // std::cout<<"new_res_outer = " << new_res_outer <<std::endl;

    lineSearchCount++;
    if (lineSearchCount>35) {
      std::cout<<"! TERMINATION of line search after 30 backtracking steps (with stepwidth = " << theta <<")"
          "\n ... which is a suspicious event ... " <<std::endl;
      std::cout<<"Either we are close to a solution or something went completely wrong." <<std::endl;
      std::cout<<"Check your results (convergence?) and input data!"<<std::endl;
      nu_ = 1.0 + (nu_-1.0)*0.5;
      std::cout<<"Try smaller nu: "<<nu_<<std::endl;
      break;
    }

  }

  for (UInt i=0; i<nrParameter_; i++) {
      std::cout<<" paramter("<<i<<") = "<< parameter_[i]<<" + "<< parameterC_[i]
          << " i "<< std::endl;
    }

  parameter_old=parameter_;
  parameter_oldC=parameterC_;
  residuumParIdent_=new_res_outer;

} // end NewtonNuMethods

}

