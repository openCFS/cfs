#include <fstream>
#include <iostream>
#include <string>

#include "multiHarmonicDriver.hh"
#include "stdSolveStep.hh"

#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "Domain/domain.hh"
#include "PDE/StdPDE.hh"
#include "DataInOut/MHMaterialDataFile.hh"


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



namespace CoupledField
{

  // ========================================================================
  // ========================= piezoParamIdent - Part ===========================
  // ========================================================================

  // constructor
  // opens datafiles: measuredData.dat for input, imedCurve.dat and piezoLog.dat for output

  MultiHarmonicDriver::MultiHarmonicDriver(Domain * adomain,
                                           UInt stepOffset,
                                           Double timeOffset,
                                           std::string driverTag,
                                           Boolean isPartOfSequence)
    :SingleDriver(adomain, stepOffset, timeOffset, 
                  driverTag, isPartOfSequence){

    ENTER_FCN( "multiharmonic::multiharmonic" );

    //  ptDomain = adomain;
    ptMyPDE_ = NULL;

    StdVector<std::string> keyVec, attrVec, valVec;
    
    attrVec = "tag";
    valVec = driverTag_;
    
    // Get time stepping information from parameter object
    keyVec = "multiHarmonic", "startFreq";
    params->Get(keyVec, attrVec, valVec, startFreq_);
  
    keyVec = "multiHarmonic", "stopFreq";
    params->Get(keyVec, attrVec, valVec, stopFreq_);
    
    keyVec = "multiHarmonic", "numFreq";
    params->Get(keyVec, attrVec, valVec, numFreq_);

    keyVec = "multiHarmonic", "nrMultiHarmonics";
    params->Get(keyVec, attrVec, valVec, nrMultHarms_);

 
    Char* measuredData="measuredData.dat";
    allMeasuredData = new std::ifstream(measuredData, std::basic_ios<char>::in);

    if (!allMeasuredData)
      {
        std::cerr << "\n File measuredData.dat does not exist!" << std::endl;
        exit(1);
      }

    std::cout<<"\n Opens impedCurve.dat and piezoLog.dat ... "<<std::endl;

    std::string filename= "imped.dat";

    impedCurve = new std::ofstream(filename.c_str(),std::basic_ios<char>::out);

    if (!impedCurve)
      {
        std::cerr << "\n ImpedanceCurve.dat could not be initialized" << std::endl;
      }

    std::string filenameLog= "piezoLog.dat";

    piezoLog = new std::ofstream(filenameLog.c_str(),std::basic_ios<char>::out);
    
    if (!piezoLog)
      {
        std::cerr << "\n piezoLog.dat could not be initialized" << std::endl;
      }

    std::string filenameParLog= "parLog.dat";
    parLog = new std::ofstream(filenameParLog.c_str(),std::basic_ios<char>::out);

    if (!parLog)
      {
        std::cerr << "\n piezoLog.dat could not be initialized" << std::endl;
      }

  } // end of constructor

  // destructor
  
  MultiHarmonicDriver::~MultiHarmonicDriver(){
    ENTER_FCN( "MultiHarmonicDriver::~MultiHarmonicDriver" );
    allMeasuredData->close();
    impedCurve->close();
    piezoLog->close();
    parLog->close();
  }

  void MultiHarmonicDriver::SolveProblem() {
    ENTER_FCN( "MultiHarmonicDriver::SolveProblem" );
    
    if (! isPartOfSequence_)
      ptdomain_->PrintGrid();
    if (isPartOfSequence_ == FALSE){     
      GetMyPDEs();
    
      //! cast pointer to BasePDE * to pointer of SinglePDE *
      ptMyPDE_ = dynamic_cast<SinglePDE*>(ptPDE_);
      Info->StartProgress ("Starting to solve problem", FALSE);
    }
    ptMyPDE_->WriteGeneralPDEdefines();
    //    pdeId_ = ptMyPDE_->GetPDEId();


    //    MHAssembleMatrices();

    MHMaterialDataFile *ptMHFiles_;
    ptMHFiles_ = new MHMaterialDataFile;
    ptMHFiles_->computeIndexSet(2,2,nrMultHarms_);

    UInt maxP=5;
    
    parameterCoeff_.Resize(10,maxP);
    parameter_.Resize(10);

    parameterCoeff_[0][0]=1.0; // c_11, a_0=1, andere a_i=0, somit konstant
    parameterCoeff_[1][0]=1.0; // c_33, a_0=1
    parameterCoeff_[1][1]=0.5; // c_33, a_1=0.5
    parameterCoeff_[1][2]=0.25; // c_33, a_2=0.25 Polynom vom Grad 3
    parameterCoeff_[2][0]=1.0; // c_12, a_0=1, andere a_i=0, somit konstant
    parameterCoeff_[3][0]=1.0; // c_13, a_0=1, andere a_i=0, somit konstant
    parameterCoeff_[4][0]=1.0; // c_44, a_0=1, andere a_i=0, somit konstant
    parameterCoeff_[5][0]=1.0; // e_15, a_0=1, andere a_i=0, somit konstant
    parameterCoeff_[6][0]=1.0; // e_13, a_0=1, andere a_i=0, somit konstant
    parameterCoeff_[7][0]=1.0; // e_33
    parameterCoeff_[7][1]=0.5; // e_33 Polynom vom Grad 2
    parameterCoeff_[8][0]=1.0; // eps_11, a_0=1, andere a_i=0, somit konstant
    parameterCoeff_[9][0]=1.0; // eps_33, 
    parameterCoeff_[9][1]=1.0; // eps_33, Polynom vom Grad 2


    UInt element;
    element=5;
    Integer delta=-1;

    calcParameterCurveAtElement(parameter_, parameterCoeff_, element, nrMultHarms_, delta, maxP);
    

    Matrix<Integer> exps;
    Matrix<Integer> count;
    ptMHFiles_->getExponentArray(exps,count,3,3,-3);    
    std::cout<<"delta = -3: " <<std::endl;
    std::cout<< exps << std::endl;
    std::cout<< count << std::endl;

    Boolean reset = TRUE;
    
    if (! isPartOfSequence_)
      ptdomain_->PrintGrid();

    // if driver is not part of multiSequence Driver, get list
    // of pdes which have to be solved and intialize them
    if (isPartOfSequence_ == FALSE){     
      GetMyPDEs();
      Info->StartProgress ("Starting to solve problem", FALSE);
    }
  
    ptPDE_->WriteGeneralPDEdefines();
      
    UInt fstep;
    Double actFreq  = startFreq_;
    Double freqIncr = (stopFreq_ - startFreq_) / numFreq_;
    std::string errMsg;
  
    // branch for single PDE
    for (fstep = 1; fstep <= numFreq_; fstep++) {
      Info->WriteHarmonicStep(ptPDE_->GetName(), fstep, actFreq);
    
      ptPDE_->GetSolveStep()->SetActFreq(actFreq);
      ptPDE_->GetSolveStep()->SetActStep(fstep);

      //      std::cout<<"\n multiHarm: 1 " <<std::endl;
      ptPDE_->GetSolveStep()->PreStepHarmonic(reset);
    
      //      std::cout<<"\n multiHarm: 2"  <<std::endl;
      ptPDE_->GetSolveStep()->SolveStepHarmonic(reset);
    
      //      std::cout<<"\n multiHarm: 3 " <<std::endl;
      ptPDE_->GetSolveStep()->PostStepHarmonic(reset);
    
      ptPDE_->PostProcess();
      ptPDE_->WriteResultsInFile( fstep, actFreq);
    
      actFreq += freqIncr;
    }
  }

  void MultiHarmonicDriver::MHAssembleMatrices(){
    ENTER_FCN( "MultiHarmonicDriver::MHAssembleMatrices" );

    // StdVector<RegionIdType> MHsubdoms_;
//     Assemble * ptAssemble_;
//     Matrix<Double> elemmat;

//     MHsubdoms_ = ptMyPDE_->getPDE_subdoms();
//     ptAssemble_=ptMyPDE_->getPDE_assemble();

//     Vector<Double> harmonicVec;
//     Double dampTransform = 1.0;
    
//     // initialize reassembling "indicator" vector
//     for (Integer MHInd=-nrMultHarms_; MHInd<=nrMultHarms_;MHInd++){

//     if (ptAssemble_->firstTime_)
//        for (UInt actMat=0; actMat < ptAssemble_->nrMatrices_; actMat++)
//          ptAssemble_->reassembleMat_[actMat] = TRUE;

//     std::cout<<"MultiHarmonicDriver::MHAssembleMatrices Nr "<< MHInd <<std::endl;

//     for (UInt actDom=0; actDom < MHsubdoms_.GetSize(); actDom++) {     
//       StdVector<Elem*> elemssd;

//        ptAssemble_->ptgrid_->GetVolElems(elemssd, MHsubdoms_[actDom]);

//        for(UInt actInteg=0; actInteg < ptAssemble_->integrators_[actDom]->GetSize(); actInteg++) {

//          std::cout<<" Integrator "<< actInteg <<std::endl;
       

//          IntegratorDescriptor * actDescriptor = (*ptAssemble_->integrators_[actDom])[actInteg];

//          actDescriptor->GetIntegrator()->SetSubdomain(actDom);

//          if (ptAssemble_->alternateMaterialData_ == TRUE)
//            actDescriptor->GetIntegrator()->SetMaterial(ptAssemble_->ptMaterial_);
      
                    
//          // assemble only if nonlinear or first time
//          if (ptAssemble_->reassembleMat_[actDescriptor->DestMat()] || ptAssemble_->firstTime_) {

//            if ( actDescriptor->IsReducedInt()) {
//              //set all elements to reduced integration!!
//              ptAssemble_->SetFE2ReducedInt();
// //           }

// //              if ( actDescriptor->GetIntegrator()->IsFracDamping() ) {
// // //             // get multiplicative pre factor depending on time step size
// //                Double damp;
// //                damp = ptMyPDE_->GetFracDampMatrixCoeff(actDom);
// //                actDescriptor->GetIntegrator()->SetFactor(damp);
// //              }
//            }
//          }
       

// //           dampTransform = 1.0;
//        if ( ( actDescriptor->GetIntegrator()->IsRaylDamping() 
//               || actDescriptor->GetSecondaryMat() != NOTYPE ) 
//             && ptAssemble_->startFreq_ > 0 ) {

//              // Obtain frequency value to which the damping parameters
//              // in the material file do belong
//              StdVector<Double> freqs;
//              Double matDataFreq;
//              params->GetList( "matDataFreq", freqs, "harmonic" );
//              if ( freqs.GetSize() == 1 ) {
//                matDataFreq = freqs[0] * 2.0 * PI;
//              }
//              else {
//                matDataFreq = ptAssemble_->startFreq_;
//              }

             
// //             // get multiplicative pre factor depending on frequency
//              if ( matDataFreq > 0 && actFreq_ > 0 ) {
//                FEMatrixType destMat =
//                  actDescriptor->GetIntegrator()->GetBaseType();

//                if ( destMat == STIFFNESS ) {
//                  dampTransform = matDataFreq / actFreq_;
//                  Info->PrintF( "", " dampTransform (stiffness matrix) = %e\n",
//                                dampTransform );
//                }
//                else if ( destMat == MASS ) {
//                  dampTransform = actFreq_ / matDataFreq;
//                  Info->PrintF( "", " dampTransform (mass matrix) = %e\n",
//                                dampTransform );
//                }
//              }

//              if ( actDescriptor->GetIntegrator()->IsRaylDamping() ) {
//                actDescriptor->GetIntegrator()->SetFactor(dampTransform);
//              }
//        }
          
       
//            //put pointer to array containing the material parameter 
//            // for each element
//        //actDescriptor->GetIntegrator()->SetMaterialArray(matArray_);

//        // HIER SELBER ÜBERLEGEN; WAS MIT DEM MATERIAL PASSIERT ....

//            for (UInt actEl=0; actEl< elemssd.GetSize(); actEl++) {
//              actDescriptor->GetIntegrator()->SetElemNr(actEl);

//              BaseFE * ptEl = elemssd[actEl]->ptElem;
//              StdVector<UInt> connecth = elemssd[actEl]->connect;
         
//              Matrix<Double> ptCoord;
//              ptAssemble_->GetElemCoords(connecth, ptCoord);
                    
// //             // map connect to PDE node numbers
//              StdVector<Integer> connect_PDE1, connect_PDE2;
                    
//              ptAssemble_->ptEQN1_->Node2EQN(connecth, connect_PDE1);
//              ptAssemble_->ptEQN2_->Node2EQN(connecth, connect_PDE2);
                    
//              Matrix<Double> elSol;
                    
//              actDescriptor->GetIntegrator()->SetElemPtr(ptEl);
//              FEMatrixType destMat = actDescriptor->DestMat();
           
                
// //             // this matrix is nonlinear and, therefore, 
// //             // has to be reassembled next time
//              if (actDescriptor->IsNonLin()) {
//                ptAssemble_->oneIntIsNonlin_ = TRUE;
//                ptAssemble_->reassembleMat_[actDescriptor->DestMat()] = TRUE;
//                ptAssemble_->sol_->GetElemSolutionAsMatrix(elSol, connecth);
//                actDescriptor->GetIntegrator()->SetActElemSol(elSol);
//              }       
           
                   
// //             // ================================================================
// //             //                             assemble matrices
// //             // ================================================================

//              actDescriptor->GetIntegrator()->
//                CalcElementMatrix(ptCoord, elemmat);
//              //             std::cout<<elemmat<<std::endl;
                  
//              piezoMaterialType matType = actDescriptor->GetPiezoMaterialType();
//              actDescriptor->SetPiezoMaterialType(matType);


//              if (ptAssemble_->analysisType_ == HARMONIC) {
//                ptAssemble_->TransformMatrix2Harmonic(harmonicVec,elemmat, 
//                                         actDescriptor->GetOrigMatrixType(),
//                                         actDescriptor->GetPiezoMaterialType());

//                if (destMat== MASS){
//                  elemmat = elemmat*MHInd;
//                  std::cout<<elemmat<<std::endl;
//                  getchar();

//                  ptAssemble_->algsys_->SetElementMatrix( destMat, &harmonicVec[0], 
//                                                          ptAssemble_->pdeId1_, connect_PDE1.GetPointer(), 
//                                                          connect_PDE1.GetSize(),
//                                                          ptAssemble_->pdeId2_, connect_PDE2.GetPointer(), 
//                                                          connect_PDE2.GetSize() );
//                }else 
//                  ptAssemble_->algsys_->SetElementMatrix( destMat, &harmonicVec[0], 
//                                                          ptAssemble_->pdeId1_, connect_PDE1.GetPointer(), 
//                                                          connect_PDE1.GetSize(),
//                                                          ptAssemble_->pdeId2_, connect_PDE2.GetPointer(), 
//                                                          connect_PDE2.GetSize() );

//              }
//            }
//        }
// //             else {

// //               algsys_->SetElementMatrix( destMat, elemmat.GetDataPointer(), 
// //                                          pdeId1_, connect_PDE1.GetPointer(), 
// //                                          connect_PDE1.GetSize(), 
// //                                          pdeId2_, connect_PDE2.GetPointer(), 
// //                                          connect_PDE2.GetSize() );
// //             }
// // #ifdef DEBUG
// //             // output matrices
// //             if (destMat == STIFFNESS) {
// //               (*debug) << "Stiffness matrix of Element " 
// //                        << actEl << std::endl;
// //             }

// //             if (destMat == MASS) {
// //               (*debug) << "Mass      matrix of Element " 
// //                        << actEl << std::endl;
// //             }

// //             if (destMat == DAMPING) {
// //               (*debug) << "Damping   matrix of Element " 
// //                        << actEl << std::endl;
// //             }

// //             if (destMat == SYSTEM) {
// //               (*debug) << "System matrix of Element " 
// //                        << actEl << std::endl;
// //             }
            
// //             (*debug) << elemmat << std::endl;

// //             if ( !elemmat.IsSymmetric() ) {
// //               (*debug) << " --> Matrix is not symmetric " 
// //                        << std::endl << std::endl;
// //             }
// //             else {
// //               (*debug) << " --> Matrix is symmetric " 
// //                        << std::endl << std::endl;
// //             }
// // #endif
// //             if (actDescriptor->GetSecondaryMat() != NOTYPE) {
// //               Double damp = dampTransform * actDescriptor->GetSecMatFac();
// //               elemmat *= damp;
// //               if (analysisType_ == HARMONIC) {
// //                 TransformMatrix2Harmonic(harmonicVec,elemmat,
// //                                          actDescriptor->GetOrigSecMatrixType(),
// //                                          actDescriptor->GetPiezoMaterialType());
        
// //                 algsys_->SetElementMatrix(destMat, &harmonicVec[0], 
// //                                           pdeId1_, connect_PDE1.GetPointer(),
// //                                           connect_PDE1.GetSize(), 
// //                                           pdeId2_, connect_PDE2.GetPointer(),
// //                                           connect_PDE2.GetSize() );
// //               }
// //               else 
// //                 algsys_->SetElementMatrix(actDescriptor->GetSecondaryMat(), 
// //                                           elemmat.GetDataPointer(), 
// //                                           pdeId1_, connect_PDE1.GetPointer(), 
// //                                           connect_PDE1.GetSize(), 
// //                                           pdeId2_, connect_PDE2.GetPointer(), 
// //                                           connect_PDE2.GetSize());
// //             }
                  
// //           } //over all elements of subdomain            
                
// //         } //check, if we have to assemble
          
// //         if ( actDescriptor->IsReducedInt()) {
// //           //set all elements back to standard integration!!
// //           SetFE2StandardInt();
// //         }
          

// //       } //integrators
        
  
    //   } //subdomains
    //  }// multiharmonics 
  }// MHAssemble Matrices
    


  void MultiHarmonicDriver::updateMaterialData(Vector<Double> & parameter, MaterialData * ptMaterial){
    ENTER_FCN( "multiharmonic::updateMaterialData");
  }



  void MultiHarmonicDriver::calcParameterCurveAtElement(Vector<Complex> & parameter, 
                                                        Matrix<Double> & parameterCoeff_, UInt element,Integer N,
                                                        Integer delta, UInt PP){

    UInt j;
    Complex prod;
    prod=Complex(1,0);
    UInt p;
    Double pFac;
    Vector<UInt> exponent;
    Double binomCoeffNenner;
    Complex innerSum;
    innerSum=Complex(0,0);
    //    getExponentArray(exponent, N, PP, delta);

    Complex elementSolution;
    elementSolution=Complex(0.1,0.5);

    for (UInt i=0;i<parameter.GetSize();i++){
      //      std::cout<<" Parameter Nr = " << i <<std::endl;
      j=0;
      p=1;
      //      std::cout<<parameterCoeff_[i][j]<<std::endl;
      while(parameterCoeff_[i][j]!=0){

         binomCoeffNenner=1;
         for (UInt e=0;e<exponent.GetSize();e++)
           if (exponent[e]!=0)
             binomCoeffNenner*=exponent[e];

         pFac=1;
         for (UInt pInd=1; pInd<=p; pInd++)
           pFac*=pInd;
         
         //         getExponentArray(exponent, N, p, delta);

          for (Integer n=-N;n<=N;n++)
            prod=prod*std::pow(elementSolution,exponent[n]);
        
          prod=prod*pFac/binomCoeffNenner;

          innerSum+=prod;
         
         parameter[i]=parameter[i]+parameterCoeff_[i][j]*innerSum;
          j++;

         innerSum=Complex(0,0);
         p++;
    }
      prod=Complex(1,0);
      //      p++;

    }
    std::cout<<parameter<<std::endl;

  }
  

  //  void MultiHarmonicDriver::providePolynomialCoefficients(Complex & parameter, UInt indexOfParameter)

  
 


}// end of namespace coupled field
