#include <fstream>
#include <iostream>
#include <string>

#include "multiHarmonicDriver.hh"
#include "stdSolveStep.hh"


#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "Domain/domain.hh"
#include "PDE/StdPDE.hh"
#include "PDE/piezoPDE.hh"
#include "Materials/baseMaterial.hh"
#include "Utils/elemstoresol.hh"


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

    ptMHFiles_ = new MHMaterialDataFile;

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

    adjustDamping_ = FALSE;
    std::string damp;
    keyVec = "multiHarmonic", "adjustDamping";
    //   adjustDamping_ =  params->IsSet(keyVec, attrVec, valVec, adjustDamping_);
    adjustDamping_ = params->IsSet("adjustDamping",  "multiHarmonic");


 
    std::cout<<"\n Opens impedCurve.dat and piezoLog.dat ... "<<std::endl;

    std::string filename= "imped.dat";

    impedCurve = new std::ofstream(filename.c_str(),std::basic_ios<char>::out);

    if (!impedCurve)
      {
        std::cout << "\n ImpedanceCurve.dat could not be initialized" << std::endl;
      }

    std::string filenameLog= "piezoLog.dat";

    piezoLog = new std::ofstream(filenameLog.c_str(),std::basic_ios<char>::out);
    
    if (!piezoLog)
      {
        std::cout << "\n piezoLog.dat could not be initialized" << std::endl;
      }

    std::string filenameParLog= "parLog.dat";
    parLog = new std::ofstream(filenameParLog.c_str(),std::basic_ios<char>::out);

    if (!parLog)
      {
        std::cout << "\n piezoLog.dat could not be initialized" << std::endl;
      }

  } // end of constructor

  // destructor
  
  MultiHarmonicDriver::~MultiHarmonicDriver(){
    ENTER_FCN( "MultiHarmonicDriver::~MultiHarmonicDriver" );
    //    allMeasuredData->close();
    impedCurve->close();
    piezoLog->close();
    parLog->close();
  }

  void MultiHarmonicDriver::SolveProblem() {
    ENTER_FCN( "MultiHarmonicDriver::SolveProblem" );
    
    PiezoPDE * ptPiezoPDE;
    if (! isPartOfSequence_)
      ptdomain_->PrintGrid();
    if (isPartOfSequence_ == FALSE){     
      GetMyPDEs();
    
      //! cast pointer to BasePDE * to pointer of SinglePDE *
      ptMyPDE_ = dynamic_cast<SinglePDE*>(ptPDE_);

      
      ptPiezoPDE = dynamic_cast<PiezoPDE*>(ptPDE_);
      Info->StartProgress ("Starting to solve problem", FALSE);
    }
    ptMyPDE_->WriteGeneralPDEdefines();
    //    pdeId_ = ptMyPDE_->GetPDEId();

    std::map<RegionIdType, BaseMaterial*> ptMaterial;

    ptMaterial = ptMyPDE_->getPDEMaterialData();   // Pointer to MaterialData

    Assemble * ptAssemble_;

    ptAssemble_ = ptMyPDE_->getPDE_assemble();

    ptAssemble_->SetMaterialPointer(ptMaterial);


    // ptMHFiles_->computeIndexSet(2,2,nrMultHarms_);
    //     ptMHFiles_->computeIndexSet(3,3,nrMultHarms_);
    //     std::cout<<"    ptMHFiles_->computeIndexSet2,2,nrMultHarms_ went fine " <<std::endl;


    //  MHAssembleMatrices();

    UInt element;
    element=2;
    //    Integer delta=3;

   //  calcParameterCurveAtElement(parameter_, parameterCoeff_, element, nrMultHarms_, delta, maxP);
//     calcParameterCurveAtElement(parameter_, parameterCoeff_, element, nrMultHarms_, 1, maxP);
    
    Matrix<Integer> exps;
    Matrix<Integer> count;
   
    Boolean reset = TRUE;
        
    UInt fstep;
    Double actFreq  = startFreq_;
    Double freqIncr;
    
    if (numFreq_ > 1) {
      freqIncr = (stopFreq_ - startFreq_) / (numFreq_-1);
    }
    else {
      freqIncr = stopFreq_ - startFreq_ ;
    }
    
    std::string errMsg;
  
  
    // branch for single PDE
    for (fstep = 1; fstep <= numFreq_; fstep++) {
      Info->WriteHarmonicStep(ptPDE_->GetName(), fstep, actFreq);
    
      ptPDE_->GetSolveStep()->SetActFreq(actFreq);
      ptPDE_->GetSolveStep()->SetActStep(fstep);
      std::cout<<" Step in multiharmonicDriver =  "<< fstep <<std::endl;
      //      std::cout<<"\n multiHarm: 1 " <<std::endl;
      ptPDE_->GetSolveStep()->PreStepHarmonic(reset);   
      //      std::cout<<"\n multiHarm: 2"  <<std::endl;
      ptPDE_->GetSolveStep()->SolveStepHarmonic(reset);
      //      std::cout<<"\n multiHarm: 3 " <<std::endl;
      ptPDE_->GetSolveStep()->PostStepHarmonic(reset);

      ptMyPDE_->PostProcess();

      Grid * ptGrid;
      NodeEQN * ptNodeEqn;
      //      Assemble * ptAssemble;
      //      Domain * ptDomain;
     

      Boolean reset = TRUE;
      ptGrid = ptMyPDE_->getPDE_grid();
      ptNodeEqn = ptMyPDE_->getPDE_eqnData();


      StdVector<Elem*> elemssd;
      StdVector<RegionIdType>  subdoms = ptMyPDE_->getPDE_subdoms();
      //      ptMaterial=ptMyPDE_->getPDEMaterialData();   // Pointer to MaterialData
      ptGrid->GetVolElems(elemssd,subdoms[0]); // gets element list elemssd
      
      BaseNodeStoreSol * ptSol = ptMyPDE_->getPDESolution();
      NodeStoreSol<Complex> * ptNodeStoreSol;
      ptNodeStoreSol = dynamic_cast<NodeStoreSol<Complex>*>(ptSol);

      StdVector<Integer> connect_PDE;                 
      //      getchar();

      EfieldInZDir_.Resize(elemssd.GetSize());
    //loop over elements
      for (UInt actEl=0; actEl< elemssd.GetSize(); actEl++) {
        BaseFE * ptEl = elemssd[actEl]->ptElem;
        StdVector<UInt> connecth = elemssd[actEl]->connect;
        
        Matrix<Double> ptCoord;
        ptGrid->GetElemNodesCoord(ptCoord, connecth);

//         std::cout<<" ElementNummer = " << actEl <<std::endl;
//         std::cout<<(ptCoord)<<std::endl;        
        // map connect to PDE node numbers
        
        const UInt nrIntPts = ptEl->GetNumIntPoints(); 

        ptNodeEqn->Node2EQN(connecth, connect_PDE);
//         std::cout<<connect_PDE<<std::endl;
//         std::cout<<connecth<<std::endl;
        Vector<Complex> elSolVec; 
        ptNodeStoreSol->GetElemSolution(elSolVec,connecth);
        Vector<Complex> solVecElecPot;
        UInt dim = ptMyPDE_->getPDE_spaceDim();
        
        if (dim==3){
          solVecElecPot.Resize(elSolVec.GetSize()/4);
          UInt j=0;
          for (UInt i=3; i<elSolVec.GetSize();i=i+4){
            solVecElecPot[j] = elSolVec[i];
            j++;
          }
        }
        else if (dim==2){
          solVecElecPot.Resize(elSolVec.GetSize()/3);
          UInt j=0;
          for (UInt i=2; i<elSolVec.GetSize();i=i+3){
            solVecElecPot[j] = elSolVec[i];
            j++;
          }
        }
//         std::cout<<"solVecElecPot:"<<std::endl;
//         std::cout<<solVecElecPot<<std::endl;

        //Vector<Complex> globSolVec;

        //        ptNodeStoreSol->GetSolVectorSingleDof(ELEC_POTENTIAL, 0, globSolVec);
        //        std::cout<<" Globaler Solution Vector ElecPot : " <<std::endl;
        // std::cout<<globSolVec<<std::endl;

        Complex ElecFieldInDirZ;
        ElecFieldInDirZ = Complex(0.0,0.0);

        for ( UInt actIntPt = 1; actIntPt <= nrIntPts; actIntPt++ ) {

          Matrix<Double> xiDx;
          ptEl->GetGlobDerivShFncAtIp(xiDx, actIntPt, ptCoord);
//           std::cout<<"xiDx"<<std::endl;
//           std::cout<<xiDx[actIntPt-1][1]<<std::endl;
          ElecFieldInDirZ+=xiDx[actIntPt-1][1]*solVecElecPot[actIntPt-1];

          //          Matrix<Complex> 
        }
//         std::cout <<" ElecFieldInDirZ : " <<std::endl;
//         std::cout << ElecFieldInDirZ  <<std::endl;
       
        EfieldInZDir_[actEl]=ElecFieldInDirZ;
        
        ptNodeEqn->Node2EQN(connecth, connect_PDE);
        //        std::cout<<connect_PDE<<std::endl;
        // std::cout<<connecth<<std::endl;
        

        //        ptElemFE->GetShFncAtIp(shFnc, actIntPt);

        //      std::cout<<elSolVec<<std::endl;
        //        getchar();
        
        std::map<RegionIdType, BaseMaterial*> actSDMat = ptMaterial;
        //        Boolean isdamping=TRUE;
      }
      
      ptPiezoPDE->SetEfieldInZDir_(EfieldInZDir_);

   
    //    std::cout<<"piezoParamIdent::createAndSetRHSforJacobian 1 "<< std::endl; 
   
        ptPDE_->WriteResultsInFile( fstep, actFreq);

        //write history data
        ptPDE_->WriteHistoryInFile(fstep, actFreq);
        
        actFreq += freqIncr;
       
    }
  } // end Solve

}// end of namespace coupled field
