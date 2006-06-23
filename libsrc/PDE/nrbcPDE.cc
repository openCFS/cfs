#include <fstream>
#include <iostream>
#include <string>
#include <math.h>
#include <stdlib.h>

#include "nrbcPDE.hh" 
#include "DataInOut/Unverg/outUnverg.hh"
#include "DataInOut/GMV/outGMV.hh"
#include "Forms/abcInt.hh"
#include "Forms/nrbcInt.hh"
#include "Forms/massInt.hh"
#include "Estimator/spaceerror.hh"
#include "newmark.hh"
#include "newmarkFracDamp.hh"
#include "DataInOut/WriteInfo.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "Utils/mathfunctions.hh"
#include "Utils/nodestoresol.hh"
#include "Driver/solveStepAcoustic.hh"
#include "CoupledPDE/pdecoupling.hh"
#include "Driver/assemble.hh"

#ifdef PARALLEL
#include <mpi.h>
#endif

namespace CoupledField {

  // =========================================================================
  // set solution information
  // =========================================================================
  nrbcPDE::nrbcPDE(Grid * aptgrid, TimeFunc *aptTimeFunc, 
                   WriteResults *aptOut, std::string pdeNameWithIndex,
                   StdVector<SolutionType> localsolType)
    :SinglePDE(aptgrid,aptOut,aptTimeFunc) {

    ENTER_FCN( "nrbcPDE::nrbcPDE" );

    Error( "Not working at the moment (-> conctact Andreas)",
           __FILE__, __LINE__ );

    Double C1_, C2_, C3_, C4_, C5_;


    params->Get("absorbOrder", order_, "absorParams");
    params->Get("C1", C1_, "absorParams");
    params->Get("C2", C2_, "absorParams");
    params->Get("C3", C3_, "absorParams");
    params->Get("C4", C4_, "absorParams");
    params->Get("C5", C5_, "absorParams");

    //order_ = 5;
    //    solDofs_ = 5;

    //solDofs_=order_;
    //std::cout<<" order_ = "<< order_<<std::endl;
    
    //pdename_          = "nrbc";
    pdename_  = pdeNameWithIndex;
//     //adding a coefficient to PDE name for solving a system of similar PDEs
//     char * auxChar=new char[2];
//     auxpdeName_=new char[20];
//     sprintf(auxChar,"%i",indexofPDE_);
//     std::strcpy(auxpdeName_,pdename_.c_str());
//     std::strcat(auxpdeName_,auxChar);
    
    pdematerialclass_ = FLUID;

    isAcouCoupled_ = false;
    saveRHSval_ = false;
    saveRHSvalHist_ = false;

    // Same as for acoustic PDE
    // PDE formulation either in acoustic potential or pressure
    std::string str;
//     params->Get( "formulation", str, "pdeList", pdename_ );
//     String2Enum( str, formulation_ );
//     str = "Using * " + str + " as state variable in formulation of PDE\n";
//     Info->PrintF( pdename_, str.c_str() );
// We fix the formulation for this PDE
     String2Enum( "nrbcPhi", formulation_ );

    // class NodeStoreSol will be initialized with auxiliar function phi.
     //solTypes_ = NRBC_PHI;

    // timestepping formulation
    //params->Get( "timeSteppingFormulation", str, "pdeList", pdename_ );
    
    //effStiffMatrix is assumed for nrbcPDE!!!
    
    str="effStiffMatrix";
    
    if ( str == "effMassMatrix" ) {
      effectiveMass_ = true;
      Info->PrintF( pdename_, 
                    "      * effective mass matrix timestepping\n");
    } 
    else {
      effectiveMass_ = false;
      Info->PrintF( pdename_, 
                    "      * effective stiffness matrix timestepping\n");
    }

    // axisymmetric setup    
    isaxi_ = params->HasValue( "type", "axi", "geometry" );


    // ===========================
    // DODO GridAdaption
    // check if <specialNodes> is given, to determine whether we to write them
    // to file,
    StdVector<std::string> vecKey, vecVals, vecAttr;
    StdVector<std::string> vecNodes, vecFiles;
    vecKey = pdename_, "storeResults", "specialNodes", "saveNodes";
    vecAttr = "", "", "";
    vecVals = "", "", "";
    params->GetList(vecKey, vecAttr, vecVals, vecNodes);
    
    // specialNodes is given => set flag
    if(vecNodes.GetSize() > 0) {
      m_bWriteSpecialBCs = true;
      
      // get file to store specialNodes
      vecKey = pdename_, "storeResults", "specialNodes", "file";
      std::string strFile;
      params->Get(vecKey, strFile);
      
      // setup grid adaption object
      // it needs: the .data-file, the .xml-file and the grid
      m_pGridAdaption = new GridAdaption(strFile, aptgrid);
    }
    else {
      m_bWriteSpecialBCs = false;
    }
  }


  void nrbcPDE::DefineIntegrators() {

    ENTER_FCN( "nrbcPDE::DefineIntegrators" );

    UInt nrbcMatType=0; // 1=Rmat, 2=pmat, 3=Qmat
    Double beta,alpha;
    Double C0, Cj;
    Double coeff_Cdamp, coeff_Pmass, coeff_Rstiff, coeff_Qstiff;
    UInt actSD = 0;  // For NRBC there is only one subdomain 
    Cj =1;
    C0 =1;
    beta = ( 1/C0 + 1/Cj );
    alpha = (1/(Cj*Cj) - 1/(C0*C0));
    
    // create new entity list
    shared_ptr<ElemList> actSDList( new ElemList(ptgrid_ ) );
    actSDList->SetRegion( subdoms_[actSD] );

    //================================================== 
    //   NRBC equation
    //==================================================


    //==========================================
    // damping (mass) integrator for C matrix
    //==========================================
    coeff_Cdamp = beta;
#ifdef DEBUG
    (*debug) << std::endl << "coeff_Cdamp = "
             << coeff_Cdamp << std::endl << std::endl;
#endif
    // we use the mass integrator since it is the same as the damping one
    BaseForm * bilinear_Cdamp  = new MassInt(coeff_Cdamp, order_, isaxi_);
    BiLinFormContext * CdampContext = 
      new BiLinFormContext(bilinear_Cdamp, DAMPING );


    CdampContext->SetPtPdes(this, this);
    CdampContext->SetResults( results_[0], results_[0],
                              actSDList, actSDList );
    //      assemble_->AddBiLinearForm(CmassIntDescr );

    
    //============================================
    // nrbc stiffness integrator for R matrix
    //============================================
    coeff_Rstiff = 1;
    nrbcMatType = 1;
#ifdef DEBUG
    (*debug) << std::endl << "coeff_Rstiff = "
             << coeff_Rstiff << std::endl << std::endl;
#endif
    BaseForm * bilinear_Rstiff  = new NrbcInt(coeff_Rstiff, order_,
                                              nrbcMatType, isaxi_);
    BiLinFormContext * RstiffContext = 
      new BiLinFormContext(bilinear_Rstiff, STIFFNESS );
    RstiffContext->SetPtPdes(this, this);
    RstiffContext->SetResults( results_[0], results_[0],
                               actSDList, actSDList );


    // =========================================
    //         ALSO P and Q go to the LHS
    // =========================================
    
    //=========================================
    // nrbc mass integrator for P matrix
    //=========================================
    coeff_Pmass = alpha;
    nrbcMatType = 2;
#ifdef DEBUG
    (*debug) << std::endl << "coeff_Pmass = "
             << coeff_Pmass << std::endl << std::endl;
#endif
    BaseForm * bilinear_Pmass  = 
      new NrbcInt(coeff_Pmass, order_,
                  nrbcMatType, isaxi_);
    BiLinFormContext * PmassContext = 
      new BiLinFormContext(bilinear_Pmass, MASS );
    PmassContext->SetPtPdes(this, this);
    PmassContext->SetResults( results_[0], results_[0],
                              actSDList, actSDList );
    //==================================================
    // nrbc tangential derivatives integrator for Q matrix
    //==================================================
    coeff_Qstiff = 1.0;
    // nrbcMatType = 3 corresponds to tangential integrator
    nrbcMatType = 3;
#ifdef DEBUG
    (*debug) << std::endl << "coeff_Qdamp = "
             << coeff_Qstiff << std::endl << std::endl;
#endif
    BaseForm * bilinear_Qstiff  = 
      new NrbcInt(coeff_Qstiff, order_,
                  nrbcMatType, isaxi_);
    BiLinFormContext * QstiffContext = 
      new BiLinFormContext(bilinear_Qstiff, STIFFNESS );
    QstiffContext->SetPtPdes(this, this);
    QstiffContext->SetResults( results_[0], results_[0],
                               actSDList, actSDList );

    //Add the integrators
    assemble_->AddBiLinearForm( CdampContext );
    assemble_->AddBiLinearForm( RstiffContext );
    assemble_->AddBiLinearForm( PmassContext );
    assemble_->AddBiLinearForm( QstiffContext );
    
    
    // Give result to equation numbering class
    eqnMap_->AddResult( *results_[0], actSDList );  
    
  }
  
    void nrbcPDE::DefineSolveStep() {
      ENTER_FCN( "nrbcPDE::DefineSolveStep" );
      
      solveStep_ = new SolveStepAcoustic(*this);
    }
    


  // ======================================================
  // TIME STEPPING SECTION
  // ======================================================

  void nrbcPDE::InitTimeStepping() {
    ENTER_FCN( "nrbcPDE::InitTimeStepping" );

    // this includes rayleigh and thermoviscous damping
    if ( effectiveMass_ == false ) {
      TS_alg_ = new Newmark( algsys_ );
    }
    else {
      TS_alg_ = new NewmarkEffMass( algsys_ );
    }
    
    
  }



  // ======================================================
  // COUPLING SECTION
  // ======================================================

  void nrbcPDE::InitCoupling(PDECoupling * Coupling) {
    
    ENTER_FCN( "nrbcPDE::InitCoupling" );
    
    isIterCoupled_ = true;
    ptCoupling_   = Coupling;
    
    // Intialize the memory of the coupling values
    for (UInt i=0; i<ptCoupling_->GetNumOutputCouplings(); i++) {
      if (ptCoupling_->GetOutputQuantity(i) == NRBC_PHI)    {
        ptCoupling_->CreateCouplingVector(i,isComplex_);
        
        // now since we need a incremental formulation, 
		//  initialize some necessary vectors
        isIncrFormulation_ = true;
      }
    }
  }
  

  void nrbcPDE::CalcOutputCoupling() {

    ENTER_FCN( "nrbcPDE::CalcOutputCoupling" );

    UInt dof;
    SolutionType quantity;
    StdVector<Elem*> * couplingElems = NULL;
    StdVector<UInt> * couplingNodes = NULL;
    CFSVector * temp_values = NULL;

    StdVector<SolutionType> SolTypWithPDEindex;
  
    // at first, check if this PDE is iterative coupled
    if (isIterCoupled_ == false)
      return;

    // loop over all output coupling quantities
    for (UInt i=0; i<ptCoupling_->GetNumOutputCouplings(); i++) {
      quantity = ptCoupling_->GetOutputQuantity(i);
      ptCoupling_->GetOutputValues(i, temp_values);

      // hard coded cast, since coupling is only possible with
      // real valued entries
      Vector<Double> * values = dynamic_cast<Vector<Double>*>(temp_values);

      switch(ptCoupling_->GetOutputType(i)) {

      case NODE:
        if (quantity == NRBC_PHI) {
          ptCoupling_->GetOutputElements(i, couplingElems);
          ptCoupling_->GetOutputNodes(i, couplingNodes);
          dof = ptCoupling_->GetOutputDof(i);
          
          CalcAcouCouplingRHS(couplingElems, *couplingNodes,
                              *values, dof);                              
        }
        break;

      case ELEM:
        Error("No Element coupling output", __FILE__,__LINE__);
      }
    }
  }


  void nrbcPDE::
  CalcAcouCouplingRHS( StdVector<Elem*> * couplingElems, 
                       StdVector<UInt> & couplingNodes,
                       Vector<Double>& elemCouplingSols,
                       UInt couplingdof ) {
    
    ENTER_FCN( "nrbcPDE::CalcAcouCouplingRHS" );
    
    Error( "Not working yet", __FILE__, __LINE__ );
  //   // This function should be called only to send Phi1 to acousticPDE
//     Double density = 0.0;
//     Double sign = 0.0;
//     Integer matIndex = -1;
//     Elem * ptVolElem = NULL;
//     Matrix<Double> ptCoord, elemMat;
//     Vector<Double> sol, PHI1sol, NRBCPhi1surfelem;
    
//     elemCouplingSols.Init(0.0);
    
//     for (UInt actElem=0; actElem<couplingElems->GetSize(); actElem++) {
      
//       // Perform cast from volume element to surface element, since
//       // mech-acou coupling makes only sense on surface elements
//       SurfElem * actCoupleElem = 
//         dynamic_cast<SurfElem*> ((*couplingElems)[actElem]);
      
//       BaseFE * ptElem = actCoupleElem->ptElem;
//       StdVector<UInt> & connecth = actCoupleElem->connect;
//       GetElemCoords(connecth, ptCoord);
      
//       // We have only one vol element neighbor and the whole domain is
//       // acoustic with the same mat parameters so no need to search more
//       matIndex = subdoms_.Find(actCoupleElem->ptVolElem1->regionId);
//       ptVolElem = actCoupleElem->ptVolElem1;
//       sign = -1.0 * actCoupleElem->normalSign;

      
//       if ( matIndex == -1 && (actCoupleElem->ptVolElem2!=NULL)) {
//         (*error) << "nrbcPDE::CalcAcouCouplingRHS: The volume "
//                  << "element neighbour of surface element Nr. "
//                  << actCoupleElem->elemNum << " do not belong to my regions!";
//         Error( __FILE__, __LINE__ );
//       }
      
//       // Assign correct density
//  //      density = materialData_[matIndex].GetDensity();
//       density=1.0;
      
//       BaseForm * bilinear_mass = new MassInt(ptElem, density, isaxi_);
//       bilinear_mass->CalcElementMatrix(ptCoord, elemMat);
//       delete bilinear_mass;     
      
//       GetSolVecOfElement(sol, connecth);
// //       for(UInt k=0; k<sol.GetSize(); k++)
// //         if (abs(sol[k])<=5e-17)
// //           sol[k]=0;
//       // solution for J=1 which is the one to be sent to acouPDE
//       PHI1sol.Resize(connecth.GetSize());

//       //std::cout<<"order_= "<<order_<<std::endl;

//       for (UInt actNode=0; actNode < connecth.GetSize(); actNode++)
//         for (UInt actDof=0; actDof<1; actDof++)//only for PHI_x,1 (J=1)
//           PHI1sol[actNode] = sol[actDof + actNode*order_];
      
//       //std::cout<<"---NRBC-PDE COUPLING OUTPUT----------------"<<std::endl;
//       //   std::cout<<"PHI1sol = "<<std::endl;
//       //   std::cout<<PHI1sol<<std::endl;
     
//       NRBCPhi1surfelem = elemMat * PHI1sol; //surface elem vector

// //         std::cout<<"elemMat = "<<std::endl;
// //          std::cout<<elemMat<<std::endl;

// //          std::cout<<"NRBCPhi1surfelem = "<<std::endl;
// //          std::cout<<NRBCPhi1surfelem<<std::endl;

//       for (UInt actNode=0; actNode<ptCoord.GetSizeRow(); actNode++) {
//         UInt nodePos = 0;
        
//         while(connecth[actNode] != couplingNodes[nodePos] &&
//               nodePos < couplingNodes.GetSize()) {
//           nodePos++;
//         }
//         elemCouplingSols[nodePos] += NRBCPhi1surfelem[actNode];

// //          std::cout<<"elemCouplingNRBCPhi1 = "<<std::endl;
// //          std::cout<<elemCouplingSols<<std::endl;

//       }   
//     }

//     //std::cout<<"---------END NRBC-PDE COUPLING OUTPUT INFO---------"<<std::endl;     

  }





  bool nrbcPDE::HasOutput(SolutionType output) {
    ENTER_FCN( "nrbcPDE::HasOutput" );
    if (output == NRBC_PHI) {
      return true;
    }
    return false;
  }


  // ======================================================
  // POSTPROCESSING SECTION
  // ======================================================

  void nrbcPDE::WriteResultsInFile(const UInt kstep,
                                       const Double asteptime,
                                       UInt stepOffset,
                                       Double timeOffset) {
    ENTER_FCN( "nrbcPDE::WriteResultsInFile" );

#ifdef PARALLEL //only one thread should write the output
    int commrank;
    MPI_Comm_rank(MPI_COMM_WORLD,&commrank);
    if (!commrank) {
#endif
      NodeStoreSol<Double> solIm_mesh;
      NodeStoreSol<Double> * solTransient;
      NodeStoreSol<Complex> * solHarmonic;
      
      
      Double actTime = asteptime + timeOffset;
      UInt actStep = kstep + stepOffset;
      
      if (analysistype_==HARMONIC) {
        solHarmonic = dynamic_cast<NodeStoreSol<Complex>*>(sol_);

        if (saveSol_)
          outFile_->WriteNodeSolutionHarmonic(*solHarmonic,  actStep, 
                                              actTime, complexFormat_);
        if (saveSolHist_)
          outFile_->WriteNodeHistoryHarmonic(*solHarmonic,  actStep, 
                                             actTime, complexFormat_);
      }
      else {  
        solTransient = dynamic_cast<NodeStoreSol<Double>*>(sol_);


        if (saveSol_){
          outFile_->WriteNodeSolutionTransient(*solTransient,actStep,actTime);
         
        }
        if (saveSolHist_){
          outFile_->WriteNodeHistoryTransient(*solTransient, actStep,actTime);
        }

	// DODO here
	if(m_bWriteSpecialBCs) {
	  m_pGridAdaption->Add2DataFile(*solTransient, actTime);
	}
        
        if (saveDeriv1_) {
          solDeriv1_.SetAlgSysVector(getS1()); 
          outFile_->WriteNodeSolutionTransient(solDeriv1_, actStep, actTime);

        }
        
        if (saveDeriv1Hist_) {
          solDeriv1_.SetAlgSysVector(getS1()); 
          outFile_->WriteNodeHistoryTransient(solDeriv1_, actStep, actTime);
        }

        if (saveDeriv2_) {
          solDeriv2_.SetAlgSysVector(getS2());
          outFile_->WriteNodeSolutionTransient(solDeriv2_, actStep, actTime);
        }
        if (saveDeriv2Hist_){
          solDeriv2_.SetAlgSysVector(getS2());
          outFile_->WriteNodeHistoryTransient(solDeriv2_, actStep, actTime);
        }
        if (saveRHSval_){
          outFile_->WriteNodeSolutionTransient(rhs_, actStep, actTime);
        }
        if (saveRHSvalHist_){
          outFile_->WriteNodeHistoryTransient(rhs_, actStep, actTime); 
        }
      }


#ifdef PARALLEL
    }//!commrank
#endif
  }
  
  void nrbcPDE::WriteHistoryInFile(const UInt kstep,
                                       const Double asteptime,
                                       UInt stepOffset,
                                       Double timeOffset) {
    ENTER_FCN( "nrbcPDEPDE::WriteHistoryInFile" );

#ifdef PARALLEL //only one thread should write the output
    int commrank;
    MPI_Comm_rank(MPI_COMM_WORLD,&commrank);
    if (!commrank) {
#endif
      NodeStoreSol<Double> solIm_mesh;
      NodeStoreSol<Double> * solTransient;
      NodeStoreSol<Complex> * solHarmonic;
      
      
      Double actTime = asteptime + timeOffset;
      UInt actStep = kstep + stepOffset;
      
      if (analysistype_==HARMONIC) {
        solHarmonic = dynamic_cast<NodeStoreSol<Complex>*>(sol_);

        if (saveSolHist_)
          outFile_->WriteNodeHistoryHarmonic(*solHarmonic,  actStep, 
                                             actTime, complexFormat_);

        if (saveDeriv1Hist_) {
          // multiply solution with j * omega
        }
      }
      else {  
        solTransient = dynamic_cast<NodeStoreSol<Double>*>(sol_);

        if (saveSolHist_){
          outFile_->WriteNodeHistoryTransient(*solTransient, actStep,actTime);
        }

        if (saveDeriv1Hist_) {
          solDeriv1_.SetAlgSysVector(getS1()); 
          outFile_->WriteNodeHistoryTransient(solDeriv1_, actStep, actTime);
        }

        if (saveDeriv2Hist_){
          solDeriv2_.SetAlgSysVector(getS2());
          outFile_->WriteNodeHistoryTransient(solDeriv2_, actStep, actTime);
        }

        if (saveRHSvalHist_){
          outFile_->WriteNodeHistoryTransient(rhs_, actStep, actTime); 
        }
        
      }
 
#ifdef PARALLEL
    }//!commrank
#endif
  }
  
    
  // ***********************************************************************
  //   Obtain information on desired output quantities from parameter file
  // ***********************************************************************
  void nrbcPDE::ReadStoreResults() {  
    ENTER_FCN( "nrbcPDE::ReadStoreResults" );
    
    // Construct vectors for restricted parameter search
    StdVector<std::string> keyVec;
    StdVector<std::string> attrVec;
    StdVector<std::string> valVec;
    std::string quantity;
    
    // *****************************
    // Determine nodal results
    // ***************************** 
    StdVector<std::string> nodeValues;
    keyVec  = pdename_, "storeResults", "nodeResults", "region";
    attrVec = "", "", "type";  

    if ( formulation_ == NRBC_PHI ) {
      // --- pressure ---
      Enum2String(ACOU_PRESSURE, quantity);
      valVec = "", "", quantity;
      params->GetList( keyVec, attrVec, valVec, nodeValues);
      if (nodeValues.GetSize() > 0) {
        (*error) << "It makes no sense to have a PDE in element values "
                 << "and retrieve \n node results in pressure. "
                 << "Try element results!";
        Error( __FILE__, __LINE__ );
      }

      // --- NRBC_PHI ---
      Enum2String(NRBC_PHI, quantity);
      valVec = "", "", quantity;
      params->GetList( keyVec, attrVec, valVec, nodeValues);
      if (nodeValues.GetSize() > 0) {
        saveSol_ = true;
        hasOutput_ = true;
      }

      // --- acoustic rhsval ---
      Enum2String(ACOU_RHSVAL, quantity);
      valVec = "", "", quantity;
      params->GetList( keyVec, attrVec, valVec, nodeValues);
      if (nodeValues.GetSize() > 0) {
	saveRHSval_ = true;
	hasOutput_ = true;
      }



    }
  else if ( formulation_ == ACOU_PRESSURE ) 
    {
      // --- pressure ---
      Enum2String(ACOU_PRESSURE, quantity);
      valVec = "", "", quantity;
      params->GetList( keyVec, attrVec, valVec, nodeValues);
      if (nodeValues.GetSize() > 0) {
        std::string warnMsg;
        warnMsg = "Due to the restrictions in the .unv file format, the ";
        warnMsg += "acoustic pressure is written as acoustic (fluid) ";
        warnMsg += "potential!";
        Warning(warnMsg.c_str(), __FILE__, __LINE__);

        sol_->SetSolutionType(ACOU_PRESSURE);
        sol_->SetNumDofs(1);
        sol_->Init();
        saveSol_ = true;
        hasOutput_ = true;
      }
    }
    
    
    // *****************************
    // Determine element results
    // *****************************
    StdVector<std::string> elementValues;
    keyVec  = pdename_, "storeResults", "elementResults", "region";
    
    // --- much to do here ---

    // *****************************
    // Determine nodal history
    // *****************************
    StdVector<std::string> saveNodeHist; 
    keyVec  = pdename_, "storeResults", "nodeHistory", "saveNodes";
    attrVec = "", "", "type";

    if ( formulation_ == NRBC_PHI) {
      // pressure
      Enum2String(ACOU_PRESSURE, quantity);
      valVec = "", "", quantity;
      params->GetList( keyVec, attrVec, valVec, saveNodeHist);
      if (saveNodeHist.GetSize() > 0) {
        (*error) <<  "It makes no sense to have a PDE in acoustic potential "
                 << " and retrieve \n history node results in pressure.";
        Error( __FILE__, __LINE__ );
      }

      // --- acoustic potential  ---
      Enum2String(NRBC_PHI, quantity);
      valVec  = "", "", quantity;
      params->GetList( keyVec, attrVec, valVec, saveNodeHist );
    
      if (saveNodeHist.GetSize() > 0) {
        saveSolHist_ = true;
        hasOutput_ = true;
        Info->PrintF( pdename_, "Saving acouPotential for Nodes:\n" );
        for ( UInt k = 0; k < saveNodeHist.GetSize(); k++ ) {
          Info->PrintF( pdename_, "  %s\n", saveNodeHist[k].c_str() );
        }

      }


      // --- acoustic RHS ---
      Enum2String(ACOU_RHSVAL, quantity);
      valVec  = "", "", quantity;
      params->GetList( keyVec, attrVec, valVec, saveNodeHist );
    
      if (saveNodeHist.GetSize() > 0) {
        saveRHSvalHist_ = true;
        hasOutput_ = true;
        Info->PrintF( pdename_, "Saving acouRHSval for Nodes:\n" );
        for ( UInt k = 0; k < saveNodeHist.GetSize(); k++ ) {
          Info->PrintF( pdename_, "  %s\n", saveNodeHist[k].c_str() );
        }


      }


    }

    else if ( formulation_ == ACOU_PRESSURE ) {
      // --- NRBC_PHI ---
      Enum2String(NRBC_PHI, quantity);
      valVec = "", "", quantity;
      params->GetList( keyVec, attrVec, valVec, nodeValues);
      if (nodeValues.GetSize() > 0) {
        Error("It makes no sense to have a PDE in pressure and retrieve \n \
results in acoustic potential.", __FILE__,__LINE__);
      }

      // --- pressure ---
      Enum2String(ACOU_PRESSURE, quantity);
      valVec  = "", "", quantity;
      params->GetList( keyVec, attrVec, valVec, saveNodeHist );
    
      if (saveNodeHist.GetSize() > 0) {
        sol_->SetSolutionType(ACOU_PRESSURE);
        sol_->SetNumDofs(1);
        sol_->Init();
        saveSolHist_ = true;
        hasOutput_ = true;
        Info->PrintF( pdename_, "Saving acouPotential for Nodes:\n" );
        for ( UInt k = 0; k < saveNodeHist.GetSize(); k++ ) {
          Info->PrintF( pdename_, "  %s\n", saveNodeHist[k].c_str() );
        }
      }
    }
  }


} // end of namespace
