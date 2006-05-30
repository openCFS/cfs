#include <fstream>
#include <iostream>
#include <sstream>
#include <math.h>

#include "DataInOut/Unverg/outUnverg.hh"
#include "DataInOut/GMV/outGMV.hh"
#include "Forms/forms_header.hh"
#include "Forms/linElecInt.hh"
#include "Estimator/spaceerror.hh"
#include "DataInOut/WriteInfo.hh"
#include "Driver/assemble.hh"
#include "General/defs.hh"

#include <Matrix/matrix.hh>
#include <Utils/vector.hh>
#include <Forms/gradfieldop.hh>
#include "elecPDE.hh"
#include <General/defs.hh>
#include "DataInOut/ParamHandling/BaseParamHandler.hh" 
#include <string>
#include "Utils/StdVector.hh"
#include "Driver/solveStepElec.hh"
#include "CoupledPDE/pdecoupling.hh"

#ifdef TCL_INTERFACE
#include "DataInOut/Scripting/cfsmessenger.hh" 
#endif

#ifdef PARALLEL
#include <mpi.h>
#endif


namespace CoupledField {


  // ***************
  //   Constructor
  // ***************
  ElecPDE::ElecPDE(Grid * aptgrid,  TimeFunc *aptTimeFunc, WriteResults *aptOut)
    :SinglePDE(aptgrid, aptOut, aptTimeFunc) {

    ENTER_FCN( "ElecPDE::ElecPDE" );

    // =====================================================================
    // set solution information
    // =====================================================================
    dofspernode_ = 1;
    solTypes_ = ELEC_POTENTIAL;
    solDofs_ = 1;
    pdename_          = "electrostatic";
    pdematerialclass_ = ELECTROSTATIC;
 
    geoUpdate_ = FALSE;
    nonLin_    = FALSE;
    isAlwaysStatic_ = TRUE;
    isPiezoCoupled_ = FALSE;

    //check, if problem is axisymmetric
    if ( params->HasValue( "type", "axi", "geometry" ) ) isaxi_ = TRUE;

  }
  


  void ElecPDE::DefineIntegrators()
  {
    ENTER_FCN( "ElecPDE::DefineIntegerators" );

    BaseForm *form, *formC;

   //transform the type
    SubTensorType tensorType;

    if ( dim_ == 3 ) {
      tensorType = FULL;
    }
    else {
      if ( isaxi_ == TRUE ) {
	tensorType = AXI;
      }
      else {
	// 2d: plane case
	tensorType = PLANE_STRAIN;
      }
    }

    // if the pde is piezo-coupled, the electrostatic entries
    // have to multiplied with -1
    Double factor = 1.0;
    if ( isPiezoCoupled_ == TRUE )
      factor *= -1.0;  
  
    for ( UInt actSD = 0; actSD < subdoms_.GetSize(); actSD++ ) {
      form = new linElecInt( materials_[subdoms_[actSD]], tensorType );
      form->SetFactor( factor );

      IntegratorDescriptor * stiffIntDescr = 
        new IntegratorDescriptor(form, STIFFNESS, nonLin_);

      stiffIntDescr->SetPDEIds(this, this);
      assemble_->AddIntegrator(stiffIntDescr, subdoms_[actSD]);

      // check for complex valued material parameter
      if( params->HasValue( "type", "imagMaterialParameter", 
                            "materialDataType" ) ) {
        DataType matType = IMAG; 

        formC = new linElecInt( materials_[subdoms_[actSD]], tensorType );
	formC->SetFactor( factor );
	formC->SetMatDataType(matType);

	IntegratorDescriptor * stiffIntDescrC = 
	  new IntegratorDescriptor(formC, STIFFNESS, nonLin_);

	stiffIntDescrC->SetPDEIds(this, this);
	stiffIntDescrC->SetMatDataType(matType);
	assemble_->AddIntegrator(stiffIntDescrC, subdoms_[actSD]);
      }

    }
  }

  void ElecPDE::DefineSolveStep() {
    ENTER_FCN( "ElecPDE::DefineSolveStep" );
    solveStep_ = new SolveStepElec(*this);
  }

  void ElecPDE::ReadSpecialBCs( ) {
    ENTER_FCN( "ElecPDE::ReadSpecialBCs" );
   
    StdVector<UInt> nodes;
    std::string aux;

    // Get values from parameter file
    StdVector<std::string> keyVec, attrVec, valVec;
    StdVector<std::string> node1, node2;
    StdVector<Double> resistance, capacitance, inductance;
    attrVec = "", "tag", "";
    valVec = "", bcSequenceTag_, "";
    
    keyVec = pdename_, "bcsAndLoads", "impedance", "node1";
    params->GetList( keyVec, attrVec, valVec, node1 ) ;

    keyVec = pdename_, "bcsAndLoads", "impedance", "node2";
    params->GetList( keyVec, attrVec, valVec, node2 ) ;

    keyVec = pdename_, "bcsAndLoads", "impedance", "resistance";
    params->GetList( keyVec, attrVec, valVec, resistance ) ;

    keyVec = pdename_, "bcsAndLoads", "impedance", "capacitance";
    params->GetList( keyVec, attrVec, valVec, capacitance ) ;

    keyVec = pdename_, "bcsAndLoads", "impedance", "inductance";
    params->GetList( keyVec, attrVec, valVec, inductance ) ;

    // Check, that all vectors have the same length
    if ( node1.GetSize() != node2.GetSize() ||
         node1.GetSize() != resistance.GetSize() ||
         node1.GetSize() != capacitance.GetSize() ||
         node1.GetSize() != inductance.GetSize() ) {
      (*error) << "Inconsistent definition of impedances!";
      Error( __FILE__, __LINE__ );
    }

    // Make sure that we are in harmonic mode
    if ( node1.GetSize() != 0 && analysistype_ != HARMONIC) {
      (*error) << "Impedances are only available in harmonic simulation "
               << "at the moment!";
      Error( __FILE__, __LINE__ );
    }
    
    // Create impedances
    impedances_.Resize( node1.GetSize() );
    std::ostringstream out;
    for ( UInt i = 0; i < impedances_.GetSize(); i++ ) {

      ptgrid_->GetNodesByName( nodes, node1[i] );
      impedances_[i].node1 = nodes[0];

      ptgrid_->GetNodesByName( nodes, node2[i] );
      impedances_[i].node2 = nodes[0];
      
      impedances_[i].resistance = resistance[i];
      impedances_[i].capacitance = capacitance[i];
      impedances_[i].inductance = inductance[i];
      
      out.clear();
      out.str("");
      out << "Impedance defined between node " << impedances_[i]. node1
          << " and node " << impedances_[i].node2  << " with:\n";
      Info->PrintF( pdename_, out.str().c_str() );
      out.clear();
      out.str("");
      out << "\tR = " << impedances_[i].resistance 
          << "\tL = " << impedances_[i].inductance         
          << "\tC = " << impedances_[i].capacitance << std::endl;
      Info->PrintF( pdename_, out.str().c_str() );
      Info->PrintF( pdename_, "\n" );
    }
}





  // ======================================================
  // POSTPROCESSING SECTION
  // ======================================================


  // **********************
  //   WriteResultsInFile
  // **********************
  void ElecPDE::WriteResultsInFile( const UInt kstep,
                                    const Double asteptime,
                                    UInt stepOffset,
                                    Double timeOffset ) {

    ENTER_FCN( "ElecPDE::WriteResultsInFile" );

    UInt actStep = kstep + stepOffset;
    Double actTime = timeOffset + asteptime;

    // only first process should write output
#ifdef PARALLEL
    int commrank;
    MPI_Comm_rank( MPI_COMM_WORLD, &commrank );
    if ( commrank != 0 ) {
      return;
    }
#endif
    
    // ATTENTION:
    // The errorMap should be assigned as a StoreSolution, not as a 
    // Vector. This is only temporary!
    if ( !isComplex_ ) {

      // Down-cast
      NodeStoreSol<Double> *solConverted;
      solConverted = dynamic_cast<NodeStoreSol<Double>*>(sol_);

      // write electric potential
      if (saveSol_) {
        outFile_->WriteNodeSolutionTransient(*solConverted, actStep, actTime);
      }

      if ( calcEfield_.GetSize() != 0 ) {
        ElemStoreSol<Double> & eFieldConv = 
          dynamic_cast<ElemStoreSol<Double>& >(*E_);
        outFile_->WriteElemSolutionTransient(eFieldConv, actStep, actTime);
      }

      if ( calcCharges_.GetSize() != 0 ) {
        ElemStoreSol<Double> & chargeConv = 
          dynamic_cast<ElemStoreSol<Double>& >(*charges_);
        outFile_->WriteElemSolutionTransient(chargeConv, actStep, actTime);
      }

#ifdef ADAPTGRID
      if (flags->CalcErrorMap_ ) {

        ElemStoreSol<Double> error, error_Mesh;

        std::cout << "Do write" << std::endl;

        // this is only a temporar solution
        error.SetNumSolutions(1);
        error.SetNumElems(errorMap_.GetSize());
        error.SetSolutionType(NO_SOLUTION_TYPE);
        error.SetNumDofs(dofspernode_);
        error.Init(0);

        //Error.SetAlgSysVector(errorMap_);
          
        // ATTENTION!!
        // up to now now transformation of the Error performed,
        // since the calculation of the error is done on the global
        // element numeration
        // Error.TransformElemSolution(Error_Mesh,subdoms_,ptgrid_);
        // OutFile_->WriteElemSolution(errorMap_, laststepcalc_, time,
        // "relERR-E-Potential"); 
        outFile_->WriteElemSolutionTransient(error_Mesh, actStep, actTime); 
      }
#endif

     
    }

    else {

      // Down-cast
      NodeStoreSol<Complex> *mySol = NULL;
      mySol = dynamic_cast< NodeStoreSol<Complex>* >( sol_ );

      // Write electric potential
      if ( saveSol_ == TRUE ) {
        outFile_->WriteNodeSolutionHarmonic( *mySol, actStep, actTime,
                                             complexFormat_ );
      }
      
      if ( calcEfield_.GetSize() != 0 ) {
        ElemStoreSol<Complex> & eFieldConv = 
          dynamic_cast<ElemStoreSol<Complex>& >(*E_);
        outFile_->WriteElemSolutionHarmonic( eFieldConv, actStep, 
                                             actTime, complexFormat_ );
      }

      if ( calcCharges_.GetSize() != 0 ) {
        ElemStoreSol<Complex> & chargeConv = 
          dynamic_cast<ElemStoreSol<Complex>& >(*charges_);
        outFile_->WriteElemSolutionHarmonic( chargeConv, actStep, 
                                             actTime, complexFormat_ );
      }
    }


    // The following section was used by Gerhard to compute sum of forces over
    // different iteration and time steps. The sum was written into the .data
    // stream. Since this is not available anymore, this is commented out
#ifdef COMMENTET_OUT
    if (isIterCoupled_ == TRUE) {
      //   // TMPORARILY
      SolutionType quantity;
      StdVector<UInt> * couplingNodes     = NULL;
      CFSVector * values = 0;
      Vector<Double> sumForces(dim_);
      sumForces.Init();
    
      // loop over all output coupling quantities
      for (UInt actCoupling=0;
           actCoupling<ptCoupling_->GetNumOutputCouplings(); actCoupling++) {
        quantity = ptCoupling_->GetOutputQuantity(actCoupling);
        ptCoupling_->GetOutputValues(actCoupling, values);
        
        Vector<Double> const &temp = dynamic_cast<Vector<Double> &>(*values);
        switch(ptCoupling_->GetOutputType(actCoupling)) {
            
        case NODE:      
          ptCoupling_->GetOutputNodes(actCoupling, couplingNodes);
          if (quantity == ELEC_FORCE_VWP) {
            for (UInt iDof=0; iDof<dim_; iDof++)
              for (UInt iNode=0; iNode<couplingNodes->GetSize(); iNode++)
                sumForces[iDof] += temp[iNode*dim_ + iDof];

            *data << lasttimecalc_ << "\t";
            for (UInt i=0; i<dim_; i++)
              *data << sumForces[i]<< "\t";

            *data << std::endl;
          }
          break;

        case ELEM:
          Error( "Element input coupling not implemented for elecPDE",
                 __FILE__, __LINE__ );
        } // switch
      } // for
    }
#endif
  }


  // **********************
  //   WriteHistoryInFile
  // **********************
  void ElecPDE::WriteHistoryInFile( const UInt kstep,
                                    const Double asteptime,
                                    UInt stepOffset,
                                    Double timeOffset ) {

    ENTER_FCN( "ElecPDE::WriteHistoryInFile" );

    UInt actStep = kstep + stepOffset;
    Double actTime = timeOffset + asteptime;

    // only first process should write output
#ifdef PARALLEL
    int commrank;
    MPI_Comm_rank( MPI_COMM_WORLD, &commrank );
    if ( commrank != 0 ) {
      return;
    }
#endif
    
    if ( analysistype_ == STATIC || analysistype_ == TRANSIENT ) {

      // Down-cast
      NodeStoreSol<Double> *solConverted;
      solConverted = dynamic_cast<NodeStoreSol<Double>*>(sol_);

      if ( saveSolHist_ ) {
        outFile_->WriteNodeHistoryTransient(*solConverted, actStep, actTime);
      }

    }
    else {

      // Down-cast
      NodeStoreSol<Complex> *mySol = NULL;
      mySol = dynamic_cast< NodeStoreSol<Complex>* >( sol_ );

      // Write electric potential
      if ( saveSolHist_ == TRUE ) {
        outFile_->WriteNodeHistoryHarmonic( *mySol, actStep, actTime,
                                            complexFormat_ );
      }
    }
  }



  // ***************
  //   PostProcess
  // ***************
  void ElecPDE::PostProcess() {

    ENTER_FCN( "ElecPDE::PostProcess" );

    // *** Electric Field Intensity ***
    if (calcEfield_.GetSize() !=0 ){
      if (analysistype_ == TRANSIENT ||
          analysistype_ == STATIC ) {        
        CalcElectricField<Double>();
      } else {
        CalcElectricField<Complex>();
      }
    }
  
    // *** Electric Charges ***
    if (calcCharges_.GetSize() !=0 ) {
      if (analysistype_ == TRANSIENT ||
          analysistype_ == STATIC ) {        
        CalcCharges<Double>();
      } else {
        CalcCharges<Complex>();
      }
    }
    // *** Electric Energy ***
    if ( calcEnergy_.GetSize() != 0 ) {
      if (analysistype_ == TRANSIENT ||
          analysistype_ == STATIC ) {        
        CalcEnergy<Double>();
      } else {
        CalcEnergy<Complex>();
      }
    }     
    

#ifdef TCL_INTERFACE
    StdVector<std::string> context;
    context.Push_back( pdename_ );
    context.Push_back( GenStr(solveStep_->GetActStep() ) );
    
    if ( analysistype_ == TRANSIENT ||
         analysistype_ == STATIC ) {
      context.Push_back( GenStr(solveStep_->GetActTime() ) );
    } else {
      context.Push_back( GenStr(solveStep_->GetActFreq() ) );
    }
    messenger->TriggerEvent( CFSMessenger::CFS_PostProcess, 
                             context );
#endif   
  }
  
  template <class TYPE>
  void ElecPDE::CalcElectricField()
  {
    ENTER_FCN( "ElecPDE::PostProcess" );
  

    Vector<Double> lCoord;
    StdVector<Elem*> elemssd;
    UInt counterElems=0;
    UInt pdeElem;

    
    NodeStoreSol<TYPE> & solhelp = dynamic_cast<NodeStoreSol<TYPE>&>(*sol_);
    Vector<TYPE> tempE ;
    GradientFieldOp<TYPE> * FieldOp = 
        new GradientFieldOp<TYPE>(ptgrid_, this, eqnData_, solhelp, 
                                  ELEC_POTENTIAL, isaxi_);
      // loop over all subdomains
    for (UInt isd=0; isd<calcEfield_.GetSize(); isd++)
      {
        
        // ------ Calculation of the electric field ------
        ptgrid_->GetVolElems( elemssd, calcEfield_[isd] );
        
        // loop over elements of subdomain
        for (UInt iel=0; iel< elemssd.GetSize(); iel++,counterElems++)
          {
            elemssd[iel]->ptElem->GetCoordMidPoint(lCoord);
            FieldOp->CalcElemGradField( tempE, elemssd[iel], lCoord, 1); 
            pdeElem = eqnData_->Mesh2PDEElem(elemssd[iel]->elemNum);
            E_->SetElemResult(pdeElem-1, tempE);
          }
      }
    
    delete FieldOp;
  }

  template <class TYPE>
  void ElecPDE::CalcCharges()
  {
    ENTER_FCN( "ElecPDE::CalcCharges" );

    NodeStoreSol<TYPE> * solhelp = dynamic_cast<NodeStoreSol<TYPE>*>(sol_);
    StdVector<SurfElem*> surfElems;

    Vector<Double> lCoordSurf, lCoordVol, normal;
    Vector<TYPE> elemDField;
    Elem * ptVolElem;
    BaseFE * ptSurfElemFE, * ptVolElemFE;
    Double permittivity = 0.0;
    TYPE elemNormalD = 0.0;
    TYPE charge = 0.0;
    TYPE sumOfCharges = 0.0;
    UInt pdeElemNum = 0;
    Double normSign = 0;
 
    // Check, if charges are written on the surface elements or onto the 
    // volume elements
    std::string outputType = "surface";
    params->Get("outputType",outputType,"electrostatic","charge");

 
    // Create vector with interpolation coordinate.
    // For simplicity we only evaluate the integral
    // in coordinate origin
    lCoordSurf.Resize(dim_-1);
    lCoordSurf.Init(0);
  
    // Create operator for electric flux density and charge calculation
    GradientFieldOp<TYPE> * dFieldOp = 
      new GradientFieldOp<TYPE>(ptgrid_, this, eqnData_, *solhelp,
                                  ELEC_POTENTIAL, isaxi_);

    ElecChargeOp<TYPE> * chargeOp = 
      new ElecChargeOp<TYPE>(ptgrid_, this, eqnData_, isaxi_);
  
    // loop over all subdomains
    for (UInt iSD=0; iSD<calcCharges_.GetSize(); iSD++){
    
      // get surface and acoording volume elements
      ptgrid_->GetSurfElems( surfElems, calcCharges_[iSD] );
    
      // loop over all surface elements
      for (UInt iElem=0; iElem<surfElems.GetSize(); iElem++)
        {
        
          // Determine, which volume element is the right neighbour for the 
          // calculation
          if ( chargeNeighborRegion_.
               Find(surfElems[iElem]->ptVolElem1->regionId) != -1 ) {
            ptVolElem = surfElems[iElem]->ptVolElem1;
            normSign = -1.0;
          } else {
            ptVolElem = surfElems[iElem]->ptVolElem2;
            normSign = 1.0;
          }

          normSign *= (Double) surfElems[iElem]->normalSign;

          ptSurfElemFE = surfElems[iElem]->ptElem; 
          ptVolElemFE = ptVolElem->ptElem;
        
          const StdVector<UInt> & surfConnect = surfElems[iElem]->connect;
          const StdVector<UInt> & volConnect = ptVolElem->connect;
        
          // calculate volume integration coordinates from
          // surfe integration coordinat for evalauting the 
          // electric flux density on the surface of the volume
          // element
          ptSurfElemFE->GetCoordMidPoint(lCoordSurf);
          ptVolElemFE->GetLocalIntPoints4Surface(surfConnect, volConnect,
                                                 lCoordSurf, lCoordVol);


          // Get the right material parameter for actual volume element
//           for (UInt i=0; i<subdoms_.GetSize(); i++)
//             {
//               if (subdoms_[i] == ptVolElem->regionId)
//                 materialData_[i]->GetScalar(permittivity,ELEC_PERMITTIVITY,REAL);
//             }
          BaseMaterial * myMat = materials_[ptVolElem->regionId];
          myMat->GetScalar(permittivity,ELEC_PERMITTIVITY,REAL);

          // Calc electric flux density
          dFieldOp->CalcElemGradField(elemDField, ptVolElem, 
                                      lCoordVol,permittivity);
        
          // Calc global normal
          ptgrid_->CalcSurfNormal(normal, *surfElems[iElem]);

          normal *= normSign;

          // Since the routine CalcLineNormal always computes a normal
          // which points in the OPPOSITE direction of the volume element,
          // we have to multiply the normal with -1 to get the correct sign for
          // the charges
          elemNormalD = 0.0;
          for ( UInt iComp = 0; iComp < normal.GetSize(); iComp++ ) {
            elemNormalD =  elemDField[iComp] * normal[iComp];
          }
          
          chargeOp->CalcElemCharge(charge, surfElems[iElem], 
                                   lCoordSurf, elemNormalD);

          if ( outputType == "surface" ) {
            pdeElemNum = eqnData_->Mesh2PDEElem(surfElems[iElem]->elemNum);
          } else {
            pdeElemNum = eqnData_->Mesh2PDEElem(ptVolElem->elemNum);
          }
        
          // Create temporary vector, since SetElemResult only
          // can handle these
          Vector<TYPE> chargeVec(1);
          chargeVec[0] = charge;
          sumOfCharges +=charge;
          charges_->SetElemResult(pdeElemNum-1, chargeVec);
        
        }
    }
    std::string outstring = "Sum of electric charges:\n";
    outstring += GenStr(sumOfCharges);
    Info->PrintF(pdename_, outstring.c_str());

    delete chargeOp;
    delete dFieldOp;
  }

  template <class TYPE>
  void ElecPDE::CalcEnergy()
  {
    ENTER_FCN( "ElecPDE::CalcEnergy" );

    Matrix<Double> elemmat;  
    Matrix<Double> ptCoord;
    BaseFE         * ptElem;
    TYPE totalE = 0.0;

    SubTensorType tensorType;

    if ( dim_ == 3 ) {
      tensorType = FULL;
    }
    else {
      if ( isaxi_ == TRUE ) {
	tensorType = AXI;
      }
      else {
	// 2d: plane case
	tensorType = PLANE_STRAIN;
      }
    }


    StdVector<UInt> connecth;
    Vector<TYPE> help;
    Double factor = 1.0;
    UInt i, j;
    Vector<TYPE> energy(subdoms_.GetSize());

    for (i=0; i<subdoms_.GetSize(); i++)
      {
        
        StdVector<Elem*> elemssd;
        ptgrid_->GetVolElems( elemssd,subdoms_[i] );

        energy[i] = 0;
        for (j=0; j < elemssd.GetSize(); j++)
          {  
            ptElem=elemssd[j]->ptElem;
            BaseForm * bilinear_stiff = 
              new linElecInt(materials_[subdoms_[i]],tensorType);
            bilinear_stiff->SetFactor(factor);
            bilinear_stiff->SetElemPtr(ptElem);

            connecth=elemssd[j]->connect;
            GetElemCoords(connecth, ptCoord);
            bilinear_stiff->CalcElementMatrix(ptCoord, elemmat);

            // Mape Mesh to PDE node numbers
            //Mesh2PDENode(connect_PDE,connecth,mesh2PDENode_);

            //        EqnData_->Mesh2Eqn(Eqns,connecth);
            //        (*debug) << "Nodes:" << connecth << std::endl;
            //        (*debug) << "Eqns :" << Eqns << std::endl;


            Vector<TYPE> elpot;
            sol_->GetElemSolution(elpot, connecth);
            help =  elemmat * elpot;
            energy[i] += 0.5 * (help * elpot);

            delete bilinear_stiff;          
          }  

        totalE += energy[i];
      }

    std::string unit = "(Ws)";
    std::string resulttype = "Electric Energy";
    std::string analysis;
    Double analysisVal;
    if ( analysistype_ == HARMONIC ) {
      analysis    = "Frequency:";
      analysisVal = solveStep_->GetActFreq();
    }
    else {
      analysis    = "Time:";
      analysisVal = solveStep_->GetActTime();
    }

    StdVector<std::string> subdomNames;
    ptgrid_->RegionIdToName( subdomNames, subdoms_ );

    Info->WriteResult(pdename_,  resulttype, subdomNames, energy, unit,
                      analysis, analysisVal);


    StdVector<std::string> suball(1);
    Vector<TYPE> tmp(1);
    suball[0] = "Sum";
    tmp[0] = totalE;
    Info->WriteResult(pdename_,  resulttype, suball, tmp, unit,
                      analysis, analysisVal);

  }




  // ======================================================
  // COUPLING SECTION
  // ======================================================



  void ElecPDE::InitCoupling(PDECoupling * Coupling)
  {
    ENTER_FCN( "ElecPDE::InitCoupling" );
  
    isIterCoupled_ = TRUE;
    ptCoupling_   = Coupling;
    
    StdVector<std::string> * nRegions;
    StdVector<RegionIdType> nRegionIds;
    const UInt numCouplings = ptCoupling_->GetNumOutputCouplings();
    
    nonLin_ = FALSE;

    // Initialization of coupling helper arrays
    std::string quantity;
    StdVector<UInt> * couplingnodes = NULL;

    for (UInt actCoupling=0; actCoupling<numCouplings; actCoupling++) {
      // Initialize arrays for coupling surface elements
      if (ptCoupling_->GetOutputQuantity(actCoupling) == ELEC_FORCE_VWP
          || ptCoupling_->GetOutputQuantity(actCoupling) == 
          ELEC_INTERFACE_FORCE) {
      
        ptCoupling_->GetOutputNodes(actCoupling, couplingnodes);
        if (couplingnodes == NULL)
          std::cerr << "Couplingnodes = 0!!!!" << std::endl;
      
        // if quantity is elecFocrceVWP, get volume neighbours lying next to
        // coupling nodes, because these volume elements have to be 
        // moved 'virtually'
        if (ptCoupling_->GetOutputQuantity(actCoupling) == ELEC_FORCE_VWP) {
          NodeStoreSol<Double> * solhelp = 
            dynamic_cast<NodeStoreSol<Double> *>(sol_);
          ForceOp_ = new  ElecForceOp(ptgrid_, this, eqnData_, 
                                      *solhelp, dim_, materials_, isaxi_);

          ptCoupling_->GetOutputNeighbourRegion(actCoupling, nRegions);
          ptgrid_->RegionNameToId(nRegionIds,*nRegions);
          ForceOp_->Setup(nRegionIds, *couplingnodes);
        }
      
        else if ( ptCoupling_->GetOutputQuantity(actCoupling)
                  == ELEC_INTERFACE_FORCE) {
          Error( "Currently ELEC_INTERFACE_FORCE not supported", __FILE__,
                 __LINE__ );
        }

        else {
          Enum2String(ptCoupling_->GetOutputQuantity(actCoupling), quantity);
          std::string errMsg = "Coupling " + quantity +  " not known! ";    
          Error(errMsg.c_str(), __FILE__,__LINE__);
        }
      
        // Intialize the memory of the coupling values
        ptCoupling_->CreateCouplingVector(actCoupling,isComplex_);
      
      
      } // end for (actNode)
    } 
  }
  


  void ElecPDE::CalcOutputCoupling()
  {
    ENTER_FCN( "ElecPDE::CalcOutputCoupling" );

    SolutionType quantity;
    StdVector<UInt> * couplingNodes     = NULL;
    CFSVector * values = 0;
    UInt forcesCount = 0;

    // at first, check if this PDE is iterative coupled
    if (isIterCoupled_ == FALSE)
      return;

    // loop over all output coupling quantities
    for (UInt actCoupling=0; 
         actCoupling<ptCoupling_->GetNumOutputCouplings(); 
         actCoupling++)
      {
        quantity = ptCoupling_->GetOutputQuantity(actCoupling);
        ptCoupling_->GetOutputValues(actCoupling, values);

        Vector<Double> * temp = dynamic_cast<Vector<Double> *>(values);
      
        switch(ptCoupling_->GetOutputType(actCoupling))
          {
          
          case NODE:        
            ptCoupling_->GetOutputNodes(actCoupling, couplingNodes);
          
            if (quantity == ELEC_POTENTIAL)
              sol_->NodeSolutionToCoupling(*values, *couplingNodes);
            
            if (quantity == ELEC_FORCE_VWP)
              {
                Vector<Double> totalForce;
                ForceOp_->CalcNodeForce(*temp, totalForce);

                // write information in .info-file
                Info->PrintF(pdename_, "Sum of electrostatic force (VWM):\n");
                Info->PrintVec(totalForce);
                forcesCount++;
              }
                  
            break;
          
          case ELEM:
            Error("No Element coupling output", __FILE__,__LINE__);
          }
      }
  }


  Boolean ElecPDE::HasOutput(SolutionType output)
  {
    ENTER_FCN( "ElecPDE::HasOutput" );
  
    switch (output)
      {
      case ELEC_FORCE_VWP:
        return TRUE;
        break;
      case ELEC_POTENTIAL:
        return TRUE;
        break;
      case ELEC_FIELD_INTENSITY:
        return TRUE;
        break;
      case ELEC_INTERFACE_FORCE:
        return TRUE;
        break;
      default:
        return FALSE;
        break;
      }
    return FALSE;
  }


  void ElecPDE::SetPiezoCoupling()
  {
    ENTER_FCN( "ElecPDE::SetPiezoCoupling" );
  
    isPiezoCoupled_ = TRUE;

  }

  void ElecPDE::AssembleSpecial( ) {

    ENTER_FCN( "ElecPDE::AssembleSpecial" );

    Double omega, temp;
    Complex factorR, factorL, factorC, factorAll;

    // Iterate over all impedances
    for ( UInt i = 0; i < impedances_.GetSize(); i++ ) {
      
      // Calculate parameters
      omega = solveStep_->GetActFreq() * 2.0 *PI;

      // *** Resistance ***
      if ( impedances_[i].resistance > EPS ) {
        temp = - ( omega / (omega*omega * impedances_[i].resistance) );
      } else {
        temp = 0.0;
      }
      factorR = Complex( 0, temp );
      
      // *** Inductance ***
      if ( impedances_[i].inductance > EPS ) {
        temp = - 1.0 / ( omega * omega * impedances_[i].inductance );
      } else {
        temp = 0.0;
      }
      factorL = Complex( temp, 0.0 );
      
      // *** Capcitance ***
      factorC = Complex( impedances_[i].capacitance, 0.0 );
      
      factorAll = factorR + factorL + factorC;
      
      // Set Factors
      Integer eqn1, eqn2;
      UInt eqnDof1, eqnDof2;
      Complex oldEntry;
      eqnData_->Node2EQN( impedances_[i].node1, 1, eqn1, eqnDof1);
      eqnData_->Node2EQN( impedances_[i].node2, 1, eqn2, eqnDof2);
      
      // -- Diagonal part --
      algsys_->GetMatrixEntry( SYSTEM, pdeId_, eqn1, 1, pdeId_, eqn1, 1,
                               oldEntry);
      algsys_->SetMatrixEntry( SYSTEM, pdeId_, eqn1, 1, pdeId_, eqn1, 1,
                               oldEntry - factorAll, true);      
      
      algsys_->GetMatrixEntry( SYSTEM, pdeId_, eqn2, 1, pdeId_, eqn2, 1,
                               oldEntry);
      algsys_->SetMatrixEntry( SYSTEM, pdeId_, eqn2, 1, pdeId_, eqn2, 1,
                               oldEntry - factorAll, true);      

      // -- Off diagonal part --
      algsys_->GetMatrixEntry( SYSTEM, pdeId_, eqn1, 1, pdeId_, eqn2, 1,
                               oldEntry );
      algsys_->SetMatrixEntry( SYSTEM, pdeId_, eqn2, 1, pdeId_, eqn2, 1,
                               oldEntry + factorAll, true );
    }



  }
  
  void ElecPDE::SetupMatrixGraphSpecial( ) {
    ENTER_FCN( "ElecPDE::SetupMatrixGraphSpecial" );

    Integer eqn1, eqn2;
    UInt eqnDof1, eqnDof2;
    StdVector<Integer> connect;
    
    for ( UInt i = 0; i < impedances_.GetSize(); i++ ) {
      
      eqnData_->Node2EQN( impedances_[i].node1, 1, eqn1, eqnDof1);
      eqnData_->Node2EQN( impedances_[i].node2, 1, eqn2, eqnDof2);
      
      // Set element positions
      connect.Clear();
      connect.Push_back(eqn1);
      connect.Push_back(eqn2);
      
      algsys_->SetElementPos( pdeId_, connect.GetPointer(), 2,
                              pdeId_, connect.GetPointer(), 2,
                              true);  
    }

    
  }

  // ***********************************************************************
  //   Obtain information on desired output quantities from parameter file
  // ***********************************************************************
  void ElecPDE::ReadStoreResults() {

    ENTER_FCN( "ElecPDE::ReadStoreResults" );

    // Construct vectors for restricted parameter search
    StdVector<std::string> keyVec;
    StdVector<std::string> attrVec;
    StdVector<std::string> valVec;
    StdVector<std::string> regionNames;
    std::string quantity;

    // *****************************
    // Determine nodal results
    // *****************************

    // --- Electric Potential ---
    StdVector<std::string> nodeValues;
    Enum2String(ELEC_POTENTIAL, quantity);
    keyVec  = pdename_, "storeResults", "nodeResults", "region";
    attrVec = "", "", "type";
    valVec = "", "", quantity;
    params->GetList( keyVec, attrVec, valVec, nodeValues);
    if (nodeValues.GetSize() > 0) {
      saveSol_ = TRUE;
      hasOutput_ = TRUE;
    }
  
    // *****************************
    // Determine element results
    // *****************************
    keyVec  = pdename_, "storeResults", "elemResults", "region";
    attrVec = "", "", "type";

    // --- Electric Field Intensity ---
    Enum2String(ELEC_FIELD_INTENSITY, quantity);
    valVec  = "", "", quantity;
    params->GetList( keyVec, attrVec, valVec, regionNames );
    ptgrid_->RegionNameToId( calcEfield_, regionNames );

    // If the the symbolic name is "all" compute electric field for all regions
    if ( calcEfield_.GetSize() == 1 && calcEfield_[0] == ALL_REGIONS ) {
      calcEfield_ = subdoms_;
    }

    if ( calcEfield_.GetSize() > 0 ) {
      hasOutput_ = TRUE;
      Info->PrintF( pdename_, " Computing electric field for regions:\n" );
      for ( UInt k = 0; k < regionNames.GetSize(); k++ ) {
        Info->PrintF( pdename_, " %s\n", regionNames[k].c_str() );
      }
      if ( !isComplex_ ) {
        E_ = new ElemStoreSol<Double>;
          } else {
        E_ = new ElemStoreSol<Complex>;
      }
      E_->SetNumSolutions(1);
      E_->SetSolutionType(ELEC_FIELD_INTENSITY);
      E_->SetNumElems(numElems_);
      E_->SetNumDofs(dim_);
      E_->SetPtrEQNData(eqnData_, ptgrid_);
      E_->Init(); 
    }

    // --- Electric Energy ---
    Enum2String(ELEC_ENERGY, quantity);
    valVec  = "", "", quantity;
    params->GetList( keyVec, attrVec, valVec, regionNames );
    ptgrid_->RegionNameToId( calcEnergy_, regionNames );

    // If the the symbolic name is "all" compute energy for all regions
    if ( calcEnergy_.GetSize() == 1 && calcEnergy_[0] == ALL_REGIONS ) {
      calcEnergy_ = subdoms_;
    }

    if ( calcEnergy_.GetSize() > 0 ) {
      hasOutput_ = TRUE;
      Info->PrintF( pdename_, " Computing energy for regions:\n" );
      for ( UInt k = 0; k < regionNames.GetSize(); k++ ) {
        Info->PrintF( pdename_, " %s\n", regionNames[k].c_str() );
      }
    }
  
    // --- Electric Charges ---
    //check for charge computation
    StdVector<std::string> temp;
    params->GetList( "region", regionNames, pdename_, "charge" );
    params->GetList( "element", temp, pdename_, "charge" );

    ptgrid_->RegionNameToId( chargeNeighborRegion_, regionNames );
    ptgrid_->RegionNameToId( calcCharges_, temp );

    if (calcCharges_.GetSize() > 0)
      {
        
        // Check, if Piezo-Coupling is active
        if ( isPiezoCoupled_ == TRUE ) {
          (*warning) << "The electrostatic PDE is directly piezo-coupled.\n"
                     << "Therefore the computation of charges in the "
                     << "electrostatic part will yield not the correct result!";
            Warning( __FILE__, __LINE__ );
        }
        
        hasOutput_ = TRUE;
        Info->PrintF( pdename_,
                      " Computing electric charges for regions:\n");
        for ( UInt k = 0; k < temp.GetSize(); k++ ) {
          Info->PrintF( pdename_, " %s\n", temp[k].c_str() );
        } 
        if ( !isComplex_ ) {
          charges_ = new ElemStoreSol<Double>;
        } else {
          charges_ = new ElemStoreSol<Complex>;
        }
        // Resize solution arrays
        charges_->SetNumSolutions(1);
        charges_->SetSolutionType(ELEC_CHARGE);
        charges_->SetNumElems(numElems_);
        charges_->SetNumDofs(1);
        charges_->SetPtrEQNData(eqnData_, ptgrid_);
        charges_->Init();
      } 
  
    // *****************************
    // Determine nodal history
    // *****************************
    StdVector<std::string> saveNodeHist;
    keyVec  = pdename_, "storeResults", "nodeHistory", "saveNodes";
    attrVec = "", "", "type";

    // --- Electric Potential ---
    Enum2String(ELEC_POTENTIAL, quantity);
    valVec  = "", "", quantity;
    params->GetList( keyVec, attrVec, valVec, saveNodeHist );
  
    if (saveNodeHist.GetSize() > 0) {
      saveSolHist_ = TRUE;
      hasOutput_ = TRUE;
      Info->PrintF( pdename_, " Saving ElecPotential for Nodes:\n" );
      for ( UInt k = 0; k < saveNodeHist.GetSize(); k++ ) {
        Info->PrintF( pdename_, " %s\n", saveNodeHist[k].c_str() );
      }
    }
  
    // *****************************
    // Determine element history
    // *****************************
    StdVector<std::string> saveElemHist;
    keyVec  = pdename_, "storeResults", "elemHistory", "saveElems";
    attrVec = "", "", "";
    valVec = "", "", "";
    params->GetList(keyVec, attrVec, valVec, saveElemHist);
  
    if (saveElemHist.GetSize() > 0) {
      std::string errMsg = pdename_;
      errMsg += ": Saving history elements is not implemented yet!\n";
      errMsg += "Meanwhile you can use 'unvtool' to extract element data.";
      Error( errMsg.c_str(), __FILE__, __LINE__);
    }
  
  }


} // end of namespace

