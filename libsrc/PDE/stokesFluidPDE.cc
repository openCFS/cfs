#include "stokesFluidPDE.hh"

#include <sstream>
#include <iomanip>

#include "Forms/forms_header.hh"
#include "Forms/linElastInt.hh"
#include "Forms/massInt.hh"
#include "Forms/linPressureInt.hh"
#include "DataInOut/writeresults.hh"
#include "Driver/assemble.hh"
#include "newmark.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "Driver/solveStepStokesFluid.hh" 
#include "CoupledPDE/pdecoupling.hh"
#include "Domain/domain.hh"
#include "Utils/coordSystem.hh"

#ifdef TCL_INTERFACE
#include "DataInOut/Scripting/cfsmessenger.hh" 
#endif

namespace CoupledField
{

  StokesFluidPDE::StokesFluidPDE(Grid * aptgrid, TimeFunc *aptTimeFunc, WriteResults *aptOut)
    :SinglePDE(aptgrid, aptOut, aptTimeFunc)

  {
    ENTER_FCN( "StokesFluidPDE::StokesFluidPDE" );

    isAlwaysStatic_ = TRUE;

    pdename_          = "stokesFluid";
    pdematerialclass_ = "fluid";


    // ****************************
    // DETERMINE GEOMETRY
    // ****************************

    // Get problem geometry and PDE subtype
    params->Get( "subType", subType_, pdename_ );
    std::string probGeo;
    params->Get( "type", probGeo, "geometry" );

    // Set number of degrees of freedom and
    // ensure that subtype fits to problem geometry
    if ( subType_ == "3d" && probGeo == "3d" ) {
      dofspernode_ = 8;
      Info->PrintF("", "=== 3D PROBLEM\n");
    }
    else if ( subType_ == "axi" && probGeo == "axi" ) {
      isaxi_ = TRUE;
      dofspernode_ = 4;
      Info->PrintF("", "=== AXISYSMMETRIC PROBLEM\n");
    }
    else if ( subType_ == "plane" && probGeo == "plane" ) {
      dofspernode_ = 4;
      Info->PrintF("", "=== PLANE PROBLEM\n");
    }
    else
      {
        std::string errmsg = "Subtype '" + subType_;
        errmsg += "' of PDE '" + pdename_ + "' does not fit to problem ";
        errmsg += "geometry '" + probGeo + "'\n";
        Info->Error( errmsg, __FILE__, __LINE__ );
      }

    // =====================================================================
    // set solution information
    // =====================================================================
    solTypes_ = STOKESFLUID_VEL_PRES_VORT;
    solDofs_ = dofspernode_;
   }

  StokesFluidPDE::~StokesFluidPDE()
  {
    ENTER_FCN( "StokesFluidPDE::~StokesFluidPDE" );
  }

  
  void StokesFluidPDE::InitNonLin()
  {
    ENTER_FCN( "StokesFluidPDE::InitNonLin");

    // ==============================================================
    // NOTE: Currently we can only treat convective non-linearity and
    //       we assume that for a stokes fluid PDE all regions either
    //       are linear or non-linear!
    // ==============================================================
    StdVector<std::string> nonLinRegion;
    params->GetList( "nonLinear", nonLinRegion, pdename_, "region" );
    // Should not happen with validating parser, but beware!
    if ( nonLinRegion.GetSize() == 0 ) {
      nonLin_ = FALSE;
    }
    else {
      for ( UInt k = 1; k < nonLinRegion.GetSize(); k++ ) {
        if ( nonLinRegion[k] != nonLinRegion[0] ) {
          Info->Error( "Non-linearity should be the same for all regions!",
                       __FILE__, __LINE__ );
        }
      }
      nonLin_ = nonLinRegion[0] == "convection" ? TRUE : FALSE;
    }

    // In non-linear case determine type of line search strategy
    if ( nonLin_ == TRUE ) {

      Info->PrintF( pdename_, "Non-linearity in %d regions\n",
                    nonLinRegion.GetSize() );

      // type of line search
      params->Get( "type", lineSearch_, pdename_, "lineSearch" );

      if ( lineSearch_ == "none" ) {
        Info->PrintF( pdename_, "Performing no line search\n" );
      }
      else {
        Info->PrintF( pdename_, "Will perform line search\n" );
      }

    }

    // If no non-linearity we do not perform line search anyhow
    else {
      lineSearch_ = "none";
    }

    if( nonLin_ == TRUE )
      {
        // incremental stopping criterion
        params->Get( "incStopCrit", incStopCrit_, pdename_, "nonLinear" );

        // residual stopping criterion
        params->Get( "resStopCrit", residualStopCrit_, pdename_, "nonLinear" );
        
        // maximal number of NL-iterations
        params->Get("maxNumIters", nonLinMaxIter_, pdename_, "nonLinear");
      }
  }

  void StokesFluidPDE::DefineIntegrators()
  {
    ENTER_FCN( "StokesFluidPDE::DefineIntegerators" );

    // help variables for parameter checking
    StdVector<std::string> keyVec;
    StdVector<std::string> attrVec;
    StdVector<std::string> valVec;
    
    Double density, dynamicViscosity;

    for (UInt actSD = 0; actSD < subdoms_.GetSize(); actSD++)
      {

      density          = materialData_[actSD].GetDensity();
      dynamicViscosity = materialData_[actSD].GetCompressibility();
      
      Info->PrintF( pdename_, "density = %e\n", density);
      Info->PrintF( pdename_, "dynamic Viscosity = %e\n", dynamicViscosity);
      
      // ==============  add stiffness ======================================
      BaseForm * bilinearStiff = GetStiffIntegrator(density, dynamicViscosity);

      IntegratorDescriptor * actIntDescr =
        new IntegratorDescriptor(bilinearStiff, STIFFNESS);

	  actIntDescr->SetPDEIds(this, this);
    
      assemble_->AddIntegrator(actIntDescr, subdoms_[actSD]);

      // ==============  add nonlinear stiffness ============================
      if (nonLin_)
        {
          BaseForm *nLinPart = NULL;

          if (subType_ == "3d")
            {   
              nLinPart = 
                new nLinStokesFluid3DInt_Convective(density, dynamicViscosity);    
            }
          else if (subType_ == "plane")
            {
              nLinPart = 
                new nLinStokesFluidPlaneInt_Convective(density, dynamicViscosity);    
            }
          else if (subType_ == "axi")
            {
              Info->Error("StokesFluid with convection needs to be implemented for axisymmetry! ",__FILE__,__LINE__);
              //nLinPart =
              //  new nLinStokesFluidAxiInt_Convective(density, dynamicViscosity);    
            }
          IntegratorDescriptor * stiffNLDescr = 
	        new IntegratorDescriptor(nLinPart, STIFFNESS, nonLin_);

	      stiffNLDescr->SetPDEIds(this, this);
	      assemble_->AddIntegrator(stiffNLDescr, subdoms_[actSD]);
        }
        // ==================== RHS ===========================================
        if (nonLin_)
          {
            BaseForm * rhsSource = 
              new nLinStokesFluid_linFormInt(density, dynamicViscosity, isaxi_);
            assemble_->AddRhsIntegrator(rhsSource, subdoms_[actSD], nonLin_);
          }
      }
  }


  BaseForm *
  StokesFluidPDE::GetStiffIntegrator(Double density, Double dynamicViscosity)
  {

    ENTER_FCN( "StokesFluidPDE::GetStiffIntegrator" );
  
    BaseForm * bilinearStiff = NULL;

        if (subType_ == "plane")
          {
            bilinearStiff = new StokesFluidPlaneInt(density, dynamicViscosity);
          }
        else if (subType_ == "axi")
          {
            Info->Error("StokesFluid for axisymmetry needs to be implemented! ",__FILE__,__LINE__);
            //bilinearStiff = new StokesFluidAxiInt(density, dynamicViscosity);
          }
        else if (subType_ == "3d")
          {
            bilinearStiff = new StokesFluid3DInt(density, dynamicViscosity);
          }
        else 
          Info->Error("Unknown subtype in stokesFluid PDE! ",__FILE__,__LINE__);
    
    return bilinearStiff;
  }

 
  void StokesFluidPDE::DefineSolveStep()
  {
    ENTER_FCN( "StokesFluidPDE::DefineSolveStep" );

    solveStep_ = new SolveStepStokesFluid(*this);
  }





  // ======================================================
  // COUPLING SECTION
  // ======================================================


  void StokesFluidPDE::InitCoupling(PDECoupling * Coupling)
  {
    ENTER_FCN( "StokesFluidPDE::InitCoupling" );

    isIterCoupled_ = TRUE;
    ptCoupling_   = Coupling;

    const UInt numCouplings = ptCoupling_->GetNumOutputCouplings();

    StdVector<StdVector<UInt> > elemNodeToCouplingNode_tmp;
    elemNodeToCouplingNode_.Resize(numCouplings);
    
    for (UInt i = 0; i < numCouplings; i++) {

      // Intialize the memory of the coupling values
      ptCoupling_->CreateCouplingVector(i,isComplex_);

      if (ptCoupling_->GetOutputQuantity(i) == STOKESFLUID_FORCE) {
        // now since we need a incremental formulation, 
        //  initialize some necessary vectors
        isIncrFormulation_ = TRUE;
      }
    }
  }

  void StokesFluidPDE::GetPresSolVecOfElement( Vector<Double>& elemSol,
                                       StdVector<UInt>& connecth ) {

    ENTER_FCN( "StokesFluidPDE::GetPresSolVecOfElement" );

    // stokesFluid pressure of element nodes
    elemSol.Resize(1 * connecth.GetSize());
    elemSol.Init(0);
    Integer eqnNr; 
    UInt eqnDof;
    UInt presDof;
    UInt dofsPerEQN = eqnData_->GetNumDofsPerEQN();

    NodeStoreSol<Double> * solhelp = 
      dynamic_cast<NodeStoreSol<Double>*>(sol_);
    Vector<Double> sol = solhelp->GetAlgSysVector();

    if (subType_ == "plane" || subType_ == "axi")
      {
        presDof = 3;
      }
    else if (subType_ == "3d")
      {
        presDof = 8;
      }
    else 
      Info->Error("Unknown PDE subtype! ",__FILE__,__LINE__);

  
    for(UInt actNode=0; actNode<connecth.GetSize(); actNode++) {
      eqnData_->Node2EQN(connecth[actNode],presDof,eqnNr,eqnDof);
      if (eqnNr!= 0) {
        elemSol[actNode] = sol[eqnDof-1 + dofsPerEQN*(abs(eqnNr-1))];
      }
      else {
        elemSol[actNode] = 0.0;
      }
    }
  }

  void StokesFluidPDE::CalcOutputCoupling()
  {
    ENTER_FCN( "StokesFluidPDE::CalcOutputCoupling" );
    UInt dof;
    SolutionType quantity;
    StdVector<Elem*> * couplingElems = NULL;
    StdVector<UInt> * couplingNodes = NULL;
    CFSVector * temp_values = NULL;
    //UInt regionCount = 0;
  
    // at first, check if this PDE is iterative coupled
    if (isIterCoupled_ == FALSE)
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
        
        ptCoupling_->GetOutputNodes(i, couplingNodes);
        
        if (quantity == STOKESFLUID_FORCE) {

          ptCoupling_->GetOutputElements(i, couplingElems);
          dof = ptCoupling_->GetOutputDof(i);
          CalcMechCouplingRHS(couplingElems, *couplingNodes, *values, dof);
        }
        break;

      case ELEM:
        Error("No Element coupling output", __FILE__,__LINE__);
      }
    }
  }

  void StokesFluidPDE::CalcMechCouplingRHS( StdVector<Elem*> * couplingElems, 
                                            StdVector<UInt> & couplingNodes,
                                            Vector<Double>& elemCouplingSols,
                                            UInt couplingdof ) {
    
    ENTER_FCN( "StokesFluidPDE::CalcMechCouplingRHS" );
    
    Double density = 0.0;
    Double dynamicViscosity = 0.0;
    Double sign = 0.0;
    Integer matIndex = -1;
    Elem * ptVolElem = NULL;
    Matrix<Double> ptCoord, elemMat;
    Vector<Double> sol, forceOnElem, normal;
    
    elemCouplingSols.Init(0.0);
    
    for (UInt actElem=0; actElem<couplingElems->GetSize(); actElem++) {
      
      // Perform cast from volume element to surface element, since
      // mech-acou coupling makes only sense on surface elements
      SurfElem * actCoupleElem = 
        dynamic_cast<SurfElem*> ((*couplingElems)[actElem]);

      if (actCoupleElem == NULL) {
        Error( "No elements found for coupling!", __FILE__, __LINE__ );
      }
      
      BaseFE * ptElem = actCoupleElem->ptElem;
      StdVector<UInt> & connecth = actCoupleElem->connect;
      GetElemCoords(connecth, ptCoord);
      
      // Try to find according region for first neighbouring volume
      // element of the surface element
      matIndex = subdoms_.Find(actCoupleElem->ptVolElem1->regionId);
      
      // If first volume element does not belong to stokesFluid PDE, try the
      // second one
      if ( matIndex == -1 ) {
        matIndex = subdoms_.Find(actCoupleElem->ptVolElem2->regionId);
        ptVolElem = actCoupleElem->ptVolElem2;
        sign = actCoupleElem->normalSign;
      } else {
        ptVolElem = actCoupleElem->ptVolElem1;
        sign = -1.0 * actCoupleElem->normalSign;
      }
      
      if ( matIndex == -1) {
        (*error) << "StokesFluidPDE::CalcMechCouplingRHS: The two volume "
                 << "element neighbours of surface element Nr. "
                 << actCoupleElem->elemNum << " do not belong to my regions!";
        Error( __FILE__, __LINE__ );
      }
      
      // Assign correct density and dynamicViscosity
      density = materialData_[matIndex].GetDensity();
      dynamicViscosity = materialData_[matIndex].GetCompressibility();
      
      BaseForm * bilinear_mass = new MassInt(ptElem, 1.0, isaxi_);
      bilinear_mass->CalcElementMatrix(ptCoord, elemMat);
      delete bilinear_mass;     

      // acoustic
      //GetDerivSolVecOfElement(sol, connecth);

      //stokesFluid      
      GetPresSolVecOfElement(sol, connecth);
     
      
      forceOnElem = elemMat * sol;
      
      // force has to be added on RHS with negative sign
      forceOnElem *= - 1.0;
      
      // the normal vector points outwards of the MECHANICAL domain
      // (see. Kaltenbacher, "Num. Sim. of Mechatr. Act. & Sens." chapter 8.2)
      ptgrid_->CalcSurfNormal(normal, *actCoupleElem);
      normal *= sign;

      for (UInt actNode=0; actNode<ptCoord.GetSizeCol(); actNode++) {
        UInt nodePos = 0;
          
        while(connecth[actNode] != couplingNodes[nodePos] &&
              nodePos < couplingNodes.GetSize()) {
          nodePos++;
        }
          
        for (UInt actDof=0; actDof < couplingdof ; actDof++) {
          elemCouplingSols[nodePos*couplingdof +actDof] += 
            forceOnElem[actNode] * normal[actDof];
        }
      }
    }
  }

  Boolean StokesFluidPDE::HasOutput(SolutionType output)
  {
    ENTER_FCN( "StokesFluidPDE::HasOutput" );

    if (output == STOKESFLUID_FORCE)
      return TRUE;

    return FALSE;
  }



  // ======================================================
  // TIME STEPPING SECTION
  // ======================================================


  void StokesFluidPDE :: InitTimeStepping()
  {
    ENTER_FCN( "StokesFluidPDE::InitTimeStepping" );
  }



  // ======================================================
  // POSTPROCESSING SECTION
  // ======================================================


  void StokesFluidPDE::WriteResultsInFile(const UInt kstep,
                                   const Double asteptime,
                                   UInt stepOffset,
                                   Double timeOffset)
  {
    ENTER_FCN( "StokesFluidPDE::WriteResultsInFile" );

    NodeStoreSol<Double> * solTransient;

    Double actTime = asteptime + timeOffset;
    UInt actStep = kstep + stepOffset;
    
#ifdef WRITE_RHS
    NodeStoreSol<Double> rhs;
    rhs.SetNumSolutions(1);
    rhs.SetNumNodes(numPDENodes_);
    rhs.SetSolutionType(MECH_RHSVAL);
    rhs.SetNumDofs(dim_);
    rhs.SetPtrEQNData(eqnData_, ptgrid_);
    rhs.Init(0.0);
    
    Double *ptRHS;
    algsys_->GetRHSVal( ptRHS );
    rhs.CopyFromAlgSysDataPointer(ptRHS);
    outFile_->WriteNodeSolutionTransient(rhs, actStep, actTime);
#endif

    if ( isComplex_ == FALSE ) {
      solTransient = dynamic_cast<NodeStoreSol<Double>*>(sol_);
    
      if (saveSol_ == TRUE ) 
        outFile_->WriteNodeSolutionTransient(*solTransient, actStep, actTime);
   
    } 
  }


  void StokesFluidPDE::WriteHistoryInFile(const UInt kstep,
                                   const Double asteptime,
                                   UInt stepOffset,
                                   Double timeOffset)
  {
    ENTER_FCN( "StokesFluidPDE::WriteHistoryInFile" );

    NodeStoreSol<Double> * solTransient;

    Double actTime = asteptime + timeOffset;
    UInt actStep = kstep + stepOffset;
    
    if ( isComplex_ == FALSE ) {
      solTransient = dynamic_cast<NodeStoreSol<Double>*>(sol_);
      
      if (saveSolHist_ == TRUE)
        outFile_->WriteNodeHistoryTransient(*solTransient, actStep, actTime);
    }
  }

  // ***********************************************************************
  //   Obtain information on desired output quantities from parameter file
  // ***********************************************************************
  void StokesFluidPDE::ReadStoreResults() {

    ENTER_FCN( "StokesFluidPDE::ReadStoreResults" );

    // Construct vectors for restricted parameter search
    StdVector<std::string> keyVec;
    StdVector<std::string> attrVec;
    StdVector<std::string> valVec;
    StdVector<std::string> regionNames;
    std::string quantity;
  
    // -------------------------
    //  Determine nodal results
    // -------------------------
    StdVector<std::string> nodeValues;
    keyVec  = pdename_, "storeResults", "nodeResults", "region";
    attrVec = "", "", "type";  

    // --- StokesFluidVelocity ---
    Enum2String(STOKESFLUID_VEL_PRES_VORT, quantity);
    valVec = "", "", quantity;
    params->GetList( keyVec, attrVec, valVec, nodeValues);
    if (nodeValues.GetSize() > 0) {
      saveSol_ = TRUE;
      hasOutput_ = TRUE;
    }


    // -------------------------
    //  Determine nodal history
    // -------------------------
    StdVector<std::string> saveNodeHist; 
    keyVec  = pdename_, "storeResults", "nodeHistory", "saveNodes";
    attrVec = "", "", "type";

    // --- StokesFluidVelocity ---
    Enum2String(STOKESFLUID_VEL_PRES_VORT, quantity);
    valVec  = "", "", quantity;
    params->GetList( keyVec, attrVec, valVec, saveNodeHist );
  
    if (saveNodeHist.GetSize() > 0) {
      saveSolHist_ = TRUE;
      hasOutput_ = TRUE;
      Info->PrintF( pdename_, " Saving StokesFluidVelocity for Nodes:\n" );
      for ( UInt k = 0; k < saveNodeHist.GetSize(); k++ ) {
        Info->PrintF( pdename_, " %s\n", saveNodeHist[k].c_str() );
      }
    }
  }


  // ************************************************************
  //   PostProcess
  // ************************************************************

  void StokesFluidPDE::PostProcess() {
    ENTER_FCN( "StokesFluidPDE::PostProcess" );
  }
} // end namespace CoupledField


