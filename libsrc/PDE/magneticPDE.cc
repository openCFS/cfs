#include <fstream>

#include "magneticPDE.hh"

#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "Driver/stdSolveStep.hh"
#include "Utils/Coil.hh"
#include "Utils/SmoothSpline.hh"
#include "Utils/LinInterpolate.hh"
#include "Forms/curlfieldop.hh"
#include "Forms/forms_header.hh"
#include "trapezoidal.hh"
#include "DataInOut/writeresults.hh"
#include "CoupledPDE/pdecoupling.hh"

#ifdef TCL_INTERFACE
#include "DataInOut/Scripting/cfsmessenger.hh" 
#endif

namespace CoupledField {

  // **************
  //  Constructor
  // **************
  MagPDE::MagPDE(Grid * aptgrid, TimeFunc *aptTimeFunc, WriteResults *aptOut)
    :SinglePDE(aptgrid, aptOut, aptTimeFunc) {

    ENTER_FCN( "MagPDE::MagPDE" );

    // =====================================================================
    // set solution information
    // =====================================================================
    dofspernode_ = 1;
    solTypes_ = MAG_POTENTIAL;
    solDofs_ = 1;
    pdename_          = "magnetic";
    pdematerialclass_ = ELECTROMAGNETIC;

    // ---------------------------------------------------------
    //   Get special geometry
    // ---------------------------------------------------------
    isaxi_ = params->HasValue( "type", "axi", "geometry" )
      == TRUE ? TRUE : FALSE;

    // ---------------------------
    //   Set coupling parameters
    // ---------------------------
    deltCoords_.Resize( dim_, numPDENodes_ );
  }


  // *************
  //  Destructor
  // *************
  MagPDE::~MagPDE() {
    ENTER_FCN( "MagPDE::~MagPDE" );
    for ( UInt k = 0; k < coilDef_.GetSize(); k++ ) {
      delete coilDef_[k];
    }
  }


  // *********************************
  //  Read special boundary conitions
  // *********************************

  void MagPDE::ReadSpecialBCs() {

    ENTER_FCN( "MagPDE::ReadSpecialBCs" );

    // --------------------------------------------------------------------
    //   Get information about coils and open files for measurement coils
    // --------------------------------------------------------------------
    ReadCoils();
 
    // -----------------------------
    // Check for permanent magnets
    // -----------------------------
    ReadMagnets();
 
  }
  

  // ****************************
  //  Initialize Nonlinearities
  // ****************************
  void MagPDE::InitNonLin() {

    ENTER_FCN( "MagPDE::InitNonLin" );

    // ------------------------------------
    //   Get information on non-linearity
    // ------------------------------------
    StdVector<std::string> regionNames;

    nonLin_ = FALSE;
    params->GetList( "nonLinear", nonLinType_, pdename_, "region" );
    
    for ( UInt k = 0; k < nonLinType_.GetSize(); k++ ) {
      if ( nonLinType_[k] != "no" ) {
        nonLin_ = TRUE;
        break;
      }
    }

    if( nonLin_ ) {

      // solution method
      params->Get( "method", nonLinMethod_, pdename_, "nonLinear" );

      // perform logging?
      nonLinLogging_ = params->IsSet( "logging", pdename_, "nonLinear" );

      // type of line search
      params->Get( "type", lineSearch_, pdename_, "lineSearch" );
 
      // incremental stopping criterion
      params->Get( "incStopCrit", incStopCrit_, pdename_, "nonLinear" );
      
      // residual stopping criterion
      params->Get( "resStopCrit", residualStopCrit_, pdename_, "nonLinear" );

      // maximal number of NL-iterations
      params->Get("maxNumIters", nonLinMaxIter_, pdename_, "nonLinear");
    }
  }


  // *****************************
  //   Definition of Integrators
  // *****************************
  void MagPDE::DefineIntegrators() {

    ENTER_FCN( "MagPDE::DefineIntegerators" );

    // Loop over all regions this PDE lives on
    for ( UInt actSD = 0; actSD < subdoms_.GetSize(); actSD++ ) {

      // Get reluctivity for this domain and perform consistency check
      Double reluctivity;
      materials_[subdoms_[actSD]]->GetScalar(reluctivity,MAG_RELUCTIVITY,REAL);
 
      if ( nonLinType_[actSD] != "no" ) {

        //read in the BH-curve data and compute the approximation
        std::string nlfnc = materials_[subdoms_[actSD]]->GetNonlinFileName();
        ApproxData *nlinFnc = new SmoothSpline(nlfnc);
        //ApproxData *nlinFnc = new LinInterpolate(nlfnc);
        nlinFnc->CalcBestParameter();
        nlinFnc->CalcApproximation();

        BaseForm *curlcurl2D = new nLinCurlCurlNode2DInt(nlinFnc,reluctivity,
                                                         isaxi_);
        curlcurl2D->SetNonLinMethod(nonLinMethod_);      

	IntegratorDescriptor * stiffDescr = 
	  new IntegratorDescriptor(curlcurl2D, STIFFNESS, TRUE);
	stiffDescr->SetPDEIds(this, this);        
	assemble_->AddIntegrator(stiffDescr,  subdoms_[actSD]);

        // nonlinear RHS linearform!!
        BaseForm * rhsSource = new nLinMagNode2D_linFormInt(nlinFnc,
                                                            reluctivity,
                                                            isaxi_);
        assemble_->AddRhsIntegrator(rhsSource, subdoms_[actSD], TRUE);
      }
      else {
        BaseForm *curlcurl2D = new CurlCurlNode2DInt(reluctivity, isaxi_);
	IntegratorDescriptor * stiffDescr = 
	  new IntegratorDescriptor(curlcurl2D, STIFFNESS, FALSE);
	stiffDescr->SetPDEIds(this, this);        
	assemble_->AddIntegrator(stiffDescr,  subdoms_[actSD]);

        if (nonLin_==TRUE) {
          // for nonlinear RHS linearform we need linear and nonlinear
          // subdomains
          BaseForm * rhsSource = new nLinMagNode2D_linFormInt(reluctivity,
                                                              isaxi_);
          assemble_->AddRhsIntegrator(rhsSource, subdoms_[actSD], TRUE);
        }
      }

      // Mass matrix
      Double conductivity;
      materials_[subdoms_[actSD]]->GetScalar(conductivity,MAG_CONDUCTIVITY,REAL);
      BaseForm *bilinear_mass = new MassInt(conductivity, dofspernode_,isaxi_);

      IntegratorDescriptor * massDescr = 
	new IntegratorDescriptor(bilinear_mass, MASS, FALSE);
	massDescr->SetPDEIds(this, this);        
	assemble_->AddIntegrator(massDescr,  subdoms_[actSD]);

      // If this subdomain is a coil we have to do special things
      for ( UInt coil = 0; coil < coilDef_.GetSize(); coil++ ) {
        if ( subdoms_[actSD] == coilRegionId_[coil] ) {
          Double factor = coilDef_[coil]->value_ /
            coilDef_[coil]->windingCrossSection_;
          BaseForm *coilSource = new VolumeSrcInt( factor, isaxi_ );      
          assemble_->AddRhsSrcIntegrator( coilSource,
                                          subdoms_[actSD],
                                          coilDef_[coil]->dynamicsFile_,
                                          nonLin_ );
        }
      }

      // check, if this subdomain is a permanent magnet
      for ( UInt perm = 0; perm < magnetsDomain_.GetSize(); perm++ ) {
        if ( subdoms_[actSD] == magnetsDomain_[perm] ) {

          if ( dim_  == 3 ) {
            Error( "Permanent magnets not implemented for 3D computation",
                   __FILE__, __LINE__ );
          }

          //change direction of magnetization, so that we can use the
          //standard GlobalDerivatives and obtain: (curl w) . M
          Vector<Double> magnetization(dim_);
          if (isaxi_) {
            magnetization[0] = magnetsOriY_[perm];
            magnetization[1] = -magnetsOriX_[perm];
          }
          else {
            magnetization[0] = -magnetsOriY_[perm];
            magnetization[1] = magnetsOriX_[perm];
          }
          
          std::string fncname = "none";
          Boolean nlin = FALSE;
          BaseForm *permSource = new MagPerm2DInt(magnetization, 
                                                  reluctivity, isaxi_ );

          assemble_->AddRhsSrcIntegrator( permSource,
                                          subdoms_[actSD],
                                          fncname,
                                          nlin );
        }
      }
    }
  }

  void MagPDE::DefineSolveStep()
  {
    ENTER_FCN( "MagPDE::DefineSolveStep" );

    solveStep_ = new StdSolveStep(*this);
  }


  // ======================================================
  // TIME-STEPPING SECTION
  // ======================================================

  void MagPDE :: InitTimeStepping() {
    ENTER_FCN( "MagPDE::InitTimeStepping" );
    UInt rhsSize = eqnData_->GetNumEQNs() *
      eqnData_->GetNumDofsPerEQN();
    
    TS_alg_ = new Trapezoidal( algsys_, rhsSize );
  }



  // ======================================================
  // POSTPROCESSING SECTION
  // ======================================================

  void MagPDE::WriteResultsInFile(const UInt kstep,
                                  const Double asteptime,
                                  UInt stepOffset,
                                  Double timeOffset) {

    ENTER_FCN( "MagPDE::WriteResultsInFile" );

    NodeStoreSol<Double> * solTransient;
    NodeStoreSol<Complex> * solHarmonic;

    Double actTime = asteptime + timeOffset;
    UInt actStep = kstep + stepOffset;
    
    if (analysistype_ == STATIC ||
        analysistype_ == TRANSIENT) {
      solTransient = dynamic_cast<NodeStoreSol<Double>*>(sol_);

      if (saveSol_) {
        outFile_->WriteNodeSolutionTransient(*solTransient, 
                                             actStep, actTime);
      }

      if (calcBfield_.GetSize() !=0 ) {
        outFile_->WriteElemSolutionTransient(B_, actStep, actTime);
      }
      
      if (calcEddy_.GetSize() !=0 ) {
        outFile_->WriteElemSolutionTransient(Jeddy_, actStep, actTime);
      }

      
      if (calcForceVWP_.GetSize() > 0 )
        outFile_->WriteNodeSolutionTransient(Force_, 
                                             actStep, actTime);
        
      
    }

    else {

      if (analysistype_ == HARMONIC) {
        
        solHarmonic = dynamic_cast<NodeStoreSol<Complex>*>(sol_);
        
        if (saveSol_)
          outFile_->WriteNodeSolutionHarmonic(*solHarmonic, actStep,
                                              actTime, complexFormat_);
      }
      else {
        (*error) << "MagPDE: Only static, transient and harmonic results can "
                 << "be written";
        Error( __FILE__, __LINE__ );
      }
    }
  }


  void MagPDE::WriteHistoryInFile(const UInt kstep,
                                  const Double asteptime,
                                  UInt stepOffset,
                                  Double timeOffset) {

    ENTER_FCN( "MagPDE::WriteHistoryInFile" );

    NodeStoreSol<Double> * solTransient;
    NodeStoreSol<Complex> * solHarmonic;

    Double actTime = asteptime + timeOffset;
    UInt actStep = kstep + stepOffset;
    
    if (analysistype_ == STATIC ||
        analysistype_ == TRANSIENT) {
      solTransient = dynamic_cast<NodeStoreSol<Double>*>(sol_);

      if (saveSolHist_) {
        outFile_->WriteNodeHistoryTransient(*solTransient, 
                                            actStep, actTime);
      }     
    }
    else {
      if (analysistype_ == HARMONIC) {        
        if (saveSolHist_)
          outFile_->WriteNodeHistoryHarmonic(*solHarmonic, actStep,
                                             actTime, complexFormat_);
      }
      else {
        (*error) << "MagPDE: Only static, transient and harmonic results can "
                 << "be written";
        Error( __FILE__, __LINE__ );
      }
    }
  }


  // ***************
  //   PostProcess
  // ***************
  void MagPDE::PostProcess() {

    ENTER_FCN( "MagPDE::PostProcess" );

    NodeStoreSol<Double> * solhelp = dynamic_cast<NodeStoreSol<Double>*>(sol_);
  
    if (calcEnergy_.GetSize() !=0 )
      CalcEnergy();

    if (calcBfield_.GetSize() !=0 ) {

      CurlNodeOp * FieldOp = new CurlNodeOp(ptgrid_, this, eqnData_,
                                            *solhelp);
      FieldOp->Set2DType(isaxi_);
 
      // ------ Calculation of the magnetic flux density  ------

      Vector<Double> lCoord;
      StdVector<Elem*> elemssd;
      UInt counterElems=0;
      Vector<Double> tempB;
      UInt pdeElem;
      
      // Resize solution arrays
      B_.SetNumSolutions(1);
      B_.SetSolutionType(MAG_FLUX_DENSITY);
      B_.SetNumElems(numElems_);
      B_.SetNumDofs(dim_);
      B_.SetPtrEQNData(eqnData_, ptgrid_);
      B_.Init(0);
      
      // loop over all subdomains
      for (UInt isd=0; isd<calcBfield_.GetSize(); isd++) {

        // get vector of Elem of subdomain with color: subdoms[isd]
        ptgrid_->GetVolElems(elemssd,calcBfield_[isd]);
          
        // loop over elements of subdomain
        for (UInt iel=0; iel< elemssd.GetSize(); iel++,counterElems++) {
          pdeElem = eqnData_->Mesh2PDEElem(elemssd[iel]->elemNum);
          elemssd[iel]->ptElem->GetCoordMidPoint( lCoord );
          FieldOp->CalcElemCurlNode( tempB, elemssd[iel], lCoord); 
          B_.SetElemResult(pdeElem-1, tempB);
        }
      }
      delete FieldOp;
    }

    if (calcEddy_.GetSize() !=0 ) {

      Vector<Double> lCoord;
      StdVector<Elem*> elemssd;
      Vector<Double> shpFnc, tmp;
      Vector<Double> magVecDeriv1Elem;
      StdVector<UInt> connect;
      StdVector<Integer>  connect_PDE;
      Double conductivity = 0.0;

      UInt counterElems=0;
      UInt pdeElem;

      // dimension hard coded for .unv file!
      Vector<Double> JeddyElem(3);
      JeddyElem.Init();
      
      // Resize solution arrays
      Jeddy_.SetNumSolutions(1);
      Jeddy_.SetSolutionType(MAG_EDDY_CURRENT);
      Jeddy_.SetNumElems(numElems_);

      // dimension hard coded for .unv file!
      Jeddy_.SetNumDofs(3);  
      Jeddy_.SetPtrEQNData(eqnData_, ptgrid_);
      Jeddy_.Init(0);

      // loop over all subdomains
      for (UInt actSD=0; actSD<calcEddy_.GetSize(); actSD++) {

        // get vector of Elem of subdomain with color: subdoms[isd]
        ptgrid_->GetVolElems( elemssd, calcEddy_[actSD] );
          
        // Get the right material parameter for actual subdomain
        for (UInt iSD=0; iSD<subdoms_.GetSize(); iSD++)
          if (subdoms_[iSD] == calcEddy_[actSD])
	    materials_[subdoms_[iSD]]->GetScalar(conductivity,MAG_CONDUCTIVITY,REAL);

        // loop over elements of subdomain
        for ( UInt actEl=0; actEl< elemssd.GetSize();
              actEl++,counterElems++ ) {
          BaseFE * ptEl = elemssd[actEl]->ptElem;
          elemssd[actEl]->ptElem->GetCoordMidPoint( lCoord );
          ptEl->GetShFnc(shpFnc,lCoord);
          pdeElem = eqnData_->Mesh2PDEElem(elemssd[actEl]->elemNum);

          connect = elemssd[actEl]->connect;
          
          GetDerivSolVecOfElement(magVecDeriv1Elem,connect);
          JeddyElem[0] = magVecDeriv1Elem * shpFnc;
          JeddyElem[0] *= -conductivity;
          Jeddy_.SetElemResult(pdeElem-1, JeddyElem);
        }
      }
    }

    if (UIfilename_.size() > 0 && analysistype_ == TRANSIENT) {
      Vector<Double> uiSD;
      ComputeUI(uiSD);
      WriteUI2File(uiSD);
    }

    //check for force computation
    if ( calcForceVWP_.GetSize() > 0 ) {
      Vector<Double> totalForce;
      Vector<Double> ForceValues(dim_*ForceNodes_.GetSize());

      ForceOpVWP_->CalcNodeForce(ForceValues, totalForce);

      // put information into storesol

      UInt actPos =  0;
      Double sum = 0.0;
      for (UInt iNode = 0; iNode < ForceNodes_.GetSize(); iNode++) {
        for (UInt iDof = 0; iDof < dim_; iDof++) {
          std::cerr << "Value before: " 
                    << Force_(ForceNodes_[iNode]-1,iDof) << std::endl;
          std::cerr << "Force(" << ForceNodes_[iNode] << ", " << iDof+1
                    << ") = " << ForceValues[actPos] << std::endl;
          Force_(ForceNodes_[iNode]-1,iDof) = ForceValues[actPos++];
          std::cerr << "Value afterwards: " 
                    << Force_(ForceNodes_[iNode]-1,iDof) << std::endl;
        }
      }


      // write information in .info-file
      Info->PrintF(pdename_, "Sum of magnetic force (VWM):\n");
      Info->PrintVec(totalForce);
    }

    // Last but no least trigger postprocessing fromt within script-file
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


  // **************
  //   CalcEnergy
  // **************
  void MagPDE::CalcEnergy() {

    ENTER_FCN( "MagPDE::CalcEnergy" );

    Matrix<Double> elemmat;  
    Matrix<Double> ptCoord;
    BaseFE         * ptElem;

    StdVector<UInt> connecth;
    StdVector<Integer> Eqns;  
    Vector<double> help;
    BaseForm *bilinear_stiff = NULL;
    UInt i, j;
    Vector<Double> energy(calcEnergy_.GetSize());

    for (i=0; i<calcEnergy_.GetSize(); i++) {

      Double reluctivity;
      materials_[calcEnergy_[i]]->GetScalar(reluctivity,MAG_RELUCTIVITY,REAL);

      // find related region in subdoms_
      Integer sdIndex = subdoms_.Find( calcEnergy_[i] );
      if( sdIndex < 0 ) {
        (*error) << "Region " << ptgrid_->RegionIdToName(calcEnergy_[i])
                 << " is not contained in set of regions for magnetic PDE.";
        Error( __FILE__, __LINE__ );
      }

      // Create stiffness integrator
      if ( nonLinType_[sdIndex] != "no" ) {

        //read in the BH-curve data and compute the approximation
        std::string nlfnc = materials_[subdoms_[sdIndex]]->GetNonlinFileName();
        ApproxData *nlinFnc = new SmoothSpline(nlfnc);
        nlinFnc->CalcBestParameter();
        nlinFnc->CalcApproximation();
        bilinear_stiff = new nLinCurlCurlNode2DInt(nlinFnc,reluctivity,
                                                   isaxi_);

        
        // VERY IMPORTANT: Set nonlinear-method "fixpoint", as otherwise also
        // the Frechet part of the stiffness is calculated!
        bilinear_stiff->SetNonLinMethod( "fixPoint" );
      } else {
        bilinear_stiff = new CurlCurlNode2DInt(reluctivity, isaxi_);
      }

      StdVector<Elem*> elemssd;
      ptgrid_->GetVolElems( elemssd,calcEnergy_[i] );

      energy[i] = 0;
      for (j=0; j < elemssd.GetSize(); j++) {

        ptElem=elemssd[j]->ptElem;
        bilinear_stiff->SetElemPtr( ptElem );

        connecth=elemssd[j]->connect;
        GetElemCoords(connecth, ptCoord);
        Vector<Double> magvecpot;
        sol_->GetElemSolution(magvecpot,connecth);
        if(  nonLinType_[sdIndex] != "no" ) {
          Matrix<Double> magPotMatrix;
          sol_->GetElemSolutionAsMatrix( magPotMatrix,connecth);
          bilinear_stiff->SetActElemSol( magPotMatrix );
        }
        
        bilinear_stiff->CalcElementMatrix(ptCoord, elemmat);
        sol_->GetElemSolution(magvecpot,connecth);
        
        help = elemmat * magvecpot;
        energy[i] += 0.5 *( help * magvecpot);
        
      }  
      delete bilinear_stiff;    
    }

    std::string resulttype = "Electric Energy";
    std::string unit = "(Ws)";
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

    StdVector<std::string> regionNames;
    ptgrid_->RegionIdToName( regionNames, calcEnergy_ );
    Info->WriteResult(pdename_,  resulttype, regionNames, energy, unit,
                      analysis, analysisVal);
  }

  // *************
  //   ComputeUI
  // *************
  void MagPDE::ComputeUI(Vector<Double>& uiSD) {

    ENTER_FCN( "MagPDE::ComputeUI" );

    uiSD.Resize(coilDef_.GetSize());
  
    // loop over all subdomains
    for (UInt actSD=0; actSD<subdoms_.GetSize(); actSD++) {

      for (UInt dom=0; dom<coilDef_.GetSize(); dom++) {
        if (subdoms_[actSD] == coilRegionId_[dom]) {
           
          StdVector<Elem*> elemssd;             
          // get vector of Elem of subdomain with color: subdoms[isd]
          ptgrid_->GetVolElems( elemssd,subdoms_[actSD] );
            

          // loop over elements of subdomain        
          for (UInt actEl=0; actEl< elemssd.GetSize(); actEl++) {
            BaseFE * ptEl = elemssd[actEl]->ptElem;
                
            const UInt nrIntPts= ptEl->GetNumIntPoints();
            const Vector<Double> & intWeights = ptEl->GetIntWeights();  
            Double jacDet;
                
            StdVector<UInt> connect;
            connect = elemssd[actEl]->connect;

            Matrix<Double> ptCoord;
            GetElemCoords(connect, ptCoord);

            Vector<Double> magVecDeriv1Elem;
            GetDerivSolVecOfElement(magVecDeriv1Elem,connect);
            Double uiElem=0;
                
            for (UInt actIntPt=1; actIntPt<=nrIntPts;  actIntPt++) {
              Vector<Double> shapeFnc;
              jacDet = ptEl->CalcJacobianDetAtIp(actIntPt, ptCoord);    
              ptEl -> GetShFncAtIp(shapeFnc, actIntPt);

              uiElem += shapeFnc * magVecDeriv1Elem;
                    
              if (isaxi_) {
                Vector<Double> coordAtIP;
                 coordAtIP = ptCoord * shapeFnc;
                uiElem += shapeFnc * magVecDeriv1Elem * 2 * PI * coordAtIP[0]
                  * intWeights[actIntPt-1];
              }
              else {
                uiElem += shapeFnc * magVecDeriv1Elem * intWeights[actIntPt-1];
              }
            }
            
            uiSD[dom] += uiElem;
          }    
        }
      }
    }
  }


  // ***************************************************
  //   Store currents/voltages and inductivity in file
  // ***************************************************
  void MagPDE::WriteUI2File(Vector<Double>& uiSD) {

    ENTER_FCN( "MagPDE::WriteUI2File" );

    Vector<UInt> coilIDs;   // just positive ids
    coilIDs.Push_back(abs(coilDef_[0]->id_));

    Integer maxID = coilDef_[0]->id_;

    for (UInt dom=1; dom < coilDef_.GetSize(); dom++) {

      Boolean isInVec = FALSE;
      
      for (UInt dom2=0; dom2 < coilIDs.GetSize(); dom2++)    
        if (abs(coilDef_[dom]->id_) == coilIDs[dom2])
          isInVec = TRUE;
      
      if (!isInVec)
        coilIDs.Push_back(abs(coilDef_[dom]->id_));
      
      if (abs(coilDef_[dom]->id_) > maxID)
        maxID = coilDef_[dom]->id_;
    }

    Vector<Double> uiID(maxID);
    uiID.Init();

    for (UInt dom=0; dom < coilDef_.GetSize(); dom++) {
      Integer actCoilID = coilDef_[dom]->id_;
      uiID[abs(actCoilID)-1] += uiSD[dom] * actCoilID/abs(actCoilID);
    }

    *UIfile_ << solveStep_->GetActTime() << " \t";
    for (UInt actID=0; actID < coilIDs.GetSize(); actID++) {
      if ( coilDef_[coilIDs[actID]-1]->coilType_ == Coil::MEASUREMENT2D )
        *UIfile_ << uiID[coilIDs[actID]-1] *
          coilDef_[coilIDs[actID]-1]->windingCrossSection_ << " \t";
    }
  
    *UIfile_ << myEndl;
  
  }


  // ******************************************************
  //   Query parameter object for information about coils
  // ******************************************************
  void MagPDE::ReadCoils() {

    ENTER_FCN( "MagEdgePDE::ReadCoils" );

    StdVector<std::string> coilNamesAux, helper;
    StdVector<std::string> keyVec;
    StdVector<std::string> attrVec;
    StdVector<std::string> valVec;

    // Get list of coils for magnetic PDE
    params->GetCoilList( coilNamesAux, pdename_, bcSequenceTag_);
    ptgrid_->RegionNameToId( coilRegionId_, coilNamesAux );

    // Now get description of each coil and generate corresponding
    // coil object
    UInt nrCoils = coilRegionId_.GetSize();
    if ( nrCoils > 0 ) {
      Info->PrintF( pdename_, "Using the following coils:\n" );
      coilDef_.Reserve( nrCoils );
      for ( UInt k = 0; k < nrCoils; k++ ) {
        coilDef_.Push_back( new Coil( coilRegionId_[k], pdename_, ptgrid_) );
        Info->PrintCoil( *coilDef_[k], analysistype_ );
      }
    }
  }


  // ********************************************************
  //   Query parameter object for information about magnets
  // ********************************************************
  void MagPDE::ReadMagnets() {

    ENTER_FCN( "MagEdgePDE::ReadMagnets" );

    // get domain names of magnets
    //    params->GetList( "name", magnetsDomain_, pdename_, "magnets" );
    // vectors for parameter handling
    StdVector<std::string> keyVec;
    StdVector<std::string> attrVec;
    StdVector<std::string> valVec;
    StdVector<std::string> regionNames;

    keyVec = pdename_, "magnets", "magnet", "name";
    attrVec = "", "", "tag";
    valVec  = "", "", bcSequenceTag_;

    params->GetList(keyVec, attrVec, valVec, regionNames);
    ptgrid_->RegionNameToId(magnetsDomain_, regionNames);

    if ( magnetsDomain_.GetSize() > 0 ) {

      Info->PrintF( pdename_,
                    "Found permanent magnets in the following regions:\n" );

      Double tmpDir;

      // Construct vectors for restricted search parameter
      StdVector<std::string> keyVec;
      StdVector<std::string> attrVec;
      StdVector<std::string> valVec;
      attrVec = "", "", "name";

      // for each magnet ...
      for ( UInt k = 0; k < magnetsDomain_.GetSize(); k++ ) {

        // ... read direction of magnetisation
        valVec = "", "", regionNames[k];

        keyVec  = pdename_, "magnets", "magnet", "orientX";
        params->Get( keyVec, attrVec, valVec, tmpDir);
        magnetsOriX_.Push_back( tmpDir);

        keyVec  = pdename_, "magnets", "magnet", "orientY";
        params->Get( keyVec, attrVec, valVec, tmpDir );
        magnetsOriY_.Push_back( tmpDir );

        keyVec  = pdename_, "magnets", "magnet", "orientZ";
        params->Get( keyVec, attrVec, valVec, tmpDir );
        magnetsOriZ_.Push_back( tmpDir );

        // ... report name to logfile
        Info->PrintF( pdename_, " %s\n", regionNames[k].c_str());
      }

      Info->PrintF( "", "\n" );
    }
  }


  // ***********************************************************************
  //   Obtain information on desired output quantities from parameter file
  // ***********************************************************************
  void MagPDE::ReadStoreResults() {

    ENTER_FCN( "MagPDE::ReadStoreResults" );

    // Construct vectors for restricted parameter search
    StdVector<std::string> keyVec;
    StdVector<std::string> attrVec;
    StdVector<std::string> valVec;
    StdVector<std::string> regionNames;
    StdVector<RegionIdType> regionIds;
    std::string quantity;

    // *****************************
    // Determine nodal results
    // ***************************** 
    StdVector<std::string> nodeValues;
    Enum2String(MAG_POTENTIAL, quantity);
    keyVec  = pdename_, "storeResults", "nodeResults", "region";
    attrVec = "", "", "type";
    valVec = "", "", quantity;
    params->GetList( keyVec, attrVec, valVec, nodeValues);
    if (nodeValues.GetSize() > 0) {
      saveSol_ = TRUE;
      hasOutput_ = TRUE;
    }

    Enum2String(MAG_FORCE_VWP, quantity);
    keyVec  = pdename_, "storeResults", "nodeResults", "nodes";
    attrVec = "", "", "type";
    valVec = "", "",quantity;
    params->GetList( keyVec, attrVec, valVec, calcForceVWP_);

    if ( calcForceVWP_.GetSize() > 0 ) {

      UInt numNodes = 0;
      
      // initalize force nodestoresol
       // intialize corresponding storesolution object
      Force_.SetNumSolutions(1);
      Force_.SetNumNodes(numPDENodes_);
      Force_.SetSolutionType(MAG_FORCE_VWP);
      Force_.SetNumDofs(dim_);
      Force_.SetPtrEQNData(eqnData_, ptgrid_); 
      Force_.Init();


      // count complete number of nodes
      for ( UInt i=0; i<calcForceVWP_.GetSize(); i++ ) {
        numNodes+= ptgrid_->GetNumNodes(calcForceVWP_[i]);

        ForceNodes_.Resize(numNodes);
      
        UInt inode =0;
        StdVector<UInt> nodesConverted;
        std::list<UInt>::iterator it;

        // get for each nodeslist all nodes
        for (UInt i=0; i<calcForceVWP_.GetSize(); i++)
          {
            ptgrid_->GetNodesByName( nodesConverted, calcForceVWP_[i]);
          
            for ( UInt j=0; j<nodesConverted.GetSize(); j++, inode++)
              ForceNodes_[inode] = nodesConverted[j];
          }

      }
      //      std::cout << "ForceNodes:\n" << ForceNodes_ << std::endl;

      //initialize the force operator
      NodeStoreSol<Double> * solhelp = dynamic_cast<NodeStoreSol<Double> *>(sol_);
      ForceOpVWP_ = new  MagForceOp(ptgrid_, this, eqnData_, *solhelp, dim_, 
                                    materials_,  isaxi_);
      
      keyVec  = pdename_, "storeResults", "nodeResults", "region";
      attrVec = "", "", "type";
      valVec = "", "", quantity;
      params->GetList( keyVec, attrVec, valVec, regionNames);
      ptgrid_->RegionNameToId( regionIds, regionNames ); 
      
      ForceOpVWP_->Setup( regionIds, ForceNodes_);
    }
   

    // *****************************
    // Determine element results
    // *****************************
    StdVector<std::string> elemResults;
    keyVec  = pdename_, "storeResults", "elemResults", "region";
    attrVec = "", "", "type";  
    
    // --- Magnetic Flux Density ---
    Enum2String(MAG_FLUX_DENSITY, quantity);
    valVec  = "", "", quantity;
    params->GetList( keyVec, attrVec, valVec, regionNames );
    ptgrid_->RegionNameToId( calcBfield_, regionNames );
    
    // If the symbolic name is "all" compute magnetic Fieldfor all regions
    if ( calcBfield_.GetSize() == 1 && calcBfield_[0] == ALL_REGIONS ) {
      calcBfield_ = subdoms_;
    }
    
    // Log to info file
    if ( calcBfield_.GetSize() > 0 ) {
      hasOutput_ = TRUE;
      Info->PrintF( pdename_,
                    "Computing magFluxDensity for regions:\n" );
      for ( UInt k = 0; k < regionNames.GetSize(); k++ ) {
        Info->PrintF( pdename_, " %s\n", regionNames[k].c_str() );
      }
      Info->PrintF( "", "\n" );
    }
    
    // --- Magnetic Energy ---
    Enum2String(MAG_ENERGY, quantity);
    valVec  = "", "", quantity;
    params->GetList( keyVec, attrVec, valVec, regionNames );
    ptgrid_->RegionNameToId( calcEnergy_, regionNames );
    
    // If the symbolic name is "all" compute magnetic Fieldfor all regions
    if ( calcEnergy_.GetSize() == 1 && calcEnergy_[0] == ALL_REGIONS ) {
      calcEnergy_ = subdoms_;
    }
    
    // Log to info file
    if ( calcEnergy_.GetSize() > 0 ) {
      hasOutput_ = TRUE;    
      Info->PrintF( pdename_,
                    "Computing magEnergy for regions:\n" );
      for ( UInt k = 0; k < regionNames.GetSize(); k++ ) {
        Info->PrintF( pdename_, " %s\n", regionNames[k].c_str() );
      }
      Info->PrintF( "", "\n" );
    }

    // --- Magnetic Eddy Current ---
    Enum2String(MAG_EDDY_CURRENT, quantity);
    valVec  = "", "", quantity;
    params->GetList( keyVec, attrVec, valVec, regionNames );
    ptgrid_->RegionNameToId( calcEddy_, regionNames );
    
    // If the symbolic name is "all" compute magnetic Fieldfor all regions
    if ( calcEddy_.GetSize() == 1 && calcEddy_[0] == ALL_REGIONS ) {
      calcEddy_ = subdoms_;
    }
    
    // Log to info file
    if ( calcEddy_.GetSize() > 0 ) {
      hasOutput_ =TRUE;
      Info->PrintF( pdename_,
                    "Computing magEddyCurrent for regions:\n" );
      for ( UInt k = 0; k < regionNames.GetSize(); k++ ) {
        Info->PrintF( pdename_, " %s\n", regionNames[k].c_str() );
      }
      Info->PrintF( "", "\n" );
    }
   
    // *****************************
    // Determine nodal history
    // *****************************
    StdVector<std::string> saveNodeHist;
    keyVec  = pdename_, "storeResults", "nodeHistory", "saveNodes";
    attrVec = "", "", "type";
    
    // --- Magnetic Potential ---
    Enum2String(MAG_POTENTIAL, quantity);
    valVec  = "", "", quantity;
    params->GetList( keyVec, attrVec, valVec, saveNodeHist );
  
    if (saveNodeHist.GetSize() > 0) {
      saveSolHist_ = TRUE;
      hasOutput_ = TRUE;
      Info->PrintF( pdename_, "Saving magPotential for Nodes:\n" );
      for ( UInt k = 0; k < saveNodeHist.GetSize(); k++ ) {
        Info->PrintF( pdename_, " %s\n", saveNodeHist[k].c_str() );
      }
      Info->PrintF( "", "\n" );
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


  // ======================================================
  // COUPLING SECTION
  // ======================================================


  void MagPDE::InitCoupling(PDECoupling * Coupling) {

    ENTER_FCN( "MagPDE::InitCoupling" );
  
    isIterCoupled_ = TRUE;
    ptCoupling_   = Coupling;

    // Enable update of geometry
    const UInt numCouplings = ptCoupling_->GetNumOutputCouplings();  

    StdVector<StdVector<UInt> > elemNodeToCouplingNode_tmp;
    elemNodeToCouplingNode_.Resize(numCouplings);


    for ( UInt actCoupling = 0; actCoupling < numCouplings; actCoupling++ ){

      if (ptCoupling_->GetOutputQuantity(actCoupling) == MAG_FORCE_LORENTZ) {

        // Intialize the memory of the coupling values
        ptCoupling_->CreateCouplingVector(actCoupling,isComplex_);

        //get the element-node to coupling node matching
        StdVector<std::string> couplRegions;
        StdVector<RegionIdType> regionIds;
        ptCoupling_->GetOutputRegions(actCoupling, couplRegions);
        ptgrid_->RegionNameToId( regionIds, couplRegions );

        //Get total number of coupling elements
        UInt totalCouplingElems = ptgrid_->GetNumElems( regionIds );
        
        elemNodeToCouplingNode_tmp.Clear();
        elemNodeToCouplingNode_tmp.Resize(totalCouplingElems);

        UInt offset = 0;
        for ( UInt reg = 0; reg < couplRegions.GetSize(); reg++ ) {

          // find subdomain index
          Integer SDidx=-1; for (UInt sd=0; sd<subdoms_.GetSize(); sd++) {
            if (regionIds[reg] == subdoms_[sd]) {
              SDidx = sd;
              break;
            }
          }
              
          if (SDidx==-1) {
            (*error) << "magneticPDE: Coupling Region is not within the "
                     << "subdomains of the PDE";
            Error( __FILE__, __LINE__ );
          }

          StdVector<Elem*> elemssd;
          ptgrid_->GetVolElems(elemssd, subdoms_[SDidx]);

          StdVector<UInt> * couplingnodes = NULL;
          ptCoupling_->GetOutputNodes(actCoupling, couplingnodes);
          if ( couplingnodes == NULL ) {
            Error( "Pointer couplingnodes = NULL!", __FILE__, __LINE__ );
          }

          for (UInt actEl=0; actEl< elemssd.GetSize(); actEl++) {
            StdVector<UInt> & connecth = elemssd[actEl]->connect;
            elemNodeToCouplingNode_tmp[offset+actEl].Resize(connecth.GetSize());

            for ( UInt ielemnode = 0; ielemnode < connecth.GetSize();
                  ielemnode++ ) {
              for ( UInt cnode = 0; cnode < (*couplingnodes).GetSize();
                    cnode++ ) {
                if (connecth[ielemnode] == (*couplingnodes)[cnode] ) {
                  elemNodeToCouplingNode_tmp[offset+actEl][ielemnode] = cnode;
                  break;
                }
              }
            }
          }

          //in the case, that we have more than one coupling region!
          offset = elemssd.GetSize();
        }

        elemNodeToCouplingNode_[actCoupling]  = elemNodeToCouplingNode_tmp;
      }
    }
  }

  void MagPDE::CalcOutputCoupling() {

    ENTER_FCN( "MagPDE::CalcOutputCoupling" );

    SolutionType quantity;
    StdVector<UInt> * couplingNodes = NULL;
    CFSVector * values = NULL;
    UInt forcesCount = 0;

    // at first, check if this PDE is iterative coupled
    if (isIterCoupled_ == FALSE)
      return;

    // loop over all output coupling quantities
    for ( UInt actCoupling = 0;
          actCoupling < ptCoupling_->GetNumOutputCouplings();
          actCoupling++ ) {

      quantity = ptCoupling_->GetOutputQuantity(actCoupling);
      ptCoupling_->GetOutputValues(actCoupling, values);

      Vector<Double> *temp = dynamic_cast<Vector<Double> *>(values);
      
      switch(ptCoupling_->GetOutputType(actCoupling)) {
          
      case NODE:
        ptCoupling_->GetOutputNodes(actCoupling, couplingNodes);
          
        if (quantity == MAG_FORCE_LORENTZ) {
          CalcNodeForceLorentz(*temp, elemNodeToCouplingNode_[forcesCount],
                               actCoupling, couplingNodes->GetSize());
          forcesCount++;
        }

        break;
          
      case ELEM:
        Error( "No Element coupling output", __FILE__, __LINE__ );
      }
    }
  }

  void MagPDE::
  CalcNodeForceLorentz(Vector<Double> & force, 
                       StdVector<StdVector<UInt> > & elemNodeToCouplingNode,
                       UInt actCoupling, UInt numCouplingNodes) {

    ENTER_FCN( "MagPDE::CalcNodeForceLorentz" );

    NodeStoreSol<Double> *solhelp = dynamic_cast<NodeStoreSol<Double> *>(sol_);

    //get the coupling regions
    StdVector<std::string> couplRegions;
    StdVector<RegionIdType> regionIds;
    ptCoupling_->GetOutputRegions(actCoupling, couplRegions);
    ptgrid_->RegionNameToId( regionIds, couplRegions );
    MagLorentzForceOp *ForceOp; 

    ForceOp  = new MagLorentzForceOp( ptgrid_, this, eqnData_, *solhelp, isaxi_ );

    force.Init(0.0);

    Vector<Double> Jeddy;

    UInt offset = 0;
    for (UInt reg=0; reg<couplRegions.GetSize(); reg++) {

      //find subdomain index
      Integer SDidx=-1; for (UInt sd=0; sd<subdoms_.GetSize(); sd++) {
        if (regionIds[reg] == subdoms_[sd]) {
          SDidx = sd;
          break;
        }
      }

      Double conductivity;
      materials_[subdoms_[SDidx]]->GetScalar(conductivity,MAG_CONDUCTIVITY,REAL);      
      
      StdVector<Elem*> elemssd;
      ptgrid_->GetVolElems(elemssd, subdoms_[SDidx]);
  
      for (UInt actEl=0; actEl< elemssd.GetSize(); actEl++) {
        StdVector<UInt> connecth = elemssd[actEl]->connect;
        Matrix<Double> elemForce;
        GetDerivSolVecOfElement(Jeddy,connecth);
        Jeddy *= -conductivity;

        ForceOp->CalcElemMagLorentzForce(elemForce, Jeddy, elemssd[actEl]);

        // Add the element force to the according coupling node
        for (UInt ielemnode=0; ielemnode<connecth.GetSize(); ielemnode++) {
          for( UInt idim=0; idim<dim_; idim++) {
            force[elemNodeToCouplingNode[actEl+offset][ielemnode]*dim_+idim]
              += elemForce[ielemnode][idim];
          }
        }
      }

      //in the case, that we have more than one coupling region!
      offset = elemssd.GetSize();
    }
  }

  Boolean MagPDE::HasOutput( SolutionType output ) {
    ENTER_FCN( "MagPDE::HasOutput" );
    if (output == MAG_FORCE_LORENTZ) {
      return TRUE;
    }
    return FALSE;
  }

} // end of namespace

