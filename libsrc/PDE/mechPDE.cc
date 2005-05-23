#include "mechPDE.hh"

#include "Forms/forms_header.hh"
#include "Forms/linElastInt.hh"
#include "Forms/massInt.hh"
#include "DataInOut/writeresults.hh"
#include "Driver/assemble.hh"
#include "newmark.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "Driver/solveStepMech.hh" 
#include "CoupledPDE/pdecoupling.hh"

namespace CoupledField
{

  MechPDE::MechPDE(Grid * aptgrid, TimeFunc *aptTimeFunc, WriteResults *aptOut)
    :SinglePDE(aptgrid, aptOut, aptTimeFunc), 
     lambdaMat(NULL),
     mueMat(NULL),
     preStressVal_(0.0)

  {
    ENTER_FCN( "MechPDE::MechPDE" );

    pdename_          = "mechanic";
    pdematerialclass_ = "piezo";
   


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
      dofspernode_ = 3;
      Info->PrintF("", "=== 3D PROBLEM\n");
    }
    else if ( subType_ == "axi" && probGeo == "axi" ) {
      isaxi_ = TRUE;
      dofspernode_ = 2;
      Info->PrintF("", "=== AXISYSMMETRIC PROBLEM\n");
    }
    else if ( subType_ == "planeStrain" && probGeo == "plane" ) {
      dofspernode_ = 2;
      Info->PrintF("", "=== PLANE STRAIN PROBLEM\n");
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
    solTypes_ = MECH_DISPLACEMENT;
    solDofs_ = dofspernode_;
    
    effectiveMass_ = params->IsSet( "effMass" );
    
    // *********************************
    //  Check damping model
    // *********************************
    
    if( params->HasValue( "type", "rayleigh", pdename_, "damping" ) )
      {
        dampingType_ = RAYLEIGH;
        Info->PrintF(pdename_, " Using RAYLEIGH damping\n" );
      }
    else
      {
        dampingType_ = NONE;
      }

    if (dampingType_)
      needsDampingMatrix_ = TRUE;
      
    // *********************************
    //  Check for pressure loads
    // *********************************

    //check for pressure loads
    StdVector<std::string> regionNames;

    params->GetList( "name"    , regionNames , pdename_, "pressure" );
    params->GetList( "value"   , pressVals_ , pdename_, "pressure" );
    params->GetList( "dynamics", pressFnc_  , pdename_, "pressure" );
    
    ptgrid_->RegionNameToId( pressSurf_, regionNames );
    // Check consistency
    if ( pressSurf_.GetSize() != pressVals_.GetSize() )
      {
        std::string errmsg = "PressureLoads: ";
        errmsg += "#name = " + Info->GenStr(pressSurf_.GetSize());
        errmsg += ", #value = " + Info->GenStr(pressVals_.GetSize());
        errmsg += ", #dynamics = " + pressFnc_.GetSize() + '\n';
        Info->Error( errmsg, __FILE__, __LINE__ );
      }

    if (pressSurf_.GetSize() > 0)
      surfdoms_ = pressSurf_;
    // We need not have as many function/filenames as pressureloads!
    for ( Integer k = pressFnc_.GetSize(); k < pressSurf_.GetSize(); k++ )
      {
        pressFnc_.Push_back( "none" );
      }

    //check for prestressing
    ReadPreStressing();

  }


  MechPDE::~MechPDE()
  {
    ENTER_FCN( "MechPDE::~MechPDE" );

    if  (lambdaMat)
      delete lambdaMat;
    
    if  (mueMat)
      delete mueMat;
  }


  
  void MechPDE::InitNonLin()
  {
    ENTER_FCN( "MechPDE::InitNonLin");

    // ==============================================================
    // NOTE: Currently we can only treat geometric non-linearity and
    //       we assume that for a mechanic PDE all regions either
    //       are linear or non-linear!
    // ==============================================================
    StdVector<std::string> nonLinRegion;
    params->GetList( "nonLinear", nonLinRegion, pdename_, "region" );
    // Should not happen with validating parser, but beware!
    if ( nonLinRegion.GetSize() == 0 ) {
      nonLin_ = FALSE;
    }
    else {
      for ( Integer k = 1; k < nonLinRegion.GetSize(); k++ ) {
        if ( nonLinRegion[k] != nonLinRegion[0] ) {
          Info->Error( "Non-linearity should be the same for all regions!",
                       __FILE__, __LINE__ );
        }
      }
      nonLin_ = nonLinRegion[0] == "geo" ? TRUE : FALSE;
    }

    // In non-linear case determine type of line search strategy
    if ( nonLin_ == TRUE ) {

      Info->PrintF( pdename_, "Non-linearity in %d regions\n",
                    nonLinRegion.GetSize() );

      // type of line search
      params->Get( "type", lineSearch_, pdename_, "lineSearch" );

      if ( lineSearch_ == "no" ) {
        Info->PrintF( pdename_, "Performing no line search\n" );
      }
      else {
        Info->PrintF( pdename_, "Will perform line search\n" );
      }

    }

    // If no non-linearity we do not perform line search anyhow
    else {
      lineSearch_ = "no";
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

    // ------------------------------------------
    //   Get information on reduced integration
    // ------------------------------------------
    params->GetList( "reducedInt", reducedIntegration_, pdename_, "region" );

    if ( nonLin_ == TRUE ) {
      for ( Integer i = 0; i < reducedIntegration_.GetSize(); i++ ) {
        if ( reducedIntegration_[i] == "yes" ) {
          (*error) << "Currently we do not support non-linearity with "
                   << "reduced integration!";
          Error( __FILE__, __LINE__ );
        }
      }
    }
  }
  

  void MechPDE::DefineIntegrators()
  {
    ENTER_FCN( "MechPDE::DefineIntegerators" );

    //voulme integrators
    for (int actSD = 0; actSD < subdoms_.GetSize(); actSD++)
      {

        // ==============  add stiffness ======================================

        MaterialData actSDMat(materialData_[actSD]);

        if (reducedIntegration_[actSD] == "yes")  
          {
            // ==================  add reduced stiffness ==========================

            std::cout << "Do reduced int" << std::endl;
          
            lambdaMat = new MaterialData(actSDMat);
            mueMat = new MaterialData(actSDMat);
          
            // use a "different" material data set for reduced integration
            CalcReducedMat(*lambdaMat, *mueMat, actSDMat);      

            BaseForm * bilinearStiff1 = GetStiffIntegrator(*mueMat);
            IntegratorDescriptor * actIntDescr1 =
              new IntegratorDescriptor(bilinearStiff1, STIFFNESS);

            //check for damping
            if (dampingType_ == RAYLEIGH) {
              actIntDescr1->SetSecondaryMat(DAMPING, actSDMat.GetDampingBeta(),analysistype_);
            }
            assemble_->AddIntegrator(actIntDescr1, subdoms_[actSD]);


            BaseForm * bilinearStiff2 = GetStiffIntegrator(*lambdaMat);
            IntegratorDescriptor * actIntDescr2 =
              new IntegratorDescriptor(bilinearStiff2, STIFFNESS);

            //for this integrator, we need reduced integration
            //see Hughes, pp.219
            actIntDescr2->SetReducedInt();

            //check for damping
            if (dampingType_ == RAYLEIGH) {
              actIntDescr2->SetSecondaryMat(DAMPING, actSDMat.GetDampingBeta(),analysistype_);
            }
            assemble_->AddIntegrator(actIntDescr2, subdoms_[actSD]);

          }

        else 
          {
            // ==============  add "standard" stiffness ===========================
            BaseForm * bilinearStiff = GetStiffIntegrator(actSDMat);
            IntegratorDescriptor * actIntDescr =
              new IntegratorDescriptor(bilinearStiff, STIFFNESS);

            //check for damping
            if (dampingType_ == RAYLEIGH) {
              actIntDescr->SetSecondaryMat(DAMPING, actSDMat.GetDampingBeta(),analysistype_);
            }
            assemble_->AddIntegrator(actIntDescr, subdoms_[actSD]);


            //for prestressing
            for ( Integer preStr=0; preStr<preStressDomain_.GetSize(); preStr++ ) {
              if ( subdoms_[actSD] == preStressDomain_[preStr]) {
                Vector<Double> preStrVal(3);
                preStrVal[0] = preStressValX_[preStr];
                preStrVal[1] = preStressValY_[preStr];
                preStrVal[2] = preStressValZ_[preStr];

                BaseForm * bilinearPreStress;
                if (subType_ == "planeStrain")
                  bilinearPreStress = new PreStressIntPlaneStrain(actSDMat, preStrVal);
                else if (subType_ == "axi")
                  bilinearPreStress = new PreStressIntAxi(actSDMat, preStrVal);
                else if (subType_ == "3d")
                  bilinearPreStress = new PreStressInt3D(actSDMat, preStrVal);
                else 
                  Info->Error("Unknown subtype in mech PDE! ",__FILE__,__LINE__);               

                IntegratorDescriptor * actIntDescrPre =
                  new IntegratorDescriptor(bilinearPreStress, STIFFNESS);
                assemble_->AddIntegrator(actIntDescrPre, subdoms_[actSD]);
              }
            }

          }


        // ==============  add nonlinear stiffness ============================
        if (nonLin_)
          {
            BaseForm *nLinPart1 = NULL;
            BaseForm *nLinPart2 = NULL;
          
            if (subType_ == "3d")
              {   
                nLinPart1 = new nLinMech3dInt_BNonLin(actSDMat);    
                nLinPart2 = new nLinMech3dInt_PiolaStress(actSDMat);
              }

            else if (subType_ == "planeStrain")
              {
                nLinPart1 = new nLinMechPlaneStrainInt_BNonLin(actSDMat);    
                nLinPart2 = new nLinMechPlaneStrainInt_PiolaStress(actSDMat);

              }
            else if (subType_ == "axi")
              {
                nLinPart1 = new nLinMechAxiInt_BNonLin(actSDMat);    
                nLinPart2 = new nLinMechAxiInt_PiolaStress(actSDMat);

              }
          
            assemble_->AddIntegrator(nLinPart1, subdoms_[actSD], STIFFNESS,
                                     nonLin_);
            assemble_->AddIntegrator(nLinPart2, subdoms_[actSD], STIFFNESS,
                                     nonLin_);

            // assemble prestress, if in config-file given
            //    if (preStressVal_)
            //      AssemblePreStressMat(ptEl, connect_PDE, ptCoord,
            //      actMatData, elDisp);
          }


        // ==============  add mass ===========================================
        double density = actSDMat.GetDensity();
        BaseForm * bilinearMass  = new MassInt(density, dofspernode_, isaxi_);

        IntegratorDescriptor * actIntDescr =
          new IntegratorDescriptor(bilinearMass, MASS);

        //check for damping (mass part)
        if (dampingType_ == RAYLEIGH) {
          actIntDescr->SetSecondaryMat(DAMPING, actSDMat.GetDampingAlfa(),analysistype_);
        }

        assemble_->AddIntegrator(actIntDescr, subdoms_[actSD]);


        // ==================== RHS ===========================================
        if (nonLin_)
          {
            BaseForm * rhsSource = new nLinMech_linFormInt(actSDMat, isaxi_);
            assemble_->AddRhsIntegrator(rhsSource, subdoms_[actSD], nonLin_);
          }      
      }

    //surface integrators
    //RHS-part
    Integer nonlin = 0;
    for (Integer actSF = 0; actSF < pressSurf_.GetSize(); actSF++) {
      BaseForm * rhsSrcSurf = new PressureLinForm(pressVals_[actSF], isaxi_);
      assemble_->AddRhsSrcSurfIntegrator(rhsSrcSurf, pressSurf_[actSF], pressFnc_[actSF],
                                         nonlin);
    }
    
  }


  BaseForm *
  MechPDE::GetStiffIntegrator(MaterialData& actSDMat, Boolean reducedInt)
  {

    ENTER_FCN( "MechPDE::GetStiffIntegrator" );
  
    BaseForm * bilinearStiff = NULL;

    if (subType_ == "planeStrain")
      bilinearStiff = new mechPlainStrainInt(actSDMat);
    else if (subType_ == "axi")
      bilinearStiff = new mechAxiInt(actSDMat);
    else if (subType_ == "3d")
      bilinearStiff = new mech3DInt(actSDMat);
    else 
      Info->Error("Unknown subtype in mech PDE! ",__FILE__,__LINE__);

    return bilinearStiff;
  }

  void MechPDE::DefineSolveStep()
  {
    ENTER_FCN( "MechPDE::DefineSolveStep" );

    solveStep_ = new SolveStepMech(*this);
  }


  void MechPDE::
  GetDirection(Directions& dir, const std::string keyword)
  {
    ENTER_FCN( "MechPDE::GetDirection" );

    std::string dirString;

    params->Get( keyword, dirString, pdename_);

    if (dirString == "x")
      dir = X;
    else if (dirString == "y")
      dir = Y;
    else if (dirString == "z")
      dir = Z;
    else if (dirString == "radXY")
      dir = radXY;
    else if (dirString == "radXZ")
      dir = radXZ;
    else if (dirString == "radYZ")
      dir = radYZ;
    else
      {
        // According to the Schema definition of the parameter file this cannot
        // happen. Did the parser not perform validation?
        std::string errmsg = "Direction should be one of x, y, z, radXY, ";
        errmsg += "radXZ, radYZ\nand not" + dirString;
        Info->Error( errmsg, __FILE__, __LINE__ );
      }
  }




  void MechPDE::CalcReducedMat(MaterialData& lambdaMat, MaterialData& mueMat, MaterialData& mat)
  {
    ENTER_FCN( "MechPDE::CalcReducedMat" );

    Double lambda, mue;
    mat.GetPiezoMatrixData(1,2, lambda);
    mat.GetPiezoMatrixData(4,4, mue);
  
    Matrix<Double> * lMechMat = lambdaMat.GetMatrix();
    lMechMat -> Init();

    Matrix<Double> * mueMechMat = mueMat.GetMatrix();
    mueMechMat -> Init();


    for(Integer actRow=0; actRow<3; actRow++)
      {
        for(Integer actCol=0; actCol<3; actCol++)
          (*lMechMat)[actRow][actCol] = lambda;

        (*mueMechMat)[actRow][actRow] = 2*mue;
      }

    for(Integer actRow=3; actRow<6; actRow++)
      (*mueMechMat)[actRow][actRow] = mue;  

    std::cout << "LambadMat:\n" << *lMechMat << std::endl;
    std::cout << "MuMat:\n" << *mueMechMat << std::endl;

  }



  Integer MechPDE::GetNrBCDof(const std::string & dofStartString)
  {
    ENTER_FCN( "MechPDE::GetNrBCDof" );

    Integer nrActDof;
    
    if (dofStartString == "ux")
      nrActDof = 1;
    else 
      if (dofStartString == "uy")
        nrActDof = 2;
      else
        if (dofspernode_ == 3)
          if (dofStartString == "uz")
            nrActDof = 3;
          else
            Error("Unknown dof-type in homog. BC; substring must start with ux, uy or uz!!",__FILE__,__LINE__);
        else
          Error("Unknown dof-type in homog. BC; substring must start with ux or uy!!",__FILE__,__LINE__);
    
    return nrActDof;
  }




  // ======================================================
  // COUPLING SECTION
  // ======================================================


  void MechPDE::InitCoupling(PDECoupling * Coupling)
  {
    ENTER_FCN( "MechPDE::InitCoupling" );
  
    isIterCoupled_ = TRUE;
    ptCoupling_   = Coupling;
  
    for (Integer i=0; i<ptCoupling_->GetNumOutputCouplings(); i++)
      {
        if (ptCoupling_->GetOutputQuantity(i) == MECH_DISPLACEMENT)
          {
            // Intialize the memory of the coupling values
            ptCoupling_->CreateCouplingVector(i,isComplex_);
          }

        if (ptCoupling_->GetOutputQuantity(i) == MECH_FORCE)
          {
            // Intialize the memory of the coupling values
            ptCoupling_->CreateCouplingVector(i,isComplex_); 

            //now since we need a incremental formulation, initialize some necessary vectors
            isIncrFormulation_ = TRUE;
            solIncr_.Resize(eqnData_->GetNumEQNs() * eqnData_->GetNumDofsPerEQN());
            actSol_.Resize(eqnData_->GetNumEQNs() * eqnData_->GetNumDofsPerEQN());
            solIncr_.Init(0);
            actSol_.Init(0);
          }
      }

  }





  void MechPDE::CalcOutputCoupling()
  {
    ENTER_FCN( "MechPDE::CalcOutputCoupling" );

    Integer dof = 0;
    SolutionType quantity;
    StdVector<Integer> * couplingnodes = NULL;
    StdVector<Elem*> * couplingElems = NULL;
    CFSVector * temp_values = NULL;
    Vector<Double> * values;
    StdVector<MaterialData*> * materials = NULL;
    StdVector<std::string> outputRegions;
  

    // loop over all output coupling quantities
    for (Integer i=0; i<ptCoupling_->GetNumOutputCouplings(); i++)
      {
        quantity = ptCoupling_->GetOutputQuantity(i);
        ptCoupling_->GetOutputValues(i, temp_values);

        values = dynamic_cast<Vector<Double>*>(temp_values);
        
        switch(ptCoupling_->GetOutputType(i))
          {
          case NODE:
          
            if (quantity == MECH_DISPLACEMENT)
              {
                ptCoupling_->GetOutputNodes(i, couplingnodes);
                      
                sol_->NodeSolutionToCoupling(*values, *couplingnodes);
              }
          

            if (quantity == MECH_FORCE)
              {
                ptCoupling_->GetOutputNodes(i, couplingnodes);
                ptCoupling_->GetOutputElements(i, couplingElems);
                ptCoupling_->GetOppositeMaterials(i, materials);
                dof = ptCoupling_->GetOutputDof(i);

                CalcAcousticCouplingRHS(couplingElems, 
                                        *materials,
                                        *couplingnodes, 
                                        *values, dof);           
              } 
            break;

          case ELEM:
            Error("No Element coupling output", __FILE__,__LINE__);
          }
      }
  }


  void MechPDE::CalcAcousticCouplingRHS( StdVector<Elem*> * couplingElems, 
                                         StdVector<MaterialData*> & materials,
                                         StdVector<Integer>& couplingNodes,
                                         Vector<Double> & elemCouplingSols,
                                         Integer couplingdof )
  {
    ENTER_FCN( "MechPDE::CalcAcousticCouplingRHS" );

    Matrix<Double> ptCoord, elemMat;
    Vector<Double> normal, sol, forceOnElem, nSol;
    Elem * ptVolElem;
    Double sign = 0.0;
    Double density = 0.0;
    Integer matIndex = -1;

    elemCouplingSols.Init(0.0);
  
    for (Integer actElem=0; actElem<couplingElems->GetSize(); actElem++)
      {
        // Perform cast from volume element to surface element, since
        // mech-acou coupling makes only sense on surface elements
        SurfElem * actCoupleElem = 
          dynamic_cast<SurfElem*> ((*couplingElems)[actElem]);
        

        BaseFE * ptElem = actCoupleElem->ptElem;
        StdVector<Integer> & connecth = (*couplingElems)[actElem]->connect;
        GetElemCoords(connecth, ptCoord);
      
        // Try to find according region for first neighbouring volume
        // element of the surface element
        matIndex = subdoms_.Find(actCoupleElem->ptVolElem1->regionId);
        
        // If first volume element does not belong to acoustic PDE, try the
        // second one
        if ( matIndex == -1 ) {
          matIndex = subdoms_.Find(actCoupleElem->ptVolElem2->regionId);
          ptVolElem = actCoupleElem->ptVolElem2;
          sign = -1.0 * actCoupleElem->normalSign;
        } else {
          ptVolElem = actCoupleElem->ptVolElem1;
          sign = 1.0 * actCoupleElem->normalSign;
        }
        
        if ( matIndex == -1) {
          (*error) << "MechPDE::CalcAcousticCouplingRHS: The two volume "
                   << "element neighbours of surface element Nr. "
                   << actCoupleElem->elemNum << " do not belong to my regions!";
          Error( __FILE__, __LINE__ );
        }
      
        // Assign correct density
        density = materials[actElem]->GetDensity();
        
        // get correct density belonging to the the neighbouring element
        // in the fluid subdomain
        //density = (*couplingMaterials)[actElem]->GetDensity();
      
        BaseForm * bilinear_mass = new MassInt(ptElem, density, isaxi_);
        Matrix<Double> elemmat;
        bilinear_mass->CalcElementMatrix(ptCoord, elemmat);
        delete bilinear_mass;       

        
        GetDerivSolVecOfElement(sol, connecth);
        nSol.Resize(connecth.GetSize());   // solution in normal direction

        // the normal vector points outwards of the mechanical domain
        // (see. Kaltenbacher, "Num. Sim. of Mech. Act. & Sens." chapter 8.2)
        ptgrid_->CalcSurfNormal( normal, *actCoupleElem );
        normal *= sign;
      

        for (Integer actNode=0; actNode < connecth.GetSize(); actNode++)
          for (Integer actDof=0; actDof<dofspernode_; actDof++)
            nSol[actNode] += sol[actDof + actNode*dofspernode_] * normal[actDof];


        Vector<Double> forceOnElem = elemmat * nSol;  
      
        for (Integer actNode=0; actNode<ptCoord.GetSizeRow(); actNode++)
          {
            Integer nodePos = 0;
          
            while(connecth[actNode] != couplingNodes[nodePos] && nodePos < couplingNodes.GetSize()) 
              nodePos++;
            elemCouplingSols[nodePos] += forceOnElem[actNode];
          }      
      }
  } 



  Boolean MechPDE::HasOutput(SolutionType output)
  {
    ENTER_FCN( "MechPDE::HasOutput" );

    if (output == MECH_DISPLACEMENT || output == MECH_FORCE)
      return TRUE;

    return FALSE;
  }



  // ======================================================
  // TIME STEPPING SECTION
  // ======================================================


  void MechPDE :: InitTimeStepping()
  {
    ENTER_FCN( "MechPDE::InitTimeStepping" );

    if (effectiveMass_)  
      TS_alg_ = new NewmarkEffMass(pdename_, algsys_, eqnData_, needsDampingMatrix_);
    else
      TS_alg_ = new Newmark(pdename_, algsys_, eqnData_, needsDampingMatrix_);

  }



  // ======================================================
  // POSTPROCESSING SECTION
  // ======================================================


  void MechPDE::WriteResultsInFile(const Integer kstep,
                                   const Double asteptime, 
                                   Integer stepOffset,
                                   Double timeOffset)
  {
    ENTER_FCN( "MechPDE::WriteResultsInFile" );

    NodeStoreSol<Double> sol_der1Array, sol_der2Array;
    NodeStoreSol<Double> * solTransient;
    NodeStoreSol<Complex> * solHarmonic;

    lasttimecalc_=asteptime;
    laststepcalc_=kstep;
    Double actTime = lasttimecalc_ + timeOffset;
    Integer actStep = laststepcalc_ + stepOffset;
 
    if (analysistype_ == STATIC ||
        analysistype_ == TRANSIENT) {
      solTransient = dynamic_cast<NodeStoreSol<Double>*>(sol_);
    
      if (saveSol_ == TRUE ) 
        outFile_->WriteNodeSolutionTransient(*solTransient, actStep, actTime);
    
      if (saveSolHist_ == TRUE)
        outFile_->WriteNodeHistoryTransient(*solTransient, actStep, actTime);
    
      if (analysistype_== TRANSIENT) {
        if (saveDeriv1_ == TRUE)
          {
            solDeriv1_.SetAlgSysVector(getS1());
            //        outFile_->WriteNodeSolutionTransient(solDeriv1_, actStep, actTime);
          
            if (saveDeriv1Hist_ == TRUE)
              outFile_->WriteNodeHistoryTransient(solDeriv1_, actStep, actTime);
          }
      
        if (saveDeriv2_ == TRUE)
          {
            solDeriv2_.SetAlgSysVector(getS2());
            outFile_->WriteNodeSolutionTransient(solDeriv2_, actStep, actTime);
          }
        if (saveDeriv2Hist_ == TRUE)
          outFile_->WriteNodeHistoryTransient(solDeriv2_, actStep, actTime);
      }
    
      //element results
      if (calcStress_.GetSize() !=0 ) {
        outFile_->WriteElemSolutionTransient(Stress_, actStep, actTime);
      }
    }
    else if (analysistype_ == HARMONIC) {
      solHarmonic = dynamic_cast<NodeStoreSol<Complex>*>(sol_);

      if (saveSol_ == TRUE )
        outFile_->WriteNodeSolutionHarmonic(*solHarmonic,  actFreqStep_, 
                                            actFrequency_, complexFormat_);
      if (saveSolHist_ == TRUE)
        outFile_->WriteNodeHistoryHarmonic(*solHarmonic,  actFreqStep_, 
                                           actFrequency_, complexFormat_);
    
    } else
      Error("MechPDE: Only static, transient and harmonic results cna be written",
            __FILE__, __LINE__);
  
  }

  // ***********************************************************************
  //   Obtain information on desired output quantities from parameter file
  // ***********************************************************************
  void MechPDE::ReadStoreResults() {

    ENTER_FCN( "MechPDE::ReadStoreResults" );

    // Construct vectors for restricted parameter search
    StdVector<std::string> keyVec;
    StdVector<std::string> attrVec;
    StdVector<std::string> valVec;
    StdVector<std::string> regionNames;
    std::string quantity;
  
    // *****************************
    // Determine nodal results
    // ***************************** 
    StdVector<std::string> nodeValues;
    keyVec  = pdename_, "storeResults", "nodeResults", "region";
    attrVec = "", "", "type";  

    // --- Mechanic Displacement ---
    Enum2String(MECH_DISPLACEMENT, quantity);
    valVec = "", "", quantity;
    params->GetList( keyVec, attrVec, valVec, nodeValues);
    if (nodeValues.GetSize() > 0) {
      saveSol_ = TRUE;
      hasOutput_ = TRUE;
    }


    // --- Mechanic Velocity ---
    Enum2String(MECH_VELOCITY, quantity);
    valVec = "", "", quantity;
    params->GetList( keyVec, attrVec, valVec, nodeValues);
    if (nodeValues.GetSize() > 0) {
      saveDeriv1_ = TRUE;
      hasOutput_ = TRUE;
    
      // intialize corresponding storesolution object
      solDeriv1_.SetNumSolutions(1);
      solDeriv1_.SetNumNodes(numPDENodes_);
      solDeriv1_.SetSolutionType(MECH_VELOCITY);
      solDeriv1_.SetNumDofs(dim_);
      solDeriv1_.SetPtrEQNData(eqnData_, ptgrid_); 
      solDeriv1_.Init();
    }

    // --- Mechanic Acceleration ---
    Enum2String(MECH_ACCELERATION, quantity);
    valVec = "", "", quantity;
    params->GetList( keyVec, attrVec, valVec, nodeValues);
    if (nodeValues.GetSize() > 0) {
      saveDeriv2_ = TRUE;
      hasOutput_ = TRUE;
      
      // intialize corresponding storesolution object
      solDeriv2_.SetNumSolutions(1);
      solDeriv2_.SetNumNodes(numPDENodes_);
      solDeriv2_.SetSolutionType(MECH_ACCELERATION);
      solDeriv2_.SetNumDofs(dim_);
      solDeriv2_.SetPtrEQNData(eqnData_, ptgrid_);
      solDeriv2_.Init();
    }

    // *****************************
    // Determine element results
    // *****************************
    StdVector<std::string> elemResults;
    keyVec  = pdename_, "storeResults", "elemResults", "region";
    attrVec = "", "", "type";  

    // --- Mechanic Stress ---
    Enum2String(MECH_STRESS, quantity);
    valVec  = "", "", quantity;
    params->GetList( keyVec, attrVec, valVec, regionNames );
    ptgrid_->RegionNameToId( calcStress_, regionNames );
    
    // If the symbolic name is "all" compute electric field for all regions
    if ( calcStress_.GetSize() == 1 && calcStress_[0] == ALL_REGIONS ) {
      calcStress_ = subdoms_;
    }

    // Log to info file
    if ( calcStress_.GetSize() > 0 ) {
      hasOutput_ = TRUE;
      Info->PrintF( pdename_,
                    " Computing mechanical stress for regions:\n");
      for ( Integer k = 0; k < regionNames.GetSize(); k++ ) {
        Info->PrintF( pdename_, " %s\n", regionNames[k].c_str() );
      }
    }

    // --- Mechanic Energy ---
    Enum2String(MECH_ENERGY, quantity);
    valVec  = "", "", quantity;
    params->GetList( keyVec, attrVec, valVec, regionNames);
    ptgrid_->RegionNameToId( calcEnergy_, regionNames );

    // If the symbolic name is "all" compute electric field for all regions
    if ( calcEnergy_.GetSize() == 1 && calcEnergy_[0] == ALL_REGIONS ) {
      calcEnergy_ = subdoms_;
    }

    // Log to info file
    if ( calcEnergy_.GetSize() > 0 ) {
      hasOutput_ = TRUE;
      Info->PrintF( pdename_,
                    " Computing mechanical Energy for regions:\n");
      for ( Integer k = 0; k < regionNames.GetSize(); k++ ) {
        Info->PrintF( pdename_, " %s\n", regionNames[k].c_str() );
      }
    }

    // *****************************
    // Determine nodal history
    // *****************************
    StdVector<std::string> saveNodeHist; 
    keyVec  = pdename_, "storeResults", "nodeHistory", "saveNodes";
    attrVec = "", "", "type";

    // --- mechDisplacement ---
    Enum2String(MECH_DISPLACEMENT, quantity);
    valVec  = "", "", quantity;
    params->GetList( keyVec, attrVec, valVec, saveNodeHist );
  
    if (saveNodeHist.GetSize() > 0) {
      saveSolHist_ = TRUE;
      hasOutput_ = TRUE;
      Info->PrintF( pdename_, " Saving mechDisplacement for Nodes:\n" );
      for ( Integer k = 0; k < saveNodeHist.GetSize(); k++ ) {
        Info->PrintF( pdename_, " %s\n", saveNodeHist[k].c_str() );
      }
    }
  
    // --- mechVelocity ---
    Enum2String(MECH_VELOCITY, quantity);
    valVec  = "", "", quantity;
    params->GetList( keyVec, attrVec, valVec, saveNodeHist );
  
    if (saveNodeHist.GetSize() > 0) {
      saveDeriv1Hist_ = TRUE;
      hasOutput_ = TRUE;
      Info->PrintF( pdename_, " Saving mechVelocity for Nodes:\n" );
      for ( Integer k = 0; k < saveNodeHist.GetSize(); k++ ) {
        Info->PrintF( pdename_, " %s\n", saveNodeHist[k].c_str() );
      }
    }

    // --- mechAcceleration ---
    Enum2String(MECH_ACCELERATION, quantity);
    valVec  = "", "", quantity;
    params->GetList( keyVec, attrVec, valVec, saveNodeHist );
  
    if (saveNodeHist.GetSize() > 0) {
      saveDeriv1Hist_ = TRUE;
      hasOutput_ = TRUE;
      Info->PrintF( pdename_, " Saving mechAcceleration for Nodes:\n" );
      for ( Integer k = 0; k < saveNodeHist.GetSize(); k++ ) {
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

    if (saveElemHist.GetSize() >   0) {
      std::string errMsg = pdename_;
      errMsg += ": Saving history elements is not implemented yet!\n";
      errMsg += "Meanwhile you can use 'unvtool' to extract element data.";
      Error( errMsg.c_str(), __FILE__, __LINE__);
    }
  }


  // ************************************************************
  //   PostProcess
  // ************************************************************

  void MechPDE::PostProcess() {

    ENTER_FCN( "MechPDE::PostProcess" );

    //check for mechanical energy calculation
    if (calcEnergy_.GetSize() !=0 ) 
      CalcEnergy();
  
    if (calcStress_.GetSize() !=0 ) {

      //get the correct bilinearform
      ShortInt stressDim;
      Vector<Double> intPoint;

      if (subType_ == "planeStrain") {
        stressDim = 3;
      }

      else if (subType_ == "axi") {
        stressDim = 4;
      }
    
      else if (subType_ == "3d") {
        stressDim = 6;
      }
    
      else 
        Info->Error("StressOp: Unknown subtype in mech PDE! ",__FILE__,__LINE__);  
      
      
      // Resize solution arrays
      Stress_.SetNumSolutions(1);
      Stress_.SetSolutionType(MECH_STRESS);
      Stress_.SetNumElems(numElems_);
      Stress_.SetNumDofs(6);
      Stress_.SetPtrEQNData(eqnData_, ptgrid_);
      Stress_.Init(0);
      
      Vector<Double> elemStress, sortedStress;
      elemStress.Resize(stressDim);
      elemStress.Init(0);
      sortedStress.Resize(6);

      
      // loop over all subdomains
      for (Integer isd=0; isd<subdoms_.GetSize(); isd++) {
        
        MaterialData actSDMat(materialData_[isd]);
        MechStressStrain *stress;
        
        if (subType_ == "planeStrain") 
          stress = new MechStressStrainPlaneStrain(actSDMat);

        else if (subType_ == "axi") 
          stress = new MechStressStrainAxi(actSDMat);

        else if (subType_ == "3d") 
          stress = new MechStressStrain3D(actSDMat);

        
        // get vector of Elements of subdomains
        StdVector<Elem*> elemssd;     
        ptgrid_->GetVolElems( elemssd,subdoms_[isd] );
        
        // loop over elements of subdomain
        for (Integer iel=0; iel< elemssd.GetSize(); iel++) {
          Integer pdeElem = eqnData_->Mesh2PDEElem(elemssd[iel]->elemNum);
          elemssd[iel]->ptElem->GetCoordMidPoint(intPoint);
          //set element pointer
          BaseFE * ptEl = elemssd[iel]->ptElem;
          stress->SetElemPtr(ptEl);
          
          //set element solution        
          Matrix<Double> elSol;
          StdVector<Integer> connecth = elemssd[iel]->connect;
          sol_->GetElemSolutionAsMatrix(elSol, connecth);
          stress->SetActElemSol(elSol);
          
          //get coordinates of element
          Matrix<Double> ptCoord;
          GetElemCoords(connecth, ptCoord);
          
          Vector<Double> actStress;     
          
          //set the integration point
          stress->SetIntPoint(intPoint);

          //calculates the stress
          stress->CalcStressVec(elemStress,1,ptCoord);
          sortStresses(elemStress,sortedStress);
          Stress_.SetElemResult(pdeElem-1, sortedStress);
        }

        delete stress;
      }
    }
  }
  
  // ********************************************************
  //   Query parameter object for information about magnets
  // ********************************************************
  void MechPDE::ReadPreStressing() {

    ENTER_FCN( "MechPDE::ReadPreStressing" );

    StdVector<std::string> keyVec;
    StdVector<std::string> attrVec;
    StdVector<std::string> valVec;
    StdVector<std::string> regionNames;

    keyVec = pdename_, "preStressing", "preStress", "name";
    attrVec = "", "", "tag";
    valVec  = "", "", bcSequenceTag_;

    params->GetList(keyVec, attrVec, valVec, regionNames );
    ptgrid_->RegionNameToId( preStressDomain_, regionNames );

    if ( preStressDomain_.GetSize() > 0 ) {

      Info->PrintF( pdename_,
                    " Found prestressing in the following regions:\n" );

      Double tmpDir;

      // Construct vectors for restricted search parameter
      StdVector<std::string> keyVec;
      StdVector<std::string> attrVec;
      StdVector<std::string> valVec;
      attrVec = "", "", "name";

      // for each prestress domain ...
      for ( UInt k = 0; k < preStressDomain_.GetSize(); k++ ) {

        // ... read direction of magnetisation
        valVec = "", "", regionNames[k];

        keyVec  = pdename_, "preStressing", "preStress", "orientX";
        params->Get( keyVec, attrVec, valVec, tmpDir);
        preStressValX_.Push_back( tmpDir);

        keyVec  = pdename_, "preStressing", "preStress", "orientY";
        params->Get( keyVec, attrVec, valVec, tmpDir );
        preStressValY_.Push_back( tmpDir );

        keyVec  = pdename_, "preStressing", "preStress", "orientZ";
        params->Get( keyVec, attrVec, valVec, tmpDir );
        preStressValZ_.Push_back( tmpDir );

        // ... report name to logfile
        Info->PrintF( pdename_, "%s\n", regionNames[k].c_str());
      }
    }
  }

  void MechPDE::CalcEnergy()
  {
    ENTER_FCN( "MechPDE::CalcEnergy" );

    Matrix<Double> elemmat;  
    Matrix<Double> ptCoord;
    BaseFE         * ptElem;

    StdVector<Integer> connecth;
    Vector<double> help;

    Double totalE = 0;

    Integer i, j;
    Vector<Double> energy(subdoms_.GetSize());

    for (i=0; i<subdoms_.GetSize(); i++) {
    
      StdVector<Elem*> elemssd;
      ptgrid_->GetVolElems( elemssd,subdoms_[i] );

      //get material
      MaterialData actSDMat(materialData_[i]);
    
      energy[i] = 0;
      for (j=0; j < elemssd.GetSize(); j++) {  
        ptElem=elemssd[j]->ptElem;
        BaseForm * bilinear_stiff = GetStiffIntegrator(actSDMat);

        connecth=elemssd[j]->connect;
        GetElemCoords(connecth, ptCoord);
        bilinear_stiff->SetElemPtr(ptElem);
        bilinear_stiff->CalcElementMatrix(ptCoord, elemmat);

        Vector<Double> eldisp;
        sol_->GetElemSolution(eldisp, connecth);
        help = elemmat * eldisp;
        energy[i] += help * eldisp;

        delete bilinear_stiff;      

      }  

      totalE += energy[i];

    }

    std::string resulttype = "Mechanic Deformation Energy";
    std::string unit = "(Ws)";
    std::string analysis;
    Double analysisVal;
    if ( analysistype_ == HARMONIC ) {
      analysis    = "Frequency:";
      analysisVal = actFrequency_;
    }
    else {
      analysis    = "Time:";
      analysisVal = lasttimecalc_;
    }
    
    StdVector<std::string> regionNames;
    ptgrid_->RegionIdToName( regionNames, subdoms_ );
    Info->WriteResult(pdename_,  resulttype, regionNames, energy, unit,
                      analysis, analysisVal);

    StdVector<std::string> suball(1);
    Vector<Double> tmp(1);
    suball[0] = "Summe";
    tmp[0] = totalE;
    Info->WriteResult(pdename_,  resulttype, suball, tmp, unit,
                      analysis, analysisVal);
  }

} // end namespace CoupledField
