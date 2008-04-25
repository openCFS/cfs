// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "stokesFluidPDE.hh"

#include <sstream>
#include <iomanip>

#include "Forms/forms_header.hh"
#include "Forms/linElastInt.hh"
#include "Forms/massInt.hh"
#include "Forms/linPressureInt.hh"
#include "Driver/assemble.hh"
#include "newmark.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "Driver/solveStepStokesFluid.hh" 
#include "CoupledPDE/pdecoupling.hh"
#include "Domain/domain.hh"
#include "Utils/coordSystem.hh"
#include "Driver/assemble.hh"

#ifdef USE_SCRIPTING
#include "DataInOut/Scripting/cfsmessenger.hh" 
#endif

namespace CoupledField
{

  StokesFluidPDE::StokesFluidPDE( Grid * aptgrid, ParamNode* paramNode )
    :SinglePDE( aptgrid, paramNode )

  {

    isAlwaysStatic_ = true;

    pdename_          = "stokesFluid";
    pdematerialclass_ = FLOW;
    maxTimeDerivOrder_ = 0;


    // ****************************
    // DETERMINE GEOMETRY
    // ****************************

    // Get problem geometry and PDE subtype
    myParam_->Get("subType", subType_ );

    std::string probGeo;
    param->Get("domain")->Get("geometryType", probGeo );

    // Set number of degrees of freedom and
    // ensure that subtype fits to problem geometry
    UInt dofspernode = 0;
    if ( subType_ == "3d" && probGeo == "3d" ) {
      dofspernode = 7;
      Info->PrintF("", "=== 3D PROBLEM\n");
    }
    else if ( subType_ == "axi" && probGeo == "axi" ) {
      isaxi_ = true;
      dofspernode = 4;
      Info->PrintF("", "=== AXISYSMMETRIC PROBLEM\n");
    }
    else if ( subType_ == "plane" && probGeo == "plane" ) {
      dofspernode = 4;
      Info->PrintF("", "=== PLANE PROBLEM\n");
    }
    else
      {
        EXCEPTION( "Subtype '" << subType_ << "' of PDE '" <<  pdename_ 
                   << "' does not fit to problem geometry '"
                   << probGeo << "'" );
      }

  
   }

  StokesFluidPDE::~StokesFluidPDE()
  {
  }

  
  void StokesFluidPDE::InitNonLin()
  {

    EXCEPTION( "Method was not adapted during last change in parameterhandling interface ");

//     // ==============================================================
//     // NOTE: Currently we can only treat convective non-linearity and
//     //       we assume that for a stokes fluid PDE all regions either
//     //       are linear or non-linear!
//     // ==============================================================
//     StdVector<std::string> nonLinRegion;
//     params->GetList( "nonLinear", nonLinRegion, pdename_, "region" );
//     // Should not happen with validating parser, but beware!
//     if ( nonLinRegion.GetSize() == 0 ) {
//       nonLin_ = false;
//     }
//     else {
//       for ( UInt k = 1; k < nonLinRegion.GetSize(); k++ ) {
//         if ( nonLinRegion[k] != nonLinRegion[0] ) {
//           EXCEPTION( "Non-linearity should be the same for all regions!" );
//         }
//       }
//       nonLin_ = nonLinRegion[0] == "convection" ? true : false;
//     }

//     // In non-linear case determine type of line search strategy
//     if ( nonLin_ == true ) {

//       Info->PrintF( pdename_, "Non-linearity in %d regions\n",
//                     nonLinRegion.GetSize() );

//       // type of line search
//       params->Get( "type", lineSearch_, pdename_, "lineSearch" );

//       if ( lineSearch_ == "none" ) {
//         Info->PrintF( pdename_, "Performing no line search\n" );
//       }
//       else {
//         Info->PrintF( pdename_, "Will perform line search\n" );
//       }

//     }

//     // If no non-linearity we do not perform line search anyhow
//     else {
//       lineSearch_ = "none";
//     }

//     if( nonLin_ == true )
//       {
//         // incremental stopping criterion
//         params->Get( "incStopCrit", incStopCrit_, pdename_, "nonLinear" );

//         // residual stopping criterion
//         params->Get( "resStopCrit", residualStopCrit_, pdename_, "nonLinear" );
        
//         // maximal number of NL-iterations
//         params->Get("maxNumIters", nonLinMaxIter_, pdename_, "nonLinear");
//       }
  }

  void StokesFluidPDE::DefineIntegrators()
  {

    EXCEPTION( "Not implemented ");
//     // help variables for parameter checking
//     StdVector<std::string> keyVec;
//     StdVector<std::string> attrVec;
//     StdVector<std::string> valVec;
    
//     Double density, dynamicViscosity;

//     // ------------------------------------------
//     //   Get information on reduced integration
//     // ------------------------------------------
//     params->GetList( "reducedInt", reducedIntegration_, pdename_, "region" );

//     for (UInt actSD = 0; actSD < subdoms_.GetSize(); actSD++) {

//       // create new entity list
//       shared_ptr<ElemList> actSDList( new ElemList(ptgrid_ ) );
//       actSDList->SetRegion( subdoms_[actSD] );
      
//       BaseMaterial * actMat = materials_[subdoms_[actSD]];
//       actMat->GetScalar(density,DENSITY,REAL);
//       actMat->GetScalar(dynamicViscosity,DYNAMIC_VISCOSITY,REAL);
//       //dynamicViscosity = 1.85e-5;
      
//       Info->PrintF( pdename_, "density = %e\n", density);
//       Info->PrintF( pdename_, "dynamic Viscosity = %e\n", dynamicViscosity);
      
//       // ==============  add stiffness ======================================
//       BaseForm * bilinearStiff = GetStiffIntegrator(density, dynamicViscosity);
      
//       BiLinFormContext * actContext =
//         new BiLinFormContext(bilinearStiff, STIFFNESS );

//       actContext->SetPtPdes(this, this);
//       actContext->SetResults( results_[0], results_[0],
//                                actSDList, actSDList );
//       //actContext->SetReducedInt();    
	  
//       assemble_->AddBiLinearForm( actContext );

//       // ==============  add nonlinear stiffness ============================
//       if (nonLin_)
//         {
//           BaseForm *nLinPart = NULL;

//           if (subType_ == "3d")
//             {   
//               nLinPart = 
//                 new nLinStokesFluid3DInt_Convective(density, dynamicViscosity);    
//             }
//           else if (subType_ == "plane")
//             {
//               nLinPart = 
//                 new nLinStokesFluidPlaneInt_Convective(density, dynamicViscosity);    
//             }
//           else if (subType_ == "axi")
//             {
//               EXCEPTION("StokesFluid with convection needs to be implemented for axisymmetry! " );
//               //nLinPart =
//               //  new nLinStokesFluidAxiInt_Convective(density, dynamicViscosity);    
//             }
//           BiLinFormContext * stiffNLContext = 
//             new BiLinFormContext(nLinPart, STIFFNESS );
          
//           stiffNLContext->SetPtPdes(this, this);
//           stiffNLContext->SetResults( results_[0], results_[0],
//                                     actSDList, actSDList );
//           assemble_->AddBiLinearForm( stiffNLContext );
//         }
//         // ==================== RHS ===========================================
//         if (nonLin_) {
//           EXCEPTION( "Not working -> Gerhard has to derive the "
//                      <<  " 'nLinStokesFluid_linFormInt' from LinearForm!" );
          
//           // LinearForm * rhsSource = 
// //             new nLinStokesFluid_linFormInt(density, dynamicViscosity, isaxi_);
// //           LinearFormContext * rhsContext = 
// //             new LinearFormContext( rhsSource, nonLin_ );
// //           rhsContext->SetPtPde( this );
// //           rhsContext->SetResult( results_[0], actSDList );
// //           assemble_->AddLinearForm( rhsContext );
//         }
        
//         // Give result to equation numbering class
//       eqnMap_->AddResult( *results_[0], actSDList );
//     }
  }


  BaseForm *
  StokesFluidPDE::GetStiffIntegrator(Double density, Double dynamicViscosity)
  {

  
    BaseForm * bilinearStiff = NULL;

        if (subType_ == "plane")
          {
            bilinearStiff = new StokesFluidPlaneInt(density, dynamicViscosity);
          }
        else if (subType_ == "axi")
          {
            EXCEPTION("StokesFluid for axisymmetry needs to be implemented! " );
            //bilinearStiff = new StokesFluidAxiInt(density, dynamicViscosity);
          }
        else if (subType_ == "3d")
          {
            bilinearStiff = new StokesFluid3DInt(density, dynamicViscosity);
          }
        else 
          EXCEPTION("Unknown subtype in stokesFluid PDE! " );
    
    return bilinearStiff;
  }

 
  void StokesFluidPDE::DefineSolveStep()
  {

    solveStep_ = new SolveStepStokesFluid(*this);
  }





  // ======================================================
  // COUPLING SECTION
  // ======================================================

  void StokesFluidPDE::InitCoupling(PDECoupling * Coupling)
  {

    isIterCoupled_ = true;
    ptCoupling_   = Coupling;
    UInt velocityDofEndSpec;

    const UInt numCouplings = ptCoupling_->GetNumOutputCouplings();

    StdVector<StdVector<UInt> > elemNodeToCouplingNode_tmp;
    elemNodeToCouplingNode_.Resize(numCouplings);
    elemNodeToCouplingNode_.Init();

    if ( subType_ == "3d" ) {
      velocityDofEndSpec=3;
    }
    else if ( subType_ == "plane" ) {
      velocityDofEndSpec=2;
    }
    else {
      EXCEPTION( "Cannot handle subtype_ for ID_BC handling" );
    }
    
    // input couplings
    for (UInt i=0; i<ptCoupling_->GetNumInputCouplings(); i++)
      {
        Warning( "Not sure if this statement is needed (Andreas)!",
                 __FILE__, __LINE__ );
       //  // check for input mechanic velocity
//         if (ptCoupling_->GetInputQuantity(i) == MECH_VELOCITY)
//           {
//             numDirichletBCs_ += (velocityDofEndSpec * ptCoupling_->GetInputNumNodes(i));
//           }
      }

    for (UInt i = 0; i < numCouplings; i++) {

      // Intialize the memory of the coupling values
      ptCoupling_->CreateCouplingVector(i,isComplex_);

      if (ptCoupling_->GetOutputQuantity(i) == STOKESFLUID_FORCE) {
        // now since we need a incremental formulation, 
        //  initialize some necessary vectors
        isIncrFormulation_ = true;
      }
    }
  }

  void StokesFluidPDE::CalcInputCoupling() {

    EXCEPTION( "Method was not adapted during last change in parameterhandling interface ");

//     std::string errMsg;
//     StdVector<UInt> * nodes;
//     CFSVector * val;
//     Integer eqnNr;
//     UInt couplingDof;
//     UInt velocityDofEndSpec;
//     // Determine maximal allowed equation number for algebraic system
//     Integer maxAllowedEqn = (Integer)algsys_->GetDimension();

//     // at first, check if this PDE is iterative coupled
//     if (isIterCoupled_ == false)
//       return;

//     // Reset counter for boundary conditions
//     couplingBCsCounter_ = 0;
 

//     // Outer loop over all INPUT coupling terms
//     for (UInt i=0; i<ptCoupling_->GetNumInputCouplings(); i++) {

//       //    ptCoupling_ = &ptCoupling_[i];
//       ptCoupling_->GetInputValues(i, val);
//       couplingDof = ptCoupling_->GetInputDof(i);
    
//       // Up to now, Coupling is only possible with
//       // Real valued solutions
//       Vector<Double> const & help = dynamic_cast<Vector<Double>&>(*val);

//       switch(ptCoupling_->GetInputType(i)) {
        
//         // -------------------
//         // COORDINATE COUPLING
//         // -------------------
//       case COORD: 
        
//         // Set flag that the geometry has changed

//         ptCoupling_->GetInputNodes(i, nodes);
//         ptgrid_->SetNodeOffset( *nodes , help );
      
      
//         break;

//         // -------------------
//         // RHS COUPLING
//         // -------------------
//       case RHS:
//         ptCoupling_->GetInputNodes(i, nodes);
          
//         //for (UInt dof=0; dof<ptCoupling_->GetInputDof(i); dof++)
//         for ( UInt dof = 0; dof < couplingDof; dof++ ) {
//           for ( UInt j = 0; j < nodes->GetSize(); j++ ) {
//             eqnNr = eqnMap_->GetNodeEqn((*nodes)[j],dof+1);
//             if ( eqnNr != 0 && eqnNr <= maxAllowedEqn ) {
//               algsys_->SetNodeRHS( help[ dof + couplingDof * j ], pdeId_,
//                                    eqnNr, 1);
//             }

// #ifdef DEBUG
//             else if ( eqnNr > maxAllowedEqn ) {
//               (*debug) << "StokesFluidPDE::CalcInputCoupling: "
//                        << "(" << pdename_ << ") "
//                        << "Refused to pass "
//                        << "eqnNr = " << eqnNr << " to SetNodeRHS(), since "
//                        << "it execeeds numLastFreeDof = " << maxAllowedEqn
//                        << std::endl;
//             }
//             else if ( eqnNr == 0 ) {
//               (*debug) << "StokesFluidPDE::CalcInputCoupling: "
//                        << "(" << pdename_ << ") "
//                        << "Refused to pass "
//                        << "eqnNr = " << eqnNr << " to SetNodeRHS(), since "
//                        << "it is fixed by hom. Dirichlet BC" << std::endl;
//             }
// #endif

//           }
//         }
//         break;

//         // -----------------------
//         // InhomDirichlet COUPLING
//         // -----------------------
//       case ID_BC:

//         // How do we want to treat inhomogeneous Dirichlet boundary
//         // conditions?
//         {
//           bool usePenalty = true;
//           std::string aux;
//           StdVector<std::string> keyVec;
//           StdVector<std::string> attrVec;
//           StdVector<std::string> valVec;
//           keyVec  = "linearSystems", "system", "setup", "idbcHandling";
//           attrVec = "", "name", "";
//           valVec  = "", pdename_, "";
//           params->Get( keyVec, attrVec, valVec, aux );
//           usePenalty = aux == "penalty" ? true : false;
//           Info->PrintF( pdename_, "Treating IDBCs using '%s' approach\n",
//                         aux.c_str() );
//           if ( usePenalty == false ) {
//             EXCEPTION( "Cannot use inhom. Dirichlet coupling together with "
//                        << "IDBC elimination, since the equation numbering "
//                        << "object does currently not have the information "
//                        << "required to put those values at the end of the "
//                        << "equation number interval! Please set idbcHandling="
//                        << '"' << "penalty" << '"' << " in your xml-file" );
//           }
//         }

//         // Set flag that the boundary conditions have to be incorporated
//         updateCouplingBCs_ = true;

//         ptCoupling_->GetInputNodes(i, nodes);

//         if ( subType_ == "3d" ) {
//           velocityDofEndSpec=3;
// 	}
//         else if ( subType_ == "plane" ) {
//           velocityDofEndSpec=2;
// 	}
// 	else {
//           EXCEPTION( "Cannot handle subtype_ for ID_BC handling" );
// 	}

//         Info->PrintF( pdename_, "\nVelocity ID-BC: \n");
	
//         for ( UInt dof = 0; dof < velocityDofEndSpec; dof++ ) {
//           for ( UInt j = 0; j < nodes->GetSize();
//                 j++, couplingBCsCounter_++) {

//             eqnNr = eqnMap_->GetNodeEqn((*nodes)[j],dof+1 );

//             if (eqnNr==0) {
//               EXCEPTION( "The specified coupling node has no equation number" );
//             }

//             //Info->PrintF( pdename_, "Node=%d\tdof=%d\thelp[%d]=%e\n",
//             //                        (*nodes)[j], dof+1,
//             //                        dof+j*velocityDofEndSpec, help[dof+j*velocityDofEndSpec]);
	    
//             algsys_->SetDirichlet( couplingBCsCounter_ + 1, pdeId_, eqnNr,
//                                    help[dof+j*velocityDofEndSpec], 1 );
//           }
//         }
//         break;
          
//       case MAT:
//         EXCEPTION( "Not implemented yet" );
//         break;

//       }
//     }
  }


  void StokesFluidPDE::GetPresSolVecOfElement( Vector<Double>& elemSol,
                                       StdVector<UInt>& connecth ) {


    // stokesFluid pressure of element nodes
    elemSol.Resize(1 * connecth.GetSize());
    elemSol.Init(0);
    Integer eqnNr; 
    UInt presDof;

    NodeStoreSol<Double> * solhelp = 
      dynamic_cast<NodeStoreSol<Double>*>(sol_);
    Vector<Double> sol = solhelp->GetAlgSysVector();

    if (subType_ == "plane" || subType_ == "axi")
      {
        presDof = 3;
      }
    else if (subType_ == "3d")
      {
        presDof = 7;
      }
    else 
      EXCEPTION( "Unknown PDE subtype!" );

  
    for(UInt actNode=0; actNode<connecth.GetSize(); actNode++) {
      eqnNr = eqnMap_->GetNodeEqn(connecth[actNode],presDof);
      if (eqnNr!= 0) {
        elemSol[actNode] = sol[(abs(eqnNr-1))];
      }
      else {
        elemSol[actNode] = 0.0;
      }
    }
  }

  void StokesFluidPDE::CalcOutputCoupling()
  {
    UInt dof;
    SolutionType quantity;
    StdVector<Elem*> * couplingElems = NULL;
    StdVector<UInt> * couplingNodes = NULL;
    CFSVector * temp_values = NULL;
    //UInt regionCount = 0;
  
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
    
    
    Error( "Not working yet", __FILE__, __LINE__ );
//     Double density = 0.0;
//     Double dynamicViscosity = 0.0;
//     Double sign = 0.0;
//     Integer matIndex = -1;
//     Elem * ptVolElem = NULL;
//     Matrix<Double> ptCoord, elemMat;
//     Vector<Double> sol, forceOnElem, normal;
    
//     elemCouplingSols.Init(0.0);
    
//     for (UInt actElem=0; actElem<couplingElems->GetSize(); actElem++) {
      
//       // Perform cast from volume element to surface element, since
//       // mech-acou coupling makes only sense on surface elements
//       SurfElem * actCoupleElem = 
//         dynamic_cast<SurfElem*> ((*couplingElems)[actElem]);

//       if (actCoupleElem == NULL) {
//         Error( "No elements found for coupling!", __FILE__, __LINE__ );
//       }
      
//       BaseFE * ptElem = actCoupleElem->ptElem;
//       StdVector<UInt> & connecth = actCoupleElem->connect;
//       GetElemCoords(connecth, ptCoord);
      
//       // Try to find according region for first neighbouring volume
//       // element of the surface element
//       matIndex = subdoms_.Find(actCoupleElem->ptVolElem1->regionId);
      
//       // If first volume element does not belong to stokesFluid PDE, try the
//       // second one
//       if ( matIndex == -1 ) {
//         matIndex = subdoms_.Find(actCoupleElem->ptVolElem2->regionId);
//         ptVolElem = actCoupleElem->ptVolElem2;
//         sign = actCoupleElem->normalSign;
//       } else {
//         ptVolElem = actCoupleElem->ptVolElem1;
//         sign = -1.0 * actCoupleElem->normalSign;
//       }
      
//       if ( matIndex == -1) {
//         (*error) << "StokesFluidPDE::CalcMechCouplingRHS: The two volume "
//                  << "element neighbours of surface element Nr. "
//                  << actCoupleElem->elemNum << " do not belong to my regions!";
//         Error( __FILE__, __LINE__ );
//       }
      
//       // Assign correct density and dynamicViscosity
//       materials_[subdoms_[matIndex]]->
//         GetScalar(density,DENSITY,REAL);
//       materials_[subdoms_[matIndex]]->
//         GetScalar(dynamicViscosity,DYNAMIC_VISCOSITY,REAL);
//       //dynamicViscosity = 1.85e-5;
      
//       BaseForm * bilinear_mass = new MassInt(ptElem, 1.0, isaxi_);
//       bilinear_mass->CalcElementMatrix(ptCoord, elemMat);
//       delete bilinear_mass;     

//       // acoustic
//       //GetDerivSolVecOfElement(sol, connecth);

//       //stokesFluid      
//       GetPresSolVecOfElement(sol, connecth);
     
      
//       forceOnElem = elemMat * sol;
      
//       // force has to be added on RHS with negative sign
//       forceOnElem *= - 1.0;
      
//       // the normal vector points outwards of the MECHANICAL domain
//       // (see. Kaltenbacher, "Num. Sim. of Mechatr. Act. & Sens." chapter 8.2)
//       ptgrid_->CalcSurfNormal(normal, *actCoupleElem);
//       normal *= sign;

//       for (UInt actNode=0; actNode<ptCoord.GetSizeCol(); actNode++) {
//         UInt nodePos = 0;
          
//         while(connecth[actNode] != couplingNodes[nodePos] &&
//               nodePos < couplingNodes.GetSize()) {
//           nodePos++;
//         }
          
//         for (UInt actDof=0; actDof < couplingdof ; actDof++) {
//           elemCouplingSols[nodePos*couplingdof +actDof] += 
//             forceOnElem[actNode] * normal[actDof];
//         }
//       }
//     }
  }

  bool StokesFluidPDE::HasOutput(SolutionType output)
  {

    if (output == STOKESFLUID_FORCE)
      return true;

    return false;
  }



  // ======================================================
  // TIME STEPPING SECTION
  // ======================================================


  void StokesFluidPDE::InitTimeStepping()
  {
  }



 

  // ***********************************************************************
  //   Obtain information on desired output quantities from parameter file
  // ***********************************************************************
  void StokesFluidPDE::DefineAvailResults() {

    Error( "Implement!", __FILE__, __LINE__ );

    // =====================================================================
    // set solution information
    // =====================================================================
    
    // Create new resultDof object
    shared_ptr<ResultInfo> res1(new ResultInfo);
    shared_ptr<AnsatzFct> fct(new LagrangeFct);
    res1->resultType = STOKESFLUID_VEL_PRES_VORT;
    res1->definedOn = ResultInfo::NODE;
    res1->fctType = fct;
    if( subType_ == "3d" ) {
      res1->dofNames = "vx", "vy", "vz", "omegax", "omegay", "omegaz", "p";
    } else if ( subType_ == "axi" ) {
      res1->dofNames = "vr", "vphi", "omega", "p";
    } else if( subType_ == "plane") {
      res1->dofNames = "vx", "vy", "p", "omega";
    } else {
      (*error) << "StokesFluidPDE: subType '" << subType_ << "' not known!";
      Error( __FILE__, __LINE__ );
    }
    res1->entryType = ResultInfo::VECTOR;
    results_.Push_back( res1 );    
    availResults_.insert( res1 );


  }

  // ************************************************************
  //   CalcResults
  // ************************************************************

  void StokesFluidPDE::CalcResults( shared_ptr<BaseResult> res ) {
    Error( "Implement!", __FILE__, __LINE__ );
  }
} // end namespace CoupledField


