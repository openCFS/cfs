// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <fstream>
#include <iostream>
#include <string>
#include <math.h>

#include "acousticPDE.hh" 
#include "Forms/forms_header.hh"
#include "Forms/abcInt.hh"
#include "Forms/linNeumannInt.hh"
#include "Forms/nonConformingInt.hh"
#include "Estimator/spaceerror.hh"
#include "newmark.hh"
#include "newmarkFracDamp.hh"
#include "DataInOut/WriteInfo.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/Logging/cfslog.hh"
#include "Utils/mathfunctions.hh"
#include "Utils/nodestoresol.hh"
#include "Driver/solveStepAcoustic.hh"
#include "CoupledPDE/pdecoupling.hh"
#include "Domain/ansatzFct.hh"
#include "Driver/assemble.hh"
#include "Domain/domain.hh"

#ifdef USE_SCRIPTING
#include "DataInOut/Scripting/cfsmessenger.hh" 
#endif

#ifdef PARALLEL
#include <mpi.h>
#endif

namespace CoupledField {

  DECLARE_LOG(acoupde)
  DEFINE_LOG(acoupde, "pde.acoustic")


  // =========================================================================
  // set solution information
  // =========================================================================
  AcousticPDE::AcousticPDE( Grid* aptgrid, ParamNode* paramNode )
    : SinglePDE( aptgrid, paramNode ) {

    ENTER_FCN( "AcousticPDE::AcousticPDE" );

    pdename_          = "acoustic";
    pdematerialclass_ = FLUID;
    maxTimeDerivOrder_ = 2;
    fracDamping_ = false;

    isMechCoupled_ = false;
    isNrbcCoupled_ = false;
    isBubbleCoupled_ = false;

    mHandle_ =  domain->GetMathParser()->GetNewHandle();
//      MathParser * mParser =  domain->GetMathParser();
//      mParser->SetExpr( mHandle_, "t" );

    // PDE formulation either in acoustic potential or pressure
    std::string str = myParam_->Get("formulation")->AsString();
    String2Enum( str, formulation_ );
    str = "Using * " + str + " as state variable in formulation of PDE\n";
    Info->PrintF( pdename_, str.c_str() );


    // timestepping formulation
    str = "";
    myParam_->Get( "timeSteppingFormulation", str,  false );
    if ( str == "effMassMatrix" ) {
      effectiveMass_ = true;
      Info->PrintF( pdename_, 
                    "      * effective mass matrix timestepping\n");
    } 
    else if ( str == "diagMassMatrix" ) {
      diagMass_      = true;
      effectiveMass_ = true;
      Info->PrintF( pdename_, 
                    "      * diagonal mass matrix in explicit timestepping\n");
    } 
    else {
      effectiveMass_ = false;
      Info->PrintF( pdename_, 
                    "      * effective stiffness matrix timestepping\n");
    }
    
    //To check if acoustic is coupled with nrbc and set isNrbcCoupled
    myParam_->Get( "isCoupledNrbc", isNrbcCoupled_, false );

  }


  // *********************************************
  //   Check what type of damping should be used
  // *********************************************
  void AcousticPDE::ReadDampingInformation( ) {

    ENTER_FCN( "AcousticPDE::ReadDampingInformation" );

    fracMemory_ = 0;
    Integer firstFrac=-1;

    // try to get dampingList
    ParamNode * dampListNode = myParam_->Get( "dampingList", false );
    if( !dampListNode) return;
    
    // get specific damping nodes
    StdVector<ParamNode*> dampNodes = dampListNode->GetChildren();
    std::map<std::string, DampingType> idDampType;
    for( UInt i = 0; i < dampNodes.GetSize(); i++ ) {

      std::string dampString = dampNodes[i]->GetName();
      std::string actId = dampNodes[i]->Get("id")->AsString();

      // determine type of damping
      DampingType actType;
      String2Enum( dampString, actType );
      
      // make special things for fractional damping
      if( actType == FRACTIONAL ) {
        fracDamping_ = true;

        // Find first region containing fractional damping
        if ( firstFrac < 0 )
          firstFrac = i;
        
        // Gather additional information for fractional damping model
        std::string fracAlg = "gl";
        dampNodes[i]->Get( "fracAlg", fracAlg, false );
        
        std::string interpol = "no";
        dampNodes[i]->Get( "interpolation", interpol, false );
        
        UInt fracMem = 1;
        dampNodes[i]->Get( "fracMemory", fracMem, false );
        
        if  ( fracAlg == "gl" ) {
          if (interpol == "no" )
            actType = FRACTIONAL_GL;
          else {
            actType= FRACTIONAL_GL_INT;
          }
        }
        else if ( fracAlg == "blank" ) {
          if (interpol == "no" )
            actType = FRACTIONAL_BLANK;
          else {
            actType= FRACTIONAL_BLANK_INT;
          }
        }
        
        // up to now take maximum of specified fracMemory values
        if ( fracMem > fracMemory_ )
          fracMemory_ = fracMem;
      }
      // store damping type string
      idDampType[actId] = actType;
      
    }
    if ( fracDamping_== true ) {
      Info->PrintF(pdename_, "Memory size for fractional damping  is: %d\n",
                   fracMemory_ );
    }
    
    // Run over all region and set entry in "regionNonLinId"
    StdVector<ParamNode*> regionNodes = 
      myParam_->Get("regionList")->GetChildren();
    
    RegionIdType actRegionId;
    std::string actRegionName, actDampingId;
    
    if( regionNodes.GetSize() > 0 ) {
      Info->PrintF( pdename_, "Damping in following region(s)\n" );
    }
    
    for (UInt k = 0; k < regionNodes.GetSize(); k++) {
      regionNodes[k]->Get( "name", actRegionName );
      regionNodes[k]->Get( "dampingId", actDampingId );
      if( actDampingId == "" )
        continue;
      
      actRegionId = ptgrid_->RegionNameToId( actRegionName );
      
      // Check actDampingId was already registerd
      if( idDampType.count( actDampingId ) == 0 ) {
        EXCEPTION( "Damping with id '" << actDampingId 
                   << "' was not defined in 'dampingList'" );
      }
      
      dampingList_[actRegionId] = idDampType[actDampingId];
      
      // Log to info file
      std::string dampString;
      Enum2String( dampingList_[actRegionId], dampString );
      
      if( dampingList_[actRegionId] == FRACTIONAL_GL ) {
        dampString += "( Gruenwald-Letnikov algorithm )";
      }
      if( dampingList_[actRegionId] == FRACTIONAL_BLANK ) {
        dampString += "( Blanks algorithm )";
      }
      if( dampingList_[actRegionId] == FRACTIONAL_GL_INT ||
          dampingList_[actRegionId] == FRACTIONAL_BLANK_INT ) {
        dampString += "(linear interpol. of single past values)";
      }
      Info->PrintF( pdename_, " %s: %s\n", actRegionName.c_str(), 
                    dampString.c_str() );      
    }
    Info->PrintF( pdename_, "\n" );
  }
  
  void AcousticPDE::ReadSpecialBCs() {
    ENTER_FCN( "AcousticPDE::ReadSpecialBCs" );

    // ***************************************************************
    //   If no other damping type is specified and we have absorbing
    //   boundary conditions, then use ABCDAMP
    // ***************************************************************
    StdVector<std::string> auxVec;
    absorbingBCs_ = false;
    ParamNode * bcNode = myParam_->Get( "bcsAndLoads", false );
    if( bcNode ) {
      StdVector<ParamNode*> abcNodes = bcNode->GetList( "absorbingBCs" );
      
      RegionIdType actRegion;
      for( UInt i = 0; i < abcNodes.GetSize(); i++ ) {
        std::string regionName = abcNodes[i]->Get("name")->AsString(); 
        actRegion = ptgrid_->RegionNameToId( regionName );
        absBCs_.Push_back( actRegion );
        absorbingBCs_ = true;
        Info->PrintF( pdename_, 
                      "Apply Absorbing Boundary Conditions on surfaceRegion '%s'\n",
                      regionName.c_str() );
      }
    }
  }


  // *************************************************************
  //   Check what type of nonlinear PDE formulation should be used
  // *************************************************************
  void  AcousticPDE::InitNonLin() {
    ENTER_FCN( "AcousticPDE::InitNonLin" );
    
    nonLin_ = false; 
    // Check, if "nonLinList" is present
    ParamNode * nonLinListNode = myParam_->Get("nonLinList", false );
    if( !nonLinListNode) 
      return;
    
    // Get nonlinear types
    StdVector<ParamNode*> nonLinNodes = nonLinListNode->GetChildren();
    for( UInt i = 0; i < nonLinNodes.GetSize(); i++ ) {

      std::string actTypeString = nonLinNodes[i]->GetName();
      std::string actId = nonLinNodes[i]->Get("id")->AsString();

      NonLinType actType;
      String2Enum( actTypeString, actType );
      
      // check for specific types
      if( actType == WESTERVELT ) {
        nonLin_ = true;
        if ( formulation_ == ACOU_POTENTIAL ) 
          EXCEPTION( "Acoustic potential formulation not supported for "
                     << "Westervelt equation" );
      } else if( actType == KUZNETSOV ) {
        nonLin_ = true;
        if ( formulation_ == ACOU_PRESSURE ) {
          EXCEPTION( "Acoustic pressure formulation not supported for"
                     << "Kuznetsov equation!" );
        }
      } else if( actType == VARIABLE_SOS_CN1 ||
                 actType == VARIABLE_SOS_CN2 ||
                 actType == VARIABLE_SOS_CN2Mean ) {
        variableSpeedOfSoundCN_ = true;
      }
      nonLinIdType_[actId] = actType;
    }
    
   // Run over all region and set entry in "regionNonLinId"
    StdVector<ParamNode*> regionNodes = 
      myParam_->Get("regionList")->GetChildren();
    
    RegionIdType actRegionId;
    std::string actRegionName, actNonLinId;
    
    if( regionNodes.GetSize() > 0 ) {
      Info->PrintF( pdename_, "Non-linearity in following region(s)\n" );
    }
    for( UInt i = 0; i < regionNodes.GetSize(); i++ ) {
      
      regionNodes[i]->Get( "name", actRegionName );
      regionNodes[i]->Get( "nonLinId", actNonLinId );
      
      if( actNonLinId == "" )
        continue;
      
      actRegionId = ptgrid_->RegionNameToId( actRegionName );
      
      // Check nonLinId was already registerd
      if( nonLinIdType_.find( actNonLinId) == nonLinIdType_.end() ) {
        EXCEPTION( "NonLinearity with id '" << actNonLinId 
                   << "' was not defined in 'nonLinList'" );
      }
      
      regionNonLinId_[actRegionId] = actNonLinId;
      regionNonLinType_[actRegionId] = nonLinIdType_[actNonLinId];
      
      // Log to info file
      std::string nonLinString;
      Enum2String( nonLinIdType_[actNonLinId], nonLinString );
      Info->PrintF( pdename_, " %s: %s\n", actRegionName.c_str(), 
                    nonLinString.c_str() );
      
    }
    
    Info->PrintF( pdename_, "\n" );
  }


  void AcousticPDE::DefineIntegrators() {

    ENTER_FCN( "AcousticPDE::DefineIntegerators" );

    // =======================================================================
    //  Read flow data
    // =======================================================================
    ReadFlowData();

    Double density, compressibility, c0;
    Double coeffmass, coeffdamp;

    for (UInt actSD = 0; actSD < subdoms_.GetSize(); actSD++) {

      // create new entity list
      shared_ptr<ElemList> actSDList( new ElemList(ptgrid_ ) );
      actSDList->SetRegion( subdoms_[actSD] );

      // get current region name and get grip of paramNode
      RegionIdType actRegion = subdoms_[actSD];
      std::string actRegionName;
      actRegionName = ptgrid_->RegionIdToName( actRegion );
      ParamNode * actRegionNode = 
        myParam_->Get("regionList")->Get( "region", "name", actRegionName );

      // ********************************************************************
      //   Attention:
      //   In AcousticPDE ALL integrators are multiplied with density!
      // ********************************************************************

      materials_[actRegion]->GetScalar(density,DENSITY,REAL);
      materials_[actRegion]->GetScalar(compressibility,ACOU_BULK_MODULUS,REAL);
      c0 = sqrt(compressibility/density);

      // if pde couples with mechanic, we have to multiply the density by -1
      if ( isMechCoupled_ == true && formulation_ !=  ACOU_PRESSURE) {
        density *= -1.0;
      }

      // ********************************************************************
      //   Linear wave equation
      // ********************************************************************
      
      // Check for Perfectly matchec layers
      if ( dampingList_[actRegion] == PML ) {
        if ( analysistype_ != HARMONIC ) {
          EXCEPTION( "PML just supported for Harmonic-Analysis" );
        }
        
        //read data for PML layer
        
        //type of PML damping
        std::string dampingTypePML;
        
        // inner / outer region
        Matrix<Double> inner;
        Matrix<Double> outer;
        
        //damping factor
        Double dampPML;

        std::string id = actRegionNode->Get("dampingId")->AsString();
        ParamNode * pmlNode = myParam_->Get("dampingList")->Get("pml", "id", id);
        ReadDataPML(dampingTypePML, inner, dampPML, pmlNode );
        dampPML *= c0;
        
        GetPMLLayerData(inner, outer, actRegion);
        
        //====================================================================
        //	 stiffness integrator for PML
        //====================================================================

        std::string formsType = "laplaceInt";

        //set real part
        BaseForm * bilinearStiffReal = 
          new PMLInt(formsType, density, dampingTypePML, dampPML, isaxi_);

        bilinearStiffReal->SetPosPML(inner,outer);
        DataType matType = REAL;
        bilinearStiffReal->SetMatDataType(matType);

        BiLinFormContext * stiffContextReal = 
          new BiLinFormContext( bilinearStiffReal, STIFFNESS );

        stiffContextReal->SetPtPdes(this, this);
        stiffContextReal->SetResults( results_[0], results_[0],
                                      actSDList, actSDList );
        stiffContextReal->SetEntryType(matType);
        assemble_->AddBiLinearForm( stiffContextReal);

        //set imaginary part
        BaseForm * bilinearStiffImag = 
          new PMLInt( formsType, density, dampingTypePML, dampPML, isaxi_ );

        bilinearStiffImag->SetPosPML(inner,outer);
        matType = IMAG;
        bilinearStiffImag->SetMatDataType( matType );

        BiLinFormContext * stiffContextImag = 
          new BiLinFormContext( bilinearStiffImag, STIFFNESS );


        stiffContextImag->SetPtPdes(this, this);
        stiffContextImag->SetResults( results_[0], results_[0],
                                      actSDList, actSDList );
        stiffContextImag->SetEntryType(matType);
        assemble_->AddBiLinearForm( stiffContextImag );


        //====================================================================
        //	 mass integrator for PML
        //====================================================================
        formsType = "massInt";

        //	Double dampMass = (2.0/3.0)*density/c0;
        Double massFactor = density/(c0*c0);

        //set real part
        BaseForm * bilinearMassReal = 
          new PMLInt( formsType, massFactor, dampingTypePML, dampPML, isaxi_ );

        bilinearMassReal->SetPosPML(inner,outer);
        matType = REAL;
        bilinearMassReal->SetMatDataType(matType);

        BiLinFormContext * massContextReal = 
          new BiLinFormContext( bilinearMassReal, MASS);

        massContextReal->SetPtPdes(this, this);
        massContextReal->SetResults( results_[0], results_[0],
                                     actSDList, actSDList );
        massContextReal->SetEntryType( matType );
        assemble_->AddBiLinearForm( massContextReal );

        //set imaginary part
        BaseForm * bilinearMassImag = 
          new PMLInt( formsType, massFactor, dampingTypePML, dampPML, isaxi_ );

        bilinearMassImag->SetPosPML( inner, outer );
        matType = IMAG;
        bilinearMassImag->SetMatDataType(matType);

        BiLinFormContext * massContextImag = 
          new BiLinFormContext( bilinearMassImag, MASS );

        massContextImag->SetPtPdes( this, this );
        massContextImag->SetResults( results_[0], results_[0],
                                     actSDList, actSDList );
        massContextImag->SetEntryType( matType );
        assemble_->AddBiLinearForm( massContextImag );

        // Finally add the stiffness/mass integrators
        assemble_->AddBiLinearForm( stiffContextReal );
        assemble_->AddBiLinearForm( stiffContextImag );
        assemble_->AddBiLinearForm( massContextReal );
        assemble_->AddBiLinearForm( massContextImag );
      } // end of pml part

      else {
        // stiffness integrator 
        BaseForm * bilinearStiff = new LaplaceInt( density, isaxi_ );        
        BiLinFormContext * stiffContext = 
          new BiLinFormContext( bilinearStiff, STIFFNESS );

        stiffContext->SetResults( results_[0], results_[0],
                                  actSDList, actSDList );
        stiffContext->SetPtPdes( this, this );

        // mass integrator
        coeffmass = density / (c0*c0);

        MassInt* bilinearMass  = new MassInt(coeffmass, 1, isaxi_);
	if ( diagMass_ ) {
	  // diagonal mass matrix
	  bilinearMass->SetDiagMass();
	}

        BiLinFormContext * massContext = 
          new BiLinFormContext( bilinearMass, MASS );
        massContext->SetResults( results_[0], results_[0],
                                 actSDList, actSDList );
        massContext->SetPtPdes(this, this);

	//check  for Flow Data
        if ( regionFlowNodes_.count( actRegion) > 0 ) {
          if ( formulation_ == ACOU_PRESSURE )
	    EXCEPTION("Pierce-Equation just possible in velocity potential formulation" );
          
	  //read flow data 
          ParamNode * flowNode = regionFlowNodes_[actRegion];
	  SimpleFlow* flowData = new SimpleFlow();
	  flowData->ReadFlowData( flowNode, dim_);

	  Double coeffPierceStiff = -density / (c0*c0);
	  BaseForm * bilinearPierceStiff  = new PierceStiffInt(coeffPierceStiff, 
							       flowData,
							       isaxi_);
	  BiLinFormContext * pierceStiffContext = 
	    new BiLinFormContext( bilinearPierceStiff, STIFFNESS );
	  pierceStiffContext->SetResults( results_[0], results_[0],
					  actSDList, actSDList );
	  pierceStiffContext->SetPtPdes(this, this);
	  assemble_->AddBiLinearForm( pierceStiffContext );    


	  Double coeffPierceDamp = 2.0 * density / (c0*c0);
	  BaseForm * bilinearPierceDamp  = new PierceDampInt(coeffPierceDamp, 
							     flowData,
							     isaxi_);
	  BiLinFormContext * pierceDampContext = 
	    new BiLinFormContext( bilinearPierceDamp, DAMPING );
	  pierceDampContext->SetResults( results_[0], results_[0],
					 actSDList, actSDList );
	  pierceDampContext->SetPtPdes(this, this);
	  assemble_->AddBiLinearForm( pierceDampContext );        	
        }

        // ********************************************************************
        //   Additional terms for damping
        // ********************************************************************

        //check  for damping layer
        if ( dampingList_[actRegion] == DAMPLAYER ) {
	  
          //type of damping fnc
          std::string dampFnc;

          std::string id = actRegionNode->Get("dampingId")->AsString();
          ParamNode * dampLayerNode = myParam_->Get("dampingList")->Get("dampLayer", "id", id);

          //damping data
          Double dampFactor, dampFactorMax, startRadius, stopRadius;
          Vector<Double> mPoint;
          ReadDataDampLayer( dampFnc, mPoint, dampFactor, dampFactorMax, 
                             startRadius, stopRadius, dampLayerNode );

          //get the Rayleigh material parameters
          Double alpha, beta;
          materials_[actRegion]->GetScalar(alpha,RAYLEIGH_ALPHA,REAL);
          materials_[actRegion]->GetScalar(beta,RAYLEIGH_BETA,REAL);
    
          // stiffness part
          stiffContext->SetSecDestMat( DAMPING, beta );
          stiffContext->SetDampLayer(dampFnc, mPoint, dampFactor, 
                                     dampFactorMax, startRadius, 
                                     stopRadius);	    
          // mass part
          massContext->SetSecDestMat( DAMPING, alpha );
          massContext->SetDampLayer(dampFnc, mPoint, dampFactor, 
                                    dampFactorMax, startRadius, 
                                    stopRadius);
        }

        if ( dampingList_.size() > 0 ) {
	  
          // We check, if damping has been specified for all regions.
          if ( dampingList_.size() != subdoms_.GetSize() ) {
            (*warning) << "Mismatch between dampingList_ and subdoms_!"
                       << "Size(dampingList_): " << dampingList_.size()
                       << "Size(subdoms_): " << subdoms_.GetSize();
            Warning(__FILE__, __LINE__);  
          }
	  
          if (dampingList_[actRegion] == RAYLEIGH) {
            // This works even after assemble_->AddIntegrator() is executed
            //   because of the pointers...
            Double alpha, beta;
            materials_[actRegion]->GetScalar(alpha,RAYLEIGH_ALPHA,REAL);
            materials_[actRegion]->GetScalar(beta,RAYLEIGH_BETA,REAL);
    
            // stiffness part
            stiffContext->SetSecDestMat( DAMPING, beta );
	    
            // mass part
            massContext->SetSecDestMat( DAMPING, alpha );
          }
	  
	  
          else if ( dampingList_[actRegion] == THERMOVISCOUS ) {
            Double viscousAlpha;
            materials_[actRegion]->GetScalar(viscousAlpha, ACOU_ALPHA,REAL);

            coeffdamp  =  density * 2.0 * viscousAlpha * c0;
            BaseForm * bilinearStiff  = new LaplaceInt(coeffdamp, isaxi_);  
            BiLinFormContext * dampContext = 
              new BiLinFormContext(bilinearStiff, DAMPING );
            dampContext->SetResults( results_[0], results_[0],
                                     actSDList, actSDList );
            dampContext->SetPtPdes(this, this);   
            assemble_->AddBiLinearForm( dampContext );
          }
	  
          else if ( dampingList_[actRegion] == FRACTIONAL_GL ||
                    dampingList_[actRegion] == FRACTIONAL_GL_INT ) {

            Double fracAlpha, fracExp;
            materials_[actRegion]->GetScalar(fracAlpha,ACOU_ALPHA,REAL);
            materials_[actRegion]->GetScalar(fracExp,FRACTIONAL_EXPONENT,REAL);

            coeffdamp = - density * 2.0 * fracAlpha / c0 / sin((fracExp-1.0)*PI/2.0);

            MassInt * bilinearDamp  = 
              new MassInt(coeffdamp, 1, isaxi_);
            Double fracDampCoeff = GetFracDampMatrixCoeff( actRegion );
            bilinearDamp->SetFactor( GenStr(fracDampCoeff ) );

            // formulation using DAMPING matrix
            // adapt NewmarkFracDamp::Init and StdPDE::GetFracDampMatrixCoeff
            // IntegratorDescriptor * dampIntDescr = 
            //   new IntegratorDescriptor(bilinearDamp, DAMPING);
	    
            // two matrices formulation
            // added to STIFFNESS matrix because, because 
            //   matrix_factors[STIFFNESS] = 1.0
            BiLinFormContext * dampContext = 
              new BiLinFormContext( bilinearDamp, STIFFNESS );
            dampContext->SetResults( results_[0], results_[0],
                                     actSDList, actSDList );
            dampContext->SetPtPdes(this, this);
            assemble_->AddBiLinearForm( dampContext );
          }
	  
          else if  ( dampingList_[actRegion] == FRACTIONAL_BLANK ||
                     dampingList_[actRegion] == FRACTIONAL_BLANK_INT ) {

            Double fracAlpha, fracExp;
            materials_[actRegion]->GetScalar(fracAlpha,ACOU_ALPHA,REAL);
            materials_[actRegion]->GetScalar(fracExp,FRACTIONAL_EXPONENT,REAL);

            coeffdamp =  - density * 2.0 * fracAlpha / c0 / sin((fracExp-1.0)*PI/2.0);
            // prefactor of blank alg
            coeffdamp *= exp(-gammaln(1.0- (fracExp- 1.0)) ); 
            // weight factor of index 0
            coeffdamp *= 1.0/(1.0- (fracExp- 1.0));           
            MassInt * bilinearDamp  = 
              new MassInt(coeffdamp, 1, isaxi_);
            Double fracDampCoeff = GetFracDampMatrixCoeff( actRegion );
            bilinearDamp->SetFactor( GenStr(fracDampCoeff) );

            // formulation using DAMPING matrix
            // adapt NewmarkFracDamp::Init and StdPDE::GetFracDampMatrixCoeff
            // IntegratorDescriptor * dampIntDescr = 
            //   new IntegratorDescriptor(bilinearDamp, DAMPING);
	    
            // two matrices formulation
            // added to STIFFNEss matrix because, because 
            //   matrix_factors[STIFFNESS] = 1.0
            BiLinFormContext * dampContext = 
              new BiLinFormContext( bilinearDamp, STIFFNESS );
            dampContext->SetResults( results_[0], results_[0],
                                     actSDList, actSDList );
            dampContext->SetPtPdes(this, this);
            assemble_->AddBiLinearForm( dampContext );

          }
        }
        
        // **********************************************************************
        //   Additional terms for cavitation computation
        // **********************************************************************
        
        // check, if in the current sequenceStep a pde with name "bubble" is 
        // defined
        ParamNode * bubbleNode = 
          param->Get("sequenceStep", "index", GenStr(sequenceStep_) )
          ->Get("pdeList")->Get("bubble", false);

        if( bubbleNode ){
          isBubbleCoupled_ = true;
          Info->PrintF(  pdename_, " will use bubble-bilinearform\n");
          
          
          
	  //get paramters and starting values for the computation of coefficients
	  Double sonicVel, densityforbubble, viscosity, bubbleDensity, frequency;
          bubbleNode->Get("bubbleDynamic")->Get("density", density );
          bubbleNode->Get("bubbleDynamic")->Get("sonicVel", sonicVel );
	  bubbleNode->Get("bubbleDynamic")->Get( "viscosity", viscosity );
	  bubbleNode->Get("bubbles")->Get( "bubbleNumDensity", bubbleDensity );
          bubbleNode->Get("excitation")->Get( "frequency", frequency );

	  //create two additional matricies due to the influence of the bubbles
	  //onto the wave behaviour
	  BubbleStiffInt * bubbleStiff = new BubbleStiffInt( density,
							     densityforbubble,
							     sonicVel, viscosity,
							     bubbleDensity,
							     frequency, isaxi_);        
	  BiLinFormContext * bubbleStiffContext = 
	    new BiLinFormContext( bubbleStiff, STIFFNESS );

	  bubbleStiffContext->SetResults( results_[0], results_[0],
				    actSDList, actSDList );
	  bubbleStiffContext->SetPtPdes( this, this );
	  
	  assemble_->AddBiLinearForm( bubbleStiffContext );
	  
	  BubbleDampInt * bubbleDamp  = new BubbleDampInt(density,
							  densityforbubble,
							  sonicVel, viscosity,
							  bubbleDensity,
							  frequency, isaxi_);
	  BiLinFormContext * bubbleDampContext = 
	    new BiLinFormContext( bubbleDamp, DAMPING );

	  bubbleDampContext->SetResults( results_[0], results_[0],
					 actSDList, actSDList );
	  bubbleDampContext->SetPtPdes(this, this);

	  assemble_->AddBiLinearForm( bubbleDampContext );
	  
	  //remember the forms objects to associate them with the result
	  //vectors of the bubble dynamics
	  bubbleDampIntMap_[actRegion] = bubbleDamp;
	  bubbleStiffIntMap_[actRegion] = bubbleStiff;
        }
        
        
        
        // Finally add the stiffness/mass integrators
        assemble_->AddBiLinearForm( stiffContext );
        assemble_->AddBiLinearForm( massContext );
      }
      
      // Give result to equation numbering class
      eqnMap_->AddResult( *results_[0], actSDList );
    }
    
    // **********************************************************************
    //   inhom. Neumann boundary condition
    // **********************************************************************
    for( UInt iBc = 0; iBc < inBcs_.GetSize(); iBc++ ) {
      
      // get current Bc
      InhomNeumannBc const & actBc = *inBcs_[iBc];
      
      //BaseForm *neumannBC = new VolumeSrcInt( amplitude, isaxi_ );
      LinearSurfForm *neumannBC = new LinNeumannInt( actBc.value, actBc.phase,
                                                     DENSITY, isaxi_ );
      neumannBC->SetVoluInfo( materials_ );
      LinearFormContext * neumannContext = 
        new LinearFormContext( neumannBC );
      neumannContext->SetPtPde( this );
      neumannContext->SetResult( actBc.result, actBc.entities );
      assemble_->AddLinearForm( neumannContext );
      
      // Give result to equation numbering class
      eqnMap_->AddResult( *actBc.result, actBc.entities );
    }

    // **********************************************************************
    //   surface-integration: Absorbing boundaries
    // **********************************************************************
    if ( absorbingBCs_ == true) { // && analysistype_ != HARMONIC ) {
      for (UInt actSD = 0; actSD < absBCs_.GetSize(); actSD++) {

        // create new entity list
        shared_ptr<SurfElemList> actSDList( new SurfElemList(ptgrid_ ) );
        actSDList->SetRegion( absBCs_[actSD] );
        
        AbsorbingBCsInt * bilinear_damp = new AbsorbingBCsInt(isaxi_);
        bilinear_damp->SetFirstVoluInfo(pdename_, materials_ );

        // In the case of acou-mech coupling we have to multiply the 
        // abc-Integrator matrix with -1
        if ( isMechCoupled_ == true ) {
          bilinear_damp->SetFactor("-1.0");
        }
        // In the case of acou-nrbc coupling we have to multiply the 
        // abc-Integrator matrix with C0 and multiply it by nrbcBeta0
        if ( isNrbcCoupled_ == true ) {
          bilinear_damp->SetFactor(GenStr(sqrt( compressibility / density )));
        }
        BiLinFormContext * abcContext = 
          new BiLinFormContext( bilinear_damp, DAMPING );
        abcContext->SetPtPdes(this, this);     
        abcContext->SetResults( results_[0], results_[0],
                                actSDList, actSDList );
        assemble_->AddBiLinearForm( abcContext );

        // Give result to equation numbering class
        eqnMap_->AddResult( *results_[0], actSDList );
      }
    }

    // =======================================================================
    // Integrators for NonConforming Interfaces
    // =======================================================================
    for( UInt i = 0; i < ncIFaces_.GetSize(); i++ ) {
      
      // get regionId of Lagrangian surface
      StdVector<std::string> keyVec, attrVec, valVec;
      std::string slaveSide;
      std::string ncIfaceName = ptgrid_->RegionIdToName(ncIFaces_[i]);

      ParamNode* ncIfaceListNode;
      ncIfaceListNode = param->Get("domain")->Get("ncInterfaceList");
      slaveSide = ncIfaceListNode->Get("ncInterface",
                                       "name",
                                       ncIfaceName)->Get("slaveSide")->AsString();

      // Part 1: Define integrator M(Psi, Lambda) on
      //         non-conforming interface
      LOG_DBG2(acoupde) << "NonMatching: Defining nonconforming integrator"
                        << " for M on interface '"
                        << ptgrid_->RegionIdToName(ncIFaces_[i]) << "'.";
      shared_ptr<ElemList> actNcList( new ElemList(ptgrid_ ) );
      actNcList->SetRegion( ncIFaces_[i] );
      
      NonConformingInt * ncInt = 
        new NonConformingInt( 1, isaxi_ );

      NcBiLinFormContext * stiffIntDescr = 
	new NcBiLinFormContext( ncInt , STIFFNESS );

      // Force assembling of M(Psi, Lambda)^T
      stiffIntDescr->SetCounterPart( true );

      stiffIntDescr->SetPtPdes(this, this);
      stiffIntDescr->SetResults( results_[0], results_[1],
                                 actNcList, actNcList );
      
      assemble_->AddBiLinearForm( stiffIntDescr );


      // Part 2: Define integrator D(Psi, Lambda) on
      //         Lagrangian surface
      LOG_DBG2(acoupde) << "NonMatching: Defining mass integrator"
                        << " for D on interface '"
                        << ptgrid_->RegionIdToName(ncIFaces_[i]) << "'.";
      shared_ptr<SurfElemList> actSDList( new SurfElemList(ptgrid_ ) );      
      actSDList->SetRegion( ptgrid_->RegionNameToId( slaveSide ) );

      // D(Psi, Lambda) has the form of a standard mass
      // integrator with factor 1.0
      MassInt * dMatInt = new MassInt( 1.0, 1, isaxi_ );
      BiLinFormContext * dMatContext = 
        new BiLinFormContext( dMatInt, STIFFNESS );

      // Force assembling of D(Psi, Lambda)^T
      dMatContext->SetCounterPart( true );
      dMatContext->SetPtPdes( this, this );
      dMatContext->SetResults( results_[0], results_[1],
                               actSDList, actSDList );
      
      assemble_->AddBiLinearForm( dMatContext );

      // Give result LAGRANGE_MULT to equation numbering class
      eqnMap_->AddResult( *results_[1], actSDList );
    }
  }

   void AcousticPDE::DefineSolveStep() {
     ENTER_FCN( "AcousticPDE::DefineSolveStep" );

     solveStep_ = new SolveStepAcoustic(*this);
   }



//    // ======================================================
//    // ALGSYS SECTION (SOLVER, ...) need for acoububble coupling
//    // ======================================================
//    void AcousticPDE::DefineAlgSys() {
    
//      ENTER_FCN( "Acoustic::DefineAlgSys" );

//      //=====================================================
//      // Only if acousticPDE is iteratively coupled with
//      // the bubblePDE
//      //====================================================

//      if( isBubbleCoupled_){

//        // Process input couplings of the bubblePDE to the
//        // additional baseforms of the acousticPDE 
//        UInt numInCouplings = ptCoupling_->GetNumInputCouplings();
      
//        // bubble coupling data
//        StdVector<UInt> elemNumbers;
//        Vector<Double> * radius = NULL;
//        Vector<Double> * radiusDeriv  = NULL;
//        std::cerr<<"definealssys numOutCouplings "<<numInCouplings<< std::endl;
//        for (UInt i = 0; i < numInCouplings; i++) {
      
//  	std::cerr<<"definealgsys ptCoupling_->GetInputQuantity(i)"<<ptCoupling_->GetInputQuantity(i)<< std::endl;
//  	if (ptCoupling_->GetInputQuantity(i) == BUBBLE_RADIUS) {
	
//  	  // get coupling elements and create vector with element numbers
//  	  StdVector<Elem*> * elems = NULL;
//  	  ptCoupling_->GetInputElements(i, elems );
 	
//  	  elemNumbers.Resize( elems->GetSize() );
//  	  for (UInt iElem = 0; iElem < elems->GetSize(); iElem++ ) {
// 	    elemNumbers[iElem] = (*elems)[iElem]->elemNum;
//  	  }
	
//  	  // get buffer with coupling values 
//  	  CFSVector * temp;
//  	  ptCoupling_->GetInputValues( i, temp );    
//  	  std::cerr<<"acou definealgsys size temp"<< (* temp).GetSize()<<std::endl;

//  	  radius = dynamic_cast<Vector<Double>*>(temp);

//  	} else if (ptCoupling_->GetInputQuantity(i) == BUBBLE_RADIUS_DERIV_1 ) {
	
//  	  // get buffer with coupling values 
//  	  CFSVector * temp;
//  	  ptCoupling_->GetInputValues( i, temp );
//  	  radiusDeriv = dynamic_cast<Vector<Double>*>(temp);
	
//  	}
//  	std::cerr<<"test in acou definealgsys nach"<< std::endl;
//        }

//        // Iterate over all BubbleMass/Stiff-Integrators and pass the
//        // radius data to them

//        std::map<RegionIdType, BubbleDampInt*>::iterator dampIt;
//        std::map<RegionIdType, BubbleStiffInt*>::iterator stiffIt;

//        for( dampIt = bubbleDampIntMap_.begin();
//  	   dampIt != bubbleDampIntMap_.end(); dampIt++ ) {
//  	dampIt->second->SetValues( elemNumbers, radius, radiusDeriv);
//        }

//        for( stiffIt = bubbleStiffIntMap_.begin();
//  	   stiffIt != bubbleStiffIntMap_.end();
//  	   stiffIt++ ) {
//  	stiffIt->second->SetValues( elemNumbers, radius, radiusDeriv);
//        }
//      }
    
//      SinglePDE::DefineAlgSys();

//    }




  // ======================================================
  // TIME STEPPING SECTION
  // ======================================================

  void AcousticPDE::InitTimeStepping() {
    ENTER_FCN( "AcousticPDE::InitTimeStepping" );

    // this includes rayleigh and thermoviscous damping
    if ( fracDamping_ == false ) { 
      if ( effectiveMass_ == true ) {
        if ( diagMass_ == true ) {
	  //explicit time stepping
	  TS_alg_ = new NewmarkEffMass( algsys_, true );
	}
	else {
	  TS_alg_ = new NewmarkEffMass( algsys_ );
	}
      }
      else if ( effectiveMass_ == false ) {
        TS_alg_ = new Newmark( algsys_ );
      }
    }
    else {
      if ( effectiveMass_ == false )
        TS_alg_ = new NewmarkFracDamp( algsys_,  
                                       pdeId_, eqnMap_, ptgrid_, this,
                                       results_[0],
                                       subdoms_, dampingList_ );
      else
        EXCEPTION( "This needs to be implemented!" );
    }
  
    // Needed for nonlinear wave equation, get memory for linear part of RHS
    if ( nonLin_ ) {
      RhsLinVal_.Resize(numPDENodes_);
      RhsLinVal_.Init();
    }

  }



  // =========================================================================
  // COUPLING SECTION
  // =========================================================================

  void AcousticPDE::InitCoupling(PDECoupling * Coupling) {
    
    ENTER_FCN( "AcousticPDE::InitCoupling" );
    
    isIterCoupled_ = true;
    ptCoupling_   = Coupling;


    // Process all otput couplings
    UInt numOutCouplings = ptCoupling_->GetNumOutputCouplings();

    StdVector<StdVector<UInt> > elemNodeToCouplingNode_tmp;
    elemNodeToCouplingNode_.Resize(numOutCouplings);
    elemNodeToCouplingNode_.Init();
    
    for (UInt i = 0; i < numOutCouplings; i++) {

      if (ptCoupling_->GetOutputQuantity(i) == ACOU_FORCE) {
        // Intialize the memory of the coupling values
        ptCoupling_->CreateCouplingVector(i,isComplex_);

        // now since we need a incremental formulation, 
        //  initialize some necessary vectors
        isIncrFormulation_ = true;
      }

      else if (ptCoupling_->GetOutputQuantity(i) == ACOU_POT_NRBC)    {
        // Intialize the memory of the coupling values
        ptCoupling_->CreateCouplingVector(i,isComplex_);
        
        // now since we need a incremental formulation, 
        //  initialize some necessary vectors
        isIncrFormulation_ = true;
      }

      else if (ptCoupling_->GetOutputQuantity(i) == ACOU_PRESSURE ||
               ptCoupling_->GetOutputQuantity(i) == ACOU_PRESSURE_DERIV_1 ) {
        // Intialize the memory of the coupling values
        ptCoupling_->CreateCouplingVector(i,isComplex_);
        
      }
      
      else if (ptCoupling_->GetOutputQuantity(i) == ACOU_POWERDENSITY) {
        // Intialize the memory of the coupling values
        // DIRTY HACK ALARM!
        // =================
        // In case of transient-transient coupling coupling vector is Double.
        // In case of harmonic-transient coupling we also want a Double vector!
        ptCoupling_->CreateCouplingVector(i,false);

        //get the element-node to coupling node matching
        StdVector<std::string> couplRegions;
        ptCoupling_->GetOutputRegions(i, couplRegions);
        StdVector<RegionIdType> regionIds;
        ptgrid_->RegionNameToId( regionIds, couplRegions );

        //Get total number of coupling elements
        UInt totalCouplingElems = ptgrid_->GetNumElems( regionIds );
        
        elemNodeToCouplingNode_tmp.Resize(totalCouplingElems);
        elemNodeToCouplingNode_tmp.Init();
        
        UInt offset = 0;
        for ( UInt reg = 0; reg < couplRegions.GetSize(); reg++ ) {

          // find subdomain index = SDidx
          Integer SDidx = subdoms_.Find( regionIds[reg] );
          if (SDidx==-1) {
            EXCEPTION( "AcousticPDE: Coupling region is not within the "
                       << "subdomains of the PDE!" );
          }

          // get elements belonging to subdomain
          StdVector<Elem*> elemssd;
          ptgrid_->GetVolElems(elemssd, subdoms_[SDidx]);

          StdVector<UInt> * couplingnodes = NULL;
          ptCoupling_->GetOutputNodes(i, couplingnodes);
          if ( couplingnodes == NULL ) {
            EXCEPTION( "The pointer 'couplingnodes' is NULL!" );
          }

          for (UInt actEl=0; actEl< elemssd.GetSize(); actEl++) {
            StdVector<UInt> & connecth = elemssd[actEl]->connect;
            elemNodeToCouplingNode_tmp[offset+actEl].Resize(connecth.GetSize());
            elemNodeToCouplingNode_tmp[offset+actEl].Init();

            for ( UInt elnode = 0; elnode < connecth.GetSize(); elnode++ ) {
              for (UInt cnode=0; cnode<(*couplingnodes).GetSize(); cnode++) {

                if (connecth[elnode] == (*couplingnodes)[cnode] ) {
                  elemNodeToCouplingNode_tmp[offset+actEl][elnode] = cnode;
                  break;
                }
              }
            }
          }
          //in the case, that we have more than one coupling region!
          offset = elemssd.GetSize();
        }
        elemNodeToCouplingNode_[i]  = elemNodeToCouplingNode_tmp;
      }

    }



//     // Process input couplings
//     UInt numInCouplings = ptCoupling_->GetNumInputCouplings();

//     // bubble coupling data
//     StdVector<UInt> elemNumbers;
//     Vector<Double> * radius = NULL;
//     Vector<Double> * radiusDeriv  = NULL;
//     	  std::cerr<<"acouInitCoupl numOutCouplings "<<numOutCouplings<< std::endl;
//     for (UInt i = 0; i < numInCouplings; i++) {
      
// 	  std::cerr<<"acouInitCoupl ptCoupling_->GetInputQuantity(i)"<<ptCoupling_->GetInputQuantity(i)<< std::endl;
//       if (ptCoupling_->GetInputQuantity(i) == BUBBLE_RADIUS) {
	
//  	// get coupling elements and create vector with element numbers
// 	StdVector<Elem*> * elems = NULL;
//  	ptCoupling_->GetInputElements(i, elems );
 	
//  	elemNumbers.Resize( elems->GetSize() );
//  	for (UInt iElem = 0; iElem < elems->GetSize(); iElem++ ) {
//  	  elemNumbers[iElem] = (*elems)[iElem]->elemNum;
//  	}
	
// 	// get buffer with coupling values 
// 	CFSVector * temp;
// 	ptCoupling_->GetInputValues( i, temp );    
// 	  std::cerr<<"acou initcoupling size temp"<< (* temp).GetSize()<<std::endl;

// 	radius = dynamic_cast<Vector<Double>*>(temp);

//       } else if (ptCoupling_->GetInputQuantity(i) == BUBBLE_RADIUS_DERIV_1 ) {
	
// 	// get buffer with coupling values 
// 	CFSVector * temp;
// 	ptCoupling_->GetInputValues( i, temp );
// 	radiusDeriv = dynamic_cast<Vector<Double>*>(temp);
	
//       }
// 	  std::cerr<<"test in acou initcouplng nach"<< std::endl;
//     }

//     // Iterate over all BubbleMass/Stiff-Integrators and pass the
//     // radius data to them

//     std::map<RegionIdType, BubbleDampInt*>::iterator dampIt;
//     std::map<RegionIdType, BubbleStiffInt*>::iterator stiffIt;

//     for( dampIt = bubbleDampIntMap_.begin();
// 	 dampIt != bubbleDampIntMap_.end(); dampIt++ ) {
//       dampIt->second->SetValues( elemNumbers, radius, radiusDeriv);
//     }

//     for( stiffIt = bubbleStiffIntMap_.begin();
// 	 stiffIt != bubbleStiffIntMap_.end();
// 	 stiffIt++ ) {
//       stiffIt->second->SetValues( elemNumbers, radius, radiusDeriv);

//     }
      
  }
  

  void AcousticPDE::CalcOutputCoupling() {

    ENTER_FCN( "AcousticPDE::CalcOutputCoupling" );

    UInt dof;
    SolutionType quantity;
    StdVector<Elem*> * couplingElems = NULL;
    StdVector<UInt> * couplingNodes = NULL;
    CFSVector * temp_values = NULL;
    UInt regionCount = 0;
  
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
        if (quantity == ACOU_FORCE) {
          
          ptCoupling_->GetOutputElements(i, couplingElems);
          dof = ptCoupling_->GetOutputDof(i);
          CalcMechCouplingRHS(couplingElems, *couplingNodes, *values, dof);
        }
        else if (quantity == ACOU_POWERDENSITY) {
          if ( isComplex_ == false ) {
            CalcHeatCouplingRHS<Double>( *values, 
                                         elemNodeToCouplingNode_[regionCount],
                                         i, couplingNodes->GetSize() );
          } 
          else {
            CalcHeatCouplingRHS<Complex>(*values, 
                                         elemNodeToCouplingNode_[regionCount],
                                         i, couplingNodes->GetSize() );
          }        
          regionCount++;
        }
        else if (quantity == ACOU_POT_NRBC) {
          ptCoupling_->GetOutputElements(i, couplingElems);
          dof = ptCoupling_->GetOutputDof(i);
          
          // Here this call gives the values to phi0 in nrbcPDE
          CalcNRBCCouplingRHS(couplingElems, *couplingNodes,
                              *values, dof);                              
        }
        break;

      case ELEM:
        ptCoupling_->GetOutputElements(i, couplingElems);
        if ( quantity == ACOU_PRESSURE ||
             quantity == ACOU_PRESSURE_DERIV_1 ) {
          CalcBubblePressure( *couplingElems, *values, quantity );
        }
        break;
      }
    }
  }


  void AcousticPDE::
  CalcMechCouplingRHS( StdVector<Elem*> * couplingElems, 
                       StdVector<UInt> & couplingNodes,
                       Vector<Double>& elemCouplingSols,
                       UInt couplingdof ) {
    
    ENTER_FCN( "AcousticPDE::CalcMechCouplingRHS" );
    
    EXCEPTION( "Not working at the moment" );
    
    //     Double density = 0.0;
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
      
    //       // If first volume element does not belong to acoustic PDE, try the
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
    //         (*error) << "AcousticPDE::CalcMechCouplingRHS: The two volume "
    //                  << "element neighbours of surface element Nr. "
    //                  << actCoupleElem->elemNum << " do not belong to my regions!";
    //         Error( __FILE__, __LINE__ );
    //       }
      
    //       // Assign correct density
    //       materials_[subdoms_[matIndex]]->GetScalar(density,DENSITY,REAL);
      
    //       BaseForm * bilinear_mass = new MassInt(ptElem, density, isaxi_);
    //       bilinear_mass->CalcElementMatrix(ptCoord, elemMat);
    //       delete bilinear_mass;     
      
    //       GetDerivSolVecOfElement(sol, connecth);
      
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
  
  void AcousticPDE::
  CalcNRBCCouplingRHS( StdVector<Elem*> * couplingElems, 
                       StdVector<UInt> & couplingNodes,
                       Vector<Double>& elemCouplingSols,
                       UInt couplingdof ) {
    
    ENTER_FCN( "AcousticPDE::CalcNRBCCouplingRHS" );
    
    EXCEPTION( "Not working at the moment" );

    //     Double density = 0.0;
    //     Double sign = 0.0;
    //     Double coeff_Pmass, coeff_Qstiff;
    //     UInt nrbcMatType=0; // 1=Rmat, 2=pmat, 3=Qmat
    //     Integer matIndex = -1;
    //     Elem * ptVolElem = NULL;
    //     Matrix<Double> ptCoord, elemMat;
    //     Vector<Double> sol, deriv2sol, Qmat_x_sol, Pmat_x_2derSol, normal;
    
    //     elemCouplingSols.Init(0.0);
    
    //     for (UInt actElem=0; actElem<couplingElems->GetSize(); actElem++) {
      
    //       // Perform cast from volume element to surface element, since
    //       // nrbc-acou coupling makes only sense on surface elements
    //       SurfElem * actSurfCoupleElem = 
    //         dynamic_cast<SurfElem*> ((*couplingElems)[actElem]);

    //       if (actSurfCoupleElem == NULL) {
    //         Error( "No elements found for coupling!", __FILE__, __LINE__ );
    //       }
      
    //       BaseFE * ptElem = actSurfCoupleElem->ptElem;
    //       StdVector<UInt> & connecth = actSurfCoupleElem->connect;
    //       GetElemCoords(connecth, ptCoord);
      
    //       // Try to find according region for first neighbouring volume
    //       // element of the surface element
    //       matIndex = subdoms_.Find(actSurfCoupleElem->ptVolElem1->regionId);
      
    //       // If first volume element does not belong to acoustic PDE, try the
    //       // second one
    //       if ( matIndex == -1 ) {
    //         matIndex = subdoms_.Find(actSurfCoupleElem->ptVolElem2->regionId);
    //         ptVolElem = actSurfCoupleElem->ptVolElem2;
    //         //        sign = actSurfCoupleElem->normalSign;
    //       } else {
    //         ptVolElem = actSurfCoupleElem->ptVolElem1;
    //         //        sign = -1.0 * actSurfCoupleElem->normalSign;
    //       }
      
    //       if ( matIndex == -1) {
    //         (*error) << "AcousticPDE::CalcNRBCCouplingRHS: The two volume "
    //                  << "element neighbours of surface element Nr. "
    //                  << actSurfCoupleElem->elemNum << " do not belong to my regions!"; //         Error( __FILE__, __LINE__ );
    //       }
      
    //       // Density set to 1 since for this mass integrator no fac is needed
    //       density = 1.0;
      
    //       // !!!!!!!!!!!!!!!USE PTELEM!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //       // FIRST WE COMPUTE THE "Qmat_x_sol" RHS term
    //       // Standard "mass" matrix shape so using default constructor
    //       coeff_Qstiff=1.0;
    //       nrbcMatType=3;
    //       NrbcInt * bilinear_Qstiff  = 
    //         new NrbcInt(ptElem, coeff_Qstiff, nrbcMatType, isaxi_);
      
    //       // Now we get the Q matrix with tangential derivatives 
    //       bilinear_Qstiff->CalcElementMatrix(ptCoord, elemMat);
      
    //       //std::cout<<"%%%%%%----ACOUSTIC-PDE COUPLING OUTPUT----%%%%%"<<std::endl;
    //       // std::cout<<"Qmatrix= "<<std::endl;
    //       // std::cout<<elemMat<<std::endl;
      
    //       delete  bilinear_Qstiff;
    //       GetSolVecOfElement(sol, connecth);
    //       //          for(UInt k=0; k<sol.GetSize(); k++)
    //       //           if (abs(sol[k])<=5e-16)
    //       //             sol[k]=0;
        
    //       // std::cout<<"AcouSolution= "<<std::endl;
    //       // std::cout<<sol<<std::endl;
    //       Qmat_x_sol = elemMat * sol;
      
    //       // std::cout<<"Qmat_x_sol= "<<std::endl;
    //       // std::cout<<Qmat_x_sol<<std::endl;
    //       //  Q_x_sol has to be added on RHS of nrbcPDE1 with negative sign
      
    //       Qmat_x_sol *= - 1.0;
    //       //END OF "Qmat_x_sol" Computation   
      
    //       //NOW WE COMPUTE "Pmat_x_Deriv2Sol" RHS term 
    //       coeff_Pmass=1.0;
    //       nrbcMatType=2;
    //       NrbcInt * bilinear_Pmass  = 
    //         new NrbcInt(ptElem, coeff_Pmass, nrbcMatType, isaxi_);
      
    //       bilinear_Pmass->CalcElementMatrix(ptCoord, elemMat);
        
    //       //std::cout<<"Pmatrix= "<<std::endl;
    //       // std::cout<<elemMat<<std::endl;
      
    //       delete  bilinear_Pmass;
      
    //       GetDeriv2SolVecOfElement(deriv2sol, connecth);
      
    //       // std::cout<<"Acou2ndDerivSolution= "<<std::endl;
    //       // std::cout<<deriv2sol<<std::endl;
      
    //       Pmat_x_2derSol = elemMat * deriv2sol;
    //       // Computation of alpha factor (like in nrbcPDE)
    //       Double c0, alphaNRBC, Cj, density, compressibility;
    //       //actSD is 0 since there is only one acoustic domain
    //       materials_.begin()->second->GetScalar( density, DENSITY, REAL );
    //       density=1;
    //       materials_.begin()->second->GetScalar( compressibility, ACOU_BULK_MODULUS, REAL );
    //       c0 = sqrt(compressibility/density);
    //       Cj = 1.0;
    //       alphaNRBC = (1/(Cj*Cj) - 1/(c0*c0));
    //       Pmat_x_2derSol *= alphaNRBC;
      
    //       // std::cout<<"Pmat_x_2derSol= "<<std::endl;
    //       // std::cout<<Pmat_x_2derSol<<std::endl;
      
    //       // END OF "Pmat_x_2derSol" Computation   
      
    //       for (UInt actNode=0; actNode<ptCoord.GetSizeRow(); actNode++) {
    //         UInt nodePos = 0;
    //         while(connecth[actNode] != couplingNodes[nodePos] &&
    //               nodePos < couplingNodes.GetSize()) {
    //           nodePos++;
    //         }
        
    //         elemCouplingSols[nodePos*couplingdof] = 
    //           (Qmat_x_sol[actNode] + Pmat_x_2derSol[actNode]);
        
    //         // std::cout <<"elemCouplingSols["<<(nodePos*couplingdof)<<"]= "
    //         //           <<elemCouplingSols[nodePos*couplingdof]<<std::endl;
        
    //       }
    //       // std::cout<<"elemCouplingSols = "<<elemCouplingSols<<std::endl;
      
    //       //std::cout <<"%%%--END ACOUSTIC-PDE COUPLING OUTPUT INFO----%%%%"
    //       //          <<std::endl;                  
      
    //     }
  }

  
  template <class TYPE>
  void AcousticPDE::
  CalcHeatCouplingRHS(Vector<Double> & sourceValue, 
                      StdVector<StdVector<UInt> > & elemNodeToCouplingNode,
                      UInt actCoupling,
                      UInt numCouplingNodes) {


    ENTER_FCN( "AcousticPDE::CalcHeatCouplingRHS" );
    
    // get the coupling regions
    StdVector<std::string> couplRegions;
    ptCoupling_->GetOutputRegions(actCoupling, couplRegions);
    StdVector<RegionIdType> regionIds;
    ptgrid_->RegionNameToId( regionIds, couplRegions );

    // Operator for calculating energy density
    AcouPowerDensityOp<TYPE> *SourceOp = 
      new AcouPowerDensityOp<TYPE>( ptgrid_, this, eqnMap_, isaxi_ );

    // initialize output vector
    sourceValue.Init();

    Vector<Double> elemPowerDensity;
    
    UInt offset = 0;
    for (UInt reg=0; reg<couplRegions.GetSize(); reg++) {
      
      // find subdomain index
      Integer SDidx = subdoms_.Find( regionIds[reg] );
      Double density;
      materials_[regionIds[reg]]->GetScalar(density,DENSITY,REAL);
      
      // get elements belonging to subdomain
      ElemList actSDList(ptgrid_ );
      actSDList.SetRegion( subdoms_[SDidx] );
      
      EntityIterator it = actSDList.GetIterator();
      UInt actEl = 0;
      for ( it.Begin(); !it.IsEnd(); it++, actEl++) {

        SourceOp->CalcElemPD(elemPowerDensity, it, density);

        // Add the element energy to the according coupling node
        StdVector<UInt> const & connecth = it.GetElem()->connect;
        for (UInt elnode=0; elnode<connecth.GetSize(); elnode++) {

          sourceValue[elemNodeToCouplingNode[actEl+offset][elnode]]
            += elemPowerDensity[elnode];
        }
      }
      
      //in the case, that we have more than one coupling region!
      offset = actSDList.GetSize();
    }
  }
  void AcousticPDE::CalcBubblePressure( StdVector<Elem*> & elems,
                                        Vector<Double>& couplVals,
                                        SolutionType solType ) {

    ENTER_FCN( "AcousticPDE::CalcBubblePressure" );

//     //NUR BEI BOSCH!!!
//     mHandle_ =  domain->GetMathParser()->GetNewHandle();
//     MathParser * mParser =  domain->GetMathParser();
//     mParser->SetExpr( mHandle_, "t" );

    // Initialize coupling values
    couplVals.Init();

    Vector<Double> pressureOut(3);
    pressureOut.Init();

  //  MathParser * mParser =  domain->GetMathParser();

    if ( solType == ACOU_PRESSURE ) {

      Vector<Double> elPressure;
      for (UInt iElem = 0; iElem < elems.GetSize(); iElem++ ) {
        ElemList elemList(ptgrid_);
        elemList.SetElement( elems[iElem] );
        EntityIterator it = elemList.GetIterator();
        sol_->GetElemSolution( elPressure, it );
        Double pressure = 0.0;

        for( UInt iNode = 0; iNode < elPressure.GetSize(); iNode++ ) {
          pressure += elPressure[iNode];
        }
        pressure /= elPressure.GetSize();

        // store back to vector
        couplVals[iElem] = pressure;


//  	if( locElemNum == 2 ){
//  	  pressureOut[0]= mParser->Eval( mHandle_ );
//  	  pressureOut[1]= pressure;
//  	  //	  std::cout<<pressureOut[0]<<"   "<<pressureOut[1];
//  	}
	


      }

    } else {
      Vector<Double> elPressureDeriv;
      for (UInt iElem = 0; iElem < elems.GetSize(); iElem++ ) {
        ElemList elemList(ptgrid_);
        elemList.SetElement( elems[iElem] );
        EntityIterator it = elemList.GetIterator();
        GetDerivSolVecOfElement( elPressureDeriv, it );
        Double pressureDeriv = 0.0;

        for( UInt iNode = 0; iNode < elPressureDeriv.GetSize(); iNode++ ) {
          pressureDeriv += elPressureDeriv[iNode];
        }
        pressureDeriv /= elPressureDeriv.GetSize();

        // store back to vector
        couplVals[iElem] = pressureDeriv;


//  	locElemNum = eqnMap_->Mesh2PdeElem(elems[iElem]->elemNum );
//  	if( locElemNum == 2 ){
//  	  pressureOut[2]= pressureDeriv;
//  	  //	  std::cout<<"    "<<pressureOut[2]<<std::endl;; 

//  	}
	  
      }
    }

  }  


  bool AcousticPDE::HasOutput(SolutionType output) {
    ENTER_FCN( "AcousticPDE::HasOutput" );
    if ((output == ACOU_FORCE) || (output == ACOU_POWERDENSITY)) {
      return true;
    }
    if (output == ACOU_POT_NRBC) {
      return true;
    }
    
    if ( formulation_ == ACOU_PRESSURE ) {
      if ( output == ACOU_PRESSURE ) {
        return true;
      }
      if ( output == ACOU_PRESSURE_DERIV_1 ) {
        return true;
      }
    }
    return false;
  }



  // =========================================================================
  // POSTPROCESSING SECTION
  // =========================================================================

  

  void AcousticPDE::CalcResults( shared_ptr<BaseResult> result ) {
    ENTER_FCN( "AcousticPDE::CalcResults" );
    
    switch (result->GetResultInfo()->resultType ) {
      
    case ACOU_POTENTIAL:
      if( formulation_ == ACOU_POTENTIAL ) {
        if( isComplex_ ) {
          ExtractResult<Complex>( result, sol_ );
        } else {
          ExtractResult<Double>( result, sol_ );
        }
      } else {
        EXCEPTION( "Acoupotential only available for potential formulation!" );
      }
      break;

    case ACOU_PRESSURE:
      if( formulation_ == ACOU_POTENTIAL ) {
        if( isComplex_ ) {
          CalcElemPressure<Complex>( result );
        } else {
          CalcElemPressure<Double>( result );
        }
      } else {
        if( isComplex_ ) {
          ExtractResult<Complex>( result, sol_ );
        } else {
          ExtractResult<Double>( result, sol_ );
        }
      }
      break;
      
    case ACOU_POTENTIAL_DERIV_1:
      ExtractDerivResult( result, 1 );
      break;

    case ACOU_POTENTIAL_DERIV_2:
      ExtractDerivResult( result, 2 );
      break;

    case ACOU_PRESSURE_DERIV_1:
      ExtractDerivResult( result, 1 );
      break;

    case ACOU_PRESSURE_DERIV_2:
      ExtractDerivResult( result, 2 );
      break;

    case ACOU_RHSVAL:
      Warning( "Not yet implemented!");
      break;

    case ACOU_INTENSITY:
      if( isComplex_ ) {
        CalcAcouIntensity<Complex>( result );
      } else {
        EXCEPTION( "Acoustic intensity only computable for harmonic "
                   << "analysis!" );
      }
      break;

    case ACOU_POWER:
      if( isComplex_ ) {
        CalcAcouPower<Complex>( result );
      } else {
        EXCEPTION( "Acoustic power only computable for harmonic "
                   << "analysis!" );
      }
      break;
    case ACOU_FORCE:
      if( isComplex_ ) {
        CalcForce<Complex>( result );
      } else {
        CalcForce<Double>( result );
      }
      break;
      
    default:
      Warning( "Resulttype not computable by acoustic PDE",
               __FILE__, __LINE__ );
    }
   }
  
 
  template <class TYPE>
  void AcousticPDE::CalcForce( shared_ptr<BaseResult> vals ) {  
    ENTER_FCN( "AcousticPDE::CalcForce" );
    
    Matrix<Double> ptCoord;

    // get data from result object and resize its vector
    Result<TYPE> &  actRes = 
      dynamic_cast<Result<TYPE>&>(*vals);
    EntityIterator it = actRes.GetEntityList()->GetIterator();
    
    Vector<TYPE> & actVal = actRes.GetVector();
    actVal.Resize( actRes.GetEntityList()->GetSize()  );
    
    // loop over elements
    for (it.Begin(); !it.IsEnd(); it++ ) {

      // Perform cast from volume element to surface element,
      //  since calculation of force makes only sense on a surface
      const SurfElem * actSaveElem = it.GetSurfElem();
      BaseFE * ptElem = actSaveElem->ptElem;
      StdVector<UInt> connect = actSaveElem->connect;
      ptgrid_->GetElemNodesCoord( ptCoord, connect, false );
              
      Vector<Double> valueElem;
      if ( formulation_ == ACOU_POTENTIAL ) {

        Integer matIndex = -1;
        // Try to find according region for first neighbouring volume
        // element of the surface element
        matIndex = subdoms_.Find(actSaveElem->ptVolElem1->regionId);

        // If first volume element does not belong to acoustic PDE, try the
        // second one
        Elem * ptVolElem = NULL;
        if ( matIndex == -1 ) {
          matIndex = subdoms_.Find(actSaveElem->ptVolElem2->regionId);
          ptVolElem = actSaveElem->ptVolElem2;
        }
        else {
          ptVolElem = actSaveElem->ptVolElem1;
        }
      	
        if ( matIndex == -1) {
          EXCEPTION( "AcousticPDE::CalcForce: The two volume element"
                     << "neighbours of surface element no. "
                     << actSaveElem->elemNum << " don't belong to my region!" );
        }
    
        // Assign correct density
        Double density;
        materials_[subdoms_[matIndex]]->GetScalar(density,DENSITY,REAL);
      	
        // retrieve 1st derivative, since F = rho * dpsi/dt * A
        //ElemList elemList(ptgrid_);
        //elemList.SetElement( actSaveElem );
        //EntityIterator it = elemList.GetIterator();
        GetDerivSolVecOfElement(valueElem,  it);

        valueElem *= density;
      }
      else if ( formulation_ == ACOU_PRESSURE ) {

        // retrieve solution, since F = p * A
        //ElemList elemList(ptgrid_);
        //elemList.SetElement( actSaveElem );
        //EntityIterator it = elemList.GetIterator();
        GetSolVecOfElement(valueElem,it);
      }

      const UInt nrIntPts= ptElem->GetNumIntPoints();
      const Vector<Double> & intWeights = ptElem->GetIntWeights();

      Double forceElem = 0.0;
      Double jacDet;
      for (UInt actIntPt=1; actIntPt<=nrIntPts;  actIntPt++) {

        jacDet = ptElem->CalcJacobianDetAtIp(actIntPt, ptCoord, actSaveElem);
        Vector<Double> shapeFnc;  
        ptElem -> GetShFncAtIp(shapeFnc, actIntPt, actSaveElem);

        if (isaxi_) {
          Vector<Double> coordAtIP;
          coordAtIP = ptCoord * shapeFnc;
          forceElem +=  (intWeights[actIntPt-1] * jacDet
                         * 2 * PI) * coordAtIP[0] * (valueElem * shapeFnc);
        }
        else {
          forceElem +=  intWeights[actIntPt-1] * jacDet 
            * (shapeFnc * valueElem);
        }
      }

      actVal[it.GetPos()] = forceElem;
    }

    
    TYPE sumAcouForce = 0;
    for(UInt k=0; k<actVal.GetSize(); k++ ) {
      sumAcouForce += actVal[k];
    }
  }

  template <class TYPE>
  void AcousticPDE::CalcElemPressure( shared_ptr<BaseResult> vals) {  
    ENTER_FCN( "AcousticPDE::CalcElemPressure" );

    // get data from result object and resize its vector
    Result<TYPE> &  actRes = 
      dynamic_cast<Result<TYPE>&>(*vals);
    EntityIterator it = actRes.GetEntityList()->GetIterator();
    
    Vector<TYPE> & actVal = actRes.GetVector();
    actVal.Resize( actRes.GetEntityList()->GetSize()  );

    // Density of element in region
    Double density = 0.0;
    TYPE elemPressure = 0.0;
    
    // loop over all elements of subdomain
    for (it.Begin(); !it.IsEnd(); it++ ) {
      BaseFE * ptElem = it.GetElem()->ptElem;
      RegionIdType actRegion = it.GetElem()->regionId;
      materials_[actRegion]->GetScalar(density,DENSITY,REAL);
      
      // get element coordinates
      StdVector<UInt> const & connect = it.GetElem()->connect;
      Matrix<Double> ptCoord;
      ptgrid_->GetElemNodesCoord( ptCoord, connect, false );
      
      // get shape function at center of the element
      Vector<Double> shapeFnc;
      Vector<Double> LCoord;
      ptElem -> GetCoordMidPoint(LCoord);
      ptElem -> GetShFnc(shapeFnc,LCoord,it.GetElem());
      
      // retrieve 1st derivative and multiply with density, 
      //  since p = rho * dpsi/dt
      Vector<TYPE> valueElem;
      GetDerivSolVecOfElement(valueElem, it);
      
      elemPressure = valueElem * shapeFnc * density;
      actVal[it.GetPos()] = elemPressure;
    }
  }
  

  template <class TYPE>
  void AcousticPDE::CalcAcouIntensity( shared_ptr<BaseResult> vals ) {
    ENTER_FCN( "AcousticPDE::CalcAcouIntensity" );

    // currently we just support harmonic analysis: complex data

    //get frequency
    MathParser * parser = domain->GetMathParser();
    parser->SetExpr( mHandle_, "f" );
    Double actFreq = parser->Eval( mHandle_ );

    //check solution type and compute factor
    SolutionType solType;
    Complex multVal = 0;

    // factor 0.5 is due to the fact, that the values are peak values
    if ( formulation_ == ACOU_PRESSURE ) {
      solType = ACOU_PRESSURE;
      multVal = Complex(0.0, -0.5/ (2.0*PI*actFreq) );
    } else {
      solType = ACOU_POTENTIAL;
      multVal = Complex(0.0, -0.5*(2.0*PI*actFreq));
    }

    NodeStoreSol<TYPE> * solhelp = dynamic_cast<NodeStoreSol<TYPE>*>(sol_);

    //some help variables
    Vector<Double> lCoordSurf, lCoordVol, normal;
    Vector<TYPE> gradVal(dim_);
    Vector<TYPE> elemIntensity(dim_);

    Elem * ptVolElem;
    BaseFE * ptSurfElemFE, * ptVolElemFE;

    TYPE gradNormal = 0.0;
    Double normSign = 0;
    Double density  = 0.0;

    // Create vector with interpolation coordinate.
    // For simplicity we only evaluate the integral
    // in coordinate origin
    lCoordSurf.Resize(dim_-1);
    lCoordSurf.Init(0);

    // Create operator for gradient computation of solution
    GradientFieldOp<TYPE> * gradOp = 
      new GradientFieldOp<TYPE>(ptgrid_, this, eqnMap_, *solhelp,
                                solType, results_[0], isaxi_);

    TYPE factorI;  
    Result<TYPE> &  actRes = 
      dynamic_cast<Result<TYPE>&>(*vals);
    EntityIterator it = actRes.GetEntityList()->GetIterator();
    
    Vector<TYPE> & actVal = actRes.GetVector();
    //std::cerr << "size of entityList is " << 
    actVal.Resize( actRes.GetEntityList()->GetSize() * dim_ );
    
    // loop over all surface elements
    UInt counterElems = 0;
    for ( it.Begin(); !it.IsEnd(); it++, counterElems++ ) {
      
      const SurfElem * actSurfElem = it.GetSurfElem();
      // Determine, which volume element is the right neighbour for the 
      // calculation;
      // our normal should point out of the correct neighbor volume element!
      if ( surfNeighborRegions_[vals] ==
           actSurfElem->ptVolElem1->regionId ) {
        ptVolElem = actSurfElem->ptVolElem1;
        normSign = 1.0;
      } 
      else {
        ptVolElem = actSurfElem->ptVolElem2;
        normSign = -1.0;
      }
      
      normSign *= (Double) actSurfElem->normalSign;
      
      ptSurfElemFE = actSurfElem->ptElem; 
      ptVolElemFE = ptVolElem->ptElem;
      
      const StdVector<UInt> & surfConnect = actSurfElem->connect;
      const StdVector<UInt> & volConnect = ptVolElem->connect;
      
      // calculate volume integration coordinates from
      // surface integration coordinat for evalauting the 
      // gradient on the surface of the volume element
      ptSurfElemFE->GetCoordMidPoint(lCoordSurf);
      ptVolElemFE->GetLocalIntPoints4Surface(surfConnect, volConnect,
                                             lCoordSurf, lCoordVol);
      
      BaseMaterial * myMat = materials_[ptVolElem->regionId];
      myMat->GetScalar(density,DENSITY,REAL);
      
      Matrix<Double> CornerCoords; 
      ptgrid_->GetElemNodesCoord( CornerCoords, surfConnect, false );
      
      // Calc gradient
      ElemList elList(ptgrid_);
      elList.SetElement( ptVolElem );
      EntityIterator it2 = elList.GetIterator();
      gradOp->CalcElemGradField(gradVal, it2, 
                                lCoordVol,1.0);
      // Calc global normal
      ptgrid_->CalcSurfNormal(normal, *actSurfElem);
      
      normal    *= normSign;
      gradNormal = normal * gradVal;
      
      //get average solution
      TYPE elemSol = 0;
      TYPE nodeSol;
      for ( UInt k=0; k<surfConnect.GetSize(); k++) {
        solhelp->Get(solType, surfConnect[k]-1,0, nodeSol);
        elemSol += nodeSol;
      }
      elemSol /= (Double)surfConnect.GetSize(); 
      
      // get the conjugate complex value
      elemSol = std::conj(elemSol);
      
      if ( formulation_ == ACOU_PRESSURE ) {
        factorI   = multVal * elemSol / density;
        elemIntensity = gradVal * factorI;
      } else {
        factorI   = multVal * density * elemSol;
        elemIntensity = gradVal * factorI;
      }
      //std::cerr << "elemIntensity = " << elemIntensity.Serialize();
      // loop over dofs
      for(UInt iDim = 0; iDim < dim_; iDim++ ) {
        actVal[it.GetPos()*dim_ + iDim] = elemIntensity[iDim];
      }
    }
    
    delete gradOp;
  }
  
  template<class TYPE>
  void AcousticPDE::CalcAcouPower( shared_ptr<BaseResult> vals ) {
    ENTER_FCN( "AcousticPDE::CalcAcouPower" );
    
    // currently we just support harmonic analysis: complex data

    //get frequency
    MathParser * parser = domain->GetMathParser();
    parser->SetExpr( mHandle_, "f" );
    Double actFreq = parser->Eval( mHandle_ );

    //check solution type and compute factor
    SolutionType solType;
    Complex multVal = 0;

    // factor 0.5 is due to the fact, that the values are peak values
    if ( formulation_ == ACOU_PRESSURE ) {
      solType = ACOU_PRESSURE;
      multVal = Complex(0.0, -0.5/ (2.0*PI*actFreq) );
    }
    else {
      solType = ACOU_POTENTIAL;
      multVal = Complex(0.0, -0.5*(2.0*PI*actFreq));
    }

    NodeStoreSol<TYPE> * solhelp = dynamic_cast<NodeStoreSol<TYPE>*>(sol_);
    

    //some help variables
    Vector<Double> lCoordSurf, lCoordVol, normal;
    Vector<TYPE> gradVal(dim_);
    Vector<TYPE> elemIntensity(dim_);

    
    Elem * ptVolElem;
    BaseFE * ptSurfElemFE, * ptVolElemFE;

    TYPE gradNormal = 0.0;
    TYPE elemPower  = 0.0;
    Double normSign = 0;
    Double density  = 0.0;

    // Create vector with interpolation coordinate.
    // For simplicity we only evaluate the integral
    // in coordinate origin
    lCoordSurf.Resize(dim_-1);
    lCoordSurf.Init(0);

    // Create operator for gradient computation of solution
    GradientFieldOp<TYPE> * gradOp = 
      new GradientFieldOp<TYPE>(ptgrid_, this, eqnMap_, *solhelp,
                                solType, results_[0], isaxi_);
    
    // convert result object
    Result<TYPE> &  actRes = 
      dynamic_cast<Result<TYPE>&>(*vals);      
    EntityIterator regionIt = actRes.GetEntityList()->GetIterator();

    // resize vector
    Vector<TYPE> & actVal = actRes.GetVector();
    actVal.Resize( actRes.GetEntityList()->GetSize() );    
    actVal.Init();

    // Loop over regions
    for( regionIt.Begin(); !regionIt.IsEnd(); regionIt++ ) {
      SurfElemList actSDList(ptgrid_ );
      actSDList.SetRegion( regionIt.GetRegion() );
      EntityIterator it = actSDList.GetIterator();
      
      // Loop over all surface elements
      UInt counterElems = 0;
      for ( it.Begin(); !it.IsEnd(); it++, counterElems++ ) {
   
        const SurfElem * actSurfElem = it.GetSurfElem();

	// Determine, which volume element is the right neighbour for the 
	// calculation;
	// our normal should point out of the correct neighbor volume element!
	if (   surfNeighborRegions_[vals] ==
              actSurfElem->ptVolElem1->regionId ) {
          ptVolElem = actSurfElem->ptVolElem1;
	  normSign = 1.0;
	} 
	else {
	  ptVolElem = actSurfElem->ptVolElem2;
	  normSign = -1.0;
	}
     
	normSign *= (Double) actSurfElem->normalSign;
     
	ptSurfElemFE = actSurfElem->ptElem; 
	ptVolElemFE = ptVolElem->ptElem;
     
	const StdVector<UInt> & surfConnect = actSurfElem->connect;
	const StdVector<UInt> & volConnect = ptVolElem->connect;
     
        // get material information
        BaseMaterial * myMat = materials_[ptVolElem->regionId];
	myMat->GetScalar(density,DENSITY,REAL);

        // Pass ansatz function to surface and volume element
        ptSurfElemFE->SetAnsatzFct( results_[0]->fctType );
        ptVolElemFE->SetAnsatzFct(  results_[0]->fctType );
        
        // Get weights and points of surface element
        UInt nrIntPts= ptSurfElemFE->GetNumIntPoints();
        const Vector<Double> & intWeights = ptSurfElemFE->GetIntWeights(); 
        const Vector<Double> * intPoints = ptSurfElemFE->GetIntPoints(); 
        
        // get global coordintes of surface element
        Matrix<Double> CornerCoords; 
        ptgrid_->GetElemNodesCoord( CornerCoords, surfConnect, false );

        // Loop over integration points
        elemPower = 0.0;
        TYPE helpVal = 0.0;
        for( UInt iPt = 0; iPt < nrIntPts; iPt++ ) {
          
          // calculate volume integration coordinates from
          // surface integration coordinat for evalauting the 
          // gradient on the surface of the volume element
          ptVolElemFE->GetLocalIntPoints4Surface(surfConnect, volConnect,
                                                 intPoints[iPt], lCoordVol);

          // Calculate jacobian gradient of surface element
          Double jacDet = ptSurfElemFE->CalcJacobianDetAtIp(iPt+1, 
                                                            CornerCoords, 
                                                            actSurfElem );

          // Calc gradient
          ElemList elList(ptgrid_);
          elList.SetElement( ptVolElem );
          EntityIterator it2 = elList.GetIterator();
          gradOp->CalcElemGradField(gradVal, it2, 
                                    lCoordVol,1.0);
          // Calc global normal
          ptgrid_->CalcSurfNormal(normal, *actSurfElem);
          
          normal    *= normSign;
          gradNormal = normal * gradVal;

          // get solution of element and interpolate into integration point
          Vector<TYPE> elemSol;
          Vector<Double> shapeFnc;
          GetSolVecOfElement( elemSol, it2);
          ptVolElemFE->GetShFnc( shapeFnc, lCoordVol, ptVolElem );
          TYPE intPointSol = shapeFnc * elemSol;

          // get the conjugate complex value
          intPointSol = std::conj(intPointSol);
          
          if ( formulation_ == ACOU_PRESSURE ) {
            helpVal =  multVal * gradNormal * intPointSol * (1.0 / density );
          }
          else {
            helpVal = multVal * density  * gradNormal * intPointSol;
          }
          
          // integrate value
          elemPower += intWeights[iPt] * helpVal * jacDet; 
        }
        actVal[regionIt.GetPos()] += elemPower;
        
      }
    }

    delete gradOp;
  }





  /**  template <class TYPE>
  void AcousticPDE::CalcAcouPower( shared_ptr<BaseResult> vals ) {
    ENTER_FCN( "AcousticPDE::CalcAcouPower" );
    
    // currently we just support harmonic analysis: complex data

    //get frequency
    MathParser * parser = domain->GetMathParser();
    parser->SetExpr( mHandle_, "f" );
    Double actFreq = parser->Eval( mHandle_ );

    //check solution type and compute factor
    SolutionType solType;
    Complex multVal = 0;

    // factor 0.5 is due to the fact, that the values are peak values
    if ( formulation_ == ACOU_PRESSURE ) {
      solType = ACOU_PRESSURE;
      multVal = Complex(0.0, -0.5/ (2.0*PI*actFreq) );
    }
    else {
      solType = ACOU_POTENTIAL;
      multVal = Complex(0.0, -0.5*(2.0*PI*actFreq));
    }

    NodeStoreSol<TYPE> * solhelp = dynamic_cast<NodeStoreSol<TYPE>*>(sol_);
    

    //some help variables
    Vector<Double> lCoordSurf, lCoordVol, normal;
    Vector<TYPE> gradVal(dim_);
    Vector<TYPE> elemIntensity(dim_);

    
    Elem * ptVolElem;
    BaseFE * ptSurfElemFE, * ptVolElemFE;

    TYPE gradNormal = 0.0;
    TYPE elemPower  = 0.0;
    Double normSign = 0;
    Double density  = 0.0;
    UInt pdeElemNum = 0;

    // Create vector with interpolation coordinate.
    // For simplicity we only evaluate the integral
    // in coordinate origin
    lCoordSurf.Resize(dim_-1);
    lCoordSurf.Init(0);

    // Create operator for gradient computation of solution
    GradientFieldOp<TYPE> * gradOp = 
      new GradientFieldOp<TYPE>(ptgrid_, this, eqnMap_, *solhelp,
                                solType, results_[0], isaxi_);
    
    // convert result object
    Result<TYPE> &  actRes = 
      dynamic_cast<Result<TYPE>&>(*vals);      
    EntityIterator regionIt = actRes.GetEntityList()->GetIterator();

    // resize vector
    Vector<TYPE> & actVal = actRes.GetVector();
    actVal.Resize( actRes.GetEntityList()->GetSize() );    

    // Loop over regions
    for( regionIt.Begin(); !regionIt.IsEnd(); regionIt++ ) {
      SurfElemList actSDList(ptgrid_ );
      actSDList.SetRegion( regionIt.GetRegion() );
      EntityIterator it = actSDList.GetIterator();
      
      // Loop over all surface elements
      UInt counterElems = 0;
      for ( it.Begin(); !it.IsEnd(); it++, counterElems++ ) {
   
        const SurfElem * actSurfElem = it.GetSurfElem();
	// Determine, which volume element is the right neighbour for the 
	// calculation;
	// our normal should point out of the correct neighbor volume element!
	if (   surfNeighborRegions_[vals] ==
              actSurfElem->ptVolElem1->regionId ) {
          ptVolElem = actSurfElem->ptVolElem1;
	  normSign = 1.0;
	} 
	else {
	  ptVolElem = actSurfElem->ptVolElem2;
	  normSign = -1.0;
	}
     
	normSign *= (Double) actSurfElem->normalSign;
     
	ptSurfElemFE = actSurfElem->ptElem; 
	ptVolElemFE = ptVolElem->ptElem;
     
	const StdVector<UInt> & surfConnect = actSurfElem->connect;
	const StdVector<UInt> & volConnect = ptVolElem->connect;
     
	// calculate volume integration coordinates from
	// surface integration coordinat for evalauting the 
	// gradient on the surface of the volume element
	ptSurfElemFE->GetCoordMidPoint(lCoordSurf);
	ptVolElemFE->GetLocalIntPoints4Surface(surfConnect, volConnect,
					       lCoordSurf, lCoordVol);
     
	BaseMaterial * myMat = materials_[ptVolElem->regionId];
	myMat->GetScalar(density,DENSITY,REAL);
     
	Matrix<Double> CornerCoords; 
	ptgrid_->GetElemNodesCoord( CornerCoords, surfConnect, false );
	Double area = ptSurfElemFE->CalcVolume( CornerCoords, isaxi_);
     
	// Calc gradient
        ElemList elList(ptgrid_);
        elList.SetElement( ptVolElem );
        EntityIterator it2 = elList.GetIterator();
	gradOp->CalcElemGradField(gradVal, it2, 
				  lCoordVol,1.0);
	// Calc global normal
	ptgrid_->CalcSurfNormal(normal, *actSurfElem);
     
	normal    *= normSign;
	gradNormal = normal * gradVal;

	//get average solution
	TYPE elemSol = 0;
	TYPE nodeSol;
	for ( UInt k=0; k<surfConnect.GetSize(); k++) {
	  solhelp->Get(solType, surfConnect[k]-1,0, nodeSol);
	  elemSol += nodeSol;
	}
	elemSol /= (Double)surfConnect.GetSize(); 
     
	// get the conjugate complex value
	elemSol = std::conj(elemSol);
     
	pdeElemNum = eqnMap_->Mesh2PdeElem(ptVolElem->elemNum);
	//pdeElemNum = eqnMap_->Mesh2PdeElem(actSurfElem->elemNum);
     
	if ( formulation_ == ACOU_PRESSURE ) {
	  elemPower = multVal * gradNormal * elemSol / density;
        }
	else {
	  elemPower = multVal * density * gradNormal * elemSol;
        }
     
	//take the real part and multiply it with surface
	elemPower *= area;

        actVal[regionIt.GetPos()] += elemPower;
        
      }
      
    }

    delete gradOp;


    } **/


  // ***********************************************************************
  //   Obtain information on desired output quantities from parameter file
  // ***********************************************************************
  void AcousticPDE::DefineAvailResults() {
    ENTER_FCN( "AcousticPDE::DefineAvailResults" );

    // === NODE POTENTIAL / PRESSURE (Primary Unknown) ===
    // check if problem is lagrange or legendre
    std::string approxType = myParam_->Get("type")->AsString();
    
    // Create new resultDof object
    shared_ptr<ResultInfo> res1(new ResultInfo);
    if ( formulation_ ==  ACOU_PRESSURE) {
      res1->resultType = ACOU_PRESSURE;
      res1->dofNames = "";
      res1->unit = "Pa";
    } else {
      res1->resultType = ACOU_POTENTIAL;
      res1->dofNames = "";
      res1->unit = "m^2/s";
    }
    
    if ( approxType == "lagrange" ) {
      shared_ptr<AnsatzFct> fct(new LagrangeFct);
      res1->definedOn = ResultInfo::NODE;
      res1->fctType = fct;
    } else {
      UInt order = myParam_->Get("order")->AsUInt();
      
      // Create new resultDof object
      shared_ptr<LegendreFct> fct(new LegendreFct);
      fct->SetIsoOrder( order );
      //fct->order_ = order;
      res1->definedOn = ResultInfo::PFEM;
      res1->fctType = fct;
    }
    res1->entryType = ResultInfo::SCALAR;
    availResults_.insert( res1 );
    results_.Push_back( res1 );      

    // === PRESSURE / POTENTIAL - 1.DERIVATIVE ===
    shared_ptr<ResultInfo> deriv1(new ResultInfo);
    if( formulation_ == ACOU_POTENTIAL ) {
      deriv1->resultType = ACOU_POTENTIAL_DERIV_1;
      deriv1->dofNames = "";
      deriv1->unit = "m^2/s^2";
    } else {
      deriv1->resultType = ACOU_PRESSURE_DERIV_1;
      deriv1->dofNames = "";
      deriv1->unit = "Pa/s";
    }
    deriv1->entryType = res1->entryType;
    deriv1->definedOn = res1->definedOn;
    deriv1->fctType = res1->fctType;
    availResults_.insert( deriv1 );

    // === PRESSURE / POTENTIAL - 2.DERIVATIVE ===
    shared_ptr<ResultInfo> deriv2(new ResultInfo);
    if( formulation_ == ACOU_POTENTIAL ) {
      deriv2->resultType = ACOU_POTENTIAL_DERIV_2;
      deriv2->dofNames = "";
      deriv2->unit = "m^2/s^3";
    } else {
      deriv2->resultType = ACOU_PRESSURE_DERIV_2;
      deriv2->dofNames = "";
      deriv2->unit = "Pa/s^2";
    }
    deriv2->entryType = res1->entryType;
    deriv2->definedOn = res1->definedOn;
    deriv2->fctType = res1->fctType;
    availResults_.insert( deriv2 );

    // === PRESSURE (element postprocessing results) ===
    if( formulation_ == ACOU_POTENTIAL ) {
      shared_ptr<ResultInfo> pres(new ResultInfo);
      pres->resultType = ACOU_PRESSURE;
      pres->dofNames = "";
      pres->unit = "Pa";
      pres->entryType = ResultInfo::SCALAR;
      pres->definedOn = ResultInfo::ELEMENT;
      pres->fctType = shared_ptr<ConstFct>(new ConstFct() );
      availResults_.insert( pres );
    }

    // === RHS VALUE ===
    shared_ptr<ResultInfo> rhs(new ResultInfo);
    rhs->resultType = ACOU_RHSVAL;
    rhs->dofNames = "";
    rhs->unit = "";
    rhs->entryType = res1->entryType;
    rhs->definedOn = res1->definedOn;
    rhs->fctType = res1->fctType;
    availResults_.insert( rhs );


    // === ACOU_INTENSITY ===
    shared_ptr<ResultInfo> intens(new ResultInfo);
    intens->resultType = ACOU_INTENSITY;
    if( dim_ == 3 ) { 
      intens->dofNames = "x", "y", "z";
    } else {
      if (isaxi_ ) {
        intens->dofNames = "r", "z";
      } else {
        intens->dofNames = "x", "y";
      }
    }
    intens->unit = "W/m^2"; 
    intens->entryType = ResultInfo::VECTOR;
    intens->definedOn = ResultInfo::SURF_ELEM;
    intens->fctType = shared_ptr<ConstFct>(new ConstFct() );
    availResults_.insert( intens );
    
    // === ACOU_POWER ===
    shared_ptr<ResultInfo> power(new ResultInfo);
    power->resultType = ACOU_POWER;
    power->dofNames = "";
    power->unit = "W";
    power->entryType = ResultInfo::SCALAR;
    power->definedOn = ResultInfo::SURF_REGION;
    power->fctType = shared_ptr<ConstFct>(new ConstFct() );
    availResults_.insert( power );

    // === ACOU_FORCE ===
    shared_ptr<ResultInfo> force(new ResultInfo);
    force->resultType = ACOU_FORCE;
    force->dofNames = "";
    force->unit = "N";
    force->entryType = ResultInfo::SCALAR;
    force->definedOn = ResultInfo::SURF_ELEM;
    force->fctType = shared_ptr<ConstFct>(new ConstFct() );
    availResults_.insert( force );

    // ===================================
    // Check for non-conforming interfaces
    // ===================================
    StdVector<std::string> ncIfaceNames, ncIfaceNamesForPDE;
    StdVector<RegionIdType> ncIfaceIds;
    
    LOG_DBG2(acoupde) << "NonMatching: Checking if nonconforming "
                      << "interfaces of PDE exist in domain.";

    ParamNode* domainNCIfaceListNode;
    domainNCIfaceListNode = param->Get("domain")->Get("ncInterfaceList", false);

    if(domainNCIfaceListNode)
    {

        StdVector<ParamNode*> pdeNCIfaceNodes = 
            param->Get("sequenceStep", "index", GenStr(sequenceStep_) )
            ->Get("pdeList")->Get("acoustic")->Get("ncInterfaceList")
            ->GetList("ncInterface");

        for (UInt i = 0; i < pdeNCIfaceNodes.GetSize(); i++) {
            std::string pdeIfaceName = pdeNCIfaceNodes[i]->Get("name")->AsString();
            std::string domainIfaceName;
        
            ParamNode* domainIfaceNode = domainNCIfaceListNode->Get("ncInterface",
                                                                    "name",
                                                                    pdeIfaceName,
                                                                    false);
            if(!domainIfaceNode)
            {
                LOG_DBG2(acoupde) << "NonMatching: Nonconforming "
                                  << "interface '" << ncIfaceNames[i]
                                  << "' does not exist in domain.";
            
                EXCEPTION( "ncInterface referenced from PDE not defined in domain!");
            }

            ncIfaceNamesForPDE.Push_back(pdeIfaceName);
        }
        ptgrid_->RegionNameToId( ncIfaceIds, ncIfaceNamesForPDE );
    
        for (UInt i = 0; i < ncIfaceIds.GetSize(); i++) {
            ncIFaces_.Push_back(ncIfaceIds[i]);
        }
    
        // In the case of the presence of non-conforming interfaces,
        // a second resultdof object has to be created, which describes the 
        // Lagrange multiplier
        if( ncIFaces_.GetSize() > 0 ) {
            LOG_DBG2(acoupde) << "NonMatching: Defining new ResultDof Lagrange.";
            shared_ptr<ResultInfo> lagr ( new ResultInfo );
            lagr->resultType = LAGRANGE_MULT;
            lagr->dofNames = "l";
            lagr->fctType = results_[0]->fctType;
            lagr->definedOn = results_[0]->definedOn;
            results_.Push_back( lagr );
        }
    }
    
  }
  
  
  void AcousticPDE::ReadFlowData() {
    ENTER_FCN( "AcousticPDE::ReadFlowDaa" );
    
    // check if bcsNode is present
    ParamNode * bcsNode = myParam_->Get("bcsAndLoads", false );
    if( !bcsNode) return;
    
    // get nodes with flowData
    StdVector<ParamNode*> flowNodes = bcsNode->GetList("flowData");
    for( UInt i = 0; i < flowNodes.GetSize(); i++ ) {
      std::string regionName = flowNodes[i]->Get("region")->AsString();
      RegionIdType regionId = ptgrid_->RegionNameToId( regionName );
      
      // store information about flow
      regionFlowNodes_[regionId] = flowNodes[i];
    }
    
  }




  //   Obtain information on desired output quantities from parameter file
  // ***********************************************************************
  void AcousticPDE::ReadDataPML(std::string& dampingTypePML, 
                                Matrix<Double>& inner, 
                                Double& dampPML, 
                                ParamNode * actNode ) {
  
    ENTER_FCN( "AcousticPDE::ReadDataPML" );

    // help variables for parameter checking
    StdVector<std::string> propGeo;
    StdVector<std::string> stringVal;
    StdVector<Double> val;

    // Check, if pml node has a child "prepRegion"
    ParamNode * propRegionNode = actNode->Get( "propRegion", false );

    // If no propagation region is defined explicitly, we 
    // let the method GetPMLLayerData() extract the geometric information
    // for the propagation region
    if( propRegionNode ) {
      
      //resize data for ptopagation region
      inner.Resize(2,dim_);
      inner.Init();
      
      //xMin
      propRegionNode->Get( "xMin", inner[0][0] );
      
      //xMax
      propRegionNode->Get( "xMax", inner[1][0] );
      
      //yMin
      propRegionNode->Get( "yMin", inner[0][1] );
      
      //yMax
      propRegionNode->Get( "yMax", inner[1][1] );
      
      if ( dim_ == 3 ) {
        //zMin
        propRegionNode->Get( "zMin", inner[0][2] );
        
        //zMax     
        propRegionNode->Get( "zMax", inner[1][2] );
      }
    }
    
    //get type of damping function
    actNode ->Get( "type", dampingTypePML );

    //get factor for damping function
    actNode->Get( "dampFactor", dampPML );

  }


  // ***********************************************************************
  //   Obtain information on desired output quantities from parameter file
  // ***********************************************************************
  void AcousticPDE::GetPMLLayerData(Matrix<Double>& inner, 
                                    Matrix<Double>& outer,
                                    RegionIdType actRegion )  {  

    ENTER_FCN( "AcousticPDE::GetPMLLayerData" );

    // inner/outer:   xmin  ymin  zmin
    //                xmax  ymax  zmax

    if ( inner.GetSizeCol() != dim_ ) {
      //we have to compute it, since the user has not specified it

      inner.Resize(2,dim_);
      inner.Init();
      
      for (UInt isd = 0; isd < subdoms_.GetSize(); isd++) {
        if ( subdoms_[isd] != actRegion ) {
          StdVector<Elem*> elemssd;
          ptgrid_->GetElems(elemssd, subdoms_[isd] );
          
          for (UInt actEl=0; actEl< elemssd.GetSize(); actEl++) {
            StdVector<UInt> & connecth = elemssd[actEl]->connect;
            
            Matrix<Double> ptCoord;
            ptgrid_->GetElemNodesCoord(ptCoord, connecth,  false );
            for (UInt i=0; i< ptCoord.GetSizeCol(); i++) {
              //minInnerX
              if ( ptCoord[0][i] < inner[0][0] )
                inner[0][0] = ptCoord[0][i];
              
              //minInnerY
              if ( ptCoord[1][i] < inner[0][1] )
                inner[0][1] = ptCoord[1][i];
              
              if ( dim_ > 2 ) {
                //minInnerZ
                if ( ptCoord[2][i] < inner[0][2] )
                  inner[0][2] = ptCoord[2][i];
              }
              
              //maxInnerX
              if ( ptCoord[0][i] > inner[1][0] )
                inner[1][0] = ptCoord[0][i];
              
              //maxInnerY
              if ( ptCoord[1][i] > inner[1][1] )
                inner[1][1] = ptCoord[1][i];
              
              if ( dim_ > 2 ) {
                //maxInnerZ
                if ( ptCoord[2][i] > inner[1][2] )
                  inner[1][2] = ptCoord[2][i];
              }
            }
          }
        }
      }
      std::ostringstream out;
      out.clear();
      out << "Acoustic propagation region:\n" 
          << "   xmin = " << inner[0][0] << std::endl
          << "   xmax = " << inner[1][0] << std::endl
          << "   ymin = " << inner[0][1] << std::endl
          << "   ymax = " << inner[1][1] << std::endl;
      if ( dim_ == 3) {
        out << "   zmin = " << inner[0][2] << std::endl
            << "   zmax = " << inner[1][2] << std::endl;
      }
      out << std::endl;
      Info->PrintF( pdename_, out.str().c_str() );
    }
    
    outer.Resize(inner.GetSizeRow(),inner.GetSizeCol());
    
    outer[0][0] = outer[1][0] = inner[1][0];
    outer[0][1] = outer[1][1] = inner[1][1];
    if (inner.GetSizeCol() > 2 ) {
      outer[0][2] = outer[1][2] = inner[1][2];
    }
    
    StdVector<Elem*> elemssd;
    ptgrid_->GetElems(elemssd, actRegion );
    
    for (UInt actEl=0; actEl< elemssd.GetSize(); actEl++) {
      StdVector<UInt> & connecth = elemssd[actEl]->connect;
      
      Matrix<Double> ptCoord;
      ptgrid_->GetElemNodesCoord(ptCoord, connecth,  false );
      for (UInt i=0; i< ptCoord.GetSizeCol(); i++) {
        //minXPML
        if ( ptCoord[0][i] < outer[0][0] )
          outer[0][0] = ptCoord[0][i];

        //minYPML
        if ( ptCoord[1][i] < outer[0][1] )
          outer[0][1] = ptCoord[1][i];
        
        if (inner.GetSizeCol() > 2 ) {
          //minZPML
          if ( ptCoord[2][i] < outer[0][2] )
            outer[0][2] = ptCoord[2][i];
        }
        
        //maxXPML
        if ( ptCoord[0][i] > outer[1][0] )
          outer[1][0] = ptCoord[0][i];
        
        //maxYPML
        if ( ptCoord[1][i] > outer[1][1] )
          outer[1][1] = ptCoord[1][i];

        if (inner.GetSizeCol() > 2 ) {
          //maxZPML
          if ( ptCoord[2][i] > outer[1][2] )
            outer[1][2] = ptCoord[2][i];
        }
      }
      
    }
    
  }


  void AcousticPDE::ReadDataDampLayer( std::string& dampingTypePML, 
                                       Vector<Double>& mPoint, 
                                       Double& dampFactor, 
                                       Double& dampFactorMax, 
                                       Double& startRadius, 
                                       Double& endRadius, 
                                       ParamNode * actNode ) {
    
    ENTER_FCN( "AcousticPDE::ReadDataDampLayer" );
    
    
    
    StdVector<std::string> stringVal;
    
    // Construct vectors for restricted parameter search
    actNode->Get( "type", dampingTypePML );
    
    //get the center point
    mPoint.Resize(dim_);

    //xM, yM,
    actNode->Get( "xM", mPoint[0] );
    actNode->Get( "yM", mPoint[1] );

    if ( dim_ == 3 ) {
      //zM
      actNode->Get( "zM", mPoint[2] );
    }

    //start radius, end radius, dampFactor
    actNode->Get( "radiusStart", startRadius );
    actNode->Get( "radiusEnd", endRadius );
    actNode->Get( "dampFactor", dampFactor );

    //get maximal damping factor (saturation)
    actNode->Get( "dampFactorMax", dampFactorMax );
  }

} // end of namespace
