#include <fstream>
#include <iostream>
#include <string>
#include <math.h>

#include "acousticPDE.hh" 
#include "DataInOut/Unverg/outUnverg.hh"
#include "DataInOut/GMV/outGMV.hh"
#include "Forms/forms_header.hh"
#include "Forms/abcInt.hh"
#include "Forms/linNeumannInt.hh"
#include "Estimator/spaceerror.hh"
#include "newmark.hh"
#include "newmarkFracDamp.hh"
#include "DataInOut/WriteInfo.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"
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



  // =========================================================================
  // set solution information
  // =========================================================================
  AcousticPDE::AcousticPDE(Grid * aptgrid, TimeFunc *aptTimeFunc, 
                           WriteResults *aptOut)
    :SinglePDE(aptgrid,aptOut,aptTimeFunc) {

    ENTER_FCN( "AcousticPDE::AcousticPDE" );

    pdename_          = "acoustic";
    pdematerialclass_ = FLUID;
    fracDamping_ = false;

    isMechCoupled_ = false;
    isNrbcCoupled_ = false;
    saveRHSval_ = false;
    saveRHSvalHist_ = false;
    isBubbleCoupled_ = false;


    // PDE formulation either in acoustic potential or pressure
    std::string str;
    params->Get( "formulation", str, "pdeList", pdename_ );
    String2Enum( str, formulation_ );
    str = "Using * " + str + " as state variable in formulation of PDE\n";
    Info->PrintF( pdename_, str.c_str() );


    // timestepping formulation
    params->Get( "timeSteppingFormulation", str, "pdeList", pdename_ );
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

    // axisymmetric setup    
    isaxi_ = params->HasValue( "type", "axi", "geometry" );

    //To check if acoustic is coupled with nrbc and set isNrbcCoupled
    if( params->IsSet( "isCoupledNrbc","pdeList" ,"acoustic" ) ) {
      isNrbcCoupled_=true;
    }
    // Create new resultDof object
    shared_ptr<ResultDof> res1(new ResultDof);
    shared_ptr<AnsatzFct> fct(new LagrangeFct);
    if ( formulation_ ==  ACOU_PRESSURE) {
      res1->resultType = ACOU_PRESSURE;
      res1->dofNames = "p";
    } else {
      res1->resultType = ACOU_POTENTIAL;
      res1->dofNames = "vp";
    }
    res1->definedOn = ResultDof::NODE;
    res1->fctType = fct;
    results_.Push_back( res1 );
  }


  // *********************************************
  //   Check what type of damping should be used
  // *********************************************
  void AcousticPDE::ReadDampingInformation( ) {

    ENTER_FCN( "AcousticPDE::ReadDampingInformation" );

    fracMemory_ = 0;
    Integer firstFrac=-1;

    // Construct vectors for restricted search parameter
    StdVector<std::string> keyVec;
    StdVector<std::string> attrVec;
    StdVector<std::string> valVec;

    for (UInt k = 0; k < subdoms_.GetSize(); k++) {

      // Be aware, that the RegionId entries in subdoms_ need not be in order!
      // We follow the order in subdoms_, which has to be in acordance with
      //  AcousticPDE::DefineIntegrators()

      RegionIdType actRegion = subdoms_[k];
      std::string actRegionName;
      actRegionName = ptgrid_->RegionIdToName( actRegion );
      keyVec = "acoustic" , "region" , "damping" , "type";
      attrVec= ""         , "name"   , "";
      valVec = ""         , actRegionName, "";
      StdVector<std::string> dampInfo;
      params->GetList( keyVec, attrVec, valVec, dampInfo);

      if ( dampInfo.IsEmpty() ) {
        dampingList_[actRegion] = NONE;
        Info->PrintF( pdename_, 
                      "No information specifying damping detected!\n" );
      }
      else if (dampInfo[0] == "no") {
        dampingList_[actRegion] = NONE;
        Info->PrintF( pdename_, 
                      "      * NO damping at all for region: %s\n", 
                      actRegionName.c_str() );;
      }
      else if (dampInfo[0] == "rayleigh") {
        dampingList_[actRegion] = RAYLEIGH;
        Info->PrintF( pdename_, 
                      "      * RAYLEIGH damping for region: %s\n",
                      actRegionName.c_str() );
      }
      else if (dampInfo[0] == "thermoViscous") {
        dampingList_[actRegion] = THERMOVISCOUS;
        Info->PrintF( pdename_, 
                      "      * THERMOVISCOUS damping for region: %s\n",
                      actRegionName.c_str() );
      }
      else if (dampInfo[0] == "fractional") {
        fracDamping_ = true;
        Info->PrintF( pdename_, 
                      "      * FRACTIONAL damping for region: %s\n",
                      actRegionName.c_str() );

        // Find first region containing fractional damping
        if ( firstFrac < 0 )
          firstFrac = k;

        // Gather additional information for fractional damping model
        keyVec = "acoustic" , "region" , "damping" , "fracAlg";
        StdVector<std::string> fracAlg;
        params->GetList( keyVec, attrVec, valVec, fracAlg );

        keyVec = "acoustic" , "region" , "damping" , "fracMemory";
        StdVector<UInt> fracMem;
        params->GetList( keyVec, attrVec, valVec, fracMem );

        keyVec = "acoustic" , "region" , "damping" , "interpolation";
        StdVector<std::string> interpol;
        params->GetList(  keyVec, attrVec, valVec, interpol );


        // Include fracAlg and interpolation info in dampingList
        if( fracAlg.IsEmpty() || fracMem.IsEmpty() || interpol.IsEmpty() ) {
          (*error) << "Specify attributes fracAlg, fracMemory " 
                   << "and interpolation for fractional damping model!";
          Error( __FILE__, __LINE__ ); 
        }
        else if  ( fracAlg[0] == "gl" ) {
          Info->PrintF( "", "\t\t\t using Gruenwald-Letnikov algorithm,\n");
          if (interpol[0] == "no" )
            dampingList_[actRegion] = FRACTIONAL_GL;
          else {
            dampingList_[actRegion] = FRACTIONAL_GL_INT;
            Info->PrintF("", 
                         "\t\t\t linear interpol. of single past values\n\n");
          }
        }
        else if ( fracAlg[0] == "blank" ) {
          Info->PrintF( "", "\t\t\t using Blanks algorithm,\n");
          if (interpol[0] == "no" )
            dampingList_[actRegion] = FRACTIONAL_BLANK;
          else {
            dampingList_[actRegion] = FRACTIONAL_BLANK_INT;
            Info->PrintF("", 
                         "\t\t\t linear interpol. of single past values\n\n");
          }
        }

        // up to now take maximum of specified fracMemory values
        if ( fracMem[0] > fracMemory_ )
          fracMemory_ = fracMem[0];
      }
    }

    if ( fracDamping_== true ) {
      Info->PrintF(pdename_, "Memory size for fractional damping  is: %d\n",
                   fracMemory_ );
    }
  }

  void AcousticPDE::ReadSpecialBCs() {
    ENTER_FCN( "AcousticPDE::ReadSpecialBCs" );

    // ***************************************************************
    //   If no other damping type is specified and we have absorbing
    //   boundary conditions, then use ABCDAMP
    // ***************************************************************
    StdVector<std::string> auxVec;
    absorbingBCs_ = false;
    params->GetList( "name", auxVec, pdename_, "absorbingBCs" );
    ptgrid_->RegionNameToId( absBCs_, auxVec );

    if ( absBCs_.GetSize() ) {
      absorbingBCs_ = true;
      Info->PrintF( pdename_, "Apply Absorbing Boundary Conditions\n" );
    }
  }


  // *************************************************************
  //   Check what type of nonlinear PDE formulation should be used
  // *************************************************************
  void  AcousticPDE::InitNonLin() {
    ENTER_FCN( "AcousticPDE::InitNonLin" );

    nonLin_ = false; //declaration in StdPDE.hh
    StdVector<std::string> strVec;
    params->GetList( "nonLinear", strVec, pdename_, "region" );
    nonLinPDEName_.Resize(strVec.GetSize());
    nonLinPDEName_.Init();
    
    for ( UInt k = 0; k < strVec.GetSize(); k++ ) {
      
      if ( strVec[k] == "no" )
        ;
      else if ( strVec[k] == "westervelt" ) {
        nonLin_ = true;
        nonLinPDEName_[k] = WESTERVELT;
        Info->PrintF(pdename_, "      * Westervelt equation for region: %d\n",
                     k );
        if ( formulation_ == ACOU_POTENTIAL )
          (*error) << "Acoustic potential formulation not supported for "
                   << "Westervelt equation";
        Error( __FILE__, __LINE__ );
      }
      else if ( strVec[k] == "kuznetsov" ) {
        nonLin_ = true;
        nonLinPDEName_[k] = KUZNETSOV;
        Info->PrintF(pdename_, "      * Kuznetsov equation for region: %d\n",
                     k );
        if ( formulation_ == ACOU_PRESSURE )
          (*error) << "Acoustic pressure formulation not supported for"
                   << "Kuznetsov equation!";
        Error( __FILE__, __LINE__ );
      }
    }
    
    if( nonLin_ ) {
      // solution method
      params->Get("method", nonLinMethod_, pdename_, "nonLinear" );
      // perform logging?
      nonLinLogging_ = params->IsSet( "logging", pdename_, "nonLinear" );
      // type of line search
      params->Get("type", lineSearch_, pdename_, "lineSearch" );
      // incremental stopping criterion
      params->Get("incStopCrit", incStopCrit_, pdename_, "nonLinear" );
      // residual stopping criterion
      params->Get("resStopCrit", residualStopCrit_, pdename_, "nonLinear" );
      // maximal number of NL-iterations
      params->Get("maxNumIters", nonLinMaxIter_, pdename_, "nonLinear");
    }
  }


  void AcousticPDE::DefineIntegrators() {

    ENTER_FCN( "AcousticPDE::DefineIntegerators" );

    // help variables for parameter checking
    StdVector<std::string> keyVec;
    StdVector<std::string> attrVec;
    StdVector<std::string> valVec;

    Double density, compressibility, c0;
    Double coeffmass, coeffdamp;

    for (UInt actSD = 0; actSD < subdoms_.GetSize(); actSD++) {

      // create new entity list
      shared_ptr<ElemList> actSDList( new ElemList(ptgrid_ ) );
      actSDList->SetRegion( subdoms_[actSD] );

      // ********************************************************************
      //   Attention:
      //   In AcousticPDE ALL integrators are multiplied with density!
      // ********************************************************************

      materials_[subdoms_[actSD]]->GetScalar(density,DENSITY,REAL);
      materials_[subdoms_[actSD]]->GetScalar(compressibility,ACOU_BULK_MODULUS,REAL);
      c0 = sqrt(compressibility/density);

      // if pde couples with mechanic, we have to multiply the density by -1
      if ( isMechCoupled_ == true && formulation_ !=  ACOU_PRESSURE) {
        density *= -1.0;
      }

      // ********************************************************************
      //   Linear wave equation
      // ********************************************************************


      //check  for PML!
      RegionIdType actRegion = subdoms_[actSD];
      std::string actRegionName;
      actRegionName = ptgrid_->RegionIdToName( actRegion );
      keyVec = "acoustic" , "region" , "pml" , "type";
      attrVec= ""         , "name"   , "";
      valVec = ""         , actRegionName, "";
      StdVector<std::string> pmlInfo;
      params->GetList( keyVec, attrVec, valVec, pmlInfo);

      if (pmlInfo.GetSize() > 0) {
        if ( analysistype_ != HARMONIC ) {
          Error("PML just supported for Harmonic-Analysis",__FILE__,__LINE__);
        }

        //read data for PML layer

        //type of PML damping
        std::string dampingTypePML;

        // inner / outer region
        Matrix<Double> inner;
        Matrix<Double> outer;

        //damping factor
        Double dampPML;

        ReadDataPML(dampingTypePML, inner, dampPML, actRegion);
        dampPML *= c0;

        GetPMLLayerData(inner, outer, actSD);

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
        stiffContextReal->SetMatDataType(matType);
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
        stiffContextImag->SetMatDataType(matType);
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
        massContextReal->SetMatDataType( matType );
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
        massContextImag->SetMatDataType( matType );
        assemble_->AddBiLinearForm( massContextImag );

        // Finally add the stiffness/mass integrators
        assemble_->AddBiLinearForm( stiffContextReal );
        assemble_->AddBiLinearForm( stiffContextImag );
        assemble_->AddBiLinearForm( massContextReal );
        assemble_->AddBiLinearForm( massContextImag );
      }

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

	RegionIdType regionFlow = subdoms_[actSD];
	std::string regionName;
	regionName = ptgrid_->RegionIdToName( regionFlow );
	keyVec = "acoustic" , "region" , "flowData" , "flowDir";
	attrVec= ""         , "name"   , "";
	valVec = ""         , regionName, "";
	StdVector<std::string> flowInfo;
	params->GetList( keyVec, attrVec, valVec, flowInfo);

	if (flowInfo.GetSize() > 0) {

	  if ( formulation_ == ACOU_PRESSURE )
	    Error("Pierce-Equation just possible in velocity potential formulation",
		  __FILE__, __LINE__);

	  //read flow data 
	  SimpleFlow* flowData = new SimpleFlow();
	  flowData->ReadFlowData(pdename_, regionName, dim_);

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
        keyVec = "acoustic" , "region" , "dampLayer" , "type";
        attrVec= ""         , "name"   , "";
        valVec = ""         , actRegionName, "";
        StdVector<std::string> dampInfo;
        params->GetList( keyVec, attrVec, valVec, dampInfo);

        if (dampInfo.GetSize() > 0) {
	  
          //type of damping fnc
          std::string dampFnc;

          //damping data
          Double dampFactor, dampFactorMax, startRadius, stopRadius;
          Vector<Double> mPoint;
          ReadDataDampLayer(dampFnc, mPoint, dampFactor, dampFactorMax, 
                            startRadius, stopRadius, actRegion);

          //get the Rayleigh material parameters
          Double alpha, beta;
          materials_[subdoms_[actSD]]->GetScalar(alpha,RAYLEIGH_ALPHA,REAL);
          materials_[subdoms_[actSD]]->GetScalar(beta,RAYLEIGH_BETA,REAL);
    
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
	  
          if (dampingList_[subdoms_[actSD]] == RAYLEIGH) {
            // This works even after assemble_->AddIntegrator() is executed
            //   because of the pointers...
            Double alpha, beta;
            materials_[subdoms_[actSD]]->GetScalar(alpha,RAYLEIGH_ALPHA,REAL);
            materials_[subdoms_[actSD]]->GetScalar(beta,RAYLEIGH_BETA,REAL);
    
            // stiffness part
            stiffContext->SetSecDestMat( DAMPING, beta );
	    
            // mass part
            massContext->SetSecDestMat( DAMPING, alpha );
          }
	  
	  
          else if ( dampingList_[subdoms_[actSD]] == THERMOVISCOUS ) {
            Double viscousAlpha;
            materials_[subdoms_[actSD]]->GetScalar(viscousAlpha, ACOU_ALPHA,REAL);

            coeffdamp  =  density * 2.0 * viscousAlpha * c0;
            BaseForm * bilinearStiff  = new LaplaceInt(coeffdamp, isaxi_);  
            BiLinFormContext * dampContext = 
              new BiLinFormContext(bilinearStiff, DAMPING );
            dampContext->SetResults( results_[0], results_[0],
                                     actSDList, actSDList );
            dampContext->SetPtPdes(this, this);   
            assemble_->AddBiLinearForm( dampContext );
          }
	  
          else if ( dampingList_[subdoms_[actSD]] == FRACTIONAL_GL ||
                    dampingList_[subdoms_[actSD]] == FRACTIONAL_GL_INT ) {

            Double fracAlpha, fracExp;
            materials_[subdoms_[actSD]]->GetScalar(fracAlpha,ACOU_ALPHA,REAL);
            materials_[subdoms_[actSD]]->GetScalar(fracExp,FRACTIONAL_EXPONENT,REAL);

            coeffdamp = - density * 2.0 * fracAlpha / c0 / sin((fracExp-1.0)*PI/2.0);

            BaseForm * bilinearDamp  = 
              new MassInt(coeffdamp, 1, isaxi_);
            bilinearDamp->SetFracDamping();
            Double fracDampCoeff = GetFracDampMatrixCoeff( subdoms_[actSD] );
            bilinearDamp->SetFactor( fracDampCoeff );

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
	  
          else if  ( dampingList_[subdoms_[actSD]] == FRACTIONAL_BLANK ||
                     dampingList_[subdoms_[actSD]] == FRACTIONAL_BLANK_INT ) {

            Double fracAlpha, fracExp;
            materials_[subdoms_[actSD]]->GetScalar(fracAlpha,ACOU_ALPHA,REAL);
            materials_[subdoms_[actSD]]->GetScalar(fracExp,FRACTIONAL_EXPONENT,REAL);

            coeffdamp =  - density * 2.0 * fracAlpha / c0 / sin((fracExp-1.0)*PI/2.0);
            // prefactor of blank alg
            coeffdamp *= exp(-gammaln(1.0- (fracExp- 1.0)) ); 
            // weight factor of index 0
            coeffdamp *= 1.0/(1.0- (fracExp- 1.0));           
            BaseForm * bilinearDamp  = 
              new MassInt(coeffdamp, 1, isaxi_);
            bilinearDamp->SetFracDamping();
            Double fracDampCoeff = GetFracDampMatrixCoeff( subdoms_[actSD] );
            bilinearDamp->SetFactor( fracDampCoeff );

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
    keyVec = "multiSequence", "step", "pde","type";
    attrVec = "", "index", "";
    valVec = "", GenStr(bcSequenceIndex_), "";
    StdVector<std::string> pdeSearchBubble;
    params->GetList( keyVec, attrVec, valVec,pdeSearchBubble );    
    for (UInt i=0;i< pdeSearchBubble.GetSize();i++) {
      if (pdeSearchBubble[i] == "bubble"){
	isBubbleCoupled_ = true;
	Info->PrintF(  pdename_, " will use bubble-bilinearform\n");
      }
    }

    if( isBubbleCoupled_){

	  //get paramters and starting values for the computation of coefficients
	  Double sonicVel, densityforbubble, viscosity, bubbleDensity, frequency;
	  params->Get( "density", densityforbubble, "bubble", "bubbleDynamic" );
	  params->Get( "sonicVel", sonicVel, "bubble", "bubbleDynamic" );
	  params->Get( "viscosity", viscosity, "bubble", "bubbleDynamic" );
	  params->Get( "bubbleNumDensity", bubbleDensity, "bubble", "bubbles" );
	  params->Get( "frequency", frequency, "bubble", "excitation" );

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



	  std::cerr<<"test in acou define integrators"<< std::endl;
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
        new LinearFormContext( neumannBC, 
                               actBc.dynamics );
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
        
        SurfForm * bilinear_damp = new AbsorbingBCsInt(isaxi_);
        bilinear_damp->SetFirstVoluInfo(pdename_, materials_ );

        // In the case of acou-mech coupling we have to multiply the 
        // abc-Integrator matrix with -1
        if ( isMechCoupled_ == true ) {
          bilinear_damp->SetFactor(-1.0);
        }
        // In the case of acou-nrbc coupling we have to multiply the 
        // abc-Integrator matrix with C0 and multiply it by nrbcBeta0
        if ( isNrbcCoupled_ == true ) {
          bilinear_damp->SetFactor(sqrt( compressibility / density ));
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
  }

  void AcousticPDE::DefineSolveStep() {
    ENTER_FCN( "AcousticPDE::DefineSolveStep" );

    solveStep_ = new SolveStepAcoustic(*this);
  }



//   // ======================================================
//   // ALGSYS SECTION (SOLVER, ...) need for acoububble coupling
//   // ======================================================
//   void AcousticPDE::DefineAlgSys() {
    
//     ENTER_FCN( "Acoustic::DefineAlgSys" );

//     //=====================================================
//     // Only if acousticPDE is iteratively coupled with
//     // the bubblePDE
//     //====================================================

//     if( isBubbleCoupled_){
//       // Process input couplings of the bubblePDE to the
//       // additional baseforms of the acousticPDE 
//       UInt numInCouplings = ptCoupling_->GetNumInputCouplings();
      
//       // bubble coupling data
//       StdVector<UInt> elemNumbers;
//       Vector<Double> * radius = NULL;
//       Vector<Double> * radiusDeriv  = NULL;
//       std::cerr<<"definealssys numOutCouplings "<<numInCouplings<< std::endl;
//       for (UInt i = 0; i < numInCouplings; i++) {
      
// 	std::cerr<<"definealgsys ptCoupling_->GetInputQuantity(i)"<<ptCoupling_->GetInputQuantity(i)<< std::endl;
// 	if (ptCoupling_->GetInputQuantity(i) == BUBBLE_RADIUS) {
	
// 	  // get coupling elements and create vector with element numbers
// 	  StdVector<Elem*> * elems = NULL;
// 	  ptCoupling_->GetInputElements(i, elems );
 	
// 	  elemNumbers.Resize( elems->GetSize() );
// 	  for (UInt iElem = 0; iElem < elems->GetSize(); iElem++ ) {
// 	    elemNumbers[iElem] = (*elems)[iElem]->elemNum;
// 	  }
	
// 	  // get buffer with coupling values 
// 	  CFSVector * temp;
// 	  ptCoupling_->GetInputValues( i, temp );    
// 	  std::cerr<<"acou definealgsys size temp"<< (* temp).GetSize()<<std::endl;

// 	  radius = dynamic_cast<Vector<Double>*>(temp);

// 	} else if (ptCoupling_->GetInputQuantity(i) == BUBBLE_RADIUS_DERIV_1 ) {
	
// 	  // get buffer with coupling values 
// 	  CFSVector * temp;
// 	  ptCoupling_->GetInputValues( i, temp );
// 	  radiusDeriv = dynamic_cast<Vector<Double>*>(temp);
	
// 	}
// 	std::cerr<<"test in acou definealgsys nach"<< std::endl;
//       }

//       // Iterate over all BubbleMass/Stiff-Integrators and pass the
//       // radius data to them

//       std::map<RegionIdType, BubbleDampInt*>::iterator dampIt;
//       std::map<RegionIdType, BubbleStiffInt*>::iterator stiffIt;

//       for( dampIt = bubbleDampIntMap_.begin();
// 	   dampIt != bubbleDampIntMap_.end(); dampIt++ ) {
// 	dampIt->second->SetValues( elemNumbers, radius, radiusDeriv);
//       }

//       for( stiffIt = bubbleStiffIntMap_.begin();
// 	   stiffIt != bubbleStiffIntMap_.end();
// 	   stiffIt++ ) {
// 	stiffIt->second->SetValues( elemNumbers, radius, radiusDeriv);
//       }
//     }
    
//     SinglePDE::DefineAlgSys();

//   }




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
                                       subdoms_, dampingList_ );
      else
        Error("This needs to be implemented!",__FILE__,__LINE__);
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
            (*error) << "AcousticPDE: Coupling region is not within the "
                     << "subdomains of the PDE!";
            Error( __FILE__, __LINE__ );
          }

          // get elements belonging to subdomain
          StdVector<Elem*> elemssd;
          ptgrid_->GetVolElems(elemssd, subdoms_[SDidx]);

          StdVector<UInt> * couplingnodes = NULL;
          ptCoupling_->GetOutputNodes(i, couplingnodes);
          if ( couplingnodes == NULL ) {
            (*error) << "The pointer 'couplingnodes' is NULL!";
            Error( __FILE__, __LINE__ );
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
    
    Error( "Not working at the moment", __FILE__, __LINE__ );
    
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
    
    Error( "Not working at the moment", __FILE__, __LINE__ );

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
    //                  << actSurfCoupleElem->elemNum << " do not belong to my regions!";
    //         Error( __FILE__, __LINE__ );
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
      StdVector<Elem*> elemssd;
      ptgrid_->GetVolElems(elemssd, subdoms_[SDidx]);
      
      for (UInt actEl=0; actEl< elemssd.GetSize(); actEl++) {

        SourceOp->CalcElemPD(elemPowerDensity, elemssd[actEl], density);

        // Add the element energy to the according coupling node
        StdVector<UInt> connecth = elemssd[actEl]->connect;
        for (UInt elnode=0; elnode<connecth.GetSize(); elnode++) {

          sourceValue[elemNodeToCouplingNode[actEl+offset][elnode]]
            += elemPowerDensity[elnode];
        }
      }
      
      //in the case, that we have more than one coupling region!
      offset = elemssd.GetSize();
    }
  }
  void AcousticPDE::CalcBubblePressure( StdVector<Elem*> & elems,
                                        Vector<Double>& couplVals,
                                        SolutionType solType ) {

    ENTER_FCN( "AcousticPDE::CalcBubblePressure" );


    // Initialize coupling values
    couplVals.Init();

    if ( solType == ACOU_PRESSURE ) {

      Vector<Double> elPressure;
      for (UInt iElem = 0; iElem < elems.GetSize(); iElem++ ) {

        StdVector<UInt> & myConnect = elems[iElem]->connect;
        sol_->GetElemSolution( elPressure, myConnect );
        Double pressure = 0.0;

        for( UInt iNode = 0; iNode < elPressure.GetSize(); iNode++ ) {
          pressure += elPressure[iNode];
        }
        pressure /= elPressure.GetSize();

        // store back to vector
        couplVals[iElem] = pressure;
      }

    } else {
      Vector<Double> elPressureDeriv;
      for (UInt iElem = 0; iElem < elems.GetSize(); iElem++ ) {

        StdVector<UInt> & myConnect = elems[iElem]->connect;
        GetDerivSolVecOfElement( elPressureDeriv, myConnect );
        Double pressureDeriv = 0.0;

        for( UInt iNode = 0; iNode < elPressureDeriv.GetSize(); iNode++ ) {
          pressureDeriv += elPressureDeriv[iNode];
        }
        pressureDeriv /= elPressureDeriv.GetSize();

        // store back to vector
        couplVals[iElem] = pressureDeriv;
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

  void AcousticPDE::PostProcess() {
    ENTER_FCN( "AcousticPDE::PostProcess" );

    StdVector<Elem*> saveElems, temp;
    for ( UInt i = 0; i < saveElemForceHist_.GetSize(); i++ ) {

      ptgrid_->GetElemsByName(temp, saveElemForceHist_[i]);
      for( UInt j = 0; j < temp.GetSize(); j++ ) {
        saveElems.Push_back(temp[j]);
      }
    }

    // compute force from pressure/acoustic potential at surface element
    if (saveElemForceHist_.GetSize() != 0 && isComplex_ == false ) {
      CalcForce(saveElems);
    }

    // compute pressure from acoustic potential
    if (calcElemPressure_.GetSize() != 0) {
      if ( isComplex_ == false) {
        CalcElemPressure<Double>();
      } 
      else {
        CalcElemPressure<Complex>();
      }
    }
     
    // *** Acoustic Power ***
    if ( calcAcouPower_.GetSize() !=0 ) {
      if (analysistype_ == HARMONIC ) {
        CalcAcouPower<Complex>();
      }
    }

    // Last but no least trigger postprocessing fromt within script-file
#ifdef USE_SCRIPTING
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

  void AcousticPDE::CalcForce( StdVector<Elem*> & saveElems ) {  
    ENTER_FCN( "AcousticPDE::CalcForce" );

    Matrix<Double> ptCoord;
    Vector<Double> forceVec(saveElems.GetSize());
	
    for (UInt actEl=0; actEl < saveElems.GetSize(); actEl++) {

      // Perform cast from volume element to surface element,
      //  since calculation of force makes only sense on a surface
      SurfElem * actSaveElem = 
        dynamic_cast<SurfElem*> ((saveElems)[actEl]);

      if (actSaveElem == NULL) {
        Error( "No elements specified in storeResults section!"
               , __FILE__, __LINE__ );
      }

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
          (*error) << "AcousticPDE::CalcForce: The two volume element"
                   << "neighbours of surface element no. "
                   << actSaveElem->elemNum << " don't belong to my region!";
          Error( __FILE__, __LINE__ );
        }
      
        // Assign correct density
        Double density;
        materials_[subdoms_[matIndex]]->GetScalar(density,DENSITY,REAL);
		
        // retrieve 1st derivative, since F = rho * dpsi/dt * A
        GetDerivSolVecOfElement(valueElem, connect);

        valueElem *= density;
      }
      else if ( formulation_ == ACOU_PRESSURE ) {

        // retrieve solution, since F = p * A
        GetSolVecOfElement(valueElem,connect);
      }

      const UInt nrIntPts= ptElem->GetNumIntPoints();
      const Vector<Double> & intWeights = ptElem->GetIntWeights();

      Vector<Double> forceElem(1); // is defined as vector, because 
      // the method SetElemResult can only handle Vector<Double>
      forceElem[0] = 0;
      Double jacDet;
      for (UInt actIntPt=1; actIntPt<=nrIntPts;  actIntPt++) {

        jacDet = ptElem->CalcJacobianDetAtIp(actIntPt, ptCoord);
        Vector<Double> shapeFnc;  
        ptElem -> GetShFncAtIp(shapeFnc, actIntPt);

        if (isaxi_) {
          Vector<Double> coordAtIP;
          coordAtIP = ptCoord * shapeFnc;
          forceElem[0] +=  (intWeights[actIntPt-1] * jacDet
                            * 2 * PI) * coordAtIP[0] * (valueElem * shapeFnc);
        }
        else {
          forceElem[0] +=  intWeights[actIntPt-1] * jacDet 
            * (shapeFnc * valueElem);
        }
      }

      // map element result back in global set of results
      UInt pdeElem;
      pdeElem = eqnMap_->Mesh2PdeElem(actSaveElem->elemNum);
      acouForce_.SetElemResult(pdeElem-1,forceElem);

      forceVec[actEl] = forceElem[0];
    }

    sumAcouForce_ = 0;
    for(UInt k=0; k<saveElems.GetSize(); k++) {
      sumAcouForce_ += forceVec[k];
    }
  }

  template <class TYPE>
  void AcousticPDE::CalcElemPressure( ) {  
    ENTER_FCN( "AcousticPDE::CalcElemPressure" );

    // loop over all subdomains
    for (UInt actSD=0; actSD<calcElemPressure_.GetSize(); actSD++) {

      // get all elements belonging to subdomain
      StdVector<Elem*> elemsSD;
      ptgrid_->GetVolElems( elemsSD, calcElemPressure_[actSD] );

      // density of elements in subdomain
      Double density;
      materials_[subdoms_[actSD]]->GetScalar(density,DENSITY,REAL);

      // loop over all elements of subdomain
      for (UInt actEl=0; actEl< elemsSD.GetSize(); actEl++) {

        BaseFE * ptElem = elemsSD[actEl]->ptElem;

        // get element coordinates
        StdVector<UInt> & connect = elemsSD[actEl]->connect;
        Matrix<Double> ptCoord;
        ptgrid_->GetElemNodesCoord( ptCoord, connect, false );

        // get shape function at center of the element
        Vector<Double> shapeFnc;
        Vector<Double> LCoord;
        ptElem -> GetCoordMidPoint(LCoord);
        ptElem -> GetShFnc(shapeFnc,LCoord);

        // retrieve 1st derivative and multiply with density, 
        //  since p = rho * dpsi/dt
        Vector<TYPE> valueElem;
        GetDerivSolVecOfElement(valueElem, connect);

        valueElem *= density;
        Vector<TYPE> pressureElem(1);
        pressureElem[0] = valueElem * shapeFnc;

        // map element result back in global set of results
        UInt pdeElem;
        pdeElem = eqnMap_->Mesh2PdeElem(elemsSD[actEl]->elemNum);

        if ( isComplex_ == true ) {
          ElemStoreSol<Complex> & pressure = 
            dynamic_cast<ElemStoreSol<Complex>&>(*acouPressure_);
          pressure.SetElemResult(pdeElem-1,pressureElem);
        }
        else {
          ElemStoreSol<Double> & pressure = 
            dynamic_cast<ElemStoreSol<Double>&>(*acouPressure_);
          pressure.SetElemResult(pdeElem-1,pressureElem);
        }

      }
    }
  }


  template <class TYPE>
  void AcousticPDE::CalcAcouPower()
  {
    ENTER_FCN( "ElecPDE::CalcAcouPower" );
    
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
    StdVector<SurfElem*> surfElems;

    UInt numSD = calcAcouPower_.GetSize();

    //some help variables
    Vector<Double> lCoordSurf, lCoordVol, normal;
    Vector<TYPE> gradVal(dim_);
    Vector<TYPE> sumOfPower(numSD);
    Vector<TYPE> elemIntensity(dim_);
    sumOfPower.Init();

    Elem * ptVolElem;
    BaseFE * ptSurfElemFE, * ptVolElemFE;

    TYPE gradNormal = 0.0;
    TYPE elemPower  = 0.0;
    Double normSign = 0;
    Double density  = 0.0;
    Double elemPowerReal;
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

    TYPE factorI;  

    // loop over all subdomains
    for (UInt iSD=0; iSD<calcAcouPower_.GetSize(); iSD++){
    
      // get surface and acoording volume elements
      ptgrid_->GetSurfElems( surfElems, calcAcouPower_[iSD] );
    
      // loop over all surface elements
      for (UInt iElem=0; iElem<surfElems.GetSize(); iElem++) {
        
	// Determine, which volume element is the right neighbour for the 
	// calculation;
	// our normal should point out of the correct neighbor volume element!
	if ( acouPowerNeighborRegion_.
	     Find(surfElems[iElem]->ptVolElem1->regionId) != -1 ) {
	  ptVolElem = surfElems[iElem]->ptVolElem1;
	  normSign = 1.0;
	} 
	else {
	  ptVolElem = surfElems[iElem]->ptVolElem2;
	  normSign = -1.0;
	}
	
	normSign *= (Double) surfElems[iElem]->normalSign;
	
	ptSurfElemFE = surfElems[iElem]->ptElem; 
	ptVolElemFE = ptVolElem->ptElem;
        
	const StdVector<UInt> & surfConnect = surfElems[iElem]->connect;
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
	gradOp->CalcElemGradField(gradVal, ptVolElem, 
				  lCoordVol,1.0);
	// Calc global normal
	ptgrid_->CalcSurfNormal(normal, *surfElems[iElem]);
	
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
	//pdeElemNum = eqnMap_->Mesh2PdeElem(surfElems[iElem]->elemNum);
        
	if ( formulation_ == ACOU_PRESSURE ) {
	  elemPower = multVal * gradNormal * elemSol / density;
	  factorI   = multVal * elemSol / density;
	  elemIntensity = factorI * gradVal;
	}
	else {
	  elemPower = multVal * density * gradNormal * elemSol;
	  factorI   = multVal * density * elemSol;
	  elemIntensity = factorI * gradVal;
	}
	
	//take the real part and multiply it with surface
	elemPower *= area;
	
	// Create temporary vector, since SetElemResult only
	// can handle these
	Vector<Complex> acouPowerVec(1);
	acouPowerVec[0] = elemPower;
	sumOfPower[iSD] += elemPower;
	acouPower_->SetElemResult(pdeElemNum-1, acouPowerVec);

	acouIntensity_->SetElemResult(pdeElemNum-1, elemIntensity);
      }
      
    }
    
    StdVector<std::string> subdomNames;
    ptgrid_->RegionIdToName( subdomNames,  calcAcouPower_ );

    Info->WriteAcouPower(pdename_, subdomNames, sumOfPower);

    delete gradOp;


  }


  void AcousticPDE::WriteResultsInFile(const UInt kstep,
                                       const Double asteptime,
                                       UInt stepOffset,
                                       Double timeOffset) {
    ENTER_FCN( "AcousticPDE::WriteResultsInFile" );

#ifdef PARALLEL //only one thread should write the output
    int commrank;
    MPI_Comm_rank(MPI_COMM_WORLD,&commrank);
    if (!commrank) {
#endif

      NodeStoreSol<Complex> * solHarmonic;
      NodeStoreSol<Double> * solTransient;

      NodeStoreSol<Complex> solD1Harmonic;

      Double actTime = asteptime + timeOffset;
      UInt actStep = kstep + stepOffset;
      
      if (analysistype_==HARMONIC) {


        solHarmonic = dynamic_cast<NodeStoreSol<Complex>*>(sol_);

        // Please remember, that in harmonic case actTime means indeed
        //  the actual frequency
        if (saveSol_)
          outFile_->WriteNodeSolutionHarmonic(*solHarmonic, actStep, 
                                              actTime, complexFormat_);

        if (saveDeriv1_) {

          Vector<Complex> solD1 = solHarmonic->GetAlgSysVector();

          Double omega = 2.0 * PI * actTime;
          Complex factor = Complex(0.0,omega); // factor = j*2*PI*f
          solD1 *= factor;


          solD1Harmonic.SetNumSolutions(1);
          solD1Harmonic.SetSolutionType(ACOU_POTENTIAL_DERIV_1);
          solD1Harmonic.SetNumDofs(1);
          solD1Harmonic.SetNumNodes(numPDENodes_);
          solD1Harmonic.SetPtrEQNData(eqnMap_.get(), ptgrid_);
          solD1Harmonic.SetRegions( subdoms_ );
          solD1Harmonic.SetResult( results_[0] );
          solD1Harmonic.Init();

          solD1Harmonic.SetAlgSysVector(solD1);
          outFile_->WriteNodeSolutionHarmonic(solD1Harmonic, actStep, 
                                              actTime, complexFormat_);
        }

        if (calcElemPressure_.GetSize() != 0 ) {
          ElemStoreSol<Complex> & pressureConverted = 
            dynamic_cast<ElemStoreSol<Complex>&>(*acouPressure_);
          outFile_->WriteElemSolutionHarmonic(pressureConverted, actStep,
                                              actTime, complexFormat_);
        }

	if ( calcAcouPower_.GetSize() > 0 ) {
	  ComplexFormat complexFormat = REAL_IMAG;

          ElemStoreSol<Complex> & acouPowerConverted = 
            dynamic_cast<ElemStoreSol<Complex>&>(*acouPower_);
          outFile_->WriteElemSolutionHarmonic(acouPowerConverted, actStep,
                                              actTime, complexFormat);

          ElemStoreSol<Complex> & acouIntensityConverted = 
            dynamic_cast<ElemStoreSol<Complex>&>(*acouIntensity_);
          outFile_->WriteElemSolutionHarmonic(acouIntensityConverted, actStep,
                                              actTime, complexFormat);
	}

      }
      else {

        solTransient = dynamic_cast<NodeStoreSol<Double>*>(sol_);

        if (saveSol_){
          outFile_->WriteNodeSolutionTransient(*solTransient,actStep,actTime);
        }

        if (saveDeriv1_) {
          solDeriv1_.SetAlgSysVector(getS1()); 
          outFile_->WriteNodeSolutionTransient(solDeriv1_, actStep, actTime);
        }

        if (saveDeriv2_) {
          solDeriv2_.SetAlgSysVector(getS2());
          outFile_->WriteNodeSolutionTransient(solDeriv2_, actStep, actTime);
        }


        if (saveRHSval_){
          outFile_->WriteNodeSolutionTransient(rhs_, actStep, actTime);
          if (plotRHSVel_==true)
            outFile_->WriteNodeSolutionTransient(rhs2_, actStep, actTime);
        }

        if (calcElemPressure_.GetSize() != 0 ) {
          ElemStoreSol<Double> & pressureConverted = 
            dynamic_cast<ElemStoreSol<Double>&>(*acouPressure_);
          outFile_->WriteElemSolutionTransient(pressureConverted,actStep,actTime);
        }

      }
 

#ifdef PARALLEL
    }//!commrank
#endif
  }

  void AcousticPDE::WriteHistoryInFile(const UInt kstep,
                                       const Double asteptime,
                                       UInt stepOffset,
                                       Double timeOffset) {
    ENTER_FCN( "AcousticPDE::WriteHistoryInFile" );

#ifdef PARALLEL //only one thread should write the output
    int commrank;
    MPI_Comm_rank(MPI_COMM_WORLD,&commrank);
    if (!commrank) {
#endif
      NodeStoreSol<Double> * solTransient;
      NodeStoreSol<Complex> * solHarmonic;
      
     
      Double actTime = asteptime + timeOffset;
      UInt actStep = kstep + stepOffset;
      
      if (analysistype_==HARMONIC) {
        solHarmonic = dynamic_cast<NodeStoreSol<Complex>*>(sol_);

        if (saveSolHist_)
          outFile_->WriteNodeHistoryHarmonic(*solHarmonic,  actStep, 
                                             actTime, complexFormat_);

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
          if (plotRHSVel_==true)
            outFile_->WriteNodeHistoryTransient(rhs2_, actStep, actTime); 
        }

        if (saveElemForceHist_.GetSize() != 0) {
          outFile_->WriteElemHistoryTransient(acouForce_, actStep, actTime);

          // lazy mans way to produce output
          std::cerr << actTime << "   " << sumAcouForce_ << std::endl;
        }

        if (saveElemPressureHist_.GetSize() != 0) {
          ElemStoreSol<Double> & pressureConverted = 
            dynamic_cast<ElemStoreSol<Double>&>(*acouPressure_);
          outFile_->WriteElemHistoryTransient(pressureConverted, actStep, 
                                              actTime);
        }

      }
 
#ifdef PARALLEL
    }//!commrank
#endif
  }


  // ***********************************************************************
  //   Obtain information on desired output quantities from parameter file
  // ***********************************************************************
  void AcousticPDE::ReadStoreResults() {  
    ENTER_FCN( "AcousticPDE::ReadStoreResults" );
    
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

    if ( formulation_ == ACOU_POTENTIAL ) {
      // --- pressure ---
      Enum2String(ACOU_PRESSURE, quantity);
      valVec = "", "", quantity;
      params->GetList( keyVec, attrVec, valVec, nodeValues);
      if (nodeValues.GetSize() > 0) {
        (*error) << "It makes no sense to have a PDE in acoustic potential "
                 << "and retrieve \n node results in pressure. "
                 << "Try element results!";
        Error( __FILE__, __LINE__ );
      }

      // --- acoustic potential ---
      Enum2String(ACOU_POTENTIAL, quantity);
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


      // --- acoustic potential, 1. Deriv ---
      Enum2String(ACOU_POTENTIAL_DERIV_1, quantity);
      valVec = "", "", quantity;
      params->GetList( keyVec, attrVec, valVec, nodeValues);
      if (nodeValues.GetSize() > 0) {

        saveDeriv1_ = true;
        hasOutput_ = true;
        if (analysistype_ != HARMONIC ) {

          // intialize corresponding storesolution object
          solDeriv1_.SetNumSolutions(1);
          solDeriv1_.SetNumNodes(numPDENodes_);
          solDeriv1_.SetSolutionType(ACOU_POTENTIAL_DERIV_1);
          solDeriv1_.SetNumDofs(1);
          solDeriv1_.SetResult( results_[0] );
          solDeriv1_.SetPtrEQNData(eqnMap_.get(), ptgrid_); 
          solDeriv1_.SetRegions( subdoms_ );
          solDeriv1_.Init(0);
        }
      }
      
      // --- acoustic potential, 2. Deriv ---
      Enum2String(ACOU_POTENTIAL_DERIV_2, quantity);
      valVec = "", "", quantity;
      params->GetList( keyVec, attrVec, valVec, nodeValues);
      if (nodeValues.GetSize() > 0 && analysistype_ != HARMONIC) {
        saveDeriv2_ = true;
        hasOutput_ = true;

        solDeriv2_.SetNumSolutions(1);
        solDeriv2_.SetNumNodes(numPDENodes_);
        solDeriv2_.SetSolutionType(ACOU_POTENTIAL_DERIV_2);
        solDeriv2_.SetNumDofs(1);
        solDeriv2_.SetResult( results_[0] );
        solDeriv2_.SetPtrEQNData(eqnMap_.get(),ptgrid_); 
        solDeriv2_.SetRegions( subdoms_ );
        solDeriv2_.Init(0);
      }
    }
    else if ( formulation_ == ACOU_PRESSURE ) {
      // --- acoustic potential ---
      Enum2String(ACOU_POTENTIAL, quantity);
      valVec = "", "", quantity;
      params->GetList( keyVec, attrVec, valVec, nodeValues);
      if (nodeValues.GetSize() > 0) {
        (*error) << "In acoustic pressure formulation, \n" 
                 << "retrieving acoustic potential requires integration.\n"
                 << "So forget it!";
        Error(__FILE__,__LINE__);
      }

      // --- pressure ---
      Enum2String(ACOU_PRESSURE, quantity);
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


      // --- acoustic pressure, 1. Deriv ---
      Enum2String(ACOU_PRESSURE_DERIV_1, quantity);
      valVec = "", "", quantity;
      params->GetList( keyVec, attrVec, valVec, nodeValues);
      if (nodeValues.GetSize() > 0) {

        saveDeriv1_ = true;
        hasOutput_ = true;
        if (analysistype_ != HARMONIC ) {

          // intialize corresponding storesolution object
          solDeriv1_.SetNumSolutions(1);
          solDeriv1_.SetNumNodes(numPDENodes_);
          solDeriv1_.SetSolutionType(ACOU_PRESSURE_DERIV_1);
          solDeriv1_.SetNumDofs(1);
          solDeriv1_.SetResult( results_[0] );
          solDeriv1_.SetPtrEQNData(eqnMap_.get(), ptgrid_); 
          solDeriv1_.SetRegions( subdoms_ );
          solDeriv1_.Init(0);
        }
      }
      
      // --- acoustic pressure, 2. Deriv ---
      Enum2String(ACOU_PRESSURE_DERIV_2, quantity);
      valVec = "", "", quantity;
      params->GetList( keyVec, attrVec, valVec, nodeValues);
      if (nodeValues.GetSize() > 0 && analysistype_ != HARMONIC) {
        saveDeriv2_ = true;
        hasOutput_ = true;

        solDeriv2_.SetNumSolutions(1);
        solDeriv2_.SetNumNodes(numPDENodes_);
        solDeriv2_.SetSolutionType(ACOU_PRESSURE_DERIV_2);
        solDeriv2_.SetNumDofs(1);
        solDeriv2_.SetResult( results_[0] );
        solDeriv2_.SetPtrEQNData(eqnMap_.get(),ptgrid_); 
        solDeriv2_.SetRegions( subdoms_ );
        solDeriv2_.Init(0);
      }
    }
    
    
    // *****************************
    // Determine element results
    // *****************************
    StdVector<std::string> regionNames;
    keyVec  = pdename_, "storeResults", "elemResults", "region";
    attrVec = "", "", "type";

    if ( formulation_ == ACOU_POTENTIAL ) {
      // pressure
      Enum2String(ACOU_PRESSURE, quantity);
      valVec = "", "", quantity;
      params->GetList( keyVec, attrVec, valVec, regionNames);
      ptgrid_->RegionNameToId( calcElemPressure_, regionNames );

      // If the the symbolic name is "all" compute quatity for all regions
      if ( calcElemPressure_.GetSize() == 1 &&
           calcElemPressure_[0] == ALL_REGIONS ) {
        calcElemPressure_ = subdoms_;
      }

      if (calcElemPressure_.GetSize() > 0) {

        hasOutput_ = true;
        Info->PrintF( pdename_, " Computing acoustic pressure for regions:\n" );
        for ( UInt k = 0; k < regionNames.GetSize(); k++ ) {
          Info->PrintF( pdename_, " %s\n", regionNames[k].c_str() );
        }

        if (analysistype_ == HARMONIC ) {
          acouPressure_ = new ElemStoreSol<Complex>;
        } 
        else {
          acouPressure_ = new ElemStoreSol<Double>;
        }

        acouPressure_->SetNumSolutions(1);
        acouPressure_->SetSolutionType(ACOU_PRESSURE);
        acouPressure_->SetNumElems(numElems_);
        acouPressure_->SetNumDofs(1);
        acouPressure_->SetResult( results_[0] );
        acouPressure_->SetPtrEQNData(eqnMap_.get(), ptgrid_);
        acouPressure_->Init();
      }
    } 


    // *****************************
    // Determine nodal history
    // *****************************
    StdVector<std::string> saveNodeHist; 
    keyVec  = pdename_, "storeResults", "nodeHistory", "saveNodes";
    attrVec = "", "", "type";

    if ( formulation_ == ACOU_POTENTIAL ) {
      // pressure
      Enum2String(ACOU_PRESSURE, quantity);
      valVec = "", "", quantity;
      params->GetList( keyVec, attrVec, valVec, saveNodeHist);
      if (saveNodeHist.GetSize() > 0) {
        (*error) <<  "In acoustic potential formulation, \n"
                 << " the quantity pressure is a element result.";
        Error( __FILE__, __LINE__ );
      }

      // --- acoustic potential  ---
      Enum2String(ACOU_POTENTIAL, quantity);
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

      // --- acoustic potential, 1. Deriv ---
      Enum2String(ACOU_POTENTIAL_DERIV_1, quantity);
      valVec  = "", "", quantity;
      params->GetList( keyVec, attrVec, valVec, saveNodeHist );
      if (saveNodeHist.GetSize() > 0) {
        saveDeriv1Hist_ = true;
        hasOutput_ = true;
        Info->PrintF( pdename_, "Saving acouPotentialD1 for Nodes:\n" );
        for ( UInt k = 0; k < saveNodeHist.GetSize(); k++ ) {
          Info->PrintF( pdename_, "  %s\n", saveNodeHist[k].c_str() );
        }
        // Check if solDeriv_1 was already set by the previous nodal
        // output check
        if( saveDeriv1_ != true ) {
          // intialize corresponding storesolution object
          solDeriv1_.SetNumSolutions(1);
          solDeriv1_.SetNumNodes(numPDENodes_);
          solDeriv1_.SetSolutionType(ACOU_POTENTIAL_DERIV_1);
          solDeriv1_.SetNumDofs(1);
          solDeriv1_.SetResult( results_[0] );
          solDeriv1_.SetPtrEQNData(eqnMap_.get(), ptgrid_); 
          solDeriv1_.Init(0);
        }
      }

      // --- acoustic potential, 2. Deriv ---
      Enum2String(ACOU_POTENTIAL_DERIV_2, quantity);
      valVec  = "", "", quantity;
      params->GetList( keyVec, attrVec, valVec, saveNodeHist );
    
      if (saveNodeHist.GetSize() > 0) {
        saveDeriv2Hist_ = true;
        hasOutput_ = true;
        Info->PrintF( pdename_, "Saving acouPotetentialD2 for Nodes:\n" );
        for ( UInt k = 0; k < saveNodeHist.GetSize(); k++ ) {
          Info->PrintF( pdename_, "  %s\n", saveNodeHist[k].c_str() );
        }
        // Check if solDeriv_2 was already set by the previous nodal
        // output check
        if( saveDeriv2_ != true ) {
          solDeriv2_.SetNumSolutions(1);
          solDeriv2_.SetNumNodes(numPDENodes_);
          solDeriv2_.SetSolutionType(ACOU_POTENTIAL_DERIV_2);
          solDeriv2_.SetNumDofs(1);
          solDeriv2_.SetResult( results_[0] );
          solDeriv2_.SetPtrEQNData(eqnMap_.get(),ptgrid_); 
          solDeriv2_.Init(0);
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
      // --- acoustic potential ---
      Enum2String(ACOU_POTENTIAL, quantity);
      valVec = "", "", quantity;
      params->GetList( keyVec, attrVec, valVec, nodeValues);
      if (nodeValues.GetSize() > 0) {
        (*error) << "It acoustic pressure formulation, \n" 
                 << "retrieving acoustic potential requires integration.\n"
                 << "So forget it!";
        Error(__FILE__,__LINE__);
      }

      // --- pressure ---
      Enum2String(ACOU_PRESSURE, quantity);
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

      // --- acoustic pressure, 1. Deriv ---
      Enum2String(ACOU_PRESSURE_DERIV_1, quantity);
      valVec  = "", "", quantity;
      params->GetList( keyVec, attrVec, valVec, saveNodeHist );
      if (saveNodeHist.GetSize() > 0) {
        saveDeriv1Hist_ = true;
        hasOutput_ = true;
        Info->PrintF( pdename_, "Saving acouPressureD1 for Nodes:\n" );
        for ( UInt k = 0; k < saveNodeHist.GetSize(); k++ ) {
          Info->PrintF( pdename_, "  %s\n", saveNodeHist[k].c_str() );
        }
        // Check if solDeriv_1 was already set by the previous nodal
        // output check
        if( saveDeriv1_ != true ) {
          // intialize corresponding storesolution object
          solDeriv1_.SetNumSolutions(1);
          solDeriv1_.SetNumNodes(numPDENodes_);
          solDeriv1_.SetSolutionType(ACOU_PRESSURE_DERIV_1);
          solDeriv1_.SetNumDofs(1);
          solDeriv1_.SetResult( results_[0] );
          solDeriv1_.SetPtrEQNData(eqnMap_.get(), ptgrid_); 
          solDeriv1_.Init(0);
        }
      }

      // --- acoustic potential, 2. Deriv ---
      Enum2String(ACOU_PRESSURE_DERIV_2, quantity);
      valVec  = "", "", quantity;
      params->GetList( keyVec, attrVec, valVec, saveNodeHist );
    
      if (saveNodeHist.GetSize() > 0) {
        saveDeriv2Hist_ = true;
        hasOutput_ = true;
        Info->PrintF( pdename_, "Saving acouPressureD2 for Nodes:\n" );
        for ( UInt k = 0; k < saveNodeHist.GetSize(); k++ ) {
          Info->PrintF( pdename_, "  %s\n", saveNodeHist[k].c_str() );
        }
        // Check if solDeriv_2 was already set by the previous nodal
        // output check
        if( saveDeriv2_ != true ) {
          solDeriv2_.SetNumSolutions(1);
          solDeriv2_.SetNumNodes(numPDENodes_);
          solDeriv2_.SetSolutionType(ACOU_PRESSURE_DERIV_2);
          solDeriv2_.SetNumDofs(1);
          solDeriv2_.SetResult( results_[0] );
          solDeriv2_.SetPtrEQNData(eqnMap_.get(),ptgrid_); 
          solDeriv2_.Init(0);
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

    // *****************************
    // Determine element history
    // *****************************
    keyVec  = pdename_, "storeResults", "elemHistory", "saveElems";
    attrVec = "", "", "type";

    if ( (formulation_==ACOU_POTENTIAL) || (formulation_==ACOU_PRESSURE) ) {

      // --- force ---
      Enum2String(ACOU_FORCE, quantity);
      valVec  = "", "", quantity;
      params->GetList( keyVec, attrVec, valVec, saveElemForceHist_ );

      if (saveElemForceHist_.GetSize() > 0) {
        hasOutput_ = true;
        acouForce_.SetNumSolutions(1);
        acouForce_.SetSolutionType(ACOU_FORCE);
        acouForce_.SetNumElems(numElems_);
        acouForce_.SetNumDofs(1);
        acouForce_.SetResult( results_[0] );
        acouForce_.SetPtrEQNData(eqnMap_.get(), ptgrid_);
        acouForce_.Init();

        Info->PrintF( pdename_, "Saving acouForce for Elements:\n" );
        for ( UInt k = 0; k < saveElemForceHist_.GetSize(); k++ ) {
          Info->PrintF( pdename_, "  %s\n", saveElemForceHist_[k].c_str() );
        }
      }

  
      // --- acoustic power ---
      StdVector<std::string> temp;
      params->GetList( "region", regionNames, pdename_, "acouIntensityPower" );
      params->GetList( "element", temp, pdename_, "acouIntensityPower" );

      ptgrid_->RegionNameToId( acouPowerNeighborRegion_, regionNames );
      ptgrid_->RegionNameToId( calcAcouPower_, temp );

      if (calcAcouPower_.GetSize() > 0) {
        if ( analysistype_ == HARMONIC ) {
	  hasOutput_ = true;
	  Info->PrintF( pdename_,
			" Computing acoustic intensity/power for regions:\n");
	  for ( UInt k = 0; k < temp.GetSize(); k++ ) {
	    Info->PrintF( pdename_, " %s\n", temp[k].c_str() );
	  }

	  acouPower_     = new ElemStoreSol<Complex>;
	  acouIntensity_ = new ElemStoreSol<Complex>;

	  // Resize solution arrays
	  acouPower_->SetNumSolutions(1);
	  acouPower_->SetSolutionType(ACOU_POWER);
	  acouPower_->SetNumElems(numElems_);
	  acouPower_->SetNumDofs(1);
	  acouPower_->SetPtrEQNData( eqnMap_.get(), ptgrid_);
	  acouPower_->Init();

	  acouIntensity_->SetNumSolutions(1);
	  acouIntensity_->SetSolutionType(ACOU_INTENSITY);
	  acouIntensity_->SetNumElems(numElems_);
	  acouIntensity_->SetNumDofs(dim_);
	  acouIntensity_->SetPtrEQNData( eqnMap_.get(), ptgrid_);
	  acouIntensity_->Init();
	}
	else {
          (*warning) << "Acoustic power computation just supported for harmonic analysis.\n";
	  Warning( __FILE__, __LINE__ );
        }
      } 
      
    }
    if ( formulation_==ACOU_POTENTIAL ) {

      // --- pressure ---
      Enum2String(ACOU_PRESSURE, quantity);
      valVec  = "", "", quantity;
      params->GetList( keyVec, attrVec, valVec, saveElemPressureHist_ );

      if (saveElemPressureHist_.GetSize() > 0) {
        hasOutput_ = true;
        acouPressure_->SetNumSolutions(1);
        acouPressure_->SetSolutionType(ACOU_PRESSURE);
        acouPressure_->SetNumElems(numElems_);
        acouPressure_->SetNumDofs(1);
        acouPressure_->SetResult( results_[0] );
        acouPressure_->SetPtrEQNData(eqnMap_.get(),ptgrid_);
        acouPressure_->Init();

        Info->PrintF( pdename_, "Saving acouPressure for Elements:\n" );
        for ( UInt k = 0; k < saveElemPressureHist_.GetSize(); k++ ) {
          Info->PrintF( pdename_, "  %s\n", saveElemPressureHist_[k].c_str() );
        }
      }
    }


  }

  // ***********************************************************************
  //   Obtain information on desired output quantities from parameter file
  // ***********************************************************************
  void AcousticPDE::ReadDataPML(std::string& dampingTypePML, 
                                Matrix<Double>& inner, 
                                Double& dampPML, 
                                RegionIdType actRegion) {
  
    ENTER_FCN( "AcousticPDE::ReadDataPML" );

    inner.Resize(dim_,dim_);
    inner.Init();

    // help variables for parameter checking
    StdVector<std::string> keyVec;
    StdVector<std::string> attrVec;
    StdVector<std::string> valVec;

    StdVector<std::string> stringVal;

    // Construct vectors for restricted parameter search
    std::string actRegionName;
    actRegionName = ptgrid_->RegionIdToName( actRegion );
    keyVec = "acoustic" , "region" , "pml" , "type";
    attrVec= ""         , "name"   , "";
    valVec = ""         , actRegionName, "";
    params->GetList( keyVec, attrVec, valVec, stringVal);
    dampingTypePML = stringVal[0];

    StdVector<Double> val;
 
    //xMin
    keyVec = "acoustic" , "region" , "pml" , "xMin";
    params->GetList( keyVec, attrVec, valVec, val);
    inner[0][0] =  val[0];

    //xMax
    keyVec = "acoustic" , "region" , "pml" , "xMax";
    params->GetList( keyVec, attrVec, valVec, val);
    inner[1][0] =  val[0];

    //yMin
    keyVec = "acoustic" , "region" , "pml" , "yMin";
    params->GetList( keyVec, attrVec, valVec, val);
    inner[0][1] =  val[0];

    //yMax
    keyVec = "acoustic" , "region" , "pml" , "yMax";
    params->GetList( keyVec, attrVec, valVec, val);
    inner[1][1] =  val[0];

    if ( dim_ == 3 ) {
      //zMin
      keyVec = "acoustic" , "region" , "pml" , "zMin";
      params->GetList( keyVec, attrVec, valVec, val);
      inner[0][2] =  val[0];

      //zMax     
      keyVec = "acoustic" , "region" , "pml" , "zMax";
      params->GetList( keyVec, attrVec, valVec, val);
      inner[1][2] =  val[0];
    }

    keyVec = "acoustic" , "region" , "pml" , "dampFactor";
    params->GetList( keyVec, attrVec, valVec, val);
    dampPML = val[0];


  }


  // ***********************************************************************
  //   Obtain information on desired output quantities from parameter file
  // ***********************************************************************
  void AcousticPDE::GetPMLLayerData(Matrix<Double>& inner, 
                                    Matrix<Double>& outer,
                                    UInt actSD)  {  

    ENTER_FCN( "AcousticPDE::GetPMLLayerData" );
    
    outer.Resize(inner.GetSizeRow(),inner.GetSizeCol());

    outer[0][0] = outer[1][0] = inner[1][0];
    outer[0][1] = outer[1][1] = inner[1][1];
    if (inner.GetSizeCol() > 2 ) {
      outer[0][2] = outer[1][2] = inner[1][2];
    }

    // inner/outer:   xmin  ymin  zmin
    //                xmax  ymax  zmax


    StdVector<Elem*> elemssd;
    ptgrid_->GetElems(elemssd, subdoms_[actSD]);

    for (UInt actEl=0; actEl< elemssd.GetSize(); actEl++) {
      StdVector<UInt> connecth = elemssd[actEl]->connect;
                             
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


  void AcousticPDE::ReadDataDampLayer(std::string& dampingTypePML, 
                                      Vector<Double>& mPoint, 
                                      Double& dampFactor, 
                                      Double& dampFactorMax, 
                                      Double& startRadius, 
                                      Double& endRadius, 
                                      RegionIdType actRegion) {
  
    ENTER_FCN( "AcousticPDE::ReadDataDampLayer" );

    mPoint.Resize(dim_);

    // help variables for parameter checking
    StdVector<std::string> keyVec;
    StdVector<std::string> attrVec;
    StdVector<std::string> valVec;

    StdVector<std::string> stringVal;

    // Construct vectors for restricted parameter search
    std::string actRegionName;
    actRegionName = ptgrid_->RegionIdToName( actRegion );
    keyVec = "acoustic" , "region" , "dampLayer" , "type";
    attrVec= ""         , "name"   , "";
    valVec = ""         , actRegionName, "";
    params->GetList( keyVec, attrVec, valVec, stringVal);
    dampingTypePML = stringVal[0];

    StdVector<Double> val;
 
    //get the center point

    //xM
    keyVec = "acoustic" , "region" , "dampLayer" , "xM";
    params->GetList( keyVec, attrVec, valVec, val);
    mPoint[0] =  val[0];

    //yM
    keyVec = "acoustic" , "region" , "dampLayer" , "yM";
    params->GetList( keyVec, attrVec, valVec, val);
    mPoint[1] =  val[0];

    if ( dim_ == 3 ) {
      //zM
      keyVec = "acoustic" , "region" , "dampLayer" , "zM";
      params->GetList( keyVec, attrVec, valVec, val);
      mPoint[2] =  val[0];
    }

    //start radius
    keyVec = "acoustic" , "region" , "dampLayer" , "radiusStart";
    params->GetList( keyVec, attrVec, valVec, val);
    startRadius =  val[0];

    //end radius
    keyVec = "acoustic" , "region" , "dampLayer" , "radiusEnd";
    params->GetList( keyVec, attrVec, valVec, val);
    endRadius =  val[0];

    //get damping factor
    keyVec = "acoustic" , "region" , "dampLayer" , "dampFactor";
    params->GetList( keyVec, attrVec, valVec, val);
    dampFactor = val[0];

    //get maximal damping factor (saturation)
    keyVec = "acoustic" , "region" , "dampLayer" , "dampFactorMax";
    params->GetList( keyVec, attrVec, valVec, val);
    dampFactorMax = val[0];

  }

} // end of namespace
