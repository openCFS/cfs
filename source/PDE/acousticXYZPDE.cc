// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <fstream>
#include <iostream>
#include <string>
#include <math.h>

#include "acousticXYZPDE.hh" 
#include "Forms/forms_header.hh"
#include "Forms/abcInt.hh"
#include "Forms/linNeumannInt.hh"
#include "Estimator/spaceerror.hh"
#include "newmark.hh"
#include "newmarkFracDamp.hh"
#include "DataInOut/WriteInfo.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "Utils/mathfunctions.hh"
#include "Utils/nodestoresol.hh"
#include "Driver/solveStepAcoustic.hh"
#include "CoupledPDE/pdecoupling.hh"
#include "Driver/assemble.hh"

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
  AcousticXYZPDE::AcousticXYZPDE( Grid * aptgrid, ParamNode* paramNode )
    :SinglePDE( aptgrid, paramNode ) {


    // ****************************
    // DETERMINE GEOMETRY
    // ****************************


    pdename_          = "acousticXYZ";
    pdematerialclass_ = FLUID;

  }


  // *********************************************
  //   Check what type of damping should be used
  // *********************************************
  void AcousticXYZPDE::ReadDampingInformation( ) {

    
    // Get region parameter nodes
    StdVector<ParamNode*> regionNodes = myParam_->Get("regionList")->GetList( "region" );

    for (UInt k = 0; k < regionNodes.GetSize(); k++) {

      // get current node and dtermine regionName/Id
      ParamNode * actNode = regionNodes[k];
      std::string actRegionName = actNode->Get("name")->AsString();
      RegionIdType actRegionId = ptgrid_->RegionNameToId( actRegionName );

      // try to get node with damping information
      ParamNode * dampNode = actNode->Get( "damping", false );
      
      // check, if damping exists at all
      if( !dampNode ) {
        dampingList_[actRegionId] = NONE;
        Info->PrintF( pdename_, 
                      "No information specifying damping detected!\n" );
        continue;
      }
      
       // get specific damping type
      std::string dampString = dampNode->Get("type")->AsString();

    
      if (dampString  == "no") {
        dampingList_[actRegionId] = NONE;
        Info->PrintF( pdename_, 
                      "      * NO damping at all for region: %s\n", 
                      actRegionName.c_str() );;
      }
      else if ( dampString == "rayleigh") {
        dampingList_[actRegionId] = RAYLEIGH;
        Info->PrintF( pdename_, 
                      "      * RAYLEIGH damping for region: %s\n",
                      actRegionName.c_str() );
      }
    }

  }


  void AcousticXYZPDE::DefineIntegrators() {


    Double density, compressibility, c0, alpha, beta, BoverA;
    Double coeffmass, coeffdamp, coeffstiff;

    for (UInt actSD = 0; actSD < subdoms_.GetSize(); actSD++) {

      BaseMaterial * actMat = materials_[subdoms_[actSD]];

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
      //   In AcousticXYZPDE ALL integrators are multiplied with density!
      // ********************************************************************

      actMat->GetScalar(density,DENSITY,REAL);
      actMat->GetScalar(compressibility,ACOU_BULK_MODULUS,REAL);

      c0 = sqrt(compressibility/density);

      actMat->GetScalar(alpha,RAYLEIGH_ALPHA,REAL);
      actMat->GetScalar(beta,RAYLEIGH_BETA,REAL);
      actMat->GetScalar(BoverA,BOVERA,REAL);

      // ********************************************************************
      //   Linear wave equation
      // ********************************************************************


      //do the stiffness part;

      //check  for PML!
      ParamNode * pmlNode = actRegionNode->Get("pml", false );
      

        if ( pmlNode ) {
	//read data for PML layer

	//type of PML damping
	std::string dampingTypePML;

	// inner / outer region
	Matrix<Double> inner;
	Matrix<Double> outer;

	//damping factor
	Double dampPML;

	ReadDataPML(dampingTypePML, inner, dampPML, pmlNode );
	dampPML *= c0;

	GetPMLLayerData(inner, outer, actRegion);

	//====================================================================
	//	 stiffness integrator for PML
	//====================================================================

	std::string formsType = "stiffness";
	coeffstiff = 1.0 / (c0*c0);
	BaseForm * bilinearStiffPML = 
	  new PMLTransXYZInt(formsType, coeffstiff, dampingTypePML, dampPML, 
			     dim_, isaxi_);

	bilinearStiffPML->SetPosPML(inner,outer);

	BiLinFormContext * stiffContextPML = 
	  new BiLinFormContext(bilinearStiffPML, STIFFNESS );

	stiffContextPML->SetPtPdes( this, this );
        stiffContextPML->SetResults( results_[0], results_[0],
                                     actSDList, actSDList );
	assemble_->AddBiLinearForm( stiffContextPML );

	//====================================================================
	//	 damping integrator for PML
	//====================================================================
	formsType = "damping";
	coeffdamp = 1.0 / (c0*c0);
	BaseForm * bilinearDampPML = 
	  new PMLTransXYZInt(formsType, coeffdamp, dampingTypePML, dampPML, 
			     dim_, isaxi_);

	bilinearDampPML->SetPosPML(inner,outer);

	BiLinFormContext * dampContextPML = 
	  new BiLinFormContext(bilinearDampPML, DAMPING );

	dampContextPML->SetPtPdes( this, this );
        dampContextPML->SetResults( results_[0], results_[0],
                                    actSDList, actSDList );

	assemble_->AddBiLinearForm( dampContextPML );
      }
      
      //      else {
      
      // stiffness integrator 
      coeffstiff = 1.0;
      BaseForm * bilinearStiff = new LaplaceXYZInt(coeffstiff,dim_,isaxi_);        
      BiLinFormContext * stiffContext = 
        new BiLinFormContext( bilinearStiff, STIFFNESS );
      
      stiffContext->SetPtPdes( this, this );
      stiffContext->SetResults( results_[0], results_[0],
                                actSDList, actSDList );
      
      // mass integrator
      coeffmass = 1.0 / (c0*c0);
      
      BaseForm * bilinearMass  = new MassInt(coeffmass, dim_, isaxi_);
      BiLinFormContext * massContext = 
	new BiLinFormContext(bilinearMass, MASS );
      
      massContext->SetPtPdes(this, this);
      massContext->SetResults( results_[0], results_[0],
                               actSDList, actSDList );
      
      
      // ********************************************************************
      //   Additional terms for damping
      // ********************************************************************
      
      if ( dampingList_.size() > 0 ) {
	(*warning) << "NO Rayleigh damping supported!";
	Warning(__FILE__, __LINE__);  
      }
      
      // Finally add the stiffness/mass integrators
      assemble_->AddBiLinearForm( stiffContext );
      assemble_->AddBiLinearForm( massContext ); 
      //      }
    }
  
  }

  void AcousticXYZPDE::DefineSolveStep() {
   
    solveStep_ = new SolveStepAcoustic(*this);
   
  }


  // ======================================================
  // TIME STEPPING SECTION
  // ======================================================

  void AcousticXYZPDE::InitTimeStepping() {
    
    // timestepping formulation
    ParamNode* myLinSysNode = FindLinearSystem( pdename_ );
    
    // <system name="acoustic"/> exists
    if( myLinSysNode ) {

      std::string str = "";
      myLinSysNode->Get( "timeSteppingFormulation", str,  false );
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
    }

    if ( effectiveMass_ == false ) {
      TS_alg_ = new Newmark( algsys_, pdename_ );
    }
    else {
      TS_alg_ = new NewmarkEffMass( algsys_ , pdename_, false );
    }

  }



  // ======================================================
  // POSTPROCESSING SECTION
  // ======================================================

  void AcousticXYZPDE::CalcResults( shared_ptr<BaseResult> res ) {

    EXCEPTION( "Not implemented -> Blame Manfred and Andreas ;_)" );

  }

  void AcousticXYZPDE::DefineAvailResults() {

    // Get problem geometry and PDE subtype
    std::string probGeo = param->Get("domain")->Get("geometryType")->AsString();

    // === PRESSURE XYZ ===
    shared_ptr<ResultInfo> res1(new ResultInfo);
    shared_ptr<AnsatzFct> fct(new LagrangeFct);
    
    // Set number of degrees of freedom and
    // ensure that subtype fits to problem geometry    
    if ( probGeo == "3d" ) {
      res1->dofNames = "x", "y", "z";
      Info->PrintF("", "=== 3D PROBLEM\n");
    }
    else if (probGeo == "axi" ) {
      isaxi_ = true;
      res1->dofNames = "r", "z";
      Info->PrintF("", "=== AXISYSMMETRIC PROBLEM\n");
    }
    else if ( probGeo == "plane" ) {
      res1->dofNames = "x", "y";
      Info->PrintF("", "=== PLANE PROBLEM\n");
    }
    res1->unit = "Pa";
    res1->resultType = ACOU_PRESSUREXYZ;
    res1->definedOn = ResultInfo::NODE;
    res1->entryType = ResultInfo::VECTOR;
    res1->fctType = fct;
    
    results_.Push_back( res1 );
    availResults_.insert( res1 );

    // === STANDARD PRESSURE ===
    shared_ptr<ResultInfo> pres (new ResultInfo );
    pres->resultType = ACOU_PRESSURE;
    pres->dofNames = "";
    pres->unit = "Pa";
    pres->definedOn = ResultInfo::NODE;
    pres->entryType = ResultInfo::SCALAR;
    pres->fctType = res1->fctType;
    availResults_.insert( pres );
  }



  // ***********************************************************************
  //   Obtain information on desired output quantities from parameter file
  // ***********************************************************************
  void AcousticXYZPDE::ReadDataPML(std::string& dampingTypePML, 
                                Matrix<Double>& inner, 
                                Double& dampPML, 
                                ParamNode * actNode ) {
  

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
  void AcousticXYZPDE::GetPMLLayerData(Matrix<Double>& inner, 
                                    Matrix<Double>& outer,
                                    RegionIdType actRegion )  {  


    // inner/outer:   xmin  ymin  zmin
    //                xmax  ymax  zmax
    
    if ( inner.GetSizeCol() != dim_ ) {
      //we have to compute it, since the user has not specified it
      
      inner.Resize(2,dim_);
      inner.Init();
      
      std::cerr << "actREgion is " << actRegion << std::endl;
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

 

} // end of namespace
