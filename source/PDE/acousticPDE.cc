// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <math.h>
#include <stddef.h>
#include <algorithm>
#include <complex>
#include <iostream>
#include <string>
#include <utility>

#include "boost/tokenizer.hpp"

#include "CoupledPDE/pdecoupling.hh"
#include "DataInOut/Logging/cfslog.hh"
#include "DataInOut/Logging/log.hpp"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/WriteInfo.hh"
#include "Domain/ansatzFct.hh"
#include "Domain/domain.hh"
#include "Domain/elem.hh"
#include "Domain/entityList.hh"
#include "Domain/grid.hh"
#include "Domain/resultInfo.hh"
#include "Domain/surfElem.hh"
#include "Driver/assemble.hh"
#include "Driver/formsContext.hh"
#include "Driver/solveStepAcoustic.hh"
#include "Elements/basefe.hh"
#include "Forms/PMLInt.hh"
#include "Forms/abcInt.hh"
#include "Forms/acouPowerDensityOp.hh"
#include "Forms/acouRHSLinForm.hh"
#include "Forms/baseForm.hh"
#include "Forms/gradfieldop.hh"
#include "Forms/laplaceInt.hh"
#include "Forms/linNeumannInt.hh"
#include "Forms/linSurfForm.hh"
#include "Forms/massInt.hh"
#include "Forms/nonConformingInt.hh"
#include "Forms/pierceInt.hh"
#include "General/Enum.hh"
#include "General/exception.hh"
#include "MatVec/SingleVector.hh"
#include "MatVec/exprt/xpr1.hh"
#include "MatVec/exprt/xpr2.hh"
#include "MatVec/matrix.hh"
#include "Materials/baseMaterial.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "PDE/SinglePDE.hh"
#include "PDE/StdPDE.hh"
#include "PDE/basePDE.hh"
#include "PDE/eqnMap.hh"
#include "PDE/timestepping.hh"
#include "Utils/basenodestoresol.hh"
#include "Utils/mathfunctions.hh"
#include "Utils/nodestoresol.hh"
#include "Utils/result.hh"
#include "acousticPDE.hh"
#include "math.h"
#include "newmark.hh"
#include "newmarkFracDamp.hh"

namespace CoupledField {

  DECLARE_LOG(acoupde)
    DEFINE_LOG(acoupde, "pde.acoustic")


  // =========================================================================
  // set solution information
  // =========================================================================
    AcousticPDE::AcousticPDE( Grid* aptgrid, PtrParamNode paramNode )
      : SinglePDE( aptgrid, paramNode ) {


      pdename_          = "acoustic";
      pdematerialclass_ = FLUID;
      maxTimeDerivOrder_ = 2;
      fracDamping_ = false;

      isAPML_ = false;
      isMechCoupled_ = false;
      variableSpeedOfSoundCN_ = false;

      justInterpolate_ = false;

      mHandle_ =  domain->GetMathParser()->GetNewHandle();
      //      MathParser * mParser =  domain->GetMathParser();
      //      mParser->SetExpr( mHandle_, "t" );

      // PDE formulation either in acoustic potential or pressure
      std::string str = myParam_->Get("formulation")->As<std::string>();
      formulation_ = SolutionTypeEnum.Parse(str);
      str = "Using * " + str + " as state variable in formulation of PDE\n";
      Info->PrintF( pdename_, str.c_str() );

      //compute on deformed geometry
      str = myParam_->Get("updatedLagrange")->As<std::string>();
      if ( str == "yes" )
      {
        str = "Compute acoustic field on deformed geometry\n";
        Info->PrintF( pdename_, str.c_str() );
        updatedLagrangeForm_ = true;
      }

      // get geometry type
      std::string probGeo;
      param->Get("domain")->GetValue("geometryType", probGeo );

      // since we have no special subType as in mechanics, 
      // the subType 9 is equal to the geometry type
      subType_ = probGeo;

      // TODO: this was never initialized up to now! Don't know if this is OK. 
      saveNodalSourcesRHS_ = false;
    }


  // *********************************************
  //   Check what type of damping should be used
  // *********************************************
  void AcousticPDE::ReadDampingInformation( ) {


    fracMemory_ = 0;
    Integer firstFrac=-1;
    std::map<std::string, DampingType> idDampType;
    std::map<std::string, shared_ptr<RaylDampingData> > idRaylData;
    
    // try to get dampingList
    PtrParamNode dampListNode = myParam_->Get( "dampingList", ParamNode::PASS );
    if( dampListNode ) {

      // get specific damping nodes
      ParamNodeList dampNodes = dampListNode->GetChildren();

      for( UInt i = 0; i < dampNodes.GetSize(); i++ ) {

        std::string dampString = dampNodes[i]->GetName();
        std::string actId = dampNodes[i]->Get("id")->As<std::string>();

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
          dampNodes[i]->GetValue( "fracAlg", fracAlg, ParamNode::PASS);

          std::string interpol = "no";
          dampNodes[i]->GetValue( "interpolation", interpol, ParamNode::PASS );

          UInt fracMem = 1;
          dampNodes[i]->GetValue( "fracMemory", fracMem, ParamNode::PASS );

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

        else if( actType == RAYLEIGH ) {
          // set data for Rayleigh damping
          shared_ptr<RaylDampingData> actRaylDamp(new RaylDampingData());
          actRaylDamp->alpha = "0.0";
          actRaylDamp->beta = "0.0";
          actRaylDamp->adjustDamping = true;
          actRaylDamp->ratioDeltaF = 0.01;
          actRaylDamp->freq = 0.0;
          dampNodes[i]->GetValue( "freq", actRaylDamp->freq, ParamNode::PASS);
          dampNodes[i]->GetValue( "ratioDeltaF", actRaylDamp->ratioDeltaF,
                             ParamNode::PASS );
          dampNodes[i]->GetValue( "adjustDamping",actRaylDamp->adjustDamping,  
                                  ParamNode::PASS);
          idRaylData[actId] = actRaylDamp;        
        }

        // store damping type string
        idDampType[actId] = actType;

      }
      if ( fracDamping_== true ) {
        Info->PrintF(pdename_, "Memory size for fractional damping  is: %d\n",
                     fracMemory_ );
      }
    }

    // Run over all region and set entry in "regionNonLinId"
    ParamNodeList regionNodes =
      myParam_->Get("regionList")->GetChildren();

    RegionIdType actRegionId;
    std::string actRegionName, actDampingId;

    if( regionNodes.GetSize() > 0 ) {
      Info->PrintF( pdename_, "Damping in following region(s)\n" );
    }

    for (UInt k = 0; k < regionNodes.GetSize(); k++) {
      regionNodes[k]->GetValue( "name", actRegionName );
      regionNodes[k]->GetValue( "dampingId", actDampingId );
      if( actDampingId == "" )
        continue;

      actRegionId = ptgrid_->GetRegion().Parse(actRegionName);

      // Check actDampingId was already registerd
      if( idDampType.count( actDampingId ) == 0 ) {
        EXCEPTION( "Damping with id '" << actDampingId
                   << "' was not defined in 'dampingList'" );
      }

      dampingList_[actRegionId] = idDampType[actDampingId];
      if ( dampingList_[actRegionId] == RAYLEIGH ){
        RaylDampingData actRayl = *(idRaylData[actDampingId]);
        Double dampFreq;

        if( actRayl.freq == 0.0 ) {
          materials_[actRegionId]->GetScalar(dampFreq,RAYLEIGH_FREQUENCY,Global::REAL);
        } else { 
          dampFreq = actRayl.freq;
        }
        // Compute Rayleigh damping parameters
        materials_[actRegionId]->
        ComputeRayleighDamping( actRayl.alpha, actRayl.beta,
                                dampFreq, actRayl.ratioDeltaF, 
                                actRayl.adjustDamping, isComplex_ );
        regionRaylDamping_[actRegionId] = actRayl;
      }

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

    // ***************************************************************
    //   If no other damping type is specified and we have absorbing
    //   boundary conditions, then use ABCDAMP
    // ***************************************************************
    PtrParamNode bcNode = myParam_->Get( "bcsAndLoads", ParamNode::PASS );
    
    if( bcNode ) {
      ParamNodeList abcNodes = bcNode->GetList( "absorbingBCs" );

      for ( UInt i=0, n=abcNodes.GetSize(); i<n; ++i )
      {
        std::string regionName = abcNodes[i]->Get("name")->As<std::string>();
        absBCs_.Push_back( ptgrid_
          ->GetEntityList( EntityList::SURF_ELEM_LIST,
                           regionName, EntityList::REGION ) );
        Info->PrintF( pdename_,
                      "Apply Absorbing Boundary Conditions on surfaceRegion '%s'\n",
                      regionName.c_str() );
      }
      
      ParamNodeList impedNodes = bcNode->GetList( "impedance" );
      for ( UInt i=0, n=impedNodes.GetSize(); i<n; ++i )
      {
        shared_ptr<InhomDirichletBc> impedBC ( new InhomDirichletBc );
        impedBC->entities = ptgrid_
            ->GetEntityList( EntityList::SURF_ELEM_LIST,
                             impedNodes[i]->Get("name")->As<std::string>(),
                             EntityList::REGION );
        impedBC->value = impedNodes[i]->Get("value")->As<std::string>();
        impedBC->phase = impedNodes[i]->Get("phase", ParamNode::INSERT)
            ->As<std::string>();
        impedanceBCs_.Push_back( impedBC );
      }
    }
  }


  // *************************************************************
  //   Check what type of nonlinear PDE formulation should be used
  // *************************************************************
  void  AcousticPDE::InitNonLin() {

    SinglePDE::InitNonLin();

   //now do PDE specifics
    nonLin_ = false;

    std::map<std::string, NonLinType>::iterator it;
    for ( it=nonLinTypes_.begin() ; it != nonLinTypes_.end(); it++ ) {
      if ( (*it).second == WESTERVELT ) {
        nonLin_ = true;
        if ( formulation_ == ACOU_POTENTIAL )
          EXCEPTION( "Acoustic potential formulation not supported for "
                     << "Westervelt equation" );
      }
      else if( (*it).second == KUZNETSOV ) {
        nonLin_ = true;
        if ( formulation_ == ACOU_PRESSURE ) {
          EXCEPTION( "Acoustic pressure formulation not supported for"
                     << "Kuznetsov equation!" );
        }
      } else if( (*it).second == VARIABLE_SOS_CN1 ||
                 (*it).second == VARIABLE_SOS_CN2 ||
                 (*it).second == VARIABLE_SOS_CN2Mean ) {
        variableSpeedOfSoundCN_ = true;
      }

    }

  }


  void AcousticPDE::DefineIntegrators() {


    // =======================================================================
    //  Read flow data
    // =======================================================================
    ReadFlowData();

    Double density, compressibility, c0;
    Double coeffmass, coeffdamp;

    bool putSecondResult2EQNclass = false;

    for (UInt actSD = 0; actSD < subdoms_.GetSize(); actSD++) {

      // create new entity list
      shared_ptr<ElemList> actSDList( new ElemList(ptgrid_ ) );
      actSDList->SetRegion( subdoms_[actSD] );

      // get current region name and get grip of paramNode
      RegionIdType actRegion = subdoms_[actSD];
      std::string actRegionName;
      actRegionName = ptgrid_->GetRegion().ToString( actRegion );

      PtrParamNode actRegionNode =
        myParam_->Get("regionList")->GetByVal( "region", "name", actRegionName );
      std::string displInputId;
      actRegionNode->GetValue("displInputId", displInputId);
      if (updatedLagrangeForm_ &&  displInputId != "")
      {
        std::string solType;
        GridDisplData gridDispData;
        gridDispData.fileName4GridDisplacements_ = displInputId;
        actRegionNode->GetValue("displSolType", solType);
        gridDispData.solType = SolutionTypeEnum.Parse(solType);
        gridDisplData_[actRegion] = gridDispData;
      }

      // ********************************************************************
      //   Attention:
      //   In AcousticPDE ALL integrators are multiplied with density!
      // ********************************************************************

      materials_[actRegion]->GetScalar(density,DENSITY,Global::REAL);
      materials_[actRegion]->GetScalar(compressibility,ACOU_BULK_MODULUS,Global::REAL);
      c0 = sqrt(compressibility/density);

      // if pde couples with mechanic, we have to multiply the density by -1
      if ( isMechCoupled_ == true && formulation_ !=  ACOU_PRESSURE) {
        density *= -1.0;
      }

      // ********************************************************************
      //   Linear wave equation
      // ********************************************************************

      // Check for Perfectly matched layers
      if ( dampingList_[actRegion] == PML ) {
        
        // check if we are axi-symmetric: At the moment
        // PML damping is not working "perfectly" for 
        // axisymmetric setups
        if( isaxi_ ) {
          WARN("PML damping for axi-symmetric setup does not yet "
               << "work 'perfectly'! \nSpurious reflections can occur, "
               << "so you might consider using absorbing boundary "
               << "conditions instead!" );
              
        }
        
        //read data for PML layer
        //type of PML damping
        std::string dampingTypePML;
        
        // inner / outer region
        Matrix<Double> inner;
        Matrix<Double> outer;
        std::string coordSysId;
        //damping factor
        Double dampPML;
        
        bool isCPML = false;

        std::string id = actRegionNode->Get("dampingId")->As<std::string>();
        PtrParamNode pmlNode = myParam_->Get("dampingList")->GetByVal("pml", "id", id);
        ReadDataPML(dampingTypePML, inner, dampPML, coordSysId, pmlNode );
        dampPML *= c0;
        
        GetPMLLayerData(inner, outer, actRegion, coordSysId );

        // check for almost PML formulation
        std::string formsType;
        std::string helpStr;
        pmlNode->GetValue("aPML",helpStr);
        if ( helpStr == "yes" ) {
          //almost PML formulation: just makes since in 3D
          if ( dim_ != 3) 
            EXCEPTION( "APML formulation just makes sense in 3D");
          
          isAPML_ = true;
          std::cout << "\n DO APML\n" << std::endl;
        }

        //check for full PML without problematic L2-term
        pmlNode->GetValue("cPML",helpStr);
        if ( helpStr == "yes" ) {
          isCPML = true;
          std::cout << "\n DO cPML\n" << std::endl;
        }


        if ( analysistype_ == HARMONIC ) {     
          //====================================================================
          //	 stiffness integrator for PML
          //====================================================================
          if ( isAPML_ )
            formsType = "laplaceIntAPML";
          else
            formsType = "laplaceInt";
          
          //set real part
          BaseForm * bilinearStiffReal =
            new PMLInt(formsType, density, dampingTypePML, dampPML, isaxi_);
          
          bilinearStiffReal->SetPosPML(inner,outer, coordSysId);
          
          BiLinFormContext * stiffContextReal =
            new BiLinFormContext( bilinearStiffReal, STIFFNESS );
          
          stiffContextReal->SetPtPdes(this, this);
          stiffContextReal->SetResults( results_[0], results_[0],
                                        actSDList, actSDList );
          // stiffContextReal->SetEntryType(matType);
          assemble_->AddBiLinearForm( stiffContextReal);
          
          //====================================================================
          //	 mass integrator for PML
          //====================================================================
          
          if ( isAPML_ )
            formsType = "massIntAPML";
          else
            formsType = "massInt";

          Double massFactor = density/(c0*c0);
          
          //set real part
          BaseForm * bilinearMassReal =
            new PMLInt( formsType, massFactor, dampingTypePML, dampPML, isaxi_ );
          
          bilinearMassReal->SetPosPML(inner,outer, coordSysId);
          
          BiLinFormContext * massContextReal =
            new BiLinFormContext( bilinearMassReal, MASS);
          
          massContextReal->SetPtPdes(this, this);
          massContextReal->SetResults( results_[0], results_[0],
                                       actSDList, actSDList );
          // massContextReal->SetEntryType( matType );
          assemble_->AddBiLinearForm( massContextReal );
          
        } // end of pml part

        else if ( analysistype_ == TRANSIENT ) {
          putSecondResult2EQNclass = true;

          //====================================================================
          //	 pressure stiffness integrator for time domain PML
          //====================================================================
          Double factorPDE = density/(c0*c0);
          std::string formsType = "pressureStiff";
          
          BaseForm * bilinearPressStiff =
            new PMLTimeInt(formsType, factorPDE, dampingTypePML, dampPML, isaxi_);
          
          bilinearPressStiff->SetPosPML(inner,outer, coordSysId);
          
          BiLinFormContext * pressStiffContext =
            new BiLinFormContext( bilinearPressStiff, STIFFNESS );
          
          pressStiffContext->SetPtPdes(this, this);
          pressStiffContext->SetResults( results_[0], results_[0],
                                         actSDList, actSDList );
          assemble_->AddBiLinearForm( pressStiffContext);
          
          //====================================================================
          //	 pressure damping integrator for time domain PML
          //====================================================================
          formsType = "pressureDamp";
          BaseForm * bilinearPressDamp =
            new PMLTimeInt(formsType, factorPDE, dampingTypePML, dampPML, isaxi_);
          
          bilinearPressDamp->SetPosPML(inner,outer, coordSysId);
          
          BiLinFormContext * pressDampContext =
            new BiLinFormContext( bilinearPressDamp, DAMPING );
          
          pressDampContext->SetPtPdes(this, this);
          pressDampContext->SetResults( results_[0], results_[0],
                                        actSDList, actSDList );
          assemble_->AddBiLinearForm( pressDampContext);

          //====================================================================
          //	 pressure gradient integrator for time domain PML
          //====================================================================
          formsType = "pressureGrad";
          factorPDE = density;
          BaseForm * bilinearPressGrad =
            new PMLTimeInt(formsType, factorPDE, dampingTypePML, dampPML, isaxi_);
          
          bilinearPressGrad->SetPosPML(inner,outer, coordSysId);
          
          BiLinFormContext * pressGradContext =
            new BiLinFormContext( bilinearPressGrad, STIFFNESS );
          
          pressGradContext->SetPtPdes(this, this);
          pressGradContext->SetResults( results_[1], results_[0],
                                        actSDList, actSDList );
          assemble_->AddBiLinearForm( pressGradContext);

          //====================================================================
          //	 axuillary mass integrator for time domain PML
          //====================================================================
          factorPDE = density;
          MassInt * bilinearAuxMass  = new MassInt(factorPDE, dim_, isaxi_);

          BiLinFormContext * auxMassContext =
            new BiLinFormContext( bilinearAuxMass, DAMPING );
          auxMassContext->SetPtPdes(this, this);
          auxMassContext->SetResults( results_[1], results_[1],
                                   actSDList, actSDList );
          assemble_->AddBiLinearForm( auxMassContext );

          //====================================================================
          //	 axuillary div integrator for time domain PML
          //====================================================================
          factorPDE = density;
          formsType = "vecAuxillaryDiv";
          BaseForm * bilinearAuxDiv =
            new PMLTimeInt(formsType, factorPDE, dampingTypePML, dampPML, isaxi_);
          
          bilinearAuxDiv->SetPosPML(inner,outer, coordSysId);
          
          BiLinFormContext * auxDivContext =
            new BiLinFormContext( bilinearAuxDiv, STIFFNESS );
          
          auxDivContext->SetPtPdes(this, this);
          auxDivContext->SetResults( results_[0], results_[1],
                                        actSDList, actSDList );
          assemble_->AddBiLinearForm( auxDivContext);

          //====================================================================
          //	 axuillary stiffness integrator for time domain PML
          //====================================================================
          formsType = "vecAuxillaryStiff";
          factorPDE = density;
          BaseForm * bilinearAuxStiff =
            new PMLTimeInt(formsType, factorPDE, dampingTypePML, dampPML, isaxi_);
          
          bilinearAuxStiff->SetPosPML(inner,outer, coordSysId);
          
          BiLinFormContext * auxStiffContext =
            new BiLinFormContext( bilinearAuxStiff, STIFFNESS );
          
          auxStiffContext->SetPtPdes(this, this);
          auxStiffContext->SetResults( results_[1], results_[1],
                                       actSDList, actSDList );
          assemble_->AddBiLinearForm( auxStiffContext);

          if ( dim_ == 3 && isAPML_ == false ) {
            //3D computation: so we need in addition a scalar auxillary variable
                
            //====================================================================
            //	 scalar axuillary grad integrator for time domain PML
            //====================================================================
            if ( !isCPML ) {
              formsType = "auxGrad";
              factorPDE = -density;
              BaseForm * bilinearAuxGrad =
                new PMLTimeInt(formsType, factorPDE, dampingTypePML, dampPML, isaxi_);
              
              bilinearAuxGrad->SetPosPML(inner,outer, coordSysId);
              
              BiLinFormContext * auxGradContext =
                new BiLinFormContext( bilinearAuxGrad, STIFFNESS );
              
              auxGradContext->SetPtPdes(this, this);
              auxGradContext->SetResults( results_[1], results_[2],
                                          actSDList, actSDList );
              assemble_->AddBiLinearForm( auxGradContext);
            }

            //====================================================================
            //   	 scalar axuillary stiffness integrator for time domain PML
            //====================================================================
            factorPDE = density/(c0*c0);
            std::string formsType = "scalarAuxStiff";
            
            BaseForm * bilinearScalarAuxStiff =
              new PMLTimeInt(formsType, factorPDE, dampingTypePML, dampPML, isaxi_);
            
            bilinearScalarAuxStiff->SetPosPML(inner,outer, coordSysId);
            
            BiLinFormContext * scalarAuxStiffContext =
              new BiLinFormContext( bilinearScalarAuxStiff, STIFFNESS );
            
            scalarAuxStiffContext->SetPtPdes(this, this);
            scalarAuxStiffContext->SetResults( results_[0], results_[2],
                                               actSDList, actSDList );
            assemble_->AddBiLinearForm( scalarAuxStiffContext);
            
            //====================================================================
            //	 scalar axuillary mass integrator for time domain PML
            //====================================================================
            coeffmass = density;
            MassInt* bilinearMass  = new MassInt(coeffmass, 1, isaxi_ );
            
            BiLinFormContext * massScalarAuxContext =
              new BiLinFormContext( bilinearMass, DAMPING );
            massScalarAuxContext->SetResults( results_[2], results_[2],
                                              actSDList, actSDList );
            massScalarAuxContext->SetPtPdes(this, this);
            assemble_->AddBiLinearForm( massScalarAuxContext );        
            
            //====================================================================
            //	 pressure mass integrator for time domain PML
            //====================================================================
            coeffmass = -density;
            MassInt* bilinearMassPress  = new MassInt(coeffmass, 1, isaxi_ );
            BiLinFormContext * massPressContext =
              new BiLinFormContext( bilinearMassPress, STIFFNESS );
            massPressContext->SetResults( results_[2], results_[0],
                                          actSDList, actSDList );
            massPressContext->SetPtPdes(this, this);
            assemble_->AddBiLinearForm( massPressContext );   
          } 
          
         

          //====================================================================
          //	 standard mass and  stiffness integrator 
          //====================================================================

          BaseForm * bilinearStiff = new LaplaceInt( density, isaxi_ );        
          BiLinFormContext * stiffContext = 
            new BiLinFormContext( bilinearStiff, STIFFNESS );
          
          stiffContext->SetResults( results_[0], results_[0],
                                    actSDList, actSDList );
          stiffContext->SetPtPdes( this, this );
          assemble_->AddBiLinearForm( stiffContext );
  
          coeffmass = density / (c0*c0);
          MassInt* bilinearMass  = new MassInt(coeffmass, 1, isaxi_ );
          
          BiLinFormContext * massContext =
            new BiLinFormContext( bilinearMass, MASS );
          massContext->SetResults( results_[0], results_[0],
                                   actSDList, actSDList );
          massContext->SetPtPdes(this, this);
          assemble_->AddBiLinearForm( massContext );        
        }    
      } // end of PML part
      else {
        // stiffness integrator
        // ready for bimaterial topology optimization
        // if mech coupling case when not pressure formulation we have to do it with -1
        BaseForm::MaterialDescriptor::Type kdt = density > 0 ? BaseForm::MaterialDescriptor::SCALAR :
                                                               BaseForm::MaterialDescriptor::MINUS_SCALAR;

        BaseForm::MaterialDescriptor md1(kdt, FLUID, DENSITY,Global::REAL);
        BaseForm* bilinearStiff = new LaplaceInt(materials_[actRegion], md1, isaxi_, updatedLagrangeForm_);
        BiLinFormContext * stiffContext = 
          new BiLinFormContext( bilinearStiff, STIFFNESS );

        stiffContext->SetResults( results_[0], results_[0],
                                  actSDList, actSDList );
        stiffContext->SetPtPdes( this, this );

        // mass integrator
        // ready for bimaterial topology optimization

        // if mech coupling case when not pressure formulation we have to do it with -1
        BaseForm::MaterialDescriptor::Type mdt = density > 0 ? BaseForm::MaterialDescriptor::MAT_1_MAT_1_BY_MAT_2 :
                                                               BaseForm::MaterialDescriptor::MINUS_MAT_1_MAT_1_BY_MAT_2;

        BaseForm::MaterialDescriptor md2(mdt, FLUID, DENSITY,ACOU_BULK_MODULUS, Global::REAL);
        MassInt* bilinearMass = new MassInt(materials_[actRegion], md2, 1, isaxi_, updatedLagrangeForm_);

        if ( diagMass_ ) {
          // diagonal mass matrix
          bilinearMass->SetDiagMass();
        }

        BiLinFormContext * massContext =
          new BiLinFormContext( bilinearMass, MASS );
        massContext->SetResults( results_[0], results_[0],
                                 actSDList, actSDList );
        massContext->SetPtPdes(this, this);

        // *******************************************************************
        // Additional terms for Pierce Equation
        // *******************************************************************
        if ( regionFlowNodes_.count( actRegion) > 0 ) {
          if ( formulation_ == ACOU_PRESSURE )
            EXCEPTION("Pierce-Equation just possible in velocity potential formulation" );

          //read flow data
          PtrParamNode flowNode = regionFlowNodes_[actRegion];
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
        //   Damping Layer
        // ********************************************************************

        //check  for damping layer
        if ( dampingList_[actRegion] == DAMPLAYER ) {

          //type of damping fnc
          std::string dampFnc;

          std::string id = actRegionNode->Get("dampingId")->As<std::string>();
          PtrParamNode dampLayerNode = myParam_->Get("dampingList")
              ->GetByVal("dampLayer", "id", id);

          //damping data
          Double dampFactor, dampFactorMax, startRadius, stopRadius;
          Vector<Double> mPoint;
          ReadDataDampLayer( dampFnc, mPoint, dampFactor, dampFactorMax,
                             startRadius, stopRadius, dampLayerNode );

          //get the Rayleigh material parameters
          WARN( "The 'DampLayer' damping is not working at the moment."
                   "The usage of the Rayleigh parameters ALPHA and BETA need to"
                   "be changed" );
          Double alpha, beta, measFreq;
          
          std::string fac;
          materials_[actRegion]->GetScalar(alpha,RAYLEIGH_ALPHA,Global::REAL);
          materials_[actRegion]->GetScalar(beta,RAYLEIGH_BETA,Global::REAL);
          materials_[actRegion]->GetScalar(measFreq,RAYLEIGH_FREQUENCY,Global::REAL);

          // stiffness part
          if( isComplex_ ) {
            fac = lexical_cast<std::string>( beta * measFreq) + "/ f";
          } else {
            fac = lexical_cast<std::string>( beta );
          }

          stiffContext->SetSecDestMat( DAMPING, fac );
          stiffContext->SetDampLayer(dampFnc, mPoint, dampFactor,
                                     dampFactorMax, startRadius,
                                     stopRadius);
          // mass part
          if( isComplex_ ) {
            fac = lexical_cast<std::string>( alpha / measFreq) + "* f";
          } else {
            fac = lexical_cast<std::string>( alpha );
          }
          massContext->SetSecDestMat( DAMPING, fac );
          massContext->SetDampLayer(dampFnc, mPoint, dampFactor,
                                    dampFactorMax, startRadius,
                                    stopRadius);
        }

        // ********************************************************************
        //   Additional terms for damping
        // ********************************************************************
        if ( dampingList_.size() > 0 ) {

          // We check, if damping has been specified for all regions.
//          if ( dampingList_.size() != subdoms_.GetSize() ) {
//            WARN("Mismatch between dampingList_ and subdoms_!"
//                 << "Size(dampingList_): " << dampingList_.size()
//                 << "Size(subdoms_): " << subdoms_.GetSize());
//}

          if (dampingList_[actRegion] == RAYLEIGH) {
            // This works even after assemble_->AddIntegrator() is executed
            //   because of the pointers...
            
            RaylDampingData & actDamp = (regionRaylDamping_[actRegion]);
            stiffContext->SetSecDestMat( DAMPING, actDamp.beta );
            massContext->SetSecDestMat( DAMPING, actDamp.alpha );
          }


          else if ( dampingList_[actRegion] == THERMOVISCOUS ) {
            Double viscousAlpha;
            materials_[actRegion]->GetScalar(viscousAlpha, ACOU_ALPHA,Global::REAL);

            coeffdamp  =  density * 2.0 * viscousAlpha * c0;
            BaseForm * bilinearStiff  = new LaplaceInt(coeffdamp, isaxi_);
            BiLinFormContext * dampContext =
              new BiLinFormContext(bilinearStiff, DAMPING );
            dampContext->SetResults( results_[0], results_[0],
                                     actSDList, actSDList );
            dampContext->SetPtPdes(this, this);
            assemble_->AddBiLinearForm( dampContext );
          }

          else if ( (analysistype_ != HARMONIC) &&
                    (dampingList_[actRegion] == FRACTIONAL_GL ||
                     dampingList_[actRegion] == FRACTIONAL_GL_INT ) ) {

            Double fracAlpha, fracExp;
            materials_[actRegion]->GetScalar(fracAlpha,ACOU_ALPHA,Global::REAL);
            materials_[actRegion]->GetScalar(fracExp,FRACTIONAL_EXPONENT,Global::REAL);

            coeffdamp = - density*2.0*fracAlpha/c0/sin((fracExp-1.0)*PI/2.0);

            MassInt * bilinearDamp = new MassInt(coeffdamp, 1, isaxi_);

            Double fracDampCoeff = GetFracDampMatrixCoeff( actRegion );
            bilinearDamp->SetSecondFactor( lexical_cast<std::string>(fracDampCoeff ) );

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

          else if  ( (analysistype_ != HARMONIC) &&
                     (dampingList_[actRegion] == FRACTIONAL_BLANK ||
                      dampingList_[actRegion] == FRACTIONAL_BLANK_INT) ) {

            Double fracAlpha, fracExp;
            materials_[actRegion]->GetScalar(fracAlpha,ACOU_ALPHA,Global::REAL);
            materials_[actRegion]->GetScalar(fracExp,FRACTIONAL_EXPONENT,Global::REAL);

            coeffdamp =  - density*2.0*fracAlpha/c0/sin((fracExp-1.0)*PI/2.0);
            // prefactor of blank alg
            coeffdamp *= exp(-gammaln(1.0- (fracExp- 1.0)) );
            // weight factor of index 0
            coeffdamp *= 1.0/(1.0- (fracExp- 1.0));

            MassInt * bilinearDamp = new MassInt(coeffdamp, 1, isaxi_);

            Double fracDampCoeff = GetFracDampMatrixCoeff( actRegion );
            bilinearDamp->SetSecondFactor( lexical_cast<std::string>(fracDampCoeff) );

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

          else if ( (analysistype_ == HARMONIC) &&
                    (dampingList_[subdoms_[actSD]] == FRACTIONAL_GL ||
                     dampingList_[subdoms_[actSD]] == FRACTIONAL_GL_INT ||
                     dampingList_[subdoms_[actSD]] == FRACTIONAL_BLANK ||
                     dampingList_[subdoms_[actSD]] == FRACTIONAL_BLANK_INT) ) {

            Double fracAlpha, fracExp;
            materials_[subdoms_[actSD]]->GetScalar(fracAlpha,ACOU_ALPHA,Global::REAL);
            materials_[subdoms_[actSD]]->GetScalar(fracExp,FRACTIONAL_EXPONENT,Global::REAL);

            Double factorReal, factorImag;
            factorReal = density*2.0*fracAlpha/c0/tan((fracExp-1.0)*PI/2.0);
            factorImag = density*2.0*fracAlpha/c0;

            // set up real and imaginary part of damping matrix
            BaseForm * bilinearDampReal = new MassInt(factorReal, 1, isaxi_);
            BaseForm * bilinearDampImag = new MassInt(factorImag, 1, isaxi_);

            std::string omegaFac;
            omegaFac = "exp((" +  lexical_cast<std::string>(fracExp) + "+1)*ln(2*pi*f))";
            bilinearDampReal->SetSecondFactor( omegaFac );
            bilinearDampImag->SetSecondFactor( omegaFac );

            // Choose stiffness matrix, because in HARMONIC calculation mass and
            //  damping matrix are multiplied by multiples of omega
            //  See method Matrix2Harmonic in assemble.cc
            BiLinFormContext * dampContextReal =
              new BiLinFormContext( bilinearDampReal, STIFFNESS );
            dampContextReal->SetPtPdes(this, this);
            dampContextReal->SetResults( results_[0], results_[0],
                                         actSDList, actSDList );

            BiLinFormContext * dampContextImag =
              new BiLinFormContext( bilinearDampImag, STIFFNESS );
            dampContextImag->SetPtPdes(this, this);
            dampContextImag->SetResults( results_[0], results_[0],
                                         actSDList, actSDList );
            // set imaginary flag of matrix context
            Global::ComplexPart complexType = Global::IMAG;
            dampContextImag->SetEntryType(complexType);

            assemble_->AddBiLinearForm( dampContextReal );
            assemble_->AddBiLinearForm( dampContextImag );
          }
        }

        // Finally add the stiffness/mass integrators
        assemble_->AddBiLinearForm( stiffContext );
        assemble_->AddBiLinearForm( massContext );
      }

      // Give result to equation numbering class
      eqnMap_->AddResult( *results_[0], actSDList );

      // check if second (third) result type is active (currently just for
      // time domain PML )
      if ( putSecondResult2EQNclass ) {
        eqnMap_->AddResult( *results_[1], actSDList );
        putSecondResult2EQNclass = false;
        if ( dim_ == 3 && isAPML_ == false ) 
          eqnMap_->AddResult( *results_[2], actSDList );
      }
    }

    // **********************************************************************
    //   inhom. Neumann boundary condition
    // **********************************************************************
    for( UInt iBc = 0; iBc < inBcs_.GetSize(); iBc++ ) {

      // get current Bc
      InhomNeumannBc const & actBc = *inBcs_[iBc];

      std::string magnStr;
      if ( isMechCoupled_ == true && formulation_ !=  ACOU_PRESSURE ) {
        // reverse the sign to get a symmetric system
        magnStr = "-(" + actBc.value + ")";
      } else {
        magnStr = actBc.value;
      }

      //BaseForm *neumannBC = new VolumeSrcInt( amplitude, isaxi_ );
      LinearSurfForm *neumannBC =
          new LinNeumannInt( magnStr, actBc.phase,
                             formulation_ == ACOU_PRESSURE ?
                                 ACOU_ACCELERATION : ACOU_VELOCITY,
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
    if ( absBCs_.GetSize() > 0 )
    {
      std::string factor;
      // In the case of acou-mech coupling we have to multiply the
      // abc-Integrator matrix with -1
      if ( isMechCoupled_ == true && formulation_ !=  ACOU_PRESSURE ) {
        factor = "-1.0";
      }
      
      for (UInt actSD=0, numSD=absBCs_.GetSize(); actSD<numSD; ++actSD)
      {

        AbsorbingBCsInt * bilinear_damp =
            new AbsorbingBCsInt( true, factor, "", isaxi_);
        bilinear_damp->SetFirstVoluInfo(pdename_, materials_ );

        BiLinFormContext * abcContext =
          new BiLinFormContext( bilinear_damp, DAMPING );
        abcContext->SetPtPdes(this, this);
        abcContext->SetResults( results_[0], results_[0],
                                absBCs_[actSD], absBCs_[actSD] );
        assemble_->AddBiLinearForm( abcContext );

        // Give result to equation numbering class
        eqnMap_->AddResult( *results_[0], absBCs_[actSD] );
      }
    }

    // =======================================================================
    //  surface-integration: Impedance boundaries
    // =======================================================================
    for ( UInt i=0, n=impedanceBCs_.GetSize(); i<n; ++i )
    {
      std::string phase;
      if ( analysistype_ == HARMONIC ) {
        phase = impedanceBCs_[i]->phase;
      }
      else {
        WARN( "You are using an acoustic impedance boundary condition in a\n"
             << " transient analysis. In this case the impedance cannot depend\n"
             << " on frequency and it cannot be complex (i.e. no phase change).\n"
             << " Is this really what you want?");
      }
      
      AbsorbingBCsInt *impedInt =
          new AbsorbingBCsInt(false, impedanceBCs_[i]->value, phase, isaxi_);
      impedInt->SetFirstVoluInfo( pdename_, materials_ );

      BiLinFormContext *impedContext =
          new BiLinFormContext( impedInt, DAMPING );
      impedContext->SetPtPdes( this, this );
      impedContext->SetResults( results_[0], results_[0],
                                impedanceBCs_[i]->entities,
                                impedanceBCs_[i]->entities );
      
      assemble_->AddBiLinearForm( impedContext );
      
      eqnMap_->AddResult( *results_[0], impedanceBCs_[i]->entities );
    }
    

    // =======================================================================
    //check for flow sources: pressure formulation
    // =======================================================================
    PtrParamNode bcNode = myParam_->Get("bcsAndLoads", ParamNode::PASS);
    
    if ( bcNode ) {
      ParamNodeList flowSrcNodes = bcNode->GetList("flowPresSrc");
  
      std::string rhsFileId;
      try {
        // iterate over all parameter nodes
        for( UInt i = 0; i < flowSrcNodes.GetSize(); i++ ) {
          std::string rhsRegion(flowSrcNodes[i]->Get("region")->As<std::string>());
          rhsFileId = flowSrcNodes[i]->Get("inputId")->As<std::string>();
  
          linFlowPressureRHSInt* sourceRHSInt = 
            new linFlowPressureRHSInt( isaxi_, rhsFileId, rhsRegion );
  
          LinearFormContext* sourceRHSContext = new LinearFormContext( sourceRHSInt );
          sourceRHSContext->SetPtPde( this );
  
          shared_ptr<ElemList> rhsElemList( new ElemList(ptgrid_ ) );
          rhsElemList->SetRegion( ptgrid_->GetRegion().Parse(rhsRegion) );
          sourceRHSContext->SetResult( results_[0], rhsElemList );
          assemble_->AddLinearForm( sourceRHSContext );
        }
      }
      catch(Exception & ex) {
        RETHROW_EXCEPTION(ex, "Could not assemble RHSFlowPressureRHSInt integrator"
                          <<" in acousticPDE" );
      }
    }


    
    // =======================================================================
    // Integrators for NonConforming Interfaces
    // =======================================================================
    PtrParamNode ncIfaceListNode
        = param->Get("domain")->Get("ncInterfaceList", ParamNode::PASS);

    // Get index of LAGRANGE_MULT result, just in case... who knows...
    UInt lmResultIdx = 0;
    for(UInt i=0, n=results_.GetSize(); i<n; i++) {
      if(results_[i]->resultType == LAGRANGE_MULT) {
        lmResultIdx = i;
        break;
      }
    }
    LOG_DBG2(acoupde) << "NonMatching: Index of LAGRANGE_MULT result: "
                     << lmResultIdx;

    for( UInt i = 0; i < ncIFaces_.GetSize(); i++ ) {

      // get regionId of Lagrangian surface
      StdVector<std::string> keyVec, attrVec, valVec;
      std::string slaveSide;
      std::string ncIfaceName = ptgrid_->GetRegion().ToString(ncIFaces_[i]);

      if (!ncIfaceListNode) {
        EXCEPTION("No ncInterfaces defined in domain section.");
      }
      PtrParamNode curNciNode = ncIfaceListNode->GetByVal("ncInterface", "name",
                                                   ncIfaceName);
      slaveSide = curNciNode->Get("slaveSide")->As<std::string>();

      // Part 1: Define integrator M(Psi, Lambda) on
      //         non-conforming interface
      LOG_DBG2(acoupde) << "NonMatching: Defining nonconforming integrator"
                        << " for M on interface '"
                        << ptgrid_->GetRegion().ToString(ncIFaces_[i]) << "'.";
      shared_ptr<ElemList> actNcList( new ElemList(ptgrid_ ) );
      actNcList->SetRegion( ncIFaces_[i] );

      NonConformingInt * ncInt =
        new NonConformingInt( 1, isaxi_ );

      NcBiLinFormContext * stiffIntDescr =
        new NcBiLinFormContext( ncInt , STIFFNESS );

      // Force assembling of M(Psi, Lambda)^T
      stiffIntDescr->SetCounterPart( true );

      stiffIntDescr->SetPtPdes(this, this);
      stiffIntDescr->SetResults( results_[0], results_[lmResultIdx],
                                 actNcList, actNcList );

      assemble_->AddBiLinearForm( stiffIntDescr );


      // Part 2: Define integrator D(Psi, Lambda) on
      //         Lagrangian surface
      LOG_DBG2(acoupde) << "NonMatching: Defining mass integrator"
                        << " for D on interface '"
                        << ptgrid_->GetRegion().ToString(ncIFaces_[i]) << "'.";
      shared_ptr<SurfElemList> actSDList( new SurfElemList(ptgrid_ ) );
      actSDList->SetRegion( ptgrid_->GetRegion().Parse(slaveSide));

      // D(Psi, Lambda) has the form of a standard mass
      // integrator with factor 1.0
      MassInt * dMatInt = new MassInt( 1.0, 1, isaxi_ );
      BiLinFormContext * dMatContext =
        new BiLinFormContext( dMatInt, STIFFNESS );

      // Force assembling of D(Psi, Lambda)^T
      dMatContext->SetCounterPart( true );
      dMatContext->SetPtPdes( this, this );
      dMatContext->SetResults( results_[0], results_[lmResultIdx],
                               actSDList, actSDList );

      assemble_->AddBiLinearForm( dMatContext );

      // Give result LAGRANGE_MULT to equation numbering class
      eqnMap_->AddResult( *results_[lmResultIdx], actSDList );
    }

    // =======================================================================
    // Integrators for acoustic RHS values (i.e. Lighthill sources)
    // =======================================================================


    std::string rhsRegion;
    PtrParamNode rhsValuesNode, bcsNode;
    StdVector<PtrParamNode> regionList;
    std::set<std::string> srcRegions;

    bcsNode = myParam_->Get("bcsAndLoads", ParamNode::PASS );
    if( bcsNode )
      rhsValuesNode = bcsNode->Get("rhsValues", ParamNode::PASS );
    if(!rhsValuesNode)
      return;

    regionList = rhsValuesNode->GetList("region");
    for (StdVector<PtrParamNode>::iterator regionIter = regionList.Begin();
        regionIter != regionList.End(); ++regionIter)
    {
      shared_ptr<NodeList> acouRHSRegionNodeList( new NodeList(ptgrid_ ) );
      try
      {
        rhsRegion = (*regionIter)->Get("name")->As<std::string>();
        PtrParamNode intNode;
        if((*regionIter)->Has("interpolation"))
        {
          intNode = (*regionIter)->Get("interpolation", ParamNode::PASS);
          intNode->GetValue("justInterpolate", justInterpolate_, ParamNode::PASS);
          
          // make sure that no source region is given twice, because that
          // does not work with our interpolation algorithm
          std::set<std::string> mySrcRegions, intersection;
          boost::char_separator<char> sep(";| ");
          boost::tokenizer< char_separator<char> > tok(intNode
              ->Get("srcRegions")->Get("names")->As<std::string>(), sep);
          mySrcRegions.insert(tok.begin(), tok.end());
          std::set_intersection(srcRegions.begin(), srcRegions.end(),
                                mySrcRegions.begin(), mySrcRegions.end(),
                                std::inserter(intersection, intersection.end()));
          if ( intersection.empty() ) {
            std::copy(mySrcRegions.begin(), mySrcRegions.end(),
                      std::inserter(srcRegions, srcRegions.end()));
          } else {
            EXCEPTION("Conservative interpolation: the same source region must"
                      << " not be used with two or more destination regions");   
          }
        }
      } catch (Exception& ex)
      {
        RETHROW_EXCEPTION(ex, "Error while trying to read parameters for AcouRHSLinForm.");
      }

      if(rhsRegion != "")
      {
        acouRHSRegionNodeList->SetNodesOfRegion( ptgrid_->GetRegion().Parse(rhsRegion) );

        AcouRHSLinForm* acouRHSInt = new AcouRHSLinForm((*regionIter));

        LinearFormContext * acouRHSContext = new LinearFormContext( acouRHSInt );
        acouRHSContext->SetPtPde( this );
        acouRHSContext->SetResult( results_[0], acouRHSRegionNodeList );
        assemble_->AddLinearForm( acouRHSContext );
      }
    }
  }

  void AcousticPDE::DefineSolveStep() {
    solveStep_ = new SolveStepAcoustic(*this, justInterpolate_);
  }


  // ======================================================
  // TIME STEPPING SECTION
  // ======================================================

  void AcousticPDE::InitTimeStepping() {

	  // timestepping formulation
	  PtrParamNode myLinSysNode = FindLinearSystem( pdename_ );

	  // <system name="acoustic"/> exists
	  if( myLinSysNode ) {

	    std::string str = "";
	    myLinSysNode->GetValue( "timeSteppingFormulation", str,  ParamNode::PASS );
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
	  }

	  PtrParamNode systemNode = FindLinearSystem(pdename_);
	  
    // this includes rayleigh and thermoviscous damping
    if ( fracDamping_ == false ) {
      if ( effectiveMass_ == true ) {
        if ( diagMass_ == true ) {
          //explicit time stepping
          TS_alg_ = new NewmarkEffMass( algsys_, systemNode, true );
        }
        else {
          TS_alg_ = new NewmarkEffMass( algsys_, systemNode, false );
        }
      }
      else if ( effectiveMass_ == false ) {
        TS_alg_ = new Newmark( algsys_, systemNode );
      }
    }
    else {
      if ( effectiveMass_ == false )
        TS_alg_ = new NewmarkFracDamp( algsys_,
                                       pdeId_, eqnMap_, ptgrid_, this,
                                       results_[0],
                                       subdoms_, dampingList_, systemNode );
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
        ptgrid_->GetRegion().Parse( couplRegions, regionIds );

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
  }


  void AcousticPDE::CalcOutputCoupling() {


    UInt dof;
    SolutionType quantity;
    StdVector<Elem*> * couplingElems = NULL;
    StdVector<UInt> * couplingNodes = NULL;
    SingleVector * temp_values = NULL;
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
        break;

      case ELEM:
        break;
      }
    }
  }


  void AcousticPDE::
  CalcMechCouplingRHS( StdVector<Elem*> * couplingElems,
                       StdVector<UInt> & couplingNodes,
                       Vector<Double>& elemCouplingSols,
                       UInt couplingdof ) {


    EXCEPTION( "Not working at the moment" );

    //     Double density = 0.0;
    //     Double sign = 0.0;
    //     Integer matIndex = -1;
    //     Elem * ptVolElem = NULL;
    //     Matrix<Double> ptCoord, elemMat;
    //     Vector<Double> sol, forceOnElem, normal;

    //     elemCouplingSols.Init();

    //     for (UInt actElem=0; actElem<couplingElems->GetSize(); actElem++) {

    //       // Perform cast from volume element to surface element, since
    //       // mech-acou coupling makes only sense on surface elements
    //       SurfElem * actCoupleElem =
    //         dynamic_cast<SurfElem*> ((*couplingElems)[actElem]);

    //       if (actCoupleElem == NULL) {
    //         EXCEPTION( "No elements found for coupling!" );
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
    //         EXCEPTION( "AcousticPDE::CalcMechCouplingRHS: The two volume "
    //                  << "element neighbours of surface element Nr. "
    //                  << actCoupleElem->elemNum << " do not belong to my regions!");
    //       }

    //       // Assign correct density
    //       materials_[subdoms_[matIndex]]->GetScalar(density,DENSITY,Global::REAL);

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

    //       for (UInt actNode=0; actNode<ptCoord.GetNumCols(); actNode++) {
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


  template <class TYPE>
  void AcousticPDE::
  CalcHeatCouplingRHS(Vector<Double> & sourceValue,
                      StdVector<StdVector<UInt> > & elemNodeToCouplingNode,
                      UInt actCoupling,
                      UInt numCouplingNodes) {

    // get the coupling regions
    StdVector<std::string> couplRegions;
    ptCoupling_->GetOutputRegions(actCoupling, couplRegions);
    StdVector<RegionIdType> regionIds;
    ptgrid_->GetRegion().Parse( couplRegions, regionIds );

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
      materials_[regionIds[reg]]->GetScalar(density,DENSITY,Global::REAL);

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


  bool AcousticPDE::HasOutput(SolutionType output) {
    if ((output == ACOU_FORCE) || (output == ACOU_POWERDENSITY)) {
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
      ExtractDerivResult( result, FIRST_DERIV );
      break;

    case ACOU_POTENTIAL_DERIV_2:
      ExtractDerivResult( result, SECOND_DERIV );
      break;

    case ACOU_PRESSURE_DERIV_1:
      ExtractDerivResult( result, FIRST_DERIV );
      break;

    case ACOU_PRESSURE_DERIV_2:
      ExtractDerivResult( result, SECOND_DERIV );
      break;

    case ACOU_RHS_LOAD:
      if( isComplex_ ) {
        ExtractRhsResult<Complex>( result, results_[0] );
      } else {
        ExtractRhsResult<Double>( result, results_[0] );
      }
      break;

    case ACOU_INTENSITY:
      if( isComplex_ ) {
        CalcAcouIntensity<Complex>( result );
      } else {
        EXCEPTION( "Acoustic intensity only computable for harmonic "
                   << "analysis!" );
      }
      break;

    case ACOU_SURFINTENSITY:
      if( isComplex_ ) {
        CalcAcouSurfIntensity<Complex>( result );
      } else {
        EXCEPTION( "Kaltenbacher's intensity only computable for harmonic "
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

    case GRAD_ACOU_SOLUTION:
      if(isComplex_)
        CalcGradSolution<Complex>(result);
      else
        CalcGradSolution<Double>(result);
      break;

    case ACOU_ENERGY:
      if( isComplex_ ) {
        //        CalcAcouEnergy<Complex>( result );
        EXCEPTION( "Acoustic energy currently only computable for time "
                   << "analysis!" );
      } else {
        //first compute particle velocity
        if ( formulation_ == ACOU_PRESSURE )
          CalcVelFromPressure( result );
        CalcAcouEnergy<Double>( result );
      }
      break;

    case ACOU_PSEUDO_DENSITY:
      if(domain->GetErsatzMaterial(false) == NULL) // no excpetion
        result->Init();
      else
        domain->GetErsatzMaterial()->ExtractResults(result, isComplex_);
      break;

    default:
      WARN( "Result type not computable by acoustic PDE" );
      break;
    }
  }


  template <class TYPE>
  void AcousticPDE::CalcForce( shared_ptr<BaseResult> vals ) {

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

      Vector<TYPE> valueElem;
      if ( formulation_ == ACOU_POTENTIAL ) {

        Integer matIndex = -1;
        // Try to find according region for first neighbouring volume
        // element of the surface element
        matIndex = subdoms_.Find(actSaveElem->ptVolElem1->regionId);

        // If first volume element does not belong to acoustic PDE, try the
        // second one
        // Elem * ptVolElem = NULL; // TODO: Unused variable ptVolElem
        if ( matIndex == -1 ) {
          matIndex = subdoms_.Find(actSaveElem->ptVolElem2->regionId);
          // ptVolElem = actSaveElem->ptVolElem2;
        }
        // else {
          // ptVolElem = actSaveElem->ptVolElem1;
        // }

        if ( matIndex == -1) {
          EXCEPTION( "AcousticPDE::CalcForce: The two volume element"
                     << "neighbours of surface element no. "
                     << actSaveElem->elemNum << " don't belong to my region!" );
        }

        // Assign correct density
        Double density;
        materials_[subdoms_[matIndex]]->GetScalar(density,DENSITY,Global::REAL);

        // retrieve 1st derivative, since F = rho * dpsi/dt * A
        //ElemList elemList(ptgrid_);
        //elemList.SetElement( actSaveElem );
        //EntityIterator it = elemList.GetIterator();
        GetDerivSolVecOfElement(valueElem,  it, results_[0]);
        valueElem *= density;
      }
      else if ( formulation_ == ACOU_PRESSURE ) {

        // retrieve solution, since F = p * A
        //ElemList elemList(ptgrid_);
        //elemList.SetElement( actSaveElem );
        //EntityIterator it = elemList.GetIterator();
        GetSolVecOfElement(valueElem,it,results_[0]);
      }

      const UInt nrIntPts= ptElem->GetNumIntPoints();
      const Vector<Double> & intWeights = ptElem->GetIntWeights();

      TYPE forceElem = 0.0;
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
      materials_[actRegion]->GetScalar(density,DENSITY,Global::REAL);

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
      GetDerivSolVecOfElement(valueElem, it, results_[0]);

      elemPressure = valueElem * shapeFnc * density;
      actVal[it.GetPos()] = elemPressure;
    }
  }


  template <class TYPE>
  void AcousticPDE::CalcAcouIntensity( shared_ptr<BaseResult> vals ) {

    //get frequency
    MathParser * parser = domain->GetMathParser();
    parser->SetExpr( mHandle_, "f" );
    Double actFreq = parser->Eval( mHandle_ );

    // factor 0.5 is due to the fact, that the values are peak values
    Complex multVal = 0;
    if ( formulation_ == ACOU_PRESSURE ) {
      multVal = Complex(0.0, -0.5/ (2.0*PI*actFreq) );
    } else {
      multVal = Complex(0.0, -0.5*(2.0*PI*actFreq));
    }

    // Create operator for gradient computation of solution
    NodeStoreSol<TYPE> * solhelp = dynamic_cast<NodeStoreSol<TYPE>*>(sol_);
    GradientFieldOp<TYPE> * gradOp =
      new GradientFieldOp<TYPE>(ptgrid_, this, eqnMap_, *solhelp, results_[0]->fctType, isaxi_);

    // get data from result object and resize its vector
    Result<TYPE> & actRes = dynamic_cast<Result<TYPE>&>(*vals);
    EntityIterator it = actRes.GetEntityList()->GetIterator();
    Vector<TYPE> & actVal = actRes.GetVector();
    actVal.Resize( actRes.GetEntityList()->GetSize() * dim_  );


    Double density = 0.0;
    TYPE elemSol = 0.0;
    Vector<TYPE> gradElemSol(dim_),  elemIntensity(dim_);

    // loop over all volume elements
    for ( it.Begin(); !it.IsEnd(); it++ ) {

      // density of element
      BaseFE * ptElem = it.GetElem()->ptElem;
      RegionIdType actRegion = it.GetElem()->regionId;
      materials_[actRegion]->GetScalar(density,DENSITY,Global::REAL);

      // get element coordinates
      StdVector<UInt> const & connect = it.GetElem()->connect;
      Matrix<Double> ptCoord;
      ptgrid_->GetElemNodesCoord( ptCoord, connect, false );

      // get shape function at center of the element
      Vector<Double> shapeFnc;
      Vector<Double> LCoord;
      ptElem->GetCoordMidPoint(LCoord);
      ptElem->GetShFnc(shapeFnc,LCoord,it.GetElem());

      // solution at center of element
      Vector<TYPE> valueElem;
      GetSolVecOfElement(valueElem, it, results_[0]);
      elemSol = valueElem * shapeFnc;

      // get the conjugate complex value
      elemSol = std::conj(elemSol);

      // calculate gradient at center of element
      gradOp->CalcElemGradField(gradElemSol, it, LCoord, 1.0);

      TYPE factorI;
      if ( formulation_ == ACOU_PRESSURE ) {
        factorI   = multVal * elemSol / density;
        elemIntensity = gradElemSol * factorI;
      } else {
        factorI   = multVal * elemSol * density;
        elemIntensity = gradElemSol * factorI;
      }

      // loop over dofs
      for(UInt iDim = 0; iDim < dim_; iDim++ ) {
        actVal[it.GetPos()*dim_ + iDim] = elemIntensity[iDim];
      }
    }

    delete gradOp;
  }

  template <class TYPE>
  void AcousticPDE::CalcGradSurfaceElement( shared_ptr<BaseResult> vals, SolutionType solType,
                                            bool apply_normal, StdVector<Vector<TYPE> >& grad_out )
  {
    //some help variables
    Vector<Double> lCoordSurf, lCoordVol, normal;
    Vector<TYPE> gradVal(dim_);
    BaseFE * ptSurfElemFE, * ptVolElemFE;

    // Create vector with interpolation coordinate.
    // For simplicity we only evaluate the integral
    // in coordinate origin
    lCoordSurf.Resize(dim_-1);
    lCoordSurf.Init(0);

    // Create operator for gradient computation of solution
    NodeStoreSol<TYPE> * solhelp = dynamic_cast<NodeStoreSol<TYPE>*>(sol_);
    GradientFieldOp<TYPE> * gradOp =
      new GradientFieldOp<TYPE>(ptgrid_, this, eqnMap_, *solhelp,
                                results_[0]->fctType, isaxi_);

    Result<TYPE> & actRes = dynamic_cast<Result<TYPE>&>(*vals);
    EntityIterator it = actRes.GetEntityList()->GetIterator();

    grad_out.Resize(actRes.GetEntityList()->GetSize());

    // loop over all surface elements
    for ( it.Begin(); !it.IsEnd(); it++ )
    {
      const SurfElem * actSurfElem = it.GetSurfElem();
      // Determine, which volume element is the right neighbor for the
      // calculation;
      // our normal should point out of the correct neighbor volume element!
      bool vol_1 = surfNeighborRegions_[vals] == actSurfElem->ptVolElem1->regionId;
      Elem* ptVolElem   = vol_1 ? actSurfElem->ptVolElem1 : actSurfElem->ptVolElem2;
      double normal_sign = vol_1 ? 1.0 : -1.0; // from CalcAcouPower() -> seem to correct surfElem->normalSign

      ptSurfElemFE = actSurfElem->ptElem;
      ptVolElemFE = ptVolElem->ptElem;

      const StdVector<UInt> & surfConnect = actSurfElem->connect;
      const StdVector<UInt> & volConnect = ptVolElem->connect;

      // calculate volume integration coordinates from
      // surface integration coordinate for evaluating the
      // gradient on the surface of the volume element
      ptSurfElemFE->GetCoordMidPoint(lCoordSurf);
      ptVolElemFE->GetLocalIntPoints4Surface(surfConnect, volConnect,
                                             lCoordSurf, lCoordVol);
      // Calc gradient
      ElemList elList(ptgrid_);
      elList.SetElement( ptVolElem );
      EntityIterator it2 = elList.GetIterator();
      gradOp->CalcElemGradField(gradVal, it2, lCoordVol,1.0);

      if(apply_normal) // check also surfElem->normalSign
        gradVal *= normal_sign;

      grad_out[it.GetPos()] = gradVal;

      LOG_DBG3(acoupde) << "CGSE: se=" << actSurfElem->elemNum 
                        << " lCoordVol=" << lCoordVol.ToString() 
                        << " grad=" << gradVal.ToString();

    }

    delete gradOp;

  }




  template <class TYPE>
  void AcousticPDE::CalcAcouSurfIntensity( shared_ptr<BaseResult> vals ) {

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

    //some help variables
    Vector<Double> lCoordSurf;
    Vector<TYPE> elemIntensity(dim_);
    // BaseFE * ptSurfElemFE; // TODO: Unused variable ptSurfElemFE
    Double density  = 0.0;

    // Create vector with interpolation coordinate.
    // For simplicity we only evaluate the integral
    // in coordinate origin
    lCoordSurf.Resize(dim_-1);
    lCoordSurf.Init(0);

    TYPE factorI;
    Result<TYPE> & actRes = dynamic_cast<Result<TYPE>&>(*vals);
    EntityIterator it = actRes.GetEntityList()->GetIterator();

    Vector<TYPE> & actVal = actRes.GetVector();
    actVal.Resize( actRes.GetEntityList()->GetSize() * dim_ );

    // calculate the gradient by CalcGradSurfaceElement()
    StdVector<Vector<TYPE> > grad_out;
    // in the original code the normal was not considered??
    CalcGradSurfaceElement(vals, solType, false, grad_out);

    NodeStoreSol<TYPE> * solhelp = dynamic_cast<NodeStoreSol<TYPE>*>(sol_);

    // loop over all surface elements
    for ( it.Begin(); !it.IsEnd(); it++ ) {

      const SurfElem * actSurfElem = it.GetSurfElem();

      // Determine, which volume element is the right neighbor for the
      // calculation;
      // our normal should point out of the correct neighbor volume element!
      bool vol_1 = surfNeighborRegions_[vals] == actSurfElem->ptVolElem1->regionId;
      Elem* ptVolElem = vol_1 ? actSurfElem->ptVolElem1 : actSurfElem->ptVolElem2;

      // ptSurfElemFE = actSurfElem->ptElem;

      const StdVector<UInt> & surfConnect = actSurfElem->connect;


      BaseMaterial * myMat = materials_[ptVolElem->regionId];
      myMat->GetScalar(density,DENSITY,Global::REAL);

      Matrix<Double> CornerCoords;
      ptgrid_->GetElemNodesCoord( CornerCoords, surfConnect, false );

      Vector<TYPE>& gradVal = grad_out[it.GetPos()];

      //get average solution
      TYPE elemSol = 0;
      TYPE nodeSol;
      for ( UInt k=0; k<surfConnect.GetSize(); k++) {
        solhelp->Get(solType, surfConnect[k]-1,0, nodeSol);
        elemSol += nodeSol;
      }
      elemSol /= (Double)surfConnect.GetSize();

      LOG_DBG3(acoupde) << "CASI: se=" << actSurfElem->elemNum 
                        << " grad=" << gradVal.ToString() 
                        << " elemSol=" << elemSol;

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
  }




  template<class TYPE>
  void AcousticPDE::CalcAcouPower(shared_ptr<BaseResult> vals)
  {
    // currently we just support harmonic analysis: complex data

    //get frequency
    MathParser * parser = domain->GetMathParser();
    parser->SetExpr( mHandle_, "f" );
    Double actFreq = parser->Eval( mHandle_ );

    //check solution type and compute factor
    Complex multVal = 0;

    // factor 0.5 is due to the fact, that the values are peak values
    if ( formulation_ == ACOU_PRESSURE ) {
      // solType = ACOU_PRESSURE;
      multVal = Complex(0.0, -0.5/ (2.0*PI*actFreq) );
    }
    else {
      // solType = ACOU_POTENTIAL;
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
                                results_[0]->fctType, isaxi_);

    // convert result object
    Result<TYPE> &  actRes =
      dynamic_cast<Result<TYPE>&>(*vals);
    EntityIterator regionIt = actRes.GetEntityList()->GetIterator();

    // resize vector
    Vector<TYPE> & actVal = actRes.GetVector();
    actVal.Resize( actRes.GetEntityList()->GetSize() );
    actVal.Init();

    // Loop over regions
    for( regionIt.Begin(); !regionIt.IsEnd(); regionIt++ )
    {
      SurfElemList actSDList(ptgrid_ );
      actSDList.SetRegion( regionIt.GetRegion() );
      EntityIterator it = actSDList.GetIterator();
      RegionIdType neighborRegionId = NO_REGION_ID;
      bool warned_before = false;
      
      // Try to determine neighbor region
      if ( surfNeighborRegions_.find(vals) != surfNeighborRegions_.end() )
      {
        neighborRegionId = surfNeighborRegions_[vals];
      }
      
      // Loop over all surface elements
      UInt counterElems = 0;
      for ( it.Begin(); !it.IsEnd(); it++, counterElems++ ) {

        const SurfElem * actSurfElem = it.GetSurfElem();

        // Determine, which volume element is the right neighbor for the
        // calculation;
        // our normal should point out of the correct neighbor volume element!
        if ( (actSurfElem->ptVolElem1 != NULL)
             && (actSurfElem->ptVolElem2 != NULL)
             && (actSurfElem->ptVolElem1->regionId
                == actSurfElem->ptVolElem2->regionId)
             && (actSurfElem->ptVolElem1->regionId != NO_REGION_ID) )
        {
          /* That means the surfRegion lies inside of a volume region. Why
           * would anyone do that? This is a problem, because we cannot
           * guarantee that the surface normals always point in the same
           * direction. But we don't want to abort the simulation, so we
           * just print a warning.
           */
          if ( !warned_before ) { // don't warn for every single element
            WARN("Calculation of result 'acouPower' on surface region '"
                << ptgrid_->GetRegion().ToString(regionIt.GetRegion())
                << "' might be wrong, because the surface is embedded in region '"
                << ptgrid_->GetRegion().ToString(actSurfElem->ptVolElem1->regionId)
                << "'. Therefore the surface normal vector could not be determined.");
            warned_before = true;
          }
        }
        if ( neighborRegionId != NO_REGION_ID ) {
          if ( (actSurfElem->ptVolElem1 != NULL) &&
               (neighborRegionId == actSurfElem->ptVolElem1->regionId) ) {
            ptVolElem = actSurfElem->ptVolElem1;
            normSign = 1.0;
          }
          else if ( (actSurfElem->ptVolElem2 != NULL) && 
                    (neighborRegionId == actSurfElem->ptVolElem2->regionId) ){
            ptVolElem = actSurfElem->ptVolElem2;
            normSign = -1.0;
          }
          else {
            WARN("Cannot calculate result 'acouPower', because region '"
                << ptgrid_->GetRegion().ToString(neighborRegionId)
                << "' is not a neighbor of surface region '"
                << ptgrid_->GetRegion().ToString(regionIt.GetRegion())
                << "'");
            break;
          }
        }
        else {
          if ( actSurfElem->ptVolElem1 != NULL
              && actSurfElem->ptVolElem2 == NULL )
          {
            ptVolElem = actSurfElem->ptVolElem1;
            normSign = 1.0;            
          }
          else if ( actSurfElem->ptVolElem1 == NULL // should never happen!?
                    && actSurfElem->ptVolElem2 != NULL )
          {
            ptVolElem = actSurfElem->ptVolElem2;
            normSign = -1.0;
          }
          else if ( actSurfElem->ptVolElem1 == NULL
                    && actSurfElem->ptVolElem2 == NULL ) {
            WARN("Cannot calculate result 'acouPower', because surface region '"
                << ptgrid_->GetRegion().ToString(regionIt.GetRegion())
                << "' has no neighboring volume region. Please check your mesh.")
            break;
          }
          else {
            WARN("Cannot calculate result 'acouPower', because surface region '"
                << ptgrid_->GetRegion().ToString(regionIt.GetRegion())
                << "' has two neighboring volume regions. Please provide the "
                << "desired volume region through the attribute 'neighborRegion'.")
            break;            
          }
        }

        normSign *= (Double) actSurfElem->normalSign;

        ptSurfElemFE = actSurfElem->ptElem;
        ptVolElemFE = ptVolElem->ptElem;

        const StdVector<UInt> & surfConnect = actSurfElem->connect;
        const StdVector<UInt> & volConnect = ptVolElem->connect;

        // get material information
        BaseMaterial * myMat = materials_[ptVolElem->regionId];
        myMat->GetScalar(density,DENSITY,Global::REAL);

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
          // surface integration coordinate for evaluating the
          // gradient on the surface of the volume element
          ptVolElemFE->GetLocalIntPoints4Surface(surfConnect, volConnect,
                                                 intPoints[iPt], lCoordVol);

          // Calculate jacobian gradient of surface element
          Double jacDet = ptSurfElemFE->CalcJacobianDetAtIp(iPt+1, CornerCoords, actSurfElem );

          // Calc gradient
          ElemList elList(ptgrid_);
          elList.SetElement( ptVolElem );
          EntityIterator it2 = elList.GetIterator();
          gradOp->CalcElemGradField(gradVal, it2, lCoordVol,1.0);

          // Calc global normal
          ptgrid_->CalcSurfNormal(normal, *actSurfElem);

          normal    *= normSign;
          gradNormal = normal * gradVal;

          LOG_DBG3(acoupde) << "CAP: se=" << actSurfElem->elemNum << " ip=" << iPt << " lCoordVol=" << lCoordVol.ToString() << " grad=" << gradVal.ToString() << " grad_n=" << gradNormal;

          // get solution of element and interpolate into integration point
          Vector<TYPE> elemSol;
          Vector<Double> shapeFnc;
          GetSolVecOfElement( elemSol, it2, results_[0]);
          ptVolElemFE->GetShFnc( shapeFnc, lCoordVol, ptVolElem );
          TYPE intPointSol = shapeFnc * elemSol;

          // get the conjugate complex value
          TYPE intPointSolConj = std::conj(intPointSol);

          if ( formulation_ == ACOU_PRESSURE ) {
            helpVal =  multVal * gradNormal * intPointSolConj * (1.0 / density );
          }
          else {
            helpVal = multVal * density  * gradNormal * intPointSolConj;
          }

          // integrate value
          elemPower += intWeights[iPt] * helpVal * jacDet;
        }

        actVal[regionIt.GetPos()] += elemPower;
      }
      if ( !it.IsEnd() ) { // calculation was aborted
        actVal[regionIt.GetPos()] = 0.0;
        break;
      }
    }
    delete gradOp;
  }


  template<class TYPE>
  void AcousticPDE::CalcAcouEnergy(shared_ptr<BaseResult> vals)
  {
    //get frequency
//     MathParser * parser = domain->GetMathParser();
//     parser->SetExpr( mHandle_, "f" );
//     Double actFreq = parser->Eval( mHandle_ );

    //check solution type and compute factor
    SolutionType solType;
    if ( formulation_ == ACOU_PRESSURE ) {
      solType = ACOU_PRESSURE;
    }
    else {
      solType = ACOU_POTENTIAL;
    }

    //some help variables
    Vector<TYPE> gradVal(dim_);

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
      // get material parameters
      Double density, compressibility;
      materials_[regionIt.GetRegion()]
        ->GetScalar(density,DENSITY,Global::REAL);
      materials_[regionIt.GetRegion()]
        ->GetScalar(compressibility,ACOU_BULK_MODULUS,Global::REAL);
      Double pressFactor = 1.0 / compressibility;

      ElemList actSDList(ptgrid_ );
      actSDList.SetRegion( regionIt.GetRegion() );
      EntityIterator elemIt = actSDList.GetIterator();
      
      RegionIdType actRegion = regionIt.GetRegion();
      Vector<Double>& acouVel = acouParticleVelocity_[actRegion];

      for ( elemIt.Begin(); !elemIt.IsEnd(); elemIt++ ) {
        BaseFE * ptElem = elemIt.GetElem()->ptElem;
        ptElem->SetAnsatzFct( results_[0]->fctType );
        const UInt nrIntPts = ptElem->GetNumIntPoints();
        const Vector<Double> & intWeights = ptElem->GetIntWeights();  
        
        Matrix<Double> ptCoord;
        ptgrid_->GetElemNodesCoord( ptCoord, 
                                    elemIt.GetElem()->connect,
                                    true );

        // get solution and time derivative of it
        Vector<TYPE> elemSol, elemSolTimeDeriv;
        GetSolVecOfElement(elemSol,elemIt,results_[0]);
        GetDerivSolVecOfElement(elemSolTimeDeriv,elemIt,results_[0]);

        // Loop over integration points
        TYPE elemEnergy = 0.0;
        Matrix<Double> xiDx, xiDxTransp;
        Vector<Double> shpFnc;
        Double jacDet;

        // sum over all integration points
        for( UInt iPt = 0; iPt < nrIntPts; iPt++ ) {

          //pressure at IP
          ptElem->GetShFncAtIp(shpFnc,iPt+1,elemIt.GetElem());
          TYPE pressureAtIp;
          if ( solType == ACOU_POTENTIAL ) {
            elemSolTimeDeriv.Inner(shpFnc,pressureAtIp); 
            pressureAtIp *= density;
          }
          else {
            elemSol.Inner(shpFnc,pressureAtIp); 
          }

          //compute global derivatives and  jacobian determinant
          ptElem->GetGlobDerivShFncAtIp(xiDx, iPt+1, ptCoord, jacDet, elemIt.GetElem());
          if (isaxi_) {
            Vector<Double> CoordAtIP;
            CoordAtIP = ptCoord * shpFnc;
            jacDet *=  2 * PI * CoordAtIP[0];
          }

          //compute square of particle velocity
          xiDx.Transpose(xiDxTransp);
          Double factorGrad;
          if ( solType == ACOU_PRESSURE ) {
            for ( UInt i=0; i<dim_; i++ ) 
              gradVal[i] = acouVel[elemIt.GetPos()*dim_ + i];

            gradVal.Inner(gradVal,factorGrad);
          }
          else {
            gradVal = xiDxTransp * elemSol; 
            gradVal.Inner(gradVal,factorGrad); 
          }

          elemEnergy += 0.5*intWeights[iPt]*jacDet*
            ( density*factorGrad + pressFactor*pressureAtIp*pressureAtIp );
        }

        actVal[regionIt.GetPos()] += elemEnergy;
      }
    }
  }

  void AcousticPDE::CalcVelFromPressure( shared_ptr<BaseResult> vals )
  {
    Vector<Double> lCoord;
    Vector<Double> tempGrad;
    NodeStoreSol<Double> & solhelp = dynamic_cast<NodeStoreSol<Double>&>(*sol_);

    // Create operator for gradient computation of solution
    GradientFieldOp<Double> * gradOp =
      new GradientFieldOp<Double>(ptgrid_, this, eqnMap_, solhelp,
                                results_[0]->fctType, isaxi_);

    // convert result object
    Result<Double> &  actRes =
      dynamic_cast<Result<Double>&>(*vals);
    EntityIterator regionIt = actRes.GetEntityList()->GetIterator();
    
    Double dt = TS_alg_->GetTimeStep() ;

    // Loop over regions
    for( regionIt.Begin(); !regionIt.IsEnd(); regionIt++ ) {
      // get material parameters
      Double density, factor;
      materials_[regionIt.GetRegion()]
        ->GetScalar(density,DENSITY,Global::REAL);

      factor = dt / density;

      ElemList actSDList(ptgrid_ );
      actSDList.SetRegion( regionIt.GetRegion() );
      EntityIterator elemIt = actSDList.GetIterator();

      //get acoustic particle velocity array for current region     
      RegionIdType actRegion = regionIt.GetRegion();
      Vector<Double>& vel = acouParticleVelocity_[actRegion];

      for ( elemIt.Begin(); !elemIt.IsEnd(); elemIt++ ) {
        elemIt.GetElem()->ptElem->GetCoordMidPoint(lCoord);
        //computes -(nabla p); minus included!
        gradOp->CalcElemGradField( tempGrad, elemIt, lCoord, 1); 

        //scale the gardient
        tempGrad *= factor;

        for(UInt iDim = 0; iDim < dim_; iDim++ ) {
          vel[elemIt.GetPos()*dim_ + iDim] += tempGrad[iDim];
        }
      }
      //      std::cout << "elIdx: " << elIdx << std::endl;
    }
    //    std::cout << "elIdxEnd: " << elIdx << std::endl;
    delete gradOp;
    
  }



  // ***********************************************************************
  //   Obtain information on desired output quantities from parameter file
  // ***********************************************************************
  void AcousticPDE::DefineAvailResults()
  {
    // === NODE POTENTIAL / PRESSURE (Primary Unknown) ===
    // check if problem is lagrange or legendre
    std::string approxType = myParam_->Get("type")->As<std::string>();

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
    }else if(  approxType == "spectral" ) {
      UInt order = myParam_->Get("order")->As<UInt>();
      shared_ptr<SpectralFct> fct(new SpectralFct);
      res1->definedOn = ResultInfo::PFEM;
      fct->SetOrder(order);
      res1->fctType = fct;

    } else {
      UInt order = myParam_->Get("order")->As<UInt>();

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

    // for time domain PML we need auxillary variables
    // check for time domain PML
    if ( analysistype_ == TRANSIENT ) { 
      bool isPML = false;

      for (UInt actSD = 0; actSD < subdoms_.GetSize(); actSD++) {
        // get current region name and get grip of paramNode
        RegionIdType actRegion = subdoms_[actSD];
        // Check for Perfectly matchec layers
        if ( dampingList_[actRegion] == PML )
          isPML = true;
      }

      if ( isPML ) {
        // === vectorial auxillary variable ===
        shared_ptr<ResultInfo> res2(new ResultInfo); 
        res2->resultType = ACOU_PMLAUXVEC;
        if( subType_ == "3d" ) {
          res2->dofNames = "phix", "phiy", "phiz";
        } 
        else if ( subType_ == "axi" ) {
          res2->dofNames = "phir", "phiz";
        } 
        else if( subType_ == "plane") {
          res2->dofNames = "phix", "phiy";
        } 
        
        res2->unit = "m/s";
        res2->entryType = ResultInfo::VECTOR;
        
        std::string approxTypeV = myParam_->Get("type")->As<std::string>();
        if ( approxTypeV == "lagrange" ) {
          shared_ptr<AnsatzFct> fctV(new LagrangeFct);
          res2->fctType = fctV;
          res2->definedOn = ResultInfo::NODE;
        } 
        else {
          UInt order = myParam_->Get("order")->As<UInt>();
          
          // Create new resultDof object
          shared_ptr<LegendreFct> fct(new LegendreFct);
          fct->SetIsoOrder( order );
          //fct->order_ = order;
          res2->definedOn = ResultInfo::PFEM;
          res2->fctType = fct;
        }
        availResults_.insert( res2 );
        results_.Push_back( res2 );

        if( subType_ == "3d" && isAPML_ == false ) {
          // === scalar auxillary variable ===
          shared_ptr<ResultInfo> res3(new ResultInfo); 
          res3->resultType = ACOU_PMLAUXSCALAR;
          res3->dofNames = "psi";
          res3->unit = "";
          if ( approxType == "lagrange" ) {
            shared_ptr<AnsatzFct> fct(new LagrangeFct);
            res3->definedOn = ResultInfo::NODE;
            res3->fctType = fct;
          } else {
            UInt order = myParam_->Get("order")->As<UInt>();
            
            // Create new resultDof object
            shared_ptr<LegendreFct> fct(new LegendreFct);
            fct->SetIsoOrder( order );
            //fct->order_ = order;
            res3->definedOn = ResultInfo::PFEM;
            res3->fctType = fct;
          }
          res3->entryType = ResultInfo::SCALAR;
          availResults_.insert( res3 );
          results_.Push_back( res3 );  
        }
      }
    }

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
    rhs->resultType = ACOU_RHS_LOAD;
    rhs->dofNames = "";
    rhs->unit = "";
    rhs->entryType = res1->entryType;
    rhs->definedOn = res1->definedOn;
    rhs->fctType = res1->fctType;
    availResults_.insert( rhs );


    // === ACOU_SURFINTENSITY ===
    shared_ptr<ResultInfo> intens(new ResultInfo);
    intens->resultType = ACOU_SURFINTENSITY;
    intens->SetVectorDOFs(dim_, isaxi_);
    intens->unit = "W/m^2";
    intens->entryType = ResultInfo::VECTOR;
    intens->definedOn = ResultInfo::SURF_ELEM;
    intens->fctType = shared_ptr<ConstFct>(new ConstFct() );
    availResults_.insert( intens );

    // === ACOU_INTENSITY ===
    shared_ptr<ResultInfo> intensity(new ResultInfo);
    intensity->resultType = ACOU_INTENSITY;
    intensity->SetVectorDOFs(dim_, isaxi_);
    intensity->unit = "W/m^2";
    intensity->entryType = ResultInfo::VECTOR;
    intensity->definedOn = ResultInfo::ELEMENT;
    intensity->fctType = shared_ptr<ConstFct>(new ConstFct() );
    availResults_.insert( intensity );


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

    // ==== GRAD_ACOU_SOLUTION ===
    shared_ptr<ResultInfo> grad_sol(new ResultInfo);
    grad_sol->resultType = GRAD_ACOU_SOLUTION;
    grad_sol->SetVectorDOFs(dim_, isaxi_);
    grad_sol->unit = "";
    grad_sol->entryType = ResultInfo::VECTOR;
    grad_sol->definedOn = ResultInfo::NODE;
    grad_sol->fctType = shared_ptr<ConstFct>(new ConstFct() ); // ??? - copy and paste! !?
    availResults_.insert( grad_sol );

    // === ACOU_ENERGY ===
    shared_ptr<ResultInfo> energy(new ResultInfo);
    energy->resultType = ACOU_ENERGY;
    energy->dofNames = "";
    energy->unit = "Ws";
    energy->entryType = ResultInfo::SCALAR;
    energy->definedOn = ResultInfo::REGION;
    energy->fctType = shared_ptr<ConstFct>(new ConstFct() );
    availResults_.insert( energy );

    // === ACOU PSEUDO DENSITY for SIMP ===
    shared_ptr<ResultInfo> acouPD(new ResultInfo);
    acouPD->resultType = ACOU_PSEUDO_DENSITY;
    acouPD->dofNames = "";
    acouPD->unit = "";
    acouPD->entryType = ResultInfo::SCALAR;
    acouPD->definedOn = ResultInfo::ELEMENT;
    acouPD->fctType = shared_ptr<ConstFct>(new ConstFct() );
    availResults_.insert(acouPD);

    // ===================================
    // Check for non-conforming interfaces
    // ===================================
    StdVector<std::string> ncIfaceNames, ncIfaceNamesForPDE;
    StdVector<RegionIdType> ncIfaceIds;

    LOG_DBG2(acoupde) << "NonMatching: Checking if nonconforming "
                      << "interfaces of PDE exist in domain.";

    PtrParamNode domainNCIfaceListNode;
    domainNCIfaceListNode = param->Get("domain")->Get("ncInterfaceList", ParamNode::PASS);

    if(domainNCIfaceListNode)
    {
      PtrParamNode ncInterfaceListNode =
        param->GetByVal("sequenceStep", std::string("index"), sequenceStep_)
        ->Get("pdeList")->Get("acoustic")->Get("ncInterfaceList", ParamNode::PASS);
      ParamNodeList pdeNCIfaceNodes;

      if(ncInterfaceListNode)
      {
        pdeNCIfaceNodes = ncInterfaceListNode->GetList("ncInterface");

        for (UInt i = 0; i < pdeNCIfaceNodes.GetSize(); i++) {
          std::string pdeIfaceName = pdeNCIfaceNodes[i]->Get("name")->As<std::string>();
          std::string domainIfaceName;

          PtrParamNode domainIfaceNode = domainNCIfaceListNode->GetByVal("ncInterface",
              "name", pdeIfaceName, ParamNode::PASS);
          if(!domainIfaceNode)
          {
            LOG_DBG2(acoupde) << "NonMatching: Nonconforming "
            << "interface '" << ncIfaceNames[i]
                                             << "' does not exist in domain.";

            EXCEPTION( "ncInterface referenced from PDE not defined in domain!");
          }

          ncIfaceNamesForPDE.Push_back(pdeIfaceName);
        }
        ptgrid_->GetRegion().Parse( ncIfaceNamesForPDE, ncIfaceIds );

        for (UInt i = 0; i < ncIfaceIds.GetSize(); i++) {
          ncIFaces_.Push_back(ncIfaceIds[i]);
        }

        // In the case of the presence of non-conforming interfaces,
        // a second resultdof object has to be created, which describes the
        // Lagrange multiplier
        if( ncIFaces_.GetSize() > 0 ) {
          LOG_DBG2(acoupde) << "NonMatching: Defining new ResultDof "
                            << "Lagrange Multiplier (LM).";
          shared_ptr<ResultInfo> lagr ( new ResultInfo );
          lagr->resultType = LAGRANGE_MULT;
          lagr->dofNames = "l";
          lagr->fctType = results_[0]->fctType;
          lagr->definedOn = results_[0]->definedOn;
          results_.Push_back( lagr );
        }
      }

    }
  }

  void AcousticPDE::ReadFlowData() {

    // check if bcsNode is present
    PtrParamNode bcsNode = myParam_->Get("bcsAndLoads", ParamNode::PASS );
    if( !bcsNode) return;

    // get nodes with flowData
    ParamNodeList flowNodes = bcsNode->GetList("flowData");
    for( UInt i = 0; i < flowNodes.GetSize(); i++ ) {
      std::string regionName = flowNodes[i]->Get("region")->As<std::string>();
      RegionIdType regionId = ptgrid_->GetRegion().Parse(regionName);

      // store information about flow
      regionFlowNodes_[regionId] = flowNodes[i];
    }

  }


  void AcousticPDE::ReadDataDampLayer( std::string& dampingTypePML,
                                       Vector<Double>& mPoint,
                                       Double& dampFactor,
                                       Double& dampFactorMax,
                                       Double& startRadius,
                                       Double& endRadius,
                                       PtrParamNode actNode ) {




    StdVector<std::string> stringVal;

    // Construct vectors for restricted parameter search
    actNode->GetValue( "type", dampingTypePML );

    //get the center point
    mPoint.Resize(dim_);

    //xM, yM,
    actNode->GetValue( "xM", mPoint[0] );
    actNode->GetValue( "yM", mPoint[1] );

    if ( dim_ == 3 ) {
      //zM
      actNode->GetValue( "zM", mPoint[2] );
    }

    //start radius, end radius, dampFactor
    actNode->GetValue( "radiusStart", startRadius );
    actNode->GetValue( "radiusEnd", endRadius );
    actNode->GetValue( "dampFactor", dampFactor );

    //get maximal damping factor (saturation)
    actNode->GetValue( "dampFactorMax", dampFactorMax );
  }

  void AcousticPDE::PreparePDE4Computation()  {
    //computation of acoustic energy in case of pressure formulation

    ResultMap::iterator it = resultLists_.begin();

    bool isAcouEnergy = false;
    // iterate over all results
    for( ; it != resultLists_.end(); it++ ) {
      ResultList & actList = it->second;

      // iterate over all solutions for each result type
      for( UInt i = 0; i < actList.GetSize(); i++ ) {
        if ( actList[i]->GetResultInfo()->resultType == ACOU_ENERGY 
             && formulation_ == ACOU_PRESSURE ) {
          isAcouEnergy = true;
        }
      }
    }

    if ( isAcouEnergy ) { 
      // get memory to store acoustic particle velocity
      RegionIdType actRegion;
      std::map<RegionIdType, BaseMaterial*>::iterator it;
      for ( it = materials_.begin(); it != materials_.end(); it++ ) {
        actRegion = it->first;
        acouParticleVelocity_[actRegion].Resize(numElems_*dim_);
        acouParticleVelocity_[actRegion].Init();
      }
    }
  }

} // end of namespace

// explicit template instantiation for GCC compiler
#ifdef __GNUC__

template
void AcousticPDE::CalcGradSurfaceElement<Complex>(
    shared_ptr<BaseResult> vals, SolutionType solType,
    bool apply_normal,
    StdVector<Vector<Complex> >& grad_out );

#endif
