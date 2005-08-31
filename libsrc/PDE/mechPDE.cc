#include "mechPDE.hh"

#include <sstream>
#include <iomanip>

#include "Forms/forms_header.hh"
#include "Forms/linElastInt.hh"
#include "Forms/massInt.hh"
#include "DataInOut/writeresults.hh"
#include "Driver/assemble.hh"
#include "newmark.hh"
#include "newmarkFracDampMech.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "Driver/solveStepMech.hh" 
#include "CoupledPDE/pdecoupling.hh"
#include "Domain/domain.hh"
#include "Utils/coordSystem.hh"

namespace CoupledField
{

  MechPDE::MechPDE(Grid * aptgrid, TimeFunc *aptTimeFunc, WriteResults *aptOut)
    :SinglePDE(aptgrid, aptOut, aptTimeFunc), 
     fracDamping_(FALSE),
     preStressVal_(0.0),
     lambdaMat(NULL),
     mueMat(NULL)

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
    for ( UInt k = pressFnc_.GetSize(); k < pressSurf_.GetSize(); k++ )
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

  void MechPDE::ReadDampingInformation( )
  {
    ENTER_FCN( "MechPDE::ReadDampingInformation" );
    

    
    fracMemory_ = 0;
    Boolean identical = TRUE; // i.e. same type of damping for all regions
    Integer firstFrac=-1;

    // Construct vectors for restricted search parameter
    StdVector<std::string> keyVec;
    StdVector<std::string> attrVec;
    StdVector<std::string> valVec;

    for (UInt k = 0; k < subdoms_.GetSize(); k++) {
      
      RegionIdType actRegion = subdoms_[k];

      std::string actRegionName;
      actRegionName = ptgrid_->RegionIdToName( actRegion );
      keyVec = "mechanic" , "region" , "damping" , "type";
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
                      actRegionName.c_str() );
      }
      else if (dampInfo[0] == "rayleigh") {
        dampingList_[actRegion] = RAYLEIGH;
        Info->PrintF( pdename_, 
                      "      * RAYLEIGH damping for region: %s\n",
                      actRegionName.c_str() );
      }
      else if (dampInfo[0] == "fractional") {
        dampingList_[actRegion] = FRACTIONAL;
        fracDamping_ = TRUE;
        Info->PrintF( pdename_, 
                      "      * FRACTIONAL damping for region: %s\n", 
                      actRegionName.c_str() );
        
        // Find first region containing fractional damping
        if ( firstFrac < 0 )
          firstFrac = k;
        
        keyVec = "mechanic" , "region" , "damping" , "fracAlg";
        StdVector<std::string> fracAlg;
        params->GetList( keyVec, attrVec, valVec, fracAlg );
        
        keyVec = "mechanic" , "region" , "damping" , "fracMemory";
        StdVector<UInt> fracMem;
        params->GetList( keyVec, attrVec, valVec, fracMem );
        
        keyVec = "mechanic" , "region" , "damping" , "interpolation";
        StdVector<std::string> interpol;
        params->GetList(  keyVec, attrVec, valVec, interpol );
        
        
        if( fracAlg.IsEmpty() || fracMem.IsEmpty() || interpol.IsEmpty() ) {
          (*error) << "Specify attributes fracAlg, fracMemory " 
                   << "and interpolation!";
          Error( __FILE__, __LINE__ ); 
        }
        else if  ( fracAlg[0] == "gl" ) {
          
          // Include fracAlg and interpolation info in dampingList
          Info->PrintF( "", "\t\t\t using Gruenwald-Letnikov algorithm,\n");
          if (interpol[0] == "no" )
            dampingList_[actRegion] = FRACTIONAL_GL;
          else {
            Error("Till now no interpolation is allowed in mechanics fractional damnping!",__FILE__,__LINE__);
          }
        }
//         else if ( fracAlg[0] == "blank" ) {
          
//           Info->PrintF( "", "\t\t\t using Blanks algorithm,\n");
//           if (interpol[0] == "no" )
//             dampingList_[actRegion] = FRACTIONAL_BLANK;
//           else {
//             dampingList_[actREgion] = FRACTIONAL_BLANK_INT;
//             Info->PrintF("", 
//                          "\t\t\t linear interpol. of single past values\n\n");
//           }
//         }
        
        // up to now take maximum of fracMemory
        if ( fracMem[0] > fracMemory_ )
          fracMemory_ = fracMem[0];
      }
      
    }
    
    // Check, if all entries are identical
    for ( UInt i = 1; i < dampingList_.size(); i++ ) {
      if ( dampingList_[subdoms_[i-1]] != dampingList_[subdoms_[i]] ) {
        identical = FALSE;
        break;
      }
    }


    // Fractional damping can only be enabled, if all regions are damped
    // this way. Oterhwise an error is thrown.
    if ( fracDamping_ == TRUE ) {
      
      if ( identical == TRUE ) {
        
        fracDamping_ = TRUE;
        Info->PrintF(pdename_, "Memory size for fractional damping  is: %d\n",
                     fracMemory_ );
      } else {
        
        (*error) << "Fractional damping can only be used if it is applied for "
                 << "ALL regions of the mechanical domain!\n"
                 << "Please check your parameter file!";
        Error( __FILE__, __LINE__ );
      }
        
        
    }
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
      for ( UInt k = 1; k < nonLinRegion.GetSize(); k++ ) {
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
      for ( UInt i = 0; i < reducedIntegration_.GetSize(); i++ ) {
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
    for (UInt actSD = 0; actSD < subdoms_.GetSize(); actSD++)
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

	    actIntDescr1->SetPDEIds(this, this);


            //check for damping
            if ( dampingList_[subdoms_[actSD]] == RAYLEIGH ) {
              actIntDescr1->SetSecondaryMat(DAMPING, actSDMat.GetDampingBeta(),analysistype_);
            }
            assemble_->AddIntegrator(actIntDescr1, subdoms_[actSD]);


            BaseForm * bilinearStiff2 = GetStiffIntegrator(*lambdaMat);
            IntegratorDescriptor * actIntDescr2 =
              new IntegratorDescriptor(bilinearStiff2, STIFFNESS);

	    actIntDescr2->SetPDEIds(this, this);


            //for this integrator, we need reduced integration
            //see Hughes, pp.219
            actIntDescr2->SetReducedInt();

            //check for damping
            if ( dampingList_[subdoms_[actSD]] == RAYLEIGH ) {
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

	    actIntDescr->SetPDEIds(this, this);


            //check for damping
            if ( dampingList_[subdoms_[actSD]] == RAYLEIGH ) {
              actIntDescr->SetSecondaryMat(DAMPING, actSDMat.GetDampingBeta(),analysistype_);
            }
	    
            assemble_->AddIntegrator(actIntDescr, subdoms_[actSD]);


            //for prestressing
            for ( UInt preStr=0; preStr<preStressDomain_.GetSize(); preStr++ ) {
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
		actIntDescrPre->SetPDEIds(this, this);
	    
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

	    IntegratorDescriptor * stiffNL1Descr = 
	      new IntegratorDescriptor(nLinPart1, STIFFNESS, nonLin_);

	    stiffNL1Descr->SetPDEIds(this, this);
	    assemble_->AddIntegrator(stiffNL1Descr, subdoms_[actSD]);

	    IntegratorDescriptor * stiffNL2Descr = 
	      new IntegratorDescriptor(nLinPart2, STIFFNESS, nonLin_);

	    stiffNL2Descr->SetPDEIds(this, this);
	    assemble_->AddIntegrator(stiffNL2Descr, subdoms_[actSD]);

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
	actIntDescr->SetPDEIds(this, this);

        //check for damping (mass part)

        if ( dampingList_[subdoms_[actSD]] == RAYLEIGH ) {
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
    Boolean nonlin = FALSE;
    for (UInt actSF = 0; actSF < pressSurf_.GetSize(); actSF++) {
      BaseForm * rhsSrcSurf = new PressureLinForm(pressVals_[actSF], isaxi_);
      assemble_->AddRhsSrcSurfIntegrator(rhsSrcSurf, pressSurf_[actSF], pressFnc_[actSF],
                                         nonlin);
    }
    
    // Trigger reading of region load
    DefineRegionLoads();
    

  }


  BaseForm *
  MechPDE::GetStiffIntegrator(MaterialData& actSDMat, Boolean reducedInt)
  {

    ENTER_FCN( "MechPDE::GetStiffIntegrator" );
  
    BaseForm * bilinearStiff = NULL;

    if( fracDamping_ == FALSE )
      {
        if (subType_ == "planeStrain")
          bilinearStiff = new mechPlainStrainInt(actSDMat);
        else if (subType_ == "axi")
          bilinearStiff = new mechAxiInt(actSDMat);
        else if (subType_ == "3d")
          bilinearStiff = new mech3DInt(actSDMat);
        else 
          Info->Error("Unknown subtype in mech PDE! ",__FILE__,__LINE__);
      } else{
      
      StdVector<std::string> keyVec, attrVec, valVec;
      keyVec = "transient", "firstDt";
      attrVec = "tag";
      valVec  =  bcSequenceTag_;
      Double dt = 0.0;
      params->Get(keyVec,attrVec,valVec,dt);
      bilinearStiff = new LinViscoElastInt(actSDMat,GetSubType(), "modifiedStiffness",dt );
    }
    
    return bilinearStiff;
  }
  
  void MechPDE::DefineRegionLoads() {
    ENTER_FCN ( "MechPDE::DefineRegionLoads" );
    
    StdVector<std::string> keyVec, attrVec, valVec;
    StdVector<std::string> tempNames, names, dofs, dynamics, refCoord, type;
    StdVector<RegionIdType> regionIds;
    StdVector<UInt> vecComp;
    StdVector<Double> loadVec;
    Vector<Double> unitLoad(dim_), totLoad(dim_), tempLoad(dim_);
    MechVolForceInt * forceInt = NULL;
    Integer index = -1;
    UInt locDof = 0;
    Double volume = 0.0;
    std::ostringstream out;


    // get names of all regions with region loads
    keyVec = "mechanic", "bcsAndLoads", "regionLoad", "name";
    attrVec = "", "tag", "";
    valVec  = "", bcSequenceTag_, "";
    params->GetList(keyVec,attrVec,valVec,tempNames);

    // Now sort the names and remove double entries
    for (UInt i = 0; i < tempNames.GetSize(); i++) {
      index = names.Find(tempNames[i]);
      if ( index == -1) {
        names.Push_back(tempNames[i]);
      }
    }

    ptgrid_->RegionNameToId(regionIds, names);

    if ( regionIds.GetSize() > 0 ) {
      Info->PrintF(pdename_, "The following regions have a region load:\n\n");
      out.clear();
      out << std::setw(15) << "name" << " | " 
          << std::setw(15) << "refCoordSys" << " | "
          << std::setw(15) << "dynamics" << " | "
          << std::setw(11) << "volume" << " | "
          << std::setw(11) << "tot. load" << " | "
          << std::setw(11) << "unit load" <<std::endl;
      Info->PrintF(pdename_, out.str().c_str());
      out.str("");
      out << std::setw(90) << std::setfill('-') << "" 
          << std::setfill(' ') << std::endl;
      Info->PrintF(pdename_, out.str().c_str());
      out.str("");
    }
                   


    // loop over all regionnames
    for (UInt i = 0; i < names.GetSize(); i++) {

      // restrict search to current region name
      attrVec = "", "tag", "name";
      valVec  = "", bcSequenceTag_, names[i];

      // get value
      keyVec = "mechanic", "bcsAndLoads", "regionLoad", "value";
      params->GetList(keyVec, attrVec, valVec, loadVec);
      
      // get dynamics
      keyVec = "mechanic", "bcsAndLoads", "regionLoad", "dynamics";
      params->GetList(keyVec, attrVec, valVec, dynamics);

      // get dofs
      keyVec = "mechanic", "bcsAndLoads", "regionLoad", "dof";
      params->GetList(keyVec, attrVec, valVec, dofs);

      // get coordinate system
      keyVec = "mechanic", "bcsAndLoads", "regionLoad", "refCoordSys";
      params->GetList(keyVec, attrVec, valVec, refCoord);

      // get load type (total / unit)
      keyVec = "mechanic", "bcsAndLoads", "regionLoad", "type";
      params->GetList(keyVec, attrVec, valVec, type);

      // check if all vectors have the same length
      if ( loadVec.GetSize() != dynamics.GetSize() ||
           loadVec.GetSize() != dofs.GetSize() ||
           loadVec.GetSize() != refCoord.GetSize() || 
           loadVec.GetSize() != type.GetSize() ) {
        (*error) << "MechPDE::DefineRegionLoads: Inconsistent definition for "
                 << "region '" << names[i] <<"'!\n"
                 << "Please correct parameter file!";
        Error( __FILE__, __LINE__ );
      }
      
      // check if all entries for dynamics, refCoord and type
      // are the same
      for (UInt k=0; k<dynamics.GetSize(); k++) {
        if (dynamics[k] != dynamics[0] ||
            refCoord[k] != refCoord[0] ||
            type[k] != type[0] ) {
          (*error) << "MechPDE::DefineRegionLoads: The region load on region "
                   << names[i] << " has not for all dofs the same entry for "
                   << "dynamics, refCoord or type (total/unit)!";
          Error( __FILE__, __LINE__ );
        }
      }

      // now create local load vector
      unitLoad.Init();
      tempLoad.Init();
      for (UInt iDim=0; iDim < loadVec.GetSize(); iDim++) {
        locDof = domain->GetCoordSystem(refCoord[iDim])->
          GetVecComponent(dofs[iDim]);
        tempLoad[locDof-1] = loadVec[iDim];
      }

      // if the load is for a complete region, divide it by the
      // volume to obtain the unit load in f / m^2
      volume = ptgrid_->CalcVolumeOfRegion(regionIds[i], isaxi_);
      if ( type[0] == "total" ) {
        totLoad = tempLoad;
        unitLoad = tempLoad / volume;
      } else {
        totLoad =  tempLoad * volume;
        unitLoad = tempLoad;
      }
        
      // create linearform and add to assemble
      forceInt = new MechVolForceInt(dim_, isaxi_);
      forceInt->SetVolForceVector(unitLoad, 
                                  domain->GetCoordSystem(refCoord[0]));
      assemble_->AddRhsSrcIntegrator( forceInt, regionIds[i],
                                      dynamics[0], nonLin_ );

      // write logging information into info file
      for (UInt k=0; k<dim_; k++) {
        out.str("");
        if ( k == 0) {
          out << std::setw(15) << names[i] << " | " 
              << std::setw(15) << refCoord[0] << " | "
              << std::setw(15) << dynamics[0] << " | "
              << std::setw(11) << volume << " | ";
        } else {
          out << std::setw(15) << "" << " | " 
              << std::setw(15) << "" << " | "
              << std::setw(15) << "" << " | "
              << std::setw(11) << "" << " | ";
              
        }
        
        out << std::setw(11) << totLoad[k] << " | "
            << std::setw(11) << unitLoad[k] << std::endl;

        Info->PrintF(pdename_,out.str().c_str());
      }
      Info->PrintF(pdename_,"\n");
    }
    
  }
  
  void MechPDE::DefineSolveStep()
  {
    ENTER_FCN( "MechPDE::DefineSolveStep" );

    solveStep_ = new SolveStepMech(*this);
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


    for(UInt actRow=0; actRow<3; actRow++)
      {
        for(UInt actCol=0; actCol<3; actCol++)
          (*lMechMat)[actRow][actCol] = lambda;

        (*mueMechMat)[actRow][actRow] = 2*mue;
      }

    for(UInt actRow=3; actRow<6; actRow++)
      (*mueMechMat)[actRow][actRow] = mue;  

    std::cout << "LambadMat:\n" << *lMechMat << std::endl;
    std::cout << "MuMat:\n" << *mueMechMat << std::endl;

  }



  // ======================================================
  // COUPLING SECTION
  // ======================================================


  void MechPDE::InitCoupling(PDECoupling * Coupling)
  {
    ENTER_FCN( "MechPDE::InitCoupling" );
  
    isIterCoupled_ = TRUE;
    ptCoupling_   = Coupling;
  
    for (UInt i=0; i<ptCoupling_->GetNumOutputCouplings(); i++)
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
          }
      }

  }





  void MechPDE::CalcOutputCoupling()
  {
    ENTER_FCN( "MechPDE::CalcOutputCoupling" );

    UInt dof = 0;
    SolutionType quantity;
    StdVector<UInt> * couplingnodes = NULL;
    StdVector<Elem*> * couplingElems = NULL;
    CFSVector * temp_values = NULL;
    Vector<Double> * values;
    StdVector<MaterialData*> * materials = NULL;
    StdVector<std::string> outputRegions;
  

    // at first, check if this PDE is iterative coupled
    if (isIterCoupled_ == FALSE)
      return;

    // loop over all output coupling quantities
    for (UInt i=0; i<ptCoupling_->GetNumOutputCouplings(); i++)
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
                                         StdVector<UInt>& couplingNodes,
                                         Vector<Double> & elemCouplingSols,
                                         UInt couplingdof )
  {
    ENTER_FCN( "MechPDE::CalcAcousticCouplingRHS" );

    Matrix<Double> ptCoord, elemMat;
    Vector<Double> normal, sol, forceOnElem, nSol;
    Elem * ptVolElem;
    Double sign = 0.0;
    Double density = 0.0;
    Integer matIndex = -1;

    elemCouplingSols.Init(0.0);
  
    for (UInt actElem=0; actElem<couplingElems->GetSize(); actElem++)
      {
        // Perform cast from volume element to surface element, since
        // mech-acou coupling makes only sense on surface elements
        SurfElem * actCoupleElem = 
          dynamic_cast<SurfElem*> ((*couplingElems)[actElem]);
        

        BaseFE * ptElem = actCoupleElem->ptElem;
        StdVector<UInt> & connecth = (*couplingElems)[actElem]->connect;
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
      

        for (UInt actNode=0; actNode < connecth.GetSize(); actNode++)
          for (UInt actDof=0; actDof<dofspernode_; actDof++)
            nSol[actNode] += sol[actDof + actNode*dofspernode_] * normal[actDof];


        Vector<Double> forceOnElem = elemmat * nSol;  
      
        for (UInt actNode=0; actNode<ptCoord.GetSizeRow(); actNode++)
          {
            UInt nodePos = 0;
          
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

    UInt rhsSize = eqnData_->GetNumEQNs() *
      eqnData_->GetNumDofsPerEQN();

    if ( fracDamping_ == FALSE ) { 
      if (effectiveMass_)  
        TS_alg_ = new NewmarkEffMass( algsys_, rhsSize );
      else
        TS_alg_ = new Newmark( algsys_, rhsSize );
    }
    else {
      if (effectiveMass_ == FALSE)
        TS_alg_ = new NewmarkFracDampMech( algsys_, rhsSize, pdeId_,
                                           eqnData_, ptgrid_, this, 
                                           subdoms_, dampingList_ );
      else
        Error("This needs to be implemented!",__FILE__,__LINE__);
    }


  }



  // ======================================================
  // POSTPROCESSING SECTION
  // ======================================================


  void MechPDE::WriteResultsInFile(const UInt kstep,
                                   const Double asteptime,
                                   UInt stepOffset,
                                   Double timeOffset)
  {
    ENTER_FCN( "MechPDE::WriteResultsInFile" );

    NodeStoreSol<Double> sol_der1Array, sol_der2Array;
    NodeStoreSol<Double> * solTransient;
    NodeStoreSol<Complex> * solHarmonic;

    Double actTime = asteptime + timeOffset;
    UInt actStep = kstep + stepOffset;
    
#ifdef WRITE_RHS
    NodeStoreSol<Double> rhs;
    rhs.SetNumSolutions(1);
    rhs.SetNumNodes(numPDENodes_);
    rhs.SetSolutionType(ACOU_RHSVAL);
    rhs.SetNumDofs(dim_);
    rhs.SetPtrEQNData(eqnData_, ptgrid_);
    rhs.Init(0.0);
    
    Double *ptRHS;
    algsys_->GetRHSVal( ptRHS );
    rhs.CopyFromAlgSysDataPointer(ptRHS);
    outFile_->WriteNodeSolutionTransient(rhs, actStep, actTime);
#endif

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
        outFile_->WriteNodeSolutionHarmonic(*solHarmonic,  actStep,
                                            actTime, complexFormat_);
      if (saveSolHist_ == TRUE)
        outFile_->WriteNodeHistoryHarmonic(*solHarmonic,  actStep, 
                                           actTime, complexFormat_);
    
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
  
    // -------------------------
    //  Determine nodal results
    // -------------------------
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

    // ---------------------------
    //  Determine element results
    // ---------------------------
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
      for ( UInt k = 0; k < regionNames.GetSize(); k++ ) {
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
      for ( UInt k = 0; k < regionNames.GetSize(); k++ ) {
        Info->PrintF( pdename_, " %s\n", regionNames[k].c_str() );
      }
    }

    // -------------------------
    //  Determine nodal history
    // -------------------------
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
      for ( UInt k = 0; k < saveNodeHist.GetSize(); k++ ) {
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
      for ( UInt k = 0; k < saveNodeHist.GetSize(); k++ ) {
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
      for ( UInt k = 0; k < saveNodeHist.GetSize(); k++ ) {
        Info->PrintF( pdename_, " %s\n", saveNodeHist[k].c_str() );
      }
    }

    // ---------------------------
    //  Determine element history
    // ---------------------------
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

    // ---------------------------
    //  Determine special results
    // ---------------------------
    params->GetList( "name", regionNames, pdename_, "volumeAboveDefSurf" );
    ptgrid_->RegionNameToId( volAboveDefSurfRegions_, regionNames );
    params->GetList( "dof", volAboveDefSurfDir_, pdename_, 
		     "volumeAboveDefSurf" );

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
        Info->Error( "StressOp: Unknown subtype in mech PDE!", __FILE__,
                     __LINE__);  

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
      for (UInt isd=0; isd<subdoms_.GetSize(); isd++) {
        
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
        for (UInt iel=0; iel< elemssd.GetSize(); iel++) {
          UInt pdeElem = eqnData_->Mesh2PDEElem(elemssd[iel]->elemNum);
          elemssd[iel]->ptElem->GetCoordMidPoint(intPoint);
          //set element pointer
          BaseFE * ptEl = elemssd[iel]->ptElem;
          stress->SetElemPtr(ptEl);
          
          //set element solution        
          Matrix<Double> elSol;
          StdVector<UInt> connecth = elemssd[iel]->connect;
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

    // Compute volume above deformed surfaces
    if ( volAboveDefSurfRegions_.GetSize() > 0 ) {
      ComputeVolDefSurf( volAboveDefSurfRegions_, volAboveDefSurfDir_ );
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

    StdVector<UInt> connecth;
    Vector<double> help;

    Double totalE = 0;

    UInt i, j;
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
      analysisVal = solveStep_->GetActFreq();
    }
    else {
      analysis    = "Time:";
      analysisVal = solveStep_->GetActTime();
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
