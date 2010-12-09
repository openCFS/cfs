// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <fstream>
#include <iostream>
#include <sstream>
#include <math.h>
#include <string>

#include "magneticScalarPDE.hh"

#include "Forms/linGradBDBInt.hh"
#include "Forms/linNeumannInt.hh"
#include "Forms/mixedInt.hh"//added for surf form
#include "Forms/gradfieldop.hh"
#include "General/defs.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ParamHandling/ParamTools.hh"
#include "DataInOut/Logging/cfslog.hh"
#include "Utils/StdVector.hh"
#include "Driver/stdSolveStep.hh"
#include "CoupledPDE/pdecoupling.hh"
#include "Domain/ansatzFct.hh"
#include "Driver/assemble.hh"

#include "Forms/linearForm.hh"// for RHS

#ifdef USE_SCRIPTING
#include "DataInOut/Scripting/cfsmessenger.hh" 
#endif

namespace CoupledField {

  DECLARE_LOG(magnetopde)
  DEFINE_LOG(magnetopde, "pde.magneto")

  // ***************
  //   Constructor
  // ***************
  MagScalarPDE::MagScalarPDE( Grid* aptgrid, PtrParamNode paramNode )
    :SinglePDE( aptgrid, paramNode ) {


    // =====================================================================
    // set solution information
    // =====================================================================
    pdename_          = "magneticScalar";
    pdematerialclass_ = ELECTROMAGNETIC;
    maxTimeDerivOrder_ = 0;
 
    nonLin_    = false;
    nonLinMaterial_ = false;
    isAlwaysStatic_ = true;
    isMagnetostrictiveCoupled_ = false;
    
    

    // check if we have a Biot-Savart caculation activated or not
     isBiotSavart_ = false;
     
     // ToDO: Add excitation using Biot-Savart current sticks
     // ......
     
  }
  

  void MagScalarPDE::InitNonLin() {

    // ----------------------------------------------------------------------
    // At the moment, we do not consider any non-linear effects for this PDE.
    // In this section we would have to add e.g. the reading in and 
    // initialization of nonlinear regions, e.g. with hysteretic effects.
    // ----------------------------------------------------------------------
    
    
//    // Check, if "nonLinList" is present
//    ParamNode * nonLinListNode = myParam_->Get("nonLinList", false );
//    if( nonLinListNode ) { 
//
//      // Get nonlinear types
//      StdVector<ParamNode*> nonLinNodes = nonLinListNode->GetChildren();
//      for( UInt i = 0; i < nonLinNodes.GetSize(); i++ ) {
//
//        std::string actTypeString = nonLinNodes[i]->GetName();
//        std::string actId = nonLinNodes[i]->Get("id")->AsString();
//
//        NonLinType actType;
//        String2Enum( actTypeString, actType );
//        nonLinIdType_[actId] = actType;
//      }
//    }
//    
//    // Run over all region and set entry in "regionNonLinId"
//    StdVector<ParamNode*> regionNodes = 
//      myParam_->Get("regionList")->GetChildren();
//
//    RegionIdType actRegionId;
//    std::string actRegionName, actNonLinId;
//
//    if( regionNodes.GetSize() > 0 ) {
//      Info->PrintF( pdename_, "Non-linearity in following region(s)\n" );
//    }
//    for( UInt i = 0; i < regionNodes.GetSize(); i++ ) {
//      
//      regionNodes[i]->Get( "name", actRegionName );
//      regionNodes[i]->Get( "nonLinId", actNonLinId );
//      
//      if( actNonLinId == "" )
//        continue;
//
//      actRegionId = ptgrid_->RegionNameToId( actRegionName );
//
//      // Check nonLinId was already registerd
//      if( nonLinIdType_.find( actNonLinId) == nonLinIdType_.end() ) {
//        EXCEPTION( "NonLinearity with id '" << actNonLinId 
//                   << "' was not defined in 'nonLinList'" );
//      }
//      
//      regionNonLinId_[actRegionId] = actNonLinId;
//
//      // get related type of nonlinearity
//      NonLinType actType = nonLinIdType_[actNonLinId];
//      regionNonLinType_[actRegionId] = actType;
//
//      // check type
//      if( actType == HYSTERESIS ) {
//        isHysteresis_ = true;
//      }
//
//      if( actType == MATERIAL ) {
//        nonLin_ = true;
//        nonLinMaterial_ = true;
//      }
//
//      // Log to info file
//      std::string nonLinString;
//      Enum2String( nonLinIdType_[actNonLinId], nonLinString );
//      Info->PrintF( pdename_, " %s: %s\n", actRegionName.c_str(), 
//                    nonLinString.c_str() );
//      
//    }
//    Info->PrintF( pdename_, "\n" );

  }

  void MagScalarPDE::DefineIntegrators() {

    
    RegionIdType actRegion;
    BaseMaterial * actSDMat = NULL;
    
    // flag for updatedLagrange formulation
    bool upLagrangeForm = true;
 
    //transform the type
    SubTensorType tensorType;
    if ( dim_ == 3 ) {
      tensorType = FULL;
    }
    else {
      if ( isaxi_ == true ) {
        tensorType = AXI;
      }
      else {
        // 2d: plane case
        tensorType = PLANE_STRAIN;
      }
    }

    // if the pde is magnetostrictive-coupled, the magnetic entries
    // have to multiplied with -1
    std::string factor = "1.0";
    if ( isMagnetostrictiveCoupled_ == true )
      factor = "-1.0";  

    // Define integrators for "standard" materials
    std::map<RegionIdType, BaseMaterial*>::iterator it;
    for ( it = materials_.begin(); it != materials_.end(); it++ ) {
      
      // Set current region and material
      actRegion = it->first;
      actSDMat = it->second;

      // Get current region node
      std::string regionName = ptgrid_->GetRegion().ToString( actRegion );

      // create new entity list
      shared_ptr<ElemList> actSDList( new ElemList(ptgrid_ ) );
      actSDList->SetRegion( actRegion );

      // --- standard real-valued stiffness integrator ---
      linGradBDBInt *  linMagForm = new linGradBDBInt( actSDMat, MAG_PERMEABILITY,tensorType,upLagrangeForm);
      linMagForm->SetFactor( factor );

      BiLinFormContext * stiffIntDescr = new BiLinFormContext(linMagForm, STIFFNESS );

      stiffIntDescr->SetPtPdes(this, this);
      stiffIntDescr->SetResults( results_[0], results_[0],actSDList, actSDList );

      assemble_->AddBiLinearForm( stiffIntDescr );
      LOG_DBG2(magnetopde) <<"bilinearform add result 1";

      // Give result to equation numbering class 
      eqnMap_->AddResult( *results_[0], actSDList );
          

 
          //==============================================================================     
          // ToDO: Additional section for Biot-Savart excitation
          //==============================================================================    
//           if(isBiotSavart_){
//              LOG_DBG2(magnetopde) <<"bisa_.dim_: "<<bisa_.dim_;//called by regions             
////               LinearForm * rhsSource;
//              //should think about the 3d case
//               MagSourceInt *rhsSource = new MagSourceInt(actSDMat,bisa_,false, upLagrangeForm); 
//               // if the pde is magnetostrictive-coupled, the magnetic entries have to multiplied with -1
//               rhsSource->SetFactor(factor);
//               LinearFormContext * rhsContext = new LinearFormContext( rhsSource );
//               rhsContext->SetPtPde( this );
//               rhsContext->SetResult( results_[0], actSDList );
//               assemble_->AddLinearForm( rhsContext );   
//               // Give result to equation numbering class  //should it be here or above?
//              eqnMap_->AddResult( *results_[0], actSDList );
//           }     
              
           
//          if(isMagnetostrictiveCoupled_){          
////             LinearForm * rhsextra = new linMagnetostrictionInt(actSDMat,bisa_,false, upLagrangeForm);  
//             LinearForm * rhsextra = new MagnetoCouplingRHSInt(actSDMat,bisa_,false, upLagrangeForm); 
//             LinearFormContext * rhsContextSpecial = new LinearFormContext( rhsextra );
//             rhsContextSpecial->SetPtPde( this );
//             rhsContextSpecial->SetResult( results_[0], actSDList );
//             assemble_->AddLinearForm( rhsContextSpecial );   
//             // Give result to equation numbering class  //should it be here or above?
//            eqnMap_->AddResult( *results_[0], actSDList );
//       
//          }
//          
    }

//    // **********************************************************************
//    //  Surf bilinear Form
//    // **********************************************************************        
//      
//        
//        std::map<RegionIdType,ListInfo>::iterator  BCsIt;
//        UInt actSD=0;
//        for( BCsIt = surfinfo_.begin();BCsIt != surfinfo_.end(); BCsIt++ ) {
//         
////        for (UInt actSD = 0; actSD < conBCs_.GetSize(); actSD++) {
//
//           shared_ptr<SurfElemList> surfElems( new SurfElemList(ptgrid_ ) );             
////           shared_ptr<ElemList> volElems( new ElemList(ptgrid_ ) );
//
//           surfElems->SetRegion(BCsIt->second.surface);
////           volElems->SetRegion(actRegion );
//           
//           ContinuityBCsInt * bilinear_surf = new ContinuityBCsInt(dim_,surfElems,isaxi_);
//           bilinear_surf->SetMaterial(actSDMat);
//           std::cout<<"act regionId:"<<actRegion<<std::endl;
//           BiLinFormContext * surfContext = new BiLinFormContext( bilinear_surf, STIFFNESS,true );//for eqns assembling
//           std::cout<<"act regionId:"<<BCsIt->second.region<<std::endl;
//           surfContext->SetVolumeNeighbor(BCsIt->second.region );
//           std::cout<<"bilinearform add result 2"<<std::endl;
//           surfContext->SetPtPdes(this, this); 
//           std::cout<<"bilinearform add result 3"<<std::endl;
//           surfContext->SetResults( results_[0], results_[0],conBCs_[actSD], conBCs_[actSD] );
//           std::cout<<"bilinearform add result 4"<<std::endl;
//           assemble_->AddBiLinearForm( surfContext );
//           std::cout<<"bilinearform add result 5"<<std::endl;
//           // Give result to equation numbering class
//           eqnMap_->AddResult( *results_[0], conBCs_[actSD] );
//           actSD++;
//           
//         }
//    
//    // **********************************************************************
//             //  Surf RHS
//             // **********************************************************************
//             for( UInt iBc = 0; iBc < inBcs_.GetSize(); iBc++ ) {
//
//               // get current Bc
//               InhomNeumannBc const & actBc = *inBcs_[iBc];
//
//               ContinuityRHS *ContinuityBC = new ContinuityRHS(isaxi_ ,bisa_);
//               ContinuityBC->SetMaterial(actSDMat);
//               // if the pde is magnetostrictive-coupled, the magnetic entries have to multiplied with -1
//               ContinuityBC->SetFactor(factor);
//   //          ContinuityBC->SetVoluInfo( materials_ );
//               LinearFormContext * BcContext =new LinearFormContext( ContinuityBC );
//               BcContext->SetPtPde( this );
//   //--------------Reference----------------
//   //            pressRhs->SetResult( results_[0], pressSurf_[actSF] );
//   //           assemble_->AddLinearForm( pressRhs );
//   //           // Give entities and result to equation numbering class
//   //           // and solution class
//   //           eqnMap_->AddResult( *results_[0], pressSurf_[actSF] );
//   //            
//               BcContext->SetResult( actBc.result, actBc.entities );
//               assemble_->AddLinearForm( BcContext );
//
//               // Give result to equation numbering class
//               eqnMap_->AddResult( *actBc.result, actBc.entities );
//             }
//              //==============================================================================  
//
//    
  
  }
  

  void MagScalarPDE::DefineSolveStep() {
    solveStep_ = new StdSolveStep(*this);
  }

  void MagScalarPDE::ReadSpecialBCs( ) 
  {
//    // ***************************************************************
//    //   continuity boundary conditions
//    // ***************************************************************
//
//    continuityBC_ = false;
//    ParamNode * bcNode = myParam_->Get( "bcsAndLoads", false );
//    if( bcNode ) {
//      StdVector<ParamNode*> contiBCNodes = bcNode->GetList( "neumannInhom" );
//
//      for( UInt i = 0; i < contiBCNodes.GetSize(); i++ ) {
//        std::string interfaceName = contiBCNodes[i]->Get("name")->AsString();
//        std::string regionName = contiBCNodes[i]->Get("region")->AsString();
//        // create new surface stress definition
//        ListInfo actinfo;
//        actinfo.surface = ptgrid_->RegionNameToId(interfaceName);
//        actinfo.region = ptgrid_->RegionNameToId(regionName);
//        conBCs_.Push_back( ptgrid_->GetEntityList( EntityList::SURF_ELEM_LIST,interfaceName, EntityList::REGION ) );
//        continuityBC_= true;
//        
//        // add surface info definition
//        RegionIdType regionId =  ptgrid_->RegionNameToId( interfaceName );
//        surfinfo_[regionId] = actinfo;
//        
//        Info->PrintF( pdename_,"Apply continuity Boundary Conditions on surfaceRegion '%s'\n",interfaceName.c_str() );
//      
//      }
//    }
//    LOG_DBG2(magnetopde) << "MagScalarPDE::ReadSpecialBCs()" ;
//   
//    LOG_DBG2(magnetopde)<<"conBCs_.GetSize():"<<conBCs_.GetSize();
  }
  
  
 
  // ======================================================
  // POSTPROCESSING SECTION
  // ======================================================

  // ***************
  //   CalcResults
  // ***************
   void MagScalarPDE::CalcResults( shared_ptr<BaseResult> res ) {
     
     switch (res->GetResultInfo()->resultType ) {
     
       case MAG_SCALAR_POTENTIAL:
          if( isComplex_ ) {
            ExtractResult<Complex>( res, sol_ );
          } else {
            ExtractResult<Double>( res, sol_ );
          }
          break;
        
        case MAG_RHS_LOAD:
          if( isComplex_ ) {
//            std::cout<<"isComplex"<<std::endl;
            ExtractRhsResult<Complex>( res, results_[0] );
          } else {
            LOG_DBG2(magnetopde) <<"notComplex";
            ExtractRhsResult<Double>( res, results_[0] );
          }
          break;

        case MAG_FLUX_DENSITY: //B field
          if( isComplex_ ) {
            CalcMagneticFluxDensity<Complex>( res );
          } else {
            LOG_DBG2(magnetopde) <<"notComplex--CalcMagneticFluxDensity";
            CalcMagneticFluxDensity<Double>( res );
          }
          break;
        case MAG_HFIELD: //H field
          if( isComplex_ ) {
            CalcMagneticField<Complex>( res );
          } else {
            CalcMagneticField<Double>( res );
          }
          break;   
//        
//        case MAG_HS_FIELD: //H_s field
//           if( isComplex_ ) {
//             CalcHs<Complex>( res );
//           } else {
//             CalcHs<Double>( res );
//           }
//         break;   
//         
//        case MAG_HM_FIELD: //H_m field
//            if( isComplex_ ) {
//             CalcHm<Complex>( res );
//            } else {
//             CalcHm<Double>( res );
//            }
//        break;  
//          
//        case MAG_ENERGY:
//          if( isComplex_ ) {
//            CalcEnergy<Complex>( res );
//          } else {
//            CalcEnergy<Double>( res );
//          }
//          break;
          
        default:
          WARN("Resulttype not computable by magneticscalar PDE");
        }        
     
   }


  template <class TYPE>//---H field-----------------------------------------
  void MagScalarPDE:: CalcMagneticField( shared_ptr<BaseResult> sol )
  {     
    Matrix<Double> CornerCoords;
    Vector<Double> globCoord,Hbiot,lCoord;
    Vector<TYPE> tempH;
    tempH.Resize(dim_);
    Hbiot.Resize(dim_);
    NodeStoreSol<TYPE> & solhelp = dynamic_cast<NodeStoreSol<TYPE>&>(*sol_);
//    std::cerr<< "sol_ is " << solVec_ << std::endl;
    GradientFieldOp<TYPE> * FieldOp = 
        new GradientFieldOp<TYPE>(ptgrid_,this,eqnMap_, solhelp, 
                                  results_[0]->fctType,false);
    Result<TYPE> &  actSol = dynamic_cast<Result<TYPE>&>(*sol);
    EntityIterator it = actSol.GetEntityList()->GetIterator();
    
    Vector<TYPE> & actVal = actSol.GetVector();
    actVal.Resize( actSol.GetEntityList()->GetSize() * dim_ );
    
    // loop over elements
    for ( it.Begin(); !it.IsEnd(); it++ ) {      
      it.GetElem()->ptElem->GetCoordMidPoint(lCoord);
      
      LOG_DBG2(magnetopde) <<"---Debug Part---lCoord---:"<<lCoord;
      ptgrid_->GetElemNodesCoord( CornerCoords, it.GetElem()->connect,true );
      it.GetElem()->ptElem->Local2GlobalCoord( globCoord, lCoord, CornerCoords,it.GetElem());      
      LOG_DBG2(magnetopde) <<"---Debug Part---globCoord---:"<<globCoord;             
      FieldOp->CalcElemGradField( tempH, it, lCoord, 1);       
      LOG_DBG2(magnetopde) <<"tempH :"<<tempH<<"---end--";
      
      if(isBiotSavart_){
      WARN("Biot Savart is not yet implemented");
      }
      // loop over dofs
      for(UInt iDim = 0; iDim < dim_; iDim++ ) {
        if(isBiotSavart_) {          
//          actVal[it.GetPos()*dim_ + iDim] = Hbiot[iDim];//only biot savart source case
          actVal[it.GetPos()*dim_ + iDim] = tempH[iDim]+ Hbiot[iDim];
          
        }
        else actVal[it.GetPos()*dim_ + iDim] = tempH[iDim];
      }
    }
    delete FieldOp;
    
  }
  
//  template <class TYPE>//---H_s field-----------------------------------------
//  void MagScalarPDE::CalcHs( shared_ptr<BaseResult> sol )
//  {     
//    Matrix<Double> CornerCoords;
//    Vector<Double> globCoord,Hbiot,lCoord;
//    Result<TYPE> &  actSol = dynamic_cast<Result<TYPE>&>(*sol);
//    EntityIterator it = actSol.GetEntityList()->GetIterator();
//    
//    Vector<TYPE> & actVal = actSol.GetVector();
//    actVal.Resize( actSol.GetEntityList()->GetSize() * dim_ );
//   if(isBiotSavart_){ 
//      // loop over elements
//      for ( it.Begin(); !it.IsEnd(); it++ ) {      
//        it.GetElem()->ptElem->GetCoordMidPoint(lCoord);
//        
//        LOG_DBG2(magnetopde) <<"---Debug Part---lCoord---:"<<lCoord;
//        ptgrid_->GetElemNodesCoord( CornerCoords, it.GetElem()->connect,true );
//        it.GetElem()->ptElem->Local2GlobalCoord( globCoord, lCoord, CornerCoords,it.GetElem());      
//        LOG_DBG2(magnetopde) <<"---Debug Part---globCoord---:"<<globCoord;
//        bisa_.CalcH(Hbiot,globCoord);
//        
//        // loop over dofs
//        for(UInt iDim = 0; iDim < dim_; iDim++ ) {
//            actVal[it.GetPos()*dim_ + iDim] = Hbiot[iDim];//only biot savart source case
//        }
//      }
//   }else
//     EXCEPTION("Biot Savart is not activated, could not cacluate the magnetic field generated by biot savart");
//}
//  
//  
//  template <class TYPE>//---H_m field-----------------------------------------
//  void MagScalarPDE::CalcHm( shared_ptr<BaseResult> sol )
//  {     
//
//    Vector<Double> lCoord;
//    Vector<TYPE> tempH;
//    tempH.Resize(dim_);
//    NodeStoreSol<TYPE> & solhelp = dynamic_cast<NodeStoreSol<TYPE>&>(*sol_);
//    GradientFieldOp<TYPE> * FieldOp = new GradientFieldOp<TYPE>(ptgrid_,this,eqnMap_, solhelp, 
//                                MAG_SCALAR_POTENTIAL, results_[0],false);
//    Result<TYPE> &  actSol = dynamic_cast<Result<TYPE>&>(*sol);
//    EntityIterator it = actSol.GetEntityList()->GetIterator();
//    
//    Vector<TYPE> & actVal = actSol.GetVector();
//    actVal.Resize( actSol.GetEntityList()->GetSize() * dim_ );
//    
//    // loop over elements
//    for ( it.Begin(); !it.IsEnd(); it++ ) {      
//      it.GetElem()->ptElem->GetCoordMidPoint(lCoord);
//      
//      FieldOp->CalcElemGradField( tempH, it, lCoord, 1); 
//      
//      LOG_DBG2(magnetopde) <<"tempH :"<<tempH<<"---end--";
//      
//      // loop over dofs
//      for(UInt iDim = 0; iDim < dim_; iDim++ ) {
//         actVal[it.GetPos()*dim_ + iDim] = tempH[iDim];
//      }
//    }
//    delete FieldOp;
//    
//  }
  
  
  template <class TYPE>//---B field----------------------------------------
  void MagScalarPDE::CalcMagneticFluxDensity( shared_ptr<BaseResult> sol )
  {
    
//     Double mu0=1.25664e-06;// vacuum
    
    // first calculate magnetic field intensity (H)
      CalcMagneticField<TYPE>( sol );
      LOG_DBG2(magnetopde)<<"CalcMagneticField ends ";
    // then, run over all elements and multiply by element permeability
      Result<TYPE> &  actSol = dynamic_cast<Result<TYPE>&>(*sol);
      EntityIterator it = actSol.GetEntityList()->GetIterator();

      Vector<TYPE> & actVal = actSol.GetVector();

      Matrix<Double> permeability;
      BaseMaterial * ptMat = NULL;
      Vector<TYPE> tempH (dim_), tempB (dim_);
//      Vector<TYPE> tmp;
      SubTensorType tensorType;
      if( dim_ == 2 ) {
        tensorType = PLANE;
      } else {
        tensorType = FULL;
      }
      
      // loop over elements
      for ( it.Begin(); !it.IsEnd(); it++ ) {
            
        ptMat = materials_[it.GetElem()->regionId];
        ptMat->GetTensor( permeability, MAG_PERMEABILITY, Global::REAL, tensorType );
                
        LOG_DBG2(magnetopde) <<"permeability: "<<permeability;
        // build temporary element vector
        tempH.Init();
        for(UInt iDim = 0; iDim < dim_; iDim++ ) {
          tempH[iDim] = actVal[it.GetPos()*dim_ + iDim];
        }
        
//        tempB = tempH *mu0; //only biot savart source field

        tempB = permeability*tempH;// all together
        
        // loop over dofs
        for(UInt iDim = 0; iDim < dim_; iDim++ ) {
          actVal[it.GetPos()*dim_ + iDim] = tempB[iDim];
        }
      }

  }
  
 

//  template <class TYPE>
//  void MagScalarPDE::CalcEnergy( shared_ptr<BaseResult> vals )
//  {
//
//    Matrix<Double> permeability;  
//    SubTensorType tensorType;
//       
//    if ( dim_ == 3 ) {
//      tensorType = FULL;
//    }
//    else {
//      if ( isaxi_ == true ) {
//        tensorType = AXI;
//      }
//      else {
//        // 2d: plane case
//        tensorType = PLANE_STRAIN;
//      }
//    }   
//
//    Result<TYPE> &  actSol = dynamic_cast<Result<TYPE>&>(*vals);
//    EntityIterator regionIt = actSol.GetEntityList()->GetIterator();
//   
//    // resize vector
//    Vector<TYPE> & finalVal = actSol.GetVector();
//    finalVal.Resize( actSol.GetEntityList()->GetSize() );// vector according to the regions
//      
//    //Energy
//    Double energy;
//    
//    // Loop over regions
//      for( regionIt.Begin(); !regionIt.IsEnd(); regionIt++ ) {
//        
//      
//        ElemList actSDList(ptgrid_ );
//        actSDList.SetRegion( regionIt.GetRegion() );
//        EntityIterator elemIt = actSDList.GetIterator();
//        
//        energy = 0;//each region has a energy      
//    
//        Vector<Double> scalarpot;
//        BaseMaterial * ptMat = NULL;
//        Vector<Double> tempH;
//        tempH.Resize(dim_);
//        Matrix<Double> CornerCoords;
//        
//        NodeStoreSol<Double> & solhelp = dynamic_cast<NodeStoreSol<Double>&>(*sol_);
//        GradientFieldOp<Double> * FieldOptmp = new GradientFieldOp<Double>(ptgrid_,this,eqnMap_, solhelp, 
//                                        MAG_SCALAR_POTENTIAL, results_[0],false);
//       
//        
//      
//       // loop over elements
//       for ( elemIt.Begin(); !elemIt.IsEnd(); elemIt++ ) {
//
//         Double energytmp=0.0; 
//                 
//          ptgrid_->GetElemNodesCoord( CornerCoords, elemIt.GetElem()->connect,true );
//          elemIt.GetElem()->ptElem->SetAnsatzFct( results_[0]->fctType );//first set the ansatz function, if it is higher order
//          shared_ptr<AnsatzFct> ansatzFct1 = elemIt.GetElem()->ptElem->GetAnsatzFct();
//          const Vector<Double> & intWeights =  elemIt.GetElem()->ptElem->GetIntWeights();
//          //=======================
//          Double permeability;
//          ptMat = materials_[elemIt.GetElem()->regionId];
//          ptMat->GetScalar(permeability,MAG_PERMEABILITY,Global::REAL);
//          LOG_DBG2(magnetopde) <<"current scalar permeability : "<<permeability;
////          std::cout<<"current scalar permeability : "<<permeability<<std::endl;
//          //=======================          
//         
//          const UInt nrIntPts = elemIt.GetElem()->ptElem->GetNumIntPoints();
//          UInt numFncs = elemIt.GetElem()->ptElem->GetNumFncs( ansatzFct1 );
//          LOG_DBG2(magnetopde) <<"numFncs: "<<numFncs<<std::endl;
//          LOG_DBG2(magnetopde) <<"Vector Potential Contribution ----4------ ";
//          const Vector<Double> * intPoints = elemIt.GetElem()->ptElem->GetIntPoints();
//          
//          Vector<Double> ShpFncAtIp, helpVec, CoordAtIP,Hfield,H_field;
//          Double tmp1;
//          LOG_DBG2(magnetopde) <<"number of intergration points:"<<nrIntPts;
//          for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++) {     
//      
//            LOG_DBG2(magnetopde)<<"actual intergration point:"<<actIntPt;
//            Double jacDet = 0;//intialize the Jacobian determinate to be zero
//            jacDet=elemIt.GetElem()->ptElem->CalcJacobianDetAtIp(actIntPt,CornerCoords,elemIt.GetElem());
//            // get the Jacobian
//            LOG_DBG2(magnetopde) <<"CornerCoords :"<<CornerCoords;
//            elemIt.GetElem()->ptElem->Local2GlobalCoord( CoordAtIP, intPoints[actIntPt-1],CornerCoords, elemIt.GetElem() );
//            elemIt.GetElem()->ptElem->GetShFncAtIp(ShpFncAtIp,actIntPt,elemIt.GetElem(),bisa_.dim_);
//            LOG_DBG2(magnetopde) <<"ShpFncAtIp:"<<ShpFncAtIp;
////======================================================================================================              
//            FieldOptmp->CalcElemGradField( tempH, elemIt,intPoints[actIntPt-1], 1);//\xi
////======================================================================================================              
//            if(isBiotSavart_){  
//              bisa_.CalcH(Hfield,CoordAtIP);//have to use the global coord
//              H_field.Resize(bisa_.dim_);
//              for(UInt iDim = 0; iDim < bisa_.dim_; iDim++ ) {
//                   H_field[iDim]= Hfield[iDim]+tempH[iDim];
//              }
//              LOG_DBG2(magnetopde) <<"Hfield : "<<Hfield; 
//            }else{
//              H_field.Resize(dim_);
//              for(UInt iDim = 0; iDim < bisa_.dim_; iDim++ ) {
//                     H_field[iDim]= tempH[iDim];
//               }
//            }
//              
//              tmp1 = H_field*H_field;
//              energytmp += tmp1*jacDet*0.5*permeability*intWeights[actIntPt-1];            
//              
//       }          
//          energy+=energytmp ;
//     }
//       finalVal[regionIt.GetPos()] = energy;
//
//   }
      
//}

  // ======================================================
  // COUPLING SECTION
  // ======================================================

  void MagScalarPDE::InitCoupling(PDECoupling * Coupling)
  {
  
    isIterCoupled_ = true;
    ptCoupling_   = Coupling;
    
    nonLin_ = false;


  }
  


  void MagScalarPDE::CalcOutputCoupling()
  {

    SolutionType quantity;
    StdVector<UInt> * couplingNodes     = NULL;
    SingleVector * values = 0;
//    UInt forcesCount = 0;

    // at first, check if this PDE is iterative coupled
    if (isIterCoupled_ == false)
      return;

    // loop over all output coupling quantities
    for (UInt actCoupling=0; 
         actCoupling<ptCoupling_->GetNumOutputCouplings(); 
         actCoupling++)
      {
        quantity = ptCoupling_->GetOutputQuantity(actCoupling);
        ptCoupling_->GetOutputValues(actCoupling, values);

//        Vector<Double> * temp = dynamic_cast<Vector<Double> *>(values);
      
        switch(ptCoupling_->GetOutputType(actCoupling))
          {
          
          case NODE:        
            ptCoupling_->GetOutputNodes(actCoupling, couplingNodes);
          
            if (quantity == MAG_SCALAR_POTENTIAL)
              sol_->NodeSolutionToCoupling(*values, *couplingNodes);
                  
            break;
          
          case ELEM:
            EXCEPTION("No Element coupling output");
          }
      }
  }


  bool MagScalarPDE::HasOutput(SolutionType output)
  {
  
    switch (output)
      {
      case MAG_SCALAR_POTENTIAL:
        return true;
        break;
      case MAG_RHS_LOAD:
        return true;
        break;
      default:
        return false;
        break;
      }
    return false;
  }


//  void MagScalarPDE::SetMagnetostrictionCoupling(Bisa & pointer)
//  {
//  
//    isMagnetostrictiveCoupled_= true;
//    if(isBiotSavart_)    pointer = bisa_;
//    else{
//      LOG_DBG2(magnetopde) <<"Do nothing to the pointer";
//    }
//  }
  
  
  void MagScalarPDE::DefineAvailResults() {
    
    // Determine vectorial dofNames
    StdVector<std::string> vecDofNames;
    if( isaxi_) {
      vecDofNames = "r", "z";
    } else if (dim_ == 2) {
      vecDofNames = "x", "y";
    } else {
      vecDofNames = "x", "y", "z";
    }
    
    //=== Magnetic scalar Potential ===
    // check if problem is lagrange or legendre
    std::string approxType = myParam_->Get("type")->As<std::string>();
    LOG_DBG2(magnetopde) <<"approxType:"<<approxType;
    
    shared_ptr<ResultInfo> res1(new ResultInfo);
    res1->resultType = MAG_SCALAR_POTENTIAL;
    res1->dofNames = "";
    res1->unit = "V";
    res1->entryType = ResultInfo::SCALAR;
    
    if ( approxType == "lagrange" ) {
      shared_ptr<AnsatzFct> fct(new LagrangeFct);
      res1->definedOn = ResultInfo::NODE;
      res1->fctType = fct;
    } else if(approxType == "legendre"){
      // Create new resultDof object
      shared_ptr<LegendreFct> fct(new LegendreFct);
      if( myParam_->Get("isIsotropic")->As<bool>() ) {
        LOG_DBG2(magnetopde) <<"isIsotropic: true";
        UInt order =  myParam_->Get("order")->As<UInt>();
        fct->SetIsoOrder( order );
      } 
      res1->definedOn = ResultInfo::PFEM;
      res1->fctType = fct;
    } else {
      EXCEPTION( "Approximation type '" << approxType
                  << "' not known" );
    }
    
    results_.Push_back( res1 );
    availResults_.insert( res1 ); 
    
    
 
        
      //=== Magnetic RHS Load ===
      // create new resultDof object
      shared_ptr<ResultInfo> rhs ( new ResultInfo );
      rhs->resultType = MAG_RHS_LOAD;
      rhs->dofNames = "";
      rhs->unit = "Am";
      rhs->definedOn = results_[0]->definedOn;//==ResultInfo::NODE
      rhs->entryType = ResultInfo::SCALAR;
      availResults_.insert( rhs );
  
      // Magnetic Field Intensity
      // create new resultDof object
      shared_ptr<ResultInfo> hfield ( new ResultInfo );
      hfield->resultType = MAG_HFIELD;
      hfield->dofNames = vecDofNames;
      hfield->unit = "A/m";
      hfield->definedOn = ResultInfo::ELEMENT;
      hfield->entryType = ResultInfo::VECTOR;
      availResults_.insert( hfield );
      
//      shared_ptr<ResultInfo> hsfield ( new ResultInfo );
//      hsfield->resultType = MAG_HS_FIELD;
//      hsfield->dofNames = vecDofNames;
//      hsfield->unit = "A/m";
//      hsfield->definedOn = ResultInfo::ELEMENT;
//      hsfield->entryType = ResultInfo::VECTOR;
//      availResults_.insert( hsfield );
//      
//      shared_ptr<ResultInfo> hmfield ( new ResultInfo );
//      hmfield->resultType = MAG_HM_FIELD;
//      hmfield->dofNames = vecDofNames;
//      hmfield->unit = "A/m";
//      hmfield->definedOn = ResultInfo::ELEMENT;
//      hmfield->entryType = ResultInfo::VECTOR;
//      availResults_.insert( hmfield );
      
      // Magnetic Flux Density
      shared_ptr<ResultInfo> flux ( new ResultInfo );
      flux->resultType = MAG_FLUX_DENSITY;//B Field
      flux->dofNames = vecDofNames;
      flux->unit = "Vs/m^2";
      flux->definedOn = ResultInfo::ELEMENT;
      flux->entryType = ResultInfo::VECTOR;
      availResults_.insert( flux );
//     
//      // Magnetic energy
//      shared_ptr<ResultInfo> energy( new ResultInfo );
//      energy->resultType = MAG_ENERGY;
//      energy->definedOn = ResultInfo::REGION;
//      energy->entryType = ResultInfo::SCALAR;
//      energy->dofNames = "";
//      energy->unit = "Ws";
//      availResults_.insert ( energy );
    

  }
}
