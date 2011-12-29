// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <stddef.h>
#include <complex>
#include <fstream>
#include <map>
#include <string>
#include <utility>

#include "CoupledPDE/pdecoupling.hh"
#include "DataInOut/Logging/cfslog.hh"
#include "DataInOut/Logging/log.hpp"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/WriteInfo.hh"
#include "Domain/ansatzFct.hh"
#include "Domain/elem.hh"
#include "Domain/entityList.hh"
#include "Domain/grid.hh"
#include "Domain/resultInfo.hh"
#include "Driver/assemble.hh"
#include "Driver/formsContext.hh"
#include "Driver/stdSolveStep.hh"
#include "Elements/basefe.hh"
#include "Forms/gradfieldop.hh"
#include "Forms/linBiotSavartSurfInt.hh"
#include "Forms/linGradBDBInt.hh"
#include "Forms/linearForm.hh"// for RHS
#include "General/Enum.hh"
#include "General/defs.hh"
#include "General/exception.hh"
#include "MatVec/exprt/xpr1.hh"
#include "MatVec/exprt/xpr2.hh"
#include "MatVec/matrix.hh"
#include "MatVec/vector.hh"
#include "Materials/baseMaterial.hh"
#include "PDE/SinglePDE.hh"
#include "PDE/StdPDE.hh"
#include "PDE/basePDE.hh"
#include "PDE/eqnMap.hh"
#include "Utils/StdVector.hh"
#include "Utils/basenodestoresol.hh"
#include "Utils/biotSavart.hh"
#include "Utils/nodestoresol.hh" // IWYU pragma: keep
#include "Utils/result.hh"
#include "magneticScalarPDE.hh"
#include "math.h"

namespace CoupledField {

class SingleVector;

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
    isAlwaysStatic_ = true;
    isMagnetostrictiveCoupled_ = false;
    
    

    // check if we have a Biot-Savart caculation activated or not
     isBiotSavart_ = false;
     
     // ToDO: Add excitation using Biot-Savart current sticks
     // ......
     
  }
  

  void MagScalarPDE::SetMagStrictCoupling() {
    isMagnetostrictiveCoupled_ = true;
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

      //==============================================================================
      //  LINEAR STIFFNESS INTEGRATOR
      //==============================================================================
      linGradBDBInt *  linMagForm = new linGradBDBInt( actSDMat, MAG_PERMEABILITY,
                                                       tensorType, upLagrangeForm );
      linMagForm->SetFactor( factor );

      BiLinFormContext * stiffIntDescr = new BiLinFormContext(linMagForm, STIFFNESS );

      stiffIntDescr->SetPtPdes(this, this);
      stiffIntDescr->SetResults( results_[0], results_[0],actSDList, actSDList );

      assemble_->AddBiLinearForm( stiffIntDescr );

      // Give result to equation numbering class 
      eqnMap_->AddResult( *results_[0], actSDList );

      //==============================================================================     
      // Additional section for Biot-Savart excitation
      //==============================================================================    
      
      // 1) Standard volume integrator
      if( biotSavart_ ) {
        BiotSavartSourceInt *rhsSource  = 
            new BiotSavartSourceInt( actSDMat, tensorType, 
                                     biotSavart_,false, upLagrangeForm); 
        // if the pde is magnetostrictive-coupled, the magnetic entries 
        // have to get multiplied with -1
        rhsSource->SetFactor(factor);        

        LinearFormContext * rhsContext = new LinearFormContext( rhsSource );
        rhsContext->SetPtPde( this );
        rhsContext->SetResult( results_[0], actSDList );
        assemble_->AddLinearForm( rhsContext );   
        // Give result to equation numbering class  //should it be here or above?
        eqnMap_->AddResult( *results_[0], actSDList );
      }     
    }
    //==============================================================================     
    // Additional section for Biot-Savart surface integrator
    //==============================================================================

    if(biotSavart_) {
      for( UInt i = 0; i < biotSavartBoundary_.GetSize(); i++ ) {

        // define boundary integrator
        LinBiotSavartSurfInt * bsInt = new LinBiotSavartSurfInt( biotSavart_, isaxi_);
        
        // set correct factor in magnetostrictive coupled case
        bsInt->SetFactor(factor); 
        bsInt->SetVoluInfo(materials_);
        LinearFormContext * surfContext = new LinearFormContext( bsInt );
        surfContext->SetResult(results_[0], biotSavartBoundary_[i]);
        surfContext->SetPtPde( this );
        
        assemble_->AddLinearForm(surfContext);
        eqnMap_->AddResult( *results_[0], biotSavartBoundary_[0]);
      } // loop over biotSavartBoundary
    } // if BiotSavart
  }

  void MagScalarPDE::DefineSolveStep() {
    solveStep_ = new StdSolveStep(*this);
  }
  
  void MagScalarPDE::PreparePDE4Computation() {
     if ( biotSavart_ ) {
       PtrParamNode bsNode = myParam_->Get("biotSavart", ParamNode::PASS );
       biotSavart_->Init( bsNode, ptgrid_, eqnMap_ );
     }
   }

  void MagScalarPDE::ReadSpecialBCs( )  {
    // -----------------------------------
    //  Check for Biot-Savart formulation
    // -----------------------------------
    PtrParamNode bsNode = myParam_->Get("biotSavart", ParamNode::PASS );
    if( bsNode ) {

      // check, if we have transient simulation
      if( analysistype_ == BasePDE::TRANSIENT ) {
        EXCEPTION( "Biot-Savart is currently just implemented for static /"
            << "harmonic simulations!" );
      }
      biotSavart_ = 
          shared_ptr<BiotSavart>(new BiotSavart());
      
      // Set formulation to H-field 
      biotSavart_->SetFormulation(BiotSavart::MAG_FIELD_STRENGTH);
    
      // Warning: do NOT set the flag "isBiotSavart_", because
      // this would trigger in the solveStep class the 
      // additional contribution of H to the solution vector, which
      // just makes sense for the vector potential formulation
    }
    
    
    // ***************************************************************
    //   Biot Savart boundary conditions
    // ***************************************************************
    PtrParamNode  bcNode = myParam_->Get( "bcsAndLoads", ParamNode::PASS );
    if( bcNode ) {
      ParamNodeList boundNodes = bcNode->GetList( "biotSavartBoundary" );

      for( UInt i = 0; i < boundNodes.GetSize(); i++ ) {
        std::string interfaceName = boundNodes[i]->Get("name")->As<std::string>();
        // create new surface stress definition
        biotSavartBoundary_.Push_back( ptgrid_->GetEntityList( EntityList::SURF_ELEM_LIST,interfaceName, EntityList::REGION ) );

        // add surface info definition
        Info->PrintF( pdename_,"Apply continuity Boundary Conditions on surfaceRegion '%s'\n",
                      interfaceName.c_str() );

      }
    }
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
        case MAG_ENERGY:
          if( isComplex_ ) {
            CalcEnergy<Complex>( res );
          } else {
            CalcEnergy<Double>( res );
          }
          break;
          
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
    tempH.Init();
    Hbiot.Resize(dim_);
    Hbiot.Init();
    NodeStoreSol<TYPE> & solhelp = dynamic_cast<NodeStoreSol<TYPE>&>(*sol_);
    
    GradientFieldOp<TYPE> * FieldOp = 
        new GradientFieldOp<TYPE>(ptgrid_,this,eqnMap_, solhelp, 
                                  results_[0]->fctType,false);
    Result<TYPE> &  actSol = dynamic_cast<Result<TYPE>&>(*sol);
    EntityIterator it = actSol.GetEntityList()->GetIterator();
    
    Vector<TYPE> & actVal = actSol.GetVector();
    actVal.Resize( actSol.GetEntityList()->GetSize() * dim_ );
    
    // loop over elements
    for ( it.Begin(); !it.IsEnd(); it++ ) {      
      Hbiot.Init();
      it.GetElem()->ptElem->GetCoordMidPoint(lCoord);
      
      ptgrid_->GetElemNodesCoord( CornerCoords, it.GetElem()->connect,true );
      it.GetElem()->ptElem->Local2GlobalCoord( globCoord, lCoord, CornerCoords,it.GetElem());      
      FieldOp->CalcElemGradField( tempH, it, lCoord, 1);       

      // if Biot-Savart coils are present, add contribution
      if( biotSavart_ ) {
        biotSavart_->CalcFieldAtPoint(Hbiot, globCoord );
      }

      // loop over dofs
      for(UInt iDim = 0; iDim < dim_; iDim++ ) {
        actVal[it.GetPos()*dim_ + iDim] = tempH[iDim]+ Hbiot[iDim];
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
  
 

  template <class TYPE>
  void MagScalarPDE::CalcEnergy( shared_ptr<BaseResult> vals )
  {

    Matrix<Double> permeability;  
    SubTensorType tensorType;
    if ( dim_ == 3 ) {
      tensorType = FULL;
    }    else {
      if ( isaxi_ == true ) {
        tensorType = AXI;
      }
      else {
        // 2d: plane case
        tensorType = PLANE_STRAIN;
      }
    }   
    Result<TYPE> &  actSol = dynamic_cast<Result<TYPE>&>(*vals);
    EntityIterator regionIt = actSol.GetEntityList()->GetIterator();

    // resize vector
    Vector<TYPE> & finalVal = actSol.GetVector();
    finalVal.Resize( actSol.GetEntityList()->GetSize() );

    //Energy
    TYPE energy = 0.0;

    // Loop over regions
    for( regionIt.Begin(); !regionIt.IsEnd(); regionIt++ ) {

      ElemList actSDList(ptgrid_ );
      actSDList.SetRegion( regionIt.GetRegion() );
      EntityIterator elemIt = actSDList.GetIterator();

      energy = 0.0;      
      BaseMaterial * ptMat = NULL;
      Matrix<Double> CornerCoords;

      NodeStoreSol<TYPE> & solhelp = dynamic_cast<NodeStoreSol<TYPE>&>(*sol_);
      GradientFieldOp<TYPE> * FieldOptmp = 
          new GradientFieldOp<TYPE>(ptgrid_,this,eqnMap_, solhelp, 
          results_[0]->fctType, results_[0], true );

      // loop over elements
      for ( elemIt.Begin(); !elemIt.IsEnd(); elemIt++ ) {

        TYPE energytmp=0.0; 
        ptgrid_->GetElemNodesCoord( CornerCoords, elemIt.GetElem()->connect,true );
        elemIt.GetElem()->ptElem->SetAnsatzFct( results_[0]->fctType );//first set the ansatz function, if it is higher order
        shared_ptr<AnsatzFct> ansatzFct1 = elemIt.GetElem()->ptElem->GetAnsatzFct();

        const Vector<Double> & intWeights =  elemIt.GetElem()->ptElem->GetIntWeights();
        Matrix<Double> permeability;
        
        ptMat = materials_[elemIt.GetElem()->regionId];
        ptMat->GetTensor( permeability, MAG_PERMEABILITY, Global::REAL, tensorType );
        const UInt nrIntPts = elemIt.GetElem()->ptElem->GetNumIntPoints();
        const Vector<Double> * intPoints = elemIt.GetElem()->ptElem->GetIntPoints();
        
        Vector<Double> CoordAtIP, biotH;
        Vector<TYPE> tempH, tempB;
        TYPE tmp1;
        for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++) {     
          Double jacDet = 0;
          jacDet=elemIt.GetElem()->ptElem->CalcJacobianDetAtIp(actIntPt,CornerCoords,elemIt.GetElem());
          elemIt.GetElem()->ptElem->Local2GlobalCoord( CoordAtIP, intPoints[actIntPt-1],CornerCoords, elemIt.GetElem() );
          FieldOptmp->CalcElemGradField( tempH, elemIt,intPoints[actIntPt-1], 1);//\xi
          if(biotSavart_){  
            biotSavart_->CalcFieldAtPoint(biotH,CoordAtIP);
            tempH = tempH + biotH;
          }
          tempB =permeability*tempH; 
          tmp1 = tempB * tempH;
          energytmp += tmp1*jacDet*0.5*intWeights[actIntPt-1];            

        } // loop over integration points          
        energy+=energytmp ;
      } // loop over elements
      
      finalVal[regionIt.GetPos()] = energy;
      
      // delete field operator
      delete FieldOptmp;
      
    } // loop over regions
}

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

      // Magnetic energy
      shared_ptr<ResultInfo> energy( new ResultInfo );
      energy->resultType = MAG_ENERGY;
      energy->definedOn = ResultInfo::REGION;
      energy->entryType = ResultInfo::SCALAR;
      energy->dofNames = "";
      energy->unit = "Ws";
      availResults_.insert ( energy );
    

  }
}
