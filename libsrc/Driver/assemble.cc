#include <fstream>
#include <iostream>
#include <string>
#include <math.h>

#include "Domain/elem.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "DataInOut/WriteInfo.hh"
#include "PDE/StdPDE.hh"
#include "PDE/SinglePDE.hh"
#include "Domain/domain.hh"
#include "Utils/coordSystem.hh"


#include "olas.hh"

#include "assemble.hh"
#include "MHassemble.hh"

namespace CoupledField {
  
  Assemble::Assemble(BaseSystem * algsys, Grid * aptgrid)
    :algsys_(algsys),
     ptgrid_(aptgrid),
     integrators_(0),
     rhsIntegrators_(0),
     rhsSrcIntegrators_(0),
     pdeId1_(NO_PDE_ID),
     pdeId2_(NO_PDE_ID)
  {
    ENTER_FCN( "Assemble::Assemble" );
    firstTime_ = TRUE;
    oneIntIsNonlin_ = FALSE;
    nrMatrices_ = 5+1;  // 5 matrices, but index starts with 1!
    reassembleMat_.Resize(nrMatrices_);
    // reassembleMat_.resize(nrMatrices_);
    nonLinGeo = FALSE;
    deltaCoords_ = NULL;
    alternateMaterialData_=FALSE;

    olasParams_ = algsys_->GetOLASParams();
    olasReport_ = algsys_->GetOLASReport(); 

    startFreq_ = 0.0;
    matArray_  = NULL;

  }


  // **************
  //   Destructor
  // **************
  Assemble::~Assemble() {

    ENTER_FCN( "Assemble::~Assemble" );
    
    for ( UInt i = 0; i < integrators_.GetSize(); i++ ) {
      for ( UInt j = 0; j < integrators_[i]->GetSize(); j++ ) {
        delete (*integrators_[i])[j];
      }
      delete integrators_[i];
    }
    integrators_.Clear();

    for ( UInt i = 0; i < surfintegrators_.GetSize(); i++ ) {
      for ( UInt j = 0; j < surfintegrators_[i]->GetSize(); j++ ) {
        delete (*surfintegrators_[i])[j];
      }
      delete surfintegrators_[i];
    }
    surfintegrators_.Clear();

    for ( UInt i = 0; i < rhsIntegrators_.GetSize(); i++ ) {
      for ( UInt j = 0; j < rhsIntegrators_[i]->GetSize(); j++ ) {
        delete (*rhsIntegrators_[i])[j];
      }
      delete rhsIntegrators_[i];
    }
    rhsIntegrators_.Clear();

    for ( UInt i = 0; i < rhsSrcIntegrators_.GetSize(); i++ ) {
      for ( UInt j = 0; j < rhsSrcIntegrators_[i]->GetSize(); j++ ) {
        delete (*rhsSrcIntegrators_[i])[j];
      }
      delete rhsSrcIntegrators_[i];
    }
    rhsSrcIntegrators_.Clear();

    for ( UInt i = 0; i < rhsSrcSurfIntegrators_.GetSize(); i++ ) {
      for ( UInt j = 0; j < rhsSrcSurfIntegrators_[i]->GetSize(); j++ ) {
        delete (*rhsSrcSurfIntegrators_[i])[j];
      }
      delete rhsSrcSurfIntegrators_[i];
    }
    rhsSrcIntegrators_.Clear();
  }


  void Assemble::SetPDEId( const PdeIdType id1,
                           const PdeIdType id2 ) {
    ENTER_FCN( "Assemble::SetPdeId" );
    
    if ( id1 == NO_PDE_ID ) {
      Error( "Assemble::SetPdeId: The Id is empty!", __FILE__, __LINE__ );
    }

    pdeId1_ = id1;

    // Set identifier for second pde equal to that of the
    // first one. This indicates, that the element matrices
    // have to be assembled onto the on-diagonal positions.
    if ( id2 == NO_PDE_ID)
      pdeId2_ = pdeId1_;
    else
      pdeId2_ = id2;
  }
  
  void Assemble::SetPtr2EQNData(NodeEQN * aPtNodeEQN1,
                                NodeEQN * aPtNodeEQN2)
  {
    ENTER_FCN( "Assemble::SetPtr2EQNData" );

    if (aPtNodeEQN1 == NULL)
      Error( "The NodeEQN object is not created yet" , 
             __FILE__, __LINE__);
    
    if (! aPtNodeEQN1->IsInitialized())
      Error( "The NodeEQN object has to be initialized before assigning it to\
              an assemble object by calling 'CalcMapping()'", 
             __FILE__, __LINE__);
    ptEQN1_ = aPtNodeEQN1;

    // check if a second eqn-object is present.
    // If not, the same is used two times
    if (aPtNodeEQN2 == NULL)
      ptEQN2_ = aPtNodeEQN1;
    else
      ptEQN2_ = aPtNodeEQN2;
    
  }


  UInt Assemble::SubDomIndex(const RegionIdType subDomId)
  {
    ENTER_FCN( "Assemble::SubDomIndex" );

    for (UInt i=0; i < subdoms_.GetSize();i++)
      {
        if (subDomId == subdoms_[i]) return i;
      }
  
    (*error) <<  "SubDomain " <<  ptgrid_->RegionIdToName(subDomId) 
             << " not defined!";
    Error( __FILE__, __LINE__ );

    return 0;
  }


  // ****************
  //   SurfDomIndex
  // ****************
  UInt Assemble::SurfDomIndex( const RegionIdType surfDomId ) {

    ENTER_FCN( "Assemble::SurfDomIndex" );

    for ( UInt i = 0; i < surfdoms_.GetSize(); i++ ) {
      if ( surfDomId == surfdoms_[i] ) {
        return i;
      }
    }

    (*error) << "Surface-Domain " << ptgrid_->RegionIdToName( surfDomId )
             << " not defined!";
    Error( __FILE__, __LINE__ );
    return 0;
  }


  void Assemble::GetElemCoords(const StdVector<UInt> connect, 
                               Matrix<Double> &coordMat)
  {
    ENTER_FCN( "Assemble:GetElemCoords" );
    
    ptgrid_->GetElemNodesCoord( coordMat, connect);
  
    if (nonLinGeo)
      {
        if (deltaCoords_ == NULL) {
          (*error) << " ElecPDE: set input_coupling_terms = "
                   << " smoothdisplacement or nonlin = no";
          Error( __FILE__, __LINE__ );
        }
        StdVector<UInt> connect_PDE;
        ptEQN1_->Mesh2PDENode(connect_PDE, connect);
        Double val;
        for (UInt i=0; i<coordMat.GetSizeRow(); i++)
          for (UInt j=0; j<coordMat.GetSizeCol(); j++) 
            {
              val = (*deltaCoords_)[i][connect_PDE[j] - 1];
              coordMat(i,j) += val;
            }
      }
  }
  
  
  // do the basic assembling stuff
  void Assemble::AssembleMatrices()
  {
    ENTER_FCN( "Assemble:AssembleMatrices" );
   
    SETPROFILE("Before AssembleMatrices");

    Vector<Double> harmonicVec;
    Double dampTransform = 1.0;

    std::string analysis;
    params->Get( "type", analysis, "analysis" );

    // initialize reassembling "indicator" vector
    if (firstTime_)
      for (UInt actMat=0; actMat < nrMatrices_; actMat++)
        reassembleMat_[actMat] = FALSE;

    for (UInt actDom=0; actDom < subdoms_.GetSize(); actDom++) {     
      StdVector<Elem*> elemssd;

      ptgrid_->GetElems(elemssd, subdoms_[actDom]);

      for(UInt actInteg=0; actInteg < integrators_[actDom]->GetSize(); 
          actInteg++) {

        IntegratorDescriptor * actDescriptor = 
          (*integrators_[actDom])[actInteg];

        // get PDE-Ids
        PdeIdType pdeId1 = actDescriptor->GetPDE1()->GetPDEId();
        PdeIdType pdeId2 = actDescriptor->GetPDE2()->GetPDEId();

        // get equation data of PDEs
        NodeEQN * ptEQN1 = actDescriptor->GetPDE1()->getPDE_eqnData();
        NodeEQN * ptEQN2 = actDescriptor->GetPDE2()->getPDE_eqnData();
                    
        actDescriptor->GetIntegrator()->SetSubdomain(actDom);

        if (alternateMaterialData_ == TRUE)
          actDescriptor->GetIntegrator()->SetMaterial(ptMaterial_);

        // assemble only if nonlinear or first time
        if (reassembleMat_[actDescriptor->DestMat()] || firstTime_) {

          if ( actDescriptor->IsReducedInt()) {
            //set all elements to reduced integration!!
            SetFE2ReducedInt();
          }

          if ( actDescriptor->GetIntegrator()->IsFracDamping() ) {
            // get multiplicative pre factor depending on time step size
            Double damp;
            damp = ptPDE_->GetFracDampMatrixCoeff(actDom);
            actDescriptor->GetIntegrator()->SetFactor(damp);
          }

          dampTransform = 1.0;
          if ( ( actDescriptor->GetIntegrator()->IsRaylDamping() 
                 || actDescriptor->GetSecondaryMat() != NOTYPE ) 
               && startFreq_ > 0 ) {

            // Obtain frequency value to which the damping parameters
            // in the material file do belong
            StdVector<Double> freqs;
            Double matDataFreq;
            params->GetList( "matDataFreq", freqs, "harmonic" );
            if ( freqs.GetSize() == 1 ) {
              matDataFreq = freqs[0] * 2.0 * PI;
            }
            else {
              matDataFreq = startFreq_;
            }

            // get multiplicative pre factor depending on frequency
            if ( matDataFreq > 0 && actFreq_ > 0 ) {
              FEMatrixType destMat =
                actDescriptor->GetIntegrator()->GetBaseType();

              if ( destMat == STIFFNESS ) {
                dampTransform = matDataFreq / actFreq_;
                if(analysis != "paramIdent")
                  Info->PrintF( "", " dampTransform for STIFFNESS ... %e\n",
                                dampTransform );
              }
              else if ( destMat == MASS ) {
                dampTransform = actFreq_ / matDataFreq;
                if (analysis != "paramIdent")
                  Info->PrintF( "", " dampTransform for MASS ........ %e\n",
                                dampTransform );
              }
            }

            if ( actDescriptor->GetIntegrator()->IsRaylDamping() ) {
              actDescriptor->GetIntegrator()->SetFactor(dampTransform);
            }
          }

	  // pass frequency value to bilinearform
	  if (analysisType_ == HARMONIC) {
	    actDescriptor->GetIntegrator()->SetFrequency(actFreq_);
	  }

	  //put pointer to array containing the material parameter 
          // for each element
          actDescriptor->GetIntegrator()->SetMaterialArray(matArray_);

          for (UInt actEl=0; actEl< elemssd.GetSize(); actEl++) {
            actDescriptor->GetIntegrator()->SetElemNr(actEl);

            BaseFE * ptEl = elemssd[actEl]->ptElem;
            StdVector<UInt> connecth = elemssd[actEl]->connect;
                             
            Matrix<Double> ptCoord;
            GetElemCoords(connecth, ptCoord);

            // map connect to PDE node numbers
            StdVector<Integer> connect_PDE1, connect_PDE2;
                    
            ptEQN1->Node2EQN(connecth, connect_PDE1);
            ptEQN2->Node2EQN(connecth, connect_PDE2);                  

            Matrix<Double> elSol;
                    
            actDescriptor->GetIntegrator()->SetElemPtr(ptEl);
            FEMatrixType destMat = actDescriptor->DestMat();
                
            // this matrix is nonlinear and, therefore, 
            // has to be reassembled next time
            if (actDescriptor->IsNonLin()) {
              oneIntIsNonlin_ = TRUE;
              reassembleMat_[actDescriptor->DestMat()] = TRUE;
              sol_->GetElemSolutionAsMatrix(elSol, connecth);
              actDescriptor->GetIntegrator()->SetActElemSol(elSol);
            }       

            // ================================================================
            //                             assemble matrices
            // ================================================================

            actDescriptor->GetIntegrator()->
              CalcElementMatrix(ptCoord, elemmat);

            piezoMaterialType matType = actDescriptor->GetPiezoMaterialType();
            actDescriptor->SetPiezoMaterialType(matType);

            if (analysisType_ == HARMONIC) {
              TransformMatrix2Harmonic(harmonicVec,elemmat, 
                                       actDescriptor->GetOrigMatrixType(),
                                       actDescriptor->GetPiezoMaterialType());

              algsys_->SetElementMatrix( destMat, &harmonicVec[0], 
                                         pdeId1, connect_PDE1.GetPointer(), 
                                         connect_PDE1.GetSize(),
                                         pdeId2, connect_PDE2.GetPointer(), 
                                         connect_PDE2.GetSize(),
                                         actDescriptor->IsSetCounterPart());
            }
            else {
              algsys_->SetElementMatrix( destMat, elemmat.GetDataPointer(), 
                                         pdeId1, connect_PDE1.GetPointer(), 
                                         connect_PDE1.GetSize(), 
                                         pdeId2, connect_PDE2.GetPointer(), 
                                         connect_PDE2.GetSize(),
                                         actDescriptor->IsSetCounterPart());

            }
#ifdef DEBUG
            UInt elemNum = elemssd[actEl]->elemNum;

            // output matrices
            if (destMat == STIFFNESS) {
              (*debug) << "Stiffness matrix of Element " 
                       << elemNum << std::endl;
            }

            if (destMat == MASS) {
              (*debug) << "Mass      matrix of Element " 
                       << elemNum << std::endl;
            }

            if (destMat == DAMPING) {
              (*debug) << "Damping   matrix of Element " 
                       << elemNum << std::endl;
            }

            if (destMat == SYSTEM) {
              (*debug) << "System matrix of Element " 
                       << elemNum << std::endl;
            }
            
            (*debug) << elemmat << std::endl;

            if ( !elemmat.IsSymmetric() ) {
              (*debug) << " --> Matrix is not symmetric " 
                       << std::endl << std::endl;
            }
            else {
              (*debug) << " --> Matrix is symmetric " 
                       << std::endl << std::endl;
            }
#endif
            if (actDescriptor->GetSecondaryMat() != NOTYPE) {

              Double damp = dampTransform * actDescriptor->GetSecMatFac();
              elemmat *= damp;
              if (analysisType_ == HARMONIC) {
                TransformMatrix2Harmonic(harmonicVec,elemmat,
                                         actDescriptor->GetOrigSecMatrixType(),
                                         actDescriptor->GetPiezoMaterialType());

                algsys_->SetElementMatrix(destMat, &harmonicVec[0], 
                                          pdeId1, connect_PDE1.GetPointer(),
                                          connect_PDE1.GetSize(), 
                                          pdeId2, connect_PDE2.GetPointer(),
                                          connect_PDE2.GetSize(),
                                          actDescriptor->IsSetCounterPart());

              }
              else {
                algsys_->SetElementMatrix(actDescriptor->GetSecondaryMat(), 
                                          elemmat.GetDataPointer(), 
                                          pdeId1, connect_PDE1.GetPointer(), 
                                          connect_PDE1.GetSize(), 
                                          pdeId2, connect_PDE2.GetPointer(), 
                                          connect_PDE2.GetSize(),
                                          actDescriptor->IsSetCounterPart());
              }
            }

          } //over all elements of subdomain            
                
        } //check, if we have to assemble
          
        if ( actDescriptor->IsReducedInt()) {
          //set all elements back to standard integration!!
          SetFE2StandardInt();
        }
          
      } //integrators

    } //subdomains

      //assemble matrices for surface integrators
    for (UInt actDom=0; actDom < surfdoms_.GetSize(); actDom++) {    

      ptgrid_->RegionIdToName(surfdoms_[actDom]);
      
      //check is necessary, because the surface-integrator 
      // could also be a RHS-Src integrator
      if (surfintegrators_[actDom]->GetSize()) {        
        StdVector<SurfElem*> elemssd;
        ptgrid_->GetSurfElems(elemssd, surfdoms_[actDom]);
          
        for (UInt actEl=0; actEl< elemssd.GetSize(); actEl++) {
          BaseFE * ptEl = elemssd[actEl]->ptElem;
          StdVector<UInt> connecth = elemssd[actEl]->connect;

          Matrix<Double> ptCoord;
          GetElemCoords(connecth, ptCoord);        
              
          Matrix<Double> elSol;
          Vector<Double> normal;
              
          // this matrix is nonlinear and, therefore, 
          // has to be reassembled next time

          if (oneIntIsNonlin_ ) {
            // fetch solution at element nodes
            //GetSolOfElement(elSol, connect_PDE);
            sol_->GetElemSolutionAsMatrix(elSol, connecth);
          }
                
          // ================================================================
          //                             assemble matrices
          // ================================================================
                
          for(UInt actInteg=0; actInteg < surfintegrators_[actDom]->GetSize(); 
              actInteg++) {

            IntegratorDescriptor * actDescriptor = 
              (*surfintegrators_[actDom])[actInteg];

            SurfForm * myForm = 
              dynamic_cast<SurfForm*>(actDescriptor->GetIntegrator());


            // get PDE-Ids
            PdeIdType pdeId1 = actDescriptor->GetPDE1()->GetPDEId();
            PdeIdType pdeId2 = actDescriptor->GetPDE2()->GetPDEId();

            // get equation data of PDEs
            NodeEQN * ptEQN1 = actDescriptor->GetPDE1()->getPDE_eqnData();
            NodeEQN * ptEQN2 = actDescriptor->GetPDE2()->getPDE_eqnData();

            // map connect to PDE node numbers
            StdVector<Integer> connect_PDE1, connect_PDE2;
            ptEQN1->Node2EQN(connecth, connect_PDE1);
            ptEQN2->Node2EQN(connecth, connect_PDE2);

            // assemble only if nonlinear or first time
            if (reassembleMat_[actDescriptor->DestMat()] || firstTime_) {
          
              myForm->SetElemPtr(ptEl);
              myForm->SetSurfElem(elemssd[actEl]);
          
              // calculate normal of surface element
              ptgrid_->CalcSurfNormal(normal,*elemssd[actEl]);
              normal *= (Double) elemssd[actEl]->normalSign;
              myForm->SetFirstVoluNormal(normal);

              // calculate normal

              FEMatrixType destMat = actDescriptor->DestMat();

              
              // this matrix is nonlinear and, therefore, 
              // has to be reassembled next time
              if (actDescriptor->IsNonLin()) {
                oneIntIsNonlin_ = TRUE;
                reassembleMat_[actDescriptor->DestMat()] = TRUE;
                actDescriptor->GetIntegrator()->SetActElemSol(elSol);
              }
                      
              actDescriptor->GetIntegrator()->CalcElementMatrix(ptCoord, elemmat);

              if (analysisType_ == HARMONIC) {
                TransformMatrix2Harmonic(harmonicVec,elemmat, 
                                         actDescriptor->GetOrigMatrixType(),
                                         actDescriptor->GetPiezoMaterialType());

                algsys_->SetElementMatrix( destMat, &harmonicVec[0], 
                                           pdeId1, connect_PDE1.GetPointer(), 
                                           connect_PDE1.GetSize(),
                                           pdeId2, connect_PDE2.GetPointer(), 
                                           connect_PDE2.GetSize(),
                                           actDescriptor->IsSetCounterPart());
              }
              else {
                algsys_->SetElementMatrix( destMat, elemmat.GetDataPointer(), 
                                           pdeId1, connect_PDE1.GetPointer(), 
                                           connect_PDE1.GetSize(),
                                           pdeId2, connect_PDE2.GetPointer(), 
                                           connect_PDE2.GetSize(),
                                           actDescriptor->IsSetCounterPart());
              }

                        
              if (actDescriptor->GetSecondaryMat()  != NOTYPE ) {
                elemmat *= actDescriptor->GetSecMatFac();
                if (analysisType_ == HARMONIC) {
                  TransformMatrix2Harmonic(harmonicVec,elemmat,
                                           actDescriptor->GetOrigSecMatrixType(),
                                           actDescriptor->GetPiezoMaterialType());
                  
                  algsys_->SetElementMatrix( destMat, &harmonicVec[0], 
                                             pdeId1, connect_PDE1.GetPointer(), 
                                             connect_PDE1.GetSize(),
                                             pdeId2, connect_PDE2.GetPointer(), 
                                             connect_PDE2.GetSize(),
                                             actDescriptor->IsSetCounterPart());
                }
                else{
                  algsys_->SetElementMatrix(  actDescriptor->GetSecondaryMat(), 
                                              elemmat.GetDataPointer(), 
                                              pdeId1, connect_PDE1.GetPointer(), 
                                              connect_PDE1.GetSize(),
                                              pdeId2, connect_PDE2.GetPointer(), 
                                              connect_PDE2.GetSize(),
                                              actDescriptor->IsSetCounterPart());
                }
              }
            }           
          }
        } // elements
      } // check for surface integrator
    } // subdomains
     
    firstTime_ = FALSE;


    SETPROFILE("After AssembleMatrices");
  }


  // do the basic assembling stuff
  void Assemble::AssembleSrcRHS(const Double time)
  {
    ENTER_FCN( "Assemble:AssembleSrcRHS" );
    AssembleRHSNodalSources(time);
    AssembleRHSIntegralSources(time);
  }


  // do the basic assembling stuff
  void Assemble::AssembleRHSIntegralSources(const Double time)
  {
    ENTER_FCN( "Assemble:AssembleRHSIntegralSources" );

    Vector<Double> harmVec;

    //add the volume sources
    for (UInt actDom=0; actDom <  subdoms_.GetSize(); actDom++) { 
      if (rhsSrcIntegrators_[actDom]->GetSize())
        {
          StdVector<Elem*> elemssd;
          ptgrid_->GetElems(elemssd, subdoms_[actDom]);
            
          Double val_tfunc = 1.0;
          Double valPhase = 0.0;
          if (analysisType_ == HARMONIC||analysisType_ == MULTIHARMONIC) 
            valPhase = rhsSrcPhase_[actDom];
          else {
            if (ptTimeFunc_->GetmaxTimeFnc() > 0 )
              val_tfunc=ptTimeFunc_->TimeFuncAtTime(time,fncname_rhs_[actDom]);
          }

          for (UInt actEl=0; actEl< elemssd.GetSize(); actEl++) {        
            
            BaseFE * ptEl = elemssd[actEl]->ptElem;
            StdVector<UInt> connecth = elemssd[actEl]->connect;
                
            Matrix<Double> ptCoord;
            GetElemCoords(connecth, ptCoord);
            
            // map connect to PDE node numbers
            StdVector<Integer> connect_PDE;
            ptEQN1_->Node2EQN(connecth, connect_PDE);
            for(UInt actRhsInt=0; 
                actRhsInt < rhsSrcIntegrators_[actDom]->GetSize(); 
                actRhsInt++)
              {
                BaseIntDescriptor * actRhsID = 
                  (*rhsSrcIntegrators_[actDom])[actRhsInt];
                    
                actRhsID->GetIntegrator()->SetElemPtr(ptEl);
                    
                Vector<Double> elemVec;
                actRhsID->GetIntegrator()->CalcElemVector(ptCoord, elemVec);
                    
                if (analysisType_ == HARMONIC||analysisType_ == MULTIHARMONIC) {
                  TransformVector2Harmonic(harmVec,elemVec,valPhase);
                  algsys_->SetElementRHS(&harmVec[0], pdeId1_, 
                                         connect_PDE.GetPointer(), 
                                         connect_PDE.GetSize());
                }
                else {
                  if (val_tfunc != 1.0)
                    elemVec *= val_tfunc;
                  algsys_->SetElementRHS(&elemVec[0], pdeId1_, 
                                         connect_PDE.GetPointer(), 
                                         connect_PDE.GetSize());
                }
              }
          }
        }
    }

    //add the surface sources
    for (UInt actSurf=0; actSurf <  surfdoms_.GetSize(); actSurf++) { 

      if (rhsSrcSurfIntegrators_[actSurf]->GetSize()) {
        StdVector<SurfElem*> elemssd;
        ptgrid_->GetSurfElems(elemssd, surfdoms_[actSurf]);
            
        Double val_tfunc = 1.0;
        Double valPhase = 0.0;
        if (analysisType_ == HARMONIC||analysisType_ == MULTIHARMONIC) 
          valPhase = rhsSrcSurfPhase_[actSurf];
        else {
          if (ptTimeFunc_->GetmaxTimeFnc() > 0 )
            val_tfunc=ptTimeFunc_->TimeFuncAtTime(time,fncname_rhsSurf_[actSurf]);
        }

        for (UInt actEl=0; actEl< elemssd.GetSize(); actEl++) {        
          BaseFE * ptEl = elemssd[actEl]->ptElem;
          StdVector<UInt> connecth = elemssd[actEl]->connect;
                
          Matrix<Double> ptCoord;
          GetElemCoords(connecth, ptCoord);
            
          // map connect to PDE node numbers
          StdVector<Integer> connect_PDE;
          ptEQN1_->Node2EQN(connecth, connect_PDE);
          for(UInt actRhsInt=0; 
              actRhsInt < rhsSrcSurfIntegrators_[actSurf]->GetSize(); 
              actRhsInt++) {
                  
            BaseIntDescriptor * actRhsID = 
              (*rhsSrcSurfIntegrators_[actSurf])[actRhsInt];
             
            LinearSurfForm * myForm =  dynamic_cast<LinearSurfForm*>(actRhsID->GetIntegrator());
             
            myForm->SetElemPtr(ptEl);
            myForm->SetSurfElem(elemssd[actEl]);

            // calculate normal of surface element
            Vector<Double> normal;
            ptgrid_->CalcSurfNormal(normal,*elemssd[actEl]);
            normal *= (Double) elemssd[actEl]->normalSign;
            myForm->SetVoluNormal(normal);
              
            Vector<Double> elemVec;
            myForm->CalcElemVector(ptCoord, elemVec);
            if (analysisType_ == HARMONIC||analysisType_ == MULTIHARMONIC) {
              TransformVector2Harmonic(harmVec,elemVec,valPhase);
              algsys_->SetElementRHS(&harmVec[0], pdeId1_, connect_PDE.GetPointer(), 
                                     connect_PDE.GetSize());
            }
            else {
              if (val_tfunc != 1.0)
                elemVec *= val_tfunc;
              algsys_->SetElementRHS(&elemVec[0], pdeId1_, 
                                     connect_PDE.GetPointer(), 
                                     connect_PDE.GetSize());
            }
          }
        }
      }
    }

  }



  // do the basic assembling stuff
  void Assemble::AssembleRHSNodalSources(const Double time)
  {
    ENTER_FCN( "Assemble:AssembleRHSNodalSources" );
    
    Integer eqnNr = 0;
    UInt eqnDof = 0;
    Double phase = 0;
    
    for (UInt actDom=0; actDom < loadDom_.GetSize(); actDom++) {
      std::string doftype = loadDom_[actDom];

      UInt dof = 1;
      if (dofsPerNode_ != 1)
        dof = domain->GetCoordSystem()->GetVecComponent( loadDof_[actDom] );

      StdVector<UInt> nodes;
      ptgrid_->GetNodesByName( nodes,loadDom_[actDom]);
      Double val = loadVals_[actDom];

      Double val_tfunc = 1.0;
      if (ptTimeFunc_->GetmaxTimeFnc() > 0 )
        val_tfunc=ptTimeFunc_->TimeFuncAtTime(time,fncname_loads_[actDom]);


      for ( UInt i=0;  i<nodes.GetSize(); i++) {
        UInt node = nodes[i];
            
        val = loadVals_[actDom] * val_tfunc;

        ptEQN1_->Node2EQN(node,dof,eqnNr,eqnDof);

        if (analysisType_ == HARMONIC || analysisType_ == MULTIHARMONIC) {
          phase = loadPhase_[actDom];
          Complex complexValue( val * cos( phase / 180 * PI ),
                                val * sin( phase / 180 * PI ) );
          algsys_->SetNodeRHS(complexValue, pdeId1_, eqnNr, eqnDof);    
        } else {
          algsys_->SetNodeRHS(val, pdeId1_, eqnNr, eqnDof);    
        }
      }
    }
  }
  






  // do the basic assembling stuff
  void Assemble::AssembleNLRHS(const Double time)
  {
    ENTER_FCN( "Assemble:AssembleNLRHS" );

    Matrix<Double> elemmat;

    StdVector<Elem*> elemssd;
   
    for (UInt actDom=0; actDom < subdoms_.GetSize(); actDom++) { 
      if (rhsIntegrators_[actDom]->GetSize()) {
        ptgrid_->GetElems(elemssd, subdoms_[actDom]);
        
        for (UInt actEl=0; actEl< elemssd.GetSize(); actEl++) {
          BaseFE * ptEl = elemssd[actEl]->ptElem;
          StdVector<UInt> connecth = elemssd[actEl]->connect;


          Matrix<Double> ptCoord;
          GetElemCoords(connecth, ptCoord);


          // map connect to PDE node numbers
          StdVector<Integer> connect_PDE;
          ptEQN1_->Node2EQN(connecth, connect_PDE);
                
          Matrix<Double> elSol;
                
          sol_->GetElemSolutionAsMatrix(elSol, connecth);
                
          // ================================================================
          //                             assemble RHS
          // ================================================================
                
          for(UInt actRhsInt=0; 
              actRhsInt < rhsIntegrators_[actDom]->GetSize(); actRhsInt++) {

            BaseIntDescriptor * actRhsID = 
              (*rhsIntegrators_[actDom])[actRhsInt];
                    
            actRhsID->GetIntegrator()->SetElemPtr(ptEl);
                    
            if (actRhsID->IsNonLin())
              actRhsID->GetIntegrator()->SetActElemSol(elSol);
                    
                    
            Vector<Double> elemVec;
            actRhsID->GetIntegrator()->CalcElemVector(ptCoord, elemVec);
            algsys_->SetElementRHS(&elemVec[0], pdeId1_, 
                                   connect_PDE.GetPointer(), 
                                   connect_PDE.GetSize());
          }
        }
      }
    }
  }


  // *******************
  //   AssembleSprings
  // *******************
  void Assemble::AssembleSprings( const Double time ) {

    ENTER_FCN( "Assemble::AssembleSprings" );
    
    UInt eqnDof;
    Integer eqnNr;
    
    for (UInt actDom=0; actDom < springDom_.GetSize(); actDom++) {
      std::string doftype = springDom_[actDom];

      UInt dof = 1;
      if ( dofsPerNode_ != 1 ) {
        dof = domain->GetCoordSystem()->GetVecComponent( springDof_[actDom] );
      }

      StdVector<UInt> nodes;
      ptgrid_->GetNodesByName(nodes, springDom_[actDom]);
        
      Double massValue_ = springMassVals_[actDom];
      Double dampingValue_ = springDampVals_[actDom];
      Double stiffnessValue_ = springStiffVals_[actDom];

      Double val_tfunc = 1.0;
      if ( ptTimeFunc_->GetmaxTimeFnc() > 0 ) {
        val_tfunc = ptTimeFunc_->TimeFuncAtTime( time,
                                                 fncname_springs_[actDom] );
      }

      for ( UInt i = 0; i < nodes.GetSize(); i++ ) {

        UInt node = nodes[i];

        massValue_      = springMassVals_[actDom]  * val_tfunc;
        dampingValue_   = springDampVals_[actDom]  * val_tfunc;
        stiffnessValue_ = springStiffVals_[actDom] * val_tfunc;

        ptEQN1_->Node2EQN( node, dof, eqnNr, eqnDof );

        if ( analysisType_ == TRANSIENT ) {
          if ( abs(massValue_) > 1e-30 ) {
            Info->PrintF( "", "Adding value %e to the mass matrix\n",
                          massValue_ );
            algsys_->AddToDiagMatrixEntry( MASS, pdeId1_, eqnNr, eqnDof,
                                           &massValue_ );
          }
          if ( abs(dampingValue_) > 1e-30 ) {

            // if (!dampingMatrix_)
            // Error("The damping value of a spring can only be added to 
            // the damping matrix when there exist one! ",__FILE__,__LINE__);

            Info->PrintF( "", "Adding value %e to the damping matrix\n",
                          dampingValue_ );
            algsys_->AddToDiagMatrixEntry( DAMPING, pdeId1_, eqnNr, eqnDof,
                                           &dampingValue_ );
          }
          if ( abs(stiffnessValue_) > 1e-30 ) {
            Info->PrintF( "", "Adding value %e to the stiffness matrix\n",
                          stiffnessValue_ );
            algsys_->AddToDiagMatrixEntry( STIFFNESS, pdeId1_, eqnNr, eqnDof,
                                           &stiffnessValue_ );
          }
        }

        else if( analysisType_ == STATIC ) {
          if ( abs(stiffnessValue_) > 1e-30 ) {
            Info->PrintF( "", "Adding value %e to the system matrix\n",
                          stiffnessValue_ );
            algsys_->AddToDiagMatrixEntry( SYSTEM, pdeId1_, eqnNr, eqnDof,
                                           &stiffnessValue_ );
          }
          if ( abs(dampingValue_) > 1e-30 || abs(massValue_) > 1e-30 ) {
            (*error) << "The damping and mass value of a spring will not "
                     << "be considered in an static analysis!";
            Error( __FILE__, __LINE__ );
          }
        }
      }
    }
  }


  void Assemble::InitNonLinMatrices() {
    ENTER_FCN( "Assemble::InitNonLinMatrices" );
    
    // return, if matrices are not yet assembled
    if (!reassembleMat_.GetSize()) {
      // if ( reassembleMat_.GetSize() == 0 ) {
      algsys_->InitMatrix();
      return;
    }
    
    // Initialize matrices in order to get BCs correct
    algsys_->InitMatrix(SYSTEM);

    if ( reassembleMat_[STIFFNESS] ) 
      algsys_->InitMatrix(STIFFNESS);
    
    
    if ( reassembleMat_[DAMPING] ) 
      algsys_->InitMatrix(DAMPING);

    if ( reassembleMat_[CONVECTION] )
      algsys_->InitMatrix(CONVECTION);
    
    if ( reassembleMat_[MASS] ) 
      algsys_->InitMatrix(MASS); 
    
  }


  void Assemble::SetGeneralParams(const std::string & pdename, 
                                  const UInt dofsPerNode,
                                  const StdVector<RegionIdType> & subdoms,
                                  const StdVector<RegionIdType> & surfdoms,
                                  const std::string bcSequenceTag)
  {
    ENTER_FCN( "Assemble::SetGeneralParams" );

    pdename_       = pdename;
    dofsPerNode_   = dofsPerNode;
    subdoms_       = subdoms;
    surfdoms_      = surfdoms;
    bcSequenceTag_ = bcSequenceTag;
 


    // parameters for reading load values =========================================
    StdVector<std::string> keyVec, attrVec, valVec;

    // parameters for reading spring values =======================================
    StdVector<std::string> keyVecSpring;

    attrVec = "", "tag", "";
    valVec = "", bcSequenceTag_, "";

    loadDom_.Clear();
    loadDof_.Clear();
    loadVals_.Clear();
    fncname_loads_.Clear();

    keyVec = pdename_, "bcsAndLoads", "load", "name";
    params->GetList(keyVec, attrVec, valVec, loadDom_);

    keyVecSpring = pdename_, "bcsAndLoads", "spring", "name";
    params->GetList(keyVecSpring, attrVec, valVec, springDom_);

    if (dofsPerNode_ == 1) {
      loadDof_.Resize(loadDom_.GetSize());
      loadDof_.Init("ux");

      springDof_.Resize(springDom_.GetSize());
      springDof_.Init("ux");
    }else {
      keyVec = pdename_, "bcsAndLoads", "load", "dof";
      params->GetList(keyVec, attrVec, valVec, loadDof_);

      keyVecSpring = pdename_, "bcsAndLoads", "spring", "dof";
      params->GetList(keyVecSpring, attrVec, valVec, springDof_);

    }
    
    keyVec = pdename_, "bcsAndLoads", "load", "value";
    params->GetList(keyVec, attrVec, valVec, loadVals_);

    keyVec = pdename_, "bcsAndLoads", "load", "phase";
    params->GetList(keyVec, attrVec, valVec, loadPhase_);

    keyVecSpring = pdename_, "bcsAndLoads", "spring", "massValue";
    params->GetList(keyVecSpring, attrVec, valVec, springMassVals_);

    keyVecSpring = pdename_, "bcsAndLoads", "spring", "dampingValue";
    params->GetList(keyVecSpring, attrVec, valVec, springDampVals_);

    keyVecSpring = pdename_, "bcsAndLoads", "spring", "stiffnessValue";
    params->GetList(keyVecSpring, attrVec, valVec, springStiffVals_);

    keyVec = pdename_, "bcsAndLoads", "load", "dynamics";
    params->GetList(keyVec, attrVec, valVec, fncname_loads_);

    keyVecSpring = pdename_, "bcsAndLoads", "spring", "dynamics";
    params->GetList(keyVecSpring, attrVec, valVec, fncname_springs_);


    // check if a fncname for load was given, although
    // in section transient there is non mentioned
    if (fncname_loads_.GetSize() > 0) {
      if (ptTimeFunc_->GetmaxTimeFnc() == 0  &&
          fncname_loads_[0] != "none") {
        std::string errmsg = "Loads: There was no time data file ";
        errmsg += "specified in the section 'transient', so 'dynamics' ";
        errmsg += "for PDE '";
        errmsg += pdename;
        errmsg += "' makes no sense!";
        Error(errmsg.c_str(), __FILE__, __LINE__);
      }
    }

    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // this section must be changed so that it is possible to
    // specify two or more timeDataFile in section transient
    // and to check them here
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // check if a fncname for spring was given, although
    // in section transient there is non mentioned
    if (fncname_springs_.GetSize() > 0)
      {
        if (ptTimeFunc_->GetmaxTimeFnc() == 0  &&
            fncname_springs_[0] != "none") {
          std::string errmsg = "Springs: There was no time data file ";
          errmsg += "specified in the section 'transient', so 'dynamics' ";
          errmsg += "for PDE '";
          errmsg += pdename;
          errmsg += "' makes no sense!";
          Error(errmsg.c_str(), __FILE__, __LINE__);
        }
      }
    

    //     std::cerr << "load names = " <<  loadDom_ << std::endl;
    //     std::cerr << "load dofs = " <<  loadDof_ << std::endl;
    //     std::cerr << "load values = " <<  loadVals_ << std::endl;
    //     std::cerr << "load dynamics = " <<  fncname_loads_ << std::endl;

    // Check load consistency
    if ( loadDom_.GetSize() != loadDof_.GetSize() ||
         loadDom_.GetSize() != loadVals_.GetSize() ||  
         loadDom_.GetSize() != loadPhase_.GetSize() ) {
      std::string errmsg = "Loads: ";
      errmsg += "#name = " + Info->GenStr(loadDom_.GetSize());
      errmsg += ", #dof = " + Info->GenStr(loadDof_.GetSize());
      errmsg += ", #value = " + Info->GenStr(loadVals_.GetSize());
      errmsg += ", #phase = " + Info->GenStr(loadPhase_.GetSize());
      errmsg += ", #dynamics = " + fncname_loads_.GetSize() + '\n';
      Info->Error( errmsg, __FILE__, __LINE__ );
    }

    // Check spring consistency
    if ( springDom_.GetSize() != springDof_.GetSize() ||
         springDom_.GetSize() != springMassVals_.GetSize() ||
         springDom_.GetSize() != springDampVals_.GetSize() ||
         springDom_.GetSize() != springStiffVals_.GetSize() ) {
      std::string errmsg = "Springs: ";
      errmsg += "#name = " + Info->GenStr(springDom_.GetSize());
      errmsg += ", #dof = " + Info->GenStr(springDof_.GetSize());
      errmsg += ", #mass value = " + Info->GenStr(springMassVals_.GetSize());
      errmsg += ", #damping value = " + Info->GenStr(springDampVals_.GetSize());
      errmsg += ", #stiffness value = " + Info->GenStr(springStiffVals_.GetSize());
      errmsg += ", #dynamics = " + fncname_springs_.GetSize() + '\n';
      Info->Error( errmsg, __FILE__, __LINE__ );
    }

    // We need not have as many function/filenames as loads!
    for ( UInt k = fncname_loads_.GetSize(); k < loadDom_.GetSize(); k++ )
      {
        fncname_loads_.Push_back( "none" );
      }

    // ????????????????????????????????????????????????????????
    // We need not have as many function/filenames as springs!
    for ( UInt k = fncname_springs_.GetSize(); k < springDom_.GetSize(); k++ )
      {
        fncname_springs_.Push_back( "none" );
      }

#ifdef DEBUG
    // loads
    (*debug) << "Assemble::SetGeneralParams: We got " << loadDom_.GetSize()
             << " interfaces with loads" << std::endl;
    (*debug) << "Loads: #interfaces = " << loadDom_.GetSize()
             << ", #dof = " << loadDof_.GetSize()
             << ", #value = " << loadVals_.GetSize()
             << ", #phase = " << loadPhase_.GetSize()
             << ", #dynamics = " << fncname_loads_.GetSize() << std::endl;
    for ( UInt k = 0; k < loadDom_.GetSize(); k++ )
      {
        (*debug) << "Loads: interface = " << loadDom_[k]
                 << ", dof = " << loadDof_[k]
                 << ", value = " << loadVals_[k]
                 << ", phase = " << loadPhase_[k];
        if ( k < fncname_loads_.GetSize() )
          {
            (*debug) << ", dynamics = " << fncname_loads_[k] << std::endl;
          }
        else
          {
            (*debug) << ", dynamics = " << std::endl;
          }
      }
    // springs
    (*debug) << "Assemble::SetGeneralParams: We got " << springDom_.GetSize()
             << " interfaces with springs" << std::endl;
    (*debug) << "Springs: #interfaces = " << springDom_.GetSize()
             << ", #dof = " << springDof_.GetSize()
             << ", #mass value = " << springMassVals_.GetSize()
             << ", #damping value = " << springDampVals_.GetSize()
             << ", #stiffness value = " << springStiffVals_.GetSize()
             << ", #dynamics = " << fncname_springs_.GetSize() << std::endl;
    for ( UInt k = 0; k < springDom_.GetSize(); k++ )
      {
        (*debug) << "Springs: interface = " << springDom_[k]
                 << ", dof = " << springDof_[k]
                 << ", mass value = " << springMassVals_[k]
                 << ", damping value = " << springDampVals_[k]
                 << ", stiffness value = " << springStiffVals_[k];
        if ( k < fncname_springs_.GetSize() )
          {
            (*debug) << ", dynamics = " << fncname_springs_[k] << std::endl;
          }
        else
          {
            (*debug) << ", dynamics = " << std::endl;
          }
      }
#endif

    // for every domain, we need an own integrator list ==========
    integrators_.Resize(subdoms_.GetSize());
    rhsSrcIntegrators_.Resize(subdoms_.GetSize());
    fncname_rhs_.Resize(subdoms_.GetSize());
    rhsIntegrators_.Resize(subdoms_.GetSize());
    for (UInt i=0; i<subdoms_.GetSize();i++) {
      integrators_[i]      = new StdVector<IntegratorDescriptor *>;
      rhsSrcIntegrators_[i]     = new StdVector<BaseIntDescriptor *>;
      rhsIntegrators_[i]        = new StdVector<BaseIntDescriptor *>;
    }

    surfintegrators_.Resize(surfdoms_.GetSize());
    rhsSrcSurfIntegrators_.Resize(surfdoms_.GetSize());
    fncname_rhsSurf_.Resize(surfdoms_.GetSize());
    for (UInt i=0; i<surfdoms_.GetSize();i++) {
      surfintegrators_[i]       = new StdVector<IntegratorDescriptor *>;
      rhsSrcSurfIntegrators_[i] = new StdVector<BaseIntDescriptor *>;
    }
    
  }
  

  // ********************
  //   SetupMatrixGraph
  // ********************
  void Assemble::SetupMatrixGraph( bool insertCounterPart ) {

    ENTER_FCN( "Assemble::SetupMatrixGraph" );
    
    SETPROFILE( "Before SetupMatrixGraph" );

    // By default we currently assume that also the "transpose" coupling
    // object will be required in the long run
    if ( pdeId1_ != pdeId2_ ) {
      algsys_->AssembleInit( pdeId1_, pdeId2_, true );
    }
    else {
      algsys_->AssembleInit( pdeId2_, pdeId1_, false );
    }

    // set the graph - connectivity matrix
    BaseFE * ptElem; 
    UInt nsub, iel;
  
    StdVector<UInt> connecth;
    StdVector<Integer> connect_PDE1, connect_PDE2;

    // First, assemble connectivity for all volume
    // elements
    for (nsub=0; nsub<subdoms_.GetSize(); nsub++) {
      StdVector<Elem*> elemssd;
      ptgrid_->GetElems(elemssd,subdoms_[nsub]);
      
      for (iel=0; iel < elemssd.GetSize(); iel++) {  
        ptElem=elemssd[iel]->ptElem;
        connecth = elemssd[iel]->connect;
        ptEQN1_->Node2EQN(connecth, connect_PDE1);
        ptEQN2_->Node2EQN(connecth, connect_PDE2);
        
        //std::cerr << "SetElelemtPos" << std::endl 
        // << "----------" << std::endl;
        //std::cerr << "connect: " << std::endl << connecth << std::endl;
        //std::cerr << "connect_PDE " << std::endl;
        //std::cerr << connect_PDE << std::endl << std::endl;
       
       
        algsys_->SetElementPos( pdeId1_, connect_PDE1.GetPointer(),
                                connect_PDE1.GetSize(), 
                                pdeId2_, connect_PDE2.GetPointer(),
                                connect_PDE2.GetSize(),
                                insertCounterPart );
      }
    }
    
    for (nsub=0; nsub<surfdoms_.GetSize(); nsub++) {
      StdVector<SurfElem*> elemssd;
      ptgrid_->GetSurfElems(elemssd,surfdoms_[nsub]);
      
      for (iel=0; iel < elemssd.GetSize(); iel++) {
        ptElem=elemssd[iel]->ptElem;
        connecth = elemssd[iel]->connect;
        ptEQN1_->Node2EQN(connecth, connect_PDE1);
        ptEQN2_->Node2EQN(connecth, connect_PDE2);
        
        //std::cerr << "SetElelemtPos" << std::endl << "----------" << std::endl;
        //std::cerr << "connect: " << std::endl << connecth << std::endl;
        //std::cerr << "connect_PDE " << std::endl;
        //std::cerr << connect_PDE << std::endl << std::endl;
        
        algsys_->SetElementPos( pdeId1_, connect_PDE1.GetPointer(),
                                connect_PDE1.GetSize(), 
                                pdeId2_, connect_PDE2.GetPointer(),
                                connect_PDE2.GetSize(),
                                insertCounterPart );
      }
    }
    
    SETPROFILE( "After SetupMatrixGraph" );

    // finish assembling procedure (see note for AssembleInit(), above)
    if ( pdeId1_ != pdeId2_ ) {
      algsys_->AssembleDone( pdeId1_, pdeId2_, true );
    }
    else {
      algsys_->AssembleDone( pdeId1_, pdeId2_, false );
    }
    SETPROFILE( "After AssembleDone" );
  }


  //set all FE-Elements to reduced integration  
  void Assemble::SetFE2ReducedInt() {
    ptQ1    -> SetReducedIntegration();
    ptQ2    -> SetReducedIntegration();
    ptTet1  -> SetReducedIntegration();
    ptL1    -> SetReducedIntegration();
    ptL2    -> SetReducedIntegration();
    ptTr1   -> SetReducedIntegration();
    ptTr2   -> SetReducedIntegration();
    ptHexa1 -> SetReducedIntegration();
    ptHexa2 -> SetReducedIntegration();
    //    ptPyra1 -> SetReducedIntegration();
    //    ptWedge1-> SetReducedIntegration();
  }

  //set all FE-Elements to reduced integration  
  void Assemble::SetFE2StandardInt()
  {
    ptQ1    -> SetStandardIntegration();
    ptQ2    -> SetStandardIntegration();
    ptTet1  -> SetStandardIntegration();
    ptL1    -> SetStandardIntegration();
    ptL2    -> SetStandardIntegration();
    ptTr1   -> SetStandardIntegration();
    ptTr2   -> SetStandardIntegration();
    ptHexa1 -> SetStandardIntegration();
    ptHexa2 -> SetStandardIntegration();
    //    ptPyra1 -> SetStandardIntegration();
    //    ptWedge1-> SetStandardIntegration();
  }


  // define integrators
  void Assemble::AddRhsIntegrator(BaseForm * integrator,
                                  const RegionIdType  regionId, 
                                  const Boolean nonLin)
  {
    ENTER_FCN( "Assemble::AddRhsIntegrator" );
    BaseIntDescriptor * actRhsID = new  BaseIntDescriptor(integrator, nonLin);
    rhsIntegrators_[SubDomIndex(regionId)]->Push_back(actRhsID);
  }


  //needed for static and transient analysis
  void Assemble::AddRhsSrcIntegrator(BaseForm * integrator,
                                     const RegionIdType regionId, 
                                     const std::string fncname,
                                     const Boolean nonLin)
  {
    ENTER_FCN( "Assemble::AddRhsSrcIntegrator" );
    BaseIntDescriptor * actRhsID = new  BaseIntDescriptor(integrator, nonLin);
    rhsSrcIntegrators_[SubDomIndex(regionId)]->Push_back(actRhsID);
    fncname_rhs_[SubDomIndex(regionId)] = fncname;
  }

  //needed for harmonic analysis
  void Assemble::AddRhsSrcIntegrator(BaseForm * integrator,
                                     const RegionIdType regionId, 
                                     const Double phaseval,
                                     const Boolean nonLin)
  {
    ENTER_FCN( "Assemble::AddRhsSrcIntegrator" );
    BaseIntDescriptor * actRhsID = new  BaseIntDescriptor(integrator, nonLin);
    rhsSrcIntegrators_[SubDomIndex(regionId)]->Push_back(actRhsID);
    rhsSrcPhase_[SubDomIndex(regionId)] = phaseval;
  }


  //needed for static and transient analysis
  void Assemble::AddRhsSrcSurfIntegrator(BaseForm * integrator,
                                         const RegionIdType regionId,
                                         const std::string fncname,
                                         const Boolean nonLin)
  {
    ENTER_FCN( "Assemble::AddRhsSrcSurfIntegrator" );
    
    BaseIntDescriptor * actRhsID = new  BaseIntDescriptor(integrator, nonLin);
    rhsSrcSurfIntegrators_[SurfDomIndex(regionId)]->Push_back(actRhsID);
    fncname_rhsSurf_[SurfDomIndex(regionId)] = fncname;
  }

  //needed for harmonic analysis
  void Assemble::AddRhsSrcSurfIntegrator(BaseForm * integrator,
                                         const RegionIdType regionId, 
                                         const Double phaseval,
                                         const Boolean nonLin)
  {
    ENTER_FCN( "Assemble::AddRhsSrcSurfIntegrator" );
    BaseIntDescriptor * actRhsID = new  BaseIntDescriptor(integrator, nonLin);
    rhsSrcSurfIntegrators_[SurfDomIndex(regionId)]->Push_back(actRhsID);
    rhsSrcSurfPhase_.Resize(1);
    rhsSrcSurfPhase_[SurfDomIndex(regionId)] = phaseval;
  }


  // ==========================================================
  // STATIC ANALYSIS
  // ==========================================================


  StaticAssemble::StaticAssemble(BaseSystem * algsys, Grid * agrid)
    :Assemble(algsys, agrid)
  {
    ENTER_FCN( "StaticAssemble::StaticAssemble" );
    SetAnalysisType(STATIC);
  }
  

  // define integrators
  void StaticAssemble::AddIntegrator(IntegratorDescriptor * actID,
                                     const RegionIdType regionId) {

    ENTER_FCN( "StaticAssemble::AddIntegrator" );

    if (actID->DestMat() == STIFFNESS)
      actID->SetDestMat(SYSTEM);

    // in static analysis, neglect damping and mass contributions
    if (actID->DestMat() !=  SYSTEM)
      return;

    integrators_[SubDomIndex(regionId)]->Push_back(actID);
  }
    

  // define surface integrators
  void StaticAssemble::AddSurfIntegrator( IntegratorDescriptor * actID,
                                          const RegionIdType regionId) {

    ENTER_FCN( "StaticAssemble::AddSurfIntegrator" );
   
    if (actID->DestMat() == STIFFNESS)
      actID->SetDestMat(SYSTEM);

    // in static analysis, neglect damping and mass contributions
    if (actID->DestMat() !=  SYSTEM)
      return;

    surfintegrators_[SurfDomIndex(regionId)]->Push_back(actID);
  }


  // ==========================================================
  // TRANSIENT ANALYSIS
  // ==========================================================


  TransientAssemble::TransientAssemble(BaseSystem * algsys, Grid * agrid)
    :Assemble(algsys, agrid)
  {
    ENTER_FCN( "TransientAssemble::TransientAssemble" );
    SetAnalysisType(TRANSIENT);

  }


  // define integrators
  void TransientAssemble::AddIntegrator( IntegratorDescriptor * actID,
                                         const RegionIdType regionId ) {
    ENTER_FCN( "TransientAssemble::AddIntegrator" );

    if (actID->DestMat() == SYSTEM) {
      Info->Error("In transient assembling, no SYSTEM matrix may be defined directly", __FILE__, __LINE__);
    }
    integrators_[SubDomIndex(regionId)]->Push_back(actID);
    
    // Pass needed FE-matrix types to the algebraic system
    algsys_->SetFEMatrixType(actID->DestMat(),
                             actID->GetPDE1()->GetPDEId(),
                             actID->GetPDE2()->GetPDEId());

    algsys_->SetFEMatrixType(actID->GetSecondaryMat(),
                             actID->GetPDE1()->GetPDEId(),
                             actID->GetPDE2()->GetPDEId());
  }

  // define integrators
  void TransientAssemble::AddSurfIntegrator( IntegratorDescriptor *actID,
                                             const RegionIdType regionId ) {
    ENTER_FCN( "TransientAssemble::AddSurfIntegrator" );

    if ( actID->DestMat() == SYSTEM ) {
      (*error) << "In transient assembling, no SYSTEM matrix may be "
               << "defined directly";
      Error( __FILE__, __LINE__ );
    }

    surfintegrators_[SurfDomIndex(regionId)]->Push_back(actID);

    // Pass needed FE-matrix types to the algebraic system
    algsys_->SetFEMatrixType(actID->DestMat(),
                             actID->GetPDE1()->GetPDEId(),
                             actID->GetPDE2()->GetPDEId());
    
    algsys_->SetFEMatrixType(actID->GetSecondaryMat(),
                             actID->GetPDE1()->GetPDEId(),
                             actID->GetPDE2()->GetPDEId());
  }



  // ==========================================================
  // BASE INTEGRATOR DESCRIPTOR
  // ==========================================================

  BaseIntDescriptor::BaseIntDescriptor() : 
    integrator(NULL),
    nonLin(FALSE),
    reducedIntegration_(FALSE),
    setCounterPart_(true) {

    ENTER_FCN( "BaseIntDescriptor::BaseIntDescriptor" );

  }

  BaseIntDescriptor::~BaseIntDescriptor() {
    ENTER_FCN( "BaseIntDescriptor::~BaseIntDescriptor" );
    delete integrator;
  }

  // define integrators
  BaseIntDescriptor::BaseIntDescriptor(BaseForm * aIntegrator,
                                       const Boolean aNonLin)
    : integrator(aIntegrator), nonLin(aNonLin), 
      reducedIntegration_(FALSE), setCounterPart_(true)
  {
    ENTER_FCN( "BaseIntDescriptor::BaseIntDescriptor " );
  }
  

  // ==========================================================
  // INTEGRATOR DESCRIPTOR
  // ==========================================================

  IntegratorDescriptor::IntegratorDescriptor()
    :BaseIntDescriptor(),
     piezoMaterialType_(REALMATERIALPARAMETER),
     destinationMatrix(SYSTEM),
     secondaryMatrix(NOTYPE),
     secMatFac(0.0)
  
  {
    ENTER_FCN( "IntegratorDescriptor::IntegratorDescriptor" );
    //    piezoMaterialType_=realMaterialParameter;
  }


  // define integrators
  IntegratorDescriptor::IntegratorDescriptor(BaseForm * aIntegrator,
                                             const FEMatrixType aDestMat,
                                             const Boolean aNonLin)
    :BaseIntDescriptor(aIntegrator, aNonLin),
     piezoMaterialType_(REALMATERIALPARAMETER),
     destinationMatrix(aDestMat),
     secondaryMatrix(NOTYPE),
     secMatFac(0.0)
  {
    ENTER_FCN( "IntegratorDescriptor::IntegratorDescriptor" );
    //   piezoMaterialType_=realMaterialParameter;
  }

  // defines a secondary destination for the calculated element marices of an
  // integrator
  void IntegratorDescriptor::SetSecondaryMat( FEMatrixType aSecMat,
                                              Double aSecMatFac,
                                              AnalysisType analysisType ) {

    ENTER_FCN( "IntegratorDescriptor::SetSecondaryMat" );
    FEMatrixType MatType = aSecMat;

    if (analysisType == HARMONIC||analysisType == MULTIHARMONIC) {
      if (aSecMat == STIFFNESS || aSecMat == MASS || aSecMat == DAMPING )
        MatType = SYSTEM;
      else
        {
          std::string error_msg = "Matrix type ";
          error_msg += aSecMat + " not supported in harmonic analysis";
          Error(error_msg.c_str(), __FILE__, __LINE__ );
        }

      SetOrigSecMatrixType(aSecMat);
    }

    else if ( analysisType == TRANSIENT ) {
      secondaryMatrix = MatType;
      secMatFac = aSecMatFac;
    }

  }


  // **************
  //   Destructor
  // **************
  IntegratorDescriptor::~IntegratorDescriptor() {
    ENTER_FCN( "IntegratorDescriptor::~IntegratorDescriptor" );
  }


  // ==========================================================
  // HARMONIC ANALYSIS
  // ==========================================================


  HarmonicAssemble::HarmonicAssemble(BaseSystem * algsys, Grid * agrid)
    :Assemble(algsys, agrid)
  {
    ENTER_FCN( "HarmonicAssemble::HarmonicAssemble" );
    SetAnalysisType(HARMONIC);

    //get start frequency of harmonic analysis
    StdVector<std::string> keyVec, attrVec, valVec;
    Double freq;

    attrVec = "tag";
    valVec  = "anyTag";
  
    // Get time stepping information from parameter object
    keyVec = "harmonic", "startFreq";
    params->Get( keyVec, attrVec, valVec, freq );

    startFreq_ = 2*PI*freq;
  }

  /// set actual frequency (already multiplied by 2*pi)
  void HarmonicAssemble::SetFrequency(Double frequency)
  {
    ENTER_FCN( "HarmonicAssemble::SetFrequency" );

    actFreq_ = 2*PI*frequency;

  } 


  // define integrators
  void HarmonicAssemble::AddIntegrator(IntegratorDescriptor * actID,
                                       const RegionIdType regionId)
  {
    ENTER_FCN( "HarmonicAssemble::AddIntegrator" );

    actID->SetOrigMatrixType(actID->DestMat());
    if (actID->DestMat() == STIFFNESS 
        || actID->DestMat() == MASS 
        || actID->DestMat() == DAMPING )
      actID->SetDestMat(SYSTEM);
    else
      {
        std::string error_msg = "Matrix type ";
        error_msg += actID->DestMat() + " not supported in harmonic analysis";
        Error(error_msg.c_str(), __FILE__, __LINE__ );
      }

    integrators_[SubDomIndex(regionId)]->Push_back(actID);

  }

  // define integrators
  void HarmonicAssemble::AddSurfIntegrator(IntegratorDescriptor * actID,
                                           const RegionIdType regionId) {
    ENTER_FCN( "HarmonicAssemble::AddSurfIntegrator" );
    
    actID->SetOrigMatrixType(actID->DestMat());
    if (actID->DestMat() == STIFFNESS 
        || actID->DestMat() == MASS 
        || actID->DestMat() == DAMPING )
      actID->SetDestMat(SYSTEM);
    else
      {
        std::string error_msg = "Matrix type ";
        error_msg += actID->DestMat() + " not supported in harmonic analysis";
        Error(error_msg.c_str(), __FILE__, __LINE__ );
      }

    surfintegrators_[SurfDomIndex(regionId)]->Push_back(actID);
  }


  void HarmonicAssemble::
  TransformMatrix2Harmonic( Vector<Double>& harmMat,
                            Matrix<Double> origMat,
                            const FEMatrixType matrixType,
                            const piezoMaterialType piezoMatType ) {

    ENTER_FCN( "HarmonicAssemble::TransformMatrix2Harmonic" );

    Integer numRow = origMat.GetSizeRow();
    Integer numCol = origMat.GetSizeCol();
    harmMat.Resize(2*numRow*numCol);

    Integer k=0;

    if (piezoMatType == REALMATERIALPARAMETER) {

      if (matrixType == STIFFNESS) {
        //	std::cout << "Assemble real stiff\n" << std::endl;

        for (Integer row=0; row<numRow; row++)
          for (Integer col=0; col<numCol; col++) {
            harmMat[k] = origMat[row][col];
            k++;
          }
      }
    
      else if (matrixType == MASS)
        {
          //          std::cout<<"Assemble real_mass\n"<<std::endl;

          Double factor = -actFreq_*actFreq_;
          for (Integer row=0; row<numRow; row++)
            for (Integer col=0; col<numCol; col++) {
              harmMat[k] = factor*origMat[row][col];
              k++;
            }
        }
    
      else if (matrixType == DAMPING)
        {       
          //    std::cout<<"real_damping"<<std::endl;
          Double factor = actFreq_;
       
          k=numRow*numCol;
          for (Integer row=0; row<numRow; row++)
            for (Integer col=0; col<numCol; col++) {
              harmMat[k] = factor*origMat[row][col];
              k++;
            }
        }
    } // end, if piezoMatType == real...

    else if(piezoMatType == IMAGMATERIALPARAMETER){  // the "imaginary parts"
   
      if (matrixType == STIFFNESS)
        {
          //          std::cout<<"Assmble imag_stiff\n"<< std::endl;
          k=numRow*numCol;
          for (Integer row=0; row<numRow; row++)
            for (Integer col=0; col<numCol; col++) {
              harmMat[k] = origMat[row][col];
              k++;
            }
        }
       
    
      else if (matrixType == MASS)
        {
          //          std::cout<<"Assemble imag_mass\n"<<std::endl;

          Double factor = -actFreq_*actFreq_;
	  k=numRow*numCol;
          for (Integer row=0; row<numRow; row++)
            for (Integer col=0; col<numCol; col++) {
              harmMat[k] = factor*origMat[row][col];
              k++;
            }
        }

      else if (matrixType == DAMPING)
        {
          Double factor = actFreq_;
          //    std::cout<<"comlex_damping"<<std::endl;
        
          k=0;  
          for (Integer row=0; row<numRow; row++)
            for (Integer col=0; col<numCol; col++) {
              harmMat[k] = -factor*origMat[row][col];
              k++;
            }
        }

    } // end if piezoMatType == imag

    else {
      (*error) <<"\n piezoMaterialType" << piezoMatType 
               << "not specified "<<std::endl;
      Error( __FILE__, __LINE__ );
    }
  }



  void  HarmonicAssemble::TransformVector2Harmonic(Vector<Double>& harmVec, 
                                                   Vector<Double> origVec, 
                                                   const Double valPhase)
  {
    ENTER_FCN( "HarmonicAssemble::TransformVector2Harmonic" );

    Integer size = origVec.GetSize();
    harmVec.Resize(2*size);

    Double valReal = cos(valPhase);
    Double valImag = sin(valPhase);

    Integer k=0;
    //real part
    for (Integer i=0; i<size; i++) {
      harmVec[k] = origVec[i]*valReal;
      k++;
    }

    //imaginary part
    for (Integer i=0; i<size; i++) {
      harmVec[k] = origVec[i]*valImag;
      k++;
    }

  }


} // end namespace CoupledField
