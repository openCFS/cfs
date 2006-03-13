#include <fstream>
#include <iostream>
#include <string>
#include <math.h>

#include "Domain/elem.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "DataInOut/WriteInfo.hh"
#include "PDE/StdPDE.hh"

#include "olas.hh"
#include "assemble.hh"
#include "DataInOut/MHMaterialData.hh"
#include "MHassemble.hh"

namespace CoupledField
{
  
  MHassemble :: MHassemble(BaseSystem * algsys, Grid * agrid,
                           const std::string bcSequenceTag )
    :Assemble(algsys, agrid, bcSequenceTag )
  {
    
    ENTER_FCN( "MHassemble::MHassemble" );

    algsys_= algsys;
    //ptgrid_=agrid;
    MHpdeId1_=pdeId1_;
    MHpdeId2_=pdeId2_;

    // SetAnalysisType(HARMONIC);
    firstTime_ = TRUE;
    oneIntIsNonlin_ = FALSE;
    nrMatrices_ = 5+1;  // 5 matrices, but index starts with 1!
    reassembleMat_.Resize(nrMatrices_);
    nonLinGeo = FALSE;
    deltaCoords_ = NULL;
    alternateMaterialData_=FALSE;

    startFreq_ = 0.0;
    //    matArray_  = NULL;

    StdVector<std::string> keyVec, attrVec, valVec;
    attrVec = "tag";
    valVec = "multiHarmonic";
    keyVec = "multiHarmonic", "nrMultiHarmonics";
    params->Get(keyVec, attrVec, valVec, nrMultHarms_);

    SetAnalysisType(MULTIHARMONIC);
    //    std::cout<<" end of construction MHAssemble"<<std::endl;



  }

  MHassemble :: ~MHassemble(){
    ENTER_FCN( "MHassemble::~MHassemble" );
  }

//   void MHassemble::MHAssembleMatrices(StdVector<RegionIdType> subdoms_, 
//                                       StdVector< StdVector<IntegratorDescriptor *>* > integrators_,
//                                       NodeEQN * ptEQN1_,
//                                       NodeEQN * ptEQN2_){
//     ENTER_FCN( "MHassemble:MHAssembleMatrices" );
//     std::cout<< "MHassemble:MHAssembleMatrices" <<std::endl;



   void MHassemble::AssembleMatrices(){
     ENTER_FCN( "MHassemble:AssembleMatrices" );
     std::cout<< "MHassemble:AssembleMatrices" <<std::endl;
    
     Vector<Double> harmonicVec;
     Double dampTransform = 1.0;
     
    // initialize reassembling "indicator" vector
    for (Integer MHIndRow=-nrMultHarms_; MHIndRow<=nrMultHarms_;MHIndRow++){
      for (Integer MHIndCol=-nrMultHarms_; MHIndCol<=nrMultHarms_;MHIndCol++){
          
        if(firstTime_)
          for (UInt actMat=0; actMat < nrMatrices_; actMat++)
            reassembleMat_[actMat] = TRUE;
        
        for (UInt actDom=0; actDom < subdoms_.GetSize(); actDom++) {     
          StdVector<Elem*> elemssd;

          ptgrid_->GetVolElems(elemssd, subdoms_[actDom]);

     //      std::cout<<"integrators_ = " <<std::endl;
//           std::cout<< integrators_[actDom]->GetSize()  <<std::endl;

          for(UInt actInteg=0; actInteg < integrators_[actDom]->GetSize(); actInteg++) {
            
            IntegratorDescriptor * actDescriptor = 
              (*integrators_[actDom])[actInteg];

            FEMatrixType destMat = actDescriptor->GetIntegrator()->GetBaseType();

            if (((MHIndRow==MHIndCol) && (destMat == MASS)) || (destMat==STIFFNESS) || (destMat==DAMPING)){ 
              // Assemble only diagonal entries in bigM and for all entries in bigK
              std::cout<<"\n Destination Matrix = "<< destMat<<std::endl;
              std::cout<<"has position in SuperMatrix = "<< MHIndRow<<" x " <<MHIndCol<<std::endl;

              actDescriptor->GetIntegrator()->SetSubdomain(actDom);

              if (alternateMaterialData_ == TRUE)
                actDescriptor->GetIntegrator()->SetMaterial(ptMaterial_);
                    
              // assemble only if nonlinear or first time
              if (reassembleMat_[actDescriptor->DestMat()] ||firstTime_) {

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
                    Info->PrintF( "", " dampTransform (stiffness matrix) = %e\n",
                                  dampTransform );
                  }
                  else if ( destMat == MASS ) {
                    dampTransform = actFreq_ / matDataFreq;
                    Info->PrintF( "", " dampTransform (mass matrix) = %e\n",
                                  dampTransform );
                  }
                }

                if ( actDescriptor->GetIntegrator()->IsRaylDamping() ) {
                  actDescriptor->GetIntegrator()->SetFactor(dampTransform);
                }
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

         //        std::cout<<"connecth: "<<std::endl;
//                 std::cout<<connecth<<std::endl;
                    
                // map connect to PDE node numbers
                StdVector<Integer> connect_PDE1, connect_PDE2;
                    
                ptEQN1_->Node2EQN(connecth, connect_PDE1);
                ptEQN2_->Node2EQN(connecth, connect_PDE2);
                    
                Matrix<Double> elSol;
                    
                actDescriptor->GetIntegrator()->SetElemPtr(ptEl);
                FEMatrixType destFEMMat = actDescriptor->DestMat();
           
                
                // this matrix is nonlinear and, therefore, 
                // has to be reassembled next time
                if (actDescriptor->IsNonLin()) {
                  oneIntIsNonlin_ = TRUE;
                  reassembleMat_[actDescriptor->DestMat()] = TRUE;
                  sol_->GetElemSolutionAsMatrix(elSol, connecth);
                  actDescriptor->GetIntegrator()->SetActElemSol(elSol);
                }       
           
                //------------------------------------------------------------------------------
                // ================================================================
                //                             assemble matrices
                // ================================================================

                // spätestens hier Material updaten

              
                if (destMat==20){
                  std::cout<<" New Material will be calculated ... " <<std::endl;
                  //  getchar();

//                   //                  if (= TRUE)
                  ptMaterial_ = actDescriptor->GetIntegrator()->GetMaterial();
                  ptMHMat_ = new MHMaterialData(ptMaterial_);

//                   std::cout<<" parameter before calcParameterCurveAtEl .... "<<ptMHMat_->parameter_<<std::endl;
//                   std::cout<<ptMHMat_->parameter_<<std::endl;
//                   std::cout<<ptMHMat_->parameterCoeff_<<std::endl;

                  ptMHMat_->calcParameterCurveAtElement(ptMHMat_->parameter_,ptMHMat_->parameterCoeff_, 10, nrMultHarms_, -1, 2);
                  //                  std::cout<<ptMHMat_->parameter_<<std::endl;

//                   std::cout<<" parameter after calcParameterCurveAtEl .... "<<std::endl;
//                   std::cout<<ptMHMat_->parameter_<<std::endl;
                  
                  ptMHMat_->updateMaterialData(ptMHMat_->parameter_,ptMaterial_);
//                   //              SetAlternatingMaterial(TRUE);
                  
                }
                //                std::cout<<"geht es hier weiter ? "<<std::endl;
                //                getchar();

                // ----------------------------------------------------------------------

                actDescriptor->GetIntegrator()->
                  CalcElementMatrix(ptCoord, elemmat);
//                 std::cout<<elemmat<<std::endl;
//                 getchar();
                  
                piezoMaterialType matType = actDescriptor->GetPiezoMaterialType();
                actDescriptor->SetPiezoMaterialType(matType);

                if (TRUE){
                  TransformMatrix2Harmonic(harmonicVec,elemmat, 
                                           actDescriptor->GetOrigMatrixType(),
                                           actDescriptor->GetPiezoMaterialType());
                  //                  std::cout<<"Transformation to harmonic vec went fine "<<std::endl;
                  //                  std::cout<<harmonicVec<<std::endl;
                    
                  if (destMat== MASS){
                    //                    std::cout<<" Now ... I will set " << destMat << " = MassMatrix "<< std::endl;
                      
                    harmonicVec = harmonicVec*Double(std::abs(MHIndRow));
                    
                    
                    algsys_->SetElementMatrix( destFEMMat, &harmonicVec[0], 
                                               pdeId1_, connect_PDE1.GetPointer(), 
                                               connect_PDE1.GetSize(),
                                               pdeId2_, connect_PDE2.GetPointer(), 
                                               connect_PDE2.GetSize() );
                  }
                  else if (destMat  == STIFFNESS) {
                    // std::cout<<" Now ... I will set elementMatrix ... "<<std::endl;
//                     getchar();
                    algsys_->SetElementMatrix( destFEMMat, &harmonicVec[0], 
                                                  pdeId1_, connect_PDE1.GetPointer(), 
                                                  connect_PDE1.GetSize(),
                                                  pdeId2_, connect_PDE2.GetPointer(), 
                                                  connect_PDE2.GetSize() );
                    

                            
              
                       if (actDescriptor->GetSecondaryMat() != NOTYPE) {
                     Double damp = dampTransform * actDescriptor->GetSecMatFac();
                     elemmat *= damp;

                     TransformMatrix2Harmonic(harmonicVec,elemmat,
                                              actDescriptor->GetOrigSecMatrixType(),
                                              actDescriptor->GetPiezoMaterialType());
        
                     algsys_->SetElementMatrix(destFEMMat, &harmonicVec[0], 
                                               pdeId1_, connect_PDE1.GetPointer(),
                                               connect_PDE1.GetSize(), 
                                               pdeId2_, connect_PDE2.GetPointer(),
                                               connect_PDE2.GetSize() );
                       }

                  }
                }
                 
                
            
              } //over all elements of subdomain            
                
              // } //check, if we have to assemble
          
              if ( actDescriptor->IsReducedInt()) {
                //set all elements back to standard integration!!
                SetFE2StandardInt();
              }
            } //end of if MHIndRow==MHIndCol
       
          } //integrators
        
        } //subdomains

      }// multiharmonics rows
    }// multiharmonics cols

  }// MHAssemble Matrices
    
 //  void MHassemble::SetPtr2EQNData(NodeEQN * aPtNodeEQN1,
//                                 NodeEQN * aPtNodeEQN2)
//   {
//     ENTER_FCN( "MHassemble::SetPtr2EQNData" );

//     if (aPtNodeEQN1 == NULL)
//       Error( "The NodeEQN object is not created yet" , 
//              __FILE__, __LINE__);
    
//     if (! aPtNodeEQN1->IsInitialized())
//       Error( "The NodeEQN object has to be initialized before assigning it to\
//               an assemble object by calling 'CalcMapping()'", 
//              __FILE__, __LINE__);
//     ptEQN1_ = aPtNodeEQN1;

//     // check if a second eqn-object is present.
//     // If not, the same is used two times
//     if (aPtNodeEQN2 == NULL)
//       ptEQN2_ = aPtNodeEQN1;
//     else
//       ptEQN2_ = aPtNodeEQN2;
    
//   }


  /// set actual frequency (already multiplied by 2*pi)
  void MHassemble::SetFrequency(Double frequency)
  {
    ENTER_FCN( "MHassemble::SetFrequency" );

    actFreq_ = 2*PI*frequency;

  } 

  void  MHassemble::TransformVector2Harmonic(Vector<Double>& harmVec, 
                                                   Vector<Double> origVec, 
                                                   const Double valPhase)
  {
    ENTER_FCN( "MHassemble::TransformVector2Harmonic" );

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

  void MHassemble :: TransformMatrix2Harmonic( Vector<Double>& harmMat,
                            Matrix<Double> origMat,
                            const FEMatrixType matrixType,
                            const piezoMaterialType piezoMatType ) {

    ENTER_FCN( "MHAssemble::TransformMatrix2Harmonic" );

    Integer numRow = origMat.GetSizeRow();
    Integer numCol = origMat.GetSizeCol();
    harmMat.Resize(2*numRow*numCol);

    Integer k=0;

    if (piezoMatType == REALMATERIALPARAMETER) {

      if (matrixType == STIFFNESS) {

        // std::cout<<"real_stiff - actfreq: "<< actFreq_<<std::endl;
        for (Integer row=0; row<numRow; row++)
          for (Integer col=0; col<numCol; col++) {
            harmMat[k] = origMat[row][col];
            k++;
          }
      }
    
      else if (matrixType == MASS)
        {
          //std::cout<<"real_mass"<<std::endl;

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
          //std::cout<<"comlex_stiff"<<std::endl;
          k=numRow*numCol;
          for (Integer row=0; row<numRow; row++)
            for (Integer col=0; col<numCol; col++) {
              harmMat[k] = origMat[row][col];
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


  void MHassemble::AddIntegrator(IntegratorDescriptor * actID,
                                 const RegionIdType regionId)
  {
    ENTER_FCN( "MHassemble::AddIntegrator" );

    actID->SetOrigMatrixType(actID->DestMat());
    if (actID->DestMat() == STIFFNESS 
        || actID->DestMat() == MASS 
        || actID->DestMat() == DAMPING )
      actID->SetDestMat(SYSTEM);
    else {
      (*error) << "Matrix type " << actID->DestMat()
               << " is not supported in harmonic analysis";
      Error( __FILE__, __LINE__ );
    }

    integrators_[SubDomIndex(regionId)]->Push_back(actID);

  }


  void MHassemble::AddSurfIntegrator(IntegratorDescriptor * actID,
				     const RegionIdType regionId)
  {
    ENTER_FCN( "MHAssemble::AddSurfIntegrator" );

    actID->SetOrigMatrixType(actID->DestMat());
    if (actID->DestMat() == STIFFNESS 
        || actID->DestMat() == MASS 
        || actID->DestMat() == DAMPING ) {
      actID->SetDestMat(SYSTEM);
    }
    else {
      (*error) << "Matrix type " << actID->DestMat()
               << " is not supported in harmonic analysis";
      Error( __FILE__, __LINE__ );
    }

    surfintegrators_[SurfDomIndex(regionId)]->Push_back(actID);
  }
  
}




