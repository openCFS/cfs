// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "assemble.hh"

#include <iostream>
#include <iomanip>
#include <boost/lexical_cast.hpp>

#include "Domain/entityList.hh"
#include "Forms/baseForm.hh"
#include "Forms/linearForm.hh"
#include "PDE/SinglePDE.hh"
#include "Domain/domain.hh"
#include "Domain/grid.hh"
#include "Driver/stdSolveStep.hh"
#include "Driver/harmonicDriver.hh"
#include "DataInOut/programOptions.hh"
#include "DataInOut/Logging/cfslog.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "OLAS/algsys/algebraicSys.hh"
#include "DataInOut/Scripting/cfsmessenger.hh"
#include <boost/progress.hpp>
#include "Utils/Timer.hh"
#include "Optimization/Optimization.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "Materials/mechanicMaterial.hh"
#include "DataInOut/ResultCache.hh"

namespace CoupledField
{
  // declare logging stream
  DECLARE_LOG(assemble)
  DEFINE_LOG(assemble, "assemble")

  Assemble::Assemble( AlgebraicSys* algsys,
                      BasePDE::AnalysisType analysis,
                      UInt maxTimeDerivOrder ) : timer_(new Timer()) {

    // init general params
    algsys_ = algsys;
    analysisType_ = analysis;
    isFirstTime_ = true;
    matrixUpdated_ = false;
    maxTimeDerivOrder_ = maxTimeDerivOrder;
    
    linForms_ = new StdVector<LinearFormContext*>();

    // Calculate matrix map from general matrix types to analysis
    // specific ones
    CreateMatrixMap();

    // Initialize reassemble-map
    ResetMatrixReassembly();

    // Obtain a new mathParser handler
    mHandle_ = domain->GetMathParser()->GetNewHandle();
    

    // the timer object is used in every AssembleMatrices() call
    info->Get("analysis")->Get(ParamNode::SUMMARY)->Get("assemble/timer")->SetValue(timer_);

    // Initialize scripting interface
    RegisterFunctions();
  }

  Assemble::~Assemble() {

    // Delete bilinear contexts
    std::set<BiLinFormContext*>::iterator it = allBiLinForms_.begin();
    for( ; it != allBiLinForms_.end(); ++ it){
      delete (*it);
    }
    

    // Delete linear contexts
    for(unsigned int i = 0; i < linForms_->GetSize(); ++i){
      delete (*linForms_)[i];
    }
    
    delete linForms_;
  }

  void Assemble::SetAlgSys(AlgebraicSys * algsys)  {
    algsys_ = algsys;
  }

  void Assemble::ResetMatrixReassembly() {
    // Initialize reassemble-map
    // Note: in the beginning all matrices have to be "re-"assembled
    // which means the destination matrices have to be cleared and
    // initialized. After the first "AssmebleMatrices" call the
    // correct values for reassembling are determined
    matReassemble_[STIFFNESS] = true;
    matReassemble_[DAMPING] = true;
    matReassemble_[MASS] = true;
    matReassemble_[AUXILIARY] = true;

    // reset also flag for "firstTime"
    isFirstTime_ = true;
    matrixUpdated_ = true;
  }

  
  BiLinFormContext* Assemble::GetBiLinForm(RegionIdType regionId, StdPDE* pde1, StdPDE* pde2,  const std::string& integrator, bool silent)
  {
    REFACTOR;
     // the EntityList has the region name as name but not the id
//     std::string region = domain->GetGrid()->GetRegion().ToString(regionId);

     BiLinFormContext* result = NULL;

//     // iterate over all descriptors
//     StdVector<BiLinFormContext*>::iterator iter;
//     int counter = 0;
//     for (iter = biLinForms_->Begin(); iter != biLinForms_->End(); iter++)
//     {
//       counter++;
//       // we are wrong if the region does not match
////       if((*iter)->GetFirstEntities()->GetName() != region) continue;
//       if((*iter)->GetFirstEntities()->GetRegion() != regionId) continue;
//       // when pde1 is given we compare it by name and continue if the names are different
////       if(pde1 != NULL && (*iter)->GetFirstPde()->GetName() != pde1->GetName()) continue;
//       if((*iter)->GetFirstPde() != pde1) continue;
////       if(pde2 != NULL && (*iter)->GetSecondPde()->GetName() != pde2->GetName()) continue;
//       if((*iter)->GetSecondPde() != pde2) continue;
//       //std::cout << counter << ":" << (*iter)->GetIntegrator()->GetName() << " vs " << integrator << std::endl;
//       if((*iter)->GetIntegrator()->GetName() != integrator) continue;
//
//       // we come here because we had no contradiction - check for uniqueness
//       if(result != NULL) throw Exception("parameters not unique!");
//       // absolutley no contradiction, save result and continue to
//       // check that there is no further match
//       result = *iter;
//     }
//
//     if(result == NULL && !silent)
//     {
//       std::string region = domain->GetGrid()->GetRegion().ToString(regionId);
//       EXCEPTION("BiLinFormContext '" << integrator << "' at region '" << region << "' not found within " << counter << " integrators");
//     }
     return result;
  }

  LinearFormContext* Assemble::GetLinearForm(RegionIdType regionId, StdPDE* pde,  const std::string& integrator, bool silent)
  {
     // the EntityList has the region name as name but not the id
     std::string region = domain->GetGrid()->GetRegion().ToString(regionId);

     LinearFormContext* result = NULL;

     // iterate over all descriptors
     for(UInt i = 0; i < linForms_->GetSize(); i++)
     {
       // we are wrong if the region does not match
       LinearFormContext* lfc = (*linForms_)[i];
       if(lfc->GetEntities()->GetName() != region) continue;
       // when pde1 is given we compare it by name and continue if the names are different
       if(lfc->GetPde()->GetName() != pde->GetName()) continue;
       if(lfc->GetIntegrator()->GetName() != integrator) continue;

       // we come here because we had no contradiction - check for uniqueness
       if(result != NULL) throw Exception("parameters not unique!");
       result = lfc;
     }

     if(result == NULL && !silent)
       EXCEPTION("LinearFormContext '" << integrator << "' at region '" << region << "' not found");
     return result;
  }

  void Assemble::AddBiLinearForm( BiLinFormContext* biLinContext ) {

    LOG_DBG(assemble) << "Adding BiLinearForm '"
                      << biLinContext->GetIntegrator()->GetName()
                      << "' on '"
                      << biLinContext->GetFirstEntities()->GetName()
                      << "'";

    // assert that Integrator is set
    assert( biLinContext->GetIntegrator() != NULL );

    // assert that some entities are set
    assert( biLinContext->GetFirstEntities() != NULL );
    assert( biLinContext->GetSecondEntities() != NULL );

    // not needed anymore, as we use thef espace to map things
    // assert that the pdes are set
//    assert( biLinContext->GetFirstPde() != NULL );
//    assert( biLinContext->GetSecondPde() != NULL );

    // If the datatype of the bilinearformcontext is "COMPLEX"
    // we have to ensure that we are in an harmonic case.
    // Otherwise we issue an error
    if( (biLinContext->GetEntryType() == Global::IMAG ||
         biLinContext->GetEntryType() == Global::COMPLEX )
        && analysisType_ != BasePDE::HARMONIC ) {
      EXCEPTION( "Can not add integrator '"
                 << biLinContext->GetIntegrator()->GetName()
                 << "' with complex/imaginary entries for a "
                 << "non-harmonic analysis." );
    }

    FEMatrixType mappedFEType = matrixMap_[biLinContext->GetDestMat()];
    FEMatrixType mappedSecFEType =
      matrixMap_[biLinContext->GetSecDestMat()];

    // Check if integrator can be assembled in this type of simulation
    if( mappedFEType != NOTYPE ) {

      
      
      
      // Store bilinear form
      allBiLinForms_.insert(biLinContext);
      
      // Note: As the shared_ptr to an Entitylist is not a unique
      // unique within CFS, we have to ensure, that the names of 
      // the entity lists rather than the pointers match!
      
      // Loop over all existing bilienarforms and check, 
      // if pair (EntityList1, EntityList2) was already defined
      std::string ent1Name = biLinContext->GetFirstEntities()->GetName();
      std::string ent2Name = biLinContext->GetSecondEntities()->GetName();
      std::pair<shared_ptr<EntityList>,shared_ptr<EntityList> > pair;
      BiLinContextListType::iterator it = biLinForms_.begin();
      bool found = false;
      for( ; it != biLinForms_.end(); ++it ) {
        if( it->first.first->GetName() == ent1Name && 
            it->first.second->GetName() == ent2Name ) {
          pair = it->first;
          found = true;
          break;
        }
      }
      if(!found) {
        pair = std::pair<shared_ptr<EntityList>,shared_ptr<EntityList> >
        (biLinContext->GetFirstEntities(), biLinContext->GetSecondEntities());
      }
      biLinForms_[pair].Push_back(biLinContext);

      // Pass needed matrix type to algebraic system
      FeFctIdType id1 = biLinContext->GetFirstFeFunction()->GetFctId();
      FeFctIdType id2 = biLinContext->GetSecondFeFunction()->GetFctId();
      algsys_->SetFEMatrixType( mappedFEType, id1, id2 );

      // Check for secondary matrix type
      if( mappedSecFEType != NOTYPE ) {
        algsys_->SetFEMatrixType( mappedSecFEType, id1, id2 );
      }
    }
    else {
      delete biLinContext;
    }
  }

  void Assemble::AddLinearForm( LinearFormContext* linContext ) {

   LOG_DBG(assemble) << "AddLinearForm: " << linContext->ToString()
                     << " on " << linContext->GetEntities()->GetName();

    // assert that Integrator is set
    assert( linContext->GetIntegrator() != NULL );

    // assert that the pdes are set
    assert( linContext->GetPde() != NULL );

    // assert that some entites are set
    assert( linContext->GetEntities() != NULL );

    linForms_->Push_back( linContext );

  }

  void Assemble::AddLoads( LoadList& load ) {

    for( UInt i=0; i<load.GetSize(); i++ ) {
      loads_.Push_back( load[i] );
    }
  }

  void Assemble::SetupMatrixGraph(FeFctIdType fctId1, FeFctIdType fctId2 ) {

    StdVector<Integer> eqnVec1, eqnVec2;
    FeFctIdType id1, id2;

    // iterate over all entitylist-pairs and 
    BiLinContextListType::iterator listIt = biLinForms_.begin();
    for ( ; listIt != biLinForms_.end(); ++listIt) {
      StdVector<BiLinFormContext*> & forms = listIt->second;
      EntityList& firstEntities = *(listIt->first.first);
      EntityList& secondEntities = *(listIt->first.second);

      
      // Loop over all entities
      EntityIterator it1 = firstEntities.GetIterator();
      EntityIterator it2 = secondEntities.GetIterator();
      for( it1.Begin(); !it1.IsEnd(); it1++, it2++ ) {
      
        // Loop over all bilinearforms
        for( UInt iForm = 0; iForm < forms.GetSize(); ++iForm ) {
          
        // get integrator
        BiLinFormContext & actContext = *forms[iForm];

        FEMatrixType destMap = matrixMap_[actContext.GetDestMat()];

        // The bilinearform gets set if
        // a) the two pde-Ids match in the same order
        // b) the two pde-Ids match in the reverse order
        //    and the counter part has to be set
        bool doAssemble = false;
        bool doTranspose = true;
        bool setCounterPart = false;

        if( fctId1 == fctId2 ) {
          setCounterPart = actContext.IsSetCounterPart();
        }

        // case a)
        if( actContext.GetFirstFeFunction()->GetFctId() == fctId1 &&
            actContext.GetSecondFeFunction()->GetFctId() == fctId2 ) {
          doAssemble = true;
          doTranspose = false;
        } else

          //case b)
          if( actContext.GetFirstFeFunction()->GetFctId() == fctId1 &&
              actContext.GetSecondFeFunction()->GetFctId() == fctId2 &&
              actContext.IsSetCounterPart() ) {
            doAssemble = true;
            doTranspose = true;
          }

        // Check if fctIds of formsContext match
        if( doAssemble ) {

          try {

            // Get equation numbers
            actContext.MapEqns( it1, it2, eqnVec1, eqnVec2, id1, id2 );

            // Pass entity eqn-connectivity to algebraic system
            if( !doTranspose ) {
              algsys_-> SetElementPos( id1, eqnVec1,
                                       id2, eqnVec2,
                                       destMap,
                                       setCounterPart );
            } else {
              algsys_-> SetElementPos( id2, eqnVec2,
                                       id1, eqnVec1,
                                       destMap,
                                       setCounterPart );
            }

          } catch (Exception& e) {
            RETHROW_EXCEPTION(e, "Could not setup matrix graph for "
                              << "BiLinearForm '"
                              << actContext.GetIntegrator()->GetName() << "' on '"
                              << actContext.GetFirstEntities()->GetName()<< "'" );
          }

        } // if
        } // loop (integrators)
      } // loop over entities
    } // loop over entity-list pair
  }

  void Assemble::AssembleMatrices() {
    
    // check for static condensation:
    
    // If static condensation is enabled, we must always pass the 
    // complete system matrix per element to the algebraic system
    // in order to be able to invert the complete inner block.
    // This influences of course also the secondary matrix factors
    // and the assembly of complex entries. This has to be refactored
    // in this case. 
    // Until we have now common structure for both cases, we
    // excplicitly distinguish both cases in two different methods.
    
    if (algsys_->UseStaticCondensation() ) {
      AssembleMatrices_Cond();
    } else {
      AssembleMatrices_Std();
    }
    
  }
  
  void Assemble::AssembleMatrices_Std() {
    Matrix<Double> elemMatrix;
    Matrix<Complex> elemMatrixC;
    StdVector<Integer> eqnVec1, eqnVec2;
    FeFctIdType fctId1, fctId2;

    timer_->Start();

    // Reset for matrix update
    matrixUpdated_ = false;

    // Temporary: Check each time for non-linearities
    // On first Assembly, assemble all matrices for each BilinearForm
    CheckNonLinearities(isFirstTime_);

    // Init all matrices, which have to be reassembled
    std::map<FEMatrixType, bool>::iterator it;
    for( it = matReassemble_.begin(); it != matReassemble_.end(); it++ ) {
      if( it->second == true ) {
        LOG_DBG2(assemble) << "AssembleMatrices: init matrix " << it->first;
        algsys_->InitMatrix( matrixMap_[it->first] );
      }
    }

    // iterate over all entitylist-pairs and 
    BiLinContextListType::iterator listIt = biLinForms_.begin();
    for ( ; listIt != biLinForms_.end(); ++listIt) {
      StdVector<BiLinFormContext*> & forms = listIt->second;
      EntityList& firstEntities = *(listIt->first.first);
      EntityList& secondEntities = *(listIt->first.second);
      UInt size = firstEntities.GetSize();

      std::cout << "  - Calculating BiLinearForms on '" 
          << firstEntities.GetName() 
          << " (" << size << " elements)'\n";
      // Total work: numElement x numForms
      boost::progress_display progress( size*forms.GetSize() );


      // Loop over all entities
      EntityIterator it1 = firstEntities.GetIterator();
      EntityIterator it2 = secondEntities.GetIterator();
      for( it1.Begin(); !it1.IsEnd(); it1++, it2++ ) {
        LOG_DBG2(assemble) << "\telems are " << it1.GetElem()->elemNum 
            << " and " << it2.GetElem()->elemNum;

        // Loop over all bilinearforms
        for( UInt iForm = 0; iForm < forms.GetSize(); ++iForm ) {

          BiLinFormContext & actContext = *forms[iForm];

          LOG_DBG(assemble) << "AssembleMatrices for bilinform " << actContext.GetIntegrator()->GetName();
          LOG_DBG2(assemble) << "bilinform=" << actContext.ToString();

          // get matrix destinations
          FEMatrixType destMat = actContext.GetDestMat();
          FEMatrixType secDestMat = actContext.GetSecDestMat();

          // get secondary matrix factor string 
          Double secMatFac = actContext.EvalSecMatFac();

          // If assemble was already called and the current destination
          // matrix must not be reassembled -> continue with next iterator
          if( matReassemble_[destMat] == false ) {
            continue;
          }

          // Update flag
          matrixUpdated_ = true;

          BaseForm * form = actContext.GetIntegrator();

          try {
            ++progress;

            // Calc element matrix
            if ( form->IsComplex() ){
              form->CalcElementMatrix( elemMatrixC, it1, it2 );
            } else {
              form->CalcElementMatrix( elemMatrix, it1, it2 );
              if(actContext.IsSetNegate()== true){
                assert(!form->IsComplex());
                elemMatrix*= (-1.0);
              }
               }

   //            // info.xml logging in detailed logging case for the first element only
   //            if(i == 0 && progOpts->DoDetailedInfo())
   //            {
   //              PtrParamNode in = info->Get("PDE")->Get(actContext.GetFirstPde()->GetName())->Get(ParamNode::PROCESS)->Get("matrices");
   //              in = in->Get("tensor", ParamNode::APPEND);
   //              in->Get("form")->SetValue(form->GetName());
   //              in->Get("region")->SetValue(domain->GetGrid()->GetRegion().ToString(it1.GetElem()->regionId));
   //              in->Get("element")->SetValue(it1.GetElem()->elemNum);
   //              // it makes sense to add the analysis id here
   //              if(form->IsComplex())
   //                in->Get("matrix")->SetValue(elemMatrixC);
   //              else
   //              {
   //                in->Get("matrix")->SetValue(elemMatrix);
   //              }
   //            }

   //            LOG_DBG3(assemble) << "CalcElemMatrix " << i << " -> " << (form->IsComplex() ? elemMatrixC.ToString(0) : elemMatrix.ToString(0));
   //            LOG_DBG2(assemble) << "CalcElemMatrix " << i << " -> " << (form->IsComplex() ? elemMatrixC.ToString(1) : elemMatrix.ToString(1));

               // Map equation numbers
               actContext.MapEqns( it1, it2, eqnVec1, eqnVec2, fctId1, fctId2 );


   //#ifdef USE_SCRIPTING
   //            // Check, if current list is element list and if element matrix
   //            // should be printed
   //            if( actContext.GetFirstEntities()->GetType() == EntityList::ELEM_LIST ) {
   //              if( printElemNums_.count( it1.GetElem()->elemNum ) ) {
   //                StdVector<std::string> args;
   //                args.Push_back( form->GetName() );
   //                args.Push_back( it1.GetIdString() );
   //                args.Push_back( it2.GetIdString() );
   //                if( form->IsComplex() ) {
   //                  args.Push_back( elemMatrixC.ToString(1) );
   //                } else {
   //                  args.Push_back( elemMatrix.ToString(1) );
   //                }
   //                messenger->TriggerEvent( CFSMessenger::CFS_AssembleMat, args );
   //              }
   //            }
   //#endif

   //            assert((form->IsComplex() && 
   //                    eqnVec1.GetSize() == elemMatrixC.GetNumRows() && 
   //                    eqnVec2.GetSize() == elemMatrixC.GetNumCols()) || !form->IsComplex());
   //            assert((!form->IsComplex() && 
   //                    eqnVec1.GetSize() == elemMatrix.GetNumRows() && 
   //                    eqnVec2.GetSize() == elemMatrix.GetNumCols()) || form->IsComplex());

               // Pass element matrix to algebraic system (primary matrix)
               if ( form->IsComplex() )
                 InsertMatrix( destMat, actContext, elemMatrixC, eqnVec1, eqnVec2, fctId1, fctId2);
               else
                 InsertMatrix( destMat, actContext, elemMatrix, eqnVec1, eqnVec2, fctId1, fctId2);

               // if optimization provides Damping Parameters, we use them, and ignore everything else
               //double secMatFacOpt = 0.0;
               //          if(domain->HasErsatzMaterialDamping() && 
               //              domain->GetErsatzMaterial()->GetErsatzMaterialDampingParameterForIntegrator(it1.GetElem(), form, secMatFacOpt)){
               //            elemMatrix *= secMatFacOpt; // only in non-complex case, complex is not known in ParamMat
               //            InsertMatrix(DAMPING, actContext, elemMatrix, eqnVec1, eqnVec2, pdeId1, pdeId2);
               //          }else 
               if (secDestMat != NOTYPE ) { // Check for secondary matrix type
                 EXCEPTION("We do not want a second matrix factor");
                 Double dampFactor = 1.0;
                 //            if ( actContext.getPtDamplayer() != NULL ) {
                 //              actContext.getPtDamplayer()->CalcDampFactor(dampFactor, it1);
                 //              // std::cout << "   dampFactor: " <<  dampFactor << std::endl;
                 //            }

                 // damping with "exotic" complex material
                 if ( form->IsComplex() ) {
                   // complex damping
                   elemMatrixC *= secMatFac * dampFactor;

                   // Pass secondary matrix part to algebraic system
                   InsertMatrix(secDestMat, actContext, elemMatrixC, eqnVec1, eqnVec2, fctId1, fctId2);
                 }
                 // "standard" Rayleigh damping. Includes the standard SIMP optimization!
                 else {
                   // Rayleigh damping
                   elemMatrix *= secMatFac * dampFactor;
                   LOG_DBG3(assemble) << "AM: e1=" << it1.GetElem()->elemNum << " Rayleigh damping form=" << form->GetName() << " sMF=" << secMatFac << " df=" <<  dampFactor;
                   // Pass secondary matrix part to algebraic system
                   InsertMatrix(secDestMat, actContext, elemMatrix, eqnVec1, eqnVec2, fctId1, fctId2);
                 }
                 // optionally with do SIMP pamping (intermediate material has complex mass damping)
                 //            if(domain->GetErsatzMaterialPamping(it1.GetElem(), form, elemMatrix)) {
                 //              LOG_DBG3(assemble) << "AM: e1=" << it1.GetElem()->elemNum << " pamping form=" << form->GetName() << " add=" << elemMatrix.ToString();
                 //              InsertMatrix(secDestMat, actContext, elemMatrix, eqnVec1, eqnVec2, pdeId1, pdeId2);
                 //            }

               } // handle secDestMat != NOTYPE

               // increment iterators
             } catch (Exception& e) {
               RETHROW_EXCEPTION(e, "Could not calculate element matrix of "
                                 << "BiLinearForm '"
                                 << form->GetName() << "' on '"
                                 << actContext.GetFirstEntities()->GetName()<< "'" );
             }
           } // loop over bilinearforms
         } // loop over entities
       }// loop over entitylist pairs
       // Change flag
       isFirstTime_ = false;

       // We have assembled all matrices, we will not do so next time
       // except: CheckNonLinearities sets one of these, or Optimization does
       matReassemble_.clear();

       timer_->Stop();
    
  }
  
    void Assemble::AssembleMatrices_Cond() {

    Matrix<Double> elemMatrix;
    Matrix<Complex> elemMatrixC;
    StdVector<Integer> eqnVec1, eqnVec2;
    FeFctIdType fctId1, fctId2;

    timer_->Start();

    // Reset for matrix update
    matrixUpdated_ = false;

    // Temporary: Check each time for non-linearities
    // On first Assembly, assemble all matrices for each BilinearForm
    CheckNonLinearities(isFirstTime_);

    // Init all matrices, which have to be reassembled
    std::map<FEMatrixType, bool>::iterator it;
    for( it = matReassemble_.begin(); it != matReassemble_.end(); it++ ) {
      if( it->second == true ) {
        LOG_DBG2(assemble) << "AssembleMatrices: init matrix " << it->first;
        algsys_->InitMatrix( matrixMap_[it->first] );
      }
    }

    // temporary matrices 
    Matrix<Double> rElemMat;
    Matrix<Complex> cElemMat;
    
    // iterate over all entitylist-pairs and 
    BiLinContextListType::iterator listIt = biLinForms_.begin();
    for ( ; listIt != biLinForms_.end(); ++listIt) {
      StdVector<BiLinFormContext*> & forms = listIt->second;
      EntityList& firstEntities = *(listIt->first.first);
      EntityList& secondEntities = *(listIt->first.second);
      UInt size = firstEntities.GetSize();

      std::cout << "  - Calculating BiLinearForms on '" 
                    << firstEntities.GetName() 
                    << " (" << size << " elements)'\n";
      // Total work: numElement x numForms
      boost::progress_display progress( size*forms.GetSize() );
      
      
      // Loop over all entities
      EntityIterator it1 = firstEntities.GetIterator();
      EntityIterator it2 = secondEntities.GetIterator();
      for( it1.Begin(); !it1.IsEnd(); it1++, it2++ ) {
        LOG_DBG2(assemble) << "\telems are " << it1.GetElem()->elemNum 
                           << " and " << it2.GetElem()->elemNum;
        try {
        // Loop over all bilinearforms
        for( UInt iForm = 0; iForm < forms.GetSize(); ++iForm ) {

          BiLinFormContext & actContext = *forms[iForm];

          LOG_DBG(assemble) << "AssembleMatrices for bilinform " << actContext.GetIntegrator()->GetName();
          LOG_DBG2(assemble) << "bilinform=" << actContext.ToString();

          // get matrix destinations
          FEMatrixType destMat = actContext.GetDestMat();
//          FEMatrixType secDestMat = actContext.GetSecDestMat();

          // get secondary matrix factor string 
//          Double secMatFac = actContext.EvalSecMatFac();

          // If assemble was already called and the current destination
          // matrix must not be reassembled -> continue with next iterator
          if( matReassemble_[destMat] == false ) {
            continue;
          }

          // Update flag
          matrixUpdated_ = true;

          BaseForm * form = actContext.GetIntegrator();

         
    

            ++progress;

            // Calc element matrix
            if ( form->IsComplex() ){
              form->CalcElementMatrix( elemMatrixC, it1, it2 );
              cElemMat += elemMatrixC;
            } else {
              form->CalcElementMatrix( elemMatrix, it1, it2 );
              if(actContext.IsSetNegate()== true){
                assert(!form->IsComplex());
                elemMatrix*= (-1.0);
              }
              if(iForm == 0) {
                 rElemMat = elemMatrix;
              } else {
            
              rElemMat += elemMatrix;
              }
            }

//            // info.xml logging in detailed logging case for the first element only
//            if(i == 0 && progOpts->DoDetailedInfo())
//            {
//              PtrParamNode in = info->Get("PDE")->Get(actContext.GetFirstPde()->GetName())->Get(ParamNode::PROCESS)->Get("matrices");
//              in = in->Get("tensor", ParamNode::APPEND);
//              in->Get("form")->SetValue(form->GetName());
//              in->Get("region")->SetValue(domain->GetGrid()->GetRegion().ToString(it1.GetElem()->regionId));
//              in->Get("element")->SetValue(it1.GetElem()->elemNum);
//              // it makes sense to add the analysis id here
//              if(form->IsComplex())
//                in->Get("matrix")->SetValue(elemMatrixC);
//              else
//              {
//                in->Get("matrix")->SetValue(elemMatrix);
//              }
//            }

//            LOG_DBG3(assemble) << "CalcElemMatrix " << i << " -> " << (form->IsComplex() ? elemMatrixC.ToString(0) : elemMatrix.ToString(0));
//            LOG_DBG2(assemble) << "CalcElemMatrix " << i << " -> " << (form->IsComplex() ? elemMatrixC.ToString(1) : elemMatrix.ToString(1));

            // Map equation numbers
            actContext.MapEqns( it1, it2, eqnVec1, eqnVec2, fctId1, fctId2 );
        } // loop over bilinearforms    // increment iterators

//#ifdef USE_SCRIPTING
//            // Check, if current list is element list and if element matrix
//            // should be printed
//            if( actContext.GetFirstEntities()->GetType() == EntityList::ELEM_LIST ) {
//              if( printElemNums_.count( it1.GetElem()->elemNum ) ) {
//                StdVector<std::string> args;
//                args.Push_back( form->GetName() );
//                args.Push_back( it1.GetIdString() );
//                args.Push_back( it2.GetIdString() );
//                if( form->IsComplex() ) {
//                  args.Push_back( elemMatrixC.ToString(1) );
//                } else {
//                  args.Push_back( elemMatrix.ToString(1) );
//                }
//                messenger->TriggerEvent( CFSMessenger::CFS_AssembleMat, args );
//              }
//            }
//#endif

//            assert((form->IsComplex() && 
//                    eqnVec1.GetSize() == elemMatrixC.GetNumRows() && 
//                    eqnVec2.GetSize() == elemMatrixC.GetNumCols()) || !form->IsComplex());
//            assert((!form->IsComplex() && 
//                    eqnVec1.GetSize() == elemMatrix.GetNumRows() && 
//                    eqnVec2.GetSize() == elemMatrix.GetNumCols()) || form->IsComplex());

            // Pass element matrix to algebraic system (primary matrix)
//            if ( form->IsComplex() )
//              InsertMatrix( destMat, actContext, elemMatrixC, eqnVec1, eqnVec2, fctId1, fctId2);
//            else
        // HARD-CODED Section:
               // just take last bilinearform and ask everything from there
               BiLinFormContext &  actContext = *forms.Last();
               //BaseForm * form = actContext.GetIntegrator();
               FEMatrixType destMat = actContext.GetDestMat();
//               FEMatrixType secDestMat = actContext.GetSecDestMat();
//               Double secMatFac = actContext.EvalSecMatFac();

              InsertMatrix( destMat, actContext, rElemMat, eqnVec1, eqnVec2, fctId1, fctId2);

            // if optimization provides Damping Parameters, we use them, and ignore everything else
            //double secMatFacOpt = 0.0;
            //          if(domain->HasErsatzMaterialDamping() && 
            //              domain->GetErsatzMaterial()->GetErsatzMaterialDampingParameterForIntegrator(it1.GetElem(), form, secMatFacOpt)){
            //            elemMatrix *= secMatFacOpt; // only in non-complex case, complex is not known in ParamMat
            //            InsertMatrix(DAMPING, actContext, elemMatrix, eqnVec1, eqnVec2, pdeId1, pdeId2);
            //          }else 
//            if (secDestMat != NOTYPE ) { // Check for secondary matrix type
//              EXCEPTION("We do not want a second matrix factor");
//              Double dampFactor = 1.0;
//              //            if ( actContext.getPtDamplayer() != NULL ) {
//              //              actContext.getPtDamplayer()->CalcDampFactor(dampFactor, it1);
//              //              // std::cout << "   dampFactor: " <<  dampFactor << std::endl;
//              //            }
//
//              // damping with "exotic" complex material
//              if ( form->IsComplex() ) {
//                // complex damping
//                elemMatrixC *= secMatFac * dampFactor;
//
//                // Pass secondary matrix part to algebraic system
//                InsertMatrix(secDestMat, actContext, elemMatrixC, eqnVec1, eqnVec2, fctId1, fctId2);
//              }
//              // "standard" Rayleigh damping. Includes the standard SIMP optimization!
//              else {
//                // Rayleigh damping
//                elemMatrix *= secMatFac * dampFactor;
//                LOG_DBG3(assemble) << "AM: e1=" << it1.GetElem()->elemNum << " Rayleigh damping form=" << form->GetName() << " sMF=" << secMatFac << " df=" <<  dampFactor;
//                // Pass secondary matrix part to algebraic system
//                InsertMatrix(secDestMat, actContext, elemMatrix, eqnVec1, eqnVec2, fctId1, fctId2);
//              }
//              // optionally with do SIMP pamping (intermediate material has complex mass damping)
//              //            if(domain->GetErsatzMaterialPamping(it1.GetElem(), form, elemMatrix)) {
//              //              LOG_DBG3(assemble) << "AM: e1=" << it1.GetElem()->elemNum << " pamping form=" << form->GetName() << " add=" << elemMatrix.ToString();
//              //              InsertMatrix(secDestMat, actContext, elemMatrix, eqnVec1, eqnVec2, pdeId1, pdeId2);
//              //            }
//
//            } // handle secDestMat != NOTYPE

          
        } catch (Exception& e) {
          RETHROW_EXCEPTION(e, "Could not calculate element matrix of "
                            << "BiLinearForm"  );
        }
        
      } // loop over entities
    }// loop over entitylist pairs
    // Change flag
    isFirstTime_ = false;

    // We have assembled all matrices, we will not do so next time
    // except: CheckNonLinearities sets one of these, or Optimization does
    matReassemble_.clear();

    timer_->Stop();
  }

  void Assemble::CalcMinMaxStrain() {
    EXCEPTION("Not implemented");
//    // iterate over all descriptors
//    StdVector<BiLinFormContext*>::iterator formsIt;
//    for ( formsIt = biLinForms_->Begin(); 
//          formsIt != biLinForms_->End(); formsIt++ ) {
//      
//      // get integrator
//      BiLinFormContext & actContext = **formsIt;
//      
//      BaseForm * form = actContext.GetIntegrator();
//      
//      form->ResetMinMaxStrain();
//        
//      // get entity iterators
//      EntityIterator  it1 = actContext.GetFirstEntities()->GetIterator();
//      EntityIterator  it2 = actContext.GetSecondEntities()->GetIterator();
//      UInt size = actContext.GetFirstEntities()->GetSize();
//      it1.Begin();
//      it2.Begin();
//        
//      // iterate over all entities
//      for ( UInt i=0; i<size; i++ ) {
//          
//        form->CalcMinMaxStrain(it1, it2 );
//        
//        // increment iterators
//        it1++;
//        it2++;
//      }
//    }
  }

  
  void Assemble::AssembleLinRHS(AdjointParameters* adjointParams)
  {
    Optimization* opt = domain->GetOptimization();
    if(adjointParams != NULL && opt != NULL){
      opt->SetAdjointRhs(adjointParams);
    }else{
      AssembleRHSLinForms(false );
      AssembleRHSLoads();
      if(opt != NULL){
        opt->RhsCalculated(adjointParams);
      }
    }
  }

  void Assemble::AssembleNonLinRHS(AdjointParameters* adjointParams) {
    Optimization* opt = domain->GetOptimization();
    if(adjointParams != NULL && opt != NULL){
      domain->GetOptimization()->SetAdjointRhs(adjointParams);
    }else{
      AssembleRHSLinForms(true );
      if(opt != NULL){
        opt->RhsCalculated(adjointParams);
      }
    }
  }

  void Assemble::AssembleRHSLinForms(bool nonLin ) {

    StdVector<Integer> eqnVec;
    FeFctIdType fctId;
    StdVector<LinearFormContext*>::iterator formsIt;

    // iterate over all descriptors
    for ( formsIt = linForms_->Begin();
          formsIt != linForms_->End();
          formsIt++ ) {

      // get integrator
      LinearFormContext & actContext = **formsIt;

      // Check, if lin/non-lin type of Context matches parameter nonLin
      if( actContext.IsNonLin() != nonLin ) {
        continue;
      }

      LinearForm * form = actContext.GetIntegrator();

      try {

        // get entity iterator
        
        EntityIterator  entIt = actContext.GetEntities()->GetIterator();
        UInt size = actContext.GetEntities()->GetSize();
        
        std::cout << "  - Calculating '" << form->GetName() << "' on '" 
            << actContext.GetEntities()->GetName() 
            << " (" << size << " elements)'\n";
        
        boost::progress_display progress( size );

        if ( analysisType_ == BasePDE::HARMONIC ) {

          Vector<Complex> elemVec;
          for ( entIt.Begin(); !entIt.IsEnd(); entIt++ ) {
            ++progress;

            // Calculate complex valued element vector
            form->CalcElemVector( elemVec, entIt );

            // Map equation numbers
            actContext.MapEqns( entIt, eqnVec, fctId );

#ifdef USE_SCRIPTING
          // Check, if current list is element list and if element matrix
          // should be printed
          if( actContext.GetEntities()->GetType() == EntityList::ELEM_LIST ) {
            if( printElemNums_.count( entIt.GetElem()->elemNum ) ) {
              // Trigger event for scripting
              StdVector<std::string> args;
              args.Push_back( form->GetName() );
              args.Push_back( entIt.GetIdString() );
              args.Push_back( elemVec.ToString(' ') );
              messenger->TriggerEvent( CFSMessenger::CFS_AssembleRhs, args );
            }
          }
#endif

#ifdef USE_SCRIPTING
          // Check, if current list is element list and if element matrix
          // should be printed
          if( actContext.GetEntities()->GetType() == EntityList::ELEM_LIST ) {
            if( printElemNums_.count( entIt.GetElem()->elemNum ) ) {
              // Trigger event for scripting
              StdVector<std::string> args;
              args.Push_back( form->GetName() );
              args.Push_back( entIt.GetIdString() );
              args.Push_back( elemVec.ToString(' ') );
              messenger->TriggerEvent( CFSMessenger::CFS_AssembleRhs, args );
            }
          }
#endif


            assert(!elemVec.ContainsNaN() && !elemVec.ContainsInf());

            algsys_-> SetElementRHS( elemVec,
                                     fctId, eqnVec );
          }

        } else {

          // That should be STATIC, TRANSIENT or EIGENFREQUENCY

          Vector<Double> elemVec;
          // iterate over all entities
          for ( entIt.Begin(); !entIt.IsEnd(); entIt++ ) {
            
            ++progress;

            // Calculate real valued element vector
            form->CalcElemVector( elemVec, entIt );

            // Map equation numbers
            actContext.MapEqns( entIt, eqnVec, fctId );
#ifdef USE_SCRIPTING
          // Check, if current list is element list and if element matrix
          // should be printed
          if( actContext.GetEntities()->GetType() == EntityList::ELEM_LIST ) {
            if( printElemNums_.count( entIt.GetElem()->elemNum ) ) {
              // Trigger event for scripting
              StdVector<std::string> args;
              args.Push_back( form->GetName() );
              args.Push_back( entIt.GetIdString() );
              args.Push_back( elemVec.ToString() );
              messenger->TriggerEvent( CFSMessenger::CFS_AssembleRhs, args );
            }
          }
#endif
#ifdef USE_SCRIPTING
          // Check, if current list is element list and if element matrix
          // should be printed
          if( actContext.GetEntities()->GetType() == EntityList::ELEM_LIST ) {
            if( printElemNums_.count( entIt.GetElem()->elemNum ) ) {
              // Trigger event for scripting
              StdVector<std::string> args;
              args.Push_back( form->GetName() );
              args.Push_back( entIt.GetIdString() );
              args.Push_back( elemVec.ToString() );
              messenger->TriggerEvent( CFSMessenger::CFS_AssembleRhs, args );
            }
          }
#endif
            
            assert(!elemVec.ContainsNaN() && !elemVec.ContainsInf());

            // Pass element vector to algebraic system
            algsys_-> SetElementRHS( elemVec,
                                     fctId, eqnVec );
          }

        }
      } catch (Exception& e) {
        RETHROW_EXCEPTION(e, "Could not calculate RHS vector for LinearForm '"
                          << form->GetName() << "' on '"
                          << actContext.GetEntities()->GetName() << "'" );
      }
    }
  }

  void Assemble::AssembleRHSLoads() {
    /*
    Vector<Double> globCoord;
    Double phase = 0.0;
    Double val = 0.0;
    Complex complexValue( 0.0, 0.0 );
    
    // get global math parser and pointer to grid
    MathParser * parser = domain->GetMathParser();
    CoordSystem * coosy = domain->GetCoordSystem();
    Grid * ptGrid = domain->GetGrid();

    // iterate over all load node-lists
    for( UInt iLoad=0, numLoads=loads_.GetSize(); iLoad<numLoads; ++iLoad )
    {

      // get current load and its equation map
      LoadBc & actLoad = *loads_[iLoad];
     // EqnMap & eqnMap = *(actLoad.eqnMap);

      // create temporary set, in which we store the equation
      // numbers which are already set, to that we avoid setting multiple
      // times the RHS of a node, if the entity list is defined on
      // surface elements
      std::set<Integer> usedEqns;

      // Provide result info to ResultCache, so we can use the input function
      // in math parser expressions
      SolutionType solType = NO_SOLUTION_TYPE;
      switch ( actLoad.result->resultType )
      {
        case MECH_DISPLACEMENT:
          solType = MECH_RHS_LOAD;
          break;
        case ELEC_POTENTIAL:
          solType = ELEC_RHS_LOAD;
          break;
        case ACOU_POTENTIAL:
        case ACOU_PRESSURE:
          solType = ACOU_RHS_LOAD;
          break;
        case MAG_POTENTIAL:
        case MAG_SCALAR_POTENTIAL:
          solType = MAG_RHS_LOAD;
          break;
        case HEAT_TEMPERATURE:
          solType = HEAT_RHS_LOAD;
          break;
        default:
          EXCEPTION("Cannot determine load corresponding to result '"
              << SolutionTypeEnum.ToString(actLoad.result->resultType)
              << "'");
      }
      ResultCache::SetInfo( ResultCache::OUT_REAL, actLoad.dof,
                            actLoad.entities->GetName(), solType );
      
      try {
        // Obtain iterator
        EntityIterator it = actLoad.entities->GetIterator();

        for( it.Begin(); !it.IsEnd(); it++ ) {

          // If iterator points to a node, pass the current coordinate
          // to the parser
          if( it.GetType() == EntityList::NODE_LIST ) {

            // Get coordinates of node
            ptGrid->GetNodeCoordinate( globCoord, it.GetNode() );
            parser->SetCoordinates( mHandle_, *coosy, globCoord );
          }

          // Provide global node number to ResultCache, so we can use
          // the input function in math parser expressions
          ResultCache::SetIndex(it.GetNode());
          
          // Evaluate load value
          ResultCache::SetOutputType( (analysisType_ == BasePDE::HARMONIC ? 
              ResultCache::OUT_AMPL : ResultCache::OUT_REAL) );
          parser->SetExpr( mHandle_, actLoad.value );
          val = parser->Eval(mHandle_ );

#ifndef NDEBUG
          if ( std::isnan(val) || std::isinf(val) )
            EXCEPTION("Trying to assemble nan/inf in AssembleRHSLoads!");
#endif

          // for a harmonic simulation: evaluate phase
          if ( analysisType_ == BasePDE::HARMONIC )
          {
            ResultCache::SetOutputType(ResultCache::OUT_PHASE);
            parser->SetExpr( mHandle_, actLoad.phase  );
            phase = parser->Eval( mHandle_ );

#ifndef NDEBUG
            if ( std::isnan(phase) || std::isinf(phase) )
              EXCEPTION("Trying to assemble nan/inf in AssembleRHSLoads!");
#endif
            complexValue = Complex( val * cos( phase / 180 * PI ),
                                    val * sin( phase / 180 * PI ) );
          }

          // Obtain equation number(s)
          StdVector<Integer> eqns;
          feSpace_.GetEqns( eqns, it, actLoad.dof );

          for( UInt iEqn=0, numEqns=eqns.GetSize(); iEqn<numEqns; ++iEqn )
          {
            // check, if RHS of current eqn was already set
            if( usedEqns.find( eqns[iEqn] ) == usedEqns.end() ) {
              usedEqns.insert( eqns[iEqn] );
              if (analysisType_ == BasePDE::HARMONIC) {
                parser->SetExpr( mHandle_, actLoad.phase  );
                phase = parser->Eval( mHandle_ );

                assert(!std::isnan(val) && !std::isinf(val) && !std::isnan(phase) && !std::isinf(phase));
                Complex complexValue( val * cos( phase / 180 * PI ),
                                      val * sin( phase / 180 * PI ) );

                algsys_->SetNodeRHS(complexValue, eqnMap.GetPdeId(),
                                    eqns[iEqn] );
              } else {
                assert(!std::isnan(val) && !std::isinf(val));
                algsys_->SetNodeRHS(val, eqnMap.GetPdeId(), eqns[iEqn] );
              } // analysis
            } // usedEqn
          } // for loop

        } // nodes
      } catch (Exception& e) {
        RETHROW_EXCEPTION(e, "Could not set RHS load with value '"
                          << actLoad.value << "' on '"
                          << actLoad.entities->GetName() << "'" );
      }
    } // loads
*/
  }

  void Assemble::ToInfo(PtrParamNode in)
  {
    WARN("Adjust to new implementation");
//    PtrParamNode list = in->Get("matrixBiLinearForms");
//
//    // iterate over all descriptors
//    StdVector<BiLinFormContext*>::iterator it;
//    for ( it = biLinForms_->Begin(); it != biLinForms_->End(); it++ )
//    {
//      // get integrator
//      BiLinFormContext & context = **it;
//
//      PtrParamNode form = list->Get("bilinearForm", ParamNode::APPEND);
//      // integrator name
//      form->Get("integrator")->SetValue(context.GetIntegrator()->GetName());
//
//      // region name of entity list
//      std::string regionName;
//      if( context.GetFirstEntities()->GetType() == EntityList::ELEM_LIST )
//      {
//        shared_ptr<ElemList> list =
//          boost::dynamic_pointer_cast<ElemList,EntityList>
//          (context.GetFirstEntities());
//        regionName = list->GetName();
//      }
//      else if ( context.GetFirstEntities()->GetType() == EntityList::SURF_ELEM_LIST )
//      {
//        shared_ptr<SurfElemList> list =
//          boost::dynamic_pointer_cast<SurfElemList,EntityList>
//          (context.GetFirstEntities());
//        regionName = list->GetName();
//      }
//      form->Get("region")->SetValue(regionName);
//
//      // add information about row / column coordinate
//      PtrParamNode row = form->Get("row", ParamNode::APPEND);
//      PtrParamNode col = form->Get("column", ParamNode::APPEND);
//      
//     // associated PDEs
//      row->Get("pde")->SetValue(context.GetFirstPde()->GetName());
//      col->Get("pde")->SetValue(context.GetSecondPde()->GetName());
//
//      // associated result types
//      std::string tmp;
//      tmp = SolutionTypeEnum.ToString(context.GetFirstResultInfo()->resultType);
//      row->Get("result")->SetValue(tmp);
//      tmp = SolutionTypeEnum.ToString(context.GetSecondResultInfo()->resultType);
//      col->Get("result")->SetValue(tmp);
//      
//      // matrix destination
//      PtrParamNode dest = form->Get("destination", ParamNode::APPEND);
//      
//      // original destination matrix
//      Enum2String(context.GetDestMat(), tmp );
//      dest->Get("feMatrix")->SetValue(tmp);
//            
//      // mapped destination matrix
//      Enum2String(matrixMap_[context.GetDestMat()], tmp );
//      dest->Get("feMatrixMapped")->SetValue(tmp);
//      
//      // secondary destination matrix and factor
//      Enum2String(context.GetSecDestMat(), tmp );
//      dest->Get("feSecondMatrix")->SetValue(tmp);
//      dest->Get("feSecondMatrixFac")->SetValue(context.GetSecMatFac());
//      
//      // additional attributes
//      PtrParamNode attr = form->Get("attributes", ParamNode::APPEND);
//      
//      // entry Type (real / imag)
//      tmp = Global::complexPart.ToString(context.GetEntryType());
//      attr->Get("entryType")->SetValue( tmp );
//      
//      // flag setcounterpart
//      tmp = context.IsSetCounterPart() ? "yes" : "no";
//      attr->Get("counterPart")->SetValue( tmp );
//      
//      // issymmetric
//      tmp = context.GetIntegrator()->IsSymmetric() ? "yes" : "no";
//      attr->Get("symmetric")->SetValue( tmp );
//      
//      // isSolDependent
//      tmp = context.GetIntegrator()->IsSolDependent() ? "yes" : "no";
//      attr->Get("solutionDependent")->SetValue( tmp );
//      
//      // updated geometry
//      tmp = context.GetIntegrator()->IsCoordUpdate() ? "yes" : "no";
//      attr->Get("updatedGeo")->SetValue( tmp );
//      
//    }
//
//    list = in->Get("rhsLinearForms");
//
//    // iterate over all descriptors
//    StdVector<LinearFormContext*>::iterator linIt;
//    for (linIt = linForms_->Begin(); linIt != linForms_->End(); linIt++)
//    {
//      PtrParamNode form = list->Get("linearForm", ParamNode::APPEND);
//
//      // get integrator
//      LinearFormContext & context = **linIt;
//
//      form->Get("integrator")->SetValue(context.GetIntegrator()->GetName());
//
//      // region name of entity list
//      std::string regionName;
//      if( context.GetEntities()->GetType() == EntityList::ELEM_LIST )
//      {
//        shared_ptr<ElemList> list =
//          boost::dynamic_pointer_cast<ElemList,EntityList>
//          (context.GetEntities());
//        regionName = list->GetName();
//      }
//      else if ( context.GetEntities()->GetType() == EntityList::SURF_ELEM_LIST )
//      {
//        shared_ptr<SurfElemList> list =
//          boost::dynamic_pointer_cast<SurfElemList,EntityList>
//          (context.GetEntities());
//        regionName = list->GetName();
//      }
//
//      form->Get("region")->SetValue(regionName);
//      
//
//      // add information about row / column coordinate
//      PtrParamNode row = form->Get("row", ParamNode::APPEND);
//
//      // associated PDEs
//      row->Get("pde")->SetValue(context.GetPde()->GetName());
//
//      // associated result types
//      std::string tmp;
//      tmp = SolutionTypeEnum.ToString(context.GetResultInfo()->resultType);
//      row->Get("result")->SetValue(tmp);
//
//    }
  }


  void Assemble::CheckNonLinearities(bool setall) {
    // iterate over all bilinearforms
    std::set<BiLinFormContext*>::iterator it;

    for( it = allBiLinForms_.begin(); it != allBiLinForms_.end(); it++ )
    {
      BiLinFormContext & actContext = **it;

      if(actContext.IsNonLin() || analysisType_ == BasePDE::HARMONIC || setall)
      {
        matReassemble_[actContext.GetDestMat()] = true;
        if ( actContext.GetSecDestMat() != NOTYPE )
          matReassemble_[actContext.GetSecDestMat()] = true;
      }
    }
  }

  void Assemble::Matrix2Harmonic(Matrix<Complex>& harmMat,
                                 Matrix<Double>& origMat,
                                 FEMatrixType matrixType,
                                 Global::ComplexPart entryType,
                                 Double omega ) {

    Double factor = 0.0;
    UInt derivOrder = 0;
    // first determine factor of due to time derivative: d/dt = jOmega
    // 0: stiffness
    // 1: damping / mass
    // 2: mass
    switch( matrixType ) {
      case STIFFNESS:
        derivOrder = 0;
        factor = 1;
        break;
      case DAMPING:
        derivOrder = 1;
        factor = omega;
        break;
      case MASS:
        if( maxTimeDerivOrder_ == 2 ) {
          derivOrder = 2;
          factor = -omega*omega;
        } else {
          derivOrder = 1;
          factor = omega;
        }
        break;
      default:
        EXCEPTION("No default conversion from double entries to matrix type"
                  << matrixType << "known" );
    }
    
    // determine if destination is real / imaginary part
    // if time derivative is 0 or 2, a real part stay a real part
    // and a imaginary one a imaginary one.
    // Only if the derivative oder = 1, we have to interchange them
    Global::ComplexPart destType;
    if( derivOrder == 1 ) {
      destType = (entryType == Global::REAL) ? Global::IMAG : Global::REAL;
    } else {
      destType = entryType;
    } 
    
    harmMat.Resize( origMat.GetNumRows(), origMat.GetNumCols() );
    harmMat.SetPart( destType, origMat );
    harmMat *= factor;
  }

  void Assemble::Matrix2Harmonic(Matrix<Complex>& harmMat,
                                 Matrix<Complex>& origMat,
                                 FEMatrixType matrixType,
                                 Global::ComplexPart entryType,
                                 Double omega ) {
                                 
    Complex factor(0.0, 0.0);
    // first determine factor of due to time derivative: d/dt = jOmega
    // 0: stiffness
    // 1: damping / mass
    // 2: mass
    switch( matrixType ) {
      case STIFFNESS:
        factor = Complex(1.0, 0.0);
        break;
      case DAMPING:
        factor = Complex(0.0, omega);
        break;
      case MASS:
        if( maxTimeDerivOrder_ == 2 ) {
          factor = Complex(-omega*omega, 0.0);
        } else {
          factor = Complex(0.0, omega);
        }
        break;
      default:
        EXCEPTION("No default conversion from double entries to matrix type"
                  << matrixType << "known" );
    }
    harmMat = origMat * factor;

  }


  void Assemble::CreateMatrixMap()
  {

    // Dependent on the type of analysis, only certain matrix types
    // (SYSTEM, STIFFNESS, MASS, DAMPING, CONVECTION, ...) are present.
    switch(analysisType_)
    {
    case BasePDE::STATIC:
      matrixMap_[SYSTEM]    = SYSTEM;
      matrixMap_[STIFFNESS] = SYSTEM;
      matrixMap_[DAMPING]   = NOTYPE;
      matrixMap_[MASS]      = NOTYPE;
      matrixMap_[AUXILIARY] = NOTYPE;
      break;

    case BasePDE::TRANSIENT:
      matrixMap_[SYSTEM]    = SYSTEM;
      matrixMap_[STIFFNESS] = STIFFNESS;
      matrixMap_[DAMPING]   = DAMPING;
      matrixMap_[MASS]      = MASS;
      matrixMap_[AUXILIARY] = AUXILIARY;
      break;

    case BasePDE::HARMONIC:
      matrixMap_[SYSTEM]    = SYSTEM;
      matrixMap_[STIFFNESS] = SYSTEM;
      matrixMap_[DAMPING]   = SYSTEM;
      matrixMap_[MASS]      = SYSTEM;
      matrixMap_[AUXILIARY] = AUXILIARY; // optimization for radiation needs this
      break;

    case BasePDE::EIGENFREQUENCY:
      matrixMap_[SYSTEM]    = NOTYPE;
      matrixMap_[STIFFNESS] = STIFFNESS;
      matrixMap_[DAMPING]   = DAMPING;
      matrixMap_[MASS]      = MASS;
      matrixMap_[AUXILIARY] = NOTYPE;
      break;

    default: EXCEPTION("Analysistype '" << BasePDE::analysisType.ToString(analysisType_) << "' not known!");
    }
  }



  bool Assemble::IsFEMatSymmetric( FeFctIdType fctId1, FeFctIdType fctId2,
                                   FEMatrixType matType  ) {

    REFACTOR;
    return true;
//    // Run over all bilinearform contexts
//    std::map<FEMatrixType, bool> isSymmetric;
//    StdVector<BiLinFormContext*>::iterator it;
//
//    // Assume at the beginning that all matrices are symmetric
//    isSymmetric[SYSTEM] = true;
//    isSymmetric[MASS] = true;
//    isSymmetric[STIFFNESS] = true;
//    isSymmetric[DAMPING] = true;
//    isSymmetric[AUXILIARY] = true;
//    
//    // set collecting all result types
//    
//    // set collecting all diagonal-positions
//    std::set<shared_ptr<ResultInfo> > allResults, diagResults;
//
//
//    // iterate over all bilinear forms
//    for( it = biLinForms_->Begin(); it != biLinForms_->End(); it++ ) {
//
//      BiLinFormContext & actCt = (**it);
//
//      // Check, where bilinearform gets assembled to diagonal block
//      if( (actCt.GetFirstPde() == actCt.GetSecondPde() )
//          && (actCt.GetFirstResultInfo() == actCt.GetSecondResultInfo() )
//          && (actCt.GetFirstEntities() == actCt.GetSecondEntities() ) ) {
//
//        // Bilinearform gets assembled to main diagonal.
//        // If bilinearform is non-symmetric, so is the related FE-matrix
//        if( !actCt.GetIntegrator()->IsSymmetric()  ) {
//          FEMatrixType mappedDest = matrixMap_[actCt.GetDestMat()];
//          isSymmetric[mappedDest] = false;
//          isSymmetric[SYSTEM] = false;
//        }
//      } else {
//
//        // BiLinearform gets assembled to off-diagonal block.
//
//        // If the bilinearorm is also assembled to the transposed block
//        // we assume that the matrix still remains symmetric.
//        // Otherwise we assume, that we need a non-symmetric matrix.
//        if( !actCt.IsSetCounterPart() ) {
//          FEMatrixType mappedDest = matrixMap_[actCt.GetDestMat()];
//          isSymmetric[mappedDest] = false;
//          isSymmetric[SYSTEM] = false;
//        }
//      }
//    }
//
//    // return flag for matrix of interest
//    return isSymmetric[feType];
  }


  void Assemble::InsertMatrix( FEMatrixType dest, BiLinFormContext& context,
                               Matrix<Double>& elemMat,
                               StdVector<Integer>& eqnVec1,
                               StdVector<Integer>& eqnVec2,
                               FeFctIdType fctId1, FeFctIdType fctId2)
  {
    // map original matrix destination to analysis-dependent one
    FEMatrixType mappedDest = matrixMap_[dest];

    Matrix<Complex> harmMat;
    // if destination matrix is NOTYPE -> leave
    if( mappedDest == NOTYPE )
      return;

    assert(!elemMat.ContainsNaN() && !elemMat.ContainsInf());

    if( analysisType_ == BasePDE::TRANSIENT
        || analysisType_ == BasePDE::STATIC
        || analysisType_ == BasePDE::EIGENFREQUENCY) {
      algsys_->SetElementMatrix( mappedDest, elemMat,
                                 fctId1, eqnVec1,
                                 fctId2, eqnVec2,
                                 context.IsSetCounterPart() );

    } else {
      assert(analysisType_ == BasePDE::HARMONIC);

      Double freq = context.GetFirstPde()->GetSolveStep()->GetActFreq();
      Double omega = freq * 2 * PI;

      Matrix2Harmonic( harmMat, elemMat, dest, context.GetEntryType(), omega );
      algsys_->SetElementMatrix( mappedDest, harmMat,
                                    fctId1, eqnVec1,
                                    fctId2, eqnVec2,
                                    context.IsSetCounterPart() );
    }

    //LOG_DBG3(assemble) << "InsertMatrix dest=" << dest << " mappedDest=" << mappedDest << " data=["
    //                   << elemMat  << "]\neqnVec1=" << eqnVec1.ToString() <<  "\neqnVec2=" << eqnVec2.ToString() << std::endl;
   
  }

  void Assemble::InsertMatrix( FEMatrixType dest, BiLinFormContext& context,
                               Matrix<Complex>& elemMat,
                               StdVector<Integer>& eqnVec1,
                               StdVector<Integer>& eqnVec2,
                               FeFctIdType fctId1, FeFctIdType fctId2) {
    Matrix<Complex> harmMat;

    // map original matrix destination to analysis-dependent one
    FEMatrixType mappedDest = matrixMap_[dest];

    assert(mappedDest != NOTYPE);
    assert(analysisType_ == BasePDE::HARMONIC);

    assert(!elemMat.ContainsNaN() && !elemMat.ContainsInf());

    Double freq = context.GetFirstPde()->GetSolveStep()->GetActFreq();
    Double omega = freq * 2 * PI;

    Matrix2Harmonic( harmMat, elemMat, dest, context.GetEntryType(), omega );

    algsys_->SetElementMatrix( mappedDest, harmMat, fctId1, eqnVec1,
                               fctId2, eqnVec2, context.IsSetCounterPart() );
  }


  void Assemble::Wrap_AddPrintElemNum( ) {
    SCRIPT_GET( StdVector<UInt>, elemNums );
    printElemNums_.insert( elemNums.Begin(), elemNums.End() );
    StdVector<UInt>::iterator it = elemNums.Begin();

  }

  void Assemble::Wrap_PrintAllElems() {
    UInt numElems = domain->GetGrid()->GetNumElems();
    for( UInt i = 0; i < numElems; i++ ) {
      printElemNums_.insert( i+ 1);
    }
  } 

  void Assemble::RegisterFunctions() {
    typedef FctPointer<Assemble> FCPT;
    StdVector<ArgList> a;
    StdVector<FCPT*> pt;
    StdVector<std::string> name;

    // --- AddPrintElemNum ---
    a.Push_back();
    a.Last().RegisterParam("elemNums", ArgList::STDVEC_UINT );
    pt.Push_back( new FCPT( this, &Assemble::Wrap_AddPrintElemNum) );
    name.Push_back( "printElemMatVec" );
    
    // --- PrintAllElemNum ---
    a.Push_back();
    pt.Push_back( new FCPT( this, &Assemble::Wrap_PrintAllElems) );
    name.Push_back( "printAllElemMatVec" );
    

    // Now register all functions with scripting
    for (UInt i = 0; i < pt.GetSize(); i++ ) {
      Script_RegisterFct(name[i], pt[i], a[i] );
    }
  }
}
