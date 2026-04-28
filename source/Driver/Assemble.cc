#include "Assemble.hh"

#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <boost/lexical_cast.hpp>
#include <boost/timer/progress_display.hpp>

#include "Domain/Domain.hh"
#include "Domain/ElemMapping/EntityLists.hh"
#include "Domain/Mesh/Grid.hh"

#include "Driver/SolveSteps/StdSolveStep.hh"
#include "Driver/HarmonicDriver.hh"

#include "DataInOut/ProgramOptions.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"

#include "PDE/SinglePDE.hh"

#include "OLAS/algsys/AlgebraicSys.hh"
#include "MatVec/SBM_Matrix.hh"

#include "Utils/Timer.hh"
#include "Utils/mathParser/mathParser.hh"

#include "Forms/BiLinForms/BiLinearForm.hh"
#include "Forms/BiLinForms/BDBInt.hh"
#include "Forms/Operators/BaseBOperator.hh"
#include "Forms/LinForms/LinearForm.hh"
#include "Forms/LinForms/SingleEntryInt.hh"

#include "Optimization/Design/DesignElement.hh"
#include "Optimization/Design/DesignSpace.hh"

#include "Utils/PythonKernel.hh"

namespace CoupledField
{
  // declare logging stream
  DEFINE_LOG(assemble, "assemble")


  Assemble::Assemble( AlgebraicSys* algsys,
                      BasePDE::AnalysisType analysis,
                      MathParser* mp,
                      PtrParamNode infoNode)
  {

    // init general params
    algsys_ = algsys;
    analysisType_ = analysis;
    mp_ = mp;
    isFirstTime_ = true;
    matrixUpdated_ = false;
    printProgressBar_ = false;
    info_ = infoNode;
    skipElemAssembly_=false;
    timer_ = shared_ptr<Timer>(new Timer());

    // Calculate matrix map from general matrix types to analysis
    // specific ones
    CreateMatrixMap();

    // Set expression for omega
    mHandle_ = mp->GetNewHandle();
    mp->SetExpr(mHandle_, "2*pi*f");
    
    // the timer object is used in every AssembleMatrices() call
    info_->Get("analysis")->Get(ParamNode::SUMMARY)->Get("assemble/timer")->SetValue(timer_);

    // Sub timer
    matrixTimer_ = info_->Get("analysis")->Get(ParamNode::SUMMARY)->Get("assemble")->Get("matrices/timer")->AsTimer(timer_);
    rhsTimer_ = info_->Get("analysis")->Get(ParamNode::SUMMARY)->Get("assemble")->Get("rhs/timer")->AsTimer(timer_);
  }

  Assemble::~Assemble() {
    for(BiLinFormContext* blfctxt : allBiLinForms_)
      delete blfctxt;

    for(LinearFormContext* lfctxt: allLinForms_)
      delete lfctxt;

    mp_->ReleaseHandle(mHandle_);
  }

  void Assemble::SetAlgSys(AlgebraicSys * algsys)  {
    algsys_ = algsys;
  }

  void Assemble::SkipElemAssembly(){
    skipElemAssembly_=true;
  }

  void Assemble::ResetMatrixReassembly() {
    // Initialize reassemble-map
    // Note: in the beginning all matrices have to be "re-"assembled
    // which means the destination matrices have to be cleared and
    // initialized. After the first "AssmebleMatrices" call the
    // correct values for reassembling are determined
    matReassemble_[STIFFNESS] = true;
    matReassemble_[GEOMETRIC_STIFFNESS] = true;
    matReassemble_[DAMPING] = true;
    matReassemble_[DAMPING_AUX] = true;
    matReassemble_[MASS] = true;
    matReassemble_[STIFFNESS_UPDATE] = true;
    matReassemble_[DAMPING_UPDATE] = true;
    matReassemble_[MASS_UPDATE] = true;
    matReassemble_[AUXILIARY] = true;

    // reset also flag for "firstTime"
    isFirstTime_ = true;
    matrixUpdated_ = true;
  }

  void Assemble::SetEqnCustomMap( const std::map<Integer, Integer>& eqnMap,
                                  const std::map<FeFctIdType, FeFctIdType>& fctIdMap ) {
    customEqnMap_ = eqnMap;
    customFctIdMap_ = fctIdMap;
  }


  BiLinFormContext* Assemble::GetBiLinForm(const std::string& integrator, RegionIdType regionId, SinglePDE* pde1, SinglePDE* pde2, bool silent)
  {
    //std::cout << "pde1= " << pde1 << "pde2= " << pde2 << "silent= " << silent << std::endl;
    for(std::set<BiLinFormContext*>::iterator it = allBiLinForms_.begin(); it != allBiLinForms_.end(); ++it )
    {
      BiLinFormContext& context = **it;

      // makes only a starts-with comparison!
      if(context.GetIntegrator()->GetName().rfind(integrator, 0) == string::npos)
        continue;

      if(context.GetFirstEntities()->GetRegion() != regionId)
        continue;

      if(pde1 != NULL && pde1 != context.GetFirstPde())
        continue;

      if(pde2 != NULL && pde2 != context.GetSecondPde())
        continue;

      return &context;
    }

    if(!silent)
      EXCEPTION("cannot find integrator " << integrator << " in region " << domain->GetGrid()->GetRegion().ToString(regionId)
                << " for pde1=" << (pde1 != NULL ? pde1->GetName() : "null") << " and pde2=" << (pde2 != NULL ? pde2->GetName() : "null"));

    return NULL;
  }

  LinearFormContext* Assemble::GetLinForm(const std::string& integrator, RegionIdType regionId, StdPDE* pde, bool silent)
  {
    for(std::set<LinearFormContext*>::iterator it = allLinForms_.begin(); it != allLinForms_.end(); ++it )
    {
      LinearFormContext& context = **it;

      // makes only a starts-with comparison!
      if(context.GetIntegrator()->GetName().rfind(integrator, 0) == string::npos)
        continue;

      if(pde != NULL && pde != context.GetPde())
        continue;

      return &context;
    }

    if(!silent)
      EXCEPTION("cannot find integrator " << integrator << " in region " << domain->GetGrid()->GetRegion().ToString(regionId)
                << " for pde=" << (pde != NULL ? pde->GetName() : "null"));

    return NULL;
  }

  bool Assemble::UseRegion(RegionIdType reg)
  {
    for (BiLinContextListType::iterator listIt = biLinForms_.begin(); listIt != biLinForms_.end(); ++listIt) {
      if(listIt->first.first->GetRegion() == reg)
        return true;
    }

     return false;
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


    // If the datatype of the bilinearformcontext is "COMPLEX"
    // we have to ensure that we are in an harmonic case.
    // Otherwise we issue an error
    if( (biLinContext->GetEntryType() == Global::IMAG ||
         biLinContext->GetEntryType() == Global::COMPLEX )
        && analysisType_ != BasePDE::HARMONIC && analysisType_ != BasePDE::INVERSESOURCE && analysisType_ != BasePDE::HARMONIC25D) {
      EXCEPTION( "Can not add integrator '"
                 << biLinContext->GetIntegrator()->GetName()
                 << "' with complex/imaginary entries for a "
                 << "non-harmonic analysis." );
    }

    FEMatrixType mappedFEType = matrixMap_[biLinContext->GetDestMat()];
    FEMatrixType mappedSecFEType = matrixMap_[biLinContext->GetSecDestMat()];

    // Check if integrator can be assembled in this type of simulation
    if( mappedFEType != NOTYPE ) {

      // Store bilinear form
      allBiLinForms_.insert(biLinContext);

      // Note: As the shared_ptr to an Entitylist is not
      // unique within CFS, we have to ensure, that the names of
      // the entity lists rather than the pointers match!

      // Loop over all existing bilinearforms and check,
      // if pair (EntityList1, EntityList2) was already defined
      const std::string ent1Name = biLinContext->GetFirstEntities()->GetName();
      const std::string ent2Name = biLinContext->GetSecondEntities()->GetName();

      std::pair<shared_ptr<EntityList>, shared_ptr<EntityList>> pair;

      bool found = false;
      for (const auto& [key, _] : biLinForms_) {
        const auto& [entity1, entity2] = key;
        if (entity1->GetName() == ent1Name && entity2->GetName() == ent2Name) {
          pair = key;
          found = true;
          break;
        }
      }
      if(!found)
        pair = std::make_pair(biLinContext->GetFirstEntities(), biLinContext->GetSecondEntities());
      biLinForms_[pair].Push_back(biLinContext);
      LOG_DBG(assemble) << "ABF: new map size: " << biLinForms_.size() << " in this entities: " << biLinForms_[pair].GetSize() << ", ptr: " << biLinContext;

      // Pass needed matrix type to algebraic system
      assert(biLinContext->GetFirstFeFunction().lock());
      assert(biLinContext->GetSecondFeFunction().lock());
      FeFctIdType id1 = biLinContext->GetFirstFeFunction().lock()->GetFctId();
      FeFctIdType id2 = biLinContext->GetSecondFeFunction().lock()->GetFctId();
      ReMapFctId(id1);
      ReMapFctId(id2);

      // determine symmetry type and complex status
      bool isSym = IsFEMatSymmetric(biLinContext);
      bool isComplex = IsFEMatComplex(biLinContext);

      algsys_->SetFEMatrixType( mappedFEType, isSym, isComplex, id1, id2 );

      // Check for secondary matrix type
      if( mappedSecFEType != NOTYPE ) {
        algsys_->SetFEMatrixType( mappedSecFEType, isSym, isComplex, id1, id2 );
      }
    }
    else {
      delete biLinContext;
    }
  }

  void Assemble::AddLinearForm( LinearFormContext* linContext ) {

   LOG_DBG(assemble) << "AddLinearForm: " << linContext->ToString() << " on " << linContext->GetEntities()->GetName()
                     << " at=" << BasePDE::analysisType.ToString(analysisType_) << " pde=" << linContext->GetPde()->ToString();

    assert(linContext->GetIntegrator() != NULL);
    assert(linContext->GetPde() != NULL);
    assert(linContext->GetEntities() != NULL);
    assert(linContext->GetPde()->GetAnalysisType() == analysisType_);

    // Store linear form
    allLinForms_.insert(linContext);

    linForms_.Push_back(linContext);

  }

  void Assemble::SetupMatrixGraph(FeFctIdType fctId1, FeFctIdType fctId2 ) {

    StdVector<Integer> eqnVec1, eqnVec2;
    FeFctIdType id1, id2;

    // Perform re-mapping of FctIds
    ReMapFctId( fctId1 );
    ReMapFctId( fctId2 );

    // iterate over all entitylist-pairs and
    BiLinContextListType::iterator listIt = biLinForms_.begin();
    for ( ; listIt != biLinForms_.end(); ++listIt) {
      StdVector<BiLinFormContext*> & forms = listIt->second;
      EntityList& firstEntities = *(listIt->first.first);
      EntityList& secondEntities = *(listIt->first.second);


      EntityIterator it1 = firstEntities.GetIterator();
      EntityIterator it2 = secondEntities.GetIterator();

      // take the maximum of both lists.
      UInt size = std::max( firstEntities.GetSize(),
		                        secondEntities.GetSize() );

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

        if( fctId1 != fctId2 ) {
          setCounterPart = actContext.IsSetCounterPart();
        }

        assert(actContext.GetFirstFeFunction().lock());
        assert(actContext.GetSecondFeFunction().lock());

        // case a)
        if( actContext.GetFirstFeFunction().lock()->GetFctId() == fctId1 &&
            actContext.GetSecondFeFunction().lock()->GetFctId() == fctId2 ) {
          doAssemble = true;
          doTranspose = false;
        } else

          //case b)
          if( actContext.GetFirstFeFunction().lock()->GetFctId() == fctId1 &&
              actContext.GetSecondFeFunction().lock()->GetFctId() == fctId2 &&
              actContext.IsSetCounterPart() ) {
            doAssemble = true;
            doTranspose = true;
          }

        // Check if fctIds of formsContext match
        if( doAssemble ) {

          try {

            NcBiLinFormContext* ncContext = dynamic_cast<NcBiLinFormContext*>(forms[iForm]);

            bool full = false;
            if(ncContext){
              if(ncContext->NeedsFullMatrix())
                full = true;
            }

            if (ncContext && full) {
              // Just get all equations, so we out a dense block in the graph
              ncContext->GetEqns(eqnVec1, eqnVec2, id1, id2);

              // Perform remapping
              ReMapEquations(eqnVec1, id1);
              ReMapEquations(eqnVec2, id2);

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
            } else {
              // Loop over all entities
              it1.Begin();
              it2.Begin();
              //for( ; !(it1.IsEnd() || it2.IsEnd()); it1++, it2++ ) {
              for ( UInt i = 0; i < size; i++ ) {
//                std::cout << "i: " << i << std::endl;

                // Get equation numbers
                actContext.MapEqns( it1, it2, eqnVec1, eqnVec2, id1, id2 );

//                std::cout << "eqnVec1: " << eqnVec1.ToString() << std::endl;
//                std::cout << "eqnVec2: " << eqnVec2.ToString() << std::endl;

                // Perform remapping
                ReMapEquations(eqnVec1, id1);
                ReMapEquations(eqnVec2, id2);

                // Pass entity eqn-connectivity to algebraic system
                if( !doTranspose ) {
//                  std::cout << "No transpose" << std::endl;
                  algsys_-> SetElementPos( id1, eqnVec1,
                      id2, eqnVec2,
                      destMap,
                      setCounterPart );
                } else {
//                  std::cout << "Transpose" << std::endl;
                  algsys_-> SetElementPos( id2, eqnVec2,
                                           id1, eqnVec1,
                                           destMap,
                                           setCounterPart );
                }

                // The size of the entity lists is checked because FeSpaceConst can add single rows/columns.
                if(firstEntities.GetSize() != 1) {
                  it1++;
                }
                if(secondEntities.GetSize() != 1) {
                  it2++;
                }
              } // loop over entities

            }
          } catch (Exception& e) {
            RETHROW_EXCEPTION(e, "Could not setup matrix graph for "
                << "BiLinearForm '"
                << actContext.GetIntegrator()->GetName() << "' on '"
                << actContext.GetFirstEntities()->GetName()<< "'" );
          }

        } // if
      } // loop (integrators)
    } // loop over entity-list pair
  }

  void Assemble::AssembleMatrices(bool isNewtonPart) {

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
      AssembleMatrices_Cond(isNewtonPart);
    } else {
      AssembleMatrices_Std(isNewtonPart);
    }

  }

  void Assemble::AssembleMatrices_Std(bool isNewtonPart) {

    LOG_DBG(assemble) << "AM_Std: AssembleMatrices_Std() enter sequence=" << domain->GetDriver()->GetActSequenceStep();

    timer_->Start();
    matrixTimer_->Start();

    // Reset for matrix update
    matrixUpdated_ = false;

    // Temporary: Check each time for non-linearities
    // On first Assembly, assemble all matrices for each BiLinearForm
    CheckNonLinearities(isFirstTime_);

    // Init all matrices, which have to be reassembled
    // Just to be done, when isNewtonPart is false!
    if ( !isNewtonPart) {
      std::map<FEMatrixType, bool>::iterator it;
      for( it = matReassemble_.begin(); it != matReassemble_.end(); it++ ) {
        if( it->second == true ) {
          LOG_DBG2(assemble) << "AssembleMatrices: init matrix " << it->first;
          algsys_->InitMatrix( matrixMap_[it->first] );
        }
      }
    }

    // iterate over all entitylist-pairs and
    BiLinContextListType::iterator listIt = biLinForms_.begin();
    for ( ; listIt != biLinForms_.end(); ++listIt) {
      StdVector<BiLinFormContext*>& forms = listIt->second;
      EntityList& firstEntities = *(listIt->first.first);
      EntityList& secondEntities = *(listIt->first.second);
      UInt size = std::max(firstEntities.GetSize(), secondEntities.GetSize());

      // Total work: numElement x numForms
      std::stringstream progStream;
      boost::timer::progress_display progress( size*forms.GetSize(), progStream );

      if(printProgressBar_)
        std::cout << "  - Calculating BiLinearForms on '"  << firstEntities.GetName() << " (" << size << " elements)'\n";

      if(!isNewtonPart)
      {
        // First loop over all integrators to check, if any of them
        // gets re-assembled
        // Loop over all bilinearforms
        bool anyReassemble = false;

        for(UInt iForm = 0; iForm < forms.GetSize(); ++iForm)
        {
          BiLinFormContext & actContext = *forms[iForm];

          // get matrix destinations
          FEMatrixType destMat = actContext.GetDestMat();
          FEMatrixType secDestMat = actContext.GetSecDestMat();

          // If assemble was already called and the current destination
          // matrix must not be reassembled -> continue with next iterator
          if(!matReassemble_[destMat] && (secDestMat == NOTYPE || !matReassemble_[secDestMat]))
            continue; // check next form
          anyReassemble = true;
        } // end form loop
        if(!anyReassemble)
          continue; // assemble next bilin form
      } // end !isNewtonPart

#pragma omp parallel num_threads(CFS_NUM_THREADS)
      {
        UInt numT = CFS_NUM_THREADS;
        UInt aThread = GetThreadNum();
        StdVector<BiLinearForm *> biLinForms(forms.GetSize());
        biLinForms.Init(NULL);

        UInt chunksize = std::floor(size/numT);
        UInt start = chunksize * aThread;
        UInt end = (aThread==numT-1)? size : (chunksize * (aThread+1));

        for( UInt iForm = 0; iForm < forms.GetSize(); ++iForm ) {
          //copy bilinear forms
          // TODO: According to fsanatize, this is a huge memory leak...
          biLinForms[iForm] = forms[iForm]->GetIntegrator()->Clone();
        }

//     #pragma omp critical
//         {
//             std::cout << "Thread #" << omp_get_thread_num() << " computing entites from " << start << " to " << end << " for " << end-start << " entities" << std::endl;
//         }


        // Loop over all entities
        EntityIterator it1 = firstEntities.GetIterator();
        EntityIterator it2 = secondEntities.GetIterator();
        //LOG_DBG2(assemble) << "\telems are " << it1.GetElem()->elemNum  << " and " << it2.GetElem()->elemNum;

        it1.Begin();
        it2.Begin();
        //take account for const space
        if(firstEntities.GetSize() != 1 ) {
          it1+=start;
        }
        if( secondEntities.GetSize() != 1 ) {
          it2+=start;
        }

        Matrix<Double> elemMatrix;
        Matrix<Complex> elemMatrixC;
        StdVector<Integer> eqnVec1, eqnVec2;
        FeFctIdType fctId1, fctId2;
        for( UInt i = start;i < end; ++i  ) {

          LOG_DBG2(assemble) << "AM_Std: elems are " << it1.GetIdString() << " and " << it2.GetIdString();

          // Loop over all bilinear forms
          for( UInt iForm = 0; iForm < forms.GetSize(); ++iForm ) {

            BiLinFormContext & actContext = *forms[iForm];
            bool bilinearFormIsNewton =  actContext.IsNewtonBiLinearForm();

            //when we just want to assemble a Newton bilinear form and
            //the bilinear form is not a Newton one, then go to the next!
            if ( isNewtonPart && !bilinearFormIsNewton )
              continue;

            //if we want to assemble all bilinear forms except Newton ones
            //and the bilinear form is a Newton one, then go to next!
            if ( !isNewtonPart && bilinearFormIsNewton )
              continue;

            // get matrix destinations
            FEMatrixType destMat = actContext.GetDestMat();
            FEMatrixType secDestMat = actContext.GetSecDestMat();

            if ( !isNewtonPart) {
              // If assemble was already called and the current destination
              // matrix must not be reassembled -> continue with next iterator
              if(!matReassemble_[destMat] && (secDestMat == NOTYPE || !matReassemble_[secDestMat]))
                continue; // check next form
            }
            // Update flag
            matrixUpdated_ = true;

            BiLinearForm * form =nullptr;
            form = UseOpenMP() ? biLinForms[iForm] : actContext.GetIntegrator();

            LOG_DBG2(assemble) << "AM_Std: bilinform " << form->GetName() << " context=" << actContext.ToString() << " complex=" << form->IsComplex();
            if(!skipElemAssembly_){
              try {

                // make only output if desired
                if( printProgressBar_) {
                  ++progress;
                  std::cout << progStream.str();
                  progStream.str("");
                }

                LOG_DBG3(assemble) << feMatrixType.ToString(destMat) << "\n";
                // Calc element matrix
                if ( form->IsComplex() ){
                  form->CalcElementMatrix( elemMatrixC, it1, it2 );
                  LOG_DBG3(assemble) << "AM_Std: e=" << it1.ToString() << " cplx CEM -> " << elemMatrixC.ToString();
                } else {
                  form->CalcElementMatrix( elemMatrix, it1, it2 );
                  if(it1.IsElemType()) {
                    LOG_DBG3(assemble) << "AM_Std: e=" << it1.GetElem()->elemNum << " reg=" << it1.GetElem()->regionId;
                  }
                  LOG_DBG3(assemble) << "AM_Std: e=" << it1.ToString() << " real CEM -> " << elemMatrix.ToString();
                  if(actContext.IsSetNegate()){
                    assert(!form->IsComplex());
                    elemMatrix*= (-1.0);
                  }
                }

                // info.xml logging in detailed logging case for the first element only
                if(i == 0 && progOpts->DoDetailedInfo())
                {
                  PtrParamNode in = domain->GetInfoRoot()->Get("sequenceStep/PDE")->Get(actContext.GetFirstPde()->GetName())->Get("exemplaryLocalMatrix");

                  // works only for element integrators
                  if(it1.GetType() != EntityList::ELEM_LIST && it1.GetType() != EntityList::SURF_ELEM_LIST && it1.GetType() != EntityList::NC_ELEM_LIST)
                    continue; // no element, no region id

                  // make sure to have only one output for non-static case (e.g. optimization)
                  std::string reg_name = domain->GetGrid()->GetRegion().ToString(it1.GetElem()->regionId);
                  bool found = false;
                  ParamNodeList list = in->GetListByVal("tensor", "form", form->GetName());
                  for(unsigned int li = 0; li < list.GetSize(); li++)
                    if(list[li]->HasByVal("region", reg_name) && list[li]->HasByVal("element", it1.GetElem()->elemNum))
                      found = true;

                  if(!found)
                  {
                    in = in->Get("tensor", ParamNode::APPEND);
                    in->Get("form")->SetValue(form->GetName());
                    in->Get("region")->SetValue(reg_name);
                    in->Get("element")->SetValue(it1.GetElem()->elemNum);
                    // it makes sense to add the analysis id here
                    if(form->IsComplex())
                      in->Get("matrix")->SetValue(elemMatrixC);
                    else
                      in->Get("matrix")->SetValue(elemMatrix);
                  }
                }

                // Map equation numbers
                actContext.MapEqns( it1, it2, eqnVec1, eqnVec2, fctId1, fctId2 );
                // Perform remapping
                ReMapEquations(eqnVec1, fctId1);
                ReMapEquations(eqnVec2, fctId2);

                LOG_DBG3(assemble) << "eqnVec1 [1-based]: " << eqnVec1.ToString() << "\n\n";
                // Pass element matrix to algebraic system (primary matrix)
                if ( form->IsComplex() )
                  InsertMatrix( destMat, actContext, elemMatrixC, eqnVec1, eqnVec2, fctId1, fctId2);
                else
                  InsertMatrix( destMat, actContext, elemMatrix, eqnVec1, eqnVec2, fctId1, fctId2);

                if (secDestMat != NOTYPE ) { // Check for secondary matrix type
                  Double dampFactor = 1.0;
                  // get secondary matrix factor string
                  Double secMatFac = actContext.EvalSecMatFac();

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
                RETHROW_EXCEPTION(e, "Could not calculate element matrix of BiLinearForm '" << form->GetName() << "' on '" << actContext.GetFirstEntities()->GetName()<< "'");
              }
            } //if block to skip assembly
          }// loop over bilinearforms
          // The size of the entity lists is checked because FeSpaceConst can add single rows/columns.
          if( firstEntities.GetSize() != 1 )
            it1++;
          if( secondEntities.GetSize() != 1 )
            it2++;
        } // loop over entities
        for( UInt iForm = 0; iForm < forms.GetSize()&&UseOpenMP(); ++iForm ) {
          //delete copied bilinear forms
          delete biLinForms[iForm];
        }

      }//OMP END
    }// loop over entitylist pairs
    // Change flag
    isFirstTime_ = false;

    // We have assembled all matrices, we will not do so next time
    // except: CheckNonLinearities sets one of these, or Optimization does
    matReassemble_.clear();

    matrixTimer_->Stop();
    timer_->Stop();

       // algsys_->GetMatrix(STIFFNESS)->Export("assemble_stiff.mtx");
       // algsys_->GetMatrix(MASS)->Export("assemble_mass.mtx");
  }

  void Assemble::InitMultHarm() {
    // Reset for matrix update
    matrixUpdated_ = false;

    // Temporary: Check each time for non-linearities
    // On first Assembly, assemble all matrices for each BiLinearForm
    CheckNonLinearities(isFirstTime_);

    // Init all matrices, which have to be reassembled
    // Just to be done, when isNewtonPart is false!
    std::map<FEMatrixType, bool>::iterator it;
    for( it = matReassemble_.begin(); it != matReassemble_.end(); it++ ) {
      if( it->second == true ) {
        LOG_DBG2(assemble) << "InitMultHarm: initialize matrix " << it->first;
        algsys_->InitMatrix( matrixMap_[it->first] );
      }
    }

    timer_->Start();
  }

  void Assemble::PostAssemble(){
    // Change flag
    isFirstTime_ = false;

    // We have assembled all matrices, we will not do so next time
    // except: CheckNonLinearities sets one of these, or Optimization does
    matReassemble_.clear();
  }

  void Assemble::AssembleMatrices_MultHarm(Integer harm, UInt N, UInt M,
             const std::map<RegionIdType, StdVector<NonLinType> >& regionNonLinTypes,
             const StdVector<Double>& multHarmFreqVec) {

    LOG_DBG(assemble) << "AM_Std: AssembleMatricesMultHarm() enter sequence="
        << domain->GetDriver()->GetActSequenceStep() <<" for nu-harmonic " << harm<<std::endl;

    // The checking of nonlinearities and initializing the matrices,
    // which have to be reassembled was already carried out in the
    // StdSolveStep::AssembleMH method

    UInt size = domain->GetDriver()->GetNumFreq();
    // same as ComputeIndex method in GraphManager, here with a lambda function
    auto ComputeIndex = [size](UInt a, UInt b ) { return (size) * a + b;};

    // store the sbm-indices of the blocks, which have to be assembled
    // NOTE: If we are using the optimized version, we only consider odd harmonics,
    // therefore the column isn't anymore iRow+harmonic !!!
    // =========== This is difficult to describe here, see notes (K.Roppert 2018) =====
    StdVector<UInt> sbmInd;
    bool isFullSys = domain->GetDriver()->IsFullSystem();
    for( UInt iRow = 0; iRow < size; ++iRow ) {
      if(harm == 0){
        // sbm-blocks on the main diagonal
        sbmInd.Push_back( ComputeIndex(iRow, iRow) );
      }else{
        Integer col = (isFullSys)? (iRow - harm) : (iRow - harm/2);
        if( col < (Integer)size && col >= 0){
          // change row/columns
          sbmInd.Push_back( ComputeIndex(col, iRow) );
        }
      }
    }

    if (IS_LOG_ENABLED(assemble, dbg3)) {
      std::cout<<"================= sbmInd = "<<sbmInd.ToString()<<std::endl;
    }

    // iterate over all entitylist-pairs and
    BiLinContextListType::iterator listIt = biLinForms_.begin();
    for ( ; listIt != biLinForms_.end(); ++listIt) {
      StdVector<BiLinFormContext*> & forms = listIt->second;
      EntityList& firstEntities = *(listIt->first.first);
      EntityList& secondEntities = *(listIt->first.second);
      StdVector<NonLinType> nonLinTypes(0);
      if( firstEntities.GetType() != 2 && secondEntities.GetType() != 2 ){
        nonLinTypes = regionNonLinTypes.find(firstEntities.GetRegion())->second;
      }
      // If the region is not nonlinear, we do not have to assemble
      // off-diagonal blocks in the multiharmonic system matrix
      if( (nonLinTypes.GetSize() == 0 && harm != 0 )){
        continue;
      }

      UInt size = std::max(firstEntities.GetSize(), secondEntities.GetSize());

#ifdef USE_OPENMP
#pragma omp parallel num_threads(CFS_NUM_THREADS)
    {
      UInt numT = CFS_NUM_THREADS;
      UInt aThread = omp_get_thread_num();
      StdVector<BiLinearForm *> biLinForms(forms.GetSize());

      UInt chunksize = std::floor(size/numT);
      UInt start = chunksize * aThread;
      UInt end = (aThread==numT-1)? size : (chunksize * (aThread+1));

      for( UInt iForm = 0; iForm < forms.GetSize(); ++iForm ) {
        //copy bilinear forms
        biLinForms[iForm] = forms[iForm]->GetIntegrator()->Clone();
      }
#else
      UInt start = 0;
      UInt end = size;
#endif

      // Loop over all entities
      EntityIterator it1 = firstEntities.GetIterator();
      EntityIterator it2 = secondEntities.GetIterator();

      it1.Begin();
      it2.Begin();
      //take account for const space
      if( firstEntities.GetSize() != 1 ) {
        it1+=start;
      }
      if( secondEntities.GetSize() != 1 ) {
        it2+=start;
      }


      Matrix<Double> elemMatrix;
      Matrix<Complex> elemMatrixC;
      StdVector<Integer> eqnVec1, eqnVec2;
      FeFctIdType fctId1, fctId2;
      for( UInt i = start; i < end; ++i  ) {

        LOG_DBG2(assemble) << "\telems are " << it1.GetIdString() << " and " << it2.GetIdString();

        // Loop over all bilinear forms
        for( UInt iForm = 0; iForm < forms.GetSize(); ++iForm ) {


          BiLinFormContext & actContext = *forms[iForm];
          bool bilinearFormIsNewton =  actContext.IsNewtonBiLinearForm();
          if( bilinearFormIsNewton ){
            EXCEPTION("Assemble::AssembleMatrices_MultHarm no newton-bilinear forms"
                      "allowed in multiharmonic analysis");
          }

          // get matrix destinations
          FEMatrixType destMat = actContext.GetDestMat();
          FEMatrixType secDestMat = actContext.GetSecDestMat();


          //================================================================
          // IMPORTANT: In multiharmonic analysis, no off-diagonal
          //            sub mass matrices!
          //================================================================
          if( harm != 0 && (destMat == FEMatrixType::DAMPING || secDestMat == FEMatrixType::DAMPING) ){
            continue;
          }


          // If assemble was already called and the current destination
          // matrix must not be reassembled -> continue with next iterator
          if( matReassemble_[destMat] == false ) {
            if( matReassemble_[secDestMat] != NOTYPE ) {
              if(  matReassemble_[secDestMat] == false ) {
                continue;
              }
            } else  {
              continue;
            }
          }


          // Update flag
          matrixUpdated_ = true;
#ifdef USE_OPENMP
          BiLinearForm * form = biLinForms[iForm];
#else
          BiLinearForm * form = actContext.GetIntegrator();
#endif

          LOG_DBG2(assemble) << "AM_Std: bilinform " << form->GetName() << " context=" << actContext.ToString() << " complex=" << form->IsComplex();

          try {

            // Calc element matrix
            if ( form->IsComplex() ){
              form->CalcElementMatrix( elemMatrixC, it1, it2 );
              LOG_DBG3(assemble) << "AM_Std: cplx CEM -> " << elemMatrixC.ToString();
            } else {
              form->CalcElementMatrix( elemMatrix, it1, it2 );
              LOG_DBG3(assemble) << "AM_Std: real CEM -> " << elemMatrix.ToString();
              if(actContext.IsSetNegate()){
                assert(!form->IsComplex());
                elemMatrix*= (-1.0);
              }
            }


            // info.xml logging in detailed logging case for the first element only
            if(i == 0 && progOpts->DoDetailedInfo())
            {
              PtrParamNode in = domain->GetInfoRoot()->Get("sequenceStep/PDE")->Get(actContext.GetFirstPde()->GetName())->Get("exemplaryLocalMatrix");

              // works only for element integrators
              if(it1.GetType() != EntityList::ELEM_LIST && it1.GetType() != EntityList::SURF_ELEM_LIST && it1.GetType() != EntityList::NC_ELEM_LIST)
                continue; // no element, no region id

              // make sure to have only one output for non-static case (e.g. optimization)
              std::string reg_name = domain->GetGrid()->GetRegion().ToString(it1.GetElem()->regionId);
              bool found = false;
              ParamNodeList list = in->GetListByVal("tensor", "form", form->GetName());
              for(unsigned int li = 0; li < list.GetSize(); li++)
                if(list[li]->HasByVal("region", reg_name) && list[li]->HasByVal("element", it1.GetElem()->elemNum))
                  found = true;

              if(!found)
              {
                in = in->Get("tensor", ParamNode::APPEND);
                in->Get("form")->SetValue(form->GetName());
                in->Get("region")->SetValue(reg_name);
                in->Get("element")->SetValue(it1.GetElem()->elemNum);
                // it makes sense to add the analysis id here
                if(form->IsComplex())
                  in->Get("matrix")->SetValue(elemMatrixC);
                else
                  in->Get("matrix")->SetValue(elemMatrix);
              }
            }

            // Map equation numbers
            actContext.MapEqns( it1, it2, eqnVec1, eqnVec2, fctId1, fctId2 );
            // Perform remapping
            ReMapEquations(eqnVec1, fctId1);
            ReMapEquations(eqnVec2, fctId2);


            // in multiharmonic analysis, the mass matrix must be multiplied by the harmonic number
            // Ideally such an operation would be performed way back in the PDE-class
            // or after the assembling in algsys_->ConstructEffectiveMatrix but the fact that
            // we need to loop over every frequency and multiply the mass matrices corresponding to the
            // frequencies with different values, prevents such a ''clean'' solution
            if( harm == 0 && multHarmFreqVec.GetSize() != 0){
              if( destMat == DAMPING ){
                // Store the sbm-indices of the blocks, which correspond to harmonic 0
                // in a vector with size 1 to pass it to InsertMatrix method.
                // This is kind of a workaround
                StdVector<UInt> diagInd(1);
                for( UInt iRow = 0; iRow < domain->GetDriver()->GetNumFreq(); ++iRow ) {
                  diagInd.Init(0);
                  diagInd[0] =  sbmInd[iRow];
                  Double f = multHarmFreqVec[iRow];

                  // For f = 0, we basically solve the static problem. If we really set the
                  // mass part for harmonic 0 to zero, the solution becomes non-unique.
                  // Therefore we replace here zero with a small relaxation parameter 1e-6
                  if( std::fabs(f) <= 1e-6 ){
                    f = 1.0e-6;
                  }
                  LOG_DBG3(assemble) << "AM_Std: real CEM MASS  in "
                      "SBM block with index "<<diagInd[0]<<" -> "
                      << elemMatrix.ToString();

                  // Pass element matrix to algebraic system (primary matrix)
                  if ( form->IsComplex() ){
                    EXCEPTION("Assembling of complex mass bilinear form in multiharmonic analysis not allowed!");
                  }else{
                    InsertMatrix( destMat, actContext, elemMatrix, eqnVec1, eqnVec2,
                        fctId1, fctId2, false, diagInd, f, true);
                  }
                }
              }else{
                // Pass element matrix to algebraic system (primary matrix)
                if ( form->IsComplex() )
                  InsertMatrix( destMat, actContext, elemMatrixC, eqnVec1, eqnVec2, fctId1, fctId2, false, sbmInd);
                else
                  InsertMatrix( destMat, actContext, elemMatrix, eqnVec1, eqnVec2, fctId1, fctId2, false, sbmInd);
              }

            }else{
              if(actContext.GetDestMat() == DAMPING ){
                EXCEPTION("Assemble::AssembleMatrices_MultHarm: This should not happen");
              }
              // Pass element matrix to algebraic system (primary matrix)
              if ( form->IsComplex() )
                InsertMatrix( destMat, actContext, elemMatrixC, eqnVec1, eqnVec2, fctId1, fctId2, false, sbmInd);
              else
                InsertMatrix( destMat, actContext, elemMatrix, eqnVec1, eqnVec2, fctId1, fctId2, false, sbmInd);
            }

            // increment iterators
          } catch (Exception& e) {
            RETHROW_EXCEPTION(e, "Could not calculate element matrix of "
                << "BiLinearForm '"
                << form->GetName() << "' on '"
                << actContext.GetFirstEntities()->GetName()<< "'" );
          }
        } // loop over bilinearforms

        // The size of the entity lists is checked because FeSpaceConst can add single rows/columns.
        if( firstEntities.GetSize() != 1 ) {
          it1++;
        }
        if( secondEntities.GetSize() != 1 ) {
          it2++;
        }

      } // loop over entities
#ifdef USE_OPENMP
      for( UInt iForm = 0; iForm < forms.GetSize(); ++iForm ) {
        //delete copied bilinear forms
        delete biLinForms[iForm];
      }
    }//OMP END
#endif

    }// loop over entitylist pairs


//algsys_->GetMatrix(SYSTEM)->Export_MultHarm("sysFirstExport", BaseMatrix::MATRIX_MARKET, N, sbmInd);


  }

  void Assemble::TimerStop(){
    if( !timer_->IsRunning() ){
      EXCEPTION("Assemble::TimerStop You want to stop a non-running timer \n"
                "this method should only be called in a multiharmonic analysis");
    }
    timer_->Stop();
  }

  // NOTE: still only for singlePDE problems!
  void Assemble::AssembleMatrices_CondTrans(bool isNewtonPart,UInt currentStage,
                                            std::map<FeFctIdType, std::map<FEMatrixType,Double> > timeStepFactors ){
    // assembling of element matrices for transient simulations
    // NEW: the time integrations factors will be integrated on element level, such that
    //      the resulting matrices can be sort directly into the system matrix
    //      This summation on element level is necessary if static condensation shall be used
    //      as the inversion of the interiour degrees of freedom has to be applied to the
    //      summed up element matrix
    //      To get the correct timestepping factors we need to provide the current stage
    //      and the matMap map (from stdSolveStep)
    //
    // NOTE: still only for singlePDE problems!
    Matrix<Double> elemMatrix;
    Matrix<Complex> elemMatrixC;
    StdVector<Integer> eqnVec1, eqnVec2;
    FeFctIdType fctId1, fctId2;

    timer_->Start();
    matrixTimer_->Start();

    // Reset for matrix update
    matrixUpdated_ = false;

    // Temporary: Check each time for non-linearities
    // On first Assembly, assemble all matrices for each BiLinearForm
    CheckNonLinearities(isFirstTime_);

    // check if any matrix needs to be reassembled
    // if yes we have to reset the whole system matrix as we add the element matrices
    // directly to it
    // if no we can leave this function

    std::map<FEMatrixType, bool>::iterator it;
    for( it = matReassemble_.begin(); it != matReassemble_.end(); it++ ) {
      if( it->second == true ) {
        LOG_DBG2(assemble) << "AssembleMatrices: init matrix " << it->first;
        algsys_->InitMatrix( matrixMap_[it->first] );
        algsys_->InitMatrix( SYSTEM );
        matrixUpdated_ = true;
      }
    }

    if(matrixUpdated_ == false){
      return;
    }

    // create temporary matrices for each pair of function ids (mech-mech,mech-acou,acou-mech,...)
    // and for each of these pairs a matrix of each type (Stiffness, mass, ...)
    std::map< std::pair<FeFctIdType,FeFctIdType>, std::map<FEMatrixType, Matrix<Double> > > rElemMats;
    //std::map< std::pair<FeFctIdType,FeFctIdType>, std::map<FEMatrixType, Matrix<Complex> > > cElemMats;

    // and we also have to store for each fctIdPair the equation vectors
    std::map< std::pair<FeFctIdType,FeFctIdType>, std::pair<StdVector<Integer>,StdVector<Integer> > > eqnVecs;
    // and for the actual context
    std::map< std::pair<FeFctIdType,FeFctIdType>, bool > counterParts;

    // iterate over all entitylist-pairs
    BiLinContextListType::iterator listIt = biLinForms_.begin();
    for ( ; listIt != biLinForms_.end(); ++listIt) {

      // get vector of all formcontexts
      StdVector<BiLinFormContext*> & forms = listIt->second;

      // get lists of elements/entities
      EntityList& firstEntities = *(listIt->first.first);
      EntityList& secondEntities = *(listIt->first.second);

      // total number of elements/entities in the first entitylist
      UInt size = firstEntities.GetSize();

      if( printProgressBar_) {
        std::cout << "  - Calculating BiLinearForms on '"
            << firstEntities.GetName()
            << " (" << size << " elements)'\n";
      }
      // Total work: numElement x numForms
      std::stringstream progStream;
      boost::timer::progress_display progress(size*forms.GetSize(), progStream);

      // Loop over all elements/entities
      EntityIterator it1 = firstEntities.GetIterator();
      EntityIterator it2 = secondEntities.GetIterator();
      for( it1.Begin(); !it1.IsEnd(); it1++, it2++ ) {

        // clear all maps (assume that each element is assigned only to one entitylist!)
        eqnVecs.clear();
        rElemMats.clear();
        //cElemMats.clear();

        LOG_DBG2(assemble) << "\telems are " << it1.GetElem()->elemNum
            << " and " << it2.GetElem()->elemNum;
        try {
          // Loop over all bilinearforms
          // for each bilinearform sort element matrix not directly into the system
          // but sort it into the rElemMats resp. cElemMats
          for( UInt iForm = 0; iForm < forms.GetSize(); ++iForm ) {

            // get current context
            BiLinFormContext & actContext = *forms[iForm];

            // get matrix destinations
            FEMatrixType destMat = actContext.GetDestMat();
            FEMatrixType secDestMat = actContext.GetSecDestMat();

            // get function ids
            fctId1 = actContext.GetFirstFeFunction().lock()->GetFctId();
            fctId2 = actContext.GetSecondFeFunction().lock()->GetFctId();

            // get eqnVecs
            actContext.MapEqns( it1, it2, eqnVec1, eqnVec2, fctId1, fctId2 );

            // perform remapping
            ReMapEquations(eqnVec1, fctId1);
            ReMapEquations(eqnVec2, fctId2);

            // build pair for map
            std::pair<FeFctIdType,FeFctIdType> fctIdPair(fctId1,fctId2);

            // sort equations into the map
            std::pair<StdVector<Integer>,StdVector<Integer> > eqnPair (eqnVec1,eqnVec2);
            eqnVecs[fctIdPair] = eqnPair;

            // sort context to map
            counterParts[fctIdPair] = actContext.IsSetCounterPart();


            //          FEMatrixType secDestMat = actContext.GetSecDestMat();

            // get secondary matrix factor string
            //          Double secMatFac = actContext.EvalSecMatFac();

            // If assemble was already called and the current destination
            // matrix must not be reassembled -> continue with next iterator
            //          if( matReassemble_[destMat] == false ) {
            //            continue;
            //          }
            //
            //          // Update flag
            //          matrixUpdated_ = true;

            BiLinearForm * form = actContext.GetIntegrator();

            // make only output if desired
            if( printProgressBar_) {
              ++progress;
              std::cout << progStream.str();
              progStream.str("");
            }

            // Calc element matrix
            if ( form->IsComplex() ){
              EXCEPTION("Complex form not implemented yet");
              //form->CalcElementMatrix( elemMatrixC, it1, it2 );
              //cElemMats[fctIdPair][destMat] += elemMatrixC;
            } else {
              form->CalcElementMatrix( elemMatrix, it1, it2 );
              if(actContext.IsSetNegate()== true){
                assert(!form->IsComplex());
                elemMatrix*= (-1.0);
              }

              // insert element matrix to dest (needed later to construct rhs)
              InsertMatrix( destMat, actContext, elemMatrix,
                            eqnVec1, eqnVec2,
                            fctId1, fctId2, true);


              if(rElemMats[fctIdPair].count(destMat) == 0){
                rElemMats[fctIdPair][destMat] = elemMatrix;
              } else {
                rElemMats[fctIdPair][destMat] += elemMatrix;
              }
            }

            // calc element matrix for second destination (for damping!)
            if (secDestMat != NOTYPE ) { // Check for secondary matrix type
              Double dampFactor = 1.0;
              // get secondary matrix factor string
              Double secMatFac = actContext.EvalSecMatFac();

              // damping with "exotic" complex material
              if ( form->IsComplex() ) {
                EXCEPTION("Complex form not implemented yet");
                // complex damping
                //elemMatrixC *= secMatFac * dampFactor;

                // Pass secondary matrix part to algebraic system
                //InsertMatrix(secDestMat, actContext, elemMatrixC, eqnVec1, eqnVec2, fctId1, fctId2, true);
              }
              // "standard" Rayleigh damping. Includes the standard SIMP optimization!
              else {
                // Rayleigh damping
                elemMatrix *= secMatFac * dampFactor;
                LOG_DBG3(assemble) << "AM: e1=" << it1.GetElem()->elemNum << " Rayleigh damping form=" << form->GetName() << " sMF=" << secMatFac << " df=" <<  dampFactor;
                // Pass secondary matrix part to algebraic system
                InsertMatrix(secDestMat, actContext, elemMatrix, eqnVec1, eqnVec2, fctId1, fctId2, true);

                if(rElemMats[fctIdPair].count(secDestMat) == 0){
                  rElemMats[fctIdPair][secDestMat] = elemMatrix;
                } else {
                  rElemMats[fctIdPair][secDestMat] += elemMatrix;
                }
              }
            } // handle secDestMat != NOTYPE

          } // loop over bilinearforms
        } catch (Exception& e) {
          RETHROW_EXCEPTION(e, "Could not calculate element matrix of "
                            << "BiLinearForm"  );
        }


        // up till now we have gathered for the current element all element-matrix-parts
        // sorted by function IDs
        // in the next step we loop over all function IDs and combine the element-matrix-parts
        // with the corresponding time-stepping factors to the complete element-matrix (w.r.t. a fctId-pair)
        // These element matrices can be passed to the algsys which can apply static condensation to it
        // NOTE: static condensation only for single fields so far! for coupled fields the element matrices for
        //       each fctId-pair have to be collected in a super-element matrix to which static condensation can
        //       then be applied

        std::map< std::pair<FeFctIdType,FeFctIdType>, bool >::iterator mapIt;
        // loop over all fucntionIds pairs
        for(mapIt = counterParts.begin(); mapIt != counterParts.end(); ++mapIt){

          // get key pair
          std::pair<FeFctIdType,FeFctIdType> keyPair = mapIt->first;

          // combine timestepping factors with matrices
          Matrix<Double> rElemMatSummed;
          //Matrix<Complex> cElemMatSummed;
          UInt count = 0;

          std::map<FEMatrixType,Double>::iterator timeIt;
          for ( timeIt = timeStepFactors[keyPair.first].begin(); timeIt != timeStepFactors[keyPair.first].end(); timeIt++ ) {

            if ( rElemMats[keyPair].count(timeIt->first) == 1 && (*timeIt).second != 0.0 ) {

              if(count == 0){
                rElemMatSummed = rElemMats[keyPair][(*timeIt).first] * timeIt->second;
              }else{
                rElemMatSummed += rElemMats[keyPair][(*timeIt).first] * timeIt->second;
              }
              count++;
            }
            //              if ( cElemMats[keyPair][(*it).first] != NULL  && (*it).second != 0.0 ) {
            //                if(count == 0){
            //                  cElemMatSummed = it->second*cElemMats[keyPair][(*it).first];
            //                }else{
            //                  cElemMatSummed += it->second*cElemMats[keyPair][(*it).first];
            //                }
            //              }
          }

          // now we can sort in the summed up matrices
          algsys_->SetElementMatrix( SYSTEM, rElemMatSummed,
                                     keyPair.first, eqnVecs[keyPair].first,
                                     keyPair.second, eqnVecs[keyPair].second,
                                     mapIt->second,false,false);

        }

      } // loop over entities
    }// loop over entitylist pairs
    // Change flag
    isFirstTime_ = false;

    // We have assembled all matrices, we will not do so next time
    // except: CheckNonLinearities sets one of these, or Optimization does
    matReassemble_.clear();

    timer_->Stop();
    matrixTimer_->Stop();

  }

    void Assemble::AssembleMatrices_Cond(bool isNewtonPart) {

    Matrix<Double> elemMatrix;
    Matrix<Complex> elemMatrixC;
    StdVector<Integer> eqnVec1, eqnVec2;
    FeFctIdType fctId1, fctId2;

    timer_->Start();
    matrixTimer_->Start();

    // Reset for matrix update
    matrixUpdated_ = false;

    // Temporary: Check each time for non-linearities
    // On first Assembly, assemble all matrices for each BiLinearForm
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

      if( printProgressBar_) {
        std::cout << "  - Calculating BiLinearForms on '"
            << firstEntities.GetName()
            << " (" << size << " elements)'\n";
      }
      // Total work: numElement x numForms
      std::stringstream progStream;
      boost::timer::progress_display progress( size*forms.GetSize(), progStream );

      // Loop over all entities
      EntityIterator it1 = firstEntities.GetIterator();
      EntityIterator it2 = secondEntities.GetIterator();
      for( it1.Begin(); !it1.IsEnd(); it1++, it2++ ) {
        LOG_DBG2(assemble) << "\telems are " << it1.GetIdString()
                           << " and " << it2.GetIdString();
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

          BiLinearForm * form = actContext.GetIntegrator();

          // make only output if desired
          if( printProgressBar_) {
            ++progress;
            std::cout << progStream.str();
            progStream.str("");
          }

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
              if(iForm == 0)
                rElemMat = elemMatrix;
              else
                rElemMat += elemMatrix;
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
            // Perform remapping
            ReMapEquations(eqnVec1, fctId1);
            ReMapEquations(eqnVec2, fctId2);

        } // loop over bilinearforms    // increment iterators
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
    matrixTimer_->Stop();
  }


  void Assemble::AssembleLinRHS()
  {
    AssembleRHSLinForms(false);
  }

  void Assemble::AssembleNonLinRHS() {
    AssembleRHSLinForms(true);
  }

  void Assemble::AssembleRHSLinForms(bool nonLin ) {

    python->CallHook(PythonKernel::ASSEMBLE_RHS);

    StdVector<Integer> eqnVec;
    FeFctIdType fctId;
    StdVector<LinearFormContext*>::iterator formsIt;

    timer_->Start();
    rhsTimer_->Start();

    // iterate over all descriptors
    for(formsIt = linForms_.Begin(); formsIt != linForms_.End(); formsIt++)
    {
      // get integrator
      LinearFormContext& actContext = **formsIt;

      // Check, if lin/non-lin type of Context matches parameter nonLin
      // For multiharmonic analysis, we don't reassemble the RHS
      if( (actContext.IsNonLin() != nonLin) && !algsys_->IsMultHarm() )
        continue; //TODO: uncomment this

      LinearForm* form = actContext.GetIntegrator();

      UInt h = form->GetHarm();

      try
      {
        // get entity iterator
        EntityIterator  entIt = actContext.GetEntities()->GetIterator();
        UInt size = actContext.GetEntities()->GetSize();

        if(printProgressBar_)
          std::cout << "  - Calculating '" << form->GetName() << "' on '" << actContext.GetEntities()->GetName() << " (" << size << " elements)'\n";

        LOG_DBG(assemble) << "ARLF: form=" << form->GetName() << " on " << actContext.GetEntities()->GetName() << " (" << size << ") at=" << BasePDE::analysisType.ToString(analysisType_);
        LOG_DBG(assemble) << "ARLF: ac-pde=" << actContext.GetPde()->ToString();
        assert(analysisType_ == actContext.GetPde()->GetAnalysisType());

        std::stringstream progStream;
        boost::timer::progress_display progress( size, progStream );

        if ( analysisType_ == BasePDE::HARMONIC || analysisType_ == BasePDE::MULTIHARMONIC || analysisType_ == BasePDE::INVERSESOURCE || analysisType_ == BasePDE::HARMONIC25D) {

          Vector<Complex> elemVec;
          for ( entIt.Begin(); !entIt.IsEnd(); entIt++ ) {
            // make only output if desired
            if( printProgressBar_) {
              ++progress;
              std::cout << progStream.str();
              progStream.str("");
            }

            // Calculate complex valued element vector
            form->CalcElemVector( elemVec, entIt );

            // Map equation numbers
            actContext.MapEqns( entIt, eqnVec, fctId );
            // Perform remapping
            ReMapEquations(eqnVec, fctId);

            assert(!elemVec.ContainsNaN() && !elemVec.ContainsInf());
            algsys_-> SetElementRHS(elemVec, fctId, eqnVec, h);
            LOG_DBG3(assemble) << "ARLF: fctId=" << fctId << " elemVec=" << elemVec.ToString() << " eqnVec=" << eqnVec.ToString();
          }
        } else {
          // That should be STATIC, TRANSIENT, EIGENFREQUENCY, EIGENVALUE or BUCKLING
          Vector<Double> elemVec;
          // iterate over all entities
          for ( entIt.Begin(); !entIt.IsEnd(); entIt++ ) {
            // make only output if desired
            if( printProgressBar_) {
              ++progress;
              std::cout << progStream.str();
              progStream.str("");
            }

            LOG_DBG3(assemble) << "ARLF: ent=" << entIt.GetPos() << "/" << entIt.GetSize() << " el=" << entIt.ToString() << " fctId=" << fctId;
            LOG_DBG3(assemble) << "ARLF: linform: " << form->ToString();

            // Calculate real valued element vector
            // check if only the real part of a complex value shall be considered
            if( form->IsExtractReal() ){
            	Vector<Complex> tmp;
            	form->CalcElemVector(tmp, entIt);
            	elemVec = tmp.GetPart(Global::REAL);
            }else{
              form->CalcElemVector(elemVec, entIt);
            }
            
            // Map equation numbers. eqnVec can be empty if nodes are not in system (e.g. not simulated region)
            actContext.MapEqns(entIt, eqnVec, fctId);
            LOG_DBG3(assemble) << "ARLF: fctId=" << fctId << " map eqnVec=" << eqnVec.GetSize() << " -> " << eqnVec.ToString();

            // Perform remapping
            ReMapEquations(eqnVec, fctId);
            LOG_DBG3(assemble) << "ARLF: fctId=" << fctId << " remap eqnVec=" << eqnVec.ToString();

            // elemVec can be independent on equations, e.g. for const or expression form
            assert(!elemVec.ContainsNaN() && !elemVec.ContainsInf());

            // Pass element vector to algebraic system only when we have equations
            if(!eqnVec.IsEmpty())
              algsys_->SetElementRHS(elemVec, fctId, eqnVec, h);
            LOG_DBG3(assemble) << "ARLF: fctId=" << fctId << " elemVec=" << elemVec.ToString() << " eqnVec=" << eqnVec.ToString();
          }
        }
      } catch (Exception& e) {
        RETHROW_EXCEPTION(e, "Could not calculate RHS vector for LinearForm '"
                          << form->GetName() << "' on '"
                          << actContext.GetEntities()->GetName() << "'" );
      }
    }
    rhsTimer_->Stop();
    timer_->Stop();
  }

  void Assemble::ToInfo(PtrParamNode in)
  {
    PtrParamNode list = in->Get("matrixBiLinearForms");

    // iterate over all descriptors
    std::set<BiLinFormContext*>::iterator it;
    for ( it = allBiLinForms_.begin(); it != allBiLinForms_.end(); it++ )
    {
      // get integrator
      BiLinFormContext& context = **it;
      BiLinearForm*     form    = context.GetIntegrator();

      PtrParamNode inf = list->Get("bilinearForm", ParamNode::APPEND);
      // integrator name
      inf->Get("integrator")->SetValue(form->GetName());
      inf->Get("type")->SetValue(form->type.ToString(form->GetType()));
      inf->Get("complex")->SetValue(form->IsComplex());

      BaseBDBInt* bdb = dynamic_cast<BaseBDBInt*>(form);
      if(bdb != NULL && bdb->GetBOp() != NULL && bdb->GetBOp()->GetName() != "")
        inf->Get("BOp")->SetValue(bdb->GetBOp()->GetName());
      // coef function for value
      if(bdb != NULL) {
        // bdb->GetCoef()->GetName() is the class Name like CoefFunctionConst
        // the ToString() is either empty, the description or something useful
        string val =bdb->GetCoef()->ToString();
        LOG_DBG2(assemble) << "TI: coef=" << bdb->GetCoef()->GetName() << " val=" << val;
        if(val == bdb->GetCoef()->GetName() || val == "" || (val.size() > 60 && !progOpts->DoDetailedInfo()))
          inf->Get("coef")->SetValue(bdb->GetCoef()->GetName()); // attribute to linearForm
        else {
          inf->Get("coef/type")->SetValue(bdb->GetCoef()->GetName()); // attribute to linearForm
          inf->Get("coef/value")->SetValue(val);
        }
      }

      // region name of entity list
      std::string regionName;
      if( context.GetFirstEntities()->GetType() == EntityList::ELEM_LIST )
      {
        shared_ptr<ElemList> list = dynamic_pointer_cast<ElemList,EntityList>(context.GetFirstEntities());
        regionName = list->GetName();
      }
      else if ( context.GetFirstEntities()->GetType() == EntityList::SURF_ELEM_LIST )
      {
        shared_ptr<SurfElemList> list = dynamic_pointer_cast<SurfElemList,EntityList>(context.GetFirstEntities());
        regionName = list->GetName();
      }
      inf->Get("region")->SetValue(regionName);

      // add information about row / column coordinate
      PtrParamNode row = inf->Get("row", ParamNode::APPEND);
      PtrParamNode col = inf->Get("column", ParamNode::APPEND);

     // associated PDEs
      row->Get("pde")->SetValue(context.GetFirstPde()->GetName());
      col->Get("pde")->SetValue(context.GetSecondPde()->GetName());

      // associated FeFunctions
      assert(context.GetFirstFeFunction().lock());
      assert(context.GetSecondFeFunction().lock());
      row->Get("functionId")->SetValue(context.GetFirstFeFunction().lock()->GetFctId());
      col->Get("functionId")->SetValue(context.GetSecondFeFunction().lock()->GetFctId());

      // associated result types
      std::string tmp;
      if (context.GetFirstResultInfo()!= NULL) {
        tmp = SolutionTypeEnum.ToString(context.GetFirstResultInfo()->resultType);
        row->Get("result")->SetValue(tmp);
        tmp = SolutionTypeEnum.ToString(context.GetSecondResultInfo()->resultType);
        col->Get("result")->SetValue(tmp);
      }


      // matrix destination
      PtrParamNode dest = inf->Get("destination", ParamNode::APPEND);

      // original destination matrix
      Enum2String(context.GetDestMat(), tmp );
      dest->Get("feMatrix")->SetValue(tmp);

      // mapped destination matrix
      Enum2String(matrixMap_[context.GetDestMat()], tmp );
      dest->Get("feMatrixMapped")->SetValue(tmp);

      // secondary destination matrix and factor
      Enum2String(context.GetSecDestMat(), tmp );
      dest->Get("feSecondMatrix")->SetValue(tmp);
      dest->Get("feSecondMatrixFac")->SetValue(context.GetSecMatFac());

      // additional attributes
      PtrParamNode attr = inf->Get("attributes", ParamNode::APPEND);

      // entry Type (real / imag)
      tmp = ComplexPartEnum.ToString(context.GetEntryType());
      attr->Get("entryType")->SetValue( tmp );

      // flag setcounterpart
      tmp = context.IsSetCounterPart() ? "yes" : "no";
      attr->Get("counterPart")->SetValue( tmp );

      // issymmetric
      tmp = context.GetIntegrator()->IsSymmetric() ? "yes" : "no";
      attr->Get("symmetric")->SetValue( tmp );

      // isSolDependent
      tmp = context.GetIntegrator()->IsSolDependent() ? "yes" : "no";
      attr->Get("solutionDependent")->SetValue( tmp );

      // updated geometry
      tmp = context.GetIntegrator()->IsCoordUpdate() ? "yes" : "no";
      attr->Get("updatedGeo")->SetValue( tmp );

    }

    list = in->Get("rhsLinearForms");

    // iterate over all descriptors
    StdVector<LinearFormContext*>::iterator linIt;
    for (linIt = linForms_.Begin(); linIt != linForms_.End(); linIt++)
    {
      PtrParamNode form = list->Get("linearForm", ParamNode::APPEND);

      // get integrator
      LinearFormContext& context = **linIt;

      form->Get("integrator")->SetValue(context.GetIntegrator()->GetName());

      shared_ptr<EntityList> el = context.GetEntities();
      form->Get("entities/name")->SetValue(el->GetName());
      form->Get("entities/size")->SetValue(el->GetSize());
      form->Get("entities/region")->SetValue(el->GetRegion() != NO_REGION_ID ? domain->GetGrid()->GetRegion().ToString(el->GetRegion()) : "no_region");

      SingleEntryInt* sei = dynamic_cast<SingleEntryInt*>(context.GetIntegrator());
      if(sei) {
        string val = sei->GetCoef()->ToString();
        LOG_DBG2(assemble) << "TI: coef=" << sei->GetCoef()->GetName() << " val=" << val;
        if(val == sei->GetCoef()->GetName() || val == "" || (val.size() > 70 && !progOpts->DoDetailedInfo()))
          form->Get("coef")->SetValue(sei->GetCoef()->GetName()); // attribute to linearForm
        else {
          form->Get("coef/type")->SetValue(sei->GetCoef()->GetName()); // attribute to linearForm
          form->Get("coef/value")->SetValue(val);
        }
      }
      // region name of entity list
      std::string regionName;
      if( context.GetEntities()->GetType() == EntityList::ELEM_LIST ) {
        shared_ptr<ElemList> list = dynamic_pointer_cast<ElemList,EntityList>(context.GetEntities());
        regionName = list->GetName();
      }
      else if ( context.GetEntities()->GetType() == EntityList::SURF_ELEM_LIST ) {
        shared_ptr<SurfElemList> list = dynamic_pointer_cast<SurfElemList,EntityList>(context.GetEntities());
        regionName = list->GetName();
      }

      // add information about row / column coordinate
      PtrParamNode row = form->Get("row", ParamNode::APPEND);
      // associated PDEs
      row->Get("pde")->SetValue(context.GetPde()->GetName());

      // associated result types
      std::string tmp;
      tmp = SolutionTypeEnum.ToString(context.GetResultInfo()->resultType);
      row->Get("result")->SetValue(tmp);

      // isSolDependent
      tmp = context.GetIntegrator()->IsSolDependent() ? "yes" : "no";
      row->Get("solutionDependent")->SetValue( tmp );
    }
  }


  void Assemble::CheckNonLinearities(bool setall) {
    // iterate over all bilinearforms
    std::set<BiLinFormContext*>::iterator it;

    for(BiLinFormContext* actContext : allBiLinForms_) {
      // we set multiple times in eigenfrequency for bloch and there we need to reassemble
      if(actContext->IsNonLin() || analysisType_ == BasePDE::HARMONIC || analysisType_ == BasePDE::MULTIHARMONIC
		     || analysisType_ ==BasePDE::INVERSESOURCE || analysisType_ == BasePDE::EIGENFREQUENCY || setall)
      {
        matReassemble_[actContext->GetDestMat()] = true;
        if ( actContext->GetSecDestMat() != NOTYPE )
          matReassemble_[actContext->GetSecDestMat()] = true;
      }
    }
    // Now we know which matrices are nonlinear (e.g. due to nonlinear stiffnes integrator)
    // However, due to the secondaryMatrix-mechanism it could happen, that initially only
    // the STIFFNESS matrix is set to reassemble. Due to the secondary matrix factor of
    // the linear stiffness integrator, also the DAMPING matrix has to be re-assembled
    // (first additional loop). In a next loop, we determine, that also the MASS integrator
    // has to be re-assembled, as his secondary-matrix is the DAMPING one, which also
    // has to be re-assembled (second loop). So to be on the save side and resolve
    // all dependencies (i.e. all matrices have to be re-assembled), we perform
    // the check three times.

    for( UInt i = 0; i < 3; ++i ) {
      for(BiLinFormContext* actContext : allBiLinForms_) {
        bool oneIsNonLin = false;

        // check primary or secondary matrix is nonlinear
        if( matReassemble_[actContext->GetDestMat()] == true ||
            ( actContext->GetSecDestMat() != NOTYPE &&
                matReassemble_[actContext->GetSecDestMat()] == true) ) {
          oneIsNonLin = true;
        }
        if( oneIsNonLin ) {
          matReassemble_[actContext->GetDestMat()] = true;
          if( actContext->GetSecDestMat() != NOTYPE )
            matReassemble_[actContext->GetSecDestMat()] = true;
        }
      } // loops over integrators
    } // 3 loops
    if( IS_LOG_ENABLED( assemble, dbg) ) {
      // Finally print status of matrices
      std::map<FEMatrixType, bool>::const_iterator it2 =  matReassemble_.begin();
      LOG_DBG(assemble) <<  "Status for matrix re-assembly: ";
      for( ; it2 != matReassemble_.end(); ++it2) {
        LOG_DBG(assemble) <<  "\t" << feMatrixType.ToString(it2->first)
                    << ": " << (it2->second ? "true" : "false");
      }
    }
  }

  void Assemble::Matrix2Complex(Matrix<Complex>& complexMat,Matrix<Double>& origMat){
    complexMat.Resize( origMat.GetNumRows(), origMat.GetNumCols() );
    complexMat.SetPart( Global::REAL, origMat );
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
      case STIFFNESS_UPDATE:
        derivOrder = 0;
        factor = 1;
        break;
      case DAMPING:
        derivOrder = 1;
        factor = omega;
        break;
      case DAMPING_AUX:
        derivOrder = 1;
        factor = 1.0 / omega;
        break;
      case DAMPING_UPDATE:
        derivOrder = 1;
        factor = omega;
        break;
      case MASS:
        derivOrder = 2;
        factor = -omega*omega;
        break;
      case MASS_UPDATE:
        derivOrder = 2;
        factor = -omega*omega;
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
    case STIFFNESS_UPDATE:
      factor = Complex(1.0, 0.0);
      break;
    case DAMPING:
      factor = Complex(0.0, omega);
      break;
    case DAMPING_AUX:
          factor = Complex(0.0, 1.0 / omega);
          break;
    case DAMPING_UPDATE:
      factor = Complex(0.0, omega);
      break;
    case MASS:
      // BLOCH CHECK for 1st time derivative order!
      factor = Complex(-omega*omega, 0.0);
      break;
    case MASS_UPDATE:
      // BLOCH CHECK for 1st time derivative order!
      factor = Complex(-omega*omega, 0.0);
      break;
    default:
      EXCEPTION("No default conversion from double entries to matrix type" << matrixType << "known");
      break;
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
      matrixMap_[STIFFNESS_UPDATE] = SYSTEM;
      matrixMap_[DAMPING_UPDATE]   = NOTYPE;
      matrixMap_[MASS_UPDATE]      = NOTYPE;
      break;

    case BasePDE::TRANSIENT:
      matrixMap_[SYSTEM]    = SYSTEM;
      matrixMap_[STIFFNESS] = STIFFNESS;
      matrixMap_[DAMPING]   = DAMPING;
      matrixMap_[MASS]      = MASS;
      matrixMap_[AUXILIARY] = AUXILIARY;
      matrixMap_[STIFFNESS_UPDATE] = STIFFNESS_UPDATE;
      matrixMap_[DAMPING_UPDATE]   = DAMPING_UPDATE;
      matrixMap_[MASS_UPDATE]      = MASS_UPDATE;
      break;

    case BasePDE::HARMONIC:
      matrixMap_[SYSTEM]    = SYSTEM;
      matrixMap_[STIFFNESS] = SYSTEM;
      matrixMap_[DAMPING]   = SYSTEM;
      matrixMap_[DAMPING_AUX]   = SYSTEM;
      matrixMap_[MASS]      = SYSTEM;
      matrixMap_[AUXILIARY] = AUXILIARY; // optimization for radiation needs this
      matrixMap_[STIFFNESS_UPDATE] = SYSTEM;
      matrixMap_[DAMPING_UPDATE]   = SYSTEM;
      matrixMap_[MASS_UPDATE]      = SYSTEM;
      break;

    case BasePDE::MULTIHARMONIC:
      matrixMap_[SYSTEM]    = SYSTEM;
      matrixMap_[STIFFNESS] = SYSTEM;
      matrixMap_[DAMPING]   = SYSTEM;
      matrixMap_[MASS]      = SYSTEM;
      matrixMap_[STIFFNESS_UPDATE] = SYSTEM;
      matrixMap_[DAMPING_UPDATE]   = SYSTEM;
      matrixMap_[MASS_UPDATE]      = SYSTEM;
      break;


    case BasePDE::INVERSESOURCE:
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
      matrixMap_[STIFFNESS_UPDATE] = STIFFNESS;
      matrixMap_[DAMPING_UPDATE]   = DAMPING;
      matrixMap_[MASS_UPDATE]      = MASS;
      break;

    case BasePDE::BUCKLING:
      matrixMap_[SYSTEM]    = NOTYPE;
      matrixMap_[STIFFNESS] = STIFFNESS;
      matrixMap_[GEOMETRIC_STIFFNESS] = GEOMETRIC_STIFFNESS;
      break;

    case BasePDE::EIGENVALUE:
        matrixMap_[SYSTEM]    = NOTYPE;
        matrixMap_[STIFFNESS] = STIFFNESS;
        matrixMap_[DAMPING]   = DAMPING;
        matrixMap_[MASS]      = MASS;
        matrixMap_[AUXILIARY] = NOTYPE;
        matrixMap_[STIFFNESS_UPDATE] = STIFFNESS;
        matrixMap_[DAMPING_UPDATE]   = DAMPING;
        matrixMap_[MASS_UPDATE]      = MASS;
        break;
    
    case BasePDE::HARMONIC25D:
        matrixMap_[SYSTEM]    = SYSTEM;
        matrixMap_[STIFFNESS] = STIFFNESS;
        matrixMap_[DAMPING]   = DAMPING;
        matrixMap_[DAMPING_AUX] = DAMPING_AUX;
        matrixMap_[MASS]      = MASS;
        break;

    default:
      EXCEPTION("Analysistype '" << BasePDE::analysisType.ToString(analysisType_)
                << "' not known!");
    }
  }

  void Assemble::ReMapEquations( StdVector<Integer>&  eqns,
                                 FeFctIdType& fctId) {
    if( customEqnMap_.size() == 0 )
      return;

    UInt numEqns = eqns.GetSize();
    StdVector<Integer> tmp(numEqns);
    for( UInt i = 0; i < numEqns; ++i )  {
      tmp[i] = customEqnMap_[eqns[i]];
    }
    eqns = tmp;
    fctId = customFctIdMap_[fctId];
  }

  void Assemble::ReMapFctId(  FeFctIdType& fctId ) {
    if( customFctIdMap_.size() != 0 ) {
      std::map<FeFctIdType, FeFctIdType>::const_iterator it;
      it = customFctIdMap_.find(fctId);
      if( it != customFctIdMap_.end() ) {
        fctId = it->second;
      } else {
        EXCEPTION("Can not re-map FeFct with id '" << fctId << "'");
      }
    }
  }

  bool Assemble::IsFEMatSymmetric( BiLinFormContext* actCt  ) {

//    if(analysisType_ == BasePDE::MULTIHARMONIC){
//      WARN("===================================================================== \n"
//           "         Developer warning: Check the symmetry of the matrix \n"
//           "         in Assemble::IsFEMatSymmetric \n"
//           "===================================================================== \n");
//    }

    bool isSymmetric = true;
    // set collecting all diagonal-positions
    std::set<shared_ptr<ResultInfo> > allResults, diagResults;

    // Check, where bilinearform gets assembled to diagonal block
    if( (actCt->GetFirstFeFunction().lock() == actCt->GetSecondFeFunction().lock())
        && (actCt->GetFirstEntities() == actCt->GetSecondEntities() ) ) {

      // Bilinearform gets assembled to main diagonal.
      // If bilinearform is non-symmetric, so is the related FE-matrix
      if( !actCt->GetIntegrator()->IsSymmetric()  ) {
        isSymmetric = false;
      }
    } else {

      // BiLinearform gets assembled to off-diagonal block.

      // If the bilinearorm is also assembled to the transposed block
      // we assume that the matrix still remains symmetric.
      // Otherwise we assume, that we need a non-symmetric matrix.
      if( !actCt->IsSetCounterPart() ) {
        isSymmetric= false;
      }
    }
    // return flag for matrix of interest
    return isSymmetric;
  }

  bool Assemble::IsFEMatComplex( BiLinFormContext* actCt  ) {

    bool isComplex = false;
    if (actCt->GetIntegrator()->IsComplex() ||
        analysisType_ == BasePDE::HARMONIC || analysisType_ == BasePDE::MULTIHARMONIC || analysisType_ == BasePDE::INVERSESOURCE || analysisType_ == BasePDE::EIGENVALUE || analysisType_ == BasePDE::HARMONIC25D) {
      isComplex = true;
    }
    return isComplex;
  }

  bool Assemble::IsRhsSolDependent() {
    bool ret = false;
    UInt size = linForms_.GetSize();
    for( UInt i = 0; i < size; ++i ) {
      ret |= linForms_[i]->IsNonLin();
    }
    return ret;
  }

  void Assemble::InsertMatrix( FEMatrixType dest, BiLinFormContext& context,
                               Matrix<Double>& elemMat,
                               StdVector<Integer>& eqnVec1,
                               StdVector<Integer>& eqnVec2,
                               FeFctIdType fctId1, FeFctIdType fctId2,
                               bool preventStaticCond,
                               const StdVector<UInt>& sbmIndices,
                               const Double& f,
                               bool isMultHarmDiag)
  {
    // map original matrix destination to analysis-dependent one
    FEMatrixType mappedDest = matrixMap_[dest];

    Matrix<Complex> harmMat;
    // if destination matrix is NOTYPE -> leave
    if( mappedDest == NOTYPE )
      return;

    assert(!elemMat.ContainsNaN() && !elemMat.ContainsInf());

    if( analysisType_ == BasePDE::TRANSIENT || analysisType_ == BasePDE::STATIC || analysisType_ == BasePDE::EIGENFREQUENCY || analysisType_ == BasePDE::BUCKLING || analysisType_ == BasePDE::EIGENVALUE) {
      if ( (analysisType_ == BasePDE::EIGENFREQUENCY || analysisType_ == BasePDE::EIGENVALUE) && (algsys_->IsMatrixComplex()) ) {
        // we have an eigenvalue problem with complex system matrices (e.g. mechanics with complex stiffness tensor)
        Matrix2Harmonic( harmMat, elemMat, STIFFNESS, context.GetEntryType(), 1.0 ); // elemMat -> harmMat with omega=1 and STIFFNESS will convert REAL->COMPLEX
        algsys_->SetElementMatrix( mappedDest, harmMat, fctId1, eqnVec1, fctId2, eqnVec2, context.IsSetCounterPart(), preventStaticCond, context.isDiagonal());
      } else {
        algsys_->SetElementMatrix( mappedDest, elemMat, fctId1, eqnVec1, fctId2, eqnVec2, context.IsSetCounterPart(), preventStaticCond, context.isDiagonal());
      }
    } else {
      assert(analysisType_ == BasePDE::HARMONIC || analysisType_ == BasePDE::MULTIHARMONIC || analysisType_ == BasePDE::INVERSESOURCE || analysisType_ == BasePDE::HARMONIC25D);

      Double omega;
      if(isMultHarmDiag){
        omega = 2 * M_PI * f;
      } else {
        omega = mp_->Eval( mHandle_ );
      }

      if ( analysisType_ == BasePDE::HARMONIC25D) {
        Matrix2Complex( harmMat, elemMat);
      } else {
        Matrix2Harmonic( harmMat, elemMat, dest, context.GetEntryType(), omega );
      }
      
      if( analysisType_ == BasePDE::MULTIHARMONIC){
        algsys_->SetElementMatrix_MultHarm( mappedDest, harmMat,
                                            fctId1, eqnVec1,
                                            fctId2, eqnVec2,
                                            context.IsSetCounterPart(),
                                            sbmIndices);
      }else{
        algsys_->SetElementMatrix( mappedDest, harmMat,
                                  fctId1, eqnVec1,
                                  fctId2, eqnVec2,
                                  context.IsSetCounterPart(),
                                  preventStaticCond,
                                  context.isDiagonal());
      }
    }

  }

  void Assemble::InsertMatrix( FEMatrixType dest, BiLinFormContext& context,
                               Matrix<Complex>& elemMat,
                               StdVector<Integer>& eqnVec1,
                               StdVector<Integer>& eqnVec2,
                               FeFctIdType fctId1, FeFctIdType fctId2,
                               bool preventStaticCond,
                               const StdVector<UInt>& sbmIndices,
                               const Double& f,
                               bool isMultHarmDiag) {
    Matrix<Complex> harmMat;

    // map original matrix destination to analysis-dependent one
    FEMatrixType mappedDest = matrixMap_[dest];

    assert(mappedDest != NOTYPE);
    // bloch mode analysis is complex
    assert(analysisType_ == BasePDE::MULTIHARMONIC || analysisType_ == BasePDE::HARMONIC || analysisType_ == BasePDE::HARMONIC25D
        || analysisType_ == BasePDE::INVERSESOURCE || analysisType_ == BasePDE::EIGENFREQUENCY || analysisType_ == BasePDE::EIGENVALUE);

    assert(!elemMat.ContainsNaN() && !elemMat.ContainsInf());

    Double omega;
    if(isMultHarmDiag){
      omega = 2 * M_PI * f;
    } else {
      omega = mp_->Eval( mHandle_ );
    }

    if(domain->GetDriver()->GetAnalysisType() == BasePDE::HARMONIC || domain->GetDriver()->GetAnalysisType() == BasePDE::INVERSESOURCE ||
       domain->GetDriver()->GetAnalysisType() == BasePDE::MULTIHARMONIC) {
        Matrix2Harmonic( harmMat, elemMat, dest, context.GetEntryType(), omega);
    } else {
      harmMat = elemMat;
    }

    if(analysisType_ == BasePDE::MULTIHARMONIC){
      algsys_->SetElementMatrix_MultHarm( mappedDest, harmMat,
                                          fctId1, eqnVec1,
                                          fctId2, eqnVec2,
                                          context.IsSetCounterPart(),
                                          sbmIndices );
    } else {
      algsys_->SetElementMatrix( mappedDest, harmMat,
                                 fctId1, eqnVec1,
                                 fctId2, eqnVec2,
                                 context.IsSetCounterPart(),
                                 preventStaticCond,
                                 context.isDiagonal() );
    }
  }
}
