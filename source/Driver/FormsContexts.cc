#include "FormsContexts.hh"

#include "Utils/mathParser/mathParser.hh"
#include "Domain/ElemMapping/EntityLists.hh"
#include "PDE/SinglePDE.hh"
#include "Domain/Mesh/Grid.hh"
#include "Domain/Domain.hh"
#include "Forms/LinForms/LinearForm.hh"
#include "Forms/BiLinForms/BiLinearForm.hh"

namespace CoupledField {

Enum<BiLinearForm::Type> BiLinearForm::type;


  BiLinFormContext::BiLinFormContext( BiLinearForm* biLinForm,
                                      FEMatrixType destMat) {

    integrator_ = biLinForm;

    destMat_ = destMat;
    secDestMat_ = NOTYPE;
    mathParser_ = domain->GetMathParser();
    secMatFacHandle_ = mathParser_->GetNewHandle();
    
    //    setCounterPart_ = false;
    entryType_ = Global::REAL;

    // Note: By default, we do not set the counter part
    // of a matrix as well, i.e. if an element matrix
    // gets assembled to a main diagonal block within the
    // big FE-Matrix, this flag must not be changed
    // through 'SetCounterPart()'
    setCounterPart_ = false;
    negateEntries_ = false;
    
    ptPde1_ = NULL;
    ptPde2_ = NULL;

    useVolEqnA_ = biLinForm->GetUseVolEqnA();
    useVolEqnB_ = biLinForm->GetUseVolEqnB();
  }

  BiLinFormContext::~BiLinFormContext() {
    
    // release math parser handle
    mathParser_->ReleaseHandle(secMatFacHandle_);
    
    // delete bilinearform
    if( integrator_ != NULL ) {
      delete integrator_;
      integrator_ = NULL;
    }
    
  }
  
  void BiLinFormContext::SetEnums()
  {
    BiLinearForm::type.SetName("BiLinearForm::Type");
    BiLinearForm::type.Add(BiLinearForm::NO_BILIN_TYPE, "no bilinear form type");
    BiLinearForm::type.Add(BiLinearForm::BILIN_WRAPPED_LIN, "BiLinWarappenLinForm");
    BiLinearForm::type.Add(BiLinearForm::BDB_INT, "BDBInt");
    BiLinearForm::type.Add(BiLinearForm::BB_INT, "BBInt");
    BiLinearForm::type.Add(BiLinearForm::AB_INT, "ABInt");
    BiLinearForm::type.Add(BiLinearForm::ADB_INT, "ABDInt");
    BiLinearForm::type.Add(BiLinearForm::SINGLE_ENTRY_BILIN_INT, "SingleEntryBiLinInt");
    BiLinearForm::type.Add(BiLinearForm::IC_MODES_INT, "ICModesInt");
  }

  void BiLinFormContext::MapEqns( EntityIterator& it1, 
                                  EntityIterator& it2,
                                  StdVector<Integer>& eqnVec1, 
                                  StdVector<Integer>& eqnVec2,
                                  FeFctIdType& id1, FeFctIdType& id2 ) {

    if( useVolEqnA_ ){
      // Use eqn evaltion from volume (use case: eval surface BLF where volume information is needed --> e.g. gradient)

      // get volume element 
      const SurfElem* sElem1 = it1.GetSurfElem();
      Elem* ptVolElem1 = sElem1->ptVolElems[0];

      this->feFct1_.lock()->GetFeSpace()->GetElemEqns(eqnVec1,ptVolElem1);
       
      id1 = feFct1_.lock()->GetFctId();

    } else {
      // Standard case

      // Get equation numbers
      feFct1_.lock()->GetFeSpace()->GetEqns( eqnVec1, it1 );

      // Get PDE IDs
      id1 = feFct1_.lock()->GetFctId();
    }

    if( useVolEqnB_ ){
      // Use eqn evaltion from volume (use case: eval surface BLF where volume information is needed --> e.g. gradient)

      // get volume element 
      const SurfElem* sElem2 = it2.GetSurfElem();
      Elem* ptVolElem2 = sElem2->ptVolElems[0];

      this->feFct2_.lock()->GetFeSpace()->GetElemEqns(eqnVec2,ptVolElem2);
       
      id2 = feFct2_.lock()->GetFctId();

    } else {
      // Standard case

      // Get equation numbers
      feFct2_.lock()->GetFeSpace()->GetEqns( eqnVec2, it2 );

      // Get PDE IDs
      id2 = feFct2_.lock()->GetFctId();
    }

  }
    
  //! Set the result types and entities the bilinearform is working on
  void BiLinFormContext::SetEntities( shared_ptr<EntityList> list1, 
                                     shared_ptr<EntityList> list2 ) {

    ent1_ = list1;
    ent2_ = list2;
  }  

  void BiLinFormContext::SetFeFunctions( shared_ptr<BaseFeFunction> fct1, 
                                         shared_ptr<BaseFeFunction> fct2) {
    feFct1_ = fct1;
    feFct2_ = fct2;  
    result1_ = feFct1_.lock()->GetResultInfo();
    result2_ = feFct2_.lock()->GetResultInfo();
    ptPde1_ = feFct1_.lock()->GetPDE();
    ptPde2_ = feFct2_.lock()->GetPDE();

    assert( integrator_ != NULL);
    // THIS HAS TO BE REPLACED BY FESPACE
    integrator_->SetFeSpace( fct1->GetFeSpace(), fct2->GetFeSpace());
  }
  
  bool BiLinFormContext::IsNonLin() {

    // Return true if either
    // - bilinearform is solution dependent (true non-linearity)
    // - bilinearform has updated Lagrangian formulation and 
    //   there are coordinate updated present in the grid class
    
    Grid * grid = ent1_->GetGrid();
    if( integrator_->IsSolDependent() 
        || (integrator_->IsCoordUpdate()  
            && grid->HasNodalOffset() ) ) {
      return true;
    } else {
      return false;
    }
  }

  void BiLinFormContext::SetSecDestMat( FEMatrixType aSecMat,
                                        std::string aSecMatFac ) {
    secDestMat_ = aSecMat;
    mathParser_->SetExpr(secMatFacHandle_, aSecMatFac);  
    mathParser_->Eval(secMatFacHandle_);
  }
  std::string BiLinFormContext::GetSecMatFac() const {
    return mathParser_->GetExpr(secMatFacHandle_);
  }
  
  Double BiLinFormContext::EvalSecMatFac() const {
    return mathParser_->Eval(secMatFacHandle_);
  }
  
  
  
  std::string BiLinFormContext::ToString()
  {
    std::ostringstream os;
    os << "BiLinFormContext nonLin=" << (IsNonLin() ? "true" : "false") << " timeFreq=" << (integrator_->IsTimeFrequencyDependent() ? "true" : "false")
       << " pde1=" << ptPde1_->GetName() << " pde2=" << (ptPde2_ == NULL ? "NULL" : ptPde2_->GetName()) 
       << " integrator=" << integrator_->GetName() << " entityList1=" << ent1_->GetName(); 
    return os.str();
  }  
  




// -------------------------------------------------------------------------

  LinearFormContext::LinearFormContext( LinearForm* linearForm ) {

    integrator_ = linearForm;

    ptPde_ = NULL;

  }

  LinearFormContext::~LinearFormContext() {

    // delete linearform
    if( integrator_ != NULL ) {
      delete integrator_;
      integrator_ = NULL;
    }
  }


  void LinearFormContext::MapEqns( EntityIterator& it,
                                   StdVector<Integer>& eqnVec,
                                   FeFctIdType& id ) {
    
    // Get equation numbers
    //map_->GetEqns( eqnVec, *result_, it );
    feFct_.lock()->GetFeSpace()->GetEqns(eqnVec,it);
    
    // Get Pde Id
    id = feFct_.lock()->GetFctId();
    
  }
  
  void LinearFormContext::SetPtPde(SinglePDE* ptPde ) {
    
    ptPde_ = ptPde;
    //map_ = ptPde_->GetEqnMap();
    

  }

  void LinearFormContext::SetEntities( shared_ptr<EntityList> list ) {

    ent_ = list;
  }

  void LinearFormContext::SetFeFunction(shared_ptr<BaseFeFunction>  fct ){
    feFct_ = fct;
    result_ = feFct_.lock()->GetResultInfo();
    ptPde_ = feFct_.lock()->GetPDE();
    integrator_->SetFeSpace( fct->GetFeSpace() );
  }

  bool LinearFormContext::IsNonLin() {
   // Return true if linearform is solution-dependent
    if( integrator_->IsSolDependent() ) {
      return true;
    } else {
      return false;
    }
  }


  std::string LinearFormContext::ToString() const
  {
    std::ostringstream os;
    os << "integrator: " << integrator_->GetName() 
       << " pde: " << ptPde_->GetName() << " at=" <<BasePDE::analysisType.ToString(ptPde_->GetAnalysisType());

    return os.str(); 
  }
  
  /***************************************************************************
   * NcBiLinFormContext
   **************************************************************************/
  
  void NcBiLinFormContext::MapEqns( EntityIterator &it1, EntityIterator &it2,
                                    StdVector<Integer> &eqnVec1,
                                    StdVector<Integer> &eqnVec2,
                                    FeFctIdType &id1, FeFctIdType &id2)
  {


    const NcSurfElem* ncElem1 = it1.GetNcSurfElem();
#ifndef NDEBUG // Paranoia pays...
    const NcSurfElem* ncElem2 = it2.GetNcSurfElem();
    
    // NcBiLinFormContext requires identical EntityIterators
    assert( ncElem1 == ncElem2 );
#endif
    
    const MortarNcSurfElem *mortarElem =
        dynamic_cast<const MortarNcSurfElem*>(ncElem1);
    // NcBiLinFormContext only works with MortarNcSurfElems at the moment
    assert( mortarElem );
    
    shared_ptr<BaseFeFunction> feFunc1 = this->feFct1_.lock(),
                               feFunc2 = this->feFct2_.lock();
    SolutionType solType1 = feFunc1->GetResultInfo()->resultType,
                 solType2 = feFunc2->GetResultInfo()->resultType;
    
    if ( (solType1 == LAGRANGE_MULT) && (solType2 != LAGRANGE_MULT) ) {
      feFunc1->GetFeSpace()->GetElemEqns(eqnVec1, mortarElem->ptSlave);
      feFunc2->GetFeSpace()->GetElemEqns(eqnVec2, mortarElem->ptMaster);
    }
    else if ( (solType1 != LAGRANGE_MULT) && (solType2 == LAGRANGE_MULT) ) {
      feFunc1->GetFeSpace()->GetElemEqns(eqnVec1, mortarElem->ptMaster);
      feFunc2->GetFeSpace()->GetElemEqns(eqnVec2, mortarElem->ptSlave);
    }
    else if ( (solType1 != LAGRANGE_MULT) && (solType2 != LAGRANGE_MULT) ) {
      // no Lagrange multiplier => probably acou-mech acoupling
      
      SurfElem *surfElem1 = NULL, *surfElem2 = NULL;
      const std::set<RegionIdType> &regions1 = feFunc1->GetRegions(),
                                   &regions2 = feFunc2->GetRegions();
      
      if (regions1.find(mortarElem->ptMaster->ptVolElems[0]->regionId)
          != regions1.end()) {
        surfElem1 = mortarElem->ptMaster;
      }
      else if (regions1.find(mortarElem->ptSlave->ptVolElems[0]->regionId)
          != regions1.end()) {
        surfElem1 = mortarElem->ptSlave;
      }
      else {
        EXCEPTION("None of the parents of NcSurfElem " << mortarElem->elemNum
            << "were found in the regions of FeFunction '"
            << SolutionTypeEnum.ToString(solType1) << "'");
      }
      
      if (regions2.find(mortarElem->ptMaster->ptVolElems[0]->regionId)
          != regions2.end()) {
        surfElem2 = mortarElem->ptMaster;
      }
      else if (regions2.find(mortarElem->ptSlave->ptVolElems[0]->regionId)
          != regions2.end()) {
        surfElem2 = mortarElem->ptSlave;
      }
      else {
        EXCEPTION("None of the parents of NcSurfElem " << mortarElem->elemNum
            << "were found in the regions of FeFunction '"
            << SolutionTypeEnum.ToString(solType2) << "'");
      }

      feFunc1->GetFeSpace()->GetElemEqns(eqnVec1, surfElem1);
      feFunc2->GetFeSpace()->GetElemEqns(eqnVec2, surfElem2);
    }
    else {
      // speaking of paranoia ...
      EXCEPTION("You cannot couple two Lagrange multipliers");
    }
    
    id1 = feFunc1->GetFctId();
    id2 = feFunc2->GetFctId();
  }
  
  void NcBiLinFormContext::GetEqns( StdVector<Integer>& eqnVec1,
                StdVector<Integer>& eqnVec2,
                FeFctIdType& id1, FeFctIdType& id2 ) const
  {

    EntityIterator it = ent1_->GetIterator();
    assert(it.GetType() == EntityList::NC_ELEM_LIST);
    assert(ent2_->GetIterator().GetType() == EntityList::NC_ELEM_LIST);

    shared_ptr<ElemList> masterElems(new ElemList(ent1_->GetGrid()));
    shared_ptr<ElemList> slaveElems(new ElemList(ent1_->GetGrid()));
    
    /*for ( ; !it.IsEnd(); it++ ) {
      const NcSurfElem* ncElem = it.GetNcSurfElem();
      const MortarNcSurfElem* mortarElem =
          dynamic_cast<const MortarNcSurfElem*>(ncElem);
      assert(mortarElem);
      masterElems->AddElement(mortarElem->ptMaster);
      slaveElems->AddElement(mortarElem->ptSlave);
    }*/
    const NcSurfElem* ncElem = it.GetNcSurfElem();
    const MortarNcSurfElem* mortarElem =
        dynamic_cast<const MortarNcSurfElem*>(ncElem);
    assert(mortarElem);
    masterElems->SetRegion(mortarElem->ptMaster->regionId);
    slaveElems->SetRegion(mortarElem->ptSlave->regionId);
    
    shared_ptr<BaseFeFunction> feFunc1 = this->feFct1_.lock(),
                               feFunc2 = this->feFct2_.lock();
    SolutionType solType1 = feFunc1->GetResultInfo()->resultType,
                 solType2 = feFunc2->GetResultInfo()->resultType;
    
    if ( (solType1 == LAGRANGE_MULT) && (solType2 != LAGRANGE_MULT) ) {
      feFunc1->GetFeSpace()->GetEntityListEqns(eqnVec1, slaveElems);
      feFunc2->GetFeSpace()->GetEntityListEqns(eqnVec2, masterElems);
    }
    else if ( (solType1 != LAGRANGE_MULT) && (solType2 == LAGRANGE_MULT) ) {
      feFunc1->GetFeSpace()->GetEntityListEqns(eqnVec1, masterElems);
      feFunc2->GetFeSpace()->GetEntityListEqns(eqnVec2, slaveElems);
    }
    else if ( (solType1 != LAGRANGE_MULT) && (solType2 != LAGRANGE_MULT) ) {
      // no Lagrange multiplier => probably acou-mech acoupling
      
      shared_ptr<ElemList> elemList1, elemList2;
      const std::set<RegionIdType> &regions1 = feFunc1->GetRegions(),
                                   &regions2 = feFunc2->GetRegions();
      
      if (regions1.find(mortarElem->ptMaster->ptVolElems[0]->regionId)
          != regions1.end()) {
        elemList1 = masterElems;
      }
      else if (regions1.find(mortarElem->ptSlave->ptVolElems[0]->regionId)
          != regions1.end()) {
        elemList1 = slaveElems;
      }
      else {
        EXCEPTION("None of the parents of NcSurfElem " << mortarElem->elemNum
            << "were found in the regions of FeFunction '"
            << SolutionTypeEnum.ToString(solType1) << "'");
      }
      
      if (regions2.find(mortarElem->ptMaster->ptVolElems[0]->regionId)
          != regions2.end()) {
        elemList2 = masterElems;
      }
      else if (regions2.find(mortarElem->ptSlave->ptVolElems[0]->regionId)
          != regions2.end()) {
        elemList2 = slaveElems;
      }
      else {
        EXCEPTION("None of the parents of NcSurfElem " << mortarElem->elemNum
            << "were found in the regions of FeFunction '"
            << SolutionTypeEnum.ToString(solType2) << "'");
      }

      feFunc1->GetFeSpace()->GetEntityListEqns(eqnVec1, elemList1);
      feFunc2->GetFeSpace()->GetEntityListEqns(eqnVec2, elemList2);
    }
    else {
      // speaking of paranoia ...
      EXCEPTION("You cannot couple two Lagrange multipliers");
    }
    
    id1 = feFunc1->GetFctId();
    id2 = feFunc2->GetFctId();
  }

  SurfaceBiLinFormContext::SurfaceBiLinFormContext( BiLinearForm* biLinForm, FEMatrixType destMat, BiLinearForm::CouplingDirection currentDirection  )
                              : NcBiLinFormContext(biLinForm,destMat){

     this->currentDirection_ = currentDirection;
   }

   //! Destructor
   SurfaceBiLinFormContext::~SurfaceBiLinFormContext(){

   }


   void SurfaceBiLinFormContext::MapEqns( EntityIterator& it1,
                             EntityIterator& it2,
                             StdVector<Integer>& eqnVec1,
                             StdVector<Integer>& eqnVec2,
                             FeFctIdType& id1, FeFctIdType& id2 ){
     //this is more or less a switch case statement which terms to consider and where to aquire the
     //equation numbers
     const SurfElem* sElem1 = it1.GetSurfElem();
     const SurfElem* sElem2 = it2.GetSurfElem();
     //just to be sure, this context requires it1 and it2 to be identical
     if(sElem1->elemNum != sElem2->elemNum){
       EXCEPTION("SurfaceBiLinFormContext requires identical iterators")
     }
     //now lets try to downcast to MortarNcSurfElem
     const MortarNcSurfElem * mSe1 = dynamic_cast<const MortarNcSurfElem*>(sElem1);

     if(mSe1->ptMaster->ptVolElems[1] != NULL){
       EXCEPTION("Master surface element has not exactly one neighbor... can this be true??")
     }
     if(mSe1->ptSlave->ptVolElems[1] != NULL){
       EXCEPTION("Master surface element has not exactly one neighbor... can this be true??")
     }
     Elem* volEMaster = mSe1->ptMaster->ptVolElems[0];
     Elem* volESlave  = mSe1->ptSlave->ptVolElems[0];

     // Hamideh: In the case of having two different PDEs id1 and id2 could be different and should be assign
     // according to the coupling direction.
     // id1 and eqnVec1 --> Test function
     // id2 and eqnVec2 --> Unknown
     // Note: In the case of LinFlowMechNitsche coupling feFct1_ is the unknown from Master(Displacement-Mechanics PDE)
     // and feFct2_ is the unknown from slave(velocity-LinFlow PDE)

     switch(currentDirection_){
     case BiLinearForm::MASTER_MASTER:
       this->feFct1_.lock()->GetFeSpace()->GetElemEqns(eqnVec1,volEMaster);
       eqnVec2 = eqnVec1;
       id1 = feFct1_.lock()->GetFctId();
       id2 = feFct1_.lock()->GetFctId();
       break;
     case BiLinearForm::SLAVE_SLAVE:
       this->feFct2_.lock()->GetFeSpace()->GetElemEqns(eqnVec1,volESlave);
       eqnVec2 = eqnVec1;
       id1 = feFct2_.lock()->GetFctId();
       id2 = feFct2_.lock()->GetFctId();
       break;
     case BiLinearForm::MASTER_SLAVE:
       this->feFct1_.lock()->GetFeSpace()->GetElemEqns(eqnVec1,volEMaster);
       this->feFct2_.lock()->GetFeSpace()->GetElemEqns(eqnVec2,volESlave);
       id1 = feFct1_.lock()->GetFctId();
       id2 = feFct2_.lock()->GetFctId();
       break;
     case BiLinearForm::SLAVE_MASTER:
       this->feFct2_.lock()->GetFeSpace()->GetElemEqns(eqnVec1,volESlave);
       this->feFct1_.lock()->GetFeSpace()->GetElemEqns(eqnVec2,volEMaster);
       id1 = feFct2_.lock()->GetFctId();
       id2 = feFct1_.lock()->GetFctId();
       break;
     default:
       EXCEPTION("Undefined coupling direction");
     }
   }

   void SurfaceBiLinFormContext::GetEqns( StdVector<Integer>& eqnVec1,
                                               StdVector<Integer>& eqnVec2,
                                               FeFctIdType& id1, FeFctIdType& id2 ) const {

     EntityIterator it = ent1_->GetIterator();
     eqnVec1.Clear();
     eqnVec2.Clear();
     shared_ptr<SurfElemList> masterElems(new SurfElemList(ent1_->GetGrid()));
     shared_ptr<SurfElemList> slaveElems(new SurfElemList(ent1_->GetGrid()));

     const NcSurfElem* ncElem = it.GetNcSurfElem();
     const MortarNcSurfElem* mortarElem =
         dynamic_cast<const MortarNcSurfElem*>(ncElem);
     assert(mortarElem);

     masterElems->SetRegion(mortarElem->ptMaster->regionId);
     slaveElems->SetRegion(mortarElem->ptSlave->regionId);

     //now we gather the volume elements
     shared_ptr<ElemList> masterVolElems(new ElemList(ent1_->GetGrid()));
     shared_ptr<ElemList> slaveVolElems(new ElemList(ent1_->GetGrid()));

     EntityIterator slaveSurfsIt = slaveElems->GetIterator();
     slaveSurfsIt.Begin();
     for( ; !slaveSurfsIt.IsEnd(); slaveSurfsIt++ ) {
       const SurfElem* sElem1 = slaveSurfsIt.GetSurfElem();
       Elem* curElem = sElem1->ptVolElems[0];
       slaveVolElems->AddElement(curElem);
     }

     EntityIterator masterSurfsIt = masterElems->GetIterator();
     masterSurfsIt.Begin();
     for( ; !masterSurfsIt.IsEnd(); masterSurfsIt++ ) {
       const SurfElem* sElem1 = masterSurfsIt.GetSurfElem();
       Elem* curElem = sElem1->ptVolElems[0];
       masterVolElems->AddElement(curElem);
     }

     //now we gather all volume elements associated
     switch(currentDirection_){
     case BiLinearForm::MASTER_MASTER:
       this->feFct1_.lock()->GetFeSpace()->GetEntityListEqns(eqnVec1,masterVolElems);
       eqnVec2 = eqnVec1;
       break;
     case BiLinearForm::SLAVE_SLAVE:
       this->feFct2_.lock()->GetFeSpace()->GetEntityListEqns(eqnVec1,slaveVolElems);
       eqnVec2 = eqnVec1;
       break;
     case BiLinearForm::MASTER_SLAVE:
       this->feFct1_.lock()->GetFeSpace()->GetEntityListEqns(eqnVec1,masterVolElems);
       this->feFct2_.lock()->GetFeSpace()->GetEntityListEqns(eqnVec2,slaveVolElems);
       break;
     case BiLinearForm::SLAVE_MASTER:
       this->feFct1_.lock()->GetFeSpace()->GetEntityListEqns(eqnVec1,slaveVolElems);
       this->feFct2_.lock()->GetFeSpace()->GetEntityListEqns(eqnVec2,masterVolElems);
       break;
     default:
       EXCEPTION("Undefined coupling direction");
     }
     id1 = this->feFct1_.lock()->GetFctId();
     id2 = this->feFct2_.lock()->GetFctId();
   }


} // namespace CoupledField
