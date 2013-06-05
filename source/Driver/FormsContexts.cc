#include "FormsContexts.hh"

#include "Domain/ElemMapping/EntityLists.hh"
#include "PDE/SinglePDE.hh"
#include "Domain/Mesh/Grid.hh"
#include "Domain/Domain.hh"
#include "Forms/LinForms/LinearForm.hh"
#include "Forms/BiLinForms/BiLinearForm.hh"

namespace CoupledField {


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
  
  void BiLinFormContext::MapEqns( EntityIterator& it1, 
                                  EntityIterator& it2,
                                  StdVector<Integer>& eqnVec1, 
                                  StdVector<Integer>& eqnVec2,
                                  FeFctIdType& id1, FeFctIdType& id2 ) {

    // Get equation numbers
    feFct1_.lock()->GetFeSpace()->GetEqns( eqnVec1, it1 );
    feFct2_.lock()->GetFeSpace()->GetEqns( eqnVec2, it2 );

    // Get PDE IDs
    id1 = feFct1_.lock()->GetFctId();
    id2 = feFct2_.lock()->GetFctId();

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
    Double dummy = mathParser_->Eval(secMatFacHandle_);
    dummy +=1.0;
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
    os << "BiLinFormContext nonLin=" << (IsNonLin() ? "true" : "false") 
       << " pde1=" << ptPde1_->GetName() << " pde2=" << (ptPde2_ == NULL ? "NULL" : ptPde2_->GetName()) 
       << " integrator=" << integrator_->GetName() << " entityList1=" << ent1_->GetName(); 
    return os.str();
  }  
  


  SurfaceBiLinFormContext::SurfaceBiLinFormContext( BiLinearForm* biLinForm, FEMatrixType destMat, BiLinearForm::CouplingDirection currentDirection  )
                             : BiLinFormContext(biLinForm,destMat){

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

    switch(currentDirection_){
    case BiLinearForm::MASTER_MASTER:
      this->feFct1_.lock()->GetFeSpace()->GetElemEqns(eqnVec1,volEMaster);
      eqnVec2 = eqnVec1;
      break;
    case BiLinearForm::SLAVE_SLAVE:
      this->feFct1_.lock()->GetFeSpace()->GetElemEqns(eqnVec1,volESlave);
      eqnVec2 = eqnVec1;
      break;
    case BiLinearForm::MASTER_SLAVE:
      this->feFct1_.lock()->GetFeSpace()->GetElemEqns(eqnVec1,volEMaster);
      this->feFct2_.lock()->GetFeSpace()->GetElemEqns(eqnVec2,volESlave);
      break;
    case BiLinearForm::SLAVE_MASTER:
      this->feFct1_.lock()->GetFeSpace()->GetElemEqns(eqnVec1,volESlave);
      this->feFct2_.lock()->GetFeSpace()->GetElemEqns(eqnVec2,volEMaster);
      break;
    default:
      EXCEPTION("Undefined coupling direction");
    }
    // Get PDE IDs
    id1 = feFct1_.lock()->GetFctId();
    id2 = feFct2_.lock()->GetFctId();
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
       << " pde: " << ptPde_->GetName();

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
    const NcSurfElem* ncElem2 = it2.GetNcSurfElem();
    
    // NcBiLinFormContext requires identical EntityIterators
    assert( ncElem1 == ncElem2 );
    
    const MortarNcSurfElem *mortarElem =
        dynamic_cast<const MortarNcSurfElem*>(ncElem1);
    // NcBiLinFormContext only works with MortarNcSurfElems at the moment
    assert( mortarElem );
    
    // TODO jens: implement the general case for two different FeFunctions
    // (e.g. MechAcou coupling)
    shared_ptr<BaseFeFunction> feFuncField = this->feFct1_.lock(),
                               feFuncLM = this->feFct2_.lock();
    if ( feFuncField->GetResultInfo()->resultType == LAGRANGE_MULT ) {
      feFuncField = this->feFct2_.lock();
      feFuncLM = this->feFct1_.lock();
      if ( feFuncField->GetResultInfo()->resultType == LAGRANGE_MULT ) {
        EXCEPTION("You cannot couple two Lagrange multipliers");
      }
    }
    feFuncField->GetFeSpace()->GetElemEqns(eqnVec1, mortarElem->ptMaster);
    feFuncLM->GetFeSpace()->GetElemEqns(eqnVec2, mortarElem->ptSlave);
    
    id1 = feFuncField->GetFctId();
    id2 = feFuncLM->GetFctId();
  }
  
} // namespace CoupledField
