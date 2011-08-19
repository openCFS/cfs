// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "formsContext.hh"

#include "Forms/baseForm.hh"
#include "Forms/linearForm.hh"
#include "Domain/entityList.hh"
#include "PDE/SinglePDE.hh"
#include "PDE/eqnMap.hh"
#include "Domain/grid.hh"
#include "Domain/domain.hh"

namespace CoupledField {


  BiLinFormContext::BiLinFormContext( BaseForm* biLinForm, 
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
    

    //ptPde1_ = NULL;
    //ptPde2_ = NULL;

    //dampingLayer_ = NULL;
  }

  BiLinFormContext::~BiLinFormContext() {
    
    // release math parser handle
    mathParser_->ReleaseHandle(secMatFacHandle_);
    
    // delete bilinearform
    if( integrator_ != NULL ) {
      delete integrator_;
      integrator_ = NULL;
    }
    
//    delete dampingLayer_;
//    dampingLayer_ = NULL;
  }
  
  void BiLinFormContext::MapEqns( EntityIterator& it1, 
                                  EntityIterator& it2,
                                  StdVector<Integer>& eqnVec1, 
                                  StdVector<Integer>& eqnVec2,
                                  FeFctIdType& id1, FeFctIdType& id2 ) {

    // Get equation numbers
    feFct1_->GetFeSpace()->GetEqns( eqnVec1, it1 );
    feFct2_->GetFeSpace()->GetEqns( eqnVec2, it2 );

    // Get PDE IDs
    id1 = ptPde1_->GetPDEId();
    id2 = ptPde2_->GetPDEId();

  }
  
  void BiLinFormContext::SetPtPdes(shared_ptr<SinglePDE> aPDE1, 
                                   shared_ptr<SinglePDE> aPDE2 ) {

    //ptPde1_ = aPDE1;
    //ptPde2_ = aPDE2;
    //OBSOLETE
    //map1_ = ptPde1_->GetEqnMap();
    //map2_ = ptPde2_->GetEqnMap();
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
    result1_ = feFct1_->GetResultInfo();
    result2_ = feFct2_->GetResultInfo();
    ptPde1_ = feFct1_->GetPDE();
    ptPde2_ = feFct2_->GetPDE();

    assert( integrator_ != NULL);
    // THIS HAS TO BE REPLACED BY FESPACE
    integrator_->SetFeSpace( fct1->GetFeSpace(), fct2->GetFeSpace());
  }
  
  bool BiLinFormContext::IsNonLin() {

    // Return true if either
    // - bilinearform is solution dependent (true non-linearity)
    // - bilinearform has updated Lagrangian formulation and 
    //   there are coordinate updated present in the grid class
    Grid * grid = domain->GetGrid();
    if( integrator_->IsSolDependent() 
        || (integrator_->IsCoordUpdate()  
            && grid->HasNodalOffset() ) ) {
      return true;
    } else {
      return false;
    }
  }

//  void BiLinFormContext::SetDampLayer(std::string& typeFnc, 
//				      Vector<Double>& mPoint, 
//				      Double& dampFactor, 
//				      Double& dampFactorMax, 
//				      Double& startRadius, 
//				      Double& endRadius) {
//
//    dampingLayer_ = new DampLayer(typeFnc);
//    dampingLayer_->SetDampingParams(mPoint, dampFactor, 
//				    dampFactorMax, startRadius, 
//				    endRadius);
//
//  }
  
  void BiLinFormContext::SetSecDestMat( FEMatrixType aSecMat,
                                        std::string aSecMatFac ) {
    secDestMat_ = aSecMat;
    mathParser_->SetExpr(secMatFacHandle_, aSecMatFac);  
    Double dummy = mathParser_->Eval(secMatFacHandle_);
    dummy +=1.0;
    std::cerr << "Secondary matrix factor is " << mathParser_->GetExpr(secMatFacHandle_) << std::endl;
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
  

// -------------------------------------------------------------------------

  LinearFormContext::LinearFormContext( LinearForm* linearForm ) {

    integrator_ = linearForm;

    //ptPde_ = NULL;

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
    feFct_->GetFeSpace()->GetEqns(eqnVec,it);
    
    // Get Pde Id
    id = ptPde_->GetPDEId();
    
  }
  
  void LinearFormContext::SetPtPde(shared_ptr<SinglePDE> ptPde ) {
    
    ptPde_ = ptPde;
    //map_ = ptPde_->GetEqnMap();
    

  }

  void LinearFormContext::SetEntities( shared_ptr<EntityList> list ) {

    ent_ = list;
  }

  void LinearFormContext::SetFeFunction(shared_ptr<BaseFeFunction>  fct ){
    feFct_ = fct;
    result_ = feFct_->GetResultInfo();
    ptPde_ = feFct_->GetPDE();
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
  
  // -------------------------------------------------------------------------
  NcBiLinFormContext::NcBiLinFormContext( BaseForm* biLinForm, 
                                          FEMatrixType destMat ) 
    : BiLinFormContext( biLinForm, destMat ) {

  }
  
  NcBiLinFormContext::~NcBiLinFormContext() {
    
    // delete bilinearform
    if( integrator_ != NULL ) {
      delete integrator_;
      integrator_ = NULL;
    }
  }

  void NcBiLinFormContext::MapEqns( EntityIterator& it1, 
                                    EntityIterator& it2,
                                    StdVector<Integer>& eqnVec1, 
                                    StdVector<Integer>& eqnVec2,
                                    FeFctIdType& id1, FeFctIdType& id2 ) {

    const NCElem *elem = dynamic_cast< const NCElem* >(it1.GetElem()); 

    ElemList masterSide(domain->GetGrid() );
    ElemList slaveSide(domain->GetGrid() );
   
    masterSide.SetElement( elem->ptSurfParent );
    slaveSide.SetElement( elem->ptLagrangeParent );

    EntityIterator masterIt, slaveIt;
    masterIt = masterSide.GetIterator();
    slaveIt = slaveSide.GetIterator();

    if (ptPde1_ == ptPde2_) {
      feFct1_->GetFeSpace()->GetEqns( eqnVec1,  masterIt );
      feFct2_->GetFeSpace()->GetEqns( eqnVec2, slaveIt );
    } else {
      // check which pde is on master side (via materials)
      std::map<RegionIdType, BaseMaterial*> pdeMat = 
        ptPde1_->getPDEMaterialData();
      std::map<RegionIdType, BaseMaterial*>::iterator itMat =
        pdeMat.find(elem->ptSurfParent->ptVolElem1->regionId);
      if (itMat != pdeMat.end()) { // pde1 is master
        feFct1_->GetFeSpace()->GetEqns( eqnVec1, masterIt );
        feFct2_->GetFeSpace()->GetEqns( eqnVec2, slaveIt );
      } else {
        feFct1_->GetFeSpace()->GetEqns( eqnVec1, slaveIt );
        feFct2_->GetFeSpace()->GetEqns( eqnVec2, masterIt );
      }
    }
    
    // Get PDE IDs
    id1 = ptPde1_->GetPDEId();
    id2 = ptPde2_->GetPDEId();

  }
  

}
