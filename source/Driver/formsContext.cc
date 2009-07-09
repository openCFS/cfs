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
    secMatFac_ = "0.0";
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

    dampingLayer_ = NULL;
  }

  BiLinFormContext::~BiLinFormContext() {
    
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
                                  PdeIdType& id1, PdeIdType& id2 ) {

    // Get equation numbers
    map1_->GetEqns( eqnVec1, *result1_, it1 );
    map2_->GetEqns( eqnVec2, *result2_, it2 );

    // Get PDE IDs
    id1 = ptPde1_->GetPDEId();
    id2 = ptPde2_->GetPDEId();

  }
  
  void BiLinFormContext::SetPtPdes(SinglePDE * aPDE1, 
                                   SinglePDE * aPDE2 ) {

    ptPde1_ = aPDE1;
    ptPde2_ = aPDE2;

    map1_ = ptPde1_->GetEqnMap();
    map2_ = ptPde2_->GetEqnMap();
  }
    
  //! Set the result types and entities the bilinearform is working on
  void BiLinFormContext::SetResults( shared_ptr<ResultInfo> result1, 
                                     shared_ptr<ResultInfo> result2,
                                     shared_ptr<EntityList> list1, 
                                     shared_ptr<EntityList> list2 ) {

    result1_ = result1;
    result2_ = result2;
    ent1_ = list1;
    ent2_ = list2;

    assert( integrator_ != NULL);
    assert( result1->fctType != NULL );
    integrator_->SetAnsatzFct( result1->fctType, result2->fctType );
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

  void BiLinFormContext::SetDampLayer(std::string& typeFnc, 
				      Vector<Double>& mPoint, 
				      Double& dampFactor, 
				      Double& dampFactorMax, 
				      Double& startRadius, 
				      Double& endRadius) {

    dampingLayer_ = new DampLayer(typeFnc);
    dampingLayer_->SetDampingParams(mPoint, dampFactor, 
				    dampFactorMax, startRadius, 
				    endRadius);

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
                                   PdeIdType& id ) {
    
    // Get equation numbers
    map_->GetEqns( eqnVec, *result_, it );
    
    // Get Pde Id
    id = ptPde_->GetPDEId();
    
  }
  
  void LinearFormContext::SetPtPde(SinglePDE * ptPde ) {
    
    ptPde_ = ptPde;
    map_ = ptPde_->GetEqnMap();
    

  }

  void LinearFormContext::SetResult( shared_ptr<ResultInfo> result,
                                     shared_ptr<EntityList> list ) {

    result_ = result;
    ent_ = list;
    integrator_->SetAnsatzFct( result->fctType );
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
                                    PdeIdType& id1, PdeIdType& id2 ) {

    const NCElem *elem = dynamic_cast< const NCElem* >(it1.GetElem()); 

    ElemList masterSide(domain->GetGrid() );
    ElemList slaveSide(domain->GetGrid() );
   
    masterSide.SetElement( elem->ptSurfParent );
    slaveSide.SetElement( elem->ptLagrangeParent );

    EntityIterator masterIt, slaveIt;
    masterIt = masterSide.GetIterator();
    slaveIt = slaveSide.GetIterator();

    if (ptPde1_ == ptPde2_) {
      map1_->GetEqns( eqnVec1, *result1_, masterIt );
      map2_->GetEqns( eqnVec2, *result2_, slaveIt );
    } else {
      // check which pde is on master side (via materials)
      std::map<RegionIdType, BaseMaterial*> pdeMat = 
        ptPde1_->getPDEMaterialData();
      std::map<RegionIdType, BaseMaterial*>::iterator itMat =
        pdeMat.find(elem->ptSurfParent->ptVolElem1->regionId);
      if (itMat != pdeMat.end()) { // pde1 is master
        map1_->GetEqns( eqnVec1, *result1_, masterIt );
        map2_->GetEqns( eqnVec2, *result2_, slaveIt );
      } else {
        map1_->GetEqns( eqnVec1, *result1_, slaveIt );
        map2_->GetEqns( eqnVec2, *result2_, masterIt );
      }
    }
    
    // Get PDE IDs
    id1 = ptPde1_->GetPDEId();
    id2 = ptPde2_->GetPDEId();

  }
  

}
