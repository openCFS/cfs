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
    ENTER_FCN( "BiLinFormContext::BiLinFormContext" );

    integrator_ = biLinForm;
    
    destMat_ = destMat;
    secDestMat_ = NOTYPE;
    secMatFac_ = "0.0";
    setCounterPart_ = false;
    entryType_ = REAL;

    // Note: By default, we do not set the transposed
    // of a matrix as well, i.e. if an element matrix
    // gets assembled to a main diagonal block within the
    // big FE-Matrix, this flag must not be changed
    // through 'SetCounterPart()'
    setCounterPart_ = false;

    ptPde1_ = NULL;
    ptPde2_ = NULL;

    dampingLayer_ = NULL;
  }

  BiLinFormContext::~BiLinFormContext() {
    ENTER_FCN( "BiLinFormContext::~BiLinFormContext" );
    
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
    ENTER_FCN( "BiLinFormContext::MapEqns" );

    // Get equation numbers
    map1_->GetEqns( eqnVec1, *result1_, it1 );
    map2_->GetEqns( eqnVec2, *result2_, it2 );

    // Get PDE IDs
    id1 = ptPde1_->GetPDEId();
    id2 = ptPde2_->GetPDEId();

  }
  
  void BiLinFormContext::SetPtPdes(SinglePDE * aPDE1, 
                                   SinglePDE * aPDE2 ) {
    ENTER_FCN( "BiLinFormContext::SetPtPdes" );

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
    ENTER_FCN( "BiLinFormContext::SetResults" );

    result1_ = result1;
    result2_ = result2;
    ent1_ = list1;
    ent2_ = list2;

    assert( integrator_ != NULL);
    assert( result1->fctType != NULL );
    integrator_->SetAnsatzFct( result1->fctType, result2->fctType );
  }  

  bool BiLinFormContext::IsNonLin() {
    ENTER_FCN( "BiLinFormContext::IsNonLin" );

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
    ENTER_FCN( "BiLinFormContext::SetDampLayer" );

    dampingLayer_ = new DampLayer(typeFnc);
    dampingLayer_->SetDampingParams(mPoint, dampFactor, 
				    dampFactorMax, startRadius, 
				    endRadius);

  }
  
  void BiLinFormContext::Dump()
  {
      std::cout << "BiLinFormContext nonLin=" << (IsNonLin() ? "true" : "false") 
                << " pde1=" << ptPde1_->GetName() << " pde2=" << (ptPde2_ == NULL ? "NULL" : ptPde2_->GetName()) 
                << " integrator=" << integrator_->GetName() << " entityList1=" << ent1_->GetName() 
                 << std::endl;
  }  
  

// -------------------------------------------------------------------------

  LinearFormContext::LinearFormContext( LinearForm* linearForm ) {
    ENTER_FCN( "LinearFormContext::LinearFormContext" );

    integrator_ = linearForm;

    ptPde_ = NULL;

  }

  LinearFormContext::~LinearFormContext() {
    ENTER_FCN( "LinearFormContext::~LinearFormContext" );

    // delete linearform
    if( integrator_ != NULL ) {
      delete integrator_;
      integrator_ = NULL;
    }
  }


  void LinearFormContext::MapEqns( EntityIterator& it,
                                   StdVector<Integer>& eqnVec,
                                   PdeIdType& id ) {
    ENTER_FCN( "LinearFormContext::MapEqns" );
    
    // Get equation numbers
    map_->GetEqns( eqnVec, *result_, it );
    
    // Get Pde Id
    id = ptPde_->GetPDEId();
    
  }
  
  void LinearFormContext::SetPtPde(SinglePDE * ptPde ) {
    ENTER_FCN( "LinearFormContext::SetPtPde" );
    
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
    ENTER_FCN( "LinearFormContext::IsNonLin" );
   // Return true if linearform is solution-dependent

    if( integrator_->IsSolDependent() ) {
      return true;
    } else {
      return false;
    }
  }

  // -------------------------------------------------------------------------
  NcBiLinFormContext::NcBiLinFormContext( BaseForm* biLinForm, 
                                          FEMatrixType destMat ) 
    : BiLinFormContext( biLinForm, destMat ) {
    ENTER_FCN( "NcBiLinFormContext::BiLinFormContext" );

  }
  
  NcBiLinFormContext::~NcBiLinFormContext() {
    ENTER_FCN( "NcBiLinFormContext::~NcBiLinFormContext" );
    
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
    ENTER_FCN( "NcBiLinFormContext::MapEqns" );

    const NCElem *elem = dynamic_cast< const NCElem* >(it1.GetElem()); 

    ElemList masterSide(domain->GetGrid() );
    ElemList slaveSide(domain->GetGrid() );
   
    masterSide.SetElement( elem->ptSurfParent );
    slaveSide.SetElement( elem->ptLagrangeParent );

    EntityIterator masterIt, slaveIt;
    masterIt = masterSide.GetIterator();
    slaveIt = slaveSide.GetIterator();

    map1_->GetEqns( eqnVec1, *result1_, masterIt );
    map2_->GetEqns( eqnVec2, *result2_, slaveIt );
    
    // Get PDE IDs
    id1 = ptPde1_->GetPDEId();
    id2 = ptPde2_->GetPDEId();

  }
  

 
}
