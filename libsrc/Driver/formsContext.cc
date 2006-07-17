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
    secMatFac_ = 0.0;
    setCounterPart_ = true;
    matDataType_ = REAL;

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
  void BiLinFormContext::SetResults( shared_ptr<ResultDof> result1, 
                                     shared_ptr<ResultDof> result2,
                                     shared_ptr<EntityList> list1, 
                                     shared_ptr<EntityList> list2 ) {
    ENTER_FCN( "BiLinFormContext::SetResults" );

    result1_ = result1;
    result2_ = result2;
    ent1_ = list1;
    ent2_ = list2;
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

// -------------------------------------------------------------------------

  LinearFormContext::LinearFormContext( LinearForm* linearForm, 
                                        Double phase,
                                        const std::string &dynamics ) {
    ENTER_FCN( "LinearFormContext::LinearFormContext" );

    integrator_ = linearForm;

    ptPde_ = NULL;
    phase_ = phase;
    dynamics_ = dynamics;
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

  bool LinearFormContext::IsNonLin() {
    ENTER_FCN( "LinearFormContext::IsNonLin" );
   // Return true if linearform is solution-dependent

    if( integrator_->IsSolDependent() ) {
      return true;
    } else {
      return false;
    }
  }


 
}
