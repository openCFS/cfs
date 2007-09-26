// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "Forms/gradfieldop.hh"

#include <string>

#include "Elements/basefe.hh"
#include "Domain/elem.hh"
#include "Domain/grid.hh"
#include "General/environment.hh"
#include "Utils/vector.hh"
#include "Matrix/matrix.hh"

#include <PDE/StdPDE.hh>

namespace CoupledField
{

  template<class TYPE>
  GradientFieldOp<TYPE>::GradientFieldOp(Grid * ptGrid, 
                                         StdPDE * ptPDE,
                                         shared_ptr<EqnMap> eqnMap,
                                         NodeStoreSol<TYPE> & potential,
                                         const SolutionType solType,
                                         shared_ptr<ResultInfo> result,
                                         bool isaxi,
                                         bool coordUpdate)
    : BaseOperator(ptGrid, ptPDE, eqnMap, isaxi, coordUpdate )
  {
    this->potential_ = &potential;
    solType_ = solType;
    result_ = result;
 
  }

  template<class TYPE>
  GradientFieldOp<TYPE>::~GradientFieldOp()
  {

  }

  template <class TYPE>
  void GradientFieldOp<TYPE>::CalcElemGradField(Vector<TYPE> & elemField,
                                                const EntityIterator& ent,
                                                const Vector<Double> & lCoord,
                                                const Double factor)
  {
  
    UInt dim;

    dim = ent.GetElem()->ptElem->GetDim();
    elemField.Resize(dim);
    elemField.Init(0.0);
   
    UInt nShFnc = 0;
    ent.GetElem()->ptElem->SetAnsatzFct( result_->fctType );
    nShFnc = ent.GetElem()->ptElem->GetNumFncs( result_->fctType );
    const StdVector<UInt> & connect = ent.GetElem()->connect;
  
    Matrix<Double> CornerCoords; 
    ptGrid_->GetElemNodesCoord( CornerCoords, connect, coordUpdate_ );

    Matrix<Double> GlobalGradient;
    ent.GetElem()->ptElem->GetGlobDerivShFnc(GlobalGradient, lCoord, 
                                         CornerCoords, ent.GetElem() );

    Vector<TYPE> potential;
    potential_->GetElemSolution( potential, ent );
    
    // loop over shape functions
    for( UInt i=0; i<dim; i++ ) {
      for( UInt j=0; j<nShFnc; j++ ) {
        elemField[i] += -GlobalGradient[j][i] * potential[j] * factor;
      }
    }
    
  }


  template<class TYPE>
  void GradientFieldOp<TYPE>::CalcSDGradField(Vector<TYPE> & elemField,
                                              const StdVector<RegionIdType> & SD, 
                                              const Vector<Double> & lCoord,
                                              const Vector<Double> & factors)
  {
  
    Error( "GradientFieldOp::CalcSDElecField: Not working yet", __FILE__, __LINE__);
    UInt nShFnc = 0;
    UInt dim;
    Matrix<Double> cornerCoords;
    Matrix<Double> globalGradient;
  
    StdVector<Elem *> subDomain;
    UInt maxelem;
    maxelem = ptGrid_->GetNumElems(SD);
    dim = ptGrid_->GetDim();

    elemField.Resize(maxelem * dim);
    elemField.Init( (TYPE) 0.0);
    //elemField.SetNumSolutions(1);
    //elemField.SetNumNodes(maxelem);
    //elemField.SetNumDofs(dim);
            
    // Iterate over all subdomains
    for( UInt iSD=0; iSD<SD.GetSize(); iSD++)
      {
        ptGrid_->GetVolElems(subDomain,SD[iSD]);

        // Iterate over whole SubDomain
        for( UInt k=0; k<subDomain.GetSize(); k++) 
          {
            nShFnc = subDomain[k]->ptElem->GetNumFncs( result_->fctType );
            subDomain[k]->ptElem->SetAnsatzFct( result_->fctType );
          
            ptGrid_->GetElemNodesCoord( cornerCoords,subDomain[k]->connect, 
                                        coordUpdate_ );
          
            subDomain[k]->ptElem-> 
              GetGlobDerivShFnc(globalGradient, lCoord, 
                                cornerCoords, subDomain[k] );
          
            // loop over shape functions
            for( UInt i=0; i<dim; i++ )
              for( UInt j=0; j<nShFnc; j++ )
                {
                  // NOT WORKING YET
                  //elecEntry = (*EPotential_)(SubDomain[k]->connect[j],1);
                  ///E(k,i) -= GlobalGradient[j][i] * (*EPotential_)[(*ptMesh2PDENode_) [SubDomain[k]->connect[j]-1]-1];      
                  //E(k,i) -= GlobalGradient[j][i] * elecEntry;
                }
          
          }
      }
  }

  // explicit teplate instantiation for gcc compiler
#ifdef __GNUC__
template class GradientFieldOp<Double>;
template class GradientFieldOp<Complex>;
#endif

  // explicit template instantiation for SGI compiler
#ifdef __sgi
#pragma instantiate GradientFieldOp<Double>
#pragma instantiate GradientFieldOp<Complex>
#endif


} // namespace CoupledField
