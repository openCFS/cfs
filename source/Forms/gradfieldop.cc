// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "Forms/gradfieldop.hh"

#include <string>

#include "Elements/basefe.hh"
#include "Domain/elem.hh"
#include "Domain/grid.hh"
#include "General/environment.hh"
#include "MatVec/vector.hh"
#include "MatVec/matrix.hh"
#include "DataInOut/Logging/cfslog.hh"
#include "PDE/StdPDE.hh"

namespace CoupledField
{
  DECLARE_LOG(gfo)
  DEFINE_LOG(gfo, "gradFieldOperator")

  template<class TYPE>
  GradientFieldOp<TYPE>::GradientFieldOp(Grid * ptGrid, 
                                         StdPDE * ptPDE,
                                         shared_ptr<EqnMap> eqnMap,
                                         NodeStoreSol<TYPE> & potential,
                                         shared_ptr<AnsatzFct> fctType,
                                         bool isaxi,
                                         bool coordUpdate)
    : BaseOperator(ptGrid, ptPDE, eqnMap, isaxi, coordUpdate )
  {
    this->potential_ = &potential;
    this->fctType_ = fctType;
  }

  template<class TYPE>
  GradientFieldOp<TYPE>::GradientFieldOp(Grid * ptGrid,
                                         StdPDE * ptPDE,
                                         shared_ptr<AnsatzFct> fctType,
                                         bool isaxi,
                                         bool coordUpdate)
    : BaseOperator(ptGrid, ptPDE, ptPDE->GetEqnMap(), isaxi, coordUpdate )
  {
    this->fctType_ = fctType;
  }


  template<class TYPE>
  GradientFieldOp<TYPE>::~GradientFieldOp()
  {

  }

  template <class TYPE>
  void GradientFieldOp<TYPE>::CalcElemGradField(Vector<TYPE> & elemField,
                                                const EntityIterator& ent,
                                                const Vector<Double> & lCoord,
                                                const Double factor,
                                                const UInt dof,
                                                const SingleVector* elem_data)
  {
  
    UInt dim;
    assert(dof >= 1 && dof <= 3);

    dim = ent.GetElem()->ptElem->GetDim();
    elemField.Resize(dim);
    elemField.Init();
   
    UInt nShFnc = 0;
    ent.GetElem()->ptElem->SetAnsatzFct( fctType_ );
    nShFnc = ent.GetElem()->ptElem->GetNumFncs( fctType_ );
    const StdVector<UInt> & connect = ent.GetElem()->connect;
  
    Matrix<Double> CornerCoords; 
    ptGrid_->GetElemNodesCoord( CornerCoords, connect, coordUpdate_ );

    Matrix<Double> GlobalGradient;
    ent.GetElem()->ptElem->GetGlobDerivShFnc(GlobalGradient, lCoord, 
                                         CornerCoords, ent.GetElem() );

    // when the solution is a vector field, it is dim*nShFnc and we access the stride
    Vector<TYPE> potential;
    // optionally it is given
    if(elem_data != NULL)
      potential = *dynamic_cast<const Vector<TYPE>*>(elem_data); // unnecessary copy!
    else
      potential_->GetElemSolution( potential, ent );
    
    const UInt stride = potential.GetSize() / nShFnc;
    assert(stride >= 1 && stride <= 3);
    assert((dof >= 1 && stride >= dof) || (dof == 1 && stride == 1));

    // loop over shape functions
    for( UInt i=0; i<dim; i++ ) {
      for( UInt j=0; j<nShFnc; j++ ) {
        // when stride = 1 we access potential[j], otherwise stride and select the dof
        elemField[i] += -GlobalGradient[j][i] * potential[stride * j + (dof-1)] * factor;
        //LOG_DBG3(gfo) << "GFO::CEGF elem=" << ent.GetElem()->elemNum << " d=" << i << " j=" << j
        //               << " ef[" << i << "]+=" << (-1.0 * GlobalGradient[j][i]) << "*"
        //               << potential[stride * j + (dof-1)] << "*" << factor
        //               << " -> " << elemField[i];
      }
    }
    
  }


  template<class TYPE>
  void GradientFieldOp<TYPE>::CalcSDGradField(Vector<TYPE> & elemField,
                                              const StdVector<RegionIdType> & SD, 
                                              const Vector<Double> & lCoord,
                                              const Vector<Double> & factors)
  {
  
    EXCEPTION( "GradientFieldOp::CalcSDElecField: Not working yet");
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
            nShFnc = subDomain[k]->ptElem->GetNumFncs( fctType_ );
            subDomain[k]->ptElem->SetAnsatzFct( fctType_ );
          
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


} // namespace CoupledField
