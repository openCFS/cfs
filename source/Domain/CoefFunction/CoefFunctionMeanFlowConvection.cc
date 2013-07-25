#include "CoefFunctionMeanFlowConvection.hh"
#include "FeBasis/FeSpace.hh"

namespace CoupledField{

  template<typename T>
  CoefFunctionMeanFlowConvection<T>::CoefFunctionMeanFlowConvection(PtrCoefFct density,
                                                                    PtrCoefFct viscosity,
                                                                    BaseBOperator* opt,
                                                                    shared_ptr<BaseFeFunction> feFnc)
    : CoefFunction(),
      density_(density),  
      viscosity_(viscosity),
      bOperator_(opt),
      feFct_(feFnc)
  {
    dimType_ = TENSOR;
  }
  
  template<typename T>
  void CoefFunctionMeanFlowConvection<T>::GetTensor( Matrix<T>& tensor, 
                                                     const LocPointMapped& lpm ) {

    Vector<T> eSolV;
    const Elem* el = lpm.ptEl;
    dynamic_cast<FeFunction<T>*>(feFct_.get())->GetElemSolution( eSolV, el );
    BaseFE* ptFe =  feFct_->GetFeSpace()->GetFe(lpm.ptEl->elemNum);
    UInt numDofs = feFct_->GetVecSize();
    StdVector< Vector<T> > eSolComp(numDofs);
    StdVector< Vector<T> > VDerivs(numDofs);
    UInt eN=eSolV.GetSize()/numDofs;
    tensor.Resize(numDofs, numDofs);
    tensor.Init();
    Double dens = 0;
    density_->GetScalar(dens,lpm);


    for(UInt rdof=0; rdof < numDofs; rdof++) 
    {      
      eSolComp[rdof].Resize(eN);
      
      for(UInt eI=0; eI < eN; eI++) 
      {
        eSolComp[rdof][eI] = eSolV[eI*numDofs+rdof];
      }
      
      bOperator_->ApplyOp( VDerivs[rdof], lpm, ptFe, eSolComp[rdof] );

      for(UInt dof=0; dof < numDofs; dof++) 
      {
        tensor[rdof][dof] = dens * VDerivs[rdof][dof];
      }
    }
  }

#ifdef EXPLICIT_TEMPLATE_INSTANTIATION
  template class CoefFunctionMeanFlowConvection<Double>;
  template class CoefFunctionMeanFlowConvection<Complex>;
#endif
}// end of namespace

