#include "CoefFunctionMeanFlowConvection.hh"
#include "FeBasis/FeSpace.hh"

namespace CoupledField{

  template<typename T, UInt DOFS>
  CoefFunctionMeanFlowConvection<T,DOFS>::CoefFunctionMeanFlowConvection(PtrCoefFct density,
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
  
  template<typename T, UInt DOFS>
  void CoefFunctionMeanFlowConvection<T,DOFS>::GetTensor( Matrix<T>& tensor, 
                                                          const LocPointMapped& lpm ) {

    Vector<T> eSolV;
    const Elem* el = lpm.ptEl;
    dynamic_cast<FeFunction<T>*>(feFct_.get())->GetElemSolution( eSolV, el );
    BaseFE* ptFe =  feFct_->GetFeSpace()->GetFe(lpm.ptEl->elemNum);
    StdVector< Vector<T> > eSolComp(DOFS);
    StdVector< Vector<T> > VDerivs(DOFS);
    UInt eN=eSolV.GetSize()/DOFS;
    tensor.Resize(DOFS, DOFS);
    tensor.Init();
    Double dens = 0;
    density_->GetScalar(dens,lpm);


    for(UInt rdof=0; rdof < DOFS; rdof++) 
    {      
      eSolComp[rdof].Resize(eN);
      
      for(UInt eI=0; eI < eN; eI++) 
      {
        eSolComp[rdof][eI] = eSolV[eI*DOFS+rdof];
      }
      
      bOperator_->ApplyOp( VDerivs[rdof], lpm, ptFe, eSolComp[rdof] );

      for(UInt dof=0; dof < DOFS; dof++) 
      {
        tensor[dof][rdof] = dens * VDerivs[rdof][dof];
      }
    }
  }

#ifdef EXPLICIT_TEMPLATE_INSTANTIATION
  template class CoefFunctionMeanFlowConvection<Double,2>;
  template class CoefFunctionMeanFlowConvection<Complex,2>;

  template class CoefFunctionMeanFlowConvection<Double,3>;
  template class CoefFunctionMeanFlowConvection<Complex,3>;
#endif
}// end of namespace

