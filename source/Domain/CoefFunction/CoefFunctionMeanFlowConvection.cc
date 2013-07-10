#include "CoefFunctionMeanFlowConvection.hh"
#include "FeBasis/FeSpace.hh"

namespace CoupledField{

  CoefFunctionMeanFlowConvection::CoefFunctionMeanFlowConvection(PtrCoefFct density,
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
  
  void CoefFunctionMeanFlowConvection::GetTensor( Matrix<Double>& tensor, 
                                                  const LocPointMapped& lpm ) {

    Vector<Complex> eSolV;
    const Elem* el = lpm.ptEl;
    dynamic_cast<FeFunction<Complex>*>(feFct_.get())->GetElemSolution( eSolV, el );
    BaseFE* ptFe =  feFct_->GetFeSpace()->GetFe(lpm.ptEl->elemNum);
    UInt numDofs = 3;
    StdVector< Vector<Double> > eSolComp(numDofs);
    StdVector< Vector<Double> > VDerivs(numDofs);
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
        eSolComp[rdof][eI] = eSolV[eI*numDofs+rdof].real();
      }
      
      bOperator_->ApplyOp( VDerivs[rdof], lpm, ptFe, eSolComp[rdof] );

      for(UInt dof=0; dof < numDofs; dof++) 
      {
        tensor[rdof][dof] = dens * VDerivs[rdof][dof];
      }
    }
  }
  
}// end of namespace
