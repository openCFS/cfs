#include "CoefFunctionMeanFlowConvection.hh"
#include "FeBasis/FeSpace.hh"

#include "DataInOut/Logging/LogConfigurator.hh"

namespace CoupledField{
DEFINE_LOG(coefMeanFlowConv, "CoefFunctionMeanFlowConvection")

  template<typename T, UInt DOFS>
  CoefFunctionMeanFlowConvection<T,DOFS>::CoefFunctionMeanFlowConvection(PtrCoefFct scalarCoef,
                                                                         BaseBOperator* opt,
                                                                         shared_ptr<BaseFeFunction> feFnc)
    : CoefFunction(),
      factor_(scalarCoef),
      bOperator_(opt),
      feFct_(feFnc)
  {
    dimType_ = TENSOR;
    dim_ = DOFS;
  }

  template<typename T, UInt DOFS>
  CoefFunctionMeanFlowConvection<T,DOFS>::CoefFunctionMeanFlowConvection(BaseBOperator* opt, shared_ptr<BaseFeFunction> feFnc)
    : CoefFunction(),
      bOperator_(opt),
      feFct_(feFnc)
  {
    dimType_ = VECTOR;
    dim_ = DOFS;
    factor_ = NULL;
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
    UInt eN=eSolV.GetSize()/DOFS; // number of nodes
    tensor.Resize(DOFS, DOFS);
    tensor.Init();
    T fact;
    if (factor_) {
      factor_->GetScalar(fact,lpm);
    } else {
      fact = T(1.0);
    }

    LOG_DBG(coefMeanFlowConv) << el->elemNum << ":\n";
    for(UInt rdof=0; rdof < DOFS; rdof++) 
    {
      eSolComp[rdof].Resize(eN);

      for(UInt eI=0; eI < eN; eI++) 
      {
        eSolComp[rdof][eI] = eSolV[eI*DOFS+rdof];// extract the x,y,z components of the nodal values
      }
      LOG_DBG(coefMeanFlowConv) << " eSolComp = " << eSolComp[rdof].ToString() << "\n";
      bOperator_->ApplyOp( VDerivs[rdof], lpm, ptFe, eSolComp[rdof] );
      LOG_DBG(coefMeanFlowConv) << " VDerivs = " << VDerivs.ToString() << "\n";

      for(UInt dof=0; dof < DOFS; dof++) 
      {
        tensor[dof][rdof] = fact * VDerivs[rdof][dof];
      }
    }
    LOG_DBG(coefMeanFlowConv) <<" tensor = " << tensor.ToString() << "\n";
  }


  template<typename T, UInt DOFS>
  void CoefFunctionMeanFlowConvection<T,DOFS>::GetVector( Vector<T>& tensor, const LocPointMapped& lpm ) {
    Matrix<T> gradV;
    this->GetTensor(gradV,lpm);
    Vector<T> v;
    feFct_->GetVector(v,lpm);
    tensor.Resize(dim_);
    gradV.Mult(v, tensor);
    // tensor = gradV * v;
  }


  template class CoefFunctionMeanFlowConvection<Double,2>;
  template class CoefFunctionMeanFlowConvection<Complex,2>;

  template class CoefFunctionMeanFlowConvection<Double,3>;
  template class CoefFunctionMeanFlowConvection<Complex,3>;
}// end of namespace

