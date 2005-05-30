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
                                         NodeEQN * ptEQN,
                                         NodeStoreSol<TYPE> & potential,
                                         const SolutionType solType,
                                         Boolean isaxi)
    : BaseOperator(ptGrid, ptPDE, ptEQN, isaxi)
  {
    ENTER_FCN( "GradientFieldOp::GradientFieldOp" );  
    this->potential_ = &potential;
    solType_ = solType;
 
  }

  template<class TYPE>
  GradientFieldOp<TYPE>::~GradientFieldOp()
  {
    ENTER_FCN( "GradientFieldOp::~GradientFieldOp" );

  }

  template <class TYPE>
  void GradientFieldOp<TYPE>::CalcElemGradField(CFSVector & elemField,
                                                const Elem * ptElement,
                                                const Vector<Double> & lCoord,
                                                const Double factor)
  {
    ENTER_FCN( "GradientFieldOp::CalcElemGradField" );
  
    UInt dim;
    Vector<TYPE> potEntry(1);
    dim = ptElement->ptElem->GetDim();
 
    Vector<TYPE> & helpElemField = dynamic_cast<Vector<TYPE>&> (elemField);   

    helpElemField.Resize(dim);
    helpElemField.Init(0.0);
   
    UInt nShFnc = 0;
    nShFnc = ptElement->ptElem->GetNumNodes();
 
    const StdVector<UInt> & connect = ptElement->connect;
  
    Matrix<Double> CornerCoords; 
    ptPDE_->GetElemCoords(connect, CornerCoords);

    Matrix<Double> GlobalGradient;

    ptElement->ptElem->GetGlobDerivShFnc(GlobalGradient, lCoord, CornerCoords);

    // loop over shape functions
    for( UInt i=0; i<dim; i++ )
      for( UInt j=0; j<nShFnc; j++ )
        {
          //std::cerr << "Longing for connect = " << connect[j] << std::endl;

          potential_->Get(solType_, connect[j]-1,0, potEntry[0]);

          //istd:cerr << "elecEntry = " << elecEntry << std::endl;
          //E[i] -= GlobalGradient[j][i] * (*EPotential_)[(*ptMesh2PDENode_) [connect[j]-1]-1];

          helpElemField.AddEntry(i,-GlobalGradient[j][i]*potEntry[0]*factor);

        }
  
  }


  template<class TYPE>
  void GradientFieldOp<TYPE>::CalcSDGradField(CFSVector & elemField,
                                              const StdVector<RegionIdType> & SD, 
                                              const Vector<Double> & lCoord,
                                              const Vector<Double> & factors)
  {
    ENTER_FCN( "GradientFieldOp::CalcSDGradField" );
  
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
            nShFnc = subDomain[k]->ptElem->GetNumNodes();
          
            ptPDE_->GetElemCoords( subDomain[k]->connect, cornerCoords );
          
            subDomain[k]->ptElem->GetGlobDerivShFnc(globalGradient, lCoord, cornerCoords);
          
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

  // explicit template instantiation for SGI compiler
#ifdef __sgi
#pragma instantiate GradientFieldOp<Double>
#pragma instantiate GradientFieldOp<Complex>
#endif


} // namespace CoupledField
