#ifndef FILE_ELEMINTEGR_2004
#define FILE_ELEMINTEGR_2004

#include <vector>
#include <string>
#include <fstream>

#include <Domain/elem.hh>
#include <General/environment.hh>
#include <Forms/linearForm.hh>
#include "Elements/elements_header.hh"

namespace CoupledField
{


  class ElemIntegr
  {
  public:
    ElemIntegr(UInt elemType);
  
    virtual ~ElemIntegr();
  

    void CreatePt2Elems( UInt elemType );
  
    void PerformIntegration(const Matrix<Double> & coordMat,
                            const Matrix<Double>& NodaldTijdxj,
                            const Matrix<Double>& NodalVal,
                            Vector<Double>& elemvec);

    void ComputeFromCombustionTij(const Matrix<Double> & coordMat,
                                  const Matrix<Double>& NodalVel, 
                                  const Vector<Double>& NodalRho,
                                  Vector<Double>& elemvec);

    void ComputeFromCombustionVector(const Matrix<Double> & coordMat,
                                     const Matrix<Double>& NodalVec,
                                     Vector<Double>& elemvec);

    void ComputeFromCombustionScalar(const Matrix<Double> & coordMat,
                                     const Vector<Double>& NodalScalar,
                                     Vector<Double>& elemvec);

    void ComputeFromCombustionTijOnSurface(const Matrix<Double> & coordMat,
                                           const Matrix<Double>& NodalVel, 
                                           const Vector<Double>& NodalRho,
                                           Vector<Double>& elemvec);

    void ComputeFromCombustionVectorOnSurface(const Matrix<Double> & coordMat,
                                              const Matrix<Double>& NodalVel,
                                              Vector<Double>& elemvec);


  protected:

    Elem *ptElem_;
    LinearFlowNoiseInt * linearLoad_;

  };
}


#endif // FILE_ELEMINTEGR_2004
