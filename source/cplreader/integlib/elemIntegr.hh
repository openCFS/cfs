#ifndef FILE_ELEMINTEGR_2004
#define FILE_ELEMINTEGR_2004

#include <vector>
#include <string>
#include <fstream>

#include "Domain/elem.hh"
#include "General/environment.hh"
#include "Forms/linearForm.hh"
#include "Elements/elements_header.hh"

namespace CoupledField
{


  class ElemIntegr
  {
  public:
    ElemIntegr(UInt elemType);

    ElemIntegr();

    ElemIntegr & operator=( ElemIntegr &rhs);

    ElemIntegr & operator=(const ElemIntegr &rhs);
  
    virtual ~ElemIntegr();
  
    UInt GetElemType() const;

    void CreatePt2Elems( UInt elemType );
  
    void PerformIntegration(const Matrix<Double> & coordMat,
                            const Matrix<Double>& NodaldTijdxj,
                            const Matrix<Double>& NodalVal,
                            Vector<Double>& elemvec,
                            Vector<Double>& nodalLoadDensity,
                            Vector<Double>& divLHTensor,
                            Double density,
                            bool surfInt);
    
    /**
     * same as PerformIntegration except it takes only the velocity at the
     * centre of the element and not at each node. The advantage is, that it may
     * cancel out numerical errors, and showed improvements with quadrativ
     * elements.
     * @param coordMat The coordinates of each node
     * @param NodalVel Nodal velocity, a vector a size of "node numbers" times
     * dimension
     * @param resVec The vector in which the result is stored
     * @param nodalLoadDensity The result normed to the volume
     * @param divLHTensor The divergence of the Lighthill tensor
     * @param surfInt true in case of surface integral computation
     */
    void PerformIntegrationCentre(const Matrix<Double> & coordMat,
                                  const Matrix<Double>& NodalVel,
                                  Vector<Double>& resVec,
                                  Vector<Double>& nodalLoadDensity,
                                  Vector<Double>& divLHTensor,
                                  Double density,
                                  bool surfInt);

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

    LinearFlowNoiseInt * linearLoad_;
    Elem *ptElem_;
    UInt elemType_;

  };

  //! This struct provides the funcionality to perform source calculation in parallel
  //! Main feature to accomplish this is a copy constructor accessible by Open MP
  //! For bigger problems, this overhead is not sgnificant
  //! The motivation is within the pointer type of the original map<UInt,ElemIntegr*>
  //! Therefor this map is created temporarily and stores the dereferenced pointers
  //! The classes can then be copied by openmp
  struct IntegrationMap{
    public:
    IntegrationMap(std::map<UInt, ElemIntegr *> inMap){
      std::map<UInt, ElemIntegr *>::iterator it = inMap.begin();
      while(it!=inMap.end()){
        integMap_[it->first] = *it->second;
        it++;
      }
    }

    IntegrationMap( IntegrationMap & toCopy){
      std::map<UInt, ElemIntegr>::iterator it = toCopy.integMap_.begin();
      while(it!=toCopy.integMap_.end()){
        integMap_[it->first] = it->second;
        it++;
      }
    }

    ElemIntegr & operator[](UInt key){
      return integMap_[key];
    }

    std::map<UInt, ElemIntegr> integMap_;
  };
}


#endif // FILE_ELEMINTEGR_2004
