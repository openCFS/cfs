#include <vector>
#include <string>
#include <fstream>

#include "Elements/elements_header.hh"
#include "General/environment.hh"
#include "elemIntegr.hh"



namespace CoupledField
{


  ElemIntegr::ElemIntegr() :
    linearLoad_(NULL),
    ptElem_(NULL),
    elemType_(0),
    volElem_(0)
  { }

  ElemIntegr::ElemIntegr(UInt elemType) :
    linearLoad_(NULL),
    ptElem_(NULL),
    elemType_(0),
    volElem_(0)
  {
#ifdef TRACE
    (*trace) << " entering ElemIntegr::ElemIntegr " << std::endl;
#endif

    ptElem_ = CreatePt2Elems(elemType);

    elemType_ = elemType;
    if(!ptElem_)
      return;

    linearLoad_ = new LinearFlowNoiseInt(ptElem_->ptElem);
  }

  ElemIntegr::~ElemIntegr()
  {
    delete linearLoad_;

    if(volElem_){
      delete volElem_->ptElem;
      delete volElem_;
    }

    if(ptElem_)
      delete ptElem_->ptElem;
    delete ptElem_;

  }

  ElemIntegr & ElemIntegr::operator=( ElemIntegr &rhs){
    if (this != &rhs) // protect against invalid self-assignment
    {
      this->elemType_ = rhs.GetElemType();
      delete ptElem_;

      this->ptElem_ = CreatePt2Elems(this->elemType_);
      if(!this->ptElem_)
        return *this;
      if(linearLoad_)
        delete linearLoad_;

      linearLoad_ = new LinearFlowNoiseInt(ptElem_->ptElem);
    }
    return *this;
  }

  ElemIntegr & ElemIntegr::operator=(const ElemIntegr &rhs){
    if (this != &rhs) // protect against invalid self-assignment
    {
      this->elemType_ = rhs.GetElemType();
      delete ptElem_;

      this->ptElem_ = CreatePt2Elems(this->elemType_);
      if(!this->ptElem_)
        return *this;
      if(linearLoad_)
        delete linearLoad_;

      linearLoad_ = new LinearFlowNoiseInt(ptElem_->ptElem);
    }
    std::cout << "created" << std::endl;
    return *this;
  }

  UInt ElemIntegr::GetElemType() const{
    return elemType_;
  }

  Elem* ElemIntegr::CreatePt2Elems( UInt type)
  {
   #ifdef TRACE
      (*trace) << " entering ElemIntegr::CreatePt2Elems " << std::endl;
   #endif
    Elem * elem = new Elem();

    switch(type)
    {
    case Elem::LINE2:
      elem->ptElem = new Line1FE();
      elem->connect = 1, 2;
      elem->edges = 1;
      elem->faces = 1;
      break;

    case Elem::TRIA3:
      elem->ptElem = new Triangle1FE();
      elem->connect = 1, 2, 3;
      elem->edges = 1, 2, 3;
      elem->faces = 1;
      break;

    case Elem::QUAD4:
      elem->ptElem = new Quad1FE();
      elem->connect = 1, 2, 3, 4;
      elem->edges = 1, 2, 3, 4;
      elem->faces = 1;
      break;

    case Elem::QUAD8:
      elem->ptElem = new Quad2FE();
      elem->connect = 1, 5, 2, 6, 3, 7, 4, 8;
      elem->edges = 1, 2, 3, 4;
      elem->faces = 1;
      break;

    case Elem::TET4:
      elem->ptElem = new Tetra1FE();
      elem->connect = 1, 2, 3, 4;
      elem->edges = 1, 2, 3, 4, 5, 6;
      elem->faces = 1, 2, 3;
      break;

    case Elem::HEXA8:
      elem->ptElem = new Hexa1FE();
      elem->connect = 1, 2, 3, 4, 5, 6, 7, 8;
      elem->edges = 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12;
      elem->faces = 1, 2, 3, 4, 5, 6;

      break;

    case Elem::PYRA5:
      elem->ptElem = new Pyra1FE();
      elem->connect = 1, 2, 3, 4, 5;
      elem->edges = 1, 2, 3, 4, 5, 6, 7, 8;
      elem->faces = 1, 2, 3, 4, 5;
      break;

    case Elem::WEDGE6:
      elem->ptElem = new Wedge1FE();
      elem->connect = 1, 2, 3, 4, 5, 6;
      elem->edges = 1, 2, 3, 4, 5, 6, 7, 8, 9, 10;
      elem->faces = 1, 2, 3, 4, 5;
      break;

    default:
      std::cerr << "Element-Type not defined!" << std::endl;
      delete elem;
      elem = NULL;
      break;
    }

    if(elem)
    {
      elem->regionId = 1;
      elem->elemNum = 1;
    }

    return elem;
  }

  void ElemIntegr::PerformIntegrationLighthill(const Matrix<Double> & coordMat,
                          const Matrix<Double>& NodaldTijdxj,
                          const Matrix<Double>& NodalVal,
                          Vector<Double>& elemvec,
                          Vector<Double>& nodalLoadDensity,
                          Vector<Double>& divLHTensor,
                          Double density){
    if(!ptElem_)
      return;

    //perform volume or surface (surfInt = true) integration
      linearLoad_->CalcElemVec4QuadwithVel(coordMat,  NodalVal, elemvec,
                                           nodalLoadDensity, divLHTensor,
                                           ptElem_, density);

  }
  
  void ElemIntegr::PerformIntegrationPresD2(const Matrix<Double> & coordMat,
                                                 const Vector<Double>& NodalPresD2,
                                                 Vector<Double>& elemvec,
                                                 Vector<Double>& nodalLoadDensity){
    if(!ptElem_)
      return;

    //perform volume or surface (surfInt = true) integration
      linearLoad_->CalcElemVecWaveWithPressD2(coordMat,  NodalPresD2, elemvec,
                                           nodalLoadDensity, ptElem_);

  }

    void ElemIntegr::PerformIntegrationLighthillwithDivTij(const Matrix<Double> & coordMat,
                          const Matrix<Double>& NodaldTijdxj,
                          const Matrix<Double>& NodalVal,
                          Vector<Double>& elemvec,
                          Vector<Double>& nodalLoadDensity,
                          Vector<Double>& divLHTensor,
                          Double density){
    if(!ptElem_)
      return;

    //perform volume or surface (surfInt = true) integration
      linearLoad_->CalcElemVec4QuadwithDivTij(coordMat,  NodaldTijdxj, elemvec,
                                           nodalLoadDensity, divLHTensor,
                                           ptElem_, density);

   }

  void ElemIntegr::PerformIntegrationLHPressure(const Matrix<Double> & coordMat,
                                                const Vector<Double>& NodalPres,
                                                Vector<Double>& elemVecPres,
                                                Vector<Double>& nodalLoadDensity){
    if(!ptElem_)
      return;

      linearLoad_->CalcElemVecLHwithPress(coordMat,  NodalPres, elemVecPres,nodalLoadDensity,
                                          ptElem_);

  }

  void ElemIntegr::PerformIntegrationAPEMass(const Matrix<Double> & coordMat,
                          const Matrix<Double>& NodalVal,
                          const Vector<Double>& NodalPres,
                          const Vector<Double>& nodalPressureTDeriv,
                          Vector<Double>& elemVecPres,
                          Double density){
    if(!ptElem_)
      return;

    //perform integration for pressure values
    linearLoad_->CalcElemVec4QuadwithPress(coordMat,
                              NodalPres,
                              nodalPressureTDeriv,
                              NodalVal,
                              elemVecPres,
                              ptElem_,
                              density);

  }

  void ElemIntegr::PerformIntegrationAPEMomentum(const Matrix<Double> & coordMat,
                          const Matrix<Double>& NodalVal,
                          const Matrix<Double>& nodalMeanVel,
                          Vector<Double>& elemVecLamb,
                          Vector<Double>& elemVecLambRhs,
                          Double density){
    if(!ptElem_)
      return;

    //perform integration for lamb vector
    linearLoad_->CalcElemVecWithLamb(coordMat,
                                     NodalVal,
                                     nodalMeanVel,
                                     elemVecLambRhs,
                                     elemVecLamb,
                                     ptElem_,
                                     density);

  }

  void ElemIntegr::PerformIntegrationMechRhs(const Matrix<Double> & coordMat,
                          const Matrix<Double>& NodalForce,
                          Vector<Double>& elemvec){

    if(!ptElem_)
      return;

    //perform integration for lamb vector
    linearLoad_->CalcElemVecSurfForce(coordMat,
        NodalForce,
        elemvec,ptElem_);

  }

  void ElemIntegr::PerformIntegrationAeroAcouSrc(const Matrix<Double> & coordMat,
                          const Matrix<Double>& NodalVal,
                          const Matrix<Double>& nodalMeanVel,
                          Vector<Double>& elemVecAeroAcou){
    if(!ptElem_)
      return;

    linearLoad_->CalcElemVec4AeroAcouSrc(coordMat,
                                         NodalVal,
                                         nodalMeanVel,
                                         elemVecAeroAcou,
                                         ptElem_);

  }

  void ElemIntegr::PerformSurfaceIntegration(const UInt volElemType,
                          const Matrix<Double>& ptVolCoord,
                          const Matrix<Double>& ptSurfCoord,
                          const Matrix<Double> & volumeVel,
                          Vector<Double> & surfNormal,
                          Double density,
                          Vector<Double> & Result,
                          Vector<Double> & ResultLHTens){
#ifdef TRACE
    (*trace) << " entering ElemIntegr::PerformSurfaceIntegration" << std::endl;
#endif

    if(!volElem_){
      volElem_ = CreatePt2Elems( volElemType );
    }else{
      if(volElem_->ptElem->feType() != (Integer)volElemType){
        delete volElem_->ptElem;
        volElem_->ptElem = 0;
        delete volElem_;
        volElem_=NULL;
        volElem_ = CreatePt2Elems( volElemType );
      }
    }

    if(!volElem_)
     Exception("ElemIntegr::PerformSurfaceIntegration error creating volume element pointer.");

    linearLoad_->CalcLighthillSurfaceTermVel(volElem_,ptElem_,
        ptVolCoord,ptSurfCoord,volumeVel,surfNormal,density,Result,ResultLHTens);
  }

  void ElemIntegr::PerformIntegrationCentre(const Matrix<Double> & coordMat,
                                            const Matrix<Double>& NodalVal,
                                            const Vector<Double>& NodalPres,
                                            const Matrix<Double>& nodalMeanVel,
                                            const Vector<Double>& nodalMeanPressure,
                                            const Vector<Double>& nodalPressureTDeriv,
                                            Vector<Double>& elemvec,
                                            Vector<Double>& nodalLoadDensity,
                                            Vector<Double>& divLHTensor,
                                            Vector<Double>& elemVecLamb,
                                            Vector<Double>& elemVecLambRhs,
                                            Vector<Double>& elemVecPres,
                                            Vector<Double>& elemVecAeroAcou,
                                            Double density)
  {
#ifdef TRACE
    (*trace) << " entering ElemIntegr::PerformIntegrationCentre" << std::endl;
#endif

    if(!ptElem_)
      return;

    linearLoad_->CalcElemVec4QuadwithVelCentre(coordMat, 
                                               NodalVal, elemvec,
                                               nodalLoadDensity, divLHTensor, 
                                               ptElem_, density);
  }


  void ElemIntegr::ComputeFromCombustionTij(const Matrix<Double> & coordMat,
                                            const Matrix<Double>& NodalVel,
                                            const Vector<Double>& NodalRho,
                                            Vector<Double>& elemvec) {
#ifdef TRACE
    (*trace) << " entering ElemIntegr::ComputeFromCombustionTij" << std::endl;
#endif

    linearLoad_->CalcElemVec4CombustionTij(coordMat,
                                           NodalVel,
                                           NodalRho,
                                           elemvec,
                                           ptElem_);
  }


  void ElemIntegr::ComputeFromCombustionVector(const Matrix<Double> & coordMat,
                                               const Matrix<Double>& NodalVec,
                                               Vector<Double>& elemvec)  {
#ifdef TRACE
    (*trace) << " entering ElemIntegr::ComputeFromCombustionVector" << std::endl;
#endif

    linearLoad_->CalcElemVec4CombustionVec(coordMat,
                                           NodalVec,
                                           elemvec,
                                           ptElem_);
  }


  void ElemIntegr::ComputeFromCombustionScalar(const Matrix<Double> & coordMat,
                                               const Vector<Double>& NodalScalar,
                                               Vector<Double>& elemvec)  {
#ifdef TRACE
    (*trace) << " entering ElemIntegr::ComputeFromCombustionScalar" << std::endl;
#endif

    linearLoad_->CalcElemVec4CombustionScalar(coordMat,
                                              NodalScalar,
                                              elemvec,
                                              ptElem_);
  }


  void ElemIntegr::ComputeFromCombustionTijOnSurface(const Matrix<Double> & coordMat,
                                                     const Matrix<Double>& NodalVel,
                                                     const Vector<Double>& NodalRho,
                                                     Vector<Double>& elemvec) {
#ifdef TRACE
    (*trace) << " entering ElemIntegr::ComputeFromCombustionTijOnSurface" << std::endl;
#endif

    linearLoad_->CalcElemVec4CombustionTijOnSurface(coordMat,
                                                    NodalVel,
                                                    NodalRho,
                                                    elemvec,
                                                    ptElem_);
  }

  void ElemIntegr::ComputeFromCombustionVectorOnSurface(const Matrix<Double> & coordMat,
                                                        const Matrix<Double>& NodalVel,
                                                        Vector<Double>& elemvec) {
#ifdef TRACE
    (*trace) << " entering ElemIntegr::ComputeFromCombustionVectorOnSurface" << std::endl;
#endif

    linearLoad_->CalcElemVec4CombustionVectorOnSurface(coordMat,
                                                       NodalVel,
                                                       elemvec,
                                                       ptElem_);
  }

}
