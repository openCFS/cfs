#include <vector>
#include <string>
#include <fstream>

#include "Elements/elements_header.hh"
#include "General/environment.hh"
#include "elemIntegr.hh"



namespace CoupledField
{


  ElemIntegr::ElemIntegr(UInt elemType) :
    linearLoad_(NULL),
    ptElem_(NULL)
  {
#ifdef TRACE
    (*trace) << " entering ElemIntegr::ElemIntegr " << std::endl;
#endif

    CreatePt2Elems(elemType);

    if(!ptElem_)
      return;

    linearLoad_ = new LinearFlowNoiseInt(ptElem_->ptElem);
  }

  ElemIntegr::~ElemIntegr()
  {
    delete linearLoad_;
    if(ptElem_)
      delete ptElem_->ptElem;
    delete ptElem_;
  }

  void ElemIntegr::CreatePt2Elems( UInt type )
  {
 #ifdef TRACE
    (*trace) << " entering ElemIntegr::CreatePt2Elems " << std::endl;
 #endif
    ptElem_ = new Elem();

    switch(type)
    {
    case Elem::LINE2:
      ptElem_->ptElem = new Line1FE();
      ptElem_->connect = 1, 2;
      ptElem_->edges = 1;
      ptElem_->faces = 1;
      break;

    case Elem::TRIA3:
      ptElem_->ptElem = new Triangle1FE();
      ptElem_->connect = 1, 2, 3;
      ptElem_->edges = 1, 2, 3;
      ptElem_->faces = 1;
      break;

    case Elem::QUAD4:
      ptElem_->ptElem = new Quad1FE();
      ptElem_->connect = 1, 2, 3, 4;
      ptElem_->edges = 1, 2, 3, 4;
      ptElem_->faces = 1;
      break;

    case Elem::QUAD8:
      ptElem_->ptElem = new Quad2FE();
      ptElem_->connect = 1, 5, 2, 6, 3, 7, 4, 8;
      ptElem_->edges = 1, 2, 3, 4;
      ptElem_->faces = 1;
      break;

    case Elem::TET4:
      ptElem_->ptElem = new Tetra1FE();
      ptElem_->connect = 1, 2, 3, 4;
      ptElem_->edges = 1, 2, 3, 4, 5, 6;
      ptElem_->faces = 1, 2, 3;
      break;

    case Elem::HEXA8:
      ptElem_->ptElem = new Hexa1FE();
      ptElem_->connect = 1, 2, 3, 4, 5, 6, 7, 8;
      ptElem_->edges = 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12;
      ptElem_->faces = 1, 2, 3, 4, 5, 6;

      break;

    case Elem::PYRA5:
      ptElem_->ptElem = new Pyra1FE();
      ptElem_->connect = 1, 2, 3, 4, 5;
      ptElem_->edges = 1, 2, 3, 4, 5, 6, 7, 8;
      ptElem_->faces = 1, 2, 3, 4, 5;
      break;

    case Elem::WEDGE6:
      ptElem_->ptElem = new Wedge1FE();
      ptElem_->connect = 1, 2, 3, 4, 5, 6;
      ptElem_->edges = 1, 2, 3, 4, 5, 6, 7, 8, 9, 10;
      ptElem_->faces = 1, 2, 3, 4, 5;
      break;

    default:
      std::cerr << "Element-Type not defined!" << std::endl;
      delete ptElem_;
      ptElem_ = NULL;
      break;
    }

    if(ptElem_)
    {
      ptElem_->regionId = 1;
      ptElem_->elemNum = 1;
    }

  }

  void ElemIntegr::PerformIntegration(const Matrix<Double> & coordMat,
                                      const Matrix<Double>& NodaldTijdxj,
                                      const Matrix<Double>& NodalVal,
                                      Vector<Double>& elemvec,
                                      Vector<Double>& nodalLoadDensity,
                                      Vector<Double>& divLHTensor,
                                      Double density)
  {
#ifdef TRACE
    (*trace) << " entering ElemIntegr::PerformIntegration" << std::endl;
#endif

    if(!ptElem_)
      return;

    linearLoad_->CalcElemVec4QuadwithVel(coordMat,  NodalVal, elemvec, nodalLoadDensity, divLHTensor, ptElem_, density);
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
