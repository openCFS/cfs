// =====================================================================================
// 
//       Filename:  coefFunction.hh
// 
//    Description:  This is the base class for every coefficient
//                  It is used in the Integrator classes to obtain the material tensor
//                  (D-Matrix) as well as in the B-Operator classes to realize e.g. a
//                  change in the derivatives. 
//
//                  Example1: consider the PML
//                  Coefficient function provides a material matrix multiplied by the
//                  complex jacobian. Moreover the same coeffient function is passed to
//                  the operator to provide the daming parameters for the "streched"
//                  derivatives.
//
//                  Example2: Consider geometric non-linear mechanics
//                  The CoefFunction can provide the solution at the current integration
//                  point to enable the b-operator to compute its matrix of derivatives
//
//                  This class is just the interface class right now, the following
//                  structure should be realized:
//                  - coefFunctionAnalytic
//                    This class descibes the quantity by an analytic expression this
//                    includes (non-linear) materials, analytic flow fields, etc.
//                  - coefFunctionGrid
//                    This class provides the Information based on information read in
//                    from a grid. E.g. in the case of (interpolated) aeroacoustic source terms 
//                  - coefFucnctionSol
//                    This function provides the information based on the current
//                    solution of the problem
//
//
// 
//        Version:  1.0
//        Created:  10/19/2011 10:16:34 AM
//       Revision:  none
//       Compiler:  g++
// 
//         Author:  Andreas Hueppe (AHU), andreas.hueppe@uni-klu.ac.at
//        Company:  Universitaet Klagenfurt
// 
// =====================================================================================

#ifndef COEFFUNCTION_HH
#define COEFFUNCTION_HH

#include "General/environment.hh"
#include "MatVec/matrix.hh"
#include "Domain/elem.hh"
#include "Elements/elemShapeMap.hh"

namespace CoupledField{

class CoefFunction{
  public:

    typedef enum{ UNDEF,SCALAR,VECTOR,TENSOR } CoefType;

    CoefFunction(){
      mType_ = UNDEF;
    }

    ~CoefFunction(){
      ;
    }

    virtual void GetTensor(Matrix<Double>& CoefMat, LocPointMapped lp, const Elem* elem){
      EXCEPTION("CoefFunction::GetTensor<Double> called: This may not happen");
    }

    virtual void GetVector(Vector<Double>& CoefMat, LocPointMapped lp, const Elem* elem){
      EXCEPTION("CoefFunction::GetVector<Double> called: This may not happen");
    }

    virtual void GetScalar(Double& CoefMat, LocPointMapped lp, const Elem* elem){
      EXCEPTION("CoefFunction::GetScalar<Double> called: This may not happen");
    }

    virtual void GetTensor(Matrix<Complex>& CoefMat, LocPointMapped lp, const Elem* elem){
      EXCEPTION("CoefFunction::GetTensor<Complex> called: This may not happen");
    }

    virtual void GetVector(Vector<Complex>& CoefMat, LocPointMapped lp, const Elem* elem){
      EXCEPTION("CoefFunction::GetVector<Complex> called: This may not happen");
    }

    virtual void GetScalar(Complex& CoefMat, LocPointMapped lp, const Elem* elem){
      EXCEPTION("CoefFunction::GetScalar<Complex> called: This may not happen");
    }


    void SetCoordinateSystem(shared_ptr<CoordSystem> cSys){
      coordSys_ = cSys;
    }

    virtual CoefType GetType(){
      return mType_;
    }

    //ElemList GetSupport(){
    //  return ElementList();
    //}

  protected:

    shared_ptr<CoordSystem> coordSys_;

    CoefType mType_;


};

}
#endif
