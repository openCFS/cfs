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
//                  Coeficient function provides a material matrix multiplied by the
//                  complex jacobian. Moreover the same coeffient function is passed to
//                  the operator to provide the daming paraeters for the "streched"
//                  derivatives.
//
//                  Example2: Consider geomatric non-linear mechanics
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
    CoefFunction(){
      ;
    }

    ~CoefFunction(){
      ;
    }

    virtual void GetMatrix(Matrix<Complex>& CoefMat, LocPointMapped lp, const Elem* elem)=0;

    virtual void GetMatrix(Matrix<Double>& CoefMat, LocPointMapped lp, const Elem* elem)=0;

    void SetCoordinateSystem(shared_ptr<CoordSystem> cSys){
      coordSys_ = cSys;
    }

    virtual void SetConstMatrix(Matrix<Double>& CoefMat)=0;

    virtual void SetConstMatrix(Matrix<Complex>& CoefMat)=0;

    //ElemList GetSupport(){
    //  return ElementList();
    //}

  protected:

    shared_ptr<CoordSystem> coordSys_;


};

}
#endif
