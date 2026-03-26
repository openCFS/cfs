// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     AeroacousticBase.hh
 *       \brief    <Description>
 *
 *       \date     Jan 7, 2017
 *       \author   kroppert
 */
//================================================================================================

#pragma once


#include <Filters/MeshFilter.hh>
#include "DataInOut/SimInput.hh"


namespace CFSDat{

class AeroacousticBase : public MeshFilter{

public:

  AeroacousticBase(UInt numWorkers, CF::PtrParamNode config, shared_ptr<ResultManager> resMan);

  virtual ~AeroacousticBase();


protected:

  virtual void PrepareCalculation() = 0;

  virtual ResultIdList SetUpstreamResults() = 0;

  virtual void AdaptFilterResults() = 0;



  void OmegaVectorProductU(const Vector<Double>& inVecVel,
                           const Vector<Double>& inVecOmega,
                           Vector<Double>& retVec,
                           const UInt& dim);


  void ScalarProduct(Vector<Double>& retVec,
                    const Vector<Double>& inVec1,
                    const Vector<Double>& inVec2,
                    const UInt& numEquPerEnt,
                    const Double& scalarFactor);

  void TensorProduct(Vector<Double>& retVec,
                    const Vector<Double>& tensor1,
                    const Vector<Double>& tensor2,
                    const UInt& numEquPerEnt,
                    const Double& scalarFactor);

  void TensorProduct(Vector<Double>& retVec,
                    const Vector<Double>& tensor1,
                    const Vector<Double>& tensor2,
                    const UInt& numEquPerEnt,
                    const Double& scalarFactor,
                    const Vector<Double>& scalar);


  void Node2Cell(Vector<Double>& returnVec,
        const Vector<Double>& inVec,
        const UInt& numEquPerEnt,
        const StdVector< StdVector<CF::UInt> >& targetSourceIndex,
        const UInt& maxNumTrgEntities,
        const UInt& gridDim);


  Vector<Double> ScalarToTwoD(const Vector<Double>& inVec);

  Vector<Double> TwoDToScalar(const Vector<Double>& inVec);

};



}
