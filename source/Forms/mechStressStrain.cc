// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <iostream>
#include <fstream>

#include "mechStressStrain.hh"

namespace CoupledField
{

  // =============================================================================
  // base class for mechanical stress computation
  // =============================================================================

  // base class for calculation of mechanical stresses
  template <class TYPE>
  MechStressStrain<TYPE>::MechStressStrain(BaseMaterial* matData, SubTensorType type) 
    : linElastInt(matData, type)
  {
    name_ = "MechStressStrain";
  }

 
  template <class TYPE>
  MechStressStrain<TYPE>::~MechStressStrain()
  {
  }


  /// calculates of stresses T (vector notation)
  // T = c . S with c the material tensor snd S the linear strains
  // S = B_lin * u 
  // see Habil. M. Kaltenbacher 

  template <class TYPE>
  void MechStressStrain<TYPE>::CalcStressVec(Vector<TYPE>& stressVec, UInt ip, EntityIterator& ent) 
  {
    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent );
    
    Matrix<TYPE> dMat;
    linElastInt::calcDMat(dMat, ent.GetElem()); 

    // convert displacement of all elem nodes into one vector: 
    // (uNode1X, uNode1Y, uNode2X, uNode2Y, ...)
    Vector<TYPE> displVec;
    elemDisp_.ConvertToVec_AppendCols(displVec);
    // linear differential operator B_lin
    Matrix<Double> linBMat;    
    calcBMat( linBMat, ip, ptCoord_);

    Vector<TYPE> linStrain(linBMat.GetSizeRow());
    linStrain.Init();
    stressVec.Init();
    //Matrix<TYPE>(linBMat).Mult(displVec,linStrain);
    //Matrix<TYPE>(dMat).Mult(linStrain,stressVec);
    linStrain = linBMat * displVec;
    stressVec = dMat * linStrain;
  }

  /// calculates green-lagrangian strains (linear part, vector notation)
  // see Habil. M. Kaltenbacher 

  template <class TYPE>
  void MechStressStrain<TYPE>::CalcStrainVec(Vector<TYPE>& strainVec, UInt ip, EntityIterator& ent) 
  {
    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent );

    Matrix<Double> dMat;
    linElastInt::calcDMat(dMat, ent.GetElem());
    
    // convert displacement of all elem nodes into one vector: 
    // (uNode1X, uNode1Y, uNode2X, uNode2Y, ...)
    Vector<TYPE> displVec;
    elemDisp_.ConvertToVec_AppendCols(displVec);

    // linear differential operator B_lin
    Matrix<Double> linBMat;    
    calcBMat( linBMat, ip, ptCoord_);

    strainVec = linBMat * displVec;
  
  }


  // returns linear B - matrix

  template <class TYPE>
  void MechStressStrain<TYPE>::
  calcBMat(Matrix<Double> & bMat, UInt ip, Matrix<Double> & ptCoord)
  {
    // linear differential operator B_lin
    linElastInt::calcBMat(bMat, ip, ptCoord);
  }


#ifdef __GNUC__
  // Explicite template instantiation
  template class MechStressStrain<Double>;
  template class MechStressStrain<Complex>;
#endif

} // end namespace CoupledField
