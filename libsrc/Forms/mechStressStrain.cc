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
  MechStressStrain<TYPE>::MechStressStrain(BaseFE * aptelem, BaseMaterial* matData,
					   SubTensorType type) 
    : linElastInt(aptelem, matData, type)

  {
    ENTER_FCN( "MechStressStrain::MechStressStrain" );
  }


  template <class TYPE>
  MechStressStrain<TYPE>::MechStressStrain(BaseMaterial* matData, SubTensorType type) 
    : linElastInt(matData, type)
  {
    ENTER_FCN( "MechStressStrain::MechStressStrain" );
  }
 
  template <class TYPE>
  MechStressStrain<TYPE>::~MechStressStrain()
  {
    ENTER_FCN( "MechStressStrain::~MechStressStrain" );
  }


  /// calculates of stresses T (vector notation)
  // T = c . S with c the material tensor snd S the linear strains
  // S = B_lin * u 
  // see Habil. M. Kaltenbacher 

  template <class TYPE>
  void MechStressStrain<TYPE>::
  CalcStressVec(Vector<TYPE>& stressVec, UInt ip, Matrix<Double> & ptCoord)
  {
    ENTER_FCN( "MechStressStrain::calcPiolaStressTensor" );

    Matrix<TYPE> dMat;
    linElastInt::calcDMat(dMat);

    // convert displacement of all elem nodes into one vector: 
    // (uNode1X, uNode1Y, uNode2X, uNode2Y, ...)
    Vector<TYPE> displVec;
    elemDisp_.ConvertToVec_AppendCols(displVec);

    // linear differential operator B_lin
    Matrix<Double> linBMat;    
    calcBMat( linBMat, ip, ptCoord);

    Vector<TYPE> linStrain(linBMat.GetSizeRow());
    //Matrix<TYPE>(linBMat).Mult(displVec,linStrain);
    //Matrix<TYPE>(dMat).Mult(linStrain,stressVec);
    linStrain = linBMat * displVec;
    stressVec = dMat * linStrain;

  }

  /// calculates green-lagrangian strains (linear part, vector notation)
  // see Habil. M. Kaltenbacher 

  template <class TYPE>
  void MechStressStrain<TYPE>::
  CalcStrainVec(Vector<TYPE>& strainVec, UInt ip, Matrix<Double> & ptCoord)
  {
    ENTER_FCN( "MechStressStrain::CalcStrainVec" );

    Matrix<Double> dMat;
    linElastInt::calcDMat(dMat);
  
    // convert displacement of all elem nodes into one vector: 
    // (uNode1X, uNode1Y, uNode2X, uNode2Y, ...)
    Vector<TYPE> displVec;
    elemDisp_.ConvertToVec_AppendCols(displVec);

    // linear differential operator B_lin
    Matrix<Double> linBMat;    
    calcBMat( linBMat, ip, ptCoord);

    strainVec = linBMat * displVec;
    //strainVec.Resize(linBMat.GetSizeRow());
    //Matrix<TYPE>(linBMat).Mult(displVec,strainVec);
  
  }


  // returns linear B - matrix

  template <class TYPE>
  void MechStressStrain<TYPE>::
  calcBMat(Matrix<Double> & bMat, UInt ip, Matrix<Double> & ptCoord)
  {
    ENTER_FCN( "MechStressStrain::calcLinBMat" );

    // linear differential operator B_lin
    linElastInt::calcBMat(bMat, ip, ptCoord);
  }




#ifdef __GNUC__
  // Explicite template instantiation
  template class MechStressStrain<Double>;
  template class MechStressStrain<Complex>;
#endif

} // end namespace CoupledField
