/*
 * CoefFunctionRotation.cc
 *
 *  Created on: 02.12.2025
 *      Author: Michael Hofmarcher and Mark Socaciu Lendvai
 */

#include "CoefFunctionRotation.hh"
#include "Domain/ElemMapping/ElemShapeMap.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "MatVec/Vector.hh"
#include "MatVec/Matrix.hh"
#include "General/Exception.hh"

namespace CoupledField {
  DEFINE_LOG(coeffunctionRotation, "coeffunctionRotation")
  DEFINE_LOG(coeffunctionRotationValidation, "coeffunctionRotationValidation")

  CoefFunctionRotation::CoefFunctionRotation(PtrCoefFct matCoeff, PtrCoefFct vec1Coeff, PtrCoefFct vec2Coeff) 
  : CoefFunction(),
    matCoeff_(matCoeff),
    vec1Coeff_(vec1Coeff),
    vec2Coeff_(vec2Coeff)
  {
    dimType_ = TENSOR;
    isAnalytic_ = false;
    isComplex_  = false;
    dependType_ = SPACE;
  }
  
  void CoefFunctionRotation::GetTensor(Matrix<Double>& rotTensor, const LocPointMapped &lpm )
  {
    // Retrieve the vectors
    Vector<Double> v1, v2;
    vec1Coeff_->GetVector(v1, lpm);
    if (v1.NormL2() < 1e-16)
        throw Exception("CoefFunctionRotation: At least one provided reference vector has length zero!");
    vec2Coeff_->GetVector(v2, lpm);
        if (v2.NormL2() < 1e-16)
        throw Exception("CoefFunctionRotation: At least one provided target vector has length zero!");
    // Retrieve original tensor
    Matrix<Double> T;
    matCoeff_->GetTensor(T, lpm);
    LOG_DBG(coeffunctionRotation) << "material Tensor " << T;
    // normalized vectors
    double normv1 = v1.NormL2();
    double normv2 = v2.NormL2();
    // scalar product
    double vdot = v1 * v2;
    // calculation of the rotation angle
    double theta = acos( vdot / (normv1 * normv2) );
    LOG_DBG(coeffunctionRotation) << "theta: " << theta;
    // dimension of the rotation vectors
    int dim_v1 = v1.GetSize();
    int dim_v2 = v2.GetSize();
    // predefine rotation matrix
    Matrix<Double> R;

    // 2D rotation
    if ( T.GetNumRows() == 2 || T.GetNumCols() == 2 ) {
        if ( dim_v1 != 2 || dim_v2 != 2 )
            throw Exception("CoefFunctionRotation: 2D materialTensor needs 2D rotationVectorFields");
        // resize the rotation matrix to 2D
        R.Resize(2,2);
        // check for positive or negative rotation
        if( (v1[0]*v2[1] - v1[1]*v2[0]) > 0 ) { // positive rotation
            R(0,0) = cos(theta);   R(0,1) = - sin(theta);
            R(1,0) = sin(theta);   R(1,1) =  cos(theta);
        }
        else {  // negative rotation
            R(0,0) = cos(-theta);   R(0,1) = - sin(-theta);
            R(1,0) = sin(-theta);   R(1,1) =  cos(-theta);
        }
        LOG_DBG(coeffunctionRotation) << "Rotation Matrix 2D " << R;
    }
    else {  // 3D rotation
        if ( dim_v1 != 3 || dim_v2 != 3 )
            throw Exception("CoefFunctionRotation: 3D materialTensor needs 3D rotationVectorFields");
        // identity matrix
        Matrix<Double> I(3,3);
        I(0,0) = I(1,1) = I(2,2) = 1;
        LOG_DBG(coeffunctionRotation) << "Identity Matrix: " << I;
        // resize the rotation matrix to 3D
        R.Resize(3,3);
        // calculate the normalized rotation axis (nra)
        Vector<Double> ra;
        v1.CrossProduct(v2, ra);
        if ( v1.Collinear(v2)) {
            // v1 and v2 are collinear: cross product is zero
            // Distinguish parallel (R = I) from antiparallel (R = -I + 2*n*n^T)
            if (v1.Inner(v2) > 0) {     // parallel: no rotation needed
                R = I;
            }
            else {  // antiparallel: 180-degree rotation
                Vector<Double> v1n;
                v1n = v1 / v1.NormL2();
                // find index of smallest absolute component of v1n
                int minIdx = 0;
                double minVal = fabs(v1n[0]);
                if (abs(v1n[1]) < minVal) { minVal = fabs(v1n[1]); minIdx = 1; }
                if (abs(v1n[2]) < minVal) { minIdx = 2; }
                // create the unit vector of this direction to cross v1 with
                Vector<Double> helper(3);
                helper[minIdx] = 1.0;
                // new cross product to obtain the rotation axis
                v1n.CrossProduct(helper, ra);
                Vector<Double> nra;
                nra= nra / nra.NormL2();
                LOG_DBG(coeffunctionRotation) << "Antiparallel rotation axis: " << nra;

                Matrix<Double> nDyadic;
                nDyadic.DyadicMult(nra, nra);
                // ==== Rodrigues rotation R for antiparallel case====
                R = nDyadic * 2.0 + I * (-1.0);
                LOG_DBG(coeffunctionRotation) << "Rotation Matrix 3D antiparallel case " << R;
            }
        }
        else {
            Vector<Double> nra;
            nra = ra / ra.NormL2();
            LOG_DBG(coeffunctionRotation) << "normalized Rotation axis: " << nra;
            // Antisymmetric cross product matrix K 
            Matrix<Double> K(3,3);
            K(0,0)=0.0;      K(0,1)=-nra[2]; K(0,2)= nra[1];
            K(1,0)= nra[2];  K(1,1)=0.0;     K(1,2)=-nra[0];
            K(2,0)=-nra[1];  K(2,1)= nra[0]; K(2,2)=0.0;
            // dyadic product
            Matrix<Double> nraDyadic;
            nraDyadic.DyadicMult(nra, nra);
            LOG_DBG(coeffunctionRotation) << "Dyadic nra: " << nraDyadic;
            // ==== Rodrigues rotation R ====
            R = I * cos(theta) + K * sin(theta) + nraDyadic *(1 - cos(theta));
            LOG_DBG(coeffunctionRotation) << "Rotation Matrix 3D " << R;
        }
    }
    // rotate the material Tensor
    T.PerformRotation(R, rotTensor);
    LOG_DBG(coeffunctionRotation) << "Rotated material tensor " << rotTensor;

    LOG_DBG3(coeffunctionRotationValidation) 
    << "; " << (lpm.ptEl ? lpm.ptEl->elemNum : 0) << "; " 
    << v1.ToString() << "; " 
    << v2.ToString() << "; " 
    << T.ToString() << "; "
    << R.ToString() << "; "
    << rotTensor.ToString();

  }  
}
