// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef MECHANICMATERIAL_DATA
#define MECHANICMATERIAL_DATA

#include "BaseMaterial.hh"

namespace CoupledField {

  //! Class for Material Data
  /*! 
     Class for handling mechanic material data
  */

  class MechanicMaterial : public BaseMaterial {

  public:

    //! Default constructor
    MechanicMaterial(MathParser* mp,
                     CoordSystem * defaultCoosy);

    //! Destructor
    virtual ~MechanicMaterial(){};

    //! Trigger finalization of mataterial (calculation of rotated matrices)
    void Finalize();

    //! Return scalar-valued coefficient function (linear)
    virtual PtrCoefFct GetScalCoefFnc(MaterialType matType,
                                      Global::ComplexPart matDataType) const;

    //! get sub vector (e.g. thermal expansion tensor in Voigt notation)
    virtual PtrCoefFct GetSubVectorCoefFnc( MaterialType matType,
                                            SubTensorType tensorType,
                                            Global::ComplexPart matDataType
                                                = Global::COMPLEX) const;

    //! Return a specific sub-tensor
    virtual PtrCoefFct GetSubTensorCoefFnc( MaterialType matType, 
                                            SubTensorType tensorType,
                                            Global::ComplexPart matDataType,
                                            bool transposed = false) const;
    
    /** Computes the error to an isotropic elasticity tensor.
     * Assume isotropy and calculate E and v, construct E(E,v) and return ||E(E,v) - tensor||_1 */
    static double CalcIsotropyError(const Matrix<double>& tensor, SubTensorType stt);

    /** Compute the Young's modulus out of an isotropic tensor.
     * You might want to check via ComputeIsotropyError() first.
     * Does not check for isotropy but is brute force! */
    static double CalcIsotropicYoungsModulus(const Matrix<double>& tensor, SubTensorType subTensor);

    /** Compute the Poisson's ratio out of an isotropic tensor.
     * @see ComputeIsotropicYoungsModulus() */
    static double CalcIsotropicPoissonsRatio(const Matrix<double>& tensor, SubTensorType subTensor);

    /** @see CalcOrthotropeProperties() */
    static StdVector<std::pair<std::string, double> > CalcIsotropicProperties(const Matrix<double>& tensor, SubTensorType stt);

    /** static helper function to calculate real stiffness tensor from elasticity modulus and Poisson's ratio (3D) */
    static void CalcIsotropicStiffnessTensorFromEAndPoisson(Matrix<Double>& out, Double emod, Double poisson);
    
    /** static helper function to calculate real stiffness tensor from lame parameters (3D) */
    static void CalcIsotropicStiffnessTensorFromLame(Matrix<Double>& out, Double lambda, Double mu);
    
    /** static helper function to calculate complex stiffness tensor from lame parameters */
    static void CalcComplexIsotropicStiffnessTensor(Matrix<Complex>& out, Complex lambda, Complex mu);

    /** Computes from a given tensor the sub-type */
    template<class T>
    static void ComputeSubTensor(Matrix<T>& matMatrix, SubTensorType subTensor, const Matrix<T>& input);

    /** Calculates the orthotrope error.
     * Is only an informal heuristic. Sums up the values where there should be zero plus the balances */
    static double CalcOrthotropeError(const Matrix<double>& tensor);

    /** Calculates orthotrope material properties and gives them with a readable string.
     * @return first a description with underliner, then the value */
    static StdVector<std::pair<std::string, double> > CalcOrthotropeProperties(const Matrix<double>& tensor, BaseMaterial* mat, SubTensorType stt, double vol);

    //! Compute the constant full stiffness tensor from isotropic,
    //! transversely isotropic or orthotropic data
    static Matrix<Complex> GetFullStiffTensor(BaseMaterial::SymmetryType symType,
                                              BaseMaterial::CoefMap &coefMap);

    //! Compute the constant full viscous tensor from the Kelvin-Voigt definition
    static Matrix<Complex> GetFullViscousTensorKV(BaseMaterial::SymmetryType symType, BaseMaterial::CoefMap &coefMap);

  private:

    /** Calculates orthotrope Youngs moduli.
     * @see CalcOrthotropePoissonsRatio() for the plane strain/ plane stress stuff */
    static StdVector<double> CalcOrthotropeYoungsModulus(const Matrix<double>& tensor, BaseMaterial* mat, SubTensorType stt, double vol);

    /** Calculates orthotrope Youngs moduli
     * 3D: 6 results
     * 2D: plane strain: 6 results with E3 assumed as core_E*vol
     * 2D: plane stress: v_21 and v_12
     * @param mat only for plane strain
     * @param vol only for plane strain
     * @return a vector with 2 or 6 entries  v_21, v_12(, v_31, v_13, v_32, v_23) */
    static StdVector<double> CalcOrthotropePoissonsRatio(const Matrix<double>& tensor, BaseMaterial* mat, SubTensorType stt, double vol);


    //! Compute elasticity tensor from given parameters
    void ComputeFullStiffTensor();

    //! Compute the viscoelastic stiffness tensor from the Prony series
    void ComputeViscoStiffTensors();

    //! Compute the thermal expansion
    void ComputeThermExpTensor();
  };

} // end of namespace

#endif
