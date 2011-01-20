// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef MECHANICMATERIAL_DATA
#define MECHANICMATERIAL_DATA

#include "baseMaterial.hh"

namespace CoupledField {

  //! Class for Material Data
  /*! 
     Class for handling mechanic material data
  */

  class MechanicMaterial : public BaseMaterial {

  public:

    //! Default constructor
    MechanicMaterial();

    //! Destructor
    ~MechanicMaterial();

    //! Trigger finalization of mataterial (calculation of rotated matrices)
    void Finalize();

   //! set a scalar string material parameter
    void SetScalar(const std::string& param, MaterialType matType);

    //! set a scalar real material parameter
    void SetScalar( Double param, MaterialType matType, 
		    Global::ComplexPart dataType );

    //! set a scalar complex material parameter
    void SetScalar( Complex param, MaterialType matType, 
		    Global::ComplexPart dataType );
    
    //! set a scalar parameter in string representation (gets later parsed by MathParser
    void SetScalar(const std::string& param, MaterialType matType, Global::ComplexPart dataType );

    //! set a real vector
    virtual void SetVector(const Vector<Double>& param, MaterialType matType,
			    Global::ComplexPart dataType );

    //! set a real material tensor
    void SetTensor(const Matrix<Double>& param, MaterialType matType,
		    Global::ComplexPart dataType );

    //! set a complex material tensor
    void SetTensor(const Matrix<Complex>& param, MaterialType matType,
		    Global::ComplexPart dataType );

    //! get a scalar string material parameter
    void GetScalar( std::string& param, 
                    MaterialType matType ) const;

     //! get a scalar integer material parameter
    void GetScalar( Integer& param, 
                    MaterialType matType) const;

    //! get a scalar real material parameter
    void GetScalar( Double& param, MaterialType matType, 
		    Global::ComplexPart dataType ) const;

    //! get a scalar complex real material parameter
    void GetScalar( Complex& param, MaterialType matType, 
		    Global::ComplexPart dataType ) const;

    //! get a real vector
    void GetVector( Vector<Double>& param, MaterialType matType,
			    Global::ComplexPart dataType ) const;

    //! get a real material tensor
    void GetTensor( Matrix<Double>& param, MaterialType matType,
		    Global::ComplexPart dataType,
		    SubTensorType = FULL ) const;	

    //! get a complex material tensor
    void GetTensor( Matrix<Complex>& param, MaterialType matType,
		    Global::ComplexPart dataType,
		    SubTensorType = FULL ) const;	
    
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

  private:

    /** Calculates orthotrope Youngs moduli. Not for plane_stress!!
     * @return a vector with 3 entries E_1, E_2, E_3 where E_3 = vol * E_core in 2D */
    static StdVector<double> CalcOrthotropeYoungsModulus(const Matrix<double>& tensor, BaseMaterial* mat, SubTensorType stt, double vol);

    /** Calculates orthotrope Youngs moduli
     * in 2D the extensions to 3D are from the plane strain extensions. Works not for plane stress yet
     * @return a vector with 6 entries  v_21, v_12, v_31, v_13, v_32, v_23 */
    static StdVector<double> CalcOrthotropePoissonsRatio(const Matrix<double>& tensor, BaseMaterial* mat, SubTensorType stt, double vol);


    /** compute the correct subTensor (3D, AXI, ..) from the material */
    void ComputeSubTensor(Matrix<Complex>& matMatrix, MaterialType matType, SubTensorType subTensor) const;


    //! Compute elasticity tensor from given parameters
    void ComputeFullStiffTensor();

    
    MathParser::HandleType mHandle_;
    
    Double density_;
    Double PoissonRatio_;
    Double RayleighAlpha_;
    Double RayleighBeta_;
    Double RayleighFrequency_;
    Double lossTangens_;

    Complex scalarEmodulus_;
    Complex scalarLameLambda_;
    Complex scalarLameMu_;

    Matrix<Complex> stiffnessTensor_;

  };

} // end of namespace

#endif
