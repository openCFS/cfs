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
     * Considers only the right upper diagonal matrix
     * 2D: E11-E22 = 0, E11-E12-2E33 = 0, E12 = 0, E23 = 0
     * @param normed divides every entry by E11
     * @return always positive as we add only abs(). */
    static double CalcIsotropyError(const Matrix<double>& tensor, bool normed);

    /** Compute the Young's modulus out of an isotropic tensor.
     * You might want to check via ComputeIsotropyError() first.
     * Does not check for isotropy but is brute force! */
    static double CalcIsotropicYoungsModulus(const Matrix<double>& tensor);

    /** Compute the Poisson's ratio out of an isotropic tensor.
     * @see ComputeIsotropicYoungsModulus() */
    static double CalcIsotropicPoissonsRatio(const Matrix<double>& tensor);

    /** static helper function to calculate real stiffness tensor from elasticity modulus and poisson ratio
      * in 2D we compute the plain strain state */ 
    static void CalcIsotropicStiffnessTensorFromEAndPoisson(Matrix<Double>& out, Double emod, Double poisson, UInt dim);
    
    /** static helper function to calculate real stiffness tensor from lame parameters
     *  in 2D we compute the plain strain state */
    static void CalcIsotropicStiffnessTensorFromLame(Matrix<Double>& out, Double lambda, Double mu, UInt dim);
    
    /** static helper function to calculate complex stiffness tensor from lame parameters
      *  in 2D we compute the plain strain state */
    static void CalcComplexIsotropicStiffnessTensor(Matrix<Complex>& out,
        const Complex &lambda, const Complex &mu, UInt dim);

    /** Calculates orthotrope Youngs moduli
     * @return a vector with 2 or 3 entries E_1, E_2 (, E_3) */
    static StdVector<double> CalcOrthotropeYoungsModulus(const Matrix<double>& tensor);

    /** Calculates orthotrope Youngs moduli
     * @return a vector with 2 or 6 entries v_21, v_12 or v_32, v_31, v_21, v_23, v_13, v_12 */
    static StdVector<double> CalcOrthotropePoissonsRatio(const Matrix<double>& tensor);

    /** Calculates the orthotrope error.
     * Is only an informal heuristic. Sums up the values where there should be zero plus the balances */
    static double CalcOrthotropeError(const Matrix<double>& tensor);


    /** Calculates orthotrope material properties and gives them with a readable string.
     * @return first a description with underliner, then the value */
    static StdVector<std::pair<std::string, double> > CalcOrthotropeProperties(const Matrix<double>& tensor);

  private:


    //! compute the2 correct subTensor (3D, AXI, ..)
    void ComputeSubTensor(Matrix<Complex>& matMatrix,
			  MaterialType matType, 
			  SubTensorType subTensor) const;

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
