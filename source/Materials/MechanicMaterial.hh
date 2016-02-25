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

    typedef StdVector<PtrCoefFct> CoefFctList; // list of coefFunctions
    typedef std::map<MaterialType, CoefFctList> CoefListMap; // map for lists of CoefFunctions like needed for relaxation tensors

    //! Default constructor
    MechanicMaterial(MathParser* mp,
                     CoordSystem * defaultCoosy);

    //! Destructor
    virtual ~MechanicMaterial(){};

    //! Trigger finalization of mataterial (calculation of rotated matrices)
    void Finalize();

    //! re-set the Re/Im part of a scalar material parameter
    void SetScalar( Double param, MaterialType matType, 
		    Global::ComplexPart dataType = Global::REAL);

    //! set a scalar complex material parameter
    void SetScalar( Complex param, MaterialType matType, 
		    Global::ComplexPart dataType = Global::REAL );
    
    //! set a scalar parameter in string representation (gets later parsed by MathParser
    void SetScalar(const std::string& param, MaterialType matType, Global::ComplexPart dataType = Global::REAL );

    //! set a real vector
    virtual void SetVector(const Vector<Double>& param, MaterialType matType,
			    Global::ComplexPart dataType );

    //! set a real material tensor
    void SetTensor(const Matrix<Double>& param, MaterialType matType,
		    Global::ComplexPart dataType );

    //! set a complex material tensor
    void SetTensor(const Matrix<Complex>& param, MaterialType matType,
		    Global::ComplexPart dataType );
    
    //! set the element of a coefFunctionList
    void SetCoefFctList(MaterialType matType, CoefFctList coefList);

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
    
    //! get an element from a real valued tensor list
    void GetTensor( Matrix<Double>& param, MaterialType matType, UInt index) const; 

    //! Return a specific sub-tensor
    virtual PtrCoefFct GetSubTensorCoefFnc( MaterialType matType, 
                                            SubTensorType tensorType,
                                            bool transposed );
    
    //! Return a list of coef functions
    CoefFctList GetCoefFctList( MaterialType matType);

    //! compute the complex valued stiffness tensor
    void RecomputeComplexViscoTensor();
    
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
    
    /** static helper function to calculate complex stiffness tensor from lame parameters */
    static void CalcIsotropicStiffnessTensor(Matrix<Complex>& out, Complex lambda, Complex mu);

    /** static helper function to calculate real stiffness tensor from lame parameters (real valued version)*/
    static void CalcIsotropicStiffnessTensor(Matrix<Double>& out, Double lambda, Double mu);
    
    /** static helper function to calculate Lame parameters from Young's Modulus and Poisson's ratio*/
    static void Enu2Lame(Double& LameLambda, Double& LameMu, Double E, Double nu);
    
    /** static helper function to calculate Lame parameters from Bulk and Shear modulus*/
    static void GK2Lame(Double& LameLambda, Double& LameMu, Double G, Double K);
    
    /** Computes from a given tensor the sub-type */
    template<class T>
    static void ComputeSubTensor(Matrix<T>& matMatrix, SubTensorType subTensor, const Matrix<T>& input);

    /** Calculates the orthotrope error.
     * Is only an informal heuristic. Sums up the values where there should be zero plus the balances */
    static double CalcOrthotropeError(const Matrix<double>& tensor);

    /** Calculates orthotrope material properties and gives them with a readable string.
     * @return first a description with underliner, then the value */
    static StdVector<std::pair<std::string, double> > CalcOrthotropeProperties(const Matrix<double>& tensor, BaseMaterial* mat, SubTensorType stt, double vol);
    
    /** */
    
    
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

    /** compute the correct subTensor (3D, AXI, ..) from the material */
    void ComputeSubTensor(Matrix<Complex>& matMatrix, MaterialType matType, SubTensorType subTensor) const;


    //! Compute elasticity tensor from given parameters
    void ComputeFullStiffTensor();

    //! Computae vector for thermal expansion
    void ComputeFullThermalExpanionVector();
    
    //! Connection to math parser instance 
    boost::signals2::connection conn_;
    
    //! Handle for expression
    MathParser::HandleType mHandle_;

  protected:
    //! map for lists of CoefFunctions like needed for relaxation tensors
    CoefListMap tensorCoefLists_;
  };

} // end of namespace

#endif
