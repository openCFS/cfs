#ifndef DESIGNMATERIAL_HH_
#define DESIGNMATERIAL_HH_

#include <map>
#include <complex>
#include <cmath>

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "General/Enum.hh"
#include "General/Environment.hh"
#include "MatVec/Matrix.hh"
#include "Optimization/Design/DesignElement.hh"
#include "Optimization/Design/MaterialTensor.hh"
#include "Optimization/OptimizationMaterial.hh"
#include "Utils/ApproxData.hh"

#include "def_use_sgpp.hh"

#ifdef USE_SGPP
#include "sgpp/base/grid/Grid.hpp"
#include "sgpp/base/operation/hash/OperationEval.hpp"
#include "sgpp/base/operation/hash/OperationNaiveEval.hpp"
#include "sgpp/base/operation/hash/OperationNaiveEvalPartialDerivative.hpp"
#include "sgpp/base/datatypes/DataVector.hpp"
#endif

namespace CoupledField {

/** This implements a function from $R^n$ to $R^{d \times d}$ for transforming a vector of Parameters
 * to a material tensor.  */
//template <class TYPE> class StdVector;

class DesignSpace;
class CoefFunctionOpt;

class DesignMaterial
{
public:
  typedef enum { NO_TYPE=-1, FMO, ISOTROPIC, LAME_ISOTROPIC, TRANSVERSAL_ISOTROPIC, TRANSVERSAL_ISOTROPIC_BOXED,
    DENSITY_TIMES_TRANSVERSAL_ISOTROPIC, DENSITY_TIMES_TRANSVERSAL_ISOTROPIC_BOXED, DENSITY_TIMES_ROT_TRANSVERSAL_ISOTROPIC,
    DENSITY_TIMES_ROT_TRANSVERSAL_ISOTROPIC_BOXED, DENSITY_TIMES_ROT_PA12, ORTHOTROPIC,
    DENSITY_TIMES_ORTHOTROPIC, DENSITY_TIMES_2D_TENSOR, DENSITY_TIMES_2D_TENSOR_CONSTANT_TRACE, DENSITY_TIMES_ROTATED_2D_TENSOR,
    D_INTERP_IN718_TENSOR, D_INTERP_IN718_TENSOR_ROT, LAMINATES, D_LAMINATES,
    HOM_RECT, D_HOM_RECT, HOM_RECT_C1, HOM_ISO_C1, MSFEM_C1, SGP_MATLAB, SGP_GRADIENTCHECK, HEAT, FEATURE_MAPPING_ANISO } Type;

    /* posibilities for the isotropic plane in transversal isotropy
     * note that parameters EMODULISO, POISSONISO are used for that plane
     * EMODUL is in the orthogonal direction, POSSION is nu_io where i is in the isotropic plane, o not
     * GMODUL is G_io where i is in the isotropic plane o not (note G_io = G_jo) */
    typedef enum { TRANSISO_XY, TRANSISO_YZ, TRANSISO_XZ } TransIsoType;

//    /** Material notation. Only for FMO we assume the design to be Hill-Mandel, in LinElastInt we use Voigt. The CFS-B-operator is also Voigt */
//    typedef enum {  NO_TYPE=-1, VOIGT, HILL_MANDEL } Notation;

    /** Method used for interpolation of material tensor and volume */
    typedef enum { NOTYPE=-1, C1, SG, FULL_BSPLINE } Interpolation;

    /** Types of rotation, the strings are read from the xml file
     * The order of application is read from right to left, hence the default xyz means that
      we rotate first around z, then y and finally araound x.
      We apply conventions of extrinsic rotation, see
      https://en.wikipedia.org/wiki/Euler_angles#Conventions_by_extrinsic_rotations
      The formats ABA are Euler angles, the format ABC are Tait–Bryan angles (Kardanwinkel in German) */
    typedef enum { ZXZ, ZYZ, YZY, YXY, XYX, XZX, XYZ, YXZ, XZY, ZXY, ZYX, YZX } RotationType;

    /** constructor, reads in DesignMaterial from XML
     * @param pn pointer to PtrParamNode */ 
    DesignMaterial(PtrParamNode pn, OptimizationMaterial::System material, StdVector<DesignID>& design, DesignSpace* space);

    /** constructor for fake DesignMaterial to be used by FeatureMappingParamMat */
    DesignMaterial(Type type, DesignSpace* space);

    ~DesignMaterial();

    bool GetTensor(MaterialTensor<double>& mt, DesignElement::Type type, SubTensorType subTensor, const Elem* elem, DesignElement::Type direction, MaterialTensorNotation notation, bool core = false);

    /** Calculate the derivative tensor from the given material parameters
     * Sets the Tensor in VOIGT notation.
     * @param core if true, return the material tensor without multiplication with penalized pseudo-density
     * @param coef only for FEATURE_MAPPING_ANISO */
    bool GetMechTensor(MaterialTensor<double>& mt, SubTensorType subTensor, const Elem* elem, DesignElement::Type direction = DesignElement::NO_DERIVATIVE, bool core = false, const CoefFunctionOpt* coef = nullptr);
    bool GetMechTensor(MaterialTensor<Complex>& mt, SubTensorType subTensor, const Elem* elem, DesignElement::Type direction = DesignElement::NO_DERIVATIVE, bool core = false, const CoefFunctionOpt* coef = nullptr);

    bool GetPiezoCouplingTensor(MaterialTensor<double>& mt, const Elem* elem, DesignElement::Type direction);

    /** Calculates MSFEM element matrix for a regular grid from material catalogue*/
    bool GetErsatzElementMatrixMSFEM(Matrix<double>& mt,const Elem* elem, DesignElement::Type direction);

    /** returns the tensor with negative design variables such the design vector is still pos. definite */
    bool GetElecTensor(MaterialTensor<double>& mt,  const Elem* elem, DesignElement::Type direction);

    /** retrieve rel. mass of element (tensor trace) or direct design variable */
    double GetMechMass(const Elem* elem, DesignElement::Type direction);

    /** return whether mass is also a design (else it is calculated from tensor) */
    bool MassIsDesign(){
      return massIsDesign_;
    }

    /** return whether damping is also design */
    bool DampingIsDesign(){
      return dampingIsDesign_;
    }

    /** retrieve damping parameters for element or derivative */
    bool GetMaterialDamping(double& alpha, double& beta, DesignElement::Type direction = DesignElement::NO_DERIVATIVE);

    /** Get a parameter of the parametric material optimization. Takes it from the thread local storage. */
    double GetParameter(const DesignElement::Type p) {
      assert(HasParameter(p));
      return params_.Mine()[p];
    }

    void static SetEnums();

    Type GetType() const { return type_; };

    void SetType(Type type) {type_ = type;};

    /** shall we add the coef->orgMat when we do not do derivative? */
    bool HasBias() const { return bias_; }

    double CalcHomVolume(Vector<double>& p, DesignElement::Type direction, bool derivative);

    Interpolation GetInterpolationMethod() const { return interpolation_; };

    /** rotate elasticity tensor in Voigt notation according to the parameters, eventually calculating a derivative.
     * The order of rations is determined by rotationType_
     *  in 3d: rotates the material by ROTANGLEFIRST around the first axis, by ROTANGLESECOND around the second axis and by ROTANGLETHIRD around the third axis in this given order or rz,ry,rx
     *  in 2d: rotates the material by ROTANGLE or rx
     *  rotation is always mathematically positive, i.e. counterclockwise
     * @param mt Material Tensor which is rotated in place (or the derivative is calculated in place)
     * @param direction if one of ROTANGLEFIRST, ROTANGLESECOND, ROTANGLETHIRD, ROTANGLE calculate the derivative of the rotation w.r.t. this parameter
     * @param angles is true if rotation angles rx,ry,rz are given by parameter, otherwise false
     * @see rotationType_
     */
    void RotateTensor(MaterialTensor<double>& mt, DesignElement::Type direction, bool angles = false, double rx = 0., double ry = 0., double rz = 0.);

    /** Calculate the Isotropic tensor */
    void GetIsoMaterialTensor(MaterialTensor<double>& mt, SubTensorType subTensor, DesignElement::Type direction);

    /** little helper for GetInterpolatedHomTensor().
     * @param vector p has the values of the design variable */
    void GetHomC1Tensor(MaterialTensor<double>& mt, Vector<double>& p, DesignElement::Type direction, SubTensorType subTensor);

    /** little helper for GetInterpolatedHomTensor().
     * @param vector p has the values of the design variable */
    void GetHomIsoC1Tensor(MaterialTensor<double>& mt, Vector<double>& p, DesignElement::Type direction, SubTensorType subTensor) const;

    void ToInfo(PtrParamNode in) const;

    /** evaluates the one dimensional C1 interpolation polynomial at point p and returns function value as double
     * f(x) = a0 + a1*p + a2*p**2 + a3*p**3 */
    double EvaluateC1Interpolation1D(double p, const Matrix<double>& coeff, int& j) const;

    /** evaluates the first derivative of the one dimensional C1 interpolation polynomial at point p
     * f'(x) = a1 + 2*a2*p + 3*a3*p**2
     * direction is only passed for an assertion */
    double EvaluateC1Interpolation1D_Deriv(double p, const Matrix<double> & coeff, int & j, DesignElement::Type direction) const;

    /** evaluates the second derivative of the one dimensional C1 interpolation polynomial at point p
     * f'(x) = 2*a2 + 6*a3*p
     * direction is only passed for an assertion */
    double EvaluateC1Interpolation1D_Deriv2(double p, const Matrix<double> & coeff, int & j, DesignElement::Type direction) const;

    /** evaluates the twodimensional C1 interpolation polynomial at point p and returns function value as double */
    double EvaluateC1Interpolation2D(Vector<double>& p, const Matrix<double>& coeff, int& j, int& k, int& m, int& n) const;

    /** evaluates the derivative of the C1 interpolation polynomial at point p w.r.t. direction */
    double EvaluateC1Interpolation2D_Deriv(Vector<double>& p, const Matrix<double>& coeff, int& j, int& k, int& m, int& n, DesignElement::Type direction) const;

    /** evaluates the threedimensional C1 interpolation polynomial at point p and returns function value as double */
    double EvaluateC1Interpolation3D(Vector<double>& p, const Matrix<double>& coeff, int& j, int& k,int& l, int& m, int& n, int& o) const;

    /** evaluates the derivative of the C1 interpolation polynomial at point p w.r.t. direction */
    double EvaluateC1Interpolation3D_Deriv(Vector<double>& p, const Matrix<double>& coeff, int& j, int& k, int& l, int& m, int& n, int& o, DesignElement::Type direction) const;

    /** evaluates the derivative of an interpolator at point p w.r.t. direction*/
    double EvaluateC1Interpolation_Deriv(Vector<double>& p, ApproxData* interpolator, DesignElement::Type direction) const;

    /** Get the index of the local interpolation interval
     * if point is outside of interval, it is set to interval's bounds*/
    int GetInterpolationIndex(const Vector<double>& interval, double& point) const;

    /** Calculates a local load factor for normalized stress by evaluating
     * an interpolated function (stored in the material catalogue file).
     * Also does a quadratic extrapolation to avoid too large values for vol -> 1. */
    double GetMicroLoadFactor(double vol, DesignElement::Type direction);

    /** checks for a parameter. Checks the thread local storage. */
    inline bool HasParameter(const DesignElement::Type p) const {
      return params_.ConstMine().find(p) != params_.ConstMine().end();
    };

protected:

    /** Set a parameter for the parametric material optimization.
     * @param global set the value only thread local (element values) or global (as in constructor) on all thread versions */
    void SetParameter(const DesignElement::Type p, const double value, bool global);

    /** to access the the local copy of the map but have the checks in debug mode */
    double GetParameter(const std::map<DesignElement::Type, double>& map, const DesignElement::Type p);

    /** Sets all Material Parameters in designMaterial for given element in the current context */
    bool CollectMaterialParametersForElement(DesignSpace* space, const Elem* elem);

    /** Return a local thread local data set to ease access within a function. Don't store!! */
    const std::map<DesignElement::Type, double>& GetParameters() const;

    /** for debugging */
    void DumpParams();

    /** this are the design variables for the current element!
     * It is filled with default and optimization variables.
     * Set by CollectMaterialParametersForElement().
     * To be thread save this needs to be in a thread local storage with CFS_NUM_THREADS.
     * Better don't access the map manually!! */
    CfsTLS<std::map<DesignElement::Type, double> > params_;

    /** mass is considered an independent design */
    bool massIsDesign_ = false;

    /** damping is also optimized */
    bool dampingIsDesign_ = false;

    /** shearing is optimized */
    bool shearIsDesign_ = false;

    /** dimension of material catalogue */
    StdVector<double> catalogueSize_;

    /** multiply mass with this, can be used to scale tensor trace */
    double massFactor_ = -1.0;

    /** for density times 2d tensor with constant trace */
    double trace_ = -1.0;

    static Enum<Type> type;
    Type type_ = NO_TYPE;

    static Enum<TransIsoType> transIsoType;
    TransIsoType transIsoType_ = TRANSISO_YZ;

    unsigned int dim = 0;

#ifdef USE_SGPP
    /** Grid for SGPP interpolation */
    std::unique_ptr<sgpp::base::Grid> grid_;
#endif

    /** returns the numbers of parameters required for this material */
    unsigned int RequiredParameters(OptimizationMaterial::System material);

    /** Check whether all required designs are available */
    bool CheckRequiredDesigns(StdVector<DesignElement::Type>& design);
private:
    /* note that most of these functions are called really often, so inlining is used */

    bool ReadCoeff(PtrParamNode pn, const string& name, int nRows, int nCols, Matrix<double>& coeff) const;

    /** create a new interpolator (cubic|bicubic|tricubic) with coeff as polynomial coefficients */
    ApproxData* CreateInterpolator(StdVector<double>& a, StdVector<double>& b, StdVector<double>& c, Matrix<double>& coeff);

    /** Calculate the Lame Tensor */
    void GetLameMaterialTensor(MaterialTensor<double>& mt, SubTensorType subTensor, DesignElement::Type direction);

    /** Calculate the Trans-Iso Tensor
     * @param core if true, return the material tensor without multiplication with penalized pseudo-density */
    void GetTransIsoMaterialTensor(MaterialTensor<double>& mt, SubTensorType subTensor, DesignElement::Type direction, bool core = false);

    /* general anisotropic FMO tensor */
    void GetElasticFMOTensor(MaterialTensor<double>& mt, SubTensorType subTensor, DesignElement::Type direction);

    /** Calculate tensor for free orthotropic material formulation
     *                    |e11 e12 0 |^2
     * E(e11,e22,e33,e12)=|e12 e22 0 |  + diag(lowerEigenvalueBound)
     *                    | 0   0 e33|
     *  */
    void GetOrthotropicMaterialTensor(MaterialTensor<double>& mt, SubTensorType subTensor, DesignElement::Type direction);

    /** Calculate the Tensor for Density times Tensor
     * @param core if true, return the material tensor without multiplication with penalized pseudo-density */
    void GetDensityTimes2dTensorTensor(MaterialTensor<double>& mt, SubTensorType subTensor, DesignElement::Type direction, bool core = false);

    /** Calculate the tensor for Laminates */
    void GetLaminatesTensor(MaterialTensor<double>& mt, SubTensorType subTensor, DesignElement::Type direction);

    /** little helper for GetInterpolatedHomTensor(). We assume we are in Hill-Mandel world
     * @param shape might also be the x or y component of the derivative! */
    void GetHomRectTensor(MaterialTensor<double>& mt, const Vector<double>& shape) const;

    /** Approximates the homogenized tensor of an a-b rectangle as used by Bendsoe and Kikuchi 1988 */
    void GetInterpolatedHomTensor(MaterialTensor<double>& mt, SubTensorType subTensor,  const Elem* elem,  DesignElement::Type direction);

    /** Gives the elasticity tensor of the Nickel basis alloy IN718 with interpolation between different crystalline microstructures */
    void GetIN718Tensor(MaterialTensor<double>& mt, SubTensorType subTensor,  DesignElement::Type direction);

    void ZeroMatrix(Matrix<double>& t, SubTensorType subTensor);

    /** set values in 2 matrix "t" in Voigt notation order */
    void Set2dMatrix(Matrix<double>& t, double t11, double t22, double t33, double t23, double t13, double t12);

    /** set values in 2 matrix "t" in Voigt notation order (doesn't assume symmetry) */
    void Set2dMatrix(Matrix<double>& t, double t11, double t12, double t13, double t21, double t22, double t23, double t31, double t32, double t33);

    /** put the entries of the 3D matrix at the right places (Voigt notation order - whatever that is in 3d???) */
    void Set3dMatrix(Matrix<double>& t, SubTensorType subTensor, double,double,double,double,double,double,double,double,double,double,double,double,double,double,double,double,double,double,double,double,double);

    /** put the entries of the transversal_isotropic matrix at the "right" places */
    void SetTransIsoMatrix(Matrix<double>& t, SubTensorType subTensor, double iD, double inD, double iG, double oD, double onD, double oG);

    /** put the entries of the isotropic matrix at the "right" places */
    void SetIsoMatrix(Matrix<double>& t, SubTensorType subTensor, double D, double nD, double G);

    // rotation matrix in 2d around z axis or in 3d around chosen coordinate axis
    void SetOneAxisRotationMatrix(Matrix<double>& R, double theta, int axis, bool derivative = false);

    /** helper function to set a rotation matrix of size 3x3
     * rotation axes ares given by rotationType_ (default XYZ, i.e. first z, then y, then x)
     * @param R the place to set the rotation matrix
     * @theta1 angle for rotation around first axis
     * @theta2 angle for rotation around second axis
     * @theta3 angle for rotation around third axis
     * @direction if given direction of the derivative to be calculated */
    void SetRotationMatrix(Matrix<double>& R, double theta1, double theta2, double theta3, DesignElement::Type direction = DesignElement::NO_DERIVATIVE);

    /** This exists only in Voigt notation! */
    void RotatePiezoCouplingTensor(Matrix<double>& t, double angle, DesignElement::Type direction);

    void RotateElecTensor(MaterialTensor<double>& mt, double angle, DesignElement::Type direction);

    /** Calculate the mass isotropic case */
    double GetIsoMaterialMass(DesignElement::Type direction);

    /** Calculate the mass lame case */
    double GetLameMaterialMass(DesignElement::Type direction);

    /** Calculate the mass trans-iso case */
    double GetTransIsoMaterialMass(DesignElement::Type direction);

    /** Calculate the mass density-times-tensor case
     * This returns the scaling factor (pseudo-density) for the normal mass matrix based on the materials actual density */
    double GetDensityTimesTensorMass(DesignElement::Type direction);

    /** Get the trans-iso mass (tensor trace) out of the corresponding tensor entries */
    double GetTransIsoMass(double iD, double iG, double oD, double oG);

    /** Get the isotropic mass (tensor trace) out of the corresponding tensor entries */
    double GetIsoMass(double D, double G);

    /** fills the row in hom_rect_samples_ */
    void FillHomRectSamples(PtrParamNode homRect, unsigned int idx, const std::string& a, const std::string& b);


    /** fills the coefficient data structure for the bicubic interpolation*/
    void FillHomRectCoeff(Matrix<double> & coeff_,const char * filename);

    /** Read detailed stats from file*/
    bool ReadDetailedStats(const char * filename, Matrix<double>& ret);

#ifdef USE_SGPP
    /** little helper for GetInterpolatedHomTensor(). We assume we are in Hill-Mandel world
     * @param vector p has the values of the design variable */
    void GetHomRectSGPPTensor(Matrix<double>& E, Vector<double>& p, DesignElement::Type direction, SubTensorType subTensor);

    /** little helper for GetInterpolatedHomTensor(). We assume we are in Hill-Mandel world
     * @param vector p has the values of the design variable */
    void GetHomRectFullBsplineTensor(MaterialTensor<double>& mt, Vector<double>& p, DesignElement::Type direction, SubTensorType subTensor) const;

    void EvaluateFullGrid();

    /** Fill sparse grid with data values*/
    void FillSparseGridWithFullGridData(Matrix<double>& data);
    void FillSparseGridWithSparseGridData(Matrix<double>& data);
    void HierarchizeSparseGridCoefficients();

    /** Initialize sparse grid for interpolation*/
    void InitializeSparseGrid(const char * filename);

    /** evaluates the derivative of the sgpp interpolation at point point in direction direction*/
    double EvaluateSGPPInterpolation_Deriv(sgpp::base::DataVector& alpha, sgpp::base::DataVector& point, DesignElement::Type direction) const;
#endif //USE_SGPP

    /** what convention for order of rotations we use (xyz, ...) */
    static Enum<RotationType> rotationType;
    RotationType rotationType_ = RotationType::XYZ; // first around z, then y, finally x

    /** sampled values for a single hom-rect 9-element by the number of shape function. Notation is Hill-Mandel!
     * 9 rows and 6 columns for with TENSOR11 being the first */
    Matrix<double> hom_rect_samples_;
    /** sampled values for coefficients of the interpolation polynomial; number of sample elements is rows and 16 columns/64 columns (3D)*/
    Matrix<double> hom_rect_coeff11_;
    Matrix<double> hom_rect_coeff12_;
    Matrix<double> hom_rect_coeff13_;
    Matrix<double> hom_rect_coeff14_;
    Matrix<double> hom_rect_coeff15_;
    Matrix<double> hom_rect_coeff16_;
    Matrix<double> hom_rect_coeff22_;
    Matrix<double> hom_rect_coeff23_;
    Matrix<double> hom_rect_coeff24_;
    Matrix<double> hom_rect_coeff25_;
    Matrix<double> hom_rect_coeff26_;
    Matrix<double> hom_rect_coeff33_;
    Matrix<double> hom_rect_coeff34_;
    Matrix<double> hom_rect_coeff35_;
    Matrix<double> hom_rect_coeff36_;
    Matrix<double> hom_rect_coeff44_;
    Matrix<double> hom_rect_coeff45_;
    Matrix<double> hom_rect_coeff46_;
    Matrix<double> hom_rect_coeff55_;
    Matrix<double> hom_rect_coeff56_;
    Matrix<double> hom_rect_coeff66_;
    Vector<double> hom_rect_a_;
    Vector<double> hom_rect_b_;
    Vector<double> hom_rect_c_;

    /** sampled values for coefficients of the interpolation polynomial for a micro load factor */
    Matrix<double> hom_mlf_coeff_;
    /** micro load factors are extrapolated above a given threshold */
    double extrapolationThreshold_ = 0.0;

    /** MSFEM element matrix coefficients of the bi-/tricubic interpolation polynomial from material catalogue; number of sample elements rows and 64 columns */
    Vector<double> msfem_a_;
    Vector<double> msfem_b_;
    Vector<double> msfem_c_;
    StdVector<Matrix<double> > msfem_coeff_;

    DesignSpace* space_ = NULL;

    /** bias means, that we add the cfs region material for any non-derivative.
     * The bias material is simply CoefFunctionOpt::orgMat */
    bool bias_ = false;

    Interpolation interpolation_ = NOTYPE;
    unsigned int level_ = 0;

    ApproxData* interpolator11_ = nullptr;
    ApproxData* interpolator12_ = nullptr;
    ApproxData* interpolator13_ = nullptr;
    ApproxData* interpolator14_ = nullptr;
    ApproxData* interpolator15_ = nullptr;
    ApproxData* interpolator16_ = nullptr;
    ApproxData* interpolator22_ = nullptr;
    ApproxData* interpolator23_ = nullptr;
    ApproxData* interpolator24_ = nullptr;
    ApproxData* interpolator25_ = nullptr;
    ApproxData* interpolator26_ = nullptr;
    ApproxData* interpolator33_ = nullptr;
    ApproxData* interpolator34_ = nullptr;
    ApproxData* interpolator35_ = nullptr;
    ApproxData* interpolator36_ = nullptr;
    ApproxData* interpolator44_ = nullptr;
    ApproxData* interpolator45_ = nullptr;
    ApproxData* interpolator46_ = nullptr;
    ApproxData* interpolator55_ = nullptr;
    ApproxData* interpolator56_ = nullptr;
    ApproxData* interpolator66_ = nullptr;
    ApproxData* interpolatorMLF_ = nullptr;

#ifdef USE_SGPP
    /** members for SGPP interpolation */
    enum SGPPBasis { LINEAR, MODLINEAR, BSPLINE, MODBSPLINE } sgpp_basis_;
    unsigned int bspline_degree_;
    sgpp::base::DataVector alpha1_;
    sgpp::base::DataVector alpha2_;
    sgpp::base::DataVector alpha3_;
    sgpp::base::DataVector alpha4_;
    sgpp::base::DataVector alpha5_;
    sgpp::base::DataVector alpha6_;
    sgpp::base::DataVector alpha7_;
    Matrix<double> full_bspline_coeff11_;
    Matrix<double> full_bspline_coeff12_;
    Matrix<double> full_bspline_coeff13_;
    Matrix<double> full_bspline_coeff22_;
    Matrix<double> full_bspline_coeff23_;
    Matrix<double> full_bspline_coeff33_;
    std::unique_ptr<sgpp::base::OperationEval> op_eval_;
    std::unique_ptr<sgpp::base::OperationNaiveEval> op_naive_eval_;
    std::unique_ptr<sgpp::base::OperationNaiveEvalPartialDerivative> op_naive_eval_partial_derivative_;
#endif //USE_SGPP
};

} // namespace

#endif /*DESIGNMATERIAL_HH_*/
