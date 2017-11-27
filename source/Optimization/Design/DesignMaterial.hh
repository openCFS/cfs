#ifndef DESIGNMATERIAL_HH_
#define DESIGNMATERIAL_HH_

#include <map>
#include <complex>
#include <cmath>

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "Optimization/Design/DesignElement.hh"
#include "General/Enum.hh"
#include "General/Environment.hh"
#include "Optimization/OptimizationMaterial.hh"
#include "MatVec/Matrix.hh"
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
template <class TYPE> class StdVector;

class ErsatzMaterial;
class DesignSpace;
class TransferFunction;

  class DesignMaterial
  {
  public:
    typedef enum { FMO, ISOTROPIC, LAME_ISOTROPIC, TRANSVERSAL_ISOTROPIC, TRANSVERSAL_ISOTROPIC_BOXED,
      DENSITY_TIMES_TRANSVERSAL_ISOTROPIC, DENSITY_TIMES_TRANSVERSAL_ISOTROPIC_BOXED, DENSITY_TIMES_ROT_TRANSVERSAL_ISOTROPIC,
      DENSITY_TIMES_ROT_TRANSVERSAL_ISOTROPIC_BOXED, DENSITY_TIMES_ROT_PA12, ORTHOTROPIC, DENSITY_TIMES_ORTHOTROPIC, DENSITY_TIMES_2D_TENSOR,
      DENSITY_TIMES_2D_TENSOR_CONSTANT_TRACE, DENSITY_TIMES_ROTATED_2D_TENSOR,D_INTERP_TENSOR, D_INTERP_TENSOR_ROT, LAMINATES, D_LAMINATES, HOM_RECT, D_HOM_RECT, HOM_RECT_C1, HOM_ISO_C1, MSFEM_C1, REDBAS_PARAM, REDBAS_FREE, GREEDY_PARAM, GREEDY_FREE, GREEDY_MAPPING, REDBAS_MAPPING} Type;

    /* possibilities for the isotropic plane in transversal isotropy
    typedef enum { FMO, ISOTROPIC, LAME_ISOTROPIC, TRANSVERSAL_ISOTROPIC, TRANSVERSAL_ISOTROPIC_BOXED, DENSITY_TIMES_TRANSVERSAL_ISOTROPIC,
      DENSITY_TIMES_TRANSVERSAL_ISOTROPIC_BOXED, DENSITY_TIMES_ROT_TRANSVERSAL_ISOTROPIC_BOXED, DENSITY_TIMES_2D_TENSOR,
      DENSITY_TIMES_2D_TENSOR_CONSTANT_TRACE, DENSITY_TIMES_ROTATED_2D_TENSOR, LAMINATES, HOM_RECT,
      REDBAS_PARAM, REDBAS_FREE, GREEDY_PARAM, GREEDY_FREE, GREEDY_MAPPING, REDBAS_MAPPING } Type; */

    /* posibilities for the isotropic plane in transversal isotropy
     * note that parameters EMODULISO, POISSONISO are used for that plane
     * EMODUL is in the orthogonal direction, POSSION is nu_io where i is in the isotropic plane, o not
     * GMODUL is G_io where i is in the isotropic plane o not (note G_io = G_jo) */
    typedef enum { TRANSISO_XY, TRANSISO_YZ, TRANSISO_XZ } TransIsoType;
    
    /** Material notation. Only for FMO we assume the design to be Hill-Mandel, in LinElastInt we use Voigt. The CFS-B-operator is also Voigt, _NO_DENSITY sets topology variable to 1 in simultaneous material and top. opt. */
    typedef enum { VOIGT, HILL_MANDEL, HILL_MANDEL_NO_DENSITY } Notation;

    /** Rotation  direction. Clockwise (CW) or Counter-clockwise (CCW) */
    typedef enum { CW, CCW } Clock;

    /** Method used for interpolation of material tensor and volume */
    typedef enum { NOTYPE=-1, C1, SG, FULL_BSPLINE } Interpolation;

    /** constructor, reads in DesignMaterial from XML
     * @param pn pointer to PtrParamNode */ 
    DesignMaterial(PtrParamNode pn, OptimizationMaterial::System material, StdVector<DesignID>& design, DesignSpace* space);
    
    /** reset the parameter space */
    void ClearParameter() { params_.clear(); }

    /** Set a parameter for the parametric material optimization */
    void SetParameter(const DesignElement::Type p, const double value);

    /** Get a parameter of the parametric material optimization */
    double GetParameter(const DesignElement::Type p) { assert(HasParameter(p)); return params_[p]; }

    /** checks for a parameter */
    bool HasParameter(const DesignElement::Type p) const { return params_.find(p) != params_.end(); }

    /** small helper to return the G-Matrix from model reduction to output as special result */
    bool GetModRedGTensor(Matrix<double>& T, const Elem* elem);

    /** the general material tensor function */
    bool GetTensor(Matrix<double>& t, DesignElement::Type type, SubTensorType subTensor, const Elem* elem, DesignElement::Type direction, DesignMaterial::Notation notation);

    /** Calculate the derivative tensor from the given material parameters */
    bool GetMechTensor(Matrix<double>& t, SubTensorType subTensor, const Elem* elem, DesignElement::Type direction = DesignElement::NO_DERIVATIVE, Notation notation = VOIGT);
    bool GetMechTensor(Matrix<Complex>& t, SubTensorType subTensor, const Elem* elem, DesignElement::Type direction = DesignElement::NO_DERIVATIVE, Notation notation = VOIGT);

    bool GetPiezoCouplingTensor(Matrix<double>& t, const Elem* elem, DesignElement::Type direction);

    /** Calculates MSFEM element matrix for a regular grid from material catalogue*/
    bool GetErsatzElementMatrixMSFEM(Matrix<double>& A,const Elem* elem, DesignElement::Type direction);

    /** helper for GetModRedTensor() but also stand alone to output G Matrix from model reduction as special result */
    void GetModRedGTensor(Matrix<double>& G, DesignElement::Type direction, const bool& all_param);

    /** returns the tensor with negative design variables such the design vector is still pos. definite */
    bool GetElecTensor(Matrix<double>& t,  const Elem* elem, DesignElement::Type direction);

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

    void static SetEnums();

    Type GetType() const { return type_; };

    void SetType(Type type) {type_ = type;};

    /** the actual notation is not stored but assumed as HILL_MANDEL for FMO problems.
     * The enum is necessary for the constraint parameter notation. */
    static Enum<Notation> notation;

    const Elem* current_elem;

    double CalcHomVolume(Vector<double>& p, DesignElement::Type direction, bool derivative);

    Interpolation GetInterpolationMethod() const { return interpolation_; };

    /** rotate elasticity tensor in Voigt notation according to the parameters, eventually calculating a derivative
     *  in 3d: rotates the material by ROTANGLEZ around the z-axis by ROTANGLEY around the y-axis and by ROTANGLEX around the x-axis in this given order or rz,ry,rx
     *  in 2d: rotates the material by ROTANGLE or rx
     * @param t Material Tensor which is rotated in place (or the derivative is calculated in place)
     * @param direction if one of ROTANGLEX, ROTANGLEY, ROTANGLEZ, ROTANGLE calculate the derivative of the rotation w.r.t. this parameter
     * @param notation can be HILL_MANDEL or VOIGT notation
     * @param clock can be CCW (counter-clockwise) or CW (clockwise)
     * @param angles is true if rotation angles rx,ry,rz are given by parameter, otherwise false
     */
    void RotateTensor(Matrix<double>& t, DesignElement::Type direction, Notation notation, Clock clock, bool angles = false, double rx = 0., double ry = 0., double rz = 0.);

    /** Calculate the Isotropic tensor */
    inline void GetIsoMaterialTensor(Matrix<double>& t, SubTensorType subTensor, DesignElement::Type direction);

    /** little helper for GetHomRectTensor(). We assume we are in Hill-Mandel world
       * @param vector p has the values of the design variable */
    void ApplyHomRectC1Tensor(Matrix<double>& E, Vector<double>& p, DesignElement::Type direction, SubTensorType subTensor) const;

    /** little helper for GetHomRectTensor(). We assume we are in Hill-Mandel world
       * @param vector p has the values of the design variable */
    void ApplyHomIsoC1Tensor(Matrix<double>& E, Vector<double>& p, DesignElement::Type direction, SubTensorType subTensor) const;

  protected:

    /** for debugging */
    void DumpParams();

    /** Sets all Material Parameters in designMaterial for given element in the current context */
    bool CollectMaterialParametersForElement(DesignSpace* space, const Elem* elem);

    std::map<DesignElement::Type, double> params_;

    /** mass is considered an independent design */
    bool massIsDesign_;

    /** damping is also optimized */
    bool dampingIsDesign_;

    /** shearing is optimized */
    bool shearIsDesign_;

    /** dimension of material catalogue */
    StdVector<double> catalogueSize_;

    /** multiply mass with this, can be used to scale tensor trace */
    double massFactor_;

    /** for density times 2d tensor with constant trace */
    double trace_;

    static Enum<Type> type;
    Type type_;

    static Enum<TransIsoType> transIsoType;
    TransIsoType transIsoType_;

    unsigned int dim;

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


    /** Calculate the Lame Tensor */
    inline void GetLameMaterialTensor(Matrix<double>& t, SubTensorType subTensor, DesignElement::Type direction);

    /** Calculate the Trans-Iso Tensor */
    inline void GetTransIsoMaterialTensor(Matrix<double>& t, SubTensorType subTensor, DesignElement::Type direction, Notation notation);
    
    /* general anisotropic FMO tensor */
    inline void GetElasticFMOTensor(Matrix<double>& t, DesignElement::Type direction, Notation notation);

    /** Calculate the orthotropic material tensor */
    inline void GetOrthotropicMaterialTensor(Matrix<double>& t, SubTensorType subTensor, DesignElement::Type direction, Notation notation);

    /** Calculate the Tensor for Density times Tensor */
    inline void GetDensityTimes2dTensorTensor(Matrix<double>& t, SubTensorType subTensor, DesignElement::Type direction);

    /** Calculate the tensor for Laminates */
    inline void GetLaminatesTensor(Matrix<double>& t, SubTensorType subTensor, DesignElement::Type direction, Notation notation);

    /** little helper for GetHomRectTensor(). We assume we are in Hill-Mandel world
     * @param shape might also be the x or y component of the derivative! */
    void ApplyHomRectTensor(Matrix<double>& E, const Vector<double>& shape) const;

    /** Approximates the homogenized tensor of an a-b rectangle as used by Bendsoe and Kikuchi 1988 */
    inline void GetHomRectTensor(Matrix<double>& t, SubTensorType subTensor,  const Elem* elem,  DesignElement::Type direction, Notation notation);

    /**Computes the homogenized tensor from the reduced-order model obtaind for the homogenization formula */
    inline void GetModRedTensor(Matrix<double>& t, DesignElement::Type direction, Notation notation);

    void GetMappingTensor(Matrix<double>& E, DesignElement::Type direction, Notation notation);

    //Computes the gradient of the mapping inside the element considered
    void GetMappingGradient(Matrix<double>& G);

    //Computes the gradient of the gradient of the mapping inside the element considered with respect to all the variables considered
    void GetMappingGradient(Matrix<double>& G, DesignElement::Type direction);

    /**Computes the homogenized tensor from the reduced-order model obtained for the homogenization formula with the greedy algorithm*/
    inline void GetGreedyTensor(Matrix<double>& t, DesignElement::Type direction, Notation notation);

    /** Approximates the homogenized tensor of an a-b rectangle as used by Bendsoe and Kikuchi 1988 */
    inline void GetInterpolatedTensor(Matrix<double>& t, SubTensorType subTensor,  DesignElement::Type direction, Notation notation);

    /** does only perform orientational optimization
     * @param mc MECHANIC, PIEZO, ELECTROSTATIC */
    inline void GetRotatedTensor(Matrix<double>& t, MaterialClass mc, DesignElement::Type direction);

    /** initialize the tensor with zeros */
    inline void ZeroTensor(Matrix<double>& t, SubTensorType subTensor);
    
    /** put values from Voigt vector to correct positions in tensor */
    inline void Set2dVoigtTensor(Matrix<double>& t, double t11, double t22, double t33, double t23, double t13, double t12);
    
    /** put values from Voigt vector to correct positions in tensor (doesn't assume symmetry) */
    inline void Set2dVoigtTensor(Matrix<double>& t, double t11, double t12, double t13, double t21, double t22, double t23, double t31, double t32, double t33);

    /** put the entries of the orthotropic tensor at the right places */
    inline void SetOrthotropicTensor(Matrix<double>& t, SubTensorType subTensor, double e11, double e12, double e13, double e22,
        double e23, double e33, double e44, double e55, double e66);

    /** put the entries of the transversal_isotropic tensor at the right places */
    inline void SetTransIsoTensor(Matrix<double>& t, SubTensorType subTensor, double iD, double inD, double iG, double oD, double onD, double oG);
    
    /** put the entries of the isotropic tensor at the right places */
    inline void SetIsoTensor(Matrix<double>& t, SubTensorType subTensor, double D, double nD, double G);

    void RotateHMStiffnessTensor(Matrix<double>& t, SubTensorType subTensor, DesignElement::Type direction, double a, Notation notation = HILL_MANDEL);

    /** helper function to set a rotation matrix of size 3x3
     * the matrix (when calculating R*x) would rotate the vector x by thetaz around the z-axis by thetay around the y-axis and by thetax around the x-axis in this given order
     * @param R the place to set the rotation matrix
     * @sthetax sin(thetax)
     * @cthetax cos(thetax)
     * @sthetay sin(thetax)
     * @cthetay cos(thetax)
     * @sthetaz sin(thetax)
     * @cthetaz cos(thetax)
     * @direction if given direction of the derivative to be calculated
     */
    void SetRotationMatrix(Matrix<double>& R, double sthetax, double cthetax, double sthetay, double cthetay, double sthetaz, double cthetaz, DesignElement::Type direction = DesignElement::NO_DERIVATIVE);

    /** This exists only in Voigt notation! */
    void RotatePiezoCouplingTensor(Matrix<double>& t, double angle, DesignElement::Type direction);

    void RotateElecTensor(Matrix<double>& t, double angle, DesignElement::Type direction);

    /** Calculate the mass isotropic case */
    inline double GetIsoMaterialMass(DesignElement::Type direction);    
    
    /** Calculate the mass lame case */
    inline double GetLameMaterialMass(DesignElement::Type direction);
    
    /** Calculate the mass trans-iso case */
    inline double GetTransIsoMaterialMass(DesignElement::Type direction);

    /** Calculate the mass density-times-tensor case
     * This returns the scaling factor (pseudo-density) for the normal mass matrix based on the materials actual density */
    inline double GetDensityTimesTensorMass(DesignElement::Type direction);

    /** Get the trans-iso mass (tensor trace) out of the corresponding tensor entries */
    inline double GetTransIsoMass(double iD, double iG, double oD, double oG);

    /** Get the isotropic mass (tensor trace) out of the corresponding tensor entries */
    inline double GetIsoMass(double D, double G);

    /** fills the row in hom_rect_samples_ */
    void FillHomRectSamples(PtrParamNode homRect, unsigned int idx, const std::string& a, const std::string& b);

    /** fills the matrices in mod_red_matrices_ **/
    void FillModRedMatrices(PtrParamNode matnode, const StdVector<std::string>& tensor_comp, const int& tensor_int, const UInt& dimbas);
    void FillModRedMatrices(PtrParamNode matnode, const StdVector<std::string>& tensor_comp, const int& tensor_int, const UInt& dimbas, const UInt& dimbastot);

    /** fills the vectors in mod_red_vectors_ **/
    void FillModRedVectors(PtrParamNode vecnode, const StdVector<std::string>& tensor_comp, const int& tensor_int, const UInt& dimbas);
    void FillModRedVectors(PtrParamNode vecnode, const StdVector<std::string>& tensor_comp, const int& tensor_int, const UInt& dimbas, const UInt& dimbastot);

    //Returns the homogenized elasticity tensor associated to the matrix G
    void GetModRedHomTensor(Matrix<double>& E, const Matrix<double>& G, const StdVector<Vector<double> >& corrector_, Notation notation);

    //Returns the derivative of the homogenized elasticity tensor associated with a matrix G and its derivative Gderiv
    void GetModRedHomTensor(Matrix<double>& E, const Matrix<double>& G, const Matrix<double>& Gderiv, const StdVector<Vector<double> >& corrector_, Notation notation);

    /** gives the SVD parameters of the 2*2 matrix G ***/
    void GetSVDGTensorParameters(const Matrix<double>& G, Vector<double>& paramvec);

    //Get the values of the parameters
    void GetModRedParamVector(Vector<double>& params);

    //Solves the corrector problems using the reduced basis
    void GetRedBasCorrector(StdVector<Vector<double> >& corrector_, const Matrix<double>& G);

    //Get the Corrector for the greedy model reduction
    void GetGreedyCorrector(StdVector<Vector<double> >& corrector_, const Vector<double>& params, const bool& all_param);

    //Compute the value of angle-dependent functions in the greedy case
    double AngleGreedyCalculus(const Vector<double>& coeffs, const double& angle);

    //Compute the value of scaling-dependent functions in the greedy caseS
    double ScalingGreedyCalculus(const Vector<double>& coeffs, const double& l);

    /** fills the coefficient data structure for the bicubic interpolation*/
    void FillHomRectCoeff(Matrix<double> & coeff_,const char * filename);

    /** evaluates the C1 interpolation polynomial at point p[0],p[1] and returns function value as double */
    double EvaluateC1Interpolation(Vector<double>& p, const Matrix<double>& coeff, double& da, double& db, int& j, int& k, int& m, int& n) const;

    /** evaluates the derivative of the C1 interpolation polynomial at point p[0],p[1] in direction 0 or 1 and returns function value as double */
    double EvaluateC1Interpolation_Deriv(Vector<double>& p, const Matrix<double>& coeff, double& da, double& db, int& j, int& k, int& m, int& n, DesignElement::Type direction) const;

    double EvaluateC1Interpolation_3D(Vector<double>& p, const Matrix<double>& coeff, double& da, double& db, double& dc, int& j, int& k,int& l, int& m, int& n, int& o) const;

    /** evaluates the derivative of the C1 interpolation polynomial at point p[0],p[1],p[2] in direction 0 or 1 and returns function value as double */
    double EvaluateC1Interpolation_Deriv_3D(Vector<double>& p, const Matrix<double>& coeff, double& da, double& db,double& dc, int& j, int& k, int& l, int& m, int& n, int& o, DesignElement::Type direction) const;
    //double EvaluateC1Interpolation(Matrix<double>& E,  Vector<double>& p, const Matrix<double> & coeff, int au,int al,int bu,int bl,int j, int k,int m,int n);

    /** Get the index of the local interpolation interval*/
    int GetInterpolationIndex(Matrix<double> interval, double& point) const;

    /** Read detailed stats from file*/
    bool ReadDetailedStats(const char * filename, Matrix<double>& ret);

#ifdef USE_SGPP
    /** little helper for GetHomRectTensor(). We assume we are in Hill-Mandel world
       * @param vector p has the values of the design variable */
    void ApplyHomRectSGPPTensor(Matrix<double>& E, Vector<double>& p, DesignElement::Type direction, SubTensorType subTensor);

    /** little helper for GetHomRectTensor(). We assume we are in Hill-Mandel world
       * @param vector p has the values of the design variable */
    void ApplyHomRectFullBsplineTensor(Matrix<double>& E, Vector<double>& p, DesignElement::Type direction, SubTensorType subTensor) const;

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

    /** sampled values for a single hom-rect 9-element by the number of shape function. Notation is Hill-Mandel!
     * 9 rows and 6 columns for with TENSOR11 being the first */
    Matrix<double> hom_rect_samples_;
    /** sampled values for coefficients of the bicubic interpolation polynomial; number of sample elements rows and 16 columns/64 columns (3D)*/
    Matrix<double> hom_rect_coeff11_;
    Matrix<double> hom_rect_coeff12_;
    Matrix<double> hom_rect_coeff22_;
    Matrix<double> hom_rect_coeff33_;
    Matrix<double> hom_rect_coeff23_;
    Matrix<double> hom_rect_coeff44_;
    Matrix<double> hom_rect_coeff55_;
    Matrix<double> hom_rect_coeff66_;
    Matrix<double> hom_rect_coeff13_;
    Matrix<double> hom_rect_coeff14_;
    Matrix<double> hom_rect_coeff15_;
    Matrix<double> hom_rect_coeff16_;
    Matrix<double> hom_rect_coeff24_;
    Matrix<double> hom_rect_coeff25_;
    Matrix<double> hom_rect_coeff26_;
    Matrix<double> hom_rect_coeff34_;
    Matrix<double> hom_rect_coeff35_;
    Matrix<double> hom_rect_coeff36_;
    Matrix<double> hom_rect_coeff45_;
    Matrix<double> hom_rect_coeff46_;
    Matrix<double> hom_rect_coeff56_;
    Matrix<double> hom_rect_a_;
    Matrix<double> hom_rect_b_;
    Matrix<double> hom_rect_c_;

    /** MSFEM element matrix coefficients of the bi-/tricubic interpolation polynomial from material catalogue; number of sample elements rows and 64 columns */
    Matrix<double> msfem_a_;
    Matrix<double> msfem_b_;
    Matrix<double> msfem_rot_;
    StdVector<Matrix<double> > msfem_coeff_;

    //** Contains the matrices and vectors with the information for the model reduction case (reduced basis or greedy)
    UInt dimension_;
    UInt dimension_tot_;

    //The matrices and vectors of the reduced model should be given in Voigt notation
    StdVector<Matrix<double> > mod_red_matrices_;

    StdVector<Vector<double> > mod_red_vectors_;

    //Mean_tensor_ = [E11, E12, E33];
    Matrix<double> mean_tensor_;

    //Gives the information: do we treat the rotation angle theta just like another parameter or not
    bool all_param_;


    //Contains information for the greedy case
    UInt Na_;
    UInt Nl_;
    double lmin_;
    double lmax_;

    //Contains the infomation about the parameters for the corrector problem in the greedy case
    StdVector<Matrix<double> > matrices_param_;

    DesignSpace* space_;

    Interpolation interpolation_;
    unsigned int level_;

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
