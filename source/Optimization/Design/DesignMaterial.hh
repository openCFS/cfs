#ifndef DESIGNMATERIAL_HH_
#define DESIGNMATERIAL_HH_

#include <map>
#include <complex>

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DesignElement.hh"
#include "General/Enum.hh"
#include "General/environment.hh"
#include "Optimization/OptimizationMaterial.hh"
#include "MatVec/matrix.hh"

namespace CoupledField {

  /** This implements a function from $R^n$ to $R^{d \times d}$ for transforming a vector of Parameters
   * to a material tensor.  */
template <class TYPE> class StdVector;

  class DesignMaterial {
    
  public:
    
    typedef enum { FMO, ISOTROPIC, LAME_ISOTROPIC, TRANSVERSAL_ISOTROPIC, TRANSVERSAL_ISOTROPIC_BOXED, DENSITY_TIMES_TRANSVERSAL_ISOTROPIC,
      DENSITY_TIMES_TRANSVERSAL_ISOTROPIC_BOXED, DENSITY_TIMES_ROT_TRANSVERSAL_ISOTROPIC_BOXED, DENSITY_TIMES_2D_TENSOR,
      DENSITY_TIMES_2D_TENSOR_CONSTANT_TRACE, DENSITY_TIMES_ROTATED_2D_TENSOR, LAMINATES, HOM_RECT, MODEL_REDUCTION, MODEL_REDUCTION2 } Type;
    
    /* posibilities for the isotropic plane in transversal isotropy
     * note that parameters EMODULISO, POISSONISO are used for that plane
     * EMODUL is in the orthogonal direction, POSSION is nu_io where i is in the isotropic plane, o not
     * GMODUL is G_io where i is in the isotropic plane o not (note G_io = G_jo) */
    typedef enum { TRANSISO_XY, TRANSISO_YZ, TRANSISO_XZ } TransIsoType;
    
    /** Material notation. Only for FMO we assume the design to be Hill-Mandel, in LinElastInt we use Voigt. The CFS-B-operator is also Voigt */
    typedef enum { VOIGT, HILL_MANDEL } Notation;

    /** constructor, reads in DesignMaterial from XML
     * @param pn pointer to PtrParamNode */ 
    DesignMaterial(PtrParamNode pn, OptimizationMaterial::System material, StdVector<DesignID>& design);
    
    /** reset the parameter space */
    void ClearParameter() { params_.clear(); }

    /** Set a parameter for the parametric material optimization */
    void SetParameter(const DesignElement::Type p, const double value);

    /** Get a parameter of the parametric material optimization */
    double GetParameter(const DesignElement::Type p) { assert(HasParameter(p)); return params_[p]; }
    
    /** checks for a parameter */
    bool HasParameter(const DesignElement::Type p) const { return params_.find(p) != params_.end(); }

    /** Calculate the derivative tensor from the given material parameters */
    void GetMaterialTensor(Matrix<double>& t, SubTensorType subTensor, DesignElement::Type direction = DesignElement::NO_DERIVATIVE, Notation notation = VOIGT);

    /** helper for GetModRedTensor() but also stand alone to output G Matrix from model reduction as special result */
    void GetModRedGTensor(Matrix<double>& G, DesignElement::Type direction);

    void GetPiezoCouplingTensor(Matrix<double>& t, DesignElement::Type direction);

    /** returns the tensor with negative design variables such the design vector is still pos. definite */
    void GetDielecTensor(Matrix<double>& t, DesignElement::Type direction);

    /** retrieve rel. mass of element (tensor trace) or derivative thereof */
    double GetMaterialMass(DesignElement::Type direction);
    
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
    
    Type GetType() const { return type_; }

    /** the actual notation is not stored but assumed as HILL_MANDEL for FMO problems.
     * The enum is necessary for the constraint parameter notation. */
    static Enum<Notation> notation;

  protected:

    /** for debugging */
    void DumpParams();

    std::map<DesignElement::Type, double> params_;
   
    /** mass is considered an independent design */
    bool massIsDesign_;
    
    /** damping is also optimized */
    bool dampingIsDesign_;
    
    /** multiply mass with this, can be used to scale tensor trace */
    double massFactor_;
    
    /** for density times 2d tensor, this is the penalization for density*/
    double penalty_;
    
    /** for density times 2d tensor with constant trace */
    double trace_;
    
    static Enum<Type> type;
    Type type_;   

    static Enum<TransIsoType> transIsoType;
    TransIsoType transIsoType_;
    
    unsigned int dim;
    
    /** returns the numbers of parameters required for this material */
    unsigned int RequiredParameters(OptimizationMaterial::System material);
    
    /** Check whether all required designs are available */
    bool CheckRequiredDesigns(StdVector<DesignElement::Type>& design);
    
  private:
    /* note that most of these functions are called really often, so inlining is used */

    /** Calculate the Isotropic tensor */
    inline void GetIsoMaterialTensor(Matrix<double>& t, SubTensorType subTensor, DesignElement::Type direction);

    /** Calculate the Lame Tensor */
    inline void GetLameMaterialTensor(Matrix<double>& t, SubTensorType subTensor, DesignElement::Type direction);

    /** Calculate the Trans-Iso Tensor */
    inline void GetTransIsoMaterialTensor(Matrix<double>& t, SubTensorType subTensor, DesignElement::Type direction);
    
    /* general anisotropic FMO tensor */
    inline void GetAnisotropicTensor(Matrix<double>& t, DesignElement::Type direction, Notation notation);

    /** Calculate the Tensor for Density times Tensor */
    inline void GetDensityTimes2dTensorTensor(Matrix<double>& t, SubTensorType subTensor, DesignElement::Type direction);
    
    /** Calculate the tensor for Laminates */
    inline void GetLaminatesTensor(Matrix<double>& t, SubTensorType subTensor, DesignElement::Type direction, Notation notation);

    /** little helper for GetHomRectTensor(). We assume we are in Hill-Mandel world
     * @param shape might also be the x or y component of the derivative! */
    void ApplyHomRectTensor(Matrix<double>& E, const Vector<double>& shape) const;

    /** Approximates the homogenized tensor of an a-b rectangle as used by Bendsoe and Kikuchi 1988 */
    inline void GetHomRectTensor(Matrix<double>& t, DesignElement::Type direction, Notation notation);

    /**Computes the homogenized tensor from the reduced-order model obtaind for the homogenization formula */
    inline void GetModRedTensor(Matrix<double>& t, DesignElement::Type direction, Notation notation);

    /** initialize the tensor with zeros */
    inline void ZeroTensor(Matrix<double>& t, SubTensorType subTensor);
    
    /** put values from Voigt vector to correct positions in tensor */
    inline void Set2dVoigtTensor(Matrix<double>& t, double t11, double t22, double t33, double t23, double t13, double t12);
    
    /** put the entries of the transversal_isotropic tensor at the right places */
    inline void SetTransIsoTensor(Matrix<double>& t, SubTensorType subTensor, double iD, double inD, double iG, double oD, double onD, double oG);
    
    /** put the entries of the isotropic tensor at the right places */
    inline void SetIsoTensor(Matrix<double>& t, SubTensorType subTensor, double D, double nD, double G);
    
    /** rotate elasticity tensor in t (in Hill-Mandel notation!) by the angle a and adjust the entries back to notation to fit with CFS++ */
    inline void RotateHMStiffnessTensor(Matrix<double>& t, SubTensorType subTensor, DesignElement::Type direction, double a, Notation notation = VOIGT);


    /** Calculate the mass isotropic case */
    inline double GetIsoMaterialMass(DesignElement::Type direction);    
    
    /** Calculate the mass lame case */
    inline double GetLameMaterialMass(DesignElement::Type direction);
    
    /** Calculate the mass trans-iso case */
    inline double GetTransIsoMaterialMass(DesignElement::Type direction);
    

    /** Get the trans-iso mass (tensor trace) out of the corresponding tensor entries */
    inline double GetTransIsoMass(double iD, double iG, double oD, double oG);
    
    /** Get the isotropic mass (tensor trace) out of the corresponding tensor entries */
    inline double GetIsoMass(double D, double G);

    /** fills the row in hom_rect_samples_ */
    void FillHomRectSamples(PtrParamNode homRect, unsigned int idx, const std::string& a, const std::string& b);


    /** fills the matrices in mod_red_matrices_ **/
    void FillModRedMatrices(PtrParamNode matnode, const StdVector<std::string>& tensor_comp, const int& tensor_int);

    /** fills the vectors in mod_red_vectors_ **/
    void FillModRedVectors(PtrParamNode vecnode, const StdVector<std::string>& tensor_comp, const int& tensor_int);

    //Solves the corrector problems using the reduced basis
    void GetModRedCorrector(StdVector<Matrix<double> >& corrector_, const Matrix<double>& G);

    //Returns the homogenized elasticity tensor associated to the matrix G
    void GetModRedHomTensor(Matrix<double>& E, const Matrix<double>& G, const StdVector<Matrix<double> >& corrector_, Notation notation);

    //Returns the derivative of the homogenized elasticity tensor associated with a matrix G and its derivative Gderiv
    void GetModRedHomTensor(Matrix<double>& E, const Matrix<double>& G, const Matrix<double>& Gderiv, const StdVector<Matrix<double> >& corrector_, Notation notation);

    /** sampled values for a single hom-rect 9-element by the number of shape function. Notation is Hill-Mandel!
     * 9 rows and 6 columns for with TENSOR11 being the first */
    Matrix<double> hom_rect_samples_;

    //** Contains the matrices and vectors with the reduced basis information for the model reduction case
    int dimension_;

    //The matrices and vectors of teh reduced model shoul be given in Voigt notation
    StdVector<Matrix<double> > mod_red_matrices_;

    StdVector<Matrix<double> > mod_red_vectors_;

    //Mean_tensor_ = [E11, E12, E33];
    Matrix<double> mean_tensor_;


  };

} // namespace

#endif /*DESIGNMATERIAL_HH_*/
