#ifndef DESIGNMATERIAL_HH_
#define DESIGNMATERIAL_HH_

#include <map>

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DesignElement.hh"
#include "General/Enum.hh"
#include "General/environment.hh"

namespace CoupledField {

  /** This implements a function from $R^n$ to $R^{d \times d}$ for transforming a vector of Parameters
   * to a material tensor.  */
template <class TYPE> class Matrix;
template <class TYPE> class StdVector;

  class DesignMaterial {
    
  public:
    
    typedef enum { ISOTROPIC, LAME_ISOTROPIC, TRANSVERSAL_ISOTROPIC, TRANSVERSAL_ISOTROPIC_BOXED, DENSITY_TIMES_TRANSVERSAL_ISOTROPIC,
      DENSITY_TIMES_TRANSVERSAL_ISOTROPIC_BOXED, DENSITY_TIMES_2D_TENSOR, DENSITY_TIMES_2D_TENSOR_CONSTANT_TRACE } Type;
    
    /* posibilities for the isotropic plane in transversal isotropy
     * note that parameters EMODULISO, POISSONISO are used for that plane
     * EMODUL is in the orthogonal direction, POSSION is nu_io where i is in the isotropic plane, o not
     * GMODUL is G_io where i is in the isotropic plane o not (note G_io = G_jo) */
    typedef enum { TRANSISO_XY, TRANSISO_YZ, TRANSISO_XZ } TransIsoType;
    
    /** constructor, reads in DesignMaterial from XML
     * @param pn pointer to PtrParamNode */ 
    DesignMaterial(PtrParamNode pn, StdVector<DesignElement::Type>& design);
    
    /** Set a parameter for the parametric material optimization */
    void SetParameter(const DesignElement::Type p, const double value);

    /** Get a parameter of the parametric material optimization */
    double GetParameter(const DesignElement::Type p);
    
    /** Calculate the derivative tensor from the given material parameters */
    void GetMaterialTensor(Matrix<double>& t, SubTensorType subTensor, DesignElement::Type direction = DesignElement::NO_DERIVATIVE);
    
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
    
  protected:
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
    unsigned int RequiredParameters();    
    
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
    
    /** Calculate the Tensor for Density times Tensor */
    inline void GetDensityTimes2dTensorTensor(Matrix<double>& t, SubTensorType subTensor, DesignElement::Type direction);
    

    /** initialize the tensor with zeros */
    inline void ZeroTensor(Matrix<double>& t, SubTensorType subTensor);
    
    /** put values from Voigt vector to correct positions in tensor */
    inline void Set2dVoigtTensor(Matrix<double>& t, SubTensorType subTensor, double t11, double t22, double t33, double t23, double t13, double t12);
    
    /** put the entries of the transversal_isotropic tensor at the right places */
    inline void SetTransIsoTensor(Matrix<double>& t, SubTensorType subTensor, double iD, double inD, double iG, double oD, double onD, double oG);
    
    /** put the entries of the isotropic tensor at the right places */
    inline void SetIsoTensor(Matrix<double>& t, SubTensorType subTensor, double D, double nD, double G);
    

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

  };

} // namespace

#endif /*DESIGNMATERIAL_HH_*/
