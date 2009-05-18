#ifndef DESIGNMATERIAL_HH_
#define DESIGNMATERIAL_HH_

#include "General/Enum.hh"
#include "DesignElement.hh"
#include "Forms/baseForm.hh"
#include <map>

namespace CoupledField {

  /** This implements a function from $R^n$ to $R^{d \times d} for transforming a vector of Parameters
   * to a material tensor.  */
  class DesignMaterial {
    
  public:
    
    typedef enum { ISOTROPIC, LAME_ISOTROPIC, TRANSVERSAL_ISOTROPIC, TRANSVERSAL_ISOTROPIC_BOXED } Type;
    
    /* posibilities for the isotropic plane in transversal isotropy
     * note that parameters EMODULISO, POISSONISO are used for that plane
     * EMODUL is in the orthogonal direction, POSSION is nu_io where i is in the isotropic plane, o not
     * GMODUL is G_io where i is in the isotropic plane o not (note G_io = G_jo) */
    typedef enum { TRANSISO_XY, TRANSISO_YZ, TRANSISO_XZ } TransIsoType;
    
    /** constructor, reads in DesignMaterial from XML
     * @param pn pointer to ParamNode */ 
    DesignMaterial(ParamNode* pn);
    
    /** returns the numbers of parameters requiered for this material */
    unsigned int RequiredParameters();    
    
    /** returns the number of required designs for this material 
     * note: no checking except the number is currently performed */
    unsigned int RequiredDesigns(){
      return RequiredParameters() - nparams_;
    }
    
    /** Set a parameter for the parametric material optimization */
    void SetParameter(const DesignElement::Type p, const double value);
    
    /** Calculate the elasticity tensor from the given material parameters */
    void GetMaterialTensor(Matrix<double>& t, SubTensorType subTensor);
    
    /** Calculate the derivative tensor from the given material parameters */
    void GetMaterialTensorDerivative(Matrix<double>& t, SubTensorType subTensor, DesignElement::Type direction);
    
    /** initialize the tensor with zeros */
    void ZeroTensor(Matrix<double>& t, SubTensorType subTensor);
    
    /** put the entries of the transversal_isotropic tensor at the right places */
    void SetTransIsoTensor(Matrix<double>& t, SubTensorType subTensor, double iD, double inD, double iG, double oD, double onD, double oG);
    
    /** put the entries of the isotropic tensor at the right places */
    void SetIsoTensor(Matrix<double>& t, SubTensorType subTensor, double D, double nD, double G);
    
    /** Conveniance method */
    void GetMaterialTensor(Matrix<double>& t, SubTensorType subTensor, DesignElement::Type direction){
      if(direction == DesignElement::NO_DERIVATIVE){
        GetMaterialTensor(t, subTensor);
      }else{
        GetMaterialTensorDerivative(t, subTensor, direction);
      }
    }
    
    void static SetEnums();
    
  protected:
    std::map<DesignElement::Type, double> params_;
    unsigned int nparams_;
    
    static Enum<Type> type;
    Type type_;   
    static Enum<TransIsoType> transIsoType;
    TransIsoType transIsoType_;
  };

} // namespace

#endif /*DESIGNMATERIAL_HH_*/
