// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef COEF_FUNCTION_JILESCACHE_HH
#define COEF_FUNCTION_JILESCACHE_HH

#include <string>
#include <map>

#include "CoefFunction.hh"
#include "FeBasis/BaseFE.hh"
#include "FeBasis/FeFunctions.hh"

#include "Materials/BaseMaterial.hh"
#include "Materials/Models/Jiles.hh"
#include "Materials/Models/EBHysteresis.hh"
#include "Materials/Models/invEBHysteresis.hh"
#include "Materials/Models/Model.hh"

namespace CoupledField  {

class MathParser;

/* Idea behind this ceoffunction:
 *
 * This class provides a interface between the actual model (eg Jiles-Model) and the PDE.
 * If there is a model used for some materialparameter, this class is created inside the PDE and
 * handles the model.
 *
 * So whenever the PDE needs to evaluate the model on a given local point, it calls this object and this coeffunction
 * calls the model.
 *
 * To implement new models, one have to create a new Model in /Materials/Models and the new model has the
 * heritate from the class model.hh/cc and provide those methods. Through the XML-Input the PDE nows what
 * Model should be used (modelName_). The model should do most of the work, eg. it should save all the
 * previous values it needs to evaluate, and it should know the current timestep (mathparser!).
 *
 * In the method Init there should then be a if-statement whichs inits the given model.
 */
template <class TYPE>
class CoefFunctionMaterialModel : public CoefFunction {

public:
  //! Constructor
  CoefFunctionMaterialModel();

  //! Destructor
  virtual ~CoefFunctionMaterialModel();

  virtual string GetName() const { return "CoefFunctionMaterialModel"; }

  // Init this coeffunction. State the dependend Coeffunction and used model
  void Init(PtrCoefFct depCoef, std::string modelName, UInt dim=3);

  // Evaluates model at certain lpm. (at what thime should be implemented in the model)
  void GetScalar(Double& coefScalar, const LocPointMapped& lpm);

  // Evaluates model at certain lpm
  void GetTensor(Matrix<Double>& coefVector, const LocPointMapped& lpm);

  // Evaluates model at certain lpm
  void GetVector(Vector<Double>& coefVector, const LocPointMapped& lpm);

  //! Return row and columns size of tensor if coefficient function is a tensor
  void GetTensorSize( UInt& numRows, UInt& numCols ) const {
    numRows = 2;
    numCols = 2;
  }

  UInt GetVecSize() const {
    return 2;
  }

  // Init  the used material model
  void InitModel(std::map<std::string, double> ParameterMap, UInt numElems);

  void InitModel(std::map<std::string, double> ParameterMap, shared_ptr<ElemList> entityList);

protected:

  // Spatial dimension of the problem
  UInt spaceDim_;
  
  // object for the model
  Model* matModel_;

  // modelname
  std::string modelName_ = "none";

  //Coeffunction, which is used to evalute the model
  PtrCoefFct depCoef_ ;
};
} //end of namespace

#endif
