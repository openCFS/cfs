/*
 * Model.hh
 *
 *  Created on: Aug 11, 2021
 *      Author: alex
 */

#ifndef SOURCE_MATERIALS_MODELS_MODEL_HH_
#define SOURCE_MATERIALS_MODELS_MODEL_HH_

#include <list>

#include <boost/utility.hpp>

#include "MatVec/Vector.hh"
#include "Domain/Domain.hh"
#include "Domain/ElemMapping/EntityLists.hh"

namespace CoupledField {

class MathParser;

class Model {

public:
  //! Constructor
  Model();

  //! Destructor
  virtual ~Model();

  virtual void Init(std::map<std::string, double> ParameterMap, UInt numElems, UInt dim){
    EXCEPTION( "Not implemented in base class");
  };

  virtual void Init(std::map<std::string, double> ParameterMap, shared_ptr<ElemList> entityList, UInt dim){
    EXCEPTION( "Not implemented in base class");
  };

  virtual Double ComputeMaterialParameter(Vector<Double> E, Integer ElemNum){
    EXCEPTION( "Not implemented in base class");
  };

  virtual Matrix<Double> ComputeTensorialMaterialParameter(Vector<Double> E, Integer ElemNum){
    EXCEPTION( "Not implemented in base class");
  };
  
  virtual Vector<Double> GetFluxDensity(Vector<Double> E, Integer ElemNum){
    EXCEPTION( "Not implemented in base class");
  };
 
  virtual Vector<Double> GetFluxDensity(Vector<Double> E, Integer ElemNum,
                                        LocPointMapped lpm, PtrCoefFct stressCoef){
    EXCEPTION( "Not implemented in base class");
  };

  virtual void UpdateStates(){
    EXCEPTION( "Not implemented in base class");
  }
};

}

#endif /* SOURCE_MATERIALS_MODELS_MODEL_HH_ */
