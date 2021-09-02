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
#include "Utils/mathParser/mathParser.hh"
#include "Domain/Domain.hh"

namespace CoupledField {

class Model {

public:
  //! Constructor
  Model();

  //! Destructor
  virtual ~Model();

  virtual void Init(std::map<std::string, double> ParameterMap, UInt numElems){
    EXCEPTION( "Not implemented in base class");
  };

  virtual Double ComputeMaterialParameter(Vector<Double> E, Integer ElemNum){
    EXCEPTION( "Not implemented in base class");
  };

};

}

#endif /* SOURCE_MATERIALS_MODELS_MODEL_HH_ */
