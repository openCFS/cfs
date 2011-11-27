
#include "Domain/CoefFunction/CoefFunction.hh"

namespace CoupledField{

shared_ptr<CoefFunction> 
CoefFunction::Generate( Global::ComplexPart format, 
                        const std::string& realVal, 
                        const std::string& imagVal  ) {
  EXCEPTION( "Implement me");
  shared_ptr<CoefFunction> ret;
  return ret;

}
shared_ptr<CoefFunction> 
CoefFunction::Generate( Global::ComplexPart format, 
                        const StdVector<std::string>& realVal, 
                        const StdVector<std::string>& imagVal ) {
  EXCEPTION( "Implement me");
  shared_ptr<CoefFunction> ret;
  return ret;
}


//! Generate tensor-valued coefficient function
shared_ptr<CoefFunction> 
CoefFunction::Generate( Global::ComplexPart format,
                        UInt numRows, UInt numCols,
                        const StdVector<std::string>& realVal,
                        const StdVector<std::string>& imagVal ) {
  EXCEPTION( "Implement me");
  shared_ptr<CoefFunction> ret;
  return ret;
}
bool CoefFunction::ExprDependsOnTimeFreq(const std::string& expr) {
  EXCEPTION( "Implement me");
}

//! Returns true, if expression depends on space
bool CoefFunction::ExprDependsOnSpace(const std::string& expr) {
  EXCEPTION( "Implement me");
}
}
