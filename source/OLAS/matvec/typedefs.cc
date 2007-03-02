
#include "matvec/typedefs.hh"

namespace OLAS {

  // Scalar matrices and vectors
template<> std::string assocType<Double >::tagM = "Double";
template<> std::string assocType<Double >::tagV = "Double";
template<> std::string assocType<Double >::tagS = "Double";
template<> std::string assocType<Complex>::tagM = "Complex";
template<> std::string assocType<Complex>::tagV = "Complex";
template<> std::string assocType<Complex>::tagS = "Complex";


}

