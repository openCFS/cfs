#ifndef FILE_TETRAFE
#define FILE_TETRAFE

#include <Elements/basefe.hh>

namespace CoupledField
{

//! Class with general description of tetrahedral element
/*! This class is derived from BaseFE. 
    It stores general procedures for each type of finite element based on tetrahedral, such as information 
    about integration points and integration weights
*/

class TetraFE : public BaseFE
{
public:
   //! Constructor with type of integration rule
   TetraFE();

   //! Deconstructor
  virtual ~TetraFE(); 

  //! return FE-Type for LAS++
#ifdef USE_OLAS
  virtual FEType feType() { return TET;}
#else
  virtual Integer feType() { return 3;}
#endif
  

protected:
  //! Set integration points
  virtual void SetIntPoints();

  /// Corrects the integration weights due to volume unequal 1
  virtual void CorrectIntWeights();  
};

}
#endif //
