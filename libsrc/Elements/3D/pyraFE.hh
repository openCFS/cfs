#ifndef FILE_PYRAFE
#define FILE_PYRAFE

#include <Elements/basefe.hh>

namespace CoupledField
{

//! Class with general description of pyramidal element
/*! This class is derived from BaseFE. 
    It stores general procedures for each type of finite element based on pyramidal, such as information 
    about integration points and integration weights
*/

class PyraFE : public BaseFE
{
public:
   //! Constructor with type of integration rule
   PyraFE();

   //! Deconstructor
  virtual ~PyraFE(); 

  //! return FE-Type for LAS++
#ifdef USE_OLAS
  virtual FEType feType() { return PYR;}
#else
  virtual Integer feType() { return 5;}
#endif
  

protected:
  //! Set integration points
  virtual void SetIntPoints();

};

}
#endif //
