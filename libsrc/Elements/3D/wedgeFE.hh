#ifndef FILE_WEDGEFE
#define FILE_WEDGEFE

#include <Elements/basefe.hh>

namespace CoupledField
{

//! Class with general description of six-node wedge element
/*! This class is derived from BaseFE. 
    It stores general procedures for each type of finite element based on wedge, such as information 
    about integration points and integration weights
*/

class WedgeFE : public BaseFE
{
public:
   //! Constructor with type of integration rule
   WedgeFE();

   //! Deconstructor
  virtual ~WedgeFE(); 

  //! return FE-Type for LAS++
  virtual Integer feType() { return WED;}
  

protected:
  //! Set integration points
  virtual void SetIntPoints();

};

}
#endif //
