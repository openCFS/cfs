#ifndef FILE_HEXAFE_2003
#define FILE_HEXAFE_2003

#include <Elements/basefe.hh>

namespace CoupledField
{

//! Class with general description of hexahedral element
/*! This class is derived from BaseElem. It stores general procedures for each type of finite element based on hexahedral, such as information
  about integration points and integration weights
*/

class HexaFE : public BaseFE
{
public:
   //! Constructor with type of integration rule
   HexaFE();

   //! Deconstructor
  virtual ~HexaFE(); 

  //! return FE-Type for LAS++
#ifdef USE_OLAS
  virtual FEType feType() { return HEX;}
#else
  virtual Integer feType() { return 4;}
#endif

///////////////////////////////////////////////////////////////////////


protected:

  //! Set integration points
  virtual void SetIntPoints();

};

}
#endif //
