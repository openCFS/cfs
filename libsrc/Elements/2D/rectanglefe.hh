#ifndef FILE_RECTANGLEFE_2003
#define FILE_RECTANGLEFE_2003

#include <Elements/basefe.hh>

namespace CoupledField
{
//! Class with general procedures for quadrilateral finite elements
/*! This class is derived from BaseFE. It stores general procedures for each type of finite element on quadralateral, such as calculation of Jacobian of transformation in standart and information about integration points and integration weights */
  
class RectangleFE : public BaseFE
{
public:

  //! Constructor with type of integration rule
  RectangleFE(); 
  
  //! Deconstructor
  virtual ~RectangleFE();

  //! return FE-Type for CLA++
  virtual Integer feType() { return QUAD;}

protected:


   //! Set integration points
  virtual void SetIntPoints();

};

} // end of namespace

#endif // FILE_RECTANGLEFE_2003
