#ifndef FILE_TRIANGLEFE_2003
#define FILE_TRIANGLEFE_2003

#include <Elements/basefe.hh>

namespace CoupledField
{
//! Class with general procedures for triangle finite elements
/*! This class is derived from BaseFE. It stores general procedures for each type of finite element on quadralateral, such as calculation of Jacobian of transformation in standart and information about integration points and integration weights */
  
class TriangleFE : public BaseFE
{
public:

  //! Constructor with type of integration rule
  TriangleFE(); 
  
  //! Deconstructor
  virtual ~TriangleFE();

  //! return FE-Type for CLA++
  virtual Integer feType() { return TRIA;}


protected:


   //! Set integration points
  virtual void SetIntPoints();


  
};

} // end of namespace

#endif // FILE_RECTANGLEFE_2003
