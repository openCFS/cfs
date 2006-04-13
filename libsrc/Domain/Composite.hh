#ifndef FILE_LAMINA
#define FILE_LAMINA


#include "General/environment.hh"


namespace CoupledField {

  // forward class declaration
  class BaseMaterial;
  
  //! Represents composite materials (lamina) made of several layers
  struct Composite {
		
    //! Name of Composite Material
    std::string name;
		
    //! Relative start height of the first lamina layer
    Double zStart;
		
    //----Definition of Layers----
		
    //! Thickness of individual lamina
    StdVector<Double> thickness;
		
    //! Material of individual lamina
    StdVector<BaseMaterial*> materials;
		
    //! Orientation with respect to lamina
    StdVector<Double> orientation;
		
  };

}// end of namespace
#endif
