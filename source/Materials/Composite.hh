#ifndef FILE_COMPOSITE
#define FILE_COMPOSITE


#include "Utils/StdVector.hh"


namespace CoupledField {

  // forward class declaration
  class BaseMaterial;
  
  //! Represents composite materials (lamina) made of several layers
  struct Composite {
		
    //! Name of Composite Material
    std::string name;
		
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
