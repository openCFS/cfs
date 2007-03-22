// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_COMPOSITE
#define FILE_COMPOSITE


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
