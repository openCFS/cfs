#ifndef FILE_CFS_SURF_ELEM_HH
#define FILE_CFS_SURF_ELEM_HH

#include <boost/array.hpp>
#include "Elem.hh"

namespace CoupledField
{

  //! Class for description of a surface finite element

  //! This class describes a surface finite element, which means that its
  //! dimension is one smaller than the highest dimensional objects in the grid.
  //! A surface element has all the properties of a volume element (geometric 
  //! and computational ones), which means that it has an element number, a 
  //! region identifier, a reference finite element for computation and
  //! all other properties of a volume element.
  //!
  //! But since a surface element has at least one volume element neighbour,
  //! it also holds references to its neighbouring volume element(s).
  struct SurfElem : public Elem
  {
  public:
    
    //! Default Constructor
    SurfElem() {
      ptVolElems[0] = NULL;
      ptVolElems[1] = NULL;
    }
    
    // ======================================================
    // GEOMETRICAL INFORMATION
    // ======================================================
    //@{ \name Geometry Information
    
    //! Array with pointer to neighbouring volume elements
    boost::array<Elem*,2> ptVolElems;
    
    //@}
  };
  
  /*! \class NcSurfElem
   *   \brief Structure specially suited for non-conforming grids and DG-methods
   *   \date 02/2012
   *   \author ahueppe
   *
   *   This struct is desigend for the need of DG and mortar Element methods.
   *   One destinct feature of these schemes is that the surface elements are
   *   closely related to their volume elements thereby it is possible during creation
   *   to store driectly local coordinates associated to volume elements.
   *   Furthermore we store a vector of neighbor elements which can be surfelems or
   *   NcSurfElems again.
   */
  struct NcSurfElem : public SurfElem{
    NcSurfElem(): SurfElem(){
      //in most cases we have at least one neighbor so we reserve memory for it
      neighbors.Reserve(1);
    }

    StdVector< Vector<Double> > localCoords;

    //waring this may not be freed!!!! just for comptability with other element types...
    StdVector< shared_ptr<NcSurfElem> > neighbors;
  };

} // end of namespace
#endif
