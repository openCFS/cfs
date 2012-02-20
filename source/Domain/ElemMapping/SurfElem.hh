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
  //! As surface elements are almost always used for surface integration 
  //! purpose, a normal vector with defined sign is needed. Therefore this 
  //! class has also a flag, indicating the direction of the normal.
  //! This can be understand as follows: 
  //! A vector perpendicular to the element's surface can simply be found by
  //! calculating the (degenerated) cross-procduct of two vectors lying in 
  //! the surface ( = connecting corner points). However, the sign of the 
  //! resulting vector will be dependend on the ordering of the connectivity
  //! of the element. Therefore the flag \a normalSign defines a factor
  //! (either 1 or -1) by which the resulting vector has to be multiplied, to
  //! point in the direction OUT of the first volume element.
  //! 
  //! Here is a snipplet of example code to demonstrate the calculation of 
  //! the surface normal:
  //! 
  //! First, get the surface normal of the element without defined sign:
  //! \verbatim
  //! Vector<Double> normalUndefSign, normalDefSign;
  //! SurfElem surfElem;
  //! Grid * ptGrid;
  //!
  //! ptGrid->CalcSurfNormal(normalUndefSign, surfElem);
  //! \endverbatim
  //! Now, in order to get a defined direction of the normal,
  //! multiply the normal with the \a normalSign:
  //! \verbatim
  //! normalDefSign = normalUndefSign * surfElem.normalSign;
  //! \endverbatim
  //! Now \c normalDefSign points in the direction OUT of 
  //! \c surfElem.ptVolElem[0] !

  struct SurfElem : public Elem
  {
  public:
    
    //! Default Constructor
    SurfElem() {
      ptVolElems[0] = NULL;
      ptVolElems[0] = NULL;
      normalSign = 0;
    }
    
    // ======================================================
    // GEOMETRICAL INFORMATION
    // ======================================================
    //@{ \name Geometry Information
    
    //! Array with pointer to neighbouring volume elements
    boost::array<Elem*,2> ptVolElems;
    
    //! Flag for indicating direction of surface normal
    
    //! If this flag (interpreted as a double) gets multiplied with the
    //! undefined normal direction, obtained by Grid::CalcSurfNormal(), the
    //! resulting normal will be oriented OUT of ptVolElem[0], i.e. for
    //! surface element on the boundary (which have just one volume neighbor)
    //! the resulting normal will point out of the domain. 
    char normalSign;
  
    //@}
  };
  
} // end of namespace
#endif
