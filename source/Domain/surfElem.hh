// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_CFS_SURF_ELEM_HH
#define FILE_CFS_SURF_ELEM_HH

#include "elem.hh"

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
  //! point in the direction of the first volume element.
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
  //! Now \c normalDefSign points in the direction of
  //! \c surfElem.ptVolElem1 !

  struct SurfElem : public Elem
  {
  public:
    
    //! Default Constructor
    SurfElem() {
      ptVolElem1 = NULL;
      ptVolElem2 = NULL;
      normalSign = 0;
    }
    
    // ======================================================
    // GEOMETRICAL INFORMATION
    // ======================================================
    //@{ \name Geometrical Information
    
    //! Pointer to first volume element
    Elem * ptVolElem1;
    
    //! Pointer to second volume element
    Elem * ptVolElem2;
    
    //! Flag for indicating direction of surface normal
    char normalSign;
  
    //@}
  };
  
} // end of namespace
#endif
