// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_CFS_NC_ELEM_HH
#define FILE_CFS_NC_ELEM_HH

#include "surfElem.hh"

namespace CoupledField
{

  //! NC_SIMON: Class for description of a finite element for coupling between
  //! non conforming interfaces

  //! This class describes no surface finite element by itself but stores
  //! pointers to two parent surface elements. One of the elements is an
  //! element on one of the non-conforming interfaces (MASTER or SLAVE). The
  //! other element is a so called "lagrange multiplier element" on which
  //! additional degrees of freedom are defined. The coupling of of the DOFs
  //! for the NCElem is inter-element whereas for Elem and SurfElem it is
  //! intra-element. 
  //!
  //! \verbatim
  //!      \               /   \           \   /           /
  //!         \         /         \         / \         /
  //!            \   /               \   /       \   /
  //!             / \                 / \         / \
  //!          /       \           /       \   /       \
  //!       /             \     /           / \           \
  //!     x-----------------x x-----------x-----x-----------x 1 2 5 9 6 10
  //!     \endverbatim
  //!     The inherited array connect specifies just node numbers
  //!     of nodes from the grid that have have been newly created by the
  //!     intersection computation on the NC-interfaces. It's meaning is
  //!     strictly geometrical since the integration involves the DOFs of the
  //!     parent elements.  The pointer ptElem stores a pointer to the BaseFE
  //!     (line, triangle, quad) corresponding to the shape of a part of or
  //!     the whole intersection of the two surface elements. It is used for
  //!     getting the location of integration points.


  struct NCElem : public SurfElem
  {
  public:
    
    //! Default Constructor
    NCElem() {
      ptLagrangeParent = NULL;
      ptSurfParent = NULL;
      isOnSlaveSide = true;
    }
    
    // ======================================================
    // GEOMETRICAL INFORMATION
    // ======================================================
    //@{ \name Geometrical Information
    
    //! NC_SIMON: Pointer to a SurfElem on which a Lagrange multiplicator is
    //! defined.
    SurfElem * ptLagrangeParent;
    
    //! NC_SIMON: Pointer to a surface element on one of MASTER or SLAVE
    //! interface.
    SurfElem * ptSurfParent;

    bool isOnSlaveSide;

    //@}
  };
  
} // end of namespace
#endif
