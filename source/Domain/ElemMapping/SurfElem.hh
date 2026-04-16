#ifndef FILE_CFS_SURF_ELEM_HH
#define FILE_CFS_SURF_ELEM_HH

#include <array>
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
  //! it also holds references to its neighbouring volume element(s).  //! As surface elements are almost always used for surface integration
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
      ptVolElems[0] = nullptr;
      ptVolElems[1] = nullptr;
    }
    
    // ======================================================
    // GEOMETRICAL INFORMATION
    // ======================================================
    //@{ \name Geometry Information
    
    //! Array with pointer to neighbouring volume elements
    std::array<Elem*,2> ptVolElems;
    
    //@}
  };
  
  /*! \class NcSurfElem
   *   \brief Structure specially suited for non-conforming grids and DG-methods
   *   \date 02/2012
   *   \author ahueppe
   *
   *   This struct is designed for the need of DG and mortar Element methods.
   *   One distinct feature of these schemes is that the surface elements are
   *   closely related to their volume elements thereby it is possible during creation
   *   to store directly local coordinates associated to volume elements.
   *   Furthermore we store a vector of neighbor elements which can be SurfElems or
   *   NcSurfElems again.
   */
  struct NcSurfElem : public SurfElem{
    NcSurfElem(): SurfElem(){
      //in most cases we have at least one neighbor so we reserve memory for it
      neighbors.Reserve(1);
    }

    StdVector< Vector<Double> > localCoords;

    //for compatibility with other element types...
    StdVector< shared_ptr<NcSurfElem> > neighbors;
  };

  /*! This struct was introduced to prevent interfering with the Discontinuous
   *  Galerkin implementation. It shall be merged with NcSurfElem, once
   *  non-matching grid are working again.
   */
  struct MortarNcSurfElem : public NcSurfElem {
    MortarNcSurfElem() : NcSurfElem() {}
    //! pointers to primary/secondary elements
    SurfElem *ptPrimary = nullptr;
    SurfElem *ptSecondary = nullptr;
    //! in case of non-coplanar interfaces we need projections to compute element intersections and to 
    //! transport the integration points of intersection elements back into the original elements.
    shared_ptr<SurfElem> projectedPrimary = nullptr;
    
    // Store the information about parallel projection between the primary and secondary surface in case of
    // translational p.b.c by reference
    const Vector<Double> *transVect = nullptr;
    // in case cake-piece projection (non-parallel surfaces, but aligned in a wedge or cake-piece shape)
    // also store the rotation vector and the center of rotation
    Vector<Double> rotationCenter;
    Double rotationAngle;
  };
} // end of namespace
#endif
