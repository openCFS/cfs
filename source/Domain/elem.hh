// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_CFS_ELEM_2002
#define FILE_CFS_ELEM_2002

#include <stddef.h>
#include <bitset>
#include <map>
#include <ostream>
#include <string>
#include <utility>

#include "General/Enum.hh"
#include "General/defs.hh"
#include "General/environment.hh"
#include "Utils/Point.hh"
#include "Utils/StdVector.hh"

namespace CoupledField
{
  class BaseFE;

  //! Class for description of a volume finite element

  //! This class describes a volume finite element, where volume means the
  //! highest dimensional element entities in the current mesh.
  //! It has to be very lightweight, since this object is created many times.
  //! It relates the geometric information of an element (node numbers)
  //! with the mathematical / computational one (reference finite element).
  //! The finite element is described by:
  //! - corner node numbers
  //! - pointer to reference finite element
  //! - element number
  //! - element subdomain identifier
  //! - refinement flag / number

  struct Elem
  {
  public:

    //! Dummy constructor
    Elem();

    virtual ~Elem()
    {
      if(neighborhood) { delete neighborhood; neighborhood = NULL; }
    }


    /** Definition of supported element types */
    typedef enum
    {
      UNDEF      = 0,
      POINT      = 1,
      LINE2      = 2,
      LINE3      = 3,
      TRIA3      = 4,
      TRIA6      = 5,
      QUAD4      = 6,
      QUAD8      = 7,
      QUAD9      = 8,
      TET4       = 9,
      TET10      = 10,
      HEXA8      = 11,
      HEXA20     = 12,
      HEXA27     = 13,
      PYRA5      = 14,
      PYRA13     = 15,
      WEDGE6     = 16,
      WEDGE15    = 17
    } FEType;
    static Enum<FEType> feType;

    /** Can this feType be a surface element */
    static bool IsSurfaceElement(FEType type, int dim);

    static UInt GetNumElemNodes(FEType type);
    static UInt GetElemDim(FEType type);
    static bool GetElemQuadratic(FEType type);

    /** Expensive method - meant for asserts but not real world applications!
     * Better use the barycenter attribute!! */
    Point ExpensiveCalcBarycenter(bool updated = false) const;

    // ======================================================
    // GEOMETRICAL INFORMATION
    // ======================================================

    //@{ \name Geometrical Information
    //! global element number
    UInt elemNum;

    //! identifier for region
    RegionIdType regionId;

    //! array with node numbers
    StdVector<UInt> connect;

    //! array with edge numbers
    StdVector<Integer> edges;

    //! array with face numbers
    StdVector<Integer> faces;

    //! bitset describing the orientation of the faces (3 for each)
    StdVector<std::bitset<3> > faceFlags;

    /** This defines the neighborhood of this element.
     * The pair entries are: first = neighbor element and second = number
     * of common nodes with this element. By this one can determine
     * if it is an face, edge or node neighbor.
     * The list is completely unsorted. To be generated via grid.
     * @see Grid::FindElementNeighorhood() */
    StdVector<std::pair<Elem*, int> >* neighborhood;

    /** The barycenter of the element, Set via Grid::SetElementBarycenters().
     * The values are by for the uninitialized case zero, be careful! Check via Grid::RegionData */
    Point barycenter;

#ifdef ADAPTGRID
    //! flag for refinement
    bool refinementFlag;

    //! number of refinement for the element
    UInt refinementNumber;

#endif
    //@}

    // ======================================================
    // COMPUTATIONAl INFORMATION
    // ======================================================

    //@{ \name Computational Information

    //! pointer to reference element representation
    BaseFE * ptElem;
    //@}

    // ======================================================
    // HELPER METHODS
    // ======================================================
    //@{ \name Helper Methods

    //! overloading operator =
    Elem & operator=(const Elem& t);
    //@}

    // Fix problems due to negative Jacobian determinants
    void CorrectConnectivity();

    std::string ToString() const
    {
      std::ostringstream os;
      os << "elemNum=" << elemNum << " region=" << regionId;
      return os.str();
    }

  private:
    static std::map<FEType, UInt> numElemNodes_;
    static std::map<FEType, UInt> elemDims_;
    static std::map<FEType, UInt> elemQuadratic_;
    static void Initialize();
  };


  inline Elem & Elem::operator=(const Elem& t)
  {
    if (this!=&t) {
      ptElem=t.ptElem;
      connect=t.connect;
      regionId=t.regionId;
#ifdef ADAPTGRID
      refinementFlag=t.refinementFlag;
      refinementNumber=t.refinementNumber;
#endif
    }
    return *this;
  }

  /** Operator to print a StdVector<Elem*> via ToString() */
  std::ostream & operator<<(std::ostream &out, const Elem*& data);

  std::ostream & operator<<(std::ostream &out, const StdVector<Elem*>& data);

  /** Operator such that we can print neighbourhood via StdVector::ToString() */
  std::ostream & operator<<(std::ostream &out, const std::pair<Elem*, int>& data);





} // end of namespace
#endif
