// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_CFS_ELEM_2002
#define FILE_CFS_ELEM_2002

#include <bitset>
#include <map>

#include "General/Enum.hh"
#include "Utils/StdVector.hh"
#include "Utils/tools.hh"

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

    //! Dummy destructor
    virtual ~Elem() {;}

  public:

    // enumeration with elements types.
    // enum ElementType{Line1, Triang1, Triang2, Quadrilateral1, Quadrilateral2};

    //! Type of finite element
    //! The enumeration contains the following values:
    //! - NOFETYPE
    //! - LINE
    //! - TRIA
    //! - QUAD
    //! - TET
    //! - HEX
    //! - PYR
    //! - WED
      //  typedef enum {NOFETYPE, LINE, TRIA, QUAD, TET, HEX, PYR, WED} FEType;

    // Definition of supported element types

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

    static UInt GetNumElemNodes(FEType type);
    static UInt GetElemDim(FEType type);
    static bool GetElemQuadratic(FEType type);

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

    //! calculation of diameter of element
    Double diameter(const Point * const ptArrayOfNodes);
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

  inline Double Elem::diameter(const Point * const ptArrayOfNodes)
  {
    if (connect.GetSize()==1)
      EXCEPTION("This function is not valid for this dimension");

    Point a=ptArrayOfNodes[connect[1]];
    Point b=ptArrayOfNodes[connect[2]];
    Point c=ptArrayOfNodes[connect[3]];

    return std::max(Point::dist(a,b), Point::dist(b,c));
  }

} // end of namespace
#endif
