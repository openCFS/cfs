/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GBTYPES_HH
#define GBTYPES_HH

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include <string>

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

//! The flags for marking a triangle
typedef enum {
  NORMAL_G      =0x00000000,
  BOUNDARY_G    =0x00000001, 
  EVEN_G        =0x00000002, 
  MARK_G        =0x00000004, 
  DELETED_G     =0x00000008, 
  RED0_G        =0x00000010, 
  RED1_G        =0x00000020, 
  RED2_G        =0x00000040,
  RED_G         =0x00000070, 
  GREEN_G       =0x00000100, 
  MARK0_G       =0x00001000, 
  MARK1_G       =0x00002000, 
  MARK2_G       =0x00004000,
  MARKED_G      =0x00008000,
  SUBDIVIDED_G  =0x00010000,
  MARKED_FOR_REUSE =0x00020000,
  ALL_G         =0xFFFFFFFF
} GbGeoStatusFlag;

//! The flags for marking the vertex
typedef enum {
  NORMAL_V       =0x00000000,
  BOUNDARY_V     =0x00000001, 
  EVEN_V         =0x00000002, 
  USED_V         =0x00000004, 
  DELETED_V      =0x00000008,
  INTERIOR_V     =0x00000010,
  INFLOW_V       =0x00000020,
  OUTFLOW_V      =0x00000040,
  EULERWALL_V    =0x00000100,
  WALLISOTHERM_V =0x00000200,
  WALLADIABAT_V  =0x00000400,
  SYMMETRY_X_V   =0x00001000,
  SYMMETRY_Y_V   =0x00002000,
  SYMMETRY_Z_V   =0x00004000,
  MARKED_V       =0x00008000,
  ALL_V          =0xFFFFFFFF
} GbVertexFlag;

//! The flags for marking an edge
typedef enum {
  NORMAL_E       =0x00000000,
  BOUNDARY_E     =0x00000001, 
  EVEN_E         =0x00000002, 
  USED_E         =0x00000004, 
  DELETED_E      =0x00000008,
  INTERIOR_E     =0x00000010,
  ALL_E          =0xFFFFFFFF
} GbEdgeStatusFlag;

typedef enum { 
  UNKNOWN_OBJECTS       =0x00000000, 
  TRIANGLE_OBJECTS      =0x00000001, 
  QUADRILATERAL_OBJECTS =0x00000002,
  SURFACE_OBJECTS       =0x00000003, 
  TETRAHEDRON_OBJECTS   =0x00000010, 
  OCTAHEDRON_OBJECTS    =0x00000020, 
  PYRAMID_OBJECTS       =0x00000040, 
  PRISM_OBJECTS         =0x00000080, 
  HEXAHEDRON_OBJECTS    =0x00000100,
  VOLUME_OBJECTS        =0x000001F0,
  BOUNDARY_OBJECTS      =0x00001000,
  HEC_OBJECTS           =0x00010000,
  EC_OBJECTS            =0x00020000,
  ALL_OBJECTS           =0xFFFFFFFF
} GbContents;

typedef enum {
  PRIMAL_MESH, DUAL_MESH
} GbPartitionMethod;

typedef enum {
  PRIMAL_STYLE, CATMULL_CLARK_STYLE, SQRT3_STYLE
} GbSubdivisionStyle;

typedef enum {
  LINEAR_SCHEME, CATMULL_CLARK_SCHEME, BUTTERFLY_SCHEME, LOOP_SCHEME, 
  SQRT3_APPROXIMATING_SCHEME, SQRT3_INTERPOLATORY_SCHEME
} GbSubdivisionScheme;

typedef bool GbBool;

typedef std::string GbString;

class GbOracle 
{
public:
  virtual ~GbOracle() {}
  virtual GbBool operator()() const = 0;
  virtual GbBool operator()(int l, int id) const = 0;
};

// for use in STL map<GbKeyName,SomeType>
class GbKeyName
{
public:
  GbKeyName () { name_ = ""; }
  GbKeyName (GbString kName) { name_ = kName; }

  operator const char*() const { return name_.c_str(); }

  GbBool operator< (const GbKeyName& rkKey) const {
    return name_ < rkKey.name_;
  }
  
protected:
  GbString name_;
};

#endif // GBTYPES_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.2  2002/03/21 14:58:57  elena
| new: changes in dat-file for reading tetrahedral (bugs in element connection)
|
| Revision 1.13  2001/09/12 09:28:40  prkipfer
| introduced adaptive tet subdivision
|
| Revision 1.12  2001/06/15 08:51:06  prkipfer
| introduced hash table key type
|
| Revision 1.11  2001/05/03 08:25:37  uflabsik
| constants added
|
| Revision 1.10  2001/04/17 12:40:47  prkipfer
| added binary complement to enums
|
| Revision 1.9  2001/02/13 11:03:17  prkipfer
| introduced boundary vertices
|
| Revision 1.8  2001/01/02 15:21:33  prkipfer
| introduced cloning and support for new base classes
|
| Revision 1.7  2000/12/12 10:06:01  prkipfer
| added new element and vertex flags
|
| Revision 1.6  2000/09/07 16:57:17  prkipfer
| added subdivision
|
| Revision 1.5  2000/08/28 09:25:20  prkipfer
| added localized subdivision mechanism
|
| Revision 1.4  2000/08/22 17:15:22  prkipfer
| added possibility to query mesh contents
|
| Revision 1.3  2000/07/25 09:10:04  prkipfer
| new content enum
|
| Revision 1.2  2000/06/14 15:39:10  prkipfer
| improved base classes and added funcstruct processing for mesh
|
| Revision 1.1.1.1  2000/06/08 16:24:43  prkipfer
| Imported source tree to start the project
|
|
+---------------------------------------------------------------------*/
