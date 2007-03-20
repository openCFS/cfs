#ifndef XMDF_DOT_H
#define XMDF_DOT_H

/*
XMDF License & Copyright Notice
Copyright Notice and Statement for the eXtensible Model Data Format (XMDF) Software Library and API.
XMDF Software Library and Utilities Copyright 2004, 2005 by the Enivironmental Modeling Research Laboratory and Brigham Young University.  All rights reserved.
Redistribution and use in source and binary forms, with or without modification, are permitted for any purpose (including commercial purposes) provided that the following conditions are met:
1.  Redistributions of source code must retain the above copyright notice, this list of conditions, and the following disclaimer.
2.  Redistributions in binary form must reproduce the above copyright notice, this list of conditions, and the following disclaimer in the documentation and/or materials provided with the distribution.
3.  In addition, redistributions of modified forms of the source or binary code must carry prominent notices stating that the original code was changed and the date of the change.
4.  All publications or advertising materials mentioning features or use of this software are asked, but not required, to acknowledge that it was developed by the EMRL at Brigham Young University in partnership with the US Army Engineer Research and Development Center and to credit the contributors.
5.  Neither the name of the University nor the names of the Contributors maybe used to endorse or promote products derived from this software without specific prior written permission from the University or the Contributors, as appropriate for the name(s) to be used.
6.  THIS SOFTWARE IS PROVIDED BY THE UNIVERSITY AND THE CONTRIBUTORS "AS IS"WITH NO WARRANTY OF ANY KIND, EITHER EXPRESSED OR IMPLIED.  In no event shall the University or the Contributors be liable for any damages suffered by the users arising out of the use of this software, even if advised of the possibility of such damage.
--------------------------------------------------------------------------
Portions of XMDF utilize the HDF5 software library.  The following statement applies to HDF5:
Copyright Notice and Statement for NCSA Hierarchical Data Format (HDF) Software Library and Utilities
NCSA HDF5 (Hierarchical Data Format 5) Software Library and Utilities Copyright 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005 by the Board of Trustees of the University of Illinois.  All rights reserved.
Contributors: National Center for Supercomputing Applications (NCSA) at the University of Illinois at Urbana-Champaign (UIUC), Lawrence Livermore National Laboratory (LLNL), Sandia National Laboratories (SNL), Los Alamos National Laboratory (LANL), Jean-loup Gailly and Mark Adler (gzip library).
Redistribution and use in source and binary forms, with or without modification, are permitted for any purpose (including commercial purposes) provided that the following conditions are met:
1.  Redistributions of source code must retain the above copyright notice, this list of conditions, and the following disclaimer.
2.  Redistributions in binary form must reproduce the above copyright notice, this list of conditions, and the following disclaimer in the documentation and/or materials provided with the distribution.
3.  In addition, redistributions of modified forms of the source or binary code must carry prominent notices stating that the original code was changed and the date of the change.
4.  All publications or advertising materials mentioning features or use of this software are asked, but not required, to acknowledge that it was developed by the National Center for Supercomputing Applications at the University of Illinois at Urbana-Champaign and to credit the contributors.
5.  Neither the name of the University nor the names of the Contributors may be used to endorse or promote products derived from this software without specific prior written permission from the University or the Contributors, as appropriate for the name(s) to be used.
6.  THIS SOFTWARE IS PROVIDED BY THE UNIVERSITY AND THE CONTRIBUTORS "AS IS" WITH NO WARRANTY OF ANY KIND, EITHER EXPRESSED OR IMPLIED.  In no event shall the University or the Contributors be liable for any damages suffered by the users arising out of the use of this software, even if advised of the possibility of such damage.
--------------------------------------------------------------------------
Portions of HDF5 were developed with support from the University of California, Lawrence Livermore National Laboratory (UC LLNL).
The following statement applies to those portions of the product and must be retained in any redistribution of source code, binaries, documentation, and/or accompanying materials:
This work was partially produced at the University of California, Lawrence Livermore National Laboratory (UC LLNL) under contract no. W-7405-ENG-48 (Contract 48) between the U.S. Department of Energy (DOE) and The Regents of the University of California (University) for the operation of UC LLNL.
DISCLAIMER:
This work was prepared as an account of work sponsored by an agency of the United States Government.  Neither the United States Government nor the University of California nor any of their employees, makes any warranty, express or implied, or assumes any liability or responsibility for the accuracy, completeness, or usefulness of any information, apparatus, product, or process  disclosed, or represents that its use would not infringe privately- owned rights.  Reference herein to any specific commercial products,  process, or service by trade name, trademark, manufacturer, or  otherwise, does not necessarily constitute or imply its endorsement,  recommendation, or favoring by the United States Government or the  University of California.  The views and opinions of authors expressed herein do not necessarily state or reflect those of the United States Government or the University of California, and shall not be used for advertising or product endorsement purposes.
--------------------------------------------------------------------------
This work was partially produced at the University of California, Lawrence Livermore National Laboratory (UC LLNL) under contract no.W-7405-ENG-48 (Contract 48) between the U.S. Department of Energy (DOE) and The Regents of the University of California (University) for the operation of UC LLNL.
DISCLAIMER:This work was prepared as an account of work sponsored by an agency of the United States Government.  Neither the United States Government nor the University of California nor any of their employees, makes any warranty, express or implied, or assumes any liability or responsibility for the accuracy, completeness, or usefulness of any information, apparatus, product, or process disclosed, or represents that its use would not infringe privately-owned rights.  Reference herein to any specific commercial products, process, or service by trade name, trademark, manufacturer, or otherwise, does not necessarily constitute or imply its endorsement, recommendation, or favoring by the United States Government or the University of California.  The views and opinions of authors expressed herein do not necessarily state or reflect those of the United States Government or the University of California, and shall not be used for advertising or product endorsement purposes.
--------------------------------------------------------------------------
*/

#include <def_use_hdf5.hh>

/* The following ifdef block is the standard way of creating macros which make
   exporting from a DLL simpler. All files within this DLL are compiled with
   the XMDF_EXPORTS symbol defined on the command line. this symbol should
   not be defined on any project that uses this DLL. This way any other project
   whose source files include this file see  functions as being
   imported from a DLL, wheras this DLL sees symbols defined with this macro
   as being exported. */
#ifdef NO_XMDF_API
#define XMDF_API  
#else
#ifdef XMDF_EXPORTS
#define XMDF_API __declspec(dllexport) 
#else
#define XMDF_API __declspec(dllimport)
#endif
#endif

/* XMDF version number */
/* Single precision float */ 
#define XMDF_VERSION  1.20


#include <hdf5.h>
#include <float.h>
#include <limits.h>
#include <math.h>

/* type defines */
typedef unsigned char xmbool;
typedef hid_t xid;

#define XFALSE 0
#define XTRUE  !XFALSE

#define MAX_ID_NAME 512

typedef struct {
  hid_t       id;
  H5I_type_t  type;
  char        name[MAX_ID_NAME];
} IdentifierInfo;


/* #defines */
#define NONE -1

#ifndef MAX_FLOAT
#define MAX_FLOAT FLT_MAX
#endif

#ifndef MIN_FLOAT
#define MIN_FLOAT (-FLT_MAX)
#endif

  /* File Version */
#define FILE_VERSION "File Version"

  /* Property Types */
#define XF_TYPE_INT     1
#define XF_TYPE_FLOAT   2
#define XF_TYPE_DOUBLE  3
#define XF_TYPE_STRING  4
#define XF_TYPE_UINT    5
#define XF_TYPE_OTHER   11

  /* Time Variables */
#define TIME_UNITS "TimeUnits"
#define TS_SECONDS "Seconds"
#define TS_MINUTES "Minutes"
#define TS_HOURS   "Hours"
#define TS_DAYS    "Days"
#define TS_NOT_APPLICABLE "None"
#define TIME_UNITS_MAXLENGTH 25

#define NO_UNITS "None"
#define UNITS_MAXLENGTH 100

/* Element types */
#define EL_TYPE_JUNCTION3              3
#define EL_TYPE_JUNCTION4              4
#define EL_TYPE_JUNCTION5              5
#define EL_TYPE_JUNCTION6              6
#define EL_TYPE_JUNCTION7              7
#define EL_TYPE_JUNCTION8              8
  /* 1D elements  */
#define EL_TYPE_1D_LINEAR            100 
#define EL_TYPE_1D_QUAD              101 
#define EL_TYPE_TRANS_1D_2D          110
  /* 2D elements  */
#define EL_TYPE_TRI_LINEAR           200
#define EL_TYPE_TRI_QUAD             201
#define EL_TYPE_QUADRILATERAL_LINEAR 210
#define EL_TYPE_QUADRILATERAL_QUAD8  211
#define EL_TYPE_QUADRILATERAL_QUAD9  212
  /* 3D elements */
#define EL_TYPE_TET_LINEAR           300
#define EL_TYPE_PRISM_LINEAR         310
#define EL_TYPE_HEX_LINEAR           320
#define EL_TYPE_PYRAMID_LINEAR       330

  /* Grid types */
#define GRID_TYPE_CARTESIAN            0
#define GRID_TYPE_CURVILINEAR          1
#define GRID_TYPE_CARTESIAN_EXTRUDED   2
#define GRID_TYPE_CURVILINEAR_EXTRUDED 3
#define GRID_TYPE_MIN                  0
#define GRID_TYPE_MAX                  3

  /* extrusion types for grids  */
#define EXTRUDE_SIGMA                  0
#define EXTRUDE_CARTESIAN              1
#define EXTRUDE_CURV_AT_CORNERS        2
#define EXTRUDE_CURV_AT_CELLS          3
#define EXTRUDE_MIN                    0
#define EXTRUDE_MAX                    3

 /* Grid orientation  */
#define ORIENTATION_RIGHT_HAND         0
#define ORIENTATION_LEFT_HAND          1

 /* Reserved property names */
#define PROP_ACTIVITY                 "activity"
#define PROP_NULL_VALUE               "nullvalue"

/* Dataset value locations for a grid */
#define GRID_LOC_CENTER                0
#define GRID_LOC_CORNER                1
#define GRID_LOC_FACES                 2
#define GRID_LOC_FACE_I                3
#define GRID_LOC_FACE_J                4
#define GRID_LOC_FACE_K                5
#define GRID_LOC_COLUMN                6
#define GRID_LOC_CENTER_CORNERS_FACES  7 /* Centers, corners, and faces */

/* #defines for coordinate conversions  */
  /* Horizontal datum  */
#define HORIZ_DATUM_LOCAL               0
#define HORIZ_DATUM_GEOGRAPHIC          1
#define HORIZ_DATUM_GEOGRAPHIC_NAD27    2
#define HORIZ_DATUM_GEOGRAPHIC_NAD83    3
#define HORIZ_DATUM_GEOGRAPHIC_HPGN     4
#define HORIZ_DATUM_UTM                 5
#define HORIZ_DATUM_UTM_NAD27           6
#define HORIZ_DATUM_UTM_NAD83           7
#define HORIZ_DATUM_UTM_HPGN            8
#define HORIZ_DATUM_STATE_PLANE_NAD27   9
#define HORIZ_DATUM_STATE_PLANE_NAD83  10
#define HORIZ_DATUM_STATE_PLANE_HPGN   11
#define HORIZ_DATUM_CPP                12

  /* Horizontal or vertical units  */
#define COORD_UNITS_US_FEET             0
#define COORD_UNITS_INTERNATIONAL_FEET  1
#define COORD_UNITS_METERS              2

  /* vertical datum */
#define VERT_DATUM_LOCAL                0
#define VERT_DATUM_NGVD_29              1
#define VERT_DATUM_NGVD_88              2

  /* latitude  */
#define LATITUDE_NORTH                  0
#define LATITUDE_SOUTH                  1

  /* longitude */
#define LONGITUDE_EAST                  0
#define LONGITUDE_WEST                  1

  /* valid UTM zones */
#define UTM_ZONE_MIN                    1
#define UTM_ZONE_MAX                   60

  /* Guid names and lengths */
#define XF_GUID                    "Guid"
#define XF_GUID_STRINGLENGTH           37

  /* Material Mapping File for Classifying */
#define XMDF_MAT_MAPING_TARGET_IDS     "Target Ids"
#define XMDF_MAT_MAPING_BACKGROUND_IDS "Background Ids"
#define XMDF_MAT_MAPING_IDS            "Mapping Ids"

/* Calendar Eras */
#define ERA_IS_BCE   0
#define ERA_IS_CE    1

/* Macros */
#define Xmax(a,b)    (((a) >= (b)) ? (a) : (b))
#define Xmin(a,b)    (((a) >= (b)) ? (b) : (a))

#define MAX_DOUBLE 1.7976931348623157e+308 
#define MIN_DOUBLE -1.7976931348623157e+308

/* Property for node ids */
#define XF_IDS                      "Ids"

/* Error Retrieval */
#define XF_MAX_ERROR_MSG_SIZE         256
#define XF_MAX_ERROR_STACK            256

/* Overwrite options in the function xfSetupToWriteDatasets */
#define XF_OVERWRITE_CLEAR_FILE           1
#define XF_OVERWRITE_CLEAR_DATASET_GROUP  2
#define XF_OVERWRITE_NONE                 3

/* File Functions */
#ifdef __cplusplus
extern "C" {
#endif

/* Library versions */
XMDF_API int xfGetLibraryVersion(float *a_Version);
XMDF_API int xfGetLibraryVersionFile(xid a_File, float *a_Version);
/* Coordinates */
XMDF_API int xfCreateCoordinateGroup (xid a_ParentId, xid *a_ChildId);
XMDF_API int xfOpenCoordinateGroup (xid a_ParentId, xid *a_ChildId);

XMDF_API int xfCreateFile(const char *a_File, xid *Id, xmbool a_Overwrite);
XMDF_API int xfOpenFile(const char *a_File, xid *Id, xmbool a_ReadOnly);
XMDF_API int xfCloseFile(xid a_FileId);

/* retreiving error messages */
/* error messages are cleared whenever a file is opened or closed */
XMDF_API int xfGetNumErrorMessages(int *a_Num);
/* a_Errors must be allocated to a_Num * XF_MAX_ERROR_MSG_SIZE */
XMDF_API int xfGetErrorMessages(int a_Num, char **a_Errors);

/* functions to determine and close identifiers not yet closed */
XMDF_API int xfGetNumOpenIdentifiers(xid a_FileId, int *a_Num);
XMDF_API int xfGetOpenIdentifiersInfo(xid a_FileId, int a_Num,
                                      IdentifierInfo *a_Info);

/* Functions used to work with Properties */
XMDF_API int xfCreatePropertyGroup(xid a_ParentId, xid *Id);

XMDF_API int xfWritePropertyString(xid a_Id, const char *a_Name, 
                   int a_Number, int a_StringLength, const char *a_Properties);
XMDF_API int xfWritePropertyUnsignedInt(xid a_Id, const char *a_Name, 
                   int a_Number, const unsigned int *a_Properties, 
                   int a_Compression);
XMDF_API int xfWritePropertyInt(xid a_Id, const char *a_Name, 
                   int a_Number, const int *a_Properties,
                   int a_Compression);
XMDF_API int xfWritePropertyDouble(xid a_Id, const char *a_Name, 
                   int a_Number, const double *a_Properties,
                   int a_Compression);
XMDF_API int xfWritePropertyFloat(xid a_Id, const char *a_Name, 
                   int a_Number, const float *a_Properties,
                   int a_Compression);

XMDF_API int xfDoesPropertyWithNameExist(xid a_Id, const char *a_Name,
                                             xmbool *Exists);
XMDF_API int xfGetPropertyNumber(xid a_Id, const char *a_Name,
                                    int *a_Number);
XMDF_API int xfGetPropertyStringLength(xid a_Id, const char *a_Name, 
                                   int *a_Number, int *a_MaxLength);
XMDF_API int xfGetPropertyType(xid a_GroupId, const char *a_Name, int *a_Type);
XMDF_API int xfAllocateReadPropertyString (xid a_Id, const char *a_Name,
                                           int *a_Number, int *a_MaxLength,
                                           char **a_Properties);
XMDF_API int xfReadPropertyString(xid a_Id, const char *a_Name,
                    int a_Number, int a_MaxLength, char *a_Properties);
XMDF_API int xfReadPropertyUnsignedInt(xid a_Id, const char *a_Name,
                    int a_Number, unsigned int *a_Properties);
XMDF_API int xfReadPropertyInt(xid a_Id, const char *a_Name,
                    int a_Number, int *a_Properties);
XMDF_API int xfReadPropertyDouble(xid a_Id, const char *a_Name,
                    int a_Number, double *a_Properties);
XMDF_API int xfReadPropertyFloat(xid a_Id, const char *a_Name,
                    int a_Number, float *a_Properties);

  /* group functions */
XMDF_API int xfCreateGenericGroup(xid a_FileId, const char *a_Path,
                                     xid *a_GroupId);
XMDF_API int xfCreateGroupForMesh(xid a_FileId, const char *a_Path,
                                         xid *a_GroupId);
XMDF_API int xfCreateGroupForGrid(xid a_FileId, const char *a_Path,
                                         xid *a_GroupId);
XMDF_API int xfCreateStationForGrid(xid a_FileId, const char *a_Path,
                                         xid *a_GroupId);
XMDF_API int xfCreateGroupForXsec(xid a_FileId, const char *a_Path,
                                         xid *a_GroupId);
XMDF_API int xfCreateGroupForMatSim(xid a_FileId, const char *a_Path,
                                         xid *a_GroupId);

XMDF_API int xfOpenGroup (xid a_ParentId, const char *a_Path, 
                             xid *a_GroupId);
XMDF_API int xfCloseGroup(xid GroupId);

XMDF_API int xfGetGroupPathsSizeForMeshes(xid FileId, int *Num,int *Maxsize);
XMDF_API int xfGetGroupPathsForMeshes(xid a_FileId, int a_Num, int a_Maxsize,
                                         char *a_Paths);

XMDF_API int xfGetGroupPathsSizeForGrids(xid FileId, int *Num,int *Maxsize);
XMDF_API int xfGetGroupPathsForGrids(xid a_FileId, int a_Num, int a_Maxsize,
                                     char *a_Paths);

XMDF_API int xfGetGroupPathsSizeForXsecs(xid a_FileId, int *a_Num,int *Maxsize);
XMDF_API int xfGetGroupPathsForXsecs(xid a_FileId, int a_Num, int a_Maxsize,
                                     char *a_Paths);

XMDF_API int xfGetGroupPathsSizeForGeomPaths(xid FileId, int *Num,int *Maxsize);
XMDF_API int xfGetGroupPathsForGeomPaths(xid a_FileId, int a_Num, int a_Maxsize,
                                         char *a_Paths);

XMDF_API int xfOpenPropertyGroup(xid a_ParentId, xid *a_GroupId);

XMDF_API int xfGetGroupAbsolutePathSize(xid a_GroupId, int *a_PathLength);
XMDF_API int xfGetGroupAbsolutePath(xid a_GroupId, int a_PathLength,
                                          char *a_Path);

  // float and double write types
XMDF_API hid_t xfGetDoubleType (void);
XMDF_API hid_t xfGetFloatType (void);
XMDF_API void xfSetFloatType (int a_BigEndian);

/* writing/reading reftimes */
XMDF_API int xfWriteReftime(xid a_Id, double a_Reftime);
XMDF_API int xfUseReftime(xid a_Id, xmbool *a_bUseReftime);
XMDF_API int xfReadReftime(xid a_Id, double *a_Reftime);

  /* functions to convert between calendar dates and julian dates */
XMDF_API int xfCalendarToJulian(xmbool a_bEra, int a_yr, int a_mo, int a_day, 
                        int a_hr, int a_min, int a_sec, double *a_julian);
XMDF_API int xfJulianToCalendar (xmbool *a_bEra, int *a_yr, int *a_mo, int *a_day,
                        int *a_hr, int *a_min, int *a_sec, double a_julian);

  /* mesh functions */
    /* writing */
XMDF_API int xfSetNumberOfElements(xid a_Id, int a_nElems);
XMDF_API int xfSetAllElemsSameType(xid a_Id, int a_Type);
XMDF_API int xfWriteElemTypes(xid a_Id, int a_nElems, const int *a_Type,
                                 int Compression);

XMDF_API int xfWriteElemNodeIds(xid a_Id, int a_nElems, int a_nMaxNodes,
                                   int *a_Ids, int Compression);
XMDF_API int xfSetNumberOfNodes(xid a_Id, int a_nNodes);
  /* x locations must be written first */
XMDF_API int xfWriteXNodeLocations(xid a_Id, int a_nNodes, double *a_Locs,
                                      int Compression);
XMDF_API int xfWriteYNodeLocations(xid a_Id, int a_nNodes, double *a_Locs);
XMDF_API int xfWriteZNodeLocations(xid a_Id, int a_nNodes, double *a_Locs);

  /* reading */
XMDF_API int xfGetNumberOfElements(xid a_Id, int *a_nElems);
XMDF_API int xfAreAllElemsSameType(xid a_Id, xmbool *a_Same);
XMDF_API int xfReadElemTypesSingleValue(xid a_Id, int *a_Type);
XMDF_API int xfReadElemTypes(xid a_Id, int a_nElems, int *a_Type);

XMDF_API int xfGetMaxNodesInElem(xid a_Id, int *a_nMaxNodes);
XMDF_API int xfReadElemNodeIds(xid a_Id, int a_nElems, int a_nMaxNodes,
                                  int *a_Ids);

XMDF_API int xfGetNumberOfNodes(xid a_Id, int *a_nNodes);
XMDF_API int xfReadXNodeLocations(xid a_Id, int a_nNodes, double *a_Locs);
XMDF_API int xfReadYNodeLocations(xid a_Id, int a_nNodes, double *a_Locs);
XMDF_API int xfReadZNodeLocations(xid a_Id, int a_nNodes, double *a_Locs);

  /* Mesh property groups */
XMDF_API int xfCreateMeshPropertyGroup (xid a_Id, xid *a_PropId);
XMDF_API int xfGetMeshPropertyGroup(xid a_Id, xid *a_PropId);
XMDF_API int xfCreateMeshNodePropertyGroup (xid a_Id, xid *a_PropId);
XMDF_API int xfGetMeshNodePropertyGroup(xid a_Id, xid *a_PropId);
XMDF_API int xfCreateMeshElementPropertyGroup (xid a_Id, xid *a_PropId);
XMDF_API int xfGetMeshElementPropertyGroup(xid a_Id, xid *a_PropId);

  /* Grid functions */
XMDF_API int xfSetGridType(xid a_Id, int a_GridType);

XMDF_API int xfSetNumberOfDimensions(xid a_Id, int a_NumDimensions);
XMDF_API int xfSetExtrusionType(xid a_Id, int a_ExtrudeType);
XMDF_API int xfSetOrigin(xid a_Id, double a_x, double a_y, double a_z);
XMDF_API int xfSetOrientation(xid a_Id, int a_Orientation);
XMDF_API int xfSetBearing(xid a_Id, double a_Bearing);
XMDF_API int xfSetDip(xid a_Id, double a_Dip);
XMDF_API int xfSetRoll(xid a_Id, double a_Roll);
XMDF_API int xfSetComputationalOrigin(xid a_Id, int a_origin);
XMDF_API int xfSetUDirection(xid a_Id, int a_direction);

XMDF_API int xfSetNumberCellsInI(xid a_Id, int a_NumI);
XMDF_API int xfSetNumberCellsInJ(xid a_Id, int a_NumJ);
XMDF_API int xfSetNumberCellsInK(xid a_Id, int a_NumK);

XMDF_API int xfSetGridCoordsI(xid a_Id, int a_NumVals, double *a_iValues);
XMDF_API int xfSetGridCoordsJ(xid a_Id, int a_NumVals, double *a_jValues);
XMDF_API int xfSetGridCoordsK(xid a_Id, int a_NumVals, double *a_kValues);

XMDF_API int xfWriteExtrudeLayerData(xid a_Id, int a_NumLayers, int a_NumVals,
                                     double *a_Values);

XMDF_API int xfGetGridType(xid a_Id, int *a_GridType);
XMDF_API int xfGetExtrusionType(xid a_Id, int *a_ExtrudeType);
XMDF_API int xfGetNumberOfDimensions(xid a_Id, int *a_NumDimensions);
XMDF_API int xfOriginDefined(xid a_Id, xmbool *a_bDefined);
XMDF_API int xfGetOrigin(xid a_Id, double *a_x, double *a_y, double *a_z);
XMDF_API int xfGetOrientation(xid a_Id, int *a_Orientation);
XMDF_API int xfBearingDefined(xid a_Id, xmbool *a_bDefined);
XMDF_API int xfGetBearing(xid a_Id, double *a_bearing);
XMDF_API int xfDipDefined(xid a_Id, xmbool *a_bDefined);
XMDF_API int xfGetDip(xid a_Id, double *a_dip);
XMDF_API int xfRollDefined(xid a_Id, xmbool *a_bDefined);
XMDF_API int xfGetRoll(xid a_Id, double *a_Roll);
XMDF_API int xfComputationalOriginDefined(xid GroupId, xmbool *bDefined);
XMDF_API int xfGetComputationalOrigin(xid GroupId, int *Origin);
XMDF_API int xfGetUDirectionDefined(xid GroupId, xmbool *bDefined);
XMDF_API int xfGetUDirection (xid GroupId, int *Direction);
XMDF_API int xfGetNumberCellsInI(xid a_Id, int *a_NumI);
XMDF_API int xfGetNumberCellsInJ(xid a_Id, int *a_NumJ);
XMDF_API int xfGetNumberCellsInK(xid a_Id, int *a_NumK);

XMDF_API int xfGetGridCoordsI(xid a_Id, int a_NumVals, double *a_iValues);
XMDF_API int xfGetGridCoordsJ(xid a_Id, int a_NumVals, double *a_jValues);
XMDF_API int xfGetGridCoordsK(xid a_Id, int a_NumVals, double *a_kValues);

XMDF_API int xfGetExtrudeNumLayers(xid a_Id, int *a_NumLayers);
XMDF_API int xfGetExtrudeValues(xid a_Id, int a_NumVals, double *a_Values);

    /* Grid Property group */
XMDF_API int xfCreateGridPropertyGroup(xid a_Id, xid *a_PropId);
XMDF_API int xfGetGridPropertyGroup(xid a_Id, xid *a_PropId);
XMDF_API int xfCreateGridCellPropertyGroup(xid a_Id, xid *a_PropId);
XMDF_API int xfGetGridCellPropertyGroup(xid a_Id, xid *a_PropId);
XMDF_API int xfCreateGridNodePropertyGroup(xid a_Id, xid *a_PropId);
XMDF_API int xfGetGridNodePropertyGroup(xid a_Id, xid *a_PropId);






  /* Xsec Functions */
//wmsejj

XMDF_API int xfSetNumberOfXSects (xid a_Id, int *a_nXSects, 
                                  int a_compression);
XMDF_API int xfGetNumberOfXSects (xid a_Id, int *a_nXSects);
XMDF_API int xfSetCSID (xid a_Id, int a_NumVals, int *a_PropId, 
                        int a_compression);
XMDF_API int xfGetCSID (xid a_Id,  int *a_NumVals, int *a_PropId);
XMDF_API int xfSetCSName (xid a_Id, int a_NumVals, int a_StrLen, char *a_PropId);
XMDF_API int xfGetCSName (xid a_Id, int *a_NumVals, int *a_StrLen, char *a_PropId);
XMDF_API int xfGetCSNameLen (xid a_Id, int *a_NumVals, int *a_StrLen);
XMDF_API int xfSetReachName (xid a_Id, int a_NumVals, int a_StrLen, char *a_PropId);
XMDF_API int xfGetReachName (xid a_Id, int *a_NumVals, int *a_StrLen, char *a_PropId);
XMDF_API int xfGetReachNameLen (xid a_Id, int *a_NumVals, int *a_StrLen);
XMDF_API int xfSetTopoID (xid a_Id, int a_NumVals, int a_StrLen, char *a_PropId);
XMDF_API int xfGetTopoID (xid a_Id, int *a_NumVals, int *a_StrLen, char *a_PropId);
XMDF_API int xfGetTopoIDLen (xid a_Id, int *a_NumVals, int *a_StrLen);
XMDF_API int xfSetStation (xid a_Id, int a_NumVals, double *a_PropId, 
                           int a_compression);
XMDF_API int xfGetStation (xid a_Id, int *a_NumVals, double *a_PropId);
XMDF_API int xfSetType (xid a_Id, int a_NumVals, int *a_PropId, 
                        int a_compression);
XMDF_API int xfGetType (xid a_Id, int *a_NumVals, int *a_PropId);
XMDF_API int xfSetpType (xid a_Id, int a_NumVals, int *a_PropId, 
                         int a_compression);
XMDF_API int xfGetpType (xid a_Id, int *a_NumVals, int *a_PropId);
XMDF_API int xfSetpCSDBLink (xid a_Id, int a_NumVals, int *a_PropId, 
                             int a_compression);
XMDF_API int xfGetpCSDBLink (xid a_Id, int *a_NumVals, int *a_PropId);
XMDF_API int xfSetNote (xid a_Id, int a_NumVals, int a_StrLen, char *a_PropId);
XMDF_API int xfGetNote (xid a_Id, int *a_NumVals, int *a_StrLen, char *a_PropId);
XMDF_API int xfGetNoteLen(xid a_Id, int *a_NumVals, int *a_StrLen);

//Geometery
XMDF_API int xfSetXSectGeomX (xid a_Id, int a_index, int a_NumVals, 
                             double *a_iValues, int a_compression);
XMDF_API int xfGetXSectGeomX (xid a_Id, int a_index, int *a_NumVals, 
                             double *a_iValues);
XMDF_API int xfSetXSectGeomY (xid a_Id, int a_index, int a_NumVals, 
                             double *a_iValues, int a_compression);
XMDF_API int xfGetXSectGeomY (xid a_Id, int a_index, int *a_NumVals, 
                             double *a_iValues);
XMDF_API int xfSetXSectGeomD (xid a_Id, int a_index, int a_NumVals, 
                             double *a_iValues, int a_compression);
XMDF_API int xfGetXSectGeomD (xid a_Id, int a_index, int *a_NumVals, 
                             double *a_iValues);
XMDF_API int xfSetXSectGeomZ (xid a_Id, int a_index, int a_NumVals, 
                             double *a_iValues, int a_compression);
XMDF_API int xfGetXSectGeomZ (xid a_Id, int a_index, int *a_NumVals, 
                             double *a_iValues);

//XSecs--Line properties
XMDF_API int xfSetXSectLinePropFrom (xid a_Id, int a_index, int a_NumVals, 
                                     double *a_PropId, int a_compression);
XMDF_API int xfGetXSectLinePropFrom (xid a_Id, int a_index, int *a_NumVals, 
                                     double *a_PropId);
XMDF_API int xfSetXSectLinePropTo (xid a_Id, int a_index, int a_NumVals, 
                                   double *a_PropId, int a_compression);
XMDF_API int xfGetXSectLinePropTo (xid a_Id, int a_index, int *a_NumVals, 
                                   double *a_PropId);
XMDF_API int xfSetXSectLinePropType (xid a_Id, int a_index, int a_NumVals, 
                                     int *a_PropId, int a_compression);
XMDF_API int xfGetXSectLinePropType (xid a_Id, int a_index, int *a_NumVals, 
                                     int *a_PropId);
XMDF_API int xfSetXSectLinePropIValue (xid a_Id, int a_index, int a_NumVals, 
                                       int *a_PropId, int a_compression);
XMDF_API int xfGetXSectLinePropIValue (xid a_Id, int a_index, int *a_NumVals, 
                                       int *a_PropId);
XMDF_API int xfSetXSectLinePropFValue (xid a_Id, int a_index, int a_NumVals, 
                                       double *a_PropId, int a_compression);
XMDF_API int xfGetXSectLinePropFValue (xid a_Id, int a_index, int *a_NumVals, 
                                       double *a_PropId);

//XSec--Line properties--2
XMDF_API int xfSetXSectLinePropID (xid a_Id, int a_NumVals, int *a_PropId, 
                                   int a_compression);
XMDF_API int xfGetXSectLinePropID (xid a_Id, int *a_NumVals, int *a_PropId);
XMDF_API int xfSetXSectLinePropName (xid a_Id, int a_NumVals, int a_StrLen, 
                                     char *a_PropId);
XMDF_API int xfGetXSectLinePropName (xid a_Id, int *a_NumVals, int *a_StrLen, 
                                     char *a_PropId);
XMDF_API int xfGetXSectLinePropNameLen (xid a_Id, int *a_NumVals, int *a_StrLen);
XMDF_API int xfSetXSectLinePropDesc (xid a_Id, int a_NumVals, int a_StrLen, 
                                     char *a_PropId);
XMDF_API int xfGetXSectLinePropDesc (xid a_Id, int *a_NumVals, int *a_StrLen, 
                                     char *a_PropId);
XMDF_API int xfGetXSectLinePropDescLen (xid a_Id, int *a_NumVals, int *a_StrLen);
XMDF_API int xfSetXSectLinePropCategory (xid a_Id, int a_NumVals, int *a_PropId, 
                                         int a_compression);
XMDF_API int xfGetXSectLinePropCategory (xid a_Id, int *a_NumVals, int *a_PropId);
XMDF_API int xfSetXSectLinePropFDefault (xid a_Id, int a_NumVals, double *a_PropId, 
                                         int a_compression);
XMDF_API int xfGetXSectLinePropFDefault (xid a_Id, int *a_NumVals, double *a_PropId);
XMDF_API int xfSetXSectLinePropExclusive (xid a_Id, int a_NumVals, int *a_PropId, 
                                          int a_compression);
XMDF_API int xfGetXSectLinePropExclusive (xid a_Id, int *a_NumVals, int *a_PropId);

XMDF_API int xfSetNumberOfLinePropEnumGroup (xid a_Id, int *a_nPropNum, 
                                             int a_compression);
XMDF_API int xfGetNumberOfLinePropEnumGroup (xid a_Id, int *a_nPropNum);
XMDF_API int xfSetNumberOfLinePropEnum (xid a_Id, int a_index, int *a_nPropNum, 
                                             int a_compression);
XMDF_API int xfGetNumberOfLinePropEnum (xid a_Id, int a_index, int *a_nPropNum);

XMDF_API int xfSetXSectLinePropEnumID (xid a_Id, int a_index, int a_NumVals, 
                                       int *a_PropId, int a_compression);
XMDF_API int xfGetXSectLinePropEnumID (xid a_Id, int a_index, int *a_NumVals, 
                                       int *a_PropId);
XMDF_API int xfSetXSectLinePropEnumMatID (xid a_Id, int a_index, int a_NumVals, 
                                          int *a_PropId, int a_compression);
XMDF_API int xfGetXSectLinePropEnumMatID (xid a_Id, int a_index, int *a_NumVals, 
                                          int *a_PropId);
XMDF_API int xfSetXSectLinePropEnumName (xid a_Id, int a_index, int a_NumVals, 
                                         int a_StrLen, char *a_PropId);
XMDF_API int xfGetXSectLinePropEnumName (xid a_Id, int a_index, int *a_NumVals, 
                                         int *a_StrLen, char *a_PropId);
XMDF_API int xfGetXSectLinePropEnumNameLen (xid a_Id, int a_index, int *a_NumVals, 
                                            int *a_StrLen);
//XSec--Point properties--1
XMDF_API int xfSetXSectPointPropMeasure (xid a_Id, int a_index, int a_NumVals, 
                                         double *a_PropId, int a_compression);
XMDF_API int xfGetXSectPointPropMeasure (xid a_Id, int a_index, int *a_NumVals, 
                                         double *a_PropId);
XMDF_API int xfSetXSectPointPropType (xid a_Id, int a_index, int a_NumVals, int *a_PropId, 
                                      int a_compression);
XMDF_API int xfGetXSectPointPropType (xid a_Id, int a_index, int *a_NumVals, int *a_PropId);

//XSec--Point Properties--2
XMDF_API int xfSetXSectPointPropID (xid a_Id, int a_NumVals, 
                                          int *a_PropId, int a_compression);
XMDF_API int xfGetXSectPointPropID (xid a_Id, int *a_NumVals, 
                                          int *a_PropId);
XMDF_API int xfSetXSectPointPropName (xid a_Id, int a_NumVals, int a_StrLen, 
                                          char *a_PropId);
XMDF_API int xfGetXSectPointPropName (xid a_Id, int *a_NumVals, int *a_StrLen, 
                                          char *a_PropId);
XMDF_API int xfGetXSectPointPropNameLen (xid a_Id, int *a_NumVals, int *a_StrLen);
XMDF_API int xfSetXSectPointPropDesc (xid a_Id, int a_NumVals, int a_StrLen, 
                                          char *a_PropId);
XMDF_API int xfGetXSectPointPropDesc (xid a_Id, int *a_NumVals, int *a_StrLen, 
                                          char *a_PropId);
XMDF_API int xfGetXSectPointPropDescLen (xid a_Id, int *a_NumVals, int *a_StrLen);
XMDF_API int xfSetXSectPointPropExclusive (xid a_Id, int a_NumVals, 
                                           int *a_PropId, int a_compression);
XMDF_API int xfGetXSectPointPropExclusive (xid a_Id, int *a_NumVals, 
                                           int *a_PropId);
//XSec--Profiles
XMDF_API int xfSetXSectTopoName (xid a_Id, int a_NumVals, int a_StrLen, char *a_PropId);
XMDF_API int xfGetXSectTopoName (xid a_Id, int *a_NumVals, int *a_StrLen, char *a_PropId);
XMDF_API int xfGetXSectTopoNameLen (xid a_Id, int *a_NumVals, int *a_StrLen);
XMDF_API int xfSetXSectTopoDesc (xid a_Id, int a_NumVals, int a_StrLen, char *a_PropId);
XMDF_API int xfGetXSectTopoDesc (xid a_Id, int *a_NumVals, int *a_StrLen, char *a_PropId);
XMDF_API int xfGetXSectTopoDescLen (xid a_Id, int *a_NumVals, int *a_StrLen);

// Metadata to XML
XMDF_API int xfGetMetadataXML(xid a_Id, char **a_xml);
XMDF_API int xfGetMetadataProfileXML (xid a_Id, char **a_xml);
XMDF_API int xfGetMetadataSourceXML (xid a_Id, char **a_xml);
XMDF_API int xfGetMetadataSpatialXML (xid a_Id, char **a_xml);

  /* Geometric path group */
XMDF_API int xfCreateGeometricPathGroup(xid a_ParentId, const char *a_Path,
               const char *a_Guid, int a_Compression,
               xid *a_PathGroup, double a_NullVal);
XMDF_API int xfWriteParticleTimestep(xid a_PathGroup, int a_nDim, double a_Time,
               int a_nPaths, double *a_Locs);

XMDF_API int xfGetPathNullVal(xid GroupId, double *NullVal);
XMDF_API int xfGetNumberOfPaths(xid GroupId, int *NumPaths);
XMDF_API int xfGetNumberOfTimes(xid GroupId, int *NumTimes);
XMDF_API int xfGetPathDimensionality(xid GroupId, int *NumDims);
XMDF_API int xfGetPathTimesArray(xid GroupId, int NumTimes, double *Times);
XMDF_API int xfReadPathLocationsAtTime(xid GroupId, int TimeIndex,
                                      int FirstPathIndex,
                                      int NumIndicies, double *Locs);
XMDF_API int xfReadPathLocationsForParticle(xid GroupId, int PathIndex, 
                                   int FirstTimeIndex, int NumTimes,
                                   double *Locs);
XMDF_API int xfReadPathLocationsForParticles(xid GroupId, int NumPaths,
                                    const int *PathIndices,
                                    int FirstTimeIndex, int NumTimes,
                                    double *Locs);

  /* Function to setup evertying necessary to start writing datasets */
XMDF_API int xfSetupToWriteDatasets(const char *a_Filename,
                           const char *MultiDatasetsGroupPath, 
                           const char *PathInMultiDatasetsGroup,
                           const char *SpatialDataObjectGuid,
                           int OverwriteOptions, xid *FileId, xid *GroupId);

  /* Multi-datasets group */
XMDF_API int xfCreateMultiDatasetsGroup(xid a_Id, const char *a_Path,
                                        const char *a_Guid, xid *a_MultiId);
XMDF_API int xfGetGroupPathsSizeForMultiDatasets(xid a_Id, int *Num,
                                                 int *Maxsize);
XMDF_API int xfGetAllGroupPathsForMultiDatasets(xid a_Id, int a_Num,
                                                int a_Maxsize, char *a_Paths);
XMDF_API int xfGetDatasetsSdoGuid(xid a_MultiDatasetsGroup, char *a_GUID);
XMDF_API int xfOpenMultiDatasetsGroup(xid a_Id, xid *DatasetsGroupId);

  /* writing dataset functions */
  /* scalar */
XMDF_API int xfCreateScalarDataset(xid DatasetGroupId, const char *Path,
                    const char *Units, const char *TimeUnits, int Compression,
                    xid *DatasetId);
XMDF_API int xfCreateScalarDatasetExtendable(xid DatasetGroupId, const char *Path,
                    const char *Units, const char *TimeUnits, float aFillVal,
                    int Compression, xid *DatasetId);
XMDF_API int xfExtendScalarDataset(xid DatasetId, int aNewSize);

XMDF_API int xfWriteScalarTimestep(xid Id, double Time, 
                    int NumValues, const float *Values);
XMDF_API int xfWriteScalarTimestepMinMax(xid Id, double Time, int NumValues,
                    const float *Values, float Min, float Max);
XMDF_API int xfSetDatasetNumTimes(xid a_Id, int a_NumTimes);

  /* vector */
XMDF_API int xfCreateVectorDataset(xid DatasetGroupId, const char *Path,
                    const char *Units, const char *TimeUnits, int Compression,
                    xid *DatasetId);

XMDF_API int xfWriteVectorTimestep(xid Id, double Time,
                    int NumValues, int NumComponents, const float *Values);
XMDF_API int xfWriteVectorTimestepMinMax(xid Id, double Time,
                    int NumValues, int NumComponents, const float *Values,
                    float Min, float Max);

  /* activity */
XMDF_API int xfWriteActivityTimestep (xid a_Id, int a_NumActive,
                                         const unsigned char *a_Active);
/* writing/reading reftimes */
XMDF_API int xfDatasetReftime(xid a_Id, double a_Reftime);
XMDF_API int xfUseDatasetReftime(xid a_Id, xmbool *a_bUseReftime);
XMDF_API int xfReadDatasetReftime(xid a_Id, double *a_Reftime);

/* finding paths to datasets */
XMDF_API int xfGetScalarDatasetsInfo(xid a_Id, int *a_Number,
                                        int *a_MaxPathLength);
XMDF_API int xfGetScalarDatasetPaths(xid a_Id, int a_Number,
                                        int a_MaxLength, char *Paths);
XMDF_API int xfGetStationInfo(xid a_Id, int *a_Number,
                                        int *a_MaxPathLength);
XMDF_API int xfGetStationPaths(xid a_Id, int a_Number,
                                        int a_MaxLength, char *Paths);
XMDF_API int xfGetVectorDatasetsInfo(xid a_Id, int *a_Number,
                                        int *a_MaxPathLength);
XMDF_API int xfGetVectorDatasetPaths(xid a_Id, int a_Number,
                                        int a_MaxLength, char *Paths);

  /* reading dataset functions */
XMDF_API int xfGetDatasetNumTimes(xid a_Id, int *a_Numtimes);
XMDF_API int xfGetDatasetNumVals (xid a_Id, int *a_Numvals);
XMDF_API int xfGetDatasetNumActive(xid a_Id, int *a_NumActivevals);
/* for vector datasets */
XMDF_API int xfGetDatasetVecNumComponents(xid a_Id, int *a_NumComponents);

XMDF_API int xfGetDatasetTimeUnits (xid a_Id, char *Units); /* units 100 long */
XMDF_API int xfGetDatasetUnits     (xid a_Id, char *Units); /* units 100 long */

XMDF_API int xfGetVectorDatasetGroupId (xid a_Id);
XMDF_API int xfGetScalarDatasetGroupId (xid a_Id);

  /* Functions to get entire arrays.  Arrays must already be allocated */
XMDF_API int xfGetDatasetTimes(xid a_DatasetId, int a_NumTimes,
                               double *a_Times);
XMDF_API int xfGetDatasetMins(xid a_DatasetId, int a_NumTimes, float *a_Mins);
XMDF_API int xfGetDatasetMaxs(xid a_DatasetId, int a_NumTimes, float *a_Maxs);
XMDF_API int xfGetDatasetMinsFloat(xid a_Id,   int a_NumTimes, float  *a_Mins);
XMDF_API int xfGetDatasetMaxsFloat(xid a_Id,   int a_NumTimes, float  *a_Maxs);
XMDF_API int xfGetDatasetMinsDouble(xid a_Id,  int a_NumTimes, double *a_Mins);
XMDF_API int xfGetDatasetMaxsDouble(xid a_Id,  int a_NumTimes, double *a_Maxs);
XMDF_API int xfGetDatasetActivity(xid a_DatasetId, unsigned char *a_Active);

  /* Functions to get portions of datasets.  Arrays must already be allocated */
XMDF_API int xfReadScalarValuesTimestep(xid a_DatasetId, 
                    int a_TimestepIndex, int a_NumVals, float *a_Values);
XMDF_API int xfReadScalarValuesTimestepFloat(xid a_DatasetId, 
                    int a_TimestepIndex, int a_NumVals, float *a_Values);
XMDF_API int xfReadScalarValuesTimestepFloatPortion (xid a_Id, 
                    int a_TimestepIndex, int a_Start, int a_NumVals,
                    float *a_Values);
XMDF_API int xfReadScalarValuesTimestepDouble(xid a_Id, 
                  int a_TimestepIndex, int a_NumVals, double *a_Values);
XMDF_API int xfReadScalarValuesTimestepDoublePortion(xid a_Id, 
                  int a_TimestepIndex, int a_Start, int a_NumVals,
                  double *a_Values);
XMDF_API int xfReadScalarValuesTimestepInt(xid a_Id, 
                  int a_TimestepIndex, int a_NumVals, int *a_Values);
XMDF_API int xfReadScalarValuesAtIndex(xid a_DatasetId,
                    int a_Index, int a_FirstTimestep, int a_NumTimes,
                    float *a_Values);
XMDF_API int xfReadScalarValuesAtIndexFloat(xid a_DatasetId,
                    int a_Index, int a_FirstTimestep, int a_NumTimes,
                    float *a_Values);
XMDF_API int xfReadScalarValuesAtIndexDouble(xid a_Id, 
                    int a_Index, int a_FirstTime, int a_NumTimes,
                    double *a_Values);
XMDF_API int xfReadScalarValuesAtIndices(xid a_Id, 
                    int a_nIndices, const int *a_Indices, int a_FirstTime,
                    int a_NumTimes, float *a_Values);
XMDF_API int xfReadScalarValuesAtIndicesFloat(xid a_Id, 
                    int a_nIndices, const int *a_Indices, int a_FirstTime,
                    int a_NumTimes, float *a_Values);
XMDF_API int xfReadVectorValuesTimestep(xid a_DatasetId,
                    int a_TimestepIndex, int a_NumVals, int a_NumComponents,
                    float *a_Values);
XMDF_API int xfReadVectorValuesTimestepFloat(xid a_Id, 
                    int a_TimestepIndex, int a_NumVals, int a_NumComponents, 
                    float *a_Values);
XMDF_API int xfReadVectorValuesTimestepFloatPortion(xid a_Id, 
                    int a_TimestepIndex, int a_Start, int a_NumVals,
                    int a_NumComponents, float *a_Values);
XMDF_API int xfReadVectorValuesTimestepDouble(xid a_Id, 
                    int a_TimestepIndex, int a_NumVals, int a_NumComponents, 
                    double *a_Values);
XMDF_API int xfReadVectorValuesAtIndex(xid a_Id, int a_Index, 
                    int a_FirstTime, int a_NumTimes, int a_NumComponents,
                    float *a_Values);
XMDF_API int xfReadVectorValuesTimestepDoublePortion (xid a_Id, 
                    int a_TimestepIndex, int a_Start, int a_NumVals,
                    int a_NumComponents, double *a_Values);
XMDF_API int xfReadVectorValuesAtIndexFloat(xid a_Id, int a_Index,
                    int a_FirstTime, int a_NumTimes, int a_NumComponents,
                    float *a_Values);

XMDF_API int xfReadActivityTimestep(xid a_Id, int Index, int NumActive,
                    xmbool *Active);
XMDF_API int xfReadActivityValuesAtIndex(xid a_Id, int a_Index,
                    int a_FirstTime, int a_NumTimes, xmbool *a_Values);
 
XMDF_API hssize_t xfGetParallelNumValuesToRead();
XMDF_API void xfSetParallelNumValuesToRead(hssize_t a);
XMDF_API void xfGetParIndexFromIndex(xid a_Id, int a_TimestepIndex,
                             int a_index, int *a_parindex);

  /* reserved dataset property functions */
XMDF_API int xfScalarDataLocation(xid a_Id, int a_DataLoc);
XMDF_API int xfVector2DDataLocations(xid a_Id, int a_DataLocI, int a_DataLocJ);
XMDF_API int xfVector3DDataLocations(xid a_Id, int a_DataLocI, int a_DataLocJ,
                                     int a_DataLocK);
XMDF_API int xfGetScalarDataLocation(xid a_Id, int *a_DataLoc);
XMDF_API int xfGetVector2DDataLocations(xid a_Id, int *a_DataLocI,
                               int *a_DataLocJ);
XMDF_API int xfGetVector3DDataLocations(xid a_Id, int *a_DataLocI,
                               int *a_DataLocJ, int *a_DataLocK);
XMDF_API int xfVectorsInLocalCoords (xid a_Id);
XMDF_API int xfAreVectorsInLocalCoords (xid a_Id, int *a_LocalCoords);

/* Coordinate conversions */
XMDF_API int xfGetHorizDatum(xid a_CoordId, int *a_val);
XMDF_API int xfGetHorizUnits(xid a_CoordId, int *a_val);
XMDF_API int xfGetVertDatum(xid a_CoordId, int *a_val);
XMDF_API int xfGetVertUnits(xid a_CoordId, int *a_val);
XMDF_API int xfGetLat(xid a_CoordId, int *a_val);
XMDF_API int xfGetLon(xid a_CoordId, int *a_val);
XMDF_API int xfGetUTMZone(xid a_CoordId, int *a_val);
XMDF_API int xfGetSPCZone(xid a_CoordId, int *a_val);
XMDF_API int xfGetHPGNArea(xid a_CoordId, int *a_val);
XMDF_API int xfGetCPPLat(xid a_CoordId, double *a_val);
XMDF_API int xfGetCPPLon(xid a_CoordId, double *a_val);
XMDF_API int xfGetEllipse(xid a_CoordId, int *a_val);
XMDF_API int xfGetMajorR(xid a_CoordId, double *a_val);
XMDF_API int xfGetMinorR(xid a_CoordId, double *a_val);

XMDF_API int xfSetHorizDatum(xid a_CoordId, int a_val);
XMDF_API int xfSetHorizUnits(xid a_CoordId, int a_val);
XMDF_API int xfSetVertDatum(xid a_CoordId, int a_val);
XMDF_API int xfSetVertUnits(xid a_CoordId, int a_val);
XMDF_API int xfSetLat(xid a_CoordId, int a_val);
XMDF_API int xfSetLon(xid a_CoordId, int a_val);
XMDF_API int xfSetUTMZone(xid a_CoordId, int a_val);
XMDF_API int xfSetSPCZone(xid a_CoordId, int a_val);
XMDF_API int xfSetHPGNArea(xid a_CoordId, int a_val);
XMDF_API int xfSetCPPLat(xid a_CoordId, double a_val);
XMDF_API int xfSetCPPLon(xid a_CoordId, double a_val);
XMDF_API int xfSetEllipse(xid a_CoordId, int a_val);
XMDF_API int xfSetMajorR(xid a_CoordId, double a_val);
XMDF_API int xfSetMinorR(xid a_CoordId, double a_val);

#ifdef __cplusplus
}

#define EMRL_XMDF

/******************************************************************************
 * CLASS     ECreateGroup
 * PURPOSE   calls xfiCreateGroup, for EMRL only!
 * NOTES     
 ******************************************************************************/
class ECreateGroup
{
  public:
    ECreateGroup() {;}
    ~ECreateGroup() {;}

    int CallXfiCreateGroup(xid a_Id, const char *Path, 
                           const char *a_GroupType);
};
#endif


#endif
