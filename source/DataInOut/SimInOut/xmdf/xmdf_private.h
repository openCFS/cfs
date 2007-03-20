#ifndef XMDF_PRIVATE_DOT_H
#define XMDF_PRIVATE_DOT_H

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

/* WARNING: This file is intended to be used only by EMRL emplyees.  
           The contents and function calls may change without warning.*/

#include "Xmdf.h"

/* Constants/enums */
#define XFBitOn(flags, bit) ((flags) |= (bit))
#define XFBitOff(flags, bit) ((flags) &= ~(bit))
#define XFBitSet(flags, bit, TorF) ((TorF) ? ((flags) |= (bit)) : ((flags) &= ~(bit)))
#define XFBitTest(flags, bit) ((((flags) & (bit)) == 0) ? 0 : 1 )

#define FILE_TYPE_ID "File Type"
#define FILE_TYPE_XMDF "Xmdf"

#define FILE_CONTENTS "Contents"
#define FILE_CONTENTS_DATASETS "Datasets"

#define GROUP_TYPE "Grouptype"
#define GROUP_TYPE_GENERIC "Generic"
#define GROUP_TYPE_STATION "Station"
#define GROUP_TYPE_MESH "MESH"
#define GROUP_TYPE_GRID "GRID"
#define GROUP_TYPE_XSECS "XSECS"
#define GROUP_TYPE_GEOMETRIC_PATH "GEOMETRIC PATHS"
#define GROUP_TYPE_PROPERTIES "PROPERTIES"
#define GROUP_TYPE_MULTI_DATASETS "MULTI DATASETS"
#define GROUP_TYPE_DATASET_SCALAR "DATASET SCALAR"
#define GROUP_TYPE_DATASET_VECTOR "DATASET VECTOR"
#define GROUP_TYPE_MESH_NODES "NodeGroup"
#define GROUP_TYPE_MESH_ELEMS "ElemGroup"
#define GROUP_TYPE_COORDS "Coordinates"
#define GROUP_TYPE_XSEC_GEOM "GeomGroup"
#define GROUP_TYPE_XSEC_LINE "LineGroup"
#define GROUP_TYPE_XSEC_POINT "PointGroup"
#define GROUP_TYPE_XSEC_PROFILES "ProfilesGroup"
#define GROUP_TYPE_XSEC_XSECS "XSECSGroup"
#define GROUP_TYPE_MAT_SIM "Material Simulation"

#define DEFAULT_GROUP_SIZE 1

#define DATASET_DSET_VALUES  "Values"
#define DATASET_DSET_TIMES   "Times"
#define DATASET_DSET_MINS    "Mins"
#define DATASET_DSET_MAXS    "Maxs"
#define DATASET_DSET_ACTIVE  "Active"
#define DATASET_REFTIME "Reftime"
#define DATASET_NUMTIMES "NumTimes"
#define DATASET_DSET_INDICES "Indices"

#define DATASET_TYPE         "Data Type"
#define DATASET_TYPE_FLOAT   0
#define DATASET_TYPE_DOUBLE  1
#define DATASET_TYPE_INT     2

#define DATASET_TIME_CHUNK_SIZE 10
#define MESH_ELEM_TYPE_CHUNK_SIZE 200

/* name for multiple datasets folder under a mesh, grid, or scattered dataset */
#define MULTI_DATASET_LOCATION "Datasets"

/* Attributes used for meshes */
#define MESH_NODE_GROUP   "Nodes"
#define MESH_ELEM_GROUP   "Elements"
#define MESH_ATT_NUMELEMS "NumElems"
#define MESH_ATT_NUMNODES "NumNodes"
#define MESH_DSET_NODE_LOCS "NodeLocs"
#define MESH_DSET_NODE_BOUND "BoundNodes"
#define MESH_DSET_NODE_SHARED_SUBDOMAIN "SharedSubdomainNodes"
#define MESH_DSET_ELEM_TYPES "Types"
#define MESH_DSET_ELEM_NODE_IDS "Nodeids"

#define MESH_PROP_GROUP "PROPERTIES"
#define MESH_NODE_PROP_GROUP "Nodes/PROPERTIES"
#define MESH_ELEM_PROP_GROUP "Elements/PROPERTIES"

/* Attributes & datasets used for grids  */
#define GRID_ATT_TYPE "Gridtype"
#define GRID_ATT_NUM_DIMENSIONS "Numdimensions"
#define GRID_ATT_EXTRUDE_TYPE "Extrudetype"
#define GRID_ATT_ORIGIN "Origin"
#define GRID_ATT_ORIENTATION "Orientation"
#define GRID_ATT_BEARING "Bearing"
#define GRID_ATT_DIP "Dip"
#define GRID_ATT_ROLL "Roll"
#define GRID_ATT_COMP_ORIG "CompOrigin"
#define GRID_ATT_U_DIR "Direction"

#define GRID_ATT_NUM_I "NumI"
#define GRID_ATT_NUM_J "NumJ"
#define GRID_ATT_NUM_K "NumK"
#define GRID_ATT_EXTRUDE_LAYERS "Extrudelayers"
#define GRID_DSET_COORDS_I "CoordsI"
#define GRID_DSET_COORDS_J "CoordsJ"
#define GRID_DSET_COORDS_K "CoordsK"
#define GRID_DSET_EXTRUDE_VALS "Extrudevals"

#define GRID_PROP_GROUP "PROPERTIES"
#define GRID_CELL_PROP_GROUP "GridCellProps"
#define GRID_NODE_PROP_GROUP "GridNodeProps"

/* Attributes used for geometric paths */
#define GEOMPATH_ATT_COMPRESSION "GeompathCompression"
#define GEOMPATH_ATT_NULLVALUE "GeompathNullvalue"
#define GEOMPATH_ATT_MINS "GeompathMins"
#define GEOMPATH_ATT_MAXS "GeompathMaxs"
#define GEOMPATH_DSET_LOCS  "Locations"
#define GEOMPATH_DSET_TIMES   "Times"
#define GEOMPATH_TIME_CHUNK_SIZE 10

/* Attributes used for xsecs */
#define XMDF_XSEC                "Cross Section"
#define XMDF_XSECT               "XSect"
#define XMDF_XSECT_PROPNUM       "Line Prop. Number"
#define XMDF_XSECT_PROPNUM2      "Material Enums"
#define XMDF_XSECT_PROPENUMNUM   "Number of Enums"

#define XSEC_GEOMETRY_GROUP "CrossSectionGeometry"
#define XSEC_LINE_GROUP "LineProperties"
#define XSEC_POINT_GROUP "PointProperties"
#define XSEC_PROFILES_GROUP "Profiles"
#define XSEC_XSECS_GROUP "XSECS"
//XSects
#define XMDF_CSID                "CS ID"
#define XMDF_CSNAME              "CS Name"
#define XMDF_REACHNAME           "Reach Name"
#define XMDF_TOPO_ID             "Topo ID"
#define XMDF_STATION             "Station"
#define XMDF_PTYPE               "pType"
#define XMDF_TYPE                "Type"
#define XMDF_CSDBLINK            "CSDBlink"
#define XMDF_NOTE                "Note"
//Geometry
#define XMDF_GEOMETRY_X          "Geometry X"
#define XMDF_GEOMETRY_Y          "Geometry Y"
#define XMDF_GEOMETRY_D          "Geometry D"
#define XMDF_GEOMETRY_Z          "Geometry Z"
//Line Properties--1
#define XMDF_LP_FROM             "Line Prop. From"
#define XMDF_LP_TO               "Line Prop. To"
#define XMDF_LP_TYPE             "Line Prop. Type"
#define XMDF_LP_IVALUE           "Line Prop. I-Value"
#define XMDF_LP_FVALUE           "Line Prop. F-Value"
//Point Properties--1
#define XMDF_PP_MEASURE          "Point Prop. Measure"
#define XMDF_PP_TYPE             "Point Prop. Type"
//Topo--Profile
#define XMDF_TOPO_NAME           "Topo Name"
#define XMDF_TOPO_DESC           "Topo Desc"
//Point Properties--2
#define XMDF_PTID                "Point Prop. ID"
#define XMDF_PP_NAME             "Point Prop. Name"
#define XMDF_PP_DESC             "Point Prop. Desc"
#define XMDF_PP_EXCLUSIVE        "Point Prop. Exclusive"
//Line Properties--2
#define XMDF_LNID                "Line Prop. ID"
#define XMDF_LN_NAME             "Line Prop. Name" 
#define XMDF_LN_DESC             "Line Prop. Description"
#define XMDF_LN_CATEGORY         "Line Prop. Category" 
#define XMDF_LN_FDEFAULT         "Line Prop. F-Default"
#define XMDF_LN_EXCLUSIVE        "Line Prop. Exclusive"
#define XMDF_LN_PROPID           "Line Prop. ID"
#define XMDF_LN_PROPENUMID       "Line Prop. Enum"
#define XMDF_LN_PROPNAME         "Line Prop. Name"

/* Metadata */
#define METADATA_MAIN "Metadata"
#define METADATA_MAIN_TITLE "Title"
#define METADATA_MAIN_ABSTRACT "Abstract"
#define METADATA_MAIN_PURPOSE "Purpose"
#define METADATA_MAIN_STARTDATE "Start Date"
#define METADATA_MAIN_TOPICS "Generic Topics"
#define METADATA_MAIN_COMMENTS "Generic Comments"
#define METADATA_MAIN_PROFILE "Profile"
#define METADATA_PROFILE_ADDRESS "Address"
#define METADATA_PROFILE_ADDRESSTYPE "Address Type"
#define METADATA_PROFILE_CONTACTPERSON "Contact Person"
#define METADATA_PROFILE_CONTACTPOSITION "Contact Position"
#define METADATA_PROFILE_EMAIL "Email"
#define METADATA_PROFILE_FAX "Fax"
#define METADATA_PROFILE_ORGANIZATION "Organization"
#define METADATA_PROFILE_PHONE "Phone"
#define METADATA_MAIN_SOURCE "Source"
#define METADATA_SOURCE_USERDEFINED "User Defined"
#define METADATA_SOURCE_XMSGENERATED "XMS Generated"
#define METADATA_MAIN_SPATIAL "Spatial"
#define METADATA_SPATIAL_BOTTOMBOUND "Bottom Bound"
#define METADATA_SPATIAL_BOUNDINGPOLY "BoundingPoly"
#define METADATA_SPATIAL_LEFTBOUND "Left Bound"
#define METADATA_SPATIAL_RIGHTBOUND "Right Bound"
#define METADATA_SPATIAL_TOPBOUND "Top Bound"

/* Metadata tags */
#define METATAG_MAIN "<metadata>"
#define METATAG_MAIN_END "</metadata>"
#define METATAG_MAIN_TITLE "<title>"
#define METATAG_MAIN_TITLE_END "</title>"
#define METATAG_MAIN_ABSTRACT "<abstract>"
#define METATAG_MAIN_ABSTRACT_END "</abstract>"
#define METATAG_MAIN_PURPOSE "<purpose>"
#define METATAG_MAIN_PURPOSE_END "</purpose>"
#define METATAG_MAIN_STARTDATE "<startdate>"
#define METATAG_MAIN_STARTDATE_END "</startdate>"
#define METATAG_MAIN_TOPIC "<generictopic>"
#define METATAG_MAIN_TOPIC_END "</generictopic>"
#define METATAG_MAIN_COMMENT "<genericcomment>"
#define METATAG_MAIN_COMMENT_END "</genericcomment>"
#define METATAG_MAIN_PROFILE "<profile>"
#define METATAG_MAIN_PROFILE_END "</profile>"
#define METATAG_PROFILE_ADDRESS "<address>"
#define METATAG_PROFILE_ADDRESS_END "</address>"
#define METATAG_PROFILE_ADDRESSTYPE "<addresstype>"
#define METATAG_PROFILE_ADDRESSTYPE_END "</addresstype>"
#define METATAG_PROFILE_CONTACTPERSON "<contactperson>"
#define METATAG_PROFILE_CONTACTPERSON_END "</contactperson>"
#define METATAG_PROFILE_CONTACTPOSITION "<contactposition>"
#define METATAG_PROFILE_CONTACTPOSITION_END "</contactposition>"
#define METATAG_PROFILE_EMAIL "<email>"
#define METATAG_PROFILE_EMAIL_END "</email>"
#define METATAG_PROFILE_FAX "<fax>"
#define METATAG_PROFILE_FAX_END "</fax>"
#define METATAG_PROFILE_ORGANIZATION "<organization>"
#define METATAG_PROFILE_ORGANIZATION_END "</organization>"
#define METATAG_PROFILE_PHONE "<phone>"
#define METATAG_PROFILE_PHONE_END "</phone>"
#define METATAG_MAIN_SOURCE "<source>"
#define METATAG_MAIN_SOURCE_END "</source>"
#define METATAG_SOURCE_USERDEFINED "<userdefined>"
#define METATAG_SOURCE_USERDEFINED_END "</userdefined>"
#define METATAG_SOURCE_XMSGENERATED "<xmsgenerated>"
#define METATAG_SOURCE_XMSGENERATED_END "</xmsgenerated>"
#define METATAG_MAIN_SPATIAL "<spatial>"
#define METATAG_MAIN_SPATIAL_END "</spatial>"
#define METATAG_SPATIAL_BOTTOMBOUND "<bottombound>"
#define METATAG_SPATIAL_BOTTOMBOUND_END "</bottombound>"
#define METATAG_SPATIAL_BOUNDINGPOLY "<boundingpoly>"
#define METATAG_SPATIAL_BOUNDINGPOLY_END "</boundingpoly>"
#define METATAG_SPATIAL_LEFTBOUND "<leftbound>"
#define METATAG_SPATIAL_LEFTBOUND_END "</leftbound>"
#define METATAG_SPATIAL_RIGHTBOUND "<rightbound>"
#define METATAG_SPATIAL_RIGHTBOUND_END "</rightbound>"
#define METATAG_SPATIAL_TOPBOUND "<topbound>"
#define METATAG_SPATIAL_TOPBOUND_END "</topbound>"

/* Attributes used for datasets */
#define DATASET_ATT_UNITS "DatasetUnits"
#define DATASET_ATT_EXTENDFILL "DatasetExtendFillValue"
#define DATASET_ATT_COMPRESSION "DatasetCompression"
#define DATASET_ATT_NULLVALUE "DatasetNullvalue"
#define DATASET_ATT_DATALOCATION "DatasetLocation"
#define DATASET_ATT_DATALOCATIONI "DatasetLocationI"
#define DATASET_ATT_DATALOCATIONJ "DatasetLocationJ"
#define DATASET_ATT_DATALOCATIONK "DatasetLocationK"
#define DATASET_ATT_LOCALCOORDS   "LocalCoords"

/* Attributes used for coordinate groups */
#define COORD_ATT_HORIZ_DATUM   "HorizontalDatum"
#define COORD_ATT_HORIZ_UNITS   "HorizontalUnits"
#define COORD_ATT_VERT_DATUM    "VerticalDatum"
#define COORD_ATT_VERT_UNITS    "VerticalUnits"
#define COORD_ATT_LATITUDE      "Latitude"
#define COORD_ATT_LONGITUDE     "Longitude"
#define COORD_ATT_UTM_ZONE      "UtmZone"
#define COORD_ATT_SPC_ZONE      "SpcZone"
#define COORD_ATT_HPGN_AREA     "HpgnArea"
#define COORD_ATT_CPP_LATITUDE  "Latitude"
#define COORD_ATT_CPP_LONGITUDE "Longitude"
#define COORD_ATT_ELLIPSE       "Ellipse"
#define COORD_ATT_MAJOR_R       "MajorR"
#define COORD_ATT_MINOR_R       "MinorR"

/* Used to determine chunking for properties */
#define XF_PROP_CHUNK_NUMBER   100
#define XF_PROP_CHUNK_MIN       10
#define XF_PROP_CHUNK_MAX     1000
#define XF_DATASET_MIN_COMPRESS_SIZE  50

/* Bitwise flags associated with XDatasetParams */
#define XF_DSET_SZIP 1
#define XF_DSET_SHUFFLE 2

/* structures */
typedef struct {
  const char *GroupType;
  xmbool       bCounting; /* Set to XTRUE to count groups and max length
                            Set to XFALSE to copy paths into array */
  xmbool       bLookInSubgroups; /* Look in subgroups of folders of specific 
                                   types */
  /* When recursing through groups we determine the groups and max length first */
  int         nGroups;
  int         MaxStringLength;
  /* After obtaining the number of groups and max length we send them back in */
  /* along with the address of a 2D character array to fill in for the paths */
  int         iCurPath;
  char       *Paths;
} XGroupIteration;

/* This structure is used to make setting up and writing datasets simplier */
/* Use the internal functions to modify the class */
typedef struct {
  int     Rank;
  hsize_t *Dims;
  hsize_t *Maxdims;
  xmbool   bChunked;
  hsize_t *Chunksize;
      /* NONE or 0 - 9 */
  int   Compression;
  unsigned int   bwFlags; /* XF_DSET_SHUFFLE etc. */
      /* setup to use fill value */
  xmbool    bUseFillValue;
  hid_t    fill_type;
  void    *fill_value;
} XDatasetParams;

#ifdef __cplusplus
extern "C" {
#endif

/* errors */
XMDF_API int xfpGetNumErrorMessages (int *a_Num);
XMDF_API int xfpGetErrorMessages (int a_Num, char **a_Errors);
XMDF_API void xfpAllocateErrorStack();
XMDF_API void xfpAddXMDFError(const char * a_Error);
XMDF_API void xfpClearErrors();
XMDF_API herr_t xfpHDF5ErrorHandler (void* client_data);
XMDF_API herr_t xfpHDF5ErrorWalk_cb (int n, H5E_error_t *err_desc, void *client_data);

/* close open identifiers when closing a file */
XMDF_API int xfpCloseOpenIdentifiers(xid a_File);

/* clear all objects in an HDF5 group */
XMDF_API int xfpClearGroup(xid a_Group);

/* Double/Float little endian and big endian */
XMDF_API int  xfpGetInitialized ();
XMDF_API void xfpInitialize ();
XMDF_API void xfpDefaultWriteFloatType (int a_BigEndian);
XMDF_API hid_t xfpGetDefaultDoubleType ();
XMDF_API hid_t xfpGetDefaultFloatType ();

/* working with HDF5 datasets */
XMDF_API int xfpReadDset1D(xid a_Id, const char *a_Name, int a_Number,
                   H5T_class_t a_DsetClass, xid a_MemType, void *a_Properties);
XMDF_API int xfpReadDset1DDouble(xid a_Id, const char *a_Name, int a_Num,
                                double *a_Array);
XMDF_API int xfpReadDset1DFloat (xid a_Id, const char *a_Name, int a_Num,
                                float  *a_Array);
XMDF_API int xfpReadDsetDouble (xid a_Id, const char *a_Name, int a_Rank,
                                const hsize_t *Dims, double *a_Array);
XMDF_API int xfpReadDsetInt (xid a_Id, const char *a_Name, int a_Rank,
                                const hsize_t *Dims, int *a_Array);
  /* portions of HDF5 datasets */
XMDF_API int xfpReadDset1DIntPortion(xid a_Id, const char *a_Name,
                hssize_t a_Start1, hsize_t a_Num1, int *a_Array);
XMDF_API int xfpReadDset1DDoublePortion(xid a_Id, const char *a_Name,
                hssize_t a_Start1, hsize_t a_Num1, double *a_Array);
XMDF_API int xfpWriteDset1DFloatPortion(xid a_Id, const char *a_Name,
                int a_Start1, int a_Num1, const float *a_Array);
XMDF_API int xfpWriteDset2DFloatPortion(xid a_Id, const char *a_Name,
                int a_Start1, int a_Num1, int a_Start2, int a_Num2,
                const float *a_Array);
XMDF_API int xfpReadDset2DFloatPortion (xid a_Id, const char *a_Name, 
                hssize_t a_Start1, hsize_t a_Num1,
                hssize_t a_Start2, hsize_t a_Num2, float *a_Array);
XMDF_API int xfpReadDset2DFloatPortions (xid a_Id, const char *a_Name, 
                int nPortions, hssize_t *a_Start1, hsize_t *a_Num1,
                hssize_t *a_Start2, hsize_t *a_Num2, float *a_Array);
XMDF_API int xfpReadDset2DDoublePortion (xid a_Id, const char *a_Name, 
                hssize_t a_Start1, hsize_t a_Num1,
                hssize_t a_Start2, hsize_t a_Num2, double *a_Array);
XMDF_API int xfpReadDset2DIntPortion (xid a_Id, const char *a_Name, 
                hssize_t a_Start1, hsize_t a_Num1,
                hssize_t a_Start2, hsize_t a_Num2, int *a_Array);
XMDF_API int xfpReadDset2DUCharPortion (xid a_Id, const char *a_Name, 
                hssize_t a_Start1, hsize_t a_Num1,
                hssize_t a_Start2, hsize_t a_Num2, xmbool *a_Array);
XMDF_API int xfpReadDset3DDoublePortion (xid a_Id, const char *a_Name, 
                hssize_t a_Start1, hsize_t a_Num1,
                hssize_t a_Start2, hsize_t a_Num2, hssize_t a_Start3,
                hsize_t a_Num3, double *a_Array);
XMDF_API int xfpReadDset3DFloatPortion (xid a_Id, const char *a_Name, 
                hssize_t a_Start1, hsize_t a_Num1,
                hssize_t a_Start2, hsize_t a_Num2, hssize_t a_Start3,
                hsize_t a_Num3, float *a_Array);
XMDF_API int xfpReadDset3DDoublePortions (xid a_Id, const char *a_Name,
                int nPortions, hssize_t *a_Start1, hsize_t *a_Num1,
                hssize_t *a_Start2, hsize_t *a_Num2, hssize_t *a_Start3,
                hsize_t *a_Num3, double *a_Array);

  /* Write arrays to datasets */
XMDF_API int xfpWriteDsetInt(xid a_Id, const char *a_Name, hid_t a_Datatype,
                            const XDatasetParams *a_Params, const int *a_Array);
XMDF_API int xfpWriteDsetUInt(xid a_Id, const char *a_Name,
                            hid_t a_Datatype, const XDatasetParams *a_Params,
                            const unsigned int *a_Array);
XMDF_API int xfpWriteDsetFloat(xid a_Id, const char *a_Name,
                   const XDatasetParams *a_Params, const float *a_Array);
XMDF_API int xfpWriteDsetDouble (xid a_Id, const char *a_Name,
                         const XDatasetParams *a_Params, const double *a_Array);
XMDF_API int xfpWriteDsetUChar(xid a_Id, const char *a_Name, hid_t a_Datatype,
                   const XDatasetParams *a_Params, const unsigned char*a_Array);
XMDF_API int xfpAppendDset1DDouble(xid a_Id, const char *a_Name,
                                  hsize_t NewTimes, const double *a_Array);
XMDF_API int xfpAppendDset1DFloat(xid a_Id, const char *a_Name,
                                  hsize_t a_NumNew, const float *a_Array);
XMDF_API int xfpAppendDset1DInt(xid a_Id, const char *a_Name, hsize_t a_NumNew,
                                const int *a_Array);
XMDF_API int xfpAppendDset2DFloatFirstDim(xid a_Id, const char *a_Name, 
                                   hsize_t a_NumNew, hsize_t a_NumDim2,
                                   const float *a_Array);
XMDF_API int xfpAppendDset2DDoubleFirstDim(xid a_Id, const char *a_Name, 
                                   hsize_t a_NumNew, hsize_t a_NumDim2,
                                   const double *a_Array);
XMDF_API int xfpAppendDset2DDoubleSecondDim (xid a_Id, const char *a_Name, 
                                   hsize_t a_NumNew, hsize_t a_NumDim1,
                                   const double *a_Array);
XMDF_API int xfpAppendDset3DDoubleFirstDim(xid a_Id, const char *a_Name, 
                                   hsize_t a_NumNew, hsize_t a_NumDim2,
                                   hsize_t a_NumDim3,
                                   const double *a_Array);
XMDF_API int xfpAppendDset3DFloatFirstDim(xid a_Id, const char *a_Name, 
                                   hsize_t a_NumNew, hsize_t a_NumDim2,
                                   hsize_t a_NumDim3,
                                   const float *a_Array);
XMDF_API int xfpAppendDset2DUCharFirstDim (xid a_Id, const char *a_Name, 
                                   hsize_t a_NumNew, hsize_t a_NumDim2,
                                   const unsigned char *a_Array);


  /* string datasets */
XMDF_API int xfpReadDatasetString (xid a_Loc, const char *a_Name, char **a_Str);

XMDF_API int xfpWriteDatasetString (xid a_Loc, const char *a_Name, 
                                  const char *a_Str);
XMDF_API int xfpWriteDatasetStrings (xid a_Loc, const char *a_Name, 
                 hsize_t a_Number, size_t a_StringLength, const char *a_Str);

/* working with HDF5 attributes */
XMDF_API int xfpDoesAttributeExist (xid a_Loc, const char *a_Name,
                                  xmbool *a_bExists);

XMDF_API int xfpReadAttributeString (xid a_Loc, const char *a_Name, int a_Num,
                                   char **a_Str);
XMDF_API int xfpReadAttributeInt (xid a_Loc, const char *a_Name, int a_Num,
                                int *a_val);
XMDF_API int xfpReadAttributeFloat (xid a_Loc, const char *a_Name, int a_Num,
                                  float *a_val);
XMDF_API int xfpReadAttributeDouble (xid a_Loc, const char *a_Name, int a_Num,
                                   double *a_val);

XMDF_API int xfpWriteAttributeString (xid a_Loc, const char *a_Name, 
                                    const char *a_Str);
XMDF_API int xfpWriteAttributeInt (xid a_Loc, const char *a_Name, int a_Number,
                                 int *a_val);
XMDF_API int xfpWriteAttributeFloat (xid a_Loc, const char *a_Name, int a_Number,
                                   float *a_val);
XMDF_API int xfpWriteAttributeDouble (xid a_Loc, const char *a_Name, int a_Number,
                                    double *a_val);

/* working with HDF5 data spaces  */
XMDF_API int xfpGetSimpleDataspaceInfoFromName (xid a_Loc, const char *a_name,
                           int *a_Rank, hsize_t **a_Dims, hsize_t **a_MaxDims);
XMDF_API int xfpGetSimpleDataspaceInfo (xid a_Id, int *a_Rank, hsize_t **a_Dims, 
                               hsize_t **a_MaxDims);
XMDF_API void xfpDestroySimpleDataspaceInfo(hsize_t **a_Dims,
                                            hsize_t **a_MaxDims);
/* Group functions */
XMDF_API int xfpCreateGroup(xid a_Id, const char *Path, const char *GroupType);
XMDF_API herr_t xfpRecurseCountGroups(xid a_Id, XGroupIteration *a_CountObj);
XMDF_API int xfpNumGroupsOfType(xid a_Id, const char *a_GroupType,
                       int *a_NumGroups, int *a_MaxLength, 
                       xmbool a_bLookInSubgroups);
XMDF_API int xfpPathGroupsOfType(xid a_Id, const char *a_GroupType, 
                        int a_NumPaths, int a_MaxPathLength, char *a_Paths,
                        xmbool a_bLookInSubgroups);

XMDF_API int xfpGetGroupType (xid a_Id, char **a_GroupType);
XMDF_API xmbool xfpIsGroupOfType(xid a_Id, const char *GroupType);

XMDF_API int xfpGetAbsolutePath (xid a_Id, char **a_Path);

/* The following functions are used with the struct XDatasetParams */
  /* always call this function before using */
XMDF_API void xfpDsetParamsInit(XDatasetParams *a_Params, int a_Rank, 
                              xmbool bChunked, int a_Compression);
  /* call this function to reset existing structures */
XMDF_API void xfpDsetParamsReset(XDatasetParams *a_Params, int a_Rank, 
                               xmbool bChunked, int a_Compression);
XMDF_API void xfpDsetParamsDestroy(XDatasetParams *a_Params);
  /* use this function to set the sizes for the parameters */
  /* chunkdim ignored if no chunking */
XMDF_API int xfpDsetParamsSetSizes(XDatasetParams *a_Params, int a_Index, 
                         hsize_t a_Dim, hsize_t a_Maxdim, hsize_t a_Chunkdim);
XMDF_API int xfpDsetParamsCreateSpace (const XDatasetParams *a_Params);
XMDF_API int xfpDsetParamsCreateProps (const XDatasetParams *a_Params);

XMDF_API int xfpDsetParamsUseFillValue(XDatasetParams *a_Params, xid a_type,
                void *a_value);
XMDF_API int xfpDsetParamsSetSZip(XDatasetParams *a_Params, xmbool a_on);
XMDF_API int xfpDsetParamsSetShuffle(XDatasetParams *a_Params, xmbool a_on);

/*parallel stuff */
XMDF_API xmbool   xfpGetRunningParallel ();
XMDF_API void     xfpSetParallel (int a_rank);
XMDF_API int      xfpGetParallelRank ();
XMDF_API hssize_t xfpGetParallelNumValuesToRead();
XMDF_API void     xfpSetParallelNumValuesToRead(hssize_t a);
XMDF_API double*  xfpStridedToContinousDouble (const double *a_OldArray,
                                       const int a_NumValues, const int a_Stride,
                                       const int a_NumPerStride, const int a_StartLoc);
XMDF_API int      xfpContinousToStridedDouble(double *a_OldArray, double** a_NewArray,
                                      const int a_NumValues, const int a_Stride,
                                      const int a_NumPerStride, const int a_StartLoc);
XMDF_API int*     xfpStridedToContinousInt(const int *a_OldArray,
                                   const int a_NumValues, const int a_Stride,
                                   const int a_NumPerStride, const int a_StartLoc);
XMDF_API int      xfpContinousToStridedInt(int *a_OldArray, int** a_NewArray,
                                   const int a_NumValues, const int a_Stride,
                                   const int a_NumPerStride, const int a_StartLoc);

#ifdef __cplusplus
}
#endif

#endif
