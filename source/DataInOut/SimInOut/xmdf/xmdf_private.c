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

#include "Xmdf.h"
#include "xmdf_private.h"
#include "ErrorDefinitions.h"
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef EMRL
#include <windows.h>
#endif

#ifdef PCAPI
#include <windows.h>
#endif

/* File globals */
/* Error messages */
static int      fg_nErrors = 0;
static char   **fg_Errors = NULL;
static hid_t    fg_Doubletype = 0;
static hid_t    fg_Floattype = 0;
static int      fg_initialized = XFALSE;
static xmbool   fg_running_parallel=XFALSE;
static int      fg_parallel_rank;
static hssize_t fg_parallel_numdatavalues=2000000;

/* WARNING: This file is intended to be used only by EMRL emplyees. */
/* ---------------------------------------------------------------------------
FUNCTION  xfpGetNumErrorMessages
PURPOSE   Get the number of messages on the error stack
NOTES     
----------------------------------------------------------------------------*/
XMDF_API int xfpGetNumErrorMessages (int *a_Num)
{
  *a_Num = fg_nErrors;
  return fg_nErrors;
} /* xfpGetNumErrorMessages */
/* ---------------------------------------------------------------------------
FUNCTION  xfpGetErrorMessages
PURPOSE   
NOTES     
----------------------------------------------------------------------------*/
XMDF_API int xfpGetErrorMessages (int a_Num, char **a_Errors)
{
  int    i = 0;

  /* return in error if the number isn't correct or if the file global is NULL */
  if (a_Num != fg_nErrors || fg_Errors == NULL) {
    return -1;
  }

  for (i = 0; i < a_Num; i++) {
    if (a_Errors[i]) {
      strcpy(a_Errors[i], fg_Errors[i]);
    }
  }
  return 1;
} /* xfpGetErrorMessages */
/* ---------------------------------------------------------------------------
FUNCTION  xfpAllocateErrorStack
PURPOSE   
NOTES     
----------------------------------------------------------------------------*/
XMDF_API void xfpAllocateErrorStack ()
{
  int  i = 0;

  if (fg_Errors == NULL) {
    /* allocate space for the maximum number of errors at the max size. */
    fg_Errors = (char**)malloc(XF_MAX_ERROR_STACK*sizeof(char*));
    for (i = 0; i < XF_MAX_ERROR_STACK; i++) {
      fg_Errors[i] = (char*)malloc(XF_MAX_ERROR_MSG_SIZE*sizeof(char));
    }
    fg_nErrors = 0;
  }
} /* xfpAllocateErrorStack */
/* ---------------------------------------------------------------------------
FUNCTION  xfpAddXMDFError
PURPOSE   
NOTES     
----------------------------------------------------------------------------*/
XMDF_API void xfpAddXMDFError (const char* a_Error)
{
  if (fg_Errors == NULL) {
    xfpAllocateErrorStack();
  }
  if (fg_nErrors < XF_MAX_ERROR_STACK) {
    strcpy(fg_Errors[fg_nErrors], a_Error);
    fg_nErrors++;
  }
} /* xfpAddXMDFError */
/* ---------------------------------------------------------------------------
FUNCTION  xfpClearErrors
PURPOSE   Clear the memory used by the error stack and reset it to zero
NOTES     
----------------------------------------------------------------------------*/
XMDF_API void xfpClearErrors ()
{
  int i;

  if (fg_Errors) {
    for (i = 0; i < XF_MAX_ERROR_STACK; i++) {
      free(fg_Errors[i]);
    }
    free(fg_Errors);
    fg_Errors = NULL;
    fg_nErrors = 0;
  }
} /* xfpClearErrors */
/* ---------------------------------------------------------------------------
FUNCTION  xfpHDF5ErrorHandler
PURPOSE   
NOTES     
----------------------------------------------------------------------------*/
XMDF_API herr_t xfpHDF5ErrorHandler (void* client_data)
{
#ifndef _WIN64
  H5Ewalk(H5E_WALK_UPWARD, xfpHDF5ErrorWalk_cb, NULL);
  return 0;
#else
  return -1;
#endif
} /* xfpHDF5ErrorHandler */
/* ---------------------------------------------------------------------------
FUNCTION  xfpHDF5ErrorWalk_cb
PURPOSE   
NOTES     
----------------------------------------------------------------------------*/
XMDF_API herr_t xfpHDF5ErrorWalk_cb (int n, H5E_error_t *err_desc, 
                                     void *client_data)
{
#ifndef _WIN64
  const char    *maj_str = NULL;
  const char    *min_str = NULL;
  const int       indent = 2;
  char           error[XF_MAX_ERROR_MSG_SIZE];

  /* Get descriptions for the major and minor error numbers */
  maj_str = H5Eget_major (err_desc->maj_num);
  min_str = H5Eget_minor (err_desc->min_num);

  /* Print error message */
  sprintf (error, "%*s#%03d: %s line %u in %s(): %s\n",
       indent, "", n, err_desc->file_name, err_desc->line,
       err_desc->func_name, err_desc->desc);
  xfpAddXMDFError(error);

  sprintf (error, "%*smajor(%02d): %s\n",
       indent*2, "", err_desc->maj_num, maj_str);
  xfpAddXMDFError(error);

  sprintf (error, "%*sminor(%02d): %s\n",
       indent*2, "", err_desc->min_num, min_str);
  xfpAddXMDFError(error);

  return 0;
#else
  return -1;
#endif
} /* xfpHDF5ErrorWalk_cb */
/******************************************************************************
 * FUNCTION  xfpCloseOpenIdentifiers
 * PURPOSE   Closes all open identifiers for a file 
 * NOTES     
 ******************2***********************************************************/
XMDF_API int xfpCloseOpenIdentifiers (xid a_Id)
{
#ifndef _WIN64
  int             i, Num, status;
  IdentifierInfo *Info;

  status = xfGetNumOpenIdentifiers(a_Id, &Num);
  if (status < 0) {
    /* HDF error */
    return -2;
  }

  if (Num > 0) {
    Info = (IdentifierInfo *) malloc(Num*sizeof(IdentifierInfo));
  
    status = xfGetOpenIdentifiersInfo(a_Id, Num, Info);
    if (status < 0) {
      free(Info);
      return status;
    }

    /* loop through and close open identifiers */
    for (i = 0; i < Num; i++) {
      switch(Info[i].type) {
        case (H5I_FILE):
          H5Fclose(Info[i].id);
          break;
        case (H5I_GROUP):
          H5Gclose(Info[i].id);         
          break;
        case (H5I_DATATYPE):
          H5Tclose(Info[i].id);
          break;
        case (H5I_DATASPACE):
          H5Sclose(Info[i].id);
          break;
        case (H5I_DATASET):
          H5Dclose(Info[i].id);
          break;
        case (H5I_ATTR):
          H5Aclose(Info[i].id);
          break;
        default:
          break;
      }
    }
   
    free(Info);
  }

  /* clear the error stack */
  xfpClearErrors();
  
  /* function successful */
  return 1;
#else
  return -1;
#endif
} /* xfpCloseOpenIdentifiers */
/* ---------------------------------------------------------------------------
FUNCTION  xfpClearGroup
PURPOSE   
NOTES     
----------------------------------------------------------------------------*/
XMDF_API int xfpClearGroup (xid a_Group)
{
#ifndef _WIN64
  herr_t  status;
  hsize_t nObjs;
  ssize_t iNameLength;
  int     i;
  char    *Path = NULL;

  /* open all groups under this group */
  status = H5Gget_num_objs(a_Group, &nObjs);
  if (status < 0) {
    return status;
  }
  /* loop through all the objects in the group */
  for (i = 0; i < nObjs; i++) {
/*    iNameLength = H5Gget_objname_by_idx(a_Group, i, GarbageBuffer, 0);*/
    iNameLength = H5Gget_objname_by_idx(a_Group, i, NULL, 0);
    if (iNameLength < 0) {
      return iNameLength;
    }
    Path = (char *)malloc(iNameLength+1*sizeof(char));
    iNameLength = H5Gget_objname_by_idx(a_Group, i, Path, iNameLength+1);
    if (iNameLength < 0) {
      return iNameLength;
    }
    /* Add NULL terminator */
    Path[iNameLength] = '\0';

    /* unlink the item */
    H5Gunlink(a_Group, Path);

    free(Path);
    Path = NULL;
  }
  return 1;
#else
  return -1;
#endif
} /* xfpClearGroup */
/*-----------------------------------------------------------------------------*/
/* FUNCTION xfpGetInitialized*/
/* PURPOSE  Returns whether or not the file globals fg_Floattype and */
/*          fg_Doubletype have been initialized*/
/* NOTES*/
/*-----------------------------------------------------------------------------*/
XMDF_API int xfpGetInitialized ()
{
  return fg_initialized;
} /* xfpInitialize*/

/*-----------------------------------------------------------------------------*/
/* FUNCTION xfpInitialize                                                      */
/* PURPOSE  Defaults the file globals fg_Floattype and fg_Doubletype           */
/* NOTES    Called anytime we create or open a file                            */
/*          **Most PC's use Little Endian, while UNIX uses Big Endian          */
/*-----------------------------------------------------------------------------*/
XMDF_API void xfpInitialize ()
{
#ifndef _WIN64
  // Initialize to little endian (PC's use Little Endian)
  fg_Floattype = H5T_IEEE_F32LE;
  fg_Doubletype = H5T_IEEE_F64LE;
  fg_initialized = XTRUE;
#endif
} /* xfpInitialize*/
/*-----------------------------------------------------------------------------*/
/* FUNCTION xfpDefaultWriteFloatType*/
/* PURPOSE  Sets either big endian or little endian float and double types*/
/* NOTES    Pass in "1" for big endian and "0" for little endian*/
/*-----------------------------------------------------------------------------*/
XMDF_API void xfpDefaultWriteFloatType (int a_BigEndian)
{
#ifndef _WIN64
  if (a_BigEndian) {
    fg_Floattype = H5T_IEEE_F32BE;
    fg_Doubletype = H5T_IEEE_F64BE;
  }
  else {
    fg_Floattype = H5T_IEEE_F32LE;
    fg_Doubletype = H5T_IEEE_F64LE;
  }
  fg_initialized = XTRUE;
#endif
} /* xfpDefaultWriteFloatType*/
/*-----------------------------------------------------------------------------*/
/* FUNCTION xfpGetDefaultFloatType*/
/* PURPOSE  returns the current float type (either big endian */
/*          or little endian)*/
/* NOTES*/
/*-----------------------------------------------------------------------------*/
XMDF_API hid_t xfpGetDefaultFloatType ()
{
  return fg_Floattype;
} /* xfpGetDefaultFloatType*/
/*-----------------------------------------------------------------------------*/
/* FUNCTION xfpGetDefaultDoubleType*/
/* PURPOSE  returns the current double type (either big endian */
/*          or little endian)*/
/* NOTES*/
/*-----------------------------------------------------------------------------*/
XMDF_API hid_t xfpGetDefaultDoubleType ()
{
  return fg_Doubletype;
} /* xfpGetDefaultDoubleType*/
/******************************************************************************
 * FUNCTION  xfpReadDset1D
 * PURPOSE   Reads a 1D dataset 
 * NOTES     the variable a_Values must already be allocated to a_Number 
 ******************2***********************************************************/
XMDF_API int xfpReadDset1D (xid a_Id, const char *a_Name, int a_Number,
              H5T_class_t a_DsetClass, xid a_MemType, void *a_Properties)
{
#ifndef _WIN64
  hid_t   DsetId, DataspaceId, DataTypeId;
  herr_t  status = 1;
  htri_t  IsSimple;
  int     Rank;
  hsize_t *Dims = NULL, *MaxDims = NULL;
  H5T_class_t DatasetClass;

  DsetId = H5Dopen(a_Id, a_Name);
  if (DsetId < 0) {
   return DsetId;
  }

  /* open the datatype and make sure it is the correct type */
  DataTypeId = H5Dget_type(DsetId);
  if (DataTypeId < 0) {
    return DataTypeId;
  }

  DatasetClass = H5Tget_class(DataTypeId);
  if (DatasetClass != a_DsetClass) {
    H5Tclose(DataTypeId);
    return ERROR_INCORRECT_DATATYPE;
  }


  /* Get the dataspace */
  DataspaceId = H5Dget_space(DsetId);

  /* Dataspace must be simply defined */
  IsSimple = H5Sis_simple(DataspaceId);
  if (!IsSimple) {
    H5Sclose(DataspaceId);
    H5Dclose(DsetId);
    H5Tclose(DataTypeId);
    return ERROR_DATASET_INVALID;
  }

  /* Get the rank and dimensions for simple dataspace */
  status = xfpGetSimpleDataspaceInfo(DataspaceId, &Rank, &Dims, &MaxDims);
  if (status < 0 || Rank != 1) {
    H5Sclose(DataspaceId);
    H5Aclose(DsetId);
    H5Tclose(DataTypeId);
    return ERROR_DATASET_INVALID;
  }

  /* Make sure the number of strings matches */
  if (a_Number != Dims[0]) {
    H5Sclose(DataspaceId);
    H5Aclose(DsetId);
    H5Tclose(DataTypeId);
    return ERROR_DATASET_INVALID;
  }

  status = H5Dread(DsetId, a_MemType, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   a_Properties);

  H5Tclose(DataTypeId);
  H5Sclose(DataspaceId);
  H5Dclose(DsetId);
 
  free(Dims);
  free(MaxDims);
  return status;
#else
  return -1;
#endif
} /* xfpReadDset1D */
/******************************************************************************
 * FUNCTION  xfpReadDset1DDouble
 * PURPOSE   Read a 1D double array from an HDF5 dataset
 * NOTES     returns positive if successful, otherwise negative
 ****************2***********************************************************/
XMDF_API int xfpReadDset1DDouble (xid a_Id, const char *a_Name, int a_Number, 
                                 double *a_Array)
{
#ifndef _WIN64
  herr_t    status;
  xid       DatasetId, SpaceId;
  int       Rank;
  hsize_t   *Dims = NULL, *Maxdims = NULL;

  DatasetId = H5Dopen(a_Id, a_Name);
  if (DatasetId < 0) {
    return ERROR_DATASET_INVALID;
  }

  /* make sure the size is right */
  SpaceId = H5Dget_space(DatasetId);
  if (SpaceId < 0) {
    H5Dclose(DatasetId);
    H5Sclose(SpaceId);
    return ERROR_DATASET_INVALID;
  }

  /* Make sure that the dataspace is simple */
  if (H5Sis_simple(SpaceId) <= 0) {
    H5Dclose(DatasetId);
    H5Sclose(SpaceId);
    return ERROR_DATASET_INVALID;
  }
  status = xfpGetSimpleDataspaceInfo(SpaceId, &Rank, &Dims, &Maxdims);

  if (status < 0 || Rank != 1 || Dims[0] != a_Number) {
    H5Sclose(SpaceId);
    H5Dclose(DatasetId);
    return ERROR_DATASET_INVALID;
  }
  
  if (Dims) {
    free(Dims);
  }
  if (Maxdims) {
    free(Maxdims);
  }

  status = H5Dread(DatasetId, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, 
                   H5P_DEFAULT, a_Array);

  H5Dclose(DatasetId);
 
  return status;
#else
  return -1;
#endif
} /* xfpReadDset1DDouble */
/******************************************************************************
 * FUNCTION  xfpReadDset1DFloat
 * PURPOSE   Read a 1D float array from an HDF5 dataset
 * NOTES     returns positive if successful, otherwise negative
 ****************2***********************************************************/
XMDF_API int xfpReadDset1DFloat (xid a_Id, const char *a_Name, int a_Number, float *a_Array)
{
#ifndef _WIN64
  herr_t    status;
  xid       DatasetId, SpaceId;
  int       Rank;
  hsize_t   *Dims = NULL, *Maxdims = NULL;

  DatasetId = H5Dopen(a_Id, a_Name);
  if (DatasetId < 0) {
    return ERROR_DATASET_INVALID;
  }

  /* make sure the size is right */
  SpaceId = H5Dget_space(DatasetId);
  if (SpaceId < 0) {
    H5Dclose(DatasetId);
    H5Sclose(SpaceId);
    return ERROR_DATASET_INVALID;
  }

  /* Make sure that the dataspace is simple */
  if (H5Sis_simple(SpaceId) <= 0) {
    H5Dclose(DatasetId);
    H5Sclose(SpaceId);
    return ERROR_DATASET_INVALID;
  }
  status = xfpGetSimpleDataspaceInfo(SpaceId, &Rank, &Dims, &Maxdims);

  if (status < 0 || Rank != 1 || Dims[0] != a_Number) {
    H5Sclose(SpaceId);
    H5Dclose(DatasetId);
    return ERROR_DATASET_INVALID;
  }
  
  if (Dims) {
    free(Dims);
  }
  if (Maxdims) {
    free(Maxdims);
  }

  status = H5Dread(DatasetId, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, 
                   H5P_DEFAULT, a_Array);

  H5Dclose(DatasetId);
 
  return status;
#else
  return -1;
#endif
} /* xfpReadDset1DFloat */
/******************************************************************************
 * FUNCTION  xfpReadDsetDouble
 * PURPOSE   Read a double array from an HDF5 dataset
 * NOTES     Data must already be allocated
 ****************2***********************************************************/
XMDF_API int xfpReadDsetDouble (xid a_Id, const char *a_Name, 
                              int a_Rank, const hsize_t *a_Dims, double *a_Array)
{
#ifndef _WIN64
  herr_t    status;
  xid       DatasetId, SpaceId, DatatypeId, IndicesId;
  int       Rank, i;
  hsize_t   *Dims = NULL, *Maxdims = NULL;
  H5T_class_t Classtype;
  double    *Array, *Indices;

  DatasetId = H5Dopen(a_Id, a_Name);
  if (DatasetId < 0) {
    return ERROR_DATASET_INVALID;
  }

  /* Make sure the datatype is a float datatype */
  DatatypeId = H5Dget_type(DatasetId);
  if (DatatypeId < 0) {
    H5Dclose(DatasetId);
    return DatatypeId;
  }
  Classtype = H5Tget_class(DatatypeId);
  H5Tclose(DatatypeId);
  if (Classtype != H5T_FLOAT) {
    H5Dclose(DatasetId);
    return ERROR_DATASET_INVALID;
  }

  SpaceId = H5Dget_space(DatasetId);
  if (SpaceId < 0) {
    H5Dclose(DatasetId);
    H5Sclose(SpaceId);
    return ERROR_DATASET_INVALID;
  }

  /* Make sure that the dataspace is simple */
  if (H5Sis_simple(SpaceId) <= 0) {
    H5Dclose(DatasetId);
    H5Sclose(SpaceId);
    return ERROR_DATASET_INVALID;
  }

  /* Make sure the rank is 2 */
  status = xfpGetSimpleDataspaceInfo(SpaceId, &Rank, &Dims, &Maxdims);

  if (status < 0 || Rank != a_Rank) {
    if (Dims) {
      free(Dims);
    }
    if (Maxdims) {
      free(Maxdims);
    }
    H5Sclose(SpaceId);
    H5Dclose(DatasetId);
    return ERROR_DATASET_INVALID;
  }
  // Make sure the sizes match up
  for (i = 0; i < Rank; i++) {
    if (a_Dims[i] != Dims[i]) {
      if (Dims) {
        free(Dims);
      }
      if (Maxdims) {
        free(Maxdims);
      }
      H5Sclose(SpaceId);
      H5Dclose(DatasetId);
      return ERROR_DATASET_INVALID;
    }
  }


  /* See if the Dataset was written out in parallel */
  /* Because if it was then it was written differently and to get the data we */
  /* do something different. We know it was written out in parallel if the indices */
  /* array exists. */
  IndicesId = H5Dopen(a_Id, DATASET_DSET_INDICES);
  if (IndicesId > -1) {

    /* create array to get the data */
    Array = (double*) malloc((size_t)Dims[0]*sizeof(double));

    /* Read the data */
    status = H5Dread(DatasetId, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                    H5P_DEFAULT, Array);

    if (Dims) {
      free(Dims);
    }
    if (Maxdims) {
      free(Maxdims);
    }

    /* read the indices array to read the values in correctly */
    SpaceId = H5Dget_space(IndicesId);
    if (SpaceId < 0) {
      H5Dclose(DatasetId);
      H5Sclose(SpaceId);
      H5Sclose(IndicesId);
      return ERROR_DATASET_INVALID;
    }
    /* Get the Dimensions */
    status = xfpGetSimpleDataspaceInfo(SpaceId, &Rank, &Dims, &Maxdims);

    /* create array to get the data */
    Indices = (double*) malloc((size_t)Dims[0]*sizeof(double));
    /* Read the Index data */
    status = H5Dread(IndicesId, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                     H5P_DEFAULT, Indices);

    /* Organize the data by the index array */
    for (i=0; i<Dims[0]; i=i+1) {
      if (Indices[i] != -1) {
        a_Array[(int)Indices[i]] = Array[i];
                
      }
    }

    H5Dclose(IndicesId);
    free(Indices);
  }
  else {
    /* Read the data */
    status = H5Dread(DatasetId, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                    H5P_DEFAULT, a_Array);
  }

  if (Dims) {
    free(Dims);
  }
  if (Maxdims) {
    free(Maxdims);
  }

  H5Sclose(SpaceId);
  H5Dclose(DatasetId);
 
  return status;
#else
  return -1;
#endif
} /* xfpReadDsetDouble */
/******************************************************************************
 * FUNCTION  xfpReadDsetInt
 * PURPOSE   Read a  int array from an HDF5 dataset
 * NOTES     Data must already be allocated
 ****************2***********************************************************/
XMDF_API int xfpReadDsetInt (xid a_Id, const char *a_Name, 
                              int a_Rank, const hsize_t *a_Dims, int *a_Array)
{
#ifndef _WIN64
  herr_t    status;
  xid       DatasetId, SpaceId, DatatypeId;
  int       Rank, i;
  hsize_t   *Dims = NULL, *Maxdims = NULL;
  H5T_class_t Classtype;

  DatasetId = H5Dopen(a_Id, a_Name);
  if (DatasetId < 0) {
    return ERROR_DATASET_INVALID;
  }

  /* Make sure the datatype is an integer datatype */
  DatatypeId = H5Dget_type(DatasetId);
  if (DatatypeId < 0) {
    H5Dclose(DatasetId);
    return DatatypeId;
  }
  Classtype = H5Tget_class(DatatypeId);
  H5Tclose(DatatypeId);
  if (Classtype != H5T_INTEGER) {
    H5Dclose(DatasetId);
    return ERROR_DATASET_INVALID;
  }

  SpaceId = H5Dget_space(DatasetId);
  if (SpaceId < 0) {
    H5Dclose(DatasetId);
    H5Sclose(SpaceId);
    return ERROR_DATASET_INVALID;
  }

  /* Make sure that the dataspace is simple */
  if (H5Sis_simple(SpaceId) <= 0) {
    H5Dclose(DatasetId);
    H5Sclose(SpaceId);
    return ERROR_DATASET_INVALID;
  }

  /* Make sure the rank is 2 */
  status = xfpGetSimpleDataspaceInfo(SpaceId, &Rank, &Dims, &Maxdims);

  if (status < 0 || Rank != a_Rank) {
    if (Dims) {
      free(Dims);
    }
    if (Maxdims) {
      free(Maxdims);
    }
    H5Sclose(SpaceId);
    H5Dclose(DatasetId);
    return ERROR_DATASET_INVALID;
  }
  // Make sure the sizes match up
  for (i = 0; i < Rank; i++) {
    if (a_Dims[i] != Dims[i]) {
      if (Dims) {
        free(Dims);
      }
      if (Maxdims) {
        free(Maxdims);
      }
      H5Sclose(SpaceId);
      H5Dclose(DatasetId);
      return ERROR_DATASET_INVALID;
    }
  }
  if (Dims) {
    free(Dims);
  }
  if (Maxdims) {
    free(Maxdims);
  }

  /* Read the data */
  status = H5Dread(DatasetId, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
                   H5P_DEFAULT, a_Array);

  H5Sclose(SpaceId);
  H5Dclose(DatasetId);
 
  return status;
#else
  return -1;
#endif
} /* xfpReadDsetInt */
/******************************************************************************
 * FUNCTION  xfpReadDset1DIntPortion
 * PURPOSE   Read a portion of a 1D int array from an HDF5 dataset
 * NOTES     returns positive if successful, otherwise negative
 *           a_Start1 is the starting index for the first array index, 
 *           a_Num1 is number to read in first array index
 *           a_Start2 and a_Num2 are same as a_Start1 and a_Num1 but for index2
 *           Data is always stored in a 1D array
 *           Data must already be allocated of size a_Num1 * a_Num2
 ****************2***********************************************************/
XMDF_API int xfpReadDset1DIntPortion (xid a_Id, const char *a_Name, 
                hssize_t a_Start1, hsize_t a_Num1, int *a_Array)
{
#ifndef _WIN64
  herr_t    status;
  xid       DatasetId, SpaceId, MemspaceId, DatatypeId;
  int       Rank;
  hsize_t   Count[1], *Dims = NULL, *Maxdims = NULL;
  hssize_t  Offset[1], Memsize;
  H5T_class_t Classtype;

  DatasetId = H5Dopen(a_Id, a_Name);
  if (DatasetId < 0) {
    return ERROR_DATASET_INVALID;
  }

  /* Make sure the datatype is a float datatype */
  DatatypeId = H5Dget_type(DatasetId);
  if (DatatypeId < 0) {
    H5Dclose(DatasetId);
    return DatatypeId;
  }
  Classtype = H5Tget_class(DatatypeId);
  H5Tclose(DatatypeId);
  if (Classtype != H5T_INTEGER) {
    H5Dclose(DatasetId);
    return ERROR_DATASET_INVALID;
  }

  SpaceId = H5Dget_space(DatasetId);
  if (SpaceId < 0) {
    H5Dclose(DatasetId);
    H5Sclose(SpaceId);
    return ERROR_DATASET_INVALID;
  }

  /* Make sure that the dataspace is simple */
  if (H5Sis_simple(SpaceId) <= 0) {
    H5Dclose(DatasetId);
    H5Sclose(SpaceId);
    return ERROR_DATASET_INVALID;
  }

  /* Make sure the rank is 2 */
  status = xfpGetSimpleDataspaceInfo(SpaceId, &Rank, &Dims, &Maxdims);
  if (Dims) {
    free(Dims);
  }
  if (Maxdims) {
    free(Maxdims);
  }
  if (status < 0 || Rank != 1) {
    H5Sclose(SpaceId);
    H5Dclose(DatasetId);
    return ERROR_DATASET_INVALID;
  }

  /* select the hyperslab */
  Count[0] = a_Num1;

  Offset[0] = a_Start1;

  status = H5Sselect_hyperslab(SpaceId, H5S_SELECT_SET, Offset, NULL, Count, 
                               NULL);
  if (status < 0) {
    H5Sclose(SpaceId);
    H5Dclose(DatasetId);
    return status;
  }

  /* Make sure the portion is valid */
  if (H5Sselect_valid(SpaceId) <= 0) {
    return ERROR_DATASET_SIZE_INCORRECT;
  }

  /* Create a memory dataspace for the data */
  Memsize = H5Sget_select_npoints(SpaceId);
  MemspaceId = H5Screate_simple(1, &Memsize, NULL);
  if (MemspaceId < 0) {
    H5Sclose(SpaceId);
    H5Dclose(DatasetId);
    return MemspaceId;
  }

  /* Read the data */
  status = H5Dread(DatasetId, H5T_NATIVE_INT, MemspaceId, SpaceId,
                   H5P_DEFAULT, a_Array);

  H5Sclose(SpaceId);
  H5Dclose(DatasetId);
  H5Sclose(MemspaceId);
 
  return status;
#else
  return -1;
#endif
} /* xfpReadDset2DIntPortion */
/******************************************************************************
 * FUNCTION  xfpReadDset1DDoublePortion
 * PURPOSE   Read a portion of a 1D double array from an HDF5 dataset
 * NOTES     returns positive if successful, otherwise negative
 *           a_Start1 is the starting index for the first array index, 
 *           a_Num1 is number to read in first array index
 *           a_Start2 and a_Num2 are same as a_Start1 and a_Num1 but for index2
 *           Data is always stored in a 1D array
 *           Data must already be allocated of size a_Num1 * a_Num2
 ****************2***********************************************************/
XMDF_API int xfpReadDset1DDoublePortion (xid a_Id, const char *a_Name, 
                hssize_t a_Start1, hsize_t a_Num1, double *a_Array)
{
#ifndef _WIN64
  herr_t    status;
  xid       DatasetId, SpaceId, MemspaceId, DatatypeId;
  int       Rank;
  hsize_t   Count[1], *Dims = NULL, *Maxdims = NULL;
  hssize_t  Offset[1], Memsize;
  H5T_class_t Classtype;

  DatasetId = H5Dopen(a_Id, a_Name);
  if (DatasetId < 0) {
    return ERROR_DATASET_INVALID;
  }

  /* Make sure the datatype is a float datatype */
  DatatypeId = H5Dget_type(DatasetId);
  if (DatatypeId < 0) {
    H5Dclose(DatasetId);
    return DatatypeId;
  }
  Classtype = H5Tget_class(DatatypeId);
  H5Tclose(DatatypeId);
  if (Classtype != H5T_FLOAT) {
    H5Dclose(DatasetId);
    return ERROR_DATASET_INVALID;
  }

  SpaceId = H5Dget_space(DatasetId);
  if (SpaceId < 0) {
    H5Dclose(DatasetId);
    H5Sclose(SpaceId);
    return ERROR_DATASET_INVALID;
  }

  /* Make sure that the dataspace is simple */
  if (H5Sis_simple(SpaceId) <= 0) {
    H5Dclose(DatasetId);
    H5Sclose(SpaceId);
    return ERROR_DATASET_INVALID;
  }

  /* Make sure the rank is 2 */
  status = xfpGetSimpleDataspaceInfo(SpaceId, &Rank, &Dims, &Maxdims);
  if (Dims) {
    free(Dims);
  }
  if (Maxdims) {
    free(Maxdims);
  }
  if (status < 0 || Rank != 1) {
    H5Sclose(SpaceId);
    H5Dclose(DatasetId);
    return ERROR_DATASET_INVALID;
  }

  /* select the hyperslab */
  Count[0] = a_Num1;

  Offset[0] = a_Start1;

  status = H5Sselect_hyperslab(SpaceId, H5S_SELECT_SET, Offset, NULL, Count, 
                               NULL);
  if (status < 0) {
    H5Sclose(SpaceId);
    H5Dclose(DatasetId);
    return status;
  }

  /* Make sure the portion is valid */
  if (H5Sselect_valid(SpaceId) <= 0) {
    return ERROR_DATASET_SIZE_INCORRECT;
  }

  /* Create a memory dataspace for the data */
  Memsize = H5Sget_select_npoints(SpaceId);
  MemspaceId = H5Screate_simple(1, &Memsize, NULL);
  if (MemspaceId < 0) {
    H5Sclose(SpaceId);
    H5Dclose(DatasetId);
    return MemspaceId;
  }

  /* Read the data */
  status = H5Dread(DatasetId, H5T_NATIVE_DOUBLE, MemspaceId, SpaceId,
                   H5P_DEFAULT, a_Array);

  H5Sclose(SpaceId);
  H5Dclose(DatasetId);
  H5Sclose(MemspaceId);
 
  return status;
#else
  return -1;
#endif
} /* xfpReadDset1DDoublePortion */
/******************************************************************************
 * FUNCTION  xfpWriteDset1DFloatPortion
 * PURPOSE   Write a portion of a 1D float array to an HDF5 dataset
 * NOTES     returns positive if successful, otherwise negative
 *           a_Start1 is the starting index for the first array index, 
 *           a_Num1 is number to write in first array index
 *           Data is always stored in a 1D array
 *           Data must already be allocated of size a_Num1
 ****************2***********************************************************/
int xfpWriteDset1DFloatPortion (xid a_Id, const char *a_Name, int a_Start1,
                               int a_Num1, const float *a_Array)
{
#ifndef _WIN64
  herr_t     status = XTRUE;
  hsize_t   *OldDims = NULL, *OldMaxDims = NULL, NewDims[2], MemDims[2];
  H5T_class_t Classtype;
  xid        DatasetId, DataspaceId, DatatypeId, FilespaceId, MemspaceId;
  int        Rank;
  hssize_t   offset[2];
  hsize_t    count[2];

  DatasetId = DataspaceId = DatatypeId = FilespaceId = MemspaceId = NONE;
  DatasetId = H5Dopen(a_Id, a_Name);
  if (DatasetId < 0) {
    return ERROR_DATASET_INVALID;
  }
  DataspaceId = H5Dget_space(DatasetId);
  if (DataspaceId < 0) {
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    return ERROR_DATASET_INVALID;
  }

  /* Dataspace must be simple and have a rank of 1 */
  if (H5Sis_simple(DataspaceId) <= 0) {
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    return ERROR_DATASET_INVALID;
  }

  /* Make sure the datatype is a float datatype */
  DatatypeId = H5Dget_type(DatasetId);
  if (DatatypeId < 0) {
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    return ERROR_DATASET_INVALID;
  }
  Classtype = H5Tget_class(DatatypeId);
  H5Tclose(DatatypeId);
  if (Classtype != H5T_FLOAT) {
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    return ERROR_DATASET_INVALID;
  }

  status = xfpGetSimpleDataspaceInfo(DataspaceId, &Rank, &OldDims, &OldMaxDims);
  if (status < 0) {
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    return ERROR_DATASET_INVALID;
  }

  /* This function is only for 1D arrays and MaxDims must be larger */
  /* than the old dimension plus the new values */
  if (Rank != 1) {
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    return ERROR_DATASET_INVALID;
  }
  /* see if we need to extend the dataset */
  if (OldDims[0] < (a_Start1 + a_Num1)) {
    NewDims[0] = a_Start1 + a_Num1;
    if (OldMaxDims[0] < NewDims[0]) {
      H5Dclose(DatasetId);
      H5Sclose(DataspaceId);
      return ERROR_DATASET_INVALID;
    }

    /* extend the dataset */
    status = H5Dextend(DatasetId, NewDims);
    if (status < 0) {
      H5Dclose(DatasetId);
      H5Sclose(DataspaceId);
      return ERROR_DATASET_INVALID;
    }
  }

  /* select a hyperslab */
  FilespaceId = H5Dget_space(DatasetId);
  if (FilespaceId < 0) {
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    return ERROR_DATASET_INVALID;
  }

  /* create the memory dataspace */
  MemDims[0] = a_Num1;
  MemspaceId = H5Screate_simple(1, MemDims, NULL);
  if (MemspaceId < 0) {
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    return ERROR_DATASET_INVALID;
  }

  /* add new information */
  offset[0] = a_Start1;
  count[0] = MemDims[0];

  status = H5Sselect_hyperslab(FilespaceId, H5S_SELECT_SET, offset, NULL, 
                                count, NULL);
  if (status < 0) {
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    return ERROR_DATASET_INVALID;
  }

  /* write the data */
  status = H5Dwrite(DatasetId, H5T_NATIVE_FLOAT, MemspaceId,
                    FilespaceId, H5P_DEFAULT, a_Array);
  if (status < 0) {
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    return ERROR_DATASET_INVALID;
  }

  /* clean up */
  if (OldDims)
    free(OldDims);
  if (OldMaxDims)
    free(OldMaxDims);
  if (MemspaceId != NONE)
    H5Sclose(MemspaceId);
  if (DatasetId != NONE)
    H5Dclose(DatasetId);
  if (DataspaceId != NONE)
    H5Sclose(DataspaceId);
  if (FilespaceId != NONE)
    H5Sclose(FilespaceId);
  return status;
#else
  return -1;
#endif
} /* xwWriteDset1DFloatPortion */
/******************************************************************************
 * FUNCTION  xfpWriteDset2DFloatPortion
 * PURPOSE   Write a portion of a 2D float array to an HDF5 dataset
 * NOTES     returns positive if successful, otherwise negative
 *           a_Start1 is the starting index for the first array index, 
 *           a_Num1 is number to write in first array index
 *           a_Start2 and a_Num2 are same as a_Start1 and a_Num1 but for index2
 *           Data is always stored in a 1D array
 *           Data must already be allocated of size a_Num1 * a_Num2
 ****************2***********************************************************/
XMDF_API int xfpWriteDset2DFloatPortion (xid a_Id, const char *a_Name, int a_Start1,
                                int a_Num1, int a_Start2, int a_Num2,
                                const float *a_Array)
{
#ifndef _WIN64
  herr_t     status = XTRUE;
  hsize_t   *OldDims = NULL, *OldMaxDims = NULL, NewDims[2], MemDims[2];
  H5T_class_t Classtype;
  xid        DatasetId, DataspaceId, DatatypeId, FilespaceId, MemspaceId;
  int        Rank, i, extend;
  hssize_t   offset[2];
  hsize_t    count[2];
  hssize_t   nMemSelect;
  hssize_t   nSelect;

  DatasetId = DataspaceId = DatatypeId = FilespaceId = MemspaceId = NONE;
  DatasetId = H5Dopen(a_Id, a_Name);
  if (DatasetId < 0) {
    return ERROR_DATASET_INVALID;
  }

  DataspaceId = H5Dget_space(DatasetId);
  if (DataspaceId < 0) {
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    return ERROR_DATASET_INVALID;
  }

  /* Dataspace must be simple and have a rank of 2 */
  if (H5Sis_simple(DataspaceId) <= 0) {
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    return ERROR_DATASET_INVALID;
  }

  status = xfpGetSimpleDataspaceInfo(DataspaceId, &Rank, &OldDims, &OldMaxDims);
  if (status < 0) {
    H5Sclose(DataspaceId);
    H5Dclose(DatasetId);
    return ERROR_DATASET_INVALID;
  }

  /* This function is only for 3D arrays and MaxDims must be larger
     than the old dimension plus the new values in both directions */
  if (Rank != 2) {
    H5Sclose(DataspaceId);
    H5Dclose(DatasetId);
    return ERROR_DATASET_INVALID;
  }
  extend = XFALSE;
  for (i = 0; i < 2; i++)
    NewDims[i] = OldDims[i];
  if (OldDims[0] < (a_Start1 + a_Num1)) {
    NewDims[0] = a_Start1 + a_Num1;
    extend = XTRUE;
  }
  if (OldDims[1] < (a_Start2 + a_Num2)) {
    NewDims[1] = a_Start2 + a_Num2;
    extend = XTRUE;
  }
  if (OldMaxDims[0] < NewDims[0] || OldMaxDims[1] < NewDims[1]) {
    H5Sclose(DataspaceId);
    H5Dclose(DatasetId);
    return ERROR_DATASET_INVALID;
  }

  /* Make sure the datatype is a float datatype */
  DatatypeId = H5Dget_type(DatasetId);
  if (DatatypeId < 0) {
    H5Sclose(DataspaceId);
    H5Dclose(DatasetId);
    return ERROR_DATASET_INVALID;
  }
  Classtype = H5Tget_class(DatatypeId);
  H5Tclose(DatatypeId);
  if (Classtype != H5T_FLOAT) {
    H5Sclose(DataspaceId);
    H5Dclose(DatasetId);
    return ERROR_DATASET_INVALID;
  }

  /* extend the dataset */
  if (extend) {
    status = H5Dextend(DatasetId, NewDims);
    if (status < 0) {
      H5Sclose(DataspaceId);
      H5Dclose(DatasetId);
      return ERROR_DATASET_INVALID;
    }
  }

  /* select a hyperslab */
  FilespaceId = H5Dget_space(DatasetId);
  if (FilespaceId < 0) {
    H5Sclose(DataspaceId);
    H5Dclose(DatasetId);
    return ERROR_DATASET_INVALID;
  }

  /* create the memory dataspace */
  MemDims[0] = a_Num1;
  MemDims[1] = a_Num2;
  MemspaceId = H5Screate_simple(2, MemDims, NULL);
  if (MemspaceId < 0) {
    H5Sclose(DataspaceId);
    H5Dclose(DatasetId);
    return ERROR_DATASET_INVALID;
  }
  nMemSelect = H5Sget_select_npoints(MemspaceId);

  /* add extras */
  offset[0] = a_Start1;
  offset[1] = a_Start2;
  count[0] = a_Num1;
  count[1] = a_Num2;

  status = H5Sselect_hyperslab(FilespaceId, H5S_SELECT_SET, offset, NULL, 
                                count, NULL);
  if (status < 0) {
    H5Sclose(DataspaceId);
    H5Dclose(DatasetId);
    return ERROR_DATASET_INVALID;
  }
  nSelect = H5Sget_select_npoints(FilespaceId);

  /* write the data */
  status = H5Dwrite(DatasetId, H5T_NATIVE_FLOAT, MemspaceId,
                    FilespaceId, H5P_DEFAULT, a_Array);
  if (status < 0) {
    H5Sclose(DataspaceId);
    H5Dclose(DatasetId);
    return ERROR_DATASET_INVALID;
  }

  if (OldDims)
    free(OldDims);
  if (OldMaxDims)
    free(OldMaxDims);
  if (MemspaceId != NONE)
    H5Sclose(MemspaceId);
  if (DatasetId != NONE)
    H5Dclose(DatasetId);
  if (DataspaceId)
    H5Sclose(DataspaceId);
  if (FilespaceId)
    H5Sclose(FilespaceId);

  return (int)status;
#else
  return -1;
#endif
} /* xfpWriteDset2DFloatPortion */
/******************************************************************************
 * FUNCTION  xfpReadDset2DFloatPortion
 * PURPOSE   Read a portion of a 2D float array from an HDF5 dataset
 * NOTES     returns positive if successful, otherwise negative
 *           a_Start1 is the starting index for the first array index, 
 *           a_Num1 is number to read in first array index
 *           a_Start2 and a_Num2 are same as a_Start1 and a_Num1 but for index2
 *           Data is always stored in a 1D array
 *           Data must already be allocated of size a_Num1 * a_Num2
 ****************2***********************************************************/
XMDF_API int xfpReadDset2DFloatPortion (xid a_Id, const char *a_Name, 
                hssize_t a_Start1, hsize_t a_Num1,
                hssize_t a_Start2, hsize_t a_Num2, float *a_Array)
{
#ifndef _WIN64
  herr_t    status;
  xid       DatasetId, SpaceId, MemspaceId, DatatypeId;
  int       Rank;
  hsize_t   Count[2], *Dims = NULL, *Maxdims = NULL;
  hssize_t  Offset[2], Memsize;
  H5T_class_t Classtype;

  DatasetId = H5Dopen(a_Id, a_Name);
  if (DatasetId < 0) {
    return ERROR_DATASET_INVALID;
  }

  /* Make sure the datatype is a float datatype */
  DatatypeId = H5Dget_type(DatasetId);
  if (DatatypeId < 0) {
    H5Dclose(DatasetId);
    return DatatypeId;
  }
  Classtype = H5Tget_class(DatatypeId);
  H5Tclose(DatatypeId);
  if (Classtype != H5T_FLOAT) {
    H5Dclose(DatasetId);
    return ERROR_DATASET_INVALID;
  }

  SpaceId = H5Dget_space(DatasetId);
  if (SpaceId < 0) {
    H5Dclose(DatasetId);
    H5Sclose(SpaceId);
    return ERROR_DATASET_INVALID;
  }

  /* Make sure that the dataspace is simple */
  if (H5Sis_simple(SpaceId) <= 0) {
    H5Dclose(DatasetId);
    H5Sclose(SpaceId);
    return ERROR_DATASET_INVALID;
  }

  /* Make sure the rank is 2 */
  status = xfpGetSimpleDataspaceInfo(SpaceId, &Rank, &Dims, &Maxdims);
  if (Dims) {
    free(Dims);
  }
  if (Maxdims) {
    free(Maxdims);
  }
  if (status < 0 || Rank != 2) {
    H5Sclose(SpaceId);
    H5Dclose(DatasetId);
    return ERROR_DATASET_INVALID;
  }

  /* select the hyperslab */
  Count[0] = a_Num1;
  Count[1] = a_Num2;

  Offset[0] = a_Start1;
  Offset[1] = a_Start2;

  status = H5Sselect_hyperslab(SpaceId, H5S_SELECT_SET, Offset, NULL, Count, 
                               NULL);
  if (status < 0) {
    H5Sclose(SpaceId);
    H5Dclose(DatasetId);
    return status;
  }

  /* Make sure the portion is valid */
  if (H5Sselect_valid(SpaceId) <= 0) {
    return ERROR_DATASET_SIZE_INCORRECT;
  }

  /* Create a memory dataspace for the data */
  Memsize = H5Sget_select_npoints(SpaceId);
  MemspaceId = H5Screate_simple(1, &Memsize, NULL);
  if (MemspaceId < 0) {
    H5Sclose(SpaceId);
    H5Dclose(DatasetId);
    return MemspaceId;
  }

  /* Read the data */
  status = H5Dread(DatasetId, H5T_NATIVE_FLOAT, MemspaceId, SpaceId,
                   H5P_DEFAULT, a_Array);

  H5Sclose(SpaceId);
  H5Dclose(DatasetId);
  H5Sclose(MemspaceId);
 
  return status;
#else
  return -1;
#endif
} /* xfpReadDset2DFloatPortion */

/******************************************************************************
 * FUNCTION  xfpReadDset2DFloatPortions
 * PURPOSE   Read a portion of a 2D float array from an HDF5 dataset
 * NOTES     returns positive if successful, otherwise negative
 *           a_Start1 is the starting index for the first array index, 
 *           a_Num1 is number to read in first array index
 *           a_Start2 and a_Num2 are same as a_Start1 and a_Num1 but for index2
 *           Data is always stored in a 1D array
 *           Data must already be allocated of size a_Num1 * a_Num2
 ****************2***********************************************************/
XMDF_API int xfpReadDset2DFloatPortions (xid a_Id, const char *a_Name, 
                int a_nPortions, hssize_t *a_Start1, hsize_t *a_Num1,
                hssize_t *a_Start2, hsize_t *a_Num2, float *a_Array)
{
#ifndef _WIN64
  herr_t    status;
  xid       DatasetId, SpaceId, MemspaceId, DatatypeId;
  int       Rank, iPortion;
  hsize_t   Count[2], *Dims = NULL, *Maxdims = NULL;
  hssize_t  Offset[2], Memsize;
  H5T_class_t Classtype;

  DatasetId = H5Dopen(a_Id, a_Name);
  if (DatasetId < 0) {
    return ERROR_DATASET_INVALID;
  }

  /* Make sure the datatype is a float datatype */
  DatatypeId = H5Dget_type(DatasetId);
  if (DatatypeId < 0) {
    H5Dclose(DatasetId);
    return DatatypeId;
  }
  Classtype = H5Tget_class(DatatypeId);
  H5Tclose(DatatypeId);
  if (Classtype != H5T_FLOAT) {
    H5Dclose(DatasetId);
    return ERROR_DATASET_INVALID;
  }

  SpaceId = H5Dget_space(DatasetId);
  if (SpaceId < 0) {
    H5Dclose(DatasetId);
    H5Sclose(SpaceId);
    return ERROR_DATASET_INVALID;
  }

  /* Make sure that the dataspace is simple */
  if (H5Sis_simple(SpaceId) <= 0) {
    H5Dclose(DatasetId);
    H5Sclose(SpaceId);
    return ERROR_DATASET_INVALID;
  }

  /* Make sure the rank is 2 */
  status = xfpGetSimpleDataspaceInfo(SpaceId, &Rank, &Dims, &Maxdims);
  if (Dims) {
    free(Dims);
  }
  if (Maxdims) {
    free(Maxdims);
  }
  if (status < 0 || Rank != 2) {
    H5Sclose(SpaceId);
    H5Dclose(DatasetId);
    return ERROR_DATASET_INVALID;
  }

  H5Sselect_none(SpaceId);

  for (iPortion = 0; iPortion < a_nPortions; iPortion++) {
    /* select the hyperslabs */
    Count[0] = a_Num1[iPortion];
    Count[1] = a_Num2[iPortion];

    Offset[0] = a_Start1[iPortion];
    Offset[1] = a_Start2[iPortion];

    status = H5Sselect_hyperslab(SpaceId, H5S_SELECT_OR, Offset, NULL, Count, 
                                 NULL);
    if (status < 0) {
      H5Sclose(SpaceId);
      H5Dclose(DatasetId);
      return status;
    }
  }

  /* Make sure the portion is valid */
  if (H5Sselect_valid(SpaceId) <= 0) {
    return ERROR_DATASET_SIZE_INCORRECT;
  }

  /* Create a memory dataspace for the data */
  Memsize = H5Sget_select_npoints(SpaceId);
  MemspaceId = H5Screate_simple(1, &Memsize, NULL);
  if (MemspaceId < 0) {
    H5Sclose(SpaceId);
    H5Dclose(DatasetId);
    return MemspaceId;
  }

  /* Read the data */
  status = H5Dread(DatasetId, H5T_NATIVE_FLOAT, MemspaceId, SpaceId,
                   H5P_DEFAULT, a_Array);

  H5Sclose(SpaceId);
  H5Dclose(DatasetId);
  H5Sclose(MemspaceId);
 
  return status;
#else
  return -1;
#endif
} /* xfpReadDset2DFloatPortions */
/******************************************************************************
 * FUNCTION  xfpReadDset2DIntPortion
 * PURPOSE   Read a portion of a 2D int array from an HDF5 dataset
 * NOTES     returns positive if successful, otherwise negative
 *           a_Start1 is the starting index for the first array index, 
 *           a_Num1 is number to read in first array index
 *           a_Start2 and a_Num2 are same as a_Start1 and a_Num1 but for index2
 *           Data is always stored in a 1D array
 *           Data must already be allocated of size a_Num1 * a_Num2
 ****************2***********************************************************/
XMDF_API int xfpReadDset2DIntPortion (xid a_Id, const char *a_Name, 
                hssize_t a_Start1, hsize_t a_Num1,
                hssize_t a_Start2, hsize_t a_Num2, int *a_Array)
{
#ifndef _WIN64
  herr_t    status;
  xid       DatasetId, SpaceId, MemspaceId, DatatypeId;
  int       Rank;
  hsize_t   Count[2], *Dims = NULL, *Maxdims = NULL;
  hssize_t  Offset[2], Memsize;
  H5T_class_t Classtype;

  DatasetId = H5Dopen(a_Id, a_Name);
  if (DatasetId < 0) {
    return ERROR_DATASET_INVALID;
  }

  /* Make sure the datatype is a float datatype */
  DatatypeId = H5Dget_type(DatasetId);
  if (DatatypeId < 0) {
    H5Dclose(DatasetId);
    return DatatypeId;
  }
  Classtype = H5Tget_class(DatatypeId);
  H5Tclose(DatatypeId);
  if (Classtype != H5T_INTEGER) {
    H5Dclose(DatasetId);
    return ERROR_DATASET_INVALID;
  }

  SpaceId = H5Dget_space(DatasetId);
  if (SpaceId < 0) {
    H5Dclose(DatasetId);
    H5Sclose(SpaceId);
    return ERROR_DATASET_INVALID;
  }

  /* Make sure that the dataspace is simple */
  if (H5Sis_simple(SpaceId) <= 0) {
    H5Dclose(DatasetId);
    H5Sclose(SpaceId);
    return ERROR_DATASET_INVALID;
  }

  /* Make sure the rank is 2 */
  status = xfpGetSimpleDataspaceInfo(SpaceId, &Rank, &Dims, &Maxdims);
  if (Dims) {
    free(Dims);
  }
  if (Maxdims) {
    free(Maxdims);
  }
  if (status < 0 || Rank != 2) {
    H5Sclose(SpaceId);
    H5Dclose(DatasetId);
    return ERROR_DATASET_INVALID;
  }

  /* select the hyperslab */
  Count[0] = a_Num1;
  Count[1] = a_Num2;

  Offset[0] = a_Start1;
  Offset[1] = a_Start2;

  status = H5Sselect_hyperslab(SpaceId, H5S_SELECT_SET, Offset, NULL, Count, 
                               NULL);
  if (status < 0) {
    H5Sclose(SpaceId);
    H5Dclose(DatasetId);
    return status;
  }

  /* Make sure the portion is valid */
  if (H5Sselect_valid(SpaceId) <= 0) {
    return ERROR_DATASET_SIZE_INCORRECT;
  }

  /* Create a memory dataspace for the data */
  Memsize = H5Sget_select_npoints(SpaceId);
  MemspaceId = H5Screate_simple(1, &Memsize, NULL);
  if (MemspaceId < 0) {
    H5Sclose(SpaceId);
    H5Dclose(DatasetId);
    return MemspaceId;
  }

  /* Read the data */
  status = H5Dread(DatasetId, H5T_NATIVE_INT, MemspaceId, SpaceId,
                   H5P_DEFAULT, a_Array);

  H5Sclose(SpaceId);
  H5Dclose(DatasetId);
  H5Sclose(MemspaceId);
 
  return status;
#else
  return -1;
#endif
} /* xfpReadDset2DIntPortion */


/******************************************************************************
 * FUNCTION  xfpReadDset2DDoublePortion
 * PURPOSE   Read a portion of a 2D double array from an HDF5 dataset
 * NOTES     returns positive if successful, otherwise negative
 *           a_Start1 is the starting index for the first array index, 
 *           a_Num1 is number to read in first array index
 *           a_Start2 and a_Num2 are same as a_Start1 and a_Num1 but for index2
 *           Data is always stored in a 1D array
 *           Data must already be allocated of size a_Num1 * a_Num2
 ****************2***********************************************************/
XMDF_API int xfpReadDset2DDoublePortion (xid a_Id, const char *a_Name, 
                hssize_t a_Start1, hsize_t a_Num1,
                hssize_t a_Start2, hsize_t a_Num2, double *a_Array)
{
#ifndef _WIN64
  herr_t    status;
  xid       DatasetId, SpaceId, DatatypeId, MemspaceId;
  int       Rank;
  hsize_t   Count[2], *Dims = NULL, *Maxdims = NULL;
  hssize_t  Offset[2], Memsize;
  H5T_class_t Classtype;

  DatasetId = H5Dopen(a_Id, a_Name);
  if (DatasetId < 0) {
    return ERROR_DATASET_INVALID;
  }

  /* Make sure the datatype is a float datatype */
  DatatypeId = H5Dget_type(DatasetId);
  if (DatatypeId < 0) {
    H5Dclose(DatasetId);
    return DatatypeId;
  }
  Classtype = H5Tget_class(DatatypeId);
  H5Tclose(DatatypeId);
  if (Classtype != H5T_FLOAT) {
    H5Dclose(DatasetId);
    return ERROR_DATASET_INVALID;
  }

  SpaceId = H5Dget_space(DatasetId);
  if (SpaceId < 0) {
    H5Dclose(DatasetId);
    H5Sclose(SpaceId);
    return ERROR_DATASET_INVALID;
  }

  /* Make sure that the dataspace is simple */
  if (H5Sis_simple(SpaceId) <= 0) {
    H5Dclose(DatasetId);
    H5Sclose(SpaceId);
    return ERROR_DATASET_INVALID;
  }

  /* Make sure the rank is 2 */
  status = xfpGetSimpleDataspaceInfo(SpaceId, &Rank, &Dims, &Maxdims);

  if (Maxdims) {
    free(Maxdims);
  }
  if (status < 0 || Rank != 2) {
    H5Sclose(SpaceId);
    H5Dclose(DatasetId);
    return ERROR_DATASET_INVALID;
  }

  /* select the hyperslab */
  Count[0] = a_Num1;
  Count[1] = a_Num2;

  Offset[0] = a_Start1;
  Offset[1] = a_Start2;

  status = H5Sselect_hyperslab(SpaceId, H5S_SELECT_SET, Offset, NULL, Count, 
                              NULL);
  if (status < 0) {
    H5Sclose(SpaceId);
    H5Dclose(DatasetId);
    return status;
  }

  /* Make sure the portion is valid */
  if (H5Sselect_valid(SpaceId) <= 0) {
    return ERROR_DATASET_SIZE_INCORRECT;
  }

  /* Create a memory dataspace for the data */
  Memsize = H5Sget_select_npoints(SpaceId);
  MemspaceId = H5Screate_simple(1, &Memsize, NULL);
  if (MemspaceId < 0) {
    H5Sclose(SpaceId);
    H5Dclose(DatasetId);
    return MemspaceId;
  }

  /* Read the data */
  status = H5Dread(DatasetId, H5T_NATIVE_DOUBLE, MemspaceId, SpaceId,
                  H5P_DEFAULT, a_Array);

  if (Dims) {
    free(Dims);
  }

  H5Sclose(SpaceId);
  H5Dclose(DatasetId);
  H5Sclose(MemspaceId);
 
  return status;
#else
  return -1;
#endif
} /* xfpReadDset2DDoublePortion */
/******************************************************************************
 * FUNCTION  xfpReadDset2DUCharPortion
 * PURPOSE   Read a portion of a 2D float array from an HDF5 dataset
 * NOTES     returns positive if successful, otherwise negative
 *           a_Start1 is the starting index for the first array index, 
 *           a_Num1 is number to read in first array index
 *           a_Start2 and a_Num2 are same as a_Start1 and a_Num1 but for index2
 *           Data is always stored in a 1D array
 *           Data must already be allocated of size a_Num1 * a_Num2
 ****************2***********************************************************/
XMDF_API int xfpReadDset2DUCharPortion (xid a_Id, const char *a_Name, 
                hssize_t a_Start1, hsize_t a_Num1,
                hssize_t a_Start2, hsize_t a_Num2, xmbool *a_Array)
{
#ifndef _WIN64
  herr_t    status;
  xid       DatasetId, SpaceId, MemspaceId, DatatypeId;
  int       Rank;
  hsize_t   Count[2], *Dims = NULL, *Maxdims = NULL;
  hssize_t  Offset[2], Memsize;
  H5T_class_t Classtype;

  DatasetId = H5Dopen(a_Id, a_Name);
  if (DatasetId < 0) {
    return ERROR_DATASET_INVALID;
  }

  /* Make sure the datatype is a float datatype */
  DatatypeId = H5Dget_type(DatasetId);
  if (DatatypeId < 0) {
    H5Dclose(DatasetId);
    return DatatypeId;
  }
  Classtype = H5Tget_class(DatatypeId);
  H5Tclose(DatatypeId);
  if (Classtype != H5T_INTEGER) {
    H5Dclose(DatasetId);
    return ERROR_DATASET_INVALID;
  }

  SpaceId = H5Dget_space(DatasetId);
  if (SpaceId < 0) {
    H5Sclose(SpaceId);
    return ERROR_DATASET_INVALID;
  }

  /* Make sure that the dataspace is simple */
  if (H5Sis_simple(SpaceId) <= 0) {
    H5Dclose(DatasetId);
    H5Sclose(SpaceId);
    return ERROR_DATASET_INVALID;
  }

  /* Make sure the rank is 2 */
  status = xfpGetSimpleDataspaceInfo(SpaceId, &Rank, &Dims, &Maxdims);
  if (Dims) {
    free(Dims);
  }
  if (Maxdims) {
    free(Maxdims);
  }
  if (status < 0 || Rank != 2) {
    H5Sclose(SpaceId);
    H5Dclose(DatasetId);
    return ERROR_DATASET_INVALID;
  }

  /* select the hyperslab */
  Count[0] = a_Num1;
  Count[1] = a_Num2;

  Offset[0] = a_Start1;
  Offset[1] = a_Start2;

  status = H5Sselect_hyperslab(SpaceId, H5S_SELECT_SET, Offset, NULL, Count, 
                               NULL);
  if (status < 0) {
    H5Sclose(SpaceId);
    H5Dclose(DatasetId);
    return status;
  }

  /* Make sure the portion is valid */
  if (H5Sselect_valid(SpaceId) <= 0) {
    return ERROR_DATASET_SIZE_INCORRECT;
  }

  /* Create a memory dataspace for the data */
  Memsize = H5Sget_select_npoints(SpaceId);
  MemspaceId = H5Screate_simple(1, &Memsize, NULL);
  if (MemspaceId < 0) {
    H5Sclose(SpaceId);
    H5Dclose(DatasetId);
    return MemspaceId;
  }

  /* Read the data */
  status = H5Dread(DatasetId, H5T_NATIVE_UCHAR, MemspaceId, SpaceId,
                   H5P_DEFAULT, a_Array);

  H5Sclose(SpaceId);
  H5Dclose(DatasetId);
  H5Sclose(MemspaceId);
 
  return status;
#else
  return -1;
#endif
} /* xfpReadDset2DUCharPortion */
/******************************************************************************
 * FUNCTION  xfpReadDset3DDoublePortion
 * PURPOSE   Read a portion of a 3D double array from an HDF5 dataset
 * NOTES     returns positive if successful, otherwise negative
 *           a_Start1 is the starting index for the first array index, 
 *           a_Num1 is number to read in first array index
 *           a_Start2 and a_Num2 are same as a_Start1 and a_Num1 but for index2
 *           Data is always stored in a 1D array
 *           Data must already be allocated of size a_Num1 * a_Num2 * a_Num3
 ****************2***********************************************************/
XMDF_API int xfpReadDset3DDoublePortion (xid a_Id, const char *a_Name, 
                hssize_t a_Start1, hsize_t a_Num1,
                hssize_t a_Start2, hsize_t a_Num2, hssize_t a_Start3,
                hsize_t a_Num3, double *a_Array)
{
#ifndef _WIN64
  herr_t    status = 1;
  xid       DatasetId, SpaceId, MemspaceId, DatatypeId;
  int       Rank;
  hsize_t   Count[3], *Dims = NULL, *Maxdims = NULL;
  hssize_t  Offset[3], Memsize = 0;
  H5T_class_t Classtype;

  DatasetId = H5Dopen(a_Id, a_Name);
  if (DatasetId < 0) {
    return ERROR_DATASET_INVALID;
  }

  /* Make sure the datatype is a float datatype */
  DatatypeId = H5Dget_type(DatasetId);
  if (DatatypeId < 0) {
    H5Dclose(DatasetId);
    return DatatypeId;
  }
  Classtype = H5Tget_class(DatatypeId);
  H5Tclose(DatatypeId);
  if (Classtype != H5T_FLOAT) {
    H5Dclose(DatasetId);
    return ERROR_DATASET_INVALID;
  }

  SpaceId = H5Dget_space(DatasetId);
  if (SpaceId < 0) {
    H5Sclose(SpaceId);
    return ERROR_DATASET_INVALID;
  }

  /* Make sure that the dataspace is simple */
  if (H5Sis_simple(SpaceId) <= 0) {
    H5Dclose(DatasetId);
    H5Sclose(SpaceId);
    return ERROR_DATASET_INVALID;
  }

  /* Make sure the rank is 3 */
  status = xfpGetSimpleDataspaceInfo(SpaceId, &Rank, &Dims, &Maxdims);
  if (Dims) {
    free(Dims);
  }
  if (Maxdims) {
    free(Maxdims);
  }
  if (status < 0 || Rank != 3) {
    H5Sclose(SpaceId);
    H5Dclose(DatasetId);
    return ERROR_DATASET_INVALID;
  }

  /* select the hyperslab */
  Count[0] = a_Num1;
  Count[1] = a_Num2;
  Count[2] = a_Num3;

  Offset[0] = a_Start1;
  Offset[1] = a_Start2;
  Offset[2] = a_Start3;

  status = H5Sselect_hyperslab(SpaceId, H5S_SELECT_SET, Offset, NULL, Count, 
                               NULL);
  if (status < 0) {
    H5Sclose(SpaceId);
    H5Dclose(DatasetId);
    return status;
  }

  /* Make sure the portion is valid */
  if (H5Sselect_valid(SpaceId) <= 0) {
    return ERROR_DATASET_SIZE_INCORRECT;
  }

  /* Create a memory dataspace for the data */
  /* don't do this way because it fails in release mode */
  Memsize = H5Sget_select_npoints(SpaceId);

  Memsize = a_Num1 * a_Num2 *a_Num3;
  MemspaceId = H5Screate_simple(1, &Memsize, NULL);
  if (MemspaceId < 0) {
    H5Sclose(SpaceId);
    H5Dclose(DatasetId);
  }

  status = H5Dread(DatasetId, H5T_NATIVE_DOUBLE, MemspaceId, SpaceId,
                   H5P_DEFAULT, a_Array);

  H5Sclose(SpaceId);
  H5Sclose(MemspaceId);
  H5Dclose(DatasetId);

 
  return (int)status;
#else
  return -1;
#endif
} /* xfpReadDset3DDoublePortion */
/******************************************************************************
 * FUNCTION  xfpReadDset3DDoublePortions
 * PURPOSE   Read several portions of a 3D double array from an HDF5 dataset
 * NOTES     returns positive if successful, otherwise negative
 *           a_Start1 is the starting index for the first array index, 
 *           a_Num1 is number to read in first array index
 *           a_Start2 and a_Num2 are same as a_Start1 and a_Num1 but for index2
 *           Data is always stored in a 1D array
 *           Data must already be allocated of size a_Num1 * a_Num2 * a_Num3
 ****************2***********************************************************/
XMDF_API int xfpReadDset3DDoublePortions (xid a_Id, const char *a_Name,
                int nPortions, hssize_t *a_Start1, hsize_t *a_Num1,
                hssize_t *a_Start2, hsize_t *a_Num2, hssize_t *a_Start3,
                hsize_t *a_Num3, double *a_Array)
{
#ifndef _WIN64
  herr_t    status = 1;
  xid       DatasetId, SpaceId, MemspaceId, DatatypeId;
  int       i, Rank;
  hsize_t   Count[3], *Dims = NULL, *Maxdims = NULL;
  hssize_t  Offset[3], Memsize = 0;
  H5T_class_t Classtype;

  DatasetId = H5Dopen(a_Id, a_Name);
  if (DatasetId < 0) {
    return ERROR_DATASET_INVALID;
  }

  /* Make sure the datatype is a float datatype */
  DatatypeId = H5Dget_type(DatasetId);
  if (DatatypeId < 0) {
    H5Dclose(DatasetId);
    return DatatypeId;
  }
  Classtype = H5Tget_class(DatatypeId);
  H5Tclose(DatatypeId);
  if (Classtype != H5T_FLOAT) {
    H5Dclose(DatasetId);
    return ERROR_DATASET_INVALID;
  }

  SpaceId = H5Dget_space(DatasetId);
  if (SpaceId < 0) {
    H5Sclose(SpaceId);
    return ERROR_DATASET_INVALID;
  }

  /* Make sure that the dataspace is simple */
  if (H5Sis_simple(SpaceId) <= 0) {
    H5Dclose(DatasetId);
    H5Sclose(SpaceId);
    return ERROR_DATASET_INVALID;
  }

  /* Make sure the rank is 3 */
  status = xfpGetSimpleDataspaceInfo(SpaceId, &Rank, &Dims, &Maxdims);
  if (Dims) {
    free(Dims);
  }
  if (Maxdims) {
    free(Maxdims);
  }
  if (status < 0 || Rank != 3) {
    H5Sclose(SpaceId);
    H5Dclose(DatasetId);
    return ERROR_DATASET_INVALID;
  }

  /* select the first hyperslab */
  Count[0] = a_Num1[0];
  Count[1] = a_Num2[0];
  Count[2] = a_Num3[0];

  Offset[0] = a_Start1[0];
  Offset[1] = a_Start2[0];
  Offset[2] = a_Start3[0];

  status = H5Sselect_hyperslab(SpaceId, H5S_SELECT_SET, Offset, NULL, Count, 
                               NULL);
  if (status < 0) {
    H5Sclose(SpaceId);
    H5Dclose(DatasetId);
    return status;
  }

  /* select the successive hyperslabs */
  for (i=1; i<nPortions; ++i) {
    Count[0] = a_Num1[i];
    Count[1] = a_Num2[i];
    Count[2] = a_Num3[i];

    Offset[0] = a_Start1[i];
    Offset[1] = a_Start2[i];
    Offset[2] = a_Start3[i];

    status = H5Sselect_hyperslab(SpaceId, H5S_SELECT_OR, Offset, NULL, Count, 
                                 NULL);
    if (status < 0) {
      H5Sclose(SpaceId);
      H5Dclose(DatasetId);
      return status;
    }
  }

  /* Make sure the portion is valid */
  if (H5Sselect_valid(SpaceId) <= 0) {
    return ERROR_DATASET_SIZE_INCORRECT;
  }

  /* Create a memory dataspace for the data */
  /* don't do this way because it fails in release mode */
  Memsize = H5Sget_select_npoints(SpaceId);
  MemspaceId = H5Screate_simple(1, &Memsize, NULL);
  if (MemspaceId < 0) {
    H5Sclose(SpaceId);
    H5Dclose(DatasetId);
  }

  status = H5Dread(DatasetId, H5T_NATIVE_DOUBLE, MemspaceId, SpaceId,
                   H5P_DEFAULT, a_Array);

  H5Sclose(SpaceId);
  H5Sclose(MemspaceId);
  H5Dclose(DatasetId);

 
  return (int)status;
#else
  return -1;
#endif
} /* xfpReadDset3DDoublePortions */
/******************************************************************************
 * FUNCTION  xfpReadDset3DFloatPortion
 * PURPOSE   Read a portion of a 3D float array from an HDF5 dataset
 * NOTES     returns positive if successful, otherwise negative
 *           a_Start1 is the starting index for the first array index, 
 *           a_Num1 is number to read in first array index
 *           a_Start2 and a_Num2 are same as a_Start1 and a_Num1 but for index2
 *           Data is always stored in a 1D array
 *           Data must already be allocated of size a_Num1 * a_Num2 * a_Num3
 ****************2***********************************************************/
XMDF_API int xfpReadDset3DFloatPortion (xid a_Id, const char *a_Name, 
                hssize_t a_Start1, hsize_t a_Num1,
                hssize_t a_Start2, hsize_t a_Num2, hssize_t a_Start3,
                hsize_t a_Num3, float *a_Array)
{
#ifndef _WIN64
  herr_t    status = 1;
  xid       DatasetId, SpaceId, MemspaceId, DatatypeId;
  int       Rank;
  hsize_t   Count[3], *Dims = NULL, *Maxdims = NULL;
  hssize_t  Offset[3], Memsize = 0;
  H5T_class_t Classtype;

  DatasetId = H5Dopen(a_Id, a_Name);
  if (DatasetId < 0) {
    return ERROR_DATASET_INVALID;
  }

  /* Make sure the datatype is a float datatype */
  DatatypeId = H5Dget_type(DatasetId);
  if (DatatypeId < 0) {
    H5Dclose(DatasetId);
    return DatatypeId;
  }
  Classtype = H5Tget_class(DatatypeId);
  H5Tclose(DatatypeId);
  if (Classtype != H5T_FLOAT) {
    H5Dclose(DatasetId);
    return ERROR_DATASET_INVALID;
  }

  SpaceId = H5Dget_space(DatasetId);
  if (SpaceId < 0) {
    H5Sclose(SpaceId);
    return ERROR_DATASET_INVALID;
  }

  /* Make sure that the dataspace is simple */
  if (H5Sis_simple(SpaceId) <= 0) {
    H5Dclose(DatasetId);
    H5Sclose(SpaceId);
    return ERROR_DATASET_INVALID;
  }

  /* Make sure the rank is 3 */
  status = xfpGetSimpleDataspaceInfo(SpaceId, &Rank, &Dims, &Maxdims);
  if (Dims) {
    free(Dims);
  }
  if (Maxdims) {
    free(Maxdims);
  }
  if (status < 0 || Rank != 3) {
    H5Sclose(SpaceId);
    H5Dclose(DatasetId);
    return ERROR_DATASET_INVALID;
  }

  /* select the hyperslab */
  Count[0] = a_Num1;
  Count[1] = a_Num2;
  Count[2] = a_Num3;

  Offset[0] = a_Start1;
  Offset[1] = a_Start2;
  Offset[2] = a_Start3;

  status = H5Sselect_hyperslab(SpaceId, H5S_SELECT_SET, Offset, NULL, Count, 
                               NULL);
  if (status < 0) {
    H5Sclose(SpaceId);
    H5Dclose(DatasetId);
    return status;
  }

  /* Make sure the portion is valid */
  if (H5Sselect_valid(SpaceId) <= 0) {
    return ERROR_DATASET_SIZE_INCORRECT;
  }

  /* Create a memory dataspace for the data */
  /* don't do this way because it fails in release mode */
  Memsize = H5Sget_select_npoints(SpaceId);

  Memsize = a_Num1 * a_Num2 *a_Num3;
  MemspaceId = H5Screate_simple(1, &Memsize, NULL);
  if (MemspaceId < 0) {
    H5Sclose(SpaceId);
    H5Dclose(DatasetId);
  }

  status = H5Dread(DatasetId, H5T_NATIVE_FLOAT, MemspaceId, SpaceId,
                   H5P_DEFAULT, a_Array);

  H5Sclose(SpaceId);
  H5Sclose(MemspaceId);
  H5Dclose(DatasetId);

 
  return (int)status;
#else
  return -1;
#endif
} /* xfpReadDset3DFloatPortion */
/* ---------------------------------------------------------------------------
FUNCTION  xfpWriteDsetInt
PURPOSE   
NOTES     
----------------------------------------------------------------------------*/
XMDF_API int xfpWriteDsetInt (xid a_Id, const char *a_Name, hid_t a_Datatype,
                            const XDatasetParams *a_Params, const int *a_Array)
{
#ifndef _WIN64
  hid_t     DataspaceId, CreateParamsId, DatasetId;
  herr_t    status;

  /* create the dataspace */
  DataspaceId = xfpDsetParamsCreateSpace(a_Params);
  if (DataspaceId < 0) {
    return DataspaceId;
  }

  /* setup the dataset creation properties */
  CreateParamsId = xfpDsetParamsCreateProps(a_Params);

  /* Create the Dataset */
  DatasetId = H5Dcreate(a_Id, a_Name, a_Datatype, DataspaceId, CreateParamsId);
  if (DatasetId < 0) {
    H5Sclose(DataspaceId);
    H5Pclose(CreateParamsId);
    return DatasetId;
  }

  status = H5Dwrite(DatasetId, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, 
                    H5P_DEFAULT, a_Array);

  H5Dclose(DatasetId);
  H5Sclose(DataspaceId);
  H5Pclose(CreateParamsId);
 
  return status;
#else
  return -1;
#endif
} /* xfpWriteDsetInt */
/* ---------------------------------------------------------------------------
FUNCTION  xfpWriteDsetUInt
PURPOSE   Write a UInt array to a HDF5 dataset
NOTES     
----------------------------------------------------------------------------*/
XMDF_API int xfpWriteDsetUInt (xid a_Id, const char *a_Name, hid_t a_Datatype,
              const XDatasetParams *a_Params, const unsigned int *a_Array)
{
#ifndef _WIN64
  hid_t     DataspaceId, CreateParamsId, DatasetId;
  herr_t    status;

  /* create the dataspace */
  DataspaceId = xfpDsetParamsCreateSpace(a_Params);
  if (DataspaceId < 0) {
    return DataspaceId;
  }

  /* setup the dataset creation properties */
  CreateParamsId = xfpDsetParamsCreateProps(a_Params);

  /* Create the Dataset */
  DatasetId = H5Dcreate(a_Id, a_Name, a_Datatype, DataspaceId, CreateParamsId);
  if (DatasetId < 0) {
    H5Sclose(DataspaceId);
    H5Pclose(CreateParamsId);
    return DatasetId;
  }

  status = H5Dwrite(DatasetId, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, 
                    H5P_DEFAULT, a_Array);

  H5Dclose(DatasetId);
  H5Sclose(DataspaceId);
  H5Pclose(CreateParamsId);
 
  return status;
#else
  return -1;
#endif
} /* xfpWriteDsetUInt */
/******************************************************************************
 * FUNCTION  xfpWriteDsetFloat
 * PURPOSE   Write a float array to a HDF5 dataset
 * NOTES     returns positive if successful, otherwise negative
 *           The array must be of continuous memory
 ******************2***********************************************************/
XMDF_API int xfpWriteDsetFloat (xid a_Id, const char *a_Name,
                           const XDatasetParams *a_Params, const float *a_Array)
{
#ifndef _WIN64
  hid_t     DataspaceId, CreateParamsId, DatasetId;
  herr_t    status=0;

  /* create the dataspace */
  DataspaceId = xfpDsetParamsCreateSpace(a_Params);
  if (DataspaceId < 0) {
    return DataspaceId;
  }

  /* setup the dataset creation properties */
  CreateParamsId = xfpDsetParamsCreateProps(a_Params);

  /* Create the Dataset */
  DatasetId = H5Dcreate(a_Id, a_Name, fg_Floattype, DataspaceId, CreateParamsId);
  if (DatasetId < 0) {
    H5Sclose(DataspaceId);
    H5Pclose(CreateParamsId);
    return DatasetId;
  }

  if (!xfpGetRunningParallel() || xfpGetParallelRank() == 0) {
    status = H5Dwrite(DatasetId, H5T_NATIVE_FLOAT, H5S_ALL, DataspaceId, 
                      H5P_DEFAULT, a_Array);
  }
  else {
    status=0;
  }

  H5Dclose(DatasetId);
  H5Sclose(DataspaceId);
  H5Pclose(CreateParamsId);
 
  return status;
#else
  return -1;
#endif
} /* xfpWriteDsetFloat */
/******************************************************************************
 * FUNCTION  xfpWriteDsetDouble
 * PURPOSE   Write a double array to a HDF5 dataset
 * NOTES     returns positive if successful, otherwise negative
 *           The array must be of continuous memory
 ******************2***********************************************************/
XMDF_API int xfpWriteDsetDouble (xid a_Id, const char *a_Name,
                          const XDatasetParams *a_Params, const double *a_Array)
{
#ifndef _WIN64
  hid_t     DataspaceId, CreateParamsId, DatasetId;
  herr_t    status;

  /* create the dataspace */
  DataspaceId = xfpDsetParamsCreateSpace(a_Params);
  if (DataspaceId < 0) {
    return DataspaceId;
  }

  /* setup the dataset creation properties */
  CreateParamsId = xfpDsetParamsCreateProps(a_Params);

  /* Create the Dataset */
  DatasetId = H5Dcreate(a_Id, a_Name, fg_Doubletype, DataspaceId, CreateParamsId);
  if (DatasetId < 0) {
    H5Sclose(DataspaceId);
    H5Pclose(CreateParamsId);
    return DatasetId;
  }

  if (!xfpGetRunningParallel() || xfpGetParallelRank() == 0) {
    status = H5Dwrite(DatasetId, H5T_NATIVE_DOUBLE, H5S_ALL, DataspaceId, 
                      H5P_DEFAULT, a_Array);
  }
  else {
    status=0;
  }

  H5Dclose(DatasetId);
  H5Sclose(DataspaceId);
  H5Pclose(CreateParamsId);
 
  return status;
#else
  return -1;
#endif
} /* xfpWriteDsetDouble */
/******************************************************************************
 * FUNCTION  xfpWriteDsetUChar
 * PURPOSE   Write a unsigned char array to a HDF5 dataset
 * NOTES     returns positive if successful, otherwise negative
 *           The array must be of continuous memory
 ******************2***********************************************************/
XMDF_API int xfpWriteDsetUChar(xid a_Id, const char *a_Name, hid_t a_Datatype,
                   const XDatasetParams *a_Params, const unsigned char *a_Array)
{
#ifndef _WIN64
  hid_t     DataspaceId, CreateParamsId, DatasetId;
  herr_t    status;

  /* create the dataspace */
  DataspaceId = xfpDsetParamsCreateSpace(a_Params);
  if (DataspaceId < 0) {
    return DataspaceId;
  }

  /* setup the dataset creation properties */
  CreateParamsId = xfpDsetParamsCreateProps(a_Params);

  /* Create the Dataset */
  DatasetId = H5Dcreate(a_Id, a_Name, a_Datatype, DataspaceId, CreateParamsId);
  if (DatasetId < 0) {
    H5Sclose(DataspaceId);
    H5Pclose(CreateParamsId);
    return DatasetId;
  }

  status = H5Dwrite(DatasetId, H5T_NATIVE_UCHAR, H5S_ALL, DataspaceId, 
                    H5P_DEFAULT, a_Array);

  H5Dclose(DatasetId);
  H5Sclose(DataspaceId);
  H5Pclose(CreateParamsId);
 
  return status;
#else
  return -1;
#endif
} /* xfpWriteDsetUChar */
/******************************************************************************
 * FUNCTION  xfpAppendDset1DDouble
 * PURPOSE   
 * NOTES     
 ******************2***********************************************************/
XMDF_API int xfpAppendDset1DDouble (xid a_Id, const char *a_Name, 
                                   hsize_t a_NumNew, const double *a_Array)
{
#ifndef _WIN64
  herr_t     status = XTRUE;
  hsize_t   *OldDims = NULL, *OldMaxDims = NULL, NewDims;
  H5T_class_t Classtype;
  xid        DatasetId, DataspaceId, DatatypeId, FilespaceId, MemspaceId;
  int        Rank;
  hssize_t   offset;

  DatasetId = H5Dopen(a_Id, a_Name);
  if (DatasetId < 0) {
    return DatasetId;
  }

  DataspaceId = H5Dget_space(DatasetId);
  if (DataspaceId < 0) {
    H5Dclose(DatasetId);
    return DataspaceId;
  }

  /* Dataspace must be simple and have a rank of 1 */
  if (H5Sis_simple(DataspaceId) <= 0) {
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    return ERROR_DATASET_INVALID;
  }

  status = xfpGetSimpleDataspaceInfo(DataspaceId, &Rank, &OldDims, &OldMaxDims);
  if (status < 0) {
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    return status;
  }

  /* This function is only for 1D arrays of rank 1 and MaxDims must be larger */
  /* than the old dimension plus the new values */
  NewDims = OldDims[0] + a_NumNew;
  if (Rank != 1 || OldMaxDims[0] < NewDims) {
    free(OldDims);
    free(OldMaxDims);
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    return ERROR_DATASET_INVALID;
  }

  /* Make sure the datatype is a float datatype */
  DatatypeId = H5Dget_type(DatasetId);
  if (DatatypeId < 0) {
    free(OldDims);
    free(OldMaxDims);
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    return DatatypeId;
  }
  Classtype = H5Tget_class(DatatypeId);
  H5Tclose(DatatypeId);
  if (Classtype != H5T_FLOAT) {
    free(OldDims);
    free(OldMaxDims);
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    return ERROR_DATASET_INVALID;
  }

  /* extend the dataset */
  status = H5Dextend(DatasetId, &NewDims);
  if (status < 0) {
    free(OldDims);
    free(OldMaxDims);
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    return ERROR_DATASET_INVALID;
  }

  /* select a hyperslab */
  FilespaceId = H5Dget_space(DatasetId);
  if (FilespaceId < 0) {
    free(OldDims);
    free(OldMaxDims);
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    return FilespaceId;
  }

  offset = OldDims[0];
  status = H5Sselect_hyperslab(FilespaceId, H5S_SELECT_SET, &offset, NULL, 
                               &a_NumNew, NULL);
  if (status < 0) {
    free(OldDims);
    free(OldMaxDims);
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    H5Sclose(FilespaceId);
    return status;
  }

  /* create the memory dataspace */
  MemspaceId = H5Screate_simple(1, &a_NumNew, NULL);
  if (MemspaceId < 0) {
    free(OldDims);
    free(OldMaxDims);
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    H5Sclose(FilespaceId);
    return MemspaceId;
  }

    /* write the data */
  if (!xfpGetRunningParallel() || xfpGetParallelRank() == 0) {
    status = H5Dwrite(DatasetId, H5T_NATIVE_DOUBLE, MemspaceId, FilespaceId,
                      H5P_DEFAULT, a_Array);
  }
  else {
    status=0;
  }

  free(OldDims);
  free(OldMaxDims);
  H5Sclose(MemspaceId);
  H5Dclose(DatasetId);
  H5Sclose(DataspaceId);
  H5Sclose(FilespaceId);

  return status;
#else
  return -1;
#endif
} /* xfpAppendDset1DDouble */
/******************************************************************************
 * FUNCTION  xfpAppendDset1DFloat
 * PURPOSE   
 * NOTES     
 ******************2***********************************************************/
XMDF_API int xfpAppendDset1DFloat (xid a_Id, const char *a_Name, hsize_t a_NumNew,
                            const float *a_Array)
{
#ifndef _WIN64
  herr_t     status = XTRUE;
  hsize_t   *OldDims = NULL, *OldMaxDims = NULL, NewDims;
  H5T_class_t Classtype;
  xid        DatasetId, DataspaceId, DatatypeId, FilespaceId, MemspaceId;
  int        Rank;
  hssize_t   offset;

  DatasetId = H5Dopen(a_Id, a_Name);
  if (DatasetId < 0) {
    return DatasetId;
  }

  DataspaceId = H5Dget_space(DatasetId);
  if (DataspaceId < 0) {
    H5Dclose(DatasetId);
    return DataspaceId;
  }

  /* Dataspace must be simple and have a rank of 1 */
  if (H5Sis_simple(DataspaceId) <= 0) {
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    return ERROR_DATASET_INVALID;
  }

  status = xfpGetSimpleDataspaceInfo(DataspaceId, &Rank, &OldDims, &OldMaxDims);
  if (status < 0) {
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    return status;
  }

  /* This function is only for 1D arrays of rank 1 and MaxDims must be larger */
  /* than the old dimension plus the new values */
  NewDims = OldDims[0] + a_NumNew;
  if (Rank != 1 || OldMaxDims[0] < NewDims) {
    free(OldDims);
    free(OldMaxDims);
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    return ERROR_DATASET_INVALID;
  }

  /* Make sure the datatype is a float datatype */
  DatatypeId = H5Dget_type(DatasetId);
  if (DatatypeId < 0) {
    free(OldDims);
    free(OldMaxDims);
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    return DatatypeId;
  }
  Classtype = H5Tget_class(DatatypeId);
  H5Tclose(DatatypeId);
  if (Classtype != H5T_FLOAT) {
    free(OldDims);
    free(OldMaxDims);
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    return ERROR_DATASET_INVALID;
  }

  /* extend the dataset */
  status = H5Dextend(DatasetId, &NewDims);
  if (status < 0) {
    free(OldDims);
    free(OldMaxDims);
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    return ERROR_DATASET_INVALID;
  }

  /* select a hyperslab */
  FilespaceId = H5Dget_space(DatasetId);
  if (FilespaceId < 0) {
    free(OldDims);
    free(OldMaxDims);
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    return FilespaceId;
  }

  offset = OldDims[0];
  status = H5Sselect_hyperslab(FilespaceId, H5S_SELECT_SET, &offset, NULL, 
                               &a_NumNew, NULL);
  if (status < 0) {
    free(OldDims);
    free(OldMaxDims);
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    H5Sclose(FilespaceId);
    return status;
  }

  /* create the memory dataspace */
  MemspaceId = H5Screate_simple(1, &a_NumNew, NULL);
  if (MemspaceId < 0) {
    free(OldDims);
    free(OldMaxDims);
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    H5Sclose(FilespaceId);
    return status;
  }

  /* write the data */
  if (!xfpGetRunningParallel() || xfpGetParallelRank() == 0) {
    status = H5Dwrite(DatasetId, H5T_NATIVE_FLOAT, MemspaceId, FilespaceId,
                      H5P_DEFAULT, a_Array);
  }
  else {
    status=0;
  }
  
  free(OldDims);
  free(OldMaxDims);
  H5Sclose(MemspaceId);
  H5Dclose(DatasetId);
  H5Sclose(DataspaceId);
  H5Sclose(FilespaceId);

  return status;
#else
  return -1;
#endif
} /* xfpAppendDset1DFloat */
/******************************************************************************
 * FUNCTION  xfpAppendDset1DFloat
 * PURPOSE   
 * NOTES     
 ******************2***********************************************************/
XMDF_API int xfpAppendDset1DInt (xid a_Id, const char *a_Name, hsize_t a_NumNew,
                                 const int *a_Array)
{
#ifndef _WIN64
  herr_t     status = XTRUE;
  hsize_t   *OldDims = NULL, *OldMaxDims = NULL, NewDims;
  H5T_class_t Classtype;
  xid        DatasetId, DataspaceId, DatatypeId, FilespaceId, MemspaceId;
  int        Rank;
  hssize_t   offset;

  DatasetId = H5Dopen(a_Id, a_Name);
  if (DatasetId < 0) {
    return DatasetId;
  }

  DataspaceId = H5Dget_space(DatasetId);
  if (DataspaceId < 0) {
    H5Dclose(DatasetId);
    return DataspaceId;
  }

  /* Dataspace must be simple and have a rank of 1 */
  if (H5Sis_simple(DataspaceId) <= 0) {
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    return ERROR_DATASET_INVALID;
  }

  status = xfpGetSimpleDataspaceInfo(DataspaceId, &Rank, &OldDims, &OldMaxDims);
  if (status < 0) {
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    return status;
  }

  /* This function is only for 1D arrays of rank 1 and MaxDims must be larger */
  /* than the old dimension plus the new values */
  NewDims = OldDims[0] + a_NumNew;
  if (Rank != 1 || OldMaxDims[0] < NewDims) {
    free(OldDims);
    free(OldMaxDims);
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    return ERROR_DATASET_INVALID;
  }

  /* Make sure the datatype is a float datatype */
  DatatypeId = H5Dget_type(DatasetId);
  if (DatatypeId < 0) {
    free(OldDims);
    free(OldMaxDims);
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    return DatatypeId;
  }
  Classtype = H5Tget_class(DatatypeId);
  H5Tclose(DatatypeId);
  if (Classtype != H5T_INTEGER) {
    free(OldDims);
    free(OldMaxDims);
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    return ERROR_DATASET_INVALID;
  }

  /* extend the dataset */
  status = H5Dextend(DatasetId, &NewDims);
  if (status < 0) {
    free(OldDims);
    free(OldMaxDims);
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    return ERROR_DATASET_INVALID;
  }

  /* select a hyperslab */
  FilespaceId = H5Dget_space(DatasetId);
  if (FilespaceId < 0) {
    free(OldDims);
    free(OldMaxDims);
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    return FilespaceId;
  }

  offset = OldDims[0];
  status = H5Sselect_hyperslab(FilespaceId, H5S_SELECT_SET, &offset, NULL, 
                               &a_NumNew, NULL);
  if (status < 0) {
    free(OldDims);
    free(OldMaxDims);
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    H5Sclose(FilespaceId);
    return status;
  }

  /* create the memory dataspace */
  MemspaceId = H5Screate_simple(1, &a_NumNew, NULL);
  if (MemspaceId < 0) {
    free(OldDims);
    free(OldMaxDims);
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    H5Sclose(FilespaceId);
    return status;
  }

  /* write the data */
  status = H5Dwrite(DatasetId, H5T_NATIVE_INT, MemspaceId, FilespaceId,
                    H5P_DEFAULT, a_Array);
  free(OldDims);
  free(OldMaxDims);
  H5Sclose(MemspaceId);
  H5Dclose(DatasetId);
  H5Sclose(DataspaceId);
  H5Sclose(FilespaceId);

  return status;
#else
  return -1;
#endif
} /* xfpAppendDset1DInt */
/******************************************************************************
 * FUNCTION  xfpAppendDset2DFloatFirstDim
 * PURPOSE   
 * NOTES     
 ******************2***********************************************************/
XMDF_API int xfpAppendDset2DFloatFirstDim (xid a_Id, const char *a_Name, 
                                   hsize_t a_NumNew, hsize_t a_NumDim2,
                                   const float *a_Array)
{
#ifndef _WIN64
  herr_t     status = XTRUE;
  hsize_t   *OldDims = NULL, *OldMaxDims = NULL, NewDims[2], MemDims[2];
  H5T_class_t Classtype;
  xid        DatasetId, DataspaceId, DatatypeId, FilespaceId, MemspaceId;
  int        Rank;
  hssize_t   offset[2];
  hsize_t    count[2];

  DatasetId = H5Dopen(a_Id, a_Name);
  if (DatasetId < 0) {
    return DatasetId;
  }

  DataspaceId = H5Dget_space(DatasetId);
  if (DataspaceId < 0) {
    H5Dclose(DatasetId);
    return DataspaceId;
  }

  /* Dataspace must be simple and have a rank of 2 */
  if (H5Sis_simple(DataspaceId) <= 0) {
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    return ERROR_DATASET_INVALID;
  }

  status = xfpGetSimpleDataspaceInfo(DataspaceId, &Rank, &OldDims, &OldMaxDims);
  if (status < 0) {
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    return status;
  }

  /* This function is only for 2D arrays and MaxDims must be larger */
  /* than the old dimension plus the new values in both directions */
  if (Rank != 2) {
    free(OldDims);
    free(OldMaxDims);
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    return ERROR_DATASET_INVALID;
  }
  NewDims[0] = OldDims[0] + a_NumNew;
  NewDims[1] = OldDims[1];
  if (OldMaxDims[0] < NewDims[0] || OldDims[1] != a_NumDim2) {
    free(OldDims);
    free(OldMaxDims);
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    return ERROR_DATASET_INVALID;
  }

  /* Make sure the datatype is a float datatype */
  DatatypeId = H5Dget_type(DatasetId);
  if (DatatypeId < 0) {
    free(OldDims);
    free(OldMaxDims);
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    return DatatypeId;
  }
  Classtype = H5Tget_class(DatatypeId);
  H5Tclose(DatatypeId);
  if (Classtype != H5T_FLOAT) {
    free(OldDims);
    free(OldMaxDims);
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    return ERROR_DATASET_INVALID;
  }

  /* extend the dataset */
  status = H5Dextend(DatasetId, NewDims);
  if (status < 0) {
    free(OldDims);
    free(OldMaxDims);
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    return ERROR_DATASET_INVALID;
  }

  /* select a hyperslab */
  FilespaceId = H5Dget_space(DatasetId);
  if (FilespaceId < 0) {
    free(OldDims);
    free(OldMaxDims);
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    return FilespaceId;
  }

  /* create the memory dataspace */
  MemDims[0] = a_NumNew;
  MemDims[1] = a_NumDim2;
  MemspaceId = H5Screate_simple(2, MemDims, NULL);
  if (MemspaceId < 0) {
    free(OldDims);
    free(OldMaxDims);
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    H5Sclose(FilespaceId);
    return status;
  }

  /* add new information */
  offset[0] = OldDims[0];
  offset[1] = 0;
  count[0] = MemDims[0];
  count[1] = MemDims[1];

  status = H5Sselect_hyperslab(FilespaceId, H5S_SELECT_SET, offset, NULL, 
                               count, NULL);
  if (status >= 0) {
    /* write the data */
    status = H5Dwrite(DatasetId, H5T_NATIVE_FLOAT, MemspaceId, FilespaceId,
                      H5P_DEFAULT, a_Array);
  }

  free(OldDims);
  free(OldMaxDims);
  H5Sclose(MemspaceId);
  H5Dclose(DatasetId);
  H5Sclose(DataspaceId);
  H5Sclose(FilespaceId);

  return status;
#else
  return -1;
#endif
} /* xfpAppendDset2DFloatFirstDim */
/******************************************************************************
 * FUNCTION  xfpAppendDset2DDoubleFirstDim
 * PURPOSE   
 * NOTES     
 ******************2***********************************************************/
XMDF_API int xfpAppendDset2DDoubleFirstDim (xid a_Id, const char *a_Name, 
                                   hsize_t a_NumNew, hsize_t a_NumDim2,
                                   const double *a_Array)
{
#ifndef _WIN64
  herr_t     status = XTRUE;
  hsize_t   *OldDims = NULL, *OldMaxDims = NULL, NewDims[2], MemDims[2];
  H5T_class_t Classtype;
  xid        DatasetId, DataspaceId, DatatypeId, FilespaceId, MemspaceId;
  int        Rank;
  hssize_t   offset[2];
  hsize_t    count[2];

  DatasetId = H5Dopen(a_Id, a_Name);
  if (DatasetId < 0) {
    return DatasetId;
  }

  DataspaceId = H5Dget_space(DatasetId);
  if (DataspaceId < 0) {
    H5Dclose(DatasetId);
    return DataspaceId;
  }

  /* Dataspace must be simple and have a rank of 2 */
  if (H5Sis_simple(DataspaceId) <= 0) {
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    return ERROR_DATASET_INVALID;
  }

  status = xfpGetSimpleDataspaceInfo(DataspaceId, &Rank, &OldDims, &OldMaxDims);
  if (status < 0) {
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    return status;
  }

  /* This function is only for 2D arrays and MaxDims must be larger */
  /* than the old dimension plus the new values in both directions */
  if (Rank != 2) {
    free(OldDims);
    free(OldMaxDims);
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    return ERROR_DATASET_INVALID;
  }
  NewDims[0] = OldDims[0] + a_NumNew;
  NewDims[1] = OldDims[1];
  if (OldMaxDims[0] < NewDims[0] || OldDims[1] != a_NumDim2) {
    free(OldDims);
    free(OldMaxDims);
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    return ERROR_DATASET_INVALID;
  }

  /* Make sure the datatype is a double datatype */
  DatatypeId = H5Dget_type(DatasetId);
  if (DatatypeId < 0) {
    free(OldDims);
    free(OldMaxDims);
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    return DatatypeId;
  }
  Classtype = H5Tget_class(DatatypeId);
  H5Tclose(DatatypeId);
  if (Classtype != H5T_FLOAT) {
    free(OldDims);
    free(OldMaxDims);
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    return ERROR_DATASET_INVALID;
  }

  /* extend the dataset */
  status = H5Dextend(DatasetId, NewDims);
  if (status < 0) {
    free(OldDims);
    free(OldMaxDims);
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    return ERROR_DATASET_INVALID;
  }

  /* select a hyperslab */
  FilespaceId = H5Dget_space(DatasetId);
  if (FilespaceId < 0) {
    free(OldDims);
    free(OldMaxDims);
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    return FilespaceId;
  }

  /* create the memory dataspace */
  MemDims[0] = a_NumNew;
  MemDims[1] = a_NumDim2;
  MemspaceId = H5Screate_simple(2, MemDims, NULL);
  if (MemspaceId < 0) {
    free(OldDims);
    free(OldMaxDims);
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    H5Sclose(FilespaceId);
    return status;
  }

  /* add new information */
  offset[0] = OldDims[0];
  offset[1] = 0;
  count[0] = MemDims[0];
  count[1] = MemDims[1];

  status = H5Sselect_hyperslab(FilespaceId, H5S_SELECT_SET, offset, NULL, 
                               count, NULL);
  if (status >= 0) {
    /* write the data */
    status = H5Dwrite(DatasetId, H5T_NATIVE_DOUBLE, MemspaceId, FilespaceId,
                      H5P_DEFAULT, a_Array);
  }

  free(OldDims);
  free(OldMaxDims);
  H5Sclose(MemspaceId);
  H5Dclose(DatasetId);
  H5Sclose(DataspaceId);
  H5Sclose(FilespaceId);

  return status;
#else
  return -1;
#endif
} /* xfpAppendDset2DDoubleFirstDim */
/******************************************************************************
 * FUNCTION  xfpAppendDset2DDoubleSecondDim
 * PURPOSE   
 * NOTES     
 ******************2***********************************************************/
XMDF_API int xfpAppendDset2DDoubleSecondDim (xid a_Id, const char *a_Name, 
                                   hsize_t a_NumNew, hsize_t a_NumDim1,
                                   const double *a_Array)
{
#ifndef _WIN64
  herr_t     status = XTRUE;
  hsize_t   *OldDims = NULL, *OldMaxDims = NULL, NewDims[2], MemDims[2];
  H5T_class_t Classtype;
  xid        DatasetId, DataspaceId, DatatypeId, FilespaceId, MemspaceId;
  int        Rank;
  hssize_t   offset[2];
  hsize_t    count[2];

  DatasetId = H5Dopen(a_Id, a_Name);
  if (DatasetId < 0) {
    return DatasetId;
  }

  DataspaceId = H5Dget_space(DatasetId);
  if (DataspaceId < 0) {
    H5Dclose(DatasetId);
    return DataspaceId;
  }

  /* Dataspace must be simple and have a rank of 2 */
  if (H5Sis_simple(DataspaceId) <= 0) {
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    return ERROR_DATASET_INVALID;
  }

  status = xfpGetSimpleDataspaceInfo(DataspaceId, &Rank, &OldDims, &OldMaxDims);
  if (status < 0) {
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    return status;
  }

  /* This function is only for 2D arrays and MaxDims must be larger */
  /* than the old dimension plus the new values in both directions */
  if (Rank != 2) {
    free(OldDims);
    free(OldMaxDims);
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    return ERROR_DATASET_INVALID;
  }
  NewDims[0] = OldDims[0];
  NewDims[1] = OldDims[1] + a_NumNew;
  if (OldMaxDims[1] < NewDims[1] || OldDims[0] != a_NumDim1) {
    free(OldDims);
    free(OldMaxDims);
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    return ERROR_DATASET_INVALID;
  }

  /* Make sure the datatype is a double datatype */
  DatatypeId = H5Dget_type(DatasetId);
  if (DatatypeId < 0) {
    free(OldDims);
    free(OldMaxDims);
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    return DatatypeId;
  }
  Classtype = H5Tget_class(DatatypeId);
  H5Tclose(DatatypeId);
  if (Classtype != H5T_FLOAT) {
    free(OldDims);
    free(OldMaxDims);
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    return ERROR_DATASET_INVALID;
  }

  /* extend the dataset */
  status = H5Dextend(DatasetId, NewDims);
  if (status < 0) {
    free(OldDims);
    free(OldMaxDims);
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    return ERROR_DATASET_INVALID;
  }

  /* select a hyperslab */
  FilespaceId = H5Dget_space(DatasetId);
  if (FilespaceId < 0) {
    free(OldDims);
    free(OldMaxDims);
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    return FilespaceId;
  }

  /* create the memory dataspace */
  MemDims[0] = a_NumDim1;
  MemDims[1] = a_NumNew;
  MemspaceId = H5Screate_simple(2, MemDims, NULL);
  if (MemspaceId < 0) {
    free(OldDims);
    free(OldMaxDims);
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    H5Sclose(FilespaceId);
    return status;
  }

  /* add new information */
  offset[0] = 0;
  offset[1] = OldDims[1];
  count[0] = MemDims[0];
  count[1] = MemDims[1];

  status = H5Sselect_hyperslab(FilespaceId, H5S_SELECT_SET, offset, NULL, 
                               count, NULL);
  if (status >= 0) {
    /* write the data */
    status = H5Dwrite(DatasetId, H5T_NATIVE_DOUBLE, MemspaceId, FilespaceId,
                      H5P_DEFAULT, a_Array);
  }

  free(OldDims);
  free(OldMaxDims);
  H5Sclose(MemspaceId);
  H5Dclose(DatasetId);
  H5Sclose(DataspaceId);
  H5Sclose(FilespaceId);

  return status;
#else
  return -1;
#endif
} /* xfpAppendDset2DDoubleSecondDim */
/******************************************************************************
 * FUNCTION  xfpAppendDset3DDoubleFirstDim
 * PURPOSE   
 * NOTES     
 ******************2***********************************************************/
XMDF_API int xfpAppendDset3DDoubleFirstDim (xid a_Id, const char *a_Name, 
                                   hsize_t a_NumNew, hsize_t a_NumDim2,
                                   hsize_t a_NumDim3,
                                   const double *a_Array)
{
#ifndef _WIN64
  herr_t     status = XTRUE;
  hsize_t   *OldDims = NULL, *OldMaxDims = NULL, NewDims[3], MemDims[3];
  H5T_class_t Classtype;
  xid        DatasetId, DataspaceId, DatatypeId, FilespaceId, MemspaceId;
  int        Rank;
  hssize_t   offset[3];
  hsize_t    count[3];
  hssize_t   nMemSelect;
  hssize_t nSelect;

  DatasetId = H5Dopen(a_Id, a_Name);
  if (DatasetId < 0) {
    return DatasetId;
  }

  DataspaceId = H5Dget_space(DatasetId);
  if (DataspaceId < 0) {
    H5Dclose(DatasetId);
    return DataspaceId;
  }

  /* Dataspace must be simple and have a rank of 2 */
  if (H5Sis_simple(DataspaceId) <= 0) {
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    return ERROR_DATASET_INVALID;
  }

  status = xfpGetSimpleDataspaceInfo(DataspaceId, &Rank, &OldDims, &OldMaxDims);
  if (status < 0) {
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    return status;
  }

  /* This function is only for 3D arrays and MaxDims must be larger */
  /* than the old dimension plus the new values in both directions */
  if (Rank != 3) {
    free(OldDims);
    free(OldMaxDims);
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    return ERROR_DATASET_INVALID;
  }
  NewDims[0] = OldDims[0] + a_NumNew;
  NewDims[1] = OldDims[1];
  NewDims[2] = OldDims[2];
  if (OldMaxDims[0] < NewDims[0] || OldDims[1] != a_NumDim2 ||
      OldDims[2] != a_NumDim3) {
    free(OldDims);
    free(OldMaxDims);
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    return ERROR_DATASET_INVALID;
  }

  /* Make sure the datatype is a float datatype */
  DatatypeId = H5Dget_type(DatasetId);
  if (DatatypeId < 0) {
    free(OldDims);
    free(OldMaxDims);
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    return DatatypeId;
  }
  Classtype = H5Tget_class(DatatypeId);
  H5Tclose(DatatypeId);
  if (Classtype != H5T_FLOAT) {
    free(OldDims);
    free(OldMaxDims);
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    return ERROR_DATASET_INVALID;
  }

  /* extend the dataset */
  status = H5Dextend(DatasetId, NewDims);
  if (status < 0) {
    free(OldDims);
    free(OldMaxDims);
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    return ERROR_DATASET_INVALID;
  }

  /* select a hyperslab */
  FilespaceId = H5Dget_space(DatasetId);
  if (FilespaceId < 0) {
    free(OldDims);
    free(OldMaxDims);
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    return FilespaceId;
  }

  /* create the memory dataspace */
  MemDims[0] = a_NumNew;
  MemDims[1] = a_NumDim2;
  MemDims[2] = a_NumDim3;
  MemspaceId = H5Screate_simple(3, MemDims, NULL);
  if (MemspaceId < 0) {
    free(OldDims);
    free(OldMaxDims);
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    H5Sclose(FilespaceId);
    return status;
  }
  nMemSelect = H5Sget_select_npoints(MemspaceId);

  /* add extras */
  offset[0] = OldDims[0];
  offset[1] = 0;
  offset[2] = 0;
  count[0] = MemDims[0];
  count[1] = MemDims[1];
  count[2] = MemDims[2];

  status = H5Sselect_hyperslab(FilespaceId, H5S_SELECT_SET, offset, NULL, 
                                   count, NULL);
  nSelect = H5Sget_select_npoints(FilespaceId);

  if (status >= 0) {
    /* write the data */
    status = H5Dwrite(DatasetId, H5T_NATIVE_DOUBLE, MemspaceId, FilespaceId,
                      H5P_DEFAULT, a_Array);
  }

  free(OldDims);
  free(OldMaxDims);
  H5Sclose(MemspaceId);
  H5Dclose(DatasetId);
  H5Sclose(DataspaceId);
  H5Sclose(FilespaceId);

  return status;
#else
  return -1;
#endif
} /* xfpAppendDset3DDoubleFirstDim */
/******************************************************************************
 * FUNCTION  xfpAppendDset3DFloatFirstDim
 * PURPOSE   
 * NOTES     
 ******************2***********************************************************/
XMDF_API int xfpAppendDset3DFloatFirstDim (xid a_Id, const char *a_Name, 
                                   hsize_t a_NumNew, hsize_t a_NumDim2,
                                   hsize_t a_NumDim3,
                                   const float *a_Array)
{
#ifndef _WIN64
  herr_t     status = XTRUE;
  hsize_t   *OldDims = NULL, *OldMaxDims = NULL, NewDims[3], MemDims[3];
  H5T_class_t Classtype;
  xid        DatasetId, DataspaceId, DatatypeId, FilespaceId, MemspaceId;
  int        Rank;
  hssize_t   offset[3];
  hsize_t    count[3];
  hssize_t   nMemSelect;
  hssize_t nSelect;

  DatasetId = H5Dopen(a_Id, a_Name);
  if (DatasetId < 0) {
    return DatasetId;
  }

  DataspaceId = H5Dget_space(DatasetId);
  if (DataspaceId < 0) {
    H5Dclose(DatasetId);
    return DataspaceId;
  }

  /* Dataspace must be simple and have a rank of 2 */
  if (H5Sis_simple(DataspaceId) <= 0) {
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    return ERROR_DATASET_INVALID;
  }

  status = xfpGetSimpleDataspaceInfo(DataspaceId, &Rank, &OldDims, &OldMaxDims);
  if (status < 0) {
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    return status;
  }

  /* This function is only for 3D arrays and MaxDims must be larger */
  /* than the old dimension plus the new values in both directions */
  if (Rank != 3) {
    free(OldDims);
    free(OldMaxDims);
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    return ERROR_DATASET_INVALID;
  }
  NewDims[0] = OldDims[0] + a_NumNew;
  NewDims[1] = OldDims[1];
  NewDims[2] = OldDims[2];
  if (OldMaxDims[0] < NewDims[0] || OldDims[1] != a_NumDim2 ||
      OldDims[2] != a_NumDim3) {
    free(OldDims);
    free(OldMaxDims);
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    return ERROR_DATASET_INVALID;
  }

  /* Make sure the datatype is a float datatype */
  DatatypeId = H5Dget_type(DatasetId);
  if (DatatypeId < 0) {
    free(OldDims);
    free(OldMaxDims);
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    return DatatypeId;
  }
  Classtype = H5Tget_class(DatatypeId);
  H5Tclose(DatatypeId);
  if (Classtype != H5T_FLOAT) {
    free(OldDims);
    free(OldMaxDims);
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    return ERROR_DATASET_INVALID;
  }

  /* extend the dataset */
  status = H5Dextend(DatasetId, NewDims);
  if (status < 0) {
    free(OldDims);
    free(OldMaxDims);
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    return ERROR_DATASET_INVALID;
  }

  /* select a hyperslab */
  FilespaceId = H5Dget_space(DatasetId);
  if (FilespaceId < 0) {
    free(OldDims);
    free(OldMaxDims);
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    return FilespaceId;
  }

  /* create the memory dataspace */
  MemDims[0] = a_NumNew;
  MemDims[1] = a_NumDim2;
  MemDims[2] = a_NumDim3;
  MemspaceId = H5Screate_simple(3, MemDims, NULL);
  if (MemspaceId < 0) {
    free(OldDims);
    free(OldMaxDims);
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    H5Sclose(FilespaceId);
    return status;
  }
  nMemSelect = H5Sget_select_npoints(MemspaceId);

  /* add extras */
  offset[0] = OldDims[0];
  offset[1] = 0;
  offset[2] = 0;
  count[0] = MemDims[0];
  count[1] = MemDims[1];
  count[2] = MemDims[2];

  status = H5Sselect_hyperslab(FilespaceId, H5S_SELECT_SET, offset, NULL, 
                                   count, NULL);
  nSelect = H5Sget_select_npoints(FilespaceId);

  if (status >= 0) {
    /* write the data */
    status = H5Dwrite(DatasetId, H5T_NATIVE_FLOAT, MemspaceId, FilespaceId,
                      H5P_DEFAULT, a_Array);
  }

  free(OldDims);
  free(OldMaxDims);
  H5Sclose(MemspaceId);
  H5Dclose(DatasetId);
  H5Sclose(DataspaceId);
  H5Sclose(FilespaceId);

  return status;
#else
  return -1;
#endif
} /* xfpAppendDset3DFloatFirstDim */
/******************************************************************************
 * FUNCTION  xfpAppendDset2DUCharFirstDim
 * PURPOSE   
 * NOTES     
 ******************2***********************************************************/
XMDF_API int xfpAppendDset2DUCharFirstDim (xid a_Id, const char *a_Name, 
                                   hsize_t a_NumNew, hsize_t a_NumDim2,
                                   const unsigned char *a_Array)
{
#ifndef _WIN64
  herr_t     status = XTRUE;
  hsize_t   *OldDims = NULL, *OldMaxDims = NULL, NewDims[2], MemDims[2];
  H5T_class_t Classtype;
  xid        DatasetId, DataspaceId, DatatypeId, FilespaceId, MemspaceId;
  int        Rank;
  hssize_t   offset[2];
  hsize_t    count[2];

  DatasetId = H5Dopen(a_Id, a_Name);
  if (DatasetId < 0) {
    return DatasetId;
  }

  DataspaceId = H5Dget_space(DatasetId);
  if (DataspaceId < 0) {
    H5Dclose(DatasetId);
    return DataspaceId;
  }

  /* Dataspace must be simple and have a rank of 2 */
  if (H5Sis_simple(DataspaceId) <= 0) {
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    return ERROR_DATASET_INVALID;
  }

  status = xfpGetSimpleDataspaceInfo(DataspaceId, &Rank, &OldDims, &OldMaxDims);
  if (status < 0) {
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    return status;
  }

  /* This function is only for 2D arrays and MaxDims must be larger */
  /* than the old dimension plus the new values in both directions */
  if (Rank != 2) {
    free(OldDims);
    free(OldMaxDims);
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    return ERROR_DATASET_INVALID;
  }
  NewDims[0] = OldDims[0] + a_NumNew;
  NewDims[1] = OldDims[1];
  if (OldMaxDims[0] < NewDims[0] || OldDims[1] != a_NumDim2) {
    free(OldDims);
    free(OldMaxDims);
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    return ERROR_DATASET_INVALID;
  }

  /* Make sure the datatype is an integer datatype (unsigned char are stored as ints) */
  DatatypeId = H5Dget_type(DatasetId);
  if (DatatypeId < 0) {
    free(OldDims);
    free(OldMaxDims);
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    return DatatypeId;
  }
  Classtype = H5Tget_class(DatatypeId);
  H5Tclose(DatatypeId);
  if (Classtype != H5T_INTEGER) {
    free(OldDims);
    free(OldMaxDims);
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    return ERROR_DATASET_INVALID;
  }

  /* extend the dataset */
  status = H5Dextend(DatasetId, NewDims);
  if (status < 0) {
    free(OldDims);
    free(OldMaxDims);
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    return ERROR_DATASET_INVALID;
  }

  /* select a hyperslab */
  FilespaceId = H5Dget_space(DatasetId);
  if (FilespaceId < 0) {
    free(OldDims);
    free(OldMaxDims);
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    return FilespaceId;
  }

  /* create the memory dataspace */
  MemDims[0] = a_NumNew;
  MemDims[1] = a_NumDim2;
  MemspaceId = H5Screate_simple(2, MemDims, NULL);
  if (MemspaceId < 0) {
    free(OldDims);
    free(OldMaxDims);
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    H5Sclose(FilespaceId);
    return status;
  }

  /* add extras */
  offset[0] = OldDims[0];
  offset[1] = 0;
  count[0] = MemDims[0];
  count[1] = MemDims[1];

  status = H5Sselect_hyperslab(FilespaceId, H5S_SELECT_SET, offset, NULL, 
                               count, NULL);
  if (status >= 0) {
    /* write the data */
    status = H5Dwrite(DatasetId, H5T_NATIVE_UCHAR, MemspaceId, FilespaceId,
                      H5P_DEFAULT, a_Array);
  }

  free(OldDims);
  free(OldMaxDims);
  H5Sclose(MemspaceId);
  H5Dclose(DatasetId);
  H5Sclose(DataspaceId);
  H5Sclose(FilespaceId);

  return status;
#else
  return -1;
#endif
} /* xfpAset2DUCharFirstDim */
/******************************************************************************
 * FUNCTION  xfpReadDatasetString
 * PURPOSE   Function designed to read a dataset consisting of a single string
 * NOTES     the string is allocated to the length of the dataset using malloc
 ******************2***********************************************************/
XMDF_API int xfpReadDatasetString (xid a_Loc, const char *a_Name, char **a_Str)
{
#ifndef _WIN64
  hid_t   DsetId, DataspaceId, DataTypeId, StringTypeId;
  herr_t  status;
  htri_t  IsSimple;
  int     Rank;
  size_t  StrSize;
  hsize_t *Dims = NULL, *MaxDims = NULL;

  /* open the datatype and make sure it is a string type */
  DsetId = H5Dopen(a_Loc, a_Name);
  if (DsetId < 0) {
   return DsetId;
  }
  DataTypeId = H5Dget_type(DsetId);
  StrSize = H5Tget_size(DataTypeId);
  /* Get the dataspace */
  DataspaceId = H5Dget_space(DsetId);

  /* Dataspace must be simply defined */
  IsSimple = H5Sis_simple(DataspaceId);
  if (!IsSimple) {
    H5Tclose(DataTypeId);
    H5Sclose(DataspaceId);
    H5Dclose(DsetId);
    return ERROR_ATTRIBUTE_NOT_SUPPORTED;
  }

  /* Get the rank and dimensions for simple dataspace */
  status = xfpGetSimpleDataspaceInfo(DataspaceId, &Rank, &Dims, &MaxDims);
  if (status < 0 || Rank != 1) {
    H5Tclose(DataTypeId);
    H5Sclose(DataspaceId);
    H5Aclose(DsetId);
    return ERROR_ATTRIBUTE_NOT_SUPPORTED;
  }

  /* Make sure we have a single string */
  if (Dims[0] != 1) {
    H5Tclose(DataTypeId);
    H5Sclose(DataspaceId);
    H5Aclose(DsetId);
    return ERROR_DATASET_SIZE_INCORRECT;
  }
  
  /* create a datatype for the type of string we are using in memory */
  StringTypeId = H5Tcopy(H5T_C_S1);
  if (StringTypeId < 0) {
    return 0;
  }
  H5Tset_strpad(StringTypeId, H5T_STR_NULLTERM);
  H5Tset_size(StringTypeId, StrSize); 

  /* allocate the array */
  *a_Str = (char *)malloc(StrSize*sizeof(char ));

  status = H5Dread(DsetId, StringTypeId, H5S_ALL, H5S_ALL, H5P_DEFAULT, *a_Str);

  /* close resources */
  H5Tclose(DataTypeId);
  H5Sclose(DataspaceId);
  H5Dclose(DsetId);

  if (Dims) {
    free(Dims);
  }
  if (MaxDims) {
    free(MaxDims);
  }

  return status;
#else
  return -1;
#endif
} /* xfpReadDatasetString */
/******************************************************************************
 * FUNCTION  xfpWriteDatasetString
 * PURPOSE   Create a dataset that is a single string
 * NOTES     
 ******************2***********************************************************/
XMDF_API int xfpWriteDatasetString (xid a_Loc, const char *a_Name, 
                             const char *a_Str)
{
#ifndef _WIN64
  hid_t   DsetId, SpaceId;
  hid_t   StringType;
  hsize_t Dims;
  herr_t  status;

  /* Create the string type */
  StringType = H5Tcopy(H5T_C_S1);
  if (StringType < 0) {
    return StringType;
  }
  H5Tset_strpad(StringType, H5T_STR_NULLTERM);

  /* Set the length of the string datatype */
  status = H5Tset_size(StringType, strlen(a_Str) + 1);

  /* Create the dataspace; */
  Dims = 1;
  SpaceId = H5Screate_simple(1, &Dims, &Dims);
    
  DsetId = H5Dcreate(a_Loc, a_Name, StringType, SpaceId, H5P_DEFAULT);
  if (DsetId < 0) {
    H5Tclose(StringType);
    H5Sclose(SpaceId);
    return DsetId;
  }

  status = H5Dwrite(DsetId, StringType, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                                                                      a_Str);

  /* close resources */
  H5Sclose(SpaceId);
  H5Dclose(DsetId);
  H5Tclose(StringType);

  return status;
#else
  return -1;
#endif
} /* xfpWriteDatasetString */
/******************************************************************************
 * FUNCTION  xfpWriteDatasetStrings
 * PURPOSE   Create a dataset that is a multiple strings
 * NOTES     
 ******************2***********************************************************/
XMDF_API int xfpWriteDatasetStrings (xid a_Loc, const char *a_Name, hsize_t a_Number,
                            size_t a_StringLength, const char *a_Str)
{
#ifndef _WIN64
  hid_t   DsetId, SpaceId, PropId;
  hid_t   StringType;
  hsize_t Dims;
  herr_t  status;

  /* Create the string type */
  StringType = H5Tcopy(H5T_C_S1);
  if (StringType < 0) {
    return StringType;
  }
  H5Tset_strpad(StringType, H5T_STR_NULLTERM);

  /* Set the length of the string datatype add one space for NULL terminator */
  status = H5Tset_size(StringType, a_StringLength + (size_t)1);

  /* Create the dataspace; */
  Dims = a_Number;
  SpaceId = H5Screate_simple(1, &a_Number, &a_Number);

  PropId = H5Pcreate(H5P_DATASET_CREATE);
  H5Pset_chunk(PropId, 1, &a_Number);
  H5Pset_deflate(PropId, 9);
  DsetId = H5Dcreate(a_Loc, a_Name, StringType, SpaceId, PropId);
  if (DsetId < 0) {
    H5Sclose(SpaceId);
    return DsetId;
  }

  status = H5Dwrite(DsetId, StringType, H5S_ALL, H5S_ALL, H5P_DEFAULT, a_Str);

  /* close resources */
  H5Sclose(SpaceId);
  H5Dclose(DsetId);
  H5Tclose(StringType);

  return status;
#else
  return -1;
#endif
} /* xfpWriteDatasetStrings */
/******************************************************************************
 * FUNCTION  xifDoesAttributeExist
 * PURPOSE   Determine if an attribute exists with a given name
 * NOTES     
 ******************2***********************************************************/
XMDF_API int xfpDoesAttributeExist (xid a_Loc, const char *a_Name,
                                  xmbool *a_bExists)
{
#ifndef _WIN64
  xid    AttId;
  herr_t (*old_func)(void*);
  void    *old_client_data;

  H5Eget_auto(&old_func, &old_client_data);

  /* Turn off error handling */
  H5Eset_auto(NULL, NULL);

  /* try to open the dataset */
  AttId = H5Aopen_name(a_Loc, a_Name);

  /* Restore previous error handler */
  H5Eset_auto(old_func, old_client_data);

  if (AttId > 0) {
    H5Aclose(AttId);
    *a_bExists = XTRUE;
  }
  else {
    *a_bExists = XFALSE;
  }
  return XTRUE;
#else
  return XFALSE;
#endif
} /* xfpDoesAttributeExist */
/******************************************************************************
 * FUNCTION  xfpReadAttributeString
 * PURPOSE   Function designed to read a single attribute string
 * NOTES     the string is allocated to the length of the attribute using malloc
 ******************2***********************************************************/
XMDF_API int xfpReadAttributeString (xid a_Loc, const char *a_Name, int a_Num, char **a_Str)
{
#ifndef _WIN64
  hid_t   AttId, DataspaceId, DataTypeId, StringTypeId;
  herr_t  status;
  htri_t  IsSimple;
  int     Rank;
  size_t  StrSize;
  hsize_t *Dims = NULL, *MaxDims = NULL;
  H5T_class_t DatasetClass;
  herr_t (*old_func)(void*);
  void  *old_client_data;

  /* make our own string type */
  /* Create the string type */
  StringTypeId = H5Tcopy(H5T_C_S1);
  if (StringTypeId < 0) {
    return StringTypeId;
  }
  H5Tset_strpad(StringTypeId, H5T_STR_NULLTERM);

  H5Eget_auto(&old_func, &old_client_data);

  /* Turn off error handling */
  H5Eset_auto(NULL, NULL);

  AttId = H5Aopen_name(a_Loc, a_Name);

  /* Restore previous error handler */
  H5Eset_auto(old_func, old_client_data);

  /* return a negative number if we did not open an attribute */
  if (AttId < 0) {
    return AttId;
  }

  /* open the datatype and make sure it is a string type */
  DataTypeId = H5Aget_type(AttId);
  if (DataTypeId < 0) {
    H5Aclose(AttId);
    return DataTypeId;
  }

  /* Make sure it is a string type */
  DatasetClass = H5Tget_class(DataTypeId);
  if (DatasetClass != H5T_STRING) {
    H5Tclose(DataTypeId);
    return ERROR_INCORRECT_DATATYPE;
  }

  /* Get the size of the string */
  StrSize = H5Tget_size(DataTypeId);

    /* Set the length of the string datatype */
  status = H5Tset_size(StringTypeId, StrSize + 1);

  /* Get the dataspace */
  DataspaceId = H5Aget_space(AttId);

  /* Dataspace must be simply defined */
  IsSimple = H5Sis_simple(DataspaceId);
  if (!IsSimple) {
    H5Tclose(DataTypeId);
    H5Sclose(DataspaceId);
    H5Aclose(AttId);
    return ERROR_ATTRIBUTE_NOT_SUPPORTED;
  }

  /* Get the rank and dimensions for simple dataspace */
  status = xfpGetSimpleDataspaceInfo(DataspaceId, &Rank, &Dims, &MaxDims);
  if (status < 0 || Rank != 1) {
    H5Tclose(DataTypeId);
    H5Sclose(DataspaceId);
    H5Aclose(AttId);
    return ERROR_ATTRIBUTE_NOT_SUPPORTED;
  }

  /* Make sure we the correct number */
  if (Dims[0] != a_Num) {
    H5Sclose(DataspaceId);
    H5Aclose(AttId);
    H5Tclose(DataTypeId);
    return ERROR_DATASET_SIZE_INCORRECT;
  }
 
  /* allocate the array */
  *a_Str = (char *)malloc((StrSize+1)*a_Num*sizeof(char ));

  status = H5Aread(AttId, StringTypeId, *a_Str);

  /* close resources */
  H5Tclose(DataTypeId);
  H5Sclose(DataspaceId);
  H5Aclose(AttId);
  free(Dims);
  free(MaxDims);

  return status;
#else
  return -1;
#endif
} /* xfpReadAttributeString */
/******************************************************************************
 * FUNCTION  xfpReadAttributeInt
 * PURPOSE   Function designed to read a single integer attribute
 * NOTES     
 ******************2***********************************************************/
XMDF_API int xfpReadAttributeInt (xid a_Loc, const char *a_Name, int a_Num, int *a_val)
{
#ifndef _WIN64
  hid_t   AttId, DataspaceId, DataTypeId;
  herr_t  status;
  htri_t  IsSimple;
  int     Rank;
  hsize_t *Dims = NULL, *MaxDims = NULL;
  H5T_class_t DatasetClass;

  AttId = H5Aopen_name(a_Loc, a_Name);
  if (AttId < 0) {
    return AttId;
  }

  /* open the datatype and make sure it is an integer type */
  DataTypeId = H5Aget_type(AttId);
  if (DataTypeId < 0) {
    H5Aclose(AttId);
    return DataTypeId;
  }

  /* Make sure it is an integer type */
  DatasetClass = H5Tget_class(DataTypeId);
  if (DatasetClass != H5T_INTEGER) {
    H5Tclose(DataTypeId);
    return ERROR_INCORRECT_DATATYPE;
  }

  /* Get the dataspace */
  DataspaceId = H5Aget_space(AttId);

  /* Dataspace must be simply defined */
  IsSimple = H5Sis_simple(DataspaceId);
  if (!IsSimple) {
    H5Tclose(DataTypeId);
    H5Sclose(DataspaceId);
    H5Aclose(AttId);
    return ERROR_ATTRIBUTE_NOT_SUPPORTED;
  }

  /* Get the rank and dimensions for simple dataspace */
  status = xfpGetSimpleDataspaceInfo(DataspaceId, &Rank, &Dims, &MaxDims);
  if (status < 0 || Rank != 1) {
    H5Tclose(DataTypeId);
    H5Sclose(DataspaceId);
    H5Aclose(AttId);
    return ERROR_ATTRIBUTE_NOT_SUPPORTED;
  }

  /* Make sure we have a single integer */
  if (Dims[0] != a_Num) {
    H5Sclose(DataspaceId);
    H5Aclose(AttId);
    H5Tclose(DataTypeId);
    return ERROR_DATASET_SIZE_INCORRECT;
  }
 
  status = H5Aread(AttId, H5T_NATIVE_INT, a_val);

  /* close resources */
  H5Tclose(DataTypeId);
  H5Sclose(DataspaceId);
  H5Aclose(AttId);
  free(Dims);
  free(MaxDims);

  return status;
#else
  return -1;
#endif
} /* xfpReadAttributeInt */
/******************************************************************************
 * FUNCTION  xfpReadAttributeFloat
 * PURPOSE   Function designed to read a single float attribute
 * NOTES     
 ******************2***********************************************************/
XMDF_API int xfpReadAttributeFloat (xid a_Loc, const char * a_Name, int a_Num, float *a_val)
{
#ifndef _WIN64
  hid_t   AttId, DataspaceId, DataTypeId;
  herr_t  status;
  htri_t  IsSimple;
  int     Rank;
  hsize_t *Dims = NULL, *MaxDims = NULL;
  H5T_class_t DatasetClass;
  xmbool   exist;

  xfpDoesAttributeExist(a_Loc, a_Name, &exist);
  if (!exist) {
    return -1;
  }
  AttId = H5Aopen_name(a_Loc, a_Name);
  if (AttId < 0) {
    return AttId;
  }

  /* open the datatype and make sure it is a string type */
  DataTypeId = H5Aget_type(AttId);
  if (DataTypeId < 0) {
    H5Aclose(AttId);
    return DataTypeId;
  }

  /* Make sure it is a float type */
  DatasetClass = H5Tget_class(DataTypeId);
  if (DatasetClass != H5T_FLOAT) {
    H5Tclose(DataTypeId);
    return ERROR_INCORRECT_DATATYPE;
  }

  /* Get the dataspace */
  DataspaceId = H5Aget_space(AttId);

  /* Dataspace must be simply defined */
  IsSimple = H5Sis_simple(DataspaceId);
  if (!IsSimple) {
    H5Tclose(DataTypeId);
    H5Sclose(DataspaceId);
    H5Aclose(AttId);
    return ERROR_ATTRIBUTE_NOT_SUPPORTED;
  }

  /* Get the rank and dimensions for simple dataspace */
  status = xfpGetSimpleDataspaceInfo(DataspaceId, &Rank, &Dims, &MaxDims);
  if (status < 0 || Rank != 1) {
    H5Tclose(DataTypeId);
    H5Sclose(DataspaceId);
    H5Aclose(AttId);
    return ERROR_ATTRIBUTE_NOT_SUPPORTED;
  }

  /* Make sure we have a single integer */
  if (Dims[0] != a_Num) {
    H5Sclose(DataspaceId);
    H5Aclose(AttId);
    H5Tclose(DataTypeId);
    return ERROR_DATASET_SIZE_INCORRECT;
  }
 
  status = H5Aread(AttId, H5T_NATIVE_FLOAT, a_val);

  /* close resources */
  H5Tclose(DataTypeId);
  H5Sclose(DataspaceId);
  H5Aclose(AttId);
  free(Dims);
  free(MaxDims);

  return status;
#else
  return -1;
#endif
} /* xfpReadAttributeFloat */
/******************************************************************************
 * FUNCTION  xfpReadAttributeDouble
 * PURPOSE   Function designed to read a single double attribute
 * NOTES     
 ******************2***********************************************************/
XMDF_API int xfpReadAttributeDouble (xid a_Loc, const char *a_Name, int a_Num, double *a_val)
{
#ifndef _WIN64
  hid_t   AttId, DataspaceId, DataTypeId;
  herr_t  status;
  htri_t  IsSimple;
  int     Rank;
  hsize_t *Dims = NULL, *MaxDims = NULL;
  H5T_class_t DatasetClass;

  AttId = H5Aopen_name(a_Loc, a_Name);
  if (AttId < 0) {
    return AttId;
  }

  /* open the datatype and make sure it is a string type */
  DataTypeId = H5Aget_type(AttId);
  if (DataTypeId < 0) {
    H5Aclose(AttId);
    return DataTypeId;
  }

  /* Make sure it is a float type */
  DatasetClass = H5Tget_class(DataTypeId);
  if (DatasetClass != H5T_FLOAT) {
    H5Tclose(DataTypeId);
    return ERROR_INCORRECT_DATATYPE;
  }

  /* Get the dataspace */
  DataspaceId = H5Aget_space(AttId);

  /* Dataspace must be simply defined */
  IsSimple = H5Sis_simple(DataspaceId);
  if (!IsSimple) {
    H5Tclose(DataTypeId);
    H5Sclose(DataspaceId);
    H5Aclose(AttId);
    return ERROR_ATTRIBUTE_NOT_SUPPORTED;
  }

  /* Get the rank and dimensions for simple dataspace */
  status = xfpGetSimpleDataspaceInfo(DataspaceId, &Rank, &Dims, &MaxDims);
  if (status < 0 || Rank != 1) {
    H5Tclose(DataTypeId);
    H5Sclose(DataspaceId);
    H5Aclose(AttId);
    return ERROR_ATTRIBUTE_NOT_SUPPORTED;
  }

  /* Make sure we have a single value */
  if (Dims[0] != a_Num) {
    H5Sclose(DataspaceId);
    H5Aclose(AttId);
    H5Tclose(DataTypeId);
    return ERROR_DATASET_SIZE_INCORRECT;
  }
 
  status = H5Aread(AttId, H5T_NATIVE_DOUBLE, a_val);

  /* close resources */
  H5Tclose(DataTypeId);
  H5Sclose(DataspaceId);
  H5Aclose(AttId);
  free(Dims);
  free(MaxDims);

  return status;
#else
  return -1;
#endif
} /* xfpReadAttributeDouble */
/******************************************************************************
 * FUNCTION  xfpWriteAttributeString
 * PURPOSE   
 * NOTES     
 ******************2***********************************************************/
XMDF_API int xfpWriteAttributeString (xid a_Loc, const char * a_Name, 
                             const char * a_Str)
{
#ifndef _WIN64
  hid_t   AttId, SpaceId;
  hid_t   StringType;
  hsize_t Dims;
  herr_t  status;
  int     length;

  /* Create the string type */
  StringType = H5Tcopy(H5T_C_S1);
  if (StringType < 0) {
    return StringType;
  }
  H5Tset_strpad(StringType, H5T_STR_NULLTERM);

  /* Set the length of the string datatype */
  length = strlen(a_Str) + 1;
  status = H5Tset_size(StringType, length);

  /* Create the dataspace; */
  Dims = 1;
  SpaceId = H5Screate_simple(1, &Dims, &Dims);
    
  AttId = H5Acreate(a_Loc, a_Name, StringType, SpaceId, H5P_DEFAULT);
  if (AttId < 0) {
    H5Sclose(SpaceId);
    return AttId;
  }

  status = H5Awrite(AttId, StringType, a_Str);

  /* close resources */
  H5Sclose(SpaceId);
  H5Aclose(AttId);
  H5Tclose(StringType);

  return status;
#else
  return -1;
#endif
} /* xfpWriteAttributeString */
/******************************************************************************
 * FUNCTION  xfpWriteAttributeInt
 * PURPOSE   
 * NOTES     
 ******************2***********************************************************/
XMDF_API int xfpWriteAttributeInt (xid a_Loc, const char *a_Name, int a_Number,
                                 int *a_val)
{
#ifndef _WIN64
  hid_t   AttId, SpaceId;
  hsize_t Dims;
  herr_t  status;
  xmbool  exist;

  /* Create the dataspace; */
  Dims = a_Number;
  SpaceId = H5Screate_simple(1, &Dims, &Dims);
  AttId = -1;

  xfpDoesAttributeExist(a_Loc, a_Name, &exist);
  if (exist) {
    AttId = H5Aopen_name(a_Loc, a_Name);
  }
  if (AttId < 0) {
    AttId = H5Acreate(a_Loc, a_Name, H5T_NATIVE_INT, SpaceId, H5P_DEFAULT);
    if (AttId < 0) {
      H5Sclose(SpaceId);
      return AttId;
    }
  }

  status = H5Awrite(AttId, H5T_NATIVE_INT, a_val);

  /* close resources */
  H5Sclose(SpaceId);
  H5Aclose(AttId);

  return status;
#else
  return -1;
#endif
} /* xfpWriteAttributeInt */
/******************************************************************************
 * FUNCTION  xfpWriteAttributeFloat
 * PURPOSE   
 * NOTES     
 ******************2***********************************************************/
XMDF_API int xfpWriteAttributeFloat (xid a_Loc, const char *a_Name, 
                                   int a_Number, float *a_val)
{
#ifndef _WIN64
  hid_t   AttId, SpaceId;
  hsize_t Dims;
  herr_t  status;

  /* Create the dataspace; */
  Dims = a_Number;
  SpaceId = H5Screate_simple(1, &Dims, &Dims);
    
  AttId = H5Acreate(a_Loc, a_Name, fg_Floattype, SpaceId, H5P_DEFAULT);
  if (AttId < 0) {
    H5Sclose(SpaceId);
    return AttId;
  }

  status = H5Awrite(AttId, H5T_NATIVE_FLOAT, a_val);

  /* close resources */
  H5Sclose(SpaceId);
  H5Aclose(AttId);

  return status;
#else
  return -1;
#endif
} /* xfpWriteAttributeFloat */
/******************************************************************************
 * FUNCTION  xfpWriteAttributeDouble
 * PURPOSE   
 * NOTES     
 ******************2***********************************************************/
XMDF_API int xfpWriteAttributeDouble (xid a_Loc, const char *a_Name, 
                                    int a_Number, double *a_val)
{
#ifndef _WIN64
  hid_t   AttId, SpaceId;
  hsize_t Dims;
  herr_t  status;

  /* Create the dataspace; */
  Dims = a_Number;
  SpaceId = H5Screate_simple(1, &Dims, &Dims);
    
  AttId = H5Acreate(a_Loc, a_Name, fg_Doubletype, SpaceId, H5P_DEFAULT);
  if (AttId < 0) {
    H5Sclose(SpaceId);
    return AttId;
  }

  status = H5Awrite(AttId, H5T_NATIVE_DOUBLE, a_val);

  /* close resources */
  H5Sclose(SpaceId);
  H5Aclose(AttId);

  return status;
#else
  return -1;
#endif
} /* xfpWriteAttributeDouble */
/******************************************************************************
 * FUNCTION  xfpGetSimpleDataspaceInfoFromName
 * PURPOSE   Get the rank, size, and max size for a simple dataspace
 * NOTES     
 ******************2***********************************************************/
XMDF_API int xfpGetSimpleDataspaceInfoFromName (xid a_Loc, const char *a_Name,
                          int *a_Rank, hsize_t **a_Dims, hsize_t **a_MaxDims)
{
#ifndef _WIN64
  hid_t    Id, SpaceId;
  int      rval;

  Id = H5Dopen(a_Loc, a_Name);
  if (Id < 0) {
    return -1;
  }
  SpaceId = H5Dget_space(Id);
  if (SpaceId < 0) {
    return -1;
  }
  
  rval = xfpGetSimpleDataspaceInfo(SpaceId, a_Rank, a_Dims, a_MaxDims);

  H5Sclose(SpaceId);
  H5Dclose(Id);
  return rval;
#else
  return -1;
#endif
} /* xfpGetSimpleDataspaceInfoFromName */
/******************************************************************************
 * FUNCTION  xfpGetSimpleDataspaceInfo
 * PURPOSE   Get the rank, size, and max size for a simple dataspace
 * NOTES     
 ******************2***********************************************************/
XMDF_API int xfpGetSimpleDataspaceInfo (xid a_Id, int *a_Rank, hsize_t **a_Dims, 
                               hsize_t **a_MaxDims)
{
#ifndef _WIN64
  *a_Rank = H5Sget_simple_extent_ndims(a_Id);
  if (*a_Rank <= 0) {
    return *a_Rank;
  }  

  /* allocate the arrays */
  *a_Dims = (hsize_t *)malloc(*a_Rank * sizeof(hsize_t));
  *a_MaxDims = (hsize_t *)malloc(*a_Rank * sizeof(hsize_t));

  /* Get the dimensions */
  return H5Sget_simple_extent_dims(a_Id, *a_Dims, *a_MaxDims);
#else
  return -1;
#endif
} /* xfpGetSimpleDataspaceInfo */
/******************************************************************************
 * FUNCTION  xfpGetSimpleDataspaceInfo
 * PURPOSE   free the space allocated by xfpGetSimpleDataspaceInfo
 * NOTES     
 ******************2***********************************************************/
XMDF_API void xfpDestroySimpleDataspaceInfo (hsize_t **a_Dims,
                                             hsize_t **a_MaxDims)
{
  if (a_Dims && *a_Dims)
    free(*a_Dims);
  if (a_MaxDims && *a_MaxDims)
    free(*a_MaxDims);
} /* xfpDestroySimpleDataspaceInfo */
/******************************************************************************
 * FUNCTION  xfpCreateGroup
 * PURPOSE   create the group starting from the file or group a_Id
 * NOTES     Each group stores the absolute path until HDF5 gives us a function
 *           to obtain it.
 *           Intermediate groups are indicated using the slash symbol (/)
 *           Intermediate groups do not to be created 
 *           Returns the group id unless there was an error and then returns a
 *             negative number
 ******************2***********************************************************/
XMDF_API int xfpCreateGroup (xid a_Id, const char *a_Path, const char *a_GroupType)
{
#ifndef _WIN64
  int    status;
  int    iStringSize;
  char  *strParentPath = NULL, *strCompletePath = NULL;
  const char *strGroups = NULL;
  char  *strNextGroup = NULL, *strExistingGrouptype;
  size_t EndOfGroup;
  hid_t  ParentGroup = NONE, Group = NONE;
  herr_t (*old_func)(void*);
  void  *old_client_data;
  xmbool  bNewGroup = XFALSE; /* See if we create a group or use an existing group */

  H5Eget_auto(&old_func, &old_client_data);

  /* Turn off error handling */
  H5Eset_auto(NULL, NULL);

  /* Get the path of the parent group */
  /* Probe. Likely to fail, but that's okay */

  status = xfpGetAbsolutePath(a_Id, &strParentPath);

  /* Restore previous error handler */
  H5Eset_auto(old_func, old_client_data);
  if (status >= 0) {
    /* allocate the string */
    iStringSize = strlen(strParentPath) + strlen(a_Path) + 2;
    strCompletePath = (char *)malloc(iStringSize*sizeof(char));
    /* append the paths together for complete path */
    sprintf(strCompletePath, "%s/%s", strParentPath, a_Path);
  }
  else {
    /* must be the main file complete path is path sent in */
    iStringSize = strlen(a_Path) + 1;
    strCompletePath = (char *)malloc(iStringSize*sizeof(char));
    /* append the paths together for complete path */
    strcpy(strCompletePath, a_Path);
  }

  /* Create the intermediate and final groups */
  strGroups = a_Path;
  strNextGroup = (char *)malloc(strlen(a_Path)+1*sizeof(char));
  /* get the first substring until the end or until we hit a slash / */
  EndOfGroup = strcspn(strGroups, "/");
  strncpy(strNextGroup, strGroups, EndOfGroup);
  strNextGroup[EndOfGroup] = '\0';
  if (EndOfGroup == strlen(strGroups)) {
    strGroups = strGroups + EndOfGroup;
  }
  else {
    strGroups = strGroups + EndOfGroup + 1;
  }
  ParentGroup = a_Id;
  status = 1;
  while (strlen(strNextGroup) > 0 && status >= 1) {
    /* open the group if it exists */
    /* Turn off error handling */
    H5Eset_auto(NULL, NULL);
    Group = H5Gopen(ParentGroup, strNextGroup);
    /* Restore previous error handler */
    H5Eset_auto(old_func, old_client_data);
    if (Group < 0) {
      /* group does not exist create */
      Group = H5Gcreate(ParentGroup, strNextGroup, DEFAULT_GROUP_SIZE);
      status = Group;
      bNewGroup = XTRUE;
    }
    /* close the parent group if it is not the group sent in */
    if (ParentGroup != a_Id) {
      H5Gclose(ParentGroup);
    }
    /* reassign the parent group to be this group and recopy the strings */
    ParentGroup = Group;
    EndOfGroup = strcspn(strGroups, "/");
    strncpy(strNextGroup, strGroups, EndOfGroup);
    strNextGroup[EndOfGroup] = '\0';
    if (strlen(strNextGroup) > 0) {
      if (strGroups[EndOfGroup] == '/') {
        /* move one past the next path marker */
        strGroups = strGroups + EndOfGroup + 1;
      }
      else {
        /* We are done set strGroups to be the NULL at end. */
        strGroups = strGroups + EndOfGroup;
      }
    }
  }
  if (status < 0) {
    if (strCompletePath) {
      free (strCompletePath);
    }
    if (strParentPath) {
      free (strParentPath);
    }
    if (strNextGroup) {
      free (strNextGroup);
    }
    return status;
  }

  /* right now ParentGroup is the group we wanted to create reassign to Group */
  Group = ParentGroup;
  
  if (bNewGroup) {
    /* Add the group type identifier to the final group*/
    status = xfpWriteAttributeString(Group, GROUP_TYPE, a_GroupType);
    if (status < 0) {
      if (strCompletePath) {
        free (strCompletePath);
      }
      if (strParentPath) {
        free (strParentPath);
      }
      if (strNextGroup) {
        free (strNextGroup);
      }
      H5Gclose(Group);
      return status;
    }
  }
  else {
    /* Found an old group make sure the identifier string is the same as we want */
    status = xfpReadAttributeString(Group, GROUP_TYPE, 1, &strExistingGrouptype);
    if (status < 0 || strcmp(strExistingGrouptype, a_GroupType) != 0) {
      /* the group type is not set correctly. */
      if (strCompletePath) {
        free (strCompletePath);
      }
      if (strParentPath) {
        free (strParentPath);
      }
      if (strNextGroup) {
        free (strNextGroup);
      }
      H5Gclose(Group);
      return status;
    }
  }

  /* free resources */
  if (strCompletePath) {
    free (strCompletePath);
  }
  if (strParentPath) {
    free (strParentPath);
  }
  if (strNextGroup) {
    free (strNextGroup);
  }

  return Group;
#else
  return -1;
#endif
} /* xfpCreateGroup */
/******************************************************************************
 * FUNCTION  xfpRecurseCountGroups
 * PURPOSE   called by the HDF5 iterator function to count the number of
 *           groups and the maximum path length
 * NOTES     
 ******************2***********************************************************/
XMDF_API herr_t xfpRecurseCountGroups(xid a_Id, XGroupIteration *a_CountObj)
{
#ifndef _WIN64
  herr_t   status;
  hsize_t  nObjs;
  ssize_t  iNameLength;
  int      i, iObjectType, PathLength;
  char    *Path = NULL, *GroupType = NULL;
  hid_t    Id, ChildId;
  xmbool      bGenericGroup = XTRUE, bRootGroupOpened = XFALSE;
  H5I_type_t   IdType;

  /* see if they are trying to look in the root directory using a file id*/
  IdType = H5Iget_type(a_Id);
  if (IdType == H5I_FILE) {
    /* open the root group*/
    Id = H5Gopen(a_Id, "/");
    bRootGroupOpened = XTRUE;
  }
  else if (IdType == H5I_GROUP) {
    Id = a_Id;
    bRootGroupOpened = XFALSE;
  }
  else {
    /* not a group or file identifier... problem.*/
    return -1;
  }

  /* open all groups under this group*/
  status = H5Gget_num_objs(Id, &nObjs);
  if (status < 0) {
    return status;
  }
  /* loop through all the objects in the group*/
  for (i = 0; i < nObjs; i++) {
    /* see if the object is a group*/
    iObjectType = H5Gget_objtype_by_idx(Id, i);
    if (iObjectType == H5G_GROUP) {
      /* Get the name of the group*/
      iNameLength = H5Gget_objname_by_idx(Id, i, NULL, 0);
      if (iNameLength < 0) {
        char *myerror = NULL;
        myerror = (char *)malloc(81*sizeof(char));
        strcpy(myerror,"Error in xfpRecurseCountGroups. The length of a group's name could not be found.");
        xfpAddXMDFError(myerror);
        free (myerror);
        return iNameLength;
      }
      Path = (char *)malloc(iNameLength+1*sizeof(char));
      iNameLength = H5Gget_objname_by_idx(Id, i, Path, iNameLength+1);
      if (iNameLength < 0) {
        return iNameLength;
      }
      /* Add NULL terminator*/
      Path[iNameLength] = '\0';

      /* Open the group*/
      ChildId = H5Gopen(Id, Path);
      if (Path) {
        free(Path);
        Path = NULL;
      }
      if (ChildId > 0) {
        /* See if the group type is the one we are counting*/
        if (xfpGetGroupType(ChildId, &GroupType) >= 0) {
          if (strcmp(GroupType, GROUP_TYPE_GENERIC) != 0) {
            bGenericGroup = XFALSE;
          }
          else {
            bGenericGroup = XTRUE;
          }
          if (strcmp(GroupType, a_CountObj->GroupType) == 0) {
            status = xfpGetAbsolutePath(ChildId, &Path);
            if (status >= 0) {
              PathLength = strlen(Path) + 1;
              if (a_CountObj->bCounting) {
                a_CountObj->nGroups++;
                a_CountObj->MaxStringLength = Xmax(a_CountObj->MaxStringLength, 
                                                 PathLength);
              }
              else {
                if (PathLength > a_CountObj->MaxStringLength || 
                    a_CountObj->iCurPath >= a_CountObj->nGroups) {
                  return ERROR_OTHER;
                }
                strcpy(&a_CountObj->Paths[a_CountObj->iCurPath * 
                                          a_CountObj->MaxStringLength], Path);
                a_CountObj->iCurPath++;
              }
              free(Path);
            }
          }
          free(GroupType);
          /* always look in generic subfolders*/
          if (a_CountObj->bLookInSubgroups || bGenericGroup) {
            /* Iterate through this group*/
            xfpRecurseCountGroups(ChildId, a_CountObj);
          }
        }
        else {
            /* Iterate through this group*/
          xfpRecurseCountGroups(ChildId, a_CountObj);
        }
        H5Gclose(ChildId);
      }
    }
  }

  if (bRootGroupOpened) {
    H5Gclose(Id);
  }
#if 0
Old method
  H5G_stat_t   ObjStats;
  int          status;
  xid          a_GroupId;
  char        *Path;
  int          PathLength;
  char        *GroupType;
  boolean      bGenericGroup = XTRUE;

  H5Gget_objinfo(a_Id, a_Name, XFALSE, &ObjStats);
  switch (ObjStats.type) {
    case H5G_GROUP:
      a_GroupId = H5Gopen(a_Id, a_Name);
      if (a_GroupId < 0) {
        return a_GroupId;
      }

      if (xfpGetGroupType(a_GroupId, &GroupType) >= 0) {
        if (strcmp(GroupType, GROUP_TYPE_GENERIC) != 0) {
          bGenericGroup = XFALSE;
        }
        if (strcmp(GroupType, a_CountObj->GroupType) == 0) {
          status = xfpGetAbsolutePath(a_GroupId, &Path);
          if (status >= 0) {
            PathLength = strlen(Path) + 1;
            if (a_CountObj->bCounting) {
              a_CountObj->nGroups++;
              a_CountObj->MaxStringLength = Xmax(CountObj->MaxStringLength, 
                                               PathLength);
            }
            else {
              if (PathLength > a_CountObj->MaxStringLength || 
                  a_CountObj->iCurPath >= a_CountObj->nGroups) {
                return ERROR_OTHER;
              }
              strcpy(&a_CountObj->Paths[a_CountObj->iCurPath * 
                                        a_CountObj->MaxStringLength], Path);
              a_CountObj->iCurPath++;
            }
            free(Path);
          }
        }
        free(GroupType);
        /* always look in generic subfolders */
        if (a_CountObj->bLookInSubgroups || bGenericGroup) {
          /* Iterate through this group */
          H5Giterate(a_Id, a_Name, NULL, xfpRecurseCountGroups, a_CountObj);
        }
      }
      else {
        /* Iterate through this group */
        H5Giterate(a_Id, a_Name, NULL, xfpRecurseCountGroups, a_CountObj);
      }
      H5Gclose(a_GroupId);
      return 0;
      break;
    default:
      return 0;
      break;
  }
#endif
  /* continue */
  return 0;
#else
  return -1;
#endif
} /* xfpRecurseCountGroups */

/******************************************************************************
 * FUNCTION  xfpNumGroupsOfType
 * PURPOSE   Determine the number of groups below the group that have a certain
 *           grouptype identifier.
 * NOTES     This is a recursive function that determines all the paths for
 *           a group. 
 ******************2***********************************************************/
XMDF_API int xfpNumGroupsOfType(xid a_Id, const char *a_GroupType,
                              int *a_NumGroups, int *a_MaxLength, 
                              xmbool a_bLookInSubgroups)
{
  int               status = 1;
  XGroupIteration   CountObj;

  CountObj.bCounting = XTRUE;
  CountObj.bLookInSubgroups = a_bLookInSubgroups;
  CountObj.GroupType = a_GroupType;
  CountObj.nGroups = 0;
  CountObj.MaxStringLength = 0;
  CountObj.Paths = NULL;
  CountObj.iCurPath = 0;

  /* RDJ - Need to change to work for groups other than the root directory */
  /*status = H5Giterate(a_Id, NULL, NULL, xfpRecurseCountGroups, &CountObj); */
  xfpRecurseCountGroups(a_Id, &CountObj);

  *a_NumGroups = CountObj.nGroups;
  *a_MaxLength = CountObj.MaxStringLength;

  return status;
} /* xfpNumGroupsOfType */
/******************************************************************************
 * FUNCTION  xfpPathGroupsOfType
 * PURPOSE   Get the paths to all scalar datasets below a starting group
 * NOTES     NumPaths and MaxPathLength are used to verify that the paths are 
 *           not to long
 ******************2***********************************************************/
XMDF_API int xfpPathGroupsOfType(xid a_Id, const char *a_GroupType, 
               int a_NumPaths, int a_MaxPathLength, char *a_Paths, 
               xmbool a_bLookInSubgroups)
{
  int         status = 1;
  XGroupIteration   CountObj;
  char             *StartPath = NULL;
  int              iCurPath;
  char             *PathsCopy = NULL;

  CountObj.bCounting = XFALSE;
  CountObj.bLookInSubgroups = a_bLookInSubgroups;
  CountObj.GroupType = a_GroupType;
  CountObj.nGroups = a_NumPaths;
  CountObj.MaxStringLength = a_MaxPathLength;
  CountObj.Paths = a_Paths;
  CountObj.iCurPath = 0;

  /*status = H5Giterate(a_Id, NULL, NULL, xfpRecurseCountGroups, &CountObj); */
  status = xfpRecurseCountGroups(a_Id, &CountObj);

  /* If not calling from root directory remove the path from each path*/
  status = xfpGetAbsolutePath(a_Id, &StartPath);
  if (status >= 0) {
    PathsCopy = (char *)malloc(a_NumPaths*a_MaxPathLength*sizeof(char));
    memcpy(PathsCopy, CountObj.Paths, a_NumPaths*a_MaxPathLength*sizeof(char));
    for (iCurPath = 0; iCurPath < a_NumPaths; iCurPath++) {
      if (strstr(&PathsCopy[iCurPath*a_MaxPathLength], StartPath) != NULL) {
        /* copy everything past the root directory and the slash afterwards */
        strcpy(&a_Paths[iCurPath*a_MaxPathLength],
               &PathsCopy[iCurPath*a_MaxPathLength+strlen(StartPath) + 1]);
      }
    }
    if (StartPath) {
      free(StartPath);
    }
    free(PathsCopy);
  }
  status = 1;

  return status;
} /* xfpPathGroupsOfType */
/******************************************************************************
 * FUNCTION  xfpGetGroupType
 * PURPOSE   Read the group type string and copy into a string which must 
 * NOTES     be deallocated using free
 *           Returns XFALSE if a group type string was not read
 ******************2***********************************************************/
XMDF_API int xfpGetGroupType (xid a_Id, char **a_GroupType)
{
#ifndef _WIN64
  int   status = 1;
  herr_t (*old_func)(void*);
  void    *old_client_data;

  H5Eget_auto(&old_func, &old_client_data);

  /* Turn off error handling because we may get an error if the attribute
   * doesn't exist but it doesn't mean there was an error */
  H5Eset_auto(NULL, NULL);

  status = xfpReadAttributeString(a_Id, GROUP_TYPE, 1, a_GroupType);

  /* Restore previous error handler */
  H5Eset_auto(old_func, old_client_data);

  return status;
#else
  return -1;
#endif
} /* xfpGetGroupType */
/******************************************************************************
 * FUNCTION  xfpIsGroupOfType
 * PURPOSE   verifies that the group has the group type indicated
 * NOTES     
 ******************2***********************************************************/
 XMDF_API xmbool xfpIsGroupOfType (xid a_Id, const char *a_GroupType)
{
#ifndef _WIN64
  char  *ActualGroupType;
  int   status;
  xmbool same;
  herr_t (*old_func)(void*);
  void    *old_client_data;

  H5Eget_auto(&old_func, &old_client_data);

  /* Turn off error handling because we may get an error if the attribute
   * doesn't exist but it doesn't mean there was an error */
  H5Eset_auto(NULL, NULL);

  status = xfpReadAttributeString(a_Id, GROUP_TYPE, 1, &ActualGroupType);

  /* Restore previous error handler */
  H5Eset_auto(old_func, old_client_data);

  if (status < 0) {
    /* couldn't read group flag must not be of type */
    return XFALSE;
  }

  same = (unsigned char)(strcmp(a_GroupType, ActualGroupType) == 0);

  free(ActualGroupType);
  return same;
#else
  return -1;
#endif
} /* xfpIsGroupOfType */
/*-----------------------------------------------------------------------------*/
/* FUNCTION  xfpGetAbsolutePath*/
/* PURPOSE   calls xfGetGroupAbsolutePathSize and xfGetGroupAbsolutePath to*/
/*           obtain the absolute path for the given Id*/
/* NOTES*/
/*-----------------------------------------------------------------------------*/
XMDF_API int xfpGetAbsolutePath (xid a_Id, char **a_Path)
{
  int         status = -1;
  int         PathLength;

  xfGetGroupAbsolutePathSize(a_Id, &PathLength);
  if (PathLength > 1) {
    *a_Path = (char *)malloc((PathLength)*sizeof(char));
    xfGetGroupAbsolutePath(a_Id, PathLength, *a_Path);
    status = 0;
  }

  return status;
} /* xfpGetAbsolutePath*/

  /* always call this function before using */
/******************************************************************************
 * FUNCTION  xfpDsetParamsInit
 * PURPOSE
 * NOTES     
 ******************2***********************************************************/
 XMDF_API void xfpDsetParamsInit (XDatasetParams *a_Params, int a_Rank,
                               xmbool a_bChunked, int a_Compression)
{
  int             i;

  a_Params->Rank = a_Rank;
  a_Params->Dims = (hsize_t *)malloc(a_Rank*sizeof(hsize_t));
  a_Params->Maxdims = (hsize_t *)malloc(a_Rank*sizeof(hsize_t));

  if (a_bChunked) {
    a_Params->bChunked = XTRUE;
    a_Params->Chunksize = (hsize_t *)malloc(a_Rank*sizeof(hsize_t));
  }
  else {
    a_Params->bChunked = XFALSE;
    a_Params->Chunksize = NULL;
  }

  /* initialize dimensions to NONE */
  for (i = 0; i < a_Rank; i++) {
    a_Params->Dims[i] = NONE;
    a_Params->Maxdims[i] = NONE;
    if (a_bChunked) {
      a_Params->Chunksize[i] = NONE;
    }
  }

  /* default use fill value to false */
  a_Params->bUseFillValue = XFALSE;
  a_Params->fill_type = 0;
  a_Params->fill_value = NULL;

  a_Params->Compression = a_Compression;
  a_Params->bwFlags = 0;
} /* xfpDsetParamsInit */
/******************************************************************************
 * FUNCTION  xfpDsetParamsReset
 * PURPOSE
 * NOTES     
 ******************2***********************************************************/
 XMDF_API void xfpDsetParamsReset (XDatasetParams *a_Params, int a_Rank,
                                xmbool a_bChunked, int a_Compression)
{
  int           i;

  /* see if the size needs to change */
  if (a_Params->Rank != a_Rank) {
    a_Params->Rank = a_Rank;
    if (a_Params->Dims) {
      a_Params->Dims = (hsize_t *)realloc(a_Params->Dims, 
                                          a_Rank*sizeof(hsize_t));
    }
    else {
      a_Params->Dims = (hsize_t *)malloc(a_Rank*sizeof(hsize_t));
    }
    if (a_Params->Maxdims) {
      a_Params->Maxdims = (hsize_t *)realloc(a_Params->Maxdims, 
                                             a_Rank*sizeof(hsize_t));
    }
    else {
      a_Params->Maxdims = (hsize_t *)malloc(a_Rank*sizeof(hsize_t));
    }
    a_Params->bChunked = a_bChunked;
    if (a_bChunked) {
      if (a_Params->Chunksize) {
        a_Params->Chunksize = (hsize_t *)realloc(a_Params->Chunksize, 
                                             a_Rank*sizeof(hsize_t));
      }
      else {
        a_Params->Chunksize = (hsize_t *)malloc(a_Rank*sizeof(hsize_t));
      }
    }
    else {
      free (a_Params->Chunksize);
      a_Params->Chunksize = NULL;
    }
    /* reinitialize items to NONE */
    for (i = 0; i < a_Rank; i++) {
      a_Params->Dims[i] = NONE;
      a_Params->Maxdims[i] = NONE;
      if (a_bChunked) {
        a_Params->Chunksize[i] = NONE;
      }
    }
  }
  else {
    /* see if we changed chunking */
    if (a_Params->bChunked != a_bChunked) {
      a_Params->bChunked = a_bChunked;
      if (a_bChunked) {
        if (a_Params->Chunksize) {
          a_Params->Chunksize = (hsize_t *)realloc(a_Params->Chunksize, 
                                               a_Rank*sizeof(hsize_t));
        }
        else {
          a_Params->Chunksize = (hsize_t *)malloc(a_Rank*sizeof(hsize_t));
        }
        for (i = 0; i < a_Rank; i++) {
          a_Params->Chunksize[i] = NONE;
        }
      }
      else {
        if (a_Params->Chunksize) {
          free(a_Params->Chunksize);
          a_Params->Chunksize = NULL;
        }
      }
    }
  }

  /* default use fill value to false */
  a_Params->bUseFillValue = XFALSE;
  a_Params->fill_type = 0;
  a_Params->fill_value = NULL;

  a_Params->Compression = a_Compression;
} /* xfpDsetParamsReset */
/******************************************************************************
 * FUNCTION  xfpDsetParamsDestroy
 * PURPOSE
 * NOTES     
 ******************2***********************************************************/
 XMDF_API void xfpDsetParamsDestroy (XDatasetParams *a_Params)
{
  if (a_Params->Dims) {
    free(a_Params->Dims);
    a_Params->Dims = NULL;
  }
  if (a_Params->Maxdims) {
    free(a_Params->Maxdims);
    a_Params->Maxdims = NULL;
  }
  if (a_Params->Chunksize) {
    free(a_Params->Chunksize);
    a_Params->Chunksize = NULL;
  }
} /* xfpDsetParamsDestroy */
/******************************************************************************
 * FUNCTION  xfpDsetParamsSetSizes
 * PURPOSE
 * NOTES     
 ******************2***********************************************************/
XMDF_API int xfpDsetParamsSetSizes (XDatasetParams *a_Params, int a_Index,
                           hsize_t a_Dim, hsize_t a_Maxdim, hsize_t a_Chunkdim)
{
  if (a_Index >= a_Params->Rank) {
    return ERROR_DATASET_SIZE_INCORRECT;
  }
  a_Params->Dims[a_Index] = a_Dim;
  a_Params->Maxdims[a_Index] = a_Maxdim;
  if (a_Params->Chunksize) {
    a_Params->Chunksize[a_Index] = a_Chunkdim;
  }
  return 1;
} /* xfpDsetParamsSetSizes */
/******************************************************************************
 * FUNCTION  xfpDsetParamsCreateSpace
 * PURPOSE   Creates a dataspace from the information in a XDatasetParam
 * NOTES     
 ******************2***********************************************************/
XMDF_API int xfpDsetParamsCreateSpace (const XDatasetParams *a_Params)
{
#ifndef _WIN64
  return H5Screate_simple(a_Params->Rank, a_Params->Dims, a_Params->Maxdims);
#else
  return -1;
#endif
} /* xfpDsetParamsCreateSpace */
/******************************************************************************
 * FUNCTION  xfpDsetParamsCreateProps
 * PURPOSE   Creates dataset creation properties
 * NOTES     
 ******************2***********************************************************/
XMDF_API int xfpDsetParamsCreateProps (const XDatasetParams *a_Params)
{
#ifndef _WIN64
  hid_t  PropId;
  herr_t status;

  PropId = H5Pcreate(H5P_DATASET_CREATE);
  if (PropId > 0) {
    if (a_Params->bChunked) {
      status = H5Pset_chunk(PropId, a_Params->Rank, a_Params->Chunksize);
      if (status < 0) {
        H5Pclose(PropId);
        return (status);
      }
    }
    if ((a_Params->Compression != NONE || XFBitTest(a_Params->bwFlags, XF_DSET_SZIP)) &&
         XFBitTest(a_Params->bwFlags, XF_DSET_SHUFFLE)) {
      /* Set the shuffle*/
      status = H5Pset_shuffle(PropId);
    }
    if (a_Params->Compression != NONE) {
      status = H5Pset_deflate(PropId, a_Params->Compression);
    }
    if (a_Params->bUseFillValue) {
      status = H5Pset_fill_value(PropId, a_Params->fill_type, a_Params->fill_value);
    }
    if (XFBitTest(a_Params->bwFlags, XF_DSET_SZIP)) {
      unsigned szip_options_mask=H5_SZIP_ALLOW_K13_OPTION_MASK|H5_SZIP_NN_OPTION_MASK;
      unsigned szip_pixels_per_block=16;
      H5Pset_szip(PropId, szip_options_mask, szip_pixels_per_block);
    }
  }

  return PropId;
#else
  return -1;
#endif
} /* xfpDsetParamsCreateProps */

/******************************************************************************
 * FUNCTION  xfpDsetParamsUseFillValue
 * PURPOSE   Creates dataset creation properties
 * NOTES     
 ******************2***********************************************************/
XMDF_API int xfpDsetParamsUseFillValue (XDatasetParams *a_Params, xid a_type,
                void *a_value)
{
  a_Params->bUseFillValue = XTRUE;
  a_Params->fill_type = a_type;
  a_Params->fill_value = a_value;
  return 1;
} /* xfpDsetParamsUseFillValue */

/******************************************************************************
 * FUNCTION  xfpDsetParamsSetSZip
 * PURPOSE   
 * NOTES     
 ******************2***********************************************************/
XMDF_API int xfpDsetParamsSetSZip (XDatasetParams *a_Params, xmbool a_on)
{
  XFBitSet(a_Params->bwFlags, XF_DSET_SZIP, a_on);
  return 1;
} /* xfpDsetParamsSetSZip */

/******************************************************************************
 * FUNCTION  xfpDsetParamsSetSZip
 * PURPOSE   
 * NOTES     
 ******************2***********************************************************/
XMDF_API int xfpDsetParamsSetShuffle (XDatasetParams *a_Params, xmbool a_on)
{
  XFBitSet(a_Params->bwFlags, XF_DSET_SHUFFLE, a_on);
  return 1;
} /* xfpDsetParamsSetShuffle */
/******************************************************************************
 * FUNCTION  xfpSetParallel
 * PURPOSE   
 * NOTES     
 ******************************************************************************/
 XMDF_API void xfpSetParallel (int a_rank)

{
  fg_running_parallel=XTRUE;
  fg_parallel_rank = a_rank;

} /* xfpSetParallel */
/******************************************************************************
 * FUNCTION  xfpGetRunningParallel
 * PURPOSE   
 * NOTES     
 ******************************************************************************/
 XMDF_API xmbool xfpGetRunningParallel ()
{
  return fg_running_parallel;

} /* xfpGetRunningParallel */
/******************************************************************************
 * FUNCTION  xfpGetParallelRank
 * PURPOSE   
 * NOTES     
 ******************************************************************************/
XMDF_API int xfpGetParallelRank ()
{
  return fg_parallel_rank;

} /* xfpGetParallelRank */
/******************************************************************************
 * FUNCTION  xfpGetParallelNumValuesToRead
 * PURPOSE   This function is to return a number of values to read in for a 
 *           datset at a time so that we dont use all of the RAM
 * NOTES     
 ******************************************************************************/
 XMDF_API hssize_t xfpGetParallelNumValuesToRead ()
{
  return fg_parallel_numdatavalues;

} /* xfpGetParallelNumValuesToRead */
/******************************************************************************
 * FUNCTION  xfpSetParallelNumValuesToRead
 * PURPOSE   This function allows you to set a number of values to read in for a 
 *           datset at a time so that we dont use all of the RAM
 * NOTES     
 ******************************************************************************/
 XMDF_API void xfpSetParallelNumValuesToRead (hssize_t a)
{
  fg_parallel_numdatavalues=a;

} /* xfpSetParallelNumValuesToRead */
/******************************************************************************
 * FUNCTION  xfpStridedToContinous
 * PURPOSE   creates a new array based on the stride that is passed in
 * NOTES     must free the memory that is created from the new array
 *           a_OldArray     = old array of values
 *           a_NumValues    = total number of new vectors/scalars in the array
 *           a_Stride       = the number of bytes to skip between each double
 *           a_NumPerStride = number of vaules to read in between each of the stride jumps
 *           a_StartLoc     = Where to start reading in the strided array 
  ******************************************************************************/
  XMDF_API double*  xfpStridedToContinousDouble (const double *a_OldArray,
                                    const int a_NumValues, const int a_Stride,
                                    const int a_NumPerStride, const int a_StartLoc)
{
  char          *ptru;
  double        *ptrd;
  int            i;
  double        *NewArray;

  /*malloc the new array*/
  NewArray = (double*) malloc(a_NumValues*a_NumPerStride*sizeof(double));

  /*get a new array with the stride*/

  /*set the unsigned char pointer to the beginng of the array*/
  ptru = (char *)(a_OldArray+a_StartLoc);
  ptrd = NewArray; 

  for (i=0; i<a_NumValues; i++) {
    /*this copies a complete vector/scalar from the old array to the new array*/
    memcpy(ptrd, ptru, a_NumPerStride*sizeof(double));

    /*shift the pointers*/ 
    ptru += a_NumPerStride*sizeof(double);
    ptrd += a_NumPerStride;
    /* add the stride to the pointer*/
    /*ptru += a_Stride*sizeof(double); */
     ptru += a_Stride;
  }

  return NewArray;
}
/******************************************************************************
 * FUNCTION  xfpStridedToContinous
 * PURPOSE   creates a new array based on the stride that is passed in
 * NOTES     a_OldArray     = old array of values
 *           a_NewArray     = the new array of values MUST BE MALLOCED ALLREADY
 *           a_NumValues    = total number of new values in the array
 *           a_Stride       = the number of bytes to skip between each double
 *           a_NumPerStride = number of values to read in between each of the stride jump
 *           a_StartLoc     = Where to start writing into the strided array 
  ******************************************************************************/
XMDF_API int xfpContinousToStridedDouble (double *a_OldArray, double** a_NewArray,
                                    const int a_NumValues, const int a_Stride,
                                    const int a_NumPerStride, const int a_StartLoc)
{
  char          *ptru;
  double        *ptrd;
  int            i;

  /*set the unsigned char pointer to the beginng of the array*/
  ptru = (char *) (*a_NewArray + a_StartLoc);
  ptrd = a_OldArray;
  for (i=0; i<a_NumValues; i++) {
    /*this copies a complete vector/scalar from the old array to the new array*/
    memcpy(ptru, ptrd, a_NumPerStride*sizeof(double));

    /*shift the pointers*/ 
    ptru += a_NumPerStride*sizeof(double);
    ptrd += a_NumPerStride;

    /* add the stride to the pointer*/
    ptru += a_Stride*sizeof(double); 
  }

  return 0;
}
/******************************************************************************
 * FUNCTION  xfpStridedToContinous
 * PURPOSE   creates a new array based on the stride that is passed in
 * NOTES     must free the memory that is created from the new array
 *           a_OldArray     = old array of values
 *           a_NumValues    = total number of new values in the array
 *           a_Stride       = the number of bytes to skip between each double
 *           a_NumPerStride = number of vaules to read in between each of the stride jumps
 *           a_StartLoc     = Where to start reading in the strided array 
  ******************************************************************************/
XMDF_API int* xfpStridedToContinousInt (const int *a_OldArray,
                                const int a_NumValues, const int a_Stride,
                                const int a_NumPerStride, const int a_StartLoc)
{
    char          *ptru;
    int           *ptrd;
    int            i;
    int           *NewArray;


  /*malloc the new array*/
  NewArray = (int*) malloc(a_NumValues*a_NumPerStride*sizeof(int));

  /*get a new array with the stride*/

  /*this code replaces the previous code to initialize the pointers*/
	ptru = (char *)(a_OldArray+a_StartLoc);
  ptrd = NewArray; 

  for (i=0; i<a_NumValues-1; i++) {
    /*this copies a complete vector/scalar from the old array to the new array*/
    memcpy(ptrd, ptru, a_NumPerStride*sizeof(int));

    /*shift the pointers*/ 
    ptru += a_NumPerStride*sizeof(int);
    ptrd += a_NumPerStride;

    /* add the stride to the pointer*/
    ptru += a_Stride; 

  }

  return NewArray;
}
/******************************************************************************
 * FUNCTION  xfpStridedToContinous
 * PURPOSE   creates a new array based on the stride that is passed in
 * NOTES     a_OldArray     = old array of values
 *           a_NewArray     = the new array of values MUST BE MALLOCED ALLREADY
 *           a_NumValues    = total number of new values in the array
 *           a_Stride       = the number of bytes to skip between each double
 *           a_NumPerStride = number of values to read in between each of the stride jump
 *           a_StartLoc     = Where to start writing into the strided array 
  ******************************************************************************/
XMDF_API int xfpContinousToStridedInt (int *a_OldArray, int** a_NewArray,
                            const int a_NumValues, const int a_Stride,
                            const int a_NumPerStride, const int a_StartLoc)
{
  char          *ptru;
  int           *ptrd;
  int            i;

  /*set the unsigned char pointer to the beginng of the array*/
  ptru = (char *) (*a_NewArray + a_StartLoc);
  ptrd = a_OldArray;
  for (i=0; i<a_NumValues; i++) {
    /*this copies a complete vector/scalar from the old array to the new array*/
    memcpy(ptru, ptrd, a_NumPerStride*sizeof(int));

    /*shift the pointers*/ 
    ptru += a_NumPerStride*sizeof(int);
    ptrd += a_NumPerStride;

    /* add the stride to the pointer*/
    ptru += a_Stride*sizeof(int); 
  }

  return 0;
}
#ifdef EMRL_XMDF
/******************************************************************************
 * FUNCTION  CallXfiCreateGroup
 * PURPOSE   Calls xfpCreateGroup, for EMRL only
 * NOTES     
 ******************************************************************************/
int ECreateGroup::CallxfpCreateGroup (xid a_Id, const char *Path, 
                                      const char *a_GroupType)
{
  return (xfpCreateGroup(a_Id, Path, a_GroupType));

} /* CCreateGroup::CallxfpCreateGroup */


#endif
