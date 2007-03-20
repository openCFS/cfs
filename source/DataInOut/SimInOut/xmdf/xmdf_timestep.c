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

#include "xmdf_private.h"
#include "xmdf_timestep.h"
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

// /* WARNING: This file is intended to be used only by EMRL emplyees. */

// XMDF_API int xfpGetNumErrorMessages (int *a_Num)
// XMDF_API int xfpGetErrorMessages (int a_Num, char **a_Errors)
// XMDF_API void xfpAllocateErrorStack ()
// XMDF_API void xfpAddXMDFError (const char* a_Error)
// XMDF_API void xfpClearErrors ()
// XMDF_API herr_t xfpHDF5ErrorHandler (void* client_data)
// XMDF_API herr_t xfpHDF5ErrorWalk_cb (int n, H5E_error_t *err_desc, 
// XMDF_API int xfpCloseOpenIdentifiers (xid a_Id)
// XMDF_API int xfpClearGroup (xid a_Group)
// XMDF_API int xfpGetInitialized ()
// XMDF_API void xfpInitialize ()
// XMDF_API void xfpDefaultWriteFloatType (int a_BigEndian)
// XMDF_API hid_t xfpGetDefaultFloatType ()
// XMDF_API hid_t xfpGetDefaultDoubleType ()
// XMDF_API int xfpReadDset1D (xid a_Id, const char *a_Name, int a_Number,
/******************************************************************************
 * FUNCTION  xftReadDset1DDouble
 * PURPOSE   Read a 1D double array from an HDF5 dataset
 * NOTES     returns positive if successful, otherwise negative
 ****************2***********************************************************/
XMDF_API int xftReadDset1DDouble (xid a_Id, const char *a_Name, int a_Number, 
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


  status = xftGetSimpleDataspaceInfoFromName(a_Id, a_Name, 0, &Rank, &Dims, &Maxdims);

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

  {
    /* select the hyperslab */
    hsize_t Count[1];
    hsize_t Offset[1];

    Count[0] = a_Number;
    Offset[0] = 0;
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
  }

  status = H5Dread(DatasetId, H5T_NATIVE_DOUBLE, H5S_ALL, SpaceId, 
                   H5P_DEFAULT, a_Array);

  H5Sclose(SpaceId);
  H5Dclose(DatasetId);
 
  return status;
#else
  return -1;
#endif
} /* xfpReadDset1DDouble */
/******************************************************************************
 * FUNCTION  xftReadDset1DFloat
 * PURPOSE   Read a 1D float array from an HDF5 dataset
 * NOTES     returns positive if successful, otherwise negative
 ****************2***********************************************************/
XMDF_API int xftReadDset1DFloat (xid a_Id, const char *a_Name, int a_Number, float *a_Array)
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
  status = xftGetSimpleDataspaceInfoFromName(a_Id, a_Name, 0, &Rank, &Dims, &Maxdims);

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

  {
    /* select the hyperslab */
    hsize_t Count[1];
    hsize_t Offset[1];

    Count[0] = a_Number;
    Offset[0] = 0;
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
  }

  status = H5Dread(DatasetId, H5T_NATIVE_FLOAT, H5S_ALL, SpaceId, 
                   H5P_DEFAULT, a_Array);

  H5Sclose(SpaceId);
  H5Dclose(DatasetId);
 
  return status;
#else
  return -1;
#endif
} /* xfpReadDset1DFloat */
// XMDF_API int xfpReadDsetDouble (xid a_Id, const char *a_Name, 
// XMDF_API int xfpReadDsetInt (xid a_Id, const char *a_Name, 
// XMDF_API int xfpReadDset1DIntPortion (xid a_Id, const char *a_Name, 
/******************************************************************************
 * FUNCTION  xftWriteDset2DFloatPortion
 * PURPOSE   Write a portion of a 2D float array to an HDF5 dataset
 * NOTES     returns positive if successful, otherwise negative
 *           a_Start1 is the starting index for the first array index, 
 *           a_Num1 is number to write in first array index
 *           a_Start2 and a_Num2 are same as a_Start1 and a_Num1 but for index2
 *           Data is always stored in a 1D array
 *           Data must already be allocated of size a_Num1 * a_Num2
 ****************2***********************************************************/
XMDF_API int xftWriteDset2DFloatPortion (xid a_Id, const char *a_Name, int a_Start1,
                                int a_Num1, int a_Start2, int a_Num2,
                                float *a_Array)
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

  status = xftGetSimpleDataspaceInfoFromName(a_Id, a_Name, !0, &Rank, &OldDims, &OldMaxDims);
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
 * FUNCTION  xftReadDset2DFloatPortion
 * PURPOSE   Read a portion of a 2D float array from an HDF5 dataset
 * NOTES     returns positive if successful, otherwise negative
 *           a_Start1 is the starting index for the first array index, 
 *           a_Num1 is number to read in first array index
 *           a_Start2 and a_Num2 are same as a_Start1 and a_Num1 but for index2
 *           Data is always stored in a 1D array
 *           Data must already be allocated of size a_Num1 * a_Num2
 ****************2***********************************************************/
XMDF_API int xftReadDset2DFloatPortion (xid a_Id, const char *a_Name, 
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
  status = xftGetSimpleDataspaceInfoFromName(a_Id, a_Name, 0, &Rank, &Dims, &Maxdims);
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
 * FUNCTION  xftReadDset2DFloatPortions
 * PURPOSE   Read a portion of a 2D float array from an HDF5 dataset
 * NOTES     returns positive if successful, otherwise negative
 *           a_Start1 is the starting index for the first array index, 
 *           a_Num1 is number to read in first array index
 *           a_Start2 and a_Num2 are same as a_Start1 and a_Num1 but for index2
 *           Data is always stored in a 1D array
 *           Data must already be allocated of size a_Num1 * a_Num2
 ****************2***********************************************************/
XMDF_API int xftReadDset2DFloatPortions (xid a_Id, const char *a_Name, 
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
  status = xftGetSimpleDataspaceInfoFromName(a_Id, a_Name, 0, &Rank, &Dims, &Maxdims);
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
 * FUNCTION  xftReadDset2DIntPortion
 * PURPOSE   Read a portion of a 2D int array from an HDF5 dataset
 * NOTES     returns positive if successful, otherwise negative
 *           a_Start1 is the starting index for the first array index, 
 *           a_Num1 is number to read in first array index
 *           a_Start2 and a_Num2 are same as a_Start1 and a_Num1 but for index2
 *           Data is always stored in a 1D array
 *           Data must already be allocated of size a_Num1 * a_Num2
 ****************2***********************************************************/
XMDF_API int xftReadDset2DIntPortion (xid a_Id, const char *a_Name, 
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
  status = xftGetSimpleDataspaceInfoFromName(a_Id, a_Name, 0, &Rank, &Dims, &Maxdims);
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
 * FUNCTION  xftReadDset2DDoublePortion
 * PURPOSE   Read a portion of a 2D double array from an HDF5 dataset
 * NOTES     returns positive if successful, otherwise negative
 *           a_Start1 is the starting index for the first array index, 
 *           a_Num1 is number to read in first array index
 *           a_Start2 and a_Num2 are same as a_Start1 and a_Num1 but for index2
 *           Data is always stored in a 1D array
 *           Data must already be allocated of size a_Num1 * a_Num2
 ****************2***********************************************************/
XMDF_API int xftReadDset2DDoublePortion (xid a_Id, const char *a_Name, 
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
  status = xftGetSimpleDataspaceInfoFromName(a_Id, a_Name, 0, &Rank, &Dims, &Maxdims);

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
 * FUNCTION  xftReadDset2DUCharPortion
 * PURPOSE   Read a portion of a 2D float array from an HDF5 dataset
 * NOTES     returns positive if successful, otherwise negative
 *           a_Start1 is the starting index for the first array index, 
 *           a_Num1 is number to read in first array index
 *           a_Start2 and a_Num2 are same as a_Start1 and a_Num1 but for index2
 *           Data is always stored in a 1D array
 *           Data must already be allocated of size a_Num1 * a_Num2
 ****************2***********************************************************/
XMDF_API int xftReadDset2DUCharPortion (xid a_Id, const char *a_Name, 
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
  status = xftGetSimpleDataspaceInfoFromName(a_Id, a_Name, 0, &Rank, &Dims, &Maxdims);
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
 * FUNCTION  xftReadDset3DDoublePortion
 * PURPOSE   Read a portion of a 3D double array from an HDF5 dataset
 * NOTES     returns positive if successful, otherwise negative
 *           a_Start1 is the starting index for the first array index, 
 *           a_Num1 is number to read in first array index
 *           a_Start2 and a_Num2 are same as a_Start1 and a_Num1 but for index2
 *           Data is always stored in a 1D array
 *           Data must already be allocated of size a_Num1 * a_Num2 * a_Num3
 ****************2***********************************************************/
XMDF_API int xftReadDset3DDoublePortion (xid a_Id, const char *a_Name, 
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
  status = xftGetSimpleDataspaceInfoFromName(a_Id, a_Name, 0, &Rank, &Dims, &Maxdims);
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
// XMDF_API int xfpReadDset3DDoublePortions (xid a_Id, const char *a_Name,
/******************************************************************************
 * FUNCTION  xftReadDset3DFloatPortion
 * PURPOSE   Read a portion of a 3D float array from an HDF5 dataset
 * NOTES     returns positive if successful, otherwise negative
 *           a_Start1 is the starting index for the first array index, 
 *           a_Num1 is number to read in first array index
 *           a_Start2 and a_Num2 are same as a_Start1 and a_Num1 but for index2
 *           Data is always stored in a 1D array
 *           Data must already be allocated of size a_Num1 * a_Num2 * a_Num3
 ****************2***********************************************************/
XMDF_API int xftReadDset3DFloatPortion (xid a_Id, const char *a_Name, 
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
  status = xftGetSimpleDataspaceInfoFromName(a_Id, a_Name, 0, &Rank, &Dims, &Maxdims);
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
// XMDF_API int xfpWriteDsetInt (xid a_Id, const char *a_Name, hid_t a_Datatype,
// XMDF_API int xfpWriteDsetUInt (xid a_Id, const char *a_Name, hid_t a_Datatype,
// XMDF_API int xfpWriteDsetFloat (xid a_Id, const char *a_Name,
// XMDF_API int xfpWriteDsetDouble (xid a_Id, const char *a_Name,
// XMDF_API int xfpWriteDsetUChar(xid a_Id, const char *a_Name, hid_t a_Datatype,
/******************************************************************************
 * FUNCTION  xftAppendDset1DDouble
 * PURPOSE   
 * NOTES     
 ******************2***********************************************************/
XMDF_API int xftAppendDset1DDouble (xid a_Id, const char *a_Name, 
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

  status = xftGetSimpleDataspaceInfoFromName(a_Id, a_Name, !0, &Rank, &OldDims, &OldMaxDims);
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
 * FUNCTION  xftAppendDset1DFloat
 * PURPOSE   
 * NOTES     
 ******************2***********************************************************/
XMDF_API int xftAppendDset1DFloat (xid a_Id, const char *a_Name, hsize_t a_NumNew,
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

  status = xftGetSimpleDataspaceInfoFromName(a_Id, a_Name, !0, &Rank, &OldDims, &OldMaxDims);
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
// XMDF_API int xfpAppendDset1DInt (xid a_Id, const char *a_Name, hsize_t a_NumNew,
/******************************************************************************
 * FUNCTION  xftAppendDset2DFloatFirstDim
 * PURPOSE   
 * NOTES     
 ******************2***********************************************************/
XMDF_API int xftAppendDset2DFloatFirstDim (xid a_Id, const char *a_Name, 
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

  status = xftGetSimpleDataspaceInfoFromName(a_Id, a_Name, !0, &Rank, &OldDims, &OldMaxDims);
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
 * FUNCTION  xftAppendDset2DDoubleFirstDim
 * PURPOSE   
 * NOTES     
 ******************2***********************************************************/
XMDF_API int xftAppendDset2DDoubleFirstDim (xid a_Id, const char *a_Name, 
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

  status = xftGetSimpleDataspaceInfoFromName(a_Id, a_Name, !0, &Rank, &OldDims, &OldMaxDims);
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
// XMDF_API int xfpAppendDset2DDoubleSecondDim (xid a_Id, const char *a_Name, 
/******************************************************************************
 * FUNCTION  xftAppendDset3DDoubleFirstDim
 * PURPOSE   
 * NOTES     
 ******************2***********************************************************/
XMDF_API int xftAppendDset3DDoubleFirstDim (xid a_Id, const char *a_Name, 
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

  status = xftGetSimpleDataspaceInfoFromName(a_Id, a_Name, !0, &Rank, &OldDims, &OldMaxDims);
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
 * FUNCTION  xftAppendDset3DFloatFirstDim
 * PURPOSE   
 * NOTES     
 ******************2***********************************************************/
XMDF_API int xftAppendDset3DFloatFirstDim (xid a_Id, const char *a_Name, 
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

  status = xftGetSimpleDataspaceInfoFromName(a_Id, a_Name, !0, &Rank, &OldDims, &OldMaxDims);
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
 * FUNCTION  xftAppendDset2DUCharFirstDim
 * PURPOSE   
 * NOTES     
 ******************2***********************************************************/
XMDF_API int xftAppendDset2DUCharFirstDim (xid a_Id, const char *a_Name, 
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

  status = xftGetSimpleDataspaceInfoFromName(a_Id, a_Name, !0, &Rank, &OldDims, &OldMaxDims);
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
// XMDF_API int xfpReadDatasetString (xid a_Loc, const char *a_Name, char **a_Str)
// XMDF_API int xfpWriteDatasetString (xid a_Loc, const char *a_Name, 
// XMDF_API int xfpWriteDatasetStrings (xid a_Loc, const char *a_Name, hsize_t a_Number,
// XMDF_API int xfpDoesAttributeExist (xid a_Loc, const char *a_Name,
// XMDF_API int xfpReadAttributeString (xid a_Loc, const char *a_Name, int a_Num, char **a_Str)
// XMDF_API int xfpReadAttributeInt (xid a_Loc, const char *a_Name, int a_Num, int *a_val)
// XMDF_API int xfpReadAttributeFloat (xid a_Loc, const char * a_Name, int a_Num, float *a_val)
// XMDF_API int xfpReadAttributeDouble (xid a_Loc, const char *a_Name, int a_Num, double *a_val)
// XMDF_API int xfpWriteAttributeString (xid a_Loc, const char * a_Name, 
// XMDF_API int xfpWriteAttributeInt (xid a_Loc, const char *a_Name, int a_Number,
// XMDF_API int xfpWriteAttributeFloat (xid a_Loc, const char *a_Name, 
// XMDF_API int xfpWriteAttributeDouble (xid a_Loc, const char *a_Name, 
/******************************************************************************
 * FUNCTION  xftGetSimpleDataspaceInfoFromName
 * PURPOSE   Get the rank, size, and max size for a simple dataspace
 * NOTES     
 ******************2***********************************************************/
XMDF_API int xftGetSimpleDataspaceInfoFromName (xid a_GroupId,
						const char *a_Name,
						const xmbool a_IsWrite,
						int *a_Rank,
						hsize_t **a_Dims,
						hsize_t **a_MaxDims)
{
#ifndef _WIN64
  int status;
  int numTimes;
  hid_t    DatasetId, SpaceId;

  DatasetId = H5Dopen(a_GroupId, a_Name);
  if (DatasetId < 0) {
    return DatasetId;
  }
  SpaceId = H5Dget_space(DatasetId);
  if (SpaceId < 0) {
    return SpaceId;
  }
  
  status = xfpGetSimpleDataspaceInfo(SpaceId, a_Rank, a_Dims, a_MaxDims);

  H5Sclose(SpaceId);
  H5Dclose(DatasetId);

  // tweak existing timestep count (OldDims[0]) if proscribed by numTimes

  status = xfpReadAttributeInt(a_GroupId, DATASET_NUMTIMES, 1, &numTimes);

  // NumTimes must be incremented before any part of a timestep is written.
  // This ensures that the Active dataset and the other datasets align.
  // However, this requires NumTimes to be decremented before each write.
  if (status >= 0 && a_IsWrite)
    {
      numTimes--;
    }

  if (status >= 0 &&
      numTimes > 0 &&
      (*a_Dims)[0] >= numTimes)
    {
      (*a_Dims)[0] = numTimes; // truncate to numTimes
    }

  return 0;
#else
  return -1;
#endif
} /* xftGetSimpleDataspaceInfoFromName */
// XMDF_API int xfpGetSimpleDataspaceInfo (xid a_GroupId,
// XMDF_API void xfpDestroySimpleDataspaceInfo (hsize_t **a_Dims,
// XMDF_API int xfpCreateGroup (xid a_Id, const char *a_Path, const char *a_GroupType)
// XMDF_API herr_t xfpRecurseCountGroups(xid a_Id, XGroupIteration *a_CountObj)
// XMDF_API int xfpNumGroupsOfType(xid a_Id, const char *a_GroupType,
// XMDF_API int xfpPathGroupsOfType(xid a_Id, const char *a_GroupType, 
// XMDF_API int xfpGetGroupType (xid a_Id, char **a_GroupType)
// XMDF_API int xfpGetAbsolutePath (xid a_Id, char **a_Path)
// XMDF_API int xfpDsetParamsSetSizes (XDatasetParams *a_Params, int a_Index,
// XMDF_API int xfpDsetParamsCreateSpace (const XDatasetParams *a_Params)
// XMDF_API int xfpDsetParamsCreateProps (const XDatasetParams *a_Params)
// XMDF_API int xfpDsetParamsUseFillValue (XDatasetParams *a_Params, xid a_type,
// XMDF_API int xfpDsetParamsSetSZip (XDatasetParams *a_Params, xmbool a_on)
// XMDF_API int xfpDsetParamsSetShuffle (XDatasetParams *a_Params, xmbool a_on)
// XMDF_API int xfpGetParallelRank ()
// XMDF_API int xfpContinousToStridedDouble (double *a_OldArray, double** a_NewArray,
// XMDF_API int* xfpStridedToContinousInt (const int *a_OldArray,
// XMDF_API int xfpContinousToStridedInt (int *a_OldArray, int** a_NewArray,
/******************************************************************************
 * FUNCTION  xftGetNumTimes
 * PURPOSE   This function returns the NumTimes value in the timestep group.
 *	     If no value is present it will return -1.
 * NOTES     An error code is returned on error.
 ******************************************************************************/
int xftGetNumTimes ( xid a_Id, int* a_NumTimes )
{
  int status = xfpReadAttributeInt(a_Id, DATASET_NUMTIMES, 1, a_NumTimes);

  if (status < 0)
    {
      *a_NumTimes = -1;
    }
    
  return XTRUE;
} /* xftGetNumTimes */
/******************************************************************************
 * FUNCTION  xftIncNumTimes
 * PURPOSE   This function increments the NumTimes value in the timestep group.
 * NOTES     A negative error code is returned on error.
 ******************************************************************************/
int xftIncNumTimes (xid a_Id)
{
  int numTimes;
  int status = xftGetNumTimes( a_Id, &numTimes );

  if (status < 0)
    return status;

  if (numTimes < 0)
    return XTRUE;
    
  numTimes++;

  status = xfpWriteAttributeInt(a_Id, DATASET_NUMTIMES, 1, &numTimes);

  return XTRUE;
} /* xftIncNumTimes */
