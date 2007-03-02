#ifndef XMDF_TIMESTEP_DOT_H
#define XMDF_TIMESTEP_DOT_H

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

#ifdef __cplusplus
extern "C" {
#endif

/* working with HDF5 datasets */
XMDF_API int xftReadDset1D(xid a_Id, const char *a_Name, int a_Number,
                   H5T_class_t a_DsetClass, xid a_MemType, void *a_Properties);
XMDF_API int xftReadDset1DDouble(xid a_Id, const char *a_Name, int a_Num,
                                double *a_Array);
XMDF_API int xftReadDset1DFloat (xid a_Id, const char *a_Name, int a_Num,
                                float  *a_Array);
XMDF_API int xftReadDsetDouble (xid a_Id, const char *a_Name, int a_Rank,
                                const hsize_t *Dims, double *a_Array);
XMDF_API int xftReadDsetInt (xid a_Id, const char *a_Name, int a_Rank,
                                const hsize_t *Dims, int *a_Array);
  /* portions of HDF5 datasets */
XMDF_API int xftReadDset1DIntPortion(xid a_Id, const char *a_Name,
                hssize_t a_Start1, hsize_t a_Num1, int *a_Array);
XMDF_API int xftWriteDset1DFloatPortion(xid a_Id, const char *a_Name,
                int a_Start1, int a_Num1, float *a_Array);
XMDF_API int xftWriteDset2DFloatPortion(xid a_Id, const char *a_Name,
                int a_Start1, int a_Num1, int a_Start2, int a_Num2,
                float *a_Array);
XMDF_API int xftReadDset2DFloatPortion (xid a_Id, const char *a_Name, 
                hssize_t a_Start1, hsize_t a_Num1,
                hssize_t a_Start2, hsize_t a_Num2, float *a_Array);
XMDF_API int xftReadDset2DFloatPortions (xid a_Id, const char *a_Name, 
                int nPortions, hssize_t *a_Start1, hsize_t *a_Num1,
                hssize_t *a_Start2, hsize_t *a_Num2, float *a_Array);
XMDF_API int xftReadDset2DDoublePortion (xid a_Id, const char *a_Name, 
                hssize_t a_Start1, hsize_t a_Num1,
                hssize_t a_Start2, hsize_t a_Num2, double *a_Array);
XMDF_API int xftReadDset2DIntPortion (xid a_Id, const char *a_Name, 
                hssize_t a_Start1, hsize_t a_Num1,
                hssize_t a_Start2, hsize_t a_Num2, int *a_Array);
XMDF_API int xftReadDset2DUCharPortion (xid a_Id, const char *a_Name, 
                hssize_t a_Start1, hsize_t a_Num1,
                hssize_t a_Start2, hsize_t a_Num2, xmbool *a_Array);
XMDF_API int xftReadDset3DDoublePortion (xid a_Id, const char *a_Name, 
                hssize_t a_Start1, hsize_t a_Num1,
                hssize_t a_Start2, hsize_t a_Num2, hssize_t a_Start3,
                hsize_t a_Num3, double *a_Array);
XMDF_API int xftReadDset3DFloatPortion (xid a_Id, const char *a_Name, 
                hssize_t a_Start1, hsize_t a_Num1,
                hssize_t a_Start2, hsize_t a_Num2, hssize_t a_Start3,
                hsize_t a_Num3, float *a_Array);
XMDF_API int xftReadDset3DDoublePortions (xid a_Id, const char *a_Name,
                int nPortions, hssize_t *a_Start1, hsize_t *a_Num1,
                hssize_t *a_Start2, hsize_t *a_Num2, hssize_t *a_Start3,
                hsize_t *a_Num3, double *a_Array);

  /* Write arrays to datasets */
XMDF_API int xftAppendDset1DDouble(xid a_Id, const char *a_Name,
                                  hsize_t NewTimes, const double *a_Array);
XMDF_API int xftAppendDset1DFloat(xid a_Id, const char *a_Name,
                                  hsize_t a_NumNew, const float *a_Array);
XMDF_API int xftAppendDset1DInt(xid a_Id, const char *a_Name, hsize_t a_NumNew,
                                const int *a_Array);
XMDF_API int xftAppendDset2DFloatFirstDim(xid a_Id, const char *a_Name, 
                                   hsize_t a_NumNew, hsize_t a_NumDim2,
                                   const float *a_Array);
XMDF_API int xftAppendDset2DDoubleFirstDim(xid a_Id, const char *a_Name, 
                                   hsize_t a_NumNew, hsize_t a_NumDim2,
                                   const double *a_Array);
XMDF_API int xftAppendDset2DDoubleSecondDim (xid a_Id, const char *a_Name, 
                                   hsize_t a_NumNew, hsize_t a_NumDim1,
                                   const double *a_Array);
XMDF_API int xftAppendDset3DDoubleFirstDim(xid a_Id, const char *a_Name, 
                                   hsize_t a_NumNew, hsize_t a_NumDim2,
                                   hsize_t a_NumDim3,
                                   const double *a_Array);
XMDF_API int xftAppendDset3DFloatFirstDim(xid a_Id, const char *a_Name, 
                                   hsize_t a_NumNew, hsize_t a_NumDim2,
                                   hsize_t a_NumDim3,
                                   const float *a_Array);
XMDF_API int xftAppendDset2DUCharFirstDim (xid a_Id, const char *a_Name, 
                                   hsize_t a_NumNew, hsize_t a_NumDim2,
                                   const unsigned char *a_Array);

/* working with HDF5 data spaces  */
XMDF_API int xftGetSimpleDataspaceInfoFromName (xid a_GroupId,
						const char *a_Name,
						const xmbool a_isWrite,
						int *a_Rank,
						hsize_t **a_Dims,
						hsize_t **a_MaxDims);
XMDF_API int xftGetSimpleDataspaceInfo (xid a_GroupId,
					const char *a_Name,
					int *a_Rank,
					hsize_t **a_Dims, 
					hsize_t **a_MaxDims);

/* getting and incrementing the finalStep, i.e. NumTimes */
int xftGetNumTimes ( xid a_Id, int* a_NumTimes );
int xftIncNumTimes (xid a_Id);

#ifdef __cplusplus
}
#endif

#endif
