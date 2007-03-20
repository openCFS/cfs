#ifndef ERRORDEFINITIONS_DOT_H
#define ERRORDEFINITIONS_DOT_H

// This file contains the error codes for errors thrown only from XFORMAT
// If the error came from HDF5, you have to use their error routines

// File errors -40xx
#define ERROR_FILE_NOT_HDF5    -4001 
#define ERROR_FILE_NOT_XMDF    -4002 

// Attribute errors -41xx
#define ERROR_ATTRIBUTE_NOT_SUPPORTED -4101

// Datatype errors -42xx
#define ERROR_INCORRECT_DATATYPE -4201

// Dataset errors -43xx
#define ERROR_DATASET_SIZE_INCORRECT -4301
#define ERROR_DATASET_NO_DATA        -4302
#define ERROR_DATASET_DOES_NOT_EXIST -4303
#define ERROR_DATASET_INVALID        -4304

// Group errors -44xx
#define ERROR_GROUP_TYPE_INCONSISTENT -4401

// Mesh errors -45xx
#define ERROR_ELEMENT_NUM_INCORRECT   -4501 // inconsistent element number
#define ERROR_NODE_NUM_INCORRECT      -4502
#define ERROR_NOT_MESH_GROUP          -4503
#define ERROR_MESH_INCOMPLETE         -4504
#define ERROR_MESH_INVALID            -4505

// Grid errors
#define ERROR_GRID_TYPE_INVALID       -4601
#define ERROR_GRID_NUM_DIMENSIONS_INVALID -4602
#define ERROR_GRID_EXTRUDE_TYPE_INVALID  -4603
#define ERROR_GRID_NUMVALS_INCORRECT    -4604

// Errors that I don't want to name 
#define ERROR_OTHER                   -9901

#endif
