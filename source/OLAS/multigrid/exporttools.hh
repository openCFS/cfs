/* $Id$ */

#ifndef _EXPORTTOOLS_HH_
#define _EXPORTTOOLS_HH_

/***********************************************************
 * This file contains macros for the matrix and vector export in AMG.
 * They are intended as a personal developing tools in OLAS AMG only.
 * The export format, binary or ascii, is steared by the format switch
 * MATRIX_VECTOR_EXPORT_BINARY.
 ***********************************************************/


/***********************************************************
 * BINARY EXPORT
 * The macros for binary matrix export use the external matrix class
 * CRSparseMatrix (by Uwe Fabricius) in the file "CRSparseMatrix_T.hpp",
 * which is NOT contained in the OLAS source. If you do not have this
 * file, you cannot use the binary export feature, until a binary export
 * method for class CRS_Matrix is implemented. Then, of course these
 * macros should use this new export method.
 ***********************************************************/

#ifdef MATRIX_VECTOR_EXPORT_BINARY
#include "CRSparseMatrix_T.hpp"

//////////////////////////////////////////////////
// export a matrix in binary format using my CRS matrix class
//   mat : CRS_Matrix        : the matrix to be exported
//   fn  : const char* const : export file name
#define  EXPORT_MATRIX( mat, fn ) { \
if( !mat.IsParallel() || ((BaseParallel&)mat).GetCommRank() == 1 ) { \
    std::cout \
    << "\033[1mBINARY\033[0m \033[36mmatrix export:\033[0m "<<fn \
    << " ("<<__FILE__<<", "<<__LINE__<<")\n"; \
} \
CRSparseMatrix<T> expMatrix; \
expMatrix.setSize(  mat.GetNrows(), mat.GetNcols() ); \
for( int matr = 1; matr <= mat.GetNrows(); matr++ ) { \
    for( int matij = mat.GetRowPointer()[matr]; \
             matij < mat.GetRowPointer()[matr+1]; matij++ ) { \
        expMatrix.setEntry( matr-1, (mat.GetColPointer())[matij]-1, \
                            (mat.GetDataPointer())[matij] ); \
    } \
} \
expMatrix.saveBinary( fn, "Exported from OLAS AMG" ); \
}

//////////////////////////////////////////////////
// export a vector in binary format
//   vec : Vector            : the vector to be exported
//   fn  : const char* const : export file name
#define  EXPORT_VECTOR( vec, fn ) { \
if( !vec.IsParallel() || ((BaseParallel&)vec).GetCommRank() == 1 ) { \
    std::cout \
    << "\033[1mBINARY\033[0m \033[36mvector export:\033[0m "<<fn \
    << " ("<<__FILE__<<", "<<__LINE__<<")\n"; \
} \
FILE *vecfile = fopen( fn, "wb" ); \
if( !vecfile ) { \
    std::cerr \
    << "\033[41;31;1m ERROR: \033[0m could not open file \033[1m" \
    << fn<<"\033[0m in macro EXPORT_VECTOR ("<<__FILE__<<", line" \
    << __LINE__<<std::endl; \
} else { \
    int vecsize = vec.GetDataSize(); \
    fwrite( &vecsize,           sizeof(int), 1,       vecfile ); \
    fwrite( vec.GetPointer()+1, sizeof(T),   vecsize, vecfile ); \
    fclose( vecfile ); \
} \
}

/***********************************************************
 * ASCII EXPORT
 ***********************************************************/

#else // EXPORT_BINARY -> this is the ASCII part

//////////////////////////////////////////////////
// export a matrix in ASCII format using the OLAS export function
//   mat : CRS_Matrix        : the matrix to be exported
//   fn  : const char* const : export file name
#define  EXPORT_MATRIX( mat, fn ) { \
if( !mat.IsParallel() || ((BaseParallel&)mat).GetCommRank() == 1 ) { \
    std::cout \
    << "\033[1mASCII\033[0m \033[36mmatrix export:\033[0m "<<fn \
    << " ("<<__FILE__<<", "<<__LINE__<<")\n"; \
} \
mat.Export( fn ); \
}

//////////////////////////////////////////////////
// export a vector in ASCII format
//   vec : Vector            : the vector to be exported
//   fn  : const char* const : export file name
#define  EXPORT_VECTOR( vec, fn ) { \
if( !vec.IsParallel() || ((BaseParallel&)vec).GetCommRank() == 1 ) { \
    std::cout \
    << "\033[1mASCII\033[0m \033[36mvector export:\033[0m "<<fn \
    << " ("<<__FILE__<<", "<<__LINE__<<")\n"; \
} \
vec.Export( fn ); \
}

/***********************************************************
 * ASCII IMPORT
 ***********************************************************/

// #define  IMPORT_PRE_MATRIX( PREMATRIX, FILENAME ) { \
// std::cout \
// << "\033[1mASCII\033[0m \033[36mmatrix import:\033[0m "<<FILENAME \
// << "\n      \033[36mfrom file "<<__FILE__<<", line "<<__LINE__ \
// << "\033[0m\n"; \
// CRSparseMatrix<T> impMatrix; \
// if( impMatrix.loadASCII(FILENAME) && impMatrix.compress() ) { \
//     PREMATRIX.Reset(); \
//     PREMATRIX.SetNumRows( impMatrix.getNumRows() ); \
//     for( Integer IMPI = 0; IMPI < impMatrix.getNumRows(); IMPI++ ) { \
//         for( Integer IMPIJ = impMatrix.getCompressedRowPointer()[IMPI]; \
//                  IMPIJ < impMatrix.getCompressedRowPointer()[IMPI+1]; IMPIJ++ ) { \
//             PREMATRIX.SetEntry( IMPI+1, \
//                                 impMatrix.getCompressedColPointer()[IMPIJ]+1, \
//                                 impMatrix.getCompressedDataPointer()[IMPIJ]  ); \
//         } \
//     } \
// } \
// }

#endif // EXPORT_BINARY

/***********************************************************/

//////////////////////////////////////////////////
// exports a matrix, whereby the file name is built using
// sprintf and one additional parameter
//   MATRIX     : matrix
//   NAMEBUFF   : char buffer, that can take the file name string
//   NAMEFORMAT : format string for sprintf
//   ID         : single paramter for sprintf
#define  EXPORT_MATRIX_FN1( MATRIX, NAMEBUFF, NAMEFORMAT, ID ) { \
sprintf( NAMEBUFF, NAMEFORMAT, ID ); \
EXPORT_MATRIX( MATRIX, NAMEBUFF ); }

//////////////////////////////////////////////////
// exports a vector, whereby the file name is built using
// sprintf and one additional parameter
//   VECTOR     : vector
//   NAMEBUFF   : char buffer, that can take the file name string
//   NAMEFORMAT : format string for sprintf
//   ID         : single paramter for sprintf
#define  EXPORT_VECTOR_FN1( VECTOR, NAMEBUFF, NAMEFORMAT, ID ) { \
sprintf( NAMEBUFF, NAMEFORMAT, ID ); \
EXPORT_VECTOR( VECTOR, NAMEBUFF ); }

/***********************************************************/

#endif // _EXPORTTOOLS_HH_
