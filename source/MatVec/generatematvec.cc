// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <complex>

#include <def_use_lapack.hh>


// Include LAPACK stuff
#ifdef USE_LAPACK
#include "OLAS/external/lapack/lapack.hh"
#include "OLAS/external/lapack/lapackgbmatrix.hh"
#endif


// include source code for templated vectors
//#include "MatVec/vector.cc"
//#include "matvec/vector_specialised.cc"
#include "vector.hh"

// include source code for templated matrices
#include "crs_matrix.hh"
#include "scrs_matrix.hh"
#include "diag_matrix.hh"
#include "sbmmatrix.hh"

#include "sbmvector.hh"

// Include that prematrix stuff
// #include "multigrid/prematrix.cc"



namespace CoupledField {



  // >>>>>>>>>>>>>>>>>>>>>>>>>>>> VECTOR PART <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<



  // *********************************************
  //   Macro for generation of templated vectors
  // *********************************************
#define CREATE_VECTOR( MATRIX_ENTRY, VECTORCLASS ) \
if ( ( eType == MATRIX_ENTRY )  ) {\
   retVector = new VECTORCLASS( length );\
   ASSERTMEM( retVector, sizeof( VECTORCLASS ) );\
}
//   (*cla) << " GenerateSingleVectorObject: Generated "\
//          << MACRO2STRING(VECTORCLASS) << " of dimension " << length\
//   << std::endl;\

  // *******************************************************
  //   Dynamic generation of vector associated with matrix
  // *******************************************************
  BaseVector* GenerateVectorObject( const BaseMatrix &m ) {


    BaseVector *retVector = NULL;
    BaseMatrix::StructureType m_structuretype;
    BaseMatrix::EntryType m_entrytype;
    BaseMatrix::StorageType m_storagetype;
    Integer m_Ncols;
    const StdMatrix  *m_StdMatrix;
    const SBM_Matrix  *sbmmat;
    SBM_Vector * sbmvec;

    m_structuretype = m.GetStructureType();
    m_Ncols = m.GetNumCols();

    switch (m_structuretype) {

    case BaseMatrix::SBM_MATRIX:

      // Downcast matrix to SBM_Matrix
      sbmmat = dynamic_cast<const SBM_Matrix*>(&m);

      // Create an SBM_Vector of fitting length
      sbmvec  = new SBM_Vector( m_Ncols );

      // Create and insert subvectors
      for ( Integer i = 0; i < m_Ncols; i++ ) {
        if ( sbmmat->GetPointer(i,i) != NULL ) {
          sbmvec->SetSubVector( dynamic_cast<SingleVector*>
                                (GenerateVectorObject((*sbmmat)(i,i))), i );
        }
      }

      // Set vector entries to zero
      sbmvec->Init();

      // Perform an upcast
      retVector = dynamic_cast<BaseVector*>(sbmvec);
      (*cla) << " GenerateVectorObject: Generated an SBM_Vector" << std::endl;

      break;

    case BaseMatrix::SPARSE_MATRIX:

      m_StdMatrix = dynamic_cast<const StdMatrix*>(&m);
      m_entrytype = m_StdMatrix->GetEntryType();
      m_storagetype = m_StdMatrix->GetStorageType();

      // Delegate work to factory for sequential SingleVectors
      retVector = GenerateSingleVectorObject( m_storagetype, m_entrytype,
                                              (UInt) m_Ncols );
      break;

    default:
      EXCEPTION( "GenerateVectorObject: Unkown matrix type" );
      break;
    }

    // That's it
    return retVector;
  }


  // **************************************************
  //   Dynamic generation of vector of specified type
  // **************************************************
  BaseVector* GenerateSingleVectorObject( const BaseMatrix::StorageType sType,
                                          const BaseMatrix::EntryType eType,
                                          const UInt length ) {


    SingleVector *retVector = NULL;

    // real valued vectors
    CREATE_VECTOR( BaseMatrix::DOUBLE,  Vector<Double>      );

    // complex valued vectors
    CREATE_VECTOR( BaseMatrix::COMPLEX, Vector<Complex>     );

    // scalar vectors to go together with LAPACK matrices
#ifdef USE_LAPACK
    CREATE_VECTOR( BaseMatrix::F77REAL8    , Vector<Double>  );
    CREATE_VECTOR( BaseMatrix::F77COMPLEX16, Vector<Complex> );
#endif

    // Check, if we were able to generate a vector object. If not, complain.
    if ( retVector == NULL ) {
      EXCEPTION( "GenerateVectorObject: Failed to generate vector object" );
    }

    // Initialise vector to zero
    retVector->Init();

    // That's it
    return retVector;
  }

  // *********************************************
  //   Macro for generation of templated vectors
  // *********************************************
#define COPY_VECTOR( MATRIX_ENTRY, VECTORCLASS ) \
if ( ( eType == MATRIX_ENTRY ) ) { \
   const VECTORCLASS& auxVec = dynamic_cast<const VECTORCLASS&>(origVec); \
   retVector = new VECTORCLASS( auxVec ); \
   ASSERTMEM( retVector, sizeof( VECTORCLASS ) ); \
}
//  (*cla) << " GenerateSingleVectorObject: Generated copy of" \
//         << MACRO2STRING(VECTORCLASS) << std::endl; 
  SingleVector* CopySingleVectorObject( const SingleVector& origVec ) {

    SingleVector *retVector = NULL;

    // Determine vector information
    BaseMatrix::EntryType eType = origVec.GetEntryType();

    // real valued vectors
    COPY_VECTOR( BaseMatrix::DOUBLE,  Vector<Double>      );

    // complex valued vectors
    COPY_VECTOR( BaseMatrix::COMPLEX,  Vector<Complex>     );

    // scalar vectors to go together with LAPACK matrices
#ifdef USE_LAPACK
    COPY_VECTOR( BaseMatrix::F77REAL8    ,  Vector<Double>  );
    COPY_VECTOR( BaseMatrix::F77COMPLEX16,  Vector<Complex> );
#endif

    // Check, if we were able to generate a vector object. If not, complain.
    if ( retVector == NULL ) {
      EXCEPTION( "GenerateVectorObject: Failed to generate vector object" );
    }

    // Return vector
    return retVector;
  }

  // >>>>>>>>>>>>>>>>>>>>>>>>>>>> MATRIX PART <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<



  // ***************************************************
  //   Central macro for generation of a matrix object
  //   Its use also takes care of the instantiation in
  //   the case of the templated matrices
  // ***************************************************
#define MATRIX_OBJ(matEntry,matStore,matrix_obj_type)\
if ((etype==matEntry) && (stype==matStore) ){\
retMat = new matrix_obj_type(nrows,ncols,fill);\
ASSERTMEM( retMat, sizeof(matrix_obj_type) );\
(*cla) << " Generated matrix of type: " << MACRO2STRING(matrix_obj_type) \
       << std::endl;}


  // ****************
  //   CRS Matrices
  // ****************
  typedef CRS_Matrix<Double>      CRSRealDof1;
  typedef CRS_Matrix<Complex>     CRSComplexDof1;

  // *****************
  //   SCRS Matrices
  // *****************
  typedef SCRS_Matrix<Double>      SCRSRealDof1;
  typedef SCRS_Matrix<Complex>     SCRSComplexDof1;

  // *****************
  //   Diag Matrices
  // *****************
  typedef Diag_Matrix<Double>      DiagRealDof1;
  typedef Diag_Matrix<Complex>     DiagComplexDof1;

  // *****************************************
  //   Dynamic generation of a matrix object
  // *****************************************
  StdMatrix* GenerateStdMatrixObject( const BaseMatrix::EntryType etype,
                                      const BaseMatrix::StorageType stype,
                                      const Integer nrows,
                                      const Integer ncols,
                                      const Integer fill ) {

    StdMatrix *retMat = NULL;

    // Branch depending on desired storage layout/StdMatrix type
    switch( stype ) {

      // CRS Matrices
    case BaseMatrix::SPARSE_NONSYM:
      MATRIX_OBJ( BaseMatrix::DOUBLE,  BaseMatrix::SPARSE_NONSYM, CRSRealDof1    );
      MATRIX_OBJ( BaseMatrix::COMPLEX, BaseMatrix::SPARSE_NONSYM, CRSComplexDof1 );
      break;

      // SCRS Matrices
    case BaseMatrix::SPARSE_SYM:
      MATRIX_OBJ( BaseMatrix::DOUBLE,  BaseMatrix::SPARSE_SYM, SCRSRealDof1    );
      MATRIX_OBJ( BaseMatrix::COMPLEX, BaseMatrix::SPARSE_SYM, SCRSComplexDof1 );
      break;

      // NSK Matrices
    case BaseMatrix::SKYLINE_NONSYM:
      EXCEPTION( "GenerateStdMatrixObject: SKYLINE_NOSYM not supported, yet!" );
      break;

      // SSK Matrices
    case BaseMatrix::SKYLINE_SYM:
      EXCEPTION( "GenerateStdMatrixObject: SSK Matrix not yet supported" );
      // MATRIX_OBJ( DOUBLE,  SKYLINE_SYM, 1, SSKRealDof1    );
      // MATRIX_OBJ( COMPLEX, SKYLINE_SYM, 1, SSKComplexDof1 );
      break;

      // Diag Matrices
    case BaseMatrix::DIAG:
      MATRIX_OBJ( BaseMatrix::DOUBLE,  BaseMatrix::DIAG, DiagRealDof1    );
      MATRIX_OBJ( BaseMatrix::COMPLEX, BaseMatrix::DIAG, DiagComplexDof1 );
      break;


#ifdef USE_LAPACK
    case BaseMatrix::LAPACK_GBMATRIX:

      // real entries
      if ( etype == BaseMatrix::DOUBLE  ) {
        retMat = new LapackGBMatrix<F77real8,double>( nrows, ncols,
                                                      BaseMatrix::F77REAL8 );
        ASSERTMEM( retMat, sizeof(LapackGBMatrix<F77real8,double>) );
        (*cla) << " Generated matrix of type: LapackGBMatrix<F77Double>"
               << std::endl;
      }

      // complex entries
      else if ( etype == BaseMatrix::COMPLEX  ) {
        retMat = new LapackGBMatrix<F77complex16,Complex>( nrows, ncols,
        		BaseMatrix::F77COMPLEX16 );
        ASSERTMEM(retMat, sizeof(LapackGBMatrix<F77complex16,Complex>) );
        (*cla) << " Generated matrix of type: LapackGBMatrix<F77Complex16>"
               << std::endl;
      }

      // nothing fits
      else {
        EXCEPTION( "Unable to generate LapackGBMatrix" );
      }
      break;

#else

    case LAPACK_GBMATRIX:
      EXCEPTION( "Compile with USE_LAPACK to enable support for LapackGBMatrix" );
      break;

#endif

      // Should not be reached (hopefully :)
    default:
      EXCEPTION( "GenerateStdMatrixObject: Request for unknown MatrixStorageType" );
    }

    // Return pointer to generated matrix object
    return retMat;

  }


  // ***************************************************
  //   Central macro for generation of a matrix object
  //   by use of a copy constructor
  // ***************************************************
#define COPY_MATRIX_OBJ( matEntry, matStore, matrixObjType )\
try {\
   const matrixObjType& auxMat = dynamic_cast<const matrixObjType&>(origMat); \
   retMat = new matrixObjType( auxMat );\
   ASSERTMEM( retMat, sizeof( matrixObjType) );\
   (*cla) << " Generated copy of matrix of type: "\
          << MACRO2STRING( matrixObjType ) \
          << std::endl;\
}\
catch(...){};



  // ******************************************************
  //   Dynamic generation of a matrix object (by copying)
  // ******************************************************
  StdMatrix* CopyStdMatrixObject( const StdMatrix &origMat  ) {


    StdMatrix *retMat = NULL;

    // Determine matrix information
    BaseMatrix::EntryType eType = origMat.GetEntryType();
    BaseMatrix::StorageType sType = origMat.GetStorageType();

    // CRS_Matrix case
    COPY_MATRIX_OBJ( DOUBLE,  SPARSE_NONSYM, CRSRealDof1    );
    COPY_MATRIX_OBJ( COMPLEX, SPARSE_NONSYM, CRSComplexDof1 )

    // SCRS_Matrix case
    COPY_MATRIX_OBJ( DOUBLE,  SPARSE_SYM, SCRSRealDof1    );
    COPY_MATRIX_OBJ( COMPLEX, SPARSE_SYM, SCRSComplexDof1 );

    // DIAG_Matrix case
    COPY_MATRIX_OBJ( DOUBLE,  DIAG, SCRSRealDof1    );
    COPY_MATRIX_OBJ( COMPLEX, DIAG, SCRSComplexDof1 );

    // Invalid or unimplemented case
    if ( retMat == NULL ) {
      EXCEPTION( "GenerateStdMatrixObject: Failed to generate StdMatrix "
               << "by copying an StdMatrix of type ( sType = "
               << BaseMatrix::storageType.ToString( sType ) << " , eType = "
               << BaseMatrix::entryType.ToString( eType ) );
    }


    return retMat;
  }


  // *********************************************
  //   Dynamic generation of a SBM matrix object
  // *********************************************
  SBM_Matrix* GenerateSBMMatrixObject( BaseMatrix::EntryType eType,
                                       UInt rowDim, UInt colDim,
                                       GraphManager *graphManager,
                                       bool symmetric ) {


	BaseMatrix::StorageType sType = BaseMatrix::NOSTORAGETYPE;
    BaseGraph *graph = NULL;
    UInt scalarRows = 0;
    UInt scalarCols = 0;

    // STEP 1: Generate empty SBM_Matrix
    SBM_Matrix *retMat = NULL;
    retMat = new SBM_Matrix( rowDim, colDim, symmetric );
    if ( retMat == NULL ) {
      EXCEPTION( "GenerateSBM_MatrixObject: "
               << "This is the end my friend!\n"
               << "Generation of empty SBM_Matrix failed!" );
    }

    // STEP 2: Populate with sub-matrices
    for ( UInt i = 0; i < rowDim; i++ ) {
      for ( UInt j = 0; j < colDim; j++ ) {

        // Check for need of sub-matrix
        if ( graphManager->SubGraphExists(i,j) ) {

          // Determine storage type
          if ( symmetric == true && i == j ) {
            sType = BaseMatrix::SPARSE_SYM;
          }
          else {
            sType = BaseMatrix::SPARSE_NONSYM;
          }

          // Get hold of graph of sub-matrix
          graph = graphManager->GetGraph(i,j);

          // Determine dimension of sub-matrix
          scalarRows = graphManager->GetGraph(i,i)->GetSize();
          scalarCols = graphManager->GetGraph(j,j)->GetSize();

          // Generate empty sub-matrix
          retMat->SetSubMatrix( i, j, eType, sType,
                                scalarRows, scalarCols,
                                graph->GetNNE() );
        }
      }
    }
    return retMat;
  }
}
