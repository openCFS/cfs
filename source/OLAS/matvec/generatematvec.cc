// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <complex>

#include <def_use_lapack.hh>

#include "utils/utils.hh"

#include "matvec/generatematvec.hh"

// Include LAPACK stuff
#ifdef USE_LAPACK
#include "external/lapack/lapack.hh"
#include "external/lapack/lapackgbmatrix.cc"
#endif

#include "matvec/matvec.hh"

// include source code for templated vectors
#include "matvec/vector.cc"
//#include "matvec/vector_specialised.cc"

// include source code for templated matrices
#include "matvec/crs_matrix.cc"
#include "matvec/scrs_matrix.cc"
#include "matvec/diag_matrix.cc"

// Include that prematrix stuff
#include "multigrid/prematrix.cc"



namespace OLAS {



  // >>>>>>>>>>>>>>>>>>>>>>>>>>>> VECTOR PART <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<



  // *********************************************
  //   Macro for generation of templated vectors
  // *********************************************
#define CREATE_VECTOR( MATRIX_ENTRY, MATRIX_DOF, VECTORCLASS ) \
if ( ( eType == MATRIX_ENTRY ) && ( blockSize == MATRIX_DOF ) ) {\
   retVector = New VECTORCLASS( length );\
   AssertMem( retVector, sizeof( VECTORCLASS ) );\
   (*cla) << " GenerateSparseVectorObject: Generated "\
          << MACRO2STRING(VECTORCLASS) << " of dimension " << length\
   << std::endl;\
}


  // *******************************************************
  //   Dynamic generation of vector associated with matrix
  // *******************************************************
  BaseVector* GenerateVectorObject( const BaseMatrix &m ) {

    ENTER_FCN("GenerateVectorObject");

    BaseVector *retVector = NULL;
    MatrixStructureType m_structuretype;
    MatrixEntryType m_entrytype;
    MatrixStorageType m_storagetype;
    Integer m_blocksize, m_Ncols;
    const StdMatrix  *m_StdMatrix;
    const SBM_Matrix  *sbmmat;
    SBM_Vector * sbmvec;

    m_structuretype = m.GetStructureType();
    m_Ncols = m.GetNcols();
    
    switch (m_structuretype) {

    case SBM_MATRIX:

      // Downcast matrix to SBM_Matrix
      sbmmat = dynamic_cast<const SBM_Matrix*>(&m);

      // Create an SBM_Vector of fitting length
      sbmvec  = new SBM_Vector( m_Ncols );

      // Create and insert subvectors
      for ( Integer i = 1; i <= m_Ncols; i++ ) {
        if ( sbmmat->GetPointer(i,i) != NULL ) {
          sbmvec->SetSubVector( dynamic_cast<SparseVector*>
                                (GenerateVectorObject((*sbmmat)(i,i))), i );
        }
      }

      // Set vector entries to zero
      sbmvec->Init();

      // Perform an upcast
      retVector = dynamic_cast<BaseVector*>(sbmvec);
      (*cla) << " GenerateVectorObject: Generated an SBM_Vector" << std::endl;

      break;

    case STDMATRIX:

      m_StdMatrix = dynamic_cast<const StdMatrix*>(&m);
      m_entrytype = m_StdMatrix->GetEntryType();
      m_blocksize = m_StdMatrix->GetBlockSize();
      m_storagetype = m_StdMatrix->GetStorageType();

      // Delegate work to factory for sequential SparseVectors
      retVector = GenerateSparseVectorObject( m_storagetype, m_entrytype,
                                           m_blocksize, m_Ncols );
      break;

    default:
      Error( "GenerateVectorObject: Unkown matrix type", __FILE__, __LINE__ );
      break;
    }

    // That's it
    return retVector;  
  }


  // **************************************************
  //   Dynamic generation of vector of specified type
  // **************************************************
  BaseVector* GenerateSparseVectorObject( const MatrixStorageType sType,
                                       const MatrixEntryType eType,
                                       const Integer blockSize,
                                       const Integer length ) {

    ENTER_FCN( "GenerateSparseVectorObject" );

    BaseVector *retVector = NULL;

    // real valued vectors
    CREATE_VECTOR( DOUBLE, 1, Vector<Double>      );

    // complex valued vectors
    CREATE_VECTOR( COMPLEX, 1, Vector<Complex>     );
    
    // scalar vectors to go together with LAPACK matrices
#ifdef USE_LAPACK
    CREATE_VECTOR( F77REAL8    , 1, Vector<Double>  );
    CREATE_VECTOR( F77COMPLEX16, 1, Vector<Complex> );
#endif
  
    // Check, if we were able to generate a vector object. If not, complain.
    if ( retVector == NULL ) {
      Error( "GenerateVectorObject: Failed to generate vector object",
             __FILE__, __LINE__ );
    }

    // Force instantiation of member functions of templated vector classes.
    // Note that the InstantiatePublicMethods() method itself should never
    // actually be called, but the compiler must not notice this. Thus, the
    // akward if condition.
    bool neverTrue = ( sType == NOSTORAGETYPE );
    if ( neverTrue ) {
      Error( "GenerateVectorObject: Wrong input! sType = NOSTORAGETYPE",
             __FILE__, __LINE__ );
      SparseVector *stdVec = dynamic_cast<SparseVector*>(retVector);
      if ( stdVec != NULL ) {
        stdVec->InstantiatePublicMethods();
      }
    }

    // Initialise vector to zero
    retVector->Init();

    // That's it
    return retVector;  
  }


  // >>>>>>>>>>>>>>>>>>>>>>>>>>>> MATRIX PART <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<



  // ***************************************************
  //   Central macro for generation of a matrix object
  //   Its use also takes care of the instantiation in
  //   the case of the templated matrices
  // ***************************************************
#define MATRIX_OBJ(matEntry,matStore,numDof,matrix_obj_type)\
if ((etype==matEntry) && (stype==matStore) && (dof==numDof)){\
retMat = New matrix_obj_type(nrows,ncols,fill);\
AssertMem( retMat, sizeof(matrix_obj_type) );\
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
  StdMatrix* GenerateStdMatrixObject( const MatrixEntryType etype,
                                      const MatrixStorageType stype,
                                      const Integer dof,
                                      const Integer nrows,
                                      const Integer ncols,
                                      const Integer fill ) {

    ENTER_FCN( "GenerateStdMatrixObject" );
    StdMatrix *retMat = NULL;

    // Branch depending on desired storage layout/StdMatrix type
    switch( stype ) {

      // CRS Matrices
    case SPARSE_NONSYM:
      MATRIX_OBJ( DOUBLE,  SPARSE_NONSYM, 1, CRSRealDof1    );
      MATRIX_OBJ( COMPLEX, SPARSE_NONSYM, 1, CRSComplexDof1 );
      break;

      // SCRS Matrices
    case SPARSE_SYM:
      MATRIX_OBJ( DOUBLE,  SPARSE_SYM, 1, SCRSRealDof1    );
      MATRIX_OBJ( COMPLEX, SPARSE_SYM, 1, SCRSComplexDof1 );
      break;

      // NSK Matrices
    case SKYLINE_NONSYM:
      Error( "GenerateStdMatrixObject: SKYLINE_NOSYM not supported, yet!",
             __FILE__, __LINE__ );
      break;

      // SSK Matrices
    case SKYLINE_SYM:
      Error( "GenerateStdMatrixObject: SSK Matrix not yet supported",
             __FILE__, __LINE__ );
      // MATRIX_OBJ( DOUBLE,  SKYLINE_SYM, 1, SSKRealDof1    );
      // MATRIX_OBJ( COMPLEX, SKYLINE_SYM, 1, SSKComplexDof1 );
      break;

      // Diag Matrices
    case DIAG:
      MATRIX_OBJ( DOUBLE,  DIAG, 1, DiagRealDof1    );
      MATRIX_OBJ( COMPLEX, DIAG, 1, DiagComplexDof1 );
      break;


#ifdef USE_LAPACK
    case LAPACK_GBMATRIX:

      // real entries
      if ( etype == DOUBLE && dof == 1 ) {
        retMat = New LapackGBMatrix<F77real8,double>( nrows, ncols,
                                                      F77REAL8 );
        AssertMem( retMat, sizeof(LapackGBMatrix<F77real8,double>) );
        (*cla) << " Generated matrix of type: LapackGBMatrix<F77Double>"
               << std::endl;
      }

      // complex entries
      else if ( etype == COMPLEX && dof == 1 ) {
        retMat = New LapackGBMatrix<F77complex16,Complex>( nrows, ncols,
                                                           F77COMPLEX16 );
        AssertMem(retMat, sizeof(LapackGBMatrix<F77complex16,Complex>) );
        (*cla) << " Generated matrix of type: LapackGBMatrix<F77Complex16>"
               << std::endl;
      }

      // nothing fits
      else {
        Error( "Unable to generate LapackGBMatrix", __FILE__, __LINE__ );
      }
      break;

#else

    case LAPACK_GBMATRIX:
      Error( "Compile with USE_LAPACK to enable support for LapackGBMatrix",
             __FILE__, __LINE__ );
      break;

#endif

      // Should not be reached (hopefully :)
    default:
      Error( "GenerateStdMatrixObject: Request for unknown MatrixStorageType",
             __FILE__, __LINE__ );
    }


    // Force instantiation of member functions of templated matrix classes
    // Note that the InstantiatePublicMethods() method itself should never
    // actually be called, but the compiler must not notice this. Thus, the
    // akward if condition.
    bool neverTrue = (dof == -2);

    if ( neverTrue ) {

      // Check for templated matrix
      if ( stype == SPARSE_SYM  || stype == SPARSE_NONSYM ||
           stype == SKYLINE_SYM || stype == SKYLINE_NONSYM ||
           stype == LAPACK_GBMATRIX  || stype == DIAG ) {
        retMat->InstantiatePublicMethods();
      }
    }

    // Return pointer to generated matrix object
    return retMat;

  }


  // ***************************************************
  //   Central macro for generation of a matrix object
  //   by use of a copy constructor
  // ***************************************************
#define COPY_MATRIX_OBJ( matEntry, matStore, numDof, matrixObjType )\
try {\
   ConstRefCast( origMat, matrixObjType, auxMat );\
   retMat = New matrixObjType( auxMat );\
   AssertMem( retMat, sizeof( matrixObjType) );\
   (*cla) << " Generated copy of matrix of type: "\
          << MACRO2STRING( matrixObjType ) \
          << std::endl;\
}\
catch(...){};



  // ******************************************************
  //   Dynamic generation of a matrix object (by copying)
  // ******************************************************
  StdMatrix* CopyStdMatrixObject( const StdMatrix &origMat  ) {

    ENTER_FCN( "CopyStdMatrixObject" );

    StdMatrix *retMat = NULL;

    // Determine matrix information
    MatrixEntryType eType = origMat.GetEntryType();
    MatrixStorageType sType = origMat.GetStorageType();
    Integer bSize = origMat.GetBlockSize();

    // CRS_Matrix case
    COPY_MATRIX_OBJ( DOUBLE,  SPARSE_NONSYM, 1, CRSRealDof1    );
    COPY_MATRIX_OBJ( COMPLEX, SPARSE_NONSYM, 1, CRSComplexDof1 )

    // SCRS_Matrix case
    COPY_MATRIX_OBJ( DOUBLE,  SPARSE_SYM, 1, SCRSRealDof1    );
    COPY_MATRIX_OBJ( COMPLEX, SPARSE_SYM, 1, SCRSComplexDof1 );

    // DIAG_Matrix case
    COPY_MATRIX_OBJ( DOUBLE,  DIAG, 1, SCRSRealDof1    );
    COPY_MATRIX_OBJ( COMPLEX, DIAG, 1, SCRSComplexDof1 );

    // Invalid or unimplemented case
    if ( retMat == NULL ) {
      (*error) << "GenerateStdMatrixObject: Failed to generate StdMatrix "
               << "by copying an StdMatrix of type ( sType = "
               << Enum2String( sType ) << " , eType = "
               << Enum2String( eType ) << " , bSize = "
               << bSize;
      Error( __FILE__, __LINE__ );
    }


    return retMat;
  }


  // *********************************************
  //   Dynamic generation of a SBM matrix object
  // *********************************************
  SBM_Matrix* GenerateSBMMatrixObject( MatrixEntryType eType,
                                       UInt rowDim, UInt colDim,
                                       GraphManagerSBMMat *graphManager,
                                       bool symmetric ) {

    ENTER_FCN( "GenerateSBMMatrixObject" );

    MatrixStorageType sType = NOSTORAGETYPE;
    BaseGraph *graph = NULL;
    UInt scalarRows = 0;
    UInt scalarCols = 0;

    // STEP 1: Generate empty SBM_Matrix
    SBM_Matrix *retMat = NULL;
    retMat = New SBM_Matrix( rowDim, colDim, symmetric );
    if ( retMat == NULL ) {
      (*error) << "GenerateSBM_MatrixObject: "
               << "This is the end my friend!\n"
               << "Generation of empty SBM_Matrix failed!";
      Error( __FILE__, __LINE__ );
    }

    // STEP 2: Populate with sub-matrices
    for ( UInt i = 1; i <= rowDim; i++ ) {
      for ( UInt j = 1; j <= colDim; j++ ) {

        // Check for need of sub-matrix
        if ( graphManager->SubGraphExists(i,j) ) {

          // Determine storage type
          if ( symmetric == true && i == j ) {
            sType = SPARSE_SYM;
          }
          else {
            sType = SPARSE_NONSYM;
          }

          // Get hold of graph of sub-matrix
          graph = graphManager->GetGraph(i,j);

          // Determine dimension of sub-matrix
          scalarRows = graphManager->GetGraph(i,i)->GetSize();
          scalarCols = graphManager->GetGraph(j,j)->GetSize();

          // Generate empty sub-matrix
          retMat->SetSubMatrix( i, j, eType, sType, 1,
                                scalarRows, scalarCols,
                                graph->GetNNE() );
        }
      }
    }
    return retMat;
  }

    // Explicitly instantiate Vector<> templates
  template class Vector<Double>;
  template class Vector<Complex>;  
}
