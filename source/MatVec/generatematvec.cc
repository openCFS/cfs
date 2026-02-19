#include <complex>

// Include LAPACK stuff
#include "OLAS/external/lapack/lapack.hh"
#include "OLAS/external/lapack/LapackGBMatrix.hh"

#include "Vector.hh"

// include source code for templated matrices
#include "CRS_Matrix.hh"
#include "SCRS_Matrix.hh"
#include "Diag_Matrix.hh"
#include "VBR_Matrix.hh"
#include "SBM_Matrix.hh"

#include "SBM_Vector.hh"
#include "DataInOut/Logging/LogConfigurator.hh"

namespace CoupledField {

// define logging stream
DEFINE_LOG(genMatVec, "genMatVec")

  // >>>>>>>>>>>>>>>>>>>>>>>>>>>> VECTOR PART <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


  // *******************************************************
  //   Dynamic generation of vector associated with matrix
  // *******************************************************
BaseVector* GenerateVectorObject( const BaseMatrix &m, BaseMatrix::EntryType entrytype ) {

    BaseVector *retVector = NULL;
    BaseMatrix::StructureType m_structuretype;
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
                                (GenerateVectorObject((*sbmmat)(i,i), entrytype)), i );
        }
      }

      // Set vector entries to zero
      sbmvec->Init();

      // Perform an upcast
      retVector = dynamic_cast<BaseVector*>(sbmvec);
      LOG_DBG(genMatVec) << " GenerateVectorObject: Generated an SBM_Vector";

      break;

    case BaseMatrix::SPARSE_MATRIX:
      m_StdMatrix = dynamic_cast<const StdMatrix*>(&m);
      if(entrytype == BaseMatrix::NOENTRYTYPE) 
      {
        entrytype = m_StdMatrix->GetEntryType();
      }
      

      // Delegate work to factory for sequential SingleVectors
      retVector = GenerateSingleVectorObject( entrytype,
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
  SingleVector* GenerateSingleVectorObject( const BaseMatrix::EntryType eType,
                                            const UInt length ) {


    SingleVector *retVector = NULL;

    if(eType == BaseMatrix::DOUBLE)
      retVector = new Vector<Double>(length);
    if(eType == BaseMatrix::COMPLEX)
      retVector = new Vector<Complex>(length);

    // scalar vectors to go together with LAPACK matrices
    if(eType == BaseMatrix::F77REAL8)
      retVector = new Vector<Double>(length);
    if(eType == BaseMatrix::F77COMPLEX16)
      retVector = new Vector<Complex>(length);

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
  }\
  LOG_DBG(genMatVec) << " GenerateSingleVectorObject: Generated copy of" \
         << MACRO2STRING(VECTORCLASS);
  
  SingleVector* CopySingleVectorObject( const SingleVector& origVec ) {

    SingleVector *retVector = NULL;

    // Determine vector information
    BaseMatrix::EntryType eType = origVec.GetEntryType();

    // real valued vectors
    COPY_VECTOR( BaseMatrix::DOUBLE,  Vector<Double>      );

    // complex valued vectors
    COPY_VECTOR( BaseMatrix::COMPLEX,  Vector<Complex>     );

    // scalar vectors to go together with LAPACK matrices
    COPY_VECTOR( BaseMatrix::F77REAL8    ,  Vector<Double>  );
    COPY_VECTOR( BaseMatrix::F77COMPLEX16,  Vector<Complex> );

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
LOG_DBG(genMatVec) << " Generated matrix of type: "<< MACRO2STRING(matrix_obj_type);}


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

  // ****************
  //   VBR Matrices
  // ****************
  typedef VBR_Matrix<Double>       VBRRealDof1;
  typedef VBR_Matrix<Complex>      VBRComplexDof1;

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
      
      // VBR Matrices
    case BaseMatrix::VAR_BLOCK_ROW:
      MATRIX_OBJ( BaseMatrix::DOUBLE,  BaseMatrix::VAR_BLOCK_ROW, VBRRealDof1    );
      MATRIX_OBJ( BaseMatrix::COMPLEX, BaseMatrix::VAR_BLOCK_ROW, VBRComplexDof1 );
      break;

      // Diag Matrices
    case BaseMatrix::DIAG:
      MATRIX_OBJ( BaseMatrix::DOUBLE,  BaseMatrix::DIAG, DiagRealDof1    );
      MATRIX_OBJ( BaseMatrix::COMPLEX, BaseMatrix::DIAG, DiagComplexDof1 );
      break;

    case BaseMatrix::LAPACK_GBMATRIX:

      // real entries
      if ( etype == BaseMatrix::DOUBLE  ) {
        retMat = new LapackGBMatrix<double,double>( nrows, ncols,
                                                      BaseMatrix::F77REAL8 );
        ASSERTMEM( retMat, sizeof(LapackGBMatrix<double,double>) );
        LOG_DBG(genMatVec) << " Generated matrix of type: LapackGBMatrix<F77Double>";
      }

      // complex entries
      else if ( etype == BaseMatrix::COMPLEX  ) {
        retMat = new LapackGBMatrix<Complex,Complex>( nrows, ncols,
        		BaseMatrix::F77COMPLEX16 );
        ASSERTMEM(retMat, sizeof(LapackGBMatrix<Complex,Complex>) );
        LOG_DBG(genMatVec) << " Generated matrix of type: LapackGBMatrix<F77Complex16>";
      }

      // nothing fits
      else {
        EXCEPTION( "Unable to generate LapackGBMatrix" );
      }
      break;


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
   LOG_DBG(genMatVec) << " Generated copy of matrix of type: "\
          << MACRO2STRING( matrixObjType ); \
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
    COPY_MATRIX_OBJ( DOUBLE,  SPARSE_SYM, SCRSRealDof1      );
    COPY_MATRIX_OBJ( COMPLEX, SPARSE_SYM, SCRSComplexDof1   );
    
    // VBR_Matrix case
    COPY_MATRIX_OBJ( DOUBLE,  VAR_BLOCK_ROW, VBRRealDof1    );
    COPY_MATRIX_OBJ( COMPLEX, VAR_BLOCK_ROW, VBRComplexDof1 );
    
    // DIAG_Matrix case
    COPY_MATRIX_OBJ( DOUBLE,  DIAG, DiagRealDof1    );
    COPY_MATRIX_OBJ( COMPLEX, DIAG, DiagComplexDof1 );
    
    

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
    BaseGraph* graph = NULL;
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
