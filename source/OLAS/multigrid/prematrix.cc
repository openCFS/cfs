// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "multigrid/prematrix.hh"

/**********************************************************/
#ifdef DEBUG_TO_CERR
#ifndef DEBUG_PREMATRIX
#define DEBUG_PREMATRIX
#endif // DEBUG_PREMATRIX
#define  debug  &std::cerr
#endif // DEBUG_PREMATRIX

#ifdef PROFILING
#ifdef PROFILE_PREMATRIX
#include "utils/utils.hh"
#endif
#endif
/**********************************************************/
#define  PREMATRIX_ROW_GRANULARITY  32
/**********************************************************/

namespace OLAS {
/**********************************************************/

template <typename T>
PreMatrix<T>::PreMatrix()
            : rows_( NULL ),
              cols_( NULL ),
              rowSize_( NULL ),
#ifndef PREMATRIX_USE_MODULO
              buffSize_( NULL ),
#endif
              numRows_( 0 ),
              hasRowsSorted_(false),
              hasDiagFirst_(false)
{
}

template <typename T>
PreMatrix<T>::~PreMatrix()
{
    
    Reset();
}

/**********************************************************/

template <typename T>
inline Integer PreMatrix<T>::GetNumNonzeros() const
{
    Integer nnz = 0;
    for( Integer i = 1; i <= numRows_; i++ )  nnz += rowSize_[i];
    return nnz;
}

/**********************************************************/

template <typename T>
inline Integer PreMatrix<T>::GetRowSize( const Integer row ) const
{
#ifdef DEBUG_PREMATRIX
    if( row < 1 || row > GetNrows() ) {
        Error( __FILE__, __LINE__,
               "PreMatrix::GetRowSize(%d): %d is not a valid "
               "row index in range [1,%d]",
               row, row, GetNrows() );
    }
#endif

    return rowSize_[row];
}

/**********************************************************/

template <typename T>
inline const T* PreMatrix<T>::GetRowData( const Integer row ) const
{
#ifdef DEBUG_PREMATRIX
    if( row < 1 || row > GetNrows() ) {
        Error( __FILE__, __LINE__,
               "PreMatrix::GetRowData(%d): %d is not a valid "
               "row index in range [1,%d]",
               row, row, GetNrows() );
    }
#endif

    return rows_[row];
}

/**********************************************************/

template <typename T>
inline const Integer* PreMatrix<T>::
GetRowCols( const Integer row ) const
{
#ifdef DEBUG_PREMATRIX
    if( row < 1 || row > GetNrows() ) {
        Error( __FILE__, __LINE__,
               "PreMatrix::GetRowCols(%d): %d is not a valid "
               "row index in range [1,%d]",
               row, row, GetNrows() );
    }
#endif

    return cols_[row];
}

template <typename T>
inline Integer PreMatrix<T>::
GetNrows() const
{
    return numRows_;
}

/**********************************************************/

template <typename T>
inline T PreMatrix<T>::GetEntry( const Integer i,
                        const Integer j ) const
{
    const Integer* const cols = GetRowCols(i);
    for( Integer ij = 1; ij <= GetRowSize(i); ij++ ) {
        if( cols[ij] == j ) {
            return GetRowData(i)[ij];
        }
    }
    return (T)0;
}

/**********************************************************/

template <typename T>
void PreMatrix<T>::SetEntry( const Integer  row, const Integer col,
                             const T& value )
{

#ifdef DEBUG_PREMATRIX
    if( row < 1 || GetNrows() < row ) {
        Error( __FILE__, __LINE__,
               "PreMatrix::SetEntry(%d,%d,..): %d is not a valid "
               "row index in range [1,%d]",
               row, col, row, GetNrows() );
    }
#endif

    // find eventually allready present entry
    for( Integer ic = 1; ic <= rowSize_[row]; ic++ ) {
        if( cols_[row][ic] == col ) {
            rows_[row][ic] = value;
            return;
        }
    }

    // entry not yet present
    CheckBufferSize( row );
    // add new entry
    rows_[row][++rowSize_[row]] = value;
    cols_[row][  rowSize_[row]] = col;
}

/**********************************************************/

template <typename T>
void PreMatrix<T>::AddToEntry( const Integer  row,
                               const Integer  col,
                               const T& value )
{

#ifdef DEBUG_PREMATRIX
    if( row < 1 || row > GetNrows() ) {
        Error( __FILE__, __LINE__,
               "PreMatrix::AddToEntry(%d,%d,..): %d is not a valid "
               "row index in range [1,%d]",
               row, col, row, GetNrows() );
    }
#endif

    // find eventually allready present entry
    for( Integer ic = 1; ic <= rowSize_[row]; ic++ ) {
        if( cols_[row][ic] == col ) {
            rows_[row][ic] += value;
            return;
        }
    }

    // entry not yet present
    CheckBufferSize( row );
    // add new entry
    rows_[row][++rowSize_[row]] = value;
    cols_[row][  rowSize_[row]] = col;
}

/**********************************************************/

template <typename T>
void PreMatrix<T>::SortDiagonal()
{

    T tdat;
    for( Integer i = 1; i <= numRows_; i++ ) {
        // search the diagonal entry
        for( Integer ij = 1; ij <= rowSize_[i]; ij++ ) {
            // place the diagonal entry at leading position
            if( cols_[i][ij] == i ) {
                // swap the column indices
                cols_[i][ij] = cols_[i][1];
                cols_[i][1]  = i;
                // swap the data
                tdat         = rows_[i][ij];
                rows_[i][ij] = rows_[i][1];
                rows_[i][1]  = tdat;
                // jump to next row
                break;
            }
#ifdef  DEBUG_PREMATRIX
            // if we reach this part of the loop in the last
            // loop run, this means that we could not find the
            // diagonal entry
            if( ij == rowSize_[i] ) {
                Warning( __FILE__, __LINE__,
                         "PreMatrix::SortDiagonal: found row [%d] "
                         "without diagonal entry", i );
                         
            }
#endif
        }
    }

    hasDiagFirst_ = true;
}

/**********************************************************/

template <typename T>
void PreMatrix<T>::Sort()
{

    hasDiagFirst_ = true;

    T tdat;
    for( Integer i = 1; i <= numRows_; i++ ) {
//std::cerr
//<< "sorting row "<<i<<" / "<<numRows_<<": size = "<<rowSize_[i]<<std::endl;
        // search the diagonal entry
        Integer ij;
        for( ij = 1; ij <= rowSize_[i]; ij++ ) {
            // place the diagonal entry at leading position
            if( cols_[i][ij] == i ) {
                // swap the column indices
                cols_[i][ij] = cols_[i][1];
                cols_[i][1]  = i;
                // swap the data
                tdat         = rows_[i][ij];
                rows_[i][ij] = rows_[i][1];
                rows_[i][1]  = tdat;
                // sort the remaining row using quick sort
                QuickSort( cols_[i]+2, rows_[i]+2, rowSize_[i]-1 );
                // next row
                break;
            }
        }
        // if we could not find a diagonal entry in this row,
        // simply sort the whole row
        if( ij > rowSize_[i] ) {
            QuickSort( cols_[i]+1, rows_[i]+1, rowSize_[i] );
            hasDiagFirst_ = false;
        }
    }

    hasRowsSorted_ = true;
}

/**********************************************************/

template <typename T>
void PreMatrix<T>::Reset()
{

    // delete array with row data
    if( rows_ != NULL ) {
        for( Integer i = 1; i <= numRows_; i++ ) {
            DeleteArray( rows_[i] ); }
        DeleteArray( rows_ );  rows_ = NULL;        
    }
    // delete array with column data
    if( cols_ != NULL ) {
        for( Integer i = 1; i <= numRows_; i++ ) {
            DeleteArray( cols_[i] ); }
        DeleteArray( cols_ );  cols_ = NULL;        
    }
    // delete array with row sizes
    DeleteArray( rowSize_ );  rowSize_ = NULL;
#ifndef PREMATRIX_USE_MODULO
    DeleteArray( buffSize_ );  buffSize_ = NULL;
#endif
    // reset number of rows
    numRows_ = 0;
    
    hasDiagFirst_  = false;
    hasRowsSorted_ = false;
}

/**********************************************************/
#ifdef PREMATRIX_ENABLE_IMPORT

template <typename T>
bool PreMatrix<T>::ImportASCII( const char* const filename )
{

    FILE    *file   = fopen( filename, "rt" );
    char    *line   = NULL;
    size_t   length = 0;
    ssize_t  read   = 0;
    
    if( file == NULL ) {
        Warning( __FILE__, __LINE__,
                 "PreMatrix::importASCII: could not open file "
                 "\"%s\". Matrix stayes unchanged.",
                 filename );
        return false;
    }

    int nrows = -1,
        ncols = -1,
        nnz   = -1;

    while( -1 != (read = getline(&line, &length, file)) ) {
        // skip comments in the file header. Here we restrict to
        // recognize comments, that start with a % character.
        if( line[0] == '%' )  continue;
        // read matrix data
        if( nrows < 0 ) {
            if( 3 == sscanf(line, "%d %d %d", &nrows, &ncols, &nnz) ) {
                if( !SetNumRows(nrows) ) {
                    fclose( file );
                    return false;
                }
            } else {
                Warning( __FILE__, __LINE__,
                         "PreMatrix::importASCII: could not read matrix "
                         "header data in file \"%s\". File might not have"
                         " the appropriate format.", filename );
                fclose( file );
                return false;
            }
        // read matrix entry
        } else {
            int     row, col;
            T val;
            if( ImportEntryASCII(line, row, col, val) ) {
                SetEntry( row, col, val );
                // read las entry
                if( --nnz <= 0 )  break;
            } else {
                Warning( __FILE__, __LINE__,
                         "PreMatrix::importASCII: Could not read matrix "
                         "entry in file \"%s\". Still %d non-zero entries "
                         "are missing. Aborting import.",
                         filename, nnz );
                fclose( file );
                return false;
            }
        }
    }

    if( line )  free( line );
    fclose( file );

    return true;
}

#endif // PREMATRIX_ENABLE_IMPORT
/**********************************************************/

template <typename T>
std::ostream& PreMatrix<T>::Print( std::ostream& out ) const
{

    out << "\033[1mPreMatrix:\033[0m" << std::endl
        << "  "<<GetNrows()<<" rows " << std::endl
        << "  "<<GetNumNonzeros()<<" entries: " << std::endl;
    
    for( Integer r = 1; r <= GetNrows(); r++ ) {
        for( Integer ic = 1; ic <= GetRowSize(r); ic++ ) {
            out << "("<<r<<","<<GetRowCols(r)[ic]<<") = "
                << GetRowData(r)[ic] << std::endl;
        }
    }
    
    return out;
}

/**********************************************************/

template <typename T>
bool PreMatrix<T>::CreateArrays( const Integer num_rows )
{

    if( !num_rows ) { Reset();  return true; }

    // delete rows and the rows array, if it cannot be reused
    if( rows_ != NULL ) {
        for( Integer i = 1; i <= numRows_; i++ ) {
            DeleteArray( rows_[i] );  rows_[i] = NULL; }
        if( num_rows != numRows_ ) {
            DeleteArray( rows_ );  rows_ = NULL; }
    }
    // new rows array
    if( rows_ == NULL ) {
        NewArray( rows_, T*, num_rows );
        for( Integer i = 1; i <= num_rows; i++ )  rows_[i] = NULL;
    }

    // delete column indices and columns array, if it cannot be reused
    if( cols_ != NULL ) {
        for( Integer i = 1; i <= numRows_; i++ ) {
            DeleteArray( cols_[i] );  cols_[i] = NULL; }
        if( num_rows != numRows_ ) {
            DeleteArray( cols_ );  cols_ = NULL; }
    }
    // new column array
    if( cols_ == NULL ) {
        NewArray( cols_, Integer*, num_rows );
        for( Integer i = 1; i <= num_rows; i++ )  cols_[i] = NULL;
    }

    // new array for row sizes
    if( num_rows != numRows_ ) {
        DeleteArray( rowSize_ );
        NewArray( rowSize_, Integer, num_rows );
#ifndef PREMATRIX_USE_MODULO
        DeleteArray( buffSize_ );
        NewArray( buffSize_, Integer, num_rows );
#endif
    }
    // must be initialized or reset in any case
    for( Integer i = 1; i <= num_rows; i++ )  rowSize_[i] = 0;
#ifndef PREMATRIX_USE_MODULO
    for( Integer i = 1; i <= num_rows; i++ )  buffSize_[i] = 0;
#endif
    numRows_ = num_rows;

    return true;
}

/**********************************************************/

template <typename T>
void PreMatrix<T>::InstantiatePublicMethods()
{
    PreMatrix<T> dummyPreMatrix( 1 );
    GetNumNonzeros();
    GetNrows();
    GetRowSize( 1 );
    GetRowData( 1 );
    GetRowCols( 1 );
    GetEntry( 1, 1 );
    SetNumRows( 1 );
    SetEntry( 1, 1, (T)1 );
    AddToEntry( 1, 1, (T)1 );
    SortDiagonal();
    Sort();
    Reset();
#ifdef PREMATRIX_ENABLE_IMPORT
    ImportASCII( NULL );
#endif
    Print( std::cout );
}

/**********************************************************/

template <typename T>
inline void PreMatrix<T>::CheckBufferSize( const Integer row )
{

#ifdef PREMATRIX_USE_MODULO
    // enlarge row, if necessary
    if( rowSize_[row] % PREMATRIX_ROW_GRANULARITY == 0 ) {
#else
    if( rowSize_[row] == buffSize_[row] ) {
#endif
        // enlarge data array
        T *newRow = NULL;
#ifdef PREMATRIX_USE_MODULO
        NewArray( newRow, T,
                  rowSize_[row] + PREMATRIX_ROW_GRANULARITY );
#else
        NewArray( newRow, T,
                  buffSize_[row] += PREMATRIX_ROW_GRANULARITY );
#endif
        for( Integer i = 1; i <= rowSize_[row]; i++ ) {
            newRow[i] = rows_[row][i]; }
        DeleteArray( rows_[row] );
        rows_[row] = newRow;
        // enlarge column array
        Integer *newCols = NULL;
#ifdef PREMATRIX_USE_MODULO
        NewArray( newCols, Integer,
                  rowSize_[row] + PREMATRIX_ROW_GRANULARITY );
#else
        NewArray( newCols, Integer, buffSize_[row] );
#endif
        for( Integer i = 1; i <= rowSize_[row]; i++ ) {
            newCols[i] = cols_[row][i]; }
        DeleteArray( cols_[row] );
        cols_[row] = newCols;
    }
}

/**********************************************************/

template <typename T>
void PreMatrix<T>::QuickSort(       Integer *const cols,
                                    T *const data,
                              const Integer        length )
{

    if( length <= 1 )  return;

    const Integer last         = length-1,
                  splitter     = cols[last];
    const T splitterData = data[last];

    Integer tcol, i = 0;
    T tdat;

    // splitt array
    for( Integer j = 0; j < last; j++ ) {
        if( cols[j] < splitter ) {
            // swap the column index
            tcol    = cols[j];
            cols[j] = cols[i];
            cols[i] = tcol;
            // swap the data
            tdat    = data[j];
            data[j] = data[i];
            data[i] = tdat;
            i++;
        }
    }
    // swap the splitter
    cols[last] = cols[i];
    cols[i]    = splitter;
    // and the corresponding data
    data[last] = data[i];
    data[i]    = splitterData;

    // quicksort the two parts
    QuickSort( cols,   data,   i++        );
    QuickSort( cols+i, data+i, length - i );
}

/**********************************************************/
#ifdef PREMATRIX_ENABLE_IMPORT

template <typename T>
bool PreMatrix<T>::
ImportEntryASCII( const char* const line,
                  Integer &row, Integer &col, T &val )
{
    return 3 == sscanf( line, "%d %d %lg", &row, &col, &val );
}

#endif // PREMATRIX_ENABLE_IMPORT
/**********************************************************/
#undef  PREMATRIX_ROW_GRANULARITY
/**********************************************************/
} // namespace OLAS

/**********************************************************
 * print-out operator
 **********************************************************/

namespace std {
template <typename T> std::ostream& operator <<
( std::ostream& out, const OLAS::PreMatrix<T>& A ) {
    return A.Print( out ); }
}

/**********************************************************/
