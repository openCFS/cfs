#include "MatVec/CoordFormat.hh"
#include <cstdio>
#include <cstring>


namespace CoupledField 
{
    //! Allowed Constructor
    template<class TYPE>
    CoordFormat<TYPE>::CoordFormat( UInt nrows, UInt ncols, UInt nnz, bool symmStorage, Layout myLayout ) 
    {
      nrows_ = nrows;
      ncols_ = ncols;
      nnz_   = nnz;

      assemblyDone_ = false;
      myLayout_     = myLayout;

      // Set flag for symmetry and if symmetric storage is requested
      // make sure that matrix is square
      symmStorage_  = symmStorage;
      if ( symmStorage == true && (ncols_ != nrows_) ) {
        EXCEPTION("You requested symmetric storage for a non-square "
                 << "matrix\n"
                 << "Dimensions are " << nrows_ << " x " << ncols_);
      }

      // Since we know the number of entries to be stored in advance,
      // we can "allocate" enough memory for the internal vector
      if ( symmStorage_ == true ) {
        entries_.reserve( (nnz_+nrows)/2 );
      }
      else {
        entries_.reserve( nnz_ );
      }

    }

    //! Constructor which imports the matrix based on a matrix market format file
    template<class TYPE>
    CoordFormat<TYPE>::CoordFormat( const std::string& fname, Layout myLayout)
    {
      //ENTER_FCN( "ImportRealMM" );

         char singleChar;
         char myString[200];
         UInt nrows;
         UInt ncols;
         UInt nentries;
         UInt ridx;
         UInt cidx;
         double value;
         symmStorage_ = false;

         assemblyDone_ = false;
         myLayout_     = myLayout;

         // Open input file and check for errors

         // Open input file and check for errors
         FILE *fp = fopen( fname.c_str(), "r" );
         if ( fp == NULL ) {
           EXCEPTION( "Cannot open file '" << fname << "' for reading");
         }

         // *******************************************************
         //   Read Matrix Market header and validate it is one of
         //
         //   %%MatrixMarket matrix coordinate real general
         //   %%MatrixMarket matrix coordinate real symmetric
         // *******************************************************
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-result"
         fscanf( fp, "%c", &singleChar );
         if ( singleChar != '%' ) {
           EXCEPTION(  "Input file '" << fname << "' seems faulty.\n" //
              << " Problems reading its header!");
         }
         fscanf( fp, "%c", &singleChar );
         if ( singleChar != '%' ) {
           EXCEPTION( "Input file '" << fname << "' seems faulty.\n" //
              << " Problems reading its header!");
         }
         fscanf( fp, "%s", myString );
         if ( strcmp( myString, "MatrixMarket" ) != 0 ) {
           EXCEPTION( "Input file '" << fname << "' seems faulty.\n"
              << " Problems reading its header!");
         }
         fscanf( fp, "%s", myString );
         if ( strcmp( myString, "matrix" ) != 0 ) {
           EXCEPTION( "Input file '" << fname << "' seems faulty.\n"
              << "Contains a '" << myString << "' and not a 'matrix'");
         }
         fscanf( fp, "%s", myString );
         if ( strcmp( myString, "coordinate" ) != 0 ) {
           EXCEPTION("Input file '" << fname << "' seems faulty.\n"
              << "Format is '" << myString << "' and not 'coordinate'");
         }
         fscanf( fp, "%s", myString );
         if ( strcmp( myString, "real" ) != 0 ) {
           EXCEPTION("Input file '" << fname << "' seems faulty.\n"
              << "Entries are of type '" << myString << "' and not 'real'");
         }
         fscanf( fp, "%s", myString );
         if ( strcmp( myString, "general" ) == 0 ) {
           symmStorage_ = false;
         }
         else if ( strcmp( myString, "symmetric") == 0 ) {
           symmStorage_ = true;
         }
         else {
           EXCEPTION("Input file '" << fname << "' seems faulty.\n" //
              << "Format is '" << myString << "' and not 'general' " //
              << "or 'symmetric'");
         }

         // Finish this line (by reading characters until the '\n' appears
         char auxChar;
         do {
           fscanf( fp, "%c", &auxChar );
         } while ( auxChar != '\n' );

         // Read lines, until we have skipped comments
         do {
           fscanf( fp, "%c", &auxChar );
           if ( auxChar == '%' ) {
       do {
         fscanf( fp, "%c", &auxChar );
       } while ( auxChar != '\n' );
           }
           else {
       ungetc( auxChar, fp );
       break;
           }
         } while ( true );

         // This line (hopefully) contains the number of rows, columns and non-zeros
         fscanf( fp, "%u%u%u", &nrows, &ncols, &nentries );

         // Minimal consistency test
         if ( nentries > (nrows * ncols) ) {
           EXCEPTION("Input file '" << fname << "' seems faulty!\n" //
              << "#rows       = " << nrows << ",\n" //
              << "#cols       = " << ncols << ", but\n" //
              << "#numEntries = " << nentries);
         }
         else {
           std::cout << " Input file '" << fname << "' contains a matrix with"
            << "\n  + rows       = " << nrows
            << "\n  + cols       = " << ncols
            << "\n  + numEntries = " << nentries << std::endl;
         }


         UInt nnz = 0;
         if ( symmStorage_ == true ) {
           nnz = ( nentries - nrows ) * 2 + nrows;
           std::cout << "  + format     = symmetric\n";
         }
         else {
           nnz = nentries;
           std::cout << "  + format     = general\n";
         }
         std::cout << "  + nnz        = " << nnz << std::endl;
         this->nnz_ = nnz;
         this->nrows_ = nrows;
         this->ncols_ = ncols;
         entries_.reserve( nnz_ );

         // Read nonzero matrix entries and their indices
         for ( unsigned int k = 0; k < nentries; k++ ) {
           fscanf( fp, "%u%u%lf", &ridx, &cidx, &value );
           this->AddEntry( ridx, cidx, value );
         }
#pragma GCC diagnostic pop

         // The matrix is assembled now
         this->FinaliseAssembly();

         // *********************
         //   Close output file
         // *********************
         if ( fclose( fp ) == EOF ) {
           WARN( "Could not close file %s after reading!" << fname );
         }
    }

    //! Finalises the assembly of the matrix

    //! Call this method after all matrix entries have been set. The object
    //! will then perform some simple consistency checks.
    //! \note
    //! - It is an error to call AddEntry after FinaliseAssembly was called.
    //! - It is an error to call any one of the access methods for matrix
    //!   entries before FinaliseAssembly was called.
    template<class TYPE>
    void CoordFormat<TYPE>::FinaliseAssembly() {

      // Determine that number of stored entries is correct
      UInt numEntries = 0;
      if ( symmStorage_ == true ) {
        numEntries = (nnz_ + ncols_) / 2;
      }
      else {
        numEntries = nnz_;
      }
      if ( numEntries > entries_.size() ) {
        EXCEPTION("Matrix is stored incompletely in CoordFormat object!\n"
                 << " Expected " << numEntries << " entries, but got only "
                 << entries_.size());
      }
      else if  ( numEntries < entries_.size() ) {
    	EXCEPTION( "Matrix is stored incorrectly in CoordFormat object!\n"
                 << " Expected " << numEntries << " entries, but got "
                 << entries_.size());
      }

      // g++ reports syntax error before ;" in the for line
      // Check that all row and column indices are inside bounds
      // for ( std::vector<CoordFormat::triple>::iterator it=entries_.begin();
      //      it != entries_.end(); it++ ) {
      //  if ( it.ridx > nrows_ || it.ridx < 1 ) {
      //    EXCEPTION("Row index of entry out of bounds:\n"
      //           << "row index    = " << it.ridx
      //           << " should be in [1," << nrows_ << "]\n"
      //           << "column index = " << it.cidx
      //             << " should be in [1," << ncols_ << "]\n"
      //             << " value       = " << it.value);
      //  }
      //  if ( it.cidx > ncols_ || it.cidx < 1 ) {
      //    EXCEPTION("Column index of entry out of bounds:\n"
      //             << "row index    = " << it.ridx
      //             << " should be in [1," << nrows_ << "]\n"
      //             << "column index = " << it.cidx
      //             << " should be in [1," << ncols_ << "]\n"
      //             << " value       = " << it.value);
      //  }
      //}

      // Object passed all tests
      assemblyDone_ = true;
    }

    
    //! Comparison operator of for a pure lexicographic ordering

    //! This comparison operator can be used to sort the matrix entries, such
    //! that they are stored in a pure lexicographic ordering. If (i,j) is the
    //! index pair of matrix entry A and (l,m) that of matrix entry B, then
    //! PureLex(A,B) will give true, if either (i<l) or if (i=l) and
    //! (j<m).
   template<class TYPE>
   bool CoordFormat<TYPE>::PureLex( const triple& a, const triple& b ) 
   {
      bool retVal = false;

      if ( a.rowIndex < b.rowIndex ) {
        retVal = true;
      }
      else if ( a.rowIndex > b.rowIndex ) {
        retVal = false;
      }
      else {
        if ( a.colIndex < b.colIndex ) {
          retVal = true;
        }
        else {
          retVal = false;
        }
      }
      return retVal;
  }


// Explicit template instantiations
template class CoordFormat<double>;
template class CoordFormat<std::complex<double>>;  
};
