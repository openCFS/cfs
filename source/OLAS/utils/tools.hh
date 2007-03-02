#ifndef OLAS_TOOLS_CLA
#define OLAS_TOOLS_CLA

#include "utils/defs.hh"
#include "utils/environment.hh"

// We need to be able to abort MPI
// in parallel setting
#ifdef PARALLEL
#include "graph/baseparallel.hh"
#endif

//! \file tools.hh
//! This file contains routines for issuing error and warning messages

namespace OLAS {
  //! Prints out an error message with according filename
  //! and line number. The program execution will be
  //! stopped.
  void Error( const Char*         Text,
              const Char*   const filename,
              const Integer       numline );

  //! Write error message and terminate program execution

  //! This method can be used to output an error message and terminate program
  //! execution. The error message will be taken from the global stringstream
  //! pointer error. If error is NULL, then a default message will be
  //! substituted. Apropriate usage of this method looks like this
  //! \code
  //! (*error) << "Cannot open file '" << fname << "' for reading";
  //! Error( __FILE__, __LINE__ );
  //! \endcode
  void Error( const Char *const filename, const Integer numline );

  //! Overloaded version of the Error function, which accepts
  //! a format string, analoguously to the C-function printf.
  void Error( const Char*   const filename,
              const Integer       numline,
              const Char*   const formatstring, ... );

  //! Write warning message and continue program execution

  //! This method can be used to output a warning message and continue program
  //! execution. The warning message will be taken from the global stringstream
  //! pointer warning. If warning is NULL, then a default message will be
  //! substituted. Apropriate usage of this method looks like this
  //! \code
  //! (*warning) << "Delta is getting very small, delta = " << delta << '\n';
  //! Warning( __FILE__, __LINE__ );
  //! \endcode
  //! The warning will be written to std::cerr as well as to the standard OLAS
  //! logfile (*cla).
  void Warning( const Char *const filename, const Integer numline );

  //! Prints out a warning  message with according filename
  //! and line number. The program execution will be
  //! continued.
  void Warning( const Char*         Text,
                const Char*   const filename = NULL,
                const Integer       numline  = 0 );

  //! Overloaded version of the Warning function, which accepts
  //! a format string, analoguously to the C-function printf.
  void Warning( const Char*   const filename,
                const Integer       numline,
                const Char*   const formatstring, ... );  
} // namespace

#endif // OLAS_TOOLS_HH
