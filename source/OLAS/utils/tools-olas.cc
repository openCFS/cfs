// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;


#include <iostream>
#include <string>

// C-includes
#include <cstdio>
#include <cstdarg>

#include "utils/tools.hh"
#include "utils/output.hh"

namespace OLAS {

  // *********************
  //   Error (version 1)
  // *********************
  void Error( const Char*          Text,
              const Char*    const filename,
              const Integer        numline ) {

    // Write text to error stream an call Error (version 2)
    (*error) << Text;
    Error( filename, numline );
  }


  // *********************
  //   Error (version 2)
  // *********************
  void Error( const Char *const filename, const Integer numline ) {

    // Flush output streams
#ifdef DEBUG
    if ( debug != NULL ) {
      (*debug).flush();
    }
#endif
    std::cout << std::flush;
    std::cerr << std::flush;

    // Obtain error message and clear string stream
    std::string errMsg = "";
    if ( error != NULL ) {
      errMsg = error->str();
      error->str( "" );
    }

    // See, if we could obtain an error message
    if ( errMsg == "" ) {
      errMsg = "Sorry, but no message was provided by the coder!";
    }

    // Cry out to standard error
    std::cerr << "\n\n \033[31mERROR:\033[0m\n\n "
	      << errMsg;

    // Add name of file and line number where the error
    // message was triggered
    if ( filename ) {
      std::cerr << "\n\n This error message was brought to you by\n "
		<< filename << ", line " << numline;
    }

    // Create some space before the user prompt re-appears
    std::cerr << "\n\n";

#ifdef PARALLEL
    OLAS_MPI_Abort(OLAS_MPI_COMM_WORLD,-1);
#endif

    exit(-1);
  }


  // *********************
  //   Error (version 3)
  // *********************
  void Error( const Char*   const filename,
              const Integer       numline,
              const Char*   const formatstring, ... ) {

    Integer  buffSize = 255,
             msgSize  = buffSize;
    char    *message  = NULL;

    // Since msgSize was initialized with buffSize,
    // this loop is entered at least once.
    while( buffSize <= msgSize ) {
        // create buffer for message
        delete [] message;
        message = New char [buffSize = msgSize+1];
        // try to format the message into the buffer
        std::va_list arg;
        va_start( arg, formatstring );
        msgSize = vsnprintf( message, buffSize, formatstring, arg );
        va_end( arg );
    }

    (*error) << message;
    delete [] message;
    Error( filename, numline );
  }


  // ***********************
  //   Warning (version 1)
  // ***********************
  void Warning( const Char*         Text,
                const Char*   const filename,
                const Integer       numline ) {

    // Write text to warning stream an call Warning (version 2)
    (*warning) << Text;
    Warning( filename, numline );
  }


  // ***********************
  //   Warning (version 2)
  // ***********************
  void Warning( const Char *const filename, const Integer numline ) {

    // Obtain warning message and clear string stream
    std::string warnMsg = "";
    if ( warning != NULL ) {
      warnMsg = warning->str();
      warning->str( "" );
    }
    else {
      warnMsg  = " A warning was triggered, but no message supplied.\n";
      warnMsg += " It occured in ";
    }

    // Write warning to standard error
    std::cerr << "\n\n \033[31mWARNING:\033[0m\n\n " << warnMsg
	      << "\n (" << filename << " " << numline << ")\n"
	      << std::endl;

    // Write warning also to logfile
    (*cla) << "\n WARNING:\n\n " << warnMsg
	   << "\n (" << filename << " " << numline << ")\n"
	   << std::endl;
  }


  // ***********************
  //   Warning (version 3)
  // ***********************
  void Warning( const Char*   const filename,
                const Integer       numline,
                const Char*   const formatstring, ... )
  {
    Integer  buffSize = 255,
             msgSize  = buffSize;
    char    *message  = NULL;

    // Since msgSize was initialized with buffSize,
    // this loop is entered at least once.
    while( buffSize <= msgSize ) {
        // create buffer for message
        delete [] message;
        message = New char [buffSize = msgSize+1];
        // try to format the message into the buffer
        std::va_list arg;
        va_start( arg, formatstring );
        msgSize = vsnprintf( message, buffSize, formatstring, arg );
        va_end( arg );
    }

    (*warning) << message;
    Warning( filename, numline );
    delete [] message;
  }


} // end of namespace
