// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_CFS_EXCEPTION_HH
#define FILE_CFS_EXCEPTION_HH

#include <string>
#include <sstream>
#include <exception>

namespace CoupledField {

 #define EXCEPTION(STR){                                          \
    std::ostringstream ostr;                                      \
    ostr << STR;                                                  \
    Exception ex__(NULL, __FILE__, __LINE__, ostr.str().c_str()); \
    throw ex__;                                                   \
  }
    
 #define RETHROW_EXCEPTION(REASON, STR){                               \
    std::ostringstream ostr;                                           \
    ostr << STR;                                                       \
    Exception ex__( & REASON, __FILE__, __LINE__, ostr.str().c_str()); \
    throw ex__;                                                        \
  }

  //! Base class for exception Handling
    class Exception : public std::exception {
    
  public:
    /** Creates an exception which can be thrown and then
     * catched by the main program. If the corresponding
     * flag is set this constructor will generate a segfault.
     * This is useful for debugging because the exception
     * and therefore the segfault will be generated at the
     * right place of the calling stack.
     */
    Exception( const Exception* reason,
               const char* const fileName = "NO_FILENAME", 
               const unsigned int lineNum = 0,
               const char* const message = "NO_MESSAGE") throw ();

    /** This is the only constructor to call an Exception by hand */
    Exception(const std::string& message,
              const char* const fileName = "NO_FILENAME", 
              const unsigned int lineNum = 0) throw ();  

    //! Copy constructor
    Exception( const Exception& exc ) throw ();

    //! Destructor
        virtual ~Exception() throw();

    //! Return message with file name and line number
    virtual const char * what() const throw ();

    //! Return exception message
    std::string GetMsg() const;
    
    //! Return file name where the exception occured
    std::string GetFileName() const;

    //! Return line number where the error occured
    unsigned int GetLineNum() const;

    //! Generate a segfault or not
    static bool segfault_;

  private:
    /** common init method */
    void init( const Exception* reason,
               const char* const fileName, 
               const unsigned int lineNum,
               const char* const message) throw ();

    //! Accumulated error message
    std::string what_;

    //! Error message
    std::string msg_;

    //! File name, where the exception occured
    std::string fileName_;

    //! Line number where the exeption occured
    unsigned int lineNum_;

    //! The exception which lead to the current exception
    Exception* reason_;

  };

}

#endif
