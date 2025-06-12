#ifndef FILE_CFS_EXCEPTION_HH
#define FILE_CFS_EXCEPTION_HH

#include <boost/function.hpp>
#include <exception>
#include <mutex>
#include <sstream>
#include <string>

namespace CoupledField {

inline std::mutex &get_exception_mutex() {
  // Local static variable to ensure a single instance
  static std::mutex exception_mutex;
  return exception_mutex;
}

/** this creates an exception and add filename and line number */
#define EXCEPTION(STR)                                                         \
  {                                                                            \
    std::lock_guard<std::mutex> lock(CoupledField::get_exception_mutex());     \
    std::ostringstream ostr;                                                   \
    ostr << STR;                                                               \
    CoupledField::Exception ex__(NULL, __FILE__, __LINE__,                     \
                                 ostr.str().c_str());                          \
    throw ex__;                                                                \
  }

/** this creates an exception without filename and line number */
#define PLAIN_EXCEPTION(STR)                                                   \
  {                                                                            \
    std::lock_guard<std::mutex> lock(CoupledField::get_exception_mutex());     \
    std::ostringstream ostr;                                                   \
    ostr << STR;                                                               \
    CoupledField::Exception ex__(ostr.str().c_str());                          \
    throw ex__;                                                                \
  }

#define RETHROW_EXCEPTION(REASON, STR)                                         \
  {                                                                            \
    std::lock_guard<std::mutex> lock(CoupledField::get_exception_mutex());     \
    std::ostringstream ostr;                                                   \
    ostr << STR;                                                               \
    Exception ex__(&REASON, __FILE__, __LINE__, ostr.str().c_str());           \
    throw ex__;                                                                \
  }

//! Macro for issuing an warning message, will not terminate the program

//! This function can be used to issue a warning message. It is actually
//! only a shortcut for calling the WARN() method of the WriteInfo class.
//! \param Text     Text of the warning message
//! \param filename This is intended to contain the
//!                 name of the module/file in which the problem occured. The
//!                 __FILE__ macro should be inserted in the call. The
//!                 argument is optional.
//! \param numline  This is intended to contain the
//!                 number of the code line of the module/file in which the
//!                 problem occured. The __LINE__ macro should be inserted in
//!                 the call. The argument is optional.
#define WARN(STR)                                                              \
  {                                                                            \
    std::lock_guard<std::mutex> lock(CoupledField::get_exception_mutex());     \
    std::ostringstream ostr;                                                   \
    ostr << STR;                                                               \
    Exception ex__(NULL, __FILE__, __LINE__, ostr.str().c_str(),               \
                   Exception::WARNING);                                        \
  }

//! Base class for exception Handling
class Exception : public std::exception {
public:
  enum SeverityType { EXCEPTION, WARNING };

  /** Creates an exception which can be thrown and then
   * catched by the main program. If the corresponding
   * flag is set this constructor will generate a segfault.
   * This is useful for debugging because the exception
   * and therefore the segfault will be generated at the
   * right place of the calling stack.
   */
  Exception(const Exception *reason, const char *const fileName = "NO_FILENAME",
            const unsigned int lineNum = 0,
            const char *const message = "NO_MESSAGE",
            SeverityType severity = EXCEPTION) throw();

  /** This is the only constructor to call an Exception by hand */
  Exception(const std::string &message,
            const char *const fileName = "NO_FILENAME",
            const unsigned int lineNum = 0,
            SeverityType severity = EXCEPTION) throw();

  /** Rethrow type */
  Exception(const std::string &message, const Exception &reason) throw();

  //! Copy constructor
  Exception(const Exception &exc) throw();

  //! Destructor
  virtual ~Exception();

  //! Set callback method for exceptions
  static void SetCallbackEx(boost::function<void(Exception &x)> cb);

  //! Set callback method for warning
  static void SetCallbackWarn(boost::function<void(Exception &x)> cb);

  //! Return message with file name and line number
  virtual const char *what() const throw();

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
  void init(const Exception *reason, const char *const fileName,
            const unsigned int lineNum, const char *const message,
            SeverityType severity) throw();

  //! Accumulated error message
  std::string what_;

  //! Error message
  std::string msg_;

  //! File name, where the exception occured
  std::string fileName_;

  //! Line number where the exeption occured
  unsigned int lineNum_;

  //! callback function for exception level
  static boost::function<void(Exception &x)> exCallback_;

  //! callback functions for warning level
  static boost::function<void(Exception &x)> warnCallback_;

  //! The exception which lead to the current exception
  Exception *reason_;
};

} // namespace CoupledField

#endif
