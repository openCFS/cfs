#ifndef ENVIRONMENT_HH
#define ENVIRONMENT_HH

/** This almost everywhere included shall have common includes:
 * - don't include too much, better remove if possible
 * - very much try to avoid including cfs headers to prevent cyclic inclusions
 * - typedefs (enum) themselves are in EnvironmentTypes.hh, here are the Enum<TYPE>s.
 *   This is necessary to prevent cyclic inclusions
 * - typedefs which need a header (e.g. shared_ptr) are to be defined here
 */

#define _USE_MATH_DEFINES
#include <cmath>

#include <typeinfo>
#include <iostream>
#include <vector>

// includes for the C99 standard datatypes (e.g. uint32_t, long double)
#include <boost/cstdint.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/scoped_ptr.hpp>
#include <cfloat>
#include <any>

#include "General/defs.hh"
#include "Enum.hh"
#include "EnvironmentTypes.hh"

//! \file Environment.hh
//! This file contains some global macro, class and enumeration data type
//! definitions for openCFS.
namespace CoupledField {

  // Import Boost's namespace
  //using namespace boost;
  using boost::shared_ptr;
  using boost::weak_ptr;
  using boost::scoped_ptr;
  using std::any_cast;
  using boost::lexical_cast;
  using boost::char_separator;
  using boost::dynamic_pointer_cast;


  // Type definition for shared_ptr<CoefFunction>
  class CoefFunction;
  typedef boost::shared_ptr<CoefFunction> PtrCoefFct;
  typedef boost::weak_ptr<CoefFunction> WeakPtrCoefFct;

  //! type of data
  struct Global {
    typedef enum {INTEGER, REAL, IMAG, COMPLEX} ComplexPart;
    static Enum<ComplexPart> complexPart;
  };


  /** it makes actually not really sense to name the Enum's with Enum.
   * Better would be to name them by the type but lowercase, e.g.
   * Enum<SolutionType> solutionType as this is an object
   */
  extern Enum<SolutionType> SolutionTypeEnum;
  extern Enum<MaterialType> MaterialTypeEnum;
  extern Enum<MaterialClass> MaterialClassEnum;
  extern Enum<MaterialTensorNotation> tensorNotation;
  extern Enum<ApproxCurveType> ApproxCurveTypeEnum;
  extern Enum<NonLinMethodType> NonLinMethodTypeEnum;

  extern Enum<FEMatrixType> feMatrixType;

  /** String2Enum/Enum2String is depreciated, better use Enum<> */
  
  //! conversion from strings to enum types
  template <class TYPE>
  void String2Enum(const std::string &in, TYPE &out);

  //! conversion from enum types to strings
  template<class TYPE>
  void Enum2String(const TYPE &in, std::string &out);


  // Instantiation for all known enum types;
#define DEFINE_ENUM_CONVERSION(TYPE)                                  \
  template<typename TYPE> void String2Enum(const std::string &in, TYPE &out); \
  template<typename TYPE> void Enum2String(const TYPE &in, std::string &out);

  DEFINE_ENUM_CONVERSION(FreqSamplingType)
  DEFINE_ENUM_CONVERSION(CouplingInputType)
  DEFINE_ENUM_CONVERSION(CouplingOutputType)
  DEFINE_ENUM_CONVERSION(CouplingRegionType)
  DEFINE_ENUM_CONVERSION(NormType)
  DEFINE_ENUM_CONVERSION(ComplexFormat)
  DEFINE_ENUM_CONVERSION(EQNType)
  DEFINE_ENUM_CONVERSION(MaterialClass)
  DEFINE_ENUM_CONVERSION(IntegMethod)
  DEFINE_ENUM_CONVERSION(NonLinType)
  DEFINE_ENUM_CONVERSION(TerminalConnector)
  DEFINE_ENUM_CONVERSION(DampingType)
  DEFINE_ENUM_CONVERSION(StopCritType)
  DEFINE_ENUM_CONVERSION(FEMatrixType)
  DEFINE_ENUM_CONVERSION(IDBCType)
  DEFINE_ENUM_CONVERSION(localInversionFlag)

#undef DEFINE_ENUM_CONVERSION

  std::string MapSolTypeToUnit(SolutionType solType);

  void SetEnvironmentEnums();

  /** Adds a generic solution after the initial initialization and returns the SolutioType
   * This can be used when dealing with generic post-processing results which might get defined at a later stage. */
  void AddGenericSolution(std::string name, Domain* domain);

  /** Function that returns the solution as string */
  std::string GetSolAsString(std::string name);

  /** sets the global CFS_NUM_THREADS variable.
   * CFS_NUM_THREADS controls e.g. parallel FEM assembly and it is important
   * to no call some cfs functions like coef functions with more threads than this!
   * OMP_NUM_THREADS and MKL_NUM_THREADS are environment variables. Setting them only effects this process, not the system.
   * @param homogenize also sets environment variables OMP_NUM_THREADS and MKL_NUM_THREADS used e.g. from libs
   * @param quiet suppress command line output when OMP and MKL_NUM_THREADS are set. */
  void SetNumberOfThreads(int numThreads, bool homogenize = false, bool quiet=true);

  /** Return the environment variable for the BLAS implementation.
   * Almost all Linux and Windows use Intel's MKL (BLAS and Pardiso) und return MKL_NUM_THREADS.
   * Apple's Accelerate Framework is used as system BLAS (not bundled) uses VECLIB_MAXIMUM_THREADS
   * OpenBLAS uses as many libs OMP_NUM_THREADS but we return DUMMY_NUM_THREADS to not double OMP_NUM_THREADS
   * NETLIB is serial so wer return OMP_NUM_THREADS
   * @param full if false return MKL, VECLIB, DUMMY */
  const char* GetBlasThreadsEnvVariable(bool full = true);

  /** true for mkl and accelerate (MKL_NUM_THREADS and VECLIB_MAXIMUM_THREADS),
   * false for openblas and netlib as we have here no environment variable to set (for simplicity DUMMY_NUM_THREADS is used) */
  bool HasBlasThreadsEnvVariable();


} // end of namespace

#endif // ENVIRONMENT_HH

