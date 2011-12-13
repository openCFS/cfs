#ifndef CFS_RESULTCACHE_HH_
#define CFS_RESULTCACHE_HH_

#include <map>
#include <string>

#include "General/defs.hh"
#include "General/environment.hh"
#include "MatVec/basematrix.hh"
#include "MatVec/vector.hh"

namespace CoupledField
{

/*! Class for caching data of a result handler for use through the math parser
 *  function 'input'.
 *
 *  All member functions and variables have to be static, because the way math
 *  parser calls functions requires it.
 */
class ResultCache
{

public:

  //! Type defining the output value of the GetResult function for harmonic
  //! data: real part, imaginary part, amplitude, phase
  typedef enum { OUT_REAL, OUT_IMAG, OUT_AMPL, OUT_PHASE } OutputType;

  //! Backend for the 'input' function in math parser
  /*! This function is called by the math parser if the user gives an
   *  'input(...)' expression. Before evaluating the expression, one has to
   *  call the SetInfo function (or the other setter functions) in order to
   *  provide solution type, dof, region name, entity index, time step, and
   *  output type.
   *  \param inputId      ID of the InputReader
   *  \param sequenceStep number of the multisequence step
   */
  static Double GetResult(const char* inputId, Double sequenceStep, Double dof);

  //! Sets the desired output type used by GetResult.
  static void SetOutputType(OutputType outType);

  //! Sets the desired dof used by GetResult.
  static void SetDof(UInt dof);

  //! Sets the desired entity (region, etc.) used by GetResult.
  static void SetEntity(const std::string& entityName);

  //! Sets the index of the desired node/element used by GetResult.
  static inline void SetIndex(UInt index) {
    index_ = index;
  }

  //! Unifies all setter functions in one
  /*! Sets all the info required by GetResult in one function call. It does
   *  not include SetIndex, because this function is usually called quite
   *  often.
   *  \param outputType   desired type of output
   *  \param dof          desired dof
   *  \param entity       desired entity (region)
   *  \param solType      desired solution type
   */
  static void SetInfo(OutputType outType,
                      UInt dof,
                      const std::string& entityName,
                      SolutionType solType);

  //! Sets the desired solution type used by GetResult.
  static void SetSolution(SolutionType solType);

  //! Sets the desired step value used by GetResult.
  static void SetStepValue(Double stepValue);

private:

  //! Function to load the step values from the input file.
  static void LoadStepValues(const std::string& readerId, UInt sequenceStep);

  //! Flag that indicates, if the cached result is still valid.
  static bool resultValid_;

  //! Flag that indicates, if the step values are still valid for the current
  //! InputReader.
  static bool stepsValid_;

  //! current time/frequency step
  static Double curStepVal_;

  // previous time/frequency step
  //static Double prevStep_;

  //! requested data type
  static OutputType outType_;

  //! entry type of dataset
  static BaseMatrix::EntryType entryType_;

  //! requested DoF
  static UInt dof_;

  //! index of requested node/element
  static UInt index_;

  //! number of the most recently used step
  static UInt lastStepNum_;

  //! total number of dofs of the current solution
  static UInt numDofs_;

  //! most recently used multi-sequence step
  static UInt sequenceStep_;

  //! most recently used entity (region)
  static std::string entityName_;

  //! ID of most recently used InputReader
  static std::string readerId_;

  //! time/frequency step values of current result
  static std::map<UInt, Double> stepValues_;

  //! type of desired result
  static SolutionType solType_;

  //! cache memory for result of doubles
  static Vector<Double> resCacheDouble_;

  //! cache memory for result of complex numbers
  static Vector<Complex> resCacheComplex_;
  
  //! mapping from global node numbers to index in dataset
  static std::map<UInt, UInt> nodeNum2Index_;
};

} // namespace CoupledField

#endif /*CFS_RESULTCACHE_HH_*/
