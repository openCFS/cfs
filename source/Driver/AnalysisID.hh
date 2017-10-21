#ifndef DRIVER_ANALYISID_HH_
#define DRIVER_ANALYISID_HH_

#include <string>
#include "DataInOut/ParamHandling/ParamNode.hh"


namespace CoupledField
{


/** An anlysis ID describes the currently solved linear system by a unique ID
 It is mainly meant for the info.xml to connect outputs, e.g. for the solver with results.
 Additionally exportLinSystem can label the system uniquely.
 Most complex stuff comes from optimization. */
struct AnalysisID
{
  AnalysisID();

  /** gives a string representation of the current state.
   * @param filename make the string nice for a filename part (exportLinSys) */
  std::string ToString(bool filename = false) const;

  /** Depending on the command line detail flag we do a lot of info.xml output. Gives the appropriate action type. */
  ParamNode::ActionType GetActionType() const;

  /** current iteration for optimization if applies. Otherwise -1 */
  int iteration;

  /** current step if applies. Otherwise -1 */
  int step;

  /** current frequency if applies. Otherwise -1.0 */
  double freq;

  /** current time if applies. Otherwise -1.0 */
  double time;

  /** current optimization excitation label if applies. */
  std::string excite;

  /** adjoint solution in optimization if applies. Otherwise false. */
  bool adjoint;
};


} // end of namespac

#endif /* DRIVER_ANALYISID_HH_ */
