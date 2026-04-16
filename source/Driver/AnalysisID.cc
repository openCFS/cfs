#include "Driver/AnalysisID.hh"
#include "DataInOut/ProgramOptions.hh"
#include "Utils/ToolsFull.hh"
#include "Domain/Domain.hh"
#include "Driver/EigenFrequencyDriver.hh"

namespace CoupledField
{

AnalysisID::AnalysisID()
{
  freq = -1.0;
  time = -1.0;
  iteration = -1;
  step = -1;
  coupleIter = -1;
  excite = "";
  adjoint = false;
}


std::string AnalysisID::ToString(bool filename) const
{
  std::string assign = filename ? "_" : ":";
  std::stringstream ss;

  if(step >= 0)
    ss << "step" << assign << step;

  assert(!(step < 0 && time >= 0.0));
  if(time >= 0.0 && !filename) // we already have step
    ss << (ss.rdbuf()->in_avail() ? "_" : "")  <<  "time" << assign << time;

  assert(!(step < 0 && freq >= 0.0));
  if(freq >= 0.0 && !filename) // we already have step
    ss << (ss.rdbuf()->in_avail() ? "_" : "")  <<  "freq" << assign << freq;

  if(coupleIter >= 0)
      ss << (ss.rdbuf()->in_avail() ? "_" : "")  << "coupleIter" << assign << coupleIter;

  if(pdeName != "")
     ss << (ss.rdbuf()->in_avail() ? "_" : "")  << assign << (filename ? ConvertToFilename(pdeName) : pdeName);

  if(iteration >= 0)
    ss << (ss.rdbuf()->in_avail() ? "_" : "")  << "iter" << assign << iteration;

  if(adjoint)
    ss << (ss.rdbuf()->in_avail() ? "_" : "")  << "adjoint";

  if(excite != "")
    ss << (ss.rdbuf()->in_avail() ? "_" : "")  << "excite" << assign << (filename ? ConvertToFilename(excite) : excite);

  // can be empty in the non-optimization case
  if(ss.str().size() == 0)
    if(domain->GetSingleDriver() != NULL && domain->GetSingleDriver()->DoBlochModeEigenfrequency())
      ss << "wv_(" << dynamic_cast<EigenFrequencyDriver*>(domain->GetSingleDriver())->GetCurrentWaveVector().ToString(TS_PLAIN,",") << ")";

  return ss.str();
}

ParamNode::ActionType AnalysisID::GetActionType() const
{
  return progOpts->DoDetailedInfo() ? ParamNode::APPEND : ParamNode::DEFAULT;
}


}
