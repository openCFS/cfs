#ifndef CPLREADER_FLOWDATATYPES
#define CPLREADER_FLOWDATATYPES

#include <General/environment.hh>
#include <Domain/resultInfo.hh>

namespace CoupledField
{
  typedef struct 
  {
    bool isActive;
    std::vector<Double> data;
    ResultInfo::EntityUnknownType definedOn;
  } FlowDataPartStruct;

  typedef std::map<SolutionType, FlowDataPartStruct> FlowDataType;

} // end of namespace

#endif // CPLREADER_FLOWDATATYPES
