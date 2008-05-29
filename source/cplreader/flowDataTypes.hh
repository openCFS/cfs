#ifndef CPLREADER_FLOWDATATYPES
#define CPLREADER_FLOWDATATYPES

#include <General/environment.hh>
#include <Domain/resultInfo.hh>

namespace CoupledField
{
  typedef struct 
  {
    // Is data actually active on the current partition?
    // E.g. we do not need to know the field on some boundaries.
    // This saves disc space!
    bool isActive;
    
    // Vector for actually saving data
    std::vector<Double> data;
    
    // On what kind of entity the result is stored (nodes, elems...)
    ResultInfo::EntityUnknownType definedOn;

    // If CFS has no SolutionType store external result name
    std::string resultName;
  } FlowDataPartStruct;

  // Data type for storing a FlowDataPartStruct per solution type
  typedef std::map<SolutionType, FlowDataPartStruct> FlowDataType;

} // end of namespace

#endif // CPLREADER_FLOWDATATYPES
