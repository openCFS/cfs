// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

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

    // Names of DOFs
    std::vector<std::string> dofNames;

    // Unit
    std::string unit;
    
    // Vector for actually saving data
    std::vector<Double> data;
    
    // On what kind of entity the result is stored (nodes, elems...)
    ResultInfo::EntityUnknownType definedOn;

    // What kind of entry is result? Scalar, vector, tensor...
    ResultInfo::EntryType entryType;

    // If CFS has no SolutionType store external result name
    std::string resultName;
  } FlowDataPartStruct;

  // Data type for storing a FlowDataPartStruct per solution type
  typedef std::map<SolutionType, FlowDataPartStruct> FlowDataType;

} // end of namespace

#endif // CPLREADER_FLOWDATATYPES
