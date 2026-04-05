#include "BaseMatrix.hh"
#include "General/Enum.hh"
#include "General/Environment.hh"

namespace CoupledField {

std::string BaseMatrix::ToInfoString() const
{
  return "st=" + MatrixStructureTypeEnum.ToString(GetStructureType()) + " et=" + MatrixEntryTypeEnum.ToString(GetEntryType()) + " mu=" + std::to_string(GetMemoryUsage());
}




} // end of namespace
