#include <numeric>

#include "BaseFE.hh"
#include "Domain/ElemMapping/Elem.hh"
#include "MatVec/Matrix.hh"


namespace CoupledField
{

  BaseFE::BaseFE() {
    actNumFncs_ = 0;
    feType_ = Elem::ET_UNDEF;
    isIsotropic_ = true;
    preComputShFnc_ = true;
  }

  BaseFE::~BaseFE() {
      
  }

  UInt BaseFE::GetNumFncs( EntityType entType,
                           UInt dof ) {
    StdVector<UInt> numFncs;
    this->GetNumFncs(numFncs, entType, dof);
   
    // Sum up all entries of vector
    return std::accumulate(numFncs.Begin(), numFncs.End(), 0);
  }
  // ************************************************************************
  // ENUM INITIALIZATION
  // ************************************************************************
      
  // Definition of entity types
  static EnumTuple entityTypeTuples[] = {
    EnumTuple(BaseFE::ALL,      "ALL"), 
    EnumTuple(BaseFE::VERTEX,   "VERTEX"),
    EnumTuple(BaseFE::NODE,     "NODE"),
    EnumTuple(BaseFE::EDGE,     "EDGE"),
    EnumTuple(BaseFE::FACE,     "FACE"),
    EnumTuple(BaseFE::INTERIOR, "INTERIOR")
  };
  Enum<BaseFE::EntityType> BaseFE::entityType = \
      Enum<BaseFE::EntityType>("Finite Element Entity Types",
          sizeof(entityTypeTuples) / sizeof(EnumTuple),
          entityTypeTuples);

} // end namespace CoupledField
