// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_FILEREADER_fvuns_2011
#define FILE_FILEREADER_fvuns_2011

#include <def_cplreader.hh>
#include "cplreader/FileReader.hh"
#include "FieldViewFile.h"

namespace CoupledField {

class FileReader_fvuns : public FileReader {

public:

  // constructor
  FileReader_fvuns(
      const std::string& name,
      const UInt dim,
      const UInt numFiles,
      const UInt startIndex);

  // destructor
  virtual ~FileReader_fvuns();

  // pure virtual functions of parent class
  void Init();
  void ReadNodalCoords(std::vector<Double> & NODECOORD);
  void ReadTopology(std::vector<UInt> & TOPOLOGYDATA, std::vector<UInt> & elemTypes);

  // not implemented in parent class
  void GetRegionElements(std::vector<UInt> & regionElements, const UInt regionIdx);
  void ReadNodalValues(std::vector<FlowDataType>& nodalFlowData,
                                   const std::vector<bool>& activeParts,
                                   const UInt timeStepIdx);


private:

  // fvuns file object
  FieldViewFile fvunsFile_;

  // the available solution types
  std::vector<SolutionType> solutionTypes_;
  // their DOF index (0 for scalar, 0,1,2 for vectors)
  std::vector<UInt> dofIndices_;

};

} // end of namespace

#endif
