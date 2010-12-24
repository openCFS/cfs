#include <string>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <iomanip>
#include <sstream>

#include <boost/tokenizer.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
namespace fs=boost::filesystem;

#include "General/exception.hh"
#include "cplreader/Settings.hh"
#include "FileReader_CfsHdf5.hh"


namespace CoupledField
{

FileReader_CfsHdf5::FileReader_CfsHdf5(const std::string& name,
                                   const UInt dim,
                                   const UInt numFiles,
                                   const UInt startIndex) :
  FileReader(name, dim, numFiles),
  hdf5Reader_()
{
}


FileReader_CfsHdf5::~FileReader_CfsHdf5()
{
}


void FileReader_CfsHdf5::Init()
{
  std::map<unsigned int, unsigned int> numSteps;
  std::map<unsigned int, H5CFS::AnalysisType> analysis;

  try
  {
    hdf5Reader_.LoadFile(name_);
  } catch (std::string &strEx) {
    EXCEPTION(strEx);
  };

  dim_ = hdf5Reader_.GetDim();
  hdf5Reader_.GetRegionNamesOfDim(regionNames_, dim_);
  if (regionNames_.size() != 1)
  {
    EXCEPTION("only works for one region!");
  }

  std::vector<std::vector<unsigned int> > connectTmp;
  std::vector<H5CFS::ElemType> elemTypesTmp;
  hdf5Reader_.GetElems(elemTypesTmp, connectTmp);
  std::vector<unsigned int> elemsRegion;
  hdf5Reader_.GetElemsOfRegion(regionNames_[0], elemsRegion);
  numElemsPerRegion_.push_back(elemsRegion.size());
  hdf5Reader_.GetNumMultiSequenceSteps(analysis, numSteps);
  if (numSteps.size() != numSteps_)
  {
    std::cerr << "WARNING: more then one multistep, not implemented."
      << " Will take results from first multistep.";
  }
  if (numSteps[1] < numSteps_)
  {
    EXCEPTION("Not enough steps in files!");
  }

  std::vector<unsigned int>::const_iterator iterElemReg = elemsRegion.begin();
  if (iterElemReg != elemsRegion.end())
  {
    maxNumElemNodes_ = (UInt)H5CFS::NUM_ELEM_NODES[elemTypesTmp[*iterElemReg]];
  } else {
    EXCEPTION("no elements in region" << regionNames_[0]);
  }

  for (; iterElemReg != elemsRegion.end(); ++iterElemReg)
  {
    if (maxNumElemNodes_ != (UInt)H5CFS::NUM_ELEM_NODES[elemTypesTmp[*iterElemReg -1]])
    {
      EXCEPTION("Varying element types in region. Not implemented yet!");
    }
  }
}


UInt FileReader_CfsHdf5::GetNumRegions()
{
  return regionNames_.size();
#if 0 // DEBUG_SZOERNER
  return (UInt)hdf5Reader_.GetNumRegions();
#endif
}


void FileReader_CfsHdf5::ReadNodalCoords(std::vector<Double>& nodalCoords)
{
  std::vector<std::vector<double> > nodalCoordsTmp;
  hdf5Reader_.GetNodeCoords(nodalCoordsTmp);
  nodalCoords.resize(3 * nodalCoordsTmp.size());
  std::vector<std::vector<double> >::iterator iterNodalCoordsTmp = nodalCoordsTmp.begin();

  for (UInt i = 0; iterNodalCoordsTmp != nodalCoordsTmp.end(); ++iterNodalCoordsTmp)
  {
    for (UInt j = 0; j < 3; ++j)
    {
      nodalCoords[i] = (*iterNodalCoordsTmp)[j];
      ++i;
    }
  }
}


void FileReader_CfsHdf5::ReadTopology(std::vector<UInt>& connectivities, \
                                    std::vector<UInt>& elemTypes)
{
  std::vector<std::vector<unsigned int> > connectTmp;
  std::vector<H5CFS::ElemType> elemTypesTmp;
  std::vector<unsigned int> elemsRegion;
  hdf5Reader_.GetElems(elemTypesTmp, connectTmp);
  hdf5Reader_.GetElemsOfRegion(regionNames_[0], elemsRegion);

  std::vector<unsigned int>::const_iterator iterElemReg = elemsRegion.begin();
  std::vector<unsigned int>::const_iterator iterConnect_lvl2, iterConnect_lvl2End;
  for (UInt i = 0; iterElemReg != elemsRegion.end(); ++iterElemReg)
  {
    iterConnect_lvl2 = connectTmp[*iterElemReg -1].begin();
    iterConnect_lvl2End = connectTmp[*iterElemReg -1].end();
    for (UInt j = 0; iterConnect_lvl2 != iterConnect_lvl2End; ++iterConnect_lvl2, ++j)
    {
      connectivities.push_back(*iterConnect_lvl2);
    }
    elemTypes.push_back((UInt)elemTypesTmp[*iterElemReg -1]);
    ++i;
  }
}


std::string FileReader_CfsHdf5::GetRegionName(const UInt regionIdx)
{
  return regionNames_[regionIdx];
}


void FileReader_CfsHdf5::GetRegionElements(std::vector<UInt> & regionElements,
  const UInt regionIdx)
{
  std::vector<unsigned int> regionElementsTmp;
  hdf5Reader_.GetElemsOfRegion(regionNames_[regionIdx], regionElementsTmp);
  regionElements.resize(regionElementsTmp.size());

  for (UInt i = 0; i < regionElementsTmp.size(); ++i)
  {
    regionElements[i] = i +1;
  }
}


void FileReader_CfsHdf5::ReadNodalValues(std::vector<FlowDataType>& nodalFlowData, \
                           const std::vector<bool>& activeParts, \
                           const UInt timeStepIdx)
{
  std::vector<shared_ptr<H5CFS::ResultInfo> > infos;
  const unsigned int sequStep = 1;
  hdf5Reader_.GetResultTypes(sequStep, infos);
  std::vector<shared_ptr<H5CFS::ResultInfo> >::const_iterator iterInfo = infos.begin();

  FlowDataPartStruct* fdPtr;
  iterInfo = infos.begin();
  for (UInt actRegion = 0; actRegion <= numRegions_; ++actRegion)
  {
    for (; iterInfo != infos.end(); ++iterInfo)
    {
      if (requiredResults_[SolutionTypeEnum.Parse((*iterInfo)->name)])
      {
        FlowDataType& fd = nodalFlowData[actRegion];
        fdPtr = &fd[SolutionTypeEnum.Parse((*iterInfo)->name)];
        shared_ptr<H5CFS::Result> resultPtr(new H5CFS::Result);
        resultPtr->resultInfo = *iterInfo;

        fdPtr->isActive = activeParts[actRegion];
        fdPtr->definedOn = ResultInfo::NODE;
        fdPtr->entryType = (ResultInfo::EntryType)(*iterInfo)->entryType;
        fdPtr->dofNames = (*iterInfo)->dofNames;
        fdPtr->unit = (*iterInfo)->unit;
        fdPtr->resultName = (*iterInfo)->name;

        hdf5Reader_.GetResult(sequStep, timeStepIdx +1, resultPtr);
        fdPtr->data = resultPtr->realVals;
      }
    }
  }
}

} // end of namespace
