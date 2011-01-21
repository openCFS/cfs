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
  firststep_ = startIndex;
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
  hdf5Reader_.GetAllRegionNames(regionNames_);
  numRegions_ = regionNames_.size();

  std::vector<std::vector<unsigned int> > connectTmp;
  std::vector<H5CFS::ElemType> elemTypesTmp;
  hdf5Reader_.GetElems(elemTypesTmp, connectTmp);
  std::vector<unsigned int> elemsRegion;
  for (UInt i = 0; i < numRegions_; ++i)
  {
    hdf5Reader_.GetElemsOfRegion(regionNames_[i], elemsRegion);
    numElemsPerRegion_.push_back(elemsRegion.size());
  }

  // get and check time steps
  hdf5Reader_.GetNumMultiSequenceSteps(analysis, numSteps);
  if (numSteps.size() != 1)
  {
    std::cerr << "WARNING: more then one multistep, not implemented."
      << " Will take results from first multistep.";
  }
  if ((numSteps[1] - firststep_ +1) < numSteps_)
  {
    EXCEPTION("Not enough steps in files!" << std::endl << \
        "Last time step in file is:  " << numSteps[1] << std::endl << \
        "Please check first time step and set it with the option --firststep");
  }
  // get time step values
  std::vector<shared_ptr<H5CFS::ResultInfo> > resInfos;
  const unsigned int sequenceStep = 1;
  hdf5Reader_.GetResultTypes(sequenceStep, resInfos);
  hdf5Reader_.GetStepValues(sequenceStep, resInfos[0], timeStepValues_);

  // get an check element types
  if (elemTypesTmp.size() != 0)
  {
    maxNumElemNodes_ = (UInt)H5CFS::NUM_ELEM_NODES[elemTypesTmp[0]];
    for (UInt i = 0; i < elemTypesTmp.size(); ++i)
    {
      if (maxNumElemNodes_ < (UInt)H5CFS::NUM_ELEM_NODES[elemTypesTmp[i]])
      {
        maxNumElemNodes_ = (UInt)H5CFS::NUM_ELEM_NODES[elemTypesTmp[i]];
      }
    }
  } else {
    EXCEPTION("no elements found!");
  }
}


UInt FileReader_CfsHdf5::GetNumRegions()
{
  return regionNames_.size();
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
  hdf5Reader_.GetElems(elemTypesTmp, connectTmp);

  elemTypes.resize(elemTypesTmp.size());
  connectivities.resize(connectTmp.size() * maxNumElemNodes_);

  std::vector<std::vector<unsigned int> >::const_iterator iterConnect_lvl1 = connectTmp.begin();
  std::vector<unsigned int>::const_iterator iterConnect_lvl2, iterConnect_lvl2End;
  for (UInt i = 0, j = 0; iterConnect_lvl1 != connectTmp.end(); ++iterConnect_lvl1, ++i)
  {
    iterConnect_lvl2    = iterConnect_lvl1->begin();
    iterConnect_lvl2End = iterConnect_lvl1->end();
    for (; iterConnect_lvl2 != iterConnect_lvl2End; ++iterConnect_lvl2)
    {
      connectivities[j] = *iterConnect_lvl2;
      ++j;
    }
    while ( j % maxNumElemNodes_ != 0 )
    {
      connectivities[j] = 0;
      ++j;
    }
    elemTypes[i] = (UInt)elemTypesTmp[i];
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

  regionElements = (std::vector<UInt>)regionElementsTmp;
}


void FileReader_CfsHdf5::ReadNodalValues(std::vector<FlowDataType>& nodalFlowData, \
                           const std::vector<bool>& activeParts, \
                           const UInt timeStepIdx)
{
  UInt timeStepIdxUpdate = firststep_ + timeStepIdx;
  std::vector<shared_ptr<H5CFS::ResultInfo> > infos;
  const unsigned int sequStep = 1;
  hdf5Reader_.GetResultTypes(sequStep, infos);
  std::vector<shared_ptr<H5CFS::ResultInfo> >::const_iterator iterInfo = infos.begin();

  FlowDataPartStruct* fdPtr;
  for (UInt actRegion = 0; actRegion < numRegions_; ++actRegion)
  {
    iterInfo = infos.begin();
    for (; iterInfo != infos.end(); ++iterInfo)
    {
      if (requiredResults_[SolutionTypeEnum.Parse((*iterInfo)->name)]
          && (*iterInfo)->listName == regionNames_[actRegion])
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

        try {
          hdf5Reader_.GetResult(sequStep, timeStepIdxUpdate, resultPtr);
        } catch (std::string &strEx) {
          EXCEPTION(strEx);
        }
        fdPtr->data = resultPtr->realVals;
      }
    }
  }
}

double FileReader_CfsHdf5::GetTimeStep(UInt stepNumber)
{
  return timeStepValues_[stepNumber];
}

} // end of namespace
