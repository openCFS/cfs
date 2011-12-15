#include <string>
#include <iostream>
#include <fstream>
#include "stdio.h"
#include <iomanip>
#include <sstream>

#include "boost/tokenizer.hpp"
#include "boost/algorithm/string/trim.hpp"
#include "boost/filesystem/operations.hpp"
#include "boost/filesystem/path.hpp"
#include "DataInOut/SimInOut/hdf5/hdf5io.hh"
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
  FileReader(name, dim, numFiles, startIndex),
  hdf5Reader_()
{
}


FileReader_CfsHdf5::~FileReader_CfsHdf5()
{
}


void FileReader_CfsHdf5::Init()
{
  Settings& settings = Settings::Instance();

  std::map<unsigned int, unsigned int> numSteps;
  std::map<unsigned int, H5CFS::AnalysisType> analysis;

  try
  {
    hdf5Reader_.LoadFile(name_);
  } catch (std::string &strEx) {
    EXCEPTION(strEx);
  };

  reduce_elementOrder_ = settings.GetInt("reduceElementOrder");
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
  // get time step values
  std::vector<shared_ptr<H5CFS::ResultInfo> > resInfos;
  const unsigned int sequenceStep = 1;
  hdf5Reader_.GetResultTypes(sequenceStep, resInfos);
  hdf5Reader_.GetStepValues(sequenceStep, resInfos[0], timeStepValues_);

  if (numSteps.size() != 1)
  {
    std::cerr << "WARNING: more then one multistep, not implemented."
      << " Will take results from first multistep.";
  }
  if ((numSteps[1] - startIndex_ +1) < numSteps_)
  {
    EXCEPTION("Not enough steps in files!" << std::endl << \
        "First time step in file is:  " << numSteps[1] - timeStepValues_.size() +1 \
        << std::endl << \
        "Last  time step in file is:  " << numSteps[1] << std::endl << \
        "maximum number of numsptes is:  " << timeStepValues_.size() << std::endl << \
        "Pleas use options --firststep and --numsteps");
  }
  if ( (startIndex_ - 1) < numSteps[1] - timeStepValues_.size())
  {
    EXCEPTION("Not enough steps in files!" << std::endl << \
        "First time step in file is:  " << numSteps[1] - timeStepValues_.size() +1 \
        << std::endl << \
        "Last  time step in file is:  " << numSteps[1] << std::endl << \
        "maximum number of numsptes is:  " << timeStepValues_.size() << std::endl << \
        "Pleas use options --firststep and --numsteps");
  }

  // get an check element types
  if (elemTypesTmp.size() != 0)
  {
    elemTypes_.resize(elemTypesTmp.size());
    // std::vector<std::vector<unsigned int> >::const_iterator iterConnect_lvl1 = connectTmp.begin();
    // TODO: Unused variable iterConnect_lvl1
    std::vector<unsigned int>::const_iterator iterConnect_lvl2, iterConnect_lvl2End;
    for (UInt i = 0; i < elemTypesTmp.size(); ++i)
    {
      if (reduce_elementOrder_)
      {
        elemTypes_[i] = (UInt)elemTypesTmp[i] -1;
      } else {
        elemTypes_[i] = (UInt)elemTypesTmp[i];
      }
    }
    maxNumElemNodes_ = (UInt)H5CFS::NUM_ELEM_NODES[elemTypes_[0]];
    for (UInt i = 0; i < elemTypesTmp.size(); ++i)
    {
      if (maxNumElemNodes_ < (UInt)H5CFS::NUM_ELEM_NODES[elemTypes_[i]])
      {
        maxNumElemNodes_ = (UInt)H5CFS::NUM_ELEM_NODES[elemTypes_[i]];
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
  UInt numNodeElem;
  std::vector<std::vector<unsigned int> > connectTmp;
  std::vector<H5CFS::ElemType> elemTypesTmp;
  hdf5Reader_.GetElems(elemTypesTmp, connectTmp);

  elemTypes = elemTypes_;
  connectivities.resize(connectTmp.size() * maxNumElemNodes_);
  std::fill(connectivities.begin(),connectivities.end(), 0);

  std::vector<std::vector<unsigned int> >::const_iterator iterConnect_lvl1 = connectTmp.begin();
  std::vector<unsigned int>::const_iterator iterConnect_lvl2, iterConnect_lvl2End;
  for (UInt i = 0, j = 0; iterConnect_lvl1 != connectTmp.end(); ++iterConnect_lvl1, ++i)
  {
    iterConnect_lvl2    = iterConnect_lvl1->begin();
    iterConnect_lvl2End = iterConnect_lvl1->end();
    // this is needed if reduce_elementOrder_ is set
    numNodeElem = (UInt)H5CFS::NUM_ELEM_NODES[elemTypes_[i]];
    for (UInt k = 0; k < numNodeElem; ++k, ++j, ++iterConnect_lvl2)
    {
      connectivities[j] = *iterConnect_lvl2;
    }
    while ( j % maxNumElemNodes_ != 0 )
    {
      ++j;
    }
  }
}

void FileReader_CfsHdf5::CorrectNumbering(std::vector<Double>& nodalCoords,
                                          std::vector<UInt>& connectivities, \
                                          std::vector<UInt>& elemTypes)
{
  // elemTypes is alread corrected in Init()
  std::set<UInt> nodesSet(connectivities.begin(), connectivities.end());
  const UInt numNodes = nodesSet.size() -1;
  std::set<UInt>::iterator nodeIt, nodeItEnd;
  nodeIt = nodesSet.begin();
  ++nodeIt; // zero is no node
  nodeItEnd = nodesSet.end();
  nodesMapAsVec_.resize(numNodes +1);
  nodesMapAsVec_[0] = 0;
  //create new numbering to get a continued numbering of the nodes
  for (UInt i = 1; i <= numNodes; ++i, ++nodeIt)
  {
    nodesMap_[*nodeIt] = i;
    nodesMapAsVec_[i] = *nodeIt;
  }

  std::vector<Double> nodalCoordsCorrect;
  nodalCoordsCorrect.resize(3 * numNodes);
  UInt nodeCnt_old;
  UInt nodeCnt = 0;
  nodeIt = nodesSet.begin();
  ++nodeIt; // do not use the zero
  for (; nodeIt != nodeItEnd; ++nodeIt)
  {
    nodeCnt_old = 3 * (*nodeIt -1);
    for (UInt i = 0; i < 3; ++i)
    {
      nodalCoordsCorrect[nodeCnt] = nodalCoords[nodeCnt_old];
      nodeCnt++;
      nodeCnt_old++;
    }
  }

  // correct the connectivities
  std::vector<UInt>::iterator iterConnect_lvl1 = connectivities.begin();
  for (; iterConnect_lvl1 != connectivities.end(); ++iterConnect_lvl1)
  {
    if (*iterConnect_lvl1 != 0)
    {
      *iterConnect_lvl1 = nodesMap_[*iterConnect_lvl1];
    }
  }

  nodalCoords = nodalCoordsCorrect;
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
  UInt timeStepIdxUpdate = startIndex_ + timeStepIdx;
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
        /* this needs to be done since in the hdf5 file not the entryType is
         * stored but its integer from MapEntryType() */
        fdPtr->entryType = H5IO::MapEntryType((Integer)(*iterInfo)->entryType);
        fdPtr->dofNames = (*iterInfo)->dofNames;
        fdPtr->unit = (*iterInfo)->unit;
        fdPtr->resultName = (*iterInfo)->name;

        try {
          hdf5Reader_.GetResult(sequStep, timeStepIdxUpdate, resultPtr);
        } catch (std::string &strEx) {
          EXCEPTION(strEx);
        }
        //TODO: This fails if we have mesh dimension of 2 but 3 components in the results
        fdPtr->data = resultPtr->realVals;
      }
    }
  }
}

void FileReader_CfsHdf5::ReduceOrderOfNodalValues( \
                           std::vector<FlowDataType>& nodalFlowData, \
                           const std::map<std::string, std::vector<UInt> >& regionNodes)
{
  std::vector<Double> correctedData;
  UInt oldNodePos, newNodePos;
  std::vector<unsigned int> origRegNodes;

  // walk through all regions
  std::map<std::string, std::vector<UInt> >::const_iterator \
    actRegNodessIter = regionNodes.begin();
  for (UInt actRegion = 0; actRegNodessIter != regionNodes.end(); \
      ++actRegNodessIter, ++actRegion)
  {
    const std::vector<UInt> actRegNodes = actRegNodessIter->second;
    FlowDataType& fd = nodalFlowData[actRegion];
    // walk thorugh all soultions
    std::map<SolutionType, FlowDataPartStruct >::iterator \
      fdIter = fd.begin();
    std::map<UInt, UInt>&  mapOrigRegToNewReg = origRegToNewReg_[actRegion];

    for (; fdIter != fd.end(); ++fdIter)
    {
      FlowDataPartStruct* fdPtr = &(fdIter->second);

      // these will be calculated later on
      if (fdPtr->isActive && fdPtr->resultName != "acouRhsLoad" 
          && fdPtr->resultName != "acouRhsLoadDensity"
          && fdPtr->resultName != "acouDivLighthillTensor")
      {
        if (!mapOrigRegToNewReg.size())
        {
          hdf5Reader_.GetNodesOfRegion(actRegNodessIter->first, origRegNodes);
          if (!nodesMapAsVec_.size())
          {
            EXCEPTION("Implementation error: nodesMapAsVec_ is not set, check if" << std::endl
                << "CorrectNumbering has been called earlier");
          }
          mapRegionalOrigToReduce(mapOrigRegToNewReg, nodesMapAsVec_, \
              origRegNodes, actRegNodessIter->second);
        }

        const UInt& numDofs = fdPtr->dofNames.size();
        std::vector<Double>& oldData = fdPtr->data;
        correctedData.resize(actRegNodes.size() * numDofs);

        // copy data from original grid to reduced grid
        for (UInt i = 0; i < actRegNodes.size(); ++i)
        {
          oldNodePos = mapOrigRegToNewReg[i] * numDofs;
          newNodePos = i * numDofs;

          for (UInt dofCnt = 0; dofCnt < numDofs; ++dofCnt)
          {
            correctedData[newNodePos] = oldData[oldNodePos];
            newNodePos++;
            oldNodePos++;
          }
        } // end copy 
      }
      fdPtr->data = correctedData;
    }
  }
}

double FileReader_CfsHdf5::GetTimeStep(UInt stepNumber)
{
  return timeStepValues_[stepNumber + startIndex_];
}

} // end of namespace
