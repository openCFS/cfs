// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <string>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <iomanip>
#include <sstream>

#include <boost/tokenizer.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
namespace fs=boost::filesystem;

#include "General/Exception.hh"
#include "cplreader/Settings.hh"
#include "FileReader_ANSYS.hh"


namespace CoupledField
{

  FileReader_ANSYS::FileReader_ANSYS(const std::string& name,
                                     const UInt dim,
                                     const UInt numFiles,
                                     const UInt startIndex) :
    FileReader(name, dim, numFiles, startIndex),
    strict_(false),
    degen_(false),
    everyThingRead_(false),
    anyThingRead_(false)
  {
    // Set chunk size to 10MB
    chunkSize_ = 10 * 1024 * 1024;
    lineIt_ = lineEnd_;
    preferredOutputPath_ = ".";
  }

  FileReader_ANSYS::~FileReader_ANSYS()
  {
    Settings& settings = Settings::Instance();
    std::vector<std::string>::const_iterator it, end;


    if(settings.GetInt("deltemp") && everyThingRead_)
    {
      it = fileNames_.begin();
      end = fileNames_.end();

      for( ; it != end; it++ )
      {
        if(settings.GetInt("verbose"))
          std::cout << "Removing " << (*it) << std::endl;

        fs::remove(*it);
      }

    }
  }

  std::string FileReader_ANSYS::GetReaderType()
  {
    Settings& settings = Settings::Instance();
    std::stringstream sstr;
    bool isMKHDF5 = false;
    bool isNACS = false;
    std::string basedir = settings.GetString("basedir");
    std::string name = settings.GetString("name");

    sstr << basedir << "/" << name << ".nodes";
    isMKHDF5  = fs::exists(sstr.str());
    sstr.clear(); sstr.str("");
    sstr << basedir << "/" << name << ".elements";
    isMKHDF5 &= fs::exists(sstr.str());

    if(isMKHDF5)
      return "MKHDF5";

    sstr.clear(); sstr.str("");
    sstr << basedir << "/" << name << "_Coordinates.dat";
    isNACS  = fs::exists(sstr.str());
    sstr.clear(); sstr.str("");
    sstr << basedir << "/" << name << "_eltype.dat";
    isNACS &= fs::exists(sstr.str());
    sstr.clear(); sstr.str("");
    sstr << basedir << "/" << name << "_elconn.dat";
    isNACS &= fs::exists(sstr.str());

    if(isNACS)
      return "NACS";

    return "FAILED";
  }

  void FileReader_ANSYS::ReadNodalCoords(std::vector<Double> & NODECOORD)
  {
    NODECOORD = nodalCoords_;
  }

  void FileReader_ANSYS::ReadTopology(std::vector<UInt> & TOPOLOGYDATA,
                                      std::vector<UInt> & elemTypes)
  {
    TOPOLOGYDATA.resize(maxNumElemNodes_ * elemNumsMap_.size());
    elemTypes.resize(elemNumsMap_.size());

    std::map<UInt, UInt>::const_iterator it, end;
    UInt numElemNodes;
    UInt elemNum;
    UInt baseIdx=0;
    UInt idx=0;
    UInt nodeNum;
    UInt origNodeNum;

    it = elemNumsMap_.begin();
    end = elemNumsMap_.end();

    for( ; it != end; it++ )
    {
      elemNum = it->first;
      Exception("Elem::GetNumElemNodes() no longer definied due to refactoring");
      //numElemNodes = Elem::GetNumElemNodes((Elem::FEType)elemTypes_[elemNum]);
      elemTypes[idx] = elemTypes_[elemNum];
      numElemNodes = 0;
      for(UInt i=0; i<numElemNodes; i++)
      {
        origNodeNum = topology_[elemNum][i];
        nodeNum = nodeNumsMap_[origNodeNum];

        if(!nodeNum)
          EXCEPTION("Node " << origNodeNum << " referenced by element "
                    << elemNum << " cannot by found in nodes array.");

        TOPOLOGYDATA[baseIdx+i] = nodeNum;
      }

      baseIdx += maxNumElemNodes_;
      idx++;
    }

  }

  void FileReader_ANSYS::GetRegionElements(std::vector<UInt> & regionElements,
                                           const UInt regionIdx)
  {
    UInt numRegionElems = elemNumsOrig_[regionIdx].size();
    UInt elemNum;

    regionElements.resize(numRegionElems);

    for(UInt i=0; i<numRegionElems; i++)
    {
      elemNum = elemNumsOrig_[regionIdx][i];
      regionElements[i] = elemNumsMap_[elemNum];
    }
  }

  void FileReader_ANSYS::ReadNodalValues(std::vector<FlowDataType>& nodalFlowData,
                                         const std::vector<bool>& activeParts,
                                         const UInt timeStepIdx)
  {
    std::cerr << "File reader for ANSYS output files is not supposed to "
      "handle nodal values!" << std::endl;
  }

  std::string FileReader_ANSYS::GetRegionName(const UInt regionIdx)
  {
    return regionNames_[regionIdx];
  }

  void FileReader_ANSYS::OpenFile(std::string fn)
  {
    inFile_.clear();
    inFile_.open(fn.c_str(), std::ios::binary);

    if (!inFile_) {
      EXCEPTION("Can't open " << fn);
    }

    // Let's remember which files we opened, so that we can delete them in the end.
    fileNames_.push_back(fn);

    // Determine size of file
    inFile_.seekg(0,std::ios::end);
    fSize_ = inFile_.tellg();

    // start from the beginning
    inFile_.seekg(0,std::ios::beg);
  }

  bool FileReader_ANSYS::ReadChunk()
  {
    UInt pos = inFile_.tellg();

    if(pos >= fSize_)
    {
      anyThingRead_ = false;
      return false;
    }
    
    char* buf;
    buf = new char[fSize_+1];
    std::fill(buf, &buf[fSize_+1], 0);
    inFile_.read(buf, fSize_);
    inFile_.close();
    buf[fSize_] = 0;

    std::string data(buf);
    delete[] buf;

    typedef boost::tokenizer<char_separator<char> > Tok;
    boost::char_separator<char> sep("\n\r");
    Tok t(data, sep);
    Tok::iterator it, end;
    it = t.begin();
    end = t.end();

    lines_.resize(std::distance(it, end));
    std::copy(it, end, lines_.begin());

    lineIt_ = lines_.begin();
    lineEnd_ = lines_.end();
    anyThingRead_ = true;

    return true;
  }

  bool FileReader_ANSYS::GetNextLine(std::string& line)
  {
    if(!anyThingRead_ || lineIt_ == lineEnd_)
    {
      if(!ReadChunk())
      {
        lines_.clear();
        return false;
      }
    }

    line = *lineIt_;
    lineIt_++;
    return true;
  }

  void FileReader_ANSYS::DegenerateElement(const Elem::FEType elemTypeIn,
                                           Elem::FEType& elemTypeOut,
                                           std::vector<UInt>& elemNodes)
  {
    static std::vector<UInt> newElemNodes;
    static std::map<Elem::FEType, std::vector<UInt> > idxMap;
    elemTypeOut = elemTypeIn;

    switch(elemTypeIn)
    {
    case Elem::ET_TRIA3:
      elemTypeOut = Elem::ET_QUAD4;
      if(idxMap[elemTypeIn].empty())
      {
        UInt elemIdxMap[] = {0, 1, 2, 2};
        std::copy(elemIdxMap, elemIdxMap+4, std::back_inserter(idxMap[elemTypeIn]));
      }
      break;

    case Elem::ET_TRIA6:
      elemTypeOut = Elem::ET_QUAD8;

      if(idxMap[elemTypeIn].empty())
      {
        UInt elemIdxMap[] = {0, 1, 2, 2, 3, 4, 2, 5};
        std::copy(elemIdxMap, elemIdxMap+8, std::back_inserter(idxMap[elemTypeIn]));
      }
      break;

    case Elem::ET_TET4:
      elemTypeOut = Elem::ET_HEXA8;

      if(idxMap[elemTypeIn].empty())
      {
        UInt elemIdxMap[] = {0, 1, 2, 2, 3, 3, 3, 3};
        std::copy(elemIdxMap, elemIdxMap+8,  std::back_inserter(idxMap[elemTypeIn]));
      }
      break;

    case Elem::ET_TET10:
      elemTypeOut = Elem::ET_HEXA20;

      if(idxMap[elemTypeIn].empty())
      {
        UInt elemIdxMap[] = {0, 1, 2, 2,
                             3, 3, 3, 3,
                             4, 5, 2, 6,
                             3, 3, 3, 3,
                             7, 8, 9, 9};
        std::copy(elemIdxMap, elemIdxMap+20,  std::back_inserter(idxMap[elemTypeIn]));
      }
      break;

    case Elem::ET_PYRA5:
      elemTypeOut = Elem::ET_HEXA8;

      if(idxMap[elemTypeIn].empty())
      {
        UInt elemIdxMap[] = {0, 1, 2, 3,
                             4, 4, 4, 4};
        std::copy(elemIdxMap, elemIdxMap+8,  std::back_inserter(idxMap[elemTypeIn]));
      }
      break;

    case Elem::ET_PYRA13:
      elemTypeOut = Elem::ET_HEXA20;

      if(idxMap[elemTypeIn].empty())
      {
        UInt elemIdxMap[] = {0, 1, 2, 3,
                             4, 4, 4, 4,
                             5, 6, 7, 8,
                             4, 4, 4, 4,
                             9, 10, 11, 12};
        std::copy(elemIdxMap, elemIdxMap+20,  std::back_inserter(idxMap[elemTypeIn]));
      }
      break;

    case Elem::ET_WEDGE6:
      elemTypeOut = Elem::ET_HEXA8;

      if(idxMap[elemTypeIn].empty())
      {
        UInt elemIdxMap[] = {0, 1, 2, 2,
                             3, 4, 5, 5};
        std::copy(elemIdxMap, elemIdxMap+8,  std::back_inserter(idxMap[elemTypeIn]));
      }
      break;

    case Elem::ET_WEDGE15:
      elemTypeOut = Elem::ET_HEXA20;

      if(idxMap[elemTypeIn].empty())
      {
        UInt elemIdxMap[] = {0, 1, 2, 2,
                             3, 4, 5, 5,
                             6, 7, 2, 8,
                             9, 10, 5, 11,
                             12, 13, 14, 14};
        std::copy(elemIdxMap, elemIdxMap+20,  std::back_inserter(idxMap[elemTypeIn]));
      }
      break;

    default:
      return;
    }

    Exception("Elem::GetNumElemNodes() no longer definied due to refactoring");
    //newElemNodes.resize(Elem::GetNumElemNodes((Elem::FEType)elemTypeOut));
    std::vector<UInt>::const_iterator it, end;
    it = idxMap[elemTypeIn].begin();
    end = idxMap[elemTypeIn].end();

    for(UInt i=0 ; it != end; it++, i++)
    {
      newElemNodes[i] = elemNodes[*it];
    }

    std::copy(newElemNodes.begin(), newElemNodes.end(), elemNodes.begin());
  }

  bool FileReader_ANSYS::CompareElements(const std::vector<UInt>& elemNodes1,
                                         const std::vector<UInt>& elemNodes2)
  {
    std::vector<UInt>::const_iterator it1, it2, end1, end2;
    it1 = elemNodes1.begin();
    end1 = elemNodes1.end();
    it2 = elemNodes2.begin();
    end2 = elemNodes2.end();

    for(; it1 != end1 && it2 != end2; it1++, it2++)
    {
      if(!(*it1))
        return true;

      if(*it1 != *it2)
        return false;
    }

    return true;
  }

  void FileReader_ANSYS::ResortNodes(std::vector<UInt>& elemNodes)
  {
     std::set<UInt> nodeSet;
     UInt pos;

     nodeSet.insert(elemNodes[0]);
     pos = 1;

     for(UInt i=1, n=elemNodes.size(); i<n; i++) {
       if(nodeSet.find(elemNodes[i]) != nodeSet.end()) {
         continue;
       }

       elemNodes[pos++] = elemNodes[i];
       nodeSet.insert(elemNodes[i]);
     }
  }

} // end of namespace
