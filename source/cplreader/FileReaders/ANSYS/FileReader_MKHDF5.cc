// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

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
#include "FileReader_MKHDF5.hh"


namespace CoupledField
{

  FileReader_MKHDF5::FileReader_MKHDF5(const std::string& name,
                                     const UInt dim,
                                     const UInt numFiles,
                                     const UInt startIndex) :
    FileReader_ANSYS(name, dim, numFiles, startIndex),
    maxOrigElemNum_(0)
  {
  }

  FileReader_MKHDF5::~FileReader_MKHDF5()
  {
  }

  void FileReader_MKHDF5::Init()
  {
    Settings& settings = Settings::Instance();
    std::stringstream sstr;
    std::string line;
    std::ostringstream errMsg;
    UInt maxNodeNum = 0;
    UInt numNodes = 0;
    std::map< std::string, std::vector<UInt> > entitySets;
    
    strict_ = settings.GetInt("strict") != 0;
    degen_ = settings.GetInt("degen") != 0;

    // Read .nodes file
    sstr.clear(); sstr.str("");
    sstr << settings.GetString("basedir") << "/" << name_ << ".nodes";

    ReadNWRITEFile(sstr.str(), entitySets);

    // Read .nodebc and .savenode files
    std::vector<std::string> nodeFiles;
    sstr.clear(); sstr.str("");
    sstr << settings.GetString("basedir") << "/" << name_ << ".nodebc";
    nodeFiles.push_back(sstr.str());
    sstr.clear(); sstr.str("");
    sstr << settings.GetString("basedir") << "/" << name_ << ".savenode";
    nodeFiles.push_back(sstr.str());

    for(UInt i=0, n=nodeFiles.size(); i < n; i++)
    {
      ReadNWRITEFile(nodeFiles[i], nodeGroups_);
    }

    // Read .elements file
    sstr.clear(); sstr.str("");
    sstr << settings.GetString("basedir") << "/" << name_ << ".elements";

    ReadEWRITEFile(sstr.str(), entitySets);

    // Read .saveelem file
    std::vector<std::string> elemFiles;
    sstr.clear(); sstr.str("");
    sstr << settings.GetString("basedir") << "/" << name_ << ".saveelem";
    elemFiles.push_back(sstr.str());

    for(UInt i=0, n=elemFiles.size(); i < n; i++)
    {
      ReadEWRITEFile(elemFiles[i], elemGroups_);
    }

    // Perform some checks for holes in numbering of nodes and elements
    std::map<UInt, UInt>::const_iterator it, end;
    UInt lastElemNum;
    UInt newElemNum = 1;
    
    it = elemTypes_.begin();
    end = elemTypes_.end();
    lastElemNum = it->first;
    elemNumsMap_[lastElemNum] = newElemNum;
    it++;
      
    // Check where exactly gaps are in the element array
    // and assign correct values to elemNumsMap_.
    for( ; it != end; it++ )
    {
      if(settings.GetInt("verbose") && (lastElemNum+1) != it->first)
        std::cerr << "Gap detected between elements " << lastElemNum
                  << " and " << it->first << std::endl;
        lastElemNum = it->first;
        elemNumsMap_[lastElemNum] = ++newElemNum;
    }
      
    UInt topoSize = topology_.size();
      
    if(maxOrigElemNum_ != topoSize)
    {
      errMsg.clear(); errMsg.str("");
      errMsg << "Element array written by ANSYS seems to have gaps!\n"
             << "Maybe issuing:\n"
             << "allsel,all,elem\n"
             << "numcmp,elem\n"
             << "before writing elements will help!";
      if(strict_)
      {
        EXCEPTION(errMsg.str());
      }
      else
        if(settings.GetInt("verbose"))
          std::cerr << errMsg.str() << std::endl;
    }
    
    everyThingRead_ = true;
    
  }

  void FileReader_MKHDF5::GetNodeGroups(std::map<std::string,
                                 std::vector<UInt> >& nodeGroups)
  {
    Settings& settings = Settings::Instance();
    std::vector<std::string> nodeFiles;
    std::string line;
    std::string groupName;
    std::stringstream sstr;
    std::set<UInt> nodeSet;

    nodeGroups = nodeGroups_;
    nodeGroups_.clear();    
  }

  void FileReader_MKHDF5::GetElemGroups(std::map<std::string,
                                        std::vector<UInt> >& elemGroups)
  {
    std::map<std::string, std::vector<UInt> >::const_iterator egIt, egEnd;
    std::vector<UInt>::const_iterator it, end;
    
    egIt = elemGroups_.begin();
    egEnd = elemGroups_.end();

    // Return new element number instead of ANSYS ones
    for( ; egIt != egEnd; egIt++ ) 
    {
      std::string elemGroupName = egIt->first;
      it = egIt->second.begin();
      end = egIt->second.end();

      if(!elemGroups[elemGroupName].empty()) 
      {
        EXCEPTION("Element group " << elemGroupName
                  << " has already been defined!");
      }

      for( ; it != end; it++ ) 
      {
        elemGroups[elemGroupName].push_back(elemNumsMap_[*it]);
      }
    }
  }

  FEType FileReader_MKHDF5::ANSYSTypeToFEType(UInt type, UInt numNodes,
                                             bool& readAnotherLine)
  {
    FEType ret = ET_UNDEF;
    readAnotherLine = false;

    switch(type)
    {
    case 2:
    case 100:
      if(numNodes == 2 || numNodes == 0)
        ret = ET_LINE2;
      break;

    case 3:
    case 101:
      if(numNodes == 3 || numNodes == 0)
        ret = ET_LINE3;
      break;

    case 4:
      if(numNodes == 3 || numNodes == 0)
        ret = ET_TRIA3;
      break;

    case 5:
      if(numNodes == 6 || numNodes == 0)
        ret = ET_TRIA6;
      break;

    case 6: // rectangle
      if(numNodes == 0)
        ret = ET_QUAD4;

      switch(numNodes)
      {
      case 3:
        ret = ET_TRIA3;
        break;

      case 4:
        ret = ET_QUAD4;
        break;
      }
      break;

    case 7: // quad. rectangle
      if(numNodes == 0)
        ret = ET_QUAD8;

      switch(numNodes)
      {
      case 6:
        ret = ET_TRIA6;
        break;

      case 8:
        ret = ET_QUAD8;
        break;
      }
      break;

    case 8:
      if(numNodes == 4 || numNodes == 0)
        ret = ET_TET4;
      break;

    case 9:
      if(numNodes == 10 || numNodes == 0)
        ret = ET_TET10;
      readAnotherLine = true;
      break;

    case 10: // hexa
      if(numNodes == 0)
        ret = ET_HEXA8;

      switch(numNodes)
      {
      case 4:
        ret = ET_TET4;
        break;

      case 5:
        ret = ET_PYRA5;
        break;

      case 6:
        ret = ET_WEDGE6;
        break;

      case 8:
        ret = ET_HEXA8;
        break;
      }
      break;

    case 11: // quad. hexa
      if(numNodes == 0)
        ret = ET_HEXA20;

      switch(numNodes)
      {
      case 10:
      ret = ET_TET10;
      break;

      case 13:
      ret = ET_PYRA13;
      break;

      case 15:
      ret = ET_WEDGE15;
      break;

      case 20:
      ret = ET_HEXA20;
      break;
      }

      readAnotherLine = true;
      break;

    case 12:
      if(numNodes == 5 || numNodes == 0)
        ret = ET_PYRA5;
      break;

    case 13:
      if(numNodes == 13 || numNodes == 0)
        ret = ET_PYRA13;
      readAnotherLine = true;
      break;

    case 14:
      if(numNodes == 6 || numNodes == 0)
        ret = ET_WEDGE6;
      break;

    case 15:
      if(numNodes == 15 || numNodes == 0)
        ret = ET_WEDGE15;
      readAnotherLine = true;
      break;
    }

    return ret;
  }

  void FileReader_MKHDF5::ReadEWRITEFile(std::string fileName,
                                         std::map< std::string, std::vector<UInt> >& elemSets)
  {
    Settings& settings = Settings::Instance();

    std::string sectionName;
    std::ostringstream errMsg;
    std::stringstream sstr;
    UInt regionIdx = 0;
    bool readAnotherLine = false;
    FEType prelimElemType, actualElemType;
    UInt ansysElemType;
    UInt mat, real, secnum, esys;
    UInt elemNum;
    UInt numRegionElems = 0;
    UInt numTotalElems = 0;
    UInt numElemNodes;
    UInt elemDim;
    std::vector<UInt> elemNodes(40);
    std::set<UInt> elemNodeSet;
    UInt dim = 0;
    bool isElementsFile;
    UInt newElemNum = elemNumsMap_.size() + 1;
    std::set<UInt> elemSet;
    std::string line;
    bool firstSection = true;

    if(settings.GetInt("verbose"))
      std::cout << "Trying to read EWRITE file: " << fileName << std::endl;
    
    isElementsFile = fileName.find( ".elements") != std::string::npos;

    try
    {
      OpenFile(fileName);
      if(settings.GetInt("verbose"))
        std::cout << "OKAY!" << std::endl;
    } catch (std::exception& ex)
    {
      if(settings.GetInt("verbose"))
        std::cout << "FAILED!" << std::endl
                  << ex.what() << std::endl;
      return;
    }

    while(GetNextLine(line)) {
      std::fill(elemNodes.begin(), elemNodes.end(), 0);
      elemNodeSet.clear();

      // strip whitespaces
      boost::trim(line);

      if(line.length() == 0)
        continue;

      if(line[0] == '[')
      {
        if(!isElementsFile)
        {
          if(!elemSet.empty())
          {
            if(!elemSets[sectionName].empty()) 
            {
              EXCEPTION("Element group '" << sectionName << "' from file '"
                        << fileName << "' has already been defined!");
            }
            
            std::copy(elemSet.begin(), elemSet.end(),
                      std::back_inserter(elemSets[sectionName]));

            elemSet.clear();
          } else 
          {
            if(settings.GetInt("strict") && !firstSection)
            {
              EXCEPTION("Element group '" << sectionName << "' from file '"
                        << fileName << "' is empty!");
            }
          }
        }

        sectionName = line.substr(1, line.length()-2);
        boost::trim( sectionName );
        firstSection = false;
        
        //        std::cout << sectionName << std::endl;

        if(isElementsFile) 
        {  
          if(!regionNames_.empty() && !numRegionElems)
          {
            std::ostringstream sstr;
            
            sstr << "Region '" << (*regionNames_.rbegin()) << "' from file '"
                 << fileName << "' contains no elements!\n"
                 << "Is this what you intended to do?\n\n";
            
            if(settings.GetInt("strict")) 
            {
              EXCEPTION(sstr.str());
            }
            else
            {
              std::cerr << sstr.str();
              regionNames_.erase(regionNames_.begin()+regionNames_.size()-1);
              numRegions_--;
              numNodesPerRegion_.erase(numNodesPerRegion_.begin()+numNodesPerRegion_.size()-1);
              numElemsPerRegion_.erase(numElemsPerRegion_.begin()+numElemsPerRegion_.size()-1);
            }            
          }

          regionIdx = regionNames_.size();
          regionNames_.push_back( sectionName );

          numRegions_++;
          numNodesPerRegion_.push_back(6);
          if(regionIdx)
            numElemsPerRegion_.push_back(numRegionElems);
          numRegionElems = 0;
        }

        
        continue;
      }

      sstr.clear(); sstr.str("");
      sstr << line;
      sstr >> elemNodes[0] >> elemNodes[1] >> elemNodes[2] >> elemNodes[3]
           >> elemNodes[4] >> elemNodes[5] >> elemNodes[6] >> elemNodes[7];
      sstr >> mat >> ansysElemType >> real >> secnum >> esys >> elemNum;
      // Im Fall von NACS files bekommen wir den Elementtyp aus einer extra Datei

      prelimElemType = ANSYSTypeToFEType(ansysElemType, 0, readAnotherLine);

      if(readAnotherLine)
      {
        if(!GetNextLine(line))
          EXCEPTION("Error while trying to read second half of elem record.")

        sstr.clear(); sstr.str("");
        sstr << line;

        sstr >> elemNodes[8] >> elemNodes[9] >> elemNodes[10] >> elemNodes[11]
             >> elemNodes[12] >> elemNodes[13] >> elemNodes[14] >> elemNodes[15]
             >> elemNodes[16] >> elemNodes[17] >> elemNodes[18] >> elemNodes[19];
      }


      // Determine number of element nodes by inserting nodes into a set.
      elemNodeSet.insert(elemNodes.begin(), elemNodes.end());
      elemNodeSet.erase((UInt) 0);
      //      std::copy(elemNodeSet.begin(), elemNodeSet.end(),
      //                std::ostream_iterator< UInt >(std::cout, " "));
      //      std::cout << std::endl;
      numElemNodes = elemNodeSet.size();

      actualElemType = ANSYSTypeToFEType(ansysElemType,
                                         numElemNodes,
                                         readAnotherLine);
      if(actualElemType == ET_UNDEF)
        EXCEPTION("Found undefined element type for elem " << elemNum
                  << " (region: " << (*regionNames_.rbegin())
                  << ", type: " << ELEM_TYPE_NAMES[prelimElemType]
                  << ", num. nodes: " << numElemNodes
                  << ")");

      if(!degen_ && prelimElemType != actualElemType)
        ResortNodes(elemNodes);

      if(degen_)
      {
        DegenerateElement(prelimElemType, actualElemType, elemNodes);
        numElemNodes = NUM_ELEM_NODES[actualElemType];
      }

      // Dieser Fehler ist fatal! Er hat nichts mit strict oder relaxed zu tun.
      if(isElementsFile) 
      {  
        if(elemTypes_[elemNum] != 0)
        {
          // Check if element has really been referenced before
          if( CompareElements(elemNodes, topology_[elemNum]) )
            EXCEPTION("Element " << elemNum << " has been referenced before in region '"
                      << regionNames_[elemRegionMap_[elemNum]] << "' !\n"
                      << "Active region is '" << (*regionNames_.rbegin()) << "'.");
          
          // Add a really big number to elemNum and try to find a unique element number.
          elemNum += 2000000000;
          while(elemTypes_.find(elemNum) != elemTypes_.end())
            elemNum++;
        }
      }
      
      if(elemTypes_[elemNum] == 0)
      {
        std::copy(elemNodes.begin(),
                  elemNodes.begin()+numElemNodes,
                  std::back_inserter(topology_[elemNum]));
        elemTypes_[elemNum] = actualElemType;

        if(isElementsFile) 
        {
          // This is just for debugging output
          elemRegionMap_[elemNum] = regionIdx;

          // Noch irgendwie die Regionen ordnen
          // - elemente in ein konsekutives array einfügen
          // - Elementnummern pro Region speichern
          // - MaxNumElemNodes von Filereadern herausfinden lassen
          // - Elementconnectivity wird immer in arrays mit n*MaxNumElemNodes behandelt
          // - Elementnummern pro Region dem MpCCIExch object zugänglich machen.
          elemNumsOrig_[regionIdx].push_back(elemNum);
          numRegionElems++;
        }
        else
        {
          elemSet.insert(elemNum);
        }
        
        
        maxOrigElemNum_ = maxOrigElemNum_ < elemNum ? elemNum : maxOrigElemNum_;
        maxNumElemNodes_ = maxNumElemNodes_ < numElemNodes ? numElemNodes : maxNumElemNodes_;
        elemDim = ELEM_DIM[actualElemType];
        dim = dim < elemDim ? elemDim : dim;
      }
      else 
      {
        elemSet.insert(elemNum);
      }
    }
    
    UInt settingsDim = (UInt) settings.GetInt("dim");
    dim = dim < (UInt) settingsDim ? settingsDim : dim;
    settings.SetInt("dim", dim);

    if(isElementsFile) 
    {
      numElemsPerRegion_.push_back(numRegionElems);
    }
    else 
    {
      if(!elemSet.empty())
      {
        if(!elemSets[sectionName].empty()) 
        {
          EXCEPTION("Element group '" << sectionName << "' from file '"
                    << fileName << "' has already been defined!");
        }
        
        std::copy(elemSet.begin(), elemSet.end(),
                  std::back_inserter(elemSets[sectionName]));
      } else
      {
        if(settings.GetInt("strict") && !firstSection)
        {
          EXCEPTION("Element group '" << sectionName << "' from file '"
                    << fileName << "' is empty!");
        }
      }
    }
    

  }
  
  void FileReader_MKHDF5::ReadNWRITEFile(std::string fileName,
                                         std::map< std::string, std::vector<UInt> >& nodeSets)
  {
    Settings& settings = Settings::Instance();

    double x, y, z;
    std::string line;
    std::ostringstream errMsg;
    std::stringstream sstr;
    UInt nodeNum;
    UInt maxNodeNum = 0;
    UInt numNodes;
    std::string sectionName;
    std::set<UInt> nodeSet;
    bool isNodesFile;
    bool firstSection = true;
    
    if(settings.GetInt("verbose"))
      std::cout << "Reading NWRITE file: " << fileName << std::endl;

    isNodesFile = fileName.find(".nodes") != std::string::npos;
  
    try
    {
      OpenFile(fileName);
      if(settings.GetInt("verbose"))
        std::cout << "OKAY!" << std::endl;
    } catch (std::exception& ex)
    {
      if(settings.GetInt("verbose"))
        std::cout << "FAILED!" << std::endl
                  << ex.what() << std::endl;
      return;
    }
    
    sectionName = "";
    numNodes = nodeNumsMap_.size();
    
    while(GetNextLine(line)) {
      x = y = z = 0;

      // strip whitespaces
      boost::trim(line);
      
      if(line.length() == 0)
        continue;
      
      if(line[0] == '[')
      {
        if(!nodeSet.empty())
        {
          if(!nodeSets[sectionName].empty()) 
          {
            EXCEPTION("Node group '" << sectionName << "' from file '"
                      << fileName << "' has already been defined!");
          }
          
          std::copy(nodeSet.begin(), nodeSet.end(),
                    std::back_inserter(nodeSets[sectionName]));
          
          nodeSet.clear();
        }
        else
        {
          if(settings.GetInt("strict") && !firstSection)
          {
            EXCEPTION("Node group '" << sectionName << "' from file '"
                      << fileName << "' is empty!");
          }
        }
        
        sectionName = line.substr(1, line.length()-2);
        boost::trim( sectionName );
        firstSection = false;

        //          std::cout << sectionName << std::endl;
        continue;
      }
      
      sstr.clear(); sstr.str("");
      sstr << line;
      
      sstr >> nodeNum >> x >> y >> z;
      
      // Hier abfragen, ob die Knotennummern konsekutiv sind
      if(isNodesFile &&
         !nodeNumsMap_.empty() &&
         (nodeNumsMap_.rbegin()->first + 1) != nodeNum)
      {
        errMsg << "Nodes are not consecutive: " << nodeNumsMap_.rbegin()->first
               << " -> " << nodeNum;
        if(strict_)
        {
          EXCEPTION(errMsg.str());
        }
        else
        {
          if(settings.GetInt("verbose"))
            std::cerr << errMsg.str() << std::endl;
        }
        errMsg.clear(); errMsg.str("");
      }
      
      // Dieser Fehler ist fatal! Er hat nichts mit strict oder relaxed zu tun.
      if(isNodesFile) 
      {
        if(nodeNumsMap_[nodeNum] != 0)
          EXCEPTION("Node " << nodeNum << " has been referenced before!\n");
      }
    
      if(nodeNumsMap_[nodeNum] == 0) 
      {  
        nodalCoords_.push_back(x);
        nodalCoords_.push_back(y);
        nodalCoords_.push_back(z);
        
        maxNodeNum = maxNodeNum < nodeNum ? nodeNum : maxNodeNum;
        numNodes++;
        nodeNumsMap_[nodeNum] = numNodes;
      }
      
      nodeSet.insert(nodeNumsMap_[nodeNum]);
      
    }

    if(isNodesFile) 
    {
      errMsg << "Node array written by ANSYS seems to have gaps!";
      if(maxNodeNum != numNodes && strict_)
      {
        EXCEPTION(errMsg.str());
      }
      else
      {
        if(settings.GetInt("verbose"))
          std::cerr << errMsg.str() << std::endl;
      }
    }

    if(!nodeSet.empty())
    {
      if(!nodeSets[sectionName].empty()) 
      {
        EXCEPTION("Node group " << sectionName << " from file "
                  << fileName << " has already been defined!");
      }
      
      std::copy(nodeSet.begin(), nodeSet.end(),
                std::back_inserter(nodeSets[sectionName]));
      
      nodeSet.clear();
    } else 
    {    
      if(settings.GetInt("strict") && !firstSection)
      {
        EXCEPTION("Node group '" << sectionName << "' from file '"
                  << fileName << "' is empty!");
      }
    }
  }
  
} // end of namespace
