
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

#include "../settings.hh"
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
    double x, y, z;
    UInt nodeNum;
    UInt maxNodeNum = 0;
    UInt numNodes = 0;

    strict_ = (bool) settings.GetInt("strict");
    degen_ = (bool) settings.GetInt("degen");
    
    sstr.clear(); sstr.str("");
    sstr << settings.GetString("basedir") << "/" << name_ << ".nodes";

    OpenFile(sstr.str());
    
    while(GetNextLine(line)) {
      x = y = z = 0;
      
      // strip whitespaces
      boost::trim(line);

      if(line == "")
        continue;
      
      sstr.clear(); sstr.str("");
      sstr << line;
      
      sstr >> nodeNum >> x >> y >> z;
      
      // Hier abfragen, ob die Knotennummern konsekutiv sind
      if(!nodeNumsMap_.empty() && 
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
      if(nodeNumsMap_[nodeNum] != 0)
        EXCEPTION("Node " << nodeNum << " has been referenced before!\n");

      nodalCoords_.push_back(x);
      nodalCoords_.push_back(y);
      nodalCoords_.push_back(z);
      
      maxNodeNum = maxNodeNum < nodeNum ? nodeNum : maxNodeNum;
      numNodes++;
      nodeNumsMap_[nodeNum] = numNodes;
    }

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
    
    std::string regionName;
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
    
    sstr.clear(); sstr.str("");
    sstr << settings.GetString("basedir") << "/" << name_ << ".elements";

    OpenFile(sstr.str());
    
    while(GetNextLine(line)) {
      std::fill(elemNodes.begin(), elemNodes.end(), 0);
      elemNodeSet.clear();
      
      // strip whitespaces
      boost::trim(line);

      if(line.length() == 0)
        continue;

      if(line[0] == '[')
      {
        regionName = line.substr(1, line.length()-2);
        boost::trim( regionName );
        regionIdx = regionNames_.size();
        regionNames_.push_back( regionName );
        //        std::cout << regionName << std::endl;
        numRegions_++;
        numNodesPerRegion_.push_back(numNodes);
        if(regionIdx)
          numElemsPerRegion_.push_back(numRegionElems);
        numRegionElems = 0;

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
      
      std::copy(elemNodes.begin(),
                elemNodes.begin()+numElemNodes,
                std::back_inserter(topology_[elemNum]));
      elemTypes_[elemNum] = actualElemType;
      // This is just for debugging output
      elemRegionMap_[elemNum] = regionIdx;

      // Noch irgendwie die Regionen ordnen
      // - elemente in ein konsekutives array einfügen
      // - Elementnummern pro Region speichern
      // - MaxNumElemNodes von Filereadern herausfinden lassen
      // - Elementconnectivity wird immer in arrays mit n*MaxNumElemNodes behandelt
      // - Elementnummern pro Region dem MpCCIExch object zugänglich machen.
      elemNumsOrig_[regionIdx].push_back(elemNum);

      maxOrigElemNum_ = maxOrigElemNum_ < elemNum ? elemNum : maxOrigElemNum_;
      maxNumElemNodes_ = maxNumElemNodes_ < numElemNodes ? numElemNodes : maxNumElemNodes_;
      elemDim = ELEM_DIM[actualElemType];
      dim = dim < elemDim ? elemDim : dim;

      numRegionElems++;
      numTotalElems++;
    }
    numElemsPerRegion_.push_back(numRegionElems);
    settings.SetInt("dim", dim);

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
             << "before writing elemens will help!";
      if(strict_)
      {
        EXCEPTION(errMsg.str());
      }
      else
        if(settings.GetInt("verbose"))
          std::cerr << errMsg.str() << std::endl;
    }
    
  }

  void FileReader_MKHDF5::GetNodeGroups(std::map<std::string,
                                 std::vector<UInt> >& nodeGroups)
  {
    Settings& settings = Settings::Instance();
    std::vector<std::string> nodeFiles;
    std::string line;
    std::string groupName;
    UInt nodeNum, newNodeNum;
    Double x,y,z;
    std::stringstream sstr;
    std::set<UInt> nodeSet;
    bool readOK;
    
    sstr.clear(); sstr.str("");
    sstr << settings.GetString("basedir") << "/" << name_ << ".nodebc";
    nodeFiles.push_back(sstr.str());
    sstr.clear(); sstr.str("");
    sstr << settings.GetString("basedir") << "/" << name_ << ".savenode";
    nodeFiles.push_back(sstr.str());

    for(UInt i=0, n=nodeFiles.size(); i < n; i++)
    {
      try 
      {
        OpenFile(nodeFiles[i]);
        readOK = true;
      } catch (std::exception& ex) 
      {
        readOK = false;
        continue;
      }
      
      
      while(GetNextLine(line)) {
        // strip whitespaces
        boost::trim(line);

        if(line.length() == 0)
          continue;

        if(line[0] == '[')
        {
          if(!nodeSet.empty())
          {
            std::copy(nodeSet.begin(), nodeSet.end(),
                      std::back_inserter(nodeGroups[groupName]));

            nodeSet.clear();
          }
        
          groupName = line.substr(1, line.length()-2);
          boost::trim( groupName );
          //          std::cout << groupName << std::endl;
          continue;
        }
      
        sstr.clear(); sstr.str("");
        sstr << line;
        sstr >> nodeNum >> x >> y >> z;

        newNodeNum = nodeNumsMap_[nodeNum];

        if(!newNodeNum)
          EXCEPTION("Node number " << nodeNum << " from node group '"
                    << groupName << "' has not been defined before!");

        nodeSet.insert(newNodeNum);
      }

      std::copy(nodeSet.begin(), nodeSet.end(),
                std::back_inserter(nodeGroups[groupName]));
    }
  }

  void FileReader_MKHDF5::GetElemGroups(std::map<std::string,
                                 std::vector<UInt> >& elemGroups)
  {
    Settings& settings = Settings::Instance();

    std::vector<std::string> elemFiles;
    std::string line;
    std::string groupName;
    UInt elemNum, newElemNum;
    UInt dummy, elemType;
    std::stringstream sstr;
    std::set<UInt> elemSet;
    FEType prelimElemType;
    bool readAnotherLine;
    bool readOK;
    
    sstr.clear(); sstr.str("");
    sstr << settings.GetString("basedir") << "/" << name_ << ".saveelem";
    elemFiles.push_back(sstr.str());
    
    for(UInt i=0, n=elemFiles.size(); i < n; i++)
    {
      try 
      {
        OpenFile(elemFiles[i]);
        readOK = true;
      } catch (std::exception& ex) 
      {
        readOK = false;
        continue;
      }
    
      while(GetNextLine(line)) {
        // strip whitespaces
        boost::trim(line);

        if(line.length() == 0)
          continue;

        if(line[0] == '[')
        {
          if(!elemSet.empty())
          {
            std::copy(elemSet.begin(), elemSet.end(),
                      std::back_inserter(elemGroups[groupName]));

            elemSet.clear();
          }
        
          groupName = line.substr(1, line.length()-2);
          boost::trim( groupName );
          //          std::cout << groupName << std::endl;
          continue;
        }
      
        sstr.clear(); sstr.str("");
        sstr << line;
        sstr >> dummy >> dummy >> dummy >> dummy 
             >> dummy >> dummy >> dummy >> dummy;
        sstr >> dummy >> elemType >> dummy >> dummy >> dummy >> elemNum;
        // Im Fall von NACS files bekommen wir den Elementtyp aus einer extra Datei

        prelimElemType = ANSYSTypeToFEType(elemType, 0, readAnotherLine);
      
        if(readAnotherLine)
        {
          if(!GetNextLine(line))
            EXCEPTION("Error while trying to read second half of elem record.");
        }
        

        newElemNum = elemNumsMap_[elemNum];

        if(!newElemNum)
          EXCEPTION("Element number " << elemNum << " from element group '"
                    << groupName << "' has not been defined before!");

        elemSet.insert(newElemNum);
      }

      std::copy(elemSet.begin(), elemSet.end(),
                std::back_inserter(elemGroups[groupName]));
    }

    everyThingRead_ = true;
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

} // end of namespace
