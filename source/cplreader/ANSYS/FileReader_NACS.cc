
#include <string>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <iomanip>
#include <sstream> 

#include <boost/regex.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/filesystem/convenience.hpp>
namespace fs=boost::filesystem;

#include "../settings.hh"
#include "FileReader_NACS.hh"


namespace CoupledField
{

  FileReader_NACS::FileReader_NACS(const std::string& name,
                                     const UInt dim,
                                     const UInt numFiles,
                                     const UInt startIndex) :
    FileReader_ANSYS(name, dim, numFiles, startIndex)
  {
  }

  FileReader_NACS::~FileReader_NACS()
  {
  }

  void FileReader_NACS::Init()
  {
    Settings& settings = Settings::Instance();
    std::stringstream sstr;
    std::string line;
    std::ostringstream errMsg;
    double x, y, z;
    UInt nodeNum = 1;
    UInt maxNodeNum = 0;
    UInt numNodes = 0;

    strict_ = (bool) settings.GetInt("strict");
    degen_ = (bool) settings.GetInt("degen");
    
    sstr.clear(); sstr.str("");
    sstr << settings.GetString("basedir") << "/" << name_ << "_Coordinates.dat";

    OpenFile(sstr.str());
    
    while(GetNextLine(line)) {
      x = y = z = 0;
      
      // strip whitespaces
      boost::trim(line);

      if(line == "")
        continue;
      
      sstr.clear(); sstr.str("");
      sstr << line;
      
      sstr >> x >> y >> z;
      
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
      nodeNum++;
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

      std::cerr << "maxNodeNum " << maxNodeNum << std::endl;
      std::cerr << "numNodes " << numNodes << std::endl;
    }

    // Read element types
    UInt elemType;
    UInt elemNum = 1;
    UInt elemDim;
    UInt dim = 0;
    sstr.clear(); sstr.str("");
    sstr << settings.GetString("basedir") << "/" << name_ << "_eltype.dat";

    OpenFile(sstr.str());

    while(GetNextLine(line)) {
      // strip whitespaces
      boost::trim(line);

      sstr.clear(); sstr.str("");
      sstr << line;
      sstr >> elemType;

      elemTypes_[elemNum++] = elemType;
    }

    // Read element connectivity
    std::vector<UInt> elemNodes(40);
    std::set<UInt> elemNodeSet;
    UInt numElemNodes;
    elemNum = 1;
    
    sstr.clear(); sstr.str("");
    sstr << settings.GetString("basedir") << "/" << name_ << "_elconn.dat";

    OpenFile(sstr.str());

    while(GetNextLine(line)) {
      std::fill(elemNodes.begin(), elemNodes.end(), 0);
      elemNodeSet.clear();

      sstr.clear(); sstr.str("");
      sstr << line;
      sstr >> elemNodes[0] >> elemNodes[1] >> elemNodes[2] >> elemNodes[3] 
           >> elemNodes[4] >> elemNodes[5] >> elemNodes[6] >> elemNodes[7]
           >> elemNodes[8] >> elemNodes[9] >> elemNodes[10] >> elemNodes[11] 
           >> elemNodes[12] >> elemNodes[13] >> elemNodes[14] >> elemNodes[15]          
           >> elemNodes[16] >> elemNodes[17] >> elemNodes[18] >> elemNodes[19];

      // Determine number of element nodes by inserting nodes into a set.
      elemNodeSet.insert(elemNodes.begin(), elemNodes.end());
      elemNodeSet.erase((UInt) 0);
      numElemNodes = elemNodeSet.size();

      maxNumElemNodes_ = numElemNodes > maxNumElemNodes_ ? numElemNodes : maxNumElemNodes_;

      FEType newFEType = DegenTypeToNativeType(elemTypes_[elemNum], numElemNodes);
      
      if(degen_)
      {
        DegenerateElement((FEType)elemTypes_[elemNum], newFEType, elemNodes);
      }
      else
      {
        ResortNodes(elemNodes);
      }
      
      elemTypes_[elemNum] = newFEType;
      elemNumsMap_[elemNum] = elemNum;
      elemDim = ELEM_DIM[newFEType];
      dim = dim < elemDim ? elemDim : dim;

      std::copy(elemNodes.begin(), elemNodes.begin() + NUM_ELEM_NODES[newFEType],
                std::back_inserter(topology_[elemNum]));

      elemNum++;
    }

    settings.SetInt("dim", dim);

    GetRegionAndGroupNames();
    
    std::vector<std::string>::const_iterator it, end;
    
    it = regionNames_.begin();
    end = regionNames_.end();
    numRegions_ = std::distance(it, end);
    
    for(UInt i=0; it != end; it++, i++)
    {
      sstr.clear(); sstr.str("");
      sstr << settings.GetString("basedir") << "/" << name_
           << "_RegionElems_" << (*it) << ".dat";
      
      OpenFile(sstr.str());

      while(GetNextLine(line)) {
        boost::trim(line);
        
        sstr.clear(); sstr.str("");
        sstr << line;
        sstr >> elemNum;

        elemRegionMap_[elemNum] = i;
        elemNumsOrig_[i].push_back(elemNum);
      }
    }
  }

  void FileReader_NACS::GetNodeGroups(std::map<std::string,
                                 std::vector<UInt> >& nodeGroups)
  {
    Settings& settings = Settings::Instance();
    std::vector<std::string> nodeFiles;
    std::string line;
    std::string groupName;
    UInt nodeNum;
    std::stringstream sstr;
    std::set<UInt> nodeSet;

    std::vector<std::string>::const_iterator it, end;
    
    it = nodeGroupNames_.begin();
    end = nodeGroupNames_.end();
    numRegions_ = std::distance(it, end);
    
    for( ; it != end; it++ )
    {
      sstr.clear(); sstr.str("");

      sstr << settings.GetString("basedir") << "/" << name_
           << "_NamedNodes_" << (*it) << ".dat"; 
       
      nodeFiles.push_back(sstr.str());
    }
    
    for(UInt i=0, n=nodeFiles.size(); i < n; i++)
    {
      try 
      {
        OpenFile(nodeFiles[i]);
      } catch (std::exception& ex) 
      {
        continue;
      }
      
      groupName = nodeGroupNames_[i];

      nodeSet.clear();
      
      while(GetNextLine(line)) {
        // strip whitespaces
        boost::trim(line);

        if(line.length() == 0)
          continue;

        sstr.clear(); sstr.str("");
        sstr << line;
        sstr >> nodeNum;

        nodeSet.insert(nodeNum);
      }

      std::copy(nodeSet.begin(), nodeSet.end(),
                std::back_inserter(nodeGroups[groupName]));
    }
  }

  void FileReader_NACS::GetElemGroups(std::map<std::string,
                                 std::vector<UInt> >& elemGroups)
  {
    Settings& settings = Settings::Instance();
    std::vector<std::string> elemFiles;
    std::string line;
    std::string groupName;
    UInt elemNum;
    std::stringstream sstr;
    std::set<UInt> elemSet;

    std::vector<std::string>::const_iterator it, end;
    
    it = elemGroupNames_.begin();
    end = elemGroupNames_.end();
    numRegions_ = std::distance(it, end);
    
    for( ; it != end; it++ )
    {
      sstr.clear(); sstr.str("");

      sstr << settings.GetString("basedir") << "/" << name_
           << "_NamedElemsElems_" << (*it) << ".dat"; 
       
      elemFiles.push_back(sstr.str());
    }
    
    for(UInt i=0, n=elemFiles.size(); i < n; i++)
    {
      try 
      {
        OpenFile(elemFiles[i]);
      } catch (std::exception& ex) 
      {
        continue;
      }
      
      groupName = elemGroupNames_[i];

      elemSet.clear();
      
      while(GetNextLine(line)) {
        // strip whitespaces
        boost::trim(line);

        if(line.length() == 0)
          continue;

        sstr.clear(); sstr.str("");
        sstr << line;
        sstr >> elemNum;

        elemSet.insert(elemNum);
      }

      std::copy(elemSet.begin(), elemSet.end(),
                std::back_inserter(elemGroups[groupName]));
    }

    everyThingRead_ = true;
  }

  FEType FileReader_NACS::DegenTypeToNativeType(UInt type, UInt numNodes)
  {
    FEType ret = (FEType) type;
    
    switch((FEType) type)
    {
    case ET_QUAD4: // rectangle
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
      
    case ET_QUAD8: // quad. rectangle
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

    case ET_HEXA8: // hexa
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
      
    case ET_HEXA20: // quad. hexa
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

    default:
      break;
    }

    
    return ret;
  }

  void FileReader_NACS::GetRegionAndGroupNames()
  {
    Settings& settings = Settings::Instance();    
    std::stringstream sstr;
    std::string name;
    
    const boost::regex datExp("\\.dat");
    const boost::regex prefixExp(name_);
    const std::string name_format("");

    sstr << settings.GetString("basedir");
    fs::path dir = sstr.str();
    fs::directory_iterator end_iter;
    for ( fs::directory_iterator dir_itr( dir );
        dir_itr != end_iter;
        ++dir_itr ) 
    {

      if ( !fs::is_directory( *dir_itr ) )
      {
        std::string fn = dir_itr->leaf();

        name = fs::basename(fn);
        name = regex_replace(name, datExp, name_format, boost::match_default | boost::format_sed);

        // Just allow .dat files
        if(fn == name)
          continue;
        
        name = regex_replace(name, prefixExp, name_format, boost::match_default | boost::format_sed);

        if(name.substr(0, 17) == "_NamedElemsElems_")
        {
          name = name.substr(17, name.length());
          
          std::cout << "Named elements: " << name << std::endl;
          elemGroupNames_.push_back(name);
        }

        if(name.substr(0, 12) == "_NamedNodes_")
        {
          name = name.substr(12, name.length());
          
          std::cout << "Named nodes: " << name << std::endl;
          nodeGroupNames_.push_back(name);
        }

        if(name.substr(0, 13) == "_RegionElems_")
        {
          name = name.substr(13, name.length());
          std::cout << "Region: " << name << std::endl;
          regionNames_.push_back(name);
        }
      }
    }
  }
  
} // end of namespace
