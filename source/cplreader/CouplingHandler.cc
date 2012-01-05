// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <fstream>
#include <iostream>
#include <iterator>
#include <set>
#include <string>
#include <algorithm>
#include <cmath>

#include "boost/tokenizer.hpp"
#include "boost/filesystem/operations.hpp"
#include "boost/filesystem/path.hpp"
#include "boost/filesystem/convenience.hpp"
#include "boost/filesystem/exception.hpp"
#include "boost/algorithm/string/trim.hpp"
namespace fs=boost::filesystem;

#include "General/exception.hh"
#include "Settings.hh"
#include "General/environment.hh"
#include "integlib/elemIntegr.hh"

#include "FlowDataTypes.hh"
#include "OutputWriter.hh"
#include "CouplingHandler.hh"

#ifdef _OPENMP
#include "omp.h"
#endif

namespace CoupledField
{

  CouplingHandler::CouplingHandler(shared_ptr<FileReader> ptFileReader,
                                   OutputWriterVectorType& outputWriters)
  {
    ptFileReader_ = ptFileReader;
    outputWriters_ = outputWriters;

    if(!ptFileReader_)
      EXCEPTION("Invalid pointer to file reader!");
  }

  CouplingHandler::~CouplingHandler()
  {
  }

  void CouplingHandler::Init(int argc, char *argv[])
  {
    Settings& settings = Settings::Instance();

    ptFileReader_->Init();

    dim_ = settings.GetInt("dim");
    doIntAverageCentre_ = settings.GetInt("doIntAverageCentre");
    reduce_elementOrder_ = settings.GetInt("reduceElementOrder");

    OutputWriterVectorType::iterator it, end;
    it = outputWriters_.begin();
    end = outputWriters_.end();

    for( ; it != end; it++) 
    {
      (*it)->Init(argc, argv,
                  ptFileReader_->GetPreferredOutputPath());
    }
    
    typedef boost::tokenizer< boost::char_separator<char> > Tok;
    boost::char_separator<char> sep(";| ");

    // Initialize vector with output fields
    Tok tokenizer(settings.GetString("outputfields"), sep);
    std::copy(tokenizer.begin(), tokenizer.end(),
              std::back_inserter(outputFields_));

    // Initialize vector with active regions
    Tok actp(settings.GetString("activeparts"), sep);
    Tok::iterator tit, tend;
    tit = actp.begin();
    tend = actp.end();

    activeParts_.resize(ptFileReader_->GetNumRegions());
    if(*tit == "all")
      std::fill(activeParts_.begin(), activeParts_.end(), true);
    else
    {
      std::stringstream sstr;
      UInt partIdx;

      for( ; tit != tend; tit++)
      {
        sstr.clear(); sstr.str("");
        sstr << *tit;
        sstr >> partIdx;

        if(partIdx > 0 && partIdx <= activeParts_.size())
          activeParts_[partIdx-1] = true;
      }
    }

    if(settings.GetInt("verbose"))
    {
      for(UInt i=0; i < activeParts_.size(); i++)
        std::cout << "Partition " << (i+1) << " active: "
                  << activeParts_[i] << std::endl;
    }

    // Initialize element integrators for source term calculation
    ptElemIntegr_[Elem::LINE2]  = new ElemIntegr(Elem::LINE2);
    ptElemIntegr_[Elem::TRIA3]  = new ElemIntegr(Elem::TRIA3);
    ptElemIntegr_[Elem::QUAD4]  = new ElemIntegr(Elem::QUAD4);
    ptElemIntegr_[Elem::QUAD8]  = new ElemIntegr(Elem::QUAD8);
    ptElemIntegr_[Elem::TET4]   = new ElemIntegr(Elem::TET4);
    ptElemIntegr_[Elem::WEDGE6] = new ElemIntegr(Elem::WEDGE6);
    ptElemIntegr_[Elem::PYRA5]  = new ElemIntegr(Elem::PYRA5);
    ptElemIntegr_[Elem::HEXA8]  = new ElemIntegr(Elem::HEXA8);

  }

  void CouplingHandler::ConvertMesh()
  {
    std::cout << "========================================"
              << "========================================"
              << std::endl;
    std::cout << "                        "
              << "Entering Mesh Conversion Phase!" << std::endl;
    std::cout << "========================================"
              << "========================================"
              << std::endl;


    Settings& settings = Settings::Instance();
    std::cout << "Reading " << settings.GetString("type") << " mesh...\n";
    uint32_t dim = settings.GetInt("dim");
    int numRegions = ptFileReader_->GetNumRegions();
    std::vector<UInt>::iterator it, it2, end;
    std::set<UInt> regionNodeSet;
    UInt maxNumElemNodes = 0;
    std::vector<std::string> regionNames;

    numRegionNodes_.resize(numRegions);

    // First read everything into internal buffers
    std::vector<Double> nodalCoordsTmp;
    ptFileReader_->ReadNodalCoords(nodalCoordsTmp);
    ptFileReader_->ReadTopology(topology_, elemTypes_);
    if (reduce_elementOrder_)
    {
      ptFileReader_->CorrectNumbering(nodalCoordsTmp, topology_, elemTypes_);
    }

    // Determine the maximum number of element nodes
    maxNumElemNodes = ptFileReader_->GetMaxNumElemNodes();
    std::set<UInt> topoSet(topology_.begin(), topology_.end());
    if(*topoSet.begin() == 0)
      topoSet.erase(topoSet.begin());

    // Throw away unused nodes
    std::set<UInt>::iterator topoIt, topoEnd;
    topoIt = topoSet.begin();
    topoEnd = topoSet.end();
    
    std::map<UInt, UInt> pointMap;
    nodalCoords_.resize(topoSet.size()*3);
    for( UInt i=0; topoIt != topoEnd; topoIt++, i++ ) 
    {
      pointMap[*topoIt] = i+1;
      
      UInt idxNew=i*3;
      UInt idxOld=(*topoIt-1)*3;
      nodalCoords_[idxNew+0] = nodalCoordsTmp[idxOld+0];
      nodalCoords_[idxNew+1] = nodalCoordsTmp[idxOld+1];
      nodalCoords_[idxNew+2] = nodalCoordsTmp[idxOld+2];
    }
    // scale the nodal coordinates
    const UInt sizeNodCoords = nodalCoords_.size();
    std::stringstream geomstr;
    geomstr << settings.GetString("geomscale");
    Double geomScaleX, geomScaleY, geomScaleZ;
    geomstr >> geomScaleX >> geomScaleY >> geomScaleZ;
    if (geomScaleX != 1.0)
    {
      UInt i = 0;
      while (i < sizeNodCoords)
      {
        nodalCoords_[i] *= geomScaleX;
        i += 3;
      }
    }
    if (geomScaleY != 1.0)
    {
      UInt i = 1;
      while (i < sizeNodCoords)
      {
        nodalCoords_[i] *= geomScaleY;
        i += 3;
      }
    }
    if (geomScaleZ != 1.0)
    {
      UInt i = 2;
      while (i < sizeNodCoords)
      {
        nodalCoords_[i] *= geomScaleZ;
        i += 3;
      }
    }
    // <-- end scaling
    

    for( UInt i=0, n=topology_.size(); i<n; i++ ) 
    {
      topology_[i] = topology_[i] == 0 ? 0 : pointMap[topology_[i]];
    }

    std::cout << "Writing nodal coords... ";

    OutputWriterVectorType::iterator owIt, owEnd;
    owIt = outputWriters_.begin();
    owEnd = outputWriters_.end();

    for( ; owIt != owEnd; owIt++) 
    {
      (*owIt)->SetDim(dim);
      (*owIt)->WriteNodalCoords(nodalCoords_);
    }

    std::cout << "done.\nWriting element connectivity... ";

    owIt = outputWriters_.begin();
    owEnd = outputWriters_.end();

    for( ; owIt != owEnd; owIt++) 
    {
      (*owIt)->WriteTopology(maxNumElemNodes,
                             elemTypes_,
                             topology_);
    }

    std::cout << "done.\n";


    std::cout << "Reading mesh done.\nConverting mesh...\n";

    for (int actRegion=0; actRegion<numRegions; actRegion++)
    {
      // Fill vector with region names
      regionNames.push_back(ptFileReader_->GetRegionName(actRegion));
      const std::string& actRegionName = regionNames[actRegion];

      // Clear temporary containers
      regionNodeSet.clear();
      regionNodes_[actRegionName].clear();

      std::cout << "Writing region " << (*regionNames.rbegin())
                << "... ";

      // Get all element numbers in current region and order them
      ptFileReader_->GetRegionElements(regionElems_[actRegion], actRegion);
      std::sort(regionElems_[actRegion].begin(), regionElems_[actRegion].end());

      // Iterate over all elements of region and check if all of them
      // have the same dimension
      std::vector<UInt>::const_iterator elemIt, elemEnd;
      elemIt = regionElems_[actRegion].begin();
      elemEnd = regionElems_[actRegion].end();
      
      for( ; elemIt != elemEnd; elemIt++ ) 
      {        
        const UInt& elemNum = *elemIt;
        UInt& regionDim = regionDims_[actRegion];
        Elem::FEType elemType = (Elem::FEType) elemTypes_[elemNum -1];
        
        if (!regionDim) 
        {
          regionDim = Elem::GetElemDim(elemType);
        }
        if ( regionDim != Elem::GetElemDim(elemType) )
        {
          EXCEPTION("Elements with different dimensions have been "
              << "encountered in region '" << (*regionNames.rbegin()) << "'!\n"
              << "The error occured while examining element "
              << elemNum << ".\n"
              << "Please check your mesh file!");
        }    
      }
      
      // Put all nodes in a partition into a set to get an ordered list
      UInt idx = 0;
      UInt regionDim = 0;

      CollectElementNodes(regionElems_[actRegion], regionNodes_[actRegionName], regionDim);

      // Write region elements and nodes
      owIt = outputWriters_.begin();
      owEnd = outputWriters_.end();
      
      for( ; owIt != owEnd; owIt++) 
      {
        (*owIt)->WriteRegion(regionNames[actRegion],
                             regionElems_[actRegion],
                             regionNodes_[actRegionName],
                             regionDim);
      }
      
      it = regionNodes_[actRegionName].begin();
      end = regionNodes_[actRegionName].end();

      numRegionNodes_[actRegion] = std::distance(it, end);

      for ( idx=0; it != end; it++, idx++ )
      {
        regionNodeIndices_[actRegion][*it] = idx;
      }


      std::cout << "done.\n";
      if(settings.GetInt("verbose"))
        std::cout << "Number of nodes in region " << regionNames[actRegion]
                  << " " << numRegionNodes_[actRegion] << std::endl;

    }

    // Read node and element groups
    ptFileReader_->GetNodeGroups(nodeGroups_);
    ptFileReader_->GetElemGroups(elemGroups_);

    // Collect nodes in element groups
    std::string elemsName;
    std::map<std::string, std::vector<UInt> >::const_iterator egIt, egEnd;
    std::map<std::string, std::vector<UInt> > elemGroupNodes;
    std::map<std::string, UInt > elemGroupDims;
    egIt = elemGroups_.begin();
    egEnd = elemGroups_.end();

    for(; egIt != egEnd; egIt++ ) {
      elemsName = egIt->first;
      CollectElementNodes(egIt->second,
                          elemGroupNodes[elemsName],
                          elemGroupDims[elemsName]);
    }
    
    // Write node and element groups
    owIt = outputWriters_.begin();
    owEnd = outputWriters_.end();
    
    for( ; owIt != owEnd; owIt++) 
    {
      (*owIt)->WriteNodeGroups(nodeGroups_);
      (*owIt)->WriteElemGroups(elemGroups_,
                               elemGroupNodes,
                               elemGroupDims);
    }

    std::cout << "Converting mesh done.\n";

    std::cout << "========================================"
              << "========================================"
              << std::endl;
    std::cout << "                        "
              << "Leaving Mesh Conversion Phase!" << std::endl;
    std::cout << "========================================"
              << "========================================"
              << std::endl;
  }

  void CouplingHandler::Couple()
  {
    Settings& settings = Settings::Instance();
    bool calcSrc = settings.GetInt("calcsrc") != 0;
    UInt stepInc = settings.GetInt("stepincr");
    UInt counter = 0;
    Double stepVal = 0;
    UInt numFiles = ptFileReader_->GetNumFiles();
    UInt numRegions = ptFileReader_->GetNumRegions();
    std::vector<UInt> timeStepNumbers;
    std::vector<Double> timeStepValues;
    UInt actRegion;
    UInt stepNum = 0;
    std::vector<std::string> regionNames;
    std::vector<FlowDataType> flowData(numRegions);
    bool readOK = true;
    // variable to gather nodes on multiple regions
    bool doCalcMultiNodes = settings.GetInt("doCalcMultiNodes");

    std::map<UInt, std::map<std::string, UInt> > multiNodes;

    std::cout << "========================================"
              << "========================================"
              << std::endl;
    std::cout << "                        "
              << "Entering Time Data Conversion Phase!" << std::endl;
    std::cout << "                        ";
    std::cout << "Number of transient files: " << numFiles << std::endl;
    std::cout << "========================================"
              << "========================================"
              << std::endl;

    OutputWriterVectorType::iterator owIt, owEnd;

    if(numFiles) 
    {
      owIt = outputWriters_.begin();
      owEnd = outputWriters_.end();
      
      for( ; owIt != owEnd; owIt++) 
      {
        (*owIt)->BeginResults(&flowData);
      }
      
      // Fill vector with region names
      for (actRegion = 0; actRegion<numRegions; actRegion++)
        regionNames.push_back(ptFileReader_->GetRegionName(actRegion));
      
      //at this point we can already calculate the surface region neighbours
      PrepareSurfaceRegions();

      while ( ( counter < numFiles*stepInc ) && readOK)
      {
        stepVal = ptFileReader_->GetTimeStep(counter);
        stepNum = counter + ptFileReader_->GetStartIndex();
        timeStepValues.push_back(stepVal);
        timeStepNumbers.push_back(stepNum);
        
        std::cout << "----------------------------------------"
                  << "----------------------------------------"
                  << std::endl;
        std::cout << "                        "
                  << "Step Number: " << stepNum << "; ";
        std::cout << " Time Step: " << stepVal << std::endl;
        std::cout << "----------------------------------------"
                  << "----------------------------------------"
                  << std::endl;
        
        // Tell output writers about begin step
        owIt = outputWriters_.begin();
        owEnd = outputWriters_.end();
        
        for( ; owIt != owEnd; owIt++) 
        {
          (*owIt)->BeginStep(stepNum, stepVal);
        }
        
        // Read nodal values for all partitions
        try
        {
          ptFileReader_->ReadNodalValues(flowData, activeParts_, counter);
          if (reduce_elementOrder_)
          {
            ptFileReader_->ReduceOrderOfNodalValues(flowData, regionNodes_);
          }
          // scale the nodal values
          // following physical fields will be checked for scaling factors
          const std::string physFieldScale_str[] = {"velscale","geomscale"};
          const SolutionType solType[] = {FLUIDMECH_VELOCITY,MECH_DISPLACEMENT};
          Double scaleX, scaleY, scaleZ;
          for (int i = 0; i < 2; ++i) // 2 physical fields
          {
            std::stringstream scaleStr;
            scaleStr << settings.GetString(physFieldScale_str[i]);

            scaleStr >> scaleX;
            if (!scaleStr.fail() && scaleX != 1.0)
              scale_PhysField_(flowData, solType[i], scaleX, 0);

            scaleStr >> scaleY;
            if (!scaleStr.fail() && scaleY != 1.0)
              scale_PhysField_(flowData, solType[i], scaleY, 1);

            if ( dim_ == 3 ) {
              scaleStr >> scaleZ;
              if (!scaleStr.fail() && scaleZ != 1.0)
                scale_PhysField_(flowData, solType[i], scaleZ, 2);
            }
          }

          // <-- end scaling velocity
          
          // Override the setting of --outprec for CFX
          if(settings.GetString("type") == "CFX" && settings.GetInt("floatDataset"))
            settings.SetString("outprec", "single");
          
          readOK = true;
        } catch (std::exception& ex)
        {
          std::cerr << "CAUGHT EXCEPTION while trying to read nodal values:"
                    << std::endl << ex.what() << std::endl
                    << "Exiting read time values loop..." << std::endl;
          
          readOK = false;
          continue;
        }
        
        
        for (actRegion = 0; actRegion<numRegions && readOK; actRegion++)
        {
          FlowDataType::iterator fdIt, fdEnd;
          
          if(!activeParts_[actRegion])
          {
            FlowDataType::iterator fIt, fEnd;
            fIt = flowData[actRegion].begin();
            fEnd = flowData[actRegion].end();
            
            for( ; fIt != fEnd; fIt++ ) {
              fIt->second.isActive = false;
            }
            
            continue;
          }


          if(calcSrc && regionDims_[actRegion] == dim_)
          {
            // We need fluidMechVelocity for Lighthill source term.
            // This must be adapted for other source term formulations!!
            if ( flowData[actRegion].find(FLUIDMECH_VELOCITY)
                 != flowData[actRegion].end() )
            {
              flowData[actRegion][FLUIDMECH_VELOCITY].isActive = true;
            }
            CalculateAcouSrcs(actRegion, flowData[actRegion]);
          }
          
        }//end of for volume regions

        //loop over all requested surface regions
        std::map<UInt,bool>::iterator sIt = surfaceRegionIndices_.begin();
        for(;sIt!=surfaceRegionIndices_.end();sIt++){
          if(sIt->second)
            CalculateSurfaceIntegral(sIt->first,flowData);
        }

        /* node which live on multiple region need to accumulate the
         * ACOU_RHS_LOAD
         * Here the common nodes of active regions is found */
        if (doCalcMultiNodes)
        {
          std::map<std::string, std::vector<UInt>* > regionNodesActive;
          std::map<std::string, std::vector<UInt> >::iterator iterRegionAct = regionNodes_.begin();
          for (; iterRegionAct != regionNodes_.end(); ++iterRegionAct)
          {
            UInt regIdx = 0;
            while (regionNames[regIdx] != iterRegionAct->first )
              ++regIdx;
            if ( flowData[regIdx].find(ACOU_RHS_LOAD)
                 != flowData[regIdx].end() )
            {
              if (flowData[regIdx][ACOU_RHS_LOAD].isActive)
              {
                regionNodesActive[iterRegionAct->first] = &iterRegionAct->second;
              }
            }
          }
          findNodeMultiRegion(regionNodesActive, multiNodes);
          doCalcMultiNodes = false;
        }
        /* start the accumulation of ACOU_RHS_LOAD on common nodes */
        std::map<UInt, std::map<std::string, UInt> >::const_iterator iterMultiNodes = multiNodes.begin();
        if (calcSrc)
        {
          for (; iterMultiNodes != multiNodes.end(); ++iterMultiNodes)
          {
            // calc accumulated values
            double accumValNodes = 0;
            std::map<std::string, UInt>::const_iterator iterRegions = iterMultiNodes->second.begin();
            for (; iterRegions != iterMultiNodes->second.end(); ++iterRegions)
            {
              const std::string& regName = iterRegions->first;
              const UInt& node = iterRegions->second;
              UInt regIdx = 0;
              while (regionNames[regIdx] != regName )
                ++regIdx;
              if ( flowData[regIdx].find(ACOU_RHS_LOAD)
                  != flowData[regIdx].end() )
              {
                if (flowData[regIdx][ACOU_RHS_LOAD].isActive)
                  accumValNodes += flowData[regIdx][ACOU_RHS_LOAD].data[node];
              }
            }
            // set accumulated values
            iterRegions = iterMultiNodes->second.begin();
            for (; iterRegions != iterMultiNodes->second.end(); ++iterRegions)
            {
              const std::string& regName = iterRegions->first;
              const UInt& node = iterRegions->second;
              UInt regIdx = 0;
              while (regionNames[regIdx] != regName )
                ++regIdx;
              if ( flowData[regIdx].find(ACOU_RHS_LOAD)
                  != flowData[regIdx].end() )
              {
                if (flowData[regIdx][ACOU_RHS_LOAD].isActive)
                  flowData[regIdx][ACOU_RHS_LOAD].data[node] = accumValNodes;
              }
            }
          }
        }// end of nodes on multiple region correction

        for (actRegion = 0; actRegion < numRegions && readOK; actRegion++)
        {
          owIt = outputWriters_.begin();
          owEnd = outputWriters_.end();
          
          for( ; owIt != owEnd; owIt++) 
          {
            (*owIt)->WriteFlowData(this,
                                   actRegion,
                                   outputFields_);
          }
          
        }//end of for
        
        owIt = outputWriters_.begin();
        owEnd = outputWriters_.end();
        
        for( ; owIt != owEnd; owIt++) 
        {
          (*owIt)->EndStep();
        }
        
        // increment time step counter
        counter = counter + stepInc;
      }//end of while

      owIt = outputWriters_.begin();
      owEnd = outputWriters_.end();
      
      for( ; owIt != owEnd; owIt++) 
      {
        (*owIt)->EndResults();
      }
    }

    owIt = outputWriters_.begin();
    owEnd = outputWriters_.end();
      
    for( ; owIt != owEnd; owIt++) 
    {
      (*owIt)->DoneWriting();
    }

    std::cout << "========================================"
              << "========================================"
              << std::endl;
    std::cout << "                        "
              << "Leaving Time Data Conversion Phase!" << std::endl;
    std::cout << "========================================"
              << "========================================"
              << std::endl;
  }
  
  void CouplingHandler::Finish()
  {
    Settings& settings = Settings::Instance();

    // Fetch custom data from file reader and write it to the UserData
    // section of the HDF5 file.
    std::map<std::string, std::string> userData;
    ptFileReader_->GetUserData(userData);
    settings.DumpXML(userData["settings"]);
 
    // Write user data
    OutputWriterVectorType::iterator it, end;
    it = outputWriters_.begin();
    end = outputWriters_.end();
    
    for( ; it != end; it++) 
    {
      (*it)->WriteUserData(userData);
    }
    
    // Delete element integrators
    std::map<UInt, ElemIntegr *>::iterator intIt, intEnd;
    intIt = ptElemIntegr_.begin();
    intEnd = ptElemIntegr_.end();
    for( ; intIt != intEnd; intIt++)
      delete intIt->second;
  }

  void CouplingHandler::CollectElementNodes(const std::vector<UInt>& elems,
                                            std::vector<UInt>& nodes,
                                            UInt& dim)
  {
    std::vector<UInt>::const_iterator it, end;
    std::set<UInt> nodeSet;

    UInt topoIdx = 0;
    UInt idx = 0;
    UInt numElemNodes;
    Elem::FEType et;
    UInt maxNumElemNodes = ptFileReader_->GetMaxNumElemNodes();


    dim = 0;
    it = elems.begin();
    end = elems.end();

    for ( ; it != end; it++ )
    {
      idx = (*it-1);
      topoIdx = idx * maxNumElemNodes;
      et = (Elem::FEType) elemTypes_[idx];
      numElemNodes = Elem::GetNumElemNodes( et );

      nodeSet.insert(&topology_[topoIdx],
                     &topology_[topoIdx+numElemNodes]);

      dim = dim < Elem::GetElemDim(et) ? Elem::GetElemDim(et) : dim;
    }

    std::copy(nodeSet.begin(), nodeSet.end(),
              std::back_inserter(nodes));
  }

  void CouplingHandler::CalculateAcouSrcs(const int regionIdx,
                                          FlowDataType& flowData)
  {
    Settings& settings = Settings::Instance();

    std::string regionName = ptFileReader_->GetRegionName(regionIdx);

    if(flowData.find(FLUIDMECH_VELOCITY) == flowData.end())
    {
      std::cerr << "Cannot calculate Lighthill sources on " << regionName
                << " since no velocity field is available!" << std::endl;
      return;
    }

    // Acquire reference to velocity field
    FlowDataPartStruct& fdps = flowData[FLUIDMECH_VELOCITY];
    std::vector<Double>& velField = fdps.data;

    if(!fdps.isActive)
    {
      std::cerr << "Will not calculate Lighthill sources on " << regionName
                << " since velocity field is not active!" << std::endl;
      return;
    }

    std::cout << "Calculating Lighthill sources on " << regionName << " ";
    Double density = settings.GetDouble("density");


    // Init Lighthill structures
    FlowDataPartStruct& fdps2 = flowData[ACOU_RHS_LOAD];
    fdps2.isActive = true; // all partitions have results
    fdps2.definedOn = ResultInfo::NODE; // nodes
    if(fdps2.dofNames.empty())
      fdps2.dofNames.push_back("-");
    fdps2.unit = MapSolTypeToUnit(ACOU_RHS_LOAD);
    fdps2.resultName = "acouRhsLoad";
    fdps2.data.resize(numRegionNodes_[regionIdx]);
    fdps2.entryType = ResultInfo::SCALAR;
    std::vector<Double>& acouRhsField = fdps2.data;

    // Fill acouRhsLoad field with zeros
    std::fill(acouRhsField.begin(), acouRhsField.end(), 0);

    
    FlowDataPartStruct& fdps3 = flowData[ACOU_RHS_LOAD_DENSITY];
    fdps3.isActive = true; // all partitions have results
    fdps3.definedOn = ResultInfo::NODE; // nodes
    if(fdps3.dofNames.empty())
      fdps3.dofNames.push_back("-");
    fdps3.unit = MapSolTypeToUnit(ACOU_RHS_LOAD_DENSITY);
    fdps3.resultName = "acouRhsLoadDensity";
    fdps3.data.resize(numRegionNodes_[regionIdx]);
    fdps3.entryType = ResultInfo::SCALAR;
    std::vector<Double>& acouRhsDensityField = fdps3.data;

    // Fill acouRhsLoadDensity field with zeros
    std::fill(acouRhsDensityField.begin(), acouRhsDensityField.end(), 0);

    int nElems = ptFileReader_->GetNumElems(regionIdx);
    
    FlowDataPartStruct& fdps4 = flowData[ACOU_DIV_LH_TENSOR];
    fdps4.isActive = true; // all partitions have results
    fdps4.definedOn = ResultInfo::ELEMENT; // elements
    if(fdps4.dofNames.empty()) {
      fdps4.dofNames.push_back("x");
      fdps4.dofNames.push_back("y");
      if(dim_ == 3)
        fdps4.dofNames.push_back("z");        
    }
    fdps4.unit = MapSolTypeToUnit(ACOU_DIV_LH_TENSOR);
    fdps4.resultName = "acouDivLighthillTensor";
    fdps4.data.resize(nElems * dim_);
    fdps4.entryType = ResultInfo::VECTOR;
    std::vector<Double>& acouDivLighthillTensor = fdps4.data;

    // Fill acouDivLighthillTensor field with zeros
    std::fill(acouDivLighthillTensor.begin(), acouDivLighthillTensor.end(), 0);

    Matrix<Double> coordMat;
    Matrix<Double> nodaldTijdxj;
    Matrix<Double> nodalVel;
    Vector<Double> elemVec;
    Vector<Double> nodalLoadDensity;
    Vector<Double> divLHTensor(dim_);

    Elem::FEType elemType;
    UInt numElemNodes;
    UInt elemDim;
    UInt elemIdx;
    UInt maxNENodes = ptFileReader_->GetMaxNumElemNodes();
    UInt nodeNum;

#ifdef _OPENMP
    IntegrationMap ptElemI(ptElemIntegr_);
    std::cout << "...parallel loop over " << nElems << " elements..." ;
#else
    std::cout << "...serial loop over " << nElems << " elements..." ;
#endif
    std::cout.flush();
    UInt velIdx,topoIdx,idx;
#pragma omp parallel private(elemIdx,elemType,numElemNodes,elemDim,coordMat, \
                                 nodaldTijdxj,nodalVel,nodeNum, \
                                 elemVec,nodalLoadDensity,velIdx,topoIdx,idx) firstprivate(divLHTensor,ptElemI) 
{
    #pragma omp for nowait
    for( int i=0; i<nElems; i++)
    {
      elemIdx = regionElems_[regionIdx][i] - 1;
      elemType = (Elem::FEType) elemTypes_[elemIdx];
      numElemNodes = Elem::GetNumElemNodes(elemType);
      elemDim = dim_ ; //Elem::GetElemDim(elemType);

      //Just calculate sources for volume elements!
      if(elemDim < dim_)
         continue;

      coordMat.Resize(elemDim, numElemNodes);
      nodaldTijdxj.Resize(elemDim, numElemNodes);
      nodalVel.Resize(elemDim, numElemNodes);

      for( UInt n=0; n<numElemNodes; n++)
      {
        nodeNum = topology_[elemIdx * maxNENodes + n];
        topoIdx = (nodeNum - 1) * 3;
        velIdx = regionNodeIndices_[regionIdx][nodeNum] * dim_;

        for( UInt d=0; d<elemDim; d++)
        {
          coordMat[d][n] = nodalCoords_[topoIdx+d];
          nodalVel[d][n] = velField[velIdx+d];
        }
      }
      if(flowData.find(SMOOTH_DISPLACEMENT) != flowData.end())
      {
        FlowDataPartStruct& dataSmoothDispl = flowData[SMOOTH_DISPLACEMENT];
        std::vector<Double>& mechSmoothVec = dataSmoothDispl.data;
        for( UInt n=0; n<numElemNodes; n++)
        {
          nodeNum = topology_[elemIdx * maxNENodes + n];
          UInt smoothIdx = regionNodeIndices_[regionIdx][nodeNum] * dim_;

          for( UInt d=0; d<elemDim; d++)
          {
            coordMat[d][n] += mechSmoothVec[smoothIdx +d];
          }
        }
      }

      try {
#ifdef _OPENMP
        if (doIntAverageCentre_)
        {
          ptElemI[elemType].PerformIntegrationCentre( coordMat, nodalVel,
                                                      elemVec, 
                                                      nodalLoadDensity, 
                                                      divLHTensor, 
                                                      density);
        } else {
          ptElemI[elemType].PerformIntegration( coordMat, nodaldTijdxj, 
                                                nodalVel, elemVec, 
                                                nodalLoadDensity, 
                                                divLHTensor, density);
        }
#else
        if (doIntAverageCentre_)
        {
          ptElemIntegr_[elemType]->PerformIntegrationCentre( coordMat, nodalVel,
                                                             elemVec, nodalLoadDensity, 
                                                             divLHTensor, density);
        } else {
          ptElemIntegr_[elemType]->PerformIntegration( coordMat, nodaldTijdxj, nodalVel,
                                                       elemVec, nodalLoadDensity, 
                                                       divLHTensor, density);
        }
#endif

      } catch (CoupledField::Exception &ex)
      {
        std::cerr << "WARN: An Exception occurred during source term "
                  << "computation:\nElement " << elemIdx+1 << std::endl;

        std::cerr << ex.what()<< std::endl;

        if(settings.GetInt("verbose")) {
          UInt oldPrec = std::cerr.precision(8);

          std::cerr << "Corner coords:\n";
          for (UInt iCol = 0, numCols = coordMat.GetNumCols();
               iCol < numCols; ++iCol) {
            for (UInt iRow = 0, numRows = coordMat.GetNumRows();
                 iRow < numRows; ++iRow) {
              std::cerr << "\t" << coordMat[iRow][iCol];
            }
            std::cerr << std::endl;
          }

          std::cerr.precision(oldPrec);
          std::cerr << ex.what();
        }

        std::cerr << "Setting contribution to acousrc to zero!\n\n";
        elemVec.Resize(numElemNodes);
        elemVec.Init();
      }

      // Add contributions of all element nodes
      for( UInt n=0; n<numElemNodes; n++)
      {
        nodeNum = topology_[elemIdx * maxNENodes + n];
        idx = regionNodeIndices_[regionIdx][nodeNum];

#ifndef NDEBUG
        if (std::isnan(elemVec[n]) || std::isinf(elemVec[n])) {
          EXCEPTION("Source term calculated on element " << i+1
                    << " is Inf or Nan.");
        }
#endif

#pragma omp atomic 
        acouRhsField[idx] -= elemVec[n];
#pragma omp atomic 
        acouRhsDensityField[idx] -= nodalLoadDensity[n];
      }
      
      // Add contributions of elements
      for( UInt n=0; n < dim_; n++)
      {
        acouDivLighthillTensor[i*dim_ + n] = divLHTensor[n];
      }      
    }
}//end of parallel region

    std::cout << "done." << std::endl;
  }

  void CouplingHandler::CalculateSurfaceIntegral(const int surfRegionIdx,
                                            std::vector<FlowDataType> & flowData){
    //1. Do standard stuff, meaning we initialize the results which make sense for surfaces
    //   Right now, this is only acouRHSLoad
    //2. loop over all elements in this region
    // 2a. obtain velocity field and coordinates for associated volume element
    // 2b. call surface integrator
    // 2c. write result to surfaceRegion and subtract it from volume region

    Settings& settings = Settings::Instance();

    std::string regionName = ptFileReader_->GetRegionName(surfRegionIdx);
    UInt maxNENodes = ptFileReader_->GetMaxNumElemNodes();
    UInt nElems = ptFileReader_->GetNumElems(surfRegionIdx);

    std::cout << "Calculating Lighthill surface sources on " << regionName << " ";
    Double density = settings.GetDouble("density");

    // Init Lighthill structures
    FlowDataPartStruct& fdps2 = flowData[surfRegionIdx][ACOU_RHS_LOAD];
    fdps2.isActive = true; // all partitions have results
    fdps2.definedOn = ResultInfo::NODE; // nodes
    if(fdps2.dofNames.empty())
      fdps2.dofNames.push_back("-");
    fdps2.unit = MapSolTypeToUnit(ACOU_RHS_LOAD);
    fdps2.resultName = "acouRhsLoad";
    fdps2.data.resize(numRegionNodes_[surfRegionIdx]);
    fdps2.entryType = ResultInfo::SCALAR;
    std::vector<Double>& SurfTermField = fdps2.data;

    // Fill acouRhsLoad field with zeros
    std::fill(SurfTermField.begin(), SurfTermField.end(), 0);

    FlowDataPartStruct& fdps4 = flowData[surfRegionIdx][ACOU_DIV_LH_TENSOR];
    fdps4.isActive = true; // all partitions have results
    fdps4.definedOn = ResultInfo::ELEMENT; // elements
    if(fdps4.dofNames.empty()) {
      fdps4.dofNames.push_back("x");
      fdps4.dofNames.push_back("y");
      if(dim_ == 3)
        fdps4.dofNames.push_back("z");
    }
    fdps4.unit = MapSolTypeToUnit(ACOU_DIV_LH_TENSOR);
    fdps4.resultName = "acouDivLighthillTensor";
    fdps4.data.resize(nElems * dim_);
    fdps4.entryType = ResultInfo::VECTOR;
    std::vector<Double>& acouDivLighthillTensor = fdps4.data;

    // Fill acouDivLighthillTensor field with zeros
    std::fill(acouDivLighthillTensor.begin(), acouDivLighthillTensor.end(), 0);

    //Define variables most of which are required for surface as well
    //as for the volume element
    Elem::FEType surfElemType,volElemType;
    UInt numSurfElemNodes,numVolElemNodes;
    UInt surfElemDim,volElemDim;
    UInt surfElemIdx,volElemIdx;
    Matrix<Double> surfCoordMat,volCoordMat;
    Matrix<Double> volNodalVel;
    Vector<Double> surfElemVec;
    Vector<Double> DivLHElemVec;
    UInt volRegionIdx,volElemIdxLocal;

    UInt velIdx,topoIdx,nodeNum;

    std::cout << "...serial loop over " << nElems << " elements..." ;
    std::cout.flush();
    for( UInt i=0; i<nElems; i++)
    {
      surfElemIdx = regionElems_[surfRegionIdx][i] - 1;

      if(surfaceNeighbors_.find(surfElemIdx)==surfaceNeighbors_.end()){
        WARN("COuld not find neighbors for a element! ignoring it...");
        continue;
      }

      //first check if we have a velocity for the current volume element
      volRegionIdx = surfaceNeighbors_[surfElemIdx].first;
      FlowDataPartStruct& volVelocityReg = flowData[volRegionIdx][FLUIDMECH_VELOCITY];
      FlowDataPartStruct& volRHSReg = flowData[volRegionIdx][ACOU_RHS_LOAD];
      FlowDataPartStruct& volDivLHReg = flowData[volRegionIdx][ACOU_DIV_LH_TENSOR];

      //ignore this element if associated volume has no data
      if(!volVelocityReg.isActive)// || !volRHSReg.isActive)
        continue;

      const std::vector<Double>& volVelField = volVelocityReg.data;

      //obtain information about volume and surface
      volElemIdxLocal = surfaceNeighbors_[surfElemIdx].second;
      volElemIdx = regionElems_[volRegionIdx][volElemIdxLocal]-1;
      volElemType = (Elem::FEType) elemTypes_[volElemIdx];
      numVolElemNodes = Elem::GetNumElemNodes(volElemType);
      volElemDim = Elem::GetElemDim(volElemType);

      surfElemType = (Elem::FEType) elemTypes_[surfElemIdx];
      numSurfElemNodes = Elem::GetNumElemNodes(surfElemType);
      surfElemDim =  Elem::GetElemDim(surfElemType);

      //resize structures
      volCoordMat.Resize(volElemDim, numVolElemNodes);
      volNodalVel.Resize(volElemDim, numVolElemNodes);
      surfCoordMat.Resize(volElemDim,numSurfElemNodes);

      //obtain information from volume element
      for( UInt n=0; n<numVolElemNodes; n++)
      {
        nodeNum = topology_[volElemIdx * maxNENodes + n];
        topoIdx = (nodeNum - 1) * 3;
        velIdx = regionNodeIndices_[volRegionIdx][nodeNum] * dim_;
        for( UInt d=0; d<volElemDim; d++)
        {
          volCoordMat[d][n] = nodalCoords_[topoIdx+d];
          volNodalVel[d][n] = volVelField[velIdx+d];
          if(volNodalVel[d][n] > 9999){
            std::cout << volVelField.size() << std::endl;
            std::cout << "problem HERE" << std::endl;
          }
        }
      }
      //obtain information for surface element i.e. coordinates
      for( UInt n=0; n<numSurfElemNodes; n++)
      {
        nodeNum = topology_[surfElemIdx * maxNENodes + n];
        topoIdx = (nodeNum - 1) * 3;
        for( UInt d=0; d<volElemDim; d++)
          surfCoordMat[d][n] = nodalCoords_[topoIdx+d];
      }

      if(flowData[volRegionIdx].find(SMOOTH_DISPLACEMENT) != flowData[volRegionIdx].end())
      {
        WARN("Smooth displacements are not supported right now for surface Elements!");
      }

      //PERFORM INTEGRATION
      try {
        if (doIntAverageCentre_)
        {
          ptElemIntegr_[surfElemType]->PerformSurfaceIntegrationCenter( volElemType, volCoordMat,
              surfCoordMat, volNodalVel, surfaceNormals_[surfElemIdx], density,surfElemVec,DivLHElemVec);
        } else {
          ptElemIntegr_[surfElemType]->PerformSurfaceIntegration(volElemType, volCoordMat,
              surfCoordMat, volNodalVel, surfaceNormals_[surfElemIdx], density,surfElemVec,DivLHElemVec);
        }

      } catch (CoupledField::Exception &ex)
      {
        std::cerr << "WARN: An Exception occurred during surface source term "
                  << "computation:\nElement " << surfElemIdx+1 << std::endl;

        std::cerr << ex.what()<< std::endl;

        if(settings.GetInt("verbose")) {
          UInt oldPrec = std::cerr.precision(8);

          std::cerr << "Corner coords:\n";
          for (UInt iCol = 0, numCols = surfCoordMat.GetNumCols();
               iCol < numCols; ++iCol) {
            for (UInt iRow = 0, numRows = surfCoordMat.GetNumRows();
                 iRow < numRows; ++iRow) {
              std::cerr << "\t" << surfCoordMat[iRow][iCol];
            }
            std::cerr << std::endl;
          }

          std::cerr.precision(oldPrec);
          std::cerr << ex.what();
        }
        std::cerr << "Setting contribution to acousrc to zero!\n\n";
        surfElemVec.Resize(numSurfElemNodes);
        surfElemVec.Init();
      }

      // Add contributions of all element nodes
      for( UInt n=0; n<numSurfElemNodes; n++)
      {
          nodeNum = topology_[surfElemIdx * maxNENodes + n];
          UInt volIdx = regionNodeIndices_[volRegionIdx][nodeNum];
          UInt surfIdx = regionNodeIndices_[surfRegionIdx][nodeNum];

          if(volRHSReg.isActive)
            volRHSReg.data[volIdx] += surfElemVec[n];

          SurfTermField[surfIdx] += surfElemVec[n];
      }
      // Add contributions of elements
      for( UInt n=0; n < dim_; n++)
      {
        acouDivLighthillTensor[i*dim_ + n] = DivLHElemVec[n];
        if(volDivLHReg.isActive){
          volDivLHReg.data[volElemIdxLocal*dim_ + n] -= DivLHElemVec[n];
        }
      }
    }

    std::cout << "done." << std::endl;
  }

  void CouplingHandler::ShrinkNodalVector(const UInt partitionIdx,
                                          const UInt numDOFs,
                                          const std::vector<Double>& input,
                                          std::vector<Double>& output)
  {
    UInt numInputNodes = input.size() / numDOFs;
    UInt numNodes = numRegionNodes_[partitionIdx];
    std::map<UInt, UInt>::const_iterator it, end;
    UInt idxInput = 0;
    UInt idxOutput = 0;

    if ( numInputNodes == numNodes )
    {
      output = input;
      return;
    }
    else if ( numInputNodes < numNodes )
    {
      EXCEPTION( "Cannot shrink vector (length " << numInputNodes
                 << "), because there are too much nodes ("
                 << numNodes << ") in region " << partitionIdx );
    }

    it = regionNodeIndices_[partitionIdx].begin();
    end = regionNodeIndices_[partitionIdx].end();
    output.resize(numNodes * numDOFs);

    for( ; it != end; it++ )
    {
      idxInput = (it->first - 1) * numDOFs;
      idxOutput = it->second * numDOFs;

      for(UInt dof=0; dof < numDOFs; dof++ )
      {
        output[idxOutput+dof] = input[idxInput+dof];
      }
    }
  }

  void CouplingHandler::GetNeighborElements(std::map<UInt,std::pair<UInt,UInt> >& volNeighbors,UInt surfRegionIdx){
    std::string regionName = ptFileReader_->GetRegionName(surfRegionIdx);
    std::cout << "Calculating Neighbor element information for region " << regionName << " ....";
    std::cout.flush();
    UInt numElems = regionElems_[surfRegionIdx].size();
    UInt maxNENodes = ptFileReader_->GetMaxNumElemNodes();
    std::vector<UInt> matches;
    matches.clear();
    UInt curNodeNum = 0;
    UInt curElemNum = 0;

    for(UInt i = 0; i< numElems; i++){
      matches.clear();
      //get number of element
      curElemNum = regionElems_[surfRegionIdx][i] - 1;
      //get number of first node
      curNodeNum = topology_[curElemNum * maxNENodes];
      //now we search through the topology vector and determine
      //all occurences of the nodeNumber and store their position
      std::vector<UInt>::iterator i = topology_.begin(), end = topology_.end();
      while(true) {
        i = std::find(i, topology_.end(), curNodeNum );
        if (i == end)
            break;
        matches.push_back(i - topology_.begin());
        i++;
      }
      //now we have a vector containing all positions of the first node of the
      //element of interest now we need to determine the element indices
      //and store only those connected to volume elements, i.e. dim_ == elemDim
      std::vector<UInt> globElemIdx;
      globElemIdx.clear();
      for(UInt j=0;j<matches.size();j++){
        UInt eNum = std::floor((Double)matches[j] / (Double)maxNENodes);
        Elem::FEType elemType = (Elem::FEType) elemTypes_[eNum];
        if(Elem::GetElemDim(elemType) == dim_ && CalcSurfaceNormalOutVolume(eNum,curElemNum))
          globElemIdx.push_back(eNum);
      }
      //if we have more than one element in globalElemIdx we assume an internal surface
      //so we set the surface region to false and return without doing anything
      if(globElemIdx.size() != 1){
        surfaceRegionIndices_[surfRegionIdx] = false;
        WARN("Found more than one neighbor for surface element in surface region we set it to false");
        return;
      }

      //now we have the candidates now we determine the real element
      //to do this unfortunately we have to search through every regionelem again
      //to obtain the region local element number
      for(UInt curReg=0;curReg<volumeRegionIndices_.size();curReg++){
        //ok now we search for the stuff
        std::vector<UInt>::iterator g = regionElems_[curReg].begin(), end = regionElems_[curReg].end();
        g = std::find(g, regionElems_[curReg].end(), globElemIdx[0]+1 );
        if (g == end)
           continue;
        else{
          UInt localVolEIdx = g-regionElems_[curReg].begin();
          volNeighbors[curElemNum] = std::pair<UInt,UInt>(curReg,localVolEIdx);
          break;
        }
      }
    }
    std::cout << "done" << std::endl;
  }

  void CouplingHandler::PrepareSurfaceRegions(){
    UInt numRegions = ptFileReader_->GetNumRegions();

    //obtain the list of surface regions or the all word
    std::vector<std::string> surfaceRegions;

    typedef boost::tokenizer< boost::char_separator<char> > Tok;
    boost::char_separator<char> sep(";| ");

    // Initialize vector with output fields
    Settings& settings = Settings::Instance();
    Tok tokenizer(settings.GetString("calcSurfaceTerms"), sep);
    std::copy(tokenizer.begin(), tokenizer.end(),
              std::back_inserter(surfaceRegions));

    //check for all and none keyword
    bool checkNONE = (std::find(surfaceRegions.begin(),surfaceRegions.end(),"none") != surfaceRegions.end());
    bool checkALL = (std::find(surfaceRegions.begin(),surfaceRegions.end(),"all") != surfaceRegions.end());
    bool doIt = false;
    if(checkNONE && checkALL){
      Exception("You supplied the all and the none keyword for surface term. confused and aborting....");
    }else if(checkNONE){
      doIt = false;
    }else{
      doIt = true;
    }

    for (UInt actRegion = 0; actRegion<numRegions; actRegion++){
      //first we determine if this region is volume or surface by
      //checking the dimension of its first element
      //WARNING this may not be true but right now we have no other choice
      Elem::FEType elemType = (Elem::FEType) elemTypes_[regionElems_[actRegion][0]];
      if(Elem::GetElemDim(elemType) == dim_-1){
        //get region name and check if its contained in user list
        surfaceRegionIndices_[actRegion] = false;
      }else{
        volumeRegionIndices_.push_back(actRegion);
      }
    }
    if(!doIt)
      return;


    std::cout << "========================================"
              << "========================================"
              << std::endl;
    std::cout << "                        "
              << "Calculating Surface Neighbor information" << std::endl;
    std::cout << "========================================"
              << "========================================"
              << std::endl;
    std::cout.flush();

    //now we have the information so we iteratre over the surface regions
    //and fill the structures
    std::map<UInt,bool>::iterator it =  surfaceRegionIndices_.begin();
    while(it != surfaceRegionIndices_.end()){
      std::string regionName = ptFileReader_->GetRegionName(it->first);
      bool curregion = (std::find(surfaceRegions.begin(),surfaceRegions.end(),regionName) != surfaceRegions.end());
      if(curregion || checkALL){
        it->second = true;
        GetNeighborElements(surfaceNeighbors_,it->first);
      }
      it++;
    }
    std::cout << "========================================"
              << "========================================"
              << std::endl;
    std::cout << "                        "
              << "Leaving Surface Neighbor calculation" << std::endl;
    std::cout << "========================================"
              << "========================================"
              << std::endl;
    std::cout.flush();
  }


  bool CouplingHandler::CalcSurfaceNormalOutVolume(UInt globalVolElemNum, UInt globalSurfaceElemNum){
    std::vector<UInt> surfaceNodes,volumeNodes;
    std::set<UInt> surfaceSet,volumeSet;
    std::set<UInt> intersectSet;
    std::set<UInt> diffSet;

    UInt maxNENodes = ptFileReader_->GetMaxNumElemNodes();
    Elem::FEType surfaceElemType = (Elem::FEType) elemTypes_[globalSurfaceElemNum];
    Elem::FEType volumeElemType = (Elem::FEType) elemTypes_[globalVolElemNum];

    UInt numVolElemNodes = Elem::GetNumElemNodes(volumeElemType);
    UInt numSurfElemNodes = Elem::GetNumElemNodes(surfaceElemType);
    surfaceNodes.resize(numSurfElemNodes);
    volumeNodes.resize(numVolElemNodes);


    for(UInt n=0; n<numSurfElemNodes; n++){
      surfaceNodes[n] = topology_[globalSurfaceElemNum * maxNENodes + n];
    }
    for(UInt n=0; n<numVolElemNodes; n++){
      volumeNodes[n] = topology_[globalVolElemNum * maxNENodes + n];
    }
    surfaceSet = std::set<UInt>(surfaceNodes.begin(), surfaceNodes.end());
    volumeSet = std::set<UInt>(volumeNodes.begin(), volumeNodes.end());

    std::set_intersection(surfaceSet.begin(), surfaceSet.end(),
                          volumeSet.begin(), volumeSet.end(),
                          std::inserter(intersectSet,intersectSet.begin()));

    if(intersectSet.size() == surfaceNodes.size()){
       //ok, now we determine the normal along with its sign and we store it
      Vector<Double> & curNormal = surfaceNormals_[globalSurfaceElemNum];
      CalcSurfaceNormal(surfaceNodes,curNormal);
      std::set_difference( volumeSet.begin(), volumeSet.end(),
                           surfaceSet.begin(), surfaceSet.end(),
                           std::inserter(diffSet,diffSet.begin()) );

      Vector<Double> diffVec( dim_ );
      Double scalarProd = 0.0;
      std::set<UInt>::const_iterator it = diffSet.begin();
      UInt coordIdxCommon = ((*intersectSet.begin()) -1)*3;
      const Double EPS = 1e-30;
      for( ; it != diffSet.end(); it++ ) {

        UInt node = *it;

        // calculate difference vector (pointing out of the volume)
        for( UInt iDim = 0; iDim < dim_; ++iDim ) {
          UInt coordIdxVol = (node -1)*3;

          Double curCoordVol = nodalCoords_[coordIdxVol+iDim];
          Double curCoordSurf = nodalCoords_[coordIdxCommon+iDim];

          diffVec[iDim] =  curCoordSurf - curCoordVol;
        }
        // normalize difference vector to 1.0
        diffVec /= diffVec.NormL2();

        //scalar product
        scalarProd = diffVec * curNormal;

        if( std::abs(scalarProd) < EPS ) {
          // we have to find another node
          continue;
        } else {
          if( scalarProd < 0.0 ) {
            curNormal *= -1.0;
          }
          break;
        }
      }

      return true;
    }else{
      return false;
    }
  }

  void CouplingHandler::CalcSurfaceNormal( const std::vector<UInt>& nodeNums,
                                 Vector<Double>& nVec){
    UInt numNodes = nodeNums.size();
    UInt coordIdx;
    Vector<Double> vec1(dim_);
    Vector<Double> vec2(dim_);
    nVec.Resize(dim_);
    nVec.Init();
    //create matrix of coordinates
    Matrix<Double> ptCoord(dim_,numNodes);
    for(UInt i=0;i<numNodes;i++){
      coordIdx = (nodeNums[i] -1)*3;
      for(UInt d=0;d<dim_;d++){
        ptCoord[d][i] = nodalCoords_[coordIdx+d];
      }
    }

    if ( dim_ == 3) {
      UInt idx=0;
      if ( ptCoord.GetNumCols() == 4) {
        idx = 3;
      }
      else if ( ptCoord.GetNumCols() == 3) {
        idx = 2;
      }
      else {
        EXCEPTION("Surface Element is no quadrilateral or triangle even though we have 3D");
      }

      for ( UInt i=0; i<dim_; i++ ) {
        vec1[i] = ptCoord[i][1]   - ptCoord[i][0];
        vec2[i] = ptCoord[i][idx] - ptCoord[i][0];
      }

      // cross product: vec1 x vec 2
      nVec[0] = vec1[1]*vec2[2] - vec1[2]*vec2[1];
      nVec[1] = vec1[2]*vec2[0] - vec1[0]*vec2[2];
      nVec[2] = vec1[0]*vec2[1] - vec1[1]*vec2[0];
    }
    else {
      nVec[0] =  ptCoord[1][1] - ptCoord[1][0];
      nVec[1] = -ptCoord[0][1] + ptCoord[0][0];
    }


    Double invLength = 1.0/nVec.NormL2();
    nVec *= invLength;

  }
} // end of namespace
