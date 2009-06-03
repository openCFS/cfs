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

#include <boost/tokenizer.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/convenience.hpp>
#include <boost/filesystem/exception.hpp>
#include <boost/algorithm/string/trim.hpp>
namespace fs=boost::filesystem;

#include "General/exception.hh"
#include "Settings.hh"
#include "General/environment.hh"
#include "integlib/elemIntegr.hh"

#include "FlowDataTypes.hh"
#include "OutputWriter.hh"
#include "CouplingHandler.hh"

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
    std::vector<UInt> regionNodes;
    UInt maxNumElemNodes = 0;
    std::vector<std::string> regionNames;

    numRegionNodes_.resize(numRegions);

    // First read everything into internal buffers
    ptFileReader_->ReadNodalCoords(nodalCoords_);
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
    ptFileReader_->ReadTopology(topology_,
                                elemTypes_);

    // Determine the maximum number of element nodes
    maxNumElemNodes = ptFileReader_->GetMaxNumElemNodes();
    std::set<UInt> topoSet(topology_.begin(), topology_.end());
    if(*topoSet.begin() == 0)
      topoSet.erase(topoSet.begin());

    std::set<UInt>::iterator topoIt, topoEnd;
    topoIt = topoSet.begin();
    topoEnd = topoSet.end();
    
    std::map<UInt, UInt> pointMap;
    for( UInt i=0; topoIt != topoEnd; topoIt++, i++ ) 
    {
      pointMap[*topoIt] = i+1;
      std::cout << (*topoIt) << " -> " << (pointMap[*topoIt]) << std::endl;
      
      UInt idxNew=i*3;
      UInt idxOld=(*topoIt-1)*3;
      nodalCoords_[idxNew+0] = nodalCoords_[idxOld+0];
      nodalCoords_[idxNew+1] = nodalCoords_[idxOld+1];
      nodalCoords_[idxNew+2] = nodalCoords_[idxOld+2];
    }
    nodalCoords_.resize(topoSet.size()*3);
    

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
    std::map<RegionIdType, UInt > regionDims;

    for (int actRegion=0; actRegion<numRegions; actRegion++)
    {
      // Fill vector with region names
      regionNames.push_back(ptFileReader_->GetRegionName(actRegion));

      // Clear temporary containers
      regionNodeSet.clear();
      regionNodes.clear();

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
        UInt elemNum = *elemIt;
        UInt regionDim = regionDims[actRegion];
        Elem::FEType elemType = (Elem::FEType) elemTypes_[elemNum-1];
        
        if(!regionDim) 
        {
          regionDims[actRegion] = Elem::GetElemDim(elemType);
        }
        else
        {
          if( regionDim != Elem::GetElemDim(elemType) )
          {
            EXCEPTION("Elements with different dimensions have been "
                      << "encountered in region '" << (*regionNames.rbegin()) << "'!\n"
                      << "The error occured while examining element "
                      << (*elemIt) << ".\n"
                      << "Please check your mesh file!");
          }    
        }
      }
      
      // Put all nodes in a partition into a set to get an ordered list
      UInt idx = 0;
      UInt regionDim = 0;

      CollectElementNodes(regionElems_[actRegion], regionNodes, regionDim);

      // Write region elements and nodes
      owIt = outputWriters_.begin();
      owEnd = outputWriters_.end();
      
      for( ; owIt != owEnd; owIt++) 
      {
        (*owIt)->WriteRegion(regionNames[actRegion],
                             regionElems_[actRegion],
                             regionNodes,
                             regionDim);
      }
      
      it = regionNodes.begin();
      end = regionNodes.end();

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
      
      while ( ( counter < numFiles ) && readOK)
      {
        stepVal = ptFileReader_->GetTimeStep(counter);
        stepNum = counter + 1;
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
          // scale the nodal values
          std::stringstream velstr;
          velstr << settings.GetString("velscale");
          Double velScaleX, velScaleY, velScaleZ;
          velstr >> velScaleX >> velScaleY >> velScaleZ;
          if (velScaleX != 1.0)
          {
            velScale_(flowData, velScaleX, 0);
          }
          if (velScaleY != 1.0)
          {
            velScale_(flowData, velScaleY, 1);
          }
          if (velScaleZ != 1.0)
          {
            velScale_(flowData, velScaleZ, 2);
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
                    << "Exiting read time values loop...";
          
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
          
          // If the user requests the calculation of the Lighthill
          // source term, follow his order!
          if(calcSrc)
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
        counter++;
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

    int nElems = ptFileReader_->GetNumElems(regionIdx);

    Matrix<Double> coordMat;
    Matrix<Double> nodaldTijdxj;
    Matrix<Double> nodalVel;
    Vector<Double> elemVec;

    Elem::FEType elemType;
    UInt numElemNodes;
    UInt elemDim;
    UInt elemIdx;
    UInt maxNENodes = ptFileReader_->GetMaxNumElemNodes();
    UInt nodeNum;

    for( int i=0; i<nElems; i++)
    {
      elemIdx = regionElems_[regionIdx][i] - 1;
      elemType = (Elem::FEType) elemTypes_[elemIdx];
      numElemNodes = Elem::GetNumElemNodes(elemType);
      elemDim = Elem::GetElemDim(elemType);

      // Just calculate sources for volume elements!
      if(elemDim < dim_)
        continue;

      coordMat.Resize(elemDim, numElemNodes);
      nodaldTijdxj.Resize(elemDim, numElemNodes);
      nodalVel.Resize(elemDim, numElemNodes);

      for( UInt n=0; n<numElemNodes; n++)
      {
        nodeNum = topology_[elemIdx * maxNENodes + n];
        UInt topoIdx = (nodeNum - 1) * 3;
        UInt velIdx = regionNodeIndices_[regionIdx][nodeNum] * dim_;

        for( UInt d=0; d<elemDim; d++)
        {
          coordMat[d][n] = nodalCoords_[topoIdx+d];
          nodalVel[d][n] = velField[velIdx+d];
        }
      }

      try {
        ptElemIntegr_[elemType]->PerformIntegration( coordMat, nodaldTijdxj,
                                                     nodalVel, elemVec, density);
      } catch (CoupledField::Exception &ex)
      {
        std::cerr << "Warning: An Exception occurred during source term "
                  << "computation:\nElement " << i+1 << " of partition "
                  << ptFileReader_->GetRegionName(regionIdx) << std::endl;

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
        UInt idx;
        nodeNum = topology_[elemIdx * maxNENodes + n];
        idx = regionNodeIndices_[regionIdx][nodeNum];

#ifndef NDEBUG
        if (std::isnan(elemVec[n]) || std::isinf(elemVec[n])) {
          EXCEPTION("Source term calculated on element " << i+1
                    << " is Inf or Nan.");
        }
#endif

        acouRhsField[idx] -= elemVec[n];
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

    if(numInputNodes == numNodes)
    {
      output = input;
      return;
    }

    it = regionNodeIndices_[partitionIdx].begin();
    end = regionNodeIndices_[partitionIdx].end();
    output.resize(numNodes * numDOFs);

    for( ; it != end; it++ )
    {
      idxInput = (it->first - 1) * numDOFs;
      idxOutput = it->second * numDOFs;

      for(UInt dof=0; dof < numDOFs; dof++ )
        output[idxOutput+dof] = input[idxInput+dof];
    }
  }

} // end of namespace
