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
#include <boost/progress.hpp>

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
              << "The error occurred while examining element "
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

      //check if we need acouLambVec, acouMixedMassLoad or explicitly request the
      //mean flow field
      if(std::find(outputFields_.begin(),outputFields_.end(),"acouLambRhs")       != outputFields_.end() ||
         std::find(outputFields_.begin(),outputFields_.end(),"acouMixedMassLoad") != outputFields_.end() ||
         std::find(outputFields_.begin(),outputFields_.end(),"aeroAcouSourceRhs") != outputFields_.end() ||
         std::find(outputFields_.begin(),outputFields_.end(),"all")               != outputFields_.end() ||
         settings.GetInt("calcMeanPresField") ||
         settings.GetInt("calcMeanVelField")){
        ComputeMeanValues(flowData);
      }

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
          const std::string physFieldScale_str[] = {"velscale","geomscale","presscale"};
          const SolutionType solType[] = {FLUIDMECH_VELOCITY,MECH_DISPLACEMENT,FLUIDMECH_PRESSURE};
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
          //HARD CODED FOR FLUIDMECH PRESSURE NEEDS TO BE CLEANER
          std::stringstream scaleStr;
          scaleStr << settings.GetString(physFieldScale_str[2]);
          scaleStr >> scaleX;
          if (!scaleStr.fail() && scaleX != 1.0){
            scale_PhysField_(flowData, solType[2], scaleX, 0);
          }

          // <-- end scaling velocity
          
          // Override the setting of --outprec for CFX
          if ( settings.GetString("type") == "CFX"
               && settings.GetString("cfxSinglePrecision") != "auto"
               && settings.GetInt("cfxSinglePrecision") )
          {
            settings.SetString("outprec", "single");
          }
          
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

          //check for surface tasks
          if(calcSrc &&  regionDims_[actRegion] == dim_-1){
            if ( flowData[actRegion].find(FLUIDMECH_FORCE)
                        != flowData[actRegion].end() )
             {
                flowData[actRegion][FLUIDMECH_FORCE].isActive = true;
             }else{
               flowData[actRegion][FLUIDMECH_FORCE].isActive = false;
             }
            CalculateMechSurfaceForce(actRegion,flowData[actRegion]);
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
            if(flowData[actRegion].find(FLUIDMECH_PRESSURE)
                 != flowData[actRegion].end() ){
              flowData[actRegion][FLUIDMECH_PRESSURE].isActive = true;
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
        // NOTE THIS: we check if there is any aeroacoustic result definied
        //    This is not that much of a hack as the situtation in which we
        //    want Lamb vector on region 1 and lighthill on region two is
        //    not very common
        if (doCalcMultiNodes)
        {
          std::map<std::string, std::vector<UInt>* > regionNodesActive;
          std::map<std::string, std::vector<UInt> >::iterator iterRegionAct = regionNodes_.begin();
          for (; iterRegionAct != regionNodes_.end(); ++iterRegionAct)
          {
            UInt regIdx = 0;
            while (regionNames[regIdx] != iterRegionAct->first ){
              ++regIdx;
            }

            if ( (flowData[regIdx].find(ACOU_RHS_LOAD) != flowData[regIdx].end()) ||
                 (flowData[regIdx].find(ACOUMIXED_MASS_LOAD) != flowData[regIdx].end()) ||
                 (flowData[regIdx].find(ACOU_LAMB_RHS) != flowData[regIdx].end()) )
            {
              if (flowData[regIdx][ACOU_RHS_LOAD].isActive ||
                  flowData[regIdx][ACOUMIXED_MASS_LOAD].isActive ||
                  flowData[regIdx][ACOU_LAMB_RHS].isActive)
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
            double accumValNodesMass = 0;
            double accumAerAcouSrc = 0;
            std::vector<double> accumValLambNodes(dim_,0.0);

            std::map<std::string, UInt>::const_iterator iterRegions = iterMultiNodes->second.begin();
            for (; iterRegions != iterMultiNodes->second.end(); ++iterRegions)
            {
              const std::string& regName = iterRegions->first;
              const UInt& node = iterRegions->second;
              UInt regIdx = 0;
              while (regionNames[regIdx] != regName ){
                ++regIdx;
              }

              if ( flowData[regIdx].find(ACOU_RHS_LOAD)
                  != flowData[regIdx].end() )
              {
                if (flowData[regIdx][ACOU_RHS_LOAD].isActive)
                  accumValNodes += flowData[regIdx][ACOU_RHS_LOAD].data[node];
              }

              if ( flowData[regIdx].find(ACOUMIXED_MASS_LOAD)
                  != flowData[regIdx].end() )
              {
                if (flowData[regIdx][ACOUMIXED_MASS_LOAD].isActive)
                  accumValNodesMass += flowData[regIdx][ACOUMIXED_MASS_LOAD].data[node];
              }

              if ( flowData[regIdx].find(AERO_ACOU_SRC_RHS)
                  != flowData[regIdx].end() )
              {
                if (flowData[regIdx][AERO_ACOU_SRC_RHS].isActive)
                  accumAerAcouSrc += flowData[regIdx][AERO_ACOU_SRC_RHS].data[node];
              }

              if ( flowData[regIdx].find(ACOU_LAMB_RHS)
                  != flowData[regIdx].end() )
              {
                if (flowData[regIdx][ACOU_LAMB_RHS].isActive){
                  for(UInt d = 0; d < dim_ ; d++){
                    accumValLambNodes[d] += flowData[regIdx][ACOU_LAMB_RHS].data[(node*dim_)+d];
                  }
                }
              }
            }
            // set accumulated values
            iterRegions = iterMultiNodes->second.begin();
            for (; iterRegions != iterMultiNodes->second.end(); ++iterRegions)
            {
              const std::string& regName = iterRegions->first;
              const UInt& node = iterRegions->second;
              UInt regIdx = 0;
              while (regionNames[regIdx] != regName ){
                  ++regIdx;
              }
              if ( flowData[regIdx].find(ACOU_RHS_LOAD)
                  != flowData[regIdx].end() )
              {
                if (flowData[regIdx][ACOU_RHS_LOAD].isActive)
                  flowData[regIdx][ACOU_RHS_LOAD].data[node] = accumValNodes;
              }

              if ( flowData[regIdx].find(ACOUMIXED_MASS_LOAD)
                  != flowData[regIdx].end() )
              {
                if (flowData[regIdx][ACOUMIXED_MASS_LOAD].isActive)
                  flowData[regIdx][ACOUMIXED_MASS_LOAD].data[node] = accumValNodesMass;
              }

              if ( flowData[regIdx].find(AERO_ACOU_SRC_RHS)
                  != flowData[regIdx].end() )
              {
                if (flowData[regIdx][AERO_ACOU_SRC_RHS].isActive)
                  flowData[regIdx][AERO_ACOU_SRC_RHS].data[node] = accumAerAcouSrc;
              }

              if ( flowData[regIdx].find(ACOU_LAMB_RHS)
                  != flowData[regIdx].end() )
              {
                if (flowData[regIdx][ACOU_LAMB_RHS].isActive){
                  for(UInt d = 0; d < dim_ ; d++){
                    flowData[regIdx][ACOU_LAMB_RHS].data[(node*dim_)+d] = accumValLambNodes[d];
                  }
                }

              }
            }
          }
        }// end of nodes on multiple region correction
        //deactive mean field for resutls file
        //if(counter > 3){
        //  for (actRegion = 0; actRegion<numRegions; actRegion++){
        //    flowData[actRegion][MEAN_FLUIDMECH_PRESSURE].isActive = false;
        //    flowData[actRegion][MEAN_FLUIDMECH_VELOCITY].isActive = false;
        //  }
        //}
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
        //reactiiveate for next step
        //if(counter > 3){
        //  for (actRegion = 0; actRegion<numRegions; actRegion++){
        //    flowData[actRegion][MEAN_FLUIDMECH_PRESSURE].isActive = true;
        //    flowData[actRegion][MEAN_FLUIDMECH_VELOCITY].isActive = true;
        //  }
        //}
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
  
  void CouplingHandler::ComputeMeanValues(std::vector<FlowDataType> & flowData){
    //ok we read in each timestep and store the values
    Settings& settings = Settings::Instance();
    Double stepVal = 0;

    UInt counter = 0;
    UInt numFiles = ptFileReader_->GetNumFiles();

    UInt userNumSteps = settings.GetInt("numStepsForMeanField");
    if(userNumSteps == 0){
      userNumSteps = numFiles;
    }
    UInt stepInc = settings.GetInt("stepincr");
    UInt numRegions = ptFileReader_->GetNumRegions();
    Double oldStepVal = 0;
    Double actDt = 0.0;
    Double startTime = 0;
    Double endTime = 0;
    std::map<UInt, std::vector<Double> > tmpMeanVelField;
    std::map<UInt, std::vector<Double> > tmpMeanPresField;

    bool velFieldRequest = false;
    bool presFieldRequest = false;

    bool meanPresRequestedByUser = settings.GetInt("calcMeanPresField");
    bool meanVelRequestedByUser = settings.GetInt("calcMeanVelField");

    if(std::find(outputFields_.begin(),outputFields_.end(),"acouLambRhs") != outputFields_.end() ||
       std::find(outputFields_.begin(),outputFields_.end(),"all") != outputFields_.end() ||
       std::find(outputFields_.begin(),outputFields_.end(),"aeroAcouSourceRhs") != outputFields_.end() ||
       meanVelRequestedByUser ){
      velFieldRequest = true;
    }

    if(std::find(outputFields_.begin(),outputFields_.end(),"acouMixedMassLoad") != outputFields_.end() ||
       std::find(outputFields_.begin(),outputFields_.end(),"aeroAcouSourceRhs") != outputFields_.end() ||
       std::find(outputFields_.begin(),outputFields_.end(),"all") != outputFields_.end() ||
       meanPresRequestedByUser ){
      presFieldRequest = true;
    }

    bool readOK = true;

    std::cout << "========================================"
              << "========================================"
              << std::endl;
    std::cout << "                        "
              << "Computing Mean Flow Fields" << std::endl;
    std::cout << "                        ";
    std::cout << "Number of transient files: " << userNumSteps << std::endl;
    std::cout << "========================================"
              << "========================================"
              << std::endl;

    boost::progress_display progress( userNumSteps );

    actDt = settings.GetDouble("timestep");
    while ( ( counter < userNumSteps*stepInc ) && readOK)
    {
      stepVal = ptFileReader_->GetTimeStep(counter);
      if(counter == 0)
        startTime = stepVal;

      if(counter == (userNumSteps-1)*stepInc){
        endTime = stepVal;
      }

      //this is now conceptually difficult
      //but user input is more significant...
      // so we try to get the time step from command line
      // if not provided there we compute it ourselves
      if(counter == 0){
        actDt = 0;
        oldStepVal = stepVal;
      }else if (settings.GetDouble("timestep") != -1.0){
        actDt = settings.GetDouble("timestep");
      }else{
        actDt = stepVal - oldStepVal;
        oldStepVal = stepVal;
      }
      //check for safty reasons for the timestep and throw error if -1
      if(actDt == -1.0){
        EXCEPTION("Error in computing mean field values: \n The timestep was calculated to be -1.0. Please make sure to specify it via the --timestep option.");
      }

      try{
       ptFileReader_->ReadNodalValues(flowData, activeParts_, counter);
       if (reduce_elementOrder_){
         ptFileReader_->ReduceOrderOfNodalValues(flowData, regionNodes_);
       }
      // Override the setting of --outprec for CFX
       if ( settings.GetString("type") == "CFX"
            && settings.GetString("cfxSinglePrecision") != "auto"
            && settings.GetInt("cfxSinglePrecision") )
       {
         settings.SetString("outprec", "single");
       }

      readOK = true;
      } catch (std::exception& ex)
      {
        std::cerr << "CAUGHT EXCEPTION while trying to read nodal values:"
                  << std::endl << ex.what() << std::endl
                  << "Exiting read time values loop..." << std::endl;

        readOK = false;
        continue;
      }
      for (UInt actRegion = 0; actRegion<numRegions && readOK; actRegion++)
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

        if(regionDims_[actRegion] == dim_){
          bool velField = true,presField = true;

          if(flowData[actRegion].find(FLUIDMECH_VELOCITY) == flowData[actRegion].end())
          {
           // std::string regionName = ptFileReader_->GetRegionName(actRegion);
           // std::cerr << "Cannot calculate mean velocity field on " << regionName
           //           << " since no velocity field is available!" << std::endl;
            velField = false;
          }
          if(flowData[actRegion].find(FLUIDMECH_PRESSURE) == flowData[actRegion].end()){
           // std::string regionName = ptFileReader_->GetRegionName(actRegion);
           // std::cerr << "Cannot calculate mean pressure field on " << regionName
           //           << " since no pressure field is available!" << std::endl;
            presField = false;
          }

          if(velField && velFieldRequest){
            FlowDataPartStruct& fdps = flowData[actRegion][FLUIDMECH_VELOCITY];
            const std::vector<Double>& velocityFieldVector = fdps.data;

            if(tmpMeanVelField.find(actRegion) == tmpMeanVelField.end()) {
              tmpMeanVelField[actRegion].resize(velocityFieldVector.size());
              std::fill(tmpMeanVelField[actRegion].begin(), tmpMeanVelField[actRegion].end(), 0);
            }

            UInt numEs = velocityFieldVector.size();
#pragma omp parallel for
            for(UInt i = 0; i<numEs;++i){
              tmpMeanVelField[actRegion][i] += velocityFieldVector[i] * actDt;
            }
          }

          if(presField && presFieldRequest){
            FlowDataPartStruct& fdps = flowData[actRegion][FLUIDMECH_PRESSURE];
            const std::vector<Double>& presureFieldVector = fdps.data;

            if(tmpMeanPresField.find(actRegion) == tmpMeanPresField.end()) {
              tmpMeanPresField[actRegion].resize(presureFieldVector.size());
              std::fill(tmpMeanPresField[actRegion].begin(), tmpMeanPresField[actRegion].end(), 0);
            }
            UInt numEs = presureFieldVector.size();
#pragma omp parallel for
            for(UInt i = 0; i<numEs;++i){
              tmpMeanPresField[actRegion][i] += presureFieldVector[i] * actDt;
            }
          }
        }
      }

      counter = counter + stepInc;
      ++progress;
    }

    Double simTime = endTime - startTime;

    std::cout << "   Averaging quantities ...";
    std::cout.flush();
    for (UInt actRegion = 0; actRegion<numRegions && readOK; actRegion++){
      if(regionDims_[actRegion] == dim_ && tmpMeanPresField.find(actRegion) != tmpMeanPresField.end()){
        FlowDataPartStruct& fdps1 = flowData[actRegion][MEAN_FLUIDMECH_PRESSURE];
        fdps1.isActive = true; // all partitions have results
        fdps1.definedOn = ResultInfo::NODE; // nodes
        if(fdps1.dofNames.empty())
          fdps1.dofNames.push_back("-");
        fdps1.unit = MapSolTypeToUnit(MEAN_FLUIDMECH_PRESSURE);
        fdps1.resultName = "meanFluidMechPressure";

        fdps1.entryType = ResultInfo::SCALAR;
        fdps1.isActive = true;

        fdps1.data.resize(tmpMeanPresField[actRegion].size());
        std::fill(fdps1.data.begin(), fdps1.data.end(), 0);

        UInt numE = tmpMeanPresField[actRegion].size();
#pragma omp parallel for
        for(UInt i = 0; i<numE;++i){
          fdps1.data[i] = tmpMeanPresField[actRegion][i] / (simTime);
        }
      }

      if(regionDims_[actRegion] == dim_ && tmpMeanVelField.find(actRegion) != tmpMeanVelField.end()){
        FlowDataPartStruct& fdps2 = flowData[actRegion][MEAN_FLUIDMECH_VELOCITY];
        fdps2.isActive = true; // all partitions have results
        fdps2.definedOn = ResultInfo::NODE; // nodes
        if(fdps2.dofNames.empty()) {
          fdps2.dofNames.push_back("x");
          fdps2.dofNames.push_back("y");
          if(dim_ == 3)
            fdps2.dofNames.push_back("z");
        }
        fdps2.unit = MapSolTypeToUnit(MEAN_FLUIDMECH_VELOCITY);
        fdps2.resultName = "meanFluidMechVelocity";
        fdps2.entryType = ResultInfo::VECTOR;
        fdps2.isActive = true;

        fdps2.data.resize(tmpMeanPresField[actRegion].size()*dim_);

        std::fill(fdps2.data.begin(), fdps2.data.end(), 0);

#pragma omp parallel for
        for(UInt i = 0; i<tmpMeanVelField[actRegion].size();++i){
          fdps2.data[i] = tmpMeanVelField[actRegion][i] / (simTime);
        }
      }
    }
    std::cout << "done" << std::endl;
    std::cout.flush();


  }

  void CouplingHandler::Finish(){
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

    //OK, we have so many different source formulations lets determine what the user wants and
    //what is available

    //lets see what is wanted
    bool computeWaveRhs = ( std::find(outputFields_.begin(),outputFields_.end(),"acouRhsLoad") != outputFields_.end() ||
                       std::find(outputFields_.begin(),outputFields_.end(),"acouDivLighthillTensor") != outputFields_.end() ||
                       std::find(outputFields_.begin(),outputFields_.end(),"acouRhsLoadDensity") != outputFields_.end() ||
                       std::find(outputFields_.begin(),outputFields_.end(),"all") != outputFields_.end());

    //prepare some booleans to determine which quantity should be used for the Wave equation
    bool useDivLHT = false;
    bool computeLHP = false;
    bool computeLHV = false;
    bool usePresD2 = false;

    if(computeWaveRhs){
      //now we determine what to do depending on the quantity the user specified
      //get the solution type string
      std::string lhSrcQuant = settings.GetString("quantityForAcouRhsLoad");
      //parse to solution type
      SolutionType lhType = SolutionTypeEnum.Parse(lhSrcQuant);
      //now we set one of our darn flags
      //BIG TODO: This needs to be more failsafe, the big variety of different formulations will require
      // a refactoring of cplreader...
      switch(lhType){
      case FLUIDMECH_VELOCITY:
        if(flowData.find(FLUIDMECH_VELOCITY)== flowData.end()){
          EXCEPTION("Trying to compute wave equation RHS with nodal velocity, but no this field was found in input data!")
        }
        computeLHV = true;
        break;
      case FLUIDMECH_PRESSURE:
        if(flowData.find(FLUIDMECH_PRESSURE)== flowData.end()){
          EXCEPTION("Trying to compute wave equation RHS with nodal pressure, but no this field was found in input data!")
        }
        computeLHP = true;
        break;
      case FLUIDMECH_DIV_LH_T:
        if(flowData.find(ACOU_DIV_LH_TENSOR_NODAL)== flowData.end()){
          EXCEPTION("Trying to compute wave equation RHS with nodal divergence of LH tensor, but this field was found in input data!")
        }
        useDivLHT = true;
        break;
      case FLUIDMECH_PRESSURE_DERIV_2:
        if(flowData.find(FLUIDMECH_PRESSURE_DERIV_2)== flowData.end()){
          EXCEPTION("Trying to compute wave equation RHS with nodal Laplacian of pressure, but this field was found in input data!")
        }
        usePresD2 = true;
        break;
      default:
        EXCEPTION("Specified a quantity for computing Wave equation RHS which is not supported");
        break;
      }
    }


    //Here we come to the perturbation equation stuff, lets see what our user wants
    bool computeAPEMass = ( std::find(outputFields_.begin(),outputFields_.end(),"acouMixedMassLoad") != outputFields_.end() ||
                            std::find(outputFields_.begin(),outputFields_.end(),"all") != outputFields_.end());

    bool computeAPEMomentum = ( std::find(outputFields_.begin(),outputFields_.end(),"acouLambRhs") != outputFields_.end() ||
                              std::find(outputFields_.begin(),outputFields_.end(),"acouLambVec") != outputFields_.end() ||
                              std::find(outputFields_.begin(),outputFields_.end(),"all") != outputFields_.end());

    bool computeAeroAcouSrc = ( std::find(outputFields_.begin(),outputFields_.end(),"aeroAcouSourceRhs") != outputFields_.end() ||
                                std::find(outputFields_.begin(),outputFields_.end(),"all") != outputFields_.end());

    //security check for PE
    if(computeAPEMass){
      if(flowData.find(FLUIDMECH_PRESSURE) == flowData.end() || flowData.find(FLUIDMECH_VELOCITY) == flowData.end() ){
        EXCEPTION("Cannot compute acouMixedMassLoad since pressure or velocity field is not available");
      }
    }
    if(computeAPEMomentum || computeAeroAcouSrc){
      if(flowData.find(FLUIDMECH_VELOCITY) == flowData.end()){
        EXCEPTION("Cannot compute acouMixedMomentumLoad or aeroAcouSrcRhs since velocity field is not available");
      }
    }

    
    //obtain the fields provided by the file reader.
    //If there would be a quantity which is requested but not provided, the checks above should already break execution
    // first we check for the velocity field
    FlowDataPartStruct& velocityStruct = flowData[FLUIDMECH_VELOCITY];
    std::vector<Double>& velField = velocityStruct.data;

    FlowDataPartStruct& meanVelocityStruct = flowData[MEAN_FLUIDMECH_VELOCITY];
    std::vector<Double>& meanVelocityField = meanVelocityStruct.data;
    
    // then we check for the divergence of Lighthill tensor field
    FlowDataPartStruct& divlhtStruct = flowData[ACOU_DIV_LH_TENSOR_NODAL];
    std::vector<Double>& divlhtField = divlhtStruct.data;

    // lets check for second derivative of fluid pressure
    FlowDataPartStruct& presD2Struct = flowData[FLUIDMECH_PRESSURE_DERIV_2];
    std::vector<Double>& presD2Field = presD2Struct.data;


    //THESE CHECKS SHOULD BE OBSOLETE NOW as inactive quantities will not be accessed later
    //BUT BE WARNED: everytime we access the map, those entries will be available! Anyhow, the aboves checks should
    // work in the first timestep which should be enough for now. It is a dirty hack anyhow....

    //now we check if we get the right this to do what is requested
    //if(!useDivLHT) {
    //
    //  if(!velocityStruct.isActive)
    //  {
    //    if(computeLHV || computeAPEMomentum || computeAeroAcouSrc){
    //      std::cerr << "Will not calculate velocity based sources on " << regionName
    //                << " since velocity field is not active!" << std::endl;
    //    }
    //    flowData.erase(FLUIDMECH_VELOCITY);
    //    flowData.erase(MEAN_FLUIDMECH_VELOCITY);
    //    computeLHV = false;
    //    computeAPEMomentum = false;
    //    computeAeroAcouSrc = false;
    //  }
    //
    //} else {
    //
    //  if(!divlhtStruct.isActive)
    //  {
    //    if(computeLHV || computeAPEMomentum || computeAeroAcouSrc){
    //      std::cerr << "Will not calculate divLHT based sources on " << regionName
    //                << " since divLHT field is not active!" << std::endl;
    //    }
    //    flowData.erase(ACOU_DIV_LH_TENSOR_NODAL);
    //    computeLHV = false;
    //    computeAPEMomentum = false;
    //    computeAeroAcouSrc = false;
    //  }
    //
    //}


    // now we turn to pressure field
    bool presFieldAvailable = true;
    FlowDataPartStruct& meanPressureStruct = flowData[MEAN_FLUIDMECH_PRESSURE];
    std::vector<Double>& meanPressureField = meanPressureStruct.data;


    FlowDataPartStruct& pressureStruct = flowData[FLUIDMECH_PRESSURE];
    std::vector<Double>& pressureField = pressureStruct.data;

    if(!pressureStruct.isActive){
      if(computeAPEMass || computeLHP){
        std::cerr << "Cannot calculate Pressure sources on " << regionName
                  << " since no pressure field is available!" << std::endl;
      }
      presFieldAvailable = false;
      computeAPEMass = false;
      computeLHP = false;
      flowData.erase(FLUIDMECH_PRESSURE);
      flowData.erase(MEAN_FLUIDMECH_PRESSURE);
    }

    // ok now will double check if the mean fields are available
    bool meanVelFieldAvail = false;

    meanVelFieldAvail = (meanVelocityField.size()>0);
//    meanPresFieldAvail = (meanPressureField.size() > 0);

    //Well thats it, now we should be sure to set the flags right to only compute things the user wants and the flow data
    //features. The flags will be checked if needed. we can go on...


    std::cout << "Calculating Aeroacoustic sources on " << regionName << " ";
    Double density = settings.GetDouble("density");
    if(density == 1.0 && computeLHV){
      WARN("The density is set to 1.0 and Lighthill sources should be computed. Check if this is what you want...")
    }

    //before initializing all structures we check for perturbed quantities in case of APEMassLoad
    if( computeAPEMass ){
      ComputePerturbedPressureField(pressureField,regionIdx,flowData[MEAN_FLUIDMECH_PRESSURE]);
      UpdatePressureTimeDeriv(perturbedPressureField_[regionIdx] , regionIdx);
    }


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
    std::fill(acouRhsField.begin(), acouRhsField.end(), 0.0);

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
    std::fill(acouRhsDensityField.begin(), acouRhsDensityField.end(), 0.0);

    int nElems = ptFileReader_->GetNumElems(regionIdx);

    FlowDataPartStruct& fdps4 = flowData[ACOU_DIV_LH_TENSOR];
    fdps4.isActive = true; // all partitions have results
    fdps4.definedOn = ResultInfo::ELEMENT; // elements
    if (fdps4.dofNames.empty()) {
      fdps4.dofNames.push_back("x");
      fdps4.dofNames.push_back("y");
      if (dim_ == 3)
        fdps4.dofNames.push_back("z");
    }
    fdps4.unit = MapSolTypeToUnit(ACOU_DIV_LH_TENSOR);
    fdps4.resultName = "acouDivLighthillTensor";
    fdps4.data.resize(nElems * dim_);
    fdps4.entryType = ResultInfo::VECTOR;
    std::vector<Double>& acouDivLighthillTensor = fdps4.data;

    // Fill acouDivLighthillTensor field with zeros
    std::fill(acouDivLighthillTensor.begin(), acouDivLighthillTensor.end(), 0.0);



    // Init the structures for AcouMixedMassLoad
    FlowDataPartStruct& fdps5 = flowData[ACOUMIXED_MASS_LOAD];
    fdps5.isActive = true; // all partitions have results
    fdps5.definedOn = ResultInfo::NODE; // nodes
    if(fdps5.dofNames.empty())
      fdps5.dofNames.push_back("-");
    fdps5.unit = MapSolTypeToUnit(ACOUMIXED_MASS_LOAD);
    fdps5.resultName = "acouMixedMassLoad";
    fdps5.data.resize(numRegionNodes_[regionIdx]);
    fdps5.entryType = ResultInfo::SCALAR;
    std::vector<Double>& acouMixedMassRhsField = fdps5.data;

    // Fill acouRhsLoad field with zeros
    std::fill(acouMixedMassRhsField.begin(), acouMixedMassRhsField.end(), 0);

    // PREPARE THE LAMB VECTOR
    FlowDataPartStruct& fdps6 = flowData[ACOU_LAMB_VEC];
    fdps6.isActive = true; // all partitions have results
    fdps6.definedOn = ResultInfo::ELEMENT; // elements
    if(fdps6.dofNames.empty()) {
      fdps6.dofNames.push_back("x");
      fdps6.dofNames.push_back("y");
      if(dim_ == 3)
        fdps6.dofNames.push_back("z");
    }
    fdps6.unit = MapSolTypeToUnit(ACOU_LAMB_VEC);
    fdps6.resultName = "acouLambVec";
    fdps6.data.resize(nElems * dim_);
    fdps6.entryType = ResultInfo::VECTOR;
    std::vector<Double>& acouLambVec = fdps6.data;

    // Fill acouLambVector field with zeros
    std::fill(acouLambVec.begin(), acouLambVec.end(), 0);

    // PREPARE THE MOMENTUM RHS VECTOR i.e. based on Lamb Vector
    FlowDataPartStruct& fdps7 = flowData[ACOU_LAMB_RHS];
    fdps7.isActive = true; // all partitions have results
    fdps7.definedOn = ResultInfo::NODE; // elements
    if(fdps7.dofNames.empty()) {
      fdps7.dofNames.push_back("x");
      fdps7.dofNames.push_back("y");
      if(dim_ == 3)
        fdps7.dofNames.push_back("z");
    }
    fdps7.unit = MapSolTypeToUnit(ACOU_LAMB_RHS);
    fdps7.resultName = "acouLambRhs";
    fdps7.data.resize(numRegionNodes_[regionIdx] * dim_);
    fdps7.entryType = ResultInfo::VECTOR;
    std::vector<Double>& acouLambVecRhsField = fdps7.data;

    // Fill acouDivLighthillTensor field with zeros
    std::fill(acouLambVecRhsField.begin(), acouLambVecRhsField.end(), 0);

    // Init the structures for AeroAcousSource i.e. extraction of fludiMechPressure
    FlowDataPartStruct& fdps8 = flowData[AERO_ACOU_SRC_RHS];
    fdps8.isActive = true; // all partitions have results
    fdps8.definedOn = ResultInfo::NODE; // nodes
    if(fdps8.dofNames.empty())
      fdps8.dofNames.push_back("-");
    fdps8.unit = MapSolTypeToUnit(AERO_ACOU_SRC_RHS);
    fdps8.resultName = "aeroAcouSourceRhs";
    fdps8.data.resize(numRegionNodes_[regionIdx]);
    fdps8.entryType = ResultInfo::SCALAR;
    std::vector<Double>& aeroAcouRhsField = fdps8.data;

    // Fill acouRhsLoad field with zeros
    std::fill(aeroAcouRhsField.begin(), aeroAcouRhsField.end(), 0);


#ifdef _OPENMP
    //now we loop over all elements serial or parallel
    int tnum = 1;
#pragma omp parallel
{
    tnum = omp_get_num_threads();
}
    IntegrationMap ptElemI(ptElemIntegr_);
    std::cout << "... " << tnum << " threads, parallel loop over " << nElems << " elements..." ;
#else
    std::cout << "...serial loop over " << nElems << " elements..." ;
#endif
    std::cout.flush();

#pragma omp parallel firstprivate(ptElemI)
{
    UInt velIdx,topoIdx,idx;
    Matrix<Double> coordMat;
    //nodal values for pressure and velocity
    Matrix<Double> nodalVel;
    Vector<Double> nodalPressure;
    Vector<Double> nodalPressureD2;


    // LIGHTHILL Terms
    Matrix<Double> nodaldTijdxj;
    Vector<Double> elemVecLH;
    Vector<Double> nodalLoadDensity;
    Vector<Double> divLHTensor(dim_);

    //VortexTerms
    Matrix<Double> nodalMeanVel;
    Vector<Double> elemVecLamb;
    Vector<Double> elemVecLambRhs;

    //Pressure source Terms
    Vector<Double> nodalMeanPressure;
    Vector<Double> nodalPressureTDeriv;
    Vector<Double> nodalPerturbedPressure;
    Vector<Double> elemVecPres;


    //aeroacoustic Source
    Vector<Double> elemVecAeroAcou;


    Elem::FEType elemType;
    UInt numElemNodes;
    UInt elemDim;
    UInt elemIdx;
    UInt maxNENodes = ptFileReader_->GetMaxNumElemNodes();
    UInt nodeNum;

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
      nodalMeanVel.Resize(elemDim, numElemNodes);
      nodalMeanPressure.Resize(numElemNodes,0.0);
      nodalPressureTDeriv.Resize(numElemNodes,0.0);
      nodalPressure.Resize(numElemNodes,0.0);
      nodalPressureD2.Resize(numElemNodes,0.0);
      nodalPerturbedPressure.Resize(numElemNodes,0.0);

      for( UInt n=0; n<numElemNodes; n++)
      {
        nodeNum = topology_[elemIdx * maxNENodes + n];
        topoIdx = (nodeNum - 1) * 3;
        velIdx = regionNodeIndices_[regionIdx][nodeNum] * dim_;

        //for the PE mass equation load we need perturbed pressure as well as
        // its time derivative
        if(computeAPEMass){
          nodalMeanPressure[n] = meanPressureField[regionNodeIndices_[regionIdx][nodeNum]];
          nodalPressureTDeriv[n] =  pressureTimeDeriv_[regionIdx][regionNodeIndices_[regionIdx][nodeNum]];
          nodalPerturbedPressure[n] = perturbedPressureField_[regionIdx][regionNodeIndices_[regionIdx][nodeNum]];
        }

        //for pressure based wave equation we just need the fluid pressure
        if(computeLHP){
          nodalPressure[n] = pressureField[regionNodeIndices_[regionIdx][nodeNum]];
        }
        if(usePresD2){
          nodalPressureD2[n] = presD2Field[regionNodeIndices_[regionIdx][nodeNum]];
        }

        for( UInt d=0; d<elemDim; d++)
        {
          coordMat[d][n] = nodalCoords_[topoIdx+d];

          if(computeLHV || computeAPEMomentum || computeAeroAcouSrc)
          {
            if (!useDivLHT) {
              nodalVel[d][n] = velField[velIdx+d];
            } else {
              nodaldTijdxj[d][n] = divlhtField[velIdx+d];
            }
            if(meanVelFieldAvail)
              nodalMeanVel[d][n] = meanVelocityField[velIdx+d];
          }
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
          ptElemI[elemType].PerformIntegrationCentre( coordMat,
                                                      nodalVel,
                                                      nodalPressure,
                                                      nodalMeanVel,
                                                      nodalMeanPressure,
                                                      nodalPressureTDeriv,
                                                      elemVecLH,
                                                      nodalLoadDensity,
                                                      divLHTensor,
                                                      elemVecLamb,
                                                      elemVecLambRhs,
                                                      elemVecPres,
                                                      elemVecAeroAcou,
                    density);
          } else {
            if (computeWaveRhs) {

              if (computeLHV) {

                ptElemI[elemType].PerformIntegrationLighthill(coordMat,
                        nodaldTijdxj,
                        nodalVel,
                        elemVecLH,
                        nodalLoadDensity,
                        divLHTensor,
                        density);
              } else if(useDivLHT) {
                ptElemI[elemType].PerformIntegrationLighthillwithDivTij(coordMat,
                        nodaldTijdxj,
                        nodalVel,
                        elemVecLH,
                        nodalLoadDensity,
                        divLHTensor,
                        density);
              } else if(computeLHP) {
                 ptElemI[elemType].PerformIntegrationLHPressure(coordMat,
                         nodalPressure,
                         elemVecLH,
                         nodalLoadDensity);
              } else if(usePresD2){
                 ptElemI[elemType].PerformIntegrationPresD2(coordMat,
                                                            nodalPressureD2,
                                                            elemVecLH,
                                                            nodalLoadDensity);
              }

          }

          if(computeAPEMass){
            ptElemI[elemType].PerformIntegrationAPEMass(coordMat,
                                                        nodalVel,
                                                        nodalPerturbedPressure,
                                                        nodalPressureTDeriv,
                                                        elemVecPres,
                                                        density);
          }

          if(computeAPEMomentum){
            ptElemI[elemType].PerformIntegrationAPEMomentum(coordMat,
                                                            nodalVel,
                                                            nodalMeanVel,
                                                            elemVecLamb,
                                                            elemVecLambRhs,
                                                            density);

          }

          if(computeAeroAcouSrc){
            ptElemI[elemType].PerformIntegrationAeroAcouSrc(coordMat,
                                                            nodalVel,
                                                            nodalMeanVel,
                                                            elemVecAeroAcou);
          }
        }
#else
        if (doIntAverageCentre_)
        {
          ptElemIntegr_[elemType]->PerformIntegrationCentre( coordMat,
                                                             nodalVel,
                                                             nodalPressure,
                                                             nodalMeanVel,
                                                             nodalMeanPressure,
                                                             nodalPressureTDeriv,
                                                             elemVecLH,
                                                             nodalLoadDensity,
                                                             divLHTensor,
                                                             elemVecLamb,
                                                             elemVecLambRhs,
                                                             elemVecPres,
                                                             elemVecAeroAcou,
                                                             density);
          } else {
            if (computeWaveRhs) {

              if (computeLHV) {

                ptElemIntegr_[elemType]->PerformIntegrationLighthill(coordMat,
                        nodaldTijdxj,
                        nodalVel,
                        elemVecLH,
                        nodalLoadDensity,
                        divLHTensor,
                        density);
              } else if(useDivLHT) {
                ptElemIntegr_[elemType]->PerformIntegrationLighthillwithDivTij(coordMat,
                        nodaldTijdxj,
                        nodalVel,
                        elemVecLH,
                        nodalLoadDensity,
                        divLHTensor,
                        density);
              }else if(computeLHP) {
                ptElemIntegr_[elemType]->PerformIntegrationLHPressure(coordMat,
                                                                 nodalPressure,
                                                                 elemVecLH,
                                                                 nodalLoadDensity);
              } else if(usePresD2){
                ptElemIntegr_[elemType]->PerformIntegrationPresD2(coordMat,
                                                                  nodalPressureD2,
                                                                  elemVecLH,
                                                                  nodalLoadDensity);

              }
            }

          if(computeAPEMass){
            ptElemIntegr_[elemType]->PerformIntegrationAPEMass(coordMat,
                                                               nodalVel,
                                                               nodalPerturbedPressure,
                                                               nodalPressureTDeriv,
                                                               elemVecPres,
                                                               density);
          }

          if(computeAPEMomentum){
            ptElemIntegr_[elemType]->PerformIntegrationAPEMomentum(coordMat,
                                                                   nodalVel,
                                                                   nodalMeanVel,
                                                                   elemVecLamb,
                                                                   elemVecLambRhs,
                                                                   density);

          }

          if(computeAeroAcouSrc){
            ptElemIntegr_[elemType]->PerformIntegrationAeroAcouSrc(coordMat,
                                                                   nodalVel,
                                                                   nodalMeanVel,
                                                                   elemVecAeroAcou);
          }
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
        elemVecLH.Resize(numElemNodes);
        elemVecLH.Init();
      }

      // Add contributions of all element nodes
      for( UInt n=0; n<numElemNodes; n++)
      {
        nodeNum = topology_[elemIdx * maxNENodes + n];
        idx = regionNodeIndices_[regionIdx][nodeNum];

#ifndef NDEBUG
        if (std::isnan(elemVecLH[n]) || std::isinf(elemVecLH[n])) {
          EXCEPTION("Source term calculated on element " << i+1
                    << " is Inf or Nan.");
        }
#endif

        if(computeLHV || computeLHP || useDivLHT || usePresD2){
#pragma omp atomic 
          acouRhsField[idx] -= elemVecLH[n];
#pragma omp atomic 
          acouRhsDensityField[idx] -= nodalLoadDensity[n];
        }

        if(computeAeroAcouSrc){
#pragma omp atomic
          aeroAcouRhsField[idx] -= elemVecAeroAcou[n];
        }

        if(computeAPEMass){
          if(presFieldAvailable){
#pragma omp atomic
            acouMixedMassRhsField[idx] -= elemVecPres[n];
          }
        }

        if(computeAPEMomentum){
          for(UInt d =0;d<dim_;++d){
#pragma omp atomic
            acouLambVecRhsField[idx*dim_ + d] -= elemVecLambRhs[n*dim_+d];
          }
        }
      }
      
      // Add contributions of elements
      for( UInt n=0; n < dim_; n++)
      {
        if(computeLHV)
          acouDivLighthillTensor[i*dim_ + n] = divLHTensor[n];
        if(computeAPEMomentum)
          acouLambVec[i*dim_+n] = elemVecLamb[n];
      }      
    }
}//end of parallel region

    std::cout << "done." << std::endl;
  }

  //! This method calculates a surface mechanical force for FSI one-way coupling
  void CouplingHandler::CalculateMechSurfaceForce(const UInt surfRegionIdx,
      FlowDataType& flowData){

    //we go as follows:
    /* 1. Check if FLUIDMECH_FORCE is available in the flowData struct
     * 2. Check if we really hava a surface region
     * 3. loop over all surface elements and compute a mass integrator on the surface elements
     * 4. add the result to the result struct
     */
    std::string regionName = ptFileReader_->GetRegionName(surfRegionIdx);
    UInt nElems = ptFileReader_->GetNumElems(surfRegionIdx);
    FlowDataPartStruct& forceStruct = flowData[FLUIDMECH_FORCE];
    std::vector<Double>& forceField = forceStruct.data;

    bool computeForce = ( std::find(outputFields_.begin(),outputFields_.end(),"mechRhsLoad") != outputFields_.end() ||
                                std::find(outputFields_.begin(),outputFields_.end(),"all") != outputFields_.end());


    if(computeForce && !forceStruct.isActive){
      EXCEPTION("the force field is not available. Going to abort");
    }

    if(!forceStruct.isActive){
      flowData.erase(FLUIDMECH_FORCE);
      return;
    }

    // Init Mechanic rhs load structures
    FlowDataPartStruct& fdps2 = flowData[MECH_RHS_LOAD];
    fdps2.isActive = true; // all partitions have results
    fdps2.definedOn = ResultInfo::NODE; // nodes
    if (fdps2.dofNames.empty()) {
      fdps2.dofNames.push_back("x");
      fdps2.dofNames.push_back("y");
      if (dim_ == 3)
        fdps2.dofNames.push_back("z");
    }
    fdps2.unit = MapSolTypeToUnit(MECH_RHS_LOAD);
    fdps2.resultName = "mechRhsLoad";
    fdps2.data.resize(numRegionNodes_[surfRegionIdx] * dim_);
    fdps2.entryType = ResultInfo::VECTOR;
    std::vector<Double>& mechRhsField = fdps2.data;
    std::fill(mechRhsField.begin(), mechRhsField.end(), 0.0);


    UInt forceIdx,topoIdx,idx;
    Matrix<Double> coordMat;
    //nodal values for force
    Matrix<Double> nodalForce;
    Vector<Double> elemVecForce;
    Elem::FEType elemType;
    UInt numElemNodes;
    //UInt elemDim;
    UInt elemIdx;
    UInt maxNENodes = ptFileReader_->GetMaxNumElemNodes();
    UInt nodeNum;

    std::cout << "Calculating mechanical surface sources on " << regionName << " ";
    std::cout << "...serial loop over " << nElems << " elements..." ;
    std::cout.flush();

    for(UInt i=0; i<nElems; i++)
    {
      elemIdx = regionElems_[surfRegionIdx][i] - 1;
      elemType = (Elem::FEType) elemTypes_[elemIdx];
      numElemNodes = Elem::GetNumElemNodes(elemType);
      //elemDim = Elem::GetElemDim(elemType);

      coordMat.Resize(dim_, numElemNodes);
      nodalForce.Resize(dim_, numElemNodes);
      for( UInt n=0; n<numElemNodes; n++)
      {
        nodeNum = topology_[elemIdx * maxNENodes + n];
        topoIdx = (nodeNum - 1) * 3;
        forceIdx = regionNodeIndices_[surfRegionIdx][nodeNum] * dim_;
        for( UInt d=0; d<dim_; d++)
        {
          coordMat[d][n] = nodalCoords_[topoIdx+d];
          nodalForce[d][n] = forceField[forceIdx+d];

        }
      }
      try{
      ptElemIntegr_[elemType]->PerformIntegrationMechRhs(coordMat,
          nodalForce,
          elemVecForce);
      }catch(CoupledField::Exception &ex){
        std::cerr << "WARN: An Exception occurred during mechanical source term "
                          << "computation:\nElement " << elemIdx+1 << std::endl;

       std::cerr << ex.what()<< std::endl;
      }
      for( UInt n=0; n<numElemNodes; n++)
      {
        nodeNum = topology_[elemIdx * maxNENodes + n];
        idx = regionNodeIndices_[surfRegionIdx][nodeNum];
        for(UInt d =0;d<dim_;++d){
          mechRhsField[idx*dim_ + d] += elemVecForce[n*dim_+d];
        }
      }
    }
    std::cout << "done!" << std::endl;

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
    // UInt surfElemDim;
    UInt volElemDim;
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
      // surfElemDim =  Elem::GetElemDim(surfElemType);

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

  void CouplingHandler::ComputePerturbedPressureField(const std::vector<Double> & actPresField,const int regIdx,
                                                      FlowDataPartStruct& fdps){
    std::vector<Double>& meanPressureField = fdps.data;
    const UInt size = meanPressureField.size();
    std::vector<Double> & pertPres = perturbedPressureField_[regIdx];
    if(pertPres.size() == 0){
      pertPres.resize(size);
    }
    int iSize = (int)size;
#pragma omp parallel for
    for(int i=0;i<iSize;++i){
      pertPres[i] = actPresField[i] - meanPressureField[i];
    }
  }

  void CouplingHandler::UpdatePressureTimeDeriv(const std::vector<Double> & actPresField, const int regIdx){

    if(actPresField.size() == 0){
      return;
    }
    Settings& settings = Settings::Instance();
    if(settings.GetDouble("timestep") == 0){
      std::cerr << "Got zero timestep value from fielreader. Aborting derivative calculation... " << std::endl;
    }
    const Double dt = settings.GetDouble("timestep");
    std::vector<Double> & myVec = pressureTimeDeriv_[regIdx];
    std::vector<Double> & oldVec = oldPressureField_n_1_[regIdx];
    std::vector<Double> & olderVec = oldPressureField_n_2_[regIdx];
    const UInt size = actPresField.size();
    bool firstStep = false;
    bool secondStep = false;
    if(oldVec.size() == 0){
      oldVec.resize(size);
      myVec.resize(size);
      olderVec.resize(size);
      firstStep=true;
    }

    if(!firstStep && olderVec.size() == 0){
      olderVec.resize(size);
      secondStep = true;
      olderVec.assign(oldVec.begin(),oldVec.end());
    }

    if(firstStep || secondStep){
      //ok, this is only first order accurate for the first two timesteps

      int iSize = (int)size;
  #pragma omp parallel for
      for(int i = 0; i < iSize; ++i){
        myVec[i] = (actPresField[i] - oldVec[i]) / dt;
        oldVec[i] = actPresField[i];
      }
    }else{
      const Double c1 = 3.0;
      const Double c2 = 4.0;
      const Double c3 = 1.0;
      const Double iDt = 1.0 / (2.0*dt);

      int iSize = (int)size;
#pragma omp parallel for
      for(int i = 0; i < iSize; ++i){
        myVec[i] = (c1 * actPresField[i] - c2 * oldVec[i] + c3*olderVec[i]) * iDt;
        olderVec[i] = oldVec[i];
        oldVec[i] = actPresField[i];
      }
    }
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
                 << "), because there are too many nodes ("
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
        UInt eNum = (UInt) std::floor((Double)matches[j] / (Double)maxNENodes);
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
