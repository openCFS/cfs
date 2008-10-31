
#include <string>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <iomanip>
#include <sstream>

#include "General/exception.hh"
#include "cplreader/Settings.hh"
#include "FileReader_FASTEST.hh"


namespace CoupledField
{

    FileReader_FASTEST::FileReader_FASTEST(const std::string& name,
                           const UInt dim,
                           const UInt numFiles,
                           const UInt startIndex) :
        FileReader(name, dim, numFiles),
        startIndex_(startIndex),
        numResults_(0)
    {
      partFmtStr_ = "";
      timeFmtStr_ = "";
    }

    FileReader_FASTEST::~FileReader_FASTEST()
    {
    }

    void FileReader_FASTEST::Init()
    {
      Settings& settings = Settings::Instance();
        std::string filename;
        char buf[128];
        UInt dummy;
        UInt i;


        if(settings.GetDouble("timestep") < 0)
          EXCEPTION("No proper time step has been specified! Use --timestep X.");

        // Let's first determine the format of filenames
        std::stringstream sstr;
        std::string regionFormatStr;
        std::string timeStepFormatStr;

        if (settings.GetInt("verbose")) {
          std::cout << "Trying to determine file name format for "
                    << ".coord/.node files... ";
        }

        for(i=1; i<10; i++)
        {
          filename = baseName_;
          sstr.str("");
          sstr << "%0" << i << "i";
          regionFormatStr = sstr.str();
          sprintf(buf, regionFormatStr.c_str(), 1);
          filename+= buf;
          filename+= ".coord";

          inFile_.clear();
          inFile_.open(filename.c_str());
          if (!inFile_) {
            //std::cerr << "Can't open " << filename << std::endl;
          }
          else {
            partFmtStr_ = regionFormatStr;
            break;
          }
        }

        if( i == 10 ) {
          EXCEPTION("Could not determine file name format for .coord/.node files");
        }
        else if (settings.GetInt("verbose")) {
          std::cout << "OK" << std::endl;
        }

        inFile_ >> dummy;
        inFile_ >> numRegions_;
        inFile_ >> dummy;
        numSteps_ = dummy < numSteps_ ? dummy : numSteps_;
        if (settings.GetInt("numsteps"))
        {
          const UInt tmp = (UInt) settings.GetInt("numsteps");
          /* only take argument if tmo does not exceed the maximal number of timesteps possible */
          if (tmp < numSteps_)
          {
            numSteps_ = tmp;
          }
        }

        inFile_.close();
        inFile_.clear();

        if (settings.GetInt("verbose")) {
          std::cout << "Trying to determine file name format for .dat files... ";
        }

        sstr.str("");
        sstr << startIndex_;

        for(i=sstr.str().length(); i<=10; i++)
        {
          filename = baseName_;
          sprintf(buf, regionFormatStr.c_str(), 1);
          filename+= buf;
          filename+= "_";
          sstr.str("");
          sstr << "%0" << i << "i";
          timeStepFormatStr = sstr.str();
          sprintf(buf, timeStepFormatStr.c_str(), startIndex_);
          filename+= buf;
          filename+= ".dat";

          inFile_.clear();
          inFile_.open(filename.c_str());
          if (!inFile_) {
            //std::cerr << "Can't open " << filename << std::endl;
          }
          else {
            timeFmtStr_ = timeStepFormatStr;
            break;
          }
        }

        if( i == 10 ) {
          EXCEPTION("Could not determine file name format for .dat files");
        }
        else if (settings.GetInt("verbose")) {
          std::cout << "OK" << std::endl;
        }

        inFile_ >> numResults_;

        inFile_.close();
        inFile_.clear();

        // Initialize mapping between columns in FASTEST result files and
        // internal variables.
        Integer dataCol;
        std::vector<std::string> fastestDOFs;
        fastestDOFs.push_back("lhsrc");
        fastestDOFs.push_back("vx");
        fastestDOFs.push_back("vy");
        fastestDOFs.push_back("vz");
        fastestDOFs.push_back("pres");

        solutionTypes_.push_back(ACOU_RHS_LOAD);
        solutionTypes_.push_back(FLUIDMECH_VELOCITY);
        solutionTypes_.push_back(FLUIDMECH_VELOCITY);
        solutionTypes_.push_back(FLUIDMECH_VELOCITY);
        solutionTypes_.push_back(FLUIDMECH_PRESSURE);

        dofIndices_.push_back(0);
        dofIndices_.push_back(0);
        dofIndices_.push_back(1);
        dofIndices_.push_back(2);
        dofIndices_.push_back(0);

        for(UInt i=0; i<fastestDOFs.size(); i++)
        {
          std::string colStr = settings.GetString(fastestDOFs[i]);

          if(colStr != "")
          {
            sstr.str("");
            sstr.clear();
            colStr = colStr.substr(3);
            sstr << colStr;
            sstr >> dataCol;

            if(sstr.fail())
            {
              EXCEPTION("Error while trying to init FASTEST column mapping");
            }

            dataColumns_.push_back(dataCol-1);
          }
          else
          {
            dataColumns_.push_back(-1);
          }
        }
        if ((dataColumns_[3] > -1) && (dim_ == 2))
          dataColumns_[3] = -1;

        std::cout << "\n Name:\t\t\t\t" << name_ << std::endl
                  << " Dimension:\t\t\t" << dim_ << std::endl
                  << " Files:\t\t\t\t" << numSteps_ << std::endl
                  << " Partitions:\t\t\t" << numRegions_ << std::endl
                  << " Results:\t\t\t" << numResults_ << std::endl
                  << " timeStep:\t\t\t" << settings.GetDouble("timestep") << std::endl
                  << " Lighthill source term column:\t" << dataColumns_[0]+1 << std::endl
                  << " vx column:\t\t\t" << dataColumns_[1]+1 << std::endl
                  << " vy column:\t\t\t" << dataColumns_[2]+1 << std::endl
                  << " vz column:\t\t\t" << dataColumns_[3]+1 << std::endl
                  << " pressure column:\t\t" << dataColumns_[4]+1 << std::endl << std::endl;

        elsize_.resize(numRegions_);
        numNodesPerRegion_.resize(numRegions_);
        numElemsPerRegion_.resize(numRegions_);

        for(UInt i=0; i<numRegions_; i++)
        {
            filename = baseName_;
            sprintf(buf, "%02i", i+1);
            filename+= buf;
            filename+= ".coord";

            inFile_.clear();
            inFile_.open(filename.c_str());
            if (!inFile_) {
              EXCEPTION("Can't open " << filename);
            }

            inFile_ >> dummy;
            inFile_ >> dummy;
            inFile_ >> dummy;
            inFile_ >> numNodesPerRegion_[i];

            inFile_.close();
            inFile_.clear();


            filename = baseName_;
            sprintf(buf, "%02i", i+1);
            filename+= buf;
            filename+= ".topo";

            inFile_.open(filename.c_str());
            if (!inFile_) {
              EXCEPTION("Can't open " << filename);
            }

            inFile_ >> dummy;
            inFile_ >> dummy;
            inFile_ >> dummy;
            inFile_ >> numElemsPerRegion_[i];
            inFile_ >> elsize_[i];

            inFile_.close();
            inFile_.clear();

            maxNumElemNodes_ = elsize_[i] > maxNumElemNodes_ ? elsize_[i] : maxNumElemNodes_;
            /*if (settings.GetInt("verbose")) {
              std::cout << "Partition " << (i+1)
                        << " nodes: " << numNodesPerRegion_[i]
                        << " elems: " << numElemsPerRegion_[i]
                        << " elsize: " << elsize_[i]
                        << std::endl;
            }*/

        }
    }


    void FileReader_FASTEST::ReadNodalCoords(std::vector<Double> & NODECOORD)
    {
        std::string filename;
        char buf[128];
        UInt dummy;
        UInt regionIdx;
        UInt numNodes = 0;
        UInt nodeOffset = 0;
        UInt nodeIdx = 0;

        for(regionIdx=0; regionIdx < numRegions_; regionIdx++)
          numNodes += numNodesPerRegion_[regionIdx];

        NODECOORD.resize(3 * numNodes);

        for(regionIdx=0; regionIdx < numRegions_; regionIdx++)
        {
          filename = baseName_;
          sprintf(buf, partFmtStr_.c_str(), regionIdx+1);
          filename+= buf;
          filename+= ".coord";

          inFile_.clear();
          inFile_.open(filename.c_str());
          if (!inFile_) {
            EXCEPTION("Can't open " << filename);
          }

          /* Set pointer to beginning of file: */
          std::string::size_type pos=0;
          inFile_.seekg(pos, std::ios::beg);

          inFile_ >> dummy;
          inFile_ >> dummy;
          inFile_ >> dummy;
          inFile_ >> dummy;

          std::vector<double> max, min;
          max.resize(3);
          min.resize(3);

          nodeIdx = nodeOffset * 3;
          inFile_ >> NODECOORD[nodeIdx + 0]
                  >> NODECOORD[nodeIdx + 1]
                  >> NODECOORD[nodeIdx + 2];
          min[0] = max[0] = NODECOORD[nodeIdx + 0];
          min[1] = max[1] = NODECOORD[nodeIdx + 1];
          min[2] = max[2] = NODECOORD[nodeIdx + 2];

          for (UInt i=1; i < numNodesPerRegion_[regionIdx]; i++)
          {
            nodeIdx += 3;

            for (UInt j=0; j < 3; j++)
            {
              inFile_ >> NODECOORD[nodeIdx + j];
              if(NODECOORD[nodeIdx + j] > max[j])
                max[j] = NODECOORD[nodeIdx + j];
              if(NODECOORD[nodeIdx + j] < min[j])
                min[j] = NODECOORD[nodeIdx + j];
              //                if (i<=10) std::cout<<"NODECOORD["<<3*i+j<<"]: "<<NODECOORD[3*i+j]<<std::endl;
            }
          }

          /*if (settings.GetInt("verbose")) {
          std::cout << "MAX COORD: ( "<< max[0] << ", " << max[1] << ", "<< max[2] << ")"<< std::endl;
          std::cout << "MIN COORD: ( "<< min[0] << ", " << min[1] << ", "<< min[2] << ")"<< std::endl;
        }*/


          inFile_.close();
          nodeOffset += numNodesPerRegion_[regionIdx];
        }
    }

    void FileReader_FASTEST::ReadTopology(std::vector<UInt> & TOPOLOGYDATA,
                                               std::vector<UInt> & elemTypes)
    {
      std::string filename;
      char buf[128];
      UInt dummy, numNodes, elemType = ET_UNDEF;
      UInt regionIdx;
      UInt numElems = 0;
      UInt elemIdx = 0;
      UInt nodeOffset = 0;

      for(regionIdx=0; regionIdx < numRegions_; regionIdx++)
        numElems += numElemsPerRegion_[regionIdx];

      TOPOLOGYDATA.resize(maxNumElemNodes_ * numElems);
      elemTypes.resize(numElems);

      for(regionIdx = 0; regionIdx < numRegions_; regionIdx++) {
        filename = baseName_;
        sprintf(buf, partFmtStr_.c_str(), regionIdx+1);
        filename+= buf;
        filename+= ".topo";

        inFile_.clear();
        inFile_.open(filename.c_str());
        if (!inFile_) {
          EXCEPTION("Can't open " << filename);
        }

        /* Set pointer to beginning of file: */
        std::string::size_type pos=0;
        inFile_.seekg(pos, std::ios::beg);

        inFile_ >> dummy >> dummy >> dummy >> dummy >> dummy;

        numNodes = elsize_[regionIdx];
        switch(numNodes)
        {
        case 2:
          elemType = ET_LINE2;
          break;
        case 3:
          elemType = ET_TRIA3;
          break;
        case 4:
          if(dim_ == 3)
            elemType = ET_TET4;
          else
            elemType = ET_QUAD4;
          break;
        case 8:
          if (dim_ == 3)
            elemType = ET_HEXA8;
          else
            elemType = ET_QUAD8;
          break;
        case 5:
          elemType = ET_PYRA5;
          break;
        case 6:
          if (dim_ == 3)
            elemType = ET_WEDGE6;
          else
            elemType = ET_TRIA6;
          break;
        }

        for (UInt i=0; i < numElemsPerRegion_[regionIdx]; i++)
        {
          for (UInt j=0; j < elsize_[regionIdx]; j++)
          {
            inFile_ >> TOPOLOGYDATA[(elemIdx * maxNumElemNodes_) + j];
            TOPOLOGYDATA[(elemIdx * maxNumElemNodes_) + j] += nodeOffset;
          }

          elemTypes[elemIdx++] = elemType;
        }

        inFile_.close();
        nodeOffset += numNodesPerRegion_[regionIdx];
      }
    }

  void FileReader_FASTEST::GetRegionElements(std::vector<UInt> & regionElements,
                                   const UInt regionIdx)
  {
    UInt elemOffset = 0;

    for(UInt i=0; i < regionIdx; i++)
      elemOffset += numElemsPerRegion_[i];

    regionElements.resize(numElemsPerRegion_[regionIdx]);

    for(UInt i=0; i < numElemsPerRegion_[regionIdx]; i++)
      regionElements[i] = elemOffset + i + 1;
  }

  //! get nodal values from the corresponding fluid datafile the new way
  void FileReader_FASTEST::ReadNodalValues(std::vector<FlowDataType>& nodalFlowData,
                                           const std::vector<bool>& activeParts,
                                           const UInt timeStepIdx)
  {
    std::string filename;
    char buf[128];
    UInt dummy;
    UInt numDOFs;
    static std::vector<Double> tempVec;

    tempVec.resize(numResults_);

    for(UInt actRegion=0; actRegion < numRegions_; actRegion++)
    {
      //      std::cout << "FASTEST: actRegion " << actRegion << std::endl;

      // Initialize flow data structures if necessary
      FlowDataType& fd = nodalFlowData[actRegion];

      for(UInt j=0; j<dataColumns_.size(); j++)
      {
        if(dataColumns_[j] > -1)
        {
          if(fd.find(solutionTypes_[j]) == fd.end())
          {
            FlowDataPartStruct& fdps = fd[solutionTypes_[j]];
            fdps.isActive = true; // all partitions have results
            fdps.definedOn = ResultInfo::NODE; // nodes
            Enum2String(solutionTypes_[j], fdps.resultName);

            switch(solutionTypes_[j])
            {
            case ACOU_RHS_LOAD:
              fdps.dofNames.push_back("-");
              fdps.unit = MapSolTypeToUnit(ACOU_RHS_LOAD);
              fdps.entryType = ResultInfo::SCALAR;
              break;
            case FLUIDMECH_VELOCITY:
              fdps.dofNames.push_back("x");
              fdps.dofNames.push_back("y");
              if (dataColumns_[3] != -1)
                fdps.dofNames.push_back("z");
              fdps.unit = MapSolTypeToUnit(FLUIDMECH_VELOCITY);
              fdps.entryType = ResultInfo::VECTOR;
              break;
            case FLUIDMECH_PRESSURE:
              fdps.dofNames.push_back("-");
              fdps.unit = MapSolTypeToUnit(FLUIDMECH_PRESSURE);
              fdps.entryType = ResultInfo::SCALAR;
              break;
            default:
              break;
            }

            numDOFs = fdps.dofNames.size();
            fdps.data.resize(numDOFs * numNodesPerRegion_[actRegion]);

            //            std::cout << "FASTEST: fdps.resultName " << fdps.resultName << std::endl;
            //            std::cout << "FASTEST: fdps.unit " << fdps.unit << std::endl;
            //            std::cout << "FASTEST: numDOFs " << numDOFs << std::endl;

          }
        }
      }

      // Open data file
      filename = baseName_;
      sprintf(buf, partFmtStr_.c_str(), actRegion+1);
      filename+= buf;
      filename+= "_";
      sprintf(buf, timeFmtStr_.c_str(), timeStepIdx+startIndex_);
      filename+= buf;
      filename+= ".dat";

      inFile_.clear();
      inFile_.open(filename.c_str());
      if (!inFile_) {
        EXCEPTION("Can't open " << filename);
      }

      /* Set pointer to beginning of file: */
      std::string::size_type pos=0;
      inFile_.seekg(pos, std::ios::beg);

      inFile_ >> dummy;

      for (UInt i=0; i < numNodesPerRegion_[actRegion]; i++)
      {
        // read all data columns from input file
        for(UInt j=0; j<numResults_; j++)
          inFile_ >> tempVec[j];

        for(UInt j=0; j<dataColumns_.size(); j++)
        {
          if(dataColumns_[j] > -1)
          {
            numDOFs = nodalFlowData[actRegion][solutionTypes_[j]].dofNames.size();
            std::vector<Double>& data = nodalFlowData[actRegion][solutionTypes_[j]].data;

            data[i*numDOFs+dofIndices_[j]] = tempVec[dataColumns_[j]];
          }
        }
      }

      inFile_.close();
    }
  }

} // end of namespace
