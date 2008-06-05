#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>

#include "../params.hh"
#include "../settings.hh"
#include "FileReader_CFXexport.hh"

namespace CoupledField {

  FileReader_CFXexport::FileReader_CFXexport(const std::string& name,
                                             const UInt dim,
                                             const UInt numFiles,
                                             const UInt startIndex) :
    FileReader(name, dim, numFiles), startIndex_(startIndex)
  {
    colX_ = 0;
    colY_ = 0;
    colZ_ = 0;
    colDens_ = 0;
    colPres_ = 0;
    colVelX_ = 0;
    colVelY_ = 0;
    colVelZ_ = 0;
    
    basename_ = "./" + name_ + "_data/results_";
  }
  
  FileReader_CFXexport::~FileReader_CFXexport()
  {
  }

  void FileReader_CFXexport::Init() {
    
    Settings& settings = Settings::Instance();
    
    if(settings.GetDouble("timeStep") < 0)
      EXCEPTION("No proper time step has been specified! Use --timestep X.");

    // put filename together
    std::ostringstream sFilename;
    sFilename.str("");
    sFilename << basename_ << startIndex_ << ".csv"; 
    
    // try to open file
    infile.clear();
    infile.open(sFilename.str().c_str(), std::ios::in);
    if (!infile) {
      EXCEPTION("Cannot open file '" << sFilename.str() << "'.");
    }
    infile.seekg(0, std::ios::beg);
    
    numPartitions_ = 1;
    elsize_.clear(); elsize_.push_back(4);
    MpCCInodes_.resize(1);
    MpCCIelems_.resize(1);
    
    // find [Data] section
    std::string lineBuf;
    while (infile.good()) {
      std::getline(infile, lineBuf);
      if (lineBuf == "[Data]")
          break;
    }

    if ( infile.bad() || infile.eof() ) {
      EXCEPTION("Cannot find [Data] section in file '" << sFilename.str()
                << "'.");
    }
    
    // get line with column headers
    while (infile.good()) {
      std::getline(infile, lineBuf);
      if (lineBuf.length() >0 && lineBuf[0] != '#')
        break;
    }
    infile.close();

    UInt curCol = 1;
    std::string sColName;
    std::istringstream strLine;
    strLine.str(lineBuf);
    numResults_ = 0;

    while (std::getline(strLine, sColName, ',')) {
      if (sColName == " X [ m ]") {
        if (colX_ == 0)
          colX_ = curCol;
      }
      else if (sColName == " Y [ m ]") {
        if (colY_ == 0)
          colY_ = curCol;
      }
      else if (sColName == " Z [ m ]") {
        if (colZ_ == 0)
          colZ_ = curCol;
      }
      else if (sColName == " Density [ kg m^-3 ]") {
        if (colDens_ == 0) {
          colDens_ = curCol;
          col2Quan_[curCol] = DENSITY;
          ++numResults_;
        }
      }
      else if (sColName == " Pressure [ Pa ]") {
        if (colPres_ == 0) {
          colPres_ = curCol;
          col2Quan_[curCol] = PRESSURE;
          ++numResults_;
        }
      }
      else if (sColName == " Velocity u [ m s^-1 ]") {
        if (colVelX_ == 0) {
          colVelX_ = curCol;
          col2Quan_[curCol] = VELOCITY_X;
          ++numResults_;
        }
      }
      else if (sColName == " Velocity v [ m s^-1 ]") {
        if (colVelY_ == 0) {
          colVelY_ = curCol;
          col2Quan_[curCol] = VELOCITY_Y;
        }
      }
      else if (sColName == " Velocity w [ m s^-1 ]") {
        if (colVelZ_ == 0) {
          colVelZ_ = curCol;
          col2Quan_[curCol] = VELOCITY_Z;
        }
      }

      ++curCol;
    }

    if (colX_ == 0 || colY_ == 0 || (dim_ == 3 && colZ_ == 0)) {
      EXCEPTION("Cannot find columns of node coordinates in file '"
                << sFilename.str() << "'.");
    }

    if (colVelX_ + colVelY_ + colVelZ_ > 0) {
      if (colVelX_ == 0 || colVelY_ == 0 || (dim_ == 3 && colVelZ_ == 0)) {
        EXCEPTION("Cannot find all columns of quantity 'fluidMechVelocity' in "
                  << "file '" << sFilename.str() << "'.");
      }
    }

    if (numResults_ == 0) {
      EXCEPTION("No known quantities found in file '"
                << sFilename.str() << "'.");
    }

    std::cout << " Name:      \t" << name_ << std::endl
              << " Dimension: \t" << dim_ << std::endl
              << " Files:     \t" << numFiles_ << std::endl
              << " Partitions:\t" << numPartitions_ << std::endl
              << " Time step: \t" << settings.GetDouble("timeStep") << std::endl
              << " Results:   \t" << numResults_ << std::endl;
    if (colDens_ > 0)
      std::cout << "   Density: \tcolumn " << colDens_ << std::endl;
    if (colPres_ > 0)
      std::cout << "   Pressure:\tcolumn " << colPres_ << std::endl;
    if (colVelX_ > 0) {
      std::cout << "   Velocity:\tcolumns " << colVelX_ << "," << colVelY_;
      if (colVelZ_ > 0)
        std::cout << "," << colVelZ_;
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }
  
  void FileReader_CFXexport::ReadNodalCoords(std::vector<Double> & NODECOORD,
                                             const UInt partitionIdx)
  {
    if (partitionIdx >= numPartitions_) {
      EXCEPTION("Invalid partition index: " << partitionIdx);
    }

    // put filename together
    std::ostringstream sFilename;
    sFilename.str("");
    sFilename << basename_ << startIndex_ << ".csv"; 
    
    // try to open file
    infile.clear();
    infile.open(sFilename.str().c_str(), std::ios::in);
    if (!infile) {
      EXCEPTION("Cannot open file '" << sFilename.str() << "'.");
    }
    infile.seekg(0, std::ios::beg);
    
    // find [Data] section
    std::string lineBuf;
    while (infile.good()) {
      std::getline(infile, lineBuf);
      if (lineBuf == "[Data]")
          break;
    }

    if ( infile.bad() || infile.eof() ) {
      EXCEPTION("Cannot find [Data] section in file '" << sFilename.str()
                << "'.");
    }
    
    // get line with column headers
    while (infile.good()) {
      std::getline(infile, lineBuf);
      if (lineBuf.length() > 0 && lineBuf[0] != '#')
        break;
    }

    UInt lineNum = 0, fieldNum, maxField = (dim_ == 3 ? colZ_ : colY_);
    std::istringstream strLine;
    std::string fieldBuf;
    NODECOORD.clear();

    // cycle through all data lines
    while (infile.good()) {
      std::getline(infile, lineBuf);

      if (lineBuf.length() > 0 && !infile.eof()) {
        if (lineBuf[0] == '[') {
          break;
        }
        else if (lineBuf[0] != '#') {
          strLine.clear();
          strLine.str(lineBuf);
          fieldNum = 1;

          NODECOORD.push_back(0.0);
          NODECOORD.push_back(0.0);
          NODECOORD.push_back(0.0);

          while (strLine.good()) {
            std::getline(strLine, fieldBuf, ',');

            if (fieldNum == colX_)
              NODECOORD[lineNum*3  ] = std::atof(fieldBuf.c_str());
            else if (fieldNum == colY_)
              NODECOORD[lineNum*3+1] = std::atof(fieldBuf.c_str());
            else if (fieldNum == colZ_)
              NODECOORD[lineNum*3+2] = std::atof(fieldBuf.c_str());

            if (++fieldNum > maxField)
              break;
          }

          ++lineNum;
        }
      }
    }

    infile.close();

    MpCCInodes_[partitionIdx] = lineNum;
  }

  void FileReader_CFXexport::ReadTopology(std::vector<UInt> & TOPOLOGYDATA,
                                          std::vector<UInt> & numNodesPerElem,
                                          std::vector<UInt> & elemTypes,
                                          const UInt partitionIdx)
  {
    if (partitionIdx >= numPartitions_) {
      EXCEPTION("Invalid partition index: " << partitionIdx);
    }

    
    // put filename together
    std::ostringstream sFilename;
    sFilename.str("");
    sFilename << basename_ << startIndex_ << ".csv"; 
    
    // try to open file
    infile.clear();
    infile.open(sFilename.str().c_str(), std::ios::in);
    if (!infile) {
      EXCEPTION("Cannot open file '" << sFilename.str() << "'.");
    }
    infile.seekg(0, std::ios::beg);
    
    // find [Faces] section
    std::string lineBuf;
    while (infile.good()) {
      std::getline(infile, lineBuf);
      if (lineBuf == "[Faces]")
          break;
    }

    if ( infile.bad() || infile.eof() ) {
      EXCEPTION("Cannot find [Faces] section in file '" << sFilename.str()
                << "'.");
    }
    
    UInt lineNum = 0;
    std::istringstream strLine;
    std::string fieldBuf;
    std::vector<UInt> connect;
    
    TOPOLOGYDATA.clear();
    numNodesPerElem.clear();
    elemTypes.clear();

    // cycle through all data lines
    while (infile.good()) {
      std::getline(infile, lineBuf);

      if (lineBuf.length() > 0 && !infile.eof()) {
        if (lineBuf[0] == '[') {
          break;
        }
        else if (lineBuf[0] != '#') {
          strLine.clear();
          strLine.str(lineBuf);
          connect.clear();

          while (strLine.good()) {
            std::getline(strLine, fieldBuf, ',');
            connect.push_back(std::atoi(fieldBuf.c_str())+1);
          }

          switch (connect.size()) {
            case 4:
              numNodesPerElem.push_back(4);
              elemTypes.push_back(ET_QUAD4);
              ++lineNum;
              break;
            case 3:
              numNodesPerElem.push_back(3);
              elemTypes.push_back(ET_TRIA3);
              ++lineNum;
              break;
            case 0:
              break;
            default:
              std::cerr << "WARNING: Element with " << connect.size()
                        << " nodes was not converted.\n";
              connect.clear();
              break;
          }
          
          for (UInt i = 0, iEnd = connect.size(); i < iEnd; ++i) {
            TOPOLOGYDATA.push_back(connect[i]);
          }
        }
      }
    }

    infile.close();

    MpCCIelems_[partitionIdx] = lineNum;
  }

  void FileReader_CFXexport::ReadNodalValues(std::vector<FlowDataType>
                                               &nodalFlowdata,
                                             const UInt timeStepIdx)
  {
    // put filename together
    std::ostringstream sFilename;
    sFilename.str("");
    sFilename << basename_ << timeStepIdx << ".csv"; 
    
    // try to open file
    infile.clear();
    infile.open(sFilename.str().c_str(), std::ios::in);
    if (!infile) {
      EXCEPTION("Cannot open file '" << sFilename.str() << "'.");
    }
    infile.seekg(0, std::ios::beg);
    
    // find [Data] section
    std::string lineBuf;
    while (infile.good()) {
      std::getline(infile, lineBuf);
      if (lineBuf == "[Data]")
          break;
    }

    if ( infile.bad() || infile.eof() ) {
      EXCEPTION("Cannot find [Data] section in file '"
                << sFilename.str() << "'.");
    }
    
    // get line with column headers
    while (infile.good()) {
      std::getline(infile, lineBuf);
      if (lineBuf.length() > 0 && lineBuf[0] != '#')
        break;
    }

    UInt lineNum = 0, fieldNum, numVelDofs = (colVelZ_ > 0 ? 3 : 2);
    std::istringstream strLine;
    std::string fieldBuf;
    std::map<UInt, Quantity>::const_iterator itQuan;
    std::map<UInt, Quantity>::const_iterator itEnd = col2Quan_.end();
    FlowDataType &flowdata = nodalFlowdata[0];
    
    // set up flow data struct
    for (itQuan = col2Quan_.begin(); itQuan != itEnd; ++itQuan) {
      FlowDataPartStruct *fdps = NULL;

      switch (itQuan->second) {
      case DENSITY:
        fdps = &flowdata[FLUIDMECH_DENSITY];
        fdps->isActive = true;
        fdps->definedOn = ResultInfo::NODE;
        fdps->dofNames.push_back("-");
        fdps->unit = MapSolTypeToUnit(FLUIDMECH_DENSITY);
        fdps->data.resize(MpCCInodes_[0]);
        Enum2String(FLUIDMECH_DENSITY, fdps->resultName);
        break;
      case PRESSURE:
        fdps = &flowdata[FLUIDMECH_PRESSURE];
        fdps->isActive = true;
        fdps->definedOn = ResultInfo::NODE;
        fdps->dofNames.push_back("-");
        fdps->unit = MapSolTypeToUnit(FLUIDMECH_PRESSURE);
        Enum2String(FLUIDMECH_PRESSURE, fdps->resultName);
        fdps->data.resize(MpCCInodes_[0]);
        break;
      case VELOCITY_X:
        fdps = &flowdata[FLUIDMECH_VELOCITY];
        fdps->isActive = true;
        fdps->definedOn = ResultInfo::NODE;
        fdps->dofNames.push_back("x");
        fdps->dofNames.push_back("y");
        if (numVelDofs == 3)
          fdps->dofNames.push_back("z");
        fdps->unit = MapSolTypeToUnit(FLUIDMECH_VELOCITY);
        Enum2String(FLUIDMECH_VELOCITY, fdps->resultName);
        fdps->data.resize(numVelDofs * MpCCInodes_[0]);
        break;
      default:
        break;
      }
    }
    
    // cycle through all data lines
    while (infile.good()) {
      std::getline(infile, lineBuf);

      if (lineBuf.length() > 0 && !infile.eof()) {
        if (lineBuf[0] == '[') {
          break;
        }
        else if (lineBuf[0] != '#') {
          strLine.clear();
          strLine.str(lineBuf);
          fieldNum = 1;
          itQuan = col2Quan_.begin();

          while (strLine.good()) {
            std::getline(strLine, fieldBuf, ',');

            if (fieldNum == itQuan->first) {
              FlowDataPartStruct *fdps;
              switch (itQuan->second) {
              case DENSITY:
                fdps = &flowdata[FLUIDMECH_DENSITY];
                fdps->data[lineNum] = atof(fieldBuf.c_str());
                break;
              case PRESSURE:
                fdps = &flowdata[FLUIDMECH_PRESSURE];
                fdps->data[lineNum] = atof(fieldBuf.c_str());
                break;
              case VELOCITY_X:
                fdps = &flowdata[FLUIDMECH_VELOCITY];
                fdps->data[lineNum*numVelDofs] = atof(fieldBuf.c_str());
                break;
              case VELOCITY_Y:
                fdps = &flowdata[FLUIDMECH_VELOCITY];
                fdps->data[lineNum*numVelDofs+1] = atof(fieldBuf.c_str());
                break;
              case VELOCITY_Z:
                fdps = &flowdata[FLUIDMECH_VELOCITY];
               fdps->data[lineNum*numVelDofs+2] = atof(fieldBuf.c_str());
                break;
              }
              ++itQuan;
            }

            if (itQuan == itEnd)
              break;
            
            ++fieldNum;
          }

          ++lineNum;
        }
      }
    }

    infile.close();
  }

  
} // namespace CoupledField
