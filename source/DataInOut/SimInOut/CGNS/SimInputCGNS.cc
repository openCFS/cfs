// =====================================================================================
// 
//       Filename:  SimInputCGNS.cc
// 
//    Description:  This is to read in a CGNS File
// 
//        Authors:  Simon Triebenbacher, simon.triebenbacher@tuwien.ac.at
//                  Andreas Hueppe (AHU), andreas.hueppe@uni-klu.ac.at
//                  Jens Grabinger, jens.grabinger@simetris.de
//        Company:  TU Wien, Universitaet Klagenfurt, SIMetris GmbH
//
// =====================================================================================

#include <sstream>
#include <iostream>
#include <vector>
#include <string>
#include <unordered_set>

#include <filesystem>
#include <boost/regex.hpp>
#include <boost/tokenizer.hpp>
namespace fs=std::filesystem;
#include <boost/algorithm/string/predicate.hpp>
namespace algo=boost::algorithm;

#include "SimInputCGNS.hh"

#include "DataInOut/Logging/LogConfigurator.hh"
#include "Utils/tools.hh"


// Size of string buffers (32 characters + terminating null byte)
// see https://cgns.github.io/CGNS_docs_current/midlevel/general.html#character
#define STRBUFSIZE 33

// Macro for checking return value of CGNS library functions that throws an
// exception on error
#define CGNS_CHECK_EX(call) \
  if (call != CG_OK) { \
    EXCEPTION("CGNS file '" << fileName_ << "':\n" << cg_get_error()); \
  }

// Macro for checking return value of CGNS library functions that displays a
// warning on error
#define CGNS_CHECK_WARN(call) \
  if (call != CG_OK) { \
    WARN("CGNS file '" << fileName_ << "':\n" << cg_get_error()); \
  }

// Macro for checking return value of CGNS library functions that displays a
// warning on error only in debug mode
#ifdef NDEBUG
# define CGNS_CHECK_DBG(call) call;
#else
# define CGNS_CHECK_DBG(call) \
    if (call != CG_OK) { \
      WARN("CGNS file '" << fileName_ << ":\n" << cg_get_error()); \
    }
#endif

namespace CoupledField {

  // static variables
  std::map<Integer, Elem::FEType> SimInputCGNS::elemTypeMap_;
  
  // declare logging stream
  DEFINE_LOG(simInputCGNS, "simInput.CGNS")

  //! Constructor
  SimInputCGNS::SimInputCGNS( std::string fileName,
                              PtrParamNode inputNode,
                              PtrParamNode infoNode ) :
    SimInput(fileName, inputNode, infoNode ),
    fileHandle_(0),
    cellDim_(0),
    numBases_(0),
    baseIndex_(1),
    numZones_(0),
    numNodes_(0),
    numElems_(0),
    numRegions_(0),
    nodeOffset_(0),
    elemOffset_(0),
    analysis_(BasePDE::NO_ANALYSIS)
  {
    capabilities_.insert( SimInput::MESH );
    capabilities_.insert( SimInput::MESH_RESULTS );

    InitElemTypeMap();
    
    numElemsOfDim_.resize(3);
    std::fill(numElemsOfDim_.begin(), numElemsOfDim_.end(), 0);
  }
  
  SimInputCGNS::~SimInputCGNS(){
    if (fileHandle_) {
      CGNS_CHECK_WARN( cg_close(fileHandle_) );
    }
  }
  
  // Error callback routine called by CGNS library
  void SimInputCGNS::ErrorCallback(int isError, char *errMsg) {
    if (isError == -1) {
      std::cerr << errMsg << std::endl;
    }
#ifndef NDEBUG
    else {
      WARN(errMsg);
    }
#endif
  }
  
  void SimInputCGNS::InitModule(){
    LOG_TRACE(simInputCGNS) << "Entering SimInputCGNS::InitModule" << std::endl;

    // Set our custom error callback routine
    CGNS_CHECK_DBG( cg_error_handler(&ErrorCallback) );
    
    // Try to find files that match the given filename pattern
    std::string baseName;
    try {
      // Get the directory part of the path
      fs::path absPath = fileName_.empty() ? fs::current_path() : fs::absolute(fileName_);
      fs::path dirPath = absPath.parent_path();

      // Get the filename path of the path and turn it into a regex
      fs::path filePath = absPath.filename();
      if (!filePath.has_extension()) {
        filePath.replace_extension("cgns");
      }
      baseName = filePath.string();
      std::string regexStr = PathPatternToRegEx(baseName);
      boost::regex fileRegex(regexStr, boost::regex::extended);
      
      std::string entry;
      fs::directory_iterator dirIt(dirPath), dirEnd;
      for ( ; dirIt != dirEnd; ++dirIt) {
        // Skip anything but regular files
        if (!fs::is_regular_file(dirIt->path())) continue;
        
        // test against the regex
        entry = dirIt->path().filename().string();
        if (boost::regex_match(entry, fileRegex)) {
          try {
            ReadFileContents(dirIt->path().string());
          }
          catch (Exception &ex) {
            WARN(ex.GetMsg());
          }
        }
      }
    }
    catch (fs::filesystem_error& ex) {
      EXCEPTION(ex.what());
      return;
    }
    
    if (fileNames_.size() == 1) {
      analysis_ = BasePDE::STATIC;
    }
    else if (fileNames_.size() == 0) {
      EXCEPTION("Could not find any CGNS files matching the pattern '"
          << baseName << "' which contain usable data.");
    }
    
    // Now create the global step map
    std::set<Double> steps;
    
    std::map< shared_ptr<ResultInfo>, std::set<Double> >::iterator
        resIt = resultSteps_.begin(),
        resEnd = resultSteps_.end();
    for ( ; resIt != resEnd; ++resIt) {
      steps.insert(resIt->second.begin(), resIt->second.end());
    }
    
    UInt step = 1;
    std::set<Double>::iterator stepIt = steps.begin(), stepEnd = steps.end();
    for ( ; stepIt != stepEnd; ++stepIt, ++step) {
      stepMap_[step] = *stepIt;
    }
  }

  void SimInputCGNS::OpenFile(const std::string &fileName) {
    if (fileHandle_) {
      if (fileName == fileName_) return; // nothing to do
      
      CGNS_CHECK_WARN( cg_close(fileHandle_) );
      fileHandle_ = 0;
      fileName_.clear();
    }
    
    if (!fs::exists(fileName)) {
      EXCEPTION("File not found: " << fileName);
    }
    
    Integer fileType;
    if (cg_is_cgns(fileName.c_str(), &fileType) != CG_OK) {
      EXCEPTION("Not a valid CGNS file: " << fileName);
    }
    
    if (cg_open(fileName.c_str(), CG_MODE_READ, &fileHandle_) != CG_OK) {
      EXCEPTION("Failed to open file: " << fileName << ",\n" << cg_get_error());
    }
    
    fileName_ = fileName;
    
    LOG_TRACE(simInputCGNS) << "Successfully opened file " << fileName;

    // Check CGNS version of opened file
    if (cg_version(fileHandle_, &fileVersion_) != CG_OK) {
      EXCEPTION("Failed to read CGNS version of file: " << fileName << ",\n" << cg_get_error());
    }
    LOG_TRACE(simInputCGNS) << "CGNS version: " << fileVersion_;
    if (fileVersion_ <= CGNS_VERSION/1000.0)
    {
      WARN("Old CGNS version: " << fileVersion_ << " (current: " << CGNS_VERSION/1000.0 << "), of file: " << fileName << std::endl
      << "If you encounter issues, use 'cgnsupdate' to update input file to supported version!");
    }
  }
  
  void SimInputCGNS::ReadFileContents(const std::string &fileName) {
    static bool warnedZoneType = false;
    static bool warnedMultiSol = false;
    static bool warnedGridLoc = false;
    
    OpenFile(fileName);

    // Check base datasets
    CGNS_CHECK_EX( cg_nbases(fileHandle_, &numBases_) );

    if (numBases_ == 0) {
      EXCEPTION("No datasets found in CGNS file '" << fileName << "'.");
    }
    if (numBases_ != 1) {
      WARN("Found " << numBases_ << " bases in CGNS file '" << fileName
           << "'.\nNACS can utilize only the first base." );
    }
    
    // Get name and dimensions of first base
    char firstBaseName[STRBUFSIZE];
    Integer physDim;
    CGNS_CHECK_EX( cg_base_read(fileHandle_, baseIndex_, firstBaseName,
                                &cellDim_, &physDim) );
    dim_ = (UInt) physDim;
    
    // Get time steps in the file
    Integer numSteps = 0;
    char baseIterName[STRBUFSIZE];
    CGNS_CHECK_DBG( cg_biter_read(fileHandle_, baseIndex_, baseIterName, &numSteps) );

    DataType_t dataType;
    StdVector<Double> steps(numSteps);
    
    if (numSteps > 0) {
      // Navigate to base iterative data node
      CGNS_CHECK_WARN( cg_goto(fileHandle_, baseIndex_, baseIterName, 0, "end") );
      
      Integer numArrays = 0;
      CGNS_CHECK_WARN( cg_narrays(&numArrays) );
      
      char arrayName[STRBUFSIZE];
      Integer a, dataDim;
      cgsize_t dimVec[12];
      for (a = 1; a <= numArrays; ++a) {
        CGNS_CHECK_WARN( cg_array_info(a, arrayName, &dataType, &dataDim, dimVec) );
        if (strcmp(arrayName, "TimeValues") == 0) {
          if (dataDim == 1 && dimVec[0] == numSteps) {
            CGNS_CHECK_EX( cg_array_read_as(a, RealDouble, (void*) steps.GetPointer()) );
            break;
          }
        }
      }
     
      if (a == numArrays) { // TimeValues array was not found
        numSteps = 0;
      }
    }
    
    // If there were transient data in other files before, we ignore files
    // without time step data
    if (analysis_ == BasePDE::TRANSIENT && numSteps == 0) {
      return;
    }
    
    // Just clear the variables that would accumulate with every time step
    numNodes_ = 0;
    numElems_ = 0;
    regionNames_.Clear(true);
    namedElems_.Clear(true);
    
    // Get number of zones
    CGNS_CHECK_EX( cg_nzones(fileHandle_, baseIndex_, &numZones_) );

    if (numZones_ == 0) {
      EXCEPTION("No grid found in CGNS file '" << fileName << "'.");
    }
    
    // Declare variables outside of loops
    bool gridFound = false, resultsFound = false;
    ZoneType_t zoneType;
    cgsize_t zoneSize[9];
    Integer numSols = 0, numFields = 0;
    UInt res;
    char zoneName[STRBUFSIZE], solName[STRBUFSIZE], fieldName[STRBUFSIZE];
    GridLocation_t solLocation;
    ResultInfo::EntityUnknownType definedOn;
    SolutionType solType;
    ResultInfo::EntryType entryType;
    
    for (Integer zone = 1; zone <= numZones_; ++zone) {
      
      // Check type of zone grid
      CGNS_CHECK_EX( cg_zone_type(fileHandle_, baseIndex_, zone, &zoneType) );

      if (zoneType != Unstructured) {
        if (!warnedZoneType) {
          WARN("Only unstructured grids in CGNS files are supported. Grid read from '"
              << fileName << "' may be incomplete.");
          warnedZoneType = true;
        }
        continue;
      }
      gridFound = true;
      
      // Get name and size of zone
      CGNS_CHECK_EX( cg_zone_read(fileHandle_, baseIndex_, zone, zoneName, zoneSize) );
      
      // Update number of nodes
      numNodes_ += zoneSize[0];
      
      // Get info about elements, regions and groups
      Integer nSections = 0;
      CGNS_CHECK_EX( cg_nsections(fileHandle_, baseIndex_, zone, &nSections) );
      
      char sectionName[STRBUFSIZE];
      ElementType_t eType;
      cgsize_t start, end;
      Integer nBoundary, parentFlag, numSecElems, numPolyFaces = 0;
      bool hasPolyhedra = false;
      StdVector<std::string> groupNames;
      std::string polyFacesName;
      
      for (Integer sec = 1; sec <= nSections; ++sec) {
        CGNS_CHECK_EX(cg_section_read(fileHandle_, baseIndex_, zone, sec,
                      sectionName, &eType, &start, &end, &nBoundary, &parentFlag));
        numSecElems = end - start + 1;

        switch (eType) {
          case NFACE_n:
            hasPolyhedra = true;
            break;
          case NGON_n:
            if (numSecElems > numPolyFaces) {
              numPolyFaces = numSecElems;
              polyFacesName = sectionName;
            }
            break;
          case ElementTypeNull:
          case ElementTypeUserDefined:
            continue;
          default:
            break;
        }

        numElems_ += numSecElems;
        
        if (strcmp(sectionName, zoneName) == 0) {
          regionNames_.push_back(zoneName);
        }
        else {
          groupNames.push_back(sectionName);
        }
      }

      // Copy group names to class variable, but leave out the polygon faces
      // if there are polyhedra
      nSections = groupNames.size();
      namedElems_.Reserve(nSections);
      for (Integer sec = 0; sec < nSections; ++sec) {
        if (hasPolyhedra && groupNames[sec] == polyFacesName) {
          numElems_ -= numPolyFaces;
          continue;
        }
        namedElems_.push_back(groupNames[sec]);
      }
      
      // Get the number of solutions
      CGNS_CHECK_EX( cg_nsols(fileHandle_, baseIndex_, zone, &numSols) );

      LOG_TRACE(simInputCGNS) << "Found " << numSols << " solutions in zone "
          << zone << std::endl;

      if (numSols < 1) continue;
      if (numSols > 1 && !warnedMultiSol) {
        WARN("NACS currently does not support reading CGNS data with more than "
            << "one time step per file. Data may be incomplete.");
        warnedMultiSol = true;
      }

      // Get name of solution (irrelevant) and location of data points
      CGNS_CHECK_EX( cg_sol_info(fileHandle_, baseIndex_, zone, 1, solName, &solLocation) );
      
      LOG_TRACE(simInputCGNS) << "Found solution '" << solName << "'" << std::endl;

      switch (solLocation) {
        case Vertex:
          definedOn = ResultInfo::NODE;
          break;
        case CellCenter:
          definedOn = ResultInfo::ELEMENT;
          break;
        default:
          if (!warnedGridLoc) {
            WARN("NACS currently does not support reading CGNS data defined on "
                << "entities other than nodes or elements. Data may be incomplete.");
            warnedGridLoc = true;
          }
          continue;
      }
      
      // Get number of physical fields
      CGNS_CHECK_EX( cg_nfields(fileHandle_, baseIndex_, zone, 1, &numFields) );
      
      // Loop over fields and check if we support it
      for (Integer field = 1; field <= numFields; ++field) {
        // Get field name
        CGNS_CHECK_EX( cg_field_info(fileHandle_, baseIndex_, zone, 1, field,
                                     &dataType, fieldName) );

        LOG_TRACE(simInputCGNS) << "field 1, name:" << fieldName << std::endl;
        
        // Convert field name to solution type
        solType = GetFieldType(fieldName, entryType);
        
        // Field is unsupported if we get NO_SOLUTION_TYPE
        if (solType != NO_SOLUTION_TYPE) {
          resultsFound = true;
          
          // Create a new ResultInfo struct
          shared_ptr<ResultInfo> resInfo(new ResultInfo());
          resInfo->resultType = solType;
          resInfo->resultName = SolutionTypeEnum.ToString(solType);
          resInfo->unit = MapSolTypeToUnit(solType);
          resInfo->definedOn = definedOn;
          resInfo->entryType = entryType;
          if (entryType == ResultInfo::VECTOR) {
            resInfo->SetVectorDOFs(dim_, false);
          }
          else {
            resInfo->dofNames = "";
          }
          
          for (res = 0; res < results_.size(); ++res) {
            if (*resInfo == *results_[res]) {
              break;
            }
          }
          if (res == results_.size()) {
            // Add to list of results only if it doesn't exist already
            results_.push_back(resInfo);
          }
          else {
            // We have to reset the variable, so we don't generate duplicates
            resInfo = results_[res];
          }
          
          // Add this zone to the list of entities for this result
          resultEntities_[resInfo].insert(zoneName);
          
          // Add step values for this result
          for (Integer s = 0; s < numSteps; ++s) {
            resultSteps_[resInfo].insert(steps[s]);
          }
        }
      }
    }
    
    // If there is no grid, there is nothing to do with this file
    if (!gridFound) {
      return;
    }
    
    // If there are no results, probably we have found the mesh file
    if (!resultsFound) {
      if (fileNames_.size() == 0) {
        // store the mesh file only if there are no other files
        fileNames_[0.0] = fileName;
      }
      return;
    }
    
    if (numSteps > 0) {
      if (analysis_ == BasePDE::NO_ANALYSIS) {
        // This clears the mesh file if it happened to be the first one found
        fileNames_.clear();
      }
      analysis_ = BasePDE::TRANSIENT;
      for (Integer s = 0; s < numSteps; ++s) {
        fileNames_[steps[s]] = fileName;
      }
    }
    else if (analysis_ == BasePDE::NO_ANALYSIS) {
      analysis_ = BasePDE::STATIC;
      fileNames_[0.0] = fileName;
    }
  }
  
  // Convert CGNS field name to NACS solution type
  SolutionType SimInputCGNS::GetFieldType(const char* field,
                                          ResultInfo::EntryType & entryType)
  {
    entryType = ResultInfo::SCALAR;
    if (strcmp(field, "Density") == 0) {
      return FLUIDMECH_DENSITY;
    }
    else if (strcmp(field, "Pressure") == 0) {
      return FLUIDMECH_PRESSURE;
    }
    else if (strcmp(field, "Temperature") == 0) {
      return HEAT_TEMPERATURE;
    }
    else if (strcmp(field, "VelocityX") == 0 ||
             strcmp(field, "VelocityY") == 0 ||
             strcmp(field, "VelocityZ") == 0)
    {
      entryType = ResultInfo::VECTOR;
      return FLUIDMECH_VELOCITY;
    }
    else if (strcmp(field, "VelocityNormal") == 0) {
      return MECH_NORMAL_VELOCITY;
    }
    else if (strcmp(field, "AcousticPressure") == 0) {
      return ACOU_PRESSURE;
    }
    else if (strcmp(field, "VelocitySound") == 0) {
      return ACOU_ELEM_SPEED_OF_SOUND;
    }
    else if (strcmp(field, "SoundIntensity") == 0) {
      return ACOU_NORMAL_INTENSITY;
    }
    else if (strcmp(field, "Voltage") == 0) {
      return ELEC_POTENTIAL;
    }
    else if (strcmp(field, "ElectricFieldX") == 0 ||
             strcmp(field, "ElectricFieldY") == 0 ||
             strcmp(field, "ElectricFieldZ") == 0)
    {
      entryType = ResultInfo::VECTOR;
      return ELEC_FIELD_INTENSITY;
    }
    else if (strcmp(field, "MagneticFieldX") == 0 ||
             strcmp(field, "MagneticFieldY") == 0 ||
             strcmp(field, "MagneticFieldZ") == 0)
    {
      entryType = ResultInfo::VECTOR;
      return MAG_FIELD_INTENSITY;
    }
    else if (strcmp(field, "CurrentDensityX") == 0 ||
             strcmp(field, "CurrentDensityY") == 0 ||
             strcmp(field, "CurrentDensityZ") == 0)
    {
      entryType = ResultInfo::VECTOR;
      return MAG_TOTAL_CURRENT_DENSITY;
    }
    else if (strcmp(field, "LorentzForceX") == 0 ||
             strcmp(field, "LorentzForceY") == 0 ||
             strcmp(field, "LorentzForceZ") == 0)
    {
      entryType = ResultInfo::VECTOR;
      return MAG_FORCE_LORENTZ;
    }
    else if (strcmp(field, "ForceX") == 0 ||
             strcmp(field, "ForceY") == 0 ||
             strcmp(field, "ForceZ") == 0)
    {
      entryType = ResultInfo::VECTOR;
      return MECH_FORCE;
    }
    else {
      return NO_SOLUTION_TYPE;
    }
  }
  
  void SimInputCGNS::GetFieldNames(SolutionType solType,
                                  StdVector<std::string> &fieldNames)
  {
    switch (solType) {
      case FLUIDMECH_DENSITY:
        fieldNames = "Density";
        break;
      case ACOU_PRESSURE:
        fieldNames = "AcousticPressure";
        break;
      case FLUIDMECH_PRESSURE:
        fieldNames = "Pressure";
        break;
      case HEAT_TEMPERATURE:
        fieldNames = "Temperature";
        break;
      case FLUIDMECH_VELOCITY:
        fieldNames = "VelocityX", "VelocityY", "VelocityZ";
        break;
      case MECH_NORMAL_VELOCITY:
        fieldNames = "VelocityNormal";
        break;
      case ACOU_ELEM_SPEED_OF_SOUND:
        fieldNames = "VelocitySound";
        break;
      case ACOU_NORMAL_INTENSITY:
        fieldNames = "SoundIntensity";
        break;
      case ELEC_POTENTIAL:
        fieldNames = "Voltage";
        break;
      case ELEC_FIELD_INTENSITY:
        fieldNames = "ElectricFieldX", "ElectricFieldY", "ElectricFieldZ";
        break;
      case MAG_FIELD_INTENSITY:
        fieldNames = "MagneticFieldX", "MagneticFieldY", "MagneticFieldZ";
        break;
      case MAG_TOTAL_CURRENT_DENSITY:
        fieldNames = "CurrentDensityX", "CurrentDensityY", "CurrentDensityZ";
        break;
      case MAG_FORCE_LORENTZ:
        fieldNames = "LorentzForceX", "LorentzForceY", "LorentzForceZ";
        break;
      case MECH_FORCE:
        fieldNames = "ForceX", "ForceY", "ForceZ";
        break;
      default:
        EXCEPTION("Result '" << SolutionTypeEnum.ToString(solType)
            << "' is not supported in CGNS files.");
    }
  }

  void SimInputCGNS::ReadMesh(Grid* ptGrid) 
  {
    LOG_TRACE(simInputCGNS) << "Entering SimInputCGNS::ReadGrid" << std::endl;

    mi_ = ptGrid;
    numNodes_ = 0;
    numElems_ = 0;
    nodeOffset_ = mi_->GetNumNodes();
    elemOffset_ = mi_->GetNumElems();
    elemToCellIdx_.Resize(numZones_);
    
    ZoneType_t zoneType;

    for (Integer zone = 1; zone <= numZones_; ++zone) {
      CGNS_CHECK_EX( cg_zone_type(fileHandle_, baseIndex_, zone, &zoneType) );
      
      switch (zoneType) {
        case Unstructured:
          ReadUnstructuredZone(baseIndex_, zone);
          break;
        default:
          EXCEPTION("Zone type " << zoneType << " is unsupported.");
          break;
      }
    }
  }
  
  void SimInputCGNS::ReadUnstructuredZone(Integer base, Integer zone) {
    // Read zone name and size
    char zoneName[STRBUFSIZE];
    cgsize_t zoneSize[3];
    
    CGNS_CHECK_EX( cg_zone_read(fileHandle_, base, zone, zoneName, zoneSize) );
    
    Integer nVertex = zoneSize[0];
    
    // How many coordinate datasets are there?
    Integer nGrids = 0;

    CGNS_CHECK_EX( cg_ngrids(fileHandle_, base, zone, &nGrids) );
    if (nGrids == 0) {
      WARN("CGNS file '" << fileName_ << "':\nZone '" << zoneName
          << "' contains no grid.");
      return;
    }
    
    // Find the main coordinate dataset (there may be other coordinate datasets
    // due to grid motion or refinement)
    Integer gridIdx = 0;
    char gridCoordName[STRBUFSIZE];

    for (gridIdx = 1; gridIdx <= nGrids; ++gridIdx) {
      CGNS_CHECK_EX( cg_grid_read(fileHandle_, base, zone, gridIdx, gridCoordName) );
      if (strncmp(gridCoordName, "GridCoordinates", STRBUFSIZE-1) == 0) {
        break;
      }
    }
    if (gridIdx > nGrids) {
      gridIdx = 1;
    }
    

    //==================================================================
    // READ COORDINATES
    //==================================================================

    // How many coordinate arrays are there?
    Integer nCoords = 0;
    
    CGNS_CHECK_EX( cg_ncoords(fileHandle_, base, zone, &nCoords) );
    
    if (nCoords != (Integer) dim_){
      EXCEPTION("Physical grid dimension " << dim_ << " does not match " <<
                "number of coordinate arrays: " << nCoords); 
    }

    // temporary storage for the coordinates
    Matrix<Double> coordMat(nCoords, nVertex);

    StdVector<std::string> coordNames;
    coordNames = "CoordinateX", "CoordinateY", "CoordinateZ";
    cgsize_t range_min[3] = {1, 1, 1};
    cgsize_t range_max[3] = {nVertex, 1, 1};

    for (Integer dof = 0; dof < nCoords ; ++dof) {
      //read in coordinates
      CGNS_CHECK_EX(cg_coord_read(fileHandle_, base, zone, coordNames[dof].c_str(),
                    RealDouble, range_min, range_max, (void*) coordMat[dof]));
    }

    // Add nodes to grid
    UInt nodeOffset = mi_->GetNumNodes();
    mi_->AddNodes(nVertex);
    
    Vector<Double> nodeCoord(nCoords);
    for (Integer node = 0; node < nVertex; ++node) {
      nodeCoord[0] = coordMat[0][node];
      nodeCoord[1] = coordMat[1][node];
      if (nCoords == 3) {
        nodeCoord[2] = coordMat[2][node];
      }
      mi_->SetNodeCoordinate(nodeOffset+node+1, nodeCoord);
    }

    numNodes_ += nVertex;
    
    
    //==================================================================
    // READ CONNECTIVITY 
    //==================================================================
    
    // Prepare element to cell number map 
    std::unordered_map<UInt, cgsize_t> & elem2cell = elemToCellIdx_[zone-1];
    elem2cell.reserve(zoneSize[1]);
    
    // How many connectivity datasets are there?
    Integer nSections = 0;
    CGNS_CHECK_EX( cg_nsections(fileHandle_, base, zone, &nSections) );
    LOG_TRACE(simInputCGNS) << "Found " << nSections
                            <<" (Surface) Regions in the Grid"
                            << std::endl;

    // Variables used for reading section info
    char sectionName[STRBUFSIZE];
    ElementType_t eType;
    cgsize_t start, end;
    Integer nBoundary, parentFlag;
    
    // First loop over all section and see if we have to deal with polyhedra
    bool hasPolyhedra = false;
    numPolyFaces_ = 0;
    cgsize_t numPolygons = 0;
    StdVector<Integer> gonSec, faceSec, stdSec;
    
    for (Integer sec = 1; sec <= nSections; ++sec) {
      CGNS_CHECK_EX(cg_section_read(fileHandle_, base, zone, sec, sectionName,
                                    &eType, &start, &end, &nBoundary, &parentFlag));
      switch (eType) {
        case NFACE_n:
          hasPolyhedra = true;
          faceSec.push_back(sec);
          break;
        case NGON_n:
          numPolygons += end - start + 1;
          if (end - start + 1 > numPolyFaces_) {
            numPolyFaces_ = end - start + 1;
          }
          gonSec.push_back(sec);
          break;
        case ElementTypeNull:
        case ElementTypeUserDefined:
          WARN("CGNS file '" << fileName_ << "':\nSkipping element section '"
              << sectionName << "' because of unsupported element type.");
          break;
        default:
          stdSec.push_back(sec);
          break;
      }
    }
    
    // Now determine the order in which to read the datasets:
    // First normal sections, then polygons, then polyhedra
    StdVector<Integer> secOrder;
    secOrder.Reserve(nSections);
    std::copy(stdSec.begin(), stdSec.end(), std::back_inserter(secOrder));
    std::copy(gonSec.begin(), gonSec.end(), std::back_inserter(secOrder));
    std::copy(faceSec.begin(), faceSec.end(), std::back_inserter(secOrder));
    
    // prepare map for polygon connectivity
    std::unordered_map<cgsize_t, StdVector< StdVector<UInt> > > polyConn;
    polyConn.reserve(numPolygons);
    
    // Variables used for reading element connectivity
    cgsize_t elemArraySize;
    Integer numElems, vertsPerElem;
    StdVector<cgsize_t> elems;
    StdVector<cgsize_t> connectOffset;    
    RegionIdType regionId;
    UInt elemOffset;
    Elem::FEType feType;
    StdVector<UInt> connect, groupElems;
    StdVector< StdVector<UInt> > polyElems;
    std::unordered_set<UInt> sectionNodes;
    
    for (Integer sec = 0; sec < nSections; ++sec) {
      CGNS_CHECK_EX( cg_section_read(fileHandle_, base, zone, secOrder[sec],
                                     sectionName, &eType, &start, &end,
                                     &nBoundary, &parentFlag) );
      nameToZoneIdx_[sectionName] = zone;
      
      CGNS_CHECK_EX( cg_ElementDataSize(fileHandle_, base, zone, secOrder[sec],
                                        &elemArraySize) );

      numElems = end - start + 1;
      numElemsPerName_[sectionName] = numElems;
      cellRange_[sectionName] = std::make_pair(start, end);
      sectionNodes.clear();
      
      LOG_TRACE(simInputCGNS) << "Reading section " << sectionName << std::endl;
      LOG_TRACE(simInputCGNS) << "------------------------------------------------------" << std::endl;
      LOG_TRACE(simInputCGNS) << "Element Type: ";
      PrintElementType(eType);
      LOG_TRACE(simInputCGNS) << std::endl;
      LOG_TRACE(simInputCGNS) << "Number of elements: " << numElems << std::endl;
      LOG_TRACE(simInputCGNS) << "Element Index Range: " << start << " to " << end << std::endl;
      LOG_TRACE(simInputCGNS) << "Boundary number: " << nBoundary << std::endl;
      LOG_TRACE(simInputCGNS) << "ParentFlag: " << parentFlag << std::endl << std::endl;
      
      // Resize connectivity and offset buffer
      elems.Resize(elemArraySize);
      connectOffset.Resize(elemArraySize);
      
      // Ignore parental Data field for now and give NULL to function
      CGNS_CHECK_EX( cg_poly_elements_read(fileHandle_, base, zone, secOrder[sec],
                                      elems.GetPointer(), connectOffset.GetPointer(), NULL) );

      if (Elem::shapes[elemTypeMap_[eType]].dim == dim_) {
        regionId = mi_->AddRegion(std::string(sectionName));
      }
      else {
        regionId = mi_->AddRegion(std::string(sectionName));
        groupElems.Resize(numElems);
      }
      elemOffset = mi_->GetNumElems();
      if (eType != NGON_n || !hasPolyhedra || numElems != numPolyFaces_) {
        mi_->AddElems((UInt)numElems);
        numElems_ += numElems;
      }
      
      if (eType == MIXED) {
        // In this special case, the element connectivity array
        // also stores the element type within it,
        // so we loop over and determine the type on the fly.

        UInt eTypeIdx = 0;

        for (Integer e = 0; e < numElems; ++e) {
          // First number is the element type
          eType = (ElementType_t) elems[eTypeIdx];
          
          // Get number of nodes for this element type
          CGNS_CHECK_EX( cg_npe(eType, &vertsPerElem) );

          // Convert to NACS element type and re-order connectivity
          feType = elemTypeMap_[eType];
          TranslateConnectivity(feType, &elems[eTypeIdx+1], connect, nodeOffset);

          // Save element data in the grid object
          mi_->SetElemData(elemOffset+e+1, feType, regionId, connect.GetPointer());
          elem2cell[elemOffset+e+1] = start+e;

          // If this is not a region, we need to add this element to a group
          if (regionId == NO_REGION_ID) {
            groupElems[e] = elemOffset + e + 1;
          }
          
          // Count elements per dimension
          ++numElemsOfDim_[Elem::shapes[feType].dim];
          
          // Add nodes of this element to set of section nodes
          sectionNodes.insert(connect.begin(), connect.end());
          
          // Calculate index of next element 
          eTypeIdx += vertsPerElem + 1;            
        }
      }
      else if (eType == NGON_n) {
        // In this special case, the elements are given as polygon with their
        // number of nodes, so we loop over and determine the type on the fly.

        UInt nNodesIdx = 0, elemNum = elemOffset + 1, numTrias;

        for (Integer e = 0; e < numElems; ++e) {
          // First number is the number of nodes
          vertsPerElem = elems[nNodesIdx];
          
          // Convert to NACS element type and re-order connectivity
          feType = PolygonToFE(vertsPerElem, &elems[nNodesIdx+1],
                               polyElems, nodeOffset);
          
          numTrias = polyElems.size();
          
          if (hasPolyhedra) {
            polyConn[start+e] = polyElems;
          }
          if (!hasPolyhedra || numElems != numPolyFaces_) {
            // If we generate more than one triangle, we need to adjust the
            // total number of elements
            if (numTrias > 1) {
              mi_->AddElems(numTrias - 1);
              numElems_ += numTrias - 1;
              numElemsPerName_[sectionName] += numTrias - 1;
            }
            numElemsOfDim_[1] += numTrias;

            for (UInt t = 0; t < numTrias; ++t) {
              // Save element data in the grid object
              mi_->SetElemData(elemNum, feType, regionId, polyElems[t].GetPointer());
              elem2cell[elemNum] = start+e;

              // If this is not a region, we need to add this element to a group
              if (regionId == NO_REGION_ID) {
                groupElems[e] = elemNum;
              }

              // Add nodes of this element to set of section nodes
              sectionNodes.insert(polyElems[t].begin(), polyElems[t].end());
              
              ++elemNum;
            }
          }
          
          // Calculate index of next element 
          nNodesIdx += vertsPerElem + 1;            
        }
      }
      else if (eType == NFACE_n) {
        // In this special case, the elements are given as polyhedra with their
        // defining faces, so we loop over and determine the type on the fly.

        UInt nFaceIdx = 0, numTets = 0;

        // estimate the number of tetrahedra produced by polyhedra
        for (Integer e = 0; e < numElems; ++e) {
          // First number is the number of faces
          vertsPerElem = elems[nFaceIdx];
          // Check if polyhedron needs to be tetrahedralized
          if (vertsPerElem == 5 || vertsPerElem == 7 || vertsPerElem > 9) {
            numTets += vertsPerElem - 3;
          }
          // Calculate index of next element 
          nFaceIdx += vertsPerElem + 1;            
        }

        // Reserve memory for additional elements due to polyhedra
        // This speeds up subsequent calls to mi_->AddElems a lot!
        mi_->ReserveElems(numTets);
        
        // Now actually read the elements
        nFaceIdx = 0;
        UInt elemNum = elemOffset + 1, numRejected = 0;

        for (Integer e = 0; e < numElems; ++e) {
          // First number is the number of faces
          vertsPerElem = elems[nFaceIdx];
          
          // Convert to NACS element type and re-order connectivity
          feType = PolyhedronToFE(vertsPerElem, &elems[nFaceIdx+1], polyConn,
                                  polyElems, nodeOffset);
          numTets = polyElems.size();
          
          if (feType != Elem::ET_UNDEF) {
            // If we generate more than one triangle, we need to adjust the
            // total number of elements
            if (numTets > 1) {
              mi_->AddElems(numTets - 1);
              numElems_ += numTets - 1;
              numElemsPerName_[sectionName] += numTets - 1;
            }
            numElemsOfDim_[2] += numTets;

            for (UInt t = 0; t < numTets; ++t) {
              // Save element data in the grid object
              mi_->SetElemData(elemNum, feType, regionId, polyElems[t].GetPointer());
              elem2cell[elemNum] = start+e;

              // If this is not a region, we need to add this element to a group
              if (regionId == NO_REGION_ID) {
                groupElems[e] = elemNum;
              }
  
              // Add nodes of this element to set of section nodes
              sectionNodes.insert(polyElems[t].begin(), polyElems[t].end());
              
              ++elemNum;
            }
          }
          else {
            ++numRejected; 
          }
          
          // Calculate index of next element 
          nFaceIdx += vertsPerElem + 1;            
        }
        
        if (numRejected > 0) {
          WARN("Rejected " << numRejected << " polyhedra due to invalid topology.");
        }
      }
      else {
        // Determine number of vertices for current element
        CGNS_CHECK_EX( cg_npe(eType, &vertsPerElem) );

        // Convert to NACS element type
        feType = elemTypeMap_[eType];

        for (Integer e = 0; e < numElems; ++e) {
          // Re-order connectivity
          TranslateConnectivity(feType, &elems[e*vertsPerElem], connect, nodeOffset);

          // Save element data in the grid object
          mi_->SetElemData(elemOffset+e+1, feType, regionId, connect.GetPointer());
          elem2cell[elemOffset+e+1] = start+e;

          // If this is not a region, we need to add this element to a group
          if (regionId == NO_REGION_ID) {
            groupElems[e] = elemOffset + e + 1;
          }

          // Count elements per dimension
          ++numElemsOfDim_[Elem::shapes[feType].dim];
          
          // Add nodes of this element to set of section nodes
          sectionNodes.insert(connect.begin(), connect.end());
        }
      }
      
      // If this is not a region add an element group
      if (regionId == NO_REGION_ID &&
          (eType != NGON_n || !hasPolyhedra || numElems != numPolyFaces_))
      {
        mi_->AddNamedElems(std::string(sectionName), groupElems);
      }
      numNodesPerName_[sectionName] = sectionNodes.size();
    }
  }
  

  void SimInputCGNS::PrintElementType(ElementType_t eType){
   
    ElementType_t test = eType;
    switch(test){
      case CGNSLIB_H::ElementTypeNull:
        LOG_TRACE(simInputCGNS) << "ElementTypeNull";
        break;
      case CGNSLIB_H::ElementTypeUserDefined:
        LOG_TRACE(simInputCGNS) << "ElementTypeUserDefined";
        break;
      case CGNSLIB_H::NODE:
        LOG_TRACE(simInputCGNS) << "NODE";
        break;
      case CGNSLIB_H::BAR_2:
        LOG_TRACE(simInputCGNS) << "BAR_2";
        break;
      case CGNSLIB_H::BAR_3:
        LOG_TRACE(simInputCGNS) << "BAR_3";
        break;
      case CGNSLIB_H::TRI_3:
        LOG_TRACE(simInputCGNS) << "TRI_3";
        break;
      case CGNSLIB_H::TRI_6:
        LOG_TRACE(simInputCGNS) << "TRI_6";
        break;
      case CGNSLIB_H::QUAD_4:
        LOG_TRACE(simInputCGNS) << "QUAD_4";
        break;
      case CGNSLIB_H::QUAD_8:
        LOG_TRACE(simInputCGNS) << "QUAD_8";
        break;
      case CGNSLIB_H::QUAD_9:
        LOG_TRACE(simInputCGNS) << "QUAD_9";
        break;
      case CGNSLIB_H::TETRA_4:
        LOG_TRACE(simInputCGNS) << "TETRA_4";
        break;
      case CGNSLIB_H::TETRA_10:
        LOG_TRACE(simInputCGNS) << "TETRA_10";
        break;
      case CGNSLIB_H::PYRA_5:
        LOG_TRACE(simInputCGNS) << "PYRA_5";
        break;
#if CGNS_VERSION > 2550
      case CGNSLIB_H::PYRA_13:
        LOG_TRACE(simInputCGNS) << "PYRA_13";
        break;
#endif
      case CGNSLIB_H::PYRA_14:
        LOG_TRACE(simInputCGNS) << "PYRA_14";
        break;
      case CGNSLIB_H::PENTA_6:
        LOG_TRACE(simInputCGNS) << "PENTA_6";
        break;
      case CGNSLIB_H::PENTA_15:
        LOG_TRACE(simInputCGNS) << "PENTA_15";
        break;
      case CGNSLIB_H::PENTA_18:
        LOG_TRACE(simInputCGNS) << "PENTA_18";
        break;
      case CGNSLIB_H::HEXA_8:
        LOG_TRACE(simInputCGNS) << "HEXA_8";
        break;
      case CGNSLIB_H::HEXA_20:
        LOG_TRACE(simInputCGNS) << "HEXA_20";
        break;
      case CGNSLIB_H::HEXA_27:
        LOG_TRACE(simInputCGNS) << "HEXA_27";
        break;
      case CGNSLIB_H::MIXED:
        LOG_TRACE(simInputCGNS) << "MIXED";
        break;
#if CGNS_VERSION > 2550
      case CGNSLIB_H::NGON_n:
        LOG_TRACE(simInputCGNS) << "NGON_n";
        break;
      case CGNSLIB_H::NFACE_n:
        LOG_TRACE(simInputCGNS) << "NFACE_n";
        break;
#endif
      default:
        LOG_TRACE(simInputCGNS) << "Unknown element Type detected";
        break;
    }
  }

  void SimInputCGNS::InitElemTypeMap(){
    if (!elemTypeMap_.empty()) return;
    
    elemTypeMap_[CGNSLIB_H::ElementTypeNull] = Elem::ET_UNDEF;
    elemTypeMap_[CGNSLIB_H::ElementTypeUserDefined] = Elem::ET_UNDEF;
    elemTypeMap_[CGNSLIB_H::NODE] = Elem::ET_POINT;
    elemTypeMap_[CGNSLIB_H::BAR_2] = Elem::ET_LINE2;
    elemTypeMap_[CGNSLIB_H::BAR_3] = Elem::ET_LINE3;
    elemTypeMap_[CGNSLIB_H::TRI_3] = Elem::ET_TRIA3;
    elemTypeMap_[CGNSLIB_H::TRI_6] = Elem::ET_TRIA6;
    elemTypeMap_[CGNSLIB_H::QUAD_4] = Elem::ET_QUAD4;
    elemTypeMap_[CGNSLIB_H::QUAD_8] = Elem::ET_QUAD8;
    elemTypeMap_[CGNSLIB_H::QUAD_9] = Elem::ET_QUAD9;
    elemTypeMap_[CGNSLIB_H::TETRA_4] = Elem::ET_TET4;
    elemTypeMap_[CGNSLIB_H::TETRA_10] = Elem::ET_TET10;
    elemTypeMap_[CGNSLIB_H::PYRA_5] = Elem::ET_PYRA5;

#if CGNS_VERSION > 2550
    elemTypeMap_[CGNSLIB_H::PYRA_13] = Elem::ET_PYRA13;
#endif
    elemTypeMap_[CGNSLIB_H::PYRA_14] = Elem::ET_PYRA14;
    elemTypeMap_[CGNSLIB_H::PENTA_6] = Elem::ET_WEDGE6;
    elemTypeMap_[CGNSLIB_H::PENTA_15] = Elem::ET_WEDGE15;
    elemTypeMap_[CGNSLIB_H::PENTA_18] = Elem::ET_WEDGE18;
    elemTypeMap_[CGNSLIB_H::HEXA_8] = Elem::ET_HEXA8;
    elemTypeMap_[CGNSLIB_H::HEXA_20] = Elem::ET_HEXA20;
    elemTypeMap_[CGNSLIB_H::HEXA_27] = Elem::ET_HEXA27;
    elemTypeMap_[CGNSLIB_H::MIXED] = Elem::ET_UNDEF;
#if CGNS_VERSION > 2550
    elemTypeMap_[CGNSLIB_H::NGON_n] = Elem::ET_UNDEF;
    elemTypeMap_[CGNSLIB_H::NFACE_n] = Elem::ET_UNDEF;
#endif
  }

  void SimInputCGNS::TranslateConnectivity(Elem::FEType feType,
                                           cgsize_t* cgnsConn,
                                           StdVector<UInt>& connect,
                                           UInt offset)
  {
    UInt numElemNodes = Elem::shapes[feType].numNodes;

    static const int trDefault[27] = {
      0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
      20, 21, 22, 23, 24, 25, 26
    };

    // Map from a CGNS HEXA_20 connectivity
    static const int trHEX20[20] = {
      0, 1, 2, 3,
      4, 5, 6, 7,
      8, 9, 10, 11,
      16, 17, 18, 19,
      12, 13, 14, 15
    };
    // Map from a CGNS HEXA_27 connectivity
    static const int trHEX27[27] = {
      0, 1, 2, 3,
      4, 5, 6, 7,
      8, 9, 10, 11,
      16, 17, 18, 19,
      12, 13, 14, 15,
      21, 22, 23, 24,
      20, 25,
      26
    };
    // Map from a CGNS PENTA_15 connectivity
    static const int trPRI15[15] = {
      0, 1, 2, 3, 4, 5, 6, 7, 8,
      12, 13, 14,
      9, 10, 11 
    };
    // Map from a CGNS PENTA_18 connectivity
    static const int trPRI18[18] = {
      0, 1, 2, 3, 4, 5,
      6, 7, 8, 12, 13, 14,
      9, 10, 11, 15, 16, 17
    };
    
    const int *tr;
    switch(feType) {
    case Elem::ET_HEXA20:
      tr = trHEX20;
      break;
    case Elem::ET_HEXA27:
      tr = trHEX27;
      break;
    case Elem::ET_WEDGE15:
      tr = trPRI15;
      break;
    case Elem::ET_WEDGE18:
      tr = trPRI18;
      break;
    default:
      tr = trDefault;
      break;
    }
    
    connect.Resize(numElemNodes);
    for (UInt n = 0; n < numElemNodes; ++n) {
      connect[n] = cgnsConn[tr[n]] + offset;
    }
  }

  // Try to guess the element type from its number of nodes
  Elem::FEType SimInputCGNS::PolygonToFE(cgsize_t numNodes, cgsize_t *nodes,
                                            StdVector< StdVector<UInt> > &elems,
                                            UInt offset)
  {
    Elem::FEType feType;
    
    switch (numNodes) {
      case 2:
        feType = Elem::ET_LINE2;
        break;
      case 4:
        feType = Elem::ET_QUAD4;
        break;
      default:
        feType = Elem::ET_TRIA3;
        break;
    }
    
    if (numNodes <= 4) {
      elems.Resize(1);
      elems[0].Resize(numNodes);
      for (Integer i = 0; i < numNodes; ++i) {
        elems[0][i] = (UInt) nodes[i] + offset;
      }
    }
    else {
      elems.Clear(true);
      Triangulate(numNodes, nodes, elems, offset);
    }
    
    return feType;
  }
  
  template<typename T>
  void SimInputCGNS::Triangulate(cgsize_t numNodes, T* nodes,
                                 StdVector< StdVector<UInt> > &elems,
                                 UInt offset)
  {
    UInt prevSize = elems.size();
    elems.Resize(prevSize + numNodes - 2);
    
    for (Integer e = 0; e < numNodes - 2; ++e) {
      elems[prevSize+e].Resize(3);
      elems[prevSize+e][0] = nodes[0] + offset;
      elems[prevSize+e][1] = nodes[e+1] + offset;
      elems[prevSize+e][2] = nodes[e+2] + offset;
      assert(elems[prevSize+e].Find(offset) == -1);
    }
  }
  
  // Convert a polyhedron into a proper element
  Elem::FEType SimInputCGNS::PolyhedronToFE(cgsize_t numFaces, cgsize_t *faces,
      std::unordered_map<cgsize_t, StdVector< StdVector<UInt> > > & polyConn,
      StdVector< StdVector<UInt> > &connect, UInt offset)
  {
    Elem::FEType feType = Elem::ET_UNDEF;
    
    // Determine types of faces
    UInt numTrias = 0, numQuads = 0, numPolys = 0;
    StdVector<Elem::FEType> faceTypes(numFaces);
    
    for (Integer f = 0; f < numFaces; ++f) {
      if (polyConn.find(abs(faces[f])) == polyConn.end()) {
#ifndef RELEASE_TO_CUSTOMER
        EXCEPTION("Could not find polyhedron face #" << abs(faces[f]));
#endif
        return feType;
      }
      
      // If a face is a true polygon, we can skip the rest of the faces
      if (polyConn[abs(faces[f])].size() > 1) {
        ++numPolys;
        break;
      }
      
      switch (polyConn[abs(faces[f])][0].size()) {
        case 3:
          ++numTrias;
          faceTypes[f] = Elem::ET_TRIA3;
          break;
        case 4:
          ++numQuads;
          faceTypes[f] = Elem::ET_QUAD4;
          break;
        case 6:
          ++numTrias;
          faceTypes[f] = Elem::ET_TRIA6;
          //order2 = true;
          break;
        case 8:
          ++numQuads;
          faceTypes[f] = Elem::ET_QUAD8;
          //order2 = true;
          break;
        case 9:
          ++numQuads;
          faceTypes[f] = Elem::ET_QUAD9;
          //order2 = true;
          break;
        default:
          ++numPolys;
          break;
      }
    }
    
    connect.Clear(true);
    
    // Simple check for supported volume elements
    if (numFaces >= 4 && numFaces <= 6 && numPolys == 0) {
      
      // Assume there will be one element
      connect.Resize(1);
      
      if (numTrias == 4 && numQuads == 0) {
        // =============
        //  TETRAHEDRON
        // =============
        // 4 triangles make a tetrahedron
        
        //feType = order2 ? Elem::ET_TET10 : Elem::ET_TET4;
        feType = Elem::ET_TET4;
        connect[0].Resize(4, 0);
        
        // take first face as the base
        connect[0][0] = polyConn[abs(faces[0])][0][0] + offset;
        connect[0][1] = polyConn[abs(faces[0])][0][1] + offset;
        connect[0][2] = polyConn[abs(faces[0])][0][2] + offset;
        
        // second face must contain 4th missing node
        for (Integer n = 0; n < 3; ++n) {
          if (connect[0].Find(polyConn[abs(faces[1])][0][n]) == -1) {
            connect[0][3] = polyConn[abs(faces[1])][0][n] + offset;
            break;
          }
        }
      }
      else if (numTrias == 4 && numQuads == 1) {
        // =========
        //  PYRAMID
        // =========
        // 4 triangles and 1 quadrilateral make a pyramid
        
        //feType = order2 ? Elem::ET_PYRA13 : Elem::ET_PYRA5;
        feType = Elem::ET_PYRA5;
        connect[0].Resize(5, 0);
        
        // take quadrilateral as the base
        for (Integer f = 0; f < numFaces; ++f) {
          if (polyConn[abs(faces[f])][0].size() == 4) {
            connect[0][0] = polyConn[abs(faces[f])][0][0] + offset;
            connect[0][1] = polyConn[abs(faces[f])][0][1] + offset;
            connect[0][2] = polyConn[abs(faces[f])][0][2] + offset;
            connect[0][3] = polyConn[abs(faces[f])][0][3] + offset;
            break;
          }
        }
        
        // first triangle must contain the 5th missing node
        for (Integer f = 0; f < numFaces; ++f) {
          if (polyConn[abs(faces[f])][0].size() == 3) {
            for (Integer n = 0; n < 3; ++n) {
              if (connect[0].Find(polyConn[abs(faces[f])][0][n]) == -1) {
                connect[0][4] = polyConn[abs(faces[f])][0][n] + offset;
              }
            }
            break;
          }
        }
      }
      else if (numTrias == 2 && numQuads == 3) {
        // =======
        //  WEDGE
        // =======
        // 2 triangles and 3 quadrilaterals make a wegde
        
        //feType = order2 ? Elem::ET_WEDGE15 : Elem::ET_WEDGE6;
        feType = Elem::ET_WEDGE6;
        connect[0].Resize(6, 0);
        
        // look for the first triangle and use it as the base
        Integer face;
        for (face = 0; face < numFaces; ++face) {
          if (polyConn[abs(faces[face])][0].size() == 3) {
            break;
          }
        }
        StdVector<UInt> & base = polyConn[abs(faces[face])][0]; 
        connect[0][0] = base[0] + offset;
        connect[0][1] = base[1] + offset;
        connect[0][2] = base[2] + offset;
  
        // find the faces that share the first and second edges with the base
        Integer edge, node;
        for (edge = 0; edge < 2; ++edge) {
          for (face = 0; face < numFaces; ++face) {
            StdVector<UInt> & fcon = polyConn[abs(faces[face])][0];
            if (fcon.size() == 3) continue;
            
            for (node = 0; node < 4; ++node) {
              if (fcon[0+node] == base[0+edge] && fcon[(1+node)%4] == base[1+edge]) {
                connect[0][3+edge] = fcon[(3+node)%4] + offset;
                connect[0][4+edge] = fcon[(2+node)%4] + offset;
                break;
              }
              if (fcon[0+node] == base[1+edge] && fcon[(1+node)%4] == base[0+edge]) {
                connect[0][3+edge] = fcon[(2+node)%4] + offset;
                connect[0][4+edge] = fcon[(3+node)%4] + offset;
                break;
              }
            }
  
            if (connect[0][4+edge] > 0) break;
          }
        }
      }
      else if (numTrias == 0 && numQuads == 6) {
        // ============
        //  HEXAHEDRON
        // ============
        // 6 quadrilaterals make a hexahedron

        //feType = order2 ? Elem::ET_HEXA20 : Elem::ET_HEXA8;
        feType = Elem::ET_HEXA8;
        connect[0].Resize(8, 0);
        
        // take first face as the base
        StdVector<UInt> & base = polyConn[abs(faces[0])][0];
        connect[0][0] = base[0] + offset;
        connect[0][1] = base[1] + offset;
        connect[0][2] = base[2] + offset;
        connect[0][3] = base[3] + offset;
  
        // find the faces that share the first and third edges with the base
        Integer edge, face, node;
        for (edge = 0; edge < 3; edge += 2) {
          for (face = 1; face < numFaces; ++face) {
            StdVector<UInt> & fcon = polyConn[abs(faces[face])][0];
  
            for (node = 0; node < 4; ++node) {
              if (fcon[0+node] == base[0+edge] && fcon[(1+node)%4] == base[1+edge]) {
                connect[0][4+edge] = fcon[(3+node)%4] + offset;
                connect[0][5+edge] = fcon[(2+node)%4] + offset;
                break;
              }
              if (fcon[0+node] == base[1+edge] && fcon[(1+node)%4] == base[0+edge]) {
                connect[0][4+edge] = fcon[(2+node)%4] + offset;
                connect[0][5+edge] = fcon[(3+node)%4] + offset;
                break;
              }
            }
  
            if (connect[0][4+edge] > 0) break;
          }
        }
      }
    }
    
    if (feType == Elem::ET_UNDEF || connect[0].Find(0) != -1) {
      // ============
      //  POLYHEDRON
      // ============
      
      // Polyhedron will be split into tetrahedra
      feType = Elem::ET_TET4;
      
      // First divide all faces into triangles
      StdVector< StdVector<UInt> > trias;
      for (Integer f = 0; f < numFaces; ++f) {
        StdVector< StdVector<UInt> > & fcon = polyConn[abs(faces[f])];
        for (UInt i = 0; i < fcon.size(); ++i) {
          Triangulate(fcon[i].size(), fcon[i].GetPointer(), trias);
        }
      }

      UInt numTrias = trias.size();
      connect.Clear(true);
      connect.Reserve(numTrias);

      // Pick the first node n0 and create tetrahedra with all triangles
      // that do not contain n0.
      UInt numTets = 0;
      UInt n0 = trias[0][0];
      for (UInt t = 1; t < numTrias; ++t) {
        if (trias[t].Find(n0) == -1) {
          connect.Resize(numTets+1);
          connect[numTets].Resize(4);
          connect[numTets][0] = trias[t][0];
          connect[numTets][1] = trias[t][1];
          connect[numTets][2] = trias[t][2];
          connect[numTets][3] = n0;
          ++numTets;
        }
      }
    }
    
    return feType;
  }

  // =========================================================================
  //  GENERAL SOLUTION INFORMATION
  // =========================================================================

  UInt SimInputCGNS::GetNumElems(Integer dim) {
    switch (dim) {
      case -1:
        return numElems_;
      case 1:
      case 2:
      case 3:
        return numElemsOfDim_[dim-1];
      default:
        EXCEPTION("Invalid dimension: " << dim);
    }
  }
  
  void SimInputCGNS::GetNumNodesElemsByName(const std::string &name,
                                            UInt &numNodes, UInt &numElems)
  {
    if (numElemsPerName_.find(name) == numElemsPerName_.end()) {
      EXCEPTION("Unknown region or group '" << name << "'");
    }
    
    numNodes = numNodesPerName_[name];
    numElems = numElemsPerName_[name];
  }
  
  void SimInputCGNS::
  GetNumMultiSequenceSteps( std::map<UInt,BasePDE::AnalysisType>& analysis,
                            std::map<UInt, UInt>& numSteps,
                            bool isHistory )
  {
    analysis.clear();
    numSteps.clear();
    
    if (!isHistory && analysis_ != BasePDE::NO_ANALYSIS) {
      analysis[1] = analysis_;
      numSteps[1] = fileNames_.size();
    }
  }


  void SimInputCGNS::
  GetStepValues( UInt sequenceStep,
                 shared_ptr<ResultInfo> info,
                 std::map<UInt, Double>& steps,
                 bool isHistory )
  {
    if (sequenceStep != 1) {
      if (numBases_ == 1) {
        EXCEPTION("Analysis step no. " << sequenceStep << " is invalid." <<
            "There is only one analysis step in CGNS file '" << fileName_ << "'");
      }
      else {
        EXCEPTION("NACS currently supports only one analysis step per CGNS file.");
      }
    }

    steps.clear();
    if (isHistory) {
      return;
    }
    
    if (analysis_ == BasePDE::STATIC) {
      steps[1] = 0.0;
      return;
    }
    
    std::map< shared_ptr<ResultInfo>, std::set<Double> >::iterator
        resIt = resultSteps_.begin(),
        resEnd = resultSteps_.end();
    for ( ; resIt != resEnd; ++resIt) {
      if (*info == *(resIt->first)) {
        break;
      }
    }
    
    if (resIt == resEnd) {
      return;
    }
    
    std::map<UInt, Double>::iterator stepIt, stepEnd;
    std::set<Double>::iterator valIt = resIt->second.begin(),
                               valEnd = resIt->second.end();
    for ( ; valIt != valEnd; ++valIt) {
      for (stepIt = stepMap_.begin(), stepEnd = stepMap_.end();
           stepIt != stepEnd;
           ++stepIt)
      {
        if (*valIt == stepIt->second) {
          steps[stepIt->first] = *valIt;
          break;
        }
      }
    }
  }


  void SimInputCGNS::
  GetResultTypes( UInt sequenceStep,
                  StdVector<shared_ptr<ResultInfo> >& infos,
                  bool isHistory )
  {
    if (sequenceStep != 1) {
      if (numBases_ == 1) {
        EXCEPTION("Analysis step no. " << sequenceStep << " is invalid." <<
            "There is only one analysis step in CGNS file '" << fileName_ << "'");
      }
      else {
        EXCEPTION("NACS currently supports only one analysis step per CGNS file.");
      }
    }
    
    infos.Clear(true);
    if (!isHistory) {
      infos = results_;
    }
  }

  void SimInputCGNS::GetResultEntityNames( UInt sequenceStep,
                                           shared_ptr<ResultInfo> info,
                                           StdVector<std::string> &names,
                                           bool isHistory)
  {
    if (sequenceStep != 1) {
      if (numBases_ == 1) {
        EXCEPTION("Analysis step no. " << sequenceStep << " is invalid." <<
            "There is only one analysis step in CGNS file '" << fileName_ << "'");
      }
      else {
        EXCEPTION("NACS currently supports only one analysis step per CGNS file.");
      }
    }

    names.Clear();
    if (isHistory)  {
      return;
    }
    
    std::map< shared_ptr<ResultInfo>, std::set<std::string> >::iterator
        resIt = resultEntities_.begin(),
        resEnd = resultEntities_.end();
    for ( ; resIt != resEnd; ++resIt) {
      if (*info == *(resIt->first)) {
        break;
      }
    }
    
    if (resIt == resEnd) {
      return;
    }
    
    std::set<std::string> &nameSet = resIt->second;
    names.Reserve(nameSet.size());
    std::copy(nameSet.begin(), nameSet.end(), std::back_inserter(names));
  }

  void SimInputCGNS::
  GetResultEntities( UInt sequenceStep,
                     shared_ptr<ResultInfo> info,
                     StdVector<shared_ptr<EntityList> >& list,
                     bool isHistory )
  {
    if (sequenceStep != 1) {
      if (numBases_ == 1) {
        EXCEPTION("Analysis step no. " << sequenceStep << " is invalid." <<
            "There is only one analysis step in CGNS file '" << fileName_ << "'");
      }
      else {
        EXCEPTION("NACS currently supports only one analysis step per CGNS file.");
      }
    }

    list.Clear();
    if (isHistory)  {
      return;
    }
    
    std::map< shared_ptr<ResultInfo>, std::set<std::string> >::iterator
        resIt = resultEntities_.begin(),
        resEnd = resultEntities_.end();
    for ( ; resIt != resEnd; ++resIt) {
      if (*info == *(resIt->first)) {
        break;
      }
    }
    
    if (resIt == resEnd) {
      return;
    }
    
    std::set<std::string> &names = resIt->second;
    UInt numRegions = names.size();
    list.Reserve(numRegions);
    EntityList::ListType listType = resIt->first->definedOn == ResultInfo::NODE ?
                                    EntityList::NODE_LIST : EntityList::ELEM_LIST;
    
    std::set<std::string>::iterator nameIt = names.begin(),
                                    nameEnd = names.end();
    for ( ; nameIt != nameEnd; ++nameIt) {
      list.push_back(mi_->GetEntityList(listType, *nameIt));
    }
  }

  void SimInputCGNS::GetResult( UInt sequenceStep,
                                UInt stepNum,
                                shared_ptr<BaseResult> result,
                                bool isHistory )
  {
    if (sequenceStep != 1) {
      if (numBases_ == 1) {
        EXCEPTION("Analysis step no. " << sequenceStep << " is invalid." <<
            "There is only one analysis step in CGNS file '" << fileName_ << "'");
      }
      else {
        EXCEPTION("NACS currently supports only one analysis step per CGNS file.");
      }
    }

    if (stepMap_.find(stepNum) == stepMap_.end()) {
      EXCEPTION("Step " << stepNum << " does not exist in CGNS file '"
          << fileName_ << "'.");
    }

    // Open the correct file for this step
    OpenFile(fileNames_[stepMap_[stepNum]]);
    
    shared_ptr<ResultInfo> resInfo = result->GetResultInfo();
    
    // determine zone for this results
    std::string regionName = result->GetEntityList()->GetName();
    if (nameToZoneIdx_.find(regionName) == nameToZoneIdx_.end()) {
      EXCEPTION("Entity '" << regionName << "' not found in CGNS file '"
          << fileName_ << "'.");
    }
    cgsize_t zone = nameToZoneIdx_[regionName];
    
    char solName[STRBUFSIZE];
    GridLocation_t solLoc;
    
    // Check that requested result lives on the same location
    CGNS_CHECK_EX( cg_sol_info(fileHandle_, baseIndex_, zone, 1, solName , &solLoc) );
    
    if ( (resInfo->definedOn == ResultInfo::ELEMENT && solLoc != CellCenter) ||
         (resInfo->definedOn == ResultInfo::NODE && solLoc != Vertex) )
    {
      EXCEPTION("Cannot read result '"
          << SolutionTypeEnum.ToString(resInfo->resultType)
          << "' from CGNS file '" << fileName_
          << "':\nresult is defined per cell and cannot be read per vertex.");
    }

    // Prepare the result vector
    UInt numDofs = resInfo->dofNames.size();
    UInt numEntities = result->GetEntityList()->GetSize();
    UInt resVecSize =  numEntities * numDofs;
    assert(numEntities == numElemsPerName_[regionName]);
    
    Vector<Double> & resVec = dynamic_cast<Result<Double>& >(*result).GetVector();
    resVec.Resize(resVecSize);
    resVec.Init();

    // Check if there is only a subset of the zone stored in the file
    PointSetType_t ptSetType;
    cgsize_t nPoints, ptRange[2];
    CGNS_CHECK_EX( cg_sol_ptset_info(fileHandle_, baseIndex_, zone, 1,
                                     &ptSetType, &nPoints) );
    if (ptSetType == PointRange) {
      CGNS_CHECK_EX( cg_sol_ptset_read(fileHandle_, baseIndex_, zone, 1, ptRange) );
    }
    else if (ptSetType == PointSetTypeNull) {
      ptRange[0] = cellRange_[regionName].first;
      ptRange[1] = cellRange_[regionName].second;
    }
    else {
      EXCEPTION("NACS currently does not support storage schemes in CGNS other "
          << "than full zones or point ranges.");
    }
    
    // Now we can set the index range needed for the current region
    cgsize_t range_min[3] = {cellRange_[regionName].first, 1, 1};
    cgsize_t range_max[3] = {cellRange_[regionName].second, 1, 1};
    
    // If our range appears to be shifted from the range in the CGNS file, that
    // is probably due to StarCCM+ not counting polygons in the cell numbering.
    // So we try to fix that here.
    UInt cellOffset = 0;
    if (range_min[0] > ptRange[0] && range_max[0] > ptRange[1]) {
      if (range_min[0] - numPolyFaces_ > 0) {
        range_min[0] -= numPolyFaces_;
        range_max[0] -= numPolyFaces_;
        cellOffset = numPolyFaces_;
      }
    }
    
    // Make sure that we don't attempt to read more data than available
    if (range_min[0] < ptRange[0]) {
      range_min[0] = ptRange[0];
    }
    if (range_max[0] > ptRange[1]) {
      range_max[0] = ptRange[1];
    }
    assert(range_min[0] < range_max[0]);
    Vector<Double> curSol(range_max[0] - range_min[0] + 1);
    
    // Get the field names for the requested result
    StdVector<std::string> fieldNames;
    GetFieldNames(resInfo->resultType, fieldNames);
    
    for (UInt f = 0, numFields = fieldNames.size(); f < numFields; ++f) {

      if ( cg_field_read(fileHandle_, baseIndex_, zone, 1, fieldNames[f].c_str(),
                         RealDouble, range_min, range_max,
                         (void *)curSol.GetPointer()) != CG_OK )
      {
        WARN("Error reading dataset '" << fieldNames[f] << "' in file '"
            << fileName_ << "':\n" << cg_get_error() << "\nData will be incomplete.");
        continue;
      }

      EntityIterator it = result->GetEntityList()->GetIterator();        
      for (UInt i = 0; i < numEntities; ++i, it++) {
        int idx = elemToCellIdx_[zone-1][it.GetElem()->elemNum - elemOffset_] - cellOffset;
        if (idx < range_min[0] || idx > range_max[0]) continue;
        resVec[i * numDofs + f] = curSol[idx - range_min[0]];
      }        
    }
  }

}
