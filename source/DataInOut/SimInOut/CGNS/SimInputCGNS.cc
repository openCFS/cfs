// =====================================================================================
// 
//       Filename:  SimInputCGNS.cc
// 
//    Description:  This is to read in a CGNS File
// 
//        Version:  1.0
//        Created:  02/22/2011 12:50:24 PM
//       Revision:  none
//       Compiler:  g++
// 
//        Authors:  Simon Triebenbacher, simon.triebenbacher@tuwien.ac.at
//                  Andreas Hueppe (AHU), andreas.hueppe@uni-klu.ac.at
//        Company:  TU Wien, Universitaet Klagenfurt
// 
// =====================================================================================

#include <sstream>
#include <iostream>
#include <vector>
#include <string>

#include <boost/regex.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/convenience.hpp>
#include <boost/filesystem/exception.hpp>
namespace fs=boost::filesystem;
#include <boost/algorithm/string/predicate.hpp>
namespace algo=boost::algorithm;

#include "DataInOut/Logging/LogConfigurator.hh"
#include "SimInputCGNS.hh"

namespace CoupledField{

  // declare logging stream
  DECLARE_LOG(simInputCGNS)
  DEFINE_LOG(simInputCGNS, "simInput.CGNS")

  //! Constructor
  SimInputCGNS::SimInputCGNS( std::string fileName, PtrParamNode inputNode,
                              PtrParamNode infoNode ) :
    SimInput(fileName, inputNode, infoNode ),
    numElems_(0),
    maxNumElemNodes_(27),
    dim_(0),
    physDim_(0)
  {
    capabilities_.insert( SimInput::MESH );
    //    capabilities_.insert( SimInput::MESH_RESULTS );

    gridInitialized_ = false;
  }
  
  SimInputCGNS::~SimInputCGNS(){
  }
  
  void SimInputCGNS::InitModule(){
    std::cout << "Entering SimInputCGNS::Init" << std::endl;

    std::string baseName;
    try
    {
      fs::path fn = fs::system_complete(fileName_);
      fn.normalize();
      fileDir_ = fn.branch_path().string();
      baseName = (fs::change_extension( fn.leaf(), "" )).string();
      if(fs::extension(fn) == "")
      {
        fn = fs::change_extension( fn, ".cgns" );
      }
      fileName_ = fn.string();
    } catch (fs::filesystem_error& ex)
    {
      EXCEPTION("Received exception: " << ex.what());
      return;
    }

    LOG_TRACE(simInputCGNS) << "fileName_: " << fileName_;
    LOG_TRACE(simInputCGNS) << "fileDir_: " << fileDir_;
    LOG_TRACE(simInputCGNS) << "baseName: " << baseName;

    InitElemTypeMap();

    ReadCGNSDirectory(fileDir_, fileNames_);

    //open the first cgns file
    Integer fn = GetFileHandle(fileName_);
    //Check if the file is what we expect
    //i.e. 1 base, 1 zone, 1 grid
    CheckFileValidity(fn);

    char firstBaseName[200];
    char firstZoneName[200];
    ZoneType_t zoneType;
    
    cg_base_read(fn,1,firstBaseName , &dim_ , &physDim_ );
    cg_zone_type(fn,1,1,&zoneType);

    switch(zoneType)
    {
    case ZoneTypeNull:
      EXCEPTION("Meshes with ZoneTypeNull are not supported.");
      break;
    case ZoneTypeUserDefined:
      EXCEPTION("Meshes with ZoneTypeUserDefined are not supported.");
      break;
    case Structured:
      EXCEPTION("Meshes with Structured Meshes are not supported.");
      break;
    case Unstructured:
      break;
    default:
      EXCEPTION("Unknown zone type.");
      break;
    }

    std::fill(vertSize_, &vertSize_[8], 0);
    cg_zone_read(fn,1,1,firstZoneName,(cgsize_t*)vertSize_);
  }

  void SimInputCGNS::ReadMesh(Grid* ptGrid) 
  {
    std::cout << "Entering SimInputCGNS::ReadGrid" << std::endl;

    mi_ = ptGrid;

    Integer fn = GetFileHandle(fileName_);

    std::cout << "Going to read " << vertSize_[0] << " Vertices, building " <<
      vertSize_[1] << " Cells ";
    std::cout << "of a " << dim_ << "D, unstructured mesh ...." << std::endl;
    ReadUnstructuredGrid(fn,physDim_,(cgsize_t*)vertSize_);

    cg_close(fn);
  }
  
  
#if 0
  //! get nodal values from the corresponding fluid datafile the new way
  void SimInputCGNS::ReadNodalValues(std::vector<FlowDataType>& nodalFlowData,
                                        const std::vector<bool>& activeParts,
                                        const UInt timeStepIdx){

    //UNSTRUCTURED VERSION!!!!!
     Integer numsols = 0; 
     char solname[33];
     GridLocation_t solLocation;
     //get grip of first file name
     std::map<Double, std::string>::iterator fIter = fileNames_.begin();
     //advance iterator to current timestep
     for(UInt i =0 ; i <timeStepIdx ; i++){
       fIter++;
     }
     if(fIter == fileNames_.end()){
       std::cerr << "Request for invalid stepIndex" << std::endl;
       exit(1);
     }
     //std::cout << "Going to read values from file " << fIter->second << std::endl;

     std::string firstFile = fileDir_ + "/" + fIter->second;
     //open the cgns file
     Integer fn = GetFileHandle(firstFile);

     cg_nsols(fn,1,1,&numsols);
     //std::cout << "found " << numsols << " solutions" << std::endl;
     cg_sol_info(fn, 1, 1, numsols, solname , &solLocation );
     //std::cout << "found " << solname << " solution" << std::endl;
     if(solLocation != Vertex){
       std::cout << "solution is definied on the cell center. mapping to nodes not supported right now" << std::endl;
       exit(1);
     }
     Integer nfields = 0;
     char fieldName[33];
     DataType_t datatype;
     std::vector< std::vector<Double> > solution;
     std::vector<Double> solutionPressure;
     std::vector< std::vector<Double> > solutionSkinFriction;
     std::vector< std::vector<Double> > solutionForce;
     solution.resize(3,std::vector<Double>(numVertices_));
     solutionSkinFriction.resize(3,std::vector<Double>(numVertices_));
     solutionForce.resize(3,std::vector<Double>(numVertices_));

     cg_nfields(fn, 1, 1, numsols, &nfields ); 
     bool velocityFound = false;
     bool pressureFound = false;
     bool frictionFound = false;
     bool fluidForceFound = false;
     for(UInt aSol = 1 ; aSol <= (UInt) nfields; aSol++){
        cg_field_info(fn, 1, 1, numsols, aSol, &datatype , fieldName );
        UInt spaceIdx = MapVelocityIndex(fieldName);
        if(spaceIdx != 9999){
          cgsize_t range_min[3] = {1,1,1};
          cgsize_t range_max[3] = {numVertices_,1,1};
          Double * curSol = new Double[numVertices_];
          cg_field_read(fn,1,1,1, fieldName, RealDouble , range_min, range_max, (void *)curSol );
          solution[spaceIdx].resize(numVertices_,0);
          solution[spaceIdx].assign(curSol,curSol +numVertices_);
          velocityFound = true;
          delete curSol;
          curSol = NULL;
        }
        spaceIdx = MapFrictionIndex(fieldName);
        if(spaceIdx != 9999){
          cgsize_t range_min[3] = {1,1,1};
          cgsize_t range_max[3] = {numVertices_,1,1};
          Double * curSol = new Double[numVertices_];
          cg_field_read(fn,1,1,1, fieldName, RealDouble , range_min, range_max, (void *)curSol );
          solutionSkinFriction[spaceIdx].resize(numVertices_,0);
          solutionSkinFriction[spaceIdx].assign(curSol,curSol +numVertices_);
          frictionFound = true;
          delete curSol;
          curSol = NULL;
        }
        spaceIdx = MapForceIndex(fieldName);
        if(spaceIdx != 9999){
          cgsize_t range_min[3] = {1,1,1};
          cgsize_t range_max[3] = {numVertices_,1,1};
          Double * curSol = new Double[numVertices_];
          cg_field_read(fn,1,1,1, fieldName, RealDouble , range_min, range_max, (void *)curSol );
          solutionForce[spaceIdx].resize(numVertices_,0);
          solutionForce[spaceIdx].assign(curSol,curSol +numVertices_);
          fluidForceFound = true;
          delete curSol;
          curSol = NULL;
        }

        if(strcmp(fieldName,"Pressure") == 0){
          cgsize_t range_min[3] = {1,1,1};
          cgsize_t range_max[3] = {numVertices_,1,1};
          Double * curSol = new Double[numVertices_];
          cg_field_read(fn,1,1,1, fieldName, RealDouble , range_min, range_max, (void *)curSol );
          solutionPressure.resize(numVertices_,0);
          solutionPressure.assign(curSol,curSol +numVertices_);
          pressureFound = true;
          delete curSol;
          curSol = NULL;
        }
     }
     if(!velocityFound){
       WARN("Did not find velocity field this could cause errors!");
     }
     //Global solution vector created
     //now fill the stupid array
     //Obtain the information
     FlowDataType curFlowData;

     //If the cgns file has information about the solution it would have "arrays"
     //but we assume that this is not hte case
     char sectionName[33];
     ElementType_t eType;
     cgsize_t start,end = 0;
     Integer nboundary = 0;
     Integer parentFlag = 0;

     std::map<Integer,std::set<UInt> >::iterator nodeIt = nodesPerRegionMap_.begin();
     for(UInt curReg = 0;curReg < nodalFlowData.size(); curReg++ ){
       FlowDataType& curType = nodalFlowData[curReg];



       cg_section_read(fn,1,1,curReg+1,sectionName,&eType,&start,&end,&nboundary,&parentFlag);
       //do the following only for volume regions i.e. nboundary == 0

        if ( (requiredResults_[FLUIDMECH_VELOCITY] ||
             requiredResults_[NO_SOLUTION_TYPE]) && velocityFound ){
            /* copy the fluid velocity values */
          FlowDataPartStruct & tmpSolStruct = curType[FLUIDMECH_VELOCITY];
          
          //          if(nboundary == 0){
            tmpSolStruct.isActive = true;
            if (tmpSolStruct.dofNames.empty()) {
              tmpSolStruct.unit = MapSolTypeToUnit(FLUIDMECH_VELOCITY);
              tmpSolStruct.dofNames.push_back("x");
              tmpSolStruct.dofNames.push_back("y");
              if(dim_ == 3){
                tmpSolStruct.dofNames.push_back("z");
              }
            
              tmpSolStruct.resultName = SolutionTypeEnum.ToString(FLUIDMECH_VELOCITY);
              tmpSolStruct.definedOn = ResultInfo::NODE;
              tmpSolStruct.entryType = ResultInfo::VECTOR;
            }
            std::set<UInt> curNodes = nodeIt->second;
            std::set<UInt>::iterator curNodeIt = curNodes.begin();
            UInt numDOFs = tmpSolStruct.dofNames.size();
            tmpSolStruct.data.resize(numDOFs * curNodes.size());
            for(UInt n = 0; n < curNodes.size() ;n++){
              UInt idx = *curNodeIt;
              for(UInt i = 0; i < dim_; i++){
                tmpSolStruct.data[n*numDOFs+i] = solution[i][idx-1];
              }
              curNodeIt++;
            }
//           }else{
//             tmpSolStruct->isActive = false;
//           }

        }

        if (  pressureFound ){
            /* copy the fluid velocity values */
          FlowDataPartStruct & tmpSolStruct = curType[FLUIDMECH_PRESSURE];
          
          //          if(nboundary == 0){
            tmpSolStruct.isActive = true;
            if (tmpSolStruct.dofNames.empty()) {
              tmpSolStruct.unit = MapSolTypeToUnit(FLUIDMECH_PRESSURE);
              tmpSolStruct.dofNames.resize(1);
            
              tmpSolStruct.resultName = SolutionTypeEnum.ToString(FLUIDMECH_PRESSURE);
              tmpSolStruct.definedOn = ResultInfo::NODE;
              tmpSolStruct.entryType = ResultInfo::SCALAR;
            }
            std::set<UInt> curNodes = nodeIt->second;
            std::set<UInt>::iterator curNodeIt = curNodes.begin();
            UInt numDOFs = tmpSolStruct.dofNames.size();
            tmpSolStruct.data.resize(numDOFs * curNodes.size());
            for(UInt n = 0; n < curNodes.size() ;n++){
              UInt idx = *curNodeIt;
              tmpSolStruct.data[n] = solutionPressure[idx-1];
              curNodeIt++;
            }
        }

        if ( frictionFound){
            /* copy the fluid velocity values */
          FlowDataPartStruct & tmpSolStruct = curType[FLUIDMECH_SKINFRICTION];
          
          //          if(nboundary == 0){
            tmpSolStruct.isActive = true;
            if (tmpSolStruct.dofNames.empty()) {
              tmpSolStruct.unit = MapSolTypeToUnit(FLUIDMECH_SKINFRICTION);
              tmpSolStruct.dofNames.push_back("x");
              tmpSolStruct.dofNames.push_back("y");
              if(dim_ == 3){
                tmpSolStruct.dofNames.push_back("z");
              }
            
              tmpSolStruct.resultName = SolutionTypeEnum.ToString(FLUIDMECH_SKINFRICTION);
              tmpSolStruct.definedOn = ResultInfo::NODE;
              tmpSolStruct.entryType = ResultInfo::VECTOR;
            }
            std::set<UInt> curNodes = nodeIt->second;
            std::set<UInt>::iterator curNodeIt = curNodes.begin();
            UInt numDOFs = tmpSolStruct.dofNames.size();
            tmpSolStruct.data.resize(numDOFs * curNodes.size());
            for(UInt n = 0; n < curNodes.size() ;n++){
              UInt idx = *curNodeIt;
              for(UInt i = 0; i < dim_; i++){
                tmpSolStruct.data[n*numDOFs+i] = solutionSkinFriction[i][idx-1];
              }
              curNodeIt++;
            }
        }

        if ( fluidForceFound){
            /* copy the fluid velocity values */
          FlowDataPartStruct & tmpSolStruct = curType[FLUIDMECH_FORCE];

          //          if(nboundary == 0){
            tmpSolStruct.isActive = true;
            if (tmpSolStruct.dofNames.empty()) {
              tmpSolStruct.unit = MapSolTypeToUnit(FLUIDMECH_FORCE);
              tmpSolStruct.dofNames.push_back("x");
              tmpSolStruct.dofNames.push_back("y");
              if(dim_ == 3){
                tmpSolStruct.dofNames.push_back("z");
              }

              tmpSolStruct.resultName = SolutionTypeEnum.ToString(FLUIDMECH_FORCE);
              tmpSolStruct.definedOn = ResultInfo::NODE;
              tmpSolStruct.entryType = ResultInfo::VECTOR;
            }
            std::set<UInt> curNodes = nodeIt->second;
            std::set<UInt>::iterator curNodeIt = curNodes.begin();
            UInt numDOFs = tmpSolStruct.dofNames.size();
            tmpSolStruct.data.resize(numDOFs * curNodes.size());
            for(UInt n = 0; n < curNodes.size() ;n++){
              UInt idx = *curNodeIt;
              for(UInt i = 0; i < dim_; i++){
                tmpSolStruct.data[n*numDOFs+i] = solutionForce[i][idx-1];
              }
              curNodeIt++;
            }
        }

      nodeIt++;
      }
       
    // std::string firstFile = fileDir_ + "/" + fileNames_[timeStepIdx]->second;
    // //open the first cgns file
    // Integer fn = GetFileHandle(firstFile);
    cg_close(fn);
  }
  
  Double SimInputCGNS::GetTimeStep(UInt stepNumber){
    Double timestep = settings.GetDouble("timestep");
    if (timestep != -1)
    {
      return timestep * stepNumber;
    }
    return -1;
  }


  std::string SimInputCGNS::GetRegionName(const UInt partitionIdx){
    std::map<Integer,std::string>::iterator rIter = regionIndexToNameMap_.find(partitionIdx+1);
    if(rIter == regionIndexToNameMap_.end()){
      std::cout << "requested unknown region index " << partitionIdx << "... going to exit" << std::endl;
      exit(1);
    }

    return rIter->second;
  }
  
  //! get user data from file reader
  void SimInputCGNS::GetUserData(std::map<std::string, std::string>& userData){
  }
  
  void SimInputCGNS::GetRegionElements(std::vector<UInt> & regionElements,
                                          const UInt regionIdx){
    std::map<Integer,StdVector<CGNSElem> >::iterator mIter = elemRegionMap_.find(regionIdx+1);
    if(mIter == elemRegionMap_.end()){
      std::cout << "requested unknown region index while getting elmeents... going to exit" << std::endl;
      exit(1);
    }
    UInt numElems = mIter->second.GetSize();
    regionElements.resize(numElems);
    for(UInt iElem = 0; iElem < numElems; iElem++){
      regionElements[iElem] = mIter->second[iElem].elemNum;
    }
  }
#endif

  void SimInputCGNS::ReadCGNSDirectory(std::string dirname, std::map<Double, std::string> & fileNames){
     std::string::size_type lastSlashPos=dirname.find_last_not_of('/');


     if (lastSlashPos != std::string::npos)
       dirname = dirname.substr( 0, lastSlashPos+1 );


     fs::path trnDir( dirname );
     //here we need to strap

     std::cout << trnDir.generic_string() << std::endl;

     fs::directory_iterator end_iter;

     //clear all contents of the map
     fileNames.clear();

     std::set<UInt>::const_iterator it, end;
     // UInt stepNum; // TODO: Unused variable stepNum
     // stepNum = 0;
     std::string fn;
     

     //new idea
     // 1. first fill an array with all cgns files
     // 2. check for the justmesh option
     // 2a. if true, and i found multiple files, ask the user which one to convert
     // 2b. if false we can continue with the standard

     //make a one based file counter
     UInt counter = 0;
     UInt justmesh = 1;// settings.GetInt("justmesh");
     UInt calcsrc = 0;//settings.GetInt("calcsrc");

     if(justmesh == 1){
       std::cout << "Going to convert a CGNS file to hdf5..." << std::endl;
        for ( fs::directory_iterator dir_itr( trnDir );
              dir_itr != end_iter;
              ++dir_itr ) {
           if ( !fs::is_directory( *dir_itr ) ) { 
             fn = dir_itr->path().filename().string();
             if(algo::ends_with(fn, ".cgns")) {
               fileNames[++counter] = fn;
             }
           }
        }
     }else if(calcsrc == 1){
       std::cout << "Going to iterate through the CGNS directory to obtain all source files" << std::endl;
       //give a wanring that the user knows whats going on
       std::cout << "Assuming that the step numbers are the last thing in the filename before the extension.\n"<<
             "If this is not the case the algorithm will fail!" << std::endl;
       UInt firstStep = 1;//settings.GetInt("firststep");
       UInt numSteps = 1;//settings.GetInt("numsteps");
       for ( fs::directory_iterator dir_itr( trnDir );
                  dir_itr != end_iter;
                  ++dir_itr ) {
         if ( !fs::is_directory( *dir_itr ) ) {
           fn = dir_itr->path().filename().string();
           std::cout << fn << std::endl;
           boost::tokenizer<> tok(fn);
           if(algo::ends_with(fn, ".cgns")) {

             //for(boost::tokenizer<>::iterator beg=tok.begin(); beg!=tok.end();++beg){
              //ASSUME THE element before the last token to be the number
              std::vector<std::string> parts( tok.begin(), tok.end() ) ;
              try{
                std::string doub =  parts[parts.size()-3] + "." + parts[parts.size()-2];
                Double number = boost::lexical_cast< Double >( doub  );
                std::cout.precision(16);
                fileNames[number] = fn;
              }catch( const boost::bad_lexical_cast & ){
                ///std::cout << "Cannot cast to Double. Maybe your files have a different name convention than expected.";
                std::cerr << "Cannot cast to double. Trying with int..." <<  std::endl;
                try{
                  Double number = boost::lexical_cast< UInt >( parts[parts.size()-2]  );
                  fileNames[number] = fn;
                }catch( const boost::bad_lexical_cast & ){
                  std::cout<< "Cannot cast to integer. Maybe your files have a different name convention than expected.";
                }
              }
            }else{
              std::cout << "Ignoring File " << fn << std::endl;
            }
          }
       }
       //now erase the number of files according to startstep parameter
       for(UInt i = 1;i<firstStep;i++){
         fileNames.erase(fileNames.begin());
       }
       if(numSteps == 0){
         numSteps = fileNames.size();
       }
       if(numSteps<fileNames.size()){
         UInt diff = fileNames.size() - numSteps;
         for(UInt i = 1;i<=diff;i++){
          std::map<Double, std::string>::iterator fiter = fileNames.end();
          fileNames.erase(--fiter);
         }
       }
       
       std::cout << "Found " << fileNames.size() << " files/steps which appear to be valid" << std::endl;
       //std::map<Double, std::string>::iterator iter = fileNames.begin();
       //for(UInt i = 1;i<=fileNames.size();i++){
       //  std::cout << iter->first << " with name " << iter->second << std::endl;
       //  iter++;
       //}
     }else{
       std::cerr << "Please specify either calcsrc or justmesh option" << std::endl;
       exit(0);
     }

     numSteps_ = 1;//fileNames.size();
  }
  
  Integer SimInputCGNS::GetFileHandle(std::string fName){
    Integer fn = -1;
    if(fs::exists(fName)){
      //#ifdef CG_MODE_READ
      if(cg_open(fName.c_str(), CG_MODE_READ, &fn) != CG_OK)
        //#else
        // Still support pre 2.5 CGNS
        //    if(cg_open(fName.c_str(), MODE_READ, &fn) != CG_OK)
        //#endif
     {
       EXCEPTION("File open failed: " << fName << "..." << std::endl <<
                 cg_get_error());
     }
    } else {
      EXCEPTION("File does not exists: " << fName << "...");
    }
    return fn;
  }

  void SimInputCGNS::CheckFileValidity(Integer fileHandle){
     Integer nbases,nzones,ngrids = 0;

     cg_nbases(fileHandle, &nbases);
     if(nbases != 1){
       cg_close(fileHandle);
       EXCEPTION("Found " << nbases << " Bases in the dataset." << std::endl <<
                 "Found invalid number of bases, expected 1...");
     }
     
     cg_nzones(fileHandle,nbases,&nzones);
     if(nzones != 1){
       cg_close(fileHandle);
       EXCEPTION("Found " << nzones << " Zones in the dataset" << std::endl <<
                 "Found invalid number of  zones, expected 1...");
     }

     cg_ngrids(fileHandle, nbases, nzones, &ngrids );
     if(ngrids != 1){
       cg_close(fileHandle);
       EXCEPTION("Found " << ngrids << " Grids in the dataset" << std::endl <<
                 "Found invalid number of  zones, expected 1...");
     }
  }

  UInt SimInputCGNS::MapCoordinateIndex(char* coordName){
    UInt coordinateIndex = 9999;
    //check if the name fulfils the naming convention
    //if not, we throw an error
    if(strcmp(coordName,"CoordinateX") == 0){
       coordinateIndex = 0;
    }else if(strcmp(coordName,"CoordinateY") == 0){
       coordinateIndex = 1;
    }else if(strcmp(coordName,"CoordinateZ") == 0){
       coordinateIndex = 2;
    }else{
      EXCEPTION("Coordinate name is not recognized.\n" <<
                "Invalid file or unsupported gridType " <<
                "(e.g. cylindrical coord system)?");
    }
    return coordinateIndex;
  }

  UInt SimInputCGNS::MapVelocityIndex(char* coordName){
    UInt coordinateIndex = 9999;
    //check if the name fulfils the naming convention
    //if not, we throw an error
    if(strcmp(coordName,"VelocityX") == 0 || strcmp(coordName,"Velocity") == 0){
       coordinateIndex = 0;
    }else if(strcmp(coordName,"VelocityY") == 0){
       coordinateIndex = 1;
    }else if(strcmp(coordName,"VelocityZ") == 0){
       coordinateIndex = 2;
    }
    return coordinateIndex;
  }

  UInt SimInputCGNS::MapFrictionIndex(char* coordName){
    UInt coordinateIndex = 9999;
    //check if the name fulfils the naming convention
    //if not, we throw an error
    if(strcmp(coordName,"SkinFrictionX") == 0 || strcmp(coordName,"SkinFriction") == 0){
       coordinateIndex = 0;
    }else if(strcmp(coordName,"SkinFrictionY") == 0){
       coordinateIndex = 1;
    }else if(strcmp(coordName,"SkinFrictionZ") == 0){
       coordinateIndex = 2;
    }
    return coordinateIndex;
  }

  UInt SimInputCGNS::MapForceIndex(char* coordName){
    UInt coordinateIndex = 9999;
    //check if the name fulfils the naming convention
    //if not, we throw an error
    if(strcmp(coordName,"ForceX") == 0 || strcmp(coordName,"Force") == 0){
       coordinateIndex = 0;
    }else if(strcmp(coordName,"ForceY") == 0){
       coordinateIndex = 1;
    }else if(strcmp(coordName,"ForceZ") == 0){
       coordinateIndex = 2;
    }
    return coordinateIndex;
  }

  void SimInputCGNS::ReadUnstructuredGrid(Integer fileHandle, Integer dim, cgsize_t* size){
    char gridCoordName[33];
    char curCoordName[33];
    Integer ncoords = 0;
    DataType_t dataType;
    
    numVertices_ = size[0];
    
    cg_grid_read(fileHandle, 1, 1, 1, gridCoordName );
    cg_ncoords(fileHandle, 1, 1, &ncoords );
    
    if(ncoords != dim){
      EXCEPTION("Physical grid dimension " << dim << " does not match " <<
                "number of coordinate arrays: " << ncoords); 
    }
    dim_ = dim;

    //==================================================================
    // READ IN COORDINATES
    //==================================================================
    cgsize_t range_min[3] = {1,1,1};
    cgsize_t range_max[3] = {numVertices_,1,1};
    
    nodeCoords_.Resize(numVertices_);
    for(UInt j = 0; j<numVertices_; j++){
      nodeCoords_[j].Resize(dim_);
    }      

    Double * curCoord = new Double[numVertices_];
    for(Integer i = 1; i<=ncoords ; i++){
      //get the first coordinate name
      cg_coord_info(fileHandle,1,1,i, &dataType , curCoordName );
      //map the index x y z
      UInt idx = MapCoordinateIndex(curCoordName);
      //read in coordinates
      cg_coord_read(fileHandle, 1, 1, curCoordName, RealDouble , range_min, range_max, (void*)curCoord );
      
      for(UInt j = 0; j<numVertices_; j++){
        nodeCoords_[j][idx] = curCoord[j];
      }      
    }
    delete [] curCoord;

    // Add nodes to grid
    mi_->AddNodes(numVertices_);
    for(UInt i = 0; i<numVertices_; i++){
      mi_->SetNodeCoordinate( i+1, nodeCoords_[i] );
    }

    //==================================================================
    // READ IN Topology 
    //==================================================================
    
    Integer nsections = 0;
    numElems_ = 0;
    cg_nsections(fileHandle,1,1,&nsections);
    std::cout << "Found " << nsections << " (Surface) Regions in the Grid"
              << std::endl;
    
    for(Integer curReg = 1; curReg<=nsections; curReg++){
      char sectionName[33];
      ElementType_t eType;
      cgsize_t start,end = 0;
      Integer nboundary = 0;
      Integer parentFlag = 0;
      cgsize_t elemArraySize= 0;
      Integer numEs = 0;
      Integer vertsPerElem= 0;
      
      cg_section_read(fileHandle,1,1,curReg,sectionName,&eType,&start,&end,
                      &nboundary,&parentFlag);
      
      //elemArraySize = numCells*NodesPerElem
      cg_ElementDataSize(fileHandle,1,1,curReg,&elemArraySize);

      numEs = (end-start+1);

      if(eType==ElementTypeNull || eType==ElementTypeUserDefined) {
        EXCEPTION("CGNS reader cannot handle nullType or user defined " <<
                  "element types.");
      }
      
      numElems_ += numEs;

      std::cout << "Reading section " << sectionName << std::endl;
      std::cout << "------------------------------------------------------" << std::endl;
      std::cout << "Element Type: ";
      PrintElementType(eType);
      std::cout << std::endl;
      std::cout << "Number of elements: " << numEs << std::endl;
      std::cout << "Element Index Range: " << start << " to " << end << std::endl;
      std::cout << "Boundary number: " << nboundary << std::endl;
      std::cout << "ParentFlag: " << parentFlag << std::endl << std::endl;
      
      //create the elements array
      StdVector<cgsize_t> curElems(elemArraySize);

      //Ignore parental Data field for now and give NULL to function
      cg_elements_read(fileHandle,1,1,curReg,&curElems[0],NULL);

      UInt elemOffset = mi_->GetNumElems();
      StdVector<UInt> connect(1024);
      RegionIdType regionId = -2;

      switch(eType) {
      case MIXED:
        {
          //ok in this special case, the element connectivity array
          //also stores the element type within it
          //so we loop over and determine the type on the fly

          mi_->AddElems(numEs);          
          // Determine element type of first element in region.
          Elem::FEType feType = elemTypeMap_[(ElementType_t)curElems[0]];
          
          AddRegionToGrid(regionId, feType, sectionName);

          //should point to the index storing the current eType
          //afterwards there will be the node numbers
          UInt eTypeIdx = 0;
          for(Integer i = 0;i<numEs;i++){
            //read
            eType = (ElementType_t)curElems[eTypeIdx];
            cg_npe(eType,&vertsPerElem);
            
            feType = elemTypeMap_[eType];
            TranslateConnectivity(feType,
                                  &curElems[eTypeIdx+1],
                                  connect);
            
            mi_->SetElemData(elemOffset+i+1,
                             elemTypeMap_[eType],
                             regionId,
                             &connect[0]);

            eTypeIdx += vertsPerElem+1;            
          }
        }
        break;
      case NGON_n:
      case NFACE_n:
        // We do not know how to handle polygonal or polyhedral cells.
        break;
      case NODE:
        // Nodes are currently not handled as elements by CFS++.
        break;
      default:
        {
          mi_->AddElems(numEs);

          //determine number of vertices for current element
          cg_npe(eType,&vertsPerElem);

          // Determine element type of first element in region.
          Elem::FEType feType = elemTypeMap_[eType];

          AddRegionToGrid(regionId, feType, sectionName);

          for(Integer i = 0;i<numEs;i++){
            TranslateConnectivity(feType,
                                  &curElems[i*vertsPerElem],
                                  connect);
            
            mi_->SetElemData(elemOffset+i+1,
                             elemTypeMap_[eType],
                             regionId,
                             &connect[0]);
          }
        }
        break;        
      }
    }
  }
  
  void SimInputCGNS::AddRegionToGrid(RegionIdType& regionId,
                                     const Elem::FEType feType,
                                     const std::string& sectionName)
  {
    regionId = -2;
    UInt elemDim = Elem::shapes[feType].dim;

    switch(elemDim) 
    {
    case 1:
      if(dim_ == 2) 
      {              
        regionId = mi_->AddRegion(std::string(sectionName));
        mi_->regionData[regionId].type = Grid::SURFACE_REGION;
      }
      break;            
    case 2:
      if(dim_ == 2) 
      {              
        regionId = mi_->AddRegion(std::string(sectionName));
        mi_->regionData[regionId].type = Grid::VOLUME_REGION;
      }
      if(dim_ == 3) 
      {              
        regionId = mi_->AddRegion(std::string(sectionName));
        mi_->regionData[regionId].type = Grid::SURFACE_REGION;
      }
      break;            
    case 3:
      if(dim_ == 3) 
      {              
        regionId = mi_->AddRegion(std::string(sectionName));
        mi_->regionData[regionId].type = Grid::VOLUME_REGION;
      }
      break;
    }
    
    if(regionId < 0) 
    {
      EXCEPTION("Could not create region '" << sectionName << "' for " <<
                "specified dimension!\nDimension of grid is " << dim_ <<
                " and elements in region are of dimension " << elemDim << ".");
    }
  }
  

  void SimInputCGNS::PrintElementType(ElementType_t eType){
   
    ElementType_t test = eType;
    switch(test){
      case CGNSLIB_H::ElementTypeNull:
        std::cout << "ElementTypeNull";
        break;
      case CGNSLIB_H::ElementTypeUserDefined:
        std::cout << "ElementTypeUserDefined";
        break;
      case CGNSLIB_H::NODE:
        std::cout << "NODE";
        break;
      case CGNSLIB_H::BAR_2:
        std::cout << "BAR_2";
        break;
      case CGNSLIB_H::BAR_3:
        std::cout << "BAR_3";
        break;
      case CGNSLIB_H::TRI_3:
        std::cout << "TRI_3";
        break;
      case CGNSLIB_H::TRI_6:
        std::cout << "TRI_6";
        break;
      case CGNSLIB_H::QUAD_4:
        std::cout << "QUAD_4";
        break;
      case CGNSLIB_H::QUAD_8:
        std::cout << "QUAD_8";
        break;
      case CGNSLIB_H::QUAD_9:
        std::cout << "QUAD_9";
        break;
      case CGNSLIB_H::TETRA_4:
        std::cout << "TETRA_4";
        break;
      case CGNSLIB_H::TETRA_10:
        std::cout << "TETRA_10";
        break;
      case CGNSLIB_H::PYRA_5:
        std::cout << "PYRA_5";
        break;
      case CGNSLIB_H::PYRA_13:
        std::cout << "PYRA_13";
        break;
      case CGNSLIB_H::PYRA_14:
        std::cout << "PYRA_14";
        break;
      case CGNSLIB_H::PENTA_6:
        std::cout << "PENTA_6";
        break;
      case CGNSLIB_H::PENTA_15:
        std::cout << "PENTA_15";
        break;
      case CGNSLIB_H::PENTA_18:
        std::cout << "PENTA_18";
        break;
      case CGNSLIB_H::HEXA_8:
        std::cout << "HEXA_8";
        break;
      case CGNSLIB_H::HEXA_20:
        std::cout << "HEXA_20";
        break;
      case CGNSLIB_H::HEXA_27:
        std::cout << "HEXA_27";
        break;
      case CGNSLIB_H::MIXED:
        std::cout << "MIXED";
        break;
      case CGNSLIB_H::NGON_n:
        std::cout << "NGON_n";
        break;
      case CGNSLIB_H::NFACE_n:
        std::cout << "NFACE_n";
        break;
      default:
        std::cout << "Unknown element Type detected";
        break;
    }
  }

  void SimInputCGNS::InitElemTypeMap(){
    elemTypeMap_.clear();
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
    elemTypeMap_[CGNSLIB_H::PYRA_13] = Elem::ET_PYRA13;
    elemTypeMap_[CGNSLIB_H::PYRA_14] = Elem::ET_PYRA14;
    elemTypeMap_[CGNSLIB_H::PENTA_6] = Elem::ET_WEDGE6;
    elemTypeMap_[CGNSLIB_H::PENTA_15] = Elem::ET_WEDGE15;
    elemTypeMap_[CGNSLIB_H::PENTA_18] = Elem::ET_WEDGE18;
    elemTypeMap_[CGNSLIB_H::HEXA_8] = Elem::ET_HEXA8;
    elemTypeMap_[CGNSLIB_H::HEXA_20] = Elem::ET_HEXA20;
    elemTypeMap_[CGNSLIB_H::HEXA_27] = Elem::ET_HEXA27;
    elemTypeMap_[CGNSLIB_H::MIXED] = Elem::ET_UNDEF;
    elemTypeMap_[CGNSLIB_H::NGON_n] = Elem::ET_UNDEF;
  }

  void SimInputCGNS::TranslateConnectivity(Elem::FEType feType,
                                           cgsize_t* cgnsConn,
                                           StdVector<UInt>& connect)
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

    for(UInt n = 0; n<numElemNodes; n++ ) {
      connect[n] = cgnsConn[tr[n]];
    }
  }
}
