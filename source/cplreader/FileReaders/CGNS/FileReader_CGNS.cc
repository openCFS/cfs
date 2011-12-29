// =====================================================================================
// 
//       Filename:  FileReader_CGNS.cc
// 
//    Description:  This is to read in a CGNS File
// 
//        Version:  1.0
//        Created:  02/22/2011 12:50:24 PM
//       Revision:  none
//       Compiler:  g++
// 
//         Author:  Andreas Hueppe (AHU), andreas.hueppe@uni-klu.ac.at
//        Company:  Universitaet Klagenfurt
// 
// =====================================================================================


#include "FileReader_CGNS.hh"
#include "cplreader/Settings.hh"

#include <sstream>
#include <iostream>
#include <vector>

#include "boost/regex.hpp"
#include "boost/filesystem/operations.hpp"
#include "boost/filesystem/path.hpp"
#include "boost/filesystem/convenience.hpp"
#include "boost/filesystem/exception.hpp"
namespace fs=boost::filesystem;
#include "boost/algorithm/string/predicate.hpp"
namespace algo=boost::algorithm;

namespace CoupledField{

  //! Constructor
  FileReader_CGNS::FileReader_CGNS(  const std::string& name,
                                     const UInt dim,
                                     const UInt numFiles) :
    FileReader(name, dim, numFiles, 1),
    numElems_(0){

     // std::cout << name_ << std::endl;
     // std::cout << baseName_ << std::endl;

    Settings& settings = Settings::Instance();
    //TODO Check if directory is valid!!
    fileDir_ = settings.GetString("basedir") + "/" + name_;
    gridInitialized_ = false;
    maxNumElemNodes_ = 0;


  }
  
  FileReader_CGNS::~FileReader_CGNS(){
  }
  
  void FileReader_CGNS::Init(){
    std::cout << "Entering FileReader_CGNS::Init" << std::endl;
    InitElemTypeMap();

    //std::string cgnsFileName = "v-rho-5/v-rho-5.000030.cgns";
    //std::cout << cgnsFileName << std::endl;
    std::string curFName;

    ReadCGNSDirectory(fileDir_, fileNames_);
    //=============================================
    //Read the mesh information from the first file
    //============================================
    ReadGrid();

    //exit(0);
  }
  
  //! get node coordinates from the corresponding file
  void FileReader_CGNS::ReadNodalCoords(std::vector<Double> & NODECOORD){
    //always 3D
    NODECOORD.resize((numVertices_)*3,0);

    for(UInt coord = 0;coord < nodeCoords_.GetSize();coord++){
      for(UInt node = 0;node < nodeCoords_[coord].GetSize();node++){
        NODECOORD[(3*node)+coord] = nodeCoords_[coord][node];
      }
    }
    return;
  }

  void FileReader_CGNS::ReadGrid(){
    std::cout << "Entering FileReader_CGNS::ReadGrid" << std::endl;

     std::string firstFile = fileDir_ + "/" + fileNames_.begin()->second;
     //open the first cgns file
     Integer fn = GetFileHandle(firstFile);
     //Check if the file is what we expect
     //i.e. 1 base, 1 zone, 1 grid
     CheckFileValidity(fn);

     char firstBaseName[200];
     char firstZoneName[200];
     int vertSize[9];
     Integer dim = 0;
     Integer physDim = 0;
     ZoneType_t gridType;

     //memset(vertSize,0,9*sizeof(cgsize_t));

     cg_base_read(fn,1,firstBaseName , &dim , &physDim );
     cg_zone_type(fn,1,1,&gridType);
     cg_zone_read(fn,1,1,firstZoneName,(cgsize_t*)vertSize);

     switch(gridType)
     {
      case ZoneTypeNull:
          std::cout << "Meshes with ZoneTypeNull are not supported." << std::endl;
          exit(1);
          break;
      case ZoneTypeUserDefined:
          std::cout << "Meshes with ZoneTypeUserDefinied are not supported." << std::endl;
          exit(1);
          break;
      case Structured:
          std::cout << "Meshes with Structured Meshes are not supported." << std::endl;
          exit(1);
          break;
      case Unstructured:
          std::cout << "Going to read " << vertSize[0] << " Vertices, building " << vertSize[1] << " Cells ";
          std::cout << "of a " << dim << "D, unstructured mesh ...." << std::endl;
          ReadUnstructuredGrid(fn,physDim,(cgsize_t*)vertSize);

          break;
     }
     CalcNumNodesPerRegion();
     cg_close(fn);
  }
  
  
  //! get topology information from the corresponding topology file
  void FileReader_CGNS::ReadTopology(std::vector<UInt> & TOPOLOGYDATA,
                                     std::vector<UInt> & elemTypes){
    TOPOLOGYDATA.resize(numElems_ * maxNumElemNodes_,0);
    elemTypes.resize(numElems_,0);
    std::map<Integer,StdVector<CGNSElem> >::iterator regIter = elemRegionMap_.begin();
    while(regIter != elemRegionMap_.end()){
      StdVector<CGNSElem> regElems = regIter->second;
      for(UInt curE = 0;curE<regElems.GetSize();curE++){
        Integer eIdx = regElems[curE].elemNum-1;
        elemTypes[eIdx] = regElems[curE].eType;
        for(UInt con = 0 ; con< regElems[curE].connect.GetSize(); con++){
          TOPOLOGYDATA[(eIdx*maxNumElemNodes_)+con] = regElems[curE].connect[con];
        }
      }
      regIter++;
    }
  }
  
  //! get nodal values from the corresponding fluid datafile the new way
  void FileReader_CGNS::ReadNodalValues(std::vector<FlowDataType>& nodalFlowData,
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
     solution.resize(3,std::vector<Double>(numVertices_));
     solutionSkinFriction.resize(3,std::vector<Double>(numVertices_));

     cg_nfields(fn, 1, 1, numsols, &nfields ); 
     bool velocityFound = false;
     bool pressureFound = false;
     bool frictionFound = false;
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
        }

        if(strcmp(fieldName,"Pressure") == 0){
          cgsize_t range_min[3] = {1,1,1};
          cgsize_t range_max[3] = {numVertices_,1,1};
          Double * curSol = new Double[numVertices_];
          cg_field_read(fn,1,1,1, fieldName, RealDouble , range_min, range_max, (void *)curSol );
          solutionPressure.resize(numVertices_,0);
          solutionPressure.assign(curSol,curSol +numVertices_);
          pressureFound = true;
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

       FlowDataPartStruct * tmpSolStruct = new FlowDataPartStruct();

       cg_section_read(fn,1,1,curReg+1,sectionName,&eType,&start,&end,&nboundary,&parentFlag);
       //do the following only for volume regions i.e. nboundary == 0

        if ( (requiredResults_[FLUIDMECH_VELOCITY] ||
             requiredResults_[NO_SOLUTION_TYPE]) && velocityFound ){
            /* copy the fluid velocity values */
          tmpSolStruct = &curType[FLUIDMECH_VELOCITY];
          
          //          if(nboundary == 0){
            tmpSolStruct->isActive = true;
            if (tmpSolStruct->dofNames.empty()) {
              tmpSolStruct->unit = MapSolTypeToUnit(FLUIDMECH_VELOCITY);
              tmpSolStruct->dofNames.push_back("x");
              tmpSolStruct->dofNames.push_back("y");
              if(dim_ == 3){
                tmpSolStruct->dofNames.push_back("z");
              }
            
              tmpSolStruct->resultName = SolutionTypeEnum.ToString(FLUIDMECH_VELOCITY);
              tmpSolStruct->definedOn = ResultInfo::NODE;
              tmpSolStruct->entryType = ResultInfo::VECTOR;
            }
            std::set<UInt> curNodes = nodeIt->second;
            std::set<UInt>::iterator curNodeIt = curNodes.begin();
            UInt numDOFs = tmpSolStruct->dofNames.size();
            tmpSolStruct->data.resize(numDOFs * curNodes.size());
            for(UInt n = 0; n < curNodes.size() ;n++){
              UInt idx = *curNodeIt;
              for(UInt i = 0; i < dim_; i++){
                tmpSolStruct->data[n*numDOFs+i] = solution[i][idx-1];
              }
              curNodeIt++;
            }
//           }else{
//             tmpSolStruct->isActive = false;
//           }

        }

        if (  pressureFound ){
            /* copy the fluid velocity values */
          tmpSolStruct = &curType[FLUIDMECH_PRESSURE];
          
          //          if(nboundary == 0){
            tmpSolStruct->isActive = true;
            if (tmpSolStruct->dofNames.empty()) {
              tmpSolStruct->unit = MapSolTypeToUnit(FLUIDMECH_PRESSURE);
              tmpSolStruct->dofNames.resize(1);
            
              tmpSolStruct->resultName = SolutionTypeEnum.ToString(FLUIDMECH_PRESSURE);
              tmpSolStruct->definedOn = ResultInfo::NODE;
              tmpSolStruct->entryType = ResultInfo::SCALAR;
            }
            std::set<UInt> curNodes = nodeIt->second;
            std::set<UInt>::iterator curNodeIt = curNodes.begin();
            UInt numDOFs = tmpSolStruct->dofNames.size();
            tmpSolStruct->data.resize(numDOFs * curNodes.size());
            for(UInt n = 0; n < curNodes.size() ;n++){
              UInt idx = *curNodeIt;
              tmpSolStruct->data[n] = solutionPressure[idx-1];
              curNodeIt++;
            }
        }

        if ( frictionFound){
            /* copy the fluid velocity values */
          tmpSolStruct = &curType[FLUIDMECH_SKINFRICTION];
          
          //          if(nboundary == 0){
            tmpSolStruct->isActive = true;
            if (tmpSolStruct->dofNames.empty()) {
              tmpSolStruct->unit = MapSolTypeToUnit(FLUIDMECH_SKINFRICTION);
              tmpSolStruct->dofNames.push_back("x");
              tmpSolStruct->dofNames.push_back("y");
              if(dim_ == 3){
                tmpSolStruct->dofNames.push_back("z");
              }
            
              tmpSolStruct->resultName = SolutionTypeEnum.ToString(FLUIDMECH_SKINFRICTION);
              tmpSolStruct->definedOn = ResultInfo::NODE;
              tmpSolStruct->entryType = ResultInfo::VECTOR;
            }
            std::set<UInt> curNodes = nodeIt->second;
            std::set<UInt>::iterator curNodeIt = curNodes.begin();
            UInt numDOFs = tmpSolStruct->dofNames.size();
            tmpSolStruct->data.resize(numDOFs * curNodes.size());
            for(UInt n = 0; n < curNodes.size() ;n++){
              UInt idx = *curNodeIt;
              for(UInt i = 0; i < dim_; i++){
                tmpSolStruct->data[n*numDOFs+i] = solutionSkinFriction[i][idx-1];
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
  
  Double FileReader_CGNS::GetTimeStep(UInt stepNumber){
    Settings& settings = Settings::Instance();
    Double timestep = settings.GetDouble("timestep");
    if (timestep != -1)
    {
      return timestep * stepNumber;
    }
    return -1;
  }


  std::string FileReader_CGNS::GetRegionName(const UInt partitionIdx){
    std::map<Integer,std::string>::iterator rIter = regionIndexToNameMap_.find(partitionIdx+1);
    if(rIter == regionIndexToNameMap_.end()){
      std::cout << "requested unknown region index " << partitionIdx << "... going to exit" << std::endl;
      exit(1);
    }

    return rIter->second;
  }
  
  //! get user data from file reader
  void FileReader_CGNS::GetUserData(std::map<std::string, std::string>& userData){
  }
  
  void FileReader_CGNS::GetRegionElements(std::vector<UInt> & regionElements,
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

  void FileReader_CGNS::ReadCGNSDirectory(std::string dirname, std::map<Double, std::string> & fileNames){
     Settings& settings = Settings::Instance();
     

     fs::path trnDir( dirname );
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
     UInt justmesh = settings.GetInt("justmesh");
     UInt calcsrc = settings.GetInt("calcsrc");

     if(justmesh == 1){
       std::cout << "Going to convert a CGNS file to hdf5..." << std::endl;
        for ( fs::directory_iterator dir_itr( trnDir );
              dir_itr != end_iter;
              ++dir_itr ) {
           if ( !fs::is_directory( *dir_itr ) ) { 
             fn = dir_itr->leaf();
             if(algo::ends_with(fn, ".cgns")) {
               fileNames[++counter] = fn;
             }
           }
        }

       if(counter > 1){
         std::string filename;
         std::cout << "Found more than one CGNS file";
         if(counter > 10){
           std::cout << " please type in the filename (case sensitive): " << std::endl;
           std::cin >> filename;
         }else{
           //print out each file with its number
           std::map<Double, std::string>::iterator mIter = fileNames.begin();
           std::cout <<std::endl;
           while(mIter!=fileNames.end()){
             std::cout << mIter->first << ") \t" << mIter->second << std::endl;
             mIter++;
           }
           UInt fNum = 0;
           while(true){
            if(fNum<1 || fNum > counter){
              std::cout << std::endl << "Please give the number of the file you want to convert:";
              std::cin >> fNum;
            }else{
              break;
            }
           }
           filename = fileNames[fNum];
         }

         //check if the file is there
         std::string ftest = name_ + "/" + filename;
         if(std::fopen(ftest.c_str(),"r")==0){
           std::cerr << "File not found " << filename << std::endl;
           exit(0);
         }else{
           fileNames.clear();
           fileNames[0] = filename;
           std::cout << "converting file " << filename << std::endl;
         }
       }

     }else if(calcsrc == 1){
       std::cout << "Going to iterate through the CGNS directory to obtain all source files" << std::endl;
       //give a wanring that the user knows whats going on
       std::cout << "Assuming that the step numbers are the last thing in the filename before the extension.\n"<<
             "If this is not the case the algorithm will fail!" << std::endl;
       UInt firstStep = settings.GetInt("firststep");
       UInt numSteps = settings.GetInt("numsteps");
       for ( fs::directory_iterator dir_itr( trnDir );
                  dir_itr != end_iter;
                  ++dir_itr ) {
         if ( !fs::is_directory( *dir_itr ) ) {
           fn = dir_itr->leaf();
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
       for(UInt i = 1;i<=firstStep;i++){
         fileNames.erase(fileNames.begin());
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
       std::cerr << "Please specify either calcsrc ot justmesh option" << std::endl;
       exit(0);
     }

     numSteps_ = fileNames.size();
  }
  
  Integer FileReader_CGNS::GetFileHandle(std::string fName){
    Integer fn = -1;
    if(fs::exists(fName)){
//#ifdef CG_MODE_READ
    if(cg_open(fName.c_str(), CG_MODE_READ, &fn) != CG_OK)
//#else
     // Still support pre 2.5 CGNS
 //    if(cg_open(fName.c_str(), MODE_READ, &fn) != CG_OK)
//#endif
     {
        std::cerr << "File open failed: " << fName << "... Going to exit" << std::endl; 
        EXCEPTION(cg_get_error());
        exit(1);
     }
    }else{
      std::cerr << "File does not exists: " << fName << "... Going to exit" << std::endl; 
      exit(1);
    }
    //std::cout << "file opened" << std::endl;
    return fn;
  }
  void FileReader_CGNS::CheckFileValidity(Integer fileHandle){
     Integer nbases,nzones,ngrids = 0;

     cg_nbases(fileHandle, &nbases);
     if(nbases != 1){
       std::cout << "ERROR: Found " << nbases << " Bases in the dataset" << std::endl;
       std::cout << "Found invalid number of bases, expected 1... Aborting" << std::endl;
       cg_close(fileHandle);
       exit(1);
     }
     
     cg_nzones(fileHandle,nbases,&nzones);
     if(nzones != 1){
       std::cout << "ERROR: Found " << nzones << " Zones in the dataset" << std::endl;
       std::cout << "Found invalid number of  zones, expected 1... Aborting" << std::endl;
       cg_close(fileHandle);
       exit(1);
     }
     cg_ngrids(fileHandle, nbases, nzones, &ngrids );
     if(ngrids != 1){
       std::cout << "ERROR: Found " << ngrids << " Grids in the dataset" << std::endl;
       std::cout << "Found invalid number of  zones, expected 1... Aborting" << std::endl;
       cg_close(fileHandle);
       exit(1);
     }
  }
  UInt FileReader_CGNS::MapCoordinateIndex(char* coordName){
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
       std::cerr << "Coordinate name is not recognized. invalid file or unsupported gridType (e.g.cylindical coord system)?" << std::endl;
       exit(1);
    }
    return coordinateIndex;
  }

  UInt FileReader_CGNS::MapVelocityIndex(char* coordName){
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

  UInt FileReader_CGNS::MapFrictionIndex(char* coordName){
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

  void FileReader_CGNS::ReadUnstructuredGrid(Integer fileHandle,Integer dim,cgsize_t* size){
     char gridCoordName[33];
     char curCoordName[33];
     Integer ncoords = 0;
     DataType_t dataType;

     numVertices_ = size[0];

     cg_grid_read(fileHandle, 1, 1, 1, gridCoordName );
     cg_ncoords(fileHandle, 1, 1, &ncoords );

     if(ncoords != dim){
       std::cerr << "Physical grid dimension " << dim << " does not match number of coordinate arrays: " << ncoords << std::endl; 
       exit(1);
     }
     dim_ = dim;

     //==================================================================
     // READ IN COORDINATES
     //==================================================================
     cgsize_t range_min[3] = {1,1,1};
     cgsize_t range_max[3] = {numVertices_,1,1};

     nodeCoords_.Resize(ncoords);
     Double * curCoord = new Double[numVertices_];
     for(Integer i = 1; i<=ncoords ; i++){
       //get the first coordinate name
       cg_coord_info(fileHandle,1,1,i, &dataType , curCoordName );
       //map the index x y z
       UInt idx = MapCoordinateIndex(curCoordName);
       //read in coordinates
       cg_coord_read(fileHandle, 1, 1, curCoordName, RealDouble , range_min, range_max, (void*)curCoord );

       //copy them into the datastructure
       std::vector<Double> tmp;
       tmp.resize(numVertices_,0);
       tmp.assign(curCoord,curCoord+numVertices_);
       nodeCoords_[idx] = StdVector<Double>(tmp);
     }
     delete [] curCoord;

     //==================================================================
     // READ IN Topology 
     //==================================================================
     
     Integer nsections = 0;
     numElems_ = 0;
     cg_nsections(fileHandle,1,1,&nsections);
     std::cout << "Found " << nsections << " (Surface) Regions in the Grid" << std::endl;

     //clear some class variables
     regionIndexToNameMap_.clear();
     elemRegionMap_.clear(); 
     numNodesPerRegion_.resize(nsections);
     numElemsPerRegion_.resize(nsections);

     numRegions_ = nsections;

     for(Integer curReg = 1; curReg<=nsections;curReg++){
       char sectionName[33];
       ElementType_t eType;
       cgsize_t start,end = 0;
       Integer nboundary = 0;
       Integer parentFlag = 0;
       cgsize_t elemArraySize= 0;
       Integer numEs = 0;
       Integer vertsPerElem= 0;

       cg_section_read(fileHandle,1,1,curReg,sectionName,&eType,&start,&end,&nboundary,&parentFlag);

       //elemArraySize = numCells*NodesPerElem
       cg_ElementDataSize(fileHandle,1,1,curReg,&elemArraySize);
       cg_npe(eType,&vertsPerElem);

       numEs = elemArraySize / vertsPerElem;
       numElems_ += numEs;

       maxNumElemNodes_ = ((UInt)vertsPerElem>maxNumElemNodes_)? vertsPerElem : maxNumElemNodes_;

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
       cgsize_t * curElems = new cgsize_t[elemArraySize];

       regionIndexToNameMap_[curReg] = std::string(sectionName);
       elemRegionMap_[curReg] = StdVector<CGNSElem>(numEs);
       elemTypeToRegionMap_[curReg] = elemTypeMap_[eType];
       //Ignore parental Data field for now and give NULL to function
       cg_elements_read(fileHandle,1,1,curReg,curElems,NULL);
       numElemsPerRegion_[curReg-1] = numEs;

       for(Integer i = 0;i<numEs;i++){
         CGNSElem curElem;
         curElem.elemNum = start+i;
         curElem.eType = elemTypeMap_[eType];
         curElem.connect.Resize(vertsPerElem);
         curElem.connect.Init();
         for(Integer j = 0;j < vertsPerElem ; j++){
           curElem.connect[j] = (UInt)curElems[(i*vertsPerElem)+j];
         }
         elemRegionMap_[curReg][i] = curElem;
       }
       delete [] curElems;
     }
  }

  void FileReader_CGNS::PrintElementType(ElementType_t eType){
   
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
      default:
        std::cout << "Unknown element Type detected";
        break;
    }
  }
  void FileReader_CGNS::InitElemTypeMap(){
    elemTypeMap_.clear();
    elemTypeMap_[CGNSLIB_H::ElementTypeNull] = Elem::UNDEF;
    elemTypeMap_[CGNSLIB_H::ElementTypeUserDefined] = Elem::UNDEF;
    elemTypeMap_[CGNSLIB_H::NODE] = Elem::POINT;
    elemTypeMap_[CGNSLIB_H::BAR_2] = Elem::LINE2;
    elemTypeMap_[CGNSLIB_H::BAR_3] = Elem::LINE3;
    elemTypeMap_[CGNSLIB_H::TRI_3] = Elem::TRIA3;
    elemTypeMap_[CGNSLIB_H::TRI_6] = Elem::TRIA6;
    elemTypeMap_[CGNSLIB_H::QUAD_4] = Elem::QUAD4;
    elemTypeMap_[CGNSLIB_H::QUAD_8] = Elem::QUAD8;
    elemTypeMap_[CGNSLIB_H::QUAD_9] = Elem::QUAD9;
    elemTypeMap_[CGNSLIB_H::TETRA_4] = Elem::TET4;
    elemTypeMap_[CGNSLIB_H::TETRA_10] = Elem::TET10;
    elemTypeMap_[CGNSLIB_H::PYRA_5] = Elem::PYRA5;
    elemTypeMap_[CGNSLIB_H::PYRA_14] = Elem::PYRA13;
    elemTypeMap_[CGNSLIB_H::PENTA_6] = Elem::UNDEF;
    elemTypeMap_[CGNSLIB_H::PENTA_15] = Elem::UNDEF;
    elemTypeMap_[CGNSLIB_H::PENTA_18] = Elem::UNDEF;
    elemTypeMap_[CGNSLIB_H::HEXA_8] = Elem::HEXA8;
    elemTypeMap_[CGNSLIB_H::HEXA_20] = Elem::HEXA20;
    elemTypeMap_[CGNSLIB_H::HEXA_27] = Elem::HEXA27;
    elemTypeMap_[CGNSLIB_H::MIXED] = Elem::UNDEF;
    elemTypeMap_[CGNSLIB_H::NGON_n] = Elem::UNDEF;
  }

  void FileReader_CGNS::CalcNumNodesPerRegion(){
    std::map<Integer,StdVector<CGNSElem> >::iterator mIter = elemRegionMap_.begin();
    numNodesPerRegion_.resize(elemRegionMap_.size());
    while(mIter !=elemRegionMap_.end()){
      StdVector<CGNSElem> curElems = mIter->second;
      std::set<UInt> tmpset;
     
      for(UInt i = 0 ; i<curElems.GetSize();i++){
        StdVector<UInt> curConn = curElems[i].connect;
        tmpset.insert(curConn.Begin(),curConn.End());
      }
      std::set<UInt>::iterator it;

      numNodesPerRegion_[mIter->first-1] = tmpset.size();
      nodesPerRegionMap_[mIter->first] = tmpset;
      std::cout << "Calculated Number of nodes or region " << mIter->first;
      std::cout << "\t result is " << tmpset.size() << std::endl;

      mIter++;
    }
  }
}
