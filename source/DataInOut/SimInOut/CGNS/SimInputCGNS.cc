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
  SimInputCGNS::SimInputCGNS( std::string fileName,
                              PtrParamNode inputNode,
                              PtrParamNode infoNode ) :
    SimInput(fileName, inputNode, infoNode ),
    numElems_(0),
    numRegions_(0),
    dim_(0),
    physDim_(0),
    nodeOffset_(0),
    cgVersion_(-1)
  {
    capabilities_.insert( SimInput::MESH );
    capabilities_.insert( SimInput::MESH_RESULTS );
  }
  
  SimInputCGNS::~SimInputCGNS(){
  }
  
  void SimInputCGNS::InitModule(){
    LOG_TRACE(simInputCGNS) << "Entering SimInputCGNS::Init" << std::endl;

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


    // ReadCGNSDirectory(fileDir_, fileNames_);

    //open the first cgns file
    Integer fn = GetFileHandle(fileName_);

    // Determine CGNS version used to write the file.
    cg_version(fn, &cgVersion_);

    InitElemTypeMap();

    //Check if the file is what we expect
    //i.e. 1 base, 1 zone, 1 grid
    CheckFileValidity(fn);

    char firstBaseName[200];
    char firstZoneName[200];
    ZoneType_t zoneType;
    
    cg_base_read(fn,1,firstBaseName , &dim_ , &physDim_ );
    cg_zone_type(fn,1,1,&zoneType);

//    std::cout << "\n Name: " << firstBaseName << "   readDIM: " <<
//    		     dim_ << "  PhysDIM: " << physDim_ << std::endl;
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
    cg_zone_read(fn,1,1,firstZoneName,vertSize_);
  }

  void SimInputCGNS::ReadMesh(Grid* ptGrid) 
  {
    LOG_TRACE(simInputCGNS) << "Entering SimInputCGNS::ReadGrid" << std::endl;

    mi_ = ptGrid;

    Integer fn = GetFileHandle(fileName_);

    LOG_TRACE(simInputCGNS) << "Going to read " << vertSize_[0]
                            << " Vertices, building " <<  vertSize_[1]
                            << " Cells " << "of a " << dim_
                            << "D, unstructured mesh ...." << std::endl;

    ReadUnstructuredGrid(fn,physDim_,(cgsize_t*)vertSize_);

    cg_close(fn);
  }
  
  void SimInputCGNS::ReadCGNSDirectory(std::string dirname,
                                       std::map<Double, std::string> & fileNames)
  {
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
    
    if(justmesh == 1) {
      LOG_TRACE(simInputCGNS) << "Going to read a single CGNS file..." << std::endl;
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

    //==================================================================
    // READ IN COORDINATES
    //==================================================================
    cgsize_t range_min[3] = {1,1,1};
    cgsize_t range_max[3] = {cgsize_t(numVertices_),1,1};

    //storing the coordinates
    StdVector< Vector<Double> > nodeCoords;
    
    nodeCoords.Resize(numVertices_);
    for(UInt j = 0; j<numVertices_; j++){
      nodeCoords[j].Resize(dim_);
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
        nodeCoords[j][idx] = curCoord[j];
      }      
    }
    delete [] curCoord;

    // Add nodes to grid
    nodeOffset_ = mi_->GetNumNodes();
    mi_->AddNodes(numVertices_);
    for(UInt i = 0; i<numVertices_; i++){
      mi_->SetNodeCoordinate( i+1, nodeCoords[i] );
    }

    //==================================================================
    // READ IN Topology 
    //==================================================================
    
    Integer nsections = 0;
    numElems_ = 0;
    cg_nsections(fileHandle,1,1,&nsections);
    LOG_TRACE(simInputCGNS) << "Found " << nsections
                            <<" (Surface) Regions in the Grid"
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

      LOG_TRACE(simInputCGNS) << "Reading section " << sectionName << std::endl;
      LOG_TRACE(simInputCGNS) << "------------------------------------------------------" << std::endl;
      LOG_TRACE(simInputCGNS) << "Element Type: ";
      PrintElementType(eType);
      LOG_TRACE(simInputCGNS) << std::endl;
      LOG_TRACE(simInputCGNS) << "Number of elements: " << numEs << std::endl;
      LOG_TRACE(simInputCGNS) << "Element Index Range: " << start << " to " << end << std::endl;
      LOG_TRACE(simInputCGNS) << "Boundary number: " << nboundary << std::endl;
      LOG_TRACE(simInputCGNS) << "ParentFlag: " << parentFlag << std::endl << std::endl;
      
      //create the elements array
      StdVector<cgsize_t> curElems(elemArraySize);

      //Ignore parental Data field for now and give NULL to function
      cg_elements_read(fileHandle,1,1,curReg,&curElems[0],NULL);

      UInt elemOffset = mi_->GetNumElems();
      StdVector<UInt> connect(1024);
      RegionIdType regionId = -2;

      switch(eType) 
      {
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
#if CGNS_VERSION > 2550
      case NGON_n:
      case NFACE_n:
        // Do nothing!!
        break;
#endif
      default:
        {
          mi_->AddElems(numEs);
          
          //determine number of vertices for current element
          cg_npe(eType,&vertsPerElem);

          Elem::FEType feType = elemTypeMap_[eType];
          
          AddRegionToGrid(regionId, feType, sectionName);
          
          for(Integer i = 0;i<numEs;i++){
            TranslateConnectivity(feType,
                                  &curElems[i*vertsPerElem],
                                  connect);
          
            mi_->SetElemData(elemOffset+i+1,
                             feType,
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

//    std::cout << "\nDIM: " << dim_ << "  elemDim: " << elemDim
//    		  << "   sectionName: " << sectionName << std::endl;

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

    numRegions_++;
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

  // =========================================================================
  //  GENERAL SOLUTION INFORMATION
  // =========================================================================

  void SimInputCGNS::
  GetNumMultiSequenceSteps( std::map<UInt,BasePDE::AnalysisType>& analysis,
                            std::map<UInt, UInt>& numSteps,
                            bool isHistory )
  {
    analysis.clear();
    numSteps.clear();
    
    // At the moment we just can read RANS solutions
    analysis[1] = BasePDE::STATIC;
    numSteps[1] = 1;
  }


  void SimInputCGNS::
  GetStepValues( UInt sequenceStep,
                 shared_ptr<ResultInfo> info,
                 std::map<UInt, Double>& steps,
                 bool isHistory )
  {
    if(sequenceStep != 1) 
    {
      EXCEPTION("Only one static (RANS) sequence step supported!\n" <<
                "Note, that this is no limitation of CGNS but " <<
                "of the current implementation.");
    }
    
    steps[1] = 0.0;
  }


  void SimInputCGNS::
  GetResultTypes( UInt sequenceStep,
                  StdVector<shared_ptr<ResultInfo> >& infos,
                  bool isHistory )
  {
    if(sequenceStep != 1) 
    {
      EXCEPTION("Only one static (RANS) sequence step supported!\n" <<
                "Note, that this is no limitation of CGNS but " <<
                "of the current implementation.");
    }

    //UNSTRUCTURED VERSION!!!!!
    Integer numsols = 0; 
    char solname[33];
    GridLocation_t solLocation;
    
    //open the cgns file
    Integer fn = GetFileHandle(fileName_);
    
    cg_nsols(fn,1,1,&numsols);
    if(!numsols){
      cg_close(fn);
      return;
    }
    
    LOG_TRACE(simInputCGNS) << "found " << numsols << " solutions" << std::endl;
    cg_sol_info(fn, 1, 1, numsols, solname , &solLocation );
    LOG_TRACE(simInputCGNS) << "found " << solname << " solution" << std::endl;
    if(solLocation != Vertex){
      cg_close(fn);
      WARN("Solution is defined on the cell center.\n" <<
           "Mapping to nodes not supported right now.");
      return;
    }

    Integer nfields = 0;
    char fieldName[33];
    DataType_t datatype;

    bool velocityFound = false;
    bool pressureFound = false;

    cg_nfields(fn, 1, 1, numsols, &nfields ); 
    for(UInt aSol = 1 ; aSol <= (UInt) nfields; aSol++){
      cg_field_info(fn, 1, 1, numsols, aSol, &datatype , fieldName );

      LOG_TRACE(simInputCGNS) << "field name " << fieldName << std::endl;

      UInt spaceIdx = MapVelocityIndex(fieldName);
      if(spaceIdx != 9999){
        velocityFound = true;
      }
      
      if(strcmp(fieldName,"Pressure") == 0){
        pressureFound = true;
      }      
    }

    cg_close(fn);

    if(velocityFound) 
    {
      // Check for subType
      StdVector<std::string> velDofNames;
      if ( dim_ == 3 ) {
        velDofNames = "x", "y", "z";
      } else { 
        velDofNames = "x", "y";
      }

      shared_ptr<ResultInfo> velocity( new ResultInfo );
      velocity->resultType = FLUIDMECH_VELOCITY;
      velocity->resultName = SolutionTypeEnum.ToString(velocity->resultType);
      velocity->dofNames = velDofNames;
      velocity->unit = MapSolTypeToUnit(velocity->resultType);
      
      velocity->definedOn = ResultInfo::NODE;
      velocity->entryType = ResultInfo::VECTOR;
      infos.Push_back( velocity );
    }

    if(pressureFound) 
    {
      shared_ptr<ResultInfo> pressure( new ResultInfo );
      pressure->resultType = FLUIDMECH_PRESSURE;
      pressure->resultName = SolutionTypeEnum.ToString(pressure->resultType);
      pressure->dofNames = "";
      pressure->unit = MapSolTypeToUnit(pressure->resultType);
      
      pressure->definedOn = ResultInfo::NODE;
      pressure->entryType = ResultInfo::SCALAR;
      infos.Push_back( pressure );
    }    
  }

  void SimInputCGNS::
  GetResultEntities( UInt sequenceStep,
                     shared_ptr<ResultInfo> info,
                     StdVector<shared_ptr<EntityList> >& list,
                     bool isHistory )
  {
    if(sequenceStep != 1) 
    {
      EXCEPTION("Only one static (RANS) sequence step supported!");
    }

    for(UInt i=0, n=mi_->regionData.GetSize(); i<n; i++)
    {
      if(mi_->regionData[i].type == Grid::VOLUME_REGION) 
      {
        list.Push_back(mi_->GetEntityList( EntityList::NODE_LIST,
                                           mi_->regionData[i].name ) );
      }
    }
  }

  void SimInputCGNS::GetResult( UInt sequenceStep,
                                UInt stepNum,
                                shared_ptr<BaseResult> result,
                                bool isHistory )
  {
    if(sequenceStep != 1) 
    {
      EXCEPTION("Only one static (RANS) sequence step supported!\n" <<
                "Note, that this is no limitation of CGNS but " <<
                "of the current implementation.");
    }

    if(stepNum != 1) 
    {
      EXCEPTION("Only one time step supported!\n" <<
                "Note, that this is no limitation of CGNS but " <<
                "of the current implementation.");
    }

    // determine region for this results
    std::string regionName =  result->GetEntityList()->GetName();

    switch( result->GetResultInfo()->definedOn ) {
    case ResultInfo::NODE:
      break;
    default:
      EXCEPTION( "Currently only results on nodes "
                 << "can be read in from a CGNS file ");
    }

    //UNSTRUCTURED VERSION!!!!!
    Integer numsols = 0; 
    char solname[33];
    GridLocation_t solLocation;
    
    //open the cgns file
    Integer fn = GetFileHandle(fileName_);
    
    cg_nsols(fn,1,1,&numsols);
    cg_sol_info(fn, 1, 1, numsols, solname , &solLocation );
    if(solLocation != Vertex){
      EXCEPTION("Solution is defined on the cell center.\n" <<
                "Mapping to nodes not supported right now.");
    }

    UInt numDofs = result->GetResultInfo()->dofNames.GetSize();
    UInt numEntities = result->GetEntityList()->GetSize();
    UInt resVecSize =  numEntities * numDofs;

    Vector<Double> & resVec = dynamic_cast<Result<Double>& >(*result).GetVector();
    resVec.Resize( resVecSize );

    Integer nfields = 0;
    char fieldName[33];
    DataType_t datatype;

    cg_nfields(fn, 1, 1, numsols, &nfields ); 
    for(UInt aSol = 1 ; aSol <= (UInt) nfields; aSol++){
      cg_field_info(fn, 1, 1, numsols, aSol, &datatype , fieldName );

      UInt spaceIdx = MapVelocityIndex(fieldName);
      if(spaceIdx != 9999 &&
         result->GetResultInfo()->resultType == FLUIDMECH_VELOCITY ) {
        cgsize_t range_min[3] = {1,1,1};
        cgsize_t range_max[3] = {cgsize_t(numVertices_),1,1};
        Vector<Double> curSol(numVertices_);
        cg_field_read(fn,1,1,1, fieldName, RealDouble , range_min, range_max, (void *)&curSol[0] );
        
        EntityIterator it = result->GetEntityList()->GetIterator();        
        for( UInt i = 0; i < numEntities; i++, it++ ) {
          UInt idx = it.GetNode()-1-nodeOffset_;
          resVec[i*numDofs+spaceIdx] = curSol[idx];
        }        
      }
      
      if(strcmp(fieldName,"Pressure") == 0 &&
         result->GetResultInfo()->resultType == FLUIDMECH_PRESSURE){
        cgsize_t range_min[3] = {1,1,1};
        cgsize_t range_max[3] = {cgsize_t(numVertices_),1,1};
        Vector<Double> curSol(numVertices_);
        cg_field_read(fn,1,1,1, fieldName, RealDouble , range_min, range_max, (void *)&curSol[0] );

        EntityIterator it = result->GetEntityList()->GetIterator();        
        for( UInt i = 0; i < numEntities; i++, it++ ) {
          UInt idx = it.GetNode()-1-nodeOffset_;
          resVec[i] = curSol[idx];
        }
      }    
    }
    
    cg_close(fn);

  }

}
