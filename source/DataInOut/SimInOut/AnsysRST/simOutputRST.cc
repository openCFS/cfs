// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <complex>

#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/exception.hpp>
namespace fs = boost::filesystem;

#include "DataInOut/WriteInfo.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "General/environment.hh"
#include "Utils/StdVector.hh"
#include "Domain/elem.hh"
#include "DataInOut/Logging/cfslog.hh"
#include "DataInOut/CommandLine/BaseCommandLineHandler.hh"

#include "simOutputRST.hh"
#include "ResWrProtos.hh"


namespace CoupledField {
  
  // declare logging stream
  DECLARE_LOG(simOutputRST)
  DEFINE_LOG(simOutputRST, "SimOutputRST")
 
    SimOutputRST::SimOutputRST( const std::string& fileName,
                                ParamNode * outputNode )
      : SimOutput( fileName, outputNode ) {

    ENTER_FCN( "SimOutputRST::SimOutputRST" );

    // Initialize variables
    formatName_ = "rst";
    fileName_ = fileName;
    dirName_ = "simoutput_rst";
    
    capabilities_.insert( MESH );
    capabilities_.insert( MESH_RESULTS );

    try 
    {
      fs::create_directory( dirName_ );
    } catch (std::exception &ex)
    {
      EXCEPTION(ex.what());
    }

    elemData_[0] = 1;  // material number
    elemData_[1] = 1;  // type number
    elemData_[2] = 1;  // real set number
    elemData_[3] = 1;  // section number
    elemData_[4] = 0;  // element coordinate system
    elemData_[5] = 0;  // birth-death flag  0, alive  1, dead
    elemData_[6] = 0;  // solid model reference
    elemData_[7] = 0;  // coded shape key
    elemData_[8] = 0;  // p-method exclude key
    elemData_[9] = 0;  // not used

    // linear element types
    ansysType2TypeId_[42] = 1;
    ansysType2TypeId_[45] = 2;

    // quadratic element types
    ansysType2TypeId_[82] = 1;
    ansysType2TypeId_[95] = 2;

    ansysType2NumNodes_[42] = 4;
    ansysType2NumNodes_[45] = 8;
    ansysType2NumNodes_[82] = 8;
    ansysType2NumNodes_[95] = 20;

  }


  SimOutputRST::~SimOutputRST() {
    ENTER_FCN( "SimOutputRST::~SimOutputRST" );
  }

  // ********
  //   Init
  // ********
  void SimOutputRST::Init( Grid* ptGrid) {
    
    ENTER_FCN( "SimOutputRST::OpenFile" );

    ptGrid_ = ptGrid;

    std::string filename;
    std::ostringstream strBuffer;

    try 
    {

      // Generate basename for output file
      filename.append( dirName_ );
      std::string pathsep = fs::path("/").native_directory_string();
      filename.append( pathsep );
      filename.append( fileName_ );
      filename.append( ".rst" );

      if(fs::exists(filename))
      {
        fs::remove( filename );
      }
    
    } catch (std::exception &ex)
    {
      EXCEPTION(ex.what());
    }

    Integer ret;
    Integer Nunit = 12;
    Integer Lunit = 6;
    char Fname[1024];
    Integer ncFname;
    char Title[80*2];
    char JobName[8];
    Integer Units = 1; // use SI units
    Integer NumDOF = 2;
    Integer DOF[] = {1, 2};
    Integer UserCode = 10;
    Integer MaxNode = ptGrid_->GetNumNodes();
    Integer NumNode = ptGrid_->GetNumNodes();
    Integer MaxElem = ptGrid_->GetNumElems();
    Integer NumElem = ptGrid_->GetNumElems();
    Integer MaxResultSet = 1;
    Integer lenTitle;
    Integer lenJobName;
    std::string tit = "Ein prima Testbeispiel!";

    std::fill(Fname, Fname+sizeof(Fname), 0);
    sprintf(Fname, "%s", filename.c_str());
    ncFname = filename.length();
    std::fill(Title, Title+sizeof(Title), ' ');
    std::copy(tit.begin(), tit.end(), Title);
    lenTitle = sizeof(Title);
    std::fill(JobName, JobName+sizeof(JobName), ' ');
    std::copy(filename.begin(), filename.begin()+sizeof(JobName), JobName);
    lenJobName = sizeof(JobName);
    
    ret = reswrbegin_(&Nunit,
                      &Lunit,
                      Fname,
                      &ncFname,
                      Title,
                      JobName,
                      &Units,
                      &NumDOF,
                      DOF,
                      &UserCode,
                      &MaxNode,
                      &NumNode,
                      &MaxElem,
                      &NumElem,
                      &MaxResultSet,
                      ncFname,
                      lenTitle,
                      lenJobName);

    if(ret != 0)
    {
      EXCEPTION("Could not open ANSYS RST file.");
    }

  }

  void SimOutputRST::WriteGrid() {
    ENTER_FCN( "SimOutputRST::WriteGrid" );
    
    LOG_TRACE(simOutputRST) << "Writing mesh";

    AnsysElementType ansysElementType[2];
    Integer num = 1;
    Integer MAXTYPE = 2;
    Integer NumType = 2;
    Integer MAXREAL = 0;
    Integer NumReal = 0;
    Integer MAXCSYS = 0;
    Integer NumCsys = 0;
    Double cSys[24];
    Double rCon;
    

    if(!ptGrid_->IsQuadratic())
    {
      //      ansysElementType[0] = 153; // SURF153
      SetElementType42(ansysElementType[0], 1); // PLANE42
      SetElementType45(ansysElementType[1], 2); // PLANE45
    }
    else
    {
      //      ansysElementType[0] = 153; // SURF153
      SetElementType82(ansysElementType[1], 1); // SOLID82
      SetElementType95(ansysElementType[1], 2); // SOLID95
    }

    reswrgeombegin_(&MAXTYPE, &NumType,
                    &MAXREAL, &NumReal,
                    &MAXCSYS, &NumCsys);

    // *****  element types  *****
    reswrtypebegin_();
    for(UInt i=0; i<NumType; i++)
    {
      reswrtype_(&ansysElementType[i].elementtypid,
                 ansysElementType[i].ielc);
    }
    reswrtypeend_();

    // *****  real constants  *****
    //    reswrrealbegin_();
    //    num = 1; reswrreal_(&num, &num, &rCon);
    //    reswrrealend_();

    // *****  Coordinate systems  *****
    //    reswrcsysbegin_();
    //    num = 1; reswrcsys_(&num, cSys);
    //    reswrcsysend_();

    // write nodes
    WriteNodes();
    
    // write elements
    WriteElements();

    reswrgeomend_();

    if ( commandLine->GetPrintGrid() == true ) {
      // Close result file
      reswrend_();
    }
  }

  void SimOutputRST::WriteNodes() {
    ENTER_FCN( "SimOutputRST::WriteNodes" );

    LOG_TRACE(simOutputRST) << "Writing nodes";
    
    // get number of nodes
    UInt numNodes = ptGrid_->GetNumNodes();
    Integer i;
    
    reswrnodebegin_();

    ansysNode_[3] = 0;
    ansysNode_[4] = 0;
    ansysNode_[5] = 0;

    Point point;
    for ( i = 1; i <= numNodes; i++ ) {
      ptGrid_->GetNodeCoordinate(point,i);
      ansysNode_[0] = point[0];
      ansysNode_[1] = point[1];
      ansysNode_[2] = point[2];

      reswrnode_(&i, ansysNode_);
    } 
    
    reswrnodeend_();

    LOG_TRACE(simOutputRST) << "Finished writing nodes" << std::endl;
  }

  void SimOutputRST::WriteElements() {
    ENTER_FCN( "SimOutputRST::WriteElements" );
   
    LOG_TRACE(simOutputRST) << "Writing elements";
    
    // iterate over all regions
    StdVector<RegionIdType> regionIds;
    StdVector<Elem*> elemVec;
    Elem * ptEl;
    Integer ansysElemMaterial;
    Integer ansysElemType;
    Integer elemNum;
    Integer numElemNodes;
    UInt iElem;
    FEType eType;

    ptGrid_->GetRegionIds( regionIds );
    
    reswrelembegin_();
    elemNum = 1;
    
    for ( UInt iReg = 0; iReg < regionIds.GetSize(); iReg++ ) {

      // get region id which serves as ANSYS element material
      ansysElemMaterial = regionIds[iReg]+1;

      ptGrid_->GetElems(elemVec,regionIds[iReg]);

      // write all element declarations
      for ( iElem = 0; iElem < elemVec.GetSize(); iElem++) {

        ptEl = elemVec[iElem];
        eType = ptEl->ptElem->feType();


        ansysElemType = ReorderConnect4Ansys(eType,
                                             ptEl->connect,
                                             connectANSYS_);
        numElemNodes = ansysType2NumNodes_[ansysElemType];

        elemData_[0] = ansysElemMaterial;
        elemData_[1] = ansysType2TypeId_[ansysElemType];
        
        reswrelem_(&elemNum, &numElemNodes,
                   connectANSYS_,
                   elemData_);

        elemNum++;
      }

    }

    reswrelemend_();
  }


  void SimOutputRST::RegisterResult( shared_ptr<BaseResult> sol,
                                     UInt saveBegin, UInt saveInc,
                                     UInt saveEnd ) {
    ENTER_FCN( "SimOutputRST::RegisterResult" );
  }

  void SimOutputRST::BeginStep( UInt stepNum, Double stepVal ) {
    ENTER_FCN( "SimOutputRST::BeginStep" );
  }
 
  void SimOutputRST::AddResult( shared_ptr<BaseResult> sol ) {
    ENTER_FCN( " SimOutputRST::AddResult" );
  }
  
  void SimOutputRST::FinishStep( ) {
    ENTER_FCN( "SimOutputRST::FinishStep" );
  }
  

  void SimOutputRST::Finalize() 
  {
    // Close result file
    reswrend_();
  }


  void SimOutputRST::
  WriteNodeElemDataTrans( const Vector<Double> & var, 
                          const StdVector<std::string> & dofNames,
                          const std::string& name, 
                          ResultInfo::EntryType entryType,
                          ResultInfo::EntityUnknownType entityType,
                          Double time ) {
    ENTER_FCN ( "SimOutputRST::WriteNodeElemDataTrans" );
  }     
  
  void SimOutputRST::
  WriteNodeElemDataHarm( const Vector<Complex> & var, 
                         const StdVector<std::string> & dofNames,
                         const std::string name, 
                         const ResultInfo::EntryType entryType,
                         ResultInfo::EntityUnknownType entityType,
                         const Double freq, 
                         const ComplexFormat outputFormat ) {
    ENTER_FCN ( "SimOutputRST::WriteNodeElemDataHarm" );
  }
  

  // change node order for right access in ANSYS
  Integer SimOutputRST::ReorderConnect4Ansys(FEType eType,
                                             StdVector<UInt>& connect,
                                             Integer* connectANSYS)
  {
    UInt numElemNodes = NUM_ELEM_NODES[eType];
    Integer ansysElemType;
    const Integer PLANE42 = 42; // linear quad
    const Integer PLANE82 = 82; // quadratic quad
    const Integer SOLID45 = 45; // linear hexa
    const Integer SOLID95 = 95; // quadratic hexa
   
    switch(eType)
    {
    case ET_TRIA3:
      std::copy(&connect[0], &connect[0]+numElemNodes, &connectANSYS[0]);
      connectANSYS[3] = connect[2];
      ansysElemType = PLANE42;
      break;

    case ET_TRIA6:
      connectANSYS[0] = connect[0];
      connectANSYS[1] = connect[1];
      connectANSYS[2] = connect[2];
      connectANSYS[3] = connect[2];
      connectANSYS[4] = connect[3];
      connectANSYS[5] = connect[4];
      connectANSYS[6] = connect[2];
      connectANSYS[7] = connect[5];
      ansysElemType = PLANE82;
      break;

    case ET_QUAD4:
      std::copy(&connect[0], &connect[0]+numElemNodes, &connectANSYS[0]);
      ansysElemType = PLANE42;
      break;

    case ET_QUAD8:
      std::copy(&connect[0], &connect[0]+numElemNodes, &connectANSYS[0]);
      ansysElemType = PLANE82;
      break;

    case ET_QUAD9:
      std::copy(&connect[0], &connect[0]+numElemNodes, &connectANSYS[0]);
      ansysElemType = PLANE82;
      break;

    case ET_TET4:
      connectANSYS[0] = connect[0]; // I
      connectANSYS[1] = connect[1]; // J
      connectANSYS[2] = connect[2]; // K
      connectANSYS[3] = connect[2]; // L
      connectANSYS[4] = connect[3]; // M
      connectANSYS[5] = connect[3]; // N
      connectANSYS[6] = connect[3]; // O
      connectANSYS[7] = connect[3]; // P
      ansysElemType = SOLID45;      
      break;                        
                                    
    case ET_TET10:                  
      connectANSYS[0] = connect[0];  // I
      connectANSYS[1] = connect[1];  // J
      connectANSYS[2] = connect[2];  // K
      connectANSYS[3] = connect[2];  // L
      connectANSYS[4] = connect[3];  // M
      connectANSYS[5] = connect[3];  // N
      connectANSYS[6] = connect[3];  // O
      connectANSYS[7] = connect[3];  // P
      connectANSYS[8] = connect[4];  // Q
      connectANSYS[9] = connect[5];  // R
      connectANSYS[10] = connect[2]; // S
      connectANSYS[11] = connect[6]; // T
      connectANSYS[12] = connect[3]; // U
      connectANSYS[13] = connect[3]; // V
      connectANSYS[14] = connect[3]; // W
      connectANSYS[15] = connect[3]; // X
      connectANSYS[16] = connect[7]; // Y
      connectANSYS[17] = connect[8]; // Z
      connectANSYS[18] = connect[9]; // A
      connectANSYS[19] = connect[9]; // B
      ansysElemType = SOLID95;
      break;

    case ET_WEDGE6:
      connectANSYS[0] = connect[0];  // I
      connectANSYS[1] = connect[1];  // J
      connectANSYS[2] = connect[2];  // K
      connectANSYS[3] = connect[2];  // L
      connectANSYS[4] = connect[3];  // M
      connectANSYS[5] = connect[4];  // N
      connectANSYS[6] = connect[5];  // O
      connectANSYS[7] = connect[5];  // P
      ansysElemType = SOLID45;
      break;

    case ET_WEDGE15:
      connectANSYS[0] = connect[0];    // I
      connectANSYS[1] = connect[1];    // J
      connectANSYS[2] = connect[2];    // K
      connectANSYS[3] = connect[2];    // L
      connectANSYS[4] = connect[3];    // M
      connectANSYS[5] = connect[4];    // N
      connectANSYS[6] = connect[5];    // O
      connectANSYS[7] = connect[5];    // P
      connectANSYS[8] = connect[6];    // Q
      connectANSYS[9] = connect[7];    // R
      connectANSYS[10] = connect[2];   // S
      connectANSYS[11] = connect[8];   // T
      connectANSYS[12] = connect[9];   // U
      connectANSYS[13] = connect[10];  // V
      connectANSYS[14] = connect[5];   // W
      connectANSYS[15] = connect[11];  // X
      connectANSYS[16] = connect[12];  // Y
      connectANSYS[17] = connect[13];  // Z
      connectANSYS[18] = connect[14];  // A
      connectANSYS[19] = connect[14];  // B
      ansysElemType = SOLID95;
      break;

    case ET_PYRA5:
      connectANSYS[0] = connect[0];  // I
      connectANSYS[1] = connect[1];  // J
      connectANSYS[2] = connect[2];  // K
      connectANSYS[3] = connect[3];  // L
      connectANSYS[4] = connect[4];  // M
      connectANSYS[5] = connect[4];  // N
      connectANSYS[6] = connect[4];  // O
      connectANSYS[7] = connect[4];  // P
      ansysElemType = SOLID45;           
      break;

    case ET_PYRA13:
      connectANSYS[0] = connect[0];    // I 
      connectANSYS[1] = connect[1];    // J 
      connectANSYS[2] = connect[2];    // K 
      connectANSYS[3] = connect[3];    // L 
      connectANSYS[4] = connect[4];    // M 
      connectANSYS[5] = connect[4];    // N 
      connectANSYS[6] = connect[4];    // O 
      connectANSYS[7] = connect[4];    // P 
      connectANSYS[8] = connect[5];    // Q 
      connectANSYS[9] = connect[6];    // R 
      connectANSYS[10] = connect[7];   // S 
      connectANSYS[11] = connect[8];   // T 
      connectANSYS[12] = connect[4];   // U 
      connectANSYS[13] = connect[4];   // V 
      connectANSYS[14] = connect[4];   // W 
      connectANSYS[15] = connect[4];   // X 
      connectANSYS[16] = connect[9];   // Y 
      connectANSYS[17] = connect[10];  // Z 
      connectANSYS[18] = connect[11];  // A 
      connectANSYS[19] = connect[12];  // B 
      ansysElemType = SOLID95;
      break;

    case ET_HEXA8:
      std::copy(&connect[0], &connect[0]+numElemNodes, &connectANSYS[0]);
      ansysElemType = SOLID45;
      break;

    case ET_HEXA20:
      std::copy(&connect[0], &connect[0]+numElemNodes, &connectANSYS[0]);
      ansysElemType = SOLID95;
      break;

    case ET_HEXA27:
      std::copy(&connect[0], &connect[0]+numElemNodes, &connectANSYS[0]);
      ansysElemType = SOLID95;
      break;

    default:
      ansysElemType = -1;
      break;
    }

    return ansysElemType;
  }


  void SimOutputRST::SetElementType42(AnsysElementType& eltype,
                                      Integer TypeNumber)
  {
    int i;

    eltype.elementtypid = TypeNumber;
    std::fill(eltype.ielc, eltype.ielc+IELCSZ, 0);

    eltype.ielc[IETYP-1] = TypeNumber;

    eltype.ielc[JETYP-1] = 42;
    eltype.ielc[21] = 4;
    eltype.ielc[25] = 1;
    eltype.ielc[27] = 1;
    eltype.ielc[33] = 3;
    eltype.ielc[34] = 1;
    eltype.ielc[46] = 1;
    eltype.ielc[49] = 1;
    eltype.ielc[51] = 4;
    eltype.ielc[55] = 2 ;
    eltype.ielc[60] = 4;
    eltype.ielc[61] = 4;
    eltype.ielc[60] = 	4;
    eltype.ielc[61] = 	4;
    eltype.ielc[62] = 	4;
    eltype.ielc[63] = 	2;
    eltype.ielc[64] = 	4;
    eltype.ielc[65] = 	22;
    eltype.ielc[66] = 	4;
    eltype.ielc[73] = 	2;
    eltype.ielc[74] = 	4;
    eltype.ielc[82] = 	4;
    eltype.ielc[83] = 	4;
    eltype.ielc[93] = 	4;
    eltype.ielc[94] = 	1;
    eltype.ielc[105] = 	1;
    eltype.ielc[108] = 	8;
    eltype.ielc[109] = 	25;
    eltype.ielc[110] = 	12;
    eltype.ielc[111] = 	139;
    eltype.ielc[112] = 	147;
    eltype.ielc[116] = 	1;
    eltype.ielc[121] = 	122;
    eltype.ielc[125] = 	152;
    eltype.ielc[126] = 	24;
    eltype.ielc[127] = 	108;
    eltype.ielc[133] = 	1347174734;
    eltype.ielc[134] = 	1161048608;
    eltype.ielc[136] = 	1;
    eltype.ielc[138] = 	1;
    eltype.ielc[140] = 	55;
    eltype.ielc[144] = 	1;


  }

  void SimOutputRST::SetElementType45(AnsysElementType& eltype,
                                      Integer TypeNumber)
  {
    int i;
    eltype.elementtypid = TypeNumber;
    std::fill(eltype.ielc, eltype.ielc+IELCSZ, 0);

    eltype.ielc[IETYP-1] = TypeNumber;

    eltype.ielc[JETYP-1] = 45;
    eltype.ielc[20] = 3;
    eltype.ielc[21] = 6;
    eltype.ielc[22] = 3;
    eltype.ielc[25] = 1;
    eltype.ielc[27] = 1;
    eltype.ielc[33] = 7;
    eltype.ielc[34] = 1;
    eltype.ielc[40] = 1;
    eltype.ielc[46] = 1;
    eltype.ielc[51] = 6;
    eltype.ielc[55] = 2 ;
    eltype.ielc[60] = 8;
    eltype.ielc[61] = 8;
    eltype.ielc[62] = 8;
    eltype.ielc[63] = 3;
    eltype.ielc[64] = 8;

	
    eltype.ielc[65] = 	22;
    eltype.ielc[66] = 	8;
    eltype.ielc[73] = 	4;
    eltype.ielc[74] = 	6;
    eltype.ielc[82] = 	8;
    eltype.ielc[83] = 	8;
    eltype.ielc[93] = 	8;
    eltype.ielc[94] = 	1;

    eltype.ielc[105] = 1;
    eltype.ielc[108] = 24;
    eltype.ielc[109] = 48;
    eltype.ielc[110] = 24;
    eltype.ielc[111] = 600;
    eltype.ielc[112] = 356;
    eltype.ielc[116] = 1;
    eltype.ielc[121] = 415;

    eltype.ielc[125] = 234;
    eltype.ielc[126] = 64;
    eltype.ielc[127] = 312;
    eltype.ielc[131] = 1;
    eltype.ielc[133] = 1397705801;
    eltype.ielc[134] = 1144272160;
    eltype.ielc[136] = 1;
    eltype.ielc[139] = 1;
    eltype.ielc[140] = 70;
    eltype.ielc[144] = 1;


  }

  void SimOutputRST::SetElementType82(AnsysElementType& eltype,
                                      Integer TypeNumber)
  {
    eltype.elementtypid = TypeNumber;

    std::fill(eltype.ielc, eltype.ielc+IELCSZ, 0);

    eltype.ielc[IETYP-1] = TypeNumber;

    eltype.ielc[JETYP-1] = 82;
    eltype.ielc[ISHAP -1] = 4;
    eltype.ielc[MNODE -1] = 2;
    eltype.ielc[KELSTO-1] = 1;
    eltype.ielc[KDOFS -1] = 3;
    eltype.ielc[NMVCTF-1]= 1;
    eltype.ielc[MATRQD-1] = 1;
    eltype.ielc[NCOMPN-1] = 4;
    eltype.ielc[MTSYM -1] = 2;
    eltype.ielc[NMNDMX-1] = 8;
    eltype.ielc[NMNDMN-1] = 4;
    eltype.ielc[NMNDST-1] = 8;
    eltype.ielc[NMDFPN-1]=2;
    eltype.ielc[NMNDAC-1]=8;
    eltype.ielc[MATRXS-1] = 22;
    eltype.ielc[NMNDNE-1] = 8;
    eltype.ielc[NMPTSF-1] = 2;
    eltype.ielc[NMPRES-1] = 4;
    eltype.ielc[NMTEMP-1] = 8;
    eltype.ielc[NMFLNC-1] = 8;
    eltype.ielc[NMNDNO-1] = 4;
    eltype.ielc[KCONIT-1] = 1;
    eltype.ielc[KNORM -1] = 1;
    eltype.ielc[NMSMIS-1] = 8;
    eltype.ielc[NMNMIS-1] = 29;
    eltype.ielc[NMNMUP-1] = 24;
    eltype.ielc[NCPTM1-1] = 156;
    eltype.ielc[NCPTM2-1] = 173;
    eltype.ielc[JSTPR -1] = 1;
    eltype.ielc[NMSSVR-1] =  54;
    eltype.ielc[NMNSVR-1] = 136;
    eltype.ielc[NMPSVR-1] = 24;
    eltype.ielc[NMCSVR-1] = 164;
    eltype.ielc[NAMRF1-1] = 1347174734;
    eltype.ielc[NAMRF2-1] = 1161310752;
    eltype.ielc[ADADES-1] = 1;
    eltype.ielc[NLSSAD-1] = 1;
    eltype.ielc[KSWTTS-1] = 77;

  }

  void SimOutputRST::SetElementType95(AnsysElementType& eltype,
                                      Integer TypeNumber)
  {
    eltype.elementtypid = TypeNumber;

    std::fill(eltype.ielc, eltype.ielc+IELCSZ, 0);

    eltype.ielc[IETYP-1] = TypeNumber;

    eltype.ielc[JETYP-1] = 95;

    eltype.ielc[20] = 3;
    eltype.ielc[21] = 6;
    eltype.ielc[22] = 3;
    eltype.ielc[24] = 2;
    eltype.ielc[27] = 1;
    eltype.ielc[33] = 7;
    eltype.ielc[34] = 1;
    eltype.ielc[40] = 0;
    eltype.ielc[46] = 1;
    eltype.ielc[51] = 6;
    eltype.ielc[55] = 2 ;
    eltype.ielc[60] = 20;
    eltype.ielc[61] = 8;
    eltype.ielc[62] = 20;
    eltype.ielc[63] = 3;
    eltype.ielc[64] = 20;

	
    eltype.ielc[65] = 	22;
    eltype.ielc[66] = 	20;
    eltype.ielc[73] = 	4;
    eltype.ielc[74] = 	6;
    eltype.ielc[82] = 	20;
    eltype.ielc[83] = 	20;
    eltype.ielc[93] = 	8;
    eltype.ielc[94] = 	1;

    eltype.ielc[105] = 1;
    eltype.ielc[108] = 24;
    eltype.ielc[109] = 60;
    eltype.ielc[110] = 60;
    eltype.ielc[111] = 2364;
    eltype.ielc[112] = 972;
    eltype.ielc[116] = 1;
    eltype.ielc[121] = 0;

    eltype.ielc[125] = 402;
    eltype.ielc[126] = 112;
    eltype.ielc[127] = 602;
    eltype.ielc[131] = 1;
    eltype.ielc[133] = 1397705801;
    eltype.ielc[134] = 1144599840;
    eltype.ielc[136] = 1;
    eltype.ielc[138] = 1;
    eltype.ielc[140] = 90;
    eltype.ielc[144] = 1;


  }

} // end of namespace
