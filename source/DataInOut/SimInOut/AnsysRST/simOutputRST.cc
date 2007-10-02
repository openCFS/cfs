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


    // Initialize variables
    formatName_ = "rst";
    fileName_ = fileName;
    dirName_ = "simoutput_rst";

    numSteps_ = 0;
    msStep_ = 0;
    stepNumOffset_ = 0;
    stepValOffset_ = 0;
    
    capabilities_.insert( MESH );
    capabilities_.insert( MESH_RESULTS );

    try 
    {
      fs::create_directory( dirName_ );
    } catch (std::exception &ex)
    {
      EXCEPTION(ex.what());
    }

    printGrid_ = commandLine->GetPrintGrid();
    printNow_ = false;

    // Initialize element data array.
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

    // Map RST elem type id to Ansys element types.
    // linear element types
    ansysType2TypeId_[PLANE42] = 1;
    ansysType2TypeId_[SOLID45] = 2;

    // quadratic element types
    ansysType2TypeId_[PLANE82] = 1;
    ansysType2TypeId_[SOLID95] = 2;

    // Store the number of nodes of the used elem types in a map.
    ansysType2NumNodes_[PLANE42] = 4;
    ansysType2NumNodes_[SOLID45] = 8;
    ansysType2NumNodes_[PLANE82] = 8;
    ansysType2NumNodes_[SOLID95] = 20;

    // Map internal element types to ANSYS element types.
    elemType2AnsysType_[ET_TRIA3]   = PLANE42;
    elemType2AnsysType_[ET_TRIA6]   = PLANE82;
    elemType2AnsysType_[ET_QUAD4]   = PLANE42;
    elemType2AnsysType_[ET_QUAD8]   = PLANE82;
    elemType2AnsysType_[ET_QUAD9]   = PLANE82;
    elemType2AnsysType_[ET_TET4]    = SOLID45;
    elemType2AnsysType_[ET_TET10]   = SOLID95;                  
    elemType2AnsysType_[ET_WEDGE6]  = SOLID45;
    elemType2AnsysType_[ET_WEDGE15] = SOLID95;
    elemType2AnsysType_[ET_PYRA5]   = SOLID45;
    elemType2AnsysType_[ET_PYRA13]  = SOLID95;
    elemType2AnsysType_[ET_HEXA8]   = SOLID45;
    elemType2AnsysType_[ET_HEXA20]  = SOLID95;
    elemType2AnsysType_[ET_HEXA27]  = SOLID95;

    // Ansys nodal dof labels.
    dofLabels_.push_back("UX  ");
    dofLabels_.push_back("UY  ");
    dofLabels_.push_back("UZ  ");
    dofLabels_.push_back("ROTX");
    dofLabels_.push_back("ROTY");
    dofLabels_.push_back("ROTZ");
    dofLabels_.push_back("AX  ");
    dofLabels_.push_back("AY  ");
    dofLabels_.push_back("AZ  ");
    dofLabels_.push_back("VX  ");
    dofLabels_.push_back("VY  ");
    dofLabels_.push_back("VZ  ");
    dofLabels_.push_back("RS13");
    dofLabels_.push_back("RS14");
    dofLabels_.push_back("RS15");
    dofLabels_.push_back("RS16");
    dofLabels_.push_back("RS17");
    dofLabels_.push_back("RS18");
    dofLabels_.push_back("PRES");
    dofLabels_.push_back("TEMP");
    dofLabels_.push_back("VOLT");
    dofLabels_.push_back("MAG ");
    dofLabels_.push_back("ENKE");
    dofLabels_.push_back("ENDS");
    dofLabels_.push_back("EMF ");
    dofLabels_.push_back("CURR");
    dofLabels_.push_back("SP01");
    dofLabels_.push_back("SP02");
    dofLabels_.push_back("SP03");
    dofLabels_.push_back("SP04");
    dofLabels_.push_back("SP05");
    dofLabels_.push_back("SP06");

    // Number of DOFs for Ansys element results
    numAnsysElemDOFs_[ENS] = 11;
    numAnsysElemDOFs_[EEL] = 7;
    numAnsysElemDOFs_[EPL] = 7;
    numAnsysElemDOFs_[ECR] = 7;
    numAnsysElemDOFs_[ETH] = 7;    
  }


  SimOutputRST::~SimOutputRST() {
  }

  void SimOutputRST::Init( Grid* ptGrid, bool printGridOnly ) {
    

    ptGrid_ = ptGrid;
    WriteGrid();
  }

  void SimOutputRST::WriteGrid() {
    
    LOG_TRACE(simOutputRST) << "Writing mesh";

    // variables for global ANSYS file header
    Integer ret;
    Integer Nunit = 12; // has to be set for .rst files
    Integer Lunit = 6;
    char Fname[1024];
    Integer ncFname;
    char Title[80*2];
    char JobName[8];
    Integer Units = 1; // use SI units
    std::vector<Integer> DOF;
    Integer UserCode = 10;
    Integer MaxNode = ptGrid_->GetNumNodes();
    Integer NumNode = ptGrid_->GetNumNodes();
    Integer MaxElem = ptGrid_->GetNumElems();
    Integer NumElem = ptGrid_->GetNumElems();
    Integer MaxResultSet; // Should be num timesteps
    Integer kan;
    Integer lenTitle;
    Integer lenJobName;

    // variables for element information header
    AnsysElementType ansysElementType[2];
    Integer num = 1;
    Integer MAXTYPE = 2;
    Integer NumType = 2;
    Integer MAXREAL = 0;
    Integer NumReal = 0;
    Integer MAXCSYS = 0;
    Integer NumCsys = 0;

    // Only write grid if corresponding flags are set.    
    if(printGrid_ || printNow_)
    {
      // Get number of DOFs
      std::map <AnsysNodalDof, UInt>::const_iterator dofIdxIt, dofIdxEnd;
      dofIdxIt = ansysNodeDof2Idx_.begin();
      dofIdxEnd = ansysNodeDof2Idx_.end();
      
      // Build up array with dofs.
      if(printNow_)
        for( ; dofIdxIt != dofIdxEnd; dofIdxIt++ )
        {
          DOF.push_back(dofIdxIt->first);
        }

      // For some reason writing the grid needs one DOF.
      if(printGrid_)
      {
        numNodeDOFs_ = 1;
        DOF.push_back(SP01);
      }
      
      // Variables for construction of file name
      std::string filename;
      std::ostringstream strBuffer;
      std::stringstream fnstream;
      std::string pathsep;

      // Generate basename for output file
      // If just the grid is written the file will end with
      // "_grid.rst". If we are in multistep N the file will end
      // with "_msN.rst".
      try 
      {
        pathsep = fs::path("/").native_directory_string();
      } catch (std::exception &ex)
      {
        EXCEPTION(ex.what());
      }
      fnstream << dirName_ << pathsep << fileName_;
      if(msStep_ != 0)
        fnstream << "_ms" << msStep_;
      else
        fnstream << "_grid";
      fnstream << ".rst";

      filename = fnstream.str();

      try 
      {
        // If file already exists, delete it.
        // ANSYS API can not handle this case.
        if(fs::exists(filename))
        {
          fs::remove( filename );
        }
      } catch (std::exception &ex)
      {
        EXCEPTION(ex.what());
      }

      // Copy file name and job name to buffers.
      std::fill(Fname, Fname+sizeof(Fname), 0);
      sprintf(Fname, "%s", filename.c_str());
      ncFname = filename.length();
      std::fill(Title, Title+sizeof(Title), ' ');
      std::copy(fileName_.begin(), fileName_.end(), Title);
      lenTitle = sizeof(Title);
      std::fill(JobName, JobName+sizeof(JobName), ' ');
      std::copy(filename.begin(), filename.begin()+sizeof(JobName), JobName);
      lenJobName = sizeof(JobName);
    
      // Determine maximum number of result sets, if results are real/complex,
      // and set ANSYS analysis type according to internal analysis type. 
      switch(analysisType_)
      {
      case STATIC:
        MaxResultSet = ((numSteps_/100)+1)*100;
        complexResults_ = false;
        timeFreq_ = "time";
        kan = 0;
        break;
      case TRANSIENT:
        MaxResultSet = ((numSteps_/100)+1)*100;
        complexResults_ = false;
        timeFreq_ = "time";
        kan = 4;
        break;
      case HARMONIC:
        MaxResultSet = ((2*numSteps_/100)+1)*100;
        complexResults_ = true;
        timeFreq_ = "frequency";
        kan = 3;
        break;
      case EIGENFREQUENCY:
        MaxResultSet = ((numSteps_/100)+1)*100;
        complexResults_ = false;
        timeFreq_ = "frequency";
        kan = 2;
        break;

      default:
        // We only care about the analysis type if we don't just
        // want to output the grid. 
        if(!printGrid_)
        {
          EXCEPTION("ANSYS RST writer does not support current analyis type.");
        }
        else
        {
          kan = 0;
          MaxResultSet = 0;
        }
        
        break;
      }

      // Open ANSYS result file.
      ret = reswrbegin_(&Nunit,
                        &Lunit,
                        Fname,
                        &ncFname,
                        Title,
                        JobName,
                        &Units,
                        &numNodeDOFs_,
                        &DOF[0],
                        &UserCode,
                        &MaxNode,
                        &NumNode,
                        &MaxElem,
                        &NumElem,
                        &MaxResultSet,
                        &kan,
                        ncFname,
                        lenTitle,
                        lenJobName);

      if(ret != 0)
      {
        EXCEPTION("Could not open ANSYS RST file.");
      }

      // Initialize element headers.
      if(!ptGrid_->IsQuadratic())
      {
        //      ansysElementType[0] = 153; // SURF153
        SetElementType42(ansysElementType[0], 1); // PLANE42
        SetElementType45(ansysElementType[1], 2); // PLANE45
      }
      else
      {
        //      ansysElementType[0] = 153; // SURF153
        SetElementType82(ansysElementType[0], 1); // SOLID82
        SetElementType95(ansysElementType[1], 2); // SOLID95
      }

      // Write geometry to result file.
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

      // write nodes
      WriteNodes();
    
      // write elements
      WriteElements();

      // Finish writing geometry
      reswrgeomend_();
    
      if ( printGrid_ ) {
        // Close result file
        reswrend_();
      }
      else
      {
        // Initialize temporary element result arrays.
        tmpElemResultVecsReal_[ENS].resize(NumElem*numAnsysElemDOFs_[ENS]);
        tmpElemResultVecsReal_[EEL].resize(NumElem*numAnsysElemDOFs_[EEL]);
        tmpElemResultVecsReal_[EPL].resize(NumElem*numAnsysElemDOFs_[EPL]);
        tmpElemResultVecsReal_[ECR].resize(NumElem*numAnsysElemDOFs_[ECR]);
        tmpElemResultVecsReal_[ETH].resize(NumElem*numAnsysElemDOFs_[ETH]);

        tmpElemResultVecsImag_[ENS].resize(NumElem*numAnsysElemDOFs_[ENS]);
        tmpElemResultVecsImag_[EEL].resize(NumElem*numAnsysElemDOFs_[EEL]);
        tmpElemResultVecsImag_[EPL].resize(NumElem*numAnsysElemDOFs_[EPL]);
        tmpElemResultVecsImag_[ECR].resize(NumElem*numAnsysElemDOFs_[ECR]);
        tmpElemResultVecsImag_[ETH].resize(NumElem*numAnsysElemDOFs_[ETH]);
      }
    }
  }

  void SimOutputRST::WriteNodes() {

    LOG_TRACE(simOutputRST) << "Writing nodes";
    
    // get number of nodes
    UInt numNodes = ptGrid_->GetNumNodes();
    Integer i;
    
    reswrnodebegin_();

    // Set node rotations to zero.
    ansysNode_[3] = 0;
    ansysNode_[4] = 0;
    ansysNode_[5] = 0;

    Point point;
    for ( i = 1; i <= numNodes; i++ ) {
      // Set coordinates for node i
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

        // Get element from grid and determine element type.
        ptEl = elemVec[iElem];
        eType = ptEl->ptElem->feType();

        // Degenerate ANSYS mechanical elements and get ANSYS element type.
        ansysElemType = ReorderConnect4Ansys(eType,
                                             ptEl->connect,
                                             connectANSYS_);
        numElemNodes = ansysType2NumNodes_[ansysElemType];

        // Set element data (region -> material number & type id)
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

  void SimOutputRST::BeginMultiSequenceStep( UInt step,
                                             AnalysisType type,
                                             UInt numSteps)
  {

    printNow_ = false;

    numSteps_ = numSteps;
    analysisType_ = type;
    msStep_ = step;
  }

  void SimOutputRST::RegisterResult( shared_ptr<BaseResult> sol,
                                     UInt saveBegin, UInt saveInc,
                                     UInt saveEnd ) {
    
    ResultInfo & resInfo = *sol->GetResultInfo();

    MapInternal2ANSYSNodeDof(resInfo.resultType);
    MapInternal2ANSYSElemDof(resInfo.resultType);
  }

  void SimOutputRST::BeginStep( UInt stepNum, Double stepVal ) {

    stepNum_ = stepNum + stepNumOffset_;
    stepVal_ = stepVal + stepValOffset_;

    // At this stage we exactly know the number of nodal DOFs
    // and can therefore write the grid.
    numNodeDOFs_ = ansysNodeDof2Idx_.size();

    if(!printNow_)
    {
      printNow_ = true;
      WriteGrid();
    }

  }
 
  void SimOutputRST::AddResult( shared_ptr<BaseResult> sol ) {

    ResultInfo & actDof = *(sol->GetResultInfo());

    LOG_DBG(simOutputRST) << "Adding result '" 
                          << actDof.resultName  << "'";
      
    resultMap_[sol->GetResultInfo()->resultName].Push_back( sol );

  }
  
  void SimOutputRST::FinishStep( ) {

    LOG_TRACE(simOutputRST) << "Starting to finish Step";
    
    std::string title;
    UInt numDOFs;
    std::vector<Double> nodeValuesReal;
    std::vector<Double> nodeValuesImag;
    std::vector<Double> elemValuesReal;
    std::vector<Double> elemValuesImag;
    UInt numNodes = ptGrid_->GetNumNodes();
    UInt numElems = ptGrid_->GetNumElems();
    nodeValuesReal.resize(numNodes*numNodeDOFs_);
    nodeValuesImag.resize(numNodes*numNodeDOFs_);
    elemValuesReal.resize(numElems*256);
    elemValuesImag.resize(numElems*256);
    Integer node;
    Integer elem;
    
    // iterate over all result types
    ResultMapType::iterator it = resultMap_.begin();
    for( ; it != resultMap_.end(); it++ ) {
      
      // get result info object and results for current result type
      ResultInfo & actInfo = *(it->second[0]->GetResultInfo());
      const StdVector<shared_ptr<BaseResult> > actResults =
        it->second;

      title = actInfo.resultName;
      if( actInfo.unit.size() != 0 ) {
        title += " (";
        title += actInfo.unit;
        title += ")";
      }
      ResultInfo::EntryType entryType =  actInfo.entryType;
      ResultInfo::EntityUnknownType entityType = actInfo.definedOn;

      // check if result is defined on nodes or elements
      StdVector<std::string> & dofNames = actInfo.dofNames;  
      if(  actInfo.definedOn != ResultInfo::NODE &&
           actInfo.definedOn != ResultInfo::ELEMENT ) 
      {
        Warning( "RST can only write results on elements and nodes",
                 __FILE__, __LINE__ );
        continue;
      }

      LOG_DBG(simOutputRST) << "Writing result '" << title << "'";

      numDOFs = dofNames.GetSize();
      
      if( actResults[0]->GetEntryType() == EntryType::DOUBLE ) {
        Vector<Double> gSol;
        FillGlobalVec<Double>(gSol, actResults, entityType );
        
        if(actInfo.definedOn == ResultInfo::NODE)
          PrepareNodeData( gSol, nodeValuesReal,
                           actInfo.resultType,
                           numDOFs, entryType);

        if(actInfo.definedOn == ResultInfo::ELEMENT)
          PrepareElemData( gSol, actInfo.resultType,
                           numDOFs, entryType);
      } else {
        Vector<Complex> gSol;
        FillGlobalVec<Complex>(gSol, actResults, entityType );

        if(actInfo.definedOn == ResultInfo::NODE)
          PrepareNodeData( gSol, nodeValuesReal, nodeValuesImag,
                           actInfo.resultType,
                           numDOFs, entryType);
        
        if(actInfo.definedOn == ResultInfo::ELEMENT)
          PrepareElemData( gSol, actInfo.resultType,
                           numDOFs, entryType);        
      }
    }

    // Generate title for current result set.
    char steptitle[400];
    std::string dummy;
    std::stringstream titstr;
    std::fill(steptitle, steptitle+sizeof(steptitle), ' ');
    titstr << fileName_ << " step: " << stepNum_ << ", "
           << timeFreq_ << ": " << stepVal_;
    if(complexResults_)
      titstr << ", real";
    
    dummy = titstr.str();
    std::copy(dummy.begin(),
              dummy.end(),
              steptitle);

    std::string dofLab;
    dofLab.resize(numNodeDOFs_*4);

    std::map <AnsysNodalDof, UInt>::const_iterator dofIdxIt, dofIdxEnd;
    dofIdxIt = ansysNodeDof2Idx_.begin();
    dofIdxEnd = ansysNodeDof2Idx_.end();
    UInt idx = 0;
    
    // Build up array with dof labels.
    for( ; dofIdxIt != dofIdxEnd; dofIdxIt++ )
    {
      char* pt = &dofLabels_[dofIdxIt->first-1][0];
      
      idx = dofIdxIt->second*4;

      std::copy(pt, pt+4, &dofLab[idx]);
    }

    Integer kcmplx = 0;
    Integer subStep = stepNum_;
    reswrsolbegin_((Integer*) &stepNum_,
                   &subStep,
                   (Integer*) &stepNum_,
                   &stepVal_,
                   &kcmplx,
                   &steptitle[0], &dofLab[0], dofLab.length());


    // Write the vector of nodal values to the RST file.
    // The functions are called *disp* for displacement but
    // that does mean nothing since we previously specified
    // that we have all 32 ANSYS dofs.
    //    reswrdispbegin_();
    //    for (node = 1; node <= numNodes; node++ ) {
    //      reswrdisp_(&node, &nodeValues[(node-1)*numNodeDOFs_]);
    //    }
    //    reswrdispend_();
    reswrnodaldata_(&nodeValuesReal[0]);

    reswreresbegin_();

    std::map <SolutionType, AnsysElemDof>::const_iterator eResIt, eResEnd;
    eResEnd = internal2AnsysElemDofMap_.end();
    StdVector<Elem*> elemVec;
    ptGrid_->GetElems(elemVec, ALL_REGIONS);

    for ( elem = 1; elem <= numElems; elem++ ) {

      eResIt = internal2AnsysElemDofMap_.begin();
      FEType eType = elemVec[elem-1]->ptElem->feType();
      numNodes = ansysType2NumNodes_[elemType2AnsysType_[eType]];

      for ( ; eResIt != eResEnd; eResIt++ ) {
        AnsysElemDof ansysDOF = eResIt->second;
        std::vector<Double>& resVec = tmpElemResultVecsReal_[ansysDOF];
        UInt numAnsysElemDofs = numAnsysElemDOFs_[ansysDOF];
        UInt idx = (elem-1)*numAnsysElemDofs;

        for ( UInt i = 0; i <numNodes; i++ ) {
          std::copy(&resVec[idx],
                    &resVec[idx+numAnsysElemDofs],
                    &elemValuesReal[i*numAnsysElemDofs]);
        }

        numAnsysElemDofs *=numNodes;
        reswrestrbegin_(&elem);
        Integer tmp = ansysDOF;
        reswrestr_(&tmp, (Integer*)&numAnsysElemDofs, &elemValuesReal[0]);
        reswrestrend_();
      }
    }

    reswreresend_();
    reswrsolend_();

    // Write imaginary part of result if one is present.
    if(complexResults_)
    {
      std::fill(steptitle, steptitle+sizeof(steptitle), ' ');
      titstr.str("");
      titstr << fileName_ << " step: " << stepNum_ << ", "
             << timeFreq_ << ": " << stepVal_ << ", imag";
      
      dummy = titstr.str();
      std::copy(dummy.begin(),
                dummy.end(),
                steptitle);

      kcmplx = 1;
      reswrsolbegin_((Integer*) &stepNum_,
                     &subStep,
                     (Integer*) &stepNum_,
                     &stepVal_,
                     &kcmplx,
                     &steptitle[0], &dofLab[0], dofLab.length());

      reswrnodaldata_(&nodeValuesImag[0]);
      
      reswreresbegin_();
      eResIt, eResEnd;
      eResEnd = internal2AnsysElemDofMap_.end();
      elemVec;
      ptGrid_->GetElems(elemVec, ALL_REGIONS);

      for ( elem = 1; elem <= numElems; elem++ ) {

        eResIt = internal2AnsysElemDofMap_.begin();
        FEType eType = elemVec[elem-1]->ptElem->feType();
        numNodes = ansysType2NumNodes_[elemType2AnsysType_[eType]];

        for ( ; eResIt != eResEnd; eResIt++ ) {
          AnsysElemDof ansysDOF = eResIt->second;
          std::vector<Double>& resVec = tmpElemResultVecsImag_[ansysDOF];
          UInt numAnsysElemDofs = numAnsysElemDOFs_[ansysDOF];
          UInt idx = (elem-1)*numAnsysElemDofs;

          for ( UInt i = 0; i <numNodes; i++ ) {
            std::copy(&resVec[idx],
                      &resVec[idx+numAnsysElemDofs],
                      &elemValuesImag[i*numAnsysElemDofs]);
          }

          numAnsysElemDofs *=numNodes;
          reswrestrbegin_(&elem);
          Integer tmp = ansysDOF;
          reswrestr_(&tmp, (Integer*)&numAnsysElemDofs, &elemValuesReal[0]);
          reswrestrend_();
        }
      }
      reswreresend_();

      reswrsolend_();      
    }
  }
  
  void SimOutputRST::FinishMultiSequenceStep( )
  {
    // set offset for step value and number to last values
    stepNumOffset_ = stepNum_;
    stepValOffset_ = stepVal_;
    
    // Close result file
    reswrend_();

    internal2AnsysNodeDofMap_.clear();
    ansysNodeDof2Idx_.clear();
    internal2AnsysElemDofMap_.clear();
  }
  

  void SimOutputRST::Finalize() 
  {
  }


  void SimOutputRST::PrepareNodeData( const Vector<Double> & var,
                                    std::vector<Double> & nodeValuesReal,
                                    SolutionType solType,
                                    UInt numDOFs,
                                    const ResultInfo::EntryType entryType) {

    AnsysNodalDof ansysDOF = internal2AnsysNodeDofMap_[solType];
    UInt numNodes = ptGrid_->GetNumNodes();
    Integer idx;

    // Solution type has not been previously mapped.
    if(ansysNodeDof2Idx_.find(ansysDOF) == ansysNodeDof2Idx_.end())
      return;
    
    for ( UInt node = 0; node < numNodes; node++ ) {
      for(UInt dof = 0; dof < numDOFs; dof++)
      {
        idx = ansysNodeDof2Idx_[(AnsysNodalDof)(ansysDOF+dof)];

        nodeValuesReal[node*numNodeDOFs_+idx] = var[node*numDOFs+dof];
      }      
    }
  }     
  
  void SimOutputRST::PrepareElemData( const Vector<Double> & var, 
                                    SolutionType solType,
                                    UInt numDOFs,
                                    ResultInfo::EntryType entryType)
  {

    UInt numElems = ptGrid_->GetNumElems();
    Integer ansIdx, intIdx;
    AnsysElemDof ansysDof = internal2AnsysElemDofMap_[solType];
    Integer elem;
    UInt numAnsysElemDofs;
    
    // Solution type has not been previously mapped.
    if( ansysDof < 0 )
      return;

    std::vector<Double>& resVec = tmpElemResultVecsReal_[ansysDof];
    numAnsysElemDofs = numAnsysElemDOFs_[ansysDof];

    for ( elem = 1; elem <= numElems; elem++ ) {
      intIdx = (elem-1)*numDOFs;
      ansIdx = (elem-1)*numAnsysElemDofs;

      for(UInt i=0; i<numDOFs; i++)
        resVec[ansIdx+i] = var[intIdx+i];
    }

  }
  
  void SimOutputRST::PrepareNodeData( const Vector<Complex> & var,
                                    std::vector<Double> & nodeValuesReal,
                                    std::vector<Double> & nodeValuesImag,
                                    SolutionType solType,
                                    UInt numDOFs,
                                    const ResultInfo::EntryType entryType) {

    AnsysNodalDof ansysDOF = internal2AnsysNodeDofMap_[solType];
    UInt numNodes = ptGrid_->GetNumNodes();
    Integer idx;

    // Solution type has not been previously mapped.
    if(ansysNodeDof2Idx_.find(ansysDOF) == ansysNodeDof2Idx_.end())
      return;
    
    for ( UInt node = 0; node < numNodes; node++ ) {
      for(UInt dof = 0; dof < numDOFs; dof++)
      {
        idx = ansysNodeDof2Idx_[(AnsysNodalDof)(ansysDOF+dof)];

        nodeValuesReal[node*numNodeDOFs_+idx] = var[node*numDOFs+dof].real();
        nodeValuesImag[node*numNodeDOFs_+idx] = var[node*numDOFs+dof].imag();
      }      
    }
  }     
  
  void SimOutputRST::PrepareElemData( const Vector<Complex> & var, 
                                    SolutionType solType,
                                    UInt numDOFs,
                                    ResultInfo::EntryType entryType)
  {

    UInt numElems = ptGrid_->GetNumElems();
    Integer ansIdx, intIdx;
    AnsysElemDof ansysDof = internal2AnsysElemDofMap_[solType];
    Integer elem;
    UInt numAnsysElemDofs;
    
    // Solution type has not been previously mapped.
    if( ansysDof < 0 )
      return;

    std::vector<Double>& resVecReal = tmpElemResultVecsReal_[ansysDof];
    std::vector<Double>& resVecImag = tmpElemResultVecsImag_[ansysDof];
    numAnsysElemDofs = numAnsysElemDOFs_[ansysDof];

    for ( elem = 1; elem <= numElems; elem++ ) {
      intIdx = (elem-1)*numDOFs;
      ansIdx = (elem-1)*numAnsysElemDofs;

      for(UInt i=0; i<numDOFs; i++)
      {
        resVecReal[ansIdx+i] = var[intIdx+i].real();
        resVecImag[ansIdx+i] = var[intIdx+i].imag();
      }
    }

  }
  

  // change node order for right access in ANSYS
  Integer SimOutputRST::ReorderConnect4Ansys(FEType eType,
                                             StdVector<UInt>& connect,
                                             Integer* connectANSYS)
  {
    UInt numElemNodes = NUM_ELEM_NODES[eType];
    Integer ansysElemType = elemType2AnsysType_[eType];
   
    switch(eType)
    {
    case ET_TRIA3:
      std::copy(&connect[0], &connect[0]+numElemNodes, &connectANSYS[0]);
      connectANSYS[3] = connect[2];
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
      break;

    case ET_QUAD4:
      std::copy(&connect[0], &connect[0]+numElemNodes, &connectANSYS[0]);
      break;

    case ET_QUAD8:
      std::copy(&connect[0], &connect[0]+numElemNodes, &connectANSYS[0]);
      break;

    case ET_QUAD9:
      std::copy(&connect[0], &connect[0]+numElemNodes, &connectANSYS[0]);
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
      break;

    case ET_HEXA8:
      std::copy(&connect[0], &connect[0]+numElemNodes, &connectANSYS[0]);
      break;

    case ET_HEXA20:
      std::copy(&connect[0], &connect[0]+numElemNodes, &connectANSYS[0]);
      break;

    case ET_HEXA27:
      std::copy(&connect[0], &connect[0]+numElemNodes, &connectANSYS[0]);
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
    eltype.elementtypid = TypeNumber;
    std::fill(eltype.ielc, eltype.ielc+IELCSZ, 0);

    eltype.ielc[IETYP-1] = TypeNumber;

    eltype.ielc[JETYP-1] = PLANE42;
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
    eltype.elementtypid = TypeNumber;
    std::fill(eltype.ielc, eltype.ielc+IELCSZ, 0);

    eltype.ielc[IETYP-1] = TypeNumber;

    eltype.ielc[JETYP-1] = SOLID45;
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

    eltype.ielc[JETYP-1] = PLANE82;
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

    eltype.ielc[JETYP-1] = SOLID95;

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

  void SimOutputRST::MapInternal2ANSYSNodeDof(SolutionType solType)
  {
    std::string solName;
    Enum2String(solType, solName);
    UInt idx = ansysNodeDof2Idx_.size();

    // Check if solution type has already been added.
    if(internal2AnsysNodeDofMap_.find(solType) !=
       internal2AnsysNodeDofMap_.end())
      return;

    switch(solType)
    {
    case MECH_DISPLACEMENT:
      internal2AnsysNodeDofMap_[MECH_DISPLACEMENT] = UX;
      ansysNodeDof2Idx_[UX] = idx + 0;
      ansysNodeDof2Idx_[UY] = idx + 1;
      ansysNodeDof2Idx_[UZ] = idx + 2;
      LOG_TRACE(simOutputRST) << "Mapped " << solName
                              << " to UX UY UZ";
      break;
    case MECH_ACCELERATION:
      internal2AnsysNodeDofMap_[MECH_ACCELERATION] = AX;
      ansysNodeDof2Idx_[AX] = idx + 0;
      ansysNodeDof2Idx_[AY] = idx + 1;
      ansysNodeDof2Idx_[AZ] = idx + 2;
      LOG_TRACE(simOutputRST) << "Mapped " << solName
                              << " to AX AY AZ";
      break;
    case MECH_VELOCITY:
      internal2AnsysNodeDofMap_[MECH_VELOCITY] = VX;
      ansysNodeDof2Idx_[VX] = idx + 0;
      ansysNodeDof2Idx_[VY] = idx + 1;
      ansysNodeDof2Idx_[VZ] = idx + 2;
      LOG_TRACE(simOutputRST) << "Mapped " << solName
                              << " to VX VY VZ";
      break;
    case ELEC_POTENTIAL:
      internal2AnsysNodeDofMap_[ELEC_POTENTIAL] = VOLT;
      ansysNodeDof2Idx_[VOLT] = idx;
      LOG_TRACE(simOutputRST) << "Mapped " << solName
                              << " to VOLT";
      break;
    case ELEC_CHARGE:
      internal2AnsysNodeDofMap_[ELEC_CHARGE] = EMF;
      ansysNodeDof2Idx_[EMF] = idx;
      LOG_TRACE(simOutputRST) << "Mapped " << solName
                              << " to EMF";
      break;
    case ACOU_POTENTIAL:
      internal2AnsysNodeDofMap_[ACOU_POTENTIAL] = PRES;
      ansysNodeDof2Idx_[PRES] = idx;
      LOG_TRACE(simOutputRST) << "Mapped " << solName
                              << " to PRES";
      break;
    case ACOU_PRESSURE:
      internal2AnsysNodeDofMap_[ACOU_PRESSURE] = PRES;
      ansysNodeDof2Idx_[PRES] = idx;
      LOG_TRACE(simOutputRST) << "Mapped " << solName
                              << " to PRES";
      break;
    case ACOU_PRESSURE_DERIV_1:
      internal2AnsysNodeDofMap_[ACOU_PRESSURE_DERIV_1] = SP01;
      ansysNodeDof2Idx_[SP01] = idx;
      LOG_TRACE(simOutputRST) << "Mapped " << solName
                              << " to SP01";
      break;
    case ACOU_POTENTIAL_DERIV_1:
      internal2AnsysNodeDofMap_[ACOU_POTENTIAL_DERIV_1] = SP01;
      ansysNodeDof2Idx_[SP01] = idx;
      LOG_TRACE(simOutputRST) << "Mapped " << solName
                              << " to SP01";
      break;
    case MAG_POTENTIAL:
      internal2AnsysNodeDofMap_[MAG_POTENTIAL] = SP04;
      ansysNodeDof2Idx_[SP04] = idx;
      LOG_TRACE(simOutputRST) << "Mapped " << solName
                              << " to SP04 - 06";
      break;
    case HEAT_TEMPERATURE:
      internal2AnsysNodeDofMap_[HEAT_TEMPERATURE] = TEMP;
      ansysNodeDof2Idx_[TEMP] = idx;
      LOG_TRACE(simOutputRST) << "Mapped " << solName
                              << " to TEMP";
      break;
    case STOKESFLUID_VELOCITY:
      internal2AnsysNodeDofMap_[STOKESFLUID_VELOCITY] = VX;
      ansysNodeDof2Idx_[VX] = idx + 0;
      ansysNodeDof2Idx_[VY] = idx + 1;
      ansysNodeDof2Idx_[VZ] = idx + 2;
      LOG_TRACE(simOutputRST) << "Mapped " << solName
                              << " to VX VY VZ";
      break;
    case STOKESFLUID_PRESSURE:
      internal2AnsysNodeDofMap_[STOKESFLUID_PRESSURE] = PRES;
      ansysNodeDof2Idx_[PRES] = idx;
      LOG_TRACE(simOutputRST) << "Mapped " << solName
                              << "to PRES";
      break;
    default:
      LOG_TRACE(simOutputRST) << "Nodal result " << solName
                              << " not mapped to ANSYS any result type.";

      internal2AnsysNodeDofMap_[solType] = UNKN;
      break;
    }
  }
  
  void SimOutputRST::MapInternal2ANSYSElemDof(SolutionType solType)
  {
    std::string solName;
    Enum2String(solType, solName);

    switch(solType)
    {
    case MECH_STRESS:
      internal2AnsysElemDofMap_[MECH_STRESS] = ENS;
      LOG_TRACE(simOutputRST) << "Mapped " << solName
                              << " to stresses (ENS).";
      break;
    case MECH_STRAIN:
      internal2AnsysElemDofMap_[MECH_STRAIN] = EEL;
      LOG_TRACE(simOutputRST) << "Mapped " << solName
                              << " to elastic strains (EEL).";
      break;
    case ELEC_FIELD_INTENSITY:
      internal2AnsysElemDofMap_[ELEC_FIELD_INTENSITY] = EPL;
      LOG_TRACE(simOutputRST) << "Mapped " << solName
                              << " to plastic strains (EPL).";
      break;
    case MAG_FLUX_DENSITY:
      internal2AnsysElemDofMap_[MAG_FLUX_DENSITY] = ECR;
      LOG_TRACE(simOutputRST) << "Mapped " << solName
                              << " to creep strains (ECR).";
      break;
    case MAG_EDDY_CURRENT:
      internal2AnsysElemDofMap_[MAG_EDDY_CURRENT] = ETH;
      LOG_TRACE(simOutputRST) << "Mapped " << solName
                              << " to thermal strains (ETH).";
      break;
    default:
      LOG_TRACE(simOutputRST) << "Element result " << solName
                              << " not mapped to ANSYS any result type.";

      internal2AnsysElemDofMap_[solType] = UNKE;
      break;
    }
  }
  
} // end of namespace
