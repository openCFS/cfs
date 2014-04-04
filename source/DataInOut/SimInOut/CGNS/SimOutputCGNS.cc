// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <string>
#include <iomanip>

#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/exception.hpp>
namespace fs = boost::filesystem;

#include "DataInOut/Logging/LogConfigurator.hh"

#include "cgnslib.h"
#if CGNS_VERSION < 3100
# define cgsize_t int
#else
# if CG_BUILD_SCOPE
#  error enumeration scoping needs to be off
# endif
#endif

#include "SimOutputCGNS.hh"

namespace CoupledField {


  // ***************
  //   Constructor
  // ***************
  SimOutputCGNS::SimOutputCGNS( std::string& fileName,
                                PtrParamNode outputNode,
                                PtrParamNode infoNode,
                                bool isRestart ) 
    : SimOutput ( fileName, outputNode, infoNode, isRestart ) {
    
    formatName_ = "cgns";
    stepNumOffset_ = 0;
    stepValOffset_ = 0.0;
    
    capabilities_.insert( MESH );
    capabilities_.insert( MESH_RESULTS );

    std::string dirString = "results_" + formatName_; 
    myParam_->GetValue("directory", dirString, ParamNode::PASS );
    dirName_ = dirString; 
  }


  // **************
  //   Destructor
  // **************
  SimOutputCGNS::~SimOutputCGNS() {
    delete output;
  }


  // *************
  //   WriteGrid
  // *************
  void SimOutputCGNS::WriteGrid() {


    if ( !outputFileOK_) {
      EXCEPTION( "File for CGNS output results is not initialized" );
    }

    if (!ptGrid_) {
      EXCEPTION("Grid pointer is not initialized" );
    }

    int geodim=ptGrid_->GetDim(), physdim=geodim;
    strcpy(baseName_,"CFS_Simulation");
    cg_base_write(indexFile_,baseName_,geodim,physdim,&indexBase_);

    writeNodesAndElements();

  }

  void  SimOutputCGNS::writeNodesAndElements() {

    if (!ptGrid_)
      EXCEPTION("ptGrid_ is not initialized" );

    UInt geodim=ptGrid_->GetDim();
    UInt numNodes=ptGrid_->GetNumNodes();
    // only 3D meshes at the moment
    if (geodim != 3)
       EXCEPTION("only 3 dimensional models supported at the moment" );

    StdVector<RegionIdType> subdoms;
    StdVector<Elem*> elemssd;
    UInt numElems=ptGrid_->GetNumElems();
    //numElems = 0;
    //std::cout << "NumElems = " << numElems << std::endl;
    //for (i=0; i<subdoms.GetSize(); i++) {
    //  numElems += elemssd.GetSize();
    //  std::cout << "NumElems = " << numElems << std::endl;
    //}

    const Elem* ptEl;

    //std::cout << "NumElems1 = " << numElems << std::endl;
    UInt nelTET4=0, nelTRIA3=0, nelHEXA8=0, nelQUAD4=0;
    for (UInt iElem=0; iElem<numElems; iElem++) {
      ptEl = ptGrid_->GetElem(iElem+1);
      if (ptEl->type == Elem::ET_TET4)
        nelTET4++;
      if (ptEl->type == Elem::ET_TRIA3)
        nelTRIA3++;
      if (ptEl->type == Elem::ET_HEXA8)
        nelHEXA8++;
      if (ptEl->type == Elem::ET_QUAD4)
        nelQUAD4++;
    }
    // only tet models at the moment
    if (nelQUAD4 != 0 || nelHEXA8 != 0)
       EXCEPTION("only tet models supported at the moment" );

    numElems = nelTET4 + nelHEXA8;
    //std::cout << "NumElems2 = " << numElems << std::endl;

    // create zone
    strcpy(zoneName_,"CFS_Mesh");
    cgsize_t isize[3][1];
    isize[0][0] = numNodes;
    isize[1][0] = numElems;
    isize[2][0] = 0;
    cg_zone_write(indexFile_,indexBase_,zoneName_,isize[0],Unstructured,&indexZone_);

    // coordinates
    double *xCoord = new double[numNodes];
    double *yCoord = new double[numNodes];
    double *zCoord = new double[numNodes];

    Vector<Double> point;
    for ( UInt i = 0; i < numNodes; i++ ) {
      ptGrid_->GetNodeCoordinate(point,i+1);
      if ( geodim == 3 ) {
        xCoord[i] = point[0];
        yCoord[i] = point[1];
        zCoord[i] = point[2];
      } else if ( geodim == 2 ) {
        xCoord[i] = point[0];
        yCoord[i] = point[1];
      }
    }

    int indexCoord;
    // write grid coordinates
    cg_coord_write(indexFile_,indexBase_,indexZone_,RealDouble,"CoordinateX",xCoord,&indexCoord);
    cg_coord_write(indexFile_,indexBase_,indexZone_,RealDouble,"CoordinateY",yCoord,&indexCoord);
    cg_coord_write(indexFile_,indexBase_,indexZone_,RealDouble,"CoordinateZ",zCoord,&indexCoord);

    // for tet4 elements
    cgsize_t *elemData = new cgsize_t [numElems*4];

    StdVector<UInt> connect;
    std::string errMsg;

    UInt k = 0;
    UInt count = 0;
    for (UInt iElem=0; iElem<numElems; iElem++) {

      ptEl = ptGrid_->GetElem(iElem+1);
      if (ptEl->type == Elem::ET_TET4) {

        connect=ptEl->connect;

        for (UInt ii=0; ii < connect.GetSize(); ii++) { 
          elemData[count] = connect[ii]; count++;
        }
        k++;

      }
    }
    if (k != numElems) {
       EXCEPTION("Internal error in generation of cgns zone: element mismatch!");
    }

    int indexSection;
    int nelemStart = 1, nelemEnd = numElems, nbdyElem=0;
    cg_section_write(indexFile_,indexBase_,indexZone_,"Elements",TETRA_4,nelemStart,
                     nelemEnd,nbdyElem,elemData,&indexSection);


    if (0 == 0) {
    std::string solname = "DummyPressure1";
    double *press = new double[numNodes];
    for (UInt node=0; node<numNodes; node++) {
        press[node] = 1.0+xCoord[node]*(1.0 - yCoord[node]*yCoord[node])*exp(1.0-zCoord[node]);
    }
    int indexFlow,indexField;
    cg_sol_write(indexFile_,indexBase_,indexZone_,solname.c_str(),Vertex,&indexFlow);
    cg_field_write(indexFile_,indexBase_,indexZone_,indexFlow,RealDouble,
            "Pressure",press,&indexField);
    for (UInt node=0; node<numNodes; node++) {
        press[node] = 2.0+xCoord[node]*(1.0 - yCoord[node]*yCoord[node])*exp(1.0-zCoord[node]);
    }
    solname = "DummyPressure2";
    cg_sol_write(indexFile_,indexBase_,indexZone_,solname.c_str(),Vertex,&indexFlow);
    cg_field_write(indexFile_,indexBase_,indexZone_,indexFlow,RealDouble,
            "Pressure",press,&indexField);
    for (UInt node=0; node<numNodes; node++) {
        press[node] = 3.0+xCoord[node]*(1.0 - yCoord[node]*yCoord[node])*exp(1.0-zCoord[node]);
    }
    solname = "DummyPressure3";
    cg_sol_write(indexFile_,indexBase_,indexZone_,solname.c_str(),Vertex,&indexFlow);
    cg_field_write(indexFile_,indexBase_,indexZone_,indexFlow,RealDouble,
            "Pressure",press,&indexField);
// create BaseIterativeData
    int nsteps=3;
    cg_biter_write(indexFile_,indexBase_,"TimeIterValues",nsteps);
// go to BaseIterativeData level and write time values
    cg_goto(indexFile_,indexBase_,"BaseIterativeData_t",1,"end");
    cgsize_t nuse=3;
    double time[3];
    time[0] = 1.0; time[1] = 2.0; time[2] = 3.0;
    cg_array_write("TimeValues",RealDouble,1,&nuse,&time);
// create ZoneIterativeData
    cg_ziter_write(indexFile_,indexBase_,indexZone_,"ZoneIterativeData");
// go to ZoneIterativeData level and give info telling which
// flow solution corresponds with which time (solname(1) corresponds
// with time(1), solname(2) with time(2), and solname(3) with time(3))
    cg_goto(indexFile_,indexBase_,"Zone_t",indexZone_,"ZoneIterativeData_t",1,"end");
    cgsize_t idata[2];
    idata[0]=32;
    idata[1]=3;
    solname="DummyPressure1      ";
    solname+="DummyPressure2      ";
    solname+="DummyPressure3      ";
    char tmp3[3][33];
    char tmp[97];
    strcpy(tmp3[0],"DummyPressure1");
    strcpy(tmp3[1],"DummyPressure2");
    strcpy(tmp3[2],"DummyPressure3");
    sprintf(tmp, "%-32s%-32s%-32s",tmp3[0],tmp3[1],tmp3[2]);
    std::cout << "tmp = " << tmp << std::endl;
    cg_array_write("FlowSolutionPointers",Character,2,idata,tmp);
// add SimulationType
    cg_simulation_type_write(indexFile_,indexBase_,TimeAccurate);
    }

    cg_close(indexFile_);
  }

  void  SimOutputCGNS::NodeElemDataTransient(const UInt dataSetNr,
                                                  const std::string & title, 
                                                  const Vector<Double> & x, 
                                                  const UInt step, 
                                                  const Double time, 
                                                  const UInt nrNodes,
                                                  const UInt nrDofs)
  {
    //
    if (!ptGrid_)
      EXCEPTION("ptGrid_ is not initialized" );

    return;

    (*output) << std::setw(6) << -1 << std::endl 
              << std::setw(6) << dataSetNr << std::endl;

    (*output).setf(std::ios::scientific);
    (*output).precision(5);
    (*output).setf(std::ios::uppercase);

    // Standard for scalar values
    UInt dataCharac = 1;
    UInt valsPerNode = 1;
 
    // needed for undocumented value of
    // Dataset 55/56 in record8 , field4
    UInt specDataCharac = 0;

    // Vector type
    if (nrDofs > 1 && nrDofs <= 3) 
      {
        valsPerNode = 3;
        dataCharac = 2;
      } 
    // Tensor
    else if(nrDofs == 6) 
      { 
        valsPerNode = 6;
        specDataCharac = 2;
        dataCharac = 4; // symmetric tensor
        //      dataCharac = 3; // vector with 6 components
        //      dataCharac = 5; // unsymmetric tensor
      }
     
 
    (*output) << " " << title << " step" << std::setw(6) << step <<
      " time   " << time << std::endl;  
    (*output) << std::endl << std::endl << std::endl << std::endl;
    (*output) << std::setw(10) << 1 << std::setw(10) << 4 << std::setw(10) 
              << dataCharac  << std::setw(10) << specDataCharac
              << std::setw(10) << 2 << std::setw(10) << valsPerNode << std::endl;
    (*output) << std::setw(10) << 2 << std::setw(10) << 1 << std::setw(10) << 1 
              << std::setw(10) <<
      step << std::endl;
    (*output) << " " << time << std::endl;       

    UInt i,j,n;
    n=nrNodes;;  
    for (i=0; i<n; i++)
      {
     
        (*output) << std::setw(10) << i+1;
        if (dataSetNr == 56)
          (*output) << std::setw(10) << valsPerNode;

        (*output) << std::endl;
     
        // in the universal file either one or three results datas must exist
        if (nrDofs == 2)
          (*output) << std::setw(UNV_WIDTH) << 0.0;

        for (j=0; j<nrDofs; j++)
          {
            //std::cerr << "trying to write " << i << ", " << j << std::endl;
            (*output) << std::setw(UNV_WIDTH) << x[i*nrDofs +j];
          }
     
        (*output) << std::endl;
      }    
    (*output) << std::setw(6) << -1 << std::endl;
  }  

  void SimOutputCGNS::NodeElemDataHarmonic(const UInt dataSetNr,
                                                const std::string & title, 
                                                const Vector<Complex> & x, 
                                                const UInt step,
                                                const Double frequency, 
                                                const ComplexFormat format,
                                                const UInt nrNodes,
                                                const UInt nrDofs)
  {
  
    UInt dataCharact = 1;
    if (!ptGrid_)
      EXCEPTION("ptGrid_ is not initialized" );

    return;
  
    (*output) << std::setw(6) << -1 << std::endl 
              << std::setw(6) << dataSetNr << std::endl;
  
    (*output).setf(std::ios::scientific);
    (*output).precision(5);
    (*output).setf(std::ios::uppercase);
  
    UInt valsPerNode = 1;
    if (nrDofs > 1)
      {
        dataCharact = 2;
        valsPerNode = 3;
      }
 

    if (format == REAL_IMAG)
      {
        // write out realpart
        (*output) << " " << title << " cw realpart" << std::setw(6) <<" frequency   " << frequency << std::endl;  
        (*output) << std::endl << std::endl << std::endl << std::endl;
        (*output) << std::setw(10) << 1 << std::setw(10) << 5 << std::setw(10) << dataCharact << std::setw(10) << 0
                  << std::setw(10) << 2 << std::setw(10) << valsPerNode << std::endl;
        (*output) << std::setw(10) << 2 << std::setw(10) << 1 << std::setw(10) << -2 << std::setw(10) <<
          step << std::endl;
        (*output) << " " << frequency << std::endl;       
      
        UInt i,j,n;
        n=nrNodes;
        for (i=0; i<n; i++)
          {
            (*output) << std::setw(10) << i+1 << std::endl;
          
            // in the universal file either one or three results datas must exist
            if (nrDofs == 2)
              (*output) << std::setw(UNV_WIDTH) << 0.0;
          
            for (j=0; j<nrDofs; j++)
              {
                //std::cerr << "trying to write " << i << ", " << j << std::endl;
                (*output) << std::setw(UNV_WIDTH) << x[i*nrDofs +j].real();
              }
          
            (*output) << std::endl;
          }    
        (*output) << std::setw(6) << -1 << std::endl;

        // write out imag part
        (*output) << std::setw(6) << -1 << std::endl << std::setw(6) << 55 << std::endl;
        (*output) << " " << title << " cw imagpart" << std::setw(6) <<" frequency   " << frequency << std::endl;  
        (*output) << std::endl << std::endl << std::endl << std::endl;
        (*output) << std::setw(10) << 1 << std::setw(10) << 5 << std::setw(10) << dataCharact << std::setw(10) << 0
                  << std::setw(10) << 2 << std::setw(10) << valsPerNode << std::endl;
        (*output) << std::setw(10) << 2 << std::setw(10) << 1 << std::setw(10) << -12 << std::setw(10) <<
          step << std::endl;
        (*output) << " " << frequency << std::endl;       
      
        for (i=0; i<n; i++)
          {
            (*output) << std::setw(10) << i+1 << std::endl;
          
            // in the universal file either one or three results datas must exist
            if (nrDofs == 2)
              (*output) << std::setw(UNV_WIDTH) << 0.0;
          
            for (j=0; j<nrDofs; j++)
              {
                //std::cerr << "trying to write " << i << ", " << j << std::endl;
                (*output) << std::setw(UNV_WIDTH) << x[i*nrDofs +j].imag();
              }
          
            (*output) << std::endl;
          }    
        (*output) << std::setw(6) << -1 << std::endl;
      }
  
    else if (format == AMPLITUDE_PHASE) {
      // write out amplitude
      (*output) << " " << title << " cw amplitude" << std::setw(6) <<" frequency   " << frequency << std::endl;  
      (*output) << std::endl << std::endl << std::endl << std::endl;
      (*output) << std::setw(10) << 1 << std::setw(10) << 5 << std::setw(10) << dataCharact << std::setw(10) << 0
                << std::setw(10) << 2 << std::setw(10) << valsPerNode << std::endl;
      (*output) << std::setw(10) << 2 << std::setw(10) << 1 << std::setw(10) << -1 << std::setw(10) <<
        step << std::endl;
      (*output) << " " << frequency << std::endl;       
    
      UInt i,j,n;
      n=nrNodes;
      for (i=0; i<n; i++)
        {
          (*output) << std::setw(10) << i+1 << std::endl;
        
          // in the universal file either one or three results datas must exist
          if (nrDofs == 2)
            (*output) << std::setw(UNV_WIDTH) << 0.0;
        
          for (j=0; j<nrDofs; j++)
            {
              //std::cerr << "trying to write " << i << ", " << j << std::endl;
              (*output) << std::setw(UNV_WIDTH) << std::abs(x[i*nrDofs +j]);
            }
        
          (*output) << std::endl;
        }    
      (*output) << std::setw(6) << -1 << std::endl;
    
      // write out phase
      (*output) << std::setw(6) << -1 << std::endl << std::setw(6) << 55 << std::endl;
      (*output) << " " << title << " cw phase" << std::setw(6) <<" frequency   " << frequency << std::endl;  
      (*output) << std::endl << std::endl << std::endl << std::endl;
      (*output) << std::setw(10) << 1 << std::setw(10) << 5 << std::setw(10) << dataCharact << std::setw(10) << 0
                << std::setw(10) << 2 << std::setw(10) << valsPerNode << std::endl;
      (*output) << std::setw(10) << 2 << std::setw(10) << 1 << std::setw(10) << -11 << std::setw(10) <<
        step << std::endl;
      (*output) << " " << frequency << std::endl;       
    
      for (i=0; i<n; i++)
        {
          (*output) << std::setw(10) << i+1 << std::endl;
        
          // in the universal file either one or three results datas must exist
          if (nrDofs == 2)
            (*output) << std::setw(UNV_WIDTH) << 0.0;
        
          for (j=0; j<nrDofs; j++)
            {
              if (abs(x[i*nrDofs +j].imag()) > 1e-16)
                (*output) << std::setw(UNV_WIDTH) << CPhase(x[i*nrDofs +j]);
              else 
                (*output) << std::setw(UNV_WIDTH) << 0.0;
            }
        
          (*output) << std::endl;
        }    
      (*output) << std::setw(6) << -1 << std::endl;
    }
  
  }

  void SimOutputCGNS::Init(Grid * ptGrid, bool printGridOnly )
  {
    ptGrid_=ptGrid;

    // concatenate output file name
    try {
      fs::create_directory( dirName_ );
    } catch (std::exception &ex) {
      EXCEPTION(ex.what());
    }
    
    fs::path fileName (dirName_);
    fileName = fileName / (fileName_ + std::string(".cgns"));
    
    // Make sure a previous file with the same name gets removed.
    try {
      if(fs::exists(fileName))
      {
        fs::remove( fileName );
      }
    } catch (std::exception &ex) {
      EXCEPTION(ex.what());
    }

    // Set type of mid level library to use.
    if(myParam_->Has("mll")) 
    {
      int file_type = CG_FILE_HDF5;

      std::string mll = myParam_->Get("mll")->As<std::string>();

      if(mll == "adf") 
      {
        file_type = CG_FILE_ADF;
      }
      if(mll == "adf")
      {
        file_type = CG_FILE_ADF2;
      }
      
      cg_set_file_type(file_type);
    }    

    //     open CGNS file for write
    outputFileOK_ = false;
    if (cg_open(fileName.c_str(),CG_MODE_WRITE,&indexFile_)) cg_error_exit();
    outputFileOK_ = true;

    WriteGrid();
  }

  void SimOutputCGNS::RegisterResult( shared_ptr<BaseResult> sol,
                                     UInt saveBegin, UInt saveInc,
                                     UInt saveEnd,
                                     bool isHistory ) {
    
  }

  void SimOutputCGNS::BeginStep( UInt stepNum, Double stepVal ) {
    
    actStep_ = stepNum + stepNumOffset_;
    actStepVal_ = stepVal + stepValOffset_;
    resultMap_.clear();
  }

  void SimOutputCGNS::AddResult( shared_ptr<BaseResult>  sol ) {
    resultMap_[sol->GetResultInfo()->resultName].Push_back( sol );
  }
  
  void SimOutputCGNS::FinishStep( ) {

    // iterate over all result types
    ResultMapType::iterator it = resultMap_.begin();
    for( ; it != resultMap_.end(); it++ ) {
      
      // check if result is defined on nodes or elements
      ResultInfo & actInfo = *(it->second[0]->GetResultInfo());
      
      const StdVector<shared_ptr<BaseResult> > actResults =
        it->second;

      
      if(  actInfo.definedOn != ResultInfo::NODE &&
           actInfo.definedOn != ResultInfo::ELEMENT &&
           actInfo.definedOn != ResultInfo::SURF_ELEM ) {
        WARN( "Unv can only write results on element and nodes" );
        continue;
      }
      
      ResultInfo::EntityUnknownType entityType = actInfo.definedOn;
      std::string title =  SolutionTypeToString( actInfo.resultType );
      
      // if result could not be mapped, omit it
      if( title == "") {
        std::stringstream warning;
        warning  <<  "Result '" << actInfo.resultName
                 << "' coul not be mappted to unv result type and is "
                 << "omitted!";
        WARN( warning.str().c_str() );
        continue;
      }
      
      StdVector<std::string> & dofNames = actInfo.dofNames;
      UInt numDofs = dofNames.GetSize();
      
      // Determine type of entity the result is defined on
      UInt mapTo = 0;
      UInt numEntities = 0;
      if( actInfo.definedOn == ResultInfo::NODE ) {
        mapTo = 55;
        numEntities = ptGrid_->GetNumNodes();;
      } else {
        mapTo = 56;
        numEntities = ptGrid_->GetNumVolElems();
      }
      
      if( actResults[0]->GetEntryType() == BaseMatrix::DOUBLE ) {
        Vector<Double> gSol;
        FillGlobalVec<Double>(gSol, actResults, entityType );
        
        // Special case: If result is mechStress, restort entries 
        // according to capa notation
        if( actInfo.resultType == MECH_STRESS ) { 
          SortStresses( gSol, dofNames);
          numDofs = 6;
        }
        
        NodeElemDataTransient( mapTo, title, gSol, actStep_, 
                               actStepVal_ , numEntities ,
                               numDofs );
      } else {
        Vector<Complex> gSol;
        FillGlobalVec<Complex>(gSol, actResults, entityType );

        // Special case: If result is mechStress, restort entries 
        // according to capa notation
        if( actInfo.resultType == MECH_STRESS ) { 
          SortStresses( gSol, dofNames);
          numDofs = 6;
        }
        
        NodeElemDataHarmonic( mapTo, title, gSol, actStep_, 
                              actStepVal_, actInfo.complexFormat,
                              numEntities, numDofs );
      }
    } // over all result types
  }
  
  
  void SimOutputCGNS::FinishMultiSequenceStep() {

    // set offset for step value and number to last values
    stepNumOffset_ = actStep_;
    stepValOffset_ = actStepVal_;
  }
  
  
  template<class TYPE>
  void SimOutputCGNS::SortStresses( Vector<TYPE> & vec,

                                   const StdVector<std::string>& dofNames ){


    //soring according to capa (unv) notation
    //our notation is Voigt: Txx Tyy Tzz Tyz Txz Txy


    // check size of dofNames:
    // 6: 3D
    // 4: axi
    // 3: plane 

    // ensure, that number of entries is multiples of number of dofs
    assert( vec.GetSize() % dofNames.GetSize() == 0 );
    UInt numDofs = dofNames.GetSize();
    UInt numEntities = (UInt) vec.GetSize() / dofNames.GetSize();

    // CAPA always expects a stress vector with 6 entries
    Vector<TYPE> sorted;
    sorted.Resize( numEntities * 6 );
    
    // a) 3D-case (6 entries for stress vector)
    if( numDofs == 6 ) {
      
      //Txx Txy Tyy Txz Tyz Tzz
      UInt j = 0;
      for( UInt i = 0; i < vec.GetSize(); i+=6 ) {
        sorted[j++] = vec[i];
        sorted[j++] = vec[i+5];
        sorted[j++] = vec[i+1];
        sorted[j++] = vec[i+4];
        sorted[j++] = vec[i+3];
        sorted[j++] = vec[i+2];
      }
    }
  
    // b) axisymmetric-case (4 entries for stress vector)  
    else if( numDofs == 4 ) {
      //Tphiphi 0 Trr 0 Trz Tzz
      
      
      UInt j = 0;
      for( UInt i = 0; i < vec.GetSize(); i+=4 ) {
        sorted[j++] = vec[i+3];
        sorted[j++] = 0.0;
        sorted[j++] = vec[i+0];
        sorted[j++] = 0.0;
        sorted[j++] = vec[i+2];
        sorted[j++] = vec[i+1];
      }
    }

    // c) planestress/strain-case (3 entries for stress vector)  
    else if( numDofs == 3 ) {
      //0 0 Tyy 0 Tyz Tzz
    
      UInt j = 0;
      for( UInt i = 0; i < vec.GetSize(); i+=3 ) {
        sorted[j++] = 0.0;
        sorted[j++] = 0.0;
        sorted[j++] = vec[i+0];
        sorted[j++] = 0.0;
        sorted[j++] = vec[i+2];
        sorted[j++] = vec[i+1];
      } 

    } else {
      EXCEPTION( "Can not resort a stress vector with " << numDofs << " entries" );
    }
   
    // store sorted vector back to original one
    vec = sorted;
  }
  
  std::string SimOutputCGNS::SolutionTypeToString(const SolutionType type) const
  {

    //   std::string warnMsg;

    switch (type)
      {
      case MECH_DISPLACEMENT:
        return "displacement";
        break;
      case MECH_ACCELERATION:
        return "acceleration";
        break;
      case MECH_VELOCITY:
        return "velocity";
        break;
      case MECH_FORCE:
        break;
      case MECH_STRESS:
        return "stress";
        break;
      case MECH_STRAIN:
        EXCEPTION("Not implemented" );
        break;
      case MECH_RHS_LOAD:
        return "mechRhsLoad";
         break;
      case ELEC_POTENTIAL:
        return "electric potential";
        break;
      case ELEC_FIELD_INTENSITY:
        return "electric field";
        break;
      case ELEC_FORCE_VWP: 
        EXCEPTION("Not implemented" );
        break;
      case ELEC_CHARGE:
        return "electric charge";
        break;
      case ELEC_FLUX_DENSITY:
        EXCEPTION("Not implemented");
        break; 
      case ELEC_ENERGY:
        EXCEPTION("Not implemented");
        break;
      case ELEC_RHS_LOAD:
        return "elecRhsLoad";
        break;
      case SMOOTH_DISPLACEMENT:
        return "displacement";
        break;
      case ACOU_POTENTIAL:
        return "fluid potential";
        break;

      case ACOU_RHS_LOAD:
        return "fluid potential";
        break;
      case ACOU_PRESSURE:
        //       warnMsg = "Due to the restrictions in the .unv file format, the ";
        //       warnMsg += "acoustic pressure is written as acoustic (fluid) potential!";
        //       WARN(warnMsg.c_str(), __FILE__, __LINE__);
        return "fluid potential";
        break;
      case ACOU_FORCE:
        EXCEPTION("Not implemented" );
        break;
      case ACOU_POTENTIAL_DERIV_1:
        return "fluid potential, 1st deriv.";
        break;
      case ACOU_POTENTIAL_DERIV_2:
        return "fluid potential, 2nd deriv.";
        break;
      case MAG_POTENTIAL:
        return "mag. vector potential";
        break;
      case MAG_FLUX_DENSITY:
        return "mag. flux density";
        break;
      case MAG_FORCE_VWP:
        EXCEPTION("Not implemented");
        break;
      case MAG_FORCE_LORENTZ:
        EXCEPTION("Not implemented");
        break;
      case MAG_ENERGY:
        EXCEPTION("Not implemented");
        break;
      case MAG_RHS_LOAD:
        return "magRhsLoad";
        break;

      case HEAT_TEMPERATURE:
        return "temperature";
        break;
      case HEAT_RHS_LOAD:
        return "temperatureRhsLoad";
        break;

      default:
        WARN( "Wrong type of solution or 'SolutionType2String'" );
        return "UNKNOWN";
        break;
      }
    return std::string();
  }

}
