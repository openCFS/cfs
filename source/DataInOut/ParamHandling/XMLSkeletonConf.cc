// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <stdio.h>

#include <def_use_xerces.hh>

#include "Utils/tools.hh"
#include "General/environment.hh"
#include "Domain/GridCFS/grid_cfs.hh"
#include "Elements/elements_header.hh"
#include "DataInOut/programOptions.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ParamHandling/SkeletonConf.hh"
#include "DataInOut/WriteInfo.hh"

#define WL( MSG ) \
(*out_) << MSG << std::endl;

namespace CoupledField {


  // ***************
  //   Constructor
  // ***************
  SkeletonConf::SkeletonConf(SimInput* simInput) 
      : simInput_(simInput)
  {

    
    level_ = 0;

    std::string xmlFile = progOpts->GetParamFileStr();

    // just test, if config-file already exists  
    std::ifstream testfile( xmlFile.c_str() );

    std::string dummy;
    testfile >> dummy;
    if ( !dummy.empty() ) {
      EXCEPTION("XML-file is not empty: please change the name of your "
                << "current conf-file" << std::endl
                << "\t \t \t \t before calling cfs with option -skel");
    }

    // open the conf-file
    out_ = new std::ofstream( xmlFile.c_str() );
    if ( out_ == NULL ) {
      EXCEPTION("Could not open XML-file '" << xmlFile
               << "' for writing!");
    }

  }


  // **************
  //   Destructor
  // **************
  SkeletonConf::~SkeletonConf() {


    out_->close();

    std::cerr << std::endl;
    std::cerr << "\t Please complete the file before starting the simulation"
              << std::endl << std::endl;
  }

  void SkeletonConf::WriteConf ()
  {

    Info->StartProgress("Writing skeleton file to disc", false);

    WL("<?xml version=\"1.0\"?>");
    WL( "<cfsSimulation xmlns=\"http://www.cfs++.org\">");
    WL( "" );
    WriteGeneral();
    WriteSubdomains(); 

    // Start Writing of sequence step
    WL( "" );
    WL( "" );
    WL( Quote("Definition of one analysis/sequence step (static/transient/...)") );
    WL( Indent() << "<sequenceStep>" );
    WL( "" );
    WriteAnalysisTypes();
        
    WL( "" )
    WL( Quote( "PDE specific parameters" ) );
    WL( Indent() << "<pdeList>" );
    WL( "" )
    WritePDE();
    WL( Indent(-1) << "</pdeList>" );
    WL( Indent() << "" );
    WriteCouplingList();

    WL( Indent(-1) << "</sequenceStep>" );
    WL( Indent(-1)<< "</cfsSimulation>" );

  }



  void SkeletonConf::WriteGeneral ()
  {

    UInt dim=simInput_->GetDim();
    std::string geomType;

    if( dim == 2) {
      geomType = "plane";
    } else {
      geomType = "3d";
    }
    
  
    std::string headerLine
      = "=============================================================";

    WL( Quote(headerLine,1));
    WL( Quote("  SKELETON-XML-FILE: PLEASE REPLACE ALL \"XXX\" AND FILL OUT! ") );
    WL( Quote(headerLine) << std::endl );
    WL( Quote(" Define I/O file and format definitions ") );
    WL( Indent()   << "<fileFormats>" );
    WL( Indent(1)  << "<output>");
    WL( Indent(1)  << "<gid/>" );
    WL( Indent(-1) << "</output>" );
    WL( Indent()   << "<materialData file=\"mat.xml\" format=\"xml\"/>" );
    WL( Indent(-1) << "</fileFormats>" << std::endl );
    
  }




  // *******************
  //   WriteSubdomains
  // *******************
  void SkeletonConf::WriteSubdomains() {


    // Close the skeleton-config-file
    out_->close();

    // Generate parser and parse XML defaults file
    std::string cfsSchema = progOpts->GetSchemaPathStr();
    cfsSchema += "/CFS-Simulation/CFS.xsd";
    
    std::string cfsDefaults = progOpts->GetSchemaPathStr();
    cfsDefaults += "/CFS-Simulation/Defaults/CFS++Defaults.xml";

    // Generate dummy parameter node
    // Note: This parameter node is queried by the finite elements themselves
    // for getting the integration parameters.
    // This should be changed in the future
    param = PtrParamNode(new ParamNode());
    
    ptQ1   = new Quad1FE();
    ptQ2   = new Quad2FE();
    ptTet1 = new Tetra1FE();
    ptTet2 = new Tetra2FE();
    ptL1   = new Line1FE();
    ptL2   = new Line2FE();
    ptTr1  = new Triangle1FE();
    ptTr2  = new Triangle2FE();
    ptHexa1 = new Hexa1FE();
    ptHexa2 = new Hexa2FE();
    ptPyra1 = new Pyra1FE();
    ptPyra2 = new Pyra2FE();
    ptWedge1 = new Wedge1FE();
    ptWedge2 = new Wedge2FE();

    // now we can delete conf-object already
    param.reset();
    
    // Reopen skeleton-conf file
    std::string xmlFile = progOpts->GetSimName() + ".xml";
    out_->clear();
    out_->open( xmlFile.c_str(), std::ios_base::app );

    StdVector<std::string> regionNames;
    
    UInt dim=simInput_->GetDim();

    if (dim == 3)
      {
        //subdomains consists of 3d elements
        if (simInput_->GetNumElems(3) == 0)
          EXCEPTION( "3D-Problem specified, but no 3D-Elements in mesh-File" );

        simInput_->GetRegionNamesOfDim(regionNames, 3);
      }

    else if (dim == 2)
      {
        //subdomains consists of 2d elements
        if (simInput_->GetNumElems(2) == 0)
          EXCEPTION( "2D-Problem specified, but no 2D-Elements in mesh-File" );
        simInput_->GetRegionNamesOfDim(regionNames, 2);
      }
    else
      EXCEPTION("Dimension of Problem not supported");


    WL( Quote("Specify geometry and material information") );
    WL( Quote("'geometryType' must be one of plane/axi/3d" ) );
    WL( Indent() << "<domain geometryType=\"XXX\">" );
    WL( "" );
    WL( Quote("List of regions",1 ) );
    WL( Indent() << "<regionList>" );
    level_++;
    for (UInt i=0; i<regionNames.GetSize(); i++) {
      WL( Indent() << "<region name=\"" 
          << regionNames[i] << "\" material=\"XXX\"/>" );
    }
    WL( Indent(-1) << "</regionList>");

    WriteLists();
    WL( Indent(-1) << "</domain>");
  }






  void SkeletonConf::WriteLists ()
  {

    StdVector<std::string> nodeNames, surfRegionNames, elemNames;
    UInt dim = simInput_-> GetDim();

    // Print surface elements
    if (dim == 3){

      //check for 2D-interface elements
      if (simInput_->GetNumElems(2) != 0)
        simInput_->GetRegionNamesOfDim(surfRegionNames, 2);
   
    } else if (dim == 2) {
    
      //check for 1D-interface elements
      if (simInput_->GetNumElems(1) != 0)
        simInput_->GetRegionNamesOfDim(surfRegionNames, 1);
    }

    // === print surface regions ===
    if (surfRegionNames.GetSize()) {
      WL( "" );
      WL( Quote("List of surface / boundary regions" ) );
      WL( Indent() << "<surfRegionList>" );
      level_++;
      for (UInt i=0; i<surfRegionNames.GetSize(); i++) {
        WL( Indent() 
            << "<surfRegion name=\"" << surfRegionNames[i] 
            << "\"/>" );
      }
      WL( Indent(-1) << "</surfRegionList>" );
    }

    // === print named elements ===
    if (simInput_->GetNumNamedElems() != 0) {
      simInput_->GetElemNames(elemNames);
      WL( "" );
      WL( Quote("List of named elements" ) );
      WL( Indent() << "<elemList>" );
      level_++;
      for (UInt i=0; i<elemNames.GetSize(); i++) {
        WL( Indent() << "<elems name=\"" << elemNames[i] 
            << "\"/>" );
      }
      WL( Indent(-1) << "</elemList>" );
    }
    
    // === print named nodes ===
     if (simInput_->GetNumNamedNodes() != 0) {
       simInput_->GetNodeNames(nodeNames);
       WL( "" );
       WL( Quote("List of named nodes" ) );
       WL( Indent() << "<nodeList>" );
       level_++;
       for (UInt i=0; i<nodeNames.GetSize(); i++) {
         WL( Indent() << "<nodes name=\"" << nodeNames[i] 
             << "\"/>" );
       }
       WL( Indent(-1) << "</nodeList>" );
     }
  }

  
  void SkeletonConf::WritePDE () {
    
    WL( Quote( "Type of PDE (mechanic/electrostatic/...)",1) );
    WL( Indent(0) << "<XXX>" );
    WL( "" );

    WL( Quote( "List of regions the PDE is defined on",1 ) );
    WL( Indent(0) << "<regionList>" );
    WL( Quote( "Normal region definition",1 ) );
    WL( Indent(0) << "<region name=\"XXX\"/>" );
    WL( Quote( "Region with damping and nonlinearity" ) );
    WL( Indent(0) << "<region name=\"XXX\" dampingId=\"XXX\" nonLinId=\"XXX\"/>" );
    WL( Indent(-1) << "</regionList>" );
    WL( "" );


    WL( Quote( "In case of nonlinearities uncomment, the following list" ) );
    WL( Indent() << "<!-- <nonLinList>" );
    WL( Indent(1) << "<geometric id=\"XXX\"/>" );
    WL( Indent(-1) << "</nonLinList> -->" );
    WL( "" );

    WL( Quote( "In case of damping unqote uncomment, the following list" ) );
    WL( Indent() << "<!-- <dampingList>" );
    WL( Indent(1) << "<rayleigh id=\"XXX\"/>" );
    WL( Indent(-1) << "</dampingList> -->" );
    WL( "" );
    
    WL( Quote( "Definition of boundary conditions and loads" ) );
    WL( Indent() << "<bcsAndLoads>" );
    WL( Indent(1) << "<dirichletHom   name=\"XXX\" dof=\"XXX\"/>" );
    WL( Indent(0) << "<dirichletInhom name=\"XXX\" dof=\"XXX\""
        << " value=\"XXX\" phase=\"XXX\"/>" );
    WL( Indent(0) << "<load           name=\"XXX\" dof=\"XXX\""
        << " value=\"XXX\" phase=\"XXX\"/>" ); 
    WL( Indent(-1) <<  "</bcsAndLoads>" );
    WL( "" );

    WL( Quote( "Storing of results" ) );
    WL( Indent() << "<storeResults>" );

    WL( "" );
    WL( Indent(1) << "<nodeResult type=\"XXX\">" );
    WL( Indent(1) << "<allRegions/>" );
    WL( Indent(0) << "<nodeList>");
    WL( Indent(1) << "<nodes name=\"XXX\"/>");
    WL( Indent(-1) << "</nodeList>");
    WL( Indent(-1) << "</nodeResult>" );

    WL( "" );
    WL( Indent(0) << "<elemResult type=\"XXX\">" );
    WL( Indent(1) << "<allRegions/>" );
    WL( Indent(-1) << "</elemResult>" );
    
    WL( "" );
    WL( Indent(0) << "<regionResult type=\"XXX\">" );
    WL( Indent(1) << "<allRegions/>" );
    WL( Indent(0) << "<regionList>");
    WL( Indent(1) << "<region name=\"XXX\"/>");
    WL( Indent(-1) << "</regionList>");
    WL( Indent(-1) << "</regionResult>" );

    WL( Indent(-1) << "</storeResults>" );
    WL( Indent(-1) << "</XXX>" );

  }

  void SkeletonConf::WriteCouplingList () {

    WL( Quote( "For coupled simulation, uncomment the following lines" ) );
    WL( Indent() << "<!--<couplingList>" );
    WL( Indent(1) << "<direct>" );
    WL( Indent(1) << "<XXXDirect>" );
    WL( Indent(1) << "<regionList>" );
    WL( Indent(1) << "<region name=\"XXX\"/>" );
    WL( Indent(-1) << "</regionList>" );
    WL( Indent(0) << "<storeResults/>" );
    WL( Indent(-1) << "</XXXDirect>" );
    WL( Indent(-1) << "</direct>" );
    WL( "" );
    WL( Indent(0) << "<iterative>" );
    WL( Indent(1) << "<XXX1XXX2 method=\"RHS\">" );
    WL( Indent(1) << "<XXX1>" );
    WL( Indent(1) << "<coupling type=\"XXX\" quantity=\"XXX\" name=\"XXX\"/>" );
    WL( Indent(-1) << "</XXX1>" );
    WL( Indent(0) << "<XXX2>" );
    WL( Indent(1) << "<coupling type=\"XXX\" quantity=\"XXX\" name=\"XXX\"/>" );
    WL( Indent(-1) << "</XXX2>" );
    WL( Indent(-1) << "</XXX1XXX2>");
    WL( "");
    WL( Indent(0) << "<nonLinear logging=\"yes\">" );
    WL( Indent(1) << "<stopCrit value=\"1e-3\" quantity=\"XXX\" l2Norm=\"rel\"/>" );
    WL( Indent(0) << "<maxNumIters> 10 </maxNumIters> " );
    WL( Indent(-1) << "</nonLinear>" );
    WL( Indent(-1) << "</iterative>" );
    WL( "" );
    WL( Indent(-1) << "</couplingList> -->" );
    WL( "" );

    WL( Quote( "For special solver definitions, uncomment the following section" ) );
    WL( Indent(0) << "<!--<linearSystems>" );
    WL( Indent(1) << "<system name=\"XXX\">" );
    WL( Indent(1) << "<solver type=\"XXX\" precond=\"XXX\"/>" );
    WL( Indent(-1) << "</system>");
    WL( Indent(-1) << "</linearSystems> -->" );
    WL( "" );
  }

  void SkeletonConf::WriteAnalysisTypes(){

    WL( Quote( "Define type of analysis", 1 ) );
    WL( Quote( "Please uncomment the analysis of choice ") );
    WL( Indent(0) << "<analysis>" );
    
    // Static Analysis
    WL( Quote("static analysis",1 ) );
    WL( Indent(0) << "<!-- <static/> -->" );
    WL( Indent(-1) << "" );

    // Transient Analysis
    WL( Quote("transient analysis",1 ) );
    WL( Indent()   << "<!-- <transient>" );
    WL( Indent(1)  << "<numSteps> 1 </numSteps>" );
    WL( Indent()   << "<deltaT> 1 </deltaT>" );
    WL( Indent(-1) << "</transient> -->" );
    WL( Indent(-1) << "" );

    // Harmonic analysis
    WL( Quote("harmonic analysis",1 ) );
    WL( Indent()   << "<!-- <harmonic>" );
    WL( Indent(1)  << "<numFreq>   1 </numFreq>" );
    WL( Indent()   << "<startFreq> 1 </startFreq>" );
    WL( Indent()   << "<stopFreq>  1 </stopFreq>" );
    WL( Indent()   << "<sampling>  1 </sampling>" );
    WL( Indent(-1) << "</harmonic> -->" );
    WL( Indent(-1) << "" );

    // Eigenfrequency analysis
    WL( Quote("eigenFrequency analysis",1 ) );
    WL( Indent() << "<!-- <eigenFrequency>" );
    WL( Indent(1) << "<isQuadratic> no </isQuadratic>" );
    WL( Indent() << "<numModes>    10 </numModes>" );
    WL( Indent() << "<freqShift>   1  </freqShift>" );
    WL( Indent() << "<writeModes> yes </writeModes>" );
    WL( Indent(-1) << "</eigenFrequency> -->") ;
    WL( Indent(-1) << "</analysis>" );

  }

}


