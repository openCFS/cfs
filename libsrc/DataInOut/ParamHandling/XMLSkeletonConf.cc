#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <stdio.h>

#include "Utils/tools.hh"
#include "General/environment.hh"
#include "DataInOut/AnsysFile/ansysfile.hh"
#include "Domain/GridCFS/grid_cfs.hh"
#include "Elements/elements_header.hh"
#include "DataInOut/CommandLine/BaseCommandLineHandler.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "DataInOut/ParamHandling/PlainXMLParamHandler.hh"
#include "DataInOut/ParamHandling/XMLParamHandler.hh"
#include "DataInOut/ParamHandling/SkeletonConf.hh"

namespace CoupledField {


  // ***************
  //   Constructor
  // ***************
  SkeletonConf::SkeletonConf() {

    ENTER_FCN( "SkeletonConf::SkeletonConf" );

    std::string xmlFile = commandLine->GetParamFile();
    std::string meshFile = commandLine->GetMeshFile();

    // just test, if config-file already exists  
    std::ifstream testfile( xmlFile.c_str() );

    std::string dummy;
    testfile >> dummy;
    if ( !dummy.empty() ) {
      std::cerr << std::endl << "  \033[31mError\033[0m: " //<< std::endl 
                << "conf-File is not empty: please change the name of your "
                << "current conf-file" << std::endl
                << "\t \t \t \t before calling cfs with option -skel"
                << std::endl<< std::endl ;
      exit(1);
    }

    // open the conf-file
    skelfile_ = new std::ofstream( xmlFile.c_str() );
    if ( skelfile_ == NULL ) {
      (*error) << "Could not open XML-file '" << xmlFile
               << "' for writing!\n";
      Error( __FILE__, __LINE__ );
    }

    Info->StartProgress("Reading in the mesh");

    // Open the mesh-file
    meshfile_ = new AnsysFile( meshFile.c_str() );
    if ( meshfile_ == NULL ) {
      (*error) << "Could not open mesh-file '" << xmlFile
               << "' for reading!\n";
      Error( __FILE__, __LINE__ );
    }

    Info->FinishProgress();
  }


  // **************
  //   Destructor
  // **************
  SkeletonConf::~SkeletonConf() {

    ENTER_FCN( "SkeletonConf::~SkeletonConf" );

    skelfile_->close();

    std::cerr << std::endl;
    std::cerr << "\t Please complete the file before starting the simulation"
              << std::endl << std::endl;
  }

  void SkeletonConf::WriteConf ()
  {
    ENTER_FCN("SkeletonConf::WriteConf");

    Info->StartProgress("Writing skeleton file to disc", FALSE);

    WriteGeneral();
    WriteSubdomains(); 
    //  WriteLists(); --> has to be done in WriteSubdomains!!

    (*skelfile_) <<myendl<< "   <!--  PDE SPECIFIC PARAMETERS -->" << myendl
		 << "   <pdeList>" << myendl << myendl;
    WritePDE();
    (*skelfile_)  << "   </pdeList>" << std::endl << myendl;

    WriteCouplingList();

    WriteAnalysisTypes();

    (*skelfile_)  << "</cfsSimulation>" << myendl;

  }



  void SkeletonConf::WriteGeneral ()
  {
    ENTER_FCN("SkeletonConf::WriteGeneral");

    (*skelfile_)  << "<?xml version=\"1.0\"?>" << myendl
                  << "<cfsSimulation xmlns=\"http://www.cfs++.org\">"
                  << myendl << myendl;    

    (*skelfile_)  << "   <!-- ==============================================="
                  << "============== -->" << std::endl 
                  << "   <!--   SKELETON-CONF-FILE: PLEASE REPLACE ALL \"XXX\""
                  << " AND FILL OUT!  -->" << std::endl
                  << "   <!-- ==============================================="
                  << "============== -->" << std::endl
                  << myendl;

    (*skelfile_)  << "   <!--  ANALYSIS (static, transient, harmonic, "
                  << "multiSequence) -->" << std::endl
                  << "   <analysis type=\"XXX\"/>" << std::endl << std::endl;

    (*skelfile_)  << "   <!--  DEFINE GEOMETRY TYPE (plane, axi, 3d)-->"
                  << std::endl
                  << "   <geometry type=\"plane\"/>" << std::endl << std::endl;

    (*skelfile_)  << "   <!--  NAME OF MATERIAL FILE -->" << std::endl
                  << "   <materialData file=\"mat.dat\"/>" << std::endl
                  << std::endl;
  }




  // *******************
  //   WriteSubdomains
  // *******************
  void SkeletonConf:: WriteSubdomains() {

    ENTER_FCN( "SkeletonConf::WriteSubdomains" );

    // Close the skeleton-config-file
    skelfile_->close();

    // Generate parser and parse XML defaults file
    std::string cfsDefaults = XMLSCHEMA;
    cfsDefaults += "/Defaults/CFS++Defaults.xml";
#ifdef USE_XERCES
    params = new XMLParamHandler( cfsDefaults.c_str() );
#else
    params = new PlainXMLParamHandler( cfsDefaults.c_str() );
#endif


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
    delete params;

    // Reopen skeleton-conf file
    std::string xmlFile = commandLine->GetSimName() + ".xml";
    skelfile_->clear();
    skelfile_->open( xmlFile.c_str(), std::ios_base::app );

    StdVector<std::string> regionNames;
    
    UInt dim=meshfile_->GetDim();

    if (dim == 3)
      {
        //subdomains consists of 3d elements
        if (meshfile_->GetNumElems(3) == 0)
          Error( "3D-Problem specified, but no 3D-Elements in mesh-File",
                 __FILE__,__LINE__);

        meshfile_->GetRegionNamesOfDim(regionNames, 3);
      }

    else if (dim == 2)
      {
        //subdomains consists of 2d elements
        if (meshfile_->GetNumElems(2) == 0)
          Error( "2D-Problem specified, but no 2D-Elements in mesh-File",
                 __FILE__,__LINE__);
        meshfile_->GetRegionNamesOfDim(regionNames, 2);
      }
    else
      Error("Dimension of Problem not supported",__FILE__,__LINE__);


    (*skelfile_) << "   <domain>" << myEndl;
    (*skelfile_) << "      <!-- LIST OF SUBDOMAINS -->"<< std::endl;

    for (UInt i=0; i<regionNames.GetSize(); i++)
      (*skelfile_) << "      <region name=\"" << regionNames[i]
                   << "\" material=\"XXX\"/>" << std::endl;

    (*skelfile_) << std::endl;

    WriteLists();
    
    (*skelfile_) << "   </domain>" << myEndl;
  }






  void SkeletonConf::WriteLists ()
  {
    ENTER_FCN("SkeletonConf::WriteLists");

    StdVector<std::string> nodeNames, surfRegionNames;
    UInt dim = meshfile_-> GetDim();

   


    //check for node-list
    if (meshfile_->GetNumNamedNodes() != 0)
      {
        meshfile_->GetNodeNames(nodeNames);

        if (nodeNames.GetSize())
          (*skelfile_) << "      <!-- LIST OF NODES --> " << myendl;
        
        for (UInt i=0; i<nodeNames.GetSize(); i++)
          (*skelfile_) << "      <nodes name=\"" << nodeNames[i] 
                       << "\"/>" << myendl;
        (*skelfile_) << myendl;
      }


    // Print surface elements
    if (dim == 3){

      //check for 2D-interface elements
      if (meshfile_->GetNumElems(2) != 0)
        meshfile_->GetRegionNamesOfDim(surfRegionNames, 2);
   
    } else if (dim == 2) {
    
      //check for 1D-interface elements
      if (meshfile_->GetNumElems(1) != 0)
        meshfile_->GetRegionNamesOfDim(surfRegionNames, 1);
    }
    if (surfRegionNames.GetSize()) {
      (*skelfile_) << "      <!--  LIST OF SURFACE ELEMENTS -->" << std::endl;
      
      for (UInt i=0; i<surfRegionNames.GetSize(); i++)
        (*skelfile_) << "      <elements name=\"" << surfRegionNames[i] 
                     << "\"/>" << myendl;
      (*skelfile_) << myendl;
    }
  }

  
  void SkeletonConf::WritePDE ()
  {
    ENTER_FCN("SkeletonConf::WritePDE");

    (*skelfile_) << "      <!-- name of pde -->" << std::endl
      << "      <XXX>" << std::endl 
      << std::endl
      << "         <region name=\"XXX\"/>" << std::endl 
      << std::endl
      << "         <!-- boundary conditions -->" << std::endl
      << "         <bcsAndLoads>" << std::endl
      << "            <dirichletHom   name=\"XXX\"/>" << std::endl
      << "            <dirichletInhom name=\"XXX\" value=\"1\" />"
      << std::endl
      << "            <load           name=\"XXX\" value=\"1\" />"
      << std::endl 
      << "         </bcsAndLoads>" << std::endl 
      << myendl
      << "         <!-- storing of results -->" << std::endl
      << "         <storeResults>" << std::endl
      << "           <!-- <nodeResults type=\"XXX\"/>                 -->"  
      << std::endl
      << "           <!-- <nodeHistory type=\"XXX\" saveNodes=\"XXX\"/> -->"
      << std::endl
      << "           <!-- <elemResults type=\"XXX\" region=\"XXX\"/>    -->" 
      << std::endl
      << "         </storeResults>" << std::endl
      << "      </XXX>" << std::endl << std::endl;
  }

  void SkeletonConf::WriteCouplingList ()
  {

    ENTER_FCN("SekeltonConf::WriteCoupling" );
    (*skelfile_) 
      << "   <!-- For coupled simulation, uncomment the following lines -->";
    (*skelfile_) 
      << myendl  
      << "   <!--<couplingList>                                     " 
      << myendl
      << "     <iterative>                                      " 
      << myendl 
      << "       <XXX1XXX2 method=\"RHS\">                        " 
      << myendl
      << "         <XXX1>                                       " 
      << myendl
      << "           <coupling type=\"XXX\" quantity=\"XXX\"        " 
      << myendl
      << "                     name=\"XXX\"/>                     " 
      << myendl
      << "         </XXX1>                                      " 
      << myendl
      << "         <XXX2>                                       " 
      << myendl
      << "           <coupling type=\"XXX\" quantity=\"XXX\"        " 
      << myendl
      << "                     name=\"XXX\"/>                    " 
      << myendl
      << "         </XXX2>                                      " 
      << myendl 
      << "       </XXX1XXX2>                                    " 
      << myendl << myendl
      << "       <nonLinear logging=\"yes\">                      " 
      << myendl
      << "         <stopCrit value=\"1e-3\" quantity=\"XXX\"        " 
      << myendl
      << "                   l2Norm=\"rel\"/                      " 
      << myendl << myendl 
      << "         <maxNumIters> 10 </maxNumIters>          " 
      << myendl
      << "       </nonLinear>                                   " 
      << myendl
      << "     </iterative>                                     " 
      << myendl
      << "   </couplingList>                                 -->" 
      << myendl << myendl << myendl;
  }

  void SkeletonConf::WriteAnalysisTypes(){
    ENTER_FCN( "SkeletonConf::WriteAnalysisTypes" );



    // Transient Analysis
    (*skelfile_)  << "   <!--In case of transient analysis, uncomment "
                  << "following lines -->" << std::endl
                  << "   <!--<transient>    " << std::endl
                  << "     <numSteps>    1    </numSteps>       " << std::endl
                  << "     <firstDt>     1e-6 </firstDt>        " << std::endl
                  << "     <stepSaveBeg> 1    </stepSaveBeg>    " << std::endl
                  << "     <stepSaveEnd> 1    </stepSaveEnd>    " << std::endl
                  << "     <stepSaveInc> 1    </stepSaveInc>    " << std::endl
                  << "     <timeDataFile name=\"XXX.dat\"/>     " << std::endl
                  << "   </transient>                        -->" << std::endl 
                  << myEndl;

    // Harmonic Analysis
    (*skelfile_)  << "   <!--In case of harmonic analysis, uncomment "
                  << "following lines -->" << std::endl
                  << "   <!--<harmonic>    " << std::endl
                  << "      <numFreq>    1   </numFreq>       " << std::endl
                  << "      <startFreq>  1e3 </startFreq>     " << std::endl
                  << "      <stopFreq>   1   </stopFreq>    " << std::endl
                  << "   </harmonic>                         -->" << std::endl 
                  << myEndl;

    // MultiSequence Analysis
    (*skelfile_)  << "   <!--In case of multiSequence analysis, uncomment "
                  << "following lines -->" << std::endl
                  << "   <!--<multiSequence>       " << std::endl
                  << "     <step index=\"1\"> " << std::endl
                  << "       <pde refTag=\"XXX\" type=\"XXX\" " 
                  << "analysis=\"XXX\"/>" << std::endl
                  << "     </step>               " << std::endl
                  << "     <step index=\"2\"> " << std::endl
                  << "       <pde refTag=\"XXX\" type=\"XXX\" " 
                  << "analysis=\"XXX\"/>" << std::endl
                  << "     </step>               " << std::endl
         << "   </multiSequence>                            -->" << std::endl 
                  << myEndl;


  }

}


