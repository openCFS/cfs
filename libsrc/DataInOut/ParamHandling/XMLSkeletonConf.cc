#ifdef XMLPARAMS

#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <stdio.h>

#include "Utils/tools.hh"
#include "General/environment.hh"
#include "DataInOut/AnsysFile/ansysfile.hh"
#include "DataInOut/ParamHandling/ConfFile.hh"
#include "Domain/bcs.hh"
#include "Domain/GridCFS/grid_cfs.hh"
#include "Elements/elements_header.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "DataInOut/ParamHandling/PlainXMLParamHandler.hh"
#include "DataInOut/ParamHandling/XMLParamHandler.hh"
#include "DataInOut/ParamHandling/SkeletonConf.hh"

namespace CoupledField
{

  SkeletonConf::SkeletonConf (const Char * aname)
  {
    ENTER_FCN("Entering SkeletonConf::SkeletonConf");

    name_=new Char[100];
    strcpy(name_, aname);

    Char * filename = new Char[100];
    strcpy(filename, aname);
    //just test, if config-file already exists  
    std::ifstream testfile(strcat(filename,".xml"));

    std::string dummy;
    testfile >> dummy;
    if (!dummy.empty())
      {
	std::cerr << std::endl << "  \033[31mError\033[0m: " //<< std::endl 
		  << "conf-File is not empty: please change the name of your current conf-file" 
		  << std::endl
		  << "\t \t \t \t before calling cfs with option -skel"
		  << std::endl<< std::endl ;
	exit(1);

      }

    //open the conf-file
    strcpy(filename, aname);
    skelfile_ = new std::ofstream(strcat(filename,".xml"));
    if (!skelfile_) 
      Error("Can't open conf-file",__FILE__,__LINE__);

    std::cerr << std::endl << "  \033[31mCFS++\033[0m: " << std::endl;
    std::cerr << "\t " << filename << " is currently written to disc" << std::endl;

    //open the mesh-file
    strcpy(filename, aname);
    meshfile_ = new AnsysFile(filename);
    if (!meshfile_) Error("Can't open mesh-file");

  }


  SkeletonConf::~SkeletonConf ()
  {
    ENTER_FCN("Entering SkeletonConf::~SkeletonConf");

    skelfile_->close();


    std::cerr << "\t Please complete the file before starting the simulation" << std::endl << std::endl;

    delete [] name_;
  }

  void SkeletonConf::WriteConf ()
  {
    ENTER_FCN("Entering SkeletonConf::WriteConf");

    WriteGeneral();
    WriteSubdomains(); 
    //  WriteLists(); --> has to be done in WriteSubdomains!!


    (*skelfile_) << myendl << "   <!--  PDE SPECIFIC PARAMETERS -->" << myendl 
		 << "   <pdeList>" << myendl << myendl;
    WritePDE();
    (*skelfile_)  << "   </pdeList>" << std::endl << myendl;


    (*skelfile_)  << "   <!--In case of transient analysis, uncomment following lines -->" << std::endl
		  << "   <!--<transient>  -->" << std::endl
		  << "   <!--   <numsteps>    1    </numsteps>    -->" << std::endl
		  << "   <!--   <firstdt>     1e-6 </firstdt>     -->" << std::endl
		  << "   <!--   <stepsavebeg> 1    </stepsavebeg> -->" << std::endl
		  << "   <!--   <stepsaveend> 1    </stepsaveend> -->" << std::endl
		  << "   <!--   <stepsaveinc> 1    </stepsaveinc> -->" << std::endl
		  << "   <!--   <timeDataFile name=\"XXX.dat\"/>  -->" << std::endl
		  << "   <!--</transient>                         -->" << std::endl 
		  << myEndl;

    (*skelfile_)  << "</cfsSimulation>" << myendl;
  }


//     analysis     & plain   & required\\
//     geometry     & plain   & required\\
//     input        & plain   & optional\\
//     output       & plain   & optional\\
//     materialData & plain   & optional\\
//     domain       & section & required \\
//     pdeList      & section & required \\
//     transient    & section & required for transient analysis \\




  void SkeletonConf::WriteGeneral ()
  {
    ENTER_FCN("Entering SkeletonConf::WriteGeneral");

    (*skelfile_)  << "<?xml version=\"1.0\"?>" << myendl
		  << "<cfsSimulation xmlns=\"http://www.cfs++.org\">" << myendl << myendl;    

    (*skelfile_)  << "   <!-- ============================================================= -->" << std::endl 
		  << "   <!--  SKELETON-CONF-FILE: PLEASE REPLACE ALL \"XXX\" AND FILL OUT! -->" << std::endl
		  << "   <!-- ============================================================= -->" << std::endl
		  << myendl;

    (*skelfile_)  << "   <!--  ANALYSIS (static, transient, harmonic) -->" << std::endl
		  << "   <analysis type=\"XXX\"/>" << std::endl << std::endl;

    (*skelfile_)  << "   <!--  DEFINE GEOMETRY TYPE -->" << std::endl
		  << "   <geometry type=\"plane\"/>" << std::endl << std::endl;

    (*skelfile_)  << "   <!--  NAME OF MATERIAL FILE -->" << std::endl
		  << "   <materialData file=\"mat.dat\"/>" << std::endl << std::endl;
  }




  void SkeletonConf:: WriteSubdomains()
  {
    ENTER_FCN("Entering SkeletonConf::WriteSubdomains");

    //close the skeleton-config-file
    skelfile_->close();

    //we need to construct the conf-file object (needed by the FE-elements)
    Char * filename = new Char[100];
    strcpy(filename, name_);

    // Generate parser and parse XML defaults file
    std::string cfsDefaults = CVSEXTERNAL;
    cfsDefaults += "/CFS++XML/Defaults/CFS++Defaults.xml";

#ifdef USE_XERCES
    params = new XMLParamHandler( cfsDefaults.c_str() );
#else
    params = new PlainXMLParamHandler( cfsDefaults.c_str() );
#endif


    ptQ   = new Quad1FE();
    ptTet = new Tetra1FE();
    ptL1  = new Line1FE();
    ptL2  = new Line2FE();
    ptTr1 = new Triangle1FE();
    ptTr2 = new Triangle2FE();
    ptHexa = new Hexa1FE();
    ptPyra = new Pyra1FE();
    ptWedge = new Wedge1FE();
    // now we can delete conf-object already
    delete params;

    //reopen skeleton-conf file
    skelfile_->clear();
    skelfile_->open(strcat(filename,".xml"),std::ios_base::app);

    std::vector<std::string> sd;
    Integer dim = meshfile_-> ReadDim();
    if (dim == 3)
      {
	//subdomains consists of 3d elements
	if (meshfile_->GetNum3DElems() == 0)
	  Error("3D-Problem specified, but no 3D-Elements in mesh-File",__FILE__,__LINE__);

	meshfile_->ReadEl3dConf(sd);
      }

    else if (dim == 2)
      {
	//subdomains consists of 2d elements
	if (meshfile_->GetNum2DElems() == 0)
	  Error("2D-Problem specified, but no 2D-Elements in mesh-File",__FILE__,__LINE__);

	meshfile_->ReadEl2dConf(sd);
      }
    else
      Error("Dimension of Problem not supported",__FILE__,__LINE__);


    (*skelfile_) << "   <domain>" << myEndl;
    (*skelfile_) << "      <!-- LIST OF SUBDOMAINS -->"<< std::endl;

    for (Integer i=0; i<sd.size(); i++)
      (*skelfile_) << "      <region name=\"" << sd[i] << "\" material=\"XXX\"/>" << std::endl;

    (*skelfile_) << std::endl;

    WriteLists();
    
    (*skelfile_) << "   </domain>" << myEndl;
  }






  void SkeletonConf::WriteLists ()
  {
    ENTER_FCN("Entering SkeletonConf::WriteLists");

    std::vector<std::string> sd;
    sd.clear();
    Integer dim = meshfile_-> ReadDim();

    if (dim == 3)

      //check for 2D-interface elements
      if (meshfile_->GetNum2DElems() != 0)
	meshfile_->ReadEl2dConf(sd);
    
    else if (dim == 2)    
      //check for 1D-interface elements
      if (meshfile_->GetNum1DElems() != 0)
	meshfile_->ReadEl1dConf(sd);

    if (sd.size())
      (*skelfile_) << "      <!--  LIST OF FACES -->" << std::endl;

    for (Integer i=0; i<sd.size(); i++)
      (*skelfile_) << "      <elements name=\"" << sd[i] << "\"/>" << myendl;


    //check for node-list
    if (meshfile_->GetNumBCs() != 0)
      {
	sd.clear();
	meshfile_->ReadBCsConf(sd);

	if (sd.size())
	  (*skelfile_) << "      <!-- LIST OF NODES FOR BCs  --> " << myendl;
	
	for (Integer i=0; i<sd.size(); i++)
	  (*skelfile_) << "      <nodes name=\"" << sd[i] << "\"/>" << myendl;
      }


    if (meshfile_->GetNumSaveNodes() )
      {
	sd.clear();
	meshfile_->ReadLevelOfSaveNodes(sd);
	if (sd.size())
	  (*skelfile_) << "      <!-- LIST OF SAVE NODES --> " << std::endl;

	for (Integer i=0; i<sd.size(); i++)
	  (*skelfile_) << "      <nodes name=\"" << sd[i] << "\"/>" << myendl;
      }
  }

  
  void SkeletonConf::WritePDE ()
  {
    ENTER_FCN("Entering SkeletonConf::WritePDE");

    (*skelfile_) << "      <!-- name of pde -->" << std::endl;
    (*skelfile_) << "      <XXX>" << std::endl 
		 << std::endl;       
    (*skelfile_) << "         <region name=\"XXX\"/>" << std::endl 
		 << std::endl;       
    (*skelfile_) << "         <!-- boundary conditions -->" << std::endl;
    (*skelfile_) << "         <bcsAndLoads>" << std::endl; 
    (*skelfile_) << "            <dirichletHom   name=\"XXX\"/>" << std::endl;    
    (*skelfile_) << "            <dirichletInhom name=\"XXX\" value=\"1\" />" << std::endl; 
    (*skelfile_) << "            <load           name=\"XXX\" value=\"1\" />" << std::endl; 
    (*skelfile_) << "         </bcsAndLoads>" << std::endl 
		 << myendl;       
    (*skelfile_) << "      </XXX>" << std::endl << std::endl;
  }
}

#endif


