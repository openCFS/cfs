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
    ENTER_FCN("SkeletonConf::SkeletonConf");

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

    Info->StartProgress("Reading in the mesh");
    //open the mesh-file
    strcpy(filename, aname);
    meshfile_ = new AnsysFile(filename);
    if (!meshfile_) Error("Can't open mesh-file");
    Info->FinishProgress();

  }


  SkeletonConf::~SkeletonConf ()
  {
    ENTER_FCN("SkeletonConf::~SkeletonConf");

    skelfile_->close();

    std::cerr << std::endl;
    std::cerr << "\t Please complete the file before starting the simulation" << std::endl << std::endl;

    delete [] name_;
  }

  void SkeletonConf::WriteConf ()
  {
    ENTER_FCN("SkeletonConf::WriteConf");

    Info->StartProgress("Writing skeleton file to disc");

    WriteGeneral();
    WriteSubdomains(); 
    //  WriteLists(); --> has to be done in WriteSubdomains!!

    (*skelfile_) << myendl << "   <!--  PDE SPECIFIC PARAMETERS -->" << myendl 
		 << "   <pdeList>" << myendl << myendl;
    WritePDE();
    (*skelfile_)  << "   </pdeList>" << std::endl << myendl;

    WriteCouplingList();

    (*skelfile_)  << "   <!--In case of transient analysis, uncomment following lines -->" << std::endl
		  << "   <!--<transient>  -->" << std::endl
		  << "   <!--   <numSteps>    1    </numSteps>    -->" << std::endl
		  << "   <!--   <firstDt>     1e-6 </firstDt>     -->" << std::endl
		  << "   <!--   <stepSaveBeg> 1    </stepSaveBeg> -->" << std::endl
		  << "   <!--   <stepSaveEnd> 1    </stepSaveEnd> -->" << std::endl
		  << "   <!--   <stepSaveInc> 1    </stepSaveInc> -->" << std::endl
		  << "   <!--   <timeDataFile name=\"XXX.dat\"/>    -->" << std::endl
		  << "   <!--</transient>                         -->" << std::endl 
		  << myEndl;

    (*skelfile_)  << "</cfsSimulation>" << myendl;

    Info->FinishProgress();
  }



  void SkeletonConf::WriteGeneral ()
  {
    ENTER_FCN("SkeletonConf::WriteGeneral");

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
    ENTER_FCN("SkeletonConf::WriteSubdomains");

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

    StdVector<std::string> sd;
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

    for (Integer i=0; i<sd.GetSize(); i++)
      (*skelfile_) << "      <region name=\"" << sd[i] << "\" material=\"XXX\"/>" << std::endl;

    (*skelfile_) << std::endl;

    WriteLists();
    
    (*skelfile_) << "   </domain>" << myEndl;
  }






  void SkeletonConf::WriteLists ()
  {
    ENTER_FCN("SkeletonConf::WriteLists");

    StdVector<std::string> sd;
    sd.Clear();
    Integer dim = meshfile_-> ReadDim();

    if (dim == 3){

      //check for 2D-interface elements
      if (meshfile_->GetNum2DElems() != 0)
	meshfile_->ReadEl2dConf(sd);
    
    } else if (dim == 2) {
    
      //check for 1D-interface elements
      if (meshfile_->GetNum1DElems() != 0) {
	meshfile_->ReadEl1dConf(sd);
      }
    }
    
    if (sd.GetSize()) {
      (*skelfile_) << "      <!--  LIST OF FACES -->" << std::endl;
      
      for (Integer i=0; i<sd.GetSize(); i++)
	(*skelfile_) << "      <elements name=\"" << sd[i] << "\"/>" << myendl;
      (*skelfile_) << myendl;
    }


    //check for node-list
    if (meshfile_->GetNumBCs() != 0)
      {
	sd.Clear();
	meshfile_->ReadBCsConf(sd);

	if (sd.GetSize())
	  (*skelfile_) << "      <!-- LIST OF NODES FOR BCs  --> " << myendl;
	
	for (Integer i=0; i<sd.GetSize(); i++)
	  (*skelfile_) << "      <nodes name=\"" << sd[i] << "\"/>" << myendl;
	(*skelfile_) << myendl;
      }


    if (meshfile_->GetNumSaveNodes() )
      {
	sd.Clear();
	meshfile_->ReadLevelOfSaveNodes(sd);
	if (sd.GetSize())
	  (*skelfile_) << "      <!-- LIST OF SAVE NODES --> " << std::endl;

	for (Integer i=0; i<sd.GetSize(); i++)
	  (*skelfile_) << "      <nodes name=\"" << sd[i] << "\"/>" << myendl;
      }
  }

  
  void SkeletonConf::WritePDE ()
  {
    ENTER_FCN("SkeletonConf::WritePDE");

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

  void SkeletonConf::WriteCouplingList ()
  {

    ENTER_FCN("SekeltonConf::WriteCoupling" );
    (*skelfile_) 
      << "   <!-- For coupled simulation, uncomment the following lines -->";
    (*skelfile_) 
      << myendl  
      << "   <!--<couplingList>                                  -->" 
      << myendl
      << "   <!--  <iterative>                                   -->" 
      << myendl 
      << "   <!--    <!-- pairwise coupling definition           -->" 
      << myendl
      << "   <!--    <XXX1XXX2 method=\"RHS\">                     -->" 
      << myendl
      << "   <!--      <XXX1>                                    -->" 
      << myendl
      << "   <!--        <coupling type=\"XXX\" quantity=\"XXX\"     -->" 
      << myendl
      << "   <!--                  name=\"XXX\"/>                  -->" 
      << myendl
      << "   <!--      </XXX1>                                   -->" 
      << myendl
      << "   <!--      <XXX2>                                    -->" 
      << myendl
      << "   <!--        <coupling type=\"XXX\" quantity=\"XXX\"     -->" 
      << myendl
      << "   <!--                   name=\"XXX\"/>                 -->" 
      << myendl
      << "   <!--      </XXX2>                                   -->" 
      << myendl 
      << "   <!--    </XXX1XXX2>                                 -->" 
      << myendl << myendl
      << "   <!--    <nonLinear logging=\"yes\">                   -->" 
      << myendl
      << "   <!--      <stopCrit value=\"1e-3\" quantity=\"XXX\"     -->" 
      << myendl
      << "   <!--                l2Norm=\"rel\"/                   -->" 
      << myendl << myendl 
      << "   <!--          <maxNumIters> 10 </maxNumIters>       -->" 
      << myendl
      << "   <!--    </nonLinear>                                -->" 
      << myendl
      << "   <!--  </iterative>                                  -->" 
      << myendl
      << "   <!--</couplingList>                                 -->" 
      << myendl << myendl;
  }
}

#endif


