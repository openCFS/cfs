#ifndef XMLPARAMS

#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <stdio.h>
#include <Utils/tools.hh>
#include <General/environment.hh>
#include <DataInOut/AnsysFile/ansysfile.hh>
#include <DataInOut/conffile.hh>
#include <Domain/bcs.hh>
#include <Domain/GridCFS/grid_cfs.hh>
#include <Elements/elements_header.hh>

#include "SkeletonConf.hh"


namespace CoupledField
{

  SkeletonConf::SkeletonConf (const Char * aname)
  {
#ifdef TRACE
    (*trace) << "Entering SkeletonConf::SkeletonConf" << std::endl;
#endif

    name_=new Char[100];
    strcpy(name_, aname);

    Char * filename = new Char[100];
    strcpy(filename, aname);
    //just test, if config-file already exists  
    std::ifstream testfile(strcat(filename,".conf"));

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
    skelfile_ = new std::ofstream(strcat(filename,".conf"));
    if (!skelfile_) 
      Error("Can't open conf-file",__FILE__,__LINE__);

    std::cerr << std::endl << "  \033[31mCFS++\033[0m: " << std::endl;
    //    strcat(filename, ".conf");
    std::cerr << "\t " << filename << " is currently written to disc" << std::endl;

    //open the mesh-file
    strcpy(filename, aname);
    meshfile_ = new AnsysFile(filename);
    if (!meshfile_) Error("Can't open mesh-file");

  }


  SkeletonConf::~SkeletonConf ()
  {
#ifdef TRACE
    (*trace) << "Entering SkeletonConf::~SkeletonConf" << std::endl;
#endif

    skelfile_->close();


    std::cerr << "\t Please complete the file before starting the simulation" << std::endl << std::endl;

    delete [] name_;
  }

  void SkeletonConf::WriteConf ()
  {
#ifdef TRACE
    (*trace) << "Entering SkeletonConf::WriteConf" << std::endl;
#endif

    WriteGeneral();
    WriteSubdomains();
    WriteLists();
    WritePDE();
  }

  void SkeletonConf::WriteGeneral ()
  {
#ifdef TRACE
    (*trace) << "Entering SkeletonConf::WriteGeneral" << std::endl;
#endif

    (*skelfile_)  << "#" << std::endl 
		  << "# -- SKELETON-CONF-FILE: PLEASE REPLACE ALL \"XXX\" AND FILL OUT!" << std::endl
		  << "#" << std::endl << std::endl;

    (*skelfile_)  << "# -- NAME OF MATERIAL FILE" << std::endl
		  << "material_file = mat.dat" << std::endl << std::endl;
    (*skelfile_)  << "# -- ANALYSIS (static, transient, harmonic)" << std::endl
		  << "analysis = XXX" << std::endl << std::endl;
    (*skelfile_)  << "# -- TIME DATA FILE (if transient analysis, uncomment following line)" << std::endl
		  << "#time_data_files = XXX non" << std::endl << std::endl;

  }

  void SkeletonConf:: WriteSubdomains()
  {
#ifdef TRACE
    (*trace) << "Entering SkeletonConf::WriteSubdomains" << std::endl;
#endif

    //close the skeleton-config-file
    skelfile_->close();

    //we need to construct the conf-file object (needed by the FE-elements)
    Char * filename = new Char[100];
    strcpy(filename, name_);
    conf = new ConfFile(filename);
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
    delete conf;

    //reopen skeleton-conf file
    skelfile_->clear();
    skelfile_->open(strcat(filename,".conf"),std::ios_base::app);

    std::vector<std::string> sd;
    Integer dim = meshfile_-> ReadDim();
    if (dim == 3)
      {
	//subdomains consists of 3d elements
	if (meshfile_->GetNum3DElems() == 0)
	  Error("3D-Problem specified, but no 3D-Elements in mesh-File",__FILE__,__LINE__);

	meshfile_->ReadEl3dConf(sd);
	(*skelfile_) << "# -- SUBDOMAINS" 
		     << std::endl << "subdomains = " << sd.size() << std::endl << std::endl;

	(*skelfile_) << "# -- LIST OF SUBDOMAINS (for each subdomain, specify the material name)" 
		     << std::endl;
	(*skelfile_) << "list_subdomains: " << std::endl;
	for (Integer i=0; i<sd.size(); i++)
	  (*skelfile_) << "\t" << sd[i] << " = XXX" << std::endl;
	(*skelfile_) << std::endl;
      }

    else if (dim == 2)
      {
	//subdomains consists of 2d elements
 	std::vector<std::string> sd;
	if (meshfile_->GetNum2DElems() == 0)
	  Error("2D-Problem specified, but no 2D-Elements in mesh-File",__FILE__,__LINE__);

	meshfile_->ReadEl2dConf(sd);
	(*skelfile_) << "# -- SUBDOMAINS" 
		     << std::endl << "subdomains = " << sd.size() << std::endl << std::endl;

	(*skelfile_) << "# -- LIST OF SUBDOMAINS (SPECIFY FOR EACH SUBDOMAIN THE MATERIAL NAME)" 
		     << std::endl;
	(*skelfile_) << "list_subdomains: " << std::endl;
	for (Integer i=0; i<sd.size(); i++)
	  (*skelfile_) << "\t" << sd[i] << " = XXX" << std::endl;
	(*skelfile_) << std::endl;
      }
    else
      Error("Dimension of Problem not supported",__FILE__,__LINE__);
  }






  void SkeletonConf::WriteLists ()
  {
#ifdef TRACE
    (*trace) << "Entering SkeletonConf::WriteLists" << std::endl;
#endif

    std::vector<std::string> sd;
    Integer dim = meshfile_-> ReadDim();
    if (dim == 3)
      {
	//check for 2D-interface elements
	if (meshfile_->GetNum2DElems() != 0)
	  {
	    sd.clear();
	    meshfile_->ReadEl2dConf(sd);
	    (*skelfile_) << "# -- LIST OF FACES" 
			 << std::endl << "list_faces = " ;
	    for (Integer i=0; i<sd.size(); i++)
	      (*skelfile_) << sd[i] << " ";
	    (*skelfile_) << "non" << std::endl << std::endl;
	  }
      }
    
    else if (dim == 2)
      {
	//check for 1D-interface elements
	if (meshfile_->GetNum1DElems() != 0)
	  {
	    sd.clear();
	    meshfile_->ReadEl1dConf(sd);
	    (*skelfile_) << "# -- LIST OF EDGES" 
			 << std::endl << "list_edges = " ;
	    for (Integer i=0; i<sd.size(); i++)
	      (*skelfile_) << sd[i] << " ";
	    (*skelfile_) << "non" << std::endl << std::endl;
	  }
      }

    //check for node-list
    if (meshfile_->GetNumBCs() != 0)
      {
	sd.clear();
	meshfile_->ReadBCsConf(sd);
	(*skelfile_) << "# -- LIST OF NODES (has to be finished by \"non\") " 
		     << std::endl << "list_nodes = " ;
	for (Integer i=0; i<sd.size(); i++)
	  (*skelfile_) << sd[i] << " ";
	(*skelfile_) << "non" << std::endl << std::endl;
      }

    //history nodes
    (*skelfile_) << "# -- SPECIFY HISTORY NODES (finish list with \"-1\") " << std::endl;
    (*skelfile_) << "history_node = -1" << std::endl << std::endl;

    if (meshfile_->GetNumSaveNodes() )
      {
	sd.clear();
	meshfile_->ReadLevelOfSaveNodes(sd);
	(*skelfile_) << "# -- LIST OF SAVE NODE LEVELS (has to be finished by \"non\") "
		     << std::endl << "save_nodes = " ;
	for (Integer i=0; i<sd.size(); i++)
	  (*skelfile_) << sd[i] << " ";
	(*skelfile_) << "non" << std::endl << std::endl;	  
      }
    
  }

  
  void SkeletonConf::WritePDE ()
  {
#ifdef TRACE
    (*trace) << "Entering SkeletonConf::WritePDE" << std::endl;
#endif

    (*skelfile_) << "# -- PDE LIST (PLEASE FILL OUT)" << std::endl;
    (*skelfile_) << "list_pdes = XXX non" << std::endl << std::endl;

    (*skelfile_) << "# -- PDE SPECIFIC PARAMETERS" << std::endl;
    (*skelfile_) << "# -- name of pde:" << std::endl;
    (*skelfile_) << "XXX:" << std::endl << std::endl;

    (*skelfile_) << "\t subdomains = XXX non" << std::endl << std::endl;

    (*skelfile_) << " \t bc_conditions:" << std::endl;
    (*skelfile_) << " \t \t homogeneous_dirichlet = XXX non" << std::endl; 
    (*skelfile_) << "#\t \t If pde has more than one dof (e.g. mechanic), uncomment next line and substitute XXX with dof (e.g. ux):" << std::endl;
    (*skelfile_) << "#\t \t homoBCDof = XXX non" << std::endl << std::endl;
    
    (*skelfile_) << " \t \t inhomogeneous_dirichlet = XXX non" << std::endl; 
    (*skelfile_) << "#\t \t If you have an inhomogeneous_dirichlet BC , uncomment the following line and substitute IHD" << std::endl; 
    (*skelfile_) << "#\t \t \t XXX = XXXvalue [, XXXfile]" << std::endl;
    (*skelfile_) << "#\t \t If pde has more than one dof (e.g. mechanic), uncomment next line and substitute XXX with dof (e.g. ux):" << std::endl;
    (*skelfile_) << "#\t \t inhomoBCDof = XXX non" << std::endl << myEndl;
    

    (*skelfile_) << "\t loads = XXX non" << std::endl;
    (*skelfile_) << "#\t \t XXX = XXXValue [, XXXFile]" << std::endl;
    (*skelfile_) << "#\t \t If pde has more than one dof (e.g. mechanic), uncomment next line and substitute XXX with dof (e.g. ux):" << std::endl;
    (*skelfile_) << "#\t \t loadDof = XXX non" << std::endl << myEndl;
    

    (*skelfile_) << "# --\t ABSORBING BCs: (if yes, uncomment  bnd_for_absBCs and specify subdomain)" << std::endl; 
    (*skelfile_) << "#\t bnd_for_absBCs = non " << std::endl << std::endl;
  }
}





#endif
