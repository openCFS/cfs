#include <fstream>
#include <iostream>
#include <string>
#include <math.h>

#include <Estimator/actimeerror.hh>
#include <DataInOut/Unverg/outUnverg.hh>
#include <DataInOut/GMV/outGMV.hh>
#include <Forms/forms_header.hh>
#include <Estimator/spaceerror.hh>
#include <multigrid.hh>

#include <DataInOut/LoadMaterialData.hh>
#include <DataInOut/MaterialData.hh>
#include "mech2dPDE.hh" 

namespace CoupledField
{

  PlainStrainPDE::PlainStrainPDE(Grid * aptgrid, BCs *aptbcs, TimeFunc *aptTimeFunc, FileType *aptFileType, 
				 WriteResults *aptOut)
    :MechPDE(aptgrid,aptbcs,aptTimeFunc,aptFileType,aptOut)
  {
#ifdef TRACE
    (*trace) << "entering PlainStrainPDE::PlainStrainPDE " << std::endl;
#endif

    pdename_ = "plainStrain";
    pdematerialclass_ = "piezo";

    dofspernode_ = 2;

    size_=ptgrid_->GetMaxnumnodes(0)*dofspernode_;
    disp_.Resize(size_);
    disp_.Init(0);

    conf->getsubdompde(subdoms_,pdename_);
    ReadBCs(pdename_);

  }


  void PlainStrainPDE::SetupMatrices(const Integer level)
  {
#ifdef TRACE
    (*trace) << "entering PlainStrainPDE::SetupMatrices" << std::endl;
#endif

    Matrix<Double> elemmat;
    // This is a smaller matrix since it is just for linear 1D elements.
    Matrix<Double> elemmatbnd;
    BaseFE * ptEl;
    Vector<Integer> connecth;
    std::vector<Elem*> elemssd;
    Integer i, j;

    // read materialfile name from configfile
    std::string matFileName;
    conf->get("material_file", matFileName);
    char * charMatFileName = c_string(matFileName);
    
    // define space for material data
    LoadMaterialData loadMat(charMatFileName);
    MaterialData* matData = new MaterialData[subdoms_.size()];



    for (i=0; i<subdoms_.size(); i++)
      {
	// load material data into array "matData"
	std::vector<std::string> matName;
	conf->getliststr(subdoms_[i], matName, "list_subdomains");	
	loadMat.GetMaterial(matData[i], matName[i], "piezo");

	ptgrid_->GetElemSD(elemssd,subdoms_[i],level);

	for (j=0; j< elemssd.size(); j++)
	  {
	    ptEl = elemssd[j]->ptElem;

	    BaseForm * bilinear_mass  = new MassInt(ptEl);
	    BaseForm * bilinear_stiff = new mechPlainStrainInt(ptEl, &matData[i]);

	    connecth=elemssd[j]->connect;

	    Matrix<Double> ptCoord;
	    ptgrid_->GetCoordNodesElemMat(connecth, ptCoord, level); 

	    // stiffness part
	    bilinear_stiff->CalcElementMatrix(ptCoord, elemmat);
	    //	    elemmat *= coeffst[i];

#ifdef DEBUG
	    (*debug) << "Connection array  " << std::endl;
	    (*debug)  << connecth  << std::endl;
	    (*debug) << "Stiffnessmatrix, ElementNumber  " <<   i << std::endl;
	    (*debug) << elemmat << std::endl;
#endif     

	    algsys_->SetElementMatrix(elemmat.getinarray(), connecth.get(), connecth.size(), STIFFNESS);

	    // mass part
	    bilinear_mass->CalcElementMatrix(ptCoord, elemmat);
	    //	    elemmat *= coeffm[i];

#ifdef DEBUG
	    (*debug) << "Massmatrix, ElementNumber  " << i << std::endl;
	    (*debug) << elemmat << std::endl;
#endif
      
	    algsys_->SetElementMatrix(elemmat.getinarray(), connecth.get(), connecth.size(), MASS);
  
	    delete bilinear_stiff;
	    delete bilinear_mass;
	  }
      }

#ifdef TRACE
    (*trace) << "Leaving PlainStrainPDE::SetupMatrices" << std::endl;
#endif
  }


  void PlainStrainPDE::SetBCs(const Integer level, const Integer update, const Double atime)
{
#ifdef TRACE
  (*trace) << "entering Mech2dPDE::SetBCs" << std::endl;
#endif

  Integer node,dof;
  Double val;

  std::string dof_ux = "ux";
  std::string dof_uy = "uy";

  // set dirichlet boundary conditions
  Integer i,j=0;
  std::list<Integer> nodes;
  Integer sizebc=bcs_hd_.size();

  // homogeneous boundary conditions
  val = 0;
  for (i=0; i< bcs_hd_.size(); i++)
    {
      std::string doftype = bcs_hd_[i]; 
      if (doftype.substr(0,2)==dof_ux)
	dof = 1;
      else 
	if (doftype.substr(0,2)==dof_uy)
	dof = 2;
      else
	Error("Unknown dof-type in homog. BC; substring is:must start with ux or uy!!",__FILE__,__LINE__);
      
#ifdef DEBUG
      (*debug) << " Homog. Dirichlet BC" << std::endl;
#endif
      nodes=ptBCs_->GetNodesLevel(bcs_hd_[i]);
   
      for (std::list<Integer>::const_iterator p=nodes.begin(); p!=nodes.end(); p++, j++)
	{
	  node=*p;
	  if (update==1)
	    {
#ifdef DEBUG
	      (*debug) << " node: " << node << " dof:" << dof << " val: " << val << std::endl;
#endif
      
      algsys_->UpdateDirichlet(i+1, val, SYSTEM);

	    }
	  else
	    { 
#ifdef DEBUG
	      (*debug) << " node: " << node << " dof:" << dof << " val: " << val << std::endl;
#endif

      algsys_->SetDirichlet(j+1, node, val, dof, SYSTEM);

	    }
	}  
    }

  //inhomogeneous boundary conditions
  for (i=0; i<bcs_id_.size(); i++)
    {
      std::string doftype = bcs_hd_[i]; 
      if (doftype.substr(0,2)==dof_ux)
	dof = 1;
      else 
	if (doftype.substr(0,2)==dof_uy)
	dof = 2;
	else
	  Error("Unknown dof-type in inhomog. BC; substring is:must start with ux or uy!!",__FILE__,__LINE__);
      
#ifdef DEBUG
      (*debug) << " Inhomog. Dirichlet BC" << std::endl;
#endif

      nodes=ptBCs_->GetNodesLevel(bcs_id_[i], level);

      val=val_id_[i];

      for (std::list<Integer>::const_iterator p=nodes.begin(); p!=nodes.end(); p++, j++)
	{
	  node=*p;
          if (update==1)
            {

#ifdef DEBUG
	      (*debug) << " node: " << node << " dof:" << dof << " val: " << val << std::endl;
#endif

              algsys_->UpdateDirichlet(j+1, val, SYSTEM);
            }
          else
            {

#ifdef DEBUG
	      (*debug) << " node: " << node << " dof:" << dof << " val: " << val << std::endl;
#endif
              algsys_->SetDirichlet(j+1, node, val, dof, SYSTEM);
            }
	}
    }

#ifdef TRACE
  (*trace) << "leaving Mech2dPDE::SetBCs" << std::endl;
#endif
}


}


