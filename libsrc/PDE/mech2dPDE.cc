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
#include <Forms/bdbInt.hh>

namespace CoupledField
{

  Mech2dPDE::Mech2dPDE(Grid * aptgrid, BCs *aptbcs, TimeFunc *aptTimeFunc, FileType *aptFileType, 
		       WriteResults *aptOut)
    :MechPDE(aptgrid,aptbcs,aptTimeFunc,aptFileType,aptOut)
  {
#ifdef TRACE
    (*trace) << "entering Mech2dPDE::Mech2dPDE " << std::endl;
#endif

    pdename_ = "mech";
    pdematerialclass_ = "piezo";

    conf->getstr("subtype", subType_, pdename_ );

    if (subType_ == "3d")
      dofspernode_ = 3;
    else
      dofspernode_ = 2;

    size_=ptgrid_->GetMaxnumnodes(0)*dofspernode_;
    disp_.Resize(size_);
    disp_.Init(0);

    conf->getsubdompde(subdoms_,pdename_);
    ReadBCs(pdename_);
  }


  void Mech2dPDE::SetupMatrices(const Integer level)
  {
#ifdef TRACE
    (*trace) << "entering Mech2dPDE::SetupMatrices" << std::endl;
#endif

    Matrix<Double> elemmat;
    // This is a smaller matrix since it is just for linear 1D elements.
    Matrix<Double> elemmatbnd;
    BaseFE * ptEl;
    Vector<Integer> connecth;
    std::vector<Elem*> elemssd;
    Integer i, j;

    for (i=0; i<subdoms_.size(); i++)
      {	
	ptgrid_->GetElemSD(elemssd,subdoms_[i],level);

	for (j=0; j< elemssd.size(); j++)
	  {
	    ptEl = elemssd[j]->ptElem;

	    BaseForm * bilinear_mass  = new MassInt(ptEl, materialData_[i]);
	    BaseForm * bilinear_stiff;
	    
	    if (subType_ == "plainStrain")
	      bilinear_stiff = new mechPlainStrainInt(ptEl, materialData_[i]);
	    else if (subType_ == "3d")
	      bilinear_stiff = new mech3DInt(ptEl, materialData_[i]);
	    else 
	      {
		std::string errMessg;
		errMessg = "Unknown subtype \"" + subType_ + "\" in mech2d PDE!";		
		Error(errMessg.c_str(),__FILE__,__LINE__);
	      }


	    connecth=elemssd[j]->connect;

	    Matrix<Double> ptCoord;
	    ptgrid_->GetCoordNodesElemMat(connecth, ptCoord, level); 

	    // stiffness part
	    bilinear_stiff->CalcElementMatrix(ptCoord, elemmat);

	    algsys_->SetElementMatrix(elemmat.getinarray(), connecth.get(), connecth.size(), SYSTEM);


	    (*debug) << "Mech3d elemmat: " << std::endl << elemmat << std::endl;

	    // 	    // mass part
	    // 	    bilinear_mass->CalcElementMatrix(ptCoord, elemmat);

	    // 	    Matrix <Double> elemMatMultDof;
	    // 	    MassMultiDof(elemMatMultDof, elemmat, dofspernode_);

	    // 	    algsys_->SetElementMatrix(elemMatMultDof.getinarray(), connecth.get(), connecth.size(), MASS);

  
	    delete bilinear_stiff;
	    delete bilinear_mass;
	  }
      }

#ifdef DEBUG
    algsys_->Print(SYSTEM);
#endif
    
    
#ifdef TRACE
    (*trace) << "Leaving Mech2dPDE::SetupMatrices" << std::endl;
#endif
  }






  void Mech2dPDE::SetBCs(const Integer level, const Integer update, const Double atime)
  {
#ifdef TRACE
    (*trace) << "entering Mech2dPDE::SetBCs" << std::endl;
#endif
    
    Integer node,dof;
    Double val;
    
    // set dirichlet boundary conditions
    Integer i,j=0;
    std::list<Integer> nodes;
    Integer sizebc=bcs_hd_.size();
    
    std::vector<Elem*> edgesBC;  // vector of 1D-elements from mesh-file
    Vector<Integer> connecth;
	  

    // homogeneous boundary conditions
    val = 0;
    for (i=0; i< bcs_hd_.size(); i++)
      {
	std::string doftype = bcs_hd_[i]; 
	dof = GetNrBCDof (doftype.substr(0,2));
      
#ifdef DEBUG
	(*debug) << " Homog. Dirichlet BC" << std::endl;
#endif
	nodes=ptBCs_->GetNodesLevel(bcs_hd_[i]);
   
	for (std::list<Integer>::const_iterator p=nodes.begin(); p!=nodes.end(); p++, j++)
	  {
	    node=*p;

#ifdef DEBUG
	    (*debug) << " node: " << node << " dof:" << dof << " val: " << val << std::endl;
#endif
	    if (update==1)
	      algsys_->UpdateDirichlet(i+1, val, SYSTEM);
	    else
	      algsys_->SetDirichlet(j+1, node, val, dof, SYSTEM);
	  }  
      }

    //inhomogeneous boundary conditions
    for (i=0; i<bcs_id_.size(); i++)
      {
	std::string doftype = bcs_id_[i]; 
	dof = GetNrBCDof (doftype.substr(0,2));
      
#ifdef DEBUG
	(*debug) << " Inhomog. Dirichlet BC" << std::endl;
#endif

	nodes = ptBCs_->GetNodesLevel(bcs_id_[i], level);

	val=val_id_[i];

	for (std::list<Integer>::const_iterator p=nodes.begin(); p!=nodes.end(); p++, j++)
	  {
	    node=*p;
#ifdef DEBUG
	    (*debug) << " node: " << node << " dof:" << dof << " val: " << val << std::endl;
#endif
	  
	    if (update==1)
	      algsys_->UpdateDirichlet(j+1, val, SYSTEM);
	    else	    
	      algsys_->SetDirichlet(j+1, node, val, dof, SYSTEM);
	  }
      }


    // load boundary conditions
    for (i=0; i<bcs_loads_.size(); i++)
      {
	std::string doftype = bcs_loads_[i]; 

	dof = GetNrBCDof (doftype.substr(0,2));

#ifdef DEBUG
	(*debug) << " Load BC" << std::endl;
#endif

	nodes = ptBCs_->GetNodesLevel(bcs_loads_[i], level);

	val = val_loads_[i];

	for (std::list<Integer>::const_iterator p=nodes.begin(); p!=nodes.end(); p++, j++)
	  {
	    node=*p;

#ifdef DEBUG
	    (*debug) << " node: " << node << " dof:" << dof << " val: " << val << std::endl;
#endif
	  
	    val = val_loads_[i];
	    algsys_->SetNodeRHS(val, node, dof);
	  
	  }
      }    
#ifdef TRACE
    (*trace) << "leaving Mech2dPDE::SetBCs" << std::endl;
#endif
  }



  Integer Mech2dPDE::GetNrBCDof(const std::string & dofStartString)
  {
#ifdef TRACE
    (*trace) << "leaving Mech2dPDE::GetNrBCDof" << std::endl;
#endif

    Integer nrActDof;
    
    if (dofStartString == "ux")
      nrActDof = 1;
    else 
      if (dofStartString == "uy")
	nrActDof = 2;
      else
	if (dofspernode_ == 3)
	  if (dofStartString == "uz")
	    nrActDof = 3;
	  else
	    Error("Unknown dof-type in homog. BC; substring must start with ux, uy or uz!!",__FILE__,__LINE__);
	else
	  Error("Unknown dof-type in homog. BC; substring must start with ux or uy!!",__FILE__,__LINE__);
    
    return nrActDof;
  }


} // end namespace CoupledField



