// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <stddef.h>
#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <boost/algorithm/string/trim.hpp>

#include "CoupledPDE/pdecoupling.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/WriteInfo.hh"
#include "Domain/ansatzFct.hh"
#include "Domain/elem.hh"
#include "Domain/domain.hh"
#include "Domain/entityList.hh"
#include "Domain/grid.hh"
#include "Domain/resultInfo.hh"
#include "Driver/assemble.hh"
#include "Driver/formsContext.hh"
#include "Driver/solveStepSmooth.hh"
#include "Elements/basefe.hh"
#include "Forms/baseForm.hh"
#include "Forms/linElastInt.hh"
#include "General/exception.hh"
#include "PDE/SinglePDE.hh"
#include "PDE/StdPDE.hh"
#include "PDE/eqnMap.hh"
#include "PDE/timestepping.hh"
#include "Utils/StdVector.hh"
#include "Utils/baseelemstoresol.hh"
#include "Utils/result.hh"
#include "math.h"
#include "pseudoTS.hh" 
#include "extLBMPDE.hh"

namespace CoupledField {

class BaseMaterial;
class SingleVector;

ExtLBMPDE::ExtLBMPDE(Grid * aptgrid, PtrParamNode paramNode )
:SinglePDE(aptgrid, paramNode ) {

  pdename_ = "externLBM";
  pdematerialclass_ = MECHANIC;
  firstTurn_ = true;
  nonLin_ = false;

  method_ = "mechanic";
  // for 2D n_z=1
  StdVector<UInt> n = aptgrid[0].GetBoundaries(0);
  n_x = n[0];
  n_y = n[1];
  n_z = n[2];
  n_elems = n_x*n_y*n_z;

  //Initializing storage for PDFs
  pdfs = (Double *)malloc(sizeof(Double) * n_elems * _QN_);
  sim_type = myParam_->Get("LBM")->Get("type")->As<std::string>();
}

void ExtLBMPDE::DefineIntegrators() {
  // code taken over from mechPDE

  // help variables for parameter checking
  RegionIdType actRegion;
  BaseMaterial* actSDMat = NULL;

  std::map<RegionIdType, BaseMaterial*>::iterator it;
  for ( it = materials_.begin(); it != materials_.end(); it++ ) {
    // Set current region and material
    actRegion = it->first;
    actSDMat = it->second;

    // create new entity list
    shared_ptr<ElemList> actSDList( new ElemList(ptgrid_ ) );
    actSDList->SetRegion( actRegion );

    // get current region name and get grip of paramNode
    std::string actRegionName;
    actRegionName = ptgrid_->GetRegion().ToString( actRegion );

    // ==============  add "standard" linear stiffness ===========================
    //transform the type
    SubTensorType tensorType = PLANE_STRAIN;
    BaseForm * bilinearStiff = new linElastInt(actSDMat, tensorType);

    BiLinFormContext * actIntDescrStiff =
        new BiLinFormContext(bilinearStiff, STIFFNESS );

    actIntDescrStiff->SetPtPdes(this, this);
    actIntDescrStiff->SetResults( results_[0], results_[0],
        actSDList, actSDList );
    assemble_->AddBiLinearForm( actIntDescrStiff );
    // Give entities and result to equation numbering class
    // and solution class
    eqnMap_->AddResult( *results_[0], actSDList );
  }
}

void ExtLBMPDE::InitCoupling(PDECoupling * coupling){ }

void ExtLBMPDE::DefineSolveStep() {
  solveStep_ = new StdSolveStep(*this);
  if (sim_type == "external") {

    executable = myParam_->Get("LBM")->Get("lbm")->As<std::string>();

    this->writeData("CFS2LBM.dat");

    std::fstream f(executable.c_str());
    if (!f.is_open())
      EXCEPTION("Could not find executable '" + executable + "', might be not in path");
    f.close();

    std::cout << "++ Calling extern LBM solver .. \n" << std::endl;
    std::string command = executable + " -o LBM2CFS.dat CFS2LBM.dat";
    system(command.c_str());
//    int err = system(command.c_str());
//    if (err)
//      EXCEPTION("LBM simulation failed, no outputs available! \n");
//
//    std::cout << "\n++ Extern LBM simulation finished" << std::endl;

    this->readData("LBM2CFS.dat");
  }
  else
    assert(false);
}

void ExtLBMPDE::InitTimeStepping()
{
  // timestepping formulation
  TS_alg_ = new PseudoTS( algsys_);
}

void ExtLBMPDE::CalcOutputCoupling() {

  SolutionType quantity;
  StdVector<UInt> * couplingnodes;
  SingleVector * values;

  // at first, check if this PDE is iterative coupled
  if (isIterCoupled_ == false)
    return;

  // loop over all output coupling quantities
  for (UInt i=0; i<ptCoupling_->GetNumOutputCouplings(); i++) {
    quantity = ptCoupling_->GetOutputQuantity(i);

    switch(ptCoupling_->GetOutputType(i)) {

    case NODE:

      ptCoupling_->GetOutputNodes(i, couplingnodes);
      ptCoupling_->GetOutputValues(i, values);

      if (quantity == EXTERNLBM_VELOCITY) {
        solDeriv1_.SetAlgSysVector(getDeriv(FIRST_DERIV));
        solDeriv1_.NodeSolutionToCoupling((*values),*couplingnodes);
      }
      break;

    case ELEM:
      EXCEPTION("No Element coupling output");
    }

  }
}

bool ExtLBMPDE::HasOutput(SolutionType output) {

  if (output == EXTERNLBM_VELOCITY || output == EXTERNLBM_DENSITY)
    return true;
  return false;
}

void ExtLBMPDE::CalcDensities( shared_ptr<BaseResult> res ) {
  Result<double>& result = dynamic_cast<Result<double>&>(*res);
  EntityIterator it = result.GetEntityList()->GetIterator();
  assert(it.GetType() == EntityList::ELEM_LIST);
  Vector<double>& dens = result.GetVector();
  dens.Resize(n_elems);
  dens.Init(0);

  for (UInt i = 0; i < n_elems; ++i) {
    for (int j = 0; j < 9; ++j) {
          dens[i] += PDF_IDX(i,j);
        }
  }
}

void ExtLBMPDE::CalcResults( shared_ptr<BaseResult> res ) {
  SolutionType solType = res->GetResultInfo()->resultType ;
  switch (solType) {
  case EXTERNLBM_DENSITY:
    CalcDensities(res);
    break;
  case EXTERNLBM_VELOCITY:
    CalcVelocities(res);
    break;
  default:
    WARN( "Result type not computable by externLBM PDE" );
    break;
  }
}

void ExtLBMPDE::CalcVelocities( shared_ptr<BaseResult> res ) {
  Result<double>& result = dynamic_cast<Result<double>&>(*res);
  EntityIterator it = result.GetEntityList()->GetIterator();
  assert(it.GetType() == EntityList::ELEM_LIST);
  Vector<double>& velo = result.GetVector();
  velo.Resize(n_elems*dim_);

//  shared_ptr<ResultInfo> resinfo = GetResultInfo(EXTERNLBM_DENSITY);
  ElemList single_elem(domain->GetGrid());
//  Vector<double> elem_sol;
  double lux, luy;
  double density;
  int elemId;
  // traverse the elements
  for(it.Begin(); !it.IsEnd(); it++)
  {
    density = 0;
    const Elem* elem = it.GetElem();
    elemId = elem->elemNum;
    single_elem.SetElement(elem);
    for (int i = 0; i < 9; ++i)
      density += PDF_IDX(elemId,i);

//    GetSolVecOfElement(elem_sol, single_elem.GetIterator(), resinfo);
//    assert(elem_sol.GetSize() == n_elems);
//    lux   = (PDF_IDX(elemId, 1) + PDF_IDX(elemId, 5) + PDF_IDX(elemId, 8)
//            - PDF_IDX(elemId, 3) - PDF_IDX(elemId, 6) - PDF_IDX(elemId, 7)) / elem_sol[elemId];
//    luy   = (PDF_IDX(elemId, 2) + PDF_IDX(elemId, 5) + PDF_IDX(elemId, 6)
//            - PDF_IDX(elemId, 4) - PDF_IDX(elemId, 7) - PDF_IDX(elemId, 8)) / elem_sol[elemId];
    lux   = (PDF_IDX(elemId, 1) + PDF_IDX(elemId, 5) + PDF_IDX(elemId, 8)
        - PDF_IDX(elemId, 3) - PDF_IDX(elemId, 6) - PDF_IDX(elemId, 7)) / density;
    luy   = (PDF_IDX(elemId, 2) + PDF_IDX(elemId, 5) + PDF_IDX(elemId, 6)
        - PDF_IDX(elemId, 4) - PDF_IDX(elemId, 7) - PDF_IDX(elemId, 8)) / density;
    velo[elemId * dim_] = lux;
    velo[elemId * dim_ + 1] = luy;
  }
}

// ***********************************************************************
//   Obtain information on desired output quantities from parameter file
// ***********************************************************************
void ExtLBMPDE::ReadSpecialResults() {
}


void ExtLBMPDE::DefineAvailResults() {
  // =====================================================================
  // set solution information
  // =====================================================================
  PtrParamNode resultsNode = myParam_->Get("storeResults", ParamNode::PASS );

  // === MECHANIC DISPLACEMENT (dummy result type)===
  shared_ptr<ResultInfo> disp(new ResultInfo);
  StdVector<std::string> dispDofNames;
  dispDofNames = "x", "y";
  disp->resultType = MECH_DISPLACEMENT;
  disp->dofNames = dispDofNames;
  disp->unit = "m";
  if( subType_ != "flatShell" ) {
    disp->entryType = ResultInfo::VECTOR;
  } else {
    disp->entryType = ResultInfo::TENSOR;
  }
  shared_ptr<AnsatzFct> fct(new LagrangeFct);
  disp->fctType = fct;
  disp->definedOn = ResultInfo::NODE;
  results_.Push_back( disp );
  availResults_.insert( disp);

  // === macroscopic LBM density ===
  shared_ptr<ResultInfo> dens(new ResultInfo);
  dens->resultType = EXTERNLBM_DENSITY;
  dens->unit =  "kg/m^3";
  dens->dofNames = "";
  dens->entryType = ResultInfo::SCALAR;
  dens->definedOn = ResultInfo::ELEMENT;
  dens->fctType = shared_ptr<ConstFct>(new ConstFct() );
//  results_.Push_back( dens );
  availResults_.insert( dens );

  // === macroscopic LBM velocity ===
  StdVector<std::string> velDofNames;
  if (dim_ == 2)
    velDofNames = "x", "y";
  else
    velDofNames = "x", "y", "z";
  shared_ptr<ResultInfo> velo(new ResultInfo);
  velo->resultType = EXTERNLBM_VELOCITY;
  velo->dofNames = velDofNames;
  velo->unit =  "m/s";
  velo->entryType = ResultInfo::VECTOR;
  velo->definedOn = ResultInfo::ELEMENT;
  velo->fctType = shared_ptr<ConstFct>(new ConstFct() );
  availResults_.insert( velo );

}

// ***********************************************************************
//   Read LBM velocities from LBM output file
// ***********************************************************************
void ExtLBMPDE::readData(const char * filename) {
  FILE * file;
  char line[1000];
  file = fopen(filename,"r");

  if (file == NULL) {
    std::string ex = "Error in opening file" + std::string(filename);
    EXCEPTION(ex);
  }

  int i = 0;
  while (fgets(line, sizeof line, file)){

    if (*line == '#')
      continue;

    if (sscanf(line, "%15le%15le%15le%15le%15le%15le%15le%15le%15le",
        &PDF_IDX(i,0),&PDF_IDX(i,1),&PDF_IDX(i,2),&PDF_IDX(i,3),&PDF_IDX(i,4),
        &PDF_IDX(i,5),&PDF_IDX(i,6),&PDF_IDX(i,7),&PDF_IDX(i,8)) != 9)
        EXCEPTION("There's something wrong with data in LBM2CFS.dat!")
    ++i;
  }
}

// ***********************************************************************
//  Write data file for extern LBM solver according to interface schema
// ***********************************************************************
void ExtLBMPDE::writeData(const char * filename){

  FILE * file;
  file = fopen(filename,"w");

  if (file == NULL) {
    std::string ex = "Error in opening file" + std::string(filename);
    EXCEPTION(ex);
  }

  fprintf(file,"# 1.0 # interface version \n");
  fprintf(file,"# optional comment line: the cell resolution \n");
  fprintf(file,"%u \n",n_x);
  fprintf(file,"%u \n",n_y);
  fprintf(file,"%u \n",n_z);
  fprintf(file,"# omega: relaxation parameter  \n");
  fprintf(file,"%.15e \n",myParam_->Get("LBM")->Get("omega")->As<Double>());
  fprintf(file,"# maximal walltime [sec] \n");
  fprintf(file,"%.15e \n", myParam_->Get("LBM")->Get("maxWallTime")->As<Double>());
  fprintf(file,"# maximal number of iterations \n");
  fprintf(file,"%u \n",myParam_->Get("LBM")->Get("maxIter")->As<Integer>());
  fprintf(file,"# convergence tolerance \n");
  fprintf(file,"%.15e \n",myParam_->Get("LBM")->Get("epsilon")->As<Double>());

  // find out which elements are inlet, outlet, boundary and inner ones
  StdVector<Elem*> elems;
  StdVector<Elem*> boundaries;
  StdVector<Elem*> inlets;
  StdVector<Elem*> outlets;
  // auxiliary vector
  StdVector<std::string> elemsClass(n_elems);
  // vector initialized with porosity value of inner cells
  elemsClass.Init("0.5");
  Grid* grd = domain->GetGrid();
  StdVector<std::string> regionNames;
  grd->GetRegionNames(regionNames);
  //specify regionId for boundary and inner region
  int boundId = -1;
  int innerId = -1;

  for (UInt i = 0; i<regionNames.GetSize(); ++i) {
    if (regionNames[i] == "boundary")
      boundId = i;
    else
      innerId = i;
  }
  if (boundId == -1 || innerId == -1)
    EXCEPTION("Mesh must contain 2 regions: 'boundary' and another region");

  grd->GetElems(boundaries,boundId);
  grd->GetElems(elems,innerId);
  grd->GetElemsByName(inlets,"inlet");
  grd->GetElemsByName(outlets,"outlet");
  for (UInt i = 0; i < boundaries.GetSize(); ++i) {
    elemsClass[boundaries[i]->elemNum - 1] = "-1.0";
  }
  for (UInt i = 0; i < inlets.GetSize(); ++i) {
    elemsClass[inlets[i]->elemNum - 1] = "-2.0";
  }
  for (UInt i = 0; i < outlets.GetSize(); ++i) {
    elemsClass[outlets[i]->elemNum - 1] = "-3.0";
  }

  // write out type of element corresponding to element id
  fprintf(file,"# domain: nx * ny * nz entries with element id \n # -3.0: outlet \n # -2.0: inlet \n # -1.0: bounce-back \n");
  for (UInt i = 0; i < n_elems; ++i) {
    fprintf(file,"%s \n",elemsClass[i].c_str());
  }

  fprintf(file,"# only for the inlet cells the space separated vectors  \n");
  // estimate u_x, u_y and u_z from xml file
  ParamNodeList inletNodes = myParam_->Get("bcsAndLoads")->GetList("inlet");
  double u_x(0.0), u_y(0.0), u_z(0.0);
  for (UInt i = 0; i < inletNodes.GetSize(); ++i) {
    std::string dof = inletNodes[i]->Get("dof")->As<std::string>();
    if (dof == "x")
      u_x = inletNodes[i]->Get("value")->As<Double>();
    else if (dof == "y")
      u_y = inletNodes[i]->Get("value")->As<Double>();
    if (dim_ == 3)
      u_z = inletNodes[i]->Get("value")->As<Double>();
  }
  for (UInt i = 0; i < inlets.GetSize(); ++i) {
    fprintf(file,"%.15e %.15e %.15e \n", u_x, u_y, u_z);
  }
  fclose(file);

  std::cout << "++ CFS2LBM.dat created" << std::endl;
}

}
