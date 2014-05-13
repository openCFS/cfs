// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <stddef.h>
#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/predicate.hpp>

#include "CoupledPDE/pdecoupling.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/WriteInfo.hh"
#include "DataInOut/Logging/cfslog.hh"
#include "DataInOut/Logging/log.hpp"
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

DECLARE_LOG(lbm)
DEFINE_LOG(lbm, "lbm")

ExtLBMPDE::ExtLBMPDE(Grid* grid, PtrParamNode pn)
:SinglePDE(grid, pn) {

  pdename_ = "externLBM";
  pdematerialclass_ = MECHANIC;
  firstTurn_ = true;
  nonLin_ = false;

  method_ = "mechanic";

  // LBM parametes
  omega_       = myParam_->Get("LBM/omega")->As<double>();
  maxWallTime_ = myParam_->Get("LBM/maxWallTime")->As<double>();
  maxIter_     = myParam_->Get("LBM/maxIter")->As<unsigned int>();
  convergence_     = myParam_->Get("LBM/convergence")->As<double>();

  PtrParamNode bcsl = myParam_->Get("bcsAndLoads");
  // only one inlet region
  u_x_ = bcsl->HasByVal("inlet", "dof", "x") ? bcsl->GetByVal("inlet", "dof", "x")->Get("value")->As<double>() : 0.0;
  u_y_ = bcsl->HasByVal("inlet", "dof", "y") ? bcsl->GetByVal("inlet", "dof", "y")->Get("value")->As<double>() : 0.0;
  u_z_ = bcsl->HasByVal("inlet", "dof", "z") ? bcsl->GetByVal("inlet", "dof", "z")->Get("value")->As<double>() : 0.0;

  iface.SetName("ExtLBMPDE::Iface");
  iface.Add(INTERNAL, "internal");
  iface.Add(EXT_MATLAB, "external_matlab_iface");
  iface.Add(EXT_CFSxLBM, "external_single_file_iface");
  iface_ = iface.Parse(pn->Get("LBM/solver")->As<std::string>());

  InitRegions(pn, grid);

  // for 2D n_z_=1
  StdVector<UInt> n = grid->GetBoundaries(boundary_reg_);
  n_x_ = n[0];
  n_y_ = n[1];
  n_z_ = n[2];
  n_elems = n_x_ * n_y_ * n_z_;

  SetupElementMapping(grid);

  //Initializing storage for PDFs
  pdfs.Resize(n_elems * 9);
}

void ExtLBMPDE::InitRegions(PtrParamNode pn, Grid* grid)
{
  // find the unique boundary region id
  ParamNodeList regions = pn->Get("regionList")->GetChildren();
  if(regions.GetSize() < 2)
    throw Exception("externLBM requires at least two regions where one has the boundary attribute set");

  boundary_reg_ = -1;
  design_reg_.Reserve(regions.GetSize() - 1); // ommit boundary region

  for(unsigned int i = 0; i < regions.GetSize(); i++)
  {
    RegionIdType reg = grid->GetRegion().Parse(regions[i]->Get("name")->As<std::string>());

    if(regions[i]->Get("boundary")->As<bool>())
    {
      if(boundary_reg_ != -1)
        throw Exception("only a single region my have the boundary attribute set");
      else
        boundary_reg_ = reg;
    }
    else
      design_reg_.Push_back(reg);
  }
  if(boundary_reg_ == -1)
    throw Exception("externLBM requires a region with boundary attribute set");
}

void ExtLBMPDE::SetupElementMapping(Grid* grid)
{
  assert(!design_reg_.IsEmpty());
  assert(n_elems == n_x_ * n_y_ * n_z_);

  // TODO here we assume that the whole mesh is for LBM and the mesh is of lexicographic ordering.
  // To be good this needs to handled by element neighbours!
  if(grid->GetNumElems() != n_elems)
    EXCEPTION("the current implementation assums the whole mesh to used for LBM and lexicographic ordered. Mesh has " << grid->GetNumElems() << " but we assume " << n_elems);

  idx_to_elem.Resize(n_elems);
  for(unsigned int i = 0; i < n_elems; i++)
    idx_to_elem[i] = i+1; // one based

  elem_to_idx.Resize(n_elems + 1); // one-based elem_nr
  for(unsigned int i = 0, n = elem_to_idx.GetSize(); i < n; i++)
    elem_to_idx[i] = i-1; // -1 for non appropriate idx
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

void ExtLBMPDE::InitCoupling(PDECoupling * coupling)
{
}

void ExtLBMPDE::DefineSolveStep()
{
  solveStep_ = new StdSolveStep(*this);

  switch(iface_)
  {
  case EXT_MATLAB:
  case EXT_CFSxLBM:
  {
    executable = myParam_->Get("LBM")->Get("lbm")->As<std::string>();

    ExportExternalSolverFiles();

    std::fstream f(executable.c_str());
    if (!f.is_open())
      EXCEPTION("Could not find executable '" + executable + "', might be not in path");
    f.close();

    std::cout << "++ Calling extern LBM solver .. \n" << std::endl;
    std::string command = executable + " -o LBM2CFS.dat CFS2LBM.dat";
    int err = system(command.c_str());
    if (err)
      EXCEPTION("LBM simulation failed, no outputs available! \n");

    if(iface_ == EXT_CFSxLBM)
      ReadData("LBM2CFS.dat");
    else
      ReadData("node_steady.dat");
    break;
  }
  case INTERNAL:
    assert(false);
  }
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

void ExtLBMPDE::CalcResults( shared_ptr<BaseResult> res )
{
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

double ExtLBMPDE::CalcLBMDensity(unsigned int idx) const
{

  double density = 0.0;
  for(unsigned int h = 0; h < 9; h++)
    density += pdf(idx,h);

  return density;
}


void ExtLBMPDE::CalcDensities( shared_ptr<BaseResult> base_result )
{
  Result<double>& res = dynamic_cast<Result<double>&>(*base_result);

  EntityIterator it = res.GetEntityList()->GetIterator();
  assert(it.GetType() == EntityList::ELEM_LIST);

  Vector<double>& val = res.GetVector();
  val.Resize(res.GetEntityList()->GetSize());

  // traverse the elements
  for(it.Begin(); !it.IsEnd(); it++)
    val[it.GetPos()] = CalcLBMDensity(elem_to_idx[it.GetElem()->elemNum]);
}


void ExtLBMPDE::CalcVelocities( shared_ptr<BaseResult> base_result )
{
  Result<double>& res = dynamic_cast<Result<double>&>(*base_result);

  EntityIterator it = res.GetEntityList()->GetIterator();
  assert(it.GetType() == EntityList::ELEM_LIST);

  Vector<double>& val = res.GetVector();
  val.Resize(res.GetEntityList()->GetSize() * dim_);

  // traverse the elements
  for(it.Begin(); !it.IsEnd(); it++)
  {
    unsigned int idx = elem_to_idx[it.GetElem()->elemNum];
    double density = CalcLBMDensity(idx);

    double ux   = (pdf(idx, 1) + pdf(idx, 5) + pdf(idx, 8) - pdf(idx, 3) - pdf(idx, 6) - pdf(idx, 7)) / density;
    double uy   = (pdf(idx, 2) + pdf(idx, 5) + pdf(idx, 6) - pdf(idx, 4) - pdf(idx, 7) - pdf(idx, 8)) / density;

    val[it.GetPos() * dim_] = ux;
    val[it.GetPos() * dim_ + 1] = uy;
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
void ExtLBMPDE::ReadData(const std::string& filename)
{
  std::ifstream file(filename.c_str());
  if(file.fail())
    throw Exception("cannot open lbm result file " + std::string(filename));

  std::string line;
  unsigned int i = 0;
  while(std::getline(file, line))
  {
    if(!boost::starts_with(line, "#"))
    {
      std::istringstream ss(line);
      //if(!(ss >>  PDF_IDX(i,0) >>  PDF_IDX(i,1) >>  PDF_IDX(i,2) >>  PDF_IDX(i,3) >>  PDF_IDX(i,4) >>  PDF_IDX(i,5) >>  PDF_IDX(i,6) >>  PDF_IDX(i,7) >>  PDF_IDX(i,8)))
      if(!(ss >>  pdf(i,0) >>  pdf(i,1) >>  pdf(i,2) >>  pdf(i,3) >>  pdf(i,4) >>  pdf(i,5) >>  pdf(i,6) >>  pdf(i,7) >>  pdf(i,8)))
        EXCEPTION("error reading nine values in line " << (i+1) << " of file " << filename);
      i++;
    }
  }

  if(i != n_elems)
    EXCEPTION("read " << i << " lines in " << filename << " but excpected " << n_elems);


  file.close();

  // LOG_DBG3(lbm) << "RD(" << filename << ") -> " << StdVector<double>::ToString(n_elems * 9, pdfs, 0);
}


void ExtLBMPDE::ExportCFS2LBM(const StdVector<double>& elements, const StdVector<Elem*>& inlets)
{
  // the single new interface file
  std::ofstream file("CFS2LBM.dat");

  file << "# 1.0 # interface version \n";
  file << "# optional comment line: the cell resolution \n";
  file << n_x_ << std::endl;
  file << n_y_ << std::endl;
  file << n_z_ << std::endl;
  file << "# omega: relaxation parameter  \n";
  file << omega_ << std::endl;
  file << "# maximal walltime [sec] \n";
  file << maxWallTime_ << std::endl;
  file << "# maximal number of iterations \n";
  file << maxIter_ << std::endl;
  file << "# convergence tolerance \n";
  file << convergence_ << std::endl;

  // write out type of element corresponding to element id
  file << "# domain: nx * ny * nz entries with element id \n# -3.0: outlet \n# -2.0: inlet \n# -1.0: bounce-back \n";
  for (UInt i = 0; i < n_elems; ++i)
    file << elements[i] << std::endl;

  file << "# only for the inlet cells the space separated vectors  \n";

  for (UInt i = 0; i < inlets.GetSize(); ++i)
    file <<  u_x_ << " " << u_y_ << " " << u_z_ << std::endl;

  file.close();

  std::cout << "++ CFS2LBM.dat created" << std::endl;
}

void ExtLBMPDE::ExportMultipleFiles(const StdVector<double>& elements, const StdVector<Elem*>& inlets)
{

  int num_non_sing = 712;

  std::ofstream data("data.dat");
  data << n_x_ << std::endl;
  data << n_y_ << std::endl;
  data << 3.0 << std::endl; // penalty parameter
  data << u_x_ << std::endl;
  data << u_y_ << std::endl;
  data << 1.0 << std::endl;  // density at the outlet
  data << omega_ << std::endl; // relaxation parameter for collision step
  data << maxIter_ << std::endl;
  data << convergence_ << std::endl; // lbm convergence tolerance
  data << num_non_sing << std::endl; // number of lines in non_sing.dat
  data << 1 << std::endl; // id of objective (1=pressure drop)
  data.close();

  std::ofstream por("por.dat");
  assert(n_elems == elements.GetSize());
  for(unsigned int i = 0, n = elements.GetSize(); i < n; i++)
    por << std::max(0.0, elements[i]) << std::endl; // boundary, inlet, outlet is all 0.0 in the old interface
  por.close();

  std::ofstream obst("obst.dat");
  // assume to have the mesh ordered from lower left starting to the right (lexicographic)
  // the OLD! interface assumes origen left upper and x downwards and y to the right!
  // TODO! Do by neighbours!
  Matrix<int> obst_m(n_y_, n_x_); // FIXME check order!
  for(unsigned int x = 0; x < n_x_; x++)
  {
    for(unsigned int y = 0; y < n_y_; y++)
    {
      // org: -1 bb, -2 inlet, -3 outlet, 0 ... 1 porisity
      // out: 1 bb1, 2 inlet, 3 outlet, 0 for porosity
      int v = (int) -1.0 * elements[n_x_ * y + x];
      obst << std::max(0, v) << " "; // org porosity would be -1
      obst_m(y, x) = v;  // origin of the matrix is left upper!
    }
    obst << std::endl;
  }
  obst.close();

  LOG_DBG2(lbm) << "EMF: obst matrix (!= obst.dat):\n" << obst_m.ToString(0, true);

  StdVector<unsigned int> non_sing;
  IdentifyNonSingualrityIndices(obst_m, non_sing);

  std::ofstream ns("non_sing.dat");
  // add one to be one based
  for(unsigned int i = 0, n = non_sing.GetSize(); i < n; i++)
    ns << (non_sing[i] + 1) << " " << std::endl; // space to allow diff with matlab's non_sing.dat
  ns.close();

}


void ExtLBMPDE::ExportExternalSolverFiles()
{
  // find out which elements are inlet, outlet, boundary and inner ones
  StdVector<Elem*> elems;
  StdVector<Elem*> boundaries;
  StdVector<Elem*> inlets;
  StdVector<Elem*> outlets;
  // auxiliary vector
  StdVector<double> elements(n_elems);
  // vector initialized with porosity value of inner cells
  elements.Init(0.5); // TODO: initial design
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
  // FIXME: check if index by elemNum ist nonsense!!
  for (UInt i = 0; i < boundaries.GetSize(); ++i)
    elements[boundaries[i]->elemNum - 1] = -1.0;
  for (UInt i = 0; i < inlets.GetSize(); ++i)
    elements[inlets[i]->elemNum - 1] = -2.0;
  for (UInt i = 0; i < outlets.GetSize(); ++i)
    elements[outlets[i]->elemNum - 1] = -3.0;

  // estimate u_x, u_y and u_z from xml file
  ParamNodeList inletNodes = myParam_->Get("bcsAndLoads")->GetList("inlet");


  if(iface_ == EXT_CFSxLBM)
    ExportCFS2LBM(elements, inlets);
  else
    ExportMultipleFiles(elements, inlets);

}

void ExtLBMPDE::IdentifyNonSingualrityIndices(Matrix<int>& obst, StdVector<unsigned int>& indices)
{
  // reimplemented matlab code
  // we have the vector to export it as non_sing.dat for the old interface
  indices.Resize(9 * n_x_ * n_y_, 0);

  // FIXME again check the proper oder.
  unsigned int k = 0;
  for(unsigned int j = 0; j < n_y_; j++)
  {
    for(unsigned int i = 0; i < n_x_; i++)
    {
      int val = obst[j][i];
      unsigned int base = j * n_x_ * 9 + i * 9; // 2D
      LOG_DBG3(lbm) << "INSI: test j=" << j << " i=" << i << " k=" << k << " base=" << base << " val -> " << val;
      if(val  != 1) // no bounce back
      {
        assert(val == 0 || val == 2 || val == 3);
        // fluid node distribution functions are not deleted
        for(unsigned int h = 0; h < 9; h++, k++) // 2D
          indices[k] = base + h;
      }
      else
      {
        // outpointing boundary distributions are not inserted and
        // distributions which point towards a bounce-back boundary point

        // test, if the node in direction f_1 (x+1) is not boundary but interior, inlet or outlet
        if(i+1 < n_x_ && obst[j][i+1] != 1) // if (inew<=lx && obst(inew,j)~=1)
        {
          indices[k] = base + 1; // x(k)=(j-1)*lx*9+(i-1)*9+2;
          LOG_DBG3(lbm) << "INSI: case 1 (x+1): j=" << j << " i=" << i << " k=" << k << " next=" << obst[j][i+1] << " -> " << indices[k];
          k++;
        }
        if(j+1 < n_y_ && obst[j+1][i] != 1) // j+1 = f_2_
        {
          indices[k] = base + 2;
          k++;
        }
        if(i > 0 && obst[j][i-1] != 1)
        {
          indices[k] = base + 3;
          k++;
        }
        if(j > 0 && obst[j-1][i] != 1)
        {
          indices[k] = base + 4;
          k++;
        }
        if(i+1 < n_x_ && j+1 < n_y_ && obst[j+1][i+1] != 1)
        {
          indices[k] = base + 5;
          k++;
        }
        if(i > 0 && j+1 < n_y_ && obst[j+1][i-1] != 1)
        {
          indices[k] = base + 6;
          k++;
        }

        if(i > 0 && j > 0 && obst[j-1][i-1] != 1)
        {
          indices[k] = base + 7;
          k++;
        }
        if(i+1 < n_x_ && j > 0 && obst[j-1][i+1] != 1)
        {
          indices[k] = base + 8;
          k++;
        }
      }
    }

  }


  // restrict indices
  indices.Resize(k);

}



}
